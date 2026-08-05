const std = @import("std");
const gas_transport = @import("gas_transport.zig");
const grid_module = @import("grid.zig");
const root_exchange = @import("plant_root_gas_exchange.zig");
const root_system = @import("plant_root_system.zig");
const plant_water = @import("plant_water_balance.zig");
const soil_properties = @import("soil_solver_properties.zig");

pub const Settings = struct {
    liquid_tortuosity_coefficient: f64,
    absolute_tolerance_g: f64,
    relative_tolerance: f64,
    maximum_iterations: u16,
};

pub fn advance(
    roots: *root_system.State,
    water: *const plant_water.Workspace,
    grid: *const grid_module.GridState,
    properties: *const soil_properties.State,
    soil_gas: *gas_transport.State,
    parameters: root_exchange.RuntimeParameters,
    atmospheric_concentration_g_per_m3: [gas_transport.species_count]f64,
    biological_domain_count_by_plant: []const u8,
    active_by_plant: []const bool,
    settings: Settings,
) !void {
    try validate(roots, water, grid, properties, soil_gas, biological_domain_count_by_plant, active_by_plant, settings);
    for (0..grid.cell_count) |cell| {
        const plant_first = cell * (roots.plant_count / grid.cell_count);
        const species_per_cell = roots.plant_count / grid.cell_count;
        for (0..grid.active_soil_layer_count[cell]) |layer| {
            const soil = try grid.layerIndex(cell, layer);
            const soil_water_m3 = grid.matrix_liquid_water_m3[soil];
            if (soil_water_m3 <= 0) continue;
            const bulk_volume_m3 = properties.matrix_bulk_volume_m3[soil];
            const water_fraction = if (bulk_volume_m3 > 0) std.math.clamp(soil_water_m3 / bulk_volume_m3, 0, 1) else 0;
            const tortuosity = settings.liquid_tortuosity_coefficient * water_fraction * water_fraction;
            const film_m = try root_exchange.soilWaterFilmThicknessM(parameters, grid.matric_potential_mpa[soil]);
            for (0..species_per_cell) |species| {
                const plant = plant_first + species;
                if (!active_by_plant[plant] or roots.roots_dead[plant]) continue;
                const population = water.plant_population_count[plant];
                if (population <= 0) continue;
                const domain_count = biological_domain_count_by_plant[plant];
                for (0..domain_count) |domain| {
                    const root = try roots.layerIndex(plant, domain, layer);
                    if (roots.aqueous_volume_m3[root] <= settings.absolute_tolerance_g or
                        roots.root_surface_area_m2_per_plant[root] <= settings.absolute_tolerance_g)
                        continue;
                    inline for (non_oxygen_gases) |gas| {
                        if (domain == 0 or gas == .carbon_dioxide) {
                            try exchangeSoilAqueous(
                                roots,
                                water,
                                grid,
                                soil_gas,
                                parameters,
                                root,
                                soil,
                                gas,
                                soil_water_m3,
                                tortuosity,
                                film_m,
                                population,
                                false,
                            );
                            if (gas == .ammonia) try exchangeSoilAqueous(
                                roots,
                                water,
                                grid,
                                soil_gas,
                                parameters,
                                root,
                                soil,
                                gas,
                                soil_water_m3,
                                tortuosity,
                                film_m,
                                population,
                                true,
                            );
                        }
                    }
                    // UPTAKE RCO2A enters the root aqueous pool after the
                    // soil-root concentration-gradient transaction and before
                    // root phase exchange.
                    const respiration_g_c = roots.actual_respiration_g_c_per_h[root];
                    const next_aqueous_co2_g_c = roots.aqueous_carbon_dioxide_g_c[root] + respiration_g_c;
                    const next_reaction_g_c = roots.aqueous_carbon_dioxide_reaction_g_c_per_h[root] + respiration_g_c;
                    if (!std.math.isFinite(respiration_g_c) or respiration_g_c < 0 or
                        !std.math.isFinite(next_aqueous_co2_g_c) or !std.math.isFinite(next_reaction_g_c))
                        return error.NonFiniteRootRespirationGasSource;
                    roots.aqueous_carbon_dioxide_g_c[root] = next_aqueous_co2_g_c;
                    roots.aqueous_carbon_dioxide_reaction_g_c_per_h[root] = next_reaction_g_c;
                }

                // The source permits root internal gas transport only in the
                // primary biological domain.
                const root = try roots.layerIndex(plant, 0, layer);
                if (roots.gaseous_volume_m3[root] <= settings.absolute_tolerance_g) continue;
                const topology = try roots.layerTopology(
                    plant,
                    0,
                    layer,
                    population,
                    properties.layer_thickness_m[soil],
                    water.woody_root_fraction[plant],
                    roots.average_secondary_length_m[root],
                );
                const cross_section_m = try root_exchange.rootGasCrossSectionPerLengthM(
                    population,
                    topology.primary_axis_count,
                    topology.secondary_axis_count,
                    water.primary_root_radius_m[plant],
                    water.secondary_root_radius_m[plant],
                    properties.layer_thickness_m[soil],
                    topology.average_secondary_length_m,
                    water.root_biome_fraction[root],
                );
                inline for (non_oxygen_gases) |gas| {
                    const environment = try root_exchange.gasEnvironment(parameters, gas, grid.soil_temperature_k[soil], 0);
                    const gas_species = transportSpecies(gas);
                    const gas_slot = @intFromEnum(gas);
                    const root_gas = rootGaseousPool(roots, gas);
                    const root_aqueous = rootAqueousPool(roots, gas);
                    const result = try root_exchange.solveRootPhaseAtmosphereExchangeG(.{
                        .gaseous_mass_g = root_gas[root],
                        .aqueous_mass_g = root_aqueous[root],
                        .gaseous_volume_m3 = roots.gaseous_volume_m3[root],
                        .aqueous_volume_m3 = roots.aqueous_volume_m3[root],
                        .water_to_air_mass_solubility_ratio = environment.water_to_air_mass_solubility_ratio,
                        .atmosphere_concentration_g_per_m3 = atmospheric_concentration_g_per_m3[@intFromEnum(gas_species)],
                        .atmosphere_conductance_m3_per_h = environment.gaseous_diffusivity_m2_per_h * cross_section_m,
                        .phase_equilibration_fraction = std.math.clamp(
                            roots.current_porosity_fraction_by_domain[try roots.domainIndex(plant, 0)],
                            0,
                            1,
                        ),
                        .maximum_iterations = settings.maximum_iterations,
                        .absolute_tolerance_g = settings.absolute_tolerance_g,
                        .relative_tolerance = settings.relative_tolerance,
                    });
                    const ledger = root * root_system.transported_root_gas_count + gas_slot;
                    const next_phase = roots.aqueous_to_gaseous_root_exchange_g_per_h[ledger] + result.aqueous_to_gaseous_exchange_g_per_h;
                    const next_atmosphere = roots.atmosphere_to_root_gas_exchange_g_per_h[ledger] + result.atmosphere_to_root_exchange_g_per_h;
                    if (!std.math.isFinite(next_phase) or !std.math.isFinite(next_atmosphere)) return error.NonFiniteRootGasLedger;
                    root_gas[root] = result.final_gaseous_mass_g;
                    root_aqueous[root] = result.final_aqueous_mass_g;
                    roots.aqueous_to_gaseous_root_exchange_g_per_h[ledger] = next_phase;
                    roots.atmosphere_to_root_gas_exchange_g_per_h[ledger] = next_atmosphere;
                }
            }
        }
    }
    try roots.validateFinite();
}

const non_oxygen_gases = [_]root_exchange.TransportedGas{
    .carbon_dioxide,
    .methane,
    .nitrous_oxide,
    .ammonia,
    .hydrogen,
};

fn exchangeSoilAqueous(
    roots: *root_system.State,
    water: *const plant_water.Workspace,
    grid: *const grid_module.GridState,
    soil_gas: *gas_transport.State,
    parameters: root_exchange.RuntimeParameters,
    root: usize,
    soil: usize,
    gas: root_exchange.TransportedGas,
    soil_water_m3: f64,
    tortuosity: f64,
    film_m: f64,
    population: f64,
    band: bool,
) !void {
    const species = transportSpecies(gas);
    const component = try gas_transport.massIndex(soil, species, soil_gas.cell_count);
    const soil_pool = if (band) &soil_gas.band_dissolved_mass_g[component] else &soil_gas.dissolved_mass_g[component];
    const root_pool = &rootAqueousPool(roots, gas)[root];
    const biome_fraction = water.root_biome_fraction[root];
    if (biome_fraction <= 0 or soil_pool.* <= 0 and root_pool.* <= 0) return;
    const environment = try root_exchange.gasEnvironment(parameters, gas, grid.soil_temperature_k[soil], 0);
    const conductance = try root_exchange.radialAqueousConductanceM3PerH(
        environment.aqueous_diffusivity_m2_per_h,
        tortuosity,
        water.root_surface_area_per_radius_m[root],
        water.root_cylinder_radius_m[root],
        film_m,
    );
    const exchange_g = try root_exchange.soilRootAqueousExchangeG(.{
        .soil_dissolved_mass_g = soil_pool.* * biome_fraction,
        .root_dissolved_mass_g = root_pool.*,
        .soil_water_volume_m3 = soil_water_m3 * biome_fraction,
        .root_water_volume_m3 = roots.aqueous_volume_m3[root],
        .water_advection_m3_per_step = @max(0, -roots.water_uptake_m3_per_h[root]) / population,
        .aqueous_diffusive_conductance_m3_per_step = conductance,
        .plant_population_count = population,
        .equilibration_fraction = 1,
    });
    const ledger = root * root_system.transported_root_gas_count + @intFromEnum(gas);
    try root_exchange.commitSoilRootAqueousExchangeWithLedgerG(
        soil_pool,
        root_pool,
        &roots.soil_to_root_gas_exchange_g_per_h[ledger],
        exchange_g,
    );
    if (gas == .ammonia) {
        if (band) {
            roots.ammonia_band_soil_exchange_g_n_per_h[root] += exchange_g;
        } else {
            roots.ammonia_nonband_soil_exchange_g_n_per_h[root] += exchange_g;
        }
    }
}

fn transportSpecies(gas: root_exchange.TransportedGas) gas_transport.Species {
    return switch (gas) {
        .carbon_dioxide => .carbon_dioxide,
        .methane => .methane,
        .nitrous_oxide => .nitrous_oxide,
        .ammonia => .ammonia,
        .hydrogen => .hydrogen,
        .oxygen => .oxygen,
    };
}

fn rootGaseousPool(roots: *root_system.State, gas: root_exchange.TransportedGas) []f64 {
    return switch (gas) {
        .carbon_dioxide => roots.gaseous_carbon_dioxide_g_c,
        .methane => roots.gaseous_methane_g_c,
        .nitrous_oxide => roots.gaseous_nitrous_oxide_g_n,
        .ammonia => roots.gaseous_ammonia_g_n,
        .hydrogen => roots.gaseous_hydrogen_g_h,
        .oxygen => roots.gaseous_oxygen_g_o,
    };
}

fn rootAqueousPool(roots: *root_system.State, gas: root_exchange.TransportedGas) []f64 {
    return switch (gas) {
        .carbon_dioxide => roots.aqueous_carbon_dioxide_g_c,
        .methane => roots.aqueous_methane_g_c,
        .nitrous_oxide => roots.aqueous_nitrous_oxide_g_n,
        .ammonia => roots.aqueous_ammonia_g_n,
        .hydrogen => roots.aqueous_hydrogen_g_h,
        .oxygen => roots.aqueous_oxygen_g_o,
    };
}

fn validate(
    roots: *const root_system.State,
    water: *const plant_water.Workspace,
    grid: *const grid_module.GridState,
    properties: *const soil_properties.State,
    soil_gas: *const gas_transport.State,
    domains: []const u8,
    active: []const bool,
    settings: Settings,
) !void {
    if (grid.cell_count == 0 or roots.plant_count % grid.cell_count != 0 or water.plant_count != roots.plant_count or
        roots.soil_layer_count != grid.soil_layer_capacity or soil_gas.cell_count != grid.layer_count or
        properties.layer_thickness_m.len != grid.layer_count or domains.len != roots.plant_count or active.len != roots.plant_count)
        return error.RootGasTransportDimensionMismatch;
    if (!std.math.isFinite(settings.liquid_tortuosity_coefficient) or settings.liquid_tortuosity_coefficient < 0 or
        !std.math.isFinite(settings.absolute_tolerance_g) or settings.absolute_tolerance_g <= 0 or
        !std.math.isFinite(settings.relative_tolerance) or settings.relative_tolerance <= 0 or settings.maximum_iterations == 0)
        return error.InvalidRootGasTransportSettings;
    for (domains) |count| if (count == 0 or count > root_system.biological_domain_count) return error.RootGasTransportDimensionMismatch;
}
