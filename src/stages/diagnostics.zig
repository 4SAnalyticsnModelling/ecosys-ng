//! Mass-balance reconstruction and debug census helpers.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
pub fn reconstructLandscapeMassBalance(context: anytype) !ecosys.mass_balance_audit.Totals {
    return ecosys.landscape_mass_balance_runtime.reconstruct(.{
        .grid = context.grid,
        .plants = context.plants,
        .snow = context.snow_transport,
        .soil_thermal = context.soil_thermal,
        .soil_properties = context.soil_solver_properties,
        .soil_gas = context.gas_transport,
        .root_gas = if (context.plant_roots.*) |*roots| roots else null,
        .soil_organic = context.soil_organic,
        .soil_organic_transport = context.soil_organic_transport,
        .surface_organic = context.surface_organic,
        .mineral_nitrogen = context.mineral_nitrogen_transport,
        .soil_chemistry = context.soil_chemistry,
        .nitrogen_fertilizer = context.soil_fertilizer_inventory,
        .mineral_fertilizer = context.mineral_fertilizer_inventory,
        .micropore_solutes = context.micropore_solute_state,
        .macropore_solutes = context.macropore_solute_state,
        .surface_chemistry = context.surface_litter_chemistry,
        .surface_fertilizer = context.surface_litter_fertilizer,
        .surface_denitrification_nitrite_g_n = context.surface_denitrification.nitrite_g_n,
        .surface = context.surface_precipitation,
        .surface_ice_water_equivalent_m3 = context.surface_litter_ice_m3,
        .surface_gas = context.litter_gas_transport,
        .surface_litter_dry_mass_megagrams = context.surface_litter_geometry.dry_mass_megagrams,
        .canopy_retention = if (context.canopy_precipitation_retention.*) |*value| value else null,
        .cell_area_m2 = context.canopy_cell_area_m2,
        .soil_mass_megagrams_scratch = context.landscape_soil_mass_megagrams_scratch,
        .parameters = .{
            .snow_ice_density_megagrams_per_m3 = context.runscript.snow_ice_density_megagrams_per_m3,
            .carbon_g_per_mol = 12,
            .nitrogen_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .phosphorus_g_per_mol = context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
            .phosphate_zone_fractions = .{
                .ammonium_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .ammonium_band = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .nitrate_band = context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                .phosphate_band = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
            },
            .surface_physical = .{
                .dry_organic_heat_capacity_megajoules_per_g_c_k = context.runscript.surface_pond_dry_organic_heat_capacity_megajoules_per_g_c_k,
                .liquid_water_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
                .ice_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.ice_heat_capacity_megajoules_per_m3_k,
                .water_molar_mass_g_per_mol = context.runscript.soil_gas_transport_parameters.water_molar_mass_g_per_mol,
                .liquid_water_density_g_per_m3 = context.runscript.soil_gas_transport_parameters.water_density_g_per_m3,
            },
        },
    }, context.landscape_boundary_ledger);
}

pub fn diagnosticStoredNitrogen_g(context: anytype) !f64 {
    const totals = try reconstructLandscapeMassBalance(context);
    return totals.residue_nitrogen_g + totals.organic_nitrogen_g +
        totals.dinitrogen_nitrogen_g + totals.ammonium_nitrogen_g +
        totals.nitrate_nitrogen_g;
}

pub fn diagnosticStoredPhosphorus_g(context: anytype) !f64 {
    const totals = try reconstructLandscapeMassBalance(context);
    return totals.residue_phosphorus_g + totals.organic_phosphorus_g +
        totals.phosphate_phosphorus_g;
}

pub fn diagnosticPhosphorusOwners_g(context: anytype) ![3]f64 {
    const totals = try reconstructLandscapeMassBalance(context);
    return .{ totals.residue_phosphorus_g, totals.organic_phosphorus_g, totals.phosphate_phosphorus_g };
}

pub fn diagnosticRelayerPhosphateOwners_g(context: anytype) [3]f64 {
    const p_mass = context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol;
    const band_fraction = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction;
    const fractions = [2]f64{ 1 - band_fraction, band_fraction };
    var result: [3]f64 = @splat(0);
    for (0..context.grid.layer_count) |layer| {
        const water = context.grid.matrix_liquid_water_m3[layer];
        const soil_mass = context.soil_solver_properties.matrix_bulk_volume_m3[layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[layer];
        const zones = [2]@TypeOf(context.soil_chemistry.non_band_phosphate[0]){ context.soil_chemistry.non_band_phosphate[layer], context.soil_chemistry.band_phosphate[layer] };
        for (zones, fractions) |zone, fraction| {
            result[0] += fraction * water * (zone.dissolved_hpo4_mol_p_per_m3 + zone.dissolved_h2po4_mol_p_per_m3) * p_mass;
            result[1] += fraction * soil_mass * (zone.adsorbed_hpo4_mol_p_per_megagram + zone.adsorbed_h2po4_mol_p_per_megagram) * p_mass;
            result[2] += fraction * water * (zone.aluminum_phosphate_solid_mol_per_m3 + zone.iron_phosphate_solid_mol_per_m3 + zone.dicalcium_phosphate_solid_mol_per_m3 + 3 * zone.hydroxyapatite_solid_mol_per_m3 + 2 * zone.monocalcium_phosphate_solid_mol_per_m3) * p_mass;
        }
    }
    return result;
}

pub fn diagnosticAmmoniumOwners_g_n(context: anytype) ![6]f64 {
    const molar_mass = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol;
    var result: [6]f64 = @splat(0);
    for (0..context.grid.cell_count) |cell| {
        const surface = context.surface_litter_chemistry.cells[cell];
        result[0] += context.surface_precipitation.litter_water_m3[cell] * (surface.ammonium_mol_per_m3 + surface.ammonia_mol_per_m3) * molar_mass;
        result[1] += context.surface_litter_geometry.dry_mass_megagrams[cell] * surface.exchange.ammonium_mol_per_megagram * molar_mass;
        const fertilizer = context.surface_litter_fertilizer.cells[cell];
        result[2] += (fertilizer.ammonium_mol_n + fertilizer.ammonia_mol_n + fertilizer.urea_mol_n) * molar_mass;
    }
    for (0..context.grid.layer_count) |layer| {
        const matrix = try context.mineral_nitrogen_transport.matrix.cellAmountsConst(layer);
        const macro = try context.mineral_nitrogen_transport.macropore.cellAmountsConst(layer);
        inline for ([_]ecosys.mineral_nitrogen_transport.Species{ .ammonium_non_band, .ammonium_band, .ammonia_non_band, .ammonia_band }) |species| result[3] += (matrix[@intFromEnum(species)] + macro[@intFromEnum(species)]) * molar_mass;
        const exchange = context.soil_chemistry.cation_exchange_mol_per_megagram[layer];
        result[4] += context.soil_solver_properties.matrix_bulk_volume_m3[layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[layer] * (exchange.ammonium_non_band + exchange.ammonium_band) * molar_mass;
        const fertilizer = context.soil_fertilizer_inventory.soil[layer];
        result[5] += (fertilizer.broadcast_ammonium_mol_n + fertilizer.broadcast_ammonia_mol_n + fertilizer.broadcast_urea_mol_n + fertilizer.banded_ammonium_mol_n + fertilizer.banded_ammonia_mol_n + fertilizer.banded_urea_mol_n) * molar_mass;
    }
    return result;
}

pub fn diagnosticSnowNitrogen_g(storage: ecosys.landscape_mass_inventory.Storage) f64 {
    return storage.dinitrogen_nitrogen_g + storage.ammonium_nitrogen_g + storage.nitrate_nitrogen_g;
}

pub fn diagnosticSnowSpeciesNitrogen_g(amounts: []const f64) f64 {
    var total_g_n: f64 = 0;
    const species_count = ecosys.snow_solute_transport.species_count;
    var first: usize = 0;
    while (first < amounts.len) : (first += species_count) {
        inline for ([_]ecosys.snow_solute_transport.Species{
            .dinitrogen_nitrogen,
            .nitrous_oxide_nitrogen,
            .ammonium_nitrogen,
            .ammonia_nitrogen,
            .nitrate_nitrogen,
        }) |species| total_g_n += amounts[first + @intFromEnum(species)];
    }
    return total_g_n;
}

pub fn diagnosticSnowDischargeNitrogen_g(discharge: []const ecosys.snow_solute_transport.SurfaceDischarge) f64 {
    var total_g_n: f64 = 0;
    for (discharge) |cell| {
        total_g_n += diagnosticSnowSpeciesNitrogen_g(&cell.litter_g);
        total_g_n += diagnosticSnowSpeciesNitrogen_g(&cell.soil_nonband_g);
        total_g_n += diagnosticSnowSpeciesNitrogen_g(&cell.soil_band_g);
    }
    return total_g_n;
}

fn debugCarbonComponents(context: anytype, totals: ecosys.mass_balance_audit.Totals, label: []const u8) !void {
    const surface_organic_storage = try ecosys.landscape_mass_inventory.aggregateSurfaceOrganic(
        context.surface_organic,
    );
    const soil_organic_storage = try ecosys.landscape_mass_inventory.aggregateSoilOrganic(
        context.soil_organic,
        context.grid,
    );
    const area = totals.landscape_area_m2;
    std.log.debug("{s} components: surf_resid={e} soil_resid={e} soil_org={e} co2={e} cum_co2in={e} cum_cout={e}", .{
        label,
        surface_organic_storage.residue_carbon_g / area,
        soil_organic_storage.residue_carbon_g / area,
        soil_organic_storage.organic_carbon_g / area,
        totals.carbon_dioxide_carbon_g / area,
        totals.cumulative_carbon_dioxide_input_g / area,
        totals.cumulative_carbon_output_g / area,
    });
    // Fine-grained surface sub-pools (cell=0 only)
    const so = context.surface_organic;
    const msub = ecosys.soil_organic_initialization.microbial_substrate_count;
    const mpop = ecosys.soil_organic_initialization.microbial_population_count;
    const mfrac = ecosys.soil_organic_initialization.kinetic_fraction_count;
    const sub = ecosys.soil_organic_initialization.substrate_count;
    const rfrac = ecosys.soil_organic_initialization.residue_fraction_count;
    const sfrac = ecosys.soil_organic_initialization.structural_fraction_count;
    var surf_microbial: f64 = 0;
    for (0..msub) |s| {
        if (s == 4) continue;
        const first = (s * mpop) * mfrac;
        for (so.microbial[first .. first + mpop * mfrac]) |p| surf_microbial += p.carbon_g_c;
    }
    var surf_residue: f64 = 0;
    for (0..3) |s| {
        const first = s * rfrac;
        for (so.residue[first .. first + rfrac]) |p| surf_residue += p.carbon_g_c;
        surf_residue += so.dissolved[s].carbon_g_c + so.adsorbed[s].carbon_g_c;
        surf_residue += so.dissolved_acetate_carbon_g_c[s] + so.adsorbed_acetate_carbon_g_c[s];
    }
    var surf_structural: f64 = 0;
    for (so.structural[0 .. sub * sfrac]) |p| surf_structural += p.carbon_g_c;
    // Also log the EXCLUDED surface pools (substrate 3-4 residue/dissolved/adsorbed, substrate 4 microbial)
    var surf_excluded_resid: f64 = 0;
    for (3..sub) |s| { // substrates 3 and 4
        const first = s * rfrac;
        for (so.residue[first .. first + rfrac]) |p| surf_excluded_resid += p.carbon_g_c;
        surf_excluded_resid += so.dissolved[s].carbon_g_c + so.adsorbed[s].carbon_g_c;
        surf_excluded_resid += so.dissolved_acetate_carbon_g_c[s] + so.adsorbed_acetate_carbon_g_c[s];
    }
    var surf_microbial_sub4: f64 = 0;
    {
        const first = (4 * mpop) * mfrac;
        for (so.microbial[first .. first + mpop * mfrac]) |p| surf_microbial_sub4 += p.carbon_g_c;
    }
    std.log.debug("{s} surf subpools: microbial={e} residue_fracs={e} structural={e} excluded_resid={e} excluded_micro4={e}", .{
        label,                      surf_microbial / area,      surf_residue / area, surf_structural / area,
        surf_excluded_resid / area, surf_microbial_sub4 / area,
    });
    // Fine-grained soil microbial by substrate category (active layers only)
    const grid = context.grid;
    var soil_microbial_resid: f64 = 0;
    var soil_microbial_org: f64 = 0;
    for (0..grid.cell_count) |cell| {
        const active = grid.active_soil_layer_count[cell];
        for (0..active) |layer| {
            const layer_cell = cell * grid.soil_layer_capacity + layer;
            const first_m = layer_cell * msub * mpop * mfrac;
            for (0..msub) |s| {
                const m_first = first_m + s * mpop * mfrac;
                for (context.soil_organic.microbial[m_first .. m_first + mpop * mfrac]) |p| {
                    if (s == 4) soil_microbial_org += p.carbon_g_c else soil_microbial_resid += p.carbon_g_c;
                }
            }
        }
    }
    std.log.debug("{s} soil microbial: resid={e} org={e}", .{
        label, soil_microbial_resid / area, soil_microbial_org / area,
    });
    // Separate gas-CO2 from bicarbonate/carbonate to diagnose balance drift.
    var gas_co2_g: f64 = 0;
    var gas_ch4_g: f64 = 0;
    const gas_state_dbg = context.gas_transport;
    for (0..grid.cell_count) |cell| {
        const active = grid.active_soil_layer_count[cell];
        for (0..active) |layer| {
            const layer_cell = cell * grid.soil_layer_capacity + layer;
            const first = layer_cell * ecosys.gas_transport.species_count;
            const end = first + ecosys.gas_transport.species_count;
            const gaseous = gas_state_dbg.gaseous_mass_g[first..end];
            const dissolved = gas_state_dbg.dissolved_mass_g[first..end];
            const macropore = gas_state_dbg.macropore_dissolved_mass_g[first..end];
            const band = gas_state_dbg.band_dissolved_mass_g[first..end];
            gas_co2_g += gaseous[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] +
                dissolved[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] +
                macropore[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] +
                band[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)];
            gas_ch4_g += gaseous[@intFromEnum(ecosys.gas_transport.Species.methane)] +
                dissolved[@intFromEnum(ecosys.gas_transport.Species.methane)] +
                macropore[@intFromEnum(ecosys.gas_transport.Species.methane)] +
                band[@intFromEnum(ecosys.gas_transport.Species.methane)];
        }
    }
    // Bicarbonate/carbonate from micropore + macropore solute amounts.
    var bicarbonate_g: f64 = 0;
    const micro = context.micropore_solute_state;
    const macro = context.macropore_solute_state;
    for (0..grid.cell_count) |cell| {
        const active = grid.active_soil_layer_count[cell];
        for (0..active) |layer| {
            const layer_cell = cell * grid.soil_layer_capacity + layer;
            const micro_amounts = try micro.cellAmountsConst(layer_cell);
            const macro_amounts = try macro.cellAmountsConst(layer_cell);
            inline for (@typeInfo(ecosys.solute_transport_species.AqueousSpecies).@"enum".fields) |field| {
                const species: ecosys.solute_transport_species.AqueousSpecies = @enumFromInt(field.value);
                const is_carbonate_carrier = switch (species) {
                    .carbonate,
                    .bicarbonate,
                    .calcium_carbonate,
                    .calcium_bicarbonate,
                    .magnesium_carbonate,
                    .magnesium_bicarbonate,
                    .sodium_carbonate,
                    => true,
                    else => false,
                };
                if (is_carbonate_carrier) {
                    bicarbonate_g += (micro_amounts[field.value] + macro_amounts[field.value]) * 12.0;
                }
            }
        }
    }
    // Chemistry CO2 and bicarbonate — CO2 is synced to gas_transport.dissolved after h-step,
    // bicarbonate is exported to micropore solute after h-step. Both are counted in EXEC
    // AFTER the first h-step, but at baseline (before any h-step) they are NOT yet counted.
    var chem_co2_g: f64 = 0;
    var chem_bicarb_g: f64 = 0;
    for (0..grid.cell_count) |cell| {
        const active = grid.active_soil_layer_count[cell];
        for (0..active) |layer| {
            const layer_cell = cell * grid.soil_layer_capacity + layer;
            const water_m3 = grid.matrix_liquid_water_m3[layer_cell];
            chem_co2_g += context.soil_chemistry.aqueous[layer_cell].carbon_dioxide * 12.0 * water_m3;
            chem_bicarb_g += context.soil_chemistry.aqueous[layer_cell].bicarbonate * 12.0 * water_m3;
        }
    }
    // Litter chemistry CO2 — not synced to litter_gas_transport, not in EXEC balance.
    var litter_chem_co2_g: f64 = 0;
    var litter_chem_bicarb_g: f64 = 0;
    var litter_gas_co2_g: f64 = 0;
    {
        const litter_chem = context.surface_litter_chemistry;
        const litter_water = context.surface_precipitation.litter_water_m3;
        for (0..grid.cell_count) |cell| {
            litter_chem_co2_g += litter_chem.cells[cell].carbon_dioxide_mol_per_m3 * 12.0 * litter_water[cell];
            litter_chem_bicarb_g += litter_chem.cells[cell].bicarbonate_mol_per_m3 * 12.0 * litter_water[cell];
        }
    }
    {
        const litter_gas = context.litter_gas_transport;
        const co2_species = @intFromEnum(ecosys.gas_transport.Species.carbon_dioxide);
        for (0..grid.cell_count) |cell| {
            const base = cell * ecosys.gas_transport.species_count;
            litter_gas_co2_g += litter_gas.gaseous_mass_g[base + co2_species] +
                litter_gas.dissolved_mass_g[base + co2_species];
        }
    }
    std.log.debug("{s} co2 split: gas_co2={e} gas_ch4={e} bicarbonate={e} chem_co2={e} chem_bicarb={e} litter_chem_co2={e} litter_chem_bicarb={e} litter_gas_co2={e}", .{
        label,                    gas_co2_g / area,            gas_ch4_g / area,        bicarbonate_g / area, chem_co2_g / area, chem_bicarb_g / area,
        litter_chem_co2_g / area, litter_chem_bicarb_g / area, litter_gas_co2_g / area,
    });
}
