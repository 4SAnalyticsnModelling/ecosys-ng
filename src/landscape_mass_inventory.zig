const std = @import("std");
const snow = @import("snow_solute_transport.zig");
const grid_module = @import("grid.zig");
const gas = @import("gas_transport.zig");
const organic = @import("soil_organic_initialization.zig");
const organic_transport = @import("soil_organic_transport.zig");
const litter_chemistry = @import("surface_litter_chemistry.zig");
const litter_fertilizer = @import("surface_litter_fertilizer.zig");
const audit = @import("mass_balance_audit.zig");
const surface_precipitation = @import("surface_precipitation.zig");
const canopy_retention = @import("canopy_precipitation_retention.zig");
const mineral_nitrogen = @import("mineral_nitrogen_transport.zig");
const nitrogen_fertilizer = @import("fertilizer_nitrogen_inventory.zig");
const mineral_fertilizer = @import("mineral_fertilizer_inventory.zig");
const soil_chemistry = @import("solute_chemistry_state.zig");
const solute_transport = @import("solute_transport.zig");
const solute_species = @import("solute_transport_species.zig");
const zone_classification = @import("solute_charge_classification.zig");
const plant_roots = @import("plant_root_system.zig");

/// Authoritative storage side of the seven EXEC conservation equations.
/// Boundary additions/removals remain in their process-owned cumulative
/// ledgers and are combined by `mass_balance_audit`.
pub const Storage = struct {
    water_m3: f64 = 0,
    heat_megajoules: f64 = 0,
    oxygen_g: f64 = 0,
    residue_carbon_g: f64 = 0,
    organic_carbon_g: f64 = 0,
    carbon_dioxide_carbon_g: f64 = 0,
    residue_nitrogen_g: f64 = 0,
    organic_nitrogen_g: f64 = 0,
    dinitrogen_nitrogen_g: f64 = 0,
    ammonium_nitrogen_g: f64 = 0,
    nitrate_nitrogen_g: f64 = 0,
    residue_phosphorus_g: f64 = 0,
    organic_phosphorus_g: f64 = 0,
    phosphate_phosphorus_g: f64 = 0,
    ion_inventory_mol: f64 = 0,

    pub fn add(self: *Storage, contribution: Storage) !void {
        inline for (std.meta.fields(Storage)) |field| {
            const value = @field(self, field.name) + @field(contribution, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteLandscapeInventory;
            @field(self, field.name) = value;
        }
    }

    pub fn validate(self: Storage) !void {
        inline for (std.meta.fields(Storage)) |field| {
            const value = @field(self, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteLandscapeInventory;
            // `heat_megajoules` is an enthalpy measured against liquid water at
            // 0 K (HEAT-001 resolution A), so a sufficiently frozen carrier may
            // legitimately be negative. Every other field is a mass or volume
            // and must remain nonnegative.
            if (comptime !std.mem.eql(u8, field.name, "heat_megajoules"))
                if (value < 0) return error.NegativeLandscapeInventory;
        }
    }
};

/// Publishes only authoritative storage fields into EXEC totals. Cumulative
/// boundary ledgers and landscape area remain owned by the caller.
/// HEAT-001 instrumentation. Per-carrier frozen water equivalent accumulated
/// by the most recent aggregation pass, printed once per census so the
/// opening and closing daily censuses can be compared carrier by carrier.
var diagnostic_snow_solid_water_equivalent_m3: f64 = 0;
var diagnostic_snow_ice_volume_water_equivalent_m3: f64 = 0;
var diagnostic_soil_matrix_ice_water_equivalent_m3: f64 = 0;
var diagnostic_soil_macropore_ice_water_equivalent_m3: f64 = 0;
var diagnostic_surface_ice_water_equivalent_m3: f64 = 0;

pub fn publishStorage(
    totals: *audit.Totals,
    storage: Storage,
) !void {
    try storage.validate();
    std.log.debug("heat latent census: snow_swe_m3={e} snow_ice_we_m3={e} soil_matrix_ice_m3={e} soil_macropore_ice_m3={e} surface_ice_m3={e}", .{
        diagnostic_snow_solid_water_equivalent_m3,
        diagnostic_snow_ice_volume_water_equivalent_m3,
        diagnostic_soil_matrix_ice_water_equivalent_m3,
        diagnostic_soil_macropore_ice_water_equivalent_m3,
        diagnostic_surface_ice_water_equivalent_m3,
    });
    totals.water_storage_m3 = storage.water_m3;
    totals.heat_storage_megajoules = storage.heat_megajoules;
    totals.oxygen_storage_g = storage.oxygen_g;
    totals.residue_carbon_g = storage.residue_carbon_g;
    totals.organic_carbon_g = storage.organic_carbon_g;
    totals.carbon_dioxide_carbon_g = storage.carbon_dioxide_carbon_g;
    totals.residue_nitrogen_g = storage.residue_nitrogen_g;
    totals.organic_nitrogen_g = storage.organic_nitrogen_g;
    totals.dinitrogen_nitrogen_g = storage.dinitrogen_nitrogen_g;
    totals.ammonium_nitrogen_g = storage.ammonium_nitrogen_g;
    totals.nitrate_nitrogen_g = storage.nitrate_nitrogen_g;
    totals.residue_phosphorus_g = storage.residue_phosphorus_g;
    totals.organic_phosphorus_g = storage.organic_phosphorus_g;
    totals.phosphate_phosphorus_g = storage.phosphate_phosphorus_g;
    totals.ion_inventory_mol = storage.ion_inventory_mol;
}

/// REDIST snowpack inventory:
/// WS = VOLSSL + VOLWSL + VOLVSL + VOLISL*DENSI
/// ENGYW = VHCPW*TKW, with gas and nutrient carriers summed in their tracked
/// element units. The modern snow owner retains only the eight transported
/// ion carriers, so their authoritative mole inventory is reconstructed from
/// the explicitly tracked element mass rather than inventing absent species.
///
/// HEAT-001 resolution A. The audited heat quantity is an *enthalpy*, not a
/// sensible heat, so every frozen carrier also contributes its latent heat of
/// fusion. The reference state is liquid water at 0 K: liquid and vapor carry
/// only `C*T`, and frozen water carries `C*T - L*V_water_equivalent`. Freeze
/// and thaw are then internal conversions that cancel exactly in the census,
/// and no boundary booking of latent heat of fusion is required or permitted.
/// Default latent heat of fusion in MJ m-3, used only by the deprecated
/// two-argument `aggregateSnow` wrapper below. See
/// `docs/binding_requests/heat_001_landscape_enthalpy.md`.
pub const default_latent_heat_of_fusion_megajoules_per_m3: f64 = 333;

/// HEAT-001 second layer. The pure-water melting temperature that separates
/// the ice branch of the enthalpy curve from the liquid branch.
///
/// The reference state is liquid water at 0 K, so a cubic metre of frozen
/// water carries
///
///     C_liquid*Tm - L + C_ice*(T - Tm)
///
/// and NOT `C_ice*T - L`. The two differ by `(C_liquid - C_ice)*Tm`, which is
/// `(4.19 - 1.9274)*273.15 = 617.5 MJ m-3`. That is nearly twice the latent
/// heat itself, so the omission is not a rounding matter: every cubic metre
/// that freezes or thaws leaked `617.5 MJ` out of the audit.
///
/// This is exactly the term that made the census disagree with
/// `soil_enthalpy_balance.stateAtTemperature`, which conserves
/// `C*(T - Tm) + L*liquid`. Subtracting the two forms leaves
/// `C_liquid*Tm*W_total` with `W_total` the total water, and the pieces that
/// depend on the liquid/ice split cancel only when this offset is present.
/// Measured on Ottawa day one: `(C_ice - C_liquid)*Tm*dV_ice` over the
/// day-one soil ice gain of `1.1168e5 m3` predicts `-6.902e7 MJ`, that is
/// `-7.918e-1 MJ m-2` against a recorded deviation of `-8.128e-1 MJ m-2`,
/// reproducing 97.4 % of it.
///
/// The default exists for the same reason as the latent-heat defaults above:
/// `src/ecosys_ng.zig` is owned by the Integrator lane. See
/// `docs/binding_requests/heat_001_landscape_enthalpy.md`.
pub const default_pure_water_melting_temperature_k: f64 = 273.15;

/// Enthalpy of one cubic metre of frozen water (water-equivalent volume)
/// relative to liquid water at 0 K. Single definition shared by every ice
/// carrier so the five carriers cannot drift apart.
///
/// HEAT-001 third layer. This is `pub` because the *boundary* must use the
/// identical definition, not merely an equivalent one. Snowfall crosses the
/// landscape boundary already frozen, so
/// `landscape_boundary_ledger.accumulateAcceptedPrecipitationHeat` has to
/// credit exactly the enthalpy that `aggregateSnowEnthalpy` will then store
/// for the same water equivalent. Measured on Ottawa day one: the boundary
/// booked the superseded `C_ice*T - L` form while the census stored the
/// ice-branch form, and the `5.753e5 m3` of day-one snowfall times the
/// `617.5 MJ m-3` difference is `+4.075 MJ m-2` of the `+4.468 MJ m-2`
/// deviation. A duplicated expression here would reopen exactly that gap the
/// next time either side changed.
pub fn frozenWaterEnthalpyPerM3(
    temperature_k: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    ice_heat_capacity_megajoules_per_m3_k: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
    pure_water_melting_temperature_k: f64,
) f64 {
    return liquid_water_heat_capacity_megajoules_per_m3_k *
        pure_water_melting_temperature_k -
        latent_heat_of_fusion_megajoules_per_m3 +
        ice_heat_capacity_megajoules_per_m3_k *
            (temperature_k - pure_water_melting_temperature_k);
}

/// Deprecated. Retained only because `src/ecosys_ng.zig` calls this arity
/// from a *nitrogen* diagnostic, where the heat field is discarded, and that
/// file is owned by the Integrator lane. Do not use it for anything that
/// reads `heat_megajoules`.
pub fn aggregateSnow(
    state: *const snow.State,
    ice_density_megagrams_per_m3: f64,
) !Storage {
    return aggregateSnowEnthalpy(
        state,
        ice_density_megagrams_per_m3,
        default_latent_heat_of_fusion_megajoules_per_m3,
        4.19,
        1.9274,
    );
}

pub fn aggregateSnowEnthalpy(
    state: *const snow.State,
    ice_density_megagrams_per_m3: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
    /// HEAT-001 second layer. Needed to re-base the snow owner's combined
    /// sensible capacity onto the shared ice-branch enthalpy definition.
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    ice_heat_capacity_megajoules_per_m3_k: f64,
) !Storage {
    inline for (.{ liquid_water_heat_capacity_megajoules_per_m3_k, ice_heat_capacity_megajoules_per_m3_k }) |value|
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSnowInventoryHeatCapacity;
    if (!std.math.isFinite(ice_density_megagrams_per_m3) or
        ice_density_megagrams_per_m3 <= 0 or
        ice_density_megagrams_per_m3 > 1)
        return error.InvalidSnowIceDensity;
    if (!std.math.isFinite(latent_heat_of_fusion_megajoules_per_m3) or
        latent_heat_of_fusion_megajoules_per_m3 <= 0)
        return error.InvalidLatentHeatOfFusion;
    const layer_count = try std.math.mul(
        usize,
        state.cell_count,
        state.layer_capacity,
    );
    if (state.active.len != layer_count or
        state.solid_snow_water_equivalent_m3.len != layer_count or
        state.liquid_water_volume_m3.len != layer_count or
        state.vapor_water_equivalent_m3.len != layer_count or
        state.ice_volume_m3.len != layer_count or
        state.temperature_k.len != layer_count or
        state.heat_capacity_megajoules_per_k.len != layer_count or
        state.amount_g.len != try std.math.mul(usize, layer_count, snow.species_count))
        return error.SnowInventoryDimensionMismatch;

    var result: Storage = .{};
    diagnostic_snow_solid_water_equivalent_m3 = 0;
    diagnostic_snow_ice_volume_water_equivalent_m3 = 0;
    for (0..layer_count) |layer| {
        const solid = state.solid_snow_water_equivalent_m3[layer];
        const liquid = state.liquid_water_volume_m3[layer];
        const vapor = state.vapor_water_equivalent_m3[layer];
        const ice = state.ice_volume_m3[layer];
        const temperature = state.temperature_k[layer];
        const heat_capacity = state.heat_capacity_megajoules_per_k[layer];
        inline for (.{ solid, liquid, vapor, ice, temperature, heat_capacity }) |value|
            if (!std.math.isFinite(value)) return error.NonFiniteSnowInventory;
        if (solid < 0 or liquid < 0 or vapor < 0 or ice < 0 or
            temperature < 0 or heat_capacity < 0)
            return error.NegativeSnowInventory;

        result.water_m3 += solid + liquid + vapor + ice * ice_density_megagrams_per_m3;
        const frozen_water_equivalent_m3 =
            solid + ice * ice_density_megagrams_per_m3;
        diagnostic_snow_solid_water_equivalent_m3 += solid;
        diagnostic_snow_ice_volume_water_equivalent_m3 += ice * ice_density_megagrams_per_m3;
        // HEAT-001 second layer. The snow owner publishes one combined
        // `heat_capacity_megajoules_per_k`, so the frozen part cannot simply be
        // dropped from it the way the soil and surface aggregators do. Instead
        // the sensible product is taken as published and the frozen carriers
        // are then *re-based* onto the shared ice-branch definition, by
        // removing the `C_i*T` the owner already counted for them and adding
        // the full `C_l*Tm - L + C_i*(T - Tm)` in its place.
        //
        // The net correction per cubic metre is `(C_l - C_i)*Tm`, identical to
        // the soil and surface corrections, which is the property that keeps
        // snow-to-soil and snow-to-surface transfers internal.
        const frozen_enthalpy_correction_megajoules_per_m3 =
            frozenWaterEnthalpyPerM3(
                temperature,
                liquid_water_heat_capacity_megajoules_per_m3_k,
                ice_heat_capacity_megajoules_per_m3_k,
                latent_heat_of_fusion_megajoules_per_m3,
                default_pure_water_melting_temperature_k,
            ) - ice_heat_capacity_megajoules_per_m3_k * temperature;
        result.heat_megajoules += heat_capacity * temperature +
            frozen_enthalpy_correction_megajoules_per_m3 *
                frozen_water_equivalent_m3;
        const amounts = state.amount_g[layer * snow.species_count .. (layer + 1) * snow.species_count];
        for (amounts) |amount| {
            if (!std.math.isFinite(amount)) return error.NonFiniteSnowInventory;
            if (amount < 0) return error.NegativeSnowInventory;
        }
        result.carbon_dioxide_carbon_g +=
            speciesAmount(amounts, .carbon_dioxide_carbon) +
            speciesAmount(amounts, .methane_carbon);
        result.oxygen_g += speciesAmount(amounts, .oxygen);
        result.dinitrogen_nitrogen_g +=
            speciesAmount(amounts, .dinitrogen_nitrogen) +
            speciesAmount(amounts, .nitrous_oxide_nitrogen);
        result.ammonium_nitrogen_g +=
            speciesAmount(amounts, .ammonium_nitrogen) +
            speciesAmount(amounts, .ammonia_nitrogen);
        result.nitrate_nitrogen_g += speciesAmount(amounts, .nitrate_nitrogen);
        result.phosphate_phosphorus_g +=
            speciesAmount(amounts, .hydrogen_phosphate_phosphorus) +
            speciesAmount(amounts, .dihydrogen_phosphate_phosphorus);
        result.ion_inventory_mol +=
            speciesAmount(amounts, .aluminum) / 26.9815385 +
            speciesAmount(amounts, .iron) / 55.845 +
            speciesAmount(amounts, .calcium) / 40.078 +
            speciesAmount(amounts, .magnesium) / 24.305 +
            speciesAmount(amounts, .sodium) / 22.98976928 +
            speciesAmount(amounts, .potassium) / 39.0983 +
            speciesAmount(amounts, .sulfate_sulfur) / 32.065 +
            speciesAmount(amounts, .chloride) / 35.453;
    }
    try result.validate();
    return result;
}

/// REDIST profile physical and gas inventory. The current grid stores ice in
/// water-equivalent m3 already, so no second density conversion is applied.
/// Root gas pools are intentionally not accepted here: their authoritative
/// owner is aggregated separately before the final EXEC reduction.
pub fn aggregateSoilPhysicalAndGas(
    grid: *const grid_module.GridState,
    dry_solid_heat_capacity_megajoules_per_m3_k: []const f64,
    layer_volume_m3: []const f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    ice_heat_capacity_megajoules_per_m3_k: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
    gas_state: *const gas.State,
) !Storage {
    if (grid.layer_count !=
        try std.math.mul(usize, grid.cell_count, grid.soil_layer_capacity) or
        grid.active_soil_layer_count.len != grid.cell_count or
        dry_solid_heat_capacity_megajoules_per_m3_k.len != grid.layer_count or
        layer_volume_m3.len != grid.layer_count or
        gas_state.cell_count != grid.layer_count)
        return error.SoilInventoryDimensionMismatch;
    inline for (.{ liquid_water_heat_capacity_megajoules_per_m3_k, ice_heat_capacity_megajoules_per_m3_k }) |value|
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSoilInventoryHeatCapacity;
    if (!std.math.isFinite(latent_heat_of_fusion_megajoules_per_m3) or
        latent_heat_of_fusion_megajoules_per_m3 <= 0)
        return error.InvalidLatentHeatOfFusion;
    inline for (.{
        grid.matrix_liquid_water_m3.len,
        grid.macropore_liquid_water_m3.len,
        grid.matrix_ice_water_m3.len,
        grid.macropore_ice_water_m3.len,
        grid.water_vapor_volume_m3.len,
        grid.soil_temperature_k.len,
    }) |length| if (length != grid.layer_count)
        return error.SoilInventoryDimensionMismatch;

    var result: Storage = .{};
    diagnostic_soil_matrix_ice_water_equivalent_m3 = 0;
    diagnostic_soil_macropore_ice_water_equivalent_m3 = 0;
    for (0..grid.cell_count) |cell| {
        const active_layers = grid.active_soil_layer_count[cell];
        if (active_layers > grid.soil_layer_capacity)
            return error.InvalidActiveSoilLayerCount;
        for (0..active_layers) |layer| {
            const index = cell * grid.soil_layer_capacity + layer;
            const matrix_liquid = grid.matrix_liquid_water_m3[index];
            const macropore_liquid = grid.macropore_liquid_water_m3[index];
            const matrix_ice_water_equivalent = grid.matrix_ice_water_m3[index];
            const macropore_ice_water_equivalent = grid.macropore_ice_water_m3[index];
            const vapor_water_equivalent = grid.water_vapor_volume_m3[index];
            const temperature = grid.soil_temperature_k[index];
            const dry_solid_heat_capacity = dry_solid_heat_capacity_megajoules_per_m3_k[index];
            const volume = layer_volume_m3[index];
            inline for (.{
                matrix_liquid,
                macropore_liquid,
                matrix_ice_water_equivalent,
                macropore_ice_water_equivalent,
                vapor_water_equivalent,
                temperature,
                dry_solid_heat_capacity,
                volume,
            }) |value| {
                if (!std.math.isFinite(value)) return error.NonFiniteSoilInventory;
                if (value < 0) return error.NegativeSoilInventory;
            }
            result.water_m3 +=
                matrix_liquid +
                macropore_liquid +
                matrix_ice_water_equivalent +
                macropore_ice_water_equivalent +
                vapor_water_equivalent;
            // Reconstruct from authoritative live carriers. Using the cached
            // soil-thermal total here made EXEC heat storage depend on whether
            // its derived table had been refreshed after water and phase
            // transactions.
            // HEAT-001 resolution A: enthalpy, not sensible heat, with the
            // reference state liquid water at 0 K.
            //
            // Frozen water does NOT sit exactly one latent heat below liquid
            // water at the same temperature. Following the liquid branch up to
            // the melting point and the ice branch back down gives
            // `C_l*Tm - L + C_i*(T - Tm)` per cubic metre. The earlier form
            // `C_i*T - L` omitted `(C_l - C_i)*Tm = 617.5 MJ m-3`, which is
            // nearly twice `L`, so freeze/thaw did not cancel in the census:
            // it leaked that offset per cubic metre converted.
            //
            // The sensible term below therefore carries only the *liquid*
            // carriers, and each frozen carrier carries its whole enthalpy.
            // Written this way the census agrees term for term with
            // `soil_enthalpy_balance.stateAtTemperature`, which is the function
            // the soil solver actually conserves.
            const frozen_water_equivalent_m3 =
                matrix_ice_water_equivalent + macropore_ice_water_equivalent;
            const liquid_extensive_heat_capacity_megajoules_per_k =
                dry_solid_heat_capacity * volume +
                liquid_water_heat_capacity_megajoules_per_m3_k *
                    (matrix_liquid + macropore_liquid + vapor_water_equivalent);
            result.heat_megajoules +=
                liquid_extensive_heat_capacity_megajoules_per_k * temperature +
                frozenWaterEnthalpyPerM3(
                    temperature,
                    liquid_water_heat_capacity_megajoules_per_m3_k,
                    ice_heat_capacity_megajoules_per_m3_k,
                    latent_heat_of_fusion_megajoules_per_m3,
                    default_pure_water_melting_temperature_k,
                ) * frozen_water_equivalent_m3;
            diagnostic_soil_matrix_ice_water_equivalent_m3 += matrix_ice_water_equivalent;
            diagnostic_soil_macropore_ice_water_equivalent_m3 += macropore_ice_water_equivalent;

            const first = index * gas.species_count;
            const end = first + gas.species_count;
            const gaseous = gas_state.gaseous_mass_g[first..end];
            const dissolved = gas_state.dissolved_mass_g[first..end];
            const macropore = gas_state.macropore_dissolved_mass_g[first..end];
            const band = gas_state.band_dissolved_mass_g[first..end];
            inline for (0..gas.species_count) |species_index| {
                inline for (.{
                    gaseous[species_index],
                    dissolved[species_index],
                    macropore[species_index],
                    band[species_index],
                }) |value| {
                    if (!std.math.isFinite(value)) return error.NonFiniteSoilInventory;
                    if (value < 0) return error.NegativeSoilInventory;
                }
            }
            result.carbon_dioxide_carbon_g +=
                fourPhaseGas(gaseous, dissolved, macropore, band, .carbon_dioxide) +
                fourPhaseGas(gaseous, dissolved, macropore, band, .methane);
            result.oxygen_g +=
                fourPhaseGas(gaseous, dissolved, macropore, band, .oxygen);
            result.dinitrogen_nitrogen_g +=
                fourPhaseGas(gaseous, dissolved, macropore, band, .nitrogen) +
                fourPhaseGas(gaseous, dissolved, macropore, band, .nitrous_oxide);
            result.ammonium_nitrogen_g +=
                fourPhaseGas(gaseous, dissolved, macropore, band, .ammonia);
        }
    }
    try result.validate();
    return result;
}

/// REDIST `TLCO2P/TLOXYP/TLCH4P/TLN2OP/TLNH3P`: gaseous and aqueous
/// root/mycorrhizal gas storage. Root carbon biomass is a plant owner, but gas
/// already produced into these pools is part of the EXEC soil-system census.
pub fn aggregateRootGas(state: *const plant_roots.State) !Storage {
    const expected = try product(&.{
        state.plant_count,
        plant_roots.biological_domain_count,
        state.soil_layer_count,
    });
    if (expected == 0 or
        state.gaseous_carbon_dioxide_g_c.len != expected or
        state.aqueous_carbon_dioxide_g_c.len != expected or
        state.gaseous_oxygen_g_o.len != expected or
        state.aqueous_oxygen_g_o.len != expected or
        state.gaseous_methane_g_c.len != expected or
        state.aqueous_methane_g_c.len != expected or
        state.gaseous_nitrous_oxide_g_n.len != expected or
        state.aqueous_nitrous_oxide_g_n.len != expected or
        state.gaseous_ammonia_g_n.len != expected or
        state.aqueous_ammonia_g_n.len != expected)
        return error.RootGasInventoryDimensionMismatch;

    var result: Storage = .{};
    for (0..expected) |index| {
        inline for (.{
            state.gaseous_carbon_dioxide_g_c[index],
            state.aqueous_carbon_dioxide_g_c[index],
            state.gaseous_oxygen_g_o[index],
            state.aqueous_oxygen_g_o[index],
            state.gaseous_methane_g_c[index],
            state.aqueous_methane_g_c[index],
            state.gaseous_nitrous_oxide_g_n[index],
            state.aqueous_nitrous_oxide_g_n[index],
            state.gaseous_ammonia_g_n[index],
            state.aqueous_ammonia_g_n[index],
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootGasInventory;
        result.carbon_dioxide_carbon_g +=
            state.gaseous_carbon_dioxide_g_c[index] +
            state.aqueous_carbon_dioxide_g_c[index] +
            state.gaseous_methane_g_c[index] +
            state.aqueous_methane_g_c[index];
        result.oxygen_g += state.gaseous_oxygen_g_o[index] +
            state.aqueous_oxygen_g_o[index];
        result.dinitrogen_nitrogen_g +=
            state.gaseous_nitrous_oxide_g_n[index] +
            state.aqueous_nitrous_oxide_g_n[index];
        result.ammonium_nitrogen_g += state.gaseous_ammonia_g_n[index] +
            state.aqueous_ammonia_g_n[index];
    }
    try result.validate();
    return result;
}

/// REDIST surface `DC/DN/DP + DCC/DNC/DPC` reconstruction. The resulting
/// values all belong to the EXEC residue category (`TLRSD*`), including
/// surface microbial biomass and fire-derived charcoal.
pub fn aggregateSurfaceOrganic(
    state: *const organic.State,
) !Storage {
    const expected_microbial = try product(&.{
        state.layer_count,
        organic.microbial_substrate_count,
        organic.microbial_population_count,
        organic.kinetic_fraction_count,
    });
    const expected_residue = try product(&.{
        state.layer_count,
        organic.substrate_count,
        organic.residue_fraction_count,
    });
    const expected_mobile = try product(&.{
        state.layer_count,
        organic.substrate_count,
    });
    const expected_structural = try product(&.{
        state.layer_count,
        organic.substrate_count,
        organic.structural_fraction_count,
    });
    if (state.microbial.len != expected_microbial or
        state.residue.len != expected_residue or
        state.dissolved.len != expected_mobile or
        state.adsorbed.len != expected_mobile or
        state.dissolved_acetate_carbon_g_c.len != expected_mobile or
        state.adsorbed_acetate_carbon_g_c.len != expected_mobile or
        state.structural.len != expected_structural or
        state.colonized_structural_carbon_g_c.len != expected_structural)
        return error.SurfaceOrganicInventoryDimensionMismatch;

    var result: Storage = .{};
    for (0..state.layer_count) |cell| {
        // OMC: every microbial substrate except humus K=4.
        for (0..organic.microbial_substrate_count) |substrate| {
            if (substrate == 4) continue;
            const first =
                (cell * organic.microbial_substrate_count + substrate) *
                organic.microbial_population_count *
                organic.kinetic_fraction_count;
            const count =
                organic.microbial_population_count *
                organic.kinetic_fraction_count;
            for (state.microbial[first .. first + count]) |pool|
                try addResiduePool(&result, pool);
        }

        // ORC/OQ/OH: source surface balance includes K=0..2.
        for (0..3) |substrate| {
            const residue_first =
                (cell * organic.substrate_count + substrate) *
                organic.residue_fraction_count;
            for (state.residue[residue_first .. residue_first + organic.residue_fraction_count]) |pool| try addResiduePool(&result, pool);
            const mobile = cell * organic.substrate_count + substrate;
            try addResiduePool(&result, state.dissolved[mobile]);
            try addResiduePool(&result, state.adsorbed[mobile]);
            try addResidueCarbon(
                &result,
                state.dissolved_acetate_carbon_g_c[mobile],
            );
            try addResidueCarbon(
                &result,
                state.adsorbed_acetate_carbon_g_c[mobile],
            );
        }

        // OSC M=1..5, K=0..4. Fraction 5 is charcoal but remains TLRSD*.
        const structural_first =
            cell * organic.substrate_count * organic.structural_fraction_count;
        const structural_count =
            organic.substrate_count * organic.structural_fraction_count;
        for (state.structural[structural_first .. structural_first + structural_count]) |pool| try addResiduePool(&result, pool);
    }
    try result.validate();
    return result;
}

/// REDIST mineral-profile `DC/OC` split. Substrate K=4 is humus and enters
/// `TLORG*`; all other substrates enter `TLRSD*`. Unlike the surface block,
/// profile residue/mobile pools span K=0..4. Only runtime-active layers are
/// accepted as authoritative profile storage.
pub fn aggregateSoilOrganic(
    state: *const organic.State,
    grid: *const grid_module.GridState,
) !Storage {
    if (state.layer_count != grid.layer_count or
        grid.layer_count !=
            try std.math.mul(usize, grid.cell_count, grid.soil_layer_capacity) or
        grid.active_soil_layer_count.len != grid.cell_count)
        return error.SoilOrganicInventoryDimensionMismatch;
    try validateOrganicDimensions(state, error.SoilOrganicInventoryDimensionMismatch);

    var result: Storage = .{};
    for (0..grid.cell_count) |cell| {
        const active_layers = grid.active_soil_layer_count[cell];
        if (active_layers > grid.soil_layer_capacity)
            return error.InvalidActiveSoilLayerCount;
        for (0..active_layers) |layer| {
            const layer_cell = cell * grid.soil_layer_capacity + layer;
            for (0..organic.microbial_substrate_count) |substrate| {
                const first =
                    (layer_cell * organic.microbial_substrate_count + substrate) *
                    organic.microbial_population_count *
                    organic.kinetic_fraction_count;
                const count =
                    organic.microbial_population_count *
                    organic.kinetic_fraction_count;
                for (state.microbial[first .. first + count]) |pool| {
                    if (substrate == 4)
                        try addOrganicPool(&result, pool)
                    else
                        try addResiduePool(&result, pool);
                }
            }
            for (0..organic.substrate_count) |substrate| {
                const is_humus = substrate == 4;
                const residue_first =
                    (layer_cell * organic.substrate_count + substrate) *
                    organic.residue_fraction_count;
                for (state.residue[residue_first .. residue_first + organic.residue_fraction_count]) |pool| {
                    if (is_humus)
                        try addOrganicPool(&result, pool)
                    else
                        try addResiduePool(&result, pool);
                }
                const mobile = layer_cell * organic.substrate_count + substrate;
                if (is_humus) {
                    try addOrganicPool(&result, state.dissolved[mobile]);
                    try addOrganicPool(&result, state.adsorbed[mobile]);
                    try addOrganicCarbon(
                        &result,
                        state.dissolved_acetate_carbon_g_c[mobile],
                    );
                    try addOrganicCarbon(
                        &result,
                        state.adsorbed_acetate_carbon_g_c[mobile],
                    );
                } else {
                    try addResiduePool(&result, state.dissolved[mobile]);
                    try addResiduePool(&result, state.adsorbed[mobile]);
                    try addResidueCarbon(
                        &result,
                        state.dissolved_acetate_carbon_g_c[mobile],
                    );
                    try addResidueCarbon(
                        &result,
                        state.adsorbed_acetate_carbon_g_c[mobile],
                    );
                }
                const structural_first =
                    (layer_cell * organic.substrate_count + substrate) *
                    organic.structural_fraction_count;
                for (state.structural[structural_first .. structural_first + organic.structural_fraction_count]) |pool| {
                    if (is_humus)
                        try addOrganicPool(&result, pool)
                    else
                        try addResiduePool(&result, pool);
                }
            }
        }
    }
    try result.validate();
    return result;
}

/// TRNSFR macropore DOC/acetate that persists between transport steps.
/// After `importMicroporeIntoProfile` the micropore amounts are in
/// `profile.dissolved` and counted by `aggregateSoilOrganic`. The macropore
/// amounts remain in the transport state and must be counted here to close the
/// EXEC carbon (and N/P) balance.
pub fn aggregateSoilOrganicTransportMacropore(
    transport: *const organic_transport.State,
    grid: *const grid_module.GridState,
) !Storage {
    const layer_count = try std.math.mul(usize, grid.cell_count, grid.soil_layer_capacity);
    if (transport.layer_count != layer_count or
        grid.active_soil_layer_count.len != grid.cell_count or
        transport.macropore_amount_g.len !=
            try std.math.mul(usize, layer_count, organic_transport.component_count))
        return error.SoilOrganicTransportInventoryDimensionMismatch;

    var result: Storage = .{};
    for (0..grid.cell_count) |cell| {
        const active_layers = grid.active_soil_layer_count[cell];
        if (active_layers > grid.soil_layer_capacity)
            return error.InvalidActiveSoilLayerCount;
        for (0..active_layers) |local_layer| {
            const flat_layer = cell * grid.soil_layer_capacity + local_layer;
            for (0..organic.substrate_count) |substrate| {
                const base =
                    flat_layer * organic_transport.component_count +
                    substrate * organic_transport.components_per_substrate;
                const doc_g_c = transport.macropore_amount_g[base + @intFromEnum(organic_transport.Component.dissolved_organic_carbon)];
                const don_g_n = transport.macropore_amount_g[base + @intFromEnum(organic_transport.Component.dissolved_organic_nitrogen)];
                const dop_g_p = transport.macropore_amount_g[base + @intFromEnum(organic_transport.Component.dissolved_organic_phosphorus)];
                const acetate_g_c = transport.macropore_amount_g[base + @intFromEnum(organic_transport.Component.dissolved_acetate_carbon)];
                inline for (.{ doc_g_c, don_g_n, dop_g_p, acetate_g_c }) |value| {
                    if (!std.math.isFinite(value))
                        return error.NonFiniteSoilOrganicTransportInventory;
                    if (value < 0) return error.NegativeSoilOrganicTransportInventory;
                }
                if (substrate == organic.substrate_count - 1) {
                    result.organic_carbon_g += doc_g_c + acetate_g_c;
                    result.organic_nitrogen_g += don_g_n;
                    result.organic_phosphorus_g += dop_g_p;
                } else {
                    result.residue_carbon_g += doc_g_c + acetate_g_c;
                    result.residue_nitrogen_g += don_g_n;
                    result.residue_phosphorus_g += dop_g_p;
                }
            }
        }
    }
    try result.validate();
    return result;
}

/// REDIST profile mineral-N inventory. Aqueous matrix and macropore amounts
/// are transport-owned extensive mol N; exchangeable ammonium is chemistry-
/// owned mol N/Mg soil; and undissolved fertilizer remains in its runtime
/// inventory until dissolution. Chemistry aqueous concentrations and plant-
/// available nutrient mirrors are intentionally not counted a second time.
pub fn aggregateProfileMineralNitrogen(
    grid: *const grid_module.GridState,
    transport: *const mineral_nitrogen.State,
    chemistry: *const soil_chemistry.State,
    fertilizer: *const nitrogen_fertilizer.State,
    soil_mass_megagrams: []const f64,
    nitrogen_g_per_mol: f64,
) !Storage {
    if (grid.layer_count !=
        try std.math.mul(usize, grid.cell_count, grid.soil_layer_capacity) or
        grid.active_soil_layer_count.len != grid.cell_count or
        transport.cell_count != grid.layer_count or
        chemistry.cell_count != grid.layer_count or
        fertilizer.cell_count != grid.cell_count or
        fertilizer.layer_capacity != grid.soil_layer_capacity or
        fertilizer.soil.len != grid.layer_count or
        soil_mass_megagrams.len != grid.layer_count)
        return error.ProfileMineralNitrogenInventoryDimensionMismatch;
    if (!std.math.isFinite(nitrogen_g_per_mol) or nitrogen_g_per_mol <= 0)
        return error.InvalidNitrogenMolarMass;
    try transport.validate();

    var result: Storage = .{};
    for (0..grid.cell_count) |cell| {
        const active_layers = grid.active_soil_layer_count[cell];
        if (active_layers > grid.soil_layer_capacity)
            return error.InvalidActiveSoilLayerCount;
        for (0..active_layers) |layer| {
            const profile_cell = cell * grid.soil_layer_capacity + layer;
            const mass_megagrams = soil_mass_megagrams[profile_cell];
            if (!std.math.isFinite(mass_megagrams) or mass_megagrams < 0)
                return error.InvalidSoilMass;

            const matrix = try transport.matrix.cellAmountsConst(profile_cell);
            const macropore =
                try transport.macropore.cellAmountsConst(profile_cell);
            const ammonium_mol_n =
                nitrogenAmount(matrix, macropore, .ammonium_non_band) +
                nitrogenAmount(matrix, macropore, .ammonium_band) +
                nitrogenAmount(matrix, macropore, .ammonia_non_band) +
                nitrogenAmount(matrix, macropore, .ammonia_band);
            const nitrate_mol_n =
                nitrogenAmount(matrix, macropore, .nitrate_non_band) +
                nitrogenAmount(matrix, macropore, .nitrate_band) +
                nitrogenAmount(matrix, macropore, .nitrite_non_band) +
                nitrogenAmount(matrix, macropore, .nitrite_band);

            const exchange = chemistry.cation_exchange_mol_per_megagram[profile_cell];
            try validateFiniteNonnegativeStruct(exchange);
            const exchange_ammonium_mol_n =
                (exchange.ammonium_non_band + exchange.ammonium_band) * mass_megagrams;

            const dry = fertilizer.soil[profile_cell];
            try validateFiniteNonnegativeStruct(dry);
            const dry_ammonium_mol_n =
                dry.broadcast_ammonium_mol_n +
                dry.broadcast_ammonia_mol_n +
                dry.broadcast_urea_mol_n +
                dry.banded_ammonium_mol_n +
                dry.banded_ammonia_mol_n +
                dry.banded_urea_mol_n;
            const dry_nitrate_mol_n =
                dry.broadcast_nitrate_mol_n + dry.banded_nitrate_mol_n;

            result.ammonium_nitrogen_g +=
                (ammonium_mol_n + exchange_ammonium_mol_n +
                    dry_ammonium_mol_n) *
                nitrogen_g_per_mol;
            result.nitrate_nitrogen_g +=
                (nitrate_mol_n + dry_nitrate_mol_n) * nitrogen_g_per_mol;

            // Legacy TION atom-count convention: exchange NH4 and dry NH4
            // carry two atoms; dry NH3, urea, and nitrate carry one.
            result.ion_inventory_mol +=
                2 * exchange_ammonium_mol_n +
                2 * (dry.broadcast_ammonium_mol_n +
                    dry.banded_ammonium_mol_n) +
                dry.broadcast_ammonia_mol_n +
                dry.broadcast_urea_mol_n +
                dry.broadcast_nitrate_mol_n +
                dry.banded_ammonia_mol_n +
                dry.banded_urea_mol_n +
                dry.banded_nitrate_mol_n;
        }
    }
    try result.validate();
    return result;
}

fn nitrogenAmount(
    matrix: []const f64,
    macropore: []const f64,
    species: mineral_nitrogen.Species,
) f64 {
    const species_index = @intFromEnum(species);
    return matrix[species_index] + macropore[species_index];
}

fn validateFiniteNonnegativeStruct(value: anytype) !void {
    inline for (std.meta.fields(@TypeOf(value))) |field| {
        const number = @field(value, field.name);
        if (!std.math.isFinite(number) or number < 0)
            return error.InvalidProfileMineralNitrogenState;
    }
}

/// Remaining REDIST profile phosphate and salt inventory. Mobile matrix and
/// macropore solutes are transport-owned extensive mol. Immobile exchange,
/// phosphate surfaces, precipitates, and dry fertilizer are added from their
/// distinct runtime owners. Mineral-N fertilizer and exchangeable ammonium
/// are excluded here because `aggregateProfileMineralNitrogen` owns their
/// TION contribution.
pub fn aggregateProfilePhosphorusAndIons(
    grid: *const grid_module.GridState,
    micropore: *const solute_transport.State,
    macropore: *const solute_transport.State,
    chemistry: *const soil_chemistry.State,
    pending_fertilizer: *const mineral_fertilizer.State,
    soil_water_m3: []const f64,
    soil_mass_megagrams: []const f64,
    fractions: zone_classification.ZoneFractions,
    carbon_g_per_mol: f64,
    phosphorus_g_per_mol: f64,
) !Storage {
    const species_count = solute_species.AqueousSpecies.count;
    if (grid.layer_count !=
        try std.math.mul(usize, grid.cell_count, grid.soil_layer_capacity) or
        grid.active_soil_layer_count.len != grid.cell_count or
        micropore.cell_count != grid.layer_count or
        macropore.cell_count != grid.layer_count or
        micropore.species_count != species_count or
        macropore.species_count != species_count or
        chemistry.cell_count != grid.layer_count or
        pending_fertilizer.cell_count != grid.cell_count or
        pending_fertilizer.layer_capacity != grid.soil_layer_capacity or
        pending_fertilizer.soil.len != grid.layer_count or
        soil_water_m3.len != grid.layer_count or
        soil_mass_megagrams.len != grid.layer_count)
        return error.ProfilePhosphorusIonInventoryDimensionMismatch;
    if (!std.math.isFinite(carbon_g_per_mol) or carbon_g_per_mol <= 0 or
        !std.math.isFinite(phosphorus_g_per_mol) or
        phosphorus_g_per_mol <= 0)
        return error.InvalidElementMolarMass;
    try validatePhosphateFractions(fractions);
    try micropore.validateFinite();
    try macropore.validateFinite();

    var result: Storage = .{};
    for (0..grid.cell_count) |cell| {
        const active_layers = grid.active_soil_layer_count[cell];
        if (active_layers > grid.soil_layer_capacity)
            return error.InvalidActiveSoilLayerCount;
        for (0..active_layers) |layer| {
            const profile_cell = cell * grid.soil_layer_capacity + layer;
            const water_m3 = soil_water_m3[profile_cell];
            const mass_megagrams = soil_mass_megagrams[profile_cell];
            if (!std.math.isFinite(water_m3) or water_m3 < 0 or
                !std.math.isFinite(mass_megagrams) or mass_megagrams < 0)
                return error.InvalidProfileChemistryGeometry;

            const matrix_amounts =
                try micropore.cellAmountsConst(profile_cell);
            const macro_amounts =
                try macropore.cellAmountsConst(profile_cell);
            inline for (@typeInfo(solute_species.AqueousSpecies).@"enum".fields) |field| {
                const species: solute_species.AqueousSpecies =
                    @enumFromInt(field.value);
                const amount_mol =
                    matrix_amounts[field.value] + macro_amounts[field.value];
                result.ion_inventory_mol +=
                    amount_mol * aqueousIonAtomCount(species);
                if (isPhosphateCarrier(species))
                    result.phosphate_phosphorus_g +=
                        amount_mol * phosphorus_g_per_mol;
                if (isCarbonateCarrier(species))
                    result.carbon_dioxide_carbon_g +=
                        amount_mol * carbon_g_per_mol;
            }

            const non_band = chemistry.non_band_phosphate[profile_cell];
            const band = chemistry.band_phosphate[profile_cell];
            try validateFiniteNonnegativeStruct(non_band);
            try validateFiniteNonnegativeStruct(band);
            const immobile_non_band =
                phosphateImmobileInventory(non_band, water_m3, mass_megagrams);
            const immobile_band =
                phosphateImmobileInventory(band, water_m3, mass_megagrams);
            result.phosphate_phosphorus_g += phosphorus_g_per_mol *
                (fractions.phosphate_non_band *
                    immobile_non_band.phosphorus_mol +
                    fractions.phosphate_band * immobile_band.phosphorus_mol);
            result.ion_inventory_mol +=
                fractions.phosphate_non_band * immobile_non_band.ion_mol +
                fractions.phosphate_band * immobile_band.ion_mol;

            const exchange = chemistry.cation_exchange_mol_per_megagram[profile_cell];
            try validateFiniteNonnegativeStruct(exchange);
            const carboxyl_hydrogen =
                chemistry.carboxyl_bound_hydrogen_mol_per_megagram[profile_cell];
            if (!std.math.isFinite(carboxyl_hydrogen) or
                carboxyl_hydrogen < 0)
                return error.InvalidProfileMineralIonState;
            result.ion_inventory_mol += mass_megagrams *
                (exchange.hydrogen + exchange.aluminum + exchange.iron +
                    exchange.calcium + exchange.magnesium + exchange.sodium +
                    exchange.potassium + carboxyl_hydrogen);

            const solids = chemistry.geochemistry_solids[profile_cell];
            try validateFiniteNonnegativeStruct(solids);
            result.ion_inventory_mol +=
                water_m3 * geochemistrySolidIonAtoms(solids);
            result.carbon_dioxide_carbon_g += water_m3 *
                solids.calcite_solid_mol_per_m3 * carbon_g_per_mol;

            const pending = pending_fertilizer.soil[profile_cell];
            try validateFiniteNonnegativeStruct(pending);
            result.phosphate_phosphorus_g += phosphorus_g_per_mol *
                (2 * (pending.broadcast_monocalcium_phosphate_mol +
                    pending.banded_monocalcium_phosphate_mol) +
                    3 * pending.hydroxyapatite_mol);
            result.ion_inventory_mol += pendingMineralIonAtoms(pending);
            result.carbon_dioxide_carbon_g +=
                pending.calcite_mol * carbon_g_per_mol;
        }
    }
    try result.validate();
    return result;
}

const ImmobilePhosphateInventory = struct {
    phosphorus_mol: f64,
    ion_mol: f64,
};

fn phosphateImmobileInventory(
    state: anytype,
    water_m3: f64,
    soil_mass_megagrams: f64,
) ImmobilePhosphateInventory {
    const adsorbed_hpo4 =
        state.adsorbed_hpo4_mol_p_per_megagram * soil_mass_megagrams;
    const adsorbed_h2po4 =
        state.adsorbed_h2po4_mol_p_per_megagram * soil_mass_megagrams;
    // These reaction intermediates are not represented by generic transport
    // carrier slots; they therefore remain authoritative in chemistry.
    const dissolved_hpo4 =
        state.dissolved_hpo4_mol_p_per_m3 * water_m3;
    const dissolved_h2po4 =
        state.dissolved_h2po4_mol_p_per_m3 * water_m3;
    const aluminum_phosphate =
        state.aluminum_phosphate_solid_mol_per_m3 * water_m3;
    const iron_phosphate =
        state.iron_phosphate_solid_mol_per_m3 * water_m3;
    const dicalcium_phosphate =
        state.dicalcium_phosphate_solid_mol_per_m3 * water_m3;
    const hydroxyapatite =
        state.hydroxyapatite_solid_mol_per_m3 * water_m3;
    const monocalcium_phosphate =
        state.monocalcium_phosphate_solid_mol_per_m3 * water_m3;
    return .{
        .phosphorus_mol = dissolved_hpo4 + dissolved_h2po4 +
            adsorbed_hpo4 + adsorbed_h2po4 +
            aluminum_phosphate + iron_phosphate + dicalcium_phosphate +
            3 * hydroxyapatite + 2 * monocalcium_phosphate,
        .ion_mol = 3 * dissolved_hpo4 + 4 * dissolved_h2po4 +
            (state.deprotonated_site_mol_per_megagram +
                2 * state.hydroxyl_site_mol_per_megagram +
                3 * state.protonated_site_mol_per_megagram) *
                soil_mass_megagrams +
            3 * adsorbed_hpo4 + 4 * adsorbed_h2po4 +
            2 * (aluminum_phosphate + iron_phosphate) +
            3 * dicalcium_phosphate +
            9 * hydroxyapatite + 7 * monocalcium_phosphate,
    };
}

fn geochemistrySolidIonAtoms(state: anytype) f64 {
    return 4 * (state.gibbsite_solid_mol_per_m3 +
        state.iron_hydroxide_solid_mol_per_m3) +
        2 * (state.calcite_solid_mol_per_m3 +
            state.gypsum_solid_mol_per_m3) +
        state.aluminum_natural_silicate_mol_per_m3 +
        state.aluminum_ground_silicate_mol_per_m3 +
        state.iron_natural_silicate_mol_per_m3 +
        state.iron_ground_silicate_mol_per_m3 +
        state.calcium_natural_silicate_mol_per_m3 +
        state.calcium_ground_silicate_mol_per_m3 +
        state.magnesium_natural_silicate_mol_per_m3 +
        state.magnesium_ground_silicate_mol_per_m3 +
        state.sodium_natural_silicate_mol_per_m3 +
        state.sodium_ground_silicate_mol_per_m3 +
        state.potassium_natural_silicate_mol_per_m3 +
        state.potassium_ground_silicate_mol_per_m3;
}

fn pendingMineralIonAtoms(state: mineral_fertilizer.Inventory) f64 {
    return 7 * (state.broadcast_monocalcium_phosphate_mol +
        state.banded_monocalcium_phosphate_mol) +
        9 * state.hydroxyapatite_mol +
        2 * (state.calcite_mol + state.gypsum_mol) +
        state.aluminum_ground_silicate_mol +
        state.iron_ground_silicate_mol +
        state.calcium_ground_silicate_mol +
        state.magnesium_ground_silicate_mol +
        state.sodium_ground_silicate_mol +
        state.potassium_ground_silicate_mol;
}

/// Dry mineral fertilizer routed to litter remains outside the concentration
/// chemistry owner until water is available. It is nevertheless part of the
/// whole-landscape EXEC storage and must be counted exactly once.
pub fn aggregatePendingSurfaceMinerals(
    state: *const mineral_fertilizer.State,
    carbon_g_per_mol: f64,
    phosphorus_g_per_mol: f64,
) !Storage {
    if (state.cell_count == 0 or state.surface.len != state.cell_count)
        return error.PendingSurfaceMineralInventoryDimensionMismatch;
    inline for (.{ carbon_g_per_mol, phosphorus_g_per_mol }) |molar_mass|
        if (!std.math.isFinite(molar_mass) or molar_mass <= 0)
            return error.InvalidElementMolarMass;

    var result: Storage = .{};
    for (state.surface) |pending| {
        try validateFiniteNonnegativeStruct(pending);
        result.phosphate_phosphorus_g += phosphorus_g_per_mol *
            (2 * (pending.broadcast_monocalcium_phosphate_mol +
                pending.banded_monocalcium_phosphate_mol) +
                3 * pending.hydroxyapatite_mol);
        result.carbon_dioxide_carbon_g +=
            pending.calcite_mol * carbon_g_per_mol;
        result.ion_inventory_mol += pendingMineralIonAtoms(pending);
    }
    try result.validate();
    return result;
}

fn validatePhosphateFractions(
    fractions: zone_classification.ZoneFractions,
) !void {
    inline for (@typeInfo(zone_classification.ZoneFractions).@"struct".fields) |field| {
        const fraction = @field(fractions, field.name);
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidChemistryZoneFraction;
    }
    if (@abs(fractions.phosphate_non_band + fractions.phosphate_band - 1) >
        1e-12)
        return error.InvalidChemistryZoneFraction;
}

fn isPhosphateCarrier(species: solute_species.AqueousSpecies) bool {
    return solute_species.diffusivityClass(species) == .phosphate;
}

fn isCarbonateCarrier(species: solute_species.AqueousSpecies) bool {
    return switch (species) {
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
}

fn aqueousIonAtomCount(species: solute_species.AqueousSpecies) f64 {
    return switch (species) {
        .aluminum,
        .iron,
        .hydrogen,
        .calcium,
        .magnesium,
        .sodium,
        .potassium,
        .hydroxide,
        .sulfate,
        .chloride,
        .carbonate,
        .non_band_phosphate,
        .band_phosphate,
        => 1,
        .bicarbonate,
        .aluminum_hydroxide_1,
        .aluminum_sulfate,
        .iron_hydroxide_1,
        .iron_sulfate,
        .calcium_hydroxide,
        .calcium_carbonate,
        .calcium_sulfate,
        .magnesium_hydroxide,
        .magnesium_carbonate,
        .magnesium_sulfate,
        .sodium_carbonate,
        .sodium_sulfate,
        .potassium_sulfate,
        .non_band_calcium_phosphate,
        .band_calcium_phosphate,
        => 2,
        .aluminum_hydroxide_2,
        .iron_hydroxide_2,
        .calcium_bicarbonate,
        .magnesium_bicarbonate,
        .non_band_iron_hpo4,
        .non_band_calcium_hpo4,
        .non_band_magnesium_hpo4,
        .band_iron_hpo4,
        .band_calcium_hpo4,
        .band_magnesium_hpo4,
        => 3,
        .aluminum_hydroxide_3,
        .iron_hydroxide_3,
        .hydrogen_silicate,
        .non_band_phosphoric_acid,
        .non_band_iron_h2po4,
        .non_band_calcium_h2po4,
        .band_phosphoric_acid,
        .band_iron_h2po4,
        .band_calcium_h2po4,
        => 4,
        .aluminum_hydroxide_4,
        .iron_hydroxide_4,
        => 5,
    };
}

/// REDIST surface mineral N/P and `SST=SSS+SSF+SSX+SSP` inventory expressed
/// through the modern free-ion, exchange-site, fertilizer, and mineral owners.
/// Concentrations are made extensive with litter water; exchange values use
/// dry litter mass. Molecular coefficients retain the source's atom-count
/// convention for TION rather than charge equivalents.
pub fn aggregateSurfaceChemistry(
    chemistry: *const litter_chemistry.State,
    fertilizer: *const litter_fertilizer.State,
    denitrification_nitrite_g_n: []const f64,
    litter_water_m3: []const f64,
    litter_dry_mass_megagrams: []const f64,
    carbon_g_per_mol: f64,
    nitrogen_g_per_mol: f64,
    phosphorus_g_per_mol: f64,
) !Storage {
    const cells = chemistry.cells.len;
    if (cells == 0 or fertilizer.cells.len != cells or
        fertilizer.formulation.len != cells or litter_water_m3.len != cells or
        litter_dry_mass_megagrams.len != cells or
        denitrification_nitrite_g_n.len != cells)
        return error.SurfaceChemistryInventoryDimensionMismatch;
    inline for (.{ carbon_g_per_mol, nitrogen_g_per_mol, phosphorus_g_per_mol }) |molar_mass|
        if (!std.math.isFinite(molar_mass) or molar_mass <= 0)
            return error.InvalidInventoryMolarMass;

    var result: Storage = .{};
    for (0..cells) |cell_index| {
        const water = litter_water_m3[cell_index];
        const dry_mass = litter_dry_mass_megagrams[cell_index];
        const nitrite_g_n = denitrification_nitrite_g_n[cell_index];
        if (!std.math.isFinite(water) or !std.math.isFinite(dry_mass))
            return error.NonFiniteSurfaceChemistryInventory;
        if (!std.math.isFinite(nitrite_g_n))
            return error.NonFiniteSurfaceChemistryInventory;
        if (water < 0 or dry_mass < 0 or nitrite_g_n < 0)
            return error.NegativeSurfaceChemistryInventory;
        const cell = chemistry.cells[cell_index];
        try validateNumericStruct(cell);
        const solid_fertilizer = fertilizer.cells[cell_index];
        inline for (std.meta.fields(litter_fertilizer.Inventory)) |field| {
            const value = @field(solid_fertilizer, field.name);
            if (!std.math.isFinite(value))
                return error.NonFiniteSurfaceChemistryInventory;
            if (value < 0) return error.NegativeSurfaceChemistryInventory;
        }

        result.ammonium_nitrogen_g += nitrogen_g_per_mol * (water * (cell.ammonium_mol_per_m3 + cell.ammonia_mol_per_m3) +
            dry_mass * cell.exchange.ammonium_mol_per_megagram +
            solid_fertilizer.ammonium_mol_n +
            solid_fertilizer.ammonia_mol_n +
            solid_fertilizer.urea_mol_n);
        result.nitrate_nitrogen_g += nitrogen_g_per_mol * (water * cell.nitrate_mol_per_m3 +
            solid_fertilizer.nitrate_mol_n) + nitrite_g_n;
        result.carbon_dioxide_carbon_g += carbon_g_per_mol * water *
            (cell.carbonate_mol_per_m3 + cell.bicarbonate_mol_per_m3 +
                cell.salt_minerals.calcite_mol_per_m3);
        result.phosphate_phosphorus_g += phosphorus_g_per_mol * (water * (cell.hpo4_mol_p_per_m3 +
            cell.h2po4_mol_p_per_m3 +
            cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 +
            cell.phosphate_minerals.iron_phosphate_mol_per_m3 +
            cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3 +
            2 * cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 +
            3 * cell.phosphate_minerals.hydroxyapatite_mol_per_m3) +
            dry_mass * (cell.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram +
                cell.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram));

        const dissolved_ion_atoms_mol_per_m3 =
            cell.aluminum_mol_per_m3 +
            cell.iron_mol_per_m3 +
            cell.hydrogen_mol_per_m3 +
            cell.calcium_mol_per_m3 +
            cell.magnesium_mol_per_m3 +
            cell.sodium_mol_per_m3 +
            cell.potassium_mol_per_m3 +
            cell.hydroxide_mol_per_m3 +
            cell.sulfate_mol_per_m3 +
            cell.chloride_mol_per_m3 +
            cell.carbonate_mol_per_m3 +
            2 * cell.bicarbonate_mol_per_m3 +
            3 * cell.hpo4_mol_p_per_m3 +
            4 * cell.h2po4_mol_p_per_m3;
        const exchange_ion_atoms_mol_per_megagram =
            cell.exchange.hydrogen_mol_per_megagram +
            cell.exchange.aluminum_mol_per_megagram +
            cell.exchange.iron_mol_per_megagram +
            cell.exchange.calcium_mol_per_megagram +
            cell.exchange.magnesium_mol_per_megagram +
            cell.exchange.sodium_mol_per_megagram +
            cell.exchange.potassium_mol_per_megagram +
            2 * cell.exchange.ammonium_mol_per_megagram +
            cell.phosphate_surface.deprotonated_site_mol_per_megagram +
            2 * cell.phosphate_surface.hydroxyl_site_mol_per_megagram +
            3 * cell.phosphate_surface.protonated_site_mol_per_megagram +
            3 * cell.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram +
            4 * cell.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram;
        const mineral_ion_atoms_mol_per_m3 =
            2 * (cell.salt_minerals.calcite_mol_per_m3 +
                cell.salt_minerals.gypsum_mol_per_m3 +
                cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 +
                cell.phosphate_minerals.iron_phosphate_mol_per_m3) +
            3 * cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3 +
            4 * (cell.salt_minerals.gibbsite_mol_per_m3 +
                cell.salt_minerals.iron_hydroxide_mol_per_m3) +
            7 * cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 +
            9 * cell.phosphate_minerals.hydroxyapatite_mol_per_m3;
        const fertilizer_ion_atoms_mol =
            2 * solid_fertilizer.ammonium_mol_n +
            solid_fertilizer.ammonia_mol_n +
            solid_fertilizer.urea_mol_n +
            solid_fertilizer.nitrate_mol_n;
        result.ion_inventory_mol +=
            water * (dissolved_ion_atoms_mol_per_m3 + mineral_ion_atoms_mol_per_m3) +
            dry_mass * exchange_ion_atoms_mol_per_megagram +
            fertilizer_ion_atoms_mol;
    }
    try result.validate();
    return result;
}

pub const SurfacePhysicalParameters = struct {
    dry_organic_heat_capacity_megajoules_per_g_c_k: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    ice_heat_capacity_megajoules_per_m3_k: f64,
    /// HEAT-001 resolution A. See
    /// `docs/binding_requests/heat_001_landscape_enthalpy.md`: the default is
    /// a temporary bridge because `src/ecosys_ng.zig` builds this struct and
    /// is owned by the Integrator lane. It equals the value every shipped
    /// runscript carries and must be deleted once the binding request lands.
    latent_heat_of_fusion_megajoules_per_m3: f64 = 333,
    water_molar_mass_g_per_mol: f64,
    liquid_water_density_g_per_m3: f64,
};

/// REDIST surface `WSS`, `TENGYC`, and litter gas inventory. The modern ice
/// carrier is explicitly water-equivalent m3, so it is added directly to
/// water storage. Water vapor is held as mol by the gas owner and converted
/// to its liquid-water-equivalent volume for both storage and heat capacity.
pub fn aggregateSurfacePhysicalAndGas(
    surface: *const surface_precipitation.RuntimeState,
    surface_ice_water_equivalent_m3: []const f64,
    grid: *const grid_module.GridState,
    gas_state: *const gas.State,
    surface_organic: *const organic.State,
    parameters: SurfacePhysicalParameters,
) !Storage {
    const cells = surface.cell_count;
    if (cells == 0 or surface.litter_water_m3.len != cells or
        surface_ice_water_equivalent_m3.len != cells or
        grid.cell_count != cells or grid.surface_temperature_k.len != cells or
        gas_state.cell_count != cells or surface_organic.layer_count != cells)
        return error.SurfacePhysicalInventoryDimensionMismatch;
    inline for (std.meta.fields(SurfacePhysicalParameters)) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSurfacePhysicalInventoryParameter;
    }

    var result: Storage = .{};
    diagnostic_surface_ice_water_equivalent_m3 = 0;
    for (0..cells) |cell| {
        const liquid = surface.litter_water_m3[cell];
        const ice_water_equivalent = surface_ice_water_equivalent_m3[cell];
        const temperature = grid.surface_temperature_k[cell];
        const vapor_mol = gas_state.water_vapor_mol[cell];
        inline for (.{ liquid, ice_water_equivalent, temperature, vapor_mol }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteSurfacePhysicalInventory;
            if (value < 0) return error.NegativeSurfacePhysicalInventory;
        }
        if (temperature <= 0) return error.InvalidSurfaceTemperature;
        const vapor_water_equivalent_m3 =
            vapor_mol * parameters.water_molar_mass_g_per_mol /
            parameters.liquid_water_density_g_per_m3;
        const organic_carbon_g_c = try surface_organic.totalCarbon_g_c(cell);
        const liquid_heat_capacity_megajoules_per_k =
            parameters.dry_organic_heat_capacity_megajoules_per_g_c_k *
            organic_carbon_g_c +
            parameters.liquid_water_heat_capacity_megajoules_per_m3_k *
                (liquid + vapor_water_equivalent_m3);
        result.water_m3 +=
            liquid + vapor_water_equivalent_m3 + ice_water_equivalent;
        // HEAT-001 resolution A. `surface_ice_water_equivalent_m3` is the one
        // authoritative carrier for both surface litter ice and pond ice, so a
        // single latent term covers both. Its latent contribution makes the
        // surface solver's freeze/thaw repartition internal to the census.
        //
        // HEAT-001 second layer: the frozen enthalpy is
        // `C_l*Tm - L + C_i*(T - Tm)`, not `C_i*T - L`. Same correction and
        // same reasoning as the soil carriers above; see
        // `frozenWaterEnthalpyPerM3`. Applying it to only some carriers would
        // make the pond_domain_transaction surface-to-soil ice transfer stop
        // cancelling, and that transfer is measured at `1.34e4 m3` on Ottawa
        // day one, so all carriers must use the one definition.
        result.heat_megajoules +=
            liquid_heat_capacity_megajoules_per_k * temperature +
            frozenWaterEnthalpyPerM3(
                temperature,
                parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
                parameters.ice_heat_capacity_megajoules_per_m3_k,
                parameters.latent_heat_of_fusion_megajoules_per_m3,
                default_pure_water_melting_temperature_k,
            ) * ice_water_equivalent;
        diagnostic_surface_ice_water_equivalent_m3 += ice_water_equivalent;

        const first = cell * gas.species_count;
        const end = first + gas.species_count;
        const gaseous = gas_state.gaseous_mass_g[first..end];
        const dissolved = gas_state.dissolved_mass_g[first..end];
        const macropore = gas_state.macropore_dissolved_mass_g[first..end];
        const band = gas_state.band_dissolved_mass_g[first..end];
        inline for (0..gas.species_count) |species_index| {
            inline for (.{
                gaseous[species_index],
                dissolved[species_index],
                macropore[species_index],
                band[species_index],
            }) |value| {
                if (!std.math.isFinite(value))
                    return error.NonFiniteSurfacePhysicalInventory;
                if (value < 0) return error.NegativeSurfacePhysicalInventory;
            }
        }
        result.carbon_dioxide_carbon_g +=
            fourPhaseGas(gaseous, dissolved, macropore, band, .carbon_dioxide) +
            fourPhaseGas(gaseous, dissolved, macropore, band, .methane);
        result.oxygen_g +=
            fourPhaseGas(gaseous, dissolved, macropore, band, .oxygen);
        result.dinitrogen_nitrogen_g +=
            fourPhaseGas(gaseous, dissolved, macropore, band, .nitrogen) +
            fourPhaseGas(gaseous, dissolved, macropore, band, .nitrous_oxide);
        result.ammonium_nitrogen_g +=
            fourPhaseGas(gaseous, dissolved, macropore, band, .ammonia);
    }
    try result.validate();
    return result;
}

/// EXTRACT/REDIST `TVOLWP + TVOLWC` and `TENGYC`. Internal canopy water is
/// stored as depth per cell area; intercepted living/dead water is already
/// extensive. `previous_water_energy_megajoules` is the authoritative ENGYX carrier
/// after the accepted canopy surface transaction.
pub fn aggregateCanopyWaterAndHeat(
    plants: *const grid_module.PlantState,
    retention: *const canopy_retention.State,
    cell_area_m2: []const f64,
) !Storage {
    if (plants.cell_count == 0 or plants.species_count == 0 or
        retention.cell_count != plants.cell_count or
        retention.species_count != plants.species_count or
        cell_area_m2.len != plants.cell_count)
        return error.CanopyInventoryDimensionMismatch;
    const plant_count = try std.math.mul(
        usize,
        plants.cell_count,
        plants.species_count,
    );
    inline for (.{
        plants.canopy_water_storage_m_per_m2.len,
        retention.living_surface_water_m3.len,
        retention.standing_dead_surface_water_m3.len,
        retention.previous_water_energy_megajoules.len,
    }) |length| if (length != plant_count)
        return error.CanopyInventoryDimensionMismatch;

    var result: Storage = .{};
    for (0..plants.cell_count) |cell| {
        const area = cell_area_m2[cell];
        if (!std.math.isFinite(area) or area <= 0)
            return error.InvalidCanopyInventoryCellArea;
        for (0..plants.species_count) |species| {
            const plant = cell * plants.species_count + species;
            const internal_water_depth_m =
                plants.canopy_water_storage_m_per_m2[plant];
            const living_surface_water =
                retention.living_surface_water_m3[plant];
            const dead_surface_water =
                retention.standing_dead_surface_water_m3[plant];
            const water_energy = retention.previous_water_energy_megajoules[plant];
            inline for (.{
                internal_water_depth_m,
                living_surface_water,
                dead_surface_water,
                water_energy,
            }) |value| {
                if (!std.math.isFinite(value))
                    return error.NonFiniteCanopyInventory;
                if (value < 0) return error.NegativeCanopyInventory;
            }
            result.water_m3 +=
                internal_water_depth_m * area +
                living_surface_water +
                dead_surface_water;
            result.heat_megajoules += water_energy;
        }
    }
    try result.validate();
    return result;
}

fn validateNumericStruct(value: anytype) !void {
    const T = @TypeOf(value);
    inline for (std.meta.fields(T)) |field| switch (@typeInfo(field.type)) {
        .float => {
            const number = @field(value, field.name);
            if (!std.math.isFinite(number))
                return error.NonFiniteSurfaceChemistryInventory;
            if (number < 0) return error.NegativeSurfaceChemistryInventory;
        },
        .@"struct" => try validateNumericStruct(@field(value, field.name)),
        else => @compileError("surface chemistry inventory must be numeric"),
    };
}

fn validateOrganicDimensions(state: *const organic.State, failure: anyerror) !void {
    const expected_microbial = try product(&.{
        state.layer_count,
        organic.microbial_substrate_count,
        organic.microbial_population_count,
        organic.kinetic_fraction_count,
    });
    const expected_residue = try product(&.{
        state.layer_count,
        organic.substrate_count,
        organic.residue_fraction_count,
    });
    const expected_mobile = try product(&.{
        state.layer_count,
        organic.substrate_count,
    });
    const expected_structural = try product(&.{
        state.layer_count,
        organic.substrate_count,
        organic.structural_fraction_count,
    });
    if (state.microbial.len != expected_microbial or
        state.residue.len != expected_residue or
        state.dissolved.len != expected_mobile or
        state.adsorbed.len != expected_mobile or
        state.dissolved_acetate_carbon_g_c.len != expected_mobile or
        state.adsorbed_acetate_carbon_g_c.len != expected_mobile or
        state.structural.len != expected_structural or
        state.colonized_structural_carbon_g_c.len != expected_structural)
        return failure;
}

fn addResiduePool(
    result: *Storage,
    pool: organic.ElementPool,
) !void {
    inline for (.{
        pool.carbon_g_c,
        pool.nitrogen_g_n,
        pool.phosphorus_g_p,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOrganicInventory;
        if (value < 0) return error.NegativeSurfaceOrganicInventory;
    }
    result.residue_carbon_g += pool.carbon_g_c;
    result.residue_nitrogen_g += pool.nitrogen_g_n;
    result.residue_phosphorus_g += pool.phosphorus_g_p;
}

fn addOrganicPool(result: *Storage, pool: organic.ElementPool) !void {
    inline for (.{
        pool.carbon_g_c,
        pool.nitrogen_g_n,
        pool.phosphorus_g_p,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteSoilOrganicInventory;
        if (value < 0) return error.NegativeSoilOrganicInventory;
    }
    result.organic_carbon_g += pool.carbon_g_c;
    result.organic_nitrogen_g += pool.nitrogen_g_n;
    result.organic_phosphorus_g += pool.phosphorus_g_p;
}

fn addResidueCarbon(result: *Storage, carbon_g_c: f64) !void {
    if (!std.math.isFinite(carbon_g_c))
        return error.NonFiniteSurfaceOrganicInventory;
    if (carbon_g_c < 0) return error.NegativeSurfaceOrganicInventory;
    result.residue_carbon_g += carbon_g_c;
}

fn addOrganicCarbon(result: *Storage, carbon_g_c: f64) !void {
    if (!std.math.isFinite(carbon_g_c))
        return error.NonFiniteSoilOrganicInventory;
    if (carbon_g_c < 0) return error.NegativeSoilOrganicInventory;
    result.organic_carbon_g += carbon_g_c;
}

fn product(values: []const usize) !usize {
    var result: usize = 1;
    for (values) |value| result = try std.math.mul(usize, result, value);
    return result;
}

fn fourPhaseGas(
    gaseous: []const f64,
    dissolved: []const f64,
    macropore: []const f64,
    band: []const f64,
    species: gas.Species,
) f64 {
    const index = @intFromEnum(species);
    return gaseous[index] + dissolved[index] + macropore[index] + band[index];
}

fn speciesAmount(amounts: []const f64, species: snow.Species) f64 {
    return amounts[@intFromEnum(species)];
}

test "REDIST snow inventory sums every runtime cell and layer" {
    var state = try snow.State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    for (0..6) |layer| {
        state.active[layer] = layer != 5;
        state.solid_snow_water_equivalent_m3[layer] = 1;
        state.liquid_water_volume_m3[layer] = 2;
        state.vapor_water_equivalent_m3[layer] = 3;
        state.ice_volume_m3[layer] = 4;
        state.temperature_k[layer] = 250;
        state.heat_capacity_megajoules_per_k[layer] = 0.5;
        const values = try state.amounts(layer / 3, layer % 3);
        values[@intFromEnum(snow.Species.carbon_dioxide_carbon)] = 1;
        values[@intFromEnum(snow.Species.methane_carbon)] = 2;
        values[@intFromEnum(snow.Species.oxygen)] = 3;
        values[@intFromEnum(snow.Species.dinitrogen_nitrogen)] = 4;
        values[@intFromEnum(snow.Species.nitrous_oxide_nitrogen)] = 5;
        values[@intFromEnum(snow.Species.ammonium_nitrogen)] = 6;
        values[@intFromEnum(snow.Species.ammonia_nitrogen)] = 7;
        values[@intFromEnum(snow.Species.nitrate_nitrogen)] = 8;
        values[@intFromEnum(snow.Species.hydrogen_phosphate_phosphorus)] = 9;
        values[@intFromEnum(snow.Species.dihydrogen_phosphate_phosphorus)] = 10;
        values[@intFromEnum(snow.Species.calcium)] = 40.078;
    }
    const inventory = try aggregateSnow(&state, 0.92);
    try std.testing.expectApproxEqAbs(6 * (1 + 2 + 3 + 4 * 0.92), inventory.water_m3, 1e-12);
    // HEAT-001 resolution A, second layer: enthalpy relative to LIQUID water at
    // 0 K. The frozen carriers here are the solid snow water equivalent and the
    // ice volume in water equivalent. Each carries
    // `C_l*Tm - L + C_i*(T - Tm)` per cubic metre, not `C_i*T - L`.
    //
    // The owner publishes one combined capacity (0.5 MJ/K per layer at 250 K),
    // so the aggregator keeps that sensible product and re-bases the frozen
    // part by `(C_l - C_i)*Tm - L` per cubic metre. That is the whole content
    // of the fix: this expectation moved from `-8600.64` to `+8753.6196552`,
    // and the difference is exactly `(4.19 - 1.9274)*273.15*28.08 = 17354.26`,
    // the offset the previous form omitted.
    const frozen_water_equivalent_m3 = 6 * (1.0 + 4 * 0.92);
    try std.testing.expectApproxEqAbs(
        6 * (0.5 * 250) +
            (4.19 * 273.15 - default_latent_heat_of_fusion_megajoules_per_m3 +
                1.9274 * (250 - 273.15) - 1.9274 * 250) * frozen_water_equivalent_m3,
        inventory.heat_megajoules,
        1e-9,
    );
    try std.testing.expectEqual(@as(f64, 18), inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(@as(f64, 18), inventory.oxygen_g);
    try std.testing.expectEqual(@as(f64, 54), inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 78), inventory.ammonium_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 48), inventory.nitrate_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 114), inventory.phosphate_phosphorus_g);
    try std.testing.expectApproxEqAbs(6, inventory.ion_inventory_mol, 1e-12);
}

test "snow inventory rejects hidden non-finite inactive-layer mass" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.amount_g[snow.species_count] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSnowInventory,
        aggregateSnow(&state, 0.92),
    );
}

test "REDIST soil inventory uses runtime active layers and all gas phases" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 3,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 2;
    grid.active_soil_layer_count[1] = 1;
    for (0..grid.layer_count) |index| {
        grid.matrix_liquid_water_m3[index] = 1;
        grid.macropore_liquid_water_m3[index] = 2;
        grid.matrix_ice_water_m3[index] = 3;
        grid.macropore_ice_water_m3[index] = 4;
        grid.water_vapor_volume_m3[index] = 5;
        grid.soil_temperature_k[index] = 250;
    }
    var gas_state = try gas.State.init(std.testing.allocator, grid.layer_count);
    defer gas_state.deinit();
    for (0..grid.layer_count) |layer| {
        const first = layer * gas.species_count;
        inline for (.{ gas_state.gaseous_mass_g, gas_state.dissolved_mass_g, gas_state.macropore_dissolved_mass_g, gas_state.band_dissolved_mass_g }) |phase| {
            phase[first + @intFromEnum(gas.Species.carbon_dioxide)] = 1;
            phase[first + @intFromEnum(gas.Species.methane)] = 2;
            phase[first + @intFromEnum(gas.Species.oxygen)] = 3;
            phase[first + @intFromEnum(gas.Species.nitrogen)] = 4;
            phase[first + @intFromEnum(gas.Species.nitrous_oxide)] = 5;
            phase[first + @intFromEnum(gas.Species.ammonia)] = 6;
        }
    }
    const dry_heat_capacity = [_]f64{2} ** 6;
    const volume = [_]f64{0.5} ** 6;
    const inventory = try aggregateSoilPhysicalAndGas(
        &grid,
        &dry_heat_capacity,
        &volume,
        4,
        1.5,
        333,
        &gas_state,
    );
    const active: f64 = 3;
    try std.testing.expectEqual(active * 15, inventory.water_m3);
    // HEAT-001 resolution A, second layer: enthalpy relative to LIQUID water at
    // 0 K. Frozen carriers do not sit exactly `L` below liquid at the same
    // temperature; they follow the liquid branch to the melting point and the
    // ice branch back down, giving `C_l*Tm - L + C_i*(T - Tm)` per cubic metre.
    //
    // Note this test's heat capacities are the deliberately non-physical
    // `C_l = 4`, `C_i = 1.5`, which is useful here: the correction
    // `(C_l - C_i)*Tm` is then `2.5*273.15 = 682.875` per cubic metre rather
    // than the production `617.5`, so the expectation would not accidentally
    // agree if the code hardcoded production constants instead of using the
    // passed-in ones.
    const frozen_water_equivalent_m3: f64 = 3 + 4;
    try std.testing.expectEqual(
        active * ((2 * 0.5 + 4 * (1 + 2 + 5)) * 250 +
            (4 * 273.15 - 333 + 1.5 * (250 - 273.15)) * frozen_water_equivalent_m3),
        inventory.heat_megajoules,
    );
    try std.testing.expectEqual(active * 4 * 3, inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(active * 4 * 3, inventory.oxygen_g);
    try std.testing.expectEqual(active * 4 * 9, inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(active * 4 * 6, inventory.ammonium_nitrogen_g);
}

test "REDIST root gas inventory includes both root phases" {
    var roots = try plant_roots.State.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.gaseous_carbon_dioxide_g_c[0] = 1;
    roots.aqueous_carbon_dioxide_g_c[0] = 2;
    roots.gaseous_methane_g_c[0] = 3;
    roots.aqueous_methane_g_c[0] = 4;
    roots.gaseous_oxygen_g_o[0] = 5;
    roots.aqueous_oxygen_g_o[0] = 6;
    roots.gaseous_nitrous_oxide_g_n[0] = 7;
    roots.aqueous_nitrous_oxide_g_n[0] = 8;
    roots.gaseous_ammonia_g_n[0] = 9;
    roots.aqueous_ammonia_g_n[0] = 10;

    const inventory = try aggregateRootGas(&roots);
    try std.testing.expectEqual(@as(f64, 10), inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(@as(f64, 11), inventory.oxygen_g);
    try std.testing.expectEqual(@as(f64, 15), inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 19), inventory.ammonium_nitrogen_g);
}

test "REDIST surface organic category split excludes humus microbial pool" {
    var state = try organic.State.init(std.testing.allocator, 2);
    defer state.deinit();
    // One counted microbial pool and one excluded K=4 humus microbial pool.
    state.microbial[0] = .{
        .carbon_g_c = 1,
        .nitrogen_g_n = 2,
        .phosphorus_g_p = 3,
    };
    const humus_microbial =
        4 * organic.microbial_population_count * organic.kinetic_fraction_count;
    state.microbial[humus_microbial] = .{
        .carbon_g_c = 100,
        .nitrogen_g_n = 100,
        .phosphorus_g_p = 100,
    };
    state.residue[0] = .{
        .carbon_g_c = 4,
        .nitrogen_g_n = 5,
        .phosphorus_g_p = 6,
    };
    state.dissolved[0] = .{
        .carbon_g_c = 7,
        .nitrogen_g_n = 8,
        .phosphorus_g_p = 9,
    };
    state.adsorbed[0] = .{
        .carbon_g_c = 10,
        .nitrogen_g_n = 11,
        .phosphorus_g_p = 12,
    };
    state.dissolved_acetate_carbon_g_c[0] = 13;
    state.adsorbed_acetate_carbon_g_c[0] = 14;
    // Fifth structural fraction is persistent charcoal and must be included.
    state.structural[organic.structural_fraction_count - 1] = .{
        .carbon_g_c = 15,
        .nitrogen_g_n = 16,
        .phosphorus_g_p = 17,
    };
    state.colonized_structural_carbon_g_c[
        organic.structural_fraction_count - 1
    ] = 9; // subset diagnostic; must not be counted a second time.
    const inventory = try aggregateSurfaceOrganic(&state);
    try std.testing.expectEqual(@as(f64, 64), inventory.residue_carbon_g);
    try std.testing.expectEqual(@as(f64, 42), inventory.residue_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 47), inventory.residue_phosphorus_g);
    try std.testing.expectEqual(@as(f64, 0), inventory.organic_carbon_g);
}

test "REDIST soil organic split assigns only K=4 to humus" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 1,
            .lat_count = 1,
            .soil_layers = 2,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 1;
    var state = try organic.State.init(std.testing.allocator, grid.layer_count);
    defer state.deinit();

    state.microbial[0] = .{
        .carbon_g_c = 1,
        .nitrogen_g_n = 2,
        .phosphorus_g_p = 3,
    };
    const humus_microbial =
        4 * organic.microbial_population_count * organic.kinetic_fraction_count;
    state.microbial[humus_microbial] = .{
        .carbon_g_c = 4,
        .nitrogen_g_n = 5,
        .phosphorus_g_p = 6,
    };
    const humus_mobile = 4;
    state.dissolved[humus_mobile] = .{
        .carbon_g_c = 7,
        .nitrogen_g_n = 8,
        .phosphorus_g_p = 9,
    };
    state.dissolved_acetate_carbon_g_c[humus_mobile] = 10;
    state.structural[organic.structural_fraction_count - 1] = .{
        .carbon_g_c = 11,
        .nitrogen_g_n = 12,
        .phosphorus_g_p = 13,
    };
    const humus_charcoal =
        4 * organic.structural_fraction_count +
        organic.structural_fraction_count - 1;
    state.structural[humus_charcoal] = .{
        .carbon_g_c = 14,
        .nitrogen_g_n = 15,
        .phosphorus_g_p = 16,
    };
    // Inactive capacity must not enter the profile inventory.
    const inactive_first =
        organic.microbial_substrate_count *
        organic.microbial_population_count *
        organic.kinetic_fraction_count;
    state.microbial[inactive_first].carbon_g_c = 1000;

    const inventory = try aggregateSoilOrganic(&state, &grid);
    try std.testing.expectEqual(@as(f64, 12), inventory.residue_carbon_g);
    try std.testing.expectEqual(@as(f64, 14), inventory.residue_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 16), inventory.residue_phosphorus_g);
    try std.testing.expectEqual(@as(f64, 35), inventory.organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 28), inventory.organic_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 31), inventory.organic_phosphorus_g);
}

test "REDIST surface chemistry retains N P and TION stoichiometry" {
    var chemistry = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    var fertilizer = try litter_fertilizer.State.init(std.testing.allocator, 1);
    defer fertilizer.deinit();
    const cell = &chemistry.cells[0];
    cell.ammonium_mol_per_m3 = 1;
    cell.ammonia_mol_per_m3 = 2;
    cell.nitrate_mol_per_m3 = 3;
    cell.hpo4_mol_p_per_m3 = 4;
    cell.h2po4_mol_p_per_m3 = 5;
    cell.calcium_mol_per_m3 = 6;
    cell.bicarbonate_mol_per_m3 = 7;
    cell.exchange.ammonium_mol_per_megagram = 8;
    cell.exchange.calcium_mol_per_megagram = 9;
    cell.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram = 10;
    cell.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram = 11;
    cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 = 12;
    cell.phosphate_minerals.hydroxyapatite_mol_per_m3 = 13;
    cell.salt_minerals.gypsum_mol_per_m3 = 14;
    fertilizer.cells[0].ammonium_mol_n = 15;
    fertilizer.cells[0].ammonia_mol_n = 16;
    fertilizer.cells[0].urea_mol_n = 17;
    fertilizer.cells[0].nitrate_mol_n = 18;

    const inventory = try aggregateSurfaceChemistry(
        &chemistry,
        &fertilizer,
        &.{19},
        &.{2},
        &.{3},
        12,
        14,
        31,
    );
    try std.testing.expectEqual(
        @as(f64, 14 * (2 * (1 + 2) + 3 * 8 + 15 + 16 + 17)),
        inventory.ammonium_nitrogen_g,
    );
    try std.testing.expectEqual(
        @as(f64, 14 * (2 * 3 + 18) + 19),
        inventory.nitrate_nitrogen_g,
    );
    try std.testing.expectEqual(
        @as(f64, 31 * (2 * (4 + 5 + 2 * 12 + 3 * 13) + 3 * (10 + 11))),
        inventory.phosphate_phosphorus_g,
    );
    try std.testing.expectEqual(
        @as(f64, 12 * 2 * 7),
        inventory.carbon_dioxide_carbon_g,
    );
    const expected_ions =
        2 * ((6 + 2 * 7 + 3 * 4 + 4 * 5) +
            (2 * 14 + 7 * 12 + 9 * 13)) +
        3 * (2 * 8 + 9 + 3 * 10 + 4 * 11) +
        (2 * 15 + 16 + 17 + 18);
    try std.testing.expectEqual(
        @as(f64, expected_ions),
        inventory.ion_inventory_mol,
    );
}

test "storage publication cannot overwrite cumulative EXEC ledgers" {
    var totals = std.mem.zeroes(audit.Totals);
    totals.landscape_area_m2 = 10;
    totals.cumulative_rain_m3 = 17;
    totals.cumulative_nitrogen_output_g = 19;
    try publishStorage(&totals, .{
        .water_m3 = 1,
        .organic_carbon_g = 2,
        .ion_inventory_mol = 3,
    });
    try std.testing.expectEqual(@as(f64, 1), totals.water_storage_m3);
    try std.testing.expectEqual(@as(f64, 2), totals.organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 3), totals.ion_inventory_mol);
    try std.testing.expectEqual(@as(f64, 17), totals.cumulative_rain_m3);
    try std.testing.expectEqual(
        @as(f64, 19),
        totals.cumulative_nitrogen_output_g,
    );
}

test "REDIST surface physical inventory includes vapor and four gas phases" {
    var surface = try surface_precipitation.RuntimeState.init(
        std.testing.allocator,
        2,
    );
    defer surface.deinit();
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var gas_state = try gas.State.init(std.testing.allocator, 2);
    defer gas_state.deinit();
    var organic_state = try organic.State.init(std.testing.allocator, 2);
    defer organic_state.deinit();
    surface.litter_water_m3[0] = 1;
    surface.litter_water_m3[1] = 2;
    grid.surface_temperature_k[0] = 250;
    grid.surface_temperature_k[1] = 300;
    gas_state.water_vapor_mol[0] = 10;
    gas_state.water_vapor_mol[1] = 20;
    organic_state.structural[0].carbon_g_c = 100;
    organic_state.structural[
        organic.substrate_count * organic.structural_fraction_count
    ].carbon_g_c = 200;
    for (0..2) |cell| {
        const first = cell * gas.species_count;
        inline for (.{ gas_state.gaseous_mass_g, gas_state.dissolved_mass_g, gas_state.macropore_dissolved_mass_g, gas_state.band_dissolved_mass_g }) |phase| {
            phase[first + @intFromEnum(gas.Species.carbon_dioxide)] = 1;
            phase[first + @intFromEnum(gas.Species.methane)] = 2;
            phase[first + @intFromEnum(gas.Species.oxygen)] = 3;
            phase[first + @intFromEnum(gas.Species.nitrogen)] = 4;
            phase[first + @intFromEnum(gas.Species.nitrous_oxide)] = 5;
            phase[first + @intFromEnum(gas.Species.ammonia)] = 6;
        }
    }
    const ice = [_]f64{ 0.5, 0.25 };
    const parameters: SurfacePhysicalParameters = .{
        .dry_organic_heat_capacity_megajoules_per_g_c_k = 2.5e-6,
        .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
        .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
        .water_molar_mass_g_per_mol = 18,
        .liquid_water_density_g_per_m3 = 1e6,
    };
    const inventory = try aggregateSurfacePhysicalAndGas(
        &surface,
        &ice,
        &grid,
        &gas_state,
        &organic_state,
        parameters,
    );
    const vapor0 = 10.0 * 18.0 / 1e6;
    const vapor1 = 20.0 * 18.0 / 1e6;
    try std.testing.expectApproxEqAbs(
        1 + 2 + 0.5 + 0.25 + vapor0 + vapor1,
        inventory.water_m3,
        1e-14,
    );
    // HEAT-001 resolution A, second layer. Surface litter/pond ice is a frozen
    // carrier and carries `C_l*Tm - L + C_i*(T - Tm)` per cubic metre, not
    // `C_i*T - L`. The two cells sit at 250 K and 300 K, so this expectation
    // also pins that the ice branch is evaluated at each cell's own
    // temperature rather than at the melting point.
    const frozen_enthalpy_at = struct {
        fn f(temperature_k: f64) f64 {
            return 4.19 * 273.15 - 333 + 1.9274 * (temperature_k - 273.15);
        }
    }.f;
    const expected_heat =
        (2.5e-6 * 100 + 4.19 * (1 + vapor0)) * 250 +
        (2.5e-6 * 200 + 4.19 * (2 + vapor1)) * 300 +
        frozen_enthalpy_at(250) * 0.5 +
        frozen_enthalpy_at(300) * 0.25;
    try std.testing.expectApproxEqAbs(expected_heat, inventory.heat_megajoules, 1e-10);
    try std.testing.expectEqual(@as(f64, 24), inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(@as(f64, 24), inventory.oxygen_g);
    try std.testing.expectEqual(@as(f64, 72), inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 48), inventory.ammonium_nitrogen_g);
}

test "EXTRACT canopy inventory supports more than five runtime species" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 7,
        },
        .{ .worker_threads = 3, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var plants = try grid_module.PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var retention = try canopy_retention.State.init(
        std.testing.allocator,
        2,
        7,
    );
    defer retention.deinit();
    for (0..14) |plant| {
        plants.canopy_water_storage_m_per_m2[plant] = 0.001;
        retention.living_surface_water_m3[plant] = 0.01;
        retention.standing_dead_surface_water_m3[plant] = 0.02;
        retention.previous_water_energy_megajoules[plant] = 3;
    }
    const inventory = try aggregateCanopyWaterAndHeat(
        &plants,
        &retention,
        &.{ 10, 20 },
    );
    try std.testing.expectApproxEqAbs(
        7 * (0.001 * 10 + 0.01 + 0.02) +
            7 * (0.001 * 20 + 0.01 + 0.02),
        inventory.water_m3,
        1e-14,
    );
    try std.testing.expectEqual(@as(f64, 42), inventory.heat_megajoules);
}

test "REDIST profile mineral nitrogen counts each runtime owner once" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 2,
            .plant_populations = 1,
        },
        .{ .worker_threads = 2, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 1;
    grid.active_soil_layer_count[1] = 1;

    var transport = try mineral_nitrogen.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer transport.deinit();
    var chemistry = try soil_chemistry.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer chemistry.deinit();
    var fertilizer = try nitrogen_fertilizer.State.init(
        std.testing.allocator,
        grid.cell_count,
        grid.soil_layer_capacity,
    );
    defer fertilizer.deinit();

    // Active profile cells are 0 and 2. Matrix and macropore are both
    // authoritative and must be included.
    for ([_]usize{ 0, 2 }) |profile_cell| {
        const matrix = try transport.matrix.cellAmounts(profile_cell);
        const macropore = try transport.macropore.cellAmounts(profile_cell);
        matrix[@intFromEnum(mineral_nitrogen.Species.ammonium_non_band)] = 1;
        macropore[@intFromEnum(mineral_nitrogen.Species.ammonia_band)] = 2;
        matrix[@intFromEnum(mineral_nitrogen.Species.nitrate_non_band)] = 3;
        macropore[@intFromEnum(mineral_nitrogen.Species.nitrite_band)] = 4;
        chemistry.cation_exchange_mol_per_megagram[profile_cell]
            .ammonium_non_band = 0.5;
        fertilizer.soil[profile_cell].broadcast_ammonium_mol_n = 5;
        fertilizer.soil[profile_cell].banded_urea_mol_n = 6;
        fertilizer.soil[profile_cell].broadcast_nitrate_mol_n = 7;
    }
    // Inactive allocated capacity must never enter an authoritative total.
    (try transport.matrix.cellAmounts(1))[0] = 1_000_000;
    chemistry.cation_exchange_mol_per_megagram[3].ammonium_band = 1_000_000;
    fertilizer.soil[1].broadcast_ammonium_mol_n = 1_000_000;

    const inventory = try aggregateProfileMineralNitrogen(
        &grid,
        &transport,
        &chemistry,
        &fertilizer,
        &.{ 2, 2, 2, 2 },
        14,
    );
    // Per active cell: NHx=(1+2)+(0.5*2)+(5+6)=15 mol N;
    // NOx=(3+4)+7=14 mol N.
    try std.testing.expectEqual(@as(f64, 2 * 15 * 14), inventory.ammonium_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 2 * 14 * 14), inventory.nitrate_nitrogen_g);
    // Per cell: exchange NH4 2*1 + dry NH4 2*5 + urea 6 + nitrate 7.
    try std.testing.expectEqual(@as(f64, 2 * 25), inventory.ion_inventory_mol);
}

test "REDIST profile phosphorus and ions include matrix macropore and immobile owners" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 2,
            .plant_populations = 1,
        },
        .{ .worker_threads = 2, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 1;
    grid.active_soil_layer_count[1] = 1;

    var micropore = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        solute_species.AqueousSpecies.count,
    );
    defer micropore.deinit();
    var macropore = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        solute_species.AqueousSpecies.count,
    );
    defer macropore.deinit();
    var chemistry = try soil_chemistry.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer chemistry.deinit();
    var pending = try mineral_fertilizer.State.init(
        std.testing.allocator,
        grid.cell_count,
        grid.soil_layer_capacity,
    );
    defer pending.deinit();

    for ([_]usize{ 0, 2 }) |profile_cell| {
        const matrix = try micropore.cellAmounts(profile_cell);
        const macro = try macropore.cellAmounts(profile_cell);
        matrix[solute_species.index(.aluminum)] = 1;
        matrix[solute_species.index(.bicarbonate)] = 2;
        matrix[solute_species.index(.aluminum_hydroxide_2)] = 3;
        matrix[solute_species.index(.aluminum_hydroxide_3)] = 4;
        matrix[solute_species.index(.aluminum_hydroxide_4)] = 5;
        matrix[solute_species.index(.non_band_phosphate)] = 6;
        matrix[solute_species.index(.band_iron_hpo4)] = 7;
        macro[solute_species.index(.band_phosphoric_acid)] = 8;

        const non_band = &chemistry.non_band_phosphate[profile_cell];
        non_band.dissolved_hpo4_mol_p_per_m3 = 1;
        non_band.adsorbed_hpo4_mol_p_per_megagram = 2;
        non_band.deprotonated_site_mol_per_megagram = 1;
        non_band.monocalcium_phosphate_solid_mol_per_m3 = 1;
        const band = &chemistry.band_phosphate[profile_cell];
        band.dissolved_h2po4_mol_p_per_m3 = 2;
        band.adsorbed_h2po4_mol_p_per_megagram = 1;
        band.hydroxyl_site_mol_per_megagram = 1;
        band.hydroxyapatite_solid_mol_per_m3 = 2;

        chemistry.cation_exchange_mol_per_megagram[profile_cell].hydrogen = 1;
        chemistry.cation_exchange_mol_per_megagram[profile_cell].calcium = 2;
        chemistry.carboxyl_bound_hydrogen_mol_per_megagram[profile_cell] = 0.5;
        chemistry.geochemistry_solids[profile_cell]
            .calcite_solid_mol_per_m3 = 1;
        chemistry.geochemistry_solids[profile_cell]
            .gibbsite_solid_mol_per_m3 = 2;
        chemistry.geochemistry_solids[profile_cell]
            .aluminum_natural_silicate_mol_per_m3 = 3;

        pending.soil[profile_cell].broadcast_monocalcium_phosphate_mol = 1;
        pending.soil[profile_cell].hydroxyapatite_mol = 2;
        pending.soil[profile_cell].calcite_mol = 3;
        pending.soil[profile_cell].aluminum_ground_silicate_mol = 4;
    }
    // Capacity layers 1 and 3 are not part of the runtime profile.
    (try micropore.cellAmounts(1))[solute_species.index(.band_phosphate)] =
        1_000_000;
    chemistry.non_band_phosphate[3]
        .hydroxyapatite_solid_mol_per_m3 = 1_000_000;
    pending.soil[1].hydroxyapatite_mol = 1_000_000;

    const inventory = try aggregateProfilePhosphorusAndIons(
        &grid,
        &micropore,
        &macropore,
        &chemistry,
        &pending,
        &.{ 3, 3, 3, 3 },
        &.{ 2, 2, 2, 2 },
        .{
            .ammonium_non_band = 0.75,
            .ammonium_band = 0.25,
            .nitrate_non_band = 0.75,
            .nitrate_band = 0.25,
            .phosphate_non_band = 0.75,
            .phosphate_band = 0.25,
        },
        12,
        31,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * 45.25 * 31),
        inventory.phosphate_phosphorus_g,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * 250.5),
        inventory.ion_inventory_mol,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * 96),
        inventory.carbon_dioxide_carbon_g,
    );
}

test "REDIST dry surface mineral fertilizer remains in landscape storage" {
    var state = try mineral_fertilizer.State.init(
        std.testing.allocator,
        2,
        3,
    );
    defer state.deinit();
    state.surface[0].broadcast_monocalcium_phosphate_mol = 2;
    state.surface[0].banded_monocalcium_phosphate_mol = 3;
    state.surface[0].hydroxyapatite_mol = 4;
    state.surface[0].calcite_mol = 5;
    state.surface[0].gypsum_mol = 6;
    state.surface[0].aluminum_ground_silicate_mol = 7;
    state.surface[1].potassium_ground_silicate_mol = 8;

    const inventory = try aggregatePendingSurfaceMinerals(&state, 12, 31);
    try std.testing.expectEqual(
        @as(f64, (2 * (2 + 3) + 3 * 4) * 31),
        inventory.phosphate_phosphorus_g,
    );
    try std.testing.expectEqual(
        @as(f64, 5 * 12),
        inventory.carbon_dioxide_carbon_g,
    );
    try std.testing.expectEqual(
        @as(f64, 7 * (2 + 3) + 9 * 4 + 2 * (5 + 6) + 7 + 8),
        inventory.ion_inventory_mol,
    );
}
