const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const ion_species_count = @typeInfo(IonSpecies).@"enum".fields.len;

pub const IonSpecies = enum(u8) {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_monohydroxide,
    aluminum_dihydroxide,
    aluminum_trihydroxide,
    aluminum_tetrahydroxide,
    aluminum_sulfate,
    iron_monohydroxide,
    iron_dihydroxide,
    iron_trihydroxide,
    iron_tetrahydroxide,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    phosphate,
    phosphoric_acid,
    iron_monophosphate,
    iron_diphosphate,
    calcium_phosphate,
    calcium_hydrogen_phosphate,
    calcium_dihydrogen_phosphate,
    magnesium_hydrogen_phosphate,
};

pub const PrecipitationConcentrations = struct {
    carbon_dioxide_mol_per_m3: f64,
    methane_mol_per_m3: f64,
    oxygen_mol_per_m3: f64,
    nitrous_oxide_mol_n_per_m3: f64,
    dinitrogen_mol_n_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
    nitrate_mol_n_per_m3: f64,
    hydrogen_phosphate_mol_p_per_m3: f64,
    organic_phosphorus_mol_p_per_m3: f64,
    ions_mol_per_m3: [ion_species_count]f64,
};

pub const Parameters = struct {
    ice_density_Mg_per_m3: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    phosphorus_molar_mass_g_per_mol: f64,
};

pub const LayerPhysicalState = struct {
    heat_capacity_mj_per_k: f64,
    active_heat_capacity_threshold_mj_per_k: f64,
    liquid_water_m3: f64,
    solid_snow_water_equivalent_m3: f64,
    ice_m3: f64,
};

pub const PrimaryInventories = struct {
    carbon_dioxide_mol: f64,
    methane_mol: f64,
    oxygen_mol: f64,
    nitrous_oxide_mol_n: f64,
    dinitrogen_mol_n: f64,
    ammonium_g_n: f64,
    ammonia_g_n: f64,
    nitrate_g_n: f64,
    hydrogen_phosphate_g_p: f64,
    organic_phosphorus_g_p: f64,
};

pub const LayerState = struct {
    primary: PrimaryInventories,
    ions_mol: [ion_species_count]f64,
};

fn validatePrecipitation(concentrations: PrecipitationConcentrations) !void {
    inline for (std.meta.fields(PrecipitationConcentrations)) |field| {
        if (comptime std.mem.eql(u8, field.name, "ions_mol_per_m3")) continue;
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSnowpackPrecipitationConcentration;
    }
    for (concentrations.ions_mol_per_m3) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSnowpackPrecipitationConcentration;
}

fn equivalentWaterVolume(
    physical: LayerPhysicalState,
    ice_density_Mg_per_m3: f64,
) f64 {
    return physical.liquid_water_m3 +
        physical.solid_snow_water_equivalent_m3 +
        physical.ice_m3 * ice_density_Mg_per_m3;
}

/// Direct translation of STARTE lines 2072--2200 over a runtime snow-layer
/// slice. Active layers inherit precipitation chemistry in their total
/// water-equivalent volume; inactive layers receive zero primary inventory.
/// STARTE touches ion inventories only when dynamic salt chemistry is enabled,
/// so static mode deliberately retains their caller-initialized values.
pub fn initialize(
    layers: []LayerState,
    physical_layers: []const LayerPhysicalState,
    concentrations: PrecipitationConcentrations,
    mode: SaltEquilibriumMode,
    parameters: Parameters,
) !void {
    if (layers.len == 0 or layers.len != physical_layers.len)
        return error.SnowpackChemicalDimensionMismatch;
    inline for (.{
        parameters.ice_density_Mg_per_m3,
        parameters.nitrogen_molar_mass_g_per_mol,
        parameters.phosphorus_molar_mass_g_per_mol,
    }) |value| if (!std.math.isFinite(value) or value <= 0)
        return error.InvalidSnowpackChemicalParameter;
    try validatePrecipitation(concentrations);

    for (physical_layers) |physical| {
        inline for (std.meta.fields(LayerPhysicalState)) |field| {
            const value = @field(physical, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidSnowpackPhysicalState;
        }
        const volume_m3 = equivalentWaterVolume(
            physical,
            parameters.ice_density_Mg_per_m3,
        );
        if (!std.math.isFinite(volume_m3))
            return error.NonFiniteSnowpackEquivalentWaterVolume;
        if (physical.heat_capacity_mj_per_k >
            physical.active_heat_capacity_threshold_mj_per_k)
        {
            const checks = [_]f64{
                volume_m3 * concentrations.carbon_dioxide_mol_per_m3,
                volume_m3 * concentrations.methane_mol_per_m3,
                volume_m3 * concentrations.oxygen_mol_per_m3,
                volume_m3 * concentrations.nitrous_oxide_mol_n_per_m3,
                volume_m3 * concentrations.dinitrogen_mol_n_per_m3,
                volume_m3 * concentrations.ammonium_mol_n_per_m3 *
                    parameters.nitrogen_molar_mass_g_per_mol,
                volume_m3 * concentrations.ammonia_mol_n_per_m3 *
                    parameters.nitrogen_molar_mass_g_per_mol,
                volume_m3 * concentrations.nitrate_mol_n_per_m3 *
                    parameters.nitrogen_molar_mass_g_per_mol,
                volume_m3 * concentrations.hydrogen_phosphate_mol_p_per_m3 *
                    parameters.phosphorus_molar_mass_g_per_mol,
                volume_m3 * concentrations.organic_phosphorus_mol_p_per_m3 *
                    parameters.phosphorus_molar_mass_g_per_mol,
            };
            for (checks) |value| if (!std.math.isFinite(value))
                return error.NonFiniteSnowpackChemicalInventory;
            if (mode == .dynamic) {
                for (concentrations.ions_mol_per_m3) |value|
                    if (!std.math.isFinite(volume_m3 * value))
                        return error.NonFiniteSnowpackChemicalInventory;
            }
        }
    }

    for (layers, physical_layers) |*layer, physical| {
        const active = physical.heat_capacity_mj_per_k >
            physical.active_heat_capacity_threshold_mj_per_k;
        if (!active) {
            layer.primary = std.mem.zeroes(PrimaryInventories);
            if (mode == .dynamic) layer.ions_mol = [_]f64{0} ** ion_species_count;
            continue;
        }
        const volume_m3 = equivalentWaterVolume(
            physical,
            parameters.ice_density_Mg_per_m3,
        );
        layer.primary = .{
            .carbon_dioxide_mol = volume_m3 * concentrations.carbon_dioxide_mol_per_m3,
            .methane_mol = volume_m3 * concentrations.methane_mol_per_m3,
            .oxygen_mol = volume_m3 * concentrations.oxygen_mol_per_m3,
            .nitrous_oxide_mol_n = volume_m3 * concentrations.nitrous_oxide_mol_n_per_m3,
            .dinitrogen_mol_n = volume_m3 * concentrations.dinitrogen_mol_n_per_m3,
            .ammonium_g_n = volume_m3 * concentrations.ammonium_mol_n_per_m3 *
                parameters.nitrogen_molar_mass_g_per_mol,
            .ammonia_g_n = volume_m3 * concentrations.ammonia_mol_n_per_m3 *
                parameters.nitrogen_molar_mass_g_per_mol,
            .nitrate_g_n = volume_m3 * concentrations.nitrate_mol_n_per_m3 *
                parameters.nitrogen_molar_mass_g_per_mol,
            .hydrogen_phosphate_g_p = volume_m3 * concentrations.hydrogen_phosphate_mol_p_per_m3 *
                parameters.phosphorus_molar_mass_g_per_mol,
            .organic_phosphorus_g_p = volume_m3 * concentrations.organic_phosphorus_mol_p_per_m3 *
                parameters.phosphorus_molar_mass_g_per_mol,
        };
        if (mode == .dynamic) {
            for (
                concentrations.ions_mol_per_m3,
                &layer.ions_mol,
            ) |concentration, *inventory| inventory.* = volume_m3 * concentration;
        }
    }
}

fn testConcentrations() PrecipitationConcentrations {
    return .{
        .carbon_dioxide_mol_per_m3 = 1,
        .methane_mol_per_m3 = 2,
        .oxygen_mol_per_m3 = 3,
        .nitrous_oxide_mol_n_per_m3 = 4,
        .dinitrogen_mol_n_per_m3 = 5,
        .ammonium_mol_n_per_m3 = 6,
        .ammonia_mol_n_per_m3 = 7,
        .nitrate_mol_n_per_m3 = 8,
        .hydrogen_phosphate_mol_p_per_m3 = 9,
        .organic_phosphorus_mol_p_per_m3 = 10,
        .ions_mol_per_m3 = [_]f64{2} ** ion_species_count,
    };
}

const test_parameters: Parameters = .{
    .ice_density_Mg_per_m3 = 0.9,
    .nitrogen_molar_mass_g_per_mol = 14,
    .phosphorus_molar_mass_g_per_mol = 31,
};

test "STARTE active snow layer uses total water-equivalent volume" {
    var layers = [_]LayerState{undefined};
    try initialize(&layers, &.{.{
        .heat_capacity_mj_per_k = 2,
        .active_heat_capacity_threshold_mj_per_k = 1,
        .liquid_water_m3 = 1,
        .solid_snow_water_equivalent_m3 = 2,
        .ice_m3 = 10,
    }}, testConcentrations(), .dynamic, test_parameters);
    const volume_m3: f64 = 12;
    try std.testing.expectEqual(volume_m3, layers[0].primary.carbon_dioxide_mol);
    try std.testing.expectEqual(volume_m3 * 6 * 14, layers[0].primary.ammonium_g_n);
    try std.testing.expectEqual(volume_m3 * 9 * 31, layers[0].primary.hydrogen_phosphate_g_p);
    try std.testing.expectEqual(volume_m3 * 2, layers[0].ions_mol[@intFromEnum(IonSpecies.potassium_sulfate)]);
}

test "STARTE inactive layer zeros tracked chemistry in dynamic mode" {
    var layers = [_]LayerState{undefined};
    try initialize(&layers, &.{.{
        .heat_capacity_mj_per_k = 1,
        .active_heat_capacity_threshold_mj_per_k = 1,
        .liquid_water_m3 = 4,
        .solid_snow_water_equivalent_m3 = 5,
        .ice_m3 = 6,
    }}, testConcentrations(), .dynamic, test_parameters);
    try std.testing.expectEqual(@as(f64, 0), layers[0].primary.nitrate_g_n);
    try std.testing.expectEqual(@as(f64, 0), layers[0].ions_mol[0]);
}

test "STARTE static salt mode retains caller ion inventories" {
    var layers = [_]LayerState{.{
        .primary = undefined,
        .ions_mol = [_]f64{7} ** ion_species_count,
    }};
    try initialize(&layers, &.{.{
        .heat_capacity_mj_per_k = 0,
        .active_heat_capacity_threshold_mj_per_k = 1,
        .liquid_water_m3 = 0,
        .solid_snow_water_equivalent_m3 = 0,
        .ice_m3 = 0,
    }}, testConcentrations(), .static, test_parameters);
    try std.testing.expectEqual(@as(f64, 0), layers[0].primary.oxygen_mol);
    try std.testing.expectEqual(@as(f64, 7), layers[0].ions_mol[0]);
}
