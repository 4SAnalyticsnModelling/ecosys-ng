const std = @import("std");
const properties_module = @import("soil_solver_properties.zig");

const extensive_fields = .{
    "sand_mass_Mg",
    "silt_mass_Mg",
    "clay_mass_Mg",
    "cation_exchange_capacity_mol",
    "anion_exchange_capacity_mol",
};

/// Exact REDIST sediment gate: SAND/SILT/CLAY/XCEC/XAEC move only for FX=1.
/// Resulting intensive property views use the caller's accepted post-transfer
/// total soil masses.
pub fn transferPondedLayerIfComplete(
    properties: *properties_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    source_soil_mass_after_Mg: f64,
    destination_soil_mass_after_Mg: f64,
) !void {
    try validatePondedLayerIfComplete(properties, source, destination, fraction, source_soil_mass_after_Mg, destination_soil_mass_after_Mg);
    if (fraction != 1) return;
    inline for (extensive_fields) |field_name| {
        const values = @field(properties, field_name);
        values[destination] += values[source];
        values[source] = 0;
    }
    refresh(properties, source, source_soil_mass_after_Mg);
    refresh(properties, destination, destination_soil_mass_after_Mg);
}

pub fn validatePondedLayerIfComplete(
    properties: *const properties_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    source_soil_mass_after_Mg: f64,
    destination_soil_mass_after_Mg: f64,
) !void {
    if (source >= properties.layer_count or destination >= properties.layer_count or source == destination) return error.MineralLayerRemapIndexOutOfBounds;
    inline for (.{ fraction, source_soil_mass_after_Mg, destination_soil_mass_after_Mg }) |value| if (!std.math.isFinite(value)) return error.InvalidMineralLayerRemapInput;
    if (fraction < 0 or fraction > 1 or source_soil_mass_after_Mg < 0 or destination_soil_mass_after_Mg <= 0) return error.InvalidMineralLayerRemapInput;
    if (fraction != 1) return;
    inline for (extensive_fields) |field_name| {
        const values = @field(properties, field_name);
        if (values.len != properties.layer_count) return error.MineralLayerRemapDimensionMismatch;
        try validatePair(values[source], values[destination]);
    }
}

fn validatePair(source: f64, destination: f64) !void {
    if (!std.math.isFinite(source) or source < 0 or !std.math.isFinite(destination) or destination < 0 or !std.math.isFinite(source + destination)) return error.InvalidMineralLayerRemapState;
}

fn refresh(properties: *properties_module.State, layer: usize, soil_mass_Mg: f64) void {
    if (soil_mass_Mg > 0) {
        properties.sand_mass_fraction[layer] = properties.sand_mass_Mg[layer] / soil_mass_Mg;
        properties.clay_mass_fraction[layer] = properties.clay_mass_Mg[layer] / soil_mass_Mg;
        properties.cation_exchange_capacity_mol_per_Mg[layer] = properties.cation_exchange_capacity_mol[layer] / soil_mass_Mg;
        properties.anion_exchange_capacity_mol_per_Mg[layer] = properties.anion_exchange_capacity_mol[layer] / soil_mass_Mg;
    } else {
        properties.sand_mass_fraction[layer] = 0;
        properties.clay_mass_fraction[layer] = 0;
        properties.cation_exchange_capacity_mol_per_Mg[layer] = 0;
        properties.anion_exchange_capacity_mol_per_Mg[layer] = 0;
    }
}

test "REDIST mineral sediment transfers only for complete ponding" {
    var properties: properties_module.State = undefined;
    properties.layer_count = 2;
    var sand = [_]f64{ 6, 2 };
    var silt = [_]f64{ 3, 5 };
    var clay = [_]f64{ 1, 3 };
    var cec_mol = [_]f64{ 100, 200 };
    var aec_mol = [_]f64{ 10, 20 };
    var sand_fraction = [_]f64{ 0.6, 0.2 };
    var clay_fraction = [_]f64{ 0.1, 0.3 };
    var cec = [_]f64{ 10, 20 };
    var aec = [_]f64{ 1, 2 };
    properties.sand_mass_Mg = &sand;
    properties.silt_mass_Mg = &silt;
    properties.clay_mass_Mg = &clay;
    properties.cation_exchange_capacity_mol = &cec_mol;
    properties.anion_exchange_capacity_mol = &aec_mol;
    properties.sand_mass_fraction = &sand_fraction;
    properties.clay_mass_fraction = &clay_fraction;
    properties.cation_exchange_capacity_mol_per_Mg = &cec;
    properties.anion_exchange_capacity_mol_per_Mg = &aec;
    try transferPondedLayerIfComplete(&properties, 0, 1, 0.5, 5, 15);
    try std.testing.expectEqual(@as(f64, 6), sand[0]);
    try transferPondedLayerIfComplete(&properties, 0, 1, 1, 0, 20);
    try std.testing.expectEqual(@as(f64, 0), sand[0]);
    try std.testing.expectEqual(@as(f64, 8), sand[1]);
    try std.testing.expectEqual(@as(f64, 0.4), sand_fraction[1]);
    try std.testing.expectEqual(@as(f64, 15), cec[1]);
}
