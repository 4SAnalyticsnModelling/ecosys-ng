const std = @import("std");
const properties_module = @import("../water/solver_properties.zig");

const extensive_fields = .{
    "sand_mass_megagrams",
    "silt_mass_megagrams",
    "clay_mass_megagrams",
    "cation_exchange_capacity_mol",
    "anion_exchange_capacity_mol",
};

/// REDIST soil branch (legacy REDIST.F 9616 onward): transfer the requested
/// fraction of sediment and exchange capacity, then rebuild their intensive
/// views on the accepted post-transfer dry-soil masses.
pub fn transferLayerFraction(
    properties: *properties_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    source_soil_mass_after_megagrams: f64,
    destination_soil_mass_after_megagrams: f64,
) !void {
    try validateLayerFraction(properties, source, destination, fraction, source_soil_mass_after_megagrams, destination_soil_mass_after_megagrams);
    inline for (extensive_fields) |field_name| {
        const values = @field(properties, field_name);
        const moved = fraction * values[source];
        values[destination] += moved;
        values[source] -= moved;
    }
    refresh(properties, source, source_soil_mass_after_megagrams);
    refresh(properties, destination, destination_soil_mass_after_megagrams);
}

pub fn validateLayerFraction(
    properties: *const properties_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    source_soil_mass_after_megagrams: f64,
    destination_soil_mass_after_megagrams: f64,
) !void {
    if (source >= properties.layer_count or destination >= properties.layer_count or source == destination) return error.MineralLayerRemapIndexOutOfBounds;
    inline for (.{ fraction, source_soil_mass_after_megagrams, destination_soil_mass_after_megagrams }) |value| if (!std.math.isFinite(value)) return error.InvalidMineralLayerRemapInput;
    if (fraction < 0 or fraction > 1 or source_soil_mass_after_megagrams < 0 or destination_soil_mass_after_megagrams <= 0) return error.InvalidMineralLayerRemapInput;
    inline for (extensive_fields) |field_name| {
        const values = @field(properties, field_name);
        if (values.len != properties.layer_count) return error.MineralLayerRemapDimensionMismatch;
        try validatePair(values[source], values[destination]);
    }
}

fn validatePair(source: f64, destination: f64) !void {
    if (!std.math.isFinite(source) or source < 0 or !std.math.isFinite(destination) or destination < 0 or !std.math.isFinite(source + destination)) return error.InvalidMineralLayerRemapState;
}

fn refresh(properties: *properties_module.State, layer: usize, soil_mass_megagrams: f64) void {
    if (soil_mass_megagrams > 0) {
        properties.sand_mass_fraction[layer] = properties.sand_mass_megagrams[layer] / soil_mass_megagrams;
        properties.clay_mass_fraction[layer] = properties.clay_mass_megagrams[layer] / soil_mass_megagrams;
        properties.cation_exchange_capacity_mol_per_megagram[layer] = properties.cation_exchange_capacity_mol[layer] / soil_mass_megagrams;
        properties.anion_exchange_capacity_mol_per_megagram[layer] = properties.anion_exchange_capacity_mol[layer] / soil_mass_megagrams;
    } else {
        properties.sand_mass_fraction[layer] = 0;
        properties.clay_mass_fraction[layer] = 0;
        properties.cation_exchange_capacity_mol_per_megagram[layer] = 0;
        properties.anion_exchange_capacity_mol_per_megagram[layer] = 0;
    }
}

test "REDIST mineral sediment transfers every requested soil fraction" {
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
    properties.sand_mass_megagrams = &sand;
    properties.silt_mass_megagrams = &silt;
    properties.clay_mass_megagrams = &clay;
    properties.cation_exchange_capacity_mol = &cec_mol;
    properties.anion_exchange_capacity_mol = &aec_mol;
    properties.sand_mass_fraction = &sand_fraction;
    properties.clay_mass_fraction = &clay_fraction;
    properties.cation_exchange_capacity_mol_per_megagram = &cec;
    properties.anion_exchange_capacity_mol_per_megagram = &aec;
    try transferLayerFraction(&properties, 0, 1, 0.5, 5, 15);
    try std.testing.expectEqual(@as(f64, 3), sand[0]);
    try std.testing.expectEqual(@as(f64, 5), sand[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0 / 3.0), sand_fraction[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 50.0 / 3.0), cec[1], 1e-14);
}
