const std = @import("std");

pub const Decision = struct {
    apply_soil_transfer: bool,
    apply_band_geometry_transfer: bool,
};

pub const Inputs = struct {
    source_layer_fortran: usize,
    source_bulk_density_megagrams_per_m3: f64,
    destination_bulk_density_megagrams_per_m3: f64,
    zero_tolerance_megagrams_per_m3: f64,
    source_layer_thickness_m: f64,
    destination_layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
};

/// Exact REDIST.F soil-transfer branch at line 9439 and nested band-geometry
/// gates at lines 9490--9492. The latter execute only for a non-surface source
/// and when both source and destination depths strictly exceed DLYRM.
pub fn decide(inputs: Inputs) !Decision {
    const density_values = [_]f64{
        inputs.source_bulk_density_megagrams_per_m3,
        inputs.destination_bulk_density_megagrams_per_m3,
        inputs.zero_tolerance_megagrams_per_m3,
    };
    for (density_values) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilLayerTransferGateInput;
    if (inputs.zero_tolerance_megagrams_per_m3 < 0)
        return error.InvalidSoilLayerTransferGateThreshold;

    const apply_soil = inputs.source_bulk_density_megagrams_per_m3 > inputs.zero_tolerance_megagrams_per_m3 and
        inputs.destination_bulk_density_megagrams_per_m3 > inputs.zero_tolerance_megagrams_per_m3;
    if (!apply_soil or inputs.source_layer_fortran == 0) return .{
        .apply_soil_transfer = apply_soil,
        .apply_band_geometry_transfer = false,
    };
    const thickness_values = [_]f64{
        inputs.source_layer_thickness_m,
        inputs.destination_layer_thickness_m,
        inputs.minimum_layer_thickness_m,
    };
    for (thickness_values) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilLayerTransferGateInput;
    if (inputs.minimum_layer_thickness_m < 0)
        return error.InvalidSoilLayerTransferGateThreshold;
    return .{
        .apply_soil_transfer = apply_soil,
        .apply_band_geometry_transfer = inputs.destination_layer_thickness_m > inputs.minimum_layer_thickness_m and
            inputs.source_layer_thickness_m > inputs.minimum_layer_thickness_m,
    };
}

test "REDIST soil branch requires both positive bulk densities" {
    const active = try decide(.{
        .source_layer_fortran = 2,
        .source_bulk_density_megagrams_per_m3 = 1.2,
        .destination_bulk_density_megagrams_per_m3 = 1.3,
        .zero_tolerance_megagrams_per_m3 = 0,
        .source_layer_thickness_m = 0.2,
        .destination_layer_thickness_m = 0.3,
        .minimum_layer_thickness_m = 0.01,
    });
    try std.testing.expect(active.apply_soil_transfer);
    try std.testing.expect(active.apply_band_geometry_transfer);
}

test "surface source keeps soil branch but skips band geometry" {
    const decision = try decide(.{
        .source_layer_fortran = 0,
        .source_bulk_density_megagrams_per_m3 = 1.2,
        .destination_bulk_density_megagrams_per_m3 = 1.3,
        .zero_tolerance_megagrams_per_m3 = 0,
        .source_layer_thickness_m = std.math.nan(f64),
        .destination_layer_thickness_m = std.math.nan(f64),
        .minimum_layer_thickness_m = std.math.nan(f64),
    });
    try std.testing.expect(decision.apply_soil_transfer);
    try std.testing.expect(!decision.apply_band_geometry_transfer);
}

test "band geometry uses strict thickness comparisons" {
    var inputs = Inputs{
        .source_layer_fortran = 1,
        .source_bulk_density_megagrams_per_m3 = 1,
        .destination_bulk_density_megagrams_per_m3 = 1,
        .zero_tolerance_megagrams_per_m3 = 0,
        .source_layer_thickness_m = 0.01,
        .destination_layer_thickness_m = 0.02,
        .minimum_layer_thickness_m = 0.01,
    };
    try std.testing.expect(!(try decide(inputs)).apply_band_geometry_transfer);
    inputs.source_layer_thickness_m = 0.02;
    inputs.destination_layer_thickness_m = 0.01;
    try std.testing.expect(!(try decide(inputs)).apply_band_geometry_transfer);
}

test "bulk density equal to ZERO does not enter soil branch" {
    const decision = try decide(.{
        .source_layer_fortran = 1,
        .source_bulk_density_megagrams_per_m3 = 0.001,
        .destination_bulk_density_megagrams_per_m3 = 1,
        .zero_tolerance_megagrams_per_m3 = 0.001,
        .source_layer_thickness_m = std.math.nan(f64),
        .destination_layer_thickness_m = std.math.nan(f64),
        .minimum_layer_thickness_m = std.math.nan(f64),
    });
    try std.testing.expect(!decision.apply_soil_transfer);
    try std.testing.expect(!decision.apply_band_geometry_transfer);
}

test "invalid soil transfer gate values fail immediately" {
    var inputs = Inputs{
        .source_layer_fortran = 1,
        .source_bulk_density_megagrams_per_m3 = std.math.nan(f64),
        .destination_bulk_density_megagrams_per_m3 = 1,
        .zero_tolerance_megagrams_per_m3 = 0,
        .source_layer_thickness_m = 1,
        .destination_layer_thickness_m = 1,
        .minimum_layer_thickness_m = 0,
    };
    try std.testing.expectError(error.NonFiniteSoilLayerTransferGateInput, decide(inputs));
    inputs.source_bulk_density_megagrams_per_m3 = 1;
    inputs.minimum_layer_thickness_m = -1;
    try std.testing.expectError(error.InvalidSoilLayerTransferGateThreshold, decide(inputs));
}
