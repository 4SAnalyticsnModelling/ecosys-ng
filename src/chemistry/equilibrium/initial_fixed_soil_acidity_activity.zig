const std = @import("std");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    layer_count_including_surface: usize,
};

pub const Result = struct {
    hydrogen_activity_mol_per_m3: f64, // AHY1
    hydroxide_activity_mol_per_m3: f64, // AOH1
};

/// Direct translation of `starte.f` lines 439--440. Soil pH is cell-major,
/// then layer, including surface layer zero. The enclosing source loop does
/// not guard these assignments: each source uses the active soil layer's
/// fixed acidity while the initialization equilibrium is solved.
pub fn calculate(
    dimensions: Dimensions,
    row: usize,
    column: usize,
    layer: usize,
    soil_ph_by_cell_layer: []const f64,
    water_dissociation_mol2_per_m6: f64,
) !Result {
    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.layer_count_including_surface == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count or
        layer >= dimensions.layer_count_including_surface)
        return error.InvalidInitialAcidityDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidInitialAcidityDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.layer_count_including_surface) catch
        return error.InvalidInitialAcidityDimensions;
    if (soil_ph_by_cell_layer.len != value_count)
        return error.InvalidInitialAcidityDimensions;
    if (!std.math.isFinite(water_dissociation_mol2_per_m6) or
        water_dissociation_mol2_per_m6 <= 0)
        return error.InvalidInitialAcidityConstant;

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.layer_count_including_surface + layer;
    const soil_ph = soil_ph_by_cell_layer[value_index];
    if (!std.math.isFinite(soil_ph) or soil_ph < 0 or soil_ph > 14)
        return error.InvalidInitialSoilPh;
    const hydrogen_activity = std.math.pow(f64, 10.0, -(soil_ph - 3.0));
    const result: Result = .{
        .hydrogen_activity_mol_per_m3 = hydrogen_activity,
        .hydroxide_activity_mol_per_m3 = water_dissociation_mol2_per_m6 /
            hydrogen_activity,
    };
    if (!std.math.isFinite(result.hydrogen_activity_mol_per_m3) or
        !std.math.isFinite(result.hydroxide_activity_mol_per_m3))
        return error.NonFiniteInitialAcidityActivity;
    return result;
}

test "STARTE fixed acidity uses runtime cell-layer soil pH topology" {
    const ph = [_]f64{ 4, 5, 6, 7 };
    const result = try calculate(
        .{ .column_count = 2, .row_count = 1, .layer_count_including_surface = 2 },
        0,
        1,
        1,
        &ph,
        1.0e-8,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0e-4),
        result.hydrogen_activity_mol_per_m3,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0e-4),
        result.hydroxide_activity_mol_per_m3,
        1.0e-16,
    );
}

test "STARTE fixed acidity uses runtime water dissociation constant" {
    const ph = [_]f64{6};
    const result = try calculate(
        .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 },
        0,
        0,
        0,
        &ph,
        2.0e-7,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 2.0e-4),
        result.hydroxide_activity_mol_per_m3,
        1.0e-16,
    );
}

test "STARTE fixed acidity rejects dimensions and late invalid pH atomically" {
    const ph = [_]f64{ 7, std.math.nan(f64) };
    try std.testing.expectError(error.InvalidInitialSoilPh, calculate(
        .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 2 },
        0,
        0,
        1,
        &ph,
        1.0e-8,
    ));
    try std.testing.expectError(error.InvalidInitialAcidityDimensions, calculate(
        .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 },
        0,
        0,
        0,
        &ph,
        1.0e-8,
    ));
}
