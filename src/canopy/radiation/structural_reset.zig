const std = @import("std");

pub const Inputs = struct {
    initial_clumping_factor: []const f64,
    projected_leaf_area_m2: []const f64,
    cell_area_m2: f64,
};

pub const Outputs = struct {
    absorbed_shortwave_megajoules_m2_h: []f64,
    downward_shortwave_megajoules_m2_h: []f64,
    absorbed_par_umol_m2_s: []f64,
    downward_par_umol_m2_s: []f64,
    effective_clumping_factor: []f64,
};

pub const ResetTotals = struct {
    canopy_shortwave_megajoules_m2_h: f64 = 0,
    canopy_par_umol_m2_s: f64 = 0,
};

/// `hour1.f` lines 1023--1032. Resets cell/species canopy radiation
/// accumulators and computes CFX for a runtime number of plant species.
pub fn apply(inputs: Inputs, outputs: Outputs) !ResetTotals {
    try validate(inputs, outputs);
    @memset(outputs.absorbed_shortwave_megajoules_m2_h, 0);
    @memset(outputs.downward_shortwave_megajoules_m2_h, 0);
    @memset(outputs.absorbed_par_umol_m2_s, 0);
    @memset(outputs.downward_par_umol_m2_s, 0);
    for (inputs.initial_clumping_factor, inputs.projected_leaf_area_m2, 0..) |
        initial_clumping,
        projected_leaf_area_m2,
        species,
    | {
        const density_correction =
            1.0 - 0.025 * projected_leaf_area_m2 / inputs.cell_area_m2;
        if (density_correction < 0)
            return error.NegativeCanopyClumpingCorrection;
        outputs.effective_clumping_factor[species] =
            initial_clumping * density_correction;
    }
    return .{};
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    const species_count = inputs.initial_clumping_factor.len;
    if (species_count == 0 or
        inputs.projected_leaf_area_m2.len != species_count or
        outputs.absorbed_shortwave_megajoules_m2_h.len != species_count or
        outputs.downward_shortwave_megajoules_m2_h.len != species_count or
        outputs.absorbed_par_umol_m2_s.len != species_count or
        outputs.downward_par_umol_m2_s.len != species_count or
        outputs.effective_clumping_factor.len != species_count)
        return error.CanopyStructuralResetDimensionMismatch;
    if (!std.math.isFinite(inputs.cell_area_m2) or inputs.cell_area_m2 <= 0)
        return error.InvalidCanopyCellArea;
    for (inputs.initial_clumping_factor, inputs.projected_leaf_area_m2) |
        initial_clumping,
        projected_leaf_area_m2,
    | {
        if (!std.math.isFinite(initial_clumping) or
            !std.math.isFinite(projected_leaf_area_m2))
            return error.NonFiniteCanopyStructuralInput;
        if (initial_clumping < 0 or projected_leaf_area_m2 < 0)
            return error.InvalidCanopyStructuralInput;
        if (1.0 - 0.025 * projected_leaf_area_m2 / inputs.cell_area_m2 < 0)
            return error.NegativeCanopyClumpingCorrection;
    }
}

test "runtime species radiation accumulators reset and CFX follows source" {
    var absorbed_sw = [_]f64{ 1, 2, 3 };
    var downward_sw = [_]f64{ 4, 5, 6 };
    var absorbed_par = [_]f64{ 7, 8, 9 };
    var downward_par = [_]f64{ 10, 11, 12 };
    var clumping: [3]f64 = undefined;
    const totals = try apply(.{
        .initial_clumping_factor = &.{ 0.8, 0.9, 1.0 },
        .projected_leaf_area_m2 = &.{ 10, 20, 30 },
        .cell_area_m2 = 100,
    }, .{
        .absorbed_shortwave_megajoules_m2_h = &absorbed_sw,
        .downward_shortwave_megajoules_m2_h = &downward_sw,
        .absorbed_par_umol_m2_s = &absorbed_par,
        .downward_par_umol_m2_s = &downward_par,
        .effective_clumping_factor = &clumping,
    });
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &absorbed_sw);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &downward_par);
    try std.testing.expectApproxEqAbs(@as(f64, 0.798), clumping[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8955), clumping[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), totals.canopy_shortwave_megajoules_m2_h);
    try std.testing.expectEqual(@as(f64, 0), totals.canopy_par_umol_m2_s);
}

test "invalid clumping correction leaves all output unchanged" {
    var absorbed_sw = [_]f64{42};
    var downward_sw = [_]f64{43};
    var absorbed_par = [_]f64{44};
    var downward_par = [_]f64{45};
    var clumping = [_]f64{46};
    try std.testing.expectError(error.NegativeCanopyClumpingCorrection, apply(.{
        .initial_clumping_factor = &.{1},
        .projected_leaf_area_m2 = &.{41},
        .cell_area_m2 = 1,
    }, .{
        .absorbed_shortwave_megajoules_m2_h = &absorbed_sw,
        .downward_shortwave_megajoules_m2_h = &downward_sw,
        .absorbed_par_umol_m2_s = &absorbed_par,
        .downward_par_umol_m2_s = &downward_par,
        .effective_clumping_factor = &clumping,
    }));
    try std.testing.expectEqual(@as(f64, 42), absorbed_sw[0]);
    try std.testing.expectEqual(@as(f64, 46), clumping[0]);
}
