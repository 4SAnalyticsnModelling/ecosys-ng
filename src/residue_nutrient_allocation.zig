const std = @import("std");

pub const Parameters = struct {
    /// Fraction-fastest `[residue_class][biochemical_fraction]`,
    /// in g N (g C)-1.
    nitrogen_to_carbon_ratio: []const f64,
    /// Fraction-fastest `[residue_class][biochemical_fraction]`,
    /// in g P (g C)-1.
    phosphorus_to_carbon_ratio: []const f64,
    residue_class_count: usize,
    biochemical_fraction_count: usize,
};

pub const State = struct {
    /// CNOFC owner, fraction-fastest, in g N (g C)-1.
    nitrogen_to_carbon_ratio: []f64,
    /// CPOFC owner, fraction-fastest, in g P (g C)-1.
    phosphorus_to_carbon_ratio: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 190--213.
///
/// The three source residue classes and four biochemical fractions are
/// generalized to runtime extents. Within each residue class all N:C values
/// are assigned before its P:C values, preserving the Fortran statement order.
pub fn initialize(state: State, parameters: Parameters) !void {
    if (parameters.residue_class_count == 0 or
        parameters.biochemical_fraction_count == 0)
        return error.InvalidResidueNutrientDimensions;
    const value_count = std.math.mul(
        usize,
        parameters.residue_class_count,
        parameters.biochemical_fraction_count,
    ) catch return error.DimensionOverflow;
    if (parameters.nitrogen_to_carbon_ratio.len != value_count or
        parameters.phosphorus_to_carbon_ratio.len != value_count or
        state.nitrogen_to_carbon_ratio.len != value_count or
        state.phosphorus_to_carbon_ratio.len != value_count)
    {
        return error.ResidueNutrientDimensionMismatch;
    }

    for (parameters.nitrogen_to_carbon_ratio) |ratio| {
        if (!std.math.isFinite(ratio))
            return error.NonFiniteResidueNutrientRatio;
        if (ratio < 0) return error.InvalidResidueNutrientRatio;
    }
    for (parameters.phosphorus_to_carbon_ratio) |ratio| {
        if (!std.math.isFinite(ratio))
            return error.NonFiniteResidueNutrientRatio;
        if (ratio < 0) return error.InvalidResidueNutrientRatio;
    }

    for (0..parameters.residue_class_count) |residue_index| {
        const start =
            residue_index * parameters.biochemical_fraction_count;
        const end = start + parameters.biochemical_fraction_count;
        @memcpy(
            state.nitrogen_to_carbon_ratio[start..end],
            parameters.nitrogen_to_carbon_ratio[start..end],
        );
        @memcpy(
            state.phosphorus_to_carbon_ratio[start..end],
            parameters.phosphorus_to_carbon_ratio[start..end],
        );
    }
}

test "STARTS residue nutrient allocation reproduces lines 190 through 213" {
    const source_nitrogen = [_]f64{
        0.005, 0.005, 0.005, 0.020,
        0.020, 0.020, 0.020, 0.020,
        0.020, 0.020, 0.020, 0.020,
    };
    const source_phosphorus = [_]f64{
        0.0005, 0.0005, 0.0005, 0.0020,
        0.0020, 0.0020, 0.0020, 0.0020,
        0.0020, 0.0020, 0.0020, 0.0020,
    };
    var nitrogen = [_]f64{-1.0} ** 12;
    var phosphorus = [_]f64{-1.0} ** 12;

    try initialize(.{
        .nitrogen_to_carbon_ratio = &nitrogen,
        .phosphorus_to_carbon_ratio = &phosphorus,
    }, .{
        .nitrogen_to_carbon_ratio = &source_nitrogen,
        .phosphorus_to_carbon_ratio = &source_phosphorus,
        .residue_class_count = 3,
        .biochemical_fraction_count = 4,
    });

    try std.testing.expectEqualSlices(f64, &source_nitrogen, &nitrogen);
    try std.testing.expectEqualSlices(f64, &source_phosphorus, &phosphorus);
}

test "runtime residue and biochemical extents are not fixed to source limits" {
    const nitrogen_input = [_]f64{
        0.01, 0.02,
        0.03, 0.04,
        0.05, 0.06,
        0.07, 0.08,
        0.09, 0.10,
    };
    const phosphorus_input = [_]f64{
        0.001, 0.002,
        0.003, 0.004,
        0.005, 0.006,
        0.007, 0.008,
        0.009, 0.010,
    };
    var nitrogen = [_]f64{0.0} ** 10;
    var phosphorus = [_]f64{0.0} ** 10;

    try initialize(.{
        .nitrogen_to_carbon_ratio = &nitrogen,
        .phosphorus_to_carbon_ratio = &phosphorus,
    }, .{
        .nitrogen_to_carbon_ratio = &nitrogen_input,
        .phosphorus_to_carbon_ratio = &phosphorus_input,
        .residue_class_count = 5,
        .biochemical_fraction_count = 2,
    });

    try std.testing.expectEqualSlices(f64, &nitrogen_input, &nitrogen);
    try std.testing.expectEqualSlices(f64, &phosphorus_input, &phosphorus);
}

test "invalid late ratio leaves both destination arrays unchanged" {
    var nitrogen = [_]f64{ 7.0, 8.0, 9.0, 10.0 };
    var phosphorus = [_]f64{ 11.0, 12.0, 13.0, 14.0 };
    const nitrogen_before = nitrogen;
    const phosphorus_before = phosphorus;

    try std.testing.expectError(
        error.NonFiniteResidueNutrientRatio,
        initialize(.{
            .nitrogen_to_carbon_ratio = &nitrogen,
            .phosphorus_to_carbon_ratio = &phosphorus,
        }, .{
            .nitrogen_to_carbon_ratio = &.{ 0.1, 0.2, 0.3, 0.4 },
            .phosphorus_to_carbon_ratio = &.{ 0.01, 0.02, 0.03, std.math.nan(f64) },
            .residue_class_count = 2,
            .biochemical_fraction_count = 2,
        }),
    );

    try std.testing.expectEqualSlices(f64, &nitrogen_before, &nitrogen);
    try std.testing.expectEqualSlices(f64, &phosphorus_before, &phosphorus);
}
