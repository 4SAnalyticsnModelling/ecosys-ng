const std = @import("std");

pub const Inputs = struct {
    soil_organic_carbon_g_c_per_Mg: []const f64,
    particulate_organic_carbon_g_c_per_Mg: []const f64,
    dry_soil_mass_Mg: []const f64,
    horizontal_area_m2: f64,
    reference_layer_index: usize,
};

pub const Parameters = struct {
    reference_fraction: f64,
    maximum_reference_humus_g_c_per_m2: f64,
};

pub const Result = struct {
    midpoint_cumulative_humus_g_c_per_m2: []f64,
    reference_humus_g_c_per_m2: *f64,
};

/// Exact source-order translation of legacy `STARTS` lines 760--766.
pub fn derive(
    result: Result,
    inputs: Inputs,
    parameters: Parameters,
) !void {
    const layer_count = inputs.soil_organic_carbon_g_c_per_Mg.len;
    if (layer_count == 0 or
        inputs.particulate_organic_carbon_g_c_per_Mg.len != layer_count or
        inputs.dry_soil_mass_Mg.len != layer_count or
        result.midpoint_cumulative_humus_g_c_per_m2.len != layer_count or
        inputs.reference_layer_index >= layer_count)
        return error.HumusProfileDimensionMismatch;
    inline for (.{
        inputs.horizontal_area_m2,
        parameters.reference_fraction,
        parameters.maximum_reference_humus_g_c_per_m2,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteHumusProfileParameter;
    }
    if (inputs.horizontal_area_m2 <= 0 or
        parameters.reference_fraction < 0 or
        parameters.maximum_reference_humus_g_c_per_m2 < 0)
        return error.InvalidHumusProfileParameter;

    var cumulative_humus_g_c_per_m2: f64 = 0.0;
    var reference_midpoint_g_c_per_m2: f64 = 0.0;
    for (0..layer_count) |layer| {
        inline for (.{
            inputs.soil_organic_carbon_g_c_per_Mg[layer],
            inputs.particulate_organic_carbon_g_c_per_Mg[layer],
            inputs.dry_soil_mass_Mg[layer],
        }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteHumusProfileInput;
            if (value < 0) return error.InvalidHumusProfileInput;
        }
        const humus_concentration_g_c_per_Mg = @max(
            0.0,
            inputs.soil_organic_carbon_g_c_per_Mg[layer] -
                inputs.particulate_organic_carbon_g_c_per_Mg[layer],
        );
        const layer_humus_g_c_per_m2 =
            humus_concentration_g_c_per_Mg *
            inputs.dry_soil_mass_Mg[layer] /
            inputs.horizontal_area_m2;
        const midpoint_g_c_per_m2 =
            cumulative_humus_g_c_per_m2 +
            layer_humus_g_c_per_m2 * 0.5;
        const next_cumulative_g_c_per_m2 =
            cumulative_humus_g_c_per_m2 + layer_humus_g_c_per_m2;
        inline for (.{
            layer_humus_g_c_per_m2,
            midpoint_g_c_per_m2,
            next_cumulative_g_c_per_m2,
        }) |candidate| if (!std.math.isFinite(candidate))
            return error.HumusProfileOverflow;
        if (layer == inputs.reference_layer_index)
            reference_midpoint_g_c_per_m2 = midpoint_g_c_per_m2;
        cumulative_humus_g_c_per_m2 = next_cumulative_g_c_per_m2;
    }
    const reference_humus_g_c_per_m2 = @min(
        parameters.maximum_reference_humus_g_c_per_m2,
        parameters.reference_fraction * reference_midpoint_g_c_per_m2,
    );
    if (!std.math.isFinite(reference_humus_g_c_per_m2))
        return error.HumusProfileOverflow;

    cumulative_humus_g_c_per_m2 = 0.0;
    for (0..layer_count) |layer| {
        const humus_concentration_g_c_per_Mg = @max(
            0.0,
            inputs.soil_organic_carbon_g_c_per_Mg[layer] -
                inputs.particulate_organic_carbon_g_c_per_Mg[layer],
        );
        const layer_humus_g_c_per_m2 =
            humus_concentration_g_c_per_Mg *
            inputs.dry_soil_mass_Mg[layer] /
            inputs.horizontal_area_m2;
        result.midpoint_cumulative_humus_g_c_per_m2[layer] =
            cumulative_humus_g_c_per_m2 +
            layer_humus_g_c_per_m2 * 0.5;
        cumulative_humus_g_c_per_m2 += layer_humus_g_c_per_m2;
    }
    result.reference_humus_g_c_per_m2.* =
        reference_humus_g_c_per_m2;
}

test "STARTS humus profile preserves midpoint cumulative operation order" {
    var midpoint = [_]f64{0.0} ** 3;
    var reference: f64 = 0.0;
    try derive(.{
        .midpoint_cumulative_humus_g_c_per_m2 = &midpoint,
        .reference_humus_g_c_per_m2 = &reference,
    }, .{
        .soil_organic_carbon_g_c_per_Mg = &.{ 100, 80, 40 },
        .particulate_organic_carbon_g_c_per_Mg = &.{ 20, 100, 10 },
        .dry_soil_mass_Mg = &.{ 10, 20, 30 },
        .horizontal_area_m2 = 100,
        .reference_layer_index = 2,
    }, .{
        .reference_fraction = 0.25,
        .maximum_reference_humus_g_c_per_m2 = 5000,
    });

    try std.testing.expectEqualSlices(f64, &.{ 4.0, 8.0, 12.5 }, &midpoint);
    try std.testing.expectEqual(@as(f64, 3.125), reference);
}

test "reference accumulation obeys runtime cap" {
    var midpoint = [_]f64{0.0};
    var reference: f64 = 0.0;
    try derive(.{
        .midpoint_cumulative_humus_g_c_per_m2 = &midpoint,
        .reference_humus_g_c_per_m2 = &reference,
    }, .{
        .soil_organic_carbon_g_c_per_Mg = &.{10_000},
        .particulate_organic_carbon_g_c_per_Mg = &.{0},
        .dry_soil_mass_Mg = &.{1000},
        .horizontal_area_m2 = 1,
        .reference_layer_index = 0,
    }, .{
        .reference_fraction = 0.25,
        .maximum_reference_humus_g_c_per_m2 = 5000,
    });
    try std.testing.expectEqual(@as(f64, 5000), reference);
}

test "late invalid layer preserves profile and reference" {
    var midpoint = [_]f64{ 7, 8 };
    var reference: f64 = 9;
    const before = midpoint;
    try std.testing.expectError(
        error.NonFiniteHumusProfileInput,
        derive(.{
            .midpoint_cumulative_humus_g_c_per_m2 = &midpoint,
            .reference_humus_g_c_per_m2 = &reference,
        }, .{
            .soil_organic_carbon_g_c_per_Mg = &.{ 100, std.math.nan(f64) },
            .particulate_organic_carbon_g_c_per_Mg = &.{ 20, 20 },
            .dry_soil_mass_Mg = &.{ 10, 10 },
            .horizontal_area_m2 = 100,
            .reference_layer_index = 1,
        }, .{
            .reference_fraction = 0.25,
            .maximum_reference_humus_g_c_per_m2 = 5000,
        }),
    );
    try std.testing.expectEqualSlices(f64, &before, &midpoint);
    try std.testing.expectEqual(@as(f64, 9), reference);
}
