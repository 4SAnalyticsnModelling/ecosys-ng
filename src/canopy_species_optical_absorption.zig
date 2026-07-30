const std = @import("std");

pub const Inputs = struct {
    direct_shortwave_mj_m2_h: f64,
    direct_par_umol_m2_s: f64,
    leaf_shortwave_absorptivity: []const f64,
    leaf_par_absorptivity: []const f64,
    stalk_shortwave_absorptivity: f64,
    standing_dead_shortwave_absorptivity: f64,
    stalk_par_absorptivity: f64,
    standing_dead_par_absorptivity: f64,
};

pub const Outputs = struct {
    leaf_shortwave_mj_m2_h: []f64,
    stalk_shortwave_mj_m2_h: []f64,
    standing_dead_shortwave_mj_m2_h: []f64,
    leaf_par_umol_m2_s: []f64,
    stalk_par_umol_m2_s: []f64,
    standing_dead_par_umol_m2_s: []f64,
};

/// HOUR1 lines 1059--1066. Computes radiation normal to leaf, stalk, and
/// standing-dead surfaces for an arbitrary runtime species count.
pub fn apply(inputs: Inputs, outputs: Outputs) !void {
    try validate(inputs, outputs);
    for (inputs.leaf_shortwave_absorptivity, inputs.leaf_par_absorptivity, 0..) |
        leaf_shortwave_absorptivity,
        leaf_par_absorptivity,
        species,
    | {
        outputs.leaf_shortwave_mj_m2_h[species] =
            inputs.direct_shortwave_mj_m2_h * leaf_shortwave_absorptivity;
        outputs.stalk_shortwave_mj_m2_h[species] =
            inputs.direct_shortwave_mj_m2_h * inputs.stalk_shortwave_absorptivity;
        outputs.standing_dead_shortwave_mj_m2_h[species] =
            inputs.direct_shortwave_mj_m2_h *
            inputs.standing_dead_shortwave_absorptivity;
        outputs.leaf_par_umol_m2_s[species] =
            inputs.direct_par_umol_m2_s * leaf_par_absorptivity;
        outputs.stalk_par_umol_m2_s[species] =
            inputs.direct_par_umol_m2_s * inputs.stalk_par_absorptivity;
        outputs.standing_dead_par_umol_m2_s[species] =
            inputs.direct_par_umol_m2_s * inputs.standing_dead_par_absorptivity;
    }
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    const species_count = inputs.leaf_shortwave_absorptivity.len;
    if (species_count == 0 or
        inputs.leaf_par_absorptivity.len != species_count or
        outputs.leaf_shortwave_mj_m2_h.len != species_count or
        outputs.stalk_shortwave_mj_m2_h.len != species_count or
        outputs.standing_dead_shortwave_mj_m2_h.len != species_count or
        outputs.leaf_par_umol_m2_s.len != species_count or
        outputs.stalk_par_umol_m2_s.len != species_count or
        outputs.standing_dead_par_umol_m2_s.len != species_count)
        return error.CanopyOpticalAbsorptionDimensionMismatch;
    if (!std.math.isFinite(inputs.direct_shortwave_mj_m2_h) or
        !std.math.isFinite(inputs.direct_par_umol_m2_s))
        return error.NonFiniteCanopyOpticalAbsorptionInput;
    if (inputs.direct_shortwave_mj_m2_h < 0 or inputs.direct_par_umol_m2_s < 0)
        return error.InvalidCanopyOpticalAbsorptionInput;
    inline for (.{
        inputs.leaf_shortwave_absorptivity,
        inputs.leaf_par_absorptivity,
    }) |values| for (values) |value| try validateCoefficient(value);
    inline for (.{
        inputs.stalk_shortwave_absorptivity,
        inputs.standing_dead_shortwave_absorptivity,
        inputs.stalk_par_absorptivity,
        inputs.standing_dead_par_absorptivity,
    }) |value| try validateCoefficient(value);
}

fn validateCoefficient(value: f64) !void {
    if (!std.math.isFinite(value))
        return error.NonFiniteCanopyOpticalAbsorptionInput;
    if (value < 0 or value > 1)
        return error.InvalidCanopyOpticalAbsorptionInput;
}

test "runtime species optical absorption preserves six assignment order" {
    var leaf_sw: [3]f64 = undefined;
    var stalk_sw: [3]f64 = undefined;
    var dead_sw: [3]f64 = undefined;
    var leaf_par: [3]f64 = undefined;
    var stalk_par: [3]f64 = undefined;
    var dead_par: [3]f64 = undefined;
    try apply(.{
        .direct_shortwave_mj_m2_h = 2,
        .direct_par_umol_m2_s = 1000,
        .leaf_shortwave_absorptivity = &.{ 0.5, 0.6, 0.7 },
        .leaf_par_absorptivity = &.{ 0.8, 0.9, 1.0 },
        .stalk_shortwave_absorptivity = 0.4,
        .standing_dead_shortwave_absorptivity = 0.3,
        .stalk_par_absorptivity = 0.2,
        .standing_dead_par_absorptivity = 0.1,
    }, .{
        .leaf_shortwave_mj_m2_h = &leaf_sw,
        .stalk_shortwave_mj_m2_h = &stalk_sw,
        .standing_dead_shortwave_mj_m2_h = &dead_sw,
        .leaf_par_umol_m2_s = &leaf_par,
        .stalk_par_umol_m2_s = &stalk_par,
        .standing_dead_par_umol_m2_s = &dead_par,
    });
    try std.testing.expectEqualSlices(f64, &.{ 1, 1.2, 1.4 }, &leaf_sw);
    try std.testing.expectEqualSlices(f64, &.{ 0.8, 0.8, 0.8 }, &stalk_sw);
    try std.testing.expectEqualSlices(f64, &.{ 800, 900, 1000 }, &leaf_par);
    try std.testing.expectEqualSlices(f64, &.{ 100, 100, 100 }, &dead_par);
}

test "dimension mismatch leaves caller output unchanged" {
    var one = [_]f64{42};
    try std.testing.expectError(error.CanopyOpticalAbsorptionDimensionMismatch, apply(.{
        .direct_shortwave_mj_m2_h = 1,
        .direct_par_umol_m2_s = 1,
        .leaf_shortwave_absorptivity = &.{ 0.5, 0.5 },
        .leaf_par_absorptivity = &.{ 0.5, 0.5 },
        .stalk_shortwave_absorptivity = 0.5,
        .standing_dead_shortwave_absorptivity = 0.5,
        .stalk_par_absorptivity = 0.5,
        .standing_dead_par_absorptivity = 0.5,
    }, .{
        .leaf_shortwave_mj_m2_h = &one,
        .stalk_shortwave_mj_m2_h = &one,
        .standing_dead_shortwave_mj_m2_h = &one,
        .leaf_par_umol_m2_s = &one,
        .stalk_par_umol_m2_s = &one,
        .standing_dead_par_umol_m2_s = &one,
    }));
    try std.testing.expectEqual(@as(f64, 42), one[0]);
}
