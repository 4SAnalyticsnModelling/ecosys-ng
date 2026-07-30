const std = @import("std");

pub const Inputs = struct {
    solar_angle_sine: f64,
    leaf_area_m2: f64,
    leaf_area_presence_threshold_m2: f64,
};

/// STOMATE.F 57--58. Photosynthesis is admitted only when both source
/// comparisons are strictly true, in solar-angle then leaf-area order.
pub fn compute(inputs: Inputs) !bool {
    if (!std.math.isFinite(inputs.solar_angle_sine) or
        !std.math.isFinite(inputs.leaf_area_m2) or
        !std.math.isFinite(inputs.leaf_area_presence_threshold_m2))
        return error.NonFiniteCanopyPhotosynthesisAdmissionInput;
    if (inputs.solar_angle_sine < -1 or inputs.solar_angle_sine > 1 or
        inputs.leaf_area_m2 < 0 or inputs.leaf_area_presence_threshold_m2 < 0)
        return error.InvalidCanopyPhotosynthesisAdmissionInput;
    return inputs.solar_angle_sine > 0 and
        inputs.leaf_area_m2 > inputs.leaf_area_presence_threshold_m2;
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    scratch: []bool,
    destination: []bool,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.CanopyPhotosynthesisAdmissionDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs);
    @memcpy(destination, scratch);
}

test "positive sun and leaf area above threshold admit photosynthesis" {
    try std.testing.expect(try compute(.{
        .solar_angle_sine = 0.5,
        .leaf_area_m2 = 2,
        .leaf_area_presence_threshold_m2 = 1e-12,
    }));
}

test "source comparisons remain strict at both boundaries" {
    try std.testing.expect(!try compute(.{
        .solar_angle_sine = 0,
        .leaf_area_m2 = 2,
        .leaf_area_presence_threshold_m2 = 1,
    }));
    try std.testing.expect(!try compute(.{
        .solar_angle_sine = 0.5,
        .leaf_area_m2 = 1,
        .leaf_area_presence_threshold_m2 = 1,
    }));
}

test "runtime axes support arbitrary cell and species dimensions" {
    const inputs = [_]Inputs{
        .{
            .solar_angle_sine = 0.5,
            .leaf_area_m2 = 2,
            .leaf_area_presence_threshold_m2 = 1,
        },
        .{
            .solar_angle_sine = -0.5,
            .leaf_area_m2 = 2,
            .leaf_area_presence_threshold_m2 = 1,
        },
        .{
            .solar_angle_sine = 0.5,
            .leaf_area_m2 = 0,
            .leaf_area_presence_threshold_m2 = 0,
        },
    };
    var scratch: [3]bool = undefined;
    var destination: [3]bool = undefined;
    try computeRuntimeAxes(&inputs, &scratch, &destination);
    try std.testing.expectEqualSlices(bool, &[_]bool{ true, false, false }, &destination);
}

test "later invalid axis leaves destination unchanged" {
    const inputs = [_]Inputs{
        .{
            .solar_angle_sine = 0.5,
            .leaf_area_m2 = 2,
            .leaf_area_presence_threshold_m2 = 1,
        },
        .{
            .solar_angle_sine = std.math.nan(f64),
            .leaf_area_m2 = 2,
            .leaf_area_presence_threshold_m2 = 1,
        },
    };
    var scratch: [2]bool = undefined;
    var destination = [_]bool{ false, true };
    try std.testing.expectError(
        error.NonFiniteCanopyPhotosynthesisAdmissionInput,
        computeRuntimeAxes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqualSlices(bool, &[_]bool{ false, true }, &destination);
}
