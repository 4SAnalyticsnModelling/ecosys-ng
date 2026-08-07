const std = @import("std");

pub const UreaRelease = enum { fast, normal, slow };

pub const ActivityProfiles = struct {
    hydrolysis_initial: []f64,
    hydrolysis_current: []f64,
    nitrification_initial: []f64,
    nitrification_current: []f64,
};

pub const Inputs = struct {
    urea_fertilizer_g_n: f64,
    additional_urea_fertilizer_g_n: f64,
    formulation_code: i32,
    placement_layer: usize,
};

/// `hour1.f` lines 921--951. Initializes runtime layer profiles for urea
/// hydrolysis and nitrification inhibition without fixed NL dimensions.
pub fn apply(
    inputs: Inputs,
    urea_release: *UreaRelease,
    profiles: ActivityProfiles,
) !void {
    try validate(inputs, profiles);
    if (inputs.urea_fertilizer_g_n > 0 or inputs.additional_urea_fertilizer_g_n > 0) {
        urea_release.* = switch (inputs.formulation_code) {
            0 => .fast,
            1, 3 => .normal,
            else => .slow,
        };
        for (profiles.hydrolysis_initial, profiles.hydrolysis_current, 0..) |
            *initial,
            *current,
            layer,
        | {
            const activity: f64 = if (layer == inputs.placement_layer) 1 else 0;
            initial.* = activity;
            current.* = activity;
        }
    }
    if (inputs.formulation_code == 3 or inputs.formulation_code == 4) {
        for (profiles.nitrification_initial, profiles.nitrification_current, 0..) |
            *initial,
            *current,
            layer,
        | {
            const activity: f64 = if (layer == inputs.placement_layer) 1 else 0;
            initial.* = activity;
            current.* = activity;
        }
    }
}

fn validate(inputs: Inputs, profiles: ActivityProfiles) !void {
    if (!std.math.isFinite(inputs.urea_fertilizer_g_n) or
        !std.math.isFinite(inputs.additional_urea_fertilizer_g_n))
        return error.NonFiniteUreaFertilizerInput;
    if (inputs.urea_fertilizer_g_n < 0 or inputs.additional_urea_fertilizer_g_n < 0)
        return error.InvalidUreaFertilizerInput;
    const layer_count = profiles.hydrolysis_initial.len;
    if (layer_count == 0 or
        profiles.hydrolysis_current.len != layer_count or
        profiles.nitrification_initial.len != layer_count or
        profiles.nitrification_current.len != layer_count)
        return error.FertilizerInhibitorLayerDimensionMismatch;
    if (inputs.placement_layer >= layer_count)
        return error.FertilizerPlacementLayerOutOfRange;
}

test "urea formulation and hydrolysis activity follow source mapping" {
    var h0 = [_]f64{ 9, 9, 9 };
    var hi = h0;
    var n0 = h0;
    var ni = h0;
    var release: UreaRelease = .fast;
    try apply(.{
        .urea_fertilizer_g_n = 1,
        .additional_urea_fertilizer_g_n = 0,
        .formulation_code = 3,
        .placement_layer = 1,
    }, &release, .{
        .hydrolysis_initial = &h0,
        .hydrolysis_current = &hi,
        .nitrification_initial = &n0,
        .nitrification_current = &ni,
    });
    try std.testing.expectEqual(UreaRelease.normal, release);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1, 0 }, &h0);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1, 0 }, &hi);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1, 0 }, &n0);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1, 0 }, &ni);
}

test "nitrification inhibitor gate is independent of urea presence" {
    var h0 = [_]f64{ 8, 8 };
    var hi = h0;
    var n0 = [_]f64{ 9, 9 };
    var ni = n0;
    var release: UreaRelease = .normal;
    try apply(.{
        .urea_fertilizer_g_n = 0,
        .additional_urea_fertilizer_g_n = 0,
        .formulation_code = 4,
        .placement_layer = 0,
    }, &release, .{
        .hydrolysis_initial = &h0,
        .hydrolysis_current = &hi,
        .nitrification_initial = &n0,
        .nitrification_current = &ni,
    });
    try std.testing.expectEqual(UreaRelease.normal, release);
    try std.testing.expectEqualSlices(f64, &.{ 8, 8 }, &h0);
    try std.testing.expectEqualSlices(f64, &.{ 1, 0 }, &n0);
}

test "inactive gates preserve all caller state" {
    var h0 = [_]f64{ 1, 2 };
    var hi = [_]f64{ 3, 4 };
    var n0 = [_]f64{ 5, 6 };
    var ni = [_]f64{ 7, 8 };
    var release: UreaRelease = .slow;
    try apply(.{
        .urea_fertilizer_g_n = 0,
        .additional_urea_fertilizer_g_n = 0,
        .formulation_code = 2,
        .placement_layer = 1,
    }, &release, .{
        .hydrolysis_initial = &h0,
        .hydrolysis_current = &hi,
        .nitrification_initial = &n0,
        .nitrification_current = &ni,
    });
    try std.testing.expectEqual(UreaRelease.slow, release);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2 }, &h0);
    try std.testing.expectEqualSlices(f64, &.{ 5, 6 }, &n0);
}
