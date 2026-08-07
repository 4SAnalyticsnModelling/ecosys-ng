const std = @import("std");

pub const Parameters = struct {
    ammonium_release_fraction: f64,
};

/// Direct translation of `starte.f` lines 1693--1694. Initial manure N and
/// non-band ammonium are extensive g N per layer. All layers are validated
/// before either inventory is changed.
pub fn apply(
    non_band_ammonium_g_n: []f64,
    initial_manure_nitrogen_g_n: []f64,
    parameters: Parameters,
) !void {
    if (non_band_ammonium_g_n.len != initial_manure_nitrogen_g_n.len or
        non_band_ammonium_g_n.len == 0)
        return error.InitialManureNitrogenDimensionMismatch;
    if (!std.math.isFinite(parameters.ammonium_release_fraction))
        return error.NonFiniteInitialManureNitrogenParameter;
    if (parameters.ammonium_release_fraction < 0 or
        parameters.ammonium_release_fraction > 1)
        return error.InvalidInitialManureNitrogenParameter;
    for (non_band_ammonium_g_n, initial_manure_nitrogen_g_n) |ammonium, manure| {
        inline for (.{ ammonium, manure }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteInitialManureNitrogenInput;
            if (value < 0) return error.InvalidInitialManureNitrogenInput;
        }
        const released = parameters.ammonium_release_fraction * manure;
        const amended_ammonium = ammonium + released;
        const retained_manure = manure - released;
        inline for (.{ released, amended_ammonium, retained_manure }) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteInitialManureNitrogenResult;
    }
    for (non_band_ammonium_g_n, initial_manure_nitrogen_g_n) |*ammonium, *manure| {
        const released = parameters.ammonium_release_fraction * manure.*;
        ammonium.* = ammonium.* + released;
        manure.* = manure.* - released;
    }
}

fn sourceParameters() Parameters {
    return .{ .ammonium_release_fraction = 0.5 };
}

test "STARTE initial manure partition is conservative across runtime layers" {
    var ammonium_g_n = [_]f64{ 1, 2, 3, 4, 5, 6, 7 };
    var manure_g_n = [_]f64{ 2, 4, 6, 8, 10, 12, 14 };
    var before_g_n: f64 = 0;
    for (ammonium_g_n, manure_g_n) |ammonium, manure|
        before_g_n += ammonium + manure;
    try apply(&ammonium_g_n, &manure_g_n, sourceParameters());
    var after_g_n: f64 = 0;
    for (ammonium_g_n, manure_g_n, 0..) |ammonium, manure, layer| {
        const original_manure = @as(f64, @floatFromInt(2 * (layer + 1)));
        try std.testing.expectEqual(original_manure * 0.5, manure);
        after_g_n += ammonium + manure;
    }
    try std.testing.expectEqual(before_g_n, after_g_n);
}

test "STARTE initial manure release fraction is runtime configurable" {
    var ammonium_g_n = [_]f64{1};
    var manure_g_n = [_]f64{8};
    try apply(
        &ammonium_g_n,
        &manure_g_n,
        .{ .ammonium_release_fraction = 0.25 },
    );
    try std.testing.expectEqual(@as(f64, 3), ammonium_g_n[0]);
    try std.testing.expectEqual(@as(f64, 6), manure_g_n[0]);
}

test "STARTE initial manure preflight prevents partial mutation" {
    var ammonium_g_n = [_]f64{ 1, 2 };
    var manure_g_n = [_]f64{ 3, std.math.nan(f64) };
    try std.testing.expectError(
        error.NonFiniteInitialManureNitrogenInput,
        apply(&ammonium_g_n, &manure_g_n, sourceParameters()),
    );
    try std.testing.expectEqualSlices(f64, &.{ 1, 2 }, &ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 3), manure_g_n[0]);
    try std.testing.expect(std.math.isNan(manure_g_n[1]));
}
