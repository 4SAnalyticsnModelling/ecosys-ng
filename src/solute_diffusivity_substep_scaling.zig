const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const diffusivity_class_count = @typeInfo(species_registry.DiffusivityClass).@"enum".fields.len;

/// Exact compatibility translation of TRNSFRS.F lines 1477--1490.
///
/// Input diffusivities are m2 h-1 and outputs are m2 substep-1. Class order is
/// PO4, Al, Fe, H, Ca, Mg, Na, K, OH, SO4, Cl, CO3, HCO3, H4SiO4.
pub fn scale(
    hourly_diffusivity_m2_per_h: []const f64,
    substep_diffusivity_m2: []f64,
    substep_duration_h: f64,
) !void {
    if (hourly_diffusivity_m2_per_h.len != diffusivity_class_count or
        substep_diffusivity_m2.len != diffusivity_class_count)
        return error.SoluteDiffusivityClassDimensionMismatch;
    if (!std.math.isFinite(substep_duration_h) or substep_duration_h < 0 or substep_duration_h > 1)
        return error.InvalidSoluteSubstepDuration;

    for (hourly_diffusivity_m2_per_h) |diffusivity| {
        if (!std.math.isFinite(diffusivity) or diffusivity < 0)
            return error.InvalidHourlySoluteDiffusivity;
        if (!std.math.isFinite(diffusivity * substep_duration_h))
            return error.NonFiniteSubstepSoluteDiffusivity;
    }
    for (hourly_diffusivity_m2_per_h, substep_diffusivity_m2) |diffusivity, *scaled|
        scaled.* = diffusivity * substep_duration_h;
}

test "TRNSFRS diffusivity scaling preserves all 14 source classes" {
    var hourly: [diffusivity_class_count]f64 = undefined;
    var substep = [_]f64{-1} ** diffusivity_class_count;
    for (&hourly, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    try scale(&hourly, &substep, 0.25);
    for (hourly, substep) |hourly_value, scaled_value|
        try std.testing.expectEqual(hourly_value * 0.25, scaled_value);
}

test "TRNSFRS diffusivity class endpoints retain phosphate and silicate order" {
    var hourly = [_]f64{0} ** diffusivity_class_count;
    var substep = [_]f64{-1} ** diffusivity_class_count;
    hourly[@intFromEnum(species_registry.DiffusivityClass.phosphate)] = 2;
    hourly[@intFromEnum(species_registry.DiffusivityClass.hydrogen_silicate)] = 6;
    try scale(&hourly, &substep, 0.5);
    try std.testing.expectEqual(@as(f64, 1), substep[0]);
    try std.testing.expectEqual(@as(f64, 3), substep[diffusivity_class_count - 1]);
}

test "solute diffusivity late failure is atomic" {
    var hourly = [_]f64{2} ** diffusivity_class_count;
    var substep = [_]f64{7} ** diffusivity_class_count;
    hourly[diffusivity_class_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidHourlySoluteDiffusivity,
        scale(&hourly, &substep, 0.5),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{7} ** diffusivity_class_count), &substep);
}

test "solute diffusivity topology is exact" {
    const hourly = [_]f64{2} ** (diffusivity_class_count - 1);
    var substep = [_]f64{7} ** diffusivity_class_count;
    try std.testing.expectError(
        error.SoluteDiffusivityClassDimensionMismatch,
        scale(&hourly, &substep, 0.5),
    );
}
