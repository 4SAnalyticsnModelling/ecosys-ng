const std = @import("std");

pub const SaltMode = enum { static_equilibrium, dynamic_transport };

/// Source order of HOUR1 XQR* runoff species, lines 2709-2750.
pub const RunoffSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_3,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
    hydrogen_sulfate,
};

/// Source order of HOUR1 XQS* snow species, lines 2751-2792.
pub const SnowSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_1,
    hydrogen_phosphate_3,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
};

/// Fluxes for one runtime-selected lateral direction and boundary side.
/// Values retain the legacy species mass-per-timestep units.
pub const SurfaceSaltFluxes = struct {
    runoff: []f64,
    snow: []f64,
};

pub const ResetError = error{SpeciesCountMismatch};

/// Translates the surface portion of `hour1.f` lines 2704--2794. Static salt
/// equilibrium deliberately leaves the transport arrays untouched.
pub fn reset(
    mode: SaltMode,
    boundaries: []SurfaceSaltFluxes,
) ResetError!void {
    if (mode == .static_equilibrium) return;
    const runoff_count = std.meta.fields(RunoffSpecies).len;
    const snow_count = std.meta.fields(SnowSpecies).len;
    for (boundaries) |boundary| {
        if (boundary.runoff.len != runoff_count or boundary.snow.len != snow_count) {
            return error.SpeciesCountMismatch;
        }
    }
    for (boundaries) |boundary| {
        for (boundary.runoff) |*flux| flux.* = 0.0;
        for (boundary.snow) |*flux| flux.* = 0.0;
    }
}

test "dynamic salt mode resets every runtime boundary and named species" {
    const allocator = std.testing.allocator;
    const runoff_count = std.meta.fields(RunoffSpecies).len;
    const snow_count = std.meta.fields(SnowSpecies).len;
    const boundaries = try allocator.alloc(SurfaceSaltFluxes, 5);
    defer allocator.free(boundaries);
    for (boundaries) |*boundary| {
        boundary.runoff = try allocator.alloc(f64, runoff_count);
        boundary.snow = try allocator.alloc(f64, snow_count);
        @memset(boundary.runoff, 6.0);
        @memset(boundary.snow, 7.0);
    }
    defer for (boundaries) |boundary| {
        allocator.free(boundary.runoff);
        allocator.free(boundary.snow);
    };

    try reset(.dynamic_transport, boundaries);

    for (boundaries) |boundary| {
        for (boundary.runoff) |flux| try std.testing.expectEqual(@as(f64, 0.0), flux);
        for (boundary.snow) |flux| try std.testing.expectEqual(@as(f64, 0.0), flux);
    }
}

test "static salt mode preserves transport arrays" {
    var runoff = [_]f64{9.0};
    var snow = [_]f64{8.0};
    var boundaries = [_]SurfaceSaltFluxes{.{ .runoff = &runoff, .snow = &snow }};
    try reset(.static_equilibrium, &boundaries);
    try std.testing.expectEqual(@as(f64, 9.0), runoff[0]);
    try std.testing.expectEqual(@as(f64, 8.0), snow[0]);
}
