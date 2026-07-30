const std = @import("std");

pub const CombustionHeat = struct {
    living_canopy_mj_per_timestep: []f64,
    standing_dead_mj_per_timestep: []f64,
};

pub const InitializationError = error{
    SpeciesCountMismatch,
};

/// Translates STARTQ lines 912-915 for runtime active plant species.
pub fn initialize(heat: CombustionHeat) InitializationError!void {
    if (heat.living_canopy_mj_per_timestep.len !=
        heat.standing_dead_mj_per_timestep.len)
    {
        return error.SpeciesCountMismatch;
    }

    for (0..heat.living_canopy_mj_per_timestep.len) |species| {
        heat.living_canopy_mj_per_timestep[species] = 0.0;
        heat.standing_dead_mj_per_timestep[species] = 0.0;
    }
}

test "runtime active species combustion heat resets in STARTQ order" {
    var living_canopy = [_]f64{ 1.0, 2.0, 3.0 };
    var standing_dead = [_]f64{ 4.0, 5.0, 6.0 };

    try initialize(.{
        .living_canopy_mj_per_timestep = &living_canopy,
        .standing_dead_mj_per_timestep = &standing_dead,
    });

    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0 }, &living_canopy);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0 }, &standing_dead);
}

test "species mismatch fails before either heat ledger mutates" {
    var living_canopy = [_]f64{ 1.0, 2.0 };
    var standing_dead = [_]f64{3.0};

    try std.testing.expectError(error.SpeciesCountMismatch, initialize(.{
        .living_canopy_mj_per_timestep = &living_canopy,
        .standing_dead_mj_per_timestep = &standing_dead,
    }));
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0 }, &living_canopy);
    try std.testing.expectEqualSlices(f64, &.{3.0}, &standing_dead);
}
