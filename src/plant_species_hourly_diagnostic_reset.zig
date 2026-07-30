const std = @import("std");

/// Runtime-owned diagnostics for one user-defined plant species.
pub const SpeciesDiagnostics = struct {
    living_canopy_height_m: f64, // ZC
    standing_dead_height_m: f64, // ZG
    net_carbon_exchange_g_c_timestep: f64, // CNET
    living_canopy_combustion_g_c_timestep: f64, // RCGCKZ
    standing_dead_combustion_g_c_timestep: f64, // RCGDKZ
    canopy_snow_carbon_g: f64, // HCSNC
    canopy_snow_nitrogen_g: f64, // HZSNC
    canopy_snow_phosphorus_g: f64, // HPSNC
    previous_living_combustion_heat_mj_timestep: f64, // HCBFCY
    current_living_combustion_heat_mj_timestep: f64, // HCBFCZ
    previous_dead_combustion_heat_mj_timestep: f64, // HCBFDY
    current_dead_combustion_heat_mj_timestep: f64, // HCBFDZ
    root_combustion_g_c_timestep: []f64, // RCGSKR over runtime rooted layers
};

pub const ResetError = error{
    NonFiniteCanopyHeight,
    NegativeCanopyHeight,
    NonFiniteCombustionHeat,
};

/// Translates HOUR1 lines 3065-3079. Species count and rooted-layer count
/// are both determined by the runtime slices.
pub fn reset(species: []SpeciesDiagnostics) ResetError!f64 {
    for (species) |item| {
        if (!std.math.isFinite(item.living_canopy_height_m) or
            !std.math.isFinite(item.standing_dead_height_m))
        {
            return error.NonFiniteCanopyHeight;
        }
        if (item.living_canopy_height_m < 0.0 or item.standing_dead_height_m < 0.0) {
            return error.NegativeCanopyHeight;
        }
        if (!std.math.isFinite(item.current_living_combustion_heat_mj_timestep) or
            !std.math.isFinite(item.current_dead_combustion_heat_mj_timestep))
        {
            return error.NonFiniteCombustionHeat;
        }
    }

    var maximum_canopy_height_m: f64 = 0.0;
    for (species) |*item| {
        maximum_canopy_height_m = @max(
            maximum_canopy_height_m,
            item.living_canopy_height_m,
            item.standing_dead_height_m,
        );
        item.net_carbon_exchange_g_c_timestep = 0.0;
        item.living_canopy_combustion_g_c_timestep = 0.0;
        item.standing_dead_combustion_g_c_timestep = 0.0;
        item.canopy_snow_carbon_g = 0.0;
        item.canopy_snow_nitrogen_g = 0.0;
        item.canopy_snow_phosphorus_g = 0.0;
        item.previous_living_combustion_heat_mj_timestep =
            item.current_living_combustion_heat_mj_timestep;
        item.current_living_combustion_heat_mj_timestep = 0.0;
        item.previous_dead_combustion_heat_mj_timestep =
            item.current_dead_combustion_heat_mj_timestep;
        item.current_dead_combustion_heat_mj_timestep = 0.0;
        for (item.root_combustion_g_c_timestep) |*root_combustion| {
            root_combustion.* = 0.0;
        }
    }
    return maximum_canopy_height_m;
}

test "runtime species and root layers reset with heat-state rollover" {
    const allocator = std.testing.allocator;
    const species = try allocator.alloc(SpeciesDiagnostics, 8);
    defer allocator.free(species);
    for (species, 0..) |*item, index| {
        item.* = .{
            .living_canopy_height_m = @as(f64, @floatFromInt(index)) + 0.5,
            .standing_dead_height_m = @as(f64, @floatFromInt(index)) + 1.0,
            .net_carbon_exchange_g_c_timestep = 9.0,
            .living_canopy_combustion_g_c_timestep = 9.0,
            .standing_dead_combustion_g_c_timestep = 9.0,
            .canopy_snow_carbon_g = 9.0,
            .canopy_snow_nitrogen_g = 9.0,
            .canopy_snow_phosphorus_g = 9.0,
            .previous_living_combustion_heat_mj_timestep = 1.0,
            .current_living_combustion_heat_mj_timestep = 2.0,
            .previous_dead_combustion_heat_mj_timestep = 3.0,
            .current_dead_combustion_heat_mj_timestep = 4.0,
            .root_combustion_g_c_timestep = try allocator.alloc(f64, index + 1),
        };
        @memset(item.root_combustion_g_c_timestep, 9.0);
    }
    defer for (species) |item| allocator.free(item.root_combustion_g_c_timestep);

    const maximum_height_m = try reset(species);

    try std.testing.expectEqual(@as(f64, 8.0), maximum_height_m);
    for (species) |item| {
        try std.testing.expectEqual(@as(f64, 0.0), item.net_carbon_exchange_g_c_timestep);
        try std.testing.expectEqual(@as(f64, 2.0), item.previous_living_combustion_heat_mj_timestep);
        try std.testing.expectEqual(@as(f64, 0.0), item.current_living_combustion_heat_mj_timestep);
        try std.testing.expectEqual(@as(f64, 4.0), item.previous_dead_combustion_heat_mj_timestep);
        try std.testing.expectEqual(@as(f64, 0.0), item.current_dead_combustion_heat_mj_timestep);
        for (item.root_combustion_g_c_timestep) |value| {
            try std.testing.expectEqual(@as(f64, 0.0), value);
        }
    }
}

test "invalid species height fails before any mutation" {
    var first_roots = [_]f64{7.0};
    var second_roots = [_]f64{8.0};
    var species = [_]SpeciesDiagnostics{
        std.mem.zeroInit(SpeciesDiagnostics, .{
            .living_canopy_height_m = 1.0,
            .standing_dead_height_m = 1.0,
            .net_carbon_exchange_g_c_timestep = 5.0,
            .root_combustion_g_c_timestep = &first_roots,
        }),
        std.mem.zeroInit(SpeciesDiagnostics, .{
            .living_canopy_height_m = -1.0,
            .standing_dead_height_m = 1.0,
            .root_combustion_g_c_timestep = &second_roots,
        }),
    };
    try std.testing.expectError(error.NegativeCanopyHeight, reset(&species));
    try std.testing.expectEqual(@as(f64, 5.0), species[0].net_carbon_exchange_g_c_timestep);
    try std.testing.expectEqual(@as(f64, 7.0), first_roots[0]);
}
