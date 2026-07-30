const std = @import("std");

pub const OrganicFluxes = struct {
    carbon_g_timestep: []f64, // XOCQRS
    nitrogen_g_timestep: []f64, // XONQRS
    phosphorus_g_timestep: []f64, // XOPQRS
    acetate_g_timestep: []f64, // XOAQRS
};

pub const RunoffWater = struct {
    runoff_m3_timestep: f64, // QR
    runoff_heat_j_timestep: f64, // HQR
    snow_m3_timestep: f64, // QS
    water_m3_timestep: f64, // QW
    ice_m3_timestep: f64, // QI
    snow_heat_j_timestep: f64, // HQS
};

pub const RunoffGases = struct {
    carbon_dioxide_g_timestep: f64, // XCOQRS
    methane_g_timestep: f64, // XCHQRS
    oxygen_g_timestep: f64, // XOXQRS
    nitrogen_g_timestep: f64, // XNGQRS
    nitrous_oxide_g_timestep: f64, // XN2QRS
    hydrogen_g_timestep: f64, // XHGQRS
};

pub const RunoffSolutes = struct {
    ammonium_g_n_timestep: f64, // XN4QRW
    ammonia_g_n_timestep: f64, // XN3QRW
    nitrate_g_n_timestep: f64, // XNOQRW
    other_nitrogen_g_n_timestep: f64, // XNXQRS
    phosphate_h2po4_g_p_timestep: f64, // XP1QRW
    phosphate_hpo4_g_p_timestep: f64, // XP4QRW
};

pub const SnowSolutes = struct {
    carbon_dioxide_g_timestep: f64, // XCOQSS
    methane_g_timestep: f64, // XCHQSS
    oxygen_g_timestep: f64, // XOXQSS
    nitrogen_g_timestep: f64, // XNGQSS
    nitrous_oxide_g_timestep: f64, // XN2QSS
    ammonium_g_n_timestep: f64, // XN4QSS
    ammonia_g_n_timestep: f64, // XN3QSS
    nitrate_g_n_timestep: f64, // XNOQSS
    phosphate_h2po4_g_p_timestep: f64, // XP1QSS
    phosphate_hpo4_g_p_timestep: f64, // XP4QSS
};

pub const SurfaceBoundaryAccumulators = struct {
    water: RunoffWater,
    organic: OrganicFluxes,
    gases: RunoffGases,
    runoff_solutes: RunoffSolutes,
    snow_solutes: SnowSolutes,
};

pub const InitializationError = error{OrganicPoolLengthMismatch};

/// Translates the surface runoff portion of HOUR1 lines 2499-2535.
/// Both lateral direction and boundary-side extents are runtime allocated.
pub fn resetSurfaceRunoff(
    boundaries: []SurfaceBoundaryAccumulators,
    organic_pool_count: usize,
) InitializationError!void {
    for (boundaries) |*boundary| {
        if (boundary.organic.carbon_g_timestep.len != organic_pool_count or
            boundary.organic.nitrogen_g_timestep.len != organic_pool_count or
            boundary.organic.phosphorus_g_timestep.len != organic_pool_count or
            boundary.organic.acetate_g_timestep.len != organic_pool_count)
        {
            return error.OrganicPoolLengthMismatch;
        }

        boundary.water.runoff_m3_timestep = 0.0;
        boundary.water.runoff_heat_j_timestep = 0.0;
        for (0..organic_pool_count) |pool_index| {
            boundary.organic.carbon_g_timestep[pool_index] = 0.0;
            boundary.organic.nitrogen_g_timestep[pool_index] = 0.0;
            boundary.organic.phosphorus_g_timestep[pool_index] = 0.0;
            boundary.organic.acetate_g_timestep[pool_index] = 0.0;
        }
        inline for (std.meta.fields(RunoffGases)) |field| {
            @field(boundary.gases, field.name) = 0.0;
        }
        inline for (std.meta.fields(RunoffSolutes)) |field| {
            @field(boundary.runoff_solutes, field.name) = 0.0;
        }
        boundary.water.snow_m3_timestep = 0.0;
        boundary.water.water_m3_timestep = 0.0;
        boundary.water.ice_m3_timestep = 0.0;
        boundary.water.snow_heat_j_timestep = 0.0;
        inline for (std.meta.fields(SnowSolutes)) |field| {
            @field(boundary.snow_solutes, field.name) = 0.0;
        }
    }
}

test "surface runoff reset covers runtime boundary and organic pool extents" {
    const allocator = std.testing.allocator;
    var organic_storage = try allocator.alloc(f64, 4 * 7);
    defer allocator.free(organic_storage);
    @memset(organic_storage, 9.0);

    const boundaries = try allocator.alloc(SurfaceBoundaryAccumulators, 7);
    defer allocator.free(boundaries);
    for (boundaries, 0..) |*boundary, index| {
        const offset = index * 4;
        boundary.* = .{
            .water = std.mem.zeroes(RunoffWater),
            .organic = .{
                .carbon_g_timestep = organic_storage[offset .. offset + 1],
                .nitrogen_g_timestep = organic_storage[offset + 1 .. offset + 2],
                .phosphorus_g_timestep = organic_storage[offset + 2 .. offset + 3],
                .acetate_g_timestep = organic_storage[offset + 3 .. offset + 4],
            },
            .gases = std.mem.zeroes(RunoffGases),
            .runoff_solutes = std.mem.zeroes(RunoffSolutes),
            .snow_solutes = std.mem.zeroes(SnowSolutes),
        };
        inline for (std.meta.fields(RunoffWater)) |field| @field(boundary.water, field.name) = 9.0;
        inline for (std.meta.fields(RunoffGases)) |field| @field(boundary.gases, field.name) = 9.0;
        inline for (std.meta.fields(RunoffSolutes)) |field| @field(boundary.runoff_solutes, field.name) = 9.0;
        inline for (std.meta.fields(SnowSolutes)) |field| @field(boundary.snow_solutes, field.name) = 9.0;
    }

    try resetSurfaceRunoff(boundaries, 1);

    for (organic_storage) |value| try std.testing.expectEqual(@as(f64, 0.0), value);
    for (boundaries) |boundary| {
        inline for (std.meta.fields(RunoffWater)) |field| {
            try std.testing.expectEqual(@as(f64, 0.0), @field(boundary.water, field.name));
        }
        inline for (std.meta.fields(RunoffGases)) |field| {
            try std.testing.expectEqual(@as(f64, 0.0), @field(boundary.gases, field.name));
        }
        inline for (std.meta.fields(RunoffSolutes)) |field| {
            try std.testing.expectEqual(@as(f64, 0.0), @field(boundary.runoff_solutes, field.name));
        }
        inline for (std.meta.fields(SnowSolutes)) |field| {
            try std.testing.expectEqual(@as(f64, 0.0), @field(boundary.snow_solutes, field.name));
        }
    }
}
