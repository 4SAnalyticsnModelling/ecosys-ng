const std = @import("std");

pub const PlantActivity = enum {
    inactive,
    active,
};

pub const MortalityState = enum {
    alive,
    dead,
};

pub const Result = enum {
    skipped_inactive,
    skipped_fully_dead,
    reset,
};

pub const State = struct {
    living_leaf_area_m2_by_canopy_layer: []f64,
    living_leaf_carbon_g_c_by_canopy_layer: []f64,
    living_stalk_area_m2_by_canopy_layer: []f64,
};

pub const Conditions = struct {
    activity: PlantActivity,
    shoot_mortality: MortalityState,
    root_mortality: MortalityState,
};

/// Exact GROSUB lines 359--365 guarded hourly canopy-work reset.
///
/// The source executes for an active plant while either its shoot or root is
/// alive. All canopy-layer slices are caller allocated at runtime.
pub fn resetIfRequired(state: State, conditions: Conditions) !Result {
    if (conditions.activity != .active) return .skipped_inactive;
    if (conditions.shoot_mortality == .dead and
        conditions.root_mortality == .dead)
        return .skipped_fully_dead;
    const layer_count = state.living_leaf_area_m2_by_canopy_layer.len;
    if (layer_count == 0 or
        state.living_leaf_carbon_g_c_by_canopy_layer.len != layer_count or
        state.living_stalk_area_m2_by_canopy_layer.len != layer_count)
        return error.CanopyWorkspaceLayerDimensionMismatch;

    @memset(state.living_leaf_area_m2_by_canopy_layer, 0);
    @memset(state.living_leaf_carbon_g_c_by_canopy_layer, 0);
    @memset(state.living_stalk_area_m2_by_canopy_layer, 0);
    return .reset;
}

test "GROSUB resets runtime canopy layers when either organ system lives" {
    const allocator = std.testing.allocator;
    const layer_count = 29;
    const leaf_area = try allocator.alloc(f64, layer_count);
    defer allocator.free(leaf_area);
    const leaf_carbon = try allocator.alloc(f64, layer_count);
    defer allocator.free(leaf_carbon);
    const stalk_area = try allocator.alloc(f64, layer_count);
    defer allocator.free(stalk_area);
    @memset(leaf_area, 1);
    @memset(leaf_carbon, 2);
    @memset(stalk_area, 3);

    const result = try resetIfRequired(.{
        .living_leaf_area_m2_by_canopy_layer = leaf_area,
        .living_leaf_carbon_g_c_by_canopy_layer = leaf_carbon,
        .living_stalk_area_m2_by_canopy_layer = stalk_area,
    }, .{
        .activity = .active,
        .shoot_mortality = .dead,
        .root_mortality = .alive,
    });

    try std.testing.expectEqual(Result.reset, result);
    for (leaf_area, leaf_carbon, stalk_area) |area, carbon, stalk|
        inline for (.{ area, carbon, stalk }) |value|
            try std.testing.expectEqual(@as(f64, 0), value);
}

test "inactive and fully dead guards preserve canopy workspace" {
    inline for (.{
        Conditions{ .activity = .inactive, .shoot_mortality = .alive, .root_mortality = .alive },
        Conditions{ .activity = .active, .shoot_mortality = .dead, .root_mortality = .dead },
    }) |conditions| {
        var leaf_area = [_]f64{1};
        var leaf_carbon = [_]f64{2};
        var stalk_area = [_]f64{3};
        const result = try resetIfRequired(.{
            .living_leaf_area_m2_by_canopy_layer = &leaf_area,
            .living_leaf_carbon_g_c_by_canopy_layer = &leaf_carbon,
            .living_stalk_area_m2_by_canopy_layer = &stalk_area,
        }, conditions);
        try std.testing.expect(result != .reset);
        try std.testing.expectEqual(@as(f64, 1), leaf_area[0]);
        try std.testing.expectEqual(@as(f64, 2), leaf_carbon[0]);
        try std.testing.expectEqual(@as(f64, 3), stalk_area[0]);
    }
}

test "dimension failure occurs before any canopy mutation" {
    var leaf_area = [_]f64{ 1, 2 };
    var leaf_carbon = [_]f64{3};
    var stalk_area = [_]f64{ 4, 5 };
    const before = leaf_area;
    try std.testing.expectError(
        error.CanopyWorkspaceLayerDimensionMismatch,
        resetIfRequired(.{
            .living_leaf_area_m2_by_canopy_layer = &leaf_area,
            .living_leaf_carbon_g_c_by_canopy_layer = &leaf_carbon,
            .living_stalk_area_m2_by_canopy_layer = &stalk_area,
        }, .{
            .activity = .active,
            .shoot_mortality = .alive,
            .root_mortality = .dead,
        }),
    );
    try std.testing.expectEqualSlices(f64, &before, &leaf_area);
}
