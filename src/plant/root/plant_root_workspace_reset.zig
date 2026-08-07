const std = @import("std");
const canopy_guard = @import("../../canopy/state/live_plant_workspace_reset.zig");

pub const Dimensions = struct {
    biological_domain_count: usize,
    root_axis_count: usize,
    soil_layer_capacity: usize,
    first_active_soil_layer: usize,
    active_soil_layer_count: usize,
};

pub const State = struct {
    root_index_by_domain_axis: []u32,
    root_checked_by_domain_axis: []bool,
    total_root_length_m_by_domain: []f64,
    shoot_root_transfer_g_c_by_domain_layer: []f64,
    primary_root_count_by_domain_layer: []f64,
    living_root_count_by_domain_layer: []f64,
    maintenance_respiration_g_c_by_domain_layer: []f64,
    nitrogen_respiration_g_c_by_domain_layer: []f64,
    actual_respiration_g_c_by_domain_layer: []f64,
    root_length_m_by_domain_layer: []f64,
    primary_sink_m_by_domain_layer_axis: []f64,
    secondary_sink_m_by_domain_layer_axis: []f64,
};

/// Exact grosub.f lines 366--385 root scratch reset under the inherited
/// active-plant and partial-survival guards from lines 359--360.
///
/// Storage is domain-major and then layer/axis major. Only the runtime active
/// soil-layer interval is cleared; inactive capacity remains untouched.
pub fn resetIfRequired(
    state: State,
    dimensions: Dimensions,
    conditions: canopy_guard.Conditions,
) !canopy_guard.Result {
    if (conditions.activity != .active) return .skipped_inactive;
    if (conditions.shoot_mortality == .dead and
        conditions.root_mortality == .dead)
        return .skipped_fully_dead;
    try validate(state, dimensions);

    // Source order: NR outer, N inner.
    for (0..dimensions.root_axis_count) |axis| {
        for (0..dimensions.biological_domain_count) |domain| {
            const index = domain * dimensions.root_axis_count + axis;
            state.root_index_by_domain_axis[index] = 0;
            state.root_checked_by_domain_axis[index] = false;
        }
    }
    // Source order: N outer, then active L, then NR for sink strengths.
    for (0..dimensions.biological_domain_count) |domain| {
        state.total_root_length_m_by_domain[domain] = 0;
        for (0..dimensions.active_soil_layer_count) |local_layer| {
            const layer = dimensions.first_active_soil_layer + local_layer;
            const domain_layer = domain * dimensions.soil_layer_capacity + layer;
            inline for (.{
                state.shoot_root_transfer_g_c_by_domain_layer,
                state.primary_root_count_by_domain_layer,
                state.living_root_count_by_domain_layer,
                state.maintenance_respiration_g_c_by_domain_layer,
                state.nitrogen_respiration_g_c_by_domain_layer,
                state.actual_respiration_g_c_by_domain_layer,
                state.root_length_m_by_domain_layer,
            }) |values| values[domain_layer] = 0;
            for (0..dimensions.root_axis_count) |axis| {
                const index = domain_layer * dimensions.root_axis_count + axis;
                state.primary_sink_m_by_domain_layer_axis[index] = 0;
                state.secondary_sink_m_by_domain_layer_axis[index] = 0;
            }
        }
    }
    return .reset;
}

fn validate(state: State, dimensions: Dimensions) !void {
    if (dimensions.biological_domain_count == 0 or
        dimensions.root_axis_count == 0 or
        dimensions.soil_layer_capacity == 0 or
        dimensions.active_soil_layer_count == 0)
        return error.EmptyPlantRootWorkspace;
    const active_end = std.math.add(
        usize,
        dimensions.first_active_soil_layer,
        dimensions.active_soil_layer_count,
    ) catch return error.PlantRootWorkspaceExtentOverflow;
    if (active_end > dimensions.soil_layer_capacity)
        return error.PlantRootWorkspaceActiveLayerOutOfBounds;
    const domain_axis = try product(
        dimensions.biological_domain_count,
        dimensions.root_axis_count,
    );
    const domain_layer = try product(
        dimensions.biological_domain_count,
        dimensions.soil_layer_capacity,
    );
    const domain_layer_axis = try product(
        domain_layer,
        dimensions.root_axis_count,
    );
    if (state.root_index_by_domain_axis.len != domain_axis or
        state.root_checked_by_domain_axis.len != domain_axis or
        state.total_root_length_m_by_domain.len !=
            dimensions.biological_domain_count)
        return error.PlantRootWorkspaceDimensionMismatch;
    inline for (.{
        state.shoot_root_transfer_g_c_by_domain_layer,
        state.primary_root_count_by_domain_layer,
        state.living_root_count_by_domain_layer,
        state.maintenance_respiration_g_c_by_domain_layer,
        state.nitrogen_respiration_g_c_by_domain_layer,
        state.actual_respiration_g_c_by_domain_layer,
        state.root_length_m_by_domain_layer,
    }) |values| if (values.len != domain_layer)
        return error.PlantRootWorkspaceDimensionMismatch;
    if (state.primary_sink_m_by_domain_layer_axis.len != domain_layer_axis or
        state.secondary_sink_m_by_domain_layer_axis.len != domain_layer_axis)
        return error.PlantRootWorkspaceDimensionMismatch;
}

fn product(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        error.PlantRootWorkspaceExtentOverflow;
}

test "GROSUB clears runtime domains axes and active layers only" {
    const dimensions: Dimensions = .{
        .biological_domain_count = 3,
        .root_axis_count = 13,
        .soil_layer_capacity = 7,
        .first_active_soil_layer = 2,
        .active_soil_layer_count = 4,
    };
    var fixture = try Fixture.init(std.testing.allocator, dimensions);
    defer fixture.deinit();
    const result = try resetIfRequired(fixture.state, dimensions, live());
    try std.testing.expectEqual(canopy_guard.Result.reset, result);
    for (fixture.state.root_index_by_domain_axis) |value|
        try std.testing.expectEqual(@as(u32, 0), value);
    for (0..dimensions.biological_domain_count) |domain| {
        for (0..dimensions.soil_layer_capacity) |layer| {
            const value = fixture.state.root_length_m_by_domain_layer[
                domain * dimensions.soil_layer_capacity + layer
            ];
            try std.testing.expectEqual(
                @as(f64, if (layer >= 2 and layer < 6) 0 else 9),
                value,
            );
        }
    }
}

test "inherited inactive and fully dead guards preserve invalid storage" {
    var empty: [0]f64 = .{};
    const state: State = .{
        .root_index_by_domain_axis = &.{},
        .root_checked_by_domain_axis = &.{},
        .total_root_length_m_by_domain = &empty,
        .shoot_root_transfer_g_c_by_domain_layer = &empty,
        .primary_root_count_by_domain_layer = &empty,
        .living_root_count_by_domain_layer = &empty,
        .maintenance_respiration_g_c_by_domain_layer = &empty,
        .nitrogen_respiration_g_c_by_domain_layer = &empty,
        .actual_respiration_g_c_by_domain_layer = &empty,
        .root_length_m_by_domain_layer = &empty,
        .primary_sink_m_by_domain_layer_axis = &empty,
        .secondary_sink_m_by_domain_layer_axis = &empty,
    };
    try std.testing.expectEqual(
        canopy_guard.Result.skipped_inactive,
        try resetIfRequired(state, std.mem.zeroes(Dimensions), .{
            .activity = .inactive,
            .shoot_mortality = .alive,
            .root_mortality = .alive,
        }),
    );
}

test "late dimension mismatch leaves every preceding owner unchanged" {
    const dimensions: Dimensions = .{
        .biological_domain_count = 2,
        .root_axis_count = 3,
        .soil_layer_capacity = 4,
        .first_active_soil_layer = 1,
        .active_soil_layer_count = 2,
    };
    var fixture = try Fixture.init(std.testing.allocator, dimensions);
    defer fixture.deinit();
    const before = fixture.state.root_index_by_domain_axis[0];
    fixture.state.secondary_sink_m_by_domain_layer_axis =
        fixture.state.secondary_sink_m_by_domain_layer_axis[0..2];
    try std.testing.expectError(
        error.PlantRootWorkspaceDimensionMismatch,
        resetIfRequired(fixture.state, dimensions, live()),
    );
    try std.testing.expectEqual(before, fixture.state.root_index_by_domain_axis[0]);
}

fn live() canopy_guard.Conditions {
    return .{
        .activity = .active,
        .shoot_mortality = .alive,
        .root_mortality = .alive,
    };
}

const Fixture = struct {
    allocator: std.mem.Allocator,
    u32_values: []u32,
    bool_values: []bool,
    f64_values: []f64,
    state: State,

    fn init(allocator: std.mem.Allocator, dimensions: Dimensions) !Fixture {
        const domain_axis = dimensions.biological_domain_count * dimensions.root_axis_count;
        const domain_layer = dimensions.biological_domain_count * dimensions.soil_layer_capacity;
        const domain_layer_axis = domain_layer * dimensions.root_axis_count;
        const f64_count = dimensions.biological_domain_count + 7 * domain_layer + 2 * domain_layer_axis;
        const u32_values = try allocator.alloc(u32, domain_axis);
        errdefer allocator.free(u32_values);
        const bool_values = try allocator.alloc(bool, domain_axis);
        errdefer allocator.free(bool_values);
        const f64_values = try allocator.alloc(f64, f64_count);
        @memset(u32_values, 8);
        @memset(bool_values, true);
        @memset(f64_values, 9);
        var cursor: usize = 0;
        const total = take(f64_values, &cursor, dimensions.biological_domain_count);
        const transfer = take(f64_values, &cursor, domain_layer);
        const primary_count = take(f64_values, &cursor, domain_layer);
        const living_count = take(f64_values, &cursor, domain_layer);
        const maintenance = take(f64_values, &cursor, domain_layer);
        const nitrogen = take(f64_values, &cursor, domain_layer);
        const actual = take(f64_values, &cursor, domain_layer);
        const length = take(f64_values, &cursor, domain_layer);
        const primary_sink = take(f64_values, &cursor, domain_layer_axis);
        const secondary_sink = take(f64_values, &cursor, domain_layer_axis);
        return .{
            .allocator = allocator,
            .u32_values = u32_values,
            .bool_values = bool_values,
            .f64_values = f64_values,
            .state = .{
                .root_index_by_domain_axis = u32_values,
                .root_checked_by_domain_axis = bool_values,
                .total_root_length_m_by_domain = total,
                .shoot_root_transfer_g_c_by_domain_layer = transfer,
                .primary_root_count_by_domain_layer = primary_count,
                .living_root_count_by_domain_layer = living_count,
                .maintenance_respiration_g_c_by_domain_layer = maintenance,
                .nitrogen_respiration_g_c_by_domain_layer = nitrogen,
                .actual_respiration_g_c_by_domain_layer = actual,
                .root_length_m_by_domain_layer = length,
                .primary_sink_m_by_domain_layer_axis = primary_sink,
                .secondary_sink_m_by_domain_layer_axis = secondary_sink,
            },
        };
    }

    fn deinit(self: *Fixture) void {
        self.allocator.free(self.f64_values);
        self.allocator.free(self.bool_values);
        self.allocator.free(self.u32_values);
        self.* = undefined;
    }
};

fn take(values: []f64, cursor: *usize, count: usize) []f64 {
    const result = values[cursor.* .. cursor.* + count];
    cursor.* += count;
    return result;
}
