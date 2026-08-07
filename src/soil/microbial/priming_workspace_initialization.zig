const std = @import("std");

pub const Inputs = struct {
    substrate_unlimited_respiration_g_c: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    complex_count: usize,
    population_count: usize,
    component_count: usize,
    total_respiration_g_c: []f64,
    activity_transfer_g_c: []f64,
    dissolved_carbon_transfer_g_c: []f64,
    dissolved_nitrogen_transfer_g_n: []f64,
    dissolved_phosphorus_transfer_g_p: []f64,
    acetate_transfer_g_c: []f64,
    microbial_carbon_change_g_c: []f64,
    microbial_nitrogen_change_g_n: []f64,
    microbial_phosphorus_change_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        complex_count: usize,
        population_count: usize,
        component_count: usize,
    ) !State {
        if (complex_count == 0 or population_count == 0 or component_count == 0)
            return error.InvalidPrimingWorkspaceDimensions;
        const population_items = std.math.mul(
            usize,
            complex_count,
            population_count,
        ) catch return error.InvalidPrimingWorkspaceDimensions;
        const component_items = std.math.mul(
            usize,
            population_items,
            component_count,
        ) catch return error.InvalidPrimingWorkspaceDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.complex_count = complex_count;
        state.population_count = population_count;
        state.component_count = component_count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const count = if (comptime std.mem.startsWith(
                u8,
                field.name,
                "microbial_",
            )) component_items else complex_count;
            @field(state, field.name) = try allocator.alloc(f64, count);
            @memset(@field(state, field.name), 0);
            allocated += 1;
        };
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// Exact NITRO.F 2974--3001 pre-priming aggregation and workspace reset.
pub fn initialize(state: *State, inputs: Inputs) !void {
    const population_items = std.math.mul(
        usize,
        state.complex_count,
        state.population_count,
    ) catch return error.InvalidPrimingWorkspaceDimensions;
    if (inputs.substrate_unlimited_respiration_g_c.len != population_items)
        return error.InvalidPrimingWorkspaceDimensions;
    for (inputs.substrate_unlimited_respiration_g_c) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPrimingWorkspaceInput;

    const totals = try state.allocator.alloc(f64, state.complex_count);
    defer state.allocator.free(totals);
    for (0..state.complex_count) |complex| {
        var total: f64 = 0;
        for (0..state.population_count) |population| {
            const item = complex * state.population_count + population;
            total += inputs.substrate_unlimited_respiration_g_c[item];
            if (!std.math.isFinite(total))
                return error.NonFinitePrimingWorkspaceResult;
        }
        totals[complex] = total;
    }

    @memcpy(state.total_respiration_g_c, totals);
    @memset(state.activity_transfer_g_c, 0);
    @memset(state.dissolved_carbon_transfer_g_c, 0);
    @memset(state.dissolved_nitrogen_transfer_g_n, 0);
    @memset(state.dissolved_phosphorus_transfer_g_p, 0);
    @memset(state.acetate_transfer_g_c, 0);
    @memset(state.microbial_carbon_change_g_c, 0);
    @memset(state.microbial_nitrogen_change_g_n, 0);
    @memset(state.microbial_phosphorus_change_g_p, 0);
}

test "NITRO 2974-3001 aggregates populations then resets every accumulator" {
    var state = try State.init(std.testing.allocator, 2, 3, 3);
    defer state.deinit();
    @memset(state.activity_transfer_g_c, 7);
    @memset(state.dissolved_carbon_transfer_g_c, 7);
    @memset(state.dissolved_nitrogen_transfer_g_n, 7);
    @memset(state.dissolved_phosphorus_transfer_g_p, 7);
    @memset(state.acetate_transfer_g_c, 7);
    @memset(state.microbial_carbon_change_g_c, 7);
    @memset(state.microbial_nitrogen_change_g_n, 7);
    @memset(state.microbial_phosphorus_change_g_p, 7);
    try initialize(&state, .{
        .substrate_unlimited_respiration_g_c = &.{ 1, 2, 3, 4, 5, 6 },
    });
    try std.testing.expectEqual(6, state.total_respiration_g_c[0]);
    try std.testing.expectEqual(15, state.total_respiration_g_c[1]);
    inline for (.{
        state.activity_transfer_g_c,
        state.dissolved_carbon_transfer_g_c,
        state.dissolved_nitrogen_transfer_g_n,
        state.dissolved_phosphorus_transfer_g_p,
        state.acetate_transfer_g_c,
        state.microbial_carbon_change_g_c,
        state.microbial_nitrogen_change_g_n,
        state.microbial_phosphorus_change_g_p,
    }) |values| for (values) |value| try std.testing.expectEqual(0, value);
}

test "NITRO 2974-3001 overflow preserves entire priming workspace" {
    var state = try State.init(std.testing.allocator, 1, 2, 3);
    defer state.deinit();
    state.total_respiration_g_c[0] = 7;
    @memset(state.activity_transfer_g_c, 8);
    try std.testing.expectError(
        error.NonFinitePrimingWorkspaceResult,
        initialize(&state, .{
            .substrate_unlimited_respiration_g_c = &.{ std.math.floatMax(f64), std.math.floatMax(f64) },
        }),
    );
    try std.testing.expectEqual(7, state.total_respiration_g_c[0]);
    try std.testing.expectEqual(8, state.activity_transfer_g_c[0]);
}

test "NITRO priming workspace invalid dimensions are atomic" {
    var state = try State.init(std.testing.allocator, 1, 2, 3);
    defer state.deinit();
    state.microbial_carbon_change_g_c[0] = 7;
    try std.testing.expectError(
        error.InvalidPrimingWorkspaceDimensions,
        initialize(&state, .{ .substrate_unlimited_respiration_g_c = &.{1} }),
    );
    try std.testing.expectEqual(7, state.microbial_carbon_change_g_c[0]);
}
