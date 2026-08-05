const std = @import("std");

/// UPTAKE.F 88--317 active behavior. The source's lateral-transfer equations
/// are commented out; every active statement clears the corresponding
/// directional and cell-total flux owner.
pub const State = struct {
    allocator: std.mem.Allocator,
    row_count: usize,
    column_count: usize,
    directional_longwave_radiation_megajoules_per_step: []f64,
    directional_sensible_heat_megajoules_per_step: []f64,
    directional_water_vapor_m3_per_step: []f64,
    directional_carbon_dioxide_g_c_per_step: []f64,
    directional_methane_g_c_per_step: []f64,
    directional_oxygen_g_o_per_step: []f64,
    total_longwave_radiation_megajoules_per_step: []f64,
    total_sensible_heat_megajoules_per_step: []f64,
    total_water_vapor_m3_per_step: []f64,
    total_carbon_dioxide_g_c_per_step: []f64,
    total_methane_g_c_per_step: []f64,
    total_oxygen_g_o_per_step: []f64,

    pub fn init(allocator: std.mem.Allocator, row_count: usize, column_count: usize) !State {
        if (row_count == 0 or column_count == 0)
            return error.InvalidCanopyLateralExchangeDimensions;
        const cell_count = std.math.mul(usize, row_count, column_count) catch
            return error.InvalidCanopyLateralExchangeDimensions;
        const directional_count = std.math.mul(usize, 2, cell_count) catch
            return error.InvalidCanopyLateralExchangeDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.row_count = row_count;
        state.column_count = column_count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                const count = if (std.mem.startsWith(u8, field.name, "directional_"))
                    directional_count
                else
                    cell_count;
                @field(state, field.name) = try allocator.alloc(f64, count);
                @memset(@field(state, field.name), 0);
                allocated += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn clear(self: *State) !void {
        try self.validateDimensions();
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) @memset(@field(self, field.name), 0);
    }

    fn validateDimensions(self: *const State) !void {
        if (self.row_count == 0 or self.column_count == 0)
            return error.InvalidCanopyLateralExchangeDimensions;
        const cell_count = std.math.mul(usize, self.row_count, self.column_count) catch
            return error.InvalidCanopyLateralExchangeDimensions;
        const directional_count = std.math.mul(usize, 2, cell_count) catch
            return error.InvalidCanopyLateralExchangeDimensions;
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                const expected = if (std.mem.startsWith(u8, field.name, "directional_"))
                    directional_count
                else
                    cell_count;
                if (@field(self, field.name).len != expected)
                    return error.InvalidCanopyLateralExchangeDimensions;
            }
        }
    }
};

test "UPTAKE active lateral exchange block clears every runtime cell and face" {
    var state = try State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 7);

    try state.clear();

    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state, field.name)) |value|
            try std.testing.expectEqual(@as(f64, 0), value);
}

test "UPTAKE lateral exchange dimensions are runtime validated" {
    try std.testing.expectError(
        error.InvalidCanopyLateralExchangeDimensions,
        State.init(std.testing.allocator, 0, 3),
    );
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    const full_oxygen = state.total_oxygen_g_o_per_step;
    state.total_oxygen_g_o_per_step = state.total_oxygen_g_o_per_step[0..3];
    try std.testing.expectError(
        error.InvalidCanopyLateralExchangeDimensions,
        state.clear(),
    );
    state.total_oxygen_g_o_per_step = full_oxygen;
}

test "UPTAKE lateral exchange uses two directions for arbitrary grid shape" {
    var state = try State.init(std.testing.allocator, 3, 5);
    defer state.deinit();
    try std.testing.expectEqual(@as(usize, 15), state.total_methane_g_c_per_step.len);
    try std.testing.expectEqual(
        @as(usize, 30),
        state.directional_methane_g_c_per_step.len,
    );
}
