const std = @import("std");

pub const gas_count: usize = 6;

/// EXTRACT order: CO2-C, O2-O, CH4-C, N2O-N, NH3-N, H2-H.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    loss_g_element_per_h_by_gas_and_cell: [gas_count][]f64,

    pub fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        species_count: usize,
    ) !State {
        if (cell_count == 0 or species_count == 0)
            return error.InvalidRootGasWithdrawalPublicationDimensions;
        const values = try allocator.alloc(
            f64,
            try std.math.mul(usize, gas_count, cell_count),
        );
        @memset(values, 0);
        var slices: [gas_count][]f64 = undefined;
        for (&slices, 0..) |*slice, gas| {
            const first = gas * cell_count;
            slice.* = values[first .. first + cell_count];
        }
        return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .species_count = species_count,
            .loss_g_element_per_h_by_gas_and_cell = slices,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(
            self.loss_g_element_per_h_by_gas_and_cell[0]
                .ptr[0 .. gas_count * self.cell_count],
        );
        self.* = undefined;
    }
};

pub const Inputs = struct {
    loss_g_element_per_h_by_gas_and_plant: [gas_count][]const f64,
};

/// Exact EXTRACT lines 952–957. Source-signed withdrawal ledgers from every
/// configured plant species are preflighted and published atomically.
pub fn refresh(state: *State, inputs: Inputs) !void {
    const plant_count = try std.math.mul(
        usize,
        state.cell_count,
        state.species_count,
    );
    inline for (inputs.loss_g_element_per_h_by_gas_and_plant) |values|
        if (values.len != plant_count)
            return error.InvalidRootGasWithdrawalPublicationDimensions;

    for (0..state.cell_count) |cell| {
        for (0..gas_count) |gas|
            _ = try totalFor(state, inputs, cell, gas);
    }

    inline for (state.loss_g_element_per_h_by_gas_and_cell) |values|
        @memset(values, 0);
    for (0..state.cell_count) |cell| {
        for (0..gas_count) |gas|
            state.loss_g_element_per_h_by_gas_and_cell[gas][cell] =
                totalFor(state, inputs, cell, gas) catch unreachable;
    }
}

fn totalFor(
    state: *const State,
    inputs: Inputs,
    cell: usize,
    gas: usize,
) !f64 {
    var total: f64 = 0;
    for (0..state.species_count) |species| {
        const plant = cell * state.species_count + species;
        const loss = inputs.loss_g_element_per_h_by_gas_and_plant[gas][plant];
        if (!std.math.isFinite(loss))
            return error.NonFiniteRootGasWithdrawalPublicationInput;
        total += loss;
    }
    if (!std.math.isFinite(total))
        return error.NonFiniteRootGasWithdrawalPublication;
    return total;
}

test "root gas withdrawal publication preserves gas order and source signs" {
    var state = try State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    const gas0 = [_]f64{ -1, -2, 1, -3, -4, 2 };
    const gas1 = [_]f64{ -2, -4, 2, -6, -8, 4 };
    const gas2 = [_]f64{ -3, -6, 3, -9, -12, 6 };
    const gas3 = [_]f64{ -4, -8, 4, -12, -16, 8 };
    const gas4 = [_]f64{ -5, -10, 5, -15, -20, 10 };
    const gas5 = [_]f64{ -6, -12, 6, -18, -24, 12 };
    try refresh(&state, .{
        .loss_g_element_per_h_by_gas_and_plant = .{ &gas0, &gas1, &gas2, &gas3, &gas4, &gas5 },
    });
    const expected_first = [_]f64{ -2, -4, -6, -8, -10, -12 };
    const expected_second = [_]f64{ -5, -10, -15, -20, -25, -30 };
    for (
        state.loss_g_element_per_h_by_gas_and_cell,
        expected_first,
        expected_second,
    ) |values, first, second| {
        try std.testing.expectEqual(first, values[0]);
        try std.testing.expectEqual(second, values[1]);
    }
}

test "late invalid root gas withdrawal preserves every published array" {
    var state = try State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    for (state.loss_g_element_per_h_by_gas_and_cell, 0..) |values, gas|
        @memset(values, @as(f64, @floatFromInt(gas + 1)));
    try std.testing.expectError(
        error.NonFiniteRootGasWithdrawalPublicationInput,
        refresh(&state, .{
            .loss_g_element_per_h_by_gas_and_plant = .{
                &.{ 1, 1 },
                &.{ 1, 1 },
                &.{ 1, 1 },
                &.{ 1, 1 },
                &.{ 1, 1 },
                &.{ 1, std.math.nan(f64) },
            },
        }),
    );
    for (
        state.loss_g_element_per_h_by_gas_and_cell,
        0..,
    ) |values, gas| for (values) |value|
        try std.testing.expectEqual(
            @as(f64, @floatFromInt(gas + 1)),
            value,
        );
}
