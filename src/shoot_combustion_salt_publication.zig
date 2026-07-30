const std = @import("std");

pub const salt_species_count: usize = 8;

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    plant_species_per_cell: usize,
    released_mol_per_h_by_cell_and_salt: []f64,
    total_released_mol_per_h_by_cell: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        plant_species_per_cell: usize,
    ) !State {
        if (cell_count == 0 or plant_species_per_cell == 0)
            return error.InvalidShootCombustionSaltPublicationDimensions;
        const salt_values = try std.math.mul(usize, cell_count, salt_species_count);
        const values = try allocator.alloc(
            f64,
            try std.math.add(usize, salt_values, cell_count),
        );
        @memset(values, 0);
        return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .plant_species_per_cell = plant_species_per_cell,
            .released_mol_per_h_by_cell_and_salt = values[0..salt_values],
            .total_released_mol_per_h_by_cell = values[salt_values..],
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(
            self.released_mol_per_h_by_cell_and_salt.ptr[0 .. self.released_mol_per_h_by_cell_and_salt.len + self.cell_count],
        );
        self.* = undefined;
    }
};

pub const Inputs = struct {
    dynamic_salts_enabled: bool,
    fire_active_by_cell: []const bool,
    plant_branch_offsets: []const usize,
    released_mol_per_h_by_branch_and_salt: []const f64,
};

/// EXTRACT lines 298–315. Salt order is Al, Fe, Ca, Mg, Na, K, SO4, Cl.
/// Branches are traversed inside runtime plant-species order.
pub fn refresh(state: *State, inputs: Inputs) !void {
    const plant_count = try std.math.mul(
        usize,
        state.cell_count,
        state.plant_species_per_cell,
    );
    if (inputs.fire_active_by_cell.len != state.cell_count or
        inputs.plant_branch_offsets.len != plant_count + 1 or
        inputs.plant_branch_offsets[0] != 0)
        return error.InvalidShootCombustionSaltPublicationDimensions;
    const branch_count = inputs.plant_branch_offsets[plant_count];
    if (inputs.released_mol_per_h_by_branch_and_salt.len !=
        try std.math.mul(usize, branch_count, salt_species_count))
        return error.InvalidShootCombustionSaltPublicationDimensions;
    try validateOffsets(inputs.plant_branch_offsets, branch_count);

    if (!inputs.dynamic_salts_enabled) {
        @memset(state.released_mol_per_h_by_cell_and_salt, 0);
        @memset(state.total_released_mol_per_h_by_cell, 0);
        return;
    }
    for (0..state.cell_count) |cell| _ = try totalsForCell(state, inputs, cell);

    @memset(state.released_mol_per_h_by_cell_and_salt, 0);
    @memset(state.total_released_mol_per_h_by_cell, 0);
    for (0..state.cell_count) |cell| {
        const totals = totalsForCell(state, inputs, cell) catch unreachable;
        @memcpy(
            state.released_mol_per_h_by_cell_and_salt[cell * salt_species_count .. (cell + 1) * salt_species_count],
            &totals.by_salt,
        );
        state.total_released_mol_per_h_by_cell[cell] = totals.total;
    }
}

const Totals = struct {
    by_salt: [salt_species_count]f64 = @splat(0),
    total: f64 = 0,
};

fn totalsForCell(state: *const State, inputs: Inputs, cell: usize) !Totals {
    var result: Totals = .{};
    if (!inputs.fire_active_by_cell[cell]) return result;
    for (0..state.plant_species_per_cell) |species| {
        const plant = cell * state.plant_species_per_cell + species;
        for (inputs.plant_branch_offsets[plant]..inputs.plant_branch_offsets[plant + 1]) |branch| {
            for (0..salt_species_count) |salt| {
                const value = inputs.released_mol_per_h_by_branch_and_salt[
                    branch * salt_species_count + salt
                ];
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidShootCombustionSaltPublicationInput;
                result.by_salt[salt] += value;
                result.total += value;
                if (!std.math.isFinite(result.by_salt[salt]) or
                    !std.math.isFinite(result.total))
                    return error.NonFiniteShootCombustionSaltPublication;
            }
        }
    }
    return result;
}

fn validateOffsets(offsets: []const usize, branch_count: usize) !void {
    var previous: usize = 0;
    for (offsets) |offset| {
        if (offset < previous or offset > branch_count)
            return error.InvalidShootCombustionSaltPublicationTopology;
        previous = offset;
    }
}

test "shoot combustion salts preserve runtime plant branch and salt order" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    const released = [_]f64{
        1,   2,   3,   4,   5,   6,   7,   8,
        10,  20,  30,  40,  50,  60,  70,  80,
        100, 200, 300, 400, 500, 600, 700, 800,
    };
    try refresh(&state, .{
        .dynamic_salts_enabled = true,
        .fire_active_by_cell = &.{true},
        .plant_branch_offsets = &.{ 0, 2, 3 },
        .released_mol_per_h_by_branch_and_salt = &released,
    });
    try std.testing.expectEqualSlices(
        f64,
        &.{ 111, 222, 333, 444, 555, 666, 777, 888 },
        state.released_mol_per_h_by_cell_and_salt,
    );
    try std.testing.expectEqual(@as(f64, 3996), state.total_released_mol_per_h_by_cell[0]);
}

test "inactive fire ignores stale branch salt and publishes zero" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    const invalid = [_]f64{std.math.nan(f64)} ** salt_species_count;
    try refresh(&state, .{
        .dynamic_salts_enabled = true,
        .fire_active_by_cell = &.{false},
        .plant_branch_offsets = &.{ 0, 1 },
        .released_mol_per_h_by_branch_and_salt = &invalid,
    });
    for (state.released_mol_per_h_by_cell_and_salt) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqual(@as(f64, 0), state.total_released_mol_per_h_by_cell[0]);
}

test "late invalid shoot salt preserves complete publication" {
    var state = try State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    @memset(state.released_mol_per_h_by_cell_and_salt, 7);
    @memset(state.total_released_mol_per_h_by_cell, 7);
    var released = [_]f64{1} ** (2 * salt_species_count);
    released[released.len - 1] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidShootCombustionSaltPublicationInput,
        refresh(&state, .{
            .dynamic_salts_enabled = true,
            .fire_active_by_cell = &.{ true, true },
            .plant_branch_offsets = &.{ 0, 1, 2 },
            .released_mol_per_h_by_branch_and_salt = &released,
        }),
    );
    for (state.released_mol_per_h_by_cell_and_salt) |value|
        try std.testing.expectEqual(@as(f64, 7), value);
    for (state.total_released_mol_per_h_by_cell) |value|
        try std.testing.expectEqual(@as(f64, 7), value);
}
