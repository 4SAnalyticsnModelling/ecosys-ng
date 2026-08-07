const std = @import("std");

pub const removal_category_count: usize = 5;

pub const Termination = enum {
    retain_or_terminate,
    terminate_and_reseed,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    cell_count: usize,
    exported_carbon_g_c_by_plant: []f64,
    exported_nitrogen_g_n_by_plant: []f64,
    exported_phosphorus_g_p_by_plant: []f64,
    seasonal_storage_carbon_g_c_by_plant: []f64,
    seasonal_storage_nitrogen_g_n_by_plant: []f64,
    seasonal_storage_phosphorus_g_p_by_plant: []f64,
    exported_carbon_g_c_by_cell: []f64,
    exported_nitrogen_g_n_by_cell: []f64,
    exported_phosphorus_g_p_by_cell: []f64,
    net_biome_productivity_g_c_by_cell: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize, cell_count: usize) !State {
        if (plant_count == 0 or cell_count == 0) return error.InvalidNongrazingDisturbanceLedgerDimensions;
        const plant_values = try allocator.alloc(f64, try std.math.mul(usize, plant_count, 6));
        errdefer allocator.free(plant_values);
        const cell_values = try allocator.alloc(f64, try std.math.mul(usize, cell_count, 4));
        @memset(plant_values, 0);
        @memset(cell_values, 0);
        return .{
            .allocator = allocator,
            .plant_count = plant_count,
            .cell_count = cell_count,
            .exported_carbon_g_c_by_plant = plant_values[0 * plant_count .. 1 * plant_count],
            .exported_nitrogen_g_n_by_plant = plant_values[1 * plant_count .. 2 * plant_count],
            .exported_phosphorus_g_p_by_plant = plant_values[2 * plant_count .. 3 * plant_count],
            .seasonal_storage_carbon_g_c_by_plant = plant_values[3 * plant_count .. 4 * plant_count],
            .seasonal_storage_nitrogen_g_n_by_plant = plant_values[4 * plant_count .. 5 * plant_count],
            .seasonal_storage_phosphorus_g_p_by_plant = plant_values[5 * plant_count .. 6 * plant_count],
            .exported_carbon_g_c_by_cell = cell_values[0 * cell_count .. 1 * cell_count],
            .exported_nitrogen_g_n_by_cell = cell_values[1 * cell_count .. 2 * cell_count],
            .exported_phosphorus_g_p_by_cell = cell_values[2 * cell_count .. 3 * cell_count],
            .net_biome_productivity_g_c_by_cell = cell_values[3 * cell_count .. 4 * cell_count],
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.exported_carbon_g_c_by_plant.ptr[0 .. self.plant_count * 6]);
        self.allocator.free(self.exported_carbon_g_c_by_cell.ptr[0 .. self.cell_count * 4]);
        self.* = undefined;
    }
};

pub const Inputs = struct {
    plant_offsets_by_cell: []const usize,
    active_by_plant: []const bool,
    termination_by_plant: []const Termination,
    removed_carbon_g_c_by_plant_and_category: []const f64,
    removed_nitrogen_g_n_by_plant_and_category: []const f64,
    removed_phosphorus_g_p_by_plant_and_category: []const f64,
    returned_carbon_g_c_by_plant_and_category: []const f64,
    returned_nitrogen_g_n_by_plant_and_category: []const f64,
    returned_phosphorus_g_p_by_plant_and_category: []const f64,
};

/// GROSUB 10679–10722 non-grazing branch. Category order is nonstructural,
/// leaf, fine/nonleaf, woody, then standing dead.
pub fn accumulate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    for (0..state.cell_count) |cell| {
        var cell_next = try currentCell(state, cell);
        for (inputs.plant_offsets_by_cell[cell]..inputs.plant_offsets_by_cell[cell + 1]) |plant| {
            if (!inputs.active_by_plant[plant]) continue;
            const transaction = try calculatePlant(state, inputs, plant);
            if (inputs.termination_by_plant[plant] == .retain_or_terminate) {
                cell_next.exported_c += transaction.net_c;
                cell_next.exported_n += transaction.net_n;
                cell_next.exported_p += transaction.net_p;
                cell_next.net_biome_productivity -= transaction.net_c;
                try validateCell(cell_next);
            }
        }
    }
    for (0..state.cell_count) |cell| {
        var cell_next = currentCell(state, cell) catch unreachable;
        for (inputs.plant_offsets_by_cell[cell]..inputs.plant_offsets_by_cell[cell + 1]) |plant| {
            if (!inputs.active_by_plant[plant]) continue;
            const transaction = calculatePlant(state, inputs, plant) catch unreachable;
            state.exported_carbon_g_c_by_plant[plant] = transaction.exported_carbon_g_c_by_plant;
            state.exported_nitrogen_g_n_by_plant[plant] = transaction.exported_nitrogen_g_n_by_plant;
            state.exported_phosphorus_g_p_by_plant[plant] = transaction.exported_phosphorus_g_p_by_plant;
            state.seasonal_storage_carbon_g_c_by_plant[plant] = transaction.seasonal_storage_carbon_g_c_by_plant;
            state.seasonal_storage_nitrogen_g_n_by_plant[plant] = transaction.seasonal_storage_nitrogen_g_n_by_plant;
            state.seasonal_storage_phosphorus_g_p_by_plant[plant] = transaction.seasonal_storage_phosphorus_g_p_by_plant;
            if (inputs.termination_by_plant[plant] == .retain_or_terminate) {
                cell_next.exported_c += transaction.net_c;
                cell_next.exported_n += transaction.net_n;
                cell_next.exported_p += transaction.net_p;
                cell_next.net_biome_productivity -= transaction.net_c;
            }
        }
        state.exported_carbon_g_c_by_cell[cell] = cell_next.exported_c;
        state.exported_nitrogen_g_n_by_cell[cell] = cell_next.exported_n;
        state.exported_phosphorus_g_p_by_cell[cell] = cell_next.exported_p;
        state.net_biome_productivity_g_c_by_cell[cell] = cell_next.net_biome_productivity;
    }
}

const PlantNext = struct {
    net_c: f64,
    net_n: f64,
    net_p: f64,
    exported_carbon_g_c_by_plant: f64,
    exported_nitrogen_g_n_by_plant: f64,
    exported_phosphorus_g_p_by_plant: f64,
    seasonal_storage_carbon_g_c_by_plant: f64,
    seasonal_storage_nitrogen_g_n_by_plant: f64,
    seasonal_storage_phosphorus_g_p_by_plant: f64,
};

fn calculatePlant(state: *const State, inputs: Inputs, plant: usize) !PlantNext {
    const first = plant * removal_category_count;
    var removed = [3]f64{ 0, 0, 0 };
    var returned = [3]f64{ 0, 0, 0 };
    for (0..removal_category_count) |category| {
        const index = first + category;
        removed[0] += inputs.removed_carbon_g_c_by_plant_and_category[index];
        returned[0] += inputs.returned_carbon_g_c_by_plant_and_category[index];
        removed[1] += inputs.removed_nitrogen_g_n_by_plant_and_category[index];
        returned[1] += inputs.returned_nitrogen_g_n_by_plant_and_category[index];
        removed[2] += inputs.removed_phosphorus_g_p_by_plant_and_category[index];
        returned[2] += inputs.returned_phosphorus_g_p_by_plant_and_category[index];
    }
    const net = [3]f64{ removed[0] - returned[0], removed[1] - returned[1], removed[2] - returned[2] };
    for (net) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidNongrazingDisturbanceLedgerInput;
    const reseed = inputs.termination_by_plant[plant] == .terminate_and_reseed;
    const next: PlantNext = .{
        .net_c = net[0],
        .net_n = net[1],
        .net_p = net[2],
        .exported_carbon_g_c_by_plant = state.exported_carbon_g_c_by_plant[plant] + if (reseed) 0 else net[0],
        .exported_nitrogen_g_n_by_plant = state.exported_nitrogen_g_n_by_plant[plant] + if (reseed) 0 else net[1],
        .exported_phosphorus_g_p_by_plant = state.exported_phosphorus_g_p_by_plant[plant] + if (reseed) 0 else net[2],
        .seasonal_storage_carbon_g_c_by_plant = state.seasonal_storage_carbon_g_c_by_plant[plant] + if (reseed) net[0] else 0,
        .seasonal_storage_nitrogen_g_n_by_plant = state.seasonal_storage_nitrogen_g_n_by_plant[plant] + if (reseed) net[1] else 0,
        .seasonal_storage_phosphorus_g_p_by_plant = state.seasonal_storage_phosphorus_g_p_by_plant[plant] + if (reseed) net[2] else 0,
    };
    inline for (@typeInfo(PlantNext).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name))) return error.NonFiniteNongrazingDisturbanceLedger;
    return next;
}

const CellNext = struct { exported_c: f64, exported_n: f64, exported_p: f64, net_biome_productivity: f64 };

fn currentCell(state: *const State, cell: usize) !CellNext {
    const value: CellNext = .{ .exported_c = state.exported_carbon_g_c_by_cell[cell], .exported_n = state.exported_nitrogen_g_n_by_cell[cell], .exported_p = state.exported_phosphorus_g_p_by_cell[cell], .net_biome_productivity = state.net_biome_productivity_g_c_by_cell[cell] };
    try validateCell(value);
    return value;
}

fn validateCell(value: CellNext) !void {
    inline for (@typeInfo(CellNext).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return error.NonFiniteNongrazingDisturbanceLedger;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.plant_offsets_by_cell.len != state.cell_count + 1 or inputs.plant_offsets_by_cell[0] != 0 or inputs.plant_offsets_by_cell[state.cell_count] != state.plant_count or inputs.active_by_plant.len != state.plant_count or inputs.termination_by_plant.len != state.plant_count) return error.InvalidNongrazingDisturbanceLedgerDimensions;
    var previous: usize = 0;
    for (inputs.plant_offsets_by_cell) |offset| {
        if (offset < previous or offset > state.plant_count) return error.InvalidNongrazingDisturbanceLedgerTopology;
        previous = offset;
    }
    const pool_count = state.plant_count * removal_category_count;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != pool_count) return error.InvalidNongrazingDisturbanceLedgerDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidNongrazingDisturbanceLedgerInput;
    };
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const expected = if (std.mem.endsWith(u8, field.name, "_by_plant")) state.plant_count else state.cell_count;
        if (@field(state, field.name).len != expected) return error.InvalidNongrazingDisturbanceLedgerDimensions;
        for (@field(state, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFiniteNongrazingDisturbanceLedger;
    };
}

test "ordinary disturbance exports five-category net mass and subtracts carbon from TNBP" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    const removed = [_]f64{ 1, 2, 3, 4, 5 };
    const returned = [_]f64{ 0, 1, 1, 2, 1 };
    try accumulate(&state, .{ .plant_offsets_by_cell = &.{ 0, 1 }, .active_by_plant = &.{true}, .termination_by_plant = &.{.retain_or_terminate}, .removed_carbon_g_c_by_plant_and_category = &removed, .removed_nitrogen_g_n_by_plant_and_category = &removed, .removed_phosphorus_g_p_by_plant_and_category = &removed, .returned_carbon_g_c_by_plant_and_category = &returned, .returned_nitrogen_g_n_by_plant_and_category = &returned, .returned_phosphorus_g_p_by_plant_and_category = &returned });
    try std.testing.expectEqual(@as(f64, 10), state.exported_carbon_g_c_by_plant[0]);
    try std.testing.expectEqual(@as(f64, 10), state.exported_carbon_g_c_by_cell[0]);
    try std.testing.expectEqual(@as(f64, -10), state.net_biome_productivity_g_c_by_cell[0]);
}

test "terminate and reseed diverts net removal to seasonal storage only" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    const removed = [_]f64{ 1, 2, 3, 4, 5 };
    const returned = [_]f64{0} ** 5;
    try accumulate(&state, .{ .plant_offsets_by_cell = &.{ 0, 1 }, .active_by_plant = &.{true}, .termination_by_plant = &.{.terminate_and_reseed}, .removed_carbon_g_c_by_plant_and_category = &removed, .removed_nitrogen_g_n_by_plant_and_category = &removed, .removed_phosphorus_g_p_by_plant_and_category = &removed, .returned_carbon_g_c_by_plant_and_category = &returned, .returned_nitrogen_g_n_by_plant_and_category = &returned, .returned_phosphorus_g_p_by_plant_and_category = &returned });
    try std.testing.expectEqual(@as(f64, 15), state.seasonal_storage_carbon_g_c_by_plant[0]);
    try std.testing.expectEqual(@as(f64, 0), state.exported_carbon_g_c_by_cell[0]);
    try std.testing.expectEqual(@as(f64, 0), state.net_biome_productivity_g_c_by_cell[0]);
}

test "late cell overflow preserves plant cell and storage ledgers" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) @memset(@field(state, field.name), 7);
    state.exported_carbon_g_c_by_cell[1] = std.math.floatMax(f64);
    const removed = [_]f64{ 1, 0, 0, 0, 0, std.math.floatMax(f64), 0, 0, 0, 0 };
    const returned = [_]f64{0} ** 10;
    try std.testing.expectError(error.NonFiniteNongrazingDisturbanceLedger, accumulate(&state, .{ .plant_offsets_by_cell = &.{ 0, 1, 2 }, .active_by_plant = &.{ true, true }, .termination_by_plant = &.{ .retain_or_terminate, .retain_or_terminate }, .removed_carbon_g_c_by_plant_and_category = &removed, .removed_nitrogen_g_n_by_plant_and_category = &removed, .removed_phosphorus_g_p_by_plant_and_category = &removed, .returned_carbon_g_c_by_plant_and_category = &returned, .returned_nitrogen_g_n_by_plant_and_category = &returned, .returned_phosphorus_g_p_by_plant_and_category = &returned }));
    try std.testing.expectEqual(@as(f64, 7), state.exported_carbon_g_c_by_plant[0]);
    try std.testing.expectEqual(std.math.floatMax(f64), state.exported_carbon_g_c_by_cell[1]);
}
