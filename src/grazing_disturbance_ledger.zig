const std = @import("std");

pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    cell_count: usize,
    exported_carbon_g_c_by_plant: []f64,
    exported_nitrogen_g_n_by_plant: []f64,
    exported_phosphorus_g_p_by_plant: []f64,
    cumulative_respiration_carbon_g_c_by_plant: []f64,
    cumulative_aboveground_respiration_carbon_g_c_by_plant: []f64,
    exported_carbon_g_c_by_cell: []f64,
    exported_nitrogen_g_n_by_cell: []f64,
    exported_phosphorus_g_p_by_cell: []f64,
    ecosystem_respiration_carbon_g_c_by_cell: []f64,
    autotrophic_respiration_carbon_g_c_by_cell: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize, cell_count: usize) !State {
        if (plant_count == 0 or cell_count == 0) return error.InvalidGrazingDisturbanceLedgerDimensions;
        const plant_values = try allocator.alloc(f64, try std.math.mul(usize, plant_count, 5));
        errdefer allocator.free(plant_values);
        const cell_values = try allocator.alloc(f64, try std.math.mul(usize, cell_count, 5));
        @memset(plant_values, 0);
        @memset(cell_values, 0);
        return .{
            .allocator = allocator,
            .plant_count = plant_count,
            .cell_count = cell_count,
            .exported_carbon_g_c_by_plant = plant_values[0 * plant_count .. 1 * plant_count],
            .exported_nitrogen_g_n_by_plant = plant_values[1 * plant_count .. 2 * plant_count],
            .exported_phosphorus_g_p_by_plant = plant_values[2 * plant_count .. 3 * plant_count],
            .cumulative_respiration_carbon_g_c_by_plant = plant_values[3 * plant_count .. 4 * plant_count],
            .cumulative_aboveground_respiration_carbon_g_c_by_plant = plant_values[4 * plant_count .. 5 * plant_count],
            .exported_carbon_g_c_by_cell = cell_values[0 * cell_count .. 1 * cell_count],
            .exported_nitrogen_g_n_by_cell = cell_values[1 * cell_count .. 2 * cell_count],
            .exported_phosphorus_g_p_by_cell = cell_values[2 * cell_count .. 3 * cell_count],
            .ecosystem_respiration_carbon_g_c_by_cell = cell_values[3 * cell_count .. 4 * cell_count],
            .autotrophic_respiration_carbon_g_c_by_cell = cell_values[4 * cell_count .. 5 * cell_count],
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.exported_carbon_g_c_by_plant.ptr[0 .. self.plant_count * 5]);
        self.allocator.free(self.exported_carbon_g_c_by_cell.ptr[0 .. self.cell_count * 5]);
        self.* = undefined;
    }
};

pub const Inputs = struct {
    plant_offsets_by_cell: []const usize,
    grazing_active_by_plant: []const bool,
    removed_carbon_g_c_by_plant: []const f64,
    removed_nitrogen_g_n_by_plant: []const f64,
    removed_phosphorus_g_p_by_plant: []const f64,
    returned_carbon_g_c_by_plant: []const f64,
    returned_nitrogen_g_n_by_plant: []const f64,
    returned_phosphorus_g_p_by_plant: []const f64,
    grazer_growth_yield_by_plant: []const f64,
    grazer_respired_fraction_by_plant: []const f64,
};

/// GROSUB lines 10724–10749. Carbon export uses `GY`; N/P export does not.
/// Respired carbon retains the source negative sign in all four ledgers.
pub fn accumulate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    for (0..state.cell_count) |cell| {
        var cell_next = try currentCell(state, cell);
        for (inputs.plant_offsets_by_cell[cell]..inputs.plant_offsets_by_cell[cell + 1]) |plant| {
            if (!inputs.grazing_active_by_plant[plant]) continue;
            const transaction = try calculatePlant(state, inputs, plant);
            cell_next.exported_c += transaction.exported_c;
            cell_next.exported_n += transaction.exported_n;
            cell_next.exported_p += transaction.exported_p;
            cell_next.ecosystem_respiration -= transaction.respired_c;
            cell_next.autotrophic_respiration -= transaction.respired_c;
            try validateCell(cell_next);
        }
    }
    for (0..state.cell_count) |cell| {
        var cell_next = currentCell(state, cell) catch unreachable;
        for (inputs.plant_offsets_by_cell[cell]..inputs.plant_offsets_by_cell[cell + 1]) |plant| {
            if (!inputs.grazing_active_by_plant[plant]) continue;
            const transaction = calculatePlant(state, inputs, plant) catch unreachable;
            state.exported_carbon_g_c_by_plant[plant] = transaction.next_exported_c;
            state.exported_nitrogen_g_n_by_plant[plant] = transaction.next_exported_n;
            state.exported_phosphorus_g_p_by_plant[plant] = transaction.next_exported_p;
            state.cumulative_respiration_carbon_g_c_by_plant[plant] = transaction.next_respiration;
            state.cumulative_aboveground_respiration_carbon_g_c_by_plant[plant] = transaction.next_aboveground_respiration;
            cell_next.exported_c += transaction.exported_c;
            cell_next.exported_n += transaction.exported_n;
            cell_next.exported_p += transaction.exported_p;
            cell_next.ecosystem_respiration -= transaction.respired_c;
            cell_next.autotrophic_respiration -= transaction.respired_c;
        }
        state.exported_carbon_g_c_by_cell[cell] = cell_next.exported_c;
        state.exported_nitrogen_g_n_by_cell[cell] = cell_next.exported_n;
        state.exported_phosphorus_g_p_by_cell[cell] = cell_next.exported_p;
        state.ecosystem_respiration_carbon_g_c_by_cell[cell] = cell_next.ecosystem_respiration;
        state.autotrophic_respiration_carbon_g_c_by_cell[cell] = cell_next.autotrophic_respiration;
    }
}

const PlantTransaction = struct {
    exported_c: f64,
    exported_n: f64,
    exported_p: f64,
    respired_c: f64,
    next_exported_c: f64,
    next_exported_n: f64,
    next_exported_p: f64,
    next_respiration: f64,
    next_aboveground_respiration: f64,
};

fn calculatePlant(state: *const State, inputs: Inputs, plant: usize) !PlantTransaction {
    const net_c = inputs.removed_carbon_g_c_by_plant[plant] - inputs.returned_carbon_g_c_by_plant[plant];
    const net_n = inputs.removed_nitrogen_g_n_by_plant[plant] - inputs.returned_nitrogen_g_n_by_plant[plant];
    const net_p = inputs.removed_phosphorus_g_p_by_plant[plant] - inputs.returned_phosphorus_g_p_by_plant[plant];
    const exported_c = inputs.grazer_growth_yield_by_plant[plant] * net_c;
    const respired_c = inputs.grazer_respired_fraction_by_plant[plant] * net_c;
    const transaction: PlantTransaction = .{
        .exported_c = exported_c,
        .exported_n = net_n,
        .exported_p = net_p,
        .respired_c = respired_c,
        .next_exported_c = state.exported_carbon_g_c_by_plant[plant] + exported_c,
        .next_exported_n = state.exported_nitrogen_g_n_by_plant[plant] + net_n,
        .next_exported_p = state.exported_phosphorus_g_p_by_plant[plant] + net_p,
        .next_respiration = state.cumulative_respiration_carbon_g_c_by_plant[plant] - respired_c,
        .next_aboveground_respiration = state.cumulative_aboveground_respiration_carbon_g_c_by_plant[plant] - respired_c,
    };
    inline for (@typeInfo(PlantTransaction).@"struct".fields) |field|
        if (!std.math.isFinite(@field(transaction, field.name)))
            return error.NonFiniteGrazingDisturbanceLedger;
    return transaction;
}

const CellNext = struct {
    exported_c: f64,
    exported_n: f64,
    exported_p: f64,
    ecosystem_respiration: f64,
    autotrophic_respiration: f64,
};

fn currentCell(state: *const State, cell: usize) !CellNext {
    const result: CellNext = .{
        .exported_c = state.exported_carbon_g_c_by_cell[cell],
        .exported_n = state.exported_nitrogen_g_n_by_cell[cell],
        .exported_p = state.exported_phosphorus_g_p_by_cell[cell],
        .ecosystem_respiration = state.ecosystem_respiration_carbon_g_c_by_cell[cell],
        .autotrophic_respiration = state.autotrophic_respiration_carbon_g_c_by_cell[cell],
    };
    try validateCell(result);
    return result;
}

fn validateCell(value: CellNext) !void {
    inline for (@typeInfo(CellNext).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteGrazingDisturbanceLedger;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.plant_offsets_by_cell.len != state.cell_count + 1 or inputs.plant_offsets_by_cell[0] != 0 or inputs.plant_offsets_by_cell[state.cell_count] != state.plant_count) return error.InvalidGrazingDisturbanceLedgerDimensions;
    var previous: usize = 0;
    for (inputs.plant_offsets_by_cell) |offset| {
        if (offset < previous or offset > state.plant_count) return error.InvalidGrazingDisturbanceLedgerTopology;
        previous = offset;
    }
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != state.plant_count) return error.InvalidGrazingDisturbanceLedgerDimensions;
        for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteGrazingDisturbanceLedgerInput;
    };
    if (inputs.grazing_active_by_plant.len != state.plant_count) return error.InvalidGrazingDisturbanceLedgerDimensions;
    for (0..state.plant_count) |plant| {
        if (inputs.removed_carbon_g_c_by_plant[plant] < 0 or inputs.removed_nitrogen_g_n_by_plant[plant] < 0 or inputs.removed_phosphorus_g_p_by_plant[plant] < 0 or inputs.returned_carbon_g_c_by_plant[plant] < 0 or inputs.returned_nitrogen_g_n_by_plant[plant] < 0 or inputs.returned_phosphorus_g_p_by_plant[plant] < 0 or inputs.removed_carbon_g_c_by_plant[plant] < inputs.returned_carbon_g_c_by_plant[plant] or inputs.removed_nitrogen_g_n_by_plant[plant] < inputs.returned_nitrogen_g_n_by_plant[plant] or inputs.removed_phosphorus_g_p_by_plant[plant] < inputs.returned_phosphorus_g_p_by_plant[plant] or inputs.grazer_growth_yield_by_plant[plant] < 0 or inputs.grazer_growth_yield_by_plant[plant] > 1 or inputs.grazer_respired_fraction_by_plant[plant] < 0 or inputs.grazer_respired_fraction_by_plant[plant] > 1) return error.InvalidGrazingDisturbanceLedgerInput;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const expected = if (std.mem.endsWith(u8, field.name, "_by_plant")) state.plant_count else state.cell_count;
        if (@field(state, field.name).len != expected) return error.InvalidGrazingDisturbanceLedgerDimensions;
        for (@field(state, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFiniteGrazingDisturbanceLedgerState;
    };
}

test "grazing carbon applies growth yield while N P export full net removal" {
    var state = try State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    try accumulate(&state, .{
        .plant_offsets_by_cell = &.{ 0, 2 },
        .grazing_active_by_plant = &.{ true, true },
        .removed_carbon_g_c_by_plant = &.{ 10, 20 },
        .removed_nitrogen_g_n_by_plant = &.{ 2, 4 },
        .removed_phosphorus_g_p_by_plant = &.{ 1, 2 },
        .returned_carbon_g_c_by_plant = &.{ 2, 5 },
        .returned_nitrogen_g_n_by_plant = &.{ 0.5, 1 },
        .returned_phosphorus_g_p_by_plant = &.{ 0.25, 0.5 },
        .grazer_growth_yield_by_plant = &.{ 0.25, 0.4 },
        .grazer_respired_fraction_by_plant = &.{ 0.75, 0.6 },
    });
    try std.testing.expectEqualSlices(f64, &.{ 2, 6 }, state.exported_carbon_g_c_by_plant);
    try std.testing.expectEqual(@as(f64, 8), state.exported_carbon_g_c_by_cell[0]);
    try std.testing.expectEqual(@as(f64, 4.5), state.exported_nitrogen_g_n_by_cell[0]);
    try std.testing.expectEqual(@as(f64, -15), state.ecosystem_respiration_carbon_g_c_by_cell[0]);
}

test "inactive plants preserve every plant and cell ledger" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) @memset(@field(state, field.name), 7);
    try accumulate(&state, .{
        .plant_offsets_by_cell = &.{ 0, 1 },
        .grazing_active_by_plant = &.{false},
        .removed_carbon_g_c_by_plant = &.{10},
        .removed_nitrogen_g_n_by_plant = &.{2},
        .removed_phosphorus_g_p_by_plant = &.{1},
        .returned_carbon_g_c_by_plant = &.{2},
        .returned_nitrogen_g_n_by_plant = &.{0.5},
        .returned_phosphorus_g_p_by_plant = &.{0.25},
        .grazer_growth_yield_by_plant = &.{0.25},
        .grazer_respired_fraction_by_plant = &.{0.75},
    });
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(state, field.name)) |value| try std.testing.expectEqual(@as(f64, 7), value);
}

test "late cell overflow preserves complete plant and cell ledgers" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) @memset(@field(state, field.name), 7);
    state.exported_carbon_g_c_by_cell[1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteGrazingDisturbanceLedger, accumulate(&state, .{
        .plant_offsets_by_cell = &.{ 0, 1, 2 },
        .grazing_active_by_plant = &.{ true, true },
        .removed_carbon_g_c_by_plant = &.{ 2, std.math.floatMax(f64) },
        .removed_nitrogen_g_n_by_plant = &.{ 2, 2 },
        .removed_phosphorus_g_p_by_plant = &.{ 2, 2 },
        .returned_carbon_g_c_by_plant = &.{ 0, 0 },
        .returned_nitrogen_g_n_by_plant = &.{ 0, 0 },
        .returned_phosphorus_g_p_by_plant = &.{ 0, 0 },
        .grazer_growth_yield_by_plant = &.{ 1, 1 },
        .grazer_respired_fraction_by_plant = &.{ 0, 0 },
    }));
    try std.testing.expectEqual(@as(f64, 7), state.exported_carbon_g_c_by_plant[0]);
    try std.testing.expectEqual(std.math.floatMax(f64), state.exported_carbon_g_c_by_cell[1]);
}
