const std = @import("std");

pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    combusted_carbon_g_c_by_plant: []f64,
    combusted_nitrogen_g_n_by_plant: []f64,
    combusted_phosphorus_g_p_by_plant: []f64,
    signed_carbon_loss_g_c_by_plant: []f64,
    signed_nitrogen_loss_g_n_by_plant: []f64,
    signed_phosphorus_loss_g_p_by_plant: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize) !State {
        if (plant_count == 0) return error.InvalidStandingDeadCharcoalDimensions;
        const values = try allocator.alloc(f64, try std.math.mul(usize, plant_count, 6));
        @memset(values, 0);
        return .{
            .allocator = allocator,
            .plant_count = plant_count,
            .combusted_carbon_g_c_by_plant = values[0 * plant_count .. 1 * plant_count],
            .combusted_nitrogen_g_n_by_plant = values[1 * plant_count .. 2 * plant_count],
            .combusted_phosphorus_g_p_by_plant = values[2 * plant_count .. 3 * plant_count],
            .signed_carbon_loss_g_c_by_plant = values[3 * plant_count .. 4 * plant_count],
            .signed_nitrogen_loss_g_n_by_plant = values[4 * plant_count .. 5 * plant_count],
            .signed_phosphorus_loss_g_p_by_plant = values[5 * plant_count .. 6 * plant_count],
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.combusted_carbon_g_c_by_plant.ptr[0 .. self.plant_count * 6]);
        self.* = undefined;
    }
};

pub const Parameters = struct {
    gas_constant_j_per_mol_k: f64,
    minimum_combustion_temperature_k: f64,
    maximum_temperature_response: f64,
    arrhenius_intercept: f64,
    activation_energy_j_per_mol: f64,
    specific_combustion_g_c_per_m2_h: f64,
};

pub const Inputs = struct {
    plant_offsets_by_cell: []const usize,
    fire_active_by_cell: []const bool,
    canopy_temperature_k_by_cell: []const f64,
    horizontal_area_m2_by_cell: []const f64,
    timestep_h: f64,
    negligible_standing_dead_g_c_by_cell: []const f64,
    total_standing_dead_carbon_g_c_by_cell: []const f64,
    parameters: Parameters,
};

pub const Pools = struct {
    charcoal_carbon_g_c_by_plant: []f64,
    charcoal_nitrogen_g_n_by_plant: []f64,
    charcoal_phosphorus_g_p_by_plant: []f64,
};

/// grosub.f lines 12140–12179. The combustion fraction is cell-wide because
/// legacy `TWTSTG` is the total standing-dead C of every PFT and fraction.
pub fn apply(state: *State, pools: Pools, inputs: Inputs) !void {
    try validate(state, pools, inputs);
    const cell_count = inputs.fire_active_by_cell.len;
    for (0..cell_count) |cell| {
        const fraction = try fractionForCell(inputs, cell);
        for (inputs.plant_offsets_by_cell[cell]..inputs.plant_offsets_by_cell[cell + 1]) |plant| {
            _ = try removal(pools, plant, fraction);
        }
    }
    @memset(state.combusted_carbon_g_c_by_plant, 0);
    @memset(state.combusted_nitrogen_g_n_by_plant, 0);
    @memset(state.combusted_phosphorus_g_p_by_plant, 0);
    @memset(state.signed_carbon_loss_g_c_by_plant, 0);
    @memset(state.signed_nitrogen_loss_g_n_by_plant, 0);
    @memset(state.signed_phosphorus_loss_g_p_by_plant, 0);
    for (0..cell_count) |cell| {
        const fraction = fractionForCell(inputs, cell) catch unreachable;
        for (inputs.plant_offsets_by_cell[cell]..inputs.plant_offsets_by_cell[cell + 1]) |plant| {
            const removed = removal(pools, plant, fraction) catch unreachable;
            pools.charcoal_carbon_g_c_by_plant[plant] -= removed[0];
            pools.charcoal_nitrogen_g_n_by_plant[plant] -= removed[1];
            pools.charcoal_phosphorus_g_p_by_plant[plant] -= removed[2];
            state.combusted_carbon_g_c_by_plant[plant] = removed[0];
            state.combusted_nitrogen_g_n_by_plant[plant] = removed[1];
            state.combusted_phosphorus_g_p_by_plant[plant] = removed[2];
            state.signed_carbon_loss_g_c_by_plant[plant] = -removed[0];
            state.signed_nitrogen_loss_g_n_by_plant[plant] = -removed[1];
            state.signed_phosphorus_loss_g_p_by_plant[plant] = -removed[2];
        }
    }
}

fn fractionForCell(inputs: Inputs, cell: usize) !f64 {
    if (!inputs.fire_active_by_cell[cell]) return 0;
    const temperature_k = inputs.canopy_temperature_k_by_cell[cell];
    const area_m2 = inputs.horizontal_area_m2_by_cell[cell];
    const total_g_c = inputs.total_standing_dead_carbon_g_c_by_cell[cell];
    const negligible_g_c = inputs.negligible_standing_dead_g_c_by_cell[cell];
    const p = inputs.parameters;
    if (temperature_k <= p.minimum_combustion_temperature_k) return 0;
    const response = @min(
        p.maximum_temperature_response,
        @exp(p.arrhenius_intercept - p.activation_energy_j_per_mol /
            (p.gas_constant_j_per_mol_k * temperature_k)),
    );
    const capacity_g_c = p.specific_combustion_g_c_per_m2_h *
        response * area_m2 * inputs.timestep_h;
    if (!std.math.isFinite(response) or !std.math.isFinite(capacity_g_c))
        return error.NonFiniteStandingDeadCharcoalCombustion;
    return if (total_g_c > negligible_g_c) @min(1, capacity_g_c / total_g_c) else 0;
}

fn removal(pools: Pools, plant: usize, fraction: f64) ![3]f64 {
    const result = [3]f64{
        pools.charcoal_carbon_g_c_by_plant[plant] * fraction,
        pools.charcoal_nitrogen_g_n_by_plant[plant] * fraction,
        pools.charcoal_phosphorus_g_p_by_plant[plant] * fraction,
    };
    for (result) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidStandingDeadCharcoalPool;
    return result;
}

fn validate(state: *const State, pools: Pools, inputs: Inputs) !void {
    const cells = inputs.fire_active_by_cell.len;
    if (cells == 0 or inputs.plant_offsets_by_cell.len != cells + 1 or
        inputs.plant_offsets_by_cell[0] != 0 or
        inputs.plant_offsets_by_cell[cells] != state.plant_count or
        inputs.canopy_temperature_k_by_cell.len != cells or
        inputs.horizontal_area_m2_by_cell.len != cells or
        inputs.negligible_standing_dead_g_c_by_cell.len != cells or
        inputs.total_standing_dead_carbon_g_c_by_cell.len != cells or
        pools.charcoal_carbon_g_c_by_plant.len != state.plant_count or
        pools.charcoal_nitrogen_g_n_by_plant.len != state.plant_count or
        pools.charcoal_phosphorus_g_p_by_plant.len != state.plant_count)
        return error.InvalidStandingDeadCharcoalDimensions;
    var previous: usize = 0;
    for (inputs.plant_offsets_by_cell) |offset| {
        if (offset < previous or offset > state.plant_count)
            return error.InvalidStandingDeadCharcoalTopology;
        previous = offset;
    }
    const p = inputs.parameters;
    inline for (.{ inputs.timestep_h, p.gas_constant_j_per_mol_k, p.minimum_combustion_temperature_k, p.maximum_temperature_response, p.arrhenius_intercept, p.activation_energy_j_per_mol, p.specific_combustion_g_c_per_m2_h }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteStandingDeadCharcoalInput;
    if (inputs.timestep_h < 0 or p.gas_constant_j_per_mol_k <= 0 or
        p.minimum_combustion_temperature_k <= 0 or p.maximum_temperature_response < 0 or
        p.activation_energy_j_per_mol < 0 or
        p.specific_combustion_g_c_per_m2_h < 0)
        return error.InvalidStandingDeadCharcoalInput;
    for (0..cells) |cell| {
        inline for (.{ inputs.canopy_temperature_k_by_cell[cell], inputs.horizontal_area_m2_by_cell[cell], inputs.negligible_standing_dead_g_c_by_cell[cell], inputs.total_standing_dead_carbon_g_c_by_cell[cell] }) |value|
            if (!std.math.isFinite(value)) return error.NonFiniteStandingDeadCharcoalInput;
        if (inputs.canopy_temperature_k_by_cell[cell] <= 0 or
            inputs.horizontal_area_m2_by_cell[cell] < 0 or
            inputs.negligible_standing_dead_g_c_by_cell[cell] < 0 or
            inputs.total_standing_dead_carbon_g_c_by_cell[cell] < 0)
            return error.InvalidStandingDeadCharcoalInput;
    }
}

fn sourceParameters() Parameters {
    return .{
        .gas_constant_j_per_mol_k = 8.3143,
        .minimum_combustion_temperature_k = 473.15,
        .maximum_temperature_response = 1,
        .arrhenius_intercept = 20.620,
        .activation_energy_j_per_mol = 120_000,
        .specific_combustion_g_c_per_m2_h = 1,
    };
}

test "cell-wide charcoal fraction preserves runtime plant order and C N P" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    var carbon = [_]f64{ 4, 6, 10 };
    var nitrogen = [_]f64{ 0.4, 0.6, 1 };
    var phosphorus = [_]f64{ 0.04, 0.06, 0.1 };
    try apply(&state, .{
        .charcoal_carbon_g_c_by_plant = &carbon,
        .charcoal_nitrogen_g_n_by_plant = &nitrogen,
        .charcoal_phosphorus_g_p_by_plant = &phosphorus,
    }, .{
        .plant_offsets_by_cell = &.{ 0, 2, 3 },
        .fire_active_by_cell = &.{ true, false },
        .canopy_temperature_k_by_cell = &.{ 700, 700 },
        .horizontal_area_m2_by_cell = &.{ 5, 5 },
        .timestep_h = 1,
        .negligible_standing_dead_g_c_by_cell = &.{ 0, 0 },
        .total_standing_dead_carbon_g_c_by_cell = &.{ 20, 10 },
        .parameters = sourceParameters(),
    });
    const fraction = 5.0 / 20.0;
    try std.testing.expectApproxEqAbs(4 * fraction, state.combusted_carbon_g_c_by_plant[0], 1e-14);
    try std.testing.expectApproxEqAbs(6 * fraction, state.combusted_carbon_g_c_by_plant[1], 1e-14);
    try std.testing.expectEqual(@as(f64, 10), carbon[2]);
    try std.testing.expectApproxEqAbs(-0.4 * fraction, state.signed_nitrogen_loss_g_n_by_plant[0], 1e-14);
    try std.testing.expectApproxEqAbs(0.04 * (1 - fraction), phosphorus[0], 1e-14);
}

test "negligible total and capacity cap reproduce source branches" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var carbon = [_]f64{ 2, 3 };
    var nitrogen = [_]f64{ 0.2, 0.3 };
    var phosphorus = [_]f64{ 0.02, 0.03 };
    var parameters = sourceParameters();
    parameters.specific_combustion_g_c_per_m2_h = 100;
    try apply(&state, .{ .charcoal_carbon_g_c_by_plant = &carbon, .charcoal_nitrogen_g_n_by_plant = &nitrogen, .charcoal_phosphorus_g_p_by_plant = &phosphorus }, .{
        .plant_offsets_by_cell = &.{ 0, 1, 2 },
        .fire_active_by_cell = &.{ true, true },
        .canopy_temperature_k_by_cell = &.{ 700, 700 },
        .horizontal_area_m2_by_cell = &.{ 1, 1 },
        .timestep_h = 1,
        .negligible_standing_dead_g_c_by_cell = &.{ 2, 0 },
        .total_standing_dead_carbon_g_c_by_cell = &.{ 2, 3 },
        .parameters = parameters,
    });
    try std.testing.expectEqual(@as(f64, 2), carbon[0]);
    try std.testing.expectEqual(@as(f64, 0), carbon[1]);
    try std.testing.expectEqual(@as(f64, -3), state.signed_carbon_loss_g_c_by_plant[1]);
}

test "strict canopy temperature threshold disables charcoal combustion" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var carbon = [_]f64{2};
    var nitrogen = [_]f64{0.2};
    var phosphorus = [_]f64{0.02};
    var parameters = sourceParameters();
    parameters.minimum_combustion_temperature_k = 700;
    try apply(&state, .{ .charcoal_carbon_g_c_by_plant = &carbon, .charcoal_nitrogen_g_n_by_plant = &nitrogen, .charcoal_phosphorus_g_p_by_plant = &phosphorus }, .{
        .plant_offsets_by_cell = &.{ 0, 1 },
        .fire_active_by_cell = &.{true},
        .canopy_temperature_k_by_cell = &.{700},
        .horizontal_area_m2_by_cell = &.{100},
        .timestep_h = 1,
        .negligible_standing_dead_g_c_by_cell = &.{0},
        .total_standing_dead_carbon_g_c_by_cell = &.{2},
        .parameters = parameters,
    });
    try std.testing.expectEqual(@as(f64, 2), carbon[0]);
    try std.testing.expectEqual(@as(f64, 0), state.combusted_carbon_g_c_by_plant[0]);
}

test "late invalid plant preserves pools and publications" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 7);
    var carbon = [_]f64{ 2, 3 };
    var nitrogen = [_]f64{ 0.2, -0.3 };
    var phosphorus = [_]f64{ 0.02, 0.03 };
    try std.testing.expectError(error.InvalidStandingDeadCharcoalPool, apply(&state, .{ .charcoal_carbon_g_c_by_plant = &carbon, .charcoal_nitrogen_g_n_by_plant = &nitrogen, .charcoal_phosphorus_g_p_by_plant = &phosphorus }, .{
        .plant_offsets_by_cell = &.{ 0, 2 },
        .fire_active_by_cell = &.{true},
        .canopy_temperature_k_by_cell = &.{700},
        .horizontal_area_m2_by_cell = &.{1},
        .timestep_h = 1,
        .negligible_standing_dead_g_c_by_cell = &.{0},
        .total_standing_dead_carbon_g_c_by_cell = &.{5},
        .parameters = sourceParameters(),
    }));
    try std.testing.expectEqualSlices(f64, &.{ 2, 3 }, &carbon);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            for (@field(state, field.name)) |value| try std.testing.expectEqual(@as(f64, 7), value);
}
