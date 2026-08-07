const std = @import("std");
const growth_stages = @import("../../plant/lifecycle/growth_stages.zig");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchState = struct {
    host_nonstructural: Pool,
    symbiont_nonstructural: Pool,
    host_leaf_and_petiole_carbon_g_c: f64,
    symbiont_structural_carbon_g_c: f64,
    physiological_maturity_set: bool,
};

pub const Inputs = struct {
    fixation_type: u8,
    growth_habit: growth_stages.GrowthHabit,
    host_carbon_presence_threshold_g_c: f64,
    host_tissue_presence_threshold_g_c: f64,
    combined_carbon_presence_threshold_g_c: f64,
    initial_bacterial_carbon_g_c_per_m2: f64,
    horizontal_cell_area_m2: f64,
    exchange_rate_per_h: f64,
    timestep_h: f64,
};

/// Exact grosub.f lines 5662--5703 canopy host-symbiont mobile C:N:P
/// equilibration in ascending NB order. Signed transfers are positive from
/// host to symbiont. N and P gradients intentionally use the post-carbon
/// pools assigned by source lines 5674--5675.
pub fn applyAll(states: []BranchState, host_to_symbiont_by_branch: []Pool, inputs: Inputs) !void {
    if (states.len == 0 or states.len != host_to_symbiont_by_branch.len)
        return error.CanopyHostSymbiontExchangeDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    try validateInputs(inputs);
    for (states) |state| _ = try nextState(state, inputs);
    for (states, host_to_symbiont_by_branch) |*state, *transfer| {
        const result = try nextState(state.*, inputs);
        state.* = result.state;
        transfer.* = result.host_to_symbiont;
    }
}

const ExchangeResult = struct {
    state: BranchState,
    host_to_symbiont: Pool,
};

fn nextState(state: BranchState, inputs: Inputs) !ExchangeResult {
    try validatePool(state.host_nonstructural);
    try validatePool(state.symbiont_nonstructural);
    inline for (.{ state.host_leaf_and_petiole_carbon_g_c, state.symbiont_structural_carbon_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyHostSymbiontExchangeState;

    const zero = Pool{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
    if (!(state.host_nonstructural.carbon_g_c > inputs.host_carbon_presence_threshold_g_c and
        state.host_leaf_and_petiole_carbon_g_c > inputs.host_tissue_presence_threshold_g_c and
        (inputs.growth_habit == .perennial or !state.physiological_maturity_set)))
        return .{ .state = state, .host_to_symbiont = zero };

    const host_tissue_carbon_g_c = state.host_leaf_and_petiole_carbon_g_c;
    const effective_symbiont_carbon_g_c = @min(
        state.host_leaf_and_petiole_carbon_g_c,
        @max(
            inputs.initial_bacterial_carbon_g_c_per_m2 * inputs.horizontal_cell_area_m2,
            state.symbiont_structural_carbon_g_c,
        ),
    );
    const tissue_total_g_c = host_tissue_carbon_g_c + effective_symbiont_carbon_g_c;
    if (!(tissue_total_g_c > inputs.host_carbon_presence_threshold_g_c))
        return .{ .state = state, .host_to_symbiont = zero };

    var next = state;
    var transfer = zero;
    const carbon_difference_g_c = (state.host_nonstructural.carbon_g_c * effective_symbiont_carbon_g_c -
        state.symbiont_nonstructural.carbon_g_c * host_tissue_carbon_g_c) / tissue_total_g_c;
    transfer.carbon_g_c = inputs.exchange_rate_per_h * carbon_difference_g_c * inputs.timestep_h;
    next.host_nonstructural.carbon_g_c = state.host_nonstructural.carbon_g_c - transfer.carbon_g_c;
    next.symbiont_nonstructural.carbon_g_c = state.symbiont_nonstructural.carbon_g_c + transfer.carbon_g_c;

    const mobile_carbon_total_g_c = next.host_nonstructural.carbon_g_c + next.symbiont_nonstructural.carbon_g_c;
    if (mobile_carbon_total_g_c > inputs.combined_carbon_presence_threshold_g_c) {
        const nitrogen_difference_g_n = (state.host_nonstructural.nitrogen_g_n * next.symbiont_nonstructural.carbon_g_c -
            state.symbiont_nonstructural.nitrogen_g_n * next.host_nonstructural.carbon_g_c) / mobile_carbon_total_g_c;
        transfer.nitrogen_g_n = inputs.exchange_rate_per_h * nitrogen_difference_g_n * inputs.timestep_h;
        const phosphorus_difference_g_p = (state.host_nonstructural.phosphorus_g_p * next.symbiont_nonstructural.carbon_g_c -
            state.symbiont_nonstructural.phosphorus_g_p * next.host_nonstructural.carbon_g_c) / mobile_carbon_total_g_c;
        transfer.phosphorus_g_p = inputs.exchange_rate_per_h * phosphorus_difference_g_p * inputs.timestep_h;
        next.host_nonstructural.nitrogen_g_n = state.host_nonstructural.nitrogen_g_n - transfer.nitrogen_g_n;
        next.host_nonstructural.phosphorus_g_p = state.host_nonstructural.phosphorus_g_p - transfer.phosphorus_g_p;
        next.symbiont_nonstructural.nitrogen_g_n = state.symbiont_nonstructural.nitrogen_g_n + transfer.nitrogen_g_n;
        next.symbiont_nonstructural.phosphorus_g_p = state.symbiont_nonstructural.phosphorus_g_p + transfer.phosphorus_g_p;
    }
    try validatePool(next.host_nonstructural);
    try validatePool(next.symbiont_nonstructural);
    inline for (.{ transfer.carbon_g_c, transfer.nitrogen_g_n, transfer.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteCanopyHostSymbiontExchange;
    return .{ .state = next, .host_to_symbiont = transfer };
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.host_carbon_presence_threshold_g_c,
        inputs.host_tissue_presence_threshold_g_c,
        inputs.combined_carbon_presence_threshold_g_c,
        inputs.initial_bacterial_carbon_g_c_per_m2,
        inputs.horizontal_cell_area_m2,
        inputs.exchange_rate_per_h,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyHostSymbiontExchangeInput;
    if (inputs.horizontal_cell_area_m2 == 0 or inputs.timestep_h == 0)
        return error.InvalidCanopyHostSymbiontExchangeInput;
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyHostSymbiontExchangePool;
}

fn testInputs() Inputs {
    return .{
        .fixation_type = 4,
        .growth_habit = .perennial,
        .host_carbon_presence_threshold_g_c = 1e-12,
        .host_tissue_presence_threshold_g_c = 1e-10,
        .combined_carbon_presence_threshold_g_c = 1e-12,
        .initial_bacterial_carbon_g_c_per_m2 = 1e-4,
        .horizontal_cell_area_m2 = 100,
        .exchange_rate_per_h = 0.2,
        .timestep_h = 0.5,
    };
}

fn testState() BranchState {
    return .{
        .host_nonstructural = .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .symbiont_nonstructural = .{ .carbon_g_c = 0.2, .nitrogen_g_n = 0.002, .phosphorus_g_p = 0.0002 },
        .host_leaf_and_petiole_carbon_g_c = 20,
        .symbiont_structural_carbon_g_c = 2,
        .physiological_maturity_set = false,
    };
}

test "GROSUB FXRN exchange conserves combined mobile C N P" {
    var states = [_]BranchState{testState()};
    const before = states[0];
    var transfers: [1]Pool = undefined;
    try applyAll(&states, &transfers, testInputs());
    inline for (.{ "carbon_g_c", "nitrogen_g_n", "phosphorus_g_p" }) |field_name|
        try std.testing.expectApproxEqAbs(
            @field(before.host_nonstructural, field_name) + @field(before.symbiont_nonstructural, field_name),
            @field(states[0].host_nonstructural, field_name) + @field(states[0].symbiont_nonstructural, field_name),
            1e-14,
        );
    try std.testing.expect(transfers[0].carbon_g_c > 0);
}

test "GROSUB N and P gradients use post-carbon pools" {
    var states = [_]BranchState{testState()};
    const before = states[0];
    var transfers: [1]Pool = undefined;
    const inputs = testInputs();
    try applyAll(&states, &transfers, inputs);
    const carbon_total = states[0].host_nonstructural.carbon_g_c + states[0].symbiont_nonstructural.carbon_g_c;
    const expected_n = inputs.exchange_rate_per_h *
        (before.host_nonstructural.nitrogen_g_n * states[0].symbiont_nonstructural.carbon_g_c - before.symbiont_nonstructural.nitrogen_g_n * states[0].host_nonstructural.carbon_g_c) /
        carbon_total * inputs.timestep_h;
    try std.testing.expectApproxEqAbs(expected_n, transfers[0].nitrogen_g_n, 1e-15);
}

test "GROSUB annual physiological maturity disables host symbiont exchange" {
    var states = [_]BranchState{testState()};
    states[0].physiological_maturity_set = true;
    const before = states;
    var transfers: [1]Pool = undefined;
    var inputs = testInputs();
    inputs.growth_habit = .annual;
    try applyAll(&states, &transfers, inputs);
    try std.testing.expectEqualDeep(before, states);
    try std.testing.expectEqual(Pool{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }, transfers[0]);
}

test "GROSUB exchange runtime sweep is atomic on invalid late branch" {
    var states = [_]BranchState{testState()} ** 41;
    var transfers: [41]Pool = undefined;
    try applyAll(&states, &transfers, testInputs());
    const before_states = states;
    const before_transfers = transfers;
    states[40].host_nonstructural.nitrogen_g_n = std.math.nan(f64);
    try std.testing.expectError(error.InvalidCanopyHostSymbiontExchangePool, applyAll(&states, &transfers, testInputs()));
    // The invalid value preceded the attempted sweep; all other state and every
    // transfer remain untouched by the failed call.
    try std.testing.expectEqualDeep(before_states[0..40], states[0..40]);
    try std.testing.expectEqualDeep(before_transfers, transfers);
}

test "GROSUB root fixation type leaves canopy exchange state unread" {
    var states = [_]BranchState{testState()};
    states[0].host_nonstructural.carbon_g_c = std.math.nan(f64);
    var transfers = [_]Pool{.{ .carbon_g_c = std.math.nan(f64), .nitrogen_g_n = 0, .phosphorus_g_p = 0 }};
    var inputs = testInputs();
    inputs.fixation_type = 3;
    try applyAll(&states, &transfers, inputs);
    try std.testing.expect(std.math.isNan(states[0].host_nonstructural.carbon_g_c));
}
