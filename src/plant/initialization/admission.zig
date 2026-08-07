const std = @import("std");

pub const BalanceInitializationOption = enum {
    initialize,
    preserve,
};

pub const PlantStateOrigin = enum {
    fresh_input,
    restored_checkpoint,
};

pub const Action = enum {
    skip,
    perform,
};

pub const Plan = struct {
    reset_mass_balance_ledgers: Action,
    initialize_standing_dead: Action,
    initialize_stalk_geometry: Action,
    initialize_heat_and_water: Action,
};

pub const AdmissionError = error{
    InvalidBalanceInitializationOption,
    NegativeContinuationIndex,
};

/// Parses the STARTQ `DATA(20)` option without case sensitivity.
///
/// The legacy `NO` branch initializes fresh balance state; `YES` preserves
/// state supplied by a continuation source.
pub fn parseBalanceInitializationOption(
    value: []const u8,
) AdmissionError!BalanceInitializationOption {
    if (std.ascii.eqlIgnoreCase(value, "no")) return .initialize;
    if (std.ascii.eqlIgnoreCase(value, "yes")) return .preserve;
    return error.InvalidBalanceInitializationOption;
}

/// Translates the nested admission conditions at `startq.f` lines 637--733.
///
/// Existing balance, standing-dead, stalk, and plant heat/water owners execute
/// the returned actions; this coordinator deliberately owns no duplicate state.
pub fn determinePlan(
    balance_option: BalanceInitializationOption,
    continuation_index: i32,
    state_origin: PlantStateOrigin,
) AdmissionError!Plan {
    if (continuation_index < 0) return error.NegativeContinuationIndex;

    const initialize_balance =
        balance_option == .initialize and continuation_index == 0;
    if (!initialize_balance) {
        return .{
            .reset_mass_balance_ledgers = .skip,
            .initialize_standing_dead = .skip,
            .initialize_stalk_geometry = .skip,
            .initialize_heat_and_water = .skip,
        };
    }

    const initialize_plant_state = state_origin == .fresh_input;
    return .{
        .reset_mass_balance_ledgers = .perform,
        .initialize_standing_dead = if (initialize_plant_state) .perform else .skip,
        .initialize_stalk_geometry = if (initialize_plant_state) .perform else .skip,
        .initialize_heat_and_water = if (initialize_plant_state) .perform else .skip,
    };
}

test "fresh initial run admits every STARTQ initialization owner" {
    const option = try parseBalanceInitializationOption("No");
    const plan = try determinePlan(option, 0, .fresh_input);

    try std.testing.expectEqual(Action.perform, plan.reset_mass_balance_ledgers);
    try std.testing.expectEqual(Action.perform, plan.initialize_standing_dead);
    try std.testing.expectEqual(Action.perform, plan.initialize_stalk_geometry);
    try std.testing.expectEqual(Action.perform, plan.initialize_heat_and_water);
}

test "checkpoint state resets balances but preserves plant state in nested branch" {
    const plan = try determinePlan(.initialize, 0, .restored_checkpoint);
    try std.testing.expectEqual(Action.perform, plan.reset_mass_balance_ledgers);
    try std.testing.expectEqual(Action.skip, plan.initialize_standing_dead);
    try std.testing.expectEqual(Action.skip, plan.initialize_heat_and_water);
}

test "continuation and preserve modes skip the entire outer block" {
    const continuation = try determinePlan(.initialize, 1, .fresh_input);
    const preserve = try determinePlan(
        try parseBalanceInitializationOption("YES"),
        0,
        .fresh_input,
    );
    try std.testing.expectEqual(Action.skip, continuation.reset_mass_balance_ledgers);
    try std.testing.expectEqual(Action.skip, preserve.reset_mass_balance_ledgers);
    try std.testing.expectEqual(Action.skip, preserve.initialize_standing_dead);
}
