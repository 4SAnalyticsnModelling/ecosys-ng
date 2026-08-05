const std = @import("std");

pub const Fractions = struct {
    mineral_nitrogen_fraction: f64, // FCOMN
    mineral_phosphorus_fraction: f64, // FCOMP
};

pub const Inputs = struct {
    combusted_nitrogen_g_n: f64, // RCMBN/RWTSTDN...
    combusted_phosphorus_g_p: f64, // RCMBP/RWTSTDP...
};

pub const State = struct {
    gaseous_nitrogen_oxide_g_n: f64, // ZOX
    gaseous_phosphorus_oxide_g_p: f64, // POX
    mineral_ammonium_g_n: f64, // Z4M
    mineral_dihydrogen_phosphate_g_p: f64, // P4M
};

/// EXTRACT lines 281--284, 351--354, 387--390, 472--475, 484--487, 498--501,
/// and 513--516. Applies one source-equation partition step to running fire
/// nutrient products in source assignment order.
pub fn publishStep(state: *State, inputs: Inputs, fractions: Fractions) !void {
    inline for (std.meta.fields(Inputs)) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFireNutrientPartitionInput;
    }
    inline for (std.meta.fields(Fractions)) |field| {
        const value = @field(fractions, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidFireNutrientPartitionFraction;
    }
    inline for (std.meta.fields(State)) |field| {
        const value = @field(state.*, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteFireNutrientPartitionState;
    }

    const gaseous_nitrogen = inputs.combusted_nitrogen_g_n *
        (1.0 - fractions.mineral_nitrogen_fraction);
    const gaseous_phosphorus = inputs.combusted_phosphorus_g_p *
        (1.0 - fractions.mineral_phosphorus_fraction);
    const mineral_nitrogen = inputs.combusted_nitrogen_g_n *
        fractions.mineral_nitrogen_fraction;
    const mineral_phosphorus = inputs.combusted_phosphorus_g_p *
        fractions.mineral_phosphorus_fraction;

    const next: State = .{
        .gaseous_nitrogen_oxide_g_n = state.gaseous_nitrogen_oxide_g_n +
            gaseous_nitrogen,
        .gaseous_phosphorus_oxide_g_p = state.gaseous_phosphorus_oxide_g_p + gaseous_phosphorus,
        .mineral_ammonium_g_n = state.mineral_ammonium_g_n +
            mineral_nitrogen,
        .mineral_dihydrogen_phosphate_g_p = state.mineral_dihydrogen_phosphate_g_p + mineral_phosphorus,
    };
    inline for (std.meta.fields(State)) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteFireNutrientPartitionResult;
    state.* = next;
}

test "fire nutrient partition preserves EXTRACT gaseous and mineral splits" {
    var state = std.mem.zeroes(State);
    try publishStep(&state, .{
        .combusted_nitrogen_g_n = 10,
        .combusted_phosphorus_g_p = 6,
    }, .{
        .mineral_nitrogen_fraction = 0.2,
        .mineral_phosphorus_fraction = 0.25,
    });
    try std.testing.expectEqual(@as(f64, 8), state.gaseous_nitrogen_oxide_g_n);
    try std.testing.expectEqual(@as(f64, 4.5), state.gaseous_phosphorus_oxide_g_p);
    try std.testing.expectEqual(@as(f64, 2), state.mineral_ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 1.5), state.mineral_dihydrogen_phosphate_g_p);
}

test "fire nutrient partition accumulates and is atomic on invalid input" {
    var state = std.mem.zeroes(State);
    try publishStep(&state, .{
        .combusted_nitrogen_g_n = 1,
        .combusted_phosphorus_g_p = 2,
    }, .{
        .mineral_nitrogen_fraction = 0.5,
        .mineral_phosphorus_fraction = 0.5,
    });
    const before = state;
    try std.testing.expectError(
        error.InvalidFireNutrientPartitionFraction,
        publishStep(&state, .{
            .combusted_nitrogen_g_n = 1,
            .combusted_phosphorus_g_p = 2,
        }, .{
            .mineral_nitrogen_fraction = 2,
            .mineral_phosphorus_fraction = 0.5,
        }),
    );
    try std.testing.expectEqual(before, state);
}
