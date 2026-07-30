const std = @import("std");

pub const GrowthHabit = enum {
    herbaceous,
    woody,
};

pub const Parameters = struct {
    reference_leaf_appearance_rate_h_inv: f64,
    woody_maximum_concurrent_nodes: usize,
    first_group_upper_bound: i32,
    second_group_upper_bound: i32,
    first_group_maximum_concurrent_nodes: usize,
    second_group_maximum_concurrent_nodes: usize,
    later_group_maximum_concurrent_nodes: usize,
};

pub const Result = struct {
    perennial_node_scale: f64,
    maximum_concurrent_nodes: usize,
    main_stalk_diameter_m: f64,
};

pub const InitializationError = error{
    NonFiniteInput,
    NonPositiveLeafAppearanceRate,
    InvalidParameters,
    NonFiniteResult,
};

/// Translates STARTQ lines 292-305 for one plant species.
pub fn initialize(
    growth_habit: GrowthHabit,
    plant_group: i32,
    leaf_appearance_rate_h_inv: f64,
    parameters: Parameters,
) InitializationError!Result {
    if (!std.math.isFinite(leaf_appearance_rate_h_inv) or
        !std.math.isFinite(parameters.reference_leaf_appearance_rate_h_inv))
    {
        return error.NonFiniteInput;
    }
    if (parameters.reference_leaf_appearance_rate_h_inv <= 0.0 or
        parameters.woody_maximum_concurrent_nodes == 0 or
        parameters.first_group_maximum_concurrent_nodes == 0 or
        parameters.second_group_maximum_concurrent_nodes == 0 or
        parameters.later_group_maximum_concurrent_nodes == 0 or
        parameters.first_group_upper_bound > parameters.second_group_upper_bound)
    {
        return error.InvalidParameters;
    }

    const result = switch (growth_habit) {
        .woody => blk: {
            if (leaf_appearance_rate_h_inv <= 0.0) {
                return error.NonPositiveLeafAppearanceRate;
            }
            break :blk Result{
                .perennial_node_scale = @max(
                    1.0,
                    parameters.reference_leaf_appearance_rate_h_inv /
                        leaf_appearance_rate_h_inv,
                ),
                .maximum_concurrent_nodes = parameters.woody_maximum_concurrent_nodes,
                .main_stalk_diameter_m = 0.0,
            };
        },
        .herbaceous => Result{
            .perennial_node_scale = 1.0,
            .maximum_concurrent_nodes = if (plant_group <= parameters.first_group_upper_bound)
                parameters.first_group_maximum_concurrent_nodes
            else if (plant_group <= parameters.second_group_upper_bound)
                parameters.second_group_maximum_concurrent_nodes
            else
                parameters.later_group_maximum_concurrent_nodes,
            .main_stalk_diameter_m = 0.0,
        },
    };
    if (!std.math.isFinite(result.perennial_node_scale)) return error.NonFiniteResult;
    return result;
}

fn legacyParameters() Parameters {
    return .{
        .reference_leaf_appearance_rate_h_inv = 0.005,
        .woody_maximum_concurrent_nodes = 24,
        .first_group_upper_bound = 10,
        .second_group_upper_bound = 15,
        .first_group_maximum_concurrent_nodes = 3,
        .second_group_maximum_concurrent_nodes = 4,
        .later_group_maximum_concurrent_nodes = 5,
    };
}

test "woody node scale preserves STARTQ rate ratio and lower bound" {
    const slow = try initialize(.woody, 0, 0.0025, legacyParameters());
    try std.testing.expectEqual(@as(f64, 2.0), slow.perennial_node_scale);
    try std.testing.expectEqual(@as(usize, 24), slow.maximum_concurrent_nodes);

    const fast = try initialize(.woody, 0, 0.01, legacyParameters());
    try std.testing.expectEqual(@as(f64, 1.0), fast.perennial_node_scale);
}

test "herbaceous group boundaries select runtime node limits" {
    const parameters = legacyParameters();
    const first = try initialize(.herbaceous, 10, 0.0, parameters);
    const second = try initialize(.herbaceous, 15, 0.0, parameters);
    const later = try initialize(.herbaceous, 16, 0.0, parameters);

    try std.testing.expectEqual(@as(usize, 3), first.maximum_concurrent_nodes);
    try std.testing.expectEqual(@as(usize, 4), second.maximum_concurrent_nodes);
    try std.testing.expectEqual(@as(usize, 5), later.maximum_concurrent_nodes);
    try std.testing.expectEqual(@as(f64, 0.0), later.main_stalk_diameter_m);
}
