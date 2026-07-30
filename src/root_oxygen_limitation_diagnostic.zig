const std = @import("std");
const root_porosity = @import("plant_root_porosity.zig");

pub const Inputs = struct {
    oxygen_uptake_g_o_per_step: f64,
    oxygen_demand_g_o_per_step: f64,
    presence_threshold_g_o_per_step: f64,
};

/// UPTAKE.F 3858--3862. Return the plant-wide oxygen-limited uptake divided
/// by oxygen-unlimited demand when demand is strictly above ZEROP; otherwise
/// return zero.
pub fn compute(inputs: Inputs) !f64 {
    return root_porosity.sourceOxygenSatisfaction(
        inputs.oxygen_uptake_g_o_per_step,
        inputs.oxygen_demand_g_o_per_step,
        inputs.presence_threshold_g_o_per_step,
    );
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    scratch: []f64,
    destination: []f64,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootOxygenLimitationDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs);
    @memcpy(destination, scratch);
}

test "oxygen limitation preserves source guarded ratio" {
    try std.testing.expectEqual(
        @as(f64, 0.25),
        try compute(.{
            .oxygen_uptake_g_o_per_step = 2,
            .oxygen_demand_g_o_per_step = 8,
            .presence_threshold_g_o_per_step = 1e-12,
        }),
    );
}

test "demand equal to threshold selects exact source zero branch" {
    try std.testing.expectEqual(
        @as(f64, 0),
        try compute(.{
            .oxygen_uptake_g_o_per_step = 0,
            .oxygen_demand_g_o_per_step = 1e-12,
            .presence_threshold_g_o_per_step = 1e-12,
        }),
    );
}

test "physically impossible uptake above demand fails explicitly" {
    try std.testing.expectError(
        error.InvalidRootPorosityOxygenResult,
        compute(.{
            .oxygen_uptake_g_o_per_step = 2,
            .oxygen_demand_g_o_per_step = 1,
            .presence_threshold_g_o_per_step = 0,
        }),
    );
}

test "later nonfinite runtime input leaves destination unchanged" {
    const inputs = [_]Inputs{
        .{
            .oxygen_uptake_g_o_per_step = 1,
            .oxygen_demand_g_o_per_step = 2,
            .presence_threshold_g_o_per_step = 0,
        },
        .{
            .oxygen_uptake_g_o_per_step = std.math.nan(f64),
            .oxygen_demand_g_o_per_step = 2,
            .presence_threshold_g_o_per_step = 0,
        },
    };
    var scratch: [2]f64 = undefined;
    var destination = [_]f64{ 41, 42 };
    try std.testing.expectError(
        error.NonFiniteRootPorosityOxygenInput,
        computeRuntimeAxes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqualSlices(f64, &[_]f64{ 41, 42 }, &destination);
}
