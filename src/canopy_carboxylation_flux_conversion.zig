const std = @import("std");

pub const Inputs = struct {
    production_admitted: bool,
    timestep_h: f64,
    total_co2_fixation_umol_per_s: f64,
    total_carbohydrate_production_umol_per_s: f64,
    c3_fixation_umol_per_s_by_node: []const f64,
    c4_fixation_umol_per_s_by_node: []const f64,
};

pub const Destination = struct {
    c3_fixation_g_c_per_timestep_by_node: []f64,
    c4_fixation_g_c_per_timestep_by_node: []f64,
};

pub const Totals = struct {
    co2_fixation_g_c_per_timestep: f64,
    carbohydrate_production_g_c_per_timestep: f64,
};

/// Exact GROSUB lines 1649--1684 conversion or gated reset.
///
/// The source coefficient is 12e-6 g C/umol CO2 * 3600 s/h = 0.0432.
/// Fixed 25-node work arrays are runtime-sized here. Caller-provided scratch
/// makes publication atomic without allocating inside the numerical kernel.
pub fn convertOrClear(
    inputs: Inputs,
    scratch: Destination,
    destination: Destination,
) !Totals {
    const node_count = inputs.c3_fixation_umol_per_s_by_node.len;
    if (node_count == 0 or
        inputs.c4_fixation_umol_per_s_by_node.len != node_count or
        scratch.c3_fixation_g_c_per_timestep_by_node.len != node_count or
        scratch.c4_fixation_g_c_per_timestep_by_node.len != node_count or
        destination.c3_fixation_g_c_per_timestep_by_node.len != node_count or
        destination.c4_fixation_g_c_per_timestep_by_node.len != node_count)
        return error.CanopyCarboxylationFluxDimensionMismatch;
    inline for (.{
        inputs.timestep_h,
        inputs.total_co2_fixation_umol_per_s,
        inputs.total_carbohydrate_production_umol_per_s,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteCanopyCarboxylationFlux;
    if (inputs.timestep_h <= 0 or
        inputs.total_co2_fixation_umol_per_s < 0 or
        inputs.total_carbohydrate_production_umol_per_s < 0)
        return error.InvalidCanopyCarboxylationFlux;
    for (inputs.c3_fixation_umol_per_s_by_node, inputs.c4_fixation_umol_per_s_by_node) |c3, c4|
        inline for (.{ c3, c4 }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteCanopyCarboxylationFlux;
            if (value < 0) return error.InvalidCanopyCarboxylationFlux;
        };

    if (!inputs.production_admitted) {
        @memset(scratch.c3_fixation_g_c_per_timestep_by_node, 0);
        @memset(scratch.c4_fixation_g_c_per_timestep_by_node, 0);
        @memcpy(
            destination.c3_fixation_g_c_per_timestep_by_node,
            scratch.c3_fixation_g_c_per_timestep_by_node,
        );
        @memcpy(
            destination.c4_fixation_g_c_per_timestep_by_node,
            scratch.c4_fixation_g_c_per_timestep_by_node,
        );
        return .{
            .co2_fixation_g_c_per_timestep = 0,
            .carbohydrate_production_g_c_per_timestep = 0,
        };
    }

    const totals: Totals = .{
        .co2_fixation_g_c_per_timestep = inputs.total_co2_fixation_umol_per_s * 0.0432 * inputs.timestep_h,
        .carbohydrate_production_g_c_per_timestep = inputs.total_carbohydrate_production_umol_per_s * 0.0432 * inputs.timestep_h,
    };
    inline for (.{
        totals.co2_fixation_g_c_per_timestep,
        totals.carbohydrate_production_g_c_per_timestep,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteCanopyCarboxylationFluxResult;
    for (
        inputs.c3_fixation_umol_per_s_by_node,
        inputs.c4_fixation_umol_per_s_by_node,
        scratch.c3_fixation_g_c_per_timestep_by_node,
        scratch.c4_fixation_g_c_per_timestep_by_node,
    ) |c3_umol_per_s, c4_umol_per_s, *c3_g_c, *c4_g_c| {
        c3_g_c.* = c3_umol_per_s * 0.0432 * inputs.timestep_h;
        c4_g_c.* = c4_umol_per_s * 0.0432 * inputs.timestep_h;
        if (!std.math.isFinite(c3_g_c.*) or !std.math.isFinite(c4_g_c.*))
            return error.NonFiniteCanopyCarboxylationFluxResult;
    }
    @memcpy(
        destination.c3_fixation_g_c_per_timestep_by_node,
        scratch.c3_fixation_g_c_per_timestep_by_node,
    );
    @memcpy(
        destination.c4_fixation_g_c_per_timestep_by_node,
        scratch.c4_fixation_g_c_per_timestep_by_node,
    );
    return totals;
}

test "GROSUB converts runtime node fixation in source K order" {
    const allocator = std.testing.allocator;
    const node_count = 31;
    const c3 = try allocator.alloc(f64, node_count);
    defer allocator.free(c3);
    const c4 = try allocator.alloc(f64, node_count);
    defer allocator.free(c4);
    const scratch_c3 = try allocator.alloc(f64, node_count);
    defer allocator.free(scratch_c3);
    const scratch_c4 = try allocator.alloc(f64, node_count);
    defer allocator.free(scratch_c4);
    const result_c3 = try allocator.alloc(f64, node_count);
    defer allocator.free(result_c3);
    const result_c4 = try allocator.alloc(f64, node_count);
    defer allocator.free(result_c4);
    for (c3, c4, 0..) |*c3_rate, *c4_rate, node| {
        c3_rate.* = @floatFromInt(node);
        c4_rate.* = @as(f64, @floatFromInt(node)) + 0.5;
    }
    const totals = try convertOrClear(.{
        .production_admitted = true,
        .timestep_h = 0.5,
        .total_co2_fixation_umol_per_s = 10,
        .total_carbohydrate_production_umol_per_s = 8,
        .c3_fixation_umol_per_s_by_node = c3,
        .c4_fixation_umol_per_s_by_node = c4,
    }, .{
        .c3_fixation_g_c_per_timestep_by_node = scratch_c3,
        .c4_fixation_g_c_per_timestep_by_node = scratch_c4,
    }, .{
        .c3_fixation_g_c_per_timestep_by_node = result_c3,
        .c4_fixation_g_c_per_timestep_by_node = result_c4,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 10 * 0.0432 * 0.5), totals.co2_fixation_g_c_per_timestep, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 8 * 0.0432 * 0.5), totals.carbohydrate_production_g_c_per_timestep, 1.0e-15);
    for (result_c3, result_c4, 0..) |c3_g_c, c4_g_c, node| {
        try std.testing.expectApproxEqAbs(@as(f64, @floatFromInt(node)) * 0.0432 * 0.5, c3_g_c, 1.0e-15);
        try std.testing.expectApproxEqAbs((@as(f64, @floatFromInt(node)) + 0.5) * 0.0432 * 0.5, c4_g_c, 1.0e-15);
    }
}

test "failed production gates clear all source work values" {
    var scratch_c3 = [_]f64{ 9, 9 };
    var scratch_c4 = [_]f64{ 9, 9 };
    var result_c3 = [_]f64{ 7, 7 };
    var result_c4 = [_]f64{ 7, 7 };
    const totals = try convertOrClear(.{
        .production_admitted = false,
        .timestep_h = 1,
        .total_co2_fixation_umol_per_s = 10,
        .total_carbohydrate_production_umol_per_s = 8,
        .c3_fixation_umol_per_s_by_node = &.{ 1, 2 },
        .c4_fixation_umol_per_s_by_node = &.{ 3, 4 },
    }, .{
        .c3_fixation_g_c_per_timestep_by_node = &scratch_c3,
        .c4_fixation_g_c_per_timestep_by_node = &scratch_c4,
    }, .{
        .c3_fixation_g_c_per_timestep_by_node = &result_c3,
        .c4_fixation_g_c_per_timestep_by_node = &result_c4,
    });
    try std.testing.expectEqual(@as(f64, 0), totals.co2_fixation_g_c_per_timestep);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &result_c3);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &result_c4);
}

test "invalid late node leaves published conversion unchanged" {
    var scratch_c3 = [_]f64{ 0, 0 };
    var scratch_c4 = [_]f64{ 0, 0 };
    var result_c3 = [_]f64{ 7, 8 };
    var result_c4 = [_]f64{ 9, 10 };
    try std.testing.expectError(
        error.NonFiniteCanopyCarboxylationFlux,
        convertOrClear(.{
            .production_admitted = true,
            .timestep_h = 1,
            .total_co2_fixation_umol_per_s = 1,
            .total_carbohydrate_production_umol_per_s = 1,
            .c3_fixation_umol_per_s_by_node = &.{ 1, std.math.nan(f64) },
            .c4_fixation_umol_per_s_by_node = &.{ 1, 1 },
        }, .{
            .c3_fixation_g_c_per_timestep_by_node = &scratch_c3,
            .c4_fixation_g_c_per_timestep_by_node = &scratch_c4,
        }, .{
            .c3_fixation_g_c_per_timestep_by_node = &result_c3,
            .c4_fixation_g_c_per_timestep_by_node = &result_c4,
        }),
    );
    try std.testing.expectEqualSlices(f64, &.{ 7, 8 }, &result_c3);
    try std.testing.expectEqualSlices(f64, &.{ 9, 10 }, &result_c4);
}
