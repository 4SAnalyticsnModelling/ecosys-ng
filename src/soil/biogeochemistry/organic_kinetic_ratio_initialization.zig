const std = @import("std");

pub const Inputs = struct {
    pool_count: usize,
    weighted_residue_pool_count: usize,
    kinetic_component_count: usize,
    carbon_quantity_g_c: []const f64,
    nitrogen_quantity_g_n: []const f64,
    phosphorus_quantity_g_p: []const f64,
    carbon_presence_floor_g_c: []const f64,
    kinetic_carbon_fraction: []const f64,
    nitrogen_allocation_weight: []const f64,
    phosphorus_allocation_weight: []const f64,
    default_nitrogen_to_carbon_g_n_per_g_c: []const f64,
    default_phosphorus_to_carbon_g_p_per_g_c: []const f64,
    fraction_sum_tolerance: f64,
};

pub const State = struct {
    component_nitrogen_to_carbon_g_n_per_g_c: []f64,
    component_phosphorus_to_carbon_g_p_per_g_c: []f64,
    pool_nitrogen_to_carbon_g_n_per_g_c: []f64,
    pool_phosphorus_to_carbon_g_p_per_g_c: []f64,
};

fn product(a: usize, b: usize) !usize {
    return std.math.mul(usize, a, b) catch
        return error.OrganicKineticRatioDimensionOverflow;
}

fn validateDimensions(state: State, inputs: Inputs) !usize {
    if (inputs.pool_count == 0 or inputs.kinetic_component_count == 0 or
        inputs.weighted_residue_pool_count > inputs.pool_count)
        return error.InvalidOrganicKineticRatioDimensions;
    const component_total =
        try product(inputs.pool_count, inputs.kinetic_component_count);
    inline for (.{
        inputs.carbon_quantity_g_c,
        inputs.nitrogen_quantity_g_n,
        inputs.phosphorus_quantity_g_p,
        inputs.carbon_presence_floor_g_c,
        inputs.default_nitrogen_to_carbon_g_n_per_g_c,
        inputs.default_phosphorus_to_carbon_g_p_per_g_c,
    }) |values| {
        if (values.len != inputs.pool_count)
            return error.OrganicKineticRatioDimensionMismatch;
    }
    inline for (.{
        inputs.kinetic_carbon_fraction,
        inputs.nitrogen_allocation_weight,
        inputs.phosphorus_allocation_weight,
        state.component_nitrogen_to_carbon_g_n_per_g_c,
        state.component_phosphorus_to_carbon_g_p_per_g_c,
    }) |values| {
        if (values.len != component_total)
            return error.OrganicKineticRatioDimensionMismatch;
    }
    if (state.pool_nitrogen_to_carbon_g_n_per_g_c.len != inputs.pool_count or
        state.pool_phosphorus_to_carbon_g_p_per_g_c.len != inputs.pool_count)
        return error.OrganicKineticRatioDimensionMismatch;
    return component_total;
}

fn validateInputs(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.fraction_sum_tolerance) or
        inputs.fraction_sum_tolerance < 0.0)
        return error.InvalidOrganicKineticRatioInput;
    inline for (.{
        inputs.carbon_quantity_g_c,
        inputs.nitrogen_quantity_g_n,
        inputs.phosphorus_quantity_g_p,
        inputs.carbon_presence_floor_g_c,
        inputs.kinetic_carbon_fraction,
        inputs.nitrogen_allocation_weight,
        inputs.phosphorus_allocation_weight,
        inputs.default_nitrogen_to_carbon_g_n_per_g_c,
        inputs.default_phosphorus_to_carbon_g_p_per_g_c,
    }) |values| {
        for (values) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteOrganicKineticRatioInput;
            if (value < 0.0) return error.InvalidOrganicKineticRatioInput;
        }
    }
    for (0..inputs.pool_count) |pool| {
        var fraction_sum: f64 = 0.0;
        for (0..inputs.kinetic_component_count) |component| {
            const index = pool * inputs.kinetic_component_count + component;
            fraction_sum += inputs.kinetic_carbon_fraction[index];
        }
        if (!std.math.isFinite(fraction_sum) or
            @abs(fraction_sum - 1.0) > inputs.fraction_sum_tolerance)
            return error.InvalidOrganicKineticFractionSum;
    }
}

const PoolScale = struct {
    nitrogen: f64,
    phosphorus: f64,
};

fn weightedPoolScale(inputs: Inputs, pool: usize) !PoolScale {
    var nitrogen_denominator_g_n: f64 = 0.0;
    var phosphorus_denominator_g_p: f64 = 0.0;
    for (0..inputs.kinetic_component_count) |component| {
        const index = pool * inputs.kinetic_component_count + component;
        nitrogen_denominator_g_n +=
            inputs.carbon_quantity_g_c[pool] *
            inputs.kinetic_carbon_fraction[index] *
            inputs.nitrogen_allocation_weight[index];
        phosphorus_denominator_g_p +=
            inputs.carbon_quantity_g_c[pool] *
            inputs.kinetic_carbon_fraction[index] *
            inputs.phosphorus_allocation_weight[index];
    }
    if (!std.math.isFinite(nitrogen_denominator_g_n) or
        !std.math.isFinite(phosphorus_denominator_g_p))
        return error.NonFiniteOrganicKineticRatioResult;
    if (nitrogen_denominator_g_n <= 0.0 or
        phosphorus_denominator_g_p <= 0.0)
        return error.ZeroOrganicKineticAllocationDenominator;
    const result: PoolScale = .{
        .nitrogen = inputs.nitrogen_quantity_g_n[pool] / nitrogen_denominator_g_n,
        .phosphorus = inputs.phosphorus_quantity_g_p[pool] /
            phosphorus_denominator_g_p,
    };
    if (!std.math.isFinite(result.nitrogen) or
        !std.math.isFinite(result.phosphorus))
        return error.NonFiniteOrganicKineticRatioResult;
    return result;
}

fn componentRatios(
    inputs: Inputs,
    pool: usize,
    component: usize,
) !struct { nitrogen: f64, phosphorus: f64 } {
    const index = pool * inputs.kinetic_component_count + component;
    if (inputs.carbon_quantity_g_c[pool] <=
        inputs.carbon_presence_floor_g_c[pool])
        return .{
            .nitrogen = inputs.default_nitrogen_to_carbon_g_n_per_g_c[pool],
            .phosphorus = inputs.default_phosphorus_to_carbon_g_p_per_g_c[pool],
        };
    if (pool < inputs.weighted_residue_pool_count) {
        const scale = try weightedPoolScale(inputs, pool);
        return .{
            .nitrogen = inputs.nitrogen_allocation_weight[index] * scale.nitrogen,
            .phosphorus = inputs.phosphorus_allocation_weight[index] * scale.phosphorus,
        };
    }
    return .{
        .nitrogen = inputs.nitrogen_quantity_g_n[pool] /
            inputs.carbon_quantity_g_c[pool],
        .phosphorus = inputs.phosphorus_quantity_g_p[pool] /
            inputs.carbon_quantity_g_c[pool],
    };
}

fn validateResults(inputs: Inputs) !void {
    for (0..inputs.pool_count) |pool| {
        var aggregate_nitrogen: f64 = 0.0;
        var aggregate_phosphorus: f64 = 0.0;
        for (0..inputs.kinetic_component_count) |component| {
            const ratios = try componentRatios(inputs, pool, component);
            if (!std.math.isFinite(ratios.nitrogen) or
                !std.math.isFinite(ratios.phosphorus) or
                ratios.nitrogen < 0.0 or ratios.phosphorus < 0.0)
                return error.InvalidOrganicKineticRatioResult;
            const index = pool * inputs.kinetic_component_count + component;
            aggregate_nitrogen +=
                inputs.kinetic_carbon_fraction[index] * ratios.nitrogen;
            aggregate_phosphorus +=
                inputs.kinetic_carbon_fraction[index] * ratios.phosphorus;
        }
        if (!std.math.isFinite(aggregate_nitrogen) or
            !std.math.isFinite(aggregate_phosphorus))
            return error.NonFiniteOrganicKineticRatioResult;
    }
}

/// Exact runtime-dimension translation of legacy `STARTS` lines 1232--1277.
///
/// Input quantities may be masses or concentrations, but C, N, and P must
/// share the same basis within each pool. Output ratios are g N/g C and
/// g P/g C. Pools and biochemical components are stored pool-major.
pub fn initialize(state: State, inputs: Inputs) !void {
    _ = try validateDimensions(state, inputs);
    try validateInputs(inputs);
    try validateResults(inputs);

    for (0..inputs.pool_count) |pool| {
        state.pool_nitrogen_to_carbon_g_n_per_g_c[pool] = 0.0;
        state.pool_phosphorus_to_carbon_g_p_per_g_c[pool] = 0.0;
        for (0..inputs.kinetic_component_count) |component| {
            const ratios = try componentRatios(inputs, pool, component);
            const index = pool * inputs.kinetic_component_count + component;
            state.component_nitrogen_to_carbon_g_n_per_g_c[index] =
                ratios.nitrogen;
            state.component_phosphorus_to_carbon_g_p_per_g_c[index] =
                ratios.phosphorus;
            state.pool_nitrogen_to_carbon_g_n_per_g_c[pool] +=
                inputs.kinetic_carbon_fraction[index] * ratios.nitrogen;
            state.pool_phosphorus_to_carbon_g_p_per_g_c[pool] +=
                inputs.kinetic_carbon_fraction[index] * ratios.phosphorus;
        }
    }
}

test "STARTS weighted residue and direct organic pools conserve ratios" {
    const pool_count = 3;
    const component_count = 2;
    var component_n = [_]f64{9.0} ** (pool_count * component_count);
    var component_p = [_]f64{9.0} ** (pool_count * component_count);
    var pool_n = [_]f64{9.0} ** pool_count;
    var pool_p = [_]f64{9.0} ** pool_count;
    try initialize(.{
        .component_nitrogen_to_carbon_g_n_per_g_c = &component_n,
        .component_phosphorus_to_carbon_g_p_per_g_c = &component_p,
        .pool_nitrogen_to_carbon_g_n_per_g_c = &pool_n,
        .pool_phosphorus_to_carbon_g_p_per_g_c = &pool_p,
    }, .{
        .pool_count = pool_count,
        .weighted_residue_pool_count = 1,
        .kinetic_component_count = component_count,
        .carbon_quantity_g_c = &.{ 100, 50, 0 },
        .nitrogen_quantity_g_n = &.{ 4, 5, 0 },
        .phosphorus_quantity_g_p = &.{ 0.4, 0.5, 0 },
        .carbon_presence_floor_g_c = &.{ 1e-12, 1e-12, 1e-12 },
        .kinetic_carbon_fraction = &.{ 0.25, 0.75, 0.4, 0.6, 0.5, 0.5 },
        .nitrogen_allocation_weight = &.{ 1, 2, 1, 1, 1, 1 },
        .phosphorus_allocation_weight = &.{ 1, 2, 1, 1, 1, 1 },
        .default_nitrogen_to_carbon_g_n_per_g_c = &.{ 0.03, 0.05, 0.16 },
        .default_phosphorus_to_carbon_g_p_per_g_c = &.{ 0.003, 0.005, 0.016 },
        .fraction_sum_tolerance = 1e-12,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.04), pool_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.004), pool_p[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.10), pool_n[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), pool_p[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 0.16), component_n[4]);
    try std.testing.expectEqual(@as(f64, 0.016), pool_p[2]);
}

test "zero weighted denominator fails atomically" {
    var component_n = [_]f64{9.0} ** 2;
    var component_p = [_]f64{9.0} ** 2;
    var pool_n = [_]f64{9.0};
    var pool_p = [_]f64{9.0};
    try std.testing.expectError(
        error.ZeroOrganicKineticAllocationDenominator,
        initialize(.{
            .component_nitrogen_to_carbon_g_n_per_g_c = &component_n,
            .component_phosphorus_to_carbon_g_p_per_g_c = &component_p,
            .pool_nitrogen_to_carbon_g_n_per_g_c = &pool_n,
            .pool_phosphorus_to_carbon_g_p_per_g_c = &pool_p,
        }, .{
            .pool_count = 1,
            .weighted_residue_pool_count = 1,
            .kinetic_component_count = 2,
            .carbon_quantity_g_c = &.{100},
            .nitrogen_quantity_g_n = &.{1},
            .phosphorus_quantity_g_p = &.{0.1},
            .carbon_presence_floor_g_c = &.{1e-12},
            .kinetic_carbon_fraction = &.{ 0.5, 0.5 },
            .nitrogen_allocation_weight = &.{ 0, 0 },
            .phosphorus_allocation_weight = &.{ 0, 0 },
            .default_nitrogen_to_carbon_g_n_per_g_c = &.{0.03},
            .default_phosphorus_to_carbon_g_p_per_g_c = &.{0.003},
            .fraction_sum_tolerance = 1e-12,
        }),
    );
    try std.testing.expectEqualSlices(f64, &.{ 9, 9 }, &component_n);
    try std.testing.expectEqualSlices(f64, &.{9}, &pool_n);
}
