const std = @import("std");

pub const Inputs = struct {
    canopy_surface_temperature_c: f64,
    reference_ammonia_solubility: f64,
    canopy_dry_matter_fraction_g_c_per_g: f64,
    branch_structural_carbon_g_by_branch: []const f64,
    branch_leaf_area_m2_by_branch: []const f64,
    total_canopy_leaf_area_m2: f64,
    branch_nonstructural_nitrogen_g_per_g_c_by_branch: []const f64,
    branch_nonstructural_nitrogen_g_by_branch: []const f64,
    atmospheric_ammonia_g_n_per_m3: f64,
    aerodynamic_resistance_h_per_m: f64,
    stomatal_resistance_h_per_m: f64,
    canopy_radiation_fraction: f64,
    horizontal_cell_area_m2: f64,
    biological_timestep_h_per_step: f64,
    negligible_branch_amount_g_c: f64,
    exchange_pool_limit_fraction: f64,
};

pub const SharedTerms = struct {
    canopy_ammonia_solubility: f64,
    structural_carbon_to_air_ratio: f64,
};

/// UPTAKE.F 1634--1654. Evaluates runtime branch NH3 exchange. A preflight
/// pass guarantees that a late invalid branch cannot partially mutate output.
pub fn calculate(
    inputs: Inputs,
    branch_exchange_g_n_per_step: []f64,
) !SharedTerms {
    try validate(inputs, branch_exchange_g_n_per_step.len);
    const solubility =
        inputs.reference_ammonia_solubility *
        @exp(0.513 - 0.0171 * inputs.canopy_surface_temperature_c);
    const structural_ratio =
        1.0e-4 * inputs.canopy_dry_matter_fraction_g_c_per_g;
    if (!std.math.isFinite(solubility) or solubility <= 0 or
        !std.math.isFinite(structural_ratio))
        return error.InvalidCanopyAmmoniaExchangeResult;

    for (0..branch_exchange_g_n_per_step.len) |branch|
        _ = try branchExchange(inputs, branch, solubility, structural_ratio);
    for (0..branch_exchange_g_n_per_step.len) |branch|
        branch_exchange_g_n_per_step[branch] =
            branchExchange(inputs, branch, solubility, structural_ratio) catch
                unreachable;
    return .{
        .canopy_ammonia_solubility = solubility,
        .structural_carbon_to_air_ratio = structural_ratio,
    };
}

fn branchExchange(
    inputs: Inputs,
    branch: usize,
    solubility: f64,
    structural_ratio: f64,
) !f64 {
    if (inputs.branch_structural_carbon_g_by_branch[branch] <=
        inputs.negligible_branch_amount_g_c or
        inputs.branch_leaf_area_m2_by_branch[branch] <=
            inputs.negligible_branch_amount_g_c or
        inputs.total_canopy_leaf_area_m2 <=
            inputs.negligible_branch_amount_g_c)
        return 0;
    const resistance =
        inputs.aerodynamic_resistance_h_per_m +
        inputs.stomatal_resistance_h_per_m;
    if (resistance <= 0)
        return error.SingularCanopyAmmoniaExchangeResistance;
    const canopy_ammonia = @max(
        0,
        structural_ratio *
            inputs.branch_nonstructural_nitrogen_g_per_g_c_by_branch[branch] /
            solubility,
    );
    const nitrogen_pool = @max(
        0,
        inputs.branch_nonstructural_nitrogen_g_by_branch[branch],
    );
    const unconstrained_exchange =
        (inputs.atmospheric_ammonia_g_n_per_m3 - canopy_ammonia) /
        resistance *
        inputs.canopy_radiation_fraction *
        inputs.horizontal_cell_area_m2 *
        inputs.biological_timestep_h_per_step *
        inputs.branch_leaf_area_m2_by_branch[branch] /
        inputs.total_canopy_leaf_area_m2;
    const pool_limit =
        inputs.exchange_pool_limit_fraction * nitrogen_pool;
    const exchange = @min(
        pool_limit,
        @max(unconstrained_exchange, -pool_limit),
    );
    if (!std.math.isFinite(exchange))
        return error.NonFiniteCanopyAmmoniaExchangeResult;
    return exchange;
}

fn validate(inputs: Inputs, output_count: usize) !void {
    const branch_count = inputs.branch_structural_carbon_g_by_branch.len;
    if (inputs.branch_leaf_area_m2_by_branch.len != branch_count or
        inputs.branch_nonstructural_nitrogen_g_per_g_c_by_branch.len != branch_count or
        inputs.branch_nonstructural_nitrogen_g_by_branch.len != branch_count or
        output_count != branch_count)
        return error.CanopyAmmoniaExchangeDimensionMismatch;
    inline for (.{
        inputs.canopy_surface_temperature_c,
        inputs.reference_ammonia_solubility,
        inputs.canopy_dry_matter_fraction_g_c_per_g,
        inputs.total_canopy_leaf_area_m2,
        inputs.atmospheric_ammonia_g_n_per_m3,
        inputs.aerodynamic_resistance_h_per_m,
        inputs.stomatal_resistance_h_per_m,
        inputs.canopy_radiation_fraction,
        inputs.horizontal_cell_area_m2,
        inputs.biological_timestep_h_per_step,
        inputs.negligible_branch_amount_g_c,
        inputs.exchange_pool_limit_fraction,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidCanopyAmmoniaExchangeInput;
    if (inputs.reference_ammonia_solubility <= 0 or
        inputs.canopy_dry_matter_fraction_g_c_per_g < 0 or
        inputs.total_canopy_leaf_area_m2 < 0 or
        inputs.atmospheric_ammonia_g_n_per_m3 < 0 or
        inputs.aerodynamic_resistance_h_per_m < 0 or
        inputs.stomatal_resistance_h_per_m < 0 or
        inputs.canopy_radiation_fraction < 0 or
        inputs.horizontal_cell_area_m2 <= 0 or
        inputs.biological_timestep_h_per_step < 0 or
        inputs.negligible_branch_amount_g_c < 0 or
        inputs.exchange_pool_limit_fraction < 0)
        return error.InvalidCanopyAmmoniaExchangeInput;
    for (0..branch_count) |branch| {
        inline for (.{
            inputs.branch_structural_carbon_g_by_branch[branch],
            inputs.branch_leaf_area_m2_by_branch[branch],
            inputs.branch_nonstructural_nitrogen_g_per_g_c_by_branch[branch],
            inputs.branch_nonstructural_nitrogen_g_by_branch[branch],
        }) |value|
            if (!std.math.isFinite(value))
                return error.InvalidCanopyAmmoniaExchangeInput;
        if (inputs.branch_structural_carbon_g_by_branch[branch] < 0 or
            inputs.branch_leaf_area_m2_by_branch[branch] < 0)
            return error.InvalidCanopyAmmoniaExchangeInput;
    }
}

fn sourceInputs() Inputs {
    return .{
        .canopy_surface_temperature_c = 25,
        .reference_ammonia_solubility = 1,
        .canopy_dry_matter_fraction_g_c_per_g = 0.2,
        .branch_structural_carbon_g_by_branch = &.{ 10, 10, 0 },
        .branch_leaf_area_m2_by_branch = &.{ 2, 1, 1 },
        .total_canopy_leaf_area_m2 = 3,
        .branch_nonstructural_nitrogen_g_per_g_c_by_branch = &.{ 0.01, 100_000, 1 },
        .branch_nonstructural_nitrogen_g_by_branch = &.{ 2, 1, 1 },
        .atmospheric_ammonia_g_n_per_m3 = 1,
        .aerodynamic_resistance_h_per_m = 1,
        .stomatal_resistance_h_per_m = 1,
        .canopy_radiation_fraction = 0.5,
        .horizontal_cell_area_m2 = 10,
        .biological_timestep_h_per_step = 1,
        .negligible_branch_amount_g_c = 1e-12,
        .exchange_pool_limit_fraction = 0.1,
    };
}

test "UPTAKE branch ammonia exchange preserves source gates and pool clamps" {
    const inputs = sourceInputs();
    var exchange = [_]f64{ 9, 9, 9 };
    const shared = try calculate(inputs, &exchange);
    try std.testing.expect(shared.canopy_ammonia_solubility > 0);
    try std.testing.expectEqual(@as(f64, 2.0e-5), shared.structural_carbon_to_air_ratio);
    try std.testing.expectEqual(@as(f64, 0.2), exchange[0]);
    try std.testing.expectEqual(@as(f64, -0.1), exchange[1]);
    try std.testing.expectEqual(@as(f64, 0), exchange[2]);
}

test "strict leaf area gate zeros exchange" {
    var inputs = sourceInputs();
    inputs.total_canopy_leaf_area_m2 = inputs.negligible_branch_amount_g_c;
    var exchange = [_]f64{ 9, 9, 9 };
    _ = try calculate(inputs, &exchange);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &exchange);
}

test "late invalid branch preserves output transactionally" {
    var nitrogen = [_]f64{ 0.01, std.math.nan(f64), 1 };
    var inputs = sourceInputs();
    inputs.branch_nonstructural_nitrogen_g_per_g_c_by_branch = &nitrogen;
    var exchange = [_]f64{ 9, 9, 9 };
    try std.testing.expectError(
        error.InvalidCanopyAmmoniaExchangeInput,
        calculate(inputs, &exchange),
    );
    try std.testing.expectEqualSlices(f64, &.{ 9, 9, 9 }, &exchange);
}
