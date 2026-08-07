const std = @import("std");

pub const Parameters = struct {
    maximum_mineral_nitrogen_fraction: f64,
    minimum_mineral_nitrogen_fraction: f64,
    maximum_mineral_phosphorus_fraction: f64,
    minimum_mineral_phosphorus_fraction: f64,
    charcoal_carbon_fraction: f64,
    aerobic_noncharcoal_carbon_fraction: f64,
    anaerobic_noncharcoal_carbon_fraction: f64,
};

pub const BranchLosses = struct {
    carbon_g_c: []const f64,
    nitrogen_g_n: []const f64,
    phosphorus_g_p: []const f64,
};

pub const Products = struct {
    carbon_dioxide_g_c: f64,
    methane_g_c: f64,
    charcoal_g_c: f64,
    nitrogen_oxide_g_n: f64,
    ammonium_g_n: f64,
    phosphorus_oxide_g_p: f64,
    dihydrogen_phosphate_g_p: f64,
    mineral_nitrogen_fraction: f64,
    mineral_phosphorus_fraction: f64,
};

fn validateFraction(value: f64) !void {
    if (!std.math.isFinite(value)) return error.NonFiniteShootCombustionParameter;
    if (value < 0 or value > 1) return error.InvalidShootCombustionParameter;
}

/// Direct source-order translation of EXTRACT lines 256--284. Runtime branch
/// losses are partitioned into charcoal, CO2, CH4, gaseous N/P products, and
/// litter mineral N/P. All product quantities are extensive elemental grams.
pub fn calculate(
    losses: BranchLosses,
    combustion_temperature_response: f64,
    parameters: Parameters,
) !Products {
    if (losses.carbon_g_c.len == 0 or
        losses.carbon_g_c.len != losses.nitrogen_g_n.len or
        losses.carbon_g_c.len != losses.phosphorus_g_p.len)
        return error.ShootCombustionBranchDimensionMismatch;
    try validateFraction(combustion_temperature_response);
    inline for (std.meta.fields(Parameters)) |field|
        try validateFraction(@field(parameters, field.name));
    if (parameters.minimum_mineral_nitrogen_fraction >
        parameters.maximum_mineral_nitrogen_fraction or
        parameters.minimum_mineral_phosphorus_fraction >
            parameters.maximum_mineral_phosphorus_fraction)
        return error.InvalidShootCombustionParameter;
    if (@abs(parameters.aerobic_noncharcoal_carbon_fraction +
        parameters.anaerobic_noncharcoal_carbon_fraction - 1) > 1.0e-12)
        return error.InvalidShootCombustionParameter;
    for (losses.carbon_g_c, losses.nitrogen_g_n, losses.phosphorus_g_p) |
        carbon,
        nitrogen,
        phosphorus,
    | inline for (.{ carbon, nitrogen, phosphorus }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteShootCombustionLoss;
        if (value < 0) return error.InvalidShootCombustionLoss;
    };

    const mineral_nitrogen_fraction =
        parameters.minimum_mineral_nitrogen_fraction +
        (parameters.maximum_mineral_nitrogen_fraction -
            parameters.minimum_mineral_nitrogen_fraction) *
            (1 - combustion_temperature_response);
    const gaseous_nitrogen_fraction = 1 - mineral_nitrogen_fraction;
    const mineral_phosphorus_fraction =
        parameters.minimum_mineral_phosphorus_fraction +
        (parameters.maximum_mineral_phosphorus_fraction -
            parameters.minimum_mineral_phosphorus_fraction) *
            (1 - combustion_temperature_response);
    const gaseous_phosphorus_fraction = 1 - mineral_phosphorus_fraction;

    var result: Products = .{
        .carbon_dioxide_g_c = 0,
        .methane_g_c = 0,
        .charcoal_g_c = 0,
        .nitrogen_oxide_g_n = 0,
        .ammonium_g_n = 0,
        .phosphorus_oxide_g_p = 0,
        .dihydrogen_phosphate_g_p = 0,
        .mineral_nitrogen_fraction = mineral_nitrogen_fraction,
        .mineral_phosphorus_fraction = mineral_phosphorus_fraction,
    };
    for (losses.carbon_g_c, losses.nitrogen_g_n, losses.phosphorus_g_p) |
        carbon,
        nitrogen,
        phosphorus,
    | {
        result.charcoal_g_c += carbon * parameters.charcoal_carbon_fraction;
        result.ammonium_g_n += nitrogen * mineral_nitrogen_fraction;
        result.dihydrogen_phosphate_g_p +=
            phosphorus * mineral_phosphorus_fraction;
        result.carbon_dioxide_g_c += carbon *
            (1 - parameters.charcoal_carbon_fraction) *
            parameters.aerobic_noncharcoal_carbon_fraction;
        result.methane_g_c += carbon *
            (1 - parameters.charcoal_carbon_fraction) *
            parameters.anaerobic_noncharcoal_carbon_fraction;
        result.nitrogen_oxide_g_n += nitrogen * gaseous_nitrogen_fraction;
        result.phosphorus_oxide_g_p += phosphorus * gaseous_phosphorus_fraction;
        inline for (std.meta.fields(Products)) |field|
            if (!std.math.isFinite(@field(result, field.name)))
                return error.NonFiniteShootCombustionProduct;
    }
    return result;
}

const test_parameters: Parameters = .{
    .maximum_mineral_nitrogen_fraction = 0.8,
    .minimum_mineral_nitrogen_fraction = 0.2,
    .maximum_mineral_phosphorus_fraction = 0.9,
    .minimum_mineral_phosphorus_fraction = 0.3,
    .charcoal_carbon_fraction = 0.25,
    .aerobic_noncharcoal_carbon_fraction = 0.6,
    .anaerobic_noncharcoal_carbon_fraction = 0.4,
};

test "EXTRACT shoot products preserve source partition and branch order" {
    const result = try calculate(.{
        .carbon_g_c = &.{ 10, 20, 30 },
        .nitrogen_g_n = &.{ 1, 2, 3 },
        .phosphorus_g_p = &.{ 0.5, 1, 1.5 },
    }, 0.5, test_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), result.mineral_nitrogen_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result.mineral_phosphorus_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 15), result.charcoal_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 27), result.carbon_dioxide_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 18), result.methane_g_c, 1e-15);
}

test "EXTRACT shoot product partition conserves C N and P" {
    const result = try calculate(.{
        .carbon_g_c = &.{ 4, 6 },
        .nitrogen_g_n = &.{ 2, 3 },
        .phosphorus_g_p = &.{ 0.4, 0.6 },
    }, 0.25, test_parameters);
    try std.testing.expectApproxEqAbs(
        @as(f64, 10),
        result.carbon_dioxide_g_c + result.methane_g_c + result.charcoal_g_c,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 5),
        result.nitrogen_oxide_g_n + result.ammonium_g_n,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1),
        result.phosphorus_oxide_g_p + result.dihydrogen_phosphate_g_p,
        1e-14,
    );
}

test "EXTRACT shoot product partition rejects invalid loss atomically" {
    try std.testing.expectError(
        error.NonFiniteShootCombustionLoss,
        calculate(.{
            .carbon_g_c = &.{1},
            .nitrogen_g_n = &.{std.math.nan(f64)},
            .phosphorus_g_p = &.{1},
        }, 0.5, test_parameters),
    );
}
