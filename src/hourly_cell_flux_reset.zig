const std = @import("std");

pub const ScalarFluxes = struct {
    total_nitrogen_transfer_g_n: f64,
    total_phosphorus_transfer_g_p: f64,
    net_ecosystem_carbon_g_c: f64,
    ground_evaporation_m3: f64,
    plant_evaporation_m3: f64,
    root_water_uptake_m3: f64,
    volatilization_g_n: f64,
    co2_exchange_g_c: f64,
    methane_exchange_g_c: f64,
    oxygen_exchange_g_o: f64,
    dinitrogen_exchange_g_n: f64,
    nitrous_oxide_exchange_g_n: f64,
    ammonia_exchange_g_n: f64,
    ecosystem_respiration_g_c: f64,
    net_radiation_mj: f64,
    latent_heat_mj: f64,
    sensible_heat_mj: f64,
    ground_heat_mj: f64,
    snow_net_radiation_mj: f64,
    snow_latent_heat_mj: f64,
    snow_sensible_heat_mj: f64,
    snow_ground_heat_mj: f64,
    canopy_carbon_g_c: f64,
};

pub const State = struct {
    management_carbon_transfer_g_c: []f64,
    management_nitrogen_transfer_g_n: []f64,
    management_phosphorus_transfer_g_p: []f64,
    scalar: ScalarFluxes,
};

/// Exact cell-level HOUR1 hourly reset from hour1.f:180-222.
pub fn apply(state: *State, first_subhourly_iteration: bool) !bool {
    const count = state.management_carbon_transfer_g_c.len;
    if (state.management_nitrogen_transfer_g_n.len != count or
        state.management_phosphorus_transfer_g_p.len != count)
        return error.HourlyCellFluxResetDimensionMismatch;
    for (state.management_carbon_transfer_g_c) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteHourlyCellFluxState;
    for (state.management_nitrogen_transfer_g_n) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteHourlyCellFluxState;
    for (state.management_phosphorus_transfer_g_p) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteHourlyCellFluxState;
    inline for (@typeInfo(ScalarFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state.scalar, field.name)))
            return error.NonFiniteHourlyCellFluxState;
    if (!first_subhourly_iteration) return false;

    @memset(state.management_carbon_transfer_g_c, 0);
    @memset(state.management_nitrogen_transfer_g_n, 0);
    @memset(state.management_phosphorus_transfer_g_p, 0);
    state.scalar = std.mem.zeroes(ScalarFluxes);
    return true;
}

fn filledScalars(value: f64) ScalarFluxes {
    var result: ScalarFluxes = undefined;
    inline for (@typeInfo(ScalarFluxes).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "first subhour clears runtime management classes and all cell fluxes" {
    var carbon = [_]f64{ 1, 2, 3, 4, 5, 6 };
    var nitrogen = [_]f64{ 1, 2, 3, 4, 5, 6 };
    var phosphorus = [_]f64{ 1, 2, 3, 4, 5, 6 };
    var state: State = .{
        .management_carbon_transfer_g_c = &carbon,
        .management_nitrogen_transfer_g_n = &nitrogen,
        .management_phosphorus_transfer_g_p = &phosphorus,
        .scalar = filledScalars(9),
    };
    try std.testing.expect(try apply(&state, true));
    for (carbon ++ nitrogen ++ phosphorus) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(ScalarFluxes),
        state.scalar,
    );
}

test "later nonlinear iteration retains hourly accumulators" {
    var carbon = [_]f64{1};
    var nitrogen = [_]f64{2};
    var phosphorus = [_]f64{3};
    var state: State = .{
        .management_carbon_transfer_g_c = &carbon,
        .management_nitrogen_transfer_g_n = &nitrogen,
        .management_phosphorus_transfer_g_p = &phosphorus,
        .scalar = filledScalars(4),
    };
    const before_scalar = state.scalar;
    try std.testing.expect(!try apply(&state, false));
    try std.testing.expectEqual(@as(f64, 1), carbon[0]);
    try std.testing.expectEqual(@as(f64, 2), nitrogen[0]);
    try std.testing.expectEqual(@as(f64, 3), phosphorus[0]);
    try std.testing.expectEqualDeep(before_scalar, state.scalar);
}

test "zero runtime management classes still reset scalar fluxes" {
    var empty_c: [0]f64 = .{};
    var empty_n: [0]f64 = .{};
    var empty_p: [0]f64 = .{};
    var state: State = .{
        .management_carbon_transfer_g_c = &empty_c,
        .management_nitrogen_transfer_g_n = &empty_n,
        .management_phosphorus_transfer_g_p = &empty_p,
        .scalar = filledScalars(5),
    };
    try std.testing.expect(try apply(&state, true));
    try std.testing.expectEqualDeep(
        std.mem.zeroes(ScalarFluxes),
        state.scalar,
    );
}

test "dimension mismatch leaves every ledger unchanged" {
    var carbon = [_]f64{ 1, 2 };
    var nitrogen = [_]f64{1};
    var phosphorus = [_]f64{ 3, 4 };
    var state: State = .{
        .management_carbon_transfer_g_c = &carbon,
        .management_nitrogen_transfer_g_n = &nitrogen,
        .management_phosphorus_transfer_g_p = &phosphorus,
        .scalar = filledScalars(5),
    };
    try std.testing.expectError(
        error.HourlyCellFluxResetDimensionMismatch,
        apply(&state, true),
    );
    try std.testing.expectEqual(@as(f64, 1), carbon[0]);
    try std.testing.expectEqual(@as(f64, 5), state.scalar.ground_heat_mj);
}

test "nonfinite late scalar leaves every ledger unchanged" {
    var carbon = [_]f64{1};
    var nitrogen = [_]f64{2};
    var phosphorus = [_]f64{3};
    var state: State = .{
        .management_carbon_transfer_g_c = &carbon,
        .management_nitrogen_transfer_g_n = &nitrogen,
        .management_phosphorus_transfer_g_p = &phosphorus,
        .scalar = filledScalars(5),
    };
    state.scalar.canopy_carbon_g_c = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteHourlyCellFluxState,
        apply(&state, true),
    );
    try std.testing.expectEqual(@as(f64, 1), carbon[0]);
    try std.testing.expectEqual(@as(f64, 2), nitrogen[0]);
    try std.testing.expectEqual(@as(f64, 3), phosphorus[0]);
}
