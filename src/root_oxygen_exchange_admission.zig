const std = @import("std");

pub const Inputs = struct {
    oxygen_unlimited_respiration_g_c_per_step: f64,
    oxygen_demand_g_o_per_g_c: f64,
    root_aqueous_volume_m3: f64,
    oxygen_competition_fraction: f64,
    negligible_root_amount: f64,
    negligible_oxygen_fraction: f64,
};

pub const Result = struct {
    oxygen_demand_g_o_per_step: f64,
    gas_exchange_active: bool,
};

/// UPTAKE.F 1872--1875. Calculates root oxygen demand before applying the
/// source's three strict gas-exchange admission gates.
pub fn calculate(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidRootOxygenExchangeAdmissionInput;
    if (inputs.oxygen_unlimited_respiration_g_c_per_step < 0 or
        inputs.oxygen_demand_g_o_per_g_c < 0 or
        inputs.root_aqueous_volume_m3 < 0 or
        inputs.oxygen_competition_fraction < 0 or
        inputs.negligible_root_amount < 0 or
        inputs.negligible_oxygen_fraction < 0)
        return error.InvalidRootOxygenExchangeAdmissionInput;
    const oxygen_demand =
        inputs.oxygen_demand_g_o_per_g_c *
        inputs.oxygen_unlimited_respiration_g_c_per_step;
    if (!std.math.isFinite(oxygen_demand))
        return error.NonFiniteRootOxygenDemand;
    return .{
        .oxygen_demand_g_o_per_step = oxygen_demand,
        .gas_exchange_active = inputs.oxygen_unlimited_respiration_g_c_per_step >
            inputs.negligible_root_amount and
            inputs.root_aqueous_volume_m3 >
                inputs.negligible_root_amount and
            inputs.oxygen_competition_fraction >
                inputs.negligible_oxygen_fraction,
    };
}

test "UPTAKE root oxygen demand precedes active admission" {
    const result = try calculate(.{
        .oxygen_unlimited_respiration_g_c_per_step = 3,
        .oxygen_demand_g_o_per_g_c = 2.667,
        .root_aqueous_volume_m3 = 0.2,
        .oxygen_competition_fraction = 0.4,
        .negligible_root_amount = 1e-12,
        .negligible_oxygen_fraction = 1e-12,
    });
    try std.testing.expectEqual(@as(f64, 2.667 * 3), result.oxygen_demand_g_o_per_step);
    try std.testing.expect(result.gas_exchange_active);
}

test "strict equality at each source threshold is inactive" {
    const base = Inputs{
        .oxygen_unlimited_respiration_g_c_per_step = 1,
        .oxygen_demand_g_o_per_g_c = 2.667,
        .root_aqueous_volume_m3 = 1,
        .oxygen_competition_fraction = 1,
        .negligible_root_amount = 0.1,
        .negligible_oxygen_fraction = 0.2,
    };
    var inputs = base;
    inputs.oxygen_unlimited_respiration_g_c_per_step =
        inputs.negligible_root_amount;
    try std.testing.expect(!(try calculate(inputs)).gas_exchange_active);
    inputs = base;
    inputs.root_aqueous_volume_m3 = inputs.negligible_root_amount;
    try std.testing.expect(!(try calculate(inputs)).gas_exchange_active);
    inputs = base;
    inputs.oxygen_competition_fraction =
        inputs.negligible_oxygen_fraction;
    try std.testing.expect(!(try calculate(inputs)).gas_exchange_active);
}

test "non-finite respiration fails explicitly" {
    try std.testing.expectError(
        error.InvalidRootOxygenExchangeAdmissionInput,
        calculate(.{
            .oxygen_unlimited_respiration_g_c_per_step = std.math.nan(f64),
            .oxygen_demand_g_o_per_g_c = 2.667,
            .root_aqueous_volume_m3 = 0.2,
            .oxygen_competition_fraction = 0.4,
            .negligible_root_amount = 1e-12,
            .negligible_oxygen_fraction = 1e-12,
        }),
    );
}
