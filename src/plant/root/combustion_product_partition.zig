const std = @import("std");

pub const ElementalLoss = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Parameters = struct {
    maximum_mineral_nitrogen_fraction: f64,
    minimum_mineral_nitrogen_fraction: f64,
    maximum_mineral_phosphorus_fraction: f64,
    minimum_mineral_phosphorus_fraction: f64,
    negligible_combustion_g_c: f64,
};

pub const Inputs = struct {
    oxygen_unlimited_combustion_g_c: f64,
    aerobic_combustion_g_c: f64,
    anaerobic_combustion_g_c: f64,
    combustion_temperature_response: f64,
    /// Present only in the surface soil layer (STARTQ `RWTRV*`).
    surface_storage_loss: ?ElementalLoss,
    /// Source-ordered root, symbiont, axis, and nodule losses for this layer.
    pool_losses: []const ElementalLoss,
};

pub const Products = struct {
    active: bool,
    carbon_dioxide_g_c: f64,
    methane_g_c: f64,
    charcoal_g_c: f64,
    nitrogen_oxide_g_n: f64,
    ammonium_g_n: f64,
    phosphorus_oxide_g_p: f64,
    dihydrogen_phosphate_g_p: f64,
    charcoal_carbon_fraction: f64,
    aerobic_noncharcoal_carbon_fraction: f64,
    anaerobic_noncharcoal_carbon_fraction: f64,
    mineral_nitrogen_fraction: f64,
    mineral_phosphorus_fraction: f64,
};

fn zeroProducts() Products {
    return .{
        .active = false,
        .carbon_dioxide_g_c = 0,
        .methane_g_c = 0,
        .charcoal_g_c = 0,
        .nitrogen_oxide_g_n = 0,
        .ammonium_g_n = 0,
        .phosphorus_oxide_g_p = 0,
        .dihydrogen_phosphate_g_p = 0,
        .charcoal_carbon_fraction = 0,
        .aerobic_noncharcoal_carbon_fraction = 0,
        .anaerobic_noncharcoal_carbon_fraction = 0,
        .mineral_nitrogen_fraction = 0,
        .mineral_phosphorus_fraction = 0,
    };
}

fn validateLoss(loss: ElementalLoss) !void {
    inline for (std.meta.fields(ElementalLoss)) |field| {
        const value = @field(loss, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteRootCombustionLoss;
        if (value < 0) return error.InvalidRootCombustionLoss;
    }
}

fn addLoss(result: *Products, loss: ElementalLoss) !void {
    result.charcoal_g_c += loss.carbon_g_c * result.charcoal_carbon_fraction;
    result.ammonium_g_n += loss.nitrogen_g_n * result.mineral_nitrogen_fraction;
    result.dihydrogen_phosphate_g_p +=
        loss.phosphorus_g_p * result.mineral_phosphorus_fraction;
    result.carbon_dioxide_g_c += loss.carbon_g_c *
        (1 - result.charcoal_carbon_fraction) *
        result.aerobic_noncharcoal_carbon_fraction;
    result.methane_g_c += loss.carbon_g_c *
        (1 - result.charcoal_carbon_fraction) *
        result.anaerobic_noncharcoal_carbon_fraction;
    result.nitrogen_oxide_g_n +=
        loss.nitrogen_g_n * (1 - result.mineral_nitrogen_fraction);
    result.phosphorus_oxide_g_p +=
        loss.phosphorus_g_p * (1 - result.mineral_phosphorus_fraction);
    inline for (std.meta.fields(Products)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(result.*, field.name)))
            return error.NonFiniteRootCombustionProduct;
    }
}

/// Direct translation of EXTRACT lines 448--516 for one runtime soil layer.
/// The caller supplies root losses in the exact domain/axis order established
/// by GROSUB and includes storage only for the surface layer.
pub fn calculate(inputs: Inputs, parameters: Parameters) !Products {
    inline for (std.meta.fields(Parameters)) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteRootCombustionParameter;
        if (value < 0 or value > 1 and
            !std.mem.eql(u8, field.name, "negligible_combustion_g_c"))
            return error.InvalidRootCombustionParameter;
    }
    if (parameters.minimum_mineral_nitrogen_fraction >
        parameters.maximum_mineral_nitrogen_fraction or
        parameters.minimum_mineral_phosphorus_fraction >
            parameters.maximum_mineral_phosphorus_fraction)
        return error.InvalidRootCombustionParameter;
    inline for (.{
        inputs.oxygen_unlimited_combustion_g_c,
        inputs.aerobic_combustion_g_c,
        inputs.anaerobic_combustion_g_c,
        inputs.combustion_temperature_response,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteRootCombustionInput;
    if (inputs.oxygen_unlimited_combustion_g_c < 0 or
        inputs.aerobic_combustion_g_c < 0 or
        inputs.anaerobic_combustion_g_c < 0 or
        inputs.combustion_temperature_response < 0 or
        inputs.combustion_temperature_response > 1)
        return error.InvalidRootCombustionInput;
    if (inputs.surface_storage_loss) |loss| try validateLoss(loss);
    for (inputs.pool_losses) |loss| try validateLoss(loss);

    const emitted_carbon_g_c =
        inputs.aerobic_combustion_g_c + inputs.anaerobic_combustion_g_c;
    if (inputs.oxygen_unlimited_combustion_g_c <=
        parameters.negligible_combustion_g_c or
        emitted_carbon_g_c <= parameters.negligible_combustion_g_c)
        return zeroProducts();
    if (emitted_carbon_g_c > inputs.oxygen_unlimited_combustion_g_c)
        return error.RootCombustionExceedsOxygenUnlimitedPotential;

    var result = zeroProducts();
    result.active = true;
    result.mineral_nitrogen_fraction =
        parameters.minimum_mineral_nitrogen_fraction +
        (parameters.maximum_mineral_nitrogen_fraction -
            parameters.minimum_mineral_nitrogen_fraction) *
            (1 - inputs.combustion_temperature_response);
    result.mineral_phosphorus_fraction =
        parameters.minimum_mineral_phosphorus_fraction +
        (parameters.maximum_mineral_phosphorus_fraction -
            parameters.minimum_mineral_phosphorus_fraction) *
            (1 - inputs.combustion_temperature_response);
    result.charcoal_carbon_fraction =
        (inputs.oxygen_unlimited_combustion_g_c - emitted_carbon_g_c) /
        inputs.oxygen_unlimited_combustion_g_c;
    result.aerobic_noncharcoal_carbon_fraction =
        inputs.aerobic_combustion_g_c / emitted_carbon_g_c;
    result.anaerobic_noncharcoal_carbon_fraction =
        inputs.anaerobic_combustion_g_c / emitted_carbon_g_c;

    if (inputs.surface_storage_loss) |loss| try addLoss(&result, loss);
    for (inputs.pool_losses) |loss| try addLoss(&result, loss);
    return result;
}

const test_parameters: Parameters = .{
    .maximum_mineral_nitrogen_fraction = 0.8,
    .minimum_mineral_nitrogen_fraction = 0.2,
    .maximum_mineral_phosphorus_fraction = 0.9,
    .minimum_mineral_phosphorus_fraction = 0.3,
    .negligible_combustion_g_c = 1.0e-12,
};

test "EXTRACT root products preserve oxygen and temperature partitions" {
    const result = try calculate(.{
        .oxygen_unlimited_combustion_g_c = 100,
        .aerobic_combustion_g_c = 45,
        .anaerobic_combustion_g_c = 30,
        .combustion_temperature_response = 0.5,
        .surface_storage_loss = .{
            .carbon_g_c = 10,
            .nitrogen_g_n = 1,
            .phosphorus_g_p = 0.1,
        },
        .pool_losses = &.{
            .{ .carbon_g_c = 20, .nitrogen_g_n = 2, .phosphorus_g_p = 0.2 },
            .{ .carbon_g_c = 30, .nitrogen_g_n = 3, .phosphorus_g_p = 0.3 },
        },
    }, test_parameters);
    try std.testing.expect(result.active);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), result.charcoal_carbon_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 27), result.carbon_dioxide_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 18), result.methane_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 15), result.charcoal_g_c, 1e-14);
}

test "EXTRACT root product partition conserves all three elements" {
    const result = try calculate(.{
        .oxygen_unlimited_combustion_g_c = 10,
        .aerobic_combustion_g_c = 4,
        .anaerobic_combustion_g_c = 4,
        .combustion_temperature_response = 0.25,
        .surface_storage_loss = null,
        .pool_losses = &.{
            .{ .carbon_g_c = 8, .nitrogen_g_n = 2, .phosphorus_g_p = 1 },
        },
    }, test_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 8), result.carbon_dioxide_g_c + result.methane_g_c + result.charcoal_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.nitrogen_oxide_g_n + result.ammonium_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), result.phosphorus_oxide_g_p + result.dihydrogen_phosphate_g_p, 1e-14);
}

test "EXTRACT root product gate returns exact zero state" {
    const result = try calculate(.{
        .oxygen_unlimited_combustion_g_c = 10,
        .aerobic_combustion_g_c = 0,
        .anaerobic_combustion_g_c = 0,
        .combustion_temperature_response = 0.5,
        .surface_storage_loss = null,
        .pool_losses = &.{},
    }, test_parameters);
    try std.testing.expectEqual(zeroProducts(), result);
}
