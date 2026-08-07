const std = @import("std");
const aqueous = @import("aqueous_gas_exchange.zig");
const root_gas = @import("gas_transfer.zig");

pub const SoilRootTotals = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
    methane_g_c_per_step: f64,
    nitrous_oxide_g_n_per_step: f64,
    ammonia_non_band_g_n_per_step: f64,
    ammonia_band_g_n_per_step: f64,
    hydrogen_g_h_per_step: f64,
};

pub const RootAtmosphereTotals = struct {
    phase_flux: root_gas.GasAmounts,
    atmosphere_flux: root_gas.GasAmounts,
};

pub const RootProcessTotals = struct {
    carbon_dioxide_g_c_per_step: f64,
    internal_oxygen_g_o_per_step: f64,
    soil_oxygen_g_o_per_step: f64,
};

pub const Totals = struct {
    soil_root: SoilRootTotals,
    root_atmosphere: RootAtmosphereTotals,
    root_process: RootProcessTotals,
};

pub const Contribution = struct {
    soil_root_exchange: aqueous.Result,
    root_gas_transfer: root_gas.Result,
    /// Legacy RCO2PX after line 2477 (`CO2P1 + prior RCO2PX`).
    root_carbon_dioxide_after_source_g_c_per_step: f64,
};

/// UPTAKE.F 2543--2600 source-ordered accumulations. Each runtime entry is
/// one caller-defined axis tuple, such as cell/layer/species/root/micropore.
pub fn accumulateRuntimeAxes(
    initial: []const Totals,
    contributions: []const Contribution,
    scratch: []Totals,
    destination: []Totals,
) !void {
    if (initial.len != contributions.len or
        initial.len != scratch.len or
        initial.len != destination.len)
        return error.RootGasAccumulationDimensionMismatch;
    for (initial, contributions, scratch) |before, contribution, *candidate|
        candidate.* = try accumulate(before, contribution);
    @memcpy(destination, scratch);
}

pub fn accumulate(before: Totals, contribution: Contribution) !Totals {
    var result = before;
    const soil = contribution.soil_root_exchange;
    result.soil_root.carbon_dioxide_g_c_per_step +=
        soil.carbon_dioxide_exchange_g_c_per_step;
    result.soil_root.oxygen_g_o_per_step += soil.oxygen_from_soil_g_o_per_step;
    result.soil_root.methane_g_c_per_step += soil.methane_exchange_g_c_per_step;
    result.soil_root.nitrous_oxide_g_n_per_step +=
        soil.nitrous_oxide_exchange_g_n_per_step;
    result.soil_root.ammonia_non_band_g_n_per_step +=
        soil.ammonia_non_band_exchange_g_n_per_step;
    result.soil_root.ammonia_band_g_n_per_step +=
        soil.ammonia_band_exchange_g_n_per_step;
    result.soil_root.hydrogen_g_h_per_step += soil.hydrogen_exchange_g_h_per_step;

    const gas = contribution.root_gas_transfer;
    addGasAmounts(&result.root_atmosphere.phase_flux, gas.phase_flux);
    addGasAmounts(&result.root_atmosphere.atmosphere_flux, gas.atmosphere_flux);

    result.root_process.carbon_dioxide_g_c_per_step +=
        contribution.root_carbon_dioxide_after_source_g_c_per_step +
        soil.carbon_dioxide_exchange_g_c_per_step;
    result.root_process.internal_oxygen_g_o_per_step +=
        soil.oxygen_from_root_g_o_per_step;
    result.root_process.soil_oxygen_g_o_per_step +=
        soil.oxygen_from_soil_g_o_per_step;
    try validateTotals(result);
    return result;
}

fn addGasAmounts(total: *root_gas.GasAmounts, addition: root_gas.GasAmounts) void {
    // UPTAKE 2576--2587 order: all six phase fluxes, then all six atmosphere
    // fluxes. This helper is called once for each group in that order.
    total.carbon_dioxide_g_c += addition.carbon_dioxide_g_c;
    total.oxygen_g_o += addition.oxygen_g_o;
    total.methane_g_c += addition.methane_g_c;
    total.nitrous_oxide_g_n += addition.nitrous_oxide_g_n;
    total.ammonia_g_n += addition.ammonia_g_n;
    total.hydrogen_g_h += addition.hydrogen_g_h;
}

fn validateTotals(totals: Totals) !void {
    inline for (@typeInfo(SoilRootTotals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.soil_root, field.name)))
            return error.NonFiniteRootGasAccumulation;
    inline for (@typeInfo(root_gas.GasAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.root_atmosphere.phase_flux, field.name)) or
            !std.math.isFinite(@field(totals.root_atmosphere.atmosphere_flux, field.name)))
            return error.NonFiniteRootGasAccumulation;
    }
    inline for (@typeInfo(RootProcessTotals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.root_process, field.name)))
            return error.NonFiniteRootGasAccumulation;
}

fn testContribution(value: f64) Contribution {
    var soil = std.mem.zeroes(aqueous.Result);
    soil.carbon_dioxide_exchange_g_c_per_step = value;
    soil.oxygen_from_soil_g_o_per_step = value + 1;
    soil.oxygen_from_root_g_o_per_step = value + 2;
    soil.methane_exchange_g_c_per_step = value + 3;
    soil.nitrous_oxide_exchange_g_n_per_step = value + 4;
    soil.ammonia_non_band_exchange_g_n_per_step = value + 5;
    soil.ammonia_band_exchange_g_n_per_step = value + 6;
    soil.hydrogen_exchange_g_h_per_step = value + 7;
    var gas = std.mem.zeroes(root_gas.Result);
    gas.phase_flux = allGas(value + 8);
    gas.atmosphere_flux = allGas(value + 9);
    return .{
        .soil_root_exchange = soil,
        .root_gas_transfer = gas,
        .root_carbon_dioxide_after_source_g_c_per_step = value + 10,
    };
}

fn allGas(value: f64) root_gas.GasAmounts {
    return .{
        .carbon_dioxide_g_c = value,
        .oxygen_g_o = value,
        .methane_g_c = value,
        .nitrous_oxide_g_n = value,
        .ammonia_g_n = value,
        .hydrogen_g_h = value,
    };
}

test "source-ordered accumulation increments every gas total" {
    const before = std.mem.zeroes(Totals);
    const contribution = testContribution(1);
    const result = try accumulate(before, contribution);
    try std.testing.expectEqual(@as(f64, 1), result.soil_root.carbon_dioxide_g_c_per_step);
    try std.testing.expectEqual(@as(f64, 7), result.soil_root.ammonia_band_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 9), result.root_atmosphere.phase_flux.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 10), result.root_atmosphere.atmosphere_flux.hydrogen_g_h);
    try std.testing.expectEqual(@as(f64, 12), result.root_process.carbon_dioxide_g_c_per_step);
    try std.testing.expectEqual(@as(f64, 3), result.root_process.internal_oxygen_g_o_per_step);
    try std.testing.expectEqual(@as(f64, 2), result.root_process.soil_oxygen_g_o_per_step);
}

test "runtime axes retain independent accumulation slots" {
    const initial = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    const contributions = [_]Contribution{ testContribution(1), testContribution(20) };
    var scratch: [2]Totals = undefined;
    var destination: [2]Totals = undefined;
    try accumulateRuntimeAxes(&initial, &contributions, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 1), destination[0].soil_root.carbon_dioxide_g_c_per_step);
    try std.testing.expectEqual(@as(f64, 20), destination[1].soil_root.carbon_dioxide_g_c_per_step);
}

test "later nonfinite contribution leaves destination unchanged" {
    const initial = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    var contributions = [_]Contribution{ testContribution(1), testContribution(2) };
    contributions[1].soil_root_exchange.oxygen_from_soil_g_o_per_step =
        std.math.nan(f64);
    var scratch: [2]Totals = undefined;
    var destination = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    destination[0].root_process.soil_oxygen_g_o_per_step = 41;
    destination[1].root_process.soil_oxygen_g_o_per_step = 42;
    try std.testing.expectError(
        error.NonFiniteRootGasAccumulation,
        accumulateRuntimeAxes(&initial, &contributions, &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].root_process.soil_oxygen_g_o_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].root_process.soil_oxygen_g_o_per_step,
    );
}
