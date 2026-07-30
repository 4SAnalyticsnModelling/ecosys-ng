const std = @import("std");
const transfer = @import("root_gas_transfer_transaction.zig");
const aqueous = @import("root_aqueous_gas_exchange.zig");

pub const SoilGasMasses = struct {
    carbon_dioxide_g_c: f64,
    oxygen_g_o: f64,
    methane_g_c: f64,
    nitrous_oxide_g_n: f64,
    ammonia_g_n: f64,
    hydrogen_g_h: f64,
};

pub const SoilAqueousMasses = struct {
    carbon_dioxide_g_c: f64,
    oxygen_g_o: f64,
    methane_g_c: f64,
    nitrous_oxide_g_n: f64,
    ammonia_non_band_g_n: f64,
    ammonia_band_g_n: f64,
    hydrogen_g_h: f64,
};

pub const ExternalGasFluxes = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
    preceding_aqueous_oxygen_g_o_per_step: f64,
};

pub const Inputs = struct {
    air_filled_porosity: f64,
    minimum_air_filled_porosity: f64,
    root_fraction: f64,
    gas_aqueous_exchange_rate_per_step: f64,
    root_axis_index: usize,
    negligible_volume_m3: f64,
    prepared: transfer.PreparedTransfer,
    gas: SoilGasMasses,
    aqueous_mass: SoilAqueousMasses,
    root_exchange: aqueous.Result,
    external_flux: ExternalGasFluxes,
};

pub const PhaseFluxes = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
    methane_g_c_per_step: f64,
    nitrous_oxide_g_n_per_step: f64,
    ammonia_non_band_g_n_per_step: f64,
    ammonia_band_g_n_per_step: f64,
    hydrogen_g_h_per_step: f64,
};

pub const Result = struct {
    phase_flux: PhaseFluxes,
    gas: SoilGasMasses,
    aqueous_mass: SoilAqueousMasses,
};

/// UPTAKE.F 2367--2430 soil gas-aqueous exchange and source-ordered pool
/// update. All runtime cells commit only after the complete batch succeeds.
pub fn computeCells(
    cells: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (cells.len != scratch.len or cells.len != destination.len)
        return error.SoilGasPhaseExchangeDimensionMismatch;
    for (cells, scratch) |inputs, *candidate|
        candidate.* = try compute(inputs);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var flux = std.mem.zeroes(PhaseFluxes);
    if (inputs.air_filled_porosity > inputs.minimum_air_filled_porosity) {
        const rate = inputs.root_fraction *
            inputs.gas_aqueous_exchange_rate_per_step;
        flux.carbon_dioxide_g_c_per_step = try phaseExchange(
            rate,
            inputs.gas.carbon_dioxide_g_c,
            inputs.prepared.dissolved_capacity.carbon_dioxide,
            inputs.aqueous_mass.carbon_dioxide_g_c -
                inputs.root_exchange.carbon_dioxide_exchange_g_c_per_step,
            inputs.prepared.soil_allocated_air_m3,
        );
        const root_oxygen_after_preceding_flux =
            inputs.root_exchange.oxygen_from_soil_g_o_per_step -
            inputs.external_flux.preceding_aqueous_oxygen_g_o_per_step;
        flux.oxygen_g_o_per_step = try phaseExchange(
            rate,
            inputs.gas.oxygen_g_o,
            inputs.prepared.dissolved_capacity.oxygen,
            inputs.aqueous_mass.oxygen_g_o - root_oxygen_after_preceding_flux,
            inputs.prepared.soil_allocated_air_m3,
        );
        if (inputs.root_axis_index == 1) {
            flux.methane_g_c_per_step = try phaseExchange(
                rate,
                inputs.gas.methane_g_c,
                inputs.prepared.dissolved_capacity.methane,
                inputs.aqueous_mass.methane_g_c -
                    inputs.root_exchange.methane_exchange_g_c_per_step,
                inputs.prepared.soil_allocated_air_m3,
            );
            flux.nitrous_oxide_g_n_per_step = try phaseExchange(
                rate,
                inputs.gas.nitrous_oxide_g_n,
                inputs.prepared.dissolved_capacity.nitrous_oxide,
                inputs.aqueous_mass.nitrous_oxide_g_n -
                    inputs.root_exchange.nitrous_oxide_exchange_g_n_per_step,
                inputs.prepared.soil_allocated_air_m3,
            );
            flux.ammonia_non_band_g_n_per_step = try ammoniaPhaseExchange(
                rate,
                inputs.gas.ammonia_g_n * inputs.prepared.ammonia_non_band_air_m3 /
                    inputs.prepared.soil_allocated_air_m3,
                inputs.prepared.dissolved_capacity.ammonia,
                inputs.aqueous_mass.ammonia_non_band_g_n -
                    inputs.root_exchange.ammonia_non_band_exchange_g_n_per_step,
                inputs.prepared.ammonia_non_band_air_m3,
                inputs.root_exchange.ammonia_non_band_exchange_g_n_per_step,
                inputs.negligible_volume_m3,
            );
            // UPTAKE line 2389 intentionally uses RUPNSX, not RUPNBX, as
            // the limiter for band NH3.
            flux.ammonia_band_g_n_per_step = try ammoniaPhaseExchange(
                rate,
                inputs.gas.ammonia_g_n * inputs.prepared.ammonia_band_air_m3 /
                    inputs.prepared.soil_allocated_air_m3,
                inputs.prepared.ammonia_band_dissolved_capacity_m3,
                inputs.aqueous_mass.ammonia_band_g_n -
                    inputs.root_exchange.ammonia_band_exchange_g_n_per_step,
                inputs.prepared.ammonia_band_air_m3,
                inputs.root_exchange.ammonia_non_band_exchange_g_n_per_step,
                inputs.negligible_volume_m3,
            );
            flux.hydrogen_g_h_per_step = try phaseExchange(
                rate,
                inputs.gas.hydrogen_g_h,
                inputs.prepared.dissolved_capacity.hydrogen,
                inputs.aqueous_mass.hydrogen_g_h -
                    inputs.root_exchange.hydrogen_exchange_g_h_per_step,
                inputs.prepared.soil_allocated_air_m3,
            );
        }
    }

    var result = Result{
        .phase_flux = flux,
        .gas = inputs.gas,
        .aqueous_mass = inputs.aqueous_mass,
    };
    result.gas.oxygen_g_o =
        inputs.gas.oxygen_g_o - flux.oxygen_g_o_per_step +
        inputs.external_flux.oxygen_g_o_per_step;
    result.aqueous_mass.oxygen_g_o =
        inputs.aqueous_mass.oxygen_g_o + flux.oxygen_g_o_per_step -
        inputs.root_exchange.oxygen_from_soil_g_o_per_step;
    result.gas.carbon_dioxide_g_c =
        inputs.gas.carbon_dioxide_g_c - flux.carbon_dioxide_g_c_per_step +
        inputs.external_flux.carbon_dioxide_g_c_per_step;
    result.aqueous_mass.carbon_dioxide_g_c =
        inputs.aqueous_mass.carbon_dioxide_g_c +
        flux.carbon_dioxide_g_c_per_step -
        inputs.root_exchange.carbon_dioxide_exchange_g_c_per_step;
    result.aqueous_mass.methane_g_c =
        inputs.aqueous_mass.methane_g_c + flux.methane_g_c_per_step -
        inputs.root_exchange.methane_exchange_g_c_per_step;
    result.aqueous_mass.nitrous_oxide_g_n =
        inputs.aqueous_mass.nitrous_oxide_g_n +
        flux.nitrous_oxide_g_n_per_step -
        inputs.root_exchange.nitrous_oxide_exchange_g_n_per_step;
    result.aqueous_mass.ammonia_non_band_g_n =
        inputs.aqueous_mass.ammonia_non_band_g_n +
        flux.ammonia_non_band_g_n_per_step -
        inputs.root_exchange.ammonia_non_band_exchange_g_n_per_step;
    result.aqueous_mass.ammonia_band_g_n =
        inputs.aqueous_mass.ammonia_band_g_n +
        flux.ammonia_band_g_n_per_step -
        inputs.root_exchange.ammonia_band_exchange_g_n_per_step;
    result.aqueous_mass.hydrogen_g_h =
        inputs.aqueous_mass.hydrogen_g_h + flux.hydrogen_g_h_per_step -
        inputs.root_exchange.hydrogen_exchange_g_h_per_step;
    try validateResult(result);
    return result;
}

fn phaseExchange(
    rate_per_step: f64,
    gas_mass_g: f64,
    dissolved_capacity_m3: f64,
    aqueous_mass_after_root_exchange_g: f64,
    air_volume_m3: f64,
) !f64 {
    const denominator = dissolved_capacity_m3 + air_volume_m3;
    if (denominator <= 0) return error.InvalidSoilGasPhaseCapacity;
    return rate_per_step *
        (@max(0, gas_mass_g) * dissolved_capacity_m3 -
            @max(0, aqueous_mass_after_root_exchange_g) * air_volume_m3) /
        denominator;
}

fn ammoniaPhaseExchange(
    rate_per_step: f64,
    partitioned_gas_mass_g_n: f64,
    dissolved_capacity_m3: f64,
    aqueous_mass_after_root_exchange_g_n: f64,
    air_volume_m3: f64,
    limiter_g_n_per_step: f64,
    negligible_volume_m3: f64,
) !f64 {
    if (dissolved_capacity_m3 + air_volume_m3 <= negligible_volume_m3)
        return 0;
    const unconstrained = try phaseExchange(
        rate_per_step,
        partitioned_gas_mass_g_n,
        dissolved_capacity_m3,
        aqueous_mass_after_root_exchange_g_n,
        air_volume_m3,
    );
    return @min(limiter_g_n_per_step, @max(-limiter_g_n_per_step, unconstrained));
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.air_filled_porosity,
        inputs.minimum_air_filled_porosity,
        inputs.root_fraction,
        inputs.gas_aqueous_exchange_rate_per_step,
        inputs.negligible_volume_m3,
    }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilGasPhaseExchangeInput;
    if (inputs.root_fraction > 1)
        return error.InvalidSoilGasPhaseExchangeInput;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(PhaseFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.phase_flux, field.name)))
            return error.NonFiniteSoilGasPhaseExchangeResult;
    inline for (@typeInfo(SoilGasMasses).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.gas, field.name)))
            return error.NonFiniteSoilGasPhaseExchangeResult;
    inline for (@typeInfo(SoilAqueousMasses).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.aqueous_mass, field.name)))
            return error.NonFiniteSoilGasPhaseExchangeResult;
}

fn testInputs() Inputs {
    var prepared = std.mem.zeroes(transfer.PreparedTransfer);
    prepared.soil_allocated_air_m3 = 1;
    prepared.dissolved_capacity = .{
        .carbon_dioxide = 1,
        .oxygen = 1,
        .methane = 1,
        .nitrous_oxide = 1,
        .ammonia = 0.6,
        .hydrogen = 1,
    };
    prepared.ammonia_band_dissolved_capacity_m3 = 0.4;
    prepared.ammonia_non_band_air_m3 = 0.6;
    prepared.ammonia_band_air_m3 = 0.4;
    var root_exchange = std.mem.zeroes(aqueous.Result);
    root_exchange.oxygen_from_soil_g_o_per_step = 0.2;
    root_exchange.carbon_dioxide_exchange_g_c_per_step = 0.1;
    root_exchange.methane_exchange_g_c_per_step = 0.1;
    root_exchange.nitrous_oxide_exchange_g_n_per_step = 0.1;
    root_exchange.ammonia_non_band_exchange_g_n_per_step = 0.1;
    root_exchange.ammonia_band_exchange_g_n_per_step = 0.05;
    root_exchange.hydrogen_exchange_g_h_per_step = 0.1;
    return .{
        .air_filled_porosity = 0.3,
        .minimum_air_filled_porosity = 0.1,
        .root_fraction = 0.5,
        .gas_aqueous_exchange_rate_per_step = 0.2,
        .root_axis_index = 1,
        .negligible_volume_m3 = 1e-12,
        .prepared = prepared,
        .gas = .{
            .carbon_dioxide_g_c = 2,
            .oxygen_g_o = 2,
            .methane_g_c = 2,
            .nitrous_oxide_g_n = 2,
            .ammonia_g_n = 2,
            .hydrogen_g_h = 2,
        },
        .aqueous_mass = .{
            .carbon_dioxide_g_c = 1,
            .oxygen_g_o = 1,
            .methane_g_c = 1,
            .nitrous_oxide_g_n = 1,
            .ammonia_non_band_g_n = 0.6,
            .ammonia_band_g_n = 0.4,
            .hydrogen_g_h = 1,
        },
        .root_exchange = root_exchange,
        .external_flux = .{
            .carbon_dioxide_g_c_per_step = 0.03,
            .oxygen_g_o_per_step = 0.04,
            .preceding_aqueous_oxygen_g_o_per_step = 0.01,
        },
    };
}

test "oxygen and carbon pool updates close phase-exchange balances" {
    const inputs = testInputs();
    const result = try compute(inputs);
    const initial_oxygen = inputs.gas.oxygen_g_o + inputs.aqueous_mass.oxygen_g_o;
    const final_oxygen = result.gas.oxygen_g_o + result.aqueous_mass.oxygen_g_o;
    try std.testing.expectApproxEqAbs(
        initial_oxygen + inputs.external_flux.oxygen_g_o_per_step -
            inputs.root_exchange.oxygen_from_soil_g_o_per_step,
        final_oxygen,
        1e-12,
    );
    const initial_carbon =
        inputs.gas.carbon_dioxide_g_c + inputs.aqueous_mass.carbon_dioxide_g_c;
    const final_carbon =
        result.gas.carbon_dioxide_g_c + result.aqueous_mass.carbon_dioxide_g_c;
    try std.testing.expectApproxEqAbs(
        initial_carbon + inputs.external_flux.carbon_dioxide_g_c_per_step -
            inputs.root_exchange.carbon_dioxide_exchange_g_c_per_step,
        final_carbon,
        1e-12,
    );
}

test "air-filled threshold zeros all phase fluxes" {
    var inputs = testInputs();
    inputs.air_filled_porosity = inputs.minimum_air_filled_porosity;
    const result = try compute(inputs);
    inline for (@typeInfo(PhaseFluxes).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(result.phase_flux, field.name));
}

test "nonprimary root axis retains only carbon dioxide and oxygen exchange" {
    var inputs = testInputs();
    inputs.root_axis_index = 2;
    const result = try compute(inputs);
    try std.testing.expect(result.phase_flux.carbon_dioxide_g_c_per_step != 0);
    try std.testing.expect(result.phase_flux.oxygen_g_o_per_step != 0);
    try std.testing.expectEqual(@as(f64, 0), result.phase_flux.methane_g_c_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.phase_flux.nitrous_oxide_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.phase_flux.ammonia_non_band_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.phase_flux.ammonia_band_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.phase_flux.hydrogen_g_h_per_step);
}

test "runtime batch failure does not partially commit" {
    var cells = [_]Inputs{ testInputs(), testInputs() };
    cells[1].prepared.dissolved_capacity.oxygen = -1;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].gas.oxygen_g_o = 41;
    destination[1].gas.oxygen_g_o = 42;
    try std.testing.expectError(
        error.InvalidSoilGasPhaseCapacity,
        computeCells(&cells, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].gas.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 42), destination[1].gas.oxygen_g_o);
}
