const std = @import("std");
const aqueous = @import("root_aqueous_gas_exchange.zig");

pub const GasAmounts = struct {
    carbon_dioxide_g_c: f64,
    oxygen_g_o: f64,
    methane_g_c: f64,
    nitrous_oxide_g_n: f64,
    ammonia_g_n: f64,
    hydrogen_g_h: f64,
};

pub const GasConcentrations = struct {
    carbon_dioxide_g_c_m3: f64,
    oxygen_g_o_m3: f64,
    methane_g_c_m3: f64,
    nitrous_oxide_g_n_m3: f64,
    ammonia_g_n_m3: f64,
    hydrogen_g_h_m3: f64,
};

pub const Conductances = struct {
    carbon_dioxide_m3_per_step: f64,
    oxygen_m3_per_step: f64,
    methane_m3_per_step: f64,
    nitrous_oxide_m3_per_step: f64,
    ammonia_m3_per_step: f64,
    hydrogen_m3_per_step: f64,
};

pub const DissolvedCapacities = struct {
    carbon_dioxide_m3: f64,
    oxygen_m3: f64,
    methane_m3: f64,
    nitrous_oxide_m3: f64,
    ammonia_m3: f64,
    hydrogen_m3: f64,
};

pub const Inputs = struct {
    root_axis_index: usize,
    root_gas_volume_m3: f64,
    root_water_volume_m3: f64,
    negligible_root_gas_volume_m3: f64,
    gas: GasAmounts,
    aqueous_mass: GasAmounts,
    preceding_phase_flux: GasAmounts,
    atmosphere: GasConcentrations,
    atmosphere_conductance: Conductances,
    dissolved_capacity: DissolvedCapacities,
    gas_aqueous_exchange_rate_per_step: f64,
    soil_root_exchange: aqueous.Result,
    root_carbon_dioxide_source_g_c_per_step: f64,
};

pub const Result = struct {
    root_gas_concentration: GasConcentrations,
    atmosphere_flux: GasAmounts,
    phase_flux: GasAmounts,
    gas: GasAmounts,
    aqueous_mass: GasAmounts,
};

/// UPTAKE.F 2434--2530 root gas gate, atmosphere/phase exchanges, and
/// source-ordered root pool updates.
pub fn computeCells(
    cells: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (cells.len != scratch.len or cells.len != destination.len)
        return error.RootGasTransferDimensionMismatch;
    for (cells, scratch) |inputs, *candidate|
        candidate.* = try compute(inputs);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var concentration = std.mem.zeroes(GasConcentrations);
    var atmosphere_flux = std.mem.zeroes(GasAmounts);
    var phase_flux = std.mem.zeroes(GasAmounts);
    const active = inputs.root_axis_index == 1 and
        inputs.root_gas_volume_m3 > inputs.negligible_root_gas_volume_m3;
    if (active) {
        concentration = concentrations(inputs);
        atmosphere_flux = atmosphereFluxes(inputs, concentration);
        phase_flux = try phaseFluxes(inputs, concentration, atmosphere_flux);
    }
    var gas = inputs.gas;
    var water = inputs.aqueous_mass;
    gas.carbon_dioxide_g_c +=
        -phase_flux.carbon_dioxide_g_c + atmosphere_flux.carbon_dioxide_g_c;
    gas.oxygen_g_o += -phase_flux.oxygen_g_o + atmosphere_flux.oxygen_g_o;
    gas.methane_g_c += -phase_flux.methane_g_c + atmosphere_flux.methane_g_c;
    gas.nitrous_oxide_g_n +=
        -phase_flux.nitrous_oxide_g_n + atmosphere_flux.nitrous_oxide_g_n;
    gas.ammonia_g_n += -phase_flux.ammonia_g_n + atmosphere_flux.ammonia_g_n;
    gas.hydrogen_g_h += -phase_flux.hydrogen_g_h + atmosphere_flux.hydrogen_g_h;
    water.carbon_dioxide_g_c += phase_flux.carbon_dioxide_g_c +
        inputs.soil_root_exchange.carbon_dioxide_exchange_g_c_per_step +
        inputs.aqueous_mass.carbon_dioxide_g_c +
        inputs.root_carbon_dioxide_source_g_c_per_step;
    water.oxygen_g_o += phase_flux.oxygen_g_o -
        inputs.soil_root_exchange.oxygen_from_root_g_o_per_step;
    water.methane_g_c += phase_flux.methane_g_c +
        inputs.soil_root_exchange.methane_exchange_g_c_per_step;
    water.nitrous_oxide_g_n += phase_flux.nitrous_oxide_g_n +
        inputs.soil_root_exchange.nitrous_oxide_exchange_g_n_per_step;
    water.ammonia_g_n += phase_flux.ammonia_g_n +
        inputs.soil_root_exchange.ammonia_non_band_exchange_g_n_per_step +
        inputs.soil_root_exchange.ammonia_band_exchange_g_n_per_step;
    water.hydrogen_g_h += phase_flux.hydrogen_g_h +
        inputs.soil_root_exchange.hydrogen_exchange_g_h_per_step;
    const result = Result{
        .root_gas_concentration = concentration,
        .atmosphere_flux = atmosphere_flux,
        .phase_flux = phase_flux,
        .gas = gas,
        .aqueous_mass = water,
    };
    try validateResult(result);
    return result;
}

fn concentrations(inputs: Inputs) GasConcentrations {
    const volume = inputs.root_gas_volume_m3;
    return .{
        .carbon_dioxide_g_c_m3 = @max(
            0,
            (inputs.gas.carbon_dioxide_g_c -
                inputs.preceding_phase_flux.carbon_dioxide_g_c) / volume,
        ),
        .oxygen_g_o_m3 = @max(
            0,
            (inputs.gas.oxygen_g_o - inputs.preceding_phase_flux.oxygen_g_o) / volume,
        ),
        .methane_g_c_m3 = @max(
            0,
            (inputs.gas.methane_g_c - inputs.preceding_phase_flux.methane_g_c) / volume,
        ),
        .nitrous_oxide_g_n_m3 = @max(
            0,
            (inputs.gas.nitrous_oxide_g_n -
                inputs.preceding_phase_flux.nitrous_oxide_g_n) / volume,
        ),
        .ammonia_g_n_m3 = @max(
            0,
            (inputs.gas.ammonia_g_n - inputs.preceding_phase_flux.ammonia_g_n) / volume,
        ),
        .hydrogen_g_h_m3 = @max(
            0,
            (inputs.gas.hydrogen_g_h - inputs.preceding_phase_flux.hydrogen_g_h) / volume,
        ),
    };
}

fn atmosphereFluxes(inputs: Inputs, root: GasConcentrations) GasAmounts {
    const volume = inputs.root_gas_volume_m3;
    const conductance = inputs.atmosphere_conductance;
    return .{
        .carbon_dioxide_g_c = @min(conductance.carbon_dioxide_m3_per_step, volume) *
            (inputs.atmosphere.carbon_dioxide_g_c_m3 - root.carbon_dioxide_g_c_m3),
        .oxygen_g_o = @min(conductance.oxygen_m3_per_step, volume) *
            (inputs.atmosphere.oxygen_g_o_m3 - root.oxygen_g_o_m3),
        .methane_g_c = @min(conductance.methane_m3_per_step, volume) *
            (inputs.atmosphere.methane_g_c_m3 - root.methane_g_c_m3),
        .nitrous_oxide_g_n = @min(conductance.nitrous_oxide_m3_per_step, volume) *
            (inputs.atmosphere.nitrous_oxide_g_n_m3 - root.nitrous_oxide_g_n_m3),
        .ammonia_g_n = @min(conductance.ammonia_m3_per_step, volume) *
            (inputs.atmosphere.ammonia_g_n_m3 - root.ammonia_g_n_m3),
        .hydrogen_g_h = @min(conductance.hydrogen_m3_per_step, volume) *
            (inputs.atmosphere.hydrogen_g_h_m3 - root.hydrogen_g_h_m3),
    };
}

fn phaseFluxes(
    inputs: Inputs,
    root: GasConcentrations,
    atmosphere_flux: GasAmounts,
) !GasAmounts {
    const exchange = inputs.soil_root_exchange;
    const ammonia_uptake = exchange.ammonia_non_band_exchange_g_n_per_step +
        exchange.ammonia_band_exchange_g_n_per_step;
    const carbon_dioxide_after_sources = inputs.aqueous_mass.carbon_dioxide_g_c +
        inputs.root_carbon_dioxide_source_g_c_per_step;
    const oxygen_after_uptake = inputs.aqueous_mass.oxygen_g_o -
        exchange.oxygen_from_root_g_o_per_step;
    const methane_after_uptake = inputs.aqueous_mass.methane_g_c +
        exchange.methane_exchange_g_c_per_step;
    const nitrous_oxide_after_uptake = inputs.aqueous_mass.nitrous_oxide_g_n +
        exchange.nitrous_oxide_exchange_g_n_per_step;
    const ammonia_after_uptake = inputs.aqueous_mass.ammonia_g_n + ammonia_uptake;
    const hydrogen_after_uptake = inputs.aqueous_mass.hydrogen_g_h +
        exchange.hydrogen_exchange_g_h_per_step;
    var result = GasAmounts{
        .carbon_dioxide_g_c = try boundedPhaseFlux(
            inputs,
            root.carbon_dioxide_g_c_m3,
            inputs.dissolved_capacity.carbon_dioxide_m3,
            carbon_dioxide_after_sources,
        ),
        .oxygen_g_o = try boundedPhaseFlux(
            inputs,
            root.oxygen_g_o_m3,
            inputs.dissolved_capacity.oxygen_m3,
            oxygen_after_uptake,
        ),
        .methane_g_c = try boundedPhaseFlux(
            inputs,
            root.methane_g_c_m3,
            inputs.dissolved_capacity.methane_m3,
            methane_after_uptake,
        ),
        .nitrous_oxide_g_n = try boundedPhaseFlux(
            inputs,
            root.nitrous_oxide_g_n_m3,
            inputs.dissolved_capacity.nitrous_oxide_m3,
            nitrous_oxide_after_uptake,
        ),
        .ammonia_g_n = try boundedPhaseFlux(
            inputs,
            root.ammonia_g_n_m3,
            inputs.dissolved_capacity.ammonia_m3,
            ammonia_after_uptake,
        ),
        .hydrogen_g_h = try boundedPhaseFlux(
            inputs,
            root.hydrogen_g_h_m3,
            inputs.dissolved_capacity.hydrogen_m3,
            hydrogen_after_uptake,
        ),
    };
    result.oxygen_g_o = @min(
        result.oxygen_g_o,
        atmosphere_flux.oxygen_g_o + inputs.gas.oxygen_g_o,
    );
    return result;
}

fn boundedPhaseFlux(
    inputs: Inputs,
    root_gas_concentration_g_m3: f64,
    dissolved_capacity_m3: f64,
    aqueous_mass_after_exchange_g: f64,
) !f64 {
    const denominator = dissolved_capacity_m3 + inputs.root_gas_volume_m3;
    if (denominator <= 0) return error.InvalidRootGasPhaseCapacity;
    const unconstrained = inputs.gas_aqueous_exchange_rate_per_step *
        (@max(0, root_gas_concentration_g_m3) * dissolved_capacity_m3 -
            aqueous_mass_after_exchange_g * inputs.root_gas_volume_m3) /
        denominator;
    return @max(-aqueous_mass_after_exchange_g, unconstrained);
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.root_gas_volume_m3,
        inputs.root_water_volume_m3,
        inputs.negligible_root_gas_volume_m3,
        inputs.gas_aqueous_exchange_rate_per_step,
    }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootGasTransferInput;
    inline for (@typeInfo(Conductances).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.atmosphere_conductance, field.name)) or
            @field(inputs.atmosphere_conductance, field.name) < 0)
            return error.InvalidRootGasTransferInput;
    inline for (@typeInfo(DissolvedCapacities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.dissolved_capacity, field.name)) or
            @field(inputs.dissolved_capacity, field.name) < 0)
            return error.InvalidRootGasTransferInput;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(GasAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.atmosphere_flux, field.name)) or
            !std.math.isFinite(@field(result.phase_flux, field.name)) or
            !std.math.isFinite(@field(result.gas, field.name)) or
            !std.math.isFinite(@field(result.aqueous_mass, field.name)))
            return error.NonFiniteRootGasTransferResult;
    }
    inline for (@typeInfo(GasConcentrations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.root_gas_concentration, field.name)))
            return error.NonFiniteRootGasTransferResult;
}

fn testInputs() Inputs {
    var exchange = std.mem.zeroes(aqueous.Result);
    exchange.carbon_dioxide_exchange_g_c_per_step = 0.1;
    exchange.oxygen_from_root_g_o_per_step = 0.2;
    exchange.methane_exchange_g_c_per_step = 0.1;
    exchange.nitrous_oxide_exchange_g_n_per_step = 0.1;
    exchange.ammonia_non_band_exchange_g_n_per_step = 0.06;
    exchange.ammonia_band_exchange_g_n_per_step = 0.04;
    exchange.hydrogen_exchange_g_h_per_step = 0.1;
    return .{
        .root_axis_index = 1,
        .root_gas_volume_m3 = 1,
        .root_water_volume_m3 = 1,
        .negligible_root_gas_volume_m3 = 1e-12,
        .gas = allAmounts(2),
        .aqueous_mass = allAmounts(1),
        .preceding_phase_flux = allAmounts(0),
        .atmosphere = .{
            .carbon_dioxide_g_c_m3 = 0.5,
            .oxygen_g_o_m3 = 3,
            .methane_g_c_m3 = 0.5,
            .nitrous_oxide_g_n_m3 = 0.5,
            .ammonia_g_n_m3 = 0.5,
            .hydrogen_g_h_m3 = 0.5,
        },
        .atmosphere_conductance = .{
            .carbon_dioxide_m3_per_step = 0.1,
            .oxygen_m3_per_step = 0.1,
            .methane_m3_per_step = 0.1,
            .nitrous_oxide_m3_per_step = 0.1,
            .ammonia_m3_per_step = 0.1,
            .hydrogen_m3_per_step = 0.1,
        },
        .dissolved_capacity = .{
            .carbon_dioxide_m3 = 1,
            .oxygen_m3 = 1,
            .methane_m3 = 1,
            .nitrous_oxide_m3 = 1,
            .ammonia_m3 = 1,
            .hydrogen_m3 = 1,
        },
        .gas_aqueous_exchange_rate_per_step = 0.2,
        .soil_root_exchange = exchange,
        .root_carbon_dioxide_source_g_c_per_step = 0.05,
    };
}

fn allAmounts(value: f64) GasAmounts {
    return .{
        .carbon_dioxide_g_c = value,
        .oxygen_g_o = value,
        .methane_g_c = value,
        .nitrous_oxide_g_n = value,
        .ammonia_g_n = value,
        .hydrogen_g_h = value,
    };
}

test "active root preserves oxygen balance and source-order carbon update" {
    const inputs = testInputs();
    const result = try compute(inputs);
    const initial_oxygen = inputs.gas.oxygen_g_o + inputs.aqueous_mass.oxygen_g_o;
    const final_oxygen = result.gas.oxygen_g_o + result.aqueous_mass.oxygen_g_o;
    try std.testing.expectApproxEqAbs(
        initial_oxygen + result.atmosphere_flux.oxygen_g_o -
            inputs.soil_root_exchange.oxygen_from_root_g_o_per_step,
        final_oxygen,
        1e-12,
    );
    const initial_carbon =
        inputs.gas.carbon_dioxide_g_c + inputs.aqueous_mass.carbon_dioxide_g_c;
    const final_carbon =
        result.gas.carbon_dioxide_g_c + result.aqueous_mass.carbon_dioxide_g_c;
    try std.testing.expectApproxEqAbs(
        initial_carbon + result.atmosphere_flux.carbon_dioxide_g_c +
            inputs.soil_root_exchange.carbon_dioxide_exchange_g_c_per_step +
            inputs.aqueous_mass.carbon_dioxide_g_c +
            inputs.root_carbon_dioxide_source_g_c_per_step,
        final_carbon,
        1e-12,
    );
}

test "inactive root axis zeros transfer fluxes but preserves aqueous source update" {
    var inputs = testInputs();
    inputs.root_axis_index = 2;
    const result = try compute(inputs);
    inline for (@typeInfo(GasAmounts).@"struct".fields) |field| {
        try std.testing.expectEqual(@as(f64, 0), @field(result.atmosphere_flux, field.name));
        try std.testing.expectEqual(@as(f64, 0), @field(result.phase_flux, field.name));
    }
    try std.testing.expectApproxEqAbs(
        inputs.aqueous_mass.methane_g_c +
            inputs.soil_root_exchange.methane_exchange_g_c_per_step,
        result.aqueous_mass.methane_g_c,
        1e-12,
    );
}

test "oxygen phase flux retains gas availability cap" {
    var inputs = testInputs();
    inputs.gas.oxygen_g_o = 0.01;
    inputs.atmosphere.oxygen_g_o_m3 = 0;
    inputs.aqueous_mass.oxygen_g_o = 10;
    const result = try compute(inputs);
    try std.testing.expect(result.phase_flux.oxygen_g_o <=
        result.atmosphere_flux.oxygen_g_o + inputs.gas.oxygen_g_o);
}

test "runtime batch failure is atomic" {
    var cells = [_]Inputs{ testInputs(), testInputs() };
    cells[1].dissolved_capacity.oxygen_m3 = -1;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].gas.oxygen_g_o = 41;
    destination[1].gas.oxygen_g_o = 42;
    try std.testing.expectError(
        error.InvalidRootGasTransferInput,
        computeCells(&cells, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].gas.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 42), destination[1].gas.oxygen_g_o);
}
