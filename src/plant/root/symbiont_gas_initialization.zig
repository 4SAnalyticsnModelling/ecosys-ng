const std = @import("std");

pub const CompartmentVolumes = struct {
    gas_volume_m3: f64,
    water_volume_m3: f64,
};

pub const AmbientForcing = struct {
    air_temperature_c: f64,
    carbon_dioxide_gas_concentration_g_c_m3: f64,
    carbon_dioxide_aqueous_reference_g_c_m3: f64,
    oxygen_gas_concentration_g_o_m3: f64,
    oxygen_aqueous_reference_g_o_m3: f64,
};

pub const Parameters = struct {
    carbon_dioxide_solubility_multiplier: f64,
    carbon_dioxide_solubility_intercept: f64,
    carbon_dioxide_solubility_temperature_per_c: f64,
    oxygen_solubility_multiplier: f64,
    oxygen_solubility_intercept: f64,
    oxygen_solubility_temperature_per_c: f64,
};

pub const GasState = struct {
    carbon_dioxide_gas_g_c: f64,
    carbon_dioxide_aqueous_g_c: f64,
    carbon_dioxide_air_exchange_g_c_per_timestep: f64,
    carbon_dioxide_phase_exchange_g_c_per_timestep: f64,
    soil_root_carbon_dioxide_exchange_g_c_per_timestep: f64,
    aqueous_carbon_dioxide_production_g_c_per_timestep: f64,
    oxygen_gas_g_o: f64,
    oxygen_aqueous_g_o: f64,
    methane_gas_g_c: f64,
    methane_aqueous_g_c: f64,
    nitrous_oxide_gas_g_n: f64,
    nitrous_oxide_aqueous_g_n: f64,
    ammonia_gas_g_n: f64,
    ammonia_aqueous_g_n: f64,
    hydrogen_gas_g_h: f64,
    hydrogen_aqueous_g_h: f64,
    water_fraction: f64,
};

pub const InitializationError = error{
    StateExtentMismatch,
    NonFiniteInput,
    NegativeVolume,
    NegativeConcentration,
    InvalidParameter,
    NonFiniteResult,
};

/// Translates `startq.f` lines 791--811 for runtime population-layer compartments.
pub fn initialize(
    volumes: []const CompartmentVolumes,
    forcing: AmbientForcing,
    parameters: Parameters,
    states: []GasState,
) InitializationError!void {
    if (volumes.len != states.len) return error.StateExtentMismatch;
    try validate(volumes, forcing, parameters);

    const aqueous_co2_concentration_g_c_m3 =
        parameters.carbon_dioxide_solubility_multiplier *
        @exp(parameters.carbon_dioxide_solubility_intercept +
            parameters.carbon_dioxide_solubility_temperature_per_c *
                forcing.air_temperature_c) *
        forcing.carbon_dioxide_aqueous_reference_g_c_m3;
    const aqueous_oxygen_concentration_g_o_m3 =
        parameters.oxygen_solubility_multiplier *
        @exp(parameters.oxygen_solubility_intercept +
            parameters.oxygen_solubility_temperature_per_c *
                forcing.air_temperature_c) *
        forcing.oxygen_aqueous_reference_g_o_m3;
    if (!std.math.isFinite(aqueous_co2_concentration_g_c_m3) or
        !std.math.isFinite(aqueous_oxygen_concentration_g_o_m3))
    {
        return error.NonFiniteResult;
    }

    for (volumes, states) |volume, *state| {
        state.* = std.mem.zeroes(GasState);
        state.carbon_dioxide_gas_g_c =
            forcing.carbon_dioxide_gas_concentration_g_c_m3 * volume.gas_volume_m3;
        state.carbon_dioxide_aqueous_g_c =
            aqueous_co2_concentration_g_c_m3 * volume.water_volume_m3;
        state.oxygen_gas_g_o =
            forcing.oxygen_gas_concentration_g_o_m3 * volume.gas_volume_m3;
        state.oxygen_aqueous_g_o =
            aqueous_oxygen_concentration_g_o_m3 * volume.water_volume_m3;
        state.water_fraction = 1.0;
        inline for (std.meta.fields(GasState)) |field| {
            if (!std.math.isFinite(@field(state, field.name))) {
                return error.NonFiniteResult;
            }
        }
    }
}

fn validate(
    volumes: []const CompartmentVolumes,
    forcing: AmbientForcing,
    parameters: Parameters,
) InitializationError!void {
    inline for (std.meta.fields(AmbientForcing)) |field| {
        const value = @field(forcing, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (field.name.len > 13 and std.mem.endsWith(u8, field.name, "_m3") and value < 0.0) {
            return error.NegativeConcentration;
        }
    }
    inline for (std.meta.fields(Parameters)) |field| {
        if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteInput;
    }
    if (parameters.carbon_dioxide_solubility_multiplier < 0.0 or
        parameters.oxygen_solubility_multiplier < 0.0)
    {
        return error.InvalidParameter;
    }
    for (volumes) |volume| {
        inline for (std.meta.fields(CompartmentVolumes)) |field| {
            const value = @field(volume, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteInput;
            if (value < 0.0) return error.NegativeVolume;
        }
    }
}

fn legacyParameters() Parameters {
    return .{
        .carbon_dioxide_solubility_multiplier = 0.030,
        .carbon_dioxide_solubility_intercept = -2.621,
        .carbon_dioxide_solubility_temperature_per_c = -0.0317,
        .oxygen_solubility_multiplier = 0.032,
        .oxygen_solubility_intercept = -6.175,
        .oxygen_solubility_temperature_per_c = -0.0211,
    };
}

test "gas and aqueous contents scale by runtime compartment volumes" {
    const volumes = [_]CompartmentVolumes{
        .{ .gas_volume_m3 = 2.0, .water_volume_m3 = 3.0 },
        .{ .gas_volume_m3 = 4.0, .water_volume_m3 = 5.0 },
    };
    const forcing = AmbientForcing{
        .air_temperature_c = 20.0,
        .carbon_dioxide_gas_concentration_g_c_m3 = 0.8,
        .carbon_dioxide_aqueous_reference_g_c_m3 = 1.2,
        .oxygen_gas_concentration_g_o_m3 = 200.0,
        .oxygen_aqueous_reference_g_o_m3 = 220.0,
    };
    var states: [volumes.len]GasState = undefined;
    try initialize(&volumes, forcing, legacyParameters(), &states);

    const aqueous_co2 = 0.030 * @exp(-2.621 - 0.0317 * 20.0) * 1.2;
    try std.testing.expectEqual(@as(f64, 1.6), states[0].carbon_dioxide_gas_g_c);
    try std.testing.expectApproxEqRel(
        aqueous_co2 * 3.0,
        states[0].carbon_dioxide_aqueous_g_c,
        1.0e-14,
    );
    try std.testing.expectEqual(@as(f64, 800.0), states[1].oxygen_gas_g_o);
    try std.testing.expectEqual(@as(f64, 1.0), states[1].water_fraction);
    try std.testing.expectEqual(@as(f64, 0.0), states[1].methane_gas_g_c);
}

test "extent mismatch fails before output mutation" {
    const volumes = [_]CompartmentVolumes{
        .{ .gas_volume_m3 = 1.0, .water_volume_m3 = 1.0 },
    };
    var states: [0]GasState = .{};
    try std.testing.expectError(error.StateExtentMismatch, initialize(
        &volumes,
        std.mem.zeroes(AmbientForcing),
        legacyParameters(),
        &states,
    ));
}
