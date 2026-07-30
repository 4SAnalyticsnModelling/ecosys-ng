const std = @import("std");
const phase_change = @import("soil_water_phase_change.zig");
const retention = @import("soil_water_retention.zig");

pub const SecondaryDomain = struct {
    porous_medium_volume_m3: f64,
    total_water_equivalent_m3: f64,
    unfrozen_pressure_head_m: f64,
    mualem_van_genuchten: retention.MualemVanGenuchtenParameters,
};

pub const Parameters = struct {
    porous_medium_volume_m3: f64,
    total_water_equivalent_m3: f64,
    unfrozen_pressure_head_m: f64,
    gravitational_water_potential_mpa_per_m: f64,
    pure_water_melting_temperature_k: f64,
    dry_solid_heat_capacity_mj_per_k: f64,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
    ice_water_equivalent_heat_capacity_mj_per_m3_k: f64,
    latent_heat_of_fusion_mj_per_m3: f64,
    mualem_van_genuchten: retention.MualemVanGenuchtenParameters,
    secondary_domain: ?SecondaryDomain = null,

    pub fn validate(self: Parameters) !void {
        inline for (@typeInfo(Parameters).@"struct".fields) |field| {
            if (field.type == f64 and !std.math.isFinite(@field(self, field.name)))
                return error.NonFiniteSoilEnthalpyParameter;
        }
        if (self.porous_medium_volume_m3 <= 0 or
            self.total_water_equivalent_m3 < 0 or
            self.unfrozen_pressure_head_m > 0 or
            self.gravitational_water_potential_mpa_per_m <= 0 or
            self.pure_water_melting_temperature_k <= 0 or
            self.dry_solid_heat_capacity_mj_per_k < 0 or
            self.liquid_water_heat_capacity_mj_per_m3_k <= 0 or
            self.ice_water_equivalent_heat_capacity_mj_per_m3_k <= 0 or
            self.latent_heat_of_fusion_mj_per_m3 <= 0)
            return error.InvalidSoilEnthalpyParameter;
        try self.mualem_van_genuchten.validate();
        if (self.secondary_domain) |secondary| {
            inline for (@typeInfo(SecondaryDomain).@"struct".fields) |field| {
                if (field.type == f64 and
                    !std.math.isFinite(@field(secondary, field.name)))
                    return error.NonFiniteSoilEnthalpyParameter;
            }
            if (secondary.porous_medium_volume_m3 <= 0 or
                secondary.total_water_equivalent_m3 < 0 or
                secondary.unfrozen_pressure_head_m > 0)
                return error.InvalidSoilEnthalpyParameter;
            try secondary.mualem_van_genuchten.validate();
        }
    }
};

pub const State = struct {
    temperature_k: f64,
    liquid_water_m3: f64,
    ice_water_equivalent_m3: f64,
    secondary_liquid_water_m3: f64,
    secondary_ice_water_equivalent_m3: f64,
    sensible_heat_capacity_mj_per_k: f64,
    enthalpy_mj: f64,
};

pub const SolverOptions = struct {
    max_iterations: u16,
    absolute_enthalpy_tolerance_mj: f64 = 1.0e-12,
    relative_enthalpy_tolerance: f64 = 1.0e-10,
    minimum_temperature_k: f64 = 173.15,
    maximum_temperature_k: f64 = 373.15,
    initial_temperature_k: ?f64 = null,
};

pub const SolverResult = struct {
    state: State,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

/// Conservative cell enthalpy relative to solid ice at the pure-water
/// melting temperature. The Dall'Amico liquid/ice partition is evaluated at
/// the trial temperature, so latent and sensible energy share one coordinate.
pub fn stateAtTemperature(parameters: Parameters, temperature_k: f64) !State {
    try parameters.validate();
    if (!std.math.isFinite(temperature_k) or temperature_k <= 0)
        return error.InvalidSoilEnthalpyTemperature;
    const equilibrium = try phase_change.dallAmicoEquilibrium(.{
        .temperature_k = temperature_k,
        .total_water_equivalent_m3 = parameters.total_water_equivalent_m3,
        .porous_medium_volume_m3 = parameters.porous_medium_volume_m3,
        .unfrozen_pressure_head_m = parameters.unfrozen_pressure_head_m,
        .gravitational_water_potential_mpa_per_m = parameters.gravitational_water_potential_mpa_per_m,
        .latent_heat_of_fusion_mj_per_m3 = parameters.latent_heat_of_fusion_mj_per_m3,
        .pure_water_melting_temperature_k = parameters.pure_water_melting_temperature_k,
        .mualem_van_genuchten = parameters.mualem_van_genuchten,
    });
    const secondary_equilibrium =
        if (parameters.secondary_domain) |secondary|
            try phase_change.dallAmicoEquilibrium(.{
                .temperature_k = temperature_k,
                .total_water_equivalent_m3 = secondary.total_water_equivalent_m3,
                .porous_medium_volume_m3 = secondary.porous_medium_volume_m3,
                .unfrozen_pressure_head_m = secondary.unfrozen_pressure_head_m,
                .gravitational_water_potential_mpa_per_m = parameters.gravitational_water_potential_mpa_per_m,
                .latent_heat_of_fusion_mj_per_m3 = parameters.latent_heat_of_fusion_mj_per_m3,
                .pure_water_melting_temperature_k = parameters.pure_water_melting_temperature_k,
                .mualem_van_genuchten = secondary.mualem_van_genuchten,
            })
        else
            null;
    const secondary_liquid_water_m3 =
        if (secondary_equilibrium) |state| state.liquid_water_m3 else 0;
    const secondary_ice_water_equivalent_m3 =
        if (secondary_equilibrium) |state|
            state.ice_water_equivalent_m3
        else
            0;
    const heat_capacity_mj_per_k =
        parameters.dry_solid_heat_capacity_mj_per_k +
        parameters.liquid_water_heat_capacity_mj_per_m3_k *
            (equilibrium.liquid_water_m3 +
                secondary_liquid_water_m3) +
        parameters.ice_water_equivalent_heat_capacity_mj_per_m3_k *
            (equilibrium.ice_water_equivalent_m3 +
                secondary_ice_water_equivalent_m3);
    const enthalpy_mj =
        heat_capacity_mj_per_k *
        (temperature_k - parameters.pure_water_melting_temperature_k) +
        parameters.latent_heat_of_fusion_mj_per_m3 *
            (equilibrium.liquid_water_m3 +
                secondary_liquid_water_m3);
    if (!std.math.isFinite(heat_capacity_mj_per_k) or
        heat_capacity_mj_per_k <= 0 or
        !std.math.isFinite(enthalpy_mj))
        return error.NonFiniteSoilEnthalpyState;
    return .{
        .temperature_k = temperature_k,
        .liquid_water_m3 = equilibrium.liquid_water_m3,
        .ice_water_equivalent_m3 = equilibrium.ice_water_equivalent_m3,
        .secondary_liquid_water_m3 = secondary_liquid_water_m3,
        .secondary_ice_water_equivalent_m3 = secondary_ice_water_equivalent_m3,
        .sensible_heat_capacity_mj_per_k = heat_capacity_mj_per_k,
        .enthalpy_mj = enthalpy_mj,
    };
}

fn enthalpyDerivativeMjPerK(
    parameters: Parameters,
    temperature_k: f64,
    state: State,
) !f64 {
    const pressure_head_change_m_per_k =
        parameters.latent_heat_of_fusion_mj_per_m3 /
        (parameters.gravitational_water_potential_mpa_per_m *
            temperature_k);
    var liquid_water_change_m3_per_k: f64 = 0;
    if (state.ice_water_equivalent_m3 > 0 and
        temperature_k < depressedMeltingTemperatureK(
            parameters,
            parameters.unfrozen_pressure_head_m,
        ))
    {
        liquid_water_change_m3_per_k +=
            parameters.porous_medium_volume_m3 *
            try parameters.mualem_van_genuchten.waterCapacityPerM(
                frozenPressureHeadM(
                    parameters,
                    parameters.unfrozen_pressure_head_m,
                    temperature_k,
                ),
            ) *
            pressure_head_change_m_per_k;
    }
    if (parameters.secondary_domain) |secondary| {
        const secondary_melting_temperature_k =
            depressedMeltingTemperatureK(
                parameters,
                secondary.unfrozen_pressure_head_m,
            );
        if (state.secondary_ice_water_equivalent_m3 > 0 and
            temperature_k < secondary_melting_temperature_k)
        {
            liquid_water_change_m3_per_k +=
                secondary.porous_medium_volume_m3 *
                try secondary.mualem_van_genuchten.waterCapacityPerM(
                    frozenPressureHeadM(
                        parameters,
                        secondary.unfrozen_pressure_head_m,
                        temperature_k,
                    ),
                ) *
                parameters.latent_heat_of_fusion_mj_per_m3 /
                (parameters.gravitational_water_potential_mpa_per_m *
                    temperature_k);
        }
    }
    const heat_capacity_change_mj_per_k2 =
        (parameters.liquid_water_heat_capacity_mj_per_m3_k -
            parameters.ice_water_equivalent_heat_capacity_mj_per_m3_k) *
        liquid_water_change_m3_per_k;
    const derivative_mj_per_k =
        state.sensible_heat_capacity_mj_per_k +
        heat_capacity_change_mj_per_k2 *
            (temperature_k - parameters.pure_water_melting_temperature_k) +
        parameters.latent_heat_of_fusion_mj_per_m3 *
            liquid_water_change_m3_per_k;
    if (!std.math.isFinite(derivative_mj_per_k) or derivative_mj_per_k <= 0)
        return error.NonFiniteSoilEnthalpyDerivative;
    return derivative_mj_per_k;
}

fn depressedMeltingTemperatureK(
    parameters: Parameters,
    unfrozen_pressure_head_m: f64,
) f64 {
    const exponent =
        parameters.gravitational_water_potential_mpa_per_m *
        unfrozen_pressure_head_m /
        parameters.latent_heat_of_fusion_mj_per_m3;
    const minimum_exponent =
        @log(std.math.floatMin(f64)) -
        @log(parameters.pure_water_melting_temperature_k);
    return parameters.pure_water_melting_temperature_k *
        @exp(@max(minimum_exponent, exponent));
}

fn frozenPressureHeadM(
    parameters: Parameters,
    unfrozen_pressure_head_m: f64,
    temperature_k: f64,
) f64 {
    const melting_temperature_k =
        depressedMeltingTemperatureK(parameters, unfrozen_pressure_head_m);
    return if (temperature_k < melting_temperature_k)
        unfrozen_pressure_head_m +
            parameters.latent_heat_of_fusion_mj_per_m3 /
                parameters.gravitational_water_potential_mpa_per_m *
                (@log(temperature_k) - @log(melting_temperature_k))
    else
        unfrozen_pressure_head_m;
}

/// Safeguarded Newton-Raphson/Picard inversion used by the coupled spatial
/// heat residual. The runtime iteration ceiling is a convergence limit, not
/// a number of phase or full-model substeps.
pub fn temperatureFromEnthalpy(
    parameters: Parameters,
    target_enthalpy_mj: f64,
    options: SolverOptions,
) !SolverResult {
    try parameters.validate();
    if (!std.math.isFinite(target_enthalpy_mj))
        return error.NonFiniteSoilEnthalpyTarget;
    if (options.max_iterations == 0 or
        !std.math.isFinite(options.absolute_enthalpy_tolerance_mj) or
        options.absolute_enthalpy_tolerance_mj <= 0 or
        !std.math.isFinite(options.relative_enthalpy_tolerance) or
        options.relative_enthalpy_tolerance <= 0 or
        !std.math.isFinite(options.minimum_temperature_k) or
        !std.math.isFinite(options.maximum_temperature_k) or
        options.minimum_temperature_k <= 0 or
        options.maximum_temperature_k <= options.minimum_temperature_k)
        return error.InvalidSoilEnthalpySolverOption;
    if (options.initial_temperature_k) |initial_temperature_k|
        if (!std.math.isFinite(initial_temperature_k) or
            initial_temperature_k < options.minimum_temperature_k or
            initial_temperature_k > options.maximum_temperature_k)
            return error.InvalidSoilEnthalpySolverOption;
    var lower_temperature_k = options.minimum_temperature_k;
    var upper_temperature_k = options.maximum_temperature_k;
    var lower_state = try stateAtTemperature(parameters, lower_temperature_k);
    var upper_state = try stateAtTemperature(parameters, upper_temperature_k);
    if (target_enthalpy_mj < lower_state.enthalpy_mj or
        target_enthalpy_mj > upper_state.enthalpy_mj)
        return error.SoilEnthalpyTargetOutsideTemperatureBracket;
    var temperature_k = options.initial_temperature_k orelse
        lower_temperature_k +
            (upper_temperature_k - lower_temperature_k) *
                (target_enthalpy_mj - lower_state.enthalpy_mj) /
                (upper_state.enthalpy_mj - lower_state.enthalpy_mj);
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        const state = try stateAtTemperature(parameters, temperature_k);
        const residual_mj = state.enthalpy_mj - target_enthalpy_mj;
        const tolerance_mj = options.absolute_enthalpy_tolerance_mj +
            options.relative_enthalpy_tolerance *
                @max(1.0, @abs(target_enthalpy_mj));
        if (@abs(residual_mj) <= tolerance_mj)
            return .{
                .state = state,
                .iterations = iteration + 1,
                .newton_raphson_steps = newton_steps,
                .picard_steps = picard_steps,
            };
        if (residual_mj < 0) {
            lower_temperature_k = temperature_k;
            lower_state = state;
        } else {
            upper_temperature_k = temperature_k;
            upper_state = state;
        }
        if (upper_temperature_k - lower_temperature_k <=
            64.0 * std.math.floatEps(f64) *
                @max(1.0, @abs(temperature_k)))
        {
            const accepted_state =
                if (@abs(lower_state.enthalpy_mj - target_enthalpy_mj) <=
                @abs(upper_state.enthalpy_mj - target_enthalpy_mj))
                    lower_state
                else
                    upper_state;
            return .{
                .state = accepted_state,
                .iterations = iteration + 1,
                .newton_raphson_steps = newton_steps,
                .picard_steps = picard_steps,
            };
        }
        var accepted_newton = false;
        const derivative_mj_per_k =
            try enthalpyDerivativeMjPerK(parameters, temperature_k, state);
        const candidate_temperature_k =
            temperature_k - residual_mj / derivative_mj_per_k;
        if (candidate_temperature_k > lower_temperature_k and
            candidate_temperature_k < upper_temperature_k)
        {
            temperature_k = candidate_temperature_k;
            newton_steps += 1;
            accepted_newton = true;
        }
        if (!accepted_newton) {
            // A rejected Newton step usually means it crossed a phase-change
            // kink. Midpoint Picard contraction moves into the correct
            // branch; regula falsi can remain pinned to the low-slope warm
            // endpoint and exhaust NPH without appreciably shrinking it.
            temperature_k =
                0.5 * (lower_temperature_k + upper_temperature_k);
            picard_steps += 1;
        }
    }
    const final_state = try stateAtTemperature(parameters, temperature_k);
    std.log.err(
        "soil enthalpy Newton-Picard ceiling reached: iterations={d} target_mj={e} state_mj={e} residual_mj={e} temperature_k={e} lower_k={e} upper_k={e}",
        .{
            options.max_iterations,
            target_enthalpy_mj,
            final_state.enthalpy_mj,
            final_state.enthalpy_mj - target_enthalpy_mj,
            temperature_k,
            lower_temperature_k,
            upper_temperature_k,
        },
    );
    return error.SoilEnthalpySolverDidNotConverge;
}

fn testParameters() Parameters {
    return .{
        .porous_medium_volume_m3 = 1,
        .total_water_equivalent_m3 = 0.45,
        .unfrozen_pressure_head_m = -2,
        .gravitational_water_potential_mpa_per_m = 0.00980665,
        .pure_water_melting_temperature_k = 273.15,
        .dry_solid_heat_capacity_mj_per_k = 1.5,
        .liquid_water_heat_capacity_mj_per_m3_k = 4.19,
        .ice_water_equivalent_heat_capacity_mj_per_m3_k = 1.93,
        .latent_heat_of_fusion_mj_per_m3 = 333.7,
        .mualem_van_genuchten = .{
            .residual_water_content_m3_per_m3 = 0.05,
            .saturated_water_content_m3_per_m3 = 0.45,
            .alpha_per_m = 1.6,
            .n = 1.6,
            .saturated_hydraulic_conductivity_m_per_h = 0.01,
        },
    };
}

test "Dall'Amico enthalpy inversion conserves energy and exits early" {
    const parameters = testParameters();
    const expected = try stateAtTemperature(parameters, 268);
    const solved = try temperatureFromEnthalpy(
        parameters,
        expected.enthalpy_mj,
        .{ .max_iterations = 80 },
    );
    try std.testing.expect(solved.iterations < 80);
    try std.testing.expect(solved.newton_raphson_steps +
        solved.picard_steps > 0);
    try std.testing.expectApproxEqAbs(
        expected.temperature_k,
        solved.state.temperature_k,
        1.0e-8,
    );
    try std.testing.expectApproxEqAbs(
        expected.enthalpy_mj,
        solved.state.enthalpy_mj,
        1.0e-9,
    );
    try std.testing.expectApproxEqAbs(
        parameters.total_water_equivalent_m3,
        solved.state.liquid_water_m3 +
            solved.state.ice_water_equivalent_m3,
        1.0e-14,
    );
}

test "enthalpy outside runtime temperature bracket fails explicitly" {
    try std.testing.expectError(
        error.SoilEnthalpyTargetOutsideTemperatureBracket,
        temperatureFromEnthalpy(
            testParameters(),
            -1.0e9,
            .{ .max_iterations = 80 },
        ),
    );
}
