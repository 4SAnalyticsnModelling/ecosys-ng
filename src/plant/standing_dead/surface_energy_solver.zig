const std = @import("std");
const numerics = @import("../../core/numerics.zig");

pub const SolverPolicy = struct {
    legacy_substeps_nph: u16,
    legacy_inner_iterations_mxn: u16,
    absolute_temperature_tolerance_k: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
};

pub const Inputs = struct {
    initial_surface_temperature_k: f64,
    initial_air_temperature_k: f64,
    combustion_heat_megajoules_per_step: f64,
    canopy_air_heat_capacity_megajoules_per_k: f64,
    combined_area_radiation_fraction: f64,
    intercepted_water_m3: f64,
    retained_precipitation_m3_per_step: f64,
    dry_heat_capacity_megajoules_per_k: f64,
    wet_heat_capacity_megajoules_per_k: f64,
    minimum_energy_heat_capacity_megajoules_per_k: f64,
    absorbed_shortwave_megajoules_per_step: f64,
    absorbed_sky_longwave_megajoules_per_step: f64,
    absorbed_lateral_longwave_megajoules_per_step: f64,
    emitted_longwave_coefficient_megajoules_per_step_k4: f64,
    ground_surface_temperature_k: f64,
    radiation_fraction: f64,
    air_vapor_volume_fraction: f64,
    richardson_coefficient_k: f64,
    minimum_richardson_number: f64,
    maximum_richardson_number: f64,
    minimum_surface_resistance_h_per_m: f64,
    maximum_surface_resistance_h_per_m: f64,
    isothermal_sensible_resistance_h_per_m: f64,
    additional_latent_resistance_h_per_m: f64,
    sensible_conductance_megajoules_per_m_k_step: f64,
    latent_conductance_m2_per_step: f64,
    maximum_water_removal_fraction_per_step: f64,
    latent_heat_megajoules_per_m3: f64,
};

pub const Result = struct {
    surface_temperature_k: f64,
    air_temperature_k: f64,
    intercepted_water_m3: f64,
    emitted_longwave_megajoules_per_step: f64,
    ground_longwave_megajoules_per_step: f64,
    net_radiation_megajoules_per_step: f64,
    evaporation_m3_per_step: f64,
    latent_heat_megajoules_per_step: f64,
    vapor_convective_heat_megajoules_per_step: f64,
    sensible_heat_megajoules_per_step: f64,
    storage_heat_megajoules_per_step: f64,
    final_heat_capacity_megajoules_per_k: f64,
    residual_k: f64,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

const Context = struct {
    inputs: Inputs,
    air_temperature_k: f64,
    water_after_precipitation_m3: f64,
    precipitation_heat_megajoules_per_step: f64,
};

const Fluxes = struct {
    emitted_longwave_megajoules_per_step: f64,
    ground_longwave_megajoules_per_step: f64,
    net_radiation_megajoules_per_step: f64,
    evaporation_m3_per_step: f64,
    latent_heat_megajoules_per_step: f64,
    vapor_convective_heat_megajoules_per_step: f64,
    sensible_heat_megajoules_per_step: f64,
    storage_heat_megajoules_per_step: f64,
    final_heat_capacity_megajoules_per_k: f64,
    fixed_point_temperature_k: f64,
};

/// UPTAKE.F 3990--4167. The legacy NPH full substep loop and MXN inner loop
/// become one local early-exit Newton-Raphson/Picard solve. Their product is
/// retained as the maximum iteration evidence.
pub fn solve(inputs: Inputs, policy: SolverPolicy) !Result {
    try validate(inputs, policy);
    const maximum_iterations_u32 =
        try std.math.mul(u32, policy.legacy_substeps_nph, policy.legacy_inner_iterations_mxn);
    if (maximum_iterations_u32 == 0 or maximum_iterations_u32 > std.math.maxInt(u16))
        return error.InvalidStandingDeadEnergyIterationCeiling;
    const maximum_iterations: u16 = @intCast(maximum_iterations_u32);

    var air_temperature_k = inputs.initial_air_temperature_k;
    if (inputs.canopy_air_heat_capacity_megajoules_per_k > 0 and
        inputs.combined_area_radiation_fraction > 0)
        air_temperature_k += inputs.combustion_heat_megajoules_per_step /
            (inputs.canopy_air_heat_capacity_megajoules_per_k *
                inputs.combined_area_radiation_fraction);
    const water_after_precipitation_m3 =
        inputs.intercepted_water_m3 + inputs.retained_precipitation_m3_per_step;
    const context = Context{
        .inputs = inputs,
        .air_temperature_k = air_temperature_k,
        .water_after_precipitation_m3 = water_after_precipitation_m3,
        .precipitation_heat_megajoules_per_step = inputs.retained_precipitation_m3_per_step * 4.19 *
            inputs.initial_surface_temperature_k,
    };
    const solved = try numerics.newtonPicard(
        context,
        residual,
        derivative,
        picard,
        policy.minimum_temperature_k,
        policy.maximum_temperature_k,
        inputs.initial_surface_temperature_k,
        .{
            .absolute_tolerance = policy.absolute_temperature_tolerance_k,
            .relative_tolerance = policy.relative_tolerance,
            .picard_relaxation = policy.picard_relaxation,
            .residual_scale = @max(1, @abs(inputs.initial_surface_temperature_k)),
            .max_iterations = maximum_iterations,
            .safeguard_with_bracket = true,
        },
    );
    const flux = evaluate(context, solved.root);
    try validateFluxes(flux);
    return .{
        .surface_temperature_k = solved.root,
        .air_temperature_k = air_temperature_k,
        .intercepted_water_m3 = water_after_precipitation_m3 + flux.evaporation_m3_per_step,
        .emitted_longwave_megajoules_per_step = flux.emitted_longwave_megajoules_per_step,
        .ground_longwave_megajoules_per_step = flux.ground_longwave_megajoules_per_step,
        .net_radiation_megajoules_per_step = flux.net_radiation_megajoules_per_step,
        .evaporation_m3_per_step = flux.evaporation_m3_per_step,
        .latent_heat_megajoules_per_step = flux.latent_heat_megajoules_per_step,
        .vapor_convective_heat_megajoules_per_step = flux.vapor_convective_heat_megajoules_per_step,
        .sensible_heat_megajoules_per_step = flux.sensible_heat_megajoules_per_step,
        .storage_heat_megajoules_per_step = flux.storage_heat_megajoules_per_step,
        .final_heat_capacity_megajoules_per_k = flux.final_heat_capacity_megajoules_per_k,
        .residual_k = solved.residual,
        .iterations = solved.iterations,
        .newton_raphson_steps = solved.newton_raphson_steps,
        .picard_steps = solved.picard_steps,
    };
}

pub fn solveRuntimeSpecies(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.StandingDeadEnergyDimensionMismatch;
    for (inputs, scratch) |species_inputs, *candidate|
        candidate.* = try solve(species_inputs, policy);
    @memcpy(destination, scratch);
}

fn residual(context: Context, temperature_k: f64) f64 {
    return temperature_k - evaluate(context, temperature_k).fixed_point_temperature_k;
}

fn derivative(context: Context, temperature_k: f64) f64 {
    const increment = @max(1e-5, @abs(temperature_k) * 1e-7);
    return (residual(context, temperature_k + increment) -
        residual(context, temperature_k - increment)) / (2 * increment);
}

fn picard(context: Context, temperature_k: f64) f64 {
    return evaluate(context, temperature_k).fixed_point_temperature_k;
}

fn evaluate(context: Context, temperature_k: f64) Fluxes {
    const inputs = context.inputs;
    const emitted = inputs.emitted_longwave_coefficient_megajoules_per_step_k4 *
        std.math.pow(f64, temperature_k, 4);
    const ground = inputs.emitted_longwave_coefficient_megajoules_per_step_k4 *
        (std.math.pow(f64, temperature_k, 4) -
            std.math.pow(f64, inputs.ground_surface_temperature_k, 4)) *
        inputs.radiation_fraction;
    const net_radiation = inputs.absorbed_shortwave_megajoules_per_step +
        inputs.absorbed_sky_longwave_megajoules_per_step +
        inputs.absorbed_lateral_longwave_megajoules_per_step - emitted - ground;
    const temperature_difference_k = context.air_temperature_k - temperature_k;
    const richardson = std.math.clamp(
        inputs.richardson_coefficient_k / context.air_temperature_k *
            temperature_difference_k,
        inputs.minimum_richardson_number,
        inputs.maximum_richardson_number,
    );
    const stability = 1 - 10 * richardson;
    const resistance = @min(
        inputs.maximum_surface_resistance_h_per_m,
        @max(
            inputs.minimum_surface_resistance_h_per_m,
            inputs.isothermal_sensible_resistance_h_per_m / stability,
        ),
    );
    const sensible_conductance =
        inputs.sensible_conductance_megajoules_per_m_k_step / resistance;
    const latent_conductance =
        inputs.latent_conductance_m2_per_step /
        (resistance + inputs.additional_latent_resistance_h_per_m);
    const saturated_vapor = 2.173e-3 / temperature_k * 0.61 *
        @exp(5360 * (3.661e-3 - 1 / temperature_k));
    const water_flux = latent_conductance *
        (inputs.air_vapor_volume_fraction - saturated_vapor);
    const evaporation = if (water_flux > 0)
        water_flux
    else
        @max(
            water_flux,
            -@max(
                0,
                context.water_after_precipitation_m3 *
                    inputs.maximum_water_removal_fraction_per_step,
            ),
        );
    const latent_heat = evaporation * inputs.latent_heat_megajoules_per_m3;
    const vapor_heat = evaporation * 4.19 * temperature_k;
    const sensible_heat = sensible_conductance *
        (context.air_temperature_k - temperature_k);
    const storage_heat = net_radiation + latent_heat + sensible_heat +
        vapor_heat + context.precipitation_heat_megajoules_per_step;
    const final_capacity = inputs.wet_heat_capacity_megajoules_per_k +
        4.19 * (evaporation + inputs.retained_precipitation_m3_per_step);
    const fixed_point = if (final_capacity > inputs.minimum_energy_heat_capacity_megajoules_per_k)
        (temperature_k * inputs.wet_heat_capacity_megajoules_per_k + storage_heat) /
            final_capacity
    else
        temperature_k;
    return .{
        .emitted_longwave_megajoules_per_step = emitted,
        .ground_longwave_megajoules_per_step = ground,
        .net_radiation_megajoules_per_step = net_radiation,
        .evaporation_m3_per_step = evaporation,
        .latent_heat_megajoules_per_step = latent_heat,
        .vapor_convective_heat_megajoules_per_step = vapor_heat,
        .sensible_heat_megajoules_per_step = sensible_heat,
        .storage_heat_megajoules_per_step = storage_heat,
        .final_heat_capacity_megajoules_per_k = final_capacity,
        .fixed_point_temperature_k = fixed_point,
    };
}

fn validate(inputs: Inputs, policy: SolverPolicy) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteStandingDeadEnergyInput;
    inline for (@typeInfo(SolverPolicy).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(policy, field.name)))
            return error.NonFiniteStandingDeadEnergyPolicy;
    if (inputs.initial_surface_temperature_k <= 0 or
        inputs.initial_air_temperature_k <= 0 or
        inputs.canopy_air_heat_capacity_megajoules_per_k < 0 or
        inputs.combined_area_radiation_fraction < 0 or
        inputs.intercepted_water_m3 < 0 or
        inputs.retained_precipitation_m3_per_step < 0 or
        inputs.dry_heat_capacity_megajoules_per_k < 0 or
        inputs.wet_heat_capacity_megajoules_per_k <= 0 or
        inputs.minimum_energy_heat_capacity_megajoules_per_k < 0 or
        inputs.ground_surface_temperature_k <= 0 or
        inputs.radiation_fraction < 0 or inputs.radiation_fraction > 1 or
        inputs.air_vapor_volume_fraction < 0 or
        inputs.minimum_richardson_number > inputs.maximum_richardson_number or
        inputs.minimum_surface_resistance_h_per_m <= 0 or
        inputs.maximum_surface_resistance_h_per_m <
            inputs.minimum_surface_resistance_h_per_m or
        inputs.isothermal_sensible_resistance_h_per_m <= 0 or
        inputs.additional_latent_resistance_h_per_m < 0 or
        inputs.sensible_conductance_megajoules_per_m_k_step < 0 or
        inputs.latent_conductance_m2_per_step < 0 or
        inputs.maximum_water_removal_fraction_per_step < 0 or
        inputs.maximum_water_removal_fraction_per_step > 1 or
        inputs.latent_heat_megajoules_per_m3 <= 0 or
        policy.absolute_temperature_tolerance_k <= 0 or
        policy.relative_tolerance < 0 or policy.picard_relaxation <= 0 or
        policy.picard_relaxation > 1 or policy.minimum_temperature_k <= 0 or
        policy.minimum_temperature_k >= policy.maximum_temperature_k)
        return error.InvalidStandingDeadEnergyInput;
}

fn validateFluxes(flux: Fluxes) !void {
    inline for (@typeInfo(Fluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteStandingDeadEnergyResult;
    if (flux.final_heat_capacity_megajoules_per_k <= 0)
        return error.InvalidStandingDeadEnergyResult;
}

fn testInputs() Inputs {
    return .{
        .initial_surface_temperature_k = 280,
        .initial_air_temperature_k = 282,
        .combustion_heat_megajoules_per_step = 0,
        .canopy_air_heat_capacity_megajoules_per_k = 2,
        .combined_area_radiation_fraction = 0.5,
        .intercepted_water_m3 = 0.01,
        .retained_precipitation_m3_per_step = 0,
        .dry_heat_capacity_megajoules_per_k = 1,
        .wet_heat_capacity_megajoules_per_k = 1.1,
        .minimum_energy_heat_capacity_megajoules_per_k = 0.01,
        .absorbed_shortwave_megajoules_per_step = 0.2,
        .absorbed_sky_longwave_megajoules_per_step = 0.1,
        .absorbed_lateral_longwave_megajoules_per_step = 0,
        .emitted_longwave_coefficient_megajoules_per_step_k4 = 1e-11,
        .ground_surface_temperature_k = 278,
        .radiation_fraction = 0.5,
        .air_vapor_volume_fraction = 0.01,
        .richardson_coefficient_k = 1,
        .minimum_richardson_number = -0.1,
        .maximum_richardson_number = 0.1,
        .minimum_surface_resistance_h_per_m = 0.1,
        .maximum_surface_resistance_h_per_m = 100,
        .isothermal_sensible_resistance_h_per_m = 1,
        .additional_latent_resistance_h_per_m = 1,
        .sensible_conductance_megajoules_per_m_k_step = 0.2,
        .latent_conductance_m2_per_step = 0.001,
        .maximum_water_removal_fraction_per_step = 1,
        .latent_heat_megajoules_per_m3 = 2450,
    };
}

fn testPolicy() SolverPolicy {
    return .{
        .legacy_substeps_nph = 4,
        .legacy_inner_iterations_mxn = 25,
        .absolute_temperature_tolerance_k = 1e-9,
        .relative_tolerance = 1e-10,
        .picard_relaxation = 0.5,
        .minimum_temperature_k = 200,
        .maximum_temperature_k = 350,
    };
}

test "standing-dead local hybrid exits before legacy product ceiling" {
    const result = try solve(testInputs(), testPolicy());
    try std.testing.expect(result.iterations < 100);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
    try std.testing.expect(@abs(result.residual_k) <= 1e-9 + 1e-10 * 280);
}

test "zero energy residual exits immediately without a subhour cycle" {
    var inputs = testInputs();
    inputs.absorbed_shortwave_megajoules_per_step = 0;
    inputs.absorbed_sky_longwave_megajoules_per_step = 0;
    inputs.emitted_longwave_coefficient_megajoules_per_step_k4 = 0;
    inputs.sensible_conductance_megajoules_per_m_k_step = 0;
    inputs.latent_conductance_m2_per_step = 0;
    const result = try solve(inputs, testPolicy());
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(f64, 280), result.surface_temperature_k);
}

test "runtime species solve fails atomically on later invalid state" {
    var invalid = testInputs();
    invalid.initial_surface_temperature_k = std.math.nan(f64);
    const inputs = [_]Inputs{ testInputs(), invalid };
    var scratch: [2]Result = undefined;
    var destination: [2]Result = undefined;
    destination[0] = try solve(testInputs(), testPolicy());
    destination[1] = destination[0];
    destination[0].surface_temperature_k = 241;
    destination[1].surface_temperature_k = 242;
    try std.testing.expectError(
        error.NonFiniteStandingDeadEnergyInput,
        solveRuntimeSpecies(&inputs, testPolicy(), &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 241), destination[0].surface_temperature_k);
    try std.testing.expectEqual(@as(f64, 242), destination[1].surface_temperature_k);
}
