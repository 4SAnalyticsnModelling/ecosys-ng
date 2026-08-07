const std = @import("std");
const numerics = @import("../../core/numerics.zig");

pub const DiffusionInputs = struct {
    microbial_radius_m: f64,
    water_film_thickness_m: f64,
    tortuosity: f64,
    aqueous_oxygen_diffusivity_m2_per_step: f64,
    microbial_count_per_g_c: f64,
    active_biomass_g_c: f64,
};

/// Ports NITRO DIFOX, including the spherical microbial surface factor 12.57.
pub fn uptakeConductance_m3_per_step(inputs: DiffusionInputs) !f64 {
    inline for (@typeInfo(DiffusionInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteOxygenDiffusionInput;
    if (inputs.microbial_radius_m <= 0 or inputs.water_film_thickness_m <= 0 or inputs.tortuosity < 0 or inputs.aqueous_oxygen_diffusivity_m2_per_step < 0 or inputs.microbial_count_per_g_c < 0 or inputs.active_biomass_g_c < 0) return error.InvalidOxygenDiffusionInput;
    const radial_factor_m = inputs.microbial_radius_m * (inputs.water_film_thickness_m + inputs.microbial_radius_m) / inputs.water_film_thickness_m;
    const conductance = inputs.tortuosity * inputs.aqueous_oxygen_diffusivity_m2_per_step * 12.57 * inputs.microbial_count_per_g_c * inputs.active_biomass_g_c * radial_factor_m;
    if (!std.math.isFinite(conductance)) return error.NonFiniteOxygenConductance;
    return conductance;
}

pub const Inputs = struct {
    allocated_gaseous_oxygen_g_o: f64,
    allocated_aqueous_oxygen_g_o: f64,
    allocated_gaseous_flux_g_o: f64,
    allocated_aqueous_flux_g_o: f64,
    water_volume_m3: f64,
    air_volume_m3: f64,
    population_allocation_fraction: f64,
    oxygen_solubility_water_to_air: f64,
    gas_exchange_rate_per_step: f64,
    uptake_conductance_m3_per_step: f64,
    oxygen_half_saturation_g_o_per_m3: f64,
    maximum_oxygen_uptake_g_o: f64,
    maximum_aqueous_oxygen_concentration_g_o_per_m3: f64,
};

pub const Options = struct {
    absolute_tolerance_g_o: f64,
    relative_tolerance: f64,
    derivative_floor: f64,
    picard_relaxation: f64,
    gas_max_iterations: u16,
};

pub const Result = struct {
    gaseous_oxygen_g_o: f64,
    aqueous_oxygen_g_o: f64,
    oxygen_uptake_g_o: f64,
    gas_to_water_exchange_g_o: f64,
    demand_satisfaction_fraction: f64,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    residual_g_o: f64,
};

const Context = struct {
    inputs: Inputs,
    gas_before_exchange_g_o: f64,
    aqueous_before_exchange_g_o: f64,
    total_before_uptake_g_o: f64,
};

/// Replaces NITRO's nested NPH/NPT O2 cycles with one bounded implicit
/// Newton-Raphson/Picard solve. `gas_max_iterations` is NPH*NPG from runtime
/// options, and convergence exits immediately.
pub fn solve(inputs: Inputs, options: Options) !Result {
    try validate(inputs, options);
    const gas_before = inputs.allocated_gaseous_oxygen_g_o + inputs.allocated_gaseous_flux_g_o;
    const aqueous_before = inputs.allocated_aqueous_oxygen_g_o + inputs.allocated_aqueous_flux_g_o;
    if (gas_before < 0 or aqueous_before < 0) return error.NegativeAvailableOxygen;
    const total_before = gas_before + aqueous_before;
    if (total_before == 0 or inputs.maximum_oxygen_uptake_g_o == 0) return .{ .gaseous_oxygen_g_o = gas_before, .aqueous_oxygen_g_o = aqueous_before, .oxygen_uptake_g_o = 0, .gas_to_water_exchange_g_o = 0, .demand_satisfaction_fraction = if (inputs.maximum_oxygen_uptake_g_o == 0) 1 else 0, .iterations = 0, .newton_raphson_steps = 0, .picard_steps = 0, .residual_g_o = 0 };
    // Once the requested relative mass tolerance underflows to zero, no
    // representable Newton/Picard correction can resolve this inventory at
    // the configured precision. Preserve the trace inventory exactly instead
    // of constructing a zero absolute tolerance or silently consuming it.
    const relative_inventory_tolerance_g_o = total_before * options.relative_tolerance;
    // Arithmetic on subnormal gram inventories cannot reliably distinguish a
    // residual correction from representation noise (the Arctic example
    // reaches ~4e-314 g O). Preserve that finite trace exactly; it is far
    // below any physical model resolution and must not consume NPH*NPG.
    if (total_before <= std.math.floatMin(f64) or relative_inventory_tolerance_g_o == 0) return .{ .gaseous_oxygen_g_o = gas_before, .aqueous_oxygen_g_o = aqueous_before, .oxygen_uptake_g_o = 0, .gas_to_water_exchange_g_o = 0, .demand_satisfaction_fraction = 0, .iterations = 0, .newton_raphson_steps = 0, .picard_steps = 0, .residual_g_o = 0 };
    const context: Context = .{ .inputs = inputs, .gas_before_exchange_g_o = gas_before, .aqueous_before_exchange_g_o = aqueous_before, .total_before_uptake_g_o = total_before };
    // Start from the phase-equilibrium partition rather than a dry aqueous
    // boundary. This seeds dissolution when all available O2 begins in gas
    // and avoids a Picard two-cycle at the physical lower bound.
    const water_capacity_m3 = inputs.water_volume_m3 * inputs.oxygen_solubility_water_to_air;
    const equilibrium_aqueous_g_o = total_before * water_capacity_m3 / (water_capacity_m3 + inputs.air_volume_m3);
    const initial = std.math.clamp(@max(aqueous_before, equilibrium_aqueous_g_o), 0, total_before);
    // A run-wide absolute tolerance can exceed a depleted population's
    // entire O2 allocation. Scale it down for small inventories so an
    // initially bounded aqueous mass cannot be accepted while its implied
    // gas mass is negative.
    const inventory_absolute_tolerance_g_o = @min(options.absolute_tolerance_g_o, relative_inventory_tolerance_g_o);
    const solved = try numerics.newtonPicard(&context, residual, derivative, picard, 0, total_before, initial, .{
        .absolute_tolerance = inventory_absolute_tolerance_g_o,
        .relative_tolerance = options.relative_tolerance,
        .derivative_floor = options.derivative_floor,
        .picard_relaxation = options.picard_relaxation,
        // The residual is an oxygen mass balance, so its scale is the
        // available inventory. Potential demand may be orders of magnitude
        // larger and must not loosen conservation convergence.
        .residual_scale = total_before,
        .max_iterations = options.gas_max_iterations,
    });
    const uptake = activeUptake(context, solved.root).value;
    const exchange = solved.root - aqueous_before + uptake;
    var final_gas = gas_before - exchange;
    const roundoff_tolerance = inventory_absolute_tolerance_g_o + options.relative_tolerance * total_before;
    if (final_gas < -roundoff_tolerance or solved.root < -roundoff_tolerance or uptake < -roundoff_tolerance) {
        std.log.err(
            "invalid oxygen solution: final_gas_g_o={e} aqueous_g_o={e} uptake_g_o={e} gas_before_g_o={e} aqueous_before_g_o={e} exchange_g_o={e} residual_g_o={e} tolerance_g_o={e} iterations={d}",
            .{ final_gas, solved.root, uptake, gas_before, aqueous_before, exchange, solved.residual, roundoff_tolerance, solved.iterations },
        );
        return error.InvalidOxygenSolution;
    }
    final_gas = @max(0, final_gas);
    const mass_error = final_gas + solved.root + uptake - total_before;
    if (@abs(mass_error) > roundoff_tolerance) return error.OxygenMassBalanceFailure;
    return .{
        .gaseous_oxygen_g_o = final_gas,
        .aqueous_oxygen_g_o = @max(0, solved.root),
        .oxygen_uptake_g_o = @max(0, uptake),
        .gas_to_water_exchange_g_o = exchange,
        .demand_satisfaction_fraction = std.math.clamp(uptake / inputs.maximum_oxygen_uptake_g_o, 0, 1),
        .iterations = solved.iterations,
        .newton_raphson_steps = solved.newton_raphson_steps,
        .picard_steps = solved.picard_steps,
        .residual_g_o = solved.residual,
    };
}

const UptakeValue = struct { value: f64, derivative: f64 };

fn activeUptake(context: Context, aqueous_mass_g_o: f64) UptakeValue {
    const inputs = context.inputs;
    const effective_water_m3 = inputs.water_volume_m3 * inputs.population_allocation_fraction;
    const unconstrained_concentration = @max(0, aqueous_mass_g_o / effective_water_m3);
    const concentration = @min(inputs.maximum_aqueous_oxygen_concentration_g_o_per_m3, unconstrained_concentration);
    const dc_da: f64 = if (unconstrained_concentration > 0 and unconstrained_concentration < inputs.maximum_aqueous_oxygen_concentration_g_o_per_m3) 1 / effective_water_m3 else 0;
    const x = inputs.uptake_conductance_m3_per_step * concentration;
    if (x <= 0 or aqueous_mass_g_o <= 0) return .{ .value = 0, .derivative = 0 };
    const demand = inputs.maximum_oxygen_uptake_g_o;
    const sum = demand + inputs.uptake_conductance_m3_per_step * inputs.oxygen_half_saturation_g_o_per_m3 + x;
    const discriminant = @max(0, sum * sum - 4 * x * demand);
    const root = @sqrt(discriminant);
    // Stable form of the source quadratic. The explicit sub-hour source kept
    // uptake below the current aqueous inventory through its tiny XNPG step;
    // the implicit replacement enforces that conservation bound directly.
    const unconstrained_value = if (sum + root > 0) 2 * x * demand / (sum + root) else 0;
    if (unconstrained_value >= aqueous_mass_g_o) return .{ .value = aqueous_mass_g_o, .derivative = 1 };
    const value = unconstrained_value;
    const dx_da = inputs.uptake_conductance_m3_per_step * dc_da;
    const derivative_value = if (root > std.math.floatEps(f64)) 0.5 * (dx_da - (2 * sum * dx_da - 4 * demand * dx_da) / (2 * root)) else 0;
    return .{ .value = value, .derivative = derivative_value };
}

fn exchangeAndDerivative(context: Context, aqueous_mass_g_o: f64, uptake: UptakeValue) struct { value: f64, derivative: f64 } {
    const inputs = context.inputs;
    const water_capacity_m3 = inputs.water_volume_m3 * inputs.oxygen_solubility_water_to_air;
    const combined_capacity_m3 = water_capacity_m3 + inputs.air_volume_m3;
    const gas_mass_g_o = context.total_before_uptake_g_o - aqueous_mass_g_o - uptake.value;
    const gas_derivative = -1 - uptake.derivative;
    const value = inputs.gas_exchange_rate_per_step * (@max(0, gas_mass_g_o) * water_capacity_m3 - aqueous_mass_g_o * inputs.air_volume_m3) / combined_capacity_m3;
    const derivative_value = inputs.gas_exchange_rate_per_step * ((if (gas_mass_g_o > 0) gas_derivative * water_capacity_m3 else 0) - inputs.air_volume_m3) / combined_capacity_m3;
    return .{ .value = value, .derivative = derivative_value };
}

fn residual(context: *const Context, aqueous_mass_g_o: f64) f64 {
    const uptake = activeUptake(context.*, aqueous_mass_g_o);
    const exchange = exchangeAndDerivative(context.*, aqueous_mass_g_o, uptake);
    return aqueous_mass_g_o - context.aqueous_before_exchange_g_o + uptake.value - exchange.value;
}

fn derivative(context: *const Context, aqueous_mass_g_o: f64) f64 {
    const uptake = activeUptake(context.*, aqueous_mass_g_o);
    const exchange = exchangeAndDerivative(context.*, aqueous_mass_g_o, uptake);
    return 1 + uptake.derivative - exchange.derivative;
}

fn picard(context: *const Context, aqueous_mass_g_o: f64) f64 {
    const uptake = activeUptake(context.*, aqueous_mass_g_o);
    const exchange = exchangeAndDerivative(context.*, aqueous_mass_g_o, uptake);
    return context.aqueous_before_exchange_g_o - uptake.value + exchange.value;
}

fn validate(inputs: Inputs, options: Options) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteOxygenSolverInput;
    inline for (@typeInfo(Options).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(options, field.name))) return error.NonFiniteOxygenSolverOption;
    if (inputs.allocated_gaseous_oxygen_g_o < 0 or inputs.allocated_aqueous_oxygen_g_o < 0 or inputs.water_volume_m3 <= 0 or inputs.air_volume_m3 < 0 or inputs.population_allocation_fraction <= 0 or inputs.population_allocation_fraction > 1 or inputs.oxygen_solubility_water_to_air < 0 or inputs.gas_exchange_rate_per_step < 0 or inputs.uptake_conductance_m3_per_step < 0 or inputs.oxygen_half_saturation_g_o_per_m3 <= 0 or inputs.maximum_oxygen_uptake_g_o < 0 or inputs.maximum_aqueous_oxygen_concentration_g_o_per_m3 <= 0) return error.InvalidOxygenSolverInput;
    if (inputs.water_volume_m3 * inputs.oxygen_solubility_water_to_air + inputs.air_volume_m3 <= 0 or options.absolute_tolerance_g_o <= 0 or options.relative_tolerance <= 0 or options.derivative_floor <= 0 or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or options.gas_max_iterations == 0) return error.InvalidOxygenSolverOption;
}

fn testInputs() Inputs {
    return .{ .allocated_gaseous_oxygen_g_o = 5, .allocated_aqueous_oxygen_g_o = 1, .allocated_gaseous_flux_g_o = 0, .allocated_aqueous_flux_g_o = 0, .water_volume_m3 = 2, .air_volume_m3 = 3, .population_allocation_fraction = 1, .oxygen_solubility_water_to_air = 0.03, .gas_exchange_rate_per_step = 0.5, .uptake_conductance_m3_per_step = 2, .oxygen_half_saturation_g_o_per_m3 = 0.1, .maximum_oxygen_uptake_g_o = 0.8, .maximum_aqueous_oxygen_concentration_g_o_per_m3 = 1 };
}

fn testOptions() Options {
    return .{ .absolute_tolerance_g_o = 1e-12, .relative_tolerance = 1e-10, .derivative_floor = 1e-14, .picard_relaxation = 0.5, .gas_max_iterations = 80 };
}

test "hybrid oxygen solve converges early and conserves mass" {
    const inputs = testInputs();
    const result = try solve(inputs, testOptions());
    try std.testing.expect(result.iterations < 80);
    try std.testing.expect(result.newton_raphson_steps > 0);
    try std.testing.expectApproxEqAbs(inputs.allocated_gaseous_oxygen_g_o + inputs.allocated_aqueous_oxygen_g_o, result.gaseous_oxygen_g_o + result.aqueous_oxygen_g_o + result.oxygen_uptake_g_o, 1e-9);
    try std.testing.expect(result.demand_satisfaction_fraction > 0 and result.demand_satisfaction_fraction <= 1);
}

test "zero oxygen returns without consuming iteration budget" {
    var inputs = testInputs();
    inputs.allocated_gaseous_oxygen_g_o = 0;
    inputs.allocated_aqueous_oxygen_g_o = 0;
    const result = try solve(inputs, testOptions());
    try std.testing.expectEqual(@as(u16, 0), result.iterations);
    try std.testing.expectEqual(@as(f64, 0), result.oxygen_uptake_g_o);
}

test "subnormal oxygen trace is preserved without nonlinear iterations" {
    var inputs = testInputs();
    inputs.allocated_gaseous_oxygen_g_o = 4.0e-314;
    inputs.allocated_aqueous_oxygen_g_o = 4.0e-315;
    const result = try solve(inputs, testOptions());
    try std.testing.expectEqual(@as(u16, 0), result.iterations);
    try std.testing.expectEqual(inputs.allocated_gaseous_oxygen_g_o, result.gaseous_oxygen_g_o);
    try std.testing.expectEqual(inputs.allocated_aqueous_oxygen_g_o, result.aqueous_oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 0), result.oxygen_uptake_g_o);
}

test "oxygen diffusion conductance preserves NITRO spherical factor" {
    const conductance = try uptakeConductance_m3_per_step(.{ .microbial_radius_m = 1e-6, .water_film_thickness_m = 2e-6, .tortuosity = 0.5, .aqueous_oxygen_diffusivity_m2_per_step = 1e-4, .microbial_count_per_g_c = 1e12, .active_biomass_g_c = 0.01 });
    try std.testing.expectApproxEqRel(@as(f64, 0.5 * 1e-4 * 12.57 * 1e12 * 0.01 * 1.5e-6), conductance, 1e-14);
}

test "oxygen solver falls back to Picard when Newton derivative is rejected" {
    var options = testOptions();
    options.derivative_floor = 1e9;
    options.relative_tolerance = 1e-8;
    const result = try solve(testInputs(), options);
    try std.testing.expect(result.picard_steps > 0);
    try std.testing.expectEqual(@as(u16, 0), result.newton_raphson_steps);
}

test "depleted oxygen inventory uses a mass-conserving scaled tolerance" {
    var inputs = testInputs();
    inputs.allocated_gaseous_oxygen_g_o = 7.0e-13;
    inputs.allocated_aqueous_oxygen_g_o = 1.19e-11;
    inputs.maximum_oxygen_uptake_g_o = 1.0e-9;
    var options = testOptions();
    options.absolute_tolerance_g_o = 1.0e-11;
    options.relative_tolerance = 1.0e-8;

    const result = try solve(inputs, options);
    const before = inputs.allocated_gaseous_oxygen_g_o + inputs.allocated_aqueous_oxygen_g_o;
    try std.testing.expect(result.gaseous_oxygen_g_o >= 0);
    try std.testing.expect(result.aqueous_oxygen_g_o >= 0);
    try std.testing.expectApproxEqAbs(before, result.gaseous_oxygen_g_o + result.aqueous_oxygen_g_o + result.oxygen_uptake_g_o, before * 3.0e-8);
}

test "subnormal oxygen below relative resolution is preserved exactly" {
    var inputs = testInputs();
    inputs.allocated_gaseous_oxygen_g_o = std.math.floatTrueMin(f64);
    inputs.allocated_aqueous_oxygen_g_o = 0;
    const result = try solve(inputs, testOptions());
    try std.testing.expectEqual(inputs.allocated_gaseous_oxygen_g_o, result.gaseous_oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 0), result.aqueous_oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 0), result.oxygen_uptake_g_o);
    try std.testing.expectEqual(@as(u16, 0), result.iterations);
}
