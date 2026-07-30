const std = @import("std");
const numerics = @import("numerics.zig");

pub const Inputs = struct {
    gaseous_methane_g_c: f64,
    aqueous_methane_g_c: f64,
    gaseous_methane_flux_g_c: f64,
    aqueous_methane_flux_g_c: f64,
    methanogenesis_g_c: f64,
    water_volume_m3: f64,
    air_volume_m3: f64,
    methane_solubility_water_to_air: f64,
    gas_exchange_rate_per_step: f64,
    gas_exchange_enabled: bool,
    minimum_gaseous_methane_g_c: f64,
    methane_half_saturation_g_c_per_m3: f64,
    maximum_methane_oxidation_g_c: f64,
    biomass_conversion_efficiency_g_c_per_g_c: f64,
    growth_respiration_g_c_per_g_c: f64,
};

pub const Options = struct {
    absolute_tolerance_g_c: f64,
    relative_tolerance: f64,
    derivative_floor: f64,
    picard_relaxation: f64,
    gas_max_iterations: u16,
};

pub const Result = struct {
    gaseous_methane_g_c: f64,
    aqueous_methane_g_c: f64,
    methane_oxidation_g_c: f64,
    growth_respiration_g_c: f64,
    gas_to_water_exchange_g_c: f64,
    oxygen_demand_g_o: f64,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    residual_g_c: f64,
};

const Context = struct {
    inputs: Inputs,
    gas_before_exchange_g_c: f64,
    aqueous_before_exchange_g_c: f64,
    total_before_consumption_g_c: f64,
    respiration_per_oxidation: f64,
};

const Consumption = struct { total_g_c: f64, derivative: f64 };

/// Replaces NITRO's NPH×NPT methane dissolution/oxidation cycling with a
/// single bounded implicit Newton–Raphson/Picard solve. The caller supplies
/// the runtime NPH×NPG gas ceiling and convergence exits immediately.
pub fn solve(inputs: Inputs, options: Options) !Result {
    try validate(inputs, options);
    const gas_before = inputs.gaseous_methane_g_c + inputs.gaseous_methane_flux_g_c;
    const aqueous_before = inputs.aqueous_methane_g_c + inputs.aqueous_methane_flux_g_c + inputs.methanogenesis_g_c;
    if (gas_before < 0 or aqueous_before < 0) return error.NegativeAvailableMethane;
    const total_before = gas_before + aqueous_before;
    if (total_before == 0 or inputs.maximum_methane_oxidation_g_c == 0) return .{
        .gaseous_methane_g_c = gas_before,
        .aqueous_methane_g_c = aqueous_before,
        .methane_oxidation_g_c = 0,
        .growth_respiration_g_c = 0,
        .gas_to_water_exchange_g_c = 0,
        .oxygen_demand_g_o = 0,
        .iterations = 0,
        .newton_raphson_steps = 0,
        .picard_steps = 0,
        .residual_g_c = 0,
    };
    const respiration_per_oxidation = inputs.biomass_conversion_efficiency_g_c_per_g_c * inputs.growth_respiration_g_c_per_g_c;
    const context: Context = .{ .inputs = inputs, .gas_before_exchange_g_c = gas_before, .aqueous_before_exchange_g_c = aqueous_before, .total_before_consumption_g_c = total_before, .respiration_per_oxidation = respiration_per_oxidation };
    const solved = try numerics.newtonPicard(&context, residual, derivative, picard, 0, total_before, @min(aqueous_before, total_before), .{
        .absolute_tolerance = options.absolute_tolerance_g_c,
        .relative_tolerance = options.relative_tolerance,
        .derivative_floor = options.derivative_floor,
        .picard_relaxation = options.picard_relaxation,
        .residual_scale = @max(total_before, inputs.maximum_methane_oxidation_g_c),
        .max_iterations = options.gas_max_iterations,
    });
    const consumed = consumption(context, solved.root).total_g_c;
    const oxidation = consumed / (1 + respiration_per_oxidation);
    const respiration = consumed - oxidation;
    const gas_to_water_exchange_g_c = solved.root - aqueous_before + consumed;
    var final_gas = gas_before - gas_to_water_exchange_g_c;
    const tolerance = options.absolute_tolerance_g_c + options.relative_tolerance * total_before;
    if (final_gas < -tolerance or solved.root < -tolerance or consumed < -tolerance) return error.InvalidMethaneSolution;
    final_gas = @max(0, final_gas);
    if (@abs(final_gas + solved.root + consumed - total_before) > tolerance) return error.MethaneMassBalanceFailure;
    return .{
        .gaseous_methane_g_c = final_gas,
        .aqueous_methane_g_c = @max(0, solved.root),
        .methane_oxidation_g_c = oxidation,
        .growth_respiration_g_c = respiration,
        .gas_to_water_exchange_g_c = gas_to_water_exchange_g_c,
        .oxygen_demand_g_o = 5.333 * oxidation + 2.667 * respiration,
        .iterations = solved.iterations,
        .newton_raphson_steps = solved.newton_raphson_steps,
        .picard_steps = solved.picard_steps,
        .residual_g_c = solved.residual,
    };
}

fn consumption(context: Context, aqueous_g_c: f64) Consumption {
    if (aqueous_g_c <= 0) return .{ .total_g_c = 0, .derivative = 0 };
    const concentration = aqueous_g_c / context.inputs.water_volume_m3;
    const half_saturation = context.inputs.methane_half_saturation_g_c_per_m3;
    const total_kinetic_capacity = (1 + context.respiration_per_oxidation) * context.inputs.maximum_methane_oxidation_g_c;
    const kinetic = total_kinetic_capacity * concentration / (concentration + half_saturation);
    if (aqueous_g_c <= kinetic) return .{ .total_g_c = aqueous_g_c, .derivative = 1 };
    const kinetic_derivative = total_kinetic_capacity * half_saturation / (context.inputs.water_volume_m3 * std.math.pow(f64, concentration + half_saturation, 2));
    return .{ .total_g_c = kinetic, .derivative = kinetic_derivative };
}

fn exchange(context: Context, aqueous_g_c: f64, consumed: Consumption) Consumption {
    if (!context.inputs.gas_exchange_enabled or context.inputs.gas_exchange_rate_per_step == 0) return .{ .total_g_c = 0, .derivative = 0 };
    const water_capacity_m3 = context.inputs.water_volume_m3 * context.inputs.methane_solubility_water_to_air;
    const total_capacity_m3 = water_capacity_m3 + context.inputs.air_volume_m3;
    const unconstrained_gas = context.total_before_consumption_g_c - aqueous_g_c - consumed.total_g_c;
    const gas = @max(context.inputs.minimum_gaseous_methane_g_c, unconstrained_gas);
    const gas_derivative: f64 = if (unconstrained_gas > context.inputs.minimum_gaseous_methane_g_c) -1 - consumed.derivative else 0;
    return .{
        .total_g_c = context.inputs.gas_exchange_rate_per_step * (gas * water_capacity_m3 - aqueous_g_c * context.inputs.air_volume_m3) / total_capacity_m3,
        .derivative = context.inputs.gas_exchange_rate_per_step * (gas_derivative * water_capacity_m3 - context.inputs.air_volume_m3) / total_capacity_m3,
    };
}

fn residual(context: *const Context, aqueous_g_c: f64) f64 {
    const consumed = consumption(context.*, aqueous_g_c);
    return aqueous_g_c - context.aqueous_before_exchange_g_c + consumed.total_g_c - exchange(context.*, aqueous_g_c, consumed).total_g_c;
}

fn derivative(context: *const Context, aqueous_g_c: f64) f64 {
    const consumed = consumption(context.*, aqueous_g_c);
    return 1 + consumed.derivative - exchange(context.*, aqueous_g_c, consumed).derivative;
}

fn picard(context: *const Context, aqueous_g_c: f64) f64 {
    const consumed = consumption(context.*, aqueous_g_c);
    return context.aqueous_before_exchange_g_c - consumed.total_g_c + exchange(context.*, aqueous_g_c, consumed).total_g_c;
}

fn validate(inputs: Inputs, options: Options) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteMethaneOxidationInput;
    inline for (@typeInfo(Options).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(options, field.name))) return error.NonFiniteMethaneOxidationOption;
    if (inputs.gaseous_methane_g_c < 0 or inputs.aqueous_methane_g_c < 0 or inputs.water_volume_m3 <= 0 or inputs.air_volume_m3 < 0 or inputs.methane_solubility_water_to_air < 0 or inputs.gas_exchange_rate_per_step < 0 or inputs.minimum_gaseous_methane_g_c < 0 or inputs.methane_half_saturation_g_c_per_m3 <= 0 or inputs.maximum_methane_oxidation_g_c < 0 or inputs.biomass_conversion_efficiency_g_c_per_g_c < 0 or inputs.growth_respiration_g_c_per_g_c < 0) return error.InvalidMethaneOxidationInput;
    if (inputs.gas_exchange_enabled and inputs.water_volume_m3 * inputs.methane_solubility_water_to_air + inputs.air_volume_m3 <= 0) return error.InvalidMethaneExchangeCapacity;
    if (options.absolute_tolerance_g_c <= 0 or options.relative_tolerance <= 0 or options.derivative_floor <= 0 or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or options.gas_max_iterations == 0) return error.InvalidMethaneOxidationOption;
}

fn testInputs() Inputs {
    return .{ .gaseous_methane_g_c = 4, .aqueous_methane_g_c = 1, .gaseous_methane_flux_g_c = 0.1, .aqueous_methane_flux_g_c = 0.1, .methanogenesis_g_c = 0.2, .water_volume_m3 = 2, .air_volume_m3 = 3, .methane_solubility_water_to_air = 0.03, .gas_exchange_rate_per_step = 0.5, .gas_exchange_enabled = true, .minimum_gaseous_methane_g_c = 1e-12, .methane_half_saturation_g_c_per_m3 = 0.2, .maximum_methane_oxidation_g_c = 0.3, .biomass_conversion_efficiency_g_c_per_g_c = 0.4, .growth_respiration_g_c_per_g_c = 0.5 };
}

fn testOptions() Options {
    return .{ .absolute_tolerance_g_c = 1e-12, .relative_tolerance = 1e-10, .derivative_floor = 1e-14, .picard_relaxation = 0.5, .gas_max_iterations = 80 };
}

test "hybrid methane oxidation converges early and conserves carbon" {
    const inputs = testInputs();
    const result = try solve(inputs, testOptions());
    const before = inputs.gaseous_methane_g_c + inputs.aqueous_methane_g_c + inputs.gaseous_methane_flux_g_c + inputs.aqueous_methane_flux_g_c + inputs.methanogenesis_g_c;
    try std.testing.expect(result.iterations < 80);
    try std.testing.expect(result.newton_raphson_steps > 0);
    try std.testing.expectApproxEqAbs(before, result.gaseous_methane_g_c + result.aqueous_methane_g_c + result.methane_oxidation_g_c + result.growth_respiration_g_c, 1e-9);
    try std.testing.expectApproxEqAbs(5.333 * result.methane_oxidation_g_c + 2.667 * result.growth_respiration_g_c, result.oxygen_demand_g_o, 1e-14);
}

test "zero methane exits without iteration" {
    var inputs = testInputs();
    inputs.gaseous_methane_g_c = 0;
    inputs.aqueous_methane_g_c = 0;
    inputs.gaseous_methane_flux_g_c = 0;
    inputs.aqueous_methane_flux_g_c = 0;
    inputs.methanogenesis_g_c = 0;
    const result = try solve(inputs, testOptions());
    try std.testing.expectEqual(@as(u16, 0), result.iterations);
}
