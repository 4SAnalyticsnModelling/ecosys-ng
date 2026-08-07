const std = @import("std");
const numerics = @import("../core/numerics.zig");

pub const Parameters = struct {
    minimum_richardson_number: f64,
    maximum_richardson_number: f64,
    richardson_resistance_multiplier: f64,
    minimum_aerodynamic_resistance_h_per_m: f64,
    maximum_aerodynamic_resistance_h_per_m: f64,
    volumetric_air_heat_capacity_megajoules_per_m3_k: f64,
    minimum_air_column_height_m: f64,
    sensible_heat_conductivity_megajoules_per_m_h_k: f64,
    liquid_water_latent_heat_megajoules_per_m3: f64,
    saturation_vapor_prefactor_k: f64,
    saturation_relative_humidity: f64,
    saturation_temperature_k: f64,
    saturation_reference_inverse_temperature_per_k: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    temperature_k: []f64,
    vapor_volume_fraction: []f64,
    heat_capacity_megajoules_per_k: []f64,
    air_volume_m3: []f64,
    iteration_count: []u16,

    pub fn init(allocator: std.mem.Allocator, initial_temperature_k: []const f64, initial_vapor_volume_fraction: []const f64, cell_area_m2: []const f64, reference_height_m: f64, parameters: Parameters) !State {
        const count = initial_temperature_k.len;
        if (count == 0 or initial_vapor_volume_fraction.len != count or cell_area_m2.len != count or !std.math.isFinite(reference_height_m) or reference_height_m <= 0) return error.InvalidGroundAirInitialization;
        const temperature = try allocator.dupe(f64, initial_temperature_k);
        errdefer allocator.free(temperature);
        const vapor = try allocator.dupe(f64, initial_vapor_volume_fraction);
        errdefer allocator.free(vapor);
        const capacity = try allocator.alloc(f64, count);
        errdefer allocator.free(capacity);
        const volume = try allocator.alloc(f64, count);
        errdefer allocator.free(volume);
        const iterations = try allocator.alloc(u16, count);
        errdefer allocator.free(iterations);
        @memset(iterations, 0);
        for (0..count) |cell| {
            if (!std.math.isFinite(temperature[cell]) or temperature[cell] <= 0 or !std.math.isFinite(vapor[cell]) or vapor[cell] < 0 or !std.math.isFinite(cell_area_m2[cell]) or cell_area_m2[cell] <= 0) return error.InvalidGroundAirInitialization;
            const height = @max(parameters.minimum_air_column_height_m, reference_height_m);
            volume[cell] = height * cell_area_m2[cell];
            capacity[cell] = volume[cell] * parameters.volumetric_air_heat_capacity_megajoules_per_m3_k;
        }
        return .{ .allocator = allocator, .cell_count = count, .temperature_k = temperature, .vapor_volume_fraction = vapor, .heat_capacity_megajoules_per_k = capacity, .air_volume_m3 = volume, .iteration_count = iterations };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.iteration_count);
        self.allocator.free(self.air_volume_m3);
        self.allocator.free(self.heat_capacity_megajoules_per_k);
        self.allocator.free(self.vapor_volume_fraction);
        self.allocator.free(self.temperature_k);
        self.* = undefined;
    }

    pub fn refreshGeometry(self: *State, cell_area_m2: []const f64, reference_height_m: []const f64, parameters: Parameters) !void {
        if (cell_area_m2.len != self.cell_count or reference_height_m.len != self.cell_count) return error.GroundAirDimensionMismatch;
        for (0..self.cell_count) |cell| {
            if (!std.math.isFinite(cell_area_m2[cell]) or cell_area_m2[cell] <= 0 or !std.math.isFinite(reference_height_m[cell]) or reference_height_m[cell] <= 0) return error.InvalidGroundAirGeometry;
            self.air_volume_m3[cell] = @max(parameters.minimum_air_column_height_m, reference_height_m[cell]) * cell_area_m2[cell];
            self.heat_capacity_megajoules_per_k[cell] = self.air_volume_m3[cell] * parameters.volumetric_air_heat_capacity_megajoules_per_m3_k;
        }
    }
};

pub const Inputs = struct {
    atmospheric_temperature_k: []const f64,
    atmospheric_vapor_volume_fraction: []const f64,
    cell_area_m2: []const f64,
    bulk_richardson_coefficient_k: []const f64,
    neutral_atmospheric_resistance_h_per_m: []const f64,
    canopy_resistance_h_per_m: []const f64,
    non_atmospheric_sensible_heat_megajoules_per_h: []const f64,
    non_atmospheric_vapor_flux_m3_per_h: []const f64,
    non_atmospheric_sensible_conductance_megajoules_per_h_k: []const f64,
    non_atmospheric_sensible_source_temperature_k: []const f64,
    non_atmospheric_vapor_conductance_m3_per_h: []const f64,
    non_atmospheric_vapor_source_fraction: []const f64,
};

pub fn vaporPressureKpa(vapor_volume_fraction: f64, temperature_k: f64, parameters: Parameters) !f64 {
    if (!std.math.isFinite(vapor_volume_fraction) or vapor_volume_fraction < 0 or !std.math.isFinite(temperature_k) or temperature_k <= 0) return error.InvalidGroundAirVaporState;
    return vapor_volume_fraction * temperature_k / parameters.saturation_vapor_prefactor_k;
}

pub fn vaporVolumeFraction(vapor_pressure_kpa: f64, temperature_k: f64, parameters: Parameters) !f64 {
    if (!std.math.isFinite(vapor_pressure_kpa) or vapor_pressure_kpa < 0 or !std.math.isFinite(temperature_k) or temperature_k <= 0) return error.InvalidAtmosphericVaporState;
    return vapor_pressure_kpa * parameters.saturation_vapor_prefactor_k / temperature_k;
}

pub fn deriveSurfaceSources(cell_area_m2: []const f64, surface_sensible_heat_flux_megajoules_per_m2_h: []const f64, surface_latent_heat_flux_megajoules_per_m2_h: []const f64, parameters: Parameters, sensible_heat_source_megajoules_per_h: []f64, vapor_source_m3_per_h: []f64) !void {
    const count = cell_area_m2.len;
    inline for (.{ surface_sensible_heat_flux_megajoules_per_m2_h.len, surface_latent_heat_flux_megajoules_per_m2_h.len, sensible_heat_source_megajoules_per_h.len, vapor_source_m3_per_h.len }) |length| if (length != count) return error.GroundAirDimensionMismatch;
    try validateParameters(parameters);
    for (0..count) |cell| {
        inline for (.{ cell_area_m2[cell], surface_sensible_heat_flux_megajoules_per_m2_h[cell], surface_latent_heat_flux_megajoules_per_m2_h[cell] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteGroundAirSurfaceFlux;
        if (cell_area_m2[cell] <= 0) return error.InvalidGroundAirSurfaceArea;
        sensible_heat_source_megajoules_per_h[cell] = -surface_sensible_heat_flux_megajoules_per_m2_h[cell] * cell_area_m2[cell];
        vapor_source_m3_per_h[cell] = -surface_latent_heat_flux_megajoules_per_m2_h[cell] * cell_area_m2[cell] / parameters.liquid_water_latent_heat_megajoules_per_m3;
    }
}

const TemperatureContext = struct {
    old_temperature_k: f64,
    atmospheric_temperature_k: f64,
    heat_capacity_megajoules_per_k: f64,
    cell_area_m2: f64,
    bulk_richardson_coefficient_k: f64,
    neutral_resistance_h_per_m: f64,
    canopy_resistance_h_per_m: f64,
    source_heat_megajoules_per_h: f64,
    source_conductance_megajoules_per_h_k: f64,
    source_temperature_k: f64,
    parameters: Parameters,
};

pub fn solve(state: *State, inputs: Inputs, parameters: Parameters, solver_options: numerics.SolverOptions) !void {
    inline for (.{ inputs.atmospheric_temperature_k.len, inputs.atmospheric_vapor_volume_fraction.len, inputs.cell_area_m2.len, inputs.bulk_richardson_coefficient_k.len, inputs.neutral_atmospheric_resistance_h_per_m.len, inputs.canopy_resistance_h_per_m.len, inputs.non_atmospheric_sensible_heat_megajoules_per_h.len, inputs.non_atmospheric_vapor_flux_m3_per_h.len, inputs.non_atmospheric_sensible_conductance_megajoules_per_h_k.len, inputs.non_atmospheric_sensible_source_temperature_k.len, inputs.non_atmospheric_vapor_conductance_m3_per_h.len, inputs.non_atmospheric_vapor_source_fraction.len }) |length| if (length != state.cell_count) return error.GroundAirDimensionMismatch;
    try validateParameters(parameters);
    const candidate_temperature = try state.allocator.dupe(f64, state.temperature_k);
    defer state.allocator.free(candidate_temperature);
    const candidate_vapor = try state.allocator.dupe(f64, state.vapor_volume_fraction);
    defer state.allocator.free(candidate_vapor);
    const candidate_iterations = try state.allocator.dupe(u16, state.iteration_count);
    defer state.allocator.free(candidate_iterations);
    for (0..state.cell_count) |cell| {
        inline for (.{ state.temperature_k[cell], state.vapor_volume_fraction[cell], state.heat_capacity_megajoules_per_k[cell], state.air_volume_m3[cell], inputs.atmospheric_temperature_k[cell], inputs.atmospheric_vapor_volume_fraction[cell], inputs.cell_area_m2[cell], inputs.bulk_richardson_coefficient_k[cell], inputs.neutral_atmospheric_resistance_h_per_m[cell], inputs.canopy_resistance_h_per_m[cell], inputs.non_atmospheric_sensible_heat_megajoules_per_h[cell], inputs.non_atmospheric_vapor_flux_m3_per_h[cell], inputs.non_atmospheric_sensible_conductance_megajoules_per_h_k[cell], inputs.non_atmospheric_sensible_source_temperature_k[cell], inputs.non_atmospheric_vapor_conductance_m3_per_h[cell], inputs.non_atmospheric_vapor_source_fraction[cell] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteGroundAirInput;
        if (state.temperature_k[cell] <= 0 or state.vapor_volume_fraction[cell] < 0 or state.heat_capacity_megajoules_per_k[cell] <= 0 or state.air_volume_m3[cell] <= 0 or inputs.atmospheric_temperature_k[cell] <= 0 or inputs.atmospheric_vapor_volume_fraction[cell] < 0 or inputs.cell_area_m2[cell] <= 0 or inputs.neutral_atmospheric_resistance_h_per_m[cell] < 0 or inputs.canopy_resistance_h_per_m[cell] < 0 or inputs.non_atmospheric_sensible_conductance_megajoules_per_h_k[cell] < 0 or inputs.non_atmospheric_sensible_source_temperature_k[cell] <= 0 or inputs.non_atmospheric_vapor_conductance_m3_per_h[cell] < 0 or inputs.non_atmospheric_vapor_source_fraction[cell] < 0) return error.InvalidGroundAirInput;
        const context: TemperatureContext = .{ .old_temperature_k = state.temperature_k[cell], .atmospheric_temperature_k = inputs.atmospheric_temperature_k[cell], .heat_capacity_megajoules_per_k = state.heat_capacity_megajoules_per_k[cell], .cell_area_m2 = inputs.cell_area_m2[cell], .bulk_richardson_coefficient_k = inputs.bulk_richardson_coefficient_k[cell], .neutral_resistance_h_per_m = inputs.neutral_atmospheric_resistance_h_per_m[cell], .canopy_resistance_h_per_m = inputs.canopy_resistance_h_per_m[cell], .source_heat_megajoules_per_h = inputs.non_atmospheric_sensible_heat_megajoules_per_h[cell], .source_conductance_megajoules_per_h_k = inputs.non_atmospheric_sensible_conductance_megajoules_per_h_k[cell], .source_temperature_k = inputs.non_atmospheric_sensible_source_temperature_k[cell], .parameters = parameters };
        var options = solver_options;
        options.residual_scale = @max(1, @abs(context.source_heat_megajoules_per_h) + context.heat_capacity_megajoules_per_k);
        const solved = try numerics.newtonPicard(context, temperatureResidual, temperatureDerivative, temperaturePicard, 173.15, 373.15, context.old_temperature_k, options);
        candidate_temperature[cell] = solved.root;
        candidate_iterations[cell] = solved.iterations;
        const resistance = try atmosphereResistance(context, solved.root);
        const exchange_volume_m3_per_h = inputs.cell_area_m2[cell] / resistance;
        const old_vapor_m3 = state.vapor_volume_fraction[cell] * state.air_volume_m3[cell];
        const source_vapor_exchange_m3_per_h = inputs.non_atmospheric_vapor_conductance_m3_per_h[cell];
        const unconstrained_vapor_m3 = (old_vapor_m3 + inputs.non_atmospheric_vapor_flux_m3_per_h[cell] + exchange_volume_m3_per_h * inputs.atmospheric_vapor_volume_fraction[cell] + source_vapor_exchange_m3_per_h * inputs.non_atmospheric_vapor_source_fraction[cell]) / (1 + (exchange_volume_m3_per_h + source_vapor_exchange_m3_per_h) / state.air_volume_m3[cell]);
        const saturation_fraction = parameters.saturation_vapor_prefactor_k / solved.root * parameters.saturation_relative_humidity * @exp(parameters.saturation_temperature_k * (parameters.saturation_reference_inverse_temperature_per_k - 1 / solved.root));
        candidate_vapor[cell] = std.math.clamp(unconstrained_vapor_m3 / state.air_volume_m3[cell], 0, saturation_fraction);
        if (!std.math.isFinite(candidate_vapor[cell])) return error.NonFiniteGroundAirVaporResult;
    }
    @memcpy(state.temperature_k, candidate_temperature);
    @memcpy(state.vapor_volume_fraction, candidate_vapor);
    @memcpy(state.iteration_count, candidate_iterations);
}

fn atmosphereResistance(context: TemperatureContext, temperature_k: f64) !f64 {
    const p = context.parameters;
    const richardson = std.math.clamp(context.bulk_richardson_coefficient_k / context.atmospheric_temperature_k * (context.atmospheric_temperature_k - temperature_k), p.minimum_richardson_number, p.maximum_richardson_number);
    const stability = 1 - p.richardson_resistance_multiplier * richardson;
    if (stability <= 0) return error.InvalidGroundAirStability;
    return @min(p.maximum_aerodynamic_resistance_h_per_m, @max(p.minimum_aerodynamic_resistance_h_per_m, context.neutral_resistance_h_per_m + context.canopy_resistance_h_per_m) / stability);
}

fn temperatureResidual(context: TemperatureContext, temperature_k: f64) f64 {
    const resistance = atmosphereResistance(context, temperature_k) catch return std.math.nan(f64);
    const atmospheric_heat = context.parameters.sensible_heat_conductivity_megajoules_per_m_h_k * context.cell_area_m2 / resistance * (context.atmospheric_temperature_k - temperature_k);
    const source_heat = context.source_heat_megajoules_per_h + context.source_conductance_megajoules_per_h_k * (context.source_temperature_k - temperature_k);
    return source_heat + atmospheric_heat + context.heat_capacity_megajoules_per_k * (context.old_temperature_k - temperature_k);
}

fn temperatureDerivative(context: TemperatureContext, temperature_k: f64) f64 {
    const step = @max(1e-6, @abs(temperature_k) * 1e-6);
    return (temperatureResidual(context, temperature_k + step) - temperatureResidual(context, temperature_k - step)) / (2 * step);
}

fn temperaturePicard(context: TemperatureContext, temperature_k: f64) f64 {
    const resistance = atmosphereResistance(context, temperature_k) catch return context.old_temperature_k;
    const conductance = context.parameters.sensible_heat_conductivity_megajoules_per_m_h_k * context.cell_area_m2 / resistance;
    return (context.source_heat_megajoules_per_h + context.source_conductance_megajoules_per_h_k * context.source_temperature_k + conductance * context.atmospheric_temperature_k + context.heat_capacity_megajoules_per_k * context.old_temperature_k) / (context.source_conductance_megajoules_per_h_k + conductance + context.heat_capacity_megajoules_per_k);
}

fn validateParameters(p: Parameters) !void {
    inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(p, field.name))) return error.NonFiniteGroundAirParameter;
    if (p.maximum_richardson_number < p.minimum_richardson_number or p.richardson_resistance_multiplier <= 0 or p.minimum_aerodynamic_resistance_h_per_m <= 0 or p.maximum_aerodynamic_resistance_h_per_m < p.minimum_aerodynamic_resistance_h_per_m or p.volumetric_air_heat_capacity_megajoules_per_m3_k <= 0 or p.minimum_air_column_height_m <= 0 or p.sensible_heat_conductivity_megajoules_per_m_h_k <= 0 or p.liquid_water_latent_heat_megajoules_per_m3 <= 0 or p.saturation_vapor_prefactor_k <= 0 or p.saturation_relative_humidity < 0 or p.saturation_relative_humidity > 1 or p.saturation_temperature_k <= 0 or p.saturation_reference_inverse_temperature_per_k <= 0) return error.InvalidGroundAirParameter;
}

test "ground air converges without repeating the full model" {
    const parameters: Parameters = .{ .minimum_richardson_number = -0.1, .maximum_richardson_number = 0.05, .richardson_resistance_multiplier = 10, .minimum_aerodynamic_resistance_h_per_m = 0.00139, .maximum_aerodynamic_resistance_h_per_m = 0.0139, .volumetric_air_heat_capacity_megajoules_per_m3_k = 1.25e-3, .minimum_air_column_height_m = 5, .sensible_heat_conductivity_megajoules_per_m_h_k = 1.2e-3, .liquid_water_latent_heat_megajoules_per_m3 = 2465, .saturation_vapor_prefactor_k = 2.173e-3, .saturation_relative_humidity = 0.61, .saturation_temperature_k = 5360, .saturation_reference_inverse_temperature_per_k = 3.661e-3 };
    var state = try State.init(std.testing.allocator, &.{280}, &.{0.005}, &.{10}, 5, parameters);
    defer state.deinit();
    try solve(&state, .{ .atmospheric_temperature_k = &.{285}, .atmospheric_vapor_volume_fraction = &.{0.004}, .cell_area_m2 = &.{10}, .bulk_richardson_coefficient_k = &.{20}, .neutral_atmospheric_resistance_h_per_m = &.{0.005}, .canopy_resistance_h_per_m = &.{0.002}, .non_atmospheric_sensible_heat_megajoules_per_h = &.{0.1}, .non_atmospheric_vapor_flux_m3_per_h = &.{0.001}, .non_atmospheric_sensible_conductance_megajoules_per_h_k = &.{0}, .non_atmospheric_sensible_source_temperature_k = &.{280}, .non_atmospheric_vapor_conductance_m3_per_h = &.{0}, .non_atmospheric_vapor_source_fraction = &.{0} }, parameters, .{ .absolute_tolerance = 1e-10, .relative_tolerance = 1e-8, .max_iterations = 20, .picard_relaxation = 0.5 });
    try std.testing.expect(state.temperature_k[0] > 280);
    try std.testing.expect(state.iteration_count[0] < 20);
    try std.testing.expect(state.vapor_volume_fraction[0] >= 0);
}

test "surface fluxes become equal and opposite ground air sources" {
    const parameters: Parameters = .{ .minimum_richardson_number = -0.1, .maximum_richardson_number = 0.05, .richardson_resistance_multiplier = 10, .minimum_aerodynamic_resistance_h_per_m = 0.00139, .maximum_aerodynamic_resistance_h_per_m = 0.0139, .volumetric_air_heat_capacity_megajoules_per_m3_k = 1.25e-3, .minimum_air_column_height_m = 5, .sensible_heat_conductivity_megajoules_per_m_h_k = 1.2e-3, .liquid_water_latent_heat_megajoules_per_m3 = 2465, .saturation_vapor_prefactor_k = 2.173e-3, .saturation_relative_humidity = 0.61, .saturation_temperature_k = 5360, .saturation_reference_inverse_temperature_per_k = 3.661e-3 };
    var sensible: [1]f64 = undefined;
    var vapor: [1]f64 = undefined;
    try deriveSurfaceSources(&.{10}, &.{0.2}, &.{-0.493}, parameters, &sensible, &vapor);
    try std.testing.expectApproxEqAbs(-2.0, sensible[0], 1e-15);
    try std.testing.expectApproxEqAbs(0.002, vapor[0], 1e-15);
}
