const std = @import("std");

pub const Parameters = struct {
    saturation_vapor_prefactor_k: f64 = 2.173e-3,
    saturation_relative_humidity: f64 = 0.61,
    saturation_temperature_coefficient_k: f64 = 5360.0,
    saturation_reference_inverse_temperature_per_k: f64 = 3.661e-3,
};

pub const SolverOptions = struct {
    max_iterations: u16 = 100,
    absolute_temperature_tolerance_k: f64 = 1.0e-8,
    absolute_vapor_fraction_tolerance: f64 = 1.0e-12,
    relative_tolerance: f64 = 1.0e-10,
    picard_relaxation: f64 = 0.5,
};

pub const Inputs = struct {
    initial_temperature_k: f64,
    initial_vapor_fraction: f64,
    atmospheric_temperature_k: f64,
    atmospheric_vapor_fraction: f64,
    ground_air_temperature_k: f64,
    ground_air_vapor_fraction: f64,
    heat_capacity_mj_per_k: f64,
    air_volume_m3: f64,
    atmospheric_sensible_conductance_mj_per_h_k: f64,
    atmospheric_vapor_conductance_m3_per_h: f64,
    ground_sensible_conductance_mj_per_h_k: f64,
    ground_vapor_conductance_m3_per_h: f64,
    canopy_surface_sensible_heat_flux_mj_per_h: f64,
    canopy_surface_vapor_flux_m3_per_h: f64,
    lateral_sensible_heat_flux_mj_per_h: f64,
    lateral_vapor_flux_m3_per_h: f64,
};

pub const Result = struct {
    temperature_k: f64,
    vapor_fraction: f64,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    temperature_k: []f64,
    vapor_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize) !State {
        if (cell_count == 0 or species_count == 0) return error.InvalidCanopyAirDimensions;
        const count = try std.math.mul(usize, cell_count, species_count);
        const temperature = try allocator.alloc(f64, count);
        errdefer allocator.free(temperature);
        const vapor = try allocator.alloc(f64, count);
        @memset(temperature, 0);
        @memset(vapor, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .species_count = species_count, .temperature_k = temperature, .vapor_fraction = vapor };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.temperature_k);
        self.allocator.free(self.vapor_fraction);
        self.* = undefined;
    }

    pub fn index(self: State, cell: usize, species: usize) !usize {
        if (cell >= self.cell_count or species >= self.species_count) return error.CanopyAirIndexOutOfBounds;
        return cell * self.species_count + species;
    }
};

pub fn solve(inputs: Inputs, parameters: Parameters, options: SolverOptions) !Result {
    try validate(inputs, parameters, options);
    var temperature_k = inputs.initial_temperature_k;
    var vapor_fraction = inputs.initial_vapor_fraction;
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        const targets = try target(inputs, parameters, temperature_k, vapor_fraction);
        const temperature_residual = targets.temperature_k - temperature_k;
        const vapor_residual = targets.vapor_fraction - vapor_fraction;
        if (converged(temperature_k, vapor_fraction, temperature_residual, vapor_residual, options)) return .{
            .temperature_k = temperature_k,
            .vapor_fraction = vapor_fraction,
            .iterations = iteration + 1,
            .newton_raphson_steps = newton_steps,
            .picard_steps = picard_steps,
        };

        // Away from the saturation cap both balances are linear. At the cap
        // this is a semismooth Newton derivative; the subsequent residual
        // check selects Picard damping if the active branch changes.
        const temperature_derivative = -1.0 -
            (inputs.atmospheric_sensible_conductance_mj_per_h_k + inputs.ground_sensible_conductance_mj_per_h_k) / inputs.heat_capacity_mj_per_k;
        const vapor_derivative = -1.0 -
            (inputs.atmospheric_vapor_conductance_m3_per_h + inputs.ground_vapor_conductance_m3_per_h) / inputs.air_volume_m3;
        const candidate_temperature = temperature_k - temperature_residual / temperature_derivative;
        const candidate_vapor = vapor_fraction - vapor_residual / vapor_derivative;
        if (std.math.isFinite(candidate_temperature) and candidate_temperature > 0 and std.math.isFinite(candidate_vapor) and candidate_vapor >= 0) {
            temperature_k = candidate_temperature;
            vapor_fraction = candidate_vapor;
            newton_steps += 1;
        } else {
            temperature_k += options.picard_relaxation * temperature_residual;
            vapor_fraction = @max(0.0, vapor_fraction + options.picard_relaxation * vapor_residual);
            picard_steps += 1;
        }
    }
    return error.CanopyAirSolverDidNotConverge;
}

pub fn solveInto(state: *State, cell: usize, species: usize, inputs: Inputs, parameters: Parameters, options: SolverOptions) !Result {
    const state_index = try state.index(cell, species);
    const result = try solve(inputs, parameters, options);
    state.temperature_k[state_index] = result.temperature_k;
    state.vapor_fraction[state_index] = result.vapor_fraction;
    return result;
}

const Target = struct { temperature_k: f64, vapor_fraction: f64 };

fn target(inputs: Inputs, parameters: Parameters, temperature_k: f64, vapor_fraction: f64) !Target {
    const sensible_mj_per_h =
        inputs.atmospheric_sensible_conductance_mj_per_h_k * (inputs.atmospheric_temperature_k - temperature_k) -
        inputs.ground_sensible_conductance_mj_per_h_k * (temperature_k - inputs.ground_air_temperature_k) -
        inputs.canopy_surface_sensible_heat_flux_mj_per_h +
        inputs.lateral_sensible_heat_flux_mj_per_h;
    const temperature_target_k = inputs.initial_temperature_k + sensible_mj_per_h / inputs.heat_capacity_mj_per_k;
    if (!std.math.isFinite(temperature_target_k) or temperature_target_k <= 0) return error.InvalidCanopyAirTemperatureTarget;
    const vapor_m3_per_h =
        inputs.atmospheric_vapor_conductance_m3_per_h * (inputs.atmospheric_vapor_fraction - vapor_fraction) -
        inputs.ground_vapor_conductance_m3_per_h * (vapor_fraction - inputs.ground_air_vapor_fraction) -
        inputs.canopy_surface_vapor_flux_m3_per_h +
        inputs.lateral_vapor_flux_m3_per_h;
    const unconstrained_vapor = inputs.initial_vapor_fraction + vapor_m3_per_h / inputs.air_volume_m3;
    const saturation = parameters.saturation_vapor_prefactor_k / temperature_target_k *
        parameters.saturation_relative_humidity *
        @exp(parameters.saturation_temperature_coefficient_k * (parameters.saturation_reference_inverse_temperature_per_k - 1.0 / temperature_target_k));
    if (!std.math.isFinite(unconstrained_vapor) or !std.math.isFinite(saturation) or saturation < 0) return error.InvalidCanopyAirVaporTarget;
    return .{ .temperature_k = temperature_target_k, .vapor_fraction = std.math.clamp(unconstrained_vapor, 0.0, saturation) };
}

fn converged(temperature_k: f64, vapor_fraction: f64, temperature_residual: f64, vapor_residual: f64, options: SolverOptions) bool {
    return @abs(temperature_residual) <= options.absolute_temperature_tolerance_k + options.relative_tolerance * @max(1.0, @abs(temperature_k)) and
        @abs(vapor_residual) <= options.absolute_vapor_fraction_tolerance + options.relative_tolerance * @max(1.0, @abs(vapor_fraction));
}

fn validate(inputs: Inputs, parameters: Parameters, options: SolverOptions) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteCanopyAirInput;
    inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteCanopyAirParameter;
    if (inputs.initial_temperature_k <= 0 or inputs.atmospheric_temperature_k <= 0 or inputs.ground_air_temperature_k <= 0 or inputs.initial_vapor_fraction < 0 or inputs.atmospheric_vapor_fraction < 0 or inputs.ground_air_vapor_fraction < 0 or inputs.heat_capacity_mj_per_k <= 0 or inputs.air_volume_m3 <= 0) return error.InvalidCanopyAirInput;
    inline for (.{ inputs.atmospheric_sensible_conductance_mj_per_h_k, inputs.atmospheric_vapor_conductance_m3_per_h, inputs.ground_sensible_conductance_mj_per_h_k, inputs.ground_vapor_conductance_m3_per_h }) |value| if (value < 0) return error.InvalidCanopyAirInput;
    if (options.max_iterations == 0 or !std.math.isFinite(options.absolute_temperature_tolerance_k) or options.absolute_temperature_tolerance_k <= 0 or !std.math.isFinite(options.absolute_vapor_fraction_tolerance) or options.absolute_vapor_fraction_tolerance <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1) return error.InvalidCanopyAirSolverOptions;
}

test "UPTAKE canopy air balance supports runtime species and exits early" {
    var state = try State.init(std.testing.allocator, 1, 7);
    defer state.deinit();
    const result = try solveInto(&state, 0, 6, .{
        .initial_temperature_k = 290,
        .initial_vapor_fraction = 0.005,
        .atmospheric_temperature_k = 300,
        .atmospheric_vapor_fraction = 0.01,
        .ground_air_temperature_k = 295,
        .ground_air_vapor_fraction = 0.007,
        .heat_capacity_mj_per_k = 1,
        .air_volume_m3 = 10,
        .atmospheric_sensible_conductance_mj_per_h_k = 0.1,
        .atmospheric_vapor_conductance_m3_per_h = 0.2,
        .ground_sensible_conductance_mj_per_h_k = 0.05,
        .ground_vapor_conductance_m3_per_h = 0.1,
        .canopy_surface_sensible_heat_flux_mj_per_h = 0.1,
        .canopy_surface_vapor_flux_m3_per_h = 0.001,
        .lateral_sensible_heat_flux_mj_per_h = 0,
        .lateral_vapor_flux_m3_per_h = 0,
    }, .{}, .{});
    try std.testing.expect(result.iterations < 100);
    try std.testing.expect(state.temperature_k[6] > 290);
    try std.testing.expect(state.vapor_fraction[6] >= 0);
}
