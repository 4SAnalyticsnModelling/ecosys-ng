const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const Parameters = struct {
    reference_vapor_diffusivity_m2_per_h: f64,
    reference_temperature_k: f64,
    temperature_exponent: f64,
    minimum_air_fraction: f64,
    vapor_sensible_heat_capacity_mj_per_m3_k: f64,
};

pub const Options = struct {
    timestep_h: f64,
    full_snow_cover_depth_m: f64,
    absolute_tolerance_m3: f64,
    relative_tolerance: f64,
    max_iterations: u16,
};

pub const Report = struct { iterations: u16, converged: bool, maximum_interface_flux_m3: f64 };

/// Implicit WATSUB interlayer vapor diffusion. Conductivity is air-fraction
/// squared times temperature-adjusted diffusivity; interface conductance is
/// the exact harmonic expression. Signed carrier heat uses donor temperature.
pub fn solve(allocator: std.mem.Allocator, state: *snow.State, parameters: Parameters, options: Options) !Report {
    inline for (.{ parameters.reference_vapor_diffusivity_m2_per_h, parameters.reference_temperature_k, parameters.temperature_exponent, parameters.minimum_air_fraction, parameters.vapor_sensible_heat_capacity_mj_per_m3_k, options.timestep_h, options.full_snow_cover_depth_m, options.absolute_tolerance_m3, options.relative_tolerance }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowVaporDiffusionParameter;
    if (parameters.reference_vapor_diffusivity_m2_per_h < 0 or parameters.reference_temperature_k <= 0 or parameters.temperature_exponent < 0 or parameters.minimum_air_fraction < 0 or parameters.minimum_air_fraction > 1 or parameters.vapor_sensible_heat_capacity_mj_per_m3_k <= 0 or options.timestep_h <= 0 or options.full_snow_cover_depth_m <= 0 or options.absolute_tolerance_m3 < 0 or options.relative_tolerance < 0 or options.max_iterations == 0) return error.InvalidSnowVaporDiffusionParameter;
    const vapor = try allocator.dupe(f64, state.vapor_water_equivalent_m3);
    defer allocator.free(vapor);
    const temperature = try allocator.dupe(f64, state.temperature_k);
    defer allocator.free(temperature);
    const sensible_energy_mj = try allocator.alloc(f64, state.temperature_k.len);
    defer allocator.free(sensible_energy_mj);
    for (sensible_energy_mj, state.heat_capacity_mj_per_k, state.temperature_k) |*energy, capacity, temperature_k| energy.* = capacity * temperature_k;
    const lower = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(lower);
    const diagonal = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(diagonal);
    const upper = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(upper);
    const rhs = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(rhs);
    const concentration = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(concentration);
    const conductance = try allocator.alloc(f64, state.layer_capacity -| 1);
    defer allocator.free(conductance);
    var maximum_flux: f64 = 0;
    for (0..state.cell_count) |cell| {
        const base = cell * state.layer_capacity;
        const cover = std.math.clamp(state.cumulative_depth_m[base + state.layer_capacity - 1] / options.full_snow_cover_depth_m, 0, 1);
        for (0..state.layer_capacity -| 1) |layer| {
            const first = base + layer;
            const second = first + 1;
            if (!state.active[first] or !state.active[second] or state.air_filled_volume_m3[first] <= 0 or state.air_filled_volume_m3[second] <= 0 or state.total_layer_volume_m3[first] <= 0 or state.total_layer_volume_m3[second] <= 0) {
                conductance[layer] = 0;
                continue;
            }
            const first_air_fraction = @max(parameters.minimum_air_fraction, state.air_filled_volume_m3[first] / state.total_layer_volume_m3[first]);
            const second_air_fraction = @max(parameters.minimum_air_fraction, state.air_filled_volume_m3[second] / state.total_layer_volume_m3[second]);
            const first_diffusivity = parameters.reference_vapor_diffusivity_m2_per_h * std.math.pow(f64, temperature[first] / parameters.reference_temperature_k, parameters.temperature_exponent);
            const second_diffusivity = parameters.reference_vapor_diffusivity_m2_per_h * std.math.pow(f64, temperature[second] / parameters.reference_temperature_k, parameters.temperature_exponent);
            const first_conductivity = first_air_fraction * first_air_fraction * first_diffusivity;
            const second_conductivity = second_air_fraction * second_air_fraction * second_diffusivity;
            const denominator = first_conductivity * state.layer_thickness_m[second] + second_conductivity * state.layer_thickness_m[first];
            conductance[layer] = if (denominator > 0) 2 * first_conductivity * second_conductivity / denominator * state.horizontal_area_m2[first] * cover else 0;
            if (!std.math.isFinite(conductance[layer]) or conductance[layer] < 0) return error.InvalidSnowVaporDiffusionConductance;
        }
        for (0..state.layer_capacity) |layer| {
            const index = base + layer;
            const storage = state.air_filled_volume_m3[index];
            if (!std.math.isFinite(vapor[index]) or vapor[index] < 0 or !std.math.isFinite(storage) or storage < 0) return error.InvalidSnowVaporDiffusionState;
            const incoming = if (layer > 0) conductance[layer - 1] * options.timestep_h else 0;
            const outgoing = if (layer + 1 < state.layer_capacity) conductance[layer] * options.timestep_h else 0;
            lower[layer] = -incoming;
            diagonal[layer] = storage + incoming + outgoing;
            upper[layer] = -outgoing;
            rhs[layer] = vapor[index];
            if (storage <= 0) {
                lower[layer] = 0;
                diagonal[layer] = 1;
                upper[layer] = 0;
                rhs[layer] = 0;
            }
        }
        for (1..state.layer_capacity) |layer| {
            const factor = lower[layer] / diagonal[layer - 1];
            diagonal[layer] -= factor * upper[layer - 1];
            rhs[layer] -= factor * rhs[layer - 1];
            if (!std.math.isFinite(diagonal[layer]) or diagonal[layer] <= 0) return error.SingularSnowVaporDiffusionSystem;
        }
        concentration[state.layer_capacity - 1] = rhs[state.layer_capacity - 1] / diagonal[state.layer_capacity - 1];
        var reverse = state.layer_capacity - 1;
        while (reverse > 0) {
            reverse -= 1;
            concentration[reverse] = (rhs[reverse] - upper[reverse] * concentration[reverse + 1]) / diagonal[reverse];
        }
        for (0..state.layer_capacity) |layer| {
            const index = base + layer;
            vapor[index] = concentration[layer] * state.air_filled_volume_m3[index];
            if (!std.math.isFinite(vapor[index]) or vapor[index] < -options.absolute_tolerance_m3) return error.InvalidSnowVaporDiffusionResult;
            vapor[index] = @max(0, vapor[index]);
        }
        for (0..state.layer_capacity -| 1) |layer| {
            const first = base + layer;
            const second = first + 1;
            const flux_m3 = conductance[layer] * options.timestep_h * (concentration[layer] - concentration[layer + 1]);
            maximum_flux = @max(maximum_flux, @abs(flux_m3));
            const donor = if (flux_m3 >= 0) first else second;
            const heat_mj = parameters.vapor_sensible_heat_capacity_mj_per_m3_k * temperature[donor] * flux_m3;
            sensible_energy_mj[first] -= heat_mj;
            sensible_energy_mj[second] += heat_mj;
        }
    }
    for (temperature, sensible_energy_mj, state.heat_capacity_mj_per_k) |*temperature_k, energy_mj, capacity| if (capacity > 0) {
        temperature_k.* = energy_mj / capacity;
    };
    for (temperature) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSnowVaporDiffusionTemperature;
    @memcpy(state.vapor_water_equivalent_m3, vapor);
    @memcpy(state.temperature_k, temperature);
    return .{ .iterations = 1, .converged = true, .maximum_interface_flux_m3 = maximum_flux };
}

test "interlayer vapor diffusion conserves vapor and sensible energy" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.1}, &.{1}, &.{268}, &.{ 0.05, 0.10 }, 0.1);
    state.vapor_water_equivalent_m3[0] = 1e-4;
    state.vapor_water_equivalent_m3[1] = 0;
    const vapor_before = state.vapor_water_equivalent_m3[0];
    var energy_before: f64 = 0;
    for (state.heat_capacity_mj_per_k, state.temperature_k) |capacity, temperature_k| energy_before += capacity * temperature_k;
    const report = try solve(std.testing.allocator, &state, .{ .reference_vapor_diffusivity_m2_per_h = 0.0896, .reference_temperature_k = 298.15, .temperature_exponent = 1.75, .minimum_air_fraction = 0, .vapor_sensible_heat_capacity_mj_per_m3_k = 4.19 }, .{ .timestep_h = 1, .full_snow_cover_depth_m = 0.07, .absolute_tolerance_m3 = 1e-14, .relative_tolerance = 1e-8, .max_iterations = 20 });
    var energy_after: f64 = 0;
    for (state.heat_capacity_mj_per_k, state.temperature_k) |capacity, temperature_k| energy_after += capacity * temperature_k;
    try std.testing.expect(report.maximum_interface_flux_m3 > 0);
    try std.testing.expect(state.vapor_water_equivalent_m3[1] > 0);
    try std.testing.expectApproxEqAbs(vapor_before, state.vapor_water_equivalent_m3[0] + state.vapor_water_equivalent_m3[1], 1e-12);
    try std.testing.expectApproxEqAbs(energy_before, energy_after, 1e-10);
}
