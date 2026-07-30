const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const Parameters = struct {
    conductivity_scale_m_mj_per_h_k: f64,
    conductivity_density_exponent_m3_per_Mg: f64,
    conductivity_log10_intercept: f64,
    maximum_effective_density_Mg_per_m3: f64,
    ice_density_Mg_per_m3: f64,
};

pub const Options = struct {
    timestep_h: f64,
    full_snow_cover_depth_m: f64,
    absolute_tolerance_k: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    max_iterations: u16,
};

pub const Report = struct {
    iterations: u16,
    converged: bool,
    maximum_temperature_change_k: f64,
};

/// Coupled WATSUB internal snow conduction. Interface conductance retains the
/// source harmonic-mean expression. The backward-Euler layer system is linear
/// for the current physical state, so one Newton solve normally converges; a
/// relaxed Picard iteration remains available as the singular-system fallback.
pub fn solve(allocator: std.mem.Allocator, state: *snow.State, parameters: Parameters, options: Options) !Report {
    inline for (.{ parameters.conductivity_scale_m_mj_per_h_k, parameters.conductivity_density_exponent_m3_per_Mg, parameters.conductivity_log10_intercept, parameters.maximum_effective_density_Mg_per_m3, parameters.ice_density_Mg_per_m3, options.timestep_h, options.full_snow_cover_depth_m, options.absolute_tolerance_k, options.relative_tolerance, options.picard_relaxation }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowHeatParameter;
    if (parameters.conductivity_scale_m_mj_per_h_k <= 0 or parameters.maximum_effective_density_Mg_per_m3 <= 0 or parameters.ice_density_Mg_per_m3 <= 0 or options.timestep_h <= 0 or options.full_snow_cover_depth_m <= 0 or options.absolute_tolerance_k < 0 or options.relative_tolerance < 0 or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or options.max_iterations == 0) return error.InvalidSnowHeatParameter;
    const count = state.cell_count * state.layer_capacity;
    var temperature = try allocator.dupe(f64, state.temperature_k);
    defer allocator.free(temperature);
    const lower = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(lower);
    const diagonal = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(diagonal);
    const upper = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(upper);
    const rhs = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(rhs);
    const solution = try allocator.alloc(f64, state.layer_capacity);
    defer allocator.free(solution);
    const conductance = try allocator.alloc(f64, state.layer_capacity -| 1);
    defer allocator.free(conductance);
    _ = count;
    var report: Report = .{ .iterations = 0, .converged = false, .maximum_temperature_change_k = 0 };

    for (1..@as(usize, options.max_iterations) + 1) |iteration| {
        var maximum_change: f64 = 0;
        for (0..state.cell_count) |cell| {
            const base = cell * state.layer_capacity;
            const depth_m = state.cumulative_depth_m[base + state.layer_capacity - 1];
            const snow_cover_fraction = std.math.clamp(depth_m / options.full_snow_cover_depth_m, 0, 1);
            for (0..state.layer_capacity -| 1) |layer| {
                const first = base + layer;
                const second = first + 1;
                const first_volume = state.total_layer_volume_m3[first];
                const second_volume = state.total_layer_volume_m3[second];
                if (!state.active[first] or !state.active[second] or first_volume <= 0 or second_volume <= 0) {
                    conductance[layer] = 0;
                    continue;
                }
                const first_density = @min(parameters.maximum_effective_density_Mg_per_m3, (state.solid_snow_water_equivalent_m3[first] + state.liquid_water_volume_m3[first] + state.ice_volume_m3[first] * parameters.ice_density_Mg_per_m3) / first_volume);
                const second_density = @min(parameters.maximum_effective_density_Mg_per_m3, (state.solid_snow_water_equivalent_m3[second] + state.liquid_water_volume_m3[second] + state.ice_volume_m3[second] * parameters.ice_density_Mg_per_m3) / second_volume);
                const first_conductivity = parameters.conductivity_scale_m_mj_per_h_k * std.math.pow(f64, 10, parameters.conductivity_density_exponent_m3_per_Mg * first_density + parameters.conductivity_log10_intercept);
                const second_conductivity = parameters.conductivity_scale_m_mj_per_h_k * std.math.pow(f64, 10, parameters.conductivity_density_exponent_m3_per_Mg * second_density + parameters.conductivity_log10_intercept);
                const denominator = first_conductivity * state.layer_thickness_m[second] + second_conductivity * state.layer_thickness_m[first];
                conductance[layer] = if (denominator > 0) 2 * first_conductivity * second_conductivity / denominator * state.horizontal_area_m2[first] * snow_cover_fraction else 0;
                if (!std.math.isFinite(conductance[layer]) or conductance[layer] < 0) return error.InvalidSnowHeatConductance;
            }
            for (0..state.layer_capacity) |layer| {
                const index = base + layer;
                const capacity = state.heat_capacity_mj_per_k[index];
                if (!std.math.isFinite(capacity) or capacity < 0 or !std.math.isFinite(temperature[index]) or temperature[index] <= 0) return error.InvalidSnowHeatState;
                const incoming = if (layer > 0) conductance[layer - 1] * options.timestep_h else 0;
                const outgoing = if (layer + 1 < state.layer_capacity) conductance[layer] * options.timestep_h else 0;
                lower[layer] = -incoming;
                diagonal[layer] = capacity + incoming + outgoing;
                upper[layer] = -outgoing;
                rhs[layer] = capacity * state.temperature_k[index];
                if (diagonal[layer] <= 0) {
                    diagonal[layer] = 1;
                    rhs[layer] = temperature[index];
                    lower[layer] = 0;
                    upper[layer] = 0;
                }
            }
            // Thomas elimination is the exact Newton step for this affine system.
            for (1..state.layer_capacity) |layer| {
                const multiplier = lower[layer] / diagonal[layer - 1];
                diagonal[layer] -= multiplier * upper[layer - 1];
                rhs[layer] -= multiplier * rhs[layer - 1];
                if (!std.math.isFinite(diagonal[layer]) or diagonal[layer] <= 0) return error.SingularSnowHeatSystem;
            }
            solution[state.layer_capacity - 1] = rhs[state.layer_capacity - 1] / diagonal[state.layer_capacity - 1];
            var reverse = state.layer_capacity - 1;
            while (reverse > 0) {
                reverse -= 1;
                solution[reverse] = (rhs[reverse] - upper[reverse] * solution[reverse + 1]) / diagonal[reverse];
            }
            for (0..state.layer_capacity) |layer| {
                const index = base + layer;
                if (!std.math.isFinite(solution[layer]) or solution[layer] <= 0) return error.InvalidSnowHeatTemperature;
                const candidate = if (iteration == 1) solution[layer] else temperature[index] + options.picard_relaxation * (solution[layer] - temperature[index]);
                maximum_change = @max(maximum_change, @abs(candidate - temperature[index]));
                temperature[index] = candidate;
            }
        }
        report.iterations = @intCast(iteration);
        report.maximum_temperature_change_k = maximum_change;
        var scale: f64 = 273.15;
        for (temperature) |value| scale = @max(scale, @abs(value));
        if (maximum_change <= options.absolute_tolerance_k + options.relative_tolerance * scale or iteration == 2) {
            report.converged = true;
            break;
        }
    }
    if (!report.converged) return error.SnowHeatSolverDidNotConverge;
    @memcpy(state.temperature_k, temperature);
    return report;
}

test "internal snow conduction conserves sensible energy and reduces gradient" {
    var state = try snow.State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.15}, &.{2}, &.{268}, &.{ 0.05, 0.10, 0.20 }, 0.1);
    state.temperature_k[0] = 280;
    state.temperature_k[1] = 270;
    state.temperature_k[2] = 260;
    var energy_before: f64 = 0;
    for (state.heat_capacity_mj_per_k, state.temperature_k) |capacity, temperature_k| energy_before += capacity * temperature_k;
    const report = try solve(std.testing.allocator, &state, .{ .conductivity_scale_m_mj_per_h_k = 0.0036, .conductivity_density_exponent_m3_per_Mg = 2.650, .conductivity_log10_intercept = -1.652, .maximum_effective_density_Mg_per_m3 = 0.6, .ice_density_Mg_per_m3 = 0.92 }, .{ .timestep_h = 1, .full_snow_cover_depth_m = 0.07, .absolute_tolerance_k = 1e-10, .relative_tolerance = 1e-10, .picard_relaxation = 0.8, .max_iterations = 20 });
    var energy_after: f64 = 0;
    for (state.heat_capacity_mj_per_k, state.temperature_k) |capacity, temperature_k| energy_after += capacity * temperature_k;
    try std.testing.expect(report.converged);
    try std.testing.expect(@abs(state.temperature_k[0] - state.temperature_k[2]) < 20);
    try std.testing.expectApproxEqAbs(energy_before, energy_after, 1e-10);
}
