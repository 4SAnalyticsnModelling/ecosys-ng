const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const Options = struct {
    ice_density_Mg_per_m3: f64,
    latent_heat_of_fusion_mj_per_m3: f64,
    damping_divisor: f64,
    absolute_temperature_tolerance_k: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    max_iterations: u16,
};

pub const Report = struct {
    iterations: u16,
    converged: bool,
    maximum_temperature_residual_k: f64,
    sensible_energy_change_mj: f64,
};

const melting_temperature_k: f64 = 273.15;
const solid_heat_capacity_mj_per_m3_k: f64 = 2.095;
const liquid_heat_capacity_mj_per_m3_k: f64 = 4.19;
const ice_heat_capacity_mj_per_m3_k: f64 = 1.9274;

/// Local WATSUB freeze/thaw solve. This replaces repeated whole-model snow
/// substeps with safeguarded Newton phase increments and a relaxed Picard
/// fallback. State is committed only after every runtime snow layer is valid.
pub fn solve(allocator: std.mem.Allocator, state: *snow.State, options: Options) !Report {
    if (!std.math.isFinite(options.ice_density_Mg_per_m3) or options.ice_density_Mg_per_m3 <= 0 or
        !std.math.isFinite(options.latent_heat_of_fusion_mj_per_m3) or options.latent_heat_of_fusion_mj_per_m3 <= 0 or
        !std.math.isFinite(options.damping_divisor) or options.damping_divisor <= 0 or
        !std.math.isFinite(options.absolute_temperature_tolerance_k) or options.absolute_temperature_tolerance_k < 0 or
        !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance < 0 or
        !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or
        options.max_iterations == 0) return error.InvalidSnowPhaseSolverOptions;

    const solid = try allocator.dupe(f64, state.solid_snow_water_equivalent_m3);
    defer allocator.free(solid);
    const liquid = try allocator.dupe(f64, state.liquid_water_volume_m3);
    defer allocator.free(liquid);
    const ice = try allocator.dupe(f64, state.ice_volume_m3);
    defer allocator.free(ice);
    const temperature = try allocator.dupe(f64, state.temperature_k);
    defer allocator.free(temperature);
    const heat_capacity = try allocator.dupe(f64, state.heat_capacity_mj_per_k);
    defer allocator.free(heat_capacity);

    var report: Report = .{ .iterations = 0, .converged = false, .maximum_temperature_residual_k = 0, .sensible_energy_change_mj = 0 };
    for (1..@as(usize, options.max_iterations) + 1) |iteration| {
        var maximum_residual: f64 = 0;
        var changed = false;
        for (solid, liquid, ice, temperature, heat_capacity, state.active, 0..) |*solid_m3, *liquid_m3, *ice_m3, *temperature_k, *capacity, active, layer_index| {
            inline for (.{ solid_m3.*, liquid_m3.*, ice_m3.*, temperature_k.*, capacity.* }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowPhaseState;
            if (solid_m3.* < 0 or liquid_m3.* < 0 or ice_m3.* < 0 or temperature_k.* <= 0 or capacity.* < 0) return error.InvalidSnowPhaseState;
            if (!active or capacity.* <= 0) continue;
            const frozen_water_equivalent_m3 = solid_m3.* + ice_m3.* * options.ice_density_Mg_per_m3;
            const can_thaw = temperature_k.* > melting_temperature_k and frozen_water_equivalent_m3 > 0;
            const can_freeze = temperature_k.* < melting_temperature_k and liquid_m3.* > 0;
            if (!can_thaw and !can_freeze) continue;

            const residual_k = @abs(temperature_k.* - melting_temperature_k);
            maximum_residual = @max(maximum_residual, residual_k);
            const tolerance_k = options.absolute_temperature_tolerance_k + options.relative_tolerance * melting_temperature_k;
            if (residual_k <= tolerance_k) continue;

            // Solve E + H - Tm*C(H) = 0. Within either phase branch C(H) is
            // affine, so this is the exact Newton step; inventory clipping is
            // the safeguard. WATSUB damping plus Picard relaxation is retained
            // as the fallback for a singular derivative.
            const old_sensible_energy_mj = capacity.* * temperature_k.*;
            const solid_fraction = if (frozen_water_equivalent_m3 > 0) solid_m3.* / frozen_water_equivalent_m3 else 0;
            const ice_water_equivalent_fraction = 1 - solid_fraction;
            const heat_capacity_change_per_m3 = if (can_thaw)
                liquid_heat_capacity_mj_per_m3_k - solid_heat_capacity_mj_per_m3_k * solid_fraction - ice_heat_capacity_mj_per_m3_k * ice_water_equivalent_fraction / options.ice_density_Mg_per_m3
            else
                -liquid_heat_capacity_mj_per_m3_k + ice_heat_capacity_mj_per_m3_k / options.ice_density_Mg_per_m3;
            const heat_capacity_change_per_mj = if (can_thaw) -heat_capacity_change_per_m3 / options.latent_heat_of_fusion_mj_per_m3 else heat_capacity_change_per_m3 / options.latent_heat_of_fusion_mj_per_m3;
            const derivative = 1 - melting_temperature_k * heat_capacity_change_per_mj;
            const newton_heat_mj = (melting_temperature_k * capacity.* - old_sensible_energy_mj) / derivative;
            const requested_heat_mj = if (std.math.isFinite(newton_heat_mj) and @abs(derivative) > 1e-12) newton_heat_mj else (melting_temperature_k - temperature_k.*) / options.damping_divisor * capacity.* * options.picard_relaxation;
            if (!std.math.isFinite(requested_heat_mj)) return error.NonFiniteSnowPhaseStep;
            const available_heat_mj = if (requested_heat_mj < 0)
                options.latent_heat_of_fusion_mj_per_m3 * frozen_water_equivalent_m3
            else
                options.latent_heat_of_fusion_mj_per_m3 * liquid_m3.*;
            const phase_heat_mj = std.math.clamp(requested_heat_mj, -available_heat_mj, available_heat_mj);
            if (phase_heat_mj == 0) continue;
            if (phase_heat_mj < 0) {
                const thawed_water_equivalent_m3 = -phase_heat_mj / options.latent_heat_of_fusion_mj_per_m3;
                const from_solid_m3 = thawed_water_equivalent_m3 * solid_fraction;
                const from_ice_water_equivalent_m3 = thawed_water_equivalent_m3 - from_solid_m3;
                solid_m3.* = @max(0, solid_m3.* - from_solid_m3);
                ice_m3.* = @max(0, ice_m3.* - from_ice_water_equivalent_m3 / options.ice_density_Mg_per_m3);
                liquid_m3.* += thawed_water_equivalent_m3;
            } else {
                const frozen_liquid_m3 = phase_heat_mj / options.latent_heat_of_fusion_mj_per_m3;
                liquid_m3.* = @max(0, liquid_m3.* - frozen_liquid_m3);
                ice_m3.* += frozen_liquid_m3 / options.ice_density_Mg_per_m3;
            }
            capacity.* = solid_heat_capacity_mj_per_m3_k * solid_m3.* + liquid_heat_capacity_mj_per_m3_k * (liquid_m3.* + state.vapor_water_equivalent_m3[layer_index]) + ice_heat_capacity_mj_per_m3_k * ice_m3.*;
            if (!std.math.isFinite(capacity.*) or capacity.* <= 0) return error.InvalidSnowPhaseHeatCapacity;
            temperature_k.* = (old_sensible_energy_mj + phase_heat_mj) / capacity.*;
            if (!std.math.isFinite(temperature_k.*) or temperature_k.* <= 0) return error.InvalidSnowPhaseTemperature;
            changed = true;
        }
        report.iterations = @intCast(iteration);
        report.maximum_temperature_residual_k = maximum_residual;
        if (!changed) {
            report.converged = true;
            break;
        }
    }
    if (!report.converged) return error.SnowPhaseSolverDidNotConverge;
    for (heat_capacity, temperature, state.heat_capacity_mj_per_k, state.temperature_k) |next_capacity, next_temperature, previous_capacity, previous_temperature| {
        report.sensible_energy_change_mj +=
            next_capacity * next_temperature -
            previous_capacity * previous_temperature;
        if (!std.math.isFinite(report.sensible_energy_change_mj))
            return error.NonFiniteSnowPhaseEnergyChange;
    }
    @memcpy(state.solid_snow_water_equivalent_m3, solid);
    @memcpy(state.liquid_water_volume_m3, liquid);
    @memcpy(state.ice_volume_m3, ice);
    @memcpy(state.temperature_k, temperature);
    @memcpy(state.heat_capacity_mj_per_k, heat_capacity);
    state.refreshAllGeometry();
    return report;
}

test "freeze thaw conserves water equivalent and commits atomically" {
    var state = try snow.State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    try state.initializePhysicalState(&.{ 0.1, 0.1 }, &.{ 1, 1 }, &.{ 275, 270 }, &.{0.2}, 0.05);
    state.temperature_k[0] = 275;
    state.liquid_water_volume_m3[1] = 0.002;
    state.heat_capacity_mj_per_k[1] += liquid_heat_capacity_mj_per_m3_k * 0.002;
    const before0 = state.solid_snow_water_equivalent_m3[0] + state.liquid_water_volume_m3[0] + state.ice_volume_m3[0] * 0.92;
    const before1 = state.solid_snow_water_equivalent_m3[1] + state.liquid_water_volume_m3[1] + state.ice_volume_m3[1] * 0.92;
    const report = try solve(std.testing.allocator, &state, .{ .ice_density_Mg_per_m3 = 0.92, .latent_heat_of_fusion_mj_per_m3 = 333, .damping_divisor = 2.7185, .absolute_temperature_tolerance_k = 1e-8, .relative_tolerance = 1e-8, .picard_relaxation = 1, .max_iterations = 20 });
    try std.testing.expect(report.iterations < 20);
    try std.testing.expectApproxEqAbs(before0, state.solid_snow_water_equivalent_m3[0] + state.liquid_water_volume_m3[0] + state.ice_volume_m3[0] * 0.92, 1e-12);
    try std.testing.expectApproxEqAbs(before1, state.solid_snow_water_equivalent_m3[1] + state.liquid_water_volume_m3[1] + state.ice_volume_m3[1] * 0.92, 1e-12);
    try std.testing.expect(state.liquid_water_volume_m3[0] > 0);
    try std.testing.expect(state.ice_volume_m3[1] > 0);
    try std.testing.expect(std.math.isFinite(report.sensible_energy_change_mj));
}
