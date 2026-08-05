const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const Options = struct {
    ice_density_megagrams_per_m3: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
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
    /// Signed change in the snow owner's *enthalpy* across the accepted
    /// solve, measured on the same reference state that
    /// `landscape_mass_inventory.aggregateSnow` uses (HEAT-001 resolution A):
    /// liquid water at 0 K, with frozen water one latent heat of fusion
    /// below liquid water at the same temperature.
    ///
    /// For a pure freeze/thaw solve this is zero to rounding, which is the
    /// point: the conversion is internal and no longer needs a boundary
    /// booking. It stays published because the vapor solve shares the
    /// same call site and does move enthalpy.
    enthalpy_change_megajoules: f64,
    /// Deprecated alias of `enthalpy_change_megajoules`, retained only
    /// because the composition root `src/ecosys_ng.zig` reads this name and
    /// is owned by the Integrator lane. See
    /// `docs/binding_requests/heat_001_landscape_enthalpy.md`; delete this
    /// field once the call site is switched. The name is now a misnomer:
    /// booking the sensible change alone would double-count the latent heat
    /// that the enthalpy census already carries.
    sensible_energy_change_megajoules: f64,
};

const melting_temperature_k: f64 = 273.15;
const solid_heat_capacity_megajoules_per_m3_k: f64 = 2.095;
const liquid_heat_capacity_megajoules_per_m3_k: f64 = 4.19;
const ice_heat_capacity_megajoules_per_m3_k: f64 = 1.9274;

/// Local WATSUB freeze/thaw solve. This replaces repeated whole-model snow
/// substeps with safeguarded Newton phase increments and a relaxed Picard
/// fallback. State is committed only after every runtime snow layer is valid.
pub fn solve(allocator: std.mem.Allocator, state: *snow.State, options: Options) !Report {
    if (!std.math.isFinite(options.ice_density_megagrams_per_m3) or options.ice_density_megagrams_per_m3 <= 0 or
        !std.math.isFinite(options.latent_heat_of_fusion_megajoules_per_m3) or options.latent_heat_of_fusion_megajoules_per_m3 <= 0 or
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
    const heat_capacity = try allocator.dupe(f64, state.heat_capacity_megajoules_per_k);
    defer allocator.free(heat_capacity);

    var report: Report = .{ .iterations = 0, .converged = false, .maximum_temperature_residual_k = 0, .enthalpy_change_megajoules = 0, .sensible_energy_change_megajoules = 0 };
    for (1..@as(usize, options.max_iterations) + 1) |iteration| {
        var maximum_residual: f64 = 0;
        var changed = false;
        for (solid, liquid, ice, temperature, heat_capacity, state.active, 0..) |*solid_m3, *liquid_m3, *ice_m3, *temperature_k, *capacity, active, layer_index| {
            inline for (.{ solid_m3.*, liquid_m3.*, ice_m3.*, temperature_k.*, capacity.* }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowPhaseState;
            if (solid_m3.* < 0 or liquid_m3.* < 0 or ice_m3.* < 0 or temperature_k.* <= 0 or capacity.* < 0) return error.InvalidSnowPhaseState;
            if (!active or capacity.* <= 0) continue;
            const frozen_water_equivalent_m3 = solid_m3.* + ice_m3.* * options.ice_density_megagrams_per_m3;
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
            const old_sensible_energy_megajoules = capacity.* * temperature_k.*;
            const solid_fraction = if (frozen_water_equivalent_m3 > 0) solid_m3.* / frozen_water_equivalent_m3 else 0;
            const ice_water_equivalent_fraction = 1 - solid_fraction;
            const heat_capacity_change_per_m3 = if (can_thaw)
                liquid_heat_capacity_megajoules_per_m3_k - solid_heat_capacity_megajoules_per_m3_k * solid_fraction - ice_heat_capacity_megajoules_per_m3_k * ice_water_equivalent_fraction / options.ice_density_megagrams_per_m3
            else
                -liquid_heat_capacity_megajoules_per_m3_k + ice_heat_capacity_megajoules_per_m3_k / options.ice_density_megagrams_per_m3;
            const heat_capacity_change_per_megajoule = if (can_thaw) -heat_capacity_change_per_m3 / options.latent_heat_of_fusion_megajoules_per_m3 else heat_capacity_change_per_m3 / options.latent_heat_of_fusion_megajoules_per_m3;
            const derivative = 1 - melting_temperature_k * heat_capacity_change_per_megajoule;
            const newton_heat_megajoules = (melting_temperature_k * capacity.* - old_sensible_energy_megajoules) / derivative;
            const requested_heat_megajoules = if (std.math.isFinite(newton_heat_megajoules) and @abs(derivative) > 1e-12) newton_heat_megajoules else (melting_temperature_k - temperature_k.*) / options.damping_divisor * capacity.* * options.picard_relaxation;
            if (!std.math.isFinite(requested_heat_megajoules)) return error.NonFiniteSnowPhaseStep;
            const available_heat_megajoules = if (requested_heat_megajoules < 0)
                options.latent_heat_of_fusion_megajoules_per_m3 * frozen_water_equivalent_m3
            else
                options.latent_heat_of_fusion_megajoules_per_m3 * liquid_m3.*;
            const phase_heat_megajoules = std.math.clamp(requested_heat_megajoules, -available_heat_megajoules, available_heat_megajoules);
            if (phase_heat_megajoules == 0) continue;
            if (phase_heat_megajoules < 0) {
                const thawed_water_equivalent_m3 = -phase_heat_megajoules / options.latent_heat_of_fusion_megajoules_per_m3;
                const from_solid_m3 = thawed_water_equivalent_m3 * solid_fraction;
                const from_ice_water_equivalent_m3 = thawed_water_equivalent_m3 - from_solid_m3;
                solid_m3.* = @max(0, solid_m3.* - from_solid_m3);
                ice_m3.* = @max(0, ice_m3.* - from_ice_water_equivalent_m3 / options.ice_density_megagrams_per_m3);
                liquid_m3.* += thawed_water_equivalent_m3;
            } else {
                const frozen_liquid_m3 = phase_heat_megajoules / options.latent_heat_of_fusion_megajoules_per_m3;
                liquid_m3.* = @max(0, liquid_m3.* - frozen_liquid_m3);
                ice_m3.* += frozen_liquid_m3 / options.ice_density_megagrams_per_m3;
            }
            capacity.* = solid_heat_capacity_megajoules_per_m3_k * solid_m3.* + liquid_heat_capacity_megajoules_per_m3_k * (liquid_m3.* + state.vapor_water_equivalent_m3[layer_index]) + ice_heat_capacity_megajoules_per_m3_k * ice_m3.*;
            if (!std.math.isFinite(capacity.*) or capacity.* <= 0) return error.InvalidSnowPhaseHeatCapacity;
            temperature_k.* = (old_sensible_energy_megajoules + phase_heat_megajoules) / capacity.*;
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
    for (
        heat_capacity,
        temperature,
        solid,
        ice,
        state.heat_capacity_megajoules_per_k,
        state.temperature_k,
        state.solid_snow_water_equivalent_m3,
        state.ice_volume_m3,
    ) |
        next_capacity,
        next_temperature,
        next_solid,
        next_ice,
        previous_capacity,
        previous_temperature,
        previous_solid,
        previous_ice,
    | {
        const next_frozen_water_equivalent_m3 =
            next_solid + next_ice * options.ice_density_megagrams_per_m3;
        const previous_frozen_water_equivalent_m3 =
            previous_solid + previous_ice * options.ice_density_megagrams_per_m3;
        report.enthalpy_change_megajoules +=
            next_capacity * next_temperature -
            previous_capacity * previous_temperature -
            options.latent_heat_of_fusion_megajoules_per_m3 *
                (next_frozen_water_equivalent_m3 -
                    previous_frozen_water_equivalent_m3);
        if (!std.math.isFinite(report.enthalpy_change_megajoules))
            return error.NonFiniteSnowPhaseEnergyChange;
    }
    report.sensible_energy_change_megajoules = report.enthalpy_change_megajoules;
    @memcpy(state.solid_snow_water_equivalent_m3, solid);
    @memcpy(state.liquid_water_volume_m3, liquid);
    @memcpy(state.ice_volume_m3, ice);
    @memcpy(state.temperature_k, temperature);
    @memcpy(state.heat_capacity_megajoules_per_k, heat_capacity);
    state.refreshAllGeometry();
    return report;
}

test "freeze thaw conserves water equivalent and commits atomically" {
    var state = try snow.State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    try state.initializePhysicalState(&.{ 0.1, 0.1 }, &.{ 1, 1 }, &.{ 275, 270 }, &.{0.2}, 0.05);
    state.temperature_k[0] = 275;
    state.liquid_water_volume_m3[1] = 0.002;
    state.heat_capacity_megajoules_per_k[1] += liquid_heat_capacity_megajoules_per_m3_k * 0.002;
    const before0 = state.solid_snow_water_equivalent_m3[0] + state.liquid_water_volume_m3[0] + state.ice_volume_m3[0] * 0.92;
    const before1 = state.solid_snow_water_equivalent_m3[1] + state.liquid_water_volume_m3[1] + state.ice_volume_m3[1] * 0.92;
    const report = try solve(std.testing.allocator, &state, .{ .ice_density_megagrams_per_m3 = 0.92, .latent_heat_of_fusion_megajoules_per_m3 = 333, .damping_divisor = 2.7185, .absolute_temperature_tolerance_k = 1e-8, .relative_tolerance = 1e-8, .picard_relaxation = 1, .max_iterations = 20 });
    try std.testing.expect(report.iterations < 20);
    try std.testing.expectApproxEqAbs(before0, state.solid_snow_water_equivalent_m3[0] + state.liquid_water_volume_m3[0] + state.ice_volume_m3[0] * 0.92, 1e-12);
    try std.testing.expectApproxEqAbs(before1, state.solid_snow_water_equivalent_m3[1] + state.liquid_water_volume_m3[1] + state.ice_volume_m3[1] * 0.92, 1e-12);
    try std.testing.expect(state.liquid_water_volume_m3[0] > 0);
    try std.testing.expect(state.ice_volume_m3[1] > 0);
    try std.testing.expect(std.math.isFinite(report.sensible_energy_change_megajoules));
}
