const std = @import("std");
const snow = @import("../solute/snow_solute_transport.zig");

pub const Parameters = struct {
    vapor_volume_prefactor_k: f64,
    equilibrium_relative_humidity: f64,
    clausius_clapeyron_temperature_k: f64,
    reference_inverse_temperature_per_k: f64,
    liquid_evaporation_latent_heat_megajoules_per_m3: f64,
    snow_sublimation_latent_heat_megajoules_per_m3: f64,
    /// HEAT-001 resolution A. Latent heat of fusion, needed only to express
    /// this solve's energy change on the same enthalpy reference state that
    /// `landscape_mass_inventory.aggregateSnow` uses. Sublimation moves water
    /// from the solid carrier, which the census holds one latent heat of
    /// fusion below liquid, so the frozen-mass change must be included.
    ///
    /// The default is a temporary bridge because the composition root
    /// `src/ecosys_ng.zig` builds this struct and is owned by the Integrator
    /// lane; see `docs/binding_requests/heat_001_landscape_enthalpy.md`.
    latent_heat_of_fusion_megajoules_per_m3: f64 = 333,
};

pub const Options = struct {
    absolute_tolerance_m3: f64,
    relative_tolerance: f64,
    max_iterations: u16,
};

pub const Report = struct {
    iterations: u16,
    converged: bool,
    maximum_vapor_residual_m3: f64,
    /// Signed change in the snow owner's enthalpy, on the reference state
    /// described on `Parameters.latent_heat_of_fusion_megajoules_per_m3`.
    enthalpy_change_megajoules: f64,
    /// Deprecated alias of `enthalpy_change_megajoules`, retained only for the
    /// unedited `src/ecosys_ng.zig` call site.
    sensible_energy_change_megajoules: f64,
};

const LayerCandidate = struct {
    solid_m3: f64,
    liquid_m3: f64,
    vapor_m3: f64,
    temperature_k: f64,
    heat_capacity_megajoules_per_k: f64,
    residual_m3: f64,
};

fn candidateAt(transfer_to_vapor_m3: f64, solid_m3: f64, liquid_m3: f64, vapor_m3: f64, temperature_k: f64, heat_capacity_megajoules_per_k: f64, ice_m3: f64, air_m3: f64, parameters: Parameters) !LayerCandidate {
    const liquid_change_m3 = if (transfer_to_vapor_m3 >= 0) -@min(liquid_m3, transfer_to_vapor_m3) else -transfer_to_vapor_m3;
    const solid_change_m3 = if (transfer_to_vapor_m3 > liquid_m3) -(transfer_to_vapor_m3 - liquid_m3) else 0;
    const next_solid = solid_m3 + solid_change_m3;
    const next_liquid = liquid_m3 + liquid_change_m3;
    const next_vapor = vapor_m3 + transfer_to_vapor_m3;
    if (next_solid < -1e-12 or next_liquid < -1e-12 or next_vapor < -1e-12) return error.InvalidSnowVaporCandidate;
    const latent_heat_megajoules = parameters.liquid_evaporation_latent_heat_megajoules_per_m3 * liquid_change_m3 + parameters.snow_sublimation_latent_heat_megajoules_per_m3 * solid_change_m3;
    const next_capacity = 2.095 * @max(0, next_solid) + 4.19 * (@max(0, next_liquid) + @max(0, next_vapor)) + 1.9274 * ice_m3;
    if (!std.math.isFinite(next_capacity) or next_capacity <= 0) return error.InvalidSnowVaporHeatCapacity;
    const next_temperature = (heat_capacity_megajoules_per_k * temperature_k + latent_heat_megajoules) / next_capacity;
    if (!std.math.isFinite(next_temperature) or next_temperature <= 0) return error.InvalidSnowVaporTemperature;
    const equilibrium_fraction = parameters.vapor_volume_prefactor_k / next_temperature * parameters.equilibrium_relative_humidity * std.math.exp(parameters.clausius_clapeyron_temperature_k * (parameters.reference_inverse_temperature_per_k - 1 / next_temperature));
    const residual = next_vapor - equilibrium_fraction * air_m3;
    if (!std.math.isFinite(residual)) return error.NonFiniteSnowVaporResidual;
    return .{ .solid_m3 = @max(0, next_solid), .liquid_m3 = @max(0, next_liquid), .vapor_m3 = @max(0, next_vapor), .temperature_k = next_temperature, .heat_capacity_megajoules_per_k = next_capacity, .residual_m3 = residual };
}

/// WATSUB internal-layer evaporation/condensation. Liquid water supplies a
/// vapor deficit first and solid snow supplies the remainder, exactly matching
/// WFLVW2/WFLVS2 ordering. The local equilibrium solve exits as soon as its
/// vapor residual meets tolerance and publishes atomically.
pub fn solve(allocator: std.mem.Allocator, state: *snow.State, parameters: Parameters, options: Options) !Report {
    inline for (.{ parameters.vapor_volume_prefactor_k, parameters.equilibrium_relative_humidity, parameters.clausius_clapeyron_temperature_k, parameters.reference_inverse_temperature_per_k, parameters.liquid_evaporation_latent_heat_megajoules_per_m3, parameters.snow_sublimation_latent_heat_megajoules_per_m3, options.absolute_tolerance_m3, options.relative_tolerance }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowVaporParameter;
    if (parameters.vapor_volume_prefactor_k <= 0 or parameters.equilibrium_relative_humidity < 0 or parameters.equilibrium_relative_humidity > 1 or parameters.clausius_clapeyron_temperature_k <= 0 or parameters.reference_inverse_temperature_per_k <= 0 or parameters.liquid_evaporation_latent_heat_megajoules_per_m3 <= 0 or parameters.snow_sublimation_latent_heat_megajoules_per_m3 <= 0 or options.absolute_tolerance_m3 < 0 or options.relative_tolerance < 0 or options.max_iterations == 0) return error.InvalidSnowVaporParameter;
    const solid = try allocator.dupe(f64, state.solid_snow_water_equivalent_m3);
    defer allocator.free(solid);
    const liquid = try allocator.dupe(f64, state.liquid_water_volume_m3);
    defer allocator.free(liquid);
    const vapor = try allocator.dupe(f64, state.vapor_water_equivalent_m3);
    defer allocator.free(vapor);
    const temperature = try allocator.dupe(f64, state.temperature_k);
    defer allocator.free(temperature);
    const heat_capacity = try allocator.dupe(f64, state.heat_capacity_megajoules_per_k);
    defer allocator.free(heat_capacity);
    var report: Report = .{ .iterations = 0, .converged = false, .maximum_vapor_residual_m3 = 0, .enthalpy_change_megajoules = 0, .sensible_energy_change_megajoules = 0 };
    for (0..state.cell_count) |cell| {
        // WATSUB applies this internal equilibrium only below the surface layer.
        for (1..state.layer_capacity) |layer| {
            const index = cell * state.layer_capacity + layer;
            inline for (.{ solid[index], liquid[index], vapor[index], temperature[index], heat_capacity[index], state.air_filled_volume_m3[index] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowVaporState;
            if (solid[index] < 0 or liquid[index] < 0 or vapor[index] < 0 or temperature[index] <= 0 or heat_capacity[index] < 0 or state.air_filled_volume_m3[index] < 0) return error.InvalidSnowVaporState;
            if (!state.active[index] or state.air_filled_volume_m3[index] <= 0 or heat_capacity[index] <= 0) continue;
            const lower_bound = -vapor[index];
            const upper_bound = liquid[index] + solid[index];
            var transfer_m3: f64 = 0;
            var current = try candidateAt(transfer_m3, solid[index], liquid[index], vapor[index], temperature[index], heat_capacity[index], state.ice_volume_m3[index], state.air_filled_volume_m3[index], parameters);
            var converged = false;
            for (1..@as(usize, options.max_iterations) + 1) |iteration| {
                report.iterations = @max(report.iterations, @as(u16, @intCast(iteration)));
                const tolerance_m3 = options.absolute_tolerance_m3 + options.relative_tolerance * @max(current.vapor_m3, @abs(current.vapor_m3 - current.residual_m3));
                if (@abs(current.residual_m3) <= tolerance_m3 or (transfer_m3 <= lower_bound and current.residual_m3 > 0) or (transfer_m3 >= upper_bound and current.residual_m3 < 0)) {
                    converged = true;
                    break;
                }
                const nominal_probe = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0e-9, @abs(transfer_m3));
                const probe_transfer = if (transfer_m3 + nominal_probe <= upper_bound) transfer_m3 + nominal_probe else transfer_m3 - nominal_probe;
                const probe = try candidateAt(probe_transfer, solid[index], liquid[index], vapor[index], temperature[index], heat_capacity[index], state.ice_volume_m3[index], state.air_filled_volume_m3[index], parameters);
                const derivative = (probe.residual_m3 - current.residual_m3) / (probe_transfer - transfer_m3);
                const newton_transfer = if (std.math.isFinite(derivative) and @abs(derivative) > std.math.floatEps(f64)) std.math.clamp(transfer_m3 - current.residual_m3 / derivative, lower_bound, upper_bound) else transfer_m3;
                var accepted = false;
                var fraction: f64 = 1;
                while (fraction >= 1.0e-6) : (fraction *= 0.5) {
                    const trial_transfer = std.math.clamp(transfer_m3 + fraction * (newton_transfer - transfer_m3), lower_bound, upper_bound);
                    const trial = try candidateAt(trial_transfer, solid[index], liquid[index], vapor[index], temperature[index], heat_capacity[index], state.ice_volume_m3[index], state.air_filled_volume_m3[index], parameters);
                    if (@abs(trial.residual_m3) < @abs(current.residual_m3)) {
                        transfer_m3 = trial_transfer;
                        current = trial;
                        accepted = true;
                        break;
                    }
                }
                if (!accepted) {
                    const picard_transfer = std.math.clamp(transfer_m3 - current.residual_m3, lower_bound, upper_bound);
                    if (picard_transfer == transfer_m3) break;
                    current = try candidateAt(picard_transfer, solid[index], liquid[index], vapor[index], temperature[index], heat_capacity[index], state.ice_volume_m3[index], state.air_filled_volume_m3[index], parameters);
                    transfer_m3 = picard_transfer;
                }
            }
            report.maximum_vapor_residual_m3 = @max(report.maximum_vapor_residual_m3, @abs(current.residual_m3));
            if (!converged) return error.SnowVaporSolverDidNotConverge;
            solid[index] = current.solid_m3;
            liquid[index] = current.liquid_m3;
            vapor[index] = current.vapor_m3;
            temperature[index] = current.temperature_k;
            heat_capacity[index] = current.heat_capacity_megajoules_per_k;
        }
    }
    report.converged = true;
    for (
        heat_capacity,
        temperature,
        solid,
        state.heat_capacity_megajoules_per_k,
        state.temperature_k,
        state.solid_snow_water_equivalent_m3,
    ) |
        next_capacity,
        next_temperature,
        next_solid,
        previous_capacity,
        previous_temperature,
        previous_solid,
    | {
        // `ice_volume_m3` is untouched by this solve, so only the solid snow
        // carrier contributes a latent-of-fusion term.
        report.enthalpy_change_megajoules +=
            next_capacity * next_temperature -
            previous_capacity * previous_temperature -
            parameters.latent_heat_of_fusion_megajoules_per_m3 *
                (next_solid - previous_solid);
        if (!std.math.isFinite(report.enthalpy_change_megajoules))
            return error.NonFiniteSnowVaporEnergyChange;
    }
    report.sensible_energy_change_megajoules = report.enthalpy_change_megajoules;
    @memcpy(state.solid_snow_water_equivalent_m3, solid);
    @memcpy(state.liquid_water_volume_m3, liquid);
    @memcpy(state.vapor_water_equivalent_m3, vapor);
    @memcpy(state.temperature_k, temperature);
    @memcpy(state.heat_capacity_megajoules_per_k, heat_capacity);
    state.refreshAllGeometry();
    return report;
}

test "internal vapor equilibrium conserves total water and responds below surface" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.1}, &.{1}, &.{268}, &.{ 0.05, 0.10 }, 0.1);
    state.liquid_water_volume_m3[1] = 0.001;
    state.vapor_water_equivalent_m3[1] = 0;
    state.heat_capacity_megajoules_per_k[1] += 4.19 * 0.001;
    state.refreshAllGeometry();
    const before = state.solid_snow_water_equivalent_m3[1] + state.liquid_water_volume_m3[1] + state.vapor_water_equivalent_m3[1];
    const report = try solve(std.testing.allocator, &state, .{ .vapor_volume_prefactor_k = 2.173e-3, .equilibrium_relative_humidity = 0.61, .clausius_clapeyron_temperature_k = 5360, .reference_inverse_temperature_per_k = 3.661e-3, .liquid_evaporation_latent_heat_megajoules_per_m3 = 2465, .snow_sublimation_latent_heat_megajoules_per_m3 = 2834 }, .{ .absolute_tolerance_m3 = 1e-12, .relative_tolerance = 1e-8, .max_iterations = 20 });
    try std.testing.expect(report.converged);
    try std.testing.expect(state.vapor_water_equivalent_m3[1] > 0);
    try std.testing.expectApproxEqAbs(before, state.solid_snow_water_equivalent_m3[1] + state.liquid_water_volume_m3[1] + state.vapor_water_equivalent_m3[1], 1e-14);
    try std.testing.expect(report.sensible_energy_change_megajoules < 0);
}
