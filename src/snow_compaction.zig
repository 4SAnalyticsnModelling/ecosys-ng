const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const Parameters = struct {
    maximum_temperature_metamorphism_density_megagrams_per_m3: f64,
    temperature_metamorphism_rate_per_h: f64,
    temperature_metamorphism_exponent_per_c: f64,
    viscosity_scale_megagrams_h_per_m3: f64,
    viscosity_temperature_exponent_per_c: f64,
    viscosity_density_exponent_m3_per_megagram: f64,
    minimum_snowfall_temperature_c: f64,
    maximum_snowfall_temperature_c: f64,
    snowfall_density_temperature_coefficient_megagrams_per_m3_c_pow_1_5: f64,
};

pub const Inputs = struct {
    snowfall_water_equivalent_m3: []const f64,
    atmospheric_temperature_k: []const f64,
    timestep_h: f64,
    initial_snow_density_megagrams_per_m3: f64,
    ice_density_megagrams_per_m3: f64,
};

/// REDIST snow metamorphism and overburden compaction. Density candidates are
/// validated for every runtime layer before any state is committed.
pub fn apply(allocator: std.mem.Allocator, state: *snow.State, inputs: Inputs, parameters: Parameters) !void {
    if (inputs.snowfall_water_equivalent_m3.len != state.cell_count or inputs.atmospheric_temperature_k.len != state.cell_count) return error.SnowCompactionDimensionMismatch;
    inline for (.{ inputs.timestep_h, inputs.initial_snow_density_megagrams_per_m3, inputs.ice_density_megagrams_per_m3, parameters.maximum_temperature_metamorphism_density_megagrams_per_m3, parameters.temperature_metamorphism_rate_per_h, parameters.temperature_metamorphism_exponent_per_c, parameters.viscosity_scale_megagrams_h_per_m3, parameters.viscosity_temperature_exponent_per_c, parameters.viscosity_density_exponent_m3_per_megagram, parameters.minimum_snowfall_temperature_c, parameters.maximum_snowfall_temperature_c, parameters.snowfall_density_temperature_coefficient_megagrams_per_m3_c_pow_1_5 }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowCompactionParameter;
    if (inputs.timestep_h <= 0 or inputs.initial_snow_density_megagrams_per_m3 <= 0 or inputs.ice_density_megagrams_per_m3 <= 0 or parameters.maximum_temperature_metamorphism_density_megagrams_per_m3 <= 0 or parameters.temperature_metamorphism_rate_per_h < 0 or parameters.viscosity_scale_megagrams_h_per_m3 <= 0 or parameters.maximum_snowfall_temperature_c < parameters.minimum_snowfall_temperature_c or parameters.snowfall_density_temperature_coefficient_megagrams_per_m3_c_pow_1_5 < 0) return error.InvalidSnowCompactionParameter;

    const density = try allocator.dupe(f64, state.snow_density_megagrams_per_m3);
    defer allocator.free(density);
    for (0..state.cell_count) |cell| {
        const snowfall = inputs.snowfall_water_equivalent_m3[cell];
        const atmospheric_temperature_k = inputs.atmospheric_temperature_k[cell];
        if (!std.math.isFinite(snowfall) or snowfall < 0 or !std.math.isFinite(atmospheric_temperature_k) or atmospheric_temperature_k <= 0) return error.InvalidSnowCompactionInput;
        const area_m2 = state.horizontal_area_m2[cell * state.layer_capacity];
        if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidSnowCompactionArea;
        var overburden_water_equivalent_m3: f64 = 0;
        for (0..state.layer_capacity) |layer| {
            const index = cell * state.layer_capacity + layer;
            inline for (.{ state.solid_snow_water_equivalent_m3[index], state.liquid_water_volume_m3[index], state.ice_volume_m3[index], state.temperature_k[index], density[index] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowCompactionState;
            if (state.solid_snow_water_equivalent_m3[index] < 0 or state.liquid_water_volume_m3[index] < 0 or state.ice_volume_m3[index] < 0 or state.temperature_k[index] <= 0 or density[index] <= 0) return error.InvalidSnowCompactionState;

            const layer_water_equivalent_m3 = state.solid_snow_water_equivalent_m3[index] + state.liquid_water_volume_m3[index] + state.ice_volume_m3[index] * inputs.ice_density_megagrams_per_m3;
            overburden_water_equivalent_m3 += 0.5 * layer_water_equivalent_m3;
            // REDIST enters this block only for an extant snow layer. Keep the
            // inherited density of empty/runtime-capacity layers unchanged.
            if (!state.active[index] or state.solid_snow_water_equivalent_m3[index] <= 0) {
                overburden_water_equivalent_m3 += 0.5 * layer_water_equivalent_m3;
                continue;
            }
            if (layer == 0 and snowfall > 0 and state.solid_snow_water_equivalent_m3[index] > 0) {
                const snowfall_temperature_c = std.math.clamp(atmospheric_temperature_k - 273.15, parameters.minimum_snowfall_temperature_c, parameters.maximum_snowfall_temperature_c);
                const snowfall_density_megagrams_per_m3 = inputs.initial_snow_density_megagrams_per_m3 + parameters.snowfall_density_temperature_coefficient_megagrams_per_m3_c_pow_1_5 * std.math.pow(f64, snowfall_temperature_c - parameters.minimum_snowfall_temperature_c, 1.5);
                const prior_solid = @max(0, state.solid_snow_water_equivalent_m3[index] - snowfall);
                const mixed_solid_volume_m3 = snowfall / snowfall_density_megagrams_per_m3 + prior_solid / density[index];
                if (mixed_solid_volume_m3 > 0) density[index] = state.solid_snow_water_equivalent_m3[index] / mixed_solid_volume_m3;
            }

            const temperature_c = state.temperature_k[index] - 273.15;
            const temperature_rate = if (density[index] < parameters.maximum_temperature_metamorphism_density_megagrams_per_m3)
                density[index] * parameters.temperature_metamorphism_rate_per_h * std.math.exp(parameters.temperature_metamorphism_exponent_per_c * temperature_c)
            else
                0;
            const viscosity_megagrams_h_per_m3 = parameters.viscosity_scale_megagrams_h_per_m3 * std.math.exp(parameters.viscosity_temperature_exponent_per_c * temperature_c + parameters.viscosity_density_exponent_m3_per_megagram * density[index]);
            const overburden_rate = density[index] * overburden_water_equivalent_m3 / (area_m2 * viscosity_megagrams_h_per_m3);
            const next_density = density[index] + (temperature_rate + overburden_rate) * inputs.timestep_h;
            if (!std.math.isFinite(next_density) or next_density <= 0) return error.InvalidSnowCompactionResult;
            density[index] = next_density;
            overburden_water_equivalent_m3 += 0.5 * layer_water_equivalent_m3;
        }
    }
    @memcpy(state.snow_density_megagrams_per_m3, density);
    state.refreshAllGeometry();
}

test "REDIST compaction increases density and conserves every inventory" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.1}, &.{2}, &.{268}, &.{ 0.05, 0.15 }, 0.05);
    state.liquid_water_volume_m3[0] = 0.001;
    state.amount_g[0] = 7;
    const solid_before = state.solid_snow_water_equivalent_m3[0] + state.solid_snow_water_equivalent_m3[1];
    try apply(std.testing.allocator, &state, .{ .snowfall_water_equivalent_m3 = &.{0}, .atmospheric_temperature_k = &.{268}, .timestep_h = 1, .initial_snow_density_megagrams_per_m3 = 0.05, .ice_density_megagrams_per_m3 = 0.92 }, .{ .maximum_temperature_metamorphism_density_megagrams_per_m3 = 0.25, .temperature_metamorphism_rate_per_h = 1e-5, .temperature_metamorphism_exponent_per_c = 0.04, .viscosity_scale_megagrams_h_per_m3 = 0.25, .viscosity_temperature_exponent_per_c = -0.08, .viscosity_density_exponent_m3_per_megagram = 23, .minimum_snowfall_temperature_c = -15, .maximum_snowfall_temperature_c = 2, .snowfall_density_temperature_coefficient_megagrams_per_m3_c_pow_1_5 = 1.7e-3 });
    try std.testing.expect(state.snow_density_megagrams_per_m3[0] > 0.05);
    try std.testing.expect(state.snow_density_megagrams_per_m3[1] > 0.05);
    try std.testing.expectApproxEqAbs(solid_before, state.solid_snow_water_equivalent_m3[0] + state.solid_snow_water_equivalent_m3[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 7), state.amount_g[0]);
}

test "empty capacity layers do not accumulate fictitious density" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.02}, &.{1}, &.{268}, &.{ 0.05, 0.10 }, 0.05);
    const empty_density = state.snow_density_megagrams_per_m3[1];
    try apply(std.testing.allocator, &state, .{ .snowfall_water_equivalent_m3 = &.{0}, .atmospheric_temperature_k = &.{268}, .timestep_h = 1, .initial_snow_density_megagrams_per_m3 = 0.05, .ice_density_megagrams_per_m3 = 0.92 }, .{ .maximum_temperature_metamorphism_density_megagrams_per_m3 = 0.25, .temperature_metamorphism_rate_per_h = 1e-5, .temperature_metamorphism_exponent_per_c = 0.04, .viscosity_scale_megagrams_h_per_m3 = 0.25, .viscosity_temperature_exponent_per_c = -0.08, .viscosity_density_exponent_m3_per_megagram = 23, .minimum_snowfall_temperature_c = -15, .maximum_snowfall_temperature_c = 2, .snowfall_density_temperature_coefficient_megagrams_per_m3_c_pow_1_5 = 1.7e-3 });
    try std.testing.expectEqual(empty_density, state.snow_density_megagrams_per_m3[1]);
}
