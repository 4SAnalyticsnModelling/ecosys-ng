const std = @import("std");

pub const Inputs = struct {
    standing_dead_surface_area_m2: f64,
    standing_dead_presence_threshold_m2: f64,
    standing_dead_air_temperature_k: f64,
    standing_dead_air_vapor_volume_fraction: f64,
    precipitation_retention_rate_m3_per_h: f64,
    timestep_h: f64,
    sapwood_thickness_m: f64,
    standing_dead_mass_g_c: f64,
    stalk_volume_m3_per_g_c: f64,
    surface_water_m3: f64,
    minimum_water_solver_heat_capacity_megajoules_per_m2_k: f64,
    minimum_energy_solver_heat_capacity_megajoules_per_m2_k: f64,
    cell_area_m2: f64,
    surface_temperature_k: f64,
    large_temperature_increment_k: f64,
    small_temperature_increment_k: f64,
    radiation_fraction: f64,
    radiation_presence_threshold: f64,
    absorbed_shortwave_rate_megajoules_per_h: f64,
    emissivity: f64,
    sky_longwave_rate_megajoules_per_h: f64,
    lateral_longwave_rate_megajoules_per_h: f64,
    combined_area_radiation_fraction: f64,
    latent_conductance_m2_per_h: f64,
    sensible_conductance_megajoules_per_m_k_h: f64,
};

pub const Diagnostics = struct {
    net_radiation_megajoules_per_step: f64,
    latent_heat_megajoules_per_step: f64,
    sensible_heat_megajoules_per_step: f64,
    heat_storage_megajoules_per_step: f64,
    vapor_convective_heat_megajoules_per_step: f64,
    emitted_thermal_radiation_megajoules_per_step: f64,
    evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    ground_canopy_vapor_flux_m3_per_step: f64,
    ground_canopy_sensible_heat_megajoules_per_step: f64,
};

pub const Result = struct {
    standing_dead_present: bool,
    energy_balance_admitted: bool,
    air_temperature_k: f64,
    air_vapor_volume_fraction: f64,
    retained_precipitation_m3_per_step: f64,
    dry_heat_capacity_megajoules_per_k: f64,
    wet_heat_capacity_megajoules_per_k: f64,
    minimum_water_solver_heat_capacity_megajoules_per_k: f64,
    minimum_energy_solver_heat_capacity_megajoules_per_k: f64,
    surface_temperature_k: f64,
    previous_surface_temperature_k: f64,
    temperature_increment_k: f64,
    absorbed_shortwave_megajoules_per_step: f64,
    emitted_longwave_coefficient_megajoules_per_step_k4: f64,
    absorbed_sky_longwave_megajoules_per_step: f64,
    absorbed_lateral_longwave_megajoules_per_step: f64,
    latent_conductance_m2_per_step: f64,
    sensible_conductance_megajoules_per_m_k_step: f64,
    diagnostics: Diagnostics,
};

/// UPTAKE.F 3895--3975 standing-dead canopy and energy-balance admission.
pub fn compute(inputs: Inputs) !Result {
    try validatePresenceGate(inputs);
    var result = std.mem.zeroes(Result);
    if (inputs.standing_dead_surface_area_m2 <=
        inputs.standing_dead_presence_threshold_m2)
        return result;

    try validateActiveInputs(inputs);
    result.standing_dead_present = true;
    result.air_temperature_k = inputs.standing_dead_air_temperature_k;
    result.air_vapor_volume_fraction =
        inputs.standing_dead_air_vapor_volume_fraction;
    result.retained_precipitation_m3_per_step =
        inputs.precipitation_retention_rate_m3_per_h * inputs.timestep_h;
    result.dry_heat_capacity_megajoules_per_k = 2.496 * @min(
        inputs.sapwood_thickness_m * inputs.standing_dead_surface_area_m2,
        inputs.standing_dead_mass_g_c * inputs.stalk_volume_m3_per_g_c,
    );
    result.wet_heat_capacity_megajoules_per_k =
        result.dry_heat_capacity_megajoules_per_k +
        4.19 * @max(0, inputs.surface_water_m3);
    result.minimum_water_solver_heat_capacity_megajoules_per_k =
        inputs.minimum_water_solver_heat_capacity_megajoules_per_m2_k * inputs.cell_area_m2;
    result.minimum_energy_solver_heat_capacity_megajoules_per_k =
        inputs.minimum_energy_solver_heat_capacity_megajoules_per_m2_k * inputs.cell_area_m2;
    result.surface_temperature_k = inputs.surface_temperature_k;
    result.previous_surface_temperature_k = result.surface_temperature_k;
    result.temperature_increment_k =
        if (result.dry_heat_capacity_megajoules_per_k >
        result.minimum_water_solver_heat_capacity_megajoules_per_k)
            inputs.large_temperature_increment_k
        else
            inputs.small_temperature_increment_k;

    if (result.wet_heat_capacity_megajoules_per_k <=
        result.minimum_energy_solver_heat_capacity_megajoules_per_k or
        inputs.radiation_fraction <= inputs.radiation_presence_threshold)
        return result;

    result.energy_balance_admitted = true;
    result.absorbed_shortwave_megajoules_per_step =
        inputs.absorbed_shortwave_rate_megajoules_per_h * inputs.timestep_h;
    result.emitted_longwave_coefficient_megajoules_per_step_k4 =
        inputs.emissivity * 2.04e-10 * inputs.radiation_fraction *
        inputs.cell_area_m2 * inputs.timestep_h;
    result.absorbed_sky_longwave_megajoules_per_step =
        inputs.sky_longwave_rate_megajoules_per_h * inputs.radiation_fraction *
        inputs.timestep_h;
    result.absorbed_lateral_longwave_megajoules_per_step =
        inputs.lateral_longwave_rate_megajoules_per_h *
        inputs.combined_area_radiation_fraction * inputs.timestep_h;
    result.diagnostics = std.mem.zeroes(Diagnostics);
    result.latent_conductance_m2_per_step =
        inputs.radiation_fraction * inputs.latent_conductance_m2_per_h;
    result.sensible_conductance_megajoules_per_m_k_step =
        inputs.radiation_fraction * inputs.sensible_conductance_megajoules_per_m_k_h;
    try validateResult(result);
    return result;
}

pub fn computeRuntimeSpecies(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.StandingDeadCanopyAdmissionDimensionMismatch;
    for (inputs, scratch) |species_inputs, *candidate|
        candidate.* = try compute(species_inputs);
    @memcpy(destination, scratch);
}

fn validatePresenceGate(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.standing_dead_surface_area_m2) or
        !std.math.isFinite(inputs.standing_dead_presence_threshold_m2))
        return error.NonFiniteStandingDeadCanopyAdmissionInput;
    if (inputs.standing_dead_surface_area_m2 < 0 or
        inputs.standing_dead_presence_threshold_m2 < 0)
        return error.InvalidStandingDeadCanopyAdmissionInput;
}

fn validateActiveInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteStandingDeadCanopyAdmissionInput;
    if (inputs.standing_dead_air_temperature_k <= 0 or
        inputs.standing_dead_air_vapor_volume_fraction < 0 or
        inputs.precipitation_retention_rate_m3_per_h < 0 or
        inputs.timestep_h < 0 or inputs.sapwood_thickness_m < 0 or
        inputs.standing_dead_mass_g_c < 0 or
        inputs.stalk_volume_m3_per_g_c < 0 or inputs.surface_water_m3 < 0 or
        inputs.minimum_water_solver_heat_capacity_megajoules_per_m2_k < 0 or
        inputs.minimum_energy_solver_heat_capacity_megajoules_per_m2_k < 0 or
        inputs.cell_area_m2 <= 0 or inputs.surface_temperature_k <= 0 or
        inputs.large_temperature_increment_k <= 0 or
        inputs.small_temperature_increment_k <= 0 or
        inputs.radiation_fraction < 0 or inputs.radiation_fraction > 1 or
        inputs.radiation_presence_threshold < 0 or
        inputs.absorbed_shortwave_rate_megajoules_per_h < 0 or inputs.emissivity < 0 or
        inputs.sky_longwave_rate_megajoules_per_h < 0 or
        inputs.combined_area_radiation_fraction < 0 or
        inputs.combined_area_radiation_fraction > 1 or
        inputs.latent_conductance_m2_per_h < 0 or
        inputs.sensible_conductance_megajoules_per_m_k_h < 0)
        return error.InvalidStandingDeadCanopyAdmissionInput;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteStandingDeadCanopyAdmissionResult;
    }
}

fn testInputs() Inputs {
    return .{
        .standing_dead_surface_area_m2 = 2,
        .standing_dead_presence_threshold_m2 = 1e-12,
        .standing_dead_air_temperature_k = 280,
        .standing_dead_air_vapor_volume_fraction = 0.01,
        .precipitation_retention_rate_m3_per_h = 0.4,
        .timestep_h = 0.25,
        .sapwood_thickness_m = 0.1,
        .standing_dead_mass_g_c = 100,
        .stalk_volume_m3_per_g_c = 0.001,
        .surface_water_m3 = 0.02,
        .minimum_water_solver_heat_capacity_megajoules_per_m2_k = 0.01,
        .minimum_energy_solver_heat_capacity_megajoules_per_m2_k = 0.01,
        .cell_area_m2 = 10,
        .surface_temperature_k = 275,
        .large_temperature_increment_k = 0.1,
        .small_temperature_increment_k = 1,
        .radiation_fraction = 0.5,
        .radiation_presence_threshold = 1e-6,
        .absorbed_shortwave_rate_megajoules_per_h = 4,
        .emissivity = 0.96,
        .sky_longwave_rate_megajoules_per_h = 2,
        .lateral_longwave_rate_megajoules_per_h = -0.5,
        .combined_area_radiation_fraction = 0.3,
        .latent_conductance_m2_per_h = 8,
        .sensible_conductance_megajoules_per_m_k_h = 6,
    };
}

test "standing-dead admission preserves source initialization order" {
    const result = try compute(testInputs());
    try std.testing.expect(result.standing_dead_present);
    try std.testing.expect(result.energy_balance_admitted);
    try std.testing.expectEqual(@as(f64, 0.1), result.retained_precipitation_m3_per_step);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2496), result.dry_heat_capacity_megajoules_per_k, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3334), result.wet_heat_capacity_megajoules_per_k, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.1), result.temperature_increment_k);
    try std.testing.expectEqual(@as(f64, 1), result.absorbed_shortwave_megajoules_per_step);
    try std.testing.expectEqualDeep(std.mem.zeroes(Diagnostics), result.diagnostics);
    try std.testing.expectEqual(@as(f64, 4), result.latent_conductance_m2_per_step);
}

test "presence threshold skips all downstream input use" {
    var inputs = testInputs();
    inputs.standing_dead_surface_area_m2 = inputs.standing_dead_presence_threshold_m2;
    inputs.surface_temperature_k = std.math.nan(f64);
    try std.testing.expectEqualDeep(std.mem.zeroes(Result), try compute(inputs));
}

test "low wet heat capacity rejects only energy-balance admission" {
    var inputs = testInputs();
    inputs.minimum_energy_solver_heat_capacity_megajoules_per_m2_k = 1;
    const result = try compute(inputs);
    try std.testing.expect(result.standing_dead_present);
    try std.testing.expect(!result.energy_balance_admitted);
    try std.testing.expectEqual(@as(f64, 275), result.surface_temperature_k);
    try std.testing.expectEqual(@as(f64, 0), result.absorbed_shortwave_megajoules_per_step);
}

test "dry heat-capacity threshold selects small temperature increment" {
    var inputs = testInputs();
    inputs.minimum_water_solver_heat_capacity_megajoules_per_m2_k = 1;
    const result = try compute(inputs);
    try std.testing.expectEqual(inputs.small_temperature_increment_k, result.temperature_increment_k);
}

test "later invalid runtime species leaves destination unchanged" {
    var invalid = testInputs();
    invalid.cell_area_m2 = 0;
    const inputs = [_]Inputs{ testInputs(), invalid };
    var scratch: [2]Result = undefined;
    var destination = [_]Result{ std.mem.zeroes(Result), std.mem.zeroes(Result) };
    destination[0].air_temperature_k = 41;
    destination[1].air_temperature_k = 42;
    try std.testing.expectError(
        error.InvalidStandingDeadCanopyAdmissionInput,
        computeRuntimeSpecies(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].air_temperature_k);
    try std.testing.expectEqual(@as(f64, 42), destination[1].air_temperature_k);
}
