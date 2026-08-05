const std = @import("std");

/// Runtime hourly ceilings corresponding to the scientific controls used by
/// SOLUTE. Units follow the downstream source declarations at solute.f:942,
/// 964, 1443, and 1902 rather than the dimensionally incomplete h-1 comment
/// at solute.f:98.
pub const HourlyCeilings = struct {
    phosphate_precipitation_mol_per_m3_h: f64,
    hydroxide_mineral_precipitation_mol_per_m3_h: f64,
    anion_adsorption_mol_per_m3_h: f64,
    cation_adsorption_mol_per_m3_h: f64,
    general_solute_reaction_mol_per_m3_h: f64,
    fast_solute_reaction_mol_per_m3_h: f64,
    silicate_weathering_mol_per_m2_h: f64,
};

pub const TemperatureResponseParameters = struct {
    gas_constant_j_per_mol_k: f64,
    entropy_temperature_scale_j_per_mol_k: f64,
    arrhenius_log_prefactor: f64,
    activation_energy_j_per_mol: f64,
    low_temperature_inactivation_j_per_mol: f64,
    high_temperature_inactivation_j_per_mol: f64,
};

pub const GroundSilicateConcentrations = struct {
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
};

pub const Inputs = struct {
    timestep_h: f64,
    maximum_equilibrium_iterations: usize,
    soil_temperature_k: f64,
    matrix_water_volume_m3: f64,
    net_hydrogen_source_mol_per_step: f64,
    natural_silicate_surface_area_m2_per_m3: f64,
    ground_silicate_specific_surface_area_m2_per_mol: f64,
    ground_silicate: GroundSilicateConcentrations,
    hourly: HourlyCeilings,
    temperature_response: TemperatureResponseParameters,
};

pub const Controls = struct {
    maximum_phosphate_precipitation_mol_per_m3_iteration: f64,
    /// Legacy TPDA. In the full dynamic-salt branch SOLUTE.F line 636 aliases
    /// the hydroxyapatite ceiling to TPDX rather than deriving a second rate.
    maximum_apatite_precipitation_mol_per_m3_iteration: f64,
    maximum_hydroxide_mineral_precipitation_mol_per_m3_iteration: f64,
    maximum_anion_adsorption_mol_per_m3_iteration: f64,
    maximum_cation_adsorption_mol_per_m3_iteration: f64,
    maximum_general_solute_reaction_mol_per_m3_iteration: f64,
    maximum_fast_solute_reaction_mol_per_m3_iteration: f64,
    maximum_natural_weathering_mol_per_m3_iteration: f64,
    maximum_ground_weathering_mol_per_m3_iteration: f64,
    weathering_temperature_response: f64,
    ground_silicate_surface_area_m2_per_m3: f64,
    hydrogen_source_mol_per_m3_iteration: f64,
};

/// Per-step ceilings used when the site disables the complete salt network.
/// SOLUTE.F assigns the source `TRWH` control to both weathering (`TRW`) and
/// apatite precipitation (`TPA`). That shared source control is retained
/// explicitly because the source comment does not establish two independent
/// dimensions for it.
pub const RestrictedSaltControls = struct {
    maximum_phosphate_precipitation_mol_per_m3_step: f64,
    maximum_hydroxide_mineral_precipitation_mol_per_m3_step: f64,
    maximum_anion_adsorption_mol_per_m3_step: f64,
    maximum_cation_adsorption_mol_per_m3_step: f64,
    maximum_general_solute_reaction_mol_per_m3_step: f64,
    maximum_fast_solute_reaction_mol_per_m3_step: f64,
    shared_weathering_apatite_control_per_step: f64,
};

/// Compatibility aliases initialized from the single legacy `FIONZ`
/// substrate/product limiter. Units are fraction per reaction iteration.
pub const SubstrateProductLimits = struct {
    /// Legacy FIONN, consumed by soil reactions beginning at SOLUTE.F 2991.
    soil_fraction_per_iteration: f64,
    /// Legacy FION0, consumed by surface-litter reactions beginning at
    /// SOLUTE.F 4423.
    surface_litter_fraction_per_iteration: f64,
};

/// Direct source-order translation of SOLUTE.F lines 141--142:
/// `FIONN = FIONZ`, followed by `FION0 = FIONZ`.
///
/// Production configurations may expose process-specific controls as an
/// intentional formulation extension. Compatibility initialization must use
/// this function so both domains retain the source alias relationship.
pub fn initializeSubstrateProductLimits(
    general_fraction_per_iteration: f64,
) !SubstrateProductLimits {
    if (!std.math.isFinite(general_fraction_per_iteration) or
        general_fraction_per_iteration < 0 or
        general_fraction_per_iteration > 1)
        return error.InvalidSoluteSubstrateProductLimit;

    var result: SubstrateProductLimits = undefined;
    result.soil_fraction_per_iteration = general_fraction_per_iteration;
    result.surface_litter_fraction_per_iteration =
        general_fraction_per_iteration;
    return result;
}

/// Direct source-order translation of SOLUTE.F lines 2896--2903 and the
/// identical surface-litter assignments at lines 4009--4016.
///
/// Unlike the complete salt network, this branch performs one restricted
/// reaction pass and therefore does not divide the hourly ceilings by the
/// nonlinear iteration limit.
pub fn deriveRestrictedSalt(
    hourly: HourlyCeilings,
    timestep_h: f64,
) !RestrictedSaltControls {
    if (!std.math.isFinite(timestep_h) or timestep_h <= 0)
        return error.InvalidSoluteKineticControlInput;
    inline for (@typeInfo(HourlyCeilings).@"struct".fields) |field| {
        const value = @field(hourly, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoluteKineticControlInput;
    }

    // Preserve TPD, TPZ, TADA, TADC, TSL, TSZ, TRW, then TPA = TRW.
    const phosphate = hourly.phosphate_precipitation_mol_per_m3_h * timestep_h;
    const hydroxide_mineral =
        hourly.hydroxide_mineral_precipitation_mol_per_m3_h * timestep_h;
    const anion_adsorption = hourly.anion_adsorption_mol_per_m3_h * timestep_h;
    const cation_adsorption =
        hourly.cation_adsorption_mol_per_m3_h * timestep_h;
    const general_reaction =
        hourly.general_solute_reaction_mol_per_m3_h * timestep_h;
    const fast_reaction =
        hourly.fast_solute_reaction_mol_per_m3_h * timestep_h;
    const weathering =
        hourly.silicate_weathering_mol_per_m2_h * timestep_h;
    const result = RestrictedSaltControls{
        .maximum_phosphate_precipitation_mol_per_m3_step = phosphate,
        .maximum_hydroxide_mineral_precipitation_mol_per_m3_step = hydroxide_mineral,
        .maximum_anion_adsorption_mol_per_m3_step = anion_adsorption,
        .maximum_cation_adsorption_mol_per_m3_step = cation_adsorption,
        .maximum_general_solute_reaction_mol_per_m3_step = general_reaction,
        .maximum_fast_solute_reaction_mol_per_m3_step = fast_reaction,
        .shared_weathering_apatite_control_per_step = weathering,
    };
    inline for (@typeInfo(RestrictedSaltControls).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSoluteKineticControl;
    return result;
}

/// Direct compatibility translation of SOLUTE lines 629-642 and 816-818.
/// The hybrid production solver must not use its iteration count as a
/// scientific rate parameter without an explicit compatibility-mode choice.
pub fn derive(inputs: Inputs) !Controls {
    try validate(inputs);
    const iteration_count: f64 =
        @floatFromInt(inputs.maximum_equilibrium_iterations);

    // Preserve the source assignment order: hourly coefficient / MRXN, then
    // multiplication by XNFH.
    const phosphate_precipitation =
        inputs.hourly.phosphate_precipitation_mol_per_m3_h /
        iteration_count * inputs.timestep_h;
    const hydroxide_mineral_precipitation =
        inputs.hourly.hydroxide_mineral_precipitation_mol_per_m3_h /
        iteration_count * inputs.timestep_h;
    const anion_adsorption =
        inputs.hourly.anion_adsorption_mol_per_m3_h /
        iteration_count * inputs.timestep_h;
    const cation_adsorption =
        inputs.hourly.cation_adsorption_mol_per_m3_h /
        iteration_count * inputs.timestep_h;
    const general_solute_reaction =
        inputs.hourly.general_solute_reaction_mol_per_m3_h /
        iteration_count * inputs.timestep_h;
    const fast_solute_reaction =
        inputs.hourly.fast_solute_reaction_mol_per_m3_h /
        iteration_count * inputs.timestep_h;
    const weathering_per_area =
        inputs.hourly.silicate_weathering_mol_per_m2_h /
        iteration_count * inputs.timestep_h;
    // SOLUTE.F line 636 follows all seven timestep-scaled rate assignments.
    const apatite_precipitation = phosphate_precipitation;

    const temperature = inputs.temperature_response;
    const gas_energy_j_per_mol =
        temperature.gas_constant_j_per_mol_k * inputs.soil_temperature_k;
    const entropy_energy_j_per_mol =
        temperature.entropy_temperature_scale_j_per_mol_k *
        inputs.soil_temperature_k;
    const inactivation = 1 +
        @exp((temperature.low_temperature_inactivation_j_per_mol -
            entropy_energy_j_per_mol) / gas_energy_j_per_mol) +
        @exp((entropy_energy_j_per_mol -
            temperature.high_temperature_inactivation_j_per_mol) /
            gas_energy_j_per_mol);
    const weathering_temperature_response =
        @exp(temperature.arrhenius_log_prefactor -
            temperature.activation_energy_j_per_mol /
                gas_energy_j_per_mol) /
        inactivation;

    const ground = inputs.ground_silicate;
    const total_ground_silicate_mol_per_m3 =
        ground.aluminum_mol_per_m3 +
        ground.iron_mol_per_m3 +
        ground.calcium_mol_per_m3 +
        ground.magnesium_mol_per_m3 +
        ground.sodium_mol_per_m3 +
        ground.potassium_mol_per_m3;
    const ground_surface_area_m2_per_m3 =
        inputs.ground_silicate_specific_surface_area_m2_per_mol *
        total_ground_silicate_mol_per_m3;

    const result = Controls{
        .maximum_phosphate_precipitation_mol_per_m3_iteration = phosphate_precipitation,
        .maximum_apatite_precipitation_mol_per_m3_iteration = apatite_precipitation,
        .maximum_hydroxide_mineral_precipitation_mol_per_m3_iteration = hydroxide_mineral_precipitation,
        .maximum_anion_adsorption_mol_per_m3_iteration = anion_adsorption,
        .maximum_cation_adsorption_mol_per_m3_iteration = cation_adsorption,
        .maximum_general_solute_reaction_mol_per_m3_iteration = general_solute_reaction,
        .maximum_fast_solute_reaction_mol_per_m3_iteration = fast_solute_reaction,
        .maximum_natural_weathering_mol_per_m3_iteration = weathering_per_area *
            inputs.natural_silicate_surface_area_m2_per_m3 *
            weathering_temperature_response,
        .maximum_ground_weathering_mol_per_m3_iteration = weathering_per_area * ground_surface_area_m2_per_m3 *
            weathering_temperature_response,
        .weathering_temperature_response = weathering_temperature_response,
        .ground_silicate_surface_area_m2_per_m3 = ground_surface_area_m2_per_m3,
        .hydrogen_source_mol_per_m3_iteration = inputs.net_hydrogen_source_mol_per_step /
            (iteration_count * inputs.matrix_water_volume_m3),
    };
    inline for (@typeInfo(Controls).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSoluteKineticControl;
    }
    return result;
}

fn validate(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.timestep_h) or inputs.timestep_h <= 0 or
        inputs.maximum_equilibrium_iterations == 0 or
        !std.math.isFinite(inputs.soil_temperature_k) or
        inputs.soil_temperature_k <= 0 or
        !std.math.isFinite(inputs.matrix_water_volume_m3) or
        inputs.matrix_water_volume_m3 <= 0 or
        !std.math.isFinite(inputs.net_hydrogen_source_mol_per_step) or
        !std.math.isFinite(inputs.natural_silicate_surface_area_m2_per_m3) or
        inputs.natural_silicate_surface_area_m2_per_m3 < 0 or
        !std.math.isFinite(inputs.ground_silicate_specific_surface_area_m2_per_mol) or
        inputs.ground_silicate_specific_surface_area_m2_per_mol < 0)
    {
        return error.InvalidSoluteKineticControlInput;
    }
    inline for (@typeInfo(HourlyCeilings).@"struct".fields) |field| {
        const value = @field(inputs.hourly, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoluteKineticControlInput;
    }
    inline for (@typeInfo(GroundSilicateConcentrations).@"struct".fields) |field| {
        const value = @field(inputs.ground_silicate, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoluteKineticControlInput;
    }
    inline for (@typeInfo(TemperatureResponseParameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.temperature_response, field.name)))
            return error.InvalidSoluteKineticControlInput;
    }
    if (inputs.temperature_response.gas_constant_j_per_mol_k <= 0 or
        inputs.temperature_response.entropy_temperature_scale_j_per_mol_k <= 0 or
        inputs.temperature_response.activation_energy_j_per_mol < 0)
    {
        return error.InvalidSoluteKineticControlInput;
    }
}

fn sourceInputs() Inputs {
    return .{
        .timestep_h = 1,
        .maximum_equilibrium_iterations = 60,
        .soil_temperature_k = 298.15,
        .matrix_water_volume_m3 = 2,
        .net_hydrogen_source_mol_per_step = 6,
        .natural_silicate_surface_area_m2_per_m3 = 12,
        .ground_silicate_specific_surface_area_m2_per_mol = 1000,
        .ground_silicate = .{
            .aluminum_mol_per_m3 = 1,
            .iron_mol_per_m3 = 2,
            .calcium_mol_per_m3 = 3,
            .magnesium_mol_per_m3 = 4,
            .sodium_mol_per_m3 = 5,
            .potassium_mol_per_m3 = 6,
        },
        .hourly = .{
            .phosphate_precipitation_mol_per_m3_h = 2.5e-3,
            .hydroxide_mineral_precipitation_mol_per_m3_h = 2.5e-2,
            .anion_adsorption_mol_per_m3_h = 2.5e-3,
            .cation_adsorption_mol_per_m3_h = 2.5e-3,
            .general_solute_reaction_mol_per_m3_h = 2.5e-2,
            .fast_solute_reaction_mol_per_m3_h = 2.5e-1,
            .silicate_weathering_mol_per_m2_h = 5.0e-5,
        },
        .temperature_response = .{
            .gas_constant_j_per_mol_k = 8.3143,
            .entropy_temperature_scale_j_per_mol_k = 710,
            .arrhenius_log_prefactor = 25.229,
            .activation_energy_j_per_mol = 62_500,
            .low_temperature_inactivation_j_per_mol = 197_500,
            .high_temperature_inactivation_j_per_mol = 222_500,
        },
    };
}

test "SOLUTE kinetic controls preserve source assignment order" {
    const inputs = sourceInputs();
    const result = try derive(inputs);
    const iteration_count: f64 =
        @floatFromInt(inputs.maximum_equilibrium_iterations);
    const rt = 8.3143 * inputs.soil_temperature_k;
    const st = 710 * inputs.soil_temperature_k;
    const response = @exp(25.229 - 62_500 / rt) /
        (1 + @exp((197_500 - st) / rt) +
            @exp((st - 222_500) / rt));
    const weathering_per_area = 5.0e-5 / iteration_count;

    try std.testing.expectEqual(
        2.5e-3 / iteration_count,
        result.maximum_phosphate_precipitation_mol_per_m3_iteration,
    );
    try std.testing.expectEqual(
        result.maximum_phosphate_precipitation_mol_per_m3_iteration,
        result.maximum_apatite_precipitation_mol_per_m3_iteration,
    );
    try std.testing.expectEqual(response, result.weathering_temperature_response);
    try std.testing.expectEqual(
        21_000,
        result.ground_silicate_surface_area_m2_per_m3,
    );
    try std.testing.expectEqual(
        weathering_per_area * 12 * response,
        result.maximum_natural_weathering_mol_per_m3_iteration,
    );
    try std.testing.expectEqual(
        weathering_per_area * 21_000 * response,
        result.maximum_ground_weathering_mol_per_m3_iteration,
    );
    try std.testing.expectEqual(
        6.0 / (60.0 * 2.0),
        result.hydrogen_source_mol_per_m3_iteration,
    );
}

test "SOLUTE kinetic controls use runtime timestep and iteration count" {
    var inputs = sourceInputs();
    inputs.timestep_h = 0.5;
    inputs.maximum_equilibrium_iterations = 120;
    const result = try derive(inputs);
    try std.testing.expectEqual(
        2.5e-1 / 120.0 * 0.5,
        result.maximum_fast_solute_reaction_mol_per_m3_iteration,
    );
    try std.testing.expectEqual(
        6.0 / (120.0 * 2.0),
        result.hydrogen_source_mol_per_m3_iteration,
    );
}

test "SOLUTE kinetic controls reject invalid physical domains" {
    var inputs = sourceInputs();
    inputs.maximum_equilibrium_iterations = 0;
    try std.testing.expectError(
        error.InvalidSoluteKineticControlInput,
        derive(inputs),
    );
    inputs = sourceInputs();
    inputs.matrix_water_volume_m3 = 0;
    try std.testing.expectError(
        error.InvalidSoluteKineticControlInput,
        derive(inputs),
    );
    inputs = sourceInputs();
    inputs.ground_silicate.iron_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSoluteKineticControlInput,
        derive(inputs),
    );
}

test "restricted-salt controls preserve source timestep scaling and alias" {
    const result = try deriveRestrictedSalt(.{
        .phosphate_precipitation_mol_per_m3_h = 1,
        .hydroxide_mineral_precipitation_mol_per_m3_h = 2,
        .anion_adsorption_mol_per_m3_h = 3,
        .cation_adsorption_mol_per_m3_h = 4,
        .general_solute_reaction_mol_per_m3_h = 5,
        .fast_solute_reaction_mol_per_m3_h = 6,
        .silicate_weathering_mol_per_m2_h = 7,
    }, 0.25);
    try std.testing.expectEqual(
        @as(f64, 0.25),
        result.maximum_phosphate_precipitation_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.maximum_hydroxide_mineral_precipitation_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.75),
        result.maximum_anion_adsorption_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 1),
        result.maximum_cation_adsorption_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 1.25),
        result.maximum_general_solute_reaction_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 1.5),
        result.maximum_fast_solute_reaction_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 1.75),
        result.shared_weathering_apatite_control_per_step,
    );
}

test "restricted-salt controls reject invalid runtime controls" {
    var hourly = sourceInputs().hourly;
    try std.testing.expectError(
        error.InvalidSoluteKineticControlInput,
        deriveRestrictedSalt(hourly, 0),
    );
    hourly.fast_solute_reaction_mol_per_m3_h = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSoluteKineticControlInput,
        deriveRestrictedSalt(hourly, 1),
    );
}

test "SOLUTE lines 141 and 142 preserve FIONZ alias and assignment order" {
    const limits = try initializeSubstrateProductLimits(0.20);
    try std.testing.expectEqual(
        @as(f64, 0.20),
        limits.soil_fraction_per_iteration,
    );
    try std.testing.expectEqual(
        @as(f64, 0.20),
        limits.surface_litter_fraction_per_iteration,
    );
}

test "SOLUTE substrate-product alias rejects invalid runtime fractions" {
    try std.testing.expectError(
        error.InvalidSoluteSubstrateProductLimit,
        initializeSubstrateProductLimits(-0.01),
    );
    try std.testing.expectError(
        error.InvalidSoluteSubstrateProductLimit,
        initializeSubstrateProductLimits(1.01),
    );
    try std.testing.expectError(
        error.InvalidSoluteSubstrateProductLimit,
        initializeSubstrateProductLimits(std.math.nan(f64)),
    );
}
