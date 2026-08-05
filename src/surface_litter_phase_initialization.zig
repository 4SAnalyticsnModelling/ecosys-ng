const std = @import("std");

pub const Inputs = struct {
    organic_carbon_g_c: f64,
    total_litter_volume_m3: f64,
    negligible_litter_volume_m3: f64,
};

pub const Parameters = struct {
    dry_mass_megagrams_per_g_c: f64,
    dry_solid_density_megagrams_per_m3: f64,
    retained_liquid_water_m3_per_g_c: f64,
    organic_carbon_heat_capacity_megajoules_per_g_c_k: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    ice_heat_capacity_megajoules_per_m3_k: f64,
};

pub const State = struct {
    dry_mass_megagrams: *f64,
    pore_volume_m3: *f64,
    liquid_water_m3: *f64,
    ice_m3: *f64,
    air_volume_m3: *f64,
    porosity_fraction: *f64,
    liquid_fraction_m3_per_m3: *f64,
    ice_fraction_m3_per_m3: *f64,
    wet_heat_capacity_megajoules_per_k: *f64,
    dry_solid_heat_capacity_megajoules_per_k: *f64,
    initial_matrix_pore_volume_m3: *f64,
};

const Candidate = struct {
    dry_mass_megagrams: f64,
    pore_volume_m3: f64,
    liquid_water_m3: f64,
    ice_m3: f64,
    air_volume_m3: f64,
    porosity_fraction: f64,
    liquid_fraction_m3_per_m3: f64,
    ice_fraction_m3_per_m3: f64,
    wet_heat_capacity_megajoules_per_k: f64,
    dry_solid_heat_capacity_megajoules_per_k: f64,
    initial_matrix_pore_volume_m3: f64,
};

fn validateInputs(inputs: Inputs, parameters: Parameters) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSurfaceLitterPhaseInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(parameters, field.name)))
            return error.NonFiniteSurfaceLitterPhaseInput;
    }
    if (inputs.organic_carbon_g_c < 0.0 or
        inputs.total_litter_volume_m3 < 0.0 or
        inputs.negligible_litter_volume_m3 < 0.0 or
        parameters.dry_mass_megagrams_per_g_c < 0.0 or
        parameters.dry_solid_density_megagrams_per_m3 <= 0.0 or
        parameters.retained_liquid_water_m3_per_g_c < 0.0 or
        parameters.organic_carbon_heat_capacity_megajoules_per_g_c_k < 0.0 or
        parameters.liquid_water_heat_capacity_megajoules_per_m3_k < 0.0 or
        parameters.ice_heat_capacity_megajoules_per_m3_k < 0.0)
        return error.InvalidSurfaceLitterPhaseInput;
}

fn calculate(inputs: Inputs, parameters: Parameters) Candidate {
    const dry_mass_megagrams =
        parameters.dry_mass_megagrams_per_g_c * inputs.organic_carbon_g_c;
    const pore_volume_m3 = @max(
        0.0,
        inputs.total_litter_volume_m3 -
            dry_mass_megagrams / parameters.dry_solid_density_megagrams_per_m3,
    );
    const liquid_water_m3 =
        parameters.retained_liquid_water_m3_per_g_c *
        inputs.organic_carbon_g_c;
    const ice_m3 = 0.0;
    const air_volume_m3 = @max(
        0.0,
        pore_volume_m3 - liquid_water_m3 - ice_m3,
    );
    const has_litter_volume =
        inputs.total_litter_volume_m3 > inputs.negligible_litter_volume_m3;
    const porosity_fraction =
        if (has_litter_volume)
            pore_volume_m3 / inputs.total_litter_volume_m3
        else
            1.0;
    const liquid_fraction =
        if (has_litter_volume)
            liquid_water_m3 / inputs.total_litter_volume_m3
        else
            0.0;
    return .{
        .dry_mass_megagrams = dry_mass_megagrams,
        .pore_volume_m3 = pore_volume_m3,
        .liquid_water_m3 = liquid_water_m3,
        .ice_m3 = ice_m3,
        .air_volume_m3 = air_volume_m3,
        .porosity_fraction = porosity_fraction,
        .liquid_fraction_m3_per_m3 = liquid_fraction,
        .ice_fraction_m3_per_m3 = 0.0,
        .wet_heat_capacity_megajoules_per_k = parameters.organic_carbon_heat_capacity_megajoules_per_g_c_k *
            inputs.organic_carbon_g_c +
            parameters.liquid_water_heat_capacity_megajoules_per_m3_k *
                liquid_water_m3 +
            parameters.ice_heat_capacity_megajoules_per_m3_k * ice_m3,
        .dry_solid_heat_capacity_megajoules_per_k = 0.0,
        .initial_matrix_pore_volume_m3 = 0.0,
    };
}

fn validateCandidate(candidate: Candidate) !void {
    inline for (@typeInfo(Candidate).@"struct".fields) |field| {
        const value = @field(candidate, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterPhaseResult;
        if (value < 0.0) return error.InvalidSurfaceLitterPhaseResult;
    }
    if (candidate.porosity_fraction > 1.0 or
        candidate.liquid_fraction_m3_per_m3 > candidate.porosity_fraction or
        candidate.ice_fraction_m3_per_m3 > candidate.porosity_fraction or
        candidate.liquid_water_m3 + candidate.ice_m3 >
            candidate.pore_volume_m3)
        return error.InvalidSurfaceLitterPhaseResult;
}

/// Exact surface-litter branch of legacy `STARTS` lines 1195--1213.
///
/// Mass is Mg, volumes are m3 per cell, volumetric fractions are m3/m3,
/// and heat capacities are MJ/K per cell.
pub fn initialize(
    state: State,
    inputs: Inputs,
    parameters: Parameters,
) !void {
    try validateInputs(inputs, parameters);
    const candidate = calculate(inputs, parameters);
    try validateCandidate(candidate);

    inline for (@typeInfo(Candidate).@"struct".fields) |field|
        @field(state, field.name).* = @field(candidate, field.name);
}

const legacy_parameters: Parameters = .{
    .dry_mass_megagrams_per_g_c = 1.82e-6,
    .dry_solid_density_megagrams_per_m3 = 1.30,
    .retained_liquid_water_m3_per_g_c = 8.0e-6,
    .organic_carbon_heat_capacity_megajoules_per_g_c_k = 2.496e-6,
    .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
    .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
};

fn stateFrom(values: []f64) State {
    return .{
        .dry_mass_megagrams = &values[0],
        .pore_volume_m3 = &values[1],
        .liquid_water_m3 = &values[2],
        .ice_m3 = &values[3],
        .air_volume_m3 = &values[4],
        .porosity_fraction = &values[5],
        .liquid_fraction_m3_per_m3 = &values[6],
        .ice_fraction_m3_per_m3 = &values[7],
        .wet_heat_capacity_megajoules_per_k = &values[8],
        .dry_solid_heat_capacity_megajoules_per_k = &values[9],
        .initial_matrix_pore_volume_m3 = &values[10],
    };
}

test "STARTS initializes surface litter mass water air and heat" {
    var values = [_]f64{9.0} ** 11;
    try initialize(stateFrom(&values), .{
        .organic_carbon_g_c = 100.0,
        .total_litter_volume_m3 = 0.002,
        .negligible_litter_volume_m3 = 1.0e-12,
    }, legacy_parameters);

    const expected_dry_mass_megagrams = 1.82e-4;
    const expected_pore_volume_m3 = 0.002 - expected_dry_mass_megagrams / 1.30;
    const expected_water_m3 = 8.0e-4;
    try std.testing.expectApproxEqAbs(
        expected_dry_mass_megagrams,
        values[0],
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(
        expected_pore_volume_m3,
        values[1],
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(expected_water_m3, values[2], 1.0e-18);
    try std.testing.expectApproxEqAbs(
        expected_pore_volume_m3 - expected_water_m3,
        values[4],
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(
        expected_pore_volume_m3 / 0.002,
        values[5],
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), values[6], 1.0e-15);
    try std.testing.expectApproxEqAbs(
        2.496e-4 + 4.19 * expected_water_m3,
        values[8],
        1.0e-15,
    );
}

test "STARTS empty surface litter uses unit porosity and zero water" {
    var values = [_]f64{9.0} ** 11;
    try initialize(stateFrom(&values), .{
        .organic_carbon_g_c = 0.0,
        .total_litter_volume_m3 = 0.0,
        .negligible_litter_volume_m3 = 1.0e-12,
    }, legacy_parameters);
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
        &values,
    );
}

test "physically overfilled litter fails atomically" {
    var values = [_]f64{9.0} ** 11;
    try std.testing.expectError(
        error.InvalidSurfaceLitterPhaseResult,
        initialize(stateFrom(&values), .{
            .organic_carbon_g_c = 100.0,
            .total_litter_volume_m3 = 0.0005,
            .negligible_litter_volume_m3 = 1.0e-12,
        }, legacy_parameters),
    );
    const expected = [_]f64{9.0} ** 11;
    try std.testing.expectEqualSlices(f64, &expected, &values);
}
