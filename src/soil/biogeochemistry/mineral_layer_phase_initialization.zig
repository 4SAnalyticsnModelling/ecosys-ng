const std = @import("std");

pub const InitialContent = union(enum) {
    saturation,
    field_capacity,
    wilting_point,
    minimum,
    fraction_m3_per_m3: f64,
};

pub const PhaseInitialization = enum {
    initialize,
    preserve_existing,
};

pub const MaterialInputs = struct {
    bulk_density_megagrams_per_m3: f64,
    matrix_bulk_volume_m3: f64,
    initial_total_volume_m3: f64,
    matrix_fraction: f64,
    macropore_fraction: f64,
    rock_volume_fraction: f64,
    sand_mass_fraction: f64,
    silt_mass_fraction: f64,
    clay_mass_fraction: f64,
    humus_carbon_concentration_g_c_per_megagram: []const f64,
};

pub const PhaseInputs = struct {
    mode: PhaseInitialization,
    initial_liquid: InitialContent,
    initial_ice: InitialContent,
    field_capacity_fraction_m3_per_m3: f64,
    wilting_point_fraction_m3_per_m3: f64,
};

pub const Parameters = struct {
    saturated_matric_potential_megapascal: f64,
    calculation_floor: f64,
    humus_carbon_fraction_g_c_per_g_organic_matter: f64,
    organic_particle_density_megagrams_per_m3: f64,
    mineral_particle_density_megagrams_per_m3: f64,
    organic_heat_capacity_megajoules_per_m3_k: f64,
    nonsand_mineral_heat_capacity_megajoules_per_m3_k: f64,
    sand_and_rock_heat_capacity_megajoules_per_m3_k: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    ice_heat_capacity_megajoules_per_m3_k: f64,
    minimum_liquid_fraction_m3_per_m3: f64,
};

pub const DiagnosticState = struct {
    saturated_matric_potential_megapascal: *f64,
    oxygen_flux_rate: *f64,
    carbon_dioxide_flux_rate: *f64,
    oxygen_litter_rate: *f64,
    methane_flux_rate: *f64,
    methane_litter_rate: *f64,
    organic_carbon_transfer_rate: *f64,
};

pub const MaterialState = struct {
    particle_density_megagrams_per_m3: *f64,
    porosity_fraction: *f64,
    initial_matrix_porosity_fraction: *f64,
    matrix_pore_volume_m3: *f64,
    initial_matrix_pore_volume_m3: *f64,
    macropore_volume_m3: *f64,
    sand_mass_megagrams: *f64,
    silt_mass_megagrams: *f64,
    clay_mass_megagrams: *f64,
    dry_solid_heat_capacity_megajoules_per_k: *f64,
};

pub const PhaseState = struct {
    liquid_fraction_m3_per_m3: *f64,
    ice_fraction_m3_per_m3: *f64,
    matrix_liquid_water_m3: *f64,
    initial_matrix_liquid_water_m3: *f64,
    macropore_liquid_water_m3: *f64,
    matrix_ice_m3: *f64,
    macropore_ice_m3: *f64,
    total_air_volume_m3: *f64,
    wet_heat_capacity_megajoules_per_k: *f64,
    previous_liquid_fraction_m3_per_m3: *f64,
    previous_ice_fraction_m3_per_m3: *f64,
};

const MaterialCandidate = struct {
    particle_density_megagrams_per_m3: f64,
    porosity_fraction: f64,
    initial_matrix_porosity_fraction: f64,
    matrix_pore_volume_m3: f64,
    macropore_volume_m3: f64,
    sand_mass_megagrams: f64,
    silt_mass_megagrams: f64,
    clay_mass_megagrams: f64,
    dry_solid_heat_capacity_megajoules_per_k: f64,
};

const PhaseCandidate = struct {
    liquid_fraction_m3_per_m3: f64,
    ice_fraction_m3_per_m3: f64,
    matrix_liquid_water_m3: f64,
    macropore_liquid_water_m3: f64,
    matrix_ice_m3: f64,
    macropore_ice_m3: f64,
    total_air_volume_m3: f64,
    wet_heat_capacity_megajoules_per_k: f64,
};

fn validateScalarInputs(
    material: MaterialInputs,
    phase: PhaseInputs,
    parameters: Parameters,
) !void {
    inline for (@typeInfo(MaterialInputs).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(material, field.name)))
            return error.NonFiniteMineralLayerInitialInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(parameters, field.name)))
            return error.NonFiniteMineralLayerInitialInput;
    }
    inline for (.{
        phase.field_capacity_fraction_m3_per_m3,
        phase.wilting_point_fraction_m3_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteMineralLayerInitialInput;
    }
    inline for (.{ phase.initial_liquid, phase.initial_ice }) |selection| {
        switch (selection) {
            .fraction_m3_per_m3 => |value| {
                if (!std.math.isFinite(value))
                    return error.NonFiniteMineralLayerInitialInput;
                if (value < 0.0 or value > 1.0)
                    return error.InvalidMineralLayerInitialInput;
            },
            else => {},
        }
    }
}

fn validateDomains(
    material: MaterialInputs,
    phase: PhaseInputs,
    parameters: Parameters,
) !void {
    if (material.bulk_density_megagrams_per_m3 < 0.0 or
        material.matrix_bulk_volume_m3 < 0.0 or
        material.initial_total_volume_m3 < 0.0 or
        material.matrix_fraction < 0.0 or material.matrix_fraction > 1.0 or
        material.macropore_fraction < 0.0 or
        material.macropore_fraction > 1.0 or
        material.rock_volume_fraction < 0.0 or
        material.rock_volume_fraction > 1.0 or
        material.sand_mass_fraction < 0.0 or
        material.silt_mass_fraction < 0.0 or
        material.clay_mass_fraction < 0.0 or
        phase.field_capacity_fraction_m3_per_m3 < 0.0 or
        phase.field_capacity_fraction_m3_per_m3 > 1.0 or
        phase.wilting_point_fraction_m3_per_m3 < 0.0 or
        phase.wilting_point_fraction_m3_per_m3 > 1.0 or
        parameters.calculation_floor < 0.0 or
        parameters.humus_carbon_fraction_g_c_per_g_organic_matter <= 0.0 or
        parameters.organic_particle_density_megagrams_per_m3 <= 0.0 or
        parameters.mineral_particle_density_megagrams_per_m3 <= 0.0 or
        parameters.organic_heat_capacity_megajoules_per_m3_k < 0.0 or
        parameters.nonsand_mineral_heat_capacity_megajoules_per_m3_k < 0.0 or
        parameters.sand_and_rock_heat_capacity_megajoules_per_m3_k < 0.0 or
        parameters.liquid_water_heat_capacity_megajoules_per_m3_k < 0.0 or
        parameters.ice_heat_capacity_megajoules_per_m3_k < 0.0 or
        parameters.minimum_liquid_fraction_m3_per_m3 < 0.0)
        return error.InvalidMineralLayerInitialInput;
    for (material.humus_carbon_concentration_g_c_per_megagram) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteMineralLayerInitialInput;
        if (value < 0.0) return error.InvalidMineralLayerInitialInput;
    }
}

fn calculateMaterial(
    material: MaterialInputs,
    parameters: Parameters,
) MaterialCandidate {
    var humus_carbon_g_c_per_megagram: f64 = 0.0;
    for (material.humus_carbon_concentration_g_c_per_megagram) |value|
        humus_carbon_g_c_per_megagram += value;
    const organic_matter_g_per_megagram = @min(
        1.0e6,
        humus_carbon_g_c_per_megagram /
            parameters.humus_carbon_fraction_g_c_per_g_organic_matter,
    );

    if (material.bulk_density_megagrams_per_m3 > parameters.calculation_floor) {
        const particle_density_megagrams_per_m3 = 1.0e-6 *
            (parameters.organic_particle_density_megagrams_per_m3 *
                organic_matter_g_per_megagram +
                parameters.mineral_particle_density_megagrams_per_m3 *
                    (1.0e6 - organic_matter_g_per_megagram));
        const porosity_fraction =
            1.0 - material.bulk_density_megagrams_per_m3 /
                particle_density_megagrams_per_m3;
        const organic_volume_fraction = organic_matter_g_per_megagram * 1.0e-6 *
            material.bulk_density_megagrams_per_m3 / particle_density_megagrams_per_m3;
        const nonsand_volume_fraction =
            (material.silt_mass_fraction + material.clay_mass_fraction) *
            material.bulk_density_megagrams_per_m3 / particle_density_megagrams_per_m3;
        const sand_volume_fraction = material.sand_mass_fraction *
            material.bulk_density_megagrams_per_m3 / particle_density_megagrams_per_m3;
        return .{
            .particle_density_megagrams_per_m3 = particle_density_megagrams_per_m3,
            .porosity_fraction = porosity_fraction,
            .initial_matrix_porosity_fraction = porosity_fraction * material.matrix_fraction,
            .matrix_pore_volume_m3 = porosity_fraction * material.matrix_bulk_volume_m3,
            .macropore_volume_m3 = material.macropore_fraction * material.initial_total_volume_m3,
            .sand_mass_megagrams = material.sand_mass_fraction * material.matrix_bulk_volume_m3,
            .silt_mass_megagrams = material.silt_mass_fraction * material.matrix_bulk_volume_m3,
            .clay_mass_megagrams = material.clay_mass_fraction * material.matrix_bulk_volume_m3,
            .dry_solid_heat_capacity_megajoules_per_k = ((parameters.organic_heat_capacity_megajoules_per_m3_k *
                organic_volume_fraction +
                parameters.nonsand_mineral_heat_capacity_megajoules_per_m3_k *
                    nonsand_volume_fraction +
                parameters.sand_and_rock_heat_capacity_megajoules_per_m3_k *
                    sand_volume_fraction) *
                material.matrix_fraction +
                parameters.sand_and_rock_heat_capacity_megajoules_per_m3_k *
                    material.rock_volume_fraction) *
                material.initial_total_volume_m3,
        };
    }
    return .{
        .particle_density_megagrams_per_m3 = 0.0,
        .porosity_fraction = 1.0,
        .initial_matrix_porosity_fraction = material.matrix_fraction,
        .matrix_pore_volume_m3 = material.matrix_bulk_volume_m3,
        .macropore_volume_m3 = material.macropore_fraction * material.initial_total_volume_m3,
        .sand_mass_megagrams = material.sand_mass_fraction * material.matrix_bulk_volume_m3,
        .silt_mass_megagrams = material.silt_mass_fraction * material.matrix_bulk_volume_m3,
        .clay_mass_megagrams = material.clay_mass_fraction * material.matrix_bulk_volume_m3,
        .dry_solid_heat_capacity_megajoules_per_k = 0.0,
    };
}

fn selectedLiquidFraction(
    selection: InitialContent,
    material: MaterialCandidate,
    phase: PhaseInputs,
    parameters: Parameters,
) f64 {
    return switch (selection) {
        .saturation => material.porosity_fraction,
        .field_capacity => phase.field_capacity_fraction_m3_per_m3,
        .wilting_point => phase.wilting_point_fraction_m3_per_m3,
        .minimum => parameters.minimum_liquid_fraction_m3_per_m3,
        .fraction_m3_per_m3 => |value| value,
    };
}

fn selectedIceFraction(
    selection: InitialContent,
    liquid_fraction: f64,
    material: MaterialCandidate,
    phase: PhaseInputs,
) f64 {
    const available = material.porosity_fraction - liquid_fraction;
    return switch (selection) {
        .saturation => @max(0.0, @min(
            material.porosity_fraction,
            available,
        )),
        .field_capacity => @max(0.0, @min(
            phase.field_capacity_fraction_m3_per_m3,
            available,
        )),
        .wilting_point => @max(0.0, @min(
            phase.wilting_point_fraction_m3_per_m3,
            available,
        )),
        .minimum => 0.0,
        .fraction_m3_per_m3 => |value| value,
    };
}

fn calculatePhase(
    material_inputs: MaterialInputs,
    phase_inputs: PhaseInputs,
    parameters: Parameters,
    material: MaterialCandidate,
) PhaseCandidate {
    const liquid_fraction = selectedLiquidFraction(
        phase_inputs.initial_liquid,
        material,
        phase_inputs,
        parameters,
    );
    const ice_fraction = selectedIceFraction(
        phase_inputs.initial_ice,
        liquid_fraction,
        material,
        phase_inputs,
    );
    const matrix_liquid_m3 =
        liquid_fraction * material_inputs.matrix_bulk_volume_m3;
    const macropore_liquid_m3 =
        liquid_fraction * material.macropore_volume_m3;
    const matrix_ice_m3 =
        ice_fraction * material_inputs.matrix_bulk_volume_m3;
    const macropore_ice_m3 = ice_fraction * material.macropore_volume_m3;
    const total_air_m3 =
        @max(
            0.0,
            material.matrix_pore_volume_m3 -
                matrix_liquid_m3 -
                matrix_ice_m3,
        ) +
        @max(
            0.0,
            material.macropore_volume_m3 -
                macropore_liquid_m3 -
                macropore_ice_m3,
        );
    return .{
        .liquid_fraction_m3_per_m3 = liquid_fraction,
        .ice_fraction_m3_per_m3 = ice_fraction,
        .matrix_liquid_water_m3 = matrix_liquid_m3,
        .macropore_liquid_water_m3 = macropore_liquid_m3,
        .matrix_ice_m3 = matrix_ice_m3,
        .macropore_ice_m3 = macropore_ice_m3,
        .total_air_volume_m3 = total_air_m3,
        .wet_heat_capacity_megajoules_per_k = material.dry_solid_heat_capacity_megajoules_per_k +
            parameters.liquid_water_heat_capacity_megajoules_per_m3_k *
                (matrix_liquid_m3 + macropore_liquid_m3) +
            parameters.ice_heat_capacity_megajoules_per_m3_k *
                (matrix_ice_m3 + macropore_ice_m3),
    };
}

fn validateCandidate(candidate: anytype) !void {
    inline for (@typeInfo(@TypeOf(candidate)).@"struct".fields) |field| {
        const value = @field(candidate, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteMineralLayerInitialResult;
        if (value < 0.0) return error.InvalidMineralLayerInitialResult;
    }
}

fn validatePhysicalState(
    material: MaterialCandidate,
    phase: ?PhaseCandidate,
) !void {
    if (material.porosity_fraction > 1.0 or
        material.initial_matrix_porosity_fraction > 1.0)
        return error.InvalidMineralLayerInitialResult;
    if (phase) |value| {
        if (value.liquid_fraction_m3_per_m3 > material.porosity_fraction or
            value.ice_fraction_m3_per_m3 > material.porosity_fraction or
            value.liquid_fraction_m3_per_m3 +
                value.ice_fraction_m3_per_m3 >
                material.porosity_fraction)
            return error.InvalidMineralLayerInitialResult;
    }
}

/// Exact mineral-layer translation of legacy `STARTS` lines 1098--1190.
pub fn initialize(
    diagnostics: DiagnosticState,
    material_state: MaterialState,
    phase_state: PhaseState,
    material_inputs: MaterialInputs,
    phase_inputs: PhaseInputs,
    parameters: Parameters,
) !void {
    try validateScalarInputs(material_inputs, phase_inputs, parameters);
    try validateDomains(material_inputs, phase_inputs, parameters);
    const material = calculateMaterial(material_inputs, parameters);
    try validateCandidate(material);
    const phase =
        if (phase_inputs.mode == .initialize)
            calculatePhase(material_inputs, phase_inputs, parameters, material)
        else
            null;
    if (phase) |candidate| try validateCandidate(candidate);
    try validatePhysicalState(material, phase);

    diagnostics.saturated_matric_potential_megapascal.* =
        parameters.saturated_matric_potential_megapascal;
    diagnostics.oxygen_flux_rate.* = 0.0;
    diagnostics.carbon_dioxide_flux_rate.* = 0.0;
    diagnostics.oxygen_litter_rate.* = 0.0;
    diagnostics.methane_flux_rate.* = 0.0;
    diagnostics.methane_litter_rate.* = 0.0;
    diagnostics.organic_carbon_transfer_rate.* = 0.0;

    material_state.particle_density_megagrams_per_m3.* =
        material.particle_density_megagrams_per_m3;
    material_state.porosity_fraction.* = material.porosity_fraction;
    material_state.initial_matrix_porosity_fraction.* =
        material.initial_matrix_porosity_fraction;
    material_state.matrix_pore_volume_m3.* = material.matrix_pore_volume_m3;
    material_state.initial_matrix_pore_volume_m3.* =
        material.matrix_pore_volume_m3;
    material_state.macropore_volume_m3.* = material.macropore_volume_m3;
    material_state.sand_mass_megagrams.* = material.sand_mass_megagrams;
    material_state.silt_mass_megagrams.* = material.silt_mass_megagrams;
    material_state.clay_mass_megagrams.* = material.clay_mass_megagrams;
    material_state.dry_solid_heat_capacity_megajoules_per_k.* =
        material.dry_solid_heat_capacity_megajoules_per_k;

    if (phase) |candidate| {
        phase_state.liquid_fraction_m3_per_m3.* =
            candidate.liquid_fraction_m3_per_m3;
        phase_state.ice_fraction_m3_per_m3.* =
            candidate.ice_fraction_m3_per_m3;
        phase_state.matrix_liquid_water_m3.* =
            candidate.matrix_liquid_water_m3;
        phase_state.initial_matrix_liquid_water_m3.* =
            candidate.matrix_liquid_water_m3;
        phase_state.macropore_liquid_water_m3.* =
            candidate.macropore_liquid_water_m3;
        phase_state.matrix_ice_m3.* = candidate.matrix_ice_m3;
        phase_state.macropore_ice_m3.* = candidate.macropore_ice_m3;
        phase_state.total_air_volume_m3.* = candidate.total_air_volume_m3;
        phase_state.wet_heat_capacity_megajoules_per_k.* =
            candidate.wet_heat_capacity_megajoules_per_k;
        phase_state.previous_liquid_fraction_m3_per_m3.* =
            candidate.liquid_fraction_m3_per_m3;
        phase_state.previous_ice_fraction_m3_per_m3.* =
            candidate.ice_fraction_m3_per_m3;
    }
}

fn testStates(values: []f64) struct {
    diagnostics: DiagnosticState,
    material: MaterialState,
    phase: PhaseState,
} {
    return .{
        .diagnostics = .{
            .saturated_matric_potential_megapascal = &values[0],
            .oxygen_flux_rate = &values[1],
            .carbon_dioxide_flux_rate = &values[2],
            .oxygen_litter_rate = &values[3],
            .methane_flux_rate = &values[4],
            .methane_litter_rate = &values[5],
            .organic_carbon_transfer_rate = &values[6],
        },
        .material = .{
            .particle_density_megagrams_per_m3 = &values[7],
            .porosity_fraction = &values[8],
            .initial_matrix_porosity_fraction = &values[9],
            .matrix_pore_volume_m3 = &values[10],
            .initial_matrix_pore_volume_m3 = &values[11],
            .macropore_volume_m3 = &values[12],
            .sand_mass_megagrams = &values[13],
            .silt_mass_megagrams = &values[14],
            .clay_mass_megagrams = &values[15],
            .dry_solid_heat_capacity_megajoules_per_k = &values[16],
        },
        .phase = .{
            .liquid_fraction_m3_per_m3 = &values[17],
            .ice_fraction_m3_per_m3 = &values[18],
            .matrix_liquid_water_m3 = &values[19],
            .initial_matrix_liquid_water_m3 = &values[20],
            .macropore_liquid_water_m3 = &values[21],
            .matrix_ice_m3 = &values[22],
            .macropore_ice_m3 = &values[23],
            .total_air_volume_m3 = &values[24],
            .wet_heat_capacity_megajoules_per_k = &values[25],
            .previous_liquid_fraction_m3_per_m3 = &values[26],
            .previous_ice_fraction_m3_per_m3 = &values[27],
        },
    };
}

const test_parameters: Parameters = .{
    .saturated_matric_potential_megapascal = -0.001,
    .calculation_floor = 1.0e-12,
    .humus_carbon_fraction_g_c_per_g_organic_matter = 0.55,
    .organic_particle_density_megagrams_per_m3 = 1.30,
    .mineral_particle_density_megagrams_per_m3 = 2.66,
    .organic_heat_capacity_megajoules_per_m3_k = 2.496,
    .nonsand_mineral_heat_capacity_megajoules_per_m3_k = 2.385,
    .sand_and_rock_heat_capacity_megajoules_per_m3_k = 2.128,
    .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
    .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
    .minimum_liquid_fraction_m3_per_m3 = 1.0e-3,
};

test "STARTS initializes mineral pore water ice air and heat state" {
    var values = [_]f64{9.0} ** 28;
    const state = testStates(&values);
    try initialize(
        state.diagnostics,
        state.material,
        state.phase,
        .{
            .bulk_density_megagrams_per_m3 = 1.2,
            .matrix_bulk_volume_m3 = 10.0,
            .initial_total_volume_m3 = 12.0,
            .matrix_fraction = 0.9,
            .macropore_fraction = 0.05,
            .rock_volume_fraction = 0.02,
            .sand_mass_fraction = 0.5,
            .silt_mass_fraction = 0.3,
            .clay_mass_fraction = 0.2,
            .humus_carbon_concentration_g_c_per_megagram = &.{
                10_000,
                20_000,
                5_000,
                15_000,
            },
        },
        .{
            .mode = .initialize,
            .initial_liquid = .field_capacity,
            .initial_ice = .wilting_point,
            .field_capacity_fraction_m3_per_m3 = 0.30,
            .wilting_point_fraction_m3_per_m3 = 0.10,
        },
        test_parameters,
    );
    try std.testing.expectEqual(@as(f64, -0.001), values[0]);
    try std.testing.expectEqual(@as(f64, 0.0), values[1]);
    try std.testing.expectEqual(@as(f64, 0.30), values[17]);
    try std.testing.expectEqual(@as(f64, 0.10), values[18]);
    try std.testing.expectEqual(@as(f64, 3.0), values[19]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.18), values[21], 1.0e-15);
    try std.testing.expectEqual(@as(f64, 1.0), values[22]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.06), values[23], 1.0e-15);
}

test "phase preservation retains existing phase values" {
    var values = [_]f64{9.0} ** 28;
    const state = testStates(&values);
    try initialize(
        state.diagnostics,
        state.material,
        state.phase,
        .{
            .bulk_density_megagrams_per_m3 = 0.0,
            .matrix_bulk_volume_m3 = 2.0,
            .initial_total_volume_m3 = 3.0,
            .matrix_fraction = 0.8,
            .macropore_fraction = 0.1,
            .rock_volume_fraction = 0.0,
            .sand_mass_fraction = 0.5,
            .silt_mass_fraction = 0.3,
            .clay_mass_fraction = 0.2,
            .humus_carbon_concentration_g_c_per_megagram = &.{0},
        },
        .{
            .mode = .preserve_existing,
            .initial_liquid = .minimum,
            .initial_ice = .minimum,
            .field_capacity_fraction_m3_per_m3 = 0.3,
            .wilting_point_fraction_m3_per_m3 = 0.1,
        },
        test_parameters,
    );
    try std.testing.expectEqual(@as(f64, 1.0), values[8]);
    try std.testing.expectEqual(@as(f64, 2.0), values[10]);
    const expected_phase = [_]f64{9.0} ** 11;
    try std.testing.expectEqualSlices(f64, &expected_phase, values[17..28]);
}

test "invalid explicit fraction fails before mutation" {
    var values = [_]f64{9.0} ** 28;
    const state = testStates(&values);
    try std.testing.expectError(
        error.InvalidMineralLayerInitialInput,
        initialize(
            state.diagnostics,
            state.material,
            state.phase,
            .{
                .bulk_density_megagrams_per_m3 = 1.2,
                .matrix_bulk_volume_m3 = 10.0,
                .initial_total_volume_m3 = 12.0,
                .matrix_fraction = 0.9,
                .macropore_fraction = 0.05,
                .rock_volume_fraction = 0.02,
                .sand_mass_fraction = 0.5,
                .silt_mass_fraction = 0.3,
                .clay_mass_fraction = 0.2,
                .humus_carbon_concentration_g_c_per_megagram = &.{100},
            },
            .{
                .mode = .initialize,
                .initial_liquid = .{ .fraction_m3_per_m3 = -0.1 },
                .initial_ice = .minimum,
                .field_capacity_fraction_m3_per_m3 = 0.3,
                .wilting_point_fraction_m3_per_m3 = 0.1,
            },
            test_parameters,
        ),
    );
    const expected = [_]f64{9.0} ** 28;
    try std.testing.expectEqualSlices(f64, &expected, &values);
}
