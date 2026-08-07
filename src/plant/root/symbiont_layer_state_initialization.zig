const std = @import("std");

pub const PopulationGeometry = struct {
    maximum_primary_radius_m: f64,
    maximum_secondary_radius_m: f64,
};

pub const WaterState = struct {
    uptake_m3_per_timestep: f64,
    total_potential_megapascal: f64,
    osmotic_potential_megapascal: f64,
    turgor_potential_megapascal: f64,
};

pub const NonstructuralState = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
    carbon_concentration_g_c_g_c: f64,
    nitrogen_concentration_g_n_g_c: f64,
    phosphorus_concentration_g_p_g_c: f64,
};

pub const MorphologyState = struct {
    maximum_protein_concentration_g_g_c: f64,
    total_carbon_g_c: f64,
    structural_carbon_g_c: f64,
    protein_g: f64,
    primary_root_number: f64,
    total_root_length_m: f64,
    primary_root_length_m: f64,
    root_density: f64,
    protein_volume_m3: f64,
    water_volume_m3: f64,
    primary_radius_m: f64,
    secondary_radius_m: f64,
    root_surface_area_m2: f64,
    axial_length_m: f64,
};

pub const NutrientFluxState = struct {
    ammonium_uptake_g_n_per_timestep: f64,
    nitrate_uptake_g_n_per_timestep: f64,
    dihydrogen_phosphate_uptake_g_p_per_timestep: f64,
    hydrogen_phosphate_uptake_g_p_per_timestep: f64,
    band_ammonium_uptake_g_n_per_timestep: f64,
    band_nitrate_uptake_g_n_per_timestep: f64,
    band_dihydrogen_phosphate_uptake_g_p_per_timestep: f64,
    band_hydrogen_phosphate_uptake_g_p_per_timestep: f64,
    oxygen_demand_g_o_per_timestep: f64,
    ammonium_demand_g_n_per_timestep: f64,
    band_ammonium_demand_g_n_per_timestep: f64,
    nitrate_demand_g_n_per_timestep: f64,
    band_nitrate_demand_g_n_per_timestep: f64,
    dihydrogen_phosphate_demand_g_p_per_timestep: f64,
    hydrogen_phosphate_demand_g_p_per_timestep: f64,
    band_dihydrogen_phosphate_demand_g_p_per_timestep: f64,
    band_hydrogen_phosphate_demand_g_p_per_timestep: f64,
};

pub const PopulationLayerState = struct {
    water: WaterState,
    nonstructural: NonstructuralState,
    morphology: MorphologyState,
    nutrient_flux: NutrientFluxState,
};

pub const SpeciesTotals = struct {
    active_root_count: usize,
    ammonium_uptake_g_n_per_timestep: f64,
    nitrate_uptake_g_n_per_timestep: f64,
    dihydrogen_phosphate_uptake_g_p_per_timestep: f64,
    hydrogen_phosphate_uptake_g_p_per_timestep: f64,
    nitrogen_fixation_g_n_per_timestep: f64,
};

pub const Parameters = struct {
    initial_total_water_potential_megapascal: f64,
    initial_axial_length_m: f64,
};

pub const InitializationError = error{
    StateExtentMismatch,
    NonFiniteInput,
    InvalidRadius,
    InvalidParameter,
    NonFiniteResult,
};

/// Translates `startq.f` lines 741--790 with runtime population and soil-layer
/// extents. Flattened state order is population-major, then layer.
pub fn initialize(
    populations: []const PopulationGeometry,
    layer_count: usize,
    leaf_osmotic_potential_at_zero_total_megapascal: f64,
    maximum_root_protein_concentration_g_g_c: f64,
    parameters: Parameters,
    states: []PopulationLayerState,
) InitializationError!SpeciesTotals {
    const expected_state_count = std.math.mul(usize, populations.len, layer_count) catch
        return error.StateExtentMismatch;
    if (states.len != expected_state_count) return error.StateExtentMismatch;
    try validate(
        populations,
        leaf_osmotic_potential_at_zero_total_megapascal,
        maximum_root_protein_concentration_g_g_c,
        parameters,
    );

    for (populations, 0..) |population, population_index| {
        for (0..layer_count) |layer_index| {
            const total_potential_megapascal = parameters.initial_total_water_potential_megapascal;
            const osmotic_potential_megapascal =
                leaf_osmotic_potential_at_zero_total_megapascal + total_potential_megapascal;
            const turgor_potential_megapascal =
                @max(0.0, total_potential_megapascal - osmotic_potential_megapascal);
            const state = PopulationLayerState{
                .water = .{
                    .uptake_m3_per_timestep = 0.0,
                    .total_potential_megapascal = total_potential_megapascal,
                    .osmotic_potential_megapascal = osmotic_potential_megapascal,
                    .turgor_potential_megapascal = turgor_potential_megapascal,
                },
                .nonstructural = std.mem.zeroes(NonstructuralState),
                .morphology = .{
                    .maximum_protein_concentration_g_g_c = maximum_root_protein_concentration_g_g_c,
                    .total_carbon_g_c = 0.0,
                    .structural_carbon_g_c = 0.0,
                    .protein_g = 0.0,
                    .primary_root_number = 0.0,
                    .total_root_length_m = 0.0,
                    .primary_root_length_m = 0.0,
                    .root_density = 0.0,
                    .protein_volume_m3 = 0.0,
                    .water_volume_m3 = 0.0,
                    .primary_radius_m = population.maximum_primary_radius_m,
                    .secondary_radius_m = population.maximum_secondary_radius_m,
                    .root_surface_area_m2 = 0.0,
                    .axial_length_m = parameters.initial_axial_length_m,
                },
                .nutrient_flux = std.mem.zeroes(NutrientFluxState),
            };
            inline for (std.meta.fields(WaterState)) |field| {
                if (!std.math.isFinite(@field(state.water, field.name))) {
                    return error.NonFiniteResult;
                }
            }
            states[population_index * layer_count + layer_index] = state;
        }
    }
    return std.mem.zeroes(SpeciesTotals);
}

fn validate(
    populations: []const PopulationGeometry,
    leaf_osmotic_potential_at_zero_total_megapascal: f64,
    maximum_root_protein_concentration_g_g_c: f64,
    parameters: Parameters,
) InitializationError!void {
    inline for (.{
        leaf_osmotic_potential_at_zero_total_megapascal,
        maximum_root_protein_concentration_g_g_c,
        parameters.initial_total_water_potential_megapascal,
        parameters.initial_axial_length_m,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
    }
    if (maximum_root_protein_concentration_g_g_c < 0.0 or
        parameters.initial_axial_length_m < 0.0 or
        parameters.initial_total_water_potential_megapascal > 0.0)
    {
        return error.InvalidParameter;
    }
    for (populations) |population| {
        inline for (std.meta.fields(PopulationGeometry)) |field| {
            const radius_m = @field(population, field.name);
            if (!std.math.isFinite(radius_m)) return error.NonFiniteInput;
            if (radius_m <= 0.0) return error.InvalidRadius;
        }
    }
}

test "runtime populations and layers preserve STARTQ initial states" {
    const populations = [_]PopulationGeometry{
        .{ .maximum_primary_radius_m = 1.0e-4, .maximum_secondary_radius_m = 5.0e-5 },
        .{ .maximum_primary_radius_m = 2.5e-6, .maximum_secondary_radius_m = 2.5e-6 },
    };
    var states: [6]PopulationLayerState = undefined;
    const totals = try initialize(
        &populations,
        3,
        -0.4,
        0.2,
        .{ .initial_total_water_potential_megapascal = -0.01, .initial_axial_length_m = 1.0e-3 },
        &states,
    );

    try std.testing.expectEqual(std.mem.zeroes(SpeciesTotals), totals);
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.41),
        states[0].water.osmotic_potential_megapascal,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.4),
        states[0].water.turgor_potential_megapascal,
        1.0e-15,
    );
    try std.testing.expectEqual(@as(f64, 1.0e-4), states[2].morphology.primary_radius_m);
    try std.testing.expectEqual(@as(f64, 2.5e-6), states[3].morphology.primary_radius_m);
    try std.testing.expectEqual(@as(f64, 1.0e-3), states[5].morphology.axial_length_m);
}

test "invalid extent fails before state mutation" {
    const populations = [_]PopulationGeometry{
        .{ .maximum_primary_radius_m = 1.0e-4, .maximum_secondary_radius_m = 5.0e-5 },
    };
    var states: [1]PopulationLayerState = undefined;
    @memset(std.mem.asBytes(&states), 0xff);
    try std.testing.expectError(error.StateExtentMismatch, initialize(
        &populations,
        2,
        -0.4,
        0.2,
        .{ .initial_total_water_potential_megapascal = -0.01, .initial_axial_length_m = 1.0e-3 },
        &states,
    ));
}
