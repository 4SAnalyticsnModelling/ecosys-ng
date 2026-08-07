const std = @import("std");

pub const Inputs = struct {
    ambient_air_temperature_k: f64,
    canopy_air_temperature_k: f64,
    canopy_surface_temperature_k: f64,
    minimum_richardson_number: f64,
    maximum_richardson_number: f64,
    bulk_richardson_coefficient_k: f64,
    minimum_aerodynamic_resistance_h_per_m: f64,
    maximum_aerodynamic_resistance_h_per_m: f64,
    biome_isothermal_boundary_resistance_h_per_m: f64,
    biome_below_canopy_aerodynamic_resistance_h_per_m: f64,
    species_below_canopy_aerodynamic_resistance_h_per_m: f64,
    biome_sensible_boundary_conductance_megajoules_per_m_k_step: f64,
    biome_latent_boundary_conductance_m2_per_step: f64,
    total_canopy_radiation_fraction: f64,
    species_sensible_boundary_conductance_megajoules_per_m_k_step: f64,
    species_latent_boundary_conductance_m2_per_step: f64,
    sensible_surface_resistance_h_per_m: f64,
    latent_surface_resistance_h_per_m: f64,
    richardson_response_coefficient: f64,
};

pub const Result = struct {
    ambient_to_canopy_temperature_difference_k: f64,
    ambient_richardson_number: f64,
    ambient_stability_factor: f64,
    canopy_boundary_layer_resistance_h_per_m: f64,
    above_species_aerodynamic_resistance_h_per_m: f64,
    total_canopy_aerodynamic_resistance_h_per_m: f64,
    ground_sensible_conductance_megajoules_per_m_k_step: f64,
    ground_latent_conductance_m2_per_step: f64,
    canopy_to_surface_temperature_difference_k: f64,
    surface_richardson_number: f64,
    surface_stability_factor: f64,
    adjusted_sensible_surface_resistance_h_per_m: f64,
    canopy_sensible_conductance_megajoules_per_k_step: f64,
    canopy_latent_conductance_m3_per_step: f64,
};

/// UPTAKE.F 926--941. One boundary-layer evaluation inside the runtime
/// MXN-compatible Newton/Picard canopy residual.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const ambient_difference =
        inputs.ambient_air_temperature_k - inputs.canopy_air_temperature_k;
    const ambient_richardson = @max(
        inputs.minimum_richardson_number,
        @min(
            inputs.maximum_richardson_number,
            inputs.bulk_richardson_coefficient_k /
                inputs.ambient_air_temperature_k *
                ambient_difference,
        ),
    );
    const ambient_stability =
        1 - inputs.richardson_response_coefficient * ambient_richardson;
    if (ambient_stability == 0) return error.SingularCanopyBoundaryLayerResistance;
    const boundary_resistance = @min(
        inputs.maximum_aerodynamic_resistance_h_per_m,
        @max(
            inputs.minimum_aerodynamic_resistance_h_per_m,
            inputs.biome_isothermal_boundary_resistance_h_per_m /
                ambient_stability,
        ),
    );
    const above_species_resistance = @max(
        0,
        inputs.biome_below_canopy_aerodynamic_resistance_h_per_m -
            inputs.species_below_canopy_aerodynamic_resistance_h_per_m,
    );
    const total_resistance =
        boundary_resistance + above_species_resistance;
    const ground_sensible =
        inputs.biome_sensible_boundary_conductance_megajoules_per_m_k_step /
        inputs.species_below_canopy_aerodynamic_resistance_h_per_m *
        inputs.total_canopy_radiation_fraction;
    const ground_latent =
        inputs.biome_latent_boundary_conductance_m2_per_step /
        inputs.species_below_canopy_aerodynamic_resistance_h_per_m *
        inputs.total_canopy_radiation_fraction;
    const surface_difference =
        inputs.canopy_air_temperature_k - inputs.canopy_surface_temperature_k;
    const surface_richardson = @max(
        inputs.minimum_richardson_number,
        @min(
            inputs.maximum_richardson_number,
            inputs.bulk_richardson_coefficient_k /
                inputs.canopy_air_temperature_k *
                surface_difference,
        ),
    );
    const surface_stability =
        1 - inputs.richardson_response_coefficient * surface_richardson;
    if (surface_stability == 0) return error.SingularCanopyBoundaryLayerResistance;
    const adjusted_surface_resistance = @min(
        inputs.maximum_aerodynamic_resistance_h_per_m,
        @max(
            inputs.minimum_aerodynamic_resistance_h_per_m,
            inputs.sensible_surface_resistance_h_per_m /
                surface_stability,
        ),
    );
    const canopy_sensible =
        inputs.species_sensible_boundary_conductance_megajoules_per_m_k_step /
        adjusted_surface_resistance;
    const canopy_latent =
        inputs.species_latent_boundary_conductance_m2_per_step /
        (adjusted_surface_resistance +
            inputs.latent_surface_resistance_h_per_m);
    const result = Result{
        .ambient_to_canopy_temperature_difference_k = ambient_difference,
        .ambient_richardson_number = ambient_richardson,
        .ambient_stability_factor = ambient_stability,
        .canopy_boundary_layer_resistance_h_per_m = boundary_resistance,
        .above_species_aerodynamic_resistance_h_per_m = above_species_resistance,
        .total_canopy_aerodynamic_resistance_h_per_m = total_resistance,
        .ground_sensible_conductance_megajoules_per_m_k_step = ground_sensible,
        .ground_latent_conductance_m2_per_step = ground_latent,
        .canopy_to_surface_temperature_difference_k = surface_difference,
        .surface_richardson_number = surface_richardson,
        .surface_stability_factor = surface_stability,
        .adjusted_sensible_surface_resistance_h_per_m = adjusted_surface_resistance,
        .canopy_sensible_conductance_megajoules_per_k_step = canopy_sensible,
        .canopy_latent_conductance_m3_per_step = canopy_latent,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyBoundaryLayerResistance;
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyBoundaryLayerResistanceInput;
    if (inputs.ambient_air_temperature_k <= 0 or
        inputs.canopy_air_temperature_k <= 0 or
        inputs.canopy_surface_temperature_k <= 0 or
        inputs.minimum_richardson_number >
            inputs.maximum_richardson_number or
        inputs.minimum_aerodynamic_resistance_h_per_m < 0 or
        inputs.maximum_aerodynamic_resistance_h_per_m <
            inputs.minimum_aerodynamic_resistance_h_per_m or
        inputs.species_below_canopy_aerodynamic_resistance_h_per_m <= 0 or
        inputs.latent_surface_resistance_h_per_m < 0)
        return error.InvalidCanopyBoundaryLayerResistanceInput;
}

fn sourceInputs() Inputs {
    return .{
        .ambient_air_temperature_k = 300,
        .canopy_air_temperature_k = 295,
        .canopy_surface_temperature_k = 290,
        .minimum_richardson_number = -0.1,
        .maximum_richardson_number = 0.05,
        .bulk_richardson_coefficient_k = 2,
        .minimum_aerodynamic_resistance_h_per_m = 0.01,
        .maximum_aerodynamic_resistance_h_per_m = 10,
        .biome_isothermal_boundary_resistance_h_per_m = 1,
        .biome_below_canopy_aerodynamic_resistance_h_per_m = 3,
        .species_below_canopy_aerodynamic_resistance_h_per_m = 2,
        .biome_sensible_boundary_conductance_megajoules_per_m_k_step = 8,
        .biome_latent_boundary_conductance_m2_per_step = 6,
        .total_canopy_radiation_fraction = 0.5,
        .species_sensible_boundary_conductance_megajoules_per_m_k_step = 4,
        .species_latent_boundary_conductance_m2_per_step = 5,
        .sensible_surface_resistance_h_per_m = 0.2,
        .latent_surface_resistance_h_per_m = 0.3,
        .richardson_response_coefficient = 10,
    };
}

test "UPTAKE boundary layer resistance preserves source operation order" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const ambient_ri = @max(-0.1, @min(0.05, 2.0 / 300.0 * 5.0));
    const surface_ri = @max(-0.1, @min(0.05, 2.0 / 295.0 * 5.0));
    const expected_boundary = @min(10.0, @max(0.01, 1.0 / (1 - 10 * ambient_ri)));
    const expected_surface = @min(10.0, @max(0.01, 0.2 / (1 - 10 * surface_ri)));
    try std.testing.expectEqual(ambient_ri, result.ambient_richardson_number);
    try std.testing.expectApproxEqAbs(expected_boundary, result.canopy_boundary_layer_resistance_h_per_m, 5e-16);
    try std.testing.expectEqual(@as(f64, 1), result.above_species_aerodynamic_resistance_h_per_m);
    try std.testing.expectEqual(@as(f64, 2), result.ground_sensible_conductance_megajoules_per_m_k_step);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), result.ground_latent_conductance_m2_per_step, 5e-16);
    try std.testing.expectApproxEqAbs(expected_surface, result.adjusted_sensible_surface_resistance_h_per_m, 5e-16);
    try std.testing.expectApproxEqAbs(4.0 / expected_surface, result.canopy_sensible_conductance_megajoules_per_k_step, 1e-14);
    try std.testing.expectApproxEqAbs(5.0 / (expected_surface + 0.3), result.canopy_latent_conductance_m3_per_step, 1e-14);
}

test "resistance bounds and nonnegative above-species branch are retained" {
    var inputs = sourceInputs();
    inputs.biome_below_canopy_aerodynamic_resistance_h_per_m = 1;
    inputs.species_below_canopy_aerodynamic_resistance_h_per_m = 2;
    inputs.biome_isothermal_boundary_resistance_h_per_m = 100;
    const result = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.above_species_aerodynamic_resistance_h_per_m);
    try std.testing.expectEqual(@as(f64, 10), result.canopy_boundary_layer_resistance_h_per_m);
}

test "singular Richardson stability fails explicitly" {
    var inputs = sourceInputs();
    inputs.maximum_richardson_number = 0.1;
    inputs.bulk_richardson_coefficient_k = 6;
    try std.testing.expectError(
        error.SingularCanopyBoundaryLayerResistance,
        calculate(inputs),
    );
}
