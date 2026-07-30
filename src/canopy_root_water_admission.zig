const std = @import("std");

pub const Inputs = struct {
    emergence_day: u32,
    wet_canopy_heat_capacity_mj_per_k: f64,
    minimum_canopy_heat_capacity_mj_per_k: f64,
    absorbed_radiation_fraction: f64,
    negligible_radiation_fraction: f64,
    primary_root_depth_m: f64,
    seeding_depth_m: f64,
    soil_surface_reference_elevation_m: f64,
    canopy_height_m: f64,
    canopy_total_water_potential_mpa: f64,
    wood_elastic_modulus_mpa: f64,
    hydraulic_canopy_height_fraction: f64,
    gravitational_water_potential_mpa_per_m: f64,
    conducting_radius_scale: f64,
    minimum_conducting_radius_factor: f64,
    conducting_radius_exponent: f64,
};

pub const Result = struct {
    total_soil_root_conductance_m_per_h_mpa: f64,
    hydraulic_canopy_height_m: f64,
    gravitational_canopy_water_potential_mpa: f64,
    stalk_to_root_conducting_radius_factor: f64,
};

/// UPTAKE.F 642--661. Returns null when any strict source admission condition
/// fails; otherwise returns the source-ordered hydraulic initialization.
pub fn calculate(inputs: Inputs) !?Result {
    try validate(inputs);
    if (inputs.emergence_day == 0 or
        inputs.wet_canopy_heat_capacity_mj_per_k <=
            inputs.minimum_canopy_heat_capacity_mj_per_k or
        inputs.absorbed_radiation_fraction <=
            inputs.negligible_radiation_fraction or
        inputs.primary_root_depth_m <=
            inputs.seeding_depth_m +
                inputs.soil_surface_reference_elevation_m)
        return null;

    const total_conductance: f64 = 0;
    const hydraulic_height =
        inputs.hydraulic_canopy_height_fraction * inputs.canopy_height_m;
    const gravitational_potential =
        -inputs.gravitational_water_potential_mpa_per_m * hydraulic_height;
    const radius_factor = inputs.conducting_radius_scale *
        std.math.pow(
            f64,
            @max(
                inputs.minimum_conducting_radius_factor,
                1 +
                    inputs.canopy_total_water_potential_mpa /
                        inputs.wood_elastic_modulus_mpa,
            ),
            inputs.conducting_radius_exponent,
        );
    const result = Result{
        .total_soil_root_conductance_m_per_h_mpa = total_conductance,
        .hydraulic_canopy_height_m = hydraulic_height,
        .gravitational_canopy_water_potential_mpa = gravitational_potential,
        .stalk_to_root_conducting_radius_factor = radius_factor,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyRootWaterAdmission;
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyRootWaterAdmissionInput;
    if (inputs.wet_canopy_heat_capacity_mj_per_k < 0 or
        inputs.minimum_canopy_heat_capacity_mj_per_k < 0 or
        inputs.absorbed_radiation_fraction < 0 or
        inputs.negligible_radiation_fraction < 0 or
        inputs.primary_root_depth_m < 0 or
        inputs.seeding_depth_m < 0 or
        inputs.canopy_height_m < 0 or
        inputs.wood_elastic_modulus_mpa <= 0 or
        inputs.hydraulic_canopy_height_fraction < 0 or
        inputs.gravitational_water_potential_mpa_per_m < 0 or
        inputs.conducting_radius_scale < 0 or
        inputs.minimum_conducting_radius_factor < 0 or
        inputs.conducting_radius_exponent <= 0)
        return error.InvalidCanopyRootWaterAdmissionInput;
}

fn admittedInputs() Inputs {
    return .{
        .emergence_day = 120,
        .wet_canopy_heat_capacity_mj_per_k = 2,
        .minimum_canopy_heat_capacity_mj_per_k = 1,
        .absorbed_radiation_fraction = 0.4,
        .negligible_radiation_fraction = 1e-12,
        .primary_root_depth_m = 0.6,
        .seeding_depth_m = 0.05,
        .soil_surface_reference_elevation_m = 0,
        .canopy_height_m = 2,
        .canopy_total_water_potential_mpa = -10,
        .wood_elastic_modulus_mpa = 50,
        .hydraulic_canopy_height_fraction = 0.8,
        .gravitational_water_potential_mpa_per_m = 0.0098,
        .conducting_radius_scale = 3.75e3,
        .minimum_conducting_radius_factor = 0.5,
        .conducting_radius_exponent = 4,
    };
}

test "UPTAKE admitted canopy initializes hydraulic preamble in source order" {
    const result = (try calculate(admittedInputs())).?;
    try std.testing.expectEqual(@as(f64, 0), result.total_soil_root_conductance_m_per_h_mpa);
    try std.testing.expectEqual(@as(f64, 1.6), result.hydraulic_canopy_height_m);
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.01568),
        result.gravitational_canopy_water_potential_mpa,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 3.75e3 * 0.8 * 0.8 * 0.8 * 0.8),
        result.stalk_to_root_conducting_radius_factor,
        1e-12,
    );
}

test "UPTAKE convergence admission preserves all four strict gates" {
    var inputs = admittedInputs();
    inputs.emergence_day = 0;
    try std.testing.expect((try calculate(inputs)) == null);
    inputs = admittedInputs();
    inputs.wet_canopy_heat_capacity_mj_per_k =
        inputs.minimum_canopy_heat_capacity_mj_per_k;
    try std.testing.expect((try calculate(inputs)) == null);
    inputs = admittedInputs();
    inputs.absorbed_radiation_fraction =
        inputs.negligible_radiation_fraction;
    try std.testing.expect((try calculate(inputs)) == null);
    inputs = admittedInputs();
    inputs.primary_root_depth_m =
        inputs.seeding_depth_m + inputs.soil_surface_reference_elevation_m;
    try std.testing.expect((try calculate(inputs)) == null);
}

test "conducting radius response retains source floor and fourth power" {
    var inputs = admittedInputs();
    inputs.canopy_total_water_potential_mpa = -100;
    const result = (try calculate(inputs)).?;
    try std.testing.expectEqual(
        @as(f64, 3.75e3 * 0.5 * 0.5 * 0.5 * 0.5),
        result.stalk_to_root_conducting_radius_factor,
    );
}
