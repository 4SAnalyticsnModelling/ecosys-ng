const std = @import("std");

pub const Components = struct {
    leaf_shortwave: []f64,
    stalk_shortwave: []f64,
    standing_dead_shortwave: []f64,
    leaf_par: []f64,
    stalk_par: []f64,
    standing_dead_par: []f64,
};

pub const SpeciesAccumulators = struct {
    absorbed: Components,
    backscattered: Components,
    forward_scattered: Components,
};

pub const AngularPar = struct {
    direct_par_umol_per_m2_s: []const f64,
    diffuse_par_umol_per_m2_s: []f64,
    total_par_umol_per_m2_s: []f64,
};

pub const State = struct {
    diffuse_transmittance_above: f64,
    current_layer_interception_fraction: f64,
};

/// HOUR1 lines 1482--1511. Angular PAR layout is
/// species -> inclination -> azimuth, matching source NZ -> N -> M traversal.
pub fn apply(
    preceding_interception_fraction: f64,
    inclination_count: usize,
    azimuth_count: usize,
    state: *State,
    accumulators: SpeciesAccumulators,
    angular_par: AngularPar,
) !?f64 {
    const species_count = try validate(
        preceding_interception_fraction,
        inclination_count,
        azimuth_count,
        state.*,
        accumulators,
        angular_par,
    );
    if (preceding_interception_fraction +
        state.current_layer_interception_fraction <= 1.0)
        return null;
    const normalization =
        (1.0 - preceding_interception_fraction) /
        ((1.0 - preceding_interception_fraction) -
            (1.0 - preceding_interception_fraction -
                state.current_layer_interception_fraction));
    if (!std.math.isFinite(normalization) or normalization < 0)
        return error.InvalidDiffuseInterceptionNormalization;

    state.diffuse_transmittance_above =
        state.diffuse_transmittance_above * normalization;
    state.current_layer_interception_fraction =
        state.current_layer_interception_fraction * normalization;
    for (0..species_count) |species| {
        scaleAt(accumulators.absorbed, species, normalization);
        scaleAt(accumulators.backscattered, species, normalization);
        scaleAt(accumulators.forward_scattered, species, normalization);
        for (0..inclination_count) |inclination| {
            for (0..azimuth_count) |azimuth| {
                const index = (species * inclination_count + inclination) *
                    azimuth_count + azimuth;
                angular_par.diffuse_par_umol_per_m2_s[index] =
                    angular_par.diffuse_par_umol_per_m2_s[index] * normalization;
                angular_par.total_par_umol_per_m2_s[index] =
                    angular_par.direct_par_umol_per_m2_s[index] +
                    angular_par.diffuse_par_umol_per_m2_s[index];
            }
        }
    }
    return normalization;
}

fn scaleAt(components: Components, species: usize, scale: f64) void {
    components.leaf_shortwave[species] =
        components.leaf_shortwave[species] * scale;
    components.stalk_shortwave[species] =
        components.stalk_shortwave[species] * scale;
    components.standing_dead_shortwave[species] =
        components.standing_dead_shortwave[species] * scale;
    components.leaf_par[species] = components.leaf_par[species] * scale;
    components.stalk_par[species] = components.stalk_par[species] * scale;
    components.standing_dead_par[species] =
        components.standing_dead_par[species] * scale;
}

fn validate(
    preceding: f64,
    inclination_count: usize,
    azimuth_count: usize,
    state: State,
    accumulators: SpeciesAccumulators,
    angular_par: AngularPar,
) !usize {
    inline for (.{
        preceding,
        state.diffuse_transmittance_above,
        state.current_layer_interception_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDiffuseInterceptionNormalizationInput;
    if (inclination_count == 0 or azimuth_count == 0)
        return error.ZeroDiffuseInterceptionAngularExtent;
    const species_count = accumulators.absorbed.leaf_shortwave.len;
    if (species_count == 0)
        return error.ZeroDiffuseInterceptionSpeciesExtent;
    inline for (@typeInfo(SpeciesAccumulators).@"struct".fields) |group_field| {
        const components = @field(accumulators, group_field.name);
        inline for (@typeInfo(Components).@"struct".fields) |field| {
            const values = @field(components, field.name);
            if (values.len != species_count)
                return error.DiffuseInterceptionNormalizationDimensionMismatch;
            for (values) |value| if (!std.math.isFinite(value))
                return error.InvalidDiffuseInterceptionNormalizationInput;
        }
    }
    const angle_count = try std.math.mul(usize, inclination_count, azimuth_count);
    const par_count = try std.math.mul(usize, species_count, angle_count);
    if (angular_par.direct_par_umol_per_m2_s.len != par_count or
        angular_par.diffuse_par_umol_per_m2_s.len != par_count or
        angular_par.total_par_umol_per_m2_s.len != par_count)
        return error.DiffuseInterceptionNormalizationDimensionMismatch;
    inline for (.{
        angular_par.direct_par_umol_per_m2_s,
        angular_par.diffuse_par_umol_per_m2_s,
        angular_par.total_par_umol_per_m2_s,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDiffuseInterceptionNormalizationInput;
    return species_count;
}

fn partitioned(values: []f64, species_count: usize) Components {
    return .{
        .leaf_shortwave = values[0 * species_count .. 1 * species_count],
        .stalk_shortwave = values[1 * species_count .. 2 * species_count],
        .standing_dead_shortwave = values[2 * species_count .. 3 * species_count],
        .leaf_par = values[3 * species_count .. 4 * species_count],
        .stalk_par = values[4 * species_count .. 5 * species_count],
        .standing_dead_par = values[5 * species_count .. 6 * species_count],
    };
}

test "diffuse normalization scales fields then reconstructs angular PAR" {
    var absorbed = [_]f64{2} ** 6;
    var back = [_]f64{4} ** 6;
    var forward = [_]f64{6} ** 6;
    const direct_par = [_]f64{ 10, 20, 30, 40 };
    var diffuse_par = [_]f64{ 2, 4, 6, 8 };
    var total_par = [_]f64{0} ** 4;
    var state: State = .{
        .diffuse_transmittance_above = 0.8,
        .current_layer_interception_fraction = 0.4,
    };
    const normalization = (try apply(
        0.8,
        2,
        2,
        &state,
        .{
            .absorbed = partitioned(&absorbed, 1),
            .backscattered = partitioned(&back, 1),
            .forward_scattered = partitioned(&forward, 1),
        },
        .{
            .direct_par_umol_per_m2_s = &direct_par,
            .diffuse_par_umol_per_m2_s = &diffuse_par,
            .total_par_umol_per_m2_s = &total_par,
        },
    )).?;
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), normalization, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.diffuse_transmittance_above, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), absorbed[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), diffuse_par[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 11), total_par[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 44), total_par[3], 1e-14);
}

test "nonexcess diffuse interception leaves fields unchanged" {
    var values = [_]f64{42} ** 6;
    const direct_par = [_]f64{1};
    var diffuse_par = [_]f64{2};
    var total_par = [_]f64{3};
    var state: State = .{
        .diffuse_transmittance_above = 0.9,
        .current_layer_interception_fraction = 0.2,
    };
    const result = try apply(0.7, 1, 1, &state, .{
        .absorbed = partitioned(&values, 1),
        .backscattered = partitioned(&values, 1),
        .forward_scattered = partitioned(&values, 1),
    }, .{
        .direct_par_umol_per_m2_s = &direct_par,
        .diffuse_par_umol_per_m2_s = &diffuse_par,
        .total_par_umol_per_m2_s = &total_par,
    });
    try std.testing.expect(result == null);
    try std.testing.expectEqual(@as(f64, 42), values[0]);
    try std.testing.expectEqual(@as(f64, 3), total_par[0]);
}
