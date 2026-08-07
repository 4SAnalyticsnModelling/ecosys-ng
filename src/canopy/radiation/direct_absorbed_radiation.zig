const std = @import("std");

pub const SpeciesRadiation = struct {
    leaf_shortwave_megajoules_m2_h: []const f64,
    stalk_shortwave_megajoules_m2_h: []const f64,
    standing_dead_shortwave_megajoules_m2_h: []const f64,
    leaf_par_umol_m2_s: []const f64,
    stalk_par_umol_m2_s: []const f64,
    standing_dead_par_umol_m2_s: []const f64,
};

pub const Outputs = struct {
    leaf_shortwave_megajoules_m2_h: []f64,
    stalk_shortwave_megajoules_m2_h: []f64,
    standing_dead_shortwave_megajoules_m2_h: []f64,
    leaf_par_umol_m2_s: []f64,
    stalk_par_umol_m2_s: []f64,
    standing_dead_par_umol_m2_s: []f64,
    diffuse_leaf_par_umol_m2_s: []f64,
    total_leaf_par_umol_m2_s: []f64,
};

/// `hour1.f` lines 1113--1127. `incidence` is azimuth-major/inclination-minor;
/// direct outputs add species, and layer outputs add a final layer index.
pub fn apply(
    incidence: []const f64,
    species_radiation: SpeciesRadiation,
    layer_count: usize,
    outputs: Outputs,
) !void {
    const species_count = try validate(
        incidence,
        species_radiation,
        layer_count,
        outputs,
    );
    for (incidence, 0..) |signed_incidence, angle| {
        const projected_incidence = @abs(signed_incidence);
        for (0..species_count) |species| {
            const direct_index = angle * species_count + species;
            outputs.leaf_shortwave_megajoules_m2_h[direct_index] =
                species_radiation.leaf_shortwave_megajoules_m2_h[species] *
                projected_incidence;
            outputs.stalk_shortwave_megajoules_m2_h[direct_index] =
                species_radiation.stalk_shortwave_megajoules_m2_h[species] *
                projected_incidence;
            outputs.standing_dead_shortwave_megajoules_m2_h[direct_index] =
                species_radiation.standing_dead_shortwave_megajoules_m2_h[species] *
                projected_incidence;
            outputs.leaf_par_umol_m2_s[direct_index] =
                species_radiation.leaf_par_umol_m2_s[species] *
                projected_incidence;
            outputs.stalk_par_umol_m2_s[direct_index] =
                species_radiation.stalk_par_umol_m2_s[species] *
                projected_incidence;
            outputs.standing_dead_par_umol_m2_s[direct_index] =
                species_radiation.standing_dead_par_umol_m2_s[species] *
                projected_incidence;
            for (0..layer_count) |layer| {
                const layer_index = direct_index * layer_count + layer;
                outputs.diffuse_leaf_par_umol_m2_s[layer_index] = 0;
                outputs.total_leaf_par_umol_m2_s[layer_index] =
                    outputs.leaf_par_umol_m2_s[direct_index];
            }
        }
    }
}

fn validate(
    incidence: []const f64,
    species_radiation: SpeciesRadiation,
    layer_count: usize,
    outputs: Outputs,
) !usize {
    const species_count = species_radiation.leaf_shortwave_megajoules_m2_h.len;
    if (incidence.len == 0 or species_count == 0 or layer_count == 0)
        return error.ZeroCanopyDirectRadiationExtent;
    inline for (@typeInfo(SpeciesRadiation).@"struct".fields) |field| {
        const values = @field(species_radiation, field.name);
        if (values.len != species_count)
            return error.CanopyDirectRadiationDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyDirectRadiationInput;
    }
    for (incidence) |value| if (!std.math.isFinite(value))
        return error.InvalidCanopyDirectRadiationInput;
    const direct_count = try std.math.mul(usize, incidence.len, species_count);
    const layer_value_count = try std.math.mul(usize, direct_count, layer_count);
    inline for (.{
        outputs.leaf_shortwave_megajoules_m2_h,
        outputs.stalk_shortwave_megajoules_m2_h,
        outputs.standing_dead_shortwave_megajoules_m2_h,
        outputs.leaf_par_umol_m2_s,
        outputs.stalk_par_umol_m2_s,
        outputs.standing_dead_par_umol_m2_s,
    }) |values| if (values.len != direct_count)
        return error.CanopyDirectRadiationDimensionMismatch;
    if (outputs.diffuse_leaf_par_umol_m2_s.len != layer_value_count or
        outputs.total_leaf_par_umol_m2_s.len != layer_value_count)
        return error.CanopyDirectRadiationDimensionMismatch;
    return species_count;
}

test "angle species layer traversal initializes absorbed radiation" {
    var leaf_sw: [4]f64 = undefined;
    var stalk_sw: [4]f64 = undefined;
    var dead_sw: [4]f64 = undefined;
    var leaf_par: [4]f64 = undefined;
    var stalk_par: [4]f64 = undefined;
    var dead_par: [4]f64 = undefined;
    var diffuse_par: [12]f64 = undefined;
    var total_par: [12]f64 = undefined;
    try apply(&.{ 0.5, -0.25 }, .{
        .leaf_shortwave_megajoules_m2_h = &.{ 2, 4 },
        .stalk_shortwave_megajoules_m2_h = &.{ 1, 2 },
        .standing_dead_shortwave_megajoules_m2_h = &.{ 0.5, 1 },
        .leaf_par_umol_m2_s = &.{ 100, 200 },
        .stalk_par_umol_m2_s = &.{ 50, 100 },
        .standing_dead_par_umol_m2_s = &.{ 25, 50 },
    }, 3, .{
        .leaf_shortwave_megajoules_m2_h = &leaf_sw,
        .stalk_shortwave_megajoules_m2_h = &stalk_sw,
        .standing_dead_shortwave_megajoules_m2_h = &dead_sw,
        .leaf_par_umol_m2_s = &leaf_par,
        .stalk_par_umol_m2_s = &stalk_par,
        .standing_dead_par_umol_m2_s = &dead_par,
        .diffuse_leaf_par_umol_m2_s = &diffuse_par,
        .total_leaf_par_umol_m2_s = &total_par,
    });
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 0.5, 1 }, &leaf_sw);
    try std.testing.expectEqualSlices(f64, &.{ 50, 50, 50 }, total_par[0..3]);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, diffuse_par[9..12]);
}

test "dimension mismatch leaves outputs unchanged" {
    var one = [_]f64{42};
    try std.testing.expectError(error.CanopyDirectRadiationDimensionMismatch, apply(
        &.{0.5},
        .{
            .leaf_shortwave_megajoules_m2_h = &.{1},
            .stalk_shortwave_megajoules_m2_h = &.{1},
            .standing_dead_shortwave_megajoules_m2_h = &.{1},
            .leaf_par_umol_m2_s = &.{1},
            .stalk_par_umol_m2_s = &.{1},
            .standing_dead_par_umol_m2_s = &.{1},
        },
        2,
        .{
            .leaf_shortwave_megajoules_m2_h = &one,
            .stalk_shortwave_megajoules_m2_h = &one,
            .standing_dead_shortwave_megajoules_m2_h = &one,
            .leaf_par_umol_m2_s = &one,
            .stalk_par_umol_m2_s = &one,
            .standing_dead_par_umol_m2_s = &one,
            .diffuse_leaf_par_umol_m2_s = &one,
            .total_leaf_par_umol_m2_s = &one,
        },
    ));
    try std.testing.expectEqual(@as(f64, 42), one[0]);
}
