const std = @import("std");

pub const Components = struct {
    leaf_shortwave: []const f64,
    stalk_shortwave: []const f64,
    standing_dead_shortwave: []const f64,
    leaf_par: []const f64,
    stalk_par: []const f64,
    standing_dead_par: []const f64,
};

pub const Inputs = struct {
    direct_absorbed: Components,
    diffuse_absorbed: Components,
    direct_backscattered: Components,
    diffuse_backscattered: Components,
    direct_forward_scattered: Components,
    diffuse_forward_scattered: Components,
    leaf_shortwave_transmittance: []const f64,
    leaf_par_transmittance: []const f64,
    leaf_shortwave_albedo: []const f64,
    leaf_par_albedo: []const f64,
    stalk_shortwave_albedo: f64,
    standing_dead_shortwave_albedo: f64,
    stalk_par_albedo: f64,
    standing_dead_par_albedo: f64,
    inverse_azimuth_area_per_m2: f64,
};

pub const Outputs = struct {
    forward_shortwave_mj_per_m2_h: *f64,
    forward_par_umol_per_m2_s: *f64,
    backscattered_shortwave_mj_per_m2_h: *f64,
    backscattered_par_umol_per_m2_s: *f64,
    species_live_shortwave_mj_h: []f64,
    species_dead_shortwave_mj_h: []f64,
    species_live_par_umol_s: []f64,
    species_dead_par_umol_s: []f64,
    cell_canopy_shortwave_mj_h: *f64,
    cell_canopy_par_umol_s: *f64,
};

/// HOUR1 lines 1541--1578. Preserves per-species component assembly,
/// layer scattered-flux additions, and species/cell diagnostic additions.
pub fn apply(inputs: Inputs, outputs: Outputs) !void {
    const species_count = try validate(inputs, outputs);
    for (0..species_count) |species| {
        const absorbed = sumAt(
            inputs.direct_absorbed,
            inputs.diffuse_absorbed,
            species,
        );
        const backscattered = sumAt(
            inputs.direct_backscattered,
            inputs.diffuse_backscattered,
            species,
        );
        const forward_scattered = sumAt(
            inputs.direct_forward_scattered,
            inputs.diffuse_forward_scattered,
            species,
        );
        outputs.forward_shortwave_mj_per_m2_h.* +=
            (absorbed.leaf_shortwave *
                inputs.leaf_shortwave_transmittance[species] +
                forward_scattered.leaf_shortwave *
                    inputs.leaf_shortwave_albedo[species] +
                forward_scattered.stalk_shortwave *
                    inputs.stalk_shortwave_albedo +
                forward_scattered.standing_dead_shortwave *
                    inputs.standing_dead_shortwave_albedo) *
            inputs.inverse_azimuth_area_per_m2;
        outputs.forward_par_umol_per_m2_s.* +=
            (absorbed.leaf_par * inputs.leaf_par_transmittance[species] +
                forward_scattered.leaf_par *
                    inputs.leaf_par_albedo[species] +
                forward_scattered.stalk_par * inputs.stalk_par_albedo +
                forward_scattered.standing_dead_par *
                    inputs.standing_dead_par_albedo) *
            inputs.inverse_azimuth_area_per_m2;
        outputs.backscattered_shortwave_mj_per_m2_h.* +=
            (backscattered.leaf_shortwave *
                inputs.leaf_shortwave_albedo[species] +
                backscattered.stalk_shortwave *
                    inputs.stalk_shortwave_albedo +
                backscattered.standing_dead_shortwave *
                    inputs.standing_dead_shortwave_albedo) *
            inputs.inverse_azimuth_area_per_m2;
        outputs.backscattered_par_umol_per_m2_s.* +=
            (backscattered.leaf_par * inputs.leaf_par_albedo[species] +
                backscattered.stalk_par * inputs.stalk_par_albedo +
                backscattered.standing_dead_par *
                    inputs.standing_dead_par_albedo) *
            inputs.inverse_azimuth_area_per_m2;
        outputs.species_live_shortwave_mj_h[species] +=
            absorbed.leaf_shortwave + absorbed.stalk_shortwave;
        outputs.species_dead_shortwave_mj_h[species] +=
            absorbed.standing_dead_shortwave;
        outputs.species_live_par_umol_s[species] +=
            absorbed.leaf_par + absorbed.stalk_par;
        outputs.species_dead_par_umol_s[species] += absorbed.standing_dead_par;
        outputs.cell_canopy_shortwave_mj_h.* +=
            absorbed.leaf_shortwave + absorbed.stalk_shortwave +
            absorbed.standing_dead_shortwave;
        outputs.cell_canopy_par_umol_s.* +=
            absorbed.leaf_par + absorbed.stalk_par +
            absorbed.standing_dead_par;
    }
}

fn sumAt(direct: Components, diffuse: Components, species: usize) ComponentValues {
    return .{
        .leaf_shortwave = direct.leaf_shortwave[species] + diffuse.leaf_shortwave[species],
        .stalk_shortwave = direct.stalk_shortwave[species] + diffuse.stalk_shortwave[species],
        .standing_dead_shortwave = direct.standing_dead_shortwave[species] +
            diffuse.standing_dead_shortwave[species],
        .leaf_par = direct.leaf_par[species] + diffuse.leaf_par[species],
        .stalk_par = direct.stalk_par[species] + diffuse.stalk_par[species],
        .standing_dead_par = direct.standing_dead_par[species] + diffuse.standing_dead_par[species],
    };
}

const ComponentValues = struct {
    leaf_shortwave: f64,
    stalk_shortwave: f64,
    standing_dead_shortwave: f64,
    leaf_par: f64,
    stalk_par: f64,
    standing_dead_par: f64,
};

fn validate(inputs: Inputs, outputs: Outputs) !usize {
    const species_count = inputs.direct_absorbed.leaf_shortwave.len;
    if (species_count == 0) return error.ZeroCanopyAggregationSpeciesExtent;
    inline for (.{
        inputs.direct_absorbed,
        inputs.diffuse_absorbed,
        inputs.direct_backscattered,
        inputs.diffuse_backscattered,
        inputs.direct_forward_scattered,
        inputs.diffuse_forward_scattered,
    }) |components| inline for (@typeInfo(Components).@"struct".fields) |field| {
        const values = @field(components, field.name);
        if (values.len != species_count)
            return error.CanopyAggregationDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyAggregationInput;
    };
    inline for (.{
        inputs.leaf_shortwave_transmittance,
        inputs.leaf_par_transmittance,
        inputs.leaf_shortwave_albedo,
        inputs.leaf_par_albedo,
        outputs.species_live_shortwave_mj_h,
        outputs.species_dead_shortwave_mj_h,
        outputs.species_live_par_umol_s,
        outputs.species_dead_par_umol_s,
    }) |values| if (values.len != species_count)
        return error.CanopyAggregationDimensionMismatch;
    inline for (.{
        inputs.leaf_shortwave_transmittance,
        inputs.leaf_par_transmittance,
        inputs.leaf_shortwave_albedo,
        inputs.leaf_par_albedo,
        outputs.species_live_shortwave_mj_h,
        outputs.species_dead_shortwave_mj_h,
        outputs.species_live_par_umol_s,
        outputs.species_dead_par_umol_s,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyAggregationInput;
    inline for (.{
        inputs.stalk_shortwave_albedo,
        inputs.standing_dead_shortwave_albedo,
        inputs.stalk_par_albedo,
        inputs.standing_dead_par_albedo,
        inputs.inverse_azimuth_area_per_m2,
        outputs.forward_shortwave_mj_per_m2_h.*,
        outputs.forward_par_umol_per_m2_s.*,
        outputs.backscattered_shortwave_mj_per_m2_h.*,
        outputs.backscattered_par_umol_per_m2_s.*,
        outputs.cell_canopy_shortwave_mj_h.*,
        outputs.cell_canopy_par_umol_s.*,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyAggregationInput;
    return species_count;
}

fn componentSlices(values: []const f64, species_count: usize) Components {
    return .{
        .leaf_shortwave = values[0 * species_count .. 1 * species_count],
        .stalk_shortwave = values[1 * species_count .. 2 * species_count],
        .standing_dead_shortwave = values[2 * species_count .. 3 * species_count],
        .leaf_par = values[3 * species_count .. 4 * species_count],
        .stalk_par = values[4 * species_count .. 5 * species_count],
        .standing_dead_par = values[5 * species_count .. 6 * species_count],
    };
}

test "species absorbed and scattered aggregation preserves source order" {
    const direct = [_]f64{1} ** 6;
    const diffuse = [_]f64{1} ** 6;
    var forward_sw: f64 = 0;
    var forward_par: f64 = 0;
    var back_sw: f64 = 0;
    var back_par: f64 = 0;
    var live_sw = [_]f64{0};
    var dead_sw = [_]f64{0};
    var live_par = [_]f64{0};
    var dead_par = [_]f64{0};
    var cell_sw: f64 = 0;
    var cell_par: f64 = 0;
    try apply(.{
        .direct_absorbed = componentSlices(&direct, 1),
        .diffuse_absorbed = componentSlices(&diffuse, 1),
        .direct_backscattered = componentSlices(&direct, 1),
        .diffuse_backscattered = componentSlices(&diffuse, 1),
        .direct_forward_scattered = componentSlices(&direct, 1),
        .diffuse_forward_scattered = componentSlices(&diffuse, 1),
        .leaf_shortwave_transmittance = &.{0.1},
        .leaf_par_transmittance = &.{0.1},
        .leaf_shortwave_albedo = &.{0.2},
        .leaf_par_albedo = &.{0.2},
        .stalk_shortwave_albedo = 0.3,
        .standing_dead_shortwave_albedo = 0.4,
        .stalk_par_albedo = 0.3,
        .standing_dead_par_albedo = 0.4,
        .inverse_azimuth_area_per_m2 = 0.25,
    }, .{
        .forward_shortwave_mj_per_m2_h = &forward_sw,
        .forward_par_umol_per_m2_s = &forward_par,
        .backscattered_shortwave_mj_per_m2_h = &back_sw,
        .backscattered_par_umol_per_m2_s = &back_par,
        .species_live_shortwave_mj_h = &live_sw,
        .species_dead_shortwave_mj_h = &dead_sw,
        .species_live_par_umol_s = &live_par,
        .species_dead_par_umol_s = &dead_par,
        .cell_canopy_shortwave_mj_h = &cell_sw,
        .cell_canopy_par_umol_s = &cell_par,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), forward_sw, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.45), back_sw, 1e-15);
    try std.testing.expectEqual(@as(f64, 4), live_sw[0]);
    try std.testing.expectEqual(@as(f64, 2), dead_sw[0]);
    try std.testing.expectEqual(@as(f64, 6), cell_sw);
    try std.testing.expectEqual(@as(f64, 6), cell_par);
}
