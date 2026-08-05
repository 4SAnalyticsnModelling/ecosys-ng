const std = @import("std");

pub const Inputs = struct {
    species_count: usize,
    inclination_count: usize,
    azimuth_count: usize,
    sky_zone_count: usize,
    working_shortwave_megajoules_per_m2_h: f64,
    working_par_umol_per_m2_s: f64,
    leaf_area_m2_per_m2: []const f64,
    stalk_area_m2_per_m2: []const f64,
    standing_dead_area_m2_per_m2: []const f64,
    leaf_clumping: []const f64,
    woody_clumping: f64,
    sky_incidence: []const f64,
    leaf_shortwave_absorptivity: []const f64,
    leaf_par_absorptivity: []const f64,
    stalk_shortwave_absorptivity: f64,
    dead_shortwave_absorptivity: f64,
    stalk_par_absorptivity: f64,
    dead_par_absorptivity: f64,
    leaf_shortwave_transmittance: []const f64,
    leaf_par_transmittance: []const f64,
    inverse_azimuth_area_per_m2: f64,
};

pub const Outputs = struct {
    diffuse_par_umol_per_m2_s: []f64,
    total_par_umol_per_m2_s: []f64,
    forward_shortwave_megajoules_per_m2_h: *f64,
    forward_par_umol_per_m2_s: *f64,
    species_live_shortwave_megajoules_h: []f64,
    species_dead_shortwave_megajoules_h: []f64,
    species_live_par_umol_s: []f64,
    species_dead_par_umol_s: []f64,
    cell_shortwave_megajoules_h: *f64,
    cell_par_umol_s: *f64,
};

/// HOUR1 lines 1693--1733. Preserves species -> inclination -> azimuth ->
/// sky-zone traversal and the subsequent species publication order.
pub fn apply(inputs: Inputs, outputs: Outputs) !void {
    try validate(inputs, outputs);
    for (0..inputs.species_count) |species| {
        var leaf_sw: f64 = 0;
        var stalk_sw: f64 = 0;
        var dead_sw: f64 = 0;
        var leaf_par: f64 = 0;
        var stalk_par: f64 = 0;
        var dead_par: f64 = 0;
        for (0..inputs.inclination_count) |inclination| {
            const surface = species * inputs.inclination_count + inclination;
            const leaf_surface = inputs.leaf_area_m2_per_m2[surface] *
                inputs.leaf_clumping[species];
            const stalk_surface = inputs.stalk_area_m2_per_m2[surface] *
                inputs.woody_clumping;
            const dead_surface = inputs.standing_dead_area_m2_per_m2[surface] *
                inputs.woody_clumping;
            for (0..inputs.azimuth_count) |azimuth| {
                const angle = (species * inputs.inclination_count + inclination) *
                    inputs.azimuth_count + azimuth;
                for (0..inputs.sky_zone_count) |sky| {
                    const incidence = inputs.sky_incidence[
                        (inclination * inputs.azimuth_count + azimuth) *
                            inputs.sky_zone_count + sky
                    ];
                    const leaf_sw_flux = inputs.working_shortwave_megajoules_per_m2_h *
                        incidence * inputs.leaf_shortwave_absorptivity[species];
                    const stalk_sw_flux = inputs.working_shortwave_megajoules_per_m2_h *
                        incidence * inputs.stalk_shortwave_absorptivity;
                    const dead_sw_flux = inputs.working_shortwave_megajoules_per_m2_h *
                        incidence * inputs.dead_shortwave_absorptivity;
                    const leaf_par_flux = inputs.working_par_umol_per_m2_s *
                        incidence * inputs.leaf_par_absorptivity[species];
                    const stalk_par_flux = inputs.working_par_umol_per_m2_s *
                        incidence * inputs.stalk_par_absorptivity;
                    const dead_par_flux = inputs.working_par_umol_per_m2_s *
                        incidence * inputs.dead_par_absorptivity;
                    outputs.diffuse_par_umol_per_m2_s[angle] += leaf_par_flux;
                    outputs.total_par_umol_per_m2_s[angle] += leaf_par_flux;
                    leaf_sw += leaf_surface * leaf_sw_flux;
                    stalk_sw += stalk_surface * stalk_sw_flux;
                    dead_sw += dead_surface * dead_sw_flux;
                    leaf_par += leaf_surface * leaf_par_flux;
                    stalk_par += stalk_surface * stalk_par_flux;
                    dead_par += dead_surface * dead_par_flux;
                }
            }
        }
        outputs.forward_shortwave_megajoules_per_m2_h.* += leaf_sw *
            inputs.leaf_shortwave_transmittance[species] *
            inputs.inverse_azimuth_area_per_m2;
        outputs.forward_par_umol_per_m2_s.* += leaf_par *
            inputs.leaf_par_transmittance[species] *
            inputs.inverse_azimuth_area_per_m2;
        outputs.species_live_shortwave_megajoules_h[species] += leaf_sw + stalk_sw;
        outputs.species_dead_shortwave_megajoules_h[species] += dead_sw;
        outputs.species_live_par_umol_s[species] += leaf_par + stalk_par;
        outputs.species_dead_par_umol_s[species] += dead_par;
        outputs.cell_shortwave_megajoules_h.* += leaf_sw + stalk_sw + dead_sw;
        outputs.cell_par_umol_s.* += leaf_par + stalk_par + dead_par;
    }
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    inline for (.{
        inputs.species_count,
        inputs.inclination_count,
        inputs.azimuth_count,
        inputs.sky_zone_count,
    }) |extent| if (extent == 0) return error.ZeroUpwardScatterExtent;
    const surface_count = try std.math.mul(
        usize,
        inputs.species_count,
        inputs.inclination_count,
    );
    const angle_count = try std.math.mul(usize, surface_count, inputs.azimuth_count);
    const incidence_count = try std.math.mul(
        usize,
        try std.math.mul(usize, inputs.inclination_count, inputs.azimuth_count),
        inputs.sky_zone_count,
    );
    inline for (.{ inputs.leaf_area_m2_per_m2, inputs.stalk_area_m2_per_m2, inputs.standing_dead_area_m2_per_m2 }) |values| if (values.len != surface_count)
        return error.UpwardScatterDimensionMismatch;
    inline for (.{
        inputs.leaf_clumping,
        inputs.leaf_shortwave_absorptivity,
        inputs.leaf_par_absorptivity,
        inputs.leaf_shortwave_transmittance,
        inputs.leaf_par_transmittance,
        outputs.species_live_shortwave_megajoules_h,
        outputs.species_dead_shortwave_megajoules_h,
        outputs.species_live_par_umol_s,
        outputs.species_dead_par_umol_s,
    }) |values| if (values.len != inputs.species_count)
        return error.UpwardScatterDimensionMismatch;
    if (inputs.sky_incidence.len != incidence_count or
        outputs.diffuse_par_umol_per_m2_s.len != angle_count or
        outputs.total_par_umol_per_m2_s.len != angle_count)
        return error.UpwardScatterDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (field.type == f64 and (!std.math.isFinite(value) or value < 0))
            return error.InvalidUpwardScatterInput;
    }
    inline for (.{
        inputs.leaf_area_m2_per_m2,
        inputs.stalk_area_m2_per_m2,
        inputs.standing_dead_area_m2_per_m2,
        inputs.leaf_clumping,
        inputs.sky_incidence,
        inputs.leaf_shortwave_absorptivity,
        inputs.leaf_par_absorptivity,
        inputs.leaf_shortwave_transmittance,
        inputs.leaf_par_transmittance,
        outputs.diffuse_par_umol_per_m2_s,
        outputs.total_par_umol_per_m2_s,
        outputs.species_live_shortwave_megajoules_h,
        outputs.species_dead_shortwave_megajoules_h,
        outputs.species_live_par_umol_s,
        outputs.species_dead_par_umol_s,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidUpwardScatterInput;
    inline for (.{
        outputs.forward_shortwave_megajoules_per_m2_h.*,
        outputs.forward_par_umol_per_m2_s.*,
        outputs.cell_shortwave_megajoules_h.*,
        outputs.cell_par_umol_s.*,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidUpwardScatterInput;
}

test "upward scatter absorbs and publishes one runtime species" {
    var diffuse_par = [_]f64{0};
    var total_par = [_]f64{10};
    var forward_sw: f64 = 0;
    var forward_par: f64 = 0;
    var live_sw = [_]f64{0};
    var dead_sw = [_]f64{0};
    var live_par = [_]f64{0};
    var dead_par = [_]f64{0};
    var cell_sw: f64 = 0;
    var cell_par: f64 = 0;
    try apply(.{
        .species_count = 1,
        .inclination_count = 1,
        .azimuth_count = 1,
        .sky_zone_count = 1,
        .working_shortwave_megajoules_per_m2_h = 2,
        .working_par_umol_per_m2_s = 100,
        .leaf_area_m2_per_m2 = &.{2},
        .stalk_area_m2_per_m2 = &.{1},
        .standing_dead_area_m2_per_m2 = &.{0.5},
        .leaf_clumping = &.{0.5},
        .woody_clumping = 0.5,
        .sky_incidence = &.{0.5},
        .leaf_shortwave_absorptivity = &.{0.5},
        .leaf_par_absorptivity = &.{0.5},
        .stalk_shortwave_absorptivity = 0.5,
        .dead_shortwave_absorptivity = 0.5,
        .stalk_par_absorptivity = 0.5,
        .dead_par_absorptivity = 0.5,
        .leaf_shortwave_transmittance = &.{0.2},
        .leaf_par_transmittance = &.{0.2},
        .inverse_azimuth_area_per_m2 = 0.25,
    }, .{
        .diffuse_par_umol_per_m2_s = &diffuse_par,
        .total_par_umol_per_m2_s = &total_par,
        .forward_shortwave_megajoules_per_m2_h = &forward_sw,
        .forward_par_umol_per_m2_s = &forward_par,
        .species_live_shortwave_megajoules_h = &live_sw,
        .species_dead_shortwave_megajoules_h = &dead_sw,
        .species_live_par_umol_s = &live_par,
        .species_dead_par_umol_s = &dead_par,
        .cell_shortwave_megajoules_h = &cell_sw,
        .cell_par_umol_s = &cell_par,
    });
    try std.testing.expectEqual(@as(f64, 25), diffuse_par[0]);
    try std.testing.expectEqual(@as(f64, 35), total_par[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.025), forward_sw, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.25), forward_par, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), live_sw[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.125), dead_sw[0], 1e-15);
}
