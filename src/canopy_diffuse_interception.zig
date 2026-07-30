const std = @import("std");

pub const Components = struct {
    leaf_shortwave: f64 = 0,
    stalk_shortwave: f64 = 0,
    standing_dead_shortwave: f64 = 0,
    leaf_par: f64 = 0,
    stalk_par: f64 = 0,
    standing_dead_par: f64 = 0,
};

pub const Absorptivity = struct {
    leaf_shortwave: f64,
    stalk_shortwave: f64,
    standing_dead_shortwave: f64,
    leaf_par: f64,
    stalk_par: f64,
    standing_dead_par: f64,
};

pub const SurfaceAreas = struct {
    leaf_m2_per_m2: f64,
    stalk_m2_per_m2: f64,
    standing_dead_m2_per_m2: f64,
    leaf_azimuth_per_m2: f64,
    stalk_azimuth_per_m2: f64,
    standing_dead_azimuth_per_m2: f64,
};

pub const Inputs = struct {
    working_diffuse_shortwave_mj_per_m2_h: f64,
    working_diffuse_par_umol_per_m2_s: f64,
    sky_incidence_fraction: []const f64,
    horizontal_sky_incidence_fraction: []const f64,
    is_backscattered: []const bool,
    absorptivity: Absorptivity,
    surfaces: SurfaceAreas,
};

pub const Accumulators = struct {
    diffuse_leaf_par_umol_per_m2_s: f64,
    total_leaf_par_umol_per_m2_s: f64,
    absorbed: Components,
    backscattered: Components,
    forward_scattered: Components,
    intercepted_horizontal_area_fraction: f64,
};

/// HOUR1 lines 1360--1408 for one inclination/azimuth/species/layer.
/// Runtime sky-zone traversal retains all source assignments and additions.
pub fn apply(inputs: Inputs, accumulators: *Accumulators) !void {
    try validate(inputs, accumulators.*);
    for (
        inputs.sky_incidence_fraction,
        inputs.horizontal_sky_incidence_fraction,
        inputs.is_backscattered,
    ) |incidence, horizontal_incidence, is_backscattered| {
        const leaf_shortwave = inputs.working_diffuse_shortwave_mj_per_m2_h *
            incidence * inputs.absorptivity.leaf_shortwave;
        const stalk_shortwave = inputs.working_diffuse_shortwave_mj_per_m2_h *
            incidence * inputs.absorptivity.stalk_shortwave;
        const dead_shortwave = inputs.working_diffuse_shortwave_mj_per_m2_h *
            incidence * inputs.absorptivity.standing_dead_shortwave;
        const leaf_par = inputs.working_diffuse_par_umol_per_m2_s *
            incidence * inputs.absorptivity.leaf_par;
        const stalk_par = inputs.working_diffuse_par_umol_per_m2_s *
            incidence * inputs.absorptivity.stalk_par;
        const dead_par = inputs.working_diffuse_par_umol_per_m2_s *
            incidence * inputs.absorptivity.standing_dead_par;
        accumulators.diffuse_leaf_par_umol_per_m2_s += leaf_par;
        accumulators.total_leaf_par_umol_per_m2_s += leaf_par;
        add(
            &accumulators.absorbed,
            inputs.surfaces,
            leaf_shortwave,
            stalk_shortwave,
            dead_shortwave,
            leaf_par,
            stalk_par,
            dead_par,
        );
        accumulators.intercepted_horizontal_area_fraction +=
            (inputs.surfaces.leaf_azimuth_per_m2 +
                inputs.surfaces.stalk_azimuth_per_m2 +
                inputs.surfaces.standing_dead_azimuth_per_m2) *
            horizontal_incidence;
        if (is_backscattered) {
            add(
                &accumulators.backscattered,
                inputs.surfaces,
                leaf_shortwave,
                stalk_shortwave,
                dead_shortwave,
                leaf_par,
                stalk_par,
                dead_par,
            );
        } else {
            add(
                &accumulators.forward_scattered,
                inputs.surfaces,
                leaf_shortwave,
                stalk_shortwave,
                dead_shortwave,
                leaf_par,
                stalk_par,
                dead_par,
            );
        }
    }
}

fn add(
    target: *Components,
    surfaces: SurfaceAreas,
    leaf_shortwave: f64,
    stalk_shortwave: f64,
    dead_shortwave: f64,
    leaf_par: f64,
    stalk_par: f64,
    dead_par: f64,
) void {
    target.leaf_shortwave += surfaces.leaf_m2_per_m2 * leaf_shortwave;
    target.stalk_shortwave += surfaces.stalk_m2_per_m2 * stalk_shortwave;
    target.standing_dead_shortwave += surfaces.standing_dead_m2_per_m2 * dead_shortwave;
    target.leaf_par += surfaces.leaf_m2_per_m2 * leaf_par;
    target.stalk_par += surfaces.stalk_m2_per_m2 * stalk_par;
    target.standing_dead_par += surfaces.standing_dead_m2_per_m2 * dead_par;
}

fn validate(inputs: Inputs, accumulators: Accumulators) !void {
    const zone_count = inputs.sky_incidence_fraction.len;
    if (zone_count == 0)
        return error.ZeroCanopyDiffuseSkyZoneExtent;
    if (inputs.horizontal_sky_incidence_fraction.len != zone_count or
        inputs.is_backscattered.len != zone_count)
        return error.CanopyDiffuseInterceptionDimensionMismatch;
    inline for (.{
        inputs.working_diffuse_shortwave_mj_per_m2_h,
        inputs.working_diffuse_par_umol_per_m2_s,
        accumulators.diffuse_leaf_par_umol_per_m2_s,
        accumulators.total_leaf_par_umol_per_m2_s,
        accumulators.intercepted_horizontal_area_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyDiffuseInterceptionInput;
    inline for (.{
        inputs.sky_incidence_fraction,
        inputs.horizontal_sky_incidence_fraction,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyDiffuseInterceptionInput;
    inline for (.{ inputs.absorptivity, inputs.surfaces }) |values|
        inline for (@typeInfo(@TypeOf(values)).@"struct".fields) |field| {
            const value = @field(values, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidCanopyDiffuseInterceptionInput;
        };
    inline for (.{ accumulators.absorbed, accumulators.backscattered, accumulators.forward_scattered }) |values|
        inline for (@typeInfo(Components).@"struct".fields) |field|
            if (!std.math.isFinite(@field(values, field.name)))
                return error.InvalidCanopyDiffuseInterceptionInput;
}

test "diffuse sky zones preserve absorption and scattering order" {
    var accumulators: Accumulators = .{
        .diffuse_leaf_par_umol_per_m2_s = 0,
        .total_leaf_par_umol_per_m2_s = 10,
        .absorbed = .{},
        .backscattered = .{},
        .forward_scattered = .{},
        .intercepted_horizontal_area_fraction = 0,
    };
    try apply(.{
        .working_diffuse_shortwave_mj_per_m2_h = 2,
        .working_diffuse_par_umol_per_m2_s = 100,
        .sky_incidence_fraction = &.{ 0.25, 0.75 },
        .horizontal_sky_incidence_fraction = &.{ 0.4, 0.6 },
        .is_backscattered = &.{ true, false },
        .absorptivity = .{
            .leaf_shortwave = 0.5,
            .stalk_shortwave = 0.5,
            .standing_dead_shortwave = 0.5,
            .leaf_par = 0.5,
            .stalk_par = 0.5,
            .standing_dead_par = 0.5,
        },
        .surfaces = .{
            .leaf_m2_per_m2 = 5,
            .stalk_m2_per_m2 = 1,
            .standing_dead_m2_per_m2 = 0.5,
            .leaf_azimuth_per_m2 = 0.5,
            .stalk_azimuth_per_m2 = 0.1,
            .standing_dead_azimuth_per_m2 = 0.05,
        },
    }, &accumulators);
    try std.testing.expectEqual(@as(f64, 50), accumulators.diffuse_leaf_par_umol_per_m2_s);
    try std.testing.expectEqual(@as(f64, 60), accumulators.total_leaf_par_umol_per_m2_s);
    try std.testing.expectEqual(@as(f64, 5), accumulators.absorbed.leaf_shortwave);
    try std.testing.expectEqual(@as(f64, 1.25), accumulators.backscattered.leaf_shortwave);
    try std.testing.expectEqual(@as(f64, 3.75), accumulators.forward_scattered.leaf_shortwave);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.65),
        accumulators.intercepted_horizontal_area_fraction,
        1e-15,
    );
}

test "zone mismatch leaves accumulators unchanged" {
    var accumulators: Accumulators = .{
        .diffuse_leaf_par_umol_per_m2_s = 42,
        .total_leaf_par_umol_per_m2_s = 43,
        .absorbed = .{},
        .backscattered = .{},
        .forward_scattered = .{},
        .intercepted_horizontal_area_fraction = 44,
    };
    try std.testing.expectError(
        error.CanopyDiffuseInterceptionDimensionMismatch,
        apply(.{
            .working_diffuse_shortwave_mj_per_m2_h = 1,
            .working_diffuse_par_umol_per_m2_s = 1,
            .sky_incidence_fraction = &.{ 0.5, 0.5 },
            .horizontal_sky_incidence_fraction = &.{0.5},
            .is_backscattered = &.{ true, false },
            .absorptivity = .{
                .leaf_shortwave = 1,
                .stalk_shortwave = 1,
                .standing_dead_shortwave = 1,
                .leaf_par = 1,
                .stalk_par = 1,
                .standing_dead_par = 1,
            },
            .surfaces = .{
                .leaf_m2_per_m2 = 1,
                .stalk_m2_per_m2 = 1,
                .standing_dead_m2_per_m2 = 1,
                .leaf_azimuth_per_m2 = 1,
                .stalk_azimuth_per_m2 = 1,
                .standing_dead_azimuth_per_m2 = 1,
            },
        }, &accumulators),
    );
    try std.testing.expectEqual(@as(f64, 42), accumulators.diffuse_leaf_par_umol_per_m2_s);
}
