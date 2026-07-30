const std = @import("std");

pub const ScatteringDirection = enum(u8) {
    backward = 1,
    forward = 2,
};

pub const Inputs = struct {
    leaf_inclination_sine: []const f64,
    leaf_inclination_cosine: []const f64,
    sky_azimuth_class_count: usize,
    leaf_azimuth_class_count: usize,
    sky_elevation_rad: f64,
    pi: f64,
    reflected_angle_threshold_rad: f64,
};

/// Immutable, heap-owned canopy interception geometry. Three-dimensional
/// fields preserve Fortran column-major order:
/// `[leaf_azimuth][leaf_inclination][sky_azimuth]`.
pub const Geometry = struct {
    allocator: std.mem.Allocator,
    leaf_inclination_sine: []f64,
    leaf_inclination_cosine: []f64,
    leaf_azimuth_rad: []f64,
    sky_azimuth_rad: []f64,
    sky_azimuth_sine: []f64,
    sky_azimuth_cosine: []f64,
    diffuse_leaf_projection_fraction: []f64,
    diffuse_horizontal_projection_fraction: []f64,
    scattering_direction: []ScatteringDirection,
    sky_sine_sum: f64,

    pub fn deinit(self: *Geometry) void {
        self.allocator.free(self.leaf_inclination_sine);
        self.allocator.free(self.leaf_inclination_cosine);
        self.allocator.free(self.leaf_azimuth_rad);
        self.allocator.free(self.sky_azimuth_rad);
        self.allocator.free(self.sky_azimuth_sine);
        self.allocator.free(self.sky_azimuth_cosine);
        self.allocator.free(self.diffuse_leaf_projection_fraction);
        self.allocator.free(self.diffuse_horizontal_projection_fraction);
        self.allocator.free(self.scattering_direction);
        self.* = undefined;
    }

    pub fn index(
        self: Geometry,
        sky_index: usize,
        inclination_index: usize,
        leaf_azimuth_index: usize,
    ) usize {
        return (leaf_azimuth_index * self.leaf_inclination_sine.len +
            inclination_index) * self.sky_azimuth_rad.len + sky_index;
    }
};

/// Source-order translation of legacy `STARTS` lines 142--181.
pub fn initialize(allocator: std.mem.Allocator, inputs: Inputs) !Geometry {
    const inclination_count = inputs.leaf_inclination_sine.len;
    if (inclination_count == 0 or
        inputs.leaf_inclination_cosine.len != inclination_count or
        inputs.sky_azimuth_class_count == 0 or
        inputs.leaf_azimuth_class_count == 0)
    {
        return error.InvalidCanopyGeometryDimensions;
    }
    inline for (.{
        inputs.sky_elevation_rad,
        inputs.pi,
        inputs.reflected_angle_threshold_rad,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteCanopyGeometry;
    }
    if (inputs.pi <= 0 or
        inputs.sky_elevation_rad <= 0 or
        inputs.sky_elevation_rad >= inputs.pi)
        return error.InvalidCanopyGeometryAngle;
    for (inputs.leaf_inclination_sine, inputs.leaf_inclination_cosine) |
        sine,
        cosine,
    | {
        if (!std.math.isFinite(sine) or !std.math.isFinite(cosine))
            return error.NonFiniteCanopyGeometry;
        if (@abs(sine) > 1 or @abs(cosine) > 1)
            return error.InvalidCanopyGeometryAngle;
    }

    const orientation_count = std.math.mul(
        usize,
        inclination_count,
        inputs.leaf_azimuth_class_count,
    ) catch return error.DimensionOverflow;
    const projection_count = std.math.mul(
        usize,
        orientation_count,
        inputs.sky_azimuth_class_count,
    ) catch return error.DimensionOverflow;

    var geometry: Geometry = .{
        .allocator = allocator,
        .leaf_inclination_sine = try allocator.dupe(
            f64,
            inputs.leaf_inclination_sine,
        ),
        .leaf_inclination_cosine = undefined,
        .leaf_azimuth_rad = undefined,
        .sky_azimuth_rad = undefined,
        .sky_azimuth_sine = undefined,
        .sky_azimuth_cosine = undefined,
        .diffuse_leaf_projection_fraction = undefined,
        .diffuse_horizontal_projection_fraction = undefined,
        .scattering_direction = undefined,
        .sky_sine_sum = 0.0,
    };
    errdefer allocator.free(geometry.leaf_inclination_sine);
    geometry.leaf_inclination_cosine = try allocator.dupe(
        f64,
        inputs.leaf_inclination_cosine,
    );
    errdefer allocator.free(geometry.leaf_inclination_cosine);
    geometry.leaf_azimuth_rad =
        try allocator.alloc(f64, inputs.leaf_azimuth_class_count);
    errdefer allocator.free(geometry.leaf_azimuth_rad);
    geometry.sky_azimuth_rad =
        try allocator.alloc(f64, inputs.sky_azimuth_class_count);
    errdefer allocator.free(geometry.sky_azimuth_rad);
    geometry.sky_azimuth_sine =
        try allocator.alloc(f64, inputs.sky_azimuth_class_count);
    errdefer allocator.free(geometry.sky_azimuth_sine);
    geometry.sky_azimuth_cosine =
        try allocator.alloc(f64, inputs.sky_azimuth_class_count);
    errdefer allocator.free(geometry.sky_azimuth_cosine);
    geometry.diffuse_leaf_projection_fraction =
        try allocator.alloc(f64, projection_count);
    errdefer allocator.free(geometry.diffuse_leaf_projection_fraction);
    geometry.diffuse_horizontal_projection_fraction =
        try allocator.alloc(f64, projection_count);
    errdefer allocator.free(geometry.diffuse_horizontal_projection_fraction);
    geometry.scattering_direction =
        try allocator.alloc(ScatteringDirection, projection_count);
    errdefer allocator.free(geometry.scattering_direction);

    const leaf_azimuth_count_f64: f64 =
        @floatFromInt(inputs.leaf_azimuth_class_count);
    for (0..inputs.leaf_azimuth_class_count) |leaf_azimuth_index| {
        geometry.leaf_azimuth_rad[leaf_azimuth_index] =
            (@as(f64, @floatFromInt(leaf_azimuth_index)) + 0.5) *
            inputs.pi / leaf_azimuth_count_f64;
    }

    const sky_count_f64: f64 = @floatFromInt(inputs.sky_azimuth_class_count);
    for (0..inputs.sky_azimuth_class_count) |sky_index| {
        geometry.sky_azimuth_rad[sky_index] =
            inputs.pi *
            (2.0 * @as(f64, @floatFromInt(sky_index)) + 1.0) /
            sky_count_f64;
        geometry.sky_azimuth_sine[sky_index] =
            @sin(inputs.sky_elevation_rad);
        geometry.sky_azimuth_cosine[sky_index] =
            @cos(inputs.sky_elevation_rad);
        geometry.sky_sine_sum += geometry.sky_azimuth_sine[sky_index];

        for (0..inputs.leaf_azimuth_class_count) |leaf_azimuth_index| {
            const azimuth_cosine = @cos(
                geometry.leaf_azimuth_rad[leaf_azimuth_index] -
                    geometry.sky_azimuth_rad[sky_index],
            );
            for (0..inclination_count) |inclination_index| {
                const signed_projection =
                    geometry.leaf_inclination_cosine[inclination_index] *
                    geometry.sky_azimuth_sine[sky_index] +
                    geometry.leaf_inclination_sine[inclination_index] *
                        geometry.sky_azimuth_cosine[sky_index] *
                        azimuth_cosine;
                if (signed_projection < -1 or signed_projection > 1)
                    return error.InvalidCanopyProjection;
                const index = geometry.index(
                    sky_index,
                    inclination_index,
                    leaf_azimuth_index,
                );
                geometry.diffuse_leaf_projection_fraction[index] =
                    @abs(signed_projection);
                geometry.diffuse_horizontal_projection_fraction[index] =
                    geometry.diffuse_leaf_projection_fraction[index] /
                    geometry.sky_azimuth_sine[sky_index];

                const reflected_angle_rad =
                    if (geometry.leaf_inclination_cosine[inclination_index] >
                    geometry.sky_azimuth_sine[sky_index])
                        std.math.acos(signed_projection)
                    else
                        -std.math.acos(signed_projection);
                const outgoing_angle_rad =
                    if (reflected_angle_rad >
                    inputs.reflected_angle_threshold_rad)
                        inputs.sky_elevation_rad + 2.0 * reflected_angle_rad
                    else
                        inputs.sky_elevation_rad -
                            2.0 * (inputs.pi + reflected_angle_rad);
                geometry.scattering_direction[index] =
                    if (outgoing_angle_rad > 0 and
                    outgoing_angle_rad < inputs.pi)
                        .backward
                    else
                        .forward;
            }
        }
    }
    return geometry;
}

test "STARTS four-class geometry reproduces source constants and extents" {
    var geometry = try initialize(std.testing.allocator, .{
        .leaf_inclination_sine = &.{ 0.195, 0.556, 0.831, 0.981 },
        .leaf_inclination_cosine = &.{ 0.981, 0.831, 0.556, 0.195 },
        .sky_azimuth_class_count = 4,
        .leaf_azimuth_class_count = 4,
        .sky_elevation_rad = 3.1416 / 4.0,
        .pi = 3.1416,
        .reflected_angle_threshold_rad = -1.5708,
    });
    defer geometry.deinit();

    try std.testing.expectEqual(@as(usize, 64), geometry.scattering_direction.len);
    try std.testing.expectApproxEqAbs(
        @as(f64, 3.1416 / 8.0),
        geometry.leaf_azimuth_rad[0],
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        4.0 * @sin(3.1416 / 4.0),
        geometry.sky_sine_sum,
        1.0e-15,
    );
    const index = geometry.index(0, 0, 0);
    const signed_projection =
        0.981 * @sin(3.1416 / 4.0) +
        0.195 * @cos(3.1416 / 4.0) *
            @cos(3.1416 / 8.0 - 3.1416 / 4.0);
    try std.testing.expectApproxEqAbs(
        @abs(signed_projection),
        geometry.diffuse_leaf_projection_fraction[index],
        1.0e-15,
    );
}

test "runtime class counts allocate only requested geometry" {
    var geometry = try initialize(std.testing.allocator, .{
        .leaf_inclination_sine = &.{ 0.5, 0.8 },
        .leaf_inclination_cosine = &.{ 0.866025403784, 0.6 },
        .sky_azimuth_class_count = 3,
        .leaf_azimuth_class_count = 5,
        .sky_elevation_rad = std.math.pi / 4.0,
        .pi = std.math.pi,
        .reflected_angle_threshold_rad = -std.math.pi / 2.0,
    });
    defer geometry.deinit();

    try std.testing.expectEqual(@as(usize, 30), geometry.diffuse_leaf_projection_fraction.len);
    for (geometry.diffuse_leaf_projection_fraction) |fraction| {
        try std.testing.expect(std.math.isFinite(fraction));
        try std.testing.expect(fraction >= 0 and fraction <= 1);
    }
}
