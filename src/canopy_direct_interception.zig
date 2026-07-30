const std = @import("std");

pub const Components = struct {
    leaf_shortwave: f64 = 0,
    stalk_shortwave: f64 = 0,
    standing_dead_shortwave: f64 = 0,
    leaf_par: f64 = 0,
    stalk_par: f64 = 0,
    standing_dead_par: f64 = 0,
};

pub const DirectFields = struct {
    leaf_shortwave: []const f64,
    stalk_shortwave: []const f64,
    standing_dead_shortwave: []const f64,
    leaf_par: []const f64,
    stalk_par: []const f64,
    standing_dead_par: []const f64,
};

pub const Inputs = struct {
    inclination_count: usize,
    azimuth_count: usize,
    leaf_area_by_inclination_m2_per_m2: []const f64,
    stalk_area_by_inclination_m2_per_m2: []const f64,
    standing_dead_area_by_inclination_m2_per_m2: []const f64,
    leaf_clumping_factor: f64,
    woody_clumping_factor: f64,
    inverse_cell_area_per_m2: f64,
    inverse_azimuth_area_per_m2: f64,
    direct_transmittance_above: f64,
    horizontal_relative_incidence: []const f64,
    /// Source IALBS: one is backscatter; every other value is forward.
    is_backscattered: []const bool,
    direct_fields: DirectFields,
};

pub const Accumulators = struct {
    absorbed: Components,
    backscattered: Components,
    forward_scattered: Components,
    intercepted_horizontal_area_fraction: f64,
};

/// HOUR1 lines 1293--1345 for one species and layer. Preserves inclination
/// -> azimuth traversal and every accumulator addition in source order.
pub fn apply(inputs: Inputs, accumulators: *Accumulators) !void {
    try validate(inputs, accumulators.*);
    for (0..inputs.inclination_count) |inclination| {
        const leaf_unshaded_m2_per_m2 =
            inputs.leaf_area_by_inclination_m2_per_m2[inclination] *
            inputs.leaf_clumping_factor;
        const leaf_azimuth_area_per_m2 =
            leaf_unshaded_m2_per_m2 * inputs.inverse_azimuth_area_per_m2;
        _ = leaf_azimuth_area_per_m2;
        const leaf_sunlit_m2_per_m2 =
            leaf_unshaded_m2_per_m2 * inputs.direct_transmittance_above;
        const leaf_sunlit_horizontal_fraction =
            leaf_sunlit_m2_per_m2 * inputs.inverse_cell_area_per_m2;
        const stalk_unshaded_m2_per_m2 =
            inputs.stalk_area_by_inclination_m2_per_m2[inclination] *
            inputs.woody_clumping_factor;
        const stalk_azimuth_area_per_m2 =
            stalk_unshaded_m2_per_m2 * inputs.inverse_azimuth_area_per_m2;
        _ = stalk_azimuth_area_per_m2;
        const stalk_sunlit_m2_per_m2 =
            stalk_unshaded_m2_per_m2 * inputs.direct_transmittance_above;
        const stalk_sunlit_horizontal_fraction =
            stalk_sunlit_m2_per_m2 * inputs.inverse_cell_area_per_m2;
        const dead_unshaded_m2_per_m2 =
            inputs.standing_dead_area_by_inclination_m2_per_m2[inclination] *
            inputs.woody_clumping_factor;
        const dead_azimuth_area_per_m2 =
            dead_unshaded_m2_per_m2 * inputs.inverse_azimuth_area_per_m2;
        _ = dead_azimuth_area_per_m2;
        const dead_sunlit_m2_per_m2 =
            dead_unshaded_m2_per_m2 * inputs.direct_transmittance_above;
        const dead_sunlit_horizontal_fraction =
            dead_sunlit_m2_per_m2 * inputs.inverse_cell_area_per_m2;

        for (0..inputs.azimuth_count) |azimuth| {
            const angle = azimuth * inputs.inclination_count + inclination;
            addScaled(
                &accumulators.absorbed,
                inputs.direct_fields,
                angle,
                leaf_sunlit_m2_per_m2,
                stalk_sunlit_m2_per_m2,
                dead_sunlit_m2_per_m2,
            );
            accumulators.intercepted_horizontal_area_fraction +=
                (leaf_sunlit_horizontal_fraction +
                    stalk_sunlit_horizontal_fraction +
                    dead_sunlit_horizontal_fraction) *
                inputs.horizontal_relative_incidence[angle];
            if (inputs.is_backscattered[angle]) {
                addScaled(
                    &accumulators.backscattered,
                    inputs.direct_fields,
                    angle,
                    leaf_sunlit_m2_per_m2,
                    stalk_sunlit_m2_per_m2,
                    dead_sunlit_m2_per_m2,
                );
            } else {
                addScaled(
                    &accumulators.forward_scattered,
                    inputs.direct_fields,
                    angle,
                    leaf_sunlit_m2_per_m2,
                    stalk_sunlit_m2_per_m2,
                    dead_sunlit_m2_per_m2,
                );
            }
        }
    }
}

fn addScaled(
    target: *Components,
    fields: DirectFields,
    angle: usize,
    leaf_area_m2: f64,
    stalk_area_m2: f64,
    dead_area_m2: f64,
) void {
    target.leaf_shortwave += leaf_area_m2 * fields.leaf_shortwave[angle];
    target.stalk_shortwave += stalk_area_m2 * fields.stalk_shortwave[angle];
    target.standing_dead_shortwave +=
        dead_area_m2 * fields.standing_dead_shortwave[angle];
    target.leaf_par += leaf_area_m2 * fields.leaf_par[angle];
    target.stalk_par += stalk_area_m2 * fields.stalk_par[angle];
    target.standing_dead_par += dead_area_m2 * fields.standing_dead_par[angle];
}

fn validate(inputs: Inputs, accumulators: Accumulators) !void {
    if (inputs.inclination_count == 0 or inputs.azimuth_count == 0)
        return error.ZeroCanopyDirectInterceptionExtent;
    const angle_count = try std.math.mul(
        usize,
        inputs.inclination_count,
        inputs.azimuth_count,
    );
    inline for (.{
        inputs.leaf_area_by_inclination_m2_per_m2,
        inputs.stalk_area_by_inclination_m2_per_m2,
        inputs.standing_dead_area_by_inclination_m2_per_m2,
    }) |values| if (values.len != inputs.inclination_count)
        return error.CanopyDirectInterceptionDimensionMismatch;
    if (inputs.horizontal_relative_incidence.len != angle_count or
        inputs.is_backscattered.len != angle_count)
        return error.CanopyDirectInterceptionDimensionMismatch;
    inline for (@typeInfo(DirectFields).@"struct".fields) |field|
        if (@field(inputs.direct_fields, field.name).len != angle_count)
            return error.CanopyDirectInterceptionDimensionMismatch;
    inline for (.{
        inputs.leaf_clumping_factor,
        inputs.woody_clumping_factor,
        inputs.inverse_cell_area_per_m2,
        inputs.inverse_azimuth_area_per_m2,
        inputs.direct_transmittance_above,
        accumulators.intercepted_horizontal_area_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyDirectInterceptionInput;
    inline for (.{
        inputs.leaf_area_by_inclination_m2_per_m2,
        inputs.stalk_area_by_inclination_m2_per_m2,
        inputs.standing_dead_area_by_inclination_m2_per_m2,
        inputs.horizontal_relative_incidence,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyDirectInterceptionInput;
    inline for (@typeInfo(DirectFields).@"struct".fields) |field|
        for (@field(inputs.direct_fields, field.name)) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidCanopyDirectInterceptionInput;
    inline for (.{ accumulators.absorbed, accumulators.backscattered, accumulators.forward_scattered }) |components|
        inline for (@typeInfo(Components).@"struct".fields) |field|
            if (!std.math.isFinite(@field(components, field.name)))
                return error.InvalidCanopyDirectInterceptionInput;
}

test "surface scaling and direct interception preserve N then M order" {
    const fields: DirectFields = .{
        .leaf_shortwave = &.{ 2, 2 },
        .stalk_shortwave = &.{ 3, 3 },
        .standing_dead_shortwave = &.{ 4, 4 },
        .leaf_par = &.{ 20, 20 },
        .stalk_par = &.{ 30, 30 },
        .standing_dead_par = &.{ 40, 40 },
    };
    var accumulators: Accumulators = .{
        .absorbed = .{},
        .backscattered = .{},
        .forward_scattered = .{},
        .intercepted_horizontal_area_fraction = 0,
    };
    try apply(.{
        .inclination_count = 1,
        .azimuth_count = 2,
        .leaf_area_by_inclination_m2_per_m2 = &.{10},
        .stalk_area_by_inclination_m2_per_m2 = &.{2},
        .standing_dead_area_by_inclination_m2_per_m2 = &.{1},
        .leaf_clumping_factor = 0.5,
        .woody_clumping_factor = 0.5,
        .inverse_cell_area_per_m2 = 0.1,
        .inverse_azimuth_area_per_m2 = 0.25,
        .direct_transmittance_above = 0.8,
        .horizontal_relative_incidence = &.{ 1, 1 },
        .is_backscattered = &.{ true, false },
        .direct_fields = fields,
    }, &accumulators);
    try std.testing.expectEqual(@as(f64, 16), accumulators.absorbed.leaf_shortwave);
    try std.testing.expectEqual(@as(f64, 8), accumulators.backscattered.leaf_shortwave);
    try std.testing.expectEqual(@as(f64, 8), accumulators.forward_scattered.leaf_shortwave);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.04),
        accumulators.intercepted_horizontal_area_fraction,
        1e-15,
    );
}

test "dimension mismatch leaves accumulators unchanged" {
    var accumulators: Accumulators = .{
        .absorbed = .{ .leaf_shortwave = 42 },
        .backscattered = .{},
        .forward_scattered = .{},
        .intercepted_horizontal_area_fraction = 43,
    };
    const one = [_]f64{1};
    try std.testing.expectError(
        error.CanopyDirectInterceptionDimensionMismatch,
        apply(.{
            .inclination_count = 1,
            .azimuth_count = 2,
            .leaf_area_by_inclination_m2_per_m2 = &one,
            .stalk_area_by_inclination_m2_per_m2 = &one,
            .standing_dead_area_by_inclination_m2_per_m2 = &one,
            .leaf_clumping_factor = 1,
            .woody_clumping_factor = 1,
            .inverse_cell_area_per_m2 = 1,
            .inverse_azimuth_area_per_m2 = 1,
            .direct_transmittance_above = 1,
            .horizontal_relative_incidence = &one,
            .is_backscattered = &.{true},
            .direct_fields = .{
                .leaf_shortwave = &one,
                .stalk_shortwave = &one,
                .standing_dead_shortwave = &one,
                .leaf_par = &one,
                .stalk_par = &one,
                .standing_dead_par = &one,
            },
        }, &accumulators),
    );
    try std.testing.expectEqual(@as(f64, 42), accumulators.absorbed.leaf_shortwave);
    try std.testing.expectEqual(@as(f64, 43), accumulators.intercepted_horizontal_area_fraction);
}
