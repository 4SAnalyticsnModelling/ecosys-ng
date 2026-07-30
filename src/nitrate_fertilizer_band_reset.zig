const std = @import("std");

pub const PairedPool = struct {
    non_band: []f64,
    band: []f64,
};

pub const State = struct {
    active: bool,
    row_width_m: f64,
    band_center_depth_m: f64,
    band_upper_depth_m: f64,
    band_depth_m: []f64,
    band_width_m: []f64,
    band_volume_fraction: []f64,
    non_band_volume_fraction: []f64,
    nitrate_g_n: PairedPool,
    nitrite_g_n: PairedPool,
};

pub const Inputs = struct {
    target_layer: usize,
    application_depth_m: f64,
    layer_upper_depth_m: []const f64,
    layer_thickness_m: []const f64,
    requested_row_width_m: f64,
    new_banded_ammonium_g_n_per_m2: f64,
    new_banded_ammonia_g_n_per_m2: f64,
    new_banded_urea_g_n_per_m2: f64,
    new_banded_nitrate_g_n_per_m2: f64,
    maximum_initial_band_dimension_m: f64 = 0.025,
    minimum_active_layer_thickness_m: f64,
    maximum_band_volume_fraction: f64 = 0.9999,
};

/// Exact NO3 band reset and conservative repartition from hour1.f:356-388.
pub fn apply(state: *State, inputs: Inputs) !bool {
    const count = inputs.layer_thickness_m.len;
    if (count == 0 or inputs.layer_upper_depth_m.len != count or
        inputs.target_layer >= count)
        return error.NitrateBandDimensionMismatch;
    inline for (.{
        state.band_depth_m,
        state.band_width_m,
        state.band_volume_fraction,
        state.non_band_volume_fraction,
        state.nitrate_g_n.non_band,
        state.nitrate_g_n.band,
        state.nitrite_g_n.non_band,
        state.nitrite_g_n.band,
    }) |values| if (values.len != count)
        return error.NitrateBandDimensionMismatch;
    inline for (.{
        inputs.application_depth_m,
        inputs.requested_row_width_m,
        inputs.new_banded_ammonium_g_n_per_m2,
        inputs.new_banded_ammonia_g_n_per_m2,
        inputs.new_banded_urea_g_n_per_m2,
        inputs.new_banded_nitrate_g_n_per_m2,
        inputs.maximum_initial_band_dimension_m,
        inputs.minimum_active_layer_thickness_m,
        inputs.maximum_band_volume_fraction,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteNitrateBandInput;
    if (inputs.application_depth_m < 0 or inputs.requested_row_width_m <= 0 or
        inputs.new_banded_ammonium_g_n_per_m2 < 0 or
        inputs.new_banded_ammonia_g_n_per_m2 < 0 or
        inputs.new_banded_urea_g_n_per_m2 < 0 or
        inputs.new_banded_nitrate_g_n_per_m2 < 0 or
        inputs.maximum_initial_band_dimension_m <= 0 or
        inputs.minimum_active_layer_thickness_m < 0 or
        inputs.maximum_band_volume_fraction <= 0 or
        inputs.maximum_band_volume_fraction >= 1)
        return error.InvalidNitrateBandInput;
    inline for (.{ state.nitrate_g_n, state.nitrite_g_n }) |pool|
        for (pool.non_band, pool.band) |non_band, band| {
            if (!std.math.isFinite(non_band) or non_band < 0 or
                !std.math.isFinite(band) or band < 0 or
                !std.math.isFinite(non_band + band))
                return error.InvalidNitrateBandState;
        };
    for (inputs.layer_upper_depth_m, inputs.layer_thickness_m) |upper, thickness|
        if (!std.math.isFinite(upper) or upper < 0 or
            !std.math.isFinite(thickness) or thickness <= 0)
            return error.InvalidNitrateBandLayerGeometry;

    const has_new_band =
        inputs.new_banded_ammonium_g_n_per_m2 +
        inputs.new_banded_ammonia_g_n_per_m2 +
        inputs.new_banded_urea_g_n_per_m2 +
        inputs.new_banded_nitrate_g_n_per_m2 > 0;
    const has_stranded_target_band =
        (state.nitrate_g_n.band[inputs.target_layer] > 0 or
            state.nitrite_g_n.band[inputs.target_layer] > 0) and
        !state.active;
    if (!has_new_band and !has_stranded_target_band) return false;

    const target_depth = @min(
        inputs.maximum_initial_band_dimension_m,
        inputs.application_depth_m -
            inputs.layer_upper_depth_m[inputs.target_layer],
    );
    if (!std.math.isFinite(target_depth) or target_depth < 0)
        return error.InvalidNitrateBandTargetDepth;
    const target_width = @min(
        inputs.maximum_initial_band_dimension_m,
        inputs.requested_row_width_m,
    );
    for (0..count) |layer| {
        const depth = if (layer == inputs.target_layer) target_depth else 0;
        const width = if (layer == inputs.target_layer) target_width else 0;
        const fraction = if (inputs.layer_thickness_m[layer] >
            inputs.minimum_active_layer_thickness_m)
            @min(
                inputs.maximum_band_volume_fraction,
                width / inputs.requested_row_width_m *
                    depth / inputs.layer_thickness_m[layer],
            )
        else
            0;
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction >= 1)
            return error.InvalidNitrateBandFraction;
    }

    state.active = true;
    state.row_width_m = inputs.requested_row_width_m;
    state.band_center_depth_m = inputs.application_depth_m;
    state.band_upper_depth_m =
        inputs.application_depth_m -
        inputs.maximum_initial_band_dimension_m;
    for (0..count) |layer| {
        const depth = if (layer == inputs.target_layer) target_depth else 0;
        const width = if (layer == inputs.target_layer) target_width else 0;
        const fraction = if (inputs.layer_thickness_m[layer] >
            inputs.minimum_active_layer_thickness_m)
            @min(
                inputs.maximum_band_volume_fraction,
                width / inputs.requested_row_width_m *
                    depth / inputs.layer_thickness_m[layer],
            )
        else
            0;
        state.band_depth_m[layer] = depth;
        state.band_width_m[layer] = width;
        state.band_volume_fraction[layer] = fraction;
        state.non_band_volume_fraction[layer] = 1 - fraction;
        inline for (.{ &state.nitrate_g_n, &state.nitrite_g_n }) |pool| {
            const total = pool.non_band[layer] + pool.band[layer];
            pool.non_band[layer] = total * (1 - fraction);
            pool.band[layer] = total * fraction;
        }
    }
    return true;
}

fn makeState(
    geometry: *[2]f64,
    band_fraction: *[2]f64,
    non_band_fraction: *[2]f64,
    nitrate_non: *[2]f64,
    nitrate_band: *[2]f64,
    nitrite_non: *[2]f64,
    nitrite_band: *[2]f64,
) State {
    return .{
        .active = false,
        .row_width_m = 0,
        .band_center_depth_m = 0,
        .band_upper_depth_m = 0,
        .band_depth_m = geometry,
        .band_width_m = geometry,
        .band_volume_fraction = band_fraction,
        .non_band_volume_fraction = non_band_fraction,
        .nitrate_g_n = .{ .non_band = nitrate_non, .band = nitrate_band },
        .nitrite_g_n = .{ .non_band = nitrite_non, .band = nitrite_band },
    };
}

test "any banded nitrogen initializes nitrate geometry and conserves pools" {
    var geometry = [_]f64{ 9, 9 };
    var band_fraction = [_]f64{ 9, 9 };
    var non_band_fraction = [_]f64{ 9, 9 };
    var nitrate_non = [_]f64{ 10, 18 };
    var nitrate_band = [_]f64{ 0, 2 };
    var nitrite_non = [_]f64{ 5, 8 };
    var nitrite_band = [_]f64{ 0, 2 };
    var state = makeState(&geometry, &band_fraction, &non_band_fraction, &nitrate_non, &nitrate_band, &nitrite_non, &nitrite_band);
    try std.testing.expect(try apply(&state, .{
        .target_layer = 1,
        .application_depth_m = 0.12,
        .layer_upper_depth_m = &.{ 0, 0.1 },
        .layer_thickness_m = &.{ 0.1, 0.2 },
        .requested_row_width_m = 0.05,
        .new_banded_ammonium_g_n_per_m2 = 1,
        .new_banded_ammonia_g_n_per_m2 = 0,
        .new_banded_urea_g_n_per_m2 = 0,
        .new_banded_nitrate_g_n_per_m2 = 0,
        .minimum_active_layer_thickness_m = 0.001,
    }));
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), band_fraction[1], 1e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, 20),
        nitrate_non[1] + nitrate_band[1],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 1), nitrate_band[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), nitrite_band[1], 1e-15);
}

test "no new or stranded nitrate band is a no-op" {
    var geometry = [_]f64{ 0, 0 };
    var band_fraction = [_]f64{ 0, 0 };
    var non_band_fraction = [_]f64{ 1, 1 };
    var nitrate_non = [_]f64{ 1, 1 };
    var nitrate_band = [_]f64{ 0, 0 };
    var nitrite_non = [_]f64{ 1, 1 };
    var nitrite_band = [_]f64{ 0, 0 };
    var state = makeState(&geometry, &band_fraction, &non_band_fraction, &nitrate_non, &nitrate_band, &nitrite_non, &nitrite_band);
    try std.testing.expect(!try apply(&state, .{
        .target_layer = 0,
        .application_depth_m = 0,
        .layer_upper_depth_m = &.{ 0, 0.1 },
        .layer_thickness_m = &.{ 0.1, 0.2 },
        .requested_row_width_m = 0.05,
        .new_banded_ammonium_g_n_per_m2 = 0,
        .new_banded_ammonia_g_n_per_m2 = 0,
        .new_banded_urea_g_n_per_m2 = 0,
        .new_banded_nitrate_g_n_per_m2 = 0,
        .minimum_active_layer_thickness_m = 0.001,
    }));
    try std.testing.expectEqual(@as(f64, 1), nitrate_non[0]);
}

test "invalid target geometry fails before repartition" {
    var geometry = [_]f64{ 0, 0 };
    var band_fraction = [_]f64{ 0, 0 };
    var non_band_fraction = [_]f64{ 1, 1 };
    var nitrate_non = [_]f64{ 1, 1 };
    var nitrate_band = [_]f64{ 0, 0 };
    var nitrite_non = [_]f64{ 1, 1 };
    var nitrite_band = [_]f64{ 0, 0 };
    var state = makeState(&geometry, &band_fraction, &non_band_fraction, &nitrate_non, &nitrate_band, &nitrite_non, &nitrite_band);
    try std.testing.expectError(
        error.InvalidNitrateBandTargetDepth,
        apply(&state, .{
            .target_layer = 1,
            .application_depth_m = 0.05,
            .layer_upper_depth_m = &.{ 0, 0.1 },
            .layer_thickness_m = &.{ 0.1, 0.2 },
            .requested_row_width_m = 0.05,
            .new_banded_ammonium_g_n_per_m2 = 0,
            .new_banded_ammonia_g_n_per_m2 = 0,
            .new_banded_urea_g_n_per_m2 = 0,
            .new_banded_nitrate_g_n_per_m2 = 1,
            .minimum_active_layer_thickness_m = 0.001,
        }),
    );
    try std.testing.expectEqual(@as(f64, 1), nitrate_non[0]);
}
