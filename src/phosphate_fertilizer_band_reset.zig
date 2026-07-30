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
    /// Runtime pool list binds HOUR1 H0/H1/H2/H3 phosphate, Fe/Ca/Mg
    /// complexes, anion sites, exchangeable phosphate, and five precipitates.
    pools: []PairedPool,
};

pub const Inputs = struct {
    target_layer: usize,
    dissolved_h2po4_pool_index: usize,
    application_depth_m: f64,
    layer_upper_depth_m: []const f64,
    layer_thickness_m: []const f64,
    requested_row_width_m: f64,
    new_banded_monocalcium_phosphate_g_p_per_m2: f64,
    maximum_initial_band_dimension_m: f64 = 0.025,
    minimum_active_layer_thickness_m: f64,
    maximum_band_volume_fraction: f64 = 0.9999,
};

/// Exact PO4 band reset and conservative repartition from hour1.f:410-507.
pub fn apply(state: *State, inputs: Inputs) !bool {
    const layer_count = inputs.layer_thickness_m.len;
    if (layer_count == 0 or inputs.layer_upper_depth_m.len != layer_count or
        inputs.target_layer >= layer_count or state.pools.len == 0 or
        inputs.dissolved_h2po4_pool_index >= state.pools.len)
        return error.PhosphateBandDimensionMismatch;
    inline for (.{
        state.band_depth_m,
        state.band_width_m,
        state.band_volume_fraction,
        state.non_band_volume_fraction,
    }) |values| if (values.len != layer_count)
        return error.PhosphateBandDimensionMismatch;
    for (state.pools) |pool|
        if (pool.non_band.len != layer_count or pool.band.len != layer_count)
            return error.PhosphateBandDimensionMismatch;
    inline for (.{
        inputs.application_depth_m,
        inputs.requested_row_width_m,
        inputs.new_banded_monocalcium_phosphate_g_p_per_m2,
        inputs.maximum_initial_band_dimension_m,
        inputs.minimum_active_layer_thickness_m,
        inputs.maximum_band_volume_fraction,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFinitePhosphateBandInput;
    if (inputs.application_depth_m < 0 or inputs.requested_row_width_m <= 0 or
        inputs.new_banded_monocalcium_phosphate_g_p_per_m2 < 0 or
        inputs.maximum_initial_band_dimension_m <= 0 or
        inputs.minimum_active_layer_thickness_m < 0 or
        inputs.maximum_band_volume_fraction <= 0 or
        inputs.maximum_band_volume_fraction >= 1)
        return error.InvalidPhosphateBandInput;
    for (state.pools) |pool| for (pool.non_band, pool.band) |non_band, band| {
        if (!std.math.isFinite(non_band) or non_band < 0 or
            !std.math.isFinite(band) or band < 0 or
            !std.math.isFinite(non_band + band))
            return error.InvalidPhosphateBandState;
    };
    for (inputs.layer_upper_depth_m, inputs.layer_thickness_m) |upper, thickness|
        if (!std.math.isFinite(upper) or upper < 0 or
            !std.math.isFinite(thickness) or thickness <= 0)
            return error.InvalidPhosphateBandLayerGeometry;

    const h2po4 = state.pools[inputs.dissolved_h2po4_pool_index];
    const activate =
        inputs.new_banded_monocalcium_phosphate_g_p_per_m2 > 0 or
        (h2po4.band[inputs.target_layer] > 0 and !state.active);
    if (!activate) return false;

    const target_depth = @min(
        inputs.maximum_initial_band_dimension_m,
        inputs.application_depth_m -
            inputs.layer_upper_depth_m[inputs.target_layer],
    );
    if (!std.math.isFinite(target_depth) or target_depth < 0)
        return error.InvalidPhosphateBandTargetDepth;
    const target_width = @min(
        inputs.maximum_initial_band_dimension_m,
        inputs.requested_row_width_m,
    );
    for (0..layer_count) |layer| {
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
            return error.InvalidPhosphateBandFraction;
    }

    state.active = true;
    state.row_width_m = inputs.requested_row_width_m;
    state.band_center_depth_m = inputs.application_depth_m;
    state.band_upper_depth_m =
        inputs.application_depth_m -
        inputs.maximum_initial_band_dimension_m;
    for (0..layer_count) |layer| {
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
        for (state.pools) |pool| {
            const total = pool.non_band[layer] + pool.band[layer];
            pool.non_band[layer] = total * (1 - fraction);
            pool.band[layer] = total * fraction;
        }
    }
    return true;
}

test "runtime phosphate pools all repartition conservatively" {
    var depth = [_]f64{ 0, 0 };
    var width = [_]f64{ 0, 0 };
    var band_fraction = [_]f64{ 0, 0 };
    var non_band_fraction = [_]f64{ 1, 1 };
    var h2_non = [_]f64{ 10, 18 };
    var h2_band = [_]f64{ 0, 2 };
    var exchange_non = [_]f64{ 5, 9 };
    var exchange_band = [_]f64{ 0, 1 };
    var precip_non = [_]f64{ 20, 30 };
    var precip_band = [_]f64{ 0, 10 };
    var pools = [_]PairedPool{
        .{ .non_band = &h2_non, .band = &h2_band },
        .{ .non_band = &exchange_non, .band = &exchange_band },
        .{ .non_band = &precip_non, .band = &precip_band },
    };
    var state: State = .{
        .active = false,
        .row_width_m = 0,
        .band_center_depth_m = 0,
        .band_upper_depth_m = 0,
        .band_depth_m = &depth,
        .band_width_m = &width,
        .band_volume_fraction = &band_fraction,
        .non_band_volume_fraction = &non_band_fraction,
        .pools = &pools,
    };
    try std.testing.expect(try apply(&state, .{
        .target_layer = 1,
        .dissolved_h2po4_pool_index = 0,
        .application_depth_m = 0.12,
        .layer_upper_depth_m = &.{ 0, 0.1 },
        .layer_thickness_m = &.{ 0.1, 0.2 },
        .requested_row_width_m = 0.05,
        .new_banded_monocalcium_phosphate_g_p_per_m2 = 1,
        .minimum_active_layer_thickness_m = 0.001,
    }));
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), band_fraction[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), h2_band[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), exchange_band[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), precip_band[1], 1e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, 40),
        precip_non[1] + precip_band[1],
        1e-15,
    );
}

test "stranded H2PO4 band reactivates without new fertilizer" {
    var depth = [_]f64{0};
    var width = [_]f64{0};
    var band_fraction = [_]f64{0};
    var non_band_fraction = [_]f64{1};
    var h2_non = [_]f64{0};
    var h2_band = [_]f64{2};
    var pools = [_]PairedPool{
        .{ .non_band = &h2_non, .band = &h2_band },
    };
    var state: State = .{
        .active = false,
        .row_width_m = 0,
        .band_center_depth_m = 0,
        .band_upper_depth_m = 0,
        .band_depth_m = &depth,
        .band_width_m = &width,
        .band_volume_fraction = &band_fraction,
        .non_band_volume_fraction = &non_band_fraction,
        .pools = &pools,
    };
    try std.testing.expect(try apply(&state, .{
        .target_layer = 0,
        .dissolved_h2po4_pool_index = 0,
        .application_depth_m = 0.02,
        .layer_upper_depth_m = &.{0},
        .layer_thickness_m = &.{0.1},
        .requested_row_width_m = 0.05,
        .new_banded_monocalcium_phosphate_g_p_per_m2 = 0,
        .minimum_active_layer_thickness_m = 0.001,
    }));
    try std.testing.expect(state.active);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), h2_band[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.8), h2_non[0], 1e-15);
}

test "invalid late pool prevents every geometry and chemistry mutation" {
    var depth = [_]f64{7};
    var width = [_]f64{8};
    var band_fraction = [_]f64{0};
    var non_band_fraction = [_]f64{1};
    var h2_non = [_]f64{1};
    var h2_band = [_]f64{0};
    var bad_non = [_]f64{std.math.nan(f64)};
    var bad_band = [_]f64{0};
    var pools = [_]PairedPool{
        .{ .non_band = &h2_non, .band = &h2_band },
        .{ .non_band = &bad_non, .band = &bad_band },
    };
    var state: State = .{
        .active = false,
        .row_width_m = 0,
        .band_center_depth_m = 0,
        .band_upper_depth_m = 0,
        .band_depth_m = &depth,
        .band_width_m = &width,
        .band_volume_fraction = &band_fraction,
        .non_band_volume_fraction = &non_band_fraction,
        .pools = &pools,
    };
    try std.testing.expectError(
        error.InvalidPhosphateBandState,
        apply(&state, .{
            .target_layer = 0,
            .dissolved_h2po4_pool_index = 0,
            .application_depth_m = 0.02,
            .layer_upper_depth_m = &.{0},
            .layer_thickness_m = &.{0.1},
            .requested_row_width_m = 0.05,
            .new_banded_monocalcium_phosphate_g_p_per_m2 = 1,
            .minimum_active_layer_thickness_m = 0.001,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), depth[0]);
    try std.testing.expectEqual(@as(f64, 1), h2_non[0]);
}
