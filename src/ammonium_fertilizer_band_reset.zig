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
    soluble_ammonium_g_n: PairedPool,
    soluble_ammonia_g_n: PairedPool,
    exchangeable_ammonium_mol_n: PairedPool,
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
    maximum_initial_band_dimension_m: f64 = 0.025,
    minimum_active_layer_thickness_m: f64,
    maximum_band_volume_fraction: f64 = 0.9999,
};

/// Exact NH4 band reset and conservative repartition from hour1.f:303-342.
pub fn apply(state: *State, inputs: Inputs) !bool {
    const count = inputs.layer_thickness_m.len;
    if (count == 0 or inputs.layer_upper_depth_m.len != count or
        inputs.target_layer >= count)
        return error.AmmoniumBandDimensionMismatch;
    inline for (.{
        state.band_depth_m,
        state.band_width_m,
        state.band_volume_fraction,
        state.non_band_volume_fraction,
        state.soluble_ammonium_g_n.non_band,
        state.soluble_ammonium_g_n.band,
        state.soluble_ammonia_g_n.non_band,
        state.soluble_ammonia_g_n.band,
        state.exchangeable_ammonium_mol_n.non_band,
        state.exchangeable_ammonium_mol_n.band,
    }) |values| if (values.len != count)
        return error.AmmoniumBandDimensionMismatch;
    inline for (.{
        inputs.application_depth_m,
        inputs.requested_row_width_m,
        inputs.new_banded_ammonium_g_n_per_m2,
        inputs.new_banded_ammonia_g_n_per_m2,
        inputs.new_banded_urea_g_n_per_m2,
        inputs.maximum_initial_band_dimension_m,
        inputs.minimum_active_layer_thickness_m,
        inputs.maximum_band_volume_fraction,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteAmmoniumBandInput;
    if (inputs.application_depth_m < 0 or inputs.requested_row_width_m <= 0 or
        inputs.new_banded_ammonium_g_n_per_m2 < 0 or
        inputs.new_banded_ammonia_g_n_per_m2 < 0 or
        inputs.new_banded_urea_g_n_per_m2 < 0 or
        inputs.maximum_initial_band_dimension_m <= 0 or
        inputs.minimum_active_layer_thickness_m < 0 or
        inputs.maximum_band_volume_fraction <= 0 or
        inputs.maximum_band_volume_fraction >= 1)
        return error.InvalidAmmoniumBandInput;

    inline for (.{
        state.soluble_ammonium_g_n,
        state.soluble_ammonia_g_n,
        state.exchangeable_ammonium_mol_n,
    }) |pool| for (pool.non_band, pool.band) |non_band, band| {
        if (!std.math.isFinite(non_band) or non_band < 0 or
            !std.math.isFinite(band) or band < 0 or
            !std.math.isFinite(non_band + band))
            return error.InvalidAmmoniumBandState;
    };
    for (inputs.layer_upper_depth_m, inputs.layer_thickness_m) |upper, thickness|
        if (!std.math.isFinite(upper) or upper < 0 or
            !std.math.isFinite(thickness) or thickness <= 0)
            return error.InvalidAmmoniumBandLayerGeometry;

    const has_new_band =
        inputs.new_banded_ammonium_g_n_per_m2 +
        inputs.new_banded_ammonia_g_n_per_m2 +
        inputs.new_banded_urea_g_n_per_m2 > 0;
    const has_stranded_target_band =
        (state.soluble_ammonium_g_n.band[inputs.target_layer] > 0 or
            state.soluble_ammonia_g_n.band[inputs.target_layer] > 0) and
        !state.active;
    if (!has_new_band and !has_stranded_target_band) return false;

    const target_depth = @min(
        inputs.maximum_initial_band_dimension_m,
        inputs.application_depth_m -
            inputs.layer_upper_depth_m[inputs.target_layer],
    );
    if (!std.math.isFinite(target_depth) or target_depth < 0)
        return error.InvalidAmmoniumBandTargetDepth;
    const target_width = @min(
        inputs.maximum_initial_band_dimension_m,
        inputs.requested_row_width_m,
    );

    // Validate every derived fraction and conserved total before mutation.
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
            return error.InvalidAmmoniumBandFraction;
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
        inline for (.{
            &state.soluble_ammonium_g_n,
            &state.soluble_ammonia_g_n,
            &state.exchangeable_ammonium_mol_n,
        }) |pool| {
            const total = pool.non_band[layer] + pool.band[layer];
            pool.non_band[layer] = total * (1 - fraction);
            pool.band[layer] = total * fraction;
        }
    }
    return true;
}

fn stateFor(
    depth: []f64,
    width: []f64,
    band_fraction: []f64,
    non_band_fraction: []f64,
    nh4_non: []f64,
    nh4_band: []f64,
    nh3_non: []f64,
    nh3_band: []f64,
    exchange_non: []f64,
    exchange_band: []f64,
) State {
    return .{
        .active = false,
        .row_width_m = 0,
        .band_center_depth_m = 0,
        .band_upper_depth_m = 0,
        .band_depth_m = depth,
        .band_width_m = width,
        .band_volume_fraction = band_fraction,
        .non_band_volume_fraction = non_band_fraction,
        .soluble_ammonium_g_n = .{ .non_band = nh4_non, .band = nh4_band },
        .soluble_ammonia_g_n = .{ .non_band = nh3_non, .band = nh3_band },
        .exchangeable_ammonium_mol_n = .{
            .non_band = exchange_non,
            .band = exchange_band,
        },
    };
}

test "new band geometry conservatively repartitions all ammonium pools" {
    var depth = [_]f64{ 9, 9 };
    var width = depth;
    var band_fraction = depth;
    var non_band_fraction = depth;
    var nh4_non = [_]f64{ 8, 20 };
    var nh4_band = [_]f64{ 2, 0 };
    var nh3_non = [_]f64{ 4, 10 };
    var nh3_band = [_]f64{ 1, 0 };
    var exchange_non = [_]f64{ 6, 30 };
    var exchange_band = [_]f64{ 4, 0 };
    var state = stateFor(&depth, &width, &band_fraction, &non_band_fraction, &nh4_non, &nh4_band, &nh3_non, &nh3_band, &exchange_non, &exchange_band);
    try std.testing.expect(try apply(&state, .{
        .target_layer = 1,
        .application_depth_m = 0.12,
        .layer_upper_depth_m = &.{ 0, 0.1 },
        .layer_thickness_m = &.{ 0.1, 0.2 },
        .requested_row_width_m = 0.05,
        .new_banded_ammonium_g_n_per_m2 = 1,
        .new_banded_ammonia_g_n_per_m2 = 0,
        .new_banded_urea_g_n_per_m2 = 0,
        .minimum_active_layer_thickness_m = 0.001,
    }));
    try std.testing.expectEqual(@as(f64, 0), band_fraction[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), band_fraction[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 10), nh4_band[0] + nh4_non[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 20), nh4_band[1] + nh4_non[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), nh4_band[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), nh3_band[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), exchange_band[1], 1e-15);
}

test "no new or stranded band leaves state unchanged" {
    var values = [_]f64{0};
    var state = stateFor(&values, &values, &values, &values, &values, &values, &values, &values, &values, &values);
    try std.testing.expect(!try apply(&state, .{
        .target_layer = 0,
        .application_depth_m = 0,
        .layer_upper_depth_m = &.{0},
        .layer_thickness_m = &.{0.1},
        .requested_row_width_m = 0.05,
        .new_banded_ammonium_g_n_per_m2 = 0,
        .new_banded_ammonia_g_n_per_m2 = 0,
        .new_banded_urea_g_n_per_m2 = 0,
        .minimum_active_layer_thickness_m = 0.001,
    }));
}

test "stranded soluble band reactivates when flag is false" {
    var depth = [_]f64{0};
    var width = [_]f64{0};
    var fractions = [_]f64{0};
    var non_fractions = [_]f64{1};
    var nh4_non = [_]f64{0};
    var nh4_band = [_]f64{2};
    var nh3_non = [_]f64{0};
    var nh3_band = [_]f64{0};
    var exchange_non = [_]f64{0};
    var exchange_band = [_]f64{0};
    var state = stateFor(&depth, &width, &fractions, &non_fractions, &nh4_non, &nh4_band, &nh3_non, &nh3_band, &exchange_non, &exchange_band);
    try std.testing.expect(try apply(&state, .{
        .target_layer = 0,
        .application_depth_m = 0.02,
        .layer_upper_depth_m = &.{0},
        .layer_thickness_m = &.{0.1},
        .requested_row_width_m = 0.05,
        .new_banded_ammonium_g_n_per_m2 = 0,
        .new_banded_ammonia_g_n_per_m2 = 0,
        .new_banded_urea_g_n_per_m2 = 0,
        .minimum_active_layer_thickness_m = 0.001,
    }));
    try std.testing.expect(state.active);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), nh4_band[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.8), nh4_non[0], 1e-15);
}
