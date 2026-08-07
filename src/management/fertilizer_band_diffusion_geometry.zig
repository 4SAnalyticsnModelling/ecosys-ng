const std = @import("std");

pub const LayerProperties = struct {
    top_depth_m: f64,
    bottom_depth_m: f64,
    thickness_m: f64,
    ammonium_diffusivity_m2_h: f64,
    tortuosity: f64,
};

pub const AmmoniumBandGeometry = struct {
    total_depth_m: f64,
    penetration_front_depth_m: f64,
    layer_band_depth_m: []f64,
    layer_band_width_m: []f64,
    nonband_volume_fraction: []f64,
    band_volume_fraction: []f64,
    nonband_fractional_change_per_timestep: []f64,
};

pub const AmmoniumPools = struct {
    aqueous_ammonium_nonband_g_n: []f64,
    aqueous_ammonium_band_g_n: []f64,
    aqueous_ammonia_nonband_g_n: []f64,
    aqueous_ammonia_band_g_n: []f64,
    adsorbed_ammonium_nonband_g_n: []f64,
    adsorbed_ammonium_band_g_n: []f64,
    fertilizer_ammonium_nonband_g_n: []f64,
    fertilizer_ammonium_band_g_n: []f64,
    fertilizer_ammonia_nonband_g_n: []f64,
    fertilizer_ammonia_band_g_n: []f64,
    fertilizer_urea_nonband_g_n: []f64,
    fertilizer_urea_band_g_n: []f64,
};

pub const BandApplication = enum {
    unbanded,
    banded,
};

pub const Inputs = struct {
    band_application: BandApplication,
    first_active_layer_index: usize,
    layer_index: usize,
    row_width_m: f64,
    solute_timestep_h: f64,
    depth_threshold_m: f64,
    minimum_layer_thickness_m: f64,
    absent_band_fraction_threshold: f64,
    maximum_band_volume_fraction: f64,
    layers: []const LayerProperties,
};

pub const UpdateError = error{
    LayerCountMismatch,
    LayerIndexOutOfBounds,
    NonFiniteInput,
    InvalidGeometry,
    InvalidFraction,
    NegativeTransportProperty,
    NonPositivePreviousNonbandFraction,
    NonFiniteResult,
};

/// Translates the NH4 branch of `hour1.f` lines 4888--4983 for one soil layer.
///
/// Call this once per active layer, in increasing layer order, interleaved with
/// the corresponding nitrate and phosphate layer updates to retain HOUR1 order.
pub fn updateAmmoniumLayer(
    inputs: Inputs,
    geometry: *AmmoniumBandGeometry,
    pools: *AmmoniumPools,
) UpdateError!void {
    try validate(inputs, geometry.*, pools.*);
    if (inputs.band_application != .banded or inputs.row_width_m <= 0.0) return;

    const layer_index = inputs.layer_index;
    const layer = inputs.layers[layer_index];
    const layer_is_within_band =
        layer_index == inputs.first_active_layer_index or layer.top_depth_m < geometry.total_depth_m;
    if (!layer_is_within_band) {
        amalgamateLayer(layer_index, geometry, pools);
        return;
    }

    var width_change_m: f64 = 0.0;
    if (geometry.layer_band_depth_m[layer_index] > inputs.depth_threshold_m) {
        width_change_m = 0.5 *
            @sqrt(layer.ammonium_diffusivity_m2_h * layer.tortuosity) *
            inputs.solute_timestep_h;
        geometry.layer_band_width_m[layer_index] = @min(
            inputs.row_width_m,
            geometry.layer_band_width_m[layer_index] + width_change_m,
        );
    } else {
        geometry.layer_band_width_m[layer_index] = 0.0;
    }

    if (layer.bottom_depth_m >= geometry.penetration_front_depth_m and
        layer.top_depth_m < geometry.penetration_front_depth_m)
    {
        geometry.penetration_front_depth_m =
            @max(0.0, geometry.penetration_front_depth_m - width_change_m);
        geometry.layer_band_depth_m[layer_index] = @min(
            layer.thickness_m,
            geometry.layer_band_depth_m[layer_index] + width_change_m,
        );
    }
    if (layer_index > inputs.first_active_layer_index and
        geometry.penetration_front_depth_m < layer.top_depth_m and
        geometry.band_volume_fraction[layer_index - 1] <
            inputs.absent_band_fraction_threshold)
    {
        geometry.layer_band_depth_m[layer_index - 1] =
            layer.top_depth_m - geometry.penetration_front_depth_m;
        geometry.layer_band_width_m[layer_index - 1] =
            geometry.layer_band_width_m[layer_index];
    }
    if (layer.bottom_depth_m >= geometry.total_depth_m and
        layer.top_depth_m < geometry.total_depth_m)
    {
        geometry.total_depth_m += width_change_m;
        geometry.layer_band_depth_m[layer_index] = @min(
            layer.thickness_m,
            geometry.layer_band_depth_m[layer_index] + width_change_m,
        );
    }
    if (geometry.total_depth_m > layer.bottom_depth_m and
        geometry.penetration_front_depth_m <= layer.top_depth_m)
    {
        const next_layer_index = layer_index + 1;
        if (next_layer_index >= inputs.layers.len) return error.LayerIndexOutOfBounds;
        if (geometry.band_volume_fraction[next_layer_index] <
            inputs.absent_band_fraction_threshold)
        {
            geometry.layer_band_width_m[next_layer_index] =
                geometry.layer_band_width_m[layer_index];
            geometry.layer_band_depth_m[next_layer_index] =
                geometry.total_depth_m - layer.bottom_depth_m;
        }
    }

    const previous_nonband_fraction = geometry.nonband_volume_fraction[layer_index];
    if (previous_nonband_fraction <= 0.0) {
        return error.NonPositivePreviousNonbandFraction;
    }
    if (layer.thickness_m > inputs.minimum_layer_thickness_m) {
        geometry.band_volume_fraction[layer_index] = @max(
            0.0,
            @min(
                inputs.maximum_band_volume_fraction,
                geometry.layer_band_width_m[layer_index] / inputs.row_width_m *
                    geometry.layer_band_depth_m[layer_index] / layer.thickness_m,
            ),
        );
    } else {
        geometry.band_volume_fraction[layer_index] = 0.0;
    }
    geometry.nonband_volume_fraction[layer_index] =
        1.0 - geometry.band_volume_fraction[layer_index];
    geometry.nonband_fractional_change_per_timestep[layer_index] = @min(
        0.0,
        (geometry.nonband_volume_fraction[layer_index] - previous_nonband_fraction) /
            previous_nonband_fraction,
    );

    if (!std.math.isFinite(width_change_m) or
        !std.math.isFinite(geometry.total_depth_m) or
        !std.math.isFinite(geometry.penetration_front_depth_m) or
        !std.math.isFinite(geometry.nonband_fractional_change_per_timestep[layer_index]))
    {
        return error.NonFiniteResult;
    }
}

fn amalgamateLayer(
    layer_index: usize,
    geometry: *AmmoniumBandGeometry,
    pools: *AmmoniumPools,
) void {
    geometry.nonband_fractional_change_per_timestep[layer_index] = 0.0;
    geometry.layer_band_depth_m[layer_index] = 0.0;
    geometry.layer_band_width_m[layer_index] = 0.0;
    geometry.nonband_volume_fraction[layer_index] = 1.0;
    geometry.band_volume_fraction[layer_index] = 0.0;

    inline for (.{
        .{ "aqueous_ammonium_nonband_g_n", "aqueous_ammonium_band_g_n" },
        .{ "aqueous_ammonia_nonband_g_n", "aqueous_ammonia_band_g_n" },
        .{ "adsorbed_ammonium_nonband_g_n", "adsorbed_ammonium_band_g_n" },
        .{ "fertilizer_ammonium_nonband_g_n", "fertilizer_ammonium_band_g_n" },
        .{ "fertilizer_ammonia_nonband_g_n", "fertilizer_ammonia_band_g_n" },
        .{ "fertilizer_urea_nonband_g_n", "fertilizer_urea_band_g_n" },
    }) |field_names| {
        @field(pools, field_names[0])[layer_index] +=
            @field(pools, field_names[1])[layer_index];
        @field(pools, field_names[1])[layer_index] = 0.0;
    }
}

fn validate(
    inputs: Inputs,
    geometry: AmmoniumBandGeometry,
    pools: AmmoniumPools,
) UpdateError!void {
    const layer_count = inputs.layers.len;
    inline for (std.meta.fields(AmmoniumBandGeometry)) |field| {
        if (field.type != f64 and @field(geometry, field.name).len != layer_count) {
            return error.LayerCountMismatch;
        }
    }
    inline for (std.meta.fields(AmmoniumPools)) |field| {
        if (@field(pools, field.name).len != layer_count) return error.LayerCountMismatch;
    }
    if (inputs.layer_index >= layer_count or inputs.first_active_layer_index > inputs.layer_index) {
        return error.LayerIndexOutOfBounds;
    }
    const scalar_values = [_]f64{
        inputs.row_width_m,
        inputs.solute_timestep_h,
        inputs.depth_threshold_m,
        inputs.minimum_layer_thickness_m,
        inputs.absent_band_fraction_threshold,
        inputs.maximum_band_volume_fraction,
        geometry.total_depth_m,
        geometry.penetration_front_depth_m,
    };
    for (scalar_values) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
    }
    if (inputs.row_width_m < 0.0 or inputs.solute_timestep_h < 0.0 or
        inputs.depth_threshold_m < 0.0 or inputs.minimum_layer_thickness_m < 0.0)
    {
        return error.InvalidGeometry;
    }
    if (inputs.absent_band_fraction_threshold < 0.0 or
        inputs.maximum_band_volume_fraction <= 0.0 or
        inputs.maximum_band_volume_fraction >= 1.0)
    {
        return error.InvalidFraction;
    }
    for (inputs.layers) |layer| {
        inline for (std.meta.fields(LayerProperties)) |field| {
            if (!std.math.isFinite(@field(layer, field.name))) return error.NonFiniteInput;
        }
        if (layer.top_depth_m < 0.0 or layer.bottom_depth_m < layer.top_depth_m or
            layer.thickness_m < 0.0)
        {
            return error.InvalidGeometry;
        }
        if (layer.ammonium_diffusivity_m2_h < 0.0 or layer.tortuosity < 0.0) {
            return error.NegativeTransportProperty;
        }
    }
    inline for (std.meta.fields(AmmoniumBandGeometry)) |field| {
        if (field.type != f64) {
            for (@field(geometry, field.name)) |value| {
                if (!std.math.isFinite(value)) return error.NonFiniteInput;
            }
        }
    }
    inline for (std.meta.fields(AmmoniumPools)) |field| {
        for (@field(pools, field.name)) |value| {
            if (!std.math.isFinite(value)) return error.NonFiniteInput;
        }
    }
}

fn zeroPools(storage: *[12][2]f64) AmmoniumPools {
    return .{
        .aqueous_ammonium_nonband_g_n = &storage[0],
        .aqueous_ammonium_band_g_n = &storage[1],
        .aqueous_ammonia_nonband_g_n = &storage[2],
        .aqueous_ammonia_band_g_n = &storage[3],
        .adsorbed_ammonium_nonband_g_n = &storage[4],
        .adsorbed_ammonium_band_g_n = &storage[5],
        .fertilizer_ammonium_nonband_g_n = &storage[6],
        .fertilizer_ammonium_band_g_n = &storage[7],
        .fertilizer_ammonia_nonband_g_n = &storage[8],
        .fertilizer_ammonia_band_g_n = &storage[9],
        .fertilizer_urea_nonband_g_n = &storage[10],
        .fertilizer_urea_band_g_n = &storage[11],
    };
}

test "ammonium band expands width depth and volume fraction" {
    const layers = [_]LayerProperties{
        .{
            .top_depth_m = 0.0,
            .bottom_depth_m = 0.1,
            .thickness_m = 0.1,
            .ammonium_diffusivity_m2_h = 4.0e-6,
            .tortuosity = 0.25,
        },
        .{
            .top_depth_m = 0.1,
            .bottom_depth_m = 0.2,
            .thickness_m = 0.1,
            .ammonium_diffusivity_m2_h = 4.0e-6,
            .tortuosity = 0.25,
        },
    };
    var depths = [_]f64{ 0.02, 0.0 };
    var widths = [_]f64{ 0.01, 0.0 };
    var nonband = [_]f64{ 0.99, 1.0 };
    var band = [_]f64{ 0.01, 0.0 };
    var changes = [_]f64{ 0.0, 0.0 };
    var geometry = AmmoniumBandGeometry{
        .total_depth_m = 0.05,
        .penetration_front_depth_m = 0.05,
        .layer_band_depth_m = &depths,
        .layer_band_width_m = &widths,
        .nonband_volume_fraction = &nonband,
        .band_volume_fraction = &band,
        .nonband_fractional_change_per_timestep = &changes,
    };
    var pool_storage = std.mem.zeroes([12][2]f64);
    var pools = zeroPools(&pool_storage);

    try updateAmmoniumLayer(.{
        .band_application = .banded,
        .first_active_layer_index = 0,
        .layer_index = 0,
        .row_width_m = 0.2,
        .solute_timestep_h = 1.0,
        .depth_threshold_m = 0.0,
        .minimum_layer_thickness_m = 1.0e-9,
        .absent_band_fraction_threshold = 1.0e-12,
        .maximum_band_volume_fraction = 0.9999,
        .layers = &layers,
    }, &geometry, &pools);

    const expected_change_m = 0.5 * @sqrt(4.0e-6 * 0.25);
    try std.testing.expectApproxEqRel(
        0.01 + expected_change_m,
        widths[0],
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        0.02 + 2.0 * expected_change_m,
        depths[0],
        1.0e-14,
    );
    try std.testing.expect(band[0] > 0.0 and band[0] < 1.0);
    try std.testing.expectApproxEqAbs(1.0, band[0] + nonband[0], 1.0e-15);
}

test "layer beyond ammonium band amalgamates band pools" {
    const layers = [_]LayerProperties{
        .{
            .top_depth_m = 0.0,
            .bottom_depth_m = 0.1,
            .thickness_m = 0.1,
            .ammonium_diffusivity_m2_h = 0.0,
            .tortuosity = 0.0,
        },
        .{
            .top_depth_m = 0.1,
            .bottom_depth_m = 0.2,
            .thickness_m = 0.1,
            .ammonium_diffusivity_m2_h = 0.0,
            .tortuosity = 0.0,
        },
    };
    var depths = [_]f64{ 0.0, 0.03 };
    var widths = [_]f64{ 0.0, 0.04 };
    var nonband = [_]f64{ 1.0, 0.8 };
    var band = [_]f64{ 0.0, 0.2 };
    var changes = [_]f64{ 0.0, -0.1 };
    var geometry = AmmoniumBandGeometry{
        .total_depth_m = 0.05,
        .penetration_front_depth_m = 0.05,
        .layer_band_depth_m = &depths,
        .layer_band_width_m = &widths,
        .nonband_volume_fraction = &nonband,
        .band_volume_fraction = &band,
        .nonband_fractional_change_per_timestep = &changes,
    };
    var pool_storage = std.mem.zeroes([12][2]f64);
    pool_storage[0][1] = 2.0;
    pool_storage[1][1] = 3.0;
    pool_storage[10][1] = 5.0;
    pool_storage[11][1] = 7.0;
    var pools = zeroPools(&pool_storage);

    try updateAmmoniumLayer(.{
        .band_application = .banded,
        .first_active_layer_index = 0,
        .layer_index = 1,
        .row_width_m = 0.2,
        .solute_timestep_h = 1.0,
        .depth_threshold_m = 0.0,
        .minimum_layer_thickness_m = 1.0e-9,
        .absent_band_fraction_threshold = 1.0e-12,
        .maximum_band_volume_fraction = 0.9999,
        .layers = &layers,
    }, &geometry, &pools);

    try std.testing.expectEqual(@as(f64, 5.0), pool_storage[0][1]);
    try std.testing.expectEqual(@as(f64, 0.0), pool_storage[1][1]);
    try std.testing.expectEqual(@as(f64, 12.0), pool_storage[10][1]);
    try std.testing.expectEqual(@as(f64, 0.0), pool_storage[11][1]);
    try std.testing.expectEqual(@as(f64, 1.0), nonband[1]);
    try std.testing.expectEqual(@as(f64, 0.0), band[1]);
}
