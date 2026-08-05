const std = @import("std");

pub const Context = struct {
    source_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    minimum_layer_thickness_m: f64,
    source_layer_thickness_m: f64,
    destination_layer_thickness_m: f64,
    redistribution_fraction: f64,
    band_enabled_flag: u8,
    row_spacing_m: f64,
    current_is_upper_layer: bool,
    preceding_cumulative_depth_m: f64,
    application_depth_m: f64,
    current_layer_bottom_depth_m: f64,
    current_layer_thickness_m: f64,
    current_band_width_m: f64,
    current_band_depth_m: f64,
};

pub const State = struct {
    source_band_width_m: f64,
    destination_band_width_m: f64,
    source_band_depth_m: f64,
    destination_band_depth_m: f64,
    source_banded_fraction: f64,
    destination_banded_fraction: f64,
    source_nonbanded_fraction: f64,
    destination_nonbanded_fraction: f64,
};

fn finiteStruct(value: anytype) bool {
    inline for (std.meta.fields(@TypeOf(value))) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(value, field.name))) return false;
    }
    return true;
}

/// Direct translation of the PO4 band-dimension block in REDIST 9559--9579.
pub fn redistribute(context: Context, state: *State) !void {
    if (!finiteStruct(context) or !finiteStruct(state.*) or context.redistribution_fraction < 0 or
        context.redistribution_fraction > 1 or context.minimum_layer_thickness_m < 0)
        return error.InvalidSoilPhosphateBandGeometryInput;
    if (context.source_layer == 0 or context.source_bulk_density_megagrams_m3 <= 0 or
        context.destination_bulk_density_megagrams_m3 <= 0 or
        context.destination_layer_thickness_m <= context.minimum_layer_thickness_m or
        context.source_layer_thickness_m <= context.minimum_layer_thickness_m or
        context.band_enabled_flag != 1 or context.row_spacing_m <= 0 or
        (!context.current_is_upper_layer and context.preceding_cumulative_depth_m >= context.application_depth_m)) return;

    var next = state.*;
    const current_band_area_m2 = context.current_band_width_m * context.current_layer_thickness_m;
    var source_band_area_m2 = state.source_band_width_m * context.source_layer_thickness_m;
    var destination_band_area_m2 = state.destination_band_width_m * context.destination_layer_thickness_m;
    const transferred_band_area_m2 = @min(context.redistribution_fraction * current_band_area_m2, source_band_area_m2);
    destination_band_area_m2 = destination_band_area_m2 + transferred_band_area_m2;
    source_band_area_m2 = source_band_area_m2 - transferred_band_area_m2;
    next.destination_band_width_m = destination_band_area_m2 / context.destination_layer_thickness_m;
    next.source_band_width_m = source_band_area_m2 / context.source_layer_thickness_m;

    if (context.current_layer_bottom_depth_m >= context.application_depth_m) {
        const transferred_band_depth_m = @min(context.redistribution_fraction * context.current_band_depth_m, state.source_band_depth_m);
        next.destination_band_depth_m = state.destination_band_depth_m + transferred_band_depth_m;
        next.source_band_depth_m = state.source_band_depth_m - transferred_band_depth_m;
    }
    next.destination_banded_fraction = @max(0, @min(0.9999, next.destination_band_width_m / context.row_spacing_m * next.destination_band_depth_m / context.destination_layer_thickness_m));
    next.source_banded_fraction = @max(0, @min(0.9999, next.source_band_width_m / context.row_spacing_m * next.source_band_depth_m / context.source_layer_thickness_m));
    next.destination_nonbanded_fraction = 1.0 - next.destination_banded_fraction;
    next.source_nonbanded_fraction = 1.0 - next.source_banded_fraction;
    if (!finiteStruct(next)) return error.NonFiniteSoilPhosphateBandGeometryResult;
    state.* = next;
}

test "REDIST phosphate band geometry preserves transfer and update order" {
    var state = State{ .source_band_width_m = 0.4, .destination_band_width_m = 0.1, .source_band_depth_m = 0.3, .destination_band_depth_m = 0.1, .source_banded_fraction = 0, .destination_banded_fraction = 0, .source_nonbanded_fraction = 1, .destination_nonbanded_fraction = 1 };
    try redistribute(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .minimum_layer_thickness_m = 0.01, .source_layer_thickness_m = 0.5, .destination_layer_thickness_m = 0.5, .redistribution_fraction = 0.5, .band_enabled_flag = 1, .row_spacing_m = 1, .current_is_upper_layer = true, .preceding_cumulative_depth_m = 0, .application_depth_m = 0.2, .current_layer_bottom_depth_m = 0.3, .current_layer_thickness_m = 0.4, .current_band_width_m = 0.5, .current_band_depth_m = 0.2 }, &state);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.source_band_width_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), state.destination_band_width_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.source_band_depth_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.destination_band_depth_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), state.destination_banded_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.08), state.source_banded_fraction, 1e-15);
}

test "REDIST phosphate band geometry obeys source and thickness gates" {
    var state = State{ .source_band_width_m = 1, .destination_band_width_m = 1, .source_band_depth_m = 1, .destination_band_depth_m = 1, .source_banded_fraction = 0, .destination_banded_fraction = 0, .source_nonbanded_fraction = 1, .destination_nonbanded_fraction = 1 };
    const before = state;
    try redistribute(.{ .source_layer = 0, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .minimum_layer_thickness_m = 0.1, .source_layer_thickness_m = 1, .destination_layer_thickness_m = 1, .redistribution_fraction = 1, .band_enabled_flag = 1, .row_spacing_m = 1, .current_is_upper_layer = true, .preceding_cumulative_depth_m = 0, .application_depth_m = 1, .current_layer_bottom_depth_m = 1, .current_layer_thickness_m = 1, .current_band_width_m = 1, .current_band_depth_m = 1 }, &state);
    try std.testing.expectEqualDeep(before, state);
}

test "REDIST phosphate band geometry failure is atomic" {
    var state = State{ .source_band_width_m = 1, .destination_band_width_m = 1, .source_band_depth_m = 1, .destination_band_depth_m = 1, .source_banded_fraction = 0, .destination_banded_fraction = 0, .source_nonbanded_fraction = 1, .destination_nonbanded_fraction = 1 };
    const before = state;
    try std.testing.expectError(error.InvalidSoilPhosphateBandGeometryInput, redistribute(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .minimum_layer_thickness_m = 0.1, .source_layer_thickness_m = 1, .destination_layer_thickness_m = 1, .redistribution_fraction = std.math.nan(f64), .band_enabled_flag = 1, .row_spacing_m = 1, .current_is_upper_layer = true, .preceding_cumulative_depth_m = 0, .application_depth_m = 1, .current_layer_bottom_depth_m = 1, .current_layer_thickness_m = 1, .current_band_width_m = 1, .current_band_depth_m = 1 }, &state));
    try std.testing.expectEqualDeep(before, state);
}
