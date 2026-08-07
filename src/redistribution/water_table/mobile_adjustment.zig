const std = @import("std");

pub const Inputs = struct {
    water_table_type: []const i32, // IDTBL: 2 natural mobile, 4 artificial mobile
    soil_surface_elevation_m: []const f64, // CDPTH(NU-1)
    net_boundary_water_transfer_m3_step: []const f64, // HVOLO
    cell_area_m2: []const f64, // AREA(3,NU)
};
pub const State = struct {
    natural_reference_depth_m: []const f64, // DTBLZ
    natural_current_depth_m: []f64, // DTBLX
    artificial_reference_depth_m: []const f64, // DTBLD
    artificial_current_depth_m: []f64, // DTBLY
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 11090--11100 for one runtime-indexed cell.
pub fn adjustCell(cell: usize, inputs: Inputs, state: State) !bool {
    const cell_count = inputs.water_table_type.len;
    if (cell_count == 0 or cell >= cell_count) return error.MobileWaterTableDimensionMismatch;
    inline for (.{ inputs.soil_surface_elevation_m, inputs.net_boundary_water_transfer_m3_step, inputs.cell_area_m2, state.natural_reference_depth_m, state.natural_current_depth_m, state.artificial_reference_depth_m, state.artificial_current_depth_m }) |values|
        if (values.len != cell_count) return error.MobileWaterTableDimensionMismatch;
    inline for (.{ inputs.soil_surface_elevation_m, inputs.net_boundary_water_transfer_m3_step, inputs.cell_area_m2, state.natural_reference_depth_m, state.natural_current_depth_m, state.artificial_reference_depth_m, state.artificial_current_depth_m }) |values|
        if (!finiteSlice(values)) return error.InvalidMobileWaterTableInput;
    if (inputs.cell_area_m2[cell] <= 0) return error.InvalidMobileWaterTableInput;
    const water_table_type = inputs.water_table_type[cell];
    if (water_table_type != 2 and water_table_type != 4) return false;

    var next_natural = state.natural_reference_depth_m[cell] + inputs.soil_surface_elevation_m[cell];
    next_natural = next_natural - inputs.net_boundary_water_transfer_m3_step[cell] / inputs.cell_area_m2[cell] -
        0.00167 * (next_natural - state.natural_reference_depth_m[cell] - inputs.soil_surface_elevation_m[cell]);
    var next_artificial = state.artificial_current_depth_m[cell];
    if (water_table_type == 4) {
        next_artificial = next_artificial - inputs.net_boundary_water_transfer_m3_step[cell] / inputs.cell_area_m2[cell] -
            0.00167 * (next_artificial - state.artificial_reference_depth_m[cell]);
    }
    if (!std.math.isFinite(next_natural) or !std.math.isFinite(next_artificial))
        return error.NonFiniteMobileWaterTableResult;
    state.natural_current_depth_m[cell] = next_natural;
    if (water_table_type == 4) state.artificial_current_depth_m[cell] = next_artificial;
    return true;
}

test "REDIST mobile water table preserves IDTBL gates and sequential equations" {
    const types = [_]i32{ 0, 2, 4 };
    const surface = [_]f64{ 1, 2, 3 };
    const transfer = [_]f64{ 10, 20, 30 };
    const area = [_]f64{ 10, 10, 10 };
    const natural_reference = [_]f64{ 4, 5, 6 };
    var natural_current = [_]f64{ 7, 8, 9 };
    const artificial_reference = [_]f64{ 10, 11, 12 };
    var artificial_current = [_]f64{ 13, 14, 15 };
    const inputs = Inputs{ .water_table_type = &types, .soil_surface_elevation_m = &surface, .net_boundary_water_transfer_m3_step = &transfer, .cell_area_m2 = &area };
    const state = State{ .natural_reference_depth_m = &natural_reference, .natural_current_depth_m = &natural_current, .artificial_reference_depth_m = &artificial_reference, .artificial_current_depth_m = &artificial_current };
    try std.testing.expect(!try adjustCell(0, inputs, state));
    try std.testing.expect(try adjustCell(1, inputs, state));
    try std.testing.expectEqual(@as(f64, 5), natural_current[1]);
    try std.testing.expectEqual(@as(f64, 14), artificial_current[1]);
    try std.testing.expect(try adjustCell(2, inputs, state));
    try std.testing.expectEqual(@as(f64, 6), natural_current[2]);
    try std.testing.expectApproxEqAbs(15 - 3 - 0.00167 * (15 - 12), artificial_current[2], 1e-12);
}

test "REDIST mobile artificial late overflow leaves natural update atomic" {
    const types = [_]i32{4};
    const surface = [_]f64{1};
    const transfer = [_]f64{-std.math.floatMax(f64)};
    const area = [_]f64{1};
    const natural_reference = [_]f64{2};
    var natural_current = [_]f64{3};
    const artificial_reference = [_]f64{0};
    var artificial_current = [_]f64{std.math.floatMax(f64)};
    try std.testing.expectError(error.NonFiniteMobileWaterTableResult, adjustCell(0, .{ .water_table_type = &types, .soil_surface_elevation_m = &surface, .net_boundary_water_transfer_m3_step = &transfer, .cell_area_m2 = &area }, .{ .natural_reference_depth_m = &natural_reference, .natural_current_depth_m = &natural_current, .artificial_reference_depth_m = &artificial_reference, .artificial_current_depth_m = &artificial_current }));
    try std.testing.expectEqual(@as(f64, 3), natural_current[0]);
    try std.testing.expectEqual(std.math.floatMax(f64), artificial_current[0]);
}
