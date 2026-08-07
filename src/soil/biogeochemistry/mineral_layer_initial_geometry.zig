const std = @import("std");

pub const Inputs = struct {
    cell_count: usize,
    layer_capacity: usize,
    first_active_layer: []const usize,
    active_layer_count: []const usize,
    /// Explicit `[cell][boundary]`, with `layer_capacity + 1` boundaries.
    boundary_depth_m: []const f64,
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
    matrix_volume_fraction: []const f64,
    bulk_density_megagrams_per_m3: []const f64,
    initial_bulk_density_megagrams_per_m3: []const f64,
    calculation_floor: f64,
};

pub const State = struct {
    macropore_fraction: []f64,
    layer_thickness_m: []f64,
    layer_midpoint_depth_m: []f64,
    layer_bottom_depth_from_surface_m: []f64,
    layer_midpoint_depth_from_surface_m: []f64,
    total_volume_m3: []f64,
    matrix_volume_m3: []f64,
    reference_matrix_volume_m3: []f64,
    initial_total_volume_m3: []f64,
    dry_soil_mass_megagrams: []f64,
    total_root_density_m_per_m3: []f64,
    east_west_face_area_m2: []f64,
    north_south_face_area_m2: []f64,
    surface_boundary_depth_m: []f64,
    initial_surface_boundary_depth_m: []f64,
};

/// Exact mineral-layer branch of legacy `STARTS` lines 568--586.
pub fn initialize(state: State, inputs: Inputs) !void {
    if (inputs.cell_count == 0 or inputs.layer_capacity == 0)
        return error.InvalidMineralLayerGeometryDimensions;
    const layer_count = std.math.mul(
        usize,
        inputs.cell_count,
        inputs.layer_capacity,
    ) catch return error.DimensionOverflow;
    const boundary_stride = std.math.add(
        usize,
        inputs.layer_capacity,
        1,
    ) catch return error.DimensionOverflow;
    const boundary_count = std.math.mul(
        usize,
        inputs.cell_count,
        boundary_stride,
    ) catch return error.DimensionOverflow;
    if (inputs.first_active_layer.len != inputs.cell_count or
        inputs.active_layer_count.len != inputs.cell_count or
        inputs.boundary_depth_m.len != boundary_count or
        inputs.horizontal_cell_width_m.len != inputs.cell_count or
        inputs.vertical_cell_width_m.len != inputs.cell_count or
        inputs.matrix_volume_fraction.len != layer_count or
        inputs.bulk_density_megagrams_per_m3.len != layer_count or
        inputs.initial_bulk_density_megagrams_per_m3.len != layer_count)
        return error.MineralLayerGeometryDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const expected = if (std.mem.eql(
            u8,
            field.name,
            "surface_boundary_depth_m",
        ) or std.mem.eql(
            u8,
            field.name,
            "initial_surface_boundary_depth_m",
        ))
            inputs.cell_count
        else
            layer_count;
        if (@field(state, field.name).len != expected)
            return error.MineralLayerGeometryDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.calculation_floor) or
        inputs.calculation_floor <= 0)
        return error.InvalidMineralLayerGeometryInput;

    for (0..inputs.cell_count) |cell| {
        const first = inputs.first_active_layer[cell];
        const count = inputs.active_layer_count[cell];
        if (count == 0 or first >= inputs.layer_capacity or
            count > inputs.layer_capacity - first)
            return error.InvalidActiveMineralLayerRange;
        inline for (.{
            inputs.horizontal_cell_width_m[cell],
            inputs.vertical_cell_width_m[cell],
        }) |width_m| {
            if (!std.math.isFinite(width_m))
                return error.NonFiniteMineralLayerGeometryInput;
            if (width_m <= 0) return error.InvalidMineralLayerGeometryInput;
        }
        const boundary_base = cell * boundary_stride;
        var previous_relative_bottom_m: f64 = 0.0;
        for (first..first + count) |layer| {
            const index = cell * inputs.layer_capacity + layer;
            const top_m = inputs.boundary_depth_m[boundary_base + layer];
            const bottom_m =
                inputs.boundary_depth_m[boundary_base + layer + 1];
            inline for (.{
                top_m,
                bottom_m,
                inputs.matrix_volume_fraction[index],
                inputs.bulk_density_megagrams_per_m3[index],
                inputs.initial_bulk_density_megagrams_per_m3[index],
                state.macropore_fraction[index],
            }) |value| {
                if (!std.math.isFinite(value))
                    return error.NonFiniteMineralLayerGeometryInput;
            }
            if (bottom_m <= top_m or
                inputs.matrix_volume_fraction[index] < 0 or
                inputs.matrix_volume_fraction[index] > 1 or
                inputs.bulk_density_megagrams_per_m3[index] < 0 or
                inputs.initial_bulk_density_megagrams_per_m3[index] < 0 or
                state.macropore_fraction[index] < 0 or
                state.macropore_fraction[index] > 1)
                return error.InvalidMineralLayerGeometryInput;

            const thickness_m = bottom_m - top_m;
            const area_m2 = inputs.horizontal_cell_width_m[cell] *
                inputs.vertical_cell_width_m[cell];
            const total_volume_m3 = area_m2 * thickness_m;
            const matrix_volume_m3 = total_volume_m3 *
                inputs.matrix_volume_fraction[index];
            const first_bottom_m =
                inputs.boundary_depth_m[boundary_base + first + 1];
            const first_thickness_m =
                first_bottom_m -
                inputs.boundary_depth_m[boundary_base + first];
            const relative_bottom_m =
                bottom_m - first_bottom_m + first_thickness_m;
            const relative_midpoint_m =
                0.5 * (relative_bottom_m + previous_relative_bottom_m);
            inline for (.{
                thickness_m,
                0.5 * (bottom_m + top_m),
                relative_bottom_m,
                relative_midpoint_m,
                total_volume_m3,
                matrix_volume_m3,
                inputs.bulk_density_megagrams_per_m3[index] * matrix_volume_m3,
                thickness_m * inputs.vertical_cell_width_m[cell],
                thickness_m * inputs.horizontal_cell_width_m[cell],
            }) |candidate| {
                if (!std.math.isFinite(candidate))
                    return error.MineralLayerGeometryOverflow;
            }
            previous_relative_bottom_m = relative_bottom_m;
        }
    }

    for (0..inputs.cell_count) |cell| {
        const first = inputs.first_active_layer[cell];
        const count = inputs.active_layer_count[cell];
        const boundary_base = cell * boundary_stride;
        var previous_relative_bottom_m: f64 = 0.0;
        for (first..first + count) |layer| {
            const index = cell * inputs.layer_capacity + layer;
            if (inputs.initial_bulk_density_megagrams_per_m3[index] <=
                inputs.calculation_floor)
                state.macropore_fraction[index] = 0.0;
            const top_m = inputs.boundary_depth_m[boundary_base + layer];
            const bottom_m =
                inputs.boundary_depth_m[boundary_base + layer + 1];
            state.layer_thickness_m[index] = bottom_m - top_m;
            state.layer_midpoint_depth_m[index] =
                0.5 * (bottom_m + top_m);
            const first_bottom_m =
                inputs.boundary_depth_m[boundary_base + first + 1];
            const first_thickness_m =
                first_bottom_m -
                inputs.boundary_depth_m[boundary_base + first];
            state.layer_bottom_depth_from_surface_m[index] =
                bottom_m - first_bottom_m + first_thickness_m;
            state.layer_midpoint_depth_from_surface_m[index] =
                0.5 *
                (state.layer_bottom_depth_from_surface_m[index] +
                    previous_relative_bottom_m);
            previous_relative_bottom_m =
                state.layer_bottom_depth_from_surface_m[index];
            state.total_volume_m3[index] =
                inputs.horizontal_cell_width_m[cell] *
                inputs.vertical_cell_width_m[cell] *
                state.layer_thickness_m[index];
            state.matrix_volume_m3[index] =
                state.total_volume_m3[index] *
                inputs.matrix_volume_fraction[index];
            state.reference_matrix_volume_m3[index] =
                state.matrix_volume_m3[index];
            state.initial_total_volume_m3[index] =
                state.total_volume_m3[index];
            state.dry_soil_mass_megagrams[index] =
                inputs.bulk_density_megagrams_per_m3[index] *
                state.matrix_volume_m3[index];
            state.total_root_density_m_per_m3[index] = 0.0;
            state.east_west_face_area_m2[index] =
                state.layer_thickness_m[index] *
                inputs.vertical_cell_width_m[cell];
            state.north_south_face_area_m2[index] =
                state.layer_thickness_m[index] *
                inputs.horizontal_cell_width_m[cell];
        }
        state.surface_boundary_depth_m[cell] =
            inputs.boundary_depth_m[boundary_base + first + 1] -
            state.layer_thickness_m[cell * inputs.layer_capacity + first];
        state.initial_surface_boundary_depth_m[cell] =
            state.surface_boundary_depth_m[cell];
    }
}

test "STARTS mineral geometry reproduces volume mass depth and face equations" {
    const cells = 1;
    const layers = 2;
    var layer_fields = [_]f64{0.25} ** (13 * cells * layers);
    var cell_fields = [_]f64{9.0} ** (2 * cells);
    const state: State = .{
        .macropore_fraction = layer_fields[0..2],
        .layer_thickness_m = layer_fields[2..4],
        .layer_midpoint_depth_m = layer_fields[4..6],
        .layer_bottom_depth_from_surface_m = layer_fields[6..8],
        .layer_midpoint_depth_from_surface_m = layer_fields[8..10],
        .total_volume_m3 = layer_fields[10..12],
        .matrix_volume_m3 = layer_fields[12..14],
        .reference_matrix_volume_m3 = layer_fields[14..16],
        .initial_total_volume_m3 = layer_fields[16..18],
        .dry_soil_mass_megagrams = layer_fields[18..20],
        .total_root_density_m_per_m3 = layer_fields[20..22],
        .east_west_face_area_m2 = layer_fields[22..24],
        .north_south_face_area_m2 = layer_fields[24..26],
        .surface_boundary_depth_m = cell_fields[0..1],
        .initial_surface_boundary_depth_m = cell_fields[1..2],
    };
    try initialize(state, .{
        .cell_count = cells,
        .layer_capacity = layers,
        .first_active_layer = &.{0},
        .active_layer_count = &.{2},
        .boundary_depth_m = &.{ -0.1, 0.1, 0.4 },
        .horizontal_cell_width_m = &.{10},
        .vertical_cell_width_m = &.{20},
        .matrix_volume_fraction = &.{ 0.8, 0.7 },
        .bulk_density_megagrams_per_m3 = &.{ 1.2, 1.3 },
        .initial_bulk_density_megagrams_per_m3 = &.{ 0.0, 1.3 },
        .calculation_floor = 1e-15,
    });
    try std.testing.expectEqual(@as(f64, 0.0), state.macropore_fraction[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.layer_thickness_m[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 40), state.total_volume_m3[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 32), state.matrix_volume_m3[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 38.4), state.dry_soil_mass_megagrams[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.layer_bottom_depth_from_surface_m[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.35), state.layer_midpoint_depth_from_surface_m[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.1), state.surface_boundary_depth_m[0], 1e-15);
}

test "invalid late layer leaves all geometry unchanged" {
    var layer_fields = [_]f64{7.0} ** (13 * 2);
    var cell_fields = [_]f64{7.0} ** 2;
    const before = layer_fields;
    const state: State = .{
        .macropore_fraction = layer_fields[0..2],
        .layer_thickness_m = layer_fields[2..4],
        .layer_midpoint_depth_m = layer_fields[4..6],
        .layer_bottom_depth_from_surface_m = layer_fields[6..8],
        .layer_midpoint_depth_from_surface_m = layer_fields[8..10],
        .total_volume_m3 = layer_fields[10..12],
        .matrix_volume_m3 = layer_fields[12..14],
        .reference_matrix_volume_m3 = layer_fields[14..16],
        .initial_total_volume_m3 = layer_fields[16..18],
        .dry_soil_mass_megagrams = layer_fields[18..20],
        .total_root_density_m_per_m3 = layer_fields[20..22],
        .east_west_face_area_m2 = layer_fields[22..24],
        .north_south_face_area_m2 = layer_fields[24..26],
        .surface_boundary_depth_m = cell_fields[0..1],
        .initial_surface_boundary_depth_m = cell_fields[1..2],
    };
    try std.testing.expectError(
        error.InvalidMineralLayerGeometryInput,
        initialize(state, .{
            .cell_count = 1,
            .layer_capacity = 2,
            .first_active_layer = &.{0},
            .active_layer_count = &.{2},
            .boundary_depth_m = &.{ 0, 0.1, 0.3 },
            .horizontal_cell_width_m = &.{10},
            .vertical_cell_width_m = &.{20},
            .matrix_volume_fraction = &.{ 0.8, 1.2 },
            .bulk_density_megagrams_per_m3 = &.{ 1.2, 1.3 },
            .initial_bulk_density_megagrams_per_m3 = &.{ 1.2, 1.3 },
            .calculation_floor = 1e-15,
        }),
    );
    try std.testing.expectEqualSlices(f64, &before, &layer_fields);
}
