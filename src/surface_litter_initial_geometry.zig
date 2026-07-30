const std = @import("std");

pub const Inputs = struct {
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
    /// Pool-fastest `[cell][residue_pool]`, in g C m-2.
    residue_carbon_g_c_per_m2: []const f64,
    residue_dry_bulk_density_Mg_per_m3: []const f64,
    residue_pool_count: usize,
    carbon_mass_conversion_Mg_per_g_c: f64,
};

pub const State = struct {
    horizontal_width_m: []f64,
    vertical_width_m: []f64,
    horizontal_area_m2: []f64,
    surface_depth_m: []f64,
    organic_carbon_g_c: []f64,
    initial_organic_carbon_g_c: []f64,
    dry_litter_volume_m3: []f64,
    total_volume_m3: []f64,
    matrix_volume_m3: []f64,
    reference_matrix_volume_m3: []f64,
    dry_mass_Mg: []f64,
    thickness_m: []f64,
    initial_thickness_m: []f64,
    total_surface_area_m2: *f64,
};

const Candidate = struct {
    area_m2: f64,
    organic_carbon_g_c: f64,
    dry_litter_volume_m3: f64,
    dry_mass_Mg: f64,
    thickness_m: f64,
};

fn candidate(inputs: Inputs, cell: usize) Candidate {
    const area_m2 = inputs.horizontal_cell_width_m[cell] *
        inputs.vertical_cell_width_m[cell];
    var organic_carbon_g_c_per_m2: f64 = 0.0;
    var dry_litter_volume_m3_per_m2: f64 = 0.0;
    for (0..inputs.residue_pool_count) |pool| {
        const index = cell * inputs.residue_pool_count + pool;
        organic_carbon_g_c_per_m2 +=
            inputs.residue_carbon_g_c_per_m2[index];
        dry_litter_volume_m3_per_m2 +=
            inputs.residue_carbon_g_c_per_m2[index] * 1.0e-6 /
            inputs.residue_dry_bulk_density_Mg_per_m3[pool];
    }
    const organic_carbon_g_c = organic_carbon_g_c_per_m2 * area_m2;
    const dry_litter_volume_m3 =
        dry_litter_volume_m3_per_m2 * area_m2;
    const dry_mass_Mg =
        inputs.carbon_mass_conversion_Mg_per_g_c * organic_carbon_g_c;
    return .{
        .area_m2 = area_m2,
        .organic_carbon_g_c = organic_carbon_g_c,
        .dry_litter_volume_m3 = dry_litter_volume_m3,
        .dry_mass_Mg = dry_mass_Mg,
        .thickness_m = dry_litter_volume_m3 / area_m2,
    };
}

/// Exact surface-layer branch of legacy `STARTS` lines 545--566.
pub fn initialize(state: State, inputs: Inputs) !void {
    const cell_count = inputs.horizontal_cell_width_m.len;
    if (cell_count == 0 or inputs.residue_pool_count == 0)
        return error.InvalidSurfaceLitterInitialDimensions;
    if (inputs.vertical_cell_width_m.len != cell_count or
        inputs.residue_dry_bulk_density_Mg_per_m3.len !=
            inputs.residue_pool_count)
        return error.SurfaceLitterInitialDimensionMismatch;
    const carbon_count = std.math.mul(
        usize,
        cell_count,
        inputs.residue_pool_count,
    ) catch return error.DimensionOverflow;
    if (inputs.residue_carbon_g_c_per_m2.len != carbon_count)
        return error.SurfaceLitterInitialDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == []f64 and @field(state, field.name).len != cell_count)
            return error.SurfaceLitterInitialDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.carbon_mass_conversion_Mg_per_g_c) or
        inputs.carbon_mass_conversion_Mg_per_g_c < 0)
        return error.InvalidSurfaceLitterInitialParameter;
    for (inputs.residue_dry_bulk_density_Mg_per_m3) |density| {
        if (!std.math.isFinite(density))
            return error.NonFiniteSurfaceLitterInitialInput;
        if (density <= 0) return error.InvalidSurfaceLitterInitialParameter;
    }
    for (inputs.residue_carbon_g_c_per_m2) |carbon| {
        if (!std.math.isFinite(carbon))
            return error.NonFiniteSurfaceLitterInitialInput;
        if (carbon < 0) return error.InvalidSurfaceLitterInitialInput;
    }

    var total_surface_area_m2: f64 = 0.0;
    for (0..cell_count) |cell| {
        inline for (.{
            inputs.horizontal_cell_width_m[cell],
            inputs.vertical_cell_width_m[cell],
        }) |width_m| {
            if (!std.math.isFinite(width_m))
                return error.NonFiniteSurfaceLitterInitialInput;
            if (width_m <= 0) return error.InvalidSurfaceLitterInitialInput;
        }
        const values = candidate(inputs, cell);
        inline for (@typeInfo(Candidate).@"struct".fields) |field| {
            if (!std.math.isFinite(@field(values, field.name)))
                return error.SurfaceLitterInitialGeometryOverflow;
        }
        total_surface_area_m2 += values.area_m2;
        if (!std.math.isFinite(total_surface_area_m2))
            return error.SurfaceLitterInitialGeometryOverflow;
    }

    state.total_surface_area_m2.* = 0.0;
    for (0..cell_count) |cell| {
        const values = candidate(inputs, cell);
        state.horizontal_width_m[cell] = inputs.horizontal_cell_width_m[cell];
        state.vertical_width_m[cell] = inputs.vertical_cell_width_m[cell];
        state.horizontal_area_m2[cell] = values.area_m2;
        state.total_surface_area_m2.* += values.area_m2;
        state.surface_depth_m[cell] = 0.0;
        state.organic_carbon_g_c[cell] = values.organic_carbon_g_c;
        state.initial_organic_carbon_g_c[cell] = values.organic_carbon_g_c;
        state.dry_litter_volume_m3[cell] = values.dry_litter_volume_m3;
        state.total_volume_m3[cell] = values.dry_litter_volume_m3;
        state.matrix_volume_m3[cell] = values.dry_litter_volume_m3;
        state.reference_matrix_volume_m3[cell] = values.dry_litter_volume_m3;
        state.dry_mass_Mg[cell] = values.dry_mass_Mg;
        state.thickness_m[cell] = values.thickness_m;
        state.initial_thickness_m[cell] = values.thickness_m;
    }
}

test "STARTS surface litter geometry reproduces three-pool source branch" {
    var fields = [_]f64{9.0} ** 13;
    var total_area: f64 = 9.0;
    const state: State = .{
        .horizontal_width_m = fields[0..1],
        .vertical_width_m = fields[1..2],
        .horizontal_area_m2 = fields[2..3],
        .surface_depth_m = fields[3..4],
        .organic_carbon_g_c = fields[4..5],
        .initial_organic_carbon_g_c = fields[5..6],
        .dry_litter_volume_m3 = fields[6..7],
        .total_volume_m3 = fields[7..8],
        .matrix_volume_m3 = fields[8..9],
        .reference_matrix_volume_m3 = fields[9..10],
        .dry_mass_Mg = fields[10..11],
        .thickness_m = fields[11..12],
        .initial_thickness_m = fields[12..13],
        .total_surface_area_m2 = &total_area,
    };
    try initialize(state, .{
        .horizontal_cell_width_m = &.{10},
        .vertical_cell_width_m = &.{20},
        .residue_carbon_g_c_per_m2 = &.{ 10, 20, 30 },
        .residue_dry_bulk_density_Mg_per_m3 = &.{ 0.1, 0.0125, 0.025 },
        .residue_pool_count = 3,
        .carbon_mass_conversion_Mg_per_g_c = 1.82e-6,
    });

    const expected_carbon_g_c: f64 = (10 + 20 + 30) * 200;
    const expected_volume_m3: f64 =
        (10.0e-6 / 0.1 + 20.0e-6 / 0.0125 +
            30.0e-6 / 0.025) * 200;
    try std.testing.expectEqual(@as(f64, 200), total_area);
    try std.testing.expectEqual(expected_carbon_g_c, state.organic_carbon_g_c[0]);
    try std.testing.expectApproxEqAbs(
        expected_volume_m3,
        state.dry_litter_volume_m3[0],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        expected_volume_m3 / 200,
        state.thickness_m[0],
        1e-15,
    );
}

test "runtime cell and residue-pool counts determine all extents" {
    var fields = [_]f64{0.0} ** (13 * 2);
    var total_area: f64 = 0.0;
    const state: State = .{
        .horizontal_width_m = fields[0..2],
        .vertical_width_m = fields[2..4],
        .horizontal_area_m2 = fields[4..6],
        .surface_depth_m = fields[6..8],
        .organic_carbon_g_c = fields[8..10],
        .initial_organic_carbon_g_c = fields[10..12],
        .dry_litter_volume_m3 = fields[12..14],
        .total_volume_m3 = fields[14..16],
        .matrix_volume_m3 = fields[16..18],
        .reference_matrix_volume_m3 = fields[18..20],
        .dry_mass_Mg = fields[20..22],
        .thickness_m = fields[22..24],
        .initial_thickness_m = fields[24..26],
        .total_surface_area_m2 = &total_area,
    };
    try initialize(state, .{
        .horizontal_cell_width_m = &.{ 10, 20 },
        .vertical_cell_width_m = &.{ 30, 40 },
        .residue_carbon_g_c_per_m2 = &.{ 1, 2, 3, 4, 5, 6, 7, 8 },
        .residue_dry_bulk_density_Mg_per_m3 = &.{ 0.1, 0.2, 0.3, 0.4 },
        .residue_pool_count = 4,
        .carbon_mass_conversion_Mg_per_g_c = 2e-6,
    });
    try std.testing.expectEqual(@as(f64, 1100), total_area);
    try std.testing.expectEqual(@as(f64, 26 * 800), state.organic_carbon_g_c[1]);
}
