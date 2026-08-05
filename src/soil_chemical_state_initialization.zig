const std = @import("std");

pub const LayerKind = enum { surface_litter, mineral_soil };

pub const FertilizerInventory = struct {
    ammonium_g_n: *f64,
    ammonia_g_n: *f64,
    urea_g_n: *f64,
    nitrate_g_n: *f64,
};

pub const FertilizerBandGeometry = struct {
    ammonium_width_m: *f64,
    ammonium_depth_m: *f64,
    nitrate_width_m: *f64,
    nitrate_depth_m: *f64,
    phosphate_width_m: *f64,
    phosphate_depth_m: *f64,
};

pub const FertilizerVolumeFractions = struct {
    nonband_ammonium_fraction: *f64,
    nonband_nitrate_fraction: *f64,
    nonband_phosphate_fraction: *f64,
    band_ammonium_fraction: *f64,
    band_nitrate_fraction: *f64,
    band_phosphate_fraction: *f64,
};

pub const ReactionDemandState = struct {
    oxygen_g_o_per_h: *f64,
    ammonium_g_n_per_h: *f64,
    nitrate_g_n_per_h: *f64,
    nitrite_g_n_per_h: *f64,
    nitrous_oxide_g_n_per_h: *f64,
    h2po4_g_p_per_h: *f64,
    hpo4_g_p_per_h: *f64,
    microbial_carbon_g_c_per_h: *f64,
    band_ammonium_g_n_per_h: *f64,
    band_ammonia_g_n_per_h: *f64,
    band_nitrite_g_n_per_h: *f64,
    band_h2po4_g_p_per_h: *f64,
    band_hpo4_g_p_per_h: *f64,
    band_microbial_carbon_g_c_per_h: *f64,
};

pub const SubsurfaceBoundaryOrganicState = struct {
    carbon_g_c: []f64,
    nitrogen_g_n: []f64,
    phosphorus_g_p: []f64,
    acetate_g_c: []f64,
};

pub const InhibitorActivity = struct {
    current_urea_hydrolysis_fraction: *f64,
    initial_urea_hydrolysis_fraction: *f64,
    current_nitrification_fraction: *f64,
    initial_nitrification_fraction: *f64,
};

pub const State = struct {
    broadcast_fertilizer: FertilizerInventory,
    band_fertilizer: FertilizerInventory,
    band_geometry: FertilizerBandGeometry,
    volume_fraction: FertilizerVolumeFractions,
    reaction_demand: ReactionDemandState,
    subsurface_boundary_organic: SubsurfaceBoundaryOrganicState,
    inhibitor_activity: InhibitorActivity,
    layer_combustion_heat_megajoules: *f64,
    cell_combustion_heat_megajoules: *f64,
};

fn zeroPointerFields(value: anytype) void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        @field(value, field.name).* = 0.0;
}

/// Exact source-order translation of legacy `STARTS` lines 1589--1642.
///
/// Fertilizer inventories use g N, band geometry uses m, reaction demands
/// use g element/h, boundary pools use g C/N/P, and heat uses MJ.
pub fn initialize(
    state: State,
    layer_kind: LayerKind,
    organic_pool_count: usize,
) !void {
    if (organic_pool_count == 0)
        return error.InvalidSoilChemicalInitialDimensions;
    inline for (@typeInfo(SubsurfaceBoundaryOrganicState).@"struct".fields) |field|
        if (@field(state.subsurface_boundary_organic, field.name).len !=
            organic_pool_count)
            return error.SoilChemicalInitialDimensionMismatch;

    zeroPointerFields(state.broadcast_fertilizer);
    if (layer_kind == .mineral_soil) {
        zeroPointerFields(state.band_fertilizer);
        zeroPointerFields(state.band_geometry);
    }

    state.volume_fraction.nonband_ammonium_fraction.* = 1.0;
    state.volume_fraction.nonband_nitrate_fraction.* = 1.0;
    state.volume_fraction.nonband_phosphate_fraction.* = 1.0;
    state.volume_fraction.band_ammonium_fraction.* = 0.0;
    state.volume_fraction.band_nitrate_fraction.* = 0.0;
    state.volume_fraction.band_phosphate_fraction.* = 0.0;

    zeroPointerFields(state.reaction_demand);
    if (layer_kind == .mineral_soil) {
        inline for (@typeInfo(SubsurfaceBoundaryOrganicState).@"struct".fields) |field|
            @memset(
                @field(state.subsurface_boundary_organic, field.name),
                0.0,
            );
    }
    zeroPointerFields(state.inhibitor_activity);
    state.layer_combustion_heat_megajoules.* = 0.0;
    state.cell_combustion_heat_megajoules.* = 0.0;
}

fn pointerStruct(comptime T: type, values: []f64, start: usize) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields, 0..) |field, index|
        @field(result, field.name) = &values[start + index];
    return result;
}

const fertilizer_field_count = @typeInfo(FertilizerInventory).@"struct".fields.len;
const geometry_field_count =
    @typeInfo(FertilizerBandGeometry).@"struct".fields.len;
const volume_field_count =
    @typeInfo(FertilizerVolumeFractions).@"struct".fields.len;
const demand_field_count =
    @typeInfo(ReactionDemandState).@"struct".fields.len;
const inhibitor_field_count =
    @typeInfo(InhibitorActivity).@"struct".fields.len;
const scalar_count = 2 * fertilizer_field_count + geometry_field_count +
    volume_field_count + demand_field_count + inhibitor_field_count + 2;

fn testState(values: []f64, boundary: []f64) State {
    var offset: usize = 0;
    const broadcast =
        pointerStruct(FertilizerInventory, values, offset);
    offset += fertilizer_field_count;
    const band = pointerStruct(FertilizerInventory, values, offset);
    offset += fertilizer_field_count;
    const geometry =
        pointerStruct(FertilizerBandGeometry, values, offset);
    offset += geometry_field_count;
    const volume =
        pointerStruct(FertilizerVolumeFractions, values, offset);
    offset += volume_field_count;
    const demand = pointerStruct(ReactionDemandState, values, offset);
    offset += demand_field_count;
    const inhibitor = pointerStruct(InhibitorActivity, values, offset);
    offset += inhibitor_field_count;
    return .{
        .broadcast_fertilizer = broadcast,
        .band_fertilizer = band,
        .band_geometry = geometry,
        .volume_fraction = volume,
        .reaction_demand = demand,
        .subsurface_boundary_organic = .{
            .carbon_g_c = boundary[0..2],
            .nitrogen_g_n = boundary[2..4],
            .phosphorus_g_p = boundary[4..6],
            .acetate_g_c = boundary[6..8],
        },
        .inhibitor_activity = inhibitor,
        .layer_combustion_heat_megajoules = &values[offset],
        .cell_combustion_heat_megajoules = &values[offset + 1],
    };
}

test "STARTS mineral layer resets all chemical initialization state" {
    var values = [_]f64{9.0} ** scalar_count;
    var boundary = [_]f64{9.0} ** 8;
    try initialize(testState(&values, &boundary), .mineral_soil, 2);

    for (values[0 .. 2 * fertilizer_field_count + geometry_field_count]) |value|
        try std.testing.expectEqual(@as(f64, 0.0), value);
    const volume_start =
        2 * fertilizer_field_count + geometry_field_count;
    try std.testing.expectEqualSlices(
        f64,
        &.{ 1, 1, 1, 0, 0, 0 },
        values[volume_start .. volume_start + volume_field_count],
    );
    for (values[volume_start + volume_field_count ..]) |value|
        try std.testing.expectEqual(@as(f64, 0.0), value);
    for (boundary) |value|
        try std.testing.expectEqual(@as(f64, 0.0), value);
}

test "STARTS surface layer preserves mineral-only band and boundary state" {
    var values = [_]f64{9.0} ** scalar_count;
    var boundary = [_]f64{9.0} ** 8;
    try initialize(testState(&values, &boundary), .surface_litter, 2);

    for (values[0..fertilizer_field_count]) |value|
        try std.testing.expectEqual(@as(f64, 0.0), value);
    for (values[fertilizer_field_count .. 2 * fertilizer_field_count + geometry_field_count]) |value|
        try std.testing.expectEqual(@as(f64, 9.0), value);
    for (boundary) |value|
        try std.testing.expectEqual(@as(f64, 9.0), value);
    const volume_start =
        2 * fertilizer_field_count + geometry_field_count;
    try std.testing.expectEqualSlices(
        f64,
        &.{ 1, 1, 1, 0, 0, 0 },
        values[volume_start .. volume_start + volume_field_count],
    );
}

test "boundary dimension failure is atomic" {
    var values = [_]f64{9.0} ** scalar_count;
    var boundary = [_]f64{9.0} ** 8;
    try std.testing.expectError(
        error.SoilChemicalInitialDimensionMismatch,
        initialize(testState(&values, &boundary), .mineral_soil, 3),
    );
    const expected_values = [_]f64{9.0} ** scalar_count;
    const expected_boundary = [_]f64{9.0} ** 8;
    try std.testing.expectEqualSlices(f64, &expected_values, &values);
    try std.testing.expectEqualSlices(f64, &expected_boundary, &boundary);
}
