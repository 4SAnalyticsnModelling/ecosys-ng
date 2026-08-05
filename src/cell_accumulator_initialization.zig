const std = @import("std");

pub const BoundaryLedgers = struct {
    precipitation_input_m3: []f64,
    carbon_dioxide_exchange_g_c: []f64,
    methane_exchange_g_c: []f64,
    oxygen_exchange_g_o: []f64,
    dinitrogen_exchange_g_n: []f64,
    nitrous_oxide_exchange_g_n: []f64,
    ammonia_exchange_g_n: []f64,
    soil_dinitrogen_fixation_g_n: []f64,
    organic_fertilizer_input_g_c: []f64,
    fertilizer_input_g_n: []f64,
    fertilizer_input_g_p: []f64,
    water_output_m3: []f64,
    evaporation_output_m3: []f64,
    runoff_output_m3: []f64,
    sediment_output_g: []f64,
    crop_output_g_c: []f64,
    salt_output_mol: []f64,
};

pub const DissolvedOutputs = struct {
    organic_runoff_g_c: []f64,
    organic_drainage_g_c: []f64,
    organic_runoff_g_n: []f64,
    organic_drainage_g_n: []f64,
    organic_runoff_g_p: []f64,
    organic_drainage_g_p: []f64,
    inorganic_runoff_g_c: []f64,
    inorganic_drainage_g_c: []f64,
    inorganic_runoff_g_n: []f64,
    inorganic_drainage_g_n: []f64,
    inorganic_runoff_g_p: []f64,
    inorganic_drainage_g_p: []f64,
};

pub const PlantAndDrainageState = struct {
    plant_sink_g_c: []f64,
    plant_sink_g_n: []f64,
    plant_sink_g_p: []f64,
    below_root_water_drainage_m3: []f64,
    below_root_nitrogen_drainage_g_n: []f64,
    below_root_phosphorus_drainage_g_p: []f64,
    ammonium_band_depth_m: []f64,
    nitrate_band_depth_m: []f64,
    phosphate_band_depth_m: []f64,
    phosphate_band_upper_depth_m: []f64,
};

pub const CanopyDiagnostics = struct {
    ground_radiation_fraction: []f64,
    net_radiation_megajoules_per_h: []f64,
    latent_heat_flux_megajoules_per_h: []f64,
    sensible_heat_flux_megajoules_per_h: []f64,
    storage_heat_flux_megajoules_per_h: []f64,
    evapotranspiration_m3_per_h: []f64,
    net_co2_exchange_g_c_per_h: []f64,
    surface_water_m3: []f64,
    living_leaf_area_m2: []f64,
    living_stalk_area_m2: []f64,
    standing_dead_area_m2: []f64,
    total_canopy_area_m2: []f64,
    canopy_water_retention_m3_per_h: []f64,
    biome_population_count: []f64,
    day_length_h: []f64,
    ground_albedo: []f64,
    harvest_g_c: []f64,
    harvest_g_n: []f64,
    harvest_g_p: []f64,
    rainfall_impact_energy_j: []f64,
};

pub const State = struct {
    boundary: BoundaryLedgers,
    dissolved: DissolvedOutputs,
    plant_and_drainage: PlantAndDrainageState,
    canopy: CanopyDiagnostics,
};

fn validateSlices(value: anytype, cell_count: usize) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (@field(value, field.name).len != cell_count) return false;
    }
    return true;
}

/// Exact source-order translation of legacy `STARTS` lines 414--472.
pub fn initialize(
    state: State,
    surface_albedo: []const f64,
    cell_count: usize,
) !void {
    if (cell_count == 0) return error.InvalidCellAccumulatorDimensions;
    if (surface_albedo.len != cell_count or
        !validateSlices(state.boundary, cell_count) or
        !validateSlices(state.dissolved, cell_count) or
        !validateSlices(state.plant_and_drainage, cell_count) or
        !validateSlices(state.canopy, cell_count))
    {
        return error.CellAccumulatorDimensionMismatch;
    }
    for (surface_albedo) |albedo| {
        if (!std.math.isFinite(albedo))
            return error.NonFiniteSurfaceAlbedo;
        if (albedo < 0 or albedo > 1) return error.InvalidSurfaceAlbedo;
    }

    for (0..cell_count) |cell| {
        state.boundary.precipitation_input_m3[cell] = 0.0;
        state.boundary.carbon_dioxide_exchange_g_c[cell] = 0.0;
        state.boundary.methane_exchange_g_c[cell] = 0.0;
        state.boundary.oxygen_exchange_g_o[cell] = 0.0;
        state.boundary.dinitrogen_exchange_g_n[cell] = 0.0;
        state.boundary.nitrous_oxide_exchange_g_n[cell] = 0.0;
        state.boundary.ammonia_exchange_g_n[cell] = 0.0;
        state.boundary.soil_dinitrogen_fixation_g_n[cell] = 0.0;
        state.boundary.organic_fertilizer_input_g_c[cell] = 0.0;
        state.boundary.fertilizer_input_g_n[cell] = 0.0;
        state.boundary.fertilizer_input_g_p[cell] = 0.0;
        state.boundary.water_output_m3[cell] = 0.0;
        state.boundary.evaporation_output_m3[cell] = 0.0;
        state.boundary.runoff_output_m3[cell] = 0.0;
        state.boundary.sediment_output_g[cell] = 0.0;
        state.boundary.crop_output_g_c[cell] = 0.0;
        state.dissolved.organic_runoff_g_c[cell] = 0.0;
        state.dissolved.organic_drainage_g_c[cell] = 0.0;
        state.dissolved.organic_runoff_g_n[cell] = 0.0;
        state.dissolved.organic_drainage_g_n[cell] = 0.0;
        state.dissolved.organic_runoff_g_p[cell] = 0.0;
        state.dissolved.organic_drainage_g_p[cell] = 0.0;
        state.dissolved.inorganic_runoff_g_c[cell] = 0.0;
        state.dissolved.inorganic_drainage_g_c[cell] = 0.0;
        state.dissolved.inorganic_runoff_g_n[cell] = 0.0;
        state.dissolved.inorganic_drainage_g_n[cell] = 0.0;
        state.dissolved.inorganic_runoff_g_p[cell] = 0.0;
        state.dissolved.inorganic_drainage_g_p[cell] = 0.0;
        state.boundary.salt_output_mol[cell] = 0.0;
        state.plant_and_drainage.plant_sink_g_c[cell] = 0.0;
        state.plant_and_drainage.plant_sink_g_n[cell] = 0.0;
        state.plant_and_drainage.plant_sink_g_p[cell] = 0.0;
        state.plant_and_drainage.below_root_water_drainage_m3[cell] = 0.0;
        state.plant_and_drainage.below_root_nitrogen_drainage_g_n[cell] = 0.0;
        state.plant_and_drainage.below_root_phosphorus_drainage_g_p[cell] = 0.0;
        state.plant_and_drainage.ammonium_band_depth_m[cell] = 0.0;
        state.plant_and_drainage.nitrate_band_depth_m[cell] = 0.0;
        state.plant_and_drainage.phosphate_band_depth_m[cell] = 0.0;
        state.plant_and_drainage.phosphate_band_upper_depth_m[cell] = 0.0;
        state.canopy.ground_radiation_fraction[cell] = 1.0;
        state.canopy.net_radiation_megajoules_per_h[cell] = 0.0;
        state.canopy.latent_heat_flux_megajoules_per_h[cell] = 0.0;
        state.canopy.sensible_heat_flux_megajoules_per_h[cell] = 0.0;
        state.canopy.storage_heat_flux_megajoules_per_h[cell] = 0.0;
        state.canopy.evapotranspiration_m3_per_h[cell] = 0.0;
        state.canopy.net_co2_exchange_g_c_per_h[cell] = 0.0;
        state.canopy.surface_water_m3[cell] = 0.0;
        state.canopy.living_leaf_area_m2[cell] = 0.0;
        state.canopy.living_stalk_area_m2[cell] = 0.0;
        state.canopy.standing_dead_area_m2[cell] = 0.0;
        state.canopy.total_canopy_area_m2[cell] = 0.0;
        state.canopy.canopy_water_retention_m3_per_h[cell] = 0.0;
        state.canopy.biome_population_count[cell] = 0.0;
        state.canopy.day_length_h[cell] = 12.0;
        state.canopy.ground_albedo[cell] = surface_albedo[cell];
        state.canopy.harvest_g_c[cell] = 0.0;
        state.canopy.harvest_g_n[cell] = 0.0;
        state.canopy.harvest_g_p[cell] = 0.0;
        state.canopy.rainfall_impact_energy_j[cell] = 0.0;
    }
}

fn filledSlices(comptime T: type, allocator: std.mem.Allocator, value: f64, count: usize) !T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| {
        @field(result, field.name) = try allocator.alloc(f64, count);
        @memset(@field(result, field.name), value);
    }
    return result;
}

fn freeSlices(value: anytype, allocator: std.mem.Allocator) void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        allocator.free(@field(value, field.name));
    }
}

test "STARTS initializes every per-cell accumulator with source values" {
    const allocator = std.testing.allocator;
    const cells = 3;
    const state: State = .{
        .boundary = try filledSlices(BoundaryLedgers, allocator, 7, cells),
        .dissolved = try filledSlices(DissolvedOutputs, allocator, 7, cells),
        .plant_and_drainage = try filledSlices(PlantAndDrainageState, allocator, 7, cells),
        .canopy = try filledSlices(CanopyDiagnostics, allocator, 7, cells),
    };
    defer freeSlices(state.boundary, allocator);
    defer freeSlices(state.dissolved, allocator);
    defer freeSlices(state.plant_and_drainage, allocator);
    defer freeSlices(state.canopy, allocator);

    try initialize(state, &.{ 0.1, 0.2, 0.3 }, cells);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 1 }, state.canopy.ground_radiation_fraction);
    try std.testing.expectEqualSlices(f64, &.{ 12, 12, 12 }, state.canopy.day_length_h);
    try std.testing.expectEqualSlices(f64, &.{ 0.1, 0.2, 0.3 }, state.canopy.ground_albedo);
    inline for (@typeInfo(BoundaryLedgers).@"struct".fields) |field| {
        for (@field(state.boundary, field.name)) |value|
            try std.testing.expectEqual(@as(f64, 0.0), value);
    }
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, state.canopy.rainfall_impact_energy_j);
}

test "invalid late albedo preserves every destination" {
    const allocator = std.testing.allocator;
    const cells = 2;
    const state: State = .{
        .boundary = try filledSlices(BoundaryLedgers, allocator, 7, cells),
        .dissolved = try filledSlices(DissolvedOutputs, allocator, 7, cells),
        .plant_and_drainage = try filledSlices(PlantAndDrainageState, allocator, 7, cells),
        .canopy = try filledSlices(CanopyDiagnostics, allocator, 7, cells),
    };
    defer freeSlices(state.boundary, allocator);
    defer freeSlices(state.dissolved, allocator);
    defer freeSlices(state.plant_and_drainage, allocator);
    defer freeSlices(state.canopy, allocator);

    try std.testing.expectError(
        error.NonFiniteSurfaceAlbedo,
        initialize(state, &.{ 0.2, std.math.nan(f64) }, cells),
    );
    inline for (@typeInfo(CanopyDiagnostics).@"struct".fields) |field| {
        for (@field(state.canopy, field.name)) |value|
            try std.testing.expectEqual(@as(f64, 7.0), value);
    }
    try std.testing.expectEqual(@as(f64, 7.0), state.boundary.precipitation_input_m3[0]);
}
