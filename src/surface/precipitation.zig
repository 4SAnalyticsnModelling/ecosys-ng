const std = @import("std");
const builtin = @import("builtin");
const water_flux = @import("../soil/water/flux.zig");
const surface_water_flow = @import("water_flow.zig");
const retention = @import("../soil/water/retention.zig");
const snow_cover = @import("../soil/water/snow_cover_fraction.zig");

pub const RuntimeState = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    litter_water_m3: []f64,
    litter_water_capacity_m3: []f64,
    litter_cover_fraction: []f64,
    solid_snow_water_equivalent_m3: []f64,
    rain_to_snow_m3_per_h: []f64,
    snow_to_snow_m3_per_h: []f64,
    water_to_litter_m3_per_h: []f64,
    water_to_matrix_m3_per_h: []f64,
    water_to_macropore_m3_per_h: []f64,
    heat_to_snow_megajoules_per_h: []f64,
    heat_to_litter_megajoules_per_h: []f64,
    heat_to_soil_megajoules_per_h: []f64,
    rainfall_m3_per_h: []f64,
    snowfall_water_equivalent_m3_per_h: []f64,
    intercepted_rain_m3_per_h: []f64,
    snow_cover_fraction: []f64,
    atmospheric_temperature_k: []f64,
    matrix_fraction: []f64,
    macropore_fraction: []f64,
    matrix_air_capacity_m3: []f64,
    macropore_air_capacity_m3: []f64,
    litter_absent_above_water_table: []bool,
    rainfall_impact_energy_j: []f64,
    cumulative_rainfall_impact_energy_j: []f64,
    saturated_hydraulic_conductivity_multiplier: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !RuntimeState {
        if (cell_count == 0) return error.InvalidSurfaceIngressDimensions;
        var result: RuntimeState = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(RuntimeState).@"struct".fields) |field| if ((field.type == []f64 or field.type == []bool) and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(RuntimeState).@"struct".fields) |field| if (field.type == []f64 or field.type == []bool) {
            @field(result, field.name) = try allocator.alloc(std.meta.Elem(field.type), cell_count);
            @memset(@field(result, field.name), if (field.type == []bool) false else 0);
            allocated += 1;
        };
        @memset(result.saturated_hydraulic_conductivity_multiplier, 1);
        return result;
    }

    pub fn deinit(self: *RuntimeState) void {
        inline for (@typeInfo(RuntimeState).@"struct".fields) |field| if (field.type == []f64 or field.type == []bool) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// Binds WTHR/HOUR1 extensive carriers to the runtime top-soil geometry and
/// then publishes the WATSUB surface ledgers without allocating per hour.
pub fn prepareFromModel(state: *RuntimeState, atmosphere: *const @import("../atmosphere/atmospheric_forcing.zig").State, grid: *const @import("../state/grid.zig").GridState, canopy: ?*const @import("../canopy/energy/precipitation_retention.zig").State, cell_area_m2: []const f64, snow_depth_m: []const f64, full_snow_cover_depth_m: f64, maximum_iterations: u16) !void {
    if (atmosphere.cell_count != state.cell_count or grid.cell_count != state.cell_count or cell_area_m2.len != state.cell_count or snow_depth_m.len != state.cell_count or !std.math.isFinite(full_snow_cover_depth_m) or full_snow_cover_depth_m <= 0 or maximum_iterations == 0) return error.SurfaceIngressDimensionMismatch;
    if (canopy) |value| if (value.cell_count != state.cell_count) return error.SurfaceIngressDimensionMismatch;
    for (0..state.cell_count) |cell| {
        const area = cell_area_m2[cell];
        const top = cell * grid.soil_layer_capacity;
        const total_pore_capacity = grid.matrix_pore_capacity_m3[top] + grid.macropore_pore_capacity_m3[top];
        if (!std.math.isFinite(area) or area <= 0 or !std.math.isFinite(snow_depth_m[cell]) or snow_depth_m[cell] < 0 or total_pore_capacity <= 0) return error.InvalidSurfaceIngressModelState;
        state.rainfall_m3_per_h[cell] = atmosphere.rainfall_m[cell] * area;
        state.snowfall_water_equivalent_m3_per_h[cell] = atmosphere.snowfall_water_equivalent_m[cell] * area;
        state.intercepted_rain_m3_per_h[cell] = if (canopy) |value| value.cell_retention_m3_per_h[cell] else 0;
        // WATSUB `FSNW`, `watsub.f` 386--392. Delegated to the single
        // authoritative owner. This site previously squared the depth ratio
        // where the source takes its SQUARE ROOT, and omitted the `FSNX`
        // `1.0e-3` floor, understating cover by up to 52x for thin snow. See
        // `docs/traceability/watsub_snow_cover_fraction_exponent_defect.md`.
        // Consumers derive the snow-free fraction as `1 - snow_cover_fraction`,
        // which is exact here because the owner returns the covered fraction as
        // the exact complement of the floored snow-free fraction.
        state.snow_cover_fraction[cell] = (try snow_cover.evaluate(snow_depth_m[cell], full_snow_cover_depth_m)).snow_fraction;
        state.atmospheric_temperature_k[cell] = atmosphere.air_temperature_k[cell];
        state.matrix_fraction[cell] = grid.matrix_pore_capacity_m3[top] / total_pore_capacity;
        state.macropore_fraction[cell] = grid.macropore_pore_capacity_m3[top] / total_pore_capacity;
        state.matrix_air_capacity_m3[cell] = @max(
            0,
            grid.matrix_pore_capacity_m3[top] -
                grid.matrix_liquid_water_m3[top] -
                grid.matrix_ice_water_m3[top],
        );
        state.macropore_air_capacity_m3[cell] = @max(
            0,
            grid.macropore_pore_capacity_m3[top] -
                grid.macropore_liquid_water_m3[top] -
                grid.macropore_ice_water_m3[top],
        );
    }
    try prepareRuntimeIngress(state, .{ .rainfall_m3_per_h = state.rainfall_m3_per_h, .snowfall_water_equivalent_m3_per_h = state.snowfall_water_equivalent_m3_per_h, .intercepted_rain_m3_per_h = state.intercepted_rain_m3_per_h, .snow_cover_fraction = state.snow_cover_fraction, .atmospheric_temperature_k = state.atmospheric_temperature_k, .matrix_fraction = state.matrix_fraction, .macropore_fraction = state.macropore_fraction, .matrix_air_capacity_m3 = state.matrix_air_capacity_m3, .macropore_air_capacity_m3 = state.macropore_air_capacity_m3, .litter_absent_above_water_table = state.litter_absent_above_water_table, .maximum_iterations = maximum_iterations });
}

pub const RuntimeInputs = struct {
    rainfall_m3_per_h: []const f64,
    snowfall_water_equivalent_m3_per_h: []const f64,
    intercepted_rain_m3_per_h: []const f64,
    snow_cover_fraction: []const f64,
    atmospheric_temperature_k: []const f64,
    matrix_fraction: []const f64,
    macropore_fraction: []const f64,
    matrix_air_capacity_m3: []const f64,
    macropore_air_capacity_m3: []const f64,
    litter_absent_above_water_table: []const bool,
    maximum_iterations: u16,
};

/// Prepares and atomically publishes WATSUB precipitation carriers for all
/// runtime cells. Soil/snow storage commit consumes these ledgers separately.
pub fn prepareRuntimeIngress(state: *RuntimeState, inputs: RuntimeInputs) !void {
    const count = state.cell_count;
    inline for (.{ inputs.rainfall_m3_per_h, inputs.snowfall_water_equivalent_m3_per_h, inputs.intercepted_rain_m3_per_h, inputs.snow_cover_fraction, inputs.atmospheric_temperature_k, inputs.matrix_fraction, inputs.macropore_fraction, inputs.matrix_air_capacity_m3, inputs.macropore_air_capacity_m3 }) |values| if (values.len != count) return error.SurfaceIngressDimensionMismatch;
    if (inputs.litter_absent_above_water_table.len != count or inputs.maximum_iterations == 0) return error.SurfaceIngressDimensionMismatch;
    // Validate every cell before publishing any ledger. Repeating these pure,
    // inexpensive calculations avoids both per-hour allocation and partial
    // flux publication when a later cell contains invalid data.
    for (0..count) |cell| _ = try calculateRuntimeIngressCell(state, inputs, cell);
    for (0..count) |cell| {
        const result = try calculateRuntimeIngressCell(state, inputs, cell);
        const partition = result.partition;
        const redistributed = result.redistributed;
        state.rain_to_snow_m3_per_h[cell] = partition.rain_to_snow_m3_per_h;
        state.snow_to_snow_m3_per_h[cell] = partition.snow_to_snow_m3_per_h;
        state.water_to_litter_m3_per_h[cell] = redistributed.litter_water_m3;
        state.water_to_matrix_m3_per_h[cell] = redistributed.matrix_water_m3;
        state.water_to_macropore_m3_per_h[cell] = redistributed.macropore_water_m3;
        state.heat_to_snow_megajoules_per_h[cell] = partition.heat_to_snow_megajoules_per_h;
        state.heat_to_litter_megajoules_per_h[cell] = redistributed.litter_heat_megajoules;
        state.heat_to_soil_megajoules_per_h[cell] = redistributed.soil_heat_megajoules;
    }
}

fn calculateRuntimeIngressCell(state: *const RuntimeState, inputs: RuntimeInputs, cell: usize) !struct { partition: SurfacePartition, redistributed: Redistribution } {
    const partition = try partitionAtmosphericWater(.{ .rain_and_irrigation_m3_per_h = inputs.rainfall_m3_per_h[cell], .intercepted_rain_m3_per_h = inputs.intercepted_rain_m3_per_h[cell], .snowfall_water_equivalent_m3_per_h = inputs.snowfall_water_equivalent_m3_per_h[cell], .snow_cover_fraction = inputs.snow_cover_fraction[cell], .snow_free_fraction = 1.0 - inputs.snow_cover_fraction[cell], .litter_cover_fraction = state.litter_cover_fraction[cell], .litter_water_capacity_m3 = state.litter_water_capacity_m3[cell], .litter_water_m3 = state.litter_water_m3[cell], .litter_absent_above_water_table = inputs.litter_absent_above_water_table[cell], .atmospheric_temperature_k = inputs.atmospheric_temperature_k[cell], .matrix_fraction = inputs.matrix_fraction[cell], .macropore_fraction = inputs.macropore_fraction[cell], .maximum_iterations = inputs.maximum_iterations });
    const redistributed = try redistribute(.{ .soil_surface_present = true, .rain_and_irrigation_to_matrix_m3 = partition.water_to_matrix_m3_per_h, .rain_and_irrigation_to_macropore_m3 = partition.water_to_macropore_m3_per_h, .rain_and_irrigation_to_litter_m3 = partition.water_to_litter_m3_per_h, .matrix_air_capacity_m3 = inputs.matrix_air_capacity_m3[cell], .macropore_air_capacity_m3 = inputs.macropore_air_capacity_m3[cell], .snow_free_fraction = 1.0 - inputs.snow_cover_fraction[cell], .atmospheric_temperature_k = inputs.atmospheric_temperature_k[cell], .litter_input_heat_megajoules = partition.heat_to_litter_megajoules_per_h, .soil_input_heat_megajoules = partition.heat_to_soil_megajoules_per_h });
    return .{ .partition = partition, .redistributed = redistributed };
}

pub fn commitRuntimeIngress(state: *RuntimeState, timestep_h: f64) !void {
    if (!std.math.isFinite(timestep_h) or timestep_h <= 0) return error.InvalidSurfaceIngressTimestep;
    for (0..state.cell_count) |cell| {
        const next_litter = state.litter_water_m3[cell] + state.water_to_litter_m3_per_h[cell] * timestep_h;
        const next_snow = state.solid_snow_water_equivalent_m3[cell] + (state.rain_to_snow_m3_per_h[cell] + state.snow_to_snow_m3_per_h[cell]) * timestep_h;
        if (!std.math.isFinite(next_litter) or !std.math.isFinite(next_snow) or next_litter < 0 or next_snow < 0) return error.InvalidSurfaceIngressCommit;
    }
    for (0..state.cell_count) |cell| {
        state.litter_water_m3[cell] += state.water_to_litter_m3_per_h[cell] * timestep_h;
        state.solid_snow_water_equivalent_m3[cell] += (state.rain_to_snow_m3_per_h[cell] + state.snow_to_snow_m3_per_h[cell]) * timestep_h;
    }
}

/// Atomically publishes the prepared liquid carriers into authoritative top
/// matrix/macropore storage. All cells validate before any storage changes.
pub fn commitSoilIngress(state: *const RuntimeState, grid: *@import("../state/grid.zig").GridState, hydrology: *@import("../transport/hydrology.zig").State, timestep_h: f64) !void {
    if (grid.cell_count != state.cell_count or hydrology.columns * hydrology.rows != state.cell_count or !std.math.isFinite(timestep_h) or timestep_h <= 0) return error.SurfaceIngressDimensionMismatch;
    for (0..state.cell_count) |cell| {
        const top = cell * grid.soil_layer_capacity;
        const matrix_input = state.water_to_matrix_m3_per_h[cell] * timestep_h;
        const macropore_input = state.water_to_macropore_m3_per_h[cell] * timestep_h;
        const matrix_air_m3 = grid.matrix_pore_capacity_m3[top] -
            grid.matrix_liquid_water_m3[top] -
            grid.matrix_ice_water_m3[top];
        const macropore_air_m3 = grid.macropore_pore_capacity_m3[top] -
            grid.macropore_liquid_water_m3[top] -
            grid.macropore_ice_water_m3[top];
        const matrix_roundoff_m3 =
            poreCapacityRoundoffToleranceM3(
                grid.matrix_pore_capacity_m3[top],
            );
        const macropore_roundoff_m3 =
            poreCapacityRoundoffToleranceM3(
                grid.macropore_pore_capacity_m3[top],
            );
        if (!std.math.isFinite(matrix_input) or
            !std.math.isFinite(macropore_input) or matrix_input < 0 or
            macropore_input < 0 or !std.math.isFinite(matrix_air_m3) or
            !std.math.isFinite(macropore_air_m3) or
            matrix_air_m3 < -matrix_roundoff_m3 or
            macropore_air_m3 < -macropore_roundoff_m3 or
            matrix_input >
                @max(0, matrix_air_m3) + matrix_roundoff_m3 or
            macropore_input >
                @max(0, macropore_air_m3) + macropore_roundoff_m3)
        {
            if (!builtin.is_test)
                std.log.err(
                    "invalid surface ingress commit: cell={} matrix_input_m3={} matrix_liquid_m3={} matrix_ice_m3={} matrix_pore_capacity_m3={} matrix_air_m3={} macropore_input_m3={} macropore_liquid_m3={} macropore_ice_m3={} macropore_pore_capacity_m3={} macropore_air_m3={}",
                    .{
                        cell,
                        matrix_input,
                        grid.matrix_liquid_water_m3[top],
                        grid.matrix_ice_water_m3[top],
                        grid.matrix_pore_capacity_m3[top],
                        matrix_air_m3,
                        macropore_input,
                        grid.macropore_liquid_water_m3[top],
                        grid.macropore_ice_water_m3[top],
                        grid.macropore_pore_capacity_m3[top],
                        macropore_air_m3,
                    },
                );
            return error.InvalidSurfaceIngressCommit;
        }
    }
    for (0..state.cell_count) |cell| {
        const top = cell * grid.soil_layer_capacity;
        const matrix_input = state.water_to_matrix_m3_per_h[cell] * timestep_h;
        const macropore_input = state.water_to_macropore_m3_per_h[cell] * timestep_h;
        grid.matrix_liquid_water_m3[top] += matrix_input;
        grid.macropore_liquid_water_m3[top] += macropore_input;
        grid.liquid_water_m3[top] += matrix_input + macropore_input;
        grid.matrix_air_volume_m3[top] = @max(
            0,
            grid.matrix_pore_capacity_m3[top] -
                grid.matrix_liquid_water_m3[top] -
                grid.matrix_ice_water_m3[top],
        );
        grid.macropore_air_volume_m3[top] = @max(
            0,
            grid.macropore_pore_capacity_m3[top] -
                grid.macropore_liquid_water_m3[top] -
                grid.macropore_ice_water_m3[top],
        );
        grid.air_volume_m3[top] = grid.matrix_air_volume_m3[top] + grid.macropore_air_volume_m3[top];
        hydrology.micropore_water_volume_m3[top] = grid.matrix_liquid_water_m3[top];
        hydrology.macropore_water_volume_m3[top] = grid.macropore_liquid_water_m3[top];
        hydrology.matrix_air_volume_m3[top] = grid.matrix_air_volume_m3[top];
        hydrology.macropore_air_volume_m3[top] = grid.macropore_air_volume_m3[top];
        hydrology.air_volume_m3[top] = grid.air_volume_m3[top];
    }
}

fn poreCapacityRoundoffToleranceM3(capacity_m3: f64) f64 {
    return @max(
        1.0e-12,
        64.0 * std.math.floatEps(f64) *
            @max(1.0, @abs(capacity_m3)),
    );
}

pub fn bindSoilHeatIngress(state: *const RuntimeState, grid: *const @import("../state/grid.zig").GridState, soil_heat_source_megajoules: []f64, timestep_h: f64) !void {
    if (grid.cell_count != state.cell_count or soil_heat_source_megajoules.len != grid.layer_count or !std.math.isFinite(timestep_h) or timestep_h <= 0) return error.SurfaceIngressDimensionMismatch;
    for (0..state.cell_count) |cell| {
        const heat = state.heat_to_soil_megajoules_per_h[cell] * timestep_h;
        if (!std.math.isFinite(heat)) return error.InvalidSurfaceIngressCommit;
        soil_heat_source_megajoules[cell * grid.soil_layer_capacity] += heat;
    }
}

pub const SurfacePartitionInputs = struct {
    rain_and_irrigation_m3_per_h: f64,
    intercepted_rain_m3_per_h: f64,
    snowfall_water_equivalent_m3_per_h: f64,
    snow_cover_fraction: f64,
    snow_free_fraction: f64,
    litter_cover_fraction: f64,
    litter_water_capacity_m3: f64,
    litter_water_m3: f64,
    litter_absent_above_water_table: bool,
    atmospheric_temperature_k: f64,
    matrix_fraction: f64,
    macropore_fraction: f64,
    maximum_iterations: u16,
};

pub const SurfacePartition = struct {
    rain_to_snow_m3_per_h: f64,
    snow_to_snow_m3_per_h: f64,
    heat_to_snow_megajoules_per_h: f64,
    water_to_litter_m3_per_h: f64,
    water_to_matrix_m3_per_h: f64,
    water_to_macropore_m3_per_h: f64,
    heat_to_litter_megajoules_per_h: f64,
    heat_to_soil_megajoules_per_h: f64,
};

/// WATSUB precipitation partition before the nonlinear soil solve. NPH is
/// retained only as the iteration ceiling used by the litter-capacity bound.
pub fn partitionAtmosphericWater(inputs: SurfacePartitionInputs) !SurfacePartition {
    inline for (@typeInfo(SurfacePartitionInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfacePrecipitationInput;
    // TFLWC is signed: canopy water above its current capacity drains and
    // appears as negative retention, increasing water delivered below canopy.
    if (inputs.rain_and_irrigation_m3_per_h < 0 or inputs.intercepted_rain_m3_per_h > inputs.rain_and_irrigation_m3_per_h or inputs.snowfall_water_equivalent_m3_per_h < 0 or inputs.snow_cover_fraction < 0 or inputs.snow_cover_fraction > 1 or inputs.snow_free_fraction < 0 or inputs.snow_free_fraction > 1 or inputs.litter_cover_fraction < 0 or inputs.litter_cover_fraction > 1 or inputs.litter_water_capacity_m3 < 0 or inputs.litter_water_m3 < 0 or inputs.atmospheric_temperature_k <= 0 or inputs.matrix_fraction < 0 or inputs.macropore_fraction < 0 or @abs(inputs.matrix_fraction + inputs.macropore_fraction - 1) > 1e-10 or inputs.maximum_iterations == 0) return error.InvalidSurfacePrecipitationInput;
    const rain_after_interception = inputs.rain_and_irrigation_m3_per_h - inputs.intercepted_rain_m3_per_h;
    const rain_to_snow = rain_after_interception * inputs.snow_cover_fraction;
    const surface_water = rain_after_interception * inputs.snow_free_fraction;
    const water_to_litter = if (inputs.litter_absent_above_water_table)
        surface_water
    else
        @max(0.0, @min(surface_water * inputs.litter_cover_fraction, (inputs.litter_water_capacity_m3 - inputs.litter_water_m3) * inputs.snow_free_fraction * @as(f64, @floatFromInt(inputs.maximum_iterations))));
    const water_to_soil = surface_water - water_to_litter;
    return .{
        .rain_to_snow_m3_per_h = rain_to_snow,
        .snow_to_snow_m3_per_h = inputs.snowfall_water_equivalent_m3_per_h,
        .heat_to_snow_megajoules_per_h = inputs.atmospheric_temperature_k * (2.095 * inputs.snowfall_water_equivalent_m3_per_h + 4.19 * rain_to_snow),
        .water_to_litter_m3_per_h = water_to_litter,
        .water_to_matrix_m3_per_h = water_to_soil * inputs.matrix_fraction,
        .water_to_macropore_m3_per_h = water_to_soil * inputs.macropore_fraction,
        .heat_to_litter_megajoules_per_h = 4.19 * inputs.atmospheric_temperature_k * water_to_litter,
        .heat_to_soil_megajoules_per_h = 4.19 * inputs.atmospheric_temperature_k * water_to_soil,
    };
}

pub const SolutePrecipitationRouting = struct { rain_to_litter_m3: f64, irrigation_to_litter_m3: f64, rain_to_snow_m3: f64, irrigation_to_snow_m3: f64 };

/// Exact FLQRQ/FLQRI/FLQGQ/FLQGI routing used by TRNSFR and TRNSFRS.
pub fn routePrecipitationSolutes(snowfall_m3_per_h: f64, rain_m3_per_h: f64, combined_rain_m3_per_h: f64, irrigation_m3_per_h: f64, water_to_litter_m3_per_h: f64, snow_heat_capacity_megajoules_per_k: f64, minimum_snow_heat_capacity_megajoules_per_k: f64, timestep_h: f64) !SolutePrecipitationRouting {
    inline for (.{ snowfall_m3_per_h, rain_m3_per_h, combined_rain_m3_per_h, irrigation_m3_per_h, water_to_litter_m3_per_h, snow_heat_capacity_megajoules_per_k, minimum_snow_heat_capacity_megajoules_per_k, timestep_h }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfacePrecipitationInput;
    if (snowfall_m3_per_h < 0 or rain_m3_per_h < 0 or combined_rain_m3_per_h < 0 or irrigation_m3_per_h < 0 or water_to_litter_m3_per_h < 0 or snow_heat_capacity_megajoules_per_k < 0 or minimum_snow_heat_capacity_megajoules_per_k < 0 or timestep_h <= 0) return error.InvalidSurfacePrecipitationInput;
    if (snowfall_m3_per_h > 0 or (rain_m3_per_h > 0 and snow_heat_capacity_megajoules_per_k > minimum_snow_heat_capacity_megajoules_per_k)) return .{ .rain_to_litter_m3 = 0, .irrigation_to_litter_m3 = 0, .rain_to_snow_m3 = combined_rain_m3_per_h * timestep_h, .irrigation_to_snow_m3 = irrigation_m3_per_h * timestep_h };
    const total_liquid = combined_rain_m3_per_h + irrigation_m3_per_h;
    if (total_liquid > 0 and snow_heat_capacity_megajoules_per_k <= minimum_snow_heat_capacity_megajoules_per_k) {
        const rain_to_litter = water_to_litter_m3_per_h * combined_rain_m3_per_h / total_liquid * timestep_h;
        const irrigation_to_litter = water_to_litter_m3_per_h * irrigation_m3_per_h / total_liquid * timestep_h;
        return .{ .rain_to_litter_m3 = rain_to_litter, .irrigation_to_litter_m3 = irrigation_to_litter, .rain_to_snow_m3 = combined_rain_m3_per_h * timestep_h - rain_to_litter, .irrigation_to_snow_m3 = irrigation_m3_per_h * timestep_h - irrigation_to_litter };
    }
    return .{ .rain_to_litter_m3 = 0, .irrigation_to_litter_m3 = 0, .rain_to_snow_m3 = 0, .irrigation_to_snow_m3 = 0 };
}

pub const RedistributionInputs = struct {
    soil_surface_present: bool,
    rain_and_irrigation_to_matrix_m3: f64,
    rain_and_irrigation_to_macropore_m3: f64,
    rain_and_irrigation_to_litter_m3: f64,
    matrix_air_capacity_m3: f64,
    macropore_air_capacity_m3: f64,
    snow_free_fraction: f64,
    atmospheric_temperature_k: f64,
    litter_input_heat_megajoules: f64,
    soil_input_heat_megajoules: f64,
};

pub const Redistribution = struct {
    litter_water_m3: f64,
    matrix_water_m3: f64,
    macropore_water_m3: f64,
    litter_heat_megajoules: f64,
    soil_heat_megajoules: f64,
};

pub fn redistribute(inputs: RedistributionInputs) !Redistribution {
    inline for (@typeInfo(RedistributionInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfacePrecipitationInput;
    if (inputs.rain_and_irrigation_to_matrix_m3 < 0 or inputs.rain_and_irrigation_to_macropore_m3 < 0 or inputs.rain_and_irrigation_to_litter_m3 < 0 or inputs.matrix_air_capacity_m3 < 0 or inputs.macropore_air_capacity_m3 < 0 or inputs.snow_free_fraction < 0 or inputs.snow_free_fraction > 1 or inputs.atmospheric_temperature_k <= 0) return error.InvalidSurfacePrecipitationInput;
    if (!inputs.soil_surface_present) return .{ .litter_water_m3 = inputs.rain_and_irrigation_to_litter_m3, .matrix_water_m3 = inputs.rain_and_irrigation_to_matrix_m3, .macropore_water_m3 = inputs.rain_and_irrigation_to_macropore_m3, .litter_heat_megajoules = inputs.litter_input_heat_megajoules, .soil_heat_megajoules = inputs.soil_input_heat_megajoules };
    const matrix_excess = @max(0.0, inputs.rain_and_irrigation_to_matrix_m3 - inputs.matrix_air_capacity_m3 * inputs.snow_free_fraction);
    const macro_excess = @max(0.0, inputs.rain_and_irrigation_to_macropore_m3 - inputs.macropore_air_capacity_m3 * inputs.snow_free_fraction);
    const redirected_heat = 4.19 * inputs.atmospheric_temperature_k * (matrix_excess + macro_excess);
    return .{ .litter_water_m3 = inputs.rain_and_irrigation_to_litter_m3 + matrix_excess + macro_excess, .matrix_water_m3 = inputs.rain_and_irrigation_to_matrix_m3 - matrix_excess, .macropore_water_m3 = inputs.rain_and_irrigation_to_macropore_m3 - macro_excess, .litter_heat_megajoules = inputs.litter_input_heat_megajoules + redirected_heat, .soil_heat_megajoules = inputs.soil_input_heat_megajoules - redirected_heat };
}

pub const GasExchangeParameters = struct {
    reference_time_h: f64,
    wet_exponent: f64,
    dry_exponent: f64,
    transition_water_fraction: f64,
    iteration_fraction: f64,
    aqueous_tortuosity_coefficient: f64,
};

pub const GasExchange = struct { air_water_rate_per_step: f64, aqueous_tortuosity: f64 };

pub fn litterGasExchange(total_pore_volume_m3: f64, ice_m3: f64, water_m3: f64, air_m3: f64, water_retention_capacity_m3: f64, parameters: GasExchangeParameters) !GasExchange {
    inline for (.{ total_pore_volume_m3, ice_m3, water_m3, air_m3, water_retention_capacity_m3 }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfacePrecipitationInput;
    inline for (@typeInfo(GasExchangeParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteSurfacePrecipitationInput;
    if (total_pore_volume_m3 < 0 or ice_m3 < 0 or water_m3 < 0 or air_m3 < 0 or water_retention_capacity_m3 < 0 or parameters.reference_time_h <= 0 or parameters.transition_water_fraction < 0 or parameters.transition_water_fraction > 1 or parameters.iteration_fraction < 0 or parameters.iteration_fraction > 1 or parameters.aqueous_tortuosity_coefficient < 0) return error.InvalidSurfacePrecipitationInput;
    const available_pores = total_pore_volume_m3 - ice_m3;
    var rate: f64 = 0;
    if (available_pores > 0 and air_m3 > 0) {
        const relative_water = std.math.clamp(water_m3 / available_pores, 0, 1);
        const exponent = if (relative_water > parameters.transition_water_fraction) parameters.wet_exponent else parameters.dry_exponent;
        const rate_per_h = 1.0 / ((1.0 / parameters.reference_time_h) * @exp(exponent * (relative_water - parameters.transition_water_fraction)));
        // WATSUB applied DFGS*XNPT explicitly inside repeated gas subcycles.
        // The full-hour local solve instead integrates the same first-order
        // kinetic rate exactly. This remains bounded and removes an artificial
        // donor-clamp discontinuity from the Newton/Picard residual.
        rate = -std.math.expm1(-@max(0.0, rate_per_h) * parameters.iteration_fraction);
    }
    const water_holding_fraction = if (water_retention_capacity_m3 > 0) @min(1.0, water_m3 / water_retention_capacity_m3) else 1.0;
    return .{ .air_water_rate_per_step = rate, .aqueous_tortuosity = parameters.aqueous_tortuosity_coefficient * water_holding_fraction * water_holding_fraction };
}

pub const RainfallImpactParameters = struct {
    direct_energy_intercept_j_per_mm: f64,
    direct_energy_log_coefficient_j_per_mm: f64,
    throughfall_energy_height_coefficient_j_per_mm_sqrt_m: f64,
    throughfall_energy_intercept_j_per_mm: f64,
    maximum_canopy_height_m: f64,
    ponding_attenuation_per_mm: f64,
    conductivity_damage_per_j_per_megagram_per_megagram: f64,
    conductivity_recovery_fraction_per_h: f64,
};

pub const RainfallImpactInputs = struct {
    direct_precipitation_mm_per_h: f64,
    throughfall_mm_per_h: f64,
    total_precipitation_mm_per_h: f64,
    canopy_height_m: f64,
    excess_surface_storage_m3: f64,
    ground_surface_retention_m3: f64,
    surface_area_m2: f64,
    bare_soil_fraction: f64,
    time_fraction: f64,
    surface_silt_megagrams_per_megagram: f64,
    surface_clay_megagrams_per_megagram: f64,
};

pub const RainfallImpact = struct { incremental_energy_j: f64, cumulative_energy_j: f64, saturated_conductivity_multiplier: f64 };

/// Exact WATSUB `ENGYP=ENGYP*(1-FENGYP*XNFH)` surface-conductivity
/// recovery. Energy is J m-2, the recovery coefficient is h-1, and the
/// timestep is h. Apply before adding current-timestep rainfall energy,
/// including on dry timesteps.
pub fn recoverRainfallImpactEnergy(previous_cumulative_energy_j_per_m2: f64, timestep_h: f64, recovery_fraction_per_h: f64) !f64 {
    inline for (.{ previous_cumulative_energy_j_per_m2, timestep_h, recovery_fraction_per_h }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePrecipitationInput;
    if (previous_cumulative_energy_j_per_m2 < 0 or timestep_h <= 0 or recovery_fraction_per_h < 0 or recovery_fraction_per_h * timestep_h > 1)
        return error.InvalidSurfacePrecipitationInput;
    const recovered = previous_cumulative_energy_j_per_m2 * (1.0 - recovery_fraction_per_h * timestep_h);
    if (!std.math.isFinite(recovered) or recovered < 0) return error.NonFiniteRainfallImpact;
    return recovered;
}

pub fn rainfallConductivityMultiplier(cumulative_energy_j_per_m2: f64, surface_silt_megagrams_per_megagram: f64, surface_clay_megagrams_per_megagram: f64, damage_per_j_per_megagram_per_megagram: f64) !f64 {
    inline for (.{ cumulative_energy_j_per_m2, surface_silt_megagrams_per_megagram, surface_clay_megagrams_per_megagram, damage_per_j_per_megagram_per_megagram }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePrecipitationInput;
    if (cumulative_energy_j_per_m2 < 0 or surface_silt_megagrams_per_megagram < 0 or surface_clay_megagrams_per_megagram < 0 or damage_per_j_per_megagram_per_megagram < 0)
        return error.InvalidSurfacePrecipitationInput;
    const multiplier = @exp(-damage_per_j_per_megagram_per_megagram * (surface_silt_megagrams_per_megagram + surface_clay_megagrams_per_megagram) * cumulative_energy_j_per_m2);
    if (!std.math.isFinite(multiplier)) return error.NonFiniteRainfallImpact;
    return multiplier;
}

pub fn rainfallImpact(previous_cumulative_energy_j: f64, inputs: RainfallImpactInputs, parameters: RainfallImpactParameters) !RainfallImpact {
    if (!std.math.isFinite(previous_cumulative_energy_j) or previous_cumulative_energy_j < 0) return error.InvalidSurfacePrecipitationInput;
    inline for (@typeInfo(RainfallImpactInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfacePrecipitationInput;
    inline for (@typeInfo(RainfallImpactParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteSurfacePrecipitationInput;
    if (inputs.direct_precipitation_mm_per_h < 0 or inputs.throughfall_mm_per_h < 0 or inputs.total_precipitation_mm_per_h <= 0 or inputs.canopy_height_m < 0 or inputs.excess_surface_storage_m3 < 0 or inputs.ground_surface_retention_m3 < 0 or inputs.surface_area_m2 <= 0 or inputs.bare_soil_fraction < 0 or inputs.bare_soil_fraction > 1 or inputs.time_fraction <= 0 or inputs.time_fraction > 1 or inputs.surface_silt_megagrams_per_megagram < 0 or inputs.surface_clay_megagrams_per_megagram < 0 or parameters.maximum_canopy_height_m < 0 or parameters.ponding_attenuation_per_mm < 0 or parameters.conductivity_damage_per_j_per_megagram_per_megagram < 0 or parameters.conductivity_recovery_fraction_per_h < 0) return error.InvalidSurfacePrecipitationInput;
    const direct_energy = if (inputs.direct_precipitation_mm_per_h > 0) @max(0.0, parameters.direct_energy_intercept_j_per_mm + parameters.direct_energy_log_coefficient_j_per_mm * @log(inputs.total_precipitation_mm_per_h)) else 0;
    const throughfall_energy = if (inputs.throughfall_mm_per_h > 0) @max(0.0, parameters.throughfall_energy_height_coefficient_j_per_mm_sqrt_m * @sqrt(@min(parameters.maximum_canopy_height_m, inputs.canopy_height_m)) - parameters.throughfall_energy_intercept_j_per_mm) else 0;
    const ponded_depth_mm = 1.0e3 * @max(0.0, inputs.excess_surface_storage_m3 - inputs.ground_surface_retention_m3) / inputs.surface_area_m2;
    const incremental = if (direct_energy + throughfall_energy > 0) (direct_energy * inputs.direct_precipitation_mm_per_h + throughfall_energy * inputs.throughfall_mm_per_h) * @exp(-parameters.ponding_attenuation_per_mm * ponded_depth_mm) * inputs.bare_soil_fraction * inputs.time_fraction else 0;
    const cumulative = previous_cumulative_energy_j + incremental;
    const multiplier = try rainfallConductivityMultiplier(cumulative, inputs.surface_silt_megagrams_per_megagram, inputs.surface_clay_megagrams_per_megagram, parameters.conductivity_damage_per_j_per_megagram_per_megagram);
    if (!std.math.isFinite(incremental) or !std.math.isFinite(cumulative) or !std.math.isFinite(multiplier)) return error.NonFiniteRainfallImpact;
    return .{ .incremental_energy_j = incremental, .cumulative_energy_j = cumulative, .saturated_conductivity_multiplier = multiplier };
}

test "the runtime snow cover carrier uses the WATSUB square-root relation" {
    // Regression guard for the exponent defect. `prepareFromModel` is the only
    // writer of `snow_cover_fraction`, and every consumer derives its snow-free
    // fraction as `1 - snow_cover_fraction`, so an inverted exponent here
    // silently mis-partitions rain, heat, and albedo for every thin-snow hour.
    // This test pins the delegation rather than re-deriving the formula.
    const owner = @import("../soil/water/snow_cover_fraction.zig");
    const full_cover_m = 0.07;
    const cover = try owner.evaluate(0.005, full_cover_m);
    // The value production used to produce at this depth.
    const previous_squared = std.math.pow(f64, 0.005 / full_cover_m, 2);
    try std.testing.expect(cover.snow_fraction > previous_squared * 50);
    // And the snow-free complement consumers compute stays a valid fraction.
    const snow_free = 1.0 - cover.snow_fraction;
    try std.testing.expectApproxEqAbs(cover.snow_free_fraction, snow_free, 1e-15);
    try std.testing.expect(snow_free > 0 and snow_free < 1);
}

test "precipitation exceeding pore air is redirected to litter with heat" {
    const result = try redistribute(.{ .soil_surface_present = true, .rain_and_irrigation_to_matrix_m3 = 0.3, .rain_and_irrigation_to_macropore_m3 = 0.2, .rain_and_irrigation_to_litter_m3 = 0.1, .matrix_air_capacity_m3 = 0.1, .macropore_air_capacity_m3 = 0.1, .snow_free_fraction = 1, .atmospheric_temperature_k = 280, .litter_input_heat_megajoules = 1, .soil_input_heat_megajoules = 10 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result.litter_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), result.matrix_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), result.macropore_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 11), result.litter_heat_megajoules + result.soil_heat_megajoules, 1e-12);
}

test "runtime surface ingress retains separate rain snow litter and pore carriers" {
    var state = try RuntimeState.init(std.testing.allocator, 2);
    defer state.deinit();
    state.litter_water_capacity_m3[0] = 0.2;
    state.litter_water_capacity_m3[1] = 0.2;
    state.litter_cover_fraction[0] = 0.5;
    state.litter_cover_fraction[1] = 0.5;
    try prepareRuntimeIngress(&state, .{ .rainfall_m3_per_h = &.{ 1, 0.5 }, .snowfall_water_equivalent_m3_per_h = &.{ 0.2, 0 }, .intercepted_rain_m3_per_h = &.{ 0.1, 0 }, .snow_cover_fraction = &.{ 0.25, 0 }, .atmospheric_temperature_k = &.{ 280, 285 }, .matrix_fraction = &.{ 0.8, 0.7 }, .macropore_fraction = &.{ 0.2, 0.3 }, .matrix_air_capacity_m3 = &.{ 0.1, 1 }, .macropore_air_capacity_m3 = &.{ 0.1, 1 }, .litter_absent_above_water_table = &.{ false, false }, .maximum_iterations = 4 });
    const first_liquid = state.rain_to_snow_m3_per_h[0] + state.water_to_litter_m3_per_h[0] + state.water_to_matrix_m3_per_h[0] + state.water_to_macropore_m3_per_h[0];
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), first_liquid, 1.0e-12);
    try commitRuntimeIngress(&state, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2) + state.rain_to_snow_m3_per_h[0], state.solid_snow_water_equivalent_m3[0], 1.0e-12);
    try std.testing.expect(state.litter_water_m3[0] > 0);
}

test "top-soil ingress commit validates every runtime cell before mutation" {
    const config = try @import("../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 4 });
    var grid = try @import("../state/grid.zig").GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var hydrology = try @import("../transport/hydrology.zig").State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    var state = try RuntimeState.init(std.testing.allocator, 2);
    defer state.deinit();
    @memset(grid.matrix_pore_capacity_m3, 1);
    @memset(grid.macropore_pore_capacity_m3, 1);
    @memset(grid.matrix_air_volume_m3, 1);
    @memset(grid.macropore_air_volume_m3, 1);
    @memset(grid.air_volume_m3, 2);
    state.water_to_matrix_m3_per_h[0] = 0.2;
    state.water_to_matrix_m3_per_h[1] = 2;
    try std.testing.expectError(error.InvalidSurfaceIngressCommit, commitSoilIngress(&state, &grid, &hydrology, 1));
    try std.testing.expectEqual(@as(f64, 0), grid.matrix_liquid_water_m3[0]);
    state.water_to_matrix_m3_per_h[1] = 0.3;
    state.water_to_macropore_m3_per_h[1] = 0.4;
    try commitSoilIngress(&state, &grid, &hydrology, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), grid.matrix_liquid_water_m3[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), grid.liquid_water_m3[1], 1e-15);
    try std.testing.expectApproxEqAbs(grid.matrix_liquid_water_m3[1], hydrology.micropore_water_volume_m3[1], 1e-15);
}

test "atmospheric rain snow litter and pore partition conserves liquid water and heat" {
    const result = try partitionAtmosphericWater(.{ .rain_and_irrigation_m3_per_h = 1, .intercepted_rain_m3_per_h = 0.1, .snowfall_water_equivalent_m3_per_h = 0.2, .snow_cover_fraction = 0.25, .snow_free_fraction = 0.75, .litter_cover_fraction = 0.5, .litter_water_capacity_m3 = 1, .litter_water_m3 = 0.9, .litter_absent_above_water_table = false, .atmospheric_temperature_k = 280, .matrix_fraction = 0.8, .macropore_fraction = 0.2, .maximum_iterations = 4 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), result.rain_to_snow_m3_per_h + result.water_to_litter_m3_per_h + result.water_to_matrix_m3_per_h + result.water_to_macropore_m3_per_h, 1e-12);
    try std.testing.expectApproxEqAbs(4.19 * 280 * 0.9 + 2.095 * 280 * 0.2, result.heat_to_snow_megajoules_per_h + result.heat_to_litter_megajoules_per_h + result.heat_to_soil_megajoules_per_h, 1e-10);
}

test "solute precipitation routing follows snow presence and dry-surface branches" {
    const snow = try routePrecipitationSolutes(0.1, 0.5, 0.6, 0.4, 0.2, 0, 0, 0.5);
    try std.testing.expectEqual(@as(f64, 0), snow.rain_to_litter_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), snow.rain_to_snow_m3, 1e-12);
    const surface = try routePrecipitationSolutes(0, 0.6, 0.6, 0.4, 0.2, 0, 0, 0.5);
    try std.testing.expectApproxEqAbs(@as(f64, 0.06), surface.rain_to_litter_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.04), surface.irrigation_to_litter_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), surface.rain_to_snow_m3 + surface.irrigation_to_snow_m3 + surface.rain_to_litter_m3 + surface.irrigation_to_litter_m3, 1e-12);
}

test "rainfall impact attenuates conductivity and accumulates after recovery" {
    const recovered = try recoverRainfallImpactEnergy(10, 1, 5.0e-4);
    const result = try rainfallImpact(recovered, .{ .direct_precipitation_mm_per_h = 2, .throughfall_mm_per_h = 1, .total_precipitation_mm_per_h = 3, .canopy_height_m = 2, .excess_surface_storage_m3 = 0, .ground_surface_retention_m3 = 0, .surface_area_m2 = 1, .bare_soil_fraction = 1, .time_fraction = 1, .surface_silt_megagrams_per_megagram = 0.2, .surface_clay_megagrams_per_megagram = 0.2 }, .{ .direct_energy_intercept_j_per_mm = 8.95, .direct_energy_log_coefficient_j_per_mm = 8.44, .throughfall_energy_height_coefficient_j_per_mm_sqrt_m = 15.8, .throughfall_energy_intercept_j_per_mm = 5.87, .maximum_canopy_height_m = 2.5, .ponding_attenuation_per_mm = 2, .conductivity_damage_per_j_per_megagram_per_megagram = 1e-3, .conductivity_recovery_fraction_per_h = 5.0e-4 });
    try std.testing.expectApproxEqAbs(@as(f64, 9.995), recovered, 2.0e-15);
    try std.testing.expect(result.incremental_energy_j > 0);
    try std.testing.expect(result.cumulative_energy_j > recovered);
    try std.testing.expect(result.saturated_conductivity_multiplier > 0 and result.saturated_conductivity_multiplier < 1);
}

test "WATSUB rainfall impact recovery advances on dry hours and rejects overshoot" {
    try std.testing.expectApproxEqAbs(@as(f64, 37.48125), try recoverRainfallImpactEnergy(37.5, 1, 5.0e-4), 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 37.490625), try recoverRainfallImpactEnergy(37.5, 0.5, 5.0e-4), 1.0e-14);
    try std.testing.expectError(error.InvalidSurfacePrecipitationInput, recoverRainfallImpactEnergy(1, 2, 0.6));
    try std.testing.expectError(error.NonFiniteSurfacePrecipitationInput, recoverRainfallImpactEnergy(std.math.nan(f64), 1, 5.0e-4));
}

test "litter gas exchange response is bounded on dry branch" {
    const result = try litterGasExchange(1, 0, 0.1, 0.9, 0.5, .{ .reference_time_h = 1, .wet_exponent = 2, .dry_exponent = -2, .transition_water_fraction = 0.5, .iteration_fraction = 1, .aqueous_tortuosity_coefficient = 0.7 });
    try std.testing.expect(result.air_water_rate_per_step >= 0 and result.air_water_rate_per_step <= 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.028), result.aqueous_tortuosity, 1e-12);
}

test "whole-step gas exchange integrates legacy substep kinetics exponentially" {
    const result = try litterGasExchange(1, 0, 0.8, 0.2, 1, .{
        .reference_time_h = 2,
        .wet_exponent = 0,
        .dry_exponent = 0,
        .transition_water_fraction = 0.5,
        .iteration_fraction = 1,
        .aqueous_tortuosity_coefficient = 0.7,
    });
    try std.testing.expectApproxEqAbs(-std.math.expm1(@as(f64, -2.0)), result.air_water_rate_per_step, 1e-15);
    try std.testing.expect(result.air_water_rate_per_step < 1);
}

test "litter-soil flux uses litter donor temperature when flow is into soil" {
    const litter_parameters = try retention.carselParrishDefault(.silt_loam, 0.8);
    const soil_parameters = try retention.carselParrishDefault(.loam, null);
    const result = try surface_water_flow.litterSoilFlux(.{
        .litter_water_m3 = 0.7,
        .soil_matrix_water_m3 = 0.2,
        .litter_air_m3 = 0.1,
        .soil_matrix_air_m3 = 0.5,
        .litter_volume_m3 = 1,
        .soil_matrix_bulk_volume_m3 = 1,
        // Water content alone does not order two pools by potential: the
        // litter's theta_s is 0.8, so it must be near-saturated on its own
        // curve to sit above a soil at theta = 0.2.
        .litter_water_fraction = 0.7,
        .soil_water_fraction = 0.2,
        .litter_parameters = litter_parameters,
        .soil_parameters = soil_parameters,
        .litter_external_water_potential_megapascal = 0,
        .soil_external_water_potential_megapascal = 0,
        .litter_thickness_m = 0.05,
        .soil_thickness_m = 0.05,
        .soil_face_area_m2 = 1,
        .litter_cover_fraction = 1,
        .wet_litter_cover_fraction = 1,
        .time_fraction = 1,
        .soil_excess_pore_volume_m3 = 0,
        .litter_temperature_k = 300,
        .soil_temperature_k = 280,
    });
    // The wetter litter is at the higher (less negative) matric potential on
    // its own curve, so water must move down into the drier soil.
    try std.testing.expect(result.water_m3 > 0);
    try std.testing.expectApproxEqAbs(
        4.19 * 300 * result.water_m3,
        result.convective_heat_megajoules,
        1e-12,
    );
}

test "litter-soil flux uses soil donor temperature when flow is upward" {
    const litter_parameters = try retention.carselParrishDefault(.silt_loam, 0.8);
    const soil_parameters = try retention.carselParrishDefault(.silt_loam, 0.8);
    const inputs = surface_water_flow.LitterSoilInputs{
        .litter_water_m3 = 0.2,
        .soil_matrix_water_m3 = 1,
        .litter_air_m3 = 0.5,
        .soil_matrix_air_m3 = 0.5,
        .litter_volume_m3 = 1,
        .soil_matrix_bulk_volume_m3 = 1,
        .litter_water_fraction = 0.2,
        .soil_water_fraction = 0.4,
        .litter_parameters = litter_parameters,
        .soil_parameters = soil_parameters,
        .litter_external_water_potential_megapascal = 0,
        .soil_external_water_potential_megapascal = 0,
        .litter_thickness_m = 0.05,
        .soil_thickness_m = 0.05,
        .soil_face_area_m2 = 1,
        .litter_cover_fraction = 1,
        .wet_litter_cover_fraction = 1,
        .time_fraction = 1,
        .soil_excess_pore_volume_m3 = 0,
        .litter_temperature_k = 300,
        .soil_temperature_k = 280,
    };
    const result = try surface_water_flow.litterSoilFlux(inputs);
    try std.testing.expect(result.water_m3 <= 0);
    try std.testing.expectApproxEqAbs(
        4.19 * 280 * result.water_m3,
        result.convective_heat_megajoules,
        1e-15,
    );
}
