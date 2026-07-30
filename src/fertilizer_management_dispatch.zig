const std = @import("std");
const delimited_input = @import("delimited_input.zig");
const fertilizer_schedule = @import("fertilizer_schedule.zig");
const land_management = @import("land_management.zig");
const Date = @import("options.zig").Date;
const nitrogen_inventory = @import("fertilizer_nitrogen_inventory.zig");
const surface_fertilizer = @import("surface_litter_fertilizer.zig");
const organic = @import("soil_organic_initialization.zig");
const organic_parameters = @import("soil_organic_parameters.zig");
const organic_application = @import("organic_fertilizer_application.zig");
const mineral_fertilizer = @import("mineral_fertilizer_inventory.zig");
const execution_calendar_date = @import("execution_calendar_date.zig");

/// Dense runtime lookup built once when a scene is activated. `null` denotes
/// an explicit case-insensitive NO schedule, never a missing file lookup.
pub const ScheduleMap = struct {
    allocator: std.mem.Allocator,
    catalog_index_by_cell: []?usize,

    pub fn init(allocator: std.mem.Allocator, assignments: land_management.Assignments, unit_by_cell: []const usize, catalog: fertilizer_schedule.Catalog) !ScheduleMap {
        const map = try allocator.alloc(?usize, unit_by_cell.len);
        errdefer allocator.free(map);
        for (unit_by_cell, 0..) |unit_index, cell| {
            if (unit_index >= assignments.units.len) return error.LandManagementUnitIndexOutOfBounds;
            const name = assignments.units[unit_index].fertilizer_file;
            map[cell] = if (delimited_input.isNo(name)) null else catalog.find(name) orelse return error.FertilizerScheduleMissingFromCatalog;
        }
        return .{ .allocator = allocator, .catalog_index_by_cell = map };
    }

    pub fn deinit(self: *ScheduleMap) void {
        self.allocator.free(self.catalog_index_by_cell);
        self.* = undefined;
    }
};

pub fn dispatchDate(map: ScheduleMap, catalog: fertilizer_schedule.Catalog, date: Date, context: anytype, comptime apply: fn (@TypeOf(context), usize, *const fertilizer_schedule.Event) anyerror!void) !usize {
    try validateDispatchDate(date);
    var applied: usize = 0;
    for (map.catalog_index_by_cell, 0..) |maybe_schedule, cell| {
        const schedule_index = maybe_schedule orelse continue;
        if (schedule_index >= catalog.entries.items.len) return error.FertilizerScheduleIndexOutOfBounds;
        for (catalog.entries.items[schedule_index].events) |*event| {
            if (event.date.day != date.day or event.date.month != date.month or (!event.date.isRecurring() and event.date.year != date.year)) continue;
            try apply(context, cell, event);
            applied = try std.math.add(usize, applied, 1);
        }
    }
    return applied;
}

fn validateDispatchDate(date: Date) !void {
    if (date.year == 0) return error.InvalidFertilizerDispatchDate;
    _ = execution_calendar_date.dayOfYear(.{ .day = date.day, .month = date.month, .year = date.year }) catch return error.InvalidFertilizerDispatchDate;
}

pub const NitrogenApplyContext = struct {
    soil: *nitrogen_inventory.State,
    surface: *surface_fertilizer.State,
    cell_area_m2: []const f64,
    active_soil_layer_count: []const usize,
    soil_layer_thickness_m: []const f64,
    nitrogen_molar_mass_g_per_mol: f64,
    surface_organic: *const organic.State,
    source_hour_one_through_twenty_four: u8,
    solar_noon_hour_by_cell: []const u8,
};

pub fn applyNitrogen(context: *NitrogenApplyContext, cell: usize, event: *const fertilizer_schedule.Event) !void {
    if (!try isApplicationHour(context.source_hour_one_through_twenty_four, context.solar_noon_hour_by_cell, cell)) return;
    if (cell >= context.soil.cell_count or cell >= context.cell_area_m2.len or cell >= context.active_soil_layer_count.len) return error.FertilizerDispatchCellOutOfBounds;
    const layer_count = context.active_soil_layer_count[cell];
    if (layer_count == 0 or layer_count > context.soil.layer_capacity) return error.InvalidActiveSoilLayerCount;
    const first = try std.math.mul(usize, cell, context.soil.layer_capacity);
    if (first + layer_count > context.soil_layer_thickness_m.len) return error.FertilizerDispatchLayerExtentMismatch;
    const area_m2 = context.cell_area_m2[cell];
    if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidFertilizerCellArea;
    const surface_carbon_g_c = try context.surface_organic.totalCarbon_g_c(cell);
    const cover_fraction = 1.0 - @exp(-0.008 * surface_carbon_g_c / area_m2);
    try nitrogen_inventory.applyEventNitrogen(context.soil, context.surface, cell, area_m2, context.nitrogen_molar_mass_g_per_mol, cover_fraction, context.soil_layer_thickness_m[first .. first + layer_count], event.*);
}

pub const MineralApplyContext = struct {
    inventory: *mineral_fertilizer.State,
    surface_organic: *const organic.State,
    cell_area_m2: []const f64,
    active_soil_layer_count: []const usize,
    soil_layer_thickness_m: []const f64,
    source_hour_one_through_twenty_four: u8,
    solar_noon_hour_by_cell: []const u8,
};

pub fn applyMinerals(context: *MineralApplyContext, cell: usize, event: *const fertilizer_schedule.Event) !void {
    if (!try isApplicationHour(context.source_hour_one_through_twenty_four, context.solar_noon_hour_by_cell, cell)) return;
    if (cell >= context.inventory.cell_count or cell >= context.surface_organic.layer_count or cell >= context.cell_area_m2.len or cell >= context.active_soil_layer_count.len) return error.FertilizerDispatchCellOutOfBounds;
    const active_layers = context.active_soil_layer_count[cell];
    const first = try std.math.mul(usize, cell, context.inventory.layer_capacity);
    if (active_layers == 0 or active_layers > context.inventory.layer_capacity or first + active_layers > context.soil_layer_thickness_m.len) return error.FertilizerDispatchLayerExtentMismatch;
    const area_m2 = context.cell_area_m2[cell];
    if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidFertilizerCellArea;
    const surface_carbon_g_c = try context.surface_organic.totalCarbon_g_c(cell);
    const cover_fraction = 1.0 - @exp(-0.008 * surface_carbon_g_c / area_m2);
    try mineral_fertilizer.applyEvent(context.inventory, cell, area_m2, cover_fraction, context.soil_layer_thickness_m[first .. first + active_layers], event.*);
}

pub const OrganicApplyContext = struct {
    soil: *organic.State,
    surface: *organic.State,
    parameters: *const organic_parameters.OwnedParameters,
    cell_area_m2: []const f64,
    active_soil_layer_count: []const usize,
    soil_layer_capacity: usize,
    soil_layer_thickness_m: []const f64,
    daily_organic_carbon_input_g_c: []f64,
    daily_biome_carbon_input_g_c: []f64,
    daily_organic_phosphorus_input_g_p: []f64,
    daily_organic_nitrogen_input_g_n: []f64,
    source_hour_one_through_twenty_four: u8,
    solar_noon_hour_by_cell: []const u8,
};

pub fn applyOrganic(context: *OrganicApplyContext, cell: usize, event: *const fertilizer_schedule.Event) !void {
    if (!try isApplicationHour(context.source_hour_one_through_twenty_four, context.solar_noon_hour_by_cell, cell)) return;
    if (cell >= context.cell_area_m2.len or cell >= context.active_soil_layer_count.len or cell >= context.daily_organic_carbon_input_g_c.len or cell >= context.daily_biome_carbon_input_g_c.len or cell >= context.daily_organic_phosphorus_input_g_p.len or cell >= context.daily_organic_nitrogen_input_g_n.len or cell >= context.surface.layer_count) return error.FertilizerDispatchCellOutOfBounds;
    const active_layers = context.active_soil_layer_count[cell];
    const first = try std.math.mul(usize, cell, context.soil_layer_capacity);
    if (active_layers == 0 or active_layers > context.soil_layer_capacity or first + active_layers > context.soil_layer_thickness_m.len or first + active_layers > context.soil.layer_count) return error.FertilizerDispatchLayerExtentMismatch;
    const area_m2 = context.cell_area_m2[cell];
    if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidFertilizerCellArea;
    const local_layer = try organicLayerAtDepth(context.soil_layer_thickness_m[first .. first + active_layers], event.application_depth_m);
    const destination = if (event.application_depth_m == 0) context.surface else context.soil;
    const destination_layer = if (event.application_depth_m == 0) cell else first + local_layer;
    const plant: organic.ElementPool = .{
        .carbon_g_c = event.plant_residue_g_per_m2.carbon * area_m2,
        .nitrogen_g_n = event.plant_residue_g_per_m2.nitrogen * area_m2,
        .phosphorus_g_p = event.plant_residue_g_per_m2.phosphorus * area_m2,
    };
    const manure: organic.ElementPool = .{
        .carbon_g_c = event.manure_g_per_m2.carbon * area_m2,
        .nitrogen_g_n = event.manure_g_per_m2.nitrogen * area_m2,
        .phosphorus_g_p = event.manure_g_per_m2.phosphorus * area_m2,
    };
    try organic_application.apply(destination, destination_layer, .plant_residue, event.plant_residue_type, plant, context.parameters);
    try organic_application.apply(destination, destination_layer, .manure, event.manure_type, manure, context.parameters);
    const next = context.daily_organic_carbon_input_g_c[cell] + plant.carbon_g_c + manure.carbon_g_c;
    if (!std.math.isFinite(next)) return error.OrganicFertilizerApplicationOverflow;
    context.daily_organic_carbon_input_g_c[cell] = next;
    const phosphorus_next = context.daily_organic_phosphorus_input_g_p[cell] + plant.phosphorus_g_p + manure.phosphorus_g_p;
    if (!std.math.isFinite(phosphorus_next)) return error.OrganicFertilizerApplicationOverflow;
    context.daily_organic_phosphorus_input_g_p[cell] = phosphorus_next;
    const nitrogen_next = context.daily_organic_nitrogen_input_g_n[cell] + plant.nitrogen_g_n + manure.nitrogen_g_n;
    if (!std.math.isFinite(nitrogen_next)) return error.OrganicFertilizerApplicationOverflow;
    context.daily_organic_nitrogen_input_g_n[cell] = nitrogen_next;
    if (event.manure_type < 3) {
        const biome_next = context.daily_biome_carbon_input_g_c[cell] + plant.carbon_g_c + manure.carbon_g_c;
        if (!std.math.isFinite(biome_next)) return error.OrganicFertilizerApplicationOverflow;
        context.daily_biome_carbon_input_g_c[cell] = biome_next;
    }
}

fn isApplicationHour(source_hour: u8, solar_noon_by_cell: []const u8, cell: usize) !bool {
    if (source_hour == 0 or source_hour > 24) return error.InvalidHourlyFertilizerSchedule;
    if (cell >= solar_noon_by_cell.len) return error.FertilizerDispatchCellOutOfBounds;
    const solar_noon = solar_noon_by_cell[cell];
    if (solar_noon > 24) return error.InvalidHourlyFertilizerSchedule;
    return source_hour == solar_noon;
}

fn organicLayerAtDepth(thickness_m: []const f64, depth_m: f64) !usize {
    var lower_boundary_m: f64 = 0;
    for (thickness_m, 0..) |thickness, layer| {
        lower_boundary_m += thickness;
        if (depth_m <= lower_boundary_m) return layer;
    }
    return error.FertilizerApplicationBelowSoilProfile;
}

test "fertilizer dispatch handles recurring dates and case-insensitive NO" {
    const allocator = std.testing.allocator;
    var assignments = try land_management.parse(allocator, "1 1 1 1\ntillage annual NO\n2 1 2 1\ntillage no irrigation\n");
    defer assignments.deinit();
    const units = try assignments.buildCellUnitMap(allocator, .{ .west_column = 1, .north_row = 1, .east_column = 2, .south_row = 1 });
    defer allocator.free(units);
    var catalog = fertilizer_schedule.Catalog.init(allocator);
    defer catalog.deinit();
    _ = try catalog.appendFromSource("annual", "01050000 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n");
    var map = try ScheduleMap.init(allocator, assignments, units, catalog);
    defer map.deinit();
    var count: usize = 0;
    const Callback = struct {
        fn apply(output: *usize, cell: usize, event: *const fertilizer_schedule.Event) !void {
            try std.testing.expectEqual(@as(usize, 0), cell);
            try std.testing.expectEqual(@as(f64, 1), event.nitrogen_g_per_m2.broadcast_ammonium);
            output.* += 1;
        }
    };
    try std.testing.expectEqual(@as(usize, 1), try dispatchDate(map, catalog, .{ .day = 1, .month = 5, .year = 2026 }, &count, Callback.apply));
    try std.testing.expectEqual(@as(usize, 1), count);
}

test "fertilizer dispatch date preserves DAY modulo-four chronology" {
    try validateDispatchDate(.{ .day = 29, .month = 2, .year = 1900 });
    try validateDispatchDate(.{ .day = 30, .month = 4, .year = 1901 });
    try std.testing.expectError(error.InvalidFertilizerDispatchDate, validateDispatchDate(.{ .day = 0, .month = 2, .year = 1901 }));
    try std.testing.expectError(
        error.InvalidFertilizerDispatchDate,
        validateDispatchDate(.{ .day = 1, .month = 1, .year = 0 }),
    );
}

test "HOUR1 fertilizer admission uses each grid cell solar noon" {
    const solar_noon = [_]u8{ 11, 13 };
    try std.testing.expect(try isApplicationHour(11, &solar_noon, 0));
    try std.testing.expect(!try isApplicationHour(11, &solar_noon, 1));
    try std.testing.expect(!try isApplicationHour(12, &solar_noon, 0));
    try std.testing.expect(try isApplicationHour(13, &solar_noon, 1));
    try std.testing.expectError(
        error.InvalidHourlyFertilizerSchedule,
        isApplicationHour(0, &solar_noon, 0),
    );
    try std.testing.expectError(
        error.FertilizerDispatchCellOutOfBounds,
        isApplicationHour(11, &solar_noon, 2),
    );
}

test "organic dispatch resolves runtime depth and publishes UORGF and eligible TNBP input" {
    const allocator = std.testing.allocator;
    var soil = try organic.State.init(allocator, 2);
    defer soil.deinit();
    var surface = try organic.State.init(allocator, 1);
    defer surface.deinit();
    var parameters = try organic_parameters.sourceParameters(allocator);
    defer parameters.deinit();
    var daily_organic = [_]f64{0};
    var daily_biome = [_]f64{0};
    var daily_phosphorus = [_]f64{0};
    var daily_nitrogen = [_]f64{0};
    var context: OrganicApplyContext = .{
        .soil = &soil,
        .surface = &surface,
        .parameters = &parameters,
        .cell_area_m2 = &.{10},
        .active_soil_layer_count = &.{2},
        .soil_layer_capacity = 2,
        .soil_layer_thickness_m = &.{ 0.1, 0.2 },
        .daily_organic_carbon_input_g_c = &daily_organic,
        .daily_biome_carbon_input_g_c = &daily_biome,
        .daily_organic_phosphorus_input_g_p = &daily_phosphorus,
        .daily_organic_nitrogen_input_g_n = &daily_nitrogen,
        .source_hour_one_through_twenty_four = 12,
        .solar_noon_hour_by_cell = &.{12},
    };
    const event: fertilizer_schedule.Event = .{
        .date = .{ .day = 1, .month = 5, .year = 0 },
        .nitrogen_g_per_m2 = .{ .broadcast_ammonium = 0, .broadcast_ammonia = 0, .broadcast_urea = 0, .broadcast_nitrate = 0, .banded_ammonium = 0, .banded_ammonia = 0, .banded_urea = 0, .banded_nitrate = 0 },
        .phosphorus_g_per_m2 = .{ .broadcast_monocalcium_phosphate = 0, .banded_monocalcium_phosphate = 0, .broadcast_hydroxyapatite = 0 },
        .calcium_carbonate_g_ca_per_m2 = 0,
        .calcium_sulfate_g_ca_per_m2 = 0,
        .plant_residue_g_per_m2 = .{ .carbon = 2, .nitrogen = 0.2, .phosphorus = 0.02 },
        .manure_g_per_m2 = .{ .carbon = 3, .nitrogen = 0.3, .phosphorus = 0.03 },
        .application_depth_m = 0.15,
        .band_row_width_m = 0,
        .fertilizer_formulation = 0,
        .plant_residue_type = 2,
        .manure_type = 2,
    };
    try applyOrganic(&context, 0, &event);
    try std.testing.expectEqual(@as(f64, 50), daily_organic[0]);
    try std.testing.expectEqual(@as(f64, 50), daily_biome[0]);
    try std.testing.expectEqual(@as(f64, 0.5), daily_phosphorus[0]);
    try std.testing.expectEqual(@as(f64, 5), daily_nitrogen[0]);
    var upper_dissolved_carbon_g_c: f64 = 0;
    var lower_dissolved_carbon_g_c: f64 = 0;
    for (soil.dissolved[0..organic.substrate_count]) |pool| upper_dissolved_carbon_g_c += pool.carbon_g_c;
    for (soil.dissolved[organic.substrate_count .. 2 * organic.substrate_count]) |pool| lower_dissolved_carbon_g_c += pool.carbon_g_c;
    try std.testing.expectEqual(@as(f64, 0), upper_dissolved_carbon_g_c);
    try std.testing.expect(lower_dissolved_carbon_g_c > 0);
}
