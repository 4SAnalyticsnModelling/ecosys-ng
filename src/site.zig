const std = @import("std");
const delimited_input = @import("delimited_input.zig");

pub const LateralConnectionMode = enum(u8) {
    connected = 1,
    disconnected = 3,
};

pub const Site = struct {
    allocator: std.mem.Allocator,
    latitude_degrees_north: f64,
    longitude_degrees_east: f64,
    elevation_m: f64,
    mean_annual_air_temperature_c: f64,
    water_table_mode: u8,
    atmospheric_oxygen_umol_mol: f64,
    atmospheric_nitrogen_umol_mol: f64,
    atmospheric_co2_umol_mol: f64,
    atmospheric_methane_umol_mol: f64,
    atmospheric_nitrous_oxide_umol_mol: f64,
    atmospheric_ammonia_umol_mol: f64,
    ecosystem_type: i32,
    salinity_enabled: bool,
    erosion_mode: i32,
    lateral_connection_mode: LateralConnectionMode,
    initial_water_table_depth_m: f64,
    natural_water_table_surface_slope: f64,
    artificial_water_table_depth_m: ?f64,
    artificial_water_table_surface_slope: ?f64,
    /// North, east, south, west, matching READI's record order.
    surface_runoff_boundary_fraction: [4]f64,
    natural_water_table_distance_m: [4]f64,
    natural_subsurface_exchange_fraction: [4]f64,
    lower_boundary_exchange_fraction: f64,
    artificial_water_table_distance_m: [4]f64,
    artificial_subsurface_exchange_fraction: [4]f64,
    horizontal_cell_widths_m: []f64,
    vertical_cell_widths_m: []f64,

    pub fn deinit(self: *Site) void {
        self.allocator.free(self.horizontal_cell_widths_m);
        self.allocator.free(self.vertical_cell_widths_m);
        self.* = undefined;
    }

    /// Legacy IERSNG modes 1 and 3 include erosion. Mode -1 disables all
    /// profile disturbance, 0 is freeze-thaw only, and 2 is freeze-thaw plus
    /// soil-organic-matter gain/loss.
    pub fn erosionEnabled(self: Site) bool {
        return self.erosion_mode == 1 or self.erosion_mode == 3;
    }
};

/// Parses strict, line-aware runtime site records.
pub fn parse(allocator: std.mem.Allocator, source: []const u8, east_column: usize, south_row: usize) !Site {
    if (east_column == 0 or south_row == 0) return error.InvalidSiteDimensions;
    var records = delimited_input.records(source);
    var record1 = delimited_input.recordTokens(try nextRecord(&records));
    var record2 = delimited_input.recordTokens(try nextRecord(&records));
    var record3 = delimited_input.recordTokens(try nextRecord(&records));
    var record4 = delimited_input.recordTokens(try nextRecord(&records));

    var result: Site = undefined;
    result.allocator = allocator;
    result.latitude_degrees_north = try number(f64, &record1);
    result.longitude_degrees_east = try number(f64, &record1);
    result.elevation_m = try number(f64, &record1);
    result.mean_annual_air_temperature_c = try number(f64, &record1);
    result.water_table_mode = try number(u8, &record1);
    result.atmospheric_oxygen_umol_mol = try number(f64, &record2);
    result.atmospheric_nitrogen_umol_mol = try number(f64, &record2);
    result.atmospheric_co2_umol_mol = try number(f64, &record2);
    result.atmospheric_methane_umol_mol = try number(f64, &record2);
    result.atmospheric_nitrous_oxide_umol_mol = try number(f64, &record2);
    result.atmospheric_ammonia_umol_mol = try number(f64, &record2);
    result.ecosystem_type = try number(i32, &record3);
    result.salinity_enabled = (try number(i32, &record3)) != 0;
    result.erosion_mode = try number(i32, &record3);
    result.lateral_connection_mode = switch (try number(u8, &record3)) {
        1 => .connected,
        3 => .disconnected,
        else => return error.InvalidLateralConnectionMode,
    };
    result.initial_water_table_depth_m = try number(f64, &record3);
    result.natural_water_table_surface_slope = try number(f64, &record3);
    for (&result.surface_runoff_boundary_fraction) |*value| value.* = try number(f64, &record4);
    for (&result.natural_water_table_distance_m) |*value| value.* = try number(f64, &record4);
    for (&result.natural_subsurface_exchange_fraction) |*value| value.* = try number(f64, &record4);
    result.lower_boundary_exchange_fraction = try number(f64, &record4);
    try requireEnd(&record1);
    try requireEnd(&record2);
    try requireEnd(&record3);
    try requireEnd(&record4);
    if (result.water_table_mode >= 3) {
        var artificial_table = delimited_input.recordTokens(try nextRecord(&records));
        result.artificial_water_table_depth_m = try number(f64, &artificial_table);
        result.artificial_water_table_surface_slope = try number(f64, &artificial_table);
        try requireEnd(&artificial_table);
        var artificial_boundaries = delimited_input.recordTokens(try nextRecord(&records));
        for (&result.artificial_water_table_distance_m) |*value| value.* = try number(f64, &artificial_boundaries);
        for (&result.artificial_subsurface_exchange_fraction) |*value| value.* = try number(f64, &artificial_boundaries);
        try requireEnd(&artificial_boundaries);
    } else {
        result.artificial_water_table_depth_m = null;
        result.artificial_water_table_surface_slope = null;
        result.artificial_water_table_distance_m = [_]f64{0} ** 4;
        result.artificial_subsurface_exchange_fraction = [_]f64{0} ** 4;
    }

    var horizontal_widths = delimited_input.recordTokens(try nextRecord(&records));
    var vertical_widths = delimited_input.recordTokens(try nextRecord(&records));

    result.horizontal_cell_widths_m = try allocator.alloc(f64, east_column);
    errdefer allocator.free(result.horizontal_cell_widths_m);
    for (result.horizontal_cell_widths_m) |*width| width.* = try number(f64, &horizontal_widths);
    result.vertical_cell_widths_m = try allocator.alloc(f64, south_row);
    errdefer allocator.free(result.vertical_cell_widths_m);
    for (result.vertical_cell_widths_m) |*width| width.* = try number(f64, &vertical_widths);
    try requireEnd(&horizontal_widths);
    try requireEnd(&vertical_widths);
    if (records.next() != null) return error.TrailingSiteRecord;
    try validate(result);
    return result;
}

fn validate(site: Site) !void {
    if (!std.math.isFinite(site.latitude_degrees_north) or site.latitude_degrees_north < -90 or site.latitude_degrees_north > 90) return error.InvalidLatitude;
    if (!std.math.isFinite(site.longitude_degrees_east) or site.longitude_degrees_east < -180 or site.longitude_degrees_east > 180) return error.InvalidLongitude;
    if (!std.math.isFinite(site.elevation_m)) return error.InvalidElevation;
    if (!std.math.isFinite(site.mean_annual_air_temperature_c)) return error.InvalidMeanAnnualTemperature;
    if (site.water_table_mode > 4) return error.InvalidWaterTableMode;
    if (site.erosion_mode < -1 or site.erosion_mode > 3) return error.InvalidErosionMode;
    if (!std.math.isFinite(site.natural_water_table_surface_slope) or site.natural_water_table_surface_slope < 0 or site.natural_water_table_surface_slope > 1) return error.InvalidWaterTableSlope;
    if (site.artificial_water_table_surface_slope) |slope| if (!std.math.isFinite(slope) or slope < 0 or slope > 1) return error.InvalidWaterTableSlope;
    for (site.surface_runoff_boundary_fraction) |value| if (!unitFraction(value)) return error.InvalidSurfaceRunoffBoundaryFraction;
    for (site.natural_water_table_distance_m) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidNaturalWaterTableDistance;
    for (site.natural_subsurface_exchange_fraction) |value| if (!unitFraction(value)) return error.InvalidNaturalSubsurfaceExchangeFraction;
    if (!unitFraction(site.lower_boundary_exchange_fraction)) return error.InvalidLowerBoundaryExchangeFraction;
    for (site.artificial_water_table_distance_m) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidArtificialWaterTableDistance;
    for (site.artificial_subsurface_exchange_fraction) |value| if (!unitFraction(value)) return error.InvalidArtificialSubsurfaceExchangeFraction;
    for (site.horizontal_cell_widths_m) |width| if (!std.math.isFinite(width) or width <= 0) return error.InvalidCellWidth;
    for (site.vertical_cell_widths_m) |width| if (!std.math.isFinite(width) or width <= 0) return error.InvalidCellWidth;
}

fn unitFraction(value: f64) bool {
    return std.math.isFinite(value) and value >= 0 and value <= 1;
}

fn nextRecord(records: anytype) ![]const u8 {
    const record = records.next() orelse return error.UnexpectedEndOfSiteFile;
    if (hasEmptyExplicitField(record)) return error.EmptySiteRecordValue;
    return record;
}

fn hasEmptyExplicitField(record: []const u8) bool {
    const content = if (std.mem.indexOfScalar(u8, record, '#')) |comment|
        record[0..comment]
    else
        record;
    const trimmed = std.mem.trim(u8, content, " \r");
    if (trimmed.len == 0) return false;

    var field_start: usize = 0;
    var saw_explicit_delimiter = false;
    for (trimmed, 0..) |byte, index| {
        if (byte != ',' and byte != '|' and byte != '\t') continue;
        if (std.mem.trim(u8, trimmed[field_start..index], " \r").len == 0)
            return true;
        field_start = index + 1;
        saw_explicit_delimiter = true;
    }
    return saw_explicit_delimiter and
        std.mem.trim(u8, trimmed[field_start..], " \r").len == 0;
}

fn requireEnd(tokens: anytype) !void {
    if (tokens.next() != null) return error.TrailingSiteRecordData;
}

fn number(comptime T: type, tokens: anytype) !T {
    const text = tokens.next() orelse return error.UnexpectedEndOfSiteRecord;
    return switch (@typeInfo(T)) {
        .float => std.fmt.parseFloat(T, text),
        .int => std.fmt.parseInt(T, text, 10),
        else => @compileError("unsupported site number type"),
    };
}

const test_site_source = "45.3 -75.7 92 5.4 3\n2.1E+05 7.8E+05 360.0 1.8 0.3 0.002\n33 1 3 1 1.0 0.0\n0 1 1 0 10 0 10 0 1 0 1 0 0\n1.5 1.0\n0 10 0 10 0 1 0 1\n1\n1\n";

test "parse self-contained site with Fortran record semantics" {
    const source = test_site_source;
    var site = try parse(std.testing.allocator, source, 1, 1);
    defer site.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 45.3), site.latitude_degrees_north, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -75.7), site.longitude_degrees_east, 1.0e-12);
    try std.testing.expectEqual(@as(i32, 33), site.ecosystem_type);
    try std.testing.expect(site.salinity_enabled);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), site.horizontal_cell_widths_m[0], 1.0e-12);
    try std.testing.expectEqual([_]f64{ 0, 1, 1, 0 }, site.surface_runoff_boundary_fraction);
    try std.testing.expectEqual(LateralConnectionMode.connected, site.lateral_connection_mode);
    try std.testing.expectEqual([_]f64{ 10, 0, 10, 0 }, site.natural_water_table_distance_m);
    try std.testing.expectEqual([_]f64{ 1, 0, 1, 0 }, site.natural_subsurface_exchange_fraction);
    try std.testing.expectEqual(@as(f64, 0), site.lower_boundary_exchange_fraction);
    try std.testing.expectEqual([_]f64{ 0, 10, 0, 10 }, site.artificial_water_table_distance_m);
    try std.testing.expectEqual([_]f64{ 0, 1, 0, 1 }, site.artificial_subsurface_exchange_fraction);
    try std.testing.expect(site.erosionEnabled());
}

test "erosion mode follows the complete IERSNG option domain" {
    inline for ([_]i32{ -1, 0, 1, 2, 3 }) |mode| {
        var source_buffer: [512]u8 = undefined;
        const source = try std.fmt.bufPrint(
            &source_buffer,
            "45.3 -75.7 92 5.4 0\n2.1E+05 7.8E+05 360.0 1.8 0.3 0.002\n33 1 {d} 1 1.0 0.0\n0 1 1 0 10 0 10 0 1 0 1 0 0\n1\n1\n",
            .{mode},
        );
        var site = try parse(std.testing.allocator, source, 1, 1);
        defer site.deinit();
        try std.testing.expectEqual(mode == 1 or mode == 3, site.erosionEnabled());
    }
}

test "erosion mode outside IERSNG domain fails immediately" {
    const source = "45.3 -75.7 92 5.4 0\n2.1E+05 7.8E+05 360.0 1.8 0.3 0.002\n33 1 4 1 1.0 0.0\n0 1 1 0 10 0 10 0 1 0 1 0 0\n1\n1\n";
    try std.testing.expectError(error.InvalidErosionMode, parse(std.testing.allocator, source, 1, 1));
}

test "site records reject empty explicit delimiter fields" {
    inline for (.{
        "45.3,-75.7,,92,5.4,0\n",
        "45.3|-75.7| |92|5.4|0\n",
        "45.3\t-75.7\t\t92\t5.4\t0\n",
    }) |first_record| {
        const source = first_record ++
            "2.1E+05 7.8E+05 360.0 1.8 0.3 0.002\n" ++
            "33 1 3 1 1.0 0.0\n" ++
            "0 1 1 0 10 0 10 0 1 0 1 0 0\n" ++
            "1.5 1.0\n0 10 0 10 0 1 0 1\n1\n1\n";
        try std.testing.expectError(
            error.EmptySiteRecordValue,
            parse(std.testing.allocator, source, 1, 1),
        );
    }
}

test "site empty-field check preserves valid spacing and comments" {
    try std.testing.expect(!hasEmptyExplicitField(
        "45.3  -75.7  92  5.4  0 # valid spaces",
    ));
    try std.testing.expect(!hasEmptyExplicitField(
        "45.3, -75.7 | 92\t5.4, 0 # mixed delimiters",
    ));
    try std.testing.expect(!hasEmptyExplicitField("# comment only"));
}
