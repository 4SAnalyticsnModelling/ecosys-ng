const std = @import("std");
const grid_input_files = @import("../io/input/grid_input_files.zig");
const runscript = @import("../driver/runscript.zig");
const site_catalog = @import("site_catalog.zig");
const topography_catalog = @import("topography_catalog.zig");
const topography_module = @import("topography.zig");

pub const Assignments = struct {
    allocator: std.mem.Allocator,
    column_count: usize,
    row_count: usize,
    site_catalog_index_by_cell: []usize,
    topography_catalog_index_by_cell: []usize,
    topography_unit_index_by_cell: []usize,
    horizontal_cell_width_m: []f64,
    vertical_cell_width_m: []f64,
    initial_water_table_depth_m: []f64,
    natural_water_table_surface_slope: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        domain: runscript.Domain,
        files: grid_input_files.CellFiles,
        sites: site_catalog.Catalog,
        topographies: topography_catalog.Catalog,
    ) !Assignments {
        const column_count = try domain.columns();
        const row_count = try domain.rows();
        const cell_count = try std.math.mul(usize, column_count, row_count);
        if (files.column_count != column_count or files.row_count != row_count or
            files.site_file_by_cell.len != cell_count or
            files.topography_file_by_cell.len != cell_count or
            files.soil_file_by_cell.len != cell_count)
            return error.GridEnvironmentDimensionMismatch;

        var result: Assignments = undefined;
        result.allocator = allocator;
        result.column_count = column_count;
        result.row_count = row_count;
        var allocated: usize = 0;
        errdefer result.freeAllocated(allocated);
        inline for (.{
            "site_catalog_index_by_cell",
            "topography_catalog_index_by_cell",
            "topography_unit_index_by_cell",
        }) |field_name| {
            @field(result, field_name) = try allocator.alloc(usize, cell_count);
            allocated += 1;
        }
        result.horizontal_cell_width_m = try allocator.alloc(f64, cell_count);
        allocated += 1;
        result.vertical_cell_width_m = try allocator.alloc(f64, cell_count);
        allocated += 1;
        result.initial_water_table_depth_m = try allocator.alloc(f64, cell_count);
        allocated += 1;
        result.natural_water_table_surface_slope = try allocator.alloc(f64, cell_count);
        allocated += 1;

        for (0..cell_count) |cell| {
            const local_column = cell % column_count;
            const local_row = cell / column_count;
            const global_column = try std.math.add(usize, domain.west_column, local_column);
            const global_row = try std.math.add(usize, domain.north_row, local_row);

            const site_index = sites.find(files.site_file_by_cell[cell]) orelse
                return error.SiteFileNotLoaded;
            const topography_index = topographies.find(files.topography_file_by_cell[cell]) orelse
                return error.TopographyFileNotLoaded;
            const selected_site = sites.entries.items[site_index].site;
            if (selected_site.horizontal_cell_widths_m.len != 1 or
                selected_site.vertical_cell_widths_m.len != 1)
                return error.SiteDoesNotCoverGridCell;
            const unit_index = try topographies.entries.items[topography_index].topography.unitForCell(
                global_column,
                global_row,
            );

            result.site_catalog_index_by_cell[cell] = site_index;
            result.topography_catalog_index_by_cell[cell] = topography_index;
            result.topography_unit_index_by_cell[cell] = unit_index;
            result.horizontal_cell_width_m[cell] = selected_site.horizontal_cell_widths_m[0];
            result.vertical_cell_width_m[cell] = selected_site.vertical_cell_widths_m[0];
            result.initial_water_table_depth_m[cell] = selected_site.initial_water_table_depth_m;
            result.natural_water_table_surface_slope[cell] = selected_site.natural_water_table_surface_slope;
        }
        return result;
    }

    pub fn deinit(self: *Assignments) void {
        self.freeAllocated(7);
        self.* = undefined;
    }

    /// Materializes one landscape unit per runtime cell for consumers that
    /// operate on a domain-wide topography. Terrain attributes come from the
    /// selected per-cell topography file; the soil filename comes exclusively
    /// from the compulsory grid mapping.
    pub fn buildDomainTopography(
        self: Assignments,
        domain: runscript.Domain,
        files: grid_input_files.CellFiles,
        topographies: topography_catalog.Catalog,
    ) !topography_module.Topography {
        const cell_count = try std.math.mul(usize, self.column_count, self.row_count);
        if (files.soil_file_by_cell.len != cell_count) return error.GridEnvironmentDimensionMismatch;
        const units = try self.allocator.alloc(topography_module.LandscapeUnit, cell_count);
        errdefer self.allocator.free(units);
        var initialized: usize = 0;
        errdefer for (units[0..initialized]) |unit| self.allocator.free(unit.soil_profile_file);

        for (0..cell_count) |cell| {
            const local_column = cell % self.column_count;
            const local_row = cell / self.column_count;
            const global_column = try std.math.add(usize, domain.west_column, local_column);
            const global_row = try std.math.add(usize, domain.north_row, local_row);
            const catalog_index = self.topography_catalog_index_by_cell[cell];
            if (catalog_index >= topographies.entries.items.len)
                return error.TopographyCatalogIndexOutOfRange;
            const selected_topography = topographies.entries.items[catalog_index].topography;
            const unit_index = self.topography_unit_index_by_cell[cell];
            if (unit_index >= selected_topography.units.len)
                return error.TopographyUnitIndexOutOfRange;
            const source_unit = selected_topography.units[unit_index];
            const soil_file = try self.allocator.dupe(u8, files.soil_file_by_cell[cell]);
            errdefer self.allocator.free(soil_file);
            units[cell] = .{
                .west_column = global_column,
                .north_row = global_row,
                .east_column = global_column,
                .south_row = global_row,
                .compass_aspect_degrees = source_unit.compass_aspect_degrees,
                .geometric_aspect_degrees = source_unit.geometric_aspect_degrees,
                .slope_degrees = source_unit.slope_degrees,
                .unused_slope_input = source_unit.unused_slope_input,
                .initial_snowpack_depth_m = source_unit.initial_snowpack_depth_m,
                .soil_profile_file = soil_file,
            };
            initialized += 1;
        }
        return .{ .allocator = self.allocator, .units = units };
    }

    fn freeAllocated(self: *Assignments, allocated: usize) void {
        if (allocated >= 7) self.allocator.free(self.natural_water_table_surface_slope);
        if (allocated >= 6) self.allocator.free(self.initial_water_table_depth_m);
        if (allocated >= 5) self.allocator.free(self.vertical_cell_width_m);
        if (allocated >= 4) self.allocator.free(self.horizontal_cell_width_m);
        if (allocated >= 3) self.allocator.free(self.topography_unit_index_by_cell);
        if (allocated >= 2) self.allocator.free(self.topography_catalog_index_by_cell);
        if (allocated >= 1) self.allocator.free(self.site_catalog_index_by_cell);
    }
};

test "runtime cells resolve repeated and distinct cached environments" {
    const allocator = std.testing.allocator;
    const site_source = @import("../core/test_fixtures.zig").site_source;
    var sites = site_catalog.Catalog.init(allocator);
    defer sites.deinit();
    _ = try sites.appendFromSource("site_shared", site_source, 1, 1);
    _ = try sites.appendFromSource("site_east", site_source, 1, 1);

    var topographies = topography_catalog.Catalog.init(allocator);
    defer topographies.deinit();
    _ = try topographies.appendFromSource(
        "terrain_shared",
        "1 1 2 1 90 1 0 0\nunused_soil\n",
    );
    _ = try topographies.appendFromSource(
        "terrain_east",
        "2 1 2 1 180 2 0 0\nunused_soil\n",
    );

    var files = try grid_input_files.CellFiles.parse(
        allocator,
        \\grid_cell 1 1 site_shared terrain_shared soil_west
        \\grid_cell 2 1 site_east terrain_east soil_east
    ,
        2,
        1,
    );
    defer files.deinit();
    var assignments = try Assignments.init(
        allocator,
        .{ .west_column = 1, .north_row = 1, .east_column = 2, .south_row = 1 },
        files,
        sites,
        topographies,
    );
    defer assignments.deinit();

    try std.testing.expect(assignments.site_catalog_index_by_cell[0] != assignments.site_catalog_index_by_cell[1]);
    try std.testing.expect(assignments.topography_catalog_index_by_cell[0] != assignments.topography_catalog_index_by_cell[1]);
    try std.testing.expectEqual(@as(usize, 0), assignments.topography_unit_index_by_cell[0]);
    try std.testing.expectEqual(@as(usize, 0), assignments.topography_unit_index_by_cell[1]);
}

test "environment assignment requires every referenced cache entry" {
    const allocator = std.testing.allocator;
    var sites = site_catalog.Catalog.init(allocator);
    defer sites.deinit();
    var topographies = topography_catalog.Catalog.init(allocator);
    defer topographies.deinit();
    var files = try grid_input_files.CellFiles.parse(
        allocator,
        "grid_cell 1 1 absent_site absent_topography soil\n",
        1,
        1,
    );
    defer files.deinit();
    try std.testing.expectError(
        error.SiteFileNotLoaded,
        Assignments.init(
            allocator,
            .{ .west_column = 1, .north_row = 1, .east_column = 1, .south_row = 1 },
            files,
            sites,
            topographies,
        ),
    );
}

test "domain topography uses selected terrain but explicit mapped soil" {
    const allocator = std.testing.allocator;
    const site_source = @import("../core/test_fixtures.zig").site_source;
    var sites = site_catalog.Catalog.init(allocator);
    defer sites.deinit();
    _ = try sites.appendFromSource("site", site_source, 1, 1);
    var topographies = topography_catalog.Catalog.init(allocator);
    defer topographies.deinit();
    _ = try topographies.appendFromSource(
        "terrain",
        "1 1 1 1 135 7 0 0.2\nignored_topography_soil\n",
    );
    var files = try grid_input_files.CellFiles.parse(
        allocator,
        "grid_cell 1 1 site terrain authoritative_soil\n",
        1,
        1,
    );
    defer files.deinit();
    const domain: runscript.Domain = .{
        .west_column = 1,
        .north_row = 1,
        .east_column = 1,
        .south_row = 1,
    };
    var assignments = try Assignments.init(allocator, domain, files, sites, topographies);
    defer assignments.deinit();
    var combined = try assignments.buildDomainTopography(domain, files, topographies);
    defer combined.deinit();
    try std.testing.expectEqualStrings("authoritative_soil", combined.units[0].soil_profile_file);
    try std.testing.expectEqual(@as(f64, 135), combined.units[0].compass_aspect_degrees);
    try std.testing.expectEqual(@as(f64, 7), combined.units[0].slope_degrees);
}
