const std = @import("std");
const delimited_input = @import("delimited_input.zig");

pub const CellFiles = struct {
    allocator: std.mem.Allocator,
    column_count: usize,
    row_count: usize,
    site_file_by_cell: [][]u8,
    topography_file_by_cell: [][]u8,
    soil_file_by_cell: [][]u8,

    pub fn parse(
        allocator: std.mem.Allocator,
        source: []const u8,
        column_count: usize,
        row_count: usize,
    ) !CellFiles {
        const cell_count = try validateDimensions(column_count, row_count);
        const site_files = try allocator.alloc([]u8, cell_count);
        errdefer allocator.free(site_files);
        const topography_files = try allocator.alloc([]u8, cell_count);
        errdefer allocator.free(topography_files);
        const soil_files = try allocator.alloc([]u8, cell_count);
        errdefer allocator.free(soil_files);
        const assigned = try allocator.alloc(bool, cell_count);
        defer allocator.free(assigned);
        @memset(assigned, false);

        errdefer for (assigned, 0..) |is_assigned, cell| {
            if (is_assigned) {
                allocator.free(site_files[cell]);
                allocator.free(topography_files[cell]);
                allocator.free(soil_files[cell]);
            }
        };

        var records = delimited_input.records(source);
        while (records.next()) |record| {
            if (hasEmptyExplicitField(record))
                return error.EmptyGridInputRecordValue;
            var fields = delimited_input.recordTokens(record);
            const record_name = fields.next() orelse unreachable;
            if (!std.ascii.eqlIgnoreCase(record_name, "grid_cell"))
                return error.InvalidGridCellInputRecord;
            const column = try coordinate(&fields);
            const row = try coordinate(&fields);
            if (column > column_count or row > row_count)
                return error.GridCellInputCoordinateOutOfRange;
            const cell = (row - 1) * column_count + column - 1;
            if (assigned[cell]) return error.DuplicateGridCellCoordinate;
            const site_file = try duplicateFileName(allocator, &fields);
            errdefer allocator.free(site_file);
            const topography_file = try duplicateFileName(allocator, &fields);
            errdefer allocator.free(topography_file);
            const soil_file = try duplicateFileName(allocator, &fields);
            errdefer allocator.free(soil_file);
            try requireEnd(&fields);
            site_files[cell] = site_file;
            topography_files[cell] = topography_file;
            soil_files[cell] = soil_file;
            assigned[cell] = true;
        }
        for (assigned) |is_assigned| if (!is_assigned)
            return error.MissingGridCellInput;

        return .{
            .allocator = allocator,
            .column_count = column_count,
            .row_count = row_count,
            .site_file_by_cell = site_files,
            .topography_file_by_cell = topography_files,
            .soil_file_by_cell = soil_files,
        };
    }

    pub fn deinit(self: *CellFiles) void {
        for (0..self.site_file_by_cell.len) |cell| {
            self.allocator.free(self.site_file_by_cell[cell]);
            self.allocator.free(self.topography_file_by_cell[cell]);
            self.allocator.free(self.soil_file_by_cell[cell]);
        }
        self.allocator.free(self.site_file_by_cell);
        self.allocator.free(self.topography_file_by_cell);
        self.allocator.free(self.soil_file_by_cell);
        self.* = undefined;
    }
};

pub const WeatherFiles = struct {
    allocator: std.mem.Allocator,
    column_count: usize,
    row_count: usize,
    file_by_cell: [][]u8,

    pub fn parse(
        allocator: std.mem.Allocator,
        source: []const u8,
        column_count: usize,
        row_count: usize,
    ) !WeatherFiles {
        const cell_count = try validateDimensions(column_count, row_count);
        const files = try allocator.alloc([]u8, cell_count);
        errdefer allocator.free(files);
        const assigned = try allocator.alloc(bool, cell_count);
        defer allocator.free(assigned);
        @memset(assigned, false);

        errdefer for (files, assigned) |file, is_assigned| {
            if (is_assigned) allocator.free(file);
        };

        var records = delimited_input.records(source);
        while (records.next()) |record| {
            if (hasEmptyExplicitField(record))
                return error.EmptyGridInputRecordValue;
            var fields = delimited_input.recordTokens(record);
            const record_name = fields.next() orelse unreachable;
            if (!std.ascii.eqlIgnoreCase(record_name, "weather_cell"))
                return error.InvalidWeatherCellInputRecord;
            const column = try coordinate(&fields);
            const row = try coordinate(&fields);
            if (column > column_count or row > row_count)
                return error.WeatherCellInputCoordinateOutOfRange;
            const cell = (row - 1) * column_count + column - 1;
            if (assigned[cell]) return error.DuplicateWeatherCellCoordinate;
            const file = try duplicateFileName(allocator, &fields);
            errdefer allocator.free(file);
            try requireEnd(&fields);
            files[cell] = file;
            assigned[cell] = true;
        }
        for (assigned) |is_assigned| if (!is_assigned)
            return error.MissingWeatherCellInput;
        return .{
            .allocator = allocator,
            .column_count = column_count,
            .row_count = row_count,
            .file_by_cell = files,
        };
    }

    pub fn deinit(self: *WeatherFiles) void {
        for (self.file_by_cell) |file| self.allocator.free(file);
        self.allocator.free(self.file_by_cell);
        self.* = undefined;
    }
};

fn validateDimensions(column_count: usize, row_count: usize) !usize {
    if (column_count == 0 or row_count == 0) return error.InvalidGridInputDimensions;
    return std.math.mul(usize, column_count, row_count);
}

fn coordinate(fields: *delimited_input.TokenIterator) !usize {
    const text = fields.next() orelse return error.IncompleteGridInputRecord;
    const value = std.fmt.parseUnsigned(usize, text, 10) catch
        return error.InvalidGridInputCoordinate;
    if (value == 0) return error.GridInputCoordinateIsZero;
    return value;
}

fn duplicateFileName(
    allocator: std.mem.Allocator,
    fields: *delimited_input.TokenIterator,
) ![]u8 {
    const name = fields.next() orelse return error.IncompleteGridInputRecord;
    if (name.len == 0 or std.mem.indexOfAny(u8, name, "\r\n#") != null)
        return error.InvalidGridInputFileName;
    return allocator.dupe(u8, name);
}

fn requireEnd(fields: *delimited_input.TokenIterator) !void {
    if (fields.next() != null) return error.TrailingGridInputRecordData;
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

test "runtime grid cells may use distinct site topography soil and weather files" {
    var cell_files = try CellFiles.parse(
        std.testing.allocator,
        \\# Records may repeat filenames or select distinct inputs.
        \\grid_cell,1,1,site_a,topography_a,soil_a
        \\GrId_CeLl|2|1|site_b|topography_b|soil_b # eastern cell
    ,
        2,
        1,
    );
    defer cell_files.deinit();
    try std.testing.expectEqualStrings("site_a", cell_files.site_file_by_cell[0]);
    try std.testing.expectEqualStrings("site_b", cell_files.site_file_by_cell[1]);
    try std.testing.expectEqualStrings("topography_b", cell_files.topography_file_by_cell[1]);
    try std.testing.expectEqualStrings("soil_a", cell_files.soil_file_by_cell[0]);

    var weather = try WeatherFiles.parse(
        std.testing.allocator,
        \\weather_cell 1 1 weather_shared.csv
        \\WEATHER_CELL 2 1 weather_east
    ,
        2,
        1,
    );
    defer weather.deinit();
    try std.testing.expectEqualStrings("weather_shared.csv", weather.file_by_cell[0]);
    try std.testing.expectEqualStrings("weather_east", weather.file_by_cell[1]);
}

test "adjacent cells may intentionally reuse every input filename" {
    var cell_files = try CellFiles.parse(
        std.testing.allocator,
        \\grid_cell 1 1 shared_site shared_topography shared_soil
        \\grid_cell 2 1 shared_site shared_topography shared_soil
    ,
        2,
        1,
    );
    defer cell_files.deinit();
    try std.testing.expectEqualStrings(cell_files.site_file_by_cell[0], cell_files.site_file_by_cell[1]);
    try std.testing.expectEqualStrings(cell_files.topography_file_by_cell[0], cell_files.topography_file_by_cell[1]);
    try std.testing.expectEqualStrings(cell_files.soil_file_by_cell[0], cell_files.soil_file_by_cell[1]);

    var weather = try WeatherFiles.parse(
        std.testing.allocator,
        \\weather_cell 1 1 shared_weather
        \\weather_cell 2 1 shared_weather
    ,
        2,
        1,
    );
    defer weather.deinit();
    try std.testing.expectEqualStrings(weather.file_by_cell[0], weather.file_by_cell[1]);
}

test "per-cell input mappings reject missing repeated coordinates and trailing values" {
    try std.testing.expectError(
        error.MissingGridCellInput,
        CellFiles.parse(
            std.testing.allocator,
            "grid_cell 1 1 site topography soil\n",
            2,
            1,
        ),
    );
    try std.testing.expectError(
        error.DuplicateWeatherCellCoordinate,
        WeatherFiles.parse(
            std.testing.allocator,
            "weather_cell 1 1 a\nweather_cell 1 1 b\n",
            1,
            1,
        ),
    );
    try std.testing.expectError(
        error.TrailingGridInputRecordData,
        CellFiles.parse(
            std.testing.allocator,
            "grid_cell 1 1 site topography soil extra\n",
            1,
            1,
        ),
    );
}

test "per-cell input mappings reject empty explicit delimiter fields" {
    inline for (.{
        "grid_cell,1,1,site,,soil\n",
        "grid_cell|1|1|site| |soil\n",
        "grid_cell\t1\t1\tsite\t\tsoil\n",
    }) |source| try std.testing.expectError(
        error.EmptyGridInputRecordValue,
        CellFiles.parse(std.testing.allocator, source, 1, 1),
    );

    try std.testing.expectError(
        error.EmptyGridInputRecordValue,
        WeatherFiles.parse(
            std.testing.allocator,
            "weather_cell,1,1, # missing compulsory filename\n",
            1,
            1,
        ),
    );
}

test "grid input empty-field check preserves spacing and comments" {
    try std.testing.expect(!hasEmptyExplicitField(
        "grid_cell  1  1  site  topography  soil # valid spaces",
    ));
    try std.testing.expect(!hasEmptyExplicitField(
        "grid_cell, 1 | 1\tsite, topography | soil # mixed delimiters",
    ));
    try std.testing.expect(!hasEmptyExplicitField("# comment only"));
}
