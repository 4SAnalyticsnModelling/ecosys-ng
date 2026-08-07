//! Input-path resolution, output catalogs and calendar helpers.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
/// Resolves the species input file name used to label a plant output file, so
/// output names carry a readable species instead of the source model's positional
/// digit.
///
/// The runscript declares a population capacity (five in the Ottawa example) while a
/// scene may assign fewer species, so indices beyond the assigned set are real and
/// expected. Those populations carry no species and would emit all-zero rows, so
/// this returns `null` for them and the caller skips the write entirely: an output
/// file is created only for `soil_or_eco` and for genuinely assigned species. The
/// same applies to a run with no plant assignments at all, which then writes no
/// plant files rather than a set of empty ones.
pub fn outputSpeciesLabel(
    assignments: ?ecosys.plant_assignment.Assignments,
    unit_by_cell: ?[]const usize,
    cell: usize,
    species: usize,
) ?[]const u8 {
    if (assignments) |resolved| if (unit_by_cell) |units| {
        if (cell < units.len) {
            const unit_index = units[cell];
            if (unit_index < resolved.units.len) {
                const assigned = resolved.units[unit_index].species;
                if (species < assigned.len) return assigned[species].species_file;
            }
        }
    };
    return null;
}

pub fn sameCalendarDay(left: ecosys.weather.Timestamp, right: ecosys.weather.Timestamp) bool {
    return left.year == right.year and left.day_of_year == right.day_of_year and left.month == right.month and left.day_of_month == right.day_of_month;
}

pub fn sameWeatherTimestamp(left: ecosys.weather.Timestamp, right: ecosys.weather.Timestamp) bool {
    return left.year == right.year and
        left.day_of_year == right.day_of_year and
        left.month == right.month and
        left.day_of_month == right.day_of_month and
        left.hour == right.hour and
        left.minute == right.minute;
}

pub fn dayOfYearFromTimestamp(timestamp: ecosys.weather.Timestamp) !u16 {
    if (timestamp.day_of_year) |day| return day;
    const date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
    return try (ecosys.plant_management.PackedDate{ .day = date.day, .month = date.month, .year = date.year }).dayOfYear(date.year);
}

pub fn soilOutputCatalog(allocator: std.mem.Allocator, editor_index: usize, layers: usize) !ecosys.soil_output_catalog.Catalog {
    return switch (editor_index) {
        0 => ecosys.soil_output_catalog.carbon(allocator, layers, layers, layers),
        1 => ecosys.soil_output_catalog.water(allocator, layers),
        2 => ecosys.soil_output_catalog.nitrogen(allocator, layers, layers),
        3 => ecosys.soil_output_catalog.phosphorus(allocator),
        4 => ecosys.soil_output_catalog.heat(allocator, layers),
        5 => ecosys.soil_output_catalog.dailyCarbon(allocator, layers),
        6 => ecosys.soil_output_catalog.dailyWater(allocator, layers, layers, layers),
        7 => ecosys.soil_output_catalog.dailyNitrogen(allocator, layers),
        8 => ecosys.soil_output_catalog.dailyPhosphorus(allocator, layers),
        9 => ecosys.soil_output_catalog.dailyHeat(allocator, layers, layers),
        else => error.OutputEditorIndexOutOfBounds,
    };
}

pub fn dailySoilWaterPotentialLayerCount(soil_layers: usize) usize {
    return soil_layers;
}

pub fn plantOutputCatalog(allocator: std.mem.Allocator, editor_index: usize, layers: usize) !ecosys.plant_output_catalog.Catalog {
    return switch (editor_index) {
        0 => ecosys.plant_output_catalog.carbon(allocator),
        1 => ecosys.plant_output_catalog.water(allocator, layers),
        2 => ecosys.plant_output_catalog.nitrogen(allocator, layers),
        3 => ecosys.plant_output_catalog.phosphorus(allocator, layers),
        4 => ecosys.plant_output_catalog.heat(allocator),
        5 => ecosys.plant_output_catalog.dailyCarbon(allocator, layers),
        6 => ecosys.plant_output_catalog.dailyWater(allocator),
        7 => ecosys.plant_output_catalog.dailyNitrogen(allocator),
        8 => ecosys.plant_output_catalog.dailyPhosphorus(allocator),
        9 => ecosys.plant_output_catalog.dailyDevelopment(allocator),
        else => error.OutputEditorIndexOutOfBounds,
    };
}

pub fn resolveInputPath(allocator: std.mem.Allocator, io: std.Io, runscript_directory: []const u8, name: []const u8) ![]u8 {
    const direct = try std.fs.path.join(allocator, &.{ runscript_directory, name });
    const direct_file = std.Io.Dir.cwd().openFile(io, direct, .{}) catch |err| switch (err) {
        error.FileNotFound => null,
        else => return err,
    };
    if (direct_file) |file| {
        file.close(io);
        return direct;
    }
    defer allocator.free(direct);
    const parent = std.fs.path.dirname(runscript_directory) orelse return error.InputFileNotFound;
    const inherited = try std.fs.path.join(allocator, &.{ parent, name });
    errdefer allocator.free(inherited);
    const inherited_file = try std.Io.Dir.cwd().openFile(io, inherited, .{});
    inherited_file.close(io);
    return inherited;
}

/// Opens the OS temporary directory, returning null if the path cannot be
/// determined. On Windows, uses GetTempPathW (kernel32, no libc required).
/// On POSIX, reads $TMPDIR or falls back to /tmp. The caller owns the returned
/// Dir handle and must call close(io) on it.
pub fn openOsTempDir(allocator: std.mem.Allocator, io: std.Io) !?std.Io.Dir {
    if (comptime @import("builtin").os.tag == .windows) {
        const GetTempPathW = struct {
            extern "kernel32" fn GetTempPathW(
                nBufferLength: u32,
                lpBuffer: [*]u16,
            ) callconv(.winapi) u32;
        }.GetTempPathW;
        var buf: [std.os.windows.MAX_PATH + 1]u16 = undefined;
        const len = GetTempPathW(buf.len, &buf);
        if (len == 0) return null;
        // GetTempPathW appends a trailing backslash; strip it for openDirAbsolute.
        const raw = buf[0..len];
        const path_u16 = if (raw[raw.len - 1] == '\\') raw[0 .. raw.len - 1] else raw;
        const path = try std.unicode.utf16LeToUtf8Alloc(allocator, path_u16);
        defer allocator.free(path);
        return try std.Io.Dir.openDirAbsolute(io, path, .{});
    } else {
        const env_path = std.posix.getenv("TMPDIR") orelse std.posix.getenv("TMP") orelse "/tmp";
        return try std.Io.Dir.openDirAbsolute(io, env_path, .{});
    }
}
