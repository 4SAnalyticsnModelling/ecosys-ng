const std = @import("std");

pub const Time = struct {
    day_of_year: usize, // I
    hour: usize, // J
    year: i32, // IYRC
    forcing_step: usize, // NFZ
    final_forcing_step: usize, // NFH
};

pub const Cell = struct {
    column: usize, // NX
    row: usize, // NY
};

pub const Inputs = struct {
    time: Time,
    cell: Cell,
    last_layer: usize, // NL; slices include surface index zero
    organic_carbon_g_c: []const f64, // ORGC
    charcoal_organic_carbon_g_c: []const f64, // ORGCC
    horizontal_area_m2: []const f64, // AREA(3,...)
};

fn finite(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

pub fn isAnnualDiagnosticTime(time: Time) bool {
    return time.day_of_year == 365 and time.hour == 24 and time.forcing_step == time.final_forcing_step;
}

fn writeRecord(writer: *std.Io.Writer, label: []const u8, time: Time, cell: Cell, values_g_m2: []const f64) !void {
    try writer.print("{s}\t{d}\t{d}\t{d}\t{d}\t{d}", .{ label, time.day_of_year, time.hour, time.year, cell.column, cell.row });
    for (values_g_m2) |value| try writer.print("\t{e}", .{value});
    try writer.writeByte('\n');
}

/// Direct translation of REDIST 12848--12854. The legacy fixed-width records
/// are intentionally emitted as mandatory tab-delimited text records.
pub fn writeAnnual(allocator: std.mem.Allocator, organic_writer: *std.Io.Writer, charcoal_writer: *std.Io.Writer, inputs: Inputs) !void {
    if (!isAnnualDiagnosticTime(inputs.time)) return;
    const count = inputs.last_layer + 1;
    if (inputs.organic_carbon_g_c.len != count or inputs.charcoal_organic_carbon_g_c.len != count or inputs.horizontal_area_m2.len != count) return error.AnnualOrganicDiagnosticDimensionMismatch;
    if (!finite(inputs.organic_carbon_g_c) or !finite(inputs.charcoal_organic_carbon_g_c) or !finite(inputs.horizontal_area_m2)) return error.InvalidAnnualOrganicDiagnosticInput;
    const staged = try allocator.alloc(f64, 2 * count);
    defer allocator.free(staged);
    for (0..count) |layer| {
        const area_m2 = inputs.horizontal_area_m2[layer];
        if (area_m2 <= 0) return error.InvalidAnnualOrganicDiagnosticArea;
        staged[layer] = inputs.organic_carbon_g_c[layer] / area_m2;
        staged[count + layer] = inputs.charcoal_organic_carbon_g_c[layer] / area_m2;
        if (!std.math.isFinite(staged[layer]) or !std.math.isFinite(staged[count + layer])) return error.NonFiniteAnnualOrganicDiagnosticResult;
    }
    try writeRecord(organic_writer, "ORGC", inputs.time, inputs.cell, staged[0..count]);
    try writeRecord(charcoal_writer, "ORGCC", inputs.time, inputs.cell, staged[count .. 2 * count]);
}

test "REDIST annual organic diagnostics are runtime-sized and tab-delimited" {
    var organic_bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer organic_bytes.deinit();
    var charcoal_bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer charcoal_bytes.deinit();
    const organic = [_]f64{ 10, 40, 90 };
    const charcoal = [_]f64{ 2, 8, 18 };
    const area = [_]f64{ 2, 4, 10 };
    try writeAnnual(std.testing.allocator, &organic_bytes.writer, &charcoal_bytes.writer, .{
        .time = .{ .day_of_year = 365, .hour = 24, .year = 2001, .forcing_step = 4, .final_forcing_step = 4 },
        .cell = .{ .column = 7, .row = 9 },
        .last_layer = 2,
        .organic_carbon_g_c = &organic,
        .charcoal_organic_carbon_g_c = &charcoal,
        .horizontal_area_m2 = &area,
    });
    try std.testing.expectEqualStrings("ORGC\t365\t24\t2001\t7\t9\t5e0\t1e1\t9e0\n", organic_bytes.written());
    try std.testing.expectEqualStrings("ORGCC\t365\t24\t2001\t7\t9\t1e0\t2e0\t1.8e0\n", charcoal_bytes.written());
}

test "REDIST annual diagnostic gate skips invalid inactive records" {
    var organic_bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer organic_bytes.deinit();
    var charcoal_bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer charcoal_bytes.deinit();
    try writeAnnual(std.testing.allocator, &organic_bytes.writer, &charcoal_bytes.writer, .{
        .time = .{ .day_of_year = 364, .hour = 24, .year = 2001, .forcing_step = 1, .final_forcing_step = 2 },
        .cell = .{ .column = 0, .row = 0 },
        .last_layer = 99,
        .organic_carbon_g_c = &.{},
        .charcoal_organic_carbon_g_c = &.{},
        .horizontal_area_m2 = &.{},
    });
    try std.testing.expectEqual(@as(usize, 0), organic_bytes.written().len);
    try std.testing.expectEqual(@as(usize, 0), charcoal_bytes.written().len);
}
