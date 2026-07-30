const std = @import("std");

pub const Deviation = struct {
    water_m: f64,
    heat_mj_per_m2: f64,
    oxygen_g_per_m2: f64,
    carbon_g_per_m2: f64,
    nitrogen_g_per_m2: f64,
    phosphorus_g_per_m2: f64,
    ions_mol_per_m2: f64,
};

pub const Record = struct {
    execution_day: u64,
    year: i32,
    deviation: Deviation,
};

/// ecosys-ng tab-delimited form of EXEC format 213. Scientific field order is
/// unchanged: water, heat, O2, C, N, P, ions.
pub fn writeHeader(writer: *std.Io.Writer) !void {
    try writer.writeAll(
        "execution_day\tyear\twater_balance_m\t" ++
            "heat_balance_MJ_m-2\toxygen_balance_g_m-2\t" ++
            "carbon_balance_g_m-2\tnitrogen_balance_g_m-2\t" ++
            "phosphorus_balance_g_m-2\tion_balance_mol_m-2\n",
    );
}

/// All values are validated before the first byte is offered to the writer,
/// preventing a NaN or invalid date from leaving a partial scientific row.
pub fn writeRecord(writer: *std.Io.Writer, record: Record) !void {
    if (record.execution_day == 0)
        return error.InvalidExecutionBalanceOutputDay;
    const values = [_]f64{
        record.deviation.water_m,
        record.deviation.heat_mj_per_m2,
        record.deviation.oxygen_g_per_m2,
        record.deviation.carbon_g_per_m2,
        record.deviation.nitrogen_g_per_m2,
        record.deviation.phosphorus_g_per_m2,
        record.deviation.ions_mol_per_m2,
    };
    for (values) |value| if (!std.math.isFinite(value))
        return error.NonFiniteExecutionBalanceOutput;

    try writer.print("{d}\t{d}", .{
        record.execution_day,
        record.year,
    });
    for (values) |value| try writer.print("\t{e}", .{value});
    try writer.writeByte('\n');
}

test "EXEC diagnostic writes exact seven-domain tab record order" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeHeader(&bytes.writer);
    try writeRecord(&bytes.writer, .{
        .execution_day = 81,
        .year = 1998,
        .deviation = .{
            .water_m = 1,
            .heat_mj_per_m2 = 2,
            .oxygen_g_per_m2 = 3,
            .carbon_g_per_m2 = 4,
            .nitrogen_g_per_m2 = 5,
            .phosphorus_g_per_m2 = 6,
            .ions_mol_per_m2 = 7,
        },
    });
    try std.testing.expectEqualStrings(
        "execution_day\tyear\twater_balance_m\t" ++
            "heat_balance_MJ_m-2\toxygen_balance_g_m-2\t" ++
            "carbon_balance_g_m-2\tnitrogen_balance_g_m-2\t" ++
            "phosphorus_balance_g_m-2\tion_balance_mol_m-2\n" ++
            "81\t1998\t1e0\t2e0\t3e0\t4e0\t5e0\t6e0\t7e0\n",
        bytes.written(),
    );
}

test "invalid late domain writes no partial EXEC row" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try std.testing.expectError(
        error.NonFiniteExecutionBalanceOutput,
        writeRecord(&bytes.writer, .{
            .execution_day = 2,
            .year = 2001,
            .deviation = .{
                .water_m = 1,
                .heat_mj_per_m2 = 2,
                .oxygen_g_per_m2 = 3,
                .carbon_g_per_m2 = 4,
                .nitrogen_g_per_m2 = 5,
                .phosphorus_g_per_m2 = 6,
                .ions_mol_per_m2 = std.math.nan(f64),
            },
        }),
    );
    try std.testing.expectEqual(@as(usize, 0), bytes.written().len);
}

test "invalid execution day writes no diagnostic bytes" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try std.testing.expectError(
        error.InvalidExecutionBalanceOutputDay,
        writeRecord(&bytes.writer, .{
            .execution_day = 0,
            .year = 2001,
            .deviation = .{
                .water_m = 0,
                .heat_mj_per_m2 = 0,
                .oxygen_g_per_m2 = 0,
                .carbon_g_per_m2 = 0,
                .nitrogen_g_per_m2 = 0,
                .phosphorus_g_per_m2 = 0,
                .ions_mol_per_m2 = 0,
            },
        }),
    );
    try std.testing.expectEqual(@as(usize, 0), bytes.written().len);
}
