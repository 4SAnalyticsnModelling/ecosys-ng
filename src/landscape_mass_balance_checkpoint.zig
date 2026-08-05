const std = @import("std");
const audit = @import("mass_balance_audit.zig");
const boundary = @import("landscape_boundary_ledger.zig");

const magic = "ECOSMBAL";
const version: u32 = 1;

pub const State = struct {
    boundary_ledger: boundary.State,
    monitor: ?audit.Monitor,
};

pub fn write(writer: anytype, state: State) !void {
    try validate(state);
    try writer.writeAll(magic);
    try writer.writeInt(u32, version, .little);
    inline for (std.meta.fields(boundary.Fluxes)) |field|
        try writer.writeInt(
            u64,
            @bitCast(@field(state.boundary_ledger.cumulative, field.name)),
            .little,
        );
    try writer.writeByte(@intFromBool(state.monitor != null));
    if (state.monitor) |monitor| {
        inline for (std.meta.fields(audit.Balance)) |field|
            try writer.writeInt(
                u64,
                @bitCast(@field(monitor.baseline, field.name)),
                .little,
            );
        try writer.writeInt(
            u64,
            @bitCast(monitor.tolerance_per_m2),
            .little,
        );
    }
}

pub fn read(reader: *std.Io.Reader) !State {
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic))
        return error.InvalidLandscapeMassBalanceCheckpointMagic;
    if (try reader.takeInt(u32, .little) != version)
        return error.UnsupportedLandscapeMassBalanceCheckpointVersion;
    var result: State = .{ .boundary_ledger = .{}, .monitor = null };
    inline for (std.meta.fields(boundary.Fluxes)) |field|
        @field(result.boundary_ledger.cumulative, field.name) =
            @bitCast(try reader.takeInt(u64, .little));
    result.monitor = switch (try reader.takeByte()) {
        0 => null,
        1 => monitor: {
            var baseline: audit.Balance = undefined;
            inline for (std.meta.fields(audit.Balance)) |field|
                @field(baseline, field.name) =
                    @bitCast(try reader.takeInt(u64, .little));
            break :monitor .{
                .baseline = baseline,
                .tolerance_per_m2 = @bitCast(try reader.takeInt(u64, .little)),
            };
        },
        else => return error.InvalidLandscapeMassBalanceMonitorTag,
    };
    if (reader.peekByte()) |_|
        return error.TrailingLandscapeMassBalanceCheckpointData
    else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
    try validate(result);
    return result;
}

fn validate(state: State) !void {
    inline for (std.meta.fields(boundary.Fluxes)) |field| {
        const value = @field(state.boundary_ledger.cumulative, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLandscapeMassBalanceCheckpointValue;
    }
    if (state.monitor) |monitor| {
        inline for (std.meta.fields(audit.Balance)) |field|
            if (!std.math.isFinite(@field(monitor.baseline, field.name)))
                return error.InvalidLandscapeMassBalanceCheckpointValue;
        if (!std.math.isFinite(monitor.tolerance_per_m2) or
            monitor.tolerance_per_m2 < 0)
            return error.InvalidLandscapeMassBalanceCheckpointValue;
    }
}

test "landscape mass balance checkpoint round trips cumulative history and monitor" {
    var bytes: [2048]u8 = undefined;
    var writer = std.Io.Writer.fixed(&bytes);
    var boundary_state: boundary.State = .{};
    try boundary_state.accumulateAccepted(.{
        .rain_m3 = 2.5,
        .heat_input_megajoules = 4,
        .dinitrogen_input_g_n = 0.3,
        .ion_output_mol = 0.02,
    });
    const state: State = .{
        .boundary_ledger = boundary_state,
        .monitor = .{
            .baseline = .{
                .water_m3 = -2.5,
                .heat_megajoules = -4,
                .oxygen_g = 0,
                .carbon_g = 0,
                .nitrogen_g = -0.3,
                .phosphorus_g = 0,
                .ions_mol = 0.02,
            },
            .tolerance_per_m2 = 1e-6,
        },
    };
    try write(&writer, state);
    var reader = std.Io.Reader.fixed(writer.buffered());
    const restored = try read(&reader);
    try std.testing.expectEqual(
        @as(f64, 2.5),
        restored.boundary_ledger.cumulative.rain_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.3),
        restored.boundary_ledger.cumulative.dinitrogen_input_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, -4),
        restored.monitor.?.baseline.heat_megajoules,
    );
    try std.testing.expectEqual(
        @as(f64, 1e-6),
        restored.monitor.?.tolerance_per_m2,
    );
}

test "landscape mass balance checkpoint rejects invalid and trailing state" {
    var invalid: State = .{ .boundary_ledger = .{}, .monitor = null };
    invalid.boundary_ledger.cumulative.rain_m3 = -1;
    var bytes: [2048]u8 = undefined;
    var invalid_writer = std.Io.Writer.fixed(&bytes);
    try std.testing.expectError(
        error.InvalidLandscapeMassBalanceCheckpointValue,
        write(&invalid_writer, invalid),
    );

    var writer = std.Io.Writer.fixed(&bytes);
    try write(&writer, .{ .boundary_ledger = .{}, .monitor = null });
    try writer.writeByte(99);
    var reader = std.Io.Reader.fixed(writer.buffered());
    try std.testing.expectError(
        error.TrailingLandscapeMassBalanceCheckpointData,
        read(&reader),
    );
}
