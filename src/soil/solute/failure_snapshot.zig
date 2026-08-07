const std = @import("std");
const chemistry = @import("chemistry_state.zig");
const solver = @import("reaction_solver.zig");

const magic = "ECOSSOL!";
const format_version: u32 = 1;
const checksum_seed: u64 = 0x45434f53534f4c46;

pub const Context = extern struct {
    execution_id: u64 = 0,
    scenario_id: u64 = 0,
    repeat_id: u64 = 0,
    scene_id: u64 = 0,
    scene_hour: u64 = 0,
    year: u64 = 0,
    day_of_year: u64 = 0,
    hour: u64 = 0,
    global_cell_id: u64 = 0,
    soil_layer_id: u64 = 0,
    packed_cell_index: u64 = 0,
};

pub const ReplayCase = struct {
    allocator: std.mem.Allocator,
    state: chemistry.State,
    parameters: chemistry.ReactionParameters,
    options: solver.Options,
    context: Context,

    pub fn deinit(self: *ReplayCase) void {
        self.state.deinit();
        self.* = undefined;
    }
};

pub fn capture(
    allocator: std.mem.Allocator,
    source: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: solver.Options,
    context: Context,
) !ReplayCase {
    if (cell_index >= source.cell_count)
        return error.ChemistryCellIndexOutOfBounds;
    var state = try chemistry.State.init(allocator, 1);
    errdefer state.deinit();
    const count = chemistry.State.packedComponentCount();
    const packed_values = try allocator.alloc(f64, count);
    defer allocator.free(packed_values);
    try source.packCell(cell_index, packed_values);
    try state.unpackCell(0, packed_values);
    return .{
        .allocator = allocator,
        .state = state,
        .parameters = parameters,
        .options = options,
        .context = context,
    };
}

pub fn write(
    allocator: std.mem.Allocator,
    writer: anytype,
    replay_case: *const ReplayCase,
) !void {
    var payload: std.Io.Writer.Allocating = .init(allocator);
    defer payload.deinit();
    try payload.writer.writeAll(std.mem.asBytes(&replay_case.context));
    try payload.writer.writeAll(std.mem.asBytes(&replay_case.parameters));
    try payload.writer.writeAll(std.mem.asBytes(&replay_case.options));
    var packed_values: [chemistry.State.packedComponentCount()]f64 = undefined;
    try replay_case.state.packCell(0, &packed_values);
    try payload.writer.writeAll(std.mem.asBytes(&packed_values));
    const bytes = payload.written();
    try writer.writeAll(magic);
    try writer.writeInt(u32, format_version, .little);
    try writer.writeInt(u64, @intCast(bytes.len), .little);
    try writer.writeInt(
        u64,
        std.hash.Wyhash.hash(checksum_seed, bytes),
        .little,
    );
    try writer.writeAll(bytes);
}

pub fn read(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
) !ReplayCase {
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic))
        return error.InvalidSoluteSnapshotMagic;
    if (try reader.takeInt(u32, .little) != format_version)
        return error.UnsupportedSoluteSnapshotVersion;
    const payload_size = std.math.cast(
        usize,
        try reader.takeInt(u64, .little),
    ) orelse return error.InvalidSoluteSnapshotSize;
    const checksum = try reader.takeInt(u64, .little);
    const expected_size = @sizeOf(Context) +
        @sizeOf(chemistry.ReactionParameters) +
        @sizeOf(solver.Options) +
        chemistry.State.packedComponentCount() * @sizeOf(f64);
    if (payload_size != expected_size)
        return error.InvalidSoluteSnapshotSize;
    const payload = try allocator.alloc(u8, payload_size);
    defer allocator.free(payload);
    reader.readSliceAll(payload) catch
        return error.TruncatedSoluteSnapshot;
    if (std.hash.Wyhash.hash(checksum_seed, payload) != checksum)
        return error.SoluteSnapshotChecksumMismatch;
    if (reader.peekByte()) |_|
        return error.TrailingSoluteSnapshotData
    else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }

    var context = std.mem.zeroes(Context);
    var parameters = std.mem.zeroes(chemistry.ReactionParameters);
    var options = std.mem.zeroes(solver.Options);
    var offset: usize = 0;
    @memcpy(std.mem.asBytes(&context), payload[offset..][0..@sizeOf(Context)]);
    offset += @sizeOf(Context);
    @memcpy(
        std.mem.asBytes(&parameters),
        payload[offset..][0..@sizeOf(chemistry.ReactionParameters)],
    );
    offset += @sizeOf(chemistry.ReactionParameters);
    @memcpy(
        std.mem.asBytes(&options),
        payload[offset..][0..@sizeOf(solver.Options)],
    );
    offset += @sizeOf(solver.Options);
    var packed_values: [chemistry.State.packedComponentCount()]f64 = undefined;
    @memcpy(std.mem.asBytes(&packed_values), payload[offset..]);
    var state = try chemistry.State.init(allocator, 1);
    errdefer state.deinit();
    try state.unpackCell(0, &packed_values);
    return .{
        .allocator = allocator,
        .state = state,
        .parameters = parameters,
        .options = options,
        .context = context,
    };
}

test "solute failure snapshot is self contained checksummed and replayable" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    var packed_values: [chemistry.State.packedComponentCount()]f64 =
        [_]f64{1} ** chemistry.State.packedComponentCount();
    try state.unpackCell(0, &packed_values);
    var parameters = std.mem.zeroes(chemistry.ReactionParameters);
    parameters.fractions = .{
        .ammonium_non_band = 1,
        .ammonium_band = 0,
        .nitrate_non_band = 1,
        .nitrate_band = 0,
        .phosphate_non_band = 1,
        .phosphate_band = 0,
    };
    const options: solver.Options = .{ .max_iterations = 1 };
    var captured = try capture(
        std.testing.allocator,
        &state,
        0,
        parameters,
        options,
        .{ .scene_hour = 3, .packed_cell_index = 0 },
    );
    defer captured.deinit();
    state.aqueous[0].hydroxide = 99;
    var encoded: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer encoded.deinit();
    try write(std.testing.allocator, &encoded.writer, &captured);
    var reader: std.Io.Reader = .fixed(encoded.written());
    var restored = try read(std.testing.allocator, &reader);
    defer restored.deinit();
    try std.testing.expectEqual(@as(u64, 3), restored.context.scene_hour);
    try std.testing.expectEqual(@as(f64, 1), restored.state.aqueous[0].hydroxide);
    try std.testing.expectEqual(options, restored.options);
    encoded.written()[encoded.written().len - 1] ^= 1;
    var corrupt_reader: std.Io.Reader = .fixed(encoded.written());
    try std.testing.expectError(
        error.SoluteSnapshotChecksumMismatch,
        read(std.testing.allocator, &corrupt_reader),
    );
}
