const std = @import("std");
const chemistry = @import("solute_chemistry_state.zig");
const solver = @import("solute_reaction_solver.zig");
const snapshot = @import("solute_failure_snapshot.zig");

pub const Request = struct {
    io: std.Io,
    directory: std.Io.Dir,
    file_path: []const u8,
    context: snapshot.Context,
};

pub fn captureAndWrite(
    allocator: std.mem.Allocator,
    request: Request,
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: solver.Options,
) !void {
    var replay_case = try snapshot.capture(
        allocator,
        state,
        cell_index,
        parameters,
        options,
        request.context,
    );
    defer replay_case.deinit();
    const buffer = try allocator.alloc(u8, 64 * 1024);
    defer allocator.free(buffer);
    var atomic_file = try request.directory.createFileAtomic(
        request.io,
        request.file_path,
        .{ .replace = true },
    );
    defer atomic_file.deinit(request.io);
    var writer = atomic_file.file.writerStreaming(request.io, buffer);
    try snapshot.write(allocator, &writer.interface, &replay_case);
    try writer.interface.flush();
    try atomic_file.file.sync(request.io);
    try atomic_file.replace(request.io);
}

pub fn reportPreservingSolverError(
    allocator: std.mem.Allocator,
    request: Request,
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: solver.Options,
    solver_error: anyerror,
) anyerror {
    captureAndWrite(
        allocator,
        request,
        state,
        cell_index,
        parameters,
        options,
    ) catch |reporter_error| {
        std.log.warn(
            "SOLUTE failure snapshot write failed: path={s} reporter_error={s} preserved_solver_error={s}",
            .{
                request.file_path,
                @errorName(reporter_error),
                @errorName(solver_error),
            },
        );
    };
    return solver_error;
}
