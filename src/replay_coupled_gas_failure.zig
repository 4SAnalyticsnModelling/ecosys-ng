const std = @import("std");
const ecosys = @import("ecosys_ng");

/// Replays one self-contained pre-solve coupled-gas failure snapshot.
pub fn main(init: std.process.Init) !void {
    var debug_allocator: std.heap.DebugAllocator(.{}) = .init;
    defer _ = debug_allocator.deinit();
    const allocator = debug_allocator.allocator();
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.log.err(
            "usage: ecosys_ng_replay_gas <ecosys-ng-gas-failure-*.bin>",
            .{},
        );
        return error.MissingCoupledGasFailureSnapshotPath;
    }
    const bytes = try std.Io.Dir.cwd().readFileAlloc(
        init.io,
        args[1],
        allocator,
        .limited(256 * 1024 * 1024),
    );
    defer allocator.free(bytes);
    var reader: std.Io.Reader = .fixed(bytes);
    var replay_case = try ecosys.coupled_gas_failure_snapshot.read(
        allocator,
        &reader,
        .{},
    );
    defer replay_case.deinit();
    const result = ecosys.coupled_gas_solver.solve(
        allocator,
        &replay_case.state,
        ecosys.coupled_gas_failure_reporter.replayInputs(&replay_case),
        ecosys.coupled_gas_failure_reporter.replayOptions(replay_case.options),
    ) catch |err| {
        std.log.err(
            "coupled gas replay failed reproducibly: snapshot={s} cells={d} faces={d} max_iterations={d} error={s}",
            .{
                args[1],
                replay_case.state.cell_count,
                replay_case.faces.len,
                replay_case.options.max_iterations,
                @errorName(err),
            },
        );
        return err;
    };
    std.log.info(
        "coupled gas replay converged: snapshot={s} iterations={d} newton_steps={d} picard_steps={d} maximum_scaled_residual={e}",
        .{
            args[1],
            result.iterations,
            result.newton_raphson_steps,
            result.picard_steps,
            result.maximum_scaled_residual,
        },
    );
}
