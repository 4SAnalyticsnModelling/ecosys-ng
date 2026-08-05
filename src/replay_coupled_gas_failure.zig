const std = @import("std");
const ecosys = @import("ecosys_ng");

/// Replays one self-contained pre-solve coupled-gas failure snapshot.
pub fn main(init: std.process.Init) !void {
    var debug_allocator: std.heap.DebugAllocator(.{}) = .init;
    defer _ = debug_allocator.deinit();
    const allocator = debug_allocator.allocator();
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2 and args.len != 3) {
        std.log.err(
            "usage: ecosys_ng_replay_gas <ecosys-ng-gas-failure-*.bin> [max_iterations]",
            .{},
        );
        return error.MissingCoupledGasFailureSnapshotPath;
    }
    // An explicit ceiling exists so the minimum sufficient iteration budget for
    // a real captured state can be measured. It only ever *lowers* the recorded
    // ceiling below what the snapshot carries: a replay may not manufacture
    // headroom the failing production solve did not have, because that would
    // turn a diagnosis into the very "widen the ceiling" move the migration
    // rules forbid.
    const requested_ceiling: ?u16 = if (args.len == 3)
        std.fmt.parseInt(u16, args[2], 10) catch
            return error.InvalidCoupledGasReplayIterationCeiling
    else
        null;
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
    var options = ecosys.coupled_gas_failure_reporter.replayOptions(replay_case.options);
    if (requested_ceiling) |ceiling| {
        if (ceiling == 0 or ceiling > replay_case.options.max_iterations) {
            std.log.err(
                "requested replay ceiling must be within the captured ceiling: requested={d} captured={d}",
                .{ ceiling, replay_case.options.max_iterations },
            );
            return error.InvalidCoupledGasReplayIterationCeiling;
        }
        options.max_iterations = ceiling;
    }
    // An archived version 1 snapshot did not record the REDIST bubble receiver
    // map, so the replayed system releases bubbles into the source cell. Say so
    // on every line, in both the success and failure paths, so no report can
    // cite a v1 replay as faithful evidence about bubbling without noticing.
    const bubbling_fidelity = if (replay_case.bubble_receiver_map_captured)
        "captured"
    else
        "not-captured-v1-source-cell-default";
    if (!replay_case.bubble_receiver_map_captured) {
        var any_bubbling = false;
        for (replay_case.bubbling_enabled) |enabled| any_bubbling = any_bubbling or enabled;
        if (any_bubbling) std.log.warn(
            "coupled gas replay snapshot predates bubble receiver capture: snapshot={s} bubbling_cells_present=true replayed_destination=source_cell",
            .{args[1]},
        );
    }
    const result = ecosys.coupled_gas_solver.solve(
        allocator,
        &replay_case.state,
        ecosys.coupled_gas_failure_reporter.replayInputs(&replay_case),
        options,
    ) catch |err| {
        std.log.err(
            "coupled gas replay failed reproducibly: snapshot={s} cells={d} faces={d} max_iterations={d} bubble_receiver_map={s} error={s}",
            .{
                args[1],
                replay_case.state.cell_count,
                replay_case.faces.len,
                options.max_iterations,
                bubbling_fidelity,
                @errorName(err),
            },
        );
        return err;
    };
    std.log.info(
        "coupled gas replay converged: snapshot={s} iterations={d} newton_steps={d} picard_steps={d} maximum_scaled_residual={e} max_iterations={d} bubble_receiver_map={s}",
        .{
            args[1],
            result.iterations,
            result.newton_raphson_steps,
            result.picard_steps,
            result.maximum_scaled_residual,
            options.max_iterations,
            bubbling_fidelity,
        },
    );
}
