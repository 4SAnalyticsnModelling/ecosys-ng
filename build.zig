const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const module = b.addModule("ecosys_ng", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    const executable = b.addExecutable(.{
        .name = "ecosys_ng",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/ecosys_ng.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "ecosys_ng", .module = module }},
        }),
    });
    b.installArtifact(executable);
    const gas_replay_executable = b.addExecutable(.{
        .name = "ecosys_ng_replay_gas",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/replay_coupled_gas_failure.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "ecosys_ng", .module = module }},
        }),
    });
    b.installArtifact(gas_replay_executable);
    const solute_replay_executable = b.addExecutable(.{
        .name = "ecosys_ng_replay_solute",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/replay_solute_failure.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "ecosys_ng", .module = module }},
        }),
    });
    b.installArtifact(solute_replay_executable);
    const check = b.step("check", "Compile ecosys-ng without installing it");
    check.dependOn(&executable.step);
    check.dependOn(&gas_replay_executable.step);
    check.dependOn(&solute_replay_executable.step);

    const run = b.addRunArtifact(executable);
    run.step.dependOn(b.getInstallStep());
    if (b.args) |args| run.addArgs(args);
    b.step("run", "Run ecosys-ng").dependOn(&run.step);

    const replay_gas = b.addRunArtifact(gas_replay_executable);
    replay_gas.step.dependOn(b.getInstallStep());
    if (b.args) |args| replay_gas.addArgs(args);
    b.step(
        "replay-gas",
        "Replay one coupled-gas failure snapshot",
    ).dependOn(&replay_gas.step);
    const replay_solute = b.addRunArtifact(solute_replay_executable);
    replay_solute.step.dependOn(b.getInstallStep());
    if (b.args) |args| replay_solute.addArgs(args);
    b.step(
        "replay-solute",
        "Replay one SOLUTE failure snapshot",
    ).dependOn(&replay_solute.step);

    const tests = b.addTest(.{ .root_module = module });
    const run_tests = b.addRunArtifact(tests);
    b.step("test", "Run unit tests").dependOn(&run_tests.step);

    const freezing_column_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/production_freezing_column_validation.zig"),
            .target = target,
            .optimize = .ReleaseFast,
        }),
    });
    const run_freezing_column_tests =
        b.addRunArtifact(freezing_column_tests);
    run_freezing_column_tests.addArgs(
        &.{ "--test-filter", "production dual-phase" },
    );
    b.step(
        "validate-freezing-column",
        "Run the published 500-cell production freeze-thaw validation",
    ).dependOn(&run_freezing_column_tests.step);
}
