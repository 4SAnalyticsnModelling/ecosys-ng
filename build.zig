const std = @import("std");

const TestModulePaths = [_][]const u8{
    "src/util/utils.zig",
    "src/util/input_parser.zig",
    "src/io/geo_attributes.zig",
};

fn buildImports(b: *std.Build, names: []const []const u8, modules: []const *std.Build.Module) []const std.Build.Module.Import {
    std.debug.assert(names.len == modules.len);
    var imports = b.allocator.alloc(std.Build.Module.Import, names.len) catch @panic("module names!!");
    for (names, modules, 0..) |name, module, i| {
        imports[i] = .{ .name = name, .module = module };
    }
    return imports;
}

pub fn build(b: *std.Build) void {
    const target = b.graph.host;
    const optimize = b.standardOptimizeOption(.{});

    const mods = .{
        .utils = b.createModule(.{ .root_source_file = b.path("src/util/utils.zig") }),
        .input_parser = b.createModule(.{ .root_source_file = b.path("src/util/input_parser.zig") }),
        .geo_attr = b.createModule(.{ .root_source_file = b.path("src/io/geo_attributes.zig") }),
    };

    mods.input_parser.addImport("utils", mods.utils);
    mods.geo_attr.addImport("utils", mods.utils);
    mods.geo_attr.addImport("input_parser", mods.input_parser);

    const import_names = &[_][]const u8{ "utils", "input_parser", "geo_attr" };
    const import_modules = &[_]*std.Build.Module{ mods.utils, mods.input_parser, mods.geo_attr };
    const common_imports = buildImports(b, import_names, import_modules);

    const exe = b.addExecutable(.{
        .name = "ecosys-ng",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/ecosys-ng.zig"),
            .target = target,
            .optimize = optimize,
            .imports = common_imports,
        }),
    });

    // exe.stack_size = 16 * 1024 * 1024; // increase stack size to 16 MB to accommodate large arrays

    //custom binary folder with `zig build -p .` command
    const install_exe = b.addInstallArtifact(
        exe,
        .{
            .dest_dir = .{
                .override = .{ .custom = "ecosys-ng-bin" },
            },
        },
    );

    b.getInstallStep().dependOn(&install_exe.step);

    const run_exe = b.addRunArtifact(exe);

    const test_step = b.step("test", "Run ecosys-ng code test blocks");

    inline for (TestModulePaths) |path| {
        const test_blocks = b.addTest(.{
            .root_module = b.createModule(.{
                .root_source_file = b.path(path),
                .target = target,
                .optimize = optimize,
                .imports = common_imports,
            }),
        });
        const run_test_blocks = b.addRunArtifact(test_blocks);
        test_step.dependOn(&run_test_blocks.step);
    }

    test_step.dependOn(b.getInstallStep()); //create binary along with testing

    const run_step = b.step("run", "Run ecosys-ng application");
    run_step.dependOn(&run_exe.step);
}
