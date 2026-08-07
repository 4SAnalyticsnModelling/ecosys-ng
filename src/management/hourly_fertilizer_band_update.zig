const std = @import("std");
const geometry_module = @import("hourly_fertilizer_band_geometry.zig");
const repartition_module = @import("fertilizer_band_inventory_repartition.zig");

pub const Family = enum {
    ammonium,
    nitrate,
    phosphate,
};

pub const Update = struct {
    family: Family,
    geometry_state: *geometry_module.State,
    layer_geometry: geometry_module.LayerGeometry,
    forcing: geometry_module.Forcing,
    geometry_output: geometry_module.Output,
    pools_in_source_order: []const repartition_module.Pool,
};

pub const Result = struct {
    family: Family,
    geometry_changed_layer_count: usize,
    repartition: repartition_module.Result,
};

/// Reusable transaction storage for one cell and one nutrient family.
pub const Workspace = struct {
    allocator: std.mem.Allocator,
    geometry: geometry_module.Workspace,
    original_band_depth_m: []f64,
    original_band_width_m: []f64,
    original_band_volume_fraction: []f64,
    original_non_band_volume_fraction: []f64,
    original_relative_non_band_change: []f64,
    original_band_disappeared: []bool,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !Workspace {
        var geometry = try geometry_module.Workspace.init(
            allocator,
            layer_count,
        );
        errdefer geometry.deinit();
        const depth = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(depth);
        const width = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(width);
        const band = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(band);
        const non_band = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(non_band);
        const relative = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(relative);
        const disappeared = try allocator.alloc(bool, layer_count);
        return .{
            .allocator = allocator,
            .geometry = geometry,
            .original_band_depth_m = depth,
            .original_band_width_m = width,
            .original_band_volume_fraction = band,
            .original_non_band_volume_fraction = non_band,
            .original_relative_non_band_change = relative,
            .original_band_disappeared = disappeared,
        };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.original_band_disappeared);
        self.allocator.free(self.original_relative_non_band_change);
        self.allocator.free(self.original_non_band_volume_fraction);
        self.allocator.free(self.original_band_volume_fraction);
        self.allocator.free(self.original_band_width_m);
        self.allocator.free(self.original_band_depth_m);
        self.geometry.deinit();
        self.* = undefined;
    }
};

/// Composes HOUR1 geometry and SOLUTE repartition for one nutrient family.
///
/// This is an atomic standalone integration boundary: a repartition failure
/// restores geometry and FVL outputs exactly. It does not imply that the two
/// phases may be adjacent in production; legacy reactions occur between them.
pub fn updateFamily(
    update: Update,
    workspace: *Workspace,
) !Result {
    try validateWorkspace(update, workspace);
    const original_upper_edge_depth_m =
        update.geometry_state.upper_edge_depth_m;
    const original_lower_edge_depth_m =
        update.geometry_state.lower_edge_depth_m;
    @memcpy(
        workspace.original_band_depth_m,
        update.geometry_state.band_depth_m,
    );
    @memcpy(
        workspace.original_band_width_m,
        update.geometry_state.band_width_m,
    );
    @memcpy(
        workspace.original_band_volume_fraction,
        update.geometry_state.band_volume_fraction,
    );
    @memcpy(
        workspace.original_non_band_volume_fraction,
        update.geometry_state.non_band_volume_fraction,
    );
    @memcpy(
        workspace.original_relative_non_band_change,
        update.geometry_output.relative_non_band_change,
    );
    @memcpy(
        workspace.original_band_disappeared,
        update.geometry_output.band_disappeared,
    );

    try geometry_module.update(
        update.geometry_state,
        update.layer_geometry,
        update.forcing,
        update.geometry_output,
        &workspace.geometry,
    );

    const repartition = repartition_module.repartition(.{
        .relative_non_band_change = update.geometry_output.relative_non_band_change,
        .band_active = update.geometry_state.active,
        .first_active_layer = update.layer_geometry.first_active_layer,
        .last_active_layer = update.layer_geometry.last_active_layer,
    }, update.pools_in_source_order) catch |err| {
        restoreGeometry(
            update,
            workspace,
            original_upper_edge_depth_m,
            original_lower_edge_depth_m,
        );
        return err;
    };

    var changed_layer_count: usize = 0;
    for (
        update.geometry_state.band_volume_fraction,
        workspace.original_band_volume_fraction,
    ) |current, original| {
        if (current != original) changed_layer_count += 1;
    }
    return .{
        .family = update.family,
        .geometry_changed_layer_count = changed_layer_count,
        .repartition = repartition,
    };
}

fn validateWorkspace(update: Update, workspace: *const Workspace) !void {
    const count = update.geometry_state.band_depth_m.len;
    inline for (.{
        workspace.original_band_depth_m,
        workspace.original_band_width_m,
        workspace.original_band_volume_fraction,
        workspace.original_non_band_volume_fraction,
        workspace.original_relative_non_band_change,
    }) |values| if (values.len != count)
        return error.FertilizerBandUpdateWorkspaceDimensionMismatch;
    if (workspace.original_band_disappeared.len != count)
        return error.FertilizerBandUpdateWorkspaceDimensionMismatch;
}

fn restoreGeometry(
    update: Update,
    workspace: *const Workspace,
    upper_edge_depth_m: f64,
    lower_edge_depth_m: f64,
) void {
    @memcpy(
        update.geometry_state.band_depth_m,
        workspace.original_band_depth_m,
    );
    @memcpy(
        update.geometry_state.band_width_m,
        workspace.original_band_width_m,
    );
    @memcpy(
        update.geometry_state.band_volume_fraction,
        workspace.original_band_volume_fraction,
    );
    @memcpy(
        update.geometry_state.non_band_volume_fraction,
        workspace.original_non_band_volume_fraction,
    );
    @memcpy(
        update.geometry_output.relative_non_band_change,
        workspace.original_relative_non_band_change,
    );
    @memcpy(
        update.geometry_output.band_disappeared,
        workspace.original_band_disappeared,
    );
    update.geometry_state.upper_edge_depth_m = upper_edge_depth_m;
    update.geometry_state.lower_edge_depth_m = lower_edge_depth_m;
}

fn fixtureState(
    depth: []f64,
    width: []f64,
    band: []f64,
    non_band: []f64,
) geometry_module.State {
    return .{
        .active = true,
        .row_spacing_m = 0.1,
        .upper_edge_depth_m = 0.04,
        .lower_edge_depth_m = 0.06,
        .band_depth_m = depth,
        .band_width_m = width,
        .band_volume_fraction = band,
        .non_band_volume_fraction = non_band,
    };
}

fn fixtureGeometry() geometry_module.LayerGeometry {
    return .{
        .upper_depth_m = &.{ 0, 0.1 },
        .lower_depth_m = &.{ 0.1, 0.2 },
        .thickness_m = &.{ 0.1, 0.1 },
        .first_active_layer = 0,
        .last_active_layer = 1,
        .minimum_active_thickness_m = 1e-6,
        .structural_presence_threshold = 1e-12,
    };
}

fn fixtureForcing() geometry_module.Forcing {
    return .{
        .diffusivity_m2_per_h = &.{ 0.0001, 0.0001 },
        .tortuosity = &.{ 1, 1 },
        .timestep_h = 1,
    };
}

test "geometry FVL feeds conservative SOLUTE repartition in exact local order" {
    var depth = [_]f64{ 0.02, 0 };
    var width = [_]f64{ 0.02, 0 };
    var band_fraction = [_]f64{ 0.004, 0 };
    var non_band_fraction = [_]f64{ 0.996, 1 };
    var relative = [_]f64{ 0, 0 };
    var disappeared = [_]bool{ false, false };
    var state = fixtureState(
        &depth,
        &width,
        &band_fraction,
        &non_band_fraction,
    );
    const dissolved_g_n = [_]f64{ 14, 0 };
    var non_band_change_mol = [_]f64{ 0, 0 };
    var band_change_mol = [_]f64{ 0, 0 };
    var broadcast_fertilizer_mol = [_]f64{ 10, 0 };
    var banded_fertilizer_mol = [_]f64{ 2, 0 };
    const pools = [_]repartition_module.Pool{
        .{
            .name = "dissolved ammonium",
            .inventory_unit = .grams_nitrogen,
            .storage = .{ .deferred = .{
                .non_band_inventory = &dissolved_g_n,
                .non_band_change = &non_band_change_mol,
                .band_change = &band_change_mol,
                .change_unit = .moles,
                .inventory_units_per_change_unit = 14,
            } },
        },
        .{
            .name = "undissolved ammonium",
            .inventory_unit = .moles,
            .storage = .{ .immediate = .{
                .non_band_inventory = &broadcast_fertilizer_mol,
                .band_inventory = &banded_fertilizer_mol,
            } },
        },
    };
    var workspace = try Workspace.init(std.testing.allocator, 2);
    defer workspace.deinit();

    const result = try updateFamily(.{
        .family = .ammonium,
        .geometry_state = &state,
        .layer_geometry = fixtureGeometry(),
        .forcing = fixtureForcing(),
        .geometry_output = .{
            .relative_non_band_change = &relative,
            .band_disappeared = &disappeared,
        },
        .pools_in_source_order = &pools,
    }, &workspace);

    try std.testing.expectEqual(Family.ammonium, result.family);
    try std.testing.expectEqual(@as(usize, 1), result.repartition.layers_repartitioned);
    try std.testing.expect(relative[0] < 0);
    try std.testing.expectApproxEqAbs(
        relative[0],
        non_band_change_mol[0],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        non_band_change_mol[0] + band_change_mol[0],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 12),
        broadcast_fertilizer_mol[0] + banded_fertilizer_mol[0],
        1e-14,
    );
}

test "late repartition rejection restores geometry and output atomically" {
    var depth = [_]f64{ 0.02, 0 };
    var width = [_]f64{ 0.02, 0 };
    var band_fraction = [_]f64{ 0.004, 0 };
    var non_band_fraction = [_]f64{ 0.996, 1 };
    const depth_before = depth;
    const width_before = width;
    const band_before = band_fraction;
    const non_band_before = non_band_fraction;
    var relative = [_]f64{ 7, 8 };
    var disappeared = [_]bool{ true, true };
    const relative_before = relative;
    const disappeared_before = disappeared;
    var state = fixtureState(
        &depth,
        &width,
        &band_fraction,
        &non_band_fraction,
    );
    const upper_before = state.upper_edge_depth_m;
    const lower_before = state.lower_edge_depth_m;
    var invalid_non = [_]f64{ std.math.nan(f64), 0 };
    var invalid_band = [_]f64{ 0, 0 };
    const pools = [_]repartition_module.Pool{
        .{
            .name = "invalid fertilizer",
            .inventory_unit = .moles,
            .storage = .{ .immediate = .{
                .non_band_inventory = &invalid_non,
                .band_inventory = &invalid_band,
            } },
        },
    };
    var workspace = try Workspace.init(std.testing.allocator, 2);
    defer workspace.deinit();

    try std.testing.expectError(error.InvalidBandInventoryState, updateFamily(.{
        .family = .nitrate,
        .geometry_state = &state,
        .layer_geometry = fixtureGeometry(),
        .forcing = fixtureForcing(),
        .geometry_output = .{
            .relative_non_band_change = &relative,
            .band_disappeared = &disappeared,
        },
        .pools_in_source_order = &pools,
    }, &workspace));
    try std.testing.expectEqualSlices(f64, &depth_before, &depth);
    try std.testing.expectEqualSlices(f64, &width_before, &width);
    try std.testing.expectEqualSlices(f64, &band_before, &band_fraction);
    try std.testing.expectEqualSlices(f64, &non_band_before, &non_band_fraction);
    try std.testing.expectEqualSlices(f64, &relative_before, &relative);
    try std.testing.expectEqualSlices(bool, &disappeared_before, &disappeared);
    try std.testing.expectEqual(upper_before, state.upper_edge_depth_m);
    try std.testing.expectEqual(lower_before, state.lower_edge_depth_m);
}
