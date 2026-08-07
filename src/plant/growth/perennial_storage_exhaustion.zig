const std = @import("std");

pub const State = struct {
    branch_dead: []bool,
    shoot_dead: []bool,
    root_dead: []bool,
    population_density_per_m2: []f64,
    population_count: []f64,
    canopy_water_m3: []f64,
    canopy_surface_water_m3: []f64,
    standing_dead_water_m3: []f64,
};

pub const Inputs = struct {
    is_perennial: bool,
    storage_carbon_g_c: f64,
    storage_nitrogen_g_n: f64,
    storage_phosphorus_g_p: f64,
    population_scaled_presence_threshold: f64,
};

fn validate(state: State, inputs: Inputs) !void {
    if (state.shoot_dead.len != 1 or state.root_dead.len != 1 or state.population_density_per_m2.len != 1 or state.population_count.len != 1 or state.canopy_water_m3.len != 1 or state.canopy_surface_water_m3.len != 1 or state.standing_dead_water_m3.len != 1) return error.PerennialStorageExhaustionDimensionMismatch;
    inline for (.{ inputs.storage_carbon_g_c, inputs.storage_nitrogen_g_n, inputs.storage_phosphorus_g_p }) |value| if (!std.math.isFinite(value)) return error.InvalidPerennialStorageExhaustionState;
    inline for (.{ inputs.population_scaled_presence_threshold, state.population_density_per_m2[0], state.population_count[0], state.canopy_water_m3[0], state.canopy_surface_water_m3[0], state.standing_dead_water_m3[0] }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPerennialStorageExhaustionState;
}

/// Exact GROSUB 7535--7553 DEATHR state transition. A perennial dies when any
/// seasonal C/N/P pool is at or below the population-scaled threshold.
pub fn apply(state: State, workspace: State, inputs: Inputs) !bool {
    try validate(state, inputs);
    try validate(workspace, inputs);
    @memcpy(workspace.branch_dead, state.branch_dead);
    inline for (.{ "shoot_dead", "root_dead", "population_density_per_m2", "population_count", "canopy_water_m3", "canopy_surface_water_m3", "standing_dead_water_m3" }) |field_name| @memcpy(@field(workspace, field_name), @field(state, field_name));
    const dies = inputs.is_perennial and (inputs.storage_carbon_g_c <= inputs.population_scaled_presence_threshold or inputs.storage_nitrogen_g_n <= inputs.population_scaled_presence_threshold or inputs.storage_phosphorus_g_p <= inputs.population_scaled_presence_threshold);
    if (dies) {
        @memset(workspace.branch_dead, true);
        workspace.shoot_dead[0] = true;
        workspace.root_dead[0] = true;
        workspace.population_density_per_m2[0] = 0;
        workspace.population_count[0] = 0;
        workspace.standing_dead_water_m3[0] += workspace.canopy_surface_water_m3[0] + workspace.canopy_water_m3[0];
        workspace.canopy_surface_water_m3[0] = 0;
        workspace.canopy_water_m3[0] = 0;
        if (!std.math.isFinite(workspace.standing_dead_water_m3[0])) return error.NonFinitePerennialStorageExhaustionResult;
    }
    @memcpy(state.branch_dead, workspace.branch_dead);
    inline for (.{ "shoot_dead", "root_dead", "population_density_per_m2", "population_count", "canopy_water_m3", "canopy_surface_water_m3", "standing_dead_water_m3" }) |field_name| @memcpy(@field(state, field_name), @field(workspace, field_name));
    return dies;
}

test "GROSUB perennial storage exhaustion kills all runtime branches and conserves water" {
    var branches = [_]bool{ false, false, false, false, false, false, false };
    var flags = [_]bool{false};
    var density = [_]f64{2};
    var population = [_]f64{20};
    var canopy = [_]f64{3};
    var surface = [_]f64{4};
    var dead = [_]f64{5};
    var wb = [_]bool{false} ** 7;
    var wf1 = [_]bool{false};
    var wf2 = [_]bool{false};
    var wd = [_]f64{0};
    var wp = [_]f64{0};
    var wc = [_]f64{0};
    var ws = [_]f64{0};
    var wdead = [_]f64{0};
    const state = State{ .branch_dead = &branches, .shoot_dead = &flags, .root_dead = &flags, .population_density_per_m2 = &density, .population_count = &population, .canopy_water_m3 = &canopy, .canopy_surface_water_m3 = &surface, .standing_dead_water_m3 = &dead };
    const work = State{ .branch_dead = &wb, .shoot_dead = &wf1, .root_dead = &wf2, .population_density_per_m2 = &wd, .population_count = &wp, .canopy_water_m3 = &wc, .canopy_surface_water_m3 = &ws, .standing_dead_water_m3 = &wdead };
    try std.testing.expect(try apply(state, work, .{ .is_perennial = true, .storage_carbon_g_c = 0.1, .storage_nitrogen_g_n = 0, .storage_phosphorus_g_p = 0.1, .population_scaled_presence_threshold = 0 }));
    for (branches) |dead_branch| try std.testing.expect(dead_branch);
    try std.testing.expectApproxEqAbs(12, dead[0], 1e-12);
    try std.testing.expectEqual(@as(f64, 0), canopy[0]);
    try std.testing.expectEqual(@as(f64, 0), surface[0]);
}

test "GROSUB annual does not die from exhausted seasonal storage" {
    var branches = [_]bool{false};
    var shoot = [_]bool{false};
    var root = [_]bool{false};
    var density = [_]f64{1};
    var population = [_]f64{1};
    var canopy = [_]f64{1};
    var surface = [_]f64{1};
    var dead = [_]f64{1};
    var wb = branches;
    var wsht = shoot;
    var wr = root;
    var wd = [_]f64{0};
    var wp = [_]f64{0};
    var wc = [_]f64{0};
    var wsurf = [_]f64{0};
    var wdead = [_]f64{0};
    try std.testing.expect(!try apply(.{ .branch_dead = &branches, .shoot_dead = &shoot, .root_dead = &root, .population_density_per_m2 = &density, .population_count = &population, .canopy_water_m3 = &canopy, .canopy_surface_water_m3 = &surface, .standing_dead_water_m3 = &dead }, .{ .branch_dead = &wb, .shoot_dead = &wsht, .root_dead = &wr, .population_density_per_m2 = &wd, .population_count = &wp, .canopy_water_m3 = &wc, .canopy_surface_water_m3 = &wsurf, .standing_dead_water_m3 = &wdead }, .{ .is_perennial = false, .storage_carbon_g_c = 0, .storage_nitrogen_g_n = 0, .storage_phosphorus_g_p = 0, .population_scaled_presence_threshold = 0 }));
    try std.testing.expectEqual(@as(f64, 1), population[0]);
    try std.testing.expectEqual(@as(f64, 1), dead[0]);
}

test "GROSUB negative perennial storage triggers death instead of input rejection" {
    var branches = [_]bool{false};
    var flags = [_]bool{false};
    var density = [_]f64{1};
    var population = [_]f64{1};
    var canopy = [_]f64{0};
    var surface = [_]f64{0};
    var dead = [_]f64{0};
    var wb = branches;
    var wf1 = flags;
    var wf2 = flags;
    var wd = [_]f64{0};
    var wp = [_]f64{0};
    var wc = [_]f64{0};
    var ws = [_]f64{0};
    var wdead = [_]f64{0};
    try std.testing.expect(try apply(.{ .branch_dead = &branches, .shoot_dead = &flags, .root_dead = &flags, .population_density_per_m2 = &density, .population_count = &population, .canopy_water_m3 = &canopy, .canopy_surface_water_m3 = &surface, .standing_dead_water_m3 = &dead }, .{ .branch_dead = &wb, .shoot_dead = &wf1, .root_dead = &wf2, .population_density_per_m2 = &wd, .population_count = &wp, .canopy_water_m3 = &wc, .canopy_surface_water_m3 = &ws, .standing_dead_water_m3 = &wdead }, .{ .is_perennial = true, .storage_carbon_g_c = -0.1, .storage_nitrogen_g_n = 1, .storage_phosphorus_g_p = 1, .population_scaled_presence_threshold = 0 }));
}
