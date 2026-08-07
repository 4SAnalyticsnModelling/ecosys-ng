const std = @import("std");

pub const GrowthHabit = enum { annual, perennial };
pub const SaltMode = enum { static_equilibrium, dynamic };
pub const standard_salt_species_count: usize = 8;

pub const State = struct {
    branch_mobile_carbon_g_c: []f64,
    branch_mobile_nitrogen_g_n: []f64,
    branch_mobile_phosphorus_g_p: []f64,
    root_mobile_carbon_g_c_by_layer: []f64,
    root_mobile_nitrogen_g_n_by_layer: []f64,
    root_mobile_phosphorus_g_p_by_layer: []f64,
    branch_salt_mol_by_species_branch: []f64,
    root_salt_mol_by_species_layer: []f64,
    total_leaf_petiole_carbon_g_c: []f64,
};

pub const Inputs = struct {
    branch_count: usize,
    soil_layer_count: usize,
    salt_species_count: usize = standard_salt_species_count,
    shared_planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    branch_is_alive: []const bool,
    branch_leaf_petiole_carbon_g_c: []const f64,
    active_root_carbon_g_c_by_layer: []const f64,
    root_layer_sink_strength: []const f64,
    total_root_sink_strength: f64,
    total_root_carbon_g_c: f64,
    branch_nonwoody_carbon_fraction: f64,
    root_nonwoody_carbon_fraction: f64,
    minimum_sink_ratio: f64,
    base_exchange_rate_per_h: f64,
    leaf_petiole_allocation_fraction: f64,
    growth_habit: GrowthHabit,
    salt_mode: SaltMode,
    salt_exchange_rate_per_h: f64,
    biological_timestep_h: f64,
    minimum_pool_g_c: f64,
};

fn branchSaltIndex(inputs: Inputs, species: usize, branch: usize) usize {
    return species * inputs.branch_count + branch;
}

fn rootSaltIndex(inputs: Inputs, species: usize, layer: usize) usize {
    return species * inputs.soil_layer_count + layer;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, inputs: Inputs) !void {
    const branch_salt_count = std.math.mul(usize, inputs.salt_species_count, inputs.branch_count) catch return error.ShootRootExchangeDimensionOverflow;
    const root_salt_count = std.math.mul(usize, inputs.salt_species_count, inputs.soil_layer_count) catch return error.ShootRootExchangeDimensionOverflow;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const expected = if (std.mem.startsWith(u8, field.name, "branch_mobile")) inputs.branch_count else if (std.mem.startsWith(u8, field.name, "root_mobile")) inputs.soil_layer_count else if (std.mem.startsWith(u8, field.name, "branch_salt")) branch_salt_count else if (std.mem.startsWith(u8, field.name, "root_salt")) root_salt_count else 1;
        if (values.len != expected) return error.ShootRootExchangeDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidShootRootExchangeState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.branch_count == 0 or inputs.soil_layer_count == 0 or inputs.salt_species_count == 0 or inputs.shared_planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or
        inputs.branch_is_alive.len != inputs.branch_count or inputs.branch_leaf_petiole_carbon_g_c.len != inputs.branch_count or inputs.active_root_carbon_g_c_by_layer.len != inputs.soil_layer_count or inputs.root_layer_sink_strength.len != inputs.soil_layer_count) return error.ShootRootExchangeDimensionMismatch;
    const nonnegative = [_]f64{ inputs.total_root_sink_strength, inputs.total_root_carbon_g_c, inputs.branch_nonwoody_carbon_fraction, inputs.root_nonwoody_carbon_fraction, inputs.minimum_sink_ratio, inputs.base_exchange_rate_per_h, inputs.leaf_petiole_allocation_fraction, inputs.salt_exchange_rate_per_h, inputs.biological_timestep_h, inputs.minimum_pool_g_c };
    for (nonnegative) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidShootRootExchangeInput;
    inline for (.{ inputs.branch_leaf_petiole_carbon_g_c, inputs.active_root_carbon_g_c_by_layer, inputs.root_layer_sink_strength }) |values|
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidShootRootExchangeInput;
    inline for (.{ inputs.branch_nonwoody_carbon_fraction, inputs.root_nonwoody_carbon_fraction, inputs.minimum_sink_ratio, inputs.leaf_petiole_allocation_fraction }) |value| if (value > 1.0) return error.InvalidShootRootExchangeInput;
}

fn exchangeSalt(workspace: State, inputs: Inputs, species: usize, branch: usize, layer: usize, branch_carbon_g_c: f64, root_carbon_g_c: f64, combined_carbon_g_c: f64) void {
    const branch_index = branchSaltIndex(inputs, species, branch);
    const root_index = rootSaltIndex(inputs, species, layer);
    const branch_salt_mol = @max(0.0, workspace.branch_salt_mol_by_species_branch[branch_index]);
    const root_salt_mol = @max(0.0, workspace.root_salt_mol_by_species_layer[root_index]);
    const difference_mol = (branch_salt_mol * root_carbon_g_c - root_salt_mol * branch_carbon_g_c) / combined_carbon_g_c;
    const transfer_mol = inputs.salt_exchange_rate_per_h * difference_mol * inputs.biological_timestep_h;
    workspace.branch_salt_mol_by_species_branch[branch_index] -= transfer_mol;
    workspace.root_salt_mol_by_species_layer[root_index] += transfer_mol;
}

/// Exact GROSUB 8241--8438 shoot/root sink weighting and mobile C:N:P/salt
/// exchange. Runtime topology is branches x shared NU..NI layers. Masses use
/// g C/N/P or mol salt; exchange rates are h-1 and timestep is h.
pub fn apply(allocator: std.mem.Allocator, state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    copyState(workspace, state);
    const root_layer_weight = try allocator.alloc(f64, inputs.soil_layer_count);
    defer allocator.free(root_layer_weight);
    @memset(root_layer_weight, 0.0);
    const branch_weight = try allocator.alloc(f64, inputs.branch_count);
    defer allocator.free(branch_weight);
    @memset(branch_weight, 0.0);

    const preceding_total_leaf_carbon_g_c = state.total_leaf_petiole_carbon_g_c[0];
    const canopy_weight = if (preceding_total_leaf_carbon_g_c > inputs.minimum_pool_g_c) @min(1.0, 0.667 * inputs.total_root_carbon_g_c / preceding_total_leaf_carbon_g_c) else 1.0;
    const root_weight = if (inputs.total_root_carbon_g_c > inputs.minimum_pool_g_c) @min(1.0, preceding_total_leaf_carbon_g_c / (0.667 * inputs.total_root_carbon_g_c)) else 1.0;
    for (inputs.shared_planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer|
        root_layer_weight[layer] = if (inputs.total_root_sink_strength > inputs.minimum_pool_g_c) @max(0.0, inputs.root_layer_sink_strength[layer] / inputs.total_root_sink_strength) else 1.0;

    workspace.total_leaf_petiole_carbon_g_c[0] = 0.0;
    for (inputs.branch_leaf_petiole_carbon_g_c) |carbon_g_c| workspace.total_leaf_petiole_carbon_g_c[0] += carbon_g_c;
    for (0..inputs.branch_count) |branch| {
        if (!inputs.branch_is_alive[branch]) continue;
        branch_weight[branch] = if (workspace.total_leaf_petiole_carbon_g_c[0] > inputs.minimum_pool_g_c) @max(0.0, inputs.branch_leaf_petiole_carbon_g_c[branch] / workspace.total_leaf_petiole_carbon_g_c[0]) else 1.0;
        const carbon_exchange_rate_per_h = if (inputs.growth_habit == .annual) @max(5.0e-3, inputs.base_exchange_rate_per_h * std.math.pow(f64, inputs.leaf_petiole_allocation_fraction, 0.25)) else inputs.base_exchange_rate_per_h;
        for (inputs.shared_planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            const weighted_branch_carbon_g_c = inputs.branch_leaf_petiole_carbon_g_c[branch] * inputs.branch_nonwoody_carbon_fraction * root_layer_weight[layer] * canopy_weight;
            const weighted_root_carbon_g_c = inputs.active_root_carbon_g_c_by_layer[layer] * inputs.root_nonwoody_carbon_fraction * branch_weight[branch] * root_weight;
            const effective_branch_carbon_g_c = @max(0.0, weighted_branch_carbon_g_c, inputs.minimum_sink_ratio * weighted_root_carbon_g_c);
            const effective_root_carbon_g_c = @max(0.0, weighted_root_carbon_g_c, inputs.minimum_sink_ratio * weighted_branch_carbon_g_c);
            const combined_effective_carbon_g_c = effective_branch_carbon_g_c + effective_root_carbon_g_c;
            if (combined_effective_carbon_g_c <= inputs.minimum_pool_g_c) continue;
            const branch_mobile_carbon_g_c = @max(0.0, workspace.branch_mobile_carbon_g_c[branch] * root_layer_weight[layer]);
            const root_mobile_carbon_g_c = @max(0.0, workspace.root_mobile_carbon_g_c_by_layer[layer] * branch_weight[branch]);
            const carbon_difference_g_c = (branch_mobile_carbon_g_c * effective_root_carbon_g_c - root_mobile_carbon_g_c * effective_branch_carbon_g_c) / combined_effective_carbon_g_c;
            const carbon_transfer_g_c = carbon_exchange_rate_per_h * carbon_difference_g_c * inputs.biological_timestep_h;
            workspace.branch_mobile_carbon_g_c[branch] -= carbon_transfer_g_c;
            workspace.root_mobile_carbon_g_c_by_layer[layer] += carbon_transfer_g_c;
            const combined_mobile_carbon_g_c = root_mobile_carbon_g_c + branch_mobile_carbon_g_c;
            var nitrogen_transfer_g_n: f64 = 0.0;
            var phosphorus_transfer_g_p: f64 = 0.0;
            if (combined_mobile_carbon_g_c > inputs.minimum_pool_g_c) {
                const branch_nitrogen_g_n = @max(0.0, workspace.branch_mobile_nitrogen_g_n[branch] * root_layer_weight[layer]);
                const root_nitrogen_g_n = @max(0.0, workspace.root_mobile_nitrogen_g_n_by_layer[layer] * branch_weight[branch]);
                const nitrogen_difference_g_n = (branch_nitrogen_g_n * root_mobile_carbon_g_c - root_nitrogen_g_n * branch_mobile_carbon_g_c) / combined_mobile_carbon_g_c;
                nitrogen_transfer_g_n = inputs.base_exchange_rate_per_h * nitrogen_difference_g_n * inputs.biological_timestep_h;
                const branch_phosphorus_g_p = @max(0.0, workspace.branch_mobile_phosphorus_g_p[branch] * root_layer_weight[layer]);
                const root_phosphorus_g_p = @max(0.0, workspace.root_mobile_phosphorus_g_p_by_layer[layer] * branch_weight[branch]);
                const phosphorus_difference_g_p = (branch_phosphorus_g_p * root_mobile_carbon_g_c - root_phosphorus_g_p * branch_mobile_carbon_g_c) / combined_mobile_carbon_g_c;
                phosphorus_transfer_g_p = inputs.base_exchange_rate_per_h * phosphorus_difference_g_p * inputs.biological_timestep_h;
            }
            workspace.branch_mobile_nitrogen_g_n[branch] -= nitrogen_transfer_g_n;
            workspace.root_mobile_nitrogen_g_n_by_layer[layer] += nitrogen_transfer_g_n;
            workspace.branch_mobile_phosphorus_g_p[branch] -= phosphorus_transfer_g_p;
            workspace.root_mobile_phosphorus_g_p_by_layer[layer] += phosphorus_transfer_g_p;
            if (inputs.salt_mode == .dynamic) {
                const raw_branch_carbon_g_c = @max(0.0, inputs.branch_leaf_petiole_carbon_g_c[branch]);
                const raw_root_carbon_g_c = @max(0.0, inputs.active_root_carbon_g_c_by_layer[layer]);
                const combined_raw_carbon_g_c = raw_branch_carbon_g_c + raw_root_carbon_g_c;
                if (combined_raw_carbon_g_c <= 0.0) return error.InvalidShootRootSaltExchangeDenominator;
                for (0..inputs.salt_species_count) |species| exchangeSalt(workspace, inputs, species, branch, layer, raw_branch_carbon_g_c, raw_root_carbon_g_c, combined_raw_carbon_g_c);
            }
        }
    }
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

test "GROSUB shoot-root exchange conserves C N P and salts" {
    var branch_c = [_]f64{ 8, 2 };
    var branch_n = [_]f64{ 0.8, 0.2 };
    var branch_p = [_]f64{ 0.08, 0.02 };
    var root_c = [_]f64{ 1, 3 };
    var root_n = [_]f64{ 0.1, 0.3 };
    var root_p = [_]f64{ 0.01, 0.03 };
    var branch_salts = [_]f64{ 4, 1 } ** standard_salt_species_count;
    var root_salts = [_]f64{ 1, 2 } ** standard_salt_species_count;
    var total_leaf = [_]f64{10};
    var wc = [_]f64{0} ** 2;
    var wn = [_]f64{0} ** 2;
    var wp = [_]f64{0} ** 2;
    var wrc = [_]f64{0} ** 2;
    var wrn = [_]f64{0} ** 2;
    var wrp = [_]f64{0} ** 2;
    var wbs = [_]f64{0} ** (standard_salt_species_count * 2);
    var wrs = [_]f64{0} ** (standard_salt_species_count * 2);
    var wtl = [_]f64{0};
    const state: State = .{ .branch_mobile_carbon_g_c = &branch_c, .branch_mobile_nitrogen_g_n = &branch_n, .branch_mobile_phosphorus_g_p = &branch_p, .root_mobile_carbon_g_c_by_layer = &root_c, .root_mobile_nitrogen_g_n_by_layer = &root_n, .root_mobile_phosphorus_g_p_by_layer = &root_p, .branch_salt_mol_by_species_branch = &branch_salts, .root_salt_mol_by_species_layer = &root_salts, .total_leaf_petiole_carbon_g_c = &total_leaf };
    const workspace: State = .{ .branch_mobile_carbon_g_c = &wc, .branch_mobile_nitrogen_g_n = &wn, .branch_mobile_phosphorus_g_p = &wp, .root_mobile_carbon_g_c_by_layer = &wrc, .root_mobile_nitrogen_g_n_by_layer = &wrn, .root_mobile_phosphorus_g_p_by_layer = &wrp, .branch_salt_mol_by_species_branch = &wbs, .root_salt_mol_by_species_layer = &wrs, .total_leaf_petiole_carbon_g_c = &wtl };
    const inputs: Inputs = .{ .branch_count = 2, .soil_layer_count = 2, .shared_planting_layer_index = 0, .deepest_rooted_layer_index = 1, .branch_is_alive = &.{ true, true }, .branch_leaf_petiole_carbon_g_c = &.{ 8, 2 }, .active_root_carbon_g_c_by_layer = &.{ 2, 3 }, .root_layer_sink_strength = &.{ 1, 1 }, .total_root_sink_strength = 2, .total_root_carbon_g_c = 5, .branch_nonwoody_carbon_fraction = 1, .root_nonwoody_carbon_fraction = 1, .minimum_sink_ratio = 0.1, .base_exchange_rate_per_h = 0.01, .leaf_petiole_allocation_fraction = 1, .growth_habit = .perennial, .salt_mode = .dynamic, .salt_exchange_rate_per_h = 0.01, .biological_timestep_h = 1, .minimum_pool_g_c = 1e-12 };
    const c_before = branch_c[0] + branch_c[1] + root_c[0] + root_c[1];
    const n_before = branch_n[0] + branch_n[1] + root_n[0] + root_n[1];
    const p_before = branch_p[0] + branch_p[1] + root_p[0] + root_p[1];
    try apply(std.testing.allocator, state, workspace, inputs);
    try std.testing.expectApproxEqAbs(c_before, branch_c[0] + branch_c[1] + root_c[0] + root_c[1], 1e-13);
    try std.testing.expectApproxEqAbs(n_before, branch_n[0] + branch_n[1] + root_n[0] + root_n[1], 1e-13);
    try std.testing.expectApproxEqAbs(p_before, branch_p[0] + branch_p[1] + root_p[0] + root_p[1], 1e-13);
    for (0..standard_salt_species_count) |species| {
        const total = branch_salts[species * 2] + branch_salts[species * 2 + 1] + root_salts[species * 2] + root_salts[species * 2 + 1];
        try std.testing.expectApproxEqAbs(8.0, total, 1e-13);
    }
}

test "GROSUB shoot-root exchange rolls back an excessive transfer" {
    var branch_c = [_]f64{1};
    var branch_n = [_]f64{0.1};
    var branch_p = [_]f64{0.01};
    var root_c = [_]f64{0};
    var root_n = [_]f64{0};
    var root_p = [_]f64{0};
    var branch_salts = [_]f64{0} ** standard_salt_species_count;
    var root_salts = [_]f64{0} ** standard_salt_species_count;
    var total_leaf = [_]f64{1};
    var wc = [_]f64{0};
    var wn = [_]f64{0};
    var wp = [_]f64{0};
    var wrc = [_]f64{0};
    var wrn = [_]f64{0};
    var wrp = [_]f64{0};
    var wbs = [_]f64{0} ** standard_salt_species_count;
    var wrs = [_]f64{0} ** standard_salt_species_count;
    var wtl = [_]f64{0};
    const state: State = .{ .branch_mobile_carbon_g_c = &branch_c, .branch_mobile_nitrogen_g_n = &branch_n, .branch_mobile_phosphorus_g_p = &branch_p, .root_mobile_carbon_g_c_by_layer = &root_c, .root_mobile_nitrogen_g_n_by_layer = &root_n, .root_mobile_phosphorus_g_p_by_layer = &root_p, .branch_salt_mol_by_species_branch = &branch_salts, .root_salt_mol_by_species_layer = &root_salts, .total_leaf_petiole_carbon_g_c = &total_leaf };
    const workspace: State = .{ .branch_mobile_carbon_g_c = &wc, .branch_mobile_nitrogen_g_n = &wn, .branch_mobile_phosphorus_g_p = &wp, .root_mobile_carbon_g_c_by_layer = &wrc, .root_mobile_nitrogen_g_n_by_layer = &wrn, .root_mobile_phosphorus_g_p_by_layer = &wrp, .branch_salt_mol_by_species_branch = &wbs, .root_salt_mol_by_species_layer = &wrs, .total_leaf_petiole_carbon_g_c = &wtl };
    const inputs: Inputs = .{ .branch_count = 1, .soil_layer_count = 1, .shared_planting_layer_index = 0, .deepest_rooted_layer_index = 0, .branch_is_alive = &.{true}, .branch_leaf_petiole_carbon_g_c = &.{1}, .active_root_carbon_g_c_by_layer = &.{1}, .root_layer_sink_strength = &.{1}, .total_root_sink_strength = 1, .total_root_carbon_g_c = 1, .branch_nonwoody_carbon_fraction = 1, .root_nonwoody_carbon_fraction = 1, .minimum_sink_ratio = 0, .base_exchange_rate_per_h = 3, .leaf_petiole_allocation_fraction = 1, .growth_habit = .perennial, .salt_mode = .static_equilibrium, .salt_exchange_rate_per_h = 0, .biological_timestep_h = 1, .minimum_pool_g_c = 1e-12 };
    try std.testing.expectError(error.InvalidShootRootExchangeState, apply(std.testing.allocator, state, workspace, inputs));
    try std.testing.expectEqual(@as(f64, 1), branch_c[0]);
    try std.testing.expectEqual(@as(f64, 0), root_c[0]);
}
