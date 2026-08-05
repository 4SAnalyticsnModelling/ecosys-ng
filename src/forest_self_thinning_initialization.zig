const std = @import("std");

pub const State = struct {
    harvest_height_or_leaf_removal_by_species: []f64,
    preceding_harvest_height_or_leaf_removal_by_species: []f64,
    senesced_carbon_g_c_by_species_material_pool: []f64,
    senesced_nitrogen_g_n_by_species_material_pool: []f64,
    senesced_phosphorus_g_p_by_species_material_pool: []f64,
    senesced_inorganic_nitrogen_g_n_by_species: []f64,
    senesced_inorganic_phosphorus_g_p_by_species: []f64,
    harvest_type_by_species: []i32,
    termination_reseeding_code_by_species: []i32,
    thinning_fraction_by_species: []f64,
    harvest_efficiency_by_species_destination_material_pool: []f64,
};

pub const Inputs = struct {
    plant_species_count: usize,
    harvest_material_pool_count: usize,
    hourly_substep_index: usize,
    simulation_day_index: i32,
    hour_index: i32,
    local_solar_noon_h: f64,
    branch_type_by_species: []const i32,
    growth_type_by_species: []const i32,
    stem_diameter_m_by_species: []const f64,
    living_population_m2_by_species: []const f64,
};

fn poolIndex(inputs: Inputs, species: usize, pool: usize) usize {
    return species * inputs.harvest_material_pool_count + pool;
}

fn efficiencyIndex(inputs: Inputs, species: usize, destination: usize, pool: usize) usize {
    return (species * 2 + destination) * inputs.harvest_material_pool_count + pool;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateInputs(inputs: Inputs) !void {
    // Source EHVST assignments define four scientific material categories.
    // Runtime layouts may append categories, but cannot omit these four.
    if (inputs.plant_species_count == 0 or inputs.harvest_material_pool_count < 4) return error.ForestThinningDimensionMismatch;
    inline for (.{ inputs.branch_type_by_species, inputs.growth_type_by_species }) |values| if (values.len != inputs.plant_species_count) return error.ForestThinningDimensionMismatch;
    inline for (.{ inputs.stem_diameter_m_by_species, inputs.living_population_m2_by_species }) |values| {
        if (values.len != inputs.plant_species_count) return error.ForestThinningDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidForestThinningInput;
    }
    if (!std.math.isFinite(inputs.local_solar_noon_h)) return error.InvalidForestThinningInput;
}

fn validateState(state: State, inputs: Inputs) !void {
    const pool_count = std.math.mul(usize, inputs.plant_species_count, inputs.harvest_material_pool_count) catch return error.ForestThinningDimensionOverflow;
    const efficiency_count = std.math.mul(usize, pool_count, 2) catch return error.ForestThinningDimensionOverflow;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const expected = if (std.mem.indexOf(u8, field.name, "material_pool") != null) (if (std.mem.startsWith(u8, field.name, "harvest_efficiency")) efficiency_count else pool_count) else inputs.plant_species_count;
        if (values.len != expected) return error.ForestThinningDimensionMismatch;
        if (field.type == []f64) {
            const signed = std.mem.indexOf(u8, field.name, "harvest_height_or_leaf_removal") != null;
            for (values) |value| if (!std.math.isFinite(value) or (!signed and value < 0.0)) return error.InvalidForestThinningState;
        }
    }
}

fn thinningTrigger(inputs: Inputs, state: State, species: usize) bool {
    const harvest_type = state.harvest_type_by_species[species];
    return @mod(inputs.simulation_day_index, 30) == 0 and
        inputs.hour_index == @as(i32, @intFromFloat(@trunc(inputs.local_solar_noon_h))) and
        inputs.branch_type_by_species[species] != 0 and
        inputs.growth_type_by_species[species] > 1 and
        (harvest_type < 0 or harvest_type == 4 or harvest_type == 6);
}

/// Exact GROSUB 8516--8554 first-substep disturbance reset and forest
/// self-thinning. Runtime topology is plant species x harvest material pools.
/// Population uses plants m-2, stem diameter m, and thinning is a fraction.
/// With four material pools, efficiency assignments exactly match EHVST(1:2,1:4).
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    copyState(workspace, state);
    if (inputs.hourly_substep_index != 1) return;
    for (0..inputs.plant_species_count) |species| {
        workspace.preceding_harvest_height_or_leaf_removal_by_species[species] = workspace.harvest_height_or_leaf_removal_by_species[species];
        for (0..inputs.harvest_material_pool_count) |pool| {
            workspace.senesced_carbon_g_c_by_species_material_pool[poolIndex(inputs, species, pool)] = 0.0;
            workspace.senesced_nitrogen_g_n_by_species_material_pool[poolIndex(inputs, species, pool)] = 0.0;
            workspace.senesced_phosphorus_g_p_by_species_material_pool[poolIndex(inputs, species, pool)] = 0.0;
        }
        workspace.senesced_inorganic_nitrogen_g_n_by_species[species] = 0.0;
        workspace.senesced_inorganic_phosphorus_g_p_by_species[species] = 0.0;
        if (!thinningTrigger(inputs, workspace, species)) continue;

        const stem_diameter_m = inputs.stem_diameter_m_by_species[species];
        const living_population_m2 = inputs.living_population_m2_by_species[species];
        const thinning_fraction = if (stem_diameter_m > 0.0 and living_population_m2 > 0.0) blk: {
            const target_population_m2 = 0.1 * std.math.pow(f64, stem_diameter_m / 0.25, -1.6);
            break :blk @max(0.0, 1.0e-1 * (living_population_m2 - target_population_m2) / living_population_m2);
        } else 0.0;
        workspace.harvest_type_by_species[species] = 0;
        workspace.termination_reseeding_code_by_species[species] = 0;
        workspace.harvest_height_or_leaf_removal_by_species[species] = 1000.0;
        workspace.thinning_fraction_by_species[species] = thinning_fraction;
        for (0..inputs.harvest_material_pool_count) |pool| {
            workspace.harvest_efficiency_by_species_destination_material_pool[efficiencyIndex(inputs, species, 0, pool)] = if (pool < 3) 1.0 else 0.0;
            workspace.harvest_efficiency_by_species_destination_material_pool[efficiencyIndex(inputs, species, 1, pool)] = 0.0;
        }
    }
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

fn makeState(harvest: []f64, preceding: []f64, senesced_c: []f64, senesced_n: []f64, senesced_p: []f64, inorganic_n: []f64, inorganic_p: []f64, harvest_type: []i32, termination: []i32, thinning: []f64, efficiency: []f64) State {
    return .{ .harvest_height_or_leaf_removal_by_species = harvest, .preceding_harvest_height_or_leaf_removal_by_species = preceding, .senesced_carbon_g_c_by_species_material_pool = senesced_c, .senesced_nitrogen_g_n_by_species_material_pool = senesced_n, .senesced_phosphorus_g_p_by_species_material_pool = senesced_p, .senesced_inorganic_nitrogen_g_n_by_species = inorganic_n, .senesced_inorganic_phosphorus_g_p_by_species = inorganic_p, .harvest_type_by_species = harvest_type, .termination_reseeding_code_by_species = termination, .thinning_fraction_by_species = thinning, .harvest_efficiency_by_species_destination_material_pool = efficiency };
}

test "GROSUB hourly reset and forest self-thinning preserve assignment order" {
    var harvest = [_]f64{2};
    var preceding = [_]f64{9};
    var senesced_c = [_]f64{1} ** 4;
    var senesced_n = [_]f64{2} ** 4;
    var senesced_p = [_]f64{3} ** 4;
    var inorganic_n = [_]f64{4};
    var inorganic_p = [_]f64{5};
    var harvest_type = [_]i32{4};
    var termination = [_]i32{2};
    var thinning = [_]f64{0};
    var efficiency = [_]f64{9} ** 8;
    var wh = [_]f64{0};
    var wp = [_]f64{0};
    var wc = [_]f64{0} ** 4;
    var wn = [_]f64{0} ** 4;
    var wph = [_]f64{0} ** 4;
    var win = [_]f64{0};
    var wip = [_]f64{0};
    var wht = [_]i32{0};
    var wt = [_]i32{0};
    var wthin = [_]f64{0};
    var we = [_]f64{0} ** 8;
    const state = makeState(&harvest, &preceding, &senesced_c, &senesced_n, &senesced_p, &inorganic_n, &inorganic_p, &harvest_type, &termination, &thinning, &efficiency);
    const workspace = makeState(&wh, &wp, &wc, &wn, &wph, &win, &wip, &wht, &wt, &wthin, &we);
    const inputs: Inputs = .{ .plant_species_count = 1, .harvest_material_pool_count = 4, .hourly_substep_index = 1, .simulation_day_index = 30, .hour_index = 12, .local_solar_noon_h = 12.8, .branch_type_by_species = &.{1}, .growth_type_by_species = &.{2}, .stem_diameter_m_by_species = &.{0.25}, .living_population_m2_by_species = &.{1.0} };
    try apply(state, workspace, inputs);
    try std.testing.expectEqual(@as(f64, 2), preceding[0]);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0, 0 }, &senesced_c);
    try std.testing.expectApproxEqAbs(0.09, thinning[0], 1e-14);
    try std.testing.expectEqual(@as(f64, 1000), harvest[0]);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 1, 0, 0, 0, 0, 0 }, &efficiency);
}
