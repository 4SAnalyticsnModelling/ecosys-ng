const std = @import("std");

pub const State = struct {
    living_population_m2_by_species: []f64,
    living_population_by_species: []f64,
    standing_dead_population_by_species: []f64,
    clumping_factor_by_species: []f64,
    harvest_height_m_or_leaf_area_fraction_by_species: []f64,
};

pub const Inputs = struct {
    plant_species_count: usize,
    canopy_layer_count: usize,
    harvest_type_by_species: []const i32,
    termination_reseeding_code_by_species: []const i32,
    thinning_fraction_by_species: []const f64,
    seeded_population_m2_by_species: []const f64,
    planting_layer_area_m2_by_species: []const f64,
    hour_index: i32,
    local_solar_noon_h: f64,
    total_canopy_leaf_area_m2: f64,
    canopy_layer_bottom_height_m: []const f64,
    canopy_leaf_area_m2_by_layer: []const f64,
    minimum_canopy_leaf_area_m2: f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.plant_species_count == 0 or inputs.canopy_layer_count == 0 or inputs.canopy_layer_bottom_height_m.len != inputs.canopy_layer_count + 1 or inputs.canopy_leaf_area_m2_by_layer.len != inputs.canopy_layer_count) return error.DisturbancePreparationDimensionMismatch;
    inline for (.{ inputs.harvest_type_by_species, inputs.termination_reseeding_code_by_species }) |values| if (values.len != inputs.plant_species_count) return error.DisturbancePreparationDimensionMismatch;
    inline for (.{ inputs.thinning_fraction_by_species, inputs.seeded_population_m2_by_species, inputs.planting_layer_area_m2_by_species }) |values| {
        if (values.len != inputs.plant_species_count) return error.DisturbancePreparationDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidDisturbancePreparationInput;
    }
    const scalars = [_]f64{ inputs.local_solar_noon_h, inputs.total_canopy_leaf_area_m2, inputs.minimum_canopy_leaf_area_m2 };
    for (scalars) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidDisturbancePreparationInput;
    for (inputs.canopy_layer_bottom_height_m) |value| if (!std.math.isFinite(value)) return error.InvalidDisturbancePreparationInput;
    for (inputs.canopy_leaf_area_m2_by_layer) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidDisturbancePreparationInput;
}

fn validateEventValues(state: State, inputs: Inputs) !void {
    for (inputs.harvest_type_by_species, 0..) |harvest_type, species| {
        const event_value = state.harvest_height_m_or_leaf_area_fraction_by_species[species];
        if (harvest_type <= 2 and event_value < -1.0) return error.InvalidDisturbancePreparationState;
        if ((harvest_type == 4 or harvest_type == 6) and event_value < 0.0) return error.InvalidDisturbancePreparationState;
        if (inputs.thinning_fraction_by_species[species] > 1.0) return error.InvalidDisturbancePreparationInput;
    }
}

fn validateState(state: State, species_count: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != species_count) return error.DisturbancePreparationDimensionMismatch;
        const signed = std.mem.startsWith(u8, field.name, "harvest_height");
        for (values) |value| if (!std.math.isFinite(value) or (!signed and value < 0.0)) return error.InvalidDisturbancePreparationState;
    }
}

fn eventEnabled(inputs: Inputs, harvest_type: i32) bool {
    return (harvest_type >= 0 and
        inputs.hour_index == @as(i32, @intFromFloat(@trunc(inputs.local_solar_noon_h))) and
        harvest_type != 4 and harvest_type != 6) or harvest_type == 4 or harvest_type == 6;
}

/// Exact GROSUB 8566--8625 above-ground disturbance gate, population update,
/// pruning clumping update, and fractional-LAI-to-cutting-height conversion.
/// Runtime topology is plant species x canopy layers. Population is plants m-2
/// or plants cell-1, leaf area m2, height m, and thinning is dimensionless.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs.plant_species_count);
    try validateState(workspace, inputs.plant_species_count);
    try validateEventValues(state, inputs);
    copyState(workspace, state);
    for (0..inputs.plant_species_count) |species| {
        const harvest_type = inputs.harvest_type_by_species[species];
        if (!eventEnabled(inputs, harvest_type) or harvest_type == 4 or harvest_type == 6) continue;
        const thinning_fraction = inputs.thinning_fraction_by_species[species];
        if (inputs.termination_reseeding_code_by_species[species] != 2) {
            workspace.living_population_m2_by_species[species] *= 1.0 - thinning_fraction;
            workspace.living_population_by_species[species] *= 1.0 - thinning_fraction;
            workspace.standing_dead_population_by_species[species] *= 1.0 - thinning_fraction;
        } else {
            workspace.living_population_m2_by_species[species] = inputs.seeded_population_m2_by_species[species];
            workspace.living_population_by_species[species] = workspace.living_population_m2_by_species[species] * inputs.planting_layer_area_m2_by_species[species];
            workspace.standing_dead_population_by_species[species] = workspace.living_population_m2_by_species[species] * inputs.planting_layer_area_m2_by_species[species];
        }
        if (harvest_type == 3) workspace.clumping_factor_by_species[species] *= workspace.harvest_height_m_or_leaf_area_fraction_by_species[species];
        if (harvest_type <= 2 and workspace.harvest_height_m_or_leaf_area_fraction_by_species[species] < 0.0) {
            const remaining_leaf_area_m2 = (1.0 - @abs(workspace.harvest_height_m_or_leaf_area_fraction_by_species[species])) * inputs.total_canopy_leaf_area_m2;
            var accumulated_leaf_area_m2: f64 = 0.0;
            for (0..inputs.canopy_layer_count) |layer| {
                const lower_height_m = inputs.canopy_layer_bottom_height_m[layer];
                const upper_height_m = inputs.canopy_layer_bottom_height_m[layer + 1];
                const layer_leaf_area_m2 = inputs.canopy_leaf_area_m2_by_layer[layer];
                if (upper_height_m > lower_height_m and layer_leaf_area_m2 > inputs.minimum_canopy_leaf_area_m2 and accumulated_leaf_area_m2 < remaining_leaf_area_m2) {
                    workspace.harvest_height_m_or_leaf_area_fraction_by_species[species] = if (accumulated_leaf_area_m2 + layer_leaf_area_m2 > remaining_leaf_area_m2)
                        lower_height_m + ((remaining_leaf_area_m2 - accumulated_leaf_area_m2) / layer_leaf_area_m2) * (upper_height_m - lower_height_m)
                    else
                        0.0;
                    accumulated_leaf_area_m2 += layer_leaf_area_m2;
                }
            }
        }
    }
    try validateState(workspace, inputs.plant_species_count);
    copyState(state, workspace);
}

fn makeState(population_m2: []f64, population: []f64, dead: []f64, clumping: []f64, harvest: []f64) State {
    return .{ .living_population_m2_by_species = population_m2, .living_population_by_species = population, .standing_dead_population_by_species = dead, .clumping_factor_by_species = clumping, .harvest_height_m_or_leaf_area_fraction_by_species = harvest };
}

test "GROSUB disturbance preparation scales population and resolves cutting height" {
    var population_m2 = [_]f64{ 10, 5 };
    var population = [_]f64{ 100, 50 };
    var dead = [_]f64{ 20, 10 };
    var clumping = [_]f64{ 0.8, 0.7 };
    var harvest = [_]f64{ -0.5, 2 };
    var wp2 = [_]f64{0} ** 2;
    var wp = [_]f64{0} ** 2;
    var wd = [_]f64{0} ** 2;
    var wc = [_]f64{0} ** 2;
    var wh = [_]f64{0} ** 2;
    const state = makeState(&population_m2, &population, &dead, &clumping, &harvest);
    const workspace = makeState(&wp2, &wp, &wd, &wc, &wh);
    const inputs: Inputs = .{ .plant_species_count = 2, .canopy_layer_count = 2, .harvest_type_by_species = &.{ 2, 2 }, .termination_reseeding_code_by_species = &.{ 0, 2 }, .thinning_fraction_by_species = &.{ 0.1, 0.2 }, .seeded_population_m2_by_species = &.{ 3, 4 }, .planting_layer_area_m2_by_species = &.{ 10, 10 }, .hour_index = 12, .local_solar_noon_h = 12.8, .total_canopy_leaf_area_m2 = 10, .canopy_layer_bottom_height_m = &.{ 0, 1, 2 }, .canopy_leaf_area_m2_by_layer = &.{ 4, 6 }, .minimum_canopy_leaf_area_m2 = 1e-12 };
    try apply(state, workspace, inputs);
    try std.testing.expectEqual(@as(f64, 9), population_m2[0]);
    try std.testing.expectEqual(@as(f64, 4), population_m2[1]);
    try std.testing.expectEqual(@as(f64, 40), population[1]);
    try std.testing.expectApproxEqAbs(1.1666666666666667, harvest[0], 1e-14);
}

test "GROSUB disturbance gate leaves noon-only harvest unchanged off hour" {
    var population_m2 = [_]f64{10};
    var population = [_]f64{100};
    var dead = [_]f64{20};
    var clumping = [_]f64{0.8};
    var harvest = [_]f64{1};
    var wp2 = [_]f64{0};
    var wp = [_]f64{0};
    var wd = [_]f64{0};
    var wc = [_]f64{0};
    var wh = [_]f64{0};
    const state = makeState(&population_m2, &population, &dead, &clumping, &harvest);
    const workspace = makeState(&wp2, &wp, &wd, &wc, &wh);
    const inputs: Inputs = .{ .plant_species_count = 1, .canopy_layer_count = 1, .harvest_type_by_species = &.{2}, .termination_reseeding_code_by_species = &.{0}, .thinning_fraction_by_species = &.{0.5}, .seeded_population_m2_by_species = &.{1}, .planting_layer_area_m2_by_species = &.{1}, .hour_index = 11, .local_solar_noon_h = 12, .total_canopy_leaf_area_m2 = 1, .canopy_layer_bottom_height_m = &.{ 0, 1 }, .canopy_leaf_area_m2_by_layer = &.{1}, .minimum_canopy_leaf_area_m2 = 0 };
    try apply(state, workspace, inputs);
    try std.testing.expectEqual(@as(f64, 10), population_m2[0]);
}
