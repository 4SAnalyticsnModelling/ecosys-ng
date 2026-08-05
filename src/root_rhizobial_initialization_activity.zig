const std = @import("std");

pub const State = struct {
    structural_carbon_g_c_by_layer: []f64,
    structural_nitrogen_g_n_by_layer: []f64,
    structural_phosphorus_g_p_by_layer: []f64,
    mobile_carbon_g_c_by_layer: []f64,
    mobile_nitrogen_g_n_by_layer: []f64,
    mobile_phosphorus_g_p_by_layer: []f64,
    mobile_carbon_concentration_g_c_per_g_c_by_layer: []f64,
    mobile_nitrogen_concentration_g_n_per_g_c_by_layer: []f64,
    mobile_phosphorus_concentration_g_p_per_g_c_by_layer: []f64,
    mobile_nitrogen_to_carbon_ratio_by_layer: []f64,
    mobile_nitrogen_to_phosphorus_ratio_by_layer: []f64,
    carbon_recycling_fraction_by_layer: []f64,
    nitrogen_recycling_fraction_by_layer: []f64,
    phosphorus_recycling_fraction_by_layer: []f64,
    nutrient_activity_fraction_by_layer: []f64,
    nodule_to_root_carbon_ratio_by_layer: []f64,
};

pub const Inputs = struct {
    soil_layer_count: usize,
    first_active_layer_index: usize,
    plant_deepest_rooted_layer_index: usize,
    fixation_type: i32,
    first_subhour: bool,
    fire_active: bool,
    host_root_carbon_g_c_by_layer: []const f64,
    structural_presence_threshold_g_c: f64,
    host_root_presence_threshold_g_c: f64,
    initial_bacterial_carbon_g_c_per_m2: f64,
    first_active_layer_area_m2: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
    nitrogen_inhibition_g_n_per_g_c: f64,
    phosphorus_inhibition_g_p_per_g_c: f64,
};

fn validateState(state: State, layers: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != layers) return error.RootRhizobialActivityDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRhizobialActivityState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.soil_layer_count == 0 or inputs.first_active_layer_index > inputs.plant_deepest_rooted_layer_index or inputs.plant_deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.host_root_carbon_g_c_by_layer.len != inputs.soil_layer_count) return error.RootRhizobialActivityDimensionMismatch;
    inline for (.{ inputs.structural_presence_threshold_g_c, inputs.host_root_presence_threshold_g_c, inputs.initial_bacterial_carbon_g_c_per_m2, inputs.first_active_layer_area_m2, inputs.target_nitrogen_per_carbon_g_n_per_g_c, inputs.target_phosphorus_per_carbon_g_p_per_g_c, inputs.nitrogen_inhibition_g_n_per_g_c, inputs.phosphorus_inhibition_g_p_per_g_c }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRhizobialActivityInput;
    if (inputs.target_nitrogen_per_carbon_g_n_per_g_c == 0 or inputs.target_phosphorus_per_carbon_g_p_per_g_c == 0 or inputs.nitrogen_inhibition_g_n_per_g_c == 0 or inputs.phosphorus_inhibition_g_p_per_g_c == 0) return error.InvalidRootRhizobialActivityInput;
    for (inputs.host_root_carbon_g_c_by_layer) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRhizobialActivityInput;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

/// Exact GROSUB 7560--7647 rhizobial infection and C:N:P activity-state
/// construction. Pools are plant-wide by runtime soil layer. Masses use g C,
/// g N, and g P; infection density is g C m-2 and area is m2.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateState(state, inputs.soil_layer_count);
    try validateState(workspace, inputs.soil_layer_count);
    try validateInputs(inputs);
    copyState(workspace, state);
    if (inputs.fixation_type >= 1 and inputs.fixation_type <= 3) for (inputs.first_active_layer_index..inputs.plant_deepest_rooted_layer_index + 1) |layer| {
        const host_root_g_c = inputs.host_root_carbon_g_c_by_layer[layer];
        if (host_root_g_c <= inputs.host_root_presence_threshold_g_c) continue;
        if (inputs.first_subhour and workspace.structural_carbon_g_c_by_layer[layer] <= 0 and !inputs.fire_active) {
            const infected_carbon_g_c = inputs.initial_bacterial_carbon_g_c_per_m2 * inputs.first_active_layer_area_m2;
            workspace.structural_carbon_g_c_by_layer[layer] += infected_carbon_g_c;
            workspace.structural_nitrogen_g_n_by_layer[layer] += infected_carbon_g_c * inputs.target_nitrogen_per_carbon_g_n_per_g_c;
            workspace.structural_phosphorus_g_p_by_layer[layer] += infected_carbon_g_c * inputs.target_phosphorus_per_carbon_g_p_per_g_c;
        }
        const structural_g_c = workspace.structural_carbon_g_c_by_layer[layer];
        const mobile_c = if (structural_g_c > inputs.structural_presence_threshold_g_c) @max(0, workspace.mobile_carbon_g_c_by_layer[layer] / structural_g_c) else 1;
        const mobile_n = if (structural_g_c > inputs.structural_presence_threshold_g_c) @max(0, workspace.mobile_nitrogen_g_n_by_layer[layer] / structural_g_c) else 1;
        const mobile_p = if (structural_g_c > inputs.structural_presence_threshold_g_c) @max(0, workspace.mobile_phosphorus_g_p_by_layer[layer] / structural_g_c) else 1;
        workspace.mobile_carbon_concentration_g_c_per_g_c_by_layer[layer] = mobile_c;
        workspace.mobile_nitrogen_concentration_g_n_per_g_c_by_layer[layer] = mobile_n;
        workspace.mobile_phosphorus_concentration_g_p_per_g_c_by_layer[layer] = mobile_p;
        workspace.mobile_nitrogen_to_carbon_ratio_by_layer[layer] = if (mobile_c > 0) mobile_n / mobile_c else 0;
        workspace.mobile_nitrogen_to_phosphorus_ratio_by_layer[layer] = if (mobile_p > 0) mobile_n / mobile_p else 0;
        if (mobile_c > 0) {
            workspace.carbon_recycling_fraction_by_layer[layer] = std.math.clamp(@min(mobile_n / (mobile_n + mobile_c * inputs.nitrogen_inhibition_g_n_per_g_c), mobile_p / (mobile_p + mobile_c * inputs.phosphorus_inhibition_g_p_per_g_c)), 0, 1);
            workspace.nitrogen_recycling_fraction_by_layer[layer] = std.math.clamp(mobile_c / (mobile_c + mobile_n / inputs.nitrogen_inhibition_g_n_per_g_c), 0, 1);
            workspace.phosphorus_recycling_fraction_by_layer[layer] = std.math.clamp(mobile_c / (mobile_c + mobile_p / inputs.phosphorus_inhibition_g_p_per_g_c), 0, 1);
        } else {
            workspace.carbon_recycling_fraction_by_layer[layer] = 1;
            workspace.nitrogen_recycling_fraction_by_layer[layer] = 0;
            workspace.phosphorus_recycling_fraction_by_layer[layer] = 0;
        }
        if (structural_g_c > inputs.structural_presence_threshold_g_c) {
            const nitrogen_factor = @min(1, @sqrt(workspace.structural_nitrogen_g_n_by_layer[layer] / (structural_g_c * inputs.target_nitrogen_per_carbon_g_n_per_g_c)));
            const phosphorus_factor = @min(1, @sqrt(workspace.structural_phosphorus_g_p_by_layer[layer] / (structural_g_c * inputs.target_phosphorus_per_carbon_g_p_per_g_c)));
            workspace.nutrient_activity_fraction_by_layer[layer] = @min(nitrogen_factor, phosphorus_factor);
        } else workspace.nutrient_activity_fraction_by_layer[layer] = 1;
        workspace.nodule_to_root_carbon_ratio_by_layer[layer] = if (host_root_g_c > inputs.host_root_presence_threshold_g_c) @max(0, structural_g_c / host_root_g_c) else 0;
    };
    try validateState(workspace, inputs.soil_layer_count);
    copyState(state, workspace);
}

fn makeState(values: *[16][2]f64) State {
    return .{ .structural_carbon_g_c_by_layer = &values[0], .structural_nitrogen_g_n_by_layer = &values[1], .structural_phosphorus_g_p_by_layer = &values[2], .mobile_carbon_g_c_by_layer = &values[3], .mobile_nitrogen_g_n_by_layer = &values[4], .mobile_phosphorus_g_p_by_layer = &values[5], .mobile_carbon_concentration_g_c_per_g_c_by_layer = &values[6], .mobile_nitrogen_concentration_g_n_per_g_c_by_layer = &values[7], .mobile_phosphorus_concentration_g_p_per_g_c_by_layer = &values[8], .mobile_nitrogen_to_carbon_ratio_by_layer = &values[9], .mobile_nitrogen_to_phosphorus_ratio_by_layer = &values[10], .carbon_recycling_fraction_by_layer = &values[11], .nitrogen_recycling_fraction_by_layer = &values[12], .phosphorus_recycling_fraction_by_layer = &values[13], .nutrient_activity_fraction_by_layer = &values[14], .nodule_to_root_carbon_ratio_by_layer = &values[15] };
}

test "GROSUB first-subhour rhizobial infection immediately constructs activity state" {
    var values: [16][2]f64 = std.mem.zeroes([16][2]f64);
    values[3] = .{ 2, 4 };
    values[4] = .{ 0.2, 0.4 };
    values[5] = .{ 0.1, 0.2 };
    var work: [16][2]f64 = std.mem.zeroes([16][2]f64);
    try apply(makeState(&values), makeState(&work), .{ .soil_layer_count = 2, .first_active_layer_index = 0, .plant_deepest_rooted_layer_index = 1, .fixation_type = 2, .first_subhour = true, .fire_active = false, .host_root_carbon_g_c_by_layer = &.{ 20, 10 }, .structural_presence_threshold_g_c = 1e-12, .host_root_presence_threshold_g_c = 1e-12, .initial_bacterial_carbon_g_c_per_m2 = 0.5, .first_active_layer_area_m2 = 2, .target_nitrogen_per_carbon_g_n_per_g_c = 0.1, .target_phosphorus_per_carbon_g_p_per_g_c = 0.05, .nitrogen_inhibition_g_n_per_g_c = 0.2, .phosphorus_inhibition_g_p_per_g_c = 0.1 });
    try std.testing.expectEqual(@as(f64, 1), values[0][0]);
    try std.testing.expectEqual(@as(f64, 1), values[0][1]);
    try std.testing.expectApproxEqAbs(1, values[14][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.05, values[15][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.1, values[9][0], 1e-12);
    try std.testing.expectApproxEqAbs(2, values[10][0], 1e-12);
}

test "GROSUB rhizobial activity skips unrooted layers and rejects zero inhibition atomically" {
    var values: [16][2]f64 = std.mem.zeroes([16][2]f64);
    values[0][1] = 2;
    const before = values;
    var work: [16][2]f64 = std.mem.zeroes([16][2]f64);
    try std.testing.expectError(error.InvalidRootRhizobialActivityInput, apply(makeState(&values), makeState(&work), .{ .soil_layer_count = 2, .first_active_layer_index = 0, .plant_deepest_rooted_layer_index = 1, .fixation_type = 1, .first_subhour = true, .fire_active = false, .host_root_carbon_g_c_by_layer = &.{ 0, 10 }, .structural_presence_threshold_g_c = 0, .host_root_presence_threshold_g_c = 0, .initial_bacterial_carbon_g_c_per_m2 = 1, .first_active_layer_area_m2 = 1, .target_nitrogen_per_carbon_g_n_per_g_c = 0.1, .target_phosphorus_per_carbon_g_p_per_g_c = 0.1, .nitrogen_inhibition_g_n_per_g_c = 0, .phosphorus_inhibition_g_p_per_g_c = 0.1 }));
    try std.testing.expectEqualDeep(before, values);
}
