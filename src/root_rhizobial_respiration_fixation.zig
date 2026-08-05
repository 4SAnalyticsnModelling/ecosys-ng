const std = @import("std");

pub const State = struct {
    substrate_respiration_oxygen_unlimited_g_c_by_layer: []f64,
    substrate_respiration_g_c_by_layer: []f64,
    maintenance_respiration_g_c_by_layer: []f64,
    net_growth_respiration_oxygen_unlimited_g_c_by_layer: []f64,
    net_growth_respiration_g_c_by_layer: []f64,
    growth_respiration_oxygen_unlimited_g_c_by_layer: []f64,
    growth_respiration_g_c_by_layer: []f64,
    maintenance_deficit_oxygen_unlimited_g_c_by_layer: []f64,
    maintenance_deficit_g_c_by_layer: []f64,
    fixation_respiration_requirement_g_c_by_layer: []f64,
    fixation_respiration_g_c_by_layer: []f64,
    nitrogen_fixation_g_n_by_layer: []f64,
    total_nitrogen_fixation_g_n: []f64,
};

pub const Inputs = struct {
    soil_layer_count: usize,
    first_active_layer_index: usize,
    plant_deepest_rooted_layer_index: usize,
    fixation_type: i32,
    layer_active: []const bool,
    mobile_carbon_g_c_by_layer: []const f64,
    structural_carbon_g_c_by_layer: []const f64,
    structural_nitrogen_g_n_by_layer: []const f64,
    nutrient_activity_fraction_by_layer: []const f64,
    growth_temperature_response_by_layer: []const f64,
    maintenance_temperature_response_by_layer: []const f64,
    growth_water_response_by_layer: []const f64,
    maintenance_water_response_by_layer: []const f64,
    oxygen_constraint_fraction_by_layer: []const f64,
    acidity_response_by_layer: []const f64,
    mobile_nitrogen_to_carbon_ratio_by_layer: []const f64,
    mobile_nitrogen_to_phosphorus_ratio_by_layer: []const f64,
    specific_respiration_per_h: f64,
    nitrogen_to_carbon_inhibition_ratio: f64,
    nitrogen_to_phosphorus_inhibition_ratio: f64,
    specific_maintenance_g_c_per_g_n_h: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    nitrogen_fixation_yield_g_n_per_g_c: f64,
    biological_timestep_h: f64,
    fixation_presence_threshold_g_c: f64,
};

fn validateState(state: State, layers: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (field.type == []f64 and values.len != (if (std.mem.eql(u8, field.name, "total_nitrogen_fixation_g_n")) 1 else layers)) return error.RootRhizobialRespirationDimensionMismatch;
        const signed = std.mem.startsWith(u8, field.name, "net_growth_respiration_");
        for (values) |value| if (!std.math.isFinite(value) or (!signed and value < 0)) return error.InvalidRootRhizobialRespirationState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.soil_layer_count == 0 or inputs.first_active_layer_index > inputs.plant_deepest_rooted_layer_index or inputs.plant_deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.layer_active.len != inputs.soil_layer_count) return error.RootRhizobialRespirationDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| switch (field.type) {
        []const f64 => {
            const values = @field(inputs, field.name);
            if (values.len != inputs.soil_layer_count) return error.RootRhizobialRespirationDimensionMismatch;
            for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRhizobialRespirationInput;
        },
        f64 => if (!std.math.isFinite(@field(inputs, field.name)) or @field(inputs, field.name) < 0) return error.InvalidRootRhizobialRespirationInput,
        else => {},
    };
    inline for (.{ inputs.nitrogen_to_carbon_inhibition_ratio, inputs.nitrogen_to_phosphorus_inhibition_ratio, inputs.nitrogen_fixation_yield_g_n_per_g_c }) |value| if (value == 0) return error.InvalidRootRhizobialRespirationInput;
    inline for (.{ inputs.nutrient_activity_fraction_by_layer, inputs.growth_water_response_by_layer, inputs.maintenance_water_response_by_layer, inputs.oxygen_constraint_fraction_by_layer, inputs.acidity_response_by_layer }) |values| for (values) |value| if (value > 1) return error.InvalidRootRhizobialRespirationInput;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

/// Exact GROSUB 7665--7730 rhizobial substrate, maintenance, growth and N2
/// fixation respiration. Traversal is ascending runtime root layer. Flux units
/// are g C or g N per biological timestep; rate inputs remain per hour.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateState(state, inputs.soil_layer_count);
    try validateState(workspace, inputs.soil_layer_count);
    try validateInputs(inputs);
    copyState(workspace, state);
    if (inputs.fixation_type >= 1 and inputs.fixation_type <= 3) for (inputs.first_active_layer_index..inputs.plant_deepest_rooted_layer_index + 1) |layer| {
        if (!inputs.layer_active[layer]) continue;
        const unlimited = @max(0, @min(inputs.mobile_carbon_g_c_by_layer[layer], inputs.specific_respiration_per_h * inputs.structural_carbon_g_c_by_layer[layer]) * inputs.nutrient_activity_fraction_by_layer[layer] * inputs.growth_temperature_response_by_layer[layer] / (1 + @max(inputs.mobile_nitrogen_to_carbon_ratio_by_layer[layer] / inputs.nitrogen_to_carbon_inhibition_ratio, inputs.mobile_nitrogen_to_phosphorus_ratio_by_layer[layer] / inputs.nitrogen_to_phosphorus_inhibition_ratio)) * inputs.growth_water_response_by_layer[layer] * inputs.biological_timestep_h);
        const limited = unlimited * inputs.oxygen_constraint_fraction_by_layer[layer];
        const maintenance = @max(0, inputs.specific_maintenance_g_c_per_g_n_h * inputs.maintenance_temperature_response_by_layer[layer] * inputs.structural_nitrogen_g_n_by_layer[layer] * inputs.maintenance_water_response_by_layer[layer] * inputs.acidity_response_by_layer[layer] * inputs.biological_timestep_h);
        const net_unlimited = unlimited - maintenance;
        const net_limited = limited - maintenance;
        const growth_unlimited = @max(0, net_unlimited);
        const growth_limited = @max(0, net_limited);
        workspace.substrate_respiration_oxygen_unlimited_g_c_by_layer[layer] = unlimited;
        workspace.substrate_respiration_g_c_by_layer[layer] = limited;
        workspace.maintenance_respiration_g_c_by_layer[layer] = maintenance;
        workspace.net_growth_respiration_oxygen_unlimited_g_c_by_layer[layer] = net_unlimited;
        workspace.net_growth_respiration_g_c_by_layer[layer] = net_limited;
        workspace.growth_respiration_oxygen_unlimited_g_c_by_layer[layer] = growth_unlimited;
        workspace.growth_respiration_g_c_by_layer[layer] = growth_limited;
        workspace.maintenance_deficit_oxygen_unlimited_g_c_by_layer[layer] = @max(0, -net_unlimited);
        workspace.maintenance_deficit_g_c_by_layer[layer] = @max(0, -net_limited);
        const fixation_requirement = @max(0, inputs.structural_carbon_g_c_by_layer[layer] * inputs.target_nitrogen_per_carbon_g_n_per_g_c - inputs.structural_nitrogen_g_n_by_layer[layer]) / inputs.nitrogen_fixation_yield_g_n_per_g_c;
        const fixation_respiration = if (growth_limited > inputs.fixation_presence_threshold_g_c) growth_limited * fixation_requirement / (growth_limited + fixation_requirement) else 0;
        workspace.fixation_respiration_requirement_g_c_by_layer[layer] = fixation_requirement;
        workspace.fixation_respiration_g_c_by_layer[layer] = fixation_respiration;
        workspace.nitrogen_fixation_g_n_by_layer[layer] = fixation_respiration * inputs.nitrogen_fixation_yield_g_n_per_g_c;
        workspace.total_nitrogen_fixation_g_n[0] += workspace.nitrogen_fixation_g_n_by_layer[layer];
    };
    // Net respiration differences are signed in Fortran and therefore are the
    // only state fields permitted below zero.
    inline for (.{ workspace.substrate_respiration_oxygen_unlimited_g_c_by_layer, workspace.substrate_respiration_g_c_by_layer, workspace.maintenance_respiration_g_c_by_layer, workspace.growth_respiration_oxygen_unlimited_g_c_by_layer, workspace.growth_respiration_g_c_by_layer, workspace.maintenance_deficit_oxygen_unlimited_g_c_by_layer, workspace.maintenance_deficit_g_c_by_layer, workspace.fixation_respiration_requirement_g_c_by_layer, workspace.fixation_respiration_g_c_by_layer, workspace.nitrogen_fixation_g_n_by_layer, workspace.total_nitrogen_fixation_g_n }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteRootRhizobialRespiration;
    inline for (.{ workspace.net_growth_respiration_oxygen_unlimited_g_c_by_layer, workspace.net_growth_respiration_g_c_by_layer }) |values| for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootRhizobialRespiration;
    copyState(state, workspace);
}

fn makeState(values: *[12][2]f64, total: *[1]f64) State {
    return .{ .substrate_respiration_oxygen_unlimited_g_c_by_layer = &values[0], .substrate_respiration_g_c_by_layer = &values[1], .maintenance_respiration_g_c_by_layer = &values[2], .net_growth_respiration_oxygen_unlimited_g_c_by_layer = &values[3], .net_growth_respiration_g_c_by_layer = &values[4], .growth_respiration_oxygen_unlimited_g_c_by_layer = &values[5], .growth_respiration_g_c_by_layer = &values[6], .maintenance_deficit_oxygen_unlimited_g_c_by_layer = &values[7], .maintenance_deficit_g_c_by_layer = &values[8], .fixation_respiration_requirement_g_c_by_layer = &values[9], .fixation_respiration_g_c_by_layer = &values[10], .nitrogen_fixation_g_n_by_layer = &values[11], .total_nitrogen_fixation_g_n = total };
}

fn makeInputs(oxygen: []const f64) Inputs {
    return .{ .soil_layer_count = 2, .first_active_layer_index = 0, .plant_deepest_rooted_layer_index = 1, .fixation_type = 2, .layer_active = &.{ true, true }, .mobile_carbon_g_c_by_layer = &.{ 10, 8 }, .structural_carbon_g_c_by_layer = &.{ 5, 4 }, .structural_nitrogen_g_n_by_layer = &.{ 0.2, 0.1 }, .nutrient_activity_fraction_by_layer = &.{ 1, 1 }, .growth_temperature_response_by_layer = &.{ 1, 1 }, .maintenance_temperature_response_by_layer = &.{ 1, 1 }, .growth_water_response_by_layer = &.{ 1, 1 }, .maintenance_water_response_by_layer = &.{ 1, 1 }, .oxygen_constraint_fraction_by_layer = oxygen, .acidity_response_by_layer = &.{ 1, 1 }, .mobile_nitrogen_to_carbon_ratio_by_layer = &.{ 0.1, 0.1 }, .mobile_nitrogen_to_phosphorus_ratio_by_layer = &.{ 1, 1 }, .specific_respiration_per_h = 1, .nitrogen_to_carbon_inhibition_ratio = 1, .nitrogen_to_phosphorus_inhibition_ratio = 1, .specific_maintenance_g_c_per_g_n_h = 0.5, .target_nitrogen_per_carbon_g_n_per_g_c = 0.1, .nitrogen_fixation_yield_g_n_per_g_c = 0.2, .biological_timestep_h = 1, .fixation_presence_threshold_g_c = 1e-12 };
}

test "GROSUB rhizobial respiration separates oxygen limits maintenance and fixation" {
    var values: [12][2]f64 = std.mem.zeroes([12][2]f64);
    var total = [_]f64{0};
    var work: [12][2]f64 = std.mem.zeroes([12][2]f64);
    var wt = [_]f64{0};
    try apply(makeState(&values, &total), makeState(&work, &wt), makeInputs(&.{ 0.5, 1 }));
    try std.testing.expectApproxEqAbs(values[0][0] * 0.5, values[1][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.1, values[2][0], 1e-12);
    try std.testing.expectApproxEqAbs(@max(0, values[1][0] - 0.1), values[6][0], 1e-12);
    try std.testing.expectApproxEqAbs(values[11][0] + values[11][1], total[0], 1e-12);
}

test "GROSUB rhizobial respiration rejects invalid oxygen atomically" {
    var values: [12][2]f64 = std.mem.zeroes([12][2]f64);
    const before = values;
    var total = [_]f64{0};
    var work: [12][2]f64 = std.mem.zeroes([12][2]f64);
    var wt = [_]f64{0};
    try std.testing.expectError(error.InvalidRootRhizobialRespirationInput, apply(makeState(&values, &total), makeState(&work, &wt), makeInputs(&.{ 0.5, 1.1 })));
    try std.testing.expectEqualDeep(before, values);
}
