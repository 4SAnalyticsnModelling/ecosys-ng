const std = @import("std");

pub const State = struct {
    carbon_recycling_fraction_by_layer: []f64,
    nitrogen_recycling_fraction_by_layer: []f64,
    phosphorus_recycling_fraction_by_layer: []f64,
    specific_decomposition_fraction_by_layer: []f64,
    decomposition_loss_carbon_g_c_by_layer: []f64,
    decomposition_loss_nitrogen_g_n_by_layer: []f64,
    decomposition_loss_phosphorus_g_p_by_layer: []f64,
    decomposition_recycled_carbon_g_c_by_layer: []f64,
    decomposition_recycled_nitrogen_g_n_by_layer: []f64,
    decomposition_recycled_phosphorus_g_p_by_layer: []f64,
    decomposition_litter_carbon_g_c_by_layer: []f64,
    decomposition_litter_nitrogen_g_n_by_layer: []f64,
    decomposition_litter_phosphorus_g_p_by_layer: []f64,
    growth_substrate_carbon_g_c_by_layer: []f64,
    structural_growth_carbon_g_c_by_layer: []f64,
    growth_respiration_g_c_by_layer: []f64,
    growth_nitrogen_g_n_by_layer: []f64,
    growth_phosphorus_g_p_by_layer: []f64,
    senescence_loss_carbon_g_c_by_layer: []f64,
    senescence_loss_nitrogen_g_n_by_layer: []f64,
    senescence_loss_phosphorus_g_p_by_layer: []f64,
    senescence_recycled_carbon_g_c_by_layer: []f64,
    senescence_recycled_nitrogen_g_n_by_layer: []f64,
    senescence_recycled_phosphorus_g_p_by_layer: []f64,
    senescence_litter_carbon_g_c_by_layer: []f64,
    senescence_litter_nitrogen_g_n_by_layer: []f64,
    senescence_litter_phosphorus_g_p_by_layer: []f64,
};

pub const Inputs = struct {
    soil_layer_count: usize,
    first_active_layer_index: usize,
    plant_deepest_rooted_layer_index: usize,
    fixation_type: i32,
    layer_active: []const bool,
    structural_carbon_g_c_by_layer: []const f64,
    structural_nitrogen_g_n_by_layer: []const f64,
    structural_phosphorus_g_p_by_layer: []const f64,
    mobile_carbon_g_c_by_layer: []const f64,
    mobile_nitrogen_g_n_by_layer: []const f64,
    mobile_phosphorus_g_p_by_layer: []const f64,
    mobile_nitrogen_concentration_by_layer: []const f64,
    mobile_phosphorus_concentration_by_layer: []const f64,
    carbon_balance_fraction_by_layer: []const f64,
    nitrogen_balance_fraction_by_layer: []const f64,
    phosphorus_balance_fraction_by_layer: []const f64,
    nodule_to_root_carbon_ratio_by_layer: []const f64,
    growth_temperature_response_by_layer: []const f64,
    growth_water_response_by_layer: []const f64,
    maintenance_respiration_g_c_by_layer: []const f64,
    substrate_respiration_g_c_by_layer: []const f64,
    growth_respiration_available_g_c_by_layer: []const f64,
    fixation_respiration_g_c_by_layer: []const f64,
    maintenance_deficit_g_c_by_layer: []const f64,
    minimum_carbon_recycling_fraction: f64,
    carbon_recycling_range_fraction: f64,
    maximum_nitrogen_recycling_fraction: f64,
    maximum_phosphorus_recycling_fraction: f64,
    specific_decomposition_rate_per_h: f64,
    decomposition_control_ratio: f64,
    biological_timestep_h: f64,
    growth_yield_g_c_per_g_c: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
    nitrogen_half_saturation_g_n_per_g_c: f64,
    phosphorus_half_saturation_g_p_per_g_c: f64,
    structural_presence_threshold_g_c: f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validate(state: State, inputs: Inputs) !void {
    if (inputs.soil_layer_count == 0 or inputs.first_active_layer_index > inputs.plant_deepest_rooted_layer_index or inputs.plant_deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.layer_active.len != inputs.soil_layer_count) return error.RootRhizobialRecyclingDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != inputs.soil_layer_count) return error.RootRhizobialRecyclingDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootRhizobialRecyclingState;
    }
    inline for (@typeInfo(Inputs).@"struct".fields) |field| switch (field.type) {
        []const f64 => {
            const values = @field(inputs, field.name);
            if (values.len != inputs.soil_layer_count) return error.RootRhizobialRecyclingDimensionMismatch;
            for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRhizobialRecyclingInput;
        },
        f64 => if (!std.math.isFinite(@field(inputs, field.name)) or @field(inputs, field.name) < 0) return error.InvalidRootRhizobialRecyclingInput,
        else => {},
    };
    if (inputs.decomposition_control_ratio == 0 or inputs.growth_yield_g_c_per_g_c >= 1 or inputs.nitrogen_half_saturation_g_n_per_g_c == 0 or inputs.phosphorus_half_saturation_g_p_per_g_c == 0) return error.InvalidRootRhizobialRecyclingInput;
}

/// Exact GROSUB 7758--7835 rhizobial recycling, decomposition, growth demand,
/// and maintenance-driven senescence. All outputs are per runtime soil layer;
/// masses are g C/N/P per biological timestep.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validate(state, inputs);
    try validate(workspace, inputs);
    copyState(workspace, state);
    if (inputs.fixation_type >= 1 and inputs.fixation_type <= 3) for (inputs.first_active_layer_index..inputs.plant_deepest_rooted_layer_index + 1) |layer| {
        if (!inputs.layer_active[layer]) continue;
        const recycle_c = inputs.minimum_carbon_recycling_fraction + inputs.carbon_balance_fraction_by_layer[layer] * inputs.carbon_recycling_range_fraction;
        const recycle_n = inputs.nitrogen_balance_fraction_by_layer[layer] * inputs.maximum_nitrogen_recycling_fraction;
        const recycle_p = inputs.phosphorus_balance_fraction_by_layer[layer] * inputs.maximum_phosphorus_recycling_fraction;
        const decomposition_fraction = @min(1, inputs.specific_decomposition_rate_per_h * inputs.nodule_to_root_carbon_ratio_by_layer[layer] / inputs.decomposition_control_ratio) * @sqrt(inputs.growth_temperature_response_by_layer[layer] * inputs.growth_water_response_by_layer[layer]) * inputs.biological_timestep_h;
        const loss_c = decomposition_fraction * inputs.structural_carbon_g_c_by_layer[layer];
        const loss_n = decomposition_fraction * inputs.structural_nitrogen_g_n_by_layer[layer];
        const loss_p = decomposition_fraction * inputs.structural_phosphorus_g_p_by_layer[layer];
        const recycled_c = loss_c * recycle_c;
        const recycled_n = loss_n * (recycle_n + (1 - recycle_n) * recycle_c);
        const recycled_p = loss_p * (recycle_p + (1 - recycle_p) * recycle_c);
        workspace.carbon_recycling_fraction_by_layer[layer] = recycle_c;
        workspace.nitrogen_recycling_fraction_by_layer[layer] = recycle_n;
        workspace.phosphorus_recycling_fraction_by_layer[layer] = recycle_p;
        workspace.specific_decomposition_fraction_by_layer[layer] = decomposition_fraction;
        workspace.decomposition_loss_carbon_g_c_by_layer[layer] = loss_c;
        workspace.decomposition_loss_nitrogen_g_n_by_layer[layer] = loss_n;
        workspace.decomposition_loss_phosphorus_g_p_by_layer[layer] = loss_p;
        workspace.decomposition_recycled_carbon_g_c_by_layer[layer] = recycled_c;
        workspace.decomposition_recycled_nitrogen_g_n_by_layer[layer] = recycled_n;
        workspace.decomposition_recycled_phosphorus_g_p_by_layer[layer] = recycled_p;
        workspace.decomposition_litter_carbon_g_c_by_layer[layer] = loss_c - recycled_c;
        workspace.decomposition_litter_nitrogen_g_n_by_layer[layer] = loss_n - recycled_n;
        workspace.decomposition_litter_phosphorus_g_p_by_layer[layer] = loss_p - recycled_p;

        const growth_substrate = @min(inputs.mobile_carbon_g_c_by_layer[layer] - @min(inputs.maintenance_respiration_g_c_by_layer[layer], inputs.substrate_respiration_g_c_by_layer[layer]) - inputs.fixation_respiration_g_c_by_layer[layer] + recycled_c, (inputs.growth_respiration_available_g_c_by_layer[layer] - inputs.fixation_respiration_g_c_by_layer[layer]) / (1 - inputs.growth_yield_g_c_per_g_c));
        const structural_growth = growth_substrate * inputs.growth_yield_g_c_per_g_c;
        workspace.growth_substrate_carbon_g_c_by_layer[layer] = growth_substrate;
        workspace.structural_growth_carbon_g_c_by_layer[layer] = structural_growth;
        workspace.growth_respiration_g_c_by_layer[layer] = inputs.fixation_respiration_g_c_by_layer[layer] + growth_substrate * (1 - inputs.growth_yield_g_c_per_g_c);
        workspace.growth_nitrogen_g_n_by_layer[layer] = @max(0, @min(inputs.mobile_nitrogen_g_n_by_layer[layer], structural_growth * inputs.target_nitrogen_per_carbon_g_n_per_g_c)) * inputs.mobile_nitrogen_concentration_by_layer[layer] / (inputs.mobile_nitrogen_concentration_by_layer[layer] + inputs.nitrogen_half_saturation_g_n_per_g_c);
        workspace.growth_phosphorus_g_p_by_layer[layer] = @max(0, @min(inputs.mobile_phosphorus_g_p_by_layer[layer], structural_growth * inputs.target_phosphorus_per_carbon_g_p_per_g_c)) * inputs.mobile_phosphorus_concentration_by_layer[layer] / (inputs.mobile_phosphorus_concentration_by_layer[layer] + inputs.phosphorus_half_saturation_g_p_per_g_c);

        if (inputs.maintenance_deficit_g_c_by_layer[layer] > 0 and inputs.structural_carbon_g_c_by_layer[layer] > inputs.structural_presence_threshold_g_c) {
            const senescence_c = inputs.maintenance_deficit_g_c_by_layer[layer];
            const senescence_n = senescence_c * inputs.structural_nitrogen_g_n_by_layer[layer] / inputs.structural_carbon_g_c_by_layer[layer];
            const senescence_p = senescence_c * inputs.structural_phosphorus_g_p_by_layer[layer] / inputs.structural_carbon_g_c_by_layer[layer];
            workspace.senescence_loss_carbon_g_c_by_layer[layer] = senescence_c;
            workspace.senescence_loss_nitrogen_g_n_by_layer[layer] = senescence_n;
            workspace.senescence_loss_phosphorus_g_p_by_layer[layer] = senescence_p;
            workspace.senescence_recycled_carbon_g_c_by_layer[layer] = senescence_c * recycle_c;
            workspace.senescence_recycled_nitrogen_g_n_by_layer[layer] = senescence_n * (recycle_n + (1 - recycle_n) * recycle_c);
            workspace.senescence_recycled_phosphorus_g_p_by_layer[layer] = senescence_p * (recycle_p + (1 - recycle_p) * recycle_c);
            workspace.senescence_litter_carbon_g_c_by_layer[layer] = senescence_c - workspace.senescence_recycled_carbon_g_c_by_layer[layer];
            workspace.senescence_litter_nitrogen_g_n_by_layer[layer] = senescence_n - workspace.senescence_recycled_nitrogen_g_n_by_layer[layer];
            workspace.senescence_litter_phosphorus_g_p_by_layer[layer] = senescence_p - workspace.senescence_recycled_phosphorus_g_p_by_layer[layer];
        } else {
            inline for (.{ workspace.senescence_loss_carbon_g_c_by_layer, workspace.senescence_loss_nitrogen_g_n_by_layer, workspace.senescence_loss_phosphorus_g_p_by_layer, workspace.senescence_recycled_carbon_g_c_by_layer, workspace.senescence_recycled_nitrogen_g_n_by_layer, workspace.senescence_recycled_phosphorus_g_p_by_layer, workspace.senescence_litter_carbon_g_c_by_layer, workspace.senescence_litter_nitrogen_g_n_by_layer, workspace.senescence_litter_phosphorus_g_p_by_layer }) |values| values[layer] = 0;
        }
    };
    inline for (@typeInfo(State).@"struct".fields) |field| for (@field(workspace, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootRhizobialRecyclingResult;
    copyState(state, workspace);
}

fn makeState(values: *[27][1]f64) State {
    var result: State = undefined;
    inline for (@typeInfo(State).@"struct".fields, 0..) |field, index| @field(result, field.name) = &values[index];
    return result;
}

fn testInputs() Inputs {
    const one = &[_]f64{1};
    return .{ .soil_layer_count = 1, .first_active_layer_index = 0, .plant_deepest_rooted_layer_index = 0, .fixation_type = 2, .layer_active = &.{true}, .structural_carbon_g_c_by_layer = &.{10}, .structural_nitrogen_g_n_by_layer = &.{1}, .structural_phosphorus_g_p_by_layer = &.{0.5}, .mobile_carbon_g_c_by_layer = &.{5}, .mobile_nitrogen_g_n_by_layer = &.{1}, .mobile_phosphorus_g_p_by_layer = &.{0.2}, .mobile_nitrogen_concentration_by_layer = one, .mobile_phosphorus_concentration_by_layer = one, .carbon_balance_fraction_by_layer = &.{0.5}, .nitrogen_balance_fraction_by_layer = &.{0.5}, .phosphorus_balance_fraction_by_layer = &.{0.5}, .nodule_to_root_carbon_ratio_by_layer = &.{0.2}, .growth_temperature_response_by_layer = one, .growth_water_response_by_layer = one, .maintenance_respiration_g_c_by_layer = &.{0.1}, .substrate_respiration_g_c_by_layer = &.{1}, .growth_respiration_available_g_c_by_layer = &.{2}, .fixation_respiration_g_c_by_layer = &.{0.2}, .maintenance_deficit_g_c_by_layer = &.{0.5}, .minimum_carbon_recycling_fraction = 0.1, .carbon_recycling_range_fraction = 0.4, .maximum_nitrogen_recycling_fraction = 0.8, .maximum_phosphorus_recycling_fraction = 0.6, .specific_decomposition_rate_per_h = 0.5, .decomposition_control_ratio = 0.1, .biological_timestep_h = 1, .growth_yield_g_c_per_g_c = 0.4, .target_nitrogen_per_carbon_g_n_per_g_c = 0.1, .target_phosphorus_per_carbon_g_p_per_g_c = 0.05, .nitrogen_half_saturation_g_n_per_g_c = 0.1, .phosphorus_half_saturation_g_p_per_g_c = 0.1, .structural_presence_threshold_g_c = 1e-12 };
}

test "GROSUB rhizobial decomposition and senescence partition C N P exactly" {
    var values: [27][1]f64 = std.mem.zeroes([27][1]f64);
    var work: [27][1]f64 = std.mem.zeroes([27][1]f64);
    try apply(makeState(&values), makeState(&work), testInputs());
    try std.testing.expectApproxEqAbs(values[4][0], values[7][0] + values[10][0], 1e-12);
    try std.testing.expectApproxEqAbs(values[5][0], values[8][0] + values[11][0], 1e-12);
    try std.testing.expectApproxEqAbs(values[6][0], values[9][0] + values[12][0], 1e-12);
    try std.testing.expectApproxEqAbs(values[18][0], values[21][0] + values[24][0], 1e-12);
    try std.testing.expectApproxEqAbs(values[19][0], values[22][0] + values[25][0], 1e-12);
    try std.testing.expectApproxEqAbs(values[20][0], values[23][0] + values[26][0], 1e-12);
}

test "GROSUB rhizobial recycling rejects singular yield atomically" {
    var values: [27][1]f64 = std.mem.zeroes([27][1]f64);
    const before = values;
    var work: [27][1]f64 = std.mem.zeroes([27][1]f64);
    var parameters = testInputs();
    parameters.growth_yield_g_c_per_g_c = 1;
    try std.testing.expectError(error.InvalidRootRhizobialRecyclingInput, apply(makeState(&values), makeState(&work), parameters));
    try std.testing.expectEqualDeep(before, values);
}
