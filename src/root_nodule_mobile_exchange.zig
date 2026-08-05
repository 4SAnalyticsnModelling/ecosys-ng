const std = @import("std");

pub const State = struct {
    root_mobile_carbon_g_c_by_layer: []f64,
    root_mobile_nitrogen_g_n_by_layer: []f64,
    root_mobile_phosphorus_g_p_by_layer: []f64,
    nodule_mobile_carbon_g_c_by_layer: []f64,
    nodule_mobile_nitrogen_g_n_by_layer: []f64,
    nodule_mobile_phosphorus_g_p_by_layer: []f64,
};

pub const Inputs = struct {
    soil_layer_count: usize,
    fixation_type: i32,
    /// Source `NU`, the first active soil layer (not planting layer `NG`).
    planting_layer_index: usize,
    plant_deepest_rooted_layer_index: usize,
    layer_active: []const bool,
    root_structural_carbon_g_c_by_layer: []const f64,
    nodule_structural_carbon_g_c_by_layer: []const f64,
    minimum_mobile_pool_g_c_by_layer: []const f64,
    minimum_root_structural_carbon_g_c_by_layer: []const f64,
    initial_nodule_carbon_g_c_m2: f64,
    root_zone_area_m2: f64,
    exchange_rate_h_inv: f64,
    biological_timestep_h: f64,
    site_type: i32,
    planting_day_marker: i32,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, layer_count: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != layer_count) return error.RootNoduleExchangeDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootNoduleExchangeState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.soil_layer_count == 0 or
        inputs.planting_layer_index > inputs.plant_deepest_rooted_layer_index or
        inputs.plant_deepest_rooted_layer_index >= inputs.soil_layer_count or
        inputs.layer_active.len != inputs.soil_layer_count)
        return error.RootNoduleExchangeDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != inputs.soil_layer_count) return error.RootNoduleExchangeDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootNoduleExchangeInput;
    };
    const scalars = [_]f64{
        inputs.initial_nodule_carbon_g_c_m2,
        inputs.root_zone_area_m2,
        inputs.exchange_rate_h_inv,
        inputs.biological_timestep_h,
    };
    for (scalars) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootNoduleExchangeInput;
}

/// Exact GROSUB 7940--7977 root--nodule non-structural C:N:P exchange.
/// Arrays are runtime-sized soil layers; masses are g C/N/P, area is m2,
/// rates are h-1, and the biological timestep is h. Carbon is committed to
/// the workspace before the N and P concentration differences are evaluated.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateState(state, inputs.soil_layer_count);
    try validateState(workspace, inputs.soil_layer_count);
    try validateInputs(inputs);
    copyState(workspace, state);

    const fixation_active = inputs.fixation_type >= 1 and inputs.fixation_type <= 3;
    const exchange_window_open = inputs.site_type != 0 or inputs.planting_day_marker == 0;
    if (fixation_active and exchange_window_open) for (inputs.planting_layer_index..inputs.plant_deepest_rooted_layer_index + 1) |layer| {
        if (!inputs.layer_active[layer] or
            workspace.root_mobile_carbon_g_c_by_layer[layer] <= inputs.minimum_mobile_pool_g_c_by_layer[layer] or
            inputs.root_structural_carbon_g_c_by_layer[layer] <= inputs.minimum_root_structural_carbon_g_c_by_layer[layer]) continue;

        const root_structural_carbon_g_c = inputs.root_structural_carbon_g_c_by_layer[layer];
        const nodule_structural_carbon_g_c = @min(
            root_structural_carbon_g_c,
            @max(
                inputs.initial_nodule_carbon_g_c_m2 * inputs.root_zone_area_m2,
                inputs.nodule_structural_carbon_g_c_by_layer[layer],
            ),
        );
        const combined_structural_carbon_g_c = root_structural_carbon_g_c + nodule_structural_carbon_g_c;
        if (combined_structural_carbon_g_c <= inputs.minimum_mobile_pool_g_c_by_layer[layer]) continue;

        const carbon_difference_g_c =
            (workspace.root_mobile_carbon_g_c_by_layer[layer] * nodule_structural_carbon_g_c -
                workspace.nodule_mobile_carbon_g_c_by_layer[layer] * root_structural_carbon_g_c) /
            combined_structural_carbon_g_c;
        const carbon_transfer_g_c = inputs.exchange_rate_h_inv * carbon_difference_g_c * inputs.biological_timestep_h;
        workspace.root_mobile_carbon_g_c_by_layer[layer] -= carbon_transfer_g_c;
        workspace.nodule_mobile_carbon_g_c_by_layer[layer] += carbon_transfer_g_c;

        const combined_mobile_carbon_g_c = workspace.root_mobile_carbon_g_c_by_layer[layer] +
            workspace.nodule_mobile_carbon_g_c_by_layer[layer];
        if (combined_mobile_carbon_g_c <= inputs.minimum_mobile_pool_g_c_by_layer[layer]) continue;

        const nitrogen_difference_g_n =
            (workspace.root_mobile_nitrogen_g_n_by_layer[layer] * workspace.nodule_mobile_carbon_g_c_by_layer[layer] -
                workspace.nodule_mobile_nitrogen_g_n_by_layer[layer] * workspace.root_mobile_carbon_g_c_by_layer[layer]) /
            combined_mobile_carbon_g_c;
        const nitrogen_transfer_g_n = inputs.exchange_rate_h_inv * nitrogen_difference_g_n * inputs.biological_timestep_h;
        const phosphorus_difference_g_p =
            (workspace.root_mobile_phosphorus_g_p_by_layer[layer] * workspace.nodule_mobile_carbon_g_c_by_layer[layer] -
                workspace.nodule_mobile_phosphorus_g_p_by_layer[layer] * workspace.root_mobile_carbon_g_c_by_layer[layer]) /
            combined_mobile_carbon_g_c;
        const phosphorus_transfer_g_p = inputs.exchange_rate_h_inv * phosphorus_difference_g_p * inputs.biological_timestep_h;
        workspace.root_mobile_nitrogen_g_n_by_layer[layer] -= nitrogen_transfer_g_n;
        workspace.root_mobile_phosphorus_g_p_by_layer[layer] -= phosphorus_transfer_g_p;
        workspace.nodule_mobile_nitrogen_g_n_by_layer[layer] += nitrogen_transfer_g_n;
        workspace.nodule_mobile_phosphorus_g_p_by_layer[layer] += phosphorus_transfer_g_p;
    };

    try validateState(workspace, inputs.soil_layer_count);
    copyState(state, workspace);
}

fn testState(carbon_root: []f64, nitrogen_root: []f64, phosphorus_root: []f64, carbon_nodule: []f64, nitrogen_nodule: []f64, phosphorus_nodule: []f64) State {
    return .{
        .root_mobile_carbon_g_c_by_layer = carbon_root,
        .root_mobile_nitrogen_g_n_by_layer = nitrogen_root,
        .root_mobile_phosphorus_g_p_by_layer = phosphorus_root,
        .nodule_mobile_carbon_g_c_by_layer = carbon_nodule,
        .nodule_mobile_nitrogen_g_n_by_layer = nitrogen_nodule,
        .nodule_mobile_phosphorus_g_p_by_layer = phosphorus_nodule,
    };
}

test "GROSUB root-nodule exchange conserves each mobile element" {
    var root_c = [_]f64{ 6.0, 8.0 };
    var root_n = [_]f64{ 0.6, 0.8 };
    var root_p = [_]f64{ 0.12, 0.16 };
    var nodule_c = [_]f64{ 2.0, 1.0 };
    var nodule_n = [_]f64{ 0.1, 0.05 };
    var nodule_p = [_]f64{ 0.02, 0.01 };
    var work_root_c = [_]f64{0} ** 2;
    var work_root_n = [_]f64{0} ** 2;
    var work_root_p = [_]f64{0} ** 2;
    var work_nodule_c = [_]f64{0} ** 2;
    var work_nodule_n = [_]f64{0} ** 2;
    var work_nodule_p = [_]f64{0} ** 2;
    const state = testState(&root_c, &root_n, &root_p, &nodule_c, &nodule_n, &nodule_p);
    const workspace = testState(&work_root_c, &work_root_n, &work_root_p, &work_nodule_c, &work_nodule_n, &work_nodule_p);
    var inputs: Inputs = .{
        .soil_layer_count = 2,
        .fixation_type = 2,
        .planting_layer_index = 0,
        .plant_deepest_rooted_layer_index = 1,
        .layer_active = &.{ true, true },
        .root_structural_carbon_g_c_by_layer = &.{ 10.0, 12.0 },
        .nodule_structural_carbon_g_c_by_layer = &.{ 2.0, 3.0 },
        .minimum_mobile_pool_g_c_by_layer = &.{ 1e-12, 1e-12 },
        .minimum_root_structural_carbon_g_c_by_layer = &.{ 1e-12, 1e-12 },
        .initial_nodule_carbon_g_c_m2 = 0.5,
        .root_zone_area_m2 = 1.0,
        .exchange_rate_h_inv = 0.1,
        .biological_timestep_h = 1.0,
        .site_type = 1,
        .planting_day_marker = 9,
    };
    const totals_before = [_]f64{ root_c[0] + nodule_c[0], root_n[0] + nodule_n[0], root_p[0] + nodule_p[0] };
    try apply(state, workspace, inputs);
    try std.testing.expectApproxEqAbs(totals_before[0], root_c[0] + nodule_c[0], 1e-14);
    try std.testing.expectApproxEqAbs(totals_before[1], root_n[0] + nodule_n[0], 1e-14);
    try std.testing.expectApproxEqAbs(totals_before[2], root_p[0] + nodule_p[0], 1e-14);
    try std.testing.expectApproxEqAbs(6.066666666666666, root_c[0], 1e-14);
    try std.testing.expectApproxEqAbs(1.9333333333333333, nodule_c[0], 1e-14);
    const after_fixing_exchange = root_c;
    inputs.fixation_type = 0;
    try apply(state, workspace, inputs);
    try std.testing.expectEqualDeep(after_fixing_exchange, root_c);
}

test "GROSUB root-nodule exchange is atomic when a transfer overdraws a pool" {
    var root_c = [_]f64{1.0};
    var root_n = [_]f64{0.1};
    var root_p = [_]f64{0.01};
    var nodule_c = [_]f64{0.0};
    var nodule_n = [_]f64{0.0};
    var nodule_p = [_]f64{0.0};
    var work_root_c = [_]f64{0.0};
    var work_root_n = [_]f64{0.0};
    var work_root_p = [_]f64{0.0};
    var work_nodule_c = [_]f64{0.0};
    var work_nodule_n = [_]f64{0.0};
    var work_nodule_p = [_]f64{0.0};
    const state = testState(&root_c, &root_n, &root_p, &nodule_c, &nodule_n, &nodule_p);
    const workspace = testState(&work_root_c, &work_root_n, &work_root_p, &work_nodule_c, &work_nodule_n, &work_nodule_p);
    const inputs: Inputs = .{
        .soil_layer_count = 1,
        .fixation_type = 2,
        .planting_layer_index = 0,
        .plant_deepest_rooted_layer_index = 0,
        .layer_active = &.{true},
        .root_structural_carbon_g_c_by_layer = &.{1.0},
        .nodule_structural_carbon_g_c_by_layer = &.{10.0},
        .minimum_mobile_pool_g_c_by_layer = &.{1e-12},
        .minimum_root_structural_carbon_g_c_by_layer = &.{1e-12},
        .initial_nodule_carbon_g_c_m2 = 0.0,
        .root_zone_area_m2 = 1.0,
        .exchange_rate_h_inv = 3.0,
        .biological_timestep_h = 1.0,
        .site_type = 1,
        .planting_day_marker = 1,
    };
    try std.testing.expectError(error.InvalidRootNoduleExchangeState, apply(state, workspace, inputs));
    try std.testing.expectEqual(@as(f64, 1.0), root_c[0]);
    try std.testing.expectEqual(@as(f64, 0.0), nodule_c[0]);
}
