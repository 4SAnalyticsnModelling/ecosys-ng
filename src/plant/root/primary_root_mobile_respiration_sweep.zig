const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    mobile_carbon_g_c_by_domain_layer: []f64,
    mobile_nitrogen_g_n_by_domain_layer: []f64,
    mobile_phosphorus_g_p_by_domain_layer: []f64,
    /// Source `RCO2A`; respiration is published as a negative C flux.
    net_root_carbon_flux_g_c_by_domain_layer: []f64,
    oxygen_unlimited_respiration_g_c_by_domain_layer: []f64,
    carbon_unlimited_respiration_g_c_by_domain_layer: []f64,
};

pub const Workspace = State;

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    /// Source `NI(NZ,NY,NX)`, shared across biological domains.
    deepest_rooted_layer_index: usize,
    tip_active_by_domain_layer_axis: []const bool,
    metabolism_by_domain_layer_axis: []const metabolism.SecondaryRootResult,
    senescence_by_domain_layer_axis: []const metabolism.SecondaryRootSenescence,
    nonwoody_carbon_fraction: f64,
    nonwoody_nitrogen_fraction: f64,
    root_depth_m_by_domain_axis: []const f64,
    /// Source `SDPTH(NZ,NY,NX)`.
    seeding_depth_m: f64,
    /// Source `NG(NZ,NY,NX)`.
    first_rooted_layer_index: usize,
    /// Source `NINR(NR,NZ,NY,NX)`, flattened by root axis only.
    deepest_axis_layer_index_by_axis: []const usize,
    layer_bottom_depth_m: []const f64,
    primary_root_length_m_by_domain_layer_axis: []const f64,
};

fn copyState(target: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(target, field.name), @field(source, field.name));
}

fn validateShape(state: State, domain_layers: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (@field(state, field.name).len != domain_layers) return error.PrimaryRootMobileRespirationDimensionMismatch;
}

fn validateState(state: State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| for (@field(state, field.name)) |value| {
        if (!std.math.isFinite(value)) return error.InvalidPrimaryRootMobileRespirationState;
        if (!std.mem.eql(u8, field.name, "net_root_carbon_flux_g_c_by_domain_layer") and value < 0)
            return error.InvalidPrimaryRootMobileRespirationState;
    };
}

fn updateTip(state: State, inputs: Inputs, domain: usize, layer: usize, axis: usize) !void {
    const domain_layer = domain * inputs.soil_layer_count + layer;
    const domain_axis = domain * inputs.root_axis_count + axis;
    const index = domain_layer * inputs.root_axis_count + axis;
    const process = inputs.metabolism_by_domain_layer_axis[index];
    const loss = inputs.senescence_by_domain_layer_axis[index];
    inline for (@typeInfo(metabolism.SecondaryRootResult).@"struct".fields) |field| if (!std.math.isFinite(@field(process, field.name)) or @field(process, field.name) < 0) return error.InvalidPrimaryRootMetabolismResult;
    inline for (@typeInfo(metabolism.SecondaryRootSenescence).@"struct".fields) |field| if (!std.math.isFinite(@field(loss, field.name)) or @field(loss, field.name) < 0) return error.InvalidPrimaryRootSenescenceResult;
    if (loss.senesced_fraction > 1) return error.InvalidPrimaryRootSenescenceResult;
    const recovered_c = loss.senesced_fraction * loss.recyclable_carbon_g_c * inputs.nonwoody_carbon_fraction;
    const recovered_n = loss.senesced_fraction * loss.recyclable_nitrogen_g_n * inputs.nonwoody_nitrogen_fraction;
    // Source line 6771 deliberately has no FWODRP(1) multiplier.
    const recovered_p = loss.senesced_fraction * loss.recyclable_phosphorus_g_p;
    state.mobile_carbon_g_c_by_domain_layer[domain_layer] += -@min(process.maintenance_respiration_g_c_per_h, process.substrate_respiration_actual_g_c_per_h) - process.growth_and_respiration_carbon_actual_g_c_per_h - process.nitrogen_assimilation_respiration_actual_g_c_per_h - loss.respiration_actual_g_c_per_h + recovered_c;
    state.mobile_nitrogen_g_n_by_domain_layer[domain_layer] += -process.nitrogen_growth_actual_g_n_per_h + recovered_n;
    state.mobile_phosphorus_g_p_by_domain_layer[domain_layer] += -process.phosphorus_growth_actual_g_p_per_h + recovered_p;
    const respiration = try metabolism.assemble(.{
        .maintenance_demand_g_c = process.maintenance_respiration_g_c_per_h,
        .substrate_respiration_actual_g_c = process.substrate_respiration_actual_g_c_per_h,
        .substrate_respiration_oxygen_unlimited_g_c = process.substrate_respiration_oxygen_unlimited_g_c_per_h,
        .growth_respiration_actual_g_c = process.growth_respiration_actual_g_c_per_h,
        .growth_respiration_oxygen_unlimited_g_c = process.growth_respiration_oxygen_unlimited_g_c_per_h,
        .senescence_respiration_actual_g_c = loss.respiration_actual_g_c_per_h,
        .senescence_respiration_oxygen_unlimited_g_c = loss.respiration_oxygen_unlimited_g_c_per_h,
        .nitrogen_assimilation_respiration_actual_g_c = process.nitrogen_assimilation_respiration_actual_g_c_per_h,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = process.nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h,
    });
    const first_layer = inputs.first_rooted_layer_index;
    const deepest_layer = inputs.deepest_axis_layer_index_by_axis[axis];
    if (inputs.root_depth_m_by_domain_axis[domain_axis] > inputs.layer_bottom_depth_m[first_layer]) {
        const denominator_m = inputs.root_depth_m_by_domain_axis[domain_axis] - inputs.seeding_depth_m;
        if (!(denominator_m > 0)) return error.InvalidPrimaryRootRespirationGeometry;
        var allocated: f64 = 0;
        for (first_layer..deepest_layer + 1) |respiration_layer| {
            const length_index = (domain * inputs.soil_layer_count + respiration_layer) * inputs.root_axis_count + axis;
            const fraction = if (respiration_layer < deepest_layer)
                @min(1, inputs.primary_root_length_m_by_domain_layer_axis[length_index] / denominator_m)
            else
                1 - allocated;
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidPrimaryRootRespirationFraction;
            allocated += fraction;
            const target = domain * inputs.soil_layer_count + respiration_layer;
            state.net_root_carbon_flux_g_c_by_domain_layer[target] -= respiration.actual_g_c * fraction;
            state.oxygen_unlimited_respiration_g_c_by_domain_layer[target] += respiration.oxygen_unlimited_g_c * fraction;
            state.carbon_unlimited_respiration_g_c_by_domain_layer[target] += respiration.carbon_unlimited_g_c * fraction;
        }
        if (@abs(allocated - 1) > 1e-12) return error.InvalidPrimaryRootRespirationFraction;
    } else {
        state.net_root_carbon_flux_g_c_by_domain_layer[domain_layer] -= respiration.actual_g_c;
        state.oxygen_unlimited_respiration_g_c_by_domain_layer[domain_layer] += respiration.oxygen_unlimited_g_c;
        state.carbon_unlimited_respiration_g_c_by_domain_layer[domain_layer] += respiration.carbon_unlimited_g_c;
    }
    try validateState(state);
}

/// Exact grosub.f lines 6766--6828 N,L,NR primary mobile-pool and nested LL
/// respiration transaction. C/N/P units are grams per biological hour;
/// depths and lengths are m and allocation fractions are dimensionless.
pub fn apply(state: State, workspace: Workspace, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootMobileRespirationDimensionOverflow;
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.PrimaryRootMobileRespirationDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootMobileRespirationDimensionOverflow;
    try validateShape(state, domain_layers);
    try validateShape(workspace, domain_layers);
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.tip_active_by_domain_layer_axis.len != values or
        inputs.metabolism_by_domain_layer_axis.len != values or inputs.senescence_by_domain_layer_axis.len != values or
        inputs.root_depth_m_by_domain_axis.len != domain_axes or
        inputs.deepest_axis_layer_index_by_axis.len != inputs.root_axis_count or
        inputs.layer_bottom_depth_m.len != inputs.soil_layer_count or inputs.primary_root_length_m_by_domain_layer_axis.len != values)
        return error.PrimaryRootMobileRespirationDimensionMismatch;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.first_rooted_layer_index >= inputs.soil_layer_count) return error.InvalidPrimaryRootMobileRespirationLayerRange;
    for (inputs.deepest_axis_layer_index_by_axis) |deepest| if (deepest < inputs.first_rooted_layer_index or deepest >= inputs.soil_layer_count) return error.InvalidPrimaryRootMobileRespirationLayerRange;
    inline for (.{ inputs.nonwoody_carbon_fraction, inputs.nonwoody_nitrogen_fraction }) |value| if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidPrimaryRootMobileRespirationInput;
    if (!std.math.isFinite(inputs.seeding_depth_m) or inputs.seeding_depth_m < 0) return error.InvalidPrimaryRootMobileRespirationInput;
    inline for (.{ inputs.root_depth_m_by_domain_axis, inputs.layer_bottom_depth_m, inputs.primary_root_length_m_by_domain_layer_axis }) |slice| for (slice) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootMobileRespirationInput;
    for (1..inputs.biological_domain_count) |domain| for (0..inputs.soil_layer_count) |layer| for (0..inputs.root_axis_count) |axis| if (inputs.tip_active_by_domain_layer_axis[(domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis]) return error.PrimaryRootTipOutsideHostDomain;
    try validateState(state);
    copyState(workspace, state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| for (0..inputs.root_axis_count) |axis| {
        const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
        if (inputs.tip_active_by_domain_layer_axis[index]) try updateTip(workspace, inputs, domain, layer, axis);
    };
    copyState(state, workspace);
}

fn testProcess() metabolism.SecondaryRootResult {
    return .{ .nutrient_feedback = 1, .substrate_respiration_oxygen_unlimited_g_c_per_h = 0.6, .maintenance_respiration_g_c_per_h = 0.1, .substrate_respiration_actual_g_c_per_h = 0.5, .growth_respiration_oxygen_unlimited_g_c_per_h = 0.25, .growth_respiration_actual_g_c_per_h = 0.2, .growth_and_respiration_carbon_oxygen_unlimited_g_c_per_h = 1, .growth_and_respiration_carbon_actual_g_c_per_h = 0.8, .root_growth_oxygen_unlimited_g_c_per_h = 0.75, .root_growth_actual_g_c_per_h = 0.6, .nitrogen_growth_demand_g_n_per_h = 0.075, .nitrogen_growth_actual_g_n_per_h = 0.06, .phosphorus_growth_actual_g_p_per_h = 0.006, .nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h = 0.1275, .nitrogen_assimilation_respiration_actual_g_c_per_h = 0.102 };
}

fn testLoss() metabolism.SecondaryRootSenescence {
    return .{ .respiration_oxygen_unlimited_g_c_per_h = 0.06, .respiration_actual_g_c_per_h = 0.05, .phenological_senescence_g_c_per_h = 0, .senesced_fraction = 0.1, .recyclable_carbon_g_c = 1, .recyclable_nitrogen_g_n = 0.1, .recyclable_phosphorus_g_p = 0.01 };
}

fn makeState(c: []f64, n: []f64, p: []f64, net_carbon: []f64, oxygen: []f64, carbon: []f64) State {
    return .{ .mobile_carbon_g_c_by_domain_layer = c, .mobile_nitrogen_g_n_by_domain_layer = n, .mobile_phosphorus_g_p_by_domain_layer = p, .net_root_carbon_flux_g_c_by_domain_layer = net_carbon, .oxygen_unlimited_respiration_g_c_by_domain_layer = oxygen, .carbon_unlimited_respiration_g_c_by_domain_layer = carbon };
}

fn testInputs(processes: []const metabolism.SecondaryRootResult, losses: []const metabolism.SecondaryRootSenescence) Inputs {
    return .{
        .biological_domain_count = 1,
        .soil_layer_count = 3,
        .root_axis_count = 2,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 2,
        .tip_active_by_domain_layer_axis = &.{ false, false, true, true, false, false },
        .metabolism_by_domain_layer_axis = processes,
        .senescence_by_domain_layer_axis = losses,
        .nonwoody_carbon_fraction = 0.5,
        .nonwoody_nitrogen_fraction = 0.5,
        .root_depth_m_by_domain_axis = &.{ 0.5, 0.3 },
        .seeding_depth_m = 0,
        .first_rooted_layer_index = 0,
        .deepest_axis_layer_index_by_axis = &.{ 2, 1 },
        .layer_bottom_depth_m = &.{ 0.2, 0.4, 0.6 },
        .primary_root_length_m_by_domain_layer_axis = &.{ 0.2, 0.15, 0.2, 0.15, 0.1, 0 },
    };
}

test "GROSUB primary mobile recovery and signed LL respiration preserve source order" {
    var c = [_]f64{ 5, 5, 5 };
    var n = [_]f64{ 1, 1, 1 };
    var p = [_]f64{ 0.1, 0.1, 0.1 };
    var actual = [_]f64{ 0, 0, 0 };
    var oxygen = [_]f64{ 0, 0, 0 };
    var carbon = [_]f64{ 0, 0, 0 };
    var wc = [_]f64{0} ** 3;
    var wn = [_]f64{0} ** 3;
    var wp = [_]f64{0} ** 3;
    var wa = [_]f64{0} ** 3;
    var wo = [_]f64{0} ** 3;
    var wx = [_]f64{0} ** 3;
    var processes = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)} ** 6;
    var losses = [_]metabolism.SecondaryRootSenescence{std.mem.zeroes(metabolism.SecondaryRootSenescence)} ** 6;
    processes[2] = testProcess();
    processes[3] = testProcess();
    losses[2] = testLoss();
    losses[3] = testLoss();
    const initial_c = c[1];
    const initial_n = n[1];
    const initial_p = p[1];
    try apply(makeState(&c, &n, &p, &actual, &oxygen, &carbon), makeState(&wc, &wn, &wp, &wa, &wo, &wx), testInputs(&processes, &losses));
    var total_actual: f64 = 0;
    for (actual) |value| total_actual += value;
    try std.testing.expectApproxEqAbs(@as(f64, -0.904), total_actual, 1e-15);
    try std.testing.expectApproxEqAbs(2 * (0.6 + 0.452), initial_c - c[1] + 2 * 0.05, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), initial_n - n[1] + 2 * 0.005, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.012), initial_p - p[1] + 2 * 0.001, 1e-12);
    // Line 6771 returns all recyclable P, without the 0.5 nonwoody fraction.
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), p[1], 1e-15);
}

test "GROSUB second primary tip overdraw rolls back first tip and LL respiration" {
    var c = [_]f64{ 5, 1.5, 5 };
    var n = [_]f64{ 1, 1, 1 };
    var p = [_]f64{ 0.1, 0.1, 0.1 };
    var actual = [_]f64{ 0, 0, 0 };
    var oxygen = [_]f64{ 0, 0, 0 };
    var carbon = [_]f64{ 0, 0, 0 };
    var workspace: [6][3]f64 = std.mem.zeroes([6][3]f64);
    var processes = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)} ** 6;
    var losses = [_]metabolism.SecondaryRootSenescence{std.mem.zeroes(metabolism.SecondaryRootSenescence)} ** 6;
    processes[2] = testProcess();
    processes[3] = testProcess();
    losses[2] = testLoss();
    losses[3] = testLoss();
    const before_c = c;
    const before_actual = actual;
    try std.testing.expectError(error.InvalidPrimaryRootMobileRespirationState, apply(makeState(&c, &n, &p, &actual, &oxygen, &carbon), makeState(&workspace[0], &workspace[1], &workspace[2], &workspace[3], &workspace[4], &workspace[5]), testInputs(&processes, &losses)));
    try std.testing.expectEqualDeep(before_c, c);
    try std.testing.expectEqualDeep(before_actual, actual);
}
