const std = @import("std");

/// Harvested and litter-routed leaf ledgers accumulated over canopy layers,
/// branches, and leaf nodes. Non-woody components carry the `1` source index of
/// `FWODB`/`FWODLN`/`FWODLP`; woody components carry the `0` index.
pub const Ledger = struct {
    harvested_nonwoody_carbon_g_c: f64 = 0,
    harvested_nonwoody_nitrogen_g_n: f64 = 0,
    harvested_nonwoody_phosphorus_g_p: f64 = 0,
    litter_nonwoody_carbon_g_c: f64 = 0,
    litter_nonwoody_nitrogen_g_n: f64 = 0,
    litter_nonwoody_phosphorus_g_p: f64 = 0,
    harvested_woody_carbon_g_c: f64 = 0,
    harvested_woody_nitrogen_g_n: f64 = 0,
    harvested_woody_phosphorus_g_p: f64 = 0,
    litter_woody_carbon_g_c: f64 = 0,
    litter_woody_nitrogen_g_n: f64 = 0,
    litter_woody_phosphorus_g_p: f64 = 0,
};

/// Mutable per-node leaf state plus the layer aggregates reset by the source
/// layer loop. Slice layouts are runtime sized; there is no `JC`/`K=0..25`
/// ceiling.
pub const State = struct {
    /// `WGLFL`, indexed `(layer * node_count + node) * branch_count + branch`.
    leaf_carbon_g_c_by_layer_node_branch: []f64,
    /// `WGLFLN`, same layout.
    leaf_nitrogen_g_n_by_layer_node_branch: []f64,
    /// `WGLFLP`, same layout.
    leaf_phosphorus_g_p_by_layer_node_branch: []f64,
    /// `ARLFL`, same layout, m2.
    leaf_area_m2_by_layer_node_branch: []f64,
    /// `ARSTK`, indexed `layer * branch_count + branch`, m2.
    stalk_area_m2_by_layer_branch: []f64,
    /// `ARLFV`, indexed by canopy layer, m2.
    plant_leaf_area_m2_by_layer: []f64,
    /// `WGLFV`, indexed by canopy layer, g C.
    plant_leaf_carbon_g_c_by_layer: []f64,
    /// `ARSTV`, indexed by canopy layer, m2.
    plant_stalk_area_m2_by_layer: []f64,
};

pub const WoodyPartition = struct {
    /// `FWODB(1)`, `FWODLN(1)`, `FWODLP(1)`: non-woody C/N/P fractions.
    nonwoody_carbon_fraction: f64,
    nonwoody_nitrogen_fraction: f64,
    nonwoody_phosphorus_fraction: f64,
    /// `FWODB(0)`, `FWODLN(0)`, `FWODLP(0)`: woody C/N/P fractions.
    woody_carbon_fraction: f64,
    woody_nitrogen_fraction: f64,
    woody_phosphorus_fraction: f64,
};

pub const Inputs = struct {
    branch_count: usize,
    canopy_layer_count: usize,
    leaf_node_pool_count: usize,
    /// The source stalk-area node selector `K.EQ.1`.
    stalk_area_node_index: usize = 1,
    /// True for `IHVST` 4 or 6 (animal or insect grazing).
    grazing: bool,
    /// `WTLF`: plant leaf C mass, g C.
    plant_leaf_carbon_g_c: f64,
    /// `ZEROP2`: population-scaled plant presence threshold, g C.
    plant_presence_threshold_g_c: f64,
    /// `WHVSLF`: plant structural leaf C removal demand, g C h-1.
    plant_structural_leaf_removal_g_c: f64,
    /// `WGLFBL`: branch leaf C in canopy layer, g C, layout
    /// `layer * branch_count + branch`.
    branch_leaf_carbon_g_c_by_layer_branch: []const f64,
    /// `FHVST` from the layer selector: fraction of layer mass not harvested.
    plant_retained_fraction_by_layer: []const f64,
    /// `FHVSH` from the layer selector: fraction not exported from the PFT.
    harvest_retained_fraction_by_layer: []const f64,
    woody_partition: WoodyPartition,
};

fn nodeIndex(inputs: Inputs, layer: usize, node: usize, branch: usize) usize {
    return (layer * inputs.leaf_node_pool_count + node) * inputs.branch_count + branch;
}

fn layerBranchIndex(inputs: Inputs, layer: usize, branch: usize) usize {
    return layer * inputs.branch_count + branch;
}

fn validateFraction(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0.0 or value > 1.0) return error.InvalidLeafNodeHarvestFraction;
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.branch_count == 0 or inputs.canopy_layer_count == 0 or inputs.leaf_node_pool_count == 0)
        return error.LeafNodeHarvestDimensionMismatch;
    if (inputs.stalk_area_node_index >= inputs.leaf_node_pool_count)
        return error.LeafNodeHarvestDimensionMismatch;
    const layer_branch = std.math.mul(usize, inputs.canopy_layer_count, inputs.branch_count) catch
        return error.LeafNodeHarvestDimensionOverflow;
    if (inputs.branch_leaf_carbon_g_c_by_layer_branch.len != layer_branch)
        return error.LeafNodeHarvestDimensionMismatch;
    if (inputs.plant_retained_fraction_by_layer.len != inputs.canopy_layer_count or
        inputs.harvest_retained_fraction_by_layer.len != inputs.canopy_layer_count)
        return error.LeafNodeHarvestDimensionMismatch;
    inline for (.{
        inputs.plant_leaf_carbon_g_c,
        inputs.plant_presence_threshold_g_c,
        inputs.plant_structural_leaf_removal_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidLeafNodeHarvestInput;
    for (inputs.branch_leaf_carbon_g_c_by_layer_branch) |value|
        if (!std.math.isFinite(value)) return error.InvalidLeafNodeHarvestInput;
    for (inputs.plant_retained_fraction_by_layer) |value| try validateFraction(value);
    for (inputs.harvest_retained_fraction_by_layer) |value| try validateFraction(value);
    inline for (@typeInfo(WoodyPartition).@"struct".fields) |field|
        try validateFraction(@field(inputs.woody_partition, field.name));
}

fn validateState(state: State, inputs: Inputs) !void {
    const layer_nodes = std.math.mul(usize, inputs.canopy_layer_count, inputs.leaf_node_pool_count) catch
        return error.LeafNodeHarvestDimensionOverflow;
    const node_count = std.math.mul(usize, layer_nodes, inputs.branch_count) catch
        return error.LeafNodeHarvestDimensionOverflow;
    const layer_branch = std.math.mul(usize, inputs.canopy_layer_count, inputs.branch_count) catch
        return error.LeafNodeHarvestDimensionOverflow;
    inline for (.{
        state.leaf_carbon_g_c_by_layer_node_branch,
        state.leaf_nitrogen_g_n_by_layer_node_branch,
        state.leaf_phosphorus_g_p_by_layer_node_branch,
        state.leaf_area_m2_by_layer_node_branch,
    }) |values| {
        if (values.len != node_count) return error.LeafNodeHarvestDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidLeafNodeHarvestState;
    }
    if (state.stalk_area_m2_by_layer_branch.len != layer_branch) return error.LeafNodeHarvestDimensionMismatch;
    for (state.stalk_area_m2_by_layer_branch) |value|
        if (!std.math.isFinite(value) or value < 0.0) return error.InvalidLeafNodeHarvestState;
    inline for (.{
        state.plant_leaf_area_m2_by_layer,
        state.plant_leaf_carbon_g_c_by_layer,
        state.plant_stalk_area_m2_by_layer,
    }) |values| {
        if (values.len != inputs.canopy_layer_count) return error.LeafNodeHarvestDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidLeafNodeHarvestState;
    }
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        @memcpy(@field(destination, field.name), @field(source, field.name));
}

pub const Result = struct {
    ledger: Ledger,
    /// The source `FHVST` value left in scope after the node loops. It scales
    /// `ARSTV` for every layer and carries between layers exactly as Fortran
    /// storage does.
    carried_plant_retained_fraction: f64,
};

/// Exact GROSUB 8850--8943 leaf-node harvest ledgers, remaining leaf state, and
/// layer aggregate reset. Traversal is canopy layer top-to-bottom, branch
/// ascending, leaf node descending, matching source loops 9865/9855/9845.
///
/// The trailing `ARSTV` scaling in source line 8943 uses whichever `FHVST` was
/// last assigned, including the grazing per-node assignment, so the carried
/// fraction is an explicit input and output rather than a per-layer value.
pub fn apply(
    state: State,
    workspace: State,
    inputs: Inputs,
    incoming_ledger: Ledger,
    incoming_carried_plant_retained_fraction: f64,
) !Result {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    try validateFraction(incoming_carried_plant_retained_fraction);
    inline for (@typeInfo(Ledger).@"struct".fields) |field|
        if (!std.math.isFinite(@field(incoming_ledger, field.name)) or
            @field(incoming_ledger, field.name) < 0.0) return error.InvalidLeafNodeHarvestLedger;

    copyState(workspace, state);
    var ledger = incoming_ledger;
    var carried_retained = incoming_carried_plant_retained_fraction;
    const partition = inputs.woody_partition;

    var reverse_layer = inputs.canopy_layer_count;
    while (reverse_layer > 0) {
        reverse_layer -= 1;
        const layer = reverse_layer;
        for (0..inputs.branch_count) |branch| {
            var layer_retained = inputs.plant_retained_fraction_by_layer[layer];
            var export_retained = inputs.harvest_retained_fraction_by_layer[layer];
            var remaining_demand_g_c: f64 = 0;
            if (inputs.grazing and inputs.plant_leaf_carbon_g_c > inputs.plant_presence_threshold_g_c) {
                remaining_demand_g_c = inputs.plant_structural_leaf_removal_g_c *
                    @max(0.0, inputs.branch_leaf_carbon_g_c_by_layer_branch[layerBranchIndex(inputs, layer, branch)]) /
                    inputs.plant_leaf_carbon_g_c;
                if (!std.math.isFinite(remaining_demand_g_c)) return error.NonFiniteLeafNodeHarvestDemand;
            }
            var reverse_node = inputs.leaf_node_pool_count;
            while (reverse_node > 0) {
                reverse_node -= 1;
                const node = reverse_node;
                if (!(!inputs.grazing or remaining_demand_g_c > 0.0)) continue;
                const index = nodeIndex(inputs, layer, node, branch);
                const node_carbon_g_c = workspace.leaf_carbon_g_c_by_layer_node_branch[index];
                if (inputs.grazing) {
                    if (node_carbon_g_c > remaining_demand_g_c) {
                        layer_retained = @max(0.0, @min(1.0, (node_carbon_g_c - remaining_demand_g_c) / node_carbon_g_c));
                        export_retained = layer_retained;
                    } else {
                        layer_retained = 1.0;
                        export_retained = 1.0;
                    }
                }
                carried_retained = layer_retained;
                const node_nitrogen_g_n = workspace.leaf_nitrogen_g_n_by_layer_node_branch[index];
                const node_phosphorus_g_p = workspace.leaf_phosphorus_g_p_by_layer_node_branch[index];
                const removed_fraction = 1.0 - export_retained;
                const litter_fraction = export_retained - layer_retained;

                remaining_demand_g_c -= (1.0 - layer_retained) * node_carbon_g_c;
                ledger.harvested_nonwoody_carbon_g_c += removed_fraction * node_carbon_g_c * partition.nonwoody_carbon_fraction;
                ledger.harvested_nonwoody_nitrogen_g_n += removed_fraction * node_nitrogen_g_n * partition.nonwoody_nitrogen_fraction;
                ledger.harvested_nonwoody_phosphorus_g_p += removed_fraction * node_phosphorus_g_p * partition.nonwoody_phosphorus_fraction;
                ledger.litter_nonwoody_carbon_g_c += litter_fraction * node_carbon_g_c * partition.nonwoody_carbon_fraction;
                ledger.litter_nonwoody_nitrogen_g_n += litter_fraction * node_nitrogen_g_n * partition.nonwoody_nitrogen_fraction;
                ledger.litter_nonwoody_phosphorus_g_p += litter_fraction * node_phosphorus_g_p * partition.nonwoody_phosphorus_fraction;
                ledger.harvested_woody_carbon_g_c += removed_fraction * node_carbon_g_c * partition.woody_carbon_fraction;
                ledger.harvested_woody_nitrogen_g_n += removed_fraction * node_nitrogen_g_n * partition.woody_nitrogen_fraction;
                ledger.harvested_woody_phosphorus_g_p += removed_fraction * node_phosphorus_g_p * partition.woody_phosphorus_fraction;
                ledger.litter_woody_carbon_g_c += litter_fraction * node_carbon_g_c * partition.woody_carbon_fraction;
                ledger.litter_woody_nitrogen_g_n += litter_fraction * node_nitrogen_g_n * partition.woody_nitrogen_fraction;
                ledger.litter_woody_phosphorus_g_p += litter_fraction * node_phosphorus_g_p * partition.woody_phosphorus_fraction;

                workspace.leaf_carbon_g_c_by_layer_node_branch[index] = layer_retained * node_carbon_g_c;
                workspace.leaf_nitrogen_g_n_by_layer_node_branch[index] = layer_retained * node_nitrogen_g_n;
                workspace.leaf_phosphorus_g_p_by_layer_node_branch[index] = layer_retained * node_phosphorus_g_p;
                workspace.leaf_area_m2_by_layer_node_branch[index] =
                    layer_retained * workspace.leaf_area_m2_by_layer_node_branch[index];
                if (node == inputs.stalk_area_node_index) {
                    const stalk = layerBranchIndex(inputs, layer, branch);
                    workspace.stalk_area_m2_by_layer_branch[stalk] =
                        layer_retained * workspace.stalk_area_m2_by_layer_branch[stalk];
                }
            }
        }
        workspace.plant_leaf_area_m2_by_layer[layer] = 0.0;
        workspace.plant_leaf_carbon_g_c_by_layer[layer] = 0.0;
        workspace.plant_stalk_area_m2_by_layer[layer] =
            workspace.plant_stalk_area_m2_by_layer[layer] * carried_retained;
    }

    inline for (@typeInfo(Ledger).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ledger, field.name)) or
            @field(ledger, field.name) < 0.0) return error.InvalidLeafNodeHarvestLedger;
    try validateState(workspace, inputs);
    try validateFraction(carried_retained);
    copyState(state, workspace);
    return .{ .ledger = ledger, .carried_plant_retained_fraction = carried_retained };
}

const test_partition: WoodyPartition = .{
    .nonwoody_carbon_fraction = 0.75,
    .nonwoody_nitrogen_fraction = 0.8,
    .nonwoody_phosphorus_fraction = 0.9,
    .woody_carbon_fraction = 0.25,
    .woody_nitrogen_fraction = 0.2,
    .woody_phosphorus_fraction = 0.1,
};

const TestState = struct {
    carbon: [4]f64 = .{ 10, 20, 30, 40 },
    nitrogen: [4]f64 = .{ 1, 2, 3, 4 },
    phosphorus: [4]f64 = .{ 0.1, 0.2, 0.3, 0.4 },
    area: [4]f64 = .{ 2, 4, 6, 8 },
    stalk: [2]f64 = .{ 5, 7 },
    plant_area: [2]f64 = .{ 3, 9 },
    plant_carbon: [2]f64 = .{ 11, 13 },
    plant_stalk: [2]f64 = .{ 20, 40 },

    fn view(self: *TestState) State {
        return .{
            .leaf_carbon_g_c_by_layer_node_branch = &self.carbon,
            .leaf_nitrogen_g_n_by_layer_node_branch = &self.nitrogen,
            .leaf_phosphorus_g_p_by_layer_node_branch = &self.phosphorus,
            .leaf_area_m2_by_layer_node_branch = &self.area,
            .stalk_area_m2_by_layer_branch = &self.stalk,
            .plant_leaf_area_m2_by_layer = &self.plant_area,
            .plant_leaf_carbon_g_c_by_layer = &self.plant_carbon,
            .plant_stalk_area_m2_by_layer = &self.plant_stalk,
        };
    }
};

fn testInputs(grazing: bool, retained: *const [2]f64, exported: *const [2]f64) Inputs {
    return .{
        .branch_count = 1,
        .canopy_layer_count = 2,
        .leaf_node_pool_count = 2,
        .grazing = grazing,
        .plant_leaf_carbon_g_c = 100,
        .plant_presence_threshold_g_c = 1e-6,
        .plant_structural_leaf_removal_g_c = 0,
        .branch_leaf_carbon_g_c_by_layer_branch = &.{ 30, 70 },
        .plant_retained_fraction_by_layer = retained,
        .harvest_retained_fraction_by_layer = exported,
        .woody_partition = test_partition,
    };
}

test "non-grazing leaf node harvest splits export and litter by layer fractions" {
    var state: TestState = .{};
    var workspace: TestState = .{};
    const inputs = testInputs(false, &.{ 0.5, 0.25 }, &.{ 0.75, 0.5 });
    const result = try apply(state.view(), workspace.view(), inputs, .{}, 1);

    // Layer 0 nodes retain 0.5 with export retention 0.75; layer 1 retains 0.25
    // with export retention 0.5.
    try std.testing.expectEqualSlices(f64, &.{ 5, 10, 7.5, 10 }, &state.carbon);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 1.5, 2 }, &state.area);
    // Only node index 1 scales stalk area.
    try std.testing.expectEqualSlices(f64, &.{ 2.5, 1.75 }, &state.stalk);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &state.plant_area);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &state.plant_carbon);
    // Layer 1 is scaled by its own 0.25 and layer 0 by 0.5, so 40 and 20 both
    // land on 10.
    try std.testing.expectEqualSlices(f64, &.{ 10, 10 }, &state.plant_stalk);
    try std.testing.expectEqual(@as(f64, 0.5), result.carried_plant_retained_fraction);

    const removed_carbon = 0.25 * (10.0 + 20.0) + 0.5 * (30.0 + 40.0);
    const litter_carbon = 0.25 * (10.0 + 20.0) + 0.25 * (30.0 + 40.0);
    try std.testing.expectApproxEqRel(removed_carbon * 0.75, result.ledger.harvested_nonwoody_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqRel(removed_carbon * 0.25, result.ledger.harvested_woody_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqRel(litter_carbon * 0.75, result.ledger.litter_nonwoody_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqRel(litter_carbon * 0.25, result.ledger.litter_woody_carbon_g_c, 1e-12);
}

test "non-grazing leaf carbon closes against remaining plus removed plus litter" {
    var state: TestState = .{};
    var workspace: TestState = .{};
    const initial = state.carbon;
    const inputs = testInputs(false, &.{ 0.5, 0.25 }, &.{ 0.75, 0.5 });
    const result = try apply(state.view(), workspace.view(), inputs, .{}, 1);
    var initial_total: f64 = 0;
    for (initial) |value| initial_total += value;
    var remaining_total: f64 = 0;
    for (state.carbon) |value| remaining_total += value;
    const accounted = remaining_total +
        result.ledger.harvested_nonwoody_carbon_g_c + result.ledger.harvested_woody_carbon_g_c +
        result.ledger.litter_nonwoody_carbon_g_c + result.ledger.litter_woody_carbon_g_c;
    try std.testing.expectApproxEqRel(initial_total, accounted, 1e-12);
}

test "grazing consumes branch layer demand from the top node downward" {
    var state: TestState = .{};
    var workspace: TestState = .{};
    var inputs = testInputs(true, &.{ 1, 1 }, &.{ 1, 1 });
    inputs.plant_structural_leaf_removal_g_c = 50;
    const result = try apply(state.view(), workspace.view(), inputs, .{}, 1);

    // Layer 1 demand is 50 * 70 / 100 = 35; node 1 holds 40 g C so it absorbs
    // all of it and node 0 sees no residual demand. Layer 0 demand is 15 and
    // node 1 holds 20 g C.
    try std.testing.expectApproxEqRel(@as(f64, 5), state.carbon[1], 1e-12);
    try std.testing.expectApproxEqRel(@as(f64, 5), state.carbon[3], 1e-12);
    try std.testing.expectEqual(@as(f64, 10), state.carbon[0]);
    try std.testing.expectEqual(@as(f64, 30), state.carbon[2]);
    // Grazing routes everything to export, never to litter.
    try std.testing.expectEqual(@as(f64, 0), result.ledger.litter_nonwoody_carbon_g_c);
    try std.testing.expectApproxEqRel(
        50.0 * 0.75,
        result.ledger.harvested_nonwoody_carbon_g_c,
        1e-12,
    );
}

test "grazing below the plant presence threshold leaves every node unchanged" {
    var state: TestState = .{};
    var workspace: TestState = .{};
    const initial = state.carbon;
    var inputs = testInputs(true, &.{ 0, 0 }, &.{ 0, 0 });
    inputs.plant_structural_leaf_removal_g_c = 50;
    inputs.plant_leaf_carbon_g_c = 1e-9;
    inputs.plant_presence_threshold_g_c = 1e-6;
    const result = try apply(state.view(), workspace.view(), inputs, .{}, 1);
    try std.testing.expectEqualSlices(f64, &initial, &state.carbon);
    try std.testing.expectEqual(@as(f64, 0), result.ledger.harvested_nonwoody_carbon_g_c);
    // No node assigned a fraction, so the incoming carried value survives and
    // scales the plant stalk area unchanged.
    try std.testing.expectEqual(@as(f64, 1), result.carried_plant_retained_fraction);
    try std.testing.expectEqualSlices(f64, &.{ 20, 40 }, &state.plant_stalk);
}

test "ledger accumulates onto incoming totals" {
    var state: TestState = .{};
    var workspace: TestState = .{};
    const inputs = testInputs(false, &.{ 0.5, 0.5 }, &.{ 0.5, 0.5 });
    const seeded: Ledger = .{ .harvested_nonwoody_carbon_g_c = 100 };
    const result = try apply(state.view(), workspace.view(), inputs, seeded, 1);
    try std.testing.expectApproxEqRel(
        100.0 + 0.5 * 100.0 * 0.75,
        result.ledger.harvested_nonwoody_carbon_g_c,
        1e-12,
    );
}

test "invalid dimensions and non-finite input fail before mutation" {
    var state: TestState = .{};
    var workspace: TestState = .{};
    var inputs = testInputs(false, &.{ 0.5, 0.5 }, &.{ 0.5, 0.5 });
    inputs.canopy_layer_count = 3;
    try std.testing.expectError(error.LeafNodeHarvestDimensionMismatch, apply(state.view(), workspace.view(), inputs, .{}, 1));
    inputs = testInputs(false, &.{ 0.5, 0.5 }, &.{ 0.5, 0.5 });
    inputs.plant_structural_leaf_removal_g_c = std.math.nan(f64);
    try std.testing.expectError(error.InvalidLeafNodeHarvestInput, apply(state.view(), workspace.view(), inputs, .{}, 1));
    inputs = testInputs(false, &.{ 0.5, 1.5 }, &.{ 0.5, 0.5 });
    try std.testing.expectError(error.InvalidLeafNodeHarvestFraction, apply(state.view(), workspace.view(), inputs, .{}, 1));
    try std.testing.expectEqualSlices(f64, &.{ 10, 20, 30, 40 }, &state.carbon);
}

test "decomposing the sweep by layer batches matches a single call" {
    var whole: TestState = .{};
    var whole_workspace: TestState = .{};
    var split: TestState = .{};
    var split_workspace: TestState = .{};
    const retained = [2]f64{ 0.4, 0.6 };
    const exported = [2]f64{ 0.7, 0.8 };
    const inputs = testInputs(false, &retained, &exported);
    const whole_result = try apply(whole.view(), whole_workspace.view(), inputs, .{}, 1);

    // A per-layer decomposition preserves top-down order because each layer is
    // scaled only by its own fractions.
    var upper = inputs;
    upper.plant_retained_fraction_by_layer = &.{ 1, retained[1] };
    upper.harvest_retained_fraction_by_layer = &.{ 1, exported[1] };
    const first = try apply(split.view(), split_workspace.view(), upper, .{}, 1);
    var lower = inputs;
    lower.plant_retained_fraction_by_layer = &.{ retained[0], 1 };
    lower.harvest_retained_fraction_by_layer = &.{ exported[0], 1 };
    const second = try apply(split.view(), split_workspace.view(), lower, first.ledger, first.carried_plant_retained_fraction);
    try std.testing.expectEqualSlices(f64, &whole.carbon, &split.carbon);
    try std.testing.expectApproxEqRel(
        whole_result.ledger.harvested_nonwoody_carbon_g_c,
        second.ledger.harvested_nonwoody_carbon_g_c,
        1e-12,
    );
}
