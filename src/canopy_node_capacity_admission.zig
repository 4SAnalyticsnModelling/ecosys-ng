const std = @import("std");

pub const BranchInputs = struct {
    photosynthetically_active: bool,
    branch_living: bool,
    c3_feedback_fraction: f64,
    annual_termination_fraction: f64,
};

pub const NodeInputs = struct {
    leaf_area_m2: f64,
    leaf_carbon_g_c: f64,
    leaf_protein_g: f64,
    presence_threshold: f64,
};

pub const Potentials = struct {
    c4_carboxylation_umol_per_m2_s: f64,
    c3_carboxylation_umol_per_m2_s: f64,
};

pub const BranchResult = struct {
    c3_feedback_fraction: f64,
    annual_termination_fraction: f64,
};

pub const NodeResult = struct {
    gate_evaluated: bool,
    protein_surface_density_g_per_m2: f64,
    capacity_admitted: bool,
    potentials: Potentials,
};

pub const Inputs = struct {
    species_branch_offsets: []const usize,
    branch_node_offsets: []const usize,
    branches: []const BranchInputs,
    nodes: []const NodeInputs,
    initial_potentials: []const Potentials,
};

pub const Scratch = struct {
    branches: []BranchResult,
    nodes: []NodeResult,
};

pub const Destination = struct {
    branches: []BranchResult,
    nodes: []NodeResult,
};

/// STOMATE.F 207--224 with reset counterparts at 625--637.
pub fn applyRuntimeTopology(
    inputs: Inputs,
    scratch: Scratch,
    destination: Destination,
) !void {
    try validateDimensions(inputs, scratch, destination);
    for (inputs.branches, scratch.branches, 0..) |branch, *branch_result, branch_index| {
        try validateBranch(branch);
        if (!branch.photosynthetically_active) {
            branch_result.* = .{
                .c3_feedback_fraction = 0,
                .annual_termination_fraction = 1,
            };
        } else {
            branch_result.* = .{
                .c3_feedback_fraction = branch.c3_feedback_fraction,
                .annual_termination_fraction = branch.annual_termination_fraction,
            };
        }

        const node_begin = inputs.branch_node_offsets[branch_index];
        const node_end = inputs.branch_node_offsets[branch_index + 1];
        for (node_begin..node_end) |node_index| {
            const initial = inputs.initial_potentials[node_index];
            if (!branch.photosynthetically_active) {
                scratch.nodes[node_index] = .{
                    .gate_evaluated = false,
                    .protein_surface_density_g_per_m2 = 0,
                    .capacity_admitted = false,
                    .potentials = std.mem.zeroes(Potentials),
                };
            } else if (!branch.branch_living) {
                try validatePotentials(initial);
                scratch.nodes[node_index] = .{
                    .gate_evaluated = false,
                    .protein_surface_density_g_per_m2 = 0,
                    .capacity_admitted = false,
                    .potentials = initial,
                };
            } else {
                const node = inputs.nodes[node_index];
                try validateNode(node);
                const density =
                    if (node.leaf_area_m2 > node.presence_threshold and
                    node.leaf_carbon_g_c > node.presence_threshold)
                        node.leaf_protein_g / node.leaf_area_m2
                    else
                        0;
                if (!std.math.isFinite(density) or density < 0)
                    return error.InvalidCanopyNodeProteinSurfaceDensity;
                const admitted = density > 0;
                if (admitted) try validatePotentials(initial);
                scratch.nodes[node_index] = .{
                    .gate_evaluated = true,
                    .protein_surface_density_g_per_m2 = density,
                    .capacity_admitted = admitted,
                    .potentials = if (admitted)
                        initial
                    else
                        std.mem.zeroes(Potentials),
                };
            }
        }
    }
    @memcpy(destination.branches, scratch.branches);
    @memcpy(destination.nodes, scratch.nodes);
}

fn validateDimensions(
    inputs: Inputs,
    scratch: Scratch,
    destination: Destination,
) !void {
    if (inputs.species_branch_offsets.len == 0 or
        inputs.branch_node_offsets.len != inputs.branches.len + 1 or
        inputs.nodes.len != inputs.initial_potentials.len or
        scratch.branches.len != inputs.branches.len or
        destination.branches.len != inputs.branches.len or
        scratch.nodes.len != inputs.nodes.len or
        destination.nodes.len != inputs.nodes.len)
        return error.CanopyNodeCapacityTopologyDimensionMismatch;
    if (inputs.species_branch_offsets[0] != 0 or
        inputs.species_branch_offsets[inputs.species_branch_offsets.len - 1] !=
            inputs.branches.len or
        inputs.branch_node_offsets[0] != 0 or
        inputs.branch_node_offsets[inputs.branch_node_offsets.len - 1] !=
            inputs.nodes.len)
        return error.InvalidCanopyNodeCapacityTopologyOffsets;
    for (0..inputs.species_branch_offsets.len - 1) |index|
        if (inputs.species_branch_offsets[index] >
            inputs.species_branch_offsets[index + 1])
            return error.InvalidCanopyNodeCapacityTopologyOffsets;
    for (0..inputs.branch_node_offsets.len - 1) |index|
        if (inputs.branch_node_offsets[index] >
            inputs.branch_node_offsets[index + 1])
            return error.InvalidCanopyNodeCapacityTopologyOffsets;
}

fn validateBranch(branch: BranchInputs) !void {
    if (!std.math.isFinite(branch.c3_feedback_fraction) or
        !std.math.isFinite(branch.annual_termination_fraction) or
        branch.c3_feedback_fraction < 0 or branch.c3_feedback_fraction > 1 or
        branch.annual_termination_fraction < 0 or
        branch.annual_termination_fraction > 1)
        return error.InvalidCanopyNodeCapacityBranchInput;
}

fn validateNode(node: NodeInputs) !void {
    inline for (@typeInfo(NodeInputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(node, field.name)))
            return error.NonFiniteCanopyNodeCapacityInput;
    if (node.leaf_area_m2 < 0 or node.leaf_carbon_g_c < 0 or
        node.leaf_protein_g < 0 or node.presence_threshold < 0)
        return error.InvalidCanopyNodeCapacityInput;
}

fn validatePotentials(potentials: Potentials) !void {
    if (!std.math.isFinite(potentials.c4_carboxylation_umol_per_m2_s) or
        !std.math.isFinite(potentials.c3_carboxylation_umol_per_m2_s) or
        potentials.c4_carboxylation_umol_per_m2_s < 0 or
        potentials.c3_carboxylation_umol_per_m2_s < 0)
        return error.InvalidCanopyNodeCapacityPotential;
}

test "living node admission uses strict area and leaf-carbon gates" {
    const species_offsets = [_]usize{ 0, 1 };
    const node_offsets = [_]usize{ 0, 3 };
    const branches = [_]BranchInputs{.{
        .photosynthetically_active = true,
        .branch_living = true,
        .c3_feedback_fraction = 0.5,
        .annual_termination_fraction = 1,
    }};
    const nodes = [_]NodeInputs{
        .{ .leaf_area_m2 = 2, .leaf_carbon_g_c = 4, .leaf_protein_g = 1, .presence_threshold = 1 },
        .{ .leaf_area_m2 = 1, .leaf_carbon_g_c = 4, .leaf_protein_g = 1, .presence_threshold = 1 },
        .{ .leaf_area_m2 = 2, .leaf_carbon_g_c = 1, .leaf_protein_g = 1, .presence_threshold = 1 },
    };
    const initial = [_]Potentials{
        .{ .c4_carboxylation_umol_per_m2_s = 3, .c3_carboxylation_umol_per_m2_s = 4 },
        .{ .c4_carboxylation_umol_per_m2_s = 5, .c3_carboxylation_umol_per_m2_s = 6 },
        .{ .c4_carboxylation_umol_per_m2_s = 7, .c3_carboxylation_umol_per_m2_s = 8 },
    };
    var scratch_branches: [1]BranchResult = undefined;
    var scratch_nodes: [3]NodeResult = undefined;
    var destination_branches: [1]BranchResult = undefined;
    var destination_nodes: [3]NodeResult = undefined;
    try applyRuntimeTopology(
        .{ .species_branch_offsets = &species_offsets, .branch_node_offsets = &node_offsets, .branches = &branches, .nodes = &nodes, .initial_potentials = &initial },
        .{ .branches = &scratch_branches, .nodes = &scratch_nodes },
        .{ .branches = &destination_branches, .nodes = &destination_nodes },
    );
    try std.testing.expect(destination_nodes[0].capacity_admitted);
    try std.testing.expectEqual(@as(f64, 0.5), destination_nodes[0].protein_surface_density_g_per_m2);
    try std.testing.expectEqualDeep(std.mem.zeroes(Potentials), destination_nodes[1].potentials);
    try std.testing.expectEqualDeep(std.mem.zeroes(Potentials), destination_nodes[2].potentials);
}

test "inactive phenology zeros feedback and every node potential" {
    const species_offsets = [_]usize{ 0, 1 };
    const node_offsets = [_]usize{ 0, 1 };
    const branches = [_]BranchInputs{.{ .photosynthetically_active = false, .branch_living = true, .c3_feedback_fraction = 0.8, .annual_termination_fraction = 0.2 }};
    const nodes = [_]NodeInputs{.{ .leaf_area_m2 = std.math.nan(f64), .leaf_carbon_g_c = 1, .leaf_protein_g = 1, .presence_threshold = 0 }};
    const initial = [_]Potentials{.{ .c4_carboxylation_umol_per_m2_s = std.math.nan(f64), .c3_carboxylation_umol_per_m2_s = 4 }};
    var scratch_branches: [1]BranchResult = undefined;
    var scratch_nodes: [1]NodeResult = undefined;
    var destination_branches: [1]BranchResult = undefined;
    var destination_nodes: [1]NodeResult = undefined;
    try applyRuntimeTopology(
        .{ .species_branch_offsets = &species_offsets, .branch_node_offsets = &node_offsets, .branches = &branches, .nodes = &nodes, .initial_potentials = &initial },
        .{ .branches = &scratch_branches, .nodes = &scratch_nodes },
        .{ .branches = &destination_branches, .nodes = &destination_nodes },
    );
    try std.testing.expectEqual(@as(f64, 0), destination_branches[0].c3_feedback_fraction);
    try std.testing.expectEqual(@as(f64, 1), destination_branches[0].annual_termination_fraction);
    try std.testing.expectEqualDeep(std.mem.zeroes(Potentials), destination_nodes[0].potentials);
}

test "dead branch preserves prior potentials without evaluating node gate" {
    const species_offsets = [_]usize{ 0, 1 };
    const node_offsets = [_]usize{ 0, 1 };
    const branches = [_]BranchInputs{.{ .photosynthetically_active = true, .branch_living = false, .c3_feedback_fraction = 0.5, .annual_termination_fraction = 1 }};
    const nodes = [_]NodeInputs{.{ .leaf_area_m2 = std.math.nan(f64), .leaf_carbon_g_c = 1, .leaf_protein_g = 1, .presence_threshold = 0 }};
    const initial = [_]Potentials{.{ .c4_carboxylation_umol_per_m2_s = 3, .c3_carboxylation_umol_per_m2_s = 4 }};
    var scratch_branches: [1]BranchResult = undefined;
    var scratch_nodes: [1]NodeResult = undefined;
    var destination_branches: [1]BranchResult = undefined;
    var destination_nodes: [1]NodeResult = undefined;
    try applyRuntimeTopology(
        .{ .species_branch_offsets = &species_offsets, .branch_node_offsets = &node_offsets, .branches = &branches, .nodes = &nodes, .initial_potentials = &initial },
        .{ .branches = &scratch_branches, .nodes = &scratch_nodes },
        .{ .branches = &destination_branches, .nodes = &destination_nodes },
    );
    try std.testing.expectEqualDeep(initial[0], destination_nodes[0].potentials);
    try std.testing.expect(!destination_nodes[0].gate_evaluated);
}

test "later invalid living node leaves both destination sets unchanged" {
    const species_offsets = [_]usize{ 0, 1, 2 };
    const node_offsets = [_]usize{ 0, 1, 2 };
    const branches = [_]BranchInputs{
        .{ .photosynthetically_active = true, .branch_living = true, .c3_feedback_fraction = 0.5, .annual_termination_fraction = 1 },
        .{ .photosynthetically_active = true, .branch_living = true, .c3_feedback_fraction = 0.5, .annual_termination_fraction = 1 },
    };
    const nodes = [_]NodeInputs{
        .{ .leaf_area_m2 = 2, .leaf_carbon_g_c = 2, .leaf_protein_g = 1, .presence_threshold = 0 },
        .{ .leaf_area_m2 = std.math.nan(f64), .leaf_carbon_g_c = 2, .leaf_protein_g = 1, .presence_threshold = 0 },
    };
    const initial = [_]Potentials{
        .{ .c4_carboxylation_umol_per_m2_s = 3, .c3_carboxylation_umol_per_m2_s = 4 },
        .{ .c4_carboxylation_umol_per_m2_s = 5, .c3_carboxylation_umol_per_m2_s = 6 },
    };
    var scratch_branches: [2]BranchResult = undefined;
    var scratch_nodes: [2]NodeResult = undefined;
    var destination_branches = [_]BranchResult{
        .{ .c3_feedback_fraction = 41, .annual_termination_fraction = 1 },
        .{ .c3_feedback_fraction = 42, .annual_termination_fraction = 1 },
    };
    var destination_nodes = [_]NodeResult{
        .{ .gate_evaluated = false, .protein_surface_density_g_per_m2 = 0, .capacity_admitted = false, .potentials = .{ .c4_carboxylation_umol_per_m2_s = 41, .c3_carboxylation_umol_per_m2_s = 41 } },
        .{ .gate_evaluated = false, .protein_surface_density_g_per_m2 = 0, .capacity_admitted = false, .potentials = .{ .c4_carboxylation_umol_per_m2_s = 42, .c3_carboxylation_umol_per_m2_s = 42 } },
    };
    try std.testing.expectError(
        error.NonFiniteCanopyNodeCapacityInput,
        applyRuntimeTopology(
            .{ .species_branch_offsets = &species_offsets, .branch_node_offsets = &node_offsets, .branches = &branches, .nodes = &nodes, .initial_potentials = &initial },
            .{ .branches = &scratch_branches, .nodes = &scratch_nodes },
            .{ .branches = &destination_branches, .nodes = &destination_nodes },
        ),
    );
    try std.testing.expectEqual(@as(f64, 41), destination_branches[0].c3_feedback_fraction);
    try std.testing.expectEqual(@as(f64, 42), destination_nodes[1].potentials.c3_carboxylation_umol_per_m2_s);
}
