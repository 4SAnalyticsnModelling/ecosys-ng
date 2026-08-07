const std = @import("std");

pub const State = struct {
    active_root_carbon_g_c_by_species_compartment_layer: []f64,
    actual_root_carbon_g_c_by_species_compartment_layer: []f64,
    total_plant_respiration_g_c_by_species: []f64,
    ecosystem_respiration_g_c: []f64,
    autotrophic_respiration_g_c: []f64,
};

pub const Inputs = struct {
    plant_species_count: usize,
    soil_layer_count: usize,
    maximum_association_compartment_count: usize,
    maximum_root_axis_count: usize,
    shared_planting_layer_index: usize,
    deepest_rooted_layer_index_by_species: []const usize,
    association_compartment_count_by_species: []const usize,
    root_axis_count_by_species: []const usize,
    initial_layer_index_by_species_root_axis: []const usize,
    secondary_root_carbon_g_c_by_species_compartment_layer_axis: []const f64,
    primary_root_carbon_g_c_by_species_compartment_layer_axis: []const f64,
    basal_root_carbon_g_c_by_species_compartment_axis: []const f64,
    root_respiration_g_c_by_species_compartment_layer: []const f64,
};

fn layerIndex(inputs: Inputs, species: usize, compartment: usize, layer: usize) usize {
    return (species * inputs.maximum_association_compartment_count + compartment) * inputs.soil_layer_count + layer;
}

fn axisLayerIndex(inputs: Inputs, species: usize, compartment: usize, layer: usize, axis: usize) usize {
    return layerIndex(inputs, species, compartment, layer) * inputs.maximum_root_axis_count + axis;
}

fn axisIndex(inputs: Inputs, species: usize, compartment: usize, axis: usize) usize {
    return (species * inputs.maximum_association_compartment_count + compartment) * inputs.maximum_root_axis_count + axis;
}

fn speciesAxisIndex(inputs: Inputs, species: usize, axis: usize) usize {
    return species * inputs.maximum_root_axis_count + axis;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn counts(inputs: Inputs) !struct { layers: usize, axes: usize, compartment_axes: usize, species_axes: usize } {
    const species_compartments = std.math.mul(usize, inputs.plant_species_count, inputs.maximum_association_compartment_count) catch return error.RootTotalsDimensionOverflow;
    const layers = std.math.mul(usize, species_compartments, inputs.soil_layer_count) catch return error.RootTotalsDimensionOverflow;
    return .{
        .layers = layers,
        .axes = std.math.mul(usize, layers, inputs.maximum_root_axis_count) catch return error.RootTotalsDimensionOverflow,
        .compartment_axes = std.math.mul(usize, species_compartments, inputs.maximum_root_axis_count) catch return error.RootTotalsDimensionOverflow,
        .species_axes = std.math.mul(usize, inputs.plant_species_count, inputs.maximum_root_axis_count) catch return error.RootTotalsDimensionOverflow,
    };
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.plant_species_count == 0 or inputs.soil_layer_count == 0 or inputs.maximum_association_compartment_count == 0 or inputs.maximum_root_axis_count == 0 or inputs.shared_planting_layer_index >= inputs.soil_layer_count) return error.RootTotalsDimensionMismatch;
    const expected = try counts(inputs);
    if (inputs.deepest_rooted_layer_index_by_species.len != inputs.plant_species_count or
        inputs.association_compartment_count_by_species.len != inputs.plant_species_count or
        inputs.root_axis_count_by_species.len != inputs.plant_species_count or
        inputs.initial_layer_index_by_species_root_axis.len != expected.species_axes or
        inputs.secondary_root_carbon_g_c_by_species_compartment_layer_axis.len != expected.axes or
        inputs.primary_root_carbon_g_c_by_species_compartment_layer_axis.len != expected.axes or
        inputs.basal_root_carbon_g_c_by_species_compartment_axis.len != expected.compartment_axes or
        inputs.root_respiration_g_c_by_species_compartment_layer.len != expected.layers) return error.RootTotalsDimensionMismatch;
    for (0..inputs.plant_species_count) |species| {
        if (inputs.deepest_rooted_layer_index_by_species[species] < inputs.shared_planting_layer_index or inputs.deepest_rooted_layer_index_by_species[species] >= inputs.soil_layer_count or
            inputs.association_compartment_count_by_species[species] == 0 or inputs.association_compartment_count_by_species[species] > inputs.maximum_association_compartment_count or
            inputs.root_axis_count_by_species[species] > inputs.maximum_root_axis_count) return error.RootTotalsDimensionMismatch;
        for (0..inputs.root_axis_count_by_species[species]) |axis| {
            const axis_layer = inputs.initial_layer_index_by_species_root_axis[speciesAxisIndex(inputs, species, axis)];
            if (axis_layer < inputs.shared_planting_layer_index or axis_layer > inputs.deepest_rooted_layer_index_by_species[species]) return error.RootTotalsDimensionMismatch;
        }
    }
    inline for (.{ inputs.secondary_root_carbon_g_c_by_species_compartment_layer_axis, inputs.primary_root_carbon_g_c_by_species_compartment_layer_axis, inputs.basal_root_carbon_g_c_by_species_compartment_axis }) |values|
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootTotalsInput;
    for (inputs.root_respiration_g_c_by_species_compartment_layer) |value| if (!std.math.isFinite(value)) return error.InvalidRootTotalsInput;
}

fn validateState(state: State, inputs: Inputs) !void {
    const expected = try counts(inputs);
    if (state.active_root_carbon_g_c_by_species_compartment_layer.len != expected.layers or state.actual_root_carbon_g_c_by_species_compartment_layer.len != expected.layers or state.total_plant_respiration_g_c_by_species.len != inputs.plant_species_count or state.ecosystem_respiration_g_c.len != 1 or state.autotrophic_respiration_g_c.len != 1) return error.RootTotalsDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| for (@field(state, field.name)) |value| {
        const signed = std.mem.indexOf(u8, field.name, "respiration") != null;
        if (!std.math.isFinite(value) or (!signed and value < 0.0)) return error.InvalidRootTotalsState;
    };
}

/// Exact GROSUB 8210--8228 root/nodule totals. Runtime topology is plant
/// species x association compartments x layers x root axes. Carbon and
/// respiration are g C per biological timestep.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    copyState(workspace, state);
    for (0..inputs.plant_species_count) |species| {
        for (0..inputs.association_compartment_count_by_species[species]) |compartment| {
            for (inputs.shared_planting_layer_index..inputs.deepest_rooted_layer_index_by_species[species] + 1) |layer| {
                const layer_index = layerIndex(inputs, species, compartment, layer);
                workspace.active_root_carbon_g_c_by_species_compartment_layer[layer_index] = 0.0;
                workspace.actual_root_carbon_g_c_by_species_compartment_layer[layer_index] = 0.0;
                for (0..inputs.root_axis_count_by_species[species]) |axis| {
                    const source_index = axisLayerIndex(inputs, species, compartment, layer, axis);
                    workspace.active_root_carbon_g_c_by_species_compartment_layer[layer_index] += inputs.secondary_root_carbon_g_c_by_species_compartment_layer_axis[source_index];
                    workspace.actual_root_carbon_g_c_by_species_compartment_layer[layer_index] += inputs.secondary_root_carbon_g_c_by_species_compartment_layer_axis[source_index] + inputs.primary_root_carbon_g_c_by_species_compartment_layer_axis[source_index];
                }
                const respiration_g_c = inputs.root_respiration_g_c_by_species_compartment_layer[layer_index];
                workspace.total_plant_respiration_g_c_by_species[species] += respiration_g_c;
                workspace.ecosystem_respiration_g_c[0] += respiration_g_c;
                workspace.autotrophic_respiration_g_c[0] += respiration_g_c;
            }
            for (0..inputs.root_axis_count_by_species[species]) |axis| {
                const initial_layer = inputs.initial_layer_index_by_species_root_axis[speciesAxisIndex(inputs, species, axis)];
                workspace.active_root_carbon_g_c_by_species_compartment_layer[layerIndex(inputs, species, compartment, initial_layer)] += inputs.basal_root_carbon_g_c_by_species_compartment_axis[axisIndex(inputs, species, compartment, axis)];
            }
        }
    }
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

test "GROSUB root totals retain active actual and respiration definitions" {
    const inputs: Inputs = .{
        .plant_species_count = 1,
        .soil_layer_count = 2,
        .maximum_association_compartment_count = 2,
        .maximum_root_axis_count = 2,
        .shared_planting_layer_index = 0,
        .deepest_rooted_layer_index_by_species = &.{1},
        .association_compartment_count_by_species = &.{2},
        .root_axis_count_by_species = &.{2},
        .initial_layer_index_by_species_root_axis = &.{ 0, 1 },
        .secondary_root_carbon_g_c_by_species_compartment_layer_axis = &.{ 1, 2, 3, 4, 5, 6, 7, 8 },
        .primary_root_carbon_g_c_by_species_compartment_layer_axis = &.{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 },
        .basal_root_carbon_g_c_by_species_compartment_axis = &.{ 10, 20, 30, 40 },
        .root_respiration_g_c_by_species_compartment_layer = &.{ -0.1, -0.2, -0.3, -0.4 },
    };
    var active = [_]f64{0} ** 4;
    var actual = [_]f64{0} ** 4;
    var plant_resp = [_]f64{1};
    var ecosystem = [_]f64{2};
    var autotrophic = [_]f64{3};
    var work_active = [_]f64{0} ** 4;
    var work_actual = [_]f64{0} ** 4;
    var work_plant_resp = [_]f64{0};
    var work_ecosystem = [_]f64{0};
    var work_autotrophic = [_]f64{0};
    const state: State = .{ .active_root_carbon_g_c_by_species_compartment_layer = &active, .actual_root_carbon_g_c_by_species_compartment_layer = &actual, .total_plant_respiration_g_c_by_species = &plant_resp, .ecosystem_respiration_g_c = &ecosystem, .autotrophic_respiration_g_c = &autotrophic };
    const workspace: State = .{ .active_root_carbon_g_c_by_species_compartment_layer = &work_active, .actual_root_carbon_g_c_by_species_compartment_layer = &work_actual, .total_plant_respiration_g_c_by_species = &work_plant_resp, .ecosystem_respiration_g_c = &work_ecosystem, .autotrophic_respiration_g_c = &work_autotrophic };
    try apply(state, workspace, inputs);
    try std.testing.expectEqual(@as(f64, 13), active[0]);
    try std.testing.expectApproxEqAbs(3.3, actual[0], 1e-14);
    try std.testing.expectEqual(@as(f64, 27), active[1]);
    try std.testing.expectApproxEqAbs(0.0, plant_resp[0], 1e-14);
    try std.testing.expectApproxEqAbs(1.0, ecosystem[0], 1e-14);
    try std.testing.expectApproxEqAbs(2.0, autotrophic[0], 1e-14);
}
