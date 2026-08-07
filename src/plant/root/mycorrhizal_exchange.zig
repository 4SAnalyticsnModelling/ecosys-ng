const std = @import("std");

pub const AssociationStatus = enum { absent, active };
pub const SaltMode = enum { static_equilibrium, dynamic };

pub const standard_salt_species_count: usize = 8;
pub const SaltSpecies = enum(usize) {
    aluminum,
    iron,
    calcium,
    magnesium,
    sodium,
    potassium,
    sulfate,
    chloride,
};

pub const State = struct {
    root_mobile_carbon_g_c_by_layer: []f64,
    mycorrhizal_mobile_carbon_g_c_by_layer: []f64,
    root_mobile_nitrogen_g_n_by_layer: []f64,
    mycorrhizal_mobile_nitrogen_g_n_by_layer: []f64,
    root_mobile_phosphorus_g_p_by_layer: []f64,
    mycorrhizal_mobile_phosphorus_g_p_by_layer: []f64,
    /// Layout is [salt species][soil layer][root=0,mycorrhiza=1], mol.
    salt_mol_by_species_layer_compartment: []f64,
};

pub const Inputs = struct {
    soil_layer_count: usize,
    salt_species_count: usize,
    shared_planting_layer_index: usize,
    shared_deepest_rooted_layer_index: usize,
    layer_active: []const bool,
    association_status: AssociationStatus,
    salt_mode: SaltMode,
    root_water_volume_m3_by_layer: []const f64,
    mycorrhizal_water_volume_m3_by_layer: []const f64,
    minimum_mycorrhiza_to_root_volume_ratio: f64,
    exchange_rate_per_h: f64,
    biological_timestep_h: f64,
    minimum_pool_g_c_by_layer: []const f64,
};

fn compartmentIndex(layer_count: usize, species: usize, layer: usize, compartment: usize) usize {
    return (species * layer_count + layer) * 2 + compartment;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, inputs: Inputs) !void {
    const layer_count = inputs.soil_layer_count;
    const salt_count = std.math.mul(usize, inputs.salt_species_count, layer_count) catch return error.RootMycorrhizalExchangeDimensionOverflow;
    const salt_compartment_count = std.math.mul(usize, salt_count, 2) catch return error.RootMycorrhizalExchangeDimensionOverflow;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const expected = if (std.mem.eql(u8, field.name, "salt_mol_by_species_layer_compartment")) salt_compartment_count else layer_count;
        if (values.len != expected) return error.RootMycorrhizalExchangeDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootMycorrhizalExchangeState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.soil_layer_count == 0 or inputs.salt_species_count == 0 or
        inputs.shared_planting_layer_index > inputs.shared_deepest_rooted_layer_index or
        inputs.shared_deepest_rooted_layer_index >= inputs.soil_layer_count or
        inputs.layer_active.len != inputs.soil_layer_count or
        inputs.root_water_volume_m3_by_layer.len != inputs.soil_layer_count or
        inputs.mycorrhizal_water_volume_m3_by_layer.len != inputs.soil_layer_count or
        inputs.minimum_pool_g_c_by_layer.len != inputs.soil_layer_count)
        return error.RootMycorrhizalExchangeDimensionMismatch;
    const scalars = [_]f64{ inputs.minimum_mycorrhiza_to_root_volume_ratio, inputs.exchange_rate_per_h, inputs.biological_timestep_h };
    for (scalars) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootMycorrhizalExchangeInput;
    for (inputs.root_water_volume_m3_by_layer) |value| if (!std.math.isFinite(value)) return error.InvalidRootMycorrhizalExchangeInput;
    for (inputs.mycorrhizal_water_volume_m3_by_layer) |value| if (!std.math.isFinite(value)) return error.InvalidRootMycorrhizalExchangeInput;
    for (inputs.minimum_pool_g_c_by_layer) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootMycorrhizalExchangeInput;
}

fn exchangeSalt(workspace: State, inputs: Inputs, species: usize, layer: usize, root_volume_m3: f64, mycorrhizal_volume_m3: f64, combined_volume_m3: f64) void {
    const root_index = compartmentIndex(inputs.soil_layer_count, species, layer, 0);
    const mycorrhizal_index = compartmentIndex(inputs.soil_layer_count, species, layer, 1);
    const root_salt_mol = @max(0.0, workspace.salt_mol_by_species_layer_compartment[root_index]);
    const mycorrhizal_salt_mol = @max(0.0, workspace.salt_mol_by_species_layer_compartment[mycorrhizal_index]);
    const concentration_difference_mol = (root_salt_mol * mycorrhizal_volume_m3 - mycorrhizal_salt_mol * root_volume_m3) / combined_volume_m3;
    const transfer_mol = inputs.exchange_rate_per_h * concentration_difference_mol * inputs.biological_timestep_h;
    workspace.salt_mol_by_species_layer_compartment[root_index] -= transfer_mol;
    workspace.salt_mol_by_species_layer_compartment[mycorrhizal_index] += transfer_mol;
}

/// Exact GROSUB 8048--8157 root--mycorrhizal mobile C:N:P and dynamic-salt
/// exchange. Traversal is the shared NU..NIX soil-layer range. Mobile pools
/// are g C/N/P, salts mol, water volumes m3, rate h-1, and timestep h.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    try validateInputs(inputs);
    copyState(workspace, state);
    if (inputs.association_status != .active) return;

    for (inputs.shared_planting_layer_index..inputs.shared_deepest_rooted_layer_index + 1) |layer| {
        if (!inputs.layer_active[layer]) continue;
        const root_carbon_g_c = @max(0.0, workspace.root_mobile_carbon_g_c_by_layer[layer]);
        const mycorrhizal_carbon_g_c = @max(0.0, workspace.mycorrhizal_mobile_carbon_g_c_by_layer[layer]);
        const root_water_volume_m3 = @max(0.0, inputs.root_water_volume_m3_by_layer[layer]);
        const effective_mycorrhizal_water_volume_m3 = @max(0.0, @min(
            inputs.root_water_volume_m3_by_layer[layer],
            @max(
                inputs.minimum_mycorrhiza_to_root_volume_ratio * inputs.root_water_volume_m3_by_layer[layer],
                inputs.mycorrhizal_water_volume_m3_by_layer[layer],
            ),
        ));
        const combined_effective_volume_m3 = root_water_volume_m3 + effective_mycorrhizal_water_volume_m3;
        if (combined_effective_volume_m3 <= inputs.minimum_pool_g_c_by_layer[layer]) continue;

        const carbon_difference_g_c = (root_carbon_g_c * effective_mycorrhizal_water_volume_m3 - mycorrhizal_carbon_g_c * root_water_volume_m3) / combined_effective_volume_m3;
        const carbon_transfer_g_c = inputs.exchange_rate_per_h * carbon_difference_g_c * inputs.biological_timestep_h;
        workspace.root_mobile_carbon_g_c_by_layer[layer] -= carbon_transfer_g_c;
        workspace.mycorrhizal_mobile_carbon_g_c_by_layer[layer] += carbon_transfer_g_c;

        const combined_mobile_carbon_g_c = root_carbon_g_c + mycorrhizal_carbon_g_c;
        if (combined_mobile_carbon_g_c > inputs.minimum_pool_g_c_by_layer[layer]) {
            const root_nitrogen_g_n = @max(0.0, workspace.root_mobile_nitrogen_g_n_by_layer[layer]);
            const mycorrhizal_nitrogen_g_n = @max(0.0, workspace.mycorrhizal_mobile_nitrogen_g_n_by_layer[layer]);
            const root_phosphorus_g_p = @max(0.0, workspace.root_mobile_phosphorus_g_p_by_layer[layer]);
            const mycorrhizal_phosphorus_g_p = @max(0.0, workspace.mycorrhizal_mobile_phosphorus_g_p_by_layer[layer]);
            const nitrogen_difference_g_n = (root_nitrogen_g_n * mycorrhizal_carbon_g_c - mycorrhizal_nitrogen_g_n * root_carbon_g_c) / combined_mobile_carbon_g_c;
            const nitrogen_transfer_g_n = inputs.exchange_rate_per_h * nitrogen_difference_g_n * inputs.biological_timestep_h;
            const phosphorus_difference_g_p = (root_phosphorus_g_p * mycorrhizal_carbon_g_c - mycorrhizal_phosphorus_g_p * root_carbon_g_c) / combined_mobile_carbon_g_c;
            const phosphorus_transfer_g_p = inputs.exchange_rate_per_h * phosphorus_difference_g_p * inputs.biological_timestep_h;
            workspace.root_mobile_nitrogen_g_n_by_layer[layer] -= nitrogen_transfer_g_n;
            workspace.mycorrhizal_mobile_nitrogen_g_n_by_layer[layer] += nitrogen_transfer_g_n;
            workspace.root_mobile_phosphorus_g_p_by_layer[layer] -= phosphorus_transfer_g_p;
            workspace.mycorrhizal_mobile_phosphorus_g_p_by_layer[layer] += phosphorus_transfer_g_p;

            if (inputs.salt_mode == .dynamic) {
                const raw_mycorrhizal_water_volume_m3 = @max(0.0, inputs.mycorrhizal_water_volume_m3_by_layer[layer]);
                const combined_raw_volume_m3 = root_water_volume_m3 + raw_mycorrhizal_water_volume_m3;
                for (0..inputs.salt_species_count) |species| exchangeSalt(workspace, inputs, species, layer, root_water_volume_m3, raw_mycorrhizal_water_volume_m3, combined_raw_volume_m3);
            }
        }
    }
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

fn stateFromArrays(root_c: []f64, myc_c: []f64, root_n: []f64, myc_n: []f64, root_p: []f64, myc_p: []f64, salts: []f64) State {
    return .{
        .root_mobile_carbon_g_c_by_layer = root_c,
        .mycorrhizal_mobile_carbon_g_c_by_layer = myc_c,
        .root_mobile_nitrogen_g_n_by_layer = root_n,
        .mycorrhizal_mobile_nitrogen_g_n_by_layer = myc_n,
        .root_mobile_phosphorus_g_p_by_layer = root_p,
        .mycorrhizal_mobile_phosphorus_g_p_by_layer = myc_p,
        .salt_mol_by_species_layer_compartment = salts,
    };
}

test "GROSUB root-mycorrhizal exchange conserves C N P and every salt" {
    var root_c = [_]f64{8.0};
    var myc_c = [_]f64{2.0};
    var root_n = [_]f64{0.8};
    var myc_n = [_]f64{0.1};
    var root_p = [_]f64{0.16};
    var myc_p = [_]f64{0.02};
    var salts = [_]f64{ 4.0, 1.0 } ** standard_salt_species_count;
    var work_root_c = [_]f64{0.0};
    var work_myc_c = [_]f64{0.0};
    var work_root_n = [_]f64{0.0};
    var work_myc_n = [_]f64{0.0};
    var work_root_p = [_]f64{0.0};
    var work_myc_p = [_]f64{0.0};
    var work_salts = [_]f64{0.0} ** (standard_salt_species_count * 2);
    const state = stateFromArrays(&root_c, &myc_c, &root_n, &myc_n, &root_p, &myc_p, &salts);
    const workspace = stateFromArrays(&work_root_c, &work_myc_c, &work_root_n, &work_myc_n, &work_root_p, &work_myc_p, &work_salts);
    const inputs: Inputs = .{
        .soil_layer_count = 1,
        .salt_species_count = standard_salt_species_count,
        .shared_planting_layer_index = 0,
        .shared_deepest_rooted_layer_index = 0,
        .layer_active = &.{true},
        .association_status = .active,
        .salt_mode = .dynamic,
        .root_water_volume_m3_by_layer = &.{2.0},
        .mycorrhizal_water_volume_m3_by_layer = &.{1.0},
        .minimum_mycorrhiza_to_root_volume_ratio = 0.1,
        .exchange_rate_per_h = 0.1,
        .biological_timestep_h = 1.0,
        .minimum_pool_g_c_by_layer = &.{1e-12},
    };
    const c_before = root_c[0] + myc_c[0];
    const n_before = root_n[0] + myc_n[0];
    const p_before = root_p[0] + myc_p[0];
    try apply(state, workspace, inputs);
    try std.testing.expectApproxEqAbs(c_before, root_c[0] + myc_c[0], 1e-14);
    try std.testing.expectApproxEqAbs(n_before, root_n[0] + myc_n[0], 1e-14);
    try std.testing.expectApproxEqAbs(p_before, root_p[0] + myc_p[0], 1e-14);
    for (0..standard_salt_species_count) |species| {
        const root_index = compartmentIndex(1, species, 0, 0);
        try std.testing.expectApproxEqAbs(5.0, salts[root_index] + salts[root_index + 1], 1e-14);
    }
}

test "GROSUB root-mycorrhizal exchange rolls back an excessive transfer" {
    var root_c = [_]f64{1.0};
    var myc_c = [_]f64{0.0};
    var root_n = [_]f64{0.1};
    var myc_n = [_]f64{0.0};
    var root_p = [_]f64{0.01};
    var myc_p = [_]f64{0.0};
    var salts = [_]f64{0.0} ** (standard_salt_species_count * 2);
    var work_root_c = [_]f64{0.0};
    var work_myc_c = [_]f64{0.0};
    var work_root_n = [_]f64{0.0};
    var work_myc_n = [_]f64{0.0};
    var work_root_p = [_]f64{0.0};
    var work_myc_p = [_]f64{0.0};
    var work_salts = [_]f64{0.0} ** (standard_salt_species_count * 2);
    const state = stateFromArrays(&root_c, &myc_c, &root_n, &myc_n, &root_p, &myc_p, &salts);
    const workspace = stateFromArrays(&work_root_c, &work_myc_c, &work_root_n, &work_myc_n, &work_root_p, &work_myc_p, &work_salts);
    const inputs: Inputs = .{
        .soil_layer_count = 1,
        .salt_species_count = standard_salt_species_count,
        .shared_planting_layer_index = 0,
        .shared_deepest_rooted_layer_index = 0,
        .layer_active = &.{true},
        .association_status = .active,
        .salt_mode = .static_equilibrium,
        .root_water_volume_m3_by_layer = &.{1.0},
        .mycorrhizal_water_volume_m3_by_layer = &.{1.0},
        .minimum_mycorrhiza_to_root_volume_ratio = 0.0,
        .exchange_rate_per_h = 3.0,
        .biological_timestep_h = 1.0,
        .minimum_pool_g_c_by_layer = &.{1e-12},
    };
    try std.testing.expectError(error.InvalidRootMycorrhizalExchangeState, apply(state, workspace, inputs));
    try std.testing.expectEqual(@as(f64, 1.0), root_c[0]);
    try std.testing.expectEqual(@as(f64, 0.0), myc_c[0]);
}
