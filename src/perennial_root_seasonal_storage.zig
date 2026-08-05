const std = @import("std");

pub const GrowthHabit = enum { annual, perennial };
pub const StorageTransferStatus = enum { disabled, enabled };

pub const State = struct {
    root_mobile_carbon_g_c_by_species_layer: []f64,
    root_mobile_nitrogen_g_n_by_species_layer: []f64,
    root_mobile_phosphorus_g_p_by_species_layer: []f64,
    seasonal_storage_carbon_g_c_by_species: []f64,
    seasonal_storage_nitrogen_g_n_by_species: []f64,
    seasonal_storage_phosphorus_g_p_by_species: []f64,
};

pub const Inputs = struct {
    plant_species_count: usize,
    soil_layer_count: usize,
    shared_planting_layer_index: usize,
    deepest_rooted_layer_index_by_species: []const usize,
    layer_active_by_species_layer: []const bool,
    growth_habit_by_species: []const GrowthHabit,
    storage_transfer_status_by_species: []const StorageTransferStatus,
    exchange_rate_h_inv_by_species: []const f64,
    biological_timestep_h: f64,
    minimum_nitrogen_carbon_ratio_g_n_g_c: f64,
    maximum_nitrogen_carbon_ratio_g_n_g_c: f64,
    minimum_phosphorus_carbon_ratio_g_p_g_c: f64,
    maximum_phosphorus_carbon_ratio_g_p_g_c: f64,
};

fn poolIndex(inputs: Inputs, species: usize, layer: usize) usize {
    return species * inputs.soil_layer_count + layer;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, inputs: Inputs) !void {
    const pool_count = std.math.mul(usize, inputs.plant_species_count, inputs.soil_layer_count) catch return error.SeasonalStorageDimensionOverflow;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const expected = if (std.mem.indexOf(u8, field.name, "species_layer") != null) pool_count else inputs.plant_species_count;
        if (values.len != expected) return error.SeasonalStorageDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidSeasonalStorageState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    const pool_count = std.math.mul(usize, inputs.plant_species_count, inputs.soil_layer_count) catch return error.SeasonalStorageDimensionOverflow;
    if (inputs.plant_species_count == 0 or inputs.soil_layer_count == 0 or
        inputs.shared_planting_layer_index >= inputs.soil_layer_count or
        inputs.deepest_rooted_layer_index_by_species.len != inputs.plant_species_count or
        inputs.layer_active_by_species_layer.len != pool_count or
        inputs.growth_habit_by_species.len != inputs.plant_species_count or
        inputs.storage_transfer_status_by_species.len != inputs.plant_species_count or
        inputs.exchange_rate_h_inv_by_species.len != inputs.plant_species_count)
        return error.SeasonalStorageDimensionMismatch;
    for (inputs.deepest_rooted_layer_index_by_species) |layer| if (layer < inputs.shared_planting_layer_index or layer >= inputs.soil_layer_count) return error.SeasonalStorageDimensionMismatch;
    for (inputs.exchange_rate_h_inv_by_species) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidSeasonalStorageInput;
    const positive_scalars = [_]f64{
        inputs.biological_timestep_h,
        inputs.minimum_nitrogen_carbon_ratio_g_n_g_c,
        inputs.maximum_nitrogen_carbon_ratio_g_n_g_c,
        inputs.minimum_phosphorus_carbon_ratio_g_p_g_c,
        inputs.maximum_phosphorus_carbon_ratio_g_p_g_c,
    };
    for (positive_scalars) |value| if (!std.math.isFinite(value) or value <= 0.0) return error.InvalidSeasonalStorageInput;
}

/// Exact GROSUB 8174--8199 perennial root mobile-pool transfer to seasonal
/// storage. Runtime topology is [plant species][soil layer]. Pools are g C/N/P,
/// ratios are g C per g N/P, exchange rates h-1, and timestep h.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    copyState(workspace, state);

    for (0..inputs.plant_species_count) |species| {
        if (inputs.storage_transfer_status_by_species[species] != .enabled or inputs.growth_habit_by_species[species] != .perennial) continue;
        for (inputs.shared_planting_layer_index..inputs.deepest_rooted_layer_index_by_species[species] + 1) |layer| {
            const index = poolIndex(inputs, species, layer);
            if (!inputs.layer_active_by_species_layer[index]) continue;
            const unconstrained_carbon_g_c = inputs.exchange_rate_h_inv_by_species[species] * @max(0.0, workspace.root_mobile_carbon_g_c_by_species_layer[index]) * inputs.biological_timestep_h;
            const unconstrained_nitrogen_g_n = inputs.exchange_rate_h_inv_by_species[species] * @max(0.0, workspace.root_mobile_nitrogen_g_n_by_species_layer[index]) * inputs.biological_timestep_h;
            const unconstrained_phosphorus_g_p = inputs.exchange_rate_h_inv_by_species[species] * @max(0.0, workspace.root_mobile_phosphorus_g_p_by_species_layer[index]) * inputs.biological_timestep_h;
            const carbon_transfer_g_c = @min(
                unconstrained_carbon_g_c,
                unconstrained_nitrogen_g_n / inputs.minimum_nitrogen_carbon_ratio_g_n_g_c,
                unconstrained_phosphorus_g_p / inputs.minimum_phosphorus_carbon_ratio_g_p_g_c,
            );
            const nitrogen_transfer_g_n = @min(
                unconstrained_nitrogen_g_n,
                carbon_transfer_g_c * inputs.maximum_nitrogen_carbon_ratio_g_n_g_c,
                unconstrained_phosphorus_g_p * inputs.maximum_nitrogen_carbon_ratio_g_n_g_c / inputs.minimum_phosphorus_carbon_ratio_g_p_g_c,
            );
            const phosphorus_transfer_g_p = @min(
                unconstrained_phosphorus_g_p,
                carbon_transfer_g_c * inputs.maximum_phosphorus_carbon_ratio_g_p_g_c,
                unconstrained_nitrogen_g_n * inputs.maximum_phosphorus_carbon_ratio_g_p_g_c / inputs.minimum_nitrogen_carbon_ratio_g_n_g_c,
            );
            workspace.root_mobile_carbon_g_c_by_species_layer[index] -= carbon_transfer_g_c;
            workspace.seasonal_storage_carbon_g_c_by_species[species] += carbon_transfer_g_c;
            workspace.root_mobile_nitrogen_g_n_by_species_layer[index] -= nitrogen_transfer_g_n;
            workspace.seasonal_storage_nitrogen_g_n_by_species[species] += nitrogen_transfer_g_n;
            workspace.root_mobile_phosphorus_g_p_by_species_layer[index] -= phosphorus_transfer_g_p;
            workspace.seasonal_storage_phosphorus_g_p_by_species[species] += phosphorus_transfer_g_p;
        }
    }
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

test "GROSUB perennial seasonal-storage transfer conserves C N and P" {
    var root_c = [_]f64{ 10.0, 6.0, 9.0, 7.0 };
    var root_n = [_]f64{ 2.0, 1.0, 2.0, 1.0 };
    var root_p = [_]f64{ 1.0, 0.5, 1.0, 0.5 };
    var storage_c = [_]f64{ 1.0, 2.0 };
    var storage_n = [_]f64{ 0.1, 0.2 };
    var storage_p = [_]f64{ 0.01, 0.02 };
    var work_root_c = [_]f64{0} ** 4;
    var work_root_n = [_]f64{0} ** 4;
    var work_root_p = [_]f64{0} ** 4;
    var work_storage_c = [_]f64{0} ** 2;
    var work_storage_n = [_]f64{0} ** 2;
    var work_storage_p = [_]f64{0} ** 2;
    const state: State = .{
        .root_mobile_carbon_g_c_by_species_layer = &root_c,
        .root_mobile_nitrogen_g_n_by_species_layer = &root_n,
        .root_mobile_phosphorus_g_p_by_species_layer = &root_p,
        .seasonal_storage_carbon_g_c_by_species = &storage_c,
        .seasonal_storage_nitrogen_g_n_by_species = &storage_n,
        .seasonal_storage_phosphorus_g_p_by_species = &storage_p,
    };
    const workspace: State = .{
        .root_mobile_carbon_g_c_by_species_layer = &work_root_c,
        .root_mobile_nitrogen_g_n_by_species_layer = &work_root_n,
        .root_mobile_phosphorus_g_p_by_species_layer = &work_root_p,
        .seasonal_storage_carbon_g_c_by_species = &work_storage_c,
        .seasonal_storage_nitrogen_g_n_by_species = &work_storage_n,
        .seasonal_storage_phosphorus_g_p_by_species = &work_storage_p,
    };
    const inputs: Inputs = .{
        .plant_species_count = 2,
        .soil_layer_count = 2,
        .shared_planting_layer_index = 0,
        .deepest_rooted_layer_index_by_species = &.{ 1, 1 },
        .layer_active_by_species_layer = &.{ true, true, true, true },
        .growth_habit_by_species = &.{ .perennial, .annual },
        .storage_transfer_status_by_species = &.{ .enabled, .enabled },
        .exchange_rate_h_inv_by_species = &.{ 0.1, 0.1 },
        .biological_timestep_h = 1.0,
        .minimum_nitrogen_carbon_ratio_g_n_g_c = 0.05,
        .maximum_nitrogen_carbon_ratio_g_n_g_c = 0.2,
        .minimum_phosphorus_carbon_ratio_g_p_g_c = 0.005,
        .maximum_phosphorus_carbon_ratio_g_p_g_c = 0.02,
    };
    const before_c = root_c[0] + root_c[1] + storage_c[0];
    const before_n = root_n[0] + root_n[1] + storage_n[0];
    const before_p = root_p[0] + root_p[1] + storage_p[0];
    try apply(state, workspace, inputs);
    try std.testing.expectApproxEqAbs(before_c, root_c[0] + root_c[1] + storage_c[0], 1e-14);
    try std.testing.expectApproxEqAbs(before_n, root_n[0] + root_n[1] + storage_n[0], 1e-14);
    try std.testing.expectApproxEqAbs(before_p, root_p[0] + root_p[1] + storage_p[0], 1e-14);
    try std.testing.expectEqual(@as(f64, 9.0), root_c[2]);
    try std.testing.expectEqual(@as(f64, 2.0), storage_c[1]);
}
