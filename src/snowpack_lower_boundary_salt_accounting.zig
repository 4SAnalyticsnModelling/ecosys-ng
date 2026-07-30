const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

/// Snow salt species in exact REDIST update order. Snow has no runoff-only
/// hydrogen-sulfate slot.
pub const SaltSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_monohydroxide,
    aluminum_dihydroxide,
    aluminum_trihydroxide,
    aluminum_tetrahydroxide,
    aluminum_sulfate,
    iron_monohydroxide,
    iron_dihydroxide,
    iron_trihydroxide,
    iron_tetrahydroxide,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    phosphate,
    phosphoric_acid,
    iron_hydrogen_phosphate,
    iron_dihydrogen_phosphate,
    calcium_phosphate,
    calcium_hydrogen_phosphate,
    calcium_dihydrogen_phosphate,
    magnesium_hydrogen_phosphate,
};

pub const ThreePath = enum {
    surface_residue,
    soil_matrix,
    soil_macropore,
};

pub const FivePath = enum {
    surface_residue_nonband,
    soil_matrix_nonband,
    soil_macropore_nonband,
    soil_matrix_band,
    soil_macropore_band,
};

pub const salt_species_count: usize = @typeInfo(SaltSpecies).@"enum".fields.len;
pub const non_phosphorus_species_count: usize =
    @intFromEnum(SaltSpecies.phosphate);
pub const phosphorus_species_count: usize =
    salt_species_count - non_phosphorus_species_count;
pub const three_path_count: usize = @typeInfo(ThreePath).@"enum".fields.len;
pub const five_path_count: usize = @typeInfo(FivePath).@"enum".fields.len;

pub const Inputs = struct {
    /// Current snow heat capacity [snow_layer], MJ K^-1.
    heat_capacity_mj_per_k_by_layer: []const f64,
    /// Snow-layer presence threshold, MJ K^-1.
    minimum_heat_capacity_mj_per_k: f64,
    /// `X*BLS(LS)` [snow_layer][salt_species], mol per model step.
    upper_face_mol_per_step_by_layer_species: []const f64,
    /// Non-phosphorus exits [snow_layer][33 species][ThreePath], mol/step.
    three_path_ground_mol_per_step_by_layer_species_path: []const f64,
    /// Phosphorus exits [snow_layer][8 species][FivePath], mol/step.
    five_path_ground_mol_per_step_by_layer_species_path: []const f64,
};

pub const State = struct {
    /// Accumulated `T*BLS` [snow_layer][salt_species], mol per model step.
    net_mol_per_step_by_layer_species: []f64,
};

/// Caller-owned runtime scratch makes failure atomic without allocating in the
/// per-cell transport kernel.
pub const Workspace = struct {
    net_mol_per_step_by_layer_species: []f64,
};

/// Accounts dynamic salt transport at each lowest active snow/ground interface.
///
/// Traceability: REDIST.F lines 2640--2772, within the terminal active-layer
/// branch at lines 2417, 2422, and 2587. Runtime slices replace `JS`.
/// Non-phosphorus species subtract three ground exits; the eight phosphorus
/// species subtract five non-band/band exits. Flattened arrays avoid padded
/// unused paths while retaining layer, species, and source operation order.
pub fn account(
    equilibrium_mode: SaltEquilibriumMode,
    inputs: Inputs,
    state: *State,
    workspace: Workspace,
) !void {
    if (equilibrium_mode == .static) return;

    const layer_count = try validateDimensions(inputs, state.*, workspace);
    try validateInputs(inputs, state.*, workspace);
    @memcpy(
        workspace.net_mol_per_step_by_layer_species,
        state.net_mol_per_step_by_layer_species,
    );

    for (0..layer_count) |layer| {
        if (!isTerminalActiveLayer(
            inputs.heat_capacity_mj_per_k_by_layer,
            inputs.minimum_heat_capacity_mj_per_k,
            layer,
        )) continue;
        try updateTerminalLayer(inputs, workspace, layer);
    }

    @memcpy(
        state.net_mol_per_step_by_layer_species,
        workspace.net_mol_per_step_by_layer_species,
    );
}

fn isTerminalActiveLayer(
    capacities: []const f64,
    minimum_capacity: f64,
    layer: usize,
) bool {
    if (capacities[layer] <= minimum_capacity) return false;
    return layer + 1 == capacities.len or
        capacities[layer + 1] <= minimum_capacity;
}

fn updateTerminalLayer(
    inputs: Inputs,
    workspace: Workspace,
    layer: usize,
) !void {
    const salt_offset = layer * salt_species_count;
    const three_path_offset =
        layer * non_phosphorus_species_count * three_path_count;
    for (0..non_phosphorus_species_count) |species| {
        const salt_index = salt_offset + species;
        var candidate = try add(
            workspace.net_mol_per_step_by_layer_species[salt_index],
            inputs.upper_face_mol_per_step_by_layer_species[salt_index],
        );
        const path_offset = three_path_offset + species * three_path_count;
        for (0..three_path_count) |path| {
            candidate = try subtract(
                candidate,
                inputs.three_path_ground_mol_per_step_by_layer_species_path[
                    path_offset + path
                ],
            );
        }
        workspace.net_mol_per_step_by_layer_species[salt_index] = candidate;
    }

    const five_path_offset =
        layer * phosphorus_species_count * five_path_count;
    for (0..phosphorus_species_count) |phosphorus_species| {
        const salt_index =
            salt_offset + non_phosphorus_species_count + phosphorus_species;
        var candidate = try add(
            workspace.net_mol_per_step_by_layer_species[salt_index],
            inputs.upper_face_mol_per_step_by_layer_species[salt_index],
        );
        const path_offset =
            five_path_offset + phosphorus_species * five_path_count;
        for (0..five_path_count) |path| {
            candidate = try subtract(
                candidate,
                inputs.five_path_ground_mol_per_step_by_layer_species_path[
                    path_offset + path
                ],
            );
        }
        workspace.net_mol_per_step_by_layer_species[salt_index] = candidate;
    }
}

fn validateDimensions(inputs: Inputs, state: State, workspace: Workspace) !usize {
    const layer_count = inputs.heat_capacity_mj_per_k_by_layer.len;
    if (layer_count == 0)
        return error.InvalidSnowpackBoundarySaltDimensions;
    const salt_value_count = try checkedExtent(layer_count, salt_species_count);
    const three_path_value_count = try checkedExtent(
        try checkedExtent(layer_count, non_phosphorus_species_count),
        three_path_count,
    );
    const five_path_value_count = try checkedExtent(
        try checkedExtent(layer_count, phosphorus_species_count),
        five_path_count,
    );
    if (inputs.upper_face_mol_per_step_by_layer_species.len != salt_value_count or
        inputs.three_path_ground_mol_per_step_by_layer_species_path.len !=
            three_path_value_count or
        inputs.five_path_ground_mol_per_step_by_layer_species_path.len !=
            five_path_value_count or
        state.net_mol_per_step_by_layer_species.len != salt_value_count or
        workspace.net_mol_per_step_by_layer_species.len != salt_value_count)
    {
        return error.SnowpackBoundarySaltDimensionMismatch;
    }
    return layer_count;
}

fn checkedExtent(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        error.InvalidSnowpackBoundarySaltDimensions;
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (!std.math.isFinite(inputs.minimum_heat_capacity_mj_per_k))
        return error.NonFiniteSnowpackBoundarySaltInput;
    if (inputs.minimum_heat_capacity_mj_per_k < 0)
        return error.InvalidSnowpackHeatCapacity;
    for (inputs.heat_capacity_mj_per_k_by_layer) |capacity| {
        if (!std.math.isFinite(capacity))
            return error.NonFiniteSnowpackBoundarySaltInput;
        if (capacity < 0) return error.InvalidSnowpackHeatCapacity;
    }
    try validateFinite(inputs.upper_face_mol_per_step_by_layer_species);
    try validateFinite(
        inputs.three_path_ground_mol_per_step_by_layer_species_path,
    );
    try validateFinite(
        inputs.five_path_ground_mol_per_step_by_layer_species_path,
    );
    try validateFinite(state.net_mol_per_step_by_layer_species);

    const scratch = workspace.net_mol_per_step_by_layer_species;
    if (overlap(scratch, state.net_mol_per_step_by_layer_species) or
        overlap(scratch, inputs.upper_face_mol_per_step_by_layer_species) or
        overlap(
            scratch,
            inputs.three_path_ground_mol_per_step_by_layer_species_path,
        ) or
        overlap(
            scratch,
            inputs.five_path_ground_mol_per_step_by_layer_species_path,
        ))
    {
        return error.SnowpackBoundarySaltWorkspaceOverlap;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteSnowpackBoundarySaltInput;
}

fn add(current: f64, contribution: f64) !f64 {
    const result = current + contribution;
    if (!std.math.isFinite(result))
        return error.NonFiniteSnowpackBoundarySaltResult;
    return result;
}

fn subtract(current: f64, contribution: f64) !f64 {
    const result = current - contribution;
    if (!std.math.isFinite(result))
        return error.NonFiniteSnowpackBoundarySaltResult;
    return result;
}

fn overlap(left: []const f64, right: []const f64) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(f64);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn saltLayer(values: []f64, layer: usize) []f64 {
    const start = layer * salt_species_count;
    return values[start .. start + salt_species_count];
}

fn saltLayerConst(values: []const f64, layer: usize) []const f64 {
    const start = layer * salt_species_count;
    return values[start .. start + salt_species_count];
}

fn expectSpeciesGroups(
    values: []const f64,
    layer: usize,
    non_phosphorus_expected: f64,
    phosphorus_expected: f64,
) !void {
    const layer_values = saltLayerConst(values, layer);
    for (layer_values[0..non_phosphorus_species_count]) |value|
        try std.testing.expectEqual(non_phosphorus_expected, value);
    for (layer_values[non_phosphorus_species_count..]) |value|
        try std.testing.expectEqual(phosphorus_expected, value);
}

test "dynamic terminal salt accounting retains exact path groups" {
    const layer_count = 1;
    const capacities = [_]f64{2};
    const upper = [_]f64{10} ** salt_species_count;
    const three_path_ground =
        [_]f64{2} ** (non_phosphorus_species_count * three_path_count);
    const five_path_ground =
        [_]f64{3} ** (phosphorus_species_count * five_path_count);
    var totals = [_]f64{100} ** salt_species_count;
    var scratch: [salt_species_count]f64 = undefined;
    var state = State{ .net_mol_per_step_by_layer_species = &totals };

    try account(.dynamic, .{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_mol_per_step_by_layer_species = &upper,
        .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
        .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
    }, &state, .{ .net_mol_per_step_by_layer_species = &scratch });

    try expectSpeciesGroups(&totals, 0, 104, 95);
    try std.testing.expectEqual(@as(usize, 1), layer_count);
}

test "snow change plus explicit ground paths conserves every salt exactly" {
    const capacities = [_]f64{2};
    const upper = [_]f64{20} ** salt_species_count;
    var three_path_ground =
        [_]f64{0} ** (non_phosphorus_species_count * three_path_count);
    for (0..non_phosphorus_species_count) |species| {
        for (0..three_path_count) |path|
            three_path_ground[species * three_path_count + path] =
                @as(f64, @floatFromInt(path + 1));
    }
    var five_path_ground =
        [_]f64{0} ** (phosphorus_species_count * five_path_count);
    for (0..phosphorus_species_count) |species| {
        for (0..five_path_count) |path|
            five_path_ground[species * five_path_count + path] =
                @as(f64, @floatFromInt(path + 1));
    }
    var totals = [_]f64{0} ** salt_species_count;
    var scratch: [salt_species_count]f64 = undefined;
    var state = State{ .net_mol_per_step_by_layer_species = &totals };

    try account(.dynamic, .{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_mol_per_step_by_layer_species = &upper,
        .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
        .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
    }, &state, .{ .net_mol_per_step_by_layer_species = &scratch });

    try expectSpeciesGroups(&totals, 0, 14, 5);
    for (totals[0..non_phosphorus_species_count]) |stored|
        try std.testing.expectEqual(@as(f64, 20), stored + 1 + 2 + 3);
    for (totals[non_phosphorus_species_count..]) |stored|
        try std.testing.expectEqual(@as(f64, 20), stored + 1 + 2 + 3 + 4 + 5);
}

test "species path layouts and runtime terminal layers are explicit" {
    try std.testing.expectEqual(@as(usize, 41), salt_species_count);
    try std.testing.expectEqual(@as(usize, 33), non_phosphorus_species_count);
    try std.testing.expectEqual(@as(usize, 8), phosphorus_species_count);
    try std.testing.expectEqual(@as(usize, 3), three_path_count);
    try std.testing.expectEqual(@as(usize, 5), five_path_count);
    try std.testing.expectEqual(
        @as(usize, 33),
        @intFromEnum(SaltSpecies.phosphate),
    );

    const layer_count = 7;
    const capacities = [_]f64{ 2, 2, 0, 2, 2, 2, 2 };
    var upper = [_]f64{0} ** (layer_count * salt_species_count);
    for (0..layer_count) |layer|
        @memset(saltLayer(&upper, layer), @as(f64, @floatFromInt(layer + 1)));
    const three_path_ground = [_]f64{0} **
        (layer_count * non_phosphorus_species_count * three_path_count);
    const five_path_ground = [_]f64{0} **
        (layer_count * phosphorus_species_count * five_path_count);
    var totals = [_]f64{0} ** (layer_count * salt_species_count);
    var scratch: [layer_count * salt_species_count]f64 = undefined;
    var state = State{ .net_mol_per_step_by_layer_species = &totals };

    try account(.dynamic, .{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_mol_per_step_by_layer_species = &upper,
        .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
        .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
    }, &state, .{ .net_mol_per_step_by_layer_species = &scratch });

    try expectSpeciesGroups(&totals, 0, 0, 0);
    try expectSpeciesGroups(&totals, 1, 2, 2);
    try expectSpeciesGroups(&totals, 2, 0, 0);
    try expectSpeciesGroups(&totals, 3, 0, 0);
    try expectSpeciesGroups(&totals, 6, 7, 7);
}

test "static mode bypasses dynamic dimensions and mutation" {
    var totals = [_]f64{11};
    var scratch = [_]f64{13};
    var state = State{ .net_mol_per_step_by_layer_species = &totals };
    try account(.static, .{
        .heat_capacity_mj_per_k_by_layer = &.{},
        .minimum_heat_capacity_mj_per_k = std.math.nan(f64),
        .upper_face_mol_per_step_by_layer_species = &.{},
        .three_path_ground_mol_per_step_by_layer_species_path = &.{},
        .five_path_ground_mol_per_step_by_layer_species_path = &.{},
    }, &state, .{ .net_mol_per_step_by_layer_species = &scratch });
    try std.testing.expectEqual(@as(f64, 11), totals[0]);
    try std.testing.expectEqual(@as(f64, 13), scratch[0]);
}

test "source association and atomic failures are retained" {
    const capacities = [_]f64{2};
    const upper = [_]f64{-1.0e16} ** salt_species_count;
    var three_path_ground =
        [_]f64{0} ** (non_phosphorus_species_count * three_path_count);
    for (0..non_phosphorus_species_count) |species|
        three_path_ground[species * three_path_count] = 1;
    var five_path_ground =
        [_]f64{0} ** (phosphorus_species_count * five_path_count);
    for (0..phosphorus_species_count) |species|
        five_path_ground[species * five_path_count] = 1;
    var totals = [_]f64{1.0e16} ** salt_species_count;
    var scratch: [salt_species_count]f64 = undefined;
    var state = State{ .net_mol_per_step_by_layer_species = &totals };
    const workspace = Workspace{ .net_mol_per_step_by_layer_species = &scratch };

    try account(.dynamic, .{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_mol_per_step_by_layer_species = &upper,
        .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
        .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
    }, &state, workspace);
    try expectSpeciesGroups(&totals, 0, -1, -1);

    @memset(&totals, 3);
    five_path_ground[five_path_ground.len - 1] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSnowpackBoundarySaltInput,
        account(.dynamic, .{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_mol_per_step_by_layer_species = &upper,
            .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
            .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
        }, &state, workspace),
    );
    try expectSpeciesGroups(&totals, 0, 3, 3);

    five_path_ground[five_path_ground.len - 1] = 0;
    const overflowing_upper =
        [_]f64{std.math.floatMax(f64)} ** salt_species_count;
    @memset(&totals, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteSnowpackBoundarySaltResult,
        account(.dynamic, .{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_mol_per_step_by_layer_species = &overflowing_upper,
            .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
            .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
        }, &state, workspace),
    );
    try expectSpeciesGroups(
        &totals,
        0,
        std.math.floatMax(f64),
        std.math.floatMax(f64),
    );
}

test "dimension and workspace alias failures precede mutation" {
    const capacities = [_]f64{2};
    const upper = [_]f64{1} ** salt_species_count;
    const three_path_ground =
        [_]f64{0} ** (non_phosphorus_species_count * three_path_count);
    const five_path_ground =
        [_]f64{0} ** (phosphorus_species_count * five_path_count);
    var totals = [_]f64{5} ** salt_species_count;
    var state = State{ .net_mol_per_step_by_layer_species = &totals };

    try std.testing.expectError(
        error.SnowpackBoundarySaltWorkspaceOverlap,
        account(.dynamic, .{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_mol_per_step_by_layer_species = &upper,
            .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
            .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
        }, &state, .{ .net_mol_per_step_by_layer_species = &totals }),
    );
    try expectSpeciesGroups(&totals, 0, 5, 5);

    var short_scratch: [salt_species_count - 1]f64 = undefined;
    try std.testing.expectError(
        error.SnowpackBoundarySaltDimensionMismatch,
        account(.dynamic, .{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_mol_per_step_by_layer_species = &upper,
            .three_path_ground_mol_per_step_by_layer_species_path = &three_path_ground,
            .five_path_ground_mol_per_step_by_layer_species_path = &five_path_ground,
        }, &state, .{ .net_mol_per_step_by_layer_species = &short_scratch }),
    );
    try expectSpeciesGroups(&totals, 0, 5, 5);
}
