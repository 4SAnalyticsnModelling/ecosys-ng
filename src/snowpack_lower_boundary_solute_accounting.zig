const std = @import("std");

/// Non-salt snow solutes on their elemental gram basis per model step.
pub const SoluteFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    hydrogen_phosphate_g_p_per_step: f64 = 0,
    dihydrogen_phosphate_g_p_per_step: f64 = 0,
};

/// Ground exits used by gaseous carriers in REDIST source order.
pub const ThreePathTransfer = struct {
    to_surface_residue_g_per_step: f64 = 0,
    to_soil_matrix_g_per_step: f64 = 0,
    to_soil_macropore_g_per_step: f64 = 0,
};

/// Ground exits used by nutrient and phosphate carriers in source order.
pub const FivePathTransfer = struct {
    to_surface_residue_nonband_g_per_step: f64 = 0,
    to_soil_matrix_nonband_g_per_step: f64 = 0,
    to_soil_macropore_nonband_g_per_step: f64 = 0,
    to_soil_matrix_band_g_per_step: f64 = 0,
    to_soil_macropore_band_g_per_step: f64 = 0,
};

pub const GroundTransfer = struct {
    carbon_dioxide: ThreePathTransfer = .{},
    methane: ThreePathTransfer = .{},
    oxygen: ThreePathTransfer = .{},
    dinitrogen: ThreePathTransfer = .{},
    nitrous_oxide: ThreePathTransfer = .{},
    ammonium: FivePathTransfer = .{},
    ammonia: FivePathTransfer = .{},
    nitrate: FivePathTransfer = .{},
    hydrogen_phosphate: FivePathTransfer = .{},
    dihydrogen_phosphate: FivePathTransfer = .{},
};

pub const Inputs = struct {
    /// Current snow heat capacity [snow_layer], MJ K^-1.
    heat_capacity_mj_per_k_by_layer: []const f64,
    /// Snow-layer presence threshold, MJ K^-1.
    minimum_heat_capacity_mj_per_k: f64,
    /// `X*BLS(LS)` entering each layer [snow_layer].
    upper_face_flux_by_layer: []const SoluteFlux,
    /// Explicit residue/matrix/macropore exits [snow_layer].
    ground_transfer_by_layer: []const GroundTransfer,
};

pub const State = struct {
    /// Accumulated `T*BLS` [snow_layer].
    net_flux_by_layer: []SoluteFlux,
};

/// Caller-owned runtime scratch provides atomic commit without kernel
/// allocation.
pub const Workspace = struct {
    net_flux_by_layer: []SoluteFlux,
};

/// Accounts non-salt solutes at each lowest active snow/ground interface.
///
/// Traceability: REDIST.F lines 2605--2639, within the active-current and
/// inactive-lower-layer branch at lines 2417, 2422, and 2587. Runtime slices
/// replace `JS`. Gaseous carriers subtract residue, soil-matrix, then
/// soil-macropore transfers. Nutrients and phosphates additionally retain
/// separate non-band and band matrix/macropore paths. Every scalar operation
/// remains in source order and state commits only after finite evaluation.
pub fn account(inputs: Inputs, state: *State, workspace: Workspace) !void {
    try validateDimensions(inputs, state.*, workspace);
    try validateInputs(inputs, state.*, workspace);
    @memcpy(workspace.net_flux_by_layer, state.net_flux_by_layer);

    const layer_count = inputs.heat_capacity_mj_per_k_by_layer.len;
    for (0..layer_count) |layer| {
        if (inputs.heat_capacity_mj_per_k_by_layer[layer] <=
            inputs.minimum_heat_capacity_mj_per_k)
        {
            continue;
        }
        if (layer + 1 < layer_count and
            inputs.heat_capacity_mj_per_k_by_layer[layer + 1] >
                inputs.minimum_heat_capacity_mj_per_k)
        {
            continue;
        }
        try updateBoundary(
            &workspace.net_flux_by_layer[layer],
            inputs.upper_face_flux_by_layer[layer],
            inputs.ground_transfer_by_layer[layer],
        );
    }

    @memcpy(state.net_flux_by_layer, workspace.net_flux_by_layer);
}

fn validateDimensions(inputs: Inputs, state: State, workspace: Workspace) !void {
    const layer_count = inputs.heat_capacity_mj_per_k_by_layer.len;
    if (layer_count == 0)
        return error.InvalidSnowpackBoundarySoluteDimensions;
    if (inputs.upper_face_flux_by_layer.len != layer_count or
        inputs.ground_transfer_by_layer.len != layer_count or
        state.net_flux_by_layer.len != layer_count or
        workspace.net_flux_by_layer.len != layer_count)
    {
        return error.SnowpackBoundarySoluteDimensionMismatch;
    }
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (!std.math.isFinite(inputs.minimum_heat_capacity_mj_per_k))
        return error.NonFiniteSnowpackBoundarySoluteInput;
    if (inputs.minimum_heat_capacity_mj_per_k < 0)
        return error.InvalidSnowpackHeatCapacity;
    for (inputs.heat_capacity_mj_per_k_by_layer) |capacity| {
        if (!std.math.isFinite(capacity))
            return error.NonFiniteSnowpackBoundarySoluteInput;
        if (capacity < 0) return error.InvalidSnowpackHeatCapacity;
    }
    try validateSoluteFluxes(inputs.upper_face_flux_by_layer);
    try validateGroundTransfers(inputs.ground_transfer_by_layer);
    try validateSoluteFluxes(state.net_flux_by_layer);
    if (overlap(workspace.net_flux_by_layer, state.net_flux_by_layer) or
        overlap(workspace.net_flux_by_layer, inputs.upper_face_flux_by_layer))
    {
        return error.SnowpackBoundarySoluteWorkspaceOverlap;
    }
}

fn validateSoluteFluxes(values: []const SoluteFlux) !void {
    for (values) |value| {
        inline for (@typeInfo(SoluteFlux).@"struct".fields) |field|
            if (!std.math.isFinite(@field(value, field.name)))
                return error.NonFiniteSnowpackBoundarySoluteInput;
    }
}

fn validateGroundTransfers(values: []const GroundTransfer) !void {
    for (values) |value| {
        inline for (@typeInfo(GroundTransfer).@"struct".fields) |outer_field| {
            const pathways = @field(value, outer_field.name);
            inline for (@typeInfo(@TypeOf(pathways)).@"struct".fields) |field|
                if (!std.math.isFinite(@field(pathways, field.name)))
                    return error.NonFiniteSnowpackBoundarySoluteInput;
        }
    }
}

fn updateBoundary(
    candidate: *SoluteFlux,
    upper: SoluteFlux,
    ground: GroundTransfer,
) !void {
    try updateThree(
        &candidate.carbon_dioxide_g_c_per_step,
        upper.carbon_dioxide_g_c_per_step,
        ground.carbon_dioxide,
    );
    try updateThree(
        &candidate.methane_g_c_per_step,
        upper.methane_g_c_per_step,
        ground.methane,
    );
    try updateThree(
        &candidate.oxygen_g_o_per_step,
        upper.oxygen_g_o_per_step,
        ground.oxygen,
    );
    try updateThree(
        &candidate.dinitrogen_g_n_per_step,
        upper.dinitrogen_g_n_per_step,
        ground.dinitrogen,
    );
    try updateThree(
        &candidate.nitrous_oxide_g_n_per_step,
        upper.nitrous_oxide_g_n_per_step,
        ground.nitrous_oxide,
    );
    try updateFive(
        &candidate.ammonium_g_n_per_step,
        upper.ammonium_g_n_per_step,
        ground.ammonium,
    );
    try updateFive(
        &candidate.ammonia_g_n_per_step,
        upper.ammonia_g_n_per_step,
        ground.ammonia,
    );
    try updateFive(
        &candidate.nitrate_g_n_per_step,
        upper.nitrate_g_n_per_step,
        ground.nitrate,
    );
    try updateFive(
        &candidate.hydrogen_phosphate_g_p_per_step,
        upper.hydrogen_phosphate_g_p_per_step,
        ground.hydrogen_phosphate,
    );
    try updateFive(
        &candidate.dihydrogen_phosphate_g_p_per_step,
        upper.dihydrogen_phosphate_g_p_per_step,
        ground.dihydrogen_phosphate,
    );
}

fn updateThree(
    candidate: *f64,
    upper: f64,
    ground: ThreePathTransfer,
) !void {
    candidate.* = try add(candidate.*, upper);
    candidate.* = try subtract(candidate.*, ground.to_surface_residue_g_per_step);
    candidate.* = try subtract(candidate.*, ground.to_soil_matrix_g_per_step);
    candidate.* = try subtract(candidate.*, ground.to_soil_macropore_g_per_step);
}

fn updateFive(
    candidate: *f64,
    upper: f64,
    ground: FivePathTransfer,
) !void {
    candidate.* = try add(candidate.*, upper);
    candidate.* = try subtract(
        candidate.*,
        ground.to_surface_residue_nonband_g_per_step,
    );
    candidate.* = try subtract(
        candidate.*,
        ground.to_soil_matrix_nonband_g_per_step,
    );
    candidate.* = try subtract(
        candidate.*,
        ground.to_soil_macropore_nonband_g_per_step,
    );
    candidate.* = try subtract(
        candidate.*,
        ground.to_soil_matrix_band_g_per_step,
    );
    candidate.* = try subtract(
        candidate.*,
        ground.to_soil_macropore_band_g_per_step,
    );
}

fn add(current: f64, contribution: f64) !f64 {
    const result = current + contribution;
    if (!std.math.isFinite(result))
        return error.NonFiniteSnowpackBoundarySoluteResult;
    return result;
}

fn subtract(current: f64, contribution: f64) !f64 {
    const result = current - contribution;
    if (!std.math.isFinite(result))
        return error.NonFiniteSnowpackBoundarySoluteResult;
    return result;
}

fn overlap(left: []const SoluteFlux, right: []const SoluteFlux) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(SoluteFlux);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn filledSolute(value: f64) SoluteFlux {
    var result: SoluteFlux = undefined;
    inline for (@typeInfo(SoluteFlux).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn filledThree(value: f64) ThreePathTransfer {
    return .{
        .to_surface_residue_g_per_step = value,
        .to_soil_matrix_g_per_step = value,
        .to_soil_macropore_g_per_step = value,
    };
}

fn filledFive(value: f64) FivePathTransfer {
    return .{
        .to_surface_residue_nonband_g_per_step = value,
        .to_soil_matrix_nonband_g_per_step = value,
        .to_soil_macropore_nonband_g_per_step = value,
        .to_soil_matrix_band_g_per_step = value,
        .to_soil_macropore_band_g_per_step = value,
    };
}

fn filledGround(three_value: f64, five_value: f64) GroundTransfer {
    return .{
        .carbon_dioxide = filledThree(three_value),
        .methane = filledThree(three_value),
        .oxygen = filledThree(three_value),
        .dinitrogen = filledThree(three_value),
        .nitrous_oxide = filledThree(three_value),
        .ammonium = filledFive(five_value),
        .ammonia = filledFive(five_value),
        .nitrate = filledFive(five_value),
        .hydrogen_phosphate = filledFive(five_value),
        .dihydrogen_phosphate = filledFive(five_value),
    };
}

fn expectGasValues(actual: SoluteFlux, expected: f64) !void {
    try std.testing.expectEqual(expected, actual.carbon_dioxide_g_c_per_step);
    try std.testing.expectEqual(expected, actual.methane_g_c_per_step);
    try std.testing.expectEqual(expected, actual.oxygen_g_o_per_step);
    try std.testing.expectEqual(expected, actual.dinitrogen_g_n_per_step);
    try std.testing.expectEqual(expected, actual.nitrous_oxide_g_n_per_step);
}

fn expectNutrientValues(actual: SoluteFlux, expected: f64) !void {
    try std.testing.expectEqual(expected, actual.ammonium_g_n_per_step);
    try std.testing.expectEqual(expected, actual.ammonia_g_n_per_step);
    try std.testing.expectEqual(expected, actual.nitrate_g_n_per_step);
    try std.testing.expectEqual(expected, actual.hydrogen_phosphate_g_p_per_step);
    try std.testing.expectEqual(expected, actual.dihydrogen_phosphate_g_p_per_step);
}

test "lowest active layer retains three and five pathway source order" {
    const capacities = [_]f64{2};
    const upper = [_]SoluteFlux{filledSolute(10)};
    const ground = [_]GroundTransfer{filledGround(2, 3)};
    var totals = [_]SoluteFlux{filledSolute(100)};
    var scratch: [1]SoluteFlux = undefined;
    var state = State{ .net_flux_by_layer = &totals };

    try account(.{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_flux_by_layer = &upper,
        .ground_transfer_by_layer = &ground,
    }, &state, .{ .net_flux_by_layer = &scratch });

    try expectGasValues(totals[0], 104);
    try expectNutrientValues(totals[0], 95);
}

test "snow change plus all ground exports conserves every carrier exactly" {
    const capacities = [_]f64{2};
    const upper = [_]SoluteFlux{filledSolute(20)};
    const ground = [_]GroundTransfer{.{
        .carbon_dioxide = .{
            .to_surface_residue_g_per_step = 1,
            .to_soil_matrix_g_per_step = 2,
            .to_soil_macropore_g_per_step = 3,
        },
        .methane = .{
            .to_surface_residue_g_per_step = 1,
            .to_soil_matrix_g_per_step = 2,
            .to_soil_macropore_g_per_step = 3,
        },
        .oxygen = .{
            .to_surface_residue_g_per_step = 1,
            .to_soil_matrix_g_per_step = 2,
            .to_soil_macropore_g_per_step = 3,
        },
        .dinitrogen = .{
            .to_surface_residue_g_per_step = 1,
            .to_soil_matrix_g_per_step = 2,
            .to_soil_macropore_g_per_step = 3,
        },
        .nitrous_oxide = .{
            .to_surface_residue_g_per_step = 1,
            .to_soil_matrix_g_per_step = 2,
            .to_soil_macropore_g_per_step = 3,
        },
        .ammonium = filledFive(3),
        .ammonia = filledFive(3),
        .nitrate = filledFive(3),
        .hydrogen_phosphate = filledFive(3),
        .dihydrogen_phosphate = filledFive(3),
    }};
    var totals = [_]SoluteFlux{.{}};
    var scratch: [1]SoluteFlux = undefined;
    var state = State{ .net_flux_by_layer = &totals };

    try account(.{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_flux_by_layer = &upper,
        .ground_transfer_by_layer = &ground,
    }, &state, .{ .net_flux_by_layer = &scratch });

    try expectGasValues(totals[0], 14);
    try expectNutrientValues(totals[0], 5);
    try std.testing.expectEqual(@as(f64, 20), totals[0].carbon_dioxide_g_c_per_step + 6);
    try std.testing.expectEqual(@as(f64, 20), totals[0].ammonium_g_n_per_step + 15);
    try std.testing.expectEqual(@as(f64, 20), totals[0].hydrogen_phosphate_g_p_per_step + 15);
}

test "only terminal active layers are boundary-accounted at runtime extent" {
    const capacities = [_]f64{ 2, 2, 0, 2, 2, 2, 2 };
    const upper = [_]SoluteFlux{
        filledSolute(1),
        filledSolute(2),
        filledSolute(3),
        filledSolute(4),
        filledSolute(5),
        filledSolute(6),
        filledSolute(7),
    };
    const ground = [_]GroundTransfer{.{}} ** 7;
    var totals = [_]SoluteFlux{.{}} ** 7;
    var scratch: [7]SoluteFlux = undefined;
    var state = State{ .net_flux_by_layer = &totals };

    try account(.{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_flux_by_layer = &upper,
        .ground_transfer_by_layer = &ground,
    }, &state, .{ .net_flux_by_layer = &scratch });

    try expectGasValues(totals[0], 0);
    try expectGasValues(totals[1], 2);
    try expectGasValues(totals[2], 0);
    try expectGasValues(totals[3], 0);
    try expectGasValues(totals[6], 7);
}

test "source floating point association is retained" {
    const capacities = [_]f64{2};
    const upper = [_]SoluteFlux{filledSolute(-1.0e16)};
    var first_only = filledGround(0, 0);
    first_only.carbon_dioxide.to_surface_residue_g_per_step = 1;
    first_only.methane.to_surface_residue_g_per_step = 1;
    first_only.oxygen.to_surface_residue_g_per_step = 1;
    first_only.dinitrogen.to_surface_residue_g_per_step = 1;
    first_only.nitrous_oxide.to_surface_residue_g_per_step = 1;
    first_only.ammonium.to_surface_residue_nonband_g_per_step = 1;
    first_only.ammonia.to_surface_residue_nonband_g_per_step = 1;
    first_only.nitrate.to_surface_residue_nonband_g_per_step = 1;
    first_only.hydrogen_phosphate.to_surface_residue_nonband_g_per_step = 1;
    first_only.dihydrogen_phosphate.to_surface_residue_nonband_g_per_step = 1;
    const ground = [_]GroundTransfer{first_only};
    var totals = [_]SoluteFlux{filledSolute(1.0e16)};
    var scratch: [1]SoluteFlux = undefined;
    var state = State{ .net_flux_by_layer = &totals };

    try account(.{
        .heat_capacity_mj_per_k_by_layer = &capacities,
        .minimum_heat_capacity_mj_per_k = 1,
        .upper_face_flux_by_layer = &upper,
        .ground_transfer_by_layer = &ground,
    }, &state, .{ .net_flux_by_layer = &scratch });

    try expectGasValues(totals[0], -1);
    try expectNutrientValues(totals[0], -1);
}

test "late invalid input and overflow preserve state atomically" {
    const capacities = [_]f64{ 2, 2 };
    const upper = [_]SoluteFlux{ filledSolute(1), filledSolute(1) };
    var ground = [_]GroundTransfer{ .{}, .{} };
    ground[1].dihydrogen_phosphate.to_soil_macropore_band_g_per_step =
        std.math.nan(f64);
    var totals = [_]SoluteFlux{filledSolute(3)} ** 2;
    var scratch: [2]SoluteFlux = undefined;
    var state = State{ .net_flux_by_layer = &totals };
    const workspace = Workspace{ .net_flux_by_layer = &scratch };

    try std.testing.expectError(
        error.NonFiniteSnowpackBoundarySoluteInput,
        account(.{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_flux_by_layer = &upper,
            .ground_transfer_by_layer = &ground,
        }, &state, workspace),
    );
    for (totals) |value| {
        try expectGasValues(value, 3);
        try expectNutrientValues(value, 3);
    }

    ground[1] = .{};
    var overflowing_upper = upper;
    overflowing_upper[1] = filledSolute(std.math.floatMax(f64));
    totals[1] = filledSolute(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteSnowpackBoundarySoluteResult,
        account(.{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_flux_by_layer = &overflowing_upper,
            .ground_transfer_by_layer = &ground,
        }, &state, workspace),
    );
    try expectGasValues(totals[1], std.math.floatMax(f64));
    try expectNutrientValues(totals[1], std.math.floatMax(f64));
}

test "dimension and workspace alias failures precede mutation" {
    const capacities = [_]f64{ 2, 2 };
    const upper = [_]SoluteFlux{ filledSolute(1), filledSolute(1) };
    const ground = [_]GroundTransfer{ .{}, .{} };
    var totals = [_]SoluteFlux{filledSolute(5)} ** 2;
    var state = State{ .net_flux_by_layer = &totals };

    try std.testing.expectError(
        error.SnowpackBoundarySoluteWorkspaceOverlap,
        account(.{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_flux_by_layer = &upper,
            .ground_transfer_by_layer = &ground,
        }, &state, .{ .net_flux_by_layer = &totals }),
    );
    try expectGasValues(totals[0], 5);

    var short_scratch: [1]SoluteFlux = undefined;
    try std.testing.expectError(
        error.SnowpackBoundarySoluteDimensionMismatch,
        account(.{
            .heat_capacity_mj_per_k_by_layer = &capacities,
            .minimum_heat_capacity_mj_per_k = 1,
            .upper_face_flux_by_layer = &upper,
            .ground_transfer_by_layer = &ground,
        }, &state, .{ .net_flux_by_layer = &short_scratch }),
    );
    try expectNutrientValues(totals[0], 5);
}
