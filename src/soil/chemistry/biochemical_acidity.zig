const std = @import("std");
const chemistry = @import("../solute/chemistry_state.zig");
const aggregation = @import("../biogeochemistry/transformation_aggregation.zig");

pub const Inputs = struct {
    water_volume_m3: f64,
    ammonia_oxidation_g_n: [2]f64,
    nitrate_reduction_g_n: [2]f64,
    nitrite_reduction_g_n: [2]f64,
    nitrous_oxide_reduction_g_n: f64,
    lignin_decomposition_g_c: f64,
    timestep_h: f64,
    negligible_hydrogen_mol: f64,
};

pub const StagedChange = struct {
    initial_hydrogen_mol: f64,
    biochemical_change_mol_h: f64,
    final_hydrogen_mol: f64,
    final_hydrogen_mol_per_m3: f64,
};

/// Units-safe concentration/extensive bridge for NITRO.F XZHYS (4110--4124).
/// STARTE initializes `ZHY = CHYU * FC`; SOLUTE converts it back with
/// `CHY1 = ZHY / VOLW`; REDIST commits the extensive `TRHY` amount.
pub fn stage(state: *const chemistry.State, layer: usize, inputs: Inputs) !StagedChange {
    if (layer >= state.cell_count) return error.ChemistryCellIndexOutOfBounds;
    if (!std.math.isFinite(inputs.water_volume_m3) or inputs.water_volume_m3 < 0)
        return error.InvalidBiochemicalAcidityWaterVolume;
    if (!std.math.isFinite(inputs.negligible_hydrogen_mol) or inputs.negligible_hydrogen_mol < 0)
        return error.InvalidBiochemicalAcidityTolerance;
    const concentration = state.aqueous[layer].hydrogen;
    if (!std.math.isFinite(concentration) or concentration < 0)
        return error.InvalidBiochemicalHydrogenConcentration;
    const initial_mol = concentration * inputs.water_volume_m3;
    if (!std.math.isFinite(initial_mol)) return error.InvalidBiochemicalHydrogenAmount;
    const change_mol = try aggregation.acidityChangeMolH(.{
        .ammonia_oxidation_g_n = inputs.ammonia_oxidation_g_n,
        .nitrate_reduction_g_n = inputs.nitrate_reduction_g_n,
        .nitrite_reduction_g_n = inputs.nitrite_reduction_g_n,
        .nitrous_oxide_reduction_g_n = inputs.nitrous_oxide_reduction_g_n,
        .lignin_decomposition_g_c = inputs.lignin_decomposition_g_c,
        .available_hydrogen_mol_h = initial_mol,
        .timestep_h = inputs.timestep_h,
    });
    const final_mol = initial_mol + change_mol;
    if (!std.math.isFinite(final_mol) or final_mol < -inputs.negligible_hydrogen_mol)
        return error.InvalidBiochemicalHydrogenAmount;
    if (inputs.water_volume_m3 == 0) {
        if (@abs(change_mol) > inputs.negligible_hydrogen_mol)
            return error.BiochemicalAcidityRequiresPositiveWaterVolume;
        return .{ .initial_hydrogen_mol = 0, .biochemical_change_mol_h = 0, .final_hydrogen_mol = 0, .final_hydrogen_mol_per_m3 = concentration };
    }
    return .{
        .initial_hydrogen_mol = initial_mol,
        .biochemical_change_mol_h = change_mol,
        .final_hydrogen_mol = @max(0, final_mol),
        .final_hydrogen_mol_per_m3 = @max(0, final_mol) / inputs.water_volume_m3,
    };
}

pub fn commit(state: *chemistry.State, layer: usize, staged: StagedChange) !void {
    if (layer >= state.cell_count) return error.ChemistryCellIndexOutOfBounds;
    inline for (std.meta.fields(StagedChange)) |field|
        if (!std.math.isFinite(@field(staged, field.name)))
            return error.InvalidBiochemicalHydrogenAmount;
    if (staged.initial_hydrogen_mol < 0 or staged.final_hydrogen_mol < 0 or staged.final_hydrogen_mol_per_m3 < 0)
        return error.InvalidBiochemicalHydrogenAmount;
    state.aqueous[layer].hydrogen = staged.final_hydrogen_mol_per_m3;
}

fn testInputs(water_volume_m3: f64) Inputs {
    return .{
        .water_volume_m3 = water_volume_m3,
        .ammonia_oxidation_g_n = .{ 0, 0 },
        .nitrate_reduction_g_n = .{ 0, 0 },
        .nitrite_reduction_g_n = .{ 0, 0 },
        .nitrous_oxide_reduction_g_n = 0,
        .lignin_decomposition_g_c = 0,
        .timestep_h = 1,
        .negligible_hydrogen_mol = 1e-12,
    };
}

test "wet concentration bridge conserves the staged mol H transaction" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].hydrogen = 2;
    var inputs = testInputs(0.5);
    inputs.ammonia_oxidation_g_n = .{ 2, 1 };
    inputs.nitrate_reduction_g_n = .{ 0.5, 0.5 };
    inputs.nitrite_reduction_g_n = .{ 0.25, 0.25 };
    inputs.nitrous_oxide_reduction_g_n = 0.5;
    inputs.lignin_decomposition_g_c = 2;
    const staged = try stage(&state, 0, inputs);
    try std.testing.expectApproxEqAbs(staged.final_hydrogen_mol, staged.initial_hydrogen_mol + staged.biochemical_change_mol_h, 1e-15);
    try commit(&state, 0, staged);
    try std.testing.expectApproxEqAbs(staged.final_hydrogen_mol, state.aqueous[0].hydrogen * 0.5, 1e-15);
}

test "available hydrogen bound cannot consume more extensive H than exists" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].hydrogen = 0.25;
    var inputs = testInputs(1);
    inputs.nitrate_reduction_g_n = .{ 100, 100 };
    const staged = try stage(&state, 0, inputs);
    try std.testing.expectEqual(@as(f64, -0.25), staged.biochemical_change_mol_h);
    try std.testing.expectEqual(@as(f64, 0), staged.final_hydrogen_mol);
}

test "dry nonzero source fails without mutating chemistry" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].hydrogen = 3;
    var inputs = testInputs(0);
    inputs.ammonia_oxidation_g_n = .{ 1, 0 };
    try std.testing.expectError(error.BiochemicalAcidityRequiresPositiveWaterVolume, stage(&state, 0, inputs));
    try std.testing.expectEqual(@as(f64, 3), state.aqueous[0].hydrogen);
}

test "dry zero source preserves concentration without inventing an amount" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].hydrogen = 3;
    const staged = try stage(&state, 0, testInputs(0));
    try commit(&state, 0, staged);
    try std.testing.expectEqual(@as(f64, 3), state.aqueous[0].hydrogen);
}
