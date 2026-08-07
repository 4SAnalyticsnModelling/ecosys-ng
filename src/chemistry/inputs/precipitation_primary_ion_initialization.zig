const std = @import("std");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
};

/// One-based STARTE loop coordinates (`K`, `I`, and `L`).
pub const LoopPosition = struct {
    source_index: usize,
    day_index: usize,
    layer_index: usize,
};

pub const Inputs = struct {
    ph_by_cell: []const f64, // PHR
    ammonium_g_n_per_m3_by_cell: []const f64, // CN4R
    nitrate_g_n_per_m3_by_cell: []const f64, // CNOR
    phosphate_g_p_per_m3_by_cell: []const f64, // CPOR
    aluminum_mol_per_m3_by_cell: []const f64, // CALR
    iron_mol_per_m3_by_cell: []const f64, // CFER
    calcium_mol_per_m3_by_cell: []const f64, // CCAR
    magnesium_mol_per_m3_by_cell: []const f64, // CMGR
    sodium_mol_per_m3_by_cell: []const f64, // CNAR
    potassium_mol_per_m3_by_cell: []const f64, // CKAR
    sulfate_mol_per_m3_by_cell: []const f64, // CSOR
    chloride_mol_per_m3_by_cell: []const f64, // CCLR
};

pub const Constants = struct {
    water_dissociation_mol2_per_m6: f64, // DPH2O
    nitrogen_g_per_mol_n: f64,
    phosphorus_g_per_mol_p: f64,
};

pub const State = struct {
    ph: f64 = 0, // PH1
    hydrogen_mol_per_m3: f64 = 0, // CHY1
    hydroxide_mol_per_m3: f64 = 0, // COH1
    ammonium_mol_n_per_m3: f64 = 0, // CN4Z
    nitrate_mol_n_per_m3: f64 = 0, // CNOZ
    phosphate_mol_p_per_m3: f64 = 0, // CPOZ
    aluminum_mol_per_m3: f64 = 0, // CALZ
    iron_mol_per_m3: f64 = 0, // CFEZ
    calcium_mol_per_m3: f64 = 0, // CCAZ
    magnesium_mol_per_m3: f64 = 0, // CMGZ
    sodium_mol_per_m3: f64 = 0, // CNAZ
    potassium_mol_per_m3: f64 = 0, // CKAZ
    sulfate_mol_per_m3: f64 = 0, // CSOZ
    chloride_mol_per_m3: f64 = 0, // CCLZ
};

/// Direct translation of `starte.f` lines 127--141 for precipitation (`K=1`).
/// Inactive loop positions return before inspecting dimensions or input data,
/// preserving the source guard's dormant-input behavior.
pub fn initialize(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    inputs: Inputs,
    constants: Constants,
    state: *State,
) !bool {
    if (position.source_index != 1 or position.day_index != 1 or position.layer_index != 1)
        return false;

    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidPrecipitationGridDimensions;
    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count)
        return error.InvalidPrecipitationGridDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (@field(inputs, field.name).len != cell_count)
            return error.InvalidPrecipitationInputDimensions;
    }
    if (!std.math.isFinite(constants.water_dissociation_mol2_per_m6) or
        constants.water_dissociation_mol2_per_m6 <= 0 or
        !std.math.isFinite(constants.nitrogen_g_per_mol_n) or
        constants.nitrogen_g_per_mol_n <= 0 or
        !std.math.isFinite(constants.phosphorus_g_per_mol_p) or
        constants.phosphorus_g_per_mol_p <= 0)
        return error.InvalidPrecipitationIonConstant;

    const cell = row * dimensions.column_count + column;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name)[cell];
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPrecipitationIonInput;
    }
    const ph = inputs.ph_by_cell[cell];
    if (ph > 14) return error.InvalidPrecipitationIonInput;

    const hydrogen = std.math.pow(f64, 10.0, -(ph - 3.0));
    const next: State = .{
        .ph = ph,
        .hydrogen_mol_per_m3 = hydrogen,
        .hydroxide_mol_per_m3 = constants.water_dissociation_mol2_per_m6 / hydrogen,
        .ammonium_mol_n_per_m3 = inputs.ammonium_g_n_per_m3_by_cell[cell] / constants.nitrogen_g_per_mol_n,
        .nitrate_mol_n_per_m3 = inputs.nitrate_g_n_per_m3_by_cell[cell] / constants.nitrogen_g_per_mol_n,
        .phosphate_mol_p_per_m3 = inputs.phosphate_g_p_per_m3_by_cell[cell] / constants.phosphorus_g_per_mol_p,
        .aluminum_mol_per_m3 = inputs.aluminum_mol_per_m3_by_cell[cell],
        .iron_mol_per_m3 = inputs.iron_mol_per_m3_by_cell[cell],
        .calcium_mol_per_m3 = inputs.calcium_mol_per_m3_by_cell[cell],
        .magnesium_mol_per_m3 = inputs.magnesium_mol_per_m3_by_cell[cell],
        .sodium_mol_per_m3 = inputs.sodium_mol_per_m3_by_cell[cell],
        .potassium_mol_per_m3 = inputs.potassium_mol_per_m3_by_cell[cell],
        .sulfate_mol_per_m3 = inputs.sulfate_mol_per_m3_by_cell[cell],
        .chloride_mol_per_m3 = inputs.chloride_mol_per_m3_by_cell[cell],
    };
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFinitePrecipitationIonState;
    }
    state.* = next;
    return true;
}

fn uniformInputs(values: []const f64) Inputs {
    return .{
        .ph_by_cell = values,
        .ammonium_g_n_per_m3_by_cell = values,
        .nitrate_g_n_per_m3_by_cell = values,
        .phosphate_g_p_per_m3_by_cell = values,
        .aluminum_mol_per_m3_by_cell = values,
        .iron_mol_per_m3_by_cell = values,
        .calcium_mol_per_m3_by_cell = values,
        .magnesium_mol_per_m3_by_cell = values,
        .sodium_mol_per_m3_by_cell = values,
        .potassium_mol_per_m3_by_cell = values,
        .sulfate_mol_per_m3_by_cell = values,
        .chloride_mol_per_m3_by_cell = values,
    };
}

test "STARTE precipitation primary ions preserve source order and molar bases" {
    const ph = [_]f64{ 6, 7 };
    const n4 = [_]f64{ 0, 28 };
    const no3 = [_]f64{ 0, 42 };
    const p = [_]f64{ 0, 62 };
    var inputs = uniformInputs(&ph);
    inputs.ammonium_g_n_per_m3_by_cell = &n4;
    inputs.nitrate_g_n_per_m3_by_cell = &no3;
    inputs.phosphate_g_p_per_m3_by_cell = &p;
    var state: State = .{};
    try std.testing.expect(try initialize(.{ .source_index = 1, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1 }, 0, 1, inputs, .{ .water_dissociation_mol2_per_m6 = 1.0e-8, .nitrogen_g_per_mol_n = 14, .phosphorus_g_per_mol_p = 31 }, &state));
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydrogen_mol_per_m3, 1.0e-16);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydroxide_mol_per_m3, 1.0e-16);
    try std.testing.expectEqual(@as(f64, 2), state.ammonium_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 3), state.nitrate_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 2), state.phosphate_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 7), state.chloride_mol_per_m3);
}

test "STARTE inactive precipitation guard leaves state and invalid inputs dormant" {
    var state: State = .{ .ph = 8, .calcium_mol_per_m3 = 9 };
    const before = state;
    const empty = [_]f64{};
    try std.testing.expect(!try initialize(.{ .source_index = 2, .day_index = 1, .layer_index = 1 }, .{ .column_count = 0, .row_count = 0 }, 99, 99, uniformInputs(&empty), .{ .water_dissociation_mol2_per_m6 = std.math.nan(f64), .nitrogen_g_per_mol_n = 0, .phosphorus_g_per_mol_p = 0 }, &state));
    try std.testing.expectEqualDeep(before, state);
}

test "STARTE precipitation initialization fails atomically on late invalid input" {
    const values = [_]f64{1};
    const invalid = [_]f64{std.math.nan(f64)};
    var inputs = uniformInputs(&values);
    inputs.chloride_mol_per_m3_by_cell = &invalid;
    var state: State = .{ .ph = 5, .chloride_mol_per_m3 = 6 };
    const before = state;
    try std.testing.expectError(error.InvalidPrecipitationIonInput, initialize(.{ .source_index = 1, .day_index = 1, .layer_index = 1 }, .{ .column_count = 1, .row_count = 1 }, 0, 0, inputs, .{ .water_dissociation_mol2_per_m6 = 1.0e-8, .nitrogen_g_per_mol_n = 14, .phosphorus_g_per_mol_p = 31 }, &state));
    try std.testing.expectEqualDeep(before, state);
}
