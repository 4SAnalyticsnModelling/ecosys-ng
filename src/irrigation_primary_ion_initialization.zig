const std = @import("std");
const precipitation = @import("precipitation_primary_ion_initialization.zig");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    day_count: usize,
};

/// One-based STARTE loop coordinates (`K`, `I`, and `L`).
pub const LoopPosition = precipitation.LoopPosition;
pub const Constants = precipitation.Constants;
pub const State = precipitation.State;

/// Irrigation chemistry uses Fortran day-fastest `(I,NY,NX)` storage.
pub const Inputs = struct {
    ph_by_cell_day: []const f64, // PHQ
    ammonium_g_n_per_m3_by_cell_day: []const f64, // CN4Q
    nitrate_g_n_per_m3_by_cell_day: []const f64, // CNOQ
    phosphate_g_p_per_m3_by_cell_day: []const f64, // CPOQ
    aluminum_mol_per_m3_by_cell_day: []const f64, // CALQ
    iron_mol_per_m3_by_cell_day: []const f64, // CFEQ
    calcium_mol_per_m3_by_cell_day: []const f64, // CCAQ
    magnesium_mol_per_m3_by_cell_day: []const f64, // CMGQ
    sodium_mol_per_m3_by_cell_day: []const f64, // CNAQ
    potassium_mol_per_m3_by_cell_day: []const f64, // CKAQ
    sulfate_mol_per_m3_by_cell_day: []const f64, // CSOQ
    chloride_mol_per_m3_by_cell_day: []const f64, // CCLQ
};

/// Direct translation of STARTE.F lines 149--163 for irrigation (`K=2`).
/// The source guard runs before any array access, so inactive loop positions
/// deliberately leave invalid dormant irrigation inputs unexamined.
pub fn initialize(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    inputs: Inputs,
    constants: Constants,
    state: *State,
) !bool {
    if (position.source_index != 2 or position.layer_index != 1) return false;

    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.day_count == 0 or row >= dimensions.row_count or
        column >= dimensions.column_count or position.day_index == 0 or
        position.day_index > dimensions.day_count)
        return error.InvalidIrrigationGridDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidIrrigationGridDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.day_count) catch
        return error.InvalidIrrigationGridDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (@field(inputs, field.name).len != value_count)
            return error.InvalidIrrigationInputDimensions;
    }
    if (!std.math.isFinite(constants.water_dissociation_mol2_per_m6) or
        constants.water_dissociation_mol2_per_m6 <= 0 or
        !std.math.isFinite(constants.nitrogen_g_per_mol_n) or
        constants.nitrogen_g_per_mol_n <= 0 or
        !std.math.isFinite(constants.phosphorus_g_per_mol_p) or
        constants.phosphorus_g_per_mol_p <= 0)
        return error.InvalidIrrigationIonConstant;

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.day_count + position.day_index - 1;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name)[value_index];
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidIrrigationIonInput;
    }
    const ph = inputs.ph_by_cell_day[value_index];
    if (ph > 14) return error.InvalidIrrigationIonInput;

    const hydrogen = std.math.pow(f64, 10.0, -(ph - 3.0));
    const next: State = .{
        .ph = ph,
        .hydrogen_mol_per_m3 = hydrogen,
        .hydroxide_mol_per_m3 = constants.water_dissociation_mol2_per_m6 / hydrogen,
        .ammonium_mol_n_per_m3 = inputs.ammonium_g_n_per_m3_by_cell_day[value_index] / constants.nitrogen_g_per_mol_n,
        .nitrate_mol_n_per_m3 = inputs.nitrate_g_n_per_m3_by_cell_day[value_index] / constants.nitrogen_g_per_mol_n,
        .phosphate_mol_p_per_m3 = inputs.phosphate_g_p_per_m3_by_cell_day[value_index] / constants.phosphorus_g_per_mol_p,
        .aluminum_mol_per_m3 = inputs.aluminum_mol_per_m3_by_cell_day[value_index],
        .iron_mol_per_m3 = inputs.iron_mol_per_m3_by_cell_day[value_index],
        .calcium_mol_per_m3 = inputs.calcium_mol_per_m3_by_cell_day[value_index],
        .magnesium_mol_per_m3 = inputs.magnesium_mol_per_m3_by_cell_day[value_index],
        .sodium_mol_per_m3 = inputs.sodium_mol_per_m3_by_cell_day[value_index],
        .potassium_mol_per_m3 = inputs.potassium_mol_per_m3_by_cell_day[value_index],
        .sulfate_mol_per_m3 = inputs.sulfate_mol_per_m3_by_cell_day[value_index],
        .chloride_mol_per_m3 = inputs.chloride_mol_per_m3_by_cell_day[value_index],
    };
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteIrrigationIonState;
    }
    state.* = next;
    return true;
}

fn uniformInputs(values: []const f64) Inputs {
    return .{
        .ph_by_cell_day = values,
        .ammonium_g_n_per_m3_by_cell_day = values,
        .nitrate_g_n_per_m3_by_cell_day = values,
        .phosphate_g_p_per_m3_by_cell_day = values,
        .aluminum_mol_per_m3_by_cell_day = values,
        .iron_mol_per_m3_by_cell_day = values,
        .calcium_mol_per_m3_by_cell_day = values,
        .magnesium_mol_per_m3_by_cell_day = values,
        .sodium_mol_per_m3_by_cell_day = values,
        .potassium_mol_per_m3_by_cell_day = values,
        .sulfate_mol_per_m3_by_cell_day = values,
        .chloride_mol_per_m3_by_cell_day = values,
    };
}

const source_constants: Constants = .{
    .water_dissociation_mol2_per_m6 = 1.0e-8,
    .nitrogen_g_per_mol_n = 14,
    .phosphorus_g_per_mol_p = 31,
};

test "STARTE irrigation ions use day-fastest cell storage and source molar bases" {
    const ph = [_]f64{ 4, 5, 6, 7 };
    const ammonium = [_]f64{ 0, 0, 0, 28 };
    const nitrate = [_]f64{ 0, 0, 0, 42 };
    const phosphate = [_]f64{ 0, 0, 0, 62 };
    var inputs = uniformInputs(&ph);
    inputs.ammonium_g_n_per_m3_by_cell_day = &ammonium;
    inputs.nitrate_g_n_per_m3_by_cell_day = &nitrate;
    inputs.phosphate_g_p_per_m3_by_cell_day = &phosphate;
    var state: State = .{};
    try std.testing.expect(try initialize(.{ .source_index = 2, .day_index = 2, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .day_count = 2 }, 0, 1, inputs, source_constants, &state));
    try std.testing.expectEqual(@as(f64, 7), state.ph);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydrogen_mol_per_m3, 1.0e-16);
    try std.testing.expectEqual(@as(f64, 2), state.ammonium_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 3), state.nitrate_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 2), state.phosphate_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 7), state.chloride_mol_per_m3);
}

test "STARTE inactive irrigation guard preserves state and dormant invalid data" {
    var state: State = .{ .ph = 8, .iron_mol_per_m3 = 4 };
    const before = state;
    const empty = [_]f64{};
    try std.testing.expect(!try initialize(.{ .source_index = 1, .day_index = 0, .layer_index = 1 }, .{ .column_count = 0, .row_count = 0, .day_count = 0 }, 99, 99, uniformInputs(&empty), .{ .water_dissociation_mol2_per_m6 = std.math.nan(f64), .nitrogen_g_per_mol_n = 0, .phosphorus_g_per_mol_p = 0 }, &state));
    try std.testing.expectEqualDeep(before, state);
}

test "STARTE irrigation initialization rejects late invalid input atomically" {
    const values = [_]f64{1};
    const invalid = [_]f64{std.math.nan(f64)};
    var inputs = uniformInputs(&values);
    inputs.chloride_mol_per_m3_by_cell_day = &invalid;
    var state: State = .{ .ph = 5, .chloride_mol_per_m3 = 6 };
    const before = state;
    try std.testing.expectError(error.InvalidIrrigationIonInput, initialize(.{ .source_index = 2, .day_index = 1, .layer_index = 1 }, .{ .column_count = 1, .row_count = 1, .day_count = 1 }, 0, 0, inputs, source_constants, &state));
    try std.testing.expectEqualDeep(before, state);
}
