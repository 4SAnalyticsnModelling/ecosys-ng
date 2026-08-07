const std = @import("std");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
};

pub const LoopPosition = struct {
    source_index: usize,
    day_index: usize,
    layer_index: usize,
};

pub const TrivalentMetalComplexes = struct {
    monohydroxide_mol_per_m3: f64,
    dihydroxide_mol_per_m3: f64,
    trihydroxide_mol_per_m3: f64,
    tetrahydroxide_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
};

pub const BaseCationComplexes = struct {
    calcium_hydroxide_mol_per_m3: f64,
    calcium_carbonate_mol_per_m3: f64,
    calcium_bicarbonate_mol_per_m3: f64,
    calcium_sulfate_mol_per_m3: f64,
    magnesium_hydroxide_mol_per_m3: f64,
    magnesium_carbonate_mol_per_m3: f64,
    magnesium_bicarbonate_mol_per_m3: f64,
    magnesium_sulfate_mol_per_m3: f64,
    sodium_carbonate_mol_per_m3: f64,
    sodium_sulfate_mol_per_m3: f64,
    potassium_sulfate_mol_per_m3: f64,
};

pub const PhosphateSpecies = struct {
    phosphate_mol_p_per_m3: f64,
    hydrogen_phosphate_mol_p_per_m3: f64,
    dihydrogen_phosphate_mol_p_per_m3: f64,
    phosphoric_acid_mol_p_per_m3: f64,
    iron_hydrogen_phosphate_mol_p_per_m3: f64,
    iron_dihydrogen_phosphate_mol_p_per_m3: f64,
    calcium_phosphate_mol_p_per_m3: f64,
    calcium_hydrogen_phosphate_mol_p_per_m3: f64,
    calcium_dihydrogen_phosphate_mol_p_per_m3: f64,
    magnesium_hydrogen_phosphate_mol_p_per_m3: f64,
};

pub const CellChemistry = struct {
    aluminum: TrivalentMetalComplexes,
    iron: TrivalentMetalComplexes,
    base_cations: BaseCationComplexes,
    phosphate: PhosphateSpecies,
};

pub const WorkingConcentrations = CellChemistry;

/// Direct translation of `starte.f` lines 1254--1284. All values are aqueous
/// concentrations in `mol m-3` (`mol P m-3` for phosphate-bearing species).
/// The guard precedes topology and concentration validation exactly as in the
/// enclosing precipitation publication branch.
pub fn publish(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    working: WorkingConcentrations,
    chemistry_by_cell: []CellChemistry,
) !bool {
    if (position.source_index != 1 or position.day_index != 1 or position.layer_index != 1)
        return false;

    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count)
        return error.InvalidPrecipitationComplexDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidPrecipitationComplexDimensions;
    if (chemistry_by_cell.len != cell_count)
        return error.InvalidPrecipitationComplexDimensions;
    try validateGroup(TrivalentMetalComplexes, working.aluminum);
    try validateGroup(TrivalentMetalComplexes, working.iron);
    try validateGroup(BaseCationComplexes, working.base_cations);
    try validateGroup(PhosphateSpecies, working.phosphate);

    const cell = row * dimensions.column_count + column;
    chemistry_by_cell[cell] = .{
        .aluminum = .{
            .monohydroxide_mol_per_m3 = working.aluminum.monohydroxide_mol_per_m3,
            .dihydroxide_mol_per_m3 = working.aluminum.dihydroxide_mol_per_m3,
            .trihydroxide_mol_per_m3 = working.aluminum.trihydroxide_mol_per_m3,
            .tetrahydroxide_mol_per_m3 = working.aluminum.tetrahydroxide_mol_per_m3,
            .sulfate_mol_per_m3 = working.aluminum.sulfate_mol_per_m3,
        },
        .iron = .{
            .monohydroxide_mol_per_m3 = working.iron.monohydroxide_mol_per_m3,
            .dihydroxide_mol_per_m3 = working.iron.dihydroxide_mol_per_m3,
            .trihydroxide_mol_per_m3 = working.iron.trihydroxide_mol_per_m3,
            .tetrahydroxide_mol_per_m3 = working.iron.tetrahydroxide_mol_per_m3,
            .sulfate_mol_per_m3 = working.iron.sulfate_mol_per_m3,
        },
        .base_cations = .{
            .calcium_hydroxide_mol_per_m3 = working.base_cations.calcium_hydroxide_mol_per_m3,
            .calcium_carbonate_mol_per_m3 = working.base_cations.calcium_carbonate_mol_per_m3,
            .calcium_bicarbonate_mol_per_m3 = working.base_cations.calcium_bicarbonate_mol_per_m3,
            .calcium_sulfate_mol_per_m3 = working.base_cations.calcium_sulfate_mol_per_m3,
            .magnesium_hydroxide_mol_per_m3 = working.base_cations.magnesium_hydroxide_mol_per_m3,
            .magnesium_carbonate_mol_per_m3 = working.base_cations.magnesium_carbonate_mol_per_m3,
            .magnesium_bicarbonate_mol_per_m3 = working.base_cations.magnesium_bicarbonate_mol_per_m3,
            .magnesium_sulfate_mol_per_m3 = working.base_cations.magnesium_sulfate_mol_per_m3,
            .sodium_carbonate_mol_per_m3 = working.base_cations.sodium_carbonate_mol_per_m3,
            .sodium_sulfate_mol_per_m3 = working.base_cations.sodium_sulfate_mol_per_m3,
            .potassium_sulfate_mol_per_m3 = working.base_cations.potassium_sulfate_mol_per_m3,
        },
        .phosphate = .{
            .phosphate_mol_p_per_m3 = working.phosphate.phosphate_mol_p_per_m3,
            .hydrogen_phosphate_mol_p_per_m3 = working.phosphate.hydrogen_phosphate_mol_p_per_m3,
            .dihydrogen_phosphate_mol_p_per_m3 = working.phosphate.dihydrogen_phosphate_mol_p_per_m3,
            .phosphoric_acid_mol_p_per_m3 = working.phosphate.phosphoric_acid_mol_p_per_m3,
            .iron_hydrogen_phosphate_mol_p_per_m3 = working.phosphate.iron_hydrogen_phosphate_mol_p_per_m3,
            .iron_dihydrogen_phosphate_mol_p_per_m3 = working.phosphate.iron_dihydrogen_phosphate_mol_p_per_m3,
            .calcium_phosphate_mol_p_per_m3 = working.phosphate.calcium_phosphate_mol_p_per_m3,
            .calcium_hydrogen_phosphate_mol_p_per_m3 = working.phosphate.calcium_hydrogen_phosphate_mol_p_per_m3,
            .calcium_dihydrogen_phosphate_mol_p_per_m3 = working.phosphate.calcium_dihydrogen_phosphate_mol_p_per_m3,
            .magnesium_hydrogen_phosphate_mol_p_per_m3 = working.phosphate.magnesium_hydrogen_phosphate_mol_p_per_m3,
        },
    };
    return true;
}

fn validateGroup(comptime T: type, group: T) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        const value = @field(group, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPrecipitationComplexConcentration;
    }
}

fn filled(value: f64) CellChemistry {
    return .{
        .aluminum = filledGroup(TrivalentMetalComplexes, value),
        .iron = filledGroup(TrivalentMetalComplexes, value),
        .base_cations = filledGroup(BaseCationComplexes, value),
        .phosphate = filledGroup(PhosphateSpecies, value),
    };
}

fn filledGroup(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE precipitation complex publication preserves runtime cell and source order" {
    var working = filled(1);
    working.aluminum.monohydroxide_mol_per_m3 = 2;
    working.iron.sulfate_mol_per_m3 = 3;
    working.base_cations.potassium_sulfate_mol_per_m3 = 4;
    working.phosphate.dihydrogen_phosphate_mol_p_per_m3 = 5;
    working.phosphate.magnesium_hydrogen_phosphate_mol_p_per_m3 = 6;
    var cells = [_]CellChemistry{filled(9)} ** 2;
    try std.testing.expect(try publish(.{ .source_index = 1, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1 }, 0, 1, working, &cells));
    try std.testing.expectEqual(@as(f64, 9), cells[0].aluminum.monohydroxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 2), cells[1].aluminum.monohydroxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 3), cells[1].iron.sulfate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), cells[1].base_cations.potassium_sulfate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 5), cells[1].phosphate.dihydrogen_phosphate_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 6), cells[1].phosphate.magnesium_hydrogen_phosphate_mol_p_per_m3);
}

test "STARTE inactive precipitation complex guard leaves dormant data untouched" {
    var cells = [_]CellChemistry{filled(7)};
    const before = cells;
    try std.testing.expect(!try publish(.{ .source_index = 1, .day_index = 2, .layer_index = 1 }, .{ .column_count = 0, .row_count = 0 }, 99, 99, filled(std.math.nan(f64)), &cells));
    try std.testing.expectEqualDeep(before, cells);
}

test "STARTE precipitation complex publication rejects late invalid value atomically" {
    var working = filled(1);
    working.phosphate.magnesium_hydrogen_phosphate_mol_p_per_m3 = std.math.nan(f64);
    var cells = [_]CellChemistry{filled(8)};
    const before = cells;
    try std.testing.expectError(error.InvalidPrecipitationComplexConcentration, publish(.{ .source_index = 1, .day_index = 1, .layer_index = 1 }, .{ .column_count = 1, .row_count = 1 }, 0, 0, working, &cells));
    try std.testing.expectEqualDeep(before, cells);
}
