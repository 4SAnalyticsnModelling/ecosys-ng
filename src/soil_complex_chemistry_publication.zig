const std = @import("std");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    layer_count_including_surface: usize,
};

pub const LoopPosition = struct {
    source_index: usize,
    day_index: usize,
    layer_index: usize,
};

/// Field order directly follows STARTE.F 1369--1399. Values are `mol m-3`,
/// or `mol P m-3` for phosphate-bearing species.
pub const ComplexConcentrations = struct {
    aluminum_monohydroxide_mol_per_m3: f64, // CAL1U receives CALO1
    aluminum_dihydroxide_mol_per_m3: f64, // CAL2U receives CALO2
    aluminum_trihydroxide_mol_per_m3: f64, // CAL3U receives CALO3
    aluminum_tetrahydroxide_mol_per_m3: f64, // CAL4U receives CALO4
    aluminum_sulfate_mol_per_m3: f64, // CALSU receives CALS1
    iron_monohydroxide_mol_per_m3: f64, // CFE1U receives CFEO1
    iron_dihydroxide_mol_per_m3: f64, // CFE2U receives CFEO2
    iron_trihydroxide_mol_per_m3: f64, // CFE3U receives CFEO3
    iron_tetrahydroxide_mol_per_m3: f64, // CFE4U receives CFEO4
    iron_sulfate_mol_per_m3: f64, // CFESU receives CFES1
    calcium_hydroxide_mol_per_m3: f64, // CCAOU receives CCAO1
    calcium_carbonate_mol_per_m3: f64, // CCACU receives CCAC1
    calcium_bicarbonate_mol_per_m3: f64, // CCAHU receives CCAH1
    calcium_sulfate_mol_per_m3: f64, // CCASU receives CCAS1
    magnesium_hydroxide_mol_per_m3: f64, // CMGOU receives CMGO1
    magnesium_carbonate_mol_per_m3: f64, // CMGCU receives CMGC1
    magnesium_bicarbonate_mol_per_m3: f64, // CMGHU receives CMGH1
    magnesium_sulfate_mol_per_m3: f64, // CMGSU receives CMGS1
    sodium_carbonate_mol_per_m3: f64, // CNACU receives CNAC1
    sodium_sulfate_mol_per_m3: f64, // CNASU receives CNAS1
    potassium_sulfate_mol_per_m3: f64, // CKASU receives CKAS1
    phosphate_mol_p_per_m3: f64, // CH0PU receives CH0P1
    hydrogen_phosphate_mol_p_per_m3: f64, // CH1PU receives CH1P1
    dihydrogen_phosphate_mol_p_per_m3: f64, // CH2PU receives CH2P1
    phosphoric_acid_mol_p_per_m3: f64, // CH3PU receives CH3P1
    iron_hydrogen_phosphate_mol_p_per_m3: f64, // CF1PU receives CF1P1
    iron_dihydrogen_phosphate_mol_p_per_m3: f64, // CF2PU receives CF2P1
    calcium_phosphate_mol_p_per_m3: f64, // CC0PU receives CC0P1
    calcium_hydrogen_phosphate_mol_p_per_m3: f64, // CC1PU receives CC1P1
    calcium_dihydrogen_phosphate_mol_p_per_m3: f64, // CC2PU receives CC2P1
    magnesium_hydrogen_phosphate_mol_p_per_m3: f64, // CM1PU receives CM1P1
};

/// Direct translation of STARTE.F lines 1369--1399. Storage is cell-major
/// then runtime layer, including surface layer zero.
pub fn publish(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    working: ComplexConcentrations,
    chemistry_by_cell_layer: []ComplexConcentrations,
) !bool {
    if (position.source_index != 3 or position.day_index != 1) return false;
    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.layer_count_including_surface == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count or
        position.layer_index >= dimensions.layer_count_including_surface)
        return error.InvalidSoilComplexDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidSoilComplexDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.layer_count_including_surface) catch
        return error.InvalidSoilComplexDimensions;
    if (chemistry_by_cell_layer.len != value_count)
        return error.InvalidSoilComplexDimensions;
    inline for (@typeInfo(ComplexConcentrations).@"struct".fields) |field| {
        const value = @field(working, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilComplexConcentration;
    }

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.layer_count_including_surface + position.layer_index;
    chemistry_by_cell_layer[value_index] = working;
    return true;
}

fn filled(value: f64) ComplexConcentrations {
    var result: ComplexConcentrations = undefined;
    inline for (@typeInfo(ComplexConcentrations).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "STARTE soil complex publication preserves source fields and layer topology" {
    var working = filled(1);
    working.aluminum_monohydroxide_mol_per_m3 = 2;
    working.potassium_sulfate_mol_per_m3 = 3;
    working.dihydrogen_phosphate_mol_p_per_m3 = 4;
    working.magnesium_hydrogen_phosphate_mol_p_per_m3 = 5;
    var chemistry = [_]ComplexConcentrations{filled(9)} ** 4;
    try std.testing.expect(try publish(.{ .source_index = 3, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .layer_count_including_surface = 2 }, 0, 1, working, &chemistry));
    try std.testing.expectEqual(@as(f64, 9), chemistry[2].aluminum_monohydroxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 2), chemistry[3].aluminum_monohydroxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 3), chemistry[3].potassium_sulfate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), chemistry[3].dihydrogen_phosphate_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 5), chemistry[3].magnesium_hydrogen_phosphate_mol_p_per_m3);
}

test "STARTE inactive soil complex guard leaves invalid dormant data untouched" {
    var chemistry = [_]ComplexConcentrations{filled(7)};
    const before = chemistry;
    try std.testing.expect(!try publish(.{ .source_index = 2, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .layer_count_including_surface = 0 }, 99, 99, filled(std.math.nan(f64)), &chemistry));
    try std.testing.expectEqualDeep(before, chemistry);
}

test "STARTE soil complex publication rejects late invalid field atomically" {
    var working = filled(1);
    working.magnesium_hydrogen_phosphate_mol_p_per_m3 = std.math.nan(f64);
    var chemistry = [_]ComplexConcentrations{filled(6)};
    const before = chemistry;
    try std.testing.expectError(error.InvalidSoilComplexConcentration, publish(.{ .source_index = 3, .day_index = 1, .layer_index = 0 }, .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 }, 0, 0, working, &chemistry));
    try std.testing.expectEqualDeep(before, chemistry);
}
