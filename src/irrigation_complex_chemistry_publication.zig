const std = @import("std");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    day_count: usize,
};

pub const LoopPosition = struct {
    source_index: usize,
    day_index: usize,
    layer_index: usize,
};

/// Field order directly follows STARTE.F 1311--1341. All concentrations are
/// `mol m-3`, or `mol P m-3` for phosphate-bearing species.
pub const CellDayChemistry = struct {
    aluminum_monohydroxide: f64, // CAL1Q
    aluminum_dihydroxide: f64, // CAL2Q
    aluminum_trihydroxide: f64, // CAL3Q
    aluminum_tetrahydroxide: f64, // CAL4Q
    aluminum_sulfate: f64, // CALSQ
    iron_monohydroxide: f64, // CFE1Q
    iron_dihydroxide: f64, // CFE2Q
    iron_trihydroxide: f64, // CFE3Q
    iron_tetrahydroxide: f64, // CFE4Q
    iron_sulfate: f64, // CFESQ
    calcium_hydroxide: f64, // CCAOQ
    calcium_carbonate: f64, // CCACQ
    calcium_bicarbonate: f64, // CCAHQ
    calcium_sulfate: f64, // CCASQ
    magnesium_hydroxide: f64, // CMGOQ
    magnesium_carbonate: f64, // CMGCQ
    magnesium_bicarbonate: f64, // CMGHQ
    magnesium_sulfate: f64, // CMGSQ
    sodium_carbonate: f64, // CNACQ
    sodium_sulfate: f64, // CNASQ
    potassium_sulfate: f64, // CKASQ
    phosphate: f64, // CH0PQ
    hydrogen_phosphate: f64, // CH1PQ
    dihydrogen_phosphate: f64, // CPOQ
    phosphoric_acid: f64, // CH3PQ
    iron_hydrogen_phosphate: f64, // CF1PQ
    iron_dihydrogen_phosphate: f64, // CF2PQ
    calcium_phosphate: f64, // CC0PQ
    calcium_hydrogen_phosphate: f64, // CC1PQ
    calcium_dihydrogen_phosphate: f64, // CC2PQ
    magnesium_hydrogen_phosphate: f64, // CM1PQ
};

/// Direct translation of STARTE.F lines 1311--1341. Irrigation chemistry is
/// stored with Fortran day-fastest `(I,NY,NX)` topology.
pub fn publish(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    working: CellDayChemistry,
    chemistry_by_cell_day: []CellDayChemistry,
) !bool {
    if (position.source_index != 2 or position.layer_index != 1) return false;
    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.day_count == 0 or row >= dimensions.row_count or
        column >= dimensions.column_count or position.day_index == 0 or
        position.day_index > dimensions.day_count)
        return error.InvalidIrrigationComplexDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidIrrigationComplexDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.day_count) catch
        return error.InvalidIrrigationComplexDimensions;
    if (chemistry_by_cell_day.len != value_count)
        return error.InvalidIrrigationComplexDimensions;
    inline for (@typeInfo(CellDayChemistry).@"struct".fields) |field| {
        const value = @field(working, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidIrrigationComplexConcentration;
    }

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.day_count + position.day_index - 1;
    chemistry_by_cell_day[value_index] = working;
    return true;
}

fn filled(value: f64) CellDayChemistry {
    var result: CellDayChemistry = undefined;
    inline for (@typeInfo(CellDayChemistry).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "STARTE irrigation complexes preserve source fields and day-fastest topology" {
    var working = filled(1);
    working.aluminum_monohydroxide = 2;
    working.potassium_sulfate = 3;
    working.dihydrogen_phosphate = 4;
    working.magnesium_hydrogen_phosphate = 5;
    var chemistry = [_]CellDayChemistry{filled(9)} ** 4;
    try std.testing.expect(try publish(.{ .source_index = 2, .day_index = 2, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .day_count = 2 }, 0, 1, working, &chemistry));
    try std.testing.expectEqual(@as(f64, 9), chemistry[2].aluminum_monohydroxide);
    try std.testing.expectEqual(@as(f64, 2), chemistry[3].aluminum_monohydroxide);
    try std.testing.expectEqual(@as(f64, 3), chemistry[3].potassium_sulfate);
    try std.testing.expectEqual(@as(f64, 4), chemistry[3].dihydrogen_phosphate);
    try std.testing.expectEqual(@as(f64, 5), chemistry[3].magnesium_hydrogen_phosphate);
}

test "STARTE inactive irrigation complex guard leaves invalid dormant data untouched" {
    var chemistry = [_]CellDayChemistry{filled(7)};
    const before = chemistry;
    try std.testing.expect(!try publish(.{ .source_index = 1, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .day_count = 0 }, 99, 99, filled(std.math.nan(f64)), &chemistry));
    try std.testing.expectEqualDeep(before, chemistry);
}

test "STARTE irrigation complex publication rejects late invalid field atomically" {
    var working = filled(1);
    working.magnesium_hydrogen_phosphate = std.math.nan(f64);
    var chemistry = [_]CellDayChemistry{filled(6)};
    const before = chemistry;
    try std.testing.expectError(error.InvalidIrrigationComplexConcentration, publish(.{ .source_index = 2, .day_index = 1, .layer_index = 1 }, .{ .column_count = 1, .row_count = 1, .day_count = 1 }, 0, 0, working, &chemistry));
    try std.testing.expectEqualDeep(before, chemistry);
}
