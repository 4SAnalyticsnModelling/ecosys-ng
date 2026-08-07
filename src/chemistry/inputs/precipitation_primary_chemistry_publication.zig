const std = @import("std");
const precipitation_initialization = @import("precipitation_primary_ion_initialization.zig");

pub const Dimensions = precipitation_initialization.Dimensions;
pub const LoopPosition = precipitation_initialization.LoopPosition;

/// Equilibrated working concentrations for the current precipitation source.
/// Every field is `mol m-3` (N or P elemental basis where applicable).
pub const WorkingConcentrations = struct {
    carbon_dioxide: f64, // CCO21
    methane: f64, // CCH41
    oxygen: f64, // COXY1
    dinitrogen_n: f64, // CZ2G1
    nitrous_oxide_n: f64, // CZ2O1
    ammonium_n: f64, // CN41
    ammonia_n: f64, // CN31
    aluminum: f64, // CAL1
    iron: f64, // CFE1
    hydrogen: f64, // CHY1
    calcium: f64, // CCA1
    magnesium: f64, // CMG1
    sodium: f64, // CNA1
    potassium: f64, // CKA1
    hydroxide: f64, // COH1
    sulfate: f64, // CSO41
    chloride: f64, // CCL1
    carbonate: f64, // CCO31
    bicarbonate: f64, // CHCO31
};

pub const CellChemistry = WorkingConcentrations;

/// Direct translation of `starte.f` lines 1234--1253. Runtime precipitation
/// chemistry uses one record per logical grid cell; iteration/source guards
/// run before dimensions or dormant working concentrations are inspected.
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
        return error.InvalidPrecipitationChemistryDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidPrecipitationChemistryDimensions;
    if (chemistry_by_cell.len != cell_count)
        return error.InvalidPrecipitationChemistryDimensions;
    inline for (@typeInfo(WorkingConcentrations).@"struct".fields) |field| {
        const value = @field(working, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPrecipitationWorkingConcentration;
    }

    const cell = row * dimensions.column_count + column;
    chemistry_by_cell[cell] = .{
        .carbon_dioxide = working.carbon_dioxide,
        .methane = working.methane,
        .oxygen = working.oxygen,
        .dinitrogen_n = working.dinitrogen_n,
        .nitrous_oxide_n = working.nitrous_oxide_n,
        .ammonium_n = working.ammonium_n,
        .ammonia_n = working.ammonia_n,
        .aluminum = working.aluminum,
        .iron = working.iron,
        .hydrogen = working.hydrogen,
        .calcium = working.calcium,
        .magnesium = working.magnesium,
        .sodium = working.sodium,
        .potassium = working.potassium,
        .hydroxide = working.hydroxide,
        .sulfate = working.sulfate,
        .chloride = working.chloride,
        .carbonate = working.carbonate,
        .bicarbonate = working.bicarbonate,
    };
    return true;
}

fn filled(value: f64) WorkingConcentrations {
    var result: WorkingConcentrations = undefined;
    inline for (@typeInfo(WorkingConcentrations).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "STARTE precipitation primary publication preserves cell topology and order" {
    var working = filled(1);
    working.carbon_dioxide = 2;
    working.nitrous_oxide_n = 3;
    working.hydrogen = 4;
    working.bicarbonate = 5;
    var cells = [_]CellChemistry{filled(9)} ** 2;
    try std.testing.expect(try publish(.{ .source_index = 1, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1 }, 0, 1, working, &cells));
    try std.testing.expectEqual(@as(f64, 9), cells[0].carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 2), cells[1].carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 3), cells[1].nitrous_oxide_n);
    try std.testing.expectEqual(@as(f64, 4), cells[1].hydrogen);
    try std.testing.expectEqual(@as(f64, 5), cells[1].bicarbonate);
}

test "STARTE inactive precipitation publication leaves invalid dormant data untouched" {
    var cells = [_]CellChemistry{filled(7)};
    const before = cells;
    try std.testing.expect(!try publish(.{ .source_index = 2, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0 }, 99, 99, filled(std.math.nan(f64)), &cells));
    try std.testing.expectEqualDeep(before, cells);
}

test "STARTE precipitation publication rejects late invalid field atomically" {
    var working = filled(1);
    working.bicarbonate = std.math.nan(f64);
    var cells = [_]CellChemistry{filled(6)};
    const before = cells;
    try std.testing.expectError(error.InvalidPrecipitationWorkingConcentration, publish(.{ .source_index = 1, .day_index = 1, .layer_index = 1 }, .{ .column_count = 1, .row_count = 1 }, 0, 0, working, &cells));
    try std.testing.expectEqualDeep(before, cells);
}
