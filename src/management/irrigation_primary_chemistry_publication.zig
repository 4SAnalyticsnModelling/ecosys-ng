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

pub const GasConcentrations = struct {
    carbon_dioxide_mol_per_m3: f64, // CCOQ
    methane_mol_per_m3: f64, // CCHQ
    oxygen_mol_per_m3: f64, // COXQ
    dinitrogen_mol_n_per_m3: f64, // CNNQ
    nitrous_oxide_mol_n_per_m3: f64, // CN2Q
};

pub const DissolvedConcentrations = struct {
    ammonium_mol_n_per_m3: f64, // CN4Q
    ammonia_mol_n_per_m3: f64, // CN3Q
    aluminum_mol_per_m3: f64, // CALQ
    iron_mol_per_m3: f64, // CFEQ
    hydrogen_mol_per_m3: f64, // CHYQ
    calcium_mol_per_m3: f64, // CCAQ
    magnesium_mol_per_m3: f64, // CMGQ
    sodium_mol_per_m3: f64, // CNAQ
    potassium_mol_per_m3: f64, // CKAQ
    hydroxide_mol_per_m3: f64, // COHQ
    sulfate_mol_per_m3: f64, // CSOQ
    chloride_mol_per_m3: f64, // CCLQ
    carbonate_mol_per_m3: f64, // CC3Q
    bicarbonate_mol_per_m3: f64, // CHCQ
};

pub const WorkingConcentrations = struct {
    gases: GasConcentrations,
    dissolved: DissolvedConcentrations,
};

/// Direct translation of `starte.f` lines 1291--1310. Gas values are cell-only
/// in the source, while dissolved irrigation chemistry is stored with Fortran
/// day-fastest `(I,NY,NX)` topology.
pub fn publish(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    working: WorkingConcentrations,
    gases_by_cell: []GasConcentrations,
    dissolved_by_cell_day: []DissolvedConcentrations,
) !bool {
    if (position.source_index != 2 or position.layer_index != 1) return false;
    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.day_count == 0 or row >= dimensions.row_count or
        column >= dimensions.column_count or position.day_index == 0 or
        position.day_index > dimensions.day_count)
        return error.InvalidIrrigationChemistryDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidIrrigationChemistryDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.day_count) catch
        return error.InvalidIrrigationChemistryDimensions;
    if (gases_by_cell.len != cell_count or dissolved_by_cell_day.len != value_count)
        return error.InvalidIrrigationChemistryDimensions;
    try validate(GasConcentrations, working.gases);
    try validate(DissolvedConcentrations, working.dissolved);

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.day_count + position.day_index - 1;
    gases_by_cell[cell] = working.gases;
    dissolved_by_cell_day[value_index] = working.dissolved;
    return true;
}

fn validate(comptime T: type, values: T) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        const value = @field(values, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidIrrigationWorkingConcentration;
    }
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE irrigation publication preserves cell gases and day-fastest solutes" {
    var gases = [_]GasConcentrations{filled(GasConcentrations, 9)} ** 2;
    var dissolved = [_]DissolvedConcentrations{filled(DissolvedConcentrations, 8)} ** 4;
    var working: WorkingConcentrations = .{
        .gases = filled(GasConcentrations, 1),
        .dissolved = filled(DissolvedConcentrations, 2),
    };
    working.gases.nitrous_oxide_mol_n_per_m3 = 3;
    working.dissolved.bicarbonate_mol_per_m3 = 4;
    try std.testing.expect(try publish(.{ .source_index = 2, .day_index = 2, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .day_count = 2 }, 0, 1, working, &gases, &dissolved));
    try std.testing.expectEqual(@as(f64, 9), gases[0].carbon_dioxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 3), gases[1].nitrous_oxide_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 8), dissolved[2].ammonium_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 4), dissolved[3].bicarbonate_mol_per_m3);
}

test "STARTE inactive irrigation publication leaves dormant data untouched" {
    var gases = [_]GasConcentrations{filled(GasConcentrations, 7)};
    var dissolved = [_]DissolvedConcentrations{filled(DissolvedConcentrations, 6)};
    const gases_before = gases;
    const dissolved_before = dissolved;
    const invalid: WorkingConcentrations = .{ .gases = filled(GasConcentrations, std.math.nan(f64)), .dissolved = filled(DissolvedConcentrations, std.math.nan(f64)) };
    try std.testing.expect(!try publish(.{ .source_index = 1, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .day_count = 0 }, 99, 99, invalid, &gases, &dissolved));
    try std.testing.expectEqualDeep(gases_before, gases);
    try std.testing.expectEqualDeep(dissolved_before, dissolved);
}

test "STARTE irrigation publication rejects late invalid value atomically" {
    var gases = [_]GasConcentrations{filled(GasConcentrations, 7)};
    var dissolved = [_]DissolvedConcentrations{filled(DissolvedConcentrations, 6)};
    const gases_before = gases;
    const dissolved_before = dissolved;
    var working: WorkingConcentrations = .{ .gases = filled(GasConcentrations, 1), .dissolved = filled(DissolvedConcentrations, 2) };
    working.dissolved.bicarbonate_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(error.InvalidIrrigationWorkingConcentration, publish(.{ .source_index = 2, .day_index = 1, .layer_index = 1 }, .{ .column_count = 1, .row_count = 1, .day_count = 1 }, 0, 0, working, &gases, &dissolved));
    try std.testing.expectEqualDeep(gases_before, gases);
    try std.testing.expectEqualDeep(dissolved_before, dissolved);
}
