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

pub const GasConcentrations = struct {
    carbon_dioxide_mol_per_m3: f64, // CCOU
    methane_mol_per_m3: f64, // CCHU
    oxygen_mol_per_m3: f64, // COXU; source forces zero
    dinitrogen_mol_n_per_m3: f64, // CNNU
    nitrous_oxide_mol_n_per_m3: f64, // CN2U
};

pub const DissolvedConcentrations = struct {
    ammonium_mol_n_per_source_volume: f64, // CN4U
    ammonia_mol_n_per_source_volume: f64, // CN3U
    nitrate_extract_mol_n_per_megagram: f64, // CNOU receives CNOX
    aluminum_mol_per_source_volume: f64, // CALU
    iron_mol_per_source_volume: f64, // CFEU
    hydrogen_mol_per_m3: f64, // CHYU
    calcium_mol_per_source_volume: f64, // CCAU
    magnesium_mol_per_source_volume: f64, // CMGU
    sodium_mol_per_source_volume: f64, // CNAU
    potassium_mol_per_source_volume: f64, // CKAU
    hydroxide_mol_per_m3: f64, // COHU
    sulfate_mol_s_per_source_volume: f64, // CSOU
    chloride_mol_per_source_volume: f64, // CCLU
    carbonate_mol_per_m3: f64, // CC3U
    bicarbonate_mol_per_m3: f64, // CHCU
};

pub const WorkingConcentrations = struct {
    carbon_dioxide_mol_per_m3: f64,
    methane_mol_per_m3: f64,
    dinitrogen_mol_n_per_m3: f64,
    nitrous_oxide_mol_n_per_m3: f64,
    dissolved: DissolvedConcentrations,
};

pub const Publication = struct {
    applied: bool,
    gases: ?GasConcentrations,
};

/// Direct translation of STARTE.F lines 1348--1368. Dissolved storage is
/// cell-major then runtime layer, including surface layer zero.
pub fn publish(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    working: WorkingConcentrations,
    dissolved_by_cell_layer: []DissolvedConcentrations,
) !Publication {
    if (position.source_index != 3 or position.day_index != 1)
        return .{ .applied = false, .gases = null };
    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.layer_count_including_surface == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count or
        position.layer_index >= dimensions.layer_count_including_surface)
        return error.InvalidSoilChemistryPublicationDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidSoilChemistryPublicationDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.layer_count_including_surface) catch
        return error.InvalidSoilChemistryPublicationDimensions;
    if (dissolved_by_cell_layer.len != value_count)
        return error.InvalidSoilChemistryPublicationDimensions;
    inline for (.{
        working.carbon_dioxide_mol_per_m3,
        working.methane_mol_per_m3,
        working.dinitrogen_mol_n_per_m3,
        working.nitrous_oxide_mol_n_per_m3,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSoilWorkingConcentration;
    inline for (@typeInfo(DissolvedConcentrations).@"struct".fields) |field| {
        const value = @field(working.dissolved, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilWorkingConcentration;
    }

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.layer_count_including_surface + position.layer_index;
    dissolved_by_cell_layer[value_index] = working.dissolved;
    return .{
        .applied = true,
        .gases = .{
            .carbon_dioxide_mol_per_m3 = working.carbon_dioxide_mol_per_m3,
            .methane_mol_per_m3 = working.methane_mol_per_m3,
            .oxygen_mol_per_m3 = 0.0,
            .dinitrogen_mol_n_per_m3 = working.dinitrogen_mol_n_per_m3,
            .nitrous_oxide_mol_n_per_m3 = working.nitrous_oxide_mol_n_per_m3,
        },
    };
}

fn filledDissolved(value: f64) DissolvedConcentrations {
    var result: DissolvedConcentrations = undefined;
    inline for (@typeInfo(DissolvedConcentrations).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "STARTE soil primary publication preserves layer topology and forced zero oxygen" {
    var layers = [_]DissolvedConcentrations{filledDissolved(9)} ** 4;
    var dissolved = filledDissolved(2);
    dissolved.nitrate_extract_mol_n_per_megagram = 7;
    dissolved.bicarbonate_mol_per_m3 = 8;
    const result = try publish(.{ .source_index = 3, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .layer_count_including_surface = 2 }, 0, 1, .{
        .carbon_dioxide_mol_per_m3 = 1,
        .methane_mol_per_m3 = 2,
        .dinitrogen_mol_n_per_m3 = 3,
        .nitrous_oxide_mol_n_per_m3 = 4,
        .dissolved = dissolved,
    }, &layers);
    try std.testing.expect(result.applied);
    try std.testing.expectEqual(@as(f64, 0), result.gases.?.oxygen_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 9), layers[2].nitrate_extract_mol_n_per_megagram);
    try std.testing.expectEqual(@as(f64, 7), layers[3].nitrate_extract_mol_n_per_megagram);
    try std.testing.expectEqual(@as(f64, 8), layers[3].bicarbonate_mol_per_m3);
}

test "STARTE inactive soil publication leaves invalid dormant state untouched" {
    var layers = [_]DissolvedConcentrations{filledDissolved(6)};
    const before = layers;
    const result = try publish(.{ .source_index = 2, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .layer_count_including_surface = 0 }, 99, 99, .{
        .carbon_dioxide_mol_per_m3 = std.math.nan(f64),
        .methane_mol_per_m3 = std.math.nan(f64),
        .dinitrogen_mol_n_per_m3 = std.math.nan(f64),
        .nitrous_oxide_mol_n_per_m3 = std.math.nan(f64),
        .dissolved = filledDissolved(std.math.nan(f64)),
    }, &layers);
    try std.testing.expect(!result.applied and result.gases == null);
    try std.testing.expectEqualDeep(before, layers);
}

test "STARTE soil publication rejects late invalid concentration atomically" {
    var layers = [_]DissolvedConcentrations{filledDissolved(6)};
    const before = layers;
    var dissolved = filledDissolved(1);
    dissolved.bicarbonate_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilWorkingConcentration, publish(.{ .source_index = 3, .day_index = 1, .layer_index = 0 }, .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 }, 0, 0, .{
        .carbon_dioxide_mol_per_m3 = 1,
        .methane_mol_per_m3 = 1,
        .dinitrogen_mol_n_per_m3 = 1,
        .nitrous_oxide_mol_n_per_m3 = 1,
        .dissolved = dissolved,
    }, &layers));
    try std.testing.expectEqualDeep(before, layers);
}
