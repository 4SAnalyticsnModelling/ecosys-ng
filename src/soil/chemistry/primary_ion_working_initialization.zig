const std = @import("std");
const source_initialization = @import("primary_ion_source_initialization.zig");

pub const Dimensions = source_initialization.Dimensions;
pub const LoopPosition = source_initialization.LoopPosition;

pub const Inputs = struct {
    ammonium_extract_mol_n_per_megagram_by_cell_layer: []const f64, // CN4X
    nitrate_extract_mol_n_per_megagram_by_cell_layer: []const f64, // CNOX
    phosphate_extract_mol_p_per_megagram_by_cell_layer: []const f64, // CPOX
    primary_sources_by_cell_layer: []const source_initialization.State,
};

/// STARTE working concentrations for the active soil layer. These retain the
/// soil-file mass basis (`mol Mg-1`), unlike precipitation and irrigation.
pub const State = struct {
    ammonium_mol_n_per_megagram: f64 = 0, // CN4Z
    nitrate_mol_n_per_megagram: f64 = 0, // CNOZ
    phosphate_mol_p_per_megagram: f64 = 0, // CPOZ
    aluminum_mol_per_megagram: f64 = 0, // CALZ
    iron_mol_per_megagram: f64 = 0, // CFEZ
    calcium_mol_per_megagram: f64 = 0, // CCAZ
    magnesium_mol_per_megagram: f64 = 0, // CMGZ
    sodium_mol_per_megagram: f64 = 0, // CNAZ
    potassium_mol_per_megagram: f64 = 0, // CKAZ
    sulfate_mol_s_per_megagram: f64 = 0, // CSOZ
    chloride_mol_per_megagram: f64 = 0, // CCLZ
};

/// Direct translation of `starte.f` lines 212--222 within the `K=3,I=1`
/// branch. All finite signed values pass through because later STARTE logic
/// interprets signed source values before publishing physical concentrations.
pub fn initialize(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    inputs: Inputs,
    state: *State,
) !bool {
    if (position.source_index != 3 or position.day_index != 1) return false;

    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.layer_count_including_surface == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count or
        position.layer_index >= dimensions.layer_count_including_surface)
        return error.InvalidSoilWorkingIonDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidSoilWorkingIonDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.layer_count_including_surface) catch
        return error.InvalidSoilWorkingIonDimensions;
    if (inputs.ammonium_extract_mol_n_per_megagram_by_cell_layer.len != value_count or
        inputs.nitrate_extract_mol_n_per_megagram_by_cell_layer.len != value_count or
        inputs.phosphate_extract_mol_p_per_megagram_by_cell_layer.len != value_count or
        inputs.primary_sources_by_cell_layer.len != value_count)
        return error.InvalidSoilWorkingIonDimensions;

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.layer_count_including_surface + position.layer_index;
    const source = inputs.primary_sources_by_cell_layer[value_index];
    const next: State = .{
        .ammonium_mol_n_per_megagram = inputs.ammonium_extract_mol_n_per_megagram_by_cell_layer[value_index],
        .nitrate_mol_n_per_megagram = inputs.nitrate_extract_mol_n_per_megagram_by_cell_layer[value_index],
        .phosphate_mol_p_per_megagram = inputs.phosphate_extract_mol_p_per_megagram_by_cell_layer[value_index],
        .aluminum_mol_per_megagram = source.aluminum_mol_per_megagram,
        .iron_mol_per_megagram = source.iron_mol_per_megagram,
        .calcium_mol_per_megagram = source.calcium_mol_per_megagram,
        .magnesium_mol_per_megagram = source.magnesium_mol_per_megagram,
        .sodium_mol_per_megagram = source.sodium_mol_per_megagram,
        .potassium_mol_per_megagram = source.potassium_mol_per_megagram,
        .sulfate_mol_s_per_megagram = source.sulfate_mol_s_per_megagram,
        .chloride_mol_per_megagram = source.chloride_mol_per_megagram,
    };
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteSoilWorkingIonInput;
    }
    state.* = next;
    return true;
}

fn fixtureSource(value: f64) source_initialization.State {
    return .{
        .aluminum_mol_per_megagram = value,
        .iron_mol_per_megagram = value + 1,
        .calcium_mol_per_megagram = value + 2,
        .magnesium_mol_per_megagram = value + 3,
        .sodium_mol_per_megagram = value + 4,
        .potassium_mol_per_megagram = value + 5,
        .sulfate_mol_s_per_megagram = value + 6,
        .chloride_mol_per_megagram = value + 7,
    };
}

test "STARTE soil working ions preserve source order and runtime layer topology" {
    const ammonium = [_]f64{ 1, 2, 3, 4 };
    const nitrate = [_]f64{ 5, 6, 7, 8 };
    const phosphate = [_]f64{ 9, 10, 11, 12 };
    const sources = [_]source_initialization.State{ fixtureSource(20), fixtureSource(30), fixtureSource(40), fixtureSource(-50) };
    var state: State = .{};
    try std.testing.expect(try initialize(.{ .source_index = 3, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .layer_count_including_surface = 2 }, 0, 1, .{
        .ammonium_extract_mol_n_per_megagram_by_cell_layer = &ammonium,
        .nitrate_extract_mol_n_per_megagram_by_cell_layer = &nitrate,
        .phosphate_extract_mol_p_per_megagram_by_cell_layer = &phosphate,
        .primary_sources_by_cell_layer = &sources,
    }, &state));
    try std.testing.expectEqual(@as(f64, 4), state.ammonium_mol_n_per_megagram);
    try std.testing.expectEqual(@as(f64, 8), state.nitrate_mol_n_per_megagram);
    try std.testing.expectEqual(@as(f64, 12), state.phosphate_mol_p_per_megagram);
    try std.testing.expectEqual(@as(f64, -50), state.aluminum_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -43), state.chloride_mol_per_megagram);
}

test "STARTE inactive soil working guard keeps invalid dormant data untouched" {
    const empty_f64 = [_]f64{};
    const empty_sources = [_]source_initialization.State{};
    var state: State = .{ .ammonium_mol_n_per_megagram = 3 };
    const before = state;
    try std.testing.expect(!try initialize(.{ .source_index = 3, .day_index = 2, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .layer_count_including_surface = 0 }, 99, 99, .{
        .ammonium_extract_mol_n_per_megagram_by_cell_layer = &empty_f64,
        .nitrate_extract_mol_n_per_megagram_by_cell_layer = &empty_f64,
        .phosphate_extract_mol_p_per_megagram_by_cell_layer = &empty_f64,
        .primary_sources_by_cell_layer = &empty_sources,
    }, &state));
    try std.testing.expectEqualDeep(before, state);
}

test "STARTE soil working snapshot rejects late non-finite source atomically" {
    const valid = [_]f64{1};
    var bad_source = fixtureSource(2);
    bad_source.chloride_mol_per_megagram = std.math.nan(f64);
    const sources = [_]source_initialization.State{bad_source};
    var state: State = .{ .ammonium_mol_n_per_megagram = 8, .chloride_mol_per_megagram = 9 };
    const before = state;
    try std.testing.expectError(error.NonFiniteSoilWorkingIonInput, initialize(.{ .source_index = 3, .day_index = 1, .layer_index = 0 }, .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 }, 0, 0, .{
        .ammonium_extract_mol_n_per_megagram_by_cell_layer = &valid,
        .nitrate_extract_mol_n_per_megagram_by_cell_layer = &valid,
        .phosphate_extract_mol_p_per_megagram_by_cell_layer = &valid,
        .primary_sources_by_cell_layer = &sources,
    }, &state));
    try std.testing.expectEqualDeep(before, state);
}
