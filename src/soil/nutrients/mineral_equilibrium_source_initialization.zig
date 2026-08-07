const std = @import("std");
const source_initialization = @import("../chemistry/primary_ion_source_initialization.zig");

pub const Dimensions = source_initialization.Dimensions;
pub const LoopPosition = source_initialization.LoopPosition;

pub const State = struct {
    aluminum_hydroxide_mol_per_megagram: f64 = 0, // PALOH1
    iron_hydroxide_mol_per_megagram: f64 = 0, // PFEOH1
    calcium_carbonate_mol_per_megagram: f64 = 0, // PCACO1
    calcium_sulfate_mol_per_megagram: f64 = 0, // PCASO1
    aluminum_phosphate_mol_per_megagram: f64 = 0, // PALPO1
    iron_phosphate_mol_per_megagram: f64 = 0, // PFEPO1
    calcium_hydrogen_phosphate_mol_per_megagram: f64 = 0, // PCAPD1
    apatite_mol_per_megagram: f64 = 0, // PCAPH1
};

/// Direct translation of `starte.f` lines 375--383. Source mineral values use
/// the soil-file `mol Mg-1` basis and are selected for the active runtime
/// cell-layer only when the enclosing source is soil (`K=3`).
pub fn initialize(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    sources_by_cell_layer: []const source_initialization.State,
    state: *State,
) !bool {
    if (position.source_index != 3) return false;

    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.layer_count_including_surface == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count or
        position.layer_index >= dimensions.layer_count_including_surface)
        return error.InvalidSoilMineralEquilibriumDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidSoilMineralEquilibriumDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.layer_count_including_surface) catch
        return error.InvalidSoilMineralEquilibriumDimensions;
    if (sources_by_cell_layer.len != value_count)
        return error.InvalidSoilMineralEquilibriumDimensions;

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.layer_count_including_surface + position.layer_index;
    const source = sources_by_cell_layer[value_index];
    const next: State = .{
        .aluminum_hydroxide_mol_per_megagram = source.aluminum_hydroxide_mol_per_megagram,
        .iron_hydroxide_mol_per_megagram = source.iron_hydroxide_mol_per_megagram,
        .calcium_carbonate_mol_per_megagram = source.calcium_carbonate_mol_per_megagram,
        .calcium_sulfate_mol_per_megagram = source.calcium_sulfate_mol_per_megagram,
        .aluminum_phosphate_mol_per_megagram = source.aluminum_phosphate_mol_per_megagram,
        .iron_phosphate_mol_per_megagram = source.iron_phosphate_mol_per_megagram,
        .calcium_hydrogen_phosphate_mol_per_megagram = source.calcium_hydrogen_phosphate_mol_per_megagram,
        .apatite_mol_per_megagram = source.apatite_mol_per_megagram,
    };
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteSoilMineralEquilibriumSource;
    }
    state.* = next;
    return true;
}

fn fixtureSource(start: f64) source_initialization.State {
    return .{
        .aluminum_phosphate_mol_per_megagram = start,
        .iron_phosphate_mol_per_megagram = start + 1,
        .calcium_hydrogen_phosphate_mol_per_megagram = start + 2,
        .apatite_mol_per_megagram = start + 3,
        .aluminum_hydroxide_mol_per_megagram = start + 4,
        .iron_hydroxide_mol_per_megagram = start + 5,
        .calcium_carbonate_mol_per_megagram = start + 6,
        .calcium_sulfate_mol_per_megagram = start + 7,
    };
}

test "STARTE soil mineral equilibrium snapshot preserves runtime topology and order" {
    const sources = [_]source_initialization.State{ fixtureSource(0), fixtureSource(10), fixtureSource(20), fixtureSource(-30) };
    var state: State = .{};
    try std.testing.expect(try initialize(.{ .source_index = 3, .day_index = 5, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .layer_count_including_surface = 2 }, 0, 1, &sources, &state));
    try std.testing.expectEqual(@as(f64, -26), state.aluminum_hydroxide_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -25), state.iron_hydroxide_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -24), state.calcium_carbonate_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -23), state.calcium_sulfate_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -30), state.aluminum_phosphate_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -29), state.iron_phosphate_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -28), state.calcium_hydrogen_phosphate_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, -27), state.apatite_mol_per_megagram);
}

test "STARTE non-soil mineral branch leaves dormant invalid topology untouched" {
    var state: State = .{ .calcium_carbonate_mol_per_megagram = 7 };
    const before = state;
    const empty = [_]source_initialization.State{};
    try std.testing.expect(!try initialize(.{ .source_index = 2, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .layer_count_including_surface = 0 }, 99, 99, &empty, &state));
    try std.testing.expectEqualDeep(before, state);
}

test "STARTE soil mineral snapshot rejects late non-finite value atomically" {
    var invalid = fixtureSource(1);
    invalid.apatite_mol_per_megagram = std.math.nan(f64);
    const sources = [_]source_initialization.State{invalid};
    var state: State = .{ .aluminum_hydroxide_mol_per_megagram = 8, .apatite_mol_per_megagram = 9 };
    const before = state;
    try std.testing.expectError(error.NonFiniteSoilMineralEquilibriumSource, initialize(.{ .source_index = 3, .day_index = 1, .layer_index = 0 }, .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 }, 0, 0, &sources, &state));
    try std.testing.expectEqualDeep(before, state);
}
