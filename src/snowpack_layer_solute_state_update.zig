const std = @import("std");

pub const species_count = 41;
pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    salt_simulation: SaltSimulation,
    cell_count: usize,
    snow_layer_count: usize,
    /// Cell-major, then layer-major net snowpack flux `[cell][layer][species]`.
    net_snowpack_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9425--9467.
/// Under the enclosing cell loops and ISALTG guard, every Fortran layer
/// `L=1..JS` receives its 41 ordered T*BLS flux totals.
pub fn apply(inputs: Inputs, snowpack_content_mol: []f64) !usize {
    if (inputs.salt_simulation == .disabled) return 0;
    const entry_count = try expectedEntryCount(inputs.cell_count, inputs.snow_layer_count);
    try validate(inputs, snowpack_content_mol, entry_count);

    for (0..entry_count) |index| {
        const result = snowpack_content_mol[index] + inputs.net_snowpack_flux_mol_per_step[index];
        if (!std.math.isFinite(result)) return error.NonFiniteSnowpackLayerSoluteResult;
    }
    for (0..entry_count) |index|
        snowpack_content_mol[index] = snowpack_content_mol[index] +
            inputs.net_snowpack_flux_mol_per_step[index];
    return inputs.cell_count * inputs.snow_layer_count;
}

fn expectedEntryCount(cell_count: usize, layer_count: usize) !usize {
    const cell_layers = std.math.mul(usize, cell_count, layer_count) catch
        return error.SnowpackLayerSoluteDimensionMismatch;
    return std.math.mul(usize, cell_layers, species_count) catch
        return error.SnowpackLayerSoluteDimensionMismatch;
}

fn validate(inputs: Inputs, state: []const f64, expected: usize) !void {
    if (inputs.net_snowpack_flux_mol_per_step.len != expected or state.len != expected)
        return error.SnowpackLayerSoluteDimensionMismatch;
    for (inputs.net_snowpack_flux_mol_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSnowpackLayerSoluteInput;
    for (state) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSnowpackLayerSoluteInput;
}

test "TRNSFRS updates every runtime cell and snow layer" {
    const cell_count = 2;
    const layer_count = 3;
    var flux = [_]f64{0} ** (cell_count * layer_count * species_count);
    for (0..cell_count) |cell| {
        for (0..layer_count) |layer| {
            const offset = (cell * layer_count + layer) * species_count;
            @memset(flux[offset .. offset + species_count], @floatFromInt(cell * layer_count + layer + 1));
        }
    }
    var state = [_]f64{10} ** (cell_count * layer_count * species_count);
    try std.testing.expectEqual(cell_count * layer_count, try apply(.{
        .salt_simulation = .enabled,
        .cell_count = cell_count,
        .snow_layer_count = layer_count,
        .net_snowpack_flux_mol_per_step = &flux,
    }, &state));
    try std.testing.expectEqual(@as(f64, 11), state[0]);
    try std.testing.expectEqual(@as(f64, 13), state[2 * species_count + 40]);
    try std.testing.expectEqual(@as(f64, 16), state[5 * species_count]);
}

test "disabled salt simulation preserves all layers" {
    const flux = [_]f64{};
    var state = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .salt_simulation = .disabled,
        .cell_count = std.math.maxInt(usize),
        .snow_layer_count = std.math.maxInt(usize),
        .net_snowpack_flux_mol_per_step = &flux,
    }, &state));
}

test "runtime layer count can be zero" {
    const flux = [_]f64{};
    var state = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .salt_simulation = .enabled,
        .cell_count = 4,
        .snow_layer_count = 0,
        .net_snowpack_flux_mol_per_step = &flux,
    }, &state));
}

test "runtime snow layer topology mismatch fails atomically" {
    const short_flux = [_]f64{2} ** (2 * species_count - 1);
    var state = [_]f64{10} ** (2 * species_count);
    try std.testing.expectError(error.SnowpackLayerSoluteDimensionMismatch, apply(.{
        .salt_simulation = .enabled,
        .cell_count = 1,
        .snow_layer_count = 2,
        .net_snowpack_flux_mol_per_step = &short_flux,
    }, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
}

test "late snow layer overflow leaves earlier layers atomic" {
    var flux = [_]f64{0} ** (2 * species_count);
    flux[flux.len - 1] = std.math.floatMax(f64);
    var state = [_]f64{0} ** (2 * species_count);
    state[state.len - 1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSnowpackLayerSoluteResult, apply(.{
        .salt_simulation = .enabled,
        .cell_count = 1,
        .snow_layer_count = 2,
        .net_snowpack_flux_mol_per_step = &flux,
    }, &state));
    try std.testing.expectEqual(@as(f64, 0), state[0]);
}
