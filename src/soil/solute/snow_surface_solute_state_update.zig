const std = @import("std");

pub const species_count = 41;
pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    salt_simulation: SaltSimulation,
    cell_count: usize,
    /// Cell-major net overland snow-solute flux `[cell][species]` (mol step-1).
    overland_snow_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9358--9424.
/// For every runtime grid cell under ISALTG != 0, the 41-species TQS family
/// is added to snowpack state layer one. H4SiO4 is absent from this family.
pub fn apply(inputs: Inputs, surface_snow_content_mol: []f64) !usize {
    if (inputs.salt_simulation == .disabled) return 0;
    try validate(inputs, surface_snow_content_mol);

    for (0..inputs.cell_count) |cell| {
        const offset = cell * species_count;
        for (0..species_count) |species| {
            const index = offset + species;
            const result = surface_snow_content_mol[index] + inputs.overland_snow_flux_mol_per_step[index];
            if (!std.math.isFinite(result)) return error.NonFiniteSnowSurfaceSoluteResult;
        }
    }
    for (0..inputs.cell_count) |cell| {
        const offset = cell * species_count;
        for (0..species_count) |species| {
            const index = offset + species;
            surface_snow_content_mol[index] = surface_snow_content_mol[index] +
                inputs.overland_snow_flux_mol_per_step[index];
        }
    }
    return inputs.cell_count;
}

fn validate(inputs: Inputs, state: []const f64) !void {
    const expected = std.math.mul(usize, inputs.cell_count, species_count) catch
        return error.SnowSurfaceSoluteDimensionMismatch;
    if (inputs.overland_snow_flux_mol_per_step.len != expected or state.len != expected)
        return error.SnowSurfaceSoluteDimensionMismatch;
    for (inputs.overland_snow_flux_mol_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSnowSurfaceSoluteInput;
    for (state) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSnowSurfaceSoluteInput;
}

test "TRNSFRS updates layer-one snow state for every runtime cell" {
    const cell_count = 3;
    var flux = [_]f64{0} ** (cell_count * species_count);
    for (0..cell_count) |cell| @memset(flux[cell * species_count .. (cell + 1) * species_count], @floatFromInt(cell + 1));
    var state = [_]f64{10} ** (cell_count * species_count);
    try std.testing.expectEqual(cell_count, try apply(.{
        .salt_simulation = .enabled,
        .cell_count = cell_count,
        .overland_snow_flux_mol_per_step = &flux,
    }, &state));
    try std.testing.expectEqual(@as(f64, 11), state[0]);
    try std.testing.expectEqual(@as(f64, 12), state[species_count + 40]);
    try std.testing.expectEqual(@as(f64, 13), state[2 * species_count]);
}

test "disabled salt simulation leaves all snow state unchanged" {
    const flux = [_]f64{};
    var state = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .salt_simulation = .disabled,
        .cell_count = std.math.maxInt(usize),
        .overland_snow_flux_mol_per_step = &flux,
    }, &state));
}

test "zero-cell runtime domain is valid" {
    const empty = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .salt_simulation = .enabled,
        .cell_count = 0,
        .overland_snow_flux_mol_per_step = &empty,
    }, @constCast(&empty)));
}

test "runtime cell topology mismatch fails atomically" {
    const short_flux = [_]f64{2} ** (2 * species_count - 1);
    var state = [_]f64{10} ** (2 * species_count);
    try std.testing.expectError(error.SnowSurfaceSoluteDimensionMismatch, apply(.{
        .salt_simulation = .enabled,
        .cell_count = 2,
        .overland_snow_flux_mol_per_step = &short_flux,
    }, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
}

test "late cell overflow leaves all earlier cells atomic" {
    const cell_count = 2;
    var flux = [_]f64{0} ** (cell_count * species_count);
    flux[flux.len - 1] = std.math.floatMax(f64);
    var state = [_]f64{0} ** (cell_count * species_count);
    state[state.len - 1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSnowSurfaceSoluteResult, apply(.{
        .salt_simulation = .enabled,
        .cell_count = cell_count,
        .overland_snow_flux_mol_per_step = &flux,
    }, &state));
    try std.testing.expectEqual(@as(f64, 0), state[0]);
}
