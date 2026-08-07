const std = @import("std");

pub const species_count = 42;
pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    salt_simulation: SaltSimulation,
    cell_count: usize,
    /// Cell-major TQR net overland fluxes `[cell][species]` (mol step-1).
    overland_flux_mol_per_step: []const f64,
    /// Cell-major R*FLS vertical litter-surface fluxes (mol step-1).
    vertical_surface_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9498--9539.
/// At litter layer zero, all 42 species (including H4SiO4) use exact source
/// evaluation `state + TQR + R*FLS(direction 3, layer 0)`.
pub fn apply(inputs: Inputs, residue_content_mol: []f64) !usize {
    if (inputs.salt_simulation == .disabled) return 0;
    const expected = std.math.mul(usize, inputs.cell_count, species_count) catch
        return error.SurfaceResidueSoluteDimensionMismatch;
    try validate(inputs, residue_content_mol, expected);

    for (0..expected) |index| {
        const after_overland = residue_content_mol[index] + inputs.overland_flux_mol_per_step[index];
        const result = after_overland + inputs.vertical_surface_flux_mol_per_step[index];
        if (!std.math.isFinite(after_overland) or !std.math.isFinite(result))
            return error.NonFiniteSurfaceResidueSoluteResult;
    }
    for (0..expected) |index| {
        residue_content_mol[index] = residue_content_mol[index] + inputs.overland_flux_mol_per_step[index];
        residue_content_mol[index] = residue_content_mol[index] + inputs.vertical_surface_flux_mol_per_step[index];
    }
    return inputs.cell_count;
}

fn validate(inputs: Inputs, state: []const f64, expected: usize) !void {
    if (inputs.overland_flux_mol_per_step.len != expected or
        inputs.vertical_surface_flux_mol_per_step.len != expected or state.len != expected)
        return error.SurfaceResidueSoluteDimensionMismatch;
    const slices = [_][]const f64{
        inputs.overland_flux_mol_per_step,
        inputs.vertical_surface_flux_mol_per_step,
        state,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceResidueSoluteInput;
}

test "TRNSFRS updates all 42 surface residue species for each cell" {
    const cell_count = 2;
    var overland = [_]f64{0} ** (cell_count * species_count);
    @memset(overland[0..species_count], 2);
    @memset(overland[species_count..], 3);
    const vertical = [_]f64{4} ** (cell_count * species_count);
    var state = [_]f64{10} ** (cell_count * species_count);
    try std.testing.expectEqual(cell_count, try apply(.{
        .salt_simulation = .enabled,
        .cell_count = cell_count,
        .overland_flux_mol_per_step = &overland,
        .vertical_surface_flux_mol_per_step = &vertical,
    }, &state));
    try std.testing.expectEqual(@as(f64, 16), state[0]);
    try std.testing.expectEqual(@as(f64, 16), state[41]);
    try std.testing.expectEqual(@as(f64, 17), state[species_count]);
}

test "disabled salt simulation preserves residue state" {
    const overland = [_]f64{};
    const vertical = [_]f64{};
    var state = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .salt_simulation = .disabled,
        .cell_count = std.math.maxInt(usize),
        .overland_flux_mol_per_step = &overland,
        .vertical_surface_flux_mol_per_step = &vertical,
    }, &state));
}

test "zero-cell residue domain is valid" {
    const overland = [_]f64{};
    const vertical = [_]f64{};
    var state = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .salt_simulation = .enabled,
        .cell_count = 0,
        .overland_flux_mol_per_step = &overland,
        .vertical_surface_flux_mol_per_step = &vertical,
    }, &state));
}

test "runtime residue topology mismatch fails atomically" {
    const short_overland = [_]f64{2} ** (species_count - 1);
    const vertical = [_]f64{4} ** species_count;
    var state = [_]f64{10} ** species_count;
    try std.testing.expectError(error.SurfaceResidueSoluteDimensionMismatch, apply(.{
        .salt_simulation = .enabled,
        .cell_count = 1,
        .overland_flux_mol_per_step = &short_overland,
        .vertical_surface_flux_mol_per_step = &vertical,
    }, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
}

test "late three-term overflow leaves earlier cells atomic" {
    var overland = [_]f64{0} ** species_count;
    var vertical = [_]f64{0} ** species_count;
    overland[species_count - 1] = std.math.floatMax(f64);
    vertical[species_count - 1] = std.math.floatMax(f64);
    var state = [_]f64{0} ** species_count;
    try std.testing.expectError(error.NonFiniteSurfaceResidueSoluteResult, apply(.{
        .salt_simulation = .enabled,
        .cell_count = 1,
        .overland_flux_mol_per_step = &overland,
        .vertical_surface_flux_mol_per_step = &vertical,
    }, &state));
    try std.testing.expectEqual(@as(f64, 0), state[0]);
}
