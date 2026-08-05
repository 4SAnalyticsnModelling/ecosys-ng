const std = @import("std");

/// Litter discharge has 34 salt/complex fields (including H4SiO4) followed by
/// eight non-band phosphate fields.
pub const litter_species_count = 42;
/// Soil discharge additionally carries eight band phosphate fields.
pub const soil_species_count = 50;

pub const SurfaceFluxWorkspace = struct {
    /// Named replacement for legacy integer sentinel `ICHKL`.
    surface_discharge_claimed: bool,
    /// Snowmelt solute flux routed to surface litter, mol step-1.
    litter_flux_mol_per_step: []f64,
    /// Snowmelt solute flux routed to non-band and band soil, mol step-1.
    soil_flux_mol_per_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 1945--2037.
///
/// Dimensions for both differently sized source groups are validated before
/// the sentinel or either flux workspace is mutated.
pub fn reset(workspace: *SurfaceFluxWorkspace) !void {
    if (workspace.litter_flux_mol_per_step.len != litter_species_count)
        return error.SnowmeltLitterSoluteFluxDimensionMismatch;
    if (workspace.soil_flux_mol_per_step.len != soil_species_count)
        return error.SnowmeltSoilSoluteFluxDimensionMismatch;

    workspace.surface_discharge_claimed = false;
    @memset(workspace.litter_flux_mol_per_step, 0);
    @memset(workspace.soil_flux_mol_per_step, 0);
}

test "TRNSFRS resets the discharge sentinel and both exact source groups" {
    var litter = [_]f64{3} ** litter_species_count;
    var soil = [_]f64{5} ** soil_species_count;
    var workspace = SurfaceFluxWorkspace{
        .surface_discharge_claimed = true,
        .litter_flux_mol_per_step = &litter,
        .soil_flux_mol_per_step = &soil,
    };
    try reset(&workspace);
    try std.testing.expect(!workspace.surface_discharge_claimed);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** litter_species_count), &litter);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** soil_species_count), &soil);
}

test "snowmelt surface workspaces preserve differing compact topologies" {
    try std.testing.expectEqual(@as(usize, 42), litter_species_count);
    try std.testing.expectEqual(@as(usize, 50), soil_species_count);
}

test "late soil dimension failure is atomic" {
    var litter = [_]f64{3} ** litter_species_count;
    var short_soil = [_]f64{5} ** (soil_species_count - 1);
    var workspace = SurfaceFluxWorkspace{
        .surface_discharge_claimed = true,
        .litter_flux_mol_per_step = &litter,
        .soil_flux_mol_per_step = &short_soil,
    };
    try std.testing.expectError(
        error.SnowmeltSoilSoluteFluxDimensionMismatch,
        reset(&workspace),
    );
    try std.testing.expect(workspace.surface_discharge_claimed);
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** litter_species_count), &litter);
    try std.testing.expectEqualSlices(f64, &([_]f64{5} ** (soil_species_count - 1)), &short_soil);
}

test "long litter topology fails before sentinel mutation" {
    var long_litter = [_]f64{3} ** (litter_species_count + 1);
    var soil = [_]f64{5} ** soil_species_count;
    var workspace = SurfaceFluxWorkspace{
        .surface_discharge_claimed = true,
        .litter_flux_mol_per_step = &long_litter,
        .soil_flux_mol_per_step = &soil,
    };
    try std.testing.expectError(
        error.SnowmeltLitterSoluteFluxDimensionMismatch,
        reset(&workspace),
    );
    try std.testing.expect(workspace.surface_discharge_claimed);
    try std.testing.expectEqualSlices(f64, &([_]f64{5} ** soil_species_count), &soil);
}
