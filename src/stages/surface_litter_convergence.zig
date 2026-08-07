//! Surface litter chemistry convergence and its denitrification map.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
const tile_kernels = @import("tile_kernels.zig");
pub const SurfaceDenitrificationRespirationMap = struct {
    destination_g_c: []f64,
    source_g_c: []const f64,
};

pub fn mapSurfaceDenitrificationRespiration(
    context: *SurfaceDenitrificationRespirationMap,
    cells: ecosys.compute.CellRange,
) !void {
    for (cells.first..cells.end) |cell| {
        for (0..ecosys.surface_microbial_respiration_step.litter_complex_count) |complex| {
            const compact =
                cell * ecosys.surface_microbial_respiration_step.litter_complex_count +
                complex;
            const unit =
                cell * ecosys.surface_microbial_respiration_step.unit_count_per_cell +
                complex *
                    ecosys.surface_microbial_respiration_step.source_population_count +
                ecosys.surface_denitrification_step.denitrifier_population;
            context.destination_g_c[unit] = context.source_g_c[compact];
        }
    }
}

pub fn convergeSurfaceLitterChemistry(context: anytype) !void {
    {
        const reaction_parameters = context.chemistry_reaction_parameters.*;
        context.surface_litter_chemistry_diagnostics.reset();
        var litter_chemistry_context: ecosys.surface_litter_chemistry_step.ApplyContext = .{
            .state = context.surface_litter_chemistry,
            .surface_organic = context.surface_organic,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .chemistry_parameters = reaction_parameters,
            .cation_selectivity_by_cell = context.surface_litter_cation_selectivity,
            .litter_dry_mass_megagrams = context.surface_litter_geometry.dry_mass_megagrams,
            .dynamic_salts = context.runscript.dynamic_plant_salts,
            .solver_options = .{
                .absolute_tolerance = context.config.absolute_tolerance,
                .relative_tolerance = context.config.relative_tolerance,
                .picard_relaxation = context.config.picard_relaxation,
                // Surface litter is seeded from topsoil and must close the
                // Dynamic salts use STARTE's MRXN=1000 equilibrium ceiling.
                // The fixed-pH ISALTG=0 source branch is outside that loop;
                // its already-hourly kinetic coefficients are applied once.
                .max_iterations = context.iteration_limits.initial_solute_reaction_max_iterations,
            },
            .diagnostics = context.surface_litter_chemistry_diagnostics,
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &litter_chemistry_context, ecosys.surface_litter_chemistry_step.applyTile);
        try ecosys.surface_litter_chemistry_step.publishAcceptedCarbonDioxideChanges(
            context.litter_gas_transport,
            context.surface_litter_chemistry_diagnostics,
            12,
        );
    }
}
