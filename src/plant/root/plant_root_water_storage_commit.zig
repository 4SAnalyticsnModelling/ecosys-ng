const std = @import("std");
const GridState = @import("../../state/grid.zig").GridState;
const root_system = @import("plant_root_system.zig");
const RootState = root_system.State;
const chemistry_rebase = @import("../../soil/chemistry/water_carrier_rebase.zig");
const ChemistryState = @import("../../soil/solute/chemistry_state.zig").State;

/// Commits accepted UPTAKE root-water fluxes to soil matrix storage. Root
/// fluxes use the source convention: negative removes water from soil and
/// positive returns water to soil.
pub fn commit(
    roots: *const RootState,
    grid: *GridState,
    chemistry: *ChemistryState,
    biological_domain_count_by_plant: []const u8,
) !void {
    if (grid.cell_count == 0 or roots.plant_count % grid.cell_count != 0 or
        roots.soil_layer_count != grid.soil_layer_capacity or
        biological_domain_count_by_plant.len != roots.plant_count or
        chemistry.cell_count != grid.layer_count)
        return error.PlantRootWaterStorageDimensionMismatch;
    const species_count = roots.plant_count / grid.cell_count;

    // Validate every accepted layer transaction before changing storage.
    for (0..grid.cell_count) |cell| for (0..grid.active_soil_layer_count[cell]) |layer| {
        const soil = try grid.layerIndex(cell, layer);
        const change_m3 = try layerWaterChange(
            roots,
            biological_domain_count_by_plant,
            species_count,
            cell,
            layer,
        );
        const next = grid.matrix_liquid_water_m3[soil] + change_m3;
        if (!std.math.isFinite(next) or next < 0)
            return error.PlantRootWaterUptakeExceedsSoilStorage;
    };

    for (0..grid.cell_count) |cell| for (0..grid.active_soil_layer_count[cell]) |layer| {
        const soil = try grid.layerIndex(cell, layer);
        const old_water_m3 = grid.matrix_liquid_water_m3[soil];
        const change_m3 = try layerWaterChange(
            roots,
            biological_domain_count_by_plant,
            species_count,
            cell,
            layer,
        );
        const new_water_m3 = old_water_m3 + change_m3;
        try chemistry_rebase.rebaseLayer(chemistry, soil, old_water_m3, new_water_m3);
        grid.matrix_liquid_water_m3[soil] = new_water_m3;
    };
}

fn layerWaterChange(
    roots: *const RootState,
    domain_count_by_plant: []const u8,
    species_count: usize,
    cell: usize,
    layer: usize,
) !f64 {
    var change_m3: f64 = 0;
    for (0..species_count) |species| {
        const plant = cell * species_count + species;
        const domain_count = domain_count_by_plant[plant];
        if (domain_count == 0 or domain_count > root_system.biological_domain_count)
            return error.InvalidPlantRootBiologicalDomainCount;
        for (0..domain_count) |domain| {
            const root = try roots.layerIndex(plant, domain, layer);
            const flux = roots.water_uptake_m3_per_h[root];
            if (!std.math.isFinite(flux)) return error.NonFinitePlantRootWaterFlux;
            change_m3 += flux;
        }
    }
    if (!std.math.isFinite(change_m3)) return error.NonFinitePlantRootWaterFlux;
    return change_m3;
}

test "root uptake removes soil water and preserves dissolved moles" {
    const config = @import("../../core/config.zig");
    const cfg = try config.SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 },
    );
    var grid = try GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var chemistry = try ChemistryState.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    grid.matrix_liquid_water_m3[0] = 10;
    grid.matrix_liquid_water_m3[1] = 8;
    chemistry.aqueous[0].nitrate_non_band = 2;
    roots.water_uptake_m3_per_h[try roots.layerIndex(0, 0, 0)] = -3;
    roots.water_uptake_m3_per_h[try roots.layerIndex(0, 0, 1)] = -1;
    try commit(&roots, &grid, &chemistry, &.{1});
    try std.testing.expectEqual(@as(f64, 7), grid.matrix_liquid_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 7), grid.matrix_liquid_water_m3[1]);
    try std.testing.expectApproxEqAbs(20.0 / 7.0, chemistry.aqueous[0].nitrate_non_band, 1e-14);
}
