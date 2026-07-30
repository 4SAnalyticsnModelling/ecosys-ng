const std = @import("std");
const GridState = @import("grid.zig").GridState;
const thermal_module = @import("soil_thermal.zig");

pub const Parameters = struct {
    liquid_water_heat_capacity_mj_per_m3_k: f64,
    ice_heat_capacity_mj_per_m3_k: f64,
    minimum_heat_capacity_mj_per_k: f64,
};

const Candidate = struct {
    source_matrix_liquid_water_m3: f64,
    destination_matrix_liquid_water_m3: f64,
    source_matrix_ice_water_m3: f64,
    destination_matrix_ice_water_m3: f64,
    source_matrix_pore_capacity_m3: f64,
    destination_matrix_pore_capacity_m3: f64,
    source_matrix_air_volume_m3: f64,
    destination_matrix_air_volume_m3: f64,
    source_layer_volume_m3: f64,
    destination_layer_volume_m3: f64,
    source_dry_heat_capacity_mj_per_k: f64,
    destination_dry_heat_capacity_mj_per_k: f64,
    source_total_heat_capacity_mj_per_k: f64,
    destination_total_heat_capacity_mj_per_k: f64,
    source_temperature_k: f64,
    destination_temperature_k: f64,
};

/// REDIST pond water/ice/heat remap. Only matrix-domain storage moves; the
/// macropore domain remains attached to its layer exactly as in the source.
pub fn transferLayerFraction(
    grid: *GridState,
    thermal: *thermal_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    parameters: Parameters,
) !void {
    const candidate = try calculate(grid, thermal, source, destination, fraction, parameters);

    grid.matrix_liquid_water_m3[source] = candidate.source_matrix_liquid_water_m3;
    grid.matrix_liquid_water_m3[destination] = candidate.destination_matrix_liquid_water_m3;
    grid.matrix_ice_water_m3[source] = candidate.source_matrix_ice_water_m3;
    grid.matrix_ice_water_m3[destination] = candidate.destination_matrix_ice_water_m3;
    grid.matrix_pore_capacity_m3[source] = candidate.source_matrix_pore_capacity_m3;
    grid.matrix_pore_capacity_m3[destination] = candidate.destination_matrix_pore_capacity_m3;
    grid.matrix_air_volume_m3[source] = candidate.source_matrix_air_volume_m3;
    grid.matrix_air_volume_m3[destination] = candidate.destination_matrix_air_volume_m3;
    thermal.layer_volume_m3[source] = candidate.source_layer_volume_m3;
    thermal.layer_volume_m3[destination] = candidate.destination_layer_volume_m3;

    commitDerived(grid, thermal, source, candidate.source_dry_heat_capacity_mj_per_k, candidate.source_total_heat_capacity_mj_per_k, candidate.source_temperature_k);
    commitDerived(grid, thermal, destination, candidate.destination_dry_heat_capacity_mj_per_k, candidate.destination_total_heat_capacity_mj_per_k, candidate.destination_temperature_k);
}

pub fn validateLayerFraction(
    grid: *const GridState,
    thermal: *const thermal_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    parameters: Parameters,
) !void {
    _ = try calculate(grid, thermal, source, destination, fraction, parameters);
}

fn calculate(
    grid: *const GridState,
    thermal: *const thermal_module.State,
    source: usize,
    destination: usize,
    fraction: f64,
    parameters: Parameters,
) !Candidate {
    if (source >= grid.layer_count or destination >= grid.layer_count or source >= thermal.layer_volume_m3.len or destination >= thermal.layer_volume_m3.len or source == destination) return error.WaterHeatLayerRemapIndexOutOfBounds;
    inline for (.{ fraction, parameters.liquid_water_heat_capacity_mj_per_m3_k, parameters.ice_heat_capacity_mj_per_m3_k, parameters.minimum_heat_capacity_mj_per_k }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteWaterHeatLayerRemapInput;
    }
    if (fraction < 0 or fraction > 1 or parameters.liquid_water_heat_capacity_mj_per_m3_k <= 0 or parameters.ice_heat_capacity_mj_per_m3_k <= 0 or parameters.minimum_heat_capacity_mj_per_k < 0) return error.InvalidWaterHeatLayerRemapInput;

    const remaining = 1 - fraction;
    inline for (.{
        grid.matrix_liquid_water_m3[source],                 grid.matrix_liquid_water_m3[destination],
        grid.matrix_ice_water_m3[source],                    grid.matrix_ice_water_m3[destination],
        grid.macropore_liquid_water_m3[source],              grid.macropore_liquid_water_m3[destination],
        grid.macropore_ice_water_m3[source],                 grid.macropore_ice_water_m3[destination],
        grid.matrix_pore_capacity_m3[source],                grid.matrix_pore_capacity_m3[destination],
        grid.matrix_air_volume_m3[source],                   grid.matrix_air_volume_m3[destination],
        thermal.layer_volume_m3[source],                     thermal.layer_volume_m3[destination],
        thermal.dry_solid_heat_capacity_mj_per_m3_k[source], thermal.dry_solid_heat_capacity_mj_per_m3_k[destination],
        thermal.total_heat_capacity_mj_per_m3_k[source],     thermal.total_heat_capacity_mj_per_m3_k[destination],
        grid.soil_temperature_k[source],                     grid.soil_temperature_k[destination],
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidWaterHeatLayerRemapState;

    const source_volume_before = thermal.layer_volume_m3[source];
    const destination_volume_before = thermal.layer_volume_m3[destination];
    if (source_volume_before <= 0 or destination_volume_before <= 0 or grid.soil_temperature_k[source] <= 0 or grid.soil_temperature_k[destination] <= 0) return error.InvalidWaterHeatLayerRemapState;
    const source_dry_before = thermal.dry_solid_heat_capacity_mj_per_m3_k[source] * source_volume_before;
    const destination_dry_before = thermal.dry_solid_heat_capacity_mj_per_m3_k[destination] * destination_volume_before;
    const source_heat_before = thermal.total_heat_capacity_mj_per_m3_k[source] * source_volume_before;
    const destination_heat_before = thermal.total_heat_capacity_mj_per_m3_k[destination] * destination_volume_before;
    const source_energy_before = source_heat_before * grid.soil_temperature_k[source];
    const destination_energy_before = destination_heat_before * grid.soil_temperature_k[destination];

    var result: Candidate = .{
        .source_matrix_liquid_water_m3 = remaining * grid.matrix_liquid_water_m3[source],
        .destination_matrix_liquid_water_m3 = grid.matrix_liquid_water_m3[destination] + fraction * grid.matrix_liquid_water_m3[source],
        .source_matrix_ice_water_m3 = remaining * grid.matrix_ice_water_m3[source],
        .destination_matrix_ice_water_m3 = grid.matrix_ice_water_m3[destination] + fraction * grid.matrix_ice_water_m3[source],
        .source_matrix_pore_capacity_m3 = remaining * grid.matrix_pore_capacity_m3[source],
        .destination_matrix_pore_capacity_m3 = grid.matrix_pore_capacity_m3[destination] + fraction * grid.matrix_pore_capacity_m3[source],
        .source_matrix_air_volume_m3 = remaining * grid.matrix_air_volume_m3[source],
        .destination_matrix_air_volume_m3 = grid.matrix_air_volume_m3[destination] + fraction * grid.matrix_air_volume_m3[source],
        .source_layer_volume_m3 = remaining * source_volume_before,
        .destination_layer_volume_m3 = destination_volume_before + fraction * source_volume_before,
        .source_dry_heat_capacity_mj_per_k = remaining * source_dry_before,
        .destination_dry_heat_capacity_mj_per_k = destination_dry_before + fraction * source_dry_before,
        .source_total_heat_capacity_mj_per_k = undefined,
        .destination_total_heat_capacity_mj_per_k = undefined,
        .source_temperature_k = undefined,
        .destination_temperature_k = undefined,
    };
    result.source_total_heat_capacity_mj_per_k = result.source_dry_heat_capacity_mj_per_k +
        parameters.liquid_water_heat_capacity_mj_per_m3_k * (result.source_matrix_liquid_water_m3 + grid.macropore_liquid_water_m3[source]) +
        parameters.ice_heat_capacity_mj_per_m3_k * (result.source_matrix_ice_water_m3 + grid.macropore_ice_water_m3[source]);
    result.destination_total_heat_capacity_mj_per_k = result.destination_dry_heat_capacity_mj_per_k +
        parameters.liquid_water_heat_capacity_mj_per_m3_k * (result.destination_matrix_liquid_water_m3 + grid.macropore_liquid_water_m3[destination]) +
        parameters.ice_heat_capacity_mj_per_m3_k * (result.destination_matrix_ice_water_m3 + grid.macropore_ice_water_m3[destination]);
    const source_energy_after = remaining * source_energy_before;
    const destination_energy_after = destination_energy_before + fraction * source_energy_before;
    result.destination_temperature_k = if (result.destination_total_heat_capacity_mj_per_k > parameters.minimum_heat_capacity_mj_per_k) destination_energy_after / result.destination_total_heat_capacity_mj_per_k else grid.soil_temperature_k[source];
    result.source_temperature_k = if (result.source_total_heat_capacity_mj_per_k > parameters.minimum_heat_capacity_mj_per_k) source_energy_after / result.source_total_heat_capacity_mj_per_k else result.destination_temperature_k;
    inline for (.{ result.source_total_heat_capacity_mj_per_k, result.destination_total_heat_capacity_mj_per_k, result.source_temperature_k, result.destination_temperature_k }) |value| {
        if (!std.math.isFinite(value) or value < 0) return error.InvalidWaterHeatLayerRemapResult;
    }
    if (result.source_layer_volume_m3 == 0 and result.source_total_heat_capacity_mj_per_k > parameters.minimum_heat_capacity_mj_per_k) return error.InvalidWaterHeatLayerRemapResult;
    return result;
}

fn commitDerived(grid: *GridState, thermal: *thermal_module.State, index: usize, dry_heat_capacity_mj_per_k: f64, total_heat_capacity_mj_per_k: f64, temperature_k: f64) void {
    grid.liquid_water_m3[index] = grid.matrix_liquid_water_m3[index] + grid.macropore_liquid_water_m3[index];
    grid.ice_water_m3[index] = grid.matrix_ice_water_m3[index] + grid.macropore_ice_water_m3[index];
    grid.air_volume_m3[index] = grid.matrix_air_volume_m3[index] + grid.macropore_air_volume_m3[index];
    grid.soil_temperature_k[index] = temperature_k;
    const volume = thermal.layer_volume_m3[index];
    thermal.dry_solid_heat_capacity_mj_per_m3_k[index] = if (volume > 0) dry_heat_capacity_mj_per_k / volume else 0;
    thermal.total_heat_capacity_mj_per_m3_k[index] = if (volume > 0) total_heat_capacity_mj_per_k / volume else 0;
    thermal.porosity_fraction[index] = if (volume > 0) (grid.matrix_pore_capacity_m3[index] + grid.macropore_pore_capacity_m3[index]) / volume else 0;
}

test "REDIST water heat remap conserves matrix stores and energy while leaving macropores" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 4 });
    var grid = try GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var thermal: thermal_module.State = undefined;
    var layer_volume_m3 = [_]f64{ 2, 3 };
    var dry_heat_capacity_mj_per_m3_k = [_]f64{ 1, 2 };
    var total_heat_capacity_mj_per_m3_k = [_]f64{ 2, 3 };
    var porosity_fraction = [_]f64{ 0.5, 0.5 };
    thermal.layer_volume_m3 = &layer_volume_m3;
    thermal.dry_solid_heat_capacity_mj_per_m3_k = &dry_heat_capacity_mj_per_m3_k;
    thermal.total_heat_capacity_mj_per_m3_k = &total_heat_capacity_mj_per_m3_k;
    thermal.porosity_fraction = &porosity_fraction;
    grid.matrix_liquid_water_m3[0] = 0.4;
    grid.matrix_liquid_water_m3[1] = 0.3;
    grid.macropore_liquid_water_m3[0] = 0.1;
    grid.macropore_liquid_water_m3[1] = 0.2;
    grid.matrix_ice_water_m3[0] = 0.2;
    grid.matrix_ice_water_m3[1] = 0.1;
    grid.matrix_pore_capacity_m3[0] = 1;
    grid.matrix_pore_capacity_m3[1] = 1.5;
    grid.matrix_air_volume_m3[0] = 0.4;
    grid.matrix_air_volume_m3[1] = 0.8;
    grid.soil_temperature_k[0] = 280;
    grid.soil_temperature_k[1] = 290;
    const initial_energy = 2 * 2 * 280 + 3 * 3 * 290;
    try transferLayerFraction(&grid, &thermal, 0, 1, 0.25, .{ .liquid_water_heat_capacity_mj_per_m3_k = 4.19, .ice_heat_capacity_mj_per_m3_k = 1.9274, .minimum_heat_capacity_mj_per_k = 0 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), grid.matrix_liquid_water_m3[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), grid.matrix_liquid_water_m3[1], 1e-14);
    try std.testing.expectEqual(@as(f64, 0.1), grid.macropore_liquid_water_m3[0]);
    const final_energy = thermal.total_heat_capacity_mj_per_m3_k[0] * thermal.layer_volume_m3[0] * grid.soil_temperature_k[0] +
        thermal.total_heat_capacity_mj_per_m3_k[1] * thermal.layer_volume_m3[1] * grid.soil_temperature_k[1];
    try std.testing.expectApproxEqAbs(initial_energy, final_energy, 1e-9);
}
