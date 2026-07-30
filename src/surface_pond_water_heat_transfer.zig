const std = @import("std");
const GridState = @import("grid.zig").GridState;
const thermal_module = @import("soil_thermal.zig");
const surface_geometry_module = @import("surface_litter_geometry_step.zig");

pub const Parameters = struct {
    dry_organic_heat_capacity_mj_per_g_c_k: f64,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
    ice_heat_capacity_mj_per_m3_k: f64,
    minimum_heat_capacity_mj_per_k: f64,
};

pub const Owners = struct {
    surface_liquid_water_m3: []f64,
    surface_ice_m3: []f64,
    surface_temperature_k: []f64,
    surface_geometry: *surface_geometry_module.State,
    grid: *GridState,
    soil_thermal: *thermal_module.State,
};

pub const Inputs = struct {
    cell: usize,
    destination_soil_layer: usize,
    fraction: f64,
    surface_organic_carbon_g_c: f64,
    parameters: Parameters,
};

const Candidate = struct {
    surface_liquid_water_m3: f64,
    surface_ice_m3: f64,
    surface_temperature_k: f64,
    surface_water_retention_capacity_m3: f64,
    surface_dry_litter_volume_m3: f64,
    surface_total_volume_m3: f64,
    surface_dry_mass_Mg: f64,
    surface_pore_volume_m3: f64,
    surface_air_volume_m3: f64,
    soil_matrix_liquid_water_m3: f64,
    soil_matrix_ice_water_m3: f64,
    soil_matrix_pore_capacity_m3: f64,
    soil_matrix_air_volume_m3: f64,
    soil_layer_volume_m3: f64,
    soil_dry_heat_capacity_mj_per_k: f64,
    soil_total_heat_capacity_mj_per_k: f64,
    soil_temperature_k: f64,
};

pub fn transferSurfaceFractionToSoil(owners: Owners, inputs: Inputs) !void {
    const next = try calculate(owners, inputs);
    const destination = inputs.cell * owners.grid.soil_layer_capacity + inputs.destination_soil_layer;
    owners.surface_liquid_water_m3[inputs.cell] = next.surface_liquid_water_m3;
    owners.surface_ice_m3[inputs.cell] = next.surface_ice_m3;
    owners.surface_temperature_k[inputs.cell] = next.surface_temperature_k;
    owners.surface_geometry.water_retention_capacity_m3[inputs.cell] = next.surface_water_retention_capacity_m3;
    owners.surface_geometry.dry_litter_volume_m3[inputs.cell] = next.surface_dry_litter_volume_m3;
    owners.surface_geometry.expanded_total_volume_m3[inputs.cell] = next.surface_total_volume_m3;
    owners.surface_geometry.dry_mass_Mg[inputs.cell] = next.surface_dry_mass_Mg;
    owners.surface_geometry.pore_volume_m3[inputs.cell] = next.surface_pore_volume_m3;
    owners.surface_geometry.air_volume_m3[inputs.cell] = next.surface_air_volume_m3;
    owners.grid.matrix_liquid_water_m3[destination] = next.soil_matrix_liquid_water_m3;
    owners.grid.matrix_ice_water_m3[destination] = next.soil_matrix_ice_water_m3;
    owners.grid.matrix_pore_capacity_m3[destination] = next.soil_matrix_pore_capacity_m3;
    owners.grid.matrix_air_volume_m3[destination] = next.soil_matrix_air_volume_m3;
    owners.grid.liquid_water_m3[destination] = next.soil_matrix_liquid_water_m3 + owners.grid.macropore_liquid_water_m3[destination];
    owners.grid.ice_water_m3[destination] = next.soil_matrix_ice_water_m3 + owners.grid.macropore_ice_water_m3[destination];
    owners.grid.air_volume_m3[destination] = next.soil_matrix_air_volume_m3 + owners.grid.macropore_air_volume_m3[destination];
    owners.grid.soil_temperature_k[destination] = next.soil_temperature_k;
    owners.soil_thermal.layer_volume_m3[destination] = next.soil_layer_volume_m3;
    owners.soil_thermal.dry_solid_heat_capacity_mj_per_m3_k[destination] = next.soil_dry_heat_capacity_mj_per_k / next.soil_layer_volume_m3;
    owners.soil_thermal.total_heat_capacity_mj_per_m3_k[destination] = next.soil_total_heat_capacity_mj_per_k / next.soil_layer_volume_m3;
    owners.soil_thermal.porosity_fraction[destination] = (next.soil_matrix_pore_capacity_m3 + owners.grid.macropore_pore_capacity_m3[destination]) / next.soil_layer_volume_m3;
}

pub fn validateSurfaceFractionToSoil(owners: Owners, inputs: Inputs) !void {
    _ = try calculate(owners, inputs);
}

fn calculate(owners: Owners, inputs: Inputs) !Candidate {
    if (inputs.cell >= owners.grid.cell_count or inputs.destination_soil_layer >= owners.grid.soil_layer_capacity or inputs.cell >= owners.surface_liquid_water_m3.len or inputs.cell >= owners.surface_ice_m3.len or inputs.cell >= owners.surface_temperature_k.len or inputs.cell >= owners.surface_geometry.cell_count) return error.SurfacePondWaterHeatIndexOutOfBounds;
    const destination = inputs.cell * owners.grid.soil_layer_capacity + inputs.destination_soil_layer;
    if (destination >= owners.grid.layer_count or destination >= owners.soil_thermal.layer_volume_m3.len) return error.SurfacePondWaterHeatDimensionMismatch;
    inline for (.{
        inputs.fraction,                                                  inputs.surface_organic_carbon_g_c,
        inputs.parameters.dry_organic_heat_capacity_mj_per_g_c_k,         inputs.parameters.liquid_water_heat_capacity_mj_per_m3_k,
        inputs.parameters.ice_heat_capacity_mj_per_m3_k,                  inputs.parameters.minimum_heat_capacity_mj_per_k,
        owners.surface_liquid_water_m3[inputs.cell],                      owners.surface_ice_m3[inputs.cell],
        owners.surface_temperature_k[inputs.cell],                        owners.surface_geometry.water_retention_capacity_m3[inputs.cell],
        owners.surface_geometry.dry_litter_volume_m3[inputs.cell],        owners.surface_geometry.expanded_total_volume_m3[inputs.cell],
        owners.surface_geometry.dry_mass_Mg[inputs.cell],                 owners.surface_geometry.pore_volume_m3[inputs.cell],
        owners.surface_geometry.air_volume_m3[inputs.cell],               owners.grid.matrix_liquid_water_m3[destination],
        owners.grid.matrix_ice_water_m3[destination],                     owners.grid.matrix_pore_capacity_m3[destination],
        owners.grid.matrix_air_volume_m3[destination],                    owners.grid.soil_temperature_k[destination],
        owners.soil_thermal.layer_volume_m3[destination],                 owners.soil_thermal.dry_solid_heat_capacity_mj_per_m3_k[destination],
        owners.soil_thermal.total_heat_capacity_mj_per_m3_k[destination],
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfacePondWaterHeatState;
    if (inputs.fraction < 0 or inputs.fraction > 1 or inputs.parameters.dry_organic_heat_capacity_mj_per_g_c_k <= 0 or inputs.parameters.liquid_water_heat_capacity_mj_per_m3_k <= 0 or inputs.parameters.ice_heat_capacity_mj_per_m3_k <= 0 or inputs.parameters.minimum_heat_capacity_mj_per_k < 0 or owners.surface_temperature_k[inputs.cell] <= 0 or owners.grid.soil_temperature_k[destination] <= 0 or owners.soil_thermal.layer_volume_m3[destination] <= 0) return error.InvalidSurfacePondWaterHeatState;

    const fraction = inputs.fraction;
    const remaining = 1 - fraction;
    const surface_dry_heat_capacity = inputs.parameters.dry_organic_heat_capacity_mj_per_g_c_k * inputs.surface_organic_carbon_g_c;
    const surface_total_heat_capacity = surface_dry_heat_capacity +
        inputs.parameters.liquid_water_heat_capacity_mj_per_m3_k * owners.surface_liquid_water_m3[inputs.cell] +
        inputs.parameters.ice_heat_capacity_mj_per_m3_k * owners.surface_ice_m3[inputs.cell];
    const soil_dry_heat_capacity = owners.soil_thermal.dry_solid_heat_capacity_mj_per_m3_k[destination] * owners.soil_thermal.layer_volume_m3[destination];
    const soil_heat_capacity_before = owners.soil_thermal.total_heat_capacity_mj_per_m3_k[destination] * owners.soil_thermal.layer_volume_m3[destination];
    const moved_surface_energy = fraction * surface_total_heat_capacity * owners.surface_temperature_k[inputs.cell];

    var result: Candidate = .{
        .surface_liquid_water_m3 = remaining * owners.surface_liquid_water_m3[inputs.cell],
        .surface_ice_m3 = remaining * owners.surface_ice_m3[inputs.cell],
        .surface_temperature_k = undefined,
        .surface_water_retention_capacity_m3 = remaining * owners.surface_geometry.water_retention_capacity_m3[inputs.cell],
        .surface_dry_litter_volume_m3 = remaining * owners.surface_geometry.dry_litter_volume_m3[inputs.cell],
        .surface_total_volume_m3 = remaining * owners.surface_geometry.expanded_total_volume_m3[inputs.cell],
        .surface_dry_mass_Mg = remaining * owners.surface_geometry.dry_mass_Mg[inputs.cell],
        .surface_pore_volume_m3 = remaining * owners.surface_geometry.pore_volume_m3[inputs.cell],
        .surface_air_volume_m3 = remaining * owners.surface_geometry.air_volume_m3[inputs.cell],
        .soil_matrix_liquid_water_m3 = owners.grid.matrix_liquid_water_m3[destination] + fraction * owners.surface_liquid_water_m3[inputs.cell],
        .soil_matrix_ice_water_m3 = owners.grid.matrix_ice_water_m3[destination] + fraction * owners.surface_ice_m3[inputs.cell],
        .soil_matrix_pore_capacity_m3 = undefined,
        .soil_matrix_air_volume_m3 = undefined,
        .soil_layer_volume_m3 = owners.soil_thermal.layer_volume_m3[destination] + fraction * owners.surface_geometry.expanded_total_volume_m3[inputs.cell],
        .soil_dry_heat_capacity_mj_per_k = soil_dry_heat_capacity + fraction * surface_dry_heat_capacity,
        .soil_total_heat_capacity_mj_per_k = undefined,
        .soil_temperature_k = undefined,
    };
    // A pond-to-soil transition moves all water carried by the incorporated
    // surface fraction. Wet residue can contain more water than its preceding
    // dry-litter pore estimate, so the new soil pore geometry must represent
    // at least that conserved occupancy. Otherwise the transaction silently
    // creates an over-saturated top layer that fails in the following hour.
    result.soil_matrix_pore_capacity_m3 = @max(
        owners.grid.matrix_pore_capacity_m3[destination] +
            fraction * owners.surface_geometry.pore_volume_m3[inputs.cell],
        result.soil_matrix_liquid_water_m3 +
            result.soil_matrix_ice_water_m3,
    );
    if (result.soil_matrix_pore_capacity_m3 >
        result.soil_layer_volume_m3 + 1e-12)
        return error.SurfacePondWaterExceedsTransferredLayerVolume;
    result.soil_matrix_air_volume_m3 = @max(
        0,
        result.soil_matrix_pore_capacity_m3 -
            result.soil_matrix_liquid_water_m3 -
            result.soil_matrix_ice_water_m3,
    );
    result.soil_total_heat_capacity_mj_per_k = result.soil_dry_heat_capacity_mj_per_k +
        inputs.parameters.liquid_water_heat_capacity_mj_per_m3_k * (result.soil_matrix_liquid_water_m3 + owners.grid.macropore_liquid_water_m3[destination]) +
        inputs.parameters.ice_heat_capacity_mj_per_m3_k * (result.soil_matrix_ice_water_m3 + owners.grid.macropore_ice_water_m3[destination]);
    result.soil_temperature_k = if (result.soil_total_heat_capacity_mj_per_k > inputs.parameters.minimum_heat_capacity_mj_per_k)
        (soil_heat_capacity_before * owners.grid.soil_temperature_k[destination] + moved_surface_energy) / result.soil_total_heat_capacity_mj_per_k
    else
        owners.surface_temperature_k[inputs.cell];
    const remaining_surface_heat_capacity = remaining * surface_total_heat_capacity;
    result.surface_temperature_k = if (remaining_surface_heat_capacity > inputs.parameters.minimum_heat_capacity_mj_per_k)
        (remaining * surface_total_heat_capacity * owners.surface_temperature_k[inputs.cell]) / remaining_surface_heat_capacity
    else
        result.soil_temperature_k;
    inline for (.{ result.soil_total_heat_capacity_mj_per_k, result.soil_temperature_k, result.surface_temperature_k, result.soil_layer_volume_m3 }) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSurfacePondWaterHeatResult;
    return result;
}

test "surface pond water heat and matrix geometry transfer conserve energy" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 4 });
    var grid = try GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var geometry = try surface_geometry_module.State.init(std.testing.allocator, 1);
    defer geometry.deinit();
    var thermal: thermal_module.State = undefined;
    var layer_volume = [_]f64{2};
    var dry_capacity = [_]f64{1};
    var total_capacity = [_]f64{2};
    var porosity = [_]f64{0.5};
    thermal.layer_volume_m3 = &layer_volume;
    thermal.dry_solid_heat_capacity_mj_per_m3_k = &dry_capacity;
    thermal.total_heat_capacity_mj_per_m3_k = &total_capacity;
    thermal.porosity_fraction = &porosity;
    var surface_water = [_]f64{0.4};
    var surface_ice = [_]f64{0.2};
    grid.surface_temperature_k[0] = 280;
    grid.soil_temperature_k[0] = 290;
    grid.matrix_liquid_water_m3[0] = 0.3;
    grid.matrix_ice_water_m3[0] = 0.1;
    geometry.dry_litter_volume_m3[0] = 0.2;
    geometry.expanded_total_volume_m3[0] = 1;
    geometry.water_retention_capacity_m3[0] = 0.3;
    geometry.dry_mass_Mg[0] = 0.01;
    geometry.pore_volume_m3[0] = 0.6;
    geometry.air_volume_m3[0] = 0.2;
    const parameters: Parameters = .{ .dry_organic_heat_capacity_mj_per_g_c_k = 2.496e-6, .liquid_water_heat_capacity_mj_per_m3_k = 4.19, .ice_heat_capacity_mj_per_m3_k = 1.9274, .minimum_heat_capacity_mj_per_k = 0 };
    const surface_capacity = parameters.dry_organic_heat_capacity_mj_per_g_c_k * 100 + parameters.liquid_water_heat_capacity_mj_per_m3_k * 0.4 + parameters.ice_heat_capacity_mj_per_m3_k * 0.2;
    const initial_energy = surface_capacity * 280 + 4 * 290;
    const owners: Owners = .{ .surface_liquid_water_m3 = &surface_water, .surface_ice_m3 = &surface_ice, .surface_temperature_k = grid.surface_temperature_k, .surface_geometry = &geometry, .grid = &grid, .soil_thermal = &thermal };
    try transferSurfaceFractionToSoil(owners, .{ .cell = 0, .destination_soil_layer = 0, .fraction = 0.5, .surface_organic_carbon_g_c = 100, .parameters = parameters });
    const remaining_surface_capacity = 0.5 * surface_capacity;
    const final_energy = remaining_surface_capacity * grid.surface_temperature_k[0] + thermal.total_heat_capacity_mj_per_m3_k[0] * thermal.layer_volume_m3[0] * grid.soil_temperature_k[0];
    try std.testing.expectApproxEqAbs(initial_energy, final_energy, 1e-9);
    try std.testing.expectEqual(@as(f64, 0.2), surface_water[0]);
    try std.testing.expectEqual(@as(f64, 0.5), geometry.expanded_total_volume_m3[0]);
    try std.testing.expectEqual(@as(f64, 0.15), geometry.water_retention_capacity_m3[0]);
    try std.testing.expectEqual(@as(f64, 0.005), geometry.dry_mass_Mg[0]);
}
