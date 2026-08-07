const std = @import("std");
const CellRange = @import("../../core/compute.zig").CellRange;
const GridState = @import("../../state/grid.zig").GridState;
const SoilThermalState = @import("thermal.zig").State;
const SurfaceTemperatureState = @import("../../surface/temperature_solver.zig").State;

/// Heap-backed workspace and diagnostics for conservative vertical soil heat
/// transport. Each layer slot stores the downward flux through its lower face.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    soil_layer_capacity: usize,
    downward_conductive_heat_flux_megajoules_per_m2: []f64,
    heat_balance_residual_megajoules_per_m2: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, soil_layer_capacity: usize) !State {
        const layer_count = try std.math.mul(usize, cell_count, soil_layer_capacity);
        const flux = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(flux);
        const residual = try allocator.alloc(f64, cell_count);
        @memset(flux, 0);
        @memset(residual, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .soil_layer_capacity = soil_layer_capacity, .downward_conductive_heat_flux_megajoules_per_m2 = flux, .heat_balance_residual_megajoules_per_m2 = residual };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.downward_conductive_heat_flux_megajoules_per_m2);
        self.allocator.free(self.heat_balance_residual_megajoules_per_m2);
        self.* = undefined;
    }

    pub fn validateFinite(self: State) !void {
        for (self.downward_conductive_heat_flux_megajoules_per_m2, 0..) |value, index| if (!std.math.isFinite(value)) {
            std.log.err("non-finite soil heat flux: layer_slot={d} value={e}", .{ index, value });
            return error.NonFiniteSoilHeatFlux;
        };
        for (self.heat_balance_residual_megajoules_per_m2, 0..) |value, cell| if (!std.math.isFinite(value)) {
            std.log.err("non-finite soil heat balance: cell={d} residual_megajoules_per_m2={e}", .{ cell, value });
            return error.NonFiniteSoilHeatBalance;
        };
    }
};

pub const UpdateContext = struct {
    transport: *State,
    grid: *GridState,
    thermal: *const SoilThermalState,
    surface: *const SurfaceTemperatureState,
    timestep_hours: f64,
};

/// Advances independent soil columns with conservative, flux-limited heat
/// conduction. The limiter is the modern equivalent of the WATSUB available-
/// heat bound and prevents a face flux from reversing either layer gradient.
pub fn updateTile(context: *UpdateContext, range: CellRange) !void {
    const transport = context.transport;
    if (!std.math.isFinite(context.timestep_hours) or context.timestep_hours <= 0) return error.InvalidSoilHeatTimestep;
    if (range.end > transport.cell_count or context.grid.cell_count != transport.cell_count or context.thermal.cell_count != transport.cell_count or context.surface.cell_count != transport.cell_count or context.grid.soil_layer_capacity != transport.soil_layer_capacity or context.thermal.soil_layer_capacity != transport.soil_layer_capacity) return error.SoilHeatDimensionMismatch;

    for (range.first..range.end) |cell| {
        const base = cell * transport.soil_layer_capacity;
        const active_layers = context.grid.active_soil_layer_count[cell];
        if (active_layers == 0 or active_layers > transport.soil_layer_capacity) return error.InvalidActiveSoilLayerCount;
        var initial_sensible_heat_megajoules_per_m2: f64 = 0;
        for (0..active_layers) |layer| {
            const index = base + layer;
            initial_sensible_heat_megajoules_per_m2 += arealHeatCapacity(context.thermal, index) * context.grid.soil_temperature_k[index];
        }
        const top_capacity = arealHeatCapacity(context.thermal, base);
        // Surface solver reports positive flux toward the surface; soil receives
        // the opposite flux. Limit it to the top layer's available sensible heat.
        const requested_surface_to_soil = -context.surface.conductive_heat_flux_megajoules_per_m2[cell];
        const top_difference_k = context.grid.surface_temperature_k[cell] - context.grid.soil_temperature_k[base];
        const surface_limit = @abs(top_difference_k) * top_capacity;
        const surface_to_soil = std.math.clamp(requested_surface_to_soil, -surface_limit, surface_limit);

        for (0..active_layers) |layer| transport.downward_conductive_heat_flux_megajoules_per_m2[base + layer] = 0;
        for (0..active_layers - 1) |layer| {
            const upper = base + layer;
            const lower = upper + 1;
            const upper_half_distance_m = 0.5 * context.thermal.layer_thickness_m[upper];
            const lower_half_distance_m = 0.5 * context.thermal.layer_thickness_m[lower];
            const resistance_h_k_m2_per_megajoule = upper_half_distance_m / context.thermal.thermal_conductivity_m_megajoules_per_h_k[upper] + lower_half_distance_m / context.thermal.thermal_conductivity_m_megajoules_per_h_k[lower];
            if (!std.math.isFinite(resistance_h_k_m2_per_megajoule) or resistance_h_k_m2_per_megajoule <= 0) return error.InvalidSoilThermalResistance;
            const temperature_difference_k = context.grid.soil_temperature_k[upper] - context.grid.soil_temperature_k[lower];
            const requested_flux = context.timestep_hours * temperature_difference_k / resistance_h_k_m2_per_megajoule;
            const equilibration_limit = @abs(temperature_difference_k) / (1.0 / arealHeatCapacity(context.thermal, upper) + 1.0 / arealHeatCapacity(context.thermal, lower));
            transport.downward_conductive_heat_flux_megajoules_per_m2[upper] = std.math.clamp(requested_flux, -equilibration_limit, equilibration_limit);
        }

        var incoming_megajoules_per_m2 = surface_to_soil;
        for (0..active_layers) |layer| {
            const index = base + layer;
            const outgoing_megajoules_per_m2 = transport.downward_conductive_heat_flux_megajoules_per_m2[index];
            context.grid.soil_temperature_k[index] += (incoming_megajoules_per_m2 - outgoing_megajoules_per_m2) / arealHeatCapacity(context.thermal, index);
            if (!std.math.isFinite(context.grid.soil_temperature_k[index]) or context.grid.soil_temperature_k[index] <= 0) {
                std.log.err("invalid soil temperature after heat transport: cell={d} layer={d} temperature_k={e}", .{ cell, layer, context.grid.soil_temperature_k[index] });
                return error.InvalidSoilTemperature;
            }
            incoming_megajoules_per_m2 = outgoing_megajoules_per_m2;
        }
        var final_sensible_heat_megajoules_per_m2: f64 = 0;
        for (0..active_layers) |layer| {
            const index = base + layer;
            final_sensible_heat_megajoules_per_m2 += arealHeatCapacity(context.thermal, index) * context.grid.soil_temperature_k[index];
        }
        const residual_megajoules_per_m2 = final_sensible_heat_megajoules_per_m2 - initial_sensible_heat_megajoules_per_m2 - surface_to_soil;
        transport.heat_balance_residual_megajoules_per_m2[cell] = residual_megajoules_per_m2;
        const balance_scale = @max(1.0, @max(@abs(initial_sensible_heat_megajoules_per_m2), @abs(final_sensible_heat_megajoules_per_m2)));
        if (@abs(residual_megajoules_per_m2) > 1.0e-12 * balance_scale) {
            std.log.err("soil heat conservation failure: cell={d} residual_megajoules_per_m2={e} scale_megajoules_per_m2={e}", .{ cell, residual_megajoules_per_m2, balance_scale });
            return error.SoilHeatConservationFailure;
        }
    }
}

fn arealHeatCapacity(thermal: *const SoilThermalState, index: usize) f64 {
    return thermal.total_heat_capacity_megajoules_per_m3_k[index] * thermal.layer_thickness_m[index];
}
