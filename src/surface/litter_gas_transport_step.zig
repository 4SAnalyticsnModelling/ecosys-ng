const std = @import("std");
const gas = @import("../soil/gas/transport.zig");
const atmosphere = @import("../soil/gas/atmosphere_exchange.zig");
const solver = @import("../soil/gas/coupled_gas_solver.zig");
const litter_geometry = @import("litter_geometry_step.zig");
const precipitation = @import("precipitation.zig");

pub const RuntimeParameters = struct {
    reference_temperature_k: f64 = 298.15,
    temperature_exponent: f64 = 1.75,
    free_air_diffusivity_m2_per_h: [gas.species_count]f64 = .{ 4.68e-2, 7.80e-2, 6.43e-2, 5.57e-2, 5.57e-2, 6.67e-2, 5.57e-2 },
    penman_tortuosity: f64 = 0.66,
    minimum_air_fraction: f64 = 1e-12,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    boundaries: []atmosphere.Boundary,
    water_volume_m3: []f64,
    band_water_volume_m3: []f64,
    mass_solubility_ratio: []f64,
    gas_water_exchange_rate_per_step: []f64,
    band_gas_water_exchange_rate_per_step: []f64,
    bubbling_enabled: []bool,
    atmospheric_flux_g_per_h: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroSurfaceLitterGasCellCount;
        const components = try std.math.mul(usize, cell_count, gas.species_count);
        const boundaries = try allocator.alloc(atmosphere.Boundary, cell_count);
        errdefer allocator.free(boundaries);
        const water = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(water);
        const band_water = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(band_water);
        const solubility = try allocator.alloc(f64, components);
        errdefer allocator.free(solubility);
        const exchange = try allocator.alloc(f64, components);
        errdefer allocator.free(exchange);
        const band_exchange = try allocator.alloc(f64, components);
        errdefer allocator.free(band_exchange);
        const bubbling = try allocator.alloc(bool, cell_count);
        errdefer allocator.free(bubbling);
        const atmospheric_flux = try allocator.alloc(f64, components);
        errdefer allocator.free(atmospheric_flux);
        @memset(water, 0);
        @memset(band_water, 0);
        @memset(solubility, 0);
        @memset(exchange, 0);
        @memset(band_exchange, 0);
        @memset(bubbling, false);
        @memset(atmospheric_flux, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .boundaries = boundaries, .water_volume_m3 = water, .band_water_volume_m3 = band_water, .mass_solubility_ratio = solubility, .gas_water_exchange_rate_per_step = exchange, .band_gas_water_exchange_rate_per_step = band_exchange, .bubbling_enabled = bubbling, .atmospheric_flux_g_per_h = atmospheric_flux };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.atmospheric_flux_g_per_h);
        self.allocator.free(self.bubbling_enabled);
        self.allocator.free(self.band_gas_water_exchange_rate_per_step);
        self.allocator.free(self.gas_water_exchange_rate_per_step);
        self.allocator.free(self.mass_solubility_ratio);
        self.allocator.free(self.band_water_volume_m3);
        self.allocator.free(self.water_volume_m3);
        self.allocator.free(self.boundaries);
        self.* = undefined;
    }

    pub fn advance(
        self: *State,
        gas_state: *gas.State,
        geometry: *const litter_geometry.State,
        litter_water_m3: []const f64,
        litter_ice_m3: []const f64,
        cell_area_m2: []const f64,
        atmospheric_conductance_m3_per_step: []const f64,
        atmospheric_concentration_g_per_m3: [gas.species_count]f64,
        solubility_parameters: gas.SurfaceSolubilityParameters,
        exchange_parameters: precipitation.GasExchangeParameters,
        parameters: RuntimeParameters,
        options: solver.Options,
    ) !solver.Result {
        @memset(self.atmospheric_flux_g_per_h, 0);
        if (gas_state.cell_count != self.cell_count or geometry.cell_count != self.cell_count or litter_water_m3.len != self.cell_count or litter_ice_m3.len != self.cell_count or cell_area_m2.len != self.cell_count or atmospheric_conductance_m3_per_step.len != self.cell_count) return error.SurfaceLitterGasDimensionMismatch;
        if (!std.math.isFinite(parameters.reference_temperature_k) or parameters.reference_temperature_k <= 0 or !std.math.isFinite(parameters.temperature_exponent) or parameters.temperature_exponent < 0 or !std.math.isFinite(parameters.penman_tortuosity) or parameters.penman_tortuosity < 0 or !std.math.isFinite(parameters.minimum_air_fraction) or parameters.minimum_air_fraction < 0) return error.InvalidSurfaceLitterGasParameter;
        for (0..self.cell_count) |cell| {
            const volume = geometry.expanded_total_volume_m3[cell];
            const porosity = geometry.porosity_m3_per_m3[cell];
            const air_volume = geometry.air_volume_m3[cell];
            const area = cell_area_m2[cell];
            const thickness = if (area > 0) volume / area else 0;
            if (!std.math.isFinite(volume) or volume < 0 or !std.math.isFinite(porosity) or porosity < 0 or !std.math.isFinite(air_volume) or air_volume < 0 or air_volume > volume or !std.math.isFinite(area) or area <= 0 or !std.math.isFinite(litter_water_m3[cell]) or litter_water_m3[cell] < 0 or !std.math.isFinite(litter_ice_m3[cell]) or litter_ice_m3[cell] < 0) return error.InvalidSurfaceLitterGasState;
            gas_state.air_volume_m3[cell] = air_volume;
            // Freeze-out: ice excludes dissolved gases. When all pore space is
            // occupied by ice and liquid, air_volume = 0 and the gas-water
            // exchange rate is zero (litterGasExchange requires air_m3 > 0).
            // Without this transfer, dissolved mass accumulates in the shrinking
            // liquid volume, producing physically impossible concentrations.
            // The atmospheric boundary (atmosphericDiffusiveFluxG, air_volume=0
            // branch) then immediately expels the degassed mass to atmosphere.
            if (air_volume <= 0) {
                const start = cell * gas.species_count;
                const dissolved = gas_state.dissolved_mass_g[start..][0..gas.species_count];
                const gaseous = gas_state.gaseous_mass_g[start..][0..gas.species_count];
                for (dissolved, gaseous) |*d, *g| {
                    g.* += d.*;
                    d.* = 0;
                }
            }
            self.water_volume_m3[cell] = litter_water_m3[cell];
            const solubility = try gas.surfaceSolubilityWaterToAir(gas_state.temperature_k[cell], solubility_parameters);
            var whole_step_exchange = exchange_parameters;
            whole_step_exchange.iteration_fraction = 1;
            const exchange = try precipitation.litterGasExchange(geometry.pore_volume_m3[cell], litter_ice_m3[cell], litter_water_m3[cell], air_volume, geometry.water_retention_capacity_m3[cell], whole_step_exchange);
            var interior: [gas.species_count]f64 = undefined;
            const air_fraction = if (volume > 0) air_volume / volume else 0;
            const diffusion_geometry_m = if (thickness > 0 and porosity > 0 and air_fraction > parameters.minimum_air_fraction)
                air_fraction * parameters.penman_tortuosity * air_fraction / porosity * area / thickness
            else
                0;
            const temperature_factor = std.math.pow(f64, gas_state.temperature_k[cell] / parameters.reference_temperature_k, parameters.temperature_exponent);
            for (0..gas.species_count) |species| {
                const component = cell * gas.species_count + species;
                interior[species] = diffusion_geometry_m * parameters.free_air_diffusivity_m2_per_h[species] * temperature_factor;
                self.mass_solubility_ratio[component] = solubility[species];
                self.gas_water_exchange_rate_per_step[component] = exchange.air_water_rate_per_step;
            }
            self.boundaries[cell] = .{
                .cell_index = cell,
                .aerodynamic_conductance_m3_per_step = atmospheric_conductance_m3_per_step[cell],
                .interior_conductance_m3_per_step = interior,
                .atmospheric_concentration_g_per_m3 = atmospheric_concentration_g_per_m3,
            };
        }
        return solver.solve(self.allocator, gas_state, .{
            .faces = &.{},
            .face_conductance_m3_per_step = &.{},
            .atmospheric_boundaries = self.boundaries,
            .water_volume_m3 = self.water_volume_m3,
            .band_water_volume_m3 = self.band_water_volume_m3,
            .mass_solubility_ratio = self.mass_solubility_ratio,
            .gas_water_exchange_rate_per_step = self.gas_water_exchange_rate_per_step,
            .band_gas_water_exchange_rate_per_step = self.band_gas_water_exchange_rate_per_step,
            .bubbling_enabled = self.bubbling_enabled,
            .atmospheric_flux_g_by_component = self.atmospheric_flux_g_per_h,
        }, options);
    }
};

test "PARR couples atmosphere and litter diffusivity in one local solve" {
    var geometry = try litter_geometry.State.init(std.testing.allocator, 1);
    defer geometry.deinit();
    geometry.expanded_total_volume_m3[0] = 2;
    geometry.pore_volume_m3[0] = 1;
    geometry.air_volume_m3[0] = 0.5;
    geometry.porosity_m3_per_m3[0] = 0.5;
    geometry.water_retention_capacity_m3[0] = 0.5;
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.temperature_k[0] = 298.15;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var atmospheric = [_]f64{0} ** gas.species_count;
    atmospheric[0] = 1;
    const unit_solubility = gas.SurfaceSolubilityParameters{ .reference_water_to_air = [_]f64{1} ** gas.species_count, .log_intercept = [_]f64{0} ** gas.species_count, .temperature_coefficient_per_c = [_]f64{0} ** gas.species_count };
    _ = try state.advance(&gas_state, &geometry, &.{0.5}, &.{0}, &.{2}, &.{0.1}, atmospheric, unit_solubility, .{ .reference_time_h = 1, .wet_exponent = 12, .dry_exponent = 12, .transition_water_fraction = 0.5, .iteration_fraction = 0, .aqueous_tortuosity_coefficient = 0.7 }, .{}, .{ .max_iterations = 80 });
    try std.testing.expect(state.boundaries[0].interior_conductance_m3_per_step[0] > 0);
    try std.testing.expect(state.gas_water_exchange_rate_per_step[0] > 0);
    try std.testing.expect(gas_state.gaseous_mass_g[0] + gas_state.dissolved_mass_g[0] > 0);
    try std.testing.expectApproxEqAbs(gas_state.gaseous_mass_g[0] + gas_state.dissolved_mass_g[0], state.atmospheric_flux_g_per_h[0], 1e-12);
}
