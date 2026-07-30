const std = @import("std");
const gas = @import("gas_transport.zig");
const RootSystem = @import("plant_root_system.zig").State;
const root_disturbance = @import("plant_root_disturbance.zig");
const output_binding = @import("soil_hourly_output_binding.zig");

/// Heap-owned DAY accumulators for accepted soil/litter atmospheric exchange
/// plus source-signed root withdrawal. Values retain tracked-element units.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    soil_litter_boundary_mass_g_by_cell_and_species: []f64,
    tracked_element_mass_g_by_cell_and_species: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroDailyGasFluxCells;
        const count = try std.math.mul(usize, cell_count, gas.species_count);
        const boundary = try allocator.alloc(f64, count);
        errdefer allocator.free(boundary);
        const combined = try allocator.alloc(f64, count);
        @memset(boundary, 0);
        @memset(combined, 0);
        return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .soil_litter_boundary_mass_g_by_cell_and_species = boundary,
            .tracked_element_mass_g_by_cell_and_species = combined,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.soil_litter_boundary_mass_g_by_cell_and_species);
        self.allocator.free(self.tracked_element_mass_g_by_cell_and_species);
        self.* = undefined;
    }

    pub fn reset(self: *State) void {
        @memset(self.soil_litter_boundary_mass_g_by_cell_and_species, 0);
        @memset(self.tracked_element_mass_g_by_cell_and_species, 0);
    }

    pub fn accumulateHour(
        self: *State,
        roots: ?*const RootSystem,
        plant_species_count: usize,
        soil_layer_count: usize,
        soil_atmospheric_flux_g_per_h: []const f64,
        litter_atmospheric_flux_g_per_h: []const f64,
    ) !void {
        if (soil_layer_count == 0 or soil_atmospheric_flux_g_per_h.len != self.cell_count * soil_layer_count * gas.species_count or litter_atmospheric_flux_g_per_h.len != self.cell_count * gas.species_count) return error.DailyGasFluxDimensionMismatch;
        // Validate the complete grid before changing any daily accumulator.
        for (0..self.cell_count) |cell| {
            const withdrawal = if (roots) |root_state|
                try root_disturbance.rootGasWithdrawalForCell(root_state, cell, plant_species_count)
            else
                root_disturbance.CellRootGasWithdrawal{};
            for (0..gas.species_count) |species_index| {
                const species: gas.Species = @enumFromInt(species_index);
                const root_atmosphere_exchange = if (roots) |root_state| try rootAtmosphereExchangeGPerH(root_state, cell, plant_species_count, species) else 0;
                const boundary_increment = try boundaryHourIncrement(cell, soil_layer_count, species, soil_atmospheric_flux_g_per_h, litter_atmospheric_flux_g_per_h);
                const increment = try combinedHourIncrement(species, withdrawal, root_atmosphere_exchange, boundary_increment);
                const index = cell * gas.species_count + species_index;
                const next_boundary = self.soil_litter_boundary_mass_g_by_cell_and_species[index] + boundary_increment;
                const next = self.tracked_element_mass_g_by_cell_and_species[index] + increment;
                if (!std.math.isFinite(next_boundary) or !std.math.isFinite(next)) return error.NonFiniteDailyGasFlux;
            }
        }
        for (0..self.cell_count) |cell| {
            const withdrawal = if (roots) |root_state|
                root_disturbance.rootGasWithdrawalForCell(root_state, cell, plant_species_count) catch unreachable
            else
                root_disturbance.CellRootGasWithdrawal{};
            for (0..gas.species_count) |species_index| {
                const species: gas.Species = @enumFromInt(species_index);
                const root_atmosphere_exchange = if (roots) |root_state| rootAtmosphereExchangeGPerH(root_state, cell, plant_species_count, species) catch unreachable else 0;
                const boundary_increment = boundaryHourIncrement(cell, soil_layer_count, species, soil_atmospheric_flux_g_per_h, litter_atmospheric_flux_g_per_h) catch unreachable;
                const increment = combinedHourIncrement(species, withdrawal, root_atmosphere_exchange, boundary_increment) catch unreachable;
                const index = cell * gas.species_count + species_index;
                self.soil_litter_boundary_mass_g_by_cell_and_species[index] += boundary_increment;
                self.tracked_element_mass_g_by_cell_and_species[index] += increment;
            }
        }
    }

    pub fn get(self: *const State, cell: usize, species: gas.Species) !f64 {
        if (cell >= self.cell_count) return error.DailyGasFluxCellOutOfBounds;
        return self.tracked_element_mass_g_by_cell_and_species[cell * gas.species_count + @intFromEnum(species)];
    }

    /// Exact DAY `UCO2G`, `UCH4G`, `UOXYG`, `UH2GG`, ... owner: cumulative
    /// snowpack+litter+soil exchange only. Root withdrawal remains available
    /// through `get` for consumers that explicitly require the combined
    /// ecosystem exchange.
    pub fn getSoilLitterBoundary(self: *const State, cell: usize, species: gas.Species) !f64 {
        if (cell >= self.cell_count) return error.DailyGasFluxCellOutOfBounds;
        return self.soil_litter_boundary_mass_g_by_cell_and_species[cell * gas.species_count + @intFromEnum(species)];
    }
};

fn boundaryHourIncrement(cell: usize, soil_layer_count: usize, species: gas.Species, soil_flux: []const f64, litter_flux: []const f64) !f64 {
    return output_binding.gasBoundaryExchangeG(soil_flux, litter_flux, cell * soil_layer_count, soil_layer_count, cell, species);
}

fn combinedHourIncrement(species: gas.Species, withdrawal: root_disturbance.CellRootGasWithdrawal, root_atmosphere_exchange_g_per_h: f64, boundary: f64) !f64 {
    const root_flux: f64 = switch (species) {
        .carbon_dioxide => withdrawal.carbon_dioxide_g_c_per_h,
        .oxygen => withdrawal.oxygen_g_o_per_h,
        .methane => withdrawal.methane_g_c_per_h,
        .nitrous_oxide => withdrawal.nitrous_oxide_g_n_per_h,
        .ammonia => withdrawal.ammonia_g_n_per_h,
        .hydrogen => withdrawal.hydrogen_g_h_per_h,
        .nitrogen => 0,
    };
    // Accepted root flux is positive atmosphere -> root, whereas DAY gas
    // output is positive ecosystem -> atmosphere.
    const result = boundary + root_flux - root_atmosphere_exchange_g_per_h;
    if (!std.math.isFinite(result)) return error.NonFiniteDailyGasFlux;
    return result;
}

fn rootAtmosphereExchangeGPerH(roots: *const RootSystem, cell: usize, plant_species_count: usize, species: gas.Species) !f64 {
    if (plant_species_count == 0 or roots.plant_count % plant_species_count != 0 or cell >= roots.plant_count / plant_species_count)
        return error.DailyGasFluxDimensionMismatch;
    const transaction_gas: ?usize = switch (species) {
        .carbon_dioxide => 0,
        .methane => 1,
        .nitrous_oxide => 2,
        .ammonia => 3,
        .hydrogen => 4,
        .oxygen => 5,
        .nitrogen => null,
    };
    const gas_slot = transaction_gas orelse return 0;
    var total: f64 = 0;
    const first_plant = cell * plant_species_count;
    for (first_plant..first_plant + plant_species_count) |plant|
        for (0..@import("plant_root_system.zig").biological_domain_count) |domain|
            for (0..roots.soil_layer_count) |layer| {
                const root = try roots.layerIndex(plant, domain, layer);
                total += roots.atmosphere_to_root_gas_exchange_g_per_h[root * @import("plant_root_system.zig").transported_root_gas_count + gas_slot];
            };
    if (!std.math.isFinite(total)) return error.NonFiniteDailyGasFlux;
    return total;
}

test "DAY gas accumulator combines boundaries and runtime root species" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var roots = try RootSystem.init(std.testing.allocator, 7, 2, 1);
    defer roots.deinit();
    roots.withdrawal_carbon_dioxide_loss_g_c_per_h[6] = -2;
    roots.withdrawal_hydrogen_loss_g_h_per_h[0] = -0.5;
    const carbon_root = try roots.layerIndex(6, 0, 0);
    const hydrogen_root = try roots.layerIndex(0, 0, 0);
    roots.atmosphere_to_root_gas_exchange_g_per_h[carbon_root * 6] = 0.75;
    roots.atmosphere_to_root_gas_exchange_g_per_h[hydrogen_root * 6 + 4] = 0.25;
    var soil = [_]f64{0} ** (2 * gas.species_count);
    var litter = [_]f64{0} ** gas.species_count;
    soil[@intFromEnum(gas.Species.carbon_dioxide)] = 3;
    litter[@intFromEnum(gas.Species.carbon_dioxide)] = 1;
    litter[@intFromEnum(gas.Species.hydrogen)] = 1.5;
    try state.accumulateHour(&roots, 7, 2, &soil, &litter);
    try state.accumulateHour(&roots, 7, 2, &soil, &litter);
    try std.testing.expectEqual(@as(f64, 8), try state.getSoilLitterBoundary(0, .carbon_dioxide));
    try std.testing.expectEqual(@as(f64, 3), try state.getSoilLitterBoundary(0, .hydrogen));
    try std.testing.expectEqual(@as(f64, 2.5), try state.get(0, .carbon_dioxide));
    try std.testing.expectEqual(@as(f64, 1.5), try state.get(0, .hydrogen));
    state.reset();
    try std.testing.expectEqual(@as(f64, 0), try state.getSoilLitterBoundary(0, .carbon_dioxide));
    try std.testing.expectEqual(@as(f64, 0), try state.get(0, .carbon_dioxide));
}

test "failed DAY gas accumulation leaves every cell unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.tracked_element_mass_g_by_cell_and_species[0] = 4;
    var soil = [_]f64{0} ** (2 * gas.species_count);
    var litter = [_]f64{0} ** (2 * gas.species_count);
    litter[gas.species_count] = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteGasOutputBoundaryExchange, state.accumulateHour(null, 1, 1, &soil, &litter));
    try std.testing.expectEqual(@as(f64, 4), state.tracked_element_mass_g_by_cell_and_species[0]);
}
