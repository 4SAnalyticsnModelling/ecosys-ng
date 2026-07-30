const std = @import("std");
const grid_module = @import("grid.zig");
const organic_transport = @import("soil_organic_transport.zig");
const dissolved_gas_transport = @import("soil_dissolved_gas_transport.zig");
const gaseous_transport = @import("soil_gas_transport_step.zig");
const gas = @import("gas_transport.zig");

pub const State = struct {
    allocator: std.mem.Allocator,
    dissolved_organic_carbon_runoff_g: []f64,
    dissolved_inorganic_carbon_runoff_g: []f64,
    dissolved_organic_carbon_drainage_g: []f64,
    dissolved_inorganic_carbon_drainage_g: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroDailyCarbonExportCells;
        const doc_runoff = try zeroes(allocator, cell_count);
        errdefer allocator.free(doc_runoff);
        const dic_runoff = try zeroes(allocator, cell_count);
        errdefer allocator.free(dic_runoff);
        const doc_drainage = try zeroes(allocator, cell_count);
        errdefer allocator.free(doc_drainage);
        const dic_drainage = try zeroes(allocator, cell_count);
        return .{
            .allocator = allocator,
            .dissolved_organic_carbon_runoff_g = doc_runoff,
            .dissolved_inorganic_carbon_runoff_g = dic_runoff,
            .dissolved_organic_carbon_drainage_g = doc_drainage,
            .dissolved_inorganic_carbon_drainage_g = dic_drainage,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.dissolved_inorganic_carbon_drainage_g);
        self.allocator.free(self.dissolved_organic_carbon_drainage_g);
        self.allocator.free(self.dissolved_inorganic_carbon_runoff_g);
        self.allocator.free(self.dissolved_organic_carbon_runoff_g);
        self.* = undefined;
    }

    pub fn resetDaily(self: *State) void {
        @memset(self.dissolved_organic_carbon_runoff_g, 0);
        @memset(self.dissolved_inorganic_carbon_runoff_g, 0);
        @memset(self.dissolved_organic_carbon_drainage_g, 0);
        @memset(self.dissolved_inorganic_carbon_drainage_g, 0);
    }

    /// REDIST `UDOCD`: DOC plus acetate lost through all profile boundaries.
    /// The transport ledger is source-signed (positive enters), so only its
    /// negative component is an outward loss.
    pub fn accumulateOrganicDrainageHour(self: *State, model_grid: *const grid_module.GridState, transport: *const organic_transport.State) !void {
        if (self.dissolved_organic_carbon_drainage_g.len != model_grid.cell_count or transport.layer_count != model_grid.layer_count) return error.DailyCarbonExportDimensionMismatch;
        for (0..model_grid.cell_count) |cell| {
            var outward_g_c: f64 = 0;
            for (0..model_grid.active_soil_layer_count[cell]) |local_layer| {
                const layer = try model_grid.layerIndex(cell, local_layer);
                for (0..@import("soil_organic_initialization.zig").substrate_count) |substrate| {
                    const base = layer * organic_transport.component_count + substrate * organic_transport.components_per_substrate;
                    outward_g_c += @max(0, -transport.boundary_net_flux_g[base + @intFromEnum(organic_transport.Component.dissolved_organic_carbon)]);
                    outward_g_c += @max(0, -transport.boundary_net_flux_g[base + @intFromEnum(organic_transport.Component.dissolved_acetate_carbon)]);
                }
            }
            if (!std.math.isFinite(outward_g_c)) return error.NonFiniteDailyCarbonExport;
            const after = self.dissolved_organic_carbon_drainage_g[cell] + outward_g_c;
            if (!std.math.isFinite(after)) return error.NonFiniteDailyCarbonExport;
            self.dissolved_organic_carbon_drainage_g[cell] = after;
        }
    }

    pub fn accumulateOrganicRunoffHour(self: *State, outward_g_c_by_cell: []const f64) !void {
        if (outward_g_c_by_cell.len != self.dissolved_organic_carbon_runoff_g.len) return error.DailyCarbonExportDimensionMismatch;
        for (self.dissolved_organic_carbon_runoff_g, outward_g_c_by_cell) |*daily, outward| {
            if (!std.math.isFinite(outward) or outward < 0 or !std.math.isFinite(daily.* + outward)) return error.NonFiniteDailyCarbonExport;
            daily.* += outward;
        }
    }

    pub fn accumulateInorganicRunoffHour(self: *State, outward_g_c_by_cell: []const f64) !void {
        if (outward_g_c_by_cell.len != self.dissolved_inorganic_carbon_runoff_g.len) return error.DailyCarbonExportDimensionMismatch;
        for (self.dissolved_inorganic_carbon_runoff_g, outward_g_c_by_cell) |*daily, outward| {
            if (!std.math.isFinite(outward) or outward < 0 or !std.math.isFinite(daily.* + outward)) return error.NonFiniteDailyCarbonExport;
            daily.* += outward;
        }
    }

    /// REDIST aqueous `XCOFLS/XCOFHS/XCHFLS/XCHFHS` contribution to UDICD.
    pub fn accumulateDissolvedInorganicDrainageHour(self: *State, model_grid: *const grid_module.GridState, transport: *const dissolved_gas_transport.State) !void {
        if (transport.layer_count != model_grid.layer_count or self.dissolved_inorganic_carbon_drainage_g.len != model_grid.cell_count) return error.DailyCarbonExportDimensionMismatch;
        for (0..model_grid.cell_count) |cell| {
            var outward_g_c: f64 = 0;
            for (0..model_grid.active_soil_layer_count[cell]) |local_layer| {
                const layer = try model_grid.layerIndex(cell, local_layer);
                const base = layer * gas.species_count;
                outward_g_c += @max(0, -transport.boundary_net_flux_g[base + @intFromEnum(gas.Species.carbon_dioxide)]);
                outward_g_c += @max(0, -transport.boundary_net_flux_g[base + @intFromEnum(gas.Species.methane)]);
            }
            const after = self.dissolved_inorganic_carbon_drainage_g[cell] + outward_g_c;
            if (!std.math.isFinite(after)) return error.NonFiniteDailyCarbonExport;
            self.dissolved_inorganic_carbon_drainage_g[cell] = after;
        }
    }

    /// REDIST gaseous `XCOFLG/XCHFLG` contribution to UDICD. Atmospheric
    /// top-boundary exchange is intentionally excluded; only the separately
    /// classified lateral and profile-bottom ledger is consumed here.
    pub fn accumulateGaseousInorganicDrainageHour(self: *State, model_grid: *const grid_module.GridState, transport: *const gaseous_transport.State) !void {
        if (transport.subsurface_flux_g_per_h.len != try std.math.mul(usize, model_grid.layer_count, gas.species_count) or self.dissolved_inorganic_carbon_drainage_g.len != model_grid.cell_count) return error.DailyCarbonExportDimensionMismatch;
        for (0..model_grid.cell_count) |cell| {
            var outward_g_c: f64 = 0;
            for (0..model_grid.active_soil_layer_count[cell]) |local_layer| {
                const layer = try model_grid.layerIndex(cell, local_layer);
                const base = layer * gas.species_count;
                outward_g_c += @max(0, -transport.subsurface_flux_g_per_h[base + @intFromEnum(gas.Species.carbon_dioxide)]);
                outward_g_c += @max(0, -transport.subsurface_flux_g_per_h[base + @intFromEnum(gas.Species.methane)]);
            }
            const after = self.dissolved_inorganic_carbon_drainage_g[cell] + outward_g_c;
            if (!std.math.isFinite(after)) return error.NonFiniteDailyCarbonExport;
            self.dissolved_inorganic_carbon_drainage_g[cell] = after;
        }
    }
};

fn zeroes(allocator: std.mem.Allocator, count: usize) ![]f64 {
    const values = try allocator.alloc(f64, count);
    @memset(values, 0);
    return values;
}

test "UDOCD sums DOC and acetate outward losses over runtime layers and complexes" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var model_grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer model_grid.deinit();
    var transport = try organic_transport.State.init(std.testing.allocator, model_grid.layer_count);
    defer transport.deinit();
    const first = 0;
    const second = organic_transport.component_count;
    transport.boundary_net_flux_g[first + @intFromEnum(organic_transport.Component.dissolved_organic_carbon)] = -2;
    transport.boundary_net_flux_g[first + @intFromEnum(organic_transport.Component.dissolved_acetate_carbon)] = -3;
    transport.boundary_net_flux_g[second + organic_transport.components_per_substrate + @intFromEnum(organic_transport.Component.dissolved_organic_carbon)] = -5;
    transport.boundary_net_flux_g[second + @intFromEnum(organic_transport.Component.dissolved_organic_carbon)] = 17;
    var daily = try State.init(std.testing.allocator, model_grid.cell_count);
    defer daily.deinit();
    try daily.accumulateOrganicDrainageHour(&model_grid, &transport);
    try std.testing.expectApproxEqAbs(@as(f64, 10), daily.dissolved_organic_carbon_drainage_g[0], 1e-15);
}

test "UDICD includes only outward subsurface gaseous CO2 and CH4" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var model_grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer model_grid.deinit();
    var transport = try gaseous_transport.State.init(std.testing.allocator, model_grid.layer_count);
    defer transport.deinit();
    transport.subsurface_flux_g_per_h[@intFromEnum(gas.Species.carbon_dioxide)] = -2;
    transport.subsurface_flux_g_per_h[@intFromEnum(gas.Species.methane)] = 7;
    transport.subsurface_flux_g_per_h[gas.species_count + @intFromEnum(gas.Species.methane)] = -3;
    // Surface exchange must never enter profile drainage.
    transport.atmospheric_flux_g_per_h[@intFromEnum(gas.Species.carbon_dioxide)] = -100;
    var daily = try State.init(std.testing.allocator, model_grid.cell_count);
    defer daily.deinit();
    try daily.accumulateGaseousInorganicDrainageHour(&model_grid, &transport);
    try std.testing.expectApproxEqAbs(@as(f64, 5), daily.dissolved_inorganic_carbon_drainage_g[0], 1e-15);
}
