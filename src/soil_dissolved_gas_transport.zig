const std = @import("std");
const gas = @import("gas_transport.zig");
const hydrology = @import("transport_hydrology.zig");
const transport = @import("aqueous_extensive_transport.zig");

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    boundary_net_flux_g: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.ZeroSoilDissolvedGasTransportLayers;
        const ledger = try allocator.alloc(f64, layer_count * gas.species_count);
        @memset(ledger, 0);
        return .{ .allocator = allocator, .layer_count = layer_count, .boundary_net_flux_g = ledger };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.boundary_net_flux_g);
        self.* = undefined;
    }

    pub fn validateFinite(self: *const State) !void {
        for (self.boundary_net_flux_g, 0..) |value, index| {
            if (!std.math.isFinite(value)) {
                std.log.err(
                    "non-finite dissolved-gas boundary ledger: index={d} value={e}",
                    .{ index, value },
                );
                return error.NonFiniteDissolvedGasBoundaryLedger;
            }
        }
    }
};

pub const Inputs = struct {
    faces: *const hydrology.SoilFaces,
    micropore_conductance_m3_per_step: []const f64,
    macropore_conductance_m3_per_step: []const f64,
    micropore_water_m3: []const f64,
    macropore_water_m3: []const f64,
    layer_bulk_volume_m3: []const f64,
    micropore_external_water_flux_m3_per_step: []const f64,
    macropore_external_water_flux_m3_per_step: []const f64,
    recharge_concentration_g_per_m3: []const f64,
};

/// TRNSFR aqueous transport for CO2/CH4/O2/N2/N2O/NH3/H2 in both pore
/// domains. The accepted source-signed ledger is consumed directly by REDIST.
pub fn advance(allocator: std.mem.Allocator, state: *State, gas_state: *gas.State, inputs: Inputs, options: transport.Options) !transport.Result {
    if (gas_state.cell_count != state.layer_count) return error.SoilDissolvedGasTransportDimensionMismatch;
    return transport.advance(
        allocator,
        gas_state.dissolved_mass_g,
        gas_state.macropore_dissolved_mass_g,
        state.boundary_net_flux_g,
        .{
            .species_count = gas.species_count,
            .faces = inputs.faces.micropore_faces,
            .micropore_conductance_m3_per_step = inputs.micropore_conductance_m3_per_step,
            .macropore_conductance_m3_per_step = inputs.macropore_conductance_m3_per_step,
            .micropore_water_m3 = inputs.micropore_water_m3,
            .macropore_water_m3 = inputs.macropore_water_m3,
            .layer_bulk_volume_m3 = inputs.layer_bulk_volume_m3,
            .micropore_external_water_flux_m3_per_step = inputs.micropore_external_water_flux_m3_per_step,
            .macropore_external_water_flux_m3_per_step = inputs.macropore_external_water_flux_m3_per_step,
            .recharge_concentration_per_m3 = inputs.recharge_concentration_g_per_m3,
        },
        options,
    );
}

test "CO2 and CH4 boundary losses retain tracked-carbon grams and source sign" {
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.carbon_dioxide)] = 8;
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.methane)] = 4;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var model_grid = try @import("grid.zig").GridState.init(std.testing.allocator, config);
    defer model_grid.deinit();
    @memset(model_grid.active_soil_layer_count, 1);
    var hydro = try hydrology.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer hydro.deinit();
    var faces = try hydrology.buildSoilFaces(std.testing.allocator, &hydro, &model_grid);
    defer faces.deinit();
    const zero_recharge = [_]f64{0} ** gas.species_count;
    _ = try advance(std.testing.allocator, &state, &gas_state, .{
        .faces = &faces,
        .micropore_conductance_m3_per_step = &.{},
        .macropore_conductance_m3_per_step = &.{},
        .micropore_water_m3 = &.{1},
        .macropore_water_m3 = &.{0},
        .layer_bulk_volume_m3 = &.{1},
        .micropore_external_water_flux_m3_per_step = &.{0.25},
        .macropore_external_water_flux_m3_per_step = &.{0},
        .recharge_concentration_g_per_m3 = &zero_recharge,
    }, .{ .absolute_tolerance = 1e-12, .relative_tolerance = 1e-8, .picard_relaxation = 0.5, .max_iterations = 20, .pore_exchange_fraction = 0 });
    try std.testing.expectEqual(@as(f64, -2), state.boundary_net_flux_g[@intFromEnum(gas.Species.carbon_dioxide)]);
    try std.testing.expectEqual(@as(f64, -1), state.boundary_net_flux_g[@intFromEnum(gas.Species.methane)]);
}
