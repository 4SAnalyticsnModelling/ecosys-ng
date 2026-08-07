//! Tests for `coupled_gas_solver.zig`.
//!
//! Extracted verbatim so the module beside it contains only the model
//! code. Tests that use private declarations of that module stay there,
//! since a sibling file can only reach `pub` declarations.

const atmosphere = @import("../soil/gas/atmosphere_exchange.zig");
const gas = @import("../soil/gas/transport.zig");
const numerics = @import("../core/numerics.zig");
const std = @import("std");
const coupled_gas_solver = @import("../soil/gas/coupled_gas_solver.zig");
test "coupled coupled_gas_solver.solve converges before NPH times NPG and conserves closed system" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.air_volume_m3[1] = 1;
    state.temperature_k[0] = 300;
    state.temperature_k[1] = 300;
    state.gaseous_mass_g[0] = 2;
    state.dissolved_mass_g[0] = 1;
    const n = 2 * gas.species_count;
    const conductance = [_]f64{0.01} ** gas.species_count;
    const water = [_]f64{ 1, 1 };
    const no_band_water = [_]f64{ 0, 0 };
    const solubility = [_]f64{1} ** n;
    const exchange = [_]f64{0.1} ** n;
    const no_exchange = [_]f64{0} ** n;
    const no_bubbling = [_]bool{ false, false };
    var accepted_face_flux_g = [_]f64{0} ** gas.species_count;
    const result = try coupled_gas_solver.solve(std.testing.allocator, &state, .{ .faces = &[_]gas.Face{.{ .first_cell = 0, .second_cell = 1 }}, .face_conductance_m3_per_step = &conductance, .atmospheric_boundaries = &.{}, .water_volume_m3 = &water, .band_water_volume_m3 = &no_band_water, .mass_solubility_ratio = &solubility, .gas_water_exchange_rate_per_step = &exchange, .band_gas_water_exchange_rate_per_step = &no_exchange, .bubbling_enabled = &no_bubbling, .face_flux_g_by_component = &accepted_face_flux_g }, .{ .max_iterations = 80 });
    try std.testing.expect(result.iterations < 80);
    try std.testing.expect(accepted_face_flux_g[0] > 0);
    var total: f64 = 0;
    for (state.gaseous_mass_g, state.dissolved_mass_g, state.band_dissolved_mass_g) |gaseous, dissolved, band| total += gaseous + dissolved + band;
    try std.testing.expectApproxEqAbs(@as(f64, 3), total, 1e-10);
}

test "dimensionless dense Newton resolves gram gas beside hundred-megagram dissolved pools" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.temperature_k[0] = 300;
    const oxygen = @intFromEnum(gas.Species.oxygen);
    const nitrogen = @intFromEnum(gas.Species.nitrogen);
    state.gaseous_mass_g[oxygen] = 2;
    state.dissolved_mass_g[oxygen] = 1.0e8;
    state.gaseous_mass_g[nitrogen] = 1.8;
    state.dissolved_mass_g[nitrogen] = 1.36e8;
    const initial_oxygen_g =
        state.gaseous_mass_g[oxygen] + state.dissolved_mass_g[oxygen];
    const initial_nitrogen_g =
        state.gaseous_mass_g[nitrogen] + state.dissolved_mass_g[nitrogen];
    const water = [_]f64{1.0e8};
    const no_band_water = [_]f64{0};
    const solubility = [_]f64{1} ** gas.species_count;
    const exchange = [_]f64{1} ** gas.species_count;
    const no_exchange = [_]f64{0} ** gas.species_count;
    const no_bubbling = [_]bool{false};
    const result = try coupled_gas_solver.solve(
        std.testing.allocator,
        &state,
        .{
            .faces = &.{},
            .face_conductance_m3_per_step = &.{},
            .atmospheric_boundaries = &.{},
            .water_volume_m3 = &water,
            .band_water_volume_m3 = &no_band_water,
            .mass_solubility_ratio = &solubility,
            .gas_water_exchange_rate_per_step = &exchange,
            .band_gas_water_exchange_rate_per_step = &no_exchange,
            .bubbling_enabled = &no_bubbling,
        },
        .{
            .absolute_tolerance_g = 1.0e-11,
            .relative_tolerance = 1.0e-8,
            .max_iterations = 80,
        },
    );
    try std.testing.expect(result.iterations < 80);
    try std.testing.expectApproxEqAbs(
        initial_oxygen_g,
        state.gaseous_mass_g[oxygen] + state.dissolved_mass_g[oxygen],
        1.0e-7,
    );
    try std.testing.expectApproxEqAbs(
        initial_nitrogen_g,
        state.gaseous_mass_g[nitrogen] + state.dissolved_mass_g[nitrogen],
        1.0e-7,
    );
}

test "scale-separated face converges two pressure-coupled donor species conservatively" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    @memset(state.air_volume_m3, 1);
    @memset(state.temperature_k, 300);
    const carbon_dioxide = @intFromEnum(gas.Species.carbon_dioxide);
    const methane = @intFromEnum(gas.Species.methane);
    state.gaseous_mass_g[carbon_dioxide] = 1;
    state.gaseous_mass_g[methane] = 2;
    state.gaseous_mass_g[gas.species_count + carbon_dioxide] = 1.0e8;
    state.gaseous_mass_g[gas.species_count + methane] = 2.0e8;
    const initial_carbon_dioxide =
        state.gaseous_mass_g[carbon_dioxide] +
        state.gaseous_mass_g[gas.species_count + carbon_dioxide];
    const initial_methane =
        state.gaseous_mass_g[methane] +
        state.gaseous_mass_g[gas.species_count + methane];
    var conductance = [_]f64{0} ** gas.species_count;
    conductance[carbon_dioxide] = 1.0e-9;
    conductance[methane] = 1.0e-9;
    const n = 2 * gas.species_count;
    const zero_water = [_]f64{ 0, 0 };
    const solubility = [_]f64{1} ** n;
    const no_exchange = [_]f64{0} ** n;
    const no_bubbling = [_]bool{ false, false };
    const result = try coupled_gas_solver.solve(
        std.testing.allocator,
        &state,
        .{
            .faces = &.{.{
                .first_cell = 0,
                .second_cell = 1,
            }},
            .face_conductance_m3_per_step = &conductance,
            .atmospheric_boundaries = &.{},
            .water_volume_m3 = &zero_water,
            .band_water_volume_m3 = &zero_water,
            .mass_solubility_ratio = &solubility,
            .gas_water_exchange_rate_per_step = &no_exchange,
            .band_gas_water_exchange_rate_per_step = &no_exchange,
            .bubbling_enabled = &no_bubbling,
        },
        .{ .max_iterations = 80 },
    );
    try std.testing.expect(result.iterations < 80);
    try std.testing.expectApproxEqAbs(
        initial_carbon_dioxide,
        state.gaseous_mass_g[carbon_dioxide] +
            state.gaseous_mass_g[gas.species_count + carbon_dioxide],
        1.0e-7,
    );
    try std.testing.expectApproxEqAbs(
        initial_methane,
        state.gaseous_mass_g[methane] +
            state.gaseous_mass_g[gas.species_count + methane],
        1.0e-7,
    );
}

test "scale-separated face closes a sub-inventory oxygen residual conservatively" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    @memset(state.air_volume_m3, 1);
    @memset(state.temperature_k, 300);
    const oxygen = @intFromEnum(gas.Species.oxygen);
    state.gaseous_mass_g[oxygen] = 2.9218919364755397e-1;
    state.gaseous_mass_g[gas.species_count + oxygen] =
        4.596233521006912e7;
    const initial_oxygen_g =
        state.gaseous_mass_g[oxygen] +
        state.gaseous_mass_g[gas.species_count + oxygen];
    var conductance = [_]f64{0} ** gas.species_count;
    conductance[oxygen] = 5.165116406442861e-12;
    const n = 2 * gas.species_count;
    const zero_water = [_]f64{ 0, 0 };
    const solubility = [_]f64{1} ** n;
    const no_exchange = [_]f64{0} ** n;
    const no_bubbling = [_]bool{ false, false };
    const result = try coupled_gas_solver.solve(
        std.testing.allocator,
        &state,
        .{
            .faces = &.{.{
                .first_cell = 0,
                .second_cell = 1,
            }},
            .face_conductance_m3_per_step = &conductance,
            .atmospheric_boundaries = &.{},
            .water_volume_m3 = &zero_water,
            .band_water_volume_m3 = &zero_water,
            .mass_solubility_ratio = &solubility,
            .gas_water_exchange_rate_per_step = &no_exchange,
            .band_gas_water_exchange_rate_per_step = &no_exchange,
            .bubbling_enabled = &no_bubbling,
        },
        .{
            .absolute_tolerance_g = 1.0e-12,
            .relative_tolerance = 1.0e-8,
            .max_iterations = 80,
        },
    );
    try std.testing.expect(result.iterations < 80);
    try std.testing.expectApproxEqAbs(
        initial_oxygen_g,
        state.gaseous_mass_g[oxygen] +
            state.gaseous_mass_g[gas.species_count + oxygen],
        1.0e-8,
    );
}

test "simultaneous atmospheric sources release every zero-inventory bound" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    @memset(state.air_volume_m3, 1);
    @memset(state.temperature_k, 300);
    const oxygen = @intFromEnum(gas.Species.oxygen);
    var external = [_]f64{0} ** gas.species_count;
    external[oxygen] = 0.01;
    const conductance = [_]f64{0.1} ** gas.species_count;
    const boundaries = [_]atmosphere.Boundary{
        .{ .cell_index = 0, .aerodynamic_conductance_m3_per_step = 0.1, .interior_conductance_m3_per_step = conductance, .atmospheric_concentration_g_per_m3 = external },
        .{ .cell_index = 1, .aerodynamic_conductance_m3_per_step = 0.1, .interior_conductance_m3_per_step = conductance, .atmospheric_concentration_g_per_m3 = external },
    };
    const n = 2 * gas.species_count;
    const water = [_]f64{ 0, 0 };
    const solubility = [_]f64{1} ** n;
    const no_exchange = [_]f64{0} ** n;
    const no_bubbling = [_]bool{ false, false };
    const result = try coupled_gas_solver.solve(std.testing.allocator, &state, .{
        .faces = &.{},
        .face_conductance_m3_per_step = &.{},
        .atmospheric_boundaries = &boundaries,
        .water_volume_m3 = &water,
        .band_water_volume_m3 = &water,
        .mass_solubility_ratio = &solubility,
        .gas_water_exchange_rate_per_step = &no_exchange,
        .band_gas_water_exchange_rate_per_step = &no_exchange,
        .bubbling_enabled = &no_bubbling,
    }, .{ .max_iterations = 80 });
    try std.testing.expect(result.iterations < 80);
    try std.testing.expect(state.gaseous_mass_g[oxygen] > 0);
    try std.testing.expect(state.gaseous_mass_g[gas.species_count + oxygen] > 0);
}

test "surface and subsurface boundary fluxes remain independently classified" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.temperature_k[0] = 300;
    const carbon_dioxide = @intFromEnum(gas.Species.carbon_dioxide);
    // Ideal-capacity inventory suppresses the pressure correction so this
    // test isolates the two boundary ledgers.
    state.gaseous_mass_g[carbon_dioxide] = 1.2194e4 / 300.0 * 12.0;
    const initial_mass_g = state.gaseous_mass_g[carbon_dioxide];
    var high_external = [_]f64{0} ** gas.species_count;
    high_external[carbon_dioxide] = initial_mass_g + 1;
    const zero_external = [_]f64{0} ** gas.species_count;
    const conductance = [_]f64{0.1} ** gas.species_count;
    const surface = atmosphere.Boundary{
        .cell_index = 0,
        .aerodynamic_conductance_m3_per_step = 0.1,
        .interior_conductance_m3_per_step = conductance,
        .atmospheric_concentration_g_per_m3 = high_external,
    };
    const subsurface = atmosphere.Boundary{
        .cell_index = 0,
        .aerodynamic_conductance_m3_per_step = 0.1,
        .interior_conductance_m3_per_step = conductance,
        .atmospheric_concentration_g_per_m3 = zero_external,
    };
    const n = gas.species_count;
    const water = [_]f64{0};
    const solubility = [_]f64{1} ** n;
    const no_exchange = [_]f64{0} ** n;
    const no_bubbling = [_]bool{false};
    var surface_flux = [_]f64{0} ** n;
    var subsurface_flux = [_]f64{0} ** n;
    _ = try coupled_gas_solver.solve(std.testing.allocator, &state, .{
        .faces = &.{},
        .face_conductance_m3_per_step = &.{},
        .atmospheric_boundaries = &.{surface},
        .subsurface_boundaries = &.{subsurface},
        .water_volume_m3 = &water,
        .band_water_volume_m3 = &water,
        .mass_solubility_ratio = &solubility,
        .gas_water_exchange_rate_per_step = &no_exchange,
        .band_gas_water_exchange_rate_per_step = &no_exchange,
        .bubbling_enabled = &no_bubbling,
        .atmospheric_flux_g_by_component = &surface_flux,
        .subsurface_flux_g_by_component = &subsurface_flux,
    }, .{ .max_iterations = 80 });
    try std.testing.expect(surface_flux[carbon_dioxide] > 0);
    try std.testing.expect(subsurface_flux[carbon_dioxide] < 0);
    try std.testing.expectApproxEqAbs(
        initial_mass_g + surface_flux[carbon_dioxide] + subsurface_flux[carbon_dioxide],
        state.gaseous_mass_g[carbon_dioxide],
        1e-9,
    );
}

test "TRNSFRS bubbling transfers dissolved mass to the REDIST release layer" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    @memset(state.air_volume_m3, 1);
    @memset(state.temperature_k, 300);
    const carbon_dioxide = @intFromEnum(gas.Species.carbon_dioxide);
    const source = gas.species_count + carbon_dioxide;
    state.dissolved_mass_g[source] = 600;
    const initial_total_g = state.dissolved_mass_g[source];
    const n = 2 * gas.species_count;
    const water = [_]f64{ 1, 1 };
    const no_band_water = [_]f64{ 0, 0 };
    const solubility = [_]f64{1} ** n;
    const no_exchange = [_]f64{0} ** n;
    const bubbling = [_]bool{ false, true };
    const receivers = [_]?usize{ 0, 0 };
    _ = try coupled_gas_solver.solve(std.testing.allocator, &state, .{
        .faces = &.{},
        .face_conductance_m3_per_step = &.{},
        .atmospheric_boundaries = &.{},
        .water_volume_m3 = &water,
        .band_water_volume_m3 = &no_band_water,
        .mass_solubility_ratio = &solubility,
        .gas_water_exchange_rate_per_step = &no_exchange,
        .band_gas_water_exchange_rate_per_step = &no_exchange,
        .bubbling_enabled = &bubbling,
        .bubble_receiver_cell_by_cell = &receivers,
    }, .{ .max_iterations = 200 });
    try std.testing.expect(state.gaseous_mass_g[carbon_dioxide] > 0);
    try std.testing.expectApproxEqAbs(
        initial_total_g,
        state.gaseous_mass_g[carbon_dioxide] + state.dissolved_mass_g[source],
        1e-9,
    );
}

test "bubble without a gas release layer is published as a boundary loss" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.temperature_k[0] = 300;
    const carbon_dioxide = @intFromEnum(gas.Species.carbon_dioxide);
    state.dissolved_mass_g[carbon_dioxide] = 600;
    const initial_total_g = state.dissolved_mass_g[carbon_dioxide];
    const water = [_]f64{1};
    const no_band_water = [_]f64{0};
    const solubility = [_]f64{1} ** gas.species_count;
    const no_exchange = [_]f64{0} ** gas.species_count;
    const bubbling = [_]bool{true};
    const receivers = [_]?usize{null};
    var boundary_flux = [_]f64{0} ** gas.species_count;
    _ = try coupled_gas_solver.solve(std.testing.allocator, &state, .{
        .faces = &.{},
        .face_conductance_m3_per_step = &.{},
        .atmospheric_boundaries = &.{},
        .water_volume_m3 = &water,
        .band_water_volume_m3 = &no_band_water,
        .mass_solubility_ratio = &solubility,
        .gas_water_exchange_rate_per_step = &no_exchange,
        .band_gas_water_exchange_rate_per_step = &no_exchange,
        .bubbling_enabled = &bubbling,
        .bubble_receiver_cell_by_cell = &receivers,
        .subsurface_flux_g_by_component = &boundary_flux,
    }, .{ .max_iterations = 200 });
    try std.testing.expect(boundary_flux[carbon_dioxide] < 0);
    try std.testing.expectApproxEqAbs(
        initial_total_g + boundary_flux[carbon_dioxide],
        state.dissolved_mass_g[carbon_dioxide],
        1e-9,
    );
}

test "failed coupled coupled_gas_solver.solve rolls back all three phases" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.temperature_k[0] = 300;
    state.gaseous_mass_g[0] = 2;
    const n = gas.species_count;
    const water = [_]f64{1};
    const no_band_water = [_]f64{0};
    const solubility = [_]f64{1} ** n;
    const exchange = [_]f64{0.5} ** n;
    const no_exchange = [_]f64{0} ** n;
    const no_bubbling = [_]bool{false};
    try std.testing.expectError(error.CoupledGasSolverDidNotConverge, coupled_gas_solver.solve(std.testing.allocator, &state, .{ .faces = &.{}, .face_conductance_m3_per_step = &.{}, .atmospheric_boundaries = &.{}, .water_volume_m3 = &water, .band_water_volume_m3 = &no_band_water, .mass_solubility_ratio = &solubility, .gas_water_exchange_rate_per_step = &exchange, .band_gas_water_exchange_rate_per_step = &no_exchange, .bubbling_enabled = &no_bubbling }, .{ .absolute_tolerance_g = 1e-20, .relative_tolerance = 1e-20, .max_iterations = 1 }));
    try std.testing.expectEqual(@as(f64, 2), state.gaseous_mass_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.dissolved_mass_g[0]);
}
