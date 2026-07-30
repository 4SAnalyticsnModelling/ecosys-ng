const std = @import("std");
const settling = @import("pond_particulate_settling.zig");
const inventory = @import("surface_pond_inventory_transfer.zig");
const transition = @import("surface_pond_transition_step.zig");

pub const MineralSedimentOwners = struct {
    /// Suspended mineral sediment in surface water (Mg).
    surface_sediment_Mg: []f64,
    /// Extensive receiving soil carrier (Mg).
    surface_soil_mass_Mg: []f64,
    /// Hourly REDIST TSEDSK diagnostic (Mg), overwritten by this transaction.
    settled_sediment_Mg: []f64,
};

/// Production boundary for REDIST (`redist.f` lines 333-614) in the separated
/// surface representation. A domain-wide preflight preserves atomic failure.
pub fn apply(
    owners: inventory.Owners,
    mineral_sediment: MineralSedimentOwners,
    transitions: *const transition.State,
    timestep_h: f64,
) !void {
    const fraction = try settling.settlingFraction(.{ .timestep_h = timestep_h });
    if (transitions.cell_count != owners.surface_organic.layer_count or
        mineral_sediment.surface_sediment_Mg.len != transitions.cell_count or
        mineral_sediment.surface_soil_mass_Mg.len != transitions.cell_count or
        mineral_sediment.settled_sediment_Mg.len != transitions.cell_count)
        return error.SurfacePondSettlingDimensionMismatch;

    for (0..transitions.cell_count) |cell| {
        const suspended = mineral_sediment.surface_sediment_Mg[cell];
        const soil = mineral_sediment.surface_soil_mass_Mg[cell];
        const previous_diagnostic = mineral_sediment.settled_sediment_Mg[cell];
        if (!std.math.isFinite(suspended) or suspended < 0 or
            !std.math.isFinite(soil) or soil < 0 or
            !std.math.isFinite(previous_diagnostic) or previous_diagnostic < 0)
            return error.InvalidSurfacePondSedimentInventory;
        if (!transitions.active[cell]) continue;
        try inventory.validateParticulateFractionToSoil(owners, .{
            .cell = cell,
            .destination_soil_layer = transitions.destination_soil_layer[cell],
            .fraction = fraction,
        });
        const transfer_Mg = fraction * suspended;
        if (!std.math.isFinite(transfer_Mg) or
            !std.math.isFinite(soil + transfer_Mg))
            return error.InvalidSurfacePondSedimentInventory;
    }
    @memset(mineral_sediment.settled_sediment_Mg, 0);
    for (0..transitions.cell_count) |cell| {
        if (!transitions.active[cell]) continue;
        inventory.transferParticulateFractionToSoil(owners, .{
            .cell = cell,
            .destination_soil_layer = transitions.destination_soil_layer[cell],
            .fraction = fraction,
        }) catch unreachable;
        const transfer_Mg =
            fraction * mineral_sediment.surface_sediment_Mg[cell];
        mineral_sediment.surface_sediment_Mg[cell] -= transfer_Mg;
        mineral_sediment.surface_soil_mass_Mg[cell] += transfer_Mg;
        mineral_sediment.settled_sediment_Mg[cell] = transfer_Mg;
    }
}

test "production boundary settles represented particulates and preserves dissolved pools" {
    const organic = @import("soil_organic_initialization.zig");
    const gas = @import("gas_transport.zig");
    const surface_fertilizer = @import("surface_litter_fertilizer.zig");
    const soil_fertilizer = @import("fertilizer_nitrogen_inventory.zig");
    const mineral = @import("mineral_fertilizer_inventory.zig");
    var surface_organic = try organic.State.init(std.testing.allocator, 1);
    defer surface_organic.deinit();
    var soil_organic = try organic.State.init(std.testing.allocator, 1);
    defer soil_organic.deinit();
    var surface_gas = try gas.State.init(std.testing.allocator, 1);
    defer surface_gas.deinit();
    var soil_gas = try gas.State.init(std.testing.allocator, 1);
    defer soil_gas.deinit();
    var surface_n = try surface_fertilizer.State.init(std.testing.allocator, 1);
    defer surface_n.deinit();
    var soil_n = try soil_fertilizer.State.init(std.testing.allocator, 1, 1);
    defer soil_n.deinit();
    var mineral_state = try mineral.State.init(std.testing.allocator, 1, 1);
    defer mineral_state.deinit();
    var transitions = try transition.State.init(std.testing.allocator, 1);
    defer transitions.deinit();
    transitions.active[0] = true;
    surface_organic.microbial[0].carbon_g_c = 10;
    surface_organic.dissolved[0].carbon_g_c = 7;
    surface_n.cells[0].ammonium_mol_n = 20;
    mineral_state.surface[0].gypsum_mol = 30;
    const owners: inventory.Owners = .{ .surface_organic = &surface_organic, .soil_organic = &soil_organic, .surface_gas = &surface_gas, .soil_gas = &soil_gas, .surface_nitrogen_fertilizer = &surface_n, .soil_nitrogen_fertilizer = &soil_n, .mineral_fertilizer = &mineral_state };
    var suspended_Mg = [_]f64{2};
    var soil_mass_Mg = [_]f64{10};
    var settled_Mg = [_]f64{99};
    try apply(owners, .{
        .surface_sediment_Mg = &suspended_Mg,
        .surface_soil_mass_Mg = &soil_mass_Mg,
        .settled_sediment_Mg = &settled_Mg,
    }, &transitions, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), soil_organic.microbial[0].carbon_g_c, 1e-15);
    try std.testing.expectEqual(@as(f64, 7), surface_organic.dissolved[0].carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), soil_n.soil[0].broadcast_ammonium_mol_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.03), mineral_state.soil[0].gypsum_mol, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.998), suspended_Mg[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 10.002), soil_mass_Mg[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.002), settled_Mg[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 12), suspended_Mg[0] + soil_mass_Mg[0], 1e-15);
}

test "late invalid sediment owner leaves mixed-unit particulate owners unchanged" {
    const organic = @import("soil_organic_initialization.zig");
    const gas = @import("gas_transport.zig");
    const surface_fertilizer = @import("surface_litter_fertilizer.zig");
    const soil_fertilizer = @import("fertilizer_nitrogen_inventory.zig");
    const mineral = @import("mineral_fertilizer_inventory.zig");
    var surface_organic = try organic.State.init(std.testing.allocator, 1);
    defer surface_organic.deinit();
    var soil_organic = try organic.State.init(std.testing.allocator, 1);
    defer soil_organic.deinit();
    var surface_gas = try gas.State.init(std.testing.allocator, 1);
    defer surface_gas.deinit();
    var soil_gas = try gas.State.init(std.testing.allocator, 1);
    defer soil_gas.deinit();
    var surface_n = try surface_fertilizer.State.init(std.testing.allocator, 1);
    defer surface_n.deinit();
    var soil_n = try soil_fertilizer.State.init(std.testing.allocator, 1, 1);
    defer soil_n.deinit();
    var mineral_state = try mineral.State.init(std.testing.allocator, 1, 1);
    defer mineral_state.deinit();
    var transitions = try transition.State.init(std.testing.allocator, 1);
    defer transitions.deinit();
    transitions.active[0] = true;
    surface_organic.microbial[0].carbon_g_c = 5;
    var suspended_Mg = [_]f64{2};
    var soil_mass_Mg = [_]f64{std.math.inf(f64)};
    var settled_Mg = [_]f64{3};
    const owners: inventory.Owners = .{ .surface_organic = &surface_organic, .soil_organic = &soil_organic, .surface_gas = &surface_gas, .soil_gas = &soil_gas, .surface_nitrogen_fertilizer = &surface_n, .soil_nitrogen_fertilizer = &soil_n, .mineral_fertilizer = &mineral_state };
    try std.testing.expectError(error.InvalidSurfacePondSedimentInventory, apply(
        owners,
        .{ .surface_sediment_Mg = &suspended_Mg, .surface_soil_mass_Mg = &soil_mass_Mg, .settled_sediment_Mg = &settled_Mg },
        &transitions,
        1,
    ));
    try std.testing.expectEqual(@as(f64, 5), surface_organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), suspended_Mg[0]);
    try std.testing.expectEqual(@as(f64, 3), settled_Mg[0]);
}
