const std = @import("std");
const schedule = @import("fertilizer_schedule.zig");
const surface = @import("../surface/litter_fertilizer.zig");
const soil_dissolution = @import("../soil/nutrients/fertilizer_dissolution.zig");

/// Heap-owned undissolved nitrogen fertilizer for every runtime soil layer.
/// Surface broadcast material remains in `surface_litter_fertilizer.State`;
/// banded material is retained in soil even when its requested depth is zero.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    layer_capacity: usize,
    soil: []soil_dissolution.FertilizerState,
    initial_urease_inhibition_fraction: []f64,
    current_urease_inhibition_fraction: []f64,
    formulation: []u8,
    daily_nitrogen_input_g_n: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, layer_capacity: usize) !State {
        if (cell_count == 0 or layer_capacity == 0) return error.ZeroFertilizerInventoryExtent;
        const count = try std.math.mul(usize, cell_count, layer_capacity);
        const soil = try allocator.alloc(soil_dissolution.FertilizerState, count);
        errdefer allocator.free(soil);
        const initial = try allocator.alloc(f64, count);
        errdefer allocator.free(initial);
        const current = try allocator.alloc(f64, count);
        errdefer allocator.free(current);
        const formulation = try allocator.alloc(u8, count);
        errdefer allocator.free(formulation);
        const daily_nitrogen_input_g_n = try allocator.alloc(f64, cell_count);
        @memset(soil, zeroSoilInventory());
        @memset(initial, 0);
        @memset(current, 0);
        @memset(formulation, 0);
        @memset(daily_nitrogen_input_g_n, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .layer_capacity = layer_capacity, .soil = soil, .initial_urease_inhibition_fraction = initial, .current_urease_inhibition_fraction = current, .formulation = formulation, .daily_nitrogen_input_g_n = daily_nitrogen_input_g_n };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.soil);
        self.allocator.free(self.initial_urease_inhibition_fraction);
        self.allocator.free(self.current_urease_inhibition_fraction);
        self.allocator.free(self.formulation);
        self.allocator.free(self.daily_nitrogen_input_g_n);
        self.* = undefined;
    }

    pub fn index(self: State, cell: usize, layer: usize) !usize {
        if (cell >= self.cell_count or layer >= self.layer_capacity) return error.FertilizerInventoryIndexOutOfBounds;
        return cell * self.layer_capacity + layer;
    }

    pub fn resetDaily(self: *State) void {
        @memset(self.daily_nitrogen_input_g_n, 0);
    }
};

/// Applies only the nitrogen portion of one management event. The caller must
/// dispatch P, Ca, residue, and manure through their respective mass ledgers.
pub fn applyEventNitrogen(
    state: *State,
    surface_state: *surface.State,
    cell: usize,
    cell_area_m2: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    surface_litter_cover_fraction: f64,
    active_layer_thickness_m: []const f64,
    event: schedule.Event,
) !void {
    if (cell >= state.cell_count or cell >= surface_state.cells.len) return error.FertilizerInventoryIndexOutOfBounds;
    if (active_layer_thickness_m.len == 0 or active_layer_thickness_m.len > state.layer_capacity) return error.InvalidActiveSoilLayerCount;
    if (!std.math.isFinite(cell_area_m2) or cell_area_m2 <= 0 or !std.math.isFinite(nitrogen_molar_mass_g_per_mol) or nitrogen_molar_mass_g_per_mol <= 0) return error.InvalidFertilizerMassConversion;
    if (!std.math.isFinite(event.application_depth_m) or event.application_depth_m < 0 or !std.math.isFinite(surface_litter_cover_fraction) or surface_litter_cover_fraction < 0 or surface_litter_cover_fraction > 1) return error.InvalidFertilizerApplication;
    for (active_layer_thickness_m) |thickness| if (!std.math.isFinite(thickness) or thickness <= 0) return error.InvalidSoilLayerThickness;

    const conversion_mol_per_g_per_m2 = cell_area_m2 / nitrogen_molar_mass_g_per_mol;
    const n = event.nitrogen_g_per_m2;
    const broadcast = [_]f64{ n.broadcast_ammonium, n.broadcast_ammonia, n.broadcast_urea, n.broadcast_nitrate };
    const banded = [_]f64{ n.banded_ammonium, n.banded_ammonia, n.banded_urea, n.banded_nitrate };
    for (broadcast ++ banded) |amount| if (!std.math.isFinite(amount) or amount < 0) return error.InvalidFertilizerApplication;
    const total_input_g_n = (n.broadcast_ammonium + n.broadcast_ammonia + n.broadcast_urea + n.broadcast_nitrate +
        n.banded_ammonium + n.banded_ammonia + n.banded_urea + n.banded_nitrate) * cell_area_m2;
    const daily_input_next = state.daily_nitrogen_input_g_n[cell] + total_input_g_n;
    if (!std.math.isFinite(daily_input_next)) return error.FertilizerInventoryOverflow;

    const soil_layer = try layerAtDepth(active_layer_thickness_m, event.application_depth_m);
    const soil_index = try state.index(cell, soil_layer);
    var next_surface = surface_state.cells[cell];
    var next_surface_formulation = surface_state.formulation[cell];
    var next_soil = state.soil[soil_index];
    const banded_total_g_n_per_m2 = n.banded_ammonium + n.banded_ammonia + n.banded_urea + n.banded_nitrate;
    const surface_target = event.application_depth_m == 0 and
        banded_total_g_n_per_m2 == 0 and
        event.phosphorus_g_per_m2.banded_monocalcium_phosphate == 0 and
        event.calcium_carbonate_g_ca_per_m2 == 0 and
        event.calcium_sulfate_g_ca_per_m2 == 0;
    if (surface_target) {
        // HOUR1 partitions NH4, urea, and NO3 by litter cover. Broadcast NH3
        // bypasses litter and enters the upper mineral layer in full.
        next_surface.ammonium_mol_n += n.broadcast_ammonium * conversion_mol_per_g_per_m2 * surface_litter_cover_fraction;
        next_surface.urea_mol_n += n.broadcast_urea * conversion_mol_per_g_per_m2 * surface_litter_cover_fraction;
        next_surface.nitrate_mol_n += n.broadcast_nitrate * conversion_mol_per_g_per_m2 * surface_litter_cover_fraction;
        const bare_fraction = 1.0 - surface_litter_cover_fraction;
        next_soil.broadcast_ammonium_mol_n += n.broadcast_ammonium * conversion_mol_per_g_per_m2 * bare_fraction;
        next_soil.broadcast_ammonia_mol_n += n.broadcast_ammonia * conversion_mol_per_g_per_m2;
        next_soil.broadcast_urea_mol_n += n.broadcast_urea * conversion_mol_per_g_per_m2 * bare_fraction;
        next_soil.broadcast_nitrate_mol_n += n.broadcast_nitrate * conversion_mol_per_g_per_m2 * bare_fraction;
    } else {
        next_soil.broadcast_ammonium_mol_n += n.broadcast_ammonium * conversion_mol_per_g_per_m2;
        next_soil.broadcast_ammonia_mol_n += n.broadcast_ammonia * conversion_mol_per_g_per_m2;
        next_soil.broadcast_urea_mol_n += n.broadcast_urea * conversion_mol_per_g_per_m2;
        next_soil.broadcast_nitrate_mol_n += n.broadcast_nitrate * conversion_mol_per_g_per_m2;
    }
    next_soil.banded_ammonium_mol_n += n.banded_ammonium * conversion_mol_per_g_per_m2;
    next_soil.banded_ammonia_mol_n += n.banded_ammonia * conversion_mol_per_g_per_m2;
    next_soil.banded_urea_mol_n += n.banded_urea * conversion_mol_per_g_per_m2;
    next_soil.banded_nitrate_mol_n += n.banded_nitrate * conversion_mol_per_g_per_m2;

    const surface_urea_applied = surface_target and surface_litter_cover_fraction > 0 and n.broadcast_urea > 0;
    const soil_urea_applied = n.banded_urea > 0 or (!surface_target and n.broadcast_urea > 0) or (surface_target and surface_litter_cover_fraction < 1 and n.broadcast_urea > 0);
    if (surface_urea_applied) {
        next_surface.initial_urease_inhibition_fraction = 1;
        next_surface.current_urease_inhibition_fraction = 1;
        next_surface_formulation = event.fertilizer_formulation;
    }
    if (soil_urea_applied) {
        state.initial_urease_inhibition_fraction[soil_index] = 1;
        state.current_urease_inhibition_fraction[soil_index] = 1;
        state.formulation[soil_index] = event.fertilizer_formulation;
    }
    try validateFiniteInventory(next_surface, next_soil);
    surface_state.cells[cell] = next_surface;
    surface_state.formulation[cell] = next_surface_formulation;
    state.soil[soil_index] = next_soil;
    state.daily_nitrogen_input_g_n[cell] = daily_input_next;
}

fn layerAtDepth(thickness_m: []const f64, depth_m: f64) !usize {
    var lower_boundary_m: f64 = 0;
    for (thickness_m, 0..) |thickness, layer| {
        lower_boundary_m += thickness;
        if (depth_m <= lower_boundary_m) return layer;
    }
    return error.FertilizerApplicationBelowSoilProfile;
}

fn zeroSoilInventory() soil_dissolution.FertilizerState {
    return .{ .broadcast_ammonium_mol_n = 0, .broadcast_ammonia_mol_n = 0, .broadcast_urea_mol_n = 0, .broadcast_nitrate_mol_n = 0, .banded_ammonium_mol_n = 0, .banded_ammonia_mol_n = 0, .banded_urea_mol_n = 0, .banded_nitrate_mol_n = 0 };
}

fn validateFiniteInventory(surface_inventory: surface.Inventory, soil_inventory: soil_dissolution.FertilizerState) !void {
    inline for (@typeInfo(surface.Inventory).@"struct".fields) |field| if (!std.math.isFinite(@field(surface_inventory, field.name)) or @field(surface_inventory, field.name) < 0) return error.FertilizerApplicationOverflow;
    inline for (@typeInfo(soil_dissolution.FertilizerState).@"struct".fields) |field| if (!std.math.isFinite(@field(soil_inventory, field.name)) or @field(soil_inventory, field.name) < 0) return error.FertilizerApplicationOverflow;
}

test "fertilizer event conversion conserves nitrogen and honors depth and bands" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    var litter = try surface.State.init(std.testing.allocator, 1);
    defer litter.deinit();
    const event: schedule.Event = .{
        .date = .{ .day = 1, .month = 5, .year = 2001 },
        .nitrogen_g_per_m2 = .{ .broadcast_ammonium = 14, .broadcast_ammonia = 0, .broadcast_urea = 28, .broadcast_nitrate = 42, .banded_ammonium = 14, .banded_ammonia = 28, .banded_urea = 42, .banded_nitrate = 56 },
        .phosphorus_g_per_m2 = .{ .broadcast_monocalcium_phosphate = 0, .banded_monocalcium_phosphate = 0, .broadcast_hydroxyapatite = 0 },
        .calcium_carbonate_g_ca_per_m2 = 0,
        .calcium_sulfate_g_ca_per_m2 = 0,
        .plant_residue_g_per_m2 = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .manure_g_per_m2 = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .application_depth_m = 0,
        .band_row_width_m = 0.1,
        .fertilizer_formulation = 4,
        .plant_residue_type = 0,
        .manure_type = 0,
    };
    try applyEventNitrogen(&state, &litter, 0, 2, 14, 0.25, &.{ 0.1, 0.2, 0.3 }, event);
    try std.testing.expectEqual(@as(f64, 0), litter.cells[0].ammonium_mol_n + litter.cells[0].urea_mol_n + litter.cells[0].nitrate_mol_n);
    try std.testing.expectApproxEqAbs(@as(f64, 12), state.soil[0].broadcast_ammonium_mol_n + state.soil[0].broadcast_ammonia_mol_n + state.soil[0].broadcast_urea_mol_n + state.soil[0].broadcast_nitrate_mol_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 20), state.soil[0].banded_ammonium_mol_n + state.soil[0].banded_ammonia_mol_n + state.soil[0].banded_urea_mol_n + state.soil[0].banded_nitrate_mol_n, 1e-14);
    try std.testing.expectEqual(@as(f64, 448), state.daily_nitrogen_input_g_n[0]);
    state.resetDaily();
    try std.testing.expectEqual(@as(f64, 0), state.daily_nitrogen_input_g_n[0]);
    try std.testing.expectEqual(@as(u8, 0), litter.formulation[0]);
    try std.testing.expectEqual(@as(u8, 4), state.formulation[0]);
}

test "zero-depth broadcast nitrogen follows litter cover while ammonia enters topsoil" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    var litter = try surface.State.init(std.testing.allocator, 1);
    defer litter.deinit();
    const event: schedule.Event = .{
        .date = .{ .day = 1, .month = 5, .year = 2001 },
        .nitrogen_g_per_m2 = .{ .broadcast_ammonium = 14, .broadcast_ammonia = 28, .broadcast_urea = 42, .broadcast_nitrate = 56, .banded_ammonium = 0, .banded_ammonia = 0, .banded_urea = 0, .banded_nitrate = 0 },
        .phosphorus_g_per_m2 = .{ .broadcast_monocalcium_phosphate = 0, .banded_monocalcium_phosphate = 0, .broadcast_hydroxyapatite = 0 },
        .calcium_carbonate_g_ca_per_m2 = 0,
        .calcium_sulfate_g_ca_per_m2 = 0,
        .plant_residue_g_per_m2 = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .manure_g_per_m2 = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .application_depth_m = 0,
        .band_row_width_m = 0,
        .fertilizer_formulation = 1,
        .plant_residue_type = 0,
        .manure_type = 0,
    };
    try applyEventNitrogen(&state, &litter, 0, 1, 14, 0.25, &.{0.2}, event);
    try std.testing.expectApproxEqAbs(@as(f64, 2), litter.cells[0].ammonium_mol_n + litter.cells[0].urea_mol_n + litter.cells[0].nitrate_mol_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 8), state.soil[0].broadcast_ammonium_mol_n + state.soil[0].broadcast_ammonia_mol_n + state.soil[0].broadcast_urea_mol_n + state.soil[0].broadcast_nitrate_mol_n, 1e-14);
    try std.testing.expectEqual(@as(f64, 140), state.daily_nitrogen_input_g_n[0]);
}
