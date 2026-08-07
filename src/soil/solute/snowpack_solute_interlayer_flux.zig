const std = @import("std");

/// Exact snowpack salt topology: 33 salt/complex fields followed by eight
/// phosphate fields. H4SiO4 is absent.
pub const species_per_layer = 41;

pub const InterlayerRouting = struct {
    /// Current snow heat capacity at this substep, MJ K-1 per layer.
    heat_capacity_megajoules_per_k: []const f64,
    /// Cell-specific minimum active snow heat capacity, MJ K-1.
    minimum_heat_capacity_megajoules_per_k: f64,
    /// Liquid snow water available in each source layer, m3.
    liquid_water_m3: []const f64,
    /// Water flux indexed by destination layer, m3 step-1.
    downward_water_flux_m3_per_step: []const f64,
    /// Layer-major snow solute inventory, mol layer-1.
    inventory_mol_per_layer: []const f64,
    /// Layer-major instantaneous interlayer solute flux, mol step-1.
    flux_mol_per_step: []f64,
    /// Layer-major accumulated interlayer solute flux, mol step-1.
    accumulated_flux_mol_per_step: []f64,
};

/// Exact compatibility translation of the interlayer branch in TRNSFRS.F
/// lines 2038--2172. The surface-discharge branch beginning at line 2186 is
/// intentionally caller-owned.
///
/// Legacy activity is determined by `VHCPWM > VHCPWX`, not a boolean mask.
/// Candidate arrays are runtime allocated so late validation failures leave
/// both mutable flux arrays unchanged.
pub fn route(allocator: std.mem.Allocator, routing: InterlayerRouting) !void {
    const layer_count = routing.heat_capacity_megajoules_per_k.len;
    if (routing.liquid_water_m3.len != layer_count or
        routing.downward_water_flux_m3_per_step.len != layer_count)
        return error.SnowpackInterlayerPhysicalDimensionMismatch;
    const solute_len = std.math.mul(usize, layer_count, species_per_layer) catch
        return error.SnowpackInterlayerSoluteDimensionOverflow;
    if (routing.inventory_mol_per_layer.len != solute_len or
        routing.flux_mol_per_step.len != solute_len or
        routing.accumulated_flux_mol_per_step.len != solute_len)
        return error.SnowpackInterlayerSoluteDimensionMismatch;
    if (!std.math.isFinite(routing.minimum_heat_capacity_megajoules_per_k) or
        routing.minimum_heat_capacity_megajoules_per_k < 0)
        return error.NonFiniteSnowpackInterlayerInput;

    const candidate_flux = try allocator.dupe(f64, routing.flux_mol_per_step);
    defer allocator.free(candidate_flux);
    const candidate_accumulated = try allocator.dupe(f64, routing.accumulated_flux_mol_per_step);
    defer allocator.free(candidate_accumulated);

    for (0..layer_count) |layer| {
        const heat_capacity = routing.heat_capacity_megajoules_per_k[layer];
        const water_m3 = routing.liquid_water_m3[layer];
        const water_flux_m3 = routing.downward_water_flux_m3_per_step[layer];
        if (!std.math.isFinite(heat_capacity) or heat_capacity < 0 or
            !std.math.isFinite(water_m3) or water_m3 < 0 or
            !std.math.isFinite(water_flux_m3))
            return error.NonFiniteSnowpackInterlayerInput;
        if (heat_capacity <= routing.minimum_heat_capacity_megajoules_per_k) continue;

        const destination = @min(layer_count - 1, layer + 1);
        if (layer + 1 < layer_count and
            routing.heat_capacity_megajoules_per_k[destination] > routing.minimum_heat_capacity_megajoules_per_k)
        {
            const fraction = if (water_m3 > 0)
                std.math.clamp(routing.downward_water_flux_m3_per_step[destination] / water_m3, 0, 1)
            else
                1;
            for (0..species_per_layer) |species| {
                const source_index = layer * species_per_layer + species;
                const destination_index = destination * species_per_layer + species;
                const inventory_mol = routing.inventory_mol_per_layer[source_index];
                if (!std.math.isFinite(inventory_mol)) return error.NonFiniteSnowpackInterlayerInput;
                const flux_mol = inventory_mol * fraction;
                const accumulated_mol = candidate_accumulated[destination_index] + flux_mol;
                if (!std.math.isFinite(flux_mol) or !std.math.isFinite(accumulated_mol))
                    return error.NonFiniteSnowpackInterlayerResult;
                candidate_flux[destination_index] = flux_mol;
                candidate_accumulated[destination_index] = accumulated_mol;
            }
        } else if (layer + 1 < layer_count) {
            @memset(candidate_flux[destination * species_per_layer .. (destination + 1) * species_per_layer], 0);
        }
    }
    @memcpy(routing.flux_mol_per_step, candidate_flux);
    @memcpy(routing.accumulated_flux_mol_per_step, candidate_accumulated);
}

test "TRNSFRS routes and accumulates 41 species to an active lower layer" {
    const layers = 2;
    const heat = [_]f64{ 2, 2 };
    const water = [_]f64{ 4, 1 };
    const water_flux = [_]f64{ 0, 1 };
    var inventory = [_]f64{0} ** (layers * species_per_layer);
    for (inventory[0..species_per_layer], 0..) |*value, species| value.* = @floatFromInt(species + 1);
    var flux = [_]f64{9} ** (layers * species_per_layer);
    var accumulated = [_]f64{1} ** (layers * species_per_layer);
    try route(std.testing.allocator, .{
        .heat_capacity_megajoules_per_k = &heat,
        .minimum_heat_capacity_megajoules_per_k = 1,
        .liquid_water_m3 = &water,
        .downward_water_flux_m3_per_step = &water_flux,
        .inventory_mol_per_layer = &inventory,
        .flux_mol_per_step = &flux,
        .accumulated_flux_mol_per_step = &accumulated,
    });
    try std.testing.expectEqual(@as(f64, 0.25), flux[species_per_layer]);
    try std.testing.expectEqual(@as(f64, 1.25), accumulated[species_per_layer]);
    try std.testing.expectEqual(@as(f64, 10.25), flux[2 * species_per_layer - 1]);
}

test "TRNSFRS dry source uses unit transport fraction" {
    const heat = [_]f64{ 2, 2 };
    const water = [_]f64{ 0, 0 };
    const water_flux = [_]f64{ 0, 99 };
    const inventory = [_]f64{2} ** (2 * species_per_layer);
    var flux = [_]f64{0} ** (2 * species_per_layer);
    var accumulated = [_]f64{0} ** (2 * species_per_layer);
    try route(std.testing.allocator, .{ .heat_capacity_megajoules_per_k = &heat, .minimum_heat_capacity_megajoules_per_k = 1, .liquid_water_m3 = &water, .downward_water_flux_m3_per_step = &water_flux, .inventory_mol_per_layer = &inventory, .flux_mol_per_step = &flux, .accumulated_flux_mol_per_step = &accumulated });
    try std.testing.expectEqual(@as(f64, 2), flux[species_per_layer]);
}

test "inactive lower layer clears its instantaneous destination flux" {
    const heat = [_]f64{ 2, 0 };
    const water = [_]f64{ 1, 0 };
    const water_flux = [_]f64{ 0, 1 };
    const inventory = [_]f64{2} ** (2 * species_per_layer);
    var flux = [_]f64{9} ** (2 * species_per_layer);
    var accumulated = [_]f64{7} ** (2 * species_per_layer);
    try route(std.testing.allocator, .{ .heat_capacity_megajoules_per_k = &heat, .minimum_heat_capacity_megajoules_per_k = 1, .liquid_water_m3 = &water, .downward_water_flux_m3_per_step = &water_flux, .inventory_mol_per_layer = &inventory, .flux_mol_per_step = &flux, .accumulated_flux_mol_per_step = &accumulated });
    try std.testing.expectEqual(@as(f64, 0), flux[species_per_layer]);
    try std.testing.expectEqual(@as(f64, 7), accumulated[species_per_layer]);
}

test "zero runtime snow layers preserve zero-trip behavior" {
    const empty_const: [0]f64 = .{};
    var empty: [0]f64 = .{};
    try route(std.testing.allocator, .{ .heat_capacity_megajoules_per_k = &empty_const, .minimum_heat_capacity_megajoules_per_k = 0, .liquid_water_m3 = &empty_const, .downward_water_flux_m3_per_step = &empty_const, .inventory_mol_per_layer = &empty_const, .flux_mol_per_step = &empty, .accumulated_flux_mol_per_step = &empty });
}

test "late invalid inventory leaves both outputs atomic" {
    const heat = [_]f64{ 2, 2 };
    const water = [_]f64{ 1, 1 };
    const water_flux = [_]f64{ 0, 1 };
    var inventory = [_]f64{2} ** (2 * species_per_layer);
    inventory[species_per_layer - 1] = std.math.inf(f64);
    var flux = [_]f64{9} ** (2 * species_per_layer);
    var accumulated = [_]f64{7} ** (2 * species_per_layer);
    try std.testing.expectError(error.NonFiniteSnowpackInterlayerInput, route(std.testing.allocator, .{ .heat_capacity_megajoules_per_k = &heat, .minimum_heat_capacity_megajoules_per_k = 1, .liquid_water_m3 = &water, .downward_water_flux_m3_per_step = &water_flux, .inventory_mol_per_layer = &inventory, .flux_mol_per_step = &flux, .accumulated_flux_mol_per_step = &accumulated }));
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** (2 * species_per_layer)), &flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{7} ** (2 * species_per_layer)), &accumulated);
}
