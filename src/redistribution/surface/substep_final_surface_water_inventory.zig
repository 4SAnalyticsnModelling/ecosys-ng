const std = @import("std");

pub const WaterStores = struct {
    litter_liquid_water_m3: f64,
    litter_vapor_water_equivalent_m3: f64,
    litter_ice_volume_m3: f64,
    canopy_surface_water_m3: f64,
    canopy_internal_water_m3: f64,
};

pub const WaterInventory = struct {
    landscape_water_m3: f64,
    cell_water_m3: f64,
};

pub const Contributions = struct {
    litter_water_m3: f64,
    canopy_water_m3: f64,
};

pub const Result = struct {
    inventory: WaterInventory,
    contributions: Contributions,
};

/// Exact redist.f line 5631 predicate.
pub fn isFinalSubstep(current_substep: u32, final_substep: u32) !bool {
    if (final_substep == 0 or current_substep > final_substep)
        return error.InvalidSurfaceWaterSchedule;
    return current_substep == final_substep;
}

/// Exact standalone translation of redist.f lines 5632--5635.
/// `ice_density_megagrams_m3` retains the legacy DENSI arithmetic without
/// reassociation or unit reinterpretation.
pub fn accumulate(
    initial: WaterInventory,
    stores: WaterStores,
    ice_density_megagrams_m3: f64,
) !Result {
    inline for (@typeInfo(WaterInventory).@"struct".fields) |field|
        if (!std.math.isFinite(@field(initial, field.name)))
            return error.InvalidSurfaceWaterInventory;
    inline for (@typeInfo(WaterStores).@"struct".fields) |field|
        if (!std.math.isFinite(@field(stores, field.name)))
            return error.InvalidSurfaceWaterStore;
    if (!std.math.isFinite(ice_density_megagrams_m3)) return error.InvalidIceDensity;

    const litter_water = stores.litter_liquid_water_m3 +
        stores.litter_vapor_water_equivalent_m3 +
        stores.litter_ice_volume_m3 * ice_density_megagrams_m3;
    const canopy_water = stores.canopy_surface_water_m3 + stores.canopy_internal_water_m3;
    const inventory = WaterInventory{
        .landscape_water_m3 = initial.landscape_water_m3 + litter_water + canopy_water,
        .cell_water_m3 = initial.cell_water_m3 + litter_water + canopy_water,
    };
    inline for (@typeInfo(WaterInventory).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inventory, field.name)))
            return error.NonFiniteSurfaceWaterInventory;
    if (!std.math.isFinite(litter_water) or !std.math.isFinite(canopy_water))
        return error.NonFiniteSurfaceWaterInventory;
    return .{
        .inventory = inventory,
        .contributions = .{ .litter_water_m3 = litter_water, .canopy_water_m3 = canopy_water },
    };
}

test "REDIST substep-final surface water preserves source grouping" {
    const result = try accumulate(
        .{ .landscape_water_m3 = 10.0, .cell_water_m3 = 20.0 },
        .{
            .litter_liquid_water_m3 = 1.0,
            .litter_vapor_water_equivalent_m3 = 2.0,
            .litter_ice_volume_m3 = 4.0,
            .canopy_surface_water_m3 = 5.0,
            .canopy_internal_water_m3 = 6.0,
        },
        0.5,
    );
    try std.testing.expectEqual(@as(f64, 5.0), result.contributions.litter_water_m3);
    try std.testing.expectEqual(@as(f64, 11.0), result.contributions.canopy_water_m3);
    try std.testing.expectEqual(@as(f64, 26.0), result.inventory.landscape_water_m3);
    try std.testing.expectEqual(@as(f64, 36.0), result.inventory.cell_water_m3);
}

test "REDIST substep-final water gate is independent of hour" {
    try std.testing.expect(try isFinalSubstep(4, 4));
    try std.testing.expect(!try isFinalSubstep(3, 4));
    try std.testing.expectError(error.InvalidSurfaceWaterSchedule, isFinalSubstep(5, 4));
}

test "REDIST substep-final surface water preserves signed legacy arithmetic" {
    const result = try accumulate(
        .{ .landscape_water_m3 = 0.0, .cell_water_m3 = 0.0 },
        .{
            .litter_liquid_water_m3 = 1.0,
            .litter_vapor_water_equivalent_m3 = -0.25,
            .litter_ice_volume_m3 = 0.0,
            .canopy_surface_water_m3 = 0.0,
            .canopy_internal_water_m3 = 0.0,
        },
        0.9,
    );
    try std.testing.expectEqual(@as(f64, 0.75), result.inventory.cell_water_m3);
}

test "REDIST substep-final surface water rejects invalid and overflow" {
    var stores = std.mem.zeroes(WaterStores);
    stores.canopy_internal_water_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceWaterStore,
        accumulate(.{ .landscape_water_m3 = 0.0, .cell_water_m3 = 0.0 }, stores, 0.9),
    );
    stores = std.mem.zeroes(WaterStores);
    stores.litter_liquid_water_m3 = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceWaterInventory,
        accumulate(.{ .landscape_water_m3 = std.math.floatMax(f64), .cell_water_m3 = 0.0 }, stores, 0.9),
    );
}
