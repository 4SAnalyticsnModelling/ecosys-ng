const std = @import("std");
const nitrogen_module = @import("../../management/fertilizer_nitrogen_inventory.zig");
const nitrogen_pools = @import("fertilizer_dissolution.zig");
const mineral_module = @import("../../management/mineral_fertilizer_inventory.zig");

/// REDIST ponding transaction for ZNH4FA...ZNO3FB and the mineral P/Ca/
/// ground-silicate stores. Both runtime owners validate completely before the
/// first source or destination field changes.
pub fn transferCellLayerFraction(
    nitrogen: *nitrogen_module.State,
    mineral: *mineral_module.State,
    cell: usize,
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
) !void {
    try validateCellLayerFraction(nitrogen, mineral, cell, source_layer, destination_layer, fraction);
    const source = cell * nitrogen.layer_capacity + source_layer;
    const destination = cell * nitrogen.layer_capacity + destination_layer;
    if (fraction == 0) return;
    transferStruct(nitrogen_pools.FertilizerState, &nitrogen.soil[source], &nitrogen.soil[destination], fraction);
    transferStruct(mineral_module.Inventory, &mineral.soil[source], &mineral.soil[destination], fraction);
}

pub fn validateCellLayerFraction(
    nitrogen: *const nitrogen_module.State,
    mineral: *const mineral_module.State,
    cell: usize,
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
) !void {
    if (nitrogen.cell_count != mineral.cell_count or nitrogen.layer_capacity != mineral.layer_capacity or cell >= nitrogen.cell_count or source_layer >= nitrogen.layer_capacity or destination_layer >= nitrogen.layer_capacity or source_layer == destination_layer) return error.FertilizerLayerRemapDimensionMismatch;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidFertilizerLayerRemapFraction;
    const source = cell * nitrogen.layer_capacity + source_layer;
    const destination = cell * nitrogen.layer_capacity + destination_layer;
    try validateStructTransfer(nitrogen_pools.FertilizerState, nitrogen.soil[source], nitrogen.soil[destination], fraction);
    try validateStructTransfer(mineral_module.Inventory, mineral.soil[source], mineral.soil[destination], fraction);
}

fn validateStructTransfer(comptime T: type, source: T, destination: T, fraction: f64) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        const source_value = @field(source, field.name);
        const destination_value = @field(destination, field.name);
        const moved = fraction * source_value;
        const next_source = source_value - moved;
        const next_destination = destination_value + moved;
        inline for (.{ source_value, destination_value, next_source, next_destination }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidFertilizerLayerRemapState;
    }
}

fn transferStruct(comptime T: type, source: *T, destination: *T, fraction: f64) void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        const moved = fraction * @field(source, field.name);
        @field(source, field.name) -= moved;
        @field(destination, field.name) += moved;
    }
}

test "REDIST fertilizer ponding moves every N P Ca and silicate inventory atomically" {
    var nitrogen = try nitrogen_module.State.init(std.testing.allocator, 1, 3);
    defer nitrogen.deinit();
    var mineral = try mineral_module.State.init(std.testing.allocator, 1, 3);
    defer mineral.deinit();
    nitrogen.soil[0].broadcast_ammonium_mol_n = 8;
    nitrogen.soil[0].banded_urea_mol_n = 4;
    mineral.soil[0].broadcast_monocalcium_phosphate_mol = 12;
    mineral.soil[0].calcite_mol = 16;
    mineral.soil[0].potassium_ground_silicate_mol = 20;
    nitrogen.current_urease_inhibition_fraction[0] = 0.75;
    nitrogen.formulation[0] = 11;
    try transferCellLayerFraction(&nitrogen, &mineral, 0, 0, 1, 0.25);
    try std.testing.expectEqual(@as(f64, 6), nitrogen.soil[0].broadcast_ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 2), nitrogen.soil[1].broadcast_ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 1), nitrogen.soil[1].banded_urea_mol_n);
    try std.testing.expectEqual(@as(f64, 3), mineral.soil[1].broadcast_monocalcium_phosphate_mol);
    try std.testing.expectEqual(@as(f64, 4), mineral.soil[1].calcite_mol);
    try std.testing.expectEqual(@as(f64, 5), mineral.soil[1].potassium_ground_silicate_mol);
    // REDIST does not move inhibitor/formulation control state.
    try std.testing.expectEqual(@as(f64, 0.75), nitrogen.current_urease_inhibition_fraction[0]);
    try std.testing.expectEqual(@as(f64, 0), nitrogen.current_urease_inhibition_fraction[1]);
    try std.testing.expectEqual(@as(u8, 11), nitrogen.formulation[0]);
}

test "REDIST fertilizer ponding rolls back both owners on a late invalid mineral" {
    var nitrogen = try nitrogen_module.State.init(std.testing.allocator, 1, 2);
    defer nitrogen.deinit();
    var mineral = try mineral_module.State.init(std.testing.allocator, 1, 2);
    defer mineral.deinit();
    nitrogen.soil[0].broadcast_nitrate_mol_n = 9;
    mineral.soil[0].potassium_ground_silicate_mol = std.math.nan(f64);
    try std.testing.expectError(error.InvalidFertilizerLayerRemapState, transferCellLayerFraction(&nitrogen, &mineral, 0, 0, 1, 0.5));
    try std.testing.expectEqual(@as(f64, 9), nitrogen.soil[0].broadcast_nitrate_mol_n);
    try std.testing.expectEqual(@as(f64, 0), nitrogen.soil[1].broadcast_nitrate_mol_n);
}
