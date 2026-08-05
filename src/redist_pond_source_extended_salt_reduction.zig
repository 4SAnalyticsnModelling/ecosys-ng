const std = @import("std");
const destination_transfer = @import("redist_pond_extended_salt_transfer.zig");

pub const Species = destination_transfer.Species;
pub const Unit = destination_transfer.Unit;
pub const species_count = destination_transfer.species_count;
pub const unit = destination_transfer.unit;

pub const Pools = struct {
    layer_count: usize,
    /// Species-major storage in the units returned by `unit`.
    amounts: []f64,

    fn index(self: Pools, species: Species, layer: usize) usize {
        return @intFromEnum(species) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9128--9163 (`ISALTG != 0`).
pub fn reduce(salts_enabled: bool, source_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.amounts.len != species_count * pools.layer_count)
        return error.PondSourceExtendedSaltReductionDimensionMismatch;
    if (!finiteSlice(pools.amounts) or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceExtendedSaltReductionInput;
    if (!salts_enabled) return;

    inline for (std.meta.fields(Species)) |field| {
        const species: Species = @enumFromInt(field.value);
        const index = pools.index(species, source_layer);
        if (!std.math.isFinite(remaining_fraction * pools.amounts[index]))
            return error.NonFinitePondSourceExtendedSaltReductionResult;
    }
    inline for (std.meta.fields(Species)) |field| {
        const species: Species = @enumFromInt(field.value);
        const index = pools.index(species, source_layer);
        pools.amounts[index] = remaining_fraction * pools.amounts[index];
    }
}

test "REDIST extended source salts scale in exact shared species order" {
    var amounts: [species_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    inline for (std.meta.fields(Species), 0..) |field, ordinal| {
        const species: Species = @enumFromInt(field.value);
        amounts[pools.index(species, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[pools.index(species, 1)] = 99;
    }
    try reduce(true, 0, 0.25, pools);
    inline for (std.meta.fields(Species), 0..) |field, ordinal| {
        const species: Species = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, @floatFromInt(ordinal + 1)), amounts[pools.index(species, 0)]);
        try std.testing.expectEqual(@as(f64, 99), amounts[pools.index(species, 1)]);
    }
}

test "REDIST disabled extended salts remain unchanged" {
    var amounts: [species_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts = &amounts };
    try reduce(false, 0, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST extended source salts permit layer zero and exhaustion" {
    var amounts: [species_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts = &amounts };
    try reduce(true, 0, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "REDIST extended source salt runtime dimensions are exact" {
    var amounts: [species_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondSourceExtendedSaltReductionDimensionMismatch,
        reduce(true, 0, 0.5, .{ .layer_count = 2, .amounts = &amounts }),
    );
}

test "REDIST extended source salt validation is atomic" {
    var amounts: [species_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts = &amounts };
    amounts[pools.index(.magnesium_phosphate_1, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceExtendedSaltReductionInput, reduce(true, 0, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.sulfate, 0)]);
}
