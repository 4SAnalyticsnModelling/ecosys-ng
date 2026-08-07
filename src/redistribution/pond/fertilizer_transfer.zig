const std = @import("std");

/// Exact REDIST 8967--8982 assignment order.
pub const FertilizerPool = enum(u8) {
    broadcast_ammonium,
    broadcast_ammonia,
    broadcast_urea,
    broadcast_nitrate,
    banded_ammonium,
    banded_ammonia,
    banded_urea,
    banded_nitrate,
};

pub const fertilizer_pool_count = std.meta.fields(FertilizerPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive fertilizer amounts, mol N.
    amounts_mol_n: []f64,

    fn index(self: Pools, pool: FertilizerPool, layer: usize) usize {
        return @intFromEnum(pool) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8967--8982 under the enclosing `FX == 1.0`.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_mol_n.len != fertilizer_pool_count * pools.layer_count)
        return error.PondFertilizerTransferDimensionMismatch;
    if (!finiteSlice(pools.amounts_mol_n) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondFertilizerTransferInput;
    if (fraction != 1.0) return;

    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        if (!std.math.isFinite(pools.amounts_mol_n[destination] + fraction * pools.amounts_mol_n[source]))
            return error.NonFinitePondFertilizerTransferResult;
    }
    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        pools.amounts_mol_n[destination] = pools.amounts_mol_n[destination] +
            fraction * pools.amounts_mol_n[source];
    }
}

test "REDIST full-fraction pond fertilizer transfers exact pool order" {
    var amounts: [fertilizer_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts_mol_n = &amounts };
    inline for (std.meta.fields(FertilizerPool), 0..) |field, ordinal| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        amounts[pools.index(pool, 0)] = @floatFromInt(ordinal + 1);
        amounts[pools.index(pool, 1)] = @floatFromInt((ordinal + 1) * 10);
    }
    try transfer(0, 1, 1, pools);
    inline for (std.meta.fields(FertilizerPool), 0..) |field, ordinal| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const factor: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(11 * factor, amounts[pools.index(pool, 1)]);
        try std.testing.expectEqual(factor, amounts[pools.index(pool, 0)]);
    }
}

test "REDIST pond fertilizer permits source layer zero" {
    var amounts: [fertilizer_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol_n = &amounts };
    try transfer(0, 1, 1, pools);
    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 2), amounts[pools.index(pool, 1)]);
    }
}

test "REDIST pond fertilizer requires exact full fraction" {
    var amounts: [fertilizer_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol_n = &amounts };
    try transfer(0, 1, 0.999999999, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST pond fertilizer runtime dimensions are exact" {
    var amounts: [fertilizer_pool_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondFertilizerTransferDimensionMismatch,
        transfer(0, 1, 1, .{ .layer_count = 2, .amounts_mol_n = &amounts }),
    );
}

test "REDIST pond fertilizer validation is atomic" {
    var amounts: [fertilizer_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol_n = &amounts };
    amounts[pools.index(.banded_nitrate, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondFertilizerTransferInput, transfer(0, 1, 1, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.broadcast_ammonium, 1)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.banded_nitrate, 1)]);
}
