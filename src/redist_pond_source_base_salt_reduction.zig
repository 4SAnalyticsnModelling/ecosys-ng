const std = @import("std");

/// Exact REDIST 9120--9127 order.
pub const SaltPool = enum(u8) {
    hydrogen,
    hydroxide,
    aluminum,
    iron,
    calcium,
    magnesium,
    sodium,
    potassium,
};

pub const salt_pool_count = std.meta.fields(SaltPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive salt amounts, mol.
    amounts_mol: []f64,

    fn index(self: Pools, pool: SaltPool, layer: usize) usize {
        return @intFromEnum(pool) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9120--9127: scale base source salts by FY.
pub fn reduce(source_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.amounts_mol.len != salt_pool_count * pools.layer_count)
        return error.PondSourceBaseSaltReductionDimensionMismatch;
    if (!finiteSlice(pools.amounts_mol) or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceBaseSaltReductionInput;

    inline for (std.meta.fields(SaltPool)) |field| {
        const pool: SaltPool = @enumFromInt(field.value);
        const index = pools.index(pool, source_layer);
        if (!std.math.isFinite(remaining_fraction * pools.amounts_mol[index]))
            return error.NonFinitePondSourceBaseSaltReductionResult;
    }
    inline for (std.meta.fields(SaltPool)) |field| {
        const pool: SaltPool = @enumFromInt(field.value);
        const index = pools.index(pool, source_layer);
        pools.amounts_mol[index] = remaining_fraction * pools.amounts_mol[index];
    }
}

test "REDIST source base salts scale in exact legacy order" {
    var amounts: [salt_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol = &amounts };
    inline for (std.meta.fields(SaltPool), 0..) |field, ordinal| {
        const pool: SaltPool = @enumFromInt(field.value);
        amounts[pools.index(pool, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[pools.index(pool, 1)] = 99;
    }
    try reduce(0, 0.25, pools);
    inline for (std.meta.fields(SaltPool), 0..) |field, ordinal| {
        const pool: SaltPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, @floatFromInt(ordinal + 1)), amounts[pools.index(pool, 0)]);
        try std.testing.expectEqual(@as(f64, 99), amounts[pools.index(pool, 1)]);
    }
}

test "REDIST source base salts permit layer zero and complete exhaustion" {
    var amounts: [salt_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol = &amounts };
    try reduce(0, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "REDIST source base salt runtime dimensions are exact" {
    var amounts: [salt_pool_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondSourceBaseSaltReductionDimensionMismatch,
        reduce(0, 0.5, .{ .layer_count = 2, .amounts_mol = &amounts }),
    );
}

test "REDIST source base salt validation is atomic" {
    var amounts: [salt_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol = &amounts };
    amounts[pools.index(.potassium, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceBaseSaltReductionInput, reduce(0, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.hydrogen, 0)]);
}
