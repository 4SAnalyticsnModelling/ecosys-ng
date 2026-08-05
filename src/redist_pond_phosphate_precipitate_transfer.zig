const std = @import("std");

/// Exact REDIST 9016--9035 order: non-band, then banded precipitates.
pub const PrecipitatePool = enum(u8) {
    nonband_aluminum_phosphate,
    nonband_iron_phosphate,
    nonband_dicalcium_phosphate,
    nonband_hydroxyapatite,
    nonband_monocalcium_phosphate,
    banded_aluminum_phosphate,
    banded_iron_phosphate,
    banded_dicalcium_phosphate,
    banded_hydroxyapatite,
    banded_monocalcium_phosphate,
};

pub const precipitate_pool_count = std.meta.fields(PrecipitatePool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive phosphate precipitates, mol P.
    amounts_mol_p: []f64,

    fn index(self: Pools, pool: PrecipitatePool, layer: usize) usize {
        return @intFromEnum(pool) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9016--9035 under `FX == 1.0`.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_mol_p.len != precipitate_pool_count * pools.layer_count)
        return error.PondPhosphatePrecipitateTransferDimensionMismatch;
    if (!finiteSlice(pools.amounts_mol_p) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondPhosphatePrecipitateTransferInput;
    if (fraction != 1.0) return;

    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        if (!std.math.isFinite(pools.amounts_mol_p[destination] + fraction * pools.amounts_mol_p[source]))
            return error.NonFinitePondPhosphatePrecipitateTransferResult;
    }
    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        pools.amounts_mol_p[destination] = pools.amounts_mol_p[destination] +
            fraction * pools.amounts_mol_p[source];
    }
}

test "REDIST full-fraction phosphate precipitates transfer exact pool order" {
    var amounts: [precipitate_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol_p = &amounts };
    inline for (std.meta.fields(PrecipitatePool), 0..) |field, ordinal| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        amounts[pools.index(pool, 0)] = @floatFromInt(ordinal + 1);
        amounts[pools.index(pool, 1)] = @floatFromInt((ordinal + 1) * 10);
    }
    try transfer(0, 1, 1, pools);
    inline for (std.meta.fields(PrecipitatePool), 0..) |field, ordinal| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const factor: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(11 * factor, amounts[pools.index(pool, 1)]);
        try std.testing.expectEqual(factor, amounts[pools.index(pool, 0)]);
    }
}

test "REDIST phosphate precipitates require exact full fraction" {
    var amounts: [precipitate_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol_p = &amounts };
    try transfer(0, 1, 0.999999999, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST phosphate precipitate runtime dimensions are exact" {
    var amounts: [precipitate_pool_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondPhosphatePrecipitateTransferDimensionMismatch,
        transfer(0, 1, 1, .{ .layer_count = 2, .amounts_mol_p = &amounts }),
    );
}

test "REDIST phosphate precipitate validation is atomic" {
    var amounts: [precipitate_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol_p = &amounts };
    amounts[pools.index(.banded_monocalcium_phosphate, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondPhosphatePrecipitateTransferInput, transfer(0, 1, 1, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.nonband_aluminum_phosphate, 1)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.banded_monocalcium_phosphate, 1)]);
}
