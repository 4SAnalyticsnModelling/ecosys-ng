const std = @import("std");

pub const Pools = struct {
    hydrogen_mol: []f64, // ZHY
    hydroxide_mol: []f64, // ZOH
    aluminum_mol: []f64, // ZAL
    iron_mol: []f64, // ZFE
    calcium_mol: []f64, // ZCA
    magnesium_mol: []f64, // ZMG
    sodium_mol: []f64, // ZNA
    potassium_mol: []f64, // ZKA
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8607--8622: always-active pond salts.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const len = pools.hydrogen_mol.len;
    if (len == 0 or source_layer >= len or destination_layer >= len or source_layer == destination_layer or
        pools.hydroxide_mol.len != len or pools.aluminum_mol.len != len or pools.iron_mol.len != len or
        pools.calcium_mol.len != len or pools.magnesium_mol.len != len or pools.sodium_mol.len != len or
        pools.potassium_mol.len != len)
        return error.PondNonbandBaseSaltTransferDimensionMismatch;
    inline for (.{ pools.hydrogen_mol, pools.hydroxide_mol, pools.aluminum_mol, pools.iron_mol, pools.calcium_mol, pools.magnesium_mol, pools.sodium_mol, pools.potassium_mol }) |values|
        if (!finiteSlice(values)) return error.InvalidPondNonbandBaseSaltTransferInput;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondNonbandBaseSaltTransferInput;

    const hydrogen_mol = pools.hydrogen_mol[destination_layer] + fraction * pools.hydrogen_mol[source_layer];
    const hydroxide_mol = pools.hydroxide_mol[destination_layer] + fraction * pools.hydroxide_mol[source_layer];
    const aluminum_mol = pools.aluminum_mol[destination_layer] + fraction * pools.aluminum_mol[source_layer];
    const iron_mol = pools.iron_mol[destination_layer] + fraction * pools.iron_mol[source_layer];
    const calcium_mol = pools.calcium_mol[destination_layer] + fraction * pools.calcium_mol[source_layer];
    const magnesium_mol = pools.magnesium_mol[destination_layer] + fraction * pools.magnesium_mol[source_layer];
    const sodium_mol = pools.sodium_mol[destination_layer] + fraction * pools.sodium_mol[source_layer];
    const potassium_mol = pools.potassium_mol[destination_layer] + fraction * pools.potassium_mol[source_layer];
    inline for (.{ hydrogen_mol, hydroxide_mol, aluminum_mol, iron_mol, calcium_mol, magnesium_mol, sodium_mol, potassium_mol }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondNonbandBaseSaltTransferResult;

    pools.hydrogen_mol[destination_layer] = hydrogen_mol;
    pools.hydroxide_mol[destination_layer] = hydroxide_mol;
    pools.aluminum_mol[destination_layer] = aluminum_mol;
    pools.iron_mol[destination_layer] = iron_mol;
    pools.calcium_mol[destination_layer] = calcium_mol;
    pools.magnesium_mol[destination_layer] = magnesium_mol;
    pools.sodium_mol[destination_layer] = sodium_mol;
    pools.potassium_mol[destination_layer] = potassium_mol;
}

const Fixture = struct {
    hydrogen: [3]f64 = .{ 2, 1, 0 },
    hydroxide: [3]f64 = .{ 4, 2, 0 },
    aluminum: [3]f64 = .{ 6, 3, 0 },
    iron: [3]f64 = .{ 8, 4, 0 },
    calcium: [3]f64 = .{ 10, 5, 0 },
    magnesium: [3]f64 = .{ 12, 6, 0 },
    sodium: [3]f64 = .{ 14, 7, 0 },
    potassium: [3]f64 = .{ 16, 8, 0 },

    fn pools(self: *Fixture) Pools {
        return .{
            .hydrogen_mol = &self.hydrogen,
            .hydroxide_mol = &self.hydroxide,
            .aluminum_mol = &self.aluminum,
            .iron_mol = &self.iron,
            .calcium_mol = &self.calcium,
            .magnesium_mol = &self.magnesium,
            .sodium_mol = &self.sodium,
            .potassium_mol = &self.potassium,
        };
    }
};

test "REDIST pond base salts transfer in exact source order" {
    var fixture = Fixture{};
    try transfer(0, 1, 0.25, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1.5), fixture.hydrogen[1]);
    try std.testing.expectEqual(@as(f64, 3), fixture.hydroxide[1]);
    try std.testing.expectEqual(@as(f64, 4.5), fixture.aluminum[1]);
    try std.testing.expectEqual(@as(f64, 6), fixture.iron[1]);
    try std.testing.expectEqual(@as(f64, 7.5), fixture.calcium[1]);
    try std.testing.expectEqual(@as(f64, 9), fixture.magnesium[1]);
    try std.testing.expectEqual(@as(f64, 10.5), fixture.sodium[1]);
    try std.testing.expectEqual(@as(f64, 12), fixture.potassium[1]);
    try std.testing.expectEqual(@as(f64, 2), fixture.hydrogen[0]);
}

test "REDIST zero fraction preserves pond base salts" {
    var fixture = Fixture{};
    try transfer(0, 1, 0, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.hydrogen[1]);
    try std.testing.expectEqual(@as(f64, 8), fixture.potassium[1]);
}

test "REDIST pond base salt runtime dimensions must agree" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.iron_mol = fixture.iron[0..2];
    try std.testing.expectError(
        error.PondNonbandBaseSaltTransferDimensionMismatch,
        transfer(0, 1, 0.5, pools),
    );
}

test "REDIST pond base salt validation is atomic" {
    var fixture = Fixture{};
    fixture.sodium[0] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidPondNonbandBaseSaltTransferInput,
        transfer(0, 1, 0.5, fixture.pools()),
    );
    try std.testing.expectEqual(@as(f64, 1), fixture.hydrogen[1]);
    try std.testing.expectEqual(@as(f64, 8), fixture.potassium[1]);
}
