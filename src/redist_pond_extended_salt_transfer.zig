const std = @import("std");

/// Declaration order is the exact REDIST 8624--8691 assignment order.
pub const Species = enum(u8) {
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_silicate,
    phosphate,
    phosphoric_acid,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
};

pub const species_count = std.meta.fields(Species).len;

pub const Unit = enum { mol, g_p };

/// Legacy `H0PO4` and `H3PO4` are documented as g P in HOUR1 lines 400--401;
/// salt ions and ion-pair complexes are extensive molar amounts.
pub fn unit(species: Species) Unit {
    return switch (species) {
        .phosphate, .phosphoric_acid => .g_p,
        else => .mol,
    };
}

pub const Pools = struct {
    layer_count: usize,
    /// Species-major storage: `species * layer_count + layer`.
    amounts: []f64,

    fn index(self: Pools, species: Species, layer: usize) usize {
        return @intFromEnum(species) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8623--8692 (`ISALTG != 0`).
pub fn transfer(
    salts_enabled: bool,
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
    pools: Pools,
) !void {
    if (pools.layer_count == 0 or pools.amounts.len != species_count * pools.layer_count or
        source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer)
        return error.PondExtendedSaltTransferDimensionMismatch;
    if (!finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondExtendedSaltTransferInput;
    if (!salts_enabled) return;

    // Validate every result before any pool is mutated. Recalculation below avoids
    // fixed-size scratch allocation while retaining exact species order.
    inline for (std.meta.fields(Species)) |field| {
        const species: Species = @enumFromInt(field.value);
        const destination = pools.index(species, destination_layer);
        const source = pools.index(species, source_layer);
        const result = pools.amounts[destination] + fraction * pools.amounts[source];
        if (!std.math.isFinite(result)) return error.NonFinitePondExtendedSaltTransferResult;
    }
    inline for (std.meta.fields(Species)) |field| {
        const species: Species = @enumFromInt(field.value);
        const destination = pools.index(species, destination_layer);
        const source = pools.index(species, source_layer);
        pools.amounts[destination] = pools.amounts[destination] + fraction * pools.amounts[source];
    }
}

fn setSpecies(pools: Pools, species: Species, source: f64, destination: f64) void {
    pools.amounts[pools.index(species, 0)] = source;
    pools.amounts[pools.index(species, 1)] = destination;
}

test "REDIST extended pond salts transfer all species in declaration order" {
    var amounts: [species_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    inline for (std.meta.fields(Species), 0..) |field, ordinal| {
        const species: Species = @enumFromInt(field.value);
        setSpecies(pools, species, @floatFromInt(ordinal + 1), @floatFromInt((ordinal + 1) * 10));
    }

    try transfer(true, 0, 1, 0.25, pools);
    inline for (std.meta.fields(Species), 0..) |field, ordinal| {
        const species: Species = @enumFromInt(field.value);
        const factor: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(10.25 * factor, amounts[pools.index(species, 1)]);
        try std.testing.expectEqual(factor, amounts[pools.index(species, 0)]);
    }
}

test "REDIST extended pond salts remain unchanged when salts are disabled" {
    var amounts: [species_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    try transfer(false, 0, 1, 0.5, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST extended pond salt units distinguish free phosphate mass" {
    try std.testing.expectEqual(Unit.g_p, unit(.phosphate));
    try std.testing.expectEqual(Unit.g_p, unit(.phosphoric_acid));
    try std.testing.expectEqual(Unit.mol, unit(.sulfate));
    try std.testing.expectEqual(Unit.mol, unit(.iron_phosphate_1));
}

test "REDIST extended pond salt runtime dimensions are exact" {
    var amounts: [species_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondExtendedSaltTransferDimensionMismatch,
        transfer(true, 0, 1, 0.5, .{ .layer_count = 2, .amounts = &amounts }),
    );
}

test "REDIST extended pond salt validation is atomic" {
    var amounts: [species_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    amounts[pools.index(.potassium_sulfate, 0)] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidPondExtendedSaltTransferInput,
        transfer(true, 0, 1, 0.5, pools),
    );
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.sulfate, 1)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.magnesium_phosphate_1, 1)]);
}
