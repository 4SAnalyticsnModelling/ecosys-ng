const std = @import("std");

pub const Parameters = struct {
    initial_proton_balance_concentration_mol_per_m3: f64,
};

/// Direct translation of STARTE line 1518. `ZHYSI` is an extensive proton
/// balance inventory (mol H+) initialized from field-capacity water volume.
/// All layers are validated before any destination value is changed.
pub fn initialize(
    proton_balance_inventory_mol: []f64,
    field_capacity_water_m3: []const f64,
    parameters: Parameters,
) !void {
    if (proton_balance_inventory_mol.len != field_capacity_water_m3.len or
        proton_balance_inventory_mol.len == 0)
        return error.SoilProtonBalanceDimensionMismatch;
    if (!std.math.isFinite(
        parameters.initial_proton_balance_concentration_mol_per_m3,
    ))
        return error.NonFiniteSoilProtonBalanceParameter;
    if (parameters.initial_proton_balance_concentration_mol_per_m3 < 0)
        return error.InvalidSoilProtonBalanceParameter;
    for (field_capacity_water_m3) |water_m3| {
        if (!std.math.isFinite(water_m3))
            return error.NonFiniteSoilProtonBalanceInput;
        if (water_m3 < 0) return error.InvalidSoilProtonBalanceInput;
        const inventory = parameters.initial_proton_balance_concentration_mol_per_m3 *
            water_m3;
        if (!std.math.isFinite(inventory))
            return error.NonFiniteSoilProtonBalanceResult;
    }
    for (proton_balance_inventory_mol, field_capacity_water_m3) |*inventory, water_m3|
        inventory.* = parameters.initial_proton_balance_concentration_mol_per_m3 *
            water_m3;
}

fn sourceParameters() Parameters {
    return .{ .initial_proton_balance_concentration_mol_per_m3 = 1.0e-3 };
}

test "STARTE proton balance initialization supports runtime layer count" {
    const water_m3 = [_]f64{ 0, 0.1, 0.25, 0.5, 1, 2, 4 };
    var inventory_mol = [_]f64{9} ** water_m3.len;
    try initialize(&inventory_mol, &water_m3, sourceParameters());
    for (inventory_mol, water_m3) |inventory, water|
        try std.testing.expectEqual(1.0e-3 * water, inventory);
}

test "STARTE proton balance concentration is runtime configurable" {
    var inventory_mol = [_]f64{0};
    try initialize(
        &inventory_mol,
        &.{2},
        .{ .initial_proton_balance_concentration_mol_per_m3 = 0.25 },
    );
    try std.testing.expectEqual(@as(f64, 0.5), inventory_mol[0]);
}

test "STARTE proton balance preflight preserves destination on failure" {
    var inventory_mol = [_]f64{ 7, 8 };
    try std.testing.expectError(
        error.NonFiniteSoilProtonBalanceInput,
        initialize(
            &inventory_mol,
            &.{ 1, std.math.nan(f64) },
            sourceParameters(),
        ),
    );
    try std.testing.expectEqualSlices(f64, &.{ 7, 8 }, &inventory_mol);
}
