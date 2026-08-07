const std = @import("std");

pub const MineralInventories = struct {
    gibbsite_mol: f64,
    iron_hydroxide_mol: f64,
    calcite_mol: f64,
    gypsum_mol: f64,
    aluminum_phosphate_mol: f64,
    iron_phosphate_mol: f64,
    monocalcium_phosphate_mol: f64,
    dicalcium_phosphate_mol: f64,
    hydroxyapatite_mol: f64,
};

pub const MineralConcentrations = struct {
    gibbsite_mol_per_m3: f64,
    iron_hydroxide_mol_per_m3: f64,
    calcite_mol_per_m3: f64,
    gypsum_mol_per_m3: f64,
    aluminum_phosphate_mol_per_m3: f64,
    iron_phosphate_mol_per_m3: f64,
    monocalcium_phosphate_mol_per_m3: f64,
    dicalcium_phosphate_mol_per_m3: f64,
    hydroxyapatite_mol_per_m3: f64,
};

pub const Inputs = struct {
    litter_water_volume_m3: f64,
    active_water_volume_threshold_m3: f64,
    inventories: MineralInventories,
};

pub const Result = union(enum) {
    inactive,
    wet: MineralConcentrations,
};

/// Direct source-order translation of SOLUTE.F lines 4163--4192.
///
/// The caller selects one runtime horizontal cell. This pure kernel performs
/// no allocation and does not mutate the authoritative mineral inventories.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);
    if (inputs.litter_water_volume_m3 <=
        inputs.active_water_volume_threshold_m3)
    {
        return .inactive;
    }

    const water = inputs.litter_water_volume_m3;
    const inventory = inputs.inventories;
    // SOLUTE.F 4173--4181. Preserve mineral and operation order.
    const concentrations: MineralConcentrations = .{
        .gibbsite_mol_per_m3 = @max(0.0, inventory.gibbsite_mol) / water,
        .iron_hydroxide_mol_per_m3 = @max(0.0, inventory.iron_hydroxide_mol) / water,
        .calcite_mol_per_m3 = @max(0.0, inventory.calcite_mol) / water,
        .gypsum_mol_per_m3 = @max(0.0, inventory.gypsum_mol) / water,
        .aluminum_phosphate_mol_per_m3 = @max(0.0, inventory.aluminum_phosphate_mol) / water,
        .iron_phosphate_mol_per_m3 = @max(0.0, inventory.iron_phosphate_mol) / water,
        .monocalcium_phosphate_mol_per_m3 = @max(0.0, inventory.monocalcium_phosphate_mol) / water,
        .dicalcium_phosphate_mol_per_m3 = @max(0.0, inventory.dicalcium_phosphate_mol) / water,
        .hydroxyapatite_mol_per_m3 = @max(0.0, inventory.hydroxyapatite_mol) / water,
    };
    inline for (@typeInfo(MineralConcentrations).@"struct".fields) |field| {
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterMineralConcentration;
    }
    return .{ .wet = concentrations };
}

fn validateInputs(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.litter_water_volume_m3) or
        inputs.litter_water_volume_m3 < 0 or
        !std.math.isFinite(inputs.active_water_volume_threshold_m3) or
        inputs.active_water_volume_threshold_m3 < 0)
    {
        return error.InvalidSurfaceLitterMineralConcentrationInput;
    }
    inline for (@typeInfo(MineralInventories).@"struct".fields) |field| {
        const value = @field(inputs.inventories, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterMineralConcentrationInput;
    }
}

fn testInputs() Inputs {
    return .{
        .litter_water_volume_m3 = 2,
        .active_water_volume_threshold_m3 = 1.0e-12,
        .inventories = .{
            .gibbsite_mol = 2,
            .iron_hydroxide_mol = 4,
            .calcite_mol = 6,
            .gypsum_mol = 8,
            .aluminum_phosphate_mol = 10,
            .iron_phosphate_mol = 12,
            .monocalcium_phosphate_mol = 14,
            .dicalcium_phosphate_mol = 16,
            .hydroxyapatite_mol = 18,
        },
    };
}

test "SOLUTE surface mineral reconstruction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);
    const wet = switch (result) {
        .inactive => return error.ExpectedWetSurfaceLitter,
        .wet => |value| value,
    };

    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.gibbsite_mol) /
            inputs.litter_water_volume_m3,
        wet.gibbsite_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.iron_hydroxide_mol) /
            inputs.litter_water_volume_m3,
        wet.iron_hydroxide_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.calcite_mol) /
            inputs.litter_water_volume_m3,
        wet.calcite_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.gypsum_mol) /
            inputs.litter_water_volume_m3,
        wet.gypsum_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.aluminum_phosphate_mol) /
            inputs.litter_water_volume_m3,
        wet.aluminum_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.iron_phosphate_mol) /
            inputs.litter_water_volume_m3,
        wet.iron_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.monocalcium_phosphate_mol) /
            inputs.litter_water_volume_m3,
        wet.monocalcium_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.dicalcium_phosphate_mol) /
            inputs.litter_water_volume_m3,
        wet.dicalcium_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(0.0, inputs.inventories.hydroxyapatite_mol) /
            inputs.litter_water_volume_m3,
        wet.hydroxyapatite_mol_per_m3,
    );
}

test "surface mineral concentrations reconstruct every extensive inventory" {
    const inputs = testInputs();
    const result = try calculate(inputs);
    const wet = result.wet;

    inline for (@typeInfo(MineralInventories).@"struct".fields, 0..) |field, index| {
        const concentration_field =
            @typeInfo(MineralConcentrations).@"struct".fields[index];
        try std.testing.expectApproxEqAbs(
            @field(inputs.inventories, field.name),
            @field(wet, concentration_field.name) *
                inputs.litter_water_volume_m3,
            1.0e-14,
        );
    }
}

test "surface mineral reconstruction conserves inventory as water changes" {
    const inputs = testInputs();
    const initial = (try calculate(inputs)).wet;
    var reduced_water = inputs;
    reduced_water.litter_water_volume_m3 = 0.5;
    const concentrated = (try calculate(reduced_water)).wet;

    try std.testing.expectEqual(
        initial.calcite_mol_per_m3 * inputs.litter_water_volume_m3,
        concentrated.calcite_mol_per_m3 *
            reduced_water.litter_water_volume_m3,
    );
    try std.testing.expectEqual(
        initial.hydroxyapatite_mol_per_m3 *
            inputs.litter_water_volume_m3,
        concentrated.hydroxyapatite_mol_per_m3 *
            reduced_water.litter_water_volume_m3,
    );
}

test "surface mineral water threshold is strict" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = inputs.active_water_volume_threshold_m3;
    try std.testing.expectEqual(Result.inactive, try calculate(inputs));

    inputs.litter_water_volume_m3 =
        inputs.active_water_volume_threshold_m3 * 2;
    const result = try calculate(inputs);
    try std.testing.expect(result == .wet);
}

test "surface mineral reconstruction rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.inventories.calcite_mol = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterMineralConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.inventories.gibbsite_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterMineralConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.litter_water_volume_m3 = std.math.floatMin(f64);
    inputs.active_water_volume_threshold_m3 = 0;
    inputs.inventories.gibbsite_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterMineralConcentration,
        calculate(inputs),
    );
}
