const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

/// Named anion values; units are supplied by the owning aggregate.
pub const Anions = struct {
    sulfate: f64,
    chloride: f64,
    carbonate: f64,
    bicarbonate: f64,
};

pub const AluminumComplexes = struct {
    monohydroxide: f64,
    dihydroxide: f64,
    trihydroxide: f64,
    tetrahydroxide: f64,
    sulfate: f64,
};

pub const IronComplexes = struct {
    monohydroxide: f64,
    dihydroxide: f64,
    trihydroxide: f64,
    tetrahydroxide: f64,
    sulfate: f64,
};

pub const BaseCationComplexes = struct {
    calcium_hydroxide: f64,
    calcium_carbonate: f64,
    calcium_bicarbonate: f64,
    calcium_sulfate: f64,
    magnesium_hydroxide: f64,
    magnesium_carbonate: f64,
    magnesium_bicarbonate: f64,
    magnesium_sulfate: f64,
    sodium_carbonate: f64,
    sodium_sulfate: f64,
    potassium_sulfate: f64,
};

pub const PhosphateComplexes = struct {
    phosphate: f64,
    phosphoric_acid: f64,
    iron_monophosphate: f64,
    iron_diphosphate: f64,
    calcium_phosphate: f64,
    calcium_hydrogen_phosphate: f64,
    calcium_dihydrogen_phosphate: f64,
    magnesium_hydrogen_phosphate: f64,
};

/// Every leaf value is an aqueous concentration in mol m-3 (mol P m-3 for
/// phosphate-bearing species).
pub const Concentrations = struct {
    anions: Anions,
    aluminum: AluminumComplexes,
    iron: IronComplexes,
    base_cations: BaseCationComplexes,
    phosphate: PhosphateComplexes,
};

/// Every leaf value is an extensive inventory in mol (mol P for
/// phosphate-bearing species).
pub const Inventories = struct {
    anions: Anions,
    aluminum: AluminumComplexes,
    iron: IronComplexes,
    base_cations: BaseCationComplexes,
    phosphate: PhosphateComplexes,
};

pub const State = struct {
    inventories: Inventories,
    proton_balance_mol: f64,
};

fn validateAndScale(comptime T: type, input: T, water_m3: f64) !T {
    var result: T = undefined;
    inline for (std.meta.fields(T)) |field| {
        const concentration = @field(input, field.name);
        if (!std.math.isFinite(concentration))
            return error.NonFiniteSurfaceIonComplexConcentration;
        if (concentration < 0)
            return error.InvalidSurfaceIonComplexConcentration;
        const inventory = concentration * water_m3;
        if (!std.math.isFinite(inventory))
            return error.NonFiniteSurfaceIonComplexInventory;
        @field(result, field.name) = inventory;
    }
    return result;
}

/// Direct translation of STARTE lines 1981--2050 for one surface litter cell.
/// In dynamic mode, topsoil concentrations are converted to extensive litter
/// inventories. The static branch zeros those inventories and, as in STARTE,
/// deliberately leaves the existing proton-balance inventory unchanged.
pub fn initialize(
    state: *State,
    mode: SaltEquilibriumMode,
    water_capacity_m3: f64,
    concentrations: Concentrations,
    initial_proton_balance_concentration: f64,
) !void {
    if (!std.math.isFinite(water_capacity_m3))
        return error.NonFiniteSurfaceIonComplexWaterCapacity;
    if (water_capacity_m3 < 0)
        return error.InvalidSurfaceIonComplexWaterCapacity;

    switch (mode) {
        .static => state.inventories = std.mem.zeroes(Inventories),
        .dynamic => {
            if (!std.math.isFinite(initial_proton_balance_concentration))
                return error.NonFiniteSurfaceProtonBalanceConcentration;
            if (initial_proton_balance_concentration < 0)
                return error.InvalidSurfaceProtonBalanceConcentration;
            const next: Inventories = .{
                .anions = try validateAndScale(
                    Anions,
                    concentrations.anions,
                    water_capacity_m3,
                ),
                .aluminum = try validateAndScale(
                    AluminumComplexes,
                    concentrations.aluminum,
                    water_capacity_m3,
                ),
                .iron = try validateAndScale(
                    IronComplexes,
                    concentrations.iron,
                    water_capacity_m3,
                ),
                .base_cations = try validateAndScale(
                    BaseCationComplexes,
                    concentrations.base_cations,
                    water_capacity_m3,
                ),
                .phosphate = try validateAndScale(
                    PhosphateComplexes,
                    concentrations.phosphate,
                    water_capacity_m3,
                ),
            };
            const proton_balance_mol =
                initial_proton_balance_concentration * water_capacity_m3;
            if (!std.math.isFinite(proton_balance_mol))
                return error.NonFiniteSurfaceProtonBalanceInventory;
            state.* = .{
                .inventories = next,
                .proton_balance_mol = proton_balance_mol,
            };
        },
    }
}

fn uniformConcentrations(value: f64) Concentrations {
    return .{
        .anions = .{
            .sulfate = value,
            .chloride = value,
            .carbonate = value,
            .bicarbonate = value,
        },
        .aluminum = .{
            .monohydroxide = value,
            .dihydroxide = value,
            .trihydroxide = value,
            .tetrahydroxide = value,
            .sulfate = value,
        },
        .iron = .{
            .monohydroxide = value,
            .dihydroxide = value,
            .trihydroxide = value,
            .tetrahydroxide = value,
            .sulfate = value,
        },
        .base_cations = .{
            .calcium_hydroxide = value,
            .calcium_carbonate = value,
            .calcium_bicarbonate = value,
            .calcium_sulfate = value,
            .magnesium_hydroxide = value,
            .magnesium_carbonate = value,
            .magnesium_bicarbonate = value,
            .magnesium_sulfate = value,
            .sodium_carbonate = value,
            .sodium_sulfate = value,
            .potassium_sulfate = value,
        },
        .phosphate = .{
            .phosphate = value,
            .phosphoric_acid = value,
            .iron_monophosphate = value,
            .iron_diphosphate = value,
            .calcium_phosphate = value,
            .calcium_hydrogen_phosphate = value,
            .calcium_dihydrogen_phosphate = value,
            .magnesium_hydrogen_phosphate = value,
        },
    };
}

test "STARTE dynamic surface ion complexes scale topsoil concentrations" {
    var state: State = undefined;
    try initialize(&state, .dynamic, 2.5, uniformConcentrations(4), 1.0e-3);
    try std.testing.expectEqual(
        @as(f64, 10),
        state.inventories.anions.sulfate,
    );
    try std.testing.expectEqual(
        @as(f64, 10),
        state.inventories.base_cations.potassium_sulfate,
    );
    try std.testing.expectEqual(
        @as(f64, 10),
        state.inventories.phosphate.magnesium_hydrogen_phosphate,
    );
    try std.testing.expectEqual(@as(f64, 2.5e-3), state.proton_balance_mol);
}

test "STARTE static surface salt branch zeros complexes but retains proton balance" {
    var state: State = .{
        .inventories = undefined,
        .proton_balance_mol = 7,
    };
    try initialize(&state, .static, 3, uniformConcentrations(std.math.nan(f64)), -1);
    try std.testing.expectEqual(
        @as(f64, 0),
        state.inventories.aluminum.monohydroxide,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        state.inventories.phosphate.phosphoric_acid,
    );
    try std.testing.expectEqual(@as(f64, 7), state.proton_balance_mol);
}

test "STARTE dynamic surface ion initialization is atomic on invalid input" {
    var state: State = .{
        .inventories = std.mem.zeroes(Inventories),
        .proton_balance_mol = 9,
    };
    var invalid = uniformConcentrations(1);
    invalid.iron.sulfate = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceIonComplexConcentration,
        initialize(&state, .dynamic, 2, invalid, 1.0e-3),
    );
    try std.testing.expectEqual(@as(f64, 9), state.proton_balance_mol);
    try std.testing.expectEqual(
        @as(f64, 0),
        state.inventories.anions.sulfate,
    );
}
