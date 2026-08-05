const std = @import("std");

pub const Inputs = struct {
    aluminum_source_mol_per_source_volume: f64, // CALZ
    iron_source_mol_per_source_volume: f64, // CFEZ
    hydroxide_mol_per_m3: f64, // COH1
    aluminum_hydroxide_solubility_product: f64, // SPALO
    iron_hydroxide_solubility_product: f64, // SPFEO
};

pub const Result = struct {
    aluminum_mol_per_source_volume: f64, // CAL1
    iron_mol_per_source_volume: f64, // CFE1
};

/// Direct translation of STARTE.F lines 265--274. A negative finite source
/// concentration is an initialization signal, not a physical concentration:
/// it requests the hydroxide-solubility equilibrium value for that metal.
pub fn calculate(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteInitialMetalInput;
    }
    if (inputs.hydroxide_mol_per_m3 <= 0 or
        inputs.aluminum_hydroxide_solubility_product <= 0 or
        inputs.iron_hydroxide_solubility_product <= 0)
        return error.InvalidInitialMetalEquilibriumInput;

    const aluminum_equilibrium = inputs.aluminum_hydroxide_solubility_product /
        std.math.pow(f64, inputs.hydroxide_mol_per_m3, 3.0);
    const aluminum = if (inputs.aluminum_source_mol_per_source_volume < 0)
        aluminum_equilibrium
    else
        @min(inputs.aluminum_source_mol_per_source_volume, aluminum_equilibrium);

    const iron_equilibrium = inputs.iron_hydroxide_solubility_product /
        std.math.pow(f64, inputs.hydroxide_mol_per_m3, 3.0);
    const iron = if (inputs.iron_source_mol_per_source_volume < 0)
        iron_equilibrium
    else
        @min(inputs.iron_source_mol_per_source_volume, iron_equilibrium);

    const result: Result = .{
        .aluminum_mol_per_source_volume = aluminum,
        .iron_mol_per_source_volume = iron,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteInitialMetalResult;
    }
    return result;
}

test "STARTE negative aluminum and iron signals request hydroxide equilibrium" {
    const result = try calculate(.{
        .aluminum_source_mol_per_source_volume = -1,
        .iron_source_mol_per_source_volume = -2,
        .hydroxide_mol_per_m3 = 2,
        .aluminum_hydroxide_solubility_product = 80,
        .iron_hydroxide_solubility_product = 40,
    });
    try std.testing.expectEqual(@as(f64, 10), result.aluminum_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 5), result.iron_mol_per_source_volume);
}

test "STARTE nonnegative metal sources are independently capped at equilibrium" {
    const result = try calculate(.{
        .aluminum_source_mol_per_source_volume = 4,
        .iron_source_mol_per_source_volume = 9,
        .hydroxide_mol_per_m3 = 2,
        .aluminum_hydroxide_solubility_product = 80,
        .iron_hydroxide_solubility_product = 40,
    });
    try std.testing.expectEqual(@as(f64, 4), result.aluminum_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 5), result.iron_mol_per_source_volume);
}

test "STARTE metal equilibrium rejects invalid runtime constants atomically" {
    try std.testing.expectError(error.NonFiniteInitialMetalInput, calculate(.{
        .aluminum_source_mol_per_source_volume = -1,
        .iron_source_mol_per_source_volume = 1,
        .hydroxide_mol_per_m3 = 1,
        .aluminum_hydroxide_solubility_product = 1,
        .iron_hydroxide_solubility_product = std.math.nan(f64),
    }));
    try std.testing.expectError(error.InvalidInitialMetalEquilibriumInput, calculate(.{
        .aluminum_source_mol_per_source_volume = -1,
        .iron_source_mol_per_source_volume = 1,
        .hydroxide_mol_per_m3 = 0,
        .aluminum_hydroxide_solubility_product = 1,
        .iron_hydroxide_solubility_product = 1,
    }));
}
