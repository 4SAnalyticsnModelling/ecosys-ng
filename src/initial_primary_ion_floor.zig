const std = @import("std");

pub const Concentrations = struct {
    ammonium_mol_n_per_source_volume: f64, // CN41
    ammonia_mol_n_per_source_volume: f64, // CN31
    nitrate_mol_n_per_source_volume: f64, // CNO1/CNOZ
    phosphate_mol_p_per_source_volume: f64, // CPO1/CPOZ
    aluminum_mol_per_source_volume: f64, // CAL1
    iron_mol_per_source_volume: f64, // CFE1
    calcium_mol_per_source_volume: f64, // CCA1/CCAZ
    magnesium_mol_per_source_volume: f64, // CMG1/CMGZ
    sodium_mol_per_source_volume: f64, // CNA1/CNAZ
    potassium_mol_per_source_volume: f64, // CKA1/CKAZ
    sulfate_mol_s_per_source_volume: f64, // CSO41/CSOZ
    chloride_mol_per_source_volume: f64, // CCL1/CCLZ
};

/// Direct translation of STARTE.F lines 275--286. The concentration basis is
/// selected by the enclosing runtime source: `m3` for precipitation and
/// irrigation, and the active soil source-volume basis for soil initialization.
pub fn apply(inputs: Concentrations, minimum_concentration: f64) !Concentrations {
    if (!std.math.isFinite(minimum_concentration) or minimum_concentration < 0)
        return error.InvalidInitialIonFloor;
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteInitialIonConcentration;
    }

    const result: Concentrations = .{
        .ammonium_mol_n_per_source_volume = @max(minimum_concentration, inputs.ammonium_mol_n_per_source_volume),
        .ammonia_mol_n_per_source_volume = @max(minimum_concentration, inputs.ammonia_mol_n_per_source_volume),
        .nitrate_mol_n_per_source_volume = @max(minimum_concentration, inputs.nitrate_mol_n_per_source_volume),
        .phosphate_mol_p_per_source_volume = @max(minimum_concentration, inputs.phosphate_mol_p_per_source_volume),
        .aluminum_mol_per_source_volume = @max(minimum_concentration, inputs.aluminum_mol_per_source_volume),
        .iron_mol_per_source_volume = @max(minimum_concentration, inputs.iron_mol_per_source_volume),
        .calcium_mol_per_source_volume = @max(minimum_concentration, inputs.calcium_mol_per_source_volume),
        .magnesium_mol_per_source_volume = @max(minimum_concentration, inputs.magnesium_mol_per_source_volume),
        .sodium_mol_per_source_volume = @max(minimum_concentration, inputs.sodium_mol_per_source_volume),
        .potassium_mol_per_source_volume = @max(minimum_concentration, inputs.potassium_mol_per_source_volume),
        .sulfate_mol_s_per_source_volume = @max(minimum_concentration, inputs.sulfate_mol_s_per_source_volume),
        .chloride_mol_per_source_volume = @max(minimum_concentration, inputs.chloride_mol_per_source_volume),
    };
    return result;
}

fn ascendingInputs() Concentrations {
    return .{
        .ammonium_mol_n_per_source_volume = -6,
        .ammonia_mol_n_per_source_volume = -5,
        .nitrate_mol_n_per_source_volume = -4,
        .phosphate_mol_p_per_source_volume = -3,
        .aluminum_mol_per_source_volume = -2,
        .iron_mol_per_source_volume = -1,
        .calcium_mol_per_source_volume = 0,
        .magnesium_mol_per_source_volume = 1,
        .sodium_mol_per_source_volume = 2,
        .potassium_mol_per_source_volume = 3,
        .sulfate_mol_s_per_source_volume = 4,
        .chloride_mol_per_source_volume = 5,
    };
}

test "STARTE primary ion floors preserve exact field order and runtime minimum" {
    const result = try apply(ascendingInputs(), 0.5);
    try std.testing.expectEqual(@as(f64, 0.5), result.ammonium_mol_n_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0.5), result.ammonia_mol_n_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0.5), result.nitrate_mol_n_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0.5), result.phosphate_mol_p_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0.5), result.aluminum_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0.5), result.iron_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0.5), result.calcium_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 1), result.magnesium_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 2), result.sodium_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 3), result.potassium_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 4), result.sulfate_mol_s_per_source_volume);
    try std.testing.expectEqual(@as(f64, 5), result.chloride_mol_per_source_volume);
}

test "STARTE primary ion floors accept source-scale tiny runtime minimum" {
    const result = try apply(ascendingInputs(), 1.0e-48);
    inline for (.{
        result.ammonium_mol_n_per_source_volume,
        result.ammonia_mol_n_per_source_volume,
        result.nitrate_mol_n_per_source_volume,
        result.phosphate_mol_p_per_source_volume,
        result.aluminum_mol_per_source_volume,
        result.iron_mol_per_source_volume,
        result.calcium_mol_per_source_volume,
    }) |value| try std.testing.expectEqual(@as(f64, 1.0e-48), value);
    try std.testing.expectEqual(@as(f64, 1), result.magnesium_mol_per_source_volume);
}

test "STARTE primary ion floors reject late invalid input atomically" {
    var invalid = ascendingInputs();
    invalid.chloride_mol_per_source_volume = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteInitialIonConcentration, apply(invalid, 1.0e-48));
    try std.testing.expectError(error.InvalidInitialIonFloor, apply(ascendingInputs(), -1));
}
