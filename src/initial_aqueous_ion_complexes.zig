const std = @import("std");

pub const PrimaryConcentrations = struct {
    aluminum_mol_per_source_volume: f64, // CAL1
    iron_mol_per_source_volume: f64, // CFE1
    calcium_mol_per_source_volume: f64, // CCA1
    magnesium_mol_per_source_volume: f64, // CMG1
    sodium_mol_per_source_volume: f64, // CNA1
    hydroxide_mol_per_m3: f64, // COH1
    carbonate_mol_per_m3: f64, // CCO31
    bicarbonate_mol_per_m3: f64, // CHCO31
};

/// Runtime dissociation constants use the source STARTE concentration
/// convention; their dimensions balance their corresponding result equation.
pub const DissociationConstants = struct {
    aluminum_hydroxide_1: f64, // DPAL1
    aluminum_hydroxide_2: f64, // DPAL2
    aluminum_hydroxide_3: f64, // DPAL3
    aluminum_hydroxide_4: f64, // DPAL4
    iron_hydroxide_1: f64, // DPFE1
    iron_hydroxide_2: f64, // DPFE2
    iron_hydroxide_3: f64, // DPFE3
    iron_hydroxide_4: f64, // DPFE4
    calcium_hydroxide: f64, // DPCAO
    calcium_carbonate: f64, // DPCAC
    calcium_bicarbonate: f64, // DPCAH
    magnesium_hydroxide: f64, // DPMGO
    magnesium_carbonate: f64, // DPMGC
    magnesium_bicarbonate: f64, // DPMGH
    sodium_carbonate: f64, // DPNAC
};

pub const HydroxideComplexes = struct {
    monohydroxide_mol_per_source_volume: f64,
    dihydroxide_mol_per_source_volume: f64,
    trihydroxide_mol_per_source_volume: f64,
    tetrahydroxide_mol_per_source_volume: f64,
    sulfate_mol_per_source_volume: f64,
};

pub const BaseCationComplexes = struct {
    calcium_hydroxide_mol_per_source_volume: f64,
    calcium_carbonate_mol_per_source_volume: f64,
    calcium_bicarbonate_mol_per_source_volume: f64,
    calcium_sulfate_mol_per_source_volume: f64,
    magnesium_hydroxide_mol_per_source_volume: f64,
    magnesium_carbonate_mol_per_source_volume: f64,
    magnesium_bicarbonate_mol_per_source_volume: f64,
    magnesium_sulfate_mol_per_source_volume: f64,
    sodium_carbonate_mol_per_source_volume: f64,
    sodium_sulfate_mol_per_source_volume: f64,
    potassium_sulfate_mol_per_source_volume: f64,
};

pub const PhosphateComplexSeeds = struct {
    iron_hydrogen_phosphate_mol_p_per_source_volume: f64,
    iron_dihydrogen_phosphate_mol_p_per_source_volume: f64,
    calcium_phosphate_mol_p_per_source_volume: f64,
    calcium_hydrogen_phosphate_mol_p_per_source_volume: f64,
    calcium_dihydrogen_phosphate_mol_p_per_source_volume: f64,
    magnesium_hydrogen_phosphate_mol_p_per_source_volume: f64,
};

pub const Result = struct {
    aluminum: HydroxideComplexes,
    iron: HydroxideComplexes,
    base_cations: BaseCationComplexes,
    phosphate: PhosphateComplexSeeds,
};

/// Direct translation of STARTE.F lines 287--314. The explicit zero values
/// are source initial conditions; subsequent chemistry iterations may form
/// those sulfate and phosphate pairs.
pub fn calculate(primary: PrimaryConcentrations, constants: DissociationConstants) !Result {
    inline for (@typeInfo(PrimaryConcentrations).@"struct".fields) |field| {
        const value = @field(primary, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInitialIonComplexInput;
        if (value < 0) return error.InvalidInitialIonComplexInput;
    }
    inline for (@typeInfo(DissociationConstants).@"struct".fields) |field| {
        const value = @field(constants, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInitialIonComplexConstant;
        if (value <= 0) return error.InvalidInitialIonComplexConstant;
    }

    const aluminum: HydroxideComplexes = .{
        .monohydroxide_mol_per_source_volume = primary.aluminum_mol_per_source_volume * primary.hydroxide_mol_per_m3 / constants.aluminum_hydroxide_1,
        .dihydroxide_mol_per_source_volume = primary.aluminum_mol_per_source_volume * std.math.pow(f64, primary.hydroxide_mol_per_m3, 2.0) / (constants.aluminum_hydroxide_1 * constants.aluminum_hydroxide_2),
        .trihydroxide_mol_per_source_volume = primary.aluminum_mol_per_source_volume * std.math.pow(f64, primary.hydroxide_mol_per_m3, 3.0) / (constants.aluminum_hydroxide_1 * constants.aluminum_hydroxide_2 * constants.aluminum_hydroxide_3),
        .tetrahydroxide_mol_per_source_volume = primary.aluminum_mol_per_source_volume * std.math.pow(f64, primary.hydroxide_mol_per_m3, 4.0) / (constants.aluminum_hydroxide_1 * constants.aluminum_hydroxide_2 * constants.aluminum_hydroxide_3 * constants.aluminum_hydroxide_4),
        .sulfate_mol_per_source_volume = 0.0,
    };
    const iron: HydroxideComplexes = .{
        .monohydroxide_mol_per_source_volume = primary.iron_mol_per_source_volume * primary.hydroxide_mol_per_m3 / constants.iron_hydroxide_1,
        .dihydroxide_mol_per_source_volume = primary.iron_mol_per_source_volume * std.math.pow(f64, primary.hydroxide_mol_per_m3, 2.0) / (constants.iron_hydroxide_1 * constants.iron_hydroxide_2),
        .trihydroxide_mol_per_source_volume = primary.iron_mol_per_source_volume * std.math.pow(f64, primary.hydroxide_mol_per_m3, 3.0) / (constants.iron_hydroxide_1 * constants.iron_hydroxide_2 * constants.iron_hydroxide_3),
        .tetrahydroxide_mol_per_source_volume = primary.iron_mol_per_source_volume * std.math.pow(f64, primary.hydroxide_mol_per_m3, 4.0) / (constants.iron_hydroxide_1 * constants.iron_hydroxide_2 * constants.iron_hydroxide_3 * constants.iron_hydroxide_4),
        .sulfate_mol_per_source_volume = 0.0,
    };
    const base_cations: BaseCationComplexes = .{
        .calcium_hydroxide_mol_per_source_volume = primary.calcium_mol_per_source_volume * primary.hydroxide_mol_per_m3 / constants.calcium_hydroxide,
        .calcium_carbonate_mol_per_source_volume = primary.calcium_mol_per_source_volume * std.math.pow(f64, primary.carbonate_mol_per_m3, 2.0) / constants.calcium_carbonate,
        .calcium_bicarbonate_mol_per_source_volume = primary.calcium_mol_per_source_volume * primary.bicarbonate_mol_per_m3 / constants.calcium_bicarbonate,
        .calcium_sulfate_mol_per_source_volume = 0.0,
        .magnesium_hydroxide_mol_per_source_volume = primary.magnesium_mol_per_source_volume * primary.hydroxide_mol_per_m3 / constants.magnesium_hydroxide,
        .magnesium_carbonate_mol_per_source_volume = primary.magnesium_mol_per_source_volume * std.math.pow(f64, primary.carbonate_mol_per_m3, 2.0) / constants.magnesium_carbonate,
        .magnesium_bicarbonate_mol_per_source_volume = primary.magnesium_mol_per_source_volume * primary.bicarbonate_mol_per_m3 / constants.magnesium_bicarbonate,
        .magnesium_sulfate_mol_per_source_volume = 0.0,
        .sodium_carbonate_mol_per_source_volume = primary.sodium_mol_per_source_volume * primary.carbonate_mol_per_m3 / constants.sodium_carbonate,
        .sodium_sulfate_mol_per_source_volume = 0.0,
        .potassium_sulfate_mol_per_source_volume = 0.0,
    };
    const result: Result = .{
        .aluminum = aluminum,
        .iron = iron,
        .base_cations = base_cations,
        .phosphate = .{
            .iron_hydrogen_phosphate_mol_p_per_source_volume = 0.0,
            .iron_dihydrogen_phosphate_mol_p_per_source_volume = 0.0,
            .calcium_phosphate_mol_p_per_source_volume = 0.0,
            .calcium_hydrogen_phosphate_mol_p_per_source_volume = 0.0,
            .calcium_dihydrogen_phosphate_mol_p_per_source_volume = 0.0,
            .magnesium_hydrogen_phosphate_mol_p_per_source_volume = 0.0,
        },
    };
    try validateResult(result);
    return result;
}

fn validateResult(result: Result) !void {
    inline for (.{ result.aluminum, result.iron }) |group| {
        inline for (@typeInfo(HydroxideComplexes).@"struct".fields) |field|
            if (!std.math.isFinite(@field(group, field.name))) return error.NonFiniteInitialIonComplexResult;
    }
    inline for (@typeInfo(BaseCationComplexes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.base_cations, field.name))) return error.NonFiniteInitialIonComplexResult;
}

fn unitConstants() DissociationConstants {
    var result: DissociationConstants = undefined;
    inline for (@typeInfo(DissociationConstants).@"struct".fields) |field| @field(result, field.name) = 1;
    return result;
}

test "STARTE initial ion complexes preserve equations and explicit zero seeds" {
    const result = try calculate(.{
        .aluminum_mol_per_source_volume = 2,
        .iron_mol_per_source_volume = 3,
        .calcium_mol_per_source_volume = 4,
        .magnesium_mol_per_source_volume = 5,
        .sodium_mol_per_source_volume = 6,
        .hydroxide_mol_per_m3 = 2,
        .carbonate_mol_per_m3 = 3,
        .bicarbonate_mol_per_m3 = 4,
    }, unitConstants());
    try std.testing.expectEqual(@as(f64, 4), result.aluminum.monohydroxide_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 32), result.aluminum.tetrahydroxide_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 48), result.iron.tetrahydroxide_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 36), result.base_cations.calcium_carbonate_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 45), result.base_cations.magnesium_carbonate_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 18), result.base_cations.sodium_carbonate_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0), result.aluminum.sulfate_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0), result.base_cations.potassium_sulfate_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 0), result.phosphate.magnesium_hydrogen_phosphate_mol_p_per_source_volume);
}

test "STARTE calcium and magnesium carbonate retain source squared-carbonate equation" {
    var constants = unitConstants();
    constants.calcium_carbonate = 2;
    constants.magnesium_carbonate = 4;
    const result = try calculate(.{
        .aluminum_mol_per_source_volume = 0,
        .iron_mol_per_source_volume = 0,
        .calcium_mol_per_source_volume = 2,
        .magnesium_mol_per_source_volume = 8,
        .sodium_mol_per_source_volume = 0,
        .hydroxide_mol_per_m3 = 1,
        .carbonate_mol_per_m3 = 3,
        .bicarbonate_mol_per_m3 = 1,
    }, constants);
    try std.testing.expectEqual(@as(f64, 9), result.base_cations.calcium_carbonate_mol_per_source_volume);
    try std.testing.expectEqual(@as(f64, 18), result.base_cations.magnesium_carbonate_mol_per_source_volume);
}

test "STARTE ion complex initialization validates atomically" {
    var constants = unitConstants();
    constants.sodium_carbonate = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteInitialIonComplexConstant, calculate(.{
        .aluminum_mol_per_source_volume = 1,
        .iron_mol_per_source_volume = 1,
        .calcium_mol_per_source_volume = 1,
        .magnesium_mol_per_source_volume = 1,
        .sodium_mol_per_source_volume = 1,
        .hydroxide_mol_per_m3 = 1,
        .carbonate_mol_per_m3 = 1,
        .bicarbonate_mol_per_m3 = 1,
    }, constants));
}
