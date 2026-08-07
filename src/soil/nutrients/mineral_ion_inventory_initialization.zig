const std = @import("std");

pub const Inputs = struct {
    field_capacity_water_m3: f64,
    ammonium_volume_fraction_non_band: f64,
    ammonium_volume_fraction_band: f64,
    nitrate_volume_fraction_non_band: f64,
    nitrate_volume_fraction_band: f64,
    phosphate_volume_fraction_non_band: f64,
    phosphate_volume_fraction_band: f64,
    ammonium_concentration_g_n_per_m3: f64,
    ammonia_concentration_g_n_per_m3: f64,
    nitrate_concentration_g_n_per_m3: f64,
    dihydrogen_phosphate_concentration_g_p_per_m3: f64,
    hydrogen_phosphate_concentration_g_p_per_m3: f64,
    dissolved: DissolvedConcentrations = .{},
};

pub const DissolvedConcentrations = struct {
    aluminum_mol_per_m3: f64 = 0,
    iron_mol_per_m3: f64 = 0,
    hydrogen_mol_per_m3: f64 = 0,
    calcium_mol_per_m3: f64 = 0,
    magnesium_mol_per_m3: f64 = 0,
    sodium_mol_per_m3: f64 = 0,
    potassium_mol_per_m3: f64 = 0,
    hydroxide_mol_per_m3: f64 = 0,
    sulfate_mol_per_m3: f64 = 0,
    chloride_mol_per_m3: f64 = 0,
    carbonate_mol_per_m3: f64 = 0,
    bicarbonate_mol_per_m3: f64 = 0,
    aloh1_mol_per_m3: f64 = 0,
    aloh2_mol_per_m3: f64 = 0,
    aloh3_mol_per_m3: f64 = 0,
    aloh4_mol_per_m3: f64 = 0,
    also4_mol_per_m3: f64 = 0,
    feoh1_mol_per_m3: f64 = 0,
    feoh2_mol_per_m3: f64 = 0,
    feoh3_mol_per_m3: f64 = 0,
    feoh4_mol_per_m3: f64 = 0,
    feso4_mol_per_m3: f64 = 0,
    caoh_mol_per_m3: f64 = 0,
    caco3_mol_per_m3: f64 = 0,
    cahco3_mol_per_m3: f64 = 0,
    caso4_mol_per_m3: f64 = 0,
    mgoh_mol_per_m3: f64 = 0,
    mgco3_mol_per_m3: f64 = 0,
    mghco3_mol_per_m3: f64 = 0,
    mgso4_mol_per_m3: f64 = 0,
    naco3_mol_per_m3: f64 = 0,
    naso4_mol_per_m3: f64 = 0,
    kaso4_mol_per_m3: f64 = 0,
    po4_mol_per_m3: f64 = 0,
    h3po4_mol_per_m3: f64 = 0,
    fehpo4_mol_per_m3: f64 = 0,
    feh2po4_mol_per_m3: f64 = 0,
    capo4_mol_per_m3: f64 = 0,
    cahpo4_mol_per_m3: f64 = 0,
    cah2po4_mol_per_m3: f64 = 0,
    mghpo4_mol_per_m3: f64 = 0,
};

pub const Result = struct {
    ammonium_non_band_g_n: f64,
    ammonium_band_g_n: f64,
    ammonia_non_band_g_n: f64,
    ammonia_band_g_n: f64,
    nitrate_non_band_g_n: f64,
    nitrate_band_g_n: f64,
    nitrite_non_band_g_n: f64,
    nitrite_band_g_n: f64,
    dihydrogen_phosphate_non_band_g_p: f64,
    dihydrogen_phosphate_band_g_p: f64,
    hydrogen_phosphate_non_band_g_p: f64,
    hydrogen_phosphate_band_g_p: f64,
    dissolved_ions: DissolvedInventories,
};

pub const DissolvedInventories = struct {
    aluminum_mol: f64,
    iron_mol: f64,
    hydrogen_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    hydroxide_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
    carbonate_mol: f64,
    bicarbonate_mol: f64,
    aloh1_mol: f64,
    aloh2_mol: f64,
    aloh3_mol: f64,
    aloh4_mol: f64,
    also4_mol: f64,
    feoh1_mol: f64,
    feoh2_mol: f64,
    feoh3_mol: f64,
    feoh4_mol: f64,
    feso4_mol: f64,
    caoh_mol: f64,
    caco3_mol: f64,
    cahco3_mol: f64,
    caso4_mol: f64,
    mgoh_mol: f64,
    mgco3_mol: f64,
    mghco3_mol: f64,
    mgso4_mol: f64,
    naco3_mol: f64,
    naso4_mol: f64,
    kaso4_mol: f64,
    po4_non_band_mol: f64,
    h3po4_non_band_mol: f64,
    fehpo4_non_band_mol: f64,
    feh2po4_non_band_mol: f64,
    capo4_non_band_mol: f64,
    cahpo4_non_band_mol: f64,
    cah2po4_non_band_mol: f64,
    mghpo4_non_band_mol: f64,
    po4_band_mol: f64,
    h3po4_band_mol: f64,
    fehpo4_band_mol: f64,
    feh2po4_band_mol: f64,
    capo4_band_mol: f64,
    cahpo4_band_mol: f64,
    cah2po4_band_mol: f64,
    mghpo4_band_mol: f64,
};

/// Direct translation of `starte.f` lines 1435--1517. This pure kernel initializes
/// layer inventories from concentration and zone fractions. `ZHYSI` line 1518
/// is intentionally owned by `soil_proton_balance_initialization`.
pub fn calculate(inputs: Inputs) !Result {
    inline for (.{
        inputs.field_capacity_water_m3,
        inputs.ammonium_volume_fraction_non_band,
        inputs.ammonium_volume_fraction_band,
        inputs.nitrate_volume_fraction_non_band,
        inputs.nitrate_volume_fraction_band,
        inputs.phosphate_volume_fraction_non_band,
        inputs.phosphate_volume_fraction_band,
        inputs.ammonium_concentration_g_n_per_m3,
        inputs.ammonia_concentration_g_n_per_m3,
        inputs.nitrate_concentration_g_n_per_m3,
        inputs.dihydrogen_phosphate_concentration_g_p_per_m3,
        inputs.hydrogen_phosphate_concentration_g_p_per_m3,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteSoilMineralIonInput;
        if (value < 0) return error.InvalidSoilMineralIonInput;
    }
    inline for (@typeInfo(DissolvedConcentrations).@"struct".fields) |field| {
        const value = @field(inputs.dissolved, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSoilMineralIonInput;
        if (value < 0) return error.InvalidSoilMineralIonInput;
    }

    const water = inputs.field_capacity_water_m3;
    const phosphorus_non_band = inputs.phosphate_volume_fraction_non_band;
    const phosphorus_band = inputs.phosphate_volume_fraction_band;
    const ions = inputs.dissolved;

    const result: Result = .{
        .ammonium_non_band_g_n = inputs.ammonium_concentration_g_n_per_m3 * water * inputs.ammonium_volume_fraction_non_band * 14.0,
        .ammonium_band_g_n = inputs.ammonium_concentration_g_n_per_m3 * water * inputs.ammonium_volume_fraction_band * 14.0,
        .ammonia_non_band_g_n = inputs.ammonia_concentration_g_n_per_m3 * water * inputs.ammonium_volume_fraction_non_band * 14.0,
        .ammonia_band_g_n = inputs.ammonia_concentration_g_n_per_m3 * water * inputs.ammonium_volume_fraction_band * 14.0,
        .nitrate_non_band_g_n = inputs.nitrate_concentration_g_n_per_m3 * water * inputs.nitrate_volume_fraction_non_band * 14.0,
        .nitrate_band_g_n = inputs.nitrate_concentration_g_n_per_m3 * water * inputs.nitrate_volume_fraction_band * 14.0,
        .nitrite_non_band_g_n = 0.0,
        .nitrite_band_g_n = 0.0,
        .dihydrogen_phosphate_non_band_g_p = inputs.dihydrogen_phosphate_concentration_g_p_per_m3 * water * phosphorus_non_band * 31.0,
        .dihydrogen_phosphate_band_g_p = inputs.dihydrogen_phosphate_concentration_g_p_per_m3 * water * phosphorus_band * 31.0,
        .hydrogen_phosphate_non_band_g_p = inputs.hydrogen_phosphate_concentration_g_p_per_m3 * water * phosphorus_non_band * 31.0,
        .hydrogen_phosphate_band_g_p = inputs.hydrogen_phosphate_concentration_g_p_per_m3 * water * phosphorus_band * 31.0,
        .dissolved_ions = .{
            .aluminum_mol = ions.aluminum_mol_per_m3 * water,
            .iron_mol = ions.iron_mol_per_m3 * water,
            .hydrogen_mol = ions.hydrogen_mol_per_m3 * water,
            .calcium_mol = ions.calcium_mol_per_m3 * water,
            .magnesium_mol = ions.magnesium_mol_per_m3 * water,
            .sodium_mol = ions.sodium_mol_per_m3 * water,
            .potassium_mol = ions.potassium_mol_per_m3 * water,
            .hydroxide_mol = ions.hydroxide_mol_per_m3 * water,
            .sulfate_mol = ions.sulfate_mol_per_m3 * water,
            .chloride_mol = ions.chloride_mol_per_m3 * water,
            .carbonate_mol = ions.carbonate_mol_per_m3 * water,
            .bicarbonate_mol = ions.bicarbonate_mol_per_m3 * water,
            .aloh1_mol = ions.aloh1_mol_per_m3 * water,
            .aloh2_mol = ions.aloh2_mol_per_m3 * water,
            .aloh3_mol = ions.aloh3_mol_per_m3 * water,
            .aloh4_mol = ions.aloh4_mol_per_m3 * water,
            .also4_mol = ions.also4_mol_per_m3 * water,
            .feoh1_mol = ions.feoh1_mol_per_m3 * water,
            .feoh2_mol = ions.feoh2_mol_per_m3 * water,
            .feoh3_mol = ions.feoh3_mol_per_m3 * water,
            .feoh4_mol = ions.feoh4_mol_per_m3 * water,
            .feso4_mol = ions.feso4_mol_per_m3 * water,
            .caoh_mol = ions.caoh_mol_per_m3 * water,
            .caco3_mol = ions.caco3_mol_per_m3 * water,
            .cahco3_mol = ions.cahco3_mol_per_m3 * water,
            .caso4_mol = ions.caso4_mol_per_m3 * water,
            .mgoh_mol = ions.mgoh_mol_per_m3 * water,
            .mgco3_mol = ions.mgco3_mol_per_m3 * water,
            .mghco3_mol = ions.mghco3_mol_per_m3 * water,
            .mgso4_mol = ions.mgso4_mol_per_m3 * water,
            .naco3_mol = ions.naco3_mol_per_m3 * water,
            .naso4_mol = ions.naso4_mol_per_m3 * water,
            .kaso4_mol = ions.kaso4_mol_per_m3 * water,
            .po4_non_band_mol = ions.po4_mol_per_m3 * water * phosphorus_non_band,
            .h3po4_non_band_mol = ions.h3po4_mol_per_m3 * water * phosphorus_non_band,
            .fehpo4_non_band_mol = ions.fehpo4_mol_per_m3 * water * phosphorus_non_band,
            .feh2po4_non_band_mol = ions.feh2po4_mol_per_m3 * water * phosphorus_non_band,
            .capo4_non_band_mol = ions.capo4_mol_per_m3 * water * phosphorus_non_band,
            .cahpo4_non_band_mol = ions.cahpo4_mol_per_m3 * water * phosphorus_non_band,
            .cah2po4_non_band_mol = ions.cah2po4_mol_per_m3 * water * phosphorus_non_band,
            .mghpo4_non_band_mol = ions.mghpo4_mol_per_m3 * water * phosphorus_non_band,
            .po4_band_mol = ions.po4_mol_per_m3 * water * phosphorus_band,
            .h3po4_band_mol = ions.h3po4_mol_per_m3 * water * phosphorus_band,
            .fehpo4_band_mol = ions.fehpo4_mol_per_m3 * water * phosphorus_band,
            .feh2po4_band_mol = ions.feh2po4_mol_per_m3 * water * phosphorus_band,
            .capo4_band_mol = ions.capo4_mol_per_m3 * water * phosphorus_band,
            .cahpo4_band_mol = ions.cahpo4_mol_per_m3 * water * phosphorus_band,
            .cah2po4_band_mol = ions.cah2po4_mol_per_m3 * water * phosphorus_band,
            .mghpo4_band_mol = ions.mghpo4_mol_per_m3 * water * phosphorus_band,
        },
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSoilMineralIonResult;
    inline for (@typeInfo(DissolvedInventories).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.dissolved_ions, field.name)))
            return error.NonFiniteSoilMineralIonResult;
    return result;
}

fn baseInputs() Inputs {
    return .{
        .field_capacity_water_m3 = 2.0,
        .ammonium_volume_fraction_non_band = 0.75,
        .ammonium_volume_fraction_band = 0.25,
        .nitrate_volume_fraction_non_band = 0.8,
        .nitrate_volume_fraction_band = 0.2,
        .phosphate_volume_fraction_non_band = 0.7,
        .phosphate_volume_fraction_band = 0.3,
        .ammonium_concentration_g_n_per_m3 = 1.5,
        .ammonia_concentration_g_n_per_m3 = 0.2,
        .nitrate_concentration_g_n_per_m3 = 1.0,
        .dihydrogen_phosphate_concentration_g_p_per_m3 = 0.5,
        .hydrogen_phosphate_concentration_g_p_per_m3 = 0.25,
        .dissolved = .{
            .aluminum_mol_per_m3 = 1.0,
            .sulfate_mol_per_m3 = 2.0,
            .h3po4_mol_per_m3 = 0.4,
            .cah2po4_mol_per_m3 = 0.6,
        },
    };
}

test "STARTE mineral and ion inventories preserve source initialization equations" {
    const result = try calculate(baseInputs());
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.5 * 2.0 * 0.75 * 14.0),
        result.ammonium_non_band_g_n,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.5 * 2.0 * 0.3 * 31.0),
        result.dihydrogen_phosphate_band_g_p,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), result.nitrite_non_band_g_n);
    try std.testing.expectApproxEqAbs(
        @as(f64, 2.0),
        result.dissolved_ions.aluminum_mol,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.4 * 2.0 * 0.7),
        result.dissolved_ions.h3po4_non_band_mol,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.6 * 2.0 * 0.3),
        result.dissolved_ions.cah2po4_band_mol,
        1.0e-15,
    );
}

test "STARTE runtime layer traversal supports arbitrary depth counts" {
    const inputs = [_]Inputs{
        baseInputs(),
        .{ .field_capacity_water_m3 = 0.0, .ammonium_volume_fraction_non_band = 0, .ammonium_volume_fraction_band = 1, .nitrate_volume_fraction_non_band = 1, .nitrate_volume_fraction_band = 0, .phosphate_volume_fraction_non_band = 1, .phosphate_volume_fraction_band = 0, .ammonium_concentration_g_n_per_m3 = 2, .ammonia_concentration_g_n_per_m3 = 1, .nitrate_concentration_g_n_per_m3 = 3, .dihydrogen_phosphate_concentration_g_p_per_m3 = 4, .hydrogen_phosphate_concentration_g_p_per_m3 = 5 },
        .{ .field_capacity_water_m3 = 100.0, .ammonium_volume_fraction_non_band = 0.5, .ammonium_volume_fraction_band = 0.5, .nitrate_volume_fraction_non_band = 0.5, .nitrate_volume_fraction_band = 0.5, .phosphate_volume_fraction_non_band = 0.5, .phosphate_volume_fraction_band = 0.5, .ammonium_concentration_g_n_per_m3 = 0.01, .ammonia_concentration_g_n_per_m3 = 0.02, .nitrate_concentration_g_n_per_m3 = 0.03, .dihydrogen_phosphate_concentration_g_p_per_m3 = 0.04, .hydrogen_phosphate_concentration_g_p_per_m3 = 0.05 },
    };
    for (inputs) |input| {
        const result = try calculate(input);
        try std.testing.expect(std.math.isFinite(result.ammonium_non_band_g_n));
        try std.testing.expect(std.math.isFinite(result.dissolved_ions.cah2po4_band_mol));
    }
}

test "STARTE inventory initialization rejects late invalid values" {
    var inputs = baseInputs();
    inputs.dissolved.kaso4_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSoilMineralIonInput,
        calculate(inputs),
    );
}
