const std = @import("std");

pub const Inputs = struct {
    soil_mass_megagram: f64,
    ammonium_volume_fraction_non_band: f64,
    ammonium_volume_fraction_band: f64,
    phosphate_volume_fraction_non_band: f64,
    phosphate_volume_fraction_band: f64,
    exchangeable_parameters: ExchangeableParameters,
};

pub const ExchangeableParameters = struct {
    ammonium_mol_per_megagram: f64, // XN41
    hydrogen_mol_per_megagram: f64, // XHY1
    aluminum_mol_per_megagram: f64, // XAL1
    iron_mol_per_megagram: f64, // XFE1
    calcium_mol_per_megagram: f64, // XCA1
    magnesium_mol_per_megagram: f64, // XMG1
    sodium_mol_per_megagram: f64, // XNA1
    potassium_mol_per_megagram: f64, // XKA1
    bicarbonate_mol_per_megagram: f64, // XHC1
    phosphate_hydroxide0_mol_per_megagram: f64, // XOH01
    phosphate_hydroxide1_mol_per_megagram: f64, // XOH11
    phosphate_hydroxide2_mol_per_megagram: f64, // XOH21
    phosphate_hpo4_mol_per_megagram: f64, // XH1P1
    phosphate_h2po4_mol_per_megagram: f64, // XH2P1
};

pub const State = struct {
    macropore: MacroporeInventory,
    exchangeable: ExchangeableInventory,
};

pub const MacroporeInventory = struct {
    carbon_dioxide_g: f64,
    methane_g: f64,
    oxygen_g: f64,
    nitrogen_g: f64,
    nitrous_oxide_g: f64,
    hydrogen_g: f64,
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

pub const ExchangeableInventory = struct {
    ammonium_non_band_mol: f64,
    ammonium_band_mol: f64,
    hydrogen_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    bicarbonate_mol: f64,
    hydroxide0_non_band_mol: f64,
    hydroxide1_non_band_mol: f64,
    hydroxide2_non_band_mol: f64,
    hpo4_non_band_mol: f64,
    h2po4_non_band_mol: f64,
    hydroxide0_band_mol: f64,
    hydroxide1_band_mol: f64,
    hydroxide2_band_mol: f64,
    hpo4_band_mol: f64,
    h2po4_band_mol: f64,
};

/// Direct translation of `starte.f` lines 1533--1632. Initializes all macropore
/// inventories to zero and applies exchangeable pool mass mappings from runtime
/// soil mass and band fractions.
pub fn apply(state: *State, inputs: Inputs) !void {
    inline for (.{ inputs.soil_mass_megagram, inputs.ammonium_volume_fraction_non_band, inputs.ammonium_volume_fraction_band, inputs.phosphate_volume_fraction_non_band, inputs.phosphate_volume_fraction_band }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteMacroporeExchangeableInput;
        if (value < 0) return error.InvalidMacroporeExchangeableInput;
    }
    inline for (@typeInfo(ExchangeableParameters).@"struct".fields) |field| {
        const value = @field(inputs.exchangeable_parameters, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteMacroporeExchangeableInput;
        if (value < 0) return error.InvalidMacroporeExchangeableInput;
    }

    const soil_mass = inputs.soil_mass_megagram;
    const ammonium_non_band = inputs.ammonium_volume_fraction_non_band;
    const ammonium_band = inputs.ammonium_volume_fraction_band;
    const phosphate_non_band = inputs.phosphate_volume_fraction_non_band;
    const phosphate_band = inputs.phosphate_volume_fraction_band;
    const p = inputs.exchangeable_parameters;
    const exchangeable: ExchangeableInventory = .{
        .ammonium_non_band_mol = p.ammonium_mol_per_megagram * soil_mass * ammonium_non_band,
        .ammonium_band_mol = p.ammonium_mol_per_megagram * soil_mass * ammonium_band,
        .hydrogen_mol = p.hydrogen_mol_per_megagram * soil_mass,
        .aluminum_mol = p.aluminum_mol_per_megagram * soil_mass,
        .iron_mol = p.iron_mol_per_megagram * soil_mass,
        .calcium_mol = p.calcium_mol_per_megagram * soil_mass,
        .magnesium_mol = p.magnesium_mol_per_megagram * soil_mass,
        .sodium_mol = p.sodium_mol_per_megagram * soil_mass,
        .potassium_mol = p.potassium_mol_per_megagram * soil_mass,
        .bicarbonate_mol = p.bicarbonate_mol_per_megagram * soil_mass,
        .hydroxide0_non_band_mol = p.phosphate_hydroxide0_mol_per_megagram * soil_mass * phosphate_non_band,
        .hydroxide1_non_band_mol = p.phosphate_hydroxide1_mol_per_megagram * soil_mass * phosphate_non_band,
        .hydroxide2_non_band_mol = p.phosphate_hydroxide2_mol_per_megagram * soil_mass * phosphate_non_band,
        .hpo4_non_band_mol = p.phosphate_hpo4_mol_per_megagram * soil_mass * phosphate_non_band,
        .h2po4_non_band_mol = p.phosphate_h2po4_mol_per_megagram * soil_mass * phosphate_non_band,
        .hydroxide0_band_mol = p.phosphate_hydroxide0_mol_per_megagram * soil_mass * phosphate_band,
        .hydroxide1_band_mol = p.phosphate_hydroxide1_mol_per_megagram * soil_mass * phosphate_band,
        .hydroxide2_band_mol = p.phosphate_hydroxide2_mol_per_megagram * soil_mass * phosphate_band,
        .hpo4_band_mol = p.phosphate_hpo4_mol_per_megagram * soil_mass * phosphate_band,
        .h2po4_band_mol = p.phosphate_h2po4_mol_per_megagram * soil_mass * phosphate_band,
    };
    inline for (@typeInfo(ExchangeableInventory).@"struct".fields) |field|
        if (!std.math.isFinite(@field(exchangeable, field.name)))
            return error.NonFiniteMacroporeExchangeableResult;

    // STARTE source order, including repeated H1PO4H assignment.
    state.macropore.carbon_dioxide_g = 0.0;
    state.macropore.methane_g = 0.0;
    state.macropore.oxygen_g = 0.0;
    state.macropore.nitrogen_g = 0.0;
    state.macropore.nitrous_oxide_g = 0.0;
    state.macropore.hydrogen_g = 0.0;
    state.macropore.ammonium_non_band_g_n = 0.0;
    state.macropore.ammonium_band_g_n = 0.0;
    state.macropore.ammonia_non_band_g_n = 0.0;
    state.macropore.ammonia_band_g_n = 0.0;
    state.macropore.nitrate_non_band_g_n = 0.0;
    state.macropore.nitrate_band_g_n = 0.0;
    state.macropore.nitrite_non_band_g_n = 0.0;
    state.macropore.nitrite_band_g_n = 0.0;
    state.macropore.dihydrogen_phosphate_non_band_g_p = 0.0;
    state.macropore.dihydrogen_phosphate_band_g_p = 0.0;
    state.macropore.hydrogen_phosphate_non_band_g_p = 0.0;
    state.macropore.hydrogen_phosphate_band_g_p = 0.0;
    state.macropore.aluminum_mol = 0.0;
    state.macropore.iron_mol = 0.0;
    state.macropore.hydrogen_mol = 0.0;
    state.macropore.calcium_mol = 0.0;
    state.macropore.magnesium_mol = 0.0;
    state.macropore.sodium_mol = 0.0;
    state.macropore.potassium_mol = 0.0;
    state.macropore.hydroxide_mol = 0.0;
    state.macropore.sulfate_mol = 0.0;
    state.macropore.chloride_mol = 0.0;
    state.macropore.carbonate_mol = 0.0;
    state.macropore.bicarbonate_mol = 0.0;
    state.macropore.aloh1_mol = 0.0;
    state.macropore.aloh2_mol = 0.0;
    state.macropore.aloh3_mol = 0.0;
    state.macropore.aloh4_mol = 0.0;
    state.macropore.also4_mol = 0.0;
    state.macropore.feoh1_mol = 0.0;
    state.macropore.feoh2_mol = 0.0;
    state.macropore.feoh3_mol = 0.0;
    state.macropore.feoh4_mol = 0.0;
    state.macropore.feso4_mol = 0.0;
    state.macropore.caoh_mol = 0.0;
    state.macropore.caco3_mol = 0.0;
    state.macropore.cahco3_mol = 0.0;
    state.macropore.caso4_mol = 0.0;
    state.macropore.mgoh_mol = 0.0;
    state.macropore.mgco3_mol = 0.0;
    state.macropore.mghco3_mol = 0.0;
    state.macropore.mgso4_mol = 0.0;
    state.macropore.naco3_mol = 0.0;
    state.macropore.naso4_mol = 0.0;
    state.macropore.kaso4_mol = 0.0;
    state.macropore.po4_non_band_mol = 0.0;
    state.macropore.h3po4_non_band_mol = 0.0;
    state.macropore.fehpo4_non_band_mol = 0.0;
    state.macropore.feh2po4_non_band_mol = 0.0;
    state.macropore.capo4_non_band_mol = 0.0;
    state.macropore.cahpo4_non_band_mol = 0.0;
    state.macropore.cah2po4_non_band_mol = 0.0;
    state.macropore.mghpo4_non_band_mol = 0.0;
    state.macropore.hydrogen_phosphate_non_band_g_p = 0.0;
    state.macropore.po4_band_mol = 0.0;
    state.macropore.hydrogen_phosphate_band_g_p = 0.0;
    state.macropore.h3po4_band_mol = 0.0;
    state.macropore.fehpo4_band_mol = 0.0;
    state.macropore.feh2po4_band_mol = 0.0;
    state.macropore.capo4_band_mol = 0.0;
    state.macropore.cahpo4_band_mol = 0.0;
    state.macropore.cah2po4_band_mol = 0.0;
    state.macropore.mghpo4_band_mol = 0.0;
    state.exchangeable = exchangeable;
}

fn baseState(value: f64) State {
    var state: State = undefined;
    inline for (@typeInfo(MacroporeInventory).@"struct".fields) |field|
        @field(state.macropore, field.name) = value;
    inline for (@typeInfo(ExchangeableInventory).@"struct".fields) |field|
        @field(state.exchangeable, field.name) = value;
    return state;
}

fn sourceLikeInputs() Inputs {
    return .{
        .soil_mass_megagram = 2.0,
        .ammonium_volume_fraction_non_band = 0.75,
        .ammonium_volume_fraction_band = 0.25,
        .phosphate_volume_fraction_non_band = 0.7,
        .phosphate_volume_fraction_band = 0.3,
        .exchangeable_parameters = .{
            .ammonium_mol_per_megagram = 1.0,
            .hydrogen_mol_per_megagram = 2.0,
            .aluminum_mol_per_megagram = 3.0,
            .iron_mol_per_megagram = 4.0,
            .calcium_mol_per_megagram = 5.0,
            .magnesium_mol_per_megagram = 6.0,
            .sodium_mol_per_megagram = 7.0,
            .potassium_mol_per_megagram = 8.0,
            .bicarbonate_mol_per_megagram = 9.0,
            .phosphate_hydroxide0_mol_per_megagram = 10.0,
            .phosphate_hydroxide1_mol_per_megagram = 11.0,
            .phosphate_hydroxide2_mol_per_megagram = 12.0,
            .phosphate_hpo4_mol_per_megagram = 13.0,
            .phosphate_h2po4_mol_per_megagram = 14.0,
        },
    };
}

test "STARTE macropore inventory reset and exchangeable initialization preserve source equations" {
    var state = baseState(99.0);
    try apply(&state, sourceLikeInputs());
    inline for (@typeInfo(MacroporeInventory).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(state.macropore, field.name));
    try std.testing.expectEqual(@as(f64, 1.0 * 2.0 * 0.75), state.exchangeable.ammonium_non_band_mol);
    try std.testing.expectEqual(@as(f64, 1.0 * 2.0 * 0.25), state.exchangeable.ammonium_band_mol);
    try std.testing.expectEqual(@as(f64, 10.0 * 2.0 * 0.3), state.exchangeable.hydroxide0_band_mol);
    try std.testing.expectEqual(@as(f64, 13.0 * 2.0 * 0.7), state.exchangeable.hpo4_non_band_mol);
}

test "STARTE supports runtime layer traversal by repeated per-layer apply" {
    var layers = [_]State{ baseState(1.0), baseState(2.0), baseState(3.0), baseState(4.0) };
    for (&layers) |*layer| try apply(layer, sourceLikeInputs());
    for (layers) |layer| {
        try std.testing.expectEqual(@as(f64, 0), layer.macropore.methane_g);
        try std.testing.expectEqual(@as(f64, 18.0), layer.exchangeable.bicarbonate_mol);
    }
}

test "invalid late parameter is rejected before state mutation" {
    var state = baseState(7.0);
    var inputs = sourceLikeInputs();
    inputs.exchangeable_parameters.phosphate_h2po4_mol_per_megagram = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteMacroporeExchangeableInput,
        apply(&state, inputs),
    );
    try std.testing.expectEqual(@as(f64, 7.0), state.macropore.caco3_mol);
    try std.testing.expectEqual(@as(f64, 7.0), state.exchangeable.h2po4_band_mol);
}
