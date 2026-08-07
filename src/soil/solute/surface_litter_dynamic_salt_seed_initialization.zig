const std = @import("std");

pub const TopLayerConcentrations = struct {
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    chloride_mol_per_m3: f64,
    carbonate_mol_per_m3: f64,
    bicarbonate_mol_per_m3: f64,
    aloh1_mol_per_m3: f64,
    aloh2_mol_per_m3: f64,
    aloh3_mol_per_m3: f64,
    aloh4_mol_per_m3: f64,
    also4_mol_per_m3: f64,
    feoh1_mol_per_m3: f64,
    feoh2_mol_per_m3: f64,
    feoh3_mol_per_m3: f64,
    feoh4_mol_per_m3: f64,
    feso4_mol_per_m3: f64,
    caoh_mol_per_m3: f64,
    caco3_mol_per_m3: f64,
    cahco3_mol_per_m3: f64,
    caso4_mol_per_m3: f64,
    mgoh_mol_per_m3: f64,
    mgco3_mol_per_m3: f64,
    mghco3_mol_per_m3: f64,
    mgso4_mol_per_m3: f64,
    naco3_mol_per_m3: f64,
    naso4_mol_per_m3: f64,
    kaso4_mol_per_m3: f64,
    po4_mol_per_m3: f64,
    h3po4_mol_per_m3: f64,
    fehpo4_mol_per_m3: f64,
    feh2po4_mol_per_m3: f64,
    capo4_mol_per_m3: f64,
    cahpo4_mol_per_m3: f64,
    cah2po4_mol_per_m3: f64,
    mghpo4_mol_per_m3: f64,
};

pub const Inputs = struct {
    dynamic_salt_enabled: bool, // ISALTG != 0
    surface_litter_ph: f64, // PH(0)
    surface_organic_carbon_g: f64, // ORGC(0)
    surface_soil_mass_megagram: f64, // BKVL(0)
    negligible_surface_soil_mass_megagram: f64, // ZEROS
    carboxyl_site_density_mol_per_kg_c: f64, // COOH
    water_ion_product_mol2_per_m6: f64, // DPH2O
    carboxyl_hydrogen_half_saturation_mol_per_m3: f64, // DPCOH
    aluminum_hydroxide_solubility_product: f64, // SPALO
    iron_hydroxide_solubility_product: f64, // SPFEO
    top_layer_ammonium_g_n_per_m3: f64, // CNH4(NU)
    top_layer: TopLayerConcentrations,
};

pub const State = struct {
    hydrogen_activity_mol_per_m3: f64, // CHY1
    hydroxide_activity_mol_per_m3: f64, // COH1
    total_carboxyl_exchange_sites_mol_per_megagram: f64, // XCOOH
    occupied_carboxyl_hydrogen_sites_mol_per_megagram: f64, // XHC1
    effective_cation_exchange_capacity_mol_per_megagram: f64, // CCEC0
    ammonium_mol_n_per_m3: f64, // CN41
    dissolved_aluminum_mol_per_m3: f64, // CAL1
    dissolved_iron_mol_per_m3: f64, // CFE1
    top_layer: TopLayerConcentrations,
};

/// Direct translation of `starte.f` lines 1726--1777.
/// Caller provides the already-selected top active layer concentrations.
pub fn initialize(state: *State, inputs: Inputs) !bool {
    if (!inputs.dynamic_salt_enabled) return false;
    try validateInputs(inputs);
    const candidate = try buildCandidate(inputs);
    state.* = candidate;
    return true;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.surface_litter_ph,
        inputs.surface_organic_carbon_g,
        inputs.surface_soil_mass_megagram,
        inputs.negligible_surface_soil_mass_megagram,
        inputs.carboxyl_site_density_mol_per_kg_c,
        inputs.water_ion_product_mol2_per_m6,
        inputs.carboxyl_hydrogen_half_saturation_mol_per_m3,
        inputs.aluminum_hydroxide_solubility_product,
        inputs.iron_hydroxide_solubility_product,
        inputs.top_layer_ammonium_g_n_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterDynamicSaltSeedInput;
    }
    inline for (@typeInfo(TopLayerConcentrations).@"struct".fields) |field| {
        const value = @field(inputs.top_layer, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterDynamicSaltSeedInput;
    }
    if (inputs.surface_organic_carbon_g < 0 or
        inputs.surface_soil_mass_megagram < 0 or
        inputs.negligible_surface_soil_mass_megagram < 0 or
        inputs.carboxyl_site_density_mol_per_kg_c < 0 or
        inputs.water_ion_product_mol2_per_m6 <= 0 or
        inputs.carboxyl_hydrogen_half_saturation_mol_per_m3 <= 0 or
        inputs.aluminum_hydroxide_solubility_product <= 0 or
        inputs.iron_hydroxide_solubility_product <= 0 or
        inputs.top_layer_ammonium_g_n_per_m3 < 0)
        return error.InvalidSurfaceLitterDynamicSaltSeedInput;
}

fn buildCandidate(inputs: Inputs) !State {
    const chy1 = std.math.pow(f64, 10.0, -(inputs.surface_litter_ph - 3.0));
    const coh1 = inputs.water_ion_product_mol2_per_m6 / chy1;
    const xcooh = if (inputs.surface_soil_mass_megagram >
        inputs.negligible_surface_soil_mass_megagram)
        @max(
            0.0,
            inputs.carboxyl_site_density_mol_per_kg_c *
                1.0e-06 *
                inputs.surface_organic_carbon_g /
                inputs.surface_soil_mass_megagram,
        )
    else
        0.0;
    const xhc1 = xcooh * @min(
        1.0,
        chy1 / inputs.carboxyl_hydrogen_half_saturation_mol_per_m3,
    );
    const candidate: State = .{
        .hydrogen_activity_mol_per_m3 = chy1,
        .hydroxide_activity_mol_per_m3 = coh1,
        .total_carboxyl_exchange_sites_mol_per_megagram = xcooh,
        .occupied_carboxyl_hydrogen_sites_mol_per_megagram = xhc1,
        .effective_cation_exchange_capacity_mol_per_megagram = xcooh,
        .ammonium_mol_n_per_m3 = inputs.top_layer_ammonium_g_n_per_m3 * 0.10 / 14.0,
        .dissolved_aluminum_mol_per_m3 = inputs.aluminum_hydroxide_solubility_product / std.math.pow(f64, coh1, 3.0),
        .dissolved_iron_mol_per_m3 = inputs.iron_hydroxide_solubility_product / std.math.pow(f64, coh1, 3.0),
        .top_layer = inputs.top_layer,
    };
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(candidate, field.name)))
            return error.NonFiniteSurfaceLitterDynamicSaltSeedResult;
    }
    return candidate;
}

fn testInputs() Inputs {
    return .{
        .dynamic_salt_enabled = true,
        .surface_litter_ph = 6.0,
        .surface_organic_carbon_g = 200.0,
        .surface_soil_mass_megagram = 1.0,
        .negligible_surface_soil_mass_megagram = 1.0e-12,
        .carboxyl_site_density_mol_per_kg_c = 200.0,
        .water_ion_product_mol2_per_m6 = 1.0e-14,
        .carboxyl_hydrogen_half_saturation_mol_per_m3 = 1.0e-5,
        .aluminum_hydroxide_solubility_product = 1.0e-33,
        .iron_hydroxide_solubility_product = 2.0e-33,
        .top_layer_ammonium_g_n_per_m3 = 14.0,
        .top_layer = .{
            .calcium_mol_per_m3 = 1.0,
            .magnesium_mol_per_m3 = 2.0,
            .sodium_mol_per_m3 = 3.0,
            .potassium_mol_per_m3 = 4.0,
            .sulfate_mol_per_m3 = 5.0,
            .chloride_mol_per_m3 = 6.0,
            .carbonate_mol_per_m3 = 7.0,
            .bicarbonate_mol_per_m3 = 8.0,
            .aloh1_mol_per_m3 = 9.0,
            .aloh2_mol_per_m3 = 10.0,
            .aloh3_mol_per_m3 = 11.0,
            .aloh4_mol_per_m3 = 12.0,
            .also4_mol_per_m3 = 13.0,
            .feoh1_mol_per_m3 = 14.0,
            .feoh2_mol_per_m3 = 15.0,
            .feoh3_mol_per_m3 = 16.0,
            .feoh4_mol_per_m3 = 17.0,
            .feso4_mol_per_m3 = 18.0,
            .caoh_mol_per_m3 = 19.0,
            .caco3_mol_per_m3 = 20.0,
            .cahco3_mol_per_m3 = 21.0,
            .caso4_mol_per_m3 = 22.0,
            .mgoh_mol_per_m3 = 23.0,
            .mgco3_mol_per_m3 = 24.0,
            .mghco3_mol_per_m3 = 25.0,
            .mgso4_mol_per_m3 = 26.0,
            .naco3_mol_per_m3 = 27.0,
            .naso4_mol_per_m3 = 28.0,
            .kaso4_mol_per_m3 = 29.0,
            .po4_mol_per_m3 = 30.0,
            .h3po4_mol_per_m3 = 31.0,
            .fehpo4_mol_per_m3 = 32.0,
            .feh2po4_mol_per_m3 = 33.0,
            .capo4_mol_per_m3 = 34.0,
            .cahpo4_mol_per_m3 = 35.0,
            .cah2po4_mol_per_m3 = 36.0,
            .mghpo4_mol_per_m3 = 37.0,
        },
    };
}

test "STARTE dynamic surface-litter salt seed preserves source equations and top-layer copies" {
    var state: State = undefined;
    const inputs = testInputs();
    try std.testing.expect(try initialize(&state, inputs));
    const chy1 = std.math.pow(f64, 10.0, -(inputs.surface_litter_ph - 3.0));
    const coh1 = inputs.water_ion_product_mol2_per_m6 / chy1;
    const xcooh = @max(
        0.0,
        inputs.carboxyl_site_density_mol_per_kg_c *
            1.0e-06 *
            inputs.surface_organic_carbon_g /
            inputs.surface_soil_mass_megagram,
    );
    try std.testing.expectApproxEqAbs(chy1, state.hydrogen_activity_mol_per_m3, 1.0e-20);
    try std.testing.expectApproxEqAbs(coh1, state.hydroxide_activity_mol_per_m3, 1.0e-20);
    try std.testing.expectApproxEqAbs(xcooh, state.total_carboxyl_exchange_sites_mol_per_megagram, 1.0e-20);
    try std.testing.expectApproxEqAbs(
        xcooh * @min(1.0, chy1 / inputs.carboxyl_hydrogen_half_saturation_mol_per_m3),
        state.occupied_carboxyl_hydrogen_sites_mol_per_megagram,
        1.0e-20,
    );
    try std.testing.expectEqual(@as(f64, 0.1), state.ammonium_mol_n_per_m3);
    try std.testing.expectEqualDeep(inputs.top_layer, state.top_layer);
}

test "small surface soil mass uses zero carboxyl sites branch" {
    var state: State = undefined;
    var inputs = testInputs();
    inputs.surface_soil_mass_megagram = 1.0e-14;
    inputs.negligible_surface_soil_mass_megagram = 1.0e-12;
    _ = try initialize(&state, inputs);
    try std.testing.expectEqual(@as(f64, 0), state.total_carboxyl_exchange_sites_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 0), state.occupied_carboxyl_hydrogen_sites_mol_per_megagram);
}

test "disabled dynamic salt gate is a no-op" {
    var state: State = .{
        .hydrogen_activity_mol_per_m3 = 9.0,
        .hydroxide_activity_mol_per_m3 = 9.0,
        .total_carboxyl_exchange_sites_mol_per_megagram = 9.0,
        .occupied_carboxyl_hydrogen_sites_mol_per_megagram = 9.0,
        .effective_cation_exchange_capacity_mol_per_megagram = 9.0,
        .ammonium_mol_n_per_m3 = 9.0,
        .dissolved_aluminum_mol_per_m3 = 9.0,
        .dissolved_iron_mol_per_m3 = 9.0,
        .top_layer = testInputs().top_layer,
    };
    var inputs = testInputs();
    inputs.dynamic_salt_enabled = false;
    try std.testing.expect(!(try initialize(&state, inputs)));
    try std.testing.expectEqual(@as(f64, 9.0), state.hydrogen_activity_mol_per_m3);
}

test "invalid late top-layer value fails without mutating state" {
    var state: State = .{
        .hydrogen_activity_mol_per_m3 = 8.0,
        .hydroxide_activity_mol_per_m3 = 8.0,
        .total_carboxyl_exchange_sites_mol_per_megagram = 8.0,
        .occupied_carboxyl_hydrogen_sites_mol_per_megagram = 8.0,
        .effective_cation_exchange_capacity_mol_per_megagram = 8.0,
        .ammonium_mol_n_per_m3 = 8.0,
        .dissolved_aluminum_mol_per_m3 = 8.0,
        .dissolved_iron_mol_per_m3 = 8.0,
        .top_layer = testInputs().top_layer,
    };
    var inputs = testInputs();
    inputs.top_layer.cah2po4_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterDynamicSaltSeedInput,
        initialize(&state, inputs),
    );
    try std.testing.expectEqual(@as(f64, 8.0), state.ammonium_mol_n_per_m3);
}
