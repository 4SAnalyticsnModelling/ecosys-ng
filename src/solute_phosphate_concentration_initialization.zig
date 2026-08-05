const std = @import("std");
const phosphate_network = @import("solute_phosphate_network.zig");

const minimum_concentration_mol_per_m3 = 1.0e-32;

pub const PrimaryZoneInput = struct {
    water_volume_m3: f64,
    /// Legacy BKVLPO/BKVLPB. Normally Mg, but m3 in SOLUTE's zero-BKVL
    /// fallback; retained as an explicit normalization basis.
    adsorption_normalization_basis: f64,
    net_hpo4_change_g_p_per_step: f64,
    net_h2po4_change_g_p_per_step: f64,
    root_hpo4_uptake_g_p_per_step: f64,
    root_h2po4_uptake_g_p_per_step: f64,
    soluble_hpo4_g_p: f64,
    soluble_h2po4_g_p: f64,
    unprotonated_site_mol: f64,
    hydroxyl_site_mol: f64,
    protonated_site_mol: f64,
    adsorbed_hpo4_mol_p: f64,
    adsorbed_h2po4_mol_p: f64,
    aluminum_phosphate_solid_mol: f64,
    iron_phosphate_solid_mol: f64,
    monocalcium_phosphate_solid_mol: f64,
    dicalcium_phosphate_solid_mol: f64,
    hydroxyapatite_solid_mol: f64,
};

pub const PrimaryZoneConcentrations = struct {
    hpo4_source_mol_p_per_m3_iteration: f64,
    h2po4_source_mol_p_per_m3_iteration: f64,
    dissolved_hpo4_mol_p_per_m3: f64,
    dissolved_h2po4_mol_p_per_m3: f64,
    unprotonated_site_per_normalization_basis: f64,
    hydroxyl_site_per_normalization_basis: f64,
    protonated_site_per_normalization_basis: f64,
    adsorbed_hpo4_per_normalization_basis: f64,
    adsorbed_h2po4_per_normalization_basis: f64,
    aluminum_phosphate_solid_mol_per_m3: f64,
    iron_phosphate_solid_mol_per_m3: f64,
    monocalcium_phosphate_solid_mol_per_m3: f64,
    dicalcium_phosphate_solid_mol_per_m3: f64,
    hydroxyapatite_solid_mol_per_m3: f64,
};

pub const PrimaryPairedZoneConcentrations = struct {
    non_band: PrimaryZoneConcentrations,
    band: PrimaryZoneConcentrations,
};

/// Direct source-order translation of SOLUTE.F lines 456--524 (`RH1PX`
/// through `PCAPHB`). Inputs stored in g P are converted with the runtime
/// phosphorus molar mass; reaction sources are divided across XMRXN cycles.
pub fn initializePrimaryZoneConcentrations(
    non_band: PrimaryZoneInput,
    band: PrimaryZoneInput,
    reaction_iterations: usize,
    phosphorus_molar_mass_g_per_mol: f64,
    minimum_water_volume_m3: f64,
    numerical_floor: f64,
) !PrimaryPairedZoneConcentrations {
    try validatePrimaryInputs(
        non_band,
        reaction_iterations,
        phosphorus_molar_mass_g_per_mol,
        minimum_water_volume_m3,
        numerical_floor,
    );
    try validatePrimaryInputs(
        band,
        reaction_iterations,
        phosphorus_molar_mass_g_per_mol,
        minimum_water_volume_m3,
        numerical_floor,
    );
    const iterations: f64 = @floatFromInt(reaction_iterations);
    return .{
        .non_band = initializePrimaryZone(
            non_band,
            iterations,
            phosphorus_molar_mass_g_per_mol,
            minimum_water_volume_m3,
            numerical_floor,
        ),
        .band = initializePrimaryZone(
            band,
            iterations,
            phosphorus_molar_mass_g_per_mol,
            minimum_water_volume_m3,
            numerical_floor,
        ),
    };
}

fn initializePrimaryZone(
    input: PrimaryZoneInput,
    reaction_iterations: f64,
    phosphorus_molar_mass_g_per_mol: f64,
    minimum_water_volume_m3: f64,
    numerical_floor: f64,
) PrimaryZoneConcentrations {
    if (input.water_volume_m3 <= minimum_water_volume_m3)
        return std.mem.zeroes(PrimaryZoneConcentrations);

    const phosphorus_water_volume = phosphorus_molar_mass_g_per_mol * input.water_volume_m3;
    const hpo4_source = (input.net_hpo4_change_g_p_per_step - input.root_hpo4_uptake_g_p_per_step) /
        (reaction_iterations * phosphorus_water_volume);
    const h2po4_source = (input.net_h2po4_change_g_p_per_step - input.root_h2po4_uptake_g_p_per_step) /
        (reaction_iterations * phosphorus_water_volume);
    return .{
        .hpo4_source_mol_p_per_m3_iteration = hpo4_source,
        .h2po4_source_mol_p_per_m3_iteration = h2po4_source,
        .dissolved_hpo4_mol_p_per_m3 = @max(numerical_floor, input.soluble_hpo4_g_p / phosphorus_water_volume + hpo4_source),
        .dissolved_h2po4_mol_p_per_m3 = @max(numerical_floor, input.soluble_h2po4_g_p / phosphorus_water_volume + h2po4_source),
        .unprotonated_site_per_normalization_basis = @max(numerical_floor, input.unprotonated_site_mol) / input.adsorption_normalization_basis,
        .hydroxyl_site_per_normalization_basis = @max(numerical_floor, input.hydroxyl_site_mol) / input.adsorption_normalization_basis,
        .protonated_site_per_normalization_basis = @max(numerical_floor, input.protonated_site_mol) / input.adsorption_normalization_basis,
        .adsorbed_hpo4_per_normalization_basis = @max(numerical_floor, input.adsorbed_hpo4_mol_p) / input.adsorption_normalization_basis,
        .adsorbed_h2po4_per_normalization_basis = @max(numerical_floor, input.adsorbed_h2po4_mol_p) / input.adsorption_normalization_basis,
        .aluminum_phosphate_solid_mol_per_m3 = @max(0.0, input.aluminum_phosphate_solid_mol) / input.water_volume_m3,
        .iron_phosphate_solid_mol_per_m3 = @max(0.0, input.iron_phosphate_solid_mol) / input.water_volume_m3,
        .monocalcium_phosphate_solid_mol_per_m3 = @max(0.0, input.monocalcium_phosphate_solid_mol) / input.water_volume_m3,
        .dicalcium_phosphate_solid_mol_per_m3 = @max(0.0, input.dicalcium_phosphate_solid_mol) / input.water_volume_m3,
        .hydroxyapatite_solid_mol_per_m3 = @max(0.0, input.hydroxyapatite_solid_mol) / input.water_volume_m3,
    };
}

fn validatePrimaryInputs(
    input: PrimaryZoneInput,
    reaction_iterations: usize,
    phosphorus_molar_mass_g_per_mol: f64,
    minimum_water_volume_m3: f64,
    numerical_floor: f64,
) !void {
    if (reaction_iterations == 0 or
        !std.math.isFinite(phosphorus_molar_mass_g_per_mol) or phosphorus_molar_mass_g_per_mol <= 0 or
        !std.math.isFinite(minimum_water_volume_m3) or minimum_water_volume_m3 < 0 or
        !std.math.isFinite(numerical_floor) or numerical_floor <= 0)
        return error.InvalidPrimaryPhosphateConcentrationInput;
    inline for (@typeInfo(PrimaryZoneInput).@"struct".fields) |field|
        if (!std.math.isFinite(@field(input, field.name)))
            return error.InvalidPrimaryPhosphateConcentrationInput;
    if (input.water_volume_m3 < 0 or
        (input.water_volume_m3 > minimum_water_volume_m3 and
            input.adsorption_normalization_basis <= 0))
        return error.InvalidPrimaryPhosphateConcentrationInput;
}

fn primaryTestInput(scale: f64) PrimaryZoneInput {
    return .{
        .water_volume_m3 = 2 * scale,
        .adsorption_normalization_basis = 4 * scale,
        .net_hpo4_change_g_p_per_step = 62 * scale,
        .net_h2po4_change_g_p_per_step = 31 * scale,
        .root_hpo4_uptake_g_p_per_step = 31 * scale,
        .root_h2po4_uptake_g_p_per_step = 0,
        .soluble_hpo4_g_p = 124 * scale,
        .soluble_h2po4_g_p = 62 * scale,
        .unprotonated_site_mol = 8 * scale,
        .hydroxyl_site_mol = 12 * scale,
        .protonated_site_mol = 16 * scale,
        .adsorbed_hpo4_mol_p = 20 * scale,
        .adsorbed_h2po4_mol_p = 24 * scale,
        .aluminum_phosphate_solid_mol = 2 * scale,
        .iron_phosphate_solid_mol = 4 * scale,
        .monocalcium_phosphate_solid_mol = 6 * scale,
        .dicalcium_phosphate_solid_mol = 8 * scale,
        .hydroxyapatite_solid_mol = 10 * scale,
    };
}

test "SOLUTE 456-508 preserves paired primary phosphate reconstruction" {
    const result = try initializePrimaryZoneConcentrations(
        primaryTestInput(1),
        primaryTestInput(2),
        2,
        31,
        1.0e-12,
        1.0e-32,
    );
    const expected_source = (62.0 - 31.0) / (2.0 * (31.0 * 2.0));
    try std.testing.expectEqual(expected_source, result.non_band.hpo4_source_mol_p_per_m3_iteration);
    try std.testing.expectEqual(expected_source, result.non_band.h2po4_source_mol_p_per_m3_iteration);
    try std.testing.expectEqual(124.0 / (31.0 * 2.0) + expected_source, result.non_band.dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(62.0 / (31.0 * 2.0) + expected_source, result.non_band.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 2), result.non_band.unprotonated_site_per_normalization_basis);
    try std.testing.expectEqual(@as(f64, 5), result.non_band.adsorbed_hpo4_per_normalization_basis);
    try std.testing.expectEqual(@as(f64, 1), result.non_band.aluminum_phosphate_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 5), result.non_band.hydroxyapatite_solid_mol_per_m3);
    try std.testing.expectEqual(result.non_band, result.band);
}

test "SOLUTE 477-524 dry phosphate zones reset every derived value" {
    var non_band = primaryTestInput(1);
    var band = primaryTestInput(1);
    non_band.water_volume_m3 = 1.0e-12;
    non_band.adsorption_normalization_basis = 0;
    band.water_volume_m3 = 0;
    band.adsorption_normalization_basis = 0;
    const result = try initializePrimaryZoneConcentrations(
        non_band,
        band,
        60,
        31,
        1.0e-12,
        1.0e-32,
    );
    try std.testing.expectEqualDeep(std.mem.zeroes(PrimaryZoneConcentrations), result.non_band);
    try std.testing.expectEqualDeep(std.mem.zeroes(PrimaryZoneConcentrations), result.band);
}

test "SOLUTE primary phosphate reconstruction preserves distinct floors" {
    var non_band = primaryTestInput(1);
    non_band.soluble_hpo4_g_p = -1000;
    non_band.unprotonated_site_mol = -1;
    non_band.aluminum_phosphate_solid_mol = -1;
    const result = try initializePrimaryZoneConcentrations(
        non_band,
        primaryTestInput(1),
        60,
        31,
        0,
        1.0e-32,
    );
    try std.testing.expectEqual(@as(f64, 1.0e-32), result.non_band.dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 2.5e-33), result.non_band.unprotonated_site_per_normalization_basis);
    try std.testing.expectEqual(@as(f64, 0), result.non_band.aluminum_phosphate_solid_mol_per_m3);
}

pub const DissolvedComplexAmounts = struct {
    po4_mol_p: f64,
    h3po4_mol_p: f64,
    iron_hpo4_mol_p: f64,
    iron_h2po4_mol_p: f64,
    calcium_po4_mol_p: f64,
    calcium_hpo4_mol_p: f64,
    calcium_h2po4_mol_p: f64,
    magnesium_hpo4_mol_p: f64,
};

pub const ZoneInput = struct {
    amounts: DissolvedComplexAmounts,
    water_volume_m3: f64,
};

/// Direct translation of SOLUTE lines 729--766. Both phosphate zones are
/// validated and staged before either state is published.
pub fn applyPairedZoneConcentrations(
    non_band: *phosphate_network.State,
    band: *phosphate_network.State,
    non_band_input: ZoneInput,
    band_input: ZoneInput,
    minimum_water_volume_m3: f64,
    concentration_floor_mol_p_per_m3: f64,
) !void {
    try validateZoneInput(non_band_input);
    try validateZoneInput(band_input);
    if (!std.math.isFinite(minimum_water_volume_m3) or
        minimum_water_volume_m3 < 0 or
        !std.math.isFinite(concentration_floor_mol_p_per_m3) or
        concentration_floor_mol_p_per_m3 <= 0)
    {
        return error.InvalidPhosphateConcentrationInput;
    }

    var staged_non_band = non_band.*;
    var staged_band = band.*;
    applyZone(
        &staged_non_band,
        non_band_input,
        minimum_water_volume_m3,
        concentration_floor_mol_p_per_m3,
    );
    applyZone(
        &staged_band,
        band_input,
        minimum_water_volume_m3,
        concentration_floor_mol_p_per_m3,
    );
    try validateState(staged_non_band);
    try validateState(staged_band);
    non_band.* = staged_non_band;
    band.* = staged_band;
}

fn applyZone(
    state: *phosphate_network.State,
    input: ZoneInput,
    minimum_water_volume_m3: f64,
    concentration_floor_mol_p_per_m3: f64,
) void {
    if (input.water_volume_m3 <= minimum_water_volume_m3) {
        state.dissolved_po4_mol_p_per_m3 = 0;
        state.dissolved_h3po4_mol_p_per_m3 = 0;
        state.iron_hpo4_pair_mol_per_m3 = 0;
        state.iron_h2po4_pair_mol_per_m3 = 0;
        state.calcium_po4_pair_mol_per_m3 = 0;
        state.calcium_hpo4_pair_mol_per_m3 = 0;
        state.calcium_h2po4_pair_mol_per_m3 = 0;
        state.magnesium_hpo4_pair_mol_per_m3 = 0;
        return;
    }
    const amounts = input.amounts;
    const volume = input.water_volume_m3;
    state.dissolved_po4_mol_p_per_m3 =
        concentration(amounts.po4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.dissolved_h3po4_mol_p_per_m3 =
        concentration(amounts.h3po4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.iron_hpo4_pair_mol_per_m3 =
        concentration(amounts.iron_hpo4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.iron_h2po4_pair_mol_per_m3 =
        concentration(amounts.iron_h2po4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.calcium_po4_pair_mol_per_m3 =
        concentration(amounts.calcium_po4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.calcium_hpo4_pair_mol_per_m3 =
        concentration(amounts.calcium_hpo4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.calcium_h2po4_pair_mol_per_m3 =
        concentration(amounts.calcium_h2po4_mol_p, volume, concentration_floor_mol_p_per_m3);
    state.magnesium_hpo4_pair_mol_per_m3 =
        concentration(amounts.magnesium_hpo4_mol_p, volume, concentration_floor_mol_p_per_m3);
}

fn concentration(
    amount_mol_p: f64,
    water_volume_m3: f64,
    concentration_floor_mol_p_per_m3: f64,
) f64 {
    return @max(
        concentration_floor_mol_p_per_m3,
        amount_mol_p / water_volume_m3,
    );
}

fn validateZoneInput(input: ZoneInput) !void {
    if (!std.math.isFinite(input.water_volume_m3) or
        input.water_volume_m3 < 0)
    {
        return error.InvalidPhosphateConcentrationInput;
    }
    inline for (@typeInfo(DissolvedComplexAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(input.amounts, field.name)))
            return error.InvalidPhosphateConcentrationInput;
    }
}

fn validateState(state: phosphate_network.State) !void {
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
        const value = @field(state, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPhosphateConcentrationState;
    }
}

fn scaledAmounts(scale: f64) DissolvedComplexAmounts {
    return .{
        .po4_mol_p = 1 * scale,
        .h3po4_mol_p = 2 * scale,
        .iron_hpo4_mol_p = 3 * scale,
        .iron_h2po4_mol_p = 4 * scale,
        .calcium_po4_mol_p = 5 * scale,
        .calcium_hpo4_mol_p = 6 * scale,
        .calcium_h2po4_mol_p = 7 * scale,
        .magnesium_hpo4_mol_p = 8 * scale,
    };
}

fn filledState(value: f64) phosphate_network.State {
    var state: phosphate_network.State = undefined;
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field|
        @field(state, field.name) = value;
    return state;
}

test "SOLUTE paired phosphate concentrations preserve zone conversion and owners" {
    var non_band = filledState(11);
    var band = filledState(13);
    var non_band_amounts = scaledAmounts(2);
    non_band_amounts.iron_hpo4_mol_p = -1;
    try applyPairedZoneConcentrations(
        &non_band,
        &band,
        .{ .amounts = non_band_amounts, .water_volume_m3 = 2 },
        .{ .amounts = scaledAmounts(3), .water_volume_m3 = 3 },
        0,
        2.0e-30,
    );

    try std.testing.expectEqual(@as(f64, 1), non_band.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(
        @as(f64, 2.0e-30),
        non_band.iron_hpo4_pair_mol_per_m3,
    );
    try std.testing.expectEqual(@as(f64, 8), non_band.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), band.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 8), band.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 11), non_band.dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 13), band.protonated_site_mol_per_megagram);
}

test "SOLUTE dry phosphate zone clears only reconstructed complexes" {
    var non_band = filledState(5);
    var band = filledState(7);
    try applyPairedZoneConcentrations(
        &non_band,
        &band,
        .{ .amounts = scaledAmounts(1), .water_volume_m3 = 1.0e-9 },
        .{ .amounts = scaledAmounts(1), .water_volume_m3 = 2 },
        1.0e-9,
        1.0e-32,
    );
    try std.testing.expectEqual(@as(f64, 0), non_band.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 0), non_band.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 5), non_band.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 0.5), band.dissolved_po4_mol_p_per_m3);
}

test "SOLUTE paired phosphate concentration failure is atomic" {
    var non_band = filledState(2);
    var band = filledState(3);
    const non_band_before = non_band;
    const band_before = band;
    var invalid = scaledAmounts(1);
    invalid.calcium_h2po4_mol_p = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidPhosphateConcentrationInput,
        applyPairedZoneConcentrations(
            &non_band,
            &band,
            .{ .amounts = scaledAmounts(1), .water_volume_m3 = 1 },
            .{ .amounts = invalid, .water_volume_m3 = 1 },
            0,
            1.0e-32,
        ),
    );
    try std.testing.expectEqualDeep(non_band_before, non_band);
    try std.testing.expectEqualDeep(band_before, band);
}
