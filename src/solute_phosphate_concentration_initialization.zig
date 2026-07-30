const std = @import("std");
const phosphate_network = @import("solute_phosphate_network.zig");

const minimum_concentration_mol_per_m3 = 1.0e-32;

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
) !void {
    try validateZoneInput(non_band_input);
    try validateZoneInput(band_input);
    if (!std.math.isFinite(minimum_water_volume_m3) or
        minimum_water_volume_m3 < 0)
    {
        return error.InvalidPhosphateConcentrationInput;
    }

    var staged_non_band = non_band.*;
    var staged_band = band.*;
    applyZone(
        &staged_non_band,
        non_band_input,
        minimum_water_volume_m3,
    );
    applyZone(
        &staged_band,
        band_input,
        minimum_water_volume_m3,
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
        concentration(amounts.po4_mol_p, volume);
    state.dissolved_h3po4_mol_p_per_m3 =
        concentration(amounts.h3po4_mol_p, volume);
    state.iron_hpo4_pair_mol_per_m3 =
        concentration(amounts.iron_hpo4_mol_p, volume);
    state.iron_h2po4_pair_mol_per_m3 =
        concentration(amounts.iron_h2po4_mol_p, volume);
    state.calcium_po4_pair_mol_per_m3 =
        concentration(amounts.calcium_po4_mol_p, volume);
    state.calcium_hpo4_pair_mol_per_m3 =
        concentration(amounts.calcium_hpo4_mol_p, volume);
    state.calcium_h2po4_pair_mol_per_m3 =
        concentration(amounts.calcium_h2po4_mol_p, volume);
    state.magnesium_hpo4_pair_mol_per_m3 =
        concentration(amounts.magnesium_hpo4_mol_p, volume);
}

fn concentration(amount_mol_p: f64, water_volume_m3: f64) f64 {
    return @max(
        minimum_concentration_mol_per_m3,
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
    );

    try std.testing.expectEqual(@as(f64, 1), non_band.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(
        minimum_concentration_mol_per_m3,
        non_band.iron_hpo4_pair_mol_per_m3,
    );
    try std.testing.expectEqual(@as(f64, 8), non_band.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), band.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 8), band.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 11), non_band.dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 13), band.protonated_site_mol_per_Mg);
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
        ),
    );
    try std.testing.expectEqualDeep(non_band_before, non_band);
    try std.testing.expectEqualDeep(band_before, band);
}
