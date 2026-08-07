const std = @import("std");
const phosphate_exchange = @import("phosphate_exchange.zig");

pub const ZoneState = struct {
    dissolved_h2po4_mol_p_per_m3: f64,
    dissolved_hpo4_mol_p_per_m3: f64,
    deprotonated_site_mol_per_megagram: f64,
    hydroxyl_site_mol_per_megagram: f64,
    protonated_site_mol_per_megagram: f64,
    adsorbed_h2po4_mol_p_per_megagram: f64,
    adsorbed_hpo4_mol_p_per_megagram: f64,
};

pub const ZoneUpdate = struct {
    state: *ZoneState,
    flux: phosphate_exchange.Flux,
    soil_mass_per_water_volume_megagrams_per_m3: f64,
};

/// Applies the non-band and fertilizer-band exchange reactions as one transaction.
/// A failed update leaves both zones unchanged.
pub fn commitPaired(non_band: ZoneUpdate, band: ZoneUpdate) !void {
    const staged_non_band = try stage(non_band);
    const staged_band = try stage(band);
    non_band.state.* = staged_non_band;
    band.state.* = staged_band;
}

pub fn commit(update: ZoneUpdate) !void {
    update.state.* = try stage(update);
}

fn stage(update: ZoneUpdate) !ZoneState {
    try validateState(update.state.*);
    if (!std.math.isFinite(update.soil_mass_per_water_volume_megagrams_per_m3) or
        update.soil_mass_per_water_volume_megagrams_per_m3 <= 0)
        return error.InvalidSoilWaterRatio;

    const before_phosphorus = phosphorusMolPerMg(update.state.*, update.soil_mass_per_water_volume_megagrams_per_m3);
    const before_sites = siteMolPerMg(update.state.*);
    const f = update.flux;
    inline for (@typeInfo(phosphate_exchange.Flux).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(f, field.name))) return error.NonFinitePhosphateExchangeFlux;
    }

    const water_volume_per_soil_mass_m3_per_megagram = 1.0 / update.soil_mass_per_water_volume_megagrams_per_m3;
    var next = update.state.*;
    next.dissolved_h2po4_mol_p_per_m3 -=
        (f.h2po4_with_protonated_site_mol_p_per_megagram + f.h2po4_with_hydroxyl_site_mol_p_per_megagram) /
        water_volume_per_soil_mass_m3_per_megagram;
    next.dissolved_hpo4_mol_p_per_m3 -=
        f.hpo4_with_hydroxyl_site_mol_p_per_megagram / water_volume_per_soil_mass_m3_per_megagram;
    next.deprotonated_site_mol_per_megagram -= f.hydroxyl_to_deprotonated_site_mol_per_megagram;
    next.hydroxyl_site_mol_per_megagram += f.hydroxyl_to_deprotonated_site_mol_per_megagram -
        f.protonated_to_hydroxyl_site_mol_per_megagram -
        f.h2po4_with_hydroxyl_site_mol_p_per_megagram -
        f.hpo4_with_hydroxyl_site_mol_p_per_megagram;
    next.protonated_site_mol_per_megagram += f.protonated_to_hydroxyl_site_mol_per_megagram -
        f.h2po4_with_protonated_site_mol_p_per_megagram;
    next.adsorbed_h2po4_mol_p_per_megagram +=
        f.h2po4_with_protonated_site_mol_p_per_megagram + f.h2po4_with_hydroxyl_site_mol_p_per_megagram;
    next.adsorbed_hpo4_mol_p_per_megagram += f.hpo4_with_hydroxyl_site_mol_p_per_megagram;

    try validateState(next);
    try expectConserved(before_phosphorus, phosphorusMolPerMg(next, update.soil_mass_per_water_volume_megagrams_per_m3));
    try expectConserved(before_sites, siteMolPerMg(next));
    return next;
}

fn phosphorusMolPerMg(state: ZoneState, soil_mass_per_water_volume_megagrams_per_m3: f64) f64 {
    return (state.dissolved_h2po4_mol_p_per_m3 + state.dissolved_hpo4_mol_p_per_m3) /
        soil_mass_per_water_volume_megagrams_per_m3 +
        state.adsorbed_h2po4_mol_p_per_megagram + state.adsorbed_hpo4_mol_p_per_megagram;
}

fn siteMolPerMg(state: ZoneState) f64 {
    return state.deprotonated_site_mol_per_megagram + state.hydroxyl_site_mol_per_megagram +
        state.protonated_site_mol_per_megagram + state.adsorbed_h2po4_mol_p_per_megagram +
        state.adsorbed_hpo4_mol_p_per_megagram;
}

fn validateState(state: ZoneState) !void {
    inline for (@typeInfo(ZoneState).@"struct".fields) |field| {
        const value = @field(state, field.name);
        if (!std.math.isFinite(value)) return error.NonFinitePhosphateState;
        if (value < -1e-12) return error.NegativePhosphateState;
    }
}

fn expectConserved(before: f64, after: f64) !void {
    const tolerance = 1e-11 * @max(1.0, @abs(before), @abs(after));
    if (@abs(after - before) > tolerance) return error.PhosphateExchangeConservationFailure;
}

fn sampleState() ZoneState {
    return .{
        .dissolved_h2po4_mol_p_per_m3 = 0.8,
        .dissolved_hpo4_mol_p_per_m3 = 0.4,
        .deprotonated_site_mol_per_megagram = 0.3,
        .hydroxyl_site_mol_per_megagram = 0.5,
        .protonated_site_mol_per_megagram = 0.4,
        .adsorbed_h2po4_mol_p_per_megagram = 0.1,
        .adsorbed_hpo4_mol_p_per_megagram = 0.05,
    };
}

test "phosphate exchange conserves phosphorus and adsorption sites with explicit units" {
    var state = sampleState();
    const density = 2.0;
    const phosphorus_before = phosphorusMolPerMg(state, density);
    const sites_before = siteMolPerMg(state);
    try commit(.{ .state = &state, .soil_mass_per_water_volume_megagrams_per_m3 = density, .flux = .{
        .protonated_to_hydroxyl_site_mol_per_megagram = 0.02,
        .hydroxyl_to_deprotonated_site_mol_per_megagram = -0.01,
        .h2po4_with_protonated_site_mol_p_per_megagram = 0.03,
        .h2po4_with_hydroxyl_site_mol_p_per_megagram = 0.02,
        .hpo4_with_hydroxyl_site_mol_p_per_megagram = 0.01,
    } });
    try std.testing.expectApproxEqAbs(phosphorus_before, phosphorusMolPerMg(state, density), 1e-12);
    try std.testing.expectApproxEqAbs(sites_before, siteMolPerMg(state), 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), state.dissolved_h2po4_mol_p_per_m3, 1e-12);
}

test "paired phosphate exchange rolls both zones back when either update is invalid" {
    var non_band = sampleState();
    var band = sampleState();
    const before_non_band = non_band;
    const before_band = band;
    const valid_flux: phosphate_exchange.Flux = .{
        .protonated_to_hydroxyl_site_mol_per_megagram = 0,
        .hydroxyl_to_deprotonated_site_mol_per_megagram = 0,
        .h2po4_with_protonated_site_mol_p_per_megagram = 0.01,
        .h2po4_with_hydroxyl_site_mol_p_per_megagram = 0,
        .hpo4_with_hydroxyl_site_mol_p_per_megagram = 0,
    };
    var invalid_flux = valid_flux;
    invalid_flux.h2po4_with_protonated_site_mol_p_per_megagram = 1.0;
    try std.testing.expectError(error.NegativePhosphateState, commitPaired(
        .{ .state = &non_band, .flux = valid_flux, .soil_mass_per_water_volume_megagrams_per_m3 = 1 },
        .{ .state = &band, .flux = invalid_flux, .soil_mass_per_water_volume_megagrams_per_m3 = 1 },
    ));
    try std.testing.expectEqualDeep(before_non_band, non_band);
    try std.testing.expectEqualDeep(before_band, band);
}
