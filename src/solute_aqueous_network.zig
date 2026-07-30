const std = @import("std");

pub const Fluxes = struct {
    ammonium_non_band_association: f64,
    ammonium_band_association: f64,
    carbonate_hydrogen_association: f64,
    bicarbonate_hydrogen_association: f64,
    aluminum_hydroxide_1_association: f64,
    aluminum_hydroxide_2_association: f64,
    aluminum_hydroxide_3_association: f64,
    aluminum_hydroxide_4_association: f64,
    aluminum_sulfate_association: f64,
    iron_hydroxide_1_association: f64,
    iron_hydroxide_2_association: f64,
    iron_hydroxide_3_association: f64,
    iron_hydroxide_4_association: f64,
    iron_sulfate_association: f64,
    calcium_hydroxide_association: f64,
    calcium_carbonate_association: f64,
    calcium_bicarbonate_association: f64,
    calcium_sulfate_association: f64,
    magnesium_hydroxide_association: f64,
    magnesium_carbonate_association: f64,
    magnesium_bicarbonate_association: f64,
    magnesium_sulfate_association: f64,
    sodium_carbonate_association: f64,
    sodium_sulfate_association: f64,
    potassium_sulfate_association: f64,
};

pub const AmmoniumWaterFractions = struct {
    non_band: f64,
    band: f64,
};

pub const Transformations = struct {
    ammonium_non_band: f64,
    ammonia_non_band: f64,
    ammonium_band: f64,
    ammonia_band: f64,
    hydrogen: f64,
    hydroxide: f64,
    aluminum: f64,
    aluminum_hydroxide_1: f64,
    aluminum_hydroxide_2: f64,
    aluminum_hydroxide_3: f64,
    aluminum_hydroxide_4: f64,
    aluminum_sulfate: f64,
    iron: f64,
    iron_hydroxide_1: f64,
    iron_hydroxide_2: f64,
    iron_hydroxide_3: f64,
    iron_hydroxide_4: f64,
    iron_sulfate: f64,
    calcium: f64,
    calcium_hydroxide: f64,
    calcium_carbonate: f64,
    calcium_bicarbonate: f64,
    calcium_sulfate: f64,
    magnesium: f64,
    magnesium_hydroxide: f64,
    magnesium_carbonate: f64,
    magnesium_bicarbonate: f64,
    magnesium_sulfate: f64,
    sodium: f64,
    sodium_carbonate: f64,
    sodium_sulfate: f64,
    potassium: f64,
    potassium_sulfate: f64,
    sulfate: f64,
    carbonate: f64,
    bicarbonate: f64,
    carbon_dioxide: f64,
    hydrogen_silicate: f64,
    chloride: f64,
    nitrate_non_band: f64,
    nitrate_band: f64,
};

pub const State = Transformations;

pub fn commit(state: *State, transformations: Transformations) !void {
    try validateState(state.*);
    var next = state.*;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const change = @field(transformations, field.name);
        if (!std.math.isFinite(change)) return error.NonFiniteAqueousTransformation;
        @field(next, field.name) += change;
    }
    try validateState(next);
    state.* = next;
}

/// Exact non-phosphate ion-pairing portion of SOLUTE.F's transformation block.
pub fn assemble(f: Fluxes, ammonium_water_fractions: AmmoniumWaterFractions) !Transformations {
    inline for (@typeInfo(Fluxes).@"struct".fields) |field| if (!std.math.isFinite(@field(f, field.name))) return error.NonFiniteAqueousReactionFlux;
    inline for (@typeInfo(AmmoniumWaterFractions).@"struct".fields) |field| {
        const fraction = @field(ammonium_water_fractions, field.name);
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidAmmoniumWaterFraction;
    }
    return .{
        .ammonium_non_band = f.ammonium_non_band_association,
        .ammonia_non_band = -f.ammonium_non_band_association,
        .ammonium_band = f.ammonium_band_association,
        .ammonia_band = -f.ammonium_band_association,
        .hydrogen = -f.bicarbonate_hydrogen_association - f.carbonate_hydrogen_association - f.ammonium_non_band_association * ammonium_water_fractions.non_band - f.ammonium_band_association * ammonium_water_fractions.band,
        .hydroxide = -f.calcium_hydroxide_association - f.magnesium_hydroxide_association - f.aluminum_hydroxide_1_association - f.aluminum_hydroxide_2_association - f.aluminum_hydroxide_3_association - f.aluminum_hydroxide_4_association - f.iron_hydroxide_1_association - f.iron_hydroxide_2_association - f.iron_hydroxide_3_association - f.iron_hydroxide_4_association,
        .aluminum = -f.aluminum_hydroxide_1_association - f.aluminum_sulfate_association,
        .aluminum_hydroxide_1 = f.aluminum_hydroxide_1_association - f.aluminum_hydroxide_2_association,
        .aluminum_hydroxide_2 = f.aluminum_hydroxide_2_association - f.aluminum_hydroxide_3_association,
        .aluminum_hydroxide_3 = f.aluminum_hydroxide_3_association - f.aluminum_hydroxide_4_association,
        .aluminum_hydroxide_4 = f.aluminum_hydroxide_4_association,
        .aluminum_sulfate = f.aluminum_sulfate_association,
        .iron = -f.iron_hydroxide_1_association - f.iron_sulfate_association,
        .iron_hydroxide_1 = f.iron_hydroxide_1_association - f.iron_hydroxide_2_association,
        .iron_hydroxide_2 = f.iron_hydroxide_2_association - f.iron_hydroxide_3_association,
        .iron_hydroxide_3 = f.iron_hydroxide_3_association - f.iron_hydroxide_4_association,
        .iron_hydroxide_4 = f.iron_hydroxide_4_association,
        .iron_sulfate = f.iron_sulfate_association,
        .calcium = -f.calcium_hydroxide_association - f.calcium_carbonate_association - f.calcium_bicarbonate_association - f.calcium_sulfate_association,
        .calcium_hydroxide = f.calcium_hydroxide_association,
        .calcium_carbonate = f.calcium_carbonate_association,
        .calcium_bicarbonate = f.calcium_bicarbonate_association,
        .calcium_sulfate = f.calcium_sulfate_association,
        .magnesium = -f.magnesium_hydroxide_association - f.magnesium_carbonate_association - f.magnesium_bicarbonate_association - f.magnesium_sulfate_association,
        .magnesium_hydroxide = f.magnesium_hydroxide_association,
        .magnesium_carbonate = f.magnesium_carbonate_association,
        .magnesium_bicarbonate = f.magnesium_bicarbonate_association,
        .magnesium_sulfate = f.magnesium_sulfate_association,
        .sodium = -f.sodium_carbonate_association - f.sodium_sulfate_association,
        .sodium_carbonate = f.sodium_carbonate_association,
        .sodium_sulfate = f.sodium_sulfate_association,
        .potassium = -f.potassium_sulfate_association,
        .potassium_sulfate = f.potassium_sulfate_association,
        .sulfate = -f.aluminum_sulfate_association - f.iron_sulfate_association - f.calcium_sulfate_association - f.magnesium_sulfate_association - f.sodium_sulfate_association - f.potassium_sulfate_association,
        .carbonate = -f.carbonate_hydrogen_association - f.calcium_carbonate_association - f.magnesium_carbonate_association - f.sodium_carbonate_association,
        .bicarbonate = -f.bicarbonate_hydrogen_association - f.calcium_bicarbonate_association - f.magnesium_bicarbonate_association + f.carbonate_hydrogen_association,
        .carbon_dioxide = f.bicarbonate_hydrogen_association,
        .hydrogen_silicate = 0,
        .chloride = 0,
        .nitrate_non_band = 0,
        .nitrate_band = 0,
    };
}

fn sampleFluxes() Fluxes {
    var result: Fluxes = undefined;
    var value: f64 = -0.012;
    inline for (@typeInfo(Fluxes).@"struct".fields) |field| {
        @field(result, field.name) = value;
        value += 0.001;
    }
    return result;
}

fn validateState(state: State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const value = @field(state, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteAqueousState;
        if (value < -1e-12) return error.NegativeAqueousState;
    }
}

fn filledState(value: f64) State {
    var result: State = undefined;
    inline for (@typeInfo(State).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "aqueous reaction ledger conserves every paired element" {
    const t = try assemble(sampleFluxes(), .{ .non_band = 0.8, .band = 0.2 });
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.ammonium_non_band + t.ammonia_non_band, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.ammonium_band + t.ammonia_band, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.aluminum + t.aluminum_hydroxide_1 + t.aluminum_hydroxide_2 + t.aluminum_hydroxide_3 + t.aluminum_hydroxide_4 + t.aluminum_sulfate, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.iron + t.iron_hydroxide_1 + t.iron_hydroxide_2 + t.iron_hydroxide_3 + t.iron_hydroxide_4 + t.iron_sulfate, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.calcium + t.calcium_hydroxide + t.calcium_carbonate + t.calcium_bicarbonate + t.calcium_sulfate, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.magnesium + t.magnesium_hydroxide + t.magnesium_carbonate + t.magnesium_bicarbonate + t.magnesium_sulfate, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.sodium + t.sodium_carbonate + t.sodium_sulfate, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), t.potassium + t.potassium_sulfate, 1e-14);
    const sulfur = t.sulfate + t.aluminum_sulfate + t.iron_sulfate + t.calcium_sulfate + t.magnesium_sulfate + t.sodium_sulfate + t.potassium_sulfate;
    try std.testing.expectApproxEqAbs(@as(f64, 0), sulfur, 1e-14);
    const carbonate_carbon = t.carbonate + t.bicarbonate + t.carbon_dioxide + t.calcium_carbonate + t.calcium_bicarbonate + t.magnesium_carbonate + t.magnesium_bicarbonate + t.sodium_carbonate;
    try std.testing.expectApproxEqAbs(@as(f64, 0), carbonate_carbon, 1e-14);
}

test "shared hydrogen uses runtime ammonium-zone water fractions" {
    var fluxes = sampleFluxes();
    inline for (@typeInfo(Fluxes).@"struct".fields) |field| @field(fluxes, field.name) = 0;
    fluxes.ammonium_non_band_association = 0.2;
    fluxes.ammonium_band_association = 0.4;
    const transformations = try assemble(fluxes, .{ .non_band = 0.75, .band = 0.25 });
    try std.testing.expectApproxEqAbs(@as(f64, -0.25), transformations.hydrogen, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), transformations.ammonium_non_band, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), transformations.ammonium_band, 1e-15);
}

test "aqueous state transaction rolls back all species on failure" {
    var state = filledState(1);
    const before = state;
    var transformations = filledState(0);
    transformations.calcium = -2;
    transformations.calcium_sulfate = 2;
    try std.testing.expectError(error.NegativeAqueousState, commit(&state, transformations));
    try std.testing.expectEqualDeep(before, state);
}
