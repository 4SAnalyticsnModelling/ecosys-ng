const std = @import("std");

pub const ActivityTerms = struct {
    ammonium: f64,
    hydrogen: f64,
    aluminum_cube_root: f64,
    iron_cube_root: f64,
    calcium_square_root: f64,
    magnesium_square_root: f64,
    sodium: f64,
    potassium: f64,
};

pub const Selectivity = struct {
    calcium_ammonium: f64,
    calcium_hydrogen: f64,
    calcium_trivalent: f64,
    calcium_magnesium: f64,
    calcium_sodium: f64,
    calcium_potassium: f64,
};

/// Exchange concentrations are `mol Mg-1`; trivalent and divalent values are
/// converted from charge-equivalent partitions using exact source divisors.
pub const ExchangeConcentrations = struct {
    ammonium: f64,
    hydrogen: f64,
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
};

pub const Diagnostics = struct {
    calcium_basis_mol_per_megagram: f64, // XCAX
    raw_charge_total_mol_per_megagram: f64, // XTLQ
    normalization_factor: f64, // FX
};

pub const Result = struct { exchange: ExchangeConcentrations, diagnostics: ?Diagnostics };

fn zeroResult() Result {
    return .{ .exchange = .{ .ammonium = 0, .hydrogen = 0, .aluminum = 0, .iron = 0, .calcium = 0, .magnesium = 0, .sodium = 0, .potassium = 0 }, .diagnostics = null };
}

/// Direct translation of `starte.f` lines 1822--1886. A null result means the
/// source `M == 1` branch was inactive; zero exchange is the low-organic-C branch.
pub fn partition(
    iteration: usize,
    surface_organic_carbon_g: f64,
    negligible_carbon_g: f64,
    exchange_capacity_mol_per_megagram: f64,
    activities: ActivityTerms,
    selectivity: Selectivity,
) !?Result {
    if (iteration != 1) return null;
    if (!std.math.isFinite(surface_organic_carbon_g) or surface_organic_carbon_g < 0 or
        !std.math.isFinite(negligible_carbon_g) or negligible_carbon_g < 0)
        return error.InvalidSurfaceLitterCarbon;
    if (surface_organic_carbon_g <= negligible_carbon_g) return zeroResult();
    if (!std.math.isFinite(exchange_capacity_mol_per_megagram) or exchange_capacity_mol_per_megagram < 0)
        return error.InvalidExchangeCapacity;
    inline for (@typeInfo(ActivityTerms).@"struct".fields) |field| {
        const value = @field(activities, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidGaponActivity;
    }
    if (activities.calcium_square_root <= 0) return error.InvalidGaponCalciumActivity;
    inline for (@typeInfo(Selectivity).@"struct".fields) |field| {
        const value = @field(selectivity, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidGaponSelectivity;
    }

    const calcium_basis = exchange_capacity_mol_per_megagram /
        (1.0 + selectivity.calcium_ammonium * activities.ammonium / activities.calcium_square_root +
            selectivity.calcium_hydrogen * activities.hydrogen / activities.calcium_square_root +
            selectivity.calcium_trivalent * activities.aluminum_cube_root / activities.calcium_square_root +
            selectivity.calcium_trivalent * activities.iron_cube_root / activities.calcium_square_root +
            selectivity.calcium_magnesium * activities.magnesium_square_root / activities.calcium_square_root +
            selectivity.calcium_sodium * activities.sodium / activities.calcium_square_root +
            selectivity.calcium_potassium * activities.potassium / activities.calcium_square_root);
    var ammonium = calcium_basis * activities.ammonium * selectivity.calcium_ammonium;
    var hydrogen = calcium_basis * activities.hydrogen * selectivity.calcium_hydrogen;
    var aluminum = calcium_basis * activities.aluminum_cube_root * selectivity.calcium_trivalent;
    var iron = calcium_basis * activities.iron_cube_root * selectivity.calcium_trivalent;
    var calcium = calcium_basis * activities.calcium_square_root;
    var magnesium = calcium_basis * activities.magnesium_square_root * selectivity.calcium_magnesium;
    var sodium = calcium_basis * activities.sodium * selectivity.calcium_sodium;
    var potassium = calcium_basis * activities.potassium * selectivity.calcium_potassium;
    const raw_total = ammonium + aluminum + iron + calcium + magnesium + hydrogen + sodium + potassium;
    const factor = if (raw_total > 0.0) exchange_capacity_mol_per_megagram / raw_total else 0.0;
    ammonium = factor * ammonium;
    hydrogen = factor * hydrogen;
    aluminum = factor * aluminum;
    iron = factor * iron;
    calcium = factor * calcium;
    magnesium = factor * magnesium;
    sodium = factor * sodium;
    potassium = factor * potassium;
    const result: Result = .{
        .exchange = .{ .ammonium = ammonium, .hydrogen = hydrogen, .aluminum = aluminum / 3.0, .iron = iron / 3.0, .calcium = calcium / 2.0, .magnesium = magnesium / 2.0, .sodium = sodium, .potassium = potassium },
        .diagnostics = .{ .calcium_basis_mol_per_megagram = calcium_basis, .raw_charge_total_mol_per_megagram = raw_total, .normalization_factor = factor },
    };
    inline for (@typeInfo(ExchangeConcentrations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.exchange, field.name))) return error.InvalidGaponPartition;
    return result;
}

fn ones(comptime T: type) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = 1;
    return result;
}

test "STARTE initial Gapon partition preserves denominator normalization and valence division" {
    const result = (try partition(1, 1, 0, 8, ones(ActivityTerms), ones(Selectivity))).?;
    try std.testing.expectEqual(@as(f64, 1), result.diagnostics.?.calcium_basis_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 8), result.diagnostics.?.raw_charge_total_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), result.exchange.ammonium);
    try std.testing.expectEqual(@as(f64, 1.0 / 3.0), result.exchange.aluminum);
    try std.testing.expectEqual(@as(f64, 0.5), result.exchange.calcium);
    const charge_total = result.exchange.ammonium + result.exchange.hydrogen + 3 * result.exchange.aluminum + 3 * result.exchange.iron + 2 * result.exchange.calcium + 2 * result.exchange.magnesium + result.exchange.sodium + result.exchange.potassium;
    try std.testing.expectApproxEqAbs(@as(f64, 8), charge_total, 1e-14);
}

test "STARTE low organic carbon branch zeros exchange before dormant activity validation" {
    const nan = std.math.nan(f64);
    const result = (try partition(1, 0, 0, nan, .{ .ammonium = nan, .hydrogen = nan, .aluminum_cube_root = nan, .iron_cube_root = nan, .calcium_square_root = nan, .magnesium_square_root = nan, .sodium = nan, .potassium = nan }, .{ .calcium_ammonium = nan, .calcium_hydrogen = nan, .calcium_trivalent = nan, .calcium_magnesium = nan, .calcium_sodium = nan, .calcium_potassium = nan })).?;
    try std.testing.expectEqual(@as(f64, 0), result.exchange.calcium);
    try std.testing.expect(result.diagnostics == null);
}

test "STARTE later Gapon iterations bypass all dormant inputs" {
    const nan = std.math.nan(f64);
    try std.testing.expectEqual(@as(?Result, null), try partition(2, nan, nan, nan, undefined, undefined));
}
