const std = @import("std");

pub const CationExchangeConcentrations = struct {
    ammonium_mol_n_per_megagram: f64,
    hydrogen_mol_per_megagram: f64,
    aluminum_mol_per_megagram: f64,
    iron_mol_per_megagram: f64,
    calcium_mol_per_megagram: f64,
    magnesium_mol_per_megagram: f64,
    sodium_mol_per_megagram: f64,
    potassium_mol_per_megagram: f64,
};

pub const CationExchangeInventories = struct {
    non_band_ammonium_mol_n: f64,
    band_ammonium_mol_n: f64,
    hydrogen_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
};

pub const AnionExchangeInventories = struct {
    unprotonated_hydroxyl_sites_mol: f64,
    singly_protonated_hydroxyl_sites_mol: f64,
    doubly_protonated_hydroxyl_sites_mol: f64,
    hydrogen_phosphate_sites_mol_p: f64,
    dihydrogen_phosphate_sites_mol_p: f64,
};

pub const Result = struct {
    cations: CationExchangeInventories,
    anions: AnionExchangeInventories,
};

/// Direct translation of STARTE lines 2054--2067 for one surface litter cell.
/// Exchange concentrations in mol Mg-1 are converted to extensive inventories
/// using runtime litter dry mass. Surface band and anion sites start at zero.
pub fn initialize(
    litter_dry_mass_megagrams: f64,
    concentrations: CationExchangeConcentrations,
) !Result {
    if (!std.math.isFinite(litter_dry_mass_megagrams))
        return error.NonFiniteSurfaceExchangeDryMass;
    if (litter_dry_mass_megagrams < 0)
        return error.InvalidSurfaceExchangeDryMass;
    inline for (std.meta.fields(CationExchangeConcentrations)) |field| {
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceExchangeConcentration;
        if (value < 0)
            return error.InvalidSurfaceExchangeConcentration;
    }

    const result: Result = .{
        .cations = .{
            .non_band_ammonium_mol_n = concentrations.ammonium_mol_n_per_megagram * litter_dry_mass_megagrams,
            .band_ammonium_mol_n = 0,
            .hydrogen_mol = concentrations.hydrogen_mol_per_megagram * litter_dry_mass_megagrams,
            .aluminum_mol = concentrations.aluminum_mol_per_megagram * litter_dry_mass_megagrams,
            .iron_mol = concentrations.iron_mol_per_megagram * litter_dry_mass_megagrams,
            .calcium_mol = concentrations.calcium_mol_per_megagram * litter_dry_mass_megagrams,
            .magnesium_mol = concentrations.magnesium_mol_per_megagram * litter_dry_mass_megagrams,
            .sodium_mol = concentrations.sodium_mol_per_megagram * litter_dry_mass_megagrams,
            .potassium_mol = concentrations.potassium_mol_per_megagram * litter_dry_mass_megagrams,
        },
        .anions = std.mem.zeroes(AnionExchangeInventories),
    };
    inline for (std.meta.fields(CationExchangeInventories)) |field| {
        if (!std.math.isFinite(@field(result.cations, field.name)))
            return error.NonFiniteSurfaceExchangeInventory;
    }
    return result;
}

fn sequentialConcentrations() CationExchangeConcentrations {
    return .{
        .ammonium_mol_n_per_megagram = 1,
        .hydrogen_mol_per_megagram = 2,
        .aluminum_mol_per_megagram = 3,
        .iron_mol_per_megagram = 4,
        .calcium_mol_per_megagram = 5,
        .magnesium_mol_per_megagram = 6,
        .sodium_mol_per_megagram = 7,
        .potassium_mol_per_megagram = 8,
    };
}

test "STARTE surface exchange concentrations scale by runtime litter mass" {
    const result = try initialize(2.5, sequentialConcentrations());
    try std.testing.expectEqual(@as(f64, 2.5), result.cations.non_band_ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 5), result.cations.hydrogen_mol);
    try std.testing.expectEqual(@as(f64, 20), result.cations.potassium_mol);
}

test "STARTE surface band ammonium and anion exchange sites initialize to zero" {
    const result = try initialize(1, sequentialConcentrations());
    try std.testing.expectEqual(@as(f64, 0), result.cations.band_ammonium_mol_n);
    inline for (std.meta.fields(AnionExchangeInventories)) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(result.anions, field.name));
}

test "STARTE surface exchange initialization rejects overflow" {
    var values = sequentialConcentrations();
    values.calcium_mol_per_megagram = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceExchangeInventory,
        initialize(2, values),
    );
}
