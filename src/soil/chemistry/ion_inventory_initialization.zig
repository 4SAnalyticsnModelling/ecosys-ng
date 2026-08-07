const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

/// Concentrations are `mol m-3`; sulfate is `mol S m-3`.
pub const Concentrations = struct {
    aluminum: f64,
    iron: f64,
    hydrogen: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
    hydroxide: f64,
    sulfate: f64,
    chloride: f64,
    carbonate: f64,
    bicarbonate: f64,
    aluminum_monohydroxide: f64,
    aluminum_dihydroxide: f64,
    aluminum_trihydroxide: f64,
    aluminum_tetrahydroxide: f64,
    aluminum_sulfate: f64,
    iron_monohydroxide: f64,
    iron_dihydroxide: f64,
    iron_trihydroxide: f64,
    iron_tetrahydroxide: f64,
    iron_sulfate: f64,
    calcium_hydroxide: f64,
    calcium_carbonate: f64,
    calcium_bicarbonate: f64,
    calcium_sulfate: f64,
    magnesium_hydroxide: f64,
    magnesium_carbonate: f64,
    magnesium_bicarbonate: f64,
    magnesium_sulfate: f64,
    sodium_carbonate: f64,
    sodium_sulfate: f64,
    potassium_sulfate: f64,
};

/// Inventories are moles in the modeled layer; sulfate fields are mol S.
pub const Inventories = Concentrations;

/// Direct translation of `starte.f` lines 1469--1501. The identical field
/// layout intentionally makes the one-to-one concentration/inventory mapping
/// explicit while retaining source assignment order.
pub fn initialize(control: Control, concentrations: Concentrations, aqueous_volume_m3: f64) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    if (!std.math.isFinite(aqueous_volume_m3) or aqueous_volume_m3 < 0)
        return error.InvalidAqueousVolume;
    var result: Inventories = undefined;
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| {
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilIonConcentration;
        @field(result, field.name) = value * aqueous_volume_m3;
        if (!std.math.isFinite(@field(result, field.name)))
            return error.InvalidSoilIonInventory;
    }
    return result;
}

fn filled(value: f64) Concentrations {
    var result: Concentrations = undefined;
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE soil ion inventories preserve all source mappings" {
    var concentrations = filled(1);
    concentrations.aluminum = 2;
    concentrations.bicarbonate = 3;
    concentrations.aluminum_monohydroxide = 4;
    concentrations.iron_sulfate = 5;
    concentrations.potassium_sulfate = 6;
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, concentrations, 7)).?;
    try std.testing.expectEqual(@as(f64, 14), result.aluminum);
    try std.testing.expectEqual(@as(f64, 21), result.bicarbonate);
    try std.testing.expectEqual(@as(f64, 28), result.aluminum_monohydroxide);
    try std.testing.expectEqual(@as(f64, 35), result.iron_sulfate);
    try std.testing.expectEqual(@as(f64, 42), result.potassium_sulfate);
}

test "STARTE inactive soil ion initialization ignores dormant invalid input" {
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, filled(std.math.nan(f64)), std.math.nan(f64)));
}

test "STARTE soil ion initialization rejects late invalid concentration" {
    var concentrations = filled(1);
    concentrations.potassium_sulfate = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilIonConcentration, initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, concentrations, 1));
}
