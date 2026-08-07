const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

pub const GasInventories = struct {
    carbon_dioxide: f64,
    methane: f64,
    oxygen: f64,
    dinitrogen: f64,
    nitrous_oxide: f64,
    hydrogen: f64,
};

pub const MineralNutrientInventories = struct {
    ammonium_non_band: f64,
    ammonium_band: f64,
    ammonia_non_band: f64,
    ammonia_band: f64,
    nitrate_non_band: f64,
    nitrate_band: f64,
    nitrite_non_band: f64,
    nitrite_band: f64,
    dihydrogen_phosphate_non_band: f64,
    dihydrogen_phosphate_band: f64,
    hydrogen_phosphate_non_band: f64,
    hydrogen_phosphate_band: f64,
};

pub const PrimaryIonInventories = struct {
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
};

pub const ComplexInventories = struct {
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

pub const PhosphateInventories = struct {
    phosphate: f64,
    hydrogen_phosphate: f64,
    phosphoric_acid: f64,
    iron_hydrogen_phosphate: f64,
    iron_dihydrogen_phosphate: f64,
    calcium_phosphate: f64,
    calcium_hydrogen_phosphate: f64,
    calcium_dihydrogen_phosphate: f64,
    magnesium_hydrogen_phosphate: f64,
};

/// Macropore contents retain each source variable's units (`g` or `mol`).
pub const Inventories = struct {
    gases: GasInventories,
    mineral_nutrients: MineralNutrientInventories,
    primary_ions: PrimaryIonInventories,
    complexes: ComplexInventories,
    phosphate_non_band: PhosphateInventories,
    phosphate_band: PhosphateInventories,
};

fn zeroed(comptime T: type) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = 0.0;
    return result;
}

/// Direct translation of `starte.f` lines 1533--1601 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch.
pub fn initialize(control: Control) ?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    return .{
        .gases = zeroed(GasInventories),
        .mineral_nutrients = zeroed(MineralNutrientInventories),
        .primary_ions = zeroed(PrimaryIonInventories),
        .complexes = zeroed(ComplexInventories),
        .phosphate_non_band = zeroed(PhosphateInventories),
        .phosphate_band = zeroed(PhosphateInventories),
    };
}

fn expectAllZero(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(value, field.name));
}

test "STARTE macropore initialization resets every scientific inventory group" {
    const result = initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }).?;
    try expectAllZero(result.gases);
    try expectAllZero(result.mineral_nutrients);
    try expectAllZero(result.primary_ions);
    try expectAllZero(result.complexes);
    try expectAllZero(result.phosphate_non_band);
    try expectAllZero(result.phosphate_band);
}

test "STARTE supplied macropore profile bypasses zero initialization" {
    try std.testing.expectEqual(@as(?Inventories, null), initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }));
}

test "STARTE nonzero gas initialization index bypasses macropore reset" {
    try std.testing.expectEqual(@as(?Inventories, null), initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 1 }));
}
