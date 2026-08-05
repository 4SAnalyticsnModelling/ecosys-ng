const std = @import("std");

/// Cell-scoped accumulators reset before traversing active soil layers.
pub const SoilLayerBalance = struct {
    /// TVHCP: total soil heat capacity (MJ K-1).
    total_heat_capacity_megajoules_k: f64,
    /// TVHCM: total soil-matrix heat capacity (MJ K-1).
    total_matrix_heat_capacity_megajoules_k: f64,
    /// TVOLW: total micropore liquid water (m3).
    total_matrix_water_m3: f64,
    /// TVOLV: total water-vapor equivalent (m3).
    total_water_vapor_m3: f64,
    /// TVOLWH: total macropore liquid water (m3).
    total_macropore_water_m3: f64,
    /// TVOLI: total micropore ice volume (m3).
    total_matrix_ice_m3: f64,
    /// TVOLIH: total macropore ice volume (m3).
    total_macropore_ice_m3: f64,
    /// TENGY: total soil energy (MJ).
    total_energy_megajoules: f64,
    /// TQALSI through TQKASI: primary-silicate inventories (mol).
    aluminum_silicate_mol: f64,
    iron_silicate_mol: f64,
    calcium_silicate_mol: f64,
    magnesium_silicate_mol: f64,
    sodium_silicate_mol: f64,
    potassium_silicate_mol: f64,
};

/// Exact translation of REDIST lines 5902--5915.
///
/// This creates fresh runtime state; it does not inspect overwritten legacy
/// storage. The caller subsequently traverses the runtime NU..NL layer range.
pub fn initialize() SoilLayerBalance {
    return .{
        .total_heat_capacity_megajoules_k = 0.0,
        .total_matrix_heat_capacity_megajoules_k = 0.0,
        .total_matrix_water_m3 = 0.0,
        .total_water_vapor_m3 = 0.0,
        .total_macropore_water_m3 = 0.0,
        .total_matrix_ice_m3 = 0.0,
        .total_macropore_ice_m3 = 0.0,
        .total_energy_megajoules = 0.0,
        .aluminum_silicate_mol = 0.0,
        .iron_silicate_mol = 0.0,
        .calcium_silicate_mol = 0.0,
        .magnesium_silicate_mol = 0.0,
        .sodium_silicate_mol = 0.0,
        .potassium_silicate_mol = 0.0,
    };
}

test "REDIST soil layer balance initializes all fourteen accumulators" {
    const result = initialize();
    comptime var field_count: usize = 0;
    inline for (@typeInfo(SoilLayerBalance).@"struct".fields) |field| {
        field_count += 1;
        try std.testing.expectEqual(@as(f64, 0.0), @field(result, field.name));
    }
    try std.testing.expectEqual(@as(usize, 14), field_count);
}

test "REDIST soil layer balance initialization is deterministic" {
    try std.testing.expectEqualDeep(initialize(), initialize());
}
