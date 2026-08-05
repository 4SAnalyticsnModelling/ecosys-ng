const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

pub const NitrogenPools = struct {
    soil_ammonium_g_n: f64, // ZNH4S
    manure_organic_n_g_n: f64, // OSN(1,2,...)
};

/// Direct translation of STARTE.F lines 1693--1694 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch. The second assignment uses the
/// unchanged manure value exactly as in the sequential Fortran statements.
pub fn partition(control: Control, pools: NitrogenPools) !?NitrogenPools {
    if (control.profile_source != .atmospheric_equilibrium or control.gas_initialization_index != 0) return null;
    if (!std.math.isFinite(pools.soil_ammonium_g_n) or pools.soil_ammonium_g_n < 0 or
        !std.math.isFinite(pools.manure_organic_n_g_n) or pools.manure_organic_n_g_n < 0)
        return error.InvalidInitialNitrogenPool;
    const transfer_g_n = 0.5 * pools.manure_organic_n_g_n;
    const result: NitrogenPools = .{
        .soil_ammonium_g_n = pools.soil_ammonium_g_n + transfer_g_n,
        .manure_organic_n_g_n = pools.manure_organic_n_g_n - transfer_g_n,
    };
    if (!std.math.isFinite(result.soil_ammonium_g_n) or
        !std.math.isFinite(result.manure_organic_n_g_n))
        return error.InvalidInitialNitrogenPool;
    return result;
}

test "STARTE initial manure transfers exactly half into ammonium" {
    const before: NitrogenPools = .{ .soil_ammonium_g_n = 3, .manure_organic_n_g_n = 10 };
    const after = (try partition(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, before)).?;
    try std.testing.expectEqual(@as(f64, 8), after.soil_ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 5), after.manure_organic_n_g_n);
    try std.testing.expectEqual(before.soil_ammonium_g_n + before.manure_organic_n_g_n, after.soil_ammonium_g_n + after.manure_organic_n_g_n);
}

test "STARTE zero manure leaves ammonium unchanged" {
    const after = (try partition(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, .{ .soil_ammonium_g_n = 4, .manure_organic_n_g_n = 0 })).?;
    try std.testing.expectEqual(@as(f64, 4), after.soil_ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 0), after.manure_organic_n_g_n);
}

test "STARTE inactive manure partition ignores invalid dormant pools" {
    const nan = std.math.nan(f64);
    try std.testing.expectEqual(@as(?NitrogenPools, null), try partition(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, .{ .soil_ammonium_g_n = nan, .manure_organic_n_g_n = nan }));
}
