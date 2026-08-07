const std = @import("std");

pub const MixingControl = struct {
    /// Fraction of material redistributed by the disturbance.
    mixed_fraction: f64,
    /// Fraction retained in its original layer (legacy XCORP).
    retained_fraction: f64,
    includes_all_plant_species: bool,
    excludes_primary_plant_species: bool,
};

/// Exact DAY disturbance-code interpretation from day.f:332-359.
///
/// Codes 0...10 mix every plant species. Codes 11...20 retain the
/// primary species while mixing all others. Every other code disables
/// redistribution. The source's minimum active mixing fraction of 0.05 is
/// preserved, including for code zero.
pub fn derive(disturbance_code: i32) MixingControl {
    if (disturbance_code >= 0 and disturbance_code <= 10) {
        const mixed_fraction = std.math.clamp(
            @as(f64, @floatFromInt(disturbance_code)) / 10.0,
            0.05,
            1.0,
        );
        return .{
            .mixed_fraction = mixed_fraction,
            .retained_fraction = 1.0 - mixed_fraction,
            .includes_all_plant_species = true,
            .excludes_primary_plant_species = false,
        };
    }
    if (disturbance_code >= 11 and disturbance_code <= 20) {
        const mixed_fraction = std.math.clamp(
            @as(f64, @floatFromInt(disturbance_code - 10)) / 10.0,
            0.05,
            1.0,
        );
        return .{
            .mixed_fraction = mixed_fraction,
            .retained_fraction = 1.0 - mixed_fraction,
            .includes_all_plant_species = false,
            .excludes_primary_plant_species = true,
        };
    }
    return .{
        .mixed_fraction = 0.0,
        .retained_fraction = 1.0,
        .includes_all_plant_species = false,
        .excludes_primary_plant_species = false,
    };
}

test "zero disturbance code preserves source minimum active mixing" {
    const result = derive(0);
    try std.testing.expectEqual(@as(f64, 0.05), result.mixed_fraction);
    try std.testing.expectEqual(@as(f64, 0.95), result.retained_fraction);
    try std.testing.expect(result.includes_all_plant_species);
    try std.testing.expect(!result.excludes_primary_plant_species);
}

test "all-species disturbance codes scale mixing through upper boundary" {
    var result = derive(7);
    try std.testing.expectEqual(@as(f64, 0.7), result.mixed_fraction);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.retained_fraction, 1e-15);

    result = derive(10);
    try std.testing.expectEqual(@as(f64, 1.0), result.mixed_fraction);
    try std.testing.expectEqual(@as(f64, 0.0), result.retained_fraction);
}

test "primary-species exclusion subtracts the exact code offset" {
    var result = derive(11);
    try std.testing.expectEqual(@as(f64, 0.1), result.mixed_fraction);
    try std.testing.expectEqual(@as(f64, 0.9), result.retained_fraction);
    try std.testing.expect(!result.includes_all_plant_species);
    try std.testing.expect(result.excludes_primary_plant_species);

    result = derive(20);
    try std.testing.expectEqual(@as(f64, 1.0), result.mixed_fraction);
}

test "non-tillage disturbance codes disable layer redistribution" {
    inline for (.{ -1, 21, 22, 24, 100 }) |code| {
        const result = derive(code);
        try std.testing.expectEqual(@as(f64, 0.0), result.mixed_fraction);
        try std.testing.expectEqual(@as(f64, 1.0), result.retained_fraction);
        try std.testing.expect(!result.includes_all_plant_species);
        try std.testing.expect(!result.excludes_primary_plant_species);
    }
}
