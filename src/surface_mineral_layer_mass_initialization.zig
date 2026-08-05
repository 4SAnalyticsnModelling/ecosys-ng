const std = @import("std");

pub const TextureFractions = struct {
    /// SAND(NU) in g Mg-1.
    sand_g_per_megagram: f64,
    /// SILT(NU) in g Mg-1.
    silt_g_per_megagram: f64,
    /// CLAY(NU) in g Mg-1.
    clay_g_per_megagram: f64,
};

/// Direct translation of starts.f lines 1646--1647.
///
/// Returns the total mineral mass of the reference surface layer NU
/// (legacy BKVLNM, g Mg-1) as the sum of sand, silt, and clay texture
/// fractions at that layer, clamped to zero (AMAX1(0.0, ...)).
///
/// The index NU is resolved by the caller; only the texture values at NU
/// are passed here.
pub fn calculate(texture: TextureFractions) !f64 {
    inline for (@typeInfo(TextureFractions).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(texture, field.name)))
            return error.NonFiniteMineralLayerMassInput;
    }

    const total = texture.sand_g_per_megagram +
        texture.silt_g_per_megagram +
        texture.clay_g_per_megagram;
    return @max(0.0, total);
}

test "STARTS surface mineral mass sums texture fractions" {
    const result = try calculate(.{
        .sand_g_per_megagram = 100.0,
        .silt_g_per_megagram = 200.0,
        .clay_g_per_megagram = 300.0,
    });
    try std.testing.expectEqual(@as(f64, 600.0), result);
}

test "STARTS surface mineral mass clamps negative sum to zero" {
    const result = try calculate(.{
        .sand_g_per_megagram = -100.0,
        .silt_g_per_megagram = -200.0,
        .clay_g_per_megagram = -300.0,
    });
    try std.testing.expectEqual(@as(f64, 0.0), result);
}

test "STARTS surface mineral mass rejects non-finite texture" {
    try std.testing.expectError(
        error.NonFiniteMineralLayerMassInput,
        calculate(.{
            .sand_g_per_megagram = std.math.nan(f64),
            .silt_g_per_megagram = 200.0,
            .clay_g_per_megagram = 300.0,
        }),
    );
    try std.testing.expectError(
        error.NonFiniteMineralLayerMassInput,
        calculate(.{
            .sand_g_per_megagram = 100.0,
            .silt_g_per_megagram = std.math.inf(f64),
            .clay_g_per_megagram = 300.0,
        }),
    );
}
