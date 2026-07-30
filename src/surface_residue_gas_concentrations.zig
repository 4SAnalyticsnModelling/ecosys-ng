const std = @import("std");

pub const GasMasses = struct {
    carbon_dioxide_g: f64,
    methane_g: f64,
    oxygen_g: f64,
    nitrogen_g: f64,
    nitrous_oxide_g: f64,
    ammonia_g: f64,
    hydrogen_g: f64,
};

pub const GasConcentrations = struct {
    carbon_dioxide_g_m3: f64,
    methane_g_m3: f64,
    oxygen_g_m3: f64,
    nitrogen_g_m3: f64,
    nitrous_oxide_g_m3: f64,
    ammonia_g_m3: f64,
    hydrogen_g_m3: f64,
};

pub const ConcentrationError = error{
    NonFiniteInput,
    NegativeAirVolume,
    NonFiniteResult,
};

/// Translates HOUR1 lines 4576-4592 for the surface-residue pore space.
pub fn calculate(
    masses: GasMasses,
    air_filled_volume_m3: f64,
    air_volume_threshold_m3: f64,
) ConcentrationError!GasConcentrations {
    inline for (std.meta.fields(GasMasses)) |field| {
        if (!std.math.isFinite(@field(masses, field.name))) return error.NonFiniteInput;
    }
    if (!std.math.isFinite(air_filled_volume_m3) or
        !std.math.isFinite(air_volume_threshold_m3))
    {
        return error.NonFiniteInput;
    }
    if (air_filled_volume_m3 < 0.0 or air_volume_threshold_m3 < 0.0) {
        return error.NegativeAirVolume;
    }
    if (air_filled_volume_m3 <= air_volume_threshold_m3) {
        return std.mem.zeroes(GasConcentrations);
    }

    const concentrations = GasConcentrations{
        .carbon_dioxide_g_m3 = @max(0.0, masses.carbon_dioxide_g / air_filled_volume_m3),
        .methane_g_m3 = @max(0.0, masses.methane_g / air_filled_volume_m3),
        .oxygen_g_m3 = @max(0.0, masses.oxygen_g / air_filled_volume_m3),
        .nitrogen_g_m3 = @max(0.0, masses.nitrogen_g / air_filled_volume_m3),
        .nitrous_oxide_g_m3 = @max(0.0, masses.nitrous_oxide_g / air_filled_volume_m3),
        .ammonia_g_m3 = @max(0.0, masses.ammonia_g / air_filled_volume_m3),
        .hydrogen_g_m3 = @max(0.0, masses.hydrogen_g / air_filled_volume_m3),
    };
    inline for (std.meta.fields(GasConcentrations)) |field| {
        if (!std.math.isFinite(@field(concentrations, field.name))) {
            return error.NonFiniteResult;
        }
    }
    return concentrations;
}

test "surface-residue gas masses divide by air-filled volume" {
    const concentrations = try calculate(.{
        .carbon_dioxide_g = 4.0,
        .methane_g = 2.0,
        .oxygen_g = 8.0,
        .nitrogen_g = 10.0,
        .nitrous_oxide_g = 1.0,
        .ammonia_g = -2.0,
        .hydrogen_g = 0.5,
    }, 2.0, 0.0);
    try std.testing.expectEqual(@as(f64, 2.0), concentrations.carbon_dioxide_g_m3);
    try std.testing.expectEqual(@as(f64, 0.0), concentrations.ammonia_g_m3);
    try std.testing.expectEqual(@as(f64, 0.25), concentrations.hydrogen_g_m3);
}

test "air volume at threshold returns zero concentrations" {
    const concentrations = try calculate(
        std.mem.zeroes(GasMasses),
        0.0,
        0.0,
    );
    try std.testing.expectEqual(std.mem.zeroes(GasConcentrations), concentrations);
}
