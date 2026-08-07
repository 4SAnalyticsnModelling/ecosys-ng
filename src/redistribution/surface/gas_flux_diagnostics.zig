const std = @import("std");

/// Step-integrated components used to assign the surface-litter gas fluxes.
/// Oxygen values are g O per step and methane values are g C per step.
pub const SurfaceComponents = struct {
    litter_volatilization_g: f64,
    surface_solute_flux_g: f64,
    gas_dissolution_g: f64,
    runoff_water_m3: f64,
    runoff_concentration_g_m3: f64,
    irrigation_water_m3: f64,
    irrigation_concentration_g_m3: f64,
};

pub const GasFluxDiagnostics = struct {
    /// ROXYL(0), assigned for the surface-litter layer (g O per step).
    litter_oxygen_g: f64,
    /// RCH4L(0), assigned for the surface-litter layer (g C per step).
    litter_methane_g: f64,
    /// ROXYL(NU), accumulated at the reference mineral layer (g O per step).
    mineral_oxygen_g: f64,
    /// RCH4L(NU), accumulated at the reference mineral layer (g C per step).
    mineral_methane_g: f64,
};

/// Direct translation of redist.f lines 4921--4926.
///
/// The litter diagnostics are assignments, while the mineral-layer diagnostics
/// retain their previous values and add the surface diffusion increments.
pub fn calculate(
    oxygen: SurfaceComponents,
    methane: SurfaceComponents,
    previous_mineral_oxygen_g: f64,
    previous_mineral_methane_g: f64,
    surface_diffusion_oxygen_g: f64,
    surface_diffusion_methane_g: f64,
) !GasFluxDiagnostics {
    inline for (@typeInfo(SurfaceComponents).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(oxygen, field.name)) or
            !std.math.isFinite(@field(methane, field.name)))
            return error.InvalidSurfaceGasFluxInput;
    }
    const scalars = [_]f64{
        previous_mineral_oxygen_g,
        previous_mineral_methane_g,
        surface_diffusion_oxygen_g,
        surface_diffusion_methane_g,
    };
    for (scalars) |value|
        if (!std.math.isFinite(value)) return error.InvalidSurfaceGasFluxInput;

    const result = GasFluxDiagnostics{
        .litter_oxygen_g = oxygen.litter_volatilization_g + oxygen.surface_solute_flux_g +
            oxygen.gas_dissolution_g -
            (oxygen.runoff_water_m3 * oxygen.runoff_concentration_g_m3 +
                oxygen.irrigation_water_m3 * oxygen.irrigation_concentration_g_m3),
        .litter_methane_g = methane.litter_volatilization_g + methane.surface_solute_flux_g +
            methane.gas_dissolution_g -
            (methane.runoff_water_m3 * methane.runoff_concentration_g_m3 +
                methane.irrigation_water_m3 * methane.irrigation_concentration_g_m3),
        .mineral_oxygen_g = previous_mineral_oxygen_g + surface_diffusion_oxygen_g,
        .mineral_methane_g = previous_mineral_methane_g + surface_diffusion_methane_g,
    };
    inline for (@typeInfo(GasFluxDiagnostics).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceGasFluxDiagnostic;
    return result;
}

fn zeroComponents() SurfaceComponents {
    return std.mem.zeroes(SurfaceComponents);
}

test "REDIST surface gas diagnostics preserve Fortran component order" {
    const oxygen = SurfaceComponents{
        .litter_volatilization_g = 20.0,
        .surface_solute_flux_g = 3.0,
        .gas_dissolution_g = 2.0,
        .runoff_water_m3 = 4.0,
        .runoff_concentration_g_m3 = 2.0,
        .irrigation_water_m3 = 3.0,
        .irrigation_concentration_g_m3 = 1.0,
    };
    const methane = SurfaceComponents{
        .litter_volatilization_g = 12.0,
        .surface_solute_flux_g = 4.0,
        .gas_dissolution_g = 1.0,
        .runoff_water_m3 = 2.0,
        .runoff_concentration_g_m3 = 3.0,
        .irrigation_water_m3 = 1.0,
        .irrigation_concentration_g_m3 = 2.0,
    };
    const result = try calculate(oxygen, methane, 10.0, 20.0, 0.5, -0.25);
    try std.testing.expectEqual(@as(f64, 14.0), result.litter_oxygen_g);
    try std.testing.expectEqual(@as(f64, 9.0), result.litter_methane_g);
    try std.testing.expectEqual(@as(f64, 10.5), result.mineral_oxygen_g);
    try std.testing.expectEqual(@as(f64, 19.75), result.mineral_methane_g);
}

test "REDIST surface gas diagnostics assign litter and accumulate mineral" {
    var oxygen = zeroComponents();
    oxygen.gas_dissolution_g = -2.0;
    var methane = zeroComponents();
    methane.surface_solute_flux_g = -3.0;
    const result = try calculate(oxygen, methane, 7.0, 8.0, -1.0, 2.0);
    try std.testing.expectEqual(@as(f64, -2.0), result.litter_oxygen_g);
    try std.testing.expectEqual(@as(f64, -3.0), result.litter_methane_g);
    try std.testing.expectEqual(@as(f64, 6.0), result.mineral_oxygen_g);
    try std.testing.expectEqual(@as(f64, 10.0), result.mineral_methane_g);
}

test "REDIST surface gas diagnostics reject invalid input" {
    var oxygen = zeroComponents();
    oxygen.runoff_water_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceGasFluxInput,
        calculate(oxygen, zeroComponents(), 0.0, 0.0, 0.0, 0.0),
    );
}

test "REDIST surface gas diagnostics reject arithmetic overflow" {
    var oxygen = zeroComponents();
    oxygen.litter_volatilization_g = std.math.floatMax(f64);
    oxygen.surface_solute_flux_g = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceGasFluxDiagnostic,
        calculate(oxygen, zeroComponents(), 0.0, 0.0, 0.0, 0.0),
    );
}
