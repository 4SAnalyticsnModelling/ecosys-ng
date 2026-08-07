const std = @import("std");

pub const Inputs = struct {
    /// PRECQ. Rain plus snow precipitation rate (m3 h-1).
    precipitation_rain_snow_m3_per_h: f64,
    /// PRECI. Surface irrigation rate (m3 h-1).
    surface_irrigation_m3_per_h: f64,
    /// PRECU. Subsurface irrigation rate (m3 h-1).
    subsurface_irrigation_m3_per_h: f64,
    /// TEVAPG. Total soil evaporation from watsub.f (m3 step-1).
    soil_evaporation_m3: f64,
    /// TEVAPP. Total plant transpiration from watsub.f (m3 step-1).
    plant_transpiration_m3: f64,
    /// FLW(3,NK,NY,NX). Drainage flux through lower boundary (m3 step-1).
    lower_boundary_drainage_m3: f64,
    /// XNFH. Timestep (h step-1).
    timestep_h: f64,
};

/// Per-step water flux increments for one cell. Grid-level accumulators
/// (CRAIN, CEVAP, VOLWOU) receive the same increments; cell-level
/// accumulators (URAIN, UEVAP, HVOLO, UVOLO, UDRAIN, TEVPGH, TEVPPH) are
/// written by the production binding.
pub const Increments = struct {
    /// WI = (PRECQ+PRECI)*XNFH. Precipitation volume added to CRAIN/URAIN.
    precipitation_m3: f64,
    /// WO = TEVAPG+TEVAPP. Evapotranspiration subtracted from CEVAP/UEVAP.
    evapotranspiration_m3: f64,
    /// PRECU*XNFH. Subsurface irrigation subtracted from VOLWOU/HVOLO/UVOLO.
    subsurface_irrigation_m3: f64,
    /// FLW(3,NK). Drainage added to UDRAIN.
    lower_boundary_drainage_m3: f64,
    /// TEVAPG. Soil evaporation added to TEVPGH.
    soil_evaporation_m3: f64,
    /// TEVAPP. Plant transpiration added to TEVPPH.
    plant_transpiration_m3: f64,
};

/// Direct translation of redist.f lines 4405--4416.
///
/// Computes per-cell water flux increments for the surface boundary
/// water balance. The production binding applies these to the running cell
/// and grid accumulators, preserving source-order sign conventions:
/// precipitation is positive input, evapotranspiration and subsurface
/// irrigation are negative (subtracted from grid totals).
pub fn compute(inputs: Inputs) !Increments {
    if (!std.math.isFinite(inputs.timestep_h) or inputs.timestep_h <= 0)
        return error.InvalidSurfaceWaterFluxTimestep;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (comptime !std.mem.eql(u8, field.name, "timestep_h")) {
            if (!std.math.isFinite(@field(inputs, field.name)))
                return error.InvalidSurfaceWaterFluxInput;
        }
    }
    if (inputs.precipitation_rain_snow_m3_per_h < 0 or
        inputs.surface_irrigation_m3_per_h < 0 or
        inputs.subsurface_irrigation_m3_per_h < 0)
    {
        return error.NegativeSurfaceWaterInputRate;
    }

    const result = Increments{
        .precipitation_m3 = (inputs.precipitation_rain_snow_m3_per_h +
            inputs.surface_irrigation_m3_per_h) * inputs.timestep_h,
        .evapotranspiration_m3 = inputs.soil_evaporation_m3 +
            inputs.plant_transpiration_m3,
        .subsurface_irrigation_m3 = inputs.subsurface_irrigation_m3_per_h *
            inputs.timestep_h,
        .lower_boundary_drainage_m3 = inputs.lower_boundary_drainage_m3,
        .soil_evaporation_m3 = inputs.soil_evaporation_m3,
        .plant_transpiration_m3 = inputs.plant_transpiration_m3,
    };
    inline for (@typeInfo(Increments).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceWaterFluxIncrement;
    return result;
}

test "REDIST surface water flux preserves source-order precipitation scaling" {
    const result = try compute(.{
        .precipitation_rain_snow_m3_per_h = 0.002,
        .surface_irrigation_m3_per_h = 0.001,
        .subsurface_irrigation_m3_per_h = 0.0005,
        .soil_evaporation_m3 = 0.003,
        .plant_transpiration_m3 = 0.002,
        .lower_boundary_drainage_m3 = 0.001,
        .timestep_h = 1.0,
    });
    try std.testing.expectApproxEqRel(
        @as(f64, 0.003),
        result.precipitation_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        @as(f64, 0.005),
        result.evapotranspiration_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        @as(f64, 0.0005),
        result.subsurface_irrigation_m3,
        1.0e-15,
    );
}

test "REDIST surface water flux applies runtime timestep to rate inputs" {
    const result = try compute(.{
        .precipitation_rain_snow_m3_per_h = 1.0,
        .surface_irrigation_m3_per_h = 0.0,
        .subsurface_irrigation_m3_per_h = 0.5,
        .soil_evaporation_m3 = 0.0,
        .plant_transpiration_m3 = 0.0,
        .lower_boundary_drainage_m3 = 0.0,
        .timestep_h = 0.25,
    });
    try std.testing.expectEqual(@as(f64, 0.25), result.precipitation_m3);
    try std.testing.expectEqual(@as(f64, 0.125), result.subsurface_irrigation_m3);
}

test "REDIST surface water flux passes through per-step volumes unchanged" {
    const result = try compute(.{
        .precipitation_rain_snow_m3_per_h = 0.0,
        .surface_irrigation_m3_per_h = 0.0,
        .subsurface_irrigation_m3_per_h = 0.0,
        .soil_evaporation_m3 = 0.004,
        .plant_transpiration_m3 = 0.006,
        .lower_boundary_drainage_m3 = 0.002,
        .timestep_h = 1.0,
    });
    try std.testing.expectEqual(@as(f64, 0.004), result.soil_evaporation_m3);
    try std.testing.expectEqual(@as(f64, 0.006), result.plant_transpiration_m3);
    try std.testing.expectEqual(@as(f64, 0.002), result.lower_boundary_drainage_m3);
}

test "REDIST surface water flux rejects invalid timestep" {
    try std.testing.expectError(
        error.InvalidSurfaceWaterFluxTimestep,
        compute(.{
            .precipitation_rain_snow_m3_per_h = 0.0,
            .surface_irrigation_m3_per_h = 0.0,
            .subsurface_irrigation_m3_per_h = 0.0,
            .soil_evaporation_m3 = 0.0,
            .plant_transpiration_m3 = 0.0,
            .lower_boundary_drainage_m3 = 0.0,
            .timestep_h = 0.0,
        }),
    );
}

test "REDIST surface water flux rejects negative precipitation and irrigation rates" {
    const valid = Inputs{
        .precipitation_rain_snow_m3_per_h = 0.0,
        .surface_irrigation_m3_per_h = 0.0,
        .subsurface_irrigation_m3_per_h = 0.0,
        .soil_evaporation_m3 = 0.0,
        .plant_transpiration_m3 = 0.0,
        .lower_boundary_drainage_m3 = 0.0,
        .timestep_h = 1.0,
    };
    inline for (.{
        "precipitation_rain_snow_m3_per_h",
        "surface_irrigation_m3_per_h",
        "subsurface_irrigation_m3_per_h",
    }) |field_name| {
        var inputs = valid;
        @field(inputs, field_name) = -1.0e-6;
        try std.testing.expectError(
            error.NegativeSurfaceWaterInputRate,
            compute(inputs),
        );
    }
}
