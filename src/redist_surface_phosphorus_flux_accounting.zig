const std = @import("std");

/// Phosphate concentrations in precipitation and irrigation (mol P m-3).
/// READS converts user g P m-3 to this stored basis by dividing by 31.
pub const PhosphorusConcentrations = struct {
    /// CPOR. H2PO4 in rain (mol P m-3).
    h2po4_rain_mol_p_per_m3: f64,
    /// CH1PR. HPO4 in rain (mol P m-3).
    hpo4_rain_mol_p_per_m3: f64,
    /// CPOQ. H2PO4 in irrigation (mol P m-3).
    h2po4_irrigation_mol_p_per_m3: f64,
    /// CH1PQ. HPO4 in irrigation (mol P m-3).
    hpo4_irrigation_mol_p_per_m3: f64,
};

/// P fluxes through the lower (root zone bottom) boundary (g P step-1).
pub const PhosphorusLowerBoundary = struct {
    /// XH2PFS(3,NK). H2PO4 non-band flux.
    h2po4_nonband_g: f64,
    /// XH2BFB(3,NK). H2PO4 band flux.
    h2po4_band_g: f64,
    /// XH1PFS(3,NK). HPO4 non-band flux.
    hpo4_nonband_g: f64,
    /// XH1BFB(3,NK). HPO4 band flux.
    hpo4_band_g: f64,
};

pub const Inputs = struct {
    concentrations: PhosphorusConcentrations,
    lower_boundary: PhosphorusLowerBoundary,
    /// FLQGQ+FLQRQ. Combined rain water flux to surface (m3 step-1).
    rain_water_m3: f64,
    /// FLQGI+FLQRI. Combined irrigation water flux to surface (m3 step-1).
    irrigation_water_m3: f64,
    /// PRECU. Subsurface irrigation rate (m3 h-1).
    subsurface_irrigation_m3_per_h: f64,
    /// XNFH. Timestep (h step-1).
    timestep_h: f64,
};

pub const Increments = struct {
    /// PI. Aqueous P surface input (g P step-1). Added to TPIN.
    aqueous_input_g: f64,
    /// PXB. Signed subsurface P contribution (g P step-1). Subsurface
    /// irrigation carries the source-negative sign when added to TPOU.
    subsurface_signed_flux_g: f64,
    /// PDRAIN increment. P leaching through lower boundary (g P step-1).
    lower_boundary_drain_g: f64,
};

/// Direct translation of REDIST lines 4620--4629.
///
/// 31.0 g mol-1 is phosphorus atomic weight.
pub fn compute(inp: Inputs) !Increments {
    if (!std.math.isFinite(inp.timestep_h) or inp.timestep_h <= 0 or
        !std.math.isFinite(inp.rain_water_m3) or
        !std.math.isFinite(inp.irrigation_water_m3) or
        !std.math.isFinite(inp.subsurface_irrigation_m3_per_h))
        return error.InvalidSurfacePhosphorusFluxInput;
    inline for (@typeInfo(PhosphorusConcentrations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inp.concentrations, field.name)))
            return error.InvalidSurfacePhosphorusFluxInput;
    inline for (@typeInfo(PhosphorusLowerBoundary).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inp.lower_boundary, field.name)))
            return error.InvalidSurfacePhosphorusFluxInput;
    if (inp.rain_water_m3 < 0 or inp.irrigation_water_m3 < 0 or
        inp.subsurface_irrigation_m3_per_h < 0)
    {
        return error.NegativeSurfacePhosphorusWaterInput;
    }
    inline for (@typeInfo(PhosphorusConcentrations).@"struct".fields) |field| {
        if (@field(inp.concentrations, field.name) < 0)
            return error.NegativeSurfacePhosphorusConcentration;
    }

    const c = inp.concentrations;
    const lb = inp.lower_boundary;

    // Line 4620-4623.
    const pi = 31.0 * (inp.rain_water_m3 * (c.h2po4_rain_mol_p_per_m3 + c.hpo4_rain_mol_p_per_m3) +
        inp.irrigation_water_m3 * (c.h2po4_irrigation_mol_p_per_m3 + c.hpo4_irrigation_mol_p_per_m3));

    // Line 4624.
    const pxb = -31.0 * inp.subsurface_irrigation_m3_per_h *
        (c.h2po4_irrigation_mol_p_per_m3 + c.hpo4_irrigation_mol_p_per_m3) * inp.timestep_h;

    // Lines 4627-4629.
    const pdrain = lb.h2po4_nonband_g + lb.h2po4_band_g + lb.hpo4_nonband_g + lb.hpo4_band_g;

    const result = Increments{
        .aqueous_input_g = pi,
        .subsurface_signed_flux_g = pxb,
        .lower_boundary_drain_g = pdrain,
    };
    inline for (@typeInfo(Increments).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfacePhosphorusFluxIncrement;
    return result;
}

test "REDIST surface P flux aqueous input uses 31 g/mol" {
    const result = try compute(.{
        .concentrations = .{
            .h2po4_rain_mol_p_per_m3 = 1.0,
            .hpo4_rain_mol_p_per_m3 = 0.0,
            .h2po4_irrigation_mol_p_per_m3 = 0.0,
            .hpo4_irrigation_mol_p_per_m3 = 0.0,
        },
        .lower_boundary = std.mem.zeroes(PhosphorusLowerBoundary),
        .rain_water_m3 = 0.01,
        .irrigation_water_m3 = 0.0,
        .subsurface_irrigation_m3_per_h = 0.0,
        .timestep_h = 1.0,
    });
    try std.testing.expectApproxEqRel(@as(f64, 31.0 * 0.01 * 1.0), result.aqueous_input_g, 1.0e-15);
}

test "REDIST surface P flux subsurface output negative and timestep-scaled" {
    const result = try compute(.{
        .concentrations = .{
            .h2po4_rain_mol_p_per_m3 = 0.0,
            .hpo4_rain_mol_p_per_m3 = 0.0,
            .h2po4_irrigation_mol_p_per_m3 = 2.0,
            .hpo4_irrigation_mol_p_per_m3 = 0.0,
        },
        .lower_boundary = std.mem.zeroes(PhosphorusLowerBoundary),
        .rain_water_m3 = 0.0,
        .irrigation_water_m3 = 0.0,
        .subsurface_irrigation_m3_per_h = 0.001,
        .timestep_h = 2.0,
    });
    // PXB = -31.0 * 0.001 * 2.0 * 2.0
    try std.testing.expectApproxEqRel(@as(f64, -31.0 * 0.001 * 2.0 * 2.0), result.subsurface_signed_flux_g, 1.0e-15);
}

test "REDIST surface P flux lower boundary drain sums four pathways" {
    const result = try compute(.{
        .concentrations = std.mem.zeroes(PhosphorusConcentrations),
        .lower_boundary = .{
            .h2po4_nonband_g = 1.0,
            .h2po4_band_g = 2.0,
            .hpo4_nonband_g = 3.0,
            .hpo4_band_g = 4.0,
        },
        .rain_water_m3 = 0.0,
        .irrigation_water_m3 = 0.0,
        .subsurface_irrigation_m3_per_h = 0.0,
        .timestep_h = 1.0,
    });
    try std.testing.expectEqual(@as(f64, 10.0), result.lower_boundary_drain_g);
}

test "REDIST surface P flux rejects non-finite concentration" {
    try std.testing.expectError(
        error.InvalidSurfacePhosphorusFluxInput,
        compute(.{
            .concentrations = .{
                .h2po4_rain_mol_p_per_m3 = std.math.inf(f64),
                .hpo4_rain_mol_p_per_m3 = 0.0,
                .h2po4_irrigation_mol_p_per_m3 = 0.0,
                .hpo4_irrigation_mol_p_per_m3 = 0.0,
            },
            .lower_boundary = std.mem.zeroes(PhosphorusLowerBoundary),
            .rain_water_m3 = 0.0,
            .irrigation_water_m3 = 0.0,
            .subsurface_irrigation_m3_per_h = 0.0,
            .timestep_h = 1.0,
        }),
    );
}

test "REDIST surface P flux rejects negative phosphate concentration" {
    var concentrations = std.mem.zeroes(PhosphorusConcentrations);
    concentrations.hpo4_irrigation_mol_p_per_m3 = -1.0e-6;
    try std.testing.expectError(
        error.NegativeSurfacePhosphorusConcentration,
        compute(.{
            .concentrations = concentrations,
            .lower_boundary = std.mem.zeroes(PhosphorusLowerBoundary),
            .rain_water_m3 = 0,
            .irrigation_water_m3 = 0,
            .subsurface_irrigation_m3_per_h = 0,
            .timestep_h = 1,
        }),
    );
}
