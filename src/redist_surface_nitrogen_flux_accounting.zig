const std = @import("std");

/// Nitrogen concentrations in precipitation and irrigation. Dissolved ionic
/// species are stored as mol N m-3 (legacy READS divides user g N m-3 by 14);
/// dissolved gases retain g N m-3.
pub const NitrogenConcentrations = struct {
    /// CN4R. NH4 in rain (mol N m-3).
    nh4_rain_mol_n_per_m3: f64,
    /// CN3R. NH3 in rain (mol N m-3).
    nh3_rain_mol_n_per_m3: f64,
    /// CNOR. NO3 in rain (mol N m-3).
    no3_rain_mol_n_per_m3: f64,
    /// CNNR. N2 in rain (g N m-3).
    n2_rain_g_per_m3: f64,
    /// CN2R. N2O in rain (g N m-3).
    n2o_rain_g_per_m3: f64,
    /// CN4Q. NH4 in irrigation (mol N m-3).
    nh4_irrigation_mol_n_per_m3: f64,
    /// CN3Q. NH3 in irrigation (mol N m-3).
    nh3_irrigation_mol_n_per_m3: f64,
    /// CNOQ. NO3 in irrigation (mol N m-3).
    no3_irrigation_mol_n_per_m3: f64,
    /// CNNQ. N2 in irrigation (g N m-3).
    n2_irrigation_g_per_m3: f64,
    /// CN2Q. N2O in irrigation (g N m-3).
    n2o_irrigation_g_per_m3: f64,
};

/// Gas N fluxes from trnsfr.f (g N step-1).
pub const NitrogenGasFluxes = struct {
    /// XNGDFS. N2 soil surface diffusion.
    n2_soil_diffusion_g: f64,
    /// XN2DFS. N2O soil surface diffusion.
    n2o_soil_diffusion_g: f64,
    /// XN3DFS. NH3 soil surface diffusion.
    nh3_soil_diffusion_g: f64,
    /// XNBDFS. NH4 soil surface diffusion (converted to gas).
    nh4_soil_diffusion_g: f64,
    /// XNGFLG(3,NU). N2 convective flux at mineral reference layer.
    n2_mineral_convective_g: f64,
    /// XN2FLG(3,NU). N2O convective flux at mineral reference layer.
    n2o_mineral_convective_g: f64,
    /// XN3FLG(3,NU). NH3 convective flux at mineral reference layer.
    nh3_mineral_convective_g: f64,
    /// TN2OZ. N2O bulk mass transfer.
    n2o_bulk_transfer_g: f64,
    /// TNH3Z. NH3 bulk mass transfer.
    nh3_bulk_transfer_g: f64,
    /// XN2FLG(3,0). N2O convective flux at litter.
    n2o_litter_convective_g: f64,
    /// XNGFLG(3,0). N2 convective flux at litter.
    n2_litter_convective_g: f64,
    /// XN3FLG(3,0). NH3 convective flux at litter.
    nh3_litter_convective_g: f64,
    /// XNGDFR. N2 litter volatilization.
    n2_litter_volatilization_g: f64,
    /// XN2DFR. N2O litter volatilization.
    n2o_litter_volatilization_g: f64,
    /// XN3DFR. NH3 litter volatilization.
    nh3_litter_volatilization_g: f64,
    /// XNGDFG(0). N2 soil-surface gas dissolution.
    n2_surface_dissolution_g: f64,
    /// XN2DFG(0). N2O soil-surface gas dissolution.
    n2o_surface_dissolution_g: f64,
    /// XN3DFG(0). NH3 soil-surface gas dissolution.
    nh3_surface_dissolution_g: f64,
};

/// N fluxes through the lower (root zone bottom) boundary (g N step-1).
pub const NitrogenLowerBoundary = struct {
    /// XN4FLW(3,NK). NH4 non-band water flux.
    nh4_water_g: f64,
    /// XN3FLW(3,NK). NH3 non-band water flux.
    nh3_water_g: f64,
    /// XNOFLW(3,NK). NO3 non-band water flux.
    no3_water_g: f64,
    /// XNXFLS(3,NK). NO2 non-band solute flux.
    no2_water_g: f64,
    /// XN4FLB(3,NK). NH4 band water flux.
    nh4_band_g: f64,
    /// XN3FLB(3,NK). NH3 band water flux.
    nh3_band_g: f64,
    /// XNOFLB(3,NK). NO3 band water flux.
    no3_band_g: f64,
    /// XNXFLB(3,NK). NO2 band flux.
    no2_band_g: f64,
};

pub const Inputs = struct {
    concentrations: NitrogenConcentrations,
    gas_fluxes: NitrogenGasFluxes,
    lower_boundary: NitrogenLowerBoundary,
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
    /// ZSI. Aqueous N input (g N step-1). Added to TZIN.
    aqueous_input_g: f64,
    /// ZXB. Signed subsurface N contribution (g N step-1). Irrigation input
    /// carries the source-negative sign when added to TZOU.
    subsurface_signed_flux_g: f64,
    /// ZGI. Gas N surface input (g N step-1). Added to ZN2GIN.
    gas_input_g: f64,
    /// ZDRAIN increment. N leaching through lower boundary (g N step-1).
    lower_boundary_drain_g: f64,
    /// ZNGGIN. N2 gas exchange (g N step-1).
    n2_gas_exchange_g: f64,
    /// ZN2OIN. N2O gas exchange (g N step-1).
    n2o_gas_exchange_g: f64,
    /// ZNH3IN. NH3 gas exchange (g N step-1).
    nh3_gas_exchange_g: f64,
};

/// Direct translation of REDIST lines 4558--4591.
///
/// Computes surface boundary N flux increments in source order.
pub fn compute(inp: Inputs) !Increments {
    if (!std.math.isFinite(inp.timestep_h) or inp.timestep_h <= 0 or
        !std.math.isFinite(inp.rain_water_m3) or
        !std.math.isFinite(inp.irrigation_water_m3) or
        !std.math.isFinite(inp.subsurface_irrigation_m3_per_h))
        return error.InvalidSurfaceNitrogenFluxInput;
    inline for (@typeInfo(NitrogenConcentrations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inp.concentrations, field.name)))
            return error.InvalidSurfaceNitrogenFluxInput;
    inline for (@typeInfo(NitrogenGasFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inp.gas_fluxes, field.name)))
            return error.InvalidSurfaceNitrogenFluxInput;
    inline for (@typeInfo(NitrogenLowerBoundary).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inp.lower_boundary, field.name)))
            return error.InvalidSurfaceNitrogenFluxInput;
    if (inp.rain_water_m3 < 0 or inp.irrigation_water_m3 < 0 or
        inp.subsurface_irrigation_m3_per_h < 0)
    {
        return error.NegativeSurfaceNitrogenWaterInput;
    }
    inline for (@typeInfo(NitrogenConcentrations).@"struct".fields) |field| {
        if (@field(inp.concentrations, field.name) < 0)
            return error.NegativeSurfaceNitrogenConcentration;
    }

    const c = inp.concentrations;
    const gf = inp.gas_fluxes;
    const lb = inp.lower_boundary;

    // Line 4558-4561: aqueous N input (g N).
    const zsi = (inp.rain_water_m3 * (c.nh4_rain_mol_n_per_m3 + c.nh3_rain_mol_n_per_m3 + c.no3_rain_mol_n_per_m3) +
        inp.irrigation_water_m3 * (c.nh4_irrigation_mol_n_per_m3 + c.nh3_irrigation_mol_n_per_m3 + c.no3_irrigation_mol_n_per_m3)) * 14.0;

    // Line 4562-4563: subsurface N output.
    const zxb = (-inp.subsurface_irrigation_m3_per_h * (c.n2_irrigation_g_per_m3 + c.n2o_irrigation_g_per_m3) -
        inp.subsurface_irrigation_m3_per_h *
            (c.nh4_irrigation_mol_n_per_m3 + c.nh3_irrigation_mol_n_per_m3 + c.no3_irrigation_mol_n_per_m3) * 14.0) * inp.timestep_h;

    // Lines 4566-4573: gas N input.
    const zgi = inp.rain_water_m3 * (c.n2_rain_g_per_m3 + c.n2o_rain_g_per_m3) +
        inp.irrigation_water_m3 * (c.n2_irrigation_g_per_m3 + c.n2o_irrigation_g_per_m3) +
        gf.n2_soil_diffusion_g + gf.n2o_soil_diffusion_g + gf.nh3_soil_diffusion_g +
        gf.nh4_soil_diffusion_g +
        gf.n2_mineral_convective_g + gf.n2o_mineral_convective_g + gf.nh3_mineral_convective_g +
        gf.n2o_bulk_transfer_g + gf.nh3_bulk_transfer_g +
        gf.n2o_litter_convective_g + gf.n2_litter_convective_g + gf.nh3_litter_convective_g +
        gf.n2_litter_volatilization_g + gf.n2o_litter_volatilization_g + gf.nh3_litter_volatilization_g;

    // Lines 4575-4579: N leaching through lower boundary.
    const zdrain = lb.nh4_water_g + lb.nh3_water_g + lb.no3_water_g + lb.no2_water_g +
        lb.nh4_band_g + lb.nh3_band_g + lb.no3_band_g + lb.no2_band_g;

    // Lines 4580-4583: individual gas species.
    const znggin = gf.n2_soil_diffusion_g + gf.n2_mineral_convective_g + gf.n2_surface_dissolution_g;
    const zn2oin = gf.n2o_soil_diffusion_g + gf.n2o_mineral_convective_g + gf.n2o_surface_dissolution_g;
    const znh3in = gf.nh3_soil_diffusion_g + gf.nh4_soil_diffusion_g +
        gf.nh3_mineral_convective_g + gf.nh3_surface_dissolution_g;

    const result = Increments{
        .aqueous_input_g = zsi,
        .subsurface_signed_flux_g = zxb,
        .gas_input_g = zgi,
        .lower_boundary_drain_g = zdrain,
        .n2_gas_exchange_g = znggin,
        .n2o_gas_exchange_g = zn2oin,
        .nh3_gas_exchange_g = znh3in,
    };
    inline for (@typeInfo(Increments).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceNitrogenFluxIncrement;
    return result;
}

fn zeroConcentrations() NitrogenConcentrations {
    return std.mem.zeroes(NitrogenConcentrations);
}

fn zeroGasFluxes() NitrogenGasFluxes {
    return std.mem.zeroes(NitrogenGasFluxes);
}

fn zeroLowerBoundary() NitrogenLowerBoundary {
    return std.mem.zeroes(NitrogenLowerBoundary);
}

test "REDIST surface N flux aqueous input uses 14 g/mol conversion" {
    var conc = zeroConcentrations();
    conc.nh4_rain_mol_n_per_m3 = 1.0;
    const result = try compute(.{
        .concentrations = conc,
        .gas_fluxes = zeroGasFluxes(),
        .lower_boundary = zeroLowerBoundary(),
        .rain_water_m3 = 0.01,
        .irrigation_water_m3 = 0.0,
        .subsurface_irrigation_m3_per_h = 0.0,
        .timestep_h = 1.0,
    });
    try std.testing.expectApproxEqRel(@as(f64, 0.01 * 1.0 * 14.0), result.aqueous_input_g, 1.0e-15);
}

test "REDIST surface N flux subsurface output applies timestep" {
    var conc = zeroConcentrations();
    conc.nh4_irrigation_mol_n_per_m3 = 1.0;
    const result = try compute(.{
        .concentrations = conc,
        .gas_fluxes = zeroGasFluxes(),
        .lower_boundary = zeroLowerBoundary(),
        .rain_water_m3 = 0.0,
        .irrigation_water_m3 = 0.0,
        .subsurface_irrigation_m3_per_h = 0.001,
        .timestep_h = 0.5,
    });
    // ZXB = -(0.001 * 1.0 * 14.0) * 0.5
    try std.testing.expectApproxEqRel(
        @as(f64, -0.001 * 1.0 * 14.0 * 0.5),
        result.subsurface_signed_flux_g,
        1.0e-15,
    );
}

test "REDIST surface N flux lower boundary drain sums all N pathways" {
    var lb = zeroLowerBoundary();
    lb.nh4_water_g = 1.0;
    lb.nh3_water_g = 2.0;
    lb.no3_band_g = 3.0;
    const result = try compute(.{
        .concentrations = zeroConcentrations(),
        .gas_fluxes = zeroGasFluxes(),
        .lower_boundary = lb,
        .rain_water_m3 = 0.0,
        .irrigation_water_m3 = 0.0,
        .subsurface_irrigation_m3_per_h = 0.0,
        .timestep_h = 1.0,
    });
    try std.testing.expectEqual(@as(f64, 6.0), result.lower_boundary_drain_g);
}

test "REDIST surface N flux rejects negative aqueous concentration" {
    var conc = zeroConcentrations();
    conc.no3_rain_mol_n_per_m3 = -1.0e-6;
    try std.testing.expectError(
        error.NegativeSurfaceNitrogenConcentration,
        compute(.{
            .concentrations = conc,
            .gas_fluxes = zeroGasFluxes(),
            .lower_boundary = zeroLowerBoundary(),
            .rain_water_m3 = 0,
            .irrigation_water_m3 = 0,
            .subsurface_irrigation_m3_per_h = 0,
            .timestep_h = 1,
        }),
    );
}
