const std = @import("std");

/// Atmospheric-phase litter gas pools at litter layer index 0. Each field is
/// stored as grams of its named element (C, O, N, or H), not grams of molecule.
pub const AtmGasPools = struct {
    /// CO2G(0). CO2 in litter atmosphere (g C).
    co2_g: f64,
    /// CH4G(0). CH4 in litter atmosphere (g C).
    ch4_g: f64,
    /// OXYG(0). O2 in litter atmosphere (g O).
    o2_g: f64,
    /// Z2GG(0). N2 in litter atmosphere (g N).
    n2_g: f64,
    /// Z2OG(0). N2O in litter atmosphere (g N).
    n2o_g: f64,
    /// ZNH3G(0). NH3 in litter atmosphere (g N).
    nh3_g: f64,
    /// H2GG(0). H2 in litter atmosphere (g H).
    h2_g: f64,
};

/// Signed step-integrated increments from `trnsfr.f` and `nitro.f`. Field names
/// describe the legacy term; the signs below remain exactly those in REDIST.
pub const AtmGasFluxes = struct {
    /// XCOFLG(3,0). CO2 convective flux.
    co2_convective_g: f64,
    /// XCODFG(0). CO2 volatilization-dissolution.
    co2_dissolution_g: f64,
    /// RCGOX(0). CO2 from C combustion.
    co2_c_combustion_g: f64,
    /// RC4OX(0). CO2 from CH4 combustion.
    co2_ch4_combustion_g: f64,
    /// XCHFLG(3,0). CH4 convective flux.
    ch4_convective_g: f64,
    /// XCHDFG(0). CH4 volatilization-dissolution.
    ch4_dissolution_g: f64,
    /// RCHOX(0). CH4 from methanogenesis.
    ch4_production_g: f64,
    /// XOXFLG(3,0). O2 convective flux.
    o2_convective_g: f64,
    /// XOXDFG(0). O2 volatilization-dissolution.
    o2_dissolution_g: f64,
    /// ROGOX(0). O2 consumed by C combustion.
    o2_c_combustion_g: f64,
    /// XNGFLG(3,0). N2 convective flux.
    n2_convective_g: f64,
    /// XNGDFG(0). N2 volatilization-dissolution.
    n2_dissolution_g: f64,
    /// XN2FLG(3,0). N2O convective flux.
    n2o_convective_g: f64,
    /// XN2DFG(0). N2O volatilization-dissolution.
    n2o_dissolution_g: f64,
    /// XN3FLG(3,0). NH3 convective flux.
    nh3_convective_g: f64,
    /// XN3DFG(0). NH3 volatilization-dissolution.
    nh3_dissolution_g: f64,
    /// XHGFLG(3,0). H2 convective flux.
    h2_convective_g: f64,
    /// XHGDFG(0). H2 volatilization-dissolution.
    h2_dissolution_g: f64,
};

/// Aqueous-phase litter gas and nutrient pools at layer 0, stored as grams of
/// C, O, N, H, or P according to the field.
pub const AqueousPools = struct {
    /// CO2S(0). Aqueous CO2 (g C).
    co2_g: f64,
    /// CH4S(0). Aqueous CH4 (g C).
    ch4_g: f64,
    /// OXYS(0). Aqueous O2 (g O).
    o2_g: f64,
    /// Z2GS(0). Aqueous N2 (g N).
    n2_g: f64,
    /// Z2OS(0). Aqueous N2O (g N).
    n2o_g: f64,
    /// H2GS(0). Aqueous H2 (g H).
    h2_g: f64,
    /// ZNH4S(0). NH4 in litter (g N).
    nh4_g: f64,
    /// ZNH3S(0). NH3 in litter (g N).
    nh3_g: f64,
    /// ZNO3S(0). NO3 in litter (g N).
    no3_g: f64,
    /// ZNO2S(0). NO2 in litter (g N).
    no2_g: f64,
    /// H1PO4(0). HPO4 in litter (g P).
    hpo4_g: f64,
    /// H2PO4(0). H2PO4 in litter (g P).
    h2po4_g: f64,
};

/// Transformation and transport fluxes for aqueous/nutrient pools.
pub const AqueousFluxes = struct {
    /// XCODFR(0). CO2 litter volatilization from trnsfr.
    co2_litter_vol_g: f64,
    /// XCOFLS(3,0). CO2 solute convective-diffusive flux.
    co2_solute_g: f64,
    /// RCO2O(0). CO2 respiration consumption.
    co2_respiration_g: f64,
    /// TRCO2(0). CO2 transformation from solute.f.
    co2_transformation_g: f64,
    /// XCHDFR(0). CH4 litter volatilization.
    ch4_litter_vol_g: f64,
    /// XCHFLS(3,0). CH4 solute flux.
    ch4_solute_g: f64,
    /// RCH4O(0). CH4 respiration/methanotroph consumption.
    ch4_respiration_g: f64,
    /// XOXDFR(0). O2 litter volatilization.
    o2_litter_vol_g: f64,
    /// XOXFLS(3,0). O2 solute flux.
    o2_solute_g: f64,
    /// RUPOXO(0). O2 uptake.
    o2_uptake_g: f64,
    /// XNGDFR(0). N2 litter volatilization.
    n2_litter_vol_g: f64,
    /// XNGFLS(3,0). N2 solute flux.
    n2_solute_g: f64,
    /// RN2G(0). N2 fixation/production.
    n2_transformation_g: f64,
    /// XN2GS(0). N2 soil-surface transfer.
    n2_soil_transfer_g: f64,
    /// XN2DFR(0). N2O litter volatilization.
    n2o_litter_vol_g: f64,
    /// XN2FLS(3,0). N2O solute flux.
    n2o_solute_g: f64,
    /// RN2O(0). N2O transformation.
    n2o_transformation_g: f64,
    /// XHGDFR(0). H2 litter volatilization.
    h2_litter_vol_g: f64,
    /// XHGFLS(3,0). H2 solute flux.
    h2_solute_g: f64,
    /// RH2GO(0). H2 transformation.
    h2_transformation_g: f64,
    /// XN4FLW(3,0). NH4 water flux.
    nh4_water_flux_g: f64,
    /// XNH4S(0). NH4 soil-surface transfer.
    nh4_soil_transfer_g: f64,
    /// TRN4S(0). NH4 transformation.
    nh4_transformation_g: f64,
    /// XN3DFR(0). NH3 litter volatilization.
    nh3_litter_vol_g: f64,
    /// XN3FLW(3,0). NH3 water flux.
    nh3_water_flux_g: f64,
    /// XN3DFG(0). NH3 dissolution.
    nh3_dissolution_g: f64,
    /// TRN3S(0). NH3 transformation.
    nh3_transformation_g: f64,
    /// XNOFLW(3,0). NO3 water flux.
    no3_water_flux_g: f64,
    /// XNO3S(0). NO3 soil-surface transfer.
    no3_soil_transfer_g: f64,
    /// TRNO3(0). NO3 transformation.
    no3_transformation_g: f64,
    /// XNXFLS(3,0). NO2 solute flux.
    no2_solute_g: f64,
    /// XNO2S(0). NO2 soil-surface transfer.
    no2_soil_transfer_g: f64,
    /// TRH1P(0). HPO4 transformation.
    hpo4_transformation_g: f64,
    /// XH1PFS(3,0). HPO4 solute flux.
    hpo4_solute_g: f64,
    /// XH1PS(0). HPO4 soil-surface transfer.
    hpo4_soil_transfer_g: f64,
    /// TRH2P(0). H2PO4 transformation.
    h2po4_transformation_g: f64,
    /// XH2PFS(3,0). H2PO4 solute flux.
    h2po4_solute_g: f64,
    /// XH2PS(0). H2PO4 soil-surface transfer.
    h2po4_soil_transfer_g: f64,
};

pub const Result = struct {
    atm: AtmGasPools,
    aqueous: AqueousPools,
};

/// Direct translation of REDIST lines 4832--4876.
///
/// Updates litter-layer atmospheric and aqueous gas, nutrient pools.
/// 2.667 is the stoichiometric O:C ratio for CH4 combustion (O2 g per g C).
pub fn update(
    atm: AtmGasPools,
    aq: AqueousPools,
    af: AtmGasFluxes,
    aqf: AqueousFluxes,
) !Result {
    inline for (@typeInfo(AtmGasPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(atm, field.name)))
            return error.InvalidLitterGasPoolInput;
    inline for (@typeInfo(AqueousPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(aq, field.name)))
            return error.InvalidLitterGasPoolInput;
    inline for (@typeInfo(AtmGasFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(af, field.name)))
            return error.InvalidLitterGasFluxInput;
    inline for (@typeInfo(AqueousFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(aqf, field.name)))
            return error.InvalidLitterGasFluxInput;

    // Lines 4832-4833: CO2G
    const co2g = atm.co2_g + af.co2_convective_g - af.co2_dissolution_g +
        af.co2_c_combustion_g + af.co2_ch4_combustion_g;
    // Lines 4834-4835: CH4G
    const ch4g = atm.ch4_g + af.ch4_convective_g - af.ch4_dissolution_g +
        af.ch4_production_g - af.co2_ch4_combustion_g;
    // Lines 4836-4837: OXYG  (2.667 O per C combustion)
    const oxyg = atm.o2_g + af.o2_convective_g - af.o2_dissolution_g -
        af.o2_c_combustion_g - af.co2_ch4_combustion_g * 2.667;
    // Lines 4845-4846: Z2GG
    const z2gg = atm.n2_g + af.n2_convective_g - af.n2_dissolution_g;
    // Lines 4847-4848: Z2OG
    const z2og = atm.n2o_g + af.n2o_convective_g - af.n2o_dissolution_g;
    // Lines 4849-4850: ZNH3G
    const znh3g = atm.nh3_g + af.nh3_convective_g - af.nh3_dissolution_g;
    // Lines 4851-4852: H2GG
    const h2gg = atm.h2_g + af.h2_convective_g - af.h2_dissolution_g;

    // Lines 4853-4854: CO2S
    const co2s = aq.co2_g + aqf.co2_litter_vol_g + aqf.co2_solute_g +
        af.co2_dissolution_g - aqf.co2_respiration_g + aqf.co2_transformation_g;
    // Lines 4855-4856: CH4S
    const ch4s = aq.ch4_g + aqf.ch4_litter_vol_g + aqf.ch4_solute_g +
        af.ch4_dissolution_g - aqf.ch4_respiration_g;
    // Lines 4857-4858: OXYS
    const oxys = aq.o2_g + aqf.o2_litter_vol_g + aqf.o2_solute_g +
        af.o2_dissolution_g - aqf.o2_uptake_g;
    // Lines 4859-4860: Z2GS
    const z2gs = aq.n2_g + aqf.n2_litter_vol_g + aqf.n2_solute_g +
        af.n2_dissolution_g - aqf.n2_transformation_g - aqf.n2_soil_transfer_g;
    // Lines 4861-4862: Z2OS
    const z2os = aq.n2o_g + aqf.n2o_litter_vol_g + aqf.n2o_solute_g +
        af.n2o_dissolution_g - aqf.n2o_transformation_g;
    // Lines 4863-4864: H2GS
    const h2gs = aq.h2_g + aqf.h2_litter_vol_g + aqf.h2_solute_g +
        af.h2_dissolution_g - aqf.h2_transformation_g;
    // Lines 4865-4866: ZNH4S
    const znh4s = aq.nh4_g + aqf.nh4_water_flux_g + aqf.nh4_soil_transfer_g +
        aqf.nh4_transformation_g;
    // Lines 4867-4868: ZNH3S
    const znh3s = aq.nh3_g + aqf.nh3_litter_vol_g + aqf.nh3_water_flux_g +
        aqf.nh3_dissolution_g + aqf.nh3_transformation_g;
    // Lines 4869-4870: ZNO3S
    const zno3s = aq.no3_g + aqf.no3_water_flux_g + aqf.no3_soil_transfer_g +
        aqf.no3_transformation_g;
    // Lines 4871-4872: ZNO2S
    const zno2s = aq.no2_g + aqf.no2_solute_g + aqf.no2_soil_transfer_g;
    // Lines 4873-4874: H1PO4
    const hpo4 = aq.hpo4_g + aqf.hpo4_transformation_g + aqf.hpo4_solute_g +
        aqf.hpo4_soil_transfer_g;
    // Lines 4875-4876: H2PO4
    const h2po4 = aq.h2po4_g + aqf.h2po4_transformation_g + aqf.h2po4_solute_g +
        aqf.h2po4_soil_transfer_g;

    const result = Result{
        .atm = .{
            .co2_g = co2g,
            .ch4_g = ch4g,
            .o2_g = oxyg,
            .n2_g = z2gg,
            .n2o_g = z2og,
            .nh3_g = znh3g,
            .h2_g = h2gg,
        },
        .aqueous = .{
            .co2_g = co2s,
            .ch4_g = ch4s,
            .o2_g = oxys,
            .n2_g = z2gs,
            .n2o_g = z2os,
            .h2_g = h2gs,
            .nh4_g = znh4s,
            .nh3_g = znh3s,
            .no3_g = zno3s,
            .no2_g = zno2s,
            .hpo4_g = hpo4,
            .h2po4_g = h2po4,
        },
    };
    inline for (@typeInfo(AtmGasPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.atm, field.name)))
            return error.NonFiniteLitterGasPool;
    inline for (@typeInfo(AqueousPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.aqueous, field.name)))
            return error.NonFiniteLitterGasPool;
    return result;
}

fn zeroAtm() AtmGasPools {
    return std.mem.zeroes(AtmGasPools);
}

fn zeroAq() AqueousPools {
    return std.mem.zeroes(AqueousPools);
}

fn zeroAtmFlux() AtmGasFluxes {
    return std.mem.zeroes(AtmGasFluxes);
}

fn zeroAqFlux() AqueousFluxes {
    return std.mem.zeroes(AqueousFluxes);
}

test "REDIST litter atm CO2G adds convection and combustion subtracts dissolution" {
    var af = zeroAtmFlux();
    af.co2_convective_g = 1.0;
    af.co2_c_combustion_g = 0.5;
    af.co2_dissolution_g = 0.2;
    const result = try update(zeroAtm(), zeroAq(), af, zeroAqFlux());
    try std.testing.expectApproxEqRel(@as(f64, 1.3), result.atm.co2_g, 1.0e-15);
}

test "REDIST litter atm OXYG uses 2.667 for CH4 combustion O2 stoichiometry" {
    var af = zeroAtmFlux();
    af.co2_ch4_combustion_g = 1.0;
    const result = try update(zeroAtm(), zeroAq(), af, zeroAqFlux());
    // O2G = -RC4OX * 2.667
    try std.testing.expectApproxEqRel(@as(f64, -2.667), result.atm.o2_g, 1.0e-15);
}

test "REDIST litter aqueous CO2S adds dissolution from atm gas" {
    var af = zeroAtmFlux();
    af.co2_dissolution_g = 0.3;
    const result = try update(zeroAtm(), zeroAq(), af, zeroAqFlux());
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.aqueous.co2_g, 1.0e-15);
}

test "REDIST litter aqueous NH4 sums three sources" {
    var aqf = zeroAqFlux();
    aqf.nh4_water_flux_g = 1.0;
    aqf.nh4_soil_transfer_g = 2.0;
    aqf.nh4_transformation_g = 3.0;
    const result = try update(zeroAtm(), zeroAq(), zeroAtmFlux(), aqf);
    try std.testing.expectApproxEqRel(@as(f64, 6.0), result.aqueous.nh4_g, 1.0e-15);
}

test "REDIST litter gas update is identity for zero increments" {
    const atm = AtmGasPools{
        .co2_g = 1.0,
        .ch4_g = 2.0,
        .o2_g = 3.0,
        .n2_g = 4.0,
        .n2o_g = 5.0,
        .nh3_g = 6.0,
        .h2_g = 7.0,
    };
    const aq = AqueousPools{
        .co2_g = 8.0,
        .ch4_g = 9.0,
        .o2_g = 10.0,
        .n2_g = 11.0,
        .n2o_g = 12.0,
        .h2_g = 13.0,
        .nh4_g = 14.0,
        .nh3_g = 15.0,
        .no3_g = 16.0,
        .no2_g = 17.0,
        .hpo4_g = 18.0,
        .h2po4_g = 19.0,
    };
    const result = try update(atm, aq, zeroAtmFlux(), zeroAqFlux());
    try std.testing.expectEqualDeep(atm, result.atm);
    try std.testing.expectEqualDeep(aq, result.aqueous);
}

test "REDIST litter atmospheric minor gases preserve source signs" {
    var af = zeroAtmFlux();
    af.n2_convective_g = 8.0;
    af.n2_dissolution_g = 1.0;
    af.n2o_convective_g = 7.0;
    af.n2o_dissolution_g = 2.0;
    af.nh3_convective_g = 6.0;
    af.nh3_dissolution_g = 3.0;
    af.h2_convective_g = 5.0;
    af.h2_dissolution_g = 4.0;
    const result = try update(zeroAtm(), zeroAq(), af, zeroAqFlux());
    try std.testing.expectEqual(@as(f64, 7.0), result.atm.n2_g);
    try std.testing.expectEqual(@as(f64, 5.0), result.atm.n2o_g);
    try std.testing.expectEqual(@as(f64, 3.0), result.atm.nh3_g);
    try std.testing.expectEqual(@as(f64, 1.0), result.atm.h2_g);
}

test "REDIST litter aqueous nitrogen phosphorus and gases follow source order" {
    var aqf = zeroAqFlux();
    aqf.ch4_litter_vol_g = 1.0;
    aqf.ch4_solute_g = 2.0;
    aqf.ch4_respiration_g = 3.0;
    aqf.o2_litter_vol_g = 4.0;
    aqf.o2_solute_g = 5.0;
    aqf.o2_uptake_g = 6.0;
    aqf.n2_litter_vol_g = 7.0;
    aqf.n2_solute_g = 8.0;
    aqf.n2_transformation_g = 9.0;
    aqf.n2_soil_transfer_g = 10.0;
    aqf.n2o_litter_vol_g = 11.0;
    aqf.n2o_solute_g = 12.0;
    aqf.n2o_transformation_g = 13.0;
    aqf.h2_litter_vol_g = 14.0;
    aqf.h2_solute_g = 15.0;
    aqf.h2_transformation_g = 16.0;
    aqf.nh3_litter_vol_g = 17.0;
    aqf.nh3_water_flux_g = 18.0;
    aqf.nh3_transformation_g = 19.0;
    aqf.no3_water_flux_g = 20.0;
    aqf.no3_soil_transfer_g = 21.0;
    aqf.no3_transformation_g = 22.0;
    aqf.no2_solute_g = 23.0;
    aqf.no2_soil_transfer_g = 24.0;
    aqf.hpo4_transformation_g = 25.0;
    aqf.hpo4_solute_g = 26.0;
    aqf.hpo4_soil_transfer_g = 27.0;
    aqf.h2po4_transformation_g = 28.0;
    aqf.h2po4_solute_g = 29.0;
    aqf.h2po4_soil_transfer_g = 30.0;
    const result = try update(zeroAtm(), zeroAq(), zeroAtmFlux(), aqf);
    try std.testing.expectEqual(@as(f64, 0.0), result.aqueous.ch4_g);
    try std.testing.expectEqual(@as(f64, 3.0), result.aqueous.o2_g);
    try std.testing.expectEqual(@as(f64, -4.0), result.aqueous.n2_g);
    try std.testing.expectEqual(@as(f64, 10.0), result.aqueous.n2o_g);
    try std.testing.expectEqual(@as(f64, 13.0), result.aqueous.h2_g);
    try std.testing.expectEqual(@as(f64, 54.0), result.aqueous.nh3_g);
    try std.testing.expectEqual(@as(f64, 63.0), result.aqueous.no3_g);
    try std.testing.expectEqual(@as(f64, 47.0), result.aqueous.no2_g);
    try std.testing.expectEqual(@as(f64, 78.0), result.aqueous.hpo4_g);
    try std.testing.expectEqual(@as(f64, 87.0), result.aqueous.h2po4_g);
}

test "REDIST litter gas update rejects invalid input and overflow" {
    var bad_atm = zeroAtm();
    bad_atm.n2_g = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidLitterGasPoolInput,
        update(bad_atm, zeroAq(), zeroAtmFlux(), zeroAqFlux()),
    );

    var af = zeroAtmFlux();
    af.co2_convective_g = std.math.floatMax(f64);
    var atm = zeroAtm();
    atm.co2_g = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteLitterGasPool,
        update(atm, zeroAq(), af, zeroAqFlux()),
    );
}
