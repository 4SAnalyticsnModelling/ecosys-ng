const std = @import("std");

/// Valence-weighted ion concentrations in precipitation (mol m-3).
///
/// READS converts user mass concentrations to molar concentrations before
/// these values enter STARTE speciation. Historical field suffixes are kept
/// locally for source traceability; every concentration field has mol m-3
/// semantics.
/// Species are organized by the coefficient group in the REDIST SIR formula
/// (lines 4724--4735), preserving source order exactly.
pub const IonConcentrationsRain = struct {
    // --- coefficient 1 ---
    /// CALR. Al3+ in rain.
    al_g_per_m3: f64,
    /// CFER. Fe3+ in rain.
    fe_g_per_m3: f64,
    /// CHYR. H+ in rain.
    hy_g_per_m3: f64,
    /// CCAR. Ca2+ in rain.
    ca_g_per_m3: f64,
    /// CMGR. Mg2+ in rain.
    megagrams_g_per_m3: f64,
    /// CNAR. Na+ in rain.
    na_g_per_m3: f64,
    /// CKAR. K+ in rain.
    ka_g_per_m3: f64,
    /// COHR. OH- in rain.
    oh_g_per_m3: f64,
    /// CSOR. SO42- in rain.
    so4_g_per_m3: f64,
    /// CCLR. Cl- in rain.
    cl_g_per_m3: f64,
    /// CC3R. CO32- (or HCO3-) in rain.
    co3_g_per_m3: f64,
    /// CH0PR. PO43- in rain.
    h0p_g_per_m3: f64,
    // --- coefficient 2 ---
    /// CHCR. HCO3- in rain.
    hco3_g_per_m3: f64,
    /// CAL1R. AlOH2+ in rain.
    al1_g_per_m3: f64,
    /// CALSR. AlSO4+ in rain.
    als_g_per_m3: f64,
    /// CFE1R. FeOH2+ in rain.
    fe1_g_per_m3: f64,
    /// CFESR. FeSO4+ in rain.
    fes_g_per_m3: f64,
    /// CCAOR. CaOH+ in rain.
    cao_g_per_m3: f64,
    /// CCACR. CaCO3 in rain.
    cac_g_per_m3: f64,
    /// CCASR. CaSO4 in rain.
    cas_g_per_m3: f64,
    /// CMGOR. MgOH+ in rain.
    mgo_g_per_m3: f64,
    /// CMGCR. MgCO3 in rain.
    mgc_g_per_m3: f64,
    /// CMGSR. MgSO4 in rain.
    mgs_g_per_m3: f64,
    /// CNACR. NaCO3- in rain.
    nac_g_per_m3: f64,
    /// CNASR. NaSO4- in rain.
    nas_g_per_m3: f64,
    /// CKASR. KSO4- in rain.
    kas_g_per_m3: f64,
    /// CC0PR. CaPO4- in rain.
    c0p_g_per_m3: f64,
    // --- coefficient 3 ---
    /// CAL2R. AlOH2+ in rain (ALOH2 species).
    al2_g_per_m3: f64,
    /// CFE2R. FeOH2+ in rain (FEOH2 species).
    fe2_g_per_m3: f64,
    /// CCAHR. CaHCO3+ in rain.
    cah_g_per_m3: f64,
    /// CMGHR. MgHCO3+ in rain.
    mgh_g_per_m3: f64,
    /// CF1PR. FeHPO42- in rain.
    f1p_g_per_m3: f64,
    /// CC1PR. CaHPO4 in rain.
    c1p_g_per_m3: f64,
    /// CM1PR. MgHPO4 in rain.
    m1p_g_per_m3: f64,
    // --- coefficient 4 ---
    /// CAL3R. AlOH3 in rain.
    al3_g_per_m3: f64,
    /// CFE3R. FeOH3 in rain.
    fe3_g_per_m3: f64,
    /// CH3PR. H3PO4 in rain.
    h3p_g_per_m3: f64,
    /// CF2PR. FeH2PO4- in rain.
    f2p_g_per_m3: f64,
    /// CC2PR. CaH2PO4+ in rain.
    c2p_g_per_m3: f64,
    // --- coefficient 5 ---
    /// CAL4R. AlOH4- in rain.
    al4_g_per_m3: f64,
    /// CFE4R. FeOH4- in rain.
    fe4_g_per_m3: f64,
    // Phosphate species for PIR reuse the fields already declared above:
    // h0p (PO43-), h3p (H3PO4), f1p (FeHPO42-), f2p (FeH2PO4-),
    // c0p (CaPO4-), c1p (CaHPO4), c2p (CaH2PO4+), m1p (MgHPO4).
};

/// Same 35+ species for irrigation (indexed by I). Mirrors IonConcentrationsRain
/// but drawn from irrigation concentrations C*Q(I,...).
pub const IonConcentrationsIrr = struct {
    // coefficient 1
    al_g_per_m3: f64, // CALQ
    fe_g_per_m3: f64, // CFEQ
    hy_g_per_m3: f64, // CHYQ
    ca_g_per_m3: f64, // CCAQ
    megagrams_g_per_m3: f64, // CMGQ
    na_g_per_m3: f64, // CNAQ
    ka_g_per_m3: f64, // CKAQ
    oh_g_per_m3: f64, // COHQ
    so4_g_per_m3: f64, // CSOQ
    cl_g_per_m3: f64, // CCLQ
    co3_g_per_m3: f64, // CC3Q
    h0p_g_per_m3: f64, // CH0PQ
    // coefficient 2
    hco3_g_per_m3: f64, // CHCQ
    al1_g_per_m3: f64, // CAL1Q
    als_g_per_m3: f64, // CALSQ
    fe1_g_per_m3: f64, // CFE1Q
    fes_g_per_m3: f64, // CFESQ
    cao_g_per_m3: f64, // CCAOQ
    cac_g_per_m3: f64, // CCACQ
    cas_g_per_m3: f64, // CCASQ
    mgo_g_per_m3: f64, // CMGOQ
    mgc_g_per_m3: f64, // CMGCQ
    mgs_g_per_m3: f64, // CMGSQ
    nac_g_per_m3: f64, // CNACQ
    nas_g_per_m3: f64, // CNASQ
    kas_g_per_m3: f64, // CKASQ
    c0p_g_per_m3: f64, // CC0PQ
    // coefficient 3
    al2_g_per_m3: f64, // CAL2Q
    fe2_g_per_m3: f64, // CFE2Q
    cah_g_per_m3: f64, // CCAHQ
    mgh_g_per_m3: f64, // CMGHQ
    f1p_g_per_m3: f64, // CF1PQ
    c1p_g_per_m3: f64, // CC1PQ
    m1p_g_per_m3: f64, // CM1PQ
    // coefficient 4
    al3_g_per_m3: f64, // CAL3Q
    fe3_g_per_m3: f64, // CFE3Q
    h3p_g_per_m3: f64, // CH3PQ
    f2p_g_per_m3: f64, // CF2PQ
    c2p_g_per_m3: f64, // CC2PQ
    // coefficient 5
    al4_g_per_m3: f64, // CAL4Q
    fe4_g_per_m3: f64, // CFE4Q
};

pub const Inputs = struct {
    /// ISALTG != 0 guard is caller responsibility. Kernel always computes.
    rain: IonConcentrationsRain,
    irr: IonConcentrationsIrr,
    /// PRECQ. Rain+snow rate (m3 h-1).
    rain_rate_m3_per_h: f64,
    /// PRECI. Surface irrigation rate (m3 h-1).
    surface_irr_rate_m3_per_h: f64,
    /// PRECU. Subsurface irrigation rate (m3 h-1).
    subsurface_irr_rate_m3_per_h: f64,
    /// XNFH. Timestep (h step-1).
    timestep_h: f64,
    /// Runtime phosphorus molar mass used by source factor 31 (g P mol-1).
    phosphorus_molar_mass_g_per_mol: f64,
};

pub const Increments = struct {
    /// PIR. Precipitation P flux (g P step-1). Added to TPIN.
    p_rain_g: f64,
    /// PII. Irrigation P flux (g P step-1). Added to TPIN.
    p_irr_g: f64,
    /// SIR. Precipitation ion balance increment (mol step-1). Added to TIONIN.
    ion_rain_mol: f64,
    /// SII. Irrigation ion balance increment (mol step-1). Added to TIONIN.
    ion_irr_mol: f64,
    /// SBU. Signed subsurface irrigation contribution (mol step-1). The
    /// source-negative value is added directly to TIONOU and UIONOU.
    ion_subsurface_signed_flux_mol: f64,
};

/// Valence-weighted ion sum for rain concentrations (line 4724--4735).
fn ionSumRain(c: IonConcentrationsRain) f64 {
    return
    // coefficient 1
    c.al_g_per_m3 + c.fe_g_per_m3 + c.hy_g_per_m3 + c.ca_g_per_m3 +
        c.megagrams_g_per_m3 + c.na_g_per_m3 + c.ka_g_per_m3 + c.oh_g_per_m3 +
        c.so4_g_per_m3 + c.cl_g_per_m3 + c.co3_g_per_m3 + c.h0p_g_per_m3 +
        // coefficient 2
        2.0 * (c.hco3_g_per_m3 + c.al1_g_per_m3 + c.als_g_per_m3 + c.fe1_g_per_m3 +
            c.fes_g_per_m3 + c.cao_g_per_m3 + c.cac_g_per_m3 + c.cas_g_per_m3 +
            c.mgo_g_per_m3 + c.mgc_g_per_m3 + c.mgs_g_per_m3 + c.nac_g_per_m3 +
            c.nas_g_per_m3 + c.kas_g_per_m3 + c.c0p_g_per_m3) +
        // coefficient 3
        3.0 * (c.al2_g_per_m3 + c.fe2_g_per_m3 + c.cah_g_per_m3 + c.mgh_g_per_m3 +
            c.f1p_g_per_m3 + c.c1p_g_per_m3 + c.m1p_g_per_m3) +
        // coefficient 4
        4.0 * (c.al3_g_per_m3 + c.fe3_g_per_m3 + c.h3p_g_per_m3 + c.f2p_g_per_m3 +
            c.c2p_g_per_m3) +
        // coefficient 5
        5.0 * (c.al4_g_per_m3 + c.fe4_g_per_m3);
}

/// Valence-weighted ion sum for irrigation concentrations (line 4736--4748).
fn ionSumIrr(c: IonConcentrationsIrr) f64 {
    return
    // coefficient 1
    c.al_g_per_m3 + c.fe_g_per_m3 + c.hy_g_per_m3 + c.ca_g_per_m3 +
        c.megagrams_g_per_m3 + c.na_g_per_m3 + c.ka_g_per_m3 + c.oh_g_per_m3 +
        c.so4_g_per_m3 + c.cl_g_per_m3 + c.co3_g_per_m3 + c.h0p_g_per_m3 +
        // coefficient 2
        2.0 * (c.hco3_g_per_m3 + c.al1_g_per_m3 + c.als_g_per_m3 + c.fe1_g_per_m3 +
            c.fes_g_per_m3 + c.cao_g_per_m3 + c.cac_g_per_m3 + c.cas_g_per_m3 +
            c.mgo_g_per_m3 + c.mgc_g_per_m3 + c.mgs_g_per_m3 + c.nac_g_per_m3 +
            c.nas_g_per_m3 + c.kas_g_per_m3 + c.c0p_g_per_m3) +
        // coefficient 3
        3.0 * (c.al2_g_per_m3 + c.fe2_g_per_m3 + c.cah_g_per_m3 + c.mgh_g_per_m3 +
            c.f1p_g_per_m3 + c.c1p_g_per_m3 + c.m1p_g_per_m3) +
        // coefficient 4
        4.0 * (c.al3_g_per_m3 + c.fe3_g_per_m3 + c.h3p_g_per_m3 + c.f2p_g_per_m3 +
            c.c2p_g_per_m3) +
        // coefficient 5
        5.0 * (c.al4_g_per_m3 + c.fe4_g_per_m3);
}

/// Phosphate sum for rain (used in PIR, line 4718--4720).
fn phosphateSumRain(c: IonConcentrationsRain) f64 {
    return c.h0p_g_per_m3 + c.h3p_g_per_m3 + c.f1p_g_per_m3 + c.f2p_g_per_m3 +
        c.c0p_g_per_m3 + c.c1p_g_per_m3 + c.c2p_g_per_m3 + c.m1p_g_per_m3;
}

/// Phosphate sum for irrigation (used in PII, line 4721--4723).
fn phosphateSumIrr(c: IonConcentrationsIrr) f64 {
    return c.h0p_g_per_m3 + c.h3p_g_per_m3 + c.f1p_g_per_m3 + c.f2p_g_per_m3 +
        c.c0p_g_per_m3 + c.c1p_g_per_m3 + c.c2p_g_per_m3 + c.m1p_g_per_m3;
}

/// Direct translation of redist.f lines 4717--4799 (ISALTG block body).
///
/// Caller must check ISALTG != 0 before invoking; this function always computes.
pub fn compute(inp: Inputs) !Increments {
    if (!std.math.isFinite(inp.rain_rate_m3_per_h) or inp.rain_rate_m3_per_h < 0 or
        !std.math.isFinite(inp.surface_irr_rate_m3_per_h) or inp.surface_irr_rate_m3_per_h < 0 or
        !std.math.isFinite(inp.subsurface_irr_rate_m3_per_h) or inp.subsurface_irr_rate_m3_per_h < 0 or
        !std.math.isFinite(inp.timestep_h) or inp.timestep_h <= 0 or
        !std.math.isFinite(inp.phosphorus_molar_mass_g_per_mol) or
        inp.phosphorus_molar_mass_g_per_mol <= 0)
        return error.InvalidSaltIonFluxInput;
    inline for (@typeInfo(IonConcentrationsRain).@"struct".fields) |field| {
        const concentration = @field(inp.rain, field.name);
        if (!std.math.isFinite(concentration) or concentration < 0)
            return error.InvalidSaltIonFluxInput;
    }
    inline for (@typeInfo(IonConcentrationsIrr).@"struct".fields) |field| {
        const concentration = @field(inp.irr, field.name);
        if (!std.math.isFinite(concentration) or concentration < 0)
            return error.InvalidSaltIonFluxInput;
    }

    // Lines 4718-4720: PIR
    const pir = inp.phosphorus_molar_mass_g_per_mol * inp.rain_rate_m3_per_h *
        phosphateSumRain(inp.rain) * inp.timestep_h;

    // Lines 4721-4748: PII
    const pii = inp.phosphorus_molar_mass_g_per_mol *
        inp.surface_irr_rate_m3_per_h * phosphateSumIrr(inp.irr) *
        inp.timestep_h;

    // Lines 4724-4735: SIR
    const sir = inp.rain_rate_m3_per_h * ionSumRain(inp.rain) * inp.timestep_h;

    // Lines 4736-4748: SII
    const sii = inp.surface_irr_rate_m3_per_h * ionSumIrr(inp.irr) * inp.timestep_h;

    // Lines 4783-4795: SBU
    const sbu = -inp.subsurface_irr_rate_m3_per_h * ionSumIrr(inp.irr) * inp.timestep_h;

    const result = Increments{
        .p_rain_g = pir,
        .p_irr_g = pii,
        .ion_rain_mol = sir,
        .ion_irr_mol = sii,
        .ion_subsurface_signed_flux_mol = sbu,
    };
    inline for (@typeInfo(Increments).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSaltIonFluxIncrement;
    return result;
}

fn zeroRain() IonConcentrationsRain {
    return std.mem.zeroes(IonConcentrationsRain);
}

fn zeroIrr() IonConcentrationsIrr {
    return std.mem.zeroes(IonConcentrationsIrr);
}

test "REDIST surface salt P rain uses 31 g/mol and XNFH" {
    var rain = zeroRain();
    rain.h0p_g_per_m3 = 1.0; // one PO43- species
    const result = try compute(.{
        .rain = rain,
        .irr = zeroIrr(),
        .rain_rate_m3_per_h = 0.01,
        .surface_irr_rate_m3_per_h = 0.0,
        .subsurface_irr_rate_m3_per_h = 0.0,
        .timestep_h = 2.0,
        .phosphorus_molar_mass_g_per_mol = 31.0,
    });
    // PIR = 31.0 * 0.01 * 1.0 * 2.0
    try std.testing.expectApproxEqRel(@as(f64, 31.0 * 0.01 * 1.0 * 2.0), result.p_rain_g, 1.0e-15);
    try std.testing.expectEqual(@as(f64, 0.0), result.p_irr_g);
}

test "REDIST surface salt ion rain coefficient-1 sum is direct" {
    var rain = zeroRain();
    rain.ca_g_per_m3 = 1.0;
    const result = try compute(.{
        .rain = rain,
        .irr = zeroIrr(),
        .rain_rate_m3_per_h = 0.5,
        .surface_irr_rate_m3_per_h = 0.0,
        .subsurface_irr_rate_m3_per_h = 0.0,
        .timestep_h = 1.0,
        .phosphorus_molar_mass_g_per_mol = 31.0,
    });
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.ion_rain_mol, 1.0e-15);
}

test "REDIST surface salt ion rain coefficient-2 doubles" {
    var rain = zeroRain();
    rain.hco3_g_per_m3 = 1.0; // coefficient-2 species
    const result = try compute(.{
        .rain = rain,
        .irr = zeroIrr(),
        .rain_rate_m3_per_h = 0.5,
        .surface_irr_rate_m3_per_h = 0.0,
        .subsurface_irr_rate_m3_per_h = 0.0,
        .timestep_h = 1.0,
        .phosphorus_molar_mass_g_per_mol = 31.0,
    });
    // SIR = 0.5 * 2.0 * 1.0 * 1.0
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.ion_rain_mol, 1.0e-15);
}

test "REDIST surface salt subsurface out negates irrigation sum" {
    var irr = zeroIrr();
    irr.ca_g_per_m3 = 1.0;
    const result = try compute(.{
        .rain = zeroRain(),
        .irr = irr,
        .rain_rate_m3_per_h = 0.0,
        .surface_irr_rate_m3_per_h = 0.0,
        .subsurface_irr_rate_m3_per_h = 0.1,
        .timestep_h = 0.5,
        .phosphorus_molar_mass_g_per_mol = 31.0,
    });
    // SBU = -(0.1 * 1.0 * 0.5)
    try std.testing.expectApproxEqRel(@as(f64, -0.1 * 1.0 * 0.5), result.ion_subsurface_signed_flux_mol, 1.0e-15);
}

test "REDIST surface salt accounting rejects negative molar concentration" {
    var rain = zeroRain();
    rain.al_g_per_m3 = -1.0e-6;
    try std.testing.expectError(
        error.InvalidSaltIonFluxInput,
        compute(.{
            .rain = rain,
            .irr = zeroIrr(),
            .rain_rate_m3_per_h = 0,
            .surface_irr_rate_m3_per_h = 0,
            .subsurface_irr_rate_m3_per_h = 0,
            .timestep_h = 1,
            .phosphorus_molar_mass_g_per_mol = 31,
        }),
    );
}
