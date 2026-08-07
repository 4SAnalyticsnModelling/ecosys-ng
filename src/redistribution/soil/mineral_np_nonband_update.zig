const std = @import("std");

/// Soil mineral N/P non-band micropore pools (g N or P) plus solution ions (mol).
pub const SoilNPNonbandPools = struct {
    /// ZNH3S. NH3 in non-band micropores (g N).
    nh3_g: f64,
    /// ZNH4S. NH4 in non-band micropores (g N).
    nh4_g: f64,
    /// ZNO3S. NO3 in non-band micropores (g N).
    no3_g: f64,
    /// ZNO2S. NO2 in non-band micropores (g N).
    no2_g: f64,
    /// H1PO4. HPO4 in non-band micropores (g P).
    hpo4_g: f64,
    /// H2PO4. H2PO4 in non-band micropores (g P).
    h2po4_g: f64,
};

/// Soil solution ion pools (mol): ZHY, ZOH, ZAL, ZFE, ZCA, ZMG, ZNA, ZKA.
pub const SoilSolutionIonPools = struct {
    h_mol: f64, // ZHY
    oh_mol: f64, // ZOH
    al_mol: f64, // ZAL
    fe_mol: f64, // ZFE
    ca_mol: f64, // ZCA
    megagrams_mol: f64, // ZMG
    na_mol: f64, // ZNA
    ka_mol: f64, // ZKA
};

/// N/P non-band flux terms for each of the 6 species (g N or P step-1).
pub const SoilNPNonbandFluxes = struct {
    // NH3 (line 6212)
    nh3_transport_g: f64, // TN3FLS
    nh3_dissolution_g: f64, // XN3DFG
    nh3_transformation_g: f64, // TRN3S
    nh3_root_uptake_loss_g: f64, // TUPN3S (subtracted)
    nh3_subsurface_g: f64, // RN3FLU
    nh3_pore_exchange_g: f64, // XN3FXW
    nh3_bubble_g: f64, // XN3BBL
    // NH4 (line 6215)
    nh4_transport_g: f64, // TN4FLS
    nh4_adsorption_g: f64, // XNH4S (desorption/adsorption net)
    nh4_transformation_g: f64, // TRN4S
    nh4_root_uptake_loss_g: f64, // TUPNH4 (subtracted)
    nh4_subsurface_g: f64, // RN4FLU
    nh4_pore_exchange_g: f64, // XN4FXW
    // NO3 (line 6218)
    no3_transport_g: f64, // TNOFLS
    no3_nitrification_g: f64, // XNO3S
    no3_transformation_g: f64, // TRNO3
    no3_root_uptake_loss_g: f64, // TUPNO3 (subtracted)
    no3_subsurface_g: f64, // RNOFLU
    no3_pore_exchange_g: f64, // XNOFXW
    // NO2 (line 6221)
    no2_transport_g: f64, // TNXFLS
    no2_nitrification_g: f64, // XNO2S
    no2_transformation_g: f64, // TRNO2
    no2_pore_exchange_g: f64, // XNXFXS
    // HPO4 (line 6223)
    hpo4_transport_g: f64, // TP1FLS
    hpo4_desorption_g: f64, // XH1PS
    hpo4_transformation_g: f64, // TRH1P
    hpo4_root_uptake_loss_g: f64, // TUPH1P (subtracted)
    hpo4_subsurface_g: f64, // RH1PFU
    hpo4_pore_exchange_g: f64, // XH1PXS
    // H2PO4 (line 6225)
    h2po4_transport_g: f64, // TPOFLS
    h2po4_desorption_g: f64, // XH2PS
    h2po4_transformation_g: f64, // TRH2P
    h2po4_root_uptake_loss_g: f64, // TUPH2P (subtracted)
    h2po4_subsurface_g: f64, // RH2PFU
    h2po4_pore_exchange_g: f64, // XH2PXS
};

/// Transformation increments for solution ions (mol step-1).
pub const SoilSolutionIonFluxes = struct {
    h_mol: f64, // TRHY
    oh_mol: f64, // TROH
    al_mol: f64, // TRAL
    fe_mol: f64, // TRFE
    ca_mol: f64, // TRCA
    megagrams_mol: f64, // TRMG
    na_mol: f64, // TRNA
    ka_mol: f64, // TRKA
};

pub const Result = struct {
    np: SoilNPNonbandPools,
    ions: SoilSolutionIonPools,
};

/// Direct translation of redist.f lines 6212--6267 (inner body of DO 125 L loop).
pub fn update(
    np: SoilNPNonbandPools,
    ions: SoilSolutionIonPools,
    nf: SoilNPNonbandFluxes,
    ionf: SoilSolutionIonFluxes,
) !Result {
    inline for (@typeInfo(SoilNPNonbandPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(np, field.name)))
            return error.InvalidSoilNPPool;
    inline for (@typeInfo(SoilSolutionIonPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ions, field.name)))
            return error.InvalidSoilIonPool;
    inline for (@typeInfo(SoilNPNonbandFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(nf, field.name)))
            return error.InvalidSoilNPFlux;
    inline for (@typeInfo(SoilSolutionIonFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ionf, field.name)))
            return error.InvalidSoilIonFlux;

    const new_np = SoilNPNonbandPools{
        .nh3_g = np.nh3_g + nf.nh3_transport_g + nf.nh3_dissolution_g + nf.nh3_transformation_g - nf.nh3_root_uptake_loss_g + nf.nh3_subsurface_g + nf.nh3_pore_exchange_g + nf.nh3_bubble_g,
        .nh4_g = np.nh4_g + nf.nh4_transport_g + nf.nh4_adsorption_g + nf.nh4_transformation_g - nf.nh4_root_uptake_loss_g + nf.nh4_subsurface_g + nf.nh4_pore_exchange_g,
        .no3_g = np.no3_g + nf.no3_transport_g + nf.no3_nitrification_g + nf.no3_transformation_g - nf.no3_root_uptake_loss_g + nf.no3_subsurface_g + nf.no3_pore_exchange_g,
        .no2_g = np.no2_g + nf.no2_transport_g + nf.no2_nitrification_g + nf.no2_transformation_g + nf.no2_pore_exchange_g,
        .hpo4_g = np.hpo4_g + nf.hpo4_transport_g + nf.hpo4_desorption_g + nf.hpo4_transformation_g - nf.hpo4_root_uptake_loss_g + nf.hpo4_subsurface_g + nf.hpo4_pore_exchange_g,
        .h2po4_g = np.h2po4_g + nf.h2po4_transport_g + nf.h2po4_desorption_g + nf.h2po4_transformation_g - nf.h2po4_root_uptake_loss_g + nf.h2po4_subsurface_g + nf.h2po4_pore_exchange_g,
    };
    const new_ions = SoilSolutionIonPools{
        .h_mol = ions.h_mol + ionf.h_mol,
        .oh_mol = ions.oh_mol + ionf.oh_mol,
        .al_mol = ions.al_mol + ionf.al_mol,
        .fe_mol = ions.fe_mol + ionf.fe_mol,
        .ca_mol = ions.ca_mol + ionf.ca_mol,
        .megagrams_mol = ions.megagrams_mol + ionf.megagrams_mol,
        .na_mol = ions.na_mol + ionf.na_mol,
        .ka_mol = ions.ka_mol + ionf.ka_mol,
    };

    inline for (@typeInfo(SoilNPNonbandPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(new_np, field.name)))
            return error.NonFiniteSoilNPPool;
    inline for (@typeInfo(SoilSolutionIonPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(new_ions, field.name)))
            return error.NonFiniteSoilIonPool;
    return Result{ .np = new_np, .ions = new_ions };
}

/// Applies the REDIST DO 125 body to the runtime-configured soil layer span.
pub fn updateLayers(
    np_by_layer: []SoilNPNonbandPools,
    ions_by_layer: []SoilSolutionIonPools,
    np_fluxes_by_layer: []const SoilNPNonbandFluxes,
    ion_fluxes_by_layer: []const SoilSolutionIonFluxes,
) !void {
    if (np_by_layer.len == 0 or
        ions_by_layer.len != np_by_layer.len or
        np_fluxes_by_layer.len != np_by_layer.len or
        ion_fluxes_by_layer.len != np_by_layer.len)
        return error.SoilNPDimensionMismatch;

    for (
        np_by_layer,
        ions_by_layer,
        np_fluxes_by_layer,
        ion_fluxes_by_layer,
    ) |*np, *ions, np_fluxes, ion_fluxes| {
        const result = try update(np.*, ions.*, np_fluxes, ion_fluxes);
        np.* = result.np;
        ions.* = result.ions;
    }
}

test "REDIST soil mineral N/P NH4 transport and adsorption add to pool" {
    var nf = std.mem.zeroes(SoilNPNonbandFluxes);
    nf.nh4_transport_g = 2.0;
    nf.nh4_adsorption_g = 1.0;
    const result = try update(
        std.mem.zeroes(SoilNPNonbandPools),
        std.mem.zeroes(SoilSolutionIonPools),
        nf,
        std.mem.zeroes(SoilSolutionIonFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.np.nh4_g, 1.0e-15);
}

test "REDIST soil mineral N/P root uptake subtracts from NH3 pool" {
    var nf = std.mem.zeroes(SoilNPNonbandFluxes);
    nf.nh3_dissolution_g = 5.0;
    nf.nh3_root_uptake_loss_g = 2.0;
    const result = try update(
        std.mem.zeroes(SoilNPNonbandPools),
        std.mem.zeroes(SoilSolutionIonPools),
        nf,
        std.mem.zeroes(SoilSolutionIonFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.np.nh3_g, 1.0e-15);
}

test "REDIST soil solution ions accumulate from transformations" {
    var ionf = std.mem.zeroes(SoilSolutionIonFluxes);
    ionf.ca_mol = 1.5;
    ionf.megagrams_mol = 0.5;
    const result = try update(
        std.mem.zeroes(SoilNPNonbandPools),
        std.mem.zeroes(SoilSolutionIonPools),
        std.mem.zeroes(SoilNPNonbandFluxes),
        ionf,
    );
    try std.testing.expectApproxEqRel(@as(f64, 1.5), result.ions.ca_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.ions.megagrams_mol, 1.0e-15);
}

test "REDIST soil mineral N/P rejects non-finite N flux" {
    var bad = std.mem.zeroes(SoilNPNonbandFluxes);
    bad.no3_transformation_g = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSoilNPFlux,
        update(
            std.mem.zeroes(SoilNPNonbandPools),
            std.mem.zeroes(SoilSolutionIonPools),
            bad,
            std.mem.zeroes(SoilSolutionIonFluxes),
        ),
    );
}

test "REDIST soil mineral N/P preserves source term order and signs" {
    var pools = std.mem.zeroes(SoilNPNonbandPools);
    pools.nh3_g = 10.0;
    var fluxes = std.mem.zeroes(SoilNPNonbandFluxes);
    fluxes.nh3_transport_g = 1.0;
    fluxes.nh3_dissolution_g = 2.0;
    fluxes.nh3_transformation_g = 3.0;
    fluxes.nh3_root_uptake_loss_g = 4.0;
    fluxes.nh3_subsurface_g = 5.0;
    fluxes.nh3_pore_exchange_g = 6.0;
    fluxes.nh3_bubble_g = 7.0;

    const result = try update(
        pools,
        std.mem.zeroes(SoilSolutionIonPools),
        fluxes,
        std.mem.zeroes(SoilSolutionIonFluxes),
    );
    try std.testing.expectEqual(@as(f64, 30.0), result.np.nh3_g);
}

test "REDIST soil solution ion transformation maps every species" {
    const fluxes = SoilSolutionIonFluxes{
        .h_mol = 1.0,
        .oh_mol = 2.0,
        .al_mol = 3.0,
        .fe_mol = 4.0,
        .ca_mol = 5.0,
        .megagrams_mol = 6.0,
        .na_mol = 7.0,
        .ka_mol = 8.0,
    };
    const result = try update(
        std.mem.zeroes(SoilNPNonbandPools),
        std.mem.zeroes(SoilSolutionIonPools),
        std.mem.zeroes(SoilNPNonbandFluxes),
        fluxes,
    );
    inline for (@typeInfo(SoilSolutionIonPools).@"struct".fields, 1..) |field, expected|
        try std.testing.expectEqual(@as(f64, @floatFromInt(expected)), @field(result.ions, field.name));
}

test "REDIST soil mineral N/P runtime layer traversal is independent" {
    var pools = [_]SoilNPNonbandPools{
        std.mem.zeroes(SoilNPNonbandPools),
        std.mem.zeroes(SoilNPNonbandPools),
    };
    var ions = [_]SoilSolutionIonPools{
        std.mem.zeroes(SoilSolutionIonPools),
        std.mem.zeroes(SoilSolutionIonPools),
    };
    var np_fluxes = [_]SoilNPNonbandFluxes{
        std.mem.zeroes(SoilNPNonbandFluxes),
        std.mem.zeroes(SoilNPNonbandFluxes),
    };
    var ion_fluxes = [_]SoilSolutionIonFluxes{
        std.mem.zeroes(SoilSolutionIonFluxes),
        std.mem.zeroes(SoilSolutionIonFluxes),
    };
    np_fluxes[0].no2_transport_g = 2.5;
    np_fluxes[1].h2po4_desorption_g = 4.5;
    ion_fluxes[0].al_mol = 6.5;
    ion_fluxes[1].na_mol = 8.5;

    try updateLayers(&pools, &ions, &np_fluxes, &ion_fluxes);

    try std.testing.expectEqual(@as(f64, 2.5), pools[0].no2_g);
    try std.testing.expectEqual(@as(f64, 0.0), pools[1].no2_g);
    try std.testing.expectEqual(@as(f64, 0.0), pools[0].h2po4_g);
    try std.testing.expectEqual(@as(f64, 4.5), pools[1].h2po4_g);
    try std.testing.expectEqual(@as(f64, 6.5), ions[0].al_mol);
    try std.testing.expectEqual(@as(f64, 8.5), ions[1].na_mol);
}

test "REDIST soil mineral N/P runtime traversal rejects mismatched dimensions" {
    var pools = [_]SoilNPNonbandPools{std.mem.zeroes(SoilNPNonbandPools)};
    var ions: [0]SoilSolutionIonPools = .{};
    const np_fluxes = [_]SoilNPNonbandFluxes{std.mem.zeroes(SoilNPNonbandFluxes)};
    const ion_fluxes = [_]SoilSolutionIonFluxes{std.mem.zeroes(SoilSolutionIonFluxes)};
    try std.testing.expectError(
        error.SoilNPDimensionMismatch,
        updateLayers(&pools, &ions, &np_fluxes, &ion_fluxes),
    );
}

test "REDIST soil mineral N/P rejects invalid and overflowing ion pools" {
    var invalid = std.mem.zeroes(SoilSolutionIonPools);
    invalid.h_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSoilIonPool,
        update(
            std.mem.zeroes(SoilNPNonbandPools),
            invalid,
            std.mem.zeroes(SoilNPNonbandFluxes),
            std.mem.zeroes(SoilSolutionIonFluxes),
        ),
    );

    var finite = std.mem.zeroes(SoilSolutionIonPools);
    finite.ca_mol = std.math.floatMax(f64);
    var overflowing = std.mem.zeroes(SoilSolutionIonFluxes);
    overflowing.ca_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSoilIonPool,
        update(
            std.mem.zeroes(SoilNPNonbandPools),
            finite,
            std.mem.zeroes(SoilNPNonbandFluxes),
            overflowing,
        ),
    );
}
