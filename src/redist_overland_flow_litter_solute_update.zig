const std = @import("std");

/// Litter dissolved gas and nutrient pools at layer 0 (g).
pub const LitterSolutePools = struct {
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
    /// ZNH4S(0). NH4 (g N).
    nh4_g: f64,
    /// ZNH3S(0). NH3 (g N).
    nh3_g: f64,
    /// ZNO3S(0). NO3 (g N).
    no3_g: f64,
    /// ZNO2S(0). NO2 (g N).
    no2_g: f64,
    /// H1PO4(0). HPO4 (g P).
    hpo4_g: f64,
    /// H2PO4(0). H2PO4 (g P).
    h2po4_g: f64,
};

/// Overland flow solute fluxes into litter (g step-1).
pub const OverlandSoluteFluxes = struct {
    /// TCOQRS. CO2 overland.
    co2_g: f64,
    /// TCHQRS. CH4 overland.
    ch4_g: f64,
    /// TOXQRS. O2 overland.
    o2_g: f64,
    /// TNGQRS. N2 overland.
    n2_g: f64,
    /// TN2QRS. N2O overland.
    n2o_g: f64,
    /// THGQRS. H2 overland.
    h2_g: f64,
    /// TN4QRS. NH4 overland.
    nh4_g: f64,
    /// TN3QRS. NH3 overland.
    nh3_g: f64,
    /// TNOQRS. NO3 overland.
    no3_g: f64,
    /// TNXQRS. NO2 overland.
    no2_g: f64,
    /// TP1QRS. HPO4 overland.
    hpo4_g: f64,
    /// TPOQRS. H2PO4 overland.
    h2po4_g: f64,
};

/// Direct translation of REDIST lines 5079--5090.
///
/// Caller must check `ABS(TQR) > ZEROS` before invoking.
pub fn update(pools: LitterSolutePools, fluxes: OverlandSoluteFluxes) !LitterSolutePools {
    inline for (@typeInfo(LitterSolutePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidOverlandLitterSolutePool;
    inline for (@typeInfo(OverlandSoluteFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fluxes, field.name)))
            return error.InvalidOverlandLitterSoluteFlux;

    const result = LitterSolutePools{
        .co2_g = pools.co2_g + fluxes.co2_g,
        .ch4_g = pools.ch4_g + fluxes.ch4_g,
        .o2_g = pools.o2_g + fluxes.o2_g,
        .n2_g = pools.n2_g + fluxes.n2_g,
        .n2o_g = pools.n2o_g + fluxes.n2o_g,
        .h2_g = pools.h2_g + fluxes.h2_g,
        .nh4_g = pools.nh4_g + fluxes.nh4_g,
        .nh3_g = pools.nh3_g + fluxes.nh3_g,
        .no3_g = pools.no3_g + fluxes.no3_g,
        .no2_g = pools.no2_g + fluxes.no2_g,
        .hpo4_g = pools.hpo4_g + fluxes.hpo4_g,
        .h2po4_g = pools.h2po4_g + fluxes.h2po4_g,
    };
    inline for (@typeInfo(LitterSolutePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteOverlandLitterSolutePool;
    return result;
}

test "REDIST overland litter solute updates all 12 pools independently" {
    var fluxes = std.mem.zeroes(OverlandSoluteFluxes);
    fluxes.co2_g = 1.0;
    fluxes.nh4_g = 2.0;
    fluxes.hpo4_g = 0.5;
    const result = try update(std.mem.zeroes(LitterSolutePools), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.co2_g, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.nh4_g, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.hpo4_g, 1.0e-15);
}

test "REDIST overland litter solute preserves prior pool state" {
    var pools = std.mem.zeroes(LitterSolutePools);
    pools.no3_g = 5.0;
    const result = try update(pools, std.mem.zeroes(OverlandSoluteFluxes));
    try std.testing.expectEqual(@as(f64, 5.0), result.no3_g);
}

test "REDIST overland litter solute rejects non-finite flux" {
    var bad = std.mem.zeroes(OverlandSoluteFluxes);
    bad.ch4_g = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidOverlandLitterSoluteFlux,
        update(std.mem.zeroes(LitterSolutePools), bad),
    );
}
