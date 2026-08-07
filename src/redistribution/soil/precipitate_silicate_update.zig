const std = @import("std");

/// Soil hydroxide precipitate pools (mol). 4 fields, REDIST 7175-7178.
pub const HydroxidePrecipitatePools = struct {
    aloh_mol: f64, // PALOH
    feoh_mol: f64, // PFEOH
    caco_mol: f64, // PCACO
    caso_mol: f64, // PCASO
};

/// Silicate mineral pools (mol). 12 fields, REDIST 7179-7190.
/// 6 primary silicates (QALSI/QFESI/QCASI/QMGSI/QNASI/QKASI)
/// + 6 secondary silicates (QALSIF/QFESIF/QCASIF/QMGSIF/QNASIF/QKASIF).
pub const SilicatePools = struct {
    al_mol: f64, // QALSI
    fe_mol: f64, // QFESI
    ca_mol: f64, // QCASI
    megagrams_mol: f64, // QMGSI
    na_mol: f64, // QNASI
    ka_mol: f64, // QKASI
    al2_mol: f64, // QALSIF (secondary)
    fe2_mol: f64, // QFESIF
    ca2_mol: f64, // QCASIF
    mg2_mol: f64, // QMGSIF
    na2_mol: f64, // QNASIF
    ka2_mol: f64, // QKASIF
};

/// Transformation fluxes for hydroxide precipitates (mol step-1).
pub const HydroxidePrecipitateFluxes = struct {
    aloh_mol: f64, // TRALOH
    feoh_mol: f64, // TRFEOH
    caco_mol: f64, // TRCACO
    caso_mol: f64, // TRCASO
};

/// Weathering/dissolution transformation fluxes for silicates (mol step-1).
pub const SilicateFluxes = struct {
    al_mol: f64, // TRALSI
    fe_mol: f64, // TRFESI
    ca_mol: f64, // TRCASI
    megagrams_mol: f64, // TRMGSI
    na_mol: f64, // TRNASI
    ka_mol: f64, // TRKASI
    al2_mol: f64, // TRALSIF
    fe2_mol: f64, // TRFESIF
    ca2_mol: f64, // TRCASIF
    mg2_mol: f64, // TRMGSIF
    na2_mol: f64, // TRNASIF
    ka2_mol: f64, // TRKASIF
};

pub const Result = struct {
    hydroxide: HydroxidePrecipitatePools,
    silicate: SilicatePools,
    /// UCO2S increment = 12*PCACO after update (line 7198).
    uco2s_increment_gC: f64,
};

pub const SilicateTotals = struct {
    al_mol: f64, // TQALSI
    fe_mol: f64, // TQFESI
    ca_mol: f64, // TQCASI
    megagrams_mol: f64, // TQMGSI
    na_mol: f64, // TQNASI
    ka_mol: f64, // TQKASI
    charge_equivalent_mol: f64, // TQSI
};

pub const SaltEquilibriumMode = enum { static, dynamic };

/// Direct translation of redist.f lines 7175--7198.
/// Gated on ISALTG != 0 by the caller.
pub fn update(
    hydroxide: HydroxidePrecipitatePools,
    silicate: SilicatePools,
    hf: HydroxidePrecipitateFluxes,
    sf: SilicateFluxes,
) !Result {
    inline for (@typeInfo(HydroxidePrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(hydroxide, field.name)))
            return error.InvalidHydroxidePrecipPool;
    inline for (@typeInfo(SilicatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(silicate, field.name)))
            return error.InvalidSilicatePool;
    inline for (@typeInfo(HydroxidePrecipitateFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(hf, field.name)))
            return error.InvalidHydroxidePrecipFlux;
    inline for (@typeInfo(SilicateFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(sf, field.name)))
            return error.InvalidSilicateFlux;

    const new_hydroxide = HydroxidePrecipitatePools{
        .aloh_mol = hydroxide.aloh_mol + hf.aloh_mol,
        .feoh_mol = hydroxide.feoh_mol + hf.feoh_mol,
        .caco_mol = hydroxide.caco_mol + hf.caco_mol,
        .caso_mol = hydroxide.caso_mol + hf.caso_mol,
    };
    const new_silicate = SilicatePools{
        .al_mol = silicate.al_mol + sf.al_mol,
        .fe_mol = silicate.fe_mol + sf.fe_mol,
        .ca_mol = silicate.ca_mol + sf.ca_mol,
        .megagrams_mol = silicate.megagrams_mol + sf.megagrams_mol,
        .na_mol = silicate.na_mol + sf.na_mol,
        .ka_mol = silicate.ka_mol + sf.ka_mol,
        .al2_mol = silicate.al2_mol + sf.al2_mol,
        .fe2_mol = silicate.fe2_mol + sf.fe2_mol,
        .ca2_mol = silicate.ca2_mol + sf.ca2_mol,
        .mg2_mol = silicate.mg2_mol + sf.mg2_mol,
        .na2_mol = silicate.na2_mol + sf.na2_mol,
        .ka2_mol = silicate.ka2_mol + sf.ka2_mol,
    };
    inline for (@typeInfo(HydroxidePrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(new_hydroxide, field.name)))
            return error.NonFiniteHydroxidePrecipPool;
    inline for (@typeInfo(SilicatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(new_silicate, field.name)))
            return error.NonFiniteSilicatePool;

    const uco2s_increment_gC = 12.0 * new_hydroxide.caco_mol;
    if (!std.math.isFinite(uco2s_increment_gC))
        return error.NonFinitePrecipitateCarbonDiagnostic;
    return Result{
        .hydroxide = new_hydroxide,
        .silicate = new_silicate,
        .uco2s_increment_gC = uco2s_increment_gC,
    };
}

/// Applies REDIST 7175--7198 to ascending runtime layers. Silicate totals are
/// running across layers; TQSI is reassigned after each layer as in the source.
pub fn updateLayers(
    hydroxide_by_layer: []HydroxidePrecipitatePools,
    silicate_by_layer: []SilicatePools,
    hydroxide_fluxes_by_layer: []const HydroxidePrecipitateFluxes,
    silicate_fluxes_by_layer: []const SilicateFluxes,
    totals: *SilicateTotals,
    grid_cell_carbon_g_c: *f64,
    mode: SaltEquilibriumMode,
) !void {
    if (hydroxide_by_layer.len == 0 or
        silicate_by_layer.len != hydroxide_by_layer.len or
        hydroxide_fluxes_by_layer.len != hydroxide_by_layer.len or
        silicate_fluxes_by_layer.len != hydroxide_by_layer.len)
        return error.PrecipitateSilicateDimensionMismatch;
    inline for (@typeInfo(SilicateTotals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals, field.name)))
            return error.InvalidSilicateTotal;
    if (!std.math.isFinite(grid_cell_carbon_g_c.*))
        return error.InvalidPrecipitateCarbonDiagnostic;
    if (mode == .static) return;

    var next_totals = totals.*;
    var next_carbon_g_c = grid_cell_carbon_g_c.*;
    for (
        hydroxide_by_layer,
        silicate_by_layer,
        hydroxide_fluxes_by_layer,
        silicate_fluxes_by_layer,
    ) |*hydroxide, *silicate, hydroxide_fluxes, silicate_fluxes| {
        const result = try update(hydroxide.*, silicate.*, hydroxide_fluxes, silicate_fluxes);
        next_totals.al_mol = next_totals.al_mol + result.silicate.al_mol + result.silicate.al2_mol;
        next_totals.fe_mol = next_totals.fe_mol + result.silicate.fe_mol + result.silicate.fe2_mol;
        next_totals.ca_mol = next_totals.ca_mol + result.silicate.ca_mol + result.silicate.ca2_mol;
        next_totals.megagrams_mol = next_totals.megagrams_mol + result.silicate.megagrams_mol + result.silicate.mg2_mol;
        next_totals.na_mol = next_totals.na_mol + result.silicate.na_mol + result.silicate.na2_mol;
        next_totals.ka_mol = next_totals.ka_mol + result.silicate.ka_mol + result.silicate.ka2_mol;
        next_totals.charge_equivalent_mol = 3.0 * (next_totals.al_mol + next_totals.fe_mol) +
            2.0 * (next_totals.ca_mol + next_totals.megagrams_mol) + next_totals.na_mol + next_totals.ka_mol;
        next_carbon_g_c = next_carbon_g_c + result.uco2s_increment_gC;
        inline for (@typeInfo(SilicateTotals).@"struct".fields) |field|
            if (!std.math.isFinite(@field(next_totals, field.name)))
                return error.NonFiniteSilicateTotal;
        if (!std.math.isFinite(next_carbon_g_c))
            return error.NonFinitePrecipitateCarbonDiagnostic;
        hydroxide.* = result.hydroxide;
        silicate.* = result.silicate;
    }
    totals.* = next_totals;
    grid_cell_carbon_g_c.* = next_carbon_g_c;
}

test "REDIST soil precipitate-silicate CaCO3 accumulates from dissolution" {
    var hf = std.mem.zeroes(HydroxidePrecipitateFluxes);
    hf.caco_mol = 2.0;
    const result = try update(
        std.mem.zeroes(HydroxidePrecipitatePools),
        std.mem.zeroes(SilicatePools),
        hf,
        std.mem.zeroes(SilicateFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.hydroxide.caco_mol, 1.0e-15);
}

test "REDIST soil precipitate-silicate UCO2S is 12*PCACO" {
    var hf = std.mem.zeroes(HydroxidePrecipitateFluxes);
    hf.caco_mol = 3.0;
    const result = try update(
        std.mem.zeroes(HydroxidePrecipitatePools),
        std.mem.zeroes(SilicatePools),
        hf,
        std.mem.zeroes(SilicateFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 36.0), result.uco2s_increment_gC, 1.0e-15);
}

test "REDIST soil silicate primary and secondary Al weathering are independent" {
    var sf = std.mem.zeroes(SilicateFluxes);
    sf.al_mol = 1.0;
    sf.al2_mol = 2.0;
    const result = try update(
        std.mem.zeroes(HydroxidePrecipitatePools),
        std.mem.zeroes(SilicatePools),
        std.mem.zeroes(HydroxidePrecipitateFluxes),
        sf,
    );
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.silicate.al_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.silicate.al2_mol, 1.0e-15);
}

test "REDIST soil precipitate-silicate rejects non-finite silicate flux" {
    var bad = std.mem.zeroes(SilicateFluxes);
    bad.mg2_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSilicateFlux,
        update(
            std.mem.zeroes(HydroxidePrecipitatePools),
            std.mem.zeroes(SilicatePools),
            std.mem.zeroes(HydroxidePrecipitateFluxes),
            bad,
        ),
    );
}

test "REDIST precipitate silicate runtime layers preserve running total order" {
    var hydroxide = [_]HydroxidePrecipitatePools{
        std.mem.zeroes(HydroxidePrecipitatePools),
        std.mem.zeroes(HydroxidePrecipitatePools),
    };
    var silicate = [_]SilicatePools{ std.mem.zeroes(SilicatePools), std.mem.zeroes(SilicatePools) };
    var hydroxide_fluxes = [_]HydroxidePrecipitateFluxes{
        std.mem.zeroes(HydroxidePrecipitateFluxes),
        std.mem.zeroes(HydroxidePrecipitateFluxes),
    };
    var silicate_fluxes = [_]SilicateFluxes{ std.mem.zeroes(SilicateFluxes), std.mem.zeroes(SilicateFluxes) };
    hydroxide_fluxes[0].caco_mol = 2;
    hydroxide_fluxes[1].caco_mol = 3;
    silicate_fluxes[0].al_mol = 1;
    silicate_fluxes[0].al2_mol = 2;
    silicate_fluxes[1].ca_mol = 4;
    silicate_fluxes[1].ca2_mol = 5;
    var totals = std.mem.zeroes(SilicateTotals);
    var carbon_g_c: f64 = 0;
    try updateLayers(&hydroxide, &silicate, &hydroxide_fluxes, &silicate_fluxes, &totals, &carbon_g_c, .dynamic);
    try std.testing.expectEqual(@as(f64, 3), totals.al_mol);
    try std.testing.expectEqual(@as(f64, 9), totals.ca_mol);
    try std.testing.expectEqual(@as(f64, 27), totals.charge_equivalent_mol);
    try std.testing.expectEqual(@as(f64, 60), carbon_g_c);
}

test "REDIST precipitate silicate static mode and failures" {
    var hydroxide = [_]HydroxidePrecipitatePools{std.mem.zeroes(HydroxidePrecipitatePools)};
    var silicate = [_]SilicatePools{std.mem.zeroes(SilicatePools)};
    const hydroxide_fluxes = [_]HydroxidePrecipitateFluxes{std.mem.zeroes(HydroxidePrecipitateFluxes)};
    const silicate_fluxes = [_]SilicateFluxes{std.mem.zeroes(SilicateFluxes)};
    var totals = std.mem.zeroes(SilicateTotals);
    var carbon_g_c: f64 = 2;
    try updateLayers(&hydroxide, &silicate, &hydroxide_fluxes, &silicate_fluxes, &totals, &carbon_g_c, .static);
    try std.testing.expectEqual(@as(f64, 2), carbon_g_c);

    const no_fluxes: [0]SilicateFluxes = .{};
    try std.testing.expectError(
        error.PrecipitateSilicateDimensionMismatch,
        updateLayers(&hydroxide, &silicate, &hydroxide_fluxes, &no_fluxes, &totals, &carbon_g_c, .dynamic),
    );
    totals.al_mol = std.math.floatMax(f64);
    var overflowing = silicate_fluxes;
    overflowing[0].al_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSilicateTotal,
        updateLayers(&hydroxide, &silicate, &hydroxide_fluxes, &overflowing, &totals, &carbon_g_c, .dynamic),
    );
}
