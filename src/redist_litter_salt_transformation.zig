const std = @import("std");

/// Eight base-ion litter pools (mol) at layer 0.
pub const LitterIonPools = struct {
    /// ZHY(0). H+ content (mol).
    h_mol: f64,
    /// ZOH(0). OH- content (mol).
    oh_mol: f64,
    /// ZAL(0). Al3+ content (mol).
    al_mol: f64,
    /// ZFE(0). Fe3+ content (mol).
    fe_mol: f64,
    /// ZCA(0). Ca2+ content (mol).
    ca_mol: f64,
    /// ZMG(0). Mg2+ content (mol).
    megagrams_mol: f64,
    /// ZNA(0). Na+ content (mol).
    na_mol: f64,
    /// ZKA(0). K+ content (mol).
    ka_mol: f64,
};

/// Signed, step-integrated `solute.f` transformation increments for each base
/// ion (mol per model step). Positive values add to the litter solution pool.
pub const IonTransformations = struct {
    /// TRHY(0).
    h_mol: f64,
    /// TROH(0).
    oh_mol: f64,
    /// TRAL(0).
    al_mol: f64,
    /// TRFE(0).
    fe_mol: f64,
    /// TRCA(0).
    ca_mol: f64,
    /// TRMG(0).
    megagrams_mol: f64,
    /// TRNA(0).
    na_mol: f64,
    /// TRKA(0).
    ka_mol: f64,
};

/// Direct translation of REDIST lines 4885--4892.
///
/// Applies solute.f net transformation increments to litter-layer base ions.
pub fn apply(pools: LitterIonPools, tr: IonTransformations) !LitterIonPools {
    inline for (@typeInfo(LitterIonPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidLitterIonPool;
    inline for (@typeInfo(IonTransformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(tr, field.name)))
            return error.InvalidLitterIonTransformation;

    const result = LitterIonPools{
        .h_mol = pools.h_mol + tr.h_mol,
        .oh_mol = pools.oh_mol + tr.oh_mol,
        .al_mol = pools.al_mol + tr.al_mol,
        .fe_mol = pools.fe_mol + tr.fe_mol,
        .ca_mol = pools.ca_mol + tr.ca_mol,
        .megagrams_mol = pools.megagrams_mol + tr.megagrams_mol,
        .na_mol = pools.na_mol + tr.na_mol,
        .ka_mol = pools.ka_mol + tr.ka_mol,
    };
    inline for (@typeInfo(LitterIonPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteLitterIonPool;
    return result;
}

test "REDIST litter salt transformation accumulates all eight ions" {
    const pools = LitterIonPools{ .h_mol = 1.0, .oh_mol = 2.0, .al_mol = 3.0, .fe_mol = 4.0, .ca_mol = 5.0, .megagrams_mol = 6.0, .na_mol = 7.0, .ka_mol = 8.0 };
    const tr = IonTransformations{ .h_mol = 0.1, .oh_mol = 0.2, .al_mol = 0.3, .fe_mol = 0.4, .ca_mol = 0.5, .megagrams_mol = 0.6, .na_mol = 0.7, .ka_mol = 0.8 };
    const result = try apply(pools, tr);
    try std.testing.expectApproxEqRel(@as(f64, 1.1), result.h_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.2), result.oh_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 3.3), result.al_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 4.4), result.fe_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 5.5), result.ca_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 6.6), result.megagrams_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 7.7), result.na_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 8.8), result.ka_mol, 1.0e-15);
}

test "REDIST litter salt transformation zero increments preserves state" {
    const pools = LitterIonPools{ .h_mol = 5.0, .oh_mol = 0.0, .al_mol = 0.0, .fe_mol = 0.0, .ca_mol = 0.0, .megagrams_mol = 0.0, .na_mol = 0.0, .ka_mol = 0.0 };
    const result = try apply(pools, std.mem.zeroes(IonTransformations));
    try std.testing.expectEqual(@as(f64, 5.0), result.h_mol);
}

test "REDIST litter salt transformation rejects non-finite input" {
    try std.testing.expectError(
        error.InvalidLitterIonPool,
        apply(
            .{ .h_mol = std.math.nan(f64), .oh_mol = 0.0, .al_mol = 0.0, .fe_mol = 0.0, .ca_mol = 0.0, .megagrams_mol = 0.0, .na_mol = 0.0, .ka_mol = 0.0 },
            std.mem.zeroes(IonTransformations),
        ),
    );
}

test "REDIST litter salt transformation preserves signed removal" {
    const pools = LitterIonPools{ .h_mol = 1.0, .oh_mol = 2.0, .al_mol = 3.0, .fe_mol = 4.0, .ca_mol = 5.0, .megagrams_mol = 6.0, .na_mol = 7.0, .ka_mol = 8.0 };
    const tr = IonTransformations{ .h_mol = -0.1, .oh_mol = -0.2, .al_mol = -0.3, .fe_mol = -0.4, .ca_mol = -0.5, .megagrams_mol = -0.6, .na_mol = -0.7, .ka_mol = -0.8 };
    const result = try apply(pools, tr);
    try std.testing.expectApproxEqRel(@as(f64, 0.9), result.h_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 7.2), result.ka_mol, 1.0e-15);
}

test "REDIST litter salt transformation rejects non-finite increment" {
    var tr = std.mem.zeroes(IonTransformations);
    tr.fe_mol = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidLitterIonTransformation,
        apply(std.mem.zeroes(LitterIonPools), tr),
    );
}

test "REDIST litter salt transformation rejects arithmetic overflow" {
    var pools = std.mem.zeroes(LitterIonPools);
    var tr = std.mem.zeroes(IonTransformations);
    pools.ca_mol = std.math.floatMax(f64);
    tr.ca_mol = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteLitterIonPool, apply(pools, tr));
}
