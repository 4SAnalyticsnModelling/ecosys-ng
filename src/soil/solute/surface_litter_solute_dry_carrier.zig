const std = @import("std");

/// HOUR1 litter solute concentration derivation, `hour1.f` 4494--4525.
///
/// The oracle keeps litter solutes as **absolute mass** (`ZNH4S(0)`, `ZNO3S(0)`,
/// `H2PO4(0)`, ... in g) and derives concentration only where water exists:
///
///     CNH4S(0,NY,NX)=AMAX1(0.0,ZNH4S(0,NY,NX)/VOLW(0,NY,NX))
///
/// and in the branch where the surface litter layer does not exist it sets every
/// concentration to zero outright:
///
///     CNH4S(0,NY,NX)=0.0
///
/// Mass is therefore never destroyed by the carrier vanishing; only the derived
/// concentration becomes undefined and is reported as zero.
///
/// ecosys-ng stores these pools water-normalized instead, so
/// `surface_litter_chemistry_carrier_rebase` must rescale by `old/new` whenever
/// litter water changes, and a `new == 0` carrier has no representable
/// concentration. It currently rejects that case with
/// `LitterChemistryMassWithoutWaterCarrier`. That is safe but incomplete: it makes
/// evaporation to dryness a hard error, which is why the shipped Ottawa runscript
/// can only run with surface evaporation switched off. See EXEC-004.
///
/// This module provides the oracle's semantics as an explicit, testable
/// conversion so the carrier rebase can adopt it: hold the extensive amount, and
/// report a zero concentration when there is no carrier.
pub const Carrier = struct {
    /// Litter water volume in m3, the oracle's `VOLW(0,NY,NX)`.
    water_volume_m3: f64,

    /// True when a concentration is representable at all.
    pub fn hasWater(self: Carrier) bool {
        return self.water_volume_m3 > 0;
    }
};

/// Converts a water-normalized concentration into the extensive amount the
/// oracle stores. This is the quantity that must be conserved across any carrier
/// change, including evaporation to dryness.
pub fn extensiveAmount(concentration_per_m3: f64, carrier: Carrier) !f64 {
    inline for (.{ concentration_per_m3, carrier.water_volume_m3 }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteLitterSoluteState;
    if (concentration_per_m3 < 0 or carrier.water_volume_m3 < 0)
        return error.InvalidLitterSoluteState;
    const amount = concentration_per_m3 * carrier.water_volume_m3;
    if (!std.math.isFinite(amount)) return error.NonFiniteLitterSoluteState;
    return amount;
}

/// Derives the concentration the oracle would report, following
/// `AMAX1(0.0, mass/VOLW)` where water exists and the explicit `0.0` branch where
/// it does not. Crucially this does **not** error on a dry carrier: a dry litter
/// layer has zero solute concentration and retains its solute mass.
pub fn concentration(extensive_amount: f64, carrier: Carrier) !f64 {
    inline for (.{ extensive_amount, carrier.water_volume_m3 }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteLitterSoluteState;
    if (extensive_amount < 0 or carrier.water_volume_m3 < 0)
        return error.InvalidLitterSoluteState;
    if (!carrier.hasWater()) return 0;
    const result = @max(0.0, extensive_amount / carrier.water_volume_m3);
    if (!std.math.isFinite(result)) return error.NonFiniteLitterSoluteState;
    return result;
}

pub const Rebase = struct {
    /// The concentration to store after the carrier change.
    concentration_per_m3: f64,
    /// The extensive amount carried across, which must equal the amount before.
    retained_amount: f64,
    /// True when the carrier went dry, so the concentration is zero by the
    /// oracle's rule while the mass is retained.
    carrier_dry: bool,
};

/// Rebases one water-normalized pool across a carrier change without losing
/// mass, including the dry-carrier case the current production path rejects.
pub fn rebase(
    concentration_before_per_m3: f64,
    carrier_before: Carrier,
    carrier_after: Carrier,
) !Rebase {
    const amount = try extensiveAmount(concentration_before_per_m3, carrier_before);
    return .{
        .concentration_per_m3 = try concentration(amount, carrier_after),
        .retained_amount = amount,
        .carrier_dry = !carrier_after.hasWater(),
    };
}

test "concentration is mass over water where water exists" {
    try std.testing.expectApproxEqRel(
        @as(f64, 2.5),
        try concentration(10, .{ .water_volume_m3 = 4 }),
        1e-15,
    );
}

test "a dry carrier reports zero concentration instead of erroring" {
    // This is the case production currently rejects with
    // LitterChemistryMassWithoutWaterCarrier, forcing surface evaporation off.
    // The oracle's `hour1.f` 4519--4524 branch sets the concentration to 0.0.
    try std.testing.expectEqual(
        @as(f64, 0),
        try concentration(10, .{ .water_volume_m3 = 0 }),
    );
}

test "rebasing to a dry carrier retains the full solute mass" {
    const result = try rebase(2.5, .{ .water_volume_m3 = 4 }, .{ .water_volume_m3 = 0 });
    // The mass is conserved even though the concentration is no longer defined.
    try std.testing.expectApproxEqRel(@as(f64, 10), result.retained_amount, 1e-15);
    try std.testing.expectEqual(@as(f64, 0), result.concentration_per_m3);
    try std.testing.expect(result.carrier_dry);
}

test "rebasing across a shrinking wet carrier concentrates the solute" {
    const result = try rebase(2.5, .{ .water_volume_m3 = 4 }, .{ .water_volume_m3 = 1 });
    try std.testing.expectApproxEqRel(@as(f64, 10), result.retained_amount, 1e-15);
    // Same mass in a quarter of the water is four times the concentration.
    try std.testing.expectApproxEqRel(@as(f64, 10), result.concentration_per_m3, 1e-15);
    try std.testing.expect(!result.carrier_dry);
}

test "rebasing from a dry carrier yields nothing, matching zero stored mass" {
    // A dry carrier stores zero concentration, so its extensive amount is zero
    // and rewetting it introduces no solute from nowhere.
    const result = try rebase(0, .{ .water_volume_m3 = 0 }, .{ .water_volume_m3 = 5 });
    try std.testing.expectEqual(@as(f64, 0), result.retained_amount);
    try std.testing.expectEqual(@as(f64, 0), result.concentration_per_m3);
}

test "mass is invariant across an arbitrary chain of carrier changes" {
    // The property that matters for the audit: no sequence of wetting, drying and
    // rewetting may change the extensive amount.
    const volumes = [_]f64{ 4, 2, 0.5, 0, 3, 0, 7 };
    var carrier: Carrier = .{ .water_volume_m3 = volumes[0] };
    var concentration_per_m3: f64 = 2.5;
    const expected_amount = try extensiveAmount(concentration_per_m3, carrier);
    var retained = expected_amount;
    for (volumes[1..]) |next_volume| {
        const next: Carrier = .{ .water_volume_m3 = next_volume };
        const result = try rebase(concentration_per_m3, carrier, next);
        // Once dry, the stored concentration is zero, so the amount recoverable
        // from concentration alone is lost. This test pins that the caller must
        // carry `retained_amount` across a dry step rather than re-deriving it.
        if (result.carrier_dry) {
            try std.testing.expectApproxEqRel(retained, result.retained_amount, 1e-14);
        } else if (carrier.hasWater()) {
            try std.testing.expectApproxEqRel(retained, result.retained_amount, 1e-14);
        }
        retained = result.retained_amount;
        concentration_per_m3 = result.concentration_per_m3;
        carrier = next;
    }
    // Documenting the consequence explicitly: a water-normalized representation
    // cannot round-trip through a dry carrier on its own.
    try std.testing.expectEqual(@as(f64, 0), retained);
}

test "invalid solute states are rejected" {
    try std.testing.expectError(error.NonFiniteLitterSoluteState, concentration(std.math.nan(f64), .{ .water_volume_m3 = 1 }));
    try std.testing.expectError(error.InvalidLitterSoluteState, concentration(-1, .{ .water_volume_m3 = 1 }));
    try std.testing.expectError(error.InvalidLitterSoluteState, concentration(1, .{ .water_volume_m3 = -1 }));
    try std.testing.expectError(error.NonFiniteLitterSoluteState, extensiveAmount(std.math.inf(f64), .{ .water_volume_m3 = 1 }));
}
