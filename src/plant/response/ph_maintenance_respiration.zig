const std = @import("std");

pub const Response = struct {
    hydrogen_activity_mol_m3: f64,
    maintenance_respiration_factor: f64,
};

pub const ResponseError = error{
    NonFiniteInput,
    InvalidPh,
    InvalidDoublingActivity,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 3264--3265.
///
/// `doubling_hydrogen_activity_mol_m3` is legacy PHKI: the H+ activity at
/// which the uncapped pH contribution equals one.
pub fn calculate(
    ph: f64,
    doubling_hydrogen_activity_mol_m3: f64,
) ResponseError!Response {
    if (!std.math.isFinite(ph) or
        !std.math.isFinite(doubling_hydrogen_activity_mol_m3))
    {
        return error.NonFiniteInput;
    }
    if (ph < 0.0 or ph > 14.0) return error.InvalidPh;
    if (doubling_hydrogen_activity_mol_m3 <= 0.0) {
        return error.InvalidDoublingActivity;
    }

    const hydrogen_activity_mol_m3 =
        1.0e3 * std.math.pow(f64, 10.0, -ph);
    const maintenance_respiration_factor = 1.0 +
        @min(4.0, hydrogen_activity_mol_m3 / doubling_hydrogen_activity_mol_m3);
    if (!std.math.isFinite(hydrogen_activity_mol_m3) or
        !std.math.isFinite(maintenance_respiration_factor))
    {
        return error.NonFiniteResult;
    }
    return .{
        .hydrogen_activity_mol_m3 = hydrogen_activity_mol_m3,
        .maintenance_respiration_factor = maintenance_respiration_factor,
    };
}

test "pH response preserves hydrogen activity and uncapped factor" {
    const response = try calculate(7.0, 0.001);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0e-4),
        response.hydrogen_activity_mol_m3,
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.1),
        response.maintenance_respiration_factor,
        1.0e-15,
    );
}

test "maintenance response contribution is capped at four" {
    const response = try calculate(3.0, 1.0e-6);
    try std.testing.expectEqual(@as(f64, 5.0), response.maintenance_respiration_factor);
}

test "zero doubling activity fails before division" {
    try std.testing.expectError(error.InvalidDoublingActivity, calculate(7.0, 0.0));
}
