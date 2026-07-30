const std = @import("std");

pub const Domain = enum {
    carbon,
    nitrogen,
    phosphorus,
    water,
    heat,
    oxygen,
    ions,
};

pub const Deviations = struct {
    carbon_g_per_m2: f64,
    nitrogen_g_per_m2: f64,
    phosphorus_g_per_m2: f64,
    water_m: f64,
    heat_mj_per_m2: f64,
    oxygen_g_per_m2: f64,
    ions_mol_per_m2: f64,
};

pub const Tolerances = struct {
    carbon_g_per_m2: f64 = 1e-6,
    nitrogen_g_per_m2: f64 = 1e-6,
    phosphorus_g_per_m2: f64 = 1e-6,
    water_m: f64 = 1e-6,
    heat_mj_per_m2: f64 = 1e-6,
    oxygen_g_per_m2: f64 = 1e-6,
    ions_mol_per_m2: f64 = 1e-6,
};

pub const Assessment = struct {
    failed: [7]bool,
    first_failure: ?Domain,

    pub fn hasFailure(self: Assessment) bool {
        return self.first_failure != null;
    }
};

/// Classifies every EXEC loss using the source's strict ABS(DIFF)>tolerance
/// comparisons and source priority C,N,P,water,heat,O2,ions. It does not reset
/// baselines: production callers terminate with the mapped informative error.
pub fn assess(
    deviations: Deviations,
    tolerances: Tolerances,
) !Assessment {
    const values = [_]f64{
        deviations.carbon_g_per_m2,
        deviations.nitrogen_g_per_m2,
        deviations.phosphorus_g_per_m2,
        deviations.water_m,
        deviations.heat_mj_per_m2,
        deviations.oxygen_g_per_m2,
        deviations.ions_mol_per_m2,
    };
    const limits = [_]f64{
        tolerances.carbon_g_per_m2,
        tolerances.nitrogen_g_per_m2,
        tolerances.phosphorus_g_per_m2,
        tolerances.water_m,
        tolerances.heat_mj_per_m2,
        tolerances.oxygen_g_per_m2,
        tolerances.ions_mol_per_m2,
    };
    var result: Assessment = .{
        .failed = [_]bool{false} ** 7,
        .first_failure = null,
    };
    for (values, limits, 0..) |value, limit, index| {
        if (!std.math.isFinite(value))
            return error.NonFiniteExecutionBalanceDeviation;
        if (!std.math.isFinite(limit) or limit < 0)
            return error.InvalidExecutionBalanceTolerance;
        result.failed[index] = @abs(value) > limit;
        if (result.failed[index] and result.first_failure == null)
            result.first_failure = @enumFromInt(index);
    }
    return result;
}

pub fn failureError(domain: Domain) anyerror {
    return switch (domain) {
        .carbon => error.CarbonMassBalanceLost,
        .nitrogen => error.NitrogenMassBalanceLost,
        .phosphorus => error.PhosphorusMassBalanceLost,
        .water => error.WaterMassBalanceLost,
        .heat => error.HeatMassBalanceLost,
        .oxygen => error.OxygenMassBalanceLost,
        .ions => error.IonMassBalanceLost,
    };
}

fn zeros() Deviations {
    return .{
        .carbon_g_per_m2 = 0,
        .nitrogen_g_per_m2 = 0,
        .phosphorus_g_per_m2 = 0,
        .water_m = 0,
        .heat_mj_per_m2 = 0,
        .oxygen_g_per_m2 = 0,
        .ions_mol_per_m2 = 0,
    };
}

test "EXEC threshold is strict and symmetric" {
    var deviations = zeros();
    deviations.carbon_g_per_m2 = 1e-6;
    var result = try assess(deviations, .{});
    try std.testing.expect(!result.hasFailure());
    deviations.carbon_g_per_m2 = -1.000001e-6;
    result = try assess(deviations, .{});
    try std.testing.expectEqual(Domain.carbon, result.first_failure.?);
    try std.testing.expectEqual(
        error.CarbonMassBalanceLost,
        failureError(result.first_failure.?),
    );
}

test "simultaneous failures retain every domain and source priority" {
    var deviations = zeros();
    deviations.nitrogen_g_per_m2 = 2e-6;
    deviations.water_m = -3e-6;
    deviations.ions_mol_per_m2 = 4e-6;
    const result = try assess(deviations, .{});
    try std.testing.expectEqual(Domain.nitrogen, result.first_failure.?);
    try std.testing.expect(!result.failed[@intFromEnum(Domain.carbon)]);
    try std.testing.expect(result.failed[@intFromEnum(Domain.nitrogen)]);
    try std.testing.expect(result.failed[@intFromEnum(Domain.water)]);
    try std.testing.expect(result.failed[@intFromEnum(Domain.ions)]);
}

test "unit-specific runtime tolerances do not alias domains" {
    var deviations = zeros();
    deviations.water_m = 0.5;
    deviations.heat_mj_per_m2 = 0.5;
    const result = try assess(deviations, .{
        .water_m = 1,
        .heat_mj_per_m2 = 0.1,
    });
    try std.testing.expectEqual(Domain.heat, result.first_failure.?);
    try std.testing.expect(!result.failed[@intFromEnum(Domain.water)]);
    try std.testing.expect(result.failed[@intFromEnum(Domain.heat)]);
}

test "nonfinite late domain aborts classification" {
    var deviations = zeros();
    deviations.carbon_g_per_m2 = 2e-6;
    deviations.ions_mol_per_m2 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteExecutionBalanceDeviation,
        assess(deviations, .{}),
    );
    try std.testing.expectError(
        error.InvalidExecutionBalanceTolerance,
        assess(zeros(), .{ .oxygen_g_per_m2 = -1 }),
    );
}
