const std = @import("std");

pub const ElementPool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Result = struct {
    dissolved_delta: ElementPool,
    allocated_after_dissolved: ElementPool,
};

/// HOUR1 lines 832--850. Computes DOC, DON, and DOP routed from an organic
/// amendment after microbial allocation, preserving source element order.
pub fn compute(amendment: ElementPool, allocated_before: ElementPool) !Result {
    try validatePool(amendment);
    try validatePool(allocated_before);
    if (allocated_before.carbon_g_c > amendment.carbon_g_c or
        allocated_before.nitrogen_g_n > amendment.nitrogen_g_n or
        allocated_before.phosphorus_g_p > amendment.phosphorus_g_p)
        return error.OrganicMicrobialAllocationExceedsAmendment;

    const dissolved: ElementPool = .{
        .carbon_g_c = @min(
            0.1 * allocated_before.carbon_g_c,
            amendment.carbon_g_c - allocated_before.carbon_g_c,
        ),
        .nitrogen_g_n = @min(
            0.1 * allocated_before.nitrogen_g_n,
            amendment.nitrogen_g_n - allocated_before.nitrogen_g_n,
        ),
        .phosphorus_g_p = @min(
            0.1 * allocated_before.phosphorus_g_p,
            amendment.phosphorus_g_p - allocated_before.phosphorus_g_p,
        ),
    };
    var allocated_after = allocated_before;
    allocated_after.carbon_g_c += dissolved.carbon_g_c;
    allocated_after.nitrogen_g_n += dissolved.nitrogen_g_n;
    allocated_after.phosphorus_g_p += dissolved.phosphorus_g_p;
    try validatePool(dissolved);
    try validatePool(allocated_after);
    return .{
        .dissolved_delta = dissolved,
        .allocated_after_dissolved = allocated_after,
    };
}

fn validatePool(pool: ElementPool) !void {
    inline for (@typeInfo(ElementPool).@"struct".fields) |field| {
        const value = @field(pool, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteOrganicDissolvedAllocationInput;
        if (value < 0)
            return error.InvalidOrganicDissolvedAllocationInput;
    }
}

test "dissolved allocation preserves carbon nitrogen phosphorus source order" {
    const result = try compute(
        .{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 1 },
        .{ .carbon_g_c = 20, .nitrogen_g_n = 9.5, .phosphorus_g_p = 0.2 },
    );
    try std.testing.expectEqual(@as(f64, 2), result.dissolved_delta.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0.5), result.dissolved_delta.nitrogen_g_n);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.02),
        result.dissolved_delta.phosphorus_g_p,
        1e-15,
    );
    try std.testing.expectEqual(@as(f64, 22), result.allocated_after_dissolved.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 10), result.allocated_after_dissolved.nitrogen_g_n);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.22),
        result.allocated_after_dissolved.phosphorus_g_p,
        1e-15,
    );
}

test "dissolved allocation fails when microbial allocation overdraws amendment" {
    try std.testing.expectError(
        error.OrganicMicrobialAllocationExceedsAmendment,
        compute(
            .{ .carbon_g_c = 1, .nitrogen_g_n = 1, .phosphorus_g_p = 1 },
            .{ .carbon_g_c = 1.1, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        ),
    );
}
