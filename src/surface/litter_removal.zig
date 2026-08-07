const std = @import("std");

/// Runtime-sized authoritative surface pools bound from the REDIST operation
/// 21 owners. Colonized carbon is a diagnostic subset of structural carbon,
/// so it is scaled but deliberately excluded from independent C accounting.
pub const Pools = struct {
    organic_carbon_g_c: []f64,
    organic_nitrogen_g_n: []f64,
    organic_phosphorus_g_p: []f64,
    charcoal_carbon_g_c: []f64,
    charcoal_nitrogen_g_n: []f64,
    charcoal_phosphorus_g_p: []f64,
    colonized_structural_carbon_g_c: []f64,
    counted_mineral_nitrogen_g_n: []f64,
    counted_phosphate_phosphorus_g_p: []f64,
    other_scaled_nitrogen_g_n: []f64,
};

pub const Inputs = struct {
    removal_fraction: f64,
    surface_temperature_k: f64,
    /// ORGCX-equivalent dry organic C before this removal transaction.
    dry_organic_carbon_before_g_c: f64,
    dry_organic_heat_capacity_megajoules_per_g_c_k: f64 = 2.496e-6,
};

pub const Accounting = struct {
    carbon_output_g_c: f64,
    nitrogen_output_g_n: f64,
    phosphorus_output_g_p: f64,
    heat_output_megajoules: f64,
    remaining_organic_carbon_g_c: f64,
    remaining_organic_nitrogen_g_n: f64,
    remaining_organic_phosphorus_g_p: f64,
    remaining_charcoal_carbon_g_c: f64,
    remaining_charcoal_nitrogen_g_n: f64,
    remaining_charcoal_phosphorus_g_p: f64,
};

/// Atomically applies REDIST operation 21. Callers bind every named surface
/// owner into the appropriate runtime slice. Validation and all derived
/// accounting finish before the first pool is changed.
pub fn apply(pools: Pools, inputs: Inputs) !Accounting {
    try validateInputs(inputs);
    inline for (std.meta.fields(Pools)) |field|
        try validatePool(@field(pools, field.name));

    const organic_c_before = try sumFinite(pools.organic_carbon_g_c);
    const organic_n_before = try sumFinite(pools.organic_nitrogen_g_n);
    const organic_p_before = try sumFinite(pools.organic_phosphorus_g_p);
    const charcoal_c_before = try sumFinite(pools.charcoal_carbon_g_c);
    const charcoal_n_before = try sumFinite(pools.charcoal_nitrogen_g_n);
    const charcoal_p_before = try sumFinite(pools.charcoal_phosphorus_g_p);
    const mineral_n_before =
        try sumFinite(pools.counted_mineral_nitrogen_g_n);
    const mineral_p_before =
        try sumFinite(pools.counted_phosphate_phosphorus_g_p);

    const retained_fraction = 1 - inputs.removal_fraction;
    const accounting: Accounting = .{
        .carbon_output_g_c = try productFinite(
            inputs.removal_fraction,
            try addFinite(organic_c_before, charcoal_c_before),
        ),
        .nitrogen_output_g_n = try productFinite(
            inputs.removal_fraction,
            try addFinite(
                try addFinite(organic_n_before, charcoal_n_before),
                mineral_n_before,
            ),
        ),
        .phosphorus_output_g_p = try productFinite(
            inputs.removal_fraction,
            try addFinite(
                try addFinite(organic_p_before, charcoal_p_before),
                mineral_p_before,
            ),
        ),
        .remaining_organic_carbon_g_c = try productFinite(retained_fraction, organic_c_before),
        .remaining_organic_nitrogen_g_n = try productFinite(retained_fraction, organic_n_before),
        .remaining_organic_phosphorus_g_p = try productFinite(retained_fraction, organic_p_before),
        .remaining_charcoal_carbon_g_c = try productFinite(retained_fraction, charcoal_c_before),
        .remaining_charcoal_nitrogen_g_n = try productFinite(retained_fraction, charcoal_n_before),
        .remaining_charcoal_phosphorus_g_p = try productFinite(retained_fraction, charcoal_p_before),
        .heat_output_megajoules = undefined,
    };
    const remaining_dry_carbon_g_c = try addFinite(
        accounting.remaining_organic_carbon_g_c,
        accounting.remaining_charcoal_carbon_g_c,
    );
    const dry_carbon_loss_g_c =
        inputs.dry_organic_carbon_before_g_c - remaining_dry_carbon_g_c;
    if (!std.math.isFinite(dry_carbon_loss_g_c) or dry_carbon_loss_g_c < 0)
        return error.InconsistentSurfaceOrganicCarbon;
    var result = accounting;
    result.heat_output_megajoules = try productFinite(
        try productFinite(
            inputs.dry_organic_heat_capacity_megajoules_per_g_c_k,
            dry_carbon_loss_g_c,
        ),
        inputs.surface_temperature_k,
    );

    scale(pools.organic_carbon_g_c, retained_fraction);
    scale(pools.organic_nitrogen_g_n, retained_fraction);
    scale(pools.organic_phosphorus_g_p, retained_fraction);
    scale(pools.charcoal_carbon_g_c, retained_fraction);
    scale(pools.charcoal_nitrogen_g_n, retained_fraction);
    scale(pools.charcoal_phosphorus_g_p, retained_fraction);
    scale(pools.colonized_structural_carbon_g_c, retained_fraction);
    scale(pools.counted_mineral_nitrogen_g_n, retained_fraction);
    scale(pools.counted_phosphate_phosphorus_g_p, retained_fraction);
    scale(pools.other_scaled_nitrogen_g_n, retained_fraction);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.removal_fraction,
        inputs.surface_temperature_k,
        inputs.dry_organic_carbon_before_g_c,
        inputs.dry_organic_heat_capacity_megajoules_per_g_c_k,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSurfaceLitterRemovalInput;
    if (inputs.removal_fraction < 0 or inputs.removal_fraction > 0.999)
        return error.InvalidSurfaceLitterRemovalFraction;
    if (inputs.surface_temperature_k <= 0 or
        inputs.dry_organic_carbon_before_g_c < 0 or
        inputs.dry_organic_heat_capacity_megajoules_per_g_c_k <= 0)
        return error.InvalidSurfaceLitterRemovalInput;
}

fn validatePool(pool: []const f64) !void {
    for (pool) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSurfaceLitterRemovalPool;
}

fn sumFinite(values: []const f64) !f64 {
    var total: f64 = 0;
    for (values) |value| {
        total += value;
        if (!std.math.isFinite(total))
            return error.SurfaceLitterRemovalAccountingOverflow;
    }
    return total;
}

fn addFinite(a: f64, b: f64) !f64 {
    const result = a + b;
    if (!std.math.isFinite(result))
        return error.SurfaceLitterRemovalAccountingOverflow;
    return result;
}

fn productFinite(a: f64, b: f64) !f64 {
    const result = a * b;
    if (!std.math.isFinite(result))
        return error.SurfaceLitterRemovalAccountingOverflow;
    return result;
}

fn scale(values: []f64, retained_fraction: f64) void {
    for (values) |*value| value.* *= retained_fraction;
}

test "operation 21 scales arbitrary runtime pools and accounts C N P heat" {
    var organic_c = [_]f64{ 10, 20, 30, 40 };
    var organic_n = [_]f64{ 1, 2, 3 };
    var organic_p = [_]f64{ 0.5, 1.5 };
    var charcoal_c = [_]f64{50};
    var charcoal_n = [_]f64{5};
    var charcoal_p = [_]f64{2};
    var colonized_c = [_]f64{ 8, 12 };
    var mineral_n = [_]f64{ 4, 6, 8, 2 };
    var mineral_p = [_]f64{ 3, 7 };
    var other_n = [_]f64{ 9, 11, 13, 17, 19 };
    const result = try apply(.{
        .organic_carbon_g_c = &organic_c,
        .organic_nitrogen_g_n = &organic_n,
        .organic_phosphorus_g_p = &organic_p,
        .charcoal_carbon_g_c = &charcoal_c,
        .charcoal_nitrogen_g_n = &charcoal_n,
        .charcoal_phosphorus_g_p = &charcoal_p,
        .colonized_structural_carbon_g_c = &colonized_c,
        .counted_mineral_nitrogen_g_n = &mineral_n,
        .counted_phosphate_phosphorus_g_p = &mineral_p,
        .other_scaled_nitrogen_g_n = &other_n,
    }, .{
        .removal_fraction = 0.25,
        .surface_temperature_k = 300,
        .dry_organic_carbon_before_g_c = 150,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 37.5), result.carbon_output_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7.75), result.nitrogen_output_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.5), result.phosphorus_output_g_p, 1e-12);
    try std.testing.expectApproxEqAbs(
        2.496e-6 * 37.5 * 300,
        result.heat_output_megajoules,
        1e-14,
    );
    try std.testing.expectEqual(@as(f64, 7.5), organic_c[0]);
    try std.testing.expectEqual(@as(f64, 37.5), charcoal_c[0]);
    try std.testing.expectEqual(@as(f64, 6), colonized_c[0]);
    try std.testing.expectEqual(@as(f64, 3), mineral_n[0]);
    try std.testing.expectEqual(@as(f64, 6.75), other_n[0]);
}

test "late invalid pool rolls back every operation 21 owner" {
    var organic_c = [_]f64{ 10, 20 };
    const organic_c_before = organic_c;
    var empty = [_]f64{};
    var other_n = [_]f64{ 1, std.math.nan(f64) };
    try std.testing.expectError(
        error.InvalidSurfaceLitterRemovalPool,
        apply(.{
            .organic_carbon_g_c = &organic_c,
            .organic_nitrogen_g_n = &empty,
            .organic_phosphorus_g_p = &empty,
            .charcoal_carbon_g_c = &empty,
            .charcoal_nitrogen_g_n = &empty,
            .charcoal_phosphorus_g_p = &empty,
            .colonized_structural_carbon_g_c = &empty,
            .counted_mineral_nitrogen_g_n = &empty,
            .counted_phosphate_phosphorus_g_p = &empty,
            .other_scaled_nitrogen_g_n = &other_n,
        }, .{
            .removal_fraction = 0.5,
            .surface_temperature_k = 290,
            .dry_organic_carbon_before_g_c = 30,
        }),
    );
    try std.testing.expectEqualSlices(f64, &organic_c_before, &organic_c);
    try std.testing.expectEqual(@as(f64, 1), other_n[0]);
    try std.testing.expect(std.math.isNan(other_n[1]));
}

test "inconsistent heat baseline fails before pool mutation" {
    var organic_c = [_]f64{100};
    var empty = [_]f64{};
    try std.testing.expectError(
        error.InconsistentSurfaceOrganicCarbon,
        apply(.{
            .organic_carbon_g_c = &organic_c,
            .organic_nitrogen_g_n = &empty,
            .organic_phosphorus_g_p = &empty,
            .charcoal_carbon_g_c = &empty,
            .charcoal_nitrogen_g_n = &empty,
            .charcoal_phosphorus_g_p = &empty,
            .colonized_structural_carbon_g_c = &empty,
            .counted_mineral_nitrogen_g_n = &empty,
            .counted_phosphate_phosphorus_g_p = &empty,
            .other_scaled_nitrogen_g_n = &empty,
        }, .{
            .removal_fraction = 0.1,
            .surface_temperature_k = 290,
            .dry_organic_carbon_before_g_c = 50,
        }),
    );
    try std.testing.expectEqual(@as(f64, 100), organic_c[0]);
}
