const std = @import("std");

pub const ElementPool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const NetBiomeProductivityTreatment = enum {
    include_organic_carbon_input,
    exclude_organic_carbon_input,
};

pub const Ledgers = struct {
    total_organic_fertilizer: ElementPool,
    cell_organic_inputs: ElementPool,
    net_biome_productivity_g_c: f64,
};

/// HOUR1 lines 895--903. Applies one organic amendment to routine-total and
/// cell input ledgers, then conditionally credits its carbon to TNBP.
pub fn apply(
    ledgers: *Ledgers,
    amendment: ElementPool,
    nbp_treatment: NetBiomeProductivityTreatment,
) !void {
    try validatePool(amendment);
    try validatePool(ledgers.total_organic_fertilizer);
    try validatePool(ledgers.cell_organic_inputs);
    if (!std.math.isFinite(ledgers.net_biome_productivity_g_c))
        return error.NonFiniteOrganicFertilizerLedger;

    var updated = ledgers.*;
    updated.total_organic_fertilizer.carbon_g_c += amendment.carbon_g_c;
    updated.total_organic_fertilizer.nitrogen_g_n += amendment.nitrogen_g_n;
    updated.total_organic_fertilizer.phosphorus_g_p += amendment.phosphorus_g_p;
    updated.cell_organic_inputs.carbon_g_c += amendment.carbon_g_c;
    updated.cell_organic_inputs.nitrogen_g_n += amendment.nitrogen_g_n;
    updated.cell_organic_inputs.phosphorus_g_p += amendment.phosphorus_g_p;
    if (nbp_treatment == .include_organic_carbon_input)
        updated.net_biome_productivity_g_c += amendment.carbon_g_c;
    try validatePool(updated.total_organic_fertilizer);
    try validatePool(updated.cell_organic_inputs);
    if (!std.math.isFinite(updated.net_biome_productivity_g_c))
        return error.NonFiniteOrganicFertilizerLedger;
    ledgers.* = updated;
}

fn validatePool(pool: ElementPool) !void {
    inline for (@typeInfo(ElementPool).@"struct".fields) |field| {
        const value = @field(pool, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteOrganicFertilizerLedger;
        if (value < 0)
            return error.InvalidOrganicFertilizerLedger;
    }
}

test "organic amendment updates all ledgers in source order" {
    var ledgers: Ledgers = .{
        .total_organic_fertilizer = .{
            .carbon_g_c = 1,
            .nitrogen_g_n = 2,
            .phosphorus_g_p = 3,
        },
        .cell_organic_inputs = .{
            .carbon_g_c = 4,
            .nitrogen_g_n = 5,
            .phosphorus_g_p = 6,
        },
        .net_biome_productivity_g_c = 7,
    };
    try apply(
        &ledgers,
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .include_organic_carbon_input,
    );
    try std.testing.expectEqual(@as(f64, 11), ledgers.total_organic_fertilizer.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 3), ledgers.total_organic_fertilizer.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 3.1), ledgers.total_organic_fertilizer.phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 14), ledgers.cell_organic_inputs.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 6), ledgers.cell_organic_inputs.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 6.1), ledgers.cell_organic_inputs.phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 17), ledgers.net_biome_productivity_g_c);
}

test "excluded amendment leaves net biome productivity unchanged" {
    var ledgers: Ledgers = .{
        .total_organic_fertilizer = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .cell_organic_inputs = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .net_biome_productivity_g_c = 5,
    };
    try apply(
        &ledgers,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 1, .phosphorus_g_p = 0.5 },
        .exclude_organic_carbon_input,
    );
    try std.testing.expectEqual(@as(f64, 2), ledgers.cell_organic_inputs.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 5), ledgers.net_biome_productivity_g_c);
}

test "invalid amendment leaves ledgers unchanged" {
    var ledgers: Ledgers = .{
        .total_organic_fertilizer = .{ .carbon_g_c = 1, .nitrogen_g_n = 2, .phosphorus_g_p = 3 },
        .cell_organic_inputs = .{ .carbon_g_c = 4, .nitrogen_g_n = 5, .phosphorus_g_p = 6 },
        .net_biome_productivity_g_c = 7,
    };
    try std.testing.expectError(
        error.InvalidOrganicFertilizerLedger,
        apply(
            &ledgers,
            .{ .carbon_g_c = -1, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
            .include_organic_carbon_input,
        ),
    );
    try std.testing.expectEqual(@as(f64, 1), ledgers.total_organic_fertilizer.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 7), ledgers.net_biome_productivity_g_c);
}
