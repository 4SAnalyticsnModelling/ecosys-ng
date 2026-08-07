const std = @import("std");

pub const salt_family_count = 42;
pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    simulation: SaltSimulation,
    layer: usize,
    layer_count: usize,
    incorporation_fraction: f64, // FI, dimensionless
    /// REDIST 12763--12804 source order. Each slice is indexed by soil layer.
    /// The source deliberately maps ZALOH1--4 to TZFEOH1--4GS at 12775--12778.
    layer_amounts: [salt_family_count][]f64,
    mixed_totals: [salt_family_count]f64,
};

fn finite(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of the ISALTG guard and REDIST 12762--12805.
pub fn incorporate(allocator: std.mem.Allocator, inputs: Inputs) !void {
    if (inputs.simulation == .disabled) return;
    if (inputs.layer_count == 0 or inputs.layer >= inputs.layer_count) return error.TillageSaltIncorporationDimensionMismatch;
    if (!std.math.isFinite(inputs.incorporation_fraction) or !finite(&inputs.mixed_totals)) return error.InvalidTillageSaltIncorporationInput;
    for (inputs.layer_amounts) |amounts| {
        if (amounts.len != inputs.layer_count) return error.TillageSaltIncorporationDimensionMismatch;
        if (!finite(amounts)) return error.InvalidTillageSaltIncorporationInput;
    }

    const staged = try allocator.alloc(f64, salt_family_count);
    defer allocator.free(staged);
    for (inputs.layer_amounts, inputs.mixed_totals, 0..) |amounts, total, family| {
        staged[family] = amounts[inputs.layer] + inputs.incorporation_fraction * total;
        if (!std.math.isFinite(staged[family])) return error.NonFiniteTillageSaltIncorporationResult;
    }
    for (inputs.layer_amounts, 0..) |amounts, family| amounts[inputs.layer] = staged[family];
}

test "REDIST salt incorporation follows all 42 source-ordered families" {
    var storage: [salt_family_count][2]f64 = @splat(.{ 2, 3 });
    var amounts: [salt_family_count][]f64 = undefined;
    for (0..salt_family_count) |family| amounts[family] = &storage[family];
    const totals: [salt_family_count]f64 = @splat(4);
    try incorporate(std.testing.allocator, .{
        .simulation = .enabled,
        .layer = 1,
        .layer_count = 2,
        .incorporation_fraction = 0.5,
        .layer_amounts = amounts,
        .mixed_totals = totals,
    });
    for (storage) |family| {
        try std.testing.expectEqual(@as(f64, 2), family[0]);
        try std.testing.expectEqual(@as(f64, 5), family[1]);
    }
}

test "REDIST disabled salt simulation preserves legacy skip semantics" {
    const invalid: [salt_family_count][]f64 = @splat(&.{});
    const totals: [salt_family_count]f64 = @splat(std.math.nan(f64));
    try incorporate(std.testing.allocator, .{
        .simulation = .disabled,
        .layer = 99,
        .layer_count = 0,
        .incorporation_fraction = std.math.nan(f64),
        .layer_amounts = invalid,
        .mixed_totals = totals,
    });
}

test "REDIST salt incorporation rejects overflow atomically" {
    var storage: [salt_family_count][1]f64 = @splat(.{1});
    var amounts: [salt_family_count][]f64 = undefined;
    for (0..salt_family_count) |family| amounts[family] = &storage[family];
    var totals: [salt_family_count]f64 = @splat(1);
    totals[41] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteTillageSaltIncorporationResult, incorporate(std.testing.allocator, .{
        .simulation = .enabled,
        .layer = 0,
        .layer_count = 1,
        .incorporation_fraction = 2,
        .layer_amounts = amounts,
        .mixed_totals = totals,
    }));
    for (storage) |family| try std.testing.expectEqual(@as(f64, 1), family[0]);
}
