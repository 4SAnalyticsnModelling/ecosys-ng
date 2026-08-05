const std = @import("std");

pub const Snapshot = struct {
    leaf_carbon_g_c: f64,
    leaf_nitrogen_g_n: f64,
    leaf_phosphorus_g_p: f64,
    remobilizable_carbon_g_c: f64,
    remobilizable_nitrogen_g_n: f64,
    remobilizable_phosphorus_g_p: f64,
    remobilization_fraction: f64,
};

pub const Composition = struct {
    woody_carbon_fraction: [2]f64,
    woody_nitrogen_fraction: [2]f64,
    woody_phosphorus_fraction: [2]f64,
};

pub const KineticFractions = struct {
    woody_carbon: []const f64,
    woody_nitrogen: []const f64,
    woody_phosphorus: []const f64,
    leaf_carbon: []const f64,
    leaf_nitrogen: []const f64,
    leaf_phosphorus: []const f64,
};

pub const Accumulators = struct {
    woody_carbon_g_c: []f64,
    woody_nitrogen_g_n: []f64,
    woody_phosphorus_g_p: []f64,
    leaf_carbon_g_c: []f64,
    leaf_nitrogen_g_n: []f64,
    leaf_phosphorus_g_p: []f64,
};

/// GROSUB lines 2603--2616. Publishes the selected leaf's woody and nonwoody
/// litter into a runtime number of kinetic pools. Within each pool, the exact
/// source order is woody C/N/P followed by foliar C/N/P.
pub fn publish(
    snapshot: Snapshot,
    composition: Composition,
    kinetics: KineticFractions,
    accumulators: Accumulators,
) !void {
    const kinetic_count = kinetics.woody_carbon.len;
    if (kinetic_count == 0) return error.ZeroLeafSenescenceKineticPools;
    inline for (@typeInfo(KineticFractions).@"struct".fields[1..]) |field|
        if (@field(kinetics, field.name).len != kinetic_count)
            return error.LeafSenescenceKineticDimensionMismatch;
    inline for (@typeInfo(Accumulators).@"struct".fields) |field|
        if (@field(accumulators, field.name).len != kinetic_count)
            return error.LeafSenescenceKineticDimensionMismatch;

    inline for (@typeInfo(Snapshot).@"struct".fields) |field| {
        const value = @field(snapshot, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLeafSenescenceLitterSnapshot;
    }
    if (snapshot.remobilization_fraction > 1 or
        snapshot.remobilizable_carbon_g_c > snapshot.leaf_carbon_g_c or
        snapshot.remobilizable_nitrogen_g_n > snapshot.leaf_nitrogen_g_n or
        snapshot.remobilizable_phosphorus_g_p > snapshot.leaf_phosphorus_g_p)
        return error.InvalidLeafSenescenceLitterSnapshot;
    inline for (.{
        composition.woody_carbon_fraction,
        composition.woody_nitrogen_fraction,
        composition.woody_phosphorus_fraction,
    }) |fractions| {
        inline for (fractions) |fraction|
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidLeafSenescenceComposition;
        if (@abs(fractions[0] + fractions[1] - 1) > 1.0e-12)
            return error.InvalidLeafSenescenceComposition;
    }
    try validateKinetics(kinetics, kinetic_count);

    for (0..kinetic_count) |kinetic| {
        const additions = additionsForPool(snapshot, composition, kinetics, kinetic);
        inline for (@typeInfo(Accumulators).@"struct".fields) |field| {
            const current = @field(accumulators, field.name)[kinetic];
            const updated = current + @field(additions, field.name);
            if (!std.math.isFinite(current) or current < 0 or
                !std.math.isFinite(updated) or updated < 0)
                return error.NonFiniteLeafSenescenceLitterResult;
        }
    }

    for (0..kinetic_count) |kinetic| {
        const additions = additionsForPool(snapshot, composition, kinetics, kinetic);
        accumulators.woody_carbon_g_c[kinetic] += additions.woody_carbon_g_c;
        accumulators.woody_nitrogen_g_n[kinetic] += additions.woody_nitrogen_g_n;
        accumulators.woody_phosphorus_g_p[kinetic] += additions.woody_phosphorus_g_p;
        accumulators.leaf_carbon_g_c[kinetic] += additions.leaf_carbon_g_c;
        accumulators.leaf_nitrogen_g_n[kinetic] += additions.leaf_nitrogen_g_n;
        accumulators.leaf_phosphorus_g_p[kinetic] += additions.leaf_phosphorus_g_p;
    }
}

fn validateKinetics(kinetics: KineticFractions, kinetic_count: usize) !void {
    inline for (@typeInfo(KineticFractions).@"struct".fields) |field| {
        var total: f64 = 0;
        for (@field(kinetics, field.name)) |fraction| {
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidLeafSenescenceKineticFraction;
            total += fraction;
        }
        if (!std.math.isFinite(total) or @abs(total - 1) > 1.0e-12)
            return error.InvalidLeafSenescenceKineticFraction;
    }
    _ = kinetic_count;
}

fn additionsForPool(
    snapshot: Snapshot,
    composition: Composition,
    kinetics: KineticFractions,
    kinetic: usize,
) AccumulatorValues {
    return .{
        .woody_carbon_g_c = kinetics.woody_carbon[kinetic] *
            snapshot.remobilization_fraction * snapshot.leaf_carbon_g_c *
            composition.woody_carbon_fraction[0],
        .woody_nitrogen_g_n = kinetics.woody_nitrogen[kinetic] *
            snapshot.remobilization_fraction * snapshot.leaf_nitrogen_g_n *
            composition.woody_nitrogen_fraction[0],
        .woody_phosphorus_g_p = kinetics.woody_phosphorus[kinetic] *
            snapshot.remobilization_fraction * snapshot.leaf_phosphorus_g_p *
            composition.woody_phosphorus_fraction[0],
        .leaf_carbon_g_c = kinetics.leaf_carbon[kinetic] *
            snapshot.remobilization_fraction *
            (snapshot.leaf_carbon_g_c - snapshot.remobilizable_carbon_g_c) *
            composition.woody_carbon_fraction[1],
        .leaf_nitrogen_g_n = kinetics.leaf_nitrogen[kinetic] *
            snapshot.remobilization_fraction *
            (snapshot.leaf_nitrogen_g_n - snapshot.remobilizable_nitrogen_g_n) *
            composition.woody_nitrogen_fraction[1],
        .leaf_phosphorus_g_p = kinetics.leaf_phosphorus[kinetic] *
            snapshot.remobilization_fraction *
            (snapshot.leaf_phosphorus_g_p - snapshot.remobilizable_phosphorus_g_p) *
            composition.woody_phosphorus_fraction[1],
    };
}

const AccumulatorValues = struct {
    woody_carbon_g_c: f64,
    woody_nitrogen_g_n: f64,
    woody_phosphorus_g_p: f64,
    leaf_carbon_g_c: f64,
    leaf_nitrogen_g_n: f64,
    leaf_phosphorus_g_p: f64,
};

fn uniformKinetics(values: []const f64) KineticFractions {
    return .{
        .woody_carbon = values,
        .woody_nitrogen = values,
        .woody_phosphorus = values,
        .leaf_carbon = values,
        .leaf_nitrogen = values,
        .leaf_phosphorus = values,
    };
}

fn makeAccumulators(storage: *[6][6]f64, count: usize) Accumulators {
    return .{
        .woody_carbon_g_c = storage[0][0..count],
        .woody_nitrogen_g_n = storage[1][0..count],
        .woody_phosphorus_g_p = storage[2][0..count],
        .leaf_carbon_g_c = storage[3][0..count],
        .leaf_nitrogen_g_n = storage[4][0..count],
        .leaf_phosphorus_g_p = storage[5][0..count],
    };
}

const sample_snapshot: Snapshot = .{
    .leaf_carbon_g_c = 10,
    .leaf_nitrogen_g_n = 1,
    .leaf_phosphorus_g_p = 0.2,
    .remobilizable_carbon_g_c = 5,
    .remobilizable_nitrogen_g_n = 0.8,
    .remobilizable_phosphorus_g_p = 0.17,
    .remobilization_fraction = 0.25,
};

const sample_composition: Composition = .{
    .woody_carbon_fraction = .{ 0.2, 0.8 },
    .woody_nitrogen_fraction = .{ 0.1, 0.9 },
    .woody_phosphorus_fraction = .{ 0.3, 0.7 },
};

test "GROSUB leaf litter partition supports runtime kinetic pools" {
    const fractions = [_]f64{ 0.05, 0.1, 0.15, 0.2, 0.2, 0.3 };
    var storage: [6][6]f64 = @splat(@splat(0));
    const output = makeAccumulators(&storage, fractions.len);
    try publish(sample_snapshot, sample_composition, uniformKinetics(&fractions), output);
    try std.testing.expectApproxEqAbs(0.05 * 0.25 * 10 * 0.2, output.woody_carbon_g_c[0], 1.0e-15);
    try std.testing.expectApproxEqAbs(0.3 * 0.25 * 0.2 * 0.3, output.woody_phosphorus_g_p[5], 1.0e-15);
    try std.testing.expectApproxEqAbs(0.05 * 0.25 * (10 - 5) * 0.8, output.leaf_carbon_g_c[0], 1.0e-15);
    try std.testing.expectApproxEqAbs(0.3 * 0.25 * (1 - 0.8) * 0.9, output.leaf_nitrogen_g_n[5], 1.0e-15);
}

test "partition conserves classified non-remobilized elements" {
    const fractions = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    var storage: [6][6]f64 = @splat(@splat(0));
    const output = makeAccumulators(&storage, fractions.len);
    try publish(sample_snapshot, sample_composition, uniformKinetics(&fractions), output);
    var carbon_total: f64 = 0;
    var nitrogen_total: f64 = 0;
    var phosphorus_total: f64 = 0;
    for (0..fractions.len) |kinetic| {
        carbon_total += output.woody_carbon_g_c[kinetic] + output.leaf_carbon_g_c[kinetic];
        nitrogen_total += output.woody_nitrogen_g_n[kinetic] + output.leaf_nitrogen_g_n[kinetic];
        phosphorus_total += output.woody_phosphorus_g_p[kinetic] + output.leaf_phosphorus_g_p[kinetic];
    }
    try std.testing.expectApproxEqAbs(0.25 * (10 * 0.2 + (10 - 5) * 0.8), carbon_total, 1.0e-15);
    try std.testing.expectApproxEqAbs(0.25 * (1 * 0.1 + (1 - 0.8) * 0.9), nitrogen_total, 1.0e-15);
    try std.testing.expectApproxEqAbs(0.25 * (0.2 * 0.3 + (0.2 - 0.17) * 0.7), phosphorus_total, 1.0e-15);
}

test "late accumulator overflow leaves all pools unchanged" {
    const fractions = [_]f64{ 1, 0, 0, 0, 0, 0 };
    var storage: [6][6]f64 = @splat(@splat(1));
    storage[0][0] = std.math.floatMax(f64);
    var overflowing = sample_snapshot;
    overflowing.leaf_carbon_g_c = std.math.floatMax(f64);
    overflowing.remobilizable_carbon_g_c = 0;
    overflowing.leaf_phosphorus_g_p = std.math.floatMax(f64);
    overflowing.remobilizable_phosphorus_g_p = 0;
    overflowing.remobilization_fraction = 1;
    var all_woody = sample_composition;
    all_woody.woody_carbon_fraction = .{ 1, 0 };
    try std.testing.expectError(
        error.NonFiniteLeafSenescenceLitterResult,
        publish(overflowing, all_woody, uniformKinetics(&fractions), makeAccumulators(&storage, fractions.len)),
    );
    try std.testing.expectEqual(std.math.floatMax(f64), storage[0][0]);
    try std.testing.expectEqual(@as(f64, 1), storage[1][0]);
}

test "kinetic topology and invalid snapshot fail explicitly" {
    const fractions = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var storage: [6][6]f64 = @splat(@splat(0));
    var malformed = uniformKinetics(&fractions);
    malformed.leaf_phosphorus = malformed.leaf_phosphorus[0..3];
    try std.testing.expectError(
        error.LeafSenescenceKineticDimensionMismatch,
        publish(sample_snapshot, sample_composition, malformed, makeAccumulators(&storage, fractions.len)),
    );
    var invalid = sample_snapshot;
    invalid.remobilizable_carbon_g_c = 11;
    try std.testing.expectError(
        error.InvalidLeafSenescenceLitterSnapshot,
        publish(invalid, sample_composition, uniformKinetics(&fractions), makeAccumulators(&storage, fractions.len)),
    );
}
