const std = @import("std");

pub const Inputs = struct {
    hydrogen_concentration_mol_per_m3: f64,
    hydroxide_concentration_mol_per_m3: f64,
    monovalent_activity_coefficient: f64,
    water_activity_product_mol2_per_m6: f64,
    negligible_concentration_mol_per_m3: f64,
};

pub const Result = struct {
    hydrogen_concentration_mol_per_m3: f64,
    hydroxide_concentration_mol_per_m3: f64,
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
    equal_reaction_extent_mol_per_m3: f64,
    ph: f64,
};

pub const FinalResetInputs = struct {
    initial_hydrogen_concentration_mol_per_m3: f64,
    initial_hydroxide_concentration_mol_per_m3: f64,
    accumulated_hydrogen_change_mol_per_m3: f64,
    accumulated_hydroxide_change_mol_per_m3: f64,
    monovalent_activity_coefficient: f64,
    water_activity_product_mol2_per_m6: f64,
    negligible_concentration_mol_per_m3: f64,
};

pub const SourceOrderFinalReset = struct {
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
    equal_reaction_extent_activity_mol_per_m3: f64,
    accumulated_hydrogen_change_mol_per_m3: f64,
    accumulated_hydroxide_change_mol_per_m3: f64,
    published_hydrogen_concentration_mol_per_m3: f64,
    published_hydroxide_concentration_mol_per_m3: f64,
};

pub const SourceOrderSurfaceReset = struct {
    initial_hydrogen_activity_mol_per_m3: f64,
    initial_hydroxide_activity_mol_per_m3: f64,
    activity_sum_mol_per_m3: f64,
    nonnegative_discriminant_mol2_per_m6: f64,
    equal_reaction_extent_activity_mol_per_m3: f64,
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
    hydrogen_concentration_mol_per_m3: f64,
    hydroxide_concentration_mol_per_m3: f64,
    ph: f64,
    hydrogen_floor_was_applied: bool,
    hydroxide_floor_was_applied: bool,
};

/// SOLUTE lines 834--854 H2O equilibrium reset with a cancellation-safe
/// quadratic evaluation. The reaction preserves the H+-OH- activity
/// difference. Solve the larger ion from that invariant, then obtain the
/// smaller ion by division from Kw; this remains accurate for the extreme
/// acid/base ratios that lose the smaller source-form root.
pub fn solve(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteWaterEquilibriumInput;
    if (inputs.hydrogen_concentration_mol_per_m3 < 0 or inputs.hydroxide_concentration_mol_per_m3 < 0 or inputs.monovalent_activity_coefficient <= 0 or inputs.water_activity_product_mol2_per_m6 <= 0 or inputs.negligible_concentration_mol_per_m3 < 0) return error.InvalidWaterEquilibriumInput;
    const hydrogen_activity = if (inputs.hydrogen_concentration_mol_per_m3 > inputs.negligible_concentration_mol_per_m3) inputs.hydrogen_concentration_mol_per_m3 * inputs.monovalent_activity_coefficient else inputs.hydrogen_concentration_mol_per_m3;
    const hydroxide_activity = if (inputs.hydroxide_concentration_mol_per_m3 > inputs.negligible_concentration_mol_per_m3) inputs.hydroxide_concentration_mol_per_m3 * inputs.monovalent_activity_coefficient else inputs.hydroxide_concentration_mol_per_m3;
    const activity_difference = hydrogen_activity - hydroxide_activity;
    const discriminant = std.math.hypot(
        activity_difference,
        2 * @sqrt(inputs.water_activity_product_mol2_per_m6),
    );
    const next_hydrogen_activity, const next_hydroxide_activity =
        if (activity_difference >= 0) .{
            0.5 * (activity_difference + discriminant),
            inputs.water_activity_product_mol2_per_m6 /
                (0.5 * (activity_difference + discriminant)),
        } else .{
            inputs.water_activity_product_mol2_per_m6 /
                (0.5 * (-activity_difference + discriminant)),
            0.5 * (-activity_difference + discriminant),
        };
    const extent = hydrogen_activity - next_hydrogen_activity;
    if (next_hydrogen_activity <= 0 or next_hydroxide_activity <= 0) return error.InvalidWaterEquilibriumSolution;
    const tolerance = 32 * std.math.floatEps(f64) * @max(inputs.water_activity_product_mol2_per_m6, next_hydrogen_activity * next_hydroxide_activity);
    if (@abs(next_hydrogen_activity * next_hydroxide_activity - inputs.water_activity_product_mol2_per_m6) > tolerance) return error.WaterActivityProductFailure;
    return .{
        .hydrogen_concentration_mol_per_m3 = next_hydrogen_activity / inputs.monovalent_activity_coefficient,
        .hydroxide_concentration_mol_per_m3 = next_hydroxide_activity / inputs.monovalent_activity_coefficient,
        .hydrogen_activity_mol_per_m3 = next_hydrogen_activity,
        .hydroxide_activity_mol_per_m3 = next_hydroxide_activity,
        .equal_reaction_extent_mol_per_m3 = extent,
        .ph = -@log10(next_hydrogen_activity * 1e-3),
    };
}

/// Exact equation-order diagnostic for the surface reset at SOLUTE.F
/// lines 4238--4265.
///
/// Production uses `solve`, which preserves the activity-difference
/// invariant while avoiding cancellation of the smaller quadratic root.
pub fn calculateSourceOrderSurfaceReset(
    inputs: Inputs,
) !SourceOrderSurfaceReset {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteWaterEquilibriumInput;
    }
    if (inputs.hydrogen_concentration_mol_per_m3 < 0 or
        inputs.hydroxide_concentration_mol_per_m3 < 0 or
        inputs.monovalent_activity_coefficient <= 0 or
        inputs.water_activity_product_mol2_per_m6 <= 0 or
        inputs.negligible_concentration_mol_per_m3 <= 0)
    {
        return error.InvalidWaterEquilibriumInput;
    }

    // SOLUTE.F 4248--4265. Preserve strict gates and arithmetic order.
    const initial_hydrogen_activity =
        if (inputs.hydrogen_concentration_mol_per_m3 >
        inputs.negligible_concentration_mol_per_m3)
            inputs.hydrogen_concentration_mol_per_m3 *
                inputs.monovalent_activity_coefficient
        else
            inputs.hydrogen_concentration_mol_per_m3;
    const initial_hydroxide_activity =
        if (inputs.hydroxide_concentration_mol_per_m3 >
        inputs.negligible_concentration_mol_per_m3)
            inputs.hydroxide_concentration_mol_per_m3 *
                inputs.monovalent_activity_coefficient
        else
            inputs.hydroxide_concentration_mol_per_m3;
    const activity_sum =
        initial_hydrogen_activity + initial_hydroxide_activity;
    const discriminant = @max(
        0.0,
        activity_sum * activity_sum -
            4.0 *
                (initial_hydrogen_activity *
                    initial_hydroxide_activity -
                    inputs.water_activity_product_mol2_per_m6),
    );
    const extent =
        0.5 * (activity_sum - @sqrt(discriminant));
    const unconstrained_hydrogen_activity =
        initial_hydrogen_activity - extent;
    const unconstrained_hydroxide_activity =
        initial_hydroxide_activity - extent;
    const hydrogen_activity = @max(
        inputs.negligible_concentration_mol_per_m3,
        unconstrained_hydrogen_activity,
    );
    const hydroxide_activity = @max(
        inputs.negligible_concentration_mol_per_m3,
        unconstrained_hydroxide_activity,
    );
    const hydrogen_concentration =
        hydrogen_activity / inputs.monovalent_activity_coefficient;
    const hydroxide_concentration =
        hydroxide_activity / inputs.monovalent_activity_coefficient;
    const ph = -@log10(hydrogen_activity * 1.0e-3);

    const result: SourceOrderSurfaceReset = .{
        .initial_hydrogen_activity_mol_per_m3 = initial_hydrogen_activity,
        .initial_hydroxide_activity_mol_per_m3 = initial_hydroxide_activity,
        .activity_sum_mol_per_m3 = activity_sum,
        .nonnegative_discriminant_mol2_per_m6 = discriminant,
        .equal_reaction_extent_activity_mol_per_m3 = extent,
        .hydrogen_activity_mol_per_m3 = hydrogen_activity,
        .hydroxide_activity_mol_per_m3 = hydroxide_activity,
        .hydrogen_concentration_mol_per_m3 = hydrogen_concentration,
        .hydroxide_concentration_mol_per_m3 = hydroxide_concentration,
        .ph = ph,
        .hydrogen_floor_was_applied = unconstrained_hydrogen_activity <
            inputs.negligible_concentration_mol_per_m3,
        .hydroxide_floor_was_applied = unconstrained_hydroxide_activity <
            inputs.negligible_concentration_mol_per_m3,
    };
    inline for (@typeInfo(SourceOrderSurfaceReset).@"struct".fields) |field| {
        switch (@typeInfo(field.type)) {
            .float => if (!std.math.isFinite(@field(result, field.name)))
                return error.NonFiniteWaterEquilibriumSolution,
            .bool => {},
            else => unreachable,
        }
    }
    if (result.hydrogen_activity_mol_per_m3 <= 0 or
        result.hydroxide_activity_mol_per_m3 <= 0 or
        result.hydrogen_concentration_mol_per_m3 <= 0 or
        result.hydroxide_concentration_mol_per_m3 <= 0)
    {
        return error.InvalidWaterEquilibriumSolution;
    }
    return result;
}

/// Projects a possibly signed, provisional H+/OH- pair onto the water
/// activity product while preserving its activity difference. Chemical
/// transformations are allowed to cross either free-ion inventory because
/// H+ and OH- are dependent coordinates joined by H2O dissociation; only the
/// projected thermodynamic state must be positive.
pub fn projectProvisional(
    provisional_hydrogen_concentration_mol_per_m3: f64,
    provisional_hydroxide_concentration_mol_per_m3: f64,
    monovalent_activity_coefficient: f64,
    water_activity_product_mol2_per_m6: f64,
) !Result {
    if (!std.math.isFinite(provisional_hydrogen_concentration_mol_per_m3) or
        !std.math.isFinite(provisional_hydroxide_concentration_mol_per_m3) or
        !std.math.isFinite(monovalent_activity_coefficient) or
        !std.math.isFinite(water_activity_product_mol2_per_m6))
        return error.NonFiniteWaterEquilibriumInput;
    if (monovalent_activity_coefficient <= 0 or
        water_activity_product_mol2_per_m6 <= 0)
        return error.InvalidWaterEquilibriumInput;

    const activity_difference = monovalent_activity_coefficient *
        (provisional_hydrogen_concentration_mol_per_m3 -
            provisional_hydroxide_concentration_mol_per_m3);
    const discriminant = std.math.hypot(
        activity_difference,
        2 * @sqrt(water_activity_product_mol2_per_m6),
    );
    const hydrogen_activity, const hydroxide_activity =
        if (activity_difference >= 0) .{
            0.5 * (activity_difference + discriminant),
            water_activity_product_mol2_per_m6 /
                (0.5 * (activity_difference + discriminant)),
        } else .{
            water_activity_product_mol2_per_m6 /
                (0.5 * (-activity_difference + discriminant)),
            0.5 * (-activity_difference + discriminant),
        };
    if (hydrogen_activity <= 0 or hydroxide_activity <= 0)
        return error.InvalidWaterEquilibriumSolution;
    return .{
        .hydrogen_concentration_mol_per_m3 = hydrogen_activity / monovalent_activity_coefficient,
        .hydroxide_concentration_mol_per_m3 = hydroxide_activity / monovalent_activity_coefficient,
        .hydrogen_activity_mol_per_m3 = hydrogen_activity,
        .hydroxide_activity_mol_per_m3 = hydroxide_activity,
        .equal_reaction_extent_mol_per_m3 = provisional_hydrogen_concentration_mol_per_m3 -
            hydrogen_activity / monovalent_activity_coefficient,
        .ph = -@log10(hydrogen_activity * 1e-3),
    };
}

/// Exact equation-order diagnostic for the final reset at SOLUTE.F
/// lines 2681-2725.
///
/// The source computes the equal extent from activities, then subtracts that
/// activity-space value directly from concentration-space transformation
/// totals. The returned published concentrations expose this dimensional
/// mismatch for attribution; production uses `solve`/`projectProvisional`.
pub fn calculateSourceOrderFinalReset(inputs: FinalResetInputs) !SourceOrderFinalReset {
    inline for (@typeInfo(FinalResetInputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteWaterEquilibriumInput;
    }
    if (inputs.initial_hydrogen_concentration_mol_per_m3 < 0 or
        inputs.initial_hydroxide_concentration_mol_per_m3 < 0 or
        inputs.monovalent_activity_coefficient <= 0 or
        inputs.water_activity_product_mol2_per_m6 <= 0 or
        inputs.negligible_concentration_mol_per_m3 < 0)
        return error.InvalidWaterEquilibriumInput;

    const provisional_hydrogen =
        inputs.initial_hydrogen_concentration_mol_per_m3 +
        inputs.accumulated_hydrogen_change_mol_per_m3;
    const provisional_hydroxide =
        inputs.initial_hydroxide_concentration_mol_per_m3 +
        inputs.accumulated_hydroxide_change_mol_per_m3;
    const hydrogen_activity = if (provisional_hydrogen >
        inputs.negligible_concentration_mol_per_m3)
        provisional_hydrogen * inputs.monovalent_activity_coefficient
    else
        provisional_hydrogen;
    const hydroxide_activity = if (provisional_hydroxide >
        inputs.negligible_concentration_mol_per_m3)
        provisional_hydroxide * inputs.monovalent_activity_coefficient
    else
        provisional_hydroxide;
    const activity_sum = hydrogen_activity + hydroxide_activity;
    const discriminant = @max(
        0,
        activity_sum * activity_sum -
            4 * (hydrogen_activity * hydroxide_activity -
                inputs.water_activity_product_mol2_per_m6),
    );
    const extent = 0.5 * (activity_sum - @sqrt(discriminant));
    const hydrogen_change =
        inputs.accumulated_hydrogen_change_mol_per_m3 - extent;
    const hydroxide_change =
        inputs.accumulated_hydroxide_change_mol_per_m3 - extent;
    const result = SourceOrderFinalReset{
        .hydrogen_activity_mol_per_m3 = hydrogen_activity - extent,
        .hydroxide_activity_mol_per_m3 = hydroxide_activity - extent,
        .equal_reaction_extent_activity_mol_per_m3 = extent,
        .accumulated_hydrogen_change_mol_per_m3 = hydrogen_change,
        .accumulated_hydroxide_change_mol_per_m3 = hydroxide_change,
        .published_hydrogen_concentration_mol_per_m3 = inputs.initial_hydrogen_concentration_mol_per_m3 +
            hydrogen_change,
        .published_hydroxide_concentration_mol_per_m3 = inputs.initial_hydroxide_concentration_mol_per_m3 +
            hydroxide_change,
    };
    inline for (@typeInfo(SourceOrderFinalReset).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteWaterEquilibriumSolution;
    }
    return result;
}

test "water equilibrium reaches activity product with equal transformation" {
    const inputs: Inputs = .{ .hydrogen_concentration_mol_per_m3 = 1e-3, .hydroxide_concentration_mol_per_m3 = 2e-3, .monovalent_activity_coefficient = 0.8, .water_activity_product_mol2_per_m6 = 1e-8, .negligible_concentration_mol_per_m3 = 1e-32 };
    const result = try solve(inputs);
    try std.testing.expectApproxEqRel(inputs.water_activity_product_mol2_per_m6, result.hydrogen_activity_mol_per_m3 * result.hydroxide_activity_mol_per_m3, 1e-13);
    try std.testing.expectApproxEqAbs(inputs.hydrogen_concentration_mol_per_m3 * inputs.monovalent_activity_coefficient - result.hydrogen_activity_mol_per_m3, inputs.hydroxide_concentration_mol_per_m3 * inputs.monovalent_activity_coefficient - result.hydroxide_activity_mol_per_m3, 1e-15);
    try std.testing.expect(std.math.isFinite(result.ph));
}

test "water equilibrium matches source quadratic away from cancellation" {
    const inputs: Inputs = .{
        .hydrogen_concentration_mol_per_m3 = 0.4,
        .hydroxide_concentration_mol_per_m3 = 0.1,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .negligible_concentration_mol_per_m3 = 1.0e-32,
    };
    const source_hydrogen_activity =
        inputs.hydrogen_concentration_mol_per_m3 *
        inputs.monovalent_activity_coefficient;
    const source_hydroxide_activity =
        inputs.hydroxide_concentration_mol_per_m3 *
        inputs.monovalent_activity_coefficient;
    const source_sum =
        source_hydrogen_activity + source_hydroxide_activity;
    const source_discriminant = @max(
        0,
        source_sum * source_sum -
            4 * (source_hydrogen_activity *
                source_hydroxide_activity -
                inputs.water_activity_product_mol2_per_m6),
    );
    const source_extent =
        0.5 * (source_sum - @sqrt(source_discriminant));
    const result = try solve(inputs);
    try std.testing.expectApproxEqAbs(
        source_hydrogen_activity - source_extent,
        result.hydrogen_activity_mol_per_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        source_hydroxide_activity - source_extent,
        result.hydroxide_activity_mol_per_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        source_extent,
        result.equal_reaction_extent_mol_per_m3,
        1.0e-15,
    );
}

test "SOLUTE surface pH reset preserves every source expression" {
    const inputs: Inputs = .{
        .hydrogen_concentration_mol_per_m3 = 0.4,
        .hydroxide_concentration_mol_per_m3 = 0.1,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .negligible_concentration_mol_per_m3 = 1.0e-32,
    };
    const result = try calculateSourceOrderSurfaceReset(inputs);
    const expected_hydrogen_activity =
        inputs.hydrogen_concentration_mol_per_m3 *
        inputs.monovalent_activity_coefficient;
    const expected_hydroxide_activity =
        inputs.hydroxide_concentration_mol_per_m3 *
        inputs.monovalent_activity_coefficient;
    const expected_sum =
        expected_hydrogen_activity + expected_hydroxide_activity;
    const expected_discriminant = @max(
        0.0,
        expected_sum * expected_sum -
            4.0 *
                (expected_hydrogen_activity *
                    expected_hydroxide_activity -
                    inputs.water_activity_product_mol2_per_m6),
    );
    const expected_extent =
        0.5 * (expected_sum - @sqrt(expected_discriminant));
    const expected_final_hydrogen = @max(
        inputs.negligible_concentration_mol_per_m3,
        expected_hydrogen_activity - expected_extent,
    );
    const expected_final_hydroxide = @max(
        inputs.negligible_concentration_mol_per_m3,
        expected_hydroxide_activity - expected_extent,
    );

    try std.testing.expectEqual(
        expected_hydrogen_activity,
        result.initial_hydrogen_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_hydroxide_activity,
        result.initial_hydroxide_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_sum,
        result.activity_sum_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_discriminant,
        result.nonnegative_discriminant_mol2_per_m6,
    );
    try std.testing.expectEqual(
        expected_extent,
        result.equal_reaction_extent_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_final_hydrogen,
        result.hydrogen_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_final_hydroxide,
        result.hydroxide_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_final_hydrogen /
            inputs.monovalent_activity_coefficient,
        result.hydrogen_concentration_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_final_hydroxide /
            inputs.monovalent_activity_coefficient,
        result.hydroxide_concentration_mol_per_m3,
    );
    try std.testing.expectEqual(
        -@log10(expected_final_hydrogen * 1.0e-3),
        result.ph,
    );
}

test "SOLUTE surface pH activity gates are strict" {
    const inputs: Inputs = .{
        .hydrogen_concentration_mol_per_m3 = 1.0e-4,
        .hydroxide_concentration_mol_per_m3 = 1.0e-4,
        .monovalent_activity_coefficient = 0.5,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .negligible_concentration_mol_per_m3 = 1.0e-4,
    };
    const result = try calculateSourceOrderSurfaceReset(inputs);
    try std.testing.expectEqual(
        inputs.hydrogen_concentration_mol_per_m3,
        result.initial_hydrogen_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        inputs.hydroxide_concentration_mol_per_m3,
        result.initial_hydroxide_activity_mol_per_m3,
    );
}

test "production water equilibrium preserves trace ion lost by source arithmetic" {
    const inputs: Inputs = .{
        .hydrogen_concentration_mol_per_m3 = 1.0e12,
        .hydroxide_concentration_mol_per_m3 = 1.0e-20,
        .monovalent_activity_coefficient = 0.7,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .negligible_concentration_mol_per_m3 = 1.0e-32,
    };
    const source = try calculateSourceOrderSurfaceReset(inputs);
    const production = try solve(inputs);
    const source_product =
        source.hydrogen_activity_mol_per_m3 *
        source.hydroxide_activity_mol_per_m3;
    const production_product =
        production.hydrogen_activity_mol_per_m3 *
        production.hydroxide_activity_mol_per_m3;

    try std.testing.expect(
        @abs(source_product -
            inputs.water_activity_product_mol2_per_m6) >
            0.1 * inputs.water_activity_product_mol2_per_m6,
    );
    try std.testing.expectApproxEqRel(
        inputs.water_activity_product_mol2_per_m6,
        production_product,
        4 * std.math.floatEps(f64),
    );
}

test "SOLUTE surface pH reset rejects invalid input and overflow" {
    var inputs: Inputs = .{
        .hydrogen_concentration_mol_per_m3 = 0.4,
        .hydroxide_concentration_mol_per_m3 = 0.1,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .negligible_concentration_mol_per_m3 = 1.0e-32,
    };
    inputs.negligible_concentration_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidWaterEquilibriumInput,
        calculateSourceOrderSurfaceReset(inputs),
    );

    inputs.negligible_concentration_mol_per_m3 = 1.0e-32;
    inputs.hydrogen_concentration_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteWaterEquilibriumInput,
        calculateSourceOrderSurfaceReset(inputs),
    );

    inputs.hydrogen_concentration_mol_per_m3 =
        std.math.floatMax(f64);
    inputs.monovalent_activity_coefficient = 2;
    try std.testing.expectError(
        error.NonFiniteWaterEquilibriumSolution,
        calculateSourceOrderSurfaceReset(inputs),
    );
}

test "cancellation-safe water equilibrium handles nearly balanced activities" {
    const result = try solve(.{ .hydrogen_concentration_mol_per_m3 = 1.000001e-4, .hydroxide_concentration_mol_per_m3 = 1e-4, .monovalent_activity_coefficient = 1, .water_activity_product_mol2_per_m6 = 1e-8, .negligible_concentration_mol_per_m3 = 1e-32 });
    try std.testing.expectApproxEqRel(@as(f64, 1e-8), result.hydrogen_activity_mol_per_m3 * result.hydroxide_activity_mol_per_m3, 1e-14);
}

test "water equilibrium preserves trace hydroxide under extreme acidity" {
    const result = try solve(.{
        .hydrogen_concentration_mol_per_m3 = 1e12,
        .hydroxide_concentration_mol_per_m3 = 1e-20,
        .monovalent_activity_coefficient = 0.7,
        .water_activity_product_mol2_per_m6 = 1e-8,
        .negligible_concentration_mol_per_m3 = 1e-32,
    });
    try std.testing.expect(result.hydroxide_activity_mol_per_m3 > 0);
    try std.testing.expectApproxEqRel(
        @as(f64, 1e-8),
        result.hydrogen_activity_mol_per_m3 *
            result.hydroxide_activity_mol_per_m3,
        4 * std.math.floatEps(f64),
    );
}

test "water equilibrium projects a provisional negative free ion" {
    const result = try projectProvisional(-0.2, 0.1, 0.8, 1e-8);
    try std.testing.expect(result.hydrogen_concentration_mol_per_m3 > 0);
    try std.testing.expect(result.hydroxide_concentration_mol_per_m3 > 0);
    try std.testing.expectApproxEqRel(
        @as(f64, 1e-8),
        result.hydrogen_activity_mol_per_m3 *
            result.hydroxide_activity_mol_per_m3,
        8 * std.math.floatEps(f64),
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.3),
        result.hydrogen_concentration_mol_per_m3 -
            result.hydroxide_concentration_mol_per_m3,
        1e-15,
    );
}

test "final source reset mixes activity extent into concentration totals" {
    const source = try calculateSourceOrderFinalReset(.{
        .initial_hydrogen_concentration_mol_per_m3 = 0.4,
        .initial_hydroxide_concentration_mol_per_m3 = 0.1,
        .accumulated_hydrogen_change_mol_per_m3 = 0,
        .accumulated_hydroxide_change_mol_per_m3 = 0,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1e-8,
        .negligible_concentration_mol_per_m3 = 1e-32,
    });
    try std.testing.expectApproxEqRel(
        @as(f64, 1e-8),
        source.hydrogen_activity_mol_per_m3 *
            source.hydroxide_activity_mol_per_m3,
        2e-10,
    );
    const source_published_product =
        source.published_hydrogen_concentration_mol_per_m3 *
        source.published_hydroxide_concentration_mol_per_m3 *
        0.8 * 0.8;
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.0799999583333406),
        source.equal_reaction_extent_activity_mol_per_m3,
        1e-15,
    );
    try std.testing.expect(source_published_product > 0.004);

    const corrected = try solve(.{
        .hydrogen_concentration_mol_per_m3 = 0.4,
        .hydroxide_concentration_mol_per_m3 = 0.1,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1e-8,
        .negligible_concentration_mol_per_m3 = 1e-32,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.019999989583342537),
        corrected.hydrogen_concentration_mol_per_m3 -
            source.published_hydrogen_concentration_mol_per_m3,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.01999998958334254),
        corrected.hydroxide_concentration_mol_per_m3 -
            source.published_hydroxide_concentration_mol_per_m3,
        1e-14,
    );
}
