const std = @import("std");

pub const ResourceValues = struct {
    oxygen: f64,
    ammonium_non_band: f64,
    ammonium_band: f64,
    nitrate_non_band: f64,
    nitrate_band: f64,
    phosphate_2_non_band: f64,
    phosphate_2_band: f64,
    phosphate_1_non_band: f64,
    phosphate_1_band: f64,
};

pub const Inputs = struct {
    total_competing_demand: ResourceValues,
    root_domain_demand: ResourceValues,
    minimum_competition_fraction: f64,
    fallback_root_mass_fraction: f64,
    significant_total_demand: f64,
    preceding_fraction_totals: ResourceValues,
};

pub const Result = struct {
    fractions: ResourceValues,
    fraction_totals: ResourceValues,
};

/// UPTAKE.F 1810--1863. Calculates one root/mycorrhizal domain's fractions of
/// nine competing demands, then accumulates them in exact source order.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    var fractions: ResourceValues = undefined;
    fractions.oxygen = fraction(
        inputs.total_competing_demand.oxygen,
        inputs.root_domain_demand.oxygen,
        inputs,
    );
    fractions.ammonium_non_band = fraction(
        inputs.total_competing_demand.ammonium_non_band,
        inputs.root_domain_demand.ammonium_non_band,
        inputs,
    );
    fractions.ammonium_band = fraction(
        inputs.total_competing_demand.ammonium_band,
        inputs.root_domain_demand.ammonium_band,
        inputs,
    );
    fractions.nitrate_non_band = fraction(
        inputs.total_competing_demand.nitrate_non_band,
        inputs.root_domain_demand.nitrate_non_band,
        inputs,
    );
    fractions.nitrate_band = fraction(
        inputs.total_competing_demand.nitrate_band,
        inputs.root_domain_demand.nitrate_band,
        inputs,
    );
    fractions.phosphate_2_non_band = fraction(
        inputs.total_competing_demand.phosphate_2_non_band,
        inputs.root_domain_demand.phosphate_2_non_band,
        inputs,
    );
    fractions.phosphate_2_band = fraction(
        inputs.total_competing_demand.phosphate_2_band,
        inputs.root_domain_demand.phosphate_2_band,
        inputs,
    );
    fractions.phosphate_1_non_band = fraction(
        inputs.total_competing_demand.phosphate_1_non_band,
        inputs.root_domain_demand.phosphate_1_non_band,
        inputs,
    );
    fractions.phosphate_1_band = fraction(
        inputs.total_competing_demand.phosphate_1_band,
        inputs.root_domain_demand.phosphate_1_band,
        inputs,
    );

    var totals = inputs.preceding_fraction_totals;
    totals.oxygen += fractions.oxygen;
    totals.ammonium_non_band += fractions.ammonium_non_band;
    totals.nitrate_non_band += fractions.nitrate_non_band;
    totals.phosphate_2_non_band += fractions.phosphate_2_non_band;
    totals.phosphate_1_non_band += fractions.phosphate_1_non_band;
    totals.ammonium_band += fractions.ammonium_band;
    totals.nitrate_band += fractions.nitrate_band;
    totals.phosphate_2_band += fractions.phosphate_2_band;
    totals.phosphate_1_band += fractions.phosphate_1_band;
    inline for (@typeInfo(ResourceValues).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(fractions, field.name)) or
            !std.math.isFinite(@field(totals, field.name)))
            return error.NonFiniteRootNutrientCompetitionResult;
    }
    return .{ .fractions = fractions, .fraction_totals = totals };
}

fn fraction(total: f64, individual: f64, inputs: Inputs) f64 {
    return if (total > inputs.significant_total_demand)
        @max(inputs.minimum_competition_fraction, individual / total)
    else
        inputs.fallback_root_mass_fraction;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ResourceValues).@"struct".fields) |field| {
        const total = @field(inputs.total_competing_demand, field.name);
        const individual = @field(inputs.root_domain_demand, field.name);
        const preceding = @field(inputs.preceding_fraction_totals, field.name);
        if (!std.math.isFinite(total) or !std.math.isFinite(individual) or
            !std.math.isFinite(preceding) or total < 0 or individual < 0)
            return error.InvalidRootNutrientCompetitionInput;
    }
    inline for (.{
        inputs.minimum_competition_fraction,
        inputs.fallback_root_mass_fraction,
        inputs.significant_total_demand,
    }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootNutrientCompetitionInput;
}

fn values(base: f64) ResourceValues {
    return .{
        .oxygen = base + 1,
        .ammonium_non_band = base + 2,
        .ammonium_band = base + 3,
        .nitrate_non_band = base + 4,
        .nitrate_band = base + 5,
        .phosphate_2_non_band = base + 6,
        .phosphate_2_band = base + 7,
        .phosphate_1_non_band = base + 8,
        .phosphate_1_band = base + 9,
    };
}

test "UPTAKE competition fractions preserve demand and fallback branches" {
    var totals = values(9);
    totals.ammonium_band = 0;
    totals.phosphate_1_band = 1e-14;
    const individuals = values(0);
    const result = try calculate(.{
        .total_competing_demand = totals,
        .root_domain_demand = individuals,
        .minimum_competition_fraction = 0.2,
        .fallback_root_mass_fraction = 0.3,
        .significant_total_demand = 1e-12,
        .preceding_fraction_totals = values(0),
    });
    try std.testing.expectEqual(@as(f64, 0.2), result.fractions.oxygen);
    try std.testing.expectEqual(@as(f64, 0.3), result.fractions.ammonium_band);
    try std.testing.expectEqual(@as(f64, 0.3), result.fractions.phosphate_1_band);
    try std.testing.expectEqual(
        @as(f64, 1.2),
        result.fraction_totals.oxygen,
    );
}

test "source maximum keeps individual over total fraction above one" {
    var totals = values(9);
    var individuals = values(0);
    totals.oxygen = 1;
    individuals.oxygen = 2;
    const result = try calculate(.{
        .total_competing_demand = totals,
        .root_domain_demand = individuals,
        .minimum_competition_fraction = 0.2,
        .fallback_root_mass_fraction = 0.3,
        .significant_total_demand = 1e-12,
        .preceding_fraction_totals = values(0),
    });
    try std.testing.expectEqual(@as(f64, 2), result.fractions.oxygen);
}

test "non-finite late resource fails explicitly" {
    var individuals = values(0);
    individuals.phosphate_1_band = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidRootNutrientCompetitionInput,
        calculate(.{
            .total_competing_demand = values(9),
            .root_domain_demand = individuals,
            .minimum_competition_fraction = 0.2,
            .fallback_root_mass_fraction = 0.3,
            .significant_total_demand = 1e-12,
            .preceding_fraction_totals = values(0),
        }),
    );
}
