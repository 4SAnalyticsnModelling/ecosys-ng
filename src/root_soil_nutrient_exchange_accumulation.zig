const std = @import("std");
const nutrient = @import("root_nutrient_uptake_reset.zig");
const organic = @import("root_organic_exchange.zig");

pub const ElementTotals = struct {
    carbon_g_c_per_step: f64,
    nitrogen_g_n_per_step: f64,
    phosphorus_g_p_per_step: f64,
};

pub const InorganicUptakeTotals = struct {
    ammonium_g_n_per_step: f64,
    nitrate_g_n_per_step: f64,
    dihydrogen_phosphate_g_p_per_step: f64,
    hydrogen_phosphate_g_p_per_step: f64,
};

pub const Totals = struct {
    plant_organic_exchange: ElementTotals,
    soil_dissolved_pool_change: organic.Result,
    inorganic_uptake: InorganicUptakeTotals,
};

pub const Contribution = struct {
    organic_exchange: organic.Result,
    nutrient_uptake: nutrient.State,
};

/// UPTAKE.F 3820--3842. Each cell owns a contiguous, source-ordered range of
/// runtime species/root/layer contributions described by `contribution_offsets`.
pub fn accumulateRuntimeCells(
    initial: []const Totals,
    contribution_offsets: []const usize,
    contributions: []const Contribution,
    scratch: []Totals,
    destination: []Totals,
) !void {
    if (initial.len != scratch.len or initial.len != destination.len or
        contribution_offsets.len != initial.len + 1)
        return error.RootSoilNutrientAccumulationDimensionMismatch;
    if (contribution_offsets[0] != 0 or
        contribution_offsets[contribution_offsets.len - 1] != contributions.len)
        return error.InvalidRootSoilNutrientContributionOffsets;

    for (initial, scratch, 0..) |before, *candidate, cell| {
        const begin = contribution_offsets[cell];
        const end = contribution_offsets[cell + 1];
        if (begin > end or end > contributions.len)
            return error.InvalidRootSoilNutrientContributionOffsets;
        candidate.* = try accumulateContributions(before, contributions[begin..end]);
    }
    @memcpy(destination, scratch);
}

pub fn accumulateContributions(
    before: Totals,
    contributions: []const Contribution,
) !Totals {
    try validateFinite(Totals, before);
    var result = before;
    for (contributions) |contribution| {
        try validateFinite(Contribution, contribution);
        accumulateOne(&result, contribution);
        try validateFinite(Totals, result);
    }
    return result;
}

fn accumulateOne(result: *Totals, contribution: Contribution) void {
    // UPTAKE.F 3820--3834: K-major, C/N/P plant totals, then C/N/P soil
    // changes for each of the five scientific substrate classes.
    for (0..organic.substrate_count) |substrate| {
        const carbon = contribution.organic_exchange.carbon_exchange_g_c_per_step[substrate];
        const nitrogen = contribution.organic_exchange.nitrogen_exchange_g_n_per_step[substrate];
        const phosphorus = contribution.organic_exchange.phosphorus_exchange_g_p_per_step[substrate];
        result.plant_organic_exchange.carbon_g_c_per_step += carbon;
        result.plant_organic_exchange.nitrogen_g_n_per_step += nitrogen;
        result.plant_organic_exchange.phosphorus_g_p_per_step += phosphorus;
        result.soil_dissolved_pool_change.carbon_exchange_g_c_per_step[substrate] -= carbon;
        result.soil_dissolved_pool_change.nitrogen_exchange_g_n_per_step[substrate] -= nitrogen;
        result.soil_dissolved_pool_change.phosphorus_exchange_g_p_per_step[substrate] -= phosphorus;
    }

    const uptake = contribution.nutrient_uptake;
    result.inorganic_uptake.ammonium_g_n_per_step +=
        uptake.nitrogen.ammonium_non_band.population_uptake_g_n_per_step;
    result.inorganic_uptake.ammonium_g_n_per_step +=
        uptake.nitrogen.ammonium_band.population_uptake_g_n_per_step;
    result.inorganic_uptake.nitrate_g_n_per_step +=
        uptake.nitrogen.nitrate_non_band.population_uptake_g_n_per_step;
    result.inorganic_uptake.nitrate_g_n_per_step +=
        uptake.nitrogen.nitrate_band.population_uptake_g_n_per_step;
    result.inorganic_uptake.dihydrogen_phosphate_g_p_per_step +=
        uptake.phosphorus.dihydrogen_phosphate_non_band.population_uptake_g_p_per_step;
    result.inorganic_uptake.dihydrogen_phosphate_g_p_per_step +=
        uptake.phosphorus.dihydrogen_phosphate_band.population_uptake_g_p_per_step;
    result.inorganic_uptake.hydrogen_phosphate_g_p_per_step +=
        uptake.phosphorus.hydrogen_phosphate_non_band.population_uptake_g_p_per_step;
    result.inorganic_uptake.hydrogen_phosphate_g_p_per_step +=
        uptake.phosphorus.hydrogen_phosphate_band.population_uptake_g_p_per_step;
}

fn validateFinite(comptime T: type, value: T) !void {
    switch (@typeInfo(T)) {
        .float => if (!std.math.isFinite(value))
            return error.NonFiniteRootSoilNutrientAccumulation,
        .array => |array| for (value) |element|
            try validateFinite(array.child, element),
        .@"struct" => |structure| inline for (structure.fields) |field|
            try validateFinite(field.type, @field(value, field.name)),
        else => {},
    }
}

fn testContribution(scale: f64) Contribution {
    var value = std.mem.zeroes(Contribution);
    for (0..organic.substrate_count) |substrate| {
        const index: f64 = @floatFromInt(substrate + 1);
        value.organic_exchange.carbon_exchange_g_c_per_step[substrate] = scale * index;
        value.organic_exchange.nitrogen_exchange_g_n_per_step[substrate] = scale * index * 0.1;
        value.organic_exchange.phosphorus_exchange_g_p_per_step[substrate] = scale * index * 0.01;
    }
    value.nutrient_uptake.nitrogen.ammonium_non_band.population_uptake_g_n_per_step = scale;
    value.nutrient_uptake.nitrogen.ammonium_band.population_uptake_g_n_per_step = scale * 2;
    value.nutrient_uptake.nitrogen.nitrate_non_band.population_uptake_g_n_per_step = scale * 3;
    value.nutrient_uptake.nitrogen.nitrate_band.population_uptake_g_n_per_step = scale * 4;
    value.nutrient_uptake.phosphorus.dihydrogen_phosphate_non_band.population_uptake_g_p_per_step = scale * 5;
    value.nutrient_uptake.phosphorus.dihydrogen_phosphate_band.population_uptake_g_p_per_step = scale * 6;
    value.nutrient_uptake.phosphorus.hydrogen_phosphate_non_band.population_uptake_g_p_per_step = scale * 7;
    value.nutrient_uptake.phosphorus.hydrogen_phosphate_band.population_uptake_g_p_per_step = scale * 8;
    return value;
}

test "source-ordered root-soil exchange accumulates all element groups" {
    const result = try accumulateContributions(
        std.mem.zeroes(Totals),
        &[_]Contribution{testContribution(2)},
    );
    try std.testing.expectApproxEqAbs(@as(f64, 30), result.plant_organic_exchange.carbon_g_c_per_step, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), result.plant_organic_exchange.nitrogen_g_n_per_step, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.plant_organic_exchange.phosphorus_g_p_per_step, 1e-12);
    try std.testing.expectEqual(@as(f64, 6), result.inorganic_uptake.ammonium_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 14), result.inorganic_uptake.nitrate_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 22), result.inorganic_uptake.dihydrogen_phosphate_g_p_per_step);
    try std.testing.expectEqual(@as(f64, 30), result.inorganic_uptake.hydrogen_phosphate_g_p_per_step);
}

test "organic exchange conserves carbon nitrogen and phosphorus exactly" {
    const contribution = testContribution(-3);
    const result = try accumulateContributions(
        std.mem.zeroes(Totals),
        &[_]Contribution{contribution},
    );
    var soil_carbon: f64 = 0;
    var soil_nitrogen: f64 = 0;
    var soil_phosphorus: f64 = 0;
    for (0..organic.substrate_count) |substrate| {
        soil_carbon += result.soil_dissolved_pool_change.carbon_exchange_g_c_per_step[substrate];
        soil_nitrogen += result.soil_dissolved_pool_change.nitrogen_exchange_g_n_per_step[substrate];
        soil_phosphorus += result.soil_dissolved_pool_change.phosphorus_exchange_g_p_per_step[substrate];
    }
    try std.testing.expectEqual(@as(f64, 0), result.plant_organic_exchange.carbon_g_c_per_step + soil_carbon);
    try std.testing.expectEqual(@as(f64, 0), result.plant_organic_exchange.nitrogen_g_n_per_step + soil_nitrogen);
    try std.testing.expectEqual(@as(f64, 0), result.plant_organic_exchange.phosphorus_g_p_per_step + soil_phosphorus);
}

test "runtime cells use independent source-ordered contribution ranges" {
    const initial = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    const offsets = [_]usize{ 0, 2, 3 };
    const contributions = [_]Contribution{
        testContribution(1),
        testContribution(2),
        testContribution(10),
    };
    var scratch: [2]Totals = undefined;
    var destination: [2]Totals = undefined;
    try accumulateRuntimeCells(&initial, &offsets, &contributions, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 9), destination[0].inorganic_uptake.ammonium_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 30), destination[1].inorganic_uptake.ammonium_g_n_per_step);
}

test "later nonfinite contribution leaves every destination cell unchanged" {
    const initial = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    const offsets = [_]usize{ 0, 1, 2 };
    var contributions = [_]Contribution{ testContribution(1), testContribution(2) };
    contributions[1].organic_exchange.phosphorus_exchange_g_p_per_step[4] =
        std.math.nan(f64);
    var scratch: [2]Totals = undefined;
    var destination = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    destination[0].inorganic_uptake.ammonium_g_n_per_step = 41;
    destination[1].inorganic_uptake.ammonium_g_n_per_step = 42;
    try std.testing.expectError(
        error.NonFiniteRootSoilNutrientAccumulation,
        accumulateRuntimeCells(&initial, &offsets, &contributions, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].inorganic_uptake.ammonium_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 42), destination[1].inorganic_uptake.ammonium_g_n_per_step);
}
