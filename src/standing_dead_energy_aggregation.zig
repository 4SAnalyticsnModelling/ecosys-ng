const std = @import("std");

pub const Totals = struct {
    net_radiation_mj_per_step: f64,
    latent_heat_mj_per_step: f64,
    sensible_heat_mj_per_step: f64,
    storage_heat_mj_per_step: f64,
    vapor_convective_heat_mj_per_step: f64,
    ground_longwave_mj_per_step: f64,
    evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    ground_vapor_flux_m3_per_step: f64,
    ground_sensible_heat_mj_per_step: f64,
    wet_heat_capacity_mj_per_k: f64,
    energy_balance_residual_mj_per_step: f64,
};

pub const Contribution = struct {
    net_radiation_mj_per_step: f64,
    latent_heat_mj_per_step: f64,
    sensible_heat_mj_per_step: f64,
    storage_heat_mj_per_step: f64,
    vapor_convective_heat_mj_per_step: f64,
    ground_longwave_mj_per_step: f64,
    evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    ground_vapor_flux_m3_per_step: f64,
    ground_sensible_heat_mj_per_step: f64,
    precipitation_convective_heat_mj_per_step: f64,
    intercepted_water_m3: f64,
};

pub const SpeciesInputs = struct {
    dry_heat_capacity_mj_per_k: f64,
};

/// UPTAKE.F 4240--4250. Each runtime species owns a contiguous source-ordered
/// range of local energy contributions.
pub fn accumulateRuntimeSpecies(
    initial: []const Totals,
    species_inputs: []const SpeciesInputs,
    contribution_offsets: []const usize,
    contributions: []const Contribution,
    scratch: []Totals,
    destination: []Totals,
) !void {
    if (initial.len != species_inputs.len or initial.len != scratch.len or
        initial.len != destination.len or
        contribution_offsets.len != initial.len + 1)
        return error.StandingDeadEnergyAggregationDimensionMismatch;
    if (contribution_offsets[0] != 0 or
        contribution_offsets[contribution_offsets.len - 1] != contributions.len)
        return error.InvalidStandingDeadEnergyContributionOffsets;
    for (initial, species_inputs, scratch, 0..) |before, inputs, *candidate, species| {
        const begin = contribution_offsets[species];
        const end = contribution_offsets[species + 1];
        if (begin > end or end > contributions.len)
            return error.InvalidStandingDeadEnergyContributionOffsets;
        candidate.* = try accumulate(
            before,
            inputs,
            contributions[begin..end],
        );
    }
    @memcpy(destination, scratch);
}

pub fn accumulate(
    before: Totals,
    inputs: SpeciesInputs,
    contributions: []const Contribution,
) !Totals {
    try validateFinite(Totals, before);
    if (!std.math.isFinite(inputs.dry_heat_capacity_mj_per_k) or
        inputs.dry_heat_capacity_mj_per_k < 0)
        return error.InvalidStandingDeadEnergyAggregationInput;
    var result = before;
    for (contributions) |contribution| {
        try validateContribution(contribution);
        result.net_radiation_mj_per_step += contribution.net_radiation_mj_per_step;
        result.latent_heat_mj_per_step += contribution.latent_heat_mj_per_step;
        result.sensible_heat_mj_per_step += contribution.sensible_heat_mj_per_step;
        result.storage_heat_mj_per_step += contribution.storage_heat_mj_per_step;
        result.vapor_convective_heat_mj_per_step +=
            contribution.vapor_convective_heat_mj_per_step;
        result.ground_longwave_mj_per_step +=
            contribution.ground_longwave_mj_per_step;
        result.evaporation_m3_per_step += contribution.evaporation_m3_per_step;
        result.transpiration_m3_per_step += contribution.transpiration_m3_per_step;
        result.ground_vapor_flux_m3_per_step +=
            contribution.ground_vapor_flux_m3_per_step;
        result.ground_sensible_heat_mj_per_step +=
            contribution.ground_sensible_heat_mj_per_step;
        result.wet_heat_capacity_mj_per_k =
            inputs.dry_heat_capacity_mj_per_k +
            4.19 * @max(0, contribution.intercepted_water_m3);
        result.energy_balance_residual_mj_per_step +=
            contribution.storage_heat_mj_per_step -
            (contribution.net_radiation_mj_per_step +
                contribution.latent_heat_mj_per_step +
                contribution.sensible_heat_mj_per_step +
                contribution.vapor_convective_heat_mj_per_step +
                contribution.precipitation_convective_heat_mj_per_step);
        try validateFinite(Totals, result);
    }
    return result;
}

fn validateContribution(contribution: Contribution) !void {
    try validateFinite(Contribution, contribution);
    if (contribution.intercepted_water_m3 < 0)
        return error.InvalidStandingDeadEnergyAggregationInput;
}

fn validateFinite(comptime T: type, value: T) !void {
    inline for (@typeInfo(T).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteStandingDeadEnergyAggregation;
}

fn testContribution(scale: f64, water_m3: f64) Contribution {
    return .{
        .net_radiation_mj_per_step = 5 * scale,
        .latent_heat_mj_per_step = -2 * scale,
        .sensible_heat_mj_per_step = -1 * scale,
        .storage_heat_mj_per_step = 1.5 * scale,
        .vapor_convective_heat_mj_per_step = -0.5 * scale,
        .ground_longwave_mj_per_step = 0.3 * scale,
        .evaporation_m3_per_step = -0.01 * scale,
        .transpiration_m3_per_step = 0,
        .ground_vapor_flux_m3_per_step = 0.02 * scale,
        .ground_sensible_heat_mj_per_step = 0.4 * scale,
        .precipitation_convective_heat_mj_per_step = 0,
        .intercepted_water_m3 = water_m3,
    };
}

test "standing-dead aggregation preserves ten source additions" {
    const contributions = [_]Contribution{
        testContribution(1, 0.2),
        testContribution(2, 0.1),
    };
    const result = try accumulate(
        std.mem.zeroes(Totals),
        .{ .dry_heat_capacity_mj_per_k = 2 },
        &contributions,
    );
    try std.testing.expectEqual(@as(f64, 15), result.net_radiation_mj_per_step);
    try std.testing.expectEqual(@as(f64, -6), result.latent_heat_mj_per_step);
    try std.testing.expectEqual(@as(f64, 4.5), result.storage_heat_mj_per_step);
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), result.ground_longwave_mj_per_step, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.03), result.evaporation_m3_per_step, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), result.ground_sensible_heat_mj_per_step, 1e-15);
}

test "energy closure and final wet heat-capacity invariant are exact" {
    const contribution = testContribution(1, 0.25);
    const result = try accumulate(
        std.mem.zeroes(Totals),
        .{ .dry_heat_capacity_mj_per_k = 2 },
        &[_]Contribution{contribution},
    );
    try std.testing.expectEqual(@as(f64, 0), result.energy_balance_residual_mj_per_step);
    try std.testing.expectEqual(
        @as(f64, 2) + 4.19 * @as(f64, 0.25),
        result.wet_heat_capacity_mj_per_k,
    );
}

test "runtime species ranges remain independent" {
    const initial = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    const inputs = [_]SpeciesInputs{
        .{ .dry_heat_capacity_mj_per_k = 1 },
        .{ .dry_heat_capacity_mj_per_k = 2 },
    };
    const offsets = [_]usize{ 0, 2, 3 };
    const contributions = [_]Contribution{
        testContribution(1, 0.1),
        testContribution(2, 0.2),
        testContribution(10, 0.3),
    };
    var scratch: [2]Totals = undefined;
    var destination: [2]Totals = undefined;
    try accumulateRuntimeSpecies(
        &initial,
        &inputs,
        &offsets,
        &contributions,
        &scratch,
        &destination,
    );
    try std.testing.expectEqual(@as(f64, 15), destination[0].net_radiation_mj_per_step);
    try std.testing.expectEqual(@as(f64, 50), destination[1].net_radiation_mj_per_step);
}

test "later nonfinite contribution leaves all destinations unchanged" {
    const initial = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    const inputs = [_]SpeciesInputs{
        .{ .dry_heat_capacity_mj_per_k = 1 },
        .{ .dry_heat_capacity_mj_per_k = 2 },
    };
    const offsets = [_]usize{ 0, 1, 2 };
    var contributions = [_]Contribution{
        testContribution(1, 0.1),
        testContribution(2, 0.2),
    };
    contributions[1].latent_heat_mj_per_step = std.math.nan(f64);
    var scratch: [2]Totals = undefined;
    var destination = [_]Totals{ std.mem.zeroes(Totals), std.mem.zeroes(Totals) };
    destination[0].net_radiation_mj_per_step = 41;
    destination[1].net_radiation_mj_per_step = 42;
    try std.testing.expectError(
        error.NonFiniteStandingDeadEnergyAggregation,
        accumulateRuntimeSpecies(
            &initial,
            &inputs,
            &offsets,
            &contributions,
            &scratch,
            &destination,
        ),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].net_radiation_mj_per_step);
    try std.testing.expectEqual(@as(f64, 42), destination[1].net_radiation_mj_per_step);
}
