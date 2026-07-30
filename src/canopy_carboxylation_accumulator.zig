const std = @import("std");

/// STOMATE.F line 87 and source additions at lines 347, 395, 570, and 617.
/// Each species total is reset exactly once, then receives its complete
/// caller-ordered branch/node/sample contribution range.
pub fn accumulateRuntimeSpecies(
    contribution_offsets: []const usize,
    contributions_umol_per_s: []const f64,
    scratch_umol_per_s: []f64,
    destination_umol_per_s: []f64,
) !void {
    if (scratch_umol_per_s.len != destination_umol_per_s.len or
        contribution_offsets.len != destination_umol_per_s.len + 1)
        return error.CanopyCarboxylationAccumulatorDimensionMismatch;
    if (contribution_offsets[0] != 0 or
        contribution_offsets[contribution_offsets.len - 1] !=
            contributions_umol_per_s.len)
        return error.InvalidCanopyCarboxylationContributionOffsets;

    for (scratch_umol_per_s, 0..) |*candidate, species| {
        const begin = contribution_offsets[species];
        const end = contribution_offsets[species + 1];
        if (begin > end or end > contributions_umol_per_s.len)
            return error.InvalidCanopyCarboxylationContributionOffsets;
        candidate.* = 0;
        for (contributions_umol_per_s[begin..end]) |contribution_umol_per_s| {
            if (!std.math.isFinite(contribution_umol_per_s))
                return error.NonFiniteCanopyCarboxylationContribution;
            if (contribution_umol_per_s < 0)
                return error.InvalidCanopyCarboxylationContribution;
            candidate.* += contribution_umol_per_s;
            if (!std.math.isFinite(candidate.*))
                return error.NonFiniteCanopyCarboxylationTotal;
        }
    }
    @memcpy(destination_umol_per_s, scratch_umol_per_s);
}

test "species accumulator resets once before all topology contributions" {
    const offsets = [_]usize{ 0, 4, 7 };
    const contributions = [_]f64{
        1,  2,  3,  4,
        10, 20, 30,
    };
    var scratch: [2]f64 = undefined;
    var destination = [_]f64{ 900, 900 };
    try accumulateRuntimeSpecies(
        &offsets,
        &contributions,
        &scratch,
        &destination,
    );
    try std.testing.expectEqualSlices(f64, &[_]f64{ 10, 60 }, &destination);
}

test "species with no active samples receives exact source reset zero" {
    const offsets = [_]usize{ 0, 0, 1 };
    const contributions = [_]f64{5};
    var scratch: [2]f64 = undefined;
    var destination = [_]f64{ 41, 42 };
    try accumulateRuntimeSpecies(
        &offsets,
        &contributions,
        &scratch,
        &destination,
    );
    try std.testing.expectEqualSlices(f64, &[_]f64{ 0, 5 }, &destination);
}

test "accumulation retains caller source order" {
    const offsets = [_]usize{ 0, 3 };
    const contributions = [_]f64{ 1.0e16, 1, 1 };
    var scratch: [1]f64 = undefined;
    var destination: [1]f64 = undefined;
    try accumulateRuntimeSpecies(
        &offsets,
        &contributions,
        &scratch,
        &destination,
    );
    var expected: f64 = 0;
    for (contributions) |contribution| expected += contribution;
    try std.testing.expectEqual(expected, destination[0]);
}

test "later nonfinite species contribution leaves all totals unchanged" {
    const offsets = [_]usize{ 0, 1, 2 };
    const contributions = [_]f64{ 5, std.math.nan(f64) };
    var scratch: [2]f64 = undefined;
    var destination = [_]f64{ 41, 42 };
    try std.testing.expectError(
        error.NonFiniteCanopyCarboxylationContribution,
        accumulateRuntimeSpecies(
            &offsets,
            &contributions,
            &scratch,
            &destination,
        ),
    );
    try std.testing.expectEqualSlices(f64, &[_]f64{ 41, 42 }, &destination);
}
