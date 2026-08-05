const std = @import("std");

pub const species_count = 49;

pub const Inputs = struct {
    convective_amount_per_step: []const f64,
    diffusive_amount_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 6764--6812.
/// The result uses the compact macro/micropore exchange topology: 33 salts and
/// complexes without H4SiO4, followed by 8 non-band and 8 band P species.
pub fn calculate(inputs: Inputs) ![species_count]f64 {
    if (inputs.convective_amount_per_step.len != species_count or
        inputs.diffusive_amount_per_step.len != species_count)
        return error.SoilPoreTotalSoluteExchangeDimensionMismatch;
    for (inputs.convective_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreTotalSoluteExchangeInput;
    for (inputs.diffusive_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreTotalSoluteExchangeInput;

    var total: [species_count]f64 = undefined;
    for (0..species_count) |species| {
        total[species] = inputs.convective_amount_per_step[species] + inputs.diffusive_amount_per_step[species];
        if (!std.math.isFinite(total[species])) return error.NonFiniteSoilPoreTotalSoluteExchangeResult;
    }
    return total;
}

test "TRNSFRS combines all compact pore-exchange species in source order" {
    var convective: [species_count]f64 = undefined;
    var diffusive: [species_count]f64 = undefined;
    for (0..species_count) |species| {
        convective[species] = @floatFromInt(species);
        diffusive[species] = 0.5;
    }
    const total = try calculate(.{ .convective_amount_per_step = &convective, .diffusive_amount_per_step = &diffusive });
    try std.testing.expectEqual(@as(f64, 0.5), total[0]);
    try std.testing.expectEqual(@as(f64, 32.5), total[32]);
    try std.testing.expectEqual(@as(f64, 33.5), total[33]);
    try std.testing.expectEqual(@as(f64, 48.5), total[48]);
}

test "TRNSFRS retains signed cancellation between convection and diffusion" {
    const convective = [_]f64{3} ** species_count;
    const diffusive = [_]f64{-3} ** species_count;
    const total = try calculate(.{ .convective_amount_per_step = &convective, .diffusive_amount_per_step = &diffusive });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &total);
}

test "compact topology has no H4SiO4 exchange slot" {
    const convective = [_]f64{1} ** species_count;
    const diffusive = [_]f64{2} ** species_count;
    const total = try calculate(.{ .convective_amount_per_step = &convective, .diffusive_amount_per_step = &diffusive });
    try std.testing.expectEqual(@as(usize, 49), total.len);
    try std.testing.expectEqual(@as(f64, 3), total[33]);
}

test "dimension mismatch fails before publishing a result" {
    const short = [_]f64{0} ** (species_count - 1);
    const complete = [_]f64{0} ** species_count;
    try std.testing.expectError(error.SoilPoreTotalSoluteExchangeDimensionMismatch, calculate(.{ .convective_amount_per_step = &short, .diffusive_amount_per_step = &complete }));
}

test "late non-finite diffusive input fails atomically" {
    const convective = [_]f64{1} ** species_count;
    var diffusive = [_]f64{2} ** species_count;
    diffusive[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilPoreTotalSoluteExchangeInput, calculate(.{ .convective_amount_per_step = &convective, .diffusive_amount_per_step = &diffusive }));
}
