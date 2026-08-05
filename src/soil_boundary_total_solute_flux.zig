const std = @import("std");

pub const micropore_species_count = 50;
pub const macropore_species_count = 49;

pub const Inputs = struct {
    micropore_convective_amount_per_step: []const f64,
    micropore_diffusive_amount_per_step: []const f64,
    macropore_convective_amount_per_step: []const f64,
    macropore_diffusive_amount_per_step: []const f64,
};

pub const Result = struct {
    micropore_total_amount_per_step: [micropore_species_count]f64,
    macropore_total_amount_per_step: [macropore_species_count]f64,
};

/// Compatibility translation of TRNSFRS.F lines 6061--6159.
/// The exact legacy assignment order is convective plus diffusive for each
/// micropore species first, followed by each compact macropore species.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    var result: Result = undefined;
    for (0..micropore_species_count) |species| {
        result.micropore_total_amount_per_step[species] =
            inputs.micropore_convective_amount_per_step[species] + inputs.micropore_diffusive_amount_per_step[species];
        if (!std.math.isFinite(result.micropore_total_amount_per_step[species]))
            return error.NonFiniteSoilBoundaryTotalSoluteResult;
    }
    for (0..macropore_species_count) |species| {
        result.macropore_total_amount_per_step[species] =
            inputs.macropore_convective_amount_per_step[species] + inputs.macropore_diffusive_amount_per_step[species];
        if (!std.math.isFinite(result.macropore_total_amount_per_step[species]))
            return error.NonFiniteSoilBoundaryTotalSoluteResult;
    }
    return result;
}

fn validate(inputs: Inputs) !void {
    if (inputs.micropore_convective_amount_per_step.len != micropore_species_count or
        inputs.micropore_diffusive_amount_per_step.len != micropore_species_count or
        inputs.macropore_convective_amount_per_step.len != macropore_species_count or
        inputs.macropore_diffusive_amount_per_step.len != macropore_species_count)
        return error.SoilBoundaryTotalSoluteDimensionMismatch;
    for (inputs.micropore_convective_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryTotalSoluteInput;
    for (inputs.micropore_diffusive_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryTotalSoluteInput;
    for (inputs.macropore_convective_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryTotalSoluteInput;
    for (inputs.macropore_diffusive_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryTotalSoluteInput;
}

test "TRNSFRS assembles micropore and compact macropore totals independently" {
    var micro_convective: [micropore_species_count]f64 = undefined;
    var micro_diffusive: [micropore_species_count]f64 = undefined;
    var macro_convective: [macropore_species_count]f64 = undefined;
    var macro_diffusive: [macropore_species_count]f64 = undefined;
    for (0..micropore_species_count) |species| {
        micro_convective[species] = @floatFromInt(species);
        micro_diffusive[species] = 0.5;
    }
    for (0..macropore_species_count) |species| {
        macro_convective[species] = -@as(f64, @floatFromInt(species));
        macro_diffusive[species] = 0.25;
    }
    const result = try calculate(.{ .micropore_convective_amount_per_step = &micro_convective, .micropore_diffusive_amount_per_step = &micro_diffusive, .macropore_convective_amount_per_step = &macro_convective, .macropore_diffusive_amount_per_step = &macro_diffusive });
    try std.testing.expectEqual(@as(f64, 0.5), result.micropore_total_amount_per_step[0]);
    try std.testing.expectEqual(@as(f64, 49.5), result.micropore_total_amount_per_step[49]);
    try std.testing.expectEqual(@as(f64, -47.75), result.macropore_total_amount_per_step[48]);
}

test "TRNSFRS keeps micropore H4SiO4 without adding a macropore slot" {
    var micro_convective = [_]f64{0} ** micropore_species_count;
    const micro_diffusive = [_]f64{0} ** micropore_species_count;
    const macro_convective = [_]f64{0} ** macropore_species_count;
    const macro_diffusive = [_]f64{0} ** macropore_species_count;
    micro_convective[33] = 7;
    const result = try calculate(.{ .micropore_convective_amount_per_step = &micro_convective, .micropore_diffusive_amount_per_step = &micro_diffusive, .macropore_convective_amount_per_step = &macro_convective, .macropore_diffusive_amount_per_step = &macro_diffusive });
    try std.testing.expectEqual(@as(f64, 7), result.micropore_total_amount_per_step[33]);
    try std.testing.expectEqual(@as(usize, 49), result.macropore_total_amount_per_step.len);
}

test "TRNSFRS permits opposing convection and diffusion cancellation" {
    const micro_convective = [_]f64{2} ** micropore_species_count;
    const micro_diffusive = [_]f64{-2} ** micropore_species_count;
    const macro_convective = [_]f64{-3} ** macropore_species_count;
    const macro_diffusive = [_]f64{3} ** macropore_species_count;
    const result = try calculate(.{ .micropore_convective_amount_per_step = &micro_convective, .micropore_diffusive_amount_per_step = &micro_diffusive, .macropore_convective_amount_per_step = &macro_convective, .macropore_diffusive_amount_per_step = &macro_diffusive });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** micropore_species_count), &result.micropore_total_amount_per_step);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** macropore_species_count), &result.macropore_total_amount_per_step);
}

test "dimension mismatch fails before any result is published" {
    const short_micro = [_]f64{0} ** (micropore_species_count - 1);
    const micro = [_]f64{0} ** micropore_species_count;
    const macro = [_]f64{0} ** macropore_species_count;
    try std.testing.expectError(error.SoilBoundaryTotalSoluteDimensionMismatch, calculate(.{ .micropore_convective_amount_per_step = &short_micro, .micropore_diffusive_amount_per_step = &micro, .macropore_convective_amount_per_step = &macro, .macropore_diffusive_amount_per_step = &macro }));
}

test "late non-finite input fails atomically" {
    const micro = [_]f64{0} ** micropore_species_count;
    const macro_convective = [_]f64{0} ** macropore_species_count;
    var macro_diffusive = [_]f64{0} ** macropore_species_count;
    macro_diffusive[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundaryTotalSoluteInput, calculate(.{ .micropore_convective_amount_per_step = &micro, .micropore_diffusive_amount_per_step = &micro, .macropore_convective_amount_per_step = &macro_convective, .macropore_diffusive_amount_per_step = &macro_diffusive }));
}
