const std = @import("std");

/// 33 salt/complex species (no H4SiO4), eight non-band phosphate species,
/// and eight band phosphate species.
pub const exchange_species_count = 49;

/// Exact source-order translation of TRNSFRS.F lines 3625--3673.
/// Positive mol step-1 exchange is macropore to micropore. The legacy comment
/// lists H4SiO4, but no H4SiO4 assignment exists in this compact source block.
pub fn combine(
    convective_flux_mol_per_step: []const f64,
    diffusive_flux_mol_per_step: []const f64,
    total_flux_mol_per_step: []f64,
) !void {
    if (convective_flux_mol_per_step.len != exchange_species_count or
        diffusive_flux_mol_per_step.len != exchange_species_count or
        total_flux_mol_per_step.len != exchange_species_count)
        return error.SurfacePoreTotalSoluteExchangeDimensionMismatch;
    for (0..exchange_species_count) |species| {
        const convective = convective_flux_mol_per_step[species];
        const diffusive = diffusive_flux_mol_per_step[species];
        if (!std.math.isFinite(convective) or !std.math.isFinite(diffusive))
            return error.NonFiniteSurfacePoreTotalSoluteExchangeInput;
        if (!std.math.isFinite(convective + diffusive))
            return error.NonFiniteSurfacePoreTotalSoluteExchangeResult;
    }
    for (0..exchange_species_count) |species|
        total_flux_mol_per_step[species] =
            convective_flux_mol_per_step[species] + diffusive_flux_mol_per_step[species];
}

test "TRNSFRS combines every compact pore-exchange field in source order" {
    var convective: [exchange_species_count]f64 = undefined;
    var diffusive: [exchange_species_count]f64 = undefined;
    for (&convective, &diffusive, 0..) |*convective_value, *diffusive_value, species| {
        convective_value.* = @floatFromInt(species + 1);
        diffusive_value.* = -0.25 * @as(f64, @floatFromInt(species + 1));
    }
    var total = [_]f64{0} ** exchange_species_count;
    try combine(&convective, &diffusive, &total);
    try std.testing.expectEqual(@as(f64, 0.75), total[0]);
    try std.testing.expectEqual(@as(f64, 24.75), total[32]);
    try std.testing.expectEqual(@as(f64, 25.5), total[33]);
    try std.testing.expectEqual(@as(f64, 36.75), total[48]);
}

test "compact total topology rejects a phantom H4SiO4 slot" {
    const convective = [_]f64{1} ** (exchange_species_count + 1);
    const diffusive = [_]f64{2} ** exchange_species_count;
    var total = [_]f64{9} ** exchange_species_count;
    try std.testing.expectError(
        error.SurfacePoreTotalSoluteExchangeDimensionMismatch,
        combine(&convective, &diffusive, &total),
    );
}

test "late non-finite exchange leaves total atomic" {
    const convective = [_]f64{1} ** exchange_species_count;
    var diffusive = [_]f64{2} ** exchange_species_count;
    diffusive[exchange_species_count - 1] = std.math.inf(f64);
    var total = [_]f64{9} ** exchange_species_count;
    try std.testing.expectError(
        error.NonFiniteSurfacePoreTotalSoluteExchangeInput,
        combine(&convective, &diffusive, &total),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** exchange_species_count), &total);
}

test "overflowing exchange sum leaves total atomic" {
    const convective = [_]f64{std.math.floatMax(f64)} ** exchange_species_count;
    const diffusive = [_]f64{std.math.floatMax(f64)} ** exchange_species_count;
    var total = [_]f64{9} ** exchange_species_count;
    try std.testing.expectError(
        error.NonFiniteSurfacePoreTotalSoluteExchangeResult,
        combine(&convective, &diffusive, &total),
    );
    try std.testing.expectEqual(@as(f64, 9), total[0]);
}
