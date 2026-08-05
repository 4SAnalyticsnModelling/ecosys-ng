const std = @import("std");

pub const micropore_species_count = 42;
pub const macropore_species_count = 41;
pub const band_phosphorus_species_count = 8;
pub const PoreFluxApplicability = enum { apply, skip };

pub const Inputs = struct {
    applicability: PoreFluxApplicability,
    current_layer_thickness_m: f64,
    neighbor_layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
};

pub const Totals = struct {
    micropore_mol_per_step: []f64,
    micropore_band_phosphorus_mol_per_step: []f64,
    macropore_mol_per_step: []f64,
    macropore_band_phosphorus_mol_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 9249--9347.
/// This is the ELSE branch of the dual DLYR>DLYRM test: if either layer is
/// too thin, all four exact source families are reset to zero.
pub fn reset(inputs: Inputs, totals: Totals) !bool {
    if (inputs.applicability == .skip) return false;
    try validateThickness(inputs);
    if (inputs.current_layer_thickness_m > inputs.minimum_layer_thickness_m and
        inputs.neighbor_layer_thickness_m > inputs.minimum_layer_thickness_m) return false;
    try validateActiveTotals(totals);

    @memset(totals.micropore_mol_per_step, 0);
    @memset(totals.micropore_band_phosphorus_mol_per_step, 0);
    @memset(totals.macropore_mol_per_step, 0);
    @memset(totals.macropore_band_phosphorus_mol_per_step, 0);
    return true;
}

fn validateThickness(inputs: Inputs) !void {
    const scalars = [_]f64{
        inputs.current_layer_thickness_m,
        inputs.neighbor_layer_thickness_m,
        inputs.minimum_layer_thickness_m,
    };
    for (scalars) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteThinLayerPoreSoluteInput;
}

fn validateActiveTotals(totals: Totals) !void {
    if (totals.micropore_mol_per_step.len != micropore_species_count or
        totals.micropore_band_phosphorus_mol_per_step.len != band_phosphorus_species_count or
        totals.macropore_mol_per_step.len != macropore_species_count or
        totals.macropore_band_phosphorus_mol_per_step.len != band_phosphorus_species_count)
        return error.ThinLayerPoreSoluteDimensionMismatch;
    const slices = [_][]const f64{
        totals.micropore_mol_per_step,
        totals.micropore_band_phosphorus_mol_per_step,
        totals.macropore_mol_per_step,
        totals.macropore_band_phosphorus_mol_per_step,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteThinLayerPoreSoluteInput;
}

fn totalsFor(micropore: []f64, micro_band: []f64, macropore: []f64, macro_band: []f64) Totals {
    return .{
        .micropore_mol_per_step = micropore,
        .micropore_band_phosphorus_mol_per_step = micro_band,
        .macropore_mol_per_step = macropore,
        .macropore_band_phosphorus_mol_per_step = macro_band,
    };
}

test "TRNSFRS thin current layer resets every exact pore family" {
    var micropore = [_]f64{1} ** micropore_species_count;
    var micro_band = [_]f64{2} ** band_phosphorus_species_count;
    var macropore = [_]f64{3} ** macropore_species_count;
    var macro_band = [_]f64{4} ** band_phosphorus_species_count;
    try std.testing.expect(try reset(.{
        .applicability = .apply,
        .current_layer_thickness_m = 0.1,
        .neighbor_layer_thickness_m = 0.2,
        .minimum_layer_thickness_m = 0.1,
    }, totalsFor(&micropore, &micro_band, &macropore, &macro_band)));
    try std.testing.expectEqual(@as(f64, 0), micropore[41]);
    try std.testing.expectEqual(@as(f64, 0), micro_band[7]);
    try std.testing.expectEqual(@as(f64, 0), macropore[40]);
    try std.testing.expectEqual(@as(f64, 0), macro_band[7]);
}

test "thin neighbor alone triggers reset" {
    var micropore = [_]f64{1} ** micropore_species_count;
    var micro_band = [_]f64{2} ** band_phosphorus_species_count;
    var macropore = [_]f64{3} ** macropore_species_count;
    var macro_band = [_]f64{4} ** band_phosphorus_species_count;
    _ = try reset(.{
        .applicability = .apply,
        .current_layer_thickness_m = 0.2,
        .neighbor_layer_thickness_m = 0.05,
        .minimum_layer_thickness_m = 0.1,
    }, totalsFor(&micropore, &micro_band, &macropore, &macro_band));
    try std.testing.expectEqual(@as(f64, 0), micropore[0]);
}

test "two thick layers and skipped outer block do not reset" {
    var micropore = [_]f64{1} ** micropore_species_count;
    var micro_band = [_]f64{2} ** band_phosphorus_species_count;
    var macropore = [_]f64{3} ** macropore_species_count;
    var macro_band = [_]f64{4} ** band_phosphorus_species_count;
    const totals = totalsFor(&micropore, &micro_band, &macropore, &macro_band);
    try std.testing.expect(!try reset(.{ .applicability = .apply, .current_layer_thickness_m = 0.2, .neighbor_layer_thickness_m = 0.3, .minimum_layer_thickness_m = 0.1 }, totals));
    const empty = [_]f64{};
    try std.testing.expect(!try reset(.{ .applicability = .skip, .current_layer_thickness_m = std.math.nan(f64), .neighbor_layer_thickness_m = std.math.nan(f64), .minimum_layer_thickness_m = std.math.nan(f64) }, totalsFor(@constCast(&empty), @constCast(&empty), @constCast(&empty), @constCast(&empty))));
    try std.testing.expectEqual(@as(f64, 1), micropore[0]);
}

test "runtime reset topology mismatch fails atomically" {
    var short_micropore = [_]f64{1} ** (micropore_species_count - 1);
    var micro_band = [_]f64{2} ** band_phosphorus_species_count;
    var macropore = [_]f64{3} ** macropore_species_count;
    var macro_band = [_]f64{4} ** band_phosphorus_species_count;
    try std.testing.expectError(error.ThinLayerPoreSoluteDimensionMismatch, reset(.{ .applicability = .apply, .current_layer_thickness_m = 0.01, .neighbor_layer_thickness_m = 0.2, .minimum_layer_thickness_m = 0.1 }, totalsFor(&short_micropore, &micro_band, &macropore, &macro_band)));
    try std.testing.expectEqual(@as(f64, 3), macropore[0]);
}

test "nonfinite late family fails before reset" {
    var micropore = [_]f64{1} ** micropore_species_count;
    var micro_band = [_]f64{2} ** band_phosphorus_species_count;
    var macropore = [_]f64{3} ** macropore_species_count;
    var macro_band = [_]f64{4} ** band_phosphorus_species_count;
    macro_band[7] = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteThinLayerPoreSoluteInput, reset(.{ .applicability = .apply, .current_layer_thickness_m = 0.01, .neighbor_layer_thickness_m = 0.2, .minimum_layer_thickness_m = 0.1 }, totalsFor(&micropore, &micro_band, &macropore, &macro_band)));
    try std.testing.expectEqual(@as(f64, 1), micropore[0]);
}
