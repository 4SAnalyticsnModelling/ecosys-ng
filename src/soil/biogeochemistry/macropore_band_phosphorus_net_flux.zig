const std = @import("std");

pub const species_count = 8;
pub const PoreFluxApplicability = enum { apply, skip };

pub const Inputs = struct {
    applicability: PoreFluxApplicability,
    current_layer_thickness_m: f64,
    neighbor_layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
    /// H0P,H3P,F1P,F2P,C0P,C1P,C2P,M1P macropore-band fluxes.
    local_flux_mol_per_step: []const f64,
    neighbor_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9232--9247.
/// The eight BHB terms use the containing pore block's strict layer guards
/// and exact source evaluation `T = T + R(local) - R(neighbor)`.
pub fn accumulate(inputs: Inputs, totals_mol_per_step: []f64) !bool {
    if (inputs.applicability == .skip) return false;
    try validateThickness(inputs);
    if (inputs.current_layer_thickness_m <= inputs.minimum_layer_thickness_m or
        inputs.neighbor_layer_thickness_m <= inputs.minimum_layer_thickness_m) return false;
    try validateActiveFlux(inputs, totals_mol_per_step);

    for (0..species_count) |species| {
        const after_local = totals_mol_per_step[species] + inputs.local_flux_mol_per_step[species];
        const result = after_local - inputs.neighbor_flux_mol_per_step[species];
        if (!std.math.isFinite(after_local) or !std.math.isFinite(result))
            return error.NonFiniteMacroporeBandPhosphorusResult;
    }
    for (0..species_count) |species| {
        totals_mol_per_step[species] = totals_mol_per_step[species] +
            inputs.local_flux_mol_per_step[species];
        totals_mol_per_step[species] = totals_mol_per_step[species] -
            inputs.neighbor_flux_mol_per_step[species];
    }
    return true;
}

fn validateThickness(inputs: Inputs) !void {
    const scalars = [_]f64{
        inputs.current_layer_thickness_m,
        inputs.neighbor_layer_thickness_m,
        inputs.minimum_layer_thickness_m,
    };
    for (scalars) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteMacroporeBandPhosphorusInput;
}

fn validateActiveFlux(inputs: Inputs, totals: []const f64) !void {
    if (inputs.local_flux_mol_per_step.len != species_count or
        inputs.neighbor_flux_mol_per_step.len != species_count or totals.len != species_count)
        return error.MacroporeBandPhosphorusDimensionMismatch;
    const slices = [_][]const f64{ inputs.local_flux_mol_per_step, inputs.neighbor_flux_mol_per_step, totals };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteMacroporeBandPhosphorusInput;
}

fn validInputs(local: []const f64, neighbor: []const f64) Inputs {
    return .{
        .applicability = .apply,
        .current_layer_thickness_m = 0.2,
        .neighbor_layer_thickness_m = 0.3,
        .minimum_layer_thickness_m = 0.1,
        .local_flux_mol_per_step = local,
        .neighbor_flux_mol_per_step = neighbor,
    };
}

test "TRNSFRS accumulates all eight BHB species" {
    const local = [_]f64{5} ** species_count;
    const neighbor = [_]f64{2} ** species_count;
    var totals = [_]f64{1} ** species_count;
    try std.testing.expect(try accumulate(validInputs(&local, &neighbor), &totals));
    try std.testing.expectEqualSlices(f64, &([_]f64{4} ** species_count), &totals);
}

test "BHB terms inherit strict layer guards" {
    const local = [_]f64{};
    const neighbor = [_]f64{};
    var totals = [_]f64{1} ** species_count;
    var inputs = validInputs(&local, &neighbor);
    inputs.neighbor_layer_thickness_m = inputs.minimum_layer_thickness_m;
    try std.testing.expect(!try accumulate(inputs, &totals));
    try std.testing.expectEqual(@as(f64, 1), totals[0]);
}

test "skipped pore block leaves BHB totals unchanged" {
    const local = [_]f64{};
    const neighbor = [_]f64{};
    var totals = [_]f64{};
    var inputs = validInputs(&local, &neighbor);
    inputs.applicability = .skip;
    try std.testing.expect(!try accumulate(inputs, &totals));
}

test "runtime BHB topology mismatch fails atomically" {
    const short_local = [_]f64{5} ** (species_count - 1);
    const neighbor = [_]f64{2} ** species_count;
    var totals = [_]f64{1} ** species_count;
    try std.testing.expectError(error.MacroporeBandPhosphorusDimensionMismatch, accumulate(validInputs(&short_local, &neighbor), &totals));
    try std.testing.expectEqual(@as(f64, 1), totals[0]);
}

test "late BHB overflow leaves all totals atomic" {
    var local = [_]f64{0} ** species_count;
    var neighbor = [_]f64{0} ** species_count;
    local[species_count - 1] = std.math.floatMax(f64);
    neighbor[species_count - 1] = -std.math.floatMax(f64);
    var totals = [_]f64{0} ** species_count;
    try std.testing.expectError(error.NonFiniteMacroporeBandPhosphorusResult, accumulate(validInputs(&local, &neighbor), &totals));
    try std.testing.expectEqual(@as(f64, 0), totals[0]);
}
