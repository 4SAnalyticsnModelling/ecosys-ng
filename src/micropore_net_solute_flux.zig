const std = @import("std");

pub const species_count = 42;

pub const PoreFluxApplicability = enum {
    apply,
    skip,
};

pub const Inputs = struct {
    applicability: PoreFluxApplicability,
    current_layer_thickness_m: f64,
    neighbor_layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
    /// Convective plus diffusive micropore fluxes in TRNSFRS species order.
    local_flux_mol_per_step: []const f64,
    neighbor_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9047--9133.
/// When both current and resolved-neighbor layers exceed DLYRM, all 42
/// micropore species are accumulated as `T = T + R(local) - R(neighbor)`.
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
            return error.NonFiniteMicroporeNetSoluteResult;
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
        if (!std.math.isFinite(value)) return error.NonFiniteMicroporeNetSoluteInput;
}

fn validateActiveFlux(inputs: Inputs, totals: []const f64) !void {
    if (inputs.local_flux_mol_per_step.len != species_count or
        inputs.neighbor_flux_mol_per_step.len != species_count or
        totals.len != species_count)
        return error.MicroporeNetSoluteDimensionMismatch;
    const slices = [_][]const f64{
        inputs.local_flux_mol_per_step,
        inputs.neighbor_flux_mol_per_step,
        totals,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteMicroporeNetSoluteInput;
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

test "TRNSFRS accumulates all 42 micropore species in source order" {
    const local = [_]f64{5} ** species_count;
    const neighbor = [_]f64{2} ** species_count;
    var totals = [_]f64{1} ** species_count;
    try std.testing.expect(try accumulate(validInputs(&local, &neighbor), &totals));
    try std.testing.expectEqual(@as(f64, 4), totals[0]);
    try std.testing.expectEqual(@as(f64, 4), totals[33]);
    try std.testing.expectEqual(@as(f64, 4), totals[41]);
}

test "strict current and neighbor thickness guards skip equal threshold" {
    const local = [_]f64{};
    const neighbor = [_]f64{};
    var totals = [_]f64{1} ** species_count;
    var inputs = validInputs(&local, &neighbor);
    inputs.neighbor_layer_thickness_m = inputs.minimum_layer_thickness_m;
    try std.testing.expect(!try accumulate(inputs, &totals));
    try std.testing.expectEqual(@as(f64, 1), totals[0]);
    inputs = validInputs(&local, &neighbor);
    inputs.current_layer_thickness_m = inputs.minimum_layer_thickness_m;
    try std.testing.expect(!try accumulate(inputs, &totals));
}

test "outer pore-flux applicability skips otherwise valid layers" {
    const local = [_]f64{};
    const neighbor = [_]f64{};
    var totals = [_]f64{};
    var inputs = validInputs(&local, &neighbor);
    inputs.applicability = .skip;
    try std.testing.expect(!try accumulate(inputs, &totals));
}

test "runtime micropore topology mismatch fails atomically" {
    const short_local = [_]f64{5} ** (species_count - 1);
    const neighbor = [_]f64{2} ** species_count;
    var totals = [_]f64{1} ** species_count;
    try std.testing.expectError(error.MicroporeNetSoluteDimensionMismatch, accumulate(validInputs(&short_local, &neighbor), &totals));
    try std.testing.expectEqual(@as(f64, 1), totals[0]);
}

test "late micropore overflow leaves all totals atomic" {
    var local = [_]f64{0} ** species_count;
    var neighbor = [_]f64{0} ** species_count;
    local[species_count - 1] = std.math.floatMax(f64);
    neighbor[species_count - 1] = -std.math.floatMax(f64);
    var totals = [_]f64{0} ** species_count;
    try std.testing.expectError(error.NonFiniteMicroporeNetSoluteResult, accumulate(validInputs(&local, &neighbor), &totals));
    try std.testing.expectEqual(@as(f64, 0), totals[0]);
}
