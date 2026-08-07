const std = @import("std");

pub const state_species_count = 42;
pub const transported_species_count = 41;
pub const silicic_acid_species_index = 33;
pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    salt_simulation: SaltSimulation,
    first_layer_fortran: usize,
    last_layer_fortran: usize,
    minimum_layer_thickness_m: f64,
    layer_thickness_m: []const f64,
    /// Layer-major T*FLS totals for all 42 state species (mol step-1).
    net_micropore_flux_mol_per_step: []const f64,
    /// Layer-major R*FXS and R*FLZ compact 41-species families. H4SiO4 is absent.
    macropore_exchange_flux_mol_per_step: []const f64,
    subsurface_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9569--9654.
/// Layers NU..NL update only when DLYR>DLYRM. All species use
/// `state + T*FLS + R*FXS + R*FLZ` except H4SiO4, whose source line 9638
/// deliberately uses only `state + THYSIS`.
pub fn apply(inputs: Inputs, micropore_content_mol: []f64) !usize {
    if (inputs.salt_simulation == .disabled) return 0;
    const layer_count = inputs.layer_thickness_m.len;
    try validate(inputs, micropore_content_mol, layer_count);

    var updated_layers: usize = 0;
    for (inputs.first_layer_fortran - 1..inputs.last_layer_fortran) |layer| {
        if (inputs.layer_thickness_m[layer] <= inputs.minimum_layer_thickness_m) continue;
        const state_offset = layer * state_species_count;
        const transport_offset = layer * transported_species_count;
        for (0..state_species_count) |species| {
            const state_index = state_offset + species;
            var result = micropore_content_mol[state_index] + inputs.net_micropore_flux_mol_per_step[state_index];
            if (species != silicic_acid_species_index) {
                const compact_species = if (species < silicic_acid_species_index) species else species - 1;
                const transport_index = transport_offset + compact_species;
                result += inputs.macropore_exchange_flux_mol_per_step[transport_index];
                result += inputs.subsurface_flux_mol_per_step[transport_index];
            }
            if (!std.math.isFinite(result)) return error.NonFiniteSoilMicroporeSoluteResult;
        }
    }

    for (inputs.first_layer_fortran - 1..inputs.last_layer_fortran) |layer| {
        if (inputs.layer_thickness_m[layer] <= inputs.minimum_layer_thickness_m) continue;
        const state_offset = layer * state_species_count;
        const transport_offset = layer * transported_species_count;
        for (0..state_species_count) |species| {
            const state_index = state_offset + species;
            micropore_content_mol[state_index] += inputs.net_micropore_flux_mol_per_step[state_index];
            if (species != silicic_acid_species_index) {
                const compact_species = if (species < silicic_acid_species_index) species else species - 1;
                const transport_index = transport_offset + compact_species;
                micropore_content_mol[state_index] += inputs.macropore_exchange_flux_mol_per_step[transport_index];
                micropore_content_mol[state_index] += inputs.subsurface_flux_mol_per_step[transport_index];
            }
        }
        updated_layers += 1;
    }
    return updated_layers;
}

fn validate(inputs: Inputs, state: []const f64, layer_count: usize) !void {
    if (inputs.first_layer_fortran == 0 or inputs.last_layer_fortran < inputs.first_layer_fortran or
        inputs.last_layer_fortran > layer_count)
        return error.SoilMicroporeSoluteLayerRangeMismatch;
    const state_count = std.math.mul(usize, layer_count, state_species_count) catch
        return error.SoilMicroporeSoluteDimensionMismatch;
    const transport_count = std.math.mul(usize, layer_count, transported_species_count) catch
        return error.SoilMicroporeSoluteDimensionMismatch;
    if (state.len != state_count or inputs.net_micropore_flux_mol_per_step.len != state_count or
        inputs.macropore_exchange_flux_mol_per_step.len != transport_count or
        inputs.subsurface_flux_mol_per_step.len != transport_count)
        return error.SoilMicroporeSoluteDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_layer_thickness_m))
        return error.NonFiniteSoilMicroporeSoluteInput;
    const slices = [_][]const f64{
        inputs.layer_thickness_m,
        inputs.net_micropore_flux_mol_per_step,
        inputs.macropore_exchange_flux_mol_per_step,
        inputs.subsurface_flux_mol_per_step,
        state,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilMicroporeSoluteInput;
}

fn inputsFor(thickness: []const f64, net: []const f64, exchange: []const f64, subsurface: []const f64) Inputs {
    return .{
        .salt_simulation = .enabled,
        .first_layer_fortran = 1,
        .last_layer_fortran = thickness.len,
        .minimum_layer_thickness_m = 0.1,
        .layer_thickness_m = thickness,
        .net_micropore_flux_mol_per_step = net,
        .macropore_exchange_flux_mol_per_step = exchange,
        .subsurface_flux_mol_per_step = subsurface,
    };
}

test "disabled salt simulation bypasses dormant soil topology" {
    const empty = [_]f64{};
    const inputs: Inputs = .{
        .salt_simulation = .disabled,
        .first_layer_fortran = 0,
        .last_layer_fortran = std.math.maxInt(usize),
        .minimum_layer_thickness_m = std.math.nan(f64),
        .layer_thickness_m = &empty,
        .net_micropore_flux_mol_per_step = &empty,
        .macropore_exchange_flux_mol_per_step = &empty,
        .subsurface_flux_mol_per_step = &empty,
    };
    try std.testing.expectEqual(@as(usize, 0), try apply(inputs, @constCast(&empty)));
}

test "TRNSFRS updates transported species but H4SiO4 only receives net flux" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{2} ** state_species_count;
    const exchange = [_]f64{3} ** transported_species_count;
    const subsurface = [_]f64{4} ** transported_species_count;
    var state = [_]f64{10} ** state_species_count;
    try std.testing.expectEqual(@as(usize, 1), try apply(inputsFor(&thickness, &net, &exchange, &subsurface), &state));
    try std.testing.expectEqual(@as(f64, 19), state[0]);
    try std.testing.expectEqual(@as(f64, 12), state[silicic_acid_species_index]);
    try std.testing.expectEqual(@as(f64, 19), state[41]);
}

test "runtime NU through NL range and strict thickness guard are preserved" {
    const thickness = [_]f64{ 0.2, 0.1, 0.3 };
    const net = [_]f64{2} ** (3 * state_species_count);
    const exchange = [_]f64{3} ** (3 * transported_species_count);
    const subsurface = [_]f64{4} ** (3 * transported_species_count);
    var state = [_]f64{10} ** (3 * state_species_count);
    var inputs = inputsFor(&thickness, &net, &exchange, &subsurface);
    inputs.first_layer_fortran = 2;
    try std.testing.expectEqual(@as(usize, 1), try apply(inputs, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
    try std.testing.expectEqual(@as(f64, 10), state[state_species_count]);
    try std.testing.expectEqual(@as(f64, 19), state[2 * state_species_count]);
}

test "invalid runtime layer range fails atomically" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{2} ** state_species_count;
    const exchange = [_]f64{3} ** transported_species_count;
    const subsurface = [_]f64{4} ** transported_species_count;
    var state = [_]f64{10} ** state_species_count;
    var inputs = inputsFor(&thickness, &net, &exchange, &subsurface);
    inputs.first_layer_fortran = 0;
    try std.testing.expectError(error.SoilMicroporeSoluteLayerRangeMismatch, apply(inputs, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
}

test "compact transported topology mismatch fails atomically" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{2} ** state_species_count;
    const short_exchange = [_]f64{3} ** (transported_species_count - 1);
    const subsurface = [_]f64{4} ** transported_species_count;
    var state = [_]f64{10} ** state_species_count;
    try std.testing.expectError(error.SoilMicroporeSoluteDimensionMismatch, apply(inputsFor(&thickness, &net, &short_exchange, &subsurface), &state));
}

test "late layer overflow leaves earlier state atomic" {
    const thickness = [_]f64{ 0.2, 0.2 };
    var net = [_]f64{0} ** (2 * state_species_count);
    net[net.len - 1] = std.math.floatMax(f64);
    const exchange = [_]f64{0} ** (2 * transported_species_count);
    const subsurface = [_]f64{0} ** (2 * transported_species_count);
    var state = [_]f64{0} ** (2 * state_species_count);
    state[state.len - 1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilMicroporeSoluteResult, apply(inputsFor(&thickness, &net, &exchange, &subsurface), &state));
    try std.testing.expectEqual(@as(f64, 0), state[0]);
}
