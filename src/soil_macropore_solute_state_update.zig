const std = @import("std");

pub const species_count = 41;

pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    salt_simulation: SaltSimulation,
    first_layer_fortran: usize,
    last_layer_fortran: usize,
    minimum_layer_thickness_m: f64,
    layer_thickness_m: []const f64,
    /// Layer-major T*FHS and R*FXS families (mol timestep-1).
    net_macropore_flux_mol_per_step: []const f64,
    macropore_to_micropore_exchange_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9671--9711.
/// Under inherited NU..NL and DLYR>DLYRM guards, the 41 macropore species
/// (no H4SiO4) use `state + T*FHS - R*FXS` in exact source order.
pub fn apply(inputs: Inputs, macropore_content_mol: []f64) !usize {
    if (inputs.salt_simulation == .disabled) return 0;

    const layer_count = inputs.layer_thickness_m.len;
    try validate(inputs, macropore_content_mol, layer_count);

    for (inputs.first_layer_fortran - 1..inputs.last_layer_fortran) |layer| {
        if (inputs.layer_thickness_m[layer] <= inputs.minimum_layer_thickness_m) continue;
        const offset = layer * species_count;
        for (0..species_count) |species| {
            const index = offset + species;
            const after_net = macropore_content_mol[index] + inputs.net_macropore_flux_mol_per_step[index];
            const result = after_net - inputs.macropore_to_micropore_exchange_mol_per_step[index];
            if (!std.math.isFinite(after_net) or !std.math.isFinite(result))
                return error.NonFiniteSoilMacroporeSoluteResult;
        }
    }

    var updated_layers: usize = 0;
    for (inputs.first_layer_fortran - 1..inputs.last_layer_fortran) |layer| {
        if (inputs.layer_thickness_m[layer] <= inputs.minimum_layer_thickness_m) continue;
        const offset = layer * species_count;
        for (0..species_count) |species| {
            const index = offset + species;
            macropore_content_mol[index] += inputs.net_macropore_flux_mol_per_step[index];
            macropore_content_mol[index] -= inputs.macropore_to_micropore_exchange_mol_per_step[index];
        }
        updated_layers += 1;
    }
    return updated_layers;
}

fn validate(inputs: Inputs, state: []const f64, layer_count: usize) !void {
    if (inputs.first_layer_fortran == 0 or inputs.last_layer_fortran < inputs.first_layer_fortran or
        inputs.last_layer_fortran > layer_count)
        return error.SoilMacroporeSoluteLayerRangeMismatch;
    const expected = std.math.mul(usize, layer_count, species_count) catch
        return error.SoilMacroporeSoluteDimensionMismatch;
    if (state.len != expected or inputs.net_macropore_flux_mol_per_step.len != expected or
        inputs.macropore_to_micropore_exchange_mol_per_step.len != expected)
        return error.SoilMacroporeSoluteDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_layer_thickness_m))
        return error.NonFiniteSoilMacroporeSoluteInput;
    const slices = [_][]const f64{
        inputs.layer_thickness_m,
        inputs.net_macropore_flux_mol_per_step,
        inputs.macropore_to_micropore_exchange_mol_per_step,
        state,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilMacroporeSoluteInput;
}

fn inputsFor(thickness: []const f64, net: []const f64, exchange: []const f64) Inputs {
    return .{
        .salt_simulation = .enabled,
        .first_layer_fortran = 1,
        .last_layer_fortran = thickness.len,
        .minimum_layer_thickness_m = 0.1,
        .layer_thickness_m = thickness,
        .net_macropore_flux_mol_per_step = net,
        .macropore_to_micropore_exchange_mol_per_step = exchange,
    };
}

test "disabled salt simulation bypasses dormant macropore storage" {
    const empty = [_]f64{};
    var state = [_]f64{};
    const inputs = Inputs{
        .salt_simulation = .disabled,
        .first_layer_fortran = 0,
        .last_layer_fortran = std.math.maxInt(usize),
        .minimum_layer_thickness_m = std.math.nan(f64),
        .layer_thickness_m = &empty,
        .net_macropore_flux_mol_per_step = &empty,
        .macropore_to_micropore_exchange_mol_per_step = &empty,
    };
    try std.testing.expectEqual(@as(usize, 0), try apply(inputs, &state));
}

test "TRNSFRS updates exact 41-species macropore topology" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{5} ** species_count;
    const exchange = [_]f64{2} ** species_count;
    var state = [_]f64{10} ** species_count;
    try std.testing.expectEqual(@as(usize, 1), try apply(inputsFor(&thickness, &net, &exchange), &state));
    try std.testing.expectEqual(@as(f64, 13), state[0]);
    try std.testing.expectEqual(@as(f64, 13), state[32]);
    try std.testing.expectEqual(@as(f64, 13), state[40]);
}

test "NU NL range and strict thickness guard are preserved" {
    const thickness = [_]f64{ 0.2, 0.1, 0.3 };
    const net = [_]f64{5} ** (3 * species_count);
    const exchange = [_]f64{2} ** (3 * species_count);
    var state = [_]f64{10} ** (3 * species_count);
    var inputs = inputsFor(&thickness, &net, &exchange);
    inputs.first_layer_fortran = 2;
    try std.testing.expectEqual(@as(usize, 1), try apply(inputs, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
    try std.testing.expectEqual(@as(f64, 10), state[species_count]);
    try std.testing.expectEqual(@as(f64, 13), state[2 * species_count]);
}

test "invalid macropore layer range fails atomically" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{5} ** species_count;
    const exchange = [_]f64{2} ** species_count;
    var state = [_]f64{10} ** species_count;
    var inputs = inputsFor(&thickness, &net, &exchange);
    inputs.first_layer_fortran = 0;
    try std.testing.expectError(error.SoilMacroporeSoluteLayerRangeMismatch, apply(inputs, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
}

test "runtime macropore state topology mismatch fails atomically" {
    const thickness = [_]f64{0.2};
    const short_net = [_]f64{5} ** (species_count - 1);
    const exchange = [_]f64{2} ** species_count;
    var state = [_]f64{10} ** species_count;
    try std.testing.expectError(error.SoilMacroporeSoluteDimensionMismatch, apply(inputsFor(&thickness, &short_net, &exchange), &state));
}

test "late macropore overflow leaves prior layers atomic" {
    const thickness = [_]f64{ 0.2, 0.2 };
    var net = [_]f64{0} ** (2 * species_count);
    var exchange = [_]f64{0} ** (2 * species_count);
    net[net.len - 1] = std.math.floatMax(f64);
    exchange[exchange.len - 1] = -std.math.floatMax(f64);
    var state = [_]f64{0} ** (2 * species_count);
    try std.testing.expectError(error.NonFiniteSoilMacroporeSoluteResult, apply(inputsFor(&thickness, &net, &exchange), &state));
    try std.testing.expectEqual(@as(f64, 0), state[0]);
}
