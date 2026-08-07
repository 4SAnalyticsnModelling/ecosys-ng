const std = @import("std");

pub const species_count = 8;

pub const SaltSimulation = enum { disabled, enabled };

pub const Inputs = struct {
    salt_simulation: SaltSimulation,
    first_layer_fortran: usize,
    last_layer_fortran: usize,
    minimum_layer_thickness_m: f64,
    layer_thickness_m: []const f64,
    /// Layer-major BFB, BXB, and BBZ families (mol timestep-1).
    net_micropore_band_flux_mol_per_step: []const f64,
    macropore_exchange_flux_mol_per_step: []const f64,
    band_subsurface_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 9655--9670.
/// Under the inherited NU..NL and DLYR>DLYRM guards, all eight banded
/// phosphorus states use `state + T*BFB + R*BXB + R*BBZ`.
pub fn apply(inputs: Inputs, band_content_mol: []f64) !usize {
    if (inputs.salt_simulation == .disabled) return 0;

    const layer_count = inputs.layer_thickness_m.len;
    try validate(inputs, band_content_mol, layer_count);

    for (inputs.first_layer_fortran - 1..inputs.last_layer_fortran) |layer| {
        if (inputs.layer_thickness_m[layer] <= inputs.minimum_layer_thickness_m) continue;
        const offset = layer * species_count;
        for (0..species_count) |species| {
            const index = offset + species;
            const after_net = band_content_mol[index] + inputs.net_micropore_band_flux_mol_per_step[index];
            const after_exchange = after_net + inputs.macropore_exchange_flux_mol_per_step[index];
            const result = after_exchange + inputs.band_subsurface_flux_mol_per_step[index];
            if (!std.math.isFinite(after_net) or !std.math.isFinite(after_exchange) or !std.math.isFinite(result))
                return error.NonFiniteMicroporeBandPhosphorusResult;
        }
    }

    var updated_layers: usize = 0;
    for (inputs.first_layer_fortran - 1..inputs.last_layer_fortran) |layer| {
        if (inputs.layer_thickness_m[layer] <= inputs.minimum_layer_thickness_m) continue;
        const offset = layer * species_count;
        for (0..species_count) |species| {
            const index = offset + species;
            band_content_mol[index] += inputs.net_micropore_band_flux_mol_per_step[index];
            band_content_mol[index] += inputs.macropore_exchange_flux_mol_per_step[index];
            band_content_mol[index] += inputs.band_subsurface_flux_mol_per_step[index];
        }
        updated_layers += 1;
    }
    return updated_layers;
}

fn validate(inputs: Inputs, state: []const f64, layer_count: usize) !void {
    if (inputs.first_layer_fortran == 0 or inputs.last_layer_fortran < inputs.first_layer_fortran or
        inputs.last_layer_fortran > layer_count)
        return error.MicroporeBandPhosphorusLayerRangeMismatch;
    const expected = std.math.mul(usize, layer_count, species_count) catch
        return error.MicroporeBandPhosphorusDimensionMismatch;
    if (state.len != expected or inputs.net_micropore_band_flux_mol_per_step.len != expected or
        inputs.macropore_exchange_flux_mol_per_step.len != expected or
        inputs.band_subsurface_flux_mol_per_step.len != expected)
        return error.MicroporeBandPhosphorusDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_layer_thickness_m))
        return error.NonFiniteMicroporeBandPhosphorusInput;
    const slices = [_][]const f64{
        inputs.layer_thickness_m,
        inputs.net_micropore_band_flux_mol_per_step,
        inputs.macropore_exchange_flux_mol_per_step,
        inputs.band_subsurface_flux_mol_per_step,
        state,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteMicroporeBandPhosphorusInput;
}

fn inputsFor(thickness: []const f64, net: []const f64, exchange: []const f64, subsurface: []const f64) Inputs {
    return .{
        .salt_simulation = .enabled,
        .first_layer_fortran = 1,
        .last_layer_fortran = thickness.len,
        .minimum_layer_thickness_m = 0.1,
        .layer_thickness_m = thickness,
        .net_micropore_band_flux_mol_per_step = net,
        .macropore_exchange_flux_mol_per_step = exchange,
        .band_subsurface_flux_mol_per_step = subsurface,
    };
}

test "disabled salt simulation bypasses dormant band phosphorus storage" {
    const empty = [_]f64{};
    var state = [_]f64{};
    const inputs = Inputs{
        .salt_simulation = .disabled,
        .first_layer_fortran = 0,
        .last_layer_fortran = std.math.maxInt(usize),
        .minimum_layer_thickness_m = std.math.nan(f64),
        .layer_thickness_m = &empty,
        .net_micropore_band_flux_mol_per_step = &empty,
        .macropore_exchange_flux_mol_per_step = &empty,
        .band_subsurface_flux_mol_per_step = &empty,
    };
    try std.testing.expectEqual(@as(usize, 0), try apply(inputs, &state));
}

test "TRNSFRS updates all eight banded micropore phosphorus species" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{2} ** species_count;
    const exchange = [_]f64{3} ** species_count;
    const subsurface = [_]f64{4} ** species_count;
    var state = [_]f64{10} ** species_count;
    try std.testing.expectEqual(@as(usize, 1), try apply(inputsFor(&thickness, &net, &exchange, &subsurface), &state));
    try std.testing.expectEqualSlices(f64, &([_]f64{19} ** species_count), &state);
}

test "NU NL range and strict layer thickness guard are preserved" {
    const thickness = [_]f64{ 0.2, 0.1, 0.3 };
    const net = [_]f64{2} ** (3 * species_count);
    const exchange = [_]f64{3} ** (3 * species_count);
    const subsurface = [_]f64{4} ** (3 * species_count);
    var state = [_]f64{10} ** (3 * species_count);
    var inputs = inputsFor(&thickness, &net, &exchange, &subsurface);
    inputs.first_layer_fortran = 2;
    try std.testing.expectEqual(@as(usize, 1), try apply(inputs, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
    try std.testing.expectEqual(@as(f64, 10), state[species_count]);
    try std.testing.expectEqual(@as(f64, 19), state[2 * species_count]);
}

test "invalid band layer range fails atomically" {
    const thickness = [_]f64{0.2};
    const net = [_]f64{2} ** species_count;
    const exchange = [_]f64{3} ** species_count;
    const subsurface = [_]f64{4} ** species_count;
    var state = [_]f64{10} ** species_count;
    var inputs = inputsFor(&thickness, &net, &exchange, &subsurface);
    inputs.last_layer_fortran = 2;
    try std.testing.expectError(error.MicroporeBandPhosphorusLayerRangeMismatch, apply(inputs, &state));
    try std.testing.expectEqual(@as(f64, 10), state[0]);
}

test "runtime band topology mismatch fails atomically" {
    const thickness = [_]f64{0.2};
    const short_net = [_]f64{2} ** (species_count - 1);
    const exchange = [_]f64{3} ** species_count;
    const subsurface = [_]f64{4} ** species_count;
    var state = [_]f64{10} ** species_count;
    try std.testing.expectError(error.MicroporeBandPhosphorusDimensionMismatch, apply(inputsFor(&thickness, &short_net, &exchange, &subsurface), &state));
}

test "late band overflow leaves prior layers atomic" {
    const thickness = [_]f64{ 0.2, 0.2 };
    var net = [_]f64{0} ** (2 * species_count);
    net[net.len - 1] = std.math.floatMax(f64);
    const exchange = [_]f64{0} ** (2 * species_count);
    const subsurface = [_]f64{0} ** (2 * species_count);
    var state = [_]f64{0} ** (2 * species_count);
    state[state.len - 1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteMicroporeBandPhosphorusResult, apply(inputsFor(&thickness, &net, &exchange, &subsurface), &state));
    try std.testing.expectEqual(@as(f64, 0), state[0]);
}
