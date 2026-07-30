const std = @import("std");
const Face = @import("solute_transport.zig").Face;

pub const Options = struct {
    absolute_tolerance: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    max_iterations: u16,
    maximum_convective_fraction: f64 = 1,
    pore_exchange_fraction: f64 = 1,
};

pub const Inputs = struct {
    species_count: usize,
    faces: []const Face,
    micropore_conductance_m3_per_step: []const f64,
    macropore_conductance_m3_per_step: []const f64,
    micropore_water_m3: []const f64,
    macropore_water_m3: []const f64,
    layer_bulk_volume_m3: []const f64,
    micropore_external_water_flux_m3_per_step: []const f64,
    macropore_external_water_flux_m3_per_step: []const f64,
    recharge_concentration_per_m3: []const f64,
};

pub const Result = struct {
    micropore_iterations: u16,
    macropore_iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

/// Runtime-species extensive aqueous transport. Units may be grams or moles,
/// but every amount, concentration, tolerance, and boundary ledger supplied
/// to one invocation must use the same explicit extensive unit.
pub fn advance(
    allocator: std.mem.Allocator,
    micropore_amount: []f64,
    macropore_amount: []f64,
    boundary_net_flux: []f64,
    inputs: Inputs,
    options: Options,
) !Result {
    try validate(micropore_amount, macropore_amount, boundary_net_flux, inputs, options);
    const micro_before = try allocator.dupe(f64, micropore_amount);
    defer allocator.free(micro_before);
    const macro_before = try allocator.dupe(f64, macropore_amount);
    defer allocator.free(macro_before);
    const boundary_before = try allocator.dupe(f64, boundary_net_flux);
    defer allocator.free(boundary_before);
    var committed = false;
    defer if (!committed) {
        @memcpy(micropore_amount, micro_before);
        @memcpy(macropore_amount, macro_before);
        @memcpy(boundary_net_flux, boundary_before);
    };
    const micro_result = try solve(allocator, micropore_amount, inputs.micropore_water_m3, inputs.faces, inputs.micropore_conductance_m3_per_step, inputs.species_count, options);
    const macro_result = try solve(allocator, macropore_amount, inputs.macropore_water_m3, inputs.faces, inputs.macropore_conductance_m3_per_step, inputs.species_count, options);
    @memset(boundary_net_flux, 0);
    for (0..inputs.micropore_water_m3.len) |layer| {
        const base = layer * inputs.species_count;
        try commitBoundary(micropore_amount[base..][0..inputs.species_count], inputs.micropore_water_m3[layer], inputs.micropore_external_water_flux_m3_per_step[layer], inputs.recharge_concentration_per_m3[base..][0..inputs.species_count], options.maximum_convective_fraction, boundary_net_flux[base..][0..inputs.species_count]);
        try commitBoundary(macropore_amount[base..][0..inputs.species_count], inputs.macropore_water_m3[layer], inputs.macropore_external_water_flux_m3_per_step[layer], inputs.recharge_concentration_per_m3[base..][0..inputs.species_count], options.maximum_convective_fraction, boundary_net_flux[base..][0..inputs.species_count]);
        for (0..inputs.species_count) |species| {
            const index = base + species;
            const exchange = try poreExchange(micropore_amount[index], macropore_amount[index], inputs.micropore_water_m3[layer], inputs.macropore_water_m3[layer], inputs.layer_bulk_volume_m3[layer], options.pore_exchange_fraction);
            micropore_amount[index] += exchange;
            macropore_amount[index] -= exchange;
        }
    }
    committed = true;
    return .{
        .micropore_iterations = micro_result.iterations,
        .macropore_iterations = macro_result.iterations,
        .newton_raphson_steps = micro_result.newton_steps + macro_result.newton_steps,
        .picard_steps = micro_result.picard_steps + macro_result.picard_steps,
    };
}

const SolverResult = struct { iterations: u16, newton_steps: u16, picard_steps: u16 };

fn solve(allocator: std.mem.Allocator, amounts: []f64, water: []const f64, faces: []const Face, conductance: []const f64, species_count: usize, options: Options) !SolverResult {
    const base = try allocator.dupe(f64, amounts);
    defer allocator.free(base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, amounts.len);
    defer allocator.free(residual);
    const probe = try allocator.alloc(f64, amounts.len);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, amounts.len);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, amounts.len);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, amounts.len);
    defer allocator.free(candidate_residual);
    const target = try allocator.alloc(f64, amounts.len);
    defer allocator.free(target);
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(base, current, water, faces, conductance, species_count, options.maximum_convective_fraction, target, residual);
        const norm = try scaledNorm(current, residual, options);
        if (norm <= 1) {
            @memcpy(amounts, target);
            return .{ .iterations = iteration + 1, .newton_steps = newton_steps, .picard_steps = picard_steps };
        }
        var accepted_newton = false;
        try addDirection(current, residual, 0.5, probe);
        try residualAt(base, probe, water, faces, conductance, species_count, options.maximum_convective_fraction, target, probe_residual);
        var numerator: f64 = 0;
        var denominator: f64 = 0;
        for (residual, probe_residual) |value, sampled| {
            const derivative = (sampled - value) / 0.5;
            numerator += value * derivative;
            denominator += derivative * derivative;
        }
        if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
            const fraction = std.math.clamp(-numerator / denominator, 0.05, 1.5);
            if (addDirection(current, residual, fraction, candidate)) |_| {
                if (residualAt(base, candidate, water, faces, conductance, species_count, options.maximum_convective_fraction, target, candidate_residual)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        accepted_newton = true;
                    }
                } else |_| {}
            } else |_| {}
        }
        if (accepted_newton) continue;
        try addDirection(current, residual, options.picard_relaxation, candidate);
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    return error.AqueousExtensiveTransportDidNotConverge;
}

fn residualAt(base: []const f64, trial: []const f64, water: []const f64, faces: []const Face, conductance: []const f64, species_count: usize, maximum_fraction: f64, target: []f64, residual: []f64) !void {
    @memcpy(target, base);
    for (faces, 0..) |face, face_index| for (0..species_count) |species| {
        const first = face.first_cell * species_count + species;
        const second = face.second_cell * species_count + species;
        const first_concentration = if (water[face.first_cell] > 0) trial[first] / water[face.first_cell] else 0;
        const second_concentration = if (water[face.second_cell] > 0) trial[second] / water[face.second_cell] else 0;
        const donor_fraction = if (face.water_flux_m3_per_step >= 0)
            if (water[face.first_cell] > 0) @min(maximum_fraction, face.water_flux_m3_per_step / water[face.first_cell]) else maximum_fraction
        else if (water[face.second_cell] > 0) @min(maximum_fraction, -face.water_flux_m3_per_step / water[face.second_cell]) else maximum_fraction;
        const convection = if (face.water_flux_m3_per_step >= 0) donor_fraction * trial[first] else -donor_fraction * trial[second];
        const diffusion = conductance[face_index * species_count + species] * (first_concentration - second_concentration);
        const flux = std.math.clamp(convection + diffusion, -target[second], target[first]);
        target[first] -= flux;
        target[second] += flux;
    };
    for (target, trial, residual) |fixed_point, value, *difference| {
        difference.* = fixed_point - value;
        if (!std.math.isFinite(difference.*)) return error.NonFiniteAqueousExtensiveTransport;
    }
}

fn commitBoundary(amounts: []f64, water_m3: f64, outward_water_m3: f64, recharge: []const f64, maximum_fraction: f64, ledger: []f64) !void {
    for (amounts, recharge, ledger) |*amount, concentration, *net| {
        const change = if (outward_water_m3 >= 0)
            -amount.* * (if (water_m3 > 0) @min(maximum_fraction, outward_water_m3 / water_m3) else maximum_fraction)
        else
            -outward_water_m3 * concentration;
        if (!std.math.isFinite(change) or amount.* + change < -1e-12 or !std.math.isFinite(net.* + change)) return error.InvalidAqueousExtensiveBoundaryFlux;
        amount.* = @max(0, amount.* + change);
        net.* += change;
    }
}

fn poreExchange(micro: f64, macro: f64, micro_water: f64, macro_water: f64, bulk_volume: f64, fraction: f64) !f64 {
    if (macro_water == 0) return 0;
    const exchanging_macro_water = @min(0.05 * bulk_volume, macro_water);
    const combined_water = micro_water + exchanging_macro_water;
    if (combined_water == 0) return 0;
    const exchange = fraction * (macro * micro_water - micro * exchanging_macro_water) / combined_water;
    if (!std.math.isFinite(exchange)) return error.NonFiniteAqueousExtensivePoreExchange;
    return std.math.clamp(exchange, -micro, macro);
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteAqueousExtensiveTransport;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance + options.relative_tolerance * @max(1, value)));
    }
    return maximum;
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, change, *candidate| {
        candidate.* = value + fraction * change;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidAqueousExtensiveTransportCandidate;
        candidate.* = @max(0, candidate.*);
    }
}

fn validate(micro: []const f64, macro: []const f64, boundary: []const f64, inputs: Inputs, options: Options) !void {
    if (inputs.species_count == 0 or micro.len != macro.len or micro.len != boundary.len or micro.len != inputs.micropore_water_m3.len * inputs.species_count or inputs.macropore_water_m3.len != inputs.micropore_water_m3.len or inputs.layer_bulk_volume_m3.len != inputs.micropore_water_m3.len or inputs.micropore_external_water_flux_m3_per_step.len != inputs.micropore_water_m3.len or inputs.macropore_external_water_flux_m3_per_step.len != inputs.micropore_water_m3.len or inputs.recharge_concentration_per_m3.len != micro.len or inputs.micropore_conductance_m3_per_step.len != inputs.faces.len * inputs.species_count or inputs.macropore_conductance_m3_per_step.len != inputs.faces.len * inputs.species_count) return error.AqueousExtensiveTransportDimensionMismatch;
    if (!std.math.isFinite(options.absolute_tolerance) or options.absolute_tolerance <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or options.max_iterations == 0 or !std.math.isFinite(options.maximum_convective_fraction) or options.maximum_convective_fraction < 0 or options.maximum_convective_fraction > 1 or !std.math.isFinite(options.pore_exchange_fraction) or options.pore_exchange_fraction < 0 or options.pore_exchange_fraction > 1) return error.InvalidAqueousExtensiveTransportOptions;
    for (micro, macro, inputs.recharge_concentration_per_m3) |a, b, recharge| if (!std.math.isFinite(a) or a < 0 or !std.math.isFinite(b) or b < 0 or !std.math.isFinite(recharge) or recharge < 0) return error.InvalidAqueousExtensiveTransportState;
    for (inputs.micropore_water_m3, inputs.macropore_water_m3, inputs.layer_bulk_volume_m3, inputs.micropore_external_water_flux_m3_per_step, inputs.macropore_external_water_flux_m3_per_step) |a, b, bulk, fa, fb| {
        if (!std.math.isFinite(a) or a < 0 or !std.math.isFinite(b) or b < 0 or !std.math.isFinite(bulk) or bulk <= 0 or !std.math.isFinite(fa) or !std.math.isFinite(fb)) return error.InvalidAqueousExtensiveTransportState;
    }
}

test "runtime species transport conserves internal mass and publishes boundary sign" {
    var micro = [_]f64{ 2, 4 };
    var macro = [_]f64{ 0, 0 };
    var boundary = [_]f64{ 0, 0 };
    const result = try advance(std.testing.allocator, &micro, &macro, &boundary, .{
        .species_count = 2,
        .faces = &.{},
        .micropore_conductance_m3_per_step = &.{},
        .macropore_conductance_m3_per_step = &.{},
        .micropore_water_m3 = &.{1},
        .macropore_water_m3 = &.{0},
        .layer_bulk_volume_m3 = &.{1},
        .micropore_external_water_flux_m3_per_step = &.{0.25},
        .macropore_external_water_flux_m3_per_step = &.{0},
        .recharge_concentration_per_m3 = &.{ 0, 0 },
    }, .{ .absolute_tolerance = 1e-12, .relative_tolerance = 1e-8, .picard_relaxation = 0.5, .max_iterations = 20, .pore_exchange_fraction = 0 });
    try std.testing.expectEqual(@as(u16, 1), result.micropore_iterations);
    try std.testing.expectEqualSlices(f64, &.{ 1.5, 3 }, &micro);
    try std.testing.expectEqualSlices(f64, &.{ -0.5, -1 }, &boundary);
}
