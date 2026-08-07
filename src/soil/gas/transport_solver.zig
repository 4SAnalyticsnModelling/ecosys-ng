const std = @import("std");
const gas = @import("transport.zig");

pub const Options = struct {
    absolute_tolerance_g: f64 = 1e-12,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
    /// Runtime NPG input: a convergence ceiling, not a mandatory cycle count.
    max_iterations: u16,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Coupled implicit gas diffusion solve. Unlike the old NPG full-model loop,
/// this reevaluates only the nonlinear gas kernel and exits on convergence.
pub fn solve(allocator: std.mem.Allocator, state: *gas.State, faces: []const gas.Face, conductance_m3_per_step: []const f64, options: Options) !Result {
    try validate(state, faces, conductance_m3_per_step, options);
    const base = try allocator.dupe(f64, state.gaseous_mass_g);
    defer allocator.free(base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, base.len);
    defer allocator.free(residual);
    const probe = try allocator.alloc(f64, base.len);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, base.len);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, base.len);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, base.len);
    defer allocator.free(candidate_residual);
    const target = try allocator.alloc(f64, base.len);
    defer allocator.free(target);
    var flux: [gas.species_count]f64 = undefined;
    var scratch = try gas.State.init(allocator, state.cell_count);
    defer scratch.deinit();
    @memcpy(scratch.air_volume_m3, state.air_volume_m3);
    @memcpy(scratch.temperature_k, state.temperature_k);
    @memcpy(scratch.water_vapor_mol, state.water_vapor_mol);

    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(&scratch, base, current, faces, conductance_m3_per_step, &flux, target, residual);
        const current_norm = try scaledNorm(current, residual, options);
        if (current_norm <= 1) {
            @memcpy(state.gaseous_mass_g, current);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = current_norm };
        }

        var accepted_newton = false;
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(&scratch, base, probe, faces, conductance_m3_per_step, &flux, target, probe_residual)) |_| {
                var numerator: f64 = 0;
                var denominator: f64 = 0;
                for (residual, probe_residual) |r, sampled| {
                    const derivative = (sampled - r) / options.directional_probe_fraction;
                    numerator += r * derivative;
                    denominator += derivative * derivative;
                }
                if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                    const fraction = std.math.clamp(-numerator / denominator, options.minimum_newton_fraction, options.maximum_newton_fraction);
                    if (addDirection(current, residual, fraction, candidate)) |_| {
                        if (residualAt(&scratch, base, candidate, faces, conductance_m3_per_step, &flux, target, candidate_residual)) |_| {
                            if (try scaledNorm(candidate, candidate_residual, options) < current_norm) {
                                @memcpy(current, candidate);
                                newton_steps += 1;
                                accepted_newton = true;
                            }
                        } else |_| {}
                    } else |_| {}
                }
            } else |_| {}
        } else |_| {}
        if (accepted_newton) continue;
        try addDirection(current, residual, options.picard_relaxation, candidate);
        if (maximumDifference(current, candidate) <= std.math.floatEps(f64) * @max(1.0, maximumMagnitude(current))) return error.GasTransportSolverStagnated;
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    return error.GasTransportSolverDidNotConverge;
}

fn residualAt(scratch: *gas.State, base: []const f64, trial: []const f64, faces: []const gas.Face, conductance: []const f64, flux: []f64, target: []f64, residual: []f64) !void {
    @memcpy(scratch.gaseous_mass_g, trial);
    @memcpy(target, base);
    for (faces, 0..) |face, face_index| {
        const start = face_index * gas.species_count;
        try gas.calculateFaceDiffusiveFluxesG(scratch, face, conductance[start .. start + gas.species_count], flux);
        for (flux, 0..) |value, species| {
            target[face.first_cell * gas.species_count + species] -= value;
            target[face.second_cell * gas.species_count + species] += value;
        }
    }
    for (target, trial, residual) |fixed_point, value, *r| {
        if (!std.math.isFinite(fixed_point) or fixed_point < -1e-12) return error.InvalidImplicitGasCandidate;
        r.* = fixed_point - value;
    }
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidImplicitGasCandidate;
        candidate.* = @max(0, candidate.*);
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteImplicitGasState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn validate(state: *const gas.State, faces: []const gas.Face, conductance: []const f64, options: Options) !void {
    if (conductance.len != try std.math.mul(usize, faces.len, gas.species_count)) return error.GasFaceParameterSizeMismatch;
    if (!std.math.isFinite(options.absolute_tolerance_g) or options.absolute_tolerance_g <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or options.max_iterations == 0) return error.InvalidGasTransportSolverOptions;
    for (faces) |face| if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell) return error.InvalidGasTransportFace;
}

fn maximumDifference(a: []const f64, b: []const f64) f64 {
    var maximum: f64 = 0;
    for (a, b) |left, right| maximum = @max(maximum, @abs(left - right));
    return maximum;
}

fn maximumMagnitude(values: []const f64) f64 {
    var maximum: f64 = 0;
    for (values) |value| maximum = @max(maximum, @abs(value));
    return maximum;
}

test "implicit seven-gas solve exits before NPG and conserves mass" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.air_volume_m3[1] = 1;
    (try state.gaseousMasses(0))[0] = 2;
    const conductance = [_]f64{0.1} ** gas.species_count;
    const result = try solve(std.testing.allocator, &state, &[_]gas.Face{.{ .first_cell = 0, .second_cell = 1 }}, &conductance, .{ .max_iterations = 20 });
    try std.testing.expect(result.iterations < 20);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.gaseous_mass_g[0] + state.gaseous_mass_g[gas.species_count], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0 / 3.0), state.gaseous_mass_g[0] - state.gaseous_mass_g[gas.species_count], 1e-8);
}

test "failed gas solve is atomic" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.air_volume_m3[1] = 1;
    (try state.gaseousMasses(0))[0] = 2;
    const before = try std.testing.allocator.dupe(f64, state.gaseous_mass_g);
    defer std.testing.allocator.free(before);
    const conductance = [_]f64{0.1} ** gas.species_count;
    try std.testing.expectError(error.GasTransportSolverDidNotConverge, solve(std.testing.allocator, &state, &[_]gas.Face{.{ .first_cell = 0, .second_cell = 1 }}, &conductance, .{ .max_iterations = 1 }));
    try std.testing.expectEqualSlices(f64, before, state.gaseous_mass_g);
}
