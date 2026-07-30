const std = @import("std");

pub const Options = struct {
    absolute_tolerance: f64 = 1e-11,
    relative_tolerance: f64 = 1e-8,
    derivative_floor: f64 = 1e-14,
    finite_difference_fraction: f64 = 1e-7,
    picard_relaxation: f64 = 0.5,
    /// `MRXN=60` in SOLUTE.F. This is a ceiling, not a fixed cycle count.
    max_iterations: u16 = 60,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Transactional multi-species hybrid for a single chemistry cell/zone.
/// `residualFn(context, state, residual)` and
/// `picardFn(context, state, fixed_point)` must fill caller-owned slices.
/// The supplied state is changed only after convergence.
pub fn solve(allocator: std.mem.Allocator, context: anytype, state: []f64, residualFn: anytype, picardFn: anytype, options: Options) !Result {
    try validateOptions(options, state);
    const count = state.len;
    const current = try allocator.dupe(f64, state);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, count);
    defer allocator.free(residual);
    const trial = try allocator.alloc(f64, count);
    defer allocator.free(trial);
    const newton_candidate = try allocator.alloc(f64, count);
    defer allocator.free(newton_candidate);
    const trial_residual = try allocator.alloc(f64, count);
    defer allocator.free(trial_residual);
    const fixed_point = try allocator.alloc(f64, count);
    defer allocator.free(fixed_point);

    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualFn(context, current, residual);
        const current_norm = try scaledNorm(current, residual, options);
        if (current_norm <= 1) {
            @memcpy(state, current);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = current_norm };
        }

        @memcpy(newton_candidate, current);
        var usable_newton = true;
        for (0..count) |index| {
            @memcpy(trial, current);
            const original = current[index];
            const difference = options.finite_difference_fraction * @max(1.0, @abs(original));
            trial[index] = original + difference;
            try residualFn(context, trial, trial_residual);
            const derivative = (trial_residual[index] - residual[index]) / difference;
            if (!std.math.isFinite(derivative) or @abs(derivative) <= options.derivative_floor) {
                usable_newton = false;
                break;
            }
            const candidate = original - residual[index] / derivative;
            if (!std.math.isFinite(candidate) or candidate < 0) {
                usable_newton = false;
                break;
            }
            newton_candidate[index] = candidate;
        }
        if (usable_newton) {
            try residualFn(context, newton_candidate, trial_residual);
            const trial_norm = try scaledNorm(newton_candidate, trial_residual, options);
            if (trial_norm < current_norm) {
                @memcpy(current, newton_candidate);
                newton_steps += 1;
                continue;
            }
        }

        try picardFn(context, current, fixed_point);
        var maximum_change: f64 = 0;
        for (current, fixed_point) |*value, target| {
            if (!std.math.isFinite(target) or target < 0) return error.InvalidSolutePicardState;
            const next = value.* + options.picard_relaxation * (target - value.*);
            if (!std.math.isFinite(next) or next < 0) return error.InvalidSolutePicardState;
            maximum_change = @max(maximum_change, @abs(next - value.*));
            value.* = next;
        }
        if (maximum_change <= std.math.floatEps(f64) * @max(1.0, maximumMagnitude(current)))
            return error.SoluteEquilibriumStagnated;
        picard_steps += 1;
    }
    try residualFn(context, current, residual);
    const final_norm = try scaledNorm(current, residual, options);
    std.log.err("SOLUTE equilibrium failed: iterations={d} species={d} maximum_scaled_residual={e}", .{ options.max_iterations, count, final_norm });
    return error.SoluteEquilibriumDidNotConverge;
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    if (state.len != residual.len) return error.SoluteResidualSizeMismatch;
    var maximum: f64 = 0;
    for (state, residual) |value, mismatch| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(mismatch)) return error.NonFiniteSoluteEquilibriumState;
        const tolerance = options.absolute_tolerance + options.relative_tolerance * @max(1.0, @abs(value));
        maximum = @max(maximum, @abs(mismatch) / tolerance);
    }
    return maximum;
}

fn maximumMagnitude(values: []const f64) f64 {
    var maximum: f64 = 0;
    for (values) |value| maximum = @max(maximum, @abs(value));
    return maximum;
}

fn validateOptions(options: Options, state: []const f64) !void {
    if (state.len == 0 or options.max_iterations == 0 or
        !std.math.isFinite(options.absolute_tolerance) or options.absolute_tolerance <= 0 or
        !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or
        !std.math.isFinite(options.derivative_floor) or options.derivative_floor <= 0 or
        !std.math.isFinite(options.finite_difference_fraction) or options.finite_difference_fraction <= 0 or
        !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1)
        return error.InvalidSoluteEquilibriumOptions;
    for (state) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteSoluteEquilibriumState;
}

const TestContext = struct { targets: []const f64 };

fn testResidual(context: TestContext, state: []const f64, residual: []f64) !void {
    if (state.len != context.targets.len or residual.len != state.len) return error.TestSizeMismatch;
    for (state, context.targets, residual) |value, target, *mismatch| mismatch.* = value * value - target * target;
}

fn testPicard(context: TestContext, _: []const f64, fixed_point: []f64) !void {
    @memcpy(fixed_point, context.targets);
}

test "SOLUTE hybrid uses the legacy 60 ceiling but converges early" {
    const allocator = std.testing.allocator;
    const targets = [_]f64{ 2, 3, 4 };
    var state = [_]f64{ 0.5, 0.8, 1.2 };
    const result = try solve(allocator, TestContext{ .targets = &targets }, &state, testResidual, testPicard, .{});
    try std.testing.expect(result.iterations < 60);
    try std.testing.expect(result.newton_raphson_steps > 0);
    for (state, targets) |value, target| try std.testing.expectApproxEqRel(target, value, 1e-8);
}

fn impossibleResidual(_: void, _: []const f64, residual: []f64) !void {
    @memset(residual, 1);
}

fn unchangedPicard(_: void, state: []const f64, fixed_point: []f64) !void {
    @memcpy(fixed_point, state);
}

test "failed SOLUTE solve leaves caller state unchanged" {
    var state = [_]f64{1};
    try std.testing.expectError(error.SoluteEquilibriumStagnated, solve(std.testing.allocator, {}, &state, impossibleResidual, unchangedPicard, .{ .max_iterations = 2 }));
    try std.testing.expectEqual(@as(f64, 1), state[0]);
}
