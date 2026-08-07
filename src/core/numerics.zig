const std = @import("std");

/// In-place partial-pivot Gaussian solve for small runtime Newton systems.
/// False indicates a singular or non-finite Jacobian so callers can fall back
/// to Picard without committing a candidate.
pub fn solveDenseLinearSystem(matrix: []f64, right_hand_side: []f64, dimension: usize) bool {
    if (dimension == 0 or matrix.len != dimension * dimension or right_hand_side.len != dimension) return false;
    for (0..dimension) |pivot_column| {
        var pivot_row = pivot_column;
        var pivot_magnitude = @abs(matrix[pivot_column * dimension + pivot_column]);
        for (pivot_column + 1..dimension) |row| {
            const magnitude = @abs(matrix[row * dimension + pivot_column]);
            if (magnitude > pivot_magnitude) {
                pivot_magnitude = magnitude;
                pivot_row = row;
            }
        }
        if (!std.math.isFinite(pivot_magnitude) or pivot_magnitude <= 1e-14) return false;
        if (pivot_row != pivot_column) {
            for (0..dimension) |column| std.mem.swap(f64, &matrix[pivot_column * dimension + column], &matrix[pivot_row * dimension + column]);
            std.mem.swap(f64, &right_hand_side[pivot_column], &right_hand_side[pivot_row]);
        }
        const pivot = matrix[pivot_column * dimension + pivot_column];
        for (pivot_column + 1..dimension) |row| {
            const factor = matrix[row * dimension + pivot_column] / pivot;
            if (!std.math.isFinite(factor)) return false;
            matrix[row * dimension + pivot_column] = 0;
            for (pivot_column + 1..dimension) |column| matrix[row * dimension + column] -= factor * matrix[pivot_column * dimension + column];
            right_hand_side[row] -= factor * right_hand_side[pivot_column];
        }
    }
    var row = dimension;
    while (row > 0) {
        row -= 1;
        var value = right_hand_side[row];
        for (row + 1..dimension) |column| value -= matrix[row * dimension + column] * right_hand_side[column];
        const diagonal = matrix[row * dimension + row];
        if (!std.math.isFinite(value) or !std.math.isFinite(diagonal) or @abs(diagonal) <= 1e-14) return false;
        right_hand_side[row] = value / diagonal;
        if (!std.math.isFinite(right_hand_side[row])) return false;
    }
    return true;
}

pub const SolverOptions = struct {
    absolute_tolerance: f64 = 1.0e-11,
    relative_tolerance: f64 = 1.0e-8,
    derivative_floor: f64 = 1.0e-14,
    picard_relaxation: f64 = 0.5,
    residual_scale: f64 = 1.0,
    max_iterations: u16 = 40,
    safeguard_with_bracket: bool = false,
};

pub const SolveResult = struct {
    root: f64,
    residual: f64,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

/// Newton-Raphson/Picard hybrid. A Newton-Raphson candidate is accepted only
/// when it remains inside the physical bounds and reduces the residual.
/// Otherwise the process-specific Picard fixed-point map is under-relaxed.
/// Callers with a continuous residual kink may optionally retain a sign
/// bracket and choose its midpoint when that is a safer Picard candidate.
pub fn newtonPicard(
    context: anytype,
    residualFn: anytype,
    derivativeFn: anytype,
    picardFn: anytype,
    lower_bound: f64,
    upper_bound: f64,
    initial_guess: f64,
    options: SolverOptions,
) !SolveResult {
    try validateOptions(options);
    if (!std.math.isFinite(lower_bound) or !std.math.isFinite(upper_bound) or lower_bound >= upper_bound) return error.InvalidBounds;
    var x = std.math.clamp(initial_guess, lower_bound, upper_bound);
    var iteration: u16 = 0;
    var newton_raphson_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var bracket_lower = lower_bound;
    var bracket_upper = upper_bound;
    var bracket_lower_residual: f64 = undefined;
    var bracket_upper_residual: f64 = undefined;
    var bracket_valid = false;
    if (options.safeguard_with_bracket) {
        bracket_lower_residual = residualFn(context, bracket_lower);
        bracket_upper_residual = residualFn(context, bracket_upper);
        try requireFinite("Newton-Picard lower-bound residual", bracket_lower_residual);
        try requireFinite("Newton-Picard upper-bound residual", bracket_upper_residual);
        bracket_valid =
            std.math.signbit(bracket_lower_residual) !=
            std.math.signbit(bracket_upper_residual);
    }
    while (iteration < options.max_iterations) : (iteration += 1) {
        const current_residual = residualFn(context, x);
        try requireFinite("Newton-Picard residual", current_residual);
        const residual_tolerance = options.absolute_tolerance + options.relative_tolerance * options.residual_scale;
        if (@abs(current_residual) <= residual_tolerance) return .{
            .root = x,
            .residual = current_residual,
            .iterations = iteration + 1,
            .newton_raphson_steps = newton_raphson_steps,
            .picard_steps = picard_steps,
        };
        if (bracket_valid) {
            if (std.math.signbit(current_residual) ==
                std.math.signbit(bracket_lower_residual))
            {
                bracket_lower = x;
                bracket_lower_residual = current_residual;
            } else {
                bracket_upper = x;
                bracket_upper_residual = current_residual;
            }
        }

        const derivative_value = derivativeFn(context, x);
        try requireFinite("Newton-Raphson derivative", derivative_value);
        if (@abs(derivative_value) > options.derivative_floor) {
            const candidate = x - current_residual / derivative_value;
            if (std.math.isFinite(candidate) and
                candidate >= (if (bracket_valid) bracket_lower else lower_bound) and
                candidate <= (if (bracket_valid) bracket_upper else upper_bound))
            {
                const candidate_residual = residualFn(context, candidate);
                try requireFinite("Newton-Raphson candidate residual", candidate_residual);
                if (@abs(candidate_residual) < @abs(current_residual)) {
                    x = candidate;
                    newton_raphson_steps += 1;
                    continue;
                }
            }
        }

        const fixed_point = picardFn(context, x);
        try requireFinite("Picard fixed point", fixed_point);
        const relaxed = x + options.picard_relaxation * (fixed_point - x);
        if (!std.math.isFinite(relaxed)) return error.NonFinitePicardIterate;
        var next_x = std.math.clamp(
            relaxed,
            if (bracket_valid) bracket_lower else lower_bound,
            if (bracket_valid) bracket_upper else upper_bound,
        );
        if (bracket_valid) {
            const relaxed_residual = residualFn(context, next_x);
            try requireFinite("Picard candidate residual", relaxed_residual);
            const midpoint = 0.5 * (bracket_lower + bracket_upper);
            const midpoint_residual = residualFn(context, midpoint);
            try requireFinite("bracketed Picard residual", midpoint_residual);
            if (@abs(midpoint_residual) < @abs(relaxed_residual))
                next_x = midpoint;
        }
        // Scale stagnation to the solve interval. Using a unit floor makes
        // every physically meaningful Picard step look stationary when a
        // mass inventory is much smaller than one model unit.
        const step_tolerance = 8.0 * std.math.floatEps(f64) * @max(@abs(x), upper_bound - lower_bound);
        if (@abs(next_x - x) <= step_tolerance) {
            std.log.err("Newton-Picard stagnated: iteration={d} x={e} residual={e}", .{ iteration + 1, x, current_residual });
            return error.NewtonPicardStagnated;
        }
        x = next_x;
        picard_steps += 1;
    }
    const final_residual = residualFn(context, x);
    try requireFinite("final Newton-Picard residual", final_residual);
    const residual_tolerance =
        options.absolute_tolerance +
        options.relative_tolerance * options.residual_scale;
    if (@abs(final_residual) <= residual_tolerance) return .{
        .root = x,
        .residual = final_residual,
        .iterations = options.max_iterations,
        .newton_raphson_steps = newton_raphson_steps,
        .picard_steps = picard_steps,
    };
    std.log.warn("Newton-Picard failed: iterations={d} bounds=[{e},{e}] last_x={e} residual={e}", .{ options.max_iterations, lower_bound, upper_bound, x, final_residual });
    return error.NewtonPicardDidNotConverge;
}

pub fn newtonPicardFiniteDifference(context: anytype, residualFn: anytype, picardFn: anytype, lower_bound: f64, upper_bound: f64, initial_guess: f64, options: SolverOptions) !SolveResult {
    const Adapter = struct {
        fn derivative(data: @TypeOf(context), x: f64) f64 {
            const scale = @max(@abs(x), 1.0);
            const step = @sqrt(std.math.floatEps(f64)) * scale;
            return (residualFn(data, x + step) - residualFn(data, x - step)) / (2.0 * step);
        }
    };
    return newtonPicard(context, residualFn, Adapter.derivative, picardFn, lower_bound, upper_bound, initial_guess, options);
}

pub fn requireFinite(comptime label: []const u8, value: f64) !void {
    if (!std.math.isFinite(value)) {
        std.log.err("non-finite numeric value: {s}={e}", .{ label, value });
        return error.NonFiniteNumericValue;
    }
}

fn validateOptions(options: SolverOptions) !void {
    if (!std.math.isFinite(options.absolute_tolerance) or options.absolute_tolerance <= 0 or
        !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or
        !std.math.isFinite(options.derivative_floor) or options.derivative_floor <= 0 or
        !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or
        !std.math.isFinite(options.residual_scale) or options.residual_scale <= 0 or
        options.max_iterations == 0) return error.InvalidSolverOptions;
}

fn squareMinusTwo(_: void, x: f64) f64 {
    return x * x - 2.0;
}
fn squareMinusTwoDerivative(_: void, x: f64) f64 {
    return 2.0 * x;
}
fn squareRootPicard(_: void, x: f64) f64 {
    return 2.0 / x;
}
fn unusableDerivative(_: void, _: f64) f64 {
    return 0;
}
fn tinyResidual(_: void, x: f64) f64 {
    return x - 2.0e-10;
}
fn tinyPicard(_: void, _: f64) f64 {
    return 2.0e-10;
}

test "Newton-Raphson path converges" {
    const result = try newtonPicard({}, squareMinusTwo, squareMinusTwoDerivative, squareRootPicard, 0.1, 2.0, 1.8, .{});
    try std.testing.expectApproxEqRel(@sqrt(2.0), result.root, 1.0e-10);
    try std.testing.expect(result.newton_raphson_steps > 0);
    try std.testing.expect(result.iterations < 40);
}

test "Newton-Picard exits on an already converged initial state" {
    const root = @sqrt(2.0);
    const result = try newtonPicard({}, squareMinusTwo, squareMinusTwoDerivative, squareRootPicard, 0.1, 2.0, root, .{ .max_iterations = 80 });
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(u16, 0), result.newton_raphson_steps);
    try std.testing.expectEqual(@as(u16, 0), result.picard_steps);
}

test "Picard path converges when Newton derivative is unavailable" {
    const result = try newtonPicard({}, squareMinusTwo, unusableDerivative, squareRootPicard, 0.1, 2.0, 1.8, .{ .picard_relaxation = 0.5 });
    try std.testing.expectApproxEqRel(@sqrt(2.0), result.root, 1.0e-8);
    try std.testing.expect(result.picard_steps > 0);
    try std.testing.expectEqual(@as(u16, 0), result.newton_raphson_steps);
}

test "Newton-Picard audits the iterate produced by its final allowed update" {
    const result = try newtonPicard(
        {},
        tinyResidual,
        unusableDerivative,
        tinyPicard,
        0,
        1.0e-9,
        9.0e-10,
        .{
            .absolute_tolerance = 1.0e-20,
            .relative_tolerance = 1.0e-10,
            .residual_scale = 1.0e-9,
            .picard_relaxation = 1,
            .max_iterations = 1,
        },
    );
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(u16, 1), result.picard_steps);
    try std.testing.expectApproxEqAbs(
        @as(f64, 2.0e-10),
        result.root,
        1.0e-25,
    );
}

test "Picard stagnation threshold follows a tiny physical interval" {
    const result = try newtonPicard({}, tinyResidual, unusableDerivative, tinyPicard, 0, 1.0e-9, 9.0e-10, .{
        .absolute_tolerance = 1.0e-20,
        .relative_tolerance = 1.0e-10,
        .residual_scale = 1.0e-9,
        .picard_relaxation = 0.5,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 2.0e-10), result.root, 2.0e-18);
    try std.testing.expect(result.picard_steps > 0);
}

test "finite-difference Newton-Picard converges" {
    const result = try newtonPicardFiniteDifference({}, squareMinusTwo, squareRootPicard, 0.1, 2.0, 0.7, .{});
    try std.testing.expectApproxEqRel(@sqrt(2.0), result.root, 1.0e-10);
}

test "Newton-Picard rejects unusable options" {
    try std.testing.expectError(error.InvalidSolverOptions, newtonPicard({}, squareMinusTwo, squareMinusTwoDerivative, squareRootPicard, 0.1, 2.0, 1.0, .{ .picard_relaxation = 0 }));
}
