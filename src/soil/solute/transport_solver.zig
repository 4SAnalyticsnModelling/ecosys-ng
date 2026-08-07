const std = @import("std");
const transport = @import("transport.zig");
const numerics = @import("../../core/numerics.zig");

pub const Options = struct {
    absolute_tolerance_mol: f64 = 1e-12,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
    /// Runtime NPH input: a ceiling, never a mandatory number of sweeps.
    max_iterations: u16,
    /// Optional exact accepted ledger indexed face × runtime species.
    /// Positive moves `first_cell -> second_cell`.
    face_flux_mol_by_component: ?[]f64 = null,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Solves `final = initial + transport(final)` over all supplied faces. Each
/// face contributes equal and opposite species changes, so every iterate and
/// accepted Newton direction remains globally conservative.
pub fn solve(
    allocator: std.mem.Allocator,
    state: *transport.State,
    faces: []const transport.Face,
    diffusive_conductance_m3_per_step: []const f64,
    mobility_fraction: []const f64,
    face_parameters: transport.FaceParameters,
    options: Options,
) !Result {
    try validateInputs(state, faces, diffusive_conductance_m3_per_step, mobility_fraction, options);
    const component_count = state.amount_mol.len;
    const base = try allocator.dupe(f64, state.amount_mol);
    defer allocator.free(base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, component_count);
    defer allocator.free(residual);
    const probe = try allocator.alloc(f64, component_count);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, component_count);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, component_count);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, component_count);
    defer allocator.free(candidate_residual);
    const dense_direction = try allocator.alloc(f64, component_count);
    defer allocator.free(dense_direction);
    const dense_matrix = try allocator.alloc(
        f64,
        try std.math.mul(usize, state.cell_count, state.cell_count),
    );
    defer allocator.free(dense_matrix);
    const dense_rhs = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(dense_rhs);
    const fixed_point = try allocator.alloc(f64, component_count);
    defer allocator.free(fixed_point);
    const face_flux = try allocator.alloc(f64, state.species_count);
    defer allocator.free(face_flux);
    var scratch = try transport.State.init(allocator, state.cell_count, state.species_count);
    defer scratch.deinit();
    @memcpy(scratch.water_volume_m3, state.water_volume_m3);

    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(&scratch, base, current, faces, diffusive_conductance_m3_per_step, mobility_fraction, face_parameters, face_flux, fixed_point, residual);
        const current_norm = try scaledNorm(current, residual, options);
        if (current_norm <= 1) {
            if (options.face_flux_mol_by_component) |accepted_flux| {
                try captureAcceptedFaceFlux(
                    &scratch,
                    base,
                    faces,
                    diffusive_conductance_m3_per_step,
                    mobility_fraction,
                    face_parameters,
                    face_flux,
                    residual,
                    accepted_flux,
                );
            }
            @memcpy(state.amount_mol, fixed_point);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = current_norm };
        }

        var accepted_newton = false;
        // Species do not couple in the transport equations. Build every
        // unconverged species block in one Newton iteration; solving only the
        // current worst carrier made the NPH ceiling scale with the runtime
        // carrier count and could never cover all carriers when species>NPH.
        @memset(dense_direction, 0);
        var has_dense_direction = false;
        for (0..state.species_count) |species| {
            var species_norm: f64 = 0;
            for (0..state.cell_count) |cell| {
                const component = cell * state.species_count + species;
                species_norm = @max(
                    species_norm,
                    @abs(residual[component]) /
                        (options.absolute_tolerance_mol +
                            options.relative_tolerance *
                                @max(1.0, @abs(current[component]))),
                );
            }
            if (species_norm <= 1) continue;
            if (try denseSpeciesNewtonDirection(
                &scratch,
                base,
                current,
                residual,
                faces,
                diffusive_conductance_m3_per_step,
                mobility_fraction,
                face_parameters,
                options,
                species,
                face_flux,
                fixed_point,
                probe,
                probe_residual,
                dense_matrix,
                dense_rhs,
                dense_direction,
            )) has_dense_direction = true;
        }
        if (has_dense_direction) {
            var fraction: f64 = 1;
            var search: u8 = 0;
            while (search < 20) : (search += 1) {
                if (addDirection(
                    current,
                    dense_direction,
                    fraction,
                    candidate,
                )) |_| {
                    if (residualAt(
                        &scratch,
                        base,
                        candidate,
                        faces,
                        diffusive_conductance_m3_per_step,
                        mobility_fraction,
                        face_parameters,
                        face_flux,
                        fixed_point,
                        candidate_residual,
                    )) |_| {
                        if (try scaledNorm(
                            candidate,
                            candidate_residual,
                            options,
                        ) < current_norm) {
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                            break;
                        }
                    } else |_| {}
                } else |_| {}
                fraction *= 0.5;
            }
        }
        if (accepted_newton) continue;
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(&scratch, base, probe, faces, diffusive_conductance_m3_per_step, mobility_fraction, face_parameters, face_flux, fixed_point, probe_residual)) |_| {
                var numerator: f64 = 0;
                var denominator: f64 = 0;
                for (residual, probe_residual) |base_residual, sampled_residual| {
                    const derivative = (sampled_residual - base_residual) / options.directional_probe_fraction;
                    numerator += base_residual * derivative;
                    denominator += derivative * derivative;
                }
                if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                    const fraction = std.math.clamp(-numerator / denominator, options.minimum_newton_fraction, options.maximum_newton_fraction);
                    if (addDirection(current, residual, fraction, candidate)) |_| {
                        if (residualAt(&scratch, base, candidate, faces, diffusive_conductance_m3_per_step, mobility_fraction, face_parameters, face_flux, fixed_point, candidate_residual)) |_| {
                            const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                            if (candidate_norm < current_norm) {
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
        if (vectorsEqual(current, candidate)) {
            const stagnation_index =
                try worstResidualIndex(current, residual, options);
            std.log.err(
                "solute transport Newton-Picard stagnated: iteration={d} scaled_residual={e} cell={d} species={d} amount_mol={e} residual_mol={e}",
                .{
                    iteration + 1,
                    current_norm,
                    stagnation_index / state.species_count,
                    stagnation_index % state.species_count,
                    current[stagnation_index],
                    residual[stagnation_index],
                },
            );
            return error.SoluteTransportSolverStagnated;
        }
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    try residualAt(
        &scratch,
        base,
        current,
        faces,
        diffusive_conductance_m3_per_step,
        mobility_fraction,
        face_parameters,
        face_flux,
        fixed_point,
        residual,
    );
    const final_norm = try scaledNorm(current, residual, options);
    if (final_norm <= 1) {
        if (options.face_flux_mol_by_component) |accepted_flux| {
            try captureAcceptedFaceFlux(
                &scratch,
                base,
                faces,
                diffusive_conductance_m3_per_step,
                mobility_fraction,
                face_parameters,
                face_flux,
                residual,
                accepted_flux,
            );
        }
        @memcpy(state.amount_mol, fixed_point);
        return .{
            .iterations = options.max_iterations,
            .newton_raphson_steps = newton_steps,
            .picard_steps = picard_steps,
            .maximum_scaled_residual = final_norm,
        };
    }
    const limiting_index =
        try worstResidualIndex(current, residual, options);
    std.log.warn(
        "solute transport Newton-Picard exhausted runtime ceiling: max_iterations={d} scaled_residual={e} cell={d} species={d} amount_mol={e} residual_mol={e} newton_steps={d} picard_steps={d}",
        .{
            options.max_iterations,
            final_norm,
            limiting_index / state.species_count,
            limiting_index % state.species_count,
            current[limiting_index],
            residual[limiting_index],
            newton_steps,
            picard_steps,
        },
    );
    return error.SoluteTransportSolverDidNotConverge;
}

fn denseSpeciesNewtonDirection(
    scratch: *transport.State,
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    faces: []const transport.Face,
    conductance: []const f64,
    mobility: []const f64,
    face_parameters: transport.FaceParameters,
    options: Options,
    species: usize,
    face_flux: []f64,
    fixed_point: []f64,
    sampled_state: []f64,
    sampled_residual: []f64,
    matrix: []f64,
    rhs: []f64,
    direction: []f64,
) !bool {
    if (species >= scratch.species_count or
        matrix.len != scratch.cell_count * scratch.cell_count or
        rhs.len != scratch.cell_count or
        direction.len != current.len)
        return error.TransportFaceParameterSizeMismatch;
    for (0..scratch.cell_count) |cell|
        rhs[cell] =
            -residual[cell * scratch.species_count + species];
    for (0..scratch.cell_count) |column| {
        const state_index =
            column * scratch.species_count + species;
        const epsilon = @max(
            1.0e-8,
            @sqrt(std.math.floatEps(f64)) *
                @max(1.0, @abs(current[state_index])),
        );
        @memcpy(sampled_state, current);
        sampled_state[state_index] += epsilon;
        residualAt(
            scratch,
            base,
            sampled_state,
            faces,
            conductance,
            mobility,
            face_parameters,
            face_flux,
            fixed_point,
            sampled_residual,
        ) catch return false;
        for (0..scratch.cell_count) |row| {
            const residual_index =
                row * scratch.species_count + species;
            matrix[row * scratch.cell_count + column] =
                (sampled_residual[residual_index] -
                    residual[residual_index]) /
                epsilon;
        }
    }
    if (!numerics.solveDenseLinearSystem(
        matrix,
        rhs,
        scratch.cell_count,
    )) return false;
    for (rhs, 0..) |value, cell|
        direction[cell * scratch.species_count + species] = value;
    _ = options;
    return true;
}

fn residualAt(scratch: *transport.State, base: []const f64, trial: []const f64, faces: []const transport.Face, conductance: []const f64, mobility: []const f64, face_parameters: transport.FaceParameters, face_flux: []f64, fixed_point: []f64, residual: []f64) !void {
    @memcpy(scratch.amount_mol, trial);
    @memcpy(fixed_point, base);
    for (faces, 0..) |face, face_index| {
        const start = face_index * scratch.species_count;
        try transport.calculateFaceFluxes(scratch, face, conductance[start .. start + scratch.species_count], mobility[start .. start + scratch.species_count], face_parameters, face_flux);
        for (face_flux, 0..) |flux, species_index| {
            const first_index = face.first_cell * scratch.species_count + species_index;
            const second_index = face.second_cell * scratch.species_count + species_index;
            // TRNSFRS traverses faces against the evolving substep inventory.
            // A cell shared by several faces must not donate its original
            // inventory independently to every neighbor.
            const bounded_flux = std.math.clamp(flux, -fixed_point[second_index], fixed_point[first_index]);
            fixed_point[first_index] -= bounded_flux;
            fixed_point[second_index] += bounded_flux;
        }
    }
    for (fixed_point, trial, residual) |target, value, *difference| {
        if (!std.math.isFinite(target) or target < -1e-12) return error.InvalidImplicitTransportCandidate;
        difference.* = target - value;
    }
}

fn captureAcceptedFaceFlux(
    scratch: *transport.State,
    base: []const f64,
    faces: []const transport.Face,
    conductance: []const f64,
    mobility: []const f64,
    face_parameters: transport.FaceParameters,
    face_flux: []f64,
    accumulator: []f64,
    accepted_face_flux: []f64,
) !void {
    if (accepted_face_flux.len != faces.len * scratch.species_count)
        return error.SoluteFaceFluxOutputDimensionMismatch;
    @memcpy(accumulator, base);
    @memset(accepted_face_flux, 0);
    for (faces, 0..) |face, face_index| {
        const start = face_index * scratch.species_count;
        try transport.calculateFaceFluxes(
            scratch,
            face,
            conductance[start..][0..scratch.species_count],
            mobility[start..][0..scratch.species_count],
            face_parameters,
            face_flux,
        );
        for (face_flux, 0..) |flux, species_index| {
            const first_index =
                face.first_cell * scratch.species_count + species_index;
            const second_index =
                face.second_cell * scratch.species_count + species_index;
            const bounded_flux = std.math.clamp(
                flux,
                -accumulator[second_index],
                accumulator[first_index],
            );
            accumulator[first_index] -= bounded_flux;
            accumulator[second_index] += bounded_flux;
            accepted_face_flux[start + species_index] = bounded_flux;
        }
    }
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidImplicitTransportCandidate;
        candidate.* = @max(0, candidate.*);
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteImplicitTransportState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_mol + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn validateInputs(state: *const transport.State, faces: []const transport.Face, conductance: []const f64, mobility: []const f64, options: Options) !void {
    const face_components = try std.math.mul(usize, faces.len, state.species_count);
    if (conductance.len != face_components or mobility.len != face_components) return error.TransportFaceParameterSizeMismatch;
    if (options.face_flux_mol_by_component) |fluxes|
        if (fluxes.len != face_components)
            return error.SoluteFaceFluxOutputDimensionMismatch;
    if (!std.math.isFinite(options.absolute_tolerance_mol) or options.absolute_tolerance_mol <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or options.max_iterations == 0) return error.InvalidSoluteTransportSolverOptions;
    for (faces) |face| if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell) return error.InvalidTransportFace;
}

fn vectorsEqual(a: []const f64, b: []const f64) bool {
    if (a.len != b.len) return false;
    for (a, b) |left, right| if (left != right) return false;
    return true;
}

fn worstResidualIndex(
    state: []const f64,
    residual: []const f64,
    options: Options,
) !usize {
    if (state.len == 0 or state.len != residual.len)
        return error.NonFiniteImplicitTransportState;
    var limiting_index: usize = 0;
    var limiting_norm: f64 = -1;
    for (state, residual, 0..) |value, difference, index| {
        if (!std.math.isFinite(value) or value < 0 or
            !std.math.isFinite(difference))
            return error.NonFiniteImplicitTransportState;
        const norm =
            @abs(difference) /
            (options.absolute_tolerance_mol +
                options.relative_tolerance *
                    @max(1.0, @abs(value)));
        if (norm > limiting_norm) {
            limiting_norm = norm;
            limiting_index = index;
        }
    }
    return limiting_index;
}

test "implicit transport converges before NPH and conserves species" {
    var state = try transport.State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    state.water_volume_m3[0] = 1;
    state.water_volume_m3[1] = 1;
    (try state.cellAmounts(0))[0] = 2;
    var accepted_face_flux_mol = [_]f64{0};
    const result = try solve(std.testing.allocator, &state, &[_]transport.Face{.{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0 }}, &[_]f64{0.1}, &[_]f64{1}, .{ .maximum_convective_fraction = 1 }, .{ .max_iterations = 20, .face_flux_mol_by_component = &accepted_face_flux_mol });
    try std.testing.expect(result.iterations < 20);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
    const first = (try state.cellAmountsConst(0))[0];
    const second = (try state.cellAmountsConst(1))[0];
    try std.testing.expectApproxEqAbs(
        second,
        accepted_face_flux_mol[0],
        1e-12,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 2), first + second, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0 / 3.0), first - second, 1e-8);
}

test "one Newton iteration advances every runtime species block" {
    const species_count: usize = 40;
    var state = try transport.State.init(std.testing.allocator, 2, species_count);
    defer state.deinit();
    @memset(state.water_volume_m3, 1);
    for (0..species_count) |species| {
        state.amount_mol[species] =
            @as(f64, @floatFromInt(species + 1)) * 1e-4;
        state.amount_mol[species_count + species] = 0;
    }
    const conductance = try std.testing.allocator.alloc(f64, species_count);
    defer std.testing.allocator.free(conductance);
    @memset(conductance, 0.1);
    const mobility = try std.testing.allocator.alloc(f64, species_count);
    defer std.testing.allocator.free(mobility);
    @memset(mobility, 1);
    const before = try std.testing.allocator.dupe(f64, state.amount_mol);
    defer std.testing.allocator.free(before);
    const result = try solve(
        std.testing.allocator,
        &state,
        &.{.{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0 }},
        conductance,
        mobility,
        .{ .maximum_convective_fraction = 1 },
        .{ .max_iterations = 4 },
    );
    try std.testing.expect(result.iterations <= 4);
    for (0..species_count) |species|
        try std.testing.expectApproxEqAbs(
            before[species] + before[species_count + species],
            state.amount_mol[species] +
                state.amount_mol[species_count + species],
            1e-14,
        );
}

test "shared donor is bounded across sequential runtime faces" {
    var state = try transport.State.init(std.testing.allocator, 3, 1);
    defer state.deinit();
    @memset(state.water_volume_m3, 1);
    state.amount_mol[0] = 1;
    const result = try solve(
        std.testing.allocator,
        &state,
        &.{
            .{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0 },
            .{ .first_cell = 0, .second_cell = 2, .water_flux_m3_per_step = 0 },
        },
        &.{ 10, 10 },
        &.{ 1, 1 },
        .{ .maximum_convective_fraction = 1 },
        .{ .max_iterations = 80 },
    );
    try std.testing.expect(result.iterations < 80);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.amount_mol[0] + state.amount_mol[1] + state.amount_mol[2], 1e-10);
    for (state.amount_mol) |amount| try std.testing.expect(amount >= 0);
}

test "failed implicit transport leaves state unchanged" {
    var state = try transport.State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    state.water_volume_m3[0] = 1;
    state.water_volume_m3[1] = 1;
    (try state.cellAmounts(0))[0] = 2;
    const before = try std.testing.allocator.dupe(f64, state.amount_mol);
    defer std.testing.allocator.free(before);
    var unpublished_face_flux_mol = [_]f64{123};
    try std.testing.expectError(error.SoluteTransportSolverDidNotConverge, solve(std.testing.allocator, &state, &[_]transport.Face{.{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0 }}, &[_]f64{0.1}, &[_]f64{1}, .{ .maximum_convective_fraction = 1 }, .{ .absolute_tolerance_mol = 1e-20, .relative_tolerance = 1e-20, .max_iterations = 1, .face_flux_mol_by_component = &unpublished_face_flux_mol }));
    try std.testing.expectEqualSlices(f64, before, state.amount_mol);
    try std.testing.expectEqual(@as(f64, 123), unpublished_face_flux_mol[0]);
}
