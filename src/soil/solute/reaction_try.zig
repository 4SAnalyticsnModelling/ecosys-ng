//! Split out of reaction_solver.zig by tools/split_decl_group.py.
//! Pure code motion: every decl below is an exact line slice.

const std = @import("std");
const chemistry = @import("chemistry_state.zig");
const water_equilibrium = @import("water_equilibrium.zig");
const reaction_span = @import("conservative_reaction_span.zig");
const __parent = @import("reaction_solver.zig");
const ComplementaritySearchInputs = __parent.ComplementaritySearchInputs;
const CoupledExtentReaction = __parent.CoupledExtentReaction;
const IterationDiagnostic = __parent.IterationDiagnostic;
const Options = __parent.Options;
const SolverTrace = __parent.SolverTrace;
const Workspace = __parent.Workspace;
const applyPhosphateExtent = __parent.applyPhosphateExtent;
const captureFullNetworkDirectionalComparison = __parent.captureFullNetworkDirectionalComparison;
const complementarityColumnDerivative = __parent.complementarityColumnDerivative;
const coupledExtentReactionEnabled = __parent.coupledExtentReactionEnabled;
const coupled_extent_reaction_count = __parent.coupled_extent_reaction_count;
const evaluateAt = __parent.evaluateAt;
const evaluateCandidateResidualAtFraction = __parent.evaluateCandidateResidualAtFraction;
const evaluateGlobalResidualAt = __parent.evaluateGlobalResidualAt;
const evaluatePhosphateExtentResiduals = __parent.evaluatePhosphateExtentResiduals;
const evaluateReactionSpanDerivativeColumn = __parent.evaluateReactionSpanDerivativeColumn;
const largestPhosphateZoneScaledResidual = __parent.largestPhosphateZoneScaledResidual;
const largestScaledResidualIndex = __parent.largestScaledResidualIndex;
const matrixColumnNorm = __parent.matrixColumnNorm;
const maximumDifference = __parent.maximumDifference;
const maximumMagnitude = __parent.maximumMagnitude;
const maximumPhosphateExtent = __parent.maximumPhosphateExtent;
const phosphateExtentCharacteristic = __parent.phosphateExtentCharacteristic;
const phosphateExtentControlsPackedIndex = __parent.phosphateExtentControlsPackedIndex;
const phosphateTrustRegionFraction = __parent.phosphateTrustRegionFraction;
const phosphateZoneExtentControlsPackedIndex = __parent.phosphateZoneExtentControlsPackedIndex;
const reactionSpanExtentBounds = __parent.reactionSpanExtentBounds;
const reactionSpanExtentIsSignificant = __parent.reactionSpanExtentIsSignificant;
const reactionSpanPredictedNorm = __parent.reactionSpanPredictedNorm;
const reactionSpanRowWeight = __parent.reactionSpanRowWeight;
const refineComplementarityDirectionalJacobian = __parent.refineComplementarityDirectionalJacobian;
const residualScale = __parent.residualScale;
const scaledNorm = __parent.scaledNorm;
const solveBoundedPhosphateExtents = __parent.solveBoundedPhosphateExtents;
const solveBoundedReactionSpan = __parent.solveBoundedReactionSpan;
const solveProjectedReactionSpanLeastSquares = __parent.solveProjectedReactionSpanLeastSquares;
const transformedVector = __parent.transformedVector;
const transformedVectorAdmissible = __parent.transformedVectorAdmissible;

/// Rank-revealing semismooth Newton step in the complete conservative
/// equilibrium-reaction span. Native extents are normalized by their local
/// inventory boxes before QR so mol/m3, mol/Mg, and the dimensionless Gapon
/// direction can share one numerically meaningful system.
pub fn tryFullNetworkReactionCandidate(
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    current_transformations: chemistry.CellTransformations,
    candidate_state: []f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
    diagnostic: ?*IterationDiagnostic,
    trace: ?*SolverTrace,
) !bool {
    if (diagnostic) |entry|
        entry.full_network_status = .not_attempted;
    // The generalized derivative uses the prospective active set at the
    // fixed-point endpoint and retains the same conservative ledger and
    // global merit acceptance as the realized line search.
    _ = try evaluateAt(scratch, candidate_state, parameters);
    try reaction_span.evaluateRates(
        scratch,
        0,
        parameters,
        workspace.reaction_span_rates,
    );
    const coefficients =
        try scratch.activityCoefficients(0, parameters.fractions);
    var active_count: usize = 0;
    for (workspace.reaction_span_rates, 0..) |rate, reaction| {
        if (!std.math.isFinite(rate))
            return error.NonFiniteSoluteReactionRate;
        if (rate == 0) continue;
        const bounds = try reactionSpanExtentBounds(
            scratch,
            current,
            reaction,
            current_transformations,
            parameters,
            @abs(rate),
            coefficients.monovalent_activity_coefficient,
            candidate_state,
        );
        const extent_scale = @max(
            bounds.negative_native_extent,
            bounds.positive_native_extent,
        );
        if (!std.math.isFinite(extent_scale) or extent_scale <= 0) continue;
        workspace.reaction_span_active_reactions[active_count] =
            reaction;
        workspace.reaction_span_extent_scales[active_count] =
            extent_scale;
        workspace.reaction_span_lower_bounds[active_count] =
            -options.maximum_newton_fraction *
            bounds.negative_native_extent / extent_scale;
        workspace.reaction_span_upper_bounds[active_count] =
            options.maximum_newton_fraction *
            bounds.positive_native_extent / extent_scale;
        active_count += 1;
    }
    if (diagnostic) |entry|
        entry.full_network_active_columns = active_count;
    workspace.reaction_span_active_count = active_count;
    if (active_count == 0) {
        if (diagnostic) |entry|
            entry.full_network_status = .no_active_reactions;
        return false;
    }

    const jacobian =
        workspace.reaction_span_jacobian[0 .. current.len * active_count];
    @memset(jacobian, 0);
    for (global_residual, current, 0..) |value, state_value, row| {
        const scale = residualScale(state_value, options);
        const active_weight = reactionSpanRowWeight(
            value / scale,
            current_norm,
        );
        workspace.reaction_span_rhs[row] =
            -(value / scale) * active_weight;
    }

    var populated_columns: usize = 0;
    for (0..active_count) |column| {
        const reaction =
            workspace.reaction_span_active_reactions[column];
        const rate = workspace.reaction_span_rates[reaction];
        const extent_scale =
            workspace.reaction_span_extent_scales[column];
        const available_normalized = if (rate < 0)
            -workspace.reaction_span_lower_bounds[column] /
                options.maximum_newton_fraction
        else
            workspace.reaction_span_upper_bounds[column] /
                options.maximum_newton_fraction;
        const normalized_probe_magnitude = @min(
            std.math.cbrt(std.math.floatEps(f64)),
            0.125 * available_normalized,
        );
        if (!std.math.isFinite(normalized_probe_magnitude) or
            normalized_probe_magnitude <= 64 * std.math.floatEps(f64))
        {
            continue;
        }
        const normalized_probe =
            if (rate < 0)
                -normalized_probe_magnitude
            else
                normalized_probe_magnitude;
        var probe_transformations =
            reaction_span.zeroTransformations(parameters);
        try reaction_span.addReactionExtent(
            &probe_transformations,
            reaction,
            normalized_probe * extent_scale,
            current_transformations,
            parameters,
        );
        transformedVector(
            scratch,
            current,
            probe_transformations,
            coefficients.monovalent_activity_coefficient,
            parameters.water_activity_product_mol2_per_m6,
            1,
            candidate_state,
        ) catch continue;
        evaluateGlobalResidualAt(
            scratch,
            candidate_state,
            parameters,
            residual_work,
        ) catch continue;
        var column_norm: f64 = 0;
        for (0..current.len) |row| {
            const scale = residualScale(current[row], options);
            const active_weight = reactionSpanRowWeight(
                global_residual[row] / scale,
                current_norm,
            );
            const derivative =
                (residual_work[row] - global_residual[row]) /
                normalized_probe / scale * active_weight;
            jacobian[row * active_count + column] = derivative;
            column_norm = @max(column_norm, @abs(derivative));
        }
        if (std.math.isFinite(column_norm) and column_norm > 0)
            populated_columns += 1;
    }
    if (diagnostic) |entry|
        entry.full_network_populated_columns = populated_columns;
    if (populated_columns == 0) {
        if (diagnostic) |entry|
            entry.full_network_status = .no_jacobian_columns;
        return false;
    }
    if (!solveBoundedReactionSpan(workspace, current.len, active_count)) {
        if (diagnostic) |entry| {
            entry.full_network_status = .bounded_solve_failed;
            entry.full_network_rank = workspace.reaction_span_last_rank;
        }
        return false;
    }
    if (diagnostic) |entry|
        entry.full_network_rank = workspace.reaction_span_last_rank;

    var candidate_transformations =
        reaction_span.zeroTransformations(parameters);
    var has_nonzero_extent = false;
    for (
        workspace.reaction_span_solution[0..active_count],
        0..,
    ) |normalized_extent, column| {
        if (!std.math.isFinite(normalized_extent))
            return error.NonFiniteSoluteReactionExtent;
        if (normalized_extent == 0) continue;
        try reaction_span.addReactionExtent(
            &candidate_transformations,
            workspace.reaction_span_active_reactions[column],
            normalized_extent *
                workspace.reaction_span_extent_scales[column],
            current_transformations,
            parameters,
        );
        has_nonzero_extent = true;
    }
    if (!has_nonzero_extent) {
        if (diagnostic) |entry|
            entry.full_network_status = .zero_extent;
        return false;
    }

    const inventory_fraction = transformedVectorAdmissible(
        scratch,
        current,
        candidate_transformations,
        parameters,
        1,
        candidate_state,
    ) catch {
        if (diagnostic) |entry|
            entry.full_network_status = .inventory_projection_failed;
        return false;
    };
    if (diagnostic) |entry|
        entry.full_network_inventory_fraction = inventory_fraction;
    const predicted_norm = reactionSpanPredictedNorm(
        workspace,
        current,
        global_residual,
        options,
        current_norm,
        active_count,
        inventory_fraction,
    );
    if (diagnostic) |entry|
        entry.full_network_predicted_maximum_scaled_residual =
            predicted_norm;
    if (!std.math.isFinite(predicted_norm) or
        current_norm - predicted_norm <
            1.0e-6 * @max(1.0, current_norm))
    {
        if (diagnostic) |entry|
            entry.full_network_status = .predicted_merit_rejected;
        return false;
    }
    const accepted = try tryAcceptAndersonCandidate(
        scratch,
        current,
        candidate_state,
        accepted_state,
        residual_work,
        parameters,
        options,
        current_norm,
    );
    if (!accepted) {
        if (trace) |solver_trace| {
            if (!solver_trace.full_network_comparison_valid) {
                captureFullNetworkDirectionalComparison(
                    solver_trace,
                    workspace,
                    scratch,
                    current,
                    global_residual,
                    candidate_state,
                    accepted_state,
                    residual_work,
                    parameters,
                    options,
                    current_norm,
                    active_count,
                    inventory_fraction,
                    current_transformations,
                    coefficients.monovalent_activity_coefficient,
                    diagnostic,
                );
            }
        }
        if (diagnostic) |entry|
            entry.full_network_status = .actual_merit_rejected;
        return false;
    }
    const accepted_norm =
        try scaledNorm(accepted_state, residual_work, options);
    if (current_norm - accepted_norm <
        1.0e-6 * @max(1.0, current_norm))
    {
        if (diagnostic) |entry|
            entry.full_network_status = .actual_merit_rejected;
        return false;
    }
    if (diagnostic) |entry|
        entry.full_network_status = .accepted;
    return true;
}

/// Fallback after every retained candidate has failed. Negative, zero, and
/// positive extent faces are solved together; zero uses a Clarke derivative.
/// The caller still accepts only a finite, admissible strict merit decrease.
pub fn tryRetainedComplementarityCandidate(
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    current_transformations: chemistry.CellTransformations,
    candidate_state: []f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
) !bool {
    const column_count = workspace.reaction_span_active_count;
    if (column_count == 0) return false;
    try scratch.unpackCell(0, current);
    const coefficients =
        try scratch.activityCoefficients(0, parameters.fractions);
    return tryFullNetworkComplementarityCandidate(
        workspace,
        scratch,
        current,
        global_residual,
        current_transformations,
        candidate_state,
        accepted_state,
        residual_work,
        parameters,
        options,
        current_norm,
        coefficients.monovalent_activity_coefficient,
        column_count,
    );
}

pub fn tryFullNetworkComplementarityCandidate(
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    current_transformations: chemistry.CellTransformations,
    candidate_state: []f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
    monovalent_activity_coefficient: f64,
    column_count: usize,
) !bool {
    const row_count = current.len;
    const matrix_count = row_count * column_count;
    const selected_jacobian =
        workspace.reaction_span_jacobian[0..matrix_count];
    const negative_jacobian =
        workspace.reaction_span_negative_jacobian[0..matrix_count];
    const positive_jacobian =
        workspace.reaction_span_positive_jacobian[0..matrix_count];
    @memcpy(
        workspace.reaction_span_original_lower_bounds[0..column_count],
        workspace.reaction_span_lower_bounds[0..column_count],
    );
    @memcpy(
        workspace.reaction_span_original_upper_bounds[0..column_count],
        workspace.reaction_span_upper_bounds[0..column_count],
    );
    defer {
        @memcpy(
            workspace.reaction_span_lower_bounds[0..column_count],
            workspace.reaction_span_original_lower_bounds[0..column_count],
        );
        @memcpy(
            workspace.reaction_span_upper_bounds[0..column_count],
            workspace.reaction_span_original_upper_bounds[0..column_count],
        );
        for (0..column_count) |column| {
            const reaction =
                workspace.reaction_span_active_reactions[column];
            const source = if (workspace.reaction_span_rates[reaction] < 0)
                negative_jacobian
            else
                positive_jacobian;
            for (0..row_count) |row| {
                selected_jacobian[row * column_count + column] =
                    source[row * column_count + column];
            }
        }
    }

    const inputs = ComplementaritySearchInputs{
        .scratch = scratch,
        .current = current,
        .global_residual = global_residual,
        .probe_state = candidate_state,
        .probe_residual = residual_work,
        .current_transformations = current_transformations,
        .parameters = parameters,
        .options = options,
        .current_norm = current_norm,
        .monovalent_activity_coefficient = monovalent_activity_coefficient,
        .column_count = column_count,
    };
    for (0..column_count) |column| {
        const reaction =
            workspace.reaction_span_active_reactions[column];
        const selected_negative =
            workspace.reaction_span_rates[reaction] < 0;
        for (0..row_count) |row| {
            const value =
                selected_jacobian[row * column_count + column];
            negative_jacobian[row * column_count + column] = value;
            positive_jacobian[row * column_count + column] = value;
        }
        if (workspace.reaction_span_lower_bounds[column] < 0) {
            evaluateReactionSpanDerivativeColumn(
                inputs,
                workspace,
                column,
                -1,
                negative_jacobian,
            ) catch {};
        }
        if (workspace.reaction_span_upper_bounds[column] > 0) {
            evaluateReactionSpanDerivativeColumn(
                inputs,
                workspace,
                column,
                1,
                positive_jacobian,
            ) catch {};
        }
        workspace.reaction_span_branch_states[column] =
            if (selected_negative) -1 else 1;
    }
    var branch_iteration: usize = 0;
    while (branch_iteration < 16 * column_count) : (branch_iteration += 1) {
        for (0..column_count) |column| {
            const branch =
                workspace.reaction_span_branch_states[column];
            const lower =
                workspace.reaction_span_original_lower_bounds[column];
            const upper =
                workspace.reaction_span_original_upper_bounds[column];
            workspace.reaction_span_lower_bounds[column] = switch (branch) {
                -1 => lower,
                0, 1 => @max(0, lower),
                else => unreachable,
            };
            workspace.reaction_span_upper_bounds[column] = switch (branch) {
                -1, 0 => @min(0, upper),
                1 => upper,
                else => unreachable,
            };
            for (0..row_count) |row| {
                const index = row * column_count + column;
                selected_jacobian[index] =
                    complementarityColumnDerivative(
                        branch,
                        negative_jacobian[index],
                        positive_jacobian[index],
                    );
            }
        }
        if (!solveProjectedReactionSpanLeastSquares(
            workspace,
            row_count,
            column_count,
        )) return false;

        for (0..row_count) |row| {
            var linear_residual =
                -workspace.reaction_span_rhs[row];
            for (
                workspace.reaction_span_solution[0..column_count],
                0..,
            ) |extent, column| {
                linear_residual +=
                    selected_jacobian[
                        row * column_count + column
                    ] * extent;
            }
            workspace.reaction_span_projected_rhs[row] =
                linear_residual;
        }
        const residual_norm = maximumMagnitude(
            workspace.reaction_span_projected_rhs[0..row_count],
        );
        var changed = false;
        for (
            workspace.reaction_span_solution[0..column_count],
            0..,
        ) |extent, column| {
            if (reactionSpanExtentIsSignificant(
                workspace,
                column,
                extent,
            )) continue;
            var negative_gradient: f64 = 0;
            var positive_gradient: f64 = 0;
            for (
                workspace.reaction_span_projected_rhs[0..row_count],
                0..,
            ) |linear_residual, row| {
                negative_gradient +=
                    negative_jacobian[
                        row * column_count + column
                    ] * linear_residual;
                positive_gradient +=
                    positive_jacobian[
                        row * column_count + column
                    ] * linear_residual;
            }
            const negative_available =
                workspace.reaction_span_original_lower_bounds[column] < 0;
            const positive_available =
                workspace.reaction_span_original_upper_bounds[column] > 0;
            const gradient_scale = @max(
                1.0,
                @max(
                    matrixColumnNorm(
                        negative_jacobian,
                        row_count,
                        column_count,
                        column,
                        0,
                    ),
                    matrixColumnNorm(
                        positive_jacobian,
                        row_count,
                        column_count,
                        column,
                        0,
                    ),
                ) * residual_norm,
            );
            const tolerance =
                @sqrt(std.math.floatEps(f64)) * gradient_scale;
            const negative_violation =
                if (negative_available and
                negative_gradient > tolerance)
                    negative_gradient / gradient_scale
                else
                    0;
            const positive_violation =
                if (positive_available and
                positive_gradient < -tolerance)
                    -positive_gradient / gradient_scale
                else
                    0;
            const next_branch: i8 =
                if (negative_violation == 0 and
                positive_violation == 0)
                    0
                else if (negative_violation > positive_violation)
                    -1
                else
                    1;
            if (next_branch !=
                workspace.reaction_span_branch_states[column])
            {
                workspace.reaction_span_branch_states[column] =
                    next_branch;
                changed = true;
            }
        }
        if (changed) continue;

        var candidate_transformations =
            reaction_span.zeroTransformations(parameters);
        var has_nonzero_extent = false;
        for (
            workspace.reaction_span_solution[0..column_count],
            0..,
        ) |normalized_extent, column| {
            if (!std.math.isFinite(normalized_extent))
                return error.NonFiniteSoluteReactionExtent;
            if (normalized_extent == 0) continue;
            try reaction_span.addReactionExtent(
                &candidate_transformations,
                workspace.reaction_span_active_reactions[column],
                normalized_extent *
                    workspace.reaction_span_extent_scales[column],
                current_transformations,
                parameters,
            );
            has_nonzero_extent = true;
        }
        if (!has_nonzero_extent) return false;
        const inventory_fraction = transformedVectorAdmissible(
            scratch,
            current,
            candidate_transformations,
            parameters,
            1,
            candidate_state,
        ) catch return false;
        var step_fraction = inventory_fraction;
        var predicted_norm = std.math.inf(f64);
        var backtrack: u8 = 0;
        while (backtrack < 40) : (backtrack += 1) {
            predicted_norm = reactionSpanPredictedNorm(
                workspace,
                current,
                global_residual,
                options,
                current_norm,
                column_count,
                step_fraction,
            );
            if (std.math.isFinite(predicted_norm) and
                predicted_norm < current_norm)
            {
                break;
            }
            step_fraction *= 0.5;
        }
        if (backtrack == 40) {
            return false;
        }
        _ = transformedVectorAdmissible(
            scratch,
            current,
            candidate_transformations,
            parameters,
            step_fraction,
            candidate_state,
        ) catch return false;
        if (try tryAcceptReactionExtentLineSearch(
            scratch,
            current,
            candidate_transformations,
            step_fraction,
            candidate_state,
            accepted_state,
            residual_work,
            parameters,
            options,
            current_norm,
        )) {
            return true;
        }
        if (tryAcceptExactReactionExtentCandidate(
            scratch,
            candidate_state,
            accepted_state,
            residual_work,
            parameters,
            options,
            current_norm,
        ) catch false) {
            return true;
        }
        if (!try tryAcceptAndersonCandidate(
            scratch,
            current,
            candidate_state,
            accepted_state,
            residual_work,
            parameters,
            options,
            current_norm,
        )) {
            if (refineComplementarityDirectionalJacobian(
                workspace,
                inputs,
                candidate_transformations,
                negative_jacobian,
                positive_jacobian,
                selected_jacobian,
                candidate_state,
                residual_work,
            )) {
                continue;
            }
            return false;
        }
        const accepted_norm =
            try scaledNorm(accepted_state, residual_work, options);
        return accepted_norm < current_norm;
    }
    return false;
}

pub fn tryAcceptReactionExtentLineSearch(
    scratch: *chemistry.State,
    current: []const f64,
    transformations: chemistry.CellTransformations,
    initial_fraction: f64,
    candidate_state: []f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
) !bool {
    var fraction = initial_fraction;
    var attempt: u8 = 0;
    while (attempt < 48) : (attempt += 1) {
        _ = transformedVectorAdmissible(
            scratch,
            current,
            transformations,
            parameters,
            fraction,
            candidate_state,
        ) catch {
            fraction *= 0.5;
            continue;
        };
        try evaluateGlobalResidualAt(
            scratch,
            candidate_state,
            parameters,
            residual_work,
        );
        const candidate_norm =
            try scaledNorm(candidate_state, residual_work, options);
        if (candidate_norm < current_norm) {
            @memcpy(accepted_state, candidate_state);
            return true;
        }
        fraction *= 0.5;
    }
    return false;
}

pub fn tryAcceptExactReactionExtentCandidate(
    scratch: *chemistry.State,
    candidate_state: []const f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
) !bool {
    try evaluateGlobalResidualAt(
        scratch,
        candidate_state,
        parameters,
        residual_work,
    );
    const candidate_norm =
        try scaledNorm(candidate_state, residual_work, options);
    if (candidate_norm >= current_norm) return false;
    @memcpy(accepted_state, candidate_state);
    return true;
}

pub fn tryPhosphateExtentCandidate(
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    candidate_state: []f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
) !bool {
    const limiting_index = largestScaledResidualIndex(
        current,
        global_residual,
        options,
    );
    if (!phosphateExtentControlsPackedIndex(limiting_index)) return false;
    if (!phosphateZoneExtentControlsPackedIndex(limiting_index) and
        largestPhosphateZoneScaledResidual(
            current,
            global_residual,
            options,
        ) <= 1)
    {
        return false;
    }

    const branch_residual = workspace.phosphate_extent_residual;
    try evaluatePhosphateExtentResiduals(
        scratch,
        current,
        parameters,
        branch_residual,
    );
    const jacobian = workspace.phosphate_extent_jacobian;
    @memset(jacobian, 0);
    for (global_residual, current, 0..) |value, state_value, row| {
        const scale = residualScale(state_value, options);
        const normalized_residual = value / scale;
        const active_weight = @max(
            @sqrt(std.math.floatEps(f64)),
            @abs(normalized_residual) / current_norm,
        );
        workspace.phosphate_extent_rhs[row] =
            -normalized_residual * active_weight;
    }

    for (0..coupled_extent_reaction_count) |column| {
        const reaction: CoupledExtentReaction = @enumFromInt(column);
        if (!coupledExtentReactionEnabled(reaction, limiting_index)) continue;
        var direction: f64 = if (branch_residual[column] < 0) -1 else 1;
        var available = maximumPhosphateExtent(
            current,
            reaction,
            direction,
            parameters,
        );
        if (available <= 0) {
            direction = -direction;
            available = maximumPhosphateExtent(
                current,
                reaction,
                direction,
                parameters,
            );
        }
        if (available <= 0) continue;
        const characteristic = phosphateExtentCharacteristic(
            current,
            reaction,
            parameters,
        );
        const nominal_probe =
            @sqrt(std.math.floatEps(f64)) * @max(1.0, characteristic);
        const probe_magnitude = @min(nominal_probe, 0.125 * available);
        if (probe_magnitude <=
            64 * std.math.floatEps(f64) * @max(1.0, characteristic))
        {
            continue;
        }
        @memcpy(workspace.probe_state, current);
        const signed_probe = direction * probe_magnitude;
        if (!applyPhosphateExtent(
            workspace.probe_state,
            reaction,
            signed_probe,
            parameters,
        )) continue;
        if (!tryProjectWaterPair(
            scratch,
            workspace.probe_state,
            parameters,
        )) continue;
        evaluateGlobalResidualAt(
            scratch,
            workspace.probe_state,
            parameters,
            residual_work,
        ) catch continue;
        for (0..current.len) |row| {
            const scale = residualScale(current[row], options);
            const active_weight = @max(
                @sqrt(std.math.floatEps(f64)),
                @abs(global_residual[row] / scale) / current_norm,
            );
            jacobian[row * coupled_extent_reaction_count + column] =
                (residual_work[row] - global_residual[row]) /
                signed_probe / scale * active_weight;
        }
    }

    if (solveBoundedPhosphateExtents(
        workspace,
        current,
        parameters,
        options,
        limiting_index,
    )) {
        @memcpy(candidate_state, current);
        var has_nonzero_extent = false;
        for (workspace.phosphate_extent_solution, 0..) |extent, index| {
            if (!std.math.isFinite(extent) or extent == 0) continue;
            const reaction: CoupledExtentReaction = @enumFromInt(index);
            if (!applyPhosphateExtent(
                candidate_state,
                reaction,
                extent,
                parameters,
            )) break;
            has_nonzero_extent = true;
        }
        if (has_nonzero_extent and
            try tryAcceptAndersonCandidate(
                scratch,
                current,
                candidate_state,
                accepted_state,
                residual_work,
                parameters,
                options,
                current_norm,
            ))
        {
            const accepted_norm = try scaledNorm(
                accepted_state,
                residual_work,
                options,
            );
            if (current_norm - accepted_norm >=
                1.0e-6 * @max(1.0, current_norm))
            {
                return true;
            }
        }
    }
    return false;
}

pub fn tryProjectWaterPair(
    scratch: *chemistry.State,
    vector: []f64,
    parameters: chemistry.ReactionParameters,
) bool {
    scratch.unpackCell(0, vector) catch return false;
    const coefficients =
        scratch.activityCoefficients(0, parameters.fractions) catch
            return false;
    const water = water_equilibrium.projectProvisional(
        scratch.aqueous[0].hydrogen,
        scratch.aqueous[0].hydroxide,
        coefficients.monovalent_activity_coefficient,
        parameters.water_activity_product_mol2_per_m6,
    ) catch return false;
    scratch.aqueous[0].hydrogen =
        water.hydrogen_concentration_mol_per_m3;
    scratch.aqueous[0].hydroxide =
        water.hydroxide_concentration_mol_per_m3;
    scratch.packCell(0, vector) catch return false;
    return true;
}

/// Depth-one Anderson acceleration of the complete conservative fixed-point
/// map. Both mapped endpoints share the same elemental and site totals, so
/// their affine secant candidate retains those totals. Backtracking preserves
/// nonnegative inventories; H+/OH- are then projected onto Kw.
pub fn tryAcceptAndersonCandidate(
    scratch: *chemistry.State,
    current: []const f64,
    target: []const f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
) !bool {
    var fraction = phosphateTrustRegionFraction(
        current,
        target,
        options,
        parameters.phosphate_kinetics.substrate_limit_fraction,
    );
    var attempt: u8 = 0;
    while (attempt < 32) : (attempt += 1) {
        const candidate_norm = evaluateCandidateResidualAtFraction(
            scratch,
            current,
            target,
            accepted_state,
            residual_work,
            parameters,
            options,
            fraction,
        ) catch {
            fraction *= 0.5;
            continue;
        };
        if (candidate_norm < current_norm) return true;
        fraction *= 0.5;
    }
    return false;
}

/// Proves a complete relaxed fixed-point path off to the side, then publishes
/// only its converged endpoint. No intermediate lookahead state is realized.
pub fn tryTransactionalEquilibriumLookahead(
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    maximum_iterations: u16,
    accepted_state: []f64,
    accepted_residual: []f64,
) !?u16 {
    if (maximum_iterations == 0) return null;
    @memcpy(workspace.reaction_span_best_state, current);
    var iteration: u16 = 0;
    while (iteration < maximum_iterations) : (iteration += 1) {
        const transformations = evaluateAt(
            scratch,
            workspace.reaction_span_best_state,
            parameters,
        ) catch return null;
        _ = transformedVectorAdmissible(
            scratch,
            workspace.reaction_span_best_state,
            transformations,
            parameters,
            options.picard_relaxation,
            workspace.reaction_span_best_residual,
        ) catch return null;
        evaluateGlobalResidualAt(
            scratch,
            workspace.reaction_span_best_residual,
            parameters,
            accepted_residual,
        ) catch return null;
        const norm = try scaledNorm(
            workspace.reaction_span_best_residual,
            accepted_residual,
            options,
        );
        if (norm <= 1) {
            @memcpy(
                accepted_state,
                workspace.reaction_span_best_residual,
            );
            return iteration + 1;
        }
        if (maximumDifference(
            workspace.reaction_span_best_state,
            workspace.reaction_span_best_residual,
        ) <= std.math.floatEps(f64) *
            @max(
                1.0,
                maximumMagnitude(workspace.reaction_span_best_state),
            ))
        {
            return null;
        }
        @memcpy(
            workspace.reaction_span_best_state,
            workspace.reaction_span_best_residual,
        );
    }
    return null;
}
