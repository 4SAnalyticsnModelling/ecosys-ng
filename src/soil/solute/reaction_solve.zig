//! Split out of reaction_solver.zig by tools/split_decl_group.py.
//! Pure code motion: every decl below is an exact line slice.

const std = @import("std");
const chemistry = @import("chemistry_state.zig");
const geochemistry = @import("geochemistry_network.zig");
const reaction_span = @import("conservative_reaction_span.zig");
const __parent = @import("reaction_solver.zig");
const ComplementaritySearchInputs = __parent.ComplementaritySearchInputs;
const CoupledExtentReaction = __parent.CoupledExtentReaction;
const FullNetworkReactionDiagnostic = __parent.FullNetworkReactionDiagnostic;
const IterationDiagnostic = __parent.IterationDiagnostic;
const Options = __parent.Options;
const Result = __parent.Result;
const Workspace = __parent.Workspace;
const applyKineticGeochemistryStep = __parent.applyKineticGeochemistryStep;
const combineResults = __parent.combineResults;
const coupledExtentReactionEnabled = __parent.coupledExtentReactionEnabled;
const coupled_extent_reaction_count = __parent.coupled_extent_reaction_count;
const equilibriumClosureParameters = __parent.equilibriumClosureParameters;
const evaluateAt = __parent.evaluateAt;
const exhaustsLargestResidual = __parent.exhaustsLargestResidual;
const hasKineticGeochemistry = __parent.hasKineticGeochemistry;
const largestScaledResidualIndex = __parent.largestScaledResidualIndex;
const logLargestResidual = __parent.logLargestResidual;
const logTerminalReactionDecomposition = __parent.logTerminalReactionDecomposition;
const logTerminalStagnationComponent = __parent.logTerminalStagnationComponent;
const matrixColumnNorm = __parent.matrixColumnNorm;
const maximumDifference = __parent.maximumDifference;
const maximumMagnitude = __parent.maximumMagnitude;
const phosphateExtentBounds = __parent.phosphateExtentBounds;
const rememberHistory = __parent.rememberHistory;
const residualScale = __parent.residualScale;
const scaledNorm = __parent.scaledNorm;
const selectTraceCandidate = __parent.selectTraceCandidate;
const transformedVectorAdmissible = __parent.transformedVectorAdmissible;
const tryAcceptAndersonCandidate = __parent.tryAcceptAndersonCandidate;
const tryFullNetworkReactionCandidate = __parent.tryFullNetworkReactionCandidate;
const tryPhosphateExtentCandidate = __parent.tryPhosphateExtentCandidate;
const tryRetainedComplementarityCandidate = __parent.tryRetainedComplementarityCandidate;
const tryTransactionalEquilibriumLookahead = __parent.tryTransactionalEquilibriumLookahead;
const validateAqueousMolarity = __parent.validateAqueousMolarity;
const validateOptions = __parent.validateOptions;

/// Heap-owned replay diagnostics. Production solvers pass no trace and incur
/// no allocation or logging dependency.
pub const SolverTrace = struct {
    allocator: std.mem.Allocator,
    entries: []IterationDiagnostic,
    count: usize = 0,
    full_network_comparison_valid: bool = false,
    full_network_comparison_closure: u8 = 0,
    full_network_comparison_iteration: u16 = 0,
    full_network_directional_fraction: f64 = 0,
    full_network_comparison_inventory_fraction: f64 = 0,
    full_network_first_trial_fraction: f64 = 0,
    full_network_first_trial_maximum_scaled_residual: f64 =
        std.math.inf(f64),
    full_network_limiting_component_index: usize = 0,
    full_network_reaction_count: usize = 0,
    full_network_reactions: []FullNetworkReactionDiagnostic,
    full_network_complementarity_search_attempted: bool = false,
    full_network_complementarity_ambiguous_count: usize = 0,
    full_network_complementarity_ambiguous_columns: []usize,
    full_network_complementarity_combination_count: usize = 0,
    full_network_complementarity_consistent_count: usize = 0,
    full_network_complementarity_exact_descent_count: usize = 0,
    full_network_complementarity_minimum_mismatch_count: usize =
        std.math.maxInt(usize),
    full_network_complementarity_minimum_mismatch_mask: u64 = 0,
    full_network_complementarity_mismatch_columns: []usize,
    full_network_complementarity_best_mask: u64 = 0,
    full_network_complementarity_best_predicted_merit: f64 =
        std.math.inf(f64),
    full_network_complementarity_best_exact_merit: f64 =
        std.math.inf(f64),
    full_network_complementarity_blend_attempted: bool = false,
    full_network_complementarity_blend_column: usize = 0,
    full_network_complementarity_blend_fraction: f64 =
        std.math.nan(f64),
    full_network_complementarity_blend_solution: f64 =
        std.math.nan(f64),
    full_network_complementarity_blend_predicted_merit: f64 =
        std.math.inf(f64),
    full_network_complementarity_blend_exact_merit: f64 =
        std.math.inf(f64),
    full_network_current_rates: []f64,
    full_network_directional_rates: []f64,
    full_network_base_residual: []f64,
    full_network_predicted_residual: []f64,
    full_network_realized_residual: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        maximum_entries: usize,
    ) !SolverTrace {
        if (maximum_entries == 0) return error.ZeroSoluteTraceCapacity;
        const entries = try allocator.alloc(
            IterationDiagnostic,
            maximum_entries,
        );
        errdefer allocator.free(entries);
        const component_count = chemistry.State.packedComponentCount();
        const base_residual = try allocator.alloc(f64, component_count);
        errdefer allocator.free(base_residual);
        const predicted_residual =
            try allocator.alloc(f64, component_count);
        errdefer allocator.free(predicted_residual);
        const realized_residual =
            try allocator.alloc(f64, component_count);
        errdefer allocator.free(realized_residual);
        const reactions = try allocator.alloc(
            FullNetworkReactionDiagnostic,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reactions);
        const ambiguous_columns = try allocator.alloc(
            usize,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(ambiguous_columns);
        const mismatch_columns = try allocator.alloc(
            usize,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(mismatch_columns);
        const current_rates = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(current_rates);
        const directional_rates = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        return .{
            .allocator = allocator,
            .entries = entries,
            .full_network_reactions = reactions,
            .full_network_complementarity_ambiguous_columns = ambiguous_columns,
            .full_network_complementarity_mismatch_columns = mismatch_columns,
            .full_network_current_rates = current_rates,
            .full_network_directional_rates = directional_rates,
            .full_network_base_residual = base_residual,
            .full_network_predicted_residual = predicted_residual,
            .full_network_realized_residual = realized_residual,
        };
    }

    pub fn deinit(self: *SolverTrace) void {
        self.allocator.free(self.full_network_directional_rates);
        self.allocator.free(self.full_network_current_rates);
        self.allocator.free(
            self.full_network_complementarity_ambiguous_columns,
        );
        self.allocator.free(
            self.full_network_complementarity_mismatch_columns,
        );
        self.allocator.free(self.full_network_reactions);
        self.allocator.free(self.full_network_realized_residual);
        self.allocator.free(self.full_network_predicted_residual);
        self.allocator.free(self.full_network_base_residual);
        self.allocator.free(self.entries);
        self.* = undefined;
    }

    pub fn recorded(self: *const SolverTrace) []const IterationDiagnostic {
        return self.entries[0..self.count];
    }

    fn append(
        self: *SolverTrace,
        diagnostic: IterationDiagnostic,
    ) !*IterationDiagnostic {
        if (self.count == self.entries.len)
            return error.SoluteTraceCapacityExceeded;
        self.entries[self.count] = diagnostic;
        self.count += 1;
        return &self.entries[self.count - 1];
    }
};

/// Coupled conservative Newton-Raphson/Picard solve for one runtime cell.
/// Newton acceleration changes only conservative stoichiometric reaction
/// extents, so elemental and exchange-site inventories cannot be broken by
/// independent concentration corrections.
pub fn solveCell(allocator: std.mem.Allocator, state: *chemistry.State, cell_index: usize, parameters: chemistry.ReactionParameters, options: Options) !Result {
    var workspace = try Workspace.init(allocator);
    defer workspace.deinit();
    return solveCellWithWorkspace(&workspace, state, cell_index, parameters, options);
}

/// Allocation-free solve for hourly kernels. A workspace is exclusively owned
/// by one worker and can be reused for any number of runtime soil layers.
pub fn solveCellWithWorkspace(workspace: *Workspace, state: *chemistry.State, cell_index: usize, parameters: chemistry.ReactionParameters, options: Options) !Result {
    return solveCellWithWorkspaceAndTrace(
        workspace,
        state,
        cell_index,
        parameters,
        options,
        null,
    );
}

/// Diagnostic solve with the same transactional production path and a
/// caller-owned heap trace. Intended for deterministic failure replay.
pub fn solveCellWithTrace(
    allocator: std.mem.Allocator,
    state: *chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: Options,
    trace: *SolverTrace,
) !Result {
    var workspace = try Workspace.init(allocator);
    defer workspace.deinit();
    trace.count = 0;
    trace.full_network_comparison_valid = false;
    trace.full_network_reaction_count = 0;
    trace.full_network_complementarity_search_attempted = false;
    trace.full_network_complementarity_ambiguous_count = 0;
    trace.full_network_complementarity_combination_count = 0;
    trace.full_network_complementarity_consistent_count = 0;
    trace.full_network_complementarity_exact_descent_count = 0;
    trace.full_network_complementarity_minimum_mismatch_count =
        std.math.maxInt(usize);
    trace.full_network_complementarity_minimum_mismatch_mask = 0;
    trace.full_network_complementarity_best_mask = 0;
    trace.full_network_complementarity_best_predicted_merit =
        std.math.inf(f64);
    trace.full_network_complementarity_best_exact_merit =
        std.math.inf(f64);
    trace.full_network_complementarity_blend_attempted = false;
    trace.full_network_complementarity_blend_column = 0;
    trace.full_network_complementarity_blend_fraction =
        std.math.nan(f64);
    trace.full_network_complementarity_blend_solution =
        std.math.nan(f64);
    trace.full_network_complementarity_blend_predicted_merit =
        std.math.inf(f64);
    trace.full_network_complementarity_blend_exact_merit =
        std.math.inf(f64);
    return solveCellWithWorkspaceAndTrace(
        &workspace,
        state,
        cell_index,
        parameters,
        options,
        trace,
    );
}

pub fn solveCellWithWorkspaceAndTrace(
    workspace: *Workspace,
    state: *chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: Options,
    trace: ?*SolverTrace,
) !Result {
    try validateOptions(options);
    if (cell_index >= state.cell_count)
        return error.ChemistryCellIndexOutOfBounds;
    try validateAqueousMolarity(state, cell_index);
    try state.packCell(cell_index, workspace.rollback_state);
    const original_state = workspace.rollback_state;
    errdefer state.unpackCell(cell_index, original_state) catch
        @panic("validated SOLUTE rollback state could not be restored");

    const equilibrium_parameters = equilibriumClosureParameters(parameters);
    const first = try solveEquilibriumWithWorkspace(
        workspace,
        state,
        cell_index,
        equilibrium_parameters,
        options,
        trace,
        0,
    );
    if (!hasKineticGeochemistry(parameters)) return first;

    if (first.iterations >= options.max_iterations) {
        std.log.warn(
            "SOLUTE geochemistry split exhausted its shared iteration ceiling before post-kinetic equilibrium: cell={d} max_iterations={d}",
            .{ cell_index, options.max_iterations },
        );
        return error.SoluteReactionSolverDidNotConverge;
    }
    try applyKineticGeochemistryStep(state, cell_index, parameters);
    var second_options = options;
    second_options.max_iterations -= first.iterations;
    const second = try solveEquilibriumWithWorkspace(
        workspace,
        state,
        cell_index,
        equilibrium_parameters,
        second_options,
        trace,
        1,
    );
    return combineResults(first, second);
}

pub fn solveEquilibriumWithWorkspace(
    workspace: *Workspace,
    state: *chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: Options,
    trace: ?*SolverTrace,
    closure_index: u8,
) !Result {
    try validateOptions(options);
    if (cell_index >= state.cell_count) return error.ChemistryCellIndexOutOfBounds;
    const current = workspace.current;
    const residual = workspace.residual;
    const probe_state = workspace.probe_state;
    const probe_residual = workspace.probe_residual;
    const candidate_state = workspace.candidate_state;
    const candidate_residual = workspace.candidate_residual;
    const previous_state = workspace.previous_state;
    const previous_residual = workspace.previous_residual;
    const previous_previous_state = workspace.previous_previous_state;
    const previous_previous_residual = workspace.previous_previous_residual;
    const scratch = &workspace.scratch;
    try state.packCell(cell_index, current);

    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var history_count: u2 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        const transformations = try evaluateAt(scratch, current, parameters);
        _ = try transformedVectorAdmissible(scratch, current, transformations, parameters, 1, residual);
        for (residual, current) |*change, value| change.* -= value;
        const current_norm = try scaledNorm(current, residual, options);
        const limiting_index =
            largestScaledResidualIndex(current, residual, options);
        var trace_entry: ?*IterationDiagnostic = null;
        if (trace) |solver_trace| {
            trace_entry = try solver_trace.append(.{
                .closure_index = closure_index,
                .iteration = iteration,
                .limiting_component_index = limiting_index,
                .limiting_state_value = current[limiting_index],
                .limiting_residual = residual[limiting_index],
                .current_maximum_scaled_residual = current_norm,
            });
        }
        if (current_norm <= 1) {
            selectTraceCandidate(
                trace_entry,
                .converged,
                current_norm,
            );
            try state.unpackCell(cell_index, current);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = current_norm, .converged = true };
        }
        for (candidate_state, current, residual) |*target, value, change|
            target.* = value + change;
        if (exhaustsLargestResidual(current, residual, candidate_state, options)) {
            // Activate the nonnegativity constraint before an interior Newton
            // step can keep approaching the same piecewise-reaction boundary.
            // The fixed-point map is conservative and candidate_state is the
            // already validated, substrate-limited state.
            selectTraceCandidate(
                trace_entry,
                .active_inventory_picard,
                std.math.inf(f64),
            );
            @memcpy(current, candidate_state);
            picard_steps += 1;
            continue;
        }
        const full_network_accepted = try tryFullNetworkReactionCandidate(
            workspace,
            scratch,
            current,
            residual,
            transformations,
            candidate_state,
            probe_state,
            probe_residual,
            parameters,
            options,
            current_norm,
            trace_entry,
            trace,
        );
        const full_network_norm = if (full_network_accepted) blk: {
            @memcpy(workspace.reaction_span_rhs, probe_state);
            break :blk try scaledNorm(
                probe_state,
                probe_residual,
                options,
            );
        } else std.math.inf(f64);
        const phosphate_extent_accepted = try tryPhosphateExtentCandidate(
            workspace,
            scratch,
            current,
            residual,
            candidate_state,
            probe_state,
            probe_residual,
            parameters,
            options,
            current_norm,
        );
        const phosphate_extent_norm = if (phosphate_extent_accepted)
            try scaledNorm(probe_state, probe_residual, options)
        else
            std.math.inf(f64);
        if (trace_entry) |entry| {
            entry.full_network_candidate_maximum_scaled_residual =
                full_network_norm;
            entry.phosphate_candidate_maximum_scaled_residual =
                phosphate_extent_norm;
        }
        if (full_network_accepted or phosphate_extent_accepted) {
            rememberHistory(
                current,
                residual,
                previous_state,
                previous_residual,
                previous_previous_state,
                previous_previous_residual,
                &history_count,
            );
            if (!phosphate_extent_accepted or
                full_network_norm < phosphate_extent_norm)
            {
                @memcpy(probe_state, workspace.reaction_span_rhs);
                selectTraceCandidate(
                    trace_entry,
                    .full_network_newton,
                    full_network_norm,
                );
            } else {
                selectTraceCandidate(
                    trace_entry,
                    .phosphate_extent_newton,
                    phosphate_extent_norm,
                );
            }
            @memcpy(current, probe_state);
            newton_steps += 1;
            continue;
        }
        // Continue Newton/multisecant attempts through the complete runtime
        // ceiling. Every remaining reaction in this closure is an equilibrium
        // coordinate; kinetic silicate weathering is operator-split outside
        // it. Falling back permanently to capped source-cycle Picard updates
        // can strand coupled phosphate protonation and surface exchange away
        // from equilibrium.
        var accepted_newton = false;
        if (history_count >= 2) {
            var gram_00: f64 = 0;
            var gram_01: f64 = 0;
            var gram_11: f64 = 0;
            var rhs_0: f64 = 0;
            var rhs_1: f64 = 0;
            for (
                residual,
                previous_residual,
                previous_previous_residual,
                current,
            ) |now, before, older, value| {
                const scale = residualScale(value, options);
                const difference_0 = (now - before) / scale;
                const difference_1 = (before - older) / scale;
                const normalized_now = now / scale;
                gram_00 += difference_0 * difference_0;
                gram_01 += difference_0 * difference_1;
                gram_11 += difference_1 * difference_1;
                rhs_0 += difference_0 * normalized_now;
                rhs_1 += difference_1 * normalized_now;
            }
            const determinant =
                gram_00 * gram_11 - gram_01 * gram_01;
            if (std.math.isFinite(determinant) and
                @abs(determinant) > std.math.floatEps(f64) *
                    @max(1.0, gram_00 * gram_11))
            {
                const gamma_0 =
                    (rhs_0 * gram_11 - rhs_1 * gram_01) / determinant;
                const gamma_1 =
                    (gram_00 * rhs_1 - gram_01 * rhs_0) / determinant;
                for (
                    candidate_state,
                    current,
                    residual,
                    previous_state,
                    previous_residual,
                    previous_previous_state,
                    previous_previous_residual,
                ) |*candidate, now, now_residual, before, before_residual, older, older_residual| {
                    const mapped_now = now + now_residual;
                    const mapped_before = before + before_residual;
                    const mapped_older = older + older_residual;
                    candidate.* = mapped_now -
                        gamma_0 * (mapped_now - mapped_before) -
                        gamma_1 * (mapped_before - mapped_older);
                }
                if (try tryAcceptAndersonCandidate(
                    scratch,
                    current,
                    candidate_state,
                    probe_state,
                    probe_residual,
                    parameters,
                    options,
                    current_norm,
                )) {
                    rememberHistory(
                        current,
                        residual,
                        previous_state,
                        previous_residual,
                        previous_previous_state,
                        previous_previous_residual,
                        &history_count,
                    );
                    @memcpy(current, probe_state);
                    newton_steps += 1;
                    accepted_newton = true;
                    selectTraceCandidate(
                        trace_entry,
                        .anderson_depth_two,
                        try scaledNorm(
                            probe_state,
                            probe_residual,
                            options,
                        ),
                    );
                }
            }
        }
        if (!accepted_newton and history_count >= 1) {
            var numerator: f64 = 0;
            var denominator: f64 = 0;
            for (residual, previous_residual, current) |now, before, value| {
                const scale = residualScale(value, options);
                const delta = (now - before) / scale;
                numerator += delta * (now / scale);
                denominator += delta * delta;
            }
            if (std.math.isFinite(denominator) and
                denominator > std.math.floatEps(f64))
            {
                const gamma = numerator / denominator;
                for (
                    candidate_state,
                    current,
                    residual,
                    previous_state,
                    previous_residual,
                ) |*candidate, now, now_residual, before, before_residual| {
                    const mapped_now = now + now_residual;
                    const mapped_before = before + before_residual;
                    candidate.* =
                        mapped_now - gamma * (mapped_now - mapped_before);
                }
                if (try tryAcceptAndersonCandidate(
                    scratch,
                    current,
                    candidate_state,
                    probe_state,
                    probe_residual,
                    parameters,
                    options,
                    current_norm,
                )) {
                    rememberHistory(
                        current,
                        residual,
                        previous_state,
                        previous_residual,
                        previous_previous_state,
                        previous_previous_residual,
                        &history_count,
                    );
                    @memcpy(current, probe_state);
                    newton_steps += 1;
                    accepted_newton = true;
                    selectTraceCandidate(
                        trace_entry,
                        .anderson_depth_one,
                        try scaledNorm(
                            probe_state,
                            probe_residual,
                            options,
                        ),
                    );
                }
            }
        }
        if (accepted_newton) continue;
        rememberHistory(
            current,
            residual,
            previous_state,
            previous_residual,
            previous_previous_state,
            previous_previous_residual,
            &history_count,
        );
        if (transformedVectorAdmissible(scratch, current, transformations, parameters, options.directional_probe_fraction, probe_state)) |_| {
            const probe_transformations = evaluateAt(scratch, probe_state, parameters) catch null;
            if (probe_transformations) |probe_changes| {
                if (transformedVectorAdmissible(scratch, probe_state, probe_changes, parameters, 1, probe_residual)) |_| {
                    for (probe_residual, probe_state) |*change, value| change.* -= value;
                    var numerator: f64 = 0;
                    var denominator: f64 = 0;
                    var limiting_scaled_residual: f64 = 0;
                    var limiting_residual: f64 = 0;
                    var limiting_derivative: f64 = 0;
                    for (residual, probe_residual, current) |base, probed, value| {
                        const scale = residualScale(value, options);
                        const derivative =
                            (probed - base) /
                            options.directional_probe_fraction / scale;
                        const normalized_base = base / scale;
                        numerator += normalized_base * derivative;
                        denominator += derivative * derivative;
                        if (@abs(normalized_base) >
                            limiting_scaled_residual)
                        {
                            limiting_scaled_residual =
                                @abs(normalized_base);
                            limiting_residual = normalized_base;
                            limiting_derivative = derivative;
                        }
                    }
                    if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                        const least_squares_fraction = -numerator / denominator;
                        const limiting_coordinate_fraction =
                            if (std.math.isFinite(limiting_derivative) and
                            @abs(limiting_derivative) >
                                std.math.floatEps(f64))
                                -limiting_residual / limiting_derivative
                            else
                                least_squares_fraction;
                        // The fixed-point direction is conservative for every
                        // scalar fraction. Target the coordinate that defines
                        // the infinity-norm convergence test; the subsequent
                        // global merit check rejects cross-coordinate harm.
                        const newton_fraction = std.math.clamp(
                            limiting_coordinate_fraction,
                            options.minimum_newton_fraction,
                            options.maximum_newton_fraction,
                        );
                        if (transformedVectorAdmissible(scratch, current, transformations, parameters, newton_fraction, candidate_state)) |_| {
                            if (try tryAcceptAndersonCandidate(
                                scratch,
                                current,
                                candidate_state,
                                probe_state,
                                candidate_residual,
                                parameters,
                                options,
                                current_norm,
                            )) {
                                @memcpy(current, probe_state);
                                newton_steps += 1;
                                accepted_newton = true;
                                selectTraceCandidate(
                                    trace_entry,
                                    .directional_newton,
                                    try scaledNorm(
                                        probe_state,
                                        candidate_residual,
                                        options,
                                    ),
                                );
                            }
                        } else |_| {}
                    }
                } else |_| {}
            }
        } else |_| {}
        if (accepted_newton) continue;

        if (try tryAcceptAndersonCandidate(
            scratch,
            current,
            candidate_state,
            probe_state,
            probe_residual,
            parameters,
            options,
            current_norm,
        )) {
            selectTraceCandidate(
                trace_entry,
                .full_picard,
                try scaledNorm(probe_state, probe_residual, options),
            );
            @memcpy(current, probe_state);
            picard_steps += 1;
            continue;
        }
        _ = try transformedVectorAdmissible(scratch, current, transformations, parameters, options.picard_relaxation, candidate_state);
        if (maximumDifference(current, candidate_state) <= std.math.floatEps(f64) * @max(1.0, maximumMagnitude(current))) {
            if (try tryRetainedComplementarityCandidate(
                workspace,
                scratch,
                current,
                residual,
                transformations,
                candidate_state,
                probe_state,
                probe_residual,
                parameters,
                options,
                current_norm,
            )) {
                selectTraceCandidate(
                    trace_entry,
                    .full_network_newton,
                    try scaledNorm(
                        probe_state,
                        probe_residual,
                        options,
                    ),
                );
                @memcpy(current, probe_state);
                newton_steps += 1;
                continue;
            }
            if (try tryTransactionalEquilibriumLookahead(
                workspace,
                scratch,
                current,
                parameters,
                options,
                options.max_iterations - iteration - 1,
                probe_state,
                probe_residual,
            )) |lookahead_iterations| {
                const accepted_norm = try scaledNorm(
                    probe_state,
                    probe_residual,
                    options,
                );
                selectTraceCandidate(
                    trace_entry,
                    .full_network_newton,
                    accepted_norm,
                );
                try state.unpackCell(cell_index, probe_state);
                return .{
                    .iterations = iteration + 1 +
                        lookahead_iterations,
                    .newton_raphson_steps = newton_steps + 1,
                    .picard_steps = picard_steps,
                    .maximum_scaled_residual = accepted_norm,
                    .converged = true,
                };
            }
            selectTraceCandidate(
                trace_entry,
                .stagnated,
                current_norm,
            );
            logTerminalReactionDecomposition(
                scratch,
                current,
                parameters,
            );
            logTerminalStagnationComponent(
                current,
                residual,
                options,
                limiting_index,
            );
            return error.SoluteReactionSolverStagnated;
        }
        if (try tryAcceptAndersonCandidate(
            scratch,
            current,
            candidate_state,
            probe_state,
            probe_residual,
            parameters,
            options,
            current_norm,
        )) {
            selectTraceCandidate(
                trace_entry,
                .relaxed_picard,
                try scaledNorm(probe_state, probe_residual, options),
            );
            @memcpy(current, probe_state);
            picard_steps += 1;
            continue;
        }
        if (try tryRetainedComplementarityCandidate(
            workspace,
            scratch,
            current,
            residual,
            transformations,
            candidate_state,
            probe_state,
            probe_residual,
            parameters,
            options,
            current_norm,
        )) {
            selectTraceCandidate(
                trace_entry,
                .full_network_newton,
                try scaledNorm(
                    probe_state,
                    probe_residual,
                    options,
                ),
            );
            @memcpy(current, probe_state);
            newton_steps += 1;
            continue;
        }
        if (try tryTransactionalEquilibriumLookahead(
            workspace,
            scratch,
            current,
            parameters,
            options,
            options.max_iterations - iteration - 1,
            probe_state,
            probe_residual,
        )) |lookahead_iterations| {
            const accepted_norm = try scaledNorm(
                probe_state,
                probe_residual,
                options,
            );
            selectTraceCandidate(
                trace_entry,
                .full_network_newton,
                accepted_norm,
            );
            try state.unpackCell(cell_index, probe_state);
            return .{
                .iterations = iteration + 1 +
                    lookahead_iterations,
                .newton_raphson_steps = newton_steps + 1,
                .picard_steps = picard_steps,
                .maximum_scaled_residual = accepted_norm,
                .converged = true,
            };
        }
        selectTraceCandidate(trace_entry, .stagnated, current_norm);
        logTerminalReactionDecomposition(
            scratch,
            current,
            parameters,
        );
        logTerminalStagnationComponent(
            current,
            residual,
            options,
            limiting_index,
        );
        return error.SoluteReactionSolverStagnated;
    }
    const final_changes = try evaluateAt(scratch, current, parameters);
    const final_admissible_fraction =
        try transformedVectorAdmissible(scratch, current, final_changes, parameters, 1, residual);
    for (residual, current) |*change, value| change.* -= value;
    const final_norm = try scaledNorm(current, residual, options);
    if (trace) |solver_trace| {
        const limiting_index =
            largestScaledResidualIndex(current, residual, options);
        const entry = try solver_trace.append(.{
            .closure_index = closure_index,
            .iteration = options.max_iterations,
            .limiting_component_index = limiting_index,
            .limiting_state_value = current[limiting_index],
            .limiting_residual = residual[limiting_index],
            .current_maximum_scaled_residual = final_norm,
        });
        selectTraceCandidate(entry, .iteration_ceiling, final_norm);
    }
    _ = final_admissible_fraction;
    // Legacy STARTE/SOLUTE retain their state after exactly MRXN cycles, but
    // ecosys-ng defines this routine as a convergent hybrid equilibrium solve.
    // Publishing a residual above tolerance would silently contaminate every
    // downstream chemistry pool. The live cell has not been unpacked, so this
    // failure is transactional.
    std.log.warn(
        "SOLUTE reaction network did not converge: cell={d} max_iterations={d} maximum_scaled_residual={e} newton_steps={d} picard_steps={d}",
        .{
            cell_index,
            options.max_iterations,
            final_norm,
            newton_steps,
            picard_steps,
        },
    );
    logLargestResidual(current, residual, options);
    logTerminalReactionDecomposition(scratch, current, parameters);
    return error.SoluteReactionSolverDidNotConverge;
}

pub fn solveComplementarityBlend(
    workspace: *Workspace,
    inputs: ComplementaritySearchInputs,
    selected_jacobian: []const f64,
    opposite_jacobian: []const f64,
    ambiguous_columns: []const usize,
    opposite_mask: u64,
    blend_column: usize,
    blend_fraction: f64,
) ?f64 {
    const matrix_count =
        inputs.current.len * inputs.column_count;
    @memcpy(
        workspace.reaction_span_jacobian[0..matrix_count],
        selected_jacobian,
    );
    for (ambiguous_columns, 0..) |column, bit| {
        if (opposite_mask &
            (@as(u64, 1) << @intCast(bit)) == 0)
        {
            continue;
        }
        for (0..inputs.current.len) |row|
            workspace.reaction_span_jacobian[
                row * inputs.column_count + column
            ] = opposite_jacobian[
                row * inputs.column_count + column
            ];
    }
    for (0..inputs.current.len) |row| {
        const index =
            row * inputs.column_count + blend_column;
        workspace.reaction_span_jacobian[index] =
            selected_jacobian[index] +
            blend_fraction *
                (opposite_jacobian[index] -
                    selected_jacobian[index]);
    }
    if (!solveBoundedReactionSpan(
        workspace,
        inputs.current.len,
        inputs.column_count,
    )) return null;
    return workspace.reaction_span_solution[blend_column];
}

pub fn solveBoundedPhosphateExtents(
    workspace: *Workspace,
    current: []const f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    limiting_index: usize,
) bool {
    const row_count = current.len;
    const column_count = coupled_extent_reaction_count;
    const matrix = workspace.phosphate_extent_jacobian;
    const right_hand_side = workspace.phosphate_extent_rhs;
    const solution = workspace.phosphate_extent_solution;
    const active_bounds = workspace.phosphate_extent_active_bounds;
    @memset(solution, 0);
    for (active_bounds, 0..) |*status, column| {
        const reaction: CoupledExtentReaction = @enumFromInt(column);
        status.* = if (coupledExtentReactionEnabled(
            reaction,
            limiting_index,
        )) 0 else 2;
    }

    var active_set_iteration: usize = 0;
    while (active_set_iteration < 8 * column_count) : (active_set_iteration += 1) {
        var free_count: usize = 0;
        for (active_bounds, 0..) |status, column| {
            if (status != 0) continue;
            workspace.phosphate_extent_free_columns[free_count] = column;
            free_count += 1;
        }
        @memcpy(workspace.phosphate_extent_projected_rhs, right_hand_side);
        for (0..row_count) |row| {
            for (active_bounds, solution, 0..) |status, extent, column| {
                if (status == 0) continue;
                workspace.phosphate_extent_projected_rhs[row] -=
                    matrix[row * column_count + column] * extent;
            }
            for (workspace.phosphate_extent_free_columns[0..free_count], 0..) |
                column,
                free_column,
            | {
                workspace.phosphate_extent_projected_jacobian[
                    row * free_count + free_column
                ] = matrix[row * column_count + column];
            }
        }
        if (free_count > 0) {
            if (!solvePivotedHouseholder(
                workspace.phosphate_extent_projected_jacobian[0 .. row_count * free_count],
                workspace.phosphate_extent_projected_rhs,
                workspace.phosphate_extent_residual[0..free_count],
                workspace.phosphate_extent_pivots[0..free_count],
                workspace.phosphate_extent_probe_residual[0..free_count],
                row_count,
                free_count,
                null,
            )) return false;
            for (workspace.phosphate_extent_free_columns[0..free_count], 0..) |
                column,
                free_column,
            | solution[column] =
                workspace.phosphate_extent_residual[free_column];
        }

        var violating_column: ?usize = null;
        var violating_status: i8 = 0;
        var largest_violation: f64 = 0;
        for (active_bounds, solution, 0..) |status, extent, column| {
            if (status != 0) continue;
            const bounds = phosphateExtentBounds(
                current,
                @enumFromInt(column),
                parameters,
                options,
            );
            const scale = @max(
                1.0,
                @max(
                    @abs(bounds.lower_mol_per_m3),
                    @abs(bounds.upper_mol_per_m3),
                ),
            );
            const lower_violation =
                (bounds.lower_mol_per_m3 - extent) / scale;
            const upper_violation =
                (extent - bounds.upper_mol_per_m3) / scale;
            if (lower_violation > largest_violation) {
                largest_violation = lower_violation;
                violating_column = column;
                violating_status = -1;
            }
            if (upper_violation > largest_violation) {
                largest_violation = upper_violation;
                violating_column = column;
                violating_status = 1;
            }
        }
        if (violating_column) |column| {
            const bounds = phosphateExtentBounds(
                current,
                @enumFromInt(column),
                parameters,
                options,
            );
            active_bounds[column] = violating_status;
            solution[column] = if (violating_status < 0)
                bounds.lower_mol_per_m3
            else
                bounds.upper_mol_per_m3;
            continue;
        }

        for (0..row_count) |row| {
            var linear_residual = -right_hand_side[row];
            for (solution, 0..) |extent, column|
                linear_residual +=
                    matrix[row * column_count + column] * extent;
            workspace.phosphate_extent_projected_rhs[row] = linear_residual;
        }
        const residual_norm =
            maximumMagnitude(workspace.phosphate_extent_projected_rhs);
        var release_column: ?usize = null;
        var largest_kkt_violation: f64 = 0;
        for (active_bounds, 0..) |status, column| {
            if (status == 0 or status == 2) continue;
            var gradient: f64 = 0;
            for (workspace.phosphate_extent_projected_rhs, 0..) |
                linear_residual,
                row,
            | gradient +=
                matrix[row * column_count + column] * linear_residual;
            const column_norm = matrixColumnNorm(
                matrix,
                row_count,
                column_count,
                column,
                0,
            );
            const gradient_scale = @max(1.0, column_norm * residual_norm);
            const violation = if (status < 0)
                -gradient / gradient_scale
            else
                gradient / gradient_scale;
            if (violation > @sqrt(std.math.floatEps(f64)) and
                violation > largest_kkt_violation)
            {
                largest_kkt_violation = violation;
                release_column = column;
            }
        }
        if (release_column) |column| {
            active_bounds[column] = 0;
            continue;
        }
        for (solution) |extent|
            if (!std.math.isFinite(extent)) return false;
        return true;
    }
    return false;
}

pub fn solveProjectedReactionSpanLeastSquares(
    workspace: *Workspace,
    row_count: usize,
    column_count: usize,
) bool {
    const matrix =
        workspace.reaction_span_jacobian[0 .. row_count * column_count];
    const right_hand_side =
        workspace.reaction_span_rhs[0..row_count];
    const solution =
        workspace.reaction_span_solution[0..column_count];
    const lower_bounds =
        workspace.reaction_span_lower_bounds[0..column_count];
    const upper_bounds =
        workspace.reaction_span_upper_bounds[0..column_count];
    @memset(solution, 0);
    for (right_hand_side, 0..) |value, row|
        workspace.reaction_span_projected_rhs[row] = -value;

    var sweep: usize = 0;
    while (sweep < 512) : (sweep += 1) {
        var maximum_scaled_change: f64 = 0;
        for (0..column_count) |column| {
            if (lower_bounds[column] > upper_bounds[column])
                return false;
            var column_norm_squared: f64 = 0;
            var gradient: f64 = 0;
            for (
                workspace.reaction_span_projected_rhs[0..row_count],
                0..,
            ) |linear_residual, row| {
                const coefficient =
                    matrix[row * column_count + column];
                column_norm_squared += coefficient * coefficient;
                gradient += coefficient * linear_residual;
            }
            if (!std.math.isFinite(column_norm_squared) or
                !std.math.isFinite(gradient))
            {
                return false;
            }
            if (column_norm_squared == 0) continue;
            const previous = solution[column];
            const next = std.math.clamp(
                previous - gradient / column_norm_squared,
                lower_bounds[column],
                upper_bounds[column],
            );
            const change = next - previous;
            if (change == 0) continue;
            solution[column] = next;
            for (0..row_count) |row| {
                workspace.reaction_span_projected_rhs[row] +=
                    matrix[row * column_count + column] * change;
            }
            maximum_scaled_change = @max(
                maximum_scaled_change,
                @abs(change) / @max(
                    1.0,
                    @max(
                        @abs(lower_bounds[column]),
                        @abs(upper_bounds[column]),
                    ),
                ),
            );
        }
        if (maximum_scaled_change <=
            64 * std.math.floatEps(f64))
        {
            for (solution) |extent|
                if (!std.math.isFinite(extent)) return false;
            return true;
        }
    }
    for (solution) |extent|
        if (!std.math.isFinite(extent)) return false;
    return true;
}

pub fn solveBoundedReactionSpan(
    workspace: *Workspace,
    row_count: usize,
    column_count: usize,
) bool {
    if (column_count == 0 or
        column_count > reaction_span.reaction_count)
    {
        return false;
    }
    const matrix =
        workspace.reaction_span_jacobian[0 .. row_count * column_count];
    const right_hand_side =
        workspace.reaction_span_rhs[0..row_count];
    const solution =
        workspace.reaction_span_solution[0..column_count];
    const active_bounds =
        workspace.reaction_span_active_bounds[0..column_count];
    const lower_bounds =
        workspace.reaction_span_lower_bounds[0..column_count];
    const upper_bounds =
        workspace.reaction_span_upper_bounds[0..column_count];
    @memset(solution, 0);
    @memset(active_bounds, 0);
    workspace.reaction_span_last_rank = 0;

    var active_set_iteration: usize = 0;
    while (active_set_iteration < 8 * column_count) : (active_set_iteration += 1) {
        var free_count: usize = 0;
        for (active_bounds, 0..) |status, column| {
            if (status != 0) continue;
            workspace.reaction_span_free_columns[free_count] = column;
            free_count += 1;
        }
        @memcpy(
            workspace.reaction_span_projected_rhs[0..row_count],
            right_hand_side,
        );
        for (0..row_count) |row| {
            for (active_bounds, solution, 0..) |status, extent, column| {
                if (status == 0) continue;
                workspace.reaction_span_projected_rhs[row] -=
                    matrix[row * column_count + column] * extent;
            }
            for (workspace.reaction_span_free_columns[0..free_count], 0..) |
                column,
                free_column,
            | {
                workspace.reaction_span_projected_jacobian[
                    row * free_count + free_column
                ] = matrix[row * column_count + column];
            }
        }
        if (free_count > 0) {
            if (!solvePivotedHouseholder(
                workspace.reaction_span_projected_jacobian[0 .. row_count * free_count],
                workspace.reaction_span_projected_rhs[0..row_count],
                workspace.reaction_span_residual[0..free_count],
                workspace.reaction_span_pivots[0..free_count],
                workspace.reaction_span_probe_residual[0..free_count],
                row_count,
                free_count,
                &workspace.reaction_span_last_rank,
            )) return false;
            for (workspace.reaction_span_free_columns[0..free_count], 0..) |
                column,
                free_column,
            | solution[column] =
                workspace.reaction_span_residual[free_column];
        }

        var violating_column: ?usize = null;
        var violating_status: i8 = 0;
        var largest_violation: f64 = 0;
        for (active_bounds, solution, 0..) |status, extent, column| {
            if (status != 0) continue;
            const scale = @max(
                1.0,
                @max(
                    @abs(lower_bounds[column]),
                    @abs(upper_bounds[column]),
                ),
            );
            const lower_violation =
                (lower_bounds[column] - extent) / scale;
            const upper_violation =
                (extent - upper_bounds[column]) / scale;
            if (lower_violation > largest_violation) {
                largest_violation = lower_violation;
                violating_column = column;
                violating_status = -1;
            }
            if (upper_violation > largest_violation) {
                largest_violation = upper_violation;
                violating_column = column;
                violating_status = 1;
            }
        }
        if (violating_column) |column| {
            active_bounds[column] = violating_status;
            solution[column] = if (violating_status < 0)
                lower_bounds[column]
            else
                upper_bounds[column];
            continue;
        }

        for (0..row_count) |row| {
            var linear_residual = -right_hand_side[row];
            for (solution, 0..) |extent, column|
                linear_residual +=
                    matrix[row * column_count + column] * extent;
            workspace.reaction_span_projected_rhs[row] = linear_residual;
        }
        const residual_norm = maximumMagnitude(
            workspace.reaction_span_projected_rhs[0..row_count],
        );
        var release_column: ?usize = null;
        var largest_kkt_violation: f64 = 0;
        for (active_bounds, 0..) |status, column| {
            if (status == 0 or status == 2) continue;
            var gradient: f64 = 0;
            for (
                workspace.reaction_span_projected_rhs[0..row_count],
                0..,
            ) |linear_residual, row| {
                gradient += matrix[
                    row * column_count + column
                ] * linear_residual;
            }
            const column_norm = matrixColumnNorm(
                matrix,
                row_count,
                column_count,
                column,
                0,
            );
            const gradient_scale = @max(1.0, column_norm * residual_norm);
            const violation = if (status < 0)
                -gradient / gradient_scale
            else
                gradient / gradient_scale;
            if (violation > @sqrt(std.math.floatEps(f64)) and
                violation > largest_kkt_violation)
            {
                largest_kkt_violation = violation;
                release_column = column;
            }
        }
        if (release_column) |column| {
            active_bounds[column] = 0;
            continue;
        }
        for (solution) |extent|
            if (!std.math.isFinite(extent)) return false;
        return true;
    }
    return false;
}

pub fn solvePivotedHouseholder(
    matrix: []f64,
    right_hand_side: []f64,
    solution: []f64,
    pivots: []usize,
    permutation_work: []f64,
    row_count: usize,
    column_count: usize,
    rank_output: ?*usize,
) bool {
    if (rank_output) |output| output.* = 0;
    if (row_count < column_count or
        matrix.len != row_count * column_count or
        right_hand_side.len != row_count or
        solution.len != column_count or
        pivots.len != column_count or
        permutation_work.len != column_count)
    {
        return false;
    }
    @memset(solution, 0);
    for (pivots, 0..) |*pivot, index| pivot.* = index;
    var leading_norm: f64 = 0;
    var rank: usize = 0;
    for (0..column_count) |column| {
        var pivot_column = column;
        var pivot_norm: f64 = 0;
        for (column..column_count) |candidate| {
            const norm = matrixColumnNorm(
                matrix,
                row_count,
                column_count,
                candidate,
                column,
            );
            if (norm > pivot_norm) {
                pivot_norm = norm;
                pivot_column = candidate;
            }
        }
        if (column == 0) leading_norm = pivot_norm;
        if (!std.math.isFinite(pivot_norm) or pivot_norm <=
            @sqrt(std.math.floatEps(f64)) * @max(1.0, leading_norm))
        {
            break;
        }
        if (pivot_column != column) {
            for (0..row_count) |row| {
                std.mem.swap(
                    f64,
                    &matrix[row * column_count + column],
                    &matrix[row * column_count + pivot_column],
                );
            }
            std.mem.swap(
                usize,
                &pivots[column],
                &pivots[pivot_column],
            );
        }
        const diagonal_index = column * column_count + column;
        const diagonal = matrix[diagonal_index];
        const reflected_diagonal =
            if (diagonal >= 0) -pivot_norm else pivot_norm;
        const leading_householder = diagonal - reflected_diagonal;
        if (!std.math.isFinite(leading_householder) or
            leading_householder == 0)
        {
            break;
        }
        var squared_tail: f64 = 1;
        for (column + 1..row_count) |row| {
            const index = row * column_count + column;
            matrix[index] /= leading_householder;
            squared_tail += matrix[index] * matrix[index];
        }
        const factor = 2 / squared_tail;
        for (column + 1..column_count) |target_column| {
            var dot = matrix[column * column_count + target_column];
            for (column + 1..row_count) |row| {
                dot += matrix[row * column_count + column] *
                    matrix[row * column_count + target_column];
            }
            dot *= factor;
            matrix[column * column_count + target_column] -= dot;
            for (column + 1..row_count) |row| {
                matrix[row * column_count + target_column] -=
                    matrix[row * column_count + column] * dot;
            }
        }
        var rhs_dot = right_hand_side[column];
        for (column + 1..row_count) |row|
            rhs_dot += matrix[row * column_count + column] *
                right_hand_side[row];
        rhs_dot *= factor;
        right_hand_side[column] -= rhs_dot;
        for (column + 1..row_count) |row|
            right_hand_side[row] -=
                matrix[row * column_count + column] * rhs_dot;
        matrix[diagonal_index] = reflected_diagonal;
        rank += 1;
    }
    if (rank == 0) return false;
    if (rank_output) |output| output.* = rank;
    @memset(permutation_work, 0);
    var reverse_index = rank;
    while (reverse_index > 0) {
        reverse_index -= 1;
        var value = right_hand_side[reverse_index];
        for (reverse_index + 1..rank) |column|
            value -= matrix[reverse_index * column_count + column] *
                permutation_work[column];
        const diagonal =
            matrix[reverse_index * column_count + reverse_index];
        if (!std.math.isFinite(diagonal) or diagonal == 0) return false;
        permutation_work[reverse_index] = value / diagonal;
        if (!std.math.isFinite(permutation_work[reverse_index]))
            return false;
    }
    for (0..rank) |index|
        solution[pivots[index]] = permutation_work[index];
    return true;
}
