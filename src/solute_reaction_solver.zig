const std = @import("std");
const chemistry = @import("solute_chemistry_state.zig");
const aqueous_network = @import("solute_aqueous_network.zig");
const phosphate_network = @import("solute_phosphate_network.zig");
const cation_exchange = @import("solute_cation_exchange.zig");
const geochemistry = @import("solute_geochemistry_network.zig");
const aqueous_rates = @import("solute_aqueous_reaction_rates.zig");
const phosphate_rates = @import("solute_phosphate_reaction_rates.zig");
const geochemistry_rates = @import("solute_geochemistry_reaction_rates.zig");
const water_equilibrium = @import("solute_water_equilibrium.zig");
const reaction_span = @import("solute_conservative_reaction_span.zig");

pub const Options = struct {
    absolute_tolerance: f64 = 1e-11,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
    /// `MRXN=60` in SOLUTE.F; convergence exits before this ceiling.
    max_iterations: u16 = 60,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
    converged: bool,
};

pub const CandidateKind = enum {
    none,
    converged,
    active_inventory_picard,
    full_network_newton,
    phosphate_extent_newton,
    anderson_depth_two,
    anderson_depth_one,
    directional_newton,
    full_picard,
    relaxed_picard,
    stagnated,
    iteration_ceiling,
};

pub const FullNetworkCandidateStatus = enum {
    not_attempted,
    no_active_reactions,
    no_jacobian_columns,
    bounded_solve_failed,
    zero_extent,
    inventory_projection_failed,
    predicted_merit_rejected,
    actual_merit_rejected,
    accepted,
};

pub const IterationDiagnostic = struct {
    closure_index: u8,
    iteration: u16,
    limiting_component_index: usize,
    limiting_state_value: f64,
    limiting_residual: f64,
    current_maximum_scaled_residual: f64,
    full_network_candidate_maximum_scaled_residual: f64 =
        std.math.inf(f64),
    phosphate_candidate_maximum_scaled_residual: f64 =
        std.math.inf(f64),
    full_network_status: FullNetworkCandidateStatus = .not_attempted,
    full_network_active_columns: usize = 0,
    full_network_populated_columns: usize = 0,
    full_network_rank: usize = 0,
    full_network_inventory_fraction: f64 = 0,
    full_network_predicted_maximum_scaled_residual: f64 =
        std.math.inf(f64),
    selected_candidate: CandidateKind = .none,
    selected_maximum_scaled_residual: f64 = std.math.inf(f64),
};

pub const FullNetworkReactionDiagnostic = struct {
    reaction_index: usize,
    selection_rate: f64,
    current_rate: f64,
    directional_rate: f64,
    native_extent_scale: f64,
    normalized_lower_bound: f64,
    normalized_upper_bound: f64,
    normalized_solution: f64,
    limiting_normalized_residual_derivative: f64,
    solution_side_limiting_normalized_residual_derivative: f64,
};

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

pub const DiagnosticDecomposition = struct {
    aqueous_hydrogen_mol_per_m3: f64,
    aqueous_hydroxide_mol_per_m3: f64,
    phosphate_hydrogen_mol_per_m3: f64,
    phosphate_hydroxide_mol_per_m3: f64,
    cation_exchange_hydrogen_mol_per_m3: f64,
    carboxyl_hydrogen_mol_per_m3: f64,
    geochemistry_hydrogen_mol_per_m3: f64,
    geochemistry_hydroxide_mol_per_m3: f64,
    assembled_hydrogen_mol_per_m3: f64,
    assembled_hydroxide_mol_per_m3: f64,
};

pub const Hpo4Diagnostic = struct {
    po4_protonation_mol_p_per_m3: f64,
    hpo4_protonation_mol_p_per_m3: f64,
    metal_pairing_mol_p_per_m3: f64,
    surface_adsorption_mol_p_per_m3: f64,
    assembled_change_mol_p_per_m3: f64,
};

const CoupledExtentReaction = enum {
    non_band_hpo4_protonation,
    non_band_iron_hpo4_pairing,
    non_band_iron_h2po4_pairing,
    non_band_calcium_hpo4_pairing,
    non_band_calcium_h2po4_pairing,
    non_band_h2po4_protonated_site_exchange,
    non_band_h2po4_hydroxyl_site_exchange,
    non_band_hpo4_hydroxyl_site_exchange,
    band_hpo4_protonation,
    band_iron_hpo4_pairing,
    band_iron_h2po4_pairing,
    band_calcium_hpo4_pairing,
    band_calcium_h2po4_pairing,
    band_h2po4_protonated_site_exchange,
    band_h2po4_hydroxyl_site_exchange,
    band_hpo4_hydroxyl_site_exchange,
    aqueous_calcium_hydroxide_pairing,
    aqueous_calcium_carbonate_pairing,
    aqueous_calcium_bicarbonate_pairing,
    aqueous_calcium_sulfate_pairing,
};

const coupled_extent_reaction_count =
    @typeInfo(CoupledExtentReaction).@"enum".fields.len;

pub const SecondClosureDiagnostic = struct {
    first_closure_iterations: u16,
    hpo4: Hpo4Diagnostic,
};

pub fn diagnoseSecondClosureStart(
    allocator: std.mem.Allocator,
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: Options,
) !SecondClosureDiagnostic {
    var staged = try chemistry.State.init(allocator, 1);
    defer staged.deinit();
    const packed_state = try allocator.alloc(f64, chemistry.State.packedComponentCount());
    defer allocator.free(packed_state);
    try state.packCell(cell_index, packed_state);
    try staged.unpackCell(0, packed_state);
    var workspace = try Workspace.init(allocator);
    defer workspace.deinit();
    const first = try solveEquilibriumWithWorkspace(
        &workspace,
        &staged,
        0,
        equilibriumClosureParameters(parameters),
        options,
        null,
        0,
    );
    try applyKineticGeochemistryStep(&staged, 0, parameters);
    return .{
        .first_closure_iterations = first.iterations,
        .hpo4 = try diagnoseNonBandHpo4(
            &staged,
            0,
            equilibriumClosureParameters(parameters),
        ),
    };
}

pub fn diagnoseNonBandHpo4(
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
) !Hpo4Diagnostic {
    const coefficients = try state.activityCoefficients(cell_index, parameters.fractions);
    const fluxes = try phosphate_rates.calculate(
        state.aqueous[cell_index],
        state.non_band_phosphate[cell_index],
        coefficients,
        parameters.non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3,
        parameters.phosphate_constants,
        parameters.phosphate_surface,
        parameters.phosphate_minerals,
        parameters.phosphate_kinetics,
    );
    const aqueous = fluxes.aqueous;
    const surface = fluxes.surface.hpo4_with_hydroxyl_site_mol_p_per_megagram *
        fluxes.soil_mass_per_water_volume_megagrams_per_m3;
    const pairing = aqueous.iron_hpo4_pairing_mol_p_per_m3 +
        aqueous.calcium_hpo4_pairing_mol_p_per_m3 +
        aqueous.magnesium_hpo4_pairing_mol_p_per_m3;
    return .{
        .po4_protonation_mol_p_per_m3 = aqueous.po4_hydrogen_association_mol_p_per_m3,
        .hpo4_protonation_mol_p_per_m3 = -aqueous.hpo4_hydrogen_association_mol_p_per_m3,
        .metal_pairing_mol_p_per_m3 = -pairing,
        .surface_adsorption_mol_p_per_m3 = -surface,
        .assembled_change_mol_p_per_m3 = aqueous.po4_hydrogen_association_mol_p_per_m3 -
            aqueous.hpo4_hydrogen_association_mol_p_per_m3 -
            pairing - surface,
    };
}

pub fn diagnoseCell(
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
) !DiagnosticDecomposition {
    const transformations = try state.evaluateCell(cell_index, parameters);
    const assembled = try chemistry.State.assembledAqueousChanges(
        transformations,
        state.cation_exchange_mol_per_megagram[cell_index],
    );
    const phosphate_hydrogen =
        transformations.non_band_phosphate.dissolved_hydrogen_mol_per_m3 *
        parameters.fractions.phosphate_non_band +
        transformations.band_phosphate.dissolved_hydrogen_mol_per_m3 *
            parameters.fractions.phosphate_band;
    const phosphate_hydroxide =
        transformations.non_band_phosphate.dissolved_hydroxide_mol_per_m3 *
        parameters.fractions.phosphate_non_band +
        transformations.band_phosphate.dissolved_hydroxide_mol_per_m3 *
            parameters.fractions.phosphate_band;
    return .{
        .aqueous_hydrogen_mol_per_m3 = transformations.aqueous.hydrogen,
        .aqueous_hydroxide_mol_per_m3 = transformations.aqueous.hydroxide,
        .phosphate_hydrogen_mol_per_m3 = phosphate_hydrogen,
        .phosphate_hydroxide_mol_per_m3 = phosphate_hydroxide,
        .cation_exchange_hydrogen_mol_per_m3 = -transformations.cation_adsorption_mol_per_megagram.hydrogen *
            parameters.cation_exchange_water_ratios.shared_megagrams_per_m3,
        .carboxyl_hydrogen_mol_per_m3 = -transformations.carboxyl_hydrogen_change_mol_per_megagram *
            parameters.cation_exchange_water_ratios.shared_megagrams_per_m3,
        .geochemistry_hydrogen_mol_per_m3 = transformations.geochemistry.dissolved_hydrogen_mol_per_m3,
        .geochemistry_hydroxide_mol_per_m3 = transformations.geochemistry.dissolved_hydroxide_mol_per_m3,
        .assembled_hydrogen_mol_per_m3 = assembled.hydrogen,
        .assembled_hydroxide_mol_per_m3 = assembled.hydroxide,
    };
}

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    current: []f64,
    residual: []f64,
    probe_state: []f64,
    probe_residual: []f64,
    candidate_state: []f64,
    candidate_residual: []f64,
    previous_state: []f64,
    previous_residual: []f64,
    previous_previous_state: []f64,
    previous_previous_residual: []f64,
    rollback_state: []f64,
    phosphate_extent_jacobian: []f64,
    phosphate_extent_rhs: []f64,
    phosphate_extent_projected_jacobian: []f64,
    phosphate_extent_projected_rhs: []f64,
    phosphate_extent_solution: []f64,
    phosphate_extent_residual: []f64,
    phosphate_extent_probe_residual: []f64,
    phosphate_extent_pivots: []usize,
    phosphate_extent_free_columns: []usize,
    phosphate_extent_active_bounds: []i8,
    reaction_span_jacobian: []f64,
    reaction_span_negative_jacobian: []f64,
    reaction_span_positive_jacobian: []f64,
    reaction_span_rhs: []f64,
    reaction_span_projected_jacobian: []f64,
    reaction_span_projected_rhs: []f64,
    reaction_span_solution: []f64,
    reaction_span_best_state: []f64,
    reaction_span_best_residual: []f64,
    reaction_span_residual: []f64,
    reaction_span_probe_residual: []f64,
    reaction_span_pivots: []usize,
    reaction_span_free_columns: []usize,
    reaction_span_active_bounds: []i8,
    reaction_span_active_reactions: []usize,
    reaction_span_lower_bounds: []f64,
    reaction_span_upper_bounds: []f64,
    reaction_span_original_lower_bounds: []f64,
    reaction_span_original_upper_bounds: []f64,
    reaction_span_extent_scales: []f64,
    reaction_span_rates: []f64,
    reaction_span_branch_states: []i8,
    reaction_span_active_count: usize,
    reaction_span_last_rank: usize,
    scratch: chemistry.State,

    pub fn init(allocator: std.mem.Allocator) !Workspace {
        const count = chemistry.State.packedComponentCount();
        const current = try allocator.alloc(f64, count);
        errdefer allocator.free(current);
        const residual = try allocator.alloc(f64, count);
        errdefer allocator.free(residual);
        const probe_state = try allocator.alloc(f64, count);
        errdefer allocator.free(probe_state);
        const probe_residual = try allocator.alloc(f64, count);
        errdefer allocator.free(probe_residual);
        const candidate_state = try allocator.alloc(f64, count);
        errdefer allocator.free(candidate_state);
        const candidate_residual = try allocator.alloc(f64, count);
        errdefer allocator.free(candidate_residual);
        const previous_state = try allocator.alloc(f64, count);
        errdefer allocator.free(previous_state);
        const previous_residual = try allocator.alloc(f64, count);
        errdefer allocator.free(previous_residual);
        const previous_previous_state = try allocator.alloc(f64, count);
        errdefer allocator.free(previous_previous_state);
        const previous_previous_residual = try allocator.alloc(f64, count);
        errdefer allocator.free(previous_previous_residual);
        const rollback_state = try allocator.alloc(f64, count);
        errdefer allocator.free(rollback_state);
        const phosphate_extent_jacobian = try allocator.alloc(
            f64,
            count * coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_jacobian);
        const phosphate_extent_rhs = try allocator.alloc(
            f64,
            count,
        );
        errdefer allocator.free(phosphate_extent_rhs);
        const phosphate_extent_projected_jacobian = try allocator.alloc(
            f64,
            count * coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_projected_jacobian);
        const phosphate_extent_projected_rhs = try allocator.alloc(
            f64,
            count,
        );
        errdefer allocator.free(phosphate_extent_projected_rhs);
        const phosphate_extent_solution = try allocator.alloc(
            f64,
            coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_solution);
        const phosphate_extent_residual = try allocator.alloc(
            f64,
            coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_residual);
        const phosphate_extent_probe_residual = try allocator.alloc(
            f64,
            coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_probe_residual);
        const phosphate_extent_pivots = try allocator.alloc(
            usize,
            coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_pivots);
        const phosphate_extent_free_columns = try allocator.alloc(
            usize,
            coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_free_columns);
        const phosphate_extent_active_bounds = try allocator.alloc(
            i8,
            coupled_extent_reaction_count,
        );
        errdefer allocator.free(phosphate_extent_active_bounds);
        const reaction_span_jacobian = try allocator.alloc(
            f64,
            count * reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_jacobian);
        const reaction_span_negative_jacobian = try allocator.alloc(
            f64,
            count * reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_negative_jacobian);
        const reaction_span_positive_jacobian = try allocator.alloc(
            f64,
            count * reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_positive_jacobian);
        const reaction_span_rhs = try allocator.alloc(f64, count);
        errdefer allocator.free(reaction_span_rhs);
        const reaction_span_projected_jacobian = try allocator.alloc(
            f64,
            count * reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_projected_jacobian);
        const reaction_span_projected_rhs = try allocator.alloc(f64, count);
        errdefer allocator.free(reaction_span_projected_rhs);
        const reaction_span_solution = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_solution);
        const reaction_span_best_state = try allocator.alloc(f64, count);
        errdefer allocator.free(reaction_span_best_state);
        const reaction_span_best_residual = try allocator.alloc(f64, count);
        errdefer allocator.free(reaction_span_best_residual);
        const reaction_span_residual = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_residual);
        const reaction_span_probe_residual = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_probe_residual);
        const reaction_span_pivots = try allocator.alloc(
            usize,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_pivots);
        const reaction_span_free_columns = try allocator.alloc(
            usize,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_free_columns);
        const reaction_span_active_bounds = try allocator.alloc(
            i8,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_active_bounds);
        const reaction_span_active_reactions = try allocator.alloc(
            usize,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_active_reactions);
        const reaction_span_lower_bounds = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_lower_bounds);
        const reaction_span_upper_bounds = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_upper_bounds);
        const reaction_span_original_lower_bounds = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_original_lower_bounds);
        const reaction_span_original_upper_bounds = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_original_upper_bounds);
        const reaction_span_extent_scales = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_extent_scales);
        const reaction_span_rates = try allocator.alloc(
            f64,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_rates);
        const reaction_span_branch_states = try allocator.alloc(
            i8,
            reaction_span.reaction_count,
        );
        errdefer allocator.free(reaction_span_branch_states);
        const scratch = try chemistry.State.init(allocator, 1);
        return .{
            .allocator = allocator,
            .current = current,
            .residual = residual,
            .probe_state = probe_state,
            .probe_residual = probe_residual,
            .candidate_state = candidate_state,
            .candidate_residual = candidate_residual,
            .previous_state = previous_state,
            .previous_residual = previous_residual,
            .previous_previous_state = previous_previous_state,
            .previous_previous_residual = previous_previous_residual,
            .rollback_state = rollback_state,
            .phosphate_extent_jacobian = phosphate_extent_jacobian,
            .phosphate_extent_rhs = phosphate_extent_rhs,
            .phosphate_extent_projected_jacobian = phosphate_extent_projected_jacobian,
            .phosphate_extent_projected_rhs = phosphate_extent_projected_rhs,
            .phosphate_extent_solution = phosphate_extent_solution,
            .phosphate_extent_residual = phosphate_extent_residual,
            .phosphate_extent_probe_residual = phosphate_extent_probe_residual,
            .phosphate_extent_pivots = phosphate_extent_pivots,
            .phosphate_extent_free_columns = phosphate_extent_free_columns,
            .phosphate_extent_active_bounds = phosphate_extent_active_bounds,
            .reaction_span_jacobian = reaction_span_jacobian,
            .reaction_span_negative_jacobian = reaction_span_negative_jacobian,
            .reaction_span_positive_jacobian = reaction_span_positive_jacobian,
            .reaction_span_rhs = reaction_span_rhs,
            .reaction_span_projected_jacobian = reaction_span_projected_jacobian,
            .reaction_span_projected_rhs = reaction_span_projected_rhs,
            .reaction_span_solution = reaction_span_solution,
            .reaction_span_best_state = reaction_span_best_state,
            .reaction_span_best_residual = reaction_span_best_residual,
            .reaction_span_residual = reaction_span_residual,
            .reaction_span_probe_residual = reaction_span_probe_residual,
            .reaction_span_pivots = reaction_span_pivots,
            .reaction_span_free_columns = reaction_span_free_columns,
            .reaction_span_active_bounds = reaction_span_active_bounds,
            .reaction_span_active_reactions = reaction_span_active_reactions,
            .reaction_span_lower_bounds = reaction_span_lower_bounds,
            .reaction_span_upper_bounds = reaction_span_upper_bounds,
            .reaction_span_original_lower_bounds = reaction_span_original_lower_bounds,
            .reaction_span_original_upper_bounds = reaction_span_original_upper_bounds,
            .reaction_span_extent_scales = reaction_span_extent_scales,
            .reaction_span_rates = reaction_span_rates,
            .reaction_span_branch_states = reaction_span_branch_states,
            .reaction_span_active_count = 0,
            .reaction_span_last_rank = 0,
            .scratch = scratch,
        };
    }

    pub fn deinit(self: *Workspace) void {
        self.scratch.deinit();
        self.allocator.free(self.reaction_span_branch_states);
        self.allocator.free(self.reaction_span_rates);
        self.allocator.free(self.reaction_span_extent_scales);
        self.allocator.free(self.reaction_span_original_upper_bounds);
        self.allocator.free(self.reaction_span_original_lower_bounds);
        self.allocator.free(self.reaction_span_upper_bounds);
        self.allocator.free(self.reaction_span_lower_bounds);
        self.allocator.free(self.reaction_span_active_reactions);
        self.allocator.free(self.reaction_span_active_bounds);
        self.allocator.free(self.reaction_span_free_columns);
        self.allocator.free(self.reaction_span_pivots);
        self.allocator.free(self.reaction_span_probe_residual);
        self.allocator.free(self.reaction_span_residual);
        self.allocator.free(self.reaction_span_best_residual);
        self.allocator.free(self.reaction_span_best_state);
        self.allocator.free(self.reaction_span_solution);
        self.allocator.free(self.reaction_span_projected_rhs);
        self.allocator.free(self.reaction_span_projected_jacobian);
        self.allocator.free(self.reaction_span_rhs);
        self.allocator.free(self.reaction_span_positive_jacobian);
        self.allocator.free(self.reaction_span_negative_jacobian);
        self.allocator.free(self.reaction_span_jacobian);
        self.allocator.free(self.phosphate_extent_active_bounds);
        self.allocator.free(self.phosphate_extent_free_columns);
        self.allocator.free(self.phosphate_extent_pivots);
        self.allocator.free(self.phosphate_extent_probe_residual);
        self.allocator.free(self.phosphate_extent_residual);
        self.allocator.free(self.phosphate_extent_solution);
        self.allocator.free(self.phosphate_extent_projected_rhs);
        self.allocator.free(self.phosphate_extent_projected_jacobian);
        self.allocator.free(self.phosphate_extent_rhs);
        self.allocator.free(self.phosphate_extent_jacobian);
        self.allocator.free(self.rollback_state);
        self.allocator.free(self.previous_previous_residual);
        self.allocator.free(self.previous_previous_state);
        self.allocator.free(self.previous_residual);
        self.allocator.free(self.previous_state);
        self.allocator.free(self.candidate_residual);
        self.allocator.free(self.candidate_state);
        self.allocator.free(self.probe_residual);
        self.allocator.free(self.probe_state);
        self.allocator.free(self.residual);
        self.allocator.free(self.current);
        self.* = undefined;
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

/// Evaluates the first equilibrium-closure residual without mutating the live
/// cell. `output` follows `State.packCell` ordering and native field units.
pub fn evaluateCellResidual(
    allocator: std.mem.Allocator,
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    options: Options,
    output: []f64,
) !f64 {
    try validateOptions(options);
    if (cell_index >= state.cell_count)
        return error.ChemistryCellIndexOutOfBounds;
    if (output.len != chemistry.State.packedComponentCount())
        return error.ChemistryVectorSizeMismatch;
    var workspace = try Workspace.init(allocator);
    defer workspace.deinit();
    try state.packCell(cell_index, workspace.current);
    try evaluateGlobalResidualAt(
        &workspace.scratch,
        workspace.current,
        equilibriumClosureParameters(parameters),
        output,
    );
    return scaledNorm(workspace.current, output, options);
}

fn solveCellWithWorkspaceAndTrace(
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

fn validateAqueousMolarity(
    state: *const chemistry.State,
    cell_index: usize,
) !void {
    const water_mol_per_m3 = state.water_mol_per_m3[cell_index];
    if (!std.math.isFinite(water_mol_per_m3) or water_mol_per_m3 <= 0)
        return error.InvalidChemistryWaterState;
    inline for (
        @typeInfo(aqueous_network.State).@"struct".fields,
        0..,
    ) |field, component_index| {
        const concentration = @field(state.aqueous[cell_index], field.name);
        if (!std.math.isFinite(concentration) or concentration < 0)
            return error.NonFiniteSoluteReactionState;
        if (concentration > water_mol_per_m3) {
            std.log.warn(
                "SOLUTE aqueous concentration exceeds water molarity: cell={d} packed_component={d} name=aqueous.{s} concentration_mol_per_m3={e} water_mol_per_m3={e}",
                .{
                    cell_index,
                    component_index,
                    field.name,
                    concentration,
                    water_mol_per_m3,
                },
            );
            return error.SoluteConcentrationExceedsWaterMolarity;
        }
    }
}

fn solveEquilibriumWithWorkspace(
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

fn logTerminalReactionDecomposition(
    scratch: *chemistry.State,
    current: []const f64,
    parameters: chemistry.ReactionParameters,
) void {
    scratch.unpackCell(0, current) catch |err| {
        std.log.debug(
            "SOLUTE terminal decomposition unavailable: stage=unpack error={s}",
            .{@errorName(err)},
        );
        return;
    };
    const terminal_decomposition =
        diagnoseCell(scratch, 0, parameters) catch |err| {
            std.log.debug(
                "SOLUTE terminal decomposition unavailable: stage=reaction_network error={s}",
                .{@errorName(err)},
            );
            return;
        };
    std.log.debug(
        "SOLUTE terminal hydrogen decomposition: aqueous={e} phosphate={e} exchange={e} carboxyl={e} geochemistry={e} assembled={e}",
        .{
            terminal_decomposition.aqueous_hydrogen_mol_per_m3,
            terminal_decomposition.phosphate_hydrogen_mol_per_m3,
            terminal_decomposition.cation_exchange_hydrogen_mol_per_m3,
            terminal_decomposition.carboxyl_hydrogen_mol_per_m3,
            terminal_decomposition.geochemistry_hydrogen_mol_per_m3,
            terminal_decomposition.assembled_hydrogen_mol_per_m3,
        },
    );
    const hpo4 =
        diagnoseNonBandHpo4(scratch, 0, parameters) catch |err| {
            std.log.warn(
                "SOLUTE terminal decomposition unavailable: stage=non_band_hpo4 error={s}",
                .{@errorName(err)},
            );
            return;
        };
    std.log.debug(
        "SOLUTE terminal non-band HPO4 decomposition: PO4_protonation={e} HPO4_protonation={e} metal_pairing={e} surface_adsorption={e} assembled={e}",
        .{
            hpo4.po4_protonation_mol_p_per_m3,
            hpo4.hpo4_protonation_mol_p_per_m3,
            hpo4.metal_pairing_mol_p_per_m3,
            hpo4.surface_adsorption_mol_p_per_m3,
            hpo4.assembled_change_mol_p_per_m3,
        },
    );
}

fn exhaustsLargestResidual(
    current: []const f64,
    residual: []const f64,
    candidate: []const f64,
    options: Options,
) bool {
    var worst_index: usize = 0;
    var worst_scaled: f64 = 0;
    for (current, residual, 0..) |value, change, index| {
        const scaled_residual = @abs(change) / residualScale(value, options);
        if (scaled_residual > worst_scaled) {
            worst_scaled = scaled_residual;
            worst_index = index;
        }
    }
    return residual[worst_index] < 0 and
        candidate[worst_index] <= options.absolute_tolerance;
}

test "active set advances only when the largest scaled residual exhausts its inventory" {
    const options: Options = .{
        .absolute_tolerance = 1.0e-8,
        .relative_tolerance = 0,
    };
    const current = [_]f64{ 1, 0.1, 2 };
    const residual = [_]f64{ -0.01, -0.1, 0.001 };

    try std.testing.expect(exhaustsLargestResidual(
        &current,
        &residual,
        &.{ 0.99, 0, 2.001 },
        options,
    ));
    try std.testing.expect(!exhaustsLargestResidual(
        &current,
        &residual,
        &.{ 0.99, 0.01, 2.001 },
        options,
    ));
    try std.testing.expect(!exhaustsLargestResidual(
        &current,
        &.{ -0.01, 0.1, 0.001 },
        &.{ 0.99, 0, 2.001 },
        options,
    ));
}

test "aqueous calcium extent conserves calcium ligand and charge" {
    const vector = try std.testing.allocator.alloc(
        f64,
        chemistry.State.packedComponentCount(),
    );
    defer std.testing.allocator.free(vector);
    @memset(vector, 1);
    vector[aqueousPackedIndex("calcium")] = 0.8;
    vector[aqueousPackedIndex("carbonate")] = 0.6;
    vector[aqueousPackedIndex("calcium_carbonate")] = 0.2;
    const calcium_before = vector[aqueousPackedIndex("calcium")] +
        vector[aqueousPackedIndex("calcium_carbonate")];
    const carbonate_before = vector[aqueousPackedIndex("carbonate")] +
        vector[aqueousPackedIndex("calcium_carbonate")];
    const charge_before =
        2 * vector[aqueousPackedIndex("calcium")] -
        2 * vector[aqueousPackedIndex("carbonate")];

    try std.testing.expect(applyPhosphateExtent(
        vector,
        .aqueous_calcium_carbonate_pairing,
        0.25,
        undefined,
    ));
    try std.testing.expectEqual(
        calcium_before,
        vector[aqueousPackedIndex("calcium")] +
            vector[aqueousPackedIndex("calcium_carbonate")],
    );
    try std.testing.expectEqual(
        carbonate_before,
        vector[aqueousPackedIndex("carbonate")] +
            vector[aqueousPackedIndex("calcium_carbonate")],
    );
    try std.testing.expectEqual(
        charge_before,
        2 * vector[aqueousPackedIndex("calcium")] -
            2 * vector[aqueousPackedIndex("carbonate")],
    );
    try std.testing.expectEqual(
        @as(f64, 0.35),
        maximumPhosphateExtent(
            vector,
            .aqueous_calcium_carbonate_pairing,
            1,
            undefined,
        ),
    );
}

const ReactionSpanExtentBounds = struct {
    negative_native_extent: f64,
    positive_native_extent: f64,
};

/// Rank-revealing semismooth Newton step in the complete conservative
/// equilibrium-reaction span. Native extents are normalized by their local
/// inventory boxes before QR so mol/m3, mol/Mg, and the dimensionless Gapon
/// direction can share one numerically meaningful system.
fn tryFullNetworkReactionCandidate(
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
fn tryRetainedComplementarityCandidate(
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

fn tryFullNetworkComplementarityCandidate(
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

fn complementarityColumnDerivative(
    branch: i8,
    negative_derivative: f64,
    positive_derivative: f64,
) f64 {
    return switch (branch) {
        -1 => negative_derivative,
        0 => 0.5 *
            (negative_derivative + positive_derivative),
        1 => positive_derivative,
        else => unreachable,
    };
}

test "zero-extent complementarity uses a generalized derivative" {
    try std.testing.expectEqual(
        @as(f64, -7),
        complementarityColumnDerivative(-1, -7, 11),
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        complementarityColumnDerivative(0, -7, 11),
    );
    try std.testing.expectEqual(
        @as(f64, 11),
        complementarityColumnDerivative(1, -7, 11),
    );
}

fn tryAcceptReactionExtentLineSearch(
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

fn tryAcceptExactReactionExtentCandidate(
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

fn refineComplementarityDirectionalJacobian(
    workspace: *Workspace,
    inputs: ComplementaritySearchInputs,
    transformations: chemistry.CellTransformations,
    negative_jacobian: []f64,
    positive_jacobian: []f64,
    selected_jacobian: []f64,
    probe_state: []f64,
    probe_residual: []f64,
) bool {
    var solution_norm_squared: f64 = 0;
    for (workspace.reaction_span_solution[0..inputs.column_count]) |extent|
        solution_norm_squared += extent * extent;
    if (!std.math.isFinite(solution_norm_squared) or
        solution_norm_squared <= std.math.floatEps(f64))
    {
        return false;
    }
    const directional_fraction =
        std.math.cbrt(std.math.floatEps(f64));
    _ = transformedVectorAdmissible(
        inputs.scratch,
        inputs.current,
        transformations,
        inputs.parameters,
        directional_fraction,
        probe_state,
    ) catch return false;
    evaluateGlobalResidualAt(
        inputs.scratch,
        probe_state,
        inputs.parameters,
        probe_residual,
    ) catch return false;

    var maximum_correction: f64 = 0;
    for (0..inputs.current.len) |row| {
        const scale = residualScale(
            inputs.current[row],
            inputs.options,
        );
        const weight = reactionSpanRowWeight(
            inputs.global_residual[row] / scale,
            inputs.current_norm,
        );
        const actual_directional_derivative =
            (probe_residual[row] -
                inputs.global_residual[row]) /
            directional_fraction / scale * weight;
        var predicted_directional_derivative: f64 = 0;
        for (
            workspace.reaction_span_solution[0..inputs.column_count],
            0..,
        ) |extent, column| {
            predicted_directional_derivative +=
                selected_jacobian[
                    row * inputs.column_count + column
                ] * extent;
        }
        const difference =
            actual_directional_derivative -
            predicted_directional_derivative;
        if (!std.math.isFinite(difference)) return false;
        maximum_correction =
            @max(maximum_correction, @abs(difference));
        for (
            workspace.reaction_span_solution[0..inputs.column_count],
            0..,
        ) |extent, column| {
            if (extent == 0) continue;
            const correction =
                difference * extent / solution_norm_squared;
            const index = row * inputs.column_count + column;
            selected_jacobian[index] += correction;
            switch (workspace.reaction_span_branch_states[column]) {
                -1 => negative_jacobian[index] += correction,
                1 => positive_jacobian[index] += correction,
                0 => {},
                else => unreachable,
            }
        }
    }
    return maximum_correction >
        64 * std.math.floatEps(f64);
}

fn reactionSpanRowWeight(
    normalized_residual: f64,
    current_norm: f64,
) f64 {
    return @max(
        @sqrt(std.math.floatEps(f64)),
        @abs(normalized_residual) / current_norm,
    );
}

fn reactionSpanPredictedNorm(
    workspace: *const Workspace,
    current: []const f64,
    global_residual: []const f64,
    options: Options,
    current_norm: f64,
    column_count: usize,
    step_fraction: f64,
) f64 {
    var predicted_norm: f64 = 0;
    for (global_residual, current, 0..) |value, state_value, row| {
        const scale = residualScale(state_value, options);
        const normalized_residual = value / scale;
        const weight =
            reactionSpanRowWeight(normalized_residual, current_norm);
        var predicted_weighted = normalized_residual * weight;
        for (
            workspace.reaction_span_solution[0..column_count],
            0..,
        ) |extent, column| {
            predicted_weighted +=
                workspace.reaction_span_jacobian[
                    row * column_count + column
                ] * extent * step_fraction;
        }
        predicted_norm = @max(
            predicted_norm,
            @abs(predicted_weighted / weight),
        );
    }
    return predicted_norm;
}

fn captureFullNetworkDirectionalComparison(
    trace: *SolverTrace,
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    target: []const f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
    column_count: usize,
    inventory_fraction: f64,
    current_transformations: chemistry.CellTransformations,
    monovalent_activity_coefficient: f64,
    diagnostic: ?*IterationDiagnostic,
) void {
    const first_fraction = phosphateTrustRegionFraction(
        current,
        target,
        options,
        parameters.phosphate_kinetics.substrate_limit_fraction,
    );
    const first_norm = evaluateCandidateResidualAtFraction(
        scratch,
        current,
        target,
        accepted_state,
        residual_work,
        parameters,
        options,
        first_fraction,
    ) catch std.math.inf(f64);
    const directional_fraction = @min(
        first_fraction,
        std.math.cbrt(std.math.floatEps(f64)),
    );
    _ = evaluateCandidateResidualAtFraction(
        scratch,
        current,
        target,
        accepted_state,
        residual_work,
        parameters,
        options,
        directional_fraction,
    ) catch return;

    @memset(trace.full_network_current_rates, std.math.nan(f64));
    _ = evaluateAt(scratch, current, parameters) catch {};
    reaction_span.evaluateRates(
        scratch,
        0,
        parameters,
        trace.full_network_current_rates,
    ) catch {};
    @memset(trace.full_network_directional_rates, std.math.nan(f64));
    _ = evaluateAt(scratch, accepted_state, parameters) catch {};
    reaction_span.evaluateRates(
        scratch,
        0,
        parameters,
        trace.full_network_directional_rates,
    ) catch {};

    @memcpy(trace.full_network_base_residual, global_residual);
    @memcpy(trace.full_network_realized_residual, residual_work);
    for (global_residual, current, 0..) |value, state_value, row| {
        const scale = residualScale(state_value, options);
        const normalized_residual = value / scale;
        const weight =
            reactionSpanRowWeight(normalized_residual, current_norm);
        var predicted_weighted = normalized_residual * weight;
        for (
            workspace.reaction_span_solution[0..column_count],
            0..,
        ) |extent, column| {
            predicted_weighted +=
                workspace.reaction_span_jacobian[
                    row * column_count + column
                ] * extent * inventory_fraction * directional_fraction;
        }
        trace.full_network_predicted_residual[row] =
            predicted_weighted / weight * scale;
    }
    const limiting_component_index = if (diagnostic) |entry|
        entry.limiting_component_index
    else
        largestScaledResidualIndex(current, global_residual, options);
    const limiting_scale =
        residualScale(current[limiting_component_index], options);
    const limiting_normalized_residual =
        global_residual[limiting_component_index] / limiting_scale;
    const limiting_weight = reactionSpanRowWeight(
        limiting_normalized_residual,
        current_norm,
    );
    const side_derivative_inputs = ReactionSideDerivativeInputs{
        .scratch = scratch,
        .current = current,
        .global_residual = global_residual,
        .probe_state = accepted_state,
        .probe_residual = residual_work,
        .current_transformations = current_transformations,
        .parameters = parameters,
        .options = options,
        .monovalent_activity_coefficient = monovalent_activity_coefficient,
        .limiting_component_index = limiting_component_index,
    };
    for (0..column_count) |column| {
        const reaction =
            workspace.reaction_span_active_reactions[column];
        trace.full_network_reactions[column] = .{
            .reaction_index = reaction,
            .selection_rate = workspace.reaction_span_rates[reaction],
            .current_rate = trace.full_network_current_rates[reaction],
            .directional_rate = trace.full_network_directional_rates[reaction],
            .native_extent_scale = workspace.reaction_span_extent_scales[column],
            .normalized_lower_bound = workspace.reaction_span_lower_bounds[column],
            .normalized_upper_bound = workspace.reaction_span_upper_bounds[column],
            .normalized_solution = workspace.reaction_span_solution[column],
            .limiting_normalized_residual_derivative = workspace.reaction_span_jacobian[
                limiting_component_index * column_count + column
            ] / limiting_weight,
            .solution_side_limiting_normalized_residual_derivative = reactionSpanSolutionSideLimitingDerivative(
                side_derivative_inputs,
                reaction,
                workspace.reaction_span_extent_scales[column],
                workspace.reaction_span_lower_bounds[column],
                workspace.reaction_span_upper_bounds[column],
                workspace.reaction_span_solution[column],
            ) catch std.math.nan(f64),
        };
    }
    searchReactionSpanComplementarity(
        trace,
        workspace,
        scratch,
        current,
        global_residual,
        current_transformations,
        accepted_state,
        residual_work,
        parameters,
        options,
        current_norm,
        monovalent_activity_coefficient,
        column_count,
    ) catch {};
    trace.full_network_comparison_valid = true;
    trace.full_network_limiting_component_index =
        limiting_component_index;
    trace.full_network_reaction_count = column_count;
    trace.full_network_directional_fraction = directional_fraction;
    trace.full_network_comparison_inventory_fraction =
        inventory_fraction;
    trace.full_network_first_trial_fraction = first_fraction;
    trace.full_network_first_trial_maximum_scaled_residual = first_norm;
    if (diagnostic) |entry| {
        trace.full_network_comparison_closure = entry.closure_index;
        trace.full_network_comparison_iteration = entry.iteration;
    }
}

const ComplementaritySearchInputs = struct {
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    probe_state: []f64,
    probe_residual: []f64,
    current_transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
    monovalent_activity_coefficient: f64,
    column_count: usize,
};

fn searchReactionSpanComplementarity(
    trace: *SolverTrace,
    workspace: *Workspace,
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    current_transformations: chemistry.CellTransformations,
    probe_state: []f64,
    probe_residual: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    current_norm: f64,
    monovalent_activity_coefficient: f64,
    column_count: usize,
) !void {
    trace.full_network_complementarity_search_attempted = true;
    var ambiguous_count: usize = 0;
    for (
        workspace.reaction_span_solution[0..column_count],
        0..,
    ) |solution, column| {
        if (!reactionSpanExtentIsSignificant(
            workspace,
            column,
            solution,
        )) continue;
        const reaction =
            workspace.reaction_span_active_reactions[column];
        const rate = workspace.reaction_span_rates[reaction];
        if ((solution < 0) == (rate < 0)) continue;
        trace.full_network_complementarity_ambiguous_columns[
            ambiguous_count
        ] = column;
        ambiguous_count += 1;
    }
    trace.full_network_complementarity_ambiguous_count =
        ambiguous_count;
    if (ambiguous_count == 0 or ambiguous_count > 16) return;

    const matrix_count = current.len * column_count;
    const selected_jacobian = try trace.allocator.dupe(
        f64,
        workspace.reaction_span_jacobian[0..matrix_count],
    );
    defer trace.allocator.free(selected_jacobian);
    const opposite_jacobian = try trace.allocator.dupe(
        f64,
        selected_jacobian,
    );
    defer trace.allocator.free(opposite_jacobian);
    const original_solution = try trace.allocator.dupe(
        f64,
        workspace.reaction_span_solution[0..column_count],
    );
    defer trace.allocator.free(original_solution);
    const inputs = ComplementaritySearchInputs{
        .scratch = scratch,
        .current = current,
        .global_residual = global_residual,
        .probe_state = probe_state,
        .probe_residual = probe_residual,
        .current_transformations = current_transformations,
        .parameters = parameters,
        .options = options,
        .current_norm = current_norm,
        .monovalent_activity_coefficient = monovalent_activity_coefficient,
        .column_count = column_count,
    };
    for (
        trace.full_network_complementarity_ambiguous_columns[0..ambiguous_count],
    ) |column| {
        const reaction =
            workspace.reaction_span_active_reactions[column];
        const selected_rate =
            workspace.reaction_span_rates[reaction];
        try evaluateReactionSpanDerivativeColumn(
            inputs,
            workspace,
            column,
            if (selected_rate < 0) 1 else -1,
            opposite_jacobian,
        );
    }

    const combination_count_u64 =
        @as(u64, 1) << @intCast(ambiguous_count);
    const combination_count: usize =
        @intCast(combination_count_u64);
    trace.full_network_complementarity_combination_count =
        combination_count;
    const pattern_mismatch_counts = try trace.allocator.alloc(
        usize,
        combination_count,
    );
    defer trace.allocator.free(pattern_mismatch_counts);
    @memset(pattern_mismatch_counts, std.math.maxInt(usize));
    const pattern_mismatch_columns = try trace.allocator.alloc(
        usize,
        combination_count,
    );
    defer trace.allocator.free(pattern_mismatch_columns);
    const pattern_mismatch_solutions = try trace.allocator.alloc(
        f64,
        combination_count,
    );
    defer trace.allocator.free(pattern_mismatch_solutions);
    const mismatch_work = try trace.allocator.alloc(
        usize,
        column_count,
    );
    defer trace.allocator.free(mismatch_work);
    var mask: u64 = 0;
    while (mask < combination_count_u64) : (mask += 1) {
        @memcpy(
            workspace.reaction_span_jacobian[0..matrix_count],
            selected_jacobian,
        );
        for (
            trace.full_network_complementarity_ambiguous_columns[0..ambiguous_count],
            0..,
        ) |column, bit| {
            if (mask & (@as(u64, 1) << @intCast(bit)) == 0)
                continue;
            for (0..current.len) |row|
                workspace.reaction_span_jacobian[
                    row * column_count + column
                ] = opposite_jacobian[
                    row * column_count + column
                ];
        }
        if (!solveBoundedReactionSpan(
            workspace,
            current.len,
            column_count,
        )) continue;
        const mismatch_count = reactionSpanDerivativeSideMismatchCount(
            workspace,
            column_count,
            trace.full_network_complementarity_ambiguous_columns[0..ambiguous_count],
            mask,
            mismatch_work,
        );
        const pattern_index: usize = @intCast(mask);
        pattern_mismatch_counts[pattern_index] =
            mismatch_count;
        if (mismatch_count == 1) {
            const mismatch_column = mismatch_work[0];
            pattern_mismatch_columns[pattern_index] =
                mismatch_column;
            pattern_mismatch_solutions[pattern_index] =
                workspace.reaction_span_solution[mismatch_column];
        }
        if (mismatch_count <
            trace.full_network_complementarity_minimum_mismatch_count)
        {
            trace.full_network_complementarity_minimum_mismatch_count =
                mismatch_count;
            trace.full_network_complementarity_minimum_mismatch_mask =
                mask;
            @memcpy(
                trace.full_network_complementarity_mismatch_columns[0..mismatch_count],
                mismatch_work[0..mismatch_count],
            );
        }
        if (mismatch_count != 0) continue;
        trace.full_network_complementarity_consistent_count += 1;
        const merit =
            try evaluateComplementarityCandidate(inputs, workspace) orelse
            continue;
        trace.full_network_complementarity_exact_descent_count += 1;
        if (merit.exact <
            trace.full_network_complementarity_best_exact_merit)
        {
            trace.full_network_complementarity_best_mask = mask;
            trace.full_network_complementarity_best_predicted_merit =
                merit.predicted;
            trace.full_network_complementarity_best_exact_merit =
                merit.exact;
        }
    }
    try searchComplementarityKinkBlends(
        trace,
        workspace,
        inputs,
        selected_jacobian,
        opposite_jacobian,
        pattern_mismatch_counts,
        pattern_mismatch_columns,
        pattern_mismatch_solutions,
    );
    @memcpy(
        workspace.reaction_span_jacobian[0..matrix_count],
        selected_jacobian,
    );
    @memcpy(
        workspace.reaction_span_solution[0..column_count],
        original_solution,
    );
}

const ComplementarityMerit = struct {
    predicted: f64,
    exact: f64,
};

fn evaluateComplementarityCandidate(
    inputs: ComplementaritySearchInputs,
    workspace: *Workspace,
) !?ComplementarityMerit {
    var transformations =
        reaction_span.zeroTransformations(inputs.parameters);
    var has_nonzero_extent = false;
    for (
        workspace.reaction_span_solution[0..inputs.column_count],
        0..,
    ) |normalized_extent, column| {
        if (!std.math.isFinite(normalized_extent))
            return error.NonFiniteSoluteReactionExtent;
        if (normalized_extent == 0) continue;
        try reaction_span.addReactionExtent(
            &transformations,
            workspace.reaction_span_active_reactions[column],
            normalized_extent *
                workspace.reaction_span_extent_scales[column],
            inputs.current_transformations,
            inputs.parameters,
        );
        has_nonzero_extent = true;
    }
    if (!has_nonzero_extent) return null;
    const inventory_fraction = transformedVectorAdmissible(
        inputs.scratch,
        inputs.current,
        transformations,
        inputs.parameters,
        1,
        inputs.probe_state,
    ) catch return null;
    const predicted_norm = reactionSpanPredictedNorm(
        workspace,
        inputs.current,
        inputs.global_residual,
        inputs.options,
        inputs.current_norm,
        inputs.column_count,
        inventory_fraction,
    );
    if (!std.math.isFinite(predicted_norm)) return null;
    if (!try tryAcceptAndersonCandidate(
        inputs.scratch,
        inputs.current,
        inputs.probe_state,
        workspace.rollback_state,
        inputs.probe_residual,
        inputs.parameters,
        inputs.options,
        inputs.current_norm,
    )) return null;
    const exact_norm = try scaledNorm(
        workspace.rollback_state,
        inputs.probe_residual,
        inputs.options,
    );
    if (inputs.current_norm - exact_norm <
        1.0e-6 * @max(1.0, inputs.current_norm))
    {
        return null;
    }
    return .{
        .predicted = predicted_norm,
        .exact = exact_norm,
    };
}

fn searchComplementarityKinkBlends(
    trace: *SolverTrace,
    workspace: *Workspace,
    inputs: ComplementaritySearchInputs,
    selected_jacobian: []const f64,
    opposite_jacobian: []const f64,
    pattern_mismatch_counts: []const usize,
    pattern_mismatch_columns: []const usize,
    pattern_mismatch_solutions: []const f64,
) !void {
    if (trace.full_network_complementarity_minimum_mismatch_count != 1)
        return;
    const ambiguous_columns =
        trace.full_network_complementarity_ambiguous_columns[0..trace.full_network_complementarity_ambiguous_count];
    var best_base_mask: ?u64 = null;
    var best_blend_column: usize = 0;
    var best_selected_solution: f64 = 0;
    var best_opposite_solution: f64 = 0;
    var best_endpoint_score = std.math.inf(f64);
    for (pattern_mismatch_counts, 0..) |base_count, base_index| {
        if (base_count != 1) continue;
        const base_mask: u64 = @intCast(base_index);
        for (ambiguous_columns, 0..) |blend_column, bit| {
            const bit_mask = @as(u64, 1) << @intCast(bit);
            if (base_mask & bit_mask != 0) continue;
            const other_index: usize =
                @intCast(base_mask | bit_mask);
            if (pattern_mismatch_counts[other_index] != 1 or
                pattern_mismatch_columns[base_index] != blend_column or
                pattern_mismatch_columns[other_index] != blend_column)
            {
                continue;
            }
            const endpoint_selected_solution =
                pattern_mismatch_solutions[base_index];
            const endpoint_opposite_solution =
                pattern_mismatch_solutions[other_index];
            if ((endpoint_selected_solution < 0) ==
                (endpoint_opposite_solution < 0))
            {
                continue;
            }
            const endpoint_score = @max(
                @abs(endpoint_selected_solution),
                @abs(endpoint_opposite_solution),
            );
            if (endpoint_score >= best_endpoint_score) continue;
            best_endpoint_score = endpoint_score;
            best_base_mask = base_mask;
            best_blend_column = blend_column;
            best_selected_solution = endpoint_selected_solution;
            best_opposite_solution = endpoint_opposite_solution;
        }
    }
    const base_mask = best_base_mask orelse return;
    trace.full_network_complementarity_blend_attempted = true;
    trace.full_network_complementarity_blend_column =
        best_blend_column;

    var selected_fraction: f64 = 0;
    var opposite_fraction: f64 = 1;
    var selected_solution = best_selected_solution;
    var opposite_solution = best_opposite_solution;
    var blend_fraction: f64 = 0.5;
    var blend_solution: f64 = std.math.nan(f64);
    var iteration: u8 = 0;
    while (iteration < 80) : (iteration += 1) {
        blend_fraction = if (opposite_solution != selected_solution)
            selected_fraction -
                selected_solution *
                    (opposite_fraction - selected_fraction) /
                    (opposite_solution - selected_solution)
        else
            selected_fraction +
                0.5 * (opposite_fraction - selected_fraction);
        if (!std.math.isFinite(blend_fraction) or
            blend_fraction <= selected_fraction or
            blend_fraction >= opposite_fraction)
        {
            blend_fraction = selected_fraction +
                0.5 * (opposite_fraction - selected_fraction);
        }
        blend_solution = solveComplementarityBlend(
            workspace,
            inputs,
            selected_jacobian,
            opposite_jacobian,
            ambiguous_columns,
            base_mask,
            best_blend_column,
            blend_fraction,
        ) orelse return;
        if (!reactionSpanExtentIsSignificant(
            workspace,
            best_blend_column,
            blend_solution,
        )) break;
        if ((blend_solution < 0) ==
            (selected_solution < 0))
        {
            selected_fraction = blend_fraction;
            selected_solution = blend_solution;
        } else {
            opposite_fraction = blend_fraction;
            opposite_solution = blend_solution;
        }
    }
    trace.full_network_complementarity_blend_fraction =
        blend_fraction;
    trace.full_network_complementarity_blend_solution =
        blend_solution;
    if (reactionSpanExtentIsSignificant(
        workspace,
        best_blend_column,
        blend_solution,
    )) return;
    if (reactionSpanDerivativeSideMismatchCount(
        workspace,
        inputs.column_count,
        ambiguous_columns,
        base_mask,
        null,
    ) != 0) return;
    const merit =
        try evaluateComplementarityCandidate(inputs, workspace) orelse
        return;
    trace.full_network_complementarity_blend_predicted_merit =
        merit.predicted;
    trace.full_network_complementarity_blend_exact_merit =
        merit.exact;
}

fn solveComplementarityBlend(
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

fn evaluateReactionSpanDerivativeColumn(
    inputs: ComplementaritySearchInputs,
    workspace: *const Workspace,
    column: usize,
    direction: i8,
    output_jacobian: []f64,
) !void {
    const available_normalized = if (direction < 0)
        -workspace.reaction_span_lower_bounds[column] /
            inputs.options.maximum_newton_fraction
    else
        workspace.reaction_span_upper_bounds[column] /
            inputs.options.maximum_newton_fraction;
    const probe_magnitude = @min(
        @sqrt(std.math.floatEps(f64)),
        0.125 * available_normalized,
    );
    if (!std.math.isFinite(probe_magnitude) or
        probe_magnitude <= 64 * std.math.floatEps(f64))
    {
        return error.NoAdmissibleComplementarityProbe;
    }
    const normalized_probe =
        @as(f64, @floatFromInt(direction)) * probe_magnitude;
    var transformations =
        reaction_span.zeroTransformations(inputs.parameters);
    try reaction_span.addReactionExtent(
        &transformations,
        workspace.reaction_span_active_reactions[column],
        normalized_probe *
            workspace.reaction_span_extent_scales[column],
        inputs.current_transformations,
        inputs.parameters,
    );
    try transformedVector(
        inputs.scratch,
        inputs.current,
        transformations,
        inputs.monovalent_activity_coefficient,
        inputs.parameters.water_activity_product_mol2_per_m6,
        1,
        inputs.probe_state,
    );
    try evaluateGlobalResidualAt(
        inputs.scratch,
        inputs.probe_state,
        inputs.parameters,
        inputs.probe_residual,
    );
    for (0..inputs.current.len) |row| {
        const scale = residualScale(
            inputs.current[row],
            inputs.options,
        );
        const active_weight = reactionSpanRowWeight(
            inputs.global_residual[row] / scale,
            inputs.current_norm,
        );
        output_jacobian[
            row * inputs.column_count + column
        ] = (inputs.probe_residual[row] -
            inputs.global_residual[row]) / normalized_probe /
            scale * active_weight;
    }
}

fn reactionSpanDerivativeSideMismatchCount(
    workspace: *const Workspace,
    column_count: usize,
    ambiguous_columns: []const usize,
    opposite_mask: u64,
    mismatch_columns: ?[]usize,
) usize {
    var mismatch_count: usize = 0;
    for (
        workspace.reaction_span_solution[0..column_count],
        0..,
    ) |solution, column| {
        if (!reactionSpanExtentIsSignificant(
            workspace,
            column,
            solution,
        )) continue;
        const reaction =
            workspace.reaction_span_active_reactions[column];
        const selected_rate =
            workspace.reaction_span_rates[reaction];
        var direction: i8 = if (selected_rate < 0) -1 else 1;
        for (ambiguous_columns, 0..) |ambiguous_column, bit| {
            if (column != ambiguous_column) continue;
            if (opposite_mask &
                (@as(u64, 1) << @intCast(bit)) != 0)
            {
                direction = -direction;
            }
            break;
        }
        if ((solution < 0) != (direction < 0)) {
            if (mismatch_columns) |columns|
                columns[mismatch_count] = column;
            mismatch_count += 1;
        }
    }
    return mismatch_count;
}

fn reactionSpanExtentIsSignificant(
    workspace: *const Workspace,
    column: usize,
    solution: f64,
) bool {
    const bound_scale = @max(
        1.0,
        @max(
            @abs(workspace.reaction_span_lower_bounds[column]),
            @abs(workspace.reaction_span_upper_bounds[column]),
        ),
    );
    return @abs(solution) >
        64 * std.math.floatEps(f64) * bound_scale;
}

const ReactionSideDerivativeInputs = struct {
    scratch: *chemistry.State,
    current: []const f64,
    global_residual: []const f64,
    probe_state: []f64,
    probe_residual: []f64,
    current_transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
    options: Options,
    monovalent_activity_coefficient: f64,
    limiting_component_index: usize,
};

fn reactionSpanSolutionSideLimitingDerivative(
    inputs: ReactionSideDerivativeInputs,
    reaction: usize,
    extent_scale: f64,
    normalized_lower_bound: f64,
    normalized_upper_bound: f64,
    normalized_solution: f64,
) !f64 {
    if (normalized_solution == 0)
        return std.math.nan(f64);
    const available_normalized = if (normalized_solution < 0)
        -normalized_lower_bound / inputs.options.maximum_newton_fraction
    else
        normalized_upper_bound / inputs.options.maximum_newton_fraction;
    const probe_magnitude = @min(
        std.math.cbrt(std.math.floatEps(f64)),
        0.125 * available_normalized,
    );
    if (!std.math.isFinite(probe_magnitude) or
        probe_magnitude <= 64 * std.math.floatEps(f64))
    {
        return std.math.nan(f64);
    }
    const normalized_probe =
        if (normalized_solution < 0)
            -probe_magnitude
        else
            probe_magnitude;
    var transformations =
        reaction_span.zeroTransformations(inputs.parameters);
    try reaction_span.addReactionExtent(
        &transformations,
        reaction,
        normalized_probe * extent_scale,
        inputs.current_transformations,
        inputs.parameters,
    );
    try transformedVector(
        inputs.scratch,
        inputs.current,
        transformations,
        inputs.monovalent_activity_coefficient,
        inputs.parameters.water_activity_product_mol2_per_m6,
        1,
        inputs.probe_state,
    );
    try evaluateGlobalResidualAt(
        inputs.scratch,
        inputs.probe_state,
        inputs.parameters,
        inputs.probe_residual,
    );
    const row = inputs.limiting_component_index;
    const scale = residualScale(
        inputs.current[row],
        inputs.options,
    );
    return (inputs.probe_residual[row] -
        inputs.global_residual[row]) / normalized_probe / scale;
}

fn reactionSpanExtentBounds(
    scratch: *chemistry.State,
    current: []const f64,
    reaction: usize,
    current_transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
    rate_magnitude: f64,
    monovalent_activity_coefficient: f64,
    output: []f64,
) !ReactionSpanExtentBounds {
    if (reaction == reaction_span.carboxyl_reaction_index) {
        const carboxyl_index =
            @typeInfo(aqueous_network.State).@"struct".fields.len +
            2 * @typeInfo(phosphate_network.State).@"struct".fields.len +
            @typeInfo(cation_exchange.Cations).@"struct".fields.len;
        return .{
            .negative_native_extent = current[carboxyl_index],
            .positive_native_extent = @max(
                0,
                parameters.total_carboxyl_sites_mol_per_megagram -
                    current[carboxyl_index],
            ),
        };
    }
    return .{
        .negative_native_extent = try maximumReactionSpanExtent(
            scratch,
            current,
            reaction,
            -1,
            current_transformations,
            parameters,
            rate_magnitude,
            monovalent_activity_coefficient,
            output,
        ),
        .positive_native_extent = try maximumReactionSpanExtent(
            scratch,
            current,
            reaction,
            1,
            current_transformations,
            parameters,
            rate_magnitude,
            monovalent_activity_coefficient,
            output,
        ),
    };
}

/// Locates the inventory face for one native conservative extent using the
/// same atomic state transaction as production. This intentionally treats
/// H+/OH- as dependent Kw coordinates; `transformedVector` projects them
/// exactly and all remaining negative inventories reject the trial.
fn maximumReactionSpanExtent(
    scratch: *chemistry.State,
    current: []const f64,
    reaction: usize,
    direction: f64,
    current_transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
    rate_magnitude: f64,
    monovalent_activity_coefficient: f64,
    output: []f64,
) !f64 {
    var trial = @max(
        rate_magnitude,
        64 * std.math.floatEps(f64),
    );
    var upper: ?f64 = null;
    var lower: f64 = 0;
    var contraction: u8 = 0;
    while (contraction < 48) : (contraction += 1) {
        if (reactionSpanExtentAdmissible(
            scratch,
            current,
            reaction,
            direction * trial,
            current_transformations,
            parameters,
            monovalent_activity_coefficient,
            output,
        )) {
            lower = trial;
            break;
        }
        upper = trial;
        trial *= 0.5;
    }
    if (lower == 0) return 0;

    if (upper == null) {
        var expansion: u8 = 0;
        while (expansion < 48) : (expansion += 1) {
            const next = lower * 2;
            if (!std.math.isFinite(next)) break;
            if (!reactionSpanExtentAdmissible(
                scratch,
                current,
                reaction,
                direction * next,
                current_transformations,
                parameters,
                monovalent_activity_coefficient,
                output,
            )) {
                upper = next;
                break;
            }
            lower = next;
        }
    }
    var high = upper orelse return lower;
    var search: u8 = 0;
    while (search < 36) : (search += 1) {
        const middle = lower + 0.5 * (high - lower);
        if (reactionSpanExtentAdmissible(
            scratch,
            current,
            reaction,
            direction * middle,
            current_transformations,
            parameters,
            monovalent_activity_coefficient,
            output,
        )) {
            lower = middle;
        } else {
            high = middle;
        }
        if (high - lower <= 64 * std.math.floatEps(f64) *
            @max(1.0, lower))
        {
            break;
        }
    }
    return lower;
}

fn reactionSpanExtentAdmissible(
    scratch: *chemistry.State,
    current: []const f64,
    reaction: usize,
    native_extent: f64,
    current_transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
    monovalent_activity_coefficient: f64,
    output: []f64,
) bool {
    var transformations = reaction_span.zeroTransformations(parameters);
    reaction_span.addReactionExtent(
        &transformations,
        reaction,
        native_extent,
        current_transformations,
        parameters,
    ) catch return false;
    transformedVector(
        scratch,
        current,
        transformations,
        monovalent_activity_coefficient,
        parameters.water_activity_product_mol2_per_m6,
        1,
        output,
    ) catch return false;
    return true;
}

fn tryPhosphateExtentCandidate(
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

const PhosphateExtentBounds = struct {
    lower_mol_per_m3: f64,
    upper_mol_per_m3: f64,
};

fn solveBoundedPhosphateExtents(
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

fn solveProjectedReactionSpanLeastSquares(
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

fn solveBoundedReactionSpan(
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

test "bounded reaction span rank reveals duplicate conservative directions" {
    var workspace = try Workspace.init(std.testing.allocator);
    defer workspace.deinit();
    const matrix = workspace.reaction_span_jacobian[0 .. 3 * 2];
    matrix[0] = 1;
    matrix[1] = 1;
    matrix[2] = 2;
    matrix[3] = 2;
    matrix[4] = 3;
    matrix[5] = 3;
    workspace.reaction_span_rhs[0] = 2;
    workspace.reaction_span_rhs[1] = 4;
    workspace.reaction_span_rhs[2] = 6;
    workspace.reaction_span_lower_bounds[0] = -1;
    workspace.reaction_span_lower_bounds[1] = -1;
    workspace.reaction_span_upper_bounds[0] = 1;
    workspace.reaction_span_upper_bounds[1] = 1;

    try std.testing.expect(solveBoundedReactionSpan(
        &workspace,
        3,
        2,
    ));
    try std.testing.expectApproxEqAbs(
        @as(f64, 2),
        workspace.reaction_span_solution[0] +
            workspace.reaction_span_solution[1],
        1.0e-12,
    );
    for (workspace.reaction_span_solution[0..2]) |extent| {
        try std.testing.expect(extent >= -1 and extent <= 1);
    }
}

fn phosphateExtentBounds(
    current: []const f64,
    reaction: CoupledExtentReaction,
    parameters: chemistry.ReactionParameters,
    options: Options,
) PhosphateExtentBounds {
    return .{
        .lower_mol_per_m3 = -options.maximum_newton_fraction *
            maximumPhosphateExtent(current, reaction, -1, parameters),
        .upper_mol_per_m3 = options.maximum_newton_fraction *
            maximumPhosphateExtent(current, reaction, 1, parameters),
    };
}

fn evaluateGlobalResidualAt(
    scratch: *chemistry.State,
    vector: []const f64,
    parameters: chemistry.ReactionParameters,
    output: []f64,
) !void {
    const transformations = try evaluateAt(scratch, vector, parameters);
    _ = try transformedVectorAdmissible(
        scratch,
        vector,
        transformations,
        parameters,
        1,
        output,
    );
    for (output, vector) |*change, value| change.* -= value;
}

fn evaluatePhosphateExtentResiduals(
    scratch: *chemistry.State,
    vector: []const f64,
    parameters: chemistry.ReactionParameters,
    output: []f64,
) !void {
    if (output.len != coupled_extent_reaction_count)
        return error.PhosphateExtentVectorSizeMismatch;
    try scratch.unpackCell(0, vector);
    const initial_coefficients =
        try scratch.activityCoefficients(0, parameters.fractions);
    const water = try water_equilibrium.solve(.{
        .hydrogen_concentration_mol_per_m3 = scratch.aqueous[0].hydrogen,
        .hydroxide_concentration_mol_per_m3 = scratch.aqueous[0].hydroxide,
        .monovalent_activity_coefficient = initial_coefficients.monovalent_activity_coefficient,
        .water_activity_product_mol2_per_m6 = parameters.water_activity_product_mol2_per_m6,
        .negligible_concentration_mol_per_m3 = parameters.negligible_water_ion_concentration_mol_per_m3,
    });
    scratch.aqueous[0].hydrogen =
        water.hydrogen_concentration_mol_per_m3;
    scratch.aqueous[0].hydroxide =
        water.hydroxide_concentration_mol_per_m3;
    const coefficients =
        try scratch.activityCoefficients(0, parameters.fractions);
    for (0..2) |zone_index| {
        const fraction = phosphateZoneFraction(parameters, zone_index);
        const offset = zone_index * 8;
        if (fraction == 0) {
            @memset(output[offset..][0..8], 0);
            continue;
        }
        const density = phosphateZoneDensity(parameters, zone_index);
        const zone = if (zone_index == 0)
            scratch.non_band_phosphate[0]
        else
            scratch.band_phosphate[0];
        const fluxes = try phosphate_rates.calculate(
            scratch.aqueous[0],
            zone,
            coefficients,
            density,
            parameters.phosphate_constants,
            parameters.phosphate_surface,
            parameters.phosphate_minerals,
            parameters.phosphate_kinetics,
        );
        output[offset] =
            fluxes.aqueous.hpo4_hydrogen_association_mol_p_per_m3;
        output[offset + 1] =
            fluxes.aqueous.iron_hpo4_pairing_mol_p_per_m3;
        output[offset + 2] =
            fluxes.aqueous.iron_h2po4_pairing_mol_p_per_m3;
        output[offset + 3] =
            fluxes.aqueous.calcium_hpo4_pairing_mol_p_per_m3;
        output[offset + 4] =
            fluxes.aqueous.calcium_h2po4_pairing_mol_p_per_m3;
        output[offset + 5] =
            fluxes.surface.h2po4_with_protonated_site_mol_p_per_megagram *
            density;
        output[offset + 6] =
            fluxes.surface.h2po4_with_hydroxyl_site_mol_p_per_megagram *
            density;
        output[offset + 7] =
            fluxes.surface.hpo4_with_hydroxyl_site_mol_p_per_megagram *
            density;
    }
    const aqueous_fluxes = try aqueous_rates.calculate(
        scratch.aqueous[0],
        coefficients,
        parameters.aqueous_constants,
        parameters.aqueous_kinetics,
    );
    output[
        @intFromEnum(
            CoupledExtentReaction.aqueous_calcium_hydroxide_pairing,
        )
    ] = aqueous_fluxes.calcium_hydroxide_association;
    output[
        @intFromEnum(
            CoupledExtentReaction.aqueous_calcium_carbonate_pairing,
        )
    ] = aqueous_fluxes.calcium_carbonate_association;
    output[
        @intFromEnum(
            CoupledExtentReaction.aqueous_calcium_bicarbonate_pairing,
        )
    ] = aqueous_fluxes.calcium_bicarbonate_association;
    output[
        @intFromEnum(
            CoupledExtentReaction.aqueous_calcium_sulfate_pairing,
        )
    ] = aqueous_fluxes.calcium_sulfate_association;
}

fn applyPhosphateExtent(
    vector: []f64,
    reaction: CoupledExtentReaction,
    extent_mol_per_m3: f64,
    parameters: chemistry.ReactionParameters,
) bool {
    if (!std.math.isFinite(extent_mol_per_m3)) return false;
    const ordinal: usize = @intFromEnum(reaction);
    if (ordinal >= 16) {
        switch (reaction) {
            .aqueous_calcium_hydroxide_pairing => applyAqueousCalciumExtent(
                vector,
                "hydroxide",
                "calcium_hydroxide",
                extent_mol_per_m3,
            ),
            .aqueous_calcium_carbonate_pairing => applyAqueousCalciumExtent(
                vector,
                "carbonate",
                "calcium_carbonate",
                extent_mol_per_m3,
            ),
            .aqueous_calcium_bicarbonate_pairing => applyAqueousCalciumExtent(
                vector,
                "bicarbonate",
                "calcium_bicarbonate",
                extent_mol_per_m3,
            ),
            .aqueous_calcium_sulfate_pairing => applyAqueousCalciumExtent(
                vector,
                "sulfate",
                "calcium_sulfate",
                extent_mol_per_m3,
            ),
            else => unreachable,
        }
        for (vector) |value| if (!std.math.isFinite(value)) return false;
        return true;
    }
    const zone_index = ordinal / 8;
    const family = ordinal % 8;
    const fraction = phosphateZoneFraction(parameters, zone_index);
    const density = phosphateZoneDensity(parameters, zone_index);
    if (!std.math.isFinite(fraction) or fraction <= 0 or
        !std.math.isFinite(density) or density <= 0)
    {
        return false;
    }
    switch (family) {
        0 => {
            addPacked(
                vector,
                phosphatePackedIndex(zone_index, "dissolved_hpo4_mol_p_per_m3"),
                -extent_mol_per_m3,
            );
            addPacked(
                vector,
                phosphatePackedIndex(zone_index, "dissolved_h2po4_mol_p_per_m3"),
                extent_mol_per_m3,
            );
            addPacked(
                vector,
                aqueousPackedIndex("hydrogen"),
                -fraction * extent_mol_per_m3,
            );
        },
        1 => applyPhosphatePairingExtent(
            vector,
            zone_index,
            fraction,
            "iron",
            "dissolved_hpo4_mol_p_per_m3",
            "iron_hpo4_pair_mol_per_m3",
            extent_mol_per_m3,
        ),
        2 => applyPhosphatePairingExtent(
            vector,
            zone_index,
            fraction,
            "iron",
            "dissolved_h2po4_mol_p_per_m3",
            "iron_h2po4_pair_mol_per_m3",
            extent_mol_per_m3,
        ),
        3 => applyPhosphatePairingExtent(
            vector,
            zone_index,
            fraction,
            "calcium",
            "dissolved_hpo4_mol_p_per_m3",
            "calcium_hpo4_pair_mol_per_m3",
            extent_mol_per_m3,
        ),
        4 => applyPhosphatePairingExtent(
            vector,
            zone_index,
            fraction,
            "calcium",
            "dissolved_h2po4_mol_p_per_m3",
            "calcium_h2po4_pair_mol_per_m3",
            extent_mol_per_m3,
        ),
        5 => {
            applySiteExchangeExtent(
                vector,
                zone_index,
                density,
                "dissolved_h2po4_mol_p_per_m3",
                "protonated_site_mol_per_megagram",
                "adsorbed_h2po4_mol_p_per_megagram",
                extent_mol_per_m3,
            );
            addPacked(
                vector,
                chemistry.State.packedComponentCount() - 1,
                fraction * extent_mol_per_m3,
            );
        },
        6 => {
            applySiteExchangeExtent(
                vector,
                zone_index,
                density,
                "dissolved_h2po4_mol_p_per_m3",
                "hydroxyl_site_mol_per_megagram",
                "adsorbed_h2po4_mol_p_per_megagram",
                extent_mol_per_m3,
            );
            addPacked(
                vector,
                aqueousPackedIndex("hydroxide"),
                fraction * extent_mol_per_m3,
            );
        },
        7 => {
            applySiteExchangeExtent(
                vector,
                zone_index,
                density,
                "dissolved_hpo4_mol_p_per_m3",
                "hydroxyl_site_mol_per_megagram",
                "adsorbed_hpo4_mol_p_per_megagram",
                extent_mol_per_m3,
            );
            addPacked(
                vector,
                aqueousPackedIndex("hydroxide"),
                fraction * extent_mol_per_m3,
            );
        },
        else => unreachable,
    }
    for (vector) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn applyAqueousCalciumExtent(
    vector: []f64,
    comptime ligand_name: []const u8,
    comptime pair_name: []const u8,
    extent_mol_per_m3: f64,
) void {
    applyAqueousAssociationExtent(
        vector,
        "calcium",
        ligand_name,
        pair_name,
        extent_mol_per_m3,
    );
}

fn applyAqueousAssociationExtent(
    vector: []f64,
    comptime first_name: []const u8,
    comptime second_name: []const u8,
    comptime product_name: []const u8,
    extent_mol_per_m3: f64,
) void {
    addPacked(vector, aqueousPackedIndex(first_name), -extent_mol_per_m3);
    addPacked(vector, aqueousPackedIndex(second_name), -extent_mol_per_m3);
    addPacked(vector, aqueousPackedIndex(product_name), extent_mol_per_m3);
}

fn applyPhosphatePairingExtent(
    vector: []f64,
    zone_index: usize,
    fraction: f64,
    comptime aqueous_metal_name: []const u8,
    comptime dissolved_phosphate_name: []const u8,
    comptime pair_name: []const u8,
    extent_mol_per_m3: f64,
) void {
    addPacked(
        vector,
        phosphatePackedIndex(zone_index, dissolved_phosphate_name),
        -extent_mol_per_m3,
    );
    addPacked(
        vector,
        phosphatePackedIndex(zone_index, pair_name),
        extent_mol_per_m3,
    );
    addPacked(
        vector,
        aqueousPackedIndex(aqueous_metal_name),
        -fraction * extent_mol_per_m3,
    );
}

fn applySiteExchangeExtent(
    vector: []f64,
    zone_index: usize,
    density_megagrams_per_m3: f64,
    comptime dissolved_name: []const u8,
    comptime site_name: []const u8,
    comptime adsorbed_name: []const u8,
    extent_mol_per_m3: f64,
) void {
    addPacked(
        vector,
        phosphatePackedIndex(zone_index, dissolved_name),
        -extent_mol_per_m3,
    );
    addPacked(
        vector,
        phosphatePackedIndex(zone_index, site_name),
        -extent_mol_per_m3 / density_megagrams_per_m3,
    );
    addPacked(
        vector,
        phosphatePackedIndex(zone_index, adsorbed_name),
        extent_mol_per_m3 / density_megagrams_per_m3,
    );
}

fn addPacked(vector: []f64, index: usize, change: f64) void {
    vector[index] += change;
}

fn maximumPhosphateExtent(
    vector: []const f64,
    reaction: CoupledExtentReaction,
    direction: f64,
    parameters: chemistry.ReactionParameters,
) f64 {
    const ordinal: usize = @intFromEnum(reaction);
    if (ordinal >= 16) return switch (reaction) {
        .aqueous_calcium_hydroxide_pairing => maximumAqueousCalciumExtent(
            vector,
            "hydroxide",
            "calcium_hydroxide",
            direction,
        ),
        .aqueous_calcium_carbonate_pairing => maximumAqueousCalciumExtent(
            vector,
            "carbonate",
            "calcium_carbonate",
            direction,
        ),
        .aqueous_calcium_bicarbonate_pairing => maximumAqueousCalciumExtent(
            vector,
            "bicarbonate",
            "calcium_bicarbonate",
            direction,
        ),
        .aqueous_calcium_sulfate_pairing => maximumAqueousCalciumExtent(
            vector,
            "sulfate",
            "calcium_sulfate",
            direction,
        ),
        else => unreachable,
    };
    const zone_index = ordinal / 8;
    const family = ordinal % 8;
    const fraction = phosphateZoneFraction(parameters, zone_index);
    const density = phosphateZoneDensity(parameters, zone_index);
    if (fraction <= 0 or density <= 0) return 0;
    const hpo4 = vector[
        phosphatePackedIndex(zone_index, "dissolved_hpo4_mol_p_per_m3")
    ];
    const h2po4 = vector[
        phosphatePackedIndex(zone_index, "dissolved_h2po4_mol_p_per_m3")
    ];
    if (direction > 0) return switch (family) {
        0 => @min(
            hpo4,
            vector[aqueousPackedIndex("hydrogen")] / fraction,
        ),
        1 => @min(hpo4, vector[aqueousPackedIndex("iron")] / fraction),
        2 => @min(h2po4, vector[aqueousPackedIndex("iron")] / fraction),
        3 => @min(
            hpo4,
            vector[aqueousPackedIndex("calcium")] / fraction,
        ),
        4 => @min(h2po4, vector[aqueousPackedIndex("calcium")] / fraction),
        5 => @min(
            h2po4,
            density * vector[
                phosphatePackedIndex(zone_index, "protonated_site_mol_per_megagram")
            ],
        ),
        6 => @min(
            h2po4,
            density * vector[
                phosphatePackedIndex(zone_index, "hydroxyl_site_mol_per_megagram")
            ],
        ),
        7 => @min(
            hpo4,
            density * vector[
                phosphatePackedIndex(zone_index, "hydroxyl_site_mol_per_megagram")
            ],
        ),
        else => unreachable,
    };
    return switch (family) {
        0 => h2po4,
        1 => vector[
            phosphatePackedIndex(zone_index, "iron_hpo4_pair_mol_per_m3")
        ],
        2 => vector[
            phosphatePackedIndex(zone_index, "iron_h2po4_pair_mol_per_m3")
        ],
        3 => vector[
            phosphatePackedIndex(zone_index, "calcium_hpo4_pair_mol_per_m3")
        ],
        4 => vector[
            phosphatePackedIndex(zone_index, "calcium_h2po4_pair_mol_per_m3")
        ],
        5 => @min(
            density * vector[
                phosphatePackedIndex(zone_index, "adsorbed_h2po4_mol_p_per_megagram")
            ],
            vector[chemistry.State.packedComponentCount() - 1] / fraction,
        ),
        6 => @min(
            density * vector[
                phosphatePackedIndex(zone_index, "adsorbed_h2po4_mol_p_per_megagram")
            ],
            vector[aqueousPackedIndex("hydroxide")] / fraction,
        ),
        7 => @min(
            density * vector[
                phosphatePackedIndex(zone_index, "adsorbed_hpo4_mol_p_per_megagram")
            ],
            vector[aqueousPackedIndex("hydroxide")] / fraction,
        ),
        else => unreachable,
    };
}

fn maximumAqueousCalciumExtent(
    vector: []const f64,
    comptime ligand_name: []const u8,
    comptime pair_name: []const u8,
    direction: f64,
) f64 {
    return maximumAqueousAssociationExtent(
        vector,
        "calcium",
        ligand_name,
        pair_name,
        direction,
    );
}

fn maximumAqueousAssociationExtent(
    vector: []const f64,
    comptime first_name: []const u8,
    comptime second_name: []const u8,
    comptime product_name: []const u8,
    direction: f64,
) f64 {
    if (direction > 0) return @min(
        vector[aqueousPackedIndex(first_name)],
        vector[aqueousPackedIndex(second_name)],
    );
    return vector[aqueousPackedIndex(product_name)];
}

fn phosphateExtentCharacteristic(
    vector: []const f64,
    reaction: CoupledExtentReaction,
    parameters: chemistry.ReactionParameters,
) f64 {
    return @max(
        maximumPhosphateExtent(vector, reaction, 1, parameters),
        maximumPhosphateExtent(vector, reaction, -1, parameters),
    );
}

fn phosphateZoneFraction(
    parameters: chemistry.ReactionParameters,
    zone_index: usize,
) f64 {
    return if (zone_index == 0)
        parameters.fractions.phosphate_non_band
    else
        parameters.fractions.phosphate_band;
}

fn phosphateZoneDensity(
    parameters: chemistry.ReactionParameters,
    zone_index: usize,
) f64 {
    return if (zone_index == 0)
        parameters.non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3
    else
        parameters.band_phosphate_soil_mass_per_water_volume_megagrams_per_m3;
}

fn phosphateExtentControlsPackedIndex(index: usize) bool {
    if (index == aqueousPackedIndex("hydrogen") or
        index == aqueousPackedIndex("hydroxide") or
        index == aqueousPackedIndex("iron") or
        index == aqueousPackedIndex("calcium") or
        index == aqueousPackedIndex("magnesium") or
        index == aqueousPackedIndex("carbonate") or
        index == aqueousPackedIndex("bicarbonate") or
        index == aqueousPackedIndex("sulfate") or
        index == aqueousPackedIndex("calcium_hydroxide") or
        index == aqueousPackedIndex("calcium_carbonate") or
        index == aqueousPackedIndex("calcium_bicarbonate") or
        index == aqueousPackedIndex("calcium_sulfate") or
        index == chemistry.State.packedComponentCount() - 1)
    {
        return true;
    }
    return phosphateZoneExtentControlsPackedIndex(index);
}

fn coupledExtentReactionEnabled(
    reaction: CoupledExtentReaction,
    _: usize,
) bool {
    if (@intFromEnum(reaction) < 16) return true;
    return switch (reaction) {
        .aqueous_calcium_hydroxide_pairing,
        .aqueous_calcium_carbonate_pairing,
        .aqueous_calcium_bicarbonate_pairing,
        .aqueous_calcium_sulfate_pairing,
        => true,
        else => unreachable,
    };
}

fn phosphateZoneExtentControlsPackedIndex(index: usize) bool {
    for (0..2) |zone_index| {
        inline for (.{
            "dissolved_po4_mol_p_per_m3",
            "dissolved_hpo4_mol_p_per_m3",
            "dissolved_h2po4_mol_p_per_m3",
            "dissolved_h3po4_mol_p_per_m3",
            "iron_hpo4_pair_mol_per_m3",
            "iron_h2po4_pair_mol_per_m3",
            "calcium_hpo4_pair_mol_per_m3",
            "calcium_h2po4_pair_mol_per_m3",
            "protonated_site_mol_per_megagram",
            "hydroxyl_site_mol_per_megagram",
            "adsorbed_h2po4_mol_p_per_megagram",
            "adsorbed_hpo4_mol_p_per_megagram",
        }) |field_name| {
            if (index == phosphatePackedIndex(zone_index, field_name))
                return true;
        }
    }
    return false;
}

fn largestPhosphateZoneScaledResidual(
    current: []const f64,
    residual: []const f64,
    options: Options,
) f64 {
    var largest: f64 = 0;
    for (current, residual, 0..) |value, change, index| {
        if (!phosphateZoneExtentControlsPackedIndex(index)) continue;
        largest = @max(
            largest,
            @abs(change) / residualScale(value, options),
        );
    }
    return largest;
}

fn aqueousPackedIndex(comptime field_name: []const u8) usize {
    return std.meta.fieldIndex(aqueous_network.State, field_name).?;
}

fn phosphatePackedIndex(
    zone_index: usize,
    comptime field_name: []const u8,
) usize {
    const aqueous_count =
        @typeInfo(aqueous_network.State).@"struct".fields.len;
    const phosphate_count =
        @typeInfo(phosphate_network.State).@"struct".fields.len;
    return aqueous_count + zone_index * phosphate_count +
        std.meta.fieldIndex(phosphate_network.State, field_name).?;
}

fn largestScaledResidualIndex(
    current: []const f64,
    residual: []const f64,
    options: Options,
) usize {
    var limiting_index: usize = 0;
    var limiting_scaled_residual: f64 = 0;
    for (current, residual, 0..) |value, change, index| {
        const scaled_residual = @abs(change) / residualScale(value, options);
        if (scaled_residual > limiting_scaled_residual) {
            limiting_scaled_residual = scaled_residual;
            limiting_index = index;
        }
    }
    return limiting_index;
}

/// Names the component that stagnated, on the stagnation path.
///
/// The convergence-failure path has always reported
/// `SOLUTE largest residual: ... name=...`, but the stagnation path emitted only
/// the hydrogen decomposition, so a stagnating cell could not be attributed to a
/// component at all. Lane A9 reported this as an observability gap after it had
/// to describe four stagnating examples without being able to say what stagnated.
/// `limiting_index` is already computed each iteration for the trace, so this
/// only reports a value that was being discarded.
fn logTerminalStagnationComponent(
    current: []const f64,
    residual: []const f64,
    options: Options,
    limiting_index: usize,
) void {
    if (limiting_index >= current.len or limiting_index >= residual.len) return;
    const value = current[limiting_index];
    const change = residual[limiting_index];
    std.log.debug(
        "SOLUTE stagnated component: packed_component={d} name={s} value={e} change={e} scaled={e}",
        .{
            limiting_index,
            chemistry.State.packedComponentName(limiting_index) orelse "unknown",
            value,
            change,
            @abs(change) / residualScale(value, options),
        },
    );
}

fn selectTraceCandidate(
    entry: ?*IterationDiagnostic,
    candidate: CandidateKind,
    maximum_scaled_residual: f64,
) void {
    if (entry) |diagnostic| {
        diagnostic.selected_candidate = candidate;
        diagnostic.selected_maximum_scaled_residual =
            maximum_scaled_residual;
    }
}

fn tryProjectWaterPair(
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

fn solvePivotedHouseholder(
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

fn matrixColumnNorm(
    matrix: []const f64,
    row_count: usize,
    column_count: usize,
    column: usize,
    first_row: usize,
) f64 {
    var scale: f64 = 0;
    var sum_squares: f64 = 1;
    for (first_row..row_count) |row| {
        const magnitude = @abs(matrix[row * column_count + column]);
        if (magnitude == 0) continue;
        if (scale < magnitude) {
            const ratio = scale / magnitude;
            sum_squares = 1 + sum_squares * ratio * ratio;
            scale = magnitude;
        } else {
            const ratio = magnitude / scale;
            sum_squares += ratio * ratio;
        }
    }
    return if (scale == 0) 0 else scale * @sqrt(sum_squares);
}

fn hasKineticGeochemistry(parameters: chemistry.ReactionParameters) bool {
    const kinetics = parameters.geochemistry_kinetics;
    return kinetics.maximum_hydroxide_mineral_mol_per_m3_step > 0 or
        kinetics.maximum_general_mineral_mol_per_m3_step > 0 or
        kinetics.maximum_natural_weathering_mol_per_m3_step > 0 or
        kinetics.maximum_ground_weathering_mol_per_m3_step > 0;
}

/// SOLUTE's association, ion-pairing, exchange, and mineral-saturation
/// branches calculate thermodynamic targets. Their legacy T*H limits are
/// fixed-cycle relaxation controls, not hourly source rates. Newton/Picard
/// closure therefore retains inventory-fraction bounds but removes these
/// dimensional caps. True kinetic sources are operator-split separately.
fn equilibriumClosureParameters(
    parameters: chemistry.ReactionParameters,
) chemistry.ReactionParameters {
    var result = parameters;
    // Reaction implementations retain their substrate/product fractional
    // bounds. A very large finite ceiling removes only the flat dimensional
    // limiter while satisfying finite-input validation.
    const unlimited = 1.0e100;
    if (result.carboxyl_exchange_parameters
        .maximum_exchange_mol_per_m3_per_iteration > 0)
    {
        result.carboxyl_exchange_parameters
            .maximum_exchange_mol_per_m3_per_iteration = unlimited;
    }
    if (result.aqueous_kinetics.maximum_fast_association_mol_per_m3_step > 0)
        result.aqueous_kinetics.maximum_fast_association_mol_per_m3_step =
            unlimited;
    if (result.aqueous_kinetics.maximum_slow_association_mol_per_m3_step > 0)
        result.aqueous_kinetics.maximum_slow_association_mol_per_m3_step =
            unlimited;
    if (result.phosphate_surface.maximum_exchange_mol_per_megagram_step > 0)
        result.phosphate_surface.maximum_exchange_mol_per_megagram_step = unlimited;
    if (result.phosphate_minerals) |*minerals| {
        if (minerals.maximum_phosphate_precipitation_mol_per_m3_step > 0)
            minerals.maximum_phosphate_precipitation_mol_per_m3_step =
                unlimited;
        if (minerals.maximum_apatite_precipitation_mol_per_m3_step > 0)
            minerals.maximum_apatite_precipitation_mol_per_m3_step =
                unlimited;
        if (minerals.maximum_mineral_dissolution_mol_per_m3_step > 0)
            minerals.maximum_mineral_dissolution_mol_per_m3_step = unlimited;
    }
    if (result.phosphate_kinetics.maximum_pairing_mol_per_m3_step > 0)
        result.phosphate_kinetics.maximum_pairing_mol_per_m3_step = unlimited;
    if (result.cation_exchange_parameters
        .maximum_adsorption_mol_charge_per_m3_step > 0)
    {
        result.cation_exchange_parameters
            .maximum_adsorption_mol_charge_per_m3_step = unlimited;
    }
    // Geochemical mineral transformations are bounded per-step kinetic
    // extents. They are applied once between equilibrium closures below.
    result.geochemistry_kinetics
        .maximum_hydroxide_mineral_mol_per_m3_step = 0;
    result.geochemistry_kinetics
        .maximum_general_mineral_mol_per_m3_step = 0;
    result.geochemistry_kinetics
        .maximum_natural_weathering_mol_per_m3_step = 0;
    result.geochemistry_kinetics
        .maximum_ground_weathering_mol_per_m3_step = 0;
    return result;
}

fn applyKineticGeochemistryStep(
    state: *chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
) !void {
    const changes = try state.evaluateGeochemistryTransformations(
        cell_index,
        parameters.fractions,
        parameters.geochemistry_products,
        parameters.geochemistry_kinetics,
    );
    var transformations = std.mem.zeroes(chemistry.CellTransformations);
    transformations.non_band_phosphate_water_fraction =
        parameters.fractions.phosphate_non_band;
    transformations.band_phosphate_water_fraction =
        parameters.fractions.phosphate_band;
    transformations.cation_exchange_water_ratios =
        parameters.cation_exchange_water_ratios;
    transformations.geochemistry = changes;
    transformations.carboxyl_soil_mass_per_water_volume_megagrams_per_m3 =
        parameters.cation_exchange_water_ratios.shared_megagrams_per_m3;
    try state.commitCell(cell_index, transformations);
}

fn combineResults(first: Result, second: Result) Result {
    return .{
        .iterations = std.math.add(
            u16,
            first.iterations,
            second.iterations,
        ) catch std.math.maxInt(u16),
        .newton_raphson_steps = std.math.add(
            u16,
            first.newton_raphson_steps,
            second.newton_raphson_steps,
        ) catch std.math.maxInt(u16),
        .picard_steps = std.math.add(
            u16,
            first.picard_steps,
            second.picard_steps,
        ) catch std.math.maxInt(u16),
        .maximum_scaled_residual = @max(
            first.maximum_scaled_residual,
            second.maximum_scaled_residual,
        ),
        .converged = first.converged and second.converged,
    };
}

/// Depth-one Anderson acceleration of the complete conservative fixed-point
/// map. Both mapped endpoints share the same elemental and site totals, so
/// their affine secant candidate retains those totals. Backtracking preserves
/// nonnegative inventories; H+/OH- are then projected onto Kw.
fn tryAcceptAndersonCandidate(
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

fn evaluateCandidateResidualAtFraction(
    scratch: *chemistry.State,
    current: []const f64,
    target: []const f64,
    accepted_state: []f64,
    residual_work: []f64,
    parameters: chemistry.ReactionParameters,
    options: Options,
    fraction: f64,
) !f64 {
    for (accepted_state, current, target) |*value, from, to|
        value.* = from + fraction * (to - from);
    try scratch.unpackCell(0, accepted_state);
    const coefficients =
        try scratch.activityCoefficients(0, parameters.fractions);
    const water = try water_equilibrium.projectProvisional(
        scratch.aqueous[0].hydrogen,
        scratch.aqueous[0].hydroxide,
        coefficients.monovalent_activity_coefficient,
        parameters.water_activity_product_mol2_per_m6,
    );
    scratch.aqueous[0].hydrogen =
        water.hydrogen_concentration_mol_per_m3;
    scratch.aqueous[0].hydroxide =
        water.hydroxide_concentration_mol_per_m3;
    try scratch.packCell(0, accepted_state);
    const changes = try evaluateAt(scratch, accepted_state, parameters);
    _ = try transformedVectorAdmissible(
        scratch,
        accepted_state,
        changes,
        parameters,
        1,
        residual_work,
    );
    for (residual_work, accepted_state) |*change, value|
        change.* -= value;
    return scaledNorm(accepted_state, residual_work, options);
}

fn phosphateTrustRegionFraction(
    current: []const f64,
    target: []const f64,
    options: Options,
    maximum_relative_change: f64,
) f64 {
    if (maximum_relative_change <= 0) return 1;
    const aqueous_count = @typeInfo(aqueous_network.State).@"struct".fields.len;
    const phosphate_count = @typeInfo(phosphate_network.State).@"struct".fields.len;
    // Do not let the neighboring HPO4 branch prevent a conservative trial
    // from landing exactly on an H2PO4 nonnegativity boundary. The candidate
    // still must pass the same global monotonic merit test, and the next map
    // evaluation determines whether the boundary remains active.
    for (0..2) |zone| {
        const h2po4_index = aqueous_count + zone * phosphate_count + 2;
        if (current[h2po4_index] > options.absolute_tolerance and
            target[h2po4_index] <= options.absolute_tolerance)
            return 1;
    }
    var fraction: f64 = 1;
    for (0..2) |zone| {
        const index = aqueous_count + zone * phosphate_count + 1;
        const from = current[index];
        const to = target[index];
        const change = @abs(to - from);
        if (change == 0) continue;
        const scale = residualScale(from, options);
        fraction = @min(
            fraction,
            @max(maximum_relative_change, options.picard_relaxation) *
                @max(@abs(from), scale) / change,
        );
    }
    return std.math.clamp(fraction, std.math.floatEps(f64), 1);
}

test "phosphate trust region uses runtime Picard relaxation and limits only HPO4" {
    const count = chemistry.State.packedComponentCount();
    const current = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(current);
    const target = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(target);
    @memset(current, 1);
    @memset(target, 1);
    const aqueous_count = @typeInfo(aqueous_network.State).@"struct".fields.len;
    const phosphate_count = @typeInfo(phosphate_network.State).@"struct".fields.len;
    target[aqueous_count + 1] = 3;
    target[aqueous_count + phosphate_count] = 1.0e12;
    const fraction = phosphateTrustRegionFraction(
        current,
        target,
        .{},
        0.2,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), fraction, 1e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.5),
        current[aqueous_count + 1] +
            fraction * (target[aqueous_count + 1] - current[aqueous_count + 1]),
        1e-15,
    );
    const smaller_runtime_fraction = phosphateTrustRegionFraction(
        current,
        target,
        .{ .picard_relaxation = 0.3 },
        0.2,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.15),
        smaller_runtime_fraction,
        1e-15,
    );
}

/// Proves a complete relaxed fixed-point path off to the side, then publishes
/// only its converged endpoint. No intermediate lookahead state is realized.
fn tryTransactionalEquilibriumLookahead(
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

fn rememberHistory(
    current: []const f64,
    residual: []const f64,
    previous_state: []f64,
    previous_residual: []f64,
    previous_previous_state: []f64,
    previous_previous_residual: []f64,
    history_count: *u2,
) void {
    if (history_count.* >= 1) {
        @memcpy(previous_previous_state, previous_state);
        @memcpy(previous_previous_residual, previous_residual);
    }
    @memcpy(previous_state, current);
    @memcpy(previous_residual, residual);
    history_count.* = @min(2, history_count.* + 1);
}

fn logRejectedFullStep(
    scratch: *chemistry.State,
    current: []const f64,
    transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
    diagnostic_output: []f64,
) void {
    scratch.unpackCell(0, current) catch return;
    const coefficients =
        scratch.activityCoefficients(0, parameters.fractions) catch return;
    transformedVector(
        scratch,
        current,
        transformations,
        coefficients.monovalent_activity_coefficient,
        parameters.water_activity_product_mol2_per_m6,
        1,
        diagnostic_output,
    ) catch |err| {
        std.log.warn("SOLUTE full reaction step rejected by {s}", .{@errorName(err)});
        if (err == error.NegativeAqueousState) {
            const changes =
                chemistry.State.assembledAqueousChanges(
                    transformations,
                    scratch.cation_exchange_mol_per_megagram[0],
                ) catch
                    return;
            const aqueous = scratch.aqueous[0];
            inline for (@typeInfo(aqueous_network.State).@"struct".fields, 0..) |field, index| {
                const value = @field(aqueous, field.name);
                const change = @field(changes, field.name);
                if (value + change < 0)
                    std.log.warn(
                        "SOLUTE rejected aqueous component: index={d} name={s} value={e} change={e} maximum_fraction={e}",
                        .{
                            index,
                            field.name,
                            value,
                            change,
                            value / -change,
                        },
                    );
            }
        }
    };
}

fn logLimitingInventory(
    state_values: []const f64,
    accepted_changes: []const f64,
    accepted_fraction: f64,
) void {
    if (!std.math.isFinite(accepted_fraction) or accepted_fraction <= 0)
        return;
    var limiting_index: ?usize = null;
    var limiting_fraction = std.math.inf(f64);
    for (state_values, accepted_changes, 0..) |value, accepted_change, index| {
        const raw_change = accepted_change / accepted_fraction;
        if (raw_change >= 0) continue;
        const fraction = value / -raw_change;
        if (fraction < limiting_fraction) {
            limiting_fraction = fraction;
            limiting_index = index;
        }
    }
    if (limiting_index) |index|
        std.log.warn(
            "SOLUTE active-set limit: packed_component={d} value={e} raw_change={e} maximum_fraction={e} accepted_fraction={e}",
            .{
                index,
                state_values[index],
                accepted_changes[index] / accepted_fraction,
                limiting_fraction,
                accepted_fraction,
            },
        );
}

fn logLargestResidual(
    state_values: []const f64,
    residual: []const f64,
    options: Options,
) void {
    var largest_index: usize = 0;
    var largest_scaled: f64 = 0;
    for (state_values, residual, 0..) |value, change, index| {
        const scale = options.absolute_tolerance +
            options.relative_tolerance * @max(1.0, @abs(value));
        const scaled_value = @abs(change) / scale;
        if (scaled_value > largest_scaled) {
            largest_index = index;
            largest_scaled = scaled_value;
        }
    }
    std.log.debug(
        "SOLUTE largest residual: packed_component={d} name={s} value={e} change={e} scaled={e}",
        .{
            largest_index,
            chemistry.State.packedComponentName(largest_index) orelse
                "unknown",
            state_values[largest_index],
            residual[largest_index],
            largest_scaled,
        },
    );
}

fn evaluateAt(scratch: *chemistry.State, vector: []const f64, parameters: chemistry.ReactionParameters) !chemistry.CellTransformations {
    try scratch.unpackCell(0, vector);
    const coefficients = try scratch.activityCoefficients(0, parameters.fractions);
    const water = try water_equilibrium.solve(.{
        .hydrogen_concentration_mol_per_m3 = scratch.aqueous[0].hydrogen,
        .hydroxide_concentration_mol_per_m3 = scratch.aqueous[0].hydroxide,
        .monovalent_activity_coefficient = coefficients.monovalent_activity_coefficient,
        .water_activity_product_mol2_per_m6 = parameters.water_activity_product_mol2_per_m6,
        .negligible_concentration_mol_per_m3 = parameters.negligible_water_ion_concentration_mol_per_m3,
    });
    scratch.aqueous[0].hydrogen = water.hydrogen_concentration_mol_per_m3;
    scratch.aqueous[0].hydroxide = water.hydroxide_concentration_mol_per_m3;
    return scratch.evaluateCell(0, parameters);
}

fn transformedVector(scratch: *chemistry.State, current: []const f64, transformations: chemistry.CellTransformations, monovalent_activity_coefficient: f64, water_activity_product_mol2_per_m6: f64, fraction: f64, output: []f64) !void {
    try scratch.unpackCell(0, current);
    var accepted = scaled(transformations, fraction);
    const assembled_aqueous_changes =
        try chemistry.State.assembledAqueousChanges(
            accepted,
            scratch.cation_exchange_mol_per_megagram[0],
        );
    const provisional_hydrogen =
        scratch.aqueous[0].hydrogen + assembled_aqueous_changes.hydrogen;
    const provisional_hydroxide =
        scratch.aqueous[0].hydroxide + assembled_aqueous_changes.hydroxide;
    // commitCell assembles shared proton changes from all reaction families.
    // Cancel that complete dependent-coordinate change during the conservative
    // inventory transaction; the exact Kw projection below commits it once.
    accepted.aqueous.hydrogen -= assembled_aqueous_changes.hydrogen;
    accepted.aqueous.hydroxide -= assembled_aqueous_changes.hydroxide;
    try scratch.commitCell(0, accepted);
    const water = try water_equilibrium.projectProvisional(
        provisional_hydrogen,
        provisional_hydroxide,
        monovalent_activity_coefficient,
        water_activity_product_mol2_per_m6,
    );
    scratch.aqueous[0].hydrogen = water.hydrogen_concentration_mol_per_m3;
    scratch.aqueous[0].hydroxide = water.hydroxide_concentration_mol_per_m3;
    try scratch.packCell(0, output);
    // State commit permits only sub-picomolar subtraction noise. Normalize
    // that authorized numerical zero before the next stricter charge
    // classifier sees it. H+ and OH- receive the least positive f64 because
    // thermodynamic consumers require positive activities; evaluateAt
    // immediately restores their exact Kw pair.
    for (output, 0..) |*value, index| {
        if (value.* < 0) {
            if (value.* < -1e-12)
                return error.InvalidSoluteReactionRoundoff;
            value.* = if (index == 4 or index == 5)
                std.math.floatMin(f64)
            else
                0;
        } else if ((index == 4 or index == 5) and value.* == 0) {
            value.* = std.math.floatMin(f64);
        }
    }
}

/// Finds the largest fraction no greater than the requested fraction that
/// preserves every nonnegative chemistry inventory. Individual source
/// reactions are substrate limited, but several simultaneous reactions can
/// draw from the same ion; this is the conservative line search for their
/// assembled direction.
fn transformedVectorAdmissible(scratch: *chemistry.State, current: []const f64, transformations: chemistry.CellTransformations, parameters: chemistry.ReactionParameters, requested_fraction: f64, output: []f64) !f64 {
    try scratch.unpackCell(0, current);
    const coefficients = try scratch.activityCoefficients(0, parameters.fractions);
    var fraction = requested_fraction;
    var rejected_fraction: ?f64 = null;
    var last_rejection: ?anyerror = null;
    var attempt: u8 = 0;
    while (attempt < 48) : (attempt += 1) {
        transformedVector(scratch, current, transformations, coefficients.monovalent_activity_coefficient, parameters.water_activity_product_mol2_per_m6, fraction, output) catch |err| {
            last_rejection = err;
            rejected_fraction = fraction;
            fraction *= 0.5;
            continue;
        };
        break;
    }
    // If even a 2^-48 step is inadmissible, expose the physical/domain error
    // that rejected it. Returning only a generic line-search error would hide
    // the corrupt component and violate ecosys-ng's fail-fast diagnostic
    // contract.
    if (attempt == 48) return last_rejection orelse
        error.NoAdmissibleSoluteReactionStep;
    var upper = rejected_fraction orelse return fraction;
    var lower = fraction;
    var search: u8 = 0;
    while (search < 40) : (search += 1) {
        const middle = lower + 0.5 * (upper - lower);
        transformedVector(
            scratch,
            current,
            transformations,
            coefficients.monovalent_activity_coefficient,
            parameters.water_activity_product_mol2_per_m6,
            middle,
            output,
        ) catch {
            upper = middle;
            continue;
        };
        lower = middle;
        if (upper - lower <= 16 * std.math.floatEps(f64) *
            @max(1.0, @abs(lower)))
            break;
    }
    // Land on the active boundary. The source AMAX1/AMIN1 bounds permit exact
    // exhaustion, and the next evaluation must see zero so that the active
    // reaction set can switch. H+ and OH- are projected onto Kw separately
    // and therefore cannot be exhausted by this conservative boundary.
    try transformedVector(
        scratch,
        current,
        transformations,
        coefficients.monovalent_activity_coefficient,
        parameters.water_activity_product_mol2_per_m6,
        lower,
        output,
    );
    return lower;
}

fn scaled(transformations: chemistry.CellTransformations, fraction: f64) chemistry.CellTransformations {
    var result = transformations;
    scaleStruct(aqueous_network.Transformations, &result.aqueous, fraction);
    scaleStruct(phosphate_network.Transformations, &result.non_band_phosphate, fraction);
    scaleStruct(phosphate_network.Transformations, &result.band_phosphate, fraction);
    scaleStruct(cation_exchange.Cations, &result.cation_adsorption_mol_per_megagram, fraction);
    scaleStruct(geochemistry.Transformations, &result.geochemistry, fraction);
    result.carboxyl_hydrogen_change_mol_per_megagram *= fraction;
    return result;
}

fn scaleStruct(comptime T: type, value: *T, fraction: f64) void {
    inline for (@typeInfo(T).@"struct".fields) |field| @field(value.*, field.name) *= fraction;
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, change| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(change)) return error.NonFiniteSoluteReactionState;
        maximum = @max(
            maximum,
            @abs(change) / residualScale(value, options),
        );
    }
    return maximum;
}

fn residualScale(value: f64, options: Options) f64 {
    return options.absolute_tolerance +
        options.relative_tolerance * @max(1.0, @abs(value));
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

fn validateOptions(options: Options) !void {
    if (!std.math.isFinite(options.absolute_tolerance) or options.absolute_tolerance <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or options.max_iterations == 0) return error.InvalidSoluteReactionSolverOptions;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

fn inorganicCarbonMolPerM3(state: *const chemistry.State, cell: usize) f64 {
    const aqueous = state.aqueous[cell];
    return aqueous.carbon_dioxide +
        aqueous.carbonate +
        aqueous.bicarbonate +
        aqueous.calcium_carbonate +
        aqueous.calcium_bicarbonate +
        aqueous.magnesium_carbonate +
        aqueous.magnesium_bicarbonate +
        aqueous.sodium_carbonate +
        state.geochemistry_solids[cell].calcite_solid_mol_per_m3;
}

test "dependent water ions do not limit a conservative chemistry step" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filled(aqueous_network.State, 1);
    state.non_band_phosphate[0] = filled(phosphate_network.State, 1);
    state.band_phosphate[0] = filled(phosphate_network.State, 1);
    state.cation_exchange_mol_per_megagram[0] =
        filled(cation_exchange.Cations, 1);
    state.geochemistry_solids[0] = filled(geochemistry.SolidState, 1);
    state.carboxyl_bound_hydrogen_mol_per_megagram[0] = 1;
    state.water_mol_per_m3[0] = 1;
    var scratch = try chemistry.State.init(std.testing.allocator, 1);
    defer scratch.deinit();
    var current: [chemistry.State.packedComponentCount()]f64 = undefined;
    var output: [chemistry.State.packedComponentCount()]f64 = undefined;
    try state.packCell(0, &current);
    var transformations = std.mem.zeroes(chemistry.CellTransformations);
    transformations.aqueous.hydroxide = -0.3;
    transformations.cation_exchange_water_ratios = .{
        .shared_megagrams_per_m3 = 1,
        .ammonium_non_band_megagrams_per_m3 = 1,
        .ammonium_band_megagrams_per_m3 = 1,
    };
    var test_parameters: chemistry.ReactionParameters = undefined;
    test_parameters.fractions = .{
        .ammonium_non_band = 0.8,
        .ammonium_band = 0.2,
        .nitrate_non_band = 0.6,
        .nitrate_band = 0.4,
        .phosphate_non_band = 0.7,
        .phosphate_band = 0.3,
    };
    test_parameters.water_activity_product_mol2_per_m6 = 1;
    const fraction = try transformedVectorAdmissible(
        &scratch,
        &current,
        transformations,
        test_parameters,
        4,
        &output,
    );
    try std.testing.expectEqual(@as(f64, 4), fraction);
    try std.testing.expect(output[4] > 0);
    try std.testing.expect(output[5] > 0);
    for (output) |value| try std.testing.expect(value >= 0);
}

// SOLUTE-037 regression. Extends the invariant already asserted by
// "dependent water ions do not limit a conservative chemistry step" from the
// admissibility path to the phosphate extent BOUNDS, which is where it was
// still being violated. A pH-determined OH- concentration must never cap a
// phosphate site-desorption extent, because `projectProvisional` re-derives
// OH- from the water product on every candidate, so it is not a depletable
// substrate. The source agrees: no reverse bound in `starte.f` 657--677 or
// `solute.f` 1018--1054 names `AOH1`, `COH1`, or water.
test "acidic hydroxide does not clamp phosphate site desorption" {
    var vector: [chemistry.State.packedComponentCount()]f64 = undefined;
    @memset(&vector, 0);
    var test_parameters: chemistry.ReactionParameters = undefined;
    test_parameters.fractions = .{
        .ammonium_non_band = 0.8,
        .ammonium_band = 0.2,
        .nitrate_non_band = 0.6,
        .nitrate_band = 0.4,
        .phosphate_non_band = 0.7,
        .phosphate_band = 0.3,
    };
    test_parameters
        .non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = 2;
    test_parameters
        .band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = 2;

    // `Boreal Black Spruce MB` cell 0 conditions: a strongly acidic horizon
    // where OH- is ~1e-7 while the adsorbed pool is ~1e-1, a seven-order
    // separation. Water is left at zero deliberately: under the removed
    // clamp, family 5 read the water component, so a zero there would have
    // forced a zero bound.
    vector[aqueousPackedIndex("hydroxide")] = 1.9952623149688792e-7;
    vector[phosphatePackedIndex(0, "adsorbed_h2po4_mol_p_per_megagram")] =
        1.4353423827526443e-1;
    vector[phosphatePackedIndex(0, "adsorbed_hpo4_mol_p_per_megagram")] =
        1.7756084305380726e-2;

    const density = 2.0;
    inline for (.{
        .{
            CoupledExtentReaction.non_band_h2po4_protonated_site_exchange,
            "adsorbed_h2po4_mol_p_per_megagram",
        },
        .{
            CoupledExtentReaction.non_band_h2po4_hydroxyl_site_exchange,
            "adsorbed_h2po4_mol_p_per_megagram",
        },
        .{
            CoupledExtentReaction.non_band_hpo4_hydroxyl_site_exchange,
            "adsorbed_hpo4_mol_p_per_megagram",
        },
    }) |case| {
        const reverse = maximumPhosphateExtent(
            &vector,
            case[0],
            -1,
            test_parameters,
        );
        // The bound is the adsorbed inventory, not the OH- concentration.
        try std.testing.expectApproxEqRel(
            density * vector[phosphatePackedIndex(0, case[1])],
            reverse,
            1e-15,
        );
        // And it is not squeezed to the trace OH- scale, which is what made
        // the clamp absorbing: desorption could not keep up with adsorption,
        // so the adsorbed pool ratcheted up by 30x and the solve stagnated.
        try std.testing.expect(
            reverse > 1e3 * vector[aqueousPackedIndex("hydroxide")],
        );
    }
}

test "packed residual diagnostics identify scientific components" {
    try std.testing.expectEqualStrings(
        "aqueous.ammonium_non_band",
        chemistry.State.packedComponentName(0).?,
    );
    try std.testing.expectEqualStrings(
        "water_mol_per_m3",
        chemistry.State.packedComponentName(
            chemistry.State.packedComponentCount() - 1,
        ).?,
    );
    try std.testing.expect(
        chemistry.State.packedComponentName(
            chemistry.State.packedComponentCount(),
        ) == null,
    );
}

test "day 12 snapshot rejects ammonium above water molarity" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.water_mol_per_m3[0] = 55555.906256719943;
    state.aqueous[0].ammonium_non_band = 1353591.0109655228;
    try std.testing.expectError(
        error.SoluteConcentrationExceedsWaterMolarity,
        validateAqueousMolarity(&state, 0),
    );
    state.aqueous[0].ammonium_non_band = 0.6232257391219167;
    try validateAqueousMolarity(&state, 0);
}

test "Ottawa phosphate stays accelerated and weathering applies once" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filled(aqueous_network.State, 1);
    state.non_band_phosphate[0] = filled(phosphate_network.State, 1);
    state.band_phosphate[0] = filled(phosphate_network.State, 1);
    state.cation_exchange_mol_per_megagram[0] = filled(cation_exchange.Cations, 1);
    state.geochemistry_solids[0] = filled(geochemistry.SolidState, 1);
    // Keep the synthetic chemistry within its declared physical domain.
    // Water is otherwise inert in this focused reaction-network fixture.
    state.water_mol_per_m3[0] = 100;
    const parameters: chemistry.ReactionParameters = .{
        .fractions = .{ .ammonium_non_band = 0.8, .ammonium_band = 0.2, .nitrate_non_band = 0.6, .nitrate_band = 0.4, .phosphate_non_band = 0.7, .phosphate_band = 0.3 },
        .non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = 1,
        .band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = 1,
        .cation_exchange_capacity_mol_charge_per_megagram = 10,
        .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 1, .ammonium_non_band_megagrams_per_m3 = 1, .ammonium_band_megagrams_per_m3 = 1 },
        .total_carboxyl_sites_mol_per_megagram = 0,
        .carboxyl_exchange_parameters = .{ .dissociation_constant_mol_per_m3 = 0.01, .maximum_exchange_mol_per_m3_per_iteration = 0.01, .substrate_limit_fraction_per_iteration = 0.2 },
        .aqueous_constants = filled(aqueous_rates.EquilibriumConstants, 1),
        .aqueous_kinetics = .{ .ammonium_substrate_limit_fraction = 0.2, .general_substrate_limit_fraction = 0.2, .maximum_fast_association_mol_per_m3_step = 0, .maximum_slow_association_mol_per_m3_step = 0 },
        .phosphate_constants = filled(phosphate_rates.EquilibriumConstants, 1),
        .phosphate_surface = .{ .protonated_site_equilibrium_constant = 1, .hydroxyl_site_equilibrium_constant = 1, .h2po4_exchange_equilibrium_constant = 1, .hpo4_exchange_equilibrium_constant = 1, .water_activity_product_mol2_per_m6 = 1, .h2po4_dissociation_constant = 1, .maximum_exchange_mol_per_megagram_step = 0, .substrate_limit_fraction = 0.2 },
        .phosphate_minerals = null,
        .phosphate_kinetics = .{ .substrate_limit_fraction = 0.2, .maximum_pairing_mol_per_m3_step = 0 },
        .cation_exchange_parameters = .{ .selectivity = .{ .calcium_ammonium = 1, .calcium_hydrogen = 1, .calcium_aluminum_and_iron = 1, .calcium_magnesium = 1, .calcium_sodium = 1, .calcium_potassium = 1 }, .substrate_limit_fraction = 0.2, .maximum_adsorption_mol_charge_per_m3_step = 0 },
        .geochemistry_products = filled(geochemistry_rates.SolubilityProducts, 1),
        .geochemistry_kinetics = .{ .general_substrate_limit_fraction = 0.2, .hydrogen_coupled_substrate_limit_fraction = 0.2, .maximum_hydroxide_mineral_mol_per_m3_step = 0, .maximum_general_mineral_mol_per_m3_step = 0, .calcite_hydroxide_inhibition_constant_mol_per_m3 = 1, .maximum_natural_weathering_mol_per_m3_step = 0, .maximum_ground_weathering_mol_per_m3_step = 0 },
        .water_activity_product_mol2_per_m6 = 1,
        .negligible_water_ion_concentration_mol_per_m3 = 1e-32,
    };
    // Accepted equilibrium and kinetic-calcite steps must only repartition
    // inorganic carbon among aqueous carriers and solid calcite.
    var carbon_state = try chemistry.State.init(std.testing.allocator, 1);
    defer carbon_state.deinit();
    carbon_state.aqueous[0] = filled(aqueous_network.State, 1);
    carbon_state.aqueous[0].calcium = 2;
    carbon_state.aqueous[0].carbonate = 2;
    carbon_state.non_band_phosphate[0] = filled(phosphate_network.State, 1);
    carbon_state.band_phosphate[0] = filled(phosphate_network.State, 1);
    carbon_state.cation_exchange_mol_per_megagram[0] =
        filled(cation_exchange.Cations, 1);
    carbon_state.geochemistry_solids[0] = filled(geochemistry.SolidState, 1);
    carbon_state.geochemistry_solids[0].calcite_solid_mol_per_m3 = 0.1;
    carbon_state.water_mol_per_m3[0] = 100;
    var carbon_parameters = parameters;
    carbon_parameters.geochemistry_kinetics
        .maximum_hydroxide_mineral_mol_per_m3_step = 0.01;
    const carbon_before = inorganicCarbonMolPerM3(&carbon_state, 0);
    const calcite_before =
        carbon_state.geochemistry_solids[0].calcite_solid_mol_per_m3;
    const carbon_result = try solveCell(
        std.testing.allocator,
        &carbon_state,
        0,
        carbon_parameters,
        .{ .max_iterations = 60 },
    );
    try std.testing.expect(carbon_result.converged);
    try std.testing.expect(
        carbon_state.geochemistry_solids[0].calcite_solid_mol_per_m3 !=
            calcite_before,
    );
    try std.testing.expectApproxEqAbs(
        carbon_before,
        inorganicCarbonMolPerM3(&carbon_state, 0),
        2.0e-14,
    );
    // Exact standalone reproduction of the first hourly Ottawa limiting
    // coordinate. The 0.025 and 0.00125 values are legacy fixed-cycle
    // relaxation ceilings for aqueous association and cation exchange, not
    // kinetic sources to be integrated once per model hour.
    var hourly_state = try chemistry.State.init(std.testing.allocator, 1);
    defer hourly_state.deinit();
    hourly_state.aqueous[0] = filled(aqueous_network.State, 1);
    hourly_state.aqueous[0].ammonium_non_band = 0.6025043003707329;
    hourly_state.non_band_phosphate[0] = filled(phosphate_network.State, 1);
    hourly_state.band_phosphate[0] = filled(phosphate_network.State, 1);
    hourly_state.cation_exchange_mol_per_megagram[0] =
        filled(cation_exchange.Cations, 1);
    hourly_state.geochemistry_solids[0] = filled(geochemistry.SolidState, 1);
    hourly_state.water_mol_per_m3[0] = 1;
    var hourly_parameters = parameters;
    hourly_parameters.aqueous_kinetics
        .maximum_fast_association_mol_per_m3_step = 0.025;
    hourly_parameters.cation_exchange_parameters
        .maximum_adsorption_mol_charge_per_m3_step = 0.00125;
    const non_band_phosphate_before = hourly_state.non_band_phosphate[0];
    const band_phosphate_before = hourly_state.band_phosphate[0];
    const hourly_result = try solveCell(
        std.testing.allocator,
        &hourly_state,
        0,
        hourly_parameters,
        .{ .max_iterations = 60 },
    );
    try std.testing.expect(hourly_result.converged);
    try std.testing.expect(hourly_result.iterations <= 60);
    try std.testing.expect(hourly_result.newton_raphson_steps > 0);
    try std.testing.expect(hourly_result.maximum_scaled_residual <= 1);
    try std.testing.expectEqual(
        non_band_phosphate_before,
        hourly_state.non_band_phosphate[0],
    );
    try std.testing.expectEqual(
        band_phosphate_before,
        hourly_state.band_phosphate[0],
    );

    var ottawa_ceiling_parameters = parameters;
    ottawa_ceiling_parameters.aqueous_kinetics
        .maximum_slow_association_mol_per_m3_step = 0.025;
    ottawa_ceiling_parameters.phosphate_surface
        .maximum_exchange_mol_per_megagram_step = 0.0025;
    ottawa_ceiling_parameters.geochemistry_kinetics
        .maximum_hydroxide_mineral_mol_per_m3_step = 0.0025;
    const closure_parameters = equilibriumClosureParameters(
        ottawa_ceiling_parameters,
    );
    try std.testing.expectEqual(
        @as(f64, 1.0e100),
        closure_parameters.aqueous_kinetics
            .maximum_slow_association_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        closure_parameters.geochemistry_kinetics
            .maximum_hydroxide_mineral_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 1.0e100),
        closure_parameters.phosphate_surface
            .maximum_exchange_mol_per_megagram_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        closure_parameters.phosphate_kinetics
            .maximum_pairing_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        closure_parameters.aqueous_kinetics
            .maximum_fast_association_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        ottawa_ceiling_parameters.phosphate_surface
            .substrate_limit_fraction,
        closure_parameters.phosphate_surface
            .substrate_limit_fraction,
    );
    try std.testing.expectEqual(
        ottawa_ceiling_parameters.phosphate_constants.h2po4,
        closure_parameters.phosphate_constants.h2po4,
    );
    var before: [chemistry.State.packedComponentCount()]f64 = undefined;
    try state.packCell(0, &before);
    try std.testing.expectError(
        error.SoluteReactionSolverDidNotConverge,
        solveCell(
            std.testing.allocator,
            &state,
            0,
            parameters,
            .{ .max_iterations = 1 },
        ),
    );
    var after_failure: [chemistry.State.packedComponentCount()]f64 = undefined;
    try state.packCell(0, &after_failure);
    try std.testing.expectEqualSlices(f64, &before, &after_failure);

    const result = try solveCell(
        std.testing.allocator,
        &state,
        0,
        parameters,
        .{ .max_iterations = 1000 },
    );
    try std.testing.expect(result.iterations < 1000);
    try std.testing.expect(result.converged);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
    try std.testing.expect(result.maximum_scaled_residual <= 1);

    // Standalone reproduction of the Ottawa limiting coordinate. Coupled
    // H2PO4 protonation, pairing, and surface exchange must retain access to
    // Newton/multisecant acceleration after the former 12-step cutoff.
    state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 =
        0.050459466982896495;
    var phosphate_parameters = parameters;
    phosphate_parameters.phosphate_kinetics
        .maximum_pairing_mol_per_m3_step = 0.025;
    phosphate_parameters.phosphate_surface
        .maximum_exchange_mol_per_megagram_step = 0.0025;
    const phosphate_result = try solveCell(
        std.testing.allocator,
        &state,
        0,
        phosphate_parameters,
        .{ .max_iterations = 1000 },
    );
    try std.testing.expect(phosphate_result.converged);
    try std.testing.expect(phosphate_result.maximum_scaled_residual <= 1);

    // A weathering rate is a once-per-timestep kinetic source, not an
    // equilibrium residual to be driven to zero or multiplied by MRXN.
    state.aqueous[0].hydrogen = 10;
    state.aqueous[0].aluminum = 0;
    state.aqueous[0].hydrogen_silicate = 0;
    state.aqueous[0].hydroxide = 0.00041020695203390634;
    state.geochemistry_solids[0]
        .aluminum_natural_silicate_mol_per_m3 = 1;
    state.geochemistry_solids[0].gibbsite_solid_mol_per_m3 =
        1691.7202939682704;
    var kinetic_parameters = parameters;
    kinetic_parameters.geochemistry_kinetics
        .maximum_natural_weathering_mol_per_m3_step = 5.0e-5;
    kinetic_parameters.geochemistry_kinetics
        .maximum_hydroxide_mineral_mol_per_m3_step = 0.0025;
    const family_preview = try state.evaluateGeochemistryTransformations(
        0,
        kinetic_parameters.fractions,
        kinetic_parameters.geochemistry_products,
        kinetic_parameters.geochemistry_kinetics,
    );
    try std.testing.expect(
        @abs(family_preview.gibbsite_solid_mol_per_m3) <= 0.0025,
    );
    try std.testing.expect(
        @abs(family_preview.aluminum_natural_silicate_mol_per_m3) <=
            5.0e-5,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        family_preview.dissolved_aluminum_mol_per_m3 +
            family_preview.gibbsite_solid_mol_per_m3 +
            family_preview.aluminum_natural_silicate_mol_per_m3 +
            family_preview.aluminum_ground_silicate_mol_per_m3,
        1e-15,
    );
    const rock_before = state.geochemistry_solids[0]
        .aluminum_natural_silicate_mol_per_m3;
    const gibbsite_before =
        state.geochemistry_solids[0].gibbsite_solid_mol_per_m3;
    const kinetic_result = try solveCell(
        std.testing.allocator,
        &state,
        0,
        kinetic_parameters,
        .{ .max_iterations = 1000 },
    );
    const weathered = rock_before - state.geochemistry_solids[0]
        .aluminum_natural_silicate_mol_per_m3;
    const gibbsite_change = @abs(
        gibbsite_before -
            state.geochemistry_solids[0].gibbsite_solid_mol_per_m3,
    );
    try std.testing.expect(kinetic_result.converged);
    try std.testing.expect(kinetic_result.iterations <= 1000);
    try std.testing.expect(weathered > 0);
    try std.testing.expect(weathered <= 5.0e-5 + 1e-14);
    // The authoritative solid is near 1.7e3 mol m-3, so subtracting its
    // before/after values loses several decimal ulps. Test the family extent
    // above directly; this observed-state check permits only that storage
    // roundoff, not an additional kinetic extent.
    const gibbsite_storage_roundoff =
        16 * std.math.floatEps(f64) * gibbsite_before;
    try std.testing.expect(
        gibbsite_change <= 0.0025 + gibbsite_storage_roundoff,
    );
}
