const std = @import("std");

const gas = @import("transport.zig");

const atmosphere = @import("atmosphere_exchange.zig");

const numerics = @import("../../core/numerics.zig");

pub const Inputs = struct {
    faces: []const gas.Face,
    face_conductance_m3_per_step: []const f64,
    atmospheric_boundaries: []const atmosphere.Boundary,
    subsurface_boundaries: []const atmosphere.Boundary = &.{},
    water_volume_m3: []const f64,
    band_water_volume_m3: []const f64,
    mass_solubility_ratio: []const f64,
    gas_water_exchange_rate_per_step: []const f64,
    band_gas_water_exchange_rate_per_step: []const f64,
    bubbling_enabled: []const bool,
    /// REDIST `LL=MIN(L,LG)` destination for gas released by bubbling. A
    /// null entry means no gas-phase route exists, so the released mass is a
    /// landscape boundary flux. Omitting the map conservatively releases
    /// bubbles into their source cell.
    bubble_receiver_cell_by_cell: ?[]const ?usize = null,
    /// Optional caller-owned exact accepted atmospheric flux ledger, indexed
    /// cell × species. It is written only at successful commit.
    atmospheric_flux_g_by_component: ?[]f64 = null,
    /// Separate accepted lateral/lower flux ledger; positive enters soil.
    subsurface_flux_g_by_component: ?[]f64 = null,
    /// Exact accepted internal face ledger, indexed face × gas species.
    /// Positive moves `first_cell -> second_cell`.
    face_flux_g_by_component: ?[]f64 = null,
};

pub const Options = struct {
    absolute_tolerance_g: f64 = 1e-12,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
    /// Fraction of the physical hour represented by the residual. ecosys-ng
    /// solves one full-hour equation; NPG expands only the iteration ceiling.
    transport_iteration_fraction: f64 = 1,
    /// Runtime NPH*NPG ceiling derived from the input options.
    max_iterations: u16,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Solves spatial diffusion, pressure displacement, atmospheric exchange,
/// gas-water equilibration (including band NH3), and bubbling as one nonlinear
/// system. No full model process is rerun during these iterations.
pub fn solve(allocator: std.mem.Allocator, state: *gas.State, inputs: Inputs, options: Options) !Result {
    try validate(state, inputs, options);
    const inventory_count = state.gaseous_mass_g.len;
    const unknown_count = try std.math.mul(usize, inventory_count, 3);
    const base = try allocator.alloc(f64, unknown_count);
    defer allocator.free(base);
    copyStateToVector(state, base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, unknown_count);
    defer allocator.free(residual);
    const probe = try allocator.alloc(f64, unknown_count);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, unknown_count);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, unknown_count);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, unknown_count);
    defer allocator.free(candidate_residual);
    const previous = try allocator.alloc(f64, unknown_count);
    defer allocator.free(previous);
    const previous_residual = try allocator.alloc(f64, unknown_count);
    defer allocator.free(previous_residual);
    const target = try allocator.alloc(f64, unknown_count);
    defer allocator.free(target);
    var scratch = try gas.State.init(allocator, state.cell_count);
    defer scratch.deinit();
    @memcpy(scratch.air_volume_m3, state.air_volume_m3);
    @memcpy(scratch.temperature_k, state.temperature_k);
    @memcpy(scratch.water_vapor_mol, state.water_vapor_mol);

    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    var last_scaled_residual: f64 = std.math.inf(f64);
    var has_previous = false;
    const residual_fraction = options.transport_iteration_fraction;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(allocator, &scratch, base, current, inputs, residual_fraction, target, residual);
        var norm = try scaledNorm(current, residual, options);
        last_scaled_residual = norm;
        if (norm <= 1) {
            // `target` is the conservatively assembled F(current). At the
            // accepted root it differs only within the requested tolerance,
            // while committing it preserves the transport mass balance that
            // an approximate Newton iterate alone need not preserve exactly.
            if (capturesFluxLedgers(inputs)) try residualAtCapturing(allocator, &scratch, base, current, inputs, residual_fraction, target, residual, true);
            copyVectorToState(target, state);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm };
        }
        if (iteration >= options.max_iterations / 2 and
            hasExtremeIncomingDonorResidual(
                current,
                residual,
                target,
                inputs,
                options,
                inventory_count,
            ))
            break;
        // Release a zero-inventory complementarity bound with the assembled
        // conservative target before forming a numerical Jacobian there.
        // At x=0 a positive source is physical, but a projected finite-
        // difference Newton column can be singular and waste NPH*NPG.
        var limiting_index = try worstResidualIndex(current, residual, options);
        var depleted_pair_norm =
            try scaledCoordinateResidual(
                current[limiting_index],
                residual[limiting_index],
                options,
            );
        for (0..inventory_count) |gas_index| {
            const dissolved_index = inventory_count + gas_index;
            if (isScaleDepletedGas(base[gas_index], base[dissolved_index], options) and
                current[dissolved_index] > 0)
            {
                const coordinate_norm =
                    try scaledCoordinateResidual(current[gas_index], residual[gas_index], options);
                if (coordinate_norm > depleted_pair_norm) {
                    depleted_pair_norm = coordinate_norm;
                    limiting_index = gas_index;
                }
            }
        }
        var released_bound_count: usize = 0;
        @memcpy(candidate, current);
        for (current, residual, target, candidate) |value, difference, assembled, *next| {
            if (isNumericallyAtNonnegativeBound(value, assembled) and
                difference > 0)
            {
                next.* = assembled;
                released_bound_count += 1;
            }
        }
        if (released_bound_count > 1) {
            var accepted_bound_set = false;
            if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                // This is an active-set transition, not a Newton line-search
                // step. Releasing all positive-source bounds together avoids
                // consuming one NPH*NPG iteration per tied cell/species.
                if (std.math.isFinite(candidate_norm)) {
                    accepted_bound_set = true;
                    @memcpy(previous, current);
                    @memcpy(previous_residual, residual);
                    has_previous = true;
                    @memcpy(current, candidate);
                    @memcpy(residual, candidate_residual);
                    norm = candidate_norm;
                    if (norm <= 1) {
                        if (capturesFluxLedgers(inputs)) try residualAtCapturing(allocator, &scratch, base, current, inputs, residual_fraction, target, residual, true);
                        copyVectorToState(target, state);
                        return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm };
                    }
                    limiting_index = try worstResidualIndex(current, residual, options);
                }
            } else |_| {}
            // `residualAt` writes the shared assembled target. A rejected
            // active-set probe must not leave F(candidate) paired with the
            // still-current iterate.
            if (!accepted_bound_set) {
                try residualAt(allocator, &scratch, base, current, inputs, residual_fraction, target, residual);
                norm = try scaledNorm(current, residual, options);
                limiting_index = try worstResidualIndex(current, residual, options);
            }
        }
        if (isNumericallyAtNonnegativeBound(
            current[limiting_index],
            target[limiting_index],
        ) and residual[limiting_index] > 0) {
            @memcpy(candidate, current);
            candidate[limiting_index] = target[limiting_index];
            if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                const old_coordinate = try scaledCoordinateResidual(current[limiting_index], residual[limiting_index], options);
                const new_coordinate = try scaledCoordinateResidual(candidate[limiting_index], candidate_residual[limiting_index], options);
                // Releasing a complementarity bound changes the active set.
                // A donor clamp can make the first positive branch locally
                // flat, so neither coordinate nor global line-search merit is
                // valid for this one-way transition. The next iteration's
                // bracket expansion resolves that branch.
                if (candidate[limiting_index] > current[limiting_index] and
                    std.math.isFinite(new_coordinate) and std.math.isFinite(old_coordinate))
                {
                    @memcpy(previous, current);
                    @memcpy(previous_residual, residual);
                    has_previous = true;
                    @memcpy(current, candidate);
                    // Bound activation defines the active set for this
                    // iteration; it is not itself a Newton/Picard correction.
                    try residualAt(allocator, &scratch, base, current, inputs, residual_fraction, target, residual);
                    norm = try scaledNorm(current, residual, options);
                    if (norm <= 1) {
                        if (capturesFluxLedgers(inputs)) try residualAtCapturing(allocator, &scratch, base, current, inputs, residual_fraction, target, residual, true);
                        copyVectorToState(target, state);
                        return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm };
                    }
                    limiting_index = try worstResidualIndex(current, residual, options);
                }
            } else |_| {}
        }
        // Species equations are separable, but pressure displacement couples
        // their cell coordinates. Bracket the coordinate controlling the
        // infinity norm at every scale: example validation shows donor clamps
        // can require this safeguard from the first hourly solve, while each
        // successful bracket exits immediately.
        // Bracketing only one coordinate is a robust active-set safeguard,
        // but applying it on every outer iteration serializes a tied
        // cell/species block and can consume the entire NPH*NPG ceiling.
        // Give the coupled species Newton solve the middle of the budget;
        // retain bracketing for initial branch discovery and enough final
        // polishing iterations to visit every phase of all seven gases in a
        // dry-cell active set without serializing the entire ceiling.
        const scalar_polishing_iterations: u16 = @min(options.max_iterations, 3 * gas.species_count);
        const limiting_phase_for_pair = limiting_index / inventory_count;
        const limiting_gas_index = limiting_index % inventory_count;
        const limiting_dissolved_index = inventory_count + limiting_gas_index;
        const limiting_pair_is_depleted = isScaleDepletedGas(
            base[limiting_gas_index],
            base[limiting_dissolved_index],
            options,
        );
        const phase_coupled_pair = limiting_phase_for_pair < 2 and
            limiting_pair_is_depleted and
            current[limiting_dissolved_index] > 0;
        if (residual[limiting_index] != 0 and
            (phase_coupled_pair or released_bound_count == 1 or iteration == 0 or
                options.max_iterations - iteration <= scalar_polishing_iterations))
        {
            const limiting_phase = limiting_index / inventory_count;
            if (limiting_phase < 2) {
                const gas_index = limiting_index % inventory_count;
                const dissolved_index = inventory_count + gas_index;
                const pair_total_g = target[gas_index] + target[dissolved_index];
                const pair_is_depleted =
                    isScaleDepletedGas(base[gas_index], base[dissolved_index], options);
                if (pair_is_depleted and
                    current[dissolved_index] > 0 and
                    std.math.isFinite(pair_total_g) and pair_total_g > 0)
                {
                    // The conservative transport target itself changes with
                    // the trial mixture. Follow that physical map for the
                    // depleted limiting pair before forming its local Newton
                    // block; prioritizing this coordinate prevents a later
                    // dense step from undoing its active-set release.
                    const current_pair_norm = @max(
                        try scaledCoordinateResidual(current[gas_index], residual[gas_index], options),
                        try scaledCoordinateResidual(current[dissolved_index], residual[dissolved_index], options),
                    );
                    @memcpy(candidate, current);
                    candidate[gas_index] = target[gas_index];
                    candidate[dissolved_index] = target[dissolved_index];
                    if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                        const picard_pair_norm = @max(
                            try scaledCoordinateResidual(candidate[gas_index], candidate_residual[gas_index], options),
                            try scaledCoordinateResidual(candidate[dissolved_index], candidate_residual[dissolved_index], options),
                        );
                        if (picard_pair_norm < current_pair_norm) {
                            @memcpy(previous, current);
                            @memcpy(previous_residual, residual);
                            has_previous = true;
                            @memcpy(current, candidate);
                            picard_steps += 1;
                            continue;
                        }
                    } else |_| {}
                    // Restore F(current) after a rejected Picard probe because
                    // residualAt writes the shared assembled target.
                    try residualAt(allocator, &scratch, base, current, inputs, residual_fraction, target, residual);
                    // Use a roundoff-scale absolute floor: dry-cell gas inventories can
                    // be near zero while pressure displacement has gram-scale sensitivity.
                    const gas_probe = std.math.cbrt(std.math.floatEps(f64)) *
                        @max(
                            options.absolute_tolerance_g,
                            @max(
                                @abs(current[gas_index]),
                                @abs(target[gas_index]),
                            ),
                        );
                    const dissolved_probe = std.math.cbrt(std.math.floatEps(f64)) *
                        @max(
                            options.absolute_tolerance_g,
                            @max(
                                @abs(current[dissolved_index]),
                                @abs(target[dissolved_index]),
                            ),
                        );
                    @memcpy(probe, current);
                    probe[gas_index] += gas_probe;
                    if (residualAt(allocator, &scratch, base, probe, inputs, residual_fraction, target, probe_residual)) |_| {
                        const j00 = (probe_residual[gas_index] - residual[gas_index]) / gas_probe;
                        const j10 = (probe_residual[dissolved_index] - residual[dissolved_index]) / gas_probe;
                        @memcpy(probe, current);
                        probe[dissolved_index] += dissolved_probe;
                        if (residualAt(allocator, &scratch, base, probe, inputs, residual_fraction, target, candidate_residual)) |_| {
                            const j01 = (candidate_residual[gas_index] - residual[gas_index]) / dissolved_probe;
                            const j11 = (candidate_residual[dissolved_index] - residual[dissolved_index]) / dissolved_probe;
                            const determinant = j00 * j11 - j01 * j10;
                            if (std.math.isFinite(determinant) and @abs(determinant) > std.math.floatEps(f64)) {
                                const gas_direction = (-residual[gas_index] * j11 + j01 * residual[dissolved_index]) / determinant;
                                const dissolved_direction = (j10 * residual[gas_index] - j00 * residual[dissolved_index]) / determinant;
                                var pair_fraction: f64 = 1;
                                var pair_line_search: u8 = 0;
                                var accepted_pair_newton = false;
                                while (pair_line_search < 24) : (pair_line_search += 1) {
                                    @memcpy(candidate, current);
                                    candidate[gas_index] = @max(0, current[gas_index] + pair_fraction * gas_direction);
                                    candidate[dissolved_index] = @max(0, current[dissolved_index] + pair_fraction * dissolved_direction);
                                    if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                                        const candidate_pair_norm = @max(
                                            try scaledCoordinateResidual(candidate[gas_index], candidate_residual[gas_index], options),
                                            try scaledCoordinateResidual(candidate[dissolved_index], candidate_residual[dissolved_index], options),
                                        );
                                        if (candidate_pair_norm < current_pair_norm) {
                                            @memcpy(previous, current);
                                            @memcpy(previous_residual, residual);
                                            has_previous = true;
                                            @memcpy(current, candidate);
                                            newton_steps += 1;
                                            accepted_pair_newton = true;
                                            break;
                                        }
                                    } else |_| {}
                                    pair_fraction *= 0.5;
                                }
                                if (accepted_pair_newton) continue;
                            }
                        } else |_| {}
                    } else |_| {}
                    @memcpy(candidate, current);
                    candidate[gas_index] = 0;
                    candidate[dissolved_index] = pair_total_g;
                    if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                        var lower_g: f64 = 0;
                        var lower_residual = candidate_residual[gas_index];
                        candidate[gas_index] = pair_total_g;
                        candidate[dissolved_index] = 0;
                        if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                            var upper_g = pair_total_g;
                            var upper_residual = candidate_residual[gas_index];
                            if (lower_residual == 0 or upper_residual == 0 or std.math.signbit(lower_residual) != std.math.signbit(upper_residual)) {
                                var pair_search: u8 = 0;
                                while (pair_search < 60) : (pair_search += 1) {
                                    const midpoint_g = lower_g + 0.5 * (upper_g - lower_g);
                                    @memcpy(candidate, current);
                                    candidate[gas_index] = midpoint_g;
                                    candidate[dissolved_index] = pair_total_g - midpoint_g;
                                    try residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual);
                                    const midpoint_residual = candidate_residual[gas_index];
                                    if (std.math.signbit(midpoint_residual) == std.math.signbit(lower_residual)) {
                                        lower_g = midpoint_g;
                                        lower_residual = midpoint_residual;
                                    } else {
                                        upper_g = midpoint_g;
                                        upper_residual = midpoint_residual;
                                    }
                                    if (@abs(midpoint_residual) <= options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, midpoint_g)) break;
                                }
                                const candidate_pair_norm = @max(
                                    try scaledCoordinateResidual(candidate[gas_index], candidate_residual[gas_index], options),
                                    try scaledCoordinateResidual(candidate[dissolved_index], candidate_residual[dissolved_index], options),
                                );
                                if (candidate_pair_norm < current_pair_norm and
                                    std.math.isFinite(try scaledNorm(candidate, candidate_residual, options)))
                                {
                                    @memcpy(previous, current);
                                    @memcpy(previous_residual, residual);
                                    has_previous = true;
                                    @memcpy(current, candidate);
                                    newton_steps += 1;
                                    continue;
                                }
                            }
                        } else |_| {}
                    } else |_| {}
                }
            }
            @memcpy(candidate, current);
            candidate[limiting_index] = 0;
            if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                const zero_residual = candidate_residual[limiting_index];
                if (zero_residual != 0) {
                    var lower: f64 = 0;
                    var lower_residual = zero_residual;
                    var upper = current[limiting_index];
                    var upper_residual = residual[limiting_index];
                    var bracket_search: u8 = 0;
                    while (bracket_search < 40 and std.math.signbit(lower_residual) == std.math.signbit(upper_residual)) : (bracket_search += 1) {
                        upper = @max(2 * upper, @max(target[limiting_index], @abs(residual[limiting_index])));
                        candidate[limiting_index] = upper;
                        try residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual);
                        upper_residual = candidate_residual[limiting_index];
                    }
                    if (std.math.signbit(lower_residual) != std.math.signbit(upper_residual)) {
                        var root_residual = zero_residual;
                        var root_search: u8 = 0;
                        while (root_search < 60) : (root_search += 1) {
                            const midpoint = lower + 0.5 * (upper - lower);
                            candidate[limiting_index] = midpoint;
                            try residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual);
                            root_residual = candidate_residual[limiting_index];
                            if (std.math.signbit(root_residual) == std.math.signbit(lower_residual)) {
                                lower = midpoint;
                                lower_residual = root_residual;
                            } else {
                                upper = midpoint;
                                upper_residual = root_residual;
                            }
                            if (@abs(root_residual) <= options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, midpoint)) break;
                        }
                        const old_coordinate = try scaledCoordinateResidual(current[limiting_index], residual[limiting_index], options);
                        const new_coordinate = try scaledCoordinateResidual(candidate[limiting_index], candidate_residual[limiting_index], options);
                        if (new_coordinate < old_coordinate and try scaledNorm(candidate, candidate_residual, options) <= norm) {
                            @memcpy(previous, current);
                            @memcpy(previous_residual, residual);
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            continue;
                        }
                    }
                }
            } else |_| {}
        }
        var accepted_newton = false;
        // Every rejected probe shares `target` scratch. Reassemble F(current)
        // before selecting the coupled species/face direction so active-set
        // tests never pair the current residual with a stale candidate
        // target.
        try residualAt(
            allocator,
            &scratch,
            base,
            current,
            inputs,
            residual_fraction,
            target,
            residual,
        );
        norm = try scaledNorm(current, residual, options);
        if (norm <= 1) {
            if (capturesFluxLedgers(inputs))
                try residualAtCapturing(
                    allocator,
                    &scratch,
                    base,
                    current,
                    inputs,
                    residual_fraction,
                    target,
                    residual,
                    true,
                );
            copyVectorToState(target, state);
            return .{
                .iterations = iteration + 1,
                .newton_raphson_steps = newton_steps,
                .picard_steps = picard_steps,
                .maximum_scaled_residual = norm,
            };
        }
        // Every transport and phase-exchange equation is species-separable.
        // After the first global Picard predictor, solve the species that
        // controls the infinity norm directly. Waiting for a late "tail"
        // wastes the NPH*NPG budget on unrelated CO2/O2 blocks and can make
        // otherwise convergent NPG-local solves exhaust their hourly ceiling.
        const active_species: ?usize = if (iteration >= 1)
            try worstResidualSpecies(current, residual, options, inventory_count)
        else
            null;
        // Resolve the convergence-controlling donor clamp before a
        // species-wide Picard update can report progress in unrelated cells
        // and skip this conservative face coordinate.
        if (iteration >= 1 and
            options.max_iterations - iteration <= 16 and
            try conservativeMultiSpeciesFaceNewtonStep(
                allocator,
                &scratch,
                base,
                current,
                residual,
                inputs,
                residual_fraction,
                target,
                probe,
                probe_residual,
                candidate,
                candidate_residual,
                previous,
                previous_residual,
                options,
                norm,
            ))
        {
            has_previous = true;
            newton_steps += 1;
            continue;
        }
        // Restore F(current) after rejected multi-species probes before the
        // scalar face safeguard inspects its active-set target.
        try residualAt(
            allocator,
            &scratch,
            base,
            current,
            inputs,
            residual_fraction,
            target,
            residual,
        );
        if (iteration >= 1 and try conservativeFaceNewtonStep(
            allocator,
            &scratch,
            base,
            current,
            residual,
            inputs,
            residual_fraction,
            target,
            probe,
            probe_residual,
            candidate,
            candidate_residual,
            previous,
            previous_residual,
            options,
            norm,
        )) {
            has_previous = true;
            newton_steps += 1;
            continue;
        }
        // Near the root, active transport clamps can reduce a valid dense
        // Newton direction to a slowly contracting sequence. Apply the
        // existing coupled Anderson extrapolation before accepting another
        // weak dense step; unlike scalar polishing this advances all tied
        // cells and phases in one outer iteration.
        if (iteration >= options.max_iterations / 2 and has_previous) {
            var numerator: f64 = 0;
            var denominator: f64 = 0;
            for (residual, previous_residual) |difference, old_difference| {
                const change = difference - old_difference;
                numerator += difference * change;
                denominator += change * change;
            }
            var valid_tail_anderson = std.math.isFinite(denominator) and denominator > std.math.floatEps(f64);
            if (valid_tail_anderson) {
                const mixing = numerator / denominator;
                for (current, previous, residual, previous_residual, candidate) |value, old_value, difference, old_difference, *next| {
                    const fixed_point = value + difference;
                    const old_fixed_point = old_value + old_difference;
                    next.* = fixed_point - mixing * (fixed_point - old_fixed_point);
                    if (!std.math.isFinite(next.*)) {
                        valid_tail_anderson = false;
                        break;
                    }
                    next.* = @max(0, next.*);
                }
            }
            if (valid_tail_anderson) {
                if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        @memcpy(previous, current);
                        @memcpy(previous_residual, residual);
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        continue;
                    }
                } else |_| {}
            }
        }
        if (iteration >= 1) {
            var iteration_options = options;
            iteration_options.transport_iteration_fraction = residual_fraction;
            const species_block_dimension = try std.math.mul(usize, scratch.cell_count, 3);
            // For a small complete system, the full Jacobian costs the same
            // number of residual probes as seven separately assembled
            // species blocks and additionally captures pressure displacement
            // cross derivatives. Large systems retain bounded species/Krylov
            // storage.
            // Pressure displacement couples species even though molecular
            // diffusion and phase exchange are species-local. Reserve a
            // complete-Jacobian window for a bounded system that remains
            // materially unconverged. Ordinary hours keep the cheaper
            // species blocks, while the hard tail gets enough coupled steps
            // to resolve several comparable species instead of waiting until
            // only eight iterations remain.
            const unconverged_species_count = try unconvergedSpeciesCount(
                current,
                residual,
                options,
                inventory_count,
            );
            const use_full_dense =
                unknown_count <= 256 and
                (options.max_iterations - iteration <= 8 or
                    (norm > 500 and
                        options.max_iterations - iteration <= 16) or
                    (unconverged_species_count >= 2 and
                        norm > 500 and
                        options.max_iterations - iteration <= 24) or
                    (unconverged_species_count >= 4 and
                        norm > 100 and
                        options.max_iterations - iteration <= 24));
            const use_single_species_block =
                active_species != null and
                unconverged_species_count == 1;
            const use_all_species_blocks =
                !use_single_species_block and
                !use_full_dense and
                species_block_dimension <= 256;
            @memset(probe, 0);
            var has_direction = if (use_single_species_block)
                try denseSpeciesNewtonDirection(
                    allocator,
                    &scratch,
                    base,
                    current,
                    residual,
                    inputs,
                    iteration_options,
                    target,
                    candidate_residual,
                    candidate,
                    probe,
                    active_species.?,
                )
            else if (use_full_dense)
                try denseFullNewtonDirection(allocator, &scratch, base, current, residual, inputs, iteration_options, target, candidate_residual, candidate, probe)
            else if (use_all_species_blocks)
                try denseAllSpeciesNewtonDirection(allocator, &scratch, base, current, residual, inputs, iteration_options, target, candidate_residual, candidate, probe)
            else if (active_species) |species|
                (try denseSpeciesNewtonDirection(allocator, &scratch, base, current, residual, inputs, iteration_options, target, candidate_residual, candidate, probe, species)) or
                    try krylovNewtonDirection(allocator, &scratch, base, current, residual, inputs, iteration_options, target, probe_residual, probe, active_species)
            else
                try krylovNewtonDirection(allocator, &scratch, base, current, residual, inputs, iteration_options, target, probe_residual, probe, active_species);
            if (has_direction and use_all_species_blocks) has_direction = try filterIndependentSpeciesDirections(
                allocator,
                &scratch,
                base,
                current,
                residual,
                inputs,
                iteration_options,
                target,
                candidate_residual,
                candidate,
                probe,
            );
            if (has_direction) {
                var line_fraction: f64 = 1;
                var line_search: u8 = 0;
                while (line_search < 20) : (line_search += 1) {
                    var valid_line = true;
                    for (current, probe, candidate) |value, direction, *next| {
                        next.* = value + line_fraction * direction;
                        if (!std.math.isFinite(next.*)) {
                            valid_line = false;
                            break;
                        }
                        // Gas inventories are nonnegative complementarity
                        // variables. Project the Newton trial onto the active
                        // bound instead of rejecting an otherwise useful
                        // coupled direction because one depleted pool crosses
                        // zero.
                        next.* = @max(0, next.*);
                    }
                    if (valid_line) {
                        if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                            // Positivity projection can activate a zero bound
                            // whose conservative assembled target is positive.
                            // Close that active set inside this same Newton
                            // iteration; otherwise the final projected step
                            // can consume the ceiling and leave a positive
                            // source unresolved until a nonexistent next
                            // outer iteration.
                            var released_projected_bound = false;
                            for (candidate, candidate_residual, target) |*value, difference, assembled| {
                                if (isNumericallyAtNonnegativeBound(
                                    value.*,
                                    assembled,
                                ) and difference > 0) {
                                    value.* = assembled;
                                    released_projected_bound = true;
                                }
                            }
                            if (released_projected_bound)
                                try residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual);
                            // The termination test is an infinity norm. Using
                            // an L2-only merit function can accept steps that
                            // reduce large gas blocks while leaving a small
                            // species (commonly N2O) above tolerance for the
                            // entire NPH*NPG budget.
                            const improves = if (use_all_species_blocks or use_full_dense)
                                // Dense block Newton is already protected by
                                // positivity projection and line search. Keep
                                // strong globalization away from the root,
                                // but accept strict progress once the residual
                                // is within twice the requested tolerance.
                                try newtonMeritImproves(candidate, candidate_residual, options, newtonAcceptanceTarget(norm))
                            else if (active_species != null)
                                try speciesBlockMeritImproves(current, residual, candidate, candidate_residual, options, inventory_count, active_species.?, norm)
                            else
                                try newtonMeritImproves(candidate, candidate_residual, options, norm);
                            if (improves) {
                                @memcpy(previous, current);
                                @memcpy(previous_residual, residual);
                                has_previous = true;
                                @memcpy(current, candidate);
                                newton_steps += 1;
                                accepted_newton = true;
                                break;
                            }
                        } else |_| {}
                    }
                    line_fraction *= 0.5;
                }
            }
        }
        if (accepted_newton) continue;
        if (active_species) |species| {
            @memset(probe, 0);
            var species_options = options;
            species_options.transport_iteration_fraction = residual_fraction;
            const has_species_newton_direction = try denseSpeciesNewtonDirection(allocator, &scratch, base, current, residual, inputs, species_options, target, candidate_residual, candidate, probe, species);
            if (has_species_newton_direction) {
                var newton_fraction: f64 = 1;
                var newton_search: u8 = 0;
                while (newton_search < 20) : (newton_search += 1) {
                    for (current, probe, candidate) |value, direction, *next| next.* = @max(0, value + newton_fraction * direction);
                    if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                        const current_species_norm = try scaledSpeciesNorm(current, residual, options, inventory_count, species);
                        const candidate_species_norm = try scaledSpeciesNorm(candidate, candidate_residual, options, inventory_count, species);
                        if (candidate_species_norm < newtonAcceptanceTarget(current_species_norm) and try speciesBlockMeritImproves(current, residual, candidate, candidate_residual, options, inventory_count, species, norm)) {
                            @memcpy(previous, current);
                            @memcpy(previous_residual, residual);
                            has_previous = true;
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                            break;
                        }
                    } else |_| {}
                    newton_fraction *= 0.5;
                }
            }
            if (accepted_newton) continue;
            if (has_previous) {
                var numerator: f64 = 0;
                var denominator: f64 = 0;
                for (residual, previous_residual, 0..) |difference, old_difference, index| {
                    if ((index % inventory_count) % gas.species_count != species) continue;
                    const change = difference - old_difference;
                    numerator += difference * change;
                    denominator += change * change;
                }
                if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                    const mixing = numerator / denominator;
                    @memcpy(candidate, current);
                    var valid_species_anderson = std.math.isFinite(mixing);
                    if (valid_species_anderson) for (candidate, current, previous, residual, previous_residual, 0..) |*next, value, old_value, difference, old_difference, index| {
                        if ((index % inventory_count) % gas.species_count != species) continue;
                        const fixed_point = value + difference;
                        const old_fixed_point = old_value + old_difference;
                        next.* = @max(0, fixed_point - mixing * (fixed_point - old_fixed_point));
                        if (!std.math.isFinite(next.*)) {
                            valid_species_anderson = false;
                            break;
                        }
                    };
                    if (valid_species_anderson) {
                        if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                            if (try speciesBlockMeritImproves(current, residual, candidate, candidate_residual, options, inventory_count, species, norm)) {
                                @memcpy(previous, current);
                                @memcpy(previous_residual, residual);
                                @memcpy(current, candidate);
                                newton_steps += 1;
                                accepted_newton = true;
                            }
                        } else |_| {}
                    }
                }
            }
            if (accepted_newton) continue;
            var species_fraction: f64 = 1;
            var species_search: u8 = 0;
            while (species_search < 24) : (species_search += 1) {
                @memcpy(candidate, current);
                for (candidate, residual, 0..) |*next, difference, index| {
                    if ((index % inventory_count) % gas.species_count == species)
                        next.* = @max(0, next.* + species_fraction * difference);
                }
                if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                    if (try speciesBlockMeritImproves(current, residual, candidate, candidate_residual, options, inventory_count, species, norm)) {
                        @memcpy(previous, current);
                        @memcpy(previous_residual, residual);
                        has_previous = true;
                        @memcpy(current, candidate);
                        picard_steps += 1;
                        accepted_newton = true;
                        break;
                    }
                } else |_| {}
                species_fraction *= 0.5;
            }
            // Preserve an accepted species-block correction without a
            // following global mix that could undo it.
            if (accepted_newton) continue;
        }
        // A phase bound can make the full species Jacobian singular even
        // though the one coordinate controlling the infinity norm remains
        // locally smooth. Polish that coordinate with a scalar Newton step
        // before falling back to global mixing. This stays inside the same
        // NPH*NPG iteration budget and never relaxes the requested tolerance.
        if (try coordinateNewtonStep(allocator, &scratch, base, current, residual, inputs, residual_fraction, target, probe, probe_residual, candidate, candidate_residual, previous, previous_residual, options, norm)) {
            has_previous = true;
            newton_steps += 1;
            accepted_newton = true;
        }
        if (accepted_newton) continue;
        if (has_previous) {
            var numerator: f64 = 0;
            var denominator: f64 = 0;
            for (residual, previous_residual) |difference, old_difference| {
                const change = difference - old_difference;
                numerator += difference * change;
                denominator += change * change;
            }
            var valid_anderson = std.math.isFinite(denominator) and denominator > std.math.floatEps(f64);
            if (valid_anderson) {
                const mixing = numerator / denominator;
                for (current, previous, residual, previous_residual, candidate) |value, old_value, difference, old_difference, *next| {
                    const fixed_point = value + difference;
                    const old_fixed_point = old_value + old_difference;
                    next.* = fixed_point - mixing * (fixed_point - old_fixed_point);
                    if (!std.math.isFinite(next.*)) {
                        valid_anderson = false;
                        break;
                    }
                    next.* = @max(0, next.*);
                }
            }
            if (valid_anderson) {
                if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        @memcpy(previous, current);
                        @memcpy(previous_residual, residual);
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        accepted_newton = true;
                    }
                } else |_| {}
            }
        }
        if (accepted_newton) continue;
        if (has_previous) {
            var valid_secant = true;
            for (current, previous, residual, previous_residual, candidate) |value, old_value, difference, old_difference, *next| {
                const denominator = difference - old_difference;
                if (std.math.isFinite(denominator) and @abs(denominator) > std.math.floatEps(f64) * @max(1.0, @abs(difference))) {
                    next.* = value - difference * (value - old_value) / denominator;
                } else {
                    next.* = value + options.picard_relaxation * difference;
                }
                if (!std.math.isFinite(next.*)) {
                    valid_secant = false;
                    break;
                }
                next.* = @max(0, next.*);
            }
            if (valid_secant) {
                if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        @memcpy(previous, current);
                        @memcpy(previous_residual, residual);
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        accepted_newton = true;
                    }
                } else |_| {}
            }
        }
        if (accepted_newton) continue;
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(allocator, &scratch, base, probe, inputs, residual_fraction, target, probe_residual)) |_| {
                var valid_candidate = true;
                for (current, residual, probe_residual, candidate) |value, difference, sampled_difference, *next| {
                    if (@abs(difference) <= std.math.floatEps(f64)) {
                        next.* = value;
                        continue;
                    }
                    const derivative = (sampled_difference - difference) /
                        (options.directional_probe_fraction * difference);
                    const raw_fraction = if (std.math.isFinite(derivative) and derivative < 0)
                        -1.0 / derivative
                    else
                        options.picard_relaxation;
                    const fraction = std.math.clamp(raw_fraction, options.minimum_newton_fraction, options.maximum_newton_fraction);
                    next.* = value + fraction * difference;
                    if (!std.math.isFinite(next.*) or next.* < -1e-12) {
                        valid_candidate = false;
                        break;
                    }
                    next.* = @max(0, next.*);
                }
                if (valid_candidate) {
                    if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                        if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                            @memcpy(previous, current);
                            @memcpy(previous_residual, residual);
                            has_previous = true;
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                        }
                    } else |_| {}
                }
            } else |_| {}
        } else |_| {}
        if (accepted_newton) continue;
        var picard_fraction = options.picard_relaxation;
        var accepted_picard = false;
        var backtrack: u8 = 0;
        while (backtrack < 24) : (backtrack += 1) {
            if (addDirection(current, residual, picard_fraction, candidate)) |_| {
                if (residualAt(allocator, &scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        accepted_picard = true;
                        break;
                    }
                } else |_| {}
            } else |_| {}
            picard_fraction *= 0.5;
        }
        if (!accepted_picard) {
            // A donor/positivity clamp can place the iterate exactly on a
            // nonsmooth active-set boundary where no monotone trial is
            // accepted. Take one small bounded Picard regularization step so
            // the following Krylov derivative is formed on the new branch.
            try addDirection(current, residual, @min(options.picard_relaxation, 1.0e-3), candidate);
        }
        // Do not scale stagnation by the largest inventory in the coupled
        // vector. A large CO2 pool can otherwise classify a representable,
        // convergence-controlling N2O/NH3 correction as zero. The residual is
        // already tolerance-scaled per coordinate; stagnation means that the
        // bounded correction changes no floating-point coordinate at all.
        if (vectorsEqual(current, candidate)) {
            std.log.warn("coupled gas solver stagnated: iteration={d} scaled_residual={e} newton_steps={d} picard_steps={d} unknowns={d}", .{ iteration + 1, norm, newton_steps, picard_steps, unknown_count });
            return error.CoupledGasSolverStagnated;
        }
        @memcpy(previous, current);
        @memcpy(previous_residual, residual);
        has_previous = true;
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    const semismooth_tail_budget =
        options.max_iterations - iteration;
    // Refresh diagnostics for the final accepted iterate; the loop's residual
    // otherwise describes the state before its last Newton/Picard update.
    try residualAt(allocator, &scratch, base, current, inputs, residual_fraction, target, residual);
    last_scaled_residual = try scaledNorm(current, residual, options);
    // The final projected Newton direction can activate a zero-inventory
    // bound whose conservative fixed-point target is positive. Normally the
    // next outer iteration releases that complementarity bound before
    // forming another Jacobian, but at the runtime iteration ceiling there
    // is no next iteration. Close every such bound as part of the already
    // completed final iteration, then require the same strict convergence
    // test. This is not an extra nonlinear iteration and cannot turn a
    // genuinely unconverged state into a successful solve.
    var released_final_bound = false;
    for (current, residual, target) |*value, difference, assembled| {
        if (isNumericallyAtNonnegativeBound(value.*, assembled) and
            difference > 0)
        {
            value.* = assembled;
            released_final_bound = true;
        }
    }
    if (released_final_bound) {
        try residualAt(allocator, &scratch, base, current, inputs, residual_fraction, target, residual);
        last_scaled_residual = try scaledNorm(current, residual, options);
    }
    // Reserve part of the user-provided NPH*NPG ceiling for a targeted
    // semismooth tail. Global Newton/Picard work and this face/species polish
    // therefore share one iteration budget instead of silently adding a
    // second budget after the ceiling.
    var final_active_set_closures: u16 = 0;
    while (final_active_set_closures < semismooth_tail_budget) : (final_active_set_closures += 1) {
        const closed_full_system =
            try denseFullNewtonStep(
                allocator,
                &scratch,
                base,
                current,
                residual,
                inputs,
                residual_fraction,
                target,
                probe,
                candidate,
                candidate_residual,
                previous,
                previous_residual,
                options,
                last_scaled_residual,
            );
        const closed_multi_species_face =
            if (closed_full_system)
                false
            else
                try conservativeMultiSpeciesFaceNewtonStep(
                    allocator,
                    &scratch,
                    base,
                    current,
                    residual,
                    inputs,
                    residual_fraction,
                    target,
                    probe,
                    probe_residual,
                    candidate,
                    candidate_residual,
                    previous,
                    previous_residual,
                    options,
                    last_scaled_residual,
                );
        const closed_face =
            if (closed_full_system or closed_multi_species_face)
                false
            else
                try conservativeFaceNewtonStep(
                    allocator,
                    &scratch,
                    base,
                    current,
                    residual,
                    inputs,
                    residual_fraction,
                    target,
                    probe,
                    probe_residual,
                    candidate,
                    candidate_residual,
                    previous,
                    previous_residual,
                    options,
                    last_scaled_residual,
                );
        const polished_all_species =
            if (closed_full_system or
            closed_multi_species_face or
            closed_face)
                false
            else
                try denseAllSpeciesNewtonStep(
                    allocator,
                    &scratch,
                    base,
                    current,
                    residual,
                    inputs,
                    residual_fraction,
                    target,
                    probe,
                    candidate,
                    candidate_residual,
                    previous,
                    previous_residual,
                    options,
                    last_scaled_residual,
                );
        const polished_species =
            if (closed_full_system or
            closed_multi_species_face or
            closed_face or
            polished_all_species)
                false
            else
                try denseWorstSpeciesNewtonStep(
                    allocator,
                    &scratch,
                    base,
                    current,
                    residual,
                    inputs,
                    residual_fraction,
                    target,
                    probe,
                    candidate,
                    candidate_residual,
                    previous,
                    previous_residual,
                    options,
                    last_scaled_residual,
                );
        const polished_coordinate =
            if (closed_full_system or
            closed_multi_species_face or
            closed_face or
            polished_all_species or
            polished_species)
                false
            else
                try coordinateNewtonStep(
                    allocator,
                    &scratch,
                    base,
                    current,
                    residual,
                    inputs,
                    residual_fraction,
                    target,
                    probe,
                    probe_residual,
                    candidate,
                    candidate_residual,
                    previous,
                    previous_residual,
                    options,
                    last_scaled_residual,
                );
        if (!closed_full_system and
            !closed_multi_species_face and
            !closed_face and
            !polished_all_species and
            !polished_species and
            !polished_coordinate) break;
        newton_steps += 1;
        try residualAt(
            allocator,
            &scratch,
            base,
            current,
            inputs,
            residual_fraction,
            target,
            residual,
        );
        last_scaled_residual =
            try scaledNorm(current, residual, options);
        if (last_scaled_residual <= 1) break;
    }
    if (last_scaled_residual <= 1) {
        if (capturesFluxLedgers(inputs)) try residualAtCapturing(allocator, &scratch, base, current, inputs, residual_fraction, target, residual, true);
        copyVectorToState(target, state);
        return .{ .iterations = options.max_iterations, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = last_scaled_residual };
    }
    var worst_index: usize = 0;
    var worst_scaled: f64 = 0;
    for (current, residual, 0..) |value, difference, index| {
        const scaled = @abs(difference) / (options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, @abs(value)));
        if (scaled > worst_scaled) {
            worst_scaled = scaled;
            worst_index = index;
        }
    }
    const phase = worst_index / inventory_count;
    const component = worst_index % inventory_count;
    const cell = component / gas.species_count;
    const species = component % gas.species_count;
    std.log.warn("coupled gas solver exhausted runtime ceiling: max_iterations={d} scaled_residual={e} newton_steps={d} picard_steps={d} phase={d} cell={d} species={d} mass_g={e} residual_g={e}", .{ options.max_iterations, last_scaled_residual, newton_steps, picard_steps, phase, cell, species, current[worst_index], residual[worst_index] });
    for (0..gas.species_count) |diagnostic_species| {
        std.log.warn(
            "coupled gas failure species norm: species={d} scaled_residual={e}",
            .{
                diagnostic_species,
                try scaledSpeciesNorm(
                    current,
                    residual,
                    options,
                    inventory_count,
                    diagnostic_species,
                ),
            },
        );
    }
    for (0..3) |diagnostic_phase| {
        const index = diagnostic_phase * inventory_count + component;
        std.log.warn("coupled gas failure phase state: phase={d} base_g={e} iterate_g={e} target_g={e} residual_g={e}", .{ diagnostic_phase, base[index], current[index], target[index], residual[index] });
    }
    for (inputs.faces, 0..) |face, face_index| {
        if (face.first_cell != cell and face.second_cell != cell) continue;
        const first_component = face.first_cell * gas.species_count + species;
        const second_component = face.second_cell * gas.species_count + species;
        std.log.warn("coupled gas failure face: face={d} first_cell={d} second_cell={d} conductance_m3={e} first_g={e} second_g={e}", .{
            face_index,
            face.first_cell,
            face.second_cell,
            inputs.face_conductance_m3_per_step[face_index * gas.species_count + species],
            current[first_component],
            current[second_component],
        });
    }
    return error.CoupledGasSolverDidNotConverge;
}

/// Solves all currently unconverged gaseous species on one conservative face
/// manifold. A transfer coordinate adds mass to the receiving cell and
/// removes exactly the same mass from its neighbor. The small dense Jacobian
/// therefore captures cross-species pressure displacement without changing
/// any species inventory or assembling the complete three-phase system.
fn conservativeMultiSpeciesFaceNewtonStep(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []f64,
    residual: []const f64,
    inputs: Inputs,
    residual_fraction: f64,
    target: []f64,
    probe: []f64,
    probe_residual: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    previous: []f64,
    previous_residual: []f64,
    options: Options,
    current_norm: f64,
) !bool {
    if (current_norm <= 10) return false;
    const inventory_count = scratch.gaseous_mass_g.len;
    const limiting_index =
        (try worstExtremeIncomingGasIndex(
            current,
            residual,
            target,
            inputs,
            options,
            inventory_count,
        )) orelse return false;
    if (current[limiting_index] <= 0 or
        residual[limiting_index] <= 0 or
        target[limiting_index] <= current[limiting_index])
        return false;
    const receiving_cell = limiting_index / gas.species_count;
    const limiting_species = limiting_index % gas.species_count;
    const active_species = try allocator.alloc(usize, gas.species_count);
    defer allocator.free(active_species);
    var dimension: usize = 0;
    for (0..gas.species_count) |species| {
        if (try scaledSpeciesNorm(
            current,
            residual,
            options,
            inventory_count,
            species,
        ) > 1) {
            active_species[dimension] = species;
            dimension += 1;
        }
    }
    if (dimension < 2) return false;

    for (inputs.faces) |face| {
        if (face.first_cell != receiving_cell and
            face.second_cell != receiving_cell)
            continue;
        const donor_cell =
            if (face.first_cell == receiving_cell)
                face.second_cell
            else
                face.first_cell;
        const limiting_donor =
            donor_cell * gas.species_count + limiting_species;
        if (current[limiting_donor] <=
            10.0 * current[limiting_index])
            continue;

        const matrix = try allocator.alloc(
            f64,
            try std.math.mul(usize, dimension, dimension),
        );
        defer allocator.free(matrix);
        const right_hand_side = try allocator.alloc(f64, dimension);
        defer allocator.free(right_hand_side);
        var valid_jacobian = true;
        for (active_species[0..dimension], 0..) |species, column| {
            const receiver = receiving_cell * gas.species_count + species;
            const donor = donor_cell * gas.species_count + species;
            // The donor can be many orders of magnitude larger than the
            // receiver. sqrt(epsilon)*donor then becomes a nonlocal,
            // gram-scale transfer and differentiates the wrong active-set
            // branch. Use receiver-scale truncation error, but never request
            // less than a comfortably representable subtraction at the
            // donor magnitude.
            const donor_resolution_g =
                64.0 * std.math.floatEps(f64) *
                @max(1.0, current[donor]);
            const receiver_probe_g =
                gasJacobianProbeG(current[receiver]);
            const epsilon = @min(
                current[donor],
                @max(donor_resolution_g, receiver_probe_g),
            );
            if (!std.math.isFinite(epsilon) or epsilon <= 0) {
                valid_jacobian = false;
                break;
            }
            @memcpy(probe, current);
            probe[receiver] += epsilon;
            probe[donor] -= epsilon;
            residualAt(
                allocator,
                scratch,
                base,
                probe,
                inputs,
                residual_fraction,
                target,
                probe_residual,
            ) catch {
                valid_jacobian = false;
                break;
            };
            for (active_species[0..dimension], 0..) |row_species, row| {
                const row_index =
                    receiving_cell * gas.species_count + row_species;
                matrix[row * dimension + column] =
                    (probe_residual[row_index] - residual[row_index]) /
                    epsilon;
            }
        }
        if (!valid_jacobian) continue;
        for (active_species[0..dimension], 0..) |species, row| {
            const index = receiving_cell * gas.species_count + species;
            right_hand_side[row] = -residual[index];
        }
        if (!numerics.solveDenseLinearSystem(
            matrix,
            right_hand_side,
            dimension,
        )) continue;
        @memset(probe, 0);
        for (active_species[0..dimension], right_hand_side) |species, change| {
            if (!std.math.isFinite(change)) {
                valid_jacobian = false;
                break;
            }
            const receiver = receiving_cell * gas.species_count + species;
            const donor = donor_cell * gas.species_count + species;
            probe[receiver] = change;
            probe[donor] = -change;
        }
        if (!valid_jacobian) continue;
        var fraction: f64 = 1;
        var search: u8 = 0;
        const current_limiting_norm =
            try scaledCoordinateResidual(
                current[limiting_index],
                residual[limiting_index],
                options,
            );
        while (search < 24) : (search += 1) {
            @memcpy(candidate, current);
            var nonnegative = true;
            for (candidate, probe) |*value, direction| {
                value.* += fraction * direction;
                if (!std.math.isFinite(value.*) or value.* < 0) {
                    nonnegative = false;
                    break;
                }
            }
            if (!nonnegative) {
                fraction *= 0.5;
                continue;
            }
            if (residualAt(
                allocator,
                scratch,
                base,
                candidate,
                inputs,
                residual_fraction,
                target,
                candidate_residual,
            )) |_| {
                const candidate_norm =
                    try scaledNorm(candidate, candidate_residual, options);
                const candidate_limiting_norm =
                    try scaledCoordinateResidual(
                        candidate[limiting_index],
                        candidate_residual[limiting_index],
                        options,
                    );
                if (candidate_limiting_norm < current_limiting_norm and
                    candidate_norm < current_norm)
                {
                    @memcpy(previous, current);
                    @memcpy(previous_residual, residual);
                    @memcpy(current, candidate);
                    return true;
                }
            } else |_| {}
            fraction *= 0.5;
        }
    }
    return false;
}

/// Resolves a donor-clamped gaseous face on its conservative coordinate.
/// Moving only the limiting cell can make the neighboring residual worse and
/// causes the projected scalar/full Newton paths to bounce between opposite
/// donor bounds. An equal-and-opposite pair correction follows the actual
/// face conservation manifold while the full residual/line search continues
/// to account for pressure displacement, other faces, and phase exchange.
fn conservativeFaceNewtonStep(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []f64,
    residual: []const f64,
    inputs: Inputs,
    residual_fraction: f64,
    target: []f64,
    probe: []f64,
    probe_residual: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    previous: []f64,
    previous_residual: []f64,
    options: Options,
    current_norm: f64,
) !bool {
    const inventory_count = scratch.gaseous_mass_g.len;
    const limiting_index =
        (try worstScaleSeparatedIncomingGasIndex(
            current,
            residual,
            target,
            inputs,
            options,
            inventory_count,
        )) orelse return false;
    // This scalar safeguard is reserved for an unresolved incoming pool whose
    // residual exceeds its complete receiving inventory. Less-extreme
    // scale-separated blocks stay in the coupled multi-species Newton/Picard
    // path so one face coordinate cannot preempt the global correction.
    if (current[limiting_index] <= 0 or
        residual[limiting_index] <= 0 or
        target[limiting_index] <= current[limiting_index])
        return false;
    const limiting_cell = limiting_index / gas.species_count;
    const species = limiting_index % gas.species_count;
    for (inputs.faces) |face| {
        if (face.first_cell != limiting_cell and
            face.second_cell != limiting_cell)
            continue;
        const neighbor_cell =
            if (face.first_cell == limiting_cell)
                face.second_cell
            else
                face.first_cell;
        const neighbor_index =
            neighbor_cell * gas.species_count + species;
        if (current[neighbor_index] <=
            10.0 * current[limiting_index])
            continue;
        var probe_transfer_g = gasJacobianProbeG(
            @max(current[limiting_index], current[neighbor_index]),
        );
        if (current[neighbor_index] >= probe_transfer_g) {
            // Positive transfer increases the limiting cell.
        } else if (current[limiting_index] >= probe_transfer_g) {
            probe_transfer_g = -probe_transfer_g;
        } else {
            continue;
        }
        @memcpy(probe, current);
        probe[limiting_index] += probe_transfer_g;
        probe[neighbor_index] -= probe_transfer_g;
        residualAt(
            allocator,
            scratch,
            base,
            probe,
            inputs,
            residual_fraction,
            target,
            probe_residual,
        ) catch continue;
        const derivative =
            (probe_residual[limiting_index] -
                residual[limiting_index]) /
            probe_transfer_g;
        if (!std.math.isFinite(derivative) or
            @abs(derivative) <= std.math.floatEps(f64))
            continue;
        const correction_g =
            std.math.clamp(
                -residual[limiting_index] / derivative,
                -current[limiting_index],
                current[neighbor_index],
            );
        if (!std.math.isFinite(correction_g) or correction_g == 0)
            continue;
        const requires_nonlocal_bracket =
            requiresNonlocalConservativeBracket(
                current[limiting_index],
                current[neighbor_index],
                options,
            );
        var fraction: f64 = 1;
        var search: u8 = 0;
        while (search < 24 and
            !(requires_nonlocal_bracket and correction_g < 0)) : (search += 1)
        {
            @memcpy(candidate, current);
            const accepted_transfer_g = fraction * correction_g;
            candidate[limiting_index] += accepted_transfer_g;
            candidate[neighbor_index] -= accepted_transfer_g;
            if (candidate[limiting_index] < 0 or
                candidate[neighbor_index] < 0)
            {
                fraction *= 0.5;
                continue;
            }
            if (residualAt(
                allocator,
                scratch,
                base,
                candidate,
                inputs,
                residual_fraction,
                target,
                candidate_residual,
            )) |_| {
                const current_coordinate_norm =
                    try scaledCoordinateResidual(
                        current[limiting_index],
                        residual[limiting_index],
                        options,
                    );
                const candidate_coordinate_norm =
                    try scaledCoordinateResidual(
                        candidate[limiting_index],
                        candidate_residual[limiting_index],
                        options,
                    );
                const candidate_global_norm =
                    try scaledNorm(
                        candidate,
                        candidate_residual,
                        options,
                    );
                const unresolved_positive_bound =
                    hasUnresolvedPositiveBound(
                        candidate[limiting_index],
                        target[limiting_index],
                        candidate_residual[limiting_index],
                    );
                // Changing a donor-clamped face also changes total pressure
                // and can temporarily expose a larger residual in another
                // species. Requiring monotone global infinity norm prevents
                // the semismooth active set from ever changing. Bound that
                // transient growth while still requiring strict progress in
                // the controlling species; ordinary Newton/Picard steps and
                // the final convergence test retain strict global merit.
                if (!unresolved_positive_bound and
                    candidate_coordinate_norm <
                        current_coordinate_norm and
                    candidate_global_norm <=
                        4.0 * current_norm)
                {
                    @memcpy(previous, current);
                    @memcpy(previous_residual, residual);
                    @memcpy(current, candidate);
                    return true;
                }
            } else |_| {}
            fraction *= 0.5;
        }
        // A nonlocal bracket is reserved for genuinely scale-separated
        // donor bounds. On ordinary faces the local conservative Newton
        // path above has better coupled convergence and should retain
        // control.
        if (!requires_nonlocal_bracket)
            continue;
        // The residual immediately above a donor bound can have the wrong
        // local slope: pressure displacement may initially increase the
        // receiving-cell residual even though a root exists farther along
        // the conservative face coordinate. In that semismooth case Newton
        // points back through zero and projection cycles. Expand a positive
        // transfer bracket geometrically, then bisect the same
        // equal-and-opposite coordinate. The search remains inside this one
        // outer NPH*NPG iteration and preserves face mass exactly.
        var lower_transfer_g: f64 = 0;
        var lower_residual_g = residual[limiting_index];
        var upper_transfer_g =
            @min(
                current[neighbor_index],
                @max(
                    gasJacobianProbeG(current[neighbor_index]),
                    @max(
                        target[limiting_index] -
                            current[limiting_index],
                        current[limiting_index],
                    ),
                ),
            );
        var bracket_search: u8 = 0;
        while (bracket_search < 60 and
            upper_transfer_g > lower_transfer_g) : (bracket_search += 1)
        {
            @memcpy(candidate, current);
            candidate[limiting_index] += upper_transfer_g;
            candidate[neighbor_index] -= upper_transfer_g;
            residualAt(
                allocator,
                scratch,
                base,
                candidate,
                inputs,
                residual_fraction,
                target,
                candidate_residual,
            ) catch break;
            const upper_residual_g =
                candidate_residual[limiting_index];
            if (upper_residual_g == 0 or
                std.math.signbit(upper_residual_g) !=
                    std.math.signbit(lower_residual_g))
            {
                var root_transfer_g = upper_transfer_g;
                var root_residual_g = upper_residual_g;
                var bisection: u8 = 0;
                while (bisection < 60 and root_residual_g != 0) : (bisection += 1) {
                    root_transfer_g =
                        0.5 *
                        (lower_transfer_g + upper_transfer_g);
                    @memcpy(candidate, current);
                    candidate[limiting_index] += root_transfer_g;
                    candidate[neighbor_index] -= root_transfer_g;
                    try residualAt(
                        allocator,
                        scratch,
                        base,
                        candidate,
                        inputs,
                        residual_fraction,
                        target,
                        candidate_residual,
                    );
                    root_residual_g =
                        candidate_residual[limiting_index];
                    if (std.math.signbit(root_residual_g) ==
                        std.math.signbit(lower_residual_g))
                    {
                        lower_transfer_g = root_transfer_g;
                        lower_residual_g = root_residual_g;
                    } else {
                        upper_transfer_g = root_transfer_g;
                    }
                    if (@abs(root_residual_g) <=
                        options.absolute_tolerance_g +
                            options.relative_tolerance *
                                @max(
                                    1.0,
                                    candidate[limiting_index],
                                ))
                        break;
                }
                const candidate_coordinate_norm =
                    try scaledCoordinateResidual(
                        candidate[limiting_index],
                        candidate_residual[limiting_index],
                        options,
                    );
                const current_coordinate_norm =
                    try scaledCoordinateResidual(
                        current[limiting_index],
                        residual[limiting_index],
                        options,
                    );
                const candidate_global_norm =
                    try scaledNorm(
                        candidate,
                        candidate_residual,
                        options,
                    );
                if (candidate_coordinate_norm <
                    current_coordinate_norm and
                    std.math.isFinite(candidate_global_norm))
                {
                    @memcpy(previous, current);
                    @memcpy(previous_residual, residual);
                    @memcpy(current, candidate);
                    return true;
                }
                break;
            }
            if (upper_transfer_g == current[neighbor_index])
                break;
            lower_transfer_g = upper_transfer_g;
            lower_residual_g = upper_residual_g;
            upper_transfer_g =
                @min(
                    current[neighbor_index],
                    2.0 * upper_transfer_g,
                );
        }
    }
    return false;
}

fn denseAllSpeciesNewtonStep(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []f64,
    residual: []const f64,
    inputs: Inputs,
    residual_fraction: f64,
    target: []f64,
    direction: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    previous: []f64,
    previous_residual: []f64,
    options: Options,
    current_norm: f64,
) !bool {
    @memset(direction, 0);
    var block_options = options;
    block_options.transport_iteration_fraction =
        residual_fraction;
    if (!try denseAllSpeciesNewtonDirection(
        allocator,
        scratch,
        base,
        current,
        residual,
        inputs,
        block_options,
        target,
        candidate_residual,
        candidate,
        direction,
    )) return false;
    if (!try filterIndependentSpeciesDirections(
        allocator,
        scratch,
        base,
        current,
        residual,
        inputs,
        block_options,
        target,
        candidate_residual,
        candidate,
        direction,
    )) return false;
    for (current, direction, candidate) |value, delta, *next| {
        next.* = @max(0, value + delta);
    }
    try residualAt(
        allocator,
        scratch,
        base,
        candidate,
        inputs,
        residual_fraction,
        target,
        candidate_residual,
    );
    if (try scaledNorm(
        candidate,
        candidate_residual,
        options,
    ) >= current_norm) return false;
    @memcpy(previous, current);
    @memcpy(previous_residual, residual);
    @memcpy(current, candidate);
    return true;
}

/// Applies the complete cross-species Jacobian as one bounded tail update.
/// The ordinary loop may hand a scale-separated donor active set to the
/// semismooth tail at half-budget, before its late full-dense window opens.
/// Keeping this wrapper in the shared tail budget ensures pressure-
/// displacement derivatives are still available without adding iterations.
fn denseFullNewtonStep(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []f64,
    residual: []const f64,
    inputs: Inputs,
    residual_fraction: f64,
    target: []f64,
    direction: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    previous: []f64,
    previous_residual: []f64,
    options: Options,
    current_norm: f64,
) !bool {
    @memset(direction, 0);
    var full_options = options;
    full_options.transport_iteration_fraction = residual_fraction;
    if (!try denseFullNewtonDirection(
        allocator,
        scratch,
        base,
        current,
        residual,
        inputs,
        full_options,
        target,
        candidate_residual,
        candidate,
        direction,
    )) return false;
    var fraction: f64 = 1;
    var search: u8 = 0;
    while (search < 24) : (search += 1) {
        var finite_nonnegative = true;
        for (current, direction, candidate) |value, delta, *next| {
            next.* = value + fraction * delta;
            if (!std.math.isFinite(next.*)) {
                finite_nonnegative = false;
                break;
            }
            next.* = @max(0, next.*);
        }
        if (finite_nonnegative) {
            if (residualAt(
                allocator,
                scratch,
                base,
                candidate,
                inputs,
                residual_fraction,
                target,
                candidate_residual,
            )) |_| {
                if (try scaledNorm(
                    candidate,
                    candidate_residual,
                    options,
                ) < current_norm) {
                    @memcpy(previous, current);
                    @memcpy(previous_residual, residual);
                    @memcpy(current, candidate);
                    return true;
                }
            } else |_| {}
        }
        fraction *= 0.5;
    }
    return false;
}

fn denseWorstSpeciesNewtonStep(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []f64,
    residual: []const f64,
    inputs: Inputs,
    residual_fraction: f64,
    target: []f64,
    direction: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    previous: []f64,
    previous_residual: []f64,
    options: Options,
    current_norm: f64,
) !bool {
    const inventory_count = scratch.gaseous_mass_g.len;
    const species = try worstResidualSpecies(
        current,
        residual,
        options,
        inventory_count,
    );
    @memset(direction, 0);
    var species_options = options;
    species_options.transport_iteration_fraction =
        residual_fraction;
    if (!try denseSpeciesNewtonDirection(
        allocator,
        scratch,
        base,
        current,
        residual,
        inputs,
        species_options,
        target,
        candidate_residual,
        candidate,
        direction,
        species,
    )) return false;
    var fraction: f64 = 1;
    var search: u8 = 0;
    while (search < 20) : (search += 1) {
        for (current, direction, candidate) |value, delta, *next| {
            next.* = @max(0, value + fraction * delta);
        }
        if (residualAt(
            allocator,
            scratch,
            base,
            candidate,
            inputs,
            residual_fraction,
            target,
            candidate_residual,
        )) |_| {
            if (try scaledNorm(
                candidate,
                candidate_residual,
                options,
            ) < current_norm) {
                @memcpy(previous, current);
                @memcpy(previous_residual, residual);
                @memcpy(current, candidate);
                return true;
            }
        } else |_| {}
        fraction *= 0.5;
    }
    return false;
}

fn coordinateNewtonStep(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []f64,
    residual: []const f64,
    inputs: Inputs,
    residual_fraction: f64,
    target: []f64,
    probe: []f64,
    probe_residual: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    previous: []f64,
    previous_residual: []f64,
    options: Options,
    current_norm: f64,
) !bool {
    const limiting_index = try worstResidualIndex(current, residual, options);
    // In the convergence tail, stay within the current semismooth active set;
    // a larger cbrt(epsilon) probe can cross a donor clamp and differentiate
    // the neighboring branch instead of the accepted one.
    const coordinate_scale = @max(1.0, @abs(current[limiting_index]));
    const stable_probe = @sqrt(std.math.floatEps(f64)) * coordinate_scale;
    const cancellation_floor = 64.0 * std.math.floatEps(f64) * coordinate_scale;
    // Do not shrink the probe with the residual. Near convergence that made
    // the function difference smaller than the roundoff accumulated by the
    // coupled transport/phase evaluation, producing a biased derivative and
    // dozens of tiny accepted Newton steps. sqrt(epsilon)*state_scale is the
    // standard forward-difference balance between truncation and cancellation.
    const epsilon = @max(cancellation_floor, stable_probe);
    @memcpy(probe, current);
    probe[limiting_index] += epsilon;
    var has_derivative = false;
    var derivative: f64 = 0;
    if (residualAt(allocator, scratch, base, probe, inputs, residual_fraction, target, probe_residual)) |_| {
        if (current[limiting_index] >= epsilon) {
            @memcpy(candidate, current);
            candidate[limiting_index] -= epsilon;
            if (residualAt(allocator, scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
                derivative = (probe_residual[limiting_index] - candidate_residual[limiting_index]) / (2.0 * epsilon);
                has_derivative = std.math.isFinite(derivative) and @abs(derivative) > std.math.floatEps(f64);
            } else |_| {}
        } else {
            derivative = (probe_residual[limiting_index] - residual[limiting_index]) / epsilon;
            has_derivative = std.math.isFinite(derivative) and @abs(derivative) > std.math.floatEps(f64);
        }
    } else |_| {}
    if (!has_derivative) return false;

    const correction = -residual[limiting_index] / derivative;
    const coordinate_norm = try scaledCoordinateResidual(current[limiting_index], residual[limiting_index], options);
    var fraction: f64 = 1;
    var search: u8 = 0;
    while (search < 20) : (search += 1) {
        @memcpy(candidate, current);
        candidate[limiting_index] = @max(0, current[limiting_index] + fraction * correction);
        if (residualAt(allocator, scratch, base, candidate, inputs, residual_fraction, target, candidate_residual)) |_| {
            const next_coordinate_norm = try scaledCoordinateResidual(candidate[limiting_index], candidate_residual[limiting_index], options);
            if (next_coordinate_norm < coordinate_norm and try scaledNorm(candidate, candidate_residual, options) <= current_norm) {
                @memcpy(previous, current);
                @memcpy(previous_residual, residual);
                @memcpy(current, candidate);
                return true;
            }
        } else |_| {}
        fraction *= 0.5;
    }
    return false;
}

/// A zero gas pool still participates in pressure, phase, and active-bound
/// derivatives. Probing it at `cbrt(epsilon) * 1e-8` is below cancellation
/// noise when the same residual assembly contains kilogram-scale pools.
/// A unit-scale floor and square-root relative step retain meaningful columns
/// without making the perturbation a large fraction of trace-gas inventories.
fn gasJacobianProbeG(mass_g: f64) f64 {
    return @sqrt(std.math.floatEps(f64)) *
        @max(1.0, @abs(mass_g));
}

fn residualAt(allocator: std.mem.Allocator, scratch: *gas.State, base: []const f64, trial: []const f64, inputs: Inputs, transport_iteration_fraction: f64, target: []f64, residual: []f64) !void {
    return residualAtCapturing(allocator, scratch, base, trial, inputs, transport_iteration_fraction, target, residual, false);
}

fn capturesFluxLedgers(inputs: Inputs) bool {
    return inputs.atmospheric_flux_g_by_component != null or
        inputs.subsurface_flux_g_by_component != null or
        inputs.face_flux_g_by_component != null;
}

fn residualAtCapturing(allocator: std.mem.Allocator, scratch: *gas.State, base: []const f64, trial: []const f64, inputs: Inputs, transport_iteration_fraction: f64, target: []f64, residual: []f64, capture_boundaries: bool) !void {
    copyVectorToState(trial, scratch);
    @memcpy(target, base);
    const n = scratch.gaseous_mass_g.len;
    const gas_target = target[0..n];
    const dissolved_target = target[n .. 2 * n];
    const band_target = target[2 * n .. 3 * n];
    if (capture_boundaries) {
        if (inputs.atmospheric_flux_g_by_component) |fluxes| {
            if (fluxes.len != n) return error.GasBoundaryFluxSizeMismatch;
            @memset(fluxes, 0);
        }
        if (inputs.subsurface_flux_g_by_component) |subsurface_fluxes| {
            if (subsurface_fluxes.len != n) return error.GasBoundaryFluxSizeMismatch;
            @memset(subsurface_fluxes, 0);
        }
        if (inputs.face_flux_g_by_component) |face_fluxes| {
            if (face_fluxes.len !=
                inputs.faces.len * gas.species_count)
                return error.GasFaceFluxSizeMismatch;
            @memset(face_fluxes, 0);
        }
    }
    var diffusion: [gas.species_count]f64 = undefined;
    var pressure: [gas.species_count]f64 = undefined;
    // TRNSFR pressure displacement is a sequential donor-bounded inventory
    // correction, not a constitutive equilibrium. Assemble it explicitly from
    // the conservative target before solving the differentiable diffusion and
    // phase-exchange terms.
    for (inputs.faces, 0..) |face, face_index| {
        const first_target = gas_target[face.first_cell * gas.species_count ..][0..gas.species_count];
        const second_target = gas_target[face.second_cell * gas.species_count ..][0..gas.species_count];
        try gas.adjacentPressureDrivenFluxesG(first_target, second_target, scratch.air_volume_m3[face.second_cell], scratch.temperature_k[face.second_cell], scratch.water_vapor_mol[face.second_cell], transport_iteration_fraction, &pressure);
        for (pressure, 0..) |convective, species| {
            const first_index = face.first_cell * gas.species_count + species;
            const second_index = face.second_cell * gas.species_count + species;
            const bounded = std.math.clamp(convective, -gas_target[second_index], gas_target[first_index]);
            gas_target[first_index] -= bounded;
            gas_target[second_index] += bounded;
            if (capture_boundaries) {
                if (inputs.face_flux_g_by_component) |fluxes|
                    fluxes[face_index * gas.species_count + species] += bounded;
            }
        }
    }
    for (inputs.faces, 0..) |face, face_index| {
        const face_conductance = inputs.face_conductance_m3_per_step[face_index * gas.species_count ..][0..gas.species_count];
        var iteration_conductance_m3: [gas.species_count]f64 = undefined;
        for (face_conductance, 0..) |conductance_m3, species| {
            iteration_conductance_m3[species] = conductance_m3 * transport_iteration_fraction;
        }
        try gas.calculateFaceDiffusiveFluxesG(scratch, face, &iteration_conductance_m3, &diffusion);
        for (diffusion, 0..) |diffusive, species| {
            const first_index = face.first_cell * gas.species_count + species;
            const second_index = face.second_cell * gas.species_count + species;
            const bounded = std.math.clamp(diffusive, -gas_target[second_index], gas_target[first_index]);
            gas_target[first_index] -= bounded;
            gas_target[second_index] += bounded;
            if (capture_boundaries) {
                if (inputs.face_flux_g_by_component) |fluxes|
                    fluxes[face_index * gas.species_count + species] += bounded;
            }
        }
    }
    var boundary_flux: [gas.species_count]f64 = undefined;
    for (inputs.atmospheric_boundaries) |boundary| {
        try atmosphere.calculateFluxesG(scratch, boundary, transport_iteration_fraction, &boundary_flux);
        for (boundary_flux, 0..) |flux, species| {
            const index = boundary.cell_index * gas.species_count + species;
            const accepted = @max(-gas_target[index], flux);
            gas_target[index] += accepted;
            if (capture_boundaries) {
                if (inputs.atmospheric_flux_g_by_component) |fluxes| fluxes[index] += accepted;
            }
        }
    }
    for (inputs.subsurface_boundaries) |boundary| {
        try atmosphere.calculateFluxesG(scratch, boundary, transport_iteration_fraction, &boundary_flux);
        for (boundary_flux, 0..) |flux, species| {
            const index = boundary.cell_index * gas.species_count + species;
            const accepted = @max(-gas_target[index], flux);
            gas_target[index] += accepted;
            if (capture_boundaries) {
                if (inputs.subsurface_flux_g_by_component) |fluxes| fluxes[index] += accepted;
            }
        }
    }
    var bubble: [gas.species_count]f64 = undefined;
    var band_bubble: [gas.species_count]f64 = undefined;
    for (0..scratch.cell_count) |cell| {
        const start = cell * gas.species_count;
        const solubility = inputs.mass_solubility_ratio[start..][0..gas.species_count];
        for (0..gas.species_count) |species| {
            const index = start + species;
            // Phase exchange target uses base dissolved (dissolved_target before this loop
            // body modifies it) as the Picard departure point so that the fixed-point
            // equation is target_d = base_d + rate*(D_eq - base_d), which converges to
            // D_eq for rate=1.  Using trial dissolved instead gives a fixed point of
            // (base_d + rate*D_eq)/(1+rate) = D_eq/2 when base_d=0, causing oscillation.
            const requested_phase = try gas.phaseExchangeFluxG(scratch.gaseous_mass_g[index], dissolved_target[index], scratch.air_volume_m3[cell], inputs.water_volume_m3[cell], solubility[species], inputs.gas_water_exchange_rate_per_step[index]);
            const phase = std.math.clamp(requested_phase, -dissolved_target[index], gas_target[index]);
            const requested_band_phase = try gas.phaseExchangeFluxG(scratch.gaseous_mass_g[index] - @min(scratch.gaseous_mass_g[index], @max(0, phase)), band_target[index], scratch.air_volume_m3[cell], inputs.band_water_volume_m3[cell], solubility[species], inputs.band_gas_water_exchange_rate_per_step[index]);
            const band_phase = std.math.clamp(requested_band_phase, -band_target[index], gas_target[index] - phase);
            gas_target[index] -= phase + band_phase;
            dissolved_target[index] += phase;
            band_target[index] += band_phase;
        }
        if (inputs.bubbling_enabled[cell]) {
            try gas.bubblingFluxesG(inputs.water_volume_m3[cell], scratch.temperature_k[cell], scratch.dissolved_mass_g[start..][0..gas.species_count], solubility, transport_iteration_fraction, &bubble);
            try gas.bubblingFluxesG(inputs.band_water_volume_m3[cell], scratch.temperature_k[cell], scratch.band_dissolved_mass_g[start..][0..gas.species_count], solubility, transport_iteration_fraction, &band_bubble);
            for (bubble, band_bubble, 0..) |nonband, band, species| {
                const accepted_nonband = @max(-dissolved_target[start + species], nonband);
                const accepted_band = @max(-band_target[start + species], band);
                dissolved_target[start + species] += accepted_nonband;
                band_target[start + species] += accepted_band;
                const released_g = -(accepted_nonband + accepted_band);
                if (released_g == 0) continue;
                const receiver = if (inputs.bubble_receiver_cell_by_cell) |receivers|
                    receivers[cell]
                else
                    cell;
                if (receiver) |receiver_cell| {
                    gas_target[receiver_cell * gas.species_count + species] += released_g;
                } else if (capture_boundaries) {
                    if (inputs.subsurface_flux_g_by_component) |fluxes|
                        fluxes[start + species] -= released_g;
                }
            }
        }
    }
    _ = allocator;
    for (target, trial, residual) |fixed_point, value, *difference| {
        if (!std.math.isFinite(fixed_point) or fixed_point < -1e-12) return error.InvalidCoupledGasCandidate;
        difference.* = fixed_point - value;
    }
}

fn denseFullNewtonDirection(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    inputs: Inputs,
    options: Options,
    target: []f64,
    sampled_residual: []f64,
    sampled_state: []f64,
    direction: []f64,
) !bool {
    const dimension = current.len;
    if (dimension == 0 or dimension > 256) return false;
    const matrix = try allocator.alloc(f64, try std.math.mul(usize, dimension, dimension));
    defer allocator.free(matrix);
    const rhs = try allocator.alloc(f64, dimension);
    defer allocator.free(rhs);
    const normal_matrix = try allocator.alloc(
        f64,
        try std.math.mul(usize, dimension, dimension),
    );
    defer allocator.free(normal_matrix);
    const normal_rhs = try allocator.alloc(f64, dimension);
    defer allocator.free(normal_rhs);
    const opposite_state = try allocator.alloc(f64, dimension);
    defer allocator.free(opposite_state);
    const opposite_residual = try allocator.alloc(f64, dimension);
    defer allocator.free(opposite_residual);
    for (0..dimension) |column| {
        const epsilon = gasJacobianProbeG(current[column]);
        @memcpy(sampled_state, current);
        sampled_state[column] += epsilon;
        residualAt(allocator, scratch, base, sampled_state, inputs, options.transport_iteration_fraction, target, sampled_residual) catch return false;
        const central = current[column] >= epsilon;
        if (central) {
            @memcpy(opposite_state, current);
            opposite_state[column] -= epsilon;
            residualAt(allocator, scratch, base, opposite_state, inputs, options.transport_iteration_fraction, target, opposite_residual) catch return false;
        }
        for (0..dimension) |row| matrix[row * dimension + column] = if (central)
            (sampled_residual[row] - opposite_residual[row]) / (2 * epsilon)
        else
            (sampled_residual[row] - residual[row]) / epsilon;
    }
    // Solve in dimensionless coordinates. A profile can legitimately pair a
    // 1 g gaseous inventory with more than 1e8 g dissolved mass. Row scaling
    // alone protects the merit norm but leaves the direct Jacobian severely
    // column-ill-conditioned, allowing pivoting to report success with an
    // unusable trace-gas correction. Scale residual rows by their requested
    // tolerance and state columns by their inventory magnitude before the
    // direct solve, then convert the correction back to grams.
    for (0..dimension) |row| {
        const residual_scale =
            options.absolute_tolerance_g +
            options.relative_tolerance *
                @max(1.0, @abs(current[row]));
        rhs[row] = -residual[row] / residual_scale;
        for (0..dimension) |column| {
            const state_scale = @max(1.0, @abs(current[column]));
            matrix[row * dimension + column] *=
                state_scale / residual_scale;
        }
    }
    // Preserve a least-squares fallback before the direct solve mutates its
    // matrix. Donor and phase complementarity bounds can make the semismooth
    // dimensionless Jacobian rank deficient even though a descent direction
    // exists.
    @memset(normal_matrix, 0);
    @memset(normal_rhs, 0);
    for (0..dimension) |row| {
        for (0..dimension) |column| {
            const weighted_jacobian = matrix[row * dimension + column];
            normal_rhs[column] += weighted_jacobian * rhs[row];
            for (0..dimension) |other_column| {
                normal_matrix[column * dimension + other_column] +=
                    weighted_jacobian *
                    matrix[row * dimension + other_column];
            }
        }
    }
    var maximum_normal_diagonal: f64 = 0;
    for (0..dimension) |index|
        maximum_normal_diagonal = @max(
            maximum_normal_diagonal,
            @abs(normal_matrix[index * dimension + index]),
        );
    const damping =
        1.0e-10 * @max(1.0, maximum_normal_diagonal);
    for (0..dimension) |index|
        normal_matrix[index * dimension + index] += damping;

    if (!numerics.solveDenseLinearSystem(matrix, rhs, dimension)) {
        if (!numerics.solveDenseLinearSystem(
            normal_matrix,
            normal_rhs,
            dimension,
        )) return false;
        @memcpy(rhs, normal_rhs);
    }
    for (direction, rhs, current) |*component, scaled_component, value|
        component.* = scaled_component * @max(1.0, @abs(value));
    return true;
}

fn denseSpeciesNewtonDirection(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    inputs: Inputs,
    options: Options,
    target: []f64,
    sampled_residual: []f64,
    sampled_state: []f64,
    direction: []f64,
    species: usize,
) !bool {
    const inventory_count = scratch.gaseous_mass_g.len;
    const dimension = try std.math.mul(usize, scratch.cell_count, 3);
    // Dense pivoting is a precision-oriented tail solve, not the large-grid
    // path. Above this bound retain the O(n) matrix-free Krylov storage.
    if (dimension == 0 or dimension > 256) return false;
    const matrix = try allocator.alloc(f64, try std.math.mul(usize, dimension, dimension));
    defer allocator.free(matrix);
    const rhs = try allocator.alloc(f64, dimension);
    defer allocator.free(rhs);
    const opposite_state = try allocator.alloc(f64, current.len);
    defer allocator.free(opposite_state);
    const opposite_residual = try allocator.alloc(f64, current.len);
    defer allocator.free(opposite_residual);
    for (0..dimension) |row| rhs[row] = -residual[speciesVectorIndex(row, scratch.cell_count, inventory_count, species)];
    for (0..dimension) |column| {
        const state_index = speciesVectorIndex(column, scratch.cell_count, inventory_count, species);
        const epsilon = gasJacobianProbeG(current[state_index]);
        @memcpy(sampled_state, current);
        sampled_state[state_index] += epsilon;
        residualAt(allocator, scratch, base, sampled_state, inputs, options.transport_iteration_fraction, target, sampled_residual) catch return false;
        const central = current[state_index] >= epsilon;
        if (central) {
            @memcpy(opposite_state, current);
            opposite_state[state_index] -= epsilon;
            residualAt(allocator, scratch, base, opposite_state, inputs, options.transport_iteration_fraction, target, opposite_residual) catch return false;
        }
        for (0..dimension) |row| {
            const residual_index = speciesVectorIndex(row, scratch.cell_count, inventory_count, species);
            matrix[row * dimension + column] = if (central)
                (sampled_residual[residual_index] - opposite_residual[residual_index]) / (2 * epsilon)
            else
                (sampled_residual[residual_index] - residual[residual_index]) / epsilon;
        }
    }
    if (!numerics.solveDenseLinearSystem(matrix, rhs, dimension)) return false;
    for (rhs, 0..) |value, compact_index| direction[speciesVectorIndex(compact_index, scratch.cell_count, inventory_count, species)] = value;
    return true;
}

fn denseAllSpeciesNewtonDirection(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    inputs: Inputs,
    options: Options,
    target: []f64,
    sampled_residual: []f64,
    sampled_state: []f64,
    direction: []f64,
) !bool {
    const inventory_count = scratch.gaseous_mass_g.len;
    const dimension = try std.math.mul(usize, scratch.cell_count, 3);
    if (dimension == 0 or dimension > 256) return false;
    const matrix_elements = try std.math.mul(usize, dimension, dimension);
    const matrices = try allocator.alloc(f64, try std.math.mul(usize, gas.species_count, matrix_elements));
    defer allocator.free(matrices);
    const right_hand_sides = try allocator.alloc(f64, try std.math.mul(usize, gas.species_count, dimension));
    defer allocator.free(right_hand_sides);
    @memset(direction, 0);

    // Diffusion and phase exchange are species-local, but TRNSFR pressure
    // displacement depends on total gas pressure and therefore couples the
    // seven species. Simultaneously perturbing all species aliases those
    // cross derivatives into every diagonal block, most severely for the
    // small H2 inventory. Form each diagonal block with its own perturbation;
    // the subsequent line search/global residual retains the pressure
    // coupling while avoiding the aliased Jacobian.
    for (0..gas.species_count) |species| {
        for (0..dimension) |column| {
            @memcpy(sampled_state, current);
            const state_index = speciesVectorIndex(column, scratch.cell_count, inventory_count, species);
            const epsilon = gasJacobianProbeG(current[state_index]);
            sampled_state[state_index] += epsilon;
            residualAt(allocator, scratch, base, sampled_state, inputs, options.transport_iteration_fraction, target, sampled_residual) catch return false;
            for (0..dimension) |row| {
                const residual_index = speciesVectorIndex(row, scratch.cell_count, inventory_count, species);
                matrices[species * matrix_elements + row * dimension + column] =
                    (sampled_residual[residual_index] - residual[residual_index]) / epsilon;
            }
        }
    }
    var solved_newton_block = false;
    for (0..gas.species_count) |species| {
        const matrix = matrices[species * matrix_elements ..][0..matrix_elements];
        const rhs = right_hand_sides[species * dimension ..][0..dimension];
        for (0..dimension) |row| rhs[row] = -residual[speciesVectorIndex(row, scratch.cell_count, inventory_count, species)];
        if (numerics.solveDenseLinearSystem(matrix, rhs, dimension)) {
            solved_newton_block = true;
            for (rhs, 0..) |component, compact_index| direction[speciesVectorIndex(compact_index, scratch.cell_count, inventory_count, species)] = component;
        } else {
            // A bound-inactive phase can make one numerical block singular.
            // Do not let that block throttle the shared line search; preserve
            // every solvable Newton correction and leave this species for the
            // local Newton/Anderson/Picard path.
            for (0..dimension) |compact_index| {
                const index = speciesVectorIndex(compact_index, scratch.cell_count, inventory_count, species);
                direction[index] = 0;
            }
        }
    }
    return solved_newton_block;
}

fn filterIndependentSpeciesDirections(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    inputs: Inputs,
    options: Options,
    target: []f64,
    candidate_residual: []f64,
    candidate: []f64,
    direction: []f64,
) !bool {
    const inventory_count = scratch.gaseous_mass_g.len;
    var accepted_any = false;
    for (0..gas.species_count) |species| {
        const current_norm = try scaledSpeciesNorm(current, residual, options, inventory_count, species);
        var fraction: f64 = 1;
        var accepted = false;
        var search: u8 = 0;
        while (search < 20) : (search += 1) {
            @memcpy(candidate, current);
            for (candidate, current, direction, 0..) |*next, value, component, index| {
                if ((index % inventory_count) % gas.species_count == species) next.* = @max(0, value + fraction * component);
            }
            if (residualAt(allocator, scratch, base, candidate, inputs, options.transport_iteration_fraction, target, candidate_residual)) |_| {
                if (try scaledSpeciesNorm(candidate, candidate_residual, options, inventory_count, species) < current_norm) {
                    accepted = true;
                    accepted_any = true;
                    break;
                }
            } else |_| {}
            fraction *= 0.5;
        }
        for (direction, 0..) |*component, index| {
            if ((index % inventory_count) % gas.species_count == species) component.* *= if (accepted) fraction else 0;
        }
    }
    return accepted_any;
}

fn speciesVectorIndex(compact_index: usize, cell_count: usize, inventory_count: usize, species: usize) usize {
    const phase = compact_index / cell_count;
    const cell = compact_index % cell_count;
    return phase * inventory_count + cell * gas.species_count + species;
}

fn krylovNewtonDirection(
    allocator: std.mem.Allocator,
    scratch: *gas.State,
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    inputs: Inputs,
    options: Options,
    target: []f64,
    sampled_residual: []f64,
    direction: []f64,
    active_species: ?usize,
) !bool {
    const n = current.len;
    const restart = if (active_species != null)
        @min(n, try std.math.mul(usize, scratch.cell_count, 3))
    else
        @min(@as(usize, 12), n);
    if (restart == 0) return false;
    const basis = try allocator.alloc(f64, try std.math.mul(usize, restart + 1, n));
    defer allocator.free(basis);
    const hessenberg_rows = try std.math.add(usize, restart, 1);
    const hessenberg = try allocator.alloc(f64, try std.math.mul(usize, hessenberg_rows, restart));
    defer allocator.free(hessenberg);
    const cosine = try allocator.alloc(f64, restart);
    defer allocator.free(cosine);
    const sine = try allocator.alloc(f64, restart);
    defer allocator.free(sine);
    const projected_rhs = try allocator.alloc(f64, restart + 1);
    defer allocator.free(projected_rhs);
    const work = try allocator.alloc(f64, n);
    defer allocator.free(work);
    const sampled_state = try allocator.alloc(f64, n);
    defer allocator.free(sampled_state);
    @memset(hessenberg, 0);
    @memset(projected_rhs, 0);
    @memset(direction, 0);

    var beta_squared: f64 = 0;
    const inventory_count = scratch.gaseous_mass_g.len;
    for (residual, 0..) |value, index| {
        if (active_species == null or (index % inventory_count) % gas.species_count == active_species.?) beta_squared += value * value;
    }
    const beta = @sqrt(beta_squared);
    if (!std.math.isFinite(beta) or beta <= std.math.floatEps(f64)) return false;
    const first_basis = basis[0..n];
    for (residual, first_basis, 0..) |value, *entry, index| {
        entry.* = if (active_species == null or (index % inventory_count) % gas.species_count == active_species.?) -value / beta else 0;
    }
    projected_rhs[0] = beta;
    var perturbation_scale: f64 = 1;
    for (current, 0..) |value, index| {
        if (active_species == null or (index % inventory_count) % gas.species_count == active_species.?) {
            perturbation_scale = @max(perturbation_scale, @abs(value));
        }
    }

    var used: usize = 0;
    for (0..restart) |column| {
        const vector = basis[column * n ..][0..n];
        const epsilon = jacobianProbeMagnitude(perturbation_scale, beta);
        for (current, vector, sampled_state) |value, component, *sampled| sampled.* = @max(0, value + epsilon * component);
        residualAt(allocator, scratch, base, sampled_state, inputs, options.transport_iteration_fraction, target, sampled_residual) catch return false;
        for (sampled_residual, residual, work, 0..) |sampled, original, *entry, index| {
            entry.* = if (active_species == null or (index % inventory_count) % gas.species_count == active_species.?)
                (sampled - original) / epsilon
            else
                0;
        }

        for (0..column + 1) |row| {
            const row_basis = basis[row * n ..][0..n];
            var projection: f64 = 0;
            for (work, row_basis) |value, basis_value| projection += value * basis_value;
            hessenberg[row * restart + column] = projection;
            for (work, row_basis) |*value, basis_value| value.* -= projection * basis_value;
        }
        var work_norm_squared: f64 = 0;
        for (work) |value| work_norm_squared += value * value;
        const work_norm = @sqrt(work_norm_squared);
        hessenberg[(column + 1) * restart + column] = work_norm;
        if (work_norm > std.math.floatEps(f64) and column + 1 < restart + 1) {
            const next_basis = basis[(column + 1) * n ..][0..n];
            for (work, next_basis) |value, *entry| entry.* = value / work_norm;
        }

        for (0..column) |row| {
            const upper_index = row * restart + column;
            const lower_index = (row + 1) * restart + column;
            const upper = hessenberg[upper_index];
            const lower = hessenberg[lower_index];
            hessenberg[upper_index] = cosine[row] * upper + sine[row] * lower;
            hessenberg[lower_index] = -sine[row] * upper + cosine[row] * lower;
        }
        const diagonal_index = column * restart + column;
        const subdiagonal_index = (column + 1) * restart + column;
        const diagonal = hessenberg[diagonal_index];
        const subdiagonal = hessenberg[subdiagonal_index];
        const magnitude = @sqrt(diagonal * diagonal + subdiagonal * subdiagonal);
        if (!std.math.isFinite(magnitude) or magnitude <= std.math.floatEps(f64)) break;
        cosine[column] = diagonal / magnitude;
        sine[column] = subdiagonal / magnitude;
        hessenberg[diagonal_index] = magnitude;
        hessenberg[subdiagonal_index] = 0;
        const rhs = projected_rhs[column];
        projected_rhs[column] = cosine[column] * rhs;
        projected_rhs[column + 1] = -sine[column] * rhs;
        used = column + 1;
        if (@abs(projected_rhs[column + 1]) <= 0.01 * beta) break;
    }
    if (used == 0) return false;
    var row = used;
    while (row > 0) {
        row -= 1;
        var rhs = projected_rhs[row];
        for (row + 1..used) |column| rhs -= hessenberg[row * restart + column] * projected_rhs[column];
        const diagonal = hessenberg[row * restart + row];
        if (!std.math.isFinite(diagonal) or @abs(diagonal) <= std.math.floatEps(f64)) return false;
        projected_rhs[row] = rhs / diagonal;
    }
    for (0..used) |column| {
        const vector = basis[column * n ..][0..n];
        for (direction, vector) |*value, component| value.* += projected_rhs[column] * component;
    }
    for (direction) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn validate(state: *const gas.State, inputs: Inputs, options: Options) !void {
    try state.validateFinite();
    const n = std.math.mul(usize, state.cell_count, gas.species_count) catch
        return error.CoupledGasStateSizeMismatch;
    if (inputs.face_conductance_m3_per_step.len != try std.math.mul(usize, inputs.faces.len, gas.species_count)) return error.GasFaceParameterSizeMismatch;
    if (inputs.water_volume_m3.len != state.cell_count or inputs.band_water_volume_m3.len != state.cell_count or inputs.bubbling_enabled.len != state.cell_count or inputs.mass_solubility_ratio.len != n or inputs.gas_water_exchange_rate_per_step.len != n or inputs.band_gas_water_exchange_rate_per_step.len != n) return error.CoupledGasInputSizeMismatch;
    if (inputs.bubble_receiver_cell_by_cell) |receivers| {
        if (receivers.len != state.cell_count) return error.CoupledGasInputSizeMismatch;
        for (receivers) |receiver| if (receiver) |cell|
            if (cell >= state.cell_count) return error.GasBubbleReceiverOutOfBounds;
    }
    if (inputs.atmospheric_flux_g_by_component) |fluxes| if (fluxes.len != n) return error.GasBoundaryFluxSizeMismatch;
    if (inputs.subsurface_flux_g_by_component) |fluxes| if (fluxes.len != n) return error.GasBoundaryFluxSizeMismatch;
    if (inputs.face_flux_g_by_component) |fluxes| if (fluxes.len !=
        try std.math.mul(usize, inputs.faces.len, gas.species_count))
        return error.GasFaceFluxSizeMismatch;
    if (!std.math.isFinite(options.absolute_tolerance_g) or options.absolute_tolerance_g <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or !std.math.isFinite(options.transport_iteration_fraction) or options.transport_iteration_fraction <= 0 or options.transport_iteration_fraction > 1 or options.max_iterations == 0) return error.InvalidCoupledGasSolverOptions;
    for (inputs.water_volume_m3, inputs.band_water_volume_m3) |water, band_water| if (!std.math.isFinite(water) or water < 0 or !std.math.isFinite(band_water) or band_water < 0) return error.InvalidCoupledGasInput;
    for (inputs.mass_solubility_ratio) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidCoupledGasInput;
    for (inputs.gas_water_exchange_rate_per_step, inputs.band_gas_water_exchange_rate_per_step) |rate, band_rate| if (!std.math.isFinite(rate) or rate < 0 or !std.math.isFinite(band_rate) or band_rate < 0) return error.InvalidCoupledGasInput;
    for (inputs.faces) |face| {
        if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell)
            return error.InvalidGasTransportFace;
    }
    for (inputs.face_conductance_m3_per_step) |conductance| {
        if (!std.math.isFinite(conductance) or conductance < 0) return error.InvalidCoupledGasInput;
    }
    try validateBoundaries(state.cell_count, inputs.atmospheric_boundaries);
    try validateBoundaries(state.cell_count, inputs.subsurface_boundaries);
}

fn validateBoundaries(cell_count: usize, boundaries: []const atmosphere.Boundary) !void {
    for (boundaries) |boundary| {
        if (boundary.cell_index >= cell_count or
            !std.math.isFinite(boundary.aerodynamic_conductance_m3_per_step) or
            boundary.aerodynamic_conductance_m3_per_step < 0 or
            !std.math.isFinite(boundary.pressure_exchange_fraction) or
            boundary.pressure_exchange_fraction < 0 or
            boundary.pressure_exchange_fraction > 1)
        {
            return error.InvalidGasBoundary;
        }
        for (boundary.interior_conductance_m3_per_step, boundary.atmospheric_concentration_g_per_m3) |conductance, concentration| {
            if (!std.math.isFinite(conductance) or conductance < 0 or
                !std.math.isFinite(concentration) or concentration < 0)
            {
                return error.InvalidGasBoundary;
            }
        }
    }
}

fn validationTestInputs() Inputs {
    return .{
        .faces = &.{},
        .face_conductance_m3_per_step = &.{},
        .atmospheric_boundaries = &.{},
        .water_volume_m3 = &.{0},
        .band_water_volume_m3 = &.{0},
        .mass_solubility_ratio = &([_]f64{1} ** gas.species_count),
        .gas_water_exchange_rate_per_step = &([_]f64{0} ** gas.species_count),
        .band_gas_water_exchange_rate_per_step = &([_]f64{0} ** gas.species_count),
        .bubbling_enabled = &.{false},
    };
}

test "coupled gas preflight rejects malformed and nonphysical state" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.temperature_k[0] = 293.15;
    try validate(&state, validationTestInputs(), .{ .max_iterations = 1 });

    const complete_masses = state.gaseous_mass_g;
    state.gaseous_mass_g = complete_masses[0 .. complete_masses.len - 1];
    try std.testing.expectError(error.CoupledGasStateSizeMismatch, validate(&state, validationTestInputs(), .{ .max_iterations = 1 }));
    state.gaseous_mass_g = complete_masses;

    state.gaseous_mass_g[0] = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteGasTransportState, validate(&state, validationTestInputs(), .{ .max_iterations = 1 }));
    state.gaseous_mass_g[0] = -1;
    try std.testing.expectError(error.NegativeGasTransportState, validate(&state, validationTestInputs(), .{ .max_iterations = 1 }));
    state.gaseous_mass_g[0] = 0;
    state.temperature_k[0] = 0;
    try std.testing.expectError(error.InvalidGasTransportTemperature, validate(&state, validationTestInputs(), .{ .max_iterations = 1 }));
}

test "coupled gas preflight validates topology conductance and boundaries" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.temperature_k[0] = 293.15;
    const self_face = [_]gas.Face{.{ .first_cell = 0, .second_cell = 0 }};
    const conductance = [_]f64{0} ** gas.species_count;
    var inputs = validationTestInputs();
    inputs.faces = &self_face;
    inputs.face_conductance_m3_per_step = &conductance;
    try std.testing.expectError(error.InvalidGasTransportFace, validate(&state, inputs, .{ .max_iterations = 1 }));

    const invalid_conductance = [_]f64{-1} ** gas.species_count;
    inputs.faces = &.{};
    inputs.face_conductance_m3_per_step = &.{};
    inputs.atmospheric_boundaries = &.{.{
        .cell_index = 0,
        .aerodynamic_conductance_m3_per_step = 0,
        .interior_conductance_m3_per_step = invalid_conductance,
        .atmospheric_concentration_g_per_m3 = conductance,
    }};
    try std.testing.expectError(error.InvalidGasBoundary, validate(&state, inputs, .{ .max_iterations = 1 }));
}

fn copyStateToVector(state: *const gas.State, vector: []f64) void {
    const n = state.gaseous_mass_g.len;
    @memcpy(vector[0..n], state.gaseous_mass_g);
    @memcpy(vector[n .. 2 * n], state.dissolved_mass_g);
    @memcpy(vector[2 * n .. 3 * n], state.band_dissolved_mass_g);
}

fn copyVectorToState(vector: []const f64, state: *gas.State) void {
    const n = state.gaseous_mass_g.len;
    @memcpy(state.gaseous_mass_g, vector[0..n]);
    @memcpy(state.dissolved_mass_g, vector[n .. 2 * n]);
    @memcpy(state.band_dissolved_mass_g, vector[2 * n .. 3 * n]);
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidCoupledGasCandidate;
        candidate.* = @max(0, candidate.*);
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteCoupledGasState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn newtonMeritImproves(candidate: []const f64, candidate_residual: []const f64, options: Options, current_maximum_scaled_residual: f64) !bool {
    if (!std.math.isFinite(current_maximum_scaled_residual) or current_maximum_scaled_residual < 0)
        return error.NonFiniteCoupledGasState;
    return try scaledNorm(candidate, candidate_residual, options) < current_maximum_scaled_residual;
}

fn newtonAcceptanceTarget(current_scaled_residual: f64) f64 {
    return current_scaled_residual;
}

fn speciesBlockMeritImproves(current: []const f64, current_residual: []const f64, candidate: []const f64, candidate_residual: []const f64, options: Options, inventory_count: usize, species: usize, current_maximum_scaled_residual: f64) !bool {
    const current_species = try scaledSpeciesNorm(current, current_residual, options, inventory_count, species);
    const candidate_species = try scaledSpeciesNorm(candidate, candidate_residual, options, inventory_count, species);
    const candidate_global = try scaledNorm(candidate, candidate_residual, options);
    return candidate_species < current_species and candidate_global <= current_maximum_scaled_residual;
}

fn scaledSpeciesNorm(state: []const f64, residual: []const f64, options: Options, inventory_count: usize, active_species: usize) !f64 {
    var maximum: f64 = 0;
    for (state, residual, 0..) |value, difference, index| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteCoupledGasState;
        if ((index % inventory_count) % gas.species_count != active_species) continue;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn worstResidualSpecies(state: []const f64, residual: []const f64, options: Options, inventory_count: usize) !usize {
    var worst_species: usize = 0;
    var worst_norm: f64 = -1;
    for (0..gas.species_count) |species| {
        const norm = try scaledSpeciesNorm(state, residual, options, inventory_count, species);
        if (norm > worst_norm) {
            worst_norm = norm;
            worst_species = species;
        }
    }
    return worst_species;
}

fn unconvergedSpeciesCount(
    state: []const f64,
    residual: []const f64,
    options: Options,
    inventory_count: usize,
) !usize {
    var count: usize = 0;
    for (0..gas.species_count) |species| {
        if (try scaledSpeciesNorm(
            state,
            residual,
            options,
            inventory_count,
            species,
        ) > 1) count += 1;
    }
    return count;
}

fn worstResidualIndex(state: []const f64, residual: []const f64, options: Options) !usize {
    if (state.len == 0 or state.len != residual.len) return error.NonFiniteCoupledGasState;
    var worst_index: usize = 0;
    var worst_norm: f64 = -1;
    for (state, residual, 0..) |value, difference, index| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteCoupledGasState;
        const norm = try scaledCoordinateResidual(value, difference, options);
        if (norm > worst_norm) {
            worst_norm = norm;
            worst_index = index;
        }
    }
    return worst_index;
}

fn scaledCoordinateResidual(value: f64, difference: f64, options: Options) !f64 {
    if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteCoupledGasState;
    return @abs(difference) / (options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, @abs(value)));
}

fn isScaleDepletedGas(base_gas_mass_g: f64, base_dissolved_mass_g: f64, options: Options) bool {
    // Active-set detection must occur before the gas coordinate reaches the
    // final nonlinear tolerance. Using `relative_tolerance * dissolved`
    // classified a gas pool at 6.8e-6 of its paired aqueous inventory as an
    // ordinary interior coordinate, bypassing the conservative two-phase
    // Newton block and leaving the fixed-point map on a donor-clamped cycle.
    // The square-root tolerance is the standard scale-separation threshold
    // for selecting the active set; it does not relax the root acceptance
    // criterion, which remains `scaledNorm <= 1`.
    const active_set_relative_scale =
        @sqrt(options.relative_tolerance);
    return base_gas_mass_g > 0 and
        base_dissolved_mass_g > 0 and
        base_gas_mass_g <=
            options.absolute_tolerance_g +
                active_set_relative_scale * base_dissolved_mass_g;
}

fn isNumericallyAtNonnegativeBound(value: f64, assembled_target: f64) bool {
    return value <= 64.0 * std.math.floatEps(f64) *
        @max(1.0, @abs(assembled_target));
}

fn hasUnresolvedPositiveBound(
    value: f64,
    assembled_target: f64,
    residual: f64,
) bool {
    return isNumericallyAtNonnegativeBound(value, assembled_target) and
        residual > 0;
}

fn requiresNonlocalConservativeBracket(
    receiver_mass_g: f64,
    donor_mass_g: f64,
    options: Options,
) bool {
    return donor_mass_g >
        (1.0 / @sqrt(options.relative_tolerance)) *
            receiver_mass_g;
}

fn hasExtremeIncomingDonorResidual(
    current: []const f64,
    residual: []const f64,
    target: []const f64,
    inputs: Inputs,
    options: Options,
    inventory_count: usize,
) bool {
    return (worstExtremeIncomingGasIndex(
        current,
        residual,
        target,
        inputs,
        options,
        inventory_count,
    ) catch return false) != null;
}

/// Select the largest unresolved gaseous residual that is trapped behind an
/// incoming, scale-separated donor face. The global worst coordinate may be
/// the paired dissolved phase; using it to select a face caused the
/// conservative safeguard to be skipped and spent the remaining NPH*NPG
/// budget on slow projected phase corrections.
fn worstExtremeIncomingGasIndex(
    current: []const f64,
    residual: []const f64,
    target: []const f64,
    inputs: Inputs,
    options: Options,
    inventory_count: usize,
) !?usize {
    return selectScaleSeparatedIncomingGasIndex(
        current,
        residual,
        target,
        inputs,
        options,
        inventory_count,
        true,
    );
}

/// Multi-species face Newton does not require the scalar bracket's stronger
/// residual-greater-than-inventory condition. Its dense pressure block is
/// valid for any tolerance-unconverged positive incoming coordinate.
fn worstScaleSeparatedIncomingGasIndex(
    current: []const f64,
    residual: []const f64,
    target: []const f64,
    inputs: Inputs,
    options: Options,
    inventory_count: usize,
) !?usize {
    return selectScaleSeparatedIncomingGasIndex(
        current,
        residual,
        target,
        inputs,
        options,
        inventory_count,
        false,
    );
}

fn selectScaleSeparatedIncomingGasIndex(
    current: []const f64,
    residual: []const f64,
    target: []const f64,
    inputs: Inputs,
    options: Options,
    inventory_count: usize,
    require_residual_larger_than_inventory: bool,
) !?usize {
    var worst_index: ?usize = null;
    var worst_norm: f64 = -1;
    for (inputs.faces) |face| {
        for (0..gas.species_count) |species| {
            const first =
                face.first_cell * gas.species_count + species;
            const second =
                face.second_cell * gas.species_count + species;
            if (first >= inventory_count or
                second >= inventory_count) return error.GasFaceCellOutOfRange;
            if (current[first] > 0 and residual[first] > 0 and
                (!require_residual_larger_than_inventory or
                    residual[first] > current[first]) and
                target[first] > current[first] and
                requiresNonlocalConservativeBracket(
                    current[first],
                    current[second],
                    options,
                ))
            {
                const norm = try scaledCoordinateResidual(
                    current[first],
                    residual[first],
                    options,
                );
                if (norm > worst_norm) {
                    worst_norm = norm;
                    worst_index = first;
                }
            }
            if (current[second] > 0 and residual[second] > 0 and
                (!require_residual_larger_than_inventory or
                    residual[second] > current[second]) and
                target[second] > current[second] and
                requiresNonlocalConservativeBracket(
                    current[second],
                    current[first],
                    options,
                ))
            {
                const norm = try scaledCoordinateResidual(
                    current[second],
                    residual[second],
                    options,
                );
                if (norm > worst_norm) {
                    worst_norm = norm;
                    worst_index = second;
                }
            }
        }
    }
    return worst_index;
}

fn vectorsEqual(a: []const f64, b: []const f64) bool {
    if (a.len != b.len) return false;
    for (a, b) |left, right| if (left != right) return false;
    return true;
}

fn jacobianProbeMagnitude(state_scale: f64, residual_norm: f64) f64 {
    _ = residual_norm;
    return @sqrt(std.math.floatEps(f64)) * @max(1.0, state_scale);
}

test "Newton merit cannot trade a worse limiting species for a lower L2 residual" {
    const options: Options = .{ .absolute_tolerance_g = 1, .relative_tolerance = 1e-12, .max_iterations = 2 };
    const state = [_]f64{ 0, 0 };
    const current_residual = [_]f64{ 10, 10 };
    const candidate_residual = [_]f64{ 11, 0 };
    const current_maximum = try scaledNorm(&state, &current_residual, options);
    // The candidate's L2 norm is smaller (11 versus sqrt(200)), but its
    // limiting component is worse and therefore must not be accepted.
    try std.testing.expect(!(try newtonMeritImproves(&state, &candidate_residual, options, current_maximum)));
}

test "worst residual coordinate uses the solver infinity norm scaling" {
    const options: Options = .{ .absolute_tolerance_g = 1.0e-12, .relative_tolerance = 1.0e-8, .max_iterations = 2 };
    try std.testing.expectEqual(@as(usize, 1), try worstResidualIndex(&.{ 1.0e6, 0.5, 2.0 }, &.{ 1.0e-3, 2.0e-8, 1.0e-8 }, options));
}

test "Newton globalization accepts every strict infinity-norm improvement" {
    try std.testing.expectApproxEqAbs(@as(f64, 10), newtonAcceptanceTarget(10), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), newtonAcceptanceTarget(2), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.337), newtonAcceptanceTarget(1.337), 1e-15);
}

test "coordinate Newton merit can advance one of two tied limiting residuals" {
    const options: Options = .{ .absolute_tolerance_g = 1, .relative_tolerance = 1.0e-12, .max_iterations = 2 };
    const current = [_]f64{ 0, 0 };
    const residual = [_]f64{ 2, 2 };
    const candidate = [_]f64{ 0, 0 };
    const candidate_residual = [_]f64{ 0.5, 2 };
    const current_norm = try scaledNorm(&current, &residual, options);
    try std.testing.expect(try scaledCoordinateResidual(candidate[0], candidate_residual[0], options) < try scaledCoordinateResidual(current[0], residual[0], options));
    try std.testing.expect(try scaledNorm(&candidate, &candidate_residual, options) <= current_norm);
}

test "trace-gas correction is not hidden by a large inventory" {
    const current = [_]f64{ 1.0e12, 1.0e-9 };
    const corrected = [_]f64{ 1.0e12, 1.0e-9 + 1.0e-15 };
    try std.testing.expect(!vectorsEqual(&current, &corrected));
    try std.testing.expect(vectorsEqual(&current, &current));
}

test "projected zero cannot be accepted while its conservative source is positive" {
    try std.testing.expect(hasUnresolvedPositiveBound(0, 0.016, 0.016));
    try std.testing.expect(!hasUnresolvedPositiveBound(0, 0, 0));
    try std.testing.expect(!hasUnresolvedPositiveBound(0.1, 0.2, 0.1));
    try std.testing.expect(!hasUnresolvedPositiveBound(0, 0, -0.1));
}

test "extreme donor scale selects nonlocal conservative face bracket" {
    const options: Options = .{
        .relative_tolerance = 1.0e-8,
        .max_iterations = 80,
    };
    try std.testing.expect(
        requiresNonlocalConservativeBracket(
            1.0e-4,
            1.0e8,
            options,
        ),
    );
    try std.testing.expect(
        !requiresNonlocalConservativeBracket(
            1,
            500,
            options,
        ),
    );
}

test "semismooth budget handoff requires an active incoming extreme donor" {
    const options: Options = .{
        .relative_tolerance = 1.0e-8,
        .max_iterations = 80,
    };
    const inventory_count = 2 * gas.species_count;
    var current = [_]f64{1} ** (3 * 2 * gas.species_count);
    var residual = [_]f64{0} ** (3 * 2 * gas.species_count);
    var target = current;
    current[0] = 1.0e-4;
    current[gas.species_count] = 1.0e8;
    residual[0] = 1;
    target[0] = 1.0001;
    const inputs: Inputs = .{
        .faces = &.{.{
            .first_cell = 0,
            .second_cell = 1,
        }},
        .face_conductance_m3_per_step = &([_]f64{0} ** gas.species_count),
        .atmospheric_boundaries = &.{},
        .water_volume_m3 = &.{ 0, 0 },
        .band_water_volume_m3 = &.{ 0, 0 },
        .mass_solubility_ratio = &([_]f64{1} ** (2 * gas.species_count)),
        .gas_water_exchange_rate_per_step = &([_]f64{0} ** (2 * gas.species_count)),
        .band_gas_water_exchange_rate_per_step = &([_]f64{0} ** (2 * gas.species_count)),
        .bubbling_enabled = &.{ false, false },
    };
    try std.testing.expect(hasExtremeIncomingDonorResidual(
        &current,
        &residual,
        &target,
        inputs,
        options,
        inventory_count,
    ));
    residual[0] = 0;
    target[0] = current[0];
    try std.testing.expect(!hasExtremeIncomingDonorResidual(
        &current,
        &residual,
        &target,
        inputs,
        options,
        inventory_count,
    ));
}

test "extreme gas face remains selectable when a dissolved residual is larger" {
    const options: Options = .{
        .relative_tolerance = 1.0e-8,
        .max_iterations = 80,
    };
    const inventory_count = 2 * gas.species_count;
    var current = [_]f64{1} ** (3 * 2 * gas.species_count);
    var residual = [_]f64{0} ** (3 * 2 * gas.species_count);
    var target = current;
    current[0] = 1.0e-4;
    current[gas.species_count] = 1.0e8;
    residual[0] = 1;
    target[0] = 1.0001;
    residual[inventory_count + 1] = 1.0e6;
    const inputs: Inputs = .{
        .faces = &.{.{
            .first_cell = 0,
            .second_cell = 1,
        }},
        .face_conductance_m3_per_step = &([_]f64{0} ** gas.species_count),
        .atmospheric_boundaries = &.{},
        .water_volume_m3 = &.{ 0, 0 },
        .band_water_volume_m3 = &.{ 0, 0 },
        .mass_solubility_ratio = &([_]f64{1} ** (2 * gas.species_count)),
        .gas_water_exchange_rate_per_step = &([_]f64{0} ** (2 * gas.species_count)),
        .band_gas_water_exchange_rate_per_step = &([_]f64{0} ** (2 * gas.species_count)),
        .bubbling_enabled = &.{ false, false },
    };
    try std.testing.expectEqual(
        @as(?usize, 0),
        try worstExtremeIncomingGasIndex(
            &current,
            &residual,
            &target,
            inputs,
            options,
            inventory_count,
        ),
    );
    residual[0] = 1.0e-8;
    target[0] = current[0] + residual[0];
    try std.testing.expectEqual(
        @as(?usize, null),
        try worstExtremeIncomingGasIndex(
            &current,
            &residual,
            &target,
            inputs,
            options,
            inventory_count,
        ),
    );
    try std.testing.expectEqual(
        @as(?usize, 0),
        try worstScaleSeparatedIncomingGasIndex(
            &current,
            &residual,
            &target,
            inputs,
            options,
            inventory_count,
        ),
    );
}

test "single remaining gas species is detected across all phases and cells" {
    const options: Options = .{
        .relative_tolerance = 1.0e-8,
        .max_iterations = 80,
    };
    const inventory_count = 2 * gas.species_count;
    const current = [_]f64{1} ** (3 * 2 * gas.species_count);
    var residual = [_]f64{0} ** (3 * 2 * gas.species_count);
    residual[gas.species_count + 6] = 1.0e-4;
    residual[2 * inventory_count + gas.species_count + 6] = -1.0e-4;
    try std.testing.expectEqual(
        @as(usize, 1),
        try unconvergedSpeciesCount(
            &current,
            &residual,
            options,
            inventory_count,
        ),
    );
    residual[2] = 1.0e-4;
    try std.testing.expectEqual(
        @as(usize, 2),
        try unconvergedSpeciesCount(
            &current,
            &residual,
            options,
            inventory_count,
        ),
    );
}

test "depleted phase classification is relative to its aqueous inventory" {
    const options: Options = .{
        .absolute_tolerance_g = 1.0e-14,
        .relative_tolerance = 1.0e-8,
        .max_iterations = 80,
    };
    try std.testing.expect(isScaleDepletedGas(2.3e-13, 5.42, options));
    try std.testing.expect(!isScaleDepletedGas(2.4e-5, 2.1e-5, options));
    try std.testing.expect(!isScaleDepletedGas(0, 5.42, options));
    try std.testing.expect(!isScaleDepletedGas(1.0e-13, 0, options));
}

test "matrix-free Jacobian probe retains derivative precision at gram-scale state" {
    const probe = jacobianProbeMagnitude(1, 1);
    try std.testing.expect(probe > 1.0e-9);
    try std.testing.expect(probe < 1.0e-7);
    try std.testing.expectEqual(probe * 1.0e6, jacobianProbeMagnitude(1.0e6, 1.0e6));
    try std.testing.expectEqual(probe, jacobianProbeMagnitude(1, 1.0e-8));
}

test "tail refinement selects the species controlling convergence" {
    const inventory_count = gas.species_count;
    const state = [_]f64{0} ** gas.species_count;
    var residual = [_]f64{0} ** gas.species_count;
    residual[@intFromEnum(gas.Species.nitrous_oxide)] = 3;
    residual[@intFromEnum(gas.Species.carbon_dioxide)] = 2;
    const options: Options = .{ .absolute_tolerance_g = 1, .relative_tolerance = 0, .max_iterations = 8 };
    try std.testing.expectEqual(@intFromEnum(gas.Species.nitrous_oxide), try worstResidualSpecies(&state, &residual, options, inventory_count));
}

test "species tail accepts progress when another species ties the global maximum" {
    const inventory_count = gas.species_count;
    const state = [_]f64{0} ** gas.species_count;
    var current_residual = [_]f64{0} ** gas.species_count;
    var candidate_residual = [_]f64{0} ** gas.species_count;
    const selected = @intFromEnum(gas.Species.oxygen);
    const tied = @intFromEnum(gas.Species.nitrous_oxide);
    current_residual[selected] = 5;
    current_residual[tied] = 5;
    candidate_residual[selected] = 4;
    candidate_residual[tied] = 5;
    const options: Options = .{ .absolute_tolerance_g = 1, .relative_tolerance = 0, .max_iterations = 8 };
    const global = try scaledNorm(&state, &current_residual, options);
    try std.testing.expect(try speciesBlockMeritImproves(&state, &current_residual, &state, &candidate_residual, options, inventory_count, selected, global));
}

test {
    // Keeps the extracted tests discoverable by `zig build test`,
    // which only reaches files reachable by import.
    _ = @import("../../validation/coupled_gas_solver_test.zig");
}
