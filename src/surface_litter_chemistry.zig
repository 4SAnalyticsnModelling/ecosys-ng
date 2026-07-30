const std = @import("std");
const ledger = @import("surface_litter_reaction_ledger.zig");
const activity_coefficients = @import("solute_activity_coefficients.zig");
const numerics = @import("numerics.zig");

/// All litter chemistry storage is allocated by runtime horizontal cell count.
/// Field units are inherited from the explicit reaction ledger: aqueous and
/// solid values are mol/m3; exchange values are mol/Mg.
pub const Cell = ledger.Transformations;

/// HOUR1 charge classes for surface litter. Concentration state is converted
/// back to extensive mol before the shared Debye-Huckel calculation.
pub fn activityCoefficients(cell: Cell, litter_water_volume_m3: f64) !activity_coefficients.Result {
    if (!std.math.isFinite(litter_water_volume_m3) or litter_water_volume_m3 <= 0) return error.InvalidLitterWaterVolume;
    try validateCell(cell);
    const water = litter_water_volume_m3;
    return activity_coefficients.calculate(.{
        .trivalent_cations_mol = (cell.aluminum_mol_per_m3 + cell.iron_mol_per_m3) * water,
        .trivalent_anions_mol = 0,
        .divalent_cations_mol = (cell.calcium_mol_per_m3 + cell.magnesium_mol_per_m3) * water,
        .divalent_anions_mol = (cell.sulfate_mol_per_m3 + cell.carbonate_mol_per_m3 + cell.hpo4_mol_p_per_m3) * water,
        .monovalent_cations_mol = (cell.ammonium_mol_per_m3 + cell.hydrogen_mol_per_m3 + cell.sodium_mol_per_m3 + cell.potassium_mol_per_m3) * water,
        .monovalent_anions_mol = (cell.hydroxide_mol_per_m3 + cell.nitrate_mol_per_m3 + cell.chloride_mol_per_m3 + cell.bicarbonate_mol_per_m3 + cell.h2po4_mol_p_per_m3) * water,
        .neutral_solutes_mol = (cell.ammonia_mol_per_m3 + cell.carbon_dioxide_mol_per_m3) * water,
    }, water);
}

pub const State = struct {
    allocator: std.mem.Allocator,
    cells: []Cell,
    /// Water volume represented by the stored solid-mineral concentrations.
    /// A dry cell retains its last positive reference so extensive solids
    /// remain representable without infinity.
    mineral_reference_water_m3: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroLitterChemistryCellCount;
        const cells = try allocator.alloc(Cell, cell_count);
        errdefer allocator.free(cells);
        const mineral_reference_water_m3 = try allocator.alloc(f64, cell_count);
        for (cells) |*cell| zeroValue(Cell, cell);
        @memset(mineral_reference_water_m3, 0);
        return .{
            .allocator = allocator,
            .cells = cells,
            .mineral_reference_water_m3 = mineral_reference_water_m3,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.mineral_reference_water_m3);
        self.allocator.free(self.cells);
        self.* = undefined;
    }

    pub fn bindMineralReferenceWater(self: *State, water_m3: []const f64) !void {
        if (water_m3.len != self.cells.len)
            return error.LitterMineralReferenceDimensionMismatch;
        for (water_m3) |water|
            if (!std.math.isFinite(water) or water < 0)
                return error.InvalidLitterMineralReferenceWater;
        @memcpy(self.mineral_reference_water_m3, water_m3);
    }

    /// Preserves each solid mineral's extensive mol inventory while the
    /// litter-water carrier changes.
    pub fn renormalizeMinerals(self: *State, cell_index: usize, water_m3: f64) !void {
        if (cell_index >= self.cells.len)
            return error.LitterChemistryCellIndexOutOfBounds;
        if (!std.math.isFinite(water_m3) or water_m3 < 0)
            return error.InvalidLitterMineralReferenceWater;
        const old_water_m3 = self.mineral_reference_water_m3[cell_index];
        if (!std.math.isFinite(old_water_m3) or old_water_m3 < 0)
            return error.InvalidLitterMineralReferenceWater;
        var cell = self.cells[cell_index];
        try validateMinerals(cell);
        if (old_water_m3 == 0) {
            if (hasMineralInventory(cell))
                return error.UnboundLitterMineralInventory;
            self.mineral_reference_water_m3[cell_index] = water_m3;
            return;
        }
        // A concentration is undefined when dry. Retain the last positive
        // concentration/reference pair as an exact extensive inventory.
        if (water_m3 == 0) return;
        const scale = old_water_m3 / water_m3;
        inline for (@typeInfo(@TypeOf(cell.phosphate_minerals)).@"struct".fields) |field|
            @field(cell.phosphate_minerals, field.name) *= scale;
        inline for (@typeInfo(@TypeOf(cell.salt_minerals)).@"struct".fields) |field|
            @field(cell.salt_minerals, field.name) *= scale;
        try validateMinerals(cell);
        self.cells[cell_index] = cell;
        self.mineral_reference_water_m3[cell_index] = water_m3;
    }
};

fn hasMineralInventory(cell: Cell) bool {
    inline for (@typeInfo(@TypeOf(cell.phosphate_minerals)).@"struct".fields) |field|
        if (@field(cell.phosphate_minerals, field.name) != 0) return true;
    inline for (@typeInfo(@TypeOf(cell.salt_minerals)).@"struct".fields) |field|
        if (@field(cell.salt_minerals, field.name) != 0) return true;
    return false;
}

fn validateMinerals(cell: Cell) !void {
    inline for (@typeInfo(@TypeOf(cell.phosphate_minerals)).@"struct".fields) |field| {
        const value = @field(cell.phosphate_minerals, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterMineralInventory;
    }
    inline for (@typeInfo(@TypeOf(cell.salt_minerals)).@"struct".fields) |field| {
        const value = @field(cell.salt_minerals, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterMineralInventory;
    }
}

test "solid mineral inventory survives wet dry and rewet carrier changes" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].salt_minerals.calcite_mol_per_m3 = 2;
    state.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3 = 3;
    try state.bindMineralReferenceWater(&.{1});

    try state.renormalizeMinerals(0, 0);
    try std.testing.expectEqual(@as(f64, 2), state.cells[0].salt_minerals.calcite_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), state.mineral_reference_water_m3[0]);

    try state.renormalizeMinerals(0, 0.25);
    try std.testing.expectEqual(@as(f64, 8), state.cells[0].salt_minerals.calcite_mol_per_m3);
    try std.testing.expectEqual(
        @as(f64, 12),
        state.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        state.cells[0].salt_minerals.calcite_mol_per_m3 * 0.25,
    );
}

pub const Environment = struct {
    litter_mass_per_water_volume_Mg_per_m3: f64,
    dynamic_salts: bool,
};

pub const Options = struct {
    absolute_tolerance: f64 = 1e-11,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
    /// Reaction-equilibrium ceiling `MRXN=60`; successful solves exit early.
    max_iterations: u16 = 60,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Supplies litter-specific reaction extents.  Keeping this callback free of
/// allocation lets the same solver drive CPU grid kernels and future GPU rate
/// evaluators without embedding storage policy in the science equations.
pub const Evaluator = struct {
    context: *const anyopaque,
    evaluate: *const fn (context: *const anyopaque, cell: Cell) anyerror!ledger.ReactionExtents,
    equilibrate_cation_exchange: ?*const fn (context: *const anyopaque, cell: Cell) anyerror!Cell = null,
    phosphate_mineral_equilibrium_residuals: ?*const fn (
        context: *const anyopaque,
        cell: Cell,
    ) anyerror!ledger.PhosphateMineralExtents = null,
};

/// Conservative directional Newton–Raphson with relaxed Picard fallback.
/// State is copied into local iterates and committed only after convergence.
pub fn solveCell(state: *State, cell_index: usize, environment: Environment, evaluator: Evaluator, options: Options) !Result {
    if (cell_index >= state.cells.len) return error.LitterChemistryCellIndexOutOfBounds;
    try validateOptions(options);
    if (!std.math.isFinite(environment.litter_mass_per_water_volume_Mg_per_m3) or environment.litter_mass_per_water_volume_Mg_per_m3 <= 0) return error.InvalidLitterMassWaterRatio;
    try validateCell(state.cells[cell_index]);

    var current = state.cells[cell_index];
    // SOLUTE's ISALTG=0 branch is outside `DO M=1,MRXN`: its TPDH,
    // TADAH, TADCH, TSLH, and TRWH coefficients are already hourly kinetic
    // increments scaled by XNFH. Applying them repeatedly until every rate
    // vanishes would multiply one hour of precipitation/association by the
    // equilibrium iteration ceiling. Commit exactly one conservative,
    // substrate-bounded hourly transformation. The iterative hybrid below is
    // reserved for the dynamic-salt MRXN equilibrium branch.
    if (!environment.dynamic_salts) {
        const changes = try changesAt(current, environment, evaluator);
        const maximum_fraction =
            maximumAdmissibleFraction(Cell, current, changes);
        const accepted_fraction =
            if (std.math.isFinite(maximum_fraction))
                @min(1.0, maximum_fraction)
            else
                1.0;
        if (accepted_fraction <= 0)
            return error.NoPhysicallyAdmissibleFixedPhKineticStep;
        const accepted = try applyAdmissibleFraction(
            current,
            changes,
            accepted_fraction,
        );
        state.cells[cell_index] = accepted;
        return .{
            .iterations = 1,
            .newton_raphson_steps = 0,
            .picard_steps = 0,
            // This branch integrates a bounded kinetic rate; it has no
            // unresolved nonlinear-equilibrium residual after acceptance.
            .maximum_scaled_residual = 0,
        };
    }
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        const changes = try changesAt(current, environment, evaluator);
        const current_norm = try scaledNorm(current, changes, options);
        if (current_norm <= 1) {
            state.cells[cell_index] = current;
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = current_norm };
        }
        var accepted_newton = false;
        // Capped dissolution/adsorption rates are locally constant, so their
        // directional derivative is zero until an inventory bound changes
        // the active set. Jump directly to the nearest such nonnegative
        // boundary and accept only a strict residual improvement.
        const admissible_fraction = maximumAdmissibleFraction(Cell, current, changes);
        if (std.math.isFinite(admissible_fraction) and admissible_fraction > options.maximum_newton_fraction) {
            var best = current;
            var best_norm = current_norm;
            var best_fraction: f64 = 0;
            // A capped-rate plateau can occupy most of the admissible
            // interval and hide a narrow equilibrium transition from
            // derivative probes. Sample the complete conservative segment;
            // the following local Newton step polishes the best branch.
            var sample: u8 = 1;
            while (sample <= 64) : (sample += 1) {
                const trial_fraction = admissible_fraction * @as(f64, @floatFromInt(sample)) / 64;
                if (applyFraction(current, changes, trial_fraction)) |trial| {
                    if (changesAt(trial, environment, evaluator)) |trial_changes| {
                        const trial_norm = try scaledNorm(trial, trial_changes, options);
                        if (trial_norm < best_norm or
                            (!meaningfullyImproves(current_norm, best_norm) and
                                trial_fraction > best_fraction and
                                trial_norm <= current_norm * (1 + 1e-12)))
                        {
                            best = trial;
                            best_norm = trial_norm;
                            best_fraction = trial_fraction;
                        }
                    } else |_| {}
                } else |_| {}
            }
            if (meaningfullyImproves(current_norm, best_norm)) {
                current = best;
                newton_steps += 1;
                accepted_newton = true;
            } else if (try plateauBoundaryLookahead(
                current,
                changes,
                admissible_fraction,
                environment,
                evaluator,
                options,
                current_norm,
            )) |candidate| {
                current = candidate;
                newton_steps += 1;
                accepted_newton = true;
            }
        }
        if (accepted_newton) continue;
        if (!environment.dynamic_salts) {
            if (try conservativeFixedPhosphateActiveSetNewton(
                current,
                environment,
                evaluator,
                options,
                current_norm,
            )) |candidate| {
                current = candidate;
                newton_steps += 1;
                continue;
            }
            if (try conservativeFixedPhosphateComplementarityNewton(
                current,
                environment,
                evaluator,
                options,
                current_norm,
            )) |candidate| {
                current = candidate;
                newton_steps += 1;
                continue;
            }
            if (try resolveIncompatibleFixedPhosphateSolids(
                current,
                environment,
                evaluator,
                options,
                current_norm,
            )) |candidate| {
                current = candidate;
                newton_steps += 1;
                continue;
            }
        }
        if (!environment.dynamic_salts and current_norm < 25) {
            if (try accelerateFixedPhosphateCationBounds(
                current,
                environment,
                evaluator,
                options,
                current_norm,
            )) |candidate| {
                current = candidate;
                newton_steps += 1;
                continue;
            }
        }
        // Prefer the complete conservative phosphate/exchange system.  All
        // five phosphate minerals are coordinates of this one Newton step,
        // matching SOLUTE's same-pre-update evaluation and simultaneous
        // extent accumulation.  The scalar active-set branch below is only
        // a plateau escape when the capped-rate Jacobian has no usable
        // coupled direction.
        if (!environment.dynamic_salts and current_norm < 25) {
            if (try conservativePhosphateNewton(state.allocator, current, environment, evaluator, options, current_norm)) |candidate| {
                current = candidate;
                newton_steps += 1;
                continue;
            }
        }
        if (!environment.dynamic_salts) {
            if (try conservativePhosphateActiveReactionSolve(current, environment, evaluator, options, current_norm, false)) |candidate| {
                current = candidate;
                newton_steps += 1;
                continue;
            }
        }
        if (evaluator.equilibrate_cation_exchange) |equilibrate_exchange| {
            const exchange_target = try equilibrate_exchange(evaluator.context, current);
            var exchange_fraction: f64 = 1;
            var exchange_search: u8 = 0;
            while (exchange_search < 32) : (exchange_search += 1) {
                const exchange_candidate =
                    try interpolateCell(current, exchange_target, exchange_fraction);
                if (changesAt(exchange_candidate, environment, evaluator)) |exchange_candidate_changes| {
                    const exchange_candidate_norm =
                        try scaledNorm(exchange_candidate, exchange_candidate_changes, options);
                    const exchange_candidate_extents =
                        try evaluator.evaluate(
                            evaluator.context,
                            exchange_candidate,
                        );
                    const exchange_current_extents =
                        try evaluator.evaluate(
                            evaluator.context,
                            current,
                        );
                    if (exchange_candidate_norm < current_norm or
                        (exchange_candidate_norm <=
                            current_norm * (1 + 1e-6) and
                            !exchangeResidualsConvergedAtAqueousScale(
                                current,
                                exchange_current_extents.exchange,
                                environment.litter_mass_per_water_volume_Mg_per_m3,
                                options,
                            ) and
                            exchangeResidualsConvergedAtAqueousScale(
                                exchange_candidate,
                                exchange_candidate_extents.exchange,
                                environment.litter_mass_per_water_volume_Mg_per_m3,
                                options,
                            )))
                    {
                        current = exchange_candidate;
                        newton_steps += 1;
                        accepted_newton = true;
                        break;
                    }
                } else |_| {}
                exchange_fraction *= 0.5;
            }
            if (accepted_newton) continue;
        }
        if (try conservativeAmmoniumAssociationSolve(current, environment, evaluator, options, current_norm)) |candidate| {
            current = candidate;
            newton_steps += 1;
            continue;
        }
        if (applyFraction(current, changes, options.directional_probe_fraction)) |probe| {
            if (changesAt(probe, environment, evaluator)) |probe_changes| {
                var numerator: f64 = 0;
                var denominator: f64 = 0;
                directionalProducts(Cell, changes, probe_changes, options.directional_probe_fraction, &numerator, &denominator);
                if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                    const raw_fraction = -numerator / denominator;
                    const fraction = @min(
                        @max(raw_fraction, options.minimum_newton_fraction),
                        maximumAdmissibleFraction(Cell, current, changes),
                    );
                    if (applyFraction(current, changes, fraction)) |candidate| {
                        if (changesAt(candidate, environment, evaluator)) |candidate_changes| {
                            const candidate_norm = try scaledNorm(candidate, candidate_changes, options);
                            if (meaningfullyImproves(
                                current_norm,
                                candidate_norm,
                            )) {
                                current = candidate;
                                newton_steps += 1;
                                accepted_newton = true;
                            }
                        } else |_| {}
                    } else |_| {}
                }
            } else |_| {}
        } else |_| {}
        if (accepted_newton) continue;

        const candidate =
            try boundedPicard(current, changes, options.picard_relaxation);
        // Do not scale stagnation by the largest inventory in the cell. A
        // large carbonate or exchange pool can make a representable,
        // convergence-controlling trace-ion correction look globally tiny.
        if (valuesEqual(Cell, current, candidate)) return error.LitterChemistrySolverStagnated;
        current = candidate;
        picard_steps += 1;
    }

    const final_changes = try changesAt(current, environment, evaluator);
    const final_norm = try scaledNorm(current, final_changes, options);
    // An accepted step on the last permitted iteration may itself satisfy
    // equilibrium. Check that state before reporting ceiling exhaustion.
    if (final_norm <= 1) {
        state.cells[cell_index] = current;
        return .{
            .iterations = options.max_iterations,
            .newton_raphson_steps = newton_steps,
            .picard_steps = picard_steps,
            .maximum_scaled_residual = final_norm,
        };
    }
    std.log.warn("litter chemistry failed to converge: cell={d} iterations={d} maximum_scaled_residual={e} newton_steps={d} picard_steps={d}", .{ cell_index, options.max_iterations, final_norm, newton_steps, picard_steps });
    logUnconvergedFields(Cell, "", current, final_changes, options);
    logAdmissibleLimits(Cell, "", current, final_changes, maximumAdmissibleFraction(Cell, current, final_changes));
    const final_extents = try evaluator.evaluate(evaluator.context, current);
    std.log.warn(
        "litter phosphate mineral extents: AlPO4={e} FePO4={e} CaHPO4={e} hydroxyapatite={e} monocalcium={e}",
        .{
            final_extents.phosphate_minerals.aluminum_phosphate_mol_per_m3,
            final_extents.phosphate_minerals.iron_phosphate_mol_per_m3,
            final_extents.phosphate_minerals.dicalcium_phosphate_mol_per_m3,
            final_extents.phosphate_minerals.hydroxyapatite_mol_per_m3,
            final_extents.phosphate_minerals.monocalcium_phosphate_mol_per_m3,
        },
    );
    std.log.warn(
        "litter salt mineral extents: gibbsite={e} iron_hydroxide={e} calcite={e} gypsum={e}",
        .{
            final_extents.salt_minerals.gibbsite_mol_per_m3,
            final_extents.salt_minerals.iron_hydroxide_mol_per_m3,
            final_extents.salt_minerals.calcite_mol_per_m3,
            final_extents.salt_minerals.gypsum_mol_per_m3,
        },
    );
    std.log.warn(
        "litter exchange extents mol/Mg: NH4={e} H={e} Al={e} Fe={e} Ca={e} Mg={e} Na={e} K={e}",
        .{
            final_extents.exchange.ammonium_mol_per_Mg,
            final_extents.exchange.hydrogen_mol_per_Mg,
            final_extents.exchange.aluminum_mol_per_Mg,
            final_extents.exchange.iron_mol_per_Mg,
            final_extents.exchange.calcium_mol_per_Mg,
            final_extents.exchange.magnesium_mol_per_Mg,
            final_extents.exchange.sodium_mol_per_Mg,
            final_extents.exchange.potassium_mol_per_Mg,
        },
    );
    if (evaluator.phosphate_mineral_equilibrium_residuals) |calculate| {
        const exact = try calculate(evaluator.context, current);
        std.log.warn(
            "litter phosphate saturation residuals: AlPO4={e} FePO4={e} CaHPO4={e} hydroxyapatite={e} monocalcium={e}",
            .{
                exact.aluminum_phosphate_mol_per_m3,
                exact.iron_phosphate_mol_per_m3,
                exact.dicalcium_phosphate_mol_per_m3,
                exact.hydroxyapatite_mol_per_m3,
                exact.monocalcium_phosphate_mol_per_m3,
            },
        );
    }
    return error.LitterChemistrySolverDidNotConverge;
}

/// Applies the largest representable admissible fraction at or below the
/// analytic shared-inventory bound. A mathematically exact endpoint can round
/// one of many simultaneous sinks a few ulps below zero; bisection corrects
/// only that floating-point endpoint error without changing relative rates.
fn applyAdmissibleFraction(
    current: Cell,
    changes: Cell,
    requested_fraction: f64,
) !Cell {
    if (applyFraction(current, changes, requested_fraction)) |candidate|
        return candidate
    else |err| switch (err) {
        error.NegativeLitterChemistryState => {},
        else => return err,
    }
    var lower: f64 = 0;
    var upper = requested_fraction;
    var best = current;
    var iteration: u8 = 0;
    while (iteration < 64) : (iteration += 1) {
        const middle = lower + 0.5 * (upper - lower);
        if (applyFraction(current, changes, middle)) |candidate| {
            lower = middle;
            best = candidate;
        } else |err| switch (err) {
            error.NegativeLitterChemistryState => upper = middle,
            else => return err,
        }
    }
    if (lower <= 0) return error.NoPhysicallyAdmissibleFixedPhKineticStep;
    return best;
}

/// A clipped simultaneous reaction can have a constant residual all the way
/// to an inventory boundary, with the descending active set visible only
/// after that boundary is reached. Evaluate that transition transactionally:
/// neither the neutral boundary nor any look-ahead iterate is published, and
/// return a candidate only when the complete reaction norm improves over the
/// original state. This crosses stoichiometric null plateaus without allowing
/// silent non-descent progression or a sequential mineral commit.
fn plateauBoundaryLookahead(
    current: Cell,
    changes: Cell,
    admissible_fraction: f64,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    if (!std.math.isFinite(admissible_fraction) or
        admissible_fraction <= options.maximum_newton_fraction)
        return null;
    var trial = applyFraction(
        current,
        changes,
        admissible_fraction,
    ) catch return null;
    var best: ?Cell = null;
    var best_norm = current_norm;
    const lookahead_ceiling: u8 = 64;
    var lookahead: u8 = 0;
    while (lookahead < lookahead_ceiling) : (lookahead += 1) {
        const trial_changes = changesAt(
            trial,
            environment,
            evaluator,
        ) catch return best;
        const trial_norm = try scaledNorm(trial, trial_changes, options);
        if (meaningfullyImproves(best_norm, trial_norm)) {
            best = trial;
            best_norm = trial_norm;
        }
        if (trial_norm <= 1) break;
        trial = boundedPicard(
            trial,
            trial_changes,
            options.picard_relaxation,
        ) catch return best;
    }
    return best;
}

const phosphate_coordinate_count: usize = 18;

fn exchangeResidualsConvergedAtAqueousScale(
    cell: Cell,
    extents: ledger.ExchangeAdsorption,
    density_Mg_per_m3: f64,
    options: Options,
) bool {
    inline for (@typeInfo(ledger.ExchangeAdsorption).@"struct".fields) |field| {
        const exchange_inventory = @field(cell.exchange, field.name);
        const exchange_scale = options.absolute_tolerance +
            options.relative_tolerance *
                @max(1.0, @abs(exchange_inventory));
        const aqueous_inventory = exchangeAqueousInventory(
            cell,
            field.name,
        );
        const aqueous_scale_per_Mg =
            (options.absolute_tolerance +
                options.relative_tolerance *
                    @max(1.0, @abs(aqueous_inventory))) /
            density_Mg_per_m3;
        if (@abs(@field(extents, field.name)) >
            @min(exchange_scale, aqueous_scale_per_Mg))
            return false;
    }
    return true;
}

fn exchangeAqueousInventory(cell: Cell, comptime field_name: []const u8) f64 {
    if (comptime std.mem.eql(u8, field_name, "ammonium_mol_per_Mg"))
        return cell.ammonium_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "hydrogen_mol_per_Mg"))
        return cell.hydrogen_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "aluminum_mol_per_Mg"))
        return cell.aluminum_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "iron_mol_per_Mg"))
        return cell.iron_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "calcium_mol_per_Mg"))
        return cell.calcium_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "magnesium_mol_per_Mg"))
        return cell.magnesium_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "sodium_mol_per_Mg"))
        return cell.sodium_mol_per_m3;
    if (comptime std.mem.eql(u8, field_name, "potassium_mol_per_Mg"))
        return cell.potassium_mol_per_m3;
    unreachable;
}

/// Fixed-pH variscite and strengite targets are externally buffered in
/// SOLUTE, so precipitation can terminate at the dissolved Al/Fe bounds
/// instead of at both mutually incompatible saturation targets. Accelerate
/// those two bound-active reactions together from one pre-update evaluation,
/// conserving Al, Fe, and P exactly. The caller accepts only a reduction of
/// the complete reaction norm, so this cannot hide a competing mineral.
const fixed_phosphate_coordinate_count: usize = 6;
const fixed_phosphate_mineral_count: usize = 5;

/// Solves each admissible mineral phase assemblage explicitly. An active
/// mineral is constrained to zero saturation residual; an inactive mineral
/// is fixed at zero inventory and must remain undersaturated. This avoids the
/// nearly singular Fischer-Burmeister Jacobian when several large solid
/// inventories compete for the same finite aqueous phosphate pool.
fn conservativeFixedPhosphateActiveSetNewton(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const calculate =
        evaluator.phosphate_mineral_equilibrium_residuals orelse return null;
    const density = environment.litter_mass_per_water_volume_Mg_per_m3;
    const site_total = phosphateSiteTotal(current);
    const phosphorus_total = phosphateTotal(current, density);
    const cation_totals = phosphateCationTotals(current);
    var best: ?Cell = null;
    var best_norm = current_norm;

    var mask: u8 = 0;
    while (mask < (1 << fixed_phosphate_mineral_count)) : (mask += 1) {
        var coordinates = fixedPhosphateCoordinates(current);
        for (0..fixed_phosphate_mineral_count) |mineral_index| {
            if (!mineralIsActive(mask, mineral_index))
                coordinates[mineral_index + 1] = 0;
        }
        var candidate = fixedPhosphateCandidate(
            current,
            coordinates,
            site_total,
            phosphorus_total,
            cation_totals,
            density,
        ) catch continue;

        var iteration: u8 = 0;
        while (iteration < 24) : (iteration += 1) {
            var residual: [fixed_phosphate_coordinate_count]f64 = undefined;
            const dimension = try fixedPhosphateActiveSetResidual(
                candidate,
                evaluator,
                calculate,
                mask,
                &residual,
            );
            const reduced_norm = activeSetResidualNorm(
                candidate,
                residual,
                dimension,
                mask,
                options,
            );
            if (reduced_norm <= 1) break;

            var matrix: [
                fixed_phosphate_coordinate_count *
                    fixed_phosphate_coordinate_count
            ]f64 = @splat(0);
            var coordinate_slots: [fixed_phosphate_coordinate_count]usize =
                undefined;
            coordinate_slots[0] = 0;
            var slot_count: usize = 1;
            for (0..fixed_phosphate_mineral_count) |mineral_index| {
                if (mineralIsActive(mask, mineral_index)) {
                    coordinate_slots[slot_count] = mineral_index + 1;
                    slot_count += 1;
                }
            }
            std.debug.assert(slot_count == dimension);

            var usable = true;
            for (0..dimension) |column| {
                const coordinate_index = coordinate_slots[column];
                const coordinate_scale = @max(
                    @abs(coordinates[coordinate_index]),
                    options.absolute_tolerance,
                );
                var step = std.math.cbrt(std.math.floatEps(f64)) *
                    coordinate_scale;
                var probe_residual: [fixed_phosphate_coordinate_count]f64 =
                    undefined;
                var selected_step: f64 = 0;
                var search: u8 = 0;
                while (search < 24) : (search += 1) {
                    const signed_step =
                        if (search % 2 == 0) step else -step;
                    var probe_coordinates = coordinates;
                    probe_coordinates[coordinate_index] += signed_step;
                    const probe = fixedPhosphateCandidate(
                        current,
                        probe_coordinates,
                        site_total,
                        phosphorus_total,
                        cation_totals,
                        density,
                    ) catch {
                        if (search % 2 == 1) step *= 2;
                        continue;
                    };
                    _ = fixedPhosphateActiveSetResidual(
                        probe,
                        evaluator,
                        calculate,
                        mask,
                        &probe_residual,
                    ) catch {
                        if (search % 2 == 1) step *= 2;
                        continue;
                    };
                    if (residualProbeIsInformative(
                        residual[0..dimension],
                        probe_residual[0..dimension],
                    )) {
                        selected_step = signed_step;
                        break;
                    }
                    if (search % 2 == 1) step *= 2;
                }
                if (selected_step == 0) {
                    usable = false;
                    break;
                }
                for (0..dimension) |row| {
                    const scale = activeSetResidualScale(
                        candidate,
                        row,
                        mask,
                        options,
                    );
                    matrix[row * dimension + column] =
                        (probe_residual[row] - residual[row]) /
                        selected_step * coordinate_scale / scale;
                }
            }
            if (!usable) break;

            var right_hand_side: [fixed_phosphate_coordinate_count]f64 = @splat(0);
            for (0..dimension) |row|
                right_hand_side[row] = -residual[row] /
                    activeSetResidualScale(candidate, row, mask, options);
            var solved_matrix = matrix;
            var delta = right_hand_side;
            if (!numerics.solveDenseLinearSystem(
                solved_matrix[0 .. dimension * dimension],
                delta[0..dimension],
                dimension,
            )) {
                var normal_matrix: [
                    fixed_phosphate_coordinate_count *
                        fixed_phosphate_coordinate_count
                ]f64 = @splat(0);
                var normal_right_hand_side: [fixed_phosphate_coordinate_count]f64 = @splat(0);
                if (!solveDampedLeastSquares(
                    matrix[0 .. dimension * dimension],
                    right_hand_side[0..dimension],
                    delta[0..dimension],
                    normal_matrix[0 .. dimension * dimension],
                    normal_right_hand_side[0..dimension],
                    dimension,
                )) break;
            }
            for (0..dimension) |column| {
                const coordinate_index = coordinate_slots[column];
                delta[column] *= @max(
                    @abs(coordinates[coordinate_index]),
                    options.absolute_tolerance,
                );
            }

            var accepted = false;
            var fraction: f64 = 1;
            var line_search: u8 = 0;
            while (line_search < 40) : (line_search += 1) {
                var trial_coordinates = coordinates;
                for (0..dimension) |column|
                    trial_coordinates[coordinate_slots[column]] +=
                        fraction * delta[column];
                const trial = fixedPhosphateCandidate(
                    current,
                    trial_coordinates,
                    site_total,
                    phosphorus_total,
                    cation_totals,
                    density,
                ) catch {
                    fraction *= 0.5;
                    continue;
                };
                var trial_residual: [fixed_phosphate_coordinate_count]f64 = undefined;
                _ = fixedPhosphateActiveSetResidual(
                    trial,
                    evaluator,
                    calculate,
                    mask,
                    &trial_residual,
                ) catch {
                    fraction *= 0.5;
                    continue;
                };
                if (activeSetResidualNorm(
                    trial,
                    trial_residual,
                    dimension,
                    mask,
                    options,
                ) < reduced_norm) {
                    coordinates = trial_coordinates;
                    candidate = trial;
                    accepted = true;
                    break;
                }
                fraction *= 0.5;
            }
            if (!accepted) break;
        }

        if (!try fixedPhosphateActiveSetIsAdmissible(
            candidate,
            evaluator,
            calculate,
            mask,
            options,
        )) continue;
        const changes = changesAt(candidate, environment, evaluator) catch
            continue;
        const complete_norm = try scaledNorm(candidate, changes, options);
        if (complete_norm < best_norm) {
            best = candidate;
            best_norm = complete_norm;
        }
    }
    return best;
}

fn mineralIsActive(mask: u8, mineral_index: usize) bool {
    return mask & (@as(u8, 1) << @intCast(mineral_index)) != 0;
}

fn conservativeFixedPhosphateComplementarityNewton(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const calculate =
        evaluator.phosphate_mineral_equilibrium_residuals orelse return null;
    const density = environment.litter_mass_per_water_volume_Mg_per_m3;
    const site_total = phosphateSiteTotal(current);
    const phosphorus_total = phosphateTotal(current, density);
    const cation_totals = phosphateCationTotals(current);
    const coordinates = fixedPhosphateCoordinates(current);
    var base_residual: [fixed_phosphate_coordinate_count]f64 = undefined;
    try fixedPhosphateComplementarityResidual(
        current,
        evaluator,
        calculate,
        &base_residual,
    );
    const base_reduced_norm =
        fixedPhosphateResidualNorm(current, base_residual, options);
    if (base_reduced_norm <= 1) return null;

    const coordinate_scales = [_]f64{
        @max(
            current.hpo4_mol_p_per_m3,
            current.h2po4_mol_p_per_m3,
            options.absolute_tolerance,
        ),
        @max(
            cation_totals.aluminum_mol_per_m3,
            options.absolute_tolerance,
        ),
        @max(
            cation_totals.iron_mol_per_m3,
            options.absolute_tolerance,
        ),
        @max(
            cation_totals.calcium_mol_per_m3,
            options.absolute_tolerance,
        ),
        @max(
            cation_totals.calcium_mol_per_m3 / 5,
            options.absolute_tolerance,
        ),
        @max(
            cation_totals.calcium_mol_per_m3,
            options.absolute_tolerance,
        ),
    };
    var matrix: [
        fixed_phosphate_coordinate_count *
            fixed_phosphate_coordinate_count
    ]f64 = undefined;
    for (0..fixed_phosphate_coordinate_count) |column| {
        var magnitude = std.math.cbrt(std.math.floatEps(f64)) *
            coordinate_scales[column];
        var selected_step: f64 = 0;
        var probe_residual: [fixed_phosphate_coordinate_count]f64 = undefined;
        var search: u8 = 0;
        while (search < 32) : (search += 1) {
            const step = if (search % 2 == 0) magnitude else -magnitude;
            var probe_coordinates = coordinates;
            probe_coordinates[column] += step;
            if (fixedPhosphateCandidate(
                current,
                probe_coordinates,
                site_total,
                phosphorus_total,
                cation_totals,
                density,
            )) |probe| {
                if (fixedPhosphateComplementarityResidual(
                    probe,
                    evaluator,
                    calculate,
                    &probe_residual,
                )) |_| {
                    if (residualProbeIsInformative(
                        &base_residual,
                        &probe_residual,
                    )) {
                        selected_step = step;
                        break;
                    }
                } else |_| {}
            } else |_| {}
            if (search % 2 == 1) magnitude *= 2;
        }
        if (selected_step == 0) return null;
        for (0..fixed_phosphate_coordinate_count) |row| {
            const residual_scale =
                fixedPhosphateResidualScale(current, row, options);
            matrix[row * fixed_phosphate_coordinate_count + column] =
                (probe_residual[row] - base_residual[row]) /
                selected_step * coordinate_scales[column] /
                residual_scale;
        }
    }
    var right_hand_side: [fixed_phosphate_coordinate_count]f64 = undefined;
    for (0..fixed_phosphate_coordinate_count) |row|
        right_hand_side[row] = -base_residual[row] /
            fixedPhosphateResidualScale(current, row, options);
    var solved_matrix = matrix;
    var delta = right_hand_side;
    if (!numerics.solveDenseLinearSystem(
        &solved_matrix,
        &delta,
        fixed_phosphate_coordinate_count,
    )) {
        var normal_matrix: [
            fixed_phosphate_coordinate_count *
                fixed_phosphate_coordinate_count
        ]f64 = undefined;
        var normal_right_hand_side: [fixed_phosphate_coordinate_count]f64 = undefined;
        if (!solveDampedLeastSquares(
            &matrix,
            &right_hand_side,
            &delta,
            &normal_matrix,
            &normal_right_hand_side,
            fixed_phosphate_coordinate_count,
        )) return null;
    }
    for (0..fixed_phosphate_coordinate_count) |coordinate|
        delta[coordinate] *= coordinate_scales[coordinate];

    var fraction: f64 = 1;
    var line_search: u8 = 0;
    while (line_search < 53) : (line_search += 1) {
        var candidate_coordinates = coordinates;
        for (0..fixed_phosphate_coordinate_count) |coordinate|
            candidate_coordinates[coordinate] +=
                fraction * delta[coordinate];
        if (fixedPhosphateCandidate(
            current,
            candidate_coordinates,
            site_total,
            phosphorus_total,
            cation_totals,
            density,
        )) |candidate| {
            var candidate_residual: [fixed_phosphate_coordinate_count]f64 =
                undefined;
            if (fixedPhosphateComplementarityResidual(
                candidate,
                evaluator,
                calculate,
                &candidate_residual,
            )) |_| {
                const candidate_reduced_norm = fixedPhosphateResidualNorm(
                    candidate,
                    candidate_residual,
                    options,
                );
                const candidate_changes =
                    changesAt(candidate, environment, evaluator) catch {
                        fraction *= 0.5;
                        continue;
                    };
                const candidate_norm =
                    try scaledNorm(candidate, candidate_changes, options);
                if (meaningfullyImproves(
                    base_reduced_norm,
                    candidate_reduced_norm,
                ) and candidate_norm <= 10 * current_norm)
                    return candidate;
            } else |_| {}
        } else |_| {}
        fraction *= 0.5;
    }
    return null;
}

fn fixedPhosphateCoordinates(
    cell: Cell,
) [fixed_phosphate_coordinate_count]f64 {
    const mineral = cell.phosphate_minerals;
    return .{
        cell.hpo4_mol_p_per_m3,
        mineral.aluminum_phosphate_mol_per_m3,
        mineral.iron_phosphate_mol_per_m3,
        mineral.dicalcium_phosphate_mol_per_m3,
        mineral.hydroxyapatite_mol_per_m3,
        mineral.monocalcium_phosphate_mol_per_m3,
    };
}

fn fixedPhosphateCandidate(
    current: Cell,
    coordinates: [fixed_phosphate_coordinate_count]f64,
    site_total: f64,
    phosphorus_total: f64,
    cation_totals: PhosphateCationTotals,
    density: f64,
) !Cell {
    var candidate = current;
    candidate.hpo4_mol_p_per_m3 = coordinates[0];
    candidate.phosphate_minerals = .{
        .aluminum_phosphate_mol_per_m3 = coordinates[1],
        .iron_phosphate_mol_per_m3 = coordinates[2],
        .dicalcium_phosphate_mol_per_m3 = coordinates[3],
        .hydroxyapatite_mol_per_m3 = coordinates[4],
        .monocalcium_phosphate_mol_per_m3 = coordinates[5],
    };
    try closePhosphateConservation(
        &candidate,
        site_total,
        phosphorus_total,
        cation_totals,
        density,
    );
    return candidate;
}

fn fixedPhosphateActiveSetResidual(
    cell: Cell,
    evaluator: Evaluator,
    calculate: *const fn (
        context: *const anyopaque,
        cell: Cell,
    ) anyerror!ledger.PhosphateMineralExtents,
    mask: u8,
    residual: *[fixed_phosphate_coordinate_count]f64,
) !usize {
    const extents = try evaluator.evaluate(evaluator.context, cell);
    _ = calculate;
    residual[0] = extents.h2po4_association_mol_p_per_m3;
    var row: usize = 1;
    for (0..fixed_phosphate_mineral_count) |mineral_index| {
        if (mineralIsActive(mask, mineral_index)) {
            residual[row] = phosphateMineralExtent(
                extents.phosphate_minerals,
                @enumFromInt(mineral_index),
            );
            row += 1;
        }
    }
    return row;
}

fn activeSetResidualScale(
    cell: Cell,
    row: usize,
    mask: u8,
    options: Options,
) f64 {
    if (row == 0)
        return options.absolute_tolerance +
            options.relative_tolerance *
                @max(1.0, @abs(cell.hpo4_mol_p_per_m3));
    var active_row: usize = 1;
    for (0..fixed_phosphate_mineral_count) |mineral_index| {
        if (!mineralIsActive(mask, mineral_index)) continue;
        if (active_row == row) {
            const inventory = phosphateMineralInventory(
                cell,
                @enumFromInt(mineral_index),
            );
            return options.absolute_tolerance +
                options.relative_tolerance * @max(1.0, @abs(inventory));
        }
        active_row += 1;
    }
    unreachable;
}

fn activeSetResidualNorm(
    cell: Cell,
    residual: [fixed_phosphate_coordinate_count]f64,
    dimension: usize,
    mask: u8,
    options: Options,
) f64 {
    var maximum: f64 = 0;
    for (0..dimension) |row|
        maximum = @max(
            maximum,
            @abs(residual[row]) /
                activeSetResidualScale(cell, row, mask, options),
        );
    return maximum;
}

fn fixedPhosphateActiveSetIsAdmissible(
    cell: Cell,
    evaluator: Evaluator,
    calculate: *const fn (
        context: *const anyopaque,
        cell: Cell,
    ) anyerror!ledger.PhosphateMineralExtents,
    mask: u8,
    options: Options,
) !bool {
    var residual: [fixed_phosphate_coordinate_count]f64 = undefined;
    const dimension = try fixedPhosphateActiveSetResidual(
        cell,
        evaluator,
        calculate,
        mask,
        &residual,
    );
    if (activeSetResidualNorm(
        cell,
        residual,
        dimension,
        mask,
        options,
    ) > 1) return false;

    const saturation = try calculate(evaluator.context, cell);
    for (0..fixed_phosphate_mineral_count) |mineral_index| {
        const reaction: PhosphateMineralReaction =
            @enumFromInt(mineral_index);
        const inventory = phosphateMineralInventory(cell, reaction);
        const scale = options.absolute_tolerance +
            options.relative_tolerance * @max(1.0, @abs(inventory));
        if (mineralIsActive(mask, mineral_index)) {
            if (inventory <= 0) return false;
        } else {
            if (inventory != 0) return false;
            // A missing solid is admissible only when it has no positive
            // precipitation drive. Negative values denote undersaturation.
            if (phosphateMineralExtent(saturation, reaction) > scale)
                return false;
        }
    }
    return true;
}

fn fixedPhosphateComplementarityResidual(
    cell: Cell,
    evaluator: Evaluator,
    calculate: *const fn (
        context: *const anyopaque,
        cell: Cell,
    ) anyerror!ledger.PhosphateMineralExtents,
    residual: *[fixed_phosphate_coordinate_count]f64,
) !void {
    const extents = try evaluator.evaluate(evaluator.context, cell);
    const saturation = try calculate(evaluator.context, cell);
    const mineral = cell.phosphate_minerals;
    residual.* = .{
        extents.h2po4_association_mol_p_per_m3,
        fischerBurmeister(
            mineral.aluminum_phosphate_mol_per_m3,
            -saturation.aluminum_phosphate_mol_per_m3,
        ),
        fischerBurmeister(
            mineral.iron_phosphate_mol_per_m3,
            -saturation.iron_phosphate_mol_per_m3,
        ),
        fischerBurmeister(
            mineral.dicalcium_phosphate_mol_per_m3,
            -saturation.dicalcium_phosphate_mol_per_m3,
        ),
        fischerBurmeister(
            mineral.hydroxyapatite_mol_per_m3,
            -saturation.hydroxyapatite_mol_per_m3,
        ),
        fischerBurmeister(
            mineral.monocalcium_phosphate_mol_per_m3,
            -saturation.monocalcium_phosphate_mol_per_m3,
        ),
    };
}

fn fischerBurmeister(nonnegative_left: f64, nonnegative_right: f64) f64 {
    return std.math.hypot(nonnegative_left, nonnegative_right) -
        nonnegative_left - nonnegative_right;
}

fn fixedPhosphateResidualNorm(
    cell: Cell,
    residual: [fixed_phosphate_coordinate_count]f64,
    options: Options,
) f64 {
    var maximum: f64 = 0;
    for (residual, 0..) |value, row|
        maximum = @max(
            maximum,
            @abs(value) /
                fixedPhosphateResidualScale(cell, row, options),
        );
    return maximum;
}

fn fixedPhosphateResidualScale(
    cell: Cell,
    row: usize,
    options: Options,
) f64 {
    const coordinates = fixedPhosphateCoordinates(cell);
    return options.absolute_tolerance +
        options.relative_tolerance *
            @max(1.0, @abs(coordinates[row]));
}

fn resolveIncompatibleFixedPhosphateSolids(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const calculate =
        evaluator.phosphate_mineral_equilibrium_residuals orelse return null;
    const residuals = try calculate(evaluator.context, current);
    const scale = options.absolute_tolerance +
        options.relative_tolerance;
    const aluminum_solid =
        current.phosphate_minerals.aluminum_phosphate_mol_per_m3;
    const iron_solid =
        current.phosphate_minerals.iron_phosphate_mol_per_m3;
    if (aluminum_solid <= scale or iron_solid <= scale) return null;

    if (residuals.aluminum_phosphate_mol_per_m3 < -scale and
        residuals.iron_phosphate_mol_per_m3 > scale)
        return fixedPhosphateBoundaryCandidate(
            current,
            environment,
            evaluator,
            options,
            current_norm,
            .aluminum,
        );
    if (residuals.iron_phosphate_mol_per_m3 < -scale and
        residuals.aluminum_phosphate_mol_per_m3 > scale)
        return fixedPhosphateBoundaryCandidate(
            current,
            environment,
            evaluator,
            options,
            current_norm,
            .iron,
        );
    return null;
}

const IncompatiblePhosphateSolid = enum { aluminum, iron };

fn fixedPhosphateBoundaryCandidate(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
    exhausted: IncompatiblePhosphateSolid,
) !?Cell {
    const calculate =
        evaluator.phosphate_mineral_equilibrium_residuals orelse return null;
    var base = current;
    const released = switch (exhausted) {
        .aluminum => current.phosphate_minerals
            .aluminum_phosphate_mol_per_m3,
        .iron => current.phosphate_minerals.iron_phosphate_mol_per_m3,
    };
    switch (exhausted) {
        .aluminum => {
            base.phosphate_minerals.aluminum_phosphate_mol_per_m3 = 0;
            base.aluminum_mol_per_m3 += released;
        },
        .iron => {
            base.phosphate_minerals.iron_phosphate_mol_per_m3 = 0;
            base.iron_mol_per_m3 += released;
        },
    }
    base.h2po4_mol_p_per_m3 += released;
    try validateCell(base);

    const active_residual = struct {
        fn value(
            which: IncompatiblePhosphateSolid,
            residual: ledger.PhosphateMineralExtents,
        ) f64 {
            return switch (which) {
                .aluminum => residual.iron_phosphate_mol_per_m3,
                .iron => residual.aluminum_phosphate_mol_per_m3,
            };
        }
    }.value;
    const lower_residual = active_residual(
        exhausted,
        try calculate(evaluator.context, base),
    );
    if (lower_residual <= 0) return null;
    const active_aqueous = switch (exhausted) {
        .aluminum => base.iron_mol_per_m3,
        .iron => base.aluminum_mol_per_m3,
    };
    var upper = @min(base.h2po4_mol_p_per_m3, active_aqueous);
    upper *= 1 - 64 * std.math.floatEps(f64);
    if (upper <= 0) return null;
    var upper_candidate =
        try applyCompetingFixedPhosphateExtent(base, exhausted, upper);
    const upper_residual = active_residual(
        exhausted,
        try calculate(evaluator.context, upper_candidate),
    );
    if (upper_residual > 0) return null;

    var lower: f64 = 0;
    var candidate = base;
    var iteration: u8 = 0;
    while (iteration < 96) : (iteration += 1) {
        const middle = lower + 0.5 * (upper - lower);
        const middle_candidate =
            try applyCompetingFixedPhosphateExtent(base, exhausted, middle);
        const middle_residual = active_residual(
            exhausted,
            try calculate(evaluator.context, middle_candidate),
        );
        candidate = middle_candidate;
        if (@abs(middle_residual) <=
            options.absolute_tolerance + options.relative_tolerance)
            break;
        if (middle_residual > 0) {
            lower = middle;
        } else {
            upper = middle;
            upper_candidate = middle_candidate;
        }
        if (upper - lower <= std.math.floatEps(f64) *
            @max(1.0, upper))
        {
            candidate = upper_candidate;
            break;
        }
    }
    const candidate_residuals = try calculate(evaluator.context, candidate);
    const exhausted_residual = switch (exhausted) {
        .aluminum => candidate_residuals
            .aluminum_phosphate_mol_per_m3,
        .iron => candidate_residuals.iron_phosphate_mol_per_m3,
    };
    const remaining_residual = active_residual(
        exhausted,
        candidate_residuals,
    );
    if (exhausted_residual > options.absolute_tolerance +
        options.relative_tolerance or
        @abs(remaining_residual) > options.absolute_tolerance +
            options.relative_tolerance)
        return null;
    const candidate_changes =
        try changesAt(candidate, environment, evaluator);
    const candidate_norm =
        try scaledNorm(candidate, candidate_changes, options);
    return if (candidate_norm <= 10 * current_norm) candidate else null;
}

fn applyCompetingFixedPhosphateExtent(
    base: Cell,
    exhausted: IncompatiblePhosphateSolid,
    extent: f64,
) !Cell {
    var candidate = base;
    switch (exhausted) {
        .aluminum => {
            candidate.phosphate_minerals.iron_phosphate_mol_per_m3 +=
                extent;
            candidate.iron_mol_per_m3 -= extent;
        },
        .iron => {
            candidate.phosphate_minerals.aluminum_phosphate_mol_per_m3 +=
                extent;
            candidate.aluminum_mol_per_m3 -= extent;
        },
    }
    candidate.h2po4_mol_p_per_m3 -= extent;
    normalizeRoundoffNonnegative(Cell, &candidate);
    try validateCell(candidate);
    return candidate;
}

fn accelerateFixedPhosphateCationBounds(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const extents = try evaluator.evaluate(evaluator.context, current);
    const equilibrium_residuals_callback =
        evaluator.phosphate_mineral_equilibrium_residuals orelse return null;
    const equilibrium_residuals =
        try equilibrium_residuals_callback(evaluator.context, current);
    const aluminum_extent =
        extents.phosphate_minerals.aluminum_phosphate_mol_per_m3;
    const iron_extent =
        extents.phosphate_minerals.iron_phosphate_mol_per_m3;
    var candidate = current;
    var total_phosphate_consumption: f64 = 0;
    var accelerated = false;

    inline for (.{
        .{
            "aluminum_mol_per_m3",
            "aluminum_phosphate_mol_per_m3",
            aluminum_extent,
            equilibrium_residuals.aluminum_phosphate_mol_per_m3,
        },
        .{
            "iron_mol_per_m3",
            "iron_phosphate_mol_per_m3",
            iron_extent,
            equilibrium_residuals.iron_phosphate_mol_per_m3,
        },
    }) |entry| {
        const dissolved = @field(current, entry[0]);
        const solid = @field(current.phosphate_minerals, entry[1]);
        const extent = entry[2];
        const saturation_residual = entry[3];
        const extent_scale = options.absolute_tolerance +
            options.relative_tolerance * @max(1.0, @abs(solid));
        if (dissolved > 0 and extent > extent_scale and
            saturation_residual > 0)
        {
            const rate_fraction = extent / dissolved;
            if (std.math.isFinite(rate_fraction) and rate_fraction > 0 and
                rate_fraction <= 1)
            {
                const target_dissolved =
                    @min(dissolved, 0.5 * extent_scale / rate_fraction);
                const consumed = dissolved - target_dissolved;
                if (consumed > 0) {
                    @field(candidate, entry[0]) = target_dissolved;
                    @field(candidate.phosphate_minerals, entry[1]) =
                        solid + consumed;
                    total_phosphate_consumption += consumed;
                    accelerated = true;
                }
            }
        }
    }
    if (!accelerated or
        total_phosphate_consumption > candidate.h2po4_mol_p_per_m3)
        return null;
    candidate.h2po4_mol_p_per_m3 -= total_phosphate_consumption;
    try validateCell(candidate);
    const candidate_changes = try changesAt(candidate, environment, evaluator);
    const candidate_norm = try scaledNorm(candidate, candidate_changes, options);
    // Moving directly to a reactant bound can expose the phosphate
    // association residual that the next coupled Newton step must remove.
    // Permit that bounded active-set transition only when it resolves the
    // controlling Al/Fe rate and keeps the full merit within one decade.
    const candidate_aluminum_scaled = scaledPhosphateMineralResidual(
        candidate_changes.phosphate_minerals,
        candidate,
        options,
        .aluminum_phosphate,
    );
    const candidate_iron_scaled = scaledPhosphateMineralResidual(
        candidate_changes.phosphate_minerals,
        candidate,
        options,
        .iron_phosphate,
    );
    return if ((candidate_aluminum_scaled <= 1 and
        candidate_iron_scaled <= 1) and
        candidate_norm <= 10 * current_norm)
        candidate
    else
        null;
}

fn conservativeAmmoniumAssociationSolve(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const base_extents = try evaluator.evaluate(evaluator.context, current);
    const base_extent = base_extents.ammonium_association_mol_per_m3;
    // Association changes NH4+ and NH3 by the same magnitude, so it has
    // converged only when that extent satisfies both state-variable scales.
    // Using the larger NH4+ inventory alone can silently leave a smaller
    // NH3 pool above its requested relative tolerance.
    const ammonium_scale = options.absolute_tolerance +
        options.relative_tolerance * @max(1.0, @abs(current.ammonium_mol_per_m3));
    const ammonia_scale = options.absolute_tolerance +
        options.relative_tolerance * @max(1.0, @abs(current.ammonia_mol_per_m3));
    const scale = @min(ammonium_scale, ammonia_scale);
    const exchange_scale = options.absolute_tolerance +
        options.relative_tolerance *
            @max(1.0, @abs(current.exchange.ammonium_mol_per_Mg));
    if (@abs(base_extent) <= scale and
        @abs(base_extents.exchange.ammonium_mol_per_Mg) <= exchange_scale)
        return null;
    const density = environment.litter_mass_per_water_volume_Mg_per_m3;
    const total_nitrogen = current.ammonium_mol_per_m3 +
        current.ammonia_mol_per_m3 +
        density * current.exchange.ammonium_mol_per_Mg;
    if (!std.math.isFinite(total_nitrogen) or total_nitrogen <= 0) return null;

    var lower_exchange: f64 = 0;
    var upper_exchange = total_nitrogen / density;
    const lower_candidate = try ammoniumCandidateAtExchange(
        current,
        evaluator,
        total_nitrogen,
        density,
        lower_exchange,
    );
    var lower_residual =
        (try evaluator.evaluate(evaluator.context, lower_candidate))
            .exchange.ammonium_mol_per_Mg;
    const upper_candidate = try ammoniumCandidateAtExchange(
        current,
        evaluator,
        total_nitrogen,
        density,
        upper_exchange,
    );
    const upper_residual =
        (try evaluator.evaluate(evaluator.context, upper_candidate))
            .exchange.ammonium_mol_per_Mg;
    var candidate = lower_candidate;
    if (lower_residual == 0) {
        candidate = lower_candidate;
    } else if (upper_residual == 0) {
        candidate = upper_candidate;
    } else {
        if (std.math.signbit(lower_residual) == std.math.signbit(upper_residual)) return null;
        var iteration: u8 = 0;
        while (iteration < 96) : (iteration += 1) {
            const middle_exchange =
                lower_exchange + 0.5 * (upper_exchange - lower_exchange);
            const middle_candidate = try ammoniumCandidateAtExchange(
                current,
                evaluator,
                total_nitrogen,
                density,
                middle_exchange,
            );
            const middle_residual =
                (try evaluator.evaluate(evaluator.context, middle_candidate))
                    .exchange.ammonium_mol_per_Mg;
            candidate = middle_candidate;
            if (middle_residual == 0 or
                upper_exchange - lower_exchange <=
                    8 * std.math.floatEps(f64) *
                        @max(1.0, @abs(middle_exchange)))
                break;
            if (std.math.signbit(middle_residual) == std.math.signbit(lower_residual)) {
                lower_exchange = middle_exchange;
                lower_residual = middle_residual;
            } else {
                upper_exchange = middle_exchange;
            }
        }
    }
    const candidate_extents = try evaluator.evaluate(evaluator.context, candidate);
    const candidate_extent = candidate_extents.ammonium_association_mol_per_m3;
    const candidate_exchange_extent =
        candidate_extents.exchange.ammonium_mol_per_Mg;
    const candidate_changes = try changesAt(candidate, environment, evaluator);
    const candidate_norm = try scaledNorm(candidate, candidate_changes, options);
    return if (candidate_norm < current_norm or
        (candidate_norm <= 10 * current_norm and
            @abs(candidate_extent) <= scale and
            @abs(candidate_exchange_extent) <= exchange_scale))
        candidate
    else
        null;
}

fn ammoniumCandidateAtExchange(
    reference: Cell,
    evaluator: Evaluator,
    total_nitrogen: f64,
    density: f64,
    exchange_ammonium_mol_per_Mg: f64,
) !Cell {
    var candidate = reference;
    candidate.exchange.ammonium_mol_per_Mg = exchange_ammonium_mol_per_Mg;
    const aqueous_total =
        total_nitrogen - density * exchange_ammonium_mol_per_Mg;
    if (!std.math.isFinite(aqueous_total) or aqueous_total < -1e-12)
        return error.NegativeLitterChemistryState;
    var lower: f64 = 0;
    var upper = @max(0, aqueous_total);
    var lower_residual =
        try ammoniumAssociationResidualAt(candidate, evaluator, upper, lower);
    const upper_residual =
        try ammoniumAssociationResidualAt(candidate, evaluator, upper, upper);
    if (lower_residual == 0) {
        upper = lower;
    } else if (upper_residual != 0) {
        if (std.math.signbit(lower_residual) == std.math.signbit(upper_residual))
            return error.UnbracketedLitterAmmoniumEquilibrium;
        var iteration: u8 = 0;
        while (iteration < 96) : (iteration += 1) {
            const middle = lower + 0.5 * (upper - lower);
            const middle_residual =
                try ammoniumAssociationResidualAt(candidate, evaluator, aqueous_total, middle);
            if (middle_residual == 0 or
                upper - lower <= 8 * std.math.floatEps(f64) *
                    @max(1.0, @abs(middle)))
            {
                lower = middle;
                upper = middle;
                break;
            }
            if (std.math.signbit(middle_residual) == std.math.signbit(lower_residual)) {
                lower = middle;
                lower_residual = middle_residual;
            } else {
                upper = middle;
            }
        }
    }
    candidate.ammonium_mol_per_m3 = lower + 0.5 * (upper - lower);
    candidate.ammonia_mol_per_m3 =
        @max(0, aqueous_total - candidate.ammonium_mol_per_m3);
    try validateCell(candidate);
    return candidate;
}

fn ammoniumAssociationResidualAt(
    reference: Cell,
    evaluator: Evaluator,
    aqueous_total: f64,
    ammonium_mol_per_m3: f64,
) !f64 {
    var candidate = reference;
    candidate.ammonium_mol_per_m3 = ammonium_mol_per_m3;
    candidate.ammonia_mol_per_m3 = aqueous_total - ammonium_mol_per_m3;
    try validateCell(candidate);
    const residual =
        (try evaluator.evaluate(evaluator.context, candidate)).ammonium_association_mol_per_m3;
    if (!std.math.isFinite(residual)) return error.NonFiniteLitterChemistryChange;
    return residual;
}

const PhosphateMineralReaction = enum(u8) {
    aluminum_phosphate,
    iron_phosphate,
    dicalcium_phosphate,
    hydroxyapatite,
    monocalcium_phosphate,

    fn controlledByH2po4(self: PhosphateMineralReaction) bool {
        return self != .dicalcium_phosphate;
    }
};

const PhosphateReactionResidual = union(enum) {
    association,
    mineral: PhosphateMineralReaction,
};

/// Capped mineral rates have a zero local derivative throughout most of
/// their admissible interval. Resolve that complementarity plateau in the
/// two aqueous phosphate coordinates, while aluminum phosphate is the
/// dependent conservative inventory. This retains the exact reaction
/// equations instead of regularizing rank-deficient solid coordinates.
fn conservativePhosphateActiveReactionSolve(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
    allow_non_improving_target: bool,
) !?Cell {
    const extents = try evaluator.evaluate(evaluator.context, current);
    const active_mineral = dominantPhosphateMineral(
        extents.phosphate_minerals,
        current,
        options,
    );
    const mineral_extent = phosphateMineralExtent(extents.phosphate_minerals, active_mineral);
    const mineral_scale = options.absolute_tolerance +
        options.relative_tolerance * @max(1.0, @abs(phosphateMineralInventory(current, active_mineral)));
    if (@abs(mineral_extent) <= mineral_scale) return null;

    const density = environment.litter_mass_per_water_volume_Mg_per_m3;
    const phosphorus_total = phosphateTotal(current, density);
    const aqueous_plus_slack = current.hpo4_mol_p_per_m3 +
        current.h2po4_mol_p_per_m3 +
        phosphatePerMineral(active_mineral) *
            phosphateMineralInventory(current, active_mineral);
    if (!std.math.isFinite(aqueous_plus_slack) or aqueous_plus_slack <= 0) return null;

    var candidate = current;
    if (active_mineral.controlledByH2po4()) {
        const maybe_h2po4 = try bracketPhosphateReactionRoot(
            current,
            phosphorus_total,
            density,
            evaluator,
            active_mineral,
            .{ .mineral = active_mineral },
            true,
            current.hpo4_mol_p_per_m3,
            aqueous_plus_slack - current.hpo4_mol_p_per_m3,
        );
        if (maybe_h2po4 == null) {
            return null;
        }
        const h2po4 = maybe_h2po4.?;
        const maybe_hpo4 = try bracketPhosphateReactionRoot(
            current,
            phosphorus_total,
            density,
            evaluator,
            active_mineral,
            .association,
            false,
            h2po4,
            aqueous_plus_slack - h2po4,
        );
        if (maybe_hpo4 == null) {
            return null;
        }
        const hpo4 = maybe_hpo4.?;
        candidate = try phosphateCandidateFromAqueous(
            current,
            phosphorus_total,
            density,
            active_mineral,
            hpo4,
            h2po4,
        );
    } else {
        const hpo4 = try bracketPhosphateReactionRoot(
            current,
            phosphorus_total,
            density,
            evaluator,
            active_mineral,
            .{ .mineral = active_mineral },
            false,
            current.h2po4_mol_p_per_m3,
            aqueous_plus_slack - current.h2po4_mol_p_per_m3,
        ) orelse return null;
        const h2po4 = try bracketPhosphateReactionRoot(
            current,
            phosphorus_total,
            density,
            evaluator,
            active_mineral,
            .association,
            true,
            hpo4,
            aqueous_plus_slack - hpo4,
        ) orelse return null;
        candidate = try phosphateCandidateFromAqueous(
            current,
            phosphorus_total,
            density,
            active_mineral,
            hpo4,
            h2po4,
        );
    }
    candidate = try eliminateUndersaturatedPhosphateMinerals(
        candidate,
        phosphorus_total,
        density,
        active_mineral,
        evaluator,
    );
    const candidate_changes = try changesAt(candidate, environment, evaluator);
    const candidate_norm = try scaledNorm(candidate, candidate_changes, options);
    const candidate_extents = try evaluator.evaluate(evaluator.context, candidate);
    const selected_extent_after =
        phosphateMineralExtent(candidate_extents.phosphate_minerals, active_mineral);
    return if (allow_non_improving_target or
        candidate_norm < current_norm or
        (candidate_norm <= current_norm * (1 + 1e-12) and
            @abs(selected_extent_after) <= mineral_scale))
        candidate
    else
        null;
}

fn dominantPhosphateMineral(
    extents: ledger.PhosphateMineralExtents,
    state_value: Cell,
    options: Options,
) PhosphateMineralReaction {
    var result: PhosphateMineralReaction = .aluminum_phosphate;
    var magnitude = scaledPhosphateMineralResidual(
        extents,
        state_value,
        options,
        result,
    );
    inline for ([_]PhosphateMineralReaction{
        .iron_phosphate,
        .dicalcium_phosphate,
        .hydroxyapatite,
        .monocalcium_phosphate,
    }) |reaction| {
        const candidate = scaledPhosphateMineralResidual(
            extents,
            state_value,
            options,
            reaction,
        );
        if (candidate > magnitude) {
            result = reaction;
            magnitude = candidate;
        }
    }
    return result;
}

fn scaledPhosphateMineralResidual(
    extents: ledger.PhosphateMineralExtents,
    state_value: Cell,
    options: Options,
    reaction: PhosphateMineralReaction,
) f64 {
    const inventory = phosphateMineralInventory(state_value, reaction);
    const scale = options.absolute_tolerance +
        options.relative_tolerance * @max(1.0, @abs(inventory));
    return @abs(phosphateMineralExtent(extents, reaction)) / scale;
}

fn phosphateMineralExtent(extents: ledger.PhosphateMineralExtents, reaction: PhosphateMineralReaction) f64 {
    return switch (reaction) {
        .aluminum_phosphate => extents.aluminum_phosphate_mol_per_m3,
        .iron_phosphate => extents.iron_phosphate_mol_per_m3,
        .dicalcium_phosphate => extents.dicalcium_phosphate_mol_per_m3,
        .hydroxyapatite => extents.hydroxyapatite_mol_per_m3,
        .monocalcium_phosphate => extents.monocalcium_phosphate_mol_per_m3,
    };
}

fn phosphateMineralInventory(cell: Cell, reaction: PhosphateMineralReaction) f64 {
    return phosphateMineralExtent(cell.phosphate_minerals, reaction);
}

fn setPhosphateMineralInventory(cell: *Cell, reaction: PhosphateMineralReaction, value: f64) void {
    switch (reaction) {
        .aluminum_phosphate => cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 = value,
        .iron_phosphate => cell.phosphate_minerals.iron_phosphate_mol_per_m3 = value,
        .dicalcium_phosphate => cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3 = value,
        .hydroxyapatite => cell.phosphate_minerals.hydroxyapatite_mol_per_m3 = value,
        .monocalcium_phosphate => cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 = value,
    }
}

fn phosphatePerMineral(reaction: PhosphateMineralReaction) f64 {
    return switch (reaction) {
        .aluminum_phosphate, .iron_phosphate, .dicalcium_phosphate => 1,
        .hydroxyapatite => 3,
        .monocalcium_phosphate => 2,
    };
}

fn bracketPhosphateReactionRoot(
    reference: Cell,
    phosphorus_total: f64,
    density: f64,
    evaluator: Evaluator,
    slack_mineral: PhosphateMineralReaction,
    residual_kind: PhosphateReactionResidual,
    solve_h2po4: bool,
    fixed_aqueous_phosphate: f64,
    upper_bound: f64,
) !?f64 {
    if (!std.math.isFinite(upper_bound) or upper_bound <= 0) return null;
    var lower: f64 = 0;
    var upper = upper_bound;
    var lower_residual = phosphateReactionResidualAt(
        reference,
        phosphorus_total,
        density,
        evaluator,
        slack_mineral,
        residual_kind,
        solve_h2po4,
        fixed_aqueous_phosphate,
        lower,
    ) catch |err| switch (err) {
        error.NegativeLitterChemistryState => return null,
        else => return err,
    };
    var upper_residual: f64 = undefined;
    var upper_is_admissible = false;
    var endpoint_attempt: u8 = 0;
    while (endpoint_attempt < 64) : (endpoint_attempt += 1) {
        upper_residual = phosphateReactionResidualAt(
            reference,
            phosphorus_total,
            density,
            evaluator,
            slack_mineral,
            residual_kind,
            solve_h2po4,
            fixed_aqueous_phosphate,
            upper,
        ) catch |err| switch (err) {
            error.NegativeLitterChemistryState => {
                upper = lower + 0.5 * (upper - lower);
                continue;
            },
            else => return err,
        };
        upper_is_admissible = true;
        break;
    }
    if (!upper_is_admissible) return null;
    if (lower_residual == 0) return lower;
    if (upper_residual == 0) return upper;
    if (std.math.signbit(lower_residual) == std.math.signbit(upper_residual)) return null;

    var iteration: u8 = 0;
    while (iteration < 96) : (iteration += 1) {
        const middle = lower + 0.5 * (upper - lower);
        const middle_residual = try phosphateReactionResidualAt(
            reference,
            phosphorus_total,
            density,
            evaluator,
            slack_mineral,
            residual_kind,
            solve_h2po4,
            fixed_aqueous_phosphate,
            middle,
        );
        if (middle_residual == 0 or
            upper - lower <= 8 * std.math.floatEps(f64) * @max(1.0, @abs(middle)))
            return middle;
        if (std.math.signbit(middle_residual) == std.math.signbit(lower_residual)) {
            lower = middle;
            lower_residual = middle_residual;
        } else {
            upper = middle;
        }
    }
    return lower + 0.5 * (upper - lower);
}

fn phosphateReactionResidualAt(
    reference: Cell,
    phosphorus_total: f64,
    density: f64,
    evaluator: Evaluator,
    slack_mineral: PhosphateMineralReaction,
    residual_kind: PhosphateReactionResidual,
    solve_h2po4: bool,
    fixed_aqueous_phosphate: f64,
    variable_aqueous_phosphate: f64,
) !f64 {
    const candidate = if (solve_h2po4)
        try phosphateCandidateFromAqueous(
            reference,
            phosphorus_total,
            density,
            slack_mineral,
            fixed_aqueous_phosphate,
            variable_aqueous_phosphate,
        )
    else
        try phosphateCandidateFromAqueous(
            reference,
            phosphorus_total,
            density,
            slack_mineral,
            variable_aqueous_phosphate,
            fixed_aqueous_phosphate,
        );
    const extents = try evaluator.evaluate(evaluator.context, candidate);
    const residual = switch (residual_kind) {
        .association => extents.h2po4_association_mol_p_per_m3,
        .mineral => |reaction| phosphateMineralExtent(extents.phosphate_minerals, reaction),
    };
    if (!std.math.isFinite(residual)) return error.NonFiniteLitterChemistryChange;
    return residual;
}

fn phosphateCandidateFromAqueous(
    reference: Cell,
    phosphorus_total: f64,
    density: f64,
    slack_mineral: PhosphateMineralReaction,
    hpo4_mol_p_per_m3: f64,
    h2po4_mol_p_per_m3: f64,
) !Cell {
    var candidate = reference;
    candidate.hpo4_mol_p_per_m3 = hpo4_mol_p_per_m3;
    candidate.h2po4_mol_p_per_m3 = h2po4_mol_p_per_m3;
    setPhosphateMineralInventory(&candidate, slack_mineral, 0);
    const without_slack = phosphateTotal(candidate, density);
    setPhosphateMineralInventory(
        &candidate,
        slack_mineral,
        (phosphorus_total - without_slack) / phosphatePerMineral(slack_mineral),
    );
    normalizeRoundoffNonnegative(Cell, &candidate);
    try validateCell(candidate);
    return candidate;
}

fn eliminateUndersaturatedPhosphateMinerals(
    initial: Cell,
    phosphorus_total: f64,
    density: f64,
    slack_mineral: PhosphateMineralReaction,
    evaluator: Evaluator,
) !Cell {
    var candidate = initial;
    const extents = (try evaluator.evaluate(evaluator.context, candidate)).phosphate_minerals;
    inline for ([_]PhosphateMineralReaction{
        .aluminum_phosphate,
        .iron_phosphate,
        .dicalcium_phosphate,
        .hydroxyapatite,
        .monocalcium_phosphate,
    }) |reaction| {
        if (reaction != slack_mineral and
            phosphateMineralExtent(extents, reaction) < 0 and
            phosphateMineralInventory(candidate, reaction) > 0)
            setPhosphateMineralInventory(&candidate, reaction, 0);
    }
    const slack_before = phosphateMineralInventory(candidate, slack_mineral);
    setPhosphateMineralInventory(&candidate, slack_mineral, 0);
    const without_slack = phosphateTotal(candidate, density);
    const slack_after =
        (phosphorus_total - without_slack) / phosphatePerMineral(slack_mineral);
    if (!std.math.isFinite(slack_after) or slack_after < -1e-12) {
        setPhosphateMineralInventory(&candidate, slack_mineral, slack_before);
        return initial;
    }
    setPhosphateMineralInventory(&candidate, slack_mineral, @max(0, slack_after));
    try validateCell(candidate);
    return candidate;
}

const reduced_active_coordinate_count: usize = 4;

/// Semismooth Newton step for the fixed-pH active set exposed by dry,
/// high-density litter: phosphate association, simultaneous AlPO4
/// dissolution/FePO4 precipitation, and competitive cation exchange. The
/// remaining phosphate minerals stay at their current complementarity bounds.
/// Every coordinate is reconstructed from elemental, phosphate-site, and
/// exchange-charge invariants before its residual is evaluated.
fn reducedActivePhosphateExchangeNewton(
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const equilibrate_exchange =
        evaluator.equilibrate_cation_exchange orelse return null;
    _ = evaluator.phosphate_mineral_equilibrium_residuals orelse return null;
    const density = environment.litter_mass_per_water_volume_Mg_per_m3;
    const extents = try evaluator.evaluate(evaluator.context, current);
    const mineral_scale = options.absolute_tolerance +
        options.relative_tolerance;
    if (extents.phosphate_minerals.aluminum_phosphate_mol_per_m3 >=
        -mineral_scale or
        extents.phosphate_minerals.iron_phosphate_mol_per_m3 <=
            mineral_scale or
        @abs(extents.exchange.calcium_mol_per_Mg) <=
            mineral_scale / density or
        @abs(extents.phosphate_minerals.dicalcium_phosphate_mol_per_m3) >
            mineral_scale or
        @abs(extents.phosphate_minerals.hydroxyapatite_mol_per_m3) >
            mineral_scale or
        @abs(extents.phosphate_minerals.monocalcium_phosphate_mol_per_m3) >
            mineral_scale)
        return null;

    const exchange_target =
        try equilibrate_exchange(evaluator.context, current);
    const site_total = phosphateSiteTotal(current);
    const phosphorus_total = phosphateTotal(current, density);
    const coordinates = [_]f64{
        current.hpo4_mol_p_per_m3,
        current.phosphate_minerals.aluminum_phosphate_mol_per_m3,
        current.phosphate_minerals.iron_phosphate_mol_per_m3,
        0.5,
    };
    const coordinate_scales = [_]f64{
        @max(
            current.hpo4_mol_p_per_m3,
            current.h2po4_mol_p_per_m3,
            options.absolute_tolerance,
        ),
        @max(
            current.phosphate_minerals.aluminum_phosphate_mol_per_m3,
            current.aluminum_mol_per_m3,
            options.absolute_tolerance,
        ),
        @max(
            current.phosphate_minerals.iron_phosphate_mol_per_m3,
            current.iron_mol_per_m3,
            options.absolute_tolerance,
        ),
        1,
    };
    const base_state = try reducedActiveCandidate(
        current,
        exchange_target,
        coordinates,
        site_total,
        phosphorus_total,
        density,
    );
    var base_residual: [reduced_active_coordinate_count]f64 = undefined;
    try reducedActiveResidual(
        base_state,
        evaluator,
        density,
        &base_residual,
    );
    var matrix: [
        reduced_active_coordinate_count *
            reduced_active_coordinate_count
    ]f64 = undefined;
    var probe_coordinates = coordinates;
    for (0..reduced_active_coordinate_count) |column| {
        const step = std.math.cbrt(std.math.floatEps(f64)) *
            coordinate_scales[column];
        var found_probe = false;
        var signed_step = step;
        var search: u8 = 0;
        while (search < 24) : (search += 1) {
            probe_coordinates = coordinates;
            probe_coordinates[column] += signed_step;
            if (reducedActiveCandidate(
                current,
                exchange_target,
                probe_coordinates,
                site_total,
                phosphorus_total,
                density,
            )) |probe| {
                var probe_residual: [reduced_active_coordinate_count]f64 =
                    undefined;
                if (reducedActiveResidual(
                    probe,
                    evaluator,
                    density,
                    &probe_residual,
                )) |_| {
                    for (0..reduced_active_coordinate_count) |row| {
                        matrix[row * reduced_active_coordinate_count + column] =
                            (probe_residual[row] - base_residual[row]) /
                            signed_step * coordinate_scales[column] /
                            reducedActiveResidualScale(
                                current,
                                row,
                                density,
                                options,
                            );
                    }
                    found_probe = true;
                    break;
                } else |_| {}
            } else |_| {}
            signed_step = if (signed_step > 0)
                -signed_step
            else
                -2 * signed_step;
        }
        if (!found_probe) return null;
    }
    var right_hand_side: [reduced_active_coordinate_count]f64 = undefined;
    for (0..reduced_active_coordinate_count) |row|
        right_hand_side[row] = -base_residual[row] /
            reducedActiveResidualScale(current, row, density, options);
    var solved_matrix = matrix;
    var delta = right_hand_side;
    if (!numerics.solveDenseLinearSystem(
        &solved_matrix,
        &delta,
        reduced_active_coordinate_count,
    )) {
        var normal_matrix: [
            reduced_active_coordinate_count *
                reduced_active_coordinate_count
        ]f64 = undefined;
        var normal_right_hand_side: [reduced_active_coordinate_count]f64 = undefined;
        if (!solveDampedLeastSquares(
            &matrix,
            &right_hand_side,
            &delta,
            &normal_matrix,
            &normal_right_hand_side,
            reduced_active_coordinate_count,
        )) {
            return null;
        }
    }
    for (0..reduced_active_coordinate_count) |coordinate|
        delta[coordinate] *= coordinate_scales[coordinate];

    var fraction: f64 = 1;
    var line_search: u8 = 0;
    while (line_search < 53) : (line_search += 1) {
        var candidate_coordinates = coordinates;
        for (0..reduced_active_coordinate_count) |coordinate|
            candidate_coordinates[coordinate] +=
                fraction * delta[coordinate];
        // Progress toward the exact Gapon target is a bounded homotopy
        // coordinate, never an extrapolation beyond either physical state.
        if (candidate_coordinates[3] >= 0 and
            candidate_coordinates[3] <= 1)
        {
            if (reducedActiveCandidate(
                current,
                exchange_target,
                candidate_coordinates,
                site_total,
                phosphorus_total,
                density,
            )) |candidate| {
                const candidate_changes =
                    changesAt(candidate, environment, evaluator) catch {
                        fraction *= 0.5;
                        continue;
                    };
                const candidate_norm =
                    try scaledNorm(candidate, candidate_changes, options);
                if (meaningfullyImproves(current_norm, candidate_norm))
                    return candidate;
            } else |_| {}
        }
        fraction *= 0.5;
    }
    return null;
}

fn reducedActiveCandidate(
    current: Cell,
    exchange_target: Cell,
    coordinates: [reduced_active_coordinate_count]f64,
    site_total: f64,
    phosphorus_total: f64,
    density: f64,
) !Cell {
    var candidate = try interpolateCell(
        current,
        exchange_target,
        coordinates[3],
    );
    const cation_totals = phosphateCationTotals(candidate);
    candidate.hpo4_mol_p_per_m3 = coordinates[0];
    candidate.phosphate_minerals.aluminum_phosphate_mol_per_m3 =
        coordinates[1];
    candidate.phosphate_minerals.iron_phosphate_mol_per_m3 =
        coordinates[2];
    try closePhosphateConservation(
        &candidate,
        site_total,
        phosphorus_total,
        cation_totals,
        density,
    );
    return candidate;
}

fn reducedActiveResidual(
    cell: Cell,
    evaluator: Evaluator,
    density: f64,
    residual: *[reduced_active_coordinate_count]f64,
) !void {
    const extents = try evaluator.evaluate(evaluator.context, cell);
    residual.* = .{
        extents.h2po4_association_mol_p_per_m3,
        extents.phosphate_minerals.aluminum_phosphate_mol_per_m3,
        extents.phosphate_minerals.iron_phosphate_mol_per_m3,
        density * extents.exchange.calcium_mol_per_Mg,
    };
}

fn reducedActiveResidualScale(
    cell: Cell,
    row: usize,
    density: f64,
    options: Options,
) f64 {
    const inventory = switch (row) {
        0 => @min(
            cell.hpo4_mol_p_per_m3,
            cell.h2po4_mol_p_per_m3,
        ),
        1 => @min(
            cell.aluminum_mol_per_m3,
            cell.phosphate_minerals.aluminum_phosphate_mol_per_m3,
        ),
        2 => @min(
            cell.iron_mol_per_m3,
            cell.phosphate_minerals.iron_phosphate_mol_per_m3,
        ),
        3 => @min(
            cell.calcium_mol_per_m3,
            density * cell.exchange.calcium_mol_per_Mg,
        ),
        else => unreachable,
    };
    return options.absolute_tolerance +
        options.relative_tolerance * @max(1.0, @abs(inventory));
}

fn conservativePhosphateNewton(
    allocator: std.mem.Allocator,
    current: Cell,
    environment: Environment,
    evaluator: Evaluator,
    options: Options,
    current_norm: f64,
) !?Cell {
    const matrix = try allocator.alloc(f64, phosphate_coordinate_count * phosphate_coordinate_count);
    defer allocator.free(matrix);
    const right_hand_side = try allocator.alloc(f64, phosphate_coordinate_count);
    defer allocator.free(right_hand_side);
    const base_residual = try allocator.alloc(f64, phosphate_coordinate_count);
    defer allocator.free(base_residual);
    const probe_residual = try allocator.alloc(f64, phosphate_coordinate_count);
    defer allocator.free(probe_residual);
    const delta = try allocator.alloc(f64, phosphate_coordinate_count);
    defer allocator.free(delta);
    const normal_matrix = try allocator.alloc(
        f64,
        phosphate_coordinate_count * phosphate_coordinate_count,
    );
    defer allocator.free(normal_matrix);
    const normal_right_hand_side =
        try allocator.alloc(f64, phosphate_coordinate_count);
    defer allocator.free(normal_right_hand_side);

    const site_total = phosphateSiteTotal(current);
    const phosphorus_total = phosphateTotal(current, environment.litter_mass_per_water_volume_Mg_per_m3);
    const density = environment.litter_mass_per_water_volume_Mg_per_m3;
    const coupled_totals = coupledPhosphateExchangeTotals(current, density);
    const base_changes = try changesAt(current, environment, evaluator);
    const frozen_exchange_coordinates =
        frozenExchangeCoordinates(current, base_changes, options);
    const base_exchange_target = if (evaluator.equilibrate_cation_exchange) |equilibrate|
        equilibrate(evaluator.context, current) catch null
    else
        null;
    packCoupledPhosphateResidual(
        current,
        current,
        base_changes,
        base_exchange_target,
        frozen_exchange_coordinates,
        base_residual,
    );
    try replacePhosphateMineralResiduals(
        evaluator,
        current,
        base_residual,
    );
    const base_coupled_norm =
        scaledCoordinateResidualNorm(current, base_residual, options);
    for (0..phosphate_coordinate_count) |column| {
        if (frozen_exchange_coordinates[column]) {
            for (0..phosphate_coordinate_count) |row|
                matrix[row * phosphate_coordinate_count + column] = 0;
            matrix[column * phosphate_coordinate_count + column] = 1;
            continue;
        }
        const variable_scale =
            @max(1.0, @abs(phosphateCoordinate(current, column)));
        const step = std.math.cbrt(std.math.floatEps(f64)) *
            @max(1.0, @abs(phosphateCoordinate(current, column)));
        var selected_step: f64 = 0;
        var magnitude = step;
        var search: u8 = 0;
        while (search < 40) : (search += 1) {
            const signed_step = if (search % 2 == 0)
                magnitude
            else
                -magnitude;
            var trial = current;
            setPhosphateCoordinate(&trial, column, phosphateCoordinate(current, column) + signed_step);
            if (closeCoupledPhosphateExchangeConservation(
                &trial,
                site_total,
                phosphorus_total,
                coupled_totals,
                density,
            )) |_| {
                if (changesAt(trial, environment, evaluator)) |probe_changes| {
                    const probe_exchange_target = if (evaluator.equilibrate_cation_exchange) |equilibrate|
                        equilibrate(evaluator.context, trial) catch null
                    else
                        null;
                    packCoupledPhosphateResidual(
                        current,
                        trial,
                        probe_changes,
                        probe_exchange_target,
                        frozen_exchange_coordinates,
                        probe_residual,
                    );
                    try replacePhosphateMineralResiduals(
                        evaluator,
                        trial,
                        probe_residual,
                    );
                    if (residualProbeIsInformative(
                        base_residual,
                        probe_residual,
                    )) {
                        selected_step = signed_step;
                        break;
                    }
                } else |_| {}
            } else |_| {}
            if (search % 2 == 1) magnitude *= 2;
        }
        if (selected_step == 0) return null;
        for (0..phosphate_coordinate_count) |row| {
            const scale = options.absolute_tolerance +
                options.relative_tolerance * @max(1.0, @abs(phosphateCoordinate(current, row)));
            matrix[row * phosphate_coordinate_count + column] =
                (probe_residual[row] - base_residual[row]) /
                selected_step * variable_scale / scale;
        }
    }
    for (0..phosphate_coordinate_count) |row| {
        const scale = options.absolute_tolerance +
            options.relative_tolerance * @max(1.0, @abs(phosphateCoordinate(current, row)));
        right_hand_side[row] = -base_residual[row] / scale;
    }
    @memcpy(normal_matrix, matrix);
    @memcpy(delta, right_hand_side);
    const direct_solution = numerics.solveDenseLinearSystem(
        normal_matrix,
        delta,
        phosphate_coordinate_count,
    );
    if (!direct_solution and !solveDampedLeastSquares(
        matrix,
        right_hand_side,
        delta,
        normal_matrix,
        normal_right_hand_side,
        phosphate_coordinate_count,
    )) return null;
    for (0..phosphate_coordinate_count) |coordinate|
        delta[coordinate] *=
            @max(1.0, @abs(phosphateCoordinate(current, coordinate)));

    var line_fraction: f64 = 1;
    var line_search: u8 = 0;
    while (line_search < 24) : (line_search += 1) {
        var candidate = current;
        for (0..phosphate_coordinate_count) |coordinate|
            setPhosphateCoordinate(&candidate, coordinate, phosphateCoordinate(current, coordinate) + line_fraction * delta[coordinate]);
        if (closeCoupledPhosphateExchangeConservation(
            &candidate,
            site_total,
            phosphorus_total,
            coupled_totals,
            density,
        )) |_| {
            const candidate_changes = changesAt(
                candidate,
                environment,
                evaluator,
            ) catch {
                line_fraction *= 0.5;
                continue;
            };
            const candidate_norm = try scaledNorm(candidate, candidate_changes, options);
            const candidate_exchange_target = if (evaluator.equilibrate_cation_exchange) |equilibrate|
                equilibrate(evaluator.context, candidate) catch null
            else
                null;
            packCoupledPhosphateResidual(
                current,
                candidate,
                candidate_changes,
                candidate_exchange_target,
                frozen_exchange_coordinates,
                probe_residual,
            );
            try replacePhosphateMineralResiduals(
                evaluator,
                candidate,
                probe_residual,
            );
            const candidate_coupled_norm =
                scaledCoordinateResidualNorm(current, probe_residual, options);
            // The full-cell infinity norm can be controlled by a capped
            // reaction outside this conservative coordinate block. Accept
            // a strict decrease of the coupled block while preventing a
            // material increase of that global merit; subsequent hybrid
            // iterations then resolve the remaining reaction.
            if (meaningfullyImproves(current_norm, candidate_norm) or
                (meaningfullyImproves(
                    base_coupled_norm,
                    candidate_coupled_norm,
                ) and
                    candidate_norm <= current_norm * (1 + 1e-10)))
                return candidate;
        } else |_| {}
        line_fraction *= 0.5;
    }
    return null;
}

fn meaningfullyImproves(current_norm: f64, candidate_norm: f64) bool {
    return candidate_norm <= 1 or
        candidate_norm <= current_norm * (1 - 1e-6);
}

fn replacePhosphateMineralResiduals(
    evaluator: Evaluator,
    state_value: Cell,
    residual: []f64,
) !void {
    const calculate = evaluator.phosphate_mineral_equilibrium_residuals orelse
        return;
    const mineral = try calculate(evaluator.context, state_value);
    residual[5] = mineral.aluminum_phosphate_mol_per_m3;
    residual[6] = mineral.iron_phosphate_mol_per_m3;
    residual[7] = mineral.dicalcium_phosphate_mol_per_m3;
    residual[8] = mineral.hydroxyapatite_mol_per_m3;
    residual[9] = mineral.monocalcium_phosphate_mol_per_m3;
}

fn scaledCoordinateResidualNorm(
    reference: Cell,
    residual: []const f64,
    options: Options,
) f64 {
    var result: f64 = 0;
    for (residual, 0..) |value, coordinate| {
        const scale = options.absolute_tolerance +
            options.relative_tolerance *
                @max(1.0, @abs(phosphateCoordinate(reference, coordinate)));
        result = @max(result, @abs(value) / scale);
    }
    return result;
}

fn solveDampedLeastSquares(
    jacobian: []const f64,
    right_hand_side: []const f64,
    solution: []f64,
    normal_matrix: []f64,
    normal_right_hand_side: []f64,
    dimension: usize,
) bool {
    var maximum_diagonal: f64 = 0;
    for (0..dimension) |column| {
        var projected_rhs: f64 = 0;
        for (0..dimension) |row|
            projected_rhs += jacobian[row * dimension + column] *
                right_hand_side[row];
        normal_right_hand_side[column] = projected_rhs;
        for (0..dimension) |other_column| {
            var product: f64 = 0;
            for (0..dimension) |row|
                product += jacobian[row * dimension + column] *
                    jacobian[row * dimension + other_column];
            normal_matrix[column * dimension + other_column] = product;
        }
        maximum_diagonal = @max(
            maximum_diagonal,
            normal_matrix[column * dimension + column],
        );
    }
    if (!std.math.isFinite(maximum_diagonal) or maximum_diagonal <= 0)
        return false;

    // Marquardt diagonal scaling keeps weak but physically meaningful
    // coordinates from being erased by a single very stiff activity
    // derivative. A scalar `lambda * max(diagonal)` previously collapsed
    // the aluminum-phosphate step by many orders of magnitude.
    var damping_fraction = 128 * std.math.floatEps(f64);
    var attempt: u8 = 0;
    while (attempt < 10) : (attempt += 1) {
        for (0..dimension) |row| {
            for (0..dimension) |column| {
                var product: f64 = 0;
                for (0..dimension) |source_row|
                    product += jacobian[source_row * dimension + row] *
                        jacobian[source_row * dimension + column];
                normal_matrix[row * dimension + column] = product;
            }
            normal_matrix[row * dimension + row] +=
                damping_fraction *
                @max(normal_matrix[row * dimension + row], 1e-30);
            var projected_rhs: f64 = 0;
            for (0..dimension) |source_row|
                projected_rhs += jacobian[source_row * dimension + row] *
                    right_hand_side[source_row];
            normal_right_hand_side[row] = projected_rhs;
        }
        @memcpy(solution, normal_right_hand_side);
        if (numerics.solveDenseLinearSystem(
            normal_matrix,
            solution,
            dimension,
        )) return true;
        damping_fraction *= 100;
    }
    return false;
}

fn residualProbeIsInformative(
    base_residual: []const f64,
    probe_residual: []const f64,
) bool {
    for (base_residual, probe_residual) |base, probe| {
        const threshold = 128 * std.math.floatEps(f64) *
            @max(1.0, @abs(base), @abs(probe));
        if (@abs(probe - base) > threshold) return true;
    }
    return false;
}

fn phosphateSiteTotal(cell: Cell) f64 {
    const p = cell.phosphate_surface;
    return p.deprotonated_site_mol_per_Mg + p.hydroxyl_site_mol_per_Mg +
        p.protonated_site_mol_per_Mg + p.adsorbed_hpo4_mol_p_per_Mg +
        p.adsorbed_h2po4_mol_p_per_Mg;
}

fn phosphateTotal(cell: Cell, density: f64) f64 {
    const p = cell.phosphate_minerals;
    return cell.hpo4_mol_p_per_m3 + cell.h2po4_mol_p_per_m3 +
        density * (cell.phosphate_surface.adsorbed_hpo4_mol_p_per_Mg +
            cell.phosphate_surface.adsorbed_h2po4_mol_p_per_Mg) +
        p.aluminum_phosphate_mol_per_m3 + p.iron_phosphate_mol_per_m3 +
        p.dicalcium_phosphate_mol_per_m3 + 3 * p.hydroxyapatite_mol_per_m3 +
        2 * p.monocalcium_phosphate_mol_per_m3;
}

const PhosphateCationTotals = struct {
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
};

fn phosphateCationTotals(cell: Cell) PhosphateCationTotals {
    const minerals = cell.phosphate_minerals;
    return .{
        .aluminum_mol_per_m3 = cell.aluminum_mol_per_m3 +
            minerals.aluminum_phosphate_mol_per_m3,
        .iron_mol_per_m3 = cell.iron_mol_per_m3 +
            minerals.iron_phosphate_mol_per_m3,
        .calcium_mol_per_m3 = cell.calcium_mol_per_m3 +
            minerals.dicalcium_phosphate_mol_per_m3 +
            5 * minerals.hydroxyapatite_mol_per_m3 +
            minerals.monocalcium_phosphate_mol_per_m3,
    };
}

fn closePhosphateConservation(
    cell: *Cell,
    site_total: f64,
    phosphorus_total: f64,
    cation_totals: PhosphateCationTotals,
    density: f64,
) !void {
    const surface = &cell.phosphate_surface;
    surface.adsorbed_h2po4_mol_p_per_Mg = site_total -
        surface.deprotonated_site_mol_per_Mg -
        surface.hydroxyl_site_mol_per_Mg -
        surface.protonated_site_mol_per_Mg -
        surface.adsorbed_hpo4_mol_p_per_Mg;
    const minerals = cell.phosphate_minerals;
    cell.h2po4_mol_p_per_m3 = phosphorus_total -
        cell.hpo4_mol_p_per_m3 -
        density * (surface.adsorbed_hpo4_mol_p_per_Mg + surface.adsorbed_h2po4_mol_p_per_Mg) -
        minerals.aluminum_phosphate_mol_per_m3 -
        minerals.iron_phosphate_mol_per_m3 -
        minerals.dicalcium_phosphate_mol_per_m3 -
        3 * minerals.hydroxyapatite_mol_per_m3 -
        2 * minerals.monocalcium_phosphate_mol_per_m3;
    cell.aluminum_mol_per_m3 = cation_totals.aluminum_mol_per_m3 -
        minerals.aluminum_phosphate_mol_per_m3;
    cell.iron_mol_per_m3 = cation_totals.iron_mol_per_m3 -
        minerals.iron_phosphate_mol_per_m3;
    cell.calcium_mol_per_m3 = cation_totals.calcium_mol_per_m3 -
        minerals.dicalcium_phosphate_mol_per_m3 -
        5 * minerals.hydroxyapatite_mol_per_m3 -
        minerals.monocalcium_phosphate_mol_per_m3;
    normalizeRoundoffNonnegative(Cell, cell);
    try validateCell(cell.*);
}

const CoupledPhosphateExchangeTotals = struct {
    ammoniacal_nitrogen_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    exchange_charge_mol_per_Mg: f64,
};

fn coupledPhosphateExchangeTotals(
    cell: Cell,
    density: f64,
) CoupledPhosphateExchangeTotals {
    const minerals = cell.phosphate_minerals;
    const exchange = cell.exchange;
    return .{
        .ammoniacal_nitrogen_mol_per_m3 = cell.ammonia_mol_per_m3 + cell.ammonium_mol_per_m3 +
            density * exchange.ammonium_mol_per_Mg,
        .aluminum_mol_per_m3 = cell.aluminum_mol_per_m3 +
            minerals.aluminum_phosphate_mol_per_m3 +
            density * exchange.aluminum_mol_per_Mg,
        .iron_mol_per_m3 = cell.iron_mol_per_m3 +
            minerals.iron_phosphate_mol_per_m3 +
            density * exchange.iron_mol_per_Mg,
        .calcium_mol_per_m3 = cell.calcium_mol_per_m3 +
            minerals.dicalcium_phosphate_mol_per_m3 +
            5 * minerals.hydroxyapatite_mol_per_m3 +
            minerals.monocalcium_phosphate_mol_per_m3 +
            density * exchange.calcium_mol_per_Mg,
        .magnesium_mol_per_m3 = cell.magnesium_mol_per_m3 +
            density * exchange.magnesium_mol_per_Mg,
        .sodium_mol_per_m3 = cell.sodium_mol_per_m3 +
            density * exchange.sodium_mol_per_Mg,
        .potassium_mol_per_m3 = cell.potassium_mol_per_m3 +
            density * exchange.potassium_mol_per_Mg,
        .exchange_charge_mol_per_Mg = exchange.ammonium_mol_per_Mg +
            exchange.hydrogen_mol_per_Mg +
            3 * exchange.aluminum_mol_per_Mg +
            3 * exchange.iron_mol_per_Mg +
            2 * exchange.calcium_mol_per_Mg +
            2 * exchange.magnesium_mol_per_Mg +
            exchange.sodium_mol_per_Mg +
            exchange.potassium_mol_per_Mg,
    };
}

fn closeCoupledPhosphateExchangeConservation(
    cell: *Cell,
    site_total: f64,
    phosphorus_total: f64,
    totals: CoupledPhosphateExchangeTotals,
    density: f64,
) !void {
    try closePhosphateConservation(
        cell,
        site_total,
        phosphorus_total,
        phosphateCationTotals(cell.*),
        density,
    );
    const exchange = &cell.exchange;
    exchange.calcium_mol_per_Mg =
        0.5 * (totals.exchange_charge_mol_per_Mg -
            exchange.ammonium_mol_per_Mg -
            exchange.hydrogen_mol_per_Mg -
            3 * exchange.aluminum_mol_per_Mg -
            3 * exchange.iron_mol_per_Mg -
            2 * exchange.magnesium_mol_per_Mg -
            exchange.sodium_mol_per_Mg -
            exchange.potassium_mol_per_Mg);

    const minerals = cell.phosphate_minerals;
    cell.ammonium_mol_per_m3 =
        totals.ammoniacal_nitrogen_mol_per_m3 -
        cell.ammonia_mol_per_m3 -
        density * exchange.ammonium_mol_per_Mg;
    cell.aluminum_mol_per_m3 = totals.aluminum_mol_per_m3 -
        minerals.aluminum_phosphate_mol_per_m3 -
        density * exchange.aluminum_mol_per_Mg;
    cell.iron_mol_per_m3 = totals.iron_mol_per_m3 -
        minerals.iron_phosphate_mol_per_m3 -
        density * exchange.iron_mol_per_Mg;
    cell.calcium_mol_per_m3 = totals.calcium_mol_per_m3 -
        minerals.dicalcium_phosphate_mol_per_m3 -
        5 * minerals.hydroxyapatite_mol_per_m3 -
        minerals.monocalcium_phosphate_mol_per_m3 -
        density * exchange.calcium_mol_per_Mg;
    cell.magnesium_mol_per_m3 = totals.magnesium_mol_per_m3 -
        density * exchange.magnesium_mol_per_Mg;
    cell.sodium_mol_per_m3 = totals.sodium_mol_per_m3 -
        density * exchange.sodium_mol_per_Mg;
    cell.potassium_mol_per_m3 = totals.potassium_mol_per_m3 -
        density * exchange.potassium_mol_per_Mg;
    normalizeRoundoffNonnegative(Cell, cell);
    try validateCell(cell.*);
}

fn phosphateCoordinate(cell: Cell, index: usize) f64 {
    return switch (index) {
        0 => cell.hpo4_mol_p_per_m3,
        1 => cell.phosphate_surface.deprotonated_site_mol_per_Mg,
        2 => cell.phosphate_surface.hydroxyl_site_mol_per_Mg,
        3 => cell.phosphate_surface.protonated_site_mol_per_Mg,
        4 => cell.phosphate_surface.adsorbed_hpo4_mol_p_per_Mg,
        5 => cell.phosphate_minerals.aluminum_phosphate_mol_per_m3,
        6 => cell.phosphate_minerals.iron_phosphate_mol_per_m3,
        7 => cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3,
        8 => cell.phosphate_minerals.hydroxyapatite_mol_per_m3,
        9 => cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3,
        10 => cell.exchange.ammonium_mol_per_Mg,
        11 => cell.exchange.hydrogen_mol_per_Mg,
        12 => cell.exchange.aluminum_mol_per_Mg,
        13 => cell.exchange.iron_mol_per_Mg,
        14 => cell.exchange.magnesium_mol_per_Mg,
        15 => cell.exchange.sodium_mol_per_Mg,
        16 => cell.exchange.potassium_mol_per_Mg,
        17 => cell.ammonia_mol_per_m3,
        else => unreachable,
    };
}

fn setPhosphateCoordinate(cell: *Cell, index: usize, value: f64) void {
    switch (index) {
        0 => cell.hpo4_mol_p_per_m3 = value,
        1 => cell.phosphate_surface.deprotonated_site_mol_per_Mg = value,
        2 => cell.phosphate_surface.hydroxyl_site_mol_per_Mg = value,
        3 => cell.phosphate_surface.protonated_site_mol_per_Mg = value,
        4 => cell.phosphate_surface.adsorbed_hpo4_mol_p_per_Mg = value,
        5 => cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 = value,
        6 => cell.phosphate_minerals.iron_phosphate_mol_per_m3 = value,
        7 => cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3 = value,
        8 => cell.phosphate_minerals.hydroxyapatite_mol_per_m3 = value,
        9 => cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 = value,
        10 => cell.exchange.ammonium_mol_per_Mg = value,
        11 => cell.exchange.hydrogen_mol_per_Mg = value,
        12 => cell.exchange.aluminum_mol_per_Mg = value,
        13 => cell.exchange.iron_mol_per_Mg = value,
        14 => cell.exchange.magnesium_mol_per_Mg = value,
        15 => cell.exchange.sodium_mol_per_Mg = value,
        16 => cell.exchange.potassium_mol_per_Mg = value,
        17 => cell.ammonia_mol_per_m3 = value,
        else => unreachable,
    }
}

fn packCoupledPhosphateResidual(
    reference: Cell,
    state_value: Cell,
    changes: Cell,
    equilibrium_target: ?Cell,
    frozen_exchange_coordinates: [phosphate_coordinate_count]bool,
    output: []f64,
) void {
    for (0..phosphate_coordinate_count) |index| {
        output[index] = if (frozen_exchange_coordinates[index])
            phosphateCoordinate(state_value, index) -
                phosphateCoordinate(reference, index)
        else if (index >= 10 and index <= 16 and equilibrium_target != null)
            phosphateCoordinate(equilibrium_target.?, index) -
                phosphateCoordinate(state_value, index)
        else
            phosphateCoordinate(changes, index);
    }
}

fn frozenExchangeCoordinates(
    state_value: Cell,
    changes: Cell,
    options: Options,
) [phosphate_coordinate_count]bool {
    var result = [_]bool{false} ** phosphate_coordinate_count;
    for (10..17) |index| {
        const scale = options.absolute_tolerance +
            options.relative_tolerance *
                @max(1.0, @abs(phosphateCoordinate(state_value, index)));
        result[index] =
            @abs(phosphateCoordinate(changes, index)) <= scale;
    }
    // NH3/NH4 association has its own exact conservative scalar solve.
    // Leaving NH3 in the rank-deficient phosphate normal equations lets its
    // much larger concentration scale dominate the mineral direction.
    result[17] = true;
    return result;
}

fn logAdmissibleLimits(comptime T: type, comptime prefix: []const u8, current: T, changes: T, limiting_fraction: f64) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const change = @field(changes, field.name);
            if (change < 0) {
                const fraction = @max(0, @field(current, field.name)) / -change;
                if (fraction <= limiting_fraction * (1 + 1e-10))
                    std.log.warn("litter chemistry active-set limit: field={s}{s} value={e} change={e} fraction={e}", .{ prefix, field.name, @field(current, field.name), change, fraction });
            }
        },
        .@"struct" => logAdmissibleLimits(field.type, prefix ++ field.name ++ ".", @field(current, field.name), @field(changes, field.name), limiting_fraction),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn maximumAdmissibleFraction(comptime T: type, current: T, changes: T) f64 {
    var result = std.math.inf(f64);
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const value = @field(current, field.name);
            const change = @field(changes, field.name);
            // Underflow-scale dissolution of an already extinct phase is
            // numerical noise, not a physical active-set boundary. It can
            // otherwise cap a simultaneous mineral plateau at one iteration
            // even though scaling that change by any useful finite step still
            // lies inside the roundoff normalization allowance.
            const roundoff_change = std.math.floatEps(f64) *
                @max(1.0, @abs(value));
            if (change < -roundoff_change)
                result = @min(result, @max(0, value) / -change);
        },
        .@"struct" => result = @min(result, maximumAdmissibleFraction(field.type, @field(current, field.name), @field(changes, field.name))),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
    return result;
}

fn logUnconvergedFields(comptime T: type, comptime prefix: []const u8, current: T, changes: T, options: Options) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const value = @field(current, field.name);
            const change = @field(changes, field.name);
            const scaled = @abs(change) / (options.absolute_tolerance + options.relative_tolerance * @max(1.0, @abs(value)));
            if (scaled > 1)
                std.log.warn("litter chemistry residual: field={s}{s} value={e} change={e} scaled={e}", .{ prefix, field.name, value, change, scaled });
        },
        .@"struct" => logUnconvergedFields(field.type, prefix ++ field.name ++ ".", @field(current, field.name), @field(changes, field.name), options),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn valuesEqual(comptime T: type, left: T, right: T) bool {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => if (@field(left, field.name) != @field(right, field.name)) return false,
        .@"struct" => if (!valuesEqual(field.type, @field(left, field.name), @field(right, field.name))) return false,
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
    return true;
}

fn boundedPicard(current: Cell, changes: Cell, requested_fraction: f64) !Cell {
    var fraction = requested_fraction;
    var attempts: u8 = 0;
    while (attempts < 53) : (attempts += 1) {
        if (applyFraction(current, changes, fraction)) |candidate| return candidate else |err| switch (err) {
            error.NegativeLitterChemistryState => fraction *= 0.5,
            else => return err,
        }
        if (fraction <= std.math.floatEps(f64)) break;
    }
    logInadmissibleDirections(Cell, "", current, changes);
    return error.NoPhysicallyAdmissibleLitterPicardStep;
}

fn logInadmissibleDirections(comptime T: type, comptime prefix: []const u8, current: T, changes: T) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const value = @field(current, field.name);
            const change = @field(changes, field.name);
            if (change < 0 and value <= std.math.floatEps(f64) * @max(1.0, @abs(change)))
                std.log.err("litter chemistry has no admissible direction: field={s}{s} value={e} change={e}", .{ prefix, field.name, value, change });
        },
        .@"struct" => logInadmissibleDirections(field.type, prefix ++ field.name ++ ".", @field(current, field.name), @field(changes, field.name)),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn changesAt(cell: Cell, environment: Environment, evaluator: Evaluator) !Cell {
    const extents = try evaluator.evaluate(evaluator.context, cell);
    return ledger.assemble(extents, environment.litter_mass_per_water_volume_Mg_per_m3, environment.dynamic_salts);
}

fn applyFraction(current: Cell, changes: Cell, fraction: f64) !Cell {
    var result = current;
    try addScaled(Cell, &result, changes, fraction);
    normalizeRoundoffNonnegative(Cell, &result);
    try validateCell(result);
    return result;
}

fn interpolateCell(current: Cell, target: Cell, fraction: f64) !Cell {
    var result = current;
    try interpolateValue(Cell, &result, target, fraction);
    normalizeRoundoffNonnegative(Cell, &result);
    try validateCell(result);
    return result;
}

fn interpolateValue(comptime T: type, destination: *T, target: T, fraction: f64) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const next = @field(destination.*, field.name) +
                fraction * (@field(target, field.name) -
                    @field(destination.*, field.name));
            if (!std.math.isFinite(next)) return error.NonFiniteLitterChemistryState;
            @field(destination.*, field.name) = next;
        },
        .@"struct" => try interpolateValue(
            field.type,
            &@field(destination.*, field.name),
            @field(target, field.name),
            fraction,
        ),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn normalizeRoundoffNonnegative(comptime T: type, value: *T) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const number = &@field(value.*, field.name);
            if (number.* < 0 and number.* >= -1e-12) number.* = 0;
        },
        .@"struct" => normalizeRoundoffNonnegative(field.type, &@field(value.*, field.name)),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn addScaled(comptime T: type, destination: *T, changes: T, fraction: f64) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const change = @field(changes, field.name) * fraction;
            if (!std.math.isFinite(change)) return error.NonFiniteLitterChemistryChange;
            @field(destination.*, field.name) += change;
        },
        .@"struct" => try addScaled(field.type, &@field(destination.*, field.name), @field(changes, field.name), fraction),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn validateCell(cell: Cell) !void {
    try validateValue(Cell, cell);
}

fn validateValue(comptime T: type, value: T) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const number = @field(value, field.name);
            if (!std.math.isFinite(number)) return error.NonFiniteLitterChemistryState;
            if (number < -1e-12) return error.NegativeLitterChemistryState;
        },
        .@"struct" => try validateValue(field.type, @field(value, field.name)),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn scaledNorm(state_value: Cell, changes: Cell, options: Options) !f64 {
    var maximum: f64 = 0;
    try accumulateScaledNorm(Cell, state_value, changes, options, &maximum);
    return maximum;
}

fn accumulateScaledNorm(comptime T: type, state_value: T, changes: T, options: Options, maximum: *f64) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const value = @field(state_value, field.name);
            const change = @field(changes, field.name);
            if (!std.math.isFinite(value) or !std.math.isFinite(change)) return error.NonFiniteLitterChemistryState;
            maximum.* = @max(maximum.*, @abs(change) / (options.absolute_tolerance + options.relative_tolerance * @max(1.0, @abs(value))));
        },
        .@"struct" => try accumulateScaledNorm(field.type, @field(state_value, field.name), @field(changes, field.name), options, maximum),
        else => unreachable,
    };
}

fn directionalProducts(comptime T: type, base: T, probe: T, probe_fraction: f64, numerator: *f64, denominator: *f64) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => {
            const derivative = (@field(probe, field.name) - @field(base, field.name)) / probe_fraction;
            numerator.* += @field(base, field.name) * derivative;
            denominator.* += derivative * derivative;
        },
        .@"struct" => directionalProducts(field.type, @field(base, field.name), @field(probe, field.name), probe_fraction, numerator, denominator),
        else => unreachable,
    };
}

fn maximumDifference(comptime T: type, a: T, b: T) f64 {
    var maximum: f64 = 0;
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => maximum = @max(maximum, @abs(@field(a, field.name) - @field(b, field.name))),
        .@"struct" => maximum = @max(maximum, maximumDifference(field.type, @field(a, field.name), @field(b, field.name))),
        else => unreachable,
    };
    return maximum;
}

fn maximumMagnitude(comptime T: type, value: T) f64 {
    var maximum: f64 = 0;
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => maximum = @max(maximum, @abs(@field(value, field.name))),
        .@"struct" => maximum = @max(maximum, maximumMagnitude(field.type, @field(value, field.name))),
        else => unreachable,
    };
    return maximum;
}

fn zeroValue(comptime T: type, value: *T) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => @field(value.*, field.name) = 0,
        .@"struct" => zeroValue(field.type, &@field(value.*, field.name)),
        else => @compileError("litter chemistry contains a non-numeric field"),
    };
}

fn validateOptions(options: Options) !void {
    if (!std.math.isFinite(options.absolute_tolerance) or options.absolute_tolerance <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or options.max_iterations == 0) return error.InvalidLitterChemistrySolverOptions;
}

const TestContext = struct { rate_fraction: f64 };

fn testEvaluator(raw: *const anyopaque, cell: Cell) !ledger.ReactionExtents {
    const context: *const TestContext = @ptrCast(@alignCast(raw));
    var extents: ledger.ReactionExtents = undefined;
    zeroValue(ledger.ReactionExtents, &extents);
    extents.ammonium_association_mol_per_m3 = context.rate_fraction * (cell.ammonia_mol_per_m3 - cell.ammonium_mol_per_m3);
    extents.external_hydrogen_mol_per_m3 =
        extents.ammonium_association_mol_per_m3;
    return extents;
}

fn cappedAssociationPlateauEvaluator(
    _: *const anyopaque,
    cell: Cell,
) !ledger.ReactionExtents {
    var extents = std.mem.zeroes(ledger.ReactionExtents);
    extents.ammonium_association_mol_per_m3 =
        if (cell.ammonia_mol_per_m3 > 1)
            0.25
        else
            0.5 *
                (cell.ammonia_mol_per_m3 -
                    cell.ammonium_mol_per_m3);
    return extents;
}

fn simultaneousPhosphateSinkEvaluator(
    _: *const anyopaque,
    _: Cell,
) !ledger.ReactionExtents {
    var extents = std.mem.zeroes(ledger.ReactionExtents);
    extents.phosphate_minerals.aluminum_phosphate_mol_per_m3 = 0.5;
    extents.phosphate_minerals.iron_phosphate_mol_per_m3 = 0.5;
    extents.phosphate_minerals.dicalcium_phosphate_mol_per_m3 = 0.5;
    return extents;
}

fn cationBoundPhosphateEvaluator(
    _: *const anyopaque,
    cell: Cell,
) !ledger.ReactionExtents {
    var extents = std.mem.zeroes(ledger.ReactionExtents);
    extents.phosphate_minerals.aluminum_phosphate_mol_per_m3 =
        0.2 * cell.aluminum_mol_per_m3;
    extents.phosphate_minerals.iron_phosphate_mol_per_m3 =
        0.2 * cell.iron_mol_per_m3;
    return extents;
}

fn supersaturatedAluminumAndIronPhosphate(
    _: *const anyopaque,
    _: Cell,
) !ledger.PhosphateMineralExtents {
    return .{
        .aluminum_phosphate_mol_per_m3 = 1,
        .iron_phosphate_mol_per_m3 = 1,
        .dicalcium_phosphate_mol_per_m3 = 0,
        .hydroxyapatite_mol_per_m3 = 0,
        .monocalcium_phosphate_mol_per_m3 = 0,
    };
}

test "fixed-pH cation-bound acceleration advances Al and Fe simultaneously" {
    var cell = std.mem.zeroes(Cell);
    cell.aluminum_mol_per_m3 = 2;
    cell.iron_mol_per_m3 = 3;
    cell.h2po4_mol_p_per_m3 = 10;
    cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 = 4;
    cell.phosphate_minerals.iron_phosphate_mol_per_m3 = 5;
    const context: u8 = 0;
    const evaluator = Evaluator{
        .context = &context,
        .evaluate = cationBoundPhosphateEvaluator,
        .phosphate_mineral_equilibrium_residuals = supersaturatedAluminumAndIronPhosphate,
    };
    const options = Options{};
    const changes = try changesAt(
        cell,
        .{
            .litter_mass_per_water_volume_Mg_per_m3 = 1,
            .dynamic_salts = false,
        },
        evaluator,
    );
    const norm = try scaledNorm(cell, changes, options);
    const candidate = (try accelerateFixedPhosphateCationBounds(
        cell,
        .{
            .litter_mass_per_water_volume_Mg_per_m3 = 1,
            .dynamic_salts = false,
        },
        evaluator,
        options,
        norm,
    )).?;
    try std.testing.expect(candidate.aluminum_mol_per_m3 < cell.aluminum_mol_per_m3);
    try std.testing.expect(candidate.iron_mol_per_m3 < cell.iron_mol_per_m3);
    try std.testing.expectApproxEqAbs(
        cell.aluminum_mol_per_m3 +
            cell.phosphate_minerals.aluminum_phosphate_mol_per_m3,
        candidate.aluminum_mol_per_m3 +
            candidate.phosphate_minerals.aluminum_phosphate_mol_per_m3,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        cell.iron_mol_per_m3 +
            cell.phosphate_minerals.iron_phosphate_mol_per_m3,
        candidate.iron_mol_per_m3 +
            candidate.phosphate_minerals.iron_phosphate_mol_per_m3,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        cell.h2po4_mol_p_per_m3 +
            cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 +
            cell.phosphate_minerals.iron_phosphate_mol_per_m3,
        candidate.h2po4_mol_p_per_m3 +
            candidate.phosphate_minerals.aluminum_phosphate_mol_per_m3 +
            candidate.phosphate_minerals.iron_phosphate_mol_per_m3,
        1e-14,
    );
}

test "runtime litter cell solve converges early and conserves nitrogen" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    state.cells[2].ammonia_mol_per_m3 = 2;
    const before = state.cells[2].ammonia_mol_per_m3 + state.cells[2].ammonium_mol_per_m3;
    const context = TestContext{ .rate_fraction = 0.25 };
    const result = try solveCell(&state, 2, .{ .litter_mass_per_water_volume_Mg_per_m3 = 1, .dynamic_salts = true }, .{ .context = &context, .evaluate = testEvaluator }, .{});
    try std.testing.expect(result.iterations < 60);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
    try std.testing.expectApproxEqAbs(before, state.cells[2].ammonia_mol_per_m3 + state.cells[2].ammonium_mol_per_m3, 1e-12);
    try std.testing.expectApproxEqAbs(state.cells[2].ammonia_mol_per_m3, state.cells[2].ammonium_mol_per_m3, 5e-8);
}

test "fixed-pH branch commits one hourly kinetic increment without MRXN multiplication" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].ammonia_mol_per_m3 = 2;
    const context = TestContext{ .rate_fraction = 0.25 };
    const result = try solveCell(
        &state,
        0,
        .{
            .litter_mass_per_water_volume_Mg_per_m3 = 1,
            .dynamic_salts = false,
        },
        .{ .context = &context, .evaluate = testEvaluator },
        .{ .max_iterations = 1000 },
    );
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(u16, 0), result.newton_raphson_steps);
    try std.testing.expectEqual(@as(u16, 0), result.picard_steps);
    try std.testing.expectEqual(@as(f64, 1.5), state.cells[0].ammonia_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0.5), state.cells[0].ammonium_mol_per_m3);
}

test "fixed-pH simultaneous sinks share one exact admissible substrate fraction" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].h2po4_mol_p_per_m3 = 1;
    const context: u8 = 0;
    const result = try solveCell(
        &state,
        0,
        .{
            .litter_mass_per_water_volume_Mg_per_m3 = 1,
            .dynamic_salts = false,
        },
        .{
            .context = &context,
            .evaluate = simultaneousPhosphateSinkEvaluator,
        },
        .{},
    );
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(f64, 0), state.cells[0].h2po4_mol_p_per_m3);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0 / 3.0),
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0 / 3.0),
        state.cells[0].phosphate_minerals.iron_phosphate_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0 / 3.0),
        state.cells[0].phosphate_minerals.dicalcium_phosphate_mol_per_m3,
        1e-15,
    );
}

test "transactional boundary lookahead crosses a neutral clipped plateau" {
    var current = std.mem.zeroes(Cell);
    current.ammonia_mol_per_m3 = 2;
    const environment = Environment{
        .litter_mass_per_water_volume_Mg_per_m3 = 1,
        .dynamic_salts = false,
    };
    const context: u8 = 0;
    const evaluator = Evaluator{
        .context = &context,
        .evaluate = cappedAssociationPlateauEvaluator,
    };
    const options = Options{};
    const changes = try changesAt(current, environment, evaluator);
    const current_norm = try scaledNorm(current, changes, options);
    const candidate = (try plateauBoundaryLookahead(
        current,
        changes,
        maximumAdmissibleFraction(Cell, current, changes),
        environment,
        evaluator,
        options,
        current_norm,
    )).?;
    const candidate_changes = try changesAt(candidate, environment, evaluator);
    try std.testing.expect(
        try scaledNorm(candidate, candidate_changes, options) < current_norm,
    );
    try std.testing.expectApproxEqAbs(
        current.ammonia_mol_per_m3 + current.ammonium_mol_per_m3,
        candidate.ammonia_mol_per_m3 + candidate.ammonium_mol_per_m3,
        1e-14,
    );
    try std.testing.expect(candidate.ammonia_mol_per_m3 > 0);
    try std.testing.expect(candidate.ammonium_mol_per_m3 > 0);
}

test "subnormal extinct-phase noise does not cap an admissible mineral step" {
    var current = std.mem.zeroes(Cell);
    current.h2po4_mol_p_per_m3 = 10;
    var changes = std.mem.zeroes(Cell);
    changes.h2po4_mol_p_per_m3 = -0.005;
    changes.phosphate_minerals.dicalcium_phosphate_mol_per_m3 = -3.0e-319;
    changes.phosphate_minerals.hydroxyapatite_mol_per_m3 = -9.0e-319;
    changes.phosphate_minerals.monocalcium_phosphate_mol_per_m3 = -6.0e-319;
    try std.testing.expectApproxEqAbs(
        @as(f64, 2000),
        maximumAdmissibleFraction(Cell, current, changes),
        1e-10,
    );
}

test "convergence reached on the last permitted iteration is committed" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].ammonia_mol_per_m3 = 2;
    const context = TestContext{ .rate_fraction = 0.25 };
    const result = try solveCell(&state, 0, .{ .litter_mass_per_water_volume_Mg_per_m3 = 1, .dynamic_salts = true }, .{ .context = &context, .evaluate = testEvaluator }, .{ .max_iterations = 1 });
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expect(result.maximum_scaled_residual <= 1);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.cells[0].ammonium_mol_per_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.cells[0].ammonia_mol_per_m3, 1e-12);
}

test "phosphate Newton coordinates close site phosphorus and cation conservation exactly" {
    var cell = std.mem.zeroes(Cell);
    cell.hpo4_mol_p_per_m3 = 2;
    cell.h2po4_mol_p_per_m3 = 3;
    cell.phosphate_surface = .{
        .deprotonated_site_mol_per_Mg = 0.1,
        .hydroxyl_site_mol_per_Mg = 0.2,
        .protonated_site_mol_per_Mg = 0.3,
        .adsorbed_hpo4_mol_p_per_Mg = 0.4,
        .adsorbed_h2po4_mol_p_per_Mg = 0.5,
    };
    cell.phosphate_minerals = .{
        .aluminum_phosphate_mol_per_m3 = 0.6,
        .iron_phosphate_mol_per_m3 = 0.7,
        .dicalcium_phosphate_mol_per_m3 = 0.8,
        .hydroxyapatite_mol_per_m3 = 0.9,
        .monocalcium_phosphate_mol_per_m3 = 1,
    };
    cell.aluminum_mol_per_m3 = 2;
    cell.iron_mol_per_m3 = 2;
    cell.calcium_mol_per_m3 = 8;
    const density: f64 = 2;
    const site_total = phosphateSiteTotal(cell);
    const phosphorus_total = phosphateTotal(cell, density);
    const cation_totals = phosphateCationTotals(cell);
    cell.hpo4_mol_p_per_m3 += 0.1;
    cell.phosphate_surface.deprotonated_site_mol_per_Mg += 0.01;
    cell.phosphate_surface.adsorbed_hpo4_mol_p_per_Mg -= 0.02;
    cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 += 0.03;
    cell.phosphate_minerals.hydroxyapatite_mol_per_m3 -= 0.04;
    try closePhosphateConservation(
        &cell,
        site_total,
        phosphorus_total,
        cation_totals,
        density,
    );
    try std.testing.expectApproxEqAbs(site_total, phosphateSiteTotal(cell), 1e-15);
    try std.testing.expectApproxEqAbs(phosphorus_total, phosphateTotal(cell, density), 1e-14);
    const cations_after = phosphateCationTotals(cell);
    try std.testing.expectApproxEqAbs(
        cation_totals.aluminum_mol_per_m3,
        cations_after.aluminum_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        cation_totals.iron_mol_per_m3,
        cations_after.iron_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        cation_totals.calcium_mol_per_m3,
        cations_after.calcium_mol_per_m3,
        1e-14,
    );
}

test "coupled phosphate exchange coordinates conserve elements and exchange charge" {
    var cell = std.mem.zeroes(Cell);
    cell.ammonia_mol_per_m3 = 3;
    cell.ammonium_mol_per_m3 = 5;
    cell.hpo4_mol_p_per_m3 = 4;
    cell.h2po4_mol_p_per_m3 = 6;
    cell.aluminum_mol_per_m3 = 7;
    cell.iron_mol_per_m3 = 8;
    cell.calcium_mol_per_m3 = 20;
    cell.magnesium_mol_per_m3 = 9;
    cell.sodium_mol_per_m3 = 10;
    cell.potassium_mol_per_m3 = 11;
    cell.phosphate_surface.deprotonated_site_mol_per_Mg = 0.1;
    cell.phosphate_surface.hydroxyl_site_mol_per_Mg = 0.2;
    cell.phosphate_surface.protonated_site_mol_per_Mg = 0.3;
    cell.phosphate_surface.adsorbed_hpo4_mol_p_per_Mg = 0.4;
    cell.phosphate_surface.adsorbed_h2po4_mol_p_per_Mg = 0.5;
    cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 = 0.6;
    cell.phosphate_minerals.iron_phosphate_mol_per_m3 = 0.7;
    cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3 = 0.8;
    cell.phosphate_minerals.hydroxyapatite_mol_per_m3 = 0.9;
    cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 = 1;
    cell.exchange = .{
        .ammonium_mol_per_Mg = 0.11,
        .hydrogen_mol_per_Mg = 0.12,
        .aluminum_mol_per_Mg = 0.13,
        .iron_mol_per_Mg = 0.14,
        .calcium_mol_per_Mg = 0.15,
        .magnesium_mol_per_Mg = 0.16,
        .sodium_mol_per_Mg = 0.17,
        .potassium_mol_per_Mg = 0.18,
    };
    const density: f64 = 2;
    const site_total = phosphateSiteTotal(cell);
    const phosphorus_total = phosphateTotal(cell, density);
    const totals = coupledPhosphateExchangeTotals(cell, density);

    cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 += 0.2;
    cell.phosphate_minerals.hydroxyapatite_mol_per_m3 -= 0.1;
    cell.exchange.ammonium_mol_per_Mg += 0.01;
    cell.exchange.hydrogen_mol_per_Mg -= 0.02;
    cell.exchange.aluminum_mol_per_Mg += 0.01;
    cell.exchange.magnesium_mol_per_Mg -= 0.03;
    cell.exchange.sodium_mol_per_Mg += 0.02;
    cell.ammonia_mol_per_m3 += 0.25;
    try closeCoupledPhosphateExchangeConservation(
        &cell,
        site_total,
        phosphorus_total,
        totals,
        density,
    );

    try std.testing.expectApproxEqAbs(
        phosphorus_total,
        phosphateTotal(cell, density),
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        site_total,
        phosphateSiteTotal(cell),
        1e-15,
    );
    const after = coupledPhosphateExchangeTotals(cell, density);
    inline for (@typeInfo(CoupledPhosphateExchangeTotals).@"struct".fields) |field|
        try std.testing.expectApproxEqAbs(
            @field(totals, field.name),
            @field(after, field.name),
            1e-13,
        );
}

test "damped least squares resolves a rank deficient coupled Jacobian" {
    const jacobian = [_]f64{
        1, 1,
        2, 2,
    };
    const right_hand_side = [_]f64{ 2, 4 };
    var solution = [_]f64{ 0, 0 };
    var normal_matrix = [_]f64{ 0, 0, 0, 0 };
    var normal_right_hand_side = [_]f64{ 0, 0 };
    try std.testing.expect(solveDampedLeastSquares(
        &jacobian,
        &right_hand_side,
        &solution,
        &normal_matrix,
        &normal_right_hand_side,
        2,
    ));
    try std.testing.expectApproxEqAbs(@as(f64, 1), solution[0], 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 1), solution[1], 1e-8);
}

test "Marquardt scaling preserves weak coordinates beside stiff chemistry" {
    const jacobian = [_]f64{
        1e12, 0,
        0,    1,
    };
    const right_hand_side = [_]f64{ 1e12, 1 };
    var solution = [_]f64{ 0, 0 };
    var normal_matrix = [_]f64{ 0, 0, 0, 0 };
    var normal_right_hand_side = [_]f64{ 0, 0 };
    try std.testing.expect(solveDampedLeastSquares(
        &jacobian,
        &right_hand_side,
        &solution,
        &normal_matrix,
        &normal_right_hand_side,
        2,
    ));
    try std.testing.expectApproxEqAbs(@as(f64, 1), solution[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1), solution[1], 1e-10);
}
