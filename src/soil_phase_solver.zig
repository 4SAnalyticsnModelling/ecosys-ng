const std = @import("std");
const builtin = @import("builtin");
const grid_module = @import("grid.zig");
const phase = @import("soil_water_phase_change.zig");
const numerics = @import("numerics.zig");
const retention = @import("soil_water_retention.zig");

pub const Properties = struct {
    matrix_bulk_volume_m3: []const f64,
    retention_curve: []const retention.ResolvedCurve,
    mualem_van_genuchten_parameters: []const retention.MualemVanGenuchtenParameters = &.{},
    macropore_mualem_van_genuchten_parameters: []const retention.MualemVanGenuchtenParameters = &.{},
    osmotic_potential_mpa: []const f64,
    saturation_water_potential_mpa: []const f64,
    heat_capacity_mj_per_k: []const f64,
    saturated_lateral_matrix_conductivity_m2_per_h_mpa: []const f64,
    face_area_m2: []const f64,
    macropore_spacing_m: []const f64,
    macropore_radius_m: []const f64,
    pore_exchange_enabled: []const bool,
    /// Liquid/ice equilibrium is owned exclusively by the immediately
    /// following conservative Dall'Amico heat-enthalpy residual.
    vapor: phase.VaporEquilibriumParameters,
    freeze_thaw: phase.FreezeThawParameters,
    gravitational_water_potential_mpa_per_m: f64 = 0.0098,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
    ice_heat_capacity_mj_per_m3_k: f64,
};

fn iceWaterEquivalentFactor(properties: Properties) f64 {
    _ = properties;
    return 1.0;
}

pub const Options = struct {
    max_iterations: u16,
    absolute_tolerance_m3: f64 = 1e-14,
    relative_tolerance: f64 = 1e-9,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
};

pub const Outputs = struct {
    latent_heat_mj: []f64,
    macropore_to_matrix_water_m3: []f64,
};

pub const Result = struct { iterations: u16, newton_raphson_steps: u16, picard_steps: u16, maximum_scaled_residual: f64 };

/// Simultaneously converges VOLW/VOLV/VOLI, VOLWH/VOLIH, and the cell-local
/// endpoint enthalpy temperature for every runtime cell. NPH is a maximum
/// iteration count, never a repeated model cycle.
pub fn solve(allocator: std.mem.Allocator, grid: *grid_module.GridState, properties: Properties, outputs: Outputs, options: Options) !Result {
    try validateInputs(grid, properties, outputs, options);
    const cells = grid.layer_count;
    const components_per_cell: usize = 6;
    const components = try std.math.mul(usize, cells, components_per_cell);
    const base = try allocator.alloc(f64, components);
    defer allocator.free(base);
    @memcpy(base[0 * cells .. 1 * cells], grid.matrix_liquid_water_m3);
    @memcpy(base[1 * cells .. 2 * cells], grid.water_vapor_volume_m3);
    @memcpy(base[2 * cells .. 3 * cells], grid.matrix_ice_water_m3);
    @memcpy(base[3 * cells .. 4 * cells], grid.macropore_liquid_water_m3);
    @memcpy(base[4 * cells .. 5 * cells], grid.macropore_ice_water_m3);
    @memcpy(base[5 * cells .. 6 * cells], grid.soil_temperature_k);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, components);
    defer allocator.free(residual);
    const target = try allocator.alloc(f64, components);
    defer allocator.free(target);
    const scratch = try allocator.alloc(f64, components);
    defer allocator.free(scratch);
    const trial_heat = try allocator.alloc(f64, cells);
    defer allocator.free(trial_heat);
    const trial_exchange = try allocator.alloc(f64, cells);
    defer allocator.free(trial_exchange);
    const probe = try allocator.alloc(f64, components);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, components);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, components);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, components);
    defer allocator.free(candidate_residual);
    const independent_components_per_cell: usize = 5;
    const block_jacobian = try allocator.alloc(f64, try std.math.mul(usize, cells, independent_components_per_cell * independent_components_per_cell));
    defer allocator.free(block_jacobian);
    const block_right_hand_side = try allocator.alloc(f64, try std.math.mul(usize, cells, independent_components_per_cell));
    defer allocator.free(block_right_hand_side);
    const block_delta = try allocator.alloc(f64, components);
    defer allocator.free(block_delta);
    const probe_step = try allocator.alloc(f64, cells);
    defer allocator.free(probe_step);
    const block_coordinate_fixed = try allocator.alloc(bool, try std.math.mul(usize, cells, independent_components_per_cell));
    defer allocator.free(block_coordinate_fixed);
    @memset(block_coordinate_fixed, false);
    var newton_steps: u16 = 0;
    var directional_newton_steps: u16 = 0;
    var reduced_block_newton_steps: u16 = 0;
    var diagonal_newton_steps: u16 = 0;
    var reduced_block_invalid_iterations: u16 = 0;
    var reduced_block_rejected_iterations: u16 = 0;
    var last_reduced_block_probe_error: ?anyerror = null;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(grid, properties, base, current, target, residual, scratch, trial_heat, trial_exchange);
        const norm = try scaledNorm(current, residual, options);
        if (norm <= 1) {
            try residualAt(grid, properties, base, current, target, residual, scratch, outputs.latent_heat_mj, outputs.macropore_to_matrix_water_m3);
            // `target` is the validated, pore-bounded phase transform. The
            // converged trial may differ by the allowed residual and can sit
            // microscopically above pore capacity; publishing the target
            // applies that final fixed-point correction without losing water.
            try commit(grid, target);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm };
        }
        var accepted_newton = false;
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(grid, properties, base, probe, target, probe_residual, scratch, trial_heat, trial_exchange)) |_| {
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
                        if (residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange)) |_| {
                            // This inexpensive global residual-direction probe
                            // is only a predictor. Requiring material decrease
                            // prevents tiny improvements from starving the
                            // exact cell-local Newton block near convergence.
                            // A bounded residual-direction Newton step often
                            // removes 50-75% of a freeze/thaw residual. That is
                            // substantial progress within NPH; requiring 80%
                            // removal rejected the useful step and forced six
                            // geometric Picard halvings in the Arctic case.
                            if (try scaledNorm(candidate, candidate_residual, options) <= 0.5 * norm) {
                                @memcpy(current, candidate);
                                newton_steps += 1;
                                directional_newton_steps += 1;
                                accepted_newton = true;
                            }
                        } else |_| {}
                    } else |_| {}
                }
            } else |_| {}
        } else |_| {}
        if (accepted_newton) continue;
        // Eliminate matrix liquid water from the Newton coordinates through
        // exact water-equivalent conservation. The independent coordinates
        // are vapor, matrix ice, macropore liquid, macropore ice, and endpoint
        // temperature; this removes the structural nullspace of the 6x6 block.
        var block_jacobian_valid = true;
        const independent_component = [_]usize{ 1, 2, 3, 4, 5 };
        @memset(block_coordinate_fixed, false);
        for (independent_component, 0..) |column_component, column| {
            @memcpy(probe, current);
            for (0..cells) |cell| {
                const index = column_component * cells + cell;
                const nominal = std.math.cbrt(std.math.floatEps(f64)) * @max(1e-6, @abs(current[index]));
                const density = iceWaterEquivalentFactor(properties);
                const conservation_coefficient = switch (column_component) {
                    1, 3 => 1.0,
                    2, 4 => density,
                    else => 0.0,
                };
                const matrix_air_m3 = @max(0, grid.matrix_pore_capacity_m3[cell] - current[cell] - current[2 * cells + cell]);
                const macro_air_m3 = @max(0, grid.macropore_pore_capacity_m3[cell] - current[3 * cells + cell] - current[4 * cells + cell]);
                var positive_limit = if (conservation_coefficient > 0) current[cell] / conservation_coefficient else std.math.inf(f64);
                if (column_component == 2) positive_limit = @min(positive_limit, matrix_air_m3 / (1.0 - density));
                if (column_component == 3) positive_limit = @min(positive_limit, macro_air_m3);
                if (column_component == 4) positive_limit = @min(positive_limit, macro_air_m3 / (1.0 - density));
                if (positive_limit >= nominal) {
                    probe_step[cell] = nominal;
                } else {
                    var negative_limit = current[index];
                    if (column_component == 1 or column_component == 3) negative_limit = @min(negative_limit, matrix_air_m3 / conservation_coefficient);
                    probe_step[cell] = -@min(nominal, 0.5 * negative_limit);
                    if (@abs(probe_step[cell]) <= 1e-20) probe_step[cell] = 0.5 * positive_limit;
                }
                if (!std.math.isFinite(probe_step[cell]) or @abs(probe_step[cell]) <= 1e-20) {
                    // This coordinate is fixed by simultaneous non-negativity
                    // and pore-capacity constraints in this cell. Leave the
                    // batched probe unchanged and install an identity column
                    // for this cell below instead of invalidating other cells.
                    probe_step[cell] = 0;
                    block_coordinate_fixed[cell * independent_components_per_cell + column] = true;
                    continue;
                }
                probe[index] += probe_step[cell];
                probe[cell] -= conservation_coefficient * probe_step[cell];
            }
            var probe_attempt: u8 = 0;
            while (true) : (probe_attempt += 1) {
                if (residualAt(grid, properties, base, probe, target, probe_residual, scratch, trial_heat, trial_exchange)) |_| break else |err| {
                    last_reduced_block_probe_error = err;
                    if (probe_attempt >= 15) {
                        block_jacobian_valid = false;
                        break;
                    }
                    @memcpy(probe, current);
                    const conservation_coefficient = switch (column_component) {
                        1, 3 => 1.0,
                        2, 4 => iceWaterEquivalentFactor(properties),
                        else => 0.0,
                    };
                    for (0..cells) |cell| {
                        probe_step[cell] *= 0.5;
                        probe[column_component * cells + cell] += probe_step[cell];
                        probe[cell] -= conservation_coefficient * probe_step[cell];
                    }
                }
            }
            if (!block_jacobian_valid) break;
            for (0..cells) |cell| for (independent_component, 0..) |row_component, row| {
                const residual_index = row_component * cells + cell;
                block_jacobian[cell * independent_components_per_cell * independent_components_per_cell + row * independent_components_per_cell + column] = if (probe_step[cell] == 0) @as(f64, if (row == column) 1 else 0) else (probe_residual[residual_index] - residual[residual_index]) / probe_step[cell];
            };
        }
        if (block_jacobian_valid) {
            @memset(block_delta, 0);
            for (0..cells) |cell| {
                const matrix = block_jacobian[cell * independent_components_per_cell * independent_components_per_cell ..][0 .. independent_components_per_cell * independent_components_per_cell];
                const rhs = block_right_hand_side[cell * independent_components_per_cell ..][0..independent_components_per_cell];
                for (independent_component, 0..) |component, row| rhs[row] = if (block_coordinate_fixed[cell * independent_components_per_cell + row]) 0 else -residual[component * cells + cell];
                for (0..independent_components_per_cell) |diagonal| {
                    const index = diagonal * independent_components_per_cell + diagonal;
                    matrix[index] += std.math.sqrt(std.math.floatEps(f64)) * @max(1.0, @abs(matrix[index]));
                }
                if (!numerics.solveDenseLinearSystem(matrix, rhs, independent_components_per_cell)) {
                    block_jacobian_valid = false;
                    break;
                }
                for (independent_component, 0..) |component, row| block_delta[component * cells + cell] = rhs[row];
                const ice_water_equivalent_factor = iceWaterEquivalentFactor(properties);
                block_delta[cell] = -(rhs[0] + ice_water_equivalent_factor * rhs[1] + rhs[2] + ice_water_equivalent_factor * rhs[3]);
            }
        }
        if (!block_jacobian_valid) reduced_block_invalid_iterations += 1;
        if (block_jacobian_valid) {
            var line_fraction: f64 = 1;
            var line: u8 = 0;
            while (line < 12) : (line += 1) {
                addCellFeasibleBlockDirection(
                    grid,
                    current,
                    block_delta,
                    line_fraction,
                    candidate,
                ) catch {
                    line_fraction *= 0.5;
                    continue;
                };
                if (residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        reduced_block_newton_steps += 1;
                        accepted_newton = true;
                        break;
                    }
                } else |_| {}
                line_fraction *= 0.5;
            }
            if (!accepted_newton) reduced_block_rejected_iterations += 1;
        }
        if (accepted_newton) continue;
        // Active-set switches can make the full block Jacobian singular even
        // though the limiting scalar equation has a well-defined one-sided
        // derivative. Try that semismooth diagonal Newton correction before
        // reverting to relaxed Picard iteration.
        var limiting_index: usize = 0;
        var limiting_scaled: f64 = -1;
        for (residual, current, 0..) |difference, value, index| {
            const scale = options.absolute_tolerance_m3 + options.relative_tolerance * @max(1.0, @abs(value));
            const scaled = @abs(difference) / scale;
            if (scaled > limiting_scaled) {
                limiting_scaled = scaled;
                limiting_index = index;
            }
        }
        if (limiting_index / cells == 5) {
            @memcpy(probe, current);
            const diagonal_probe_step = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0e-6, @abs(current[limiting_index]));
            probe[limiting_index] += diagonal_probe_step;
            if (residualAt(grid, properties, base, probe, target, probe_residual, scratch, trial_heat, trial_exchange)) |_| {
                var diagonal_derivative = (probe_residual[limiting_index] - residual[limiting_index]) / diagonal_probe_step;
                @memcpy(candidate, current);
                candidate[limiting_index] -= diagonal_probe_step;
                if (candidate[limiting_index] > 0) {
                    if (residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange)) |_| {
                        diagonal_derivative = (probe_residual[limiting_index] - candidate_residual[limiting_index]) / (2 * diagonal_probe_step);
                    } else |_| {}
                }
                if (std.math.isFinite(diagonal_derivative) and @abs(diagonal_derivative) > std.math.floatEps(f64)) {
                    const diagonal_delta = -residual[limiting_index] / diagonal_derivative;
                    var line_fraction: f64 = 1;
                    var line: u8 = 0;
                    while (line < 16) : (line += 1) {
                        @memcpy(candidate, current);
                        candidate[limiting_index] += line_fraction * diagonal_delta;
                        if (!std.math.isFinite(candidate[limiting_index]) or candidate[limiting_index] < 0) {
                            line_fraction *= 0.5;
                            continue;
                        }
                        if (residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange)) |_| {
                            if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                                @memcpy(current, candidate);
                                newton_steps += 1;
                                diagonal_newton_steps += 1;
                                accepted_newton = true;
                                break;
                            }
                        } else |_| {}
                        line_fraction *= 0.5;
                    }
                }
            } else |_| {}
        }
        if (accepted_newton) continue;
        // A phase coordinate can be locally active even when the simultaneous
        // five-coordinate block crosses another cell's active-set boundary.
        // Probe only the limiting cell while eliminating matrix liquid through
        // exact water-equivalent conservation. This is a local semismooth
        // Newton direction, not another model timestep.
        const limiting_component = limiting_index / cells;
        const limiting_phase_cell = limiting_index % cells;
        if (limiting_component >= 1 and limiting_component <= 4) {
            const conservation_coefficient = switch (limiting_component) {
                1, 3 => 1.0,
                2, 4 => iceWaterEquivalentFactor(properties),
                else => unreachable,
            };
            var local_probe_step = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0e-6, @abs(current[limiting_index]));
            if (current[limiting_phase_cell] < conservation_coefficient * local_probe_step) local_probe_step = -@min(local_probe_step, 0.5 * current[limiting_index]);
            if (@abs(local_probe_step) > 1.0e-20) {
                @memcpy(probe, current);
                probe[limiting_index] += local_probe_step;
                probe[limiting_phase_cell] -= conservation_coefficient * local_probe_step;
                if (probe[limiting_index] >= 0 and probe[limiting_phase_cell] >= 0) {
                    if (residualAt(grid, properties, base, probe, target, probe_residual, scratch, trial_heat, trial_exchange)) |_| {
                        const derivative = (probe_residual[limiting_index] - residual[limiting_index]) / local_probe_step;
                        if (std.math.isFinite(derivative) and @abs(derivative) > std.math.floatEps(f64)) {
                            const local_delta = -residual[limiting_index] / derivative;
                            var line_fraction: f64 = 1;
                            var line: u8 = 0;
                            while (line < 16) : (line += 1) {
                                @memcpy(candidate, current);
                                candidate[limiting_index] += line_fraction * local_delta;
                                candidate[limiting_phase_cell] -= conservation_coefficient * line_fraction * local_delta;
                                if (candidate[limiting_index] < 0 or candidate[limiting_phase_cell] < 0 or !std.math.isFinite(candidate[limiting_index]) or !std.math.isFinite(candidate[limiting_phase_cell])) {
                                    line_fraction *= 0.5;
                                    continue;
                                }
                                if (residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange)) |_| {
                                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                                        @memcpy(current, candidate);
                                        newton_steps += 1;
                                        directional_newton_steps += 1;
                                        accepted_newton = true;
                                        break;
                                    }
                                } else |_| {}
                                line_fraction *= 0.5;
                            }
                        }
                    } else |_| {}
                }
            }
        }
        if (accepted_newton) continue;
        // Try the actual fixed point first. Backtracking supplies damping when
        // an active-set transition makes it unstable; always starting at 0.5
        // needlessly halves a smooth freeze/thaw residual until NPH expires.
        var picard_fraction: f64 = 1;
        var accepted_picard = false;
        while (picard_fraction >= 1.0e-6) : (picard_fraction = if (picard_fraction == 1) @min(options.picard_relaxation, 0.5) else picard_fraction * 0.5) {
            addDirection(current, residual, picard_fraction, candidate) catch continue;
            residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange) catch continue;
            if (try scaledNorm(candidate, candidate_residual, options) >= norm) continue;
            @memcpy(current, candidate);
            accepted_picard = true;
            break;
        }
        // FINHL is clipped independently by the liquid and air stores in each
        // cell. At a clip transition, a relaxation factor that improves one
        // block can worsen another. Advance only the worst 5-equation block
        // before declaring stagnation; this is block Picard, not another
        // sub-hourly model cycle.
        if (!accepted_picard) {
            var limiting_cell: usize = 0;
            var limiting_cell_norm: f64 = -1;
            for (0..cells) |cell| {
                var cell_norm: f64 = 0;
                for (0..components_per_cell) |component| {
                    const index = component * cells + cell;
                    const scale = options.absolute_tolerance_m3 + options.relative_tolerance * @max(@abs(current[index]), @abs(current[index] + residual[index]));
                    cell_norm = @max(cell_norm, @abs(residual[index]) / scale);
                }
                if (cell_norm > limiting_cell_norm) {
                    limiting_cell_norm = cell_norm;
                    limiting_cell = cell;
                }
            }
            var block_picard_fraction: f64 = 1;
            while (block_picard_fraction >= 1.0e-6) : (block_picard_fraction = if (block_picard_fraction == 1) @min(options.picard_relaxation, 0.5) else block_picard_fraction * 0.5) {
                @memcpy(candidate, current);
                var valid = true;
                for (0..components_per_cell) |component| {
                    const index = component * cells + limiting_cell;
                    candidate[index] += block_picard_fraction * residual[index];
                    if (!std.math.isFinite(candidate[index]) or candidate[index] < 0) valid = false;
                }
                if (!valid) continue;
                residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange) catch continue;
                if (try scaledNorm(candidate, candidate_residual, options) >= norm) continue;
                @memcpy(current, candidate);
                accepted_picard = true;
                break;
            }
        }
        // At a semismooth phase or pore-capacity switch, even the six-value
        // cell block can straddle two active sets. Advance only its largest
        // currently active equation before reporting stagnation.
        if (!accepted_picard) {
            var coordinate_index: usize = 0;
            var coordinate_scaled: f64 = -1;
            for (residual, current, 0..) |difference, value, index| {
                const scale = options.absolute_tolerance_m3 + options.relative_tolerance * @max(1.0, @abs(value));
                const scaled = @abs(difference) / scale;
                if (scaled > coordinate_scaled) {
                    coordinate_scaled = scaled;
                    coordinate_index = index;
                }
            }
            var coordinate_fraction: f64 = 1;
            while (coordinate_fraction >= 1.0e-8) : (coordinate_fraction = if (coordinate_fraction == 1) @min(options.picard_relaxation, 0.5) else coordinate_fraction * 0.5) {
                @memcpy(candidate, current);
                candidate[coordinate_index] += coordinate_fraction * residual[coordinate_index];
                if (!std.math.isFinite(candidate[coordinate_index]) or candidate[coordinate_index] < 0) continue;
                residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange) catch continue;
                if (try scaledNorm(candidate, candidate_residual, options) >= norm) continue;
                @memcpy(current, candidate);
                accepted_picard = true;
                break;
            }
        }
        // Last local fallback: jointly correct the limiting phase coordinate
        // and endpoint temperature in one cell. Eliminating matrix liquid
        // preserves water exactly, while the 2x2 block captures latent-heat
        // feedback that a scalar coordinate step cannot.
        if (!accepted_picard) local_phase_temperature: {
            var phase_index: usize = 0;
            var phase_scaled: f64 = -1;
            for (1..5) |component| for (0..cells) |cell| {
                const index = component * cells + cell;
                const scale = options.absolute_tolerance_m3 + options.relative_tolerance * @max(1.0, @abs(current[index]));
                const scaled = @abs(residual[index]) / scale;
                if (scaled > phase_scaled) {
                    phase_scaled = scaled;
                    phase_index = index;
                }
            };
            const phase_component = phase_index / cells;
            const cell = phase_index % cells;
            const temperature_index = 5 * cells + cell;
            const conservation_coefficient = switch (phase_component) {
                1, 3 => 1.0,
                2, 4 => iceWaterEquivalentFactor(properties),
                else => break :local_phase_temperature,
            };
            var phase_step = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0e-6, @abs(current[phase_index]));
            if (current[cell] < conservation_coefficient * phase_step) phase_step = -@min(phase_step, 0.5 * current[phase_index]);
            if (@abs(phase_step) <= 1.0e-20) break :local_phase_temperature;
            @memcpy(probe, current);
            probe[phase_index] += phase_step;
            probe[cell] -= conservation_coefficient * phase_step;
            residualAt(grid, properties, base, probe, target, probe_residual, scratch, trial_heat, trial_exchange) catch break :local_phase_temperature;
            const local_matrix = block_jacobian[0..4];
            const local_rhs = block_right_hand_side[0..2];
            local_matrix[0] = (probe_residual[phase_index] - residual[phase_index]) / phase_step;
            local_matrix[2] = (probe_residual[temperature_index] - residual[temperature_index]) / phase_step;
            const temperature_step = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0, @abs(current[temperature_index]));
            @memcpy(probe, current);
            probe[temperature_index] += temperature_step;
            residualAt(grid, properties, base, probe, target, probe_residual, scratch, trial_heat, trial_exchange) catch break :local_phase_temperature;
            local_matrix[1] = (probe_residual[phase_index] - residual[phase_index]) / temperature_step;
            local_matrix[3] = (probe_residual[temperature_index] - residual[temperature_index]) / temperature_step;
            local_rhs[0] = -residual[phase_index];
            local_rhs[1] = -residual[temperature_index];
            if (!numerics.solveDenseLinearSystem(local_matrix, local_rhs, 2)) break :local_phase_temperature;
            var fraction: f64 = 1;
            while (fraction >= 1.0e-8) : (fraction *= 0.5) {
                @memcpy(candidate, current);
                candidate[phase_index] += fraction * local_rhs[0];
                candidate[cell] -= conservation_coefficient * fraction * local_rhs[0];
                candidate[temperature_index] += fraction * local_rhs[1];
                if (candidate[phase_index] < 0 or candidate[cell] < 0 or candidate[temperature_index] <= 0) continue;
                residualAt(grid, properties, base, candidate, target, candidate_residual, scratch, trial_heat, trial_exchange) catch continue;
                if (try scaledNorm(candidate, candidate_residual, options) >= norm) continue;
                @memcpy(current, candidate);
                accepted_picard = true;
                newton_steps += 1;
                directional_newton_steps += 1;
                break;
            }
        }
        if (!accepted_picard) {
            if (!builtin.is_test) {
                var largest_index: usize = 0;
                var largest_scaled: f64 = -1;
                for (residual, current, 0..) |difference, value, index| {
                    const scale = options.absolute_tolerance_m3 + options.relative_tolerance * @max(@abs(value), @abs(value + difference));
                    const scaled = @abs(difference) / scale;
                    if (scaled > largest_scaled) {
                        largest_scaled = scaled;
                        largest_index = index;
                    }
                }
                const component_names = [_][]const u8{ "matrix_liquid_water_m3", "water_vapor_volume_m3", "matrix_ice_water_m3", "macropore_liquid_water_m3", "macropore_ice_volume_m3", "endpoint_temperature_k" };
                std.log.err("soil phase-enthalpy Newton-Picard stagnated: iteration={d} scaled_residual={e} limiting_component={s} layer_cell={d} state={e} residual={e} target={e}", .{ iteration + 1, norm, component_names[largest_index / cells], largest_index % cells, current[largest_index], residual[largest_index], target[largest_index] });
                const limiting_cell = largest_index % cells;
                std.log.err("soil phase limiting block: temperature_k={e}->{e} matrix_water={e}->{e} vapor={e}->{e} matrix_ice={e}->{e} macropore_water={e}->{e} macropore_ice={e}->{e}", .{ current[5 * cells + limiting_cell], target[5 * cells + limiting_cell], current[limiting_cell], target[limiting_cell], current[cells + limiting_cell], target[cells + limiting_cell], current[2 * cells + limiting_cell], target[2 * cells + limiting_cell], current[3 * cells + limiting_cell], target[3 * cells + limiting_cell], current[4 * cells + limiting_cell], target[4 * cells + limiting_cell] });
            }
            return error.SoilPhaseSolverStagnated;
        }
        picard_steps += 1;
    }
    try residualAt(grid, properties, base, current, target, residual, scratch, trial_heat, trial_exchange);
    const final_norm = try scaledNorm(current, residual, options);
    if (final_norm <= 1) {
        try residualAt(grid, properties, base, current, target, residual, scratch, outputs.latent_heat_mj, outputs.macropore_to_matrix_water_m3);
        try commit(grid, target);
        return .{ .iterations = options.max_iterations, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = final_norm };
    }
    if (!builtin.is_test) {
        var limiting_index: usize = 0;
        var limiting_scaled: f64 = -1;
        for (residual, current, 0..) |difference, value, index| {
            const scale = options.absolute_tolerance_m3 + options.relative_tolerance * @max(1.0, @abs(value));
            const scaled = @abs(difference) / scale;
            if (scaled > limiting_scaled) {
                limiting_scaled = scaled;
                limiting_index = index;
            }
        }
        const component_names = [_][]const u8{ "matrix_liquid_water_m3", "water_vapor_volume_m3", "matrix_ice_water_m3", "macropore_liquid_water_m3", "macropore_ice_volume_m3", "endpoint_temperature_k" };
        std.log.err("soil phase-enthalpy iteration ceiling reached: iterations={d} directional_newton={d} reduced_block_newton={d} reduced_block_invalid={d} reduced_block_rejected={d} diagonal_newton={d} picard={d} scaled_residual={e} limiting_component={s} layer_cell={d} state={e} residual={e} target={e}", .{ options.max_iterations, directional_newton_steps, reduced_block_newton_steps, reduced_block_invalid_iterations, reduced_block_rejected_iterations, diagonal_newton_steps, picard_steps, final_norm, component_names[limiting_index / cells], limiting_index % cells, current[limiting_index], residual[limiting_index], target[limiting_index] });
        if (last_reduced_block_probe_error) |err| std.log.err("last reduced phase block probe error: {s}", .{@errorName(err)});
    }
    return error.SoilPhaseSolverDidNotConverge;
}

fn residualAt(grid: *const grid_module.GridState, properties: Properties, base: []const f64, trial: []const f64, target: []f64, residual: []f64, scratch: []f64, latent_heat: []f64, exchange_flux: []f64) !void {
    const cells = grid.layer_count;
    @memcpy(scratch, trial);
    @memcpy(target, base);
    @memset(latent_heat, 0);
    @memset(exchange_flux, 0);
    for (0..cells) |cell| {
        const matrix_water = cell;
        const vapor = cells + cell;
        const matrix_ice = 2 * cells + cell;
        const macro_water = 3 * cells + cell;
        const macro_ice = 4 * cells + cell;
        const temperature = 5 * cells + cell;
        const matrix_water_fraction = scratch[matrix_water] / properties.matrix_bulk_volume_m3[cell];
        const matrix_parameters = properties.mualem_van_genuchten_parameters[cell];
        const matric_plus_osmotic_potential_mpa =
            try matrix_parameters.pressureHeadAtWaterContent(std.math.clamp(
                matrix_water_fraction,
                matrix_parameters.residual_water_content_m3_per_m3,
                matrix_parameters.saturated_water_content_m3_per_m3,
            )) * properties.gravitational_water_potential_mpa_per_m +
            properties.osmotic_potential_mpa[cell];
        const matrix_air_before_phase_m3 = @max(0.0, grid.matrix_pore_capacity_m3[cell] - scratch[matrix_water] - scratch[matrix_ice]);
        const total_air_before_phase_m3 = matrix_air_before_phase_m3 + @max(0.0, grid.macropore_pore_capacity_m3[cell] - scratch[macro_water] - scratch[macro_ice]);
        var vapor_change = phase.vaporLiquidEquilibrium(scratch[temperature], matric_plus_osmotic_potential_mpa, scratch[vapor], total_air_before_phase_m3, scratch[matrix_water], 1, properties.vapor) catch |err| {
            if (!builtin.is_test) std.log.err("soil vapor equilibrium failed: layer_cell={d} temperature_k={e} potential_mpa={e} vapor_m3={e} air_m3={e} matrix_water_m3={e} error={s}", .{ cell, scratch[temperature], matric_plus_osmotic_potential_mpa, scratch[vapor], total_air_before_phase_m3, scratch[matrix_water], @errorName(err) });
            return err;
        };
        if (vapor_change.water_condensation_m3 > matrix_air_before_phase_m3) {
            const scale = matrix_air_before_phase_m3 / vapor_change.water_condensation_m3;
            vapor_change.water_condensation_m3 *= scale;
            vapor_change.vapor_change_m3 *= scale;
            vapor_change.latent_heat_mj *= scale;
        }
        applyChange(scratch, target, matrix_water, vapor_change.water_condensation_m3);
        applyChange(scratch, target, vapor, vapor_change.vapor_change_m3);
        const matrix_freeze = phase.FreezeThaw{
            .freezing_temperature_k = properties.freeze_thaw.pure_water_freezing_temperature_k,
            .liquid_water_change_m3 = 0,
            .ice_volume_change_m3 = 0,
            .latent_heat_mj = 0,
        };
        applyChange(scratch, target, matrix_water, matrix_freeze.liquid_water_change_m3);
        applyChange(scratch, target, matrix_ice, matrix_freeze.ice_volume_change_m3);
        const macro_freeze = phase.FreezeThaw{
            .freezing_temperature_k = properties.freeze_thaw.pure_water_freezing_temperature_k,
            .liquid_water_change_m3 = 0,
            .ice_volume_change_m3 = 0,
            .latent_heat_mj = 0,
        };
        applyChange(scratch, target, macro_water, macro_freeze.liquid_water_change_m3);
        applyChange(scratch, target, macro_ice, macro_freeze.ice_volume_change_m3);
        if (properties.pore_exchange_enabled.len != 0 and
            properties.pore_exchange_enabled[cell] and
            grid.macropore_pore_capacity_m3[cell] > 0)
        {
            const matrix_water_for_exchange = if (scratch[matrix_water] >= -1e-12) @max(0.0, scratch[matrix_water]) else return error.InvalidSoilPhaseCandidate;
            const macropore_water_for_exchange = if (scratch[macro_water] >= -1e-12) @max(0.0, scratch[macro_water]) else return error.InvalidSoilPhaseCandidate;
            const matrix_air = @max(0.0, grid.matrix_pore_capacity_m3[cell] - matrix_water_for_exchange - scratch[matrix_ice]);
            const macro_air = @max(0.0, grid.macropore_pore_capacity_m3[cell] - macropore_water_for_exchange - scratch[macro_ice]);
            const exchange_water_fraction = matrix_water_for_exchange / properties.matrix_bulk_volume_m3[cell];
            const exchange_matric_potential_mpa = try properties.retention_curve[cell].waterPotentialMpa(@max(exchange_water_fraction, std.math.floatMin(f64)));
            const exchange = phase.macroporeMatrixExchange(.{ .saturated_lateral_matrix_conductivity_m2_per_h_mpa = properties.saturated_lateral_matrix_conductivity_m2_per_h_mpa[cell], .face_area_m2 = properties.face_area_m2[cell], .saturation_water_potential_mpa = properties.saturation_water_potential_mpa[cell], .current_matric_potential_mpa = exchange_matric_potential_mpa, .macropore_spacing_m = properties.macropore_spacing_m[cell], .macropore_radius_m = properties.macropore_radius_m[cell], .time_fraction = 1, .matrix_water_m3 = matrix_water_for_exchange, .matrix_air_m3 = matrix_air, .macropore_water_m3 = macropore_water_for_exchange, .macropore_air_m3 = macro_air }) catch |err| {
                if (!builtin.is_test) std.log.err("macropore-matrix exchange failed: layer_cell={d} conductivity={e} face_area_m2={e} spacing_m={e} radius_m={e} matrix_water_m3={e} matrix_air_m3={e} macropore_water_m3={e} macropore_air_m3={e} error={s}", .{ cell, properties.saturated_lateral_matrix_conductivity_m2_per_h_mpa[cell], properties.face_area_m2[cell], properties.macropore_spacing_m[cell], properties.macropore_radius_m[cell], matrix_water_for_exchange, matrix_air, macropore_water_for_exchange, macro_air, @errorName(err) });
                return err;
            };
            applyChange(scratch, target, matrix_water, exchange);
            applyChange(scratch, target, macro_water, -exchange);
            exchange_flux[cell] = exchange;
        }
        const matrix_occupancy_m3 =
            scratch[matrix_water] + scratch[matrix_ice];
        const macropore_occupancy_m3 =
            scratch[macro_water] + scratch[macro_ice];
        if (matrix_occupancy_m3 >
            grid.matrix_pore_capacity_m3[cell] +
                poreCapacityRoundoffToleranceM3(
                    grid.matrix_pore_capacity_m3[cell],
                ) or
            macropore_occupancy_m3 >
                grid.macropore_pore_capacity_m3[cell] +
                    poreCapacityRoundoffToleranceM3(
                        grid.macropore_pore_capacity_m3[cell],
                    ))
        {
            if (!builtin.is_test) std.log.err(
                "soil phase candidate exceeds pore capacity: layer_cell={d} matrix_liquid_m3={e} matrix_ice_m3={e} matrix_capacity_m3={e} matrix_excess_m3={e} macropore_liquid_m3={e} macropore_ice_m3={e} macropore_capacity_m3={e} macropore_excess_m3={e}",
                .{
                    cell,
                    scratch[matrix_water],
                    scratch[matrix_ice],
                    grid.matrix_pore_capacity_m3[cell],
                    matrix_occupancy_m3 -
                        grid.matrix_pore_capacity_m3[cell],
                    scratch[macro_water],
                    scratch[macro_ice],
                    grid.macropore_pore_capacity_m3[cell],
                    macropore_occupancy_m3 -
                        grid.macropore_pore_capacity_m3[cell],
                },
            );
            return error.SoilPhaseCandidateExceedsPoreCapacity;
        }
        latent_heat[cell] = vapor_change.latent_heat_mj + matrix_freeze.latent_heat_mj + macro_freeze.latent_heat_mj;
        target[temperature] = try phase.endpointTemperatureFromPhaseEnthalpy(
            grid.soil_temperature_k[cell],
            properties.heat_capacity_mj_per_k[cell],
            .{ .matrix_liquid_water_m3 = base[matrix_water], .water_vapor_volume_m3 = base[vapor], .matrix_ice_volume_m3 = base[matrix_ice], .macropore_liquid_water_m3 = base[macro_water], .macropore_ice_volume_m3 = base[macro_ice] },
            .{ .matrix_liquid_water_m3 = @max(0, target[matrix_water]), .water_vapor_volume_m3 = @max(0, target[vapor]), .matrix_ice_volume_m3 = @max(0, target[matrix_ice]), .macropore_liquid_water_m3 = @max(0, target[macro_water]), .macropore_ice_volume_m3 = @max(0, target[macro_ice]) },
            0,
            .{ .liquid_heat_capacity_mj_per_m3_k = properties.liquid_water_heat_capacity_mj_per_m3_k, .ice_heat_capacity_mj_per_m3_k = properties.ice_heat_capacity_mj_per_m3_k, .vaporization_latent_heat_mj_per_m3 = properties.vapor.latent_heat_of_vaporization_mj_per_m3, .fusion_latent_heat_mj_per_m3 = properties.freeze_thaw.latent_heat_of_fusion_mj_per_m3, .ice_density_megagrams_per_m3 = iceWaterEquivalentFactor(properties) },
        );
    }
    for (target, trial, residual) |*value, trial_value, *difference| {
        if (!std.math.isFinite(value.*) or value.* < -1e-12) return error.InvalidSoilPhaseCandidate;
        // Normalize solver-scale signed zero in the authoritative target, not
        // only in its residual. The converged target is later committed.
        value.* = @max(0.0, value.*);
        difference.* = value.* - trial_value;
    }
}

fn dallAmicoChange(
    temperature_k: f64,
    liquid_water_m3: f64,
    ice_water_equivalent_m3: f64,
    porous_medium_volume_m3: f64,
    parameters: retention.MualemVanGenuchtenParameters,
    freeze_thaw: phase.FreezeThawParameters,
    gravitational_water_potential_mpa_per_m: f64,
) !phase.FreezeThaw {
    const total_water_equivalent_m3 =
        liquid_water_m3 + ice_water_equivalent_m3;
    const total_water_content =
        total_water_equivalent_m3 / porous_medium_volume_m3;
    if (total_water_content <= parameters.residual_water_content_m3_per_m3) {
        return .{
            .freezing_temperature_k = freeze_thaw.pure_water_freezing_temperature_k,
            .liquid_water_change_m3 = ice_water_equivalent_m3,
            .ice_volume_change_m3 = -ice_water_equivalent_m3,
            .latent_heat_mj = -freeze_thaw.latent_heat_of_fusion_mj_per_m3 *
                ice_water_equivalent_m3,
        };
    }
    const bounded_water_content = @min(
        total_water_content,
        parameters.saturated_water_content_m3_per_m3,
    );
    const unfrozen_pressure_head_m =
        try parameters.pressureHeadAtWaterContent(bounded_water_content);
    const equilibrium = try phase.dallAmicoEquilibrium(.{
        .temperature_k = temperature_k,
        .total_water_equivalent_m3 = total_water_equivalent_m3,
        .porous_medium_volume_m3 = porous_medium_volume_m3,
        .unfrozen_pressure_head_m = unfrozen_pressure_head_m,
        .gravitational_water_potential_mpa_per_m = gravitational_water_potential_mpa_per_m,
        .latent_heat_of_fusion_mj_per_m3 = freeze_thaw.latent_heat_of_fusion_mj_per_m3,
        .pure_water_melting_temperature_k = freeze_thaw.pure_water_freezing_temperature_k,
        .mualem_van_genuchten = parameters,
    });
    const liquid_change_m3 =
        equilibrium.liquid_water_m3 - liquid_water_m3;
    const ice_change_m3 =
        equilibrium.ice_water_equivalent_m3 - ice_water_equivalent_m3;
    return .{
        .freezing_temperature_k = equilibrium.depressed_melting_temperature_k,
        .liquid_water_change_m3 = liquid_change_m3,
        .ice_volume_change_m3 = ice_change_m3,
        .latent_heat_mj = freeze_thaw.latent_heat_of_fusion_mj_per_m3 * ice_change_m3,
    };
}

/// Ice occupies more pore volume than the liquid water that formed it. At
/// saturation WATSUB cannot create pore volume, so freezing and latent heat
/// are reduced together to the available air volume.
fn limitFreezingExpansion(change: anytype, liquid_water_m3: f64, ice_volume_m3: f64, pore_capacity_m3: f64) void {
    if (change.liquid_water_change_m3 >= 0) return;
    const volume_expansion_m3 = change.liquid_water_change_m3 + change.ice_volume_change_m3;
    if (volume_expansion_m3 <= 0) return;
    const available_air_m3 = @max(0.0, pore_capacity_m3 - liquid_water_m3 - ice_volume_m3);
    if (volume_expansion_m3 <= available_air_m3) return;
    const scale = available_air_m3 / volume_expansion_m3;
    change.liquid_water_change_m3 *= scale;
    change.ice_volume_change_m3 *= scale;
    change.latent_heat_mj *= scale;
}

fn poreCapacityRoundoffToleranceM3(capacity_m3: f64) f64 {
    return @max(
        1.0e-12,
        64.0 * std.math.floatEps(f64) *
            @max(1.0, @abs(capacity_m3)),
    );
}

fn applyChange(scratch: []f64, target: []f64, index: usize, change: f64) void {
    scratch[index] += change;
    target[index] += change;
}

fn commit(grid: *grid_module.GridState, state: []const f64) !void {
    const cells = grid.layer_count;
    @memcpy(grid.matrix_liquid_water_m3, state[0 * cells .. 1 * cells]);
    @memcpy(grid.water_vapor_volume_m3, state[1 * cells .. 2 * cells]);
    @memcpy(grid.matrix_ice_water_m3, state[2 * cells .. 3 * cells]);
    @memcpy(grid.macropore_liquid_water_m3, state[3 * cells .. 4 * cells]);
    @memcpy(grid.macropore_ice_water_m3, state[4 * cells .. 5 * cells]);
    for (0..cells) |cell| {
        grid.liquid_water_m3[cell] = grid.matrix_liquid_water_m3[cell] + grid.macropore_liquid_water_m3[cell];
        grid.ice_water_m3[cell] = grid.matrix_ice_water_m3[cell] + grid.macropore_ice_water_m3[cell];
        grid.matrix_air_volume_m3[cell] = @max(0.0, grid.matrix_pore_capacity_m3[cell] - grid.matrix_liquid_water_m3[cell] - grid.matrix_ice_water_m3[cell]);
        grid.macropore_air_volume_m3[cell] = @max(0.0, grid.macropore_pore_capacity_m3[cell] - grid.macropore_liquid_water_m3[cell] - grid.macropore_ice_water_m3[cell]);
        grid.air_volume_m3[cell] = grid.matrix_air_volume_m3[cell] + grid.macropore_air_volume_m3[cell];
    }
    try grid.validateFinite();
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidSoilPhaseCandidate;
        candidate.* = @max(0.0, candidate.*);
    }
}

/// Apply the independent cell-local Newton blocks with a separate feasible
/// fraction for each cell. Soil cells are uncoupled in this phase residual;
/// forcing every cell to use the smallest globally feasible line fraction
/// lets one saturated active set stall all other cells.
fn addCellFeasibleBlockDirection(
    grid: *const grid_module.GridState,
    current: []const f64,
    direction: []const f64,
    requested_fraction: f64,
    output: []f64,
) !void {
    const cells = grid.layer_count;
    if (current.len != 6 * cells or
        direction.len != current.len or
        output.len != current.len or
        !std.math.isFinite(requested_fraction) or
        requested_fraction <= 0)
        return error.InvalidSoilPhaseCandidate;
    @memcpy(output, current);
    for (0..cells) |cell| {
        var fraction = requested_fraction;
        for (0..6) |component| {
            const index = component * cells + cell;
            const delta = direction[index];
            if (!std.math.isFinite(delta)) return error.InvalidSoilPhaseCandidate;
            if (delta < 0)
                fraction = @min(fraction, current[index] / -delta);
        }
        const matrix_occupancy_m3 =
            current[cell] + current[2 * cells + cell];
        const matrix_occupancy_change_m3 =
            direction[cell] + direction[2 * cells + cell];
        if (matrix_occupancy_change_m3 > 0)
            fraction = @min(
                fraction,
                @max(0, grid.matrix_pore_capacity_m3[cell] -
                    matrix_occupancy_m3) /
                    matrix_occupancy_change_m3,
            );
        const macropore_occupancy_m3 =
            current[3 * cells + cell] + current[4 * cells + cell];
        const macropore_occupancy_change_m3 =
            direction[3 * cells + cell] + direction[4 * cells + cell];
        if (macropore_occupancy_change_m3 > 0)
            fraction = @min(
                fraction,
                @max(0, grid.macropore_pore_capacity_m3[cell] -
                    macropore_occupancy_m3) /
                    macropore_occupancy_change_m3,
            );
        if (!std.math.isFinite(fraction) or fraction < 0)
            return error.InvalidSoilPhaseCandidate;
        // Stay a few ulps inside active constraints so residual evaluation
        // cannot reject an otherwise feasible Newton block after rounding.
        if (fraction < requested_fraction)
            fraction *= 1.0 - 16.0 * std.math.floatEps(f64);
        for (0..6) |component| {
            const index = component * cells + cell;
            output[index] =
                @max(0, current[index] + fraction * direction[index]);
            if (!std.math.isFinite(output[index]) or
                (component == 5 and output[index] <= 0))
                return error.InvalidSoilPhaseCandidate;
        }
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteSoilPhaseSolverState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_m3 + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn validateInputs(grid: *const grid_module.GridState, properties: Properties, outputs: Outputs, options: Options) !void {
    const cells = grid.layer_count;
    inline for (@typeInfo(Properties).@"struct".fields) |field| {
        if (field.type == []const f64 and @field(properties, field.name).len != cells) return error.SoilPhaseSolverDimensionMismatch;
        if (field.type == []const bool and @field(properties, field.name).len != 0 and @field(properties, field.name).len != cells) return error.SoilPhaseSolverDimensionMismatch;
    }
    if (properties.retention_curve.len != cells) return error.SoilPhaseSolverDimensionMismatch;
    if (properties.mualem_van_genuchten_parameters.len != cells or
        properties.macropore_mualem_van_genuchten_parameters.len != cells)
    {
        return error.SoilPhaseSolverDimensionMismatch;
    }
    for (properties.mualem_van_genuchten_parameters) |parameters|
        try parameters.validate();
    for (properties.macropore_mualem_van_genuchten_parameters) |parameters|
        try parameters.validate();
    if (properties.freeze_thaw.ice_density_megagrams_per_m3 <= 0 or properties.freeze_thaw.ice_density_megagrams_per_m3 >= 1) return error.InvalidSoilPhaseInput;
    if (!std.math.isFinite(properties.gravitational_water_potential_mpa_per_m) or
        properties.gravitational_water_potential_mpa_per_m <= 0)
        return error.InvalidSoilPhaseInput;
    if (!std.math.isFinite(properties.liquid_water_heat_capacity_mj_per_m3_k) or properties.liquid_water_heat_capacity_mj_per_m3_k <= 0 or !std.math.isFinite(properties.ice_heat_capacity_mj_per_m3_k) or properties.ice_heat_capacity_mj_per_m3_k <= 0) return error.InvalidSoilPhaseInput;
    for (properties.matrix_bulk_volume_m3, properties.retention_curve) |bulk_volume_m3, curve| {
        if (!std.math.isFinite(bulk_volume_m3) or bulk_volume_m3 <= 0) return error.InvalidSoilPhaseInput;
        if (!std.math.isFinite(curve.porosity_fraction) or curve.porosity_fraction <= 0) return error.InvalidSoilPhaseInput;
    }
    if (outputs.latent_heat_mj.len != cells or outputs.macropore_to_matrix_water_m3.len != cells) return error.SoilPhaseSolverDimensionMismatch;
    if (options.max_iterations == 0 or !std.math.isFinite(options.absolute_tolerance_m3) or options.absolute_tolerance_m3 <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction) return error.InvalidSoilPhaseSolverOptions;
}

fn testProperties() Properties {
    const values = struct {
        const bulk = [_]f64{2};
        const osmotic = [_]f64{0};
        const curves = [_]retention.ResolvedCurve{.{ .porosity_fraction = 1, .curve = .{ .field_capacity_fraction = 0.6, .wilting_point_fraction = 0.2, .saturation_water_potential_mpa = -0.0005, .field_capacity_water_potential_mpa = -0.01, .wilting_point_water_potential_mpa = -1.5, .minimum_water_potential_mpa = -1.5e12, .saturation_to_field_shape = 0.5, .below_wilting_shape = 0.5 } }};
        const saturation_potential = [_]f64{-0.0005};
        // Includes a positive dry-solid contribution in addition to water.
        const capacity = [_]f64{5.0};
        const one = [_]f64{1};
        const spacing = [_]f64{1};
        const radius = [_]f64{0.01};
        const disabled = [_]bool{false};
        const matrix = [_]retention.MualemVanGenuchtenParameters{.{
            .residual_water_content_m3_per_m3 = 0.05,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 1.6,
            .n = 1.6,
            .saturated_hydraulic_conductivity_m_per_h = 0.01,
        }};
        const macropore = [_]retention.MualemVanGenuchtenParameters{.{
            .residual_water_content_m3_per_m3 = 0,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 15,
            .n = 2.68,
            .saturated_hydraulic_conductivity_m_per_h = 0.1,
        }};
    };
    return .{ .matrix_bulk_volume_m3 = &values.bulk, .retention_curve = &values.curves, .mualem_van_genuchten_parameters = &values.matrix, .macropore_mualem_van_genuchten_parameters = &values.macropore, .osmotic_potential_mpa = &values.osmotic, .saturation_water_potential_mpa = &values.saturation_potential, .heat_capacity_mj_per_k = &values.capacity, .saturated_lateral_matrix_conductivity_m2_per_h_mpa = &values.one, .face_area_m2 = &values.one, .macropore_spacing_m = &values.spacing, .macropore_radius_m = &values.radius, .pore_exchange_enabled = &values.disabled, .vapor = .{ .vapor_density_temperature_coefficient = 2.173e-3, .molecular_weight_ratio = 0.61, .clausius_clapeyron_coefficient_k = 5360, .reference_inverse_temperature_per_k = 3.661e-3, .water_molar_mass_g_per_mol = 18, .gas_constant_j_per_mol_k = 8.3143, .latent_heat_of_vaporization_mj_per_m3 = 2450 }, .freeze_thaw = .{ .freezing_potential_numerator_k_mpa = 9.0959e4, .latent_heat_of_fusion_mj_per_m3 = 333, .ice_density_megagrams_per_m3 = 0.917, .heat_capacity_temperature_feedback_per_k = 6.2913e-3, .pure_water_freezing_temperature_k = 273.15 }, .liquid_water_heat_capacity_mj_per_m3_k = 4.19, .ice_heat_capacity_mj_per_m3_k = 1.9274 };
}

fn dallAmicoTestProperties() Properties {
    const values = struct {
        const matrix = [_]retention.MualemVanGenuchtenParameters{.{
            .residual_water_content_m3_per_m3 = 0.05,
            .saturated_water_content_m3_per_m3 = 0.45,
            .alpha_per_m = 1.6,
            .n = 1.6,
            .saturated_hydraulic_conductivity_m_per_h = 0.01,
        }};
        const macropore = [_]retention.MualemVanGenuchtenParameters{.{
            .residual_water_content_m3_per_m3 = 0,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 15,
            .n = 2.68,
            .saturated_hydraulic_conductivity_m_per_h = 0.1,
        }};
    };
    var properties = testProperties();
    properties.mualem_van_genuchten_parameters = &values.matrix;
    properties.macropore_mualem_van_genuchten_parameters = &values.macropore;
    return properties;
}

test "phase hybrid converges vapor without duplicating heat-enthalpy freeze thaw" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 1;
    grid.water_vapor_volume_m3[0] = 0.01;
    grid.matrix_pore_capacity_m3[0] = 2;
    grid.soil_temperature_k[0] = 260;
    var heat_output = [_]f64{0};
    var exchange = [_]f64{0};
    const before = grid.matrix_liquid_water_m3[0] + grid.matrix_ice_water_m3[0] + grid.water_vapor_volume_m3[0];
    const result = try solve(std.testing.allocator, &grid, testProperties(), .{ .latent_heat_mj = &heat_output, .macropore_to_matrix_water_m3 = &exchange }, .{ .max_iterations = 40 });
    const after = grid.matrix_liquid_water_m3[0] + grid.matrix_ice_water_m3[0] + grid.water_vapor_volume_m3[0];
    try std.testing.expect(result.iterations < 40);
    try std.testing.expectApproxEqAbs(before, after, 1e-10);
    try std.testing.expectEqual(@as(f64, 0), grid.matrix_ice_water_m3[0]);
    try std.testing.expect(heat_output[0] != 0);
}

test "phase hybrid leaves Dall'Amico state for the coupled enthalpy solver" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    const properties = dallAmicoTestProperties();
    grid.matrix_liquid_water_m3[0] = try properties.mualem_van_genuchten_parameters[0].waterContentAtPressureHead(-2) * properties.matrix_bulk_volume_m3[0];
    grid.water_vapor_volume_m3[0] = 1.0e-5;
    grid.matrix_pore_capacity_m3[0] = properties.mualem_van_genuchten_parameters[0].saturated_water_content_m3_per_m3 * properties.matrix_bulk_volume_m3[0];
    grid.soil_temperature_k[0] = 268;
    const before = grid.matrix_liquid_water_m3[0] + grid.matrix_ice_water_m3[0] + grid.water_vapor_volume_m3[0];
    var heat_output = [_]f64{0};
    var exchange = [_]f64{0};
    const result = try solve(std.testing.allocator, &grid, properties, .{ .latent_heat_mj = &heat_output, .macropore_to_matrix_water_m3 = &exchange }, .{ .max_iterations = 80 });
    const after = grid.matrix_liquid_water_m3[0] + grid.matrix_ice_water_m3[0] + grid.water_vapor_volume_m3[0];
    try std.testing.expect(result.iterations < 80);
    try std.testing.expectApproxEqAbs(before, after, 1e-10);
    try std.testing.expectEqual(@as(f64, 0), grid.matrix_ice_water_m3[0]);
    try std.testing.expect(grid.matrix_liquid_water_m3[0] >= properties.mualem_van_genuchten_parameters[0].residual_water_content_m3_per_m3 * properties.matrix_bulk_volume_m3[0]);
}

test "phase solve cannot independently create ice at saturated pore capacity" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 1;
    grid.liquid_water_m3[0] = 1;
    grid.matrix_pore_capacity_m3[0] = 1;
    grid.soil_temperature_k[0] = 260;
    var heat_output = [_]f64{0};
    var exchange = [_]f64{0};
    _ = try solve(std.testing.allocator, &grid, testProperties(), .{ .latent_heat_mj = &heat_output, .macropore_to_matrix_water_m3 = &exchange }, .{ .max_iterations = 40 });
    try std.testing.expectEqual(@as(f64, 0), grid.matrix_ice_water_m3[0]);
    try std.testing.expect(grid.matrix_liquid_water_m3[0] + grid.matrix_ice_water_m3[0] <= grid.matrix_pore_capacity_m3[0] + 1e-12);
}

test "rejected phase solve leaves grid and outputs unchanged" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 1;
    grid.water_vapor_volume_m3[0] = 0.01;
    grid.matrix_pore_capacity_m3[0] = 2;
    grid.soil_temperature_k[0] = 260;
    var heat_output = [_]f64{99};
    var exchange = [_]f64{88};
    try std.testing.expectError(error.InvalidSoilPhaseSolverOptions, solve(std.testing.allocator, &grid, testProperties(), .{ .latent_heat_mj = &heat_output, .macropore_to_matrix_water_m3 = &exchange }, .{ .max_iterations = 0 }));
    try std.testing.expectEqual(@as(f64, 1), grid.matrix_liquid_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 99), heat_output[0]);
    try std.testing.expectEqual(@as(f64, 88), exchange[0]);
}
