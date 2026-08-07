const std = @import("std");
const builtin = @import("builtin");
const numerics = @import("../../core/numerics.zig");
const grid_module = @import("../../state/grid.zig");
const retention = @import("retention.zig");
const water_flux = @import("flux.zig");
const transport_hydrology = @import("../../transport/hydrology.zig");
const water_boundary = @import("boundary.zig");
const boundary_topology = @import("../profile/boundary_topology.zig");

pub const Axis = enum(u2) { x = 0, y = 1, z = 2 };

pub const Face = struct {
    source_cell: usize,
    destination_cell: usize,
    axis: Axis = .x,
    direction: water_flux.FaceDirection,
    source_path_length_m: f64,
    destination_path_length_m: f64,
    face_area_m2: f64,
    macropore_hydraulic_conductance_m_per_h_megapascal: f64,
};

pub const FaceGeometry = struct {
    source_path_length_m: []const f64,
    destination_path_length_m: []const f64,
    face_area_m2: []const f64,
    macropore_hydraulic_conductance_m_per_h_megapascal: []const f64,
};

/// All arrays are runtime-sized by soil cell. Directional conductivity is
/// cell-major x/y/z. Nothing in this WATSUB state has a compile-time grid cap.
pub const Properties = struct {
    matrix_bulk_volume_m3: []const f64,
    retention_curve: []const retention.ResolvedCurve,
    /// Original Mualem–van Genuchten parameters used by every Richards
    /// residual. Both domains are mandatory runtime data.
    mualem_van_genuchten_parameters: []const retention.MualemVanGenuchtenParameters = &.{},
    macropore_mualem_van_genuchten_parameters: []const retention.MualemVanGenuchtenParameters = &.{},
    macropore_spacing_m: []const f64 = &.{},
    macropore_radius_m: []const f64 = &.{},
    dual_domain_exchange_enabled: []const bool = &.{},
    dual_domain_geometry_factor: f64 = 3,
    dual_domain_scaling_coefficient: f64 = 0.4,
    frozen_hydraulic_impedance_exponent: f64 = 0,
    gravitational_water_potential_mpa_per_m: f64 = 0.00980665,
    gravitational_potential_megapascal: []const f64,
    osmotic_potential_megapascal: []const f64,
    matrix_hydraulic_conductivity_m2_per_h_megapascal: []const f64,
    rainfall_conductivity_multiplier: []const f64 = &.{},
    /// Positive extensive liquid-water source entering the matrix domain
    /// during this physical step (for example depth-routed irrigation).
    matrix_external_source_m3_per_step: []const f64 = &.{},
    hydraulic_conductivity_class_count: usize = 1,
    vertical_thickness_m: []const f64,
    osmotic_potential_multiplier: f64,
    nonlinear_time_fraction: f64 = 1,
    boundary_topology: ?*const boundary_topology.State = null,
    boundary_face_area_m2: []const f64 = &.{},
    boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal: []const f64 = &.{},
    boundary_layer_volume_m3: []const f64 = &.{},
    boundary_layer_midpoint_depth_m: []const f64 = &.{},
    boundary_layer_bottom_depth_m: []const f64 = &.{},
    /// Solver-private scratch. The caller copies it only after convergence.
    artificial_drainage_outflow_m3_per_step: ?[]f64 = null,
};

pub const Options = struct {
    /// Runtime NPH option. It is a convergence ceiling, not a sub-hour loop.
    max_iterations: u16,
    absolute_tolerance_m3: f64 = 1e-12,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Conservative whole-hour dual-domain Richards solve. Every face evaluates
/// one common nonlinear trial state, while an order-independent conservative
/// target accumulator enforces aggregate donor and receiver bounds. The
/// Newton/Picard hybrid replaces repetition of the full sub-hour model. Grid
/// storage and output fluxes are committed only after convergence.
pub fn solve(
    allocator: std.mem.Allocator,
    grid: *grid_module.GridState,
    faces: []const Face,
    properties: Properties,
    micropore_face_flux_m3_per_step: []f64,
    macropore_face_flux_m3_per_step: []f64,
    options: Options,
) !Result {
    try validateInputs(grid, faces, properties, micropore_face_flux_m3_per_step, macropore_face_flux_m3_per_step, options);
    const cells = grid.layer_count;
    const components = try std.math.mul(usize, cells, 2);
    const base = try allocator.alloc(f64, components);
    defer allocator.free(base);
    @memcpy(base[0..cells], grid.matrix_liquid_water_m3);
    @memcpy(base[cells..], grid.macropore_liquid_water_m3);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, components);
    defer allocator.free(residual);
    const target = try allocator.alloc(f64, components);
    defer allocator.free(target);
    const probe = try allocator.alloc(f64, components);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, components);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, components);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, components);
    defer allocator.free(candidate_residual);
    const best_candidate = try allocator.alloc(f64, components);
    defer allocator.free(best_candidate);
    const previous_state = try allocator.alloc(f64, components);
    defer allocator.free(previous_state);
    const previous_residual = try allocator.alloc(f64, components);
    defer allocator.free(previous_residual);
    const previous_previous_state = try allocator.alloc(f64, components);
    defer allocator.free(previous_previous_state);
    const previous_previous_residual = try allocator.alloc(f64, components);
    defer allocator.free(previous_previous_residual);
    const scratch = try allocator.alloc(f64, components);
    defer allocator.free(scratch);
    const trial_micro_flux = try allocator.alloc(f64, faces.len);
    defer allocator.free(trial_micro_flux);
    const trial_macro_flux = try allocator.alloc(f64, faces.len);
    defer allocator.free(trial_macro_flux);
    const jacobian = try allocator.alloc(f64, try std.math.mul(usize, components, components));
    defer allocator.free(jacobian);
    const newton_delta = try allocator.alloc(f64, components);
    defer allocator.free(newton_delta);
    const reduced_jacobian = try allocator.alloc(f64, try std.math.mul(usize, components, components));
    defer allocator.free(reduced_jacobian);
    const reduced_right_hand_side = try allocator.alloc(f64, components);
    defer allocator.free(reduced_right_hand_side);
    const reduced_index_by_component = try allocator.alloc(usize, components);
    defer allocator.free(reduced_index_by_component);
    const micropore_parent = try allocator.alloc(usize, cells);
    defer allocator.free(micropore_parent);
    const macropore_parent = try allocator.alloc(usize, cells);
    defer allocator.free(macropore_parent);
    const micropore_component_size = try allocator.alloc(usize, cells);
    defer allocator.free(micropore_component_size);
    const macropore_component_size = try allocator.alloc(usize, cells);
    defer allocator.free(macropore_component_size);
    for (0..cells) |cell| {
        micropore_parent[cell] = cell;
        macropore_parent[cell] = cell;
        micropore_component_size[cell] = 1;
        macropore_component_size[cell] = 1;
    }

    var newton_steps: u16 = 0;
    var dense_newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var history_count: u8 = 0;
    var final_aitken_norm: ?f64 = null;
    var final_component_secant_norm: ?f64 = null;
    var final_anderson_depth_two_norm: ?f64 = null;
    var final_anderson_depth_one_norm: ?f64 = null;
    var final_dense_candidate_norm: ?f64 = null;
    var final_dense_linear_solved = false;
    var iteration: u16 = 0;
    const active_properties = properties;
    while (iteration < options.max_iterations) : (iteration += 1) {
        // Every Newton/Picard iteration solves the actual whole-hour WATSUB
        // residual. Intermediate states are never committed.
        residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux) catch |initial_error| {
            var restored = false;
            var restoration: u8 = 0;
            while (restoration < 12) : (restoration += 1) {
                for (current, base) |*value, baseline| value.* = 0.5 * (value.* + baseline);
                if (residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                    restored = true;
                    break;
                } else |_| {}
            }
            if (!restored) {
                @memcpy(current, base);
                residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux) catch return initial_error;
            }
        };
        const current_norm = try scaledNorm(current, residual, options);
        if (current_norm <= 1) {
            try commit(grid, properties, current);
            // Re-evaluate at the converged state so published FLWM/FLWHM
            // exactly correspond to the committed water state.
            try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, micropore_face_flux_m3_per_step, macropore_face_flux_m3_per_step);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = current_norm };
        }
        var accepted_newton = false;
        if (history_count > 0 and iteration + 1 == options.max_iterations) {
            var valid_secant_candidate = true;
            var secant_component_count: usize = 0;
            for (current, residual, previous_state, previous_residual, candidate) |value, value_residual, old_value, old_residual, *next| {
                const residual_change = value_residual - old_residual;
                if (!std.math.isFinite(residual_change) or @abs(residual_change) <= std.math.floatEps(f64) * @max(1.0, @abs(value_residual))) {
                    next.* = value;
                    continue;
                }
                next.* = value - value_residual * (value - old_value) / residual_change;
                if (!std.math.isFinite(next.*) or next.* < -1.0e-12) {
                    valid_secant_candidate = false;
                    break;
                }
                next.* = @max(0, next.*);
                secant_component_count += 1;
            }
            if (valid_secant_candidate and secant_component_count > 0) {
                if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                    const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                    final_component_secant_norm = candidate_norm;
                    if (candidate_norm <= 1) {
                        rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        accepted_newton = true;
                    }
                } else |_| {}
            }
            if (accepted_newton) continue;
            if (history_count > 1) {
                var valid_aitken_candidate = true;
                var accelerated_component_count: usize = 0;
                for (current, previous_state, previous_previous_state, candidate) |value, old_value, older_value, *next| {
                    const first_difference = value - old_value;
                    const second_difference = value - 2.0 * old_value + older_value;
                    if (!std.math.isFinite(second_difference) or @abs(second_difference) <= std.math.floatEps(f64) * @max(1.0, @abs(value))) {
                        next.* = value;
                        continue;
                    }
                    next.* = value - first_difference * first_difference / second_difference;
                    if (!std.math.isFinite(next.*) or next.* < -1.0e-12) {
                        valid_aitken_candidate = false;
                        break;
                    }
                    next.* = @max(0, next.*);
                    accelerated_component_count += 1;
                }
                if (valid_aitken_candidate and accelerated_component_count > 0) {
                    if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                        const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                        final_aitken_norm = candidate_norm;
                        if (candidate_norm <= 1) {
                            rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                        }
                    } else |_| {}
                }
                if (accepted_newton) continue;
                var gram_00: f64 = 0;
                var gram_01: f64 = 0;
                var gram_11: f64 = 0;
                var rhs_0: f64 = 0;
                var rhs_1: f64 = 0;
                for (residual, previous_residual, previous_previous_residual) |current_value, previous_value, older_value| {
                    const difference_0 = current_value - previous_value;
                    const difference_1 = previous_value - older_value;
                    gram_00 += difference_0 * difference_0;
                    gram_01 += difference_0 * difference_1;
                    gram_11 += difference_1 * difference_1;
                    rhs_0 += difference_0 * current_value;
                    rhs_1 += difference_1 * current_value;
                }
                const determinant = gram_00 * gram_11 - gram_01 * gram_01;
                if (std.math.isFinite(determinant) and @abs(determinant) > std.math.floatEps(f64) * @max(1.0, gram_00 * gram_11)) {
                    const mixing_0 = (rhs_0 * gram_11 - rhs_1 * gram_01) / determinant;
                    const mixing_1 = (gram_00 * rhs_1 - gram_01 * rhs_0) / determinant;
                    var valid_candidate = std.math.isFinite(mixing_0) and std.math.isFinite(mixing_1);
                    for (current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, candidate) |value, value_residual, old_value, old_residual, older_value, older_residual, *next| {
                        const fixed_point = value + value_residual;
                        const old_fixed_point = old_value + old_residual;
                        const older_fixed_point = older_value + older_residual;
                        next.* = fixed_point - mixing_0 * (fixed_point - old_fixed_point) - mixing_1 * (old_fixed_point - older_fixed_point);
                        if (!std.math.isFinite(next.*) or next.* < -1.0e-12) valid_candidate = false else next.* = @max(0, next.*);
                    }
                    if (valid_candidate) {
                        if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                            const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                            final_anderson_depth_two_norm = candidate_norm;
                            if (candidate_norm <= 1) {
                                rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                                @memcpy(current, candidate);
                                newton_steps += 1;
                                accepted_newton = true;
                            }
                        } else |_| {}
                    }
                }
            }
            if (accepted_newton) continue;
            var numerator: f64 = 0;
            var denominator: f64 = 0;
            for (residual, previous_residual) |current_value, previous_value| {
                const change = current_value - previous_value;
                numerator += current_value * change;
                denominator += change * change;
            }
            if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                const mixing = numerator / denominator;
                var valid_candidate = std.math.isFinite(mixing);
                for (current, residual, previous_state, previous_residual, candidate) |value, value_residual, old_value, old_residual, *next| {
                    const fixed_point = value + value_residual;
                    const old_fixed_point = old_value + old_residual;
                    next.* = fixed_point - mixing * (fixed_point - old_fixed_point);
                    if (!std.math.isFinite(next.*) or next.* < -1.0e-12) valid_candidate = false else next.* = @max(0, next.*);
                }
                if (valid_candidate) {
                    if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                        const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                        final_anderson_depth_one_norm = candidate_norm;
                        if (candidate_norm <= 1) {
                            rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                        }
                    } else |_| {}
                }
            }
        }
        if (accepted_newton) continue;
        // Let bounded Picard/directional iterations settle the first
        // donor-limited active set, then spend the remaining NPH budget on
        // the spatially coupled Richards Jacobian. Large geospatial cells
        // expose wetting fronts for which a diagonal correction cannot
        // propagate information through the column quickly enough.
        var dense_jacobian_valid =
            components <= 256 and iteration + 8 >= options.max_iterations;
        if (dense_jacobian_valid) for (0..components) |column| {
            const column_cell = if (column < cells) column else column - cells;
            const pore_capacity = if (column < cells) grid.matrix_pore_capacity_m3[column_cell] else grid.macropore_pore_capacity_m3[column_cell];
            const ice_volume = if (column < cells) grid.matrix_ice_water_m3[column_cell] else grid.macropore_ice_water_m3[column_cell];
            const available_capacity = @max(0.0, pore_capacity - ice_volume);
            if (available_capacity <= 1.0e-18) {
                for (0..components) |row| jacobian[row * components + column] = if (row == column) 1 else 0;
                continue;
            }
            const nominal_perturbation = @max(
                1.0e-6,
                std.math.sqrt(std.math.floatEps(f64)) *
                    @abs(current[column]),
            );
            const upward_room = @max(0.0, available_capacity - current[column]);
            const downward_room = @max(0.0, current[column]);
            const central_perturbation = @min(
                nominal_perturbation,
                @min(0.5 * upward_room, 0.5 * downward_room),
            );
            if (central_perturbation > 1.0e-20) {
                @memcpy(probe, current);
                probe[column] += central_perturbation;
                if (residualAt(
                    grid,
                    faces,
                    active_properties,
                    base,
                    probe,
                    target,
                    probe_residual,
                    scratch,
                    trial_micro_flux,
                    trial_macro_flux,
                )) |_| {
                    @memcpy(candidate, current);
                    candidate[column] -= central_perturbation;
                    if (residualAt(
                        grid,
                        faces,
                        active_properties,
                        base,
                        candidate,
                        target,
                        candidate_residual,
                        scratch,
                        trial_micro_flux,
                        trial_macro_flux,
                    )) |_| {
                        for (0..components) |row|
                            jacobian[row * components + column] =
                                (probe_residual[row] -
                                    candidate_residual[row]) /
                                (2.0 * central_perturbation);
                        continue;
                    } else |_| {}
                } else |_| {}
            }
            const use_upward = upward_room >= downward_room;
            const perturbation =
                @min(
                    nominal_perturbation,
                    0.5 * @max(upward_room, downward_room),
                );
            if (perturbation <= 1.0e-20) {
                for (0..components) |row| jacobian[row * components + column] = if (row == column) 1 else 0;
                continue;
            }
            var signed_perturbation = if (use_upward) perturbation else -perturbation;
            var probe_attempt: u8 = 0;
            while (true) : (probe_attempt += 1) {
                @memcpy(probe, current);
                probe[column] += signed_perturbation;
                if (residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| break else |_| {
                    if (probe_attempt >= 15) {
                        dense_jacobian_valid = false;
                        break;
                    }
                    signed_perturbation *= 0.5;
                }
            }
            if (!dense_jacobian_valid) break;
            for (0..components) |row| jacobian[row * components + column] = (probe_residual[row] - residual[row]) / signed_perturbation;
        };
        if (dense_jacobian_valid) {
            for (residual, newton_delta) |value, *right_hand_side| right_hand_side.* = -value;
            // Determine components from current numerical couplings rather
            // than geometric faces: donor limiting or zero conductance can
            // split one geometric component into independent Jacobian blocks.
            for (0..cells) |cell| {
                micropore_parent[cell] = cell;
                macropore_parent[cell] = cell;
                micropore_component_size[cell] = 1;
                macropore_component_size[cell] = 1;
            }
            for (0..cells) |first| {
                for (first + 1..cells) |second| {
                    if (@abs(jacobian[first * components + second]) > 1.0e-20 or @abs(jacobian[second * components + first]) > 1.0e-20) unionComponents(micropore_parent, micropore_component_size, first, second);
                    const macro_first = cells + first;
                    const macro_second = cells + second;
                    if (@abs(jacobian[macro_first * components + macro_second]) > 1.0e-20 or @abs(jacobian[macro_second * components + macro_first]) > 1.0e-20) unionComponents(macropore_parent, macropore_component_size, first, second);
                }
            }
            // Release a component's total-storage coordinate only when the
            // current Jacobian actually contains an external source/sink.
            // A configured boundary can be inactive at a donor or water-table
            // clip and must then remain exactly mass-conserving.
            releaseIndependentStorageCoordinates(micropore_parent, micropore_component_size);
            releaseIndependentStorageCoordinates(macropore_parent, macropore_component_size);
            // Re-evaluate closed-component columns in the coordinates that the
            // linear solve actually uses.  Perturbing a store independently
            // and subtracting an independently sampled anchor column takes
            // both samples off the conserved manifold; near nonlinear
            // retention curves that difference is needlessly ill-conditioned.
            // A direct member-to-anchor transfer preserves total water during
            // the residual evaluation and gives the exact reduced derivative.
            for (0..components) |column| {
                const domain_offset = if (column < cells) @as(usize, 0) else cells;
                const column_cell = column - domain_offset;
                const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
                const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
                const root = componentRoot(parent, column_cell);
                if (component_size[root] <= 1) continue;
                const anchor = domain_offset + root;
                if (column == anchor) {
                    for (0..components) |row| jacobian[row * components + column] = 0;
                    continue;
                }
                const column_capacity = if (domain_offset == 0) grid.matrix_pore_capacity_m3[column_cell] - grid.matrix_ice_water_m3[column_cell] else grid.macropore_pore_capacity_m3[column_cell] - grid.macropore_ice_water_m3[column_cell];
                const anchor_capacity = if (domain_offset == 0) grid.matrix_pore_capacity_m3[root] - grid.matrix_ice_water_m3[root] else grid.macropore_pore_capacity_m3[root] - grid.macropore_ice_water_m3[root];
                const positive_room = @min(@max(0, column_capacity - current[column]), @max(0, current[anchor]));
                const negative_room = @min(@max(0, current[column]), @max(0, anchor_capacity - current[anchor]));
                const nominal = @max(
                    1.0e-6,
                    std.math.sqrt(std.math.floatEps(f64)) * @max(
                        @abs(current[column]),
                        @abs(current[anchor]),
                    ),
                );
                var signed_perturbation = if (positive_room >= negative_room) @min(nominal, 0.5 * positive_room) else -@min(nominal, 0.5 * negative_room);
                if (@abs(signed_perturbation) <= 1.0e-20) {
                    dense_jacobian_valid = false;
                    break;
                }
                var probe_attempt: u8 = 0;
                while (true) : (probe_attempt += 1) {
                    @memcpy(probe, current);
                    probe[column] += signed_perturbation;
                    probe[anchor] -= signed_perturbation;
                    if (residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| break else |_| {
                        if (probe_attempt >= 15) {
                            dense_jacobian_valid = false;
                            break;
                        }
                        signed_perturbation *= 0.5;
                    }
                }
                if (!dense_jacobian_valid) break;
                for (0..components) |row| jacobian[row * components + column] = (probe_residual[row] - residual[row]) / signed_perturbation;
            }
            if (!dense_jacobian_valid) continue;
            if (solveConservedNewtonSystem(jacobian, residual, newton_delta, reduced_jacobian, reduced_right_hand_side, reduced_index_by_component, components, cells, micropore_parent, micropore_component_size, macropore_parent, macropore_component_size)) {
                if (iteration + 1 == options.max_iterations) final_dense_linear_solved = true;
                var line_fraction: f64 = 1;
                var line_search: u8 = 0;
                var best_candidate_norm = current_norm;
                while (line_search < 6) : (line_search += 1) {
                    const valid_candidate = projectStorageStep(grid, current, newton_delta, line_fraction, candidate);
                    if (valid_candidate) {
                        if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                            const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                            if (iteration + 1 == options.max_iterations and (final_dense_candidate_norm == null or candidate_norm < final_dense_candidate_norm.?)) final_dense_candidate_norm = candidate_norm;
                            if (candidate_norm < best_candidate_norm) {
                                best_candidate_norm = candidate_norm;
                                @memcpy(best_candidate, candidate);
                                accepted_newton = true;
                            }
                        } else |_| {}
                    }
                    line_fraction *= 0.5;
                }
                if (accepted_newton) {
                    rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                    @memcpy(current, best_candidate);
                    newton_steps += 1;
                    dense_newton_steps += 1;
                }
            }
            if (!accepted_newton) {
                var best_trust_region_norm = current_norm;
                var damping_index: u8 = 0;
                while (damping_index < 11) : (damping_index += 1) {
                    const damping = std.math.pow(f64, 10.0, -8.0 + @as(f64, @floatFromInt(damping_index)));
                    // Solve directly in independent mass-conserving
                    // coordinates. Solving the singular full storage system
                    // and projecting afterward changes the Newton optimum and
                    // stalls when an active pore bound splits a component.
                    if (!solveConservedTrustRegionSystem(jacobian, residual, current, options, newton_delta, reduced_jacobian, reduced_right_hand_side, reduced_index_by_component, components, cells, micropore_parent, micropore_component_size, macropore_parent, macropore_component_size, damping)) continue;
                    var line_fraction: f64 = 1;
                    var line_search: u8 = 0;
                    while (line_search < 8) : (line_search += 1) {
                        const valid_candidate = projectStorageStep(grid, current, newton_delta, line_fraction, candidate);
                        if (valid_candidate) {
                            if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                                const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                                if (candidate_norm < best_trust_region_norm) {
                                    best_trust_region_norm = candidate_norm;
                                    @memcpy(best_candidate, candidate);
                                    accepted_newton = true;
                                }
                            } else |_| {}
                        }
                        line_fraction *= 0.5;
                    }
                }
                if (accepted_newton) {
                    rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                    @memcpy(current, best_candidate);
                    newton_steps += 1;
                    dense_newton_steps += 1;
                }
            }
        }
        if (accepted_newton) continue;
        // If the fully coupled matrix is ill-conditioned because a nearly
        // empty macropore coordinate is on a bound, solve the complete
        // micropore Richards column as one Newton block. This is still one
        // NPH iteration and retains the full residual (including dual-domain
        // exchange) in its line-search merit function.
        if (try domainBlockNewton(
            grid,
            faces,
            active_properties,
            options,
            base,
            current,
            target,
            residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
            0,
            false,
            best_candidate,
            candidate,
            candidate_residual,
            probe,
            probe_residual,
            reduced_jacobian,
            reduced_right_hand_side,
        )) {
            rememberIteration(
                current,
                residual,
                previous_state,
                previous_residual,
                previous_previous_state,
                previous_previous_residual,
                &history_count,
            );
            @memcpy(current, best_candidate);
            newton_steps += 1;
            continue;
        }
        // When the complete domain step is rejected at a semismooth donor
        // transition, resolve the limiting cell and two incident face rings
        // as one coupled Newton block. This consumes the current NPH
        // iteration; it is not deferred to an uncounted post-loop sweep.
        var limiting_component: usize = 0;
        var limiting_component_norm: f64 = -1;
        for (current, residual, 0..) |state, difference, component| {
            const norm =
                scaledComponentResidual(state, difference, options);
            if (norm > limiting_component_norm) {
                limiting_component_norm = norm;
                limiting_component = component;
            }
        }
        @memcpy(best_candidate, current);
        @memcpy(newton_delta, residual);
        var spatial_block_changed = try spatialDomainCellNewton(
            grid,
            faces,
            active_properties,
            base,
            current,
            target,
            residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
            limiting_component,
            candidate,
            candidate_residual,
            probe,
            probe_residual,
            reduced_index_by_component,
            reduced_jacobian,
            reduced_right_hand_side,
            options,
            true,
            false,
        );
        if (!spatial_block_changed and limiting_component < cells)
            spatial_block_changed = try spatialDomainCellNewton(
                grid,
                faces,
                active_properties,
                base,
                current,
                target,
                residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
                limiting_component,
                candidate,
                candidate_residual,
                probe,
                probe_residual,
                reduced_index_by_component,
                reduced_jacobian,
                reduced_right_hand_side,
                options,
                false,
                false,
            );
        if (spatial_block_changed) {
            @memcpy(previous_state, best_candidate);
            @memcpy(previous_residual, newton_delta);
            history_count = 1;
            newton_steps += 1;
            continue;
        }
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                var numerator: f64 = 0;
                var denominator: f64 = 0;
                for (residual, probe_residual) |base_residual, sampled_residual| {
                    const derivative = (sampled_residual - base_residual) / options.directional_probe_fraction;
                    numerator += base_residual * derivative;
                    denominator += derivative * derivative;
                }
                // Diagonal secant Newton step. Each runtime water store can
                // have a very different stiffness near saturation; a single
                // global line fraction leaves the slowest component far from
                // convergence within NPH.
                var diagonal_candidate_valid = true;
                for (current, residual, probe_residual, candidate) |value, base_residual, sampled_residual, *next| {
                    const derivative = (sampled_residual - base_residual) / options.directional_probe_fraction;
                    const fraction = if (std.math.isFinite(derivative) and @abs(derivative) > std.math.floatEps(f64)) std.math.clamp(-base_residual / derivative, options.minimum_newton_fraction, options.maximum_newton_fraction) else options.picard_relaxation;
                    next.* = value + fraction * base_residual;
                    if (!std.math.isFinite(next.*) or next.* < -1e-12) diagonal_candidate_valid = false else next.* = @max(0.0, next.*);
                }
                if (diagonal_candidate_valid) {
                    if (iteration + 1 == options.max_iterations) {
                        @memcpy(probe, candidate);
                        var best_norm = current_norm;
                        var scale_index: u8 = 0;
                        while (scale_index < 29) : (scale_index += 1) {
                            const scale = 0.5 + 0.125 * @as(f64, @floatFromInt(scale_index));
                            var scaled_candidate_valid = true;
                            for (current, probe, candidate) |value, diagonal_value, *next| {
                                next.* = value + scale * (diagonal_value - value);
                                if (!std.math.isFinite(next.*) or next.* < -1.0e-12) scaled_candidate_valid = false else next.* = @max(0, next.*);
                            }
                            if (scaled_candidate_valid) {
                                if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                                    const candidate_norm = try scaledNorm(candidate, candidate_residual, options);
                                    if (candidate_norm < best_norm) {
                                        best_norm = candidate_norm;
                                        @memcpy(best_candidate, candidate);
                                        accepted_newton = true;
                                    }
                                    var correction_index: u8 = 1;
                                    while (correction_index <= 4) : (correction_index += 1) {
                                        const correction_fraction = 0.25 * @as(f64, @floatFromInt(correction_index));
                                        var corrected_candidate_valid = true;
                                        for (candidate, candidate_residual, previous_previous_state) |value, value_residual, *corrected| {
                                            corrected.* = value + correction_fraction * value_residual;
                                            if (!std.math.isFinite(corrected.*) or corrected.* < -1.0e-12) corrected_candidate_valid = false else corrected.* = @max(0, corrected.*);
                                        }
                                        if (corrected_candidate_valid) {
                                            if (residualAt(grid, faces, active_properties, base, previous_previous_state, target, previous_previous_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                                                const corrected_norm = try scaledNorm(previous_previous_state, previous_previous_residual, options);
                                                if (corrected_norm < best_norm) {
                                                    best_norm = corrected_norm;
                                                    @memcpy(best_candidate, previous_previous_state);
                                                    accepted_newton = true;
                                                }
                                            } else |_| {}
                                        }
                                    }
                                } else |_| {}
                            }
                        }
                        if (accepted_newton) {
                            rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                            @memcpy(current, best_candidate);
                            newton_steps += 1;
                        }
                    } else if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                        if (try scaledNorm(candidate, candidate_residual, options) < current_norm) {
                            rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                        }
                    } else |_| {}
                }
                if (accepted_newton) continue;
                if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                    const fraction = std.math.clamp(-numerator / denominator, options.minimum_newton_fraction, options.maximum_newton_fraction);
                    if (addDirection(current, residual, fraction, candidate)) |_| {
                        if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                            if (try scaledNorm(candidate, candidate_residual, options) < current_norm) {
                                rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
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
        rememberIteration(current, residual, previous_state, previous_residual, previous_previous_state, previous_previous_residual, &history_count);
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux);
    const final_norm = try scaledNorm(current, residual, options);
    if (final_norm <= 1) {
        try commit(grid, properties, current);
        try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, micropore_face_flux_m3_per_step, macropore_face_flux_m3_per_step);
        return .{ .iterations = options.max_iterations, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = final_norm };
    }
    // Final bounded residual line search belongs to the last Newton iteration,
    // and is only evaluated after its primary Jacobian/secant candidates fail.
    // Residual directions conserve each closed pore domain by construction.
    var polished_norm = final_norm;
    var polish_found = false;
    var polish_index: u16 = 1;
    while (polish_index <= 256) : (polish_index += 1) {
        // Resolve the one-dimensional bounded Newton line accurately enough
        // that the answer is governed by the requested water tolerance, not
        // by a coarse line-search grid. This remains one Newton correction.
        const fraction = 0.03125 * @as(f64, @floatFromInt(polish_index));
        var valid = true;
        for (current, residual, candidate) |value, difference, *next| {
            next.* = value + fraction * difference;
            if (!std.math.isFinite(next.*) or next.* < -1e-12) valid = false else next.* = @max(0, next.*);
        }
        if (!valid) continue;
        if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
            const norm = try scaledNorm(candidate, candidate_residual, options);
            if (norm < polished_norm) {
                polished_norm = norm;
                @memcpy(best_candidate, candidate);
                polish_found = true;
            }
        } else |_| {}
    }
    // Coordinate Newton in conservative face-transfer space. Each coordinate
    // adds δ to one endpoint and removes δ from the other, so the Jacobian has
    // no closed-domain storage nullspace and no post-hoc mass projection.
    var coordinate_correction: usize = 0;
    const coordinate_correction_limit = try std.math.mul(usize, components, 4);
    while (polished_norm > 1 and coordinate_correction < coordinate_correction_limit) : (coordinate_correction += 1) {
        if (polish_found) {
            @memcpy(current, best_candidate);
            try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux);
        }
        var maximum_index: usize = 0;
        for (residual, 0..) |value, index| if (@abs(value) > @abs(residual[maximum_index])) {
            maximum_index = index;
        };
        const domain_offset = if (maximum_index < cells) @as(usize, 0) else cells;
        const cell = maximum_index - domain_offset;
        const epsilon = std.math.cbrt(std.math.floatEps(f64)) * @max(1e-6, @abs(current[maximum_index]));
        const positive_room = if (domain_offset == 0) grid.matrix_pore_capacity_m3[cell] - grid.matrix_ice_water_m3[cell] - current[maximum_index] else grid.macropore_pore_capacity_m3[cell] - grid.macropore_ice_water_m3[cell] - current[maximum_index];
        var accepted_coordinate = false;
        // The residual is target(state) - state, so a single storage
        // coordinate is nonsingular even for a closed domain: its identity
        // term restores the exact total inventory. This coordinate is also
        // essential when an active boundary source/sink changes domain mass;
        // a pairwise-only correction cannot span that mode.
        const negative_room = current[maximum_index];
        const signed_storage_probe = if (positive_room >= negative_room) @min(epsilon, 0.5 * @max(0, positive_room)) else -@min(epsilon, 0.5 * negative_room);
        if (@abs(signed_storage_probe) > 1e-20) {
            @memcpy(probe, current);
            probe[maximum_index] += signed_storage_probe;
            if (residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                var gradient: f64 = 0;
                var curvature: f64 = 0;
                for (residual, probe_residual) |value_residual, sampled_residual| {
                    const derivative = (sampled_residual - value_residual) / signed_storage_probe;
                    gradient += value_residual * derivative;
                    curvature += derivative * derivative;
                }
                if (std.math.isFinite(gradient) and std.math.isFinite(curvature) and curvature > std.math.floatEps(f64)) {
                    const correction = -gradient / curvature;
                    var line_fraction: f64 = 1;
                    var line: u8 = 0;
                    while (line < 14) : (line += 1) {
                        @memcpy(candidate, current);
                        candidate[maximum_index] += line_fraction * correction;
                        const capacity = current[maximum_index] + @max(0, positive_room);
                        if (candidate[maximum_index] >= -1e-12 and candidate[maximum_index] <= capacity + 1e-12) {
                            candidate[maximum_index] = std.math.clamp(candidate[maximum_index], 0, capacity);
                            if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                                const norm = try scaledNorm(candidate, candidate_residual, options);
                                if (norm < polished_norm) {
                                    polished_norm = norm;
                                    @memcpy(best_candidate, candidate);
                                    accepted_coordinate = true;
                                }
                            } else |_| {}
                        }
                        line_fraction *= 0.5;
                    }
                }
            } else |_| {}
        }
        for (faces) |face| {
            const neighbor = if (face.source_cell == cell)
                face.destination_cell
            else if (face.destination_cell == cell)
                face.source_cell
            else
                continue;
            const neighbor_index = domain_offset + neighbor;
            const neighbor_positive_room = if (domain_offset == 0) grid.matrix_pore_capacity_m3[neighbor] - grid.matrix_ice_water_m3[neighbor] - current[neighbor_index] else grid.macropore_pore_capacity_m3[neighbor] - grid.macropore_ice_water_m3[neighbor] - current[neighbor_index];
            const forward_room = @min(@max(0, positive_room), current[neighbor_index]);
            const backward_room = @min(current[maximum_index], @max(0, neighbor_positive_room));
            const signed_probe = if (forward_room >= backward_room) @min(epsilon, 0.5 * forward_room) else -@min(epsilon, 0.5 * backward_room);
            if (@abs(signed_probe) <= 1e-20) continue;
            @memcpy(probe, current);
            probe[maximum_index] += signed_probe;
            probe[neighbor_index] -= signed_probe;
            residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux) catch continue;
            var gradient: f64 = 0;
            var curvature: f64 = 0;
            for (residual, probe_residual) |value_residual, sampled_residual| {
                const derivative = (sampled_residual - value_residual) / signed_probe;
                gradient += value_residual * derivative;
                curvature += derivative * derivative;
            }
            if (!std.math.isFinite(gradient) or !std.math.isFinite(curvature) or curvature <= std.math.floatEps(f64)) continue;
            const transfer = -gradient / curvature;
            var line_fraction: f64 = 1;
            var line: u8 = 0;
            while (line < 10) : (line += 1) {
                @memcpy(candidate, current);
                candidate[maximum_index] += line_fraction * transfer;
                candidate[neighbor_index] -= line_fraction * transfer;
                if (candidate[maximum_index] >= -1e-12 and candidate[neighbor_index] >= -1e-12) {
                    candidate[maximum_index] = @max(0, candidate[maximum_index]);
                    candidate[neighbor_index] = @max(0, candidate[neighbor_index]);
                    if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                        const norm = try scaledNorm(candidate, candidate_residual, options);
                        if (norm < polished_norm) {
                            polished_norm = norm;
                            @memcpy(best_candidate, candidate);
                            accepted_coordinate = true;
                        }
                    } else |_| {}
                }
                line_fraction *= 0.5;
            }
        }
        if (!accepted_coordinate) break;
        polish_found = true;
    }
    // Complete the final semismooth Newton correction with nonlinear
    // Gauss-Seidel coordinates. A donor limiter can make reducing one
    // component temporarily increase its neighbor, so requiring every inner
    // coordinate to reduce the global infinity norm stalls at wetting fronts.
    // These are inner coordinates of one Newton solve, analogous to solving
    // the dense linear system, and do not consume additional NPH cycles.
    if (polish_found) {
        @memcpy(current, best_candidate);
        try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux);
    }
    var gauss_seidel_sweep: u8 = 0;
    while (polished_norm > 1 and gauss_seidel_sweep < 32) : (gauss_seidel_sweep += 1) {
        var sweep_changed = false;
        for (0..components) |component| {
            const cell = if (component < cells) component else component - cells;
            const component_tolerance = options.absolute_tolerance_m3 +
                options.relative_tolerance * @max(1.0, @abs(current[component]));
            var spatial_block_changed = false;
            if (@abs(residual[component]) > component_tolerance) {
                spatial_block_changed = try spatialDomainCellNewton(
                    grid,
                    faces,
                    active_properties,
                    base,
                    current,
                    target,
                    residual,
                    scratch,
                    trial_micro_flux,
                    trial_macro_flux,
                    component,
                    candidate,
                    candidate_residual,
                    probe,
                    probe_residual,
                    reduced_index_by_component,
                    reduced_jacobian,
                    reduced_right_hand_side,
                    options,
                    true,
                    false,
                );
                // A nearly empty macropore coordinate can make the coupled
                // star singular even though the matrix Richards block is
                // well-conditioned. Preserve the coupled attempt as the
                // primary solve, then fall back only for the matrix domain.
                if (!spatial_block_changed and component < cells)
                    spatial_block_changed = try spatialDomainCellNewton(
                        grid,
                        faces,
                        active_properties,
                        base,
                        current,
                        target,
                        residual,
                        scratch,
                        trial_micro_flux,
                        trial_macro_flux,
                        component,
                        candidate,
                        candidate_residual,
                        probe,
                        probe_residual,
                        reduced_index_by_component,
                        reduced_jacobian,
                        reduced_right_hand_side,
                        options,
                        false,
                        false,
                    );
            }
            if (spatial_block_changed) {
                sweep_changed = true;
                const norm = try scaledNorm(current, residual, options);
                if (norm < polished_norm) {
                    polished_norm = norm;
                    @memcpy(best_candidate, current);
                    polish_found = true;
                }
                if (polished_norm <= 1) break;
            }
            const capacity = if (component < cells)
                grid.matrix_pore_capacity_m3[cell] - grid.matrix_ice_water_m3[cell]
            else
                grid.macropore_pore_capacity_m3[cell] - grid.macropore_ice_water_m3[cell];
            const upward_room = @max(0.0, capacity - current[component]);
            const downward_room = @max(0.0, current[component]);
            const epsilon = std.math.cbrt(std.math.floatEps(f64)) *
                @max(1e-6, @abs(current[component]));
            if (component < cells and
                properties.dual_domain_exchange_enabled.len != 0 and
                properties.dual_domain_exchange_enabled[cell] and
                grid.macropore_pore_capacity_m3[cell] > 0 and
                try dualDomainCellNewton(
                    grid,
                    faces,
                    active_properties,
                    base,
                    current,
                    target,
                    residual,
                    scratch,
                    trial_micro_flux,
                    trial_macro_flux,
                    cell,
                    candidate,
                    candidate_residual,
                    probe,
                    probe_residual,
                ))
            {
                sweep_changed = true;
                const norm = try scaledNorm(current, residual, options);
                if (norm < polished_norm) {
                    polished_norm = norm;
                    @memcpy(best_candidate, current);
                    polish_found = true;
                }
                if (polished_norm <= 1) break;
                continue;
            }
            const signed_probe = if (upward_room >= downward_room) @min(epsilon, 0.5 * upward_room) else -@min(epsilon, 0.5 * downward_room);
            if (@abs(signed_probe) <= 1e-20) continue;
            @memcpy(probe, current);
            probe[component] += signed_probe;
            residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux) catch continue;
            const derivative = (probe_residual[component] - residual[component]) / signed_probe;
            if (!std.math.isFinite(derivative) or @abs(derivative) <= std.math.floatEps(f64)) continue;
            const correction = -residual[component] / derivative;
            var best_component_residual = @abs(residual[component]);
            var component_candidate_found = false;
            var line_fraction: f64 = 1;
            var line: u8 = 0;
            while (line < 16) : (line += 1) {
                @memcpy(candidate, current);
                candidate[component] = std.math.clamp(current[component] + line_fraction * correction, 0.0, @max(0.0, capacity));
                if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                    if (@abs(candidate_residual[component]) < best_component_residual) {
                        best_component_residual = @abs(candidate_residual[component]);
                        @memcpy(probe, candidate);
                        @memcpy(probe_residual, candidate_residual);
                        component_candidate_found = true;
                    }
                } else |_| {}
                line_fraction *= 0.5;
            }
            if (best_component_residual > component_tolerance) {
                // The donor/receiver MIN/MAX bounds make the residual
                // semismooth. Even a Newton candidate that makes modest
                // progress can remain on the wrong active branch, so retain
                // it as the incumbent while safeguarding the component over
                // its exact physical interval.
                var lower_value: f64 = 0;
                var upper_value: f64 = @max(0.0, capacity);
                var lower_residual: f64 = undefined;
                var upper_residual: f64 = undefined;
                var bracket_valid = true;
                @memcpy(candidate, current);
                candidate[component] = lower_value;
                if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                    lower_residual = candidate_residual[component];
                    if (@abs(lower_residual) < best_component_residual) {
                        best_component_residual = @abs(lower_residual);
                        @memcpy(probe, candidate);
                        @memcpy(probe_residual, candidate_residual);
                        component_candidate_found = true;
                    }
                } else |_| bracket_valid = false;
                @memcpy(candidate, current);
                candidate[component] = upper_value;
                if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                    upper_residual = candidate_residual[component];
                    if (@abs(upper_residual) < best_component_residual) {
                        best_component_residual = @abs(upper_residual);
                        @memcpy(probe, candidate);
                        @memcpy(probe_residual, candidate_residual);
                        component_candidate_found = true;
                    }
                } else |_| bracket_valid = false;
                if (bracket_valid and std.math.signbit(lower_residual) != std.math.signbit(upper_residual)) {
                    var bisection: u8 = 0;
                    while (bisection < 52) : (bisection += 1) {
                        const midpoint = 0.5 * (lower_value + upper_value);
                        @memcpy(candidate, current);
                        candidate[component] = midpoint;
                        residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux) catch break;
                        const midpoint_residual = candidate_residual[component];
                        if (@abs(midpoint_residual) < best_component_residual) {
                            best_component_residual = @abs(midpoint_residual);
                            @memcpy(probe, candidate);
                            @memcpy(probe_residual, candidate_residual);
                            component_candidate_found = true;
                        }
                        if (std.math.signbit(midpoint_residual) == std.math.signbit(lower_residual)) {
                            lower_value = midpoint;
                            lower_residual = midpoint_residual;
                        } else {
                            upper_value = midpoint;
                            upper_residual = midpoint_residual;
                        }
                    }
                }
            }
            if (!component_candidate_found) continue;
            @memcpy(current, probe);
            @memcpy(residual, probe_residual);
            sweep_changed = true;
            const norm = try scaledNorm(current, residual, options);
            if (norm < polished_norm) {
                polished_norm = norm;
                @memcpy(best_candidate, current);
                polish_found = true;
            }
            if (polished_norm <= 1) break;
        }
        if (!sweep_changed) break;
    }
    // Re-linearize at the best semismooth candidate. The Jacobian used by the
    // last outer NPH iteration can represent a different donor/receiver
    // branch, while the coordinate sweeps above may have crossed that branch.
    // A damped least-squares correction on the current active set completes
    // that final Newton solve without adding another full-model or NPH cycle.
    if (polish_found and polished_norm > 1) {
        @memcpy(current, best_candidate);
        try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux);
        var trust_refinement: u8 = 0;
        while (polished_norm > 1 and trust_refinement < 4) : (trust_refinement += 1) {
            var jacobian_valid = true;
            for (0..components) |column| {
                const column_cell = if (column < cells) column else column - cells;
                const capacity =
                    if (column < cells)
                        grid.matrix_pore_capacity_m3[column_cell] -
                            grid.matrix_ice_water_m3[column_cell]
                    else
                        grid.macropore_pore_capacity_m3[column_cell] -
                            grid.macropore_ice_water_m3[column_cell];
                const upward_room = @max(0.0, capacity - current[column]);
                const downward_room = @max(0.0, current[column]);
                const nominal = std.math.sqrt(std.math.floatEps(f64)) *
                    @max(1.0e-6, @abs(current[column]));
                const signed_step =
                    if (upward_room >= downward_room)
                        @min(nominal, 0.5 * upward_room)
                    else
                        -@min(nominal, 0.5 * downward_room);
                if (@abs(signed_step) <= 1.0e-20) {
                    for (0..components) |row|
                        jacobian[row * components + column] =
                            if (row == column) 1 else 0;
                    continue;
                }
                @memcpy(probe, current);
                probe[column] += signed_step;
                if (residualAt(grid, faces, active_properties, base, probe, target, probe_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                    for (0..components) |row|
                        jacobian[row * components + column] =
                            (probe_residual[row] - residual[row]) / signed_step;
                } else |_| {
                    jacobian_valid = false;
                    break;
                }
            }
            if (!jacobian_valid) break;
            var refinement_found = false;
            var refinement_norm = polished_norm;
            var damping_index: u8 = 0;
            while (damping_index < 15) : (damping_index += 1) {
                const damping = std.math.pow(
                    f64,
                    10.0,
                    -12.0 + @as(f64, @floatFromInt(damping_index)),
                );
                if (!solveConservedTrustRegionSystem(
                    jacobian,
                    residual,
                    current,
                    options,
                    newton_delta,
                    reduced_jacobian,
                    reduced_right_hand_side,
                    reduced_index_by_component,
                    components,
                    cells,
                    micropore_parent,
                    micropore_component_size,
                    macropore_parent,
                    macropore_component_size,
                    damping,
                )) continue;
                var fraction: f64 = 1;
                var line: u8 = 0;
                while (line < 16) : (line += 1) {
                    if (projectStorageStep(grid, current, newton_delta, fraction, candidate)) {
                        if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
                            const norm = try scaledNorm(candidate, candidate_residual, options);
                            if (norm < refinement_norm) {
                                refinement_norm = norm;
                                @memcpy(best_candidate, candidate);
                                @memcpy(probe_residual, candidate_residual);
                                refinement_found = true;
                            }
                        } else |_| {}
                    }
                    fraction *= 0.5;
                }
            }
            if (!refinement_found) break;
            polished_norm = refinement_norm;
            @memcpy(current, best_candidate);
            @memcpy(residual, probe_residual);
        }
    }
    if (polish_found and polished_norm <= 1) {
        @memcpy(current, best_candidate);
        try commit(grid, properties, current);
        try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, micropore_face_flux_m3_per_step, macropore_face_flux_m3_per_step);
        return .{ .iterations = options.max_iterations, .newton_raphson_steps = newton_steps + 1, .picard_steps = picard_steps, .maximum_scaled_residual = polished_norm };
    }
    // Report the actual best bounded candidate, not the state that preceded
    // the final conservative correction. This is diagnostic-only: failure
    // remains atomic and neither the grid nor published face fluxes change.
    if (polish_found) {
        @memcpy(current, best_candidate);
        try residualAt(grid, faces, active_properties, base, current, target, residual, scratch, trial_micro_flux, trial_macro_flux);
        polished_norm = try scaledNorm(current, residual, options);
    }
    var final_scaled_index: usize = 0;
    var final_scaled_component_norm: f64 = -1;
    for (current, residual, 0..) |state, difference, index| {
        const norm = scaledComponentResidual(state, difference, options);
        if (norm > final_scaled_component_norm) {
            final_scaled_component_norm = norm;
            final_scaled_index = index;
        }
    }
    const final_scaled_cell =
        if (final_scaled_index < cells)
            final_scaled_index
        else
            final_scaled_index - cells;
    if (polished_norm > 1 and
        active_properties.dual_domain_exchange_enabled.len != 0 and
        active_properties.dual_domain_exchange_enabled[final_scaled_cell] and
        grid.macropore_pore_capacity_m3[final_scaled_cell] > 0 and
        try dualDomainTotalPartitionNewton(
            grid,
            faces,
            active_properties,
            options,
            base,
            current,
            target,
            residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
            final_scaled_cell,
            best_candidate,
            candidate,
            candidate_residual,
            probe,
            probe_residual,
            reduced_jacobian,
            reduced_right_hand_side,
        ))
    {
        @memcpy(current, best_candidate);
        try residualAt(
            grid,
            faces,
            active_properties,
            base,
            current,
            target,
            residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
        );
        polished_norm = try scaledNorm(current, residual, options);
    }
    if (polished_norm > 1) {
        var coupled_limiting_component: usize = 0;
        var coupled_limiting_norm: f64 = -1;
        for (current, residual, 0..) |state, difference, component| {
            const norm =
                scaledComponentResidual(state, difference, options);
            if (norm > coupled_limiting_norm) {
                coupled_limiting_norm = norm;
                coupled_limiting_component = component;
            }
        }
        if (try spatialDomainCellNewton(
            grid,
            faces,
            active_properties,
            base,
            current,
            target,
            residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
            coupled_limiting_component,
            candidate,
            candidate_residual,
            probe,
            probe_residual,
            reduced_index_by_component,
            reduced_jacobian,
            reduced_right_hand_side,
            options,
            true,
            true,
        )) polished_norm = try scaledNorm(current, residual, options);
    }
    if (polished_norm <= 1) {
        try commit(grid, properties, current);
        try residualAt(
            grid,
            faces,
            active_properties,
            base,
            current,
            target,
            residual,
            scratch,
            micropore_face_flux_m3_per_step,
            macropore_face_flux_m3_per_step,
        );
        return .{
            .iterations = options.max_iterations,
            .newton_raphson_steps = newton_steps + 1,
            .picard_steps = picard_steps,
            .maximum_scaled_residual = polished_norm,
        };
    }
    // The last coupled Newton step may be limited by an inactive macropore
    // coordinate even after the matrix column is within a few tolerances.
    // Re-solve that same final linearization in the complete micropore block;
    // this is an inner refinement of the final Newton step, not another NPH
    // iteration or a sub-hour model cycle.
    if (polished_norm > 1 and try domainBlockNewton(
        grid,
        faces,
        active_properties,
        options,
        base,
        current,
        target,
        residual,
        scratch,
        trial_micro_flux,
        trial_macro_flux,
        0,
        true,
        best_candidate,
        candidate,
        candidate_residual,
        probe,
        probe_residual,
        reduced_jacobian,
        reduced_right_hand_side,
    )) {
        @memcpy(current, best_candidate);
        try residualAt(
            grid,
            faces,
            active_properties,
            base,
            current,
            target,
            residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
        );
        polished_norm = try scaledNorm(current, residual, options);
        if (polished_norm <= 1) {
            try commit(grid, properties, current);
            try residualAt(
                grid,
                faces,
                active_properties,
                base,
                current,
                target,
                residual,
                scratch,
                micropore_face_flux_m3_per_step,
                macropore_face_flux_m3_per_step,
            );
            return .{
                .iterations = options.max_iterations,
                .newton_raphson_steps = newton_steps,
                .picard_steps = picard_steps,
                .maximum_scaled_residual = polished_norm,
            };
        }
    }
    // The final semismooth Newton linearization can stop just outside the
    // tolerance when a donor/receiver limiter changes branch. Complete that
    // same final hybrid iteration with one bounded Picard substitution of the
    // limiting storage equation, then audit the complete nonlinear residual.
    // This was previously evaluated only by the failure diagnostic below,
    // even when it had already produced a converged state.
    var limiting_scaled_index: usize = 0;
    var limiting_scaled_residual: f64 = -1;
    for (residual, current, 0..) |value, state, index| {
        const tolerance = options.absolute_tolerance_m3 +
            options.relative_tolerance * @max(1.0, @abs(state));
        const scaled = @abs(value) / tolerance;
        if (scaled > limiting_scaled_residual) {
            limiting_scaled_residual = scaled;
            limiting_scaled_index = index;
        }
    }
    const limiting_scaled_cell =
        if (limiting_scaled_index < cells)
            limiting_scaled_index
        else
            limiting_scaled_index - cells;
    const limiting_scaled_capacity =
        if (limiting_scaled_index < cells)
            grid.matrix_pore_capacity_m3[limiting_scaled_cell] -
                grid.matrix_ice_water_m3[limiting_scaled_cell]
        else
            grid.macropore_pore_capacity_m3[limiting_scaled_cell] -
                grid.macropore_ice_water_m3[limiting_scaled_cell];
    const final_picard_value = std.math.clamp(
        current[limiting_scaled_index] + residual[limiting_scaled_index],
        0,
        @max(0, limiting_scaled_capacity),
    );
    @memcpy(candidate, current);
    candidate[limiting_scaled_index] = final_picard_value;
    var final_picard_norm: f64 = std.math.inf(f64);
    var final_picard_converged = false;
    if (residualAt(
        grid,
        faces,
        active_properties,
        base,
        candidate,
        target,
        candidate_residual,
        scratch,
        trial_micro_flux,
        trial_macro_flux,
    )) |_| {
        final_picard_norm = try scaledNorm(candidate, candidate_residual, options);
        if (final_picard_norm <= 1) {
            @memcpy(best_candidate, candidate);
            final_picard_converged = true;
        } else if (final_picard_value != current[limiting_scaled_index] and
            std.math.signbit(candidate_residual[limiting_scaled_index]) !=
                std.math.signbit(residual[limiting_scaled_index]))
        {
            // The complete physical interval can be non-monotone while the
            // bounded Picard image crosses the local semismooth root. Search
            // only that final limiting-coordinate bracket, and publish it
            // only if the complete nonlinear system satisfies tolerance.
            var first_value = current[limiting_scaled_index];
            var first_residual = residual[limiting_scaled_index];
            var second_value = final_picard_value;
            var bisection: u8 = 0;
            while (bisection < 52) : (bisection += 1) {
                const midpoint = 0.5 * (first_value + second_value);
                @memcpy(candidate, current);
                candidate[limiting_scaled_index] = midpoint;
                residualAt(
                    grid,
                    faces,
                    active_properties,
                    base,
                    candidate,
                    target,
                    candidate_residual,
                    scratch,
                    trial_micro_flux,
                    trial_macro_flux,
                ) catch break;
                const norm = try scaledNorm(candidate, candidate_residual, options);
                if (norm < final_picard_norm) final_picard_norm = norm;
                if (norm <= 1) {
                    @memcpy(best_candidate, candidate);
                    final_picard_converged = true;
                    final_picard_norm = norm;
                    break;
                }
                const midpoint_residual =
                    candidate_residual[limiting_scaled_index];
                if (std.math.signbit(midpoint_residual) ==
                    std.math.signbit(first_residual))
                {
                    first_value = midpoint;
                    first_residual = midpoint_residual;
                } else {
                    second_value = midpoint;
                }
            }
        }
        if (!final_picard_converged and
            active_properties.dual_domain_exchange_enabled.len != 0 and
            active_properties.dual_domain_exchange_enabled[limiting_scaled_cell] and
            grid.macropore_pore_capacity_m3[limiting_scaled_cell] > 0)
        {
            // A scalar storage change can cross a discontinuous exchange
            // donor bound. Resolve that branch on the conservative
            // matrix/macropore coordinate: the limiting store receives
            // `delta` while its paired pore domain loses the same volume.
            const paired_index =
                if (limiting_scaled_index < cells)
                    cells + limiting_scaled_cell
                else
                    limiting_scaled_cell;
            const paired_capacity =
                if (paired_index < cells)
                    grid.matrix_pore_capacity_m3[limiting_scaled_cell] -
                        grid.matrix_ice_water_m3[limiting_scaled_cell]
                else
                    grid.macropore_pore_capacity_m3[limiting_scaled_cell] -
                        grid.macropore_ice_water_m3[limiting_scaled_cell];
            const minimum_pair_delta = @max(
                -current[limiting_scaled_index],
                current[paired_index] - @max(0, paired_capacity),
            );
            const maximum_pair_delta = @min(
                @max(0, limiting_scaled_capacity) -
                    current[limiting_scaled_index],
                current[paired_index],
            );
            const pair_picard_delta = std.math.clamp(
                residual[limiting_scaled_index],
                minimum_pair_delta,
                maximum_pair_delta,
            );
            if (pair_picard_delta != 0) {
                @memcpy(candidate, current);
                candidate[limiting_scaled_index] += pair_picard_delta;
                candidate[paired_index] -= pair_picard_delta;
                if (residualAt(
                    grid,
                    faces,
                    active_properties,
                    base,
                    candidate,
                    target,
                    candidate_residual,
                    scratch,
                    trial_micro_flux,
                    trial_macro_flux,
                )) |_| {
                    const pair_norm =
                        try scaledNorm(candidate, candidate_residual, options);
                    if (pair_norm <= 1) {
                        @memcpy(best_candidate, candidate);
                        final_picard_converged = true;
                        final_picard_norm = pair_norm;
                    } else if (std.math.signbit(
                        candidate_residual[limiting_scaled_index],
                    ) != std.math.signbit(residual[limiting_scaled_index])) {
                        var first_delta: f64 = 0;
                        var first_residual = residual[limiting_scaled_index];
                        var second_delta = pair_picard_delta;
                        var pair_bisection: u8 = 0;
                        while (pair_bisection < 52) : (pair_bisection += 1) {
                            const midpoint_delta =
                                0.5 * (first_delta + second_delta);
                            @memcpy(candidate, current);
                            candidate[limiting_scaled_index] += midpoint_delta;
                            candidate[paired_index] -= midpoint_delta;
                            residualAt(
                                grid,
                                faces,
                                active_properties,
                                base,
                                candidate,
                                target,
                                candidate_residual,
                                scratch,
                                trial_micro_flux,
                                trial_macro_flux,
                            ) catch break;
                            const norm =
                                try scaledNorm(candidate, candidate_residual, options);
                            if (norm <= 1) {
                                @memcpy(best_candidate, candidate);
                                final_picard_converged = true;
                                final_picard_norm = norm;
                                break;
                            }
                            const midpoint_residual =
                                candidate_residual[limiting_scaled_index];
                            if (std.math.signbit(midpoint_residual) ==
                                std.math.signbit(first_residual))
                            {
                                first_delta = midpoint_delta;
                                first_residual = midpoint_residual;
                            } else {
                                second_delta = midpoint_delta;
                            }
                        }
                    }
                } else |_| {}
            }
        }
    } else |_| {}
    if (final_picard_converged) {
        @memcpy(current, best_candidate);
        try commit(grid, properties, current);
        try residualAt(
            grid,
            faces,
            active_properties,
            base,
            current,
            target,
            residual,
            scratch,
            micropore_face_flux_m3_per_step,
            macropore_face_flux_m3_per_step,
        );
        return .{
            .iterations = options.max_iterations,
            .newton_raphson_steps = newton_steps,
            .picard_steps = picard_steps + 1,
            .maximum_scaled_residual = final_picard_norm,
        };
    }
    if (!builtin.is_test) {
        var maximum_residual_index: usize = 0;
        for (residual, 0..) |value, index| if (@abs(value) > @abs(residual[maximum_residual_index])) {
            maximum_residual_index = index;
        };
        std.log.err("soil water Newton-Raphson/Picard did not converge: iterations={d} newton_raphson_steps={d} dense_newton_steps={d} directional_newton_steps={d} picard_steps={d} maximum_scaled_residual={e} component={d} state_m3={e} residual_m3={e}", .{ options.max_iterations, newton_steps, dense_newton_steps, newton_steps - dense_newton_steps, picard_steps, polished_norm, maximum_residual_index, current[maximum_residual_index], residual[maximum_residual_index] });
        std.log.err(
            "soil water limiting scaled component: component={d} state_m3={e} residual_m3={e} scaled_residual={e}",
            .{
                limiting_scaled_index,
                current[limiting_scaled_index],
                residual[limiting_scaled_index],
                scaledComponentResidual(
                    current[limiting_scaled_index],
                    residual[limiting_scaled_index],
                    options,
                ),
            },
        );
        const paired_component =
            if (limiting_scaled_index < cells)
                cells + limiting_scaled_cell
            else
                limiting_scaled_cell;
        std.log.err(
            "soil water limiting paired domain: component={d} state_m3={e} residual_m3={e} scaled_residual={e}",
            .{
                paired_component,
                current[paired_component],
                residual[paired_component],
                scaledComponentResidual(
                    current[paired_component],
                    residual[paired_component],
                    options,
                ),
            },
        );
        std.log.err("soil water final acceleration candidates: component_secant={?} aitken={?} anderson_depth_two={?} anderson_depth_one={?} dense_linear_solved={} dense_best={?} limiting_picard_best={e}", .{ final_component_secant_norm, final_aitken_norm, final_anderson_depth_two_norm, final_anderson_depth_one_norm, final_dense_linear_solved, final_dense_candidate_norm, final_picard_norm });
        const limiting_is_macro = maximum_residual_index >= cells;
        const limiting_cell = if (limiting_is_macro) maximum_residual_index - cells else maximum_residual_index;
        const limiting_capacity = if (limiting_is_macro) grid.macropore_pore_capacity_m3[limiting_cell] - grid.macropore_ice_water_m3[limiting_cell] else grid.matrix_pore_capacity_m3[limiting_cell] - grid.matrix_ice_water_m3[limiting_cell];
        std.log.err("soil water limiting storage: domain={s} layer_cell={d} base_m3={e} target_m3={e} capacity_after_ice_m3={e}", .{ if (limiting_is_macro) "macropore" else "matrix", limiting_cell, base[maximum_residual_index], target[maximum_residual_index], limiting_capacity });
        @memcpy(candidate, current);
        candidate[maximum_residual_index] = std.math.clamp(
            current[maximum_residual_index] + residual[maximum_residual_index],
            0,
            limiting_capacity,
        );
        if (residualAt(grid, faces, active_properties, base, candidate, target, candidate_residual, scratch, trial_micro_flux, trial_macro_flux)) |_| {
            var picard_maximum_index: usize = 0;
            for (candidate_residual, 0..) |value, index| {
                if (@abs(value) > @abs(candidate_residual[picard_maximum_index])) {
                    picard_maximum_index = index;
                }
            }
            std.log.err("soil water limiting exact-Picard probe: changed_component={d} changed_residual_m3={e} maximum_component={d} maximum_residual_m3={e}", .{ maximum_residual_index, candidate_residual[maximum_residual_index], picard_maximum_index, candidate_residual[picard_maximum_index] });
        } else |probe_error| {
            std.log.err("soil water limiting exact-Picard probe failed: error={s}", .{@errorName(probe_error)});
        }
        const limiting_domain_offset = if (limiting_is_macro) cells else 0;
        for (faces) |face| {
            if (face.source_cell != limiting_cell and face.destination_cell != limiting_cell) continue;
            const neighbor_cell = if (face.source_cell == limiting_cell) face.destination_cell else face.source_cell;
            const neighbor_index = limiting_domain_offset + neighbor_cell;
            const neighbor_capacity = if (limiting_is_macro) grid.macropore_pore_capacity_m3[neighbor_cell] - grid.macropore_ice_water_m3[neighbor_cell] else grid.matrix_pore_capacity_m3[neighbor_cell] - grid.matrix_ice_water_m3[neighbor_cell];
            std.log.err("soil water limiting neighbor: domain={s} layer_cell={d} state_m3={e} residual_m3={e} base_m3={e} target_m3={e} capacity_after_ice_m3={e}", .{ if (limiting_is_macro) "macropore" else "matrix", neighbor_cell, current[neighbor_index], residual[neighbor_index], base[neighbor_index], target[neighbor_index], neighbor_capacity });
        }
        for (faces, 0..) |face, face_index| {
            if (face.source_cell != limiting_cell and face.destination_cell != limiting_cell) continue;
            std.log.err("soil water limiting face: face={d} axis={s} source={d} destination={d} matrix_flux_m3={e} macropore_flux_m3={e}", .{ face_index, @tagName(face.axis), face.source_cell, face.destination_cell, trial_micro_flux[face_index], trial_macro_flux[face_index] });
        }
    }
    return error.SoilWaterSolverDidNotConverge;
}

/// Resolves one dual-domain layer in well-scaled physical coordinates:
/// coordinate 0 changes total layer water through the matrix store, while
/// coordinate 1 transfers water conservatively from matrix to macropore.
/// The residual rows are total-water closure and macropore partition closure.
fn dualDomainTotalPartitionNewton(
    grid: *const grid_module.GridState,
    faces: []const Face,
    properties: Properties,
    options: Options,
    base: []const f64,
    current: []const f64,
    target: []f64,
    residual: []const f64,
    scratch: []f64,
    trial_micro_flux: []f64,
    trial_macro_flux: []f64,
    cell: usize,
    accepted_state: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    probe: []f64,
    probe_residual: []f64,
    jacobian_workspace: []f64,
    delta_workspace: []f64,
) !bool {
    const cells = grid.layer_count;
    const matrix_index = cell;
    const macropore_index = cells + cell;
    const matrix_capacity =
        grid.matrix_pore_capacity_m3[cell] - grid.matrix_ice_water_m3[cell];
    const macropore_capacity =
        grid.macropore_pore_capacity_m3[cell] -
        grid.macropore_ice_water_m3[cell];
    const matrix_water = current[matrix_index];
    const macropore_water = current[macropore_index];
    const total_residual =
        residual[matrix_index] + residual[macropore_index];
    const partition_residual = residual[macropore_index];
    const jacobian = jacobian_workspace[0..4];
    const delta = delta_workspace[0..2];

    const total_positive_room = @max(0.0, matrix_capacity - matrix_water);
    const total_negative_room = @max(0.0, matrix_water);
    const total_nominal = std.math.sqrt(std.math.floatEps(f64)) *
        @max(1.0, matrix_water + macropore_water);
    const total_step = @min(
        total_nominal,
        @min(0.5 * total_positive_room, 0.5 * total_negative_room),
    );
    if (total_step <= 1.0e-20) return false;
    @memcpy(probe, current);
    probe[matrix_index] += total_step;
    residualAt(
        grid,
        faces,
        properties,
        base,
        probe,
        target,
        probe_residual,
        scratch,
        trial_micro_flux,
        trial_macro_flux,
    ) catch return false;
    @memcpy(candidate, current);
    candidate[matrix_index] -= total_step;
    residualAt(
        grid,
        faces,
        properties,
        base,
        candidate,
        target,
        candidate_residual,
        scratch,
        trial_micro_flux,
        trial_macro_flux,
    ) catch return false;
    jacobian[0] =
        (probe_residual[matrix_index] + probe_residual[macropore_index] -
            candidate_residual[matrix_index] -
            candidate_residual[macropore_index]) / (2.0 * total_step);
    jacobian[2] =
        (probe_residual[macropore_index] -
            candidate_residual[macropore_index]) / (2.0 * total_step);

    const partition_positive_room = @min(
        @max(0.0, matrix_water),
        @max(0.0, macropore_capacity - macropore_water),
    );
    const partition_negative_room = @min(
        @max(0.0, macropore_water),
        @max(0.0, matrix_capacity - matrix_water),
    );
    const partition_nominal = std.math.sqrt(std.math.floatEps(f64)) *
        @max(1.0, macropore_water);
    const partition_step = @min(
        partition_nominal,
        @min(0.5 * partition_positive_room, 0.5 * partition_negative_room),
    );
    if (partition_step <= 1.0e-20) return false;
    @memcpy(probe, current);
    probe[matrix_index] -= partition_step;
    probe[macropore_index] += partition_step;
    residualAt(
        grid,
        faces,
        properties,
        base,
        probe,
        target,
        probe_residual,
        scratch,
        trial_micro_flux,
        trial_macro_flux,
    ) catch return false;
    @memcpy(candidate, current);
    candidate[matrix_index] += partition_step;
    candidate[macropore_index] -= partition_step;
    residualAt(
        grid,
        faces,
        properties,
        base,
        candidate,
        target,
        candidate_residual,
        scratch,
        trial_micro_flux,
        trial_macro_flux,
    ) catch return false;
    jacobian[1] =
        (probe_residual[matrix_index] + probe_residual[macropore_index] -
            candidate_residual[matrix_index] -
            candidate_residual[macropore_index]) / (2.0 * partition_step);
    jacobian[3] =
        (probe_residual[macropore_index] -
            candidate_residual[macropore_index]) / (2.0 * partition_step);
    delta[0] = -total_residual;
    delta[1] = -partition_residual;
    if (!numerics.solveDenseLinearSystem(jacobian, delta, 2)) return false;

    var best_norm = try scaledNorm(current, residual, options);
    var found = false;
    var fraction: f64 = 1;
    var line: u8 = 0;
    while (line < 24) : (line += 1) {
        @memcpy(candidate, current);
        candidate[matrix_index] =
            matrix_water + fraction * (delta[0] - delta[1]);
        candidate[macropore_index] =
            macropore_water + fraction * delta[1];
        if (candidate[matrix_index] >= 0 and
            candidate[matrix_index] <= matrix_capacity and
            candidate[macropore_index] >= 0 and
            candidate[macropore_index] <= macropore_capacity)
        {
            if (residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            )) |_| {
                const norm =
                    try scaledNorm(candidate, candidate_residual, options);
                if (norm < best_norm) {
                    best_norm = norm;
                    @memcpy(accepted_state, candidate);
                    found = true;
                }
            } else |_| {}
        }
        fraction *= 0.5;
    }
    return found;
}

/// Newton correction for one complete pore domain. It avoids contaminating a
/// well-conditioned matrix Richards column with inactive, bound-constrained
/// macropore coordinates while evaluating acceptance against the full
/// dual-domain nonlinear residual.
fn domainBlockNewton(
    grid: *const grid_module.GridState,
    faces: []const Face,
    properties: Properties,
    options: Options,
    base: []const f64,
    current: []const f64,
    target: []f64,
    residual: []const f64,
    scratch: []f64,
    trial_micro_flux: []f64,
    trial_macro_flux: []f64,
    domain_offset: usize,
    fine_active_set_probe: bool,
    accepted_state: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    probe: []f64,
    probe_residual: []f64,
    jacobian_workspace: []f64,
    delta_workspace: []f64,
) !bool {
    const cells = grid.layer_count;
    if (domain_offset != 0 and domain_offset != cells)
        return error.InvalidSoilWaterDomainOffset;
    const jacobian = jacobian_workspace[0 .. cells * cells];
    const delta = delta_workspace[0..cells];
    for (0..cells) |column_cell| {
        const column = domain_offset + column_cell;
        const capacity =
            if (domain_offset == 0)
                grid.matrix_pore_capacity_m3[column_cell] -
                    grid.matrix_ice_water_m3[column_cell]
            else
                grid.macropore_pore_capacity_m3[column_cell] -
                    grid.macropore_ice_water_m3[column_cell];
        const upward_room = @max(0.0, capacity - current[column]);
        const downward_room = @max(0.0, current[column]);
        const nominal = (if (fine_active_set_probe)
            std.math.sqrt(std.math.floatEps(f64))
        else
            std.math.cbrt(std.math.floatEps(f64))) *
            @max(1.0e-6, @abs(current[column]));
        const central_step =
            @min(nominal, @min(0.5 * upward_room, 0.5 * downward_room));
        if (central_step > 1.0e-20) {
            @memcpy(probe, current);
            probe[column] += central_step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                probe,
                target,
                probe_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            @memcpy(candidate, current);
            candidate[column] -= central_step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            for (0..cells) |row_cell|
                jacobian[row_cell * cells + column_cell] =
                    (probe_residual[domain_offset + row_cell] -
                        candidate_residual[domain_offset + row_cell]) /
                    (2.0 * central_step);
        } else {
            const step =
                if (upward_room >= downward_room)
                    @min(nominal, 0.5 * upward_room)
                else
                    -@min(nominal, 0.5 * downward_room);
            if (@abs(step) <= 1.0e-20) return false;
            @memcpy(probe, current);
            probe[column] += step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                probe,
                target,
                probe_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            for (0..cells) |row_cell|
                jacobian[row_cell * cells + column_cell] =
                    (probe_residual[domain_offset + row_cell] -
                        residual[domain_offset + row_cell]) / step;
        }
    }
    for (0..cells) |row_cell|
        delta[row_cell] = -residual[domain_offset + row_cell];
    if (!numerics.solveDenseLinearSystem(jacobian, delta, cells)) return false;

    var best_norm = try scaledNorm(current, residual, options);
    var found = false;
    var fraction: f64 = 1;
    var line: u8 = 0;
    while (line < 20) : (line += 1) {
        @memcpy(candidate, current);
        var valid = true;
        for (0..cells) |cell| {
            const index = domain_offset + cell;
            const capacity =
                if (domain_offset == 0)
                    grid.matrix_pore_capacity_m3[cell] -
                        grid.matrix_ice_water_m3[cell]
                else
                    grid.macropore_pore_capacity_m3[cell] -
                        grid.macropore_ice_water_m3[cell];
            candidate[index] = current[index] + fraction * delta[cell];
            if (!std.math.isFinite(candidate[index]) or
                candidate[index] < 0 or candidate[index] > capacity)
                valid = false;
        }
        if (valid) {
            if (residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            )) |_| {
                const norm =
                    try scaledNorm(candidate, candidate_residual, options);
                if (norm < best_norm) {
                    best_norm = norm;
                    @memcpy(accepted_state, candidate);
                    found = true;
                }
            } else |_| {}
        }
        fraction *= 0.5;
    }
    return found;
}

/// Resolves the Richards coupling across every face incident on one storage
/// cell as a single bounded Newton block. This is the spatial counterpart of
/// `dualDomainCellNewton`: scalar coordinates cannot reliably cross a
/// donor/receiver active-set change because changing the centre storage also
/// changes both neighboring residual equations.
fn spatialDomainCellNewton(
    grid: *const grid_module.GridState,
    faces: []const Face,
    properties: Properties,
    base: []const f64,
    current: []f64,
    target: []f64,
    residual: []f64,
    scratch: []f64,
    trial_micro_flux: []f64,
    trial_macro_flux: []f64,
    component: usize,
    candidate: []f64,
    candidate_residual: []f64,
    probe: []f64,
    probe_residual: []f64,
    indices_workspace: []usize,
    jacobian_workspace: []f64,
    right_hand_side_workspace: []f64,
    options: Options,
    include_dual_domain: bool,
    use_scaled_merit: bool,
) !bool {
    const cells = grid.layer_count;
    const domain_offset = if (component < cells) @as(usize, 0) else cells;
    const cell = component - domain_offset;
    var local_count: usize = 1;
    indices_workspace[0] = component;
    for (faces) |face| {
        const neighbor =
            if (face.source_cell == cell)
                face.destination_cell
            else if (face.destination_cell == cell)
                face.source_cell
            else
                continue;
        const neighbor_component = domain_offset + neighbor;
        var duplicate = false;
        for (indices_workspace[0..local_count]) |index|
            if (index == neighbor_component) {
                duplicate = true;
                break;
            };
        if (!duplicate) {
            indices_workspace[local_count] = neighbor_component;
            local_count += 1;
        }
    }
    // Close one additional face ring so changing an immediate neighbor does
    // not leave its outward Richards face outside the Newton block.
    const first_ring_count = local_count;
    for (faces) |face| {
        const source_component = domain_offset + face.source_cell;
        const destination_component = domain_offset + face.destination_cell;
        var source_in_first_ring = false;
        var destination_in_first_ring = false;
        for (indices_workspace[0..first_ring_count]) |index| {
            source_in_first_ring = source_in_first_ring or index == source_component;
            destination_in_first_ring =
                destination_in_first_ring or index == destination_component;
        }
        const outward_component =
            if (source_in_first_ring and !destination_in_first_ring)
                destination_component
            else if (destination_in_first_ring and !source_in_first_ring)
                source_component
            else
                continue;
        var duplicate = false;
        for (indices_workspace[0..local_count]) |index|
            if (index == outward_component) {
                duplicate = true;
                break;
            };
        if (!duplicate) {
            indices_workspace[local_count] = outward_component;
            local_count += 1;
        }
    }
    if (include_dual_domain and
        properties.dual_domain_exchange_enabled.len != 0)
    {
        const spatial_count = local_count;
        for (indices_workspace[0..spatial_count]) |spatial_index| {
            const spatial_cell = spatial_index - domain_offset;
            if (!properties.dual_domain_exchange_enabled[spatial_cell] or
                grid.macropore_pore_capacity_m3[spatial_cell] <= 0)
                continue;
            indices_workspace[local_count] =
                if (domain_offset == 0)
                    cells + spatial_cell
                else
                    spatial_cell;
            local_count += 1;
        }
    }
    if (local_count < 2) return false;

    const indices = indices_workspace[0..local_count];
    const local_jacobian = jacobian_workspace[0 .. local_count * local_count];
    const local_delta = right_hand_side_workspace[0..local_count];
    for (indices, 0..) |column_index, column| {
        const column_is_macro = column_index >= cells;
        const column_cell =
            if (column_is_macro) column_index - cells else column_index;
        const capacity =
            if (!column_is_macro)
                grid.matrix_pore_capacity_m3[column_cell] -
                    grid.matrix_ice_water_m3[column_cell]
            else
                grid.macropore_pore_capacity_m3[column_cell] -
                    grid.macropore_ice_water_m3[column_cell];
        const upward_room = @max(0.0, capacity - current[column_index]);
        const downward_room = @max(0.0, current[column_index]);
        const nominal = @max(
            1.0e-6,
            std.math.sqrt(std.math.floatEps(f64)) *
                @abs(current[column_index]),
        );
        const central_step = @min(
            nominal,
            @min(0.5 * upward_room, 0.5 * downward_room),
        );
        if (central_step > 1.0e-20) {
            @memcpy(probe, current);
            probe[column_index] += central_step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                probe,
                target,
                probe_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            @memcpy(candidate, current);
            candidate[column_index] -= central_step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            for (indices, 0..) |row_index, row|
                local_jacobian[row * local_count + column] =
                    (probe_residual[row_index] -
                        candidate_residual[row_index]) / (2.0 * central_step);
            continue;
        }
        const step =
            if (upward_room >= downward_room)
                @min(nominal, 0.5 * upward_room)
            else
                -@min(nominal, 0.5 * downward_room);
        if (@abs(step) <= 1.0e-20) return false;
        @memcpy(probe, current);
        probe[column_index] += step;
        residualAt(
            grid,
            faces,
            properties,
            base,
            probe,
            target,
            probe_residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
        ) catch return false;
        for (indices, 0..) |row_index, row|
            local_jacobian[row * local_count + column] =
                (probe_residual[row_index] - residual[row_index]) / step;
    }
    for (indices, 0..) |index, row|
        local_delta[row] = -residual[index];
    if (!numerics.solveDenseLinearSystem(local_jacobian, local_delta, local_count))
        return false;

    var best_global_norm: f64 = if (use_scaled_merit)
        try scaledNorm(current, residual, options)
    else
        0;
    if (!use_scaled_merit) {
        for (residual) |value|
            best_global_norm = @max(best_global_norm, @abs(value));
    }
    var found = false;
    var fraction: f64 = 1;
    var line: u8 = 0;
    while (line < 20) : (line += 1) {
        @memcpy(candidate, current);
        var valid = true;
        for (indices, local_delta) |index, delta| {
            const index_is_macro = index >= cells;
            const index_cell = if (index_is_macro) index - cells else index;
            const capacity =
                if (!index_is_macro)
                    grid.matrix_pore_capacity_m3[index_cell] -
                        grid.matrix_ice_water_m3[index_cell]
                else
                    grid.macropore_pore_capacity_m3[index_cell] -
                        grid.macropore_ice_water_m3[index_cell];
            candidate[index] = current[index] + fraction * delta;
            if (!std.math.isFinite(candidate[index]) or
                candidate[index] < 0 or candidate[index] > capacity)
                valid = false;
        }
        if (valid) {
            if (residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            )) |_| {
                var global_norm: f64 = if (use_scaled_merit)
                    try scaledNorm(candidate, candidate_residual, options)
                else
                    0;
                if (!use_scaled_merit) {
                    for (candidate_residual) |value|
                        global_norm = @max(global_norm, @abs(value));
                }
                if (global_norm < best_global_norm) {
                    best_global_norm = global_norm;
                    @memcpy(probe, candidate);
                    @memcpy(probe_residual, candidate_residual);
                    found = true;
                }
            } else |_| {}
        }
        fraction *= 0.5;
    }
    if (!found) return false;
    @memcpy(current, probe);
    @memcpy(residual, probe_residual);
    return true;
}

fn dualDomainCellNewton(
    grid: *const grid_module.GridState,
    faces: []const Face,
    properties: Properties,
    base: []const f64,
    current: []f64,
    target: []f64,
    residual: []f64,
    scratch: []f64,
    trial_micro_flux: []f64,
    trial_macro_flux: []f64,
    cell: usize,
    candidate: []f64,
    candidate_residual: []f64,
    probe: []f64,
    probe_residual: []f64,
) !bool {
    const cells = grid.layer_count;
    const indices = [2]usize{ cell, cells + cell };
    var local_jacobian: [4]f64 = undefined;
    for (indices, 0..) |column_index, column| {
        const capacity =
            if (column == 0)
                grid.matrix_pore_capacity_m3[cell] -
                    grid.matrix_ice_water_m3[cell]
            else
                grid.macropore_pore_capacity_m3[cell] -
                    grid.macropore_ice_water_m3[cell];
        const upward_room = @max(0.0, capacity - current[column_index]);
        const downward_room = @max(0.0, current[column_index]);
        const nominal = std.math.cbrt(std.math.floatEps(f64)) *
            @max(1.0e-6, @abs(current[column_index]));
        const central_step = @min(
            nominal,
            @min(0.5 * upward_room, 0.5 * downward_room),
        );
        if (central_step > 1.0e-20) {
            @memcpy(probe, current);
            probe[column_index] += central_step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                probe,
                target,
                probe_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            @memcpy(candidate, current);
            candidate[column_index] -= central_step;
            residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            ) catch return false;
            for (indices, 0..) |row_index, row|
                local_jacobian[row * 2 + column] =
                    (probe_residual[row_index] -
                        candidate_residual[row_index]) / (2 * central_step);
            continue;
        }
        const step =
            if (upward_room >= downward_room)
                @min(nominal, 0.5 * upward_room)
            else
                -@min(nominal, 0.5 * downward_room);
        if (@abs(step) <= 1.0e-20) return false;
        @memcpy(probe, current);
        probe[column_index] += step;
        residualAt(
            grid,
            faces,
            properties,
            base,
            probe,
            target,
            probe_residual,
            scratch,
            trial_micro_flux,
            trial_macro_flux,
        ) catch return false;
        for (indices, 0..) |row_index, row|
            local_jacobian[row * 2 + column] =
                (probe_residual[row_index] - residual[row_index]) / step;
    }
    var local_delta = [2]f64{
        -residual[indices[0]],
        -residual[indices[1]],
    };
    if (!numerics.solveDenseLinearSystem(&local_jacobian, &local_delta, 2))
        return false;
    const current_pair_norm =
        @max(@abs(residual[indices[0]]), @abs(residual[indices[1]]));
    var best_pair_norm = current_pair_norm;
    var found = false;
    var fraction: f64 = 1;
    var line: u8 = 0;
    while (line < 20) : (line += 1) {
        @memcpy(candidate, current);
        var valid = true;
        for (indices, local_delta, 0..) |index, delta, domain| {
            const capacity =
                if (domain == 0)
                    grid.matrix_pore_capacity_m3[cell] -
                        grid.matrix_ice_water_m3[cell]
                else
                    grid.macropore_pore_capacity_m3[cell] -
                        grid.macropore_ice_water_m3[cell];
            candidate[index] = current[index] + fraction * delta;
            if (!std.math.isFinite(candidate[index]) or
                candidate[index] < 0 or candidate[index] > capacity)
                valid = false;
        }
        if (valid) {
            if (residualAt(
                grid,
                faces,
                properties,
                base,
                candidate,
                target,
                candidate_residual,
                scratch,
                trial_micro_flux,
                trial_macro_flux,
            )) |_| {
                const pair_norm = @max(
                    @abs(candidate_residual[indices[0]]),
                    @abs(candidate_residual[indices[1]]),
                );
                if (pair_norm < best_pair_norm) {
                    best_pair_norm = pair_norm;
                    @memcpy(probe, candidate);
                    @memcpy(probe_residual, candidate_residual);
                    found = true;
                }
            } else |_| {}
        }
        fraction *= 0.5;
    }
    if (!found) return false;
    @memcpy(current, probe);
    @memcpy(residual, probe_residual);
    return true;
}

/// Solves the shared runtime topology and publishes converged FLWM/FLWHM
/// directly to the water and solute transport views.
pub fn solveAndBindTransportFaces(allocator: std.mem.Allocator, grid: *grid_module.GridState, hydrology: *transport_hydrology.State, shared_faces: *transport_hydrology.SoilFaces, geometry: FaceGeometry, properties: Properties, options: Options) !Result {
    const count = shared_faces.micropore_faces.len;
    if (shared_faces.macropore_faces.len != count or shared_faces.direction_axis.len != count or geometry.source_path_length_m.len != count or geometry.destination_path_length_m.len != count or geometry.face_area_m2.len != count or geometry.macropore_hydraulic_conductance_m_per_h_megapascal.len != count) return error.SoilWaterFaceGeometryDimensionMismatch;
    const faces = try allocator.alloc(Face, count);
    defer allocator.free(faces);
    for (faces, 0..) |*face, index| {
        const axis: Axis = @enumFromInt(shared_faces.direction_axis[index]);
        face.* = .{ .source_cell = shared_faces.micropore_faces[index].first_cell, .destination_cell = shared_faces.micropore_faces[index].second_cell, .axis = axis, .direction = if (axis == .z) .vertical else .horizontal, .source_path_length_m = geometry.source_path_length_m[index], .destination_path_length_m = geometry.destination_path_length_m[index], .face_area_m2 = geometry.face_area_m2[index], .macropore_hydraulic_conductance_m_per_h_megapascal = geometry.macropore_hydraulic_conductance_m_per_h_megapascal[index] };
    }
    const artificial_drainage = try allocator.alloc(f64, grid.cell_count);
    defer allocator.free(artificial_drainage);
    @memset(artificial_drainage, 0);
    var active_properties = properties;
    active_properties.artificial_drainage_outflow_m3_per_step = artificial_drainage;
    const result = try solve(allocator, grid, faces, active_properties, shared_faces.micropore_water_flux_m3_per_step, shared_faces.macropore_water_flux_m3_per_step, options);
    @memcpy(hydrology.micropore_water_volume_m3, grid.matrix_liquid_water_m3);
    @memcpy(hydrology.macropore_water_volume_m3, grid.macropore_liquid_water_m3);
    @memcpy(hydrology.matrix_air_volume_m3, grid.matrix_air_volume_m3);
    @memcpy(hydrology.macropore_air_volume_m3, grid.macropore_air_volume_m3);
    @memcpy(hydrology.air_volume_m3, grid.air_volume_m3);
    @memcpy(hydrology.artificial_drainage_outflow_m3_per_step, artificial_drainage);
    @memset(hydrology.micropore_face_flux_m3_per_step, 0);
    @memset(hydrology.macropore_face_flux_m3_per_step, 0);
    for (shared_faces.micropore_faces, shared_faces.macropore_faces, shared_faces.direction_axis, shared_faces.micropore_water_flux_m3_per_step, shared_faces.macropore_water_flux_m3_per_step) |*micro_face, *macro_face, axis, micro, macro| {
        micro_face.water_flux_m3_per_step = micro;
        macro_face.water_flux_m3_per_step = macro;
        hydrology.micropore_face_flux_m3_per_step[micro_face.first_cell * 3 + axis] = micro;
        hydrology.macropore_face_flux_m3_per_step[macro_face.first_cell * 3 + axis] = macro;
    }
    try hydrology.validateFinite();
    return result;
}

fn componentRoot(parent: []usize, cell: usize) usize {
    var root = cell;
    while (parent[root] != root) root = parent[root];
    var current = cell;
    while (parent[current] != current) {
        const next = parent[current];
        parent[current] = root;
        current = next;
    }
    return root;
}

fn unionComponents(parent: []usize, component_size: []usize, first: usize, second: usize) void {
    var first_root = componentRoot(parent, first);
    var second_root = componentRoot(parent, second);
    if (first_root == second_root) return;
    if (component_size[first_root] < component_size[second_root]) std.mem.swap(usize, &first_root, &second_root);
    parent[second_root] = first_root;
    component_size[first_root] += component_size[second_root];
}

fn releaseIndependentStorageCoordinates(parent: []usize, component_size: []usize) void {
    // The whole-step residual is target - trial, so every storage contributes
    // a -1 identity derivative. Unlike a flux-only Laplacian, this system has
    // no mass nullspace: solving every storage coordinate is nonsingular and
    // the conservative target equations themselves enforce water balance.
    for (0..parent.len) |cell| component_size[componentRoot(parent, cell)] = 1;
}

fn solveConservedNewtonSystem(full_jacobian: []const f64, residual: []const f64, expanded_delta: []f64, reduced_jacobian: []f64, reduced_right_hand_side: []f64, reduced_index_by_component: []usize, dimension: usize, cells: usize, micropore_parent: []usize, micropore_component_size: []const usize, macropore_parent: []usize, macropore_component_size: []const usize) bool {
    const excluded = std.math.maxInt(usize);
    @memset(reduced_index_by_component, excluded);
    var reduced_dimension: usize = 0;
    for (0..dimension) |full_index| {
        const domain_offset = if (full_index < cells) @as(usize, 0) else cells;
        const cell = full_index - domain_offset;
        const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
        const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
        const root = componentRoot(parent, cell);
        if (component_size[root] > 1 and cell == root) continue;
        reduced_index_by_component[full_index] = reduced_dimension;
        reduced_dimension += 1;
    }
    if (reduced_dimension == 0) {
        @memset(expanded_delta, 0);
        return true;
    }
    for (0..dimension) |full_row| {
        const reduced_row = reduced_index_by_component[full_row];
        if (reduced_row == excluded) continue;
        reduced_right_hand_side[reduced_row] = -residual[full_row];
        for (0..dimension) |full_column| {
            const reduced_column = reduced_index_by_component[full_column];
            if (reduced_column == excluded) continue;
            const domain_offset = if (full_column < cells) @as(usize, 0) else cells;
            const cell = full_column - domain_offset;
            const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
            const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
            const root = componentRoot(parent, cell);
            const anchor_column = domain_offset + root;
            const anchor_coefficient = if (component_size[root] > 1) full_jacobian[full_row * dimension + anchor_column] else 0;
            reduced_jacobian[reduced_row * reduced_dimension + reduced_column] = full_jacobian[full_row * dimension + full_column] - anchor_coefficient;
        }
    }
    if (!numerics.solveDenseLinearSystem(reduced_jacobian[0 .. reduced_dimension * reduced_dimension], reduced_right_hand_side[0..reduced_dimension], reduced_dimension)) return false;
    @memset(expanded_delta, 0);
    for (0..dimension) |full_index| {
        const reduced_index = reduced_index_by_component[full_index];
        if (reduced_index != excluded) expanded_delta[full_index] = reduced_right_hand_side[reduced_index];
    }
    for (0..2) |domain_index| {
        const domain_offset = if (domain_index == 0) @as(usize, 0) else cells;
        const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
        const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
        for (0..cells) |cell| {
            const root = componentRoot(parent, cell);
            if (cell != root or component_size[root] <= 1) continue;
            var sum: f64 = 0;
            for (0..cells) |member| {
                if (componentRoot(parent, member) == root and member != root) sum += expanded_delta[domain_offset + member];
            }
            expanded_delta[domain_offset + root] = -sum;
        }
    }
    return true;
}

fn solveConservedTrustRegionSystem(full_jacobian: []const f64, residual: []const f64, state: []const f64, options: Options, expanded_delta: []f64, normal_matrix: []f64, reduced_right_hand_side: []f64, reduced_index_by_component: []usize, dimension: usize, cells: usize, micropore_parent: []usize, micropore_component_size: []const usize, macropore_parent: []usize, macropore_component_size: []const usize, damping: f64) bool {
    if (!std.math.isFinite(damping) or damping <= 0) return false;
    const excluded = std.math.maxInt(usize);
    @memset(reduced_index_by_component, excluded);
    var reduced_dimension: usize = 0;
    for (0..dimension) |full_index| {
        const domain_offset = if (full_index < cells) @as(usize, 0) else cells;
        const cell = full_index - domain_offset;
        const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
        const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
        const root = componentRoot(parent, cell);
        if (component_size[root] > 1 and cell == root) continue;
        reduced_index_by_component[full_index] = reduced_dimension;
        reduced_dimension += 1;
    }
    if (reduced_dimension == 0) {
        @memset(expanded_delta, 0);
        return true;
    }
    @memset(normal_matrix[0 .. reduced_dimension * reduced_dimension], 0);
    @memset(reduced_right_hand_side[0..reduced_dimension], 0);
    for (0..dimension) |full_row| {
        const residual_scale = options.absolute_tolerance_m3 +
            options.relative_tolerance *
                @max(1.0, @abs(state[full_row]));
        const scaled_residual = residual[full_row] / residual_scale;
        for (0..dimension) |first_full_column| {
            const first_reduced_column = reduced_index_by_component[first_full_column];
            if (first_reduced_column == excluded) continue;
            const first_variable_scale =
                @max(1.0, @abs(state[first_full_column]));
            const first_coefficient = conservedCoordinateCoefficient(full_jacobian, full_row, first_full_column, dimension, cells, micropore_parent, micropore_component_size, macropore_parent, macropore_component_size) *
                first_variable_scale / residual_scale;
            reduced_right_hand_side[first_reduced_column] -= first_coefficient * scaled_residual;
            for (0..dimension) |second_full_column| {
                const second_reduced_column = reduced_index_by_component[second_full_column];
                if (second_reduced_column == excluded) continue;
                const second_variable_scale =
                    @max(1.0, @abs(state[second_full_column]));
                const second_coefficient = conservedCoordinateCoefficient(full_jacobian, full_row, second_full_column, dimension, cells, micropore_parent, micropore_component_size, macropore_parent, macropore_component_size) *
                    second_variable_scale / residual_scale;
                normal_matrix[first_reduced_column * reduced_dimension + second_reduced_column] += first_coefficient * second_coefficient;
            }
        }
    }
    for (0..reduced_dimension) |index| {
        const diagonal_index = index * reduced_dimension + index;
        normal_matrix[diagonal_index] += damping * @max(1.0e-18, normal_matrix[diagonal_index]);
    }
    if (!numerics.solveDenseLinearSystem(normal_matrix[0 .. reduced_dimension * reduced_dimension], reduced_right_hand_side[0..reduced_dimension], reduced_dimension)) return false;
    expandConservedDelta(expanded_delta, reduced_right_hand_side, reduced_index_by_component, dimension, cells, micropore_parent, micropore_component_size, macropore_parent, macropore_component_size);
    for (expanded_delta, state) |*value, state_value|
        value.* *= @max(1.0, @abs(state_value));
    return true;
}

fn solveProjectedTrustRegionSystem(jacobian: []const f64, residual: []const f64, delta: []f64, normal_matrix: []f64, right_hand_side: []f64, dimension: usize, cells: usize, micropore_parent: []usize, micropore_component_size: []const usize, macropore_parent: []usize, macropore_component_size: []const usize, damping: f64) bool {
    if (!std.math.isFinite(damping) or damping <= 0) return false;
    @memset(normal_matrix[0 .. dimension * dimension], 0);
    @memset(right_hand_side[0..dimension], 0);
    for (0..dimension) |row| for (0..dimension) |first_column| {
        const first = jacobian[row * dimension + first_column];
        right_hand_side[first_column] -= first * residual[row];
        for (0..dimension) |second_column| normal_matrix[first_column * dimension + second_column] += first * jacobian[row * dimension + second_column];
    };
    for (0..dimension) |index| {
        const diagonal = index * dimension + index;
        normal_matrix[diagonal] += damping * @max(1.0e-18, normal_matrix[diagonal]);
    }
    if (!numerics.solveDenseLinearSystem(normal_matrix[0 .. dimension * dimension], right_hand_side[0..dimension], dimension)) return false;
    @memcpy(delta, right_hand_side[0..dimension]);
    // Orthogonally project onto each numerically connected domain's
    // conservative subspace rather than eliminating an anchor coordinate.
    for (0..2) |domain_index| {
        const offset = if (domain_index == 0) @as(usize, 0) else cells;
        const parent = if (domain_index == 0) micropore_parent else macropore_parent;
        const component_size = if (domain_index == 0) micropore_component_size else macropore_component_size;
        for (0..cells) |root_candidate| {
            const root = componentRoot(parent, root_candidate);
            if (root != root_candidate or component_size[root] <= 1) continue;
            var sum: f64 = 0;
            for (0..cells) |member| {
                if (componentRoot(parent, member) == root) sum += delta[offset + member];
            }
            const mean = sum / @as(f64, @floatFromInt(component_size[root]));
            for (0..cells) |member| {
                if (componentRoot(parent, member) == root) delta[offset + member] -= mean;
            }
        }
    }
    for (delta) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn conservedCoordinateCoefficient(full_jacobian: []const f64, full_row: usize, full_column: usize, dimension: usize, cells: usize, micropore_parent: []usize, micropore_component_size: []const usize, macropore_parent: []usize, macropore_component_size: []const usize) f64 {
    const domain_offset = if (full_column < cells) @as(usize, 0) else cells;
    const cell = full_column - domain_offset;
    const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
    const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
    const root = componentRoot(parent, cell);
    const anchor_coefficient = if (component_size[root] > 1) full_jacobian[full_row * dimension + domain_offset + root] else 0;
    return full_jacobian[full_row * dimension + full_column] - anchor_coefficient;
}

fn expandConservedDelta(expanded_delta: []f64, reduced_delta: []const f64, reduced_index_by_component: []const usize, dimension: usize, cells: usize, micropore_parent: []usize, micropore_component_size: []const usize, macropore_parent: []usize, macropore_component_size: []const usize) void {
    const excluded = std.math.maxInt(usize);
    @memset(expanded_delta, 0);
    for (0..dimension) |full_index| {
        const reduced_index = reduced_index_by_component[full_index];
        if (reduced_index != excluded) expanded_delta[full_index] = reduced_delta[reduced_index];
    }
    for (0..2) |domain_index| {
        const domain_offset = if (domain_index == 0) @as(usize, 0) else cells;
        const parent = if (domain_offset == 0) micropore_parent else macropore_parent;
        const component_size = if (domain_offset == 0) micropore_component_size else macropore_component_size;
        for (0..cells) |cell| {
            const root = componentRoot(parent, cell);
            if (cell != root or component_size[root] <= 1) continue;
            var sum: f64 = 0;
            for (0..cells) |member| {
                if (componentRoot(parent, member) == root and member != root) sum += expanded_delta[domain_offset + member];
            }
            expanded_delta[domain_offset + root] = -sum;
        }
    }
}

fn rememberIteration(current: []const f64, residual: []const f64, previous_state: []f64, previous_residual: []f64, previous_previous_state: []f64, previous_previous_residual: []f64, history_count: *u8) void {
    if (history_count.* > 0) {
        @memcpy(previous_previous_state, previous_state);
        @memcpy(previous_previous_residual, previous_residual);
    }
    @memcpy(previous_state, current);
    @memcpy(previous_residual, residual);
    history_count.* = @min(@as(u8, 2), history_count.* + 1);
}

fn residualAt(grid: *const grid_module.GridState, faces: []const Face, properties: Properties, base: []const f64, trial: []const f64, target: []f64, residual: []f64, scratch: []f64, micro_fluxes: []f64, macro_fluxes: []f64) !void {
    if (properties.artificial_drainage_outflow_m3_per_step) |outflow| {
        if (outflow.len != grid.cell_count) return error.ArtificialDrainageDimensionMismatch;
        @memset(outflow, 0);
    }
    const cells = grid.layer_count;
    @memcpy(scratch, trial);
    @memcpy(target, base);
    if (properties.matrix_external_source_m3_per_step.len != 0)
        for (
            target[0..cells],
            properties.matrix_external_source_m3_per_step,
        ) |*water_m3, source_m3| {
            if (!std.math.isFinite(source_m3) or source_m3 < 0)
                return error.InvalidSoilWaterExternalSource;
            water_m3.* += source_m3;
            if (!std.math.isFinite(water_m3.*))
                return error.NonFiniteSoilWaterExternalSource;
        };
    @memset(micro_fluxes, 0);
    @memset(macro_fluxes, 0);
    for (0..cells) |cell| {
        if (trial[cell] + grid.matrix_ice_water_m3[cell] > grid.matrix_pore_capacity_m3[cell] + poreCapacityRoundoffToleranceM3(grid.matrix_pore_capacity_m3[cell]) or trial[cells + cell] + grid.macropore_ice_water_m3[cell] > grid.macropore_pore_capacity_m3[cell] + poreCapacityRoundoffToleranceM3(grid.macropore_pore_capacity_m3[cell])) return error.SoilWaterCandidateExceedsPoreCapacity;
    }
    for (faces, 0..) |face, face_index| {
        const source = face.source_cell;
        const destination = face.destination_cell;
        const direction_index: usize = @intFromEnum(face.axis);
        const source_water = scratch[source];
        const destination_water = scratch[destination];
        const source_air = @max(0.0, grid.matrix_pore_capacity_m3[source] - grid.matrix_ice_water_m3[source] - source_water);
        const destination_air = @max(0.0, grid.matrix_pore_capacity_m3[destination] - grid.matrix_ice_water_m3[destination] - destination_water);
        const source_fraction = source_water / properties.matrix_bulk_volume_m3[source];
        const destination_fraction = destination_water / properties.matrix_bulk_volume_m3[destination];
        const source_matric = try matricPotentialMpaAt(properties, source, source_fraction);
        const destination_matric = try matricPotentialMpaAt(properties, destination, destination_fraction);
        const source_total = source_matric + properties.gravitational_potential_megapascal[source] + properties.osmotic_potential_multiplier * properties.osmotic_potential_megapascal[source];
        const destination_total = destination_matric + properties.gravitational_potential_megapascal[destination] + properties.osmotic_potential_multiplier * properties.osmotic_potential_megapascal[destination];
        const source_conductivity = try conductivityAt(properties, source, direction_index, source_fraction, grid.matrix_ice_water_m3[source]);
        const destination_conductivity = try conductivityAt(properties, destination, direction_index, destination_fraction, grid.matrix_ice_water_m3[destination]);
        const matrix = try water_flux.calculateMatrixFaceFlux(.{ .direction = face.direction, .source_water_m3 = source_water, .destination_water_m3 = destination_water, .source_air_m3 = source_air, .destination_air_m3 = destination_air, .source_micropore_volume_m3 = properties.matrix_bulk_volume_m3[source], .destination_micropore_volume_m3 = properties.matrix_bulk_volume_m3[destination], .source_water_fraction = source_fraction, .destination_water_fraction = destination_fraction, .source_total_water_potential_megapascal = source_total, .destination_total_water_potential_megapascal = destination_total, .source_hydraulic_conductivity_m2_per_h_megapascal = source_conductivity, .destination_hydraulic_conductivity_m2_per_h_megapascal = destination_conductivity, .source_path_length_m = face.source_path_length_m, .destination_path_length_m = face.destination_path_length_m, .face_area_m2 = face.face_area_m2, .time_fraction = properties.nonlinear_time_fraction });
        const matrix_flux_m3 = limitFluxForAssembledTarget(
            matrix.limited_water_m3,
            target[0..cells],
            source,
            destination,
            grid.matrix_pore_capacity_m3[source] -
                grid.matrix_ice_water_m3[source],
            grid.matrix_pore_capacity_m3[destination] -
                grid.matrix_ice_water_m3[destination],
        );
        applyConservativeFlux(
            target[0..cells],
            source,
            destination,
            matrix_flux_m3,
        );
        // FLWM is assigned from FLWL=FLQL in WATSUB; FLQ2 is retained by the
        // face kernel for the separate FLWLX unsaturated-water diagnostic.
        micro_fluxes[face_index] = matrix_flux_m3;

        const source_macro = scratch[cells + source];
        const destination_macro = scratch[cells + destination];
        if (grid.macropore_pore_capacity_m3[source] > 0 and grid.macropore_pore_capacity_m3[destination] > 0) {
            const source_macro_air = @max(0.0, grid.macropore_pore_capacity_m3[source] - grid.macropore_ice_water_m3[source] - source_macro);
            const destination_macro_air = @max(0.0, grid.macropore_pore_capacity_m3[destination] - grid.macropore_ice_water_m3[destination] - destination_macro);
            const calculated_macro_flux_m3 = richards: {
                const source_saturation = source_macro / grid.macropore_pore_capacity_m3[source];
                const destination_saturation = destination_macro / grid.macropore_pore_capacity_m3[destination];
                const source_parameters = properties.macropore_mualem_van_genuchten_parameters[source];
                const destination_parameters = properties.macropore_mualem_van_genuchten_parameters[destination];
                const source_head_m = try source_parameters.pressureHeadAtWaterContent(std.math.clamp(source_saturation, source_parameters.residual_water_content_m3_per_m3, 1));
                const destination_head_m = try destination_parameters.pressureHeadAtWaterContent(std.math.clamp(destination_saturation, destination_parameters.residual_water_content_m3_per_m3, 1));
                const source_macropore_conductivity = try source_parameters.hydraulicConductivityMPerH(source_head_m) / properties.gravitational_water_potential_mpa_per_m *
                    frozenHydraulicImpedance(
                        properties.frozen_hydraulic_impedance_exponent,
                        grid.macropore_ice_water_m3[source] / grid.macropore_pore_capacity_m3[source],
                        source_parameters,
                    );
                const destination_macropore_conductivity = try destination_parameters.hydraulicConductivityMPerH(destination_head_m) / properties.gravitational_water_potential_mpa_per_m *
                    frozenHydraulicImpedance(
                        properties.frozen_hydraulic_impedance_exponent,
                        grid.macropore_ice_water_m3[destination] / grid.macropore_pore_capacity_m3[destination],
                        destination_parameters,
                    );
                const flux = try water_flux.calculateMatrixFaceFlux(.{
                    .direction = face.direction,
                    .source_water_m3 = source_macro,
                    .destination_water_m3 = destination_macro,
                    .source_air_m3 = source_macro_air,
                    .destination_air_m3 = destination_macro_air,
                    .source_micropore_volume_m3 = grid.macropore_pore_capacity_m3[source],
                    .destination_micropore_volume_m3 = grid.macropore_pore_capacity_m3[destination],
                    .source_water_fraction = source_saturation,
                    .destination_water_fraction = destination_saturation,
                    .source_total_water_potential_megapascal = source_head_m * properties.gravitational_water_potential_mpa_per_m + properties.gravitational_potential_megapascal[source],
                    .destination_total_water_potential_megapascal = destination_head_m * properties.gravitational_water_potential_mpa_per_m + properties.gravitational_potential_megapascal[destination],
                    .source_hydraulic_conductivity_m2_per_h_megapascal = source_macropore_conductivity,
                    .destination_hydraulic_conductivity_m2_per_h_megapascal = destination_macropore_conductivity,
                    .source_path_length_m = face.source_path_length_m,
                    .destination_path_length_m = face.destination_path_length_m,
                    .face_area_m2 = face.face_area_m2,
                                        .time_fraction = properties.nonlinear_time_fraction,
                });
                break :richards flux.limited_water_m3;
            };
            const macro_flux_m3 = limitFluxForAssembledTarget(
                calculated_macro_flux_m3,
                target[cells..],
                source,
                destination,
                grid.macropore_pore_capacity_m3[source] -
                    grid.macropore_ice_water_m3[source],
                grid.macropore_pore_capacity_m3[destination] -
                    grid.macropore_ice_water_m3[destination],
            );
            applyConservativeFlux(
                target[cells..],
                source,
                destination,
                macro_flux_m3,
            );
            macro_fluxes[face_index] = macro_flux_m3;
        }
    }
    if (properties.dual_domain_exchange_enabled.len != 0) for (0..cells) |cell| {
        if (!properties.dual_domain_exchange_enabled[cell] or
            grid.macropore_pore_capacity_m3[cell] <= 0)
            continue;
        const matrix_water = scratch[cell];
        const macropore_water = scratch[cells + cell];
        const matrix_air = @max(0.0, grid.matrix_pore_capacity_m3[cell] -
            grid.matrix_ice_water_m3[cell] - matrix_water);
        const macropore_air = @max(0.0, grid.macropore_pore_capacity_m3[cell] -
            grid.macropore_ice_water_m3[cell] - macropore_water);
        // Potentials are evaluated at the nonlinear trial, while the flux is
        // committed to the hour-start target. Bound transfer by both views so
        // neither the trial path nor the authoritative conservative target
        // can overdraw a donor or overfill a receiver away from the root.
        const assembled_matrix_water = @max(0.0, target[cell]);
        const assembled_macropore_water =
            @max(0.0, target[cells + cell]);
        const assembled_matrix_air = @max(
            0.0,
            grid.matrix_pore_capacity_m3[cell] -
                grid.matrix_ice_water_m3[cell] -
                assembled_matrix_water,
        );
        const assembled_macropore_air = @max(
            0.0,
            grid.macropore_pore_capacity_m3[cell] -
                grid.macropore_ice_water_m3[cell] -
                assembled_macropore_water,
        );
        const matrix_parameters = properties.mualem_van_genuchten_parameters[cell];
        var exchange_matrix_parameters = matrix_parameters;
        exchange_matrix_parameters.saturated_hydraulic_conductivity_m_per_h *=
            frozenHydraulicImpedance(
                properties.frozen_hydraulic_impedance_exponent,
                grid.matrix_ice_water_m3[cell] /
                    properties.matrix_bulk_volume_m3[cell],
                matrix_parameters,
            );
        const macropore_parameters =
            properties.macropore_mualem_van_genuchten_parameters[cell];
        const matrix_fraction = matrix_water / properties.matrix_bulk_volume_m3[cell];
        const macropore_saturation =
            macropore_water / grid.macropore_pore_capacity_m3[cell];
        const matrix_pressure_head_m = try matrix_parameters.pressureHeadAtWaterContent(
            std.math.clamp(
                matrix_fraction,
                matrix_parameters.residual_water_content_m3_per_m3,
                matrix_parameters.saturated_water_content_m3_per_m3,
            ),
        );
        const macropore_pressure_head_m =
            try macropore_parameters.pressureHeadAtWaterContent(
                std.math.clamp(
                    macropore_saturation,
                    macropore_parameters.residual_water_content_m3_per_m3,
                    macropore_parameters.saturated_water_content_m3_per_m3,
                ),
            );
        const exchange_m3 = try water_flux.calculateDualDomainExchange(.{
            .matrix_parameters = exchange_matrix_parameters,
            .matrix_pressure_head_m = matrix_pressure_head_m,
            .macropore_pressure_head_m = macropore_pressure_head_m,
            .matrix_bulk_volume_m3 = properties.matrix_bulk_volume_m3[cell],
            .characteristic_matrix_length_m = properties.macropore_spacing_m[cell] -
                properties.macropore_radius_m[cell],
            .geometry_factor = properties.dual_domain_geometry_factor,
            .scaling_coefficient = properties.dual_domain_scaling_coefficient,
            .time_fraction = properties.nonlinear_time_fraction,
            .matrix_water_m3 = @min(matrix_water, assembled_matrix_water),
            .matrix_air_m3 = @min(matrix_air, assembled_matrix_air),
            .macropore_water_m3 = @min(macropore_water, assembled_macropore_water),
            .macropore_air_m3 = @min(macropore_air, assembled_macropore_air),
        });
        target[cell] += exchange_m3;
        target[cells + cell] -= exchange_m3;
    };
    // Lower boundaries participate in the same nonlinear residual as internal
    // faces. WATSUB stores an oriented face flux; multiplying by XN converts
    // it to the source-layer storage change (drainage is negative).
    if (properties.boundary_topology) |topology| for (topology.faces) |boundary_face| {
        const layer = boundary_face.layer_index;
        const matrix_water = scratch[layer];
        const macropore_water = scratch[cells + layer];
        const matrix_fraction = matrix_water / properties.matrix_bulk_volume_m3[layer];
        var matrix_change: f64 = 0;
        var macropore_change: f64 = 0;
        if (boundary_face.is_lower_boundary) {
            if (boundary_face.natural_exchange_fraction == 0) continue;
            const oriented = try water_boundary.freeDrainage(.{ .direction_sign = boundary_face.direction_sign, .slope_sine = 1, .matrix_hydraulic_conductivity_m2_per_h_megapascal = try conductivityAt(properties, layer, 2, matrix_fraction, grid.matrix_ice_water_m3[layer]), .macropore_hydraulic_conductivity_m2_per_h_megapascal = properties.boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal[layer] * macroporeFrozenHydraulicImpedance(properties, grid, layer), .face_area_m2 = properties.boundary_face_area_m2[layer], .matrix_water_available_m3 = matrix_water, .macropore_water_available_m3 = macropore_water, .recharge_frequency_divisor = 1, .recharge_time_multiplier = boundary_face.natural_exchange_fraction, .time_fraction = properties.nonlinear_time_fraction, .source_temperature_k = grid.soil_temperature_k[layer] });
            matrix_change = boundary_face.direction_sign * oriented.matrix_water_m3;
            macropore_change = boundary_face.direction_sign * oriented.macropore_water_m3;
        } else if (topology.water_table_mode[boundary_face.cell_index] != 0) {
            const axis: usize = if (boundary_face.direction == .east or boundary_face.direction == .west) 0 else 1;
            const thickness_m = properties.vertical_thickness_m[layer];
            const midpoint_m = properties.boundary_layer_midpoint_depth_m[layer];
            const bottom_m = properties.boundary_layer_bottom_depth_m[layer];
            const profile_bottom_m = properties.boundary_layer_bottom_depth_m[boundary_face.cell_index * grid.soil_layer_capacity + grid.active_soil_layer_count[boundary_face.cell_index] - 1];
            const face_area_m2 = properties.boundary_layer_volume_m3[layer] / boundary_face.directional_layer_width_m;
            const matric_megapascal = try matricPotentialMpaAt(properties, layer, matrix_fraction);
            const base_fraction = base[layer] / properties.matrix_bulk_volume_m3[layer];
            const base_matric_megapascal = try matricPotentialMpaAt(properties, layer, base_fraction);
            const current_layer_discharge_candidate = base_matric_megapascal > grid.matric_potential_megapascal[layer];
            const macro_capacity = grid.macropore_pore_capacity_m3[layer];
            const macro_water_depth_m = if (macro_capacity > 0) bottom_m - (macropore_water + grid.macropore_ice_water_m3[layer]) / macro_capacity * thickness_m else bottom_m;
            const base_macro_water_depth_m = if (macro_capacity > 0) bottom_m - (base[cells + layer] + grid.macropore_ice_water_m3[layer]) / macro_capacity * thickness_m else bottom_m;
            for ([_]bool{ false, true }) |artificial| {
                var matrix_discharge_enabled = current_layer_discharge_candidate;
                if (artificial and topology.water_table_mode[boundary_face.cell_index] < 3) continue;
                const external_depth_m = if (artificial) topology.artificial_water_table_depth_m[boundary_face.cell_index] else topology.natural_water_table_depth_m[boundary_face.cell_index];
                const table_slope = if (artificial)
                    topology.artificial_water_table_surface_slope[boundary_face.cell_index]
                else
                    topology.natural_water_table_surface_slope[boundary_face.cell_index];
                const distance_m = if (artificial) boundary_face.artificial_water_table_distance_m else boundary_face.natural_water_table_distance_m;
                const exchange_fraction = if (artificial) boundary_face.artificial_exchange_fraction else boundary_face.natural_exchange_fraction;
                if (exchange_fraction == 0) continue;
                // WATSUB IFLGU/IFLGD require every deeper layer above the
                // corresponding external table to remain wetter than its
                // previous HOUR1 matric state. Freeze this active set from the
                // hour-start state so the implicit residual remains smooth.
                if (matrix_discharge_enabled and midpoint_m < external_depth_m) {
                    const local_layer = layer % grid.soil_layer_capacity;
                    var lower_local = local_layer + 1;
                    while (lower_local < grid.active_soil_layer_count[boundary_face.cell_index]) : (lower_local += 1) {
                        const lower = boundary_face.cell_index * grid.soil_layer_capacity + lower_local;
                        if (properties.boundary_layer_midpoint_depth_m[lower] >= external_depth_m) break;
                        const lower_fraction = base[lower] / properties.matrix_bulk_volume_m3[lower];
                        const lower_matric_megapascal = try matricPotentialMpaAt(properties, lower, lower_fraction);
                        if (lower_matric_megapascal <= grid.matric_potential_megapascal[lower] or properties.boundary_layer_midpoint_depth_m[lower] > topology.active_layer_depth_m[boundary_face.cell_index]) {
                            matrix_discharge_enabled = false;
                            break;
                        }
                    }
                }
                const fraction_below = std.math.clamp((bottom_m - external_depth_m) / thickness_m, 0, 1);
                if (midpoint_m < external_depth_m and matrix_discharge_enabled) {
                    const oriented = try water_boundary.matrixDischarge(.{ .direction_sign = boundary_face.direction_sign, .slope_sine = boundary_face.slope_sine, .directional_layer_width_m = boundary_face.directional_layer_width_m, .water_table_slope = table_slope, .matric_potential_megapascal = matric_megapascal, .saturation_water_potential_megapascal = saturationMatricPotentialMpa(properties, layer), .layer_midpoint_depth_m = midpoint_m, .external_water_table_depth_m = external_depth_m, .internal_water_table_depth_m = topology.internal_water_table_depth_m[boundary_face.cell_index], .hydraulic_conductivity_m2_per_h_megapascal = try conductivityAt(properties, layer, axis, matrix_fraction, grid.matrix_ice_water_m3[layer]), .face_area_m2 = face_area_m2, .fraction_face_below_water_table = fraction_below, .recharge_frequency_divisor = distance_m, .recharge_time_multiplier = exchange_fraction, .time_fraction = properties.nonlinear_time_fraction, .source_temperature_k = grid.soil_temperature_k[layer] }, artificial);
                    const accepted_change = limitExternalStorageChange(
                        boundary_face.direction_sign * oriented.matrix_water_m3,
                        matrix_water + matrix_change,
                        target[layer] + matrix_change,
                        grid.matrix_pore_capacity_m3[layer] -
                            grid.matrix_ice_water_m3[layer],
                    );
                    matrix_change += accepted_change;
                    if (artificial) {
                        if (properties.artificial_drainage_outflow_m3_per_step) |outflow|
                            outflow[boundary_face.cell_index] +=
                                @max(0, -accepted_change);
                    }
                } else if (!matrix_discharge_enabled and midpoint_m >= external_depth_m and midpoint_m < profile_bottom_m) {
                    const oriented = try water_boundary.recharge(.{ .direction_sign = boundary_face.direction_sign, .slope_sine = boundary_face.slope_sine, .directional_layer_width_m = boundary_face.directional_layer_width_m, .water_table_slope = table_slope, .matric_potential_megapascal = matric_megapascal, .layer_or_macropore_water_depth_m = midpoint_m, .external_water_table_depth_m = external_depth_m, .hydraulic_conductivity_m2_per_h_megapascal = try conductivityAt(properties, layer, axis, properties.retention_curve[layer].porosity_fraction, grid.matrix_ice_water_m3[layer]), .face_area_m2 = face_area_m2, .fraction_face_below_water_table = fraction_below, .recharge_frequency_divisor = @max(1, distance_m), .recharge_time_multiplier = exchange_fraction, .time_fraction = properties.nonlinear_time_fraction, .available_air_volume_m3 = @max(0, grid.matrix_pore_capacity_m3[layer] - grid.matrix_ice_water_m3[layer] - matrix_water), .source_temperature_k = grid.soil_temperature_k[layer] }, false);
                    matrix_change += limitExternalStorageChange(
                        boundary_face.direction_sign * oriented.matrix_water_m3,
                        matrix_water + matrix_change,
                        target[layer] + matrix_change,
                        grid.matrix_pore_capacity_m3[layer] -
                            grid.matrix_ice_water_m3[layer],
                    );
                }
                if (macro_capacity > 0 and base_macro_water_depth_m < external_depth_m and macropore_water > 0) {
                    // The source bound is VOLWH1*XNPXX plus the signed vertical
                    // FLWHL terms. In the whole-hour implicit solve scratch
                    // already contains those vertical terms. It additionally
                    // contains horizontal transfers and earlier boundary-face
                    // withdrawals, which is the conservative generalization
                    // required when the source sub-hour sweep is removed.
                    const oriented = try water_boundary.macroporeDischarge(.{ .direction_sign = boundary_face.direction_sign, .slope_sine = boundary_face.slope_sine, .directional_layer_width_m = boundary_face.directional_layer_width_m, .water_table_slope = table_slope, .macropore_water_depth_m = macro_water_depth_m, .external_water_table_depth_m = external_depth_m, .internal_water_table_depth_m = topology.internal_water_table_depth_m[boundary_face.cell_index], .hydraulic_conductivity_m2_per_h_megapascal = properties.boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal[layer] * macroporeFrozenHydraulicImpedance(properties, grid, layer), .face_area_m2 = face_area_m2, .fraction_face_below_water_table = fraction_below, .recharge_frequency_divisor = @max(1, distance_m), .recharge_time_multiplier = exchange_fraction, .time_fraction = properties.nonlinear_time_fraction, .available_macropore_water_m3 = macropore_water, .incoming_vertical_macropore_water_m3 = 0, .outgoing_vertical_macropore_water_m3 = 0, .source_temperature_k = grid.soil_temperature_k[layer] });
                    const accepted_change = limitExternalStorageChange(
                        boundary_face.direction_sign * oriented.macropore_water_m3,
                        macropore_water + macropore_change,
                        target[cells + layer] + macropore_change,
                        macro_capacity - grid.macropore_ice_water_m3[layer],
                    );
                    macropore_change += accepted_change;
                    if (artificial) {
                        if (properties.artificial_drainage_outflow_m3_per_step) |outflow|
                            outflow[boundary_face.cell_index] +=
                                @max(0, -accepted_change);
                    }
                } else if (macro_capacity > 0 and base_macro_water_depth_m >= external_depth_m and midpoint_m < profile_bottom_m) {
                    const oriented = try water_boundary.recharge(.{ .direction_sign = boundary_face.direction_sign, .slope_sine = boundary_face.slope_sine, .directional_layer_width_m = boundary_face.directional_layer_width_m, .water_table_slope = table_slope, .matric_potential_megapascal = 0, .layer_or_macropore_water_depth_m = macro_water_depth_m, .external_water_table_depth_m = external_depth_m, .hydraulic_conductivity_m2_per_h_megapascal = properties.boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal[layer] * macroporeFrozenHydraulicImpedance(properties, grid, layer), .face_area_m2 = face_area_m2, .fraction_face_below_water_table = fraction_below, .recharge_frequency_divisor = @max(1, distance_m), .recharge_time_multiplier = exchange_fraction, .time_fraction = properties.nonlinear_time_fraction, .available_air_volume_m3 = @max(0, macro_capacity - grid.macropore_ice_water_m3[layer] - macropore_water), .source_temperature_k = grid.soil_temperature_k[layer] }, true);
                    macropore_change += limitExternalStorageChange(
                        boundary_face.direction_sign * oriented.macropore_water_m3,
                        macropore_water + macropore_change,
                        target[cells + layer] + macropore_change,
                        macro_capacity - grid.macropore_ice_water_m3[layer],
                    );
                }
            }
        }
        target[layer] += matrix_change;
        target[cells + layer] += macropore_change;
        scratch[layer] += matrix_change;
        scratch[cells + layer] += macropore_change;
    };
    for (target, trial, residual) |value, trial_value, *difference| {
        if (!std.math.isFinite(value) or value < -1e-12)
            return error.InvalidSoilWaterCandidate;
        difference.* = @max(0.0, value) - trial_value;
    }
}

fn recordArtificialDrainage(outflow_m3: []f64, cell: usize, direction_sign: f64, oriented_flux_m3: f64) void {
    outflow_m3[cell] += @max(0, -direction_sign * oriented_flux_m3);
}

fn limitExternalStorageChange(
    proposed_change_m3: f64,
    trial_water_m3: f64,
    target_water_m3: f64,
    capacity_m3: f64,
) f64 {
    const minimum_change_m3 =
        -@min(@max(0.0, trial_water_m3), @max(0.0, target_water_m3));
    const maximum_change_m3 = @min(
        @max(0.0, capacity_m3 - trial_water_m3),
        @max(0.0, capacity_m3 - target_water_m3),
    );
    return std.math.clamp(
        proposed_change_m3,
        minimum_change_m3,
        maximum_change_m3,
    );
}

fn applyConservativeFlux(
    target: []f64,
    source: usize,
    destination: usize,
    flux_m3: f64,
) void {
    target[source] -= flux_m3;
    target[destination] += flux_m3;
}

fn limitFluxForAssembledTarget(
    flux_m3: f64,
    target: []const f64,
    source: usize,
    destination: usize,
    source_capacity_m3: f64,
    destination_capacity_m3: f64,
) f64 {
    if (flux_m3 >= 0)
        return @min(
            flux_m3,
            @min(
                @max(0.0, target[source]),
                @max(0.0, destination_capacity_m3 - target[destination]),
            ),
        );
    return @max(
        flux_m3,
        -@min(
            @max(0.0, target[destination]),
            @max(0.0, source_capacity_m3 - target[source]),
        ),
    );
}

pub fn conductivityAt(properties: Properties, cell: usize, axis: usize, water_fraction: f64, ice_water_equivalent_m3: f64) !f64 {
    _ = axis;
    return unsaturatedConductivityM2PerHMpa(.{
        .parameters = properties.mualem_van_genuchten_parameters[cell],
        .water_fraction = water_fraction,
        .ice_water_equivalent_m3 = ice_water_equivalent_m3,
        .matrix_bulk_volume_m3 = properties.matrix_bulk_volume_m3[cell],
        .conductivity_multiplier = if (properties.rainfall_conductivity_multiplier.len == 0) 1 else properties.rainfall_conductivity_multiplier[cell],
        .frozen_hydraulic_impedance_exponent = properties.frozen_hydraulic_impedance_exponent,
        .gravitational_water_potential_mpa_per_m = properties.gravitational_water_potential_mpa_per_m,
    });
}

pub const UnsaturatedConductivityInputs = struct {
    parameters: retention.MualemVanGenuchtenParameters,
    water_fraction: f64,
    ice_water_equivalent_m3: f64,
    matrix_bulk_volume_m3: f64,
    conductivity_multiplier: f64 = 1,
    frozen_hydraulic_impedance_exponent: f64 = 0,
    gravitational_water_potential_mpa_per_m: f64 = 0.00980665,
};

/// The single unsaturated-conductivity evaluation in the model: pure
/// Mualem-van Genuchten K(h) with the Richards potential unit conversion, the
/// rainfall-impact multiplier, and the frozen impedance factor. Direction is
/// carried by the anisotropy already fitted into the saturated conductivity, so
/// no axis argument is needed. Every consumer (Richards faces, dual-domain
/// exchange, root uptake resistance) must route through here so that the plant
/// and the soil never disagree about K at the same water content.
pub fn unsaturatedConductivityM2PerHMpa(inputs: UnsaturatedConductivityInputs) !f64 {
    if (inputs.matrix_bulk_volume_m3 <= 0 or
        !std.math.isFinite(inputs.water_fraction) or
        !std.math.isFinite(inputs.ice_water_equivalent_m3) or
        !std.math.isFinite(inputs.conductivity_multiplier) or
        inputs.conductivity_multiplier < 0 or
        inputs.gravitational_water_potential_mpa_per_m <= 0)
        return error.InvalidUnsaturatedConductivityInput;
    const parameters = inputs.parameters;
    const bounded_water_content = std.math.clamp(
        inputs.water_fraction,
        parameters.residual_water_content_m3_per_m3,
        parameters.saturated_water_content_m3_per_m3,
    );
    const pressure_head_m = try parameters.pressureHeadAtWaterContent(bounded_water_content);
    const conductivity_m_per_h = try parameters.hydraulicConductivityMPerH(pressure_head_m);
    const ice_content_m3_per_m3 = inputs.ice_water_equivalent_m3 / inputs.matrix_bulk_volume_m3;
    const result = conductivity_m_per_h / inputs.gravitational_water_potential_mpa_per_m *
        inputs.conductivity_multiplier *
        frozenHydraulicImpedance(inputs.frozen_hydraulic_impedance_exponent, ice_content_m3_per_m3, parameters);
    if (!std.math.isFinite(result) or result < 0) return error.NonFiniteUnsaturatedConductivity;
    return result;
}

fn matricPotentialMpaAt(properties: Properties, cell: usize, water_fraction: f64) !f64 {
    const parameters = properties.mualem_van_genuchten_parameters[cell];
    const bounded_water_content = std.math.clamp(
        water_fraction,
        parameters.residual_water_content_m3_per_m3,
        parameters.saturated_water_content_m3_per_m3,
    );
    return try parameters.pressureHeadAtWaterContent(bounded_water_content) *
        properties.gravitational_water_potential_mpa_per_m;
}

fn frozenHydraulicImpedance(
    exponent: f64,
    ice_content_m3_per_m3: f64,
    parameters: retention.MualemVanGenuchtenParameters,
) f64 {
    const available_water_content =
        parameters.saturated_water_content_m3_per_m3 -
        parameters.residual_water_content_m3_per_m3;
    const fractional_ice_content = std.math.clamp(
        ice_content_m3_per_m3 / available_water_content,
        0,
        1,
    );
    return std.math.pow(f64, 10, -exponent * fractional_ice_content);
}

fn macroporeFrozenHydraulicImpedance(
    properties: Properties,
    grid: *const grid_module.GridState,
    cell: usize,
) f64 {
    if (grid.macropore_pore_capacity_m3[cell] <= 0) return 1;
    return frozenHydraulicImpedance(
        properties.frozen_hydraulic_impedance_exponent,
        grid.macropore_ice_water_m3[cell] /
            grid.macropore_pore_capacity_m3[cell],
        properties.macropore_mualem_van_genuchten_parameters[cell],
    );
}

fn saturationMatricPotentialMpa(properties: Properties, cell: usize) f64 {
    _ = properties;
    _ = cell;
    return 0;
}

fn waterPressureMpaPerM() f64 {
    return 1000.0 * 9.80665 / 1_000_000.0;
}

fn poreCapacityRoundoffToleranceM3(capacity_m3: f64) f64 {
    return @max(
        1.0e-12,
        64.0 * std.math.floatEps(f64) *
            @max(1.0, @abs(capacity_m3)),
    );
}

fn commit(grid: *grid_module.GridState, properties: Properties, state: []const f64) !void {
    const cells = grid.layer_count;
    for (0..cells) |cell| {
        const matrix = state[cell];
        const macro = state[cells + cell];
        if (!std.math.isFinite(matrix) or !std.math.isFinite(macro) or matrix < 0 or macro < 0) return error.InvalidSoilWaterCandidate;
        if (matrix + grid.matrix_ice_water_m3[cell] > grid.matrix_pore_capacity_m3[cell] + poreCapacityRoundoffToleranceM3(grid.matrix_pore_capacity_m3[cell]) or macro + grid.macropore_ice_water_m3[cell] > grid.macropore_pore_capacity_m3[cell] + poreCapacityRoundoffToleranceM3(grid.macropore_pore_capacity_m3[cell])) return error.SoilWaterCandidateExceedsPoreCapacity;
        grid.matrix_liquid_water_m3[cell] = matrix;
        grid.macropore_liquid_water_m3[cell] = macro;
        grid.liquid_water_m3[cell] = matrix + macro;
        grid.matrix_air_volume_m3[cell] = @max(0.0, grid.matrix_pore_capacity_m3[cell] - grid.matrix_ice_water_m3[cell] - matrix);
        grid.macropore_air_volume_m3[cell] = @max(0.0, grid.macropore_pore_capacity_m3[cell] - grid.macropore_ice_water_m3[cell] - macro);
        grid.air_volume_m3[cell] = grid.matrix_air_volume_m3[cell] + grid.macropore_air_volume_m3[cell];
        grid.matric_potential_megapascal[cell] = try matricPotentialMpaAt(
            properties,
            cell,
            matrix / properties.matrix_bulk_volume_m3[cell],
        );
    }
    try grid.validateFinite();
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidSoilWaterCandidate;
        candidate.* = @max(0.0, candidate.*);
    }
}

fn projectStorageStep(grid: *const grid_module.GridState, current: []const f64, direction: []const f64, fraction: f64, output: []f64) bool {
    const cells = grid.layer_count;
    if (current.len != 2 * cells or direction.len != current.len or output.len != current.len or !std.math.isFinite(fraction)) return false;
    for (current, direction, output, 0..) |value, delta, *candidate, component| {
        const cell = if (component < cells) component else component - cells;
        const capacity = if (component < cells)
            grid.matrix_pore_capacity_m3[cell] - grid.matrix_ice_water_m3[cell]
        else
            grid.macropore_pore_capacity_m3[cell] - grid.macropore_ice_water_m3[cell];
        const unconstrained = value + fraction * delta;
        if (!std.math.isFinite(unconstrained) or !std.math.isFinite(capacity) or capacity < -1.0e-12) return false;
        candidate.* = std.math.clamp(unconstrained, 0.0, @max(0.0, capacity));
    }
    return true;
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteSoilWaterSolverState;
        maximum = @max(maximum, scaledComponentResidual(value, difference, options));
    }
    return maximum;
}

fn scaledComponentResidual(state: f64, residual: f64, options: Options) f64 {
    return @abs(residual) /
        (options.absolute_tolerance_m3 +
            options.relative_tolerance * @max(1.0, @abs(state)));
}

fn validateInputs(grid: *const grid_module.GridState, faces: []const Face, properties: Properties, micro_flux: []const f64, macro_flux: []const f64, options: Options) !void {
    const cells = grid.layer_count;
    if (properties.matrix_bulk_volume_m3.len != cells or properties.retention_curve.len != cells or properties.mualem_van_genuchten_parameters.len != cells or properties.macropore_mualem_van_genuchten_parameters.len != cells or properties.dual_domain_exchange_enabled.len != cells or properties.macropore_spacing_m.len != cells or properties.macropore_radius_m.len != cells or properties.gravitational_potential_megapascal.len != cells or properties.osmotic_potential_megapascal.len != cells or properties.vertical_thickness_m.len != cells or (properties.rainfall_conductivity_multiplier.len != 0 and properties.rainfall_conductivity_multiplier.len != cells) or (properties.matrix_external_source_m3_per_step.len != 0 and properties.matrix_external_source_m3_per_step.len != cells) or micro_flux.len != faces.len or macro_flux.len != faces.len) return error.SoilWaterSolverDimensionMismatch;
    for (properties.mualem_van_genuchten_parameters) |parameters| try parameters.validate();
    for (properties.macropore_mualem_van_genuchten_parameters) |parameters| try parameters.validate();
    if (!std.math.isFinite(properties.frozen_hydraulic_impedance_exponent) or properties.frozen_hydraulic_impedance_exponent < 0 or !std.math.isFinite(properties.gravitational_water_potential_mpa_per_m) or properties.gravitational_water_potential_mpa_per_m <= 0) return error.InvalidSoilWaterSolverProperty;
    if (properties.boundary_topology != null and (properties.boundary_face_area_m2.len != cells or properties.boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal.len != cells or properties.boundary_layer_volume_m3.len != cells or properties.boundary_layer_midpoint_depth_m.len != cells or properties.boundary_layer_bottom_depth_m.len != cells)) return error.SoilWaterBoundaryDimensionMismatch;
    if (options.max_iterations == 0 or !std.math.isFinite(options.absolute_tolerance_m3) or options.absolute_tolerance_m3 <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or !std.math.isFinite(properties.osmotic_potential_multiplier) or !std.math.isFinite(properties.dual_domain_geometry_factor) or properties.dual_domain_geometry_factor <= 0 or !std.math.isFinite(properties.dual_domain_scaling_coefficient) or properties.dual_domain_scaling_coefficient <= 0) return error.InvalidSoilWaterSolverOptions;
    for (0..cells) |cell| {
        if (!std.math.isFinite(properties.matrix_bulk_volume_m3[cell]) or properties.matrix_bulk_volume_m3[cell] <= 0 or !std.math.isFinite(properties.vertical_thickness_m[cell]) or properties.vertical_thickness_m[cell] <= 0 or grid.matrix_liquid_water_m3[cell] < 0 or grid.macropore_liquid_water_m3[cell] < 0 or grid.matrix_liquid_water_m3[cell] + grid.matrix_ice_water_m3[cell] > grid.matrix_pore_capacity_m3[cell] + poreCapacityRoundoffToleranceM3(grid.matrix_pore_capacity_m3[cell]) or grid.macropore_liquid_water_m3[cell] + grid.macropore_ice_water_m3[cell] > grid.macropore_pore_capacity_m3[cell] + poreCapacityRoundoffToleranceM3(grid.macropore_pore_capacity_m3[cell])) {
            if (!builtin.is_test) std.log.err("invalid soil water input: layer_cell={d} bulk_m3={e} thickness_m={e} matrix_water_m3={e} matrix_ice_m3={e} matrix_capacity_m3={e} macropore_water_m3={e} macropore_ice_m3={e} macropore_capacity_m3={e}", .{ cell, properties.matrix_bulk_volume_m3[cell], properties.vertical_thickness_m[cell], grid.matrix_liquid_water_m3[cell], grid.matrix_ice_water_m3[cell], grid.matrix_pore_capacity_m3[cell], grid.macropore_liquid_water_m3[cell], grid.macropore_ice_water_m3[cell], grid.macropore_pore_capacity_m3[cell] });
            return error.InvalidSoilWaterSolverInput;
        }
    }
    for (faces) |face| if (face.source_cell >= cells or face.destination_cell >= cells or face.source_cell == face.destination_cell or !std.math.isFinite(face.source_path_length_m) or face.source_path_length_m <= 0 or !std.math.isFinite(face.destination_path_length_m) or face.destination_path_length_m <= 0 or !std.math.isFinite(face.face_area_m2) or face.face_area_m2 < 0 or !std.math.isFinite(face.macropore_hydraulic_conductance_m_per_h_megapascal) or face.macropore_hydraulic_conductance_m_per_h_megapascal < 0) return error.InvalidSoilWaterFace;
}

fn testCurve() retention.ResolvedCurve {
    return .{ .porosity_fraction = 0.5, .curve = .{ .field_capacity_fraction = 0.3, .wilting_point_fraction = 0.1, .saturation_water_potential_megapascal = -0.0005, .field_capacity_water_potential_megapascal = -0.01, .wilting_point_water_potential_megapascal = -1.5, .minimum_water_potential_megapascal = -1.5e12, .saturation_to_field_shape = 0.5, .below_wilting_shape = 0.5 } };
}

fn testMatrixMualemVanGenuchten() retention.MualemVanGenuchtenParameters {
    return .{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.5,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.002,
    };
}

fn testMacroporeMualemVanGenuchten() retention.MualemVanGenuchtenParameters {
    return .{
        .residual_water_content_m3_per_m3 = 0,
        .saturated_water_content_m3_per_m3 = 1,
        .alpha_per_m = 15,
        .n = 2.68,
        .saturated_hydraulic_conductivity_m_per_h = 0.1,
    };
}

test "runtime NPH hybrid water solve exits early and conserves both pore domains" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.35;
    grid.matrix_liquid_water_m3[1] = 0.15;
    grid.macropore_liquid_water_m3[0] = 0.08;
    grid.macropore_liquid_water_m3[1] = 0.02;
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    @memset(grid.macropore_pore_capacity_m3, 0.1);
    const curves = [_]retention.ResolvedCurve{ testCurve(), testCurve() };
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{ testMatrixMualemVanGenuchten(), testMatrixMualemVanGenuchten() };
    const macropore_parameters = [_]retention.MualemVanGenuchtenParameters{ testMacroporeMualemVanGenuchten(), testMacroporeMualemVanGenuchten() };
    const macropore_spacing_m = [_]f64{ 0.2, 0.2 };
    const macropore_radius_m = [_]f64{ 0.001, 0.001 };
    const dual_domain_disabled = [_]bool{ false, false };
    const bulk = [_]f64{ 1, 1 };
    const gravity = [_]f64{ 0, 0 };
    const osmotic = [_]f64{ 0, 0 };
    const conductivity = [_]f64{0.002} ** 6;
    const thickness = [_]f64{ 0.1, 0.1 };
    var micro_flux = [_]f64{0};
    var macro_flux = [_]f64{0};
    const before_micro = grid.matrix_liquid_water_m3[0] + grid.matrix_liquid_water_m3[1];
    const before_macro = grid.macropore_liquid_water_m3[0] + grid.macropore_liquid_water_m3[1];
    const result = try solve(std.testing.allocator, &grid, &.{.{ .source_cell = 0, .destination_cell = 1, .direction = .horizontal, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .macropore_hydraulic_conductance_m_per_h_megapascal = 0.01 }}, .{ .matrix_bulk_volume_m3 = &bulk, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &dual_domain_disabled, .gravitational_potential_megapascal = &gravity, .osmotic_potential_megapascal = &osmotic, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1 }, &micro_flux, &macro_flux, .{ .max_iterations = 40 });
    try std.testing.expect(result.iterations < 40);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
    try std.testing.expectApproxEqAbs(before_micro, grid.matrix_liquid_water_m3[0] + grid.matrix_liquid_water_m3[1], 1e-12);
    try std.testing.expectApproxEqAbs(before_macro, grid.macropore_liquid_water_m3[0] + grid.macropore_liquid_water_m3[1], 1e-12);
    try std.testing.expect(micro_flux[0] > 0);
    try std.testing.expect(macro_flux[0] > 0);
}

test "subsurface irrigation is an external source in the Richards residual" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 1,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-12,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.2;
    grid.matrix_pore_capacity_m3[0] = 0.5;
    grid.macropore_liquid_water_m3[0] = 0;
    grid.macropore_pore_capacity_m3[0] = 0;
    const properties: Properties = .{
        .matrix_bulk_volume_m3 = &.{1},
        .retention_curve = &.{testCurve()},
        .mualem_van_genuchten_parameters = &.{testMatrixMualemVanGenuchten()},
        .macropore_mualem_van_genuchten_parameters = &.{testMacroporeMualemVanGenuchten()},
        .macropore_spacing_m = &.{1},
        .macropore_radius_m = &.{0},
        .dual_domain_exchange_enabled = &.{false},
        .gravitational_potential_megapascal = &.{0},
        .osmotic_potential_megapascal = &.{0},
        .matrix_hydraulic_conductivity_m2_per_h_megapascal = &.{ 0, 0, 0 },
        .matrix_external_source_m3_per_step = &.{0.03},
        .vertical_thickness_m = &.{1},
        .osmotic_potential_multiplier = 1,
    };
    var matrix_flux: [0]f64 = .{};
    var macropore_flux: [0]f64 = .{};
    const result = try solve(
        std.testing.allocator,
        &grid,
        &.{},
        properties,
        &matrix_flux,
        &macropore_flux,
        .{
            .max_iterations = 20,
            .absolute_tolerance_m3 = 1e-12,
            .relative_tolerance = 1e-10,
        },
    );
    try std.testing.expect(result.iterations < 20);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.23),
        grid.matrix_liquid_water_m3[0],
        1e-11,
    );
}

test "simultaneous Richards residual is independent of face ordering" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 3,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 3 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    const initial_matrix_water_m3 = [_]f64{ 0.38, 0.24, 0.12 };
    @memcpy(grid.matrix_liquid_water_m3, &initial_matrix_water_m3);
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    @memset(grid.macropore_liquid_water_m3, 0);
    @memset(grid.macropore_pore_capacity_m3, 0);
    const curves = [_]retention.ResolvedCurve{
        testCurve(),
        testCurve(),
        testCurve(),
    };
    const matrix_parameters =
        [_]retention.MualemVanGenuchtenParameters{
            testMatrixMualemVanGenuchten(),
            testMatrixMualemVanGenuchten(),
            testMatrixMualemVanGenuchten(),
        };
    const macropore_parameters =
        [_]retention.MualemVanGenuchtenParameters{
            testMacroporeMualemVanGenuchten(),
            testMacroporeMualemVanGenuchten(),
            testMacroporeMualemVanGenuchten(),
        };
    const scalar = [_]f64{ 1, 1, 1 };
    const zero = [_]f64{ 0, 0, 0 };
    const disabled = [_]bool{ false, false, false };
    const conductivity = [_]f64{0.002} ** 9;
    const properties: Properties = .{
        .matrix_bulk_volume_m3 = &scalar,
        .retention_curve = &curves,
        .mualem_van_genuchten_parameters = &matrix_parameters,
        .macropore_mualem_van_genuchten_parameters = &macropore_parameters,
        .macropore_spacing_m = &scalar,
        .macropore_radius_m = &zero,
        .dual_domain_exchange_enabled = &disabled,
        .gravitational_potential_megapascal = &zero,
        .osmotic_potential_megapascal = &zero,
        .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity,
        .vertical_thickness_m = &scalar,
        .osmotic_potential_multiplier = 1,
    };
    const face_01: Face = .{
        .source_cell = 0,
        .destination_cell = 1,
        .direction = .horizontal,
        .source_path_length_m = 1,
        .destination_path_length_m = 1,
        .face_area_m2 = 1,
                .macropore_hydraulic_conductance_m_per_h_megapascal = 0,
    };
    const face_12: Face = .{
        .source_cell = 1,
        .destination_cell = 2,
        .direction = .horizontal,
        .source_path_length_m = 1,
        .destination_path_length_m = 1,
        .face_area_m2 = 1,
                .macropore_hydraulic_conductance_m_per_h_megapascal = 0,
    };
    const trial = [_]f64{ 0.38, 0.24, 0.12, 0, 0, 0 };
    var target_a = [_]f64{0} ** 6;
    var residual_a = [_]f64{0} ** 6;
    var scratch_a = [_]f64{0} ** 6;
    var matrix_flux_a = [_]f64{0} ** 2;
    var macropore_flux_a = [_]f64{0} ** 2;
    try residualAt(
        &grid,
        &.{ face_01, face_12 },
        properties,
        &trial,
        &trial,
        &target_a,
        &residual_a,
        &scratch_a,
        &matrix_flux_a,
        &macropore_flux_a,
    );
    var target_b = [_]f64{0} ** 6;
    var residual_b = [_]f64{0} ** 6;
    var scratch_b = [_]f64{0} ** 6;
    var matrix_flux_b = [_]f64{0} ** 2;
    var macropore_flux_b = [_]f64{0} ** 2;
    try residualAt(
        &grid,
        &.{ face_12, face_01 },
        properties,
        &trial,
        &trial,
        &target_b,
        &residual_b,
        &scratch_b,
        &matrix_flux_b,
        &macropore_flux_b,
    );
    for (residual_a, residual_b) |first, second|
        try std.testing.expectApproxEqAbs(first, second, 1e-15);
}

test "lower drainage and lateral water-table recharge converge inside NPH residual" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.4;
    grid.macropore_liquid_water_m3[0] = 0.05;
    grid.matrix_pore_capacity_m3[0] = 0.5;
    grid.macropore_pore_capacity_m3[0] = 0.1;
    var site = try @import("../../state/site.zig").parse(std.testing.allocator, @import("../../core/test_fixtures.zig").site_source, 1, 1);
    defer site.deinit();
    var topography = try @import("../../state/topography.zig").parse(std.testing.allocator, "1 1 1 1 0 0 0 0\nsoil\n");
    defer topography.deinit();
    var terrain = try @import("../../state/terrain_hydrology.zig").State.initMapped(std.testing.allocator, topography, &.{0}, &.{1}, &.{1}, 1, 1);
    defer terrain.deinit();
    var topology = try boundary_topology.State.initMapped(std.testing.allocator, &grid, &terrain, 1, 1, &.{1}, &.{1}, &.{site});
    defer topology.deinit();
    for (topology.faces) |*boundary_face| if (boundary_face.is_lower_boundary) {
        boundary_face.natural_exchange_fraction = 1;
    };
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{testMatrixMualemVanGenuchten()};
    const macropore_parameters = [_]retention.MualemVanGenuchtenParameters{testMacroporeMualemVanGenuchten()};
    const macropore_spacing_m = [_]f64{0.2};
    const macropore_radius_m = [_]f64{0.001};
    const dual_domain_disabled = [_]bool{false};
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const thickness = [_]f64{0.1};
    const conductivity = [_]f64{ 0.002, 0.002, 0.002 };
    const macro_conductivity = [_]f64{0.001};
    const before_matrix = grid.matrix_liquid_water_m3[0];
    const before_macro = grid.macropore_liquid_water_m3[0];
    _ = try solve(std.testing.allocator, &grid, &.{}, .{ .matrix_bulk_volume_m3 = &one, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &dual_domain_disabled, .gravitational_potential_megapascal = &zero, .osmotic_potential_megapascal = &zero, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1, .boundary_topology = &topology, .boundary_face_area_m2 = &one, .boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal = &macro_conductivity, .boundary_layer_volume_m3 = &one, .boundary_layer_midpoint_depth_m = &thickness, .boundary_layer_bottom_depth_m = &one }, &.{}, &.{}, .{ .max_iterations = 20 });
    try std.testing.expect(grid.matrix_liquid_water_m3[0] < before_matrix);
    try std.testing.expect(grid.macropore_liquid_water_m3[0] < before_macro);
    try std.testing.expect(grid.matrix_liquid_water_m3[0] >= 0);
    try std.testing.expect(grid.macropore_liquid_water_m3[0] >= 0);

    // Reuse the self-contained state to exercise lateral water-table recharge
    // through all four perimeter faces of the one-cell runtime grid.
    topology.natural_water_table_depth_m[0] = 0.01;
    topology.water_table_mode[0] = 1;
    for (topology.faces) |*boundary_face| {
        boundary_face.natural_exchange_fraction = if (boundary_face.is_lower_boundary) 0 else 1;
        boundary_face.artificial_exchange_fraction = 0;
    }
    grid.matrix_liquid_water_m3[0] = 0.1;
    grid.macropore_liquid_water_m3[0] = 0.01;
    _ = try solve(std.testing.allocator, &grid, &.{}, .{ .matrix_bulk_volume_m3 = &one, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &dual_domain_disabled, .gravitational_potential_megapascal = &zero, .osmotic_potential_megapascal = &zero, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1, .boundary_topology = &topology, .boundary_face_area_m2 = &one, .boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal = &macro_conductivity, .boundary_layer_volume_m3 = &one, .boundary_layer_midpoint_depth_m = &thickness, .boundary_layer_bottom_depth_m = &one }, &.{}, &.{}, .{ .max_iterations = 20 });
    try std.testing.expect(grid.matrix_liquid_water_m3[0] > 0.1);
    try std.testing.expect(grid.macropore_liquid_water_m3[0] > 0.01);

    // Artificial drainage is retained separately for OUTSD UVOLY and is
    // published only from the converged residual, never from trial iterates.
    topology.water_table_mode[0] = 3;
    topology.artificial_water_table_depth_m[0] = 0.2;
    for (topology.faces) |*boundary_face| {
        boundary_face.natural_exchange_fraction = 0;
        boundary_face.artificial_exchange_fraction = if (boundary_face.is_lower_boundary) 0 else 1;
        boundary_face.artificial_water_table_distance_m = 1;
        if (!boundary_face.is_lower_boundary) {
            boundary_face.direction_sign = 1;
            boundary_face.slope_sine = 1;
        }
    }
    grid.matrix_liquid_water_m3[0] = 0.499;
    grid.macropore_liquid_water_m3[0] = 0.05;
    grid.matric_potential_megapascal[0] = -1.0e300;
    var artificial_drainage = [_]f64{0};
    const bottom = [_]f64{0.3};
    _ = try solve(std.testing.allocator, &grid, &.{}, .{ .matrix_bulk_volume_m3 = &one, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &dual_domain_disabled, .gravitational_potential_megapascal = &zero, .osmotic_potential_megapascal = &zero, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1, .boundary_topology = &topology, .boundary_face_area_m2 = &one, .boundary_macropore_hydraulic_conductivity_m2_per_h_megapascal = &macro_conductivity, .boundary_layer_volume_m3 = &one, .boundary_layer_midpoint_depth_m = &thickness, .boundary_layer_bottom_depth_m = &bottom, .artificial_drainage_outflow_m3_per_step = &artificial_drainage }, &.{}, &.{}, .{ .max_iterations = 20 });
    try std.testing.expect(artificial_drainage[0] >= 0);
}

test "artificial drain publication preserves REDIST UVOLY loss sign" {
    var outflow = [_]f64{0};
    recordArtificialDrainage(&outflow, 0, 1, -0.25);
    recordArtificialDrainage(&outflow, 0, -1, 0.5);
    recordArtificialDrainage(&outflow, 0, 1, 0.1);
    try std.testing.expectEqual(@as(f64, 0.75), outflow[0]);
}

test "reduced Newton coordinates remove pore-domain conservation nullspaces" {
    var jacobian = [_]f64{
        -1, 1,  0,  0,
        1,  -1, 0,  0,
        0,  0,  -1, 1,
        0,  0,  1,  -1,
    };
    const residual = [_]f64{ 1, -1, 2, -2 };
    var delta = [_]f64{0} ** 4;
    var reduced_jacobian = [_]f64{0} ** 16;
    var reduced_rhs = [_]f64{0} ** 4;
    var reduced_index = [_]usize{0} ** 4;
    var micro_parent = [_]usize{ 0, 0 };
    var macro_parent = [_]usize{ 0, 0 };
    const component_size = [_]usize{ 2, 1 };
    try std.testing.expect(solveConservedNewtonSystem(&jacobian, &residual, &delta, &reduced_jacobian, &reduced_rhs, &reduced_index, 4, 2, &micro_parent, &component_size, &macro_parent, &component_size));
    try std.testing.expectApproxEqAbs(@as(f64, 0), delta[0] + delta[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0), delta[2] + delta[3], 1e-15);
    for (0..4) |row| {
        var linearized_residual = residual[row];
        for (0..4) |column| linearized_residual += jacobian[row * 4 + column] * delta[column];
        try std.testing.expectApproxEqAbs(@as(f64, 0), linearized_residual, 1e-12);
    }
    @memset(&delta, 0);
    @memset(&reduced_jacobian, 0);
    @memset(&reduced_rhs, 0);
    micro_parent = .{ 0, 0 };
    macro_parent = .{ 0, 0 };
    try std.testing.expect(solveConservedTrustRegionSystem(&jacobian, &residual, &.{ 1, 1, 1, 1 }, .{ .max_iterations = 1 }, &delta, &reduced_jacobian, &reduced_rhs, &reduced_index, 4, 2, &micro_parent, &component_size, &macro_parent, &component_size, 1e-8));
    try std.testing.expectApproxEqAbs(@as(f64, 0), delta[0] + delta[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0), delta[2] + delta[3], 1e-15);
}

test "whole-step residual retains every independent storage coordinate" {
    var parent = [_]usize{ 0, 0 };
    var size = [_]usize{ 2, 1 };
    releaseIndependentStorageCoordinates(&parent, &size);
    try std.testing.expectEqual(@as(usize, 1), size[0]);
}

test "rejected water solve leaves grid and published flux unchanged" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.35;
    grid.matrix_liquid_water_m3[1] = 0.15;
    grid.matrix_ice_water_m3[0] = 0.2;
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    const curves = [_]retention.ResolvedCurve{ testCurve(), testCurve() };
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{ testMatrixMualemVanGenuchten(), testMatrixMualemVanGenuchten() };
    const macropore_parameters = [_]retention.MualemVanGenuchtenParameters{ testMacroporeMualemVanGenuchten(), testMacroporeMualemVanGenuchten() };
    const macropore_spacing_m = [_]f64{ 0.2, 0.2 };
    const macropore_radius_m = [_]f64{ 0.001, 0.001 };
    const dual_domain_disabled = [_]bool{ false, false };
    const pair = [_]f64{ 1, 1 };
    const zeros = [_]f64{ 0, 0 };
    const conductivity = [_]f64{0.002} ** 6;
    const thickness = [_]f64{ 0.1, 0.1 };
    var micro_flux = [_]f64{99};
    var macro_flux = [_]f64{88};
    try std.testing.expectError(error.InvalidSoilWaterSolverInput, solve(std.testing.allocator, &grid, &.{.{ .source_cell = 0, .destination_cell = 1, .direction = .horizontal, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .macropore_hydraulic_conductance_m_per_h_megapascal = 0 }}, .{ .matrix_bulk_volume_m3 = &pair, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &dual_domain_disabled, .gravitational_potential_megapascal = &zeros, .osmotic_potential_megapascal = &zeros, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1 }, &micro_flux, &macro_flux, .{ .max_iterations = 1, .absolute_tolerance_m3 = 1.0e-30, .relative_tolerance = 1.0e-30 }));
    try std.testing.expectEqual(@as(f64, 0.35), grid.matrix_liquid_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 99), micro_flux[0]);
    try std.testing.expectEqual(@as(f64, 88), macro_flux[0]);
}

test "converged WATSUB flux binds to shared hydrology and solute faces" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.35;
    grid.matrix_liquid_water_m3[1] = 0.15;
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    grid.matrix_air_volume_m3[0] = 0.15;
    grid.matrix_air_volume_m3[1] = 0.35;
    @memcpy(grid.air_volume_m3, grid.matrix_air_volume_m3);
    var snow = try @import("../solute/snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try transport_hydrology.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    try hydrology.syncStorage(&grid, &snow);
    var shared_faces = try transport_hydrology.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer shared_faces.deinit();
    const curves = [_]retention.ResolvedCurve{ testCurve(), testCurve() };
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{ testMatrixMualemVanGenuchten(), testMatrixMualemVanGenuchten() };
    const macropore_parameters = [_]retention.MualemVanGenuchtenParameters{ testMacroporeMualemVanGenuchten(), testMacroporeMualemVanGenuchten() };
    const macropore_spacing_m = [_]f64{ 0.2, 0.2 };
    const macropore_radius_m = [_]f64{ 0.001, 0.001 };
    const dual_domain_disabled = [_]bool{ false, false };
    const bulk = [_]f64{ 1, 1 };
    const zeros = [_]f64{ 0, 0 };
    const conductivity = [_]f64{0.002} ** 6;
    const thickness = [_]f64{ 0.1, 0.1 };
    const one = [_]f64{1};
    const no_macro = [_]f64{0};
    _ = try solveAndBindTransportFaces(std.testing.allocator, &grid, &hydrology, &shared_faces, .{ .source_path_length_m = &one, .destination_path_length_m = &one, .face_area_m2 = &one, .macropore_hydraulic_conductance_m_per_h_megapascal = &no_macro }, .{ .matrix_bulk_volume_m3 = &bulk, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &dual_domain_disabled, .gravitational_potential_megapascal = &zeros, .osmotic_potential_megapascal = &zeros, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1 }, .{ .max_iterations = 40 });
    try std.testing.expectEqual(@as(u2, 0), shared_faces.direction_axis[0]);
    try std.testing.expect(shared_faces.micropore_faces[0].water_flux_m3_per_step > 0);
    try std.testing.expectEqual(shared_faces.micropore_faces[0].water_flux_m3_per_step, hydrology.micropore_face_flux_m3_per_step[0]);
    try std.testing.expectEqualSlices(f64, grid.matrix_liquid_water_m3, hydrology.micropore_water_volume_m3);
}

test "original Mualem conductivity is isotropic before runtime directional geometry" {
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const scalar = [_]f64{1};
    const parameters = [_]retention.MualemVanGenuchtenParameters{testMatrixMualemVanGenuchten()};
    const compatibility_table = [_]f64{ 10, 20, 30 };
    const properties: Properties = .{ .matrix_bulk_volume_m3 = &scalar, .retention_curve = &curves, .mualem_van_genuchten_parameters = &parameters, .gravitational_potential_megapascal = &scalar, .osmotic_potential_megapascal = &scalar, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &compatibility_table, .vertical_thickness_m = &scalar, .osmotic_potential_multiplier = 1 };
    const x_conductivity = try conductivityAt(properties, 0, 0, 0.3, 0);
    try std.testing.expectEqual(x_conductivity, try conductivityAt(properties, 0, 1, 0.3, 0));
    try std.testing.expectEqual(x_conductivity, try conductivityAt(properties, 0, 2, 0.3, 0));
}

test "root uptake and Richards faces share one unsaturated conductivity" {
    // Plant root uptake resistance calls `unsaturatedConductivityM2PerHMpa`
    // directly while the Richards residual calls `conductivityAt`. If those two
    // ever diverge, the plant and the soil disagree about how fast water reaches
    // a root at the same water content, and the disagreement shows up only as a
    // slow water-balance drift. Pin them to bit equality across the curve.
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const scalar = [_]f64{1};
    const multiplier = [_]f64{0.25};
    const parameters = [_]retention.MualemVanGenuchtenParameters{testMatrixMualemVanGenuchten()};
    const compatibility_table = [_]f64{ 10, 20, 30 };
    const properties: Properties = .{
        .matrix_bulk_volume_m3 = &scalar,
        .retention_curve = &curves,
        .mualem_van_genuchten_parameters = &parameters,
        .gravitational_potential_megapascal = &scalar,
        .osmotic_potential_megapascal = &scalar,
        .matrix_hydraulic_conductivity_m2_per_h_megapascal = &compatibility_table,
        .rainfall_conductivity_multiplier = &multiplier,
        .vertical_thickness_m = &scalar,
        .osmotic_potential_multiplier = 1,
        .frozen_hydraulic_impedance_exponent = 7,
    };
    for ([_]f64{ 0.05, 0.12, 0.2, 0.28, 0.35, 0.42 }) |water_fraction| {
        for ([_]f64{ 0, 0.03, 0.1 }) |ice_m3| {
            const solver_conductivity = try conductivityAt(properties, 0, 2, water_fraction, ice_m3);
            const plant_conductivity = try unsaturatedConductivityM2PerHMpa(.{
                .parameters = parameters[0],
                .water_fraction = water_fraction,
                .ice_water_equivalent_m3 = ice_m3,
                .matrix_bulk_volume_m3 = scalar[0],
                .conductivity_multiplier = multiplier[0],
                .frozen_hydraulic_impedance_exponent = properties.frozen_hydraulic_impedance_exponent,
                .gravitational_water_potential_mpa_per_m = properties.gravitational_water_potential_mpa_per_m,
            });
            try std.testing.expectEqual(solver_conductivity, plant_conductivity);
            try std.testing.expect(plant_conductivity >= 0);
        }
    }
    try std.testing.expectError(error.InvalidUnsaturatedConductivityInput, unsaturatedConductivityM2PerHMpa(.{
        .parameters = parameters[0],
        .water_fraction = 0.3,
        .ice_water_equivalent_m3 = 0,
        .matrix_bulk_volume_m3 = 0,
    }));
}

test "rainfall damage scales Mualem conductivity without changing parameters" {
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const scalar = [_]f64{1};
    const multiplier = [_]f64{0.25};
    const parameters = [_]retention.MualemVanGenuchtenParameters{testMatrixMualemVanGenuchten()};
    const compatibility_table = [_]f64{ 8, 8, 8 };
    const unscaled: Properties = .{ .matrix_bulk_volume_m3 = &scalar, .retention_curve = &curves, .mualem_van_genuchten_parameters = &parameters, .gravitational_potential_megapascal = &scalar, .osmotic_potential_megapascal = &scalar, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &compatibility_table, .vertical_thickness_m = &scalar, .osmotic_potential_multiplier = 1 };
    var scaled = unscaled;
    scaled.rainfall_conductivity_multiplier = &multiplier;
    try std.testing.expectApproxEqAbs(0.25 * try conductivityAt(unscaled, 0, 0, 0.3, 0), try conductivityAt(scaled, 0, 0, 0.3, 0), 1e-15);
    try std.testing.expectEqual(testMatrixMualemVanGenuchten(), parameters[0]);
}

test "Richards residual uses runtime Mualem van Genuchten head and conductivity" {
    const parameters: retention.MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.45,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.012,
    };
    const water_content_m3_per_m3 = 0.25;
    const pressure_head_m =
        try parameters.pressureHeadAtWaterContent(water_content_m3_per_m3);
    const expected_conductivity_m_per_h =
        try parameters.hydraulicConductivityMPerH(pressure_head_m);
    const scalar = [_]f64{1};
    const curve = [_]retention.ResolvedCurve{testCurve()};
    const parameter_slice = [_]retention.MualemVanGenuchtenParameters{parameters};
    const compatibility_table = [_]f64{ 99, 99, 99 };
    const properties: Properties = .{
        .matrix_bulk_volume_m3 = &scalar,
        .retention_curve = &curve,
        .mualem_van_genuchten_parameters = &parameter_slice,
        .gravitational_potential_megapascal = &scalar,
        .osmotic_potential_megapascal = &scalar,
        .matrix_hydraulic_conductivity_m2_per_h_megapascal = &compatibility_table,
        .vertical_thickness_m = &scalar,
        .osmotic_potential_multiplier = 1,
    };
    try std.testing.expectApproxEqAbs(
        pressure_head_m * waterPressureMpaPerM(),
        try matricPotentialMpaAt(properties, 0, water_content_m3_per_m3),
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        expected_conductivity_m_per_h,
        try conductivityAt(properties, 0, 2, water_content_m3_per_m3, 0) *
            waterPressureMpaPerM(),
        1.0e-14,
    );
}

test "Dall'Amico ice impedance scales runtime Mualem conductivity" {
    const parameters: retention.MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.45,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.012,
    };
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const scalar = [_]f64{1};
    const zero = [_]f64{0};
    const ice = [_]f64{0.2};
    const conductivity = [_]f64{1} ** 3;
    const parameters_slice = [_]retention.MualemVanGenuchtenParameters{parameters};
    const unfrozen: Properties = .{ .matrix_bulk_volume_m3 = &scalar, .retention_curve = &curves, .mualem_van_genuchten_parameters = &parameters_slice, .gravitational_potential_megapascal = &zero, .osmotic_potential_megapascal = &zero, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &scalar, .osmotic_potential_multiplier = 1, .frozen_hydraulic_impedance_exponent = 7 };
    const water_content_m3_per_m3 = 0.25;
    const unfrozen_conductivity =
        try conductivityAt(unfrozen, 0, 2, water_content_m3_per_m3, 0);
    const frozen_conductivity =
        try conductivityAt(unfrozen, 0, 2, water_content_m3_per_m3, ice[0]);
    // q = 0.2 / (0.45 - 0.05) = 0.5; impedance = 10^(-7q).
    try std.testing.expectApproxEqRel(
        unfrozen_conductivity * std.math.pow(f64, 10, -3.5),
        frozen_conductivity,
        1e-12,
    );
}

test "macropore Richards flow uses runtime shape and conserves pore water" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(
        .{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 2 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 40 },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    @memset(grid.matrix_liquid_water_m3, 0.25);
    @memset(grid.macropore_pore_capacity_m3, 0.1);
    grid.macropore_liquid_water_m3[0] = 0.08;
    grid.macropore_liquid_water_m3[1] = 0.02;
    const before = grid.macropore_liquid_water_m3[0] +
        grid.macropore_liquid_water_m3[1];
    const curves = [_]retention.ResolvedCurve{ testCurve(), testCurve() };
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{ testMatrixMualemVanGenuchten(), testMatrixMualemVanGenuchten() };
    const macro_parameters = [_]retention.MualemVanGenuchtenParameters{
        .{ .residual_water_content_m3_per_m3 = 0, .saturated_water_content_m3_per_m3 = 1, .alpha_per_m = 15, .n = 2.68, .saturated_hydraulic_conductivity_m_per_h = 0.1 },
        .{ .residual_water_content_m3_per_m3 = 0, .saturated_water_content_m3_per_m3 = 1, .alpha_per_m = 15, .n = 2.68, .saturated_hydraulic_conductivity_m_per_h = 0.1 },
    };
    const pair = [_]f64{ 1, 1 };
    const zeros = [_]f64{ 0, 0 };
    const spacing = [_]f64{ 0.2, 0.2 };
    const radius = [_]f64{ 0.001, 0.001 };
    const exchange_disabled = [_]bool{ false, false };
    const conductivity = [_]f64{0.002} ** 6;
    var matrix_flux = [_]f64{0};
    var macropore_flux = [_]f64{0};
    _ = try solve(
        std.testing.allocator,
        &grid,
        &.{.{ .source_cell = 0, .destination_cell = 1, .direction = .horizontal, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .macropore_hydraulic_conductance_m_per_h_megapascal = 999 }},
        .{ .matrix_bulk_volume_m3 = &pair, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macro_parameters, .macropore_spacing_m = &spacing, .macropore_radius_m = &radius, .dual_domain_exchange_enabled = &exchange_disabled, .gravitational_potential_megapascal = &zeros, .osmotic_potential_megapascal = &zeros, .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity, .vertical_thickness_m = &pair, .osmotic_potential_multiplier = 1 },
        &matrix_flux,
        &macropore_flux,
        .{ .max_iterations = 40 },
    );
    try std.testing.expect(macropore_flux[0] > 0);
    try std.testing.expectApproxEqAbs(
        before,
        grid.macropore_liquid_water_m3[0] + grid.macropore_liquid_water_m3[1],
        1e-12,
    );
    // The obsolete face conductance is deliberately extreme; the accepted
    // flux remains donor-bounded because the runtime Mualem parameters and
    // CNDH-derived Ksat control this branch.
    try std.testing.expect(macropore_flux[0] <= 0.08);
}

test "Gerke van Genuchten exchange converges inside the water residual" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 40 },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_pore_capacity_m3[0] = 0.45;
    grid.matrix_liquid_water_m3[0] = 0.10;
    grid.macropore_pore_capacity_m3[0] = 0.10;
    grid.macropore_liquid_water_m3[0] = 0.08;
    const before = grid.matrix_liquid_water_m3[0] +
        grid.macropore_liquid_water_m3[0];
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{.{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.45,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.01,
    }};
    const macropore_parameters = [_]retention.MualemVanGenuchtenParameters{.{
        .residual_water_content_m3_per_m3 = 0,
        .saturated_water_content_m3_per_m3 = 1,
        .alpha_per_m = 15,
        .n = 2.68,
        .saturated_hydraulic_conductivity_m_per_h = 0.1,
    }};
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const spacing = [_]f64{0.2};
    const radius = [_]f64{0.001};
    const enabled = [_]bool{true};
    const conductivity = [_]f64{ 0.002, 0.002, 0.002 };
    const result = try solve(
        std.testing.allocator,
        &grid,
        &.{},
        .{
            .matrix_bulk_volume_m3 = &one,
            .retention_curve = &curves,
            .mualem_van_genuchten_parameters = &matrix_parameters,
            .macropore_mualem_van_genuchten_parameters = &macropore_parameters,
            .macropore_spacing_m = &spacing,
            .macropore_radius_m = &radius,
            .dual_domain_exchange_enabled = &enabled,
            .gravitational_potential_megapascal = &zero,
            .osmotic_potential_megapascal = &zero,
            .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity,
            .vertical_thickness_m = &one,
            .osmotic_potential_multiplier = 1,
        },
        &.{},
        &.{},
        .{ .max_iterations = 40 },
    );
    try std.testing.expect(result.iterations < 40);
    try std.testing.expect(grid.matrix_liquid_water_m3[0] > 0.10);
    try std.testing.expect(grid.macropore_liquid_water_m3[0] < 0.08);
    try std.testing.expectApproxEqAbs(
        before,
        grid.matrix_liquid_water_m3[0] + grid.macropore_liquid_water_m3[0],
        1e-12,
    );
}

test "large geospatial dual-domain volume converges within runtime NPH" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 1,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    const matrix_bulk_volume_m3 = [_]f64{1.0e7};
    grid.matrix_pore_capacity_m3[0] = 4.5e6;
    grid.matrix_liquid_water_m3[0] = 1.0e6;
    grid.macropore_pore_capacity_m3[0] = 1.0e6;
    grid.macropore_liquid_water_m3[0] = 8.0e5;
    const water_before_m3 = grid.matrix_liquid_water_m3[0] +
        grid.macropore_liquid_water_m3[0];
    const matrix_parameters =
        [_]retention.MualemVanGenuchtenParameters{.{
            .residual_water_content_m3_per_m3 = 0.05,
            .saturated_water_content_m3_per_m3 = 0.45,
            .alpha_per_m = 3.6,
            .n = 1.56,
            .saturated_hydraulic_conductivity_m_per_h = 0.01,
        }};
    const macropore_parameters =
        [_]retention.MualemVanGenuchtenParameters{.{
            .residual_water_content_m3_per_m3 = 0,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 15,
            .n = 2.68,
            .saturated_hydraulic_conductivity_m_per_h = 0.1,
        }};
    const curves = [_]retention.ResolvedCurve{testCurve()};
    const zero = [_]f64{0};
    const one = [_]f64{1};
    const spacing_m = [_]f64{0.2};
    const radius_m = [_]f64{0.001};
    const exchange_enabled = [_]bool{true};
    const conductivity = [_]f64{ 0.002, 0.002, 0.002 };
    const result = try solve(
        std.testing.allocator,
        &grid,
        &.{},
        .{
            .matrix_bulk_volume_m3 = &matrix_bulk_volume_m3,
            .retention_curve = &curves,
            .mualem_van_genuchten_parameters = &matrix_parameters,
            .macropore_mualem_van_genuchten_parameters = &macropore_parameters,
            .macropore_spacing_m = &spacing_m,
            .macropore_radius_m = &radius_m,
            .dual_domain_exchange_enabled = &exchange_enabled,
            .gravitational_potential_megapascal = &zero,
            .osmotic_potential_megapascal = &zero,
            .matrix_hydraulic_conductivity_m2_per_h_megapascal = &conductivity,
            .vertical_thickness_m = &one,
            .osmotic_potential_multiplier = 1,
        },
        &.{},
        &.{},
        .{
            .max_iterations = 20,
            .absolute_tolerance_m3 = 1e-11,
            .relative_tolerance = 1e-8,
        },
    );
    try std.testing.expect(result.iterations <= 20);
    try std.testing.expect(result.maximum_scaled_residual <= 1);
    try std.testing.expectApproxEqAbs(
        water_before_m3,
        grid.matrix_liquid_water_m3[0] +
            grid.macropore_liquid_water_m3[0],
        1e-11 + water_before_m3 * 1e-8,
    );
}
