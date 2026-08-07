const std = @import("std");
const builtin = @import("builtin");
const numerics = @import("../../core/numerics.zig");
const grid_module = @import("../../state/grid.zig");
const water_flux = @import("../water/flux.zig");
const transport_hydrology = @import("../../transport/hydrology.zig");

pub const FaceGeometry = struct { source_path_length_m: []const f64, destination_path_length_m: []const f64, face_area_m2: []const f64 };

pub const Properties = struct {
    vapor_diffusivity_m2_per_h: []const f64,
    air_fraction: []const f64,
    porosity_fraction: []const f64,
    tortuosity: f64,
};

pub const Options = struct {
    max_iterations: u16,
    absolute_tolerance_m3: f64 = 1e-15,
    relative_tolerance: f64 = 1e-9,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
};

pub const Result = struct { iterations: u16, newton_raphson_steps: u16, picard_steps: u16, maximum_scaled_residual: f64 };

/// NPH-bounded replacement for WATSUB's repeated VOLV diffusion sweeps. The
/// converged face flux is published atomically to the shared water/heat graph.
pub fn solveAndBindTransportFaces(allocator: std.mem.Allocator, grid: *grid_module.GridState, hydrology: *transport_hydrology.State, shared_faces: *transport_hydrology.SoilFaces, geometry: FaceGeometry, properties: Properties, options: Options) !Result {
    const count = shared_faces.micropore_faces.len;
    const cells = grid.layer_count;
    if (shared_faces.vapor_flux_m3_per_step.len != count or geometry.source_path_length_m.len != count or geometry.destination_path_length_m.len != count or geometry.face_area_m2.len != count or properties.vapor_diffusivity_m2_per_h.len != cells or properties.air_fraction.len != cells or properties.porosity_fraction.len != cells or hydrology.water_vapor_volume_m3.len != cells) return error.SoilVaporSolverDimensionMismatch;
    try validateOptions(properties, options);
    const base = try allocator.dupe(f64, grid.water_vapor_volume_m3);
    defer allocator.free(base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, cells);
    defer allocator.free(residual);
    const target = try allocator.alloc(f64, cells);
    defer allocator.free(target);
    const scratch = try allocator.alloc(f64, cells);
    defer allocator.free(scratch);
    const trial_flux = try allocator.alloc(f64, count);
    defer allocator.free(trial_flux);
    const probe = try allocator.alloc(f64, cells);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, cells);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, cells);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, cells);
    defer allocator.free(candidate_residual);
    const jacobian = try allocator.alloc(f64, try std.math.mul(usize, cells, cells));
    defer allocator.free(jacobian);
    const newton_delta = try allocator.alloc(f64, cells);
    defer allocator.free(newton_delta);
    var conserved_vapor_total_m3: f64 = 0;
    for (base) |value| conserved_vapor_total_m3 += value;
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(grid, shared_faces, geometry, properties, base, current, target, residual, scratch, trial_flux);
        const norm = try scaledNorm(current, residual, options);
        if (norm <= 1) {
            try residualAt(grid, shared_faces, geometry, properties, base, current, target, residual, scratch, shared_faces.vapor_flux_m3_per_step);
            @memcpy(grid.water_vapor_volume_m3, current);
            @memcpy(hydrology.water_vapor_volume_m3, current);
            @memset(hydrology.vapor_face_flux_m3_per_step, 0);
            for (shared_faces.micropore_faces, shared_faces.direction_axis, shared_faces.vapor_flux_m3_per_step) |face, axis, flux| hydrology.vapor_face_flux_m3_per_step[face.first_cell * 3 + axis] = flux;
            try grid.validateFinite();
            try hydrology.validateFinite();
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm };
        }
        var accepted_newton = false;
        var dense_jacobian_valid = true;
        for (0..cells) |column| {
            @memcpy(probe, current);
            const perturbation = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0e-3, @abs(current[column]));
            probe[column] += perturbation;
            if (residualAt(grid, shared_faces, geometry, properties, base, probe, target, probe_residual, scratch, trial_flux)) |_| {
                for (0..cells) |row| jacobian[row * cells + column] = probe_residual[row];
                if (current[column] > perturbation) {
                    @memcpy(probe, current);
                    probe[column] -= perturbation;
                    if (residualAt(grid, shared_faces, geometry, properties, base, probe, target, candidate_residual, scratch, trial_flux)) |_| {
                        for (0..cells) |row| jacobian[row * cells + column] = (jacobian[row * cells + column] - candidate_residual[row]) / (2.0 * perturbation);
                    } else |_| {
                        dense_jacobian_valid = false;
                        break;
                    }
                } else {
                    for (0..cells) |row| jacobian[row * cells + column] = (jacobian[row * cells + column] - residual[row]) / perturbation;
                }
            } else |_| {
                dense_jacobian_valid = false;
                break;
            }
        }
        if (dense_jacobian_valid) {
            for (residual, newton_delta) |value, *rhs| rhs.* = -value;
            if (numerics.solveDenseLinearSystem(jacobian, newton_delta, cells)) {
                var line_fraction: f64 = 1;
                var line_search: u8 = 0;
                while (line_search < 8) : (line_search += 1) {
                    var valid = true;
                    for (current, newton_delta, candidate) |value, delta, *next| {
                        next.* = value + line_fraction * delta;
                        if (!std.math.isFinite(next.*) or next.* < -1e-12) valid = false else next.* = @max(0, next.*);
                    }
                    if (valid) valid = projectConservedTotal(candidate, conserved_vapor_total_m3);
                    if (valid) {
                        if (residualAt(grid, shared_faces, geometry, properties, base, candidate, target, candidate_residual, scratch, trial_flux)) |_| {
                            if (try scaledNorm(candidate, candidate_residual, options) < norm) {
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
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(grid, shared_faces, geometry, properties, base, probe, target, probe_residual, scratch, trial_flux)) |_| {
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
                        if (residualAt(grid, shared_faces, geometry, properties, base, candidate, target, candidate_residual, scratch, trial_flux)) |_| {
                            if (try scaledNorm(candidate, candidate_residual, options) < norm) {
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
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    try residualAt(grid, shared_faces, geometry, properties, base, current, target, residual, scratch, trial_flux);
    const final_norm = try scaledNorm(current, residual, options);
    if (final_norm <= 1) {
        try residualAt(grid, shared_faces, geometry, properties, base, current, target, residual, scratch, shared_faces.vapor_flux_m3_per_step);
        @memcpy(grid.water_vapor_volume_m3, current);
        @memcpy(hydrology.water_vapor_volume_m3, current);
        @memset(hydrology.vapor_face_flux_m3_per_step, 0);
        for (shared_faces.micropore_faces, shared_faces.direction_axis, shared_faces.vapor_flux_m3_per_step) |face, axis, flux| hydrology.vapor_face_flux_m3_per_step[face.first_cell * 3 + axis] = flux;
        try grid.validateFinite();
        try hydrology.validateFinite();
        return .{ .iterations = options.max_iterations, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = final_norm };
    }
    if (!builtin.is_test) std.log.err("soil vapor Newton-Raphson/Picard did not converge: iterations={d} newton_raphson_steps={d} picard_steps={d} maximum_scaled_residual={e}", .{ options.max_iterations, newton_steps, picard_steps, final_norm });
    return error.SoilVaporSolverDidNotConverge;
}

fn projectConservedTotal(values: []f64, required_total: f64) bool {
    if (values.len == 0) return false;
    var actual_total: f64 = 0;
    var largest_index: usize = 0;
    for (values, 0..) |value, index| {
        if (!std.math.isFinite(value) or value < 0) return false;
        actual_total += value;
        if (value > values[largest_index]) largest_index = index;
    }
    values[largest_index] += required_total - actual_total;
    return std.math.isFinite(values[largest_index]) and values[largest_index] >= 0;
}

fn residualAt(grid: *const grid_module.GridState, faces: *const transport_hydrology.SoilFaces, geometry: FaceGeometry, properties: Properties, base: []const f64, trial: []const f64, target: []f64, residual: []f64, scratch: []f64, output_flux: []f64) !void {
    @memcpy(scratch, trial);
    @memcpy(target, base);
    @memset(output_flux, 0);
    for (faces.micropore_faces, 0..) |face, index| {
        if (face.first_cell >= grid.layer_count or face.second_cell >= grid.layer_count or face.first_cell == face.second_cell) return error.InvalidSoilVaporFace;
        const source = face.first_cell;
        const destination = face.second_cell;
        const flux = try water_flux.calculateVaporFaceFlux(.{ .source_air_volume_m3 = grid.air_volume_m3[source], .destination_air_volume_m3 = grid.air_volume_m3[destination], .source_vapor_volume_m3 = scratch[source], .destination_vapor_volume_m3 = scratch[destination], .source_vapor_diffusivity_m2_per_h = properties.vapor_diffusivity_m2_per_h[source], .destination_vapor_diffusivity_m2_per_h = properties.vapor_diffusivity_m2_per_h[destination], .source_air_fraction = properties.air_fraction[source], .destination_air_fraction = properties.air_fraction[destination], .source_porosity_fraction = properties.porosity_fraction[source], .destination_porosity_fraction = properties.porosity_fraction[destination], .tortuosity = properties.tortuosity, .source_path_length_m = geometry.source_path_length_m[index], .destination_path_length_m = geometry.destination_path_length_m[index], .face_area_m2 = geometry.face_area_m2[index], .time_fraction = 1 });
        // In the implicit whole-hour solve, positivity is enforced on the
        // converged inventories rather than by the explicit WATSUB donor cap.
        // Using the unlimited constitutive flux keeps the residual continuous;
        // the resulting diffusion matrix is conservative and positivity
        // preserving, without repeated sub-hour sweeps.
        const implicit_flux_m3 = flux.unlimited_vapor_m3;
        output_flux[index] = implicit_flux_m3;
        target[source] -= implicit_flux_m3;
        target[destination] += implicit_flux_m3;
    }
    for (target, trial, residual) |value, trial_value, *difference| {
        if (!std.math.isFinite(value)) return error.InvalidSoilVaporCandidate;
        difference.* = value - trial_value;
    }
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidSoilVaporCandidate;
        candidate.* = @max(0.0, candidate.*);
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteSoilVaporSolverState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_m3 + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn validateOptions(properties: Properties, options: Options) !void {
    if (!std.math.isFinite(properties.tortuosity) or properties.tortuosity < 0 or options.max_iterations == 0 or !std.math.isFinite(options.absolute_tolerance_m3) or options.absolute_tolerance_m3 <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction) return error.InvalidSoilVaporSolverOptions;
}

test "vapor hybrid solve conserves VOLV and publishes shared face flux" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.air_volume_m3, 1);
    grid.water_vapor_volume_m3[0] = 0.01;
    var snow = try @import("../solute/snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try transport_hydrology.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    try hydrology.syncStorage(&grid, &snow);
    var faces = try transport_hydrology.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    const one = [_]f64{1};
    const cell_one = [_]f64{ 1, 1 };
    const porosity = [_]f64{ 0.5, 0.5 };
    const result = try solveAndBindTransportFaces(std.testing.allocator, &grid, &hydrology, &faces, .{ .source_path_length_m = &one, .destination_path_length_m = &one, .face_area_m2 = &one }, .{ .vapor_diffusivity_m2_per_h = &cell_one, .air_fraction = &porosity, .porosity_fraction = &porosity, .tortuosity = 1 }, .{ .max_iterations = 20 });
    try std.testing.expect(result.iterations < 20);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), grid.water_vapor_volume_m3[0] + grid.water_vapor_volume_m3[1], 1e-12);
    try std.testing.expect(faces.vapor_flux_m3_per_step[0] > 0);
    try std.testing.expectEqual(faces.vapor_flux_m3_per_step[0], hydrology.vapor_face_flux_m3_per_step[0]);
}

test "dense vapor Newton converges a linear face within one NPH iteration" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.air_volume_m3, 1);
    grid.water_vapor_volume_m3[0] = 0.01;
    var snow = try @import("../solute/snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try transport_hydrology.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    try hydrology.syncStorage(&grid, &snow);
    var faces = try transport_hydrology.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    const one = [_]f64{1};
    const cell_one = [_]f64{ 1, 1 };
    const porosity = [_]f64{ 0.5, 0.5 };
    const result = try solveAndBindTransportFaces(std.testing.allocator, &grid, &hydrology, &faces, .{ .source_path_length_m = &one, .destination_path_length_m = &one, .face_area_m2 = &one }, .{ .vapor_diffusivity_m2_per_h = &cell_one, .air_fraction = &porosity, .porosity_fraction = &porosity, .tortuosity = 1 }, .{ .max_iterations = 1, .minimum_newton_fraction = 0.05, .maximum_newton_fraction = 0.05 });
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), grid.water_vapor_volume_m3[0] + grid.water_vapor_volume_m3[1], 1e-12);
    try std.testing.expect(faces.vapor_flux_m3_per_step[0] > 0);
}
