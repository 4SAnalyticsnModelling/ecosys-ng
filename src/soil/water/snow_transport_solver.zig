const std = @import("std");
const snow = @import("../solute/snow_solute_transport.zig");

pub const Inputs = struct {
    atmospheric_top_input_g: []const f64,
    transport_water_volume_m3: []const f64,
    water_flux_to_lower_m3: []const f64,
    litter_water_flux_m3: []const f64,
    soil_micropore_water_flux_m3: []const f64,
    soil_macropore_water_flux_m3: []const f64,
    surface_partitions: []const snow.SurfacePartition,
};

pub const Options = struct {
    absolute_tolerance_g: f64 = 1e-12,
    relative_tolerance: f64 = 1e-8,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 0.5,
    minimum_newton_fraction: f64 = 0.05,
    maximum_newton_fraction: f64 = 1.5,
    /// Source-derived snowpack ceiling (20 unless supplied at runtime).
    max_iterations: u16 = 20,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
};

/// Implicitly advances only snow solute transport. The rest of the ecosystem
/// model is not repeated while the snow state converges.
pub fn solve(allocator: std.mem.Allocator, state: *snow.State, inputs: Inputs, options: Options, output_surface_discharge: []snow.SurfaceDischarge) !Result {
    try validate(state, inputs, options, output_surface_discharge);
    const base = try allocator.dupe(f64, state.amount_g);
    defer allocator.free(base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, base.len);
    defer allocator.free(residual);
    const probe = try allocator.alloc(f64, base.len);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, base.len);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, base.len);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, base.len);
    defer allocator.free(candidate_residual);
    const target = try allocator.alloc(f64, base.len);
    defer allocator.free(target);
    var scratch = try snow.State.init(allocator, state.cell_count, state.layer_capacity);
    defer scratch.deinit();
    @memcpy(scratch.active, state.active);
    @memcpy(scratch.liquid_water_volume_m3, inputs.transport_water_volume_m3);

    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(allocator, &scratch, base, current, inputs, target, residual, null);
        const norm = try scaledNorm(current, residual, options);
        if (norm <= 1) {
            @memcpy(state.amount_g, current);
            try residualAt(allocator, &scratch, base, current, inputs, target, residual, output_surface_discharge);
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm };
        }
        var accepted_newton = false;
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(allocator, &scratch, base, probe, inputs, target, probe_residual, null)) |_| {
                var numerator: f64 = 0;
                var denominator: f64 = 0;
                for (residual, probe_residual) |value, sampled| {
                    const derivative = (sampled - value) / options.directional_probe_fraction;
                    numerator += value * derivative;
                    denominator += derivative * derivative;
                }
                if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
                    const fraction = std.math.clamp(-numerator / denominator, options.minimum_newton_fraction, options.maximum_newton_fraction);
                    if (addDirection(current, residual, fraction, candidate)) |_| {
                        if (residualAt(allocator, &scratch, base, candidate, inputs, target, candidate_residual, null)) |_| {
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
        if (maximumDifference(current, candidate) <= std.math.floatEps(f64) * @max(1.0, maximumMagnitude(current))) return error.SnowTransportSolverStagnated;
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    return error.SnowTransportSolverDidNotConverge;
}

fn residualAt(allocator: std.mem.Allocator, scratch: *snow.State, base: []const f64, trial: []const f64, inputs: Inputs, target: []f64, residual: []f64, output_discharge: ?[]snow.SurfaceDischarge) !void {
    @memcpy(scratch.amount_g, trial);
    var fluxes = try snow.calculateFluxes(allocator, scratch, inputs.water_flux_to_lower_m3, inputs.litter_water_flux_m3, inputs.soil_micropore_water_flux_m3, inputs.soil_macropore_water_flux_m3, inputs.surface_partitions);
    defer fluxes.deinit();
    @memcpy(target, base);
    for (0..scratch.cell_count) |cell| {
        const top = try scratch.layerIndex(cell, 0);
        for (0..snow.species_count) |species| target[top * snow.species_count + species] += inputs.atmospheric_top_input_g[cell * snow.species_count + species];
        var discharged = false;
        for (0..scratch.layer_capacity) |layer| {
            const index = try scratch.layerIndex(cell, layer);
            if (!scratch.active[index]) continue;
            const has_lower = layer + 1 < scratch.layer_capacity and scratch.active[try scratch.layerIndex(cell, layer + 1)];
            for (0..snow.species_count) |species| {
                const component = index * snow.species_count + species;
                const downward = fluxes.downward_g[component];
                const actual_downward = @min(target[component], downward);
                target[component] -= actual_downward;
                if (has_lower) target[(index + 1) * snow.species_count + species] += actual_downward;
                if (!has_lower and !discharged) {
                    const requested = fluxes.surface_discharge[cell].litter_g[species] + fluxes.surface_discharge[cell].soil_nonband_g[species] + fluxes.surface_discharge[cell].soil_band_g[species];
                    target[component] -= @min(target[component], requested);
                }
            }
            if (!has_lower and !discharged) discharged = true;
        }
    }
    if (output_discharge) |output| @memcpy(output, fluxes.surface_discharge);
    for (target, trial, residual) |fixed_point, value, *difference| {
        if (!std.math.isFinite(fixed_point) or fixed_point < -1e-12) return error.InvalidImplicitSnowCandidate;
        difference.* = fixed_point - value;
    }
}

fn validate(state: *const snow.State, inputs: Inputs, options: Options, output: []snow.SurfaceDischarge) !void {
    const layers = try std.math.mul(usize, state.cell_count, state.layer_capacity);
    if (inputs.atmospheric_top_input_g.len != state.cell_count * snow.species_count or inputs.transport_water_volume_m3.len != layers or inputs.water_flux_to_lower_m3.len != layers or inputs.litter_water_flux_m3.len != state.cell_count or inputs.soil_micropore_water_flux_m3.len != state.cell_count or inputs.soil_macropore_water_flux_m3.len != state.cell_count or inputs.surface_partitions.len != state.cell_count or output.len != state.cell_count) return error.SnowTransportInputSizeMismatch;
    if (!std.math.isFinite(options.absolute_tolerance_g) or options.absolute_tolerance_g <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction or options.max_iterations == 0) return error.InvalidSnowTransportSolverOptions;
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidImplicitSnowCandidate;
        candidate.* = @max(0, candidate.*);
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteImplicitSnowState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_g + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
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

fn testPartition() snow.SurfacePartition {
    return .{ .litter_cover_fraction = 0.25, .bare_soil_fraction = 0.75, .nonband_ammonium_fraction = 0.6, .band_ammonium_fraction = 0.4, .nonband_nitrate_fraction = 0.7, .band_nitrate_fraction = 0.3, .nonband_phosphate_fraction = 0.8, .band_phosphate_fraction = 0.2 };
}

test "snow hybrid converges before legacy 20 iteration ceiling" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.active[0] = true;
    state.active[1] = true;
    state.liquid_water_volume_m3[0] = 2;
    state.liquid_water_volume_m3[1] = 2;
    @memset(try state.amounts(0, 0), 10);
    var discharge: [1]snow.SurfaceDischarge = undefined;
    const zero_input = [_]f64{0} ** snow.species_count;
    const result = try solve(std.testing.allocator, &state, .{ .atmospheric_top_input_g = &zero_input, .transport_water_volume_m3 = state.liquid_water_volume_m3, .water_flux_to_lower_m3 = &[_]f64{ 0, 0.2 }, .litter_water_flux_m3 = &[_]f64{0.1}, .soil_micropore_water_flux_m3 = &[_]f64{0.1}, .soil_macropore_water_flux_m3 = &[_]f64{0}, .surface_partitions = &[_]snow.SurfacePartition{testPartition()} }, .{}, &discharge);
    try std.testing.expect(result.iterations < 20);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
}

test "failed snow hybrid leaves caller state unchanged" {
    var state = try snow.State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    state.active[0] = true;
    state.liquid_water_volume_m3[0] = 1;
    state.amount_g[0] = 1;
    var discharge: [1]snow.SurfaceDischarge = undefined;
    const zero_input = [_]f64{0} ** snow.species_count;
    try std.testing.expectError(error.SnowTransportSolverDidNotConverge, solve(std.testing.allocator, &state, .{ .atmospheric_top_input_g = &zero_input, .transport_water_volume_m3 = state.liquid_water_volume_m3, .water_flux_to_lower_m3 = &[_]f64{0}, .litter_water_flux_m3 = &[_]f64{0.2}, .soil_micropore_water_flux_m3 = &[_]f64{0}, .soil_macropore_water_flux_m3 = &[_]f64{0}, .surface_partitions = &[_]snow.SurfacePartition{testPartition()} }, .{ .max_iterations = 1 }, &discharge));
    try std.testing.expectEqual(@as(f64, 1), state.amount_g[0]);
}
