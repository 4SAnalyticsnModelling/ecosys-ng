const std = @import("std");
const grid_module = @import("grid.zig");
const hydrology_module = @import("transport_hydrology.zig");
const organic = @import("soil_organic_initialization.zig");
const face_parameters = @import("soil_organic_face_parameters.zig");
const numerics = @import("numerics.zig");

pub const component_count: usize = face_parameters.component_count;
pub const components_per_substrate: usize = face_parameters.components_per_substrate;
pub const Component = face_parameters.Component;

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    micropore_amount_g: []f64,
    macropore_amount_g: []f64,
    boundary_net_flux_g: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.ZeroSoilOrganicTransportLayers;
        const count = try std.math.mul(usize, layer_count, component_count);
        const micropore = try allocateZero(allocator, count);
        errdefer allocator.free(micropore);
        const macropore = try allocateZero(allocator, count);
        errdefer allocator.free(macropore);
        const boundary = try allocateZero(allocator, count);
        return .{ .allocator = allocator, .layer_count = layer_count, .micropore_amount_g = micropore, .macropore_amount_g = macropore, .boundary_net_flux_g = boundary };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.boundary_net_flux_g);
        self.allocator.free(self.macropore_amount_g);
        self.allocator.free(self.micropore_amount_g);
        self.* = undefined;
    }

    /// STARTS initializes dissolved organic matter in the matrix domain.
    pub fn initializeFromProfile(self: *State, profile: *const organic.State) !void {
        try validateProfileDimensions(self, profile);
        @memset(self.macropore_amount_g, 0);
        @memset(self.boundary_net_flux_g, 0);
        try self.exportMicroporeFromProfile(profile);
    }

    /// Exports the process-owned matrix pools immediately before TRNSFR.
    pub fn exportMicroporeFromProfile(self: *State, profile: *const organic.State) !void {
        try validateProfileDimensions(self, profile);
        for (0..self.layer_count) |layer| for (0..organic.substrate_count) |substrate| {
            const pool = profile.dissolved[layer * organic.substrate_count + substrate];
            const acetate = profile.dissolved_acetate_carbon_g_c[layer * organic.substrate_count + substrate];
            try validatePool(pool, acetate);
            const base = substrateBase(layer, substrate);
            self.micropore_amount_g[base + @intFromEnum(Component.dissolved_organic_carbon)] = pool.carbon_g_c;
            self.micropore_amount_g[base + @intFromEnum(Component.dissolved_organic_nitrogen)] = pool.nitrogen_g_n;
            self.micropore_amount_g[base + @intFromEnum(Component.dissolved_organic_phosphorus)] = pool.phosphorus_g_p;
            self.micropore_amount_g[base + @intFromEnum(Component.dissolved_acetate_carbon)] = acetate;
        };
    }

    /// Imports accepted TRNSFR matrix inventories for subsequent NITRO use.
    pub fn importMicroporeIntoProfile(self: *const State, profile: *organic.State) !void {
        try validateProfileDimensions(self, profile);
        for (0..self.layer_count) |layer| for (0..organic.substrate_count) |substrate| {
            const base = substrateBase(layer, substrate);
            profile.dissolved[layer * organic.substrate_count + substrate] = .{
                .carbon_g_c = self.micropore_amount_g[base + @intFromEnum(Component.dissolved_organic_carbon)],
                .nitrogen_g_n = self.micropore_amount_g[base + @intFromEnum(Component.dissolved_organic_nitrogen)],
                .phosphorus_g_p = self.micropore_amount_g[base + @intFromEnum(Component.dissolved_organic_phosphorus)],
            };
            profile.dissolved_acetate_carbon_g_c[layer * organic.substrate_count + substrate] = self.micropore_amount_g[base + @intFromEnum(Component.dissolved_acetate_carbon)];
        };
    }
};

pub const Options = struct {
    absolute_tolerance_g: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    max_iterations: u16,
    maximum_convective_fraction: f64 = 1,
    pore_exchange_fraction: f64 = 1,
};

pub const Inputs = struct {
    micropore_conductance_m3_per_step: []const f64,
    macropore_conductance_m3_per_step: []const f64,
    matrix_water_m3: []const f64,
    macropore_water_m3: []const f64,
    layer_bulk_volume_m3: []const f64,
    micropore_external_water_flux_m3_per_step: []const f64,
    macropore_external_water_flux_m3_per_step: []const f64,
    recharge_concentration_g_per_m3: []const f64,
};

pub const Result = struct {
    micropore_iterations: u16,
    macropore_iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

/// Atomic TRNSFR organic-solute step. The profile, both pore inventories, and
/// boundary ledger are published only after every solve and sufficiency check.
pub fn advance(
    allocator: std.mem.Allocator,
    state: *State,
    profile: *organic.State,
    faces: *const hydrology_module.SoilFaces,
    inputs: Inputs,
    options: Options,
) !Result {
    try validateAdvance(state, profile, faces, inputs, options);
    const profile_dissolved_before = try allocator.dupe(organic.ElementPool, profile.dissolved);
    defer allocator.free(profile_dissolved_before);
    const profile_acetate_before = try allocator.dupe(f64, profile.dissolved_acetate_carbon_g_c);
    defer allocator.free(profile_acetate_before);
    const micropore_before = try allocator.dupe(f64, state.micropore_amount_g);
    defer allocator.free(micropore_before);
    const macropore_before = try allocator.dupe(f64, state.macropore_amount_g);
    defer allocator.free(macropore_before);
    const boundary_before = try allocator.dupe(f64, state.boundary_net_flux_g);
    defer allocator.free(boundary_before);
    var committed = false;
    defer if (!committed) {
        @memcpy(profile.dissolved, profile_dissolved_before);
        @memcpy(profile.dissolved_acetate_carbon_g_c, profile_acetate_before);
        @memcpy(state.micropore_amount_g, micropore_before);
        @memcpy(state.macropore_amount_g, macropore_before);
        @memcpy(state.boundary_net_flux_g, boundary_before);
    };

    try state.exportMicroporeFromProfile(profile);
    const micropore_result = try solveDomain(allocator, state.micropore_amount_g, inputs.matrix_water_m3, faces.micropore_faces, inputs.micropore_conductance_m3_per_step, options);
    const macropore_result = try solveDomain(allocator, state.macropore_amount_g, inputs.macropore_water_m3, faces.macropore_faces, inputs.macropore_conductance_m3_per_step, options);
    @memset(state.boundary_net_flux_g, 0);
    for (0..state.layer_count) |layer| {
        try commitBoundary(state.micropore_amount_g[layer * component_count ..][0..component_count], inputs.matrix_water_m3[layer], inputs.micropore_external_water_flux_m3_per_step[layer], inputs.recharge_concentration_g_per_m3[layer * component_count ..][0..component_count], options.maximum_convective_fraction, state.boundary_net_flux_g[layer * component_count ..][0..component_count]);
        try commitBoundary(state.macropore_amount_g[layer * component_count ..][0..component_count], inputs.macropore_water_m3[layer], inputs.macropore_external_water_flux_m3_per_step[layer], inputs.recharge_concentration_g_per_m3[layer * component_count ..][0..component_count], options.maximum_convective_fraction, state.boundary_net_flux_g[layer * component_count ..][0..component_count]);
        for (0..component_count) |component| {
            const index = layer * component_count + component;
            const exchange_g = try poreExchange(
                state.micropore_amount_g[index],
                state.macropore_amount_g[index],
                inputs.matrix_water_m3[layer],
                inputs.macropore_water_m3[layer],
                inputs.layer_bulk_volume_m3[layer],
                options.pore_exchange_fraction,
            );
            state.micropore_amount_g[index] += exchange_g;
            state.macropore_amount_g[index] -= exchange_g;
        }
    }
    try state.importMicroporeIntoProfile(profile);
    committed = true;
    return .{
        .micropore_iterations = micropore_result.iterations,
        .macropore_iterations = macropore_result.iterations,
        .newton_raphson_steps = micropore_result.newton_steps + macropore_result.newton_steps,
        .picard_steps = micropore_result.picard_steps + macropore_result.picard_steps,
    };
}

const SolverResult = struct { iterations: u16, newton_steps: u16, picard_steps: u16 };

fn solveDomain(allocator: std.mem.Allocator, amounts_g: []f64, water_m3: []const f64, faces: []const @import("solute_transport.zig").Face, conductance: []const f64, options: Options) !SolverResult {
    const base = try allocator.dupe(f64, amounts_g);
    defer allocator.free(base);
    const current = try allocator.dupe(f64, base);
    defer allocator.free(current);
    const residual = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(residual);
    const probe = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(probe);
    const probe_residual = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(probe_residual);
    const candidate = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(candidate);
    const candidate_residual = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(candidate_residual);
    const dense_direction = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(dense_direction);
    const layer_count = amounts_g.len / component_count;
    const dense_matrix = try allocator.alloc(
        f64,
        try std.math.mul(usize, layer_count, layer_count),
    );
    defer allocator.free(dense_matrix);
    const dense_rhs = try allocator.alloc(f64, layer_count);
    defer allocator.free(dense_rhs);
    const fixed_point = try allocator.alloc(f64, amounts_g.len);
    defer allocator.free(fixed_point);
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(base, current, water_m3, faces, conductance, options.maximum_convective_fraction, fixed_point, residual);
        const norm = try scaledNorm(current, residual, options);
        if (norm <= 1) {
            try enforceInternalConservation(base, current);
            @memcpy(amounts_g, current);
            return .{ .iterations = iteration + 1, .newton_steps = newton_steps, .picard_steps = picard_steps };
        }
        var accepted_newton = false;
        const limiting_index =
            try worstResidualIndex(current, residual, options);
        const limiting_component = limiting_index % component_count;
        if (try denseComponentNewtonDirection(
            base,
            current,
            residual,
            water_m3,
            faces,
            conductance,
            options,
            limiting_component,
            fixed_point,
            probe,
            probe_residual,
            dense_matrix,
            dense_rhs,
            dense_direction,
        )) {
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
                        base,
                        candidate,
                        water_m3,
                        faces,
                        conductance,
                        options.maximum_convective_fraction,
                        fixed_point,
                        candidate_residual,
                    )) |_| {
                        if (try scaledNorm(
                            candidate,
                            candidate_residual,
                            options,
                        ) < norm) {
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
        try addDirection(current, residual, 0.5, probe);
        try residualAt(base, probe, water_m3, faces, conductance, options.maximum_convective_fraction, fixed_point, probe_residual);
        var numerator: f64 = 0;
        var denominator: f64 = 0;
        for (residual, probe_residual) |value, sampled| {
            const derivative = (sampled - value) / 0.5;
            numerator += value * derivative;
            denominator += derivative * derivative;
        }
        if (std.math.isFinite(denominator) and denominator > std.math.floatEps(f64)) {
            const fraction = std.math.clamp(-numerator / denominator, 0.05, 1.5);
            if (addDirection(current, residual, fraction, candidate)) |_| {
                if (residualAt(base, candidate, water_m3, faces, conductance, options.maximum_convective_fraction, fixed_point, candidate_residual)) |_| {
                    if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                        @memcpy(current, candidate);
                        newton_steps += 1;
                        accepted_newton = true;
                    }
                } else |_| {}
            } else |_| {}
        }
        if (accepted_newton) continue;
        try addDirection(current, residual, options.picard_relaxation, candidate);
        @memcpy(current, candidate);
        picard_steps += 1;
    }
    try residualAt(
        base,
        current,
        water_m3,
        faces,
        conductance,
        options.maximum_convective_fraction,
        fixed_point,
        residual,
    );
    const final_norm = try scaledNorm(current, residual, options);
    if (final_norm <= 1) {
        try enforceInternalConservation(base, current);
        @memcpy(amounts_g, current);
        return .{
            .iterations = options.max_iterations,
            .newton_steps = newton_steps,
            .picard_steps = picard_steps,
        };
    }
    const limiting_index =
        try worstResidualIndex(current, residual, options);
    const limiting_layer = limiting_index / component_count;
    const limiting_component = limiting_index % component_count;
    std.log.err(
        "soil organic transport Newton-Picard exhausted runtime ceiling: max_iterations={d} scaled_residual={e} layer={d} substrate={d} component={d} amount_g={e} residual_g={e} newton_steps={d} picard_steps={d}",
        .{
            options.max_iterations,
            final_norm,
            limiting_layer,
            limiting_component / components_per_substrate,
            limiting_component % components_per_substrate,
            current[limiting_index],
            residual[limiting_index],
            newton_steps,
            picard_steps,
        },
    );
    return error.SoilOrganicTransportDidNotConverge;
}

/// Internal TRNSFR faces only redistribute each organic component. Newton or
/// Picard convergence controls the local equations, but a tolerance-sized sum
/// of residuals must not become a landscape source. Correct the roundoff-sized
/// closure remainder in the largest receiving pool before accepting the solve.
fn enforceInternalConservation(base: []const f64, current: []f64) !void {
    if (base.len == 0 or base.len != current.len or
        base.len % component_count != 0)
        return error.SoilOrganicTransportDimensionMismatch;
    const layer_count = base.len / component_count;
    for (0..component_count) |component| {
        var base_total: f64 = 0;
        var current_total: f64 = 0;
        var largest_index = component;
        for (0..layer_count) |layer| {
            const index = layer * component_count + component;
            const before = base[index];
            const after = current[index];
            if (!std.math.isFinite(before) or before < 0 or
                !std.math.isFinite(after) or after < 0)
                return error.NonFiniteSoilOrganicTransport;
            base_total += before;
            current_total += after;
            if (after > current[largest_index]) largest_index = index;
        }
        const correction = base_total - current_total;
        const corrected = current[largest_index] + correction;
        if (!std.math.isFinite(corrected) or corrected < 0)
            return error.SoilOrganicTransportConservationFailure;
        current[largest_index] = corrected;
    }
}

fn denseComponentNewtonDirection(
    base: []const f64,
    current: []const f64,
    residual: []const f64,
    water_m3: []const f64,
    faces: []const @import("solute_transport.zig").Face,
    conductance: []const f64,
    options: Options,
    component: usize,
    fixed_point: []f64,
    sampled_state: []f64,
    sampled_residual: []f64,
    matrix: []f64,
    rhs: []f64,
    direction: []f64,
) !bool {
    const layer_count = current.len / component_count;
    if (component >= component_count or
        current.len % component_count != 0 or
        matrix.len != layer_count * layer_count or
        rhs.len != layer_count or
        direction.len != current.len)
        return error.SoilOrganicTransportDimensionMismatch;
    @memset(direction, 0);
    for (0..layer_count) |row|
        rhs[row] = -residual[row * component_count + component];
    for (0..layer_count) |column| {
        const state_index = column * component_count + component;
        const epsilon = @max(
            1.0e-6,
            @sqrt(std.math.floatEps(f64)) *
                @max(1.0, @abs(current[state_index])),
        );
        @memcpy(sampled_state, current);
        sampled_state[state_index] += epsilon;
        residualAt(
            base,
            sampled_state,
            water_m3,
            faces,
            conductance,
            options.maximum_convective_fraction,
            fixed_point,
            sampled_residual,
        ) catch return false;
        for (0..layer_count) |row| {
            const residual_index =
                row * component_count + component;
            matrix[row * layer_count + column] =
                (sampled_residual[residual_index] -
                    residual[residual_index]) /
                epsilon;
        }
    }
    if (!numerics.solveDenseLinearSystem(
        matrix,
        rhs,
        layer_count,
    )) return false;
    for (rhs, 0..) |value, layer|
        direction[layer * component_count + component] = value;
    return true;
}

fn residualAt(base: []const f64, trial: []const f64, water: []const f64, faces: []const @import("solute_transport.zig").Face, conductance: []const f64, maximum_fraction: f64, fixed_point: []f64, residual: []f64) !void {
    @memcpy(fixed_point, base);
    for (faces, 0..) |face, face_index| for (0..component_count) |component| {
        const first = face.first_cell * component_count + component;
        const second = face.second_cell * component_count + component;
        const first_concentration = if (water[face.first_cell] > 0) trial[first] / water[face.first_cell] else 0;
        const second_concentration = if (water[face.second_cell] > 0) trial[second] / water[face.second_cell] else 0;
        const donor_fraction = if (face.water_flux_m3_per_step >= 0)
            if (water[face.first_cell] > 0) @min(maximum_fraction, face.water_flux_m3_per_step / water[face.first_cell]) else maximum_fraction
        else if (water[face.second_cell] > 0) @min(maximum_fraction, -face.water_flux_m3_per_step / water[face.second_cell]) else maximum_fraction;
        const convection = if (face.water_flux_m3_per_step >= 0) donor_fraction * trial[first] else -donor_fraction * trial[second];
        const diffusion = conductance[face_index * component_count + component] * (first_concentration - second_concentration);
        const flux = std.math.clamp(convection + diffusion, -fixed_point[second], fixed_point[first]);
        fixed_point[first] -= flux;
        fixed_point[second] += flux;
    };
    for (fixed_point, trial, residual) |target, value, *difference| {
        difference.* = target - value;
        if (!std.math.isFinite(difference.*)) return error.NonFiniteSoilOrganicTransport;
    }
}

fn commitBoundary(amounts: []f64, water_m3: f64, outward_water_m3: f64, recharge_g_per_m3: []const f64, maximum_fraction: f64, ledger: []f64) !void {
    for (amounts, recharge_g_per_m3, ledger) |*amount, recharge, *net| {
        const change = if (outward_water_m3 >= 0)
            -amount.* * (if (water_m3 > 0) @min(maximum_fraction, outward_water_m3 / water_m3) else maximum_fraction)
        else
            -outward_water_m3 * recharge;
        if (!std.math.isFinite(change) or amount.* + change < -1e-12 or !std.math.isFinite(net.* + change)) return error.InvalidSoilOrganicBoundaryFlux;
        amount.* = @max(0, amount.* + change);
        net.* += change;
    }
}

fn poreExchange(micro_g: f64, macro_g: f64, micro_water_m3: f64, macro_water_m3: f64, bulk_volume_m3: f64, fraction: f64) !f64 {
    if (macro_water_m3 == 0) return 0;
    const exchanging_macro_water_m3 = @min(0.05 * bulk_volume_m3, macro_water_m3);
    const combined_water_m3 = micro_water_m3 + exchanging_macro_water_m3;
    if (combined_water_m3 == 0) return 0;
    const exchange = fraction * (macro_g * micro_water_m3 - micro_g * exchanging_macro_water_m3) / combined_water_m3;
    if (!std.math.isFinite(exchange)) return error.NonFiniteSoilOrganicPoreExchange;
    return std.math.clamp(exchange, -micro_g, macro_g);
}

fn componentScale(value: f64, options: Options) f64 {
    // 0.1%-relative loose floor so tiny pools (< 1 g) don't need sub-nanogram convergence.
    const strict = options.absolute_tolerance_g + options.relative_tolerance * @max(1, value);
    const loose = options.relative_tolerance * 1.0e5 * value;
    return @max(strict, loose);
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value < 0 or !std.math.isFinite(difference)) return error.NonFiniteSoilOrganicTransport;
        maximum = @max(maximum, @abs(difference) / componentScale(value, options));
    }
    return maximum;
}

fn worstResidualIndex(
    state: []const f64,
    residual: []const f64,
    options: Options,
) !usize {
    if (state.len == 0 or state.len != residual.len)
        return error.NonFiniteSoilOrganicTransport;
    var limiting_index: usize = 0;
    var limiting_norm: f64 = -1;
    for (state, residual, 0..) |value, difference, index| {
        if (!std.math.isFinite(value) or value < 0 or
            !std.math.isFinite(difference))
            return error.NonFiniteSoilOrganicTransport;
        const norm = @abs(difference) / componentScale(value, options);
        if (norm > limiting_norm) {
            limiting_norm = norm;
            limiting_index = index;
        }
    }
    return limiting_index;
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, change, *candidate| {
        candidate.* = value + fraction * change;
        if (!std.math.isFinite(candidate.*) or candidate.* < -1e-12) return error.InvalidSoilOrganicTransportCandidate;
        candidate.* = @max(0, candidate.*);
    }
}

fn substrateBase(layer: usize, substrate: usize) usize {
    return layer * component_count + substrate * components_per_substrate;
}

fn validateAdvance(state: *const State, profile: *const organic.State, faces: *const hydrology_module.SoilFaces, inputs: Inputs, options: Options) !void {
    try validateProfileDimensions(state, profile);
    const layers = state.layer_count;
    const amounts = layers * component_count;
    const face_components = faces.micropore_faces.len * component_count;
    if (faces.macropore_faces.len != faces.micropore_faces.len or inputs.micropore_conductance_m3_per_step.len != face_components or inputs.macropore_conductance_m3_per_step.len != face_components or inputs.matrix_water_m3.len != layers or inputs.macropore_water_m3.len != layers or inputs.layer_bulk_volume_m3.len != layers or inputs.micropore_external_water_flux_m3_per_step.len != layers or inputs.macropore_external_water_flux_m3_per_step.len != layers or inputs.recharge_concentration_g_per_m3.len != amounts) return error.SoilOrganicTransportDimensionMismatch;
    if (!std.math.isFinite(options.absolute_tolerance_g) or options.absolute_tolerance_g <= 0 or
        !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or
        !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or
        !std.math.isFinite(options.maximum_convective_fraction) or options.maximum_convective_fraction < 0 or options.maximum_convective_fraction > 1 or
        !std.math.isFinite(options.pore_exchange_fraction) or options.pore_exchange_fraction < 0 or options.pore_exchange_fraction > 1)
        return error.InvalidSoilOrganicTransportOptions;
    if (options.max_iterations == 0) return error.InvalidSoilOrganicTransportOptions;
    for (inputs.matrix_water_m3, inputs.macropore_water_m3, inputs.layer_bulk_volume_m3, inputs.micropore_external_water_flux_m3_per_step, inputs.macropore_external_water_flux_m3_per_step) |micro, macro, bulk, micro_boundary, macro_boundary| {
        inline for (.{ micro, macro, bulk, micro_boundary, macro_boundary }) |value| if (!std.math.isFinite(value)) return error.InvalidSoilOrganicTransportInput;
        if (micro < 0 or macro < 0 or bulk <= 0) return error.InvalidSoilOrganicTransportInput;
    }
    for (inputs.recharge_concentration_g_per_m3) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilOrganicTransportInput;
}

fn validateProfileDimensions(state: *const State, profile: *const organic.State) !void {
    if (profile.layer_count != state.layer_count or profile.dissolved.len != state.layer_count * organic.substrate_count or profile.dissolved_acetate_carbon_g_c.len != profile.dissolved.len) return error.SoilOrganicTransportDimensionMismatch;
}

fn validatePool(pool: organic.ElementPool, acetate: f64) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p, acetate }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilOrganicTransportPool;
}

fn allocateZero(allocator: std.mem.Allocator, count: usize) ![]f64 {
    const values = try allocator.alloc(f64, count);
    @memset(values, 0);
    return values;
}

test "runtime organic pore transport commits exact signed boundary loss" {
    var profile = try organic.State.init(std.testing.allocator, 1);
    defer profile.deinit();
    profile.dissolved[0] = .{ .carbon_g_c = 8, .nitrogen_g_n = 4, .phosphorus_g_p = 2 };
    profile.dissolved_acetate_carbon_g_c[0] = 6;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try state.initializeFromProfile(&profile);
    state.macropore_amount_g[@intFromEnum(Component.dissolved_organic_carbon)] = 4;

    const config = try @import("config.zig").SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var model_grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer model_grid.deinit();
    @memset(model_grid.active_soil_layer_count, 1);
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer hydrology.deinit();
    var no_faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &model_grid);
    defer no_faces.deinit();
    const zero_recharge = [_]f64{0} ** component_count;
    const result = try advance(std.testing.allocator, &state, &profile, &no_faces, .{
        .micropore_conductance_m3_per_step = &.{},
        .macropore_conductance_m3_per_step = &.{},
        .matrix_water_m3 = &.{1},
        .macropore_water_m3 = &.{1},
        .layer_bulk_volume_m3 = &.{1},
        .micropore_external_water_flux_m3_per_step = &.{0.25},
        .macropore_external_water_flux_m3_per_step = &.{0},
        .recharge_concentration_g_per_m3 = &zero_recharge,
    }, .{
        .absolute_tolerance_g = 1e-12,
        .relative_tolerance = 1e-8,
        .picard_relaxation = 0.5,
        .max_iterations = 20,
        .pore_exchange_fraction = 0,
    });
    try std.testing.expectEqual(@as(u16, 1), result.micropore_iterations);
    try std.testing.expectApproxEqAbs(@as(f64, 6), profile.dissolved[0].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), profile.dissolved[0].nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), profile.dissolved[0].phosphorus_g_p, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 4.5), profile.dissolved_acetate_carbon_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -2), state.boundary_net_flux_g[@intFromEnum(Component.dissolved_organic_carbon)], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 4), state.macropore_amount_g[@intFromEnum(Component.dissolved_organic_carbon)], 1e-15);
}

test "XFRS pore exchange conserves organic mass" {
    const exchange = try poreExchange(2, 6, 4, 3, 20, 0.25);
    try std.testing.expectApproxEqAbs(@as(f64, 8), (2 + exchange) + (6 - exchange), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.1), exchange, 1e-15);
}

test "accepted internal solve removes tolerance-sized component mass drift" {
    var base = [_]f64{0} ** (2 * component_count);
    var current = [_]f64{0} ** (2 * component_count);
    const doc = @intFromEnum(Component.dissolved_organic_carbon);
    const humus_doc = components_per_substrate * 4 + doc;
    base[doc] = 8;
    base[component_count + doc] = 2;
    current[doc] = 6.000_004;
    current[component_count + doc] = 4.000_003;
    base[humus_doc] = 5;
    current[humus_doc] = 4.999_996;

    try enforceInternalConservation(&base, &current);

    try std.testing.expectEqual(@as(f64, 10), current[doc] + current[component_count + doc]);
    try std.testing.expectEqual(@as(f64, 5), current[humus_doc] + current[component_count + humus_doc]);
}
