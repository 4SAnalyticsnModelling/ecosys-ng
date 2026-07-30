const std = @import("std");
const RootState = @import("plant_root_system.zig").State;

pub const species_count: usize = 8;

pub const Parameters = struct {
    reference_temperature_k: f64,
    aqueous_diffusivity_m2_per_h_at_reference: [species_count]f64,
    aqueous_diffusivity_temperature_exponent: f64,
    root_concentration_inhibition_mol_per_m3: [species_count]f64,

    pub fn validate(self: Parameters) !void {
        if (!std.math.isFinite(self.reference_temperature_k) or self.reference_temperature_k <= 0 or
            !std.math.isFinite(self.aqueous_diffusivity_temperature_exponent) or self.aqueous_diffusivity_temperature_exponent < 0)
            return error.InvalidRootSaltParameter;
        for (self.aqueous_diffusivity_m2_per_h_at_reference) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootSaltParameter;
        for (self.root_concentration_inhibition_mol_per_m3) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidRootSaltParameter;
    }

    pub fn diffusivityM2PerH(self: Parameters, species: usize, temperature_k: f64) !f64 {
        try self.validate();
        if (species >= species_count) return error.RootSaltSpeciesOutOfBounds;
        if (!std.math.isFinite(temperature_k) or temperature_k <= 0) return error.InvalidRootSaltTemperature;
        const value = self.aqueous_diffusivity_m2_per_h_at_reference[species] *
            std.math.pow(f64, temperature_k / self.reference_temperature_k, self.aqueous_diffusivity_temperature_exponent);
        if (!std.math.isFinite(value)) return error.NonFiniteRootSaltDiffusivity;
        return value;
    }
};

pub fn compatibilityParameters() Parameters {
    return .{
        .reference_temperature_k = 298.15,
        .aqueous_diffusivity_m2_per_h_at_reference = .{5.0e-6} ** species_count,
        .aqueous_diffusivity_temperature_exponent = 6,
        .root_concentration_inhibition_mol_per_m3 = .{ 1.0e-5, 1.0e-3, 1, 1, 1.0e-3, 1, 1, 1.0e-3 },
    };
}

pub const Input = struct {
    soil_content_mol: f64,
    root_content_mol: f64,
    soil_water_volume_m3: f64,
    root_water_volume_m3: f64,
    water_advection_m3_per_step: f64,
    diffusive_conductance_m3_per_step: f64,
    plant_population_count: f64,
    equilibration_fraction: f64,
    root_concentration_inhibition_mol_per_m3: f64,
};

/// Exact dynamic-salt UPTAKE kernel shared by Al, Fe, Ca, Mg, Na, K, SO4,
/// and Cl. Positive result moves salt from soil to root; negative result is
/// release. The inhibition divisor is applied after the directional source
/// bounds, matching each RUPZ* equation.
pub fn calculateExchangeMol(input: Input) !f64 {
    inline for (@typeInfo(Input).@"struct".fields) |field| if (!std.math.isFinite(@field(input, field.name))) return error.NonFiniteRootSaltExchangeInput;
    if (input.soil_content_mol < 0 or input.root_content_mol < 0 or input.soil_water_volume_m3 <= 0 or input.root_water_volume_m3 <= 0 or input.water_advection_m3_per_step < 0 or input.diffusive_conductance_m3_per_step < 0 or input.plant_population_count <= 0 or input.equilibration_fraction < 0 or input.equilibration_fraction > 1 or input.root_concentration_inhibition_mol_per_m3 <= 0) return error.InvalidRootSaltExchangeInput;
    const soil_concentration = input.soil_content_mol / input.soil_water_volume_m3;
    const root_concentration = input.root_content_mol / input.root_water_volume_m3;
    const candidate = (input.water_advection_m3_per_step * soil_concentration +
        input.diffusive_conductance_m3_per_step * (soil_concentration - root_concentration)) * input.plant_population_count;
    const total_water = input.root_water_volume_m3 + input.soil_water_volume_m3;
    const equilibrium_extent = (input.root_water_volume_m3 * input.soil_content_mol - input.soil_water_volume_m3 * input.root_content_mol) /
        total_water * input.equilibration_fraction;
    const bounded = if (candidate > 0) @min(@max(0.0, equilibrium_extent), candidate) else @max(@min(0.0, equilibrium_extent), candidate);
    const result = bounded / (1.0 + root_concentration / input.root_concentration_inhibition_mol_per_m3);
    if (!std.math.isFinite(result)) return error.NonFiniteRootSaltExchange;
    return result;
}

pub fn commitExchangeMol(soil_content_mol: *f64, root_content_mol: *f64, accumulated_root_uptake_mol: *f64, exchange_mol: f64) !void {
    inline for (.{ soil_content_mol.*, root_content_mol.*, accumulated_root_uptake_mol.*, exchange_mol }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootSaltCommitInput;
    if (soil_content_mol.* < 0 or root_content_mol.* < 0) return error.InvalidRootSaltCommitInput;
    const next_soil = soil_content_mol.* - exchange_mol;
    const next_root = root_content_mol.* + exchange_mol;
    const next_uptake = accumulated_root_uptake_mol.* + exchange_mol;
    if (next_soil < -1.0e-12 or next_root < -1.0e-12) return error.InsufficientSaltForRootExchange;
    if (!std.math.isFinite(next_soil) or !std.math.isFinite(next_root) or !std.math.isFinite(next_uptake)) return error.NonFiniteRootSaltCommit;
    soil_content_mol.* = @max(0.0, next_soil);
    root_content_mol.* = @max(0.0, next_root);
    accumulated_root_uptake_mol.* = next_uptake;
}

pub const LayerInputs = struct {
    soil_content_mol: []f64,
    root_content_mol: []f64,
    accumulated_root_uptake_mol: []f64,
    diffusive_conductance_m3_per_step: []const f64,
    root_concentration_inhibition_mol_per_m3: []const f64,
    soil_water_volume_m3: f64,
    root_water_volume_m3: f64,
    water_advection_m3_per_step: f64,
    plant_population_count: f64,
    equilibration_fraction: f64,
};

pub const LayerCompetitor = struct {
    plant: usize,
    domain: usize,
    layer: usize,
    soil_inventory_fraction: f64,
    root_water_volume_m3: f64,
    water_advection_m3_per_step: f64,
    diffusive_geometry_m: f64,
    plant_population_count: f64,
};

pub const SolverOptions = struct {
    absolute_tolerance_mol: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    max_iterations: u16,
};

pub const SolverReport = struct {
    iterations: u16,
    newton_raphson_steps: usize,
    picard_steps: usize,
};

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    competitor_capacity: usize,
    competitors: []LayerCompetitor,
    staged_exchange_mol: []f64,
    candidate_exchange_mol: []f64,

    pub fn init(allocator: std.mem.Allocator, competitor_capacity: usize) !Workspace {
        if (competitor_capacity == 0) return error.ZeroRootSaltCompetitorCapacity;
        const competitors = try allocator.alloc(LayerCompetitor, competitor_capacity);
        errdefer allocator.free(competitors);
        const staged = try allocator.alloc(f64, try std.math.mul(usize, competitor_capacity, species_count));
        errdefer allocator.free(staged);
        const candidate = try allocator.alloc(f64, try std.math.mul(usize, competitor_capacity, species_count));
        @memset(staged, 0);
        @memset(candidate, 0);
        return .{ .allocator = allocator, .competitor_capacity = competitor_capacity, .competitors = competitors, .staged_exchange_mol = staged, .candidate_exchange_mol = candidate };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.candidate_exchange_mol);
        self.allocator.free(self.staged_exchange_mol);
        self.allocator.free(self.competitors);
        self.* = undefined;
    }

    pub fn advance(
        self: *Workspace,
        roots: *RootState,
        soil_content_mol: []f64,
        soil_water_volume_m3: f64,
        temperature_k: f64,
        parameters: Parameters,
        competitor_count: usize,
        options: SolverOptions,
    ) !SolverReport {
        if (competitor_count > self.competitor_capacity) return error.RootSaltCompetitorOutOfBounds;
        return advanceCompetingLayer(
            roots,
            soil_content_mol,
            soil_water_volume_m3,
            temperature_k,
            parameters,
            self.competitors[0..competitor_count],
            self.staged_exchange_mol[0 .. competitor_count * species_count],
            self.candidate_exchange_mol[0 .. competitor_count * species_count],
            options,
        );
    }
};

pub const GridWorkspace = struct {
    allocator: std.mem.Allocator,
    per_cell: []Workspace,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, competitor_capacity_per_cell: usize) !GridWorkspace {
        if (cell_count == 0) return error.ZeroRootSaltGridCells;
        const cells = try allocator.alloc(Workspace, cell_count);
        errdefer allocator.free(cells);
        var initialized: usize = 0;
        errdefer for (cells[0..initialized]) |*cell| cell.deinit();
        for (cells) |*cell| {
            cell.* = try Workspace.init(allocator, competitor_capacity_per_cell);
            initialized += 1;
        }
        return .{ .allocator = allocator, .per_cell = cells };
    }

    pub fn deinit(self: *GridWorkspace) void {
        for (self.per_cell) |*cell| cell.deinit();
        self.allocator.free(self.per_cell);
        self.* = undefined;
    }
};

/// Stages every runtime root competitor against one immutable layer snapshot.
/// This is the allocation-free hourly replacement for traversal-ordered RUPZ*
/// publication in UPTAKE.
pub fn advanceCompetingLayer(
    roots: *RootState,
    soil_content_mol: []f64,
    soil_water_volume_m3: f64,
    temperature_k: f64,
    parameters: Parameters,
    competitors: []const LayerCompetitor,
    staged_exchange_mol: []f64,
    candidate_exchange_mol: []f64,
    options: SolverOptions,
) !SolverReport {
    if (soil_content_mol.len != species_count or staged_exchange_mol.len != competitors.len * species_count or candidate_exchange_mol.len != staged_exchange_mol.len) return error.RootSaltSpeciesCountMismatch;
    if (!std.math.isFinite(soil_water_volume_m3) or soil_water_volume_m3 <= 0) return error.InvalidRootSaltExchangeInput;
    if (!std.math.isFinite(options.absolute_tolerance_mol) or options.absolute_tolerance_mol <= 0 or
        !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or
        !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or
        options.max_iterations == 0)
        return error.InvalidRootSaltSolverOptions;
    try parameters.validate();
    for (competitors, 0..) |competitor, competitor_index| {
        if (!std.math.isFinite(competitor.soil_inventory_fraction) or competitor.soil_inventory_fraction < 0 or competitor.soil_inventory_fraction > 1 or
            !std.math.isFinite(competitor.diffusive_geometry_m) or competitor.diffusive_geometry_m < 0)
            return error.InvalidRootSaltCompetitor;
        _ = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
        _ = competitor_index;
    }
    @memset(staged_exchange_mol, 0);
    var converged = false;
    var iterations: u16 = 0;
    var newton_raphson_steps: usize = 0;
    var picard_steps: usize = 0;
    for (0..options.max_iterations) |iteration| {
        iterations = @intCast(iteration + 1);
        var maximum_scaled_residual: f64 = 0;
        for (0..species_count) |species| {
            var total_exchange: f64 = 0;
            for (0..competitors.len) |competitor_index| total_exchange += staged_exchange_mol[competitor_index * species_count + species];
            for (competitors, 0..) |competitor, competitor_index| {
                const index = competitor_index * species_count + species;
                const proposal = try implicitExchangeProposal(roots, soil_content_mol[species], soil_water_volume_m3, temperature_k, parameters, competitor, species, staged_exchange_mol[index], total_exchange);
                const residual = proposal - staged_exchange_mol[index];
                const scale = options.absolute_tolerance_mol + options.relative_tolerance * @max(@abs(proposal), @abs(staged_exchange_mol[index]));
                maximum_scaled_residual = @max(maximum_scaled_residual, @abs(residual) / scale);
                const probe = @max(options.absolute_tolerance_mol, 1.0e-6 * @max(1, @abs(staged_exchange_mol[index])));
                const perturbed_proposal = try implicitExchangeProposal(roots, soil_content_mol[species], soil_water_volume_m3, temperature_k, parameters, competitor, species, staged_exchange_mol[index] + probe, total_exchange + probe);
                const derivative = (perturbed_proposal - proposal) / probe - 1;
                const newton = if (std.math.isFinite(derivative) and @abs(derivative) > 1.0e-12)
                    staged_exchange_mol[index] - residual / derivative
                else
                    std.math.nan(f64);
                candidate_exchange_mol[index] = if (std.math.isFinite(newton)) value: {
                    newton_raphson_steps += 1;
                    break :value newton;
                } else value: {
                    picard_steps += 1;
                    break :value staged_exchange_mol[index] + options.picard_relaxation * residual;
                };
            }
            constrainSpeciesExchange(roots, soil_content_mol[species], competitors, species, candidate_exchange_mol);
        }
        if (maximum_scaled_residual <= 1) {
            converged = true;
            break;
        }
        @memcpy(staged_exchange_mol, candidate_exchange_mol);
    }
    if (!converged) return error.RootSaltSolverDidNotConverge;

    // Candidate is the converged Newton/Picard state; the initial zero state
    // can converge without entering the copy at the end of an iteration.
    @memcpy(staged_exchange_mol, candidate_exchange_mol);
    var total_exchange = [_]f64{0} ** species_count;
    for (0..species_count) |species| for (0..competitors.len) |competitor_index| {
        total_exchange[species] += staged_exchange_mol[competitor_index * species_count + species];
    };
    for (0..species_count) |species| {
        if (!std.math.isFinite(soil_content_mol[species]) or soil_content_mol[species] < 0) return error.InvalidRootSaltCommitInput;
        if (total_exchange[species] > soil_content_mol[species] + 1.0e-12) return error.RootSaltCompetitionOverdraw;
    }
    for (competitors, 0..) |competitor, competitor_index| {
        const root_layer = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
        const base = competitor_index * species_count;
        for (0..species_count) |species| {
            const salt_index = root_layer * species_count + species;
            const next_root = roots.salt_content_mol[salt_index] + staged_exchange_mol[base + species];
            const next_uptake = roots.salt_uptake_mol_per_h[salt_index] + staged_exchange_mol[base + species];
            if (!std.math.isFinite(next_root) or !std.math.isFinite(next_uptake) or next_root < -1.0e-12) return error.InsufficientSaltForRootExchange;
        }
    }
    for (0..species_count) |species| soil_content_mol[species] = @max(0, soil_content_mol[species] - total_exchange[species]);
    for (competitors, 0..) |competitor, competitor_index| {
        const root_layer = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
        const base = competitor_index * species_count;
        for (0..species_count) |species| {
            const salt_index = root_layer * species_count + species;
            roots.salt_content_mol[salt_index] = @max(0, roots.salt_content_mol[salt_index] + staged_exchange_mol[base + species]);
            roots.salt_uptake_mol_per_h[salt_index] += staged_exchange_mol[base + species];
        }
    }
    return .{ .iterations = iterations, .newton_raphson_steps = newton_raphson_steps, .picard_steps = picard_steps };
}

fn implicitExchangeProposal(
    roots: *const RootState,
    base_soil_content_mol: f64,
    soil_water_volume_m3: f64,
    temperature_k: f64,
    parameters: Parameters,
    competitor: LayerCompetitor,
    species: usize,
    competitor_exchange_mol: f64,
    total_exchange_mol: f64,
) !f64 {
    const root_layer = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
    const root_content = roots.salt_content_mol[root_layer * species_count + species];
    const final_soil = @max(0, base_soil_content_mol - total_exchange_mol);
    const final_root = @max(0, root_content + competitor_exchange_mol);
    const soil_concentration = final_soil * competitor.soil_inventory_fraction / soil_water_volume_m3;
    const root_concentration = final_root / competitor.root_water_volume_m3;
    const conductance = try parameters.diffusivityM2PerH(species, temperature_k) * competitor.diffusive_geometry_m;
    const candidate = (competitor.water_advection_m3_per_step * soil_concentration +
        conductance * (soil_concentration - root_concentration)) * competitor.plant_population_count;
    const allocated_base_soil = base_soil_content_mol * competitor.soil_inventory_fraction;
    const equilibrium_extent = (competitor.root_water_volume_m3 * allocated_base_soil - soil_water_volume_m3 * root_content) /
        (competitor.root_water_volume_m3 + soil_water_volume_m3);
    const bounded = if (candidate > 0) @min(@max(0, equilibrium_extent), candidate) else @max(@min(0, equilibrium_extent), candidate);
    const proposal = bounded / (1 + root_concentration / parameters.root_concentration_inhibition_mol_per_m3[species]);
    if (!std.math.isFinite(proposal)) return error.NonFiniteRootSaltCompetition;
    return proposal;
}

fn constrainSpeciesExchange(roots: *const RootState, soil_content_mol: f64, competitors: []const LayerCompetitor, species: usize, candidate: []f64) void {
    var positive: f64 = 0;
    var negative: f64 = 0;
    for (competitors, 0..) |competitor, competitor_index| {
        const index = competitor_index * species_count + species;
        const root_layer = roots.layerIndex(competitor.plant, competitor.domain, competitor.layer) catch unreachable;
        candidate[index] = @max(-roots.salt_content_mol[root_layer * species_count + species], candidate[index]);
        if (candidate[index] > 0) positive += candidate[index] else negative += candidate[index];
    }
    const maximum_positive = soil_content_mol - negative;
    if (positive > maximum_positive and positive > 0) {
        const fraction = maximum_positive / positive;
        for (0..competitors.len) |competitor_index| {
            const index = competitor_index * species_count + species;
            if (candidate[index] > 0) candidate[index] *= fraction;
        }
    }
}

/// One rollback-safe UPTAKE transaction for the named eight-salt registry.
pub fn advanceLayer(inputs: LayerInputs) ![8]f64 {
    inline for (.{ inputs.soil_content_mol, inputs.root_content_mol, inputs.accumulated_root_uptake_mol, inputs.diffusive_conductance_m3_per_step, inputs.root_concentration_inhibition_mol_per_m3 }) |values| if (values.len != 8) return error.RootSaltSpeciesCountMismatch;
    var exchange: [8]f64 = undefined;
    for (0..8) |species| exchange[species] = try calculateExchangeMol(.{
        .soil_content_mol = inputs.soil_content_mol[species],
        .root_content_mol = inputs.root_content_mol[species],
        .soil_water_volume_m3 = inputs.soil_water_volume_m3,
        .root_water_volume_m3 = inputs.root_water_volume_m3,
        .water_advection_m3_per_step = inputs.water_advection_m3_per_step,
        .diffusive_conductance_m3_per_step = inputs.diffusive_conductance_m3_per_step[species],
        .plant_population_count = inputs.plant_population_count,
        .equilibration_fraction = inputs.equilibration_fraction,
        .root_concentration_inhibition_mol_per_m3 = inputs.root_concentration_inhibition_mol_per_m3[species],
    });

    var next_soil: [8]f64 = undefined;
    var next_root: [8]f64 = undefined;
    var next_uptake: [8]f64 = undefined;
    for (0..8) |species| {
        next_soil[species] = inputs.soil_content_mol[species] - exchange[species];
        next_root[species] = inputs.root_content_mol[species] + exchange[species];
        next_uptake[species] = inputs.accumulated_root_uptake_mol[species] + exchange[species];
        inline for (.{ next_soil[species], next_root[species], next_uptake[species] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootSaltCommit;
        if (next_soil[species] < -1.0e-12 or next_root[species] < -1.0e-12) return error.InsufficientSaltForRootExchange;
    }
    for (0..8) |species| {
        inputs.soil_content_mol[species] = @max(0.0, next_soil[species]);
        inputs.root_content_mol[species] = @max(0.0, next_root[species]);
        inputs.accumulated_root_uptake_mol[species] = next_uptake[species];
    }
    return exchange;
}

test "UPTAKE dynamic salt exchange preserves source inhibition and bounds" {
    const input = Input{ .soil_content_mol = 8, .root_content_mol = 1, .soil_water_volume_m3 = 2, .root_water_volume_m3 = 1, .water_advection_m3_per_step = 0.1, .diffusive_conductance_m3_per_step = 0.2, .plant_population_count = 2, .equilibration_fraction = 0.5, .root_concentration_inhibition_mol_per_m3 = 2 };
    const candidate = (0.1 * 4.0 + 0.2 * (4.0 - 1.0)) * 2.0;
    const extent = (1.0 * 8.0 - 2.0 * 1.0) / 3.0 * 0.5;
    try std.testing.expectApproxEqAbs(@min(candidate, extent) / 1.5, try calculateExchangeMol(input), 1.0e-12);
}

test "UPTAKE dynamic salt commit is conservative and atomic" {
    var soil: f64 = 2;
    var root: f64 = 1;
    var uptake: f64 = 0;
    try commitExchangeMol(&soil, &root, &uptake, 0.5);
    try std.testing.expectApproxEqAbs(@as(f64, 3), soil + root, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), uptake, 1.0e-12);
    try std.testing.expectError(error.InsufficientSaltForRootExchange, commitExchangeMol(&soil, &root, &uptake, 3));
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), soil, 1.0e-12);
}

test "UPTAKE eight-salt layer transaction conserves every named species" {
    var soil = [_]f64{8} ** 8;
    var root = [_]f64{1} ** 8;
    var uptake = [_]f64{0} ** 8;
    const conductance = [_]f64{0.2} ** 8;
    const inhibition = [_]f64{2} ** 8;
    const exchange = try advanceLayer(.{ .soil_content_mol = &soil, .root_content_mol = &root, .accumulated_root_uptake_mol = &uptake, .diffusive_conductance_m3_per_step = &conductance, .root_concentration_inhibition_mol_per_m3 = &inhibition, .soil_water_volume_m3 = 2, .root_water_volume_m3 = 1, .water_advection_m3_per_step = 0.1, .plant_population_count = 2, .equilibration_fraction = 0.5 });
    for (0..8) |species| {
        try std.testing.expectApproxEqAbs(@as(f64, 9), soil[species] + root[species], 1.0e-12);
        try std.testing.expectApproxEqAbs(exchange[species], uptake[species], 1.0e-12);
    }
}

test "dynamic salt competitors share one snapshot beyond legacy plant capacity" {
    var roots = try RootState.init(std.testing.allocator, 7, 1, 1);
    defer roots.deinit();
    var workspace = try Workspace.init(std.testing.allocator, 7);
    defer workspace.deinit();
    var soil = [_]f64{70} ** species_count;
    const before = soil;
    for (0..7) |plant| {
        const root = try roots.layerIndex(plant, 0, 0);
        roots.aqueous_volume_m3[root] = 0.1;
        workspace.competitors[plant] = .{
            .plant = plant,
            .domain = 0,
            .layer = 0,
            .soil_inventory_fraction = 1.0 / 7.0,
            .root_water_volume_m3 = 0.1,
            .water_advection_m3_per_step = 0.01,
            .diffusive_geometry_m = 10,
            .plant_population_count = 1,
        };
    }
    const report = try workspace.advance(&roots, &soil, 1, 298.15, compatibilityParameters(), 7, .{ .absolute_tolerance_mol = 1.0e-12, .relative_tolerance = 1.0e-10, .picard_relaxation = 0.5, .max_iterations = 40 });
    try std.testing.expect(report.iterations < 40);
    try std.testing.expect(report.newton_raphson_steps > 0);
    for (0..species_count) |species| {
        var root_total: f64 = 0;
        for (0..7) |plant| root_total += roots.salt_content_mol[(try roots.layerIndex(plant, 0, 0)) * species_count + species];
        try std.testing.expectApproxEqAbs(before[species], soil[species] + root_total, 1.0e-12);
    }
    const soil_before_failure = soil;
    const root_before_failure = roots.salt_content_mol[0];
    try std.testing.expectError(error.RootSaltSolverDidNotConverge, workspace.advance(&roots, &soil, 1, 298.15, compatibilityParameters(), 7, .{ .absolute_tolerance_mol = 1.0e-30, .relative_tolerance = 1.0e-30, .picard_relaxation = 0.5, .max_iterations = 1 }));
    try std.testing.expectEqualSlices(f64, &soil_before_failure, &soil);
    try std.testing.expectEqual(root_before_failure, roots.salt_content_mol[0]);
}
