const std = @import("std");
const builtin = @import("builtin");
const grid_module = @import("../../state/grid.zig");
const heat = @import("flux.zig");
const transport_hydrology = @import("../../transport/hydrology.zig");
const numerics = @import("../../core/numerics.zig");
const boundary_topology_module = @import("../profile/boundary_topology.zig");
const water_boundary = @import("../water/boundary.zig");
const enthalpy = @import("../water/enthalpy_balance.zig");
const retention = @import("../water/retention.zig");

pub const Face = struct {
    source_cell: usize,
    destination_cell: usize,
    source_path_length_m: f64,
    destination_path_length_m: f64,
    face_area_m2: f64,
};

pub const FaceGeometry = struct { source_path_length_m: []const f64, destination_path_length_m: []const f64, face_area_m2: []const f64 };

pub const GeothermalBoundary = struct {
    topology: *const boundary_topology_module.State,
    layer_bottom_depth_m: []const f64,
    lower_face_area_m2: []const f64,
    enabled_by_cell: []const bool,
    mean_annual_temperature_k_by_cell: []const f64,
    minimum_source_depth_m: f64,
    source_depth_below_profile_m: f64,
    conductivity_m_megajoules_per_h_k: f64,
    geothermal_flux_megajoules_per_m2_h: f64,
};

pub const DirichletThermalBoundaries = struct {
    cell_index: []const usize,
    temperature_k: []const f64,
    distance_from_cell_center_m: []const f64,
    face_area_m2: []const f64,
};

pub const EnthalpyCoupling = struct {
    matrix_liquid_water_m3: []const f64,
    matrix_ice_water_equivalent_m3: []const f64,
    porous_medium_volume_m3: []const f64,
    /// Current matrix pore capacity. Dynamic surface/soil geometry can change
    /// this after initialization; when supplied it is authoritative for the
    /// retention-domain volume used by the Dall'Amico closure.
    matrix_pore_capacity_m3: []const f64 = &.{},
    /// Empty derives the head from the accepted matrix liquid content and
    /// the same runtime original-MVG curve used by Richards flow.
    unfrozen_pressure_head_m: []const f64 = &.{},
    mualem_van_genuchten: []const retention.MualemVanGenuchtenParameters,
    gravitational_water_potential_mpa_per_m: f64,
    pure_water_melting_temperature_k: f64,
    ice_water_equivalent_heat_capacity_megajoules_per_m3_k: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
    solver_options: enthalpy.SolverOptions,
    macropore_liquid_water_m3: []const f64 = &.{},
    macropore_ice_water_equivalent_m3: []const f64 = &.{},
    macropore_porous_medium_volume_m3: []const f64 = &.{},
    macropore_unfrozen_pressure_head_m: []const f64 = &.{},
    macropore_mualem_van_genuchten: []const retention.MualemVanGenuchtenParameters = &.{},
};

pub const Properties = struct {
    heat_capacity_megajoules_per_k: []const f64,
    minimum_heat_capacity_megajoules_per_k: []const f64,
    bulk_density_megagrams_per_m3: []const f64,
    liquid_water_fraction: []const f64,
    ice_fraction: []const f64,
    air_fraction: []const f64,
    fraction_of_pore_volume_air_filled: []const f64,
    solid_conductivity_numerator_m_megajoules_per_h_k: []const f64,
    solid_conductivity_denominator: []const f64,
    is_top_soil_layer: []const bool,
    top_snow_heat_capacity_megajoules_per_k: []const f64,
    maximum_negligible_snow_heat_capacity_megajoules_per_k: []const f64,
    snow_storage_heat_flux_megajoules: []const f64,
    cell_heat_source_megajoules: []const f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    turbulence: heat.TurbulenceParameters,
    /// Physical duration represented by one solve. Ordinary model execution
    /// uses one hour; published fine-step validation supplies 10/3600.
    time_step_hours: f64 = 1,
    geothermal_boundary: ?GeothermalBoundary = null,
    dirichlet_thermal_boundaries: ?DirichletThermalBoundaries = null,
    enthalpy_coupling: ?EnthalpyCoupling = null,
};

pub const WaterHeatFluxes = struct {
    liquid_water_m3: []const f64,
    vapor_m3: []const f64,
    macropore_water_m3: []const f64,
};

pub const Options = struct {
    /// Runtime NPH ceiling; convergence exits without completing unused cycles.
    max_iterations: u16,
    absolute_tolerance_k: f64 = 1e-10,
    relative_tolerance: f64 = 1e-9,
    picard_relaxation: f64 = 0.5,
    directional_probe_fraction: f64 = 1.0,
    minimum_newton_fraction: f64 = 0.05,
    // Aitken/Newton acceleration of the Picard residual may legitimately
    // exceed one. Two is the stable upper bound used here; accepted steps must
    // still reduce the fully recomputed nonlinear residual.
    maximum_newton_fraction: f64 = 2.0,
    /// Runtime memory/performance tradeoff. Zero selects only the O(n)
    /// directional Newton/Picard path.
    dense_newton_max_components: usize = 256,
};

pub const Result = struct {
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
    maximum_scaled_residual: f64,
    boundary_heat_input_megajoules: f64,
    boundary_heat_output_megajoules: f64,
};

const PhaseBuffers = struct {
    matrix_liquid_m3: []f64,
    matrix_ice_m3: []f64,
    macropore_liquid_m3: []f64,
    macropore_ice_m3: []f64,
    macropore_enabled: bool,
};

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    face_count: usize,
    base: []f64,
    current: []f64,
    residual: []f64,
    target: []f64,
    scratch: []f64,
    trial_flux: []f64,
    face_buffer: []Face,
    matrix_liquid_m3: []f64,
    matrix_ice_m3: []f64,
    macropore_liquid_m3: []f64,
    macropore_ice_m3: []f64,
    probe: []f64,
    probe_residual: []f64,
    candidate: []f64,
    candidate_residual: []f64,
    jacobian: []f64,
    newton_delta: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        face_count: usize,
        dense_newton_max_components: usize,
    ) !Workspace {
        if (cell_count == 0) return error.InvalidSoilHeatWorkspaceSize;
        var result: Workspace = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        result.face_count = face_count;
        var allocated: usize = 0;
        errdefer result.freeAllocated(allocated);
        inline for (.{
            "base",
            "current",
            "residual",
            "target",
            "scratch",
            "matrix_liquid_m3",
            "matrix_ice_m3",
            "macropore_liquid_m3",
            "macropore_ice_m3",
            "probe",
            "probe_residual",
            "candidate",
            "candidate_residual",
            "newton_delta",
        }) |field_name| {
            @field(result, field_name) = try allocator.alloc(f64, cell_count);
            allocated += 1;
        }
        result.trial_flux = try allocator.alloc(f64, face_count);
        allocated += 1;
        result.face_buffer = try allocator.alloc(Face, face_count);
        allocated += 1;
        result.jacobian = try allocator.alloc(
            f64,
            if (cell_count <= dense_newton_max_components)
                try std.math.mul(usize, cell_count, cell_count)
            else
                0,
        );
        return result;
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.jacobian);
        self.allocator.free(self.face_buffer);
        self.allocator.free(self.trial_flux);
        inline for (.{
            "newton_delta",
            "candidate_residual",
            "candidate",
            "probe_residual",
            "probe",
            "macropore_ice_m3",
            "macropore_liquid_m3",
            "matrix_ice_m3",
            "matrix_liquid_m3",
            "scratch",
            "target",
            "residual",
            "current",
            "base",
        }) |field_name| self.allocator.free(@field(self, field_name));
        self.* = undefined;
    }

    fn freeAllocated(self: *Workspace, count: usize) void {
        const field_names = [_][]const u8{
            "base",
            "current",
            "residual",
            "target",
            "scratch",
            "matrix_liquid_m3",
            "matrix_ice_m3",
            "macropore_liquid_m3",
            "macropore_ice_m3",
            "probe",
            "probe_residual",
            "candidate",
            "candidate_residual",
            "newton_delta",
            "trial_flux",
            "face_buffer",
        };
        inline for (field_names, 0..) |field_name, index|
            if (index < count) self.allocator.free(@field(self, field_name));
    }
};

/// Solves the nonlinear WATSUB temperature balance over arbitrary runtime
/// faces. State and HFLWM-equivalent output remain unchanged on failure.
pub fn solve(allocator: std.mem.Allocator, grid: *grid_module.GridState, faces: []const Face, properties: Properties, water_fluxes: WaterHeatFluxes, heat_flux_megajoules: []f64, options: Options) !Result {
    var workspace = try Workspace.init(
        allocator,
        grid.layer_count,
        faces.len,
        options.dense_newton_max_components,
    );
    defer workspace.deinit();
    return solveWithWorkspace(
        &workspace,
        grid,
        faces,
        properties,
        water_fluxes,
        heat_flux_megajoules,
        options,
    );
}

pub fn solveWithWorkspace(
    workspace: *Workspace,
    grid: *grid_module.GridState,
    faces: []const Face,
    properties: Properties,
    water_fluxes: WaterHeatFluxes,
    heat_flux_megajoules: []f64,
    options: Options,
) !Result {
    try validateInputs(grid, faces, properties, water_fluxes, heat_flux_megajoules, options);
    const count = grid.layer_count;
    const use_dense_newton = count <= options.dense_newton_max_components;
    const required_jacobian_count =
        if (use_dense_newton)
            try std.math.mul(usize, count, count)
        else
            0;
    if (workspace.cell_count != count or
        workspace.face_count != faces.len or
        workspace.jacobian.len < required_jacobian_count)
        return error.SoilHeatWorkspaceDimensionMismatch;
    const base = workspace.base;
    const current = workspace.current;
    const residual = workspace.residual;
    const target = workspace.target;
    const scratch = workspace.scratch;
    const trial_flux = workspace.trial_flux;
    const probe = workspace.probe;
    const probe_residual = workspace.probe_residual;
    const candidate = workspace.candidate;
    const candidate_residual = workspace.candidate_residual;
    const jacobian = workspace.jacobian[0..required_jacobian_count];
    const newton_delta = workspace.newton_delta;
    @memcpy(base, grid.soil_temperature_k);
    @memcpy(current, base);
    const phase_buffers: PhaseBuffers = .{
        .matrix_liquid_m3 = workspace.matrix_liquid_m3,
        .matrix_ice_m3 = workspace.matrix_ice_m3,
        .macropore_liquid_m3 = workspace.macropore_liquid_m3,
        .macropore_ice_m3 = workspace.macropore_ice_m3,
        .macropore_enabled = if (properties.enthalpy_coupling) |coupling|
            coupling.macropore_mualem_van_genuchten.len != 0
        else
            false,
    };
    // Keep quadratic storage away from large/out-of-core grids. Modest
    // connected systems use the full numerical Jacobian; larger systems retain
    // the O(n) directional Newton/Picard path below.
    var newton_steps: u16 = 0;
    var picard_steps: u16 = 0;
    var iteration: u16 = 0;
    var final_norm: f64 = std.math.inf(f64);
    while (iteration < options.max_iterations) : (iteration += 1) {
        try residualAt(faces, properties, water_fluxes, base, current, target, residual, scratch, trial_flux, phase_buffers);
        const norm = try scaledNorm(current, residual, options);
        final_norm = norm;
        if (norm <= 1) {
            try residualAt(faces, properties, water_fluxes, base, current, target, residual, scratch, heat_flux_megajoules, phase_buffers);
            @memcpy(grid.soil_temperature_k, current);
            if (properties.enthalpy_coupling != null)
                try commitMatrixPhase(
                    grid,
                    phase_buffers,
                );
            try grid.validateFinite();
            const boundary_heat = try acceptedBoundaryHeat(
                properties,
                current,
                phase_buffers,
            );
            return .{ .iterations = iteration + 1, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = norm, .boundary_heat_input_megajoules = boundary_heat.input_megajoules, .boundary_heat_output_megajoules = boundary_heat.output_megajoules };
        }
        var accepted_newton = false;
        // The O(n) component/directional update is normally sufficient and is
        // much cheaper. Build the full Jacobian only for the last permitted
        // update, where it can prevent an otherwise false ceiling failure.
        if (use_dense_newton and iteration + 1 == options.max_iterations) {
            var jacobian_valid = true;
            for (0..count) |column| {
                const perturbation = std.math.cbrt(std.math.floatEps(f64)) * @max(1.0, @abs(current[column]));
                @memcpy(probe, current);
                probe[column] += perturbation;
                if (residualAt(faces, properties, water_fluxes, base, probe, target, probe_residual, scratch, trial_flux, phase_buffers)) |_| {
                    @memcpy(probe, current);
                    probe[column] -= perturbation;
                    if (probe[column] <= 0) {
                        for (0..count) |row| jacobian[row * count + column] = (probe_residual[row] - residual[row]) / perturbation;
                    } else if (residualAt(faces, properties, water_fluxes, base, probe, target, candidate_residual, scratch, trial_flux, phase_buffers)) |_| {
                        for (0..count) |row| jacobian[row * count + column] = (probe_residual[row] - candidate_residual[row]) / (2.0 * perturbation);
                    } else |_| {
                        jacobian_valid = false;
                        break;
                    }
                } else |_| {
                    jacobian_valid = false;
                    break;
                }
            }
            if (jacobian_valid) {
                for (residual, newton_delta) |value, *right_hand_side| right_hand_side.* = -value;
                if (numerics.solveDenseLinearSystem(jacobian, newton_delta, count)) {
                    var line_fraction: f64 = 1;
                    var line_search: u8 = 0;
                    while (line_search < 8) : (line_search += 1) {
                        var candidate_valid = true;
                        for (current, newton_delta, candidate) |value, delta, *next| {
                            next.* = value + line_fraction * delta;
                            if (!std.math.isFinite(next.*) or next.* <= 0) candidate_valid = false;
                        }
                        if (candidate_valid) {
                            if (residualAt(faces, properties, water_fluxes, base, candidate, target, candidate_residual, scratch, trial_flux, phase_buffers)) |_| {
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
        }
        if (accepted_newton) continue;
        if (addDirection(current, residual, options.directional_probe_fraction, probe)) |_| {
            if (residualAt(faces, properties, water_fluxes, base, probe, target, probe_residual, scratch, trial_flux, phase_buffers)) |_| {
                // A component-wise secant approximation preserves the runtime
                // memory bound while resolving the very different thermal
                // stiffnesses of snow, litter, and mineral-soil cells. The
                // former single scalar fraction forced every cell to advance
                // at the rate of the stiffest one and could miss the source
                // NPH ceiling by one otherwise unnecessary iteration.
                var component_step_valid = true;
                for (current, residual, probe_residual, candidate) |value, base_residual, sampled_residual, *next| {
                    const derivative_along_residual = (sampled_residual - base_residual) / options.directional_probe_fraction;
                    if (!std.math.isFinite(derivative_along_residual) or @abs(derivative_along_residual) <= std.math.floatEps(f64)) {
                        component_step_valid = false;
                        break;
                    }
                    const raw_fraction = -base_residual / derivative_along_residual;
                    if (!std.math.isFinite(raw_fraction) or raw_fraction <= 0) {
                        component_step_valid = false;
                        break;
                    }
                    const fraction = std.math.clamp(raw_fraction, options.minimum_newton_fraction, options.maximum_newton_fraction);
                    next.* = value + fraction * base_residual;
                    if (!std.math.isFinite(next.*) or next.* <= 0) {
                        component_step_valid = false;
                        break;
                    }
                }
                if (component_step_valid) {
                    if (residualAt(faces, properties, water_fluxes, base, candidate, target, candidate_residual, scratch, trial_flux, phase_buffers)) |_| {
                        if (try scaledNorm(candidate, candidate_residual, options) < norm) {
                            @memcpy(current, candidate);
                            newton_steps += 1;
                            accepted_newton = true;
                        }
                    } else |_| {}
                }
                if (accepted_newton) continue;
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
                        if (residualAt(faces, properties, water_fluxes, base, candidate, target, candidate_residual, scratch, trial_flux, phase_buffers)) |_| {
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
    // The last permitted Newton/Picard update must be tested before declaring
    // failure; otherwise max_iterations=N permits only N-1 useful updates.
    try residualAt(faces, properties, water_fluxes, base, current, target, residual, scratch, trial_flux, phase_buffers);
    final_norm = try scaledNorm(current, residual, options);
    if (final_norm <= 1) {
        try residualAt(faces, properties, water_fluxes, base, current, target, residual, scratch, heat_flux_megajoules, phase_buffers);
        @memcpy(grid.soil_temperature_k, current);
        if (properties.enthalpy_coupling != null)
            try commitMatrixPhase(
                grid,
                phase_buffers,
            );
        try grid.validateFinite();
        const boundary_heat = try acceptedBoundaryHeat(
            properties,
            current,
            phase_buffers,
        );
        return .{ .iterations = options.max_iterations, .newton_raphson_steps = newton_steps, .picard_steps = picard_steps, .maximum_scaled_residual = final_norm, .boundary_heat_input_megajoules = boundary_heat.input_megajoules, .boundary_heat_output_megajoules = boundary_heat.output_megajoules };
    }
    if (!builtin.is_test) std.log.err("soil heat Newton-Raphson/Picard did not converge: iterations={d} final_scaled_residual={e} newton_steps={d} picard_steps={d}", .{ options.max_iterations, final_norm, newton_steps, picard_steps });
    return error.SoilHeatSolverDidNotConverge;
}

const BoundaryHeat = struct { input_megajoules: f64, output_megajoules: f64 };

/// Re-evaluates only explicit external thermal faces at the accepted
/// temperature. `cell_heat_source_megajoules` is deliberately excluded because its
/// mapped production owner combines true sources with surface-to-soil
/// conduction; that conduction is internal to EXEC's landscape storage.
/// Internal faces and water-carried heat likewise cancel within the domain.
fn acceptedBoundaryHeat(
    properties: Properties,
    accepted_temperature_k: []const f64,
    phase_buffers: PhaseBuffers,
) !BoundaryHeat {
    var result: BoundaryHeat = .{ .input_megajoules = 0, .output_megajoules = 0 };
    if (properties.geothermal_boundary) |geothermal| {
        for (geothermal.topology.faces) |boundary_face| {
            if (!boundary_face.is_lower_boundary) continue;
            const horizontal_cell = boundary_face.cell_index;
            if (!geothermal.enabled_by_cell[horizontal_cell]) continue;
            const layer = boundary_face.layer_index;
            const lower_depth_m = geothermal.layer_bottom_depth_m[layer];
            const source_depth_m = @max(
                geothermal.minimum_source_depth_m,
                lower_depth_m + geothermal.source_depth_below_profile_m,
            );
            const deep_temperature_k =
                geothermal.mean_annual_temperature_k_by_cell[horizontal_cell] +
                geothermal.geothermal_flux_megajoules_per_m2_h * source_depth_m /
                    geothermal.conductivity_m_megajoules_per_h_k;
            const outward_heat_megajoules = try water_boundary.geothermalHeatFluxMj(
                accepted_temperature_k[layer],
                deep_temperature_k,
                geothermal.conductivity_m_megajoules_per_h_k,
                source_depth_m,
                lower_depth_m,
                geothermal.lower_face_area_m2[layer],
                properties.time_step_hours,
            );
            try accumulateSignedBoundaryHeat(&result, -outward_heat_megajoules);
        }
    }
    if (properties.dirichlet_thermal_boundaries) |boundaries| {
        for (
            boundaries.cell_index,
            boundaries.temperature_k,
            boundaries.distance_from_cell_center_m,
            boundaries.face_area_m2,
        ) |cell, boundary_temperature_k, distance_m, face_area_m2| {
            const temperature_difference_k =
                boundary_temperature_k - accepted_temperature_k[cell];
            const conductivity_m_megajoules_per_h_k =
                try heat.calculateCellConductivity(
                    cellConductivityInputs(
                        properties,
                        cell,
                        temperature_difference_k,
                        phase_buffers,
                    ),
                    properties.turbulence,
                );
            const inward_heat_megajoules =
                conductivity_m_megajoules_per_h_k * face_area_m2 *
                temperature_difference_k / distance_m *
                properties.time_step_hours;
            try accumulateSignedBoundaryHeat(&result, inward_heat_megajoules);
        }
    }
    return result;
}

fn accumulateSignedBoundaryHeat(
    result: *BoundaryHeat,
    signed_input_megajoules: f64,
) !void {
    if (!std.math.isFinite(signed_input_megajoules))
        return error.NonFiniteAcceptedBoundaryHeat;
    if (signed_input_megajoules >= 0)
        result.input_megajoules += signed_input_megajoules
    else
        result.output_megajoules -= signed_input_megajoules;
    if (!std.math.isFinite(result.input_megajoules) or
        !std.math.isFinite(result.output_megajoules))
        return error.NonFiniteAcceptedBoundaryHeat;
}

pub fn solveAndBindTransportFaces(allocator: std.mem.Allocator, grid: *grid_module.GridState, hydrology: *transport_hydrology.State, shared_faces: *transport_hydrology.SoilFaces, geometry: FaceGeometry, properties: Properties, options: Options) !Result {
    const count = shared_faces.micropore_faces.len;
    var workspace = try Workspace.init(
        allocator,
        grid.layer_count,
        count,
        options.dense_newton_max_components,
    );
    defer workspace.deinit();
    return solveAndBindTransportFacesWithWorkspace(
        &workspace,
        grid,
        hydrology,
        shared_faces,
        geometry,
        properties,
        options,
    );
}

pub fn solveAndBindTransportFacesWithWorkspace(
    workspace: *Workspace,
    grid: *grid_module.GridState,
    hydrology: *transport_hydrology.State,
    shared_faces: *transport_hydrology.SoilFaces,
    geometry: FaceGeometry,
    properties: Properties,
    options: Options,
) !Result {
    const count = shared_faces.micropore_faces.len;
    if (shared_faces.macropore_faces.len != count or shared_faces.micropore_water_flux_m3_per_step.len != count or shared_faces.macropore_water_flux_m3_per_step.len != count or shared_faces.vapor_flux_m3_per_step.len != count or shared_faces.heat_flux_megajoules_per_step.len != count or geometry.source_path_length_m.len != count or geometry.destination_path_length_m.len != count or geometry.face_area_m2.len != count) return error.SoilHeatFaceGeometryDimensionMismatch;
    if (workspace.face_count != count)
        return error.SoilHeatWorkspaceDimensionMismatch;
    const faces = workspace.face_buffer;
    for (faces, 0..) |*face, index| face.* = .{ .source_cell = shared_faces.micropore_faces[index].first_cell, .destination_cell = shared_faces.micropore_faces[index].second_cell, .source_path_length_m = geometry.source_path_length_m[index], .destination_path_length_m = geometry.destination_path_length_m[index], .face_area_m2 = geometry.face_area_m2[index] };
    const result = try solveWithWorkspace(workspace, grid, faces, properties, .{ .liquid_water_m3 = shared_faces.micropore_water_flux_m3_per_step, .vapor_m3 = shared_faces.vapor_flux_m3_per_step, .macropore_water_m3 = shared_faces.macropore_water_flux_m3_per_step }, shared_faces.heat_flux_megajoules_per_step, options);
    @memset(hydrology.heat_face_flux_megajoules_per_step, 0);
    for (faces, shared_faces.direction_axis, shared_faces.heat_flux_megajoules_per_step) |face, axis, flux| hydrology.heat_face_flux_megajoules_per_step[face.source_cell * 3 + axis] = flux;
    try hydrology.validateFinite();
    return result;
}

fn residualAt(
    faces: []const Face,
    properties: Properties,
    water_fluxes: WaterHeatFluxes,
    base: []const f64,
    trial: []const f64,
    target: []f64,
    residual: []f64,
    scratch: []f64,
    output_flux: []f64,
    phase_buffers: PhaseBuffers,
) !void {
    @memcpy(scratch, trial);
    @memcpy(target, base);
    if (properties.enthalpy_coupling) |coupling| {
        for (0..trial.len) |cell| {
            if (trial[cell] == base[cell]) {
                phase_buffers.matrix_liquid_m3[cell] =
                    coupling.matrix_liquid_water_m3[cell];
                phase_buffers.matrix_ice_m3[cell] =
                    coupling.matrix_ice_water_equivalent_m3[cell];
                if (phase_buffers.macropore_enabled) {
                    phase_buffers.macropore_liquid_m3[cell] =
                        coupling.macropore_liquid_water_m3[cell];
                    phase_buffers.macropore_ice_m3[cell] =
                        coupling.macropore_ice_water_equivalent_m3[cell];
                } else {
                    phase_buffers.macropore_liquid_m3[cell] = 0;
                    phase_buffers.macropore_ice_m3[cell] = 0;
                }
                continue;
            }
            const trial_phase = try enthalpy.stateAtTemperature(
                try enthalpyParameters(properties, coupling, cell),
                trial[cell],
            );
            phase_buffers.matrix_liquid_m3[cell] =
                trial_phase.liquid_water_m3;
            phase_buffers.matrix_ice_m3[cell] =
                trial_phase.ice_water_equivalent_m3;
            phase_buffers.macropore_liquid_m3[cell] =
                trial_phase.secondary_liquid_water_m3;
            phase_buffers.macropore_ice_m3[cell] =
                trial_phase.secondary_ice_water_equivalent_m3;
        }
    }
    for (target, properties.cell_heat_source_megajoules, properties.heat_capacity_megajoules_per_k) |*temperature, source_megajoules, capacity| temperature.* += source_megajoules / capacity;
    @memset(output_flux, 0);
    for (faces, 0..) |face, face_index| {
        const source = face.source_cell;
        const destination = face.destination_cell;
        const difference = scratch[source] - scratch[destination];
        if (difference == 0 and
            water_fluxes.liquid_water_m3[face_index] == 0 and
            water_fluxes.vapor_m3[face_index] == 0 and
            water_fluxes.macropore_water_m3[face_index] == 0)
        {
            output_flux[face_index] = 0;
            continue;
        }
        const source_conductivity = try heat.calculateCellConductivity(
            cellConductivityInputs(
                properties,
                source,
                difference,
                phase_buffers,
            ),
            properties.turbulence,
        );
        const destination_conductivity = try heat.calculateCellConductivity(
            cellConductivityInputs(
                properties,
                destination,
                difference,
                phase_buffers,
            ),
            properties.turbulence,
        );
        const flux = try heat.calculateFaceFlux(.{ .source_temperature_k = scratch[source], .destination_temperature_k = scratch[destination], .source_heat_capacity_megajoules_per_k = properties.heat_capacity_megajoules_per_k[source], .destination_heat_capacity_megajoules_per_k = properties.heat_capacity_megajoules_per_k[destination], .source_minimum_heat_capacity_megajoules_per_k = properties.minimum_heat_capacity_megajoules_per_k[source], .destination_minimum_heat_capacity_megajoules_per_k = properties.minimum_heat_capacity_megajoules_per_k[destination], .source_is_top_soil_layer = properties.is_top_soil_layer[source], .top_snow_heat_capacity_megajoules_per_k = properties.top_snow_heat_capacity_megajoules_per_k[source], .maximum_negligible_snow_heat_capacity_megajoules_per_k = properties.maximum_negligible_snow_heat_capacity_megajoules_per_k[source], .snow_storage_heat_flux_megajoules = properties.snow_storage_heat_flux_megajoules[source], .liquid_water_flux_m3 = water_fluxes.liquid_water_m3[face_index], .vapor_flux_m3 = water_fluxes.vapor_m3[face_index], .macropore_water_flux_m3 = water_fluxes.macropore_water_m3[face_index], .liquid_water_heat_capacity_megajoules_per_m3_k = properties.liquid_water_heat_capacity_megajoules_per_m3_k, .source_thermal_conductivity_m_megajoules_per_h_k = source_conductivity, .destination_thermal_conductivity_m_megajoules_per_h_k = destination_conductivity, .source_path_length_m = face.source_path_length_m, .destination_path_length_m = face.destination_path_length_m, .face_area_m2 = face.face_area_m2, .time_fraction = properties.time_step_hours });
        output_flux[face_index] = flux.total_megajoules;
        const source_delta = flux.total_megajoules / properties.heat_capacity_megajoules_per_k[source];
        const destination_delta = flux.total_megajoules / properties.heat_capacity_megajoules_per_k[destination];
        scratch[source] -= source_delta;
        scratch[destination] += destination_delta;
        target[source] -= source_delta;
        target[destination] += destination_delta;
    }
    if (properties.geothermal_boundary) |geothermal| {
        for (geothermal.topology.faces) |boundary_face| {
            if (!boundary_face.is_lower_boundary) continue;
            const horizontal_cell = boundary_face.cell_index;
            if (!geothermal.enabled_by_cell[horizontal_cell]) continue;
            const layer = boundary_face.layer_index;
            const lower_depth_m = geothermal.layer_bottom_depth_m[layer];
            const source_depth_m = @max(geothermal.minimum_source_depth_m, lower_depth_m + geothermal.source_depth_below_profile_m);
            const deep_temperature_k = geothermal.mean_annual_temperature_k_by_cell[horizontal_cell] + geothermal.geothermal_flux_megajoules_per_m2_h * source_depth_m / geothermal.conductivity_m_megajoules_per_h_k;
            const outward_heat_megajoules = try water_boundary.geothermalHeatFluxMj(scratch[layer], deep_temperature_k, geothermal.conductivity_m_megajoules_per_h_k, source_depth_m, lower_depth_m, geothermal.lower_face_area_m2[layer], properties.time_step_hours);
            target[layer] -= outward_heat_megajoules / properties.heat_capacity_megajoules_per_k[layer];
        }
    }
    if (properties.dirichlet_thermal_boundaries) |boundaries| {
        for (boundaries.cell_index, boundaries.temperature_k, boundaries.distance_from_cell_center_m, boundaries.face_area_m2) |cell, boundary_temperature_k, distance_m, face_area_m2| {
            const temperature_difference_k =
                boundary_temperature_k - scratch[cell];
            const conductivity_m_megajoules_per_h_k =
                try heat.calculateCellConductivity(
                    cellConductivityInputs(
                        properties,
                        cell,
                        temperature_difference_k,
                        phase_buffers,
                    ),
                    properties.turbulence,
                );
            const inward_heat_megajoules =
                conductivity_m_megajoules_per_h_k * face_area_m2 *
                temperature_difference_k / distance_m *
                properties.time_step_hours;
            if (!std.math.isFinite(inward_heat_megajoules))
                return error.NonFiniteDirichletSoilHeatFlux;
            target[cell] += inward_heat_megajoules /
                properties.heat_capacity_megajoules_per_k[cell];
        }
    }
    if (properties.enthalpy_coupling) |coupling| {
        for (0..target.len) |cell| {
            const non_phase_heat_megajoules =
                (target[cell] - base[cell]) *
                properties.heat_capacity_megajoules_per_k[cell];
            if (non_phase_heat_megajoules == 0) {
                target[cell] = base[cell];
                phase_buffers.matrix_liquid_m3[cell] =
                    coupling.matrix_liquid_water_m3[cell];
                phase_buffers.matrix_ice_m3[cell] =
                    coupling.matrix_ice_water_equivalent_m3[cell];
                if (phase_buffers.macropore_enabled) {
                    phase_buffers.macropore_liquid_m3[cell] =
                        coupling.macropore_liquid_water_m3[cell];
                    phase_buffers.macropore_ice_m3[cell] =
                        coupling.macropore_ice_water_equivalent_m3[cell];
                } else {
                    phase_buffers.macropore_liquid_m3[cell] = 0;
                    phase_buffers.macropore_ice_m3[cell] = 0;
                }
                continue;
            }
            const parameters = try enthalpyParameters(properties, coupling, cell);
            const base_state = enthalpy.stateAtTemperature(
                parameters,
                base[cell],
            ) catch |err| {
                std.log.err(
                    "soil enthalpy coupling rejected state: cell={d} domain=matrix temperature_k={e} total_water_equivalent_m3={e} porous_medium_volume_m3={e} residual_water_content_m3_per_m3={e} saturated_water_content_m3_per_m3={e} error={s}",
                    .{
                        cell,
                        base[cell],
                        parameters.total_water_equivalent_m3,
                        parameters.porous_medium_volume_m3,
                        parameters.mualem_van_genuchten
                            .residual_water_content_m3_per_m3,
                        parameters.mualem_van_genuchten
                            .saturated_water_content_m3_per_m3,
                        @errorName(err),
                    },
                );
                if (parameters.secondary_domain) |secondary| {
                    std.log.err(
                        "soil enthalpy coupling paired domain: cell={d} domain=macropore total_water_equivalent_m3={e} porous_medium_volume_m3={e} residual_water_content_m3_per_m3={e} saturated_water_content_m3_per_m3={e}",
                        .{
                            cell,
                            secondary.total_water_equivalent_m3,
                            secondary.porous_medium_volume_m3,
                            secondary.mualem_van_genuchten
                                .residual_water_content_m3_per_m3,
                            secondary.mualem_van_genuchten
                                .saturated_water_content_m3_per_m3,
                        },
                    );
                }
                return err;
            };
            // Retained as a validity gate on `parameters` at `base[cell]`.
            _ = base_state;
            // The base enthalpy must be valued from the liquid/ice split the
            // grid actually holds, not from the Dall'Amico equilibrium split
            // re-derived at `base[cell]`. When the incoming grid state is off
            // its equilibrium partition, the equilibrium-valued base silently
            // creates or destroys `L * (liquid_eq - liquid_actual)` of latent
            // energy every hour, which the landscape census (which values the
            // real split) then reports as unexplained heat. Measured at
            // -129237 m3 of spurious liquid over 216 h = -0.4952 MJ/m2,
            // 99.1% of the observed `spatial_heat` residual.
            const actual_liquid_water_m3 =
                coupling.matrix_liquid_water_m3[cell] +
                if (phase_buffers.macropore_enabled)
                    coupling.macropore_liquid_water_m3[cell]
                else
                    0;
            const actual_ice_water_equivalent_m3 =
                coupling.matrix_ice_water_equivalent_m3[cell] +
                if (phase_buffers.macropore_enabled)
                    coupling.macropore_ice_water_equivalent_m3[cell]
                else
                    0;
            const actual_base_enthalpy_megajoules =
                (parameters.dry_solid_heat_capacity_megajoules_per_k +
                    parameters.liquid_water_heat_capacity_megajoules_per_m3_k *
                        actual_liquid_water_m3 +
                    parameters.ice_water_equivalent_heat_capacity_megajoules_per_m3_k *
                        actual_ice_water_equivalent_m3) *
                    (base[cell] - parameters.pure_water_melting_temperature_k) +
                parameters.latent_heat_of_fusion_megajoules_per_m3 *
                    actual_liquid_water_m3;
            if (!std.math.isFinite(actual_base_enthalpy_megajoules))
                return error.NonFiniteSoilEnthalpyState;
            var enthalpy_solver_options = coupling.solver_options;
            enthalpy_solver_options.initial_temperature_k = trial[cell];
            const solved = try enthalpy.temperatureFromEnthalpy(
                parameters,
                actual_base_enthalpy_megajoules + non_phase_heat_megajoules,
                enthalpy_solver_options,
            );
            target[cell] = solved.state.temperature_k;
            phase_buffers.matrix_liquid_m3[cell] =
                solved.state.liquid_water_m3;
            phase_buffers.matrix_ice_m3[cell] =
                solved.state.ice_water_equivalent_m3;
            phase_buffers.macropore_liquid_m3[cell] =
                solved.state.secondary_liquid_water_m3;
            phase_buffers.macropore_ice_m3[cell] =
                solved.state.secondary_ice_water_equivalent_m3;
        }
    }
    for (target, trial, residual) |value, trial_value, *difference| {
        if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilHeatCandidate;
        difference.* = value - trial_value;
    }
}

/// The `unfrozen_pressure_head_m` that `soil_enthalpy_balance` and
/// `soil_water_phase_change.dallAmicoEquilibrium` require: the matric head the
/// layer's TOTAL water (liquid plus ice water equivalent) would hold on the
/// UNMODIFIED retention curve if none of it were frozen.
///
/// Both properties are load bearing, and getting either wrong is a conservation
/// defect rather than a cosmetic one.
///
/// **Total water, not liquid.** In Dall'Amico et al. (2011) the head is the
/// reference state of the freezing curve: `dallAmicoEquilibrium` depresses the
/// melting point from it via the Clapeyron exponent, walks the head down the
/// curve to the current temperature, and reads the liquid fraction back off it.
/// Pricing the LIQUID content there makes the reference state a function of the
/// answer it is supposed to determine, so as ice grows the reference dries, the
/// depressed melting point falls, and the equilibrium tracks a moving target
/// while the enthalpy that was already booked against the previous target is
/// not revisited. The pre-existing owner of this same physics,
/// `soil_phase_solver.freezeThawEquilibrium`, already prices total water.
///
/// **Unmodified curve, not an ice-shrunk one.** The head is consumed alongside
/// `mualem_van_genuchten`, which `soil_enthalpy_balance` passes through
/// unmodified to both `dallAmicoEquilibrium` and `waterCapacityPerM`. A head
/// produced from a curve with a different `saturated_water_content_m3_per_m3`
/// is not a point on the curve that then interprets it, so the enthalpy state
/// and its own temperature derivative are evaluated on two different
/// constitutive relations and the Newton iteration converges to a temperature
/// that satisfies neither.
///
/// Clamping to the saturated content is the same admissibility guard
/// `soil_phase_solver` applies. It is not a widened domain: the retention curve
/// is undefined above saturation, saturation is where the head is zero, and a
/// total water content above it is a water-balance question owned upstream, not
/// something to be absorbed by silently reshaping the curve here.
fn unfrozenPressureHeadM(
    parameters: retention.MualemVanGenuchtenParameters,
    total_water_equivalent_m3: f64,
    porous_medium_volume_m3: f64,
) !f64 {
    if (!std.math.isFinite(porous_medium_volume_m3) or porous_medium_volume_m3 <= 0)
        return error.InvalidCoupledSoilPorousMediumVolume;
    const total_water_content_m3_per_m3 =
        total_water_equivalent_m3 / porous_medium_volume_m3;
    return parameters.pressureHeadAtWaterContent(std.math.clamp(
        total_water_content_m3_per_m3,
        parameters.residual_water_content_m3_per_m3,
        parameters.saturated_water_content_m3_per_m3,
    ));
}

fn enthalpyParameters(
    properties: Properties,
    coupling: EnthalpyCoupling,
    cell: usize,
) !enthalpy.Parameters {
    const matrix_retention = coupling.mualem_van_genuchten[cell];
    const matrix_volume_m3 =
        if (coupling.matrix_pore_capacity_m3.len != 0)
            try porousMediumVolumeFromPoreCapacity(
                coupling.matrix_pore_capacity_m3[cell],
                matrix_retention.saturated_water_content_m3_per_m3,
            )
        else
            coupling.porous_medium_volume_m3[cell];
    const dry_solid_heat_capacity_megajoules_per_k =
        properties.heat_capacity_megajoules_per_k[cell] -
        properties.liquid_water_heat_capacity_megajoules_per_m3_k *
            coupling.matrix_liquid_water_m3[cell] -
        coupling.ice_water_equivalent_heat_capacity_megajoules_per_m3_k *
            coupling.matrix_ice_water_equivalent_m3[cell] -
        (if (coupling.macropore_mualem_van_genuchten.len != 0)
            properties.liquid_water_heat_capacity_megajoules_per_m3_k *
                coupling.macropore_liquid_water_m3[cell] +
                coupling.ice_water_equivalent_heat_capacity_megajoules_per_m3_k *
                    coupling.macropore_ice_water_equivalent_m3[cell]
        else
            0);
    if (!std.math.isFinite(dry_solid_heat_capacity_megajoules_per_k) or
        dry_solid_heat_capacity_megajoules_per_k <
            -64.0 * std.math.floatEps(f64) *
                @max(1.0, properties.heat_capacity_megajoules_per_k[cell]))
        return error.InvalidCoupledSoilDryHeatCapacity;
    return .{
        .porous_medium_volume_m3 = matrix_volume_m3,
        .total_water_equivalent_m3 = coupling.matrix_liquid_water_m3[cell] +
            coupling.matrix_ice_water_equivalent_m3[cell],
        .unfrozen_pressure_head_m = if (coupling.unfrozen_pressure_head_m.len != 0)
            coupling.unfrozen_pressure_head_m[cell]
        else
            try unfrozenPressureHeadM(
                matrix_retention,
                coupling.matrix_liquid_water_m3[cell] +
                    coupling.matrix_ice_water_equivalent_m3[cell],
                matrix_volume_m3,
            ),
        .gravitational_water_potential_mpa_per_m = coupling.gravitational_water_potential_mpa_per_m,
        .pure_water_melting_temperature_k = coupling.pure_water_melting_temperature_k,
        .dry_solid_heat_capacity_megajoules_per_k = @max(0, dry_solid_heat_capacity_megajoules_per_k),
        .liquid_water_heat_capacity_megajoules_per_m3_k = properties.liquid_water_heat_capacity_megajoules_per_m3_k,
        .ice_water_equivalent_heat_capacity_megajoules_per_m3_k = coupling.ice_water_equivalent_heat_capacity_megajoules_per_m3_k,
        .latent_heat_of_fusion_megajoules_per_m3 = coupling.latent_heat_of_fusion_megajoules_per_m3,
        .mualem_van_genuchten = matrix_retention,
        .secondary_domain = if (coupling.macropore_mualem_van_genuchten.len != 0 and
            coupling.macropore_porous_medium_volume_m3[cell] > 0)
            .{
                .porous_medium_volume_m3 = coupling.macropore_porous_medium_volume_m3[cell],
                .total_water_equivalent_m3 = coupling.macropore_liquid_water_m3[cell] +
                    coupling.macropore_ice_water_equivalent_m3[cell],
                .unfrozen_pressure_head_m = if (coupling.macropore_unfrozen_pressure_head_m.len != 0)
                    coupling.macropore_unfrozen_pressure_head_m[cell]
                else
                    try unfrozenPressureHeadM(
                        coupling.macropore_mualem_van_genuchten[cell],
                        coupling.macropore_liquid_water_m3[cell] +
                            coupling.macropore_ice_water_equivalent_m3[cell],
                        coupling.macropore_porous_medium_volume_m3[cell],
                    ),
                .mualem_van_genuchten = coupling.macropore_mualem_van_genuchten[cell],
            }
        else
            null,
    };
}

fn porousMediumVolumeFromPoreCapacity(
    pore_capacity_m3: f64,
    saturated_water_content_m3_per_m3: f64,
) !f64 {
    if (!std.math.isFinite(pore_capacity_m3) or pore_capacity_m3 <= 0 or
        !std.math.isFinite(saturated_water_content_m3_per_m3) or
        saturated_water_content_m3_per_m3 <= 0)
        return error.InvalidCoupledSoilPoreGeometry;
    const volume_m3 =
        pore_capacity_m3 / saturated_water_content_m3_per_m3;
    if (!std.math.isFinite(volume_m3) or volume_m3 <= 0)
        return error.InvalidCoupledSoilPoreGeometry;
    return volume_m3;
}

test "dynamic matrix pore capacity defines the Dall'Amico retention volume" {
    const pore_capacity_m3 = 588_380.0197233009;
    const saturated_water_content_m3_per_m3 = 0.5120612363148425;
    const volume_m3 = try porousMediumVolumeFromPoreCapacity(
        pore_capacity_m3,
        saturated_water_content_m3_per_m3,
    );
    try std.testing.expectApproxEqRel(
        pore_capacity_m3,
        volume_m3 * saturated_water_content_m3_per_m3,
        4 * std.math.floatEps(f64),
    );
    try std.testing.expectError(
        error.InvalidCoupledSoilPoreGeometry,
        porousMediumVolumeFromPoreCapacity(0, saturated_water_content_m3_per_m3),
    );
}

fn commitMatrixPhase(
    grid: *grid_module.GridState,
    phase_buffers: PhaseBuffers,
) !void {
    if (phase_buffers.matrix_liquid_m3.len != grid.layer_count or
        phase_buffers.matrix_ice_m3.len != grid.layer_count or
        phase_buffers.macropore_liquid_m3.len != grid.layer_count or
        phase_buffers.macropore_ice_m3.len != grid.layer_count)
        return error.SoilHeatPhaseCommitDimensionMismatch;
    @memcpy(
        grid.matrix_liquid_water_m3,
        phase_buffers.matrix_liquid_m3,
    );
    @memcpy(grid.matrix_ice_water_m3, phase_buffers.matrix_ice_m3);
    if (phase_buffers.macropore_enabled) {
        @memcpy(
            grid.macropore_liquid_water_m3,
            phase_buffers.macropore_liquid_m3,
        );
        @memcpy(
            grid.macropore_ice_water_m3,
            phase_buffers.macropore_ice_m3,
        );
    }
    for (0..grid.layer_count) |cell| {
        grid.liquid_water_m3[cell] =
            grid.matrix_liquid_water_m3[cell] +
            grid.macropore_liquid_water_m3[cell];
        grid.ice_water_m3[cell] =
            grid.matrix_ice_water_m3[cell] +
            grid.macropore_ice_water_m3[cell];
        grid.matrix_air_volume_m3[cell] = @max(
            0,
            grid.matrix_pore_capacity_m3[cell] -
                grid.matrix_liquid_water_m3[cell] -
                grid.matrix_ice_water_m3[cell],
        );
        grid.macropore_air_volume_m3[cell] = @max(
            0,
            grid.macropore_pore_capacity_m3[cell] -
                grid.macropore_liquid_water_m3[cell] -
                grid.macropore_ice_water_m3[cell],
        );
        grid.air_volume_m3[cell] =
            grid.matrix_air_volume_m3[cell] +
            grid.macropore_air_volume_m3[cell];
    }
}

fn cellConductivityInputs(
    properties: Properties,
    cell: usize,
    temperature_difference_k: f64,
    phase_buffers: PhaseBuffers,
) heat.CellConductivityInputs {
    var liquid_water_fraction = properties.liquid_water_fraction[cell];
    var ice_fraction = properties.ice_fraction[cell];
    if (properties.enthalpy_coupling) |coupling| {
        liquid_water_fraction +=
            (phase_buffers.matrix_liquid_m3[cell] -
                coupling.matrix_liquid_water_m3[cell]) /
            coupling.porous_medium_volume_m3[cell];
        ice_fraction +=
            (phase_buffers.matrix_ice_m3[cell] -
                coupling.matrix_ice_water_equivalent_m3[cell]) /
            coupling.porous_medium_volume_m3[cell];
        if (coupling.macropore_mualem_van_genuchten.len != 0 and
            coupling.macropore_porous_medium_volume_m3[cell] > 0)
        {
            liquid_water_fraction +=
                (phase_buffers.macropore_liquid_m3[cell] -
                    coupling.macropore_liquid_water_m3[cell]) /
                coupling.macropore_porous_medium_volume_m3[cell];
            ice_fraction +=
                (phase_buffers.macropore_ice_m3[cell] -
                    coupling.macropore_ice_water_equivalent_m3[cell]) /
                coupling.macropore_porous_medium_volume_m3[cell];
        }
        liquid_water_fraction = @max(0, liquid_water_fraction);
        ice_fraction = @max(0, ice_fraction);
    }
    return .{
        .bulk_density_megagrams_per_m3 = properties.bulk_density_megagrams_per_m3[cell],
        .liquid_water_fraction = liquid_water_fraction,
        .ice_fraction = ice_fraction,
        .air_fraction = properties.air_fraction[cell],
        .fraction_of_pore_volume_air_filled = properties.fraction_of_pore_volume_air_filled[cell],
        .solid_conductivity_numerator_m_megajoules_per_h_k = properties.solid_conductivity_numerator_m_megajoules_per_h_k[cell],
        .solid_conductivity_denominator = properties.solid_conductivity_denominator[cell],
        .temperature_difference_k = temperature_difference_k,
    };
}

fn addDirection(current: []const f64, direction: []const f64, fraction: f64, output: []f64) !void {
    for (current, direction, output) |value, delta, *candidate| {
        candidate.* = value + fraction * delta;
        if (!std.math.isFinite(candidate.*) or candidate.* <= 0) return error.InvalidSoilHeatCandidate;
    }
}

fn scaledNorm(state: []const f64, residual: []const f64, options: Options) !f64 {
    var maximum: f64 = 0;
    for (state, residual) |value, difference| {
        if (!std.math.isFinite(value) or value <= 0 or !std.math.isFinite(difference)) return error.NonFiniteSoilHeatSolverState;
        maximum = @max(maximum, @abs(difference) / (options.absolute_tolerance_k + options.relative_tolerance * @max(1.0, @abs(value))));
    }
    return maximum;
}

fn validateInputs(grid: *const grid_module.GridState, faces: []const Face, properties: Properties, water_fluxes: WaterHeatFluxes, heat_flux_megajoules: []const f64, options: Options) !void {
    const cells = grid.layer_count;
    inline for (@typeInfo(Properties).@"struct".fields) |field| {
        if (field.type == []const f64 and @field(properties, field.name).len != cells) return error.SoilHeatSolverDimensionMismatch;
        if (field.type == []const bool and @field(properties, field.name).len != cells) return error.SoilHeatSolverDimensionMismatch;
    }
    if (water_fluxes.liquid_water_m3.len != faces.len or water_fluxes.vapor_m3.len != faces.len or water_fluxes.macropore_water_m3.len != faces.len or heat_flux_megajoules.len != faces.len) return error.SoilHeatSolverDimensionMismatch;
    if (properties.geothermal_boundary) |geothermal| {
        if (geothermal.layer_bottom_depth_m.len != cells or geothermal.lower_face_area_m2.len != cells or geothermal.enabled_by_cell.len != geothermal.topology.water_table_mode.len or geothermal.mean_annual_temperature_k_by_cell.len != geothermal.topology.water_table_mode.len or !std.math.isFinite(geothermal.minimum_source_depth_m) or geothermal.minimum_source_depth_m <= 0 or !std.math.isFinite(geothermal.source_depth_below_profile_m) or geothermal.source_depth_below_profile_m <= 0 or !std.math.isFinite(geothermal.conductivity_m_megajoules_per_h_k) or geothermal.conductivity_m_megajoules_per_h_k <= 0 or !std.math.isFinite(geothermal.geothermal_flux_megajoules_per_m2_h)) return error.InvalidGeothermalBoundary;
        for (geothermal.mean_annual_temperature_k_by_cell) |temperature_k| if (!std.math.isFinite(temperature_k) or temperature_k <= 0) return error.InvalidGeothermalBoundary;
    }
    if (properties.dirichlet_thermal_boundaries) |boundaries| {
        const boundary_count = boundaries.cell_index.len;
        if (boundaries.temperature_k.len != boundary_count or
            boundaries.distance_from_cell_center_m.len != boundary_count or
            boundaries.face_area_m2.len != boundary_count)
            return error.SoilHeatBoundaryDimensionMismatch;
        for (boundaries.cell_index, boundaries.temperature_k, boundaries.distance_from_cell_center_m, boundaries.face_area_m2) |cell, temperature_k, distance_m, face_area_m2| {
            if (cell >= cells or
                !std.math.isFinite(temperature_k) or temperature_k <= 0 or
                !std.math.isFinite(distance_m) or distance_m <= 0 or
                !std.math.isFinite(face_area_m2) or face_area_m2 < 0)
                return error.InvalidDirichletSoilHeatBoundary;
        }
    }
    if (properties.enthalpy_coupling) |coupling| {
        inline for (.{
            coupling.matrix_liquid_water_m3,
            coupling.matrix_ice_water_equivalent_m3,
            coupling.porous_medium_volume_m3,
        }) |values| if (values.len != cells)
            return error.SoilHeatSolverDimensionMismatch;
        if (coupling.unfrozen_pressure_head_m.len != 0 and
            coupling.unfrozen_pressure_head_m.len != cells)
            return error.SoilHeatSolverDimensionMismatch;
        if (coupling.matrix_pore_capacity_m3.len != 0 and
            coupling.matrix_pore_capacity_m3.len != cells)
            return error.SoilHeatSolverDimensionMismatch;
        if (coupling.mualem_van_genuchten.len != cells)
            return error.SoilHeatSolverDimensionMismatch;
        const macropore_enabled =
            coupling.macropore_mualem_van_genuchten.len != 0;
        if (macropore_enabled) {
            inline for (.{
                coupling.macropore_liquid_water_m3,
                coupling.macropore_ice_water_equivalent_m3,
                coupling.macropore_porous_medium_volume_m3,
            }) |values| if (values.len != cells)
                return error.SoilHeatSolverDimensionMismatch;
            if (coupling.macropore_mualem_van_genuchten.len != cells or
                (coupling.macropore_unfrozen_pressure_head_m.len != 0 and
                    coupling.macropore_unfrozen_pressure_head_m.len != cells))
                return error.SoilHeatSolverDimensionMismatch;
            for (coupling.macropore_mualem_van_genuchten) |parameters|
                try parameters.validate();
        } else if (coupling.macropore_liquid_water_m3.len != 0 or
            coupling.macropore_ice_water_equivalent_m3.len != 0 or
            coupling.macropore_porous_medium_volume_m3.len != 0 or
            coupling.macropore_unfrozen_pressure_head_m.len != 0)
            return error.IncompleteSoilHeatMacroporeEnthalpyCoupling;
        inline for (.{
            coupling.gravitational_water_potential_mpa_per_m,
            coupling.pure_water_melting_temperature_k,
            coupling.ice_water_equivalent_heat_capacity_megajoules_per_m3_k,
            coupling.latent_heat_of_fusion_megajoules_per_m3,
        }) |value| if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSoilHeatEnthalpyCoupling;
        for (coupling.mualem_van_genuchten) |parameters|
            try parameters.validate();
    }
    if (!std.math.isFinite(properties.liquid_water_heat_capacity_megajoules_per_m3_k) or properties.liquid_water_heat_capacity_megajoules_per_m3_k <= 0) return error.InvalidSoilHeatCapacity;
    if (!std.math.isFinite(properties.time_step_hours) or
        properties.time_step_hours <= 0 or properties.time_step_hours > 1)
        return error.InvalidSoilHeatTimeStep;
    if (options.max_iterations == 0 or !std.math.isFinite(options.absolute_tolerance_k) or options.absolute_tolerance_k <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.directional_probe_fraction) or options.directional_probe_fraction <= 0 or !std.math.isFinite(options.minimum_newton_fraction) or options.minimum_newton_fraction <= 0 or !std.math.isFinite(options.maximum_newton_fraction) or options.maximum_newton_fraction < options.minimum_newton_fraction) return error.InvalidSoilHeatSolverOptions;
    for (0..cells) |cell| if (!std.math.isFinite(properties.heat_capacity_megajoules_per_k[cell]) or properties.heat_capacity_megajoules_per_k[cell] <= 0) return error.InvalidSoilHeatCapacity;
    for (faces) |face| if (face.source_cell >= cells or face.destination_cell >= cells or face.source_cell == face.destination_cell or !std.math.isFinite(face.source_path_length_m) or face.source_path_length_m <= 0 or !std.math.isFinite(face.destination_path_length_m) or face.destination_path_length_m <= 0 or !std.math.isFinite(face.face_area_m2) or face.face_area_m2 < 0) return error.InvalidSoilHeatFace;
}

fn testProperties() Properties {
    const values = struct {
        const capacity = [_]f64{ 2, 2 };
        const zero = [_]f64{ 0, 0 };
        const density = [_]f64{ 1, 1 };
        const liquid = [_]f64{ 0.2, 0.2 };
        const air = [_]f64{ 0.3, 0.3 };
        const numerator = [_]f64{ 0.01, 0.01 };
        const denominator = [_]f64{ 1, 1 };
        const top = [_]bool{ true, false };
    };
    return .{ .heat_capacity_megajoules_per_k = &values.capacity, .minimum_heat_capacity_megajoules_per_k = &values.zero, .bulk_density_megagrams_per_m3 = &values.density, .liquid_water_fraction = &values.liquid, .ice_fraction = &values.zero, .air_fraction = &values.air, .fraction_of_pore_volume_air_filled = &values.air, .solid_conductivity_numerator_m_megajoules_per_h_k = &values.numerator, .solid_conductivity_denominator = &values.denominator, .is_top_soil_layer = &values.top, .top_snow_heat_capacity_megajoules_per_k = &values.zero, .maximum_negligible_snow_heat_capacity_megajoules_per_k = &values.zero, .snow_storage_heat_flux_megajoules = &values.zero, .cell_heat_source_megajoules = &values.zero, .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19, .turbulence = .{ .water_fraction_threshold = 1, .air_fraction_threshold = 1, .water_rayleigh_coefficient = 0, .air_rayleigh_coefficient = 0, .water_nusselt_denominator = 1, .air_nusselt_denominator = 1 } };
}

test "hybrid heat solve exits before NPH and conserves sensible heat" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.soil_temperature_k[0] = 300;
    grid.soil_temperature_k[1] = 280;
    const zero_flux = [_]f64{0};
    var output = [_]f64{0};
    const before = 2 * grid.soil_temperature_k[0] + 2 * grid.soil_temperature_k[1];
    const result = try solve(std.testing.allocator, &grid, &.{.{ .source_cell = 0, .destination_cell = 1, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1 }}, testProperties(), .{ .liquid_water_m3 = &zero_flux, .vapor_m3 = &zero_flux, .macropore_water_m3 = &zero_flux }, &output, .{ .max_iterations = 20 });
    try std.testing.expect(result.iterations < 20);
    try std.testing.expect(result.newton_raphson_steps + result.picard_steps > 0);
    try std.testing.expectApproxEqAbs(before, 2 * grid.soil_temperature_k[0] + 2 * grid.soil_temperature_k[1], 1e-10);
    try std.testing.expect(output[0] > 0);
}

test "failed heat convergence is atomic" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.soil_temperature_k[0] = 300;
    grid.soil_temperature_k[1] = 280;
    const zero_flux = [_]f64{0};
    var output = [_]f64{99};
    try std.testing.expectError(error.SoilHeatSolverDidNotConverge, solve(std.testing.allocator, &grid, &.{.{ .source_cell = 0, .destination_cell = 1, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1 }}, testProperties(), .{ .liquid_water_m3 = &zero_flux, .vapor_m3 = &zero_flux, .macropore_water_m3 = &zero_flux }, &output, .{ .max_iterations = 1, .minimum_newton_fraction = 0.05, .maximum_newton_fraction = 0.05, .dense_newton_max_components = 0 }));
    try std.testing.expectEqual(@as(f64, 300), grid.soil_temperature_k[0]);
    try std.testing.expectEqual(@as(f64, 99), output[0]);
}

test "coupled heat residual commits Dall'Amico phase and conserves enthalpy" {
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
            .max_nonlinear_iterations = 80,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.soil_temperature_k[0] = 275;
    grid.matrix_liquid_water_m3[0] = 0.4;
    grid.macropore_liquid_water_m3[0] = 0.05;
    grid.liquid_water_m3[0] = 0.45;
    grid.matrix_pore_capacity_m3[0] = 0.5;
    grid.macropore_pore_capacity_m3[0] = 0.1;
    grid.matrix_air_volume_m3[0] = 0.1;
    grid.macropore_air_volume_m3[0] = 0.05;
    grid.air_volume_m3[0] = 0.15;
    const capacity = [_]f64{2.8855};
    const zero = [_]f64{0};
    const density = [_]f64{1};
    const liquid_fraction = [_]f64{0.45};
    const air_fraction = [_]f64{0.15};
    const numerator = [_]f64{0.01};
    const denominator = [_]f64{1};
    const top = [_]bool{true};
    const cooling_megajoules = [_]f64{-100};
    const volume = [_]f64{1};
    const macropore_volume = [_]f64{0.1};
    const curve = [_]retention.MualemVanGenuchtenParameters{.{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.5,
        .alpha_per_m = 1.6,
        .n = 1.6,
        .saturated_hydraulic_conductivity_m_per_h = 0.01,
    }};
    const macropore_curve = [_]retention.MualemVanGenuchtenParameters{.{
        .residual_water_content_m3_per_m3 = 0,
        .saturated_water_content_m3_per_m3 = 1,
        .alpha_per_m = 15,
        .n = 2.68,
        .saturated_hydraulic_conductivity_m_per_h = 0.1,
    }};
    const coupling: EnthalpyCoupling = .{
        .matrix_liquid_water_m3 = grid.matrix_liquid_water_m3,
        .matrix_ice_water_equivalent_m3 = grid.matrix_ice_water_m3,
        .porous_medium_volume_m3 = &volume,
        .mualem_van_genuchten = &curve,
        .gravitational_water_potential_mpa_per_m = 0.00980665,
        .pure_water_melting_temperature_k = 273.15,
        .ice_water_equivalent_heat_capacity_megajoules_per_m3_k = 1.93,
        .latent_heat_of_fusion_megajoules_per_m3 = 333.7,
        .solver_options = .{ .max_iterations = 80 },
        .macropore_liquid_water_m3 = grid.macropore_liquid_water_m3,
        .macropore_ice_water_equivalent_m3 = grid.macropore_ice_water_m3,
        .macropore_porous_medium_volume_m3 = &macropore_volume,
        .macropore_mualem_van_genuchten = &macropore_curve,
    };
    const properties: Properties = .{
        .heat_capacity_megajoules_per_k = &capacity,
        .minimum_heat_capacity_megajoules_per_k = &zero,
        .bulk_density_megagrams_per_m3 = &density,
        .liquid_water_fraction = &liquid_fraction,
        .ice_fraction = &zero,
        .air_fraction = &air_fraction,
        .fraction_of_pore_volume_air_filled = &air_fraction,
        .solid_conductivity_numerator_m_megajoules_per_h_k = &numerator,
        .solid_conductivity_denominator = &denominator,
        .is_top_soil_layer = &top,
        .top_snow_heat_capacity_megajoules_per_k = &zero,
        .maximum_negligible_snow_heat_capacity_megajoules_per_k = &zero,
        .snow_storage_heat_flux_megajoules = &zero,
        .cell_heat_source_megajoules = &cooling_megajoules,
        .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
        .turbulence = .{
            .water_fraction_threshold = 1,
            .air_fraction_threshold = 1,
            .water_rayleigh_coefficient = 0,
            .air_rayleigh_coefficient = 0,
            .water_nusselt_denominator = 1,
            .air_nusselt_denominator = 1,
        },
        .enthalpy_coupling = coupling,
    };
    const initial_parameters =
        try enthalpyParameters(properties, coupling, 0);
    const initial_state = try enthalpy.stateAtTemperature(
        initial_parameters,
        grid.soil_temperature_k[0],
    );
    var no_face_heat_flux: [0]f64 = .{};
    const no_water_flux: [0]f64 = .{};
    const result = try solve(
        std.testing.allocator,
        &grid,
        &.{},
        properties,
        .{
            .liquid_water_m3 = &no_water_flux,
            .vapor_m3 = &no_water_flux,
            .macropore_water_m3 = &no_water_flux,
        },
        &no_face_heat_flux,
        .{ .max_iterations = 80 },
    );
    const final_state = try enthalpy.stateAtTemperature(
        initial_parameters,
        grid.soil_temperature_k[0],
    );
    try std.testing.expect(result.iterations < 80);
    try std.testing.expect(grid.matrix_ice_water_m3[0] > 0);
    try std.testing.expect(grid.macropore_ice_water_m3[0] > 0);
    try std.testing.expectApproxEqAbs(
        initial_state.enthalpy_megajoules + cooling_megajoules[0],
        final_state.enthalpy_megajoules,
        1.0e-8,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.4),
        grid.matrix_liquid_water_m3[0] +
            grid.matrix_ice_water_m3[0],
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.05),
        grid.macropore_liquid_water_m3[0] +
            grid.macropore_ice_water_m3[0],
        1.0e-12,
    );
}

test "runtime Dirichlet thermal face uses cell distance and physical step" {
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
            .max_nonlinear_iterations = 40,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.soil_temperature_k[0] = 280;
    const capacity = [_]f64{2};
    const zero = [_]f64{0};
    const density = [_]f64{1};
    const liquid = [_]f64{0.2};
    const air = [_]f64{0.3};
    const numerator = [_]f64{0.01};
    const denominator = [_]f64{1};
    const top = [_]bool{true};
    const boundary_cell = [_]usize{0};
    const boundary_temperature = [_]f64{300};
    const boundary_distance = [_]f64{0.5};
    const boundary_area = [_]f64{1};
    const properties: Properties = .{
        .heat_capacity_megajoules_per_k = &capacity,
        .minimum_heat_capacity_megajoules_per_k = &zero,
        .bulk_density_megagrams_per_m3 = &density,
        .liquid_water_fraction = &liquid,
        .ice_fraction = &zero,
        .air_fraction = &air,
        .fraction_of_pore_volume_air_filled = &air,
        .solid_conductivity_numerator_m_megajoules_per_h_k = &numerator,
        .solid_conductivity_denominator = &denominator,
        .is_top_soil_layer = &top,
        .top_snow_heat_capacity_megajoules_per_k = &zero,
        .maximum_negligible_snow_heat_capacity_megajoules_per_k = &zero,
        .snow_storage_heat_flux_megajoules = &zero,
        .cell_heat_source_megajoules = &zero,
        .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
        .turbulence = .{
            .water_fraction_threshold = 1,
            .air_fraction_threshold = 1,
            .water_rayleigh_coefficient = 0,
            .air_rayleigh_coefficient = 0,
            .water_nusselt_denominator = 1,
            .air_nusselt_denominator = 1,
        },
        .time_step_hours = 0.25,
        .dirichlet_thermal_boundaries = .{
            .cell_index = &boundary_cell,
            .temperature_k = &boundary_temperature,
            .distance_from_cell_center_m = &boundary_distance,
            .face_area_m2 = &boundary_area,
        },
    };
    const empty_flux: [0]f64 = .{};
    var heat_flux: [0]f64 = .{};
    const conductivity = try heat.calculateCellConductivity(
        cellConductivityInputs(
            properties,
            0,
            20,
            .{
                .matrix_liquid_m3 = &heat_flux,
                .matrix_ice_m3 = &heat_flux,
                .macropore_liquid_m3 = &heat_flux,
                .macropore_ice_m3 = &heat_flux,
                .macropore_enabled = false,
            },
        ),
        properties.turbulence,
    );
    const coefficient = conductivity * boundary_area[0] *
        properties.time_step_hours /
        (boundary_distance[0] * capacity[0]);
    const expected_temperature_k =
        (280 + coefficient * boundary_temperature[0]) /
        (1 + coefficient);
    var workspace = try Workspace.init(
        std.testing.allocator,
        1,
        0,
        40,
    );
    defer workspace.deinit();
    const first_result = try solveWithWorkspace(
        &workspace,
        &grid,
        &.{},
        properties,
        .{
            .liquid_water_m3 = &empty_flux,
            .vapor_m3 = &empty_flux,
            .macropore_water_m3 = &empty_flux,
        },
        &heat_flux,
        .{ .max_iterations = 40 },
    );
    try std.testing.expectApproxEqAbs(
        expected_temperature_k,
        grid.soil_temperature_k[0],
        1.0e-9,
    );
    try std.testing.expectApproxEqAbs(
        capacity[0] * (expected_temperature_k - 280),
        first_result.boundary_heat_input_megajoules,
        1.0e-9,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        first_result.boundary_heat_output_megajoules,
    );
    const first_temperature_k = grid.soil_temperature_k[0];
    _ = try solveWithWorkspace(
        &workspace,
        &grid,
        &.{},
        properties,
        .{
            .liquid_water_m3 = &empty_flux,
            .vapor_m3 = &empty_flux,
            .macropore_water_m3 = &empty_flux,
        },
        &heat_flux,
        .{ .max_iterations = 40 },
    );
    try std.testing.expect(grid.soil_temperature_k[0] >
        first_temperature_k);
}

test "unfrozen pressure head prices total water on the unmodified retention curve" {
    // HEAT-001 regression. `1c97ba0` priced the LIQUID content against a
    // retention curve whose saturated content had been reduced by the ice
    // fraction. Both halves of that are wrong, and this test pins each one
    // independently so a future change cannot reintroduce one while the other
    // stays fixed.
    const curve = retention.MualemVanGenuchtenParameters{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.5,
        .alpha_per_m = 1.6,
        .n = 1.6,
        .saturated_hydraulic_conductivity_m_per_h = 0.01,
    };
    const volume_m3 = 1.0;

    // Property 1: the head is a function of TOTAL water, so repartitioning the
    // same total between liquid and ice cannot move it. This is what makes the
    // Dall'Amico reference state independent of the answer the equilibrium is
    // solving for. The old code made the head fall as ice grew.
    //
    // The comparison is approximate rather than exact only because
    // `total*(1-f) + total*f` is not bit-identical to `total` in binary floating
    // point; the tolerance covers the TEST's own rounding, not any latitude in
    // the property. `1e-15` is about four ulp here, while the defect this pins
    // moves the head by more than `1e-3` (asserted below).
    const total_m3 = 0.30;
    const total_head_m = try unfrozenPressureHeadM(curve, total_m3, volume_m3);
    var ice_fraction: f64 = 0;
    while (ice_fraction <= 0.9) : (ice_fraction += 0.1) {
        const liquid_m3 = total_m3 * (1.0 - ice_fraction);
        const ice_m3 = total_m3 * ice_fraction;
        try std.testing.expectApproxEqRel(
            total_head_m,
            try unfrozenPressureHeadM(curve, liquid_m3 + ice_m3, volume_m3),
            1e-15,
        );
    }

    // Property 2: the head is a point ON the curve it is consumed with, so
    // reading the water content back off that same curve returns the input.
    // An ice-shrunk curve breaks exactly this round trip, which is what left
    // `soil_enthalpy_balance` evaluating the enthalpy and its own temperature
    // derivative on two different constitutive relations.
    try std.testing.expectApproxEqRel(
        total_m3 / volume_m3,
        try curve.waterContentAtPressureHead(total_head_m),
        1e-12,
    );

    // The defect's own signature, to prove property 1 is a real constraint and
    // not a tautology a degenerate curve would also satisfy: the head the old
    // code produced for a 0.2 ice fraction is distinguishable from the correct
    // head by far more than a rounding error.
    const ice_shrunk = retention.MualemVanGenuchtenParameters{
        .residual_water_content_m3_per_m3 = curve.residual_water_content_m3_per_m3,
        .saturated_water_content_m3_per_m3 = curve.saturated_water_content_m3_per_m3 - 0.2,
        .alpha_per_m = curve.alpha_per_m,
        .n = curve.n,
        .pore_connectivity = curve.pore_connectivity,
        .saturated_hydraulic_conductivity_m_per_h = curve.saturated_hydraulic_conductivity_m_per_h,
    };
    const old_code_head_m = try ice_shrunk.pressureHeadAtWaterContent(0.1);
    try std.testing.expect(@abs(old_code_head_m - total_head_m) > 1e-3);

    // Saturation is where the head is zero, and total water at or above
    // saturation clamps there rather than leaving the retention domain. This is
    // the guard `soil_phase_solver.freezeThawEquilibrium` already applies: an
    // admissibility clamp on the INPUT, not a widened domain.
    try std.testing.expectEqual(
        @as(f64, 0),
        try unfrozenPressureHeadM(curve, 0.5, volume_m3),
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try unfrozenPressureHeadM(curve, 0.7, volume_m3),
    );

    // Degenerate geometry is an error, not a silent infinity.
    try std.testing.expectError(
        error.InvalidCoupledSoilPorousMediumVolume,
        unfrozenPressureHeadM(curve, 0.3, 0),
    );
}
