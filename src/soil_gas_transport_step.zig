const std = @import("std");
const grid_module = @import("grid.zig");
const gas = @import("gas_transport.zig");
const gas_faces = @import("gas_face_assembly.zig");
const gas_solver = @import("coupled_gas_solver.zig");
const gas_failure_reporter = @import("coupled_gas_failure_reporter.zig");
const gas_failure_snapshot = @import("coupled_gas_failure_snapshot.zig");
const gas_atmosphere = @import("gas_atmosphere_exchange.zig");
const hydrology_module = @import("transport_hydrology.zig");
const geometry_module = @import("soil_face_geometry.zig");
const surface_precipitation = @import("surface_precipitation.zig");
const boundary_topology = @import("soil_boundary_topology.zig");

pub const RuntimeParameters = struct {
    reference_temperature_k: f64 = 298.15,
    temperature_exponent: f64 = 1.75,
    free_air_diffusivity_m2_per_h: [gas.species_count]f64 = .{ 4.68e-2, 7.80e-2, 6.43e-2, 5.57e-2, 5.57e-2, 6.67e-2, 5.57e-2 },
    penman_tortuosity: f64 = 0.66,
    minimum_air_filled_porosity_m3_per_m3: f64 = 1e-12,
    water_density_g_per_m3: f64 = 1.0e6,
    water_molar_mass_g_per_mol: f64 = 18,
};

pub const SurfaceBoundaryInputs = struct {
    /// WATSUB `PARG`, already including `RATG + RAGS + RAS`.
    atmospheric_conductance_m3_per_step: []const f64,
    cell_area_m2: []const f64,
    top_layer_thickness_m: []const f64,
    atmospheric_concentration_g_per_m3: [gas.species_count]f64,
};

pub const SubsurfaceBoundaryInputs = struct {
    topology: *const boundary_topology.State,
    layer_thickness_m: []const f64,
    external_concentration_g_per_m3: [gas.species_count]f64,
};

pub const FailureReportRequest = struct {
    io: std.Io,
    directory: std.Io.Dir,
    file_path: []const u8,
    options: gas_failure_reporter.Options = .{},
};

pub const AdvanceRequest = struct {
    grid: *const grid_module.GridState,
    hydrology: *const hydrology_module.State,
    soil_faces: *const hydrology_module.SoilFaces,
    geometry: *const geometry_module.State,
    matrix_bulk_volume_m3: []const f64,
    total_porosity_fraction: []const f64,
    field_capacity_fraction: []const f64,
    gas_state: *gas.State,
    solubility_parameters: gas.SurfaceSolubilityParameters,
    exchange_parameters: surface_precipitation.GasExchangeParameters,
    ammonium_band_fraction: f64,
    surface_boundary_inputs: ?SurfaceBoundaryInputs,
    subsurface_boundary_inputs: ?SubsurfaceBoundaryInputs,
    parameters: RuntimeParameters,
    solver_options: gas_solver.Options,
    failure_report: ?FailureReportRequest = null,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    air_filled_porosity_m3_per_m3: []f64,
    free_air_diffusivity_m2_per_step: []f64,
    band_water_volume_m3: []f64,
    mass_solubility_ratio: []f64,
    gas_water_exchange_rate_per_step: []f64,
    band_gas_water_exchange_rate_per_step: []f64,
    bubbling_enabled: []bool,
    surface_boundaries: []gas_atmosphere.Boundary,
    subsurface_boundaries: []gas_atmosphere.Boundary,
    atmospheric_flux_g_per_h: []f64,
    substep_atmospheric_flux_g: []f64,
    subsurface_flux_g_per_h: []f64,
    substep_subsurface_flux_g: []f64,
    accepted_face_flux_g_per_h: []f64,
    substep_face_flux_g: []f64,
    accepted_faces: []gas.Face,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.ZeroSoilGasLayerCount;
        const components = try std.math.mul(usize, layer_count, gas.species_count);
        const air = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(air);
        const diffusivity = try allocator.alloc(f64, components);
        errdefer allocator.free(diffusivity);
        const band_water = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(band_water);
        const solubility = try allocator.alloc(f64, components);
        errdefer allocator.free(solubility);
        const exchange = try allocator.alloc(f64, components);
        errdefer allocator.free(exchange);
        const band_exchange = try allocator.alloc(f64, components);
        errdefer allocator.free(band_exchange);
        const bubbling = try allocator.alloc(bool, layer_count);
        errdefer allocator.free(bubbling);
        const boundaries = try allocator.alloc(gas_atmosphere.Boundary, layer_count);
        errdefer allocator.free(boundaries);
        const subsurface_boundaries = try allocator.alloc(gas_atmosphere.Boundary, 0);
        errdefer allocator.free(subsurface_boundaries);
        const atmospheric_flux = try allocator.alloc(f64, components);
        errdefer allocator.free(atmospheric_flux);
        const substep_flux = try allocator.alloc(f64, components);
        errdefer allocator.free(substep_flux);
        const subsurface_flux = try allocator.alloc(f64, components);
        errdefer allocator.free(subsurface_flux);
        const substep_subsurface_flux = try allocator.alloc(f64, components);
        errdefer allocator.free(substep_subsurface_flux);
        const accepted_face_flux = try allocator.alloc(f64, 0);
        errdefer allocator.free(accepted_face_flux);
        const substep_face_flux = try allocator.alloc(f64, 0);
        errdefer allocator.free(substep_face_flux);
        const accepted_faces = try allocator.alloc(gas.Face, 0);
        errdefer allocator.free(accepted_faces);
        @memset(air, 0);
        @memset(diffusivity, 0);
        @memset(band_water, 0);
        @memset(solubility, 0);
        @memset(exchange, 0);
        @memset(band_exchange, 0);
        @memset(bubbling, true);
        @memset(atmospheric_flux, 0);
        @memset(substep_flux, 0);
        @memset(subsurface_flux, 0);
        @memset(substep_subsurface_flux, 0);
        return .{ .allocator = allocator, .air_filled_porosity_m3_per_m3 = air, .free_air_diffusivity_m2_per_step = diffusivity, .band_water_volume_m3 = band_water, .mass_solubility_ratio = solubility, .gas_water_exchange_rate_per_step = exchange, .band_gas_water_exchange_rate_per_step = band_exchange, .bubbling_enabled = bubbling, .surface_boundaries = boundaries, .subsurface_boundaries = subsurface_boundaries, .atmospheric_flux_g_per_h = atmospheric_flux, .substep_atmospheric_flux_g = substep_flux, .subsurface_flux_g_per_h = subsurface_flux, .substep_subsurface_flux_g = substep_subsurface_flux, .accepted_face_flux_g_per_h = accepted_face_flux, .substep_face_flux_g = substep_face_flux, .accepted_faces = accepted_faces };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.accepted_faces);
        self.allocator.free(self.substep_face_flux_g);
        self.allocator.free(self.accepted_face_flux_g_per_h);
        self.allocator.free(self.substep_subsurface_flux_g);
        self.allocator.free(self.subsurface_flux_g_per_h);
        self.allocator.free(self.substep_atmospheric_flux_g);
        self.allocator.free(self.atmospheric_flux_g_per_h);
        self.allocator.free(self.subsurface_boundaries);
        self.allocator.free(self.surface_boundaries);
        self.allocator.free(self.bubbling_enabled);
        self.allocator.free(self.band_gas_water_exchange_rate_per_step);
        self.allocator.free(self.gas_water_exchange_rate_per_step);
        self.allocator.free(self.mass_solubility_ratio);
        self.allocator.free(self.band_water_volume_m3);
        self.allocator.free(self.free_air_diffusivity_m2_per_step);
        self.allocator.free(self.air_filled_porosity_m3_per_m3);
        self.* = undefined;
    }

    pub fn advance(self: *State, request: AdvanceRequest) !gas_solver.Result {
        const grid = request.grid;
        const hydrology = request.hydrology;
        const soil_faces = request.soil_faces;
        const geometry = request.geometry;
        const matrix_bulk_volume_m3 = request.matrix_bulk_volume_m3;
        const total_porosity_fraction = request.total_porosity_fraction;
        const field_capacity_fraction = request.field_capacity_fraction;
        const gas_state = request.gas_state;
        const solubility_parameters = request.solubility_parameters;
        const exchange_parameters = request.exchange_parameters;
        const ammonium_band_fraction = request.ammonium_band_fraction;
        const surface_boundary_inputs = request.surface_boundary_inputs;
        const subsurface_boundary_inputs = request.subsurface_boundary_inputs;
        const parameters = request.parameters;
        const solver_options = request.solver_options;
        try validate(self, grid, hydrology, matrix_bulk_volume_m3, total_porosity_fraction, field_capacity_fraction, gas_state, ammonium_band_fraction, parameters);
        // NPG enlarges only solver_options.max_iterations. Every coefficient
        // below represents the complete hour; there is no fixed gas substep.
        @memset(self.atmospheric_flux_g_per_h, 0);
        @memset(self.subsurface_flux_g_per_h, 0);
        var whole_hour_exchange = exchange_parameters;
        whole_hour_exchange.iteration_fraction = 1;
        for (0..grid.layer_count) |layer| {
            const air_fraction = std.math.clamp(grid.matrix_air_volume_m3[layer] / matrix_bulk_volume_m3[layer], 0, total_porosity_fraction[layer]);
            self.air_filled_porosity_m3_per_m3[layer] = air_fraction;
            const temperature_factor = std.math.pow(f64, grid.soil_temperature_k[layer] / parameters.reference_temperature_k, parameters.temperature_exponent);
            const solubility = try gas.surfaceSolubilityWaterToAir(grid.soil_temperature_k[layer], solubility_parameters);
            const exchange = try surface_precipitation.litterGasExchange(grid.matrix_pore_capacity_m3[layer], grid.matrix_ice_water_m3[layer], grid.matrix_liquid_water_m3[layer], grid.matrix_air_volume_m3[layer], field_capacity_fraction[layer] * matrix_bulk_volume_m3[layer], whole_hour_exchange);
            self.band_water_volume_m3[layer] = ammonium_band_fraction * grid.matrix_liquid_water_m3[layer];
            gas_state.air_volume_m3[layer] = grid.matrix_air_volume_m3[layer];
            gas_state.temperature_k[layer] = grid.soil_temperature_k[layer];
            gas_state.water_vapor_mol[layer] = grid.water_vapor_volume_m3[layer] * parameters.water_density_g_per_m3 / parameters.water_molar_mass_g_per_mol;
            for (0..gas.species_count) |species| {
                const component = layer * gas.species_count + species;
                self.free_air_diffusivity_m2_per_step[component] = parameters.free_air_diffusivity_m2_per_h[species] * temperature_factor;
                self.mass_solubility_ratio[component] = solubility[species];
                self.gas_water_exchange_rate_per_step[component] = exchange.air_water_rate_per_step;
                self.band_gas_water_exchange_rate_per_step[component] = if (species == @intFromEnum(gas.Species.ammonia)) exchange.air_water_rate_per_step else 0;
            }
        }
        var face_set = try gas_faces.buildMapped(self.allocator, soil_faces, geometry, self.air_filled_porosity_m3_per_m3, total_porosity_fraction, self.free_air_diffusivity_m2_per_step, parameters.penman_tortuosity, parameters.minimum_air_filled_porosity_m3_per_m3);
        defer face_set.deinit();
        const face_component_count = try std.math.mul(
            usize,
            face_set.faces.len,
            gas.species_count,
        );
        self.accepted_face_flux_g_per_h = try self.allocator.realloc(
            self.accepted_face_flux_g_per_h,
            face_component_count,
        );
        self.substep_face_flux_g = try self.allocator.realloc(
            self.substep_face_flux_g,
            face_component_count,
        );
        self.accepted_faces = try self.allocator.realloc(
            self.accepted_faces,
            face_set.faces.len,
        );
        const atmospheric_boundaries = if (surface_boundary_inputs) |inputs|
            try self.refreshSurfaceBoundaries(grid, total_porosity_fraction, inputs, parameters)
        else
            self.surface_boundaries[0..0];
        const subsurface_boundaries = if (subsurface_boundary_inputs) |inputs|
            try self.refreshSubsurfaceBoundaries(grid, matrix_bulk_volume_m3, total_porosity_fraction, inputs, parameters)
        else
            self.subsurface_boundaries[0..0];
        const solve_inputs: gas_solver.Inputs = .{
            .faces = face_set.faces,
            .face_conductance_m3_per_step = face_set.conductance_m3_per_step,
            .atmospheric_boundaries = atmospheric_boundaries,
            .subsurface_boundaries = subsurface_boundaries,
            .water_volume_m3 = grid.matrix_liquid_water_m3,
            .band_water_volume_m3 = self.band_water_volume_m3,
            .mass_solubility_ratio = self.mass_solubility_ratio,
            .gas_water_exchange_rate_per_step = self.gas_water_exchange_rate_per_step,
            .band_gas_water_exchange_rate_per_step = self.band_gas_water_exchange_rate_per_step,
            .bubbling_enabled = self.bubbling_enabled,
            .atmospheric_flux_g_by_component = self.substep_atmospheric_flux_g,
            .subsurface_flux_g_by_component = self.substep_subsurface_flux_g,
            .face_flux_g_by_component = self.substep_face_flux_g,
        };
        var hourly_options = solver_options;
        hourly_options.transport_iteration_fraction = 1;
        const result = gas_solver.solve(self.allocator, gas_state, solve_inputs, hourly_options) catch |err| {
            @memset(self.atmospheric_flux_g_per_h, 0);
            @memset(self.subsurface_flux_g_per_h, 0);
            @memset(self.accepted_face_flux_g_per_h, 0);
            if (request.failure_report) |report| return gas_failure_reporter.reportPreservingSolverError(
                self.allocator,
                report.io,
                report.directory,
                report.file_path,
                gas_state,
                solve_inputs,
                hourly_options,
                report.options,
                err,
            );
            return err;
        };
        @memcpy(self.atmospheric_flux_g_per_h, self.substep_atmospheric_flux_g);
        @memcpy(self.subsurface_flux_g_per_h, self.substep_subsurface_flux_g);
        @memcpy(self.accepted_face_flux_g_per_h, self.substep_face_flux_g);
        @memcpy(self.accepted_faces, face_set.faces);
        return result;
    }

    fn refreshSubsurfaceBoundaries(self: *State, grid: *const grid_module.GridState, matrix_bulk_volume_m3: []const f64, total_porosity_fraction: []const f64, inputs: SubsurfaceBoundaryInputs, parameters: RuntimeParameters) ![]const gas_atmosphere.Boundary {
        if (inputs.layer_thickness_m.len != grid.layer_count or matrix_bulk_volume_m3.len != grid.layer_count or total_porosity_fraction.len != grid.layer_count) return error.SoilGasSubsurfaceBoundaryDimensionMismatch;
        self.subsurface_boundaries = try self.allocator.realloc(self.subsurface_boundaries, inputs.topology.faces.len);
        for (inputs.topology.faces, 0..) |face, boundary_index| {
            const layer = face.layer_index;
            if (layer >= grid.layer_count) return error.SoilGasSubsurfaceBoundaryDimensionMismatch;
            const exchange_fraction = std.math.clamp(face.natural_exchange_fraction + face.artificial_exchange_fraction, 0, 1);
            const path_m = if (face.is_lower_boundary) inputs.layer_thickness_m[layer] else face.directional_layer_width_m;
            const face_area_m2 = matrix_bulk_volume_m3[layer] / path_m;
            const air_fraction = self.air_filled_porosity_m3_per_m3[layer];
            const porosity = total_porosity_fraction[layer];
            if (!std.math.isFinite(path_m) or path_m <= 0 or !std.math.isFinite(face_area_m2) or face_area_m2 <= 0 or porosity <= 0) return error.InvalidSoilGasSubsurfaceGeometry;
            const diffusion_geometry_m = exchange_fraction * air_fraction * parameters.penman_tortuosity * air_fraction / porosity * face_area_m2 / path_m;
            var interior: [gas.species_count]f64 = undefined;
            for (&interior, 0..) |*conductance, species| conductance.* = diffusion_geometry_m * self.free_air_diffusivity_m2_per_step[layer * gas.species_count + species];
            self.subsurface_boundaries[boundary_index] = .{
                .cell_index = layer,
                // The topology exchange fraction and path are already folded
                // into the soil-side conductance; a very large outer
                // conductance makes the series limit equal to that value.
                .aerodynamic_conductance_m3_per_step = std.math.floatMax(f64),
                .interior_conductance_m3_per_step = interior,
                .atmospheric_concentration_g_per_m3 = inputs.external_concentration_g_per_m3,
                .pressure_exchange_fraction = exchange_fraction,
            };
        }
        return self.subsurface_boundaries;
    }

    fn refreshSurfaceBoundaries(self: *State, grid: *const grid_module.GridState, total_porosity_fraction: []const f64, inputs: SurfaceBoundaryInputs, parameters: RuntimeParameters) ![]const gas_atmosphere.Boundary {
        const cell_count = grid.cell_count;
        if (inputs.atmospheric_conductance_m3_per_step.len != cell_count or inputs.cell_area_m2.len != cell_count or inputs.top_layer_thickness_m.len != grid.layer_count or self.surface_boundaries.len < cell_count) return error.SoilGasSurfaceBoundaryDimensionMismatch;
        for (0..cell_count) |cell| {
            const top = cell * grid.soil_layer_capacity;
            const air_fraction = self.air_filled_porosity_m3_per_m3[top];
            const porosity = total_porosity_fraction[top];
            const area = inputs.cell_area_m2[cell];
            const thickness = inputs.top_layer_thickness_m[top];
            if (!std.math.isFinite(area) or area <= 0 or !std.math.isFinite(thickness) or thickness <= 0) return error.InvalidSoilGasSurfaceGeometry;
            const diffusion_geometry_m = air_fraction * parameters.penman_tortuosity * air_fraction / porosity * area / thickness;
            var interior: [gas.species_count]f64 = undefined;
            for (&interior, 0..) |*conductance, species| conductance.* = diffusion_geometry_m * self.free_air_diffusivity_m2_per_step[top * gas.species_count + species];
            self.surface_boundaries[cell] = .{
                .cell_index = top,
                .aerodynamic_conductance_m3_per_step = inputs.atmospheric_conductance_m3_per_step[cell],
                .interior_conductance_m3_per_step = interior,
                .atmospheric_concentration_g_per_m3 = inputs.atmospheric_concentration_g_per_m3,
            };
        }
        return self.surface_boundaries[0..cell_count];
    }
};

fn validate(self: *const State, grid: *const grid_module.GridState, hydrology: *const hydrology_module.State, bulk: []const f64, porosity: []const f64, field_capacity: []const f64, gas_state: *const gas.State, band_fraction: f64, p: RuntimeParameters) !void {
    _ = hydrology;
    if (bulk.len != grid.layer_count or porosity.len != grid.layer_count or field_capacity.len != grid.layer_count or gas_state.cell_count != grid.layer_count or self.air_filled_porosity_m3_per_m3.len != grid.layer_count) return error.SoilGasDimensionMismatch;
    if (!std.math.isFinite(band_fraction) or band_fraction < 0 or band_fraction > 1) return error.InvalidSoilGasBandFraction;
    inline for (.{ p.reference_temperature_k, p.temperature_exponent, p.penman_tortuosity, p.minimum_air_filled_porosity_m3_per_m3, p.water_density_g_per_m3, p.water_molar_mass_g_per_mol }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilGasRuntimeParameter;
    if (p.reference_temperature_k == 0 or p.water_density_g_per_m3 == 0 or p.water_molar_mass_g_per_mol == 0) return error.InvalidSoilGasRuntimeParameter;
    for (p.free_air_diffusivity_m2_per_h) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilGasRuntimeParameter;
    for (bulk, porosity, field_capacity, grid.soil_temperature_k) |volume, pore, capacity, temperature| if (!std.math.isFinite(volume) or volume <= 0 or !std.math.isFinite(pore) or pore <= 0 or !std.math.isFinite(capacity) or capacity < 0 or !std.math.isFinite(temperature) or temperature <= 0) return error.InvalidSoilGasLayerState;
}

test "mapped soil gas step conserves internal gaseous inventory" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 2, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.active_soil_layer_count, 1);
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    @memset(grid.matrix_air_volume_m3, 0.25);
    @memset(grid.air_volume_m3, 0.25);
    @memset(grid.matrix_liquid_water_m3, 0.25);
    @memset(grid.soil_temperature_k, 298.15);
    var snow = try @import("snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    var geometry = try geometry_module.State.initMapped(std.testing.allocator, &grid, &faces, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 });
    defer geometry.deinit();
    var gas_state = try gas.State.init(std.testing.allocator, 2);
    defer gas_state.deinit();
    gas_state.gaseous_mass_g[0] = 2;
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    const solubility = gas.SurfaceSolubilityParameters{ .reference_water_to_air = [_]f64{1} ** gas.species_count, .log_intercept = [_]f64{0} ** gas.species_count, .temperature_coefficient_per_c = [_]f64{0} ** gas.species_count };
    _ = try state.advance(.{
        .grid = &grid,
        .hydrology = &hydrology,
        .soil_faces = &faces,
        .geometry = &geometry,
        .matrix_bulk_volume_m3 = &.{ 1, 1 },
        .total_porosity_fraction = &.{ 0.5, 0.5 },
        .field_capacity_fraction = &.{ 0.3, 0.3 },
        .gas_state = &gas_state,
        .solubility_parameters = solubility,
        .exchange_parameters = .{ .reference_time_h = 1, .wet_exponent = 12, .dry_exponent = 12, .transition_water_fraction = 0.5, .iteration_fraction = 0, .aqueous_tortuosity_coefficient = 0.7 },
        .ammonium_band_fraction = 0,
        .surface_boundary_inputs = null,
        .subsurface_boundary_inputs = null,
        .parameters = .{},
        .solver_options = .{ .max_iterations = 80 },
    });
    try std.testing.expectEqual(
        faces.micropore_faces.len * gas.species_count,
        state.accepted_face_flux_g_per_h.len,
    );
    try std.testing.expectEqual(
        faces.micropore_faces.len,
        state.accepted_faces.len,
    );
    try std.testing.expectEqual(@as(usize, 0), state.accepted_faces[0].first_cell);
    try std.testing.expectEqual(@as(usize, 1), state.accepted_faces[0].second_cell);
    try std.testing.expect(state.accepted_face_flux_g_per_h[0] > 0);
    try std.testing.expectApproxEqAbs(@as(f64, 2), gas_state.gaseous_mass_g[0] + gas_state.gaseous_mass_g[gas.species_count] + gas_state.dissolved_mass_g[0] + gas_state.dissolved_mass_g[gas.species_count], 1e-12);
    try std.testing.expect(gas_state.gaseous_mass_g[gas.species_count] > 0);
    try std.testing.expect(gas_state.dissolved_mass_g[0] + gas_state.dissolved_mass_g[gas.species_count] > 0);
}

test "surface PARG is coupled in series with runtime top-layer diffusivity" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.active_soil_layer_count, 1);
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    @memset(grid.matrix_air_volume_m3, 0.25);
    @memset(grid.air_volume_m3, 0.25);
    @memset(grid.matrix_liquid_water_m3, 0.25);
    @memset(grid.soil_temperature_k, 298.15);
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer hydrology.deinit();
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    var geometry = try geometry_module.State.initMapped(std.testing.allocator, &grid, &faces, &.{1}, &.{1}, &.{1});
    defer geometry.deinit();
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const solubility = gas.SurfaceSolubilityParameters{ .reference_water_to_air = [_]f64{1} ** gas.species_count, .log_intercept = [_]f64{0} ** gas.species_count, .temperature_coefficient_per_c = [_]f64{0} ** gas.species_count };
    var atmosphere = [_]f64{0} ** gas.species_count;
    atmosphere[0] = 1;
    _ = try state.advance(.{
        .grid = &grid,
        .hydrology = &hydrology,
        .soil_faces = &faces,
        .geometry = &geometry,
        .matrix_bulk_volume_m3 = &.{1},
        .total_porosity_fraction = &.{0.5},
        .field_capacity_fraction = &.{0.3},
        .gas_state = &gas_state,
        .solubility_parameters = solubility,
        .exchange_parameters = .{ .reference_time_h = 1, .wet_exponent = 12, .dry_exponent = 12, .transition_water_fraction = 0.5, .iteration_fraction = 0, .aqueous_tortuosity_coefficient = 0.7 },
        .ammonium_band_fraction = 0,
        .surface_boundary_inputs = .{
            .atmospheric_conductance_m3_per_step = &.{0.1},
            .cell_area_m2 = &.{2},
            .top_layer_thickness_m = &.{1},
            .atmospheric_concentration_g_per_m3 = atmosphere,
        },
        .subsurface_boundary_inputs = null,
        .parameters = .{},
        .solver_options = .{ .max_iterations = 80 },
    });
    const top_layer_conductance = 0.25 * 0.66 * 0.25 / 0.5 * 2 * 4.68e-2;
    try std.testing.expectApproxEqAbs(top_layer_conductance, state.surface_boundaries[0].interior_conductance_m3_per_step[0], 1e-15);
    try std.testing.expect(gas_state.gaseous_mass_g[0] > 0);
    try std.testing.expectApproxEqAbs(gas_state.gaseous_mass_g[0] + gas_state.dissolved_mass_g[0] + gas_state.band_dissolved_mass_g[0], state.atmospheric_flux_g_per_h[0], 1e-12);
}

test "solver catch path preserves failure and publishes a readable snapshot" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const cfg = try @import("config.zig").SimulationConfig.init(
        .{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.active_soil_layer_count, 1);
    @memset(grid.matrix_pore_capacity_m3, 1);
    @memset(grid.matrix_air_volume_m3, 1);
    @memset(grid.air_volume_m3, 1);
    @memset(grid.matrix_liquid_water_m3, 1);
    @memset(grid.soil_temperature_k, 298.15);
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer hydrology.deinit();
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    var geometry = try geometry_module.State.initMapped(
        std.testing.allocator,
        &grid,
        &faces,
        &.{1},
        &.{1},
        &.{1},
    );
    defer geometry.deinit();
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.gaseous_mass_g[0] = 2;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const solubility: gas.SurfaceSolubilityParameters = .{
        .reference_water_to_air = [_]f64{1} ** gas.species_count,
        .log_intercept = [_]f64{0} ** gas.species_count,
        .temperature_coefficient_per_c = [_]f64{0} ** gas.species_count,
    };
    try std.testing.expectError(
        error.CoupledGasSolverDidNotConverge,
        state.advance(.{
            .grid = &grid,
            .hydrology = &hydrology,
            .soil_faces = &faces,
            .geometry = &geometry,
            .matrix_bulk_volume_m3 = &.{1},
            .total_porosity_fraction = &.{1},
            .field_capacity_fraction = &.{1},
            .gas_state = &gas_state,
            .solubility_parameters = solubility,
            .exchange_parameters = .{
                .reference_time_h = 1,
                .wet_exponent = 12,
                .dry_exponent = 12,
                .transition_water_fraction = 0.5,
                .iteration_fraction = 0,
                .aqueous_tortuosity_coefficient = 0.7,
            },
            .ammonium_band_fraction = 0,
            .surface_boundary_inputs = null,
            .subsurface_boundary_inputs = null,
            .parameters = .{},
            .solver_options = .{
                .absolute_tolerance_g = 1e-20,
                .relative_tolerance = 1e-20,
                .max_iterations = 1,
            },
            .failure_report = .{
                .io = std.testing.io,
                .directory = temporary.dir,
                .file_path = "cell0-hour1-gas-failure.bin",
                .options = .{ .write_buffer_bytes = 37 },
            },
        }),
    );
    const bytes = try temporary.dir.readFileAlloc(
        std.testing.io,
        "cell0-hour1-gas-failure.bin",
        std.testing.allocator,
        .limited(1024 * 1024),
    );
    defer std.testing.allocator.free(bytes);
    var reader: std.Io.Reader = .fixed(bytes);
    var replay_case = try gas_failure_snapshot.read(
        std.testing.allocator,
        &reader,
        .{},
    );
    defer replay_case.deinit();
    try std.testing.expectEqual(@as(usize, 1), replay_case.state.cell_count);
    try std.testing.expectEqual(@as(f64, 2), replay_case.state.gaseous_mass_g[0]);
    try std.testing.expectEqual(@as(u16, 1), replay_case.options.max_iterations);
}

test "runtime perimeter face builds source geometry and exchange fraction" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.air_filled_porosity_m3_per_m3[0] = 0.25;
    @memset(state.free_air_diffusivity_m2_per_step, 0.04);
    var face = boundary_topology.Face{
        .cell_index = 0,
        .layer_index = 0,
        .direction = .east,
        .direction_sign = -1,
        .directional_layer_width_m = 2,
        .slope_sine = 0,
        .natural_water_table_distance_m = 1,
        .natural_exchange_fraction = 0.4,
        .artificial_water_table_distance_m = 0,
        .artificial_exchange_fraction = 0,
        .surface_runoff_fraction = 0,
        .is_lower_boundary = false,
    };
    var water_table_mode = [_]u8{0};
    var zero_slope = [_]f64{0};
    const topology = boundary_topology.State{
        .allocator = std.testing.allocator,
        .faces = (&face)[0..1],
        .water_table_mode = &water_table_mode,
        .natural_water_table_depth_m = &.{},
        .internal_water_table_depth_m = &.{},
        .active_layer_depth_m = &.{},
        .artificial_water_table_depth_m = &.{},
        .natural_water_table_surface_slope = &zero_slope,
        .artificial_water_table_surface_slope = &zero_slope,
    };
    const boundaries = try state.refreshSubsurfaceBoundaries(&grid, &.{8}, &.{0.5}, .{
        .topology = &topology,
        .layer_thickness_m = &.{1},
        .external_concentration_g_per_m3 = [_]f64{0} ** gas.species_count,
    }, .{});
    // AREA = bulk/path = 4 m2; conductance geometry is
    // exchange * air * 0.66 * air / porosity * AREA / path.
    const expected: f64 = 0.4 * 0.25 * 0.66 * 0.25 / 0.5 * 4.0 / 2.0 * 0.04;
    try std.testing.expectApproxEqAbs(expected, boundaries[0].interior_conductance_m3_per_step[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), boundaries[0].pressure_exchange_fraction, 1e-15);
}
