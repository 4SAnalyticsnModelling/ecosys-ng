const std = @import("std");
const grid_module = @import("grid.zig");
const hydrology_module = @import("transport_hydrology.zig");
const geometry_module = @import("soil_face_geometry.zig");
const species_module = @import("solute_transport_species.zig");

const Species = species_module.AqueousSpecies;

pub const RuntimeParameters = struct {
    reference_temperature_k: f64 = 298.15,
    temperature_exponent: f64 = 6,
    phosphate_diffusivity_m2_per_h: f64 = 3.0e-6,
    other_ion_diffusivity_m2_per_h: f64 = 5.0e-6,
    micropore_tortuosity_coefficient: f64 = 0.7,
    macropore_tortuosity_coefficient: f64 = 2.8,
    dispersivity_coefficient: f64 = 0.20,
    dispersivity_distance_exponent: f64 = 1.07,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    micropore_conductance_m3_per_step: []f64,
    macropore_conductance_m3_per_step: []f64,
    micropore_mobility_fraction: []f64,
    macropore_mobility_fraction: []f64,
    boundary_mobility_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, face_count: usize) !State {
        const count = try std.math.mul(usize, face_count, Species.count);
        const micro = try allocator.alloc(f64, count);
        errdefer allocator.free(micro);
        const macro = try allocator.alloc(f64, count);
        errdefer allocator.free(macro);
        const micro_mobile = try allocator.alloc(f64, count);
        errdefer allocator.free(micro_mobile);
        const macro_mobile = try allocator.alloc(f64, count);
        errdefer allocator.free(macro_mobile);
        const boundary_mobile = try allocator.alloc(f64, Species.count);
        errdefer allocator.free(boundary_mobile);
        @memset(micro, 0);
        @memset(macro, 0);
        @memset(micro_mobile, 1);
        @memset(macro_mobile, 1);
        @memset(boundary_mobile, 1);
        return .{ .allocator = allocator, .micropore_conductance_m3_per_step = micro, .macropore_conductance_m3_per_step = macro, .micropore_mobility_fraction = micro_mobile, .macropore_mobility_fraction = macro_mobile, .boundary_mobility_fraction = boundary_mobile };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.boundary_mobility_fraction);
        self.allocator.free(self.macropore_mobility_fraction);
        self.allocator.free(self.micropore_mobility_fraction);
        self.allocator.free(self.macropore_conductance_m3_per_step);
        self.allocator.free(self.micropore_conductance_m3_per_step);
        self.* = undefined;
    }

    /// Reconstructs HOUR1/WATSUB/STARTS/TRNSFRS face coefficients from the
    /// current runtime water state. `step_hours` is physical time, not a fixed
    /// sub-hour loop fraction; nonlinear iterations converge this one step.
    pub fn refresh(self: *State, grid: *const grid_module.GridState, faces: *const hydrology_module.SoilFaces, geometry: *const geometry_module.State, matrix_bulk_volume_m3: []const f64, phosphate_non_band_fraction: f64, phosphate_band_fraction: f64, step_hours: f64, parameters: RuntimeParameters) !void {
        try validate(grid, faces, geometry, matrix_bulk_volume_m3, phosphate_non_band_fraction, phosphate_band_fraction, step_hours, parameters, self.micropore_conductance_m3_per_step.len);
        inline for (@typeInfo(Species).@"enum".fields) |field| self.boundary_mobility_fraction[field.value] = mobility(@enumFromInt(field.value), phosphate_non_band_fraction, phosphate_band_fraction);
        for (faces.micropore_faces, 0..) |face, face_index| {
            const source = face.first_cell;
            const destination = face.second_cell;
            const path_sum_m = geometry.source_path_length_m[face_index] + geometry.destination_path_length_m[face_index];
            const mean_distance_m = 0.5 * path_sum_m;
            const area_m2 = geometry.face_area_m2[face_index];
            const source_macro_fraction = macroporeFraction(grid, source);
            const destination_macro_fraction = macroporeFraction(grid, destination);
            const source_matrix_theta = std.math.clamp(grid.matrix_liquid_water_m3[source] / matrix_bulk_volume_m3[source], 0, 1);
            const destination_matrix_theta = std.math.clamp(grid.matrix_liquid_water_m3[destination] / matrix_bulk_volume_m3[destination], 0, 1);
            const source_micro_tortuosity = parameters.micropore_tortuosity_coefficient * source_matrix_theta * source_matrix_theta * (1 - source_macro_fraction);
            const destination_micro_tortuosity = parameters.micropore_tortuosity_coefficient * destination_matrix_theta * destination_matrix_theta * (1 - destination_macro_fraction);
            const source_macro_theta = if (grid.macropore_pore_capacity_m3[source] > 0) std.math.clamp(grid.macropore_liquid_water_m3[source] / grid.macropore_pore_capacity_m3[source], 0, 1) else 0;
            const destination_macro_theta = if (grid.macropore_pore_capacity_m3[destination] > 0) std.math.clamp(grid.macropore_liquid_water_m3[destination] / grid.macropore_pore_capacity_m3[destination], 0, 1) else 0;
            const source_macro_tortuosity = @min(1.0, parameters.macropore_tortuosity_coefficient * source_macro_theta * source_macro_theta * source_macro_theta) * source_macro_fraction;
            const destination_macro_tortuosity = @min(1.0, parameters.macropore_tortuosity_coefficient * destination_macro_theta * destination_macro_theta * destination_macro_theta) * destination_macro_fraction;
            const micro_tortuosity_per_m = (source_micro_tortuosity + destination_micro_tortuosity) / path_sum_m;
            const macro_tortuosity_per_m = (source_macro_tortuosity + destination_macro_tortuosity) / path_sum_m;
            const water_velocity_m_per_step = @abs(faces.micropore_water_flux_m3_per_step[face_index]) / area_m2;
            const dispersivity_m2_per_step = parameters.dispersivity_coefficient * std.math.pow(f64, mean_distance_m, parameters.dispersivity_distance_exponent) * step_hours * @min(step_hours, water_velocity_m_per_step);
            const temperature_factor = std.math.pow(f64, grid.soil_temperature_k[destination] / parameters.reference_temperature_k, parameters.temperature_exponent);
            inline for (@typeInfo(Species).@"enum".fields) |field| {
                const species: Species = @enumFromInt(field.value);
                const component = face_index * Species.count + field.value;
                const reference_diffusivity = if (species_module.diffusivityClass(species) == .phosphate) parameters.phosphate_diffusivity_m2_per_h else parameters.other_ion_diffusivity_m2_per_h;
                const diffusivity_m2_per_step = reference_diffusivity * temperature_factor * step_hours;
                self.micropore_conductance_m3_per_step[component] = (diffusivity_m2_per_step * micro_tortuosity_per_m + dispersivity_m2_per_step) * area_m2;
                self.macropore_conductance_m3_per_step[component] = diffusivity_m2_per_step * macro_tortuosity_per_m * area_m2;
                self.micropore_mobility_fraction[component] = mobility(species, phosphate_non_band_fraction, phosphate_band_fraction);
                self.macropore_mobility_fraction[component] = self.micropore_mobility_fraction[component];
            }
        }
    }
};

fn macroporeFraction(grid: *const grid_module.GridState, cell: usize) f64 {
    const total = grid.matrix_pore_capacity_m3[cell] + grid.macropore_pore_capacity_m3[cell];
    return if (total > 0) std.math.clamp(grid.macropore_pore_capacity_m3[cell] / total, 0, 1) else 0;
}

fn mobility(species: Species, non_band: f64, band: f64) f64 {
    return switch (species) {
        .non_band_phosphate, .non_band_phosphoric_acid, .non_band_iron_hpo4, .non_band_iron_h2po4, .non_band_calcium_phosphate, .non_band_calcium_hpo4, .non_band_calcium_h2po4, .non_band_magnesium_hpo4 => non_band,
        .band_phosphate, .band_phosphoric_acid, .band_iron_hpo4, .band_iron_h2po4, .band_calcium_phosphate, .band_calcium_hpo4, .band_calcium_h2po4, .band_magnesium_hpo4 => band,
        else => 1,
    };
}

fn validate(grid: *const grid_module.GridState, faces: *const hydrology_module.SoilFaces, geometry: *const geometry_module.State, bulk: []const f64, non_band: f64, band: f64, step: f64, p: RuntimeParameters, output_len: usize) !void {
    if (bulk.len != grid.layer_count or faces.micropore_faces.len != geometry.face_area_m2.len or output_len != faces.micropore_faces.len * Species.count) return error.SoilSoluteFaceDimensionMismatch;
    inline for (.{ non_band, band }) |value| if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSoilSoluteMobilityFraction;
    if (!std.math.isFinite(non_band + band) or @abs(non_band + band - 1) > 1e-10) return error.InvalidSoilSoluteMobilityFraction;
    inline for (.{ step, p.reference_temperature_k, p.phosphate_diffusivity_m2_per_h, p.other_ion_diffusivity_m2_per_h, p.micropore_tortuosity_coefficient, p.macropore_tortuosity_coefficient, p.dispersivity_coefficient, p.dispersivity_distance_exponent }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilSoluteRuntimeParameter;
    if (step == 0 or p.reference_temperature_k == 0) return error.InvalidSoilSoluteRuntimeParameter;
    for (bulk, grid.soil_temperature_k) |volume, temperature| if (!std.math.isFinite(volume) or volume <= 0 or !std.math.isFinite(temperature) or temperature <= 0) return error.InvalidSoilSoluteLayerState;
}

test "face assembler reproduces TRNSFRS micropore and macropore forms" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 2, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    @memset(grid.active_soil_layer_count, 1);
    @memset(grid.matrix_liquid_water_m3, 0.5);
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    @memset(grid.macropore_liquid_water_m3, 0.1);
    @memset(grid.macropore_pore_capacity_m3, 0.1);
    @memset(grid.soil_temperature_k, 298.15);
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    hydrology.micropore_face_flux_m3_per_step[0] = 0.01;
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    var geometry = try geometry_module.State.initMapped(std.testing.allocator, &grid, &faces, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 });
    defer geometry.deinit();
    var state = try State.init(std.testing.allocator, faces.micropore_faces.len);
    defer state.deinit();
    try state.refresh(&grid, &faces, &geometry, &.{ 1, 1 }, 0.75, 0.25, 1, .{});
    const aluminum = species_module.index(.aluminum);
    const phosphate = species_module.index(.non_band_phosphate);
    try std.testing.expect(state.micropore_conductance_m3_per_step[aluminum] > state.micropore_conductance_m3_per_step[phosphate]);
    try std.testing.expect(state.macropore_conductance_m3_per_step[aluminum] > 0);
    try std.testing.expectEqual(@as(f64, 0.75), state.micropore_mobility_fraction[phosphate]);
}
