const std = @import("std");

pub const SimulationConfig = struct {
    grid_columns: usize,
    grid_rows: usize,
    soil_layers: usize,
    plant_populations: usize,
    worker_threads: usize,
    tile_cells: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    max_nonlinear_iterations: u16,
    picard_relaxation: f64,

    pub const Dimensions = struct {
        grid_columns: usize,
        grid_rows: usize,
        soil_layers: usize,
        plant_populations: usize,
    };

    pub const Execution = struct {
        worker_threads: usize,
        tile_cells: usize,
    };

    pub const Numerics = struct {
        relative_tolerance: f64,
        absolute_tolerance: f64,
        max_nonlinear_iterations: u16,
        picard_relaxation: f64 = 0.5,
    };

    pub fn init(dimensions: Dimensions, execution: Execution, numerics: Numerics) !SimulationConfig {
        const result = SimulationConfig{
            .grid_columns = dimensions.grid_columns,
            .grid_rows = dimensions.grid_rows,
            .soil_layers = dimensions.soil_layers,
            .plant_populations = dimensions.plant_populations,
            .worker_threads = execution.worker_threads,
            .tile_cells = execution.tile_cells,
            .relative_tolerance = numerics.relative_tolerance,
            .absolute_tolerance = numerics.absolute_tolerance,
            .max_nonlinear_iterations = numerics.max_nonlinear_iterations,
            .picard_relaxation = numerics.picard_relaxation,
        };
        try result.validate();
        return result;
    }

    pub fn validate(self: SimulationConfig) !void {
        if (self.grid_columns == 0 or self.grid_rows == 0) return error.EmptyGrid;
        if (self.soil_layers == 0) return error.NoSoilLayers;
        if (self.plant_populations == 0) return error.NoPlantSpecies;
        if (self.worker_threads == 0) return error.NoWorkerThreads;
        if (self.tile_cells == 0) return error.EmptyTile;
        if (!std.math.isFinite(self.relative_tolerance) or self.relative_tolerance <= 0) return error.InvalidRelativeTolerance;
        if (!std.math.isFinite(self.absolute_tolerance) or self.absolute_tolerance <= 0) return error.InvalidAbsoluteTolerance;
        if (self.max_nonlinear_iterations == 0) return error.NoNonlinearIterations;
        if (!std.math.isFinite(self.picard_relaxation) or self.picard_relaxation <= 0 or self.picard_relaxation > 1) return error.InvalidPicardRelaxation;
        _ = try std.math.mul(usize, self.grid_columns, self.grid_rows);
        const cells = try std.math.mul(usize, self.grid_columns, self.grid_rows);
        const cell_species = try std.math.mul(usize, cells, self.plant_populations);
        _ = try std.math.mul(usize, cell_species, self.soil_layers);
    }
};

test "explicit runtime configuration is valid" {
    _ = try SimulationConfig.init(
        .{ .grid_columns = 13, .grid_rows = 7, .soil_layers = 23, .plant_populations = 11 },
        .{ .worker_threads = 3, .tile_cells = 97 },
        .{ .relative_tolerance = 1.0e-8, .absolute_tolerance = 1.0e-11, .max_nonlinear_iterations = 40 },
    );
}
