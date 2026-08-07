const std = @import("std");

pub const SimulationConfig = struct {
    lon_count: usize,
    lat_count: usize,
    soil_layers: usize,
    plant_populations: usize,
    worker_threads: usize,
    tile_cells: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    /// Deviation permitted by the landscape mass-balance audit, in the audited
    /// species' own units (MJ/m2, g/m2, m). Distinct from `absolute_tolerance`:
    /// the audit measures conservation error of a completed step, whereas
    /// `absolute_tolerance` decides when a nonlinear solve has converged.
    /// Tying them together makes the audit untunable, because loosening it to
    /// observe a trajectory also stops the solvers from iterating.
    mass_balance_tolerance: f64,
    /// Magnitude below which a physical quantity is treated as absent (leaf
    /// area, plant mass, water volume, combustible carbon). Carries physical
    /// units, not residual units, so it is unrelated to solver convergence.
    negligible_quantity_threshold: f64,
    max_nonlinear_iterations: u16,
    picard_relaxation: f64,

    pub const Dimensions = struct {
        lon_count: usize,
        lat_count: usize,
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
        /// Defaults to `absolute_tolerance` so existing callers are unchanged.
        mass_balance_tolerance: ?f64 = null,
        /// Defaults to `absolute_tolerance` so existing callers are unchanged.
        negligible_quantity_threshold: ?f64 = null,
        max_nonlinear_iterations: u16,
        picard_relaxation: f64 = 0.5,
    };

    pub fn init(dimensions: Dimensions, execution: Execution, numerics: Numerics) !SimulationConfig {
        const result = SimulationConfig{
            .lon_count = dimensions.lon_count,
            .lat_count = dimensions.lat_count,
            .soil_layers = dimensions.soil_layers,
            .plant_populations = dimensions.plant_populations,
            .worker_threads = execution.worker_threads,
            .tile_cells = execution.tile_cells,
            .relative_tolerance = numerics.relative_tolerance,
            .absolute_tolerance = numerics.absolute_tolerance,
            .mass_balance_tolerance = numerics.mass_balance_tolerance orelse numerics.absolute_tolerance,
            .negligible_quantity_threshold = numerics.negligible_quantity_threshold orelse numerics.absolute_tolerance,
            .max_nonlinear_iterations = numerics.max_nonlinear_iterations,
            .picard_relaxation = numerics.picard_relaxation,
        };
        try result.validate();
        return result;
    }

    pub fn validate(self: SimulationConfig) !void {
        if (self.lon_count == 0 or self.lat_count == 0) return error.EmptyGrid;
        if (self.soil_layers == 0) return error.NoSoilLayers;
        if (self.plant_populations == 0) return error.NoPlantSpecies;
        if (self.worker_threads == 0) return error.NoWorkerThreads;
        if (self.tile_cells == 0) return error.EmptyTile;
        if (!std.math.isFinite(self.relative_tolerance) or self.relative_tolerance <= 0) return error.InvalidRelativeTolerance;
        if (!std.math.isFinite(self.absolute_tolerance) or self.absolute_tolerance <= 0) return error.InvalidAbsoluteTolerance;
        if (!std.math.isFinite(self.mass_balance_tolerance) or self.mass_balance_tolerance <= 0) return error.InvalidMassBalanceTolerance;
        if (!std.math.isFinite(self.negligible_quantity_threshold) or self.negligible_quantity_threshold <= 0) return error.InvalidNegligibleQuantityThreshold;
        if (self.max_nonlinear_iterations == 0) return error.NoNonlinearIterations;
        if (!std.math.isFinite(self.picard_relaxation) or self.picard_relaxation <= 0 or self.picard_relaxation > 1) return error.InvalidPicardRelaxation;
        _ = try std.math.mul(usize, self.lon_count, self.lat_count);
        const cells = try std.math.mul(usize, self.lon_count, self.lat_count);
        const cell_species = try std.math.mul(usize, cells, self.plant_populations);
        _ = try std.math.mul(usize, cell_species, self.soil_layers);
    }
};

test "explicit runtime configuration is valid" {
    _ = try SimulationConfig.init(
        .{ .lon_count = 13, .lat_count = 7, .soil_layers = 23, .plant_populations = 11 },
        .{ .worker_threads = 3, .tile_cells = 97 },
        .{ .relative_tolerance = 1.0e-8, .absolute_tolerance = 1.0e-11, .max_nonlinear_iterations = 40 },
    );
}
