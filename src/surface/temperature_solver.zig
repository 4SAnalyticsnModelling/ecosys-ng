const std = @import("std");
const CellRange = @import("../core/compute.zig").CellRange;
const freeze_thaw_energy_limit = @import("litter_freeze_thaw_energy_limit.zig");
const numerics = @import("../core/numerics.zig");
const GridState = @import("../state/grid.zig").GridState;
const AtmosphericState = @import("../atmosphere/atmospheric_forcing.zig").State;
const GroundRadiationState = @import("ground_radiation.zig").State;
const SurfaceEnergyState = @import("energy.zig").State;
const ExposureState = @import("../canopy/radiation/exposure.zig").State;
const SoilThermalState = @import("../soil/heat/thermal.zig").State;
const phase_change = @import("../soil/water/phase_change.zig");
const retention = @import("../soil/water/retention.zig");

pub const Settings = struct {
    sensible_heat_conductance_megajoules_per_m2_h_k: f64,
    latent_heat_conductance_megajoules_per_m2_h_kpa: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    latent_heat_of_vaporization_megajoules_per_m3: f64,
    surface_vapor_activity_fraction: f64,
    timestep_hours: f64,
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
    solver_options: numerics.SolverOptions,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    equilibrium_temperature_k: []f64,
    energy_residual_megajoules_per_m2: []f64,
    residual_tolerance_megajoules_per_m2: []f64,
    sensible_heat_flux_megajoules_per_m2: []f64,
    latent_heat_flux_megajoules_per_m2: []f64,
    vapor_sensible_heat_flux_megajoules_per_m2: []f64,
    conductive_heat_flux_megajoules_per_m2: []f64,
    storage_heat_flux_megajoules_per_m2: []f64,
    phase_heat_flux_megajoules_per_m2: []f64,
    liquid_water_change_m3: []f64,
    ice_water_equivalent_change_m3: []f64,
    iteration_count: []u16,
    newton_raphson_step_count: []u16,
    picard_step_count: []u16,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.EmptySurfaceTemperatureGrid;
        var result: State = .{
            .allocator = allocator,
            .cell_count = cell_count,
            .equilibrium_temperature_k = try allocator.alloc(f64, cell_count),
            .energy_residual_megajoules_per_m2 = undefined,
            .residual_tolerance_megajoules_per_m2 = undefined,
            .sensible_heat_flux_megajoules_per_m2 = undefined,
            .latent_heat_flux_megajoules_per_m2 = undefined,
            .vapor_sensible_heat_flux_megajoules_per_m2 = undefined,
            .conductive_heat_flux_megajoules_per_m2 = undefined,
            .storage_heat_flux_megajoules_per_m2 = undefined,
            .phase_heat_flux_megajoules_per_m2 = undefined,
            .liquid_water_change_m3 = undefined,
            .ice_water_equivalent_change_m3 = undefined,
            .iteration_count = undefined,
            .newton_raphson_step_count = undefined,
            .picard_step_count = undefined,
        };
        errdefer allocator.free(result.equilibrium_temperature_k);
        result.energy_residual_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.energy_residual_megajoules_per_m2);
        result.residual_tolerance_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.residual_tolerance_megajoules_per_m2);
        result.sensible_heat_flux_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.sensible_heat_flux_megajoules_per_m2);
        result.latent_heat_flux_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.latent_heat_flux_megajoules_per_m2);
        result.vapor_sensible_heat_flux_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.vapor_sensible_heat_flux_megajoules_per_m2);
        result.conductive_heat_flux_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.conductive_heat_flux_megajoules_per_m2);
        result.storage_heat_flux_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.storage_heat_flux_megajoules_per_m2);
        result.phase_heat_flux_megajoules_per_m2 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.phase_heat_flux_megajoules_per_m2);
        result.liquid_water_change_m3 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.liquid_water_change_m3);
        result.ice_water_equivalent_change_m3 = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.ice_water_equivalent_change_m3);
        result.iteration_count = try allocator.alloc(u16, cell_count);
        errdefer allocator.free(result.iteration_count);
        result.newton_raphson_step_count = try allocator.alloc(u16, cell_count);
        errdefer allocator.free(result.newton_raphson_step_count);
        result.picard_step_count = try allocator.alloc(u16, cell_count);
        @memset(result.equilibrium_temperature_k, 273.15);
        @memset(result.energy_residual_megajoules_per_m2, 0);
        @memset(result.residual_tolerance_megajoules_per_m2, 0);
        @memset(result.sensible_heat_flux_megajoules_per_m2, 0);
        @memset(result.latent_heat_flux_megajoules_per_m2, 0);
        @memset(result.vapor_sensible_heat_flux_megajoules_per_m2, 0);
        @memset(result.conductive_heat_flux_megajoules_per_m2, 0);
        @memset(result.storage_heat_flux_megajoules_per_m2, 0);
        @memset(result.phase_heat_flux_megajoules_per_m2, 0);
        @memset(result.liquid_water_change_m3, 0);
        @memset(result.ice_water_equivalent_change_m3, 0);
        @memset(result.iteration_count, 0);
        @memset(result.newton_raphson_step_count, 0);
        @memset(result.picard_step_count, 0);
        return result;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.picard_step_count);
        self.allocator.free(self.newton_raphson_step_count);
        self.allocator.free(self.iteration_count);
        self.allocator.free(self.storage_heat_flux_megajoules_per_m2);
        self.allocator.free(self.phase_heat_flux_megajoules_per_m2);
        self.allocator.free(self.liquid_water_change_m3);
        self.allocator.free(self.ice_water_equivalent_change_m3);
        self.allocator.free(self.conductive_heat_flux_megajoules_per_m2);
        self.allocator.free(self.vapor_sensible_heat_flux_megajoules_per_m2);
        self.allocator.free(self.latent_heat_flux_megajoules_per_m2);
        self.allocator.free(self.sensible_heat_flux_megajoules_per_m2);
        self.allocator.free(self.residual_tolerance_megajoules_per_m2);
        self.allocator.free(self.energy_residual_megajoules_per_m2);
        self.allocator.free(self.equilibrium_temperature_k);
        self.* = undefined;
    }

    /// Clears accepted phase changes before a potentially failing tile pass.
    /// This makes a caught nonlinear failure publish a zero change for every
    /// cell that was not successfully committed during the current pass.
    pub fn resetPhaseChangeDiagnostics(self: *State) void {
        @memset(self.liquid_water_change_m3, 0);
        @memset(self.ice_water_equivalent_change_m3, 0);
        @memset(self.phase_heat_flux_megajoules_per_m2, 0);
    }

    pub fn validateFinite(self: State) !void {
        for (self.equilibrium_temperature_k, 0..) |value, cell| if (!std.math.isFinite(value) or value <= 0) {
            std.log.err("invalid solved surface temperature: cell={d} value={e}", .{ cell, value });
            return error.InvalidSolvedSurfaceTemperature;
        };
        for (self.energy_residual_megajoules_per_m2, 0..) |value, cell| if (!std.math.isFinite(value)) {
            std.log.err("non-finite surface energy residual: cell={d} value={e}", .{ cell, value });
            return error.NonFiniteSurfaceTemperatureResidual;
        };
        inline for (.{ self.sensible_heat_flux_megajoules_per_m2, self.latent_heat_flux_megajoules_per_m2, self.vapor_sensible_heat_flux_megajoules_per_m2, self.conductive_heat_flux_megajoules_per_m2, self.storage_heat_flux_megajoules_per_m2, self.phase_heat_flux_megajoules_per_m2, self.liquid_water_change_m3, self.ice_water_equivalent_change_m3 }) |values| for (values, 0..) |value, cell| {
            if (!std.math.isFinite(value)) {
                std.log.err("non-finite surface heat/phase diagnostic: cell={d} value={e}", .{ cell, value });
                return error.NonFiniteSurfaceHeatDiagnostic;
            }
        };
    }

    pub const Diagnostics = struct {
        maximum_absolute_residual_megajoules_per_m2: f64,
        total_iterations: u64,
        total_newton_raphson_steps: u64,
        total_picard_steps: u64,
        cells_using_picard: usize,
    };

    pub fn validateConvergence(self: State) !Diagnostics {
        var diagnostics: Diagnostics = .{
            .maximum_absolute_residual_megajoules_per_m2 = 0,
            .total_iterations = 0,
            .total_newton_raphson_steps = 0,
            .total_picard_steps = 0,
            .cells_using_picard = 0,
        };
        for (0..self.cell_count) |cell| {
            const absolute_residual = @abs(self.energy_residual_megajoules_per_m2[cell]);
            const tolerance = self.residual_tolerance_megajoules_per_m2[cell];
            if (!std.math.isFinite(absolute_residual) or absolute_residual > tolerance) {
                std.log.err("surface nonlinear convergence audit failed: cell={d} residual_megajoules_per_m2={e} tolerance={e} temperature_k={e}", .{ cell, absolute_residual, tolerance, self.equilibrium_temperature_k[cell] });
                return error.SurfaceTemperatureResidualTooLarge;
            }
            diagnostics.maximum_absolute_residual_megajoules_per_m2 = @max(diagnostics.maximum_absolute_residual_megajoules_per_m2, absolute_residual);
            diagnostics.total_iterations += self.iteration_count[cell];
            diagnostics.total_newton_raphson_steps += self.newton_raphson_step_count[cell];
            diagnostics.total_picard_steps += self.picard_step_count[cell];
            if (self.picard_step_count[cell] > 0) diagnostics.cells_using_picard += 1;
        }
        return diagnostics;
    }
};

pub const ApplyContext = struct {
    result: *State,
    grid: *GridState,
    atmosphere: *const AtmosphericState,
    air_temperature_k: []const f64,
    air_vapor_pressure_kpa: []const f64,
    ground_radiation: *const GroundRadiationState,
    surface_energy: *SurfaceEnergyState,
    soil_thermal: *const SoilThermalState,
    /// Complete litter/surface heat capacity, MJ K-1 per cell.
    surface_heat_capacity_megajoules_per_k: []const f64,
    exposure: ?*const ExposureState,
    external_heat_megajoules_per_m2: []const f64,
    surface_phase: ?SurfacePhaseContext = null,
    settings: Settings,
};

pub const SurfacePhaseContext = struct {
    liquid_water_m3: []f64,
    ice_water_equivalent_m3: []f64,
    retention_capacity_m3: []const f64,
    horizontal_area_m2: []const f64,
    residual_water_content_m3_per_m3: f64,
    van_genuchten_alpha_per_m: f64,
    van_genuchten_n: f64,
    gravitational_water_potential_mpa_per_m: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
    pure_water_melting_temperature_k: f64,
};

const ResidualContext = struct {
    absorbed_shortwave_megajoules_per_m2: f64,
    atmospheric_longwave_megajoules_per_m2: f64,
    air_temperature_k: f64,
    ground_exposure_fraction: f64,
    emissivity: f64,
    sensible_heat_conductance_megajoules_per_m2_h_k: f64,
    latent_heat_conductance_megajoules_per_m2_h_kpa: f64,
    atmospheric_vapor_pressure_kpa: f64,
    surface_vapor_activity_fraction: f64,
    subsurface_temperature_k: f64,
    previous_surface_temperature_k: f64,
    conductive_heat_conductance_megajoules_per_m2_h_k: f64,
    storage_heat_conductance_megajoules_per_m2_h_k: f64,
    external_heat_megajoules_per_m2: f64 = 0,
    phase: ?CellPhaseContext = null,
};

const CellPhaseContext = struct {
    initial_liquid_water_m3: f64,
    initial_ice_water_equivalent_m3: f64,
    porous_medium_volume_m3: f64,
    horizontal_area_m2: f64,
    parameters: retention.MualemVanGenuchtenParameters,
    gravitational_water_potential_mpa_per_m: f64,
    latent_heat_of_fusion_megajoules_per_m3: f64,
    pure_water_melting_temperature_k: f64,
    /// Surface volumetric heat capacity, the oracle's `VHCPR2`. `watsub.f` 3150
    /// limits the phase change to the latent heat the temperature deficit can
    /// drive, which is proportional to this capacity.
    heat_capacity_megajoules_per_k: f64,
};

pub fn applyTile(context: *ApplyContext, range: CellRange) !void {
    try validateSettings(context.settings);
    const result = context.result;
    if (range.end > result.cell_count or context.grid.cell_count != result.cell_count or context.atmosphere.cell_count != result.cell_count or context.ground_radiation.cell_count != result.cell_count or context.surface_energy.cell_count != result.cell_count or context.soil_thermal.cell_count != result.cell_count or context.air_temperature_k.len != result.cell_count or context.air_vapor_pressure_kpa.len != result.cell_count or context.external_heat_megajoules_per_m2.len != result.cell_count or context.surface_heat_capacity_megajoules_per_k.len != result.cell_count) return error.SurfaceTemperatureDimensionMismatch;
    if (context.exposure) |exposure| if (exposure.cell_count != result.cell_count) return error.SurfaceTemperatureDimensionMismatch;
    if (context.surface_phase) |surface_phase|
        try validateSurfacePhaseContext(surface_phase, result.cell_count);
    for (range.first..range.end) |cell| {
        // A failed nonlinear solve must not expose a phase change retained
        // from an earlier timestep. The prognostic water/ice carriers are
        // committed only below, after convergence.
        result.liquid_water_change_m3[cell] = 0;
        result.ice_water_equivalent_change_m3[cell] = 0;
        result.phase_heat_flux_megajoules_per_m2[cell] = 0;
        const surface_layer_index = cell * context.grid.soil_layer_capacity;
        const layer_thickness_m = context.soil_thermal.layer_thickness_m[surface_layer_index];
        if (!std.math.isFinite(layer_thickness_m) or layer_thickness_m <= 0) return error.InvalidSurfaceSoilLayerThickness;
        const horizontal_area_m2 = if (context.surface_phase) |phase| phase.horizontal_area_m2[cell] else return error.MissingSurfaceHeatCapacityArea;
        if (!std.math.isFinite(horizontal_area_m2) or horizontal_area_m2 <= 0 or !std.math.isFinite(context.surface_heat_capacity_megajoules_per_k[cell]) or context.surface_heat_capacity_megajoules_per_k[cell] <= 0) return error.InvalidSurfaceHeatCapacity;
        const residual_context: ResidualContext = .{
            .absorbed_shortwave_megajoules_per_m2 = context.ground_radiation.absorbed_shortwave_megajoules_per_m2[cell],
            .atmospheric_longwave_megajoules_per_m2 = context.atmosphere.longwave_radiation_megajoules_per_m2[cell],
            .air_temperature_k = context.air_temperature_k[cell],
            .ground_exposure_fraction = if (context.exposure) |exposure| exposure.ground_exposure_fraction[cell] else 1,
            .emissivity = context.surface_energy.surface_emissivity[cell],
            .sensible_heat_conductance_megajoules_per_m2_h_k = context.settings.sensible_heat_conductance_megajoules_per_m2_h_k,
            .latent_heat_conductance_megajoules_per_m2_h_kpa = context.settings.latent_heat_conductance_megajoules_per_m2_h_kpa,
            .atmospheric_vapor_pressure_kpa = context.air_vapor_pressure_kpa[cell],
            .surface_vapor_activity_fraction = context.settings.surface_vapor_activity_fraction,
            .subsurface_temperature_k = context.grid.soil_temperature_k[surface_layer_index],
            .previous_surface_temperature_k = context.grid.surface_temperature_k[cell],
            // Couple the sequential surface and soil implicit solves with the
            // same accepted conductance. A conductance larger than the
            // topsoil's one-step sensible-heat capacity can demand a soil
            // temperature outside its physical bracket. Limiting the
            // conductance here (rather than clipping the committed flux later)
            // keeps the surface residual and soil source equal and opposite.
            .conductive_heat_conductance_megajoules_per_m2_h_k = @min(
                context.soil_thermal.thermal_conductivity_m_megajoules_per_h_k[surface_layer_index] /
                    (0.5 * layer_thickness_m),
                context.soil_thermal.total_heat_capacity_megajoules_per_m3_k[surface_layer_index] *
                    layer_thickness_m / context.settings.timestep_hours,
            ),
            .storage_heat_conductance_megajoules_per_m2_h_k = context.surface_heat_capacity_megajoules_per_k[cell] / horizontal_area_m2 / context.settings.timestep_hours,
            .external_heat_megajoules_per_m2 = context.external_heat_megajoules_per_m2[cell],
            .phase = if (context.surface_phase) |surface_phase| phase: {
                const total_water_equivalent_m3 =
                    surface_phase.liquid_water_m3[cell] +
                    surface_phase.ice_water_equivalent_m3[cell];
                const porous_medium_volume_m3 = @max(
                    surface_phase.retention_capacity_m3[cell],
                    total_water_equivalent_m3,
                );
                if (porous_medium_volume_m3 == 0) break :phase null;
                break :phase .{
                    .initial_liquid_water_m3 = surface_phase.liquid_water_m3[cell],
                    .initial_ice_water_equivalent_m3 = surface_phase.ice_water_equivalent_m3[cell],
                    .porous_medium_volume_m3 = porous_medium_volume_m3,
                    .horizontal_area_m2 = surface_phase.horizontal_area_m2[cell],
                    .parameters = .{
                        .residual_water_content_m3_per_m3 = surface_phase.residual_water_content_m3_per_m3,
                        .saturated_water_content_m3_per_m3 = 1,
                        .alpha_per_m = surface_phase.van_genuchten_alpha_per_m,
                        .n = surface_phase.van_genuchten_n,
                        .saturated_hydraulic_conductivity_m_per_h = 0,
                    },
                    .gravitational_water_potential_mpa_per_m = surface_phase.gravitational_water_potential_mpa_per_m,
                    .latent_heat_of_fusion_megajoules_per_m3 = surface_phase.latent_heat_of_fusion_megajoules_per_m3,
                    .pure_water_melting_temperature_k = surface_phase.pure_water_melting_temperature_k,
                    .heat_capacity_megajoules_per_k = context.surface_heat_capacity_megajoules_per_k[cell],
                };
            } else null,
        };
        if (!std.math.isFinite(residual_context.external_heat_megajoules_per_m2)) return error.NonFiniteSurfaceExternalHeat;
        var cell_solver_options = context.settings.solver_options;
        cell_solver_options.residual_scale = residualScale(residual_context, context.grid.surface_temperature_k[cell]);
        cell_solver_options.safeguard_with_bracket = true;
        const solved = numerics.newtonPicard(
            residual_context,
            residual,
            derivative,
            picard,
            context.settings.minimum_temperature_k,
            context.settings.maximum_temperature_k,
            context.grid.surface_temperature_k[cell],
            cell_solver_options,
        ) catch |err| {
            std.log.err("surface temperature solve failed: cell={d} initial_temperature_k={e}", .{ cell, context.grid.surface_temperature_k[cell] });
            return err;
        };
        result.equilibrium_temperature_k[cell] = solved.root;
        result.energy_residual_megajoules_per_m2[cell] = solved.residual;
        result.residual_tolerance_megajoules_per_m2[cell] = cell_solver_options.absolute_tolerance + cell_solver_options.relative_tolerance * cell_solver_options.residual_scale;
        result.sensible_heat_flux_megajoules_per_m2[cell] = residual_context.sensible_heat_conductance_megajoules_per_m2_h_k * (residual_context.air_temperature_k - solved.root);
        result.latent_heat_flux_megajoules_per_m2[cell] = latentHeatFlux(residual_context, solved.root);
        result.vapor_sensible_heat_flux_megajoules_per_m2[cell] =
            try vaporSensibleHeatFlux(
                result.latent_heat_flux_megajoules_per_m2[cell],
                residual_context.air_temperature_k,
                solved.root,
                context.settings.liquid_water_heat_capacity_megajoules_per_m3_k,
                context.settings.latent_heat_of_vaporization_megajoules_per_m3,
            );
        result.conductive_heat_flux_megajoules_per_m2[cell] = residual_context.conductive_heat_conductance_megajoules_per_m2_h_k * (residual_context.subsurface_temperature_k - solved.root);
        result.storage_heat_flux_megajoules_per_m2[cell] = residual_context.storage_heat_conductance_megajoules_per_m2_h_k * (residual_context.previous_surface_temperature_k - solved.root);
        const accepted_emitted_longwave_megajoules_per_m2 = residual_context.emissivity * 2.04e-10 * std.math.pow(f64, solved.root, 4) * residual_context.ground_exposure_fraction;
        context.surface_energy.emitted_sky_longwave_megajoules_per_m2[cell] = accepted_emitted_longwave_megajoules_per_m2;
        context.surface_energy.net_longwave_megajoules_per_m2[cell] = context.surface_energy.downward_sky_longwave_megajoules_per_m2[cell] - accepted_emitted_longwave_megajoules_per_m2;
        context.surface_energy.net_radiation_megajoules_per_m2[cell] = residual_context.absorbed_shortwave_megajoules_per_m2 + context.surface_energy.net_longwave_megajoules_per_m2[cell];
        if (context.surface_phase) |surface_phase| if (residual_context.phase != null) {
            const equilibrium = try surfacePhaseEquilibrium(
                residual_context.phase.?,
                solved.root,
            );
            result.liquid_water_change_m3[cell] =
                equilibrium.liquid_water_m3 -
                surface_phase.liquid_water_m3[cell];
            result.ice_water_equivalent_change_m3[cell] =
                equilibrium.ice_water_equivalent_m3 -
                surface_phase.ice_water_equivalent_m3[cell];
            result.phase_heat_flux_megajoules_per_m2[cell] =
                residual_context.phase.?.latent_heat_of_fusion_megajoules_per_m3 *
                result.ice_water_equivalent_change_m3[cell] /
                residual_context.phase.?.horizontal_area_m2;
            surface_phase.liquid_water_m3[cell] =
                equilibrium.liquid_water_m3;
            surface_phase.ice_water_equivalent_m3[cell] =
                equilibrium.ice_water_equivalent_m3;
        } else {
            result.liquid_water_change_m3[cell] = 0;
            result.ice_water_equivalent_change_m3[cell] = 0;
            result.phase_heat_flux_megajoules_per_m2[cell] = 0;
        } else {
            result.liquid_water_change_m3[cell] = 0;
            result.ice_water_equivalent_change_m3[cell] = 0;
            result.phase_heat_flux_megajoules_per_m2[cell] = 0;
        }
        result.iteration_count[cell] = solved.iterations;
        result.newton_raphson_step_count[cell] = solved.newton_raphson_steps;
        result.picard_step_count[cell] = solved.picard_steps;
        // Each tile owns its cells, so committing the converged prognostic
        // temperature is race-free under the CPU executor and GPU-ready.
        context.grid.surface_temperature_k[cell] = solved.root;
    }
}

fn residual(context: ResidualContext, temperature_k: f64) f64 {
    const downward_longwave = context.atmospheric_longwave_megajoules_per_m2 * context.ground_exposure_fraction;
    const emitted_longwave = context.emissivity * 2.04e-10 * std.math.pow(f64, temperature_k, 4) * context.ground_exposure_fraction;
    const sensible_heat = context.sensible_heat_conductance_megajoules_per_m2_h_k * (context.air_temperature_k - temperature_k);
    const conductive_heat = context.conductive_heat_conductance_megajoules_per_m2_h_k * (context.subsurface_temperature_k - temperature_k);
    const storage_heat = context.storage_heat_conductance_megajoules_per_m2_h_k * (context.previous_surface_temperature_k - temperature_k);
    return context.absorbed_shortwave_megajoules_per_m2 + context.external_heat_megajoules_per_m2 + downward_longwave - emitted_longwave + sensible_heat + latentHeatFlux(context, temperature_k) + conductive_heat + storage_heat + phaseHeatFlux(context, temperature_k);
}

fn derivative(context: ResidualContext, temperature_k: f64) f64 {
    const saturation_pressure = saturationVaporPressureKpa(temperature_k);
    const saturation_derivative = saturation_pressure * 5360.0 / (temperature_k * temperature_k);
    return -4.0 * context.emissivity * 2.04e-10 * std.math.pow(f64, temperature_k, 3) * context.ground_exposure_fraction - context.sensible_heat_conductance_megajoules_per_m2_h_k - context.latent_heat_conductance_megajoules_per_m2_h_kpa * context.surface_vapor_activity_fraction * saturation_derivative - context.conductive_heat_conductance_megajoules_per_m2_h_k - context.storage_heat_conductance_megajoules_per_m2_h_k + phaseHeatDerivative(context, temperature_k);
}

fn picard(context: ResidualContext, temperature_k: f64) f64 {
    const downward_longwave = context.atmospheric_longwave_megajoules_per_m2 * context.ground_exposure_fraction;
    const emitted_longwave = context.emissivity * 2.04e-10 * std.math.pow(f64, temperature_k, 4) * context.ground_exposure_fraction;
    const total_linear_conductance = context.sensible_heat_conductance_megajoules_per_m2_h_k + context.conductive_heat_conductance_megajoules_per_m2_h_k + context.storage_heat_conductance_megajoules_per_m2_h_k;
    return (context.sensible_heat_conductance_megajoules_per_m2_h_k * context.air_temperature_k + context.conductive_heat_conductance_megajoules_per_m2_h_k * context.subsurface_temperature_k + context.storage_heat_conductance_megajoules_per_m2_h_k * context.previous_surface_temperature_k + context.absorbed_shortwave_megajoules_per_m2 + context.external_heat_megajoules_per_m2 + downward_longwave - emitted_longwave + latentHeatFlux(context, temperature_k) + phaseHeatFlux(context, temperature_k)) / total_linear_conductance;
}

fn saturationVaporPressureKpa(temperature_k: f64) f64 {
    return 0.61 * @exp(5360.0 * (3.661e-3 - 1.0 / temperature_k));
}

fn latentHeatFlux(context: ResidualContext, temperature_k: f64) f64 {
    return context.latent_heat_conductance_megajoules_per_m2_h_kpa * (context.atmospheric_vapor_pressure_kpa - context.surface_vapor_activity_fraction * saturationVaporPressureKpa(temperature_k));
}

fn vaporSensibleHeatFlux(
    latent_heat_flux_megajoules_per_m2: f64,
    air_temperature_k: f64,
    surface_temperature_k: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    latent_heat_of_vaporization_megajoules_per_m3: f64,
) !f64 {
    inline for (.{ latent_heat_flux_megajoules_per_m2, air_temperature_k, surface_temperature_k, liquid_water_heat_capacity_megajoules_per_m3_k, latent_heat_of_vaporization_megajoules_per_m3 }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceVaporSensibleHeatInput;
    if (air_temperature_k <= 0 or surface_temperature_k <= 0 or
        liquid_water_heat_capacity_megajoules_per_m3_k <= 0 or
        latent_heat_of_vaporization_megajoules_per_m3 <= 0)
        return error.InvalidSurfaceVaporSensibleHeatInput;
    const water_flux_into_surface_m3_per_m2 =
        latent_heat_flux_megajoules_per_m2 / latent_heat_of_vaporization_megajoules_per_m3;
    const donor_temperature_k = if (water_flux_into_surface_m3_per_m2 >= 0)
        air_temperature_k
    else
        surface_temperature_k;
    const heat_flux_megajoules_per_m2 = water_flux_into_surface_m3_per_m2 *
        liquid_water_heat_capacity_megajoules_per_m3_k * donor_temperature_k;
    if (!std.math.isFinite(heat_flux_megajoules_per_m2))
        return error.NonFiniteSurfaceVaporSensibleHeatResult;
    return heat_flux_megajoules_per_m2;
}

fn surfacePhaseEquilibrium(
    context: CellPhaseContext,
    temperature_k: f64,
) !phase_change.DallAmicoEquilibrium {
    const total_water_equivalent_m3 =
        context.initial_liquid_water_m3 +
        context.initial_ice_water_equivalent_m3;
    const total_water_content =
        total_water_equivalent_m3 / context.porous_medium_volume_m3;
    if (total_water_content <=
        context.parameters.residual_water_content_m3_per_m3)
    {
        return .{
            .depressed_melting_temperature_k = context.pure_water_melting_temperature_k,
            .liquid_pressure_head_m = 0,
            .liquid_water_m3 = total_water_equivalent_m3,
            .ice_water_equivalent_m3 = 0,
        };
    }
    const unfrozen_pressure_head_m =
        try context.parameters.pressureHeadAtWaterContent(
            std.math.clamp(
                total_water_content,
                context.parameters.residual_water_content_m3_per_m3,
                context.parameters.saturated_water_content_m3_per_m3,
            ),
        );
    const equilibrium = try phase_change.dallAmicoEquilibrium(.{
        .temperature_k = temperature_k,
        .total_water_equivalent_m3 = total_water_equivalent_m3,
        .porous_medium_volume_m3 = context.porous_medium_volume_m3,
        .unfrozen_pressure_head_m = unfrozen_pressure_head_m,
        .gravitational_water_potential_mpa_per_m = context.gravitational_water_potential_mpa_per_m,
        .latent_heat_of_fusion_megajoules_per_m3 = context.latent_heat_of_fusion_megajoules_per_m3,
        .pure_water_melting_temperature_k = context.pure_water_melting_temperature_k,
        .mualem_van_genuchten = context.parameters,
    });
    return limitPhaseChangeToAvailableEnergy(context, temperature_k, total_water_equivalent_m3, equilibrium);
}

/// WATSUB 3150--3156 rate limit. The Dall'Amico equilibrium answers "what split
/// would this temperature imply at equilibrium", but the oracle asks the different
/// question "how much phase change can the temperature deficit actually drive in
/// this step", with the available phase acting only as a cap. Applying the
/// equilibrium directly let the whole litter reservoir freeze in one hour at
/// `19.17x` the energy-permitted rate, which collapsed the litter water carrier by
/// `1.1e7` and made surface evaporation impossible. See EXEC-002 and EXEC-004.
///
/// This clamps the requested ice change to the energy limit. It is applied inside
/// `surfacePhaseEquilibrium` so the Newton residual, its derivative, and the final
/// writeback all see the same limited value and the solve stays consistent.
fn limitPhaseChangeToAvailableEnergy(
    context: CellPhaseContext,
    temperature_k: f64,
    total_water_equivalent_m3: f64,
    equilibrium: phase_change.DallAmicoEquilibrium,
) !phase_change.DallAmicoEquilibrium {
    const requested_ice_change_m3 =
        equilibrium.ice_water_equivalent_m3 - context.initial_ice_water_equivalent_m3;
    if (requested_ice_change_m3 == 0) return equilibrium;
    if (!std.math.isFinite(context.heat_capacity_megajoules_per_k) or
        context.heat_capacity_megajoules_per_k <= 0) return equilibrium;

    const limit = try freeze_thaw_energy_limit.apply(.{
        .temperature_k = temperature_k,
        .heat_capacity_megajoules_per_k = context.heat_capacity_megajoules_per_k,
        .liquid_water_m3 = context.initial_liquid_water_m3,
        .ice_water_equivalent_m3 = context.initial_ice_water_equivalent_m3,
        // The Dall'Amico solve already accounts for the depressed melting point
        // through the retention curve, so express the same depression as the
        // potential the energy limit expects.
        .water_potential_megapascal = 0,
        .substep_energy_fraction = 1,
        .substep_mass_fraction = 1,
        .negligible_volume_m3 = 0,
    }, .{
        .latent_heat_of_fusion_megajoules_per_m3 = context.latent_heat_of_fusion_megajoules_per_m3,
        .freezing_point_depression_numerator = 9.0959e4,
        .freezing_temperature_coefficient = 6.2913e-3,
        .ice_density_ratio = 0.917,
    });
    // `liquid_water_change_m3` is negative when freezing, so the permitted ice
    // change is its negation.
    const permitted_ice_change_m3 = -limit.liquid_water_change_m3;
    const limited_ice_change_m3 = if (requested_ice_change_m3 > 0)
        @min(requested_ice_change_m3, @max(0.0, permitted_ice_change_m3))
    else
        @max(requested_ice_change_m3, @min(0.0, permitted_ice_change_m3));
    const limited_ice_m3 = std.math.clamp(
        context.initial_ice_water_equivalent_m3 + limited_ice_change_m3,
        0,
        total_water_equivalent_m3,
    );
    if (!std.math.isFinite(limited_ice_m3)) return equilibrium;
    return .{
        .depressed_melting_temperature_k = equilibrium.depressed_melting_temperature_k,
        .liquid_pressure_head_m = equilibrium.liquid_pressure_head_m,
        .liquid_water_m3 = total_water_equivalent_m3 - limited_ice_m3,
        .ice_water_equivalent_m3 = limited_ice_m3,
    };
}

fn phaseHeatFlux(context: ResidualContext, temperature_k: f64) f64 {
    const phase = context.phase orelse return 0;
    const equilibrium =
        surfacePhaseEquilibrium(phase, temperature_k) catch return std.math.nan(f64);
    return phase.latent_heat_of_fusion_megajoules_per_m3 *
        (equilibrium.ice_water_equivalent_m3 -
            phase.initial_ice_water_equivalent_m3) /
        phase.horizontal_area_m2;
}

fn phaseHeatDerivative(context: ResidualContext, temperature_k: f64) f64 {
    const phase = context.phase orelse return 0;
    const equilibrium =
        surfacePhaseEquilibrium(phase, temperature_k) catch
            return std.math.nan(f64);
    if (equilibrium.ice_water_equivalent_m3 <= 0 or
        temperature_k >= equilibrium.depressed_melting_temperature_k)
        return 0;
    const water_capacity_per_m =
        phase.parameters.waterCapacityPerM(
            equilibrium.liquid_pressure_head_m,
        ) catch return std.math.nan(f64);
    const liquid_water_change_m3_per_k =
        phase.porous_medium_volume_m3 * water_capacity_per_m *
        phase.latent_heat_of_fusion_megajoules_per_m3 /
        (phase.gravitational_water_potential_mpa_per_m * temperature_k);
    return -phase.latent_heat_of_fusion_megajoules_per_m3 *
        liquid_water_change_m3_per_k / phase.horizontal_area_m2;
}

fn validateSurfacePhaseContext(
    context: SurfacePhaseContext,
    cell_count: usize,
) !void {
    if (context.liquid_water_m3.len != cell_count or
        context.ice_water_equivalent_m3.len != cell_count or
        context.retention_capacity_m3.len != cell_count or
        context.horizontal_area_m2.len != cell_count)
        return error.SurfaceTemperatureDimensionMismatch;
    inline for (.{
        context.residual_water_content_m3_per_m3,
        context.van_genuchten_alpha_per_m,
        context.van_genuchten_n,
        context.gravitational_water_potential_mpa_per_m,
        context.latent_heat_of_fusion_megajoules_per_m3,
        context.pure_water_melting_temperature_k,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSurfacePhaseParameter;
    if (context.residual_water_content_m3_per_m3 < 0 or
        context.residual_water_content_m3_per_m3 >= 1 or
        context.van_genuchten_alpha_per_m <= 0 or
        context.van_genuchten_n <= 1 or
        context.gravitational_water_potential_mpa_per_m <= 0 or
        context.latent_heat_of_fusion_megajoules_per_m3 <= 0 or
        context.pure_water_melting_temperature_k <= 0)
        return error.InvalidSurfacePhaseParameter;
    for (context.liquid_water_m3, context.ice_water_equivalent_m3, context.retention_capacity_m3, context.horizontal_area_m2) |liquid, ice, capacity, area| {
        if (!std.math.isFinite(liquid) or liquid < 0 or
            !std.math.isFinite(ice) or ice < 0 or
            !std.math.isFinite(capacity) or capacity < 0 or
            !std.math.isFinite(area) or area <= 0)
            return error.InvalidSurfacePhaseState;
    }
}

fn residualScale(context: ResidualContext, temperature_k: f64) f64 {
    const downward_longwave = context.atmospheric_longwave_megajoules_per_m2 * context.ground_exposure_fraction;
    const emitted_longwave = context.emissivity * 2.04e-10 * std.math.pow(f64, temperature_k, 4) * context.ground_exposure_fraction;
    const sensible_heat = context.sensible_heat_conductance_megajoules_per_m2_h_k * (context.air_temperature_k - temperature_k);
    const latent_heat = latentHeatFlux(context, temperature_k);
    const conductive_heat = context.conductive_heat_conductance_megajoules_per_m2_h_k * (context.subsurface_temperature_k - temperature_k);
    const storage_heat = context.storage_heat_conductance_megajoules_per_m2_h_k * (context.previous_surface_temperature_k - temperature_k);
    const phase_heat = phaseHeatFlux(context, temperature_k);
    return @max(1.0e-12, @abs(context.absorbed_shortwave_megajoules_per_m2) + @abs(context.external_heat_megajoules_per_m2) + @abs(downward_longwave) + @abs(emitted_longwave) + @abs(sensible_heat) + @abs(latent_heat) + @abs(conductive_heat) + @abs(storage_heat) + @abs(phase_heat));
}

fn validateSettings(settings: Settings) !void {
    if (!std.math.isFinite(settings.sensible_heat_conductance_megajoules_per_m2_h_k) or settings.sensible_heat_conductance_megajoules_per_m2_h_k <= 0 or !std.math.isFinite(settings.latent_heat_conductance_megajoules_per_m2_h_kpa) or settings.latent_heat_conductance_megajoules_per_m2_h_kpa < 0 or !std.math.isFinite(settings.liquid_water_heat_capacity_megajoules_per_m3_k) or settings.liquid_water_heat_capacity_megajoules_per_m3_k <= 0 or !std.math.isFinite(settings.latent_heat_of_vaporization_megajoules_per_m3) or settings.latent_heat_of_vaporization_megajoules_per_m3 <= 0 or !std.math.isFinite(settings.surface_vapor_activity_fraction) or settings.surface_vapor_activity_fraction < 0 or settings.surface_vapor_activity_fraction > 1 or !std.math.isFinite(settings.timestep_hours) or settings.timestep_hours <= 0 or !std.math.isFinite(settings.minimum_temperature_k) or !std.math.isFinite(settings.maximum_temperature_k) or settings.minimum_temperature_k <= 0 or settings.minimum_temperature_k >= settings.maximum_temperature_k) return error.InvalidSurfaceTemperatureSettings;
}

test "hybrid Newton closes radiative sensible surface balance" {
    const context: ResidualContext = .{ .absorbed_shortwave_megajoules_per_m2 = 1.2, .atmospheric_longwave_megajoules_per_m2 = 0.9, .air_temperature_k = 285, .subsurface_temperature_k = 280, .previous_surface_temperature_k = 285, .ground_exposure_fraction = 1, .emissivity = 0.97, .sensible_heat_conductance_megajoules_per_m2_h_k = 0.43, .latent_heat_conductance_megajoules_per_m2_h_kpa = 0.1, .conductive_heat_conductance_megajoules_per_m2_h_k = 0.1, .storage_heat_conductance_megajoules_per_m2_h_k = 0.2, .atmospheric_vapor_pressure_kpa = 1.0, .surface_vapor_activity_fraction = 1.0 };
    const solved = try numerics.newtonPicard(context, residual, derivative, picard, 173.15, 373.15, 285, .{});
    try std.testing.expect(@abs(solved.residual) < 1.0e-7);
    try std.testing.expect(solved.root > 285);
    try std.testing.expect(solved.newton_raphson_steps > 0);
}

test "failed-pass reset removes retained phase carrier changes" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    @memset(state.liquid_water_change_m3, 1.25);
    @memset(state.ice_water_equivalent_change_m3, -1.25);
    @memset(state.phase_heat_flux_megajoules_per_m2, 416.25);

    state.resetPhaseChangeDiagnostics();

    for (state.liquid_water_change_m3) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
    for (state.ice_water_equivalent_change_m3) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
    for (state.phase_heat_flux_megajoules_per_m2) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
}

test "delayed litter combustion heat is closed inside the hybrid surface solve" {
    const baseline: ResidualContext = .{
        .absorbed_shortwave_megajoules_per_m2 = 0.8,
        .atmospheric_longwave_megajoules_per_m2 = 0.7,
        .air_temperature_k = 280,
        .subsurface_temperature_k = 278,
        .previous_surface_temperature_k = 281,
        .ground_exposure_fraction = 1,
        .emissivity = 0.97,
        .sensible_heat_conductance_megajoules_per_m2_h_k = 0.43,
        .latent_heat_conductance_megajoules_per_m2_h_kpa = 0.1,
        .conductive_heat_conductance_megajoules_per_m2_h_k = 0.1,
        .storage_heat_conductance_megajoules_per_m2_h_k = 0.2,
        .atmospheric_vapor_pressure_kpa = 0.9,
        .surface_vapor_activity_fraction = 0.8,
    };
    var heated = baseline;
    heated.external_heat_megajoules_per_m2 = 2;
    const baseline_solution = try numerics.newtonPicard(baseline, residual, derivative, picard, 173.15, 373.15, 281, .{});
    const heated_solution = try numerics.newtonPicard(heated, residual, derivative, picard, 173.15, 373.15, 281, .{});
    try std.testing.expect(heated_solution.root > baseline_solution.root);
    try std.testing.expect(@abs(heated_solution.residual) < 1e-7);
}

test "saturation vapor pressure follows the ecosys relation" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.61), saturationVaporPressureKpa(1.0 / 3.661e-3), 1.0e-14);
    try std.testing.expect(saturationVaporPressureKpa(300) > saturationVaporPressureKpa(280));
}

test "surface vapor sensible heat uses the phase donor temperature" {
    try std.testing.expectApproxEqAbs(
        0.002 * 4.19 * 285,
        try vaporSensibleHeatFlux(0.002 * 2465, 285, 280, 4.19, 2465),
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        -0.002 * 4.19 * 280,
        try vaporSensibleHeatFlux(-0.002 * 2465, 285, 280, 4.19, 2465),
        1.0e-12,
    );
}

test "analytic residual derivative agrees with finite difference" {
    const context: ResidualContext = .{ .absorbed_shortwave_megajoules_per_m2 = 0.8, .atmospheric_longwave_megajoules_per_m2 = 0.7, .air_temperature_k = 280, .subsurface_temperature_k = 278, .previous_surface_temperature_k = 281, .ground_exposure_fraction = 0.75, .emissivity = 0.97, .sensible_heat_conductance_megajoules_per_m2_h_k = 0.43, .latent_heat_conductance_megajoules_per_m2_h_kpa = 0.2, .conductive_heat_conductance_megajoules_per_m2_h_k = 0.1, .storage_heat_conductance_megajoules_per_m2_h_k = 0.2, .atmospheric_vapor_pressure_kpa = 0.9, .surface_vapor_activity_fraction = 0.85 };
    const temperature_k: f64 = 286;
    const step: f64 = 1.0e-4;
    const finite_difference = (residual(context, temperature_k + step) - residual(context, temperature_k - step)) / (2 * step);
    try std.testing.expectApproxEqRel(finite_difference, derivative(context, temperature_k), 1.0e-8);
}

test "surfacePhaseEquilibrium applies the WATSUB energy limit" {
    // Regression test for the wiring, not the limit arithmetic: reverting
    // `surfacePhaseEquilibrium` to return the raw Dall'Amico equilibrium must fail
    // here. It previously left every test passing. See EXEC-002.
    //
    // The measured Ottawa hour-2 state: 6.9736e4 m3 of litter water at 277.145 K
    // with a heat capacity of 4.19 x that volume. Unlimited, the equilibrium at a
    // below-freezing trial temperature freezes essentially the whole reservoir;
    // limited, it may only freeze what the temperature deficit supports.
    const liquid_m3 = 6.973601907996609e4;
    const phase: CellPhaseContext = .{
        .initial_liquid_water_m3 = liquid_m3,
        .initial_ice_water_equivalent_m3 = 0,
        .porous_medium_volume_m3 = liquid_m3,
        .horizontal_area_m2 = 8.717002384995762e7,
        .parameters = .{
            .residual_water_content_m3_per_m3 = 0,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 400,
            .n = 2.5,
            .saturated_hydraulic_conductivity_m_per_h = 0,
        },
        .gravitational_water_potential_mpa_per_m = 0.0098,
        .latent_heat_of_fusion_megajoules_per_m3 = 333,
        .pure_water_melting_temperature_k = 273.15,
        .heat_capacity_megajoules_per_k = 4.19 * liquid_m3,
    };

    // A trial temperature 4 K below freezing, as the hour-2 solve found.
    const limited = try surfacePhaseEquilibrium(phase, 273.15 - 4.0);

    // The energy the deficit can drive is capacity x deficit, so the ice formed
    // must not exceed that divided by the latent heat. With capacity 4.19*V and a
    // 4 K deficit that is about 0.0503 of the reservoir, far short of all of it.
    const energy_limited_ice_m3 = phase.heat_capacity_megajoules_per_k * 4.0 / 333.0;
    try std.testing.expect(limited.ice_water_equivalent_m3 <= energy_limited_ice_m3 * 1.05);
    // And it must be a small fraction of the reservoir, which is the property the
    // unlimited equilibrium violates.
    try std.testing.expect(limited.ice_water_equivalent_m3 < 0.25 * liquid_m3);
    // Water equivalent is still conserved exactly.
    try std.testing.expectApproxEqRel(
        liquid_m3,
        limited.liquid_water_m3 + limited.ice_water_equivalent_m3,
        1e-14,
    );

    // Above the freezing point nothing freezes, matching the oracle's TFREEZ gate.
    const warm = try surfacePhaseEquilibrium(phase, 277.145);
    try std.testing.expectEqual(@as(f64, 0), warm.ice_water_equivalent_m3);
}

test "surface residue Dall'Amico enthalpy converges without a sub-hour cycle" {
    const phase: CellPhaseContext = .{
        .initial_liquid_water_m3 = 0.01,
        .initial_ice_water_equivalent_m3 = 0,
        .porous_medium_volume_m3 = 0.02,
        .horizontal_area_m2 = 1,
        .parameters = .{
            .residual_water_content_m3_per_m3 = 0,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 400,
            .n = 2.5,
            .saturated_hydraulic_conductivity_m_per_h = 0,
        },
        .gravitational_water_potential_mpa_per_m = 0.0098,
        .latent_heat_of_fusion_megajoules_per_m3 = 333,
        .pure_water_melting_temperature_k = 273.15,
        // Ample capacity so the WATSUB energy limit does not bind in this test.
        .heat_capacity_megajoules_per_k = 1.0e9,
    };
    const context: ResidualContext = .{
        .absorbed_shortwave_megajoules_per_m2 = 0,
        .atmospheric_longwave_megajoules_per_m2 = 0,
        .air_temperature_k = 260,
        .subsurface_temperature_k = 260,
        .previous_surface_temperature_k = 260,
        .ground_exposure_fraction = 0,
        .emissivity = 0,
        .sensible_heat_conductance_megajoules_per_m2_h_k = 1,
        .latent_heat_conductance_megajoules_per_m2_h_kpa = 0,
        .conductive_heat_conductance_megajoules_per_m2_h_k = 1,
        .storage_heat_conductance_megajoules_per_m2_h_k = 1,
        .atmospheric_vapor_pressure_kpa = 0,
        .surface_vapor_activity_fraction = 0,
        .phase = phase,
    };
    const solved = try numerics.newtonPicard(
        context,
        residual,
        derivative,
        picard,
        240,
        300,
        260,
        .{ .max_iterations = 80 },
    );
    const equilibrium = try surfacePhaseEquilibrium(phase, solved.root);
    try std.testing.expect(solved.iterations < 80);
    try std.testing.expect(@abs(solved.residual) < 1e-7);
    try std.testing.expect(equilibrium.ice_water_equivalent_m3 > 0);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.01),
        equilibrium.liquid_water_m3 +
            equilibrium.ice_water_equivalent_m3,
        1e-14,
    );
    try std.testing.expect(solved.root > 260);
}

test "surface Dall'Amico phase derivative matches the frozen branch" {
    const phase: CellPhaseContext = .{
        .initial_liquid_water_m3 = 0.01,
        .initial_ice_water_equivalent_m3 = 0,
        .porous_medium_volume_m3 = 0.02,
        .horizontal_area_m2 = 1,
        .parameters = .{
            .residual_water_content_m3_per_m3 = 0,
            .saturated_water_content_m3_per_m3 = 1,
            .alpha_per_m = 400,
            .n = 2.5,
            .saturated_hydraulic_conductivity_m_per_h = 0,
        },
        .gravitational_water_potential_mpa_per_m = 0.0098,
        .latent_heat_of_fusion_megajoules_per_m3 = 333,
        .pure_water_melting_temperature_k = 273.15,
        // Ample capacity so the WATSUB energy limit does not bind in this test.
        .heat_capacity_megajoules_per_k = 1.0e9,
    };
    const context: ResidualContext = .{
        .absorbed_shortwave_megajoules_per_m2 = 0,
        .atmospheric_longwave_megajoules_per_m2 = 0,
        .air_temperature_k = 260,
        .subsurface_temperature_k = 260,
        .previous_surface_temperature_k = 260,
        .ground_exposure_fraction = 0,
        .emissivity = 0,
        .sensible_heat_conductance_megajoules_per_m2_h_k = 1,
        .latent_heat_conductance_megajoules_per_m2_h_kpa = 0,
        .conductive_heat_conductance_megajoules_per_m2_h_k = 1,
        .storage_heat_conductance_megajoules_per_m2_h_k = 1,
        .atmospheric_vapor_pressure_kpa = 0,
        .surface_vapor_activity_fraction = 0,
        .phase = phase,
    };
    const temperature_k: f64 = 268;
    const step_k: f64 = 1.0e-4;
    const finite_difference =
        (phaseHeatFlux(context, temperature_k + step_k) -
            phaseHeatFlux(context, temperature_k - step_k)) / (2 * step_k);
    try std.testing.expectApproxEqRel(
        finite_difference,
        phaseHeatDerivative(context, temperature_k),
        2.0e-4,
    );
}
