const std = @import("std");
const sediment_routing = @import("sediment_routing.zig");
const numerics = @import("numerics.zig");

pub const RuntimeState = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    surface_sediment_Mg: []f64,
    local_detachment_Mg: []f64,
    transportable_sediment_Mg: []f64,
    /// REDIST TSEDSK: mineral sediment settled from the pond this hour (Mg).
    pond_settled_sediment_Mg: []f64,
    surface_soil_mass_Mg: []f64,
    surface_soil_mass_initialized: []bool,
    east_boundary_open: []bool,
    west_boundary_open: []bool,
    south_boundary_open: []bool,
    north_boundary_open: []bool,
    routing: sediment_routing.RoutingState,

    pub fn init(allocator: std.mem.Allocator, columns: usize, rows: usize) !RuntimeState {
        const count = try std.math.mul(usize, columns, rows);
        if (count == 0) return error.InvalidSedimentGridDimensions;
        var result: RuntimeState = undefined;
        result.allocator = allocator;
        result.cell_count = count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(RuntimeState).@"struct".fields) |field| if ((field.type == []f64 or field.type == []bool) and allocated > 0) {
                allocated -= 1;
                allocator.free(@field(result, field.name));
            };
        }
        inline for (@typeInfo(RuntimeState).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(result, field.name) = try allocator.alloc(f64, count);
                @memset(@field(result, field.name), 0);
                allocated += 1;
            } else if (field.type == []bool) {
                @field(result, field.name) = try allocator.alloc(bool, count);
                @memset(@field(result, field.name), false);
                allocated += 1;
            }
        }
        result.routing = try sediment_routing.RoutingState.init(allocator, columns, rows);
        return result;
    }

    pub fn deinit(self: *RuntimeState) void {
        self.routing.deinit();
        inline for (@typeInfo(RuntimeState).@"struct".fields) |field| if (field.type == []f64 or field.type == []bool) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

pub const SurfacePropertyParameters = struct {
    clay_particle_diameter_um: f64 = 1,
    silt_particle_diameter_um: f64 = 10,
    sand_particle_diameter_um: f64 = 100,
    humus_particle_diameter_um: f64 = 10,
    residue_particle_diameter_um: f64 = 100,
    mineral_particle_density_Mg_per_m3: f64 = 2.66,
    organic_particle_density_Mg_per_m3: f64 = 1.30,
    reference_water_viscosity_Mg_per_m_s: f64,
    viscosity_temperature_intercept: f64,
    viscosity_temperature_coefficient_per_c: f64,
    gravitational_acceleration_m_per_s2: f64 = 9.8,
};

pub const SurfacePropertyInputs = struct {
    sand_mass_fraction: f64,
    silt_mass_fraction: f64,
    clay_mass_fraction: f64,
    humus_mass_fraction: f64,
    residue_mass_fraction: f64,
    root_length_density_m_per_m3: f64,
    surface_temperature_c: f64,
};

pub const SurfaceProperties = struct {
    mean_particle_diameter_um: f64,
    rainfall_detachability_g_per_j: f64,
    runoff_detachability: f64,
    particle_density_Mg_per_m3: f64,
    settling_velocity_m_per_h: f64,
    transport_capacity_coefficient: f64,
    transport_capacity_exponent: f64,
};

/// HOUR1 surface cohesion and EUROSEM erosion properties. These are refreshed
/// from runtime soil and organic state rather than frozen compile-time values.
pub fn deriveSurfaceProperties(inputs: SurfacePropertyInputs, parameters: SurfacePropertyParameters) !SurfaceProperties {
    inline for (@typeInfo(SurfacePropertyInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteErosionInput;
    inline for (@typeInfo(SurfacePropertyParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteErosionInput;
    if (inputs.sand_mass_fraction < 0 or inputs.silt_mass_fraction < 0 or inputs.clay_mass_fraction < 0 or inputs.humus_mass_fraction < 0 or inputs.residue_mass_fraction < 0 or inputs.root_length_density_m_per_m3 < 0 or parameters.reference_water_viscosity_Mg_per_m_s <= 0 or parameters.mineral_particle_density_Mg_per_m3 <= 1 or parameters.organic_particle_density_Mg_per_m3 <= 1 or parameters.gravitational_acceleration_m_per_s2 <= 0) return error.InvalidErosionInput;
    const diameter_um =
        parameters.clay_particle_diameter_um * inputs.clay_mass_fraction +
        parameters.silt_particle_diameter_um * inputs.silt_mass_fraction +
        parameters.sand_particle_diameter_um * inputs.sand_mass_fraction +
        parameters.humus_particle_diameter_um * inputs.humus_mass_fraction +
        parameters.residue_particle_diameter_um * inputs.residue_mass_fraction;
    const organic_fraction = std.math.clamp(inputs.humus_mass_fraction + inputs.residue_mass_fraction, 0, 1);
    const cohesion_kpa = 2 + 5 * inputs.residue_mass_fraction + 5 * inputs.humus_mass_fraction + 5 * inputs.clay_mass_fraction + 5 * (1 - @exp(-1.0e-5 * inputs.root_length_density_m_per_m3));
    const particle_density = parameters.organic_particle_density_Mg_per_m3 * organic_fraction + parameters.mineral_particle_density_Mg_per_m3 * (1 - organic_fraction);
    const viscosity = parameters.reference_water_viscosity_Mg_per_m_s * @exp(parameters.viscosity_temperature_intercept - parameters.viscosity_temperature_coefficient_per_c * inputs.surface_temperature_c);
    const settling_velocity = 3.6e3 * parameters.gravitational_acceleration_m_per_s2 * (particle_density - 1) * std.math.pow(f64, 1.0e-6 * diameter_um, 2) / (18 * viscosity);
    const result: SurfaceProperties = .{
        .mean_particle_diameter_um = diameter_um,
        .rainfall_detachability_g_per_j = 1.0e-6 * (1.5 + 2.5 * (1 - inputs.silt_mass_fraction - inputs.residue_mass_fraction)),
        .runoff_detachability = 0.79 * @exp(-0.85 * cohesion_kpa),
        .particle_density_Mg_per_m3 = particle_density,
        .settling_velocity_m_per_h = settling_velocity,
        .transport_capacity_coefficient = std.math.pow(f64, (diameter_um + 5) / 0.32, -0.6),
        .transport_capacity_exponent = std.math.pow(f64, (diameter_um + 5) / 300, 0.25),
    };
    inline for (@typeInfo(SurfaceProperties).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.NonFiniteErosionResult;
    return result;
}

pub const StepInputs = struct {
    erosion_enabled: bool,
    surface_soil_bulk_density_Mg_per_m3: f64,
    surface_soil_mass_Mg: f64,
    surface_soil_water_m3: f64,
    surface_soil_pore_volume_m3: f64,
    excess_surface_water_m3: f64,
    excess_surface_ice_m3: f64,
    surface_ponding_capacity_m3: f64,
    sediment_in_surface_water_Mg: f64,
    rainfall_kinetic_energy_j: f64,
    soil_rainfall_detachability_g_per_j: f64,
    soil_runoff_detachability: f64,
    sediment_settling_velocity_m_per_h: f64,
    grid_cell_area_m2: f64,
    soil_matrix_fraction: f64,
    snow_free_fraction: f64,
    runoff_velocity_m_per_s: f64,
    slope_sine: f64,
    surface_particle_density_Mg_per_m3: f64,
    transport_capacity_coefficient: f64,
    transport_capacity_exponent: f64,
    maximum_erodible_soil_fraction_per_step: f64,
    water_transport_timestep_h: f64,
    negligible_volume_m3: f64,
    negligible_mass_Mg: f64,
};

pub const StepResult = struct {
    rainfall_detachment_Mg: f64,
    immobile_water_deposition_Mg: f64,
    runoff_detachment_or_deposition_Mg: f64,
    net_detachment_Mg: f64,
    runoff_sediment_capacity_Mg_per_m3: f64,
    surface_sediment_concentration_Mg_per_m3: f64,
};

pub const HourlySolveOptions = struct {
    absolute_tolerance_Mg: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
    /// The runtime NPH input is a convergence ceiling, not a sub-hour loop.
    max_iterations: u16,
};

pub const HourlySolveResult = struct {
    local: StepResult,
    final_surface_sediment_Mg: f64,
    iterations: u16,
    newton_raphson_steps: u16,
    picard_steps: u16,
};

/// EROSION RERSED0: sediment available to downslope runoff after local
/// detachment/deposition, limited by both available mass and runoff volume.
pub fn calculateDownslopeTransport_Mg(surface_sediment_Mg: f64, net_local_detachment_Mg: f64, excess_surface_water_m3: f64, downslope_runoff_m3: f64, surface_soil_bulk_density_Mg_per_m3: f64, negligible_volume_m3: f64) !f64 {
    const values = [_]f64{ surface_sediment_Mg, net_local_detachment_Mg, excess_surface_water_m3, downslope_runoff_m3, surface_soil_bulk_density_Mg_per_m3, negligible_volume_m3 };
    for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteErosionInput;
    if (surface_sediment_Mg < 0 or excess_surface_water_m3 < 0 or downslope_runoff_m3 < 0 or surface_soil_bulk_density_Mg_per_m3 < 0 or negligible_volume_m3 < 0) return error.InvalidErosionInput;
    if (downslope_runoff_m3 <= 0 or surface_soil_bulk_density_Mg_per_m3 <= 0 or excess_surface_water_m3 <= negligible_volume_m3) return 0;
    const available_sediment_Mg = @max(0, surface_sediment_Mg + net_local_detachment_Mg);
    const concentration_Mg_per_m3 = available_sediment_Mg / excess_surface_water_m3;
    return @min(available_sediment_Mg, concentration_Mg_per_m3 * downslope_runoff_m3);
}

/// Ports EROSION's cell-local source terms. It is called once for each
/// converged hydrology step; the former NPH full-cycle repetition is absent.
pub fn calculateLocalStep(inputs: StepInputs) !StepResult {
    try validate(inputs);
    var result: StepResult = .{ .rainfall_detachment_Mg = 0, .immobile_water_deposition_Mg = 0, .runoff_detachment_or_deposition_Mg = 0, .net_detachment_Mg = 0, .runoff_sediment_capacity_Mg_per_m3 = 0, .surface_sediment_concentration_Mg_per_m3 = 0 };
    if (!inputs.erosion_enabled) return result;

    const water_fraction = std.math.clamp(inputs.excess_surface_water_m3 / inputs.surface_ponding_capacity_m3, 0, 1);
    const ice_fraction = std.math.clamp(inputs.excess_surface_ice_m3 / inputs.surface_ponding_capacity_m3, 0, 1);
    const mobile_surface_fraction = (1 - ice_fraction) * water_fraction;
    const soil_mass_limit_Mg = inputs.surface_soil_mass_Mg * inputs.maximum_erodible_soil_fraction_per_step;

    if (inputs.surface_soil_bulk_density_Mg_per_m3 > 0 and inputs.rainfall_kinetic_energy_j > 0 and inputs.excess_surface_water_m3 > inputs.negligible_volume_m3) {
        const water_attenuated_detachability_g_per_j = inputs.soil_rainfall_detachability_g_per_j * (1 + 2 * inputs.surface_soil_water_m3 / inputs.surface_soil_pore_volume_m3);
        result.rainfall_detachment_Mg = @min(soil_mass_limit_Mg, water_attenuated_detachability_g_per_j * inputs.rainfall_kinetic_energy_j * inputs.grid_cell_area_m2 * inputs.soil_matrix_fraction * inputs.snow_free_fraction * (1 - ice_fraction));
    }

    if (inputs.surface_soil_bulk_density_Mg_per_m3 > 0 and inputs.sediment_in_surface_water_Mg > inputs.negligible_mass_Mg and inputs.excess_surface_water_m3 > inputs.negligible_volume_m3 and mobile_surface_fraction > 0) {
        result.surface_sediment_concentration_Mg_per_m3 = @max(0, inputs.sediment_in_surface_water_Mg / inputs.excess_surface_water_m3);
        result.immobile_water_deposition_Mg = @max(-inputs.sediment_in_surface_water_Mg, -inputs.sediment_settling_velocity_m_per_h * result.surface_sediment_concentration_Mg_per_m3 * inputs.grid_cell_area_m2 * mobile_surface_fraction * inputs.soil_matrix_fraction * inputs.water_transport_timestep_h);
    }

    if (inputs.surface_soil_bulk_density_Mg_per_m3 > 0 and inputs.excess_surface_water_m3 > inputs.surface_ponding_capacity_m3 and mobile_surface_fraction > 0) {
        const stream_power_cm_per_s = 100 * inputs.runoff_velocity_m_per_s * @abs(inputs.slope_sine);
        result.runoff_sediment_capacity_Mg_per_m3 = inputs.surface_particle_density_Mg_per_m3 * inputs.transport_capacity_coefficient * std.math.pow(f64, @max(0, stream_power_cm_per_s - 0.4), inputs.transport_capacity_exponent);
        result.surface_sediment_concentration_Mg_per_m3 = @max(0, inputs.sediment_in_surface_water_Mg / inputs.excess_surface_water_m3);
        const concentration_difference = result.runoff_sediment_capacity_Mg_per_m3 - result.surface_sediment_concentration_Mg_per_m3;
        const potential_Mg = inputs.sediment_settling_velocity_m_per_h * concentration_difference * inputs.grid_cell_area_m2 * mobile_surface_fraction * inputs.soil_matrix_fraction * inputs.water_transport_timestep_h;
        if (concentration_difference > 0) {
            result.runoff_detachment_or_deposition_Mg = @min(soil_mass_limit_Mg, inputs.soil_runoff_detachability * potential_Mg);
        } else if (inputs.sediment_in_surface_water_Mg > inputs.negligible_mass_Mg) {
            result.runoff_detachment_or_deposition_Mg = @max(-inputs.sediment_in_surface_water_Mg, potential_Mg);
        }
    }
    result.net_detachment_Mg = @max(
        -inputs.sediment_in_surface_water_Mg,
        result.rainfall_detachment_Mg + result.immobile_water_deposition_Mg + result.runoff_detachment_or_deposition_Mg,
    );
    if (!std.math.isFinite(result.net_detachment_Mg)) return error.NonFiniteErosionResult;
    return result;
}

const HourlyBalance = struct {
    inputs: StepInputs,
    initial_sediment_Mg: f64,
    rainfall_detachment_Mg: f64,
    upper_sediment_Mg: f64,

    fn source(self: HourlyBalance, suspended_sediment_Mg: f64) f64 {
        const flux = concentrationDependentFluxes(
            self.inputs,
            suspended_sediment_Mg,
        );
        return self.rainfall_detachment_Mg +
            flux.immobile_water_deposition_Mg +
            flux.runoff_detachment_or_deposition_Mg;
    }
};

fn hourlyResidual(balance: HourlyBalance, suspended_sediment_Mg: f64) f64 {
    return suspended_sediment_Mg -
        balance.initial_sediment_Mg -
        balance.source(suspended_sediment_Mg);
}

fn hourlyPicard(balance: HourlyBalance, suspended_sediment_Mg: f64) f64 {
    return std.math.clamp(
        balance.initial_sediment_Mg + balance.source(suspended_sediment_Mg),
        0,
        balance.upper_sediment_Mg,
    );
}

const ConcentrationFluxes = struct {
    immobile_water_deposition_Mg: f64,
    runoff_detachment_or_deposition_Mg: f64,
    runoff_sediment_capacity_Mg_per_m3: f64,
    surface_sediment_concentration_Mg_per_m3: f64,
};

fn concentrationDependentFluxes(
    inputs: StepInputs,
    suspended_sediment_Mg: f64,
) ConcentrationFluxes {
    var result: ConcentrationFluxes = .{
        .immobile_water_deposition_Mg = 0,
        .runoff_detachment_or_deposition_Mg = 0,
        .runoff_sediment_capacity_Mg_per_m3 = 0,
        .surface_sediment_concentration_Mg_per_m3 = 0,
    };
    const water_fraction = std.math.clamp(
        inputs.excess_surface_water_m3 / inputs.surface_ponding_capacity_m3,
        0,
        1,
    );
    const ice_fraction = std.math.clamp(
        inputs.excess_surface_ice_m3 / inputs.surface_ponding_capacity_m3,
        0,
        1,
    );
    const mobile_surface_fraction = (1 - ice_fraction) * water_fraction;
    const soil_mass_limit_Mg =
        inputs.surface_soil_mass_Mg *
        inputs.maximum_erodible_soil_fraction_per_step;
    if (inputs.surface_soil_bulk_density_Mg_per_m3 > 0 and
        suspended_sediment_Mg > inputs.negligible_mass_Mg and
        inputs.excess_surface_water_m3 > inputs.negligible_volume_m3 and
        mobile_surface_fraction > 0)
    {
        result.surface_sediment_concentration_Mg_per_m3 =
            suspended_sediment_Mg / inputs.excess_surface_water_m3;
        result.immobile_water_deposition_Mg = @max(
            -suspended_sediment_Mg,
            -inputs.sediment_settling_velocity_m_per_h *
                result.surface_sediment_concentration_Mg_per_m3 *
                inputs.grid_cell_area_m2 *
                mobile_surface_fraction *
                inputs.soil_matrix_fraction *
                inputs.water_transport_timestep_h,
        );
    }
    if (inputs.surface_soil_bulk_density_Mg_per_m3 > 0 and
        inputs.excess_surface_water_m3 > inputs.surface_ponding_capacity_m3 and
        mobile_surface_fraction > 0)
    {
        const stream_power_cm_per_s =
            100 * inputs.runoff_velocity_m_per_s * @abs(inputs.slope_sine);
        result.runoff_sediment_capacity_Mg_per_m3 =
            inputs.surface_particle_density_Mg_per_m3 *
            inputs.transport_capacity_coefficient *
            std.math.pow(
                f64,
                @max(0, stream_power_cm_per_s - 0.4),
                inputs.transport_capacity_exponent,
            );
        result.surface_sediment_concentration_Mg_per_m3 =
            suspended_sediment_Mg / inputs.excess_surface_water_m3;
        const concentration_difference =
            result.runoff_sediment_capacity_Mg_per_m3 -
            result.surface_sediment_concentration_Mg_per_m3;
        const potential_Mg =
            inputs.sediment_settling_velocity_m_per_h *
            concentration_difference *
            inputs.grid_cell_area_m2 *
            mobile_surface_fraction *
            inputs.soil_matrix_fraction *
            inputs.water_transport_timestep_h;
        if (concentration_difference > 0) {
            result.runoff_detachment_or_deposition_Mg =
                @min(soil_mass_limit_Mg, inputs.soil_runoff_detachability * potential_Mg);
        } else if (suspended_sediment_Mg > inputs.negligible_mass_Mg) {
            result.runoff_detachment_or_deposition_Mg =
                @max(-suspended_sediment_Mg, potential_Mg);
        }
    }
    return result;
}

/// Replaces EROSION's repeated NPH mutation of `SED` with one bounded
/// backward-Euler balance. Newton–Raphson is attempted first and Picard is
/// the safeguarded fallback. The caller-provided NPH is only the maximum
/// number of nonlinear iterations and the solve exits as soon as its
/// extensive sediment residual converges.
pub fn calculateConvergedHourlyLocalStep(
    inputs: StepInputs,
    options: HourlySolveOptions,
) !HourlySolveResult {
    try validate(inputs);
    if (!std.math.isFinite(options.absolute_tolerance_Mg) or
        options.absolute_tolerance_Mg <= 0 or
        !std.math.isFinite(options.relative_tolerance) or
        options.relative_tolerance <= 0 or
        !std.math.isFinite(options.picard_relaxation) or
        options.picard_relaxation <= 0 or
        options.picard_relaxation > 1 or
        options.max_iterations == 0)
        return error.InvalidErosionSolverOptions;
    if (!inputs.erosion_enabled) return .{
        .local = .{
            .rainfall_detachment_Mg = 0,
            .immobile_water_deposition_Mg = 0,
            .runoff_detachment_or_deposition_Mg = 0,
            .net_detachment_Mg = 0,
            .runoff_sediment_capacity_Mg_per_m3 = 0,
            .surface_sediment_concentration_Mg_per_m3 = 0,
        },
        .final_surface_sediment_Mg = inputs.sediment_in_surface_water_Mg,
        .iterations = 0,
        .newton_raphson_steps = 0,
        .picard_steps = 0,
    };

    var rainfall_only = inputs;
    rainfall_only.sediment_in_surface_water_Mg = 0;
    rainfall_only.sediment_settling_velocity_m_per_h = 0;
    const rainfall_detachment_Mg =
        (try calculateLocalStep(rainfall_only)).rainfall_detachment_Mg;
    const soil_mass_limit_Mg =
        inputs.surface_soil_mass_Mg *
        inputs.maximum_erodible_soil_fraction_per_step;
    const upper_sediment_Mg =
        inputs.sediment_in_surface_water_Mg +
        rainfall_detachment_Mg +
        soil_mass_limit_Mg +
        @max(options.absolute_tolerance_Mg, std.math.floatEps(f64));
    if (upper_sediment_Mg <= options.absolute_tolerance_Mg) return .{
        .local = .{
            .rainfall_detachment_Mg = rainfall_detachment_Mg,
            .immobile_water_deposition_Mg = 0,
            .runoff_detachment_or_deposition_Mg = 0,
            .net_detachment_Mg = rainfall_detachment_Mg,
            .runoff_sediment_capacity_Mg_per_m3 = 0,
            .surface_sediment_concentration_Mg_per_m3 = 0,
        },
        .final_surface_sediment_Mg = rainfall_detachment_Mg,
        .iterations = 0,
        .newton_raphson_steps = 0,
        .picard_steps = 0,
    };
    const balance: HourlyBalance = .{
        .inputs = inputs,
        .initial_sediment_Mg = inputs.sediment_in_surface_water_Mg,
        .rainfall_detachment_Mg = rainfall_detachment_Mg,
        .upper_sediment_Mg = upper_sediment_Mg,
    };
    const solved = try numerics.newtonPicardFiniteDifference(
        balance,
        hourlyResidual,
        hourlyPicard,
        0,
        upper_sediment_Mg,
        std.math.clamp(
            inputs.sediment_in_surface_water_Mg + rainfall_detachment_Mg,
            0,
            upper_sediment_Mg,
        ),
        .{
            .absolute_tolerance = options.absolute_tolerance_Mg,
            .relative_tolerance = options.relative_tolerance,
            .picard_relaxation = options.picard_relaxation,
            .residual_scale = @max(
                options.absolute_tolerance_Mg,
                inputs.sediment_in_surface_water_Mg + rainfall_detachment_Mg,
            ),
            .max_iterations = options.max_iterations,
            .safeguard_with_bracket = true,
        },
    );
    const flux = concentrationDependentFluxes(inputs, solved.root);
    const net_detachment_Mg =
        solved.root - inputs.sediment_in_surface_water_Mg;
    if (!std.math.isFinite(net_detachment_Mg) or solved.root < 0)
        return error.NonFiniteErosionResult;
    return .{
        .local = .{
            .rainfall_detachment_Mg = rainfall_detachment_Mg,
            .immobile_water_deposition_Mg = flux.immobile_water_deposition_Mg,
            // Close the accepted extensive balance exactly; the residual is
            // below tolerance but must not enter the authoritative ledger.
            .runoff_detachment_or_deposition_Mg = net_detachment_Mg -
                rainfall_detachment_Mg -
                flux.immobile_water_deposition_Mg,
            .net_detachment_Mg = net_detachment_Mg,
            .runoff_sediment_capacity_Mg_per_m3 = flux.runoff_sediment_capacity_Mg_per_m3,
            .surface_sediment_concentration_Mg_per_m3 = flux.surface_sediment_concentration_Mg_per_m3,
        },
        .final_surface_sediment_Mg = solved.root,
        .iterations = solved.iterations,
        .newton_raphson_steps = solved.newton_raphson_steps,
        .picard_steps = solved.picard_steps,
    };
}

fn validate(inputs: StepInputs) !void {
    inline for (@typeInfo(StepInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteErosionInput;
    if (inputs.surface_soil_bulk_density_Mg_per_m3 < 0 or inputs.surface_soil_mass_Mg < 0 or inputs.surface_soil_water_m3 < 0 or inputs.surface_soil_pore_volume_m3 <= 0 or inputs.excess_surface_water_m3 < 0 or inputs.excess_surface_ice_m3 < 0 or inputs.surface_ponding_capacity_m3 <= 0 or inputs.sediment_in_surface_water_Mg < 0 or inputs.rainfall_kinetic_energy_j < 0 or inputs.soil_rainfall_detachability_g_per_j < 0 or inputs.soil_runoff_detachability < 0 or inputs.sediment_settling_velocity_m_per_h < 0 or inputs.grid_cell_area_m2 <= 0 or inputs.soil_matrix_fraction < 0 or inputs.soil_matrix_fraction > 1 or inputs.snow_free_fraction < 0 or inputs.snow_free_fraction > 1 or inputs.surface_particle_density_Mg_per_m3 < 0 or inputs.transport_capacity_coefficient < 0 or inputs.transport_capacity_exponent < 0 or inputs.maximum_erodible_soil_fraction_per_step < 0 or inputs.maximum_erodible_soil_fraction_per_step > 1 or inputs.water_transport_timestep_h <= 0 or inputs.negligible_volume_m3 < 0 or inputs.negligible_mass_Mg < 0) return error.InvalidErosionInput;
}

fn baseline() StepInputs {
    return .{ .erosion_enabled = true, .surface_soil_bulk_density_Mg_per_m3 = 1.2, .surface_soil_mass_Mg = 100, .surface_soil_water_m3 = 0.2, .surface_soil_pore_volume_m3 = 1, .excess_surface_water_m3 = 2, .excess_surface_ice_m3 = 0, .surface_ponding_capacity_m3 = 1, .sediment_in_surface_water_Mg = 0.2, .rainfall_kinetic_energy_j = 1, .soil_rainfall_detachability_g_per_j = 0.01, .soil_runoff_detachability = 0.5, .sediment_settling_velocity_m_per_h = 0.01, .grid_cell_area_m2 = 10, .soil_matrix_fraction = 0.8, .snow_free_fraction = 1, .runoff_velocity_m_per_s = 0.01, .slope_sine = 0.1, .surface_particle_density_Mg_per_m3 = 2.65, .transport_capacity_coefficient = 1, .transport_capacity_exponent = 1, .maximum_erodible_soil_fraction_per_step = 0.01, .water_transport_timestep_h = 0.25, .negligible_volume_m3 = 1e-12, .negligible_mass_Mg = 1e-12 };
}

test "rainfall detachment preserves EROSION water attenuation and soil mass cap" {
    var inputs = baseline();
    inputs.rainfall_kinetic_energy_j = 100;
    const result = try calculateLocalStep(inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 1), result.rainfall_detachment_Mg, 1e-12);
}

test "settling never deposits more sediment than surface water contains" {
    var inputs = baseline();
    inputs.excess_surface_water_m3 = 0.5;
    inputs.sediment_settling_velocity_m_per_h = 100;
    inputs.rainfall_kinetic_energy_j = 0;
    const result = try calculateLocalStep(inputs);
    try std.testing.expectEqual(@as(f64, -0.2), result.immobile_water_deposition_Mg);
}

test "combined deposition mechanisms cannot exceed suspended sediment" {
    var inputs = baseline();
    inputs.rainfall_kinetic_energy_j = 0;
    inputs.sediment_settling_velocity_m_per_h = 100;
    const result = try calculateLocalStep(inputs);
    try std.testing.expect(result.net_detachment_Mg >= -inputs.sediment_in_surface_water_Mg);
}

test "runoff below stream-power threshold has zero carrying capacity" {
    var inputs = baseline();
    inputs.runoff_velocity_m_per_s = 0.001;
    const result = try calculateLocalStep(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.runoff_sediment_capacity_Mg_per_m3);
    try std.testing.expect(result.runoff_detachment_or_deposition_Mg <= 0);
}

test "downslope sediment transport is limited by available sediment" {
    try std.testing.expectEqual(@as(f64, 0.6), try calculateDownslopeTransport_Mg(0.5, 0.1, 1, 2, 1.2, 1e-12));
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), try calculateDownslopeTransport_Mg(0.5, 0.1, 2, 0.5, 1.2, 1e-12), 1e-14);
}

test "HOUR1 erosion properties are derived from runtime texture and temperature" {
    const properties = try deriveSurfaceProperties(.{ .sand_mass_fraction = 0.6, .silt_mass_fraction = 0.25, .clay_mass_fraction = 0.15, .humus_mass_fraction = 0, .residue_mass_fraction = 0, .root_length_density_m_per_m3 = 0, .surface_temperature_c = 20 }, .{ .reference_water_viscosity_Mg_per_m_s = 1.0e-3, .viscosity_temperature_intercept = 0.533, .viscosity_temperature_coefficient_per_c = 0.0267 });
    try std.testing.expectApproxEqAbs(@as(f64, 62.65), properties.mean_particle_diameter_um, 1e-12);
    try std.testing.expect(properties.rainfall_detachability_g_per_j > 0);
    try std.testing.expect(properties.settling_velocity_m_per_h > 0);
}

test "hourly erosion converges before NPH without explicit sub-hour mutation" {
    var inputs = baseline();
    inputs.water_transport_timestep_h = 1;
    // The implicit solve represents the complete hour, so the NPH sequence
    // of per-pass 1/NPH caps becomes one full-hour cap.
    inputs.maximum_erodible_soil_fraction_per_step = 1;
    inputs.sediment_settling_velocity_m_per_h = 0.2;
    const solved = try calculateConvergedHourlyLocalStep(inputs, .{
        .absolute_tolerance_Mg = 1.0e-12,
        .relative_tolerance = 1.0e-10,
        .picard_relaxation = 0.5,
        .max_iterations = 20,
    });
    try std.testing.expect(solved.iterations < 20);
    try std.testing.expect(solved.newton_raphson_steps + solved.picard_steps > 0);
    try std.testing.expect(solved.final_surface_sediment_Mg >= 0);
    try std.testing.expectApproxEqAbs(
        solved.final_surface_sediment_Mg,
        inputs.sediment_in_surface_water_Mg + solved.local.net_detachment_Mg,
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        solved.local.net_detachment_Mg,
        solved.local.rainfall_detachment_Mg +
            solved.local.immobile_water_deposition_Mg +
            solved.local.runoff_detachment_or_deposition_Mg,
        1.0e-14,
    );
}

test "implicit erosion remains nonnegative under stiff settling" {
    var inputs = baseline();
    inputs.rainfall_kinetic_energy_j = 0;
    inputs.sediment_settling_velocity_m_per_h = 1.0e6;
    inputs.water_transport_timestep_h = 1;
    const solved = try calculateConvergedHourlyLocalStep(inputs, .{
        .absolute_tolerance_Mg = 1.0e-11,
        .relative_tolerance = 1.0e-8,
        .picard_relaxation = 0.5,
        .max_iterations = 20,
    });
    try std.testing.expect(solved.final_surface_sediment_Mg >= 0);
    try std.testing.expect(solved.final_surface_sediment_Mg <= inputs.sediment_in_surface_water_Mg);
}

test "disabled converged erosion returns unchanged sediment without iterations" {
    var inputs = baseline();
    inputs.erosion_enabled = false;
    const solved = try calculateConvergedHourlyLocalStep(inputs, .{
        .absolute_tolerance_Mg = 1.0e-12,
        .relative_tolerance = 1.0e-8,
        .picard_relaxation = 0.5,
        .max_iterations = 20,
    });
    try std.testing.expectEqual(@as(u16, 0), solved.iterations);
    try std.testing.expectEqual(inputs.sediment_in_surface_water_Mg, solved.final_surface_sediment_Mg);
    try std.testing.expectEqual(@as(f64, 0), solved.local.net_detachment_Mg);
}

test "invalid hourly erosion solve cannot mutate caller state" {
    const inputs = baseline();
    const before = inputs;
    try std.testing.expectError(
        error.InvalidErosionSolverOptions,
        calculateConvergedHourlyLocalStep(inputs, .{
            .absolute_tolerance_Mg = 1.0e-12,
            .relative_tolerance = 1.0e-8,
            .picard_relaxation = 0.5,
            .max_iterations = 0,
        }),
    );
    try std.testing.expectEqual(before, inputs);
}
