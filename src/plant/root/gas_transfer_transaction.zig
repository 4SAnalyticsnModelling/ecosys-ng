const std = @import("std");
const oxygen = @import("oxygen_uptake_solver.zig");

pub const WaterAndPoreVolumes = struct {
    soil_micropore_water_m3: f64,
    soil_micropore_air_m3: f64,
    root_aqueous_m3: f64,
    root_aqueous_non_band_m3: f64,
    root_aqueous_band_m3: f64,
    soil_matrix_volume_m3: f64,
};

pub const RootAllocation = struct {
    oxygen_fraction: f64,
    root_fraction: f64,
    ammonium_non_band_fraction: f64,
    ammonium_band_fraction: f64,
};

pub const RadialTransport = struct {
    tortuosity: f64,
    water_filled_threshold: f64,
    water_film_thickness_m: f64,
    root_radius_m: f64,
    root_surface_area_per_radius_m: f64,
    oxygen_diffusivity_m2_per_step: f64,
    methane_diffusivity_m2_per_step: f64,
    nitrous_oxide_diffusivity_m2_per_step: f64,
    ammonia_diffusivity_m2_per_step: f64,
    hydrogen_diffusivity_m2_per_step: f64,
};

pub const GasSolubility = struct {
    carbon_dioxide: f64,
    oxygen: f64,
    methane: f64,
    nitrous_oxide: f64,
    ammonia: f64,
    hydrogen: f64,
};

pub const SoilGasConcentration = struct {
    methane_g_m3: f64,
    nitrous_oxide_g_n_m3: f64,
    ammonia_g_n_m3: f64,
    hydrogen_g_h_m3: f64,
};

pub const AqueousGasMasses = struct {
    carbon_dioxide_g_c: f64,
    oxygen_g_o: f64,
    methane_g_c: f64,
    nitrous_oxide_g_n: f64,
    ammonia_non_band_g_n: f64,
    ammonia_band_g_n: f64,
    hydrogen_g_h: f64,
};

pub const AqueousMassPools = struct {
    soil: AqueousGasMasses,
    root: AqueousGasMasses,
};

pub const AqueousGasConcentrations = struct {
    carbon_dioxide_g_c_m3: f64,
    oxygen_g_o_m3: f64,
    methane_g_c_m3: f64,
    nitrous_oxide_g_n_m3: f64,
    ammonia_non_band_g_n_m3: f64,
    ammonia_band_g_n_m3: f64,
    hydrogen_g_h_m3: f64,
};

pub const AqueousConcentrationPair = struct {
    soil: AqueousGasConcentrations,
    root: AqueousGasConcentrations,
};

pub const ConvectiveMassFlow = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
    methane_g_c_per_step: f64,
    nitrous_oxide_g_n_per_step: f64,
    ammonia_non_band_g_n_per_step: f64,
    ammonia_band_g_n_per_step: f64,
    hydrogen_g_h_per_step: f64,
};

pub const PreparationInputs = struct {
    volumes: WaterAndPoreVolumes,
    allocation: RootAllocation,
    transport: RadialTransport,
    solubility: GasSolubility,
    gas_concentration: SoilGasConcentration,
};

pub const PreparedTransfer = struct {
    soil_oxygen_allocated_water_m3: f64,
    root_allocated_water_m3: f64,
    soil_allocated_air_m3: f64,
    combined_root_soil_water_m3: f64,
    non_band_water_m3: f64,
    band_water_m3: f64,
    soil_water_fraction: f64,
    path_length_log_ratio: f64,
    soil_to_root_oxygen_diffusivity_m3_per_step: f64,
    soil_to_root_methane_diffusivity_m3_per_step: f64,
    soil_to_root_nitrous_oxide_diffusivity_m3_per_step: f64,
    soil_to_root_ammonia_non_band_diffusivity_m3_per_step: f64,
    soil_to_root_ammonia_band_diffusivity_m3_per_step: f64,
    soil_to_root_hydrogen_diffusivity_m3_per_step: f64,
    soil_methane_g: f64,
    soil_nitrous_oxide_g_n: f64,
    soil_ammonia_g_n: f64,
    soil_hydrogen_g_h: f64,
    dissolved_capacity: GasSolubility,
    ammonia_band_dissolved_capacity_m3: f64,
    ammonia_non_band_air_m3: f64,
    ammonia_band_air_m3: f64,
};

pub const CellInputs = struct {
    preparation: PreparationInputs,
    soil_aqueous_oxygen_g_o: f64,
    preceding_soil_oxygen_flux_g_o_per_step: f64,
    root_aqueous_oxygen_g_o: f64,
    equilibrium_gas_oxygen_g_o_m3: f64,
    internal_root_diffusivity_m3_per_step: f64,
    root_water_uptake_m3_per_step: f64,
    oxygen_demand_per_plant_g_o_per_step: f64,
    oxygen_half_saturation_g_o_m3: f64,
    negligible_oxygen_g_o: f64,
    oxygen_competition_fraction: f64,
    plant_population: f64,
    other_aqueous_gas_masses: AqueousMassPools,
};

pub const CellResult = struct {
    prepared: PreparedTransfer,
    oxygen: oxygen.Result,
    aqueous_concentrations: AqueousConcentrationPair,
    convective_mass_flow: ConvectiveMassFlow,
    oxygen_balance_residual_g_o_per_plant_step: f64,
};

/// UPTAKE.F 2118--2154, preserving the legacy statement order.
pub fn prepare(inputs: PreparationInputs) !PreparedTransfer {
    try validatePreparation(inputs);
    const v = inputs.volumes;
    const a = inputs.allocation;
    const soil_oxygen_water = v.soil_micropore_water_m3 * a.oxygen_fraction;
    const root_water = v.soil_micropore_water_m3 * a.root_fraction;
    const root_air = v.soil_micropore_air_m3 * a.root_fraction;
    const combined_water = v.root_aqueous_m3 + root_water;
    const root_non_band_water = root_water * a.ammonium_non_band_fraction;
    const root_band_water = root_water * a.ammonium_band_fraction;
    const non_band_water = v.root_aqueous_non_band_m3 + root_non_band_water;
    const band_water = v.root_aqueous_band_m3 + root_band_water;
    const water_fraction = if (v.soil_matrix_volume_m3 > 0)
        @max(0, v.soil_micropore_water_m3 / v.soil_matrix_volume_m3)
    else
        0;

    var result = PreparedTransfer{
        .soil_oxygen_allocated_water_m3 = soil_oxygen_water,
        .root_allocated_water_m3 = root_water,
        .soil_allocated_air_m3 = root_air,
        .combined_root_soil_water_m3 = combined_water,
        .non_band_water_m3 = non_band_water,
        .band_water_m3 = band_water,
        .soil_water_fraction = water_fraction,
        .path_length_log_ratio = 0,
        .soil_to_root_oxygen_diffusivity_m3_per_step = 0,
        .soil_to_root_methane_diffusivity_m3_per_step = 0,
        .soil_to_root_nitrous_oxide_diffusivity_m3_per_step = 0,
        .soil_to_root_ammonia_non_band_diffusivity_m3_per_step = 0,
        .soil_to_root_ammonia_band_diffusivity_m3_per_step = 0,
        .soil_to_root_hydrogen_diffusivity_m3_per_step = 0,
        .soil_methane_g = 0,
        .soil_nitrous_oxide_g_n = 0,
        .soil_ammonia_g_n = 0,
        .soil_hydrogen_g_h = 0,
        .dissolved_capacity = .{
            .carbon_dioxide = 0,
            .oxygen = 0,
            .methane = 0,
            .nitrous_oxide = 0,
            .ammonia = 0,
            .hydrogen = 0,
        },
        .ammonia_band_dissolved_capacity_m3 = 0,
        .ammonia_non_band_air_m3 = 0,
        .ammonia_band_air_m3 = 0,
    };
    if (water_fraction <= inputs.transport.water_filled_threshold or
        a.root_fraction <= 0) return result;

    const transport = inputs.transport;
    const theta_transport = transport.tortuosity * water_fraction;
    const path = @log(
        (transport.water_film_thickness_m + transport.root_radius_m) /
            transport.root_radius_m,
    );
    if (!std.math.isFinite(path) or path <= 0)
        return error.InvalidRootGasDiffusionPath;
    const area_path_m = transport.root_surface_area_per_radius_m / path;
    result.path_length_log_ratio = path;
    result.soil_to_root_oxygen_diffusivity_m3_per_step =
        theta_transport * transport.oxygen_diffusivity_m2_per_step * area_path_m;
    result.soil_to_root_methane_diffusivity_m3_per_step =
        theta_transport * transport.methane_diffusivity_m2_per_step * area_path_m;
    result.soil_to_root_nitrous_oxide_diffusivity_m3_per_step =
        theta_transport * transport.nitrous_oxide_diffusivity_m2_per_step * area_path_m;
    result.soil_to_root_ammonia_non_band_diffusivity_m3_per_step =
        theta_transport * transport.ammonia_diffusivity_m2_per_step * area_path_m *
        a.ammonium_non_band_fraction;
    result.soil_to_root_ammonia_band_diffusivity_m3_per_step =
        theta_transport * transport.ammonia_diffusivity_m2_per_step * area_path_m *
        a.ammonium_band_fraction;
    result.soil_to_root_hydrogen_diffusivity_m3_per_step =
        theta_transport * transport.hydrogen_diffusivity_m2_per_step * area_path_m;
    const gas = inputs.gas_concentration;
    result.soil_methane_g = gas.methane_g_m3 * root_air;
    result.soil_nitrous_oxide_g_n = gas.nitrous_oxide_g_n_m3 * root_air;
    result.soil_ammonia_g_n = gas.ammonia_g_n_m3 * root_air;
    result.soil_hydrogen_g_h = gas.hydrogen_g_h_m3 * root_air;
    const s = inputs.solubility;
    result.dissolved_capacity = .{
        .carbon_dioxide = root_water * s.carbon_dioxide,
        .oxygen = root_water * s.oxygen,
        .methane = root_water * s.methane,
        .nitrous_oxide = root_water * s.nitrous_oxide,
        .ammonia = root_water * s.ammonia * a.ammonium_non_band_fraction,
        .hydrogen = root_water * s.hydrogen,
    };
    result.ammonia_band_dissolved_capacity_m3 =
        root_water * s.ammonia * a.ammonium_band_fraction;
    result.ammonia_non_band_air_m3 = root_air * a.ammonium_non_band_fraction;
    result.ammonia_band_air_m3 = root_air * a.ammonium_band_fraction;
    return result;
}

/// Computes every runtime cell into caller-owned scratch before committing any
/// result, so a failed cell cannot partially mutate the destination.
pub fn solveCells(
    cells: []const CellInputs,
    policy: oxygen.SolverPolicy,
    scratch: []CellResult,
    destination: []CellResult,
) !void {
    if (scratch.len != cells.len or destination.len != cells.len)
        return error.RootGasTransferDimensionMismatch;
    for (cells, scratch) |cell, *candidate| {
        try validateAqueousMasses(cell.other_aqueous_gas_masses);
        const prepared = try prepare(cell.preparation);
        const concentrations = try concentrationsInSourceOrder(cell, prepared);
        const mass_flow = ConvectiveMassFlow{
            .carbon_dioxide_g_c_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.carbon_dioxide_g_c_m3,
            .oxygen_g_o_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.oxygen_g_o_m3,
            .methane_g_c_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.methane_g_c_m3,
            .nitrous_oxide_g_n_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.nitrous_oxide_g_n_m3,
            .ammonia_non_band_g_n_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.ammonia_non_band_g_n_m3 *
                cell.preparation.allocation.ammonium_non_band_fraction,
            .ammonia_band_g_n_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.ammonia_band_g_n_m3 *
                cell.preparation.allocation.ammonium_band_fraction,
            .hydrogen_g_h_per_step = cell.root_water_uptake_m3_per_step *
                concentrations.soil.hydrogen_g_h_m3,
        };
        const solved = try oxygen.solve(.{
            .soil_aqueous_oxygen_g_o = cell.soil_aqueous_oxygen_g_o,
            .preceding_soil_oxygen_aqueous_flux_g_o_per_step = cell.preceding_soil_oxygen_flux_g_o_per_step,
            .soil_oxygen_allocated_water_volume_m3 = prepared.soil_oxygen_allocated_water_m3,
            .root_aqueous_oxygen_g_o = cell.root_aqueous_oxygen_g_o,
            .root_aqueous_volume_m3 = cell.preparation.volumes.root_aqueous_m3,
            .equilibrium_oxygen_concentration_g_o_per_m3 = cell.equilibrium_gas_oxygen_g_o_m3 * cell.preparation.solubility.oxygen,
            .soil_to_root_diffusivity_m3_per_step = prepared.soil_to_root_oxygen_diffusivity_m3_per_step,
            .internal_root_diffusivity_m3_per_step = cell.internal_root_diffusivity_m3_per_step,
            .root_water_uptake_m3_per_step = cell.root_water_uptake_m3_per_step,
            .oxygen_demand_per_plant_g_o_per_step = cell.oxygen_demand_per_plant_g_o_per_step,
            .oxygen_half_saturation_g_o_per_m3 = cell.oxygen_half_saturation_g_o_m3,
            .negligible_oxygen_g_o = cell.negligible_oxygen_g_o,
            .oxygen_competition_fraction = cell.oxygen_competition_fraction,
            .plant_population = cell.plant_population,
        }, policy);
        const residual = solved.soil_to_root_oxygen_flux_g_o_per_step +
            solved.internal_root_oxygen_flux_g_o_per_step -
            solved.oxygen_uptake_per_plant_g_o_per_step;
        if (!std.math.isFinite(residual))
            return error.NonFiniteRootGasTransferBalance;
        candidate.* = .{
            .prepared = prepared,
            .oxygen = solved,
            .aqueous_concentrations = concentrations,
            .convective_mass_flow = mass_flow,
            .oxygen_balance_residual_g_o_per_plant_step = residual,
        };
    }
    @memcpy(destination, scratch);
}

/// UPTAKE.F 2174--2196 concentration and convective-flow statements.
fn concentrationsInSourceOrder(
    cell: CellInputs,
    prepared: PreparedTransfer,
) !AqueousConcentrationPair {
    const soil_volume = cell.preparation.volumes.soil_micropore_water_m3 *
        cell.preparation.allocation.root_fraction;
    const root_volume = cell.preparation.volumes.root_aqueous_m3;
    if (soil_volume <= 0 or root_volume <= 0)
        return error.InvalidRootGasAqueousVolume;
    const masses = cell.other_aqueous_gas_masses;
    const equilibrium_oxygen = cell.equilibrium_gas_oxygen_g_o_m3 *
        cell.preparation.solubility.oxygen;
    const soil_oxygen_amount = cell.soil_aqueous_oxygen_g_o +
        cell.preceding_soil_oxygen_flux_g_o_per_step;
    const result = AqueousConcentrationPair{
        .soil = .{
            .carbon_dioxide_g_c_m3 = @max(0, masses.soil.carbon_dioxide_g_c / soil_volume),
            .oxygen_g_o_m3 = @min(
                equilibrium_oxygen,
                @max(0, soil_oxygen_amount / prepared.soil_oxygen_allocated_water_m3),
            ),
            .methane_g_c_m3 = @max(0, masses.soil.methane_g_c / soil_volume),
            .nitrous_oxide_g_n_m3 = @max(0, masses.soil.nitrous_oxide_g_n / soil_volume),
            .ammonia_non_band_g_n_m3 = @max(0, masses.soil.ammonia_non_band_g_n / soil_volume),
            .ammonia_band_g_n_m3 = @max(0, masses.soil.ammonia_band_g_n / soil_volume),
            .hydrogen_g_h_m3 = @max(0, masses.soil.hydrogen_g_h / soil_volume),
        },
        .root = .{
            .carbon_dioxide_g_c_m3 = @max(0, masses.root.carbon_dioxide_g_c / root_volume),
            .oxygen_g_o_m3 = @min(
                equilibrium_oxygen,
                @max(0, cell.root_aqueous_oxygen_g_o / root_volume),
            ),
            .methane_g_c_m3 = @max(0, masses.root.methane_g_c / root_volume),
            .nitrous_oxide_g_n_m3 = @max(0, masses.root.nitrous_oxide_g_n / root_volume),
            .ammonia_non_band_g_n_m3 = @max(0, masses.root.ammonia_non_band_g_n / root_volume),
            .ammonia_band_g_n_m3 = 0,
            .hydrogen_g_h_m3 = @max(0, masses.root.hydrogen_g_h / root_volume),
        },
    };
    inline for (@typeInfo(AqueousGasConcentrations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.soil, field.name)) or
            !std.math.isFinite(@field(result.root, field.name)))
            return error.NonFiniteRootGasConcentration;
    }
    return result;
}

fn validatePreparation(inputs: PreparationInputs) !void {
    inline for (@typeInfo(WaterAndPoreVolumes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.volumes, field.name)) or
            @field(inputs.volumes, field.name) < 0)
            return error.InvalidRootGasPreparationInput;
    inline for (@typeInfo(RootAllocation).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.allocation, field.name)) or
            @field(inputs.allocation, field.name) < 0 or
            @field(inputs.allocation, field.name) > 1)
            return error.InvalidRootGasPreparationInput;
    inline for (@typeInfo(RadialTransport).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.transport, field.name)) or
            @field(inputs.transport, field.name) < 0)
            return error.InvalidRootGasPreparationInput;
    inline for (@typeInfo(GasSolubility).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.solubility, field.name)) or
            @field(inputs.solubility, field.name) < 0)
            return error.InvalidRootGasPreparationInput;
    inline for (@typeInfo(SoilGasConcentration).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.gas_concentration, field.name)) or
            @field(inputs.gas_concentration, field.name) < 0)
            return error.InvalidRootGasPreparationInput;
    if (inputs.transport.root_radius_m <= 0)
        return error.InvalidRootGasPreparationInput;
}

fn validateAqueousMasses(masses: AqueousMassPools) !void {
    inline for (@typeInfo(AqueousGasMasses).@"struct".fields) |field| {
        const soil_value = @field(masses.soil, field.name);
        const root_value = @field(masses.root, field.name);
        if (!std.math.isFinite(soil_value) or soil_value < 0 or
            !std.math.isFinite(root_value) or root_value < 0)
            return error.InvalidRootGasAqueousMass;
    }
}

fn testCell() CellInputs {
    return .{
        .preparation = .{
            .volumes = .{
                .soil_micropore_water_m3 = 2,
                .soil_micropore_air_m3 = 1,
                .root_aqueous_m3 = 1,
                .root_aqueous_non_band_m3 = 0.2,
                .root_aqueous_band_m3 = 0.1,
                .soil_matrix_volume_m3 = 4,
            },
            .allocation = .{
                .oxygen_fraction = 1,
                .root_fraction = 0.5,
                .ammonium_non_band_fraction = 0.6,
                .ammonium_band_fraction = 0.4,
            },
            .transport = .{
                .tortuosity = 0.8,
                .water_filled_threshold = 0.1,
                .water_film_thickness_m = 0.001,
                .root_radius_m = 0.001,
                .root_surface_area_per_radius_m = 1,
                .oxygen_diffusivity_m2_per_step = 0.5,
                .methane_diffusivity_m2_per_step = 0.4,
                .nitrous_oxide_diffusivity_m2_per_step = 0.3,
                .ammonia_diffusivity_m2_per_step = 0.2,
                .hydrogen_diffusivity_m2_per_step = 0.1,
            },
            .solubility = .{
                .carbon_dioxide = 1,
                .oxygen = 1,
                .methane = 1,
                .nitrous_oxide = 1,
                .ammonia = 1,
                .hydrogen = 1,
            },
            .gas_concentration = .{
                .methane_g_m3 = 1,
                .nitrous_oxide_g_n_m3 = 1,
                .ammonia_g_n_m3 = 1,
                .hydrogen_g_h_m3 = 1,
            },
        },
        .soil_aqueous_oxygen_g_o = 4,
        .preceding_soil_oxygen_flux_g_o_per_step = 0,
        .root_aqueous_oxygen_g_o = 1,
        .equilibrium_gas_oxygen_g_o_m3 = 10,
        .internal_root_diffusivity_m3_per_step = 0.25,
        .root_water_uptake_m3_per_step = 0.1,
        .oxygen_demand_per_plant_g_o_per_step = 0.8,
        .oxygen_half_saturation_g_o_m3 = 0.2,
        .negligible_oxygen_g_o = 1e-12,
        .oxygen_competition_fraction = 0.4,
        .plant_population = 2,
        .other_aqueous_gas_masses = .{
            .soil = .{
                .carbon_dioxide_g_c = 2,
                .oxygen_g_o = 0,
                .methane_g_c = 1,
                .nitrous_oxide_g_n = 0.5,
                .ammonia_non_band_g_n = 0.4,
                .ammonia_band_g_n = 0.3,
                .hydrogen_g_h = 0.2,
            },
            .root = .{
                .carbon_dioxide_g_c = 1,
                .oxygen_g_o = 0,
                .methane_g_c = 0.5,
                .nitrous_oxide_g_n = 0.25,
                .ammonia_non_band_g_n = 0.2,
                .ammonia_band_g_n = 0,
                .hydrogen_g_h = 0.1,
            },
        },
    };
}

fn testPolicy(maximum_iterations: u16) oxygen.SolverPolicy {
    return .{
        .maximum_iterations = maximum_iterations,
        .absolute_tolerance_g_o_per_step = 1e-12,
        .relative_tolerance = 1e-10,
        .picard_relaxation = 0.5,
    };
}

test "runtime cell transaction closes root oxygen flux balance" {
    const cells = [_]CellInputs{testCell()};
    var scratch: [1]CellResult = undefined;
    var destination: [1]CellResult = undefined;
    try solveCells(&cells, testPolicy(30), &scratch, &destination);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        destination[0].oxygen_balance_residual_g_o_per_plant_step,
        1e-10,
    );
    try std.testing.expect(destination[0].oxygen.iterations < 30);
    try std.testing.expectApproxEqAbs(
        @as(f64, 2),
        destination[0].aqueous_concentrations.soil.carbon_dioxide_g_c_m3,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.2),
        destination[0].convective_mass_flow.carbon_dioxide_g_c_per_step,
        1e-12,
    );
}

test "batch failure leaves every destination cell unchanged" {
    var cells = [_]CellInputs{ testCell(), testCell() };
    cells[1].oxygen_demand_per_plant_g_o_per_step = std.math.nan(f64);
    var scratch: [2]CellResult = undefined;
    var destination = [_]CellResult{std.mem.zeroes(CellResult)} ** 2;
    destination[0].oxygen_balance_residual_g_o_per_plant_step = 41;
    destination[1].oxygen_balance_residual_g_o_per_plant_step = 42;
    try std.testing.expectError(
        error.InvalidRootOxygenUptakeInput,
        solveCells(&cells, testPolicy(30), &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].oxygen_balance_residual_g_o_per_plant_step);
    try std.testing.expectEqual(@as(f64, 42), destination[1].oxygen_balance_residual_g_o_per_plant_step);
}

test "runtime dimensions must match caller-owned storage" {
    const cells = [_]CellInputs{testCell()};
    var scratch: [1]CellResult = undefined;
    var destination: [0]CellResult = .{};
    try std.testing.expectError(
        error.RootGasTransferDimensionMismatch,
        solveCells(&cells, testPolicy(30), &scratch, &destination),
    );
}
