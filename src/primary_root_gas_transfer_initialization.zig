const std = @import("std");

pub const RootDomain = enum {
    primary_root,
    mycorrhizal,
};

pub const GasSpeciesValues = struct {
    carbon_dioxide: f64,
    oxygen: f64,
    methane: f64,
    nitrous_oxide: f64,
    ammonia: f64,
    hydrogen: f64,
};

pub const Inputs = struct {
    root_domain: RootDomain,
    emerged: bool,
    root_length_per_plant_m: f64,
    negligible_root_length_m: f64,
    soil_water_film_thickness_m: f64,
    root_radius_m: f64,
    root_surface_area_per_radius_m: f64,
    unallocated_oxygen_aqueous_diffusivity_m2_per_step: f64,
    root_aqueous_volume_m3: f64,
    gas_solubility: GasSpeciesValues,
    gaseous_diffusivity_m2_per_step: GasSpeciesValues,
    effective_cross_section_per_length_m: f64,
    root_porosity_m3_per_m3: f64,
    root_carbon_dioxide_flux_g_c_per_h: f64,
    gas_flux_timestep_h_per_step: f64,
};

pub const Result = struct {
    internal_transfer_active: bool,
    radial_path_log_factor: f64,
    effective_surface_area_per_path_m: f64,
    internal_oxygen_diffusivity_m3_per_step: f64,
    aqueous_capacity_m3: GasSpeciesValues,
    atmosphere_conductance_m3_per_step: GasSpeciesValues,
    gaseous_aqueous_equilibration_fraction: f64,
    root_carbon_dioxide_source_g_c_per_step: f64,
};

/// UPTAKE.F 2054--2087. Initializes primary-root internal and atmospheric gas
/// transfer, then applies the post-branch porosity and CO2-source terms.
pub fn calculate(inputs: Inputs) !Result {
    try validateCommon(inputs);
    var path_log: f64 = 0;
    var surface_per_path: f64 = 0;
    var oxygen_diffusivity: f64 = 0;
    var capacity = zeroSpecies();
    var atmosphere = zeroSpecies();
    const active =
        inputs.root_domain == .primary_root and inputs.emerged and
        inputs.root_length_per_plant_m > inputs.negligible_root_length_m;
    if (active) {
        if (inputs.root_radius_m <= 0 or
            inputs.soil_water_film_thickness_m <= 0)
            return error.InvalidPrimaryRootGasTransferGeometry;
        path_log = @log(
            (inputs.soil_water_film_thickness_m + inputs.root_radius_m) /
                inputs.root_radius_m,
        );
        if (path_log == 0)
            return error.SingularPrimaryRootGasTransferPath;
        surface_per_path =
            inputs.root_surface_area_per_radius_m / path_log;
        oxygen_diffusivity =
            inputs.unallocated_oxygen_aqueous_diffusivity_m2_per_step *
            surface_per_path;
        inline for (@typeInfo(GasSpeciesValues).@"struct".fields) |field| {
            @field(capacity, field.name) =
                inputs.root_aqueous_volume_m3 *
                @field(inputs.gas_solubility, field.name);
            @field(atmosphere, field.name) =
                @field(inputs.gaseous_diffusivity_m2_per_step, field.name) *
                inputs.effective_cross_section_per_length_m;
        }
    }
    const equilibration = @min(1, inputs.root_porosity_m3_per_m3);
    const carbon_dioxide_source =
        -inputs.root_carbon_dioxide_flux_g_c_per_h *
        inputs.gas_flux_timestep_h_per_step;
    const result = Result{
        .internal_transfer_active = active,
        .radial_path_log_factor = path_log,
        .effective_surface_area_per_path_m = surface_per_path,
        .internal_oxygen_diffusivity_m3_per_step = oxygen_diffusivity,
        .aqueous_capacity_m3 = capacity,
        .atmosphere_conductance_m3_per_step = atmosphere,
        .gaseous_aqueous_equilibration_fraction = equilibration,
        .root_carbon_dioxide_source_g_c_per_step = carbon_dioxide_source,
    };
    try validateResult(result);
    return result;
}

fn zeroSpecies() GasSpeciesValues {
    return .{
        .carbon_dioxide = 0,
        .oxygen = 0,
        .methane = 0,
        .nitrous_oxide = 0,
        .ammonia = 0,
        .hydrogen = 0,
    };
}

fn validateCommon(inputs: Inputs) !void {
    inline for (.{
        inputs.root_length_per_plant_m,
        inputs.negligible_root_length_m,
        inputs.soil_water_film_thickness_m,
        inputs.root_radius_m,
        inputs.root_surface_area_per_radius_m,
        inputs.unallocated_oxygen_aqueous_diffusivity_m2_per_step,
        inputs.root_aqueous_volume_m3,
        inputs.effective_cross_section_per_length_m,
        inputs.root_porosity_m3_per_m3,
        inputs.root_carbon_dioxide_flux_g_c_per_h,
        inputs.gas_flux_timestep_h_per_step,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidPrimaryRootGasTransferInput;
    if (inputs.root_length_per_plant_m < 0 or
        inputs.negligible_root_length_m < 0 or
        inputs.root_surface_area_per_radius_m < 0 or
        inputs.unallocated_oxygen_aqueous_diffusivity_m2_per_step < 0 or
        inputs.root_aqueous_volume_m3 < 0 or
        inputs.effective_cross_section_per_length_m < 0 or
        inputs.root_porosity_m3_per_m3 < 0 or
        inputs.gas_flux_timestep_h_per_step < 0)
        return error.InvalidPrimaryRootGasTransferInput;
    inline for (@typeInfo(GasSpeciesValues).@"struct".fields) |field| {
        const solubility = @field(inputs.gas_solubility, field.name);
        const diffusivity =
            @field(inputs.gaseous_diffusivity_m2_per_step, field.name);
        if (!std.math.isFinite(solubility) or solubility < 0 or
            !std.math.isFinite(diffusivity) or diffusivity < 0)
            return error.InvalidPrimaryRootGasTransferInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (.{
        result.radial_path_log_factor,
        result.effective_surface_area_per_path_m,
        result.internal_oxygen_diffusivity_m3_per_step,
        result.gaseous_aqueous_equilibration_fraction,
        result.root_carbon_dioxide_source_g_c_per_step,
    }) |value|
        if (!std.math.isFinite(value))
            return error.NonFinitePrimaryRootGasTransferResult;
    inline for (@typeInfo(GasSpeciesValues).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.aqueous_capacity_m3, field.name)) or
            !std.math.isFinite(
                @field(result.atmosphere_conductance_m3_per_step, field.name),
            ))
            return error.NonFinitePrimaryRootGasTransferResult;
}

fn species(base: f64) GasSpeciesValues {
    return .{
        .carbon_dioxide = base + 1,
        .oxygen = base + 2,
        .methane = base + 3,
        .nitrous_oxide = base + 4,
        .ammonia = base + 5,
        .hydrogen = base + 6,
    };
}

fn sourceInputs() Inputs {
    return .{
        .root_domain = .primary_root,
        .emerged = true,
        .root_length_per_plant_m = 2,
        .negligible_root_length_m = 1e-12,
        .soil_water_film_thickness_m = 0.01,
        .root_radius_m = 0.005,
        .root_surface_area_per_radius_m = 3,
        .unallocated_oxygen_aqueous_diffusivity_m2_per_step = 0.2,
        .root_aqueous_volume_m3 = 4,
        .gas_solubility = species(0),
        .gaseous_diffusivity_m2_per_step = species(10),
        .effective_cross_section_per_length_m = 0.5,
        .root_porosity_m3_per_m3 = 0.8,
        .root_carbon_dioxide_flux_g_c_per_h = 2,
        .gas_flux_timestep_h_per_step = 0.25,
    };
}

test "UPTAKE active primary root initializes path capacity and conductance" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const path = @log((0.01 + 0.005) / 0.005);
    try std.testing.expect(result.internal_transfer_active);
    try std.testing.expectEqual(path, result.radial_path_log_factor);
    try std.testing.expectEqual(3.0 / path, result.effective_surface_area_per_path_m);
    try std.testing.expectEqual(0.2 * (3.0 / path), result.internal_oxygen_diffusivity_m3_per_step);
    try std.testing.expectEqual(@as(f64, 4), result.aqueous_capacity_m3.carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 5.5), result.atmosphere_conductance_m3_per_step.carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 0.8), result.gaseous_aqueous_equilibration_fraction);
    try std.testing.expectEqual(@as(f64, -0.5), result.root_carbon_dioxide_source_g_c_per_step);
}

test "inactive mycorrhizal branch zeros transfer but retains post-branch terms" {
    var inputs = sourceInputs();
    inputs.root_domain = .mycorrhizal;
    inputs.root_radius_m = 0;
    inputs.soil_water_film_thickness_m = 0;
    inputs.root_porosity_m3_per_m3 = 2;
    const result = try calculate(inputs);
    try std.testing.expect(!result.internal_transfer_active);
    try std.testing.expectEqual(@as(f64, 0), result.internal_oxygen_diffusivity_m3_per_step);
    try std.testing.expectEqual(@as(f64, 1), result.gaseous_aqueous_equilibration_fraction);
    try std.testing.expectEqual(@as(f64, -0.5), result.root_carbon_dioxide_source_g_c_per_step);
}

test "active zero radial geometry fails explicitly" {
    var inputs = sourceInputs();
    inputs.soil_water_film_thickness_m = 0;
    try std.testing.expectError(
        error.InvalidPrimaryRootGasTransferGeometry,
        calculate(inputs),
    );
}
