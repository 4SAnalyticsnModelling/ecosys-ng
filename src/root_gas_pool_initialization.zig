const std = @import("std");

pub const FourPhaseAmounts = struct {
    root_gaseous: f64,
    root_aqueous: f64,
    soil_gaseous: f64,
    soil_aqueous: f64,
};

pub const ThreePhaseAmounts = struct {
    root_gaseous: f64,
    root_aqueous: f64,
    soil_aqueous: f64,
};

pub const AmmoniaAmounts = struct {
    root_gaseous: f64,
    root_aqueous: f64,
    soil_non_band_aqueous: f64,
    soil_band_aqueous: f64,
};

pub const Inputs = struct {
    carbon_dioxide: FourPhaseAmounts,
    oxygen: FourPhaseAmounts,
    methane: ThreePhaseAmounts,
    nitrous_oxide: ThreePhaseAmounts,
    ammonia: AmmoniaAmounts,
    hydrogen: ThreePhaseAmounts,
    root_mass_fraction: f64,
    oxygen_competition_fraction: f64,
    root_aqueous_volume_m3: f64,
    soil_non_band_fraction: f64,
    soil_band_fraction: f64,
    oxygen_demand_g_o_per_step: f64,
    gas_flux_timestep_h_per_step: f64,
    plant_population: f64,
    preceding_soil_oxygen_gas_flux_g_o_per_step: f64,
    preceding_soil_carbon_dioxide_gas_flux_g_c_per_step: f64,
    preceding_soil_oxygen_aqueous_flux_g_o_per_step: f64,
};

pub const Result = struct {
    carbon_dioxide: FourPhaseAmounts,
    oxygen: FourPhaseAmounts,
    methane: ThreePhaseAmounts,
    nitrous_oxide: ThreePhaseAmounts,
    ammonia: AmmoniaAmounts,
    hydrogen: ThreePhaseAmounts,
    root_non_band_aqueous_volume_m3: f64,
    root_band_aqueous_volume_m3: f64,
    oxygen_demand_per_plant_g_o_per_step: f64,
    allocated_soil_oxygen_gas_flux_g_o_per_step: f64,
    allocated_soil_carbon_dioxide_gas_flux_g_c_per_step: f64,
    allocated_soil_oxygen_aqueous_flux_g_o_per_step: f64,
};

/// UPTAKE.F 1904--1930. Initializes root pools and allocates soil pools and
/// preceding fluxes to one runtime root/mycorrhizal domain.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const result = Result{
        .carbon_dioxide = .{
            .root_gaseous = @max(0, inputs.carbon_dioxide.root_gaseous),
            .root_aqueous = @max(0, inputs.carbon_dioxide.root_aqueous),
            .soil_gaseous = @max(
                0,
                inputs.carbon_dioxide.soil_gaseous *
                    inputs.root_mass_fraction,
            ),
            .soil_aqueous = @max(
                0,
                inputs.carbon_dioxide.soil_aqueous *
                    inputs.root_mass_fraction,
            ),
        },
        .oxygen = .{
            .root_gaseous = @max(0, inputs.oxygen.root_gaseous),
            .root_aqueous = @max(0, inputs.oxygen.root_aqueous),
            .soil_gaseous = @max(
                0,
                inputs.oxygen.soil_gaseous *
                    inputs.oxygen_competition_fraction,
            ),
            .soil_aqueous = @max(
                0,
                inputs.oxygen.soil_aqueous *
                    inputs.oxygen_competition_fraction,
            ),
        },
        .methane = .{
            .root_gaseous = inputs.methane.root_gaseous,
            .root_aqueous = inputs.methane.root_aqueous,
            .soil_aqueous = inputs.methane.soil_aqueous * inputs.root_mass_fraction,
        },
        .nitrous_oxide = .{
            .root_gaseous = inputs.nitrous_oxide.root_gaseous,
            .root_aqueous = inputs.nitrous_oxide.root_aqueous,
            .soil_aqueous = inputs.nitrous_oxide.soil_aqueous * inputs.root_mass_fraction,
        },
        .ammonia = .{
            .root_gaseous = inputs.ammonia.root_gaseous,
            .root_aqueous = inputs.ammonia.root_aqueous,
            .soil_non_band_aqueous = inputs.ammonia.soil_non_band_aqueous *
                inputs.root_mass_fraction,
            .soil_band_aqueous = inputs.ammonia.soil_band_aqueous *
                inputs.root_mass_fraction,
        },
        .hydrogen = .{
            .root_gaseous = inputs.hydrogen.root_gaseous,
            .root_aqueous = inputs.hydrogen.root_aqueous,
            .soil_aqueous = inputs.hydrogen.soil_aqueous * inputs.root_mass_fraction,
        },
        .root_non_band_aqueous_volume_m3 = inputs.root_aqueous_volume_m3 *
            inputs.soil_non_band_fraction,
        .root_band_aqueous_volume_m3 = inputs.root_aqueous_volume_m3 * inputs.soil_band_fraction,
        .oxygen_demand_per_plant_g_o_per_step = inputs.oxygen_demand_g_o_per_step *
            inputs.gas_flux_timestep_h_per_step /
            inputs.plant_population,
        .allocated_soil_oxygen_gas_flux_g_o_per_step = inputs.preceding_soil_oxygen_gas_flux_g_o_per_step *
            inputs.oxygen_competition_fraction *
            inputs.gas_flux_timestep_h_per_step,
        .allocated_soil_carbon_dioxide_gas_flux_g_c_per_step = inputs.preceding_soil_carbon_dioxide_gas_flux_g_c_per_step *
            inputs.oxygen_competition_fraction *
            inputs.gas_flux_timestep_h_per_step,
        .allocated_soil_oxygen_aqueous_flux_g_o_per_step = inputs.preceding_soil_oxygen_aqueous_flux_g_o_per_step *
            inputs.oxygen_competition_fraction *
            inputs.gas_flux_timestep_h_per_step,
    };
    try validateResult(result);
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(FourPhaseAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.carbon_dioxide, field.name)) or
            !std.math.isFinite(@field(inputs.oxygen, field.name)))
            return error.InvalidRootGasPoolInitializationInput;
    }
    inline for (@typeInfo(ThreePhaseAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.methane, field.name)) or
            !std.math.isFinite(@field(inputs.nitrous_oxide, field.name)) or
            !std.math.isFinite(@field(inputs.hydrogen, field.name)))
            return error.InvalidRootGasPoolInitializationInput;
    }
    inline for (@typeInfo(AmmoniaAmounts).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.ammonia, field.name)))
            return error.InvalidRootGasPoolInitializationInput;
    inline for (.{
        inputs.root_mass_fraction,
        inputs.oxygen_competition_fraction,
        inputs.root_aqueous_volume_m3,
        inputs.soil_non_band_fraction,
        inputs.soil_band_fraction,
        inputs.oxygen_demand_g_o_per_step,
        inputs.gas_flux_timestep_h_per_step,
        inputs.plant_population,
        inputs.preceding_soil_oxygen_gas_flux_g_o_per_step,
        inputs.preceding_soil_carbon_dioxide_gas_flux_g_c_per_step,
        inputs.preceding_soil_oxygen_aqueous_flux_g_o_per_step,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidRootGasPoolInitializationInput;
    if (inputs.root_mass_fraction < 0 or
        inputs.oxygen_competition_fraction < 0 or
        inputs.root_aqueous_volume_m3 < 0 or
        inputs.soil_non_band_fraction < 0 or
        inputs.soil_band_fraction < 0 or
        inputs.oxygen_demand_g_o_per_step < 0 or
        inputs.gas_flux_timestep_h_per_step < 0 or
        inputs.plant_population <= 0)
        return error.InvalidRootGasPoolInitializationInput;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(FourPhaseAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.carbon_dioxide, field.name)) or
            !std.math.isFinite(@field(result.oxygen, field.name)))
            return error.NonFiniteRootGasPoolInitializationResult;
    }
    inline for (@typeInfo(ThreePhaseAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.methane, field.name)) or
            !std.math.isFinite(@field(result.nitrous_oxide, field.name)) or
            !std.math.isFinite(@field(result.hydrogen, field.name)))
            return error.NonFiniteRootGasPoolInitializationResult;
    }
    inline for (@typeInfo(AmmoniaAmounts).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.ammonia, field.name)))
            return error.NonFiniteRootGasPoolInitializationResult;
    inline for (.{
        result.root_non_band_aqueous_volume_m3,
        result.root_band_aqueous_volume_m3,
        result.oxygen_demand_per_plant_g_o_per_step,
        result.allocated_soil_oxygen_gas_flux_g_o_per_step,
        result.allocated_soil_carbon_dioxide_gas_flux_g_c_per_step,
        result.allocated_soil_oxygen_aqueous_flux_g_o_per_step,
    }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteRootGasPoolInitializationResult;
}

fn sourceInputs() Inputs {
    return .{
        .carbon_dioxide = .{
            .root_gaseous = -1,
            .root_aqueous = 2,
            .soil_gaseous = -3,
            .soil_aqueous = 4,
        },
        .oxygen = .{
            .root_gaseous = -5,
            .root_aqueous = 6,
            .soil_gaseous = -7,
            .soil_aqueous = 8,
        },
        .methane = .{ .root_gaseous = -1, .root_aqueous = -2, .soil_aqueous = -3 },
        .nitrous_oxide = .{ .root_gaseous = -4, .root_aqueous = -5, .soil_aqueous = -6 },
        .ammonia = .{
            .root_gaseous = -7,
            .root_aqueous = -8,
            .soil_non_band_aqueous = -9,
            .soil_band_aqueous = -10,
        },
        .hydrogen = .{ .root_gaseous = -11, .root_aqueous = -12, .soil_aqueous = -13 },
        .root_mass_fraction = 0.25,
        .oxygen_competition_fraction = 0.5,
        .root_aqueous_volume_m3 = 4,
        .soil_non_band_fraction = 0.6,
        .soil_band_fraction = 0.4,
        .oxygen_demand_g_o_per_step = 8,
        .gas_flux_timestep_h_per_step = 0.5,
        .plant_population = 2,
        .preceding_soil_oxygen_gas_flux_g_o_per_step = -2,
        .preceding_soil_carbon_dioxide_gas_flux_g_c_per_step = 4,
        .preceding_soil_oxygen_aqueous_flux_g_o_per_step = -6,
    };
}

test "UPTAKE gas pool initialization preserves selective floors and allocation" {
    const result = try calculate(sourceInputs());
    try std.testing.expectEqual(@as(f64, 0), result.carbon_dioxide.root_gaseous);
    try std.testing.expectEqual(@as(f64, 2), result.carbon_dioxide.root_aqueous);
    try std.testing.expectEqual(@as(f64, 0), result.oxygen.soil_gaseous);
    try std.testing.expectEqual(@as(f64, 4), result.oxygen.soil_aqueous);
    try std.testing.expectEqual(@as(f64, 2.4), result.root_non_band_aqueous_volume_m3);
    try std.testing.expectEqual(@as(f64, 1.6), result.root_band_aqueous_volume_m3);
    try std.testing.expectEqual(@as(f64, 2), result.oxygen_demand_per_plant_g_o_per_step);
}

test "UPTAKE signed methane nitrous ammonia hydrogen pools are not floored" {
    const result = try calculate(sourceInputs());
    try std.testing.expectEqual(@as(f64, -1), result.methane.root_gaseous);
    try std.testing.expectEqual(@as(f64, -1.5), result.nitrous_oxide.soil_aqueous);
    try std.testing.expectEqual(@as(f64, -2.5), result.ammonia.soil_band_aqueous);
    try std.testing.expectEqual(@as(f64, -3.25), result.hydrogen.soil_aqueous);
    try std.testing.expectEqual(@as(f64, -0.5), result.allocated_soil_oxygen_gas_flux_g_o_per_step);
}

test "zero plant population fails explicitly" {
    var inputs = sourceInputs();
    inputs.plant_population = 0;
    try std.testing.expectError(
        error.InvalidRootGasPoolInitializationInput,
        calculate(inputs),
    );
}
