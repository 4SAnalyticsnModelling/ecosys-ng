const std = @import("std");

pub const GasStateSource = enum {
    atmospheric_equilibrium,
    supplied_profile,
};

pub const Control = struct {
    gas_state_source: GasStateSource,
    gas_initialization_index: usize, // IGO
};

/// Atmospheric gas concentrations are `g m-3`.
pub const AtmosphericConcentrations = struct {
    carbon_dioxide_g_per_m3: f64,
    methane_g_per_m3: f64,
    oxygen_g_per_m3: f64,
    dinitrogen_g_per_m3: f64,
    nitrous_oxide_g_per_m3: f64,
    ammonia_g_per_m3: f64,
    hydrogen_g_per_m3: f64,
};

/// Solubility reference factors and exponential coefficients retain the
/// dimensions used by STARTE.F's temperature correction.
pub const SolubilityParameters = struct {
    oxygen_reference: f64, // SOXYX
    oxygen_exponential_coefficient: f64, // AOXYX
    carbon_dioxide_reference: f64, // SCO2X
    carbon_dioxide_exponential_coefficient: f64, // ACO2X
    methane_reference: f64, // SCH4X
    methane_exponential_coefficient: f64, // ACH4X
    dinitrogen_reference: f64, // SN2GX
    dinitrogen_exponential_coefficient: f64, // AN2GX
    nitrous_oxide_reference: f64, // SN2OX
    nitrous_oxide_exponential_coefficient: f64, // AN2OX
    hydrogen_reference: f64, // SH2GX
    hydrogen_exponential_coefficient: f64, // AH2GX
};

pub const LayerEnvironment = struct {
    air_volume_m3: f64, // VOLP
    aqueous_volume_m3: f64, // FC
    overlying_depth_m: f64, // CDPTH(L-1)
    water_table_depth_m: f64, // DTBLZ
    air_temperature_c: f64, // ATCA
    temperature_reference_transform: f64, // CSTR1
};

/// Gas inventories are grams in the gaseous and aqueous phases.
pub const GasInventories = struct {
    gaseous_carbon_dioxide_g: f64, // CO2G
    gaseous_methane_g: f64, // CH4G
    gaseous_oxygen_g: f64, // OXYG
    gaseous_dinitrogen_g: f64, // Z2GG
    gaseous_nitrous_oxide_g: f64, // Z2OG
    gaseous_ammonia_g: f64, // ZNH3G
    gaseous_hydrogen_g: f64, // H2GG
    aqueous_oxygen_g: f64, // OXYS
    aqueous_carbon_dioxide_g: f64, // CO2S
    aqueous_methane_g: f64, // CH4S
    aqueous_dinitrogen_g: f64, // Z2GS
    aqueous_nitrous_oxide_g: f64, // Z2OS
    aqueous_hydrogen_g: f64, // H2GS
};

fn correctedDissolvedInventory(
    atmospheric_g_per_m3: f64,
    reference_solubility: f64,
    exponential_coefficient: f64,
    temperature_reference_transform: f64,
    temperature_factor: f64,
    air_temperature_c: f64,
    aqueous_volume_m3: f64,
) f64 {
    return atmospheric_g_per_m3 * reference_solubility /
        @exp(exponential_coefficient * temperature_reference_transform) *
        @exp(temperature_factor - air_temperature_c) * aqueous_volume_m3;
}

/// Direct translation of STARTE.F lines 1410--1433. The caller supplies the
/// current soil layer; `overlying_depth_m` represents source index `L-1`.
pub fn initialize(
    control: Control,
    atmosphere: AtmosphericConcentrations,
    solubility: SolubilityParameters,
    environment: LayerEnvironment,
) !?GasInventories {
    if (control.gas_state_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    inline for (@typeInfo(AtmosphericConcentrations).@"struct".fields) |field| {
        const value = @field(atmosphere, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidAtmosphericGasConcentration;
    }
    inline for (@typeInfo(SolubilityParameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(solubility, field.name))) return error.InvalidGasSolubilityParameter;
    }
    if (!std.math.isFinite(environment.air_volume_m3) or environment.air_volume_m3 < 0 or
        !std.math.isFinite(environment.aqueous_volume_m3) or environment.aqueous_volume_m3 < 0 or
        !std.math.isFinite(environment.overlying_depth_m) or
        !std.math.isFinite(environment.water_table_depth_m) or
        !std.math.isFinite(environment.air_temperature_c) or
        !std.math.isFinite(environment.temperature_reference_transform))
        return error.InvalidSoilGasEnvironment;

    const oxygen_aqueous_g = if (environment.overlying_depth_m < environment.water_table_depth_m)
        atmosphere.oxygen_g_per_m3 * solubility.oxygen_reference /
            @exp(solubility.oxygen_exponential_coefficient * environment.temperature_reference_transform) *
            @exp(0.516 - 0.0172 * environment.air_temperature_c) * environment.aqueous_volume_m3
    else
        0.0;
    const result: GasInventories = .{
        .gaseous_carbon_dioxide_g = atmosphere.carbon_dioxide_g_per_m3 * environment.air_volume_m3,
        .gaseous_methane_g = atmosphere.methane_g_per_m3 * environment.air_volume_m3,
        .gaseous_oxygen_g = atmosphere.oxygen_g_per_m3 * environment.air_volume_m3,
        .gaseous_dinitrogen_g = atmosphere.dinitrogen_g_per_m3 * environment.air_volume_m3,
        .gaseous_nitrous_oxide_g = atmosphere.nitrous_oxide_g_per_m3 * environment.air_volume_m3,
        .gaseous_ammonia_g = atmosphere.ammonia_g_per_m3 * environment.air_volume_m3,
        .gaseous_hydrogen_g = atmosphere.hydrogen_g_per_m3 * environment.air_volume_m3,
        .aqueous_oxygen_g = oxygen_aqueous_g,
        .aqueous_carbon_dioxide_g = correctedDissolvedInventory(atmosphere.carbon_dioxide_g_per_m3, solubility.carbon_dioxide_reference, solubility.carbon_dioxide_exponential_coefficient, environment.temperature_reference_transform, 0.843, 0.0281 * environment.air_temperature_c, environment.aqueous_volume_m3),
        .aqueous_methane_g = correctedDissolvedInventory(atmosphere.methane_g_per_m3, solubility.methane_reference, solubility.methane_exponential_coefficient, environment.temperature_reference_transform, 0.597, 0.0199 * environment.air_temperature_c, environment.aqueous_volume_m3),
        .aqueous_dinitrogen_g = correctedDissolvedInventory(atmosphere.dinitrogen_g_per_m3, solubility.dinitrogen_reference, solubility.dinitrogen_exponential_coefficient, environment.temperature_reference_transform, 0.456, 0.0152 * environment.air_temperature_c, environment.aqueous_volume_m3),
        .aqueous_nitrous_oxide_g = correctedDissolvedInventory(atmosphere.nitrous_oxide_g_per_m3, solubility.nitrous_oxide_reference, solubility.nitrous_oxide_exponential_coefficient, environment.temperature_reference_transform, 0.897, 0.0299 * environment.air_temperature_c, environment.aqueous_volume_m3),
        .aqueous_hydrogen_g = correctedDissolvedInventory(atmosphere.hydrogen_g_per_m3, solubility.hydrogen_reference, solubility.hydrogen_exponential_coefficient, environment.temperature_reference_transform, 0.597, 0.0199 * environment.air_temperature_c, environment.aqueous_volume_m3),
    };
    inline for (@typeInfo(GasInventories).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidInitialSoilGasInventory;
    }
    return result;
}

fn filledAtmosphere(value: f64) AtmosphericConcentrations {
    var result: AtmosphericConcentrations = undefined;
    inline for (@typeInfo(AtmosphericConcentrations).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

fn filledSolubility(reference: f64, coefficient: f64) SolubilityParameters {
    var result: SolubilityParameters = undefined;
    inline for (@typeInfo(SolubilityParameters).@"struct".fields) |field|
        @field(result, field.name) = if (std.mem.endsWith(u8, field.name, "reference")) reference else coefficient;
    return result;
}

test "STARTE soil gas initialization preserves gas and temperature calculation order" {
    const result = (try initialize(.{ .gas_state_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filledAtmosphere(2), filledSolubility(3, 0.2), .{
        .air_volume_m3 = 4,
        .aqueous_volume_m3 = 5,
        .overlying_depth_m = 1,
        .water_table_depth_m = 2,
        .air_temperature_c = 10,
        .temperature_reference_transform = 0.5,
    })).?;
    try std.testing.expectEqual(@as(f64, 8), result.gaseous_carbon_dioxide_g);
    const expected_co2 = 2.0 * 3.0 / @exp(0.2 * 0.5) * @exp(0.843 - 0.0281 * 10.0) * 5.0;
    try std.testing.expectApproxEqRel(expected_co2, result.aqueous_carbon_dioxide_g, 1e-14);
    try std.testing.expect(result.aqueous_oxygen_g > 0);
}

test "STARTE submerged soil layer initializes aqueous oxygen to exact zero" {
    const result = (try initialize(.{ .gas_state_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filledAtmosphere(1), filledSolubility(1, 0), .{
        .air_volume_m3 = 1,
        .aqueous_volume_m3 = 1,
        .overlying_depth_m = 2,
        .water_table_depth_m = 2,
        .air_temperature_c = 0,
        .temperature_reference_transform = 0,
    })).?;
    try std.testing.expectEqual(@as(f64, 0), result.aqueous_oxygen_g);
}

test "STARTE inactive soil gas initialization does not inspect dormant invalid inputs" {
    try std.testing.expectEqual(@as(?GasInventories, null), try initialize(.{ .gas_state_source = .supplied_profile, .gas_initialization_index = 0 }, filledAtmosphere(std.math.nan(f64)), filledSolubility(std.math.nan(f64), std.math.nan(f64)), .{
        .air_volume_m3 = std.math.nan(f64),
        .aqueous_volume_m3 = std.math.nan(f64),
        .overlying_depth_m = std.math.nan(f64),
        .water_table_depth_m = std.math.nan(f64),
        .air_temperature_c = std.math.nan(f64),
        .temperature_reference_transform = std.math.nan(f64),
    }));
}
