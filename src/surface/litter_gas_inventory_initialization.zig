const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

/// Atmospheric gas concentrations are `g m-3`.
pub const Atmosphere = struct {
    carbon_dioxide: f64,
    methane: f64,
    oxygen: f64,
    dinitrogen: f64,
    nitrous_oxide: f64,
    ammonia: f64,
    hydrogen: f64,
};

pub const ReferenceSolubility = struct {
    carbon_dioxide: f64,
    methane: f64,
    oxygen: f64,
    dinitrogen: f64,
    nitrous_oxide: f64,
    hydrogen: f64,
};

/// All fields are grams in the surface-litter gaseous or aqueous phase.
pub const Inventories = struct {
    gaseous_carbon_dioxide_g: f64,
    gaseous_methane_g: f64,
    gaseous_oxygen_g: f64,
    gaseous_dinitrogen_g: f64,
    gaseous_nitrous_oxide_g: f64,
    gaseous_ammonia_g: f64,
    gaseous_hydrogen_g: f64,
    aqueous_carbon_dioxide_g: f64,
    aqueous_methane_g: f64,
    aqueous_oxygen_g: f64,
    aqueous_dinitrogen_g: f64,
    aqueous_nitrous_oxide_g: f64,
    aqueous_hydrogen_g: f64,
};

/// Direct translation of `starte.f` lines 1921--1940. Surface litter is layer
/// zero; the caller supplies its runtime air and aqueous volumes.
pub fn initialize(
    control: Control,
    atmosphere: Atmosphere,
    solubility: ReferenceSolubility,
    air_volume_m3: f64,
    aqueous_volume_m3: f64,
    air_temperature_c: f64,
) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or control.gas_initialization_index != 0) return null;
    inline for (@typeInfo(Atmosphere).@"struct".fields) |field| {
        const value = @field(atmosphere, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidAtmosphericGasConcentration;
    }
    inline for (@typeInfo(ReferenceSolubility).@"struct".fields) |field| {
        const value = @field(solubility, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidGasSolubility;
    }
    if (!std.math.isFinite(air_volume_m3) or air_volume_m3 < 0 or
        !std.math.isFinite(aqueous_volume_m3) or aqueous_volume_m3 < 0 or
        !std.math.isFinite(air_temperature_c))
        return error.InvalidSurfaceLitterGasEnvironment;
    const result: Inventories = .{
        .gaseous_carbon_dioxide_g = atmosphere.carbon_dioxide * air_volume_m3,
        .gaseous_methane_g = atmosphere.methane * air_volume_m3,
        .gaseous_oxygen_g = atmosphere.oxygen * air_volume_m3,
        .gaseous_dinitrogen_g = atmosphere.dinitrogen * air_volume_m3,
        .gaseous_nitrous_oxide_g = atmosphere.nitrous_oxide * air_volume_m3,
        .gaseous_ammonia_g = atmosphere.ammonia * air_volume_m3,
        .gaseous_hydrogen_g = atmosphere.hydrogen * air_volume_m3,
        .aqueous_carbon_dioxide_g = atmosphere.carbon_dioxide * solubility.carbon_dioxide * @exp(0.843 - 0.0281 * air_temperature_c) * aqueous_volume_m3,
        .aqueous_methane_g = atmosphere.methane * solubility.methane * @exp(0.597 - 0.0199 * air_temperature_c) * aqueous_volume_m3,
        .aqueous_oxygen_g = atmosphere.oxygen * solubility.oxygen * @exp(0.516 - 0.0172 * air_temperature_c) * aqueous_volume_m3,
        .aqueous_dinitrogen_g = atmosphere.dinitrogen * solubility.dinitrogen * @exp(0.456 - 0.0152 * air_temperature_c) * aqueous_volume_m3,
        .aqueous_nitrous_oxide_g = atmosphere.nitrous_oxide * solubility.nitrous_oxide * @exp(0.897 - 0.0299 * air_temperature_c) * aqueous_volume_m3,
        .aqueous_hydrogen_g = atmosphere.hydrogen * solubility.hydrogen * @exp(0.597 - 0.0199 * air_temperature_c) * aqueous_volume_m3,
    };
    inline for (@typeInfo(Inventories).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfaceLitterGasInventory;
    }
    return result;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE surface litter gas inventories preserve gaseous and dissolved source equations" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filled(Atmosphere, 2), filled(ReferenceSolubility, 3), 4, 5, 10)).?;
    try std.testing.expectEqual(@as(f64, 8), result.gaseous_carbon_dioxide_g);
    try std.testing.expectEqual(@as(f64, 8), result.gaseous_ammonia_g);
    const expected_co2 = 2.0 * 3.0 * @exp(0.843 - 0.0281 * 10.0) * 5.0;
    const expected_oxygen = 2.0 * 3.0 * @exp(0.516 - 0.0172 * 10.0) * 5.0;
    try std.testing.expectApproxEqRel(expected_co2, result.aqueous_carbon_dioxide_g, 1e-14);
    try std.testing.expectApproxEqRel(expected_oxygen, result.aqueous_oxygen_g, 1e-14);
}

test "STARTE surface litter dissolved gas has no soil ionic-strength divisor" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filled(Atmosphere, 1), filled(ReferenceSolubility, 1), 1, 1, 0)).?;
    try std.testing.expectApproxEqRel(@exp(0.843), result.aqueous_carbon_dioxide_g, 1e-14);
}

test "STARTE inactive surface litter gas initialization ignores invalid dormant inputs" {
    const nan = std.math.nan(f64);
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, filled(Atmosphere, nan), filled(ReferenceSolubility, nan), nan, nan, nan));
}
