const std = @import("std");

pub const source_pool_count: usize = 5;

pub const Parameters = struct {
    water_retention_m3_per_g_c: [source_pool_count]f64,
    dry_bulk_density_Mg_per_m3: [source_pool_count]f64,
    dry_mass_Mg_per_g_c: f64,
    particle_density_Mg_per_m3: f64,
    field_capacity_fraction_of_porosity: f64,
    wilting_point_fraction_of_porosity: f64,
};

pub const Inputs = struct {
    carbon_by_pool_g_c: [source_pool_count]f64,
    charcoal_carbon_g_c: f64,
    water_m3: f64,
    ice_m3: f64,
};

pub const Result = struct {
    water_retention_capacity_m3: f64,
    dry_litter_volume_m3: f64,
    expanded_total_volume_m3: f64,
    dry_mass_Mg: f64,
    pore_volume_m3: f64,
    air_volume_m3: f64,
    porosity_m3_per_m3: f64,
    field_capacity_m3_per_m3: f64,
    wilting_point_m3_per_m3: f64,
};

/// HOUR1 THETY/THETZ inverse of the log-log segment between field capacity
/// and wilting point. Used for hygroscopic and minimum-potential water limits.
pub fn waterFractionAtPotentialBelowWilting(field_capacity_fraction: f64, wilting_point_fraction: f64, field_capacity_potential_mpa: f64, wilting_point_potential_mpa: f64, target_potential_mpa: f64) !f64 {
    inline for (.{ field_capacity_fraction, wilting_point_fraction, field_capacity_potential_mpa, wilting_point_potential_mpa, target_potential_mpa }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceLitterRetentionInput;
    if (field_capacity_fraction <= wilting_point_fraction or wilting_point_fraction <= 0 or field_capacity_potential_mpa >= 0 or wilting_point_potential_mpa >= field_capacity_potential_mpa or target_potential_mpa >= wilting_point_potential_mpa) return error.InvalidSurfaceLitterRetentionInput;
    const log_field_water = @log(field_capacity_fraction);
    const log_wilting_water = @log(wilting_point_fraction);
    const log_field_potential = @log(-field_capacity_potential_mpa);
    const log_wilting_potential = @log(-wilting_point_potential_mpa);
    const log_target_potential = @log(-target_potential_mpa);
    const result = @exp((log_field_potential - log_target_potential) * (log_field_water - log_wilting_water) / (log_wilting_potential - log_field_potential) + log_field_water);
    if (!std.math.isFinite(result) or result < 0) return error.NonFiniteSurfaceLitterRetention;
    return result;
}

/// HOUR1 L=0 geometry. Pool 3 is intentionally excluded from VOLWRX/VOLR,
/// matching RC0 indices 0,1,2,4 in the source.
pub fn calculate(inputs: Inputs, parameters: Parameters) !Result {
    inline for (inputs.carbon_by_pool_g_c ++ parameters.water_retention_m3_per_g_c ++ parameters.dry_bulk_density_Mg_per_m3) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceLitterGeometryInput;
    inline for (.{ inputs.charcoal_carbon_g_c, inputs.water_m3, inputs.ice_m3, parameters.dry_mass_Mg_per_g_c, parameters.particle_density_Mg_per_m3, parameters.field_capacity_fraction_of_porosity, parameters.wilting_point_fraction_of_porosity }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceLitterGeometryInput;
    for (inputs.carbon_by_pool_g_c) |value| if (value < 0) return error.InvalidSurfaceLitterGeometryInput;
    for (parameters.water_retention_m3_per_g_c) |value| if (value < 0) return error.InvalidSurfaceLitterGeometryInput;
    for (parameters.dry_bulk_density_Mg_per_m3) |value| if (value <= 0) return error.InvalidSurfaceLitterGeometryInput;
    if (inputs.charcoal_carbon_g_c < 0 or inputs.water_m3 < 0 or inputs.ice_m3 < 0 or parameters.dry_mass_Mg_per_g_c < 0 or parameters.particle_density_Mg_per_m3 <= 0 or parameters.field_capacity_fraction_of_porosity < 0 or parameters.field_capacity_fraction_of_porosity > 1 or parameters.wilting_point_fraction_of_porosity < 0 or parameters.wilting_point_fraction_of_porosity > parameters.field_capacity_fraction_of_porosity) return error.InvalidSurfaceLitterGeometryInput;

    const included = [_]usize{ 0, 1, 2, 4 };
    var retention: f64 = 0;
    var dry_volume: f64 = 0;
    for (included) |pool| {
        retention += parameters.water_retention_m3_per_g_c[pool] * inputs.carbon_by_pool_g_c[pool];
        dry_volume += 1e-6 * inputs.carbon_by_pool_g_c[pool] / parameters.dry_bulk_density_Mg_per_m3[pool];
    }
    const excess_water_and_ice = @max(0, inputs.water_m3 + inputs.ice_m3 - retention);
    const expanded_volume = dry_volume + excess_water_and_ice;
    var total_carbon: f64 = inputs.charcoal_carbon_g_c;
    for (inputs.carbon_by_pool_g_c) |value| total_carbon += value;
    const dry_mass = parameters.dry_mass_Mg_per_g_c * total_carbon;
    const pore_volume = @max(0, dry_volume - dry_mass / parameters.particle_density_Mg_per_m3);
    const air_volume = @max(0, pore_volume - inputs.water_m3 - inputs.ice_m3);
    const porosity = if (dry_volume > 0) pore_volume / dry_volume else parameters.water_retention_m3_per_g_c[1] / parameters.dry_bulk_density_Mg_per_m3[1];
    const charcoal_volume_fraction = if (dry_volume > 0) 1e-6 * inputs.charcoal_carbon_g_c / dry_volume else 0;
    const result: Result = .{
        .water_retention_capacity_m3 = retention,
        .dry_litter_volume_m3 = dry_volume,
        .expanded_total_volume_m3 = expanded_volume,
        .dry_mass_Mg = dry_mass,
        .pore_volume_m3 = pore_volume,
        .air_volume_m3 = air_volume,
        .porosity_m3_per_m3 = porosity,
        .field_capacity_m3_per_m3 = parameters.field_capacity_fraction_of_porosity * porosity + charcoal_volume_fraction,
        .wilting_point_m3_per_m3 = parameters.wilting_point_fraction_of_porosity * porosity + charcoal_volume_fraction,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.NonFiniteSurfaceLitterGeometry;
    return result;
}

test "surface litter geometry reproduces HOUR1 pool equations" {
    const carbon = [source_pool_count]f64{ 10, 20, 30, 40, 50 };
    const parameters: Parameters = .{ .water_retention_m3_per_g_c = .{ 2e-6, 5e-6, 5e-6, 5e-6, 5e-6 }, .dry_bulk_density_Mg_per_m3 = .{ 0.1, 0.0125, 0.025, 0.025, 0.025 }, .dry_mass_Mg_per_g_c = 1.82e-6, .particle_density_Mg_per_m3 = 1.30, .field_capacity_fraction_of_porosity = 0.5, .wilting_point_fraction_of_porosity = 0.25 };
    const result = try calculate(.{ .carbon_by_pool_g_c = carbon, .charcoal_carbon_g_c = 5, .water_m3 = 0.001, .ice_m3 = 0 }, parameters);
    const retention = 2e-6 * 10 + 5e-6 * (20 + 30 + 50);
    const volume: f64 = 1e-6 * (10.0 / 0.1 + 20.0 / 0.0125 + 30.0 / 0.025 + 50.0 / 0.025);
    try std.testing.expectApproxEqAbs(retention, result.water_retention_capacity_m3, 1e-18);
    try std.testing.expectApproxEqAbs(volume, result.dry_litter_volume_m3, 1e-18);
    try std.testing.expectApproxEqAbs(volume + @max(0, 0.001 - retention), result.expanded_total_volume_m3, 1e-18);
    try std.testing.expectApproxEqAbs(@max(0, result.pore_volume_m3 - 0.001), result.air_volume_m3, 1e-18);
}

test "surface litter hygroscopic water reproduces HOUR1 THETY" {
    const value = try waterFractionAtPotentialBelowWilting(0.3, 0.1, -0.033, -1.5, -1.5e4);
    const expected = @exp((@log(0.033) - @log(1.5e4)) * (@log(0.3) - @log(0.1)) / (@log(1.5) - @log(0.033)) + @log(0.3));
    try std.testing.expectApproxEqAbs(expected, value, 1e-15);
}
