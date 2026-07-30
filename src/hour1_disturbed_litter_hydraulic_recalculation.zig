const std = @import("std");

pub const Inputs = struct {
    litter_carbon_g_by_material: []const f64,
    water_capacity_m3_g_by_material: []const f64,
    bulk_density_megagrams_m3_by_material: []const f64,
    fallback_fine_material_index: usize,
    dissolved_organic_carbon_g_c: f64,
    negligible_dry_volume_m3: f64,
    log_field_capacity_potential: f64,
    log_fc_to_wilting_potential_interval: f64,
    hygroscopic_potential_mpa: f64,
    minimum_potential_mpa: f64,
};

pub const Result = struct {
    water_holding_capacity_m3: f64,
    dry_litter_volume_m3: f64,
    porosity_m3_m3: f64,
    field_capacity_m3_m3: f64,
    wilting_point_m3_m3: f64,
    log_porosity: f64,
    log_field_capacity: f64,
    log_wilting_point: f64,
    log_porosity_to_field_capacity: f64,
    log_field_capacity_to_wilting_point: f64,
    retention_shape: f64,
    pore_interaction: f64,
    dry_end_shape: f64,
    hygroscopic_water_content_m3_m3: f64,
    minimum_water_content_m3_m3: f64,
};

/// HOUR1 lines 1927--1959. Runtime material classes replace source indices
/// 0,1,2,4 while retaining class traversal and all scalar operation order.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var water_holding_capacity_m3: f64 = 0;
    for (
        inputs.water_capacity_m3_g_by_material,
        inputs.litter_carbon_g_by_material,
    ) |specific_capacity, carbon_g|
        water_holding_capacity_m3 += specific_capacity * carbon_g;
    water_holding_capacity_m3 = @max(0.0, water_holding_capacity_m3);
    var density_weighted_volume: f64 = 0;
    for (
        inputs.litter_carbon_g_by_material,
        inputs.bulk_density_megagrams_m3_by_material,
    ) |carbon_g, bulk_density|
        density_weighted_volume += carbon_g / bulk_density;
    const dry_litter_volume_m3 =
        1.0e-6 * @max(0.0, density_weighted_volume);

    var porosity_m3_m3: f64 = undefined;
    var field_capacity_m3_m3: f64 = undefined;
    var wilting_point_m3_m3: f64 = undefined;
    if (dry_litter_volume_m3 > inputs.negligible_dry_volume_m3) {
        porosity_m3_m3 =
            water_holding_capacity_m3 / dry_litter_volume_m3;
        field_capacity_m3_m3 = 0.500 * porosity_m3_m3 +
            1.0e-6 * inputs.dissolved_organic_carbon_g_c /
                dry_litter_volume_m3;
        wilting_point_m3_m3 = 0.250 * porosity_m3_m3 +
            1.0e-6 * inputs.dissolved_organic_carbon_g_c /
                dry_litter_volume_m3;
    } else {
        const fine = inputs.fallback_fine_material_index;
        porosity_m3_m3 = inputs.water_capacity_m3_g_by_material[fine] /
            inputs.bulk_density_megagrams_m3_by_material[fine];
        field_capacity_m3_m3 = 0.500 * porosity_m3_m3;
        wilting_point_m3_m3 = 0.250 * porosity_m3_m3;
    }
    if (porosity_m3_m3 <= 0 or field_capacity_m3_m3 <= 0 or
        wilting_point_m3_m3 <= 0 or
        inputs.log_fc_to_wilting_potential_interval == 0)
        return error.InvalidDisturbedLitterHydraulicResult;
    const log_porosity = @log(porosity_m3_m3);
    const log_field_capacity = @log(field_capacity_m3_m3);
    const log_wilting_point = @log(wilting_point_m3_m3);
    const log_porosity_to_field_capacity =
        log_porosity - log_field_capacity;
    const log_field_capacity_to_wilting_point =
        log_field_capacity - log_wilting_point;
    const hygroscopic_water_content_m3_m3 = @exp(
        (inputs.log_field_capacity_potential -
            @log(-inputs.hygroscopic_potential_mpa)) *
            log_field_capacity_to_wilting_point /
            inputs.log_fc_to_wilting_potential_interval +
            log_field_capacity,
    );
    const minimum_water_content_m3_m3 = @exp(
        (inputs.log_field_capacity_potential -
            @log(-inputs.minimum_potential_mpa)) *
            log_field_capacity_to_wilting_point /
            inputs.log_fc_to_wilting_potential_interval +
            log_field_capacity,
    );
    return .{
        .water_holding_capacity_m3 = water_holding_capacity_m3,
        .dry_litter_volume_m3 = dry_litter_volume_m3,
        .porosity_m3_m3 = porosity_m3_m3,
        .field_capacity_m3_m3 = field_capacity_m3_m3,
        .wilting_point_m3_m3 = wilting_point_m3_m3,
        .log_porosity = log_porosity,
        .log_field_capacity = log_field_capacity,
        .log_wilting_point = log_wilting_point,
        .log_porosity_to_field_capacity = log_porosity_to_field_capacity,
        .log_field_capacity_to_wilting_point = log_field_capacity_to_wilting_point,
        .retention_shape = 0.50,
        .pore_interaction = 1.33,
        .dry_end_shape = 0.50,
        .hygroscopic_water_content_m3_m3 = hygroscopic_water_content_m3_m3,
        .minimum_water_content_m3_m3 = minimum_water_content_m3_m3,
    };
}

fn validate(inputs: Inputs) !void {
    const count = inputs.litter_carbon_g_by_material.len;
    if (count == 0 or inputs.water_capacity_m3_g_by_material.len != count or
        inputs.bulk_density_megagrams_m3_by_material.len != count)
        return error.DisturbedLitterHydraulicDimensionMismatch;
    if (inputs.fallback_fine_material_index >= count)
        return error.FallbackLitterMaterialOutOfRange;
    inline for (.{
        inputs.litter_carbon_g_by_material,
        inputs.water_capacity_m3_g_by_material,
        inputs.bulk_density_megagrams_m3_by_material,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDisturbedLitterHydraulicInput;
    for (inputs.bulk_density_megagrams_m3_by_material) |value|
        if (value == 0) return error.InvalidDisturbedLitterHydraulicInput;
    inline for (.{
        inputs.dissolved_organic_carbon_g_c,
        inputs.negligible_dry_volume_m3,
        inputs.log_field_capacity_potential,
        inputs.log_fc_to_wilting_potential_interval,
    }) |value| if (!std.math.isFinite(value))
        return error.InvalidDisturbedLitterHydraulicInput;
    if (inputs.dissolved_organic_carbon_g_c < 0 or
        inputs.negligible_dry_volume_m3 < 0 or
        !std.math.isFinite(inputs.hygroscopic_potential_mpa) or
        !std.math.isFinite(inputs.minimum_potential_mpa) or
        inputs.hygroscopic_potential_mpa >= 0 or
        inputs.minimum_potential_mpa >= 0)
        return error.InvalidDisturbedLitterHydraulicInput;
}

test "runtime litter materials reproduce hydrologic recalculation order" {
    const result = try compute(.{
        .litter_carbon_g_by_material = &.{100},
        .water_capacity_m3_g_by_material = &.{1e-6},
        .bulk_density_megagrams_m3_by_material = &.{0.1},
        .fallback_fine_material_index = 0,
        .dissolved_organic_carbon_g_c = 0,
        .negligible_dry_volume_m3 = 1e-12,
        .log_field_capacity_potential = 1,
        .log_fc_to_wilting_potential_interval = 2,
        .hygroscopic_potential_mpa = -10,
        .minimum_potential_mpa = -100,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.0001), result.water_holding_capacity_m3, 1e-16);
    try std.testing.expectApproxEqAbs(@as(f64, 0.001), result.dry_litter_volume_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), result.porosity_m3_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), result.field_capacity_m3_m3, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.5), result.retention_shape);
    try std.testing.expect(std.math.isFinite(result.minimum_water_content_m3_m3));
}
