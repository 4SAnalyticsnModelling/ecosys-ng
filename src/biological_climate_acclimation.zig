const std = @import("std");

pub const CellState = struct {
    mean_annual_air_temperature_c: f64,
    mean_annual_soil_temperature_c: f64,
    microbial_arrhenius_offset_c: f64,
};

pub const PlantInputs = struct {
    initial_thermal_adaptation_zone: []const f64,
    default_lower_temperature_threshold_c: []const f64,
    default_upper_temperature_threshold_c: []const f64,
    initial_floral_initiation_group: []const f64,
    planting_node_number: []const f64,
    carboxylation_type: []const u8,
    biomass_turnover_type: []const u8,
};

pub const PlantState = struct {
    thermal_adaptation_zone: []f64,
    arrhenius_offset_c: []f64,
    lower_temperature_threshold_c: []f64,
    upper_temperature_threshold_c: []f64,
    grain_number_heat_threshold_c: []f64,
    floral_initiation_node_number: []f64,
};

pub const Inputs = struct {
    incremental_climate_change: bool,
    source_hour: u8,
    daily_average_air_temperature_change_c: f64,
    initial_mean_annual_air_temperature_c: f64,
};

/// Atomic WTHR J=1 biological acclimation over arbitrary runtime plants.
pub fn apply(
    cell: *CellState,
    plants: PlantState,
    plant_inputs: PlantInputs,
    inputs: Inputs,
) !bool {
    if (inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidClimateAcclimationSourceHour;
    if (!inputs.incremental_climate_change or inputs.source_hour != 1)
        return false;
    inline for (@typeInfo(CellState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(cell.*, field.name)))
            return error.NonFiniteClimateAcclimationState;
    inline for (.{
        inputs.daily_average_air_temperature_change_c,
        inputs.initial_mean_annual_air_temperature_c,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteClimateAcclimationInput;
    const count = plant_inputs.initial_thermal_adaptation_zone.len;
    if (count == 0 or
        plant_inputs.default_lower_temperature_threshold_c.len != count or
        plant_inputs.default_upper_temperature_threshold_c.len != count or
        plant_inputs.initial_floral_initiation_group.len != count or
        plant_inputs.planting_node_number.len != count or
        plant_inputs.carboxylation_type.len != count or
        plant_inputs.biomass_turnover_type.len != count or
        plants.thermal_adaptation_zone.len != count or
        plants.arrhenius_offset_c.len != count or
        plants.lower_temperature_threshold_c.len != count or
        plants.upper_temperature_threshold_c.len != count or
        plants.grain_number_heat_threshold_c.len != count or
        plants.floral_initiation_node_number.len != count)
        return error.ClimateAcclimationPlantDimensionMismatch;

    const soil_change_c =
        0.5 * inputs.daily_average_air_temperature_change_c;
    const next_air_c = inputs.initial_mean_annual_air_temperature_c +
        inputs.daily_average_air_temperature_change_c;
    const next_soil_c = inputs.initial_mean_annual_air_temperature_c +
        soil_change_c;
    const microbial_offset_c =
        0.333 * (15 - std.math.clamp(next_soil_c, 0, 30));
    inline for (.{ soil_change_c, next_air_c, next_soil_c, microbial_offset_c }) |value| if (!std.math.isFinite(value))
        return error.ClimateAcclimationOverflow;

    // First pass validates every input and every candidate. The second pass
    // is a no-fail commit, so a bad late plant cannot partially acclimate.
    for (0..count) |plant| {
        inline for (.{
            plant_inputs.initial_thermal_adaptation_zone[plant],
            plant_inputs.default_lower_temperature_threshold_c[plant],
            plant_inputs.default_upper_temperature_threshold_c[plant],
            plant_inputs.initial_floral_initiation_group[plant],
            plant_inputs.planting_node_number[plant],
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteClimateAcclimationPlantInput;
        const candidate = derivePlant(plant_inputs, plant, inputs);
        inline for (@typeInfo(PlantCandidate).@"struct".fields) |field|
            if (!std.math.isFinite(@field(candidate, field.name)))
                return error.ClimateAcclimationOverflow;
    }

    cell.* = .{
        .mean_annual_air_temperature_c = next_air_c,
        .mean_annual_soil_temperature_c = next_soil_c,
        .microbial_arrhenius_offset_c = microbial_offset_c,
    };
    for (0..count) |plant| {
        const candidate = derivePlant(plant_inputs, plant, inputs);
        plants.thermal_adaptation_zone[plant] =
            candidate.thermal_adaptation_zone;
        plants.arrhenius_offset_c[plant] = candidate.arrhenius_offset_c;
        plants.lower_temperature_threshold_c[plant] =
            candidate.lower_temperature_threshold_c;
        plants.upper_temperature_threshold_c[plant] =
            candidate.upper_temperature_threshold_c;
        plants.grain_number_heat_threshold_c[plant] =
            candidate.grain_number_heat_threshold_c;
        plants.floral_initiation_node_number[plant] =
            candidate.floral_initiation_node_number;
    }
    return true;
}

const PlantCandidate = struct {
    thermal_adaptation_zone: f64,
    arrhenius_offset_c: f64,
    lower_temperature_threshold_c: f64,
    upper_temperature_threshold_c: f64,
    grain_number_heat_threshold_c: f64,
    floral_initiation_node_number: f64,
};

fn derivePlant(
    inputs: PlantInputs,
    plant: usize,
    climate: Inputs,
) PlantCandidate {
    const zone = inputs.initial_thermal_adaptation_zone[plant] +
        0.30 / 2.667 * climate.daily_average_air_temperature_change_c;
    const offset_c = 2.5 * (3 - zone);
    var group = inputs.initial_floral_initiation_group[plant] +
        0.30 * climate.daily_average_air_temperature_change_c;
    if (inputs.biomass_turnover_type[plant] != 0) group /= 25;
    group -= inputs.planting_node_number[plant];
    return .{
        .thermal_adaptation_zone = zone,
        .arrhenius_offset_c = offset_c,
        .lower_temperature_threshold_c = inputs.default_lower_temperature_threshold_c[plant] - offset_c,
        .upper_temperature_threshold_c = @min(
            15,
            inputs.default_upper_temperature_threshold_c[plant] - offset_c,
        ),
        .grain_number_heat_threshold_c = (if (inputs.carboxylation_type[plant] == 3)
            @as(f64, 27)
        else
            @as(f64, 30)) + 3 * zone,
        .floral_initiation_node_number = group,
    };
}

test "incremental climate acclimates arbitrary plants at source hour one" {
    const count = 7;
    var zone = [_]f64{0} ** count;
    var offset = [_]f64{0} ** count;
    var lower = [_]f64{0} ** count;
    var upper = [_]f64{0} ** count;
    var heat = [_]f64{0} ** count;
    var group = [_]f64{0} ** count;
    var cell: CellState = .{
        .mean_annual_air_temperature_c = 0,
        .mean_annual_soil_temperature_c = 0,
        .microbial_arrhenius_offset_c = 0,
    };
    const applied = try apply(&cell, .{
        .thermal_adaptation_zone = &zone,
        .arrhenius_offset_c = &offset,
        .lower_temperature_threshold_c = &lower,
        .upper_temperature_threshold_c = &upper,
        .grain_number_heat_threshold_c = &heat,
        .floral_initiation_node_number = &group,
    }, .{
        .initial_thermal_adaptation_zone = &.{ 1, 2, 3, 4, 5, 6, 7 },
        .default_lower_temperature_threshold_c = &.{ 0, 0, 0, 0, 0, 0, 0 },
        .default_upper_temperature_threshold_c = &.{ 30, 30, 30, 30, 30, 30, 30 },
        .initial_floral_initiation_group = &.{ 10, 10, 10, 10, 10, 10, 10 },
        .planting_node_number = &.{ 1, 1, 1, 1, 1, 1, 1 },
        .carboxylation_type = &.{ 3, 4, 3, 4, 3, 4, 3 },
        .biomass_turnover_type = &.{ 0, 1, 0, 1, 0, 1, 0 },
    }, .{
        .incremental_climate_change = true,
        .source_hour = 1,
        .daily_average_air_temperature_change_c = 4,
        .initial_mean_annual_air_temperature_c = 10,
    });
    try std.testing.expect(applied);
    try std.testing.expectEqual(@as(f64, 14), cell.mean_annual_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 12), cell.mean_annual_soil_temperature_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.999), cell.microbial_arrhenius_offset_c, 1e-15);
    try std.testing.expect(heat[0] < heat[1]);
    try std.testing.expect(group[0] > group[1]);
    try std.testing.expectEqual(@as(usize, count), zone.len);
}

test "upper plant temperature threshold retains source fifteen degree cap" {
    var zone = [_]f64{0};
    var offset = [_]f64{0};
    var lower = [_]f64{0};
    var upper = [_]f64{0};
    var heat = [_]f64{0};
    var group = [_]f64{0};
    var cell: CellState = .{
        .mean_annual_air_temperature_c = 0,
        .mean_annual_soil_temperature_c = 0,
        .microbial_arrhenius_offset_c = 0,
    };
    _ = try apply(&cell, .{
        .thermal_adaptation_zone = &zone,
        .arrhenius_offset_c = &offset,
        .lower_temperature_threshold_c = &lower,
        .upper_temperature_threshold_c = &upper,
        .grain_number_heat_threshold_c = &heat,
        .floral_initiation_node_number = &group,
    }, .{
        .initial_thermal_adaptation_zone = &.{3},
        .default_lower_temperature_threshold_c = &.{0},
        .default_upper_temperature_threshold_c = &.{100},
        .initial_floral_initiation_group = &.{10},
        .planting_node_number = &.{0},
        .carboxylation_type = &.{3},
        .biomass_turnover_type = &.{0},
    }, .{
        .incremental_climate_change = true,
        .source_hour = 1,
        .daily_average_air_temperature_change_c = 0,
        .initial_mean_annual_air_temperature_c = 10,
    });
    try std.testing.expectEqual(@as(f64, 15), upper[0]);
}

test "nonincremental or nonfirst hour leaves all state untouched" {
    var value = [_]f64{9};
    var cell: CellState = .{
        .mean_annual_air_temperature_c = 1,
        .mean_annual_soil_temperature_c = 2,
        .microbial_arrhenius_offset_c = 3,
    };
    const before = cell;
    const applied = try apply(&cell, .{
        .thermal_adaptation_zone = &value,
        .arrhenius_offset_c = &value,
        .lower_temperature_threshold_c = &value,
        .upper_temperature_threshold_c = &value,
        .grain_number_heat_threshold_c = &value,
        .floral_initiation_node_number = &value,
    }, .{
        .initial_thermal_adaptation_zone = &.{1},
        .default_lower_temperature_threshold_c = &.{1},
        .default_upper_temperature_threshold_c = &.{1},
        .initial_floral_initiation_group = &.{1},
        .planting_node_number = &.{1},
        .carboxylation_type = &.{3},
        .biomass_turnover_type = &.{0},
    }, .{
        .incremental_climate_change = true,
        .source_hour = 2,
        .daily_average_air_temperature_change_c = 1,
        .initial_mean_annual_air_temperature_c = 10,
    });
    try std.testing.expect(!applied);
    try std.testing.expectEqualDeep(before, cell);
    try std.testing.expectEqual(@as(f64, 9), value[0]);
}

test "invalid late plant rolls back cell and all plant outputs" {
    var zone = [_]f64{ 7, 8 };
    var offset = [_]f64{ 7, 8 };
    var lower = [_]f64{ 7, 8 };
    var upper = [_]f64{ 7, 8 };
    var heat = [_]f64{ 7, 8 };
    var group = [_]f64{ 7, 8 };
    const zone_before = zone;
    var cell: CellState = .{
        .mean_annual_air_temperature_c = 1,
        .mean_annual_soil_temperature_c = 2,
        .microbial_arrhenius_offset_c = 3,
    };
    const cell_before = cell;
    try std.testing.expectError(
        error.NonFiniteClimateAcclimationPlantInput,
        apply(&cell, .{
            .thermal_adaptation_zone = &zone,
            .arrhenius_offset_c = &offset,
            .lower_temperature_threshold_c = &lower,
            .upper_temperature_threshold_c = &upper,
            .grain_number_heat_threshold_c = &heat,
            .floral_initiation_node_number = &group,
        }, .{
            .initial_thermal_adaptation_zone = &.{ 1, std.math.nan(f64) },
            .default_lower_temperature_threshold_c = &.{ 0, 0 },
            .default_upper_temperature_threshold_c = &.{ 20, 20 },
            .initial_floral_initiation_group = &.{ 1, 1 },
            .planting_node_number = &.{ 0, 0 },
            .carboxylation_type = &.{ 3, 3 },
            .biomass_turnover_type = &.{ 0, 0 },
        }, .{
            .incremental_climate_change = true,
            .source_hour = 1,
            .daily_average_air_temperature_change_c = 1,
            .initial_mean_annual_air_temperature_c = 10,
        }),
    );
    try std.testing.expectEqualDeep(cell_before, cell);
    try std.testing.expectEqualSlices(f64, &zone_before, &zone);
}
