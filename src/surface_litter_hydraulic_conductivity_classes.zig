const std = @import("std");

pub const Inputs = struct {
    porosity_m3_m3: f64,
    field_capacity_m3_m3: f64,
    wilting_point_m3_m3: f64,
    log_porosity: f64,
    log_field_capacity: f64,
    log_wilting_point: f64,
    log_porosity_to_field_capacity: f64,
    log_field_capacity_to_wilting_point: f64,
    log_saturation_potential: f64,
    log_field_capacity_potential: f64,
    log_wilting_potential: f64,
    saturation_to_field_potential_interval: f64,
    field_to_wilting_potential_interval: f64,
    dry_end_shape: f64,
    wet_end_shape: f64,
    pore_interaction: f64,
    minimum_potential_mpa: f64,
    saturation_potential_mpa: f64,
    saturated_conductivity_m2_h_mpa: f64,
    air_entry_conductivity_fraction: f64,
};

pub const Outputs = struct {
    water_content_m3_m3: []f64,
    water_potential_mpa: []f64,
    vertical_conductivity_m2_h_mpa: []f64,
    lateral_x_conductivity_m2_h_mpa: []f64,
    lateral_y_conductivity_m2_h_mpa: []f64,
};

pub const AirEntry = struct {
    potential_mpa: f64,
    water_content_m3_m3: f64,
};

/// HOUR1 lines 1986--2034. Runtime hydraulic classes replace source JK
/// while retaining both K loops, inner M summation, and threshold crossing.
pub fn compute(inputs: Inputs, outputs: Outputs) !?AirEntry {
    const class_count = try validate(inputs, outputs);
    var denominator_sum: f64 = 0;
    for (0..class_count) |class| {
        const source_k = @as(f64, @floatFromInt(class + 1));
        const count = @as(f64, @floatFromInt(class_count));
        const x_k = @max(0.0, source_k - 0.01 * count);
        const water_content =
            inputs.porosity_m3_m3 -
            (x_k / count * inputs.porosity_m3_m3);
        outputs.water_content_m3_m3[class] = water_content;
        const potential = if (water_content < inputs.wilting_point_m3_m3)
            @max(
                inputs.minimum_potential_mpa,
                -@exp(inputs.log_wilting_potential +
                    inputs.dry_end_shape *
                        ((inputs.log_wilting_point - @log(water_content)) /
                            inputs.log_field_capacity_to_wilting_point *
                            inputs.field_to_wilting_potential_interval)),
            )
        else if (water_content < inputs.field_capacity_m3_m3)
            -@exp(inputs.log_field_capacity_potential +
                ((inputs.log_field_capacity - @log(water_content)) /
                    inputs.log_field_capacity_to_wilting_point *
                    inputs.field_to_wilting_potential_interval))
        else if (water_content < inputs.porosity_m3_m3)
            @max(
                inputs.minimum_potential_mpa,
                -@exp(inputs.log_saturation_potential +
                    std.math.pow(
                        f64,
                        @max(
                            0.0,
                            inputs.log_porosity - @log(water_content),
                        ) / inputs.log_porosity_to_field_capacity,
                        inputs.wet_end_shape,
                    ) * inputs.saturation_to_field_potential_interval),
            )
        else
            inputs.saturation_potential_mpa;
        outputs.water_potential_mpa[class] = potential;
        denominator_sum +=
            (2.0 * source_k - 1.0) / (potential * potential);
    }
    if (!std.math.isFinite(denominator_sum) or denominator_sum <= 0)
        return error.InvalidLitterHydraulicClassDenominator;

    var air_entry: ?AirEntry = null;
    for (0..class_count) |class| {
        const source_k = @as(f64, @floatFromInt(class + 1));
        const count = @as(f64, @floatFromInt(class_count));
        const x_k = @max(0.0, source_k - 0.01 * count);
        const conductivity_shape =
            std.math.pow(f64, (count - x_k) / count, inputs.pore_interaction);
        var numerator_sum: f64 = 0;
        for (class..class_count) |summed_class| {
            const source_m = @as(f64, @floatFromInt(summed_class + 1));
            const x_m = @max(0.0, source_m - 0.01 * count);
            const potential = outputs.water_potential_mpa[summed_class];
            numerator_sum +=
                (2.0 * x_m + 1.0 - 2.0 * x_k) /
                (potential * potential);
        }
        outputs.vertical_conductivity_m2_h_mpa[class] =
            inputs.saturated_conductivity_m2_h_mpa * conductivity_shape *
            numerator_sum / denominator_sum;
        outputs.lateral_x_conductivity_m2_h_mpa[class] = 0.0;
        outputs.lateral_y_conductivity_m2_h_mpa[class] = 0.0;
        if (class > 0 and
            outputs.vertical_conductivity_m2_h_mpa[class] <
                inputs.air_entry_conductivity_fraction *
                    inputs.saturated_conductivity_m2_h_mpa and
            outputs.vertical_conductivity_m2_h_mpa[class - 1] >=
                inputs.air_entry_conductivity_fraction *
                    inputs.saturated_conductivity_m2_h_mpa)
            air_entry = .{
                .potential_mpa = outputs.water_potential_mpa[class - 1],
                .water_content_m3_m3 = outputs.water_content_m3_m3[class - 1],
            };
    }
    return air_entry;
}

fn validate(inputs: Inputs, outputs: Outputs) !usize {
    const count = outputs.water_content_m3_m3.len;
    if (count == 0 or outputs.water_potential_mpa.len != count or
        outputs.vertical_conductivity_m2_h_mpa.len != count or
        outputs.lateral_x_conductivity_m2_h_mpa.len != count or
        outputs.lateral_y_conductivity_m2_h_mpa.len != count)
        return error.LitterHydraulicClassDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteLitterHydraulicClassInput;
    if (inputs.porosity_m3_m3 <= 0 or inputs.field_capacity_m3_m3 <= 0 or
        inputs.wilting_point_m3_m3 <= 0 or
        inputs.log_porosity_to_field_capacity <= 0 or
        inputs.log_field_capacity_to_wilting_point <= 0 or
        inputs.minimum_potential_mpa >= 0 or
        inputs.saturation_potential_mpa >= 0 or
        inputs.saturated_conductivity_m2_h_mpa < 0 or
        inputs.air_entry_conductivity_fraction < 0)
        return error.InvalidLitterHydraulicClassInput;
    return count;
}

test "runtime hydraulic classes produce finite vertical and zero lateral conductivity" {
    var water: [10]f64 = undefined;
    var potential: [10]f64 = undefined;
    var vertical: [10]f64 = undefined;
    var lateral_x: [10]f64 = undefined;
    var lateral_y: [10]f64 = undefined;
    _ = try compute(.{
        .porosity_m3_m3 = 0.8,
        .field_capacity_m3_m3 = 0.4,
        .wilting_point_m3_m3 = 0.2,
        .log_porosity = @log(0.8),
        .log_field_capacity = @log(0.4),
        .log_wilting_point = @log(0.2),
        .log_porosity_to_field_capacity = @log(0.8) - @log(0.4),
        .log_field_capacity_to_wilting_point = @log(0.4) - @log(0.2),
        .log_saturation_potential = @log(0.01),
        .log_field_capacity_potential = @log(0.033),
        .log_wilting_potential = @log(1.5),
        .saturation_to_field_potential_interval = 1.2,
        .field_to_wilting_potential_interval = 3.8,
        .dry_end_shape = 0.5,
        .wet_end_shape = 0.5,
        .pore_interaction = 1.33,
        .minimum_potential_mpa = -1000,
        .saturation_potential_mpa = -0.001,
        .saturated_conductivity_m2_h_mpa = 0.01,
        .air_entry_conductivity_fraction = 0.1,
    }, .{
        .water_content_m3_m3 = &water,
        .water_potential_mpa = &potential,
        .vertical_conductivity_m2_h_mpa = &vertical,
        .lateral_x_conductivity_m2_h_mpa = &lateral_x,
        .lateral_y_conductivity_m2_h_mpa = &lateral_y,
    });
    for (vertical) |value| try std.testing.expect(std.math.isFinite(value));
    const zeros = [_]f64{0} ** 10;
    try std.testing.expectEqualSlices(f64, &zeros, &lateral_x);
    try std.testing.expectEqualSlices(f64, &zeros, &lateral_y);
}
