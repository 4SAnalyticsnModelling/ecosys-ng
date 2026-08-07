const std = @import("std");

pub const RetentionInputMode = enum {
    use_supplied,
    estimate_missing,
};

pub const Inputs = struct {
    mode: RetentionInputMode,
    porosity_m3_m3: []const f64,
    field_capacity_was_missing: []const bool,
    wilting_point_was_missing: []const bool,
    sand_fraction: []const f64,
    clay_fraction: []const f64,
    organic_carbon_g_per_megagram: []const f64,
    bulk_density_megagrams_m3: []const f64,
    coarse_fragment_fraction: []const f64,
    organic_soil_threshold_g_per_megagram: f64,
};

pub const Outputs = struct {
    field_capacity_m3_m3: []f64,
    wilting_point_m3_m3: []f64,
    wet_end_shape: []f64,
    pore_interaction: []f64,
    dry_end_shape: []f64,
    log_porosity: []f64,
    log_field_capacity: []f64,
    log_wilting_point: []f64,
    log_porosity_to_field_capacity: []f64,
    log_field_capacity_to_wilting_point: []f64,
};

/// `hour1.f` lines 2039--2116. Runtime layer slices preserve the supplied-data
/// gate and source mineral/organic default FC and WP branches.
pub fn apply(inputs: Inputs, outputs: Outputs) !void {
    const layer_count = try validate(inputs, outputs);
    for (0..layer_count) |layer| {
        outputs.wet_end_shape[layer] = 0.50;
        outputs.pore_interaction[layer] = 1.33;
        outputs.dry_end_shape[layer] = 0.50;
        outputs.log_porosity[layer] = @log(inputs.porosity_m3_m3[layer]);
        if ((!inputs.field_capacity_was_missing[layer] and
            !inputs.wilting_point_was_missing[layer]) or
            inputs.mode == .use_supplied)
        {
            publishLogs(layer, outputs);
        } else {
            if (inputs.mode == .estimate_missing and
                (inputs.field_capacity_was_missing[layer] or
                    inputs.wilting_point_was_missing[layer]))
            {
                var field_capacity: f64 = if (inputs.organic_carbon_g_per_megagram[layer] <
                    inputs.organic_soil_threshold_g_per_megagram)
                    0.2576 - 0.20 * inputs.sand_fraction[layer] +
                        0.36 * inputs.clay_fraction[layer] +
                        0.60e-6 * inputs.organic_carbon_g_per_megagram[layer]
                else if (inputs.bulk_density_megagrams_m3[layer] < 0.075)
                    0.27
                else if (inputs.bulk_density_megagrams_m3[layer] < 0.195)
                    0.62
                else
                    0.71;
                field_capacity =
                    field_capacity / (1.0 - inputs.coarse_fragment_fraction[layer]);
                outputs.field_capacity_m3_m3[layer] = @min(
                    0.75 * inputs.porosity_m3_m3[layer],
                    field_capacity,
                );
                var wilting_point: f64 = if (inputs.organic_carbon_g_per_megagram[layer] <
                    inputs.organic_soil_threshold_g_per_megagram)
                    0.0260 + 0.50 * inputs.clay_fraction[layer] +
                        0.32e-6 * inputs.organic_carbon_g_per_megagram[layer]
                else if (inputs.bulk_density_megagrams_m3[layer] < 0.075)
                    0.04
                else if (inputs.bulk_density_megagrams_m3[layer] < 0.195)
                    0.15
                else
                    0.22;
                wilting_point =
                    wilting_point / (1.0 - inputs.coarse_fragment_fraction[layer]);
                outputs.wilting_point_m3_m3[layer] = @min(
                    0.75 * outputs.field_capacity_m3_m3[layer],
                    wilting_point,
                );
            }
            publishLogs(layer, outputs);
        }
    }
}

fn publishLogs(layer: usize, outputs: Outputs) void {
    outputs.log_field_capacity[layer] =
        @log(outputs.field_capacity_m3_m3[layer]);
    outputs.log_wilting_point[layer] =
        @log(outputs.wilting_point_m3_m3[layer]);
    outputs.log_porosity_to_field_capacity[layer] =
        outputs.log_porosity[layer] - outputs.log_field_capacity[layer];
    outputs.log_field_capacity_to_wilting_point[layer] =
        outputs.log_field_capacity[layer] - outputs.log_wilting_point[layer];
}

fn validate(inputs: Inputs, outputs: Outputs) !usize {
    const count = inputs.porosity_m3_m3.len;
    if (count == 0) return error.ZeroSoilRetentionLayerExtent;
    inline for (.{
        inputs.field_capacity_was_missing,
        inputs.wilting_point_was_missing,
        inputs.sand_fraction,
        inputs.clay_fraction,
        inputs.organic_carbon_g_per_megagram,
        inputs.bulk_density_megagrams_m3,
        inputs.coarse_fragment_fraction,
        outputs.field_capacity_m3_m3,
        outputs.wilting_point_m3_m3,
        outputs.wet_end_shape,
        outputs.pore_interaction,
        outputs.dry_end_shape,
        outputs.log_porosity,
        outputs.log_field_capacity,
        outputs.log_wilting_point,
        outputs.log_porosity_to_field_capacity,
        outputs.log_field_capacity_to_wilting_point,
    }) |values| if (values.len != count)
        return error.SoilRetentionLayerDimensionMismatch;
    if (!std.math.isFinite(inputs.organic_soil_threshold_g_per_megagram) or
        inputs.organic_soil_threshold_g_per_megagram < 0)
        return error.InvalidSoilRetentionLayerInput;
    for (0..count) |layer| {
        inline for (.{
            inputs.porosity_m3_m3[layer],
            inputs.sand_fraction[layer],
            inputs.clay_fraction[layer],
            inputs.organic_carbon_g_per_megagram[layer],
            inputs.bulk_density_megagrams_m3[layer],
            inputs.coarse_fragment_fraction[layer],
            outputs.field_capacity_m3_m3[layer],
            outputs.wilting_point_m3_m3[layer],
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilRetentionLayerInput;
        if (inputs.porosity_m3_m3[layer] <= 0 or
            inputs.coarse_fragment_fraction[layer] >= 1)
            return error.InvalidSoilRetentionLayerInput;
    }
    return count;
}

test "missing mineral soil retention values follow source pedotransfer order" {
    var fc = [_]f64{0.1};
    var wp = [_]f64{0.05};
    var shape: [1]f64 = undefined;
    var pore: [1]f64 = undefined;
    var dry: [1]f64 = undefined;
    var lp: [1]f64 = undefined;
    var lf: [1]f64 = undefined;
    var lw: [1]f64 = undefined;
    var pd: [1]f64 = undefined;
    var fd: [1]f64 = undefined;
    try apply(.{
        .mode = .estimate_missing,
        .porosity_m3_m3 = &.{0.5},
        .field_capacity_was_missing = &.{true},
        .wilting_point_was_missing = &.{false},
        .sand_fraction = &.{0.4},
        .clay_fraction = &.{0.2},
        .organic_carbon_g_per_megagram = &.{1000},
        .bulk_density_megagrams_m3 = &.{1.2},
        .coarse_fragment_fraction = &.{0.1},
        .organic_soil_threshold_g_per_megagram = 100_000,
    }, .{
        .field_capacity_m3_m3 = &fc,
        .wilting_point_m3_m3 = &wp,
        .wet_end_shape = &shape,
        .pore_interaction = &pore,
        .dry_end_shape = &dry,
        .log_porosity = &lp,
        .log_field_capacity = &lf,
        .log_wilting_point = &lw,
        .log_porosity_to_field_capacity = &pd,
        .log_field_capacity_to_wilting_point = &fd,
    });
    try std.testing.expect(fc[0] > wp[0]);
    try std.testing.expectEqual(@as(f64, 0.5), shape[0]);
    try std.testing.expectEqual(@as(f64, 1.33), pore[0]);
    try std.testing.expect(std.math.isFinite(fd[0]));
}
