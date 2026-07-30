const std = @import("std");

pub const fraction_count = 4;

pub const ElementPool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const StructuralDelta = struct {
    carbon_g_c: f64,
    colonized_carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Placement = enum { surface, subsurface };

pub const Inputs = struct {
    amendment: ElementPool,
    allocated: ElementPool,
    carbon_fraction: []const f64,
    base_nitrogen_to_carbon: []const f64,
    base_phosphorus_to_carbon: []const f64,
    colonized_carbon_fraction: f64,
    negligible_carbon_g_c: f64,
    residue_density_mg_c_m3: f64,
    placement: Placement,
};

pub const Result = struct {
    surface_volume_delta_m3: f64,
};

/// HOUR1 lines 851--893. Distributes the post-microbial and post-dissolved
/// amendment remainder across the four structural litter fractions.
pub fn compute(inputs: Inputs, output: []StructuralDelta) !Result {
    try validate(inputs, output.len);
    const remaining: ElementPool = .{
        .carbon_g_c = inputs.amendment.carbon_g_c - inputs.allocated.carbon_g_c,
        .nitrogen_g_n = inputs.amendment.nitrogen_g_n - inputs.allocated.nitrogen_g_n,
        .phosphorus_g_p = inputs.amendment.phosphorus_g_p - inputs.allocated.phosphorus_g_p,
    };
    var nitrogen_to_carbon = [_]f64{0} ** fraction_count;
    var phosphorus_to_carbon = [_]f64{0} ** fraction_count;
    var nitrogen_weight_total: f64 = 0;
    var phosphorus_weight_total: f64 = 0;
    if (remaining.carbon_g_c > inputs.negligible_carbon_g_c) {
        var nitrogen_requirement: f64 = 0;
        var phosphorus_requirement: f64 = 0;
        for (0..fraction_count) |fraction| {
            nitrogen_requirement += remaining.carbon_g_c *
                inputs.carbon_fraction[fraction] *
                inputs.base_nitrogen_to_carbon[fraction];
            phosphorus_requirement += remaining.carbon_g_c *
                inputs.carbon_fraction[fraction] *
                inputs.base_phosphorus_to_carbon[fraction];
        }
        if (nitrogen_requirement == 0 or phosphorus_requirement == 0)
            return error.ZeroStructuralNutrientRequirement;
        const nitrogen_scale = remaining.nitrogen_g_n / nitrogen_requirement;
        const phosphorus_scale = remaining.phosphorus_g_p / phosphorus_requirement;
        for (0..fraction_count) |fraction| {
            nitrogen_to_carbon[fraction] =
                inputs.base_nitrogen_to_carbon[fraction] * nitrogen_scale;
            phosphorus_to_carbon[fraction] =
                inputs.base_phosphorus_to_carbon[fraction] * phosphorus_scale;
            nitrogen_weight_total +=
                inputs.carbon_fraction[fraction] * nitrogen_to_carbon[fraction];
            phosphorus_weight_total +=
                inputs.carbon_fraction[fraction] * phosphorus_to_carbon[fraction];
        }
    }

    var scratch: [fraction_count]StructuralDelta = undefined;
    var surface_volume_delta_m3: f64 = 0;
    for (0..fraction_count) |fraction| {
        const carbon = inputs.carbon_fraction[fraction] * remaining.carbon_g_c;
        const nitrogen = if (nitrogen_weight_total > 0)
            inputs.carbon_fraction[fraction] * nitrogen_to_carbon[fraction] /
                nitrogen_weight_total * remaining.nitrogen_g_n
        else
            0;
        const phosphorus = if (phosphorus_weight_total > 0)
            inputs.carbon_fraction[fraction] * phosphorus_to_carbon[fraction] /
                phosphorus_weight_total * remaining.phosphorus_g_p
        else
            0;
        scratch[fraction] = .{
            .carbon_g_c = carbon,
            .colonized_carbon_g_c = carbon * inputs.colonized_carbon_fraction,
            .nitrogen_g_n = nitrogen,
            .phosphorus_g_p = phosphorus,
        };
        if (inputs.placement == .surface)
            surface_volume_delta_m3 += carbon * 1.0e-6 / inputs.residue_density_mg_c_m3;
    }
    @memcpy(output, &scratch);
    return .{ .surface_volume_delta_m3 = surface_volume_delta_m3 };
}

fn validate(inputs: Inputs, output_len: usize) !void {
    if (output_len != fraction_count or
        inputs.carbon_fraction.len != fraction_count or
        inputs.base_nitrogen_to_carbon.len != fraction_count or
        inputs.base_phosphorus_to_carbon.len != fraction_count)
        return error.StructuralAllocationDimensionMismatch;
    inline for (.{ inputs.amendment, inputs.allocated }) |pool|
        inline for (@typeInfo(ElementPool).@"struct".fields) |field| {
            const value = @field(pool, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteStructuralAllocationInput;
            if (value < 0) return error.InvalidStructuralAllocationInput;
        };
    if (inputs.allocated.carbon_g_c > inputs.amendment.carbon_g_c or
        inputs.allocated.nitrogen_g_n > inputs.amendment.nitrogen_g_n or
        inputs.allocated.phosphorus_g_p > inputs.amendment.phosphorus_g_p)
        return error.OrganicAllocationExceedsAmendment;
    inline for (.{
        inputs.carbon_fraction,
        inputs.base_nitrogen_to_carbon,
        inputs.base_phosphorus_to_carbon,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStructuralAllocationInput;
    inline for (.{ inputs.colonized_carbon_fraction, inputs.negligible_carbon_g_c }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidStructuralAllocationInput;
    if (inputs.placement == .surface and
        (!std.math.isFinite(inputs.residue_density_mg_c_m3) or
            inputs.residue_density_mg_c_m3 <= 0))
        return error.InvalidResidueDensity;
}

test "structural remainder preserves source weighting and surface volume" {
    var output: [fraction_count]StructuralDelta = undefined;
    const result = try compute(.{
        .amendment = .{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 1 },
        .allocated = .{ .carbon_g_c = 20, .nitrogen_g_n = 2, .phosphorus_g_p = 0.2 },
        .carbon_fraction = &.{ 0.1, 0.2, 0.3, 0.4 },
        .base_nitrogen_to_carbon = &.{ 0.02, 0.02, 0.02, 0.02 },
        .base_phosphorus_to_carbon = &.{ 0.002, 0.002, 0.002, 0.002 },
        .colonized_carbon_fraction = 0.05,
        .negligible_carbon_g_c = 1.0e-12,
        .residue_density_mg_c_m3 = 0.1,
        .placement = .surface,
    }, &output);
    try std.testing.expectEqual(@as(f64, 8), output[0].carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), output[0].nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.08), output[0].phosphorus_g_p, 1e-14);
    try std.testing.expectEqual(@as(f64, 0.4), output[0].colonized_carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0008), result.surface_volume_delta_m3, 1e-16);
}

test "below threshold still distributes carbon but no nutrients" {
    var output: [fraction_count]StructuralDelta = undefined;
    const result = try compute(.{
        .amendment = .{ .carbon_g_c = 1, .nitrogen_g_n = 1, .phosphorus_g_p = 1 },
        .allocated = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .carbon_fraction = &.{ 0.25, 0.25, 0.25, 0.25 },
        .base_nitrogen_to_carbon = &.{ 0, 0, 0, 0 },
        .base_phosphorus_to_carbon = &.{ 0, 0, 0, 0 },
        .colonized_carbon_fraction = 0,
        .negligible_carbon_g_c = 1,
        .residue_density_mg_c_m3 = 0,
        .placement = .subsurface,
    }, &output);
    try std.testing.expectEqual(@as(f64, 0.25), output[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), output[0].nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 0), output[0].phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 0), result.surface_volume_delta_m3);
}

test "zero nutrient requirement fails before changing output" {
    var output = [_]StructuralDelta{.{ .carbon_g_c = 42, .colonized_carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }} ** fraction_count;
    try std.testing.expectError(error.ZeroStructuralNutrientRequirement, compute(.{
        .amendment = .{ .carbon_g_c = 1, .nitrogen_g_n = 1, .phosphorus_g_p = 1 },
        .allocated = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .carbon_fraction = &.{ 0.25, 0.25, 0.25, 0.25 },
        .base_nitrogen_to_carbon = &.{ 0, 0, 0, 0 },
        .base_phosphorus_to_carbon = &.{ 0, 0, 0, 0 },
        .colonized_carbon_fraction = 0,
        .negligible_carbon_g_c = 0,
        .residue_density_mg_c_m3 = 1,
        .placement = .surface,
    }, &output));
    try std.testing.expectEqual(@as(f64, 42), output[0].carbon_g_c);
}
