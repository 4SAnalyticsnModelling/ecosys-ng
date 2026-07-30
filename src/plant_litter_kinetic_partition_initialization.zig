const std = @import("std");

pub const Traits = struct {
    growth_type: i32,
    nitrogen_fixation_type: i32,
    woody_type: i32,
};

pub const Profiles = struct {
    nonstructural: []const f64,
    nonvascular: []const f64,
    legume_foliar: []const f64,
    legume_nonfoliar: []const f64,
    herbaceous_foliar: []const f64,
    herbaceous_nonfoliar: []const f64,
    deciduous_foliar: []const f64,
    woody_nonfoliar: []const f64,
    conifer_foliar: []const f64,
    herbaceous_stalk: []const f64,
    herbaceous_root: []const f64,
    deciduous_root: []const f64,
    conifer_root: []const f64,
    coarse_wood: []const f64,
    nitrogen_weights: []const f64,
    phosphorus_weights: []const f64,
};

pub const Output = struct {
    organ_count: usize,
    component_count: usize,
    carbon_fraction: []f64,
    nitrogen_fraction: []f64,
    phosphorus_fraction: []f64,
};

pub const InitializationError = error{
    InsufficientOrganCount,
    ComponentCountMismatch,
    OutputExtentMismatch,
    NonFiniteFraction,
    NegativeFraction,
    InvalidCarbonFractionSum,
    ZeroNutrientDenominator,
    ExtentOverflow,
};

/// Translates STARTQ lines 113-285.
///
/// Organ indexes retain STARTQ order: nonstructural, foliar, non-foliar,
/// stalk, fine root, and coarse wood. Kinetic-component count is runtime.
pub fn initialize(
    traits: Traits,
    profiles: Profiles,
    output: Output,
) InitializationError!void {
    if (output.organ_count < 6) return error.InsufficientOrganCount;
    const component_count = output.component_count;
    inline for (std.meta.fields(Profiles)) |field| {
        if (@field(profiles, field.name).len != component_count) {
            return error.ComponentCountMismatch;
        }
    }
    const extent = std.math.mul(usize, output.organ_count, component_count) catch
        return error.ExtentOverflow;
    if (output.carbon_fraction.len != extent or output.nitrogen_fraction.len != extent or
        output.phosphorus_fraction.len != extent)
    {
        return error.OutputExtentMismatch;
    }
    try validateProfiles(profiles);

    copyOrgan(output, 0, profiles.nonstructural);
    if (traits.growth_type == 0) {
        copyOrgan(output, 1, profiles.nonvascular);
        copyOrgan(output, 2, profiles.nonvascular);
    } else if (traits.nitrogen_fixation_type != 0) {
        copyOrgan(output, 1, profiles.legume_foliar);
        copyOrgan(output, 2, profiles.legume_nonfoliar);
    } else if (traits.woody_type == 0 or traits.growth_type <= 1) {
        copyOrgan(output, 1, profiles.herbaceous_foliar);
        copyOrgan(output, 2, profiles.herbaceous_nonfoliar);
    } else if (traits.woody_type == 1 or traits.woody_type >= 3) {
        copyOrgan(output, 1, profiles.deciduous_foliar);
        copyOrgan(output, 2, profiles.woody_nonfoliar);
    } else {
        copyOrgan(output, 1, profiles.conifer_foliar);
        copyOrgan(output, 2, profiles.woody_nonfoliar);
    }

    if (traits.growth_type == 0) {
        copyOrgan(output, 3, profiles.nonvascular);
    } else if (traits.woody_type == 0 or traits.growth_type <= 1) {
        copyOrgan(output, 3, profiles.herbaceous_stalk);
    } else {
        copyOrgan(output, 3, profiles.coarse_wood);
    }
    if (traits.growth_type == 0) {
        copyOrgan(output, 4, profiles.nonvascular);
    } else if (traits.woody_type == 0 or traits.growth_type <= 1) {
        copyOrgan(output, 4, profiles.herbaceous_root);
    } else if (traits.woody_type == 1 or traits.woody_type >= 3) {
        copyOrgan(output, 4, profiles.deciduous_root);
    } else {
        copyOrgan(output, 4, profiles.conifer_root);
    }
    copyOrgan(output, 5, profiles.coarse_wood);

    for (0..6) |organ| {
        var nitrogen_denominator: f64 = 0.0;
        var phosphorus_denominator: f64 = 0.0;
        for (0..component_count) |component| {
            const index = organ * component_count + component;
            nitrogen_denominator +=
                output.carbon_fraction[index] * profiles.nitrogen_weights[component];
            phosphorus_denominator +=
                output.carbon_fraction[index] * profiles.phosphorus_weights[component];
        }
        if (nitrogen_denominator <= 0.0 or phosphorus_denominator <= 0.0) {
            return error.ZeroNutrientDenominator;
        }
        for (0..component_count) |component| {
            const index = organ * component_count + component;
            output.nitrogen_fraction[index] =
                output.carbon_fraction[index] * profiles.nitrogen_weights[component] /
                nitrogen_denominator;
            output.phosphorus_fraction[index] =
                output.carbon_fraction[index] * profiles.phosphorus_weights[component] /
                phosphorus_denominator;
        }
    }
}

fn copyOrgan(output: Output, organ: usize, profile: []const f64) void {
    const start = organ * output.component_count;
    @memcpy(output.carbon_fraction[start .. start + output.component_count], profile);
}

fn validateProfiles(profiles: Profiles) InitializationError!void {
    inline for (std.meta.fields(Profiles)) |field| {
        const values = @field(profiles, field.name);
        var sum: f64 = 0.0;
        for (values) |value| {
            if (!std.math.isFinite(value)) return error.NonFiniteFraction;
            if (value < 0.0) return error.NegativeFraction;
            sum += value;
        }
        if (!std.mem.eql(u8, field.name, "nitrogen_weights") and
            !std.mem.eql(u8, field.name, "phosphorus_weights") and
            @abs(sum - 1.0) > 1.0e-12)
        {
            return error.InvalidCarbonFractionSum;
        }
    }
}

test "legume profiles and nutrient normalization preserve STARTQ order" {
    const profiles = legacyProfiles();
    var carbon: [24]f64 = undefined;
    var nitrogen: [24]f64 = undefined;
    var phosphorus: [24]f64 = undefined;
    try initialize(.{
        .growth_type = 2,
        .nitrogen_fixation_type = 1,
        .woody_type = 0,
    }, profiles, .{
        .organ_count = 6,
        .component_count = 4,
        .carbon_fraction = &carbon,
        .nitrogen_fraction = &nitrogen,
        .phosphorus_fraction = &phosphorus,
    });

    try std.testing.expectEqualSlices(f64, profiles.legume_foliar, carbon[4..8]);
    try std.testing.expectEqualSlices(f64, profiles.legume_nonfoliar, carbon[8..12]);
    for (0..6) |organ| {
        var nitrogen_sum: f64 = 0.0;
        var phosphorus_sum: f64 = 0.0;
        for (0..4) |component| {
            nitrogen_sum += nitrogen[organ * 4 + component];
            phosphorus_sum += phosphorus[organ * 4 + component];
        }
        try std.testing.expectApproxEqAbs(1.0, nitrogen_sum, 1.0e-14);
        try std.testing.expectApproxEqAbs(1.0, phosphorus_sum, 1.0e-14);
    }
}

test "runtime component extent is accepted" {
    const unit = [_]f64{1.0};
    const profiles = Profiles{
        .nonstructural = &unit,
        .nonvascular = &unit,
        .legume_foliar = &unit,
        .legume_nonfoliar = &unit,
        .herbaceous_foliar = &unit,
        .herbaceous_nonfoliar = &unit,
        .deciduous_foliar = &unit,
        .woody_nonfoliar = &unit,
        .conifer_foliar = &unit,
        .herbaceous_stalk = &unit,
        .herbaceous_root = &unit,
        .deciduous_root = &unit,
        .conifer_root = &unit,
        .coarse_wood = &unit,
        .nitrogen_weights = &unit,
        .phosphorus_weights = &unit,
    };
    var carbon: [6]f64 = undefined;
    var nitrogen: [6]f64 = undefined;
    var phosphorus: [6]f64 = undefined;
    try initialize(.{ .growth_type = 0, .nitrogen_fixation_type = 0, .woody_type = 0 }, profiles, .{
        .organ_count = 6,
        .component_count = 1,
        .carbon_fraction = &carbon,
        .nitrogen_fraction = &nitrogen,
        .phosphorus_fraction = &phosphorus,
    });
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 1, 1, 1, 1 }, &carbon);
}

fn legacyProfiles() Profiles {
    return .{
        .nonstructural = &.{ 0.00, 1.00, 0.00, 0.00 },
        .nonvascular = &.{ 0.07, 0.25, 0.30, 0.38 },
        .legume_foliar = &.{ 0.16, 0.38, 0.34, 0.12 },
        .legume_nonfoliar = &.{ 0.07, 0.41, 0.37, 0.15 },
        .herbaceous_foliar = &.{ 0.08, 0.41, 0.36, 0.15 },
        .herbaceous_nonfoliar = &.{ 0.07, 0.41, 0.36, 0.16 },
        .deciduous_foliar = &.{ 0.07, 0.34, 0.36, 0.23 },
        .woody_nonfoliar = &.{ 0.000, 0.045, 0.660, 0.295 },
        .conifer_foliar = &.{ 0.07, 0.25, 0.38, 0.30 },
        .herbaceous_stalk = &.{ 0.03, 0.25, 0.57, 0.15 },
        .herbaceous_root = &.{ 0.057, 0.263, 0.542, 0.138 },
        .deciduous_root = &.{ 0.059, 0.308, 0.464, 0.169 },
        .conifer_root = &.{ 0.07, 0.25, 0.38, 0.30 },
        .coarse_wood = &.{ 0.00, 0.045, 0.660, 0.295 },
        .nitrogen_weights = &.{ 0.020, 0.010, 0.010, 0.020 },
        .phosphorus_weights = &.{ 0.0020, 0.0010, 0.0010, 0.0020 },
    };
}
