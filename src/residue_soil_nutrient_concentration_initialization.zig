const std = @import("std");

pub const LayerKind = enum {
    surface_litter,
    mineral_soil,
};

pub const Inputs = struct {
    residue_carbon_g_c_per_m2: []const f64,
    residue_nitrogen_g_n_per_m2: []const f64,
    residue_phosphorus_g_p_per_m2: []const f64,
    horizontal_area_m2: f64,
    dry_layer_mass_Mg: f64,
    total_layer_volume_m3: f64,
    area_scaled_calculation_floor: f64,
    calculation_floor: f64,
    layer_kind: LayerKind,
    soil_organic_carbon_g_c_per_Mg: f64,
    particulate_organic_carbon_g_c_per_Mg: f64,
    soil_organic_nitrogen_g_n_per_Mg: f64,
    soil_organic_phosphorus_g_p_per_Mg: f64,
    particulate_nitrogen_to_carbon_g_n_per_g_c: f64,
    particulate_phosphorus_to_carbon_g_p_per_g_c: f64,
};

pub const Result = struct {
    /// Runtime residue pools followed by particulate and humus slots.
    carbon_concentration: []f64,
    nitrogen_concentration: []f64,
    phosphorus_concentration: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 776--828.
pub fn initialize(result: Result, inputs: Inputs) !void {
    const residue_pool_count = inputs.residue_carbon_g_c_per_m2.len;
    if (residue_pool_count == 0 or
        inputs.residue_nitrogen_g_n_per_m2.len != residue_pool_count or
        inputs.residue_phosphorus_g_p_per_m2.len != residue_pool_count)
        return error.ResidueSoilConcentrationDimensionMismatch;
    const output_count = std.math.add(
        usize,
        residue_pool_count,
        2,
    ) catch return error.DimensionOverflow;
    if (result.carbon_concentration.len != output_count or
        result.nitrogen_concentration.len != output_count or
        result.phosphorus_concentration.len != output_count)
        return error.ResidueSoilConcentrationDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value))
                return error.NonFiniteResidueSoilConcentrationInput;
        }
    }
    if (inputs.horizontal_area_m2 <= 0 or
        inputs.dry_layer_mass_Mg < 0 or
        inputs.total_layer_volume_m3 < 0 or
        inputs.area_scaled_calculation_floor < 0 or
        inputs.calculation_floor < 0 or
        inputs.soil_organic_carbon_g_c_per_Mg < 0 or
        inputs.particulate_organic_carbon_g_c_per_Mg < 0 or
        inputs.soil_organic_nitrogen_g_n_per_Mg < 0 or
        inputs.soil_organic_phosphorus_g_p_per_Mg < 0 or
        inputs.particulate_nitrogen_to_carbon_g_n_per_g_c < 0 or
        inputs.particulate_phosphorus_to_carbon_g_p_per_g_c < 0)
        return error.InvalidResidueSoilConcentrationInput;
    const use_mass = inputs.dry_layer_mass_Mg >
        inputs.area_scaled_calculation_floor;
    if (!use_mass and inputs.total_layer_volume_m3 <= 0)
        return error.InvalidResidueSoilConcentrationInput;
    inline for (.{
        inputs.residue_carbon_g_c_per_m2,
        inputs.residue_nitrogen_g_n_per_m2,
        inputs.residue_phosphorus_g_p_per_m2,
    }) |values| {
        for (values) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteResidueSoilConcentrationInput;
            if (value < 0) return error.InvalidResidueSoilConcentrationInput;
            const candidate = value * inputs.horizontal_area_m2 /
                (if (use_mass)
                    inputs.dry_layer_mass_Mg
                else
                    inputs.total_layer_volume_m3);
            if (!std.math.isFinite(candidate))
                return error.ResidueSoilConcentrationOverflow;
        }
    }
    if (inputs.layer_kind == .mineral_soil and
        inputs.soil_organic_carbon_g_c_per_Mg > inputs.calculation_floor)
    {
        inline for (.{
            @min(
                inputs.particulate_nitrogen_to_carbon_g_n_per_g_c *
                    inputs.particulate_organic_carbon_g_c_per_Mg,
                inputs.soil_organic_nitrogen_g_n_per_Mg,
            ),
            @min(
                inputs.particulate_phosphorus_to_carbon_g_p_per_g_c *
                    inputs.particulate_organic_carbon_g_c_per_Mg,
                inputs.soil_organic_phosphorus_g_p_per_Mg,
            ),
        }) |candidate| if (!std.math.isFinite(candidate))
            return error.ResidueSoilConcentrationOverflow;
    }

    const denominator =
        if (use_mass)
            inputs.dry_layer_mass_Mg
        else
            inputs.total_layer_volume_m3;
    for (0..residue_pool_count) |pool| {
        result.carbon_concentration[pool] =
            inputs.residue_carbon_g_c_per_m2[pool] *
            inputs.horizontal_area_m2 / denominator;
    }
    for (0..residue_pool_count) |pool| {
        result.nitrogen_concentration[pool] =
            inputs.residue_nitrogen_g_n_per_m2[pool] *
            inputs.horizontal_area_m2 / denominator;
    }
    for (0..residue_pool_count) |pool| {
        result.phosphorus_concentration[pool] =
            inputs.residue_phosphorus_g_p_per_m2[pool] *
            inputs.horizontal_area_m2 / denominator;
    }
    const particulate = residue_pool_count;
    const humus = residue_pool_count + 1;
    if (inputs.layer_kind == .mineral_soil) {
        const soil_carbon = inputs.soil_organic_carbon_g_c_per_Mg;
        const particulate_carbon =
            inputs.particulate_organic_carbon_g_c_per_Mg;
        const soil_nitrogen = inputs.soil_organic_nitrogen_g_n_per_Mg;
        const soil_phosphorus = inputs.soil_organic_phosphorus_g_p_per_Mg;
        if (soil_carbon > inputs.calculation_floor) {
            result.carbon_concentration[particulate] = particulate_carbon;
            result.carbon_concentration[humus] =
                @max(0.0, soil_carbon -
                    result.carbon_concentration[particulate]);
            result.nitrogen_concentration[particulate] = @min(
                inputs.particulate_nitrogen_to_carbon_g_n_per_g_c *
                    result.carbon_concentration[particulate],
                soil_nitrogen,
            );
            result.nitrogen_concentration[humus] =
                @max(0.0, soil_nitrogen -
                    result.nitrogen_concentration[particulate]);
            result.phosphorus_concentration[particulate] = @min(
                inputs.particulate_phosphorus_to_carbon_g_p_per_g_c *
                    result.carbon_concentration[particulate],
                soil_phosphorus,
            );
            result.phosphorus_concentration[humus] =
                @max(0.0, soil_phosphorus -
                    result.phosphorus_concentration[particulate]);
        } else {
            result.carbon_concentration[particulate] = 0.0;
            result.carbon_concentration[humus] = 0.0;
            result.nitrogen_concentration[particulate] = 0.0;
            result.nitrogen_concentration[humus] = 0.0;
            result.phosphorus_concentration[particulate] = 0.0;
            result.phosphorus_concentration[humus] = 0.0;
        }
    } else {
        result.carbon_concentration[particulate] = 0.0;
        result.carbon_concentration[humus] = 0.0;
        result.nitrogen_concentration[particulate] = 0.0;
        result.nitrogen_concentration[humus] = 0.0;
        result.phosphorus_concentration[particulate] = 0.0;
        result.phosphorus_concentration[humus] = 0.0;
    }
}

test "STARTS mass branch converts residues then conserves mineral C N P" {
    var carbon = [_]f64{0.0} ** 5;
    var nitrogen = [_]f64{0.0} ** 5;
    var phosphorus = [_]f64{0.0} ** 5;
    try initialize(.{
        .carbon_concentration = &carbon,
        .nitrogen_concentration = &nitrogen,
        .phosphorus_concentration = &phosphorus,
    }, .{
        .residue_carbon_g_c_per_m2 = &.{ 1, 2, 3 },
        .residue_nitrogen_g_n_per_m2 = &.{ 0.1, 0.2, 0.3 },
        .residue_phosphorus_g_p_per_m2 = &.{ 0.01, 0.02, 0.03 },
        .horizontal_area_m2 = 100,
        .dry_layer_mass_Mg = 10,
        .total_layer_volume_m3 = 20,
        .area_scaled_calculation_floor = 1e-13,
        .calculation_floor = 1e-15,
        .layer_kind = .mineral_soil,
        .soil_organic_carbon_g_c_per_Mg = 100,
        .particulate_organic_carbon_g_c_per_Mg = 30,
        .soil_organic_nitrogen_g_n_per_Mg = 8,
        .soil_organic_phosphorus_g_p_per_Mg = 0.8,
        .particulate_nitrogen_to_carbon_g_n_per_g_c = 0.1,
        .particulate_phosphorus_to_carbon_g_p_per_g_c = 0.01,
    });
    try std.testing.expectEqualSlices(f64, &.{ 10, 20, 30, 30, 70 }, &carbon);
    try std.testing.expectEqual(@as(f64, 8), nitrogen[3] + nitrogen[4]);
    try std.testing.expectEqual(@as(f64, 0.8), phosphorus[3] + phosphorus[4]);
}

test "surface volume branch zeros particulate and humus slots" {
    var carbon = [_]f64{9.0} ** 4;
    var nitrogen = [_]f64{9.0} ** 4;
    var phosphorus = [_]f64{9.0} ** 4;
    try initialize(.{
        .carbon_concentration = &carbon,
        .nitrogen_concentration = &nitrogen,
        .phosphorus_concentration = &phosphorus,
    }, .{
        .residue_carbon_g_c_per_m2 = &.{ 1, 2 },
        .residue_nitrogen_g_n_per_m2 = &.{ 0.1, 0.2 },
        .residue_phosphorus_g_p_per_m2 = &.{ 0.01, 0.02 },
        .horizontal_area_m2 = 10,
        .dry_layer_mass_Mg = 0,
        .total_layer_volume_m3 = 2,
        .area_scaled_calculation_floor = 1e-14,
        .calculation_floor = 1e-15,
        .layer_kind = .surface_litter,
        .soil_organic_carbon_g_c_per_Mg = 0,
        .particulate_organic_carbon_g_c_per_Mg = 0,
        .soil_organic_nitrogen_g_n_per_Mg = 0,
        .soil_organic_phosphorus_g_p_per_Mg = 0,
        .particulate_nitrogen_to_carbon_g_n_per_g_c = 0.1,
        .particulate_phosphorus_to_carbon_g_p_per_g_c = 0.01,
    });
    try std.testing.expectEqualSlices(f64, &.{ 5, 10, 0, 0 }, &carbon);
}
