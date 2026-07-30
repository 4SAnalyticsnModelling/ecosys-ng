const std = @import("std");

pub const ratio_class_count = 3;

pub const Parameters = struct {
    substrate_class_count: usize,
    population_count: usize,
    heterotrophic_substrate_class_count: usize,
    fungal_population_index: usize,
    labile_fraction: f64,
    recalcitrant_fraction: f64,
    fungal_labile_nitrogen_to_carbon_g_n_per_g_c: f64,
    fungal_recalcitrant_nitrogen_to_carbon_g_n_per_g_c: f64,
    fungal_labile_phosphorus_to_carbon_g_p_per_g_c: f64,
    fungal_recalcitrant_phosphorus_to_carbon_g_p_per_g_c: f64,
    other_labile_nitrogen_to_carbon_g_n_per_g_c: f64,
    other_recalcitrant_nitrogen_to_carbon_g_n_per_g_c: f64,
    other_labile_phosphorus_to_carbon_g_p_per_g_c: f64,
    other_recalcitrant_phosphorus_to_carbon_g_p_per_g_c: f64,
};

pub const State = struct {
    /// `[substrate][population][labile,recalcitrant,mixture]`.
    nitrogen_to_carbon_g_n_per_g_c: []f64,
    /// `[substrate][population][labile,recalcitrant,mixture]`.
    phosphorus_to_carbon_g_p_per_g_c: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 220--249.
pub fn initialize(state: State, parameters: Parameters) !void {
    if (parameters.substrate_class_count == 0 or
        parameters.population_count == 0 or
        parameters.heterotrophic_substrate_class_count >
            parameters.substrate_class_count or
        parameters.fungal_population_index >= parameters.population_count)
    {
        return error.InvalidMicrobialNutrientDimensions;
    }
    const population_substrate_count = std.math.mul(
        usize,
        parameters.substrate_class_count,
        parameters.population_count,
    ) catch return error.DimensionOverflow;
    const value_count = std.math.mul(
        usize,
        population_substrate_count,
        ratio_class_count,
    ) catch return error.DimensionOverflow;
    if (state.nitrogen_to_carbon_g_n_per_g_c.len != value_count or
        state.phosphorus_to_carbon_g_p_per_g_c.len != value_count)
    {
        return error.MicrobialNutrientDimensionMismatch;
    }

    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(parameters, field.name);
            if (!std.math.isFinite(value))
                return error.NonFiniteMicrobialNutrientParameter;
            if (value < 0) return error.InvalidMicrobialNutrientParameter;
        }
    }
    if (@abs(parameters.labile_fraction +
        parameters.recalcitrant_fraction - 1.0) > 1.0e-12)
    {
        return error.InvalidMicrobialPartition;
    }

    for (0..parameters.substrate_class_count) |substrate_index| {
        for (0..parameters.population_count) |population_index| {
            const fungal = substrate_index <
                parameters.heterotrophic_substrate_class_count and
                population_index == parameters.fungal_population_index;
            const nitrogen_labile =
                if (fungal)
                    parameters.fungal_labile_nitrogen_to_carbon_g_n_per_g_c
                else
                    parameters.other_labile_nitrogen_to_carbon_g_n_per_g_c;
            const nitrogen_recalcitrant =
                if (fungal)
                    parameters.fungal_recalcitrant_nitrogen_to_carbon_g_n_per_g_c
                else
                    parameters.other_recalcitrant_nitrogen_to_carbon_g_n_per_g_c;
            const phosphorus_labile =
                if (fungal)
                    parameters.fungal_labile_phosphorus_to_carbon_g_p_per_g_c
                else
                    parameters.other_labile_phosphorus_to_carbon_g_p_per_g_c;
            const phosphorus_recalcitrant =
                if (fungal)
                    parameters.fungal_recalcitrant_phosphorus_to_carbon_g_p_per_g_c
                else
                    parameters.other_recalcitrant_phosphorus_to_carbon_g_p_per_g_c;

            const base =
                (substrate_index * parameters.population_count +
                    population_index) * ratio_class_count;
            state.nitrogen_to_carbon_g_n_per_g_c[base] = nitrogen_labile;
            state.nitrogen_to_carbon_g_n_per_g_c[base + 1] =
                nitrogen_recalcitrant;
            state.phosphorus_to_carbon_g_p_per_g_c[base] =
                phosphorus_labile;
            state.phosphorus_to_carbon_g_p_per_g_c[base + 1] =
                phosphorus_recalcitrant;
            state.nitrogen_to_carbon_g_n_per_g_c[base + 2] =
                parameters.labile_fraction * nitrogen_labile +
                parameters.recalcitrant_fraction * nitrogen_recalcitrant;
            state.phosphorus_to_carbon_g_p_per_g_c[base + 2] =
                parameters.labile_fraction * phosphorus_labile +
                parameters.recalcitrant_fraction * phosphorus_recalcitrant;
        }
    }
}

fn sourceParameters() Parameters {
    return .{
        .substrate_class_count = 6,
        .population_count = 7,
        .heterotrophic_substrate_class_count = 5,
        .fungal_population_index = 2,
        .labile_fraction = 0.55,
        .recalcitrant_fraction = 0.45,
        .fungal_labile_nitrogen_to_carbon_g_n_per_g_c = 0.111,
        .fungal_recalcitrant_nitrogen_to_carbon_g_n_per_g_c = 0.083,
        .fungal_labile_phosphorus_to_carbon_g_p_per_g_c = 0.0111,
        .fungal_recalcitrant_phosphorus_to_carbon_g_p_per_g_c = 0.0083,
        .other_labile_nitrogen_to_carbon_g_n_per_g_c = 0.167,
        .other_recalcitrant_nitrogen_to_carbon_g_n_per_g_c = 0.125,
        .other_labile_phosphorus_to_carbon_g_p_per_g_c = 0.0167,
        .other_recalcitrant_phosphorus_to_carbon_g_p_per_g_c = 0.0125,
    };
}

test "STARTS microbial nutrient tables reproduce fungal and other branches" {
    var nitrogen = [_]f64{-1.0} ** (6 * 7 * ratio_class_count);
    var phosphorus = [_]f64{-1.0} ** (6 * 7 * ratio_class_count);
    const parameters = sourceParameters();
    try initialize(.{
        .nitrogen_to_carbon_g_n_per_g_c = &nitrogen,
        .phosphorus_to_carbon_g_p_per_g_c = &phosphorus,
    }, parameters);

    const fungal_base = (4 * 7 + 2) * ratio_class_count;
    try std.testing.expectEqual(@as(f64, 0.111), nitrogen[fungal_base]);
    try std.testing.expectEqual(@as(f64, 0.083), nitrogen[fungal_base + 1]);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.55 * 0.111 + 0.45 * 0.083),
        nitrogen[fungal_base + 2],
        1.0e-15,
    );
    const autotroph_fungus_slot = (5 * 7 + 2) * ratio_class_count;
    try std.testing.expectEqual(
        @as(f64, 0.167),
        nitrogen[autotroph_fungus_slot],
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.55 * 0.0167 + 0.45 * 0.0125),
        phosphorus[autotroph_fungus_slot + 2],
        1.0e-15,
    );
}

test "runtime population and substrate extents determine storage" {
    var nitrogen = [_]f64{0.0} ** (3 * 4 * ratio_class_count);
    var phosphorus = [_]f64{0.0} ** (3 * 4 * ratio_class_count);
    var parameters = sourceParameters();
    parameters.substrate_class_count = 3;
    parameters.population_count = 4;
    parameters.heterotrophic_substrate_class_count = 2;
    parameters.fungal_population_index = 3;

    try initialize(.{
        .nitrogen_to_carbon_g_n_per_g_c = &nitrogen,
        .phosphorus_to_carbon_g_p_per_g_c = &phosphorus,
    }, parameters);
    const fungal_base = (1 * 4 + 3) * ratio_class_count;
    try std.testing.expectEqual(@as(f64, 0.111), nitrogen[fungal_base]);
    const non_fungal_base = (2 * 4 + 3) * ratio_class_count;
    try std.testing.expectEqual(@as(f64, 0.167), nitrogen[non_fungal_base]);
}

test "invalid partition fails without modifying either table" {
    var nitrogen = [_]f64{ 7.0, 8.0, 9.0 };
    var phosphorus = [_]f64{ 10.0, 11.0, 12.0 };
    const nitrogen_before = nitrogen;
    const phosphorus_before = phosphorus;
    var parameters = sourceParameters();
    parameters.substrate_class_count = 1;
    parameters.population_count = 1;
    parameters.heterotrophic_substrate_class_count = 1;
    parameters.fungal_population_index = 0;
    parameters.recalcitrant_fraction = 0.40;

    try std.testing.expectError(
        error.InvalidMicrobialPartition,
        initialize(.{
            .nitrogen_to_carbon_g_n_per_g_c = &nitrogen,
            .phosphorus_to_carbon_g_p_per_g_c = &phosphorus,
        }, parameters),
    );
    try std.testing.expectEqualSlices(f64, &nitrogen_before, &nitrogen);
    try std.testing.expectEqualSlices(f64, &phosphorus_before, &phosphorus);
}
