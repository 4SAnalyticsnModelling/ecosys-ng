const std = @import("std");
const activity_coefficients = @import("solute_activity_coefficients.zig");
const aqueous_network = @import("solute_aqueous_network.zig");
const phosphate_network = @import("solute_phosphate_network.zig");

pub const FreeActivities = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    carbonate_mol_per_m3: f64,
    bicarbonate_mol_per_m3: f64,
    carbon_dioxide_mol_per_m3: f64,
    ammonium_non_band_mol_per_m3: f64,
    ammonium_band_mol_per_m3: f64,
    ammonia_non_band_mol_per_m3: f64,
    ammonia_band_mol_per_m3: f64,
};

pub const MetalComplexActivities = struct {
    aluminum_hydroxide_1_mol_per_m3: f64,
    aluminum_hydroxide_2_mol_per_m3: f64,
    aluminum_hydroxide_3_mol_per_m3: f64,
    aluminum_hydroxide_4_mol_per_m3: f64,
    aluminum_sulfate_mol_per_m3: f64,
    iron_hydroxide_1_mol_per_m3: f64,
    iron_hydroxide_2_mol_per_m3: f64,
    iron_hydroxide_3_mol_per_m3: f64,
    iron_hydroxide_4_mol_per_m3: f64,
    iron_sulfate_mol_per_m3: f64,
    calcium_hydroxide_mol_per_m3: f64,
    calcium_carbonate_mol_per_m3: f64,
    calcium_bicarbonate_mol_per_m3: f64,
    calcium_sulfate_mol_per_m3: f64,
    magnesium_hydroxide_mol_per_m3: f64,
    magnesium_carbonate_mol_per_m3: f64,
    magnesium_bicarbonate_mol_per_m3: f64,
    magnesium_sulfate_mol_per_m3: f64,
    sodium_carbonate_mol_per_m3: f64,
    sodium_sulfate_mol_per_m3: f64,
    potassium_sulfate_mol_per_m3: f64,
};

pub const PhosphateActivities = struct {
    po4_mol_p_per_m3: f64,
    hpo4_mol_p_per_m3: f64,
    h2po4_mol_p_per_m3: f64,
    h3po4_mol_p_per_m3: f64,
    iron_hpo4_mol_p_per_m3: f64,
    iron_h2po4_mol_p_per_m3: f64,
    calcium_po4_mol_p_per_m3: f64,
    calcium_hpo4_mol_p_per_m3: f64,
    calcium_h2po4_mol_p_per_m3: f64,
    magnesium_hpo4_mol_p_per_m3: f64,
};

pub const Result = struct {
    free: FreeActivities,
    metal_complexes: MetalComplexActivities,
    non_band_phosphate: PhosphateActivities,
    band_phosphate: PhosphateActivities,
};

pub const SurfaceConcentrations = struct {
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    carbon_dioxide_mol_per_m3: f64,
    bicarbonate_mol_per_m3: f64,
    carbonate_mol_per_m3: f64,
    hydrogen_phosphate_mol_p_per_m3: f64,
    dihydrogen_phosphate_mol_p_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
};

pub const SurfaceActivities = SurfaceConcentrations;

pub const SurfaceActivityCoefficients = struct {
    monovalent: f64,
    divalent: f64,
    trivalent: f64,
};

/// Direct SOLUTE lines 866--922 charge-class mapping. Activity coefficients
/// are runtime values; neutral complexes retain their concentrations.
pub fn calculate(
    shared: aqueous_network.State,
    non_band: phosphate_network.State,
    band: phosphate_network.State,
    coefficients: activity_coefficients.Result,
) !Result {
    try validate(shared, non_band, band, coefficients);
    const g1 = coefficients.monovalent_activity_coefficient;
    const g2 = coefficients.divalent_activity_coefficient;
    const g3 = coefficients.trivalent_activity_coefficient;
    const result = Result{
        .free = .{
            .hydrogen_mol_per_m3 = shared.hydrogen * g1,
            .hydroxide_mol_per_m3 = shared.hydroxide * g1,
            .aluminum_mol_per_m3 = shared.aluminum * g3,
            .iron_mol_per_m3 = shared.iron * g3,
            .calcium_mol_per_m3 = shared.calcium * g2,
            .magnesium_mol_per_m3 = shared.magnesium * g2,
            .sodium_mol_per_m3 = shared.sodium * g1,
            .potassium_mol_per_m3 = shared.potassium * g1,
            .sulfate_mol_per_m3 = shared.sulfate * g2,
            .carbonate_mol_per_m3 = shared.carbonate * g2,
            .bicarbonate_mol_per_m3 = shared.bicarbonate * g1,
            .carbon_dioxide_mol_per_m3 = shared.carbon_dioxide,
            .ammonium_non_band_mol_per_m3 = shared.ammonium_non_band * g1,
            .ammonium_band_mol_per_m3 = shared.ammonium_band * g1,
            .ammonia_non_band_mol_per_m3 = shared.ammonia_non_band,
            .ammonia_band_mol_per_m3 = shared.ammonia_band,
        },
        .metal_complexes = .{
            .aluminum_hydroxide_1_mol_per_m3 = shared.aluminum_hydroxide_1 * g2,
            .aluminum_hydroxide_2_mol_per_m3 = shared.aluminum_hydroxide_2 * g1,
            .aluminum_hydroxide_3_mol_per_m3 = shared.aluminum_hydroxide_3,
            .aluminum_hydroxide_4_mol_per_m3 = shared.aluminum_hydroxide_4 * g1,
            .aluminum_sulfate_mol_per_m3 = shared.aluminum_sulfate * g1,
            .iron_hydroxide_1_mol_per_m3 = shared.iron_hydroxide_1 * g2,
            .iron_hydroxide_2_mol_per_m3 = shared.iron_hydroxide_2 * g1,
            .iron_hydroxide_3_mol_per_m3 = shared.iron_hydroxide_3,
            .iron_hydroxide_4_mol_per_m3 = shared.iron_hydroxide_4 * g1,
            .iron_sulfate_mol_per_m3 = shared.iron_sulfate * g1,
            .calcium_hydroxide_mol_per_m3 = shared.calcium_hydroxide * g1,
            .calcium_carbonate_mol_per_m3 = shared.calcium_carbonate,
            .calcium_bicarbonate_mol_per_m3 = shared.calcium_bicarbonate * g1,
            .calcium_sulfate_mol_per_m3 = shared.calcium_sulfate,
            .magnesium_hydroxide_mol_per_m3 = shared.magnesium_hydroxide * g1,
            .magnesium_carbonate_mol_per_m3 = shared.magnesium_carbonate,
            .magnesium_bicarbonate_mol_per_m3 = shared.magnesium_bicarbonate * g1,
            .magnesium_sulfate_mol_per_m3 = shared.magnesium_sulfate,
            .sodium_carbonate_mol_per_m3 = shared.sodium_carbonate * g1,
            .sodium_sulfate_mol_per_m3 = shared.sodium_sulfate * g1,
            .potassium_sulfate_mol_per_m3 = shared.potassium_sulfate * g1,
        },
        .non_band_phosphate = phosphateActivities(non_band, g1, g2, g3),
        .band_phosphate = phosphateActivities(band, g1, g2, g3),
    };
    inline for (@typeInfo(FreeActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.free, field.name)))
            return error.NonFiniteIonActivity;
    inline for (@typeInfo(MetalComplexActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.metal_complexes, field.name)))
            return error.NonFiniteIonActivity;
    inline for (@typeInfo(PhosphateActivities).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.non_band_phosphate, field.name)) or
            !std.math.isFinite(@field(result.band_phosphate, field.name)))
        {
            return error.NonFiniteIonActivity;
        }
    }
    return result;
}

/// Direct source-order translation of SOLUTE.F lines 4274--4289.
///
/// This scalar surface interface uses the same charge-class mapping as the
/// paired-zone owner above without manufacturing unused soil-zone state.
pub fn calculateSurfaceSourceOrder(
    concentrations: SurfaceConcentrations,
    coefficients: SurfaceActivityCoefficients,
) !SurfaceActivities {
    try validateSurface(concentrations, coefficients);
    const result: SurfaceActivities = .{
        .aluminum_mol_per_m3 = concentrations.aluminum_mol_per_m3 * coefficients.trivalent,
        .iron_mol_per_m3 = concentrations.iron_mol_per_m3 * coefficients.trivalent,
        .calcium_mol_per_m3 = concentrations.calcium_mol_per_m3 * coefficients.divalent,
        .magnesium_mol_per_m3 = concentrations.magnesium_mol_per_m3 * coefficients.divalent,
        .sodium_mol_per_m3 = concentrations.sodium_mol_per_m3 * coefficients.monovalent,
        .potassium_mol_per_m3 = concentrations.potassium_mol_per_m3 * coefficients.monovalent,
        .sulfate_mol_per_m3 = concentrations.sulfate_mol_per_m3 * coefficients.divalent,
        .carbon_dioxide_mol_per_m3 = concentrations.carbon_dioxide_mol_per_m3,
        .bicarbonate_mol_per_m3 = concentrations.bicarbonate_mol_per_m3 *
            coefficients.monovalent,
        .carbonate_mol_per_m3 = concentrations.carbonate_mol_per_m3 * coefficients.divalent,
        .hydrogen_phosphate_mol_p_per_m3 = concentrations.hydrogen_phosphate_mol_p_per_m3 *
            coefficients.divalent,
        .dihydrogen_phosphate_mol_p_per_m3 = concentrations.dihydrogen_phosphate_mol_p_per_m3 *
            coefficients.monovalent,
        .ammonium_mol_n_per_m3 = concentrations.ammonium_mol_n_per_m3 *
            coefficients.monovalent,
        .ammonia_mol_n_per_m3 = concentrations.ammonia_mol_n_per_m3,
    };
    inline for (@typeInfo(SurfaceActivities).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteIonActivity;
    }
    return result;
}

fn validateSurface(
    concentrations: SurfaceConcentrations,
    coefficients: SurfaceActivityCoefficients,
) !void {
    inline for (@typeInfo(SurfaceConcentrations).@"struct".fields) |field| {
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidIonActivityState;
    }
    inline for (@typeInfo(SurfaceActivityCoefficients).@"struct".fields) |field| {
        const value = @field(coefficients, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidIonActivityCoefficient;
    }
}

fn phosphateActivities(
    state: phosphate_network.State,
    monovalent: f64,
    divalent: f64,
    trivalent: f64,
) PhosphateActivities {
    return .{
        .po4_mol_p_per_m3 = state.dissolved_po4_mol_p_per_m3 * trivalent,
        .hpo4_mol_p_per_m3 = state.dissolved_hpo4_mol_p_per_m3 * divalent,
        .h2po4_mol_p_per_m3 = state.dissolved_h2po4_mol_p_per_m3 * monovalent,
        .h3po4_mol_p_per_m3 = state.dissolved_h3po4_mol_p_per_m3,
        .iron_hpo4_mol_p_per_m3 = state.iron_hpo4_pair_mol_per_m3 * divalent,
        .iron_h2po4_mol_p_per_m3 = state.iron_h2po4_pair_mol_per_m3 * monovalent,
        .calcium_po4_mol_p_per_m3 = state.calcium_po4_pair_mol_per_m3 * monovalent,
        .calcium_hpo4_mol_p_per_m3 = state.calcium_hpo4_pair_mol_per_m3,
        .calcium_h2po4_mol_p_per_m3 = state.calcium_h2po4_pair_mol_per_m3 * monovalent,
        .magnesium_hpo4_mol_p_per_m3 = state.magnesium_hpo4_pair_mol_per_m3,
    };
}

fn validate(
    shared: aqueous_network.State,
    non_band: phosphate_network.State,
    band: phosphate_network.State,
    coefficients: activity_coefficients.Result,
) !void {
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        const value = @field(shared, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidIonActivityState;
    }
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
        const non_band_value = @field(non_band, field.name);
        const band_value = @field(band, field.name);
        if (!std.math.isFinite(non_band_value) or non_band_value < 0 or
            !std.math.isFinite(band_value) or band_value < 0)
        {
            return error.InvalidIonActivityState;
        }
    }
    inline for (.{
        coefficients.monovalent_activity_coefficient,
        coefficients.divalent_activity_coefficient,
        coefficients.trivalent_activity_coefficient,
    }) |coefficient| {
        if (!std.math.isFinite(coefficient) or coefficient <= 0)
            return error.InvalidIonActivityCoefficient;
    }
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "SOLUTE ion activities preserve every source charge class" {
    const result = try calculate(
        filled(aqueous_network.State, 8),
        filled(phosphate_network.State, 8),
        filled(phosphate_network.State, 4),
        .{
            .ionic_strength_mol_per_l = 0,
            .monovalent_activity_coefficient = 0.5,
            .divalent_activity_coefficient = 0.25,
            .trivalent_activity_coefficient = 0.125,
            .total_ion_activity_mol_per_m3 = 0,
            .electrical_conductivity_dS_per_m = 0,
        },
    );
    inline for (.{
        "hydrogen_mol_per_m3",
        "hydroxide_mol_per_m3",
        "sodium_mol_per_m3",
        "potassium_mol_per_m3",
        "bicarbonate_mol_per_m3",
        "ammonium_non_band_mol_per_m3",
        "ammonium_band_mol_per_m3",
    }) |field_name|
        try std.testing.expectEqual(@as(f64, 4), @field(result.free, field_name));
    inline for (.{
        "calcium_mol_per_m3",
        "magnesium_mol_per_m3",
        "sulfate_mol_per_m3",
        "carbonate_mol_per_m3",
    }) |field_name|
        try std.testing.expectEqual(@as(f64, 2), @field(result.free, field_name));
    try std.testing.expectEqual(@as(f64, 1), result.free.aluminum_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), result.free.iron_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 8), result.free.carbon_dioxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 8), result.free.ammonia_non_band_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 2), result.metal_complexes.aluminum_hydroxide_1_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), result.metal_complexes.aluminum_hydroxide_2_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 8), result.metal_complexes.aluminum_hydroxide_3_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 8), result.metal_complexes.calcium_carbonate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), result.metal_complexes.potassium_sulfate_mol_per_m3);
}

test "SOLUTE phosphate activities distinguish both runtime zones" {
    const result = try calculate(
        filled(aqueous_network.State, 1),
        filled(phosphate_network.State, 8),
        filled(phosphate_network.State, 4),
        .{
            .ionic_strength_mol_per_l = 0,
            .monovalent_activity_coefficient = 0.5,
            .divalent_activity_coefficient = 0.25,
            .trivalent_activity_coefficient = 0.125,
            .total_ion_activity_mol_per_m3 = 0,
            .electrical_conductivity_dS_per_m = 0,
        },
    );
    try std.testing.expectEqual(
        @as(f64, 1),
        result.non_band_phosphate.po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        result.non_band_phosphate.hpo4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 4),
        result.non_band_phosphate.h2po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 8),
        result.non_band_phosphate.calcium_hpo4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.band_phosphate.po4_mol_p_per_m3,
    );
}

fn testSurfaceConcentrations() SurfaceConcentrations {
    return .{
        .aluminum_mol_per_m3 = 1,
        .iron_mol_per_m3 = 2,
        .calcium_mol_per_m3 = 3,
        .magnesium_mol_per_m3 = 4,
        .sodium_mol_per_m3 = 5,
        .potassium_mol_per_m3 = 6,
        .sulfate_mol_per_m3 = 7,
        .carbon_dioxide_mol_per_m3 = 8,
        .bicarbonate_mol_per_m3 = 9,
        .carbonate_mol_per_m3 = 10,
        .hydrogen_phosphate_mol_p_per_m3 = 11,
        .dihydrogen_phosphate_mol_p_per_m3 = 12,
        .ammonium_mol_n_per_m3 = 13,
        .ammonia_mol_n_per_m3 = 14,
    };
}

test "SOLUTE initial surface activities preserve every source expression" {
    const concentrations = testSurfaceConcentrations();
    const coefficients: SurfaceActivityCoefficients = .{
        .monovalent = 0.5,
        .divalent = 0.25,
        .trivalent = 0.125,
    };
    const result = try calculateSurfaceSourceOrder(
        concentrations,
        coefficients,
    );

    try std.testing.expectEqual(
        concentrations.aluminum_mol_per_m3 * coefficients.trivalent,
        result.aluminum_mol_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.iron_mol_per_m3 * coefficients.trivalent,
        result.iron_mol_per_m3,
    );
    inline for (.{
        "calcium_mol_per_m3",
        "magnesium_mol_per_m3",
        "sulfate_mol_per_m3",
        "carbonate_mol_per_m3",
        "hydrogen_phosphate_mol_p_per_m3",
    }) |field_name| {
        try std.testing.expectEqual(
            @field(concentrations, field_name) * coefficients.divalent,
            @field(result, field_name),
        );
    }
    inline for (.{
        "sodium_mol_per_m3",
        "potassium_mol_per_m3",
        "bicarbonate_mol_per_m3",
        "dihydrogen_phosphate_mol_p_per_m3",
        "ammonium_mol_n_per_m3",
    }) |field_name| {
        try std.testing.expectEqual(
            @field(concentrations, field_name) * coefficients.monovalent,
            @field(result, field_name),
        );
    }
    try std.testing.expectEqual(
        concentrations.carbon_dioxide_mol_per_m3,
        result.carbon_dioxide_mol_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.ammonia_mol_n_per_m3,
        result.ammonia_mol_n_per_m3,
    );
}

test "neutral surface species remain unchanged" {
    var concentrations = std.mem.zeroes(SurfaceConcentrations);
    concentrations.carbon_dioxide_mol_per_m3 = 2;
    concentrations.ammonia_mol_n_per_m3 = 3;
    const result = try calculateSurfaceSourceOrder(
        concentrations,
        .{ .monovalent = 0.7, .divalent = 0.4, .trivalent = 0.2 },
    );
    try std.testing.expectEqual(
        concentrations.carbon_dioxide_mol_per_m3,
        result.carbon_dioxide_mol_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.ammonia_mol_n_per_m3,
        result.ammonia_mol_n_per_m3,
    );
}

test "surface activity charge classes use independent runtime coefficients" {
    const concentrations = filled(SurfaceConcentrations, 8);
    const result = try calculateSurfaceSourceOrder(
        concentrations,
        .{ .monovalent = 0.5, .divalent = 0.25, .trivalent = 0.125 },
    );
    try std.testing.expectEqual(@as(f64, 1), result.aluminum_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 2), result.calcium_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), result.sodium_mol_per_m3);
    try std.testing.expectEqual(
        @as(f64, 8),
        result.carbon_dioxide_mol_per_m3,
    );
}

test "surface activities reject invalid input and overflow" {
    var concentrations = testSurfaceConcentrations();
    concentrations.calcium_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidIonActivityState,
        calculateSurfaceSourceOrder(
            concentrations,
            .{ .monovalent = 1, .divalent = 1, .trivalent = 1 },
        ),
    );

    concentrations = testSurfaceConcentrations();
    try std.testing.expectError(
        error.InvalidIonActivityCoefficient,
        calculateSurfaceSourceOrder(
            concentrations,
            .{ .monovalent = 0, .divalent = 1, .trivalent = 1 },
        ),
    );

    concentrations.aluminum_mol_per_m3 = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteIonActivity,
        calculateSurfaceSourceOrder(
            concentrations,
            .{ .monovalent = 1, .divalent = 1, .trivalent = 2 },
        ),
    );
}
