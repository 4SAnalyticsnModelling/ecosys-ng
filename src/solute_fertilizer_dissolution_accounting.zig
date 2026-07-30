const std = @import("std");
const fertilizer_dissolution = @import("soil_fertilizer_dissolution.zig");

pub const MolarMasses = struct {
    carbon_g_per_mol: f64,
    nitrogen_g_per_mol: f64,
    phosphorus_g_per_mol: f64,
};

pub const MolarTransformationTotals = struct {
    carbon_dioxide_mol_c_per_step: f64,
    gaseous_ammonia_mol_n_per_step: f64,
    non_band_ammonium_mol_n_per_step: f64,
    band_ammonium_mol_n_per_step: f64,
    non_band_ammonia_mol_n_per_step: f64,
    band_ammonia_mol_n_per_step: f64,
    non_band_nitrate_mol_n_per_step: f64,
    band_nitrate_mol_n_per_step: f64,
    non_band_nitrite_mol_n_per_step: f64,
    band_nitrite_mol_n_per_step: f64,
    non_band_hpo4_mol_p_per_step: f64,
    non_band_h2po4_mol_p_per_step: f64,
    band_hpo4_mol_p_per_step: f64,
    band_h2po4_mol_p_per_step: f64,
};

pub const MassTransformationTotals = struct {
    carbon_dioxide_g_c_per_step: f64,
    gaseous_ammonia_g_n_per_step: f64,
    non_band_ammonium_g_n_per_step: f64,
    band_ammonium_g_n_per_step: f64,
    non_band_ammonia_g_n_per_step: f64,
    band_ammonia_g_n_per_step: f64,
    non_band_nitrate_g_n_per_step: f64,
    band_nitrate_g_n_per_step: f64,
    non_band_nitrite_g_n_per_step: f64,
    band_nitrite_g_n_per_step: f64,
    non_band_hpo4_g_p_per_step: f64,
    non_band_h2po4_g_p_per_step: f64,
    band_hpo4_g_p_per_step: f64,
    band_h2po4_g_p_per_step: f64,
};

pub const Inputs = struct {
    fertilizer: fertilizer_dissolution.FertilizerState,
    dissolution: fertilizer_dissolution.DissolutionFlux,
    preceding_molar_transformations: MolarTransformationTotals,
    molar_masses: MolarMasses,
};

pub const Result = struct {
    remaining_fertilizer: fertilizer_dissolution.FertilizerState,
    molar_transformations_before_scaling: MolarTransformationTotals,
    mass_transformations: MassTransformationTotals,
};

/// Exact source-order comparator for SOLUTE.F lines 3917--3989.
///
/// The source first subtracts every dissolution flux from its undissolved
/// donor, then adds NH4 fertilizer to aqueous NH4, NH3 fertilizer to gaseous
/// NH3, urea hydrolysis to aqueous NH3, and nitrate fertilizer to aqueous
/// NO3. Only after those additions are all transformation totals converted
/// from moles to grams.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const flux = inputs.dissolution;
    var remaining = inputs.fertilizer;

    // SOLUTE.F 3940--3947. Preserve left-associated subtraction order.
    remaining.broadcast_ammonium_mol_n =
        remaining.broadcast_ammonium_mol_n -
        flux.broadcast_ammonium_non_band_mol_n -
        flux.broadcast_ammonium_band_mol_n;
    remaining.broadcast_ammonia_mol_n =
        remaining.broadcast_ammonia_mol_n -
        flux.broadcast_ammonia_non_band_mol_n -
        flux.broadcast_ammonia_band_mol_n;
    remaining.broadcast_urea_mol_n =
        remaining.broadcast_urea_mol_n -
        flux.broadcast_urea_non_band_mol_n -
        flux.broadcast_urea_band_mol_n;
    remaining.broadcast_nitrate_mol_n =
        remaining.broadcast_nitrate_mol_n -
        flux.broadcast_nitrate_non_band_mol_n -
        flux.broadcast_nitrate_band_mol_n;
    remaining.banded_ammonium_mol_n =
        remaining.banded_ammonium_mol_n -
        flux.banded_ammonium_mol_n;
    remaining.banded_ammonia_mol_n =
        remaining.banded_ammonia_mol_n -
        flux.banded_ammonia_mol_n;
    remaining.banded_urea_mol_n =
        remaining.banded_urea_mol_n -
        flux.banded_urea_mol_n;
    remaining.banded_nitrate_mol_n =
        remaining.banded_nitrate_mol_n -
        flux.banded_nitrate_mol_n;

    // SOLUTE.F 3959--3965. These recipient identities are source-specific.
    var molar = inputs.preceding_molar_transformations;
    molar.gaseous_ammonia_mol_n_per_step =
        molar.gaseous_ammonia_mol_n_per_step +
        flux.broadcast_ammonia_non_band_mol_n +
        flux.broadcast_ammonia_band_mol_n +
        flux.banded_ammonia_mol_n;
    molar.non_band_ammonium_mol_n_per_step =
        molar.non_band_ammonium_mol_n_per_step +
        flux.broadcast_ammonium_non_band_mol_n;
    molar.band_ammonium_mol_n_per_step =
        molar.band_ammonium_mol_n_per_step +
        flux.broadcast_ammonium_band_mol_n +
        flux.banded_ammonium_mol_n;
    molar.non_band_ammonia_mol_n_per_step =
        molar.non_band_ammonia_mol_n_per_step +
        flux.broadcast_urea_non_band_mol_n;
    molar.band_ammonia_mol_n_per_step =
        molar.band_ammonia_mol_n_per_step +
        flux.broadcast_urea_band_mol_n +
        flux.banded_urea_mol_n;
    molar.non_band_nitrate_mol_n_per_step =
        molar.non_band_nitrate_mol_n_per_step +
        flux.broadcast_nitrate_non_band_mol_n;
    molar.band_nitrate_mol_n_per_step =
        molar.band_nitrate_mol_n_per_step +
        flux.broadcast_nitrate_band_mol_n +
        flux.banded_nitrate_mol_n;

    const masses = inputs.molar_masses;
    // SOLUTE.F 3976--3989. Runtime masses replace PARAMETER constants.
    const mass: MassTransformationTotals = .{
        .carbon_dioxide_g_c_per_step = molar.carbon_dioxide_mol_c_per_step * masses.carbon_g_per_mol,
        .gaseous_ammonia_g_n_per_step = molar.gaseous_ammonia_mol_n_per_step * masses.nitrogen_g_per_mol,
        .non_band_ammonium_g_n_per_step = molar.non_band_ammonium_mol_n_per_step * masses.nitrogen_g_per_mol,
        .band_ammonium_g_n_per_step = molar.band_ammonium_mol_n_per_step * masses.nitrogen_g_per_mol,
        .non_band_ammonia_g_n_per_step = molar.non_band_ammonia_mol_n_per_step * masses.nitrogen_g_per_mol,
        .band_ammonia_g_n_per_step = molar.band_ammonia_mol_n_per_step * masses.nitrogen_g_per_mol,
        .non_band_nitrate_g_n_per_step = molar.non_band_nitrate_mol_n_per_step * masses.nitrogen_g_per_mol,
        .band_nitrate_g_n_per_step = molar.band_nitrate_mol_n_per_step * masses.nitrogen_g_per_mol,
        .non_band_nitrite_g_n_per_step = molar.non_band_nitrite_mol_n_per_step * masses.nitrogen_g_per_mol,
        .band_nitrite_g_n_per_step = molar.band_nitrite_mol_n_per_step * masses.nitrogen_g_per_mol,
        .non_band_hpo4_g_p_per_step = molar.non_band_hpo4_mol_p_per_step * masses.phosphorus_g_per_mol,
        .non_band_h2po4_g_p_per_step = molar.non_band_h2po4_mol_p_per_step * masses.phosphorus_g_per_mol,
        .band_hpo4_g_p_per_step = molar.band_hpo4_mol_p_per_step * masses.phosphorus_g_per_mol,
        .band_h2po4_g_p_per_step = molar.band_h2po4_mol_p_per_step * masses.phosphorus_g_per_mol,
    };
    try validateResult(remaining, molar, mass);
    return .{
        .remaining_fertilizer = remaining,
        .molar_transformations_before_scaling = molar,
        .mass_transformations = mass,
    };
}

fn validateInputs(inputs: Inputs) !void {
    try validateFiniteStruct(
        inputs.fertilizer,
        error.NonFiniteFertilizerAccountingInput,
    );
    inline for (@typeInfo(
        fertilizer_dissolution.FertilizerState,
    ).@"struct".fields) |field| {
        if (@field(inputs.fertilizer, field.name) < 0)
            return error.InvalidFertilizerAccountingInput;
    }
    try validateFiniteStruct(
        inputs.dissolution,
        error.NonFiniteFertilizerAccountingInput,
    );
    inline for (@typeInfo(
        fertilizer_dissolution.DissolutionFlux,
    ).@"struct".fields) |field| {
        if (@field(inputs.dissolution, field.name) < 0)
            return error.InvalidFertilizerAccountingInput;
    }
    try validateFiniteStruct(
        inputs.preceding_molar_transformations,
        error.NonFiniteFertilizerAccountingInput,
    );
    try validateFiniteStruct(
        inputs.molar_masses,
        error.NonFiniteFertilizerAccountingInput,
    );
    inline for (@typeInfo(MolarMasses).@"struct".fields) |field| {
        if (@field(inputs.molar_masses, field.name) <= 0)
            return error.InvalidFertilizerAccountingInput;
    }
}

fn validateResult(
    remaining: fertilizer_dissolution.FertilizerState,
    molar: MolarTransformationTotals,
    mass: MassTransformationTotals,
) !void {
    try validateFiniteStruct(
        remaining,
        error.NonFiniteFertilizerAccountingResult,
    );
    inline for (@typeInfo(
        fertilizer_dissolution.FertilizerState,
    ).@"struct".fields) |field| {
        if (@field(remaining, field.name) < 0)
            return error.FertilizerDissolutionExceedsAccountingPool;
    }
    try validateFiniteStruct(
        molar,
        error.NonFiniteFertilizerAccountingResult,
    );
    try validateFiniteStruct(
        mass,
        error.NonFiniteFertilizerAccountingResult,
    );
}

fn validateFiniteStruct(values: anytype, failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(values)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(values, field.name))) return failure;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn testFertilizer() fertilizer_dissolution.FertilizerState {
    return .{
        .broadcast_ammonium_mol_n = 10,
        .broadcast_ammonia_mol_n = 11,
        .broadcast_urea_mol_n = 12,
        .broadcast_nitrate_mol_n = 13,
        .banded_ammonium_mol_n = 14,
        .banded_ammonia_mol_n = 15,
        .banded_urea_mol_n = 16,
        .banded_nitrate_mol_n = 17,
    };
}

fn testFlux() fertilizer_dissolution.DissolutionFlux {
    return .{
        .broadcast_ammonium_non_band_mol_n = 0.1,
        .broadcast_ammonium_band_mol_n = 0.2,
        .broadcast_ammonia_non_band_mol_n = 0.3,
        .broadcast_ammonia_band_mol_n = 0.4,
        .broadcast_urea_non_band_mol_n = 0.5,
        .broadcast_urea_band_mol_n = 0.6,
        .broadcast_nitrate_non_band_mol_n = 0.7,
        .broadcast_nitrate_band_mol_n = 0.8,
        .banded_ammonium_mol_n = 0.9,
        .banded_ammonia_mol_n = 1.0,
        .banded_urea_mol_n = 1.1,
        .banded_nitrate_mol_n = 1.2,
    };
}

test "fertilizer dissolution accounting preserves every source statement" {
    const inputs: Inputs = .{
        .fertilizer = testFertilizer(),
        .dissolution = testFlux(),
        .preceding_molar_transformations = filled(MolarTransformationTotals, 0.25),
        .molar_masses = .{
            .carbon_g_per_mol = 12.01,
            .nitrogen_g_per_mol = 14.01,
            .phosphorus_g_per_mol = 30.97,
        },
    };
    const result = try calculateSourceOrder(inputs);
    const flux = inputs.dissolution;
    const remaining = result.remaining_fertilizer;
    try std.testing.expectEqual(
        inputs.fertilizer.broadcast_ammonium_mol_n -
            flux.broadcast_ammonium_non_band_mol_n -
            flux.broadcast_ammonium_band_mol_n,
        remaining.broadcast_ammonium_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.broadcast_ammonia_mol_n -
            flux.broadcast_ammonia_non_band_mol_n -
            flux.broadcast_ammonia_band_mol_n,
        remaining.broadcast_ammonia_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.broadcast_urea_mol_n -
            flux.broadcast_urea_non_band_mol_n -
            flux.broadcast_urea_band_mol_n,
        remaining.broadcast_urea_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.broadcast_nitrate_mol_n -
            flux.broadcast_nitrate_non_band_mol_n -
            flux.broadcast_nitrate_band_mol_n,
        remaining.broadcast_nitrate_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.banded_ammonium_mol_n -
            flux.banded_ammonium_mol_n,
        remaining.banded_ammonium_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.banded_ammonia_mol_n -
            flux.banded_ammonia_mol_n,
        remaining.banded_ammonia_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.banded_urea_mol_n -
            flux.banded_urea_mol_n,
        remaining.banded_urea_mol_n,
    );
    try std.testing.expectEqual(
        inputs.fertilizer.banded_nitrate_mol_n -
            flux.banded_nitrate_mol_n,
        remaining.banded_nitrate_mol_n,
    );

    const molar = result.molar_transformations_before_scaling;
    try std.testing.expectEqual(
        0.25 +
            flux.broadcast_ammonia_non_band_mol_n +
            flux.broadcast_ammonia_band_mol_n +
            flux.banded_ammonia_mol_n,
        molar.gaseous_ammonia_mol_n_per_step,
    );
    try std.testing.expectEqual(
        0.25 + flux.broadcast_ammonium_non_band_mol_n,
        molar.non_band_ammonium_mol_n_per_step,
    );
    try std.testing.expectEqual(
        0.25 +
            flux.broadcast_ammonium_band_mol_n +
            flux.banded_ammonium_mol_n,
        molar.band_ammonium_mol_n_per_step,
    );
    try std.testing.expectEqual(
        0.25 + flux.broadcast_urea_non_band_mol_n,
        molar.non_band_ammonia_mol_n_per_step,
    );
    try std.testing.expectEqual(
        0.25 +
            flux.broadcast_urea_band_mol_n +
            flux.banded_urea_mol_n,
        molar.band_ammonia_mol_n_per_step,
    );
    try std.testing.expectEqual(
        0.25 + flux.broadcast_nitrate_non_band_mol_n,
        molar.non_band_nitrate_mol_n_per_step,
    );
    try std.testing.expectEqual(
        0.25 +
            flux.broadcast_nitrate_band_mol_n +
            flux.banded_nitrate_mol_n,
        molar.band_nitrate_mol_n_per_step,
    );

    const mass = result.mass_transformations;
    inline for (@typeInfo(MassTransformationTotals).@"struct".fields, 0..) |
        field,
        index,
    | {
        const molar_field =
            @typeInfo(MolarTransformationTotals).@"struct".fields[index];
        const factor = if (index == 0)
            inputs.molar_masses.carbon_g_per_mol
        else if (index < 10)
            inputs.molar_masses.nitrogen_g_per_mol
        else
            inputs.molar_masses.phosphorus_g_per_mol;
        try std.testing.expectEqual(
            @field(molar, molar_field.name) * factor,
            @field(mass, field.name),
        );
    }
}

test "fertilizer dissolution source recipients conserve nitrogen" {
    const fertilizer = testFertilizer();
    const flux = testFlux();
    const result = try calculateSourceOrder(.{
        .fertilizer = fertilizer,
        .dissolution = flux,
        .preceding_molar_transformations = filled(MolarTransformationTotals, 0),
        .molar_masses = .{
            .carbon_g_per_mol = 12,
            .nitrogen_g_per_mol = 14,
            .phosphorus_g_per_mol = 31,
        },
    });
    const before = sumFertilizer(fertilizer);
    const remaining = sumFertilizer(result.remaining_fertilizer);
    const added = result.molar_transformations_before_scaling;
    const recipient_nitrogen =
        added.gaseous_ammonia_mol_n_per_step +
        added.non_band_ammonium_mol_n_per_step +
        added.band_ammonium_mol_n_per_step +
        added.non_band_ammonia_mol_n_per_step +
        added.band_ammonia_mol_n_per_step +
        added.non_band_nitrate_mol_n_per_step +
        added.band_nitrate_mol_n_per_step;
    try std.testing.expectApproxEqAbs(
        before - remaining,
        recipient_nitrogen,
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        recipient_nitrogen * 14,
        result.mass_transformations.gaseous_ammonia_g_n_per_step +
            result.mass_transformations.non_band_ammonium_g_n_per_step +
            result.mass_transformations.band_ammonium_g_n_per_step +
            result.mass_transformations.non_band_ammonia_g_n_per_step +
            result.mass_transformations.band_ammonia_g_n_per_step +
            result.mass_transformations.non_band_nitrate_g_n_per_step +
            result.mass_transformations.band_nitrate_g_n_per_step,
        3.0e-14,
    );
}

test "source routes fertilizer ammonia to gas and urea to aqueous ammonia" {
    var fertilizer = std.mem.zeroes(
        fertilizer_dissolution.FertilizerState,
    );
    fertilizer.broadcast_ammonia_mol_n = 1;
    fertilizer.banded_ammonia_mol_n = 1;
    fertilizer.broadcast_urea_mol_n = 1;
    fertilizer.banded_urea_mol_n = 1;
    var flux = std.mem.zeroes(fertilizer_dissolution.DissolutionFlux);
    flux.broadcast_ammonia_non_band_mol_n = 0.2;
    flux.broadcast_ammonia_band_mol_n = 0.3;
    flux.banded_ammonia_mol_n = 0.4;
    flux.broadcast_urea_non_band_mol_n = 0.1;
    flux.broadcast_urea_band_mol_n = 0.15;
    flux.banded_urea_mol_n = 0.25;
    const result = try calculateSourceOrder(.{
        .fertilizer = fertilizer,
        .dissolution = flux,
        .preceding_molar_transformations = filled(MolarTransformationTotals, 0),
        .molar_masses = .{
            .carbon_g_per_mol = 12,
            .nitrogen_g_per_mol = 14,
            .phosphorus_g_per_mol = 31,
        },
    });
    const molar = result.molar_transformations_before_scaling;
    try std.testing.expectEqual(
        @as(f64, 0.9),
        molar.gaseous_ammonia_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.1),
        molar.non_band_ammonia_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.4),
        molar.band_ammonia_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        molar.non_band_ammonium_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        molar.band_ammonium_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 12.6),
        result.mass_transformations.gaseous_ammonia_g_n_per_step,
    );
}

test "fertilizer accounting rejects invalid donors factors and overflow" {
    var inputs: Inputs = .{
        .fertilizer = testFertilizer(),
        .dissolution = testFlux(),
        .preceding_molar_transformations = filled(MolarTransformationTotals, 0),
        .molar_masses = .{
            .carbon_g_per_mol = 12,
            .nitrogen_g_per_mol = 14,
            .phosphorus_g_per_mol = 31,
        },
    };
    inputs.dissolution.broadcast_ammonium_non_band_mol_n = 20;
    try std.testing.expectError(
        error.FertilizerDissolutionExceedsAccountingPool,
        calculateSourceOrder(inputs),
    );
    inputs.dissolution = testFlux();
    inputs.molar_masses.nitrogen_g_per_mol = 0;
    try std.testing.expectError(
        error.InvalidFertilizerAccountingInput,
        calculateSourceOrder(inputs),
    );
    inputs.molar_masses.nitrogen_g_per_mol = 14;
    inputs.fertilizer.banded_nitrate_mol_n = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteFertilizerAccountingInput,
        calculateSourceOrder(inputs),
    );
    inputs.fertilizer = testFertilizer();
    inputs.preceding_molar_transformations.carbon_dioxide_mol_c_per_step =
        std.math.floatMax(f64);
    inputs.molar_masses.carbon_g_per_mol = 2;
    try std.testing.expectError(
        error.NonFiniteFertilizerAccountingResult,
        calculateSourceOrder(inputs),
    );
}

fn sumFertilizer(state: fertilizer_dissolution.FertilizerState) f64 {
    var total: f64 = 0;
    inline for (@typeInfo(
        fertilizer_dissolution.FertilizerState,
    ).@"struct".fields) |field| total += @field(state, field.name);
    return total;
}
