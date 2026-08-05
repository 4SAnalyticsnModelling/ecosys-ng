const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static_concentrations,
    dynamic_equilibria,
};

pub const AqueousTotals = struct {
    carbonate_mol_c_per_step: f64,
    carbon_dioxide_mol_c_per_step: f64,
    bicarbonate_mol_c_per_step: f64,
    sulfate_mol_s_per_step: f64,
    hydrogen_mol_per_step: f64,
    hydroxide_mol_per_step: f64,
    water_mol_per_step: f64,
    buffered_water_mol_per_step: f64,
};

pub const ExchangeAndMineralTotals = struct {
    carboxyl_hydrogen_exchange_mol_per_step: f64,
    aluminum_hydroxide_mol_mineral_per_step: f64,
    iron_hydroxide_mol_mineral_per_step: f64,
    calcite_mol_mineral_per_step: f64,
    gypsum_mol_mineral_per_step: f64,
};

pub const Totals = struct {
    aqueous: AqueousTotals,
    exchange_and_minerals: ExchangeAndMineralTotals,
};

pub const DynamicRates = struct {
    carbonate_mol_c_per_m3_step: f64,
    carbon_dioxide_mol_c_per_m3_step: f64,
    bicarbonate_mol_c_per_m3_step: f64,
    sulfate_mol_s_per_m3_step: f64,
    carboxyl_hydrogen_exchange_mol_per_megagram_step: f64,
    aluminum_hydroxide_mol_mineral_per_m3_step: f64,
    iron_hydroxide_mol_mineral_per_m3_step: f64,
    calcite_mol_mineral_per_m3_step: f64,
    gypsum_mol_mineral_per_m3_step: f64,
    hydrogen_mol_per_m3_step: f64,
    hydroxide_mol_per_m3_step: f64,
    water_mol_per_m3_step: f64,
    hydrogen_hydroxide_equilibration_mol_per_m3_step: f64,
};

pub const StaticRates = struct {
    hydrogen_mol_per_m3_step: f64,
    hydroxide_mol_per_m3_step: f64,
    hydrogen_equilibration_mol_per_m3_step: f64,
    hydroxide_equilibration_mol_per_m3_step: f64,
};

pub const Rates = union(SaltEquilibriumMode) {
    static_concentrations: StaticRates,
    dynamic_equilibria: DynamicRates,
};

pub const Inputs = struct {
    existing_totals: Totals,
    rates: Rates,
    litter_water_volume_m3: f64,
    litter_dry_mass_megagrams: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5127--5166.
///
/// Concentration changes are converted to per-cell molar totals with water
/// volume; the one exchange change is converted with litter dry mass.
pub fn calculateSourceOrder(inputs: Inputs) !Totals {
    try validateInputs(inputs);
    var totals = inputs.existing_totals;
    const water_volume = inputs.litter_water_volume_m3;

    switch (inputs.rates) {
        // SOLUTE.F 5127--5140.
        .dynamic_equilibria => |rates| {
            totals.aqueous.carbonate_mol_c_per_step +=
                rates.carbonate_mol_c_per_m3_step * water_volume;
            totals.aqueous.carbon_dioxide_mol_c_per_step +=
                rates.carbon_dioxide_mol_c_per_m3_step * water_volume;
            totals.aqueous.bicarbonate_mol_c_per_step +=
                rates.bicarbonate_mol_c_per_m3_step * water_volume;
            totals.aqueous.sulfate_mol_s_per_step +=
                rates.sulfate_mol_s_per_m3_step * water_volume;
            totals.exchange_and_minerals
                .carboxyl_hydrogen_exchange_mol_per_step +=
                rates.carboxyl_hydrogen_exchange_mol_per_megagram_step *
                inputs.litter_dry_mass_megagrams;
            totals.exchange_and_minerals
                .aluminum_hydroxide_mol_mineral_per_step +=
                rates.aluminum_hydroxide_mol_mineral_per_m3_step *
                water_volume;
            totals.exchange_and_minerals
                .iron_hydroxide_mol_mineral_per_step +=
                rates.iron_hydroxide_mol_mineral_per_m3_step *
                water_volume;
            totals.exchange_and_minerals.calcite_mol_mineral_per_step +=
                rates.calcite_mol_mineral_per_m3_step * water_volume;
            totals.exchange_and_minerals.gypsum_mol_mineral_per_step +=
                rates.gypsum_mol_mineral_per_m3_step * water_volume;
            totals.aqueous.hydrogen_mol_per_step +=
                (rates.hydrogen_mol_per_m3_step -
                    rates
                        .hydrogen_hydroxide_equilibration_mol_per_m3_step) *
                water_volume;
            totals.aqueous.hydroxide_mol_per_step +=
                (rates.hydroxide_mol_per_m3_step -
                    rates
                        .hydrogen_hydroxide_equilibration_mol_per_m3_step) *
                water_volume;
            totals.aqueous.water_mol_per_step +=
                rates.water_mol_per_m3_step * water_volume;
            totals.aqueous.buffered_water_mol_per_step +=
                rates.hydrogen_hydroxide_equilibration_mol_per_m3_step *
                water_volume;
        },

        // SOLUTE.F 5162--5166.
        .static_concentrations => |rates| {
            totals.aqueous.hydrogen_mol_per_step +=
                (rates.hydrogen_mol_per_m3_step +
                    rates.hydrogen_equilibration_mol_per_m3_step) *
                water_volume;
            totals.aqueous.hydroxide_mol_per_step +=
                (rates.hydroxide_mol_per_m3_step +
                    rates.hydroxide_equilibration_mol_per_m3_step) *
                water_volume;
            totals.aqueous.buffered_water_mol_per_step +=
                rates.hydrogen_equilibration_mol_per_m3_step *
                water_volume;
        },
    }

    try validateTotals(totals, error.NonFiniteSurfaceLitterAccumulationResult);
    return totals;
}

fn validateInputs(inputs: Inputs) !void {
    try validateTotals(
        inputs.existing_totals,
        error.InvalidSurfaceLitterAccumulationInput,
    );
    inline for (.{
        inputs.litter_water_volume_m3,
        inputs.litter_dry_mass_megagrams,
    }) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterAccumulationInput;
    }
    switch (inputs.rates) {
        inline else => |rates| try validateFiniteStruct(
            rates,
            error.InvalidSurfaceLitterAccumulationInput,
        ),
    }
}

fn validateTotals(totals: Totals, failure: anyerror) !void {
    try validateFiniteStruct(totals.aqueous, failure);
    try validateFiniteStruct(totals.exchange_and_minerals, failure);
}

fn validateFiniteStruct(value: anytype, failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(value, field.name)))
            return failure;
    }
}

fn testTotals() Totals {
    return .{
        .aqueous = .{
            .carbonate_mol_c_per_step = 1,
            .carbon_dioxide_mol_c_per_step = 2,
            .bicarbonate_mol_c_per_step = 3,
            .sulfate_mol_s_per_step = 4,
            .hydrogen_mol_per_step = 5,
            .hydroxide_mol_per_step = 6,
            .water_mol_per_step = 7,
            .buffered_water_mol_per_step = 8,
        },
        .exchange_and_minerals = .{
            .carboxyl_hydrogen_exchange_mol_per_step = 9,
            .aluminum_hydroxide_mol_mineral_per_step = 10,
            .iron_hydroxide_mol_mineral_per_step = 11,
            .calcite_mol_mineral_per_step = 12,
            .gypsum_mol_mineral_per_step = 13,
        },
    };
}

fn dynamicInputs() Inputs {
    return .{
        .existing_totals = testTotals(),
        .rates = .{ .dynamic_equilibria = .{
            .carbonate_mol_c_per_m3_step = 0.1,
            .carbon_dioxide_mol_c_per_m3_step = 0.2,
            .bicarbonate_mol_c_per_m3_step = 0.3,
            .sulfate_mol_s_per_m3_step = 0.4,
            .carboxyl_hydrogen_exchange_mol_per_megagram_step = 0.5,
            .aluminum_hydroxide_mol_mineral_per_m3_step = 0.6,
            .iron_hydroxide_mol_mineral_per_m3_step = 0.7,
            .calcite_mol_mineral_per_m3_step = 0.8,
            .gypsum_mol_mineral_per_m3_step = 0.9,
            .hydrogen_mol_per_m3_step = 1.0,
            .hydroxide_mol_per_m3_step = 1.1,
            .water_mol_per_m3_step = 1.2,
            .hydrogen_hydroxide_equilibration_mol_per_m3_step = 0.2,
        } },
        .litter_water_volume_m3 = 2,
        .litter_dry_mass_megagrams = 3,
    };
}

test "SOLUTE dynamic surface accumulation preserves every source update" {
    const inputs = dynamicInputs();
    const result = try calculateSourceOrder(inputs);
    const initial = inputs.existing_totals;
    const rates = inputs.rates.dynamic_equilibria;
    const volume = inputs.litter_water_volume_m3;

    try std.testing.expectEqual(
        initial.aqueous.carbonate_mol_c_per_step +
            rates.carbonate_mol_c_per_m3_step * volume,
        result.aqueous.carbonate_mol_c_per_step,
    );
    try std.testing.expectEqual(
        initial.aqueous.carbon_dioxide_mol_c_per_step +
            rates.carbon_dioxide_mol_c_per_m3_step * volume,
        result.aqueous.carbon_dioxide_mol_c_per_step,
    );
    try std.testing.expectEqual(
        initial.exchange_and_minerals
            .carboxyl_hydrogen_exchange_mol_per_step +
            rates.carboxyl_hydrogen_exchange_mol_per_megagram_step *
                inputs.litter_dry_mass_megagrams,
        result.exchange_and_minerals
            .carboxyl_hydrogen_exchange_mol_per_step,
    );
    try std.testing.expectEqual(
        initial.exchange_and_minerals.gypsum_mol_mineral_per_step +
            rates.gypsum_mol_mineral_per_m3_step * volume,
        result.exchange_and_minerals.gypsum_mol_mineral_per_step,
    );
    try std.testing.expectEqual(
        initial.aqueous.water_mol_per_step +
            rates.water_mol_per_m3_step * volume,
        result.aqueous.water_mol_per_step,
    );
}

test "dynamic H and OH buffering bookkeeping closes" {
    const inputs = dynamicInputs();
    const result = try calculateSourceOrder(inputs);
    const initial = inputs.existing_totals.aqueous;
    const rates = inputs.rates.dynamic_equilibria;
    const volume = inputs.litter_water_volume_m3;
    const hydrogen_delta =
        result.aqueous.hydrogen_mol_per_step -
        initial.hydrogen_mol_per_step;
    const hydroxide_delta =
        result.aqueous.hydroxide_mol_per_step -
        initial.hydroxide_mol_per_step;
    const buffer_delta =
        result.aqueous.buffered_water_mol_per_step -
        initial.buffered_water_mol_per_step;

    try std.testing.expectApproxEqAbs(
        rates.hydrogen_mol_per_m3_step * volume,
        hydrogen_delta + buffer_delta,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        rates.hydroxide_mol_per_m3_step * volume,
        hydroxide_delta + buffer_delta,
        1.0e-15,
    );
}

test "static accumulation updates only H OH and buffered water" {
    const inputs: Inputs = .{
        .existing_totals = testTotals(),
        .rates = .{ .static_concentrations = .{
            .hydrogen_mol_per_m3_step = 0.4,
            .hydroxide_mol_per_m3_step = 0.5,
            .hydrogen_equilibration_mol_per_m3_step = 0.6,
            .hydroxide_equilibration_mol_per_m3_step = 0.7,
        } },
        .litter_water_volume_m3 = 2,
        .litter_dry_mass_megagrams = 3,
    };
    const result = try calculateSourceOrder(inputs);
    const initial = inputs.existing_totals;
    const rates = inputs.rates.static_concentrations;

    try std.testing.expectEqual(
        initial.aqueous.hydrogen_mol_per_step +
            (rates.hydrogen_mol_per_m3_step +
                rates.hydrogen_equilibration_mol_per_m3_step) *
                inputs.litter_water_volume_m3,
        result.aqueous.hydrogen_mol_per_step,
    );
    try std.testing.expectEqual(
        initial.aqueous.hydroxide_mol_per_step +
            (rates.hydroxide_mol_per_m3_step +
                rates.hydroxide_equilibration_mol_per_m3_step) *
                inputs.litter_water_volume_m3,
        result.aqueous.hydroxide_mol_per_step,
    );
    try std.testing.expectEqual(
        initial.aqueous.buffered_water_mol_per_step +
            rates.hydrogen_equilibration_mol_per_m3_step *
                inputs.litter_water_volume_m3,
        result.aqueous.buffered_water_mol_per_step,
    );
    try std.testing.expectEqual(
        initial.aqueous.carbonate_mol_c_per_step,
        result.aqueous.carbonate_mol_c_per_step,
    );
    try std.testing.expectEqualDeep(
        initial.exchange_and_minerals,
        result.exchange_and_minerals,
    );
}

test "zero litter dimensions leave all totals unchanged" {
    var inputs = dynamicInputs();
    inputs.litter_water_volume_m3 = 0;
    inputs.litter_dry_mass_megagrams = 0;

    try std.testing.expectEqualDeep(
        inputs.existing_totals,
        try calculateSourceOrder(inputs),
    );
}

test "surface accumulation rejects invalid input and overflow" {
    var inputs = dynamicInputs();
    inputs.litter_water_volume_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterAccumulationInput,
        calculateSourceOrder(inputs),
    );

    inputs = dynamicInputs();
    inputs.rates.dynamic_equilibria.sulfate_mol_s_per_m3_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterAccumulationInput,
        calculateSourceOrder(inputs),
    );

    inputs = dynamicInputs();
    inputs.litter_water_volume_m3 = std.math.floatMax(f64);
    inputs.rates.dynamic_equilibria.carbonate_mol_c_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterAccumulationResult,
        calculateSourceOrder(inputs),
    );
}
