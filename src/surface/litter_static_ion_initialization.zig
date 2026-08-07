const std = @import("std");

pub const IonInventories = struct {
    hydrogen_mol: f64,
    hydroxide_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
};

pub const ExistingConcentrations = struct {
    hydrogen_phosphate_mol_p_per_m3: f64,
    dihydrogen_phosphate_mol_p_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
};

pub const ActivityCoefficients = struct {
    monovalent: f64,
    divalent: f64,
    trivalent: f64,
};

pub const CarbonateEquilibrium = struct {
    bicarbonate_dissociation_mol_per_m3: f64,
    carbonate_combined_dissociation_mol2_per_m6: f64,
};

pub const Inputs = struct {
    litter_water_volume_m3: f64,
    concentration_floor_mol_per_m3: f64,
    carbon_molar_mass_g_c_per_mol: f64,
    dissolved_co2_g_c_per_m3: f64,
    ion_inventories: IonInventories,
    existing_concentrations: ExistingConcentrations,
    activity_coefficients: ActivityCoefficients,
    carbonate_equilibrium: CarbonateEquilibrium,
};

pub const Concentrations = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    chloride_mol_per_m3: f64,
    co2_mol_c_per_m3: f64,
    bicarbonate_mol_c_per_m3: f64,
    carbonate_mol_c_per_m3: f64,
    hydrogen_phosphate_mol_p_per_m3: f64,
    dihydrogen_phosphate_mol_p_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
};

pub const Activities = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    co2_mol_c_per_m3: f64,
    bicarbonate_mol_c_per_m3: f64,
    carbonate_mol_c_per_m3: f64,
    hydrogen_phosphate_mol_p_per_m3: f64,
    dihydrogen_phosphate_mol_p_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
};

pub const Result = struct {
    concentrations: Concentrations,
    activities: Activities,
};

/// Direct source-order translation of SOLUTE.F lines 4638--4666.
///
/// This is the one-cell surface-litter initialization used when salt
/// equilibria remain static. Array index `(0, NY, NX)` is resolved by the
/// caller; this pure kernel neither allocates nor mutates authoritative state.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    const volume = inputs.litter_water_volume_m3;
    const floor = inputs.concentration_floor_mol_per_m3;
    const inventories = inputs.ion_inventories;

    // SOLUTE.F 4638--4647.
    const hydrogen = @max(floor, inventories.hydrogen_mol / volume);
    const hydroxide = @max(floor, inventories.hydroxide_mol / volume);
    const aluminum = @max(floor, inventories.aluminum_mol / volume);
    const iron = @max(floor, inventories.iron_mol / volume);
    const calcium = @max(floor, inventories.calcium_mol / volume);
    const magnesium = @max(floor, inventories.magnesium_mol / volume);
    const sodium = @max(floor, inventories.sodium_mol / volume);
    const potassium = @max(floor, inventories.potassium_mol / volume);
    const sulfate = @max(floor, inventories.sulfate_mol / volume);
    const chloride = @max(floor, inventories.chloride_mol / volume);

    // SOLUTE.F 4648--4650. CCO2S is already an aqueous concentration, so
    // unlike the ion inventories it is not divided by litter-water volume.
    const co2 = @max(
        floor,
        inputs.dissolved_co2_g_c_per_m3 /
            inputs.carbon_molar_mass_g_c_per_mol,
    );
    const bicarbonate = @max(
        floor,
        co2 *
            inputs.carbonate_equilibrium.bicarbonate_dissociation_mol_per_m3 /
            hydrogen,
    );
    const carbonate = @max(
        floor,
        co2 *
            inputs.carbonate_equilibrium
                .carbonate_combined_dissociation_mol2_per_m6 /
            std.math.pow(f64, hydrogen, 2.0),
    );

    const existing = inputs.existing_concentrations;
    const coefficients = inputs.activity_coefficients;
    const concentrations: Concentrations = .{
        .hydrogen_mol_per_m3 = hydrogen,
        .hydroxide_mol_per_m3 = hydroxide,
        .aluminum_mol_per_m3 = aluminum,
        .iron_mol_per_m3 = iron,
        .calcium_mol_per_m3 = calcium,
        .magnesium_mol_per_m3 = magnesium,
        .sodium_mol_per_m3 = sodium,
        .potassium_mol_per_m3 = potassium,
        .sulfate_mol_per_m3 = sulfate,
        .chloride_mol_per_m3 = chloride,
        .co2_mol_c_per_m3 = co2,
        .bicarbonate_mol_c_per_m3 = bicarbonate,
        .carbonate_mol_c_per_m3 = carbonate,
        .hydrogen_phosphate_mol_p_per_m3 = existing.hydrogen_phosphate_mol_p_per_m3,
        .dihydrogen_phosphate_mol_p_per_m3 = existing.dihydrogen_phosphate_mol_p_per_m3,
        .ammonium_mol_n_per_m3 = existing.ammonium_mol_n_per_m3,
        .ammonia_mol_n_per_m3 = existing.ammonia_mol_n_per_m3,
    };

    // SOLUTE.F 4651--4666.
    const activities: Activities = .{
        .hydrogen_mol_per_m3 = hydrogen * coefficients.monovalent,
        .hydroxide_mol_per_m3 = hydroxide * coefficients.monovalent,
        .aluminum_mol_per_m3 = aluminum * coefficients.trivalent,
        .iron_mol_per_m3 = iron * coefficients.trivalent,
        .calcium_mol_per_m3 = calcium * coefficients.divalent,
        .magnesium_mol_per_m3 = magnesium * coefficients.divalent,
        .sodium_mol_per_m3 = sodium * coefficients.monovalent,
        .potassium_mol_per_m3 = potassium * coefficients.monovalent,
        .sulfate_mol_per_m3 = sulfate * coefficients.divalent,
        .co2_mol_c_per_m3 = co2,
        .bicarbonate_mol_c_per_m3 = bicarbonate * coefficients.monovalent,
        .carbonate_mol_c_per_m3 = carbonate * coefficients.divalent,
        .hydrogen_phosphate_mol_p_per_m3 = existing.hydrogen_phosphate_mol_p_per_m3 *
            coefficients.divalent,
        .dihydrogen_phosphate_mol_p_per_m3 = existing.dihydrogen_phosphate_mol_p_per_m3 *
            coefficients.monovalent,
        .ammonium_mol_n_per_m3 = existing.ammonium_mol_n_per_m3 * coefficients.monovalent,
        .ammonia_mol_n_per_m3 = existing.ammonia_mol_n_per_m3,
    };

    const result: Result = .{
        .concentrations = concentrations,
        .activities = activities,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.litter_water_volume_m3,
        inputs.concentration_floor_mol_per_m3,
        inputs.carbon_molar_mass_g_c_per_mol,
        inputs.dissolved_co2_g_c_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.InvalidSurfaceLitterStaticIonInput;
    }
    if (inputs.litter_water_volume_m3 <= 0 or
        inputs.concentration_floor_mol_per_m3 <= 0 or
        inputs.carbon_molar_mass_g_c_per_mol <= 0 or
        inputs.dissolved_co2_g_c_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterStaticIonInput;
    }
    try validateNonnegativeStruct(inputs.ion_inventories);
    try validateNonnegativeStruct(inputs.existing_concentrations);
    try validatePositiveStruct(inputs.activity_coefficients);
    try validatePositiveStruct(inputs.carbonate_equilibrium);
}

fn validateNonnegativeStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value < 0)
            return error.InvalidSurfaceLitterStaticIonInput;
    }
}

fn validatePositiveStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value <= 0)
            return error.InvalidSurfaceLitterStaticIonInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| {
        const value = @field(result.concentrations, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterStaticIonResult;
        if (value < 0)
            return error.InvalidSurfaceLitterStaticIonResult;
    }
    inline for (@typeInfo(Activities).@"struct".fields) |field| {
        const value = @field(result.activities, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterStaticIonResult;
        if (value < 0)
            return error.InvalidSurfaceLitterStaticIonResult;
    }
}

fn testInputs() Inputs {
    return .{
        .litter_water_volume_m3 = 2,
        .concentration_floor_mol_per_m3 = 1.0e-32,
        .carbon_molar_mass_g_c_per_mol = 12,
        .dissolved_co2_g_c_per_m3 = 24,
        .ion_inventories = .{
            .hydrogen_mol = 0.4,
            .hydroxide_mol = 0.6,
            .aluminum_mol = 0.8,
            .iron_mol = 1.0,
            .calcium_mol = 1.2,
            .magnesium_mol = 1.4,
            .sodium_mol = 1.6,
            .potassium_mol = 1.8,
            .sulfate_mol = 2.0,
            .chloride_mol = 2.2,
        },
        .existing_concentrations = .{
            .hydrogen_phosphate_mol_p_per_m3 = 1.2,
            .dihydrogen_phosphate_mol_p_per_m3 = 1.4,
            .ammonium_mol_n_per_m3 = 1.6,
            .ammonia_mol_n_per_m3 = 1.8,
        },
        .activity_coefficients = .{
            .monovalent = 0.8,
            .divalent = 0.6,
            .trivalent = 0.4,
        },
        .carbonate_equilibrium = .{
            .bicarbonate_dissociation_mol_per_m3 = 0.05,
            .carbonate_combined_dissociation_mol2_per_m6 = 0.002,
        },
    };
}

test "SOLUTE static surface ions preserve every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const floor = inputs.concentration_floor_mol_per_m3;
    const volume = inputs.litter_water_volume_m3;
    const inventories = inputs.ion_inventories;
    const co2 = @max(
        floor,
        inputs.dissolved_co2_g_c_per_m3 /
            inputs.carbon_molar_mass_g_c_per_mol,
    );
    const hydrogen = @max(floor, inventories.hydrogen_mol / volume);

    const expected_concentrations: Concentrations = .{
        .hydrogen_mol_per_m3 = hydrogen,
        .hydroxide_mol_per_m3 = @max(floor, inventories.hydroxide_mol / volume),
        .aluminum_mol_per_m3 = @max(floor, inventories.aluminum_mol / volume),
        .iron_mol_per_m3 = @max(floor, inventories.iron_mol / volume),
        .calcium_mol_per_m3 = @max(floor, inventories.calcium_mol / volume),
        .magnesium_mol_per_m3 = @max(floor, inventories.magnesium_mol / volume),
        .sodium_mol_per_m3 = @max(floor, inventories.sodium_mol / volume),
        .potassium_mol_per_m3 = @max(floor, inventories.potassium_mol / volume),
        .sulfate_mol_per_m3 = @max(floor, inventories.sulfate_mol / volume),
        .chloride_mol_per_m3 = @max(floor, inventories.chloride_mol / volume),
        .co2_mol_c_per_m3 = co2,
        .bicarbonate_mol_c_per_m3 = @max(
            floor,
            co2 *
                inputs.carbonate_equilibrium
                    .bicarbonate_dissociation_mol_per_m3 /
                hydrogen,
        ),
        .carbonate_mol_c_per_m3 = @max(
            floor,
            co2 *
                inputs.carbonate_equilibrium
                    .carbonate_combined_dissociation_mol2_per_m6 /
                std.math.pow(f64, hydrogen, 2.0),
        ),
        .hydrogen_phosphate_mol_p_per_m3 = inputs.existing_concentrations.hydrogen_phosphate_mol_p_per_m3,
        .dihydrogen_phosphate_mol_p_per_m3 = inputs.existing_concentrations.dihydrogen_phosphate_mol_p_per_m3,
        .ammonium_mol_n_per_m3 = inputs.existing_concentrations.ammonium_mol_n_per_m3,
        .ammonia_mol_n_per_m3 = inputs.existing_concentrations.ammonia_mol_n_per_m3,
    };
    try std.testing.expectEqualDeep(
        expected_concentrations,
        result.concentrations,
    );
}

test "static surface ion activities preserve source valence mapping" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const concentrations = result.concentrations;
    const activities = result.activities;
    const coefficients = inputs.activity_coefficients;

    try std.testing.expectEqual(
        concentrations.hydrogen_mol_per_m3 * coefficients.monovalent,
        activities.hydrogen_mol_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.aluminum_mol_per_m3 * coefficients.trivalent,
        activities.aluminum_mol_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.calcium_mol_per_m3 * coefficients.divalent,
        activities.calcium_mol_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.hydrogen_phosphate_mol_p_per_m3 *
            coefficients.divalent,
        activities.hydrogen_phosphate_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.dihydrogen_phosphate_mol_p_per_m3 *
            coefficients.monovalent,
        activities.dihydrogen_phosphate_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.co2_mol_c_per_m3,
        activities.co2_mol_c_per_m3,
    );
    try std.testing.expectEqual(
        concentrations.ammonia_mol_n_per_m3,
        activities.ammonia_mol_n_per_m3,
    );
}

test "surface dissolved carbon dioxide is already a concentration" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    var larger_volume = inputs;
    larger_volume.litter_water_volume_m3 = 8;
    const larger_volume_result = try calculateSourceOrder(larger_volume);

    try std.testing.expectEqual(
        result.concentrations.co2_mol_c_per_m3,
        larger_volume_result.concentrations.co2_mol_c_per_m3,
    );
    try std.testing.expect(
        result.concentrations.calcium_mol_per_m3 !=
            larger_volume_result.concentrations.calcium_mol_per_m3,
    );
}

test "static surface ion floor prevents zero hydrogen denominator" {
    var inputs = testInputs();
    inputs.ion_inventories.hydrogen_mol = 0;
    inputs.ion_inventories.hydroxide_mol = 0;
    const result = try calculateSourceOrder(inputs);
    const floor = inputs.concentration_floor_mol_per_m3;

    try std.testing.expectEqual(
        floor,
        result.concentrations.hydrogen_mol_per_m3,
    );
    try std.testing.expectEqual(
        floor,
        result.concentrations.hydroxide_mol_per_m3,
    );
    try std.testing.expect(std.math.isFinite(
        result.concentrations.bicarbonate_mol_c_per_m3,
    ));
    try std.testing.expect(std.math.isFinite(
        result.concentrations.carbonate_mol_c_per_m3,
    ));
}

test "static surface ion initialization rejects invalid inputs" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticIonInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.concentration_floor_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticIonInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.ion_inventories.calcium_mol = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticIonInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.activity_coefficients.divalent = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticIonInput,
        calculateSourceOrder(inputs),
    );
}

test "static surface ion initialization rejects non-finite results" {
    var inputs = testInputs();
    inputs.dissolved_co2_g_c_per_m3 = std.math.floatMax(f64);
    inputs.carbonate_equilibrium.bicarbonate_dissociation_mol_per_m3 =
        std.math.floatMax(f64);

    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticIonResult,
        calculateSourceOrder(inputs),
    );
}
