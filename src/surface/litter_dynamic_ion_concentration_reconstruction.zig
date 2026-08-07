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
    bicarbonate_mol: f64,
    carbonate_mol: f64,
};

pub const IonConcentrations = struct {
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
    carbon_dioxide_mol_per_m3: f64,
    bicarbonate_mol_per_m3: f64,
    carbonate_mol_per_m3: f64,
};

pub const FloorsApplied = struct {
    aluminum: bool,
    iron: bool,
    calcium: bool,
    magnesium: bool,
    sodium: bool,
    potassium: bool,
    sulfate: bool,
    chloride: bool,
    carbon_dioxide: bool,
    bicarbonate: bool,
    carbonate: bool,
};

pub const Inputs = struct {
    litter_water_volume_m3: f64,
    external_hydrogen_production_mol_per_step: f64,
    inventories: IonInventories,
    dissolved_carbon_dioxide_g_c_per_m3: f64,
    carbon_molar_mass_g_c_per_mol: f64,
    aqueous_concentration_floor_mol_per_m3: f64,
};

pub const Result = struct {
    external_hydrogen_input_mol_per_m3_step: f64,
    concentrations: IonConcentrations,
    floors_applied: FloorsApplied,
};

/// Direct source-order translation of SOLUTE.F lines 4193--4236.
///
/// This pure kernel reconstructs the dynamic-salt aqueous state for one
/// caller-selected wet horizontal cell. It allocates and mutates no state.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);

    const water = inputs.litter_water_volume_m3;
    const floor = inputs.aqueous_concentration_floor_mol_per_m3;
    const inventory = inputs.inventories;

    // SOLUTE.F 4223--4236. Hydrogen and hydroxide are not floored here.
    const external_hydrogen =
        inputs.external_hydrogen_production_mol_per_step / water;
    const hydrogen = inventory.hydrogen_mol / water;
    const hydroxide = inventory.hydroxide_mol / water;
    const unconstrained_aluminum = inventory.aluminum_mol / water;
    const unconstrained_iron = inventory.iron_mol / water;
    const unconstrained_calcium = inventory.calcium_mol / water;
    const unconstrained_magnesium = inventory.magnesium_mol / water;
    const unconstrained_sodium = inventory.sodium_mol / water;
    const unconstrained_potassium = inventory.potassium_mol / water;
    const unconstrained_sulfate = inventory.sulfate_mol / water;
    const unconstrained_chloride = inventory.chloride_mol / water;
    const unconstrained_carbon_dioxide =
        inputs.dissolved_carbon_dioxide_g_c_per_m3 /
        inputs.carbon_molar_mass_g_c_per_mol;
    const unconstrained_bicarbonate = inventory.bicarbonate_mol / water;
    const unconstrained_carbonate = inventory.carbonate_mol / water;

    const result: Result = .{
        .external_hydrogen_input_mol_per_m3_step = external_hydrogen,
        .concentrations = .{
            .hydrogen_mol_per_m3 = hydrogen,
            .hydroxide_mol_per_m3 = hydroxide,
            .aluminum_mol_per_m3 = @max(floor, unconstrained_aluminum),
            .iron_mol_per_m3 = @max(floor, unconstrained_iron),
            .calcium_mol_per_m3 = @max(floor, unconstrained_calcium),
            .magnesium_mol_per_m3 = @max(floor, unconstrained_magnesium),
            .sodium_mol_per_m3 = @max(floor, unconstrained_sodium),
            .potassium_mol_per_m3 = @max(floor, unconstrained_potassium),
            .sulfate_mol_per_m3 = @max(floor, unconstrained_sulfate),
            .chloride_mol_per_m3 = @max(floor, unconstrained_chloride),
            .carbon_dioxide_mol_per_m3 = @max(floor, unconstrained_carbon_dioxide),
            .bicarbonate_mol_per_m3 = @max(floor, unconstrained_bicarbonate),
            .carbonate_mol_per_m3 = @max(floor, unconstrained_carbonate),
        },
        .floors_applied = .{
            .aluminum = unconstrained_aluminum < floor,
            .iron = unconstrained_iron < floor,
            .calcium = unconstrained_calcium < floor,
            .magnesium = unconstrained_magnesium < floor,
            .sodium = unconstrained_sodium < floor,
            .potassium = unconstrained_potassium < floor,
            .sulfate = unconstrained_sulfate < floor,
            .chloride = unconstrained_chloride < floor,
            .carbon_dioxide = unconstrained_carbon_dioxide < floor,
            .bicarbonate = unconstrained_bicarbonate < floor,
            .carbonate = unconstrained_carbonate < floor,
        },
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.litter_water_volume_m3,
        inputs.external_hydrogen_production_mol_per_step,
        inputs.dissolved_carbon_dioxide_g_c_per_m3,
        inputs.carbon_molar_mass_g_c_per_mol,
        inputs.aqueous_concentration_floor_mol_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.InvalidSurfaceLitterDynamicIonConcentrationInput;
    }
    if (inputs.litter_water_volume_m3 <= 0 or
        inputs.dissolved_carbon_dioxide_g_c_per_m3 < 0 or
        inputs.carbon_molar_mass_g_c_per_mol <= 0 or
        inputs.aqueous_concentration_floor_mol_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterDynamicIonConcentrationInput;
    }
    inline for (@typeInfo(IonInventories).@"struct".fields) |field| {
        const value = @field(inputs.inventories, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterDynamicIonConcentrationInput;
    }
}

fn validateResult(result: Result) !void {
    if (!std.math.isFinite(result.external_hydrogen_input_mol_per_m3_step))
        return error.NonFiniteSurfaceLitterDynamicIonConcentration;
    inline for (@typeInfo(IonConcentrations).@"struct".fields) |field| {
        const value = @field(result.concentrations, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterDynamicIonConcentration;
        if (value < 0)
            return error.InvalidSurfaceLitterDynamicIonConcentrationResult;
    }
}

fn testInputs() Inputs {
    return .{
        .litter_water_volume_m3 = 2,
        .external_hydrogen_production_mol_per_step = -0.4,
        .inventories = .{
            .hydrogen_mol = 0.2,
            .hydroxide_mol = 0.4,
            .aluminum_mol = 0.6,
            .iron_mol = 0.8,
            .calcium_mol = 1,
            .magnesium_mol = 1.2,
            .sodium_mol = 1.4,
            .potassium_mol = 1.6,
            .sulfate_mol = 1.8,
            .chloride_mol = 2,
            .bicarbonate_mol = 2.2,
            .carbonate_mol = 2.4,
        },
        .dissolved_carbon_dioxide_g_c_per_m3 = 15,
        .carbon_molar_mass_g_c_per_mol = 12,
        .aqueous_concentration_floor_mol_per_m3 = 1.0e-20,
    };
}

test "SOLUTE dynamic surface ions preserve every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);
    const inventory = inputs.inventories;
    const water = inputs.litter_water_volume_m3;
    const floor = inputs.aqueous_concentration_floor_mol_per_m3;

    try std.testing.expectEqual(
        inputs.external_hydrogen_production_mol_per_step / water,
        result.external_hydrogen_input_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        inventory.hydrogen_mol / water,
        result.concentrations.hydrogen_mol_per_m3,
    );
    try std.testing.expectEqual(
        inventory.hydroxide_mol / water,
        result.concentrations.hydroxide_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.aluminum_mol / water),
        result.concentrations.aluminum_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.iron_mol / water),
        result.concentrations.iron_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.calcium_mol / water),
        result.concentrations.calcium_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.magnesium_mol / water),
        result.concentrations.magnesium_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.sodium_mol / water),
        result.concentrations.sodium_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.potassium_mol / water),
        result.concentrations.potassium_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.sulfate_mol / water),
        result.concentrations.sulfate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.chloride_mol / water),
        result.concentrations.chloride_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(
            floor,
            inputs.dissolved_carbon_dioxide_g_c_per_m3 /
                inputs.carbon_molar_mass_g_c_per_mol,
        ),
        result.concentrations.carbon_dioxide_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.bicarbonate_mol / water),
        result.concentrations.bicarbonate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.carbonate_mol / water),
        result.concentrations.carbonate_mol_per_m3,
    );
}

test "dynamic surface ion concentrations reconstruct extensive inventories" {
    const inputs = testInputs();
    const result = try calculate(inputs);

    try std.testing.expectApproxEqAbs(
        inputs.external_hydrogen_production_mol_per_step,
        result.external_hydrogen_input_mol_per_m3_step *
            inputs.litter_water_volume_m3,
        1.0e-14,
    );
    inline for (@typeInfo(IonInventories).@"struct".fields) |field| {
        const concentration_name =
            field.name[0 .. field.name.len - "_mol".len] ++ "_mol_per_m3";
        try std.testing.expectApproxEqAbs(
            @field(inputs.inventories, field.name),
            @field(result.concentrations, concentration_name) *
                inputs.litter_water_volume_m3,
            1.0e-14,
        );
    }
    try std.testing.expectApproxEqAbs(
        inputs.dissolved_carbon_dioxide_g_c_per_m3,
        result.concentrations.carbon_dioxide_mol_per_m3 *
            inputs.carbon_molar_mass_g_c_per_mol,
        1.0e-14,
    );
}

test "dynamic surface ion reconstruction conserves inventory as water changes" {
    const inputs = testInputs();
    const initial = try calculate(inputs);
    var reduced_water = inputs;
    reduced_water.litter_water_volume_m3 = 0.5;
    const concentrated = try calculate(reduced_water);

    try std.testing.expectEqual(
        initial.concentrations.calcium_mol_per_m3 *
            inputs.litter_water_volume_m3,
        concentrated.concentrations.calcium_mol_per_m3 *
            reduced_water.litter_water_volume_m3,
    );
    try std.testing.expectEqual(
        initial.concentrations.sulfate_mol_per_m3 *
            inputs.litter_water_volume_m3,
        concentrated.concentrations.sulfate_mol_per_m3 *
            reduced_water.litter_water_volume_m3,
    );
}

test "dynamic surface ion floors are explicit by species" {
    var inputs = testInputs();
    inputs.inventories = std.mem.zeroes(IonInventories);
    inputs.dissolved_carbon_dioxide_g_c_per_m3 = 0;
    const result = try calculate(inputs);

    try std.testing.expectEqual(
        @as(f64, 0),
        result.concentrations.hydrogen_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.concentrations.hydroxide_mol_per_m3,
    );
    inline for (@typeInfo(FloorsApplied).@"struct".fields) |field|
        try std.testing.expect(@field(result.floors_applied, field.name));
}

test "dynamic surface ion reconstruction rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterDynamicIonConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.inventories.calcium_mol = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterDynamicIonConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.dissolved_carbon_dioxide_g_c_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterDynamicIonConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.litter_water_volume_m3 = std.math.floatMin(f64);
    inputs.inventories.aluminum_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterDynamicIonConcentration,
        calculate(inputs),
    );
}
