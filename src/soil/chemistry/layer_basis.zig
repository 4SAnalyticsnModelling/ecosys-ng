const std = @import("std");

pub const Parameters = struct {
    ammonium_extract_multiplier: f64,
    phosphate_extract_multiplier: f64,
    water_activity_product_mol2_per_m6: f64,
};

pub const Inputs = struct {
    soil_mass_megagrams: f64,
    pore_volume_m3: f64,
    field_capacity_water_m3: f64,
    negligible_soil_mass_megagrams: f64,
    cation_exchange_capacity_mol_charge_per_megagram: f64,
    anion_exchange_capacity_mol_charge_per_megagram: f64,
    soil_ph: f64,
    ammonium_mol_n_per_megagram: f64,
    nitrate_mol_n_per_megagram: f64,
    phosphate_mol_p_per_megagram: f64,
};

pub const Basis = struct {
    effective_soil_mass_megagrams: f64,
    soil_mass_to_field_capacity_ratio: f64,
    cation_exchange_capacity_mol_charge: f64,
    anion_exchange_capacity_mol_charge: f64,
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    ammonium_extract_mol_n_per_megagram: f64,
    nitrate_extract_mol_n_per_megagram: f64,
    phosphate_extract_mol_p_per_megagram: f64,
};

/// Direct translation of `starte.f` lines 179--191 and 208--211. This establishes
/// the layer basis used by the subsequent source-order aqueous speciation.
pub fn derive(inputs: Inputs, parameters: Parameters) !Basis {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSoilChemistryLayerInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(parameters, field.name)))
            return error.NonFiniteSoilChemistryLayerParameter;
    }
    if (inputs.soil_mass_megagrams < 0 or
        inputs.pore_volume_m3 <= 0 or
        inputs.field_capacity_water_m3 < 0 or
        inputs.negligible_soil_mass_megagrams < 0 or
        inputs.cation_exchange_capacity_mol_charge_per_megagram < 0 or
        inputs.anion_exchange_capacity_mol_charge_per_megagram < 0 or
        inputs.ammonium_mol_n_per_megagram < 0 or
        inputs.nitrate_mol_n_per_megagram < 0 or
        inputs.phosphate_mol_p_per_megagram < 0 or
        parameters.ammonium_extract_multiplier < 0 or
        parameters.phosphate_extract_multiplier < 0 or
        parameters.water_activity_product_mol2_per_m6 <= 0)
        return error.InvalidSoilChemistryLayerInput;

    const effective_soil_mass_megagrams =
        if (inputs.soil_mass_megagrams > inputs.negligible_soil_mass_megagrams)
            inputs.soil_mass_megagrams
        else
            inputs.pore_volume_m3;
    const mass_to_field_capacity =
        if (inputs.field_capacity_water_m3 > inputs.negligible_soil_mass_megagrams)
            @min(1.0, effective_soil_mass_megagrams / inputs.field_capacity_water_m3)
        else
            1.0;
    const hydrogen_mol_per_m3 =
        std.math.pow(f64, 10.0, -(inputs.soil_ph - 3.0));
    const hydroxide_mol_per_m3 =
        parameters.water_activity_product_mol2_per_m6 / hydrogen_mol_per_m3;

    const result: Basis = .{
        .effective_soil_mass_megagrams = effective_soil_mass_megagrams,
        .soil_mass_to_field_capacity_ratio = mass_to_field_capacity,
        .cation_exchange_capacity_mol_charge = inputs.cation_exchange_capacity_mol_charge_per_megagram *
            effective_soil_mass_megagrams,
        .anion_exchange_capacity_mol_charge = inputs.anion_exchange_capacity_mol_charge_per_megagram *
            effective_soil_mass_megagrams,
        .hydrogen_mol_per_m3 = hydrogen_mol_per_m3,
        .hydroxide_mol_per_m3 = hydroxide_mol_per_m3,
        .ammonium_extract_mol_n_per_megagram = parameters.ammonium_extract_multiplier *
            inputs.ammonium_mol_n_per_megagram,
        .nitrate_extract_mol_n_per_megagram = inputs.nitrate_mol_n_per_megagram,
        .phosphate_extract_mol_p_per_megagram = parameters.phosphate_extract_multiplier *
            inputs.phosphate_mol_p_per_megagram,
    };
    inline for (@typeInfo(Basis).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSoilChemistryLayerResult;
    }
    return result;
}

fn sourceParameters() Parameters {
    return .{
        .ammonium_extract_multiplier = 0.1,
        .phosphate_extract_multiplier = 0.01,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
    };
}

test "STARTE layer basis preserves source operation order and units" {
    const result = try derive(.{
        .soil_mass_megagrams = 2,
        .pore_volume_m3 = 4,
        .field_capacity_water_m3 = 5,
        .negligible_soil_mass_megagrams = 1.0e-12,
        .cation_exchange_capacity_mol_charge_per_megagram = 10,
        .anion_exchange_capacity_mol_charge_per_megagram = 3,
        .soil_ph = 7,
        .ammonium_mol_n_per_megagram = 8,
        .nitrate_mol_n_per_megagram = 6,
        .phosphate_mol_p_per_megagram = 4,
    }, sourceParameters());
    try std.testing.expectEqual(@as(f64, 2), result.effective_soil_mass_megagrams);
    try std.testing.expectEqual(@as(f64, 0.4), result.soil_mass_to_field_capacity_ratio);
    try std.testing.expectEqual(@as(f64, 20), result.cation_exchange_capacity_mol_charge);
    try std.testing.expectEqual(@as(f64, 6), result.anion_exchange_capacity_mol_charge);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), result.hydrogen_mol_per_m3, 1.0e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), result.hydroxide_mol_per_m3, 1.0e-18);
    try std.testing.expectEqual(@as(f64, 0.8), result.ammonium_extract_mol_n_per_megagram);
    try std.testing.expectEqual(@as(f64, 6), result.nitrate_extract_mol_n_per_megagram);
    try std.testing.expectEqual(@as(f64, 0.04), result.phosphate_extract_mol_p_per_megagram);
}

test "STARTE layer basis uses pore volume for a negligible soil mass" {
    const result = try derive(.{
        .soil_mass_megagrams = 0,
        .pore_volume_m3 = 0.75,
        .field_capacity_water_m3 = 0,
        .negligible_soil_mass_megagrams = 1.0e-12,
        .cation_exchange_capacity_mol_charge_per_megagram = 2,
        .anion_exchange_capacity_mol_charge_per_megagram = 1,
        .soil_ph = 6,
        .ammonium_mol_n_per_megagram = 0,
        .nitrate_mol_n_per_megagram = 0,
        .phosphate_mol_p_per_megagram = 0,
    }, sourceParameters());
    try std.testing.expectEqual(@as(f64, 0.75), result.effective_soil_mass_megagrams);
    try std.testing.expectEqual(@as(f64, 1), result.soil_mass_to_field_capacity_ratio);
    try std.testing.expectEqual(@as(f64, 1.5), result.cation_exchange_capacity_mol_charge);
}

test "STARTE layer basis rejects invalid state before publishing a result" {
    var inputs: Inputs = .{
        .soil_mass_megagrams = 1,
        .pore_volume_m3 = 1,
        .field_capacity_water_m3 = 1,
        .negligible_soil_mass_megagrams = 1.0e-12,
        .cation_exchange_capacity_mol_charge_per_megagram = 1,
        .anion_exchange_capacity_mol_charge_per_megagram = 1,
        .soil_ph = 7,
        .ammonium_mol_n_per_megagram = 1,
        .nitrate_mol_n_per_megagram = 1,
        .phosphate_mol_p_per_megagram = 1,
    };
    inputs.soil_mass_megagrams = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSoilChemistryLayerInput,
        derive(inputs, sourceParameters()),
    );
}
