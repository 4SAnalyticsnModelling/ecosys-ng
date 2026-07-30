const std = @import("std");

pub const Parameters = struct {
    carboxyl_capacity_mol_per_Mg_carbon: f64,
    carboxyl_dissociation_mol_per_m3: f64,
    negligible_site_concentration_mol_per_Mg: f64,
};

pub const Inputs = struct {
    cation_exchange_capacity_mol_charge: f64,
    effective_soil_mass_Mg: f64,
    soil_organic_carbon_g_c: f64,
    hydrogen_mol_per_m3: f64,
};

pub const Result = struct {
    cation_exchange_capacity_mol_charge_per_Mg: f64,
    total_carboxyl_sites_mol_per_Mg: f64,
    protonated_carboxyl_sites_mol_per_Mg: f64,
};

/// Direct translation of STARTE lines 400--403. The `1.0e-6` source factor
/// converts organic carbon from g C to Mg C before applying mol sites/Mg C.
pub fn calculate(inputs: Inputs, parameters: Parameters) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteCarboxylExchangeInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(parameters, field.name)))
            return error.NonFiniteCarboxylExchangeParameter;
    }
    if (inputs.cation_exchange_capacity_mol_charge < 0 or
        inputs.effective_soil_mass_Mg <= 0 or
        inputs.soil_organic_carbon_g_c < 0 or
        inputs.hydrogen_mol_per_m3 < 0 or
        parameters.carboxyl_capacity_mol_per_Mg_carbon < 0 or
        parameters.carboxyl_dissociation_mol_per_m3 <= 0 or
        parameters.negligible_site_concentration_mol_per_Mg < 0)
        return error.InvalidCarboxylExchangeInput;

    const total_carboxyl_sites_mol_per_Mg = @max(
        parameters.negligible_site_concentration_mol_per_Mg,
        parameters.carboxyl_capacity_mol_per_Mg_carbon *
            1.0e-6 *
            inputs.soil_organic_carbon_g_c /
            inputs.effective_soil_mass_Mg,
    );
    const result: Result = .{
        .cation_exchange_capacity_mol_charge_per_Mg = inputs.cation_exchange_capacity_mol_charge /
            inputs.effective_soil_mass_Mg,
        .total_carboxyl_sites_mol_per_Mg = total_carboxyl_sites_mol_per_Mg,
        .protonated_carboxyl_sites_mol_per_Mg = total_carboxyl_sites_mol_per_Mg *
            @min(
                1.0,
                inputs.hydrogen_mol_per_m3 /
                    parameters.carboxyl_dissociation_mol_per_m3,
            ),
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteCarboxylExchangeResult;
    }
    return result;
}

fn sourceParameters() Parameters {
    return .{
        .carboxyl_capacity_mol_per_Mg_carbon = 250,
        .carboxyl_dissociation_mol_per_m3 = 1.0e-2,
        .negligible_site_concentration_mol_per_Mg = 1.0e-48,
    };
}

test "STARTE soil carboxyl sites preserve source mass and occupancy order" {
    const result = try calculate(.{
        .cation_exchange_capacity_mol_charge = 60,
        .effective_soil_mass_Mg = 2,
        .soil_organic_carbon_g_c = 40_000,
        .hydrogen_mol_per_m3 = 2.5e-3,
    }, sourceParameters());
    try std.testing.expectEqual(
        @as(f64, 30),
        result.cation_exchange_capacity_mol_charge_per_Mg,
    );
    try std.testing.expectEqual(
        @as(f64, 5),
        result.total_carboxyl_sites_mol_per_Mg,
    );
    try std.testing.expectEqual(
        @as(f64, 1.25),
        result.protonated_carboxyl_sites_mol_per_Mg,
    );
}

test "STARTE soil carboxyl protonation is capped at all sites" {
    const result = try calculate(.{
        .cation_exchange_capacity_mol_charge = 1,
        .effective_soil_mass_Mg = 1,
        .soil_organic_carbon_g_c = 4000,
        .hydrogen_mol_per_m3 = 1,
    }, sourceParameters());
    try std.testing.expectEqual(
        result.total_carboxyl_sites_mol_per_Mg,
        result.protonated_carboxyl_sites_mol_per_Mg,
    );
}

test "STARTE soil carboxyl sites reject invalid mass" {
    const inputs: Inputs = .{
        .cation_exchange_capacity_mol_charge = 1,
        .effective_soil_mass_Mg = 0,
        .soil_organic_carbon_g_c = 1,
        .hydrogen_mol_per_m3 = 1,
    };
    try std.testing.expectError(
        error.InvalidCarboxylExchangeInput,
        calculate(inputs, sourceParameters()),
    );
}
