const std = @import("std");

pub const Inputs = struct {
    hydrogen_concentration_mol_per_m3: f64,
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
    h2po4_concentration_mol_p_per_m3: f64,
    h2po4_activity_mol_p_per_m3: f64,
    hpo4_concentration_mol_p_per_m3: f64,
    hpo4_activity_mol_p_per_m3: f64,
    deprotonated_site_mol_per_Mg: f64,
    hydroxyl_site_mol_per_Mg: f64,
    protonated_site_mol_per_Mg: f64,
    adsorbed_h2po4_mol_p_per_Mg: f64,
    adsorbed_hpo4_mol_p_per_Mg: f64,
    monovalent_activity_coefficient: f64,
    divalent_activity_coefficient: f64,
};

pub const Parameters = struct {
    protonated_site_equilibrium_constant: f64,
    hydroxyl_site_equilibrium_constant: f64,
    h2po4_exchange_equilibrium_constant: f64,
    hpo4_exchange_equilibrium_constant: f64,
    water_activity_product_mol2_per_m6: f64,
    h2po4_dissociation_constant: f64,
    maximum_exchange_mol_per_Mg_step: f64,
    substrate_limit_fraction: f64,
};

pub const Flux = struct {
    protonated_to_hydroxyl_site_mol_per_Mg: f64,
    hydroxyl_to_deprotonated_site_mol_per_Mg: f64,
    h2po4_with_protonated_site_mol_p_per_Mg: f64,
    h2po4_with_hydroxyl_site_mol_p_per_Mg: f64,
    hpo4_with_hydroxyl_site_mol_p_per_Mg: f64,
};

/// Direct shared kernel for SOLUTE lines 979--1054 and 1106--1174. The
/// non-band H2PO4/protonated-site expression uses explicit H and OH
/// activities; after water projection their product equals the band branch's
/// source `DPH2O` expression.
pub fn calculate(inputs: Inputs, parameters: Parameters) !Flux {
    try validate(inputs, parameters);
    const protonated_to_hydroxyl = if (inputs.protonated_site_mol_per_Mg > 0 and inputs.hydrogen_activity_mol_per_m3 > 0) bounded(
        inputs.hydroxyl_site_mol_per_Mg - parameters.protonated_site_equilibrium_constant * inputs.protonated_site_mol_per_Mg / inputs.hydrogen_activity_mol_per_m3,
        parameters.maximum_exchange_mol_per_Mg_step,
        parameters.substrate_limit_fraction * @min(inputs.hydrogen_concentration_mol_per_m3, inputs.hydroxyl_site_mol_per_Mg),
        parameters.substrate_limit_fraction * inputs.protonated_site_mol_per_Mg,
    ) else 0;
    const hydroxyl_to_deprotonated = if (inputs.hydroxyl_site_mol_per_Mg > 0 and inputs.hydrogen_activity_mol_per_m3 > 0) bounded(
        inputs.deprotonated_site_mol_per_Mg - parameters.hydroxyl_site_equilibrium_constant * inputs.hydroxyl_site_mol_per_Mg / inputs.hydrogen_activity_mol_per_m3,
        parameters.maximum_exchange_mol_per_Mg_step,
        parameters.substrate_limit_fraction * @min(inputs.hydrogen_concentration_mol_per_m3, inputs.deprotonated_site_mol_per_Mg),
        parameters.substrate_limit_fraction * inputs.hydroxyl_site_mol_per_Mg,
    ) else 0;
    const h2po4_with_protonated = if (inputs.protonated_site_mol_per_Mg > 0) blk: {
        const equilibrium = parameters.h2po4_exchange_equilibrium_constant * inputs.hydrogen_activity_mol_per_m3 * inputs.hydroxide_activity_mol_per_m3 * inputs.adsorbed_h2po4_mol_p_per_Mg / inputs.protonated_site_mol_per_Mg;
        break :blk bounded((inputs.h2po4_activity_mol_p_per_m3 - equilibrium) / inputs.monovalent_activity_coefficient, parameters.maximum_exchange_mol_per_Mg_step, parameters.substrate_limit_fraction * @min(inputs.h2po4_concentration_mol_p_per_m3, inputs.protonated_site_mol_per_Mg), parameters.substrate_limit_fraction * inputs.adsorbed_h2po4_mol_p_per_Mg);
    } else 0;
    const h2po4_with_hydroxyl = if (inputs.hydroxyl_site_mol_per_Mg > 0) blk: {
        const equilibrium = parameters.h2po4_exchange_equilibrium_constant * inputs.hydroxide_activity_mol_per_m3 * inputs.adsorbed_h2po4_mol_p_per_Mg / inputs.hydroxyl_site_mol_per_Mg;
        break :blk bounded((inputs.h2po4_activity_mol_p_per_m3 - equilibrium) / inputs.monovalent_activity_coefficient, parameters.maximum_exchange_mol_per_Mg_step, parameters.substrate_limit_fraction * @min(inputs.h2po4_concentration_mol_p_per_m3, inputs.hydroxyl_site_mol_per_Mg), parameters.substrate_limit_fraction * inputs.adsorbed_h2po4_mol_p_per_Mg);
    } else 0;
    const hpo4_with_hydroxyl = if (inputs.hydroxyl_site_mol_per_Mg > 0) blk: {
        const effective_equilibrium = parameters.hpo4_exchange_equilibrium_constant * parameters.water_activity_product_mol2_per_m6 / parameters.h2po4_dissociation_constant;
        const equilibrium = effective_equilibrium * inputs.adsorbed_hpo4_mol_p_per_Mg / inputs.hydroxyl_site_mol_per_Mg;
        break :blk bounded((inputs.hpo4_activity_mol_p_per_m3 - equilibrium) / inputs.divalent_activity_coefficient, parameters.maximum_exchange_mol_per_Mg_step, parameters.substrate_limit_fraction * @min(inputs.hpo4_concentration_mol_p_per_m3, inputs.hydroxyl_site_mol_per_Mg), parameters.substrate_limit_fraction * inputs.adsorbed_hpo4_mol_p_per_Mg);
    } else 0;
    return .{ .protonated_to_hydroxyl_site_mol_per_Mg = protonated_to_hydroxyl, .hydroxyl_to_deprotonated_site_mol_per_Mg = hydroxyl_to_deprotonated, .h2po4_with_protonated_site_mol_p_per_Mg = h2po4_with_protonated, .h2po4_with_hydroxyl_site_mol_p_per_Mg = h2po4_with_hydroxyl, .hpo4_with_hydroxyl_site_mol_p_per_Mg = hpo4_with_hydroxyl };
}

fn bounded(driving_force: f64, maximum_rate: f64, precipitation_limit: f64, dissolution_limit: f64) f64 {
    return @max(-maximum_rate, -dissolution_limit, @min(maximum_rate, precipitation_limit, driving_force));
}

fn validate(inputs: Inputs, parameters: Parameters) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name)) or @field(inputs, field.name) < 0) return error.InvalidPhosphateExchangeInput;
    inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name)) or @field(parameters, field.name) < 0) return error.InvalidPhosphateExchangeParameter;
    if (inputs.monovalent_activity_coefficient <= 0 or inputs.divalent_activity_coefficient <= 0 or parameters.water_activity_product_mol2_per_m6 <= 0 or parameters.h2po4_dissociation_constant <= 0 or parameters.substrate_limit_fraction > 1) return error.InvalidPhosphateExchangeParameter;
}

test "phosphate exchange rates obey kinetic and substrate bounds" {
    const flux = try calculate(.{ .hydrogen_concentration_mol_per_m3 = 0.1, .hydrogen_activity_mol_per_m3 = 0.08, .hydroxide_activity_mol_per_m3 = 1e-6, .h2po4_concentration_mol_p_per_m3 = 0.2, .h2po4_activity_mol_p_per_m3 = 0.16, .hpo4_concentration_mol_p_per_m3 = 0.1, .hpo4_activity_mol_p_per_m3 = 0.06, .deprotonated_site_mol_per_Mg = 0.2, .hydroxyl_site_mol_per_Mg = 0.4, .protonated_site_mol_per_Mg = 0.3, .adsorbed_h2po4_mol_p_per_Mg = 0.05, .adsorbed_hpo4_mol_p_per_Mg = 0.05, .monovalent_activity_coefficient = 0.8, .divalent_activity_coefficient = 0.6 }, .{ .protonated_site_equilibrium_constant = 0.1, .hydroxyl_site_equilibrium_constant = 0.1, .h2po4_exchange_equilibrium_constant = 0.2, .hpo4_exchange_equilibrium_constant = 0.2, .water_activity_product_mol2_per_m6 = 1e-8, .h2po4_dissociation_constant = 1e-3, .maximum_exchange_mol_per_Mg_step = 0.01, .substrate_limit_fraction = 0.2 });
    inline for (@typeInfo(Flux).@"struct".fields) |field| try std.testing.expect(@abs(@field(flux, field.name)) <= 0.01);
}

test "empty phosphate exchange sites produce finite zero branches" {
    const flux = try calculate(.{ .hydrogen_concentration_mol_per_m3 = 0.1, .hydrogen_activity_mol_per_m3 = 0.1, .hydroxide_activity_mol_per_m3 = 1e-7, .h2po4_concentration_mol_p_per_m3 = 0.1, .h2po4_activity_mol_p_per_m3 = 0.1, .hpo4_concentration_mol_p_per_m3 = 0.1, .hpo4_activity_mol_p_per_m3 = 0.1, .deprotonated_site_mol_per_Mg = 0, .hydroxyl_site_mol_per_Mg = 0, .protonated_site_mol_per_Mg = 0, .adsorbed_h2po4_mol_p_per_Mg = 0, .adsorbed_hpo4_mol_p_per_Mg = 0, .monovalent_activity_coefficient = 1, .divalent_activity_coefficient = 1 }, .{ .protonated_site_equilibrium_constant = 1, .hydroxyl_site_equilibrium_constant = 1, .h2po4_exchange_equilibrium_constant = 1, .hpo4_exchange_equilibrium_constant = 1, .water_activity_product_mol2_per_m6 = 1e-8, .h2po4_dissociation_constant = 1e-3, .maximum_exchange_mol_per_Mg_step = 0.1, .substrate_limit_fraction = 0.2 });
    try std.testing.expectEqual(@as(f64, 0), flux.h2po4_with_protonated_site_mol_p_per_Mg);
}

test "phosphate exchange matches SOLUTE non-band source equations" {
    const inputs = Inputs{
        .hydrogen_concentration_mol_per_m3 = 0.1,
        .hydrogen_activity_mol_per_m3 = 0.08,
        .hydroxide_activity_mol_per_m3 = 1.25e-7,
        .h2po4_concentration_mol_p_per_m3 = 0.2,
        .h2po4_activity_mol_p_per_m3 = 0.16,
        .hpo4_concentration_mol_p_per_m3 = 0.1,
        .hpo4_activity_mol_p_per_m3 = 0.06,
        .deprotonated_site_mol_per_Mg = 0.2,
        .hydroxyl_site_mol_per_Mg = 0.4,
        .protonated_site_mol_per_Mg = 0.3,
        .adsorbed_h2po4_mol_p_per_Mg = 0.05,
        .adsorbed_hpo4_mol_p_per_Mg = 0.05,
        .monovalent_activity_coefficient = 0.8,
        .divalent_activity_coefficient = 0.6,
    };
    const parameters = Parameters{
        .protonated_site_equilibrium_constant = 0.1,
        .hydroxyl_site_equilibrium_constant = 0.02,
        .h2po4_exchange_equilibrium_constant = 0.2,
        .hpo4_exchange_equilibrium_constant = 0.2,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .h2po4_dissociation_constant = 1.0e-3,
        .maximum_exchange_mol_per_Mg_step = 1,
        .substrate_limit_fraction = 0.5,
    };
    const result = try calculate(inputs, parameters);

    const protonated_equilibrium =
        parameters.protonated_site_equilibrium_constant *
        inputs.protonated_site_mol_per_Mg /
        inputs.hydrogen_activity_mol_per_m3;
    const hydroxyl_equilibrium =
        parameters.hydroxyl_site_equilibrium_constant *
        inputs.hydroxyl_site_mol_per_Mg /
        inputs.hydrogen_activity_mol_per_m3;
    const h2po4_protonated_equilibrium =
        parameters.h2po4_exchange_equilibrium_constant *
        inputs.hydrogen_activity_mol_per_m3 *
        inputs.hydroxide_activity_mol_per_m3 *
        inputs.adsorbed_h2po4_mol_p_per_Mg /
        inputs.protonated_site_mol_per_Mg;
    const h2po4_hydroxyl_equilibrium =
        parameters.h2po4_exchange_equilibrium_constant *
        inputs.hydroxide_activity_mol_per_m3 *
        inputs.adsorbed_h2po4_mol_p_per_Mg /
        inputs.hydroxyl_site_mol_per_Mg;
    const hpo4_equilibrium =
        parameters.hpo4_exchange_equilibrium_constant *
        parameters.water_activity_product_mol2_per_m6 /
        parameters.h2po4_dissociation_constant *
        inputs.adsorbed_hpo4_mol_p_per_Mg /
        inputs.hydroxyl_site_mol_per_Mg;

    try std.testing.expectEqual(
        bounded(
            inputs.hydroxyl_site_mol_per_Mg - protonated_equilibrium,
            1,
            0.5 * @min(
                inputs.hydrogen_concentration_mol_per_m3,
                inputs.hydroxyl_site_mol_per_Mg,
            ),
            0.5 * inputs.protonated_site_mol_per_Mg,
        ),
        result.protonated_to_hydroxyl_site_mol_per_Mg,
    );
    try std.testing.expectEqual(
        bounded(
            inputs.deprotonated_site_mol_per_Mg - hydroxyl_equilibrium,
            1,
            0.5 * @min(
                inputs.hydrogen_concentration_mol_per_m3,
                inputs.deprotonated_site_mol_per_Mg,
            ),
            0.5 * inputs.hydroxyl_site_mol_per_Mg,
        ),
        result.hydroxyl_to_deprotonated_site_mol_per_Mg,
    );
    try std.testing.expectEqual(
        bounded(
            (inputs.h2po4_activity_mol_p_per_m3 -
                h2po4_protonated_equilibrium) /
                inputs.monovalent_activity_coefficient,
            1,
            0.5 * @min(
                inputs.h2po4_concentration_mol_p_per_m3,
                inputs.protonated_site_mol_per_Mg,
            ),
            0.5 * inputs.adsorbed_h2po4_mol_p_per_Mg,
        ),
        result.h2po4_with_protonated_site_mol_p_per_Mg,
    );
    try std.testing.expectEqual(
        bounded(
            (inputs.h2po4_activity_mol_p_per_m3 -
                h2po4_hydroxyl_equilibrium) /
                inputs.monovalent_activity_coefficient,
            1,
            0.5 * @min(
                inputs.h2po4_concentration_mol_p_per_m3,
                inputs.hydroxyl_site_mol_per_Mg,
            ),
            0.5 * inputs.adsorbed_h2po4_mol_p_per_Mg,
        ),
        result.h2po4_with_hydroxyl_site_mol_p_per_Mg,
    );
    try std.testing.expectEqual(
        bounded(
            (inputs.hpo4_activity_mol_p_per_m3 - hpo4_equilibrium) /
                inputs.divalent_activity_coefficient,
            1,
            0.5 * @min(
                inputs.hpo4_concentration_mol_p_per_m3,
                inputs.hydroxyl_site_mol_per_Mg,
            ),
            0.5 * inputs.adsorbed_hpo4_mol_p_per_Mg,
        ),
        result.hpo4_with_hydroxyl_site_mol_p_per_Mg,
    );
}
