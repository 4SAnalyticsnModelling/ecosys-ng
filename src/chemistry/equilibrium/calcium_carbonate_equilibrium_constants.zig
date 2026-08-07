const std = @import("std");

pub const BaseConstants = struct {
    water_dissociation: f64,
    calcium_carbonate_solubility: f64,
    bicarbonate_dissociation: f64,
    dissolved_carbon_dioxide_dissociation: f64,
};

pub const DerivedConstants = struct {
    /// Legacy DPCO3 from SOLUTE.F line 72. Carbonate dissociation product
    /// ((mol m-3)^2) derived from dissolved CO2 and bicarbonate constants.
    carbonate_dissociation_product_mol2_per_m6: f64,
    /// Legacy SHCAC1.
    bicarbonate_proton_equilibrium: f64,
    /// Legacy SYCAC1.
    bicarbonate_hydroxide_equilibrium: f64,
    /// Legacy SHCAC2.
    carbon_dioxide_proton_equilibrium: f64,
    /// Legacy SYCAC2.
    carbon_dioxide_hydroxide_equilibrium: f64,
};

/// Exact SOLUTE carbonate and calcium-carbonate derivation from SOLUTE.F
/// lines 72 and 78--79.
///
/// Runtime inputs replace the source PARAMETER declarations. The constants
/// retain SOLUTE's concentration-based equilibrium convention and therefore
/// do not all have the same dimensions.
pub fn derive(base: BaseConstants) !DerivedConstants {
    inline for (@typeInfo(BaseConstants).@"struct".fields) |field| {
        const value = @field(base, field.name);
        if (!std.math.isFinite(value) or value <= 0.0)
            return error.InvalidCalciumCarbonateEquilibriumConstant;
    }

    // Preserve the source dependency order: DPCO3 is assigned before the
    // aluminum/iron chains and the later SHCAC/SYCAC calcium constants.
    const carbonate_dissociation_product =
        base.dissolved_carbon_dioxide_dissociation *
        base.bicarbonate_dissociation;
    const bicarbonate_proton =
        base.calcium_carbonate_solubility /
        base.bicarbonate_dissociation;
    const bicarbonate_hydroxide =
        bicarbonate_proton * base.water_dissociation;
    const carbon_dioxide_proton =
        bicarbonate_proton /
        base.dissolved_carbon_dioxide_dissociation;
    const carbon_dioxide_hydroxide =
        carbon_dioxide_proton *
        std.math.pow(f64, base.water_dissociation, 2.0);

    const result: DerivedConstants = .{
        .carbonate_dissociation_product_mol2_per_m6 = carbonate_dissociation_product,
        .bicarbonate_proton_equilibrium = bicarbonate_proton,
        .bicarbonate_hydroxide_equilibrium = bicarbonate_hydroxide,
        .carbon_dioxide_proton_equilibrium = carbon_dioxide_proton,
        .carbon_dioxide_hydroxide_equilibrium = carbon_dioxide_hydroxide,
    };
    inline for (@typeInfo(DerivedConstants).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCalciumCarbonateEquilibriumConstant;
    return result;
}

const source_defaults: BaseConstants = .{
    .water_dissociation = 1.0e-8,
    .calcium_carbonate_solubility = 3.3e-3,
    .bicarbonate_dissociation = 5.6e-8,
    .dissolved_carbon_dioxide_dissociation = 4.2e-4,
};

test "runtime derivation preserves SOLUTE calcium carbonate order" {
    const result = try derive(source_defaults);

    try std.testing.expectEqual(
        source_defaults.dissolved_carbon_dioxide_dissociation *
            source_defaults.bicarbonate_dissociation,
        result.carbonate_dissociation_product_mol2_per_m6,
    );
    try std.testing.expectApproxEqRel(
        source_defaults.calcium_carbonate_solubility /
            source_defaults.bicarbonate_dissociation,
        result.bicarbonate_proton_equilibrium,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.bicarbonate_proton_equilibrium *
            source_defaults.water_dissociation,
        result.bicarbonate_hydroxide_equilibrium,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.carbon_dioxide_proton_equilibrium *
            std.math.pow(f64, source_defaults.water_dissociation, 2.0),
        result.carbon_dioxide_hydroxide_equilibrium,
        1.0e-15,
    );
}

test "runtime derivation rejects unsafe calcium carbonate constants" {
    var base = source_defaults;
    base.bicarbonate_dissociation = 0.0;
    try std.testing.expectError(
        error.InvalidCalciumCarbonateEquilibriumConstant,
        derive(base),
    );

    base = source_defaults;
    base.calcium_carbonate_solubility = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidCalciumCarbonateEquilibriumConstant,
        derive(base),
    );
}
