const std = @import("std");

pub const BaseConstants = struct {
    water_dissociation: f64,
    calcium_hydrogen_phosphate_solubility: f64,
    hydroxyapatite_solubility: f64,
    hydrogen_phosphate_dissociation: f64,
    dihydrogen_phosphate_dissociation: f64,
};

pub const DerivedConstants = struct {
    /// Legacy SHCAD2.
    calcium_hydrogen_phosphate_proton_equilibrium: f64,
    /// Legacy SYCAD2.
    calcium_hydrogen_phosphate_hydroxide_equilibrium: f64,
    /// Legacy SHCAH1.
    hydroxyapatite_hpo4_proton_equilibrium: f64,
    /// Legacy SYCAH1.
    hydroxyapatite_hpo4_hydroxide_equilibrium: f64,
    /// Legacy SHCAH2.
    hydroxyapatite_h2po4_proton_equilibrium: f64,
    /// Legacy SYCAH2.
    hydroxyapatite_h2po4_hydroxide_equilibrium: f64,
};

/// Exact SOLUTE calcium-phosphate derivation from solute.f:93-95.
pub fn derive(base: BaseConstants) !DerivedConstants {
    inline for (@typeInfo(BaseConstants).@"struct".fields) |field| {
        const value = @field(base, field.name);
        if (!std.math.isFinite(value) or value <= 0.0)
            return error.InvalidCalciumPhosphateEquilibriumConstant;
    }

    const hydrogen_phosphate_proton =
        base.calcium_hydrogen_phosphate_solubility /
        base.dihydrogen_phosphate_dissociation;
    const hydrogen_phosphate_hydroxide =
        hydrogen_phosphate_proton * base.water_dissociation;
    const hydroxyapatite_hpo4_proton =
        base.hydroxyapatite_solubility /
        (base.water_dissociation *
            std.math.pow(f64, base.hydrogen_phosphate_dissociation, 3.0));
    const hydroxyapatite_hpo4_hydroxide =
        hydroxyapatite_hpo4_proton *
        std.math.pow(f64, base.water_dissociation, 4.0);
    const hydroxyapatite_h2po4_proton =
        hydroxyapatite_hpo4_proton /
        std.math.pow(f64, base.dihydrogen_phosphate_dissociation, 3.0);
    const hydroxyapatite_h2po4_hydroxide =
        hydroxyapatite_h2po4_proton *
        std.math.pow(f64, base.water_dissociation, 7.0);

    const result: DerivedConstants = .{
        .calcium_hydrogen_phosphate_proton_equilibrium = hydrogen_phosphate_proton,
        .calcium_hydrogen_phosphate_hydroxide_equilibrium = hydrogen_phosphate_hydroxide,
        .hydroxyapatite_hpo4_proton_equilibrium = hydroxyapatite_hpo4_proton,
        .hydroxyapatite_hpo4_hydroxide_equilibrium = hydroxyapatite_hpo4_hydroxide,
        .hydroxyapatite_h2po4_proton_equilibrium = hydroxyapatite_h2po4_proton,
        .hydroxyapatite_h2po4_hydroxide_equilibrium = hydroxyapatite_h2po4_hydroxide,
    };
    inline for (@typeInfo(DerivedConstants).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCalciumPhosphateEquilibriumConstant;
    return result;
}

const source_defaults: BaseConstants = .{
    .water_dissociation = 1.0e-8,
    .calcium_hydrogen_phosphate_solubility = 1.3e-1,
    .hydroxyapatite_solubility = 4.0e-31,
    .hydrogen_phosphate_dissociation = 4.8e-10,
    .dihydrogen_phosphate_dissociation = 6.2e-5,
};

test "runtime derivation preserves calcium phosphate source powers" {
    const result = try derive(source_defaults);
    try std.testing.expectApproxEqRel(
        result.calcium_hydrogen_phosphate_proton_equilibrium *
            source_defaults.water_dissociation,
        result.calcium_hydrogen_phosphate_hydroxide_equilibrium,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.hydroxyapatite_hpo4_proton_equilibrium /
            std.math.pow(
                f64,
                source_defaults.dihydrogen_phosphate_dissociation,
                3.0,
            ),
        result.hydroxyapatite_h2po4_proton_equilibrium,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.hydroxyapatite_h2po4_proton_equilibrium *
            std.math.pow(f64, source_defaults.water_dissociation, 7.0),
        result.hydroxyapatite_h2po4_hydroxide_equilibrium,
        1.0e-15,
    );
}

test "runtime derivation rejects unsafe calcium phosphate constants" {
    var base = source_defaults;
    base.hydrogen_phosphate_dissociation = 0.0;
    try std.testing.expectError(
        error.InvalidCalciumPhosphateEquilibriumConstant,
        derive(base),
    );

    base = source_defaults;
    base.hydroxyapatite_solubility = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidCalciumPhosphateEquilibriumConstant,
        derive(base),
    );
}
