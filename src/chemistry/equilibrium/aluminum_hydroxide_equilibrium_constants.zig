const std = @import("std");

pub const BaseConstants = struct {
    water_dissociation: f64,
    aluminum_hydroxide_solubility: f64,
    first_hydrolysis_dissociation: f64,
    second_hydrolysis_dissociation: f64,
    third_hydrolysis_dissociation: f64,
    fourth_hydrolysis_dissociation: f64,
};

pub const DerivedConstants = struct {
    /// Legacy SHALO.
    proton_based_solid_equilibrium: f64,
    /// Legacy SYAL1.
    first_hydrolysis_from_solid: f64,
    /// Legacy SHAL1.
    proton_based_first_hydrolysis: f64,
    /// Legacy SYAL2.
    second_hydrolysis_from_solid: f64,
    /// Legacy SHAL2.
    proton_based_second_hydrolysis: f64,
    /// Legacy SPAL3.
    neutral_hydroxide_equilibrium: f64,
    /// Legacy SYAL4.
    fourth_hydrolysis_from_neutral: f64,
    /// Legacy SHAL4.
    proton_based_fourth_hydrolysis: f64,
};

/// Exact SOLUTE aluminum-hydroxide constant derivation from solute.f:72-75.
///
/// The input constants are runtime scientific parameters. Their dimensions
/// follow the concentration-based equilibrium convention used by SOLUTE, so
/// derived constants do not share one common physical unit.
pub fn derive(base: BaseConstants) !DerivedConstants {
    inline for (@typeInfo(BaseConstants).@"struct".fields) |field| {
        const value = @field(base, field.name);
        if (!std.math.isFinite(value) or value <= 0.0)
            return error.InvalidAluminumEquilibriumConstant;
    }

    const proton_based_solid =
        base.aluminum_hydroxide_solubility /
        std.math.pow(f64, base.water_dissociation, 3.0);
    const first_from_solid =
        base.aluminum_hydroxide_solubility /
        base.first_hydrolysis_dissociation;
    const proton_based_first =
        first_from_solid /
        std.math.pow(f64, base.water_dissociation, 2.0);
    const second_from_solid =
        first_from_solid / base.second_hydrolysis_dissociation;
    const proton_based_second =
        second_from_solid / base.water_dissociation;
    const neutral =
        second_from_solid / base.third_hydrolysis_dissociation;
    const fourth_from_neutral =
        neutral / base.fourth_hydrolysis_dissociation;
    const proton_based_fourth =
        fourth_from_neutral * base.water_dissociation;

    const result: DerivedConstants = .{
        .proton_based_solid_equilibrium = proton_based_solid,
        .first_hydrolysis_from_solid = first_from_solid,
        .proton_based_first_hydrolysis = proton_based_first,
        .second_hydrolysis_from_solid = second_from_solid,
        .proton_based_second_hydrolysis = proton_based_second,
        .neutral_hydroxide_equilibrium = neutral,
        .fourth_hydrolysis_from_neutral = fourth_from_neutral,
        .proton_based_fourth_hydrolysis = proton_based_fourth,
    };
    inline for (@typeInfo(DerivedConstants).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteAluminumEquilibriumConstant;
    return result;
}

const source_defaults: BaseConstants = .{
    .water_dissociation = 1.0e-8,
    .aluminum_hydroxide_solubility = 1.9e-21,
    .first_hydrolysis_dissociation = 4.6e-7,
    .second_hydrolysis_dissociation = 7.3e-7,
    .third_hydrolysis_dissociation = 1.8e-5,
    .fourth_hydrolysis_dissociation = 1.2e-5,
};

test "runtime derivation preserves SOLUTE aluminum dependency chain" {
    const result = try derive(source_defaults);

    try std.testing.expectApproxEqRel(
        source_defaults.aluminum_hydroxide_solubility /
            std.math.pow(f64, source_defaults.water_dissociation, 3.0),
        result.proton_based_solid_equilibrium,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.first_hydrolysis_from_solid /
            source_defaults.second_hydrolysis_dissociation,
        result.second_hydrolysis_from_solid,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.neutral_hydroxide_equilibrium /
            source_defaults.fourth_hydrolysis_dissociation *
            source_defaults.water_dissociation,
        result.proton_based_fourth_hydrolysis,
        1.0e-15,
    );
}

test "runtime derivation rejects zero and non-finite constants" {
    var base = source_defaults;
    base.second_hydrolysis_dissociation = 0.0;
    try std.testing.expectError(
        error.InvalidAluminumEquilibriumConstant,
        derive(base),
    );

    base = source_defaults;
    base.water_dissociation = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidAluminumEquilibriumConstant,
        derive(base),
    );
}
