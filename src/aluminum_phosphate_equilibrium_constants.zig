const std = @import("std");

pub const BaseConstants = struct {
    water_dissociation: f64,
    aluminum_phosphate_solubility: f64,
    hydrogen_phosphate_dissociation: f64,
    dihydrogen_phosphate_dissociation: f64,
    aluminum_first_hydrolysis: f64,
    aluminum_second_hydrolysis: f64,
    aluminum_third_hydrolysis: f64,
    aluminum_fourth_hydrolysis: f64,
};

pub const DerivedConstants = struct {
    hpo4_solid_proton: f64, // SHA0P1
    hpo4_solid_hydroxide: f64, // SYA0P1
    hpo4_first_hydrolysis: f64, // SPA1P1
    hpo4_second_hydrolysis: f64, // SYA2P1
    hpo4_second_proton: f64, // SHA2P1
    hpo4_third_hydrolysis: f64, // SYA3P1
    hpo4_third_proton: f64, // SHA3P1
    hpo4_fourth_hydrolysis: f64, // SYA4P1
    hpo4_fourth_proton: f64, // SHA4P1
    h2po4_solid_proton: f64, // SHA0P2
    h2po4_solid_hydroxide: f64, // SYA0P2
    h2po4_first_hydrolysis: f64, // SYA1P2
    h2po4_first_proton: f64, // SHA1P2
    h2po4_second_hydrolysis: f64, // SPA2P2
    h2po4_third_hydrolysis: f64, // SYA3P2
    h2po4_third_proton: f64, // SHA3P2
    h2po4_fourth_hydrolysis: f64, // SYA4P2
    h2po4_fourth_proton: f64, // SHA4P2
};

/// Exact SOLUTE aluminum-phosphate derivation from solute.f:79-85.
pub fn derive(base: BaseConstants) !DerivedConstants {
    inline for (@typeInfo(BaseConstants).@"struct".fields) |field| {
        const value = @field(base, field.name);
        if (!std.math.isFinite(value) or value <= 0.0)
            return error.InvalidAluminumPhosphateEquilibriumConstant;
    }

    const hpo4_solid_proton =
        base.aluminum_phosphate_solubility /
        base.hydrogen_phosphate_dissociation;
    const hpo4_solid_hydroxide =
        hpo4_solid_proton * base.water_dissociation;
    const hpo4_first =
        hpo4_solid_hydroxide / base.aluminum_first_hydrolysis;
    const hpo4_second =
        hpo4_first / base.aluminum_second_hydrolysis;
    const hpo4_second_proton = hpo4_second * base.water_dissociation;
    const hpo4_third = hpo4_second / base.aluminum_third_hydrolysis;
    const hpo4_third_proton =
        hpo4_third * std.math.pow(f64, base.water_dissociation, 2.0);
    const hpo4_fourth = hpo4_third / base.aluminum_fourth_hydrolysis;
    const hpo4_fourth_proton =
        hpo4_fourth * std.math.pow(f64, base.water_dissociation, 3.0);
    const h2po4_solid_proton =
        hpo4_solid_proton / base.dihydrogen_phosphate_dissociation;
    const h2po4_solid_hydroxide =
        h2po4_solid_proton * std.math.pow(f64, base.water_dissociation, 2.0);
    const h2po4_first =
        h2po4_solid_hydroxide / base.aluminum_first_hydrolysis;
    const h2po4_first_proton = h2po4_first / base.water_dissociation;
    const h2po4_second =
        h2po4_first / base.aluminum_second_hydrolysis;
    const h2po4_third =
        h2po4_second / base.aluminum_third_hydrolysis;
    const h2po4_third_proton = h2po4_third * base.water_dissociation;
    const h2po4_fourth =
        h2po4_third / base.aluminum_fourth_hydrolysis;
    const h2po4_fourth_proton =
        h2po4_fourth * std.math.pow(f64, base.water_dissociation, 2.0);

    const result: DerivedConstants = .{
        .hpo4_solid_proton = hpo4_solid_proton,
        .hpo4_solid_hydroxide = hpo4_solid_hydroxide,
        .hpo4_first_hydrolysis = hpo4_first,
        .hpo4_second_hydrolysis = hpo4_second,
        .hpo4_second_proton = hpo4_second_proton,
        .hpo4_third_hydrolysis = hpo4_third,
        .hpo4_third_proton = hpo4_third_proton,
        .hpo4_fourth_hydrolysis = hpo4_fourth,
        .hpo4_fourth_proton = hpo4_fourth_proton,
        .h2po4_solid_proton = h2po4_solid_proton,
        .h2po4_solid_hydroxide = h2po4_solid_hydroxide,
        .h2po4_first_hydrolysis = h2po4_first,
        .h2po4_first_proton = h2po4_first_proton,
        .h2po4_second_hydrolysis = h2po4_second,
        .h2po4_third_hydrolysis = h2po4_third,
        .h2po4_third_proton = h2po4_third_proton,
        .h2po4_fourth_hydrolysis = h2po4_fourth,
        .h2po4_fourth_proton = h2po4_fourth_proton,
    };
    inline for (@typeInfo(DerivedConstants).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteAluminumPhosphateEquilibriumConstant;
    return result;
}

const source_defaults: BaseConstants = .{
    .water_dissociation = 1.0e-8,
    .aluminum_phosphate_solubility = 9.8e-15,
    .hydrogen_phosphate_dissociation = 4.8e-10,
    .dihydrogen_phosphate_dissociation = 6.2e-5,
    .aluminum_first_hydrolysis = 4.6e-7,
    .aluminum_second_hydrolysis = 7.3e-7,
    .aluminum_third_hydrolysis = 1.8e-5,
    .aluminum_fourth_hydrolysis = 1.2e-5,
};

test "runtime derivation preserves both aluminum phosphate branches" {
    const result = try derive(source_defaults);
    try std.testing.expectApproxEqRel(
        result.hpo4_solid_proton /
            source_defaults.dihydrogen_phosphate_dissociation,
        result.h2po4_solid_proton,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.hpo4_third_hydrolysis *
            std.math.pow(f64, source_defaults.water_dissociation, 2.0),
        result.hpo4_third_proton,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        result.h2po4_fourth_hydrolysis *
            std.math.pow(f64, source_defaults.water_dissociation, 2.0),
        result.h2po4_fourth_proton,
        1.0e-15,
    );
}

test "runtime derivation rejects unsafe aluminum phosphate constants" {
    var base = source_defaults;
    base.aluminum_third_hydrolysis = 0.0;
    try std.testing.expectError(
        error.InvalidAluminumPhosphateEquilibriumConstant,
        derive(base),
    );

    base = source_defaults;
    base.aluminum_phosphate_solubility = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidAluminumPhosphateEquilibriumConstant,
        derive(base),
    );
}
