const std = @import("std");

pub const MolarTotals = struct {
    carbon_dioxide_mol_c_per_step: f64,
    ammonium_mol_n_per_step: f64,
    ammonia_mol_n_per_step: f64,
    nitrate_mol_n_per_step: f64,
    hydrogen_phosphate_mol_p_per_step: f64,
    dihydrogen_phosphate_mol_p_per_step: f64,
};

pub const ElementalMolarMasses = struct {
    carbon_g_c_per_mol: f64,
    nitrogen_g_n_per_mol: f64,
    phosphorus_g_p_per_mol: f64,
};

pub const MassTotals = struct {
    carbon_dioxide_g_c_per_step: f64,
    ammonium_g_n_per_step: f64,
    ammonia_g_n_per_step: f64,
    nitrate_g_n_per_step: f64,
    hydrogen_phosphate_g_p_per_step: f64,
    dihydrogen_phosphate_g_p_per_step: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5222--5227.
///
/// Signed molar transformations are converted to elemental mass using runtime
/// molar masses. No spatial state is read or mutated by this pure kernel.
pub fn calculateSourceOrder(
    molar: MolarTotals,
    molar_mass: ElementalMolarMasses,
) !MassTotals {
    try validateInputs(molar, molar_mass);

    // SOLUTE.F 5222--5227, preserved in assignment order.
    const carbon_dioxide =
        molar.carbon_dioxide_mol_c_per_step *
        molar_mass.carbon_g_c_per_mol;
    const ammonium =
        molar.ammonium_mol_n_per_step *
        molar_mass.nitrogen_g_n_per_mol;
    const ammonia =
        molar.ammonia_mol_n_per_step *
        molar_mass.nitrogen_g_n_per_mol;
    const nitrate =
        molar.nitrate_mol_n_per_step *
        molar_mass.nitrogen_g_n_per_mol;
    const hydrogen_phosphate =
        molar.hydrogen_phosphate_mol_p_per_step *
        molar_mass.phosphorus_g_p_per_mol;
    const dihydrogen_phosphate =
        molar.dihydrogen_phosphate_mol_p_per_step *
        molar_mass.phosphorus_g_p_per_mol;

    const result: MassTotals = .{
        .carbon_dioxide_g_c_per_step = carbon_dioxide,
        .ammonium_g_n_per_step = ammonium,
        .ammonia_g_n_per_step = ammonia,
        .nitrate_g_n_per_step = nitrate,
        .hydrogen_phosphate_g_p_per_step = hydrogen_phosphate,
        .dihydrogen_phosphate_g_p_per_step = dihydrogen_phosphate,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(
    molar: MolarTotals,
    molar_mass: ElementalMolarMasses,
) !void {
    inline for (@typeInfo(MolarTotals).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(molar, field.name)))
            return error.InvalidSurfaceLitterElementalMassConversionInput;
    }
    inline for (@typeInfo(ElementalMolarMasses).@"struct".fields) |field| {
        const value = @field(molar_mass, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSurfaceLitterElementalMassConversionInput;
    }
}

fn validateResult(result: MassTotals) !void {
    inline for (@typeInfo(MassTotals).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterElementalMassConversionResult;
    }
}

fn testMolarTotals() MolarTotals {
    return .{
        .carbon_dioxide_mol_c_per_step = 0.5,
        .ammonium_mol_n_per_step = 0.6,
        .ammonia_mol_n_per_step = 0.7,
        .nitrate_mol_n_per_step = 0.8,
        .hydrogen_phosphate_mol_p_per_step = 0.9,
        .dihydrogen_phosphate_mol_p_per_step = 1.0,
    };
}

fn testMolarMasses() ElementalMolarMasses {
    return .{
        .carbon_g_c_per_mol = 12,
        .nitrogen_g_n_per_mol = 14,
        .phosphorus_g_p_per_mol = 31,
    };
}

test "SOLUTE elemental conversion preserves every source expression" {
    const molar = testMolarTotals();
    const masses = testMolarMasses();
    const result = try calculateSourceOrder(molar, masses);

    try std.testing.expectEqual(
        molar.carbon_dioxide_mol_c_per_step * masses.carbon_g_c_per_mol,
        result.carbon_dioxide_g_c_per_step,
    );
    try std.testing.expectEqual(
        molar.ammonium_mol_n_per_step * masses.nitrogen_g_n_per_mol,
        result.ammonium_g_n_per_step,
    );
    try std.testing.expectEqual(
        molar.ammonia_mol_n_per_step * masses.nitrogen_g_n_per_mol,
        result.ammonia_g_n_per_step,
    );
    try std.testing.expectEqual(
        molar.nitrate_mol_n_per_step * masses.nitrogen_g_n_per_mol,
        result.nitrate_g_n_per_step,
    );
    try std.testing.expectEqual(
        molar.hydrogen_phosphate_mol_p_per_step *
            masses.phosphorus_g_p_per_mol,
        result.hydrogen_phosphate_g_p_per_step,
    );
    try std.testing.expectEqual(
        molar.dihydrogen_phosphate_mol_p_per_step *
            masses.phosphorus_g_p_per_mol,
        result.dihydrogen_phosphate_g_p_per_step,
    );
}

test "elemental conversion preserves moles under inverse conversion" {
    const molar = testMolarTotals();
    const masses = testMolarMasses();
    const result = try calculateSourceOrder(molar, masses);

    try std.testing.expectApproxEqAbs(
        molar.carbon_dioxide_mol_c_per_step,
        result.carbon_dioxide_g_c_per_step / masses.carbon_g_c_per_mol,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        molar.ammonium_mol_n_per_step,
        result.ammonium_g_n_per_step / masses.nitrogen_g_n_per_mol,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        molar.hydrogen_phosphate_mol_p_per_step,
        result.hydrogen_phosphate_g_p_per_step /
            masses.phosphorus_g_p_per_mol,
        1.0e-16,
    );
}

test "elemental conversion preserves signed transformations" {
    var molar = testMolarTotals();
    inline for (@typeInfo(MolarTotals).@"struct".fields) |field|
        @field(molar, field.name) *= -1;
    const result = try calculateSourceOrder(molar, testMolarMasses());

    inline for (@typeInfo(MassTotals).@"struct".fields) |field|
        try std.testing.expect(@field(result, field.name) < 0);
}

test "elemental conversion rejects invalid input and overflow" {
    var molar = testMolarTotals();
    var masses = testMolarMasses();
    masses.nitrogen_g_n_per_mol = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterElementalMassConversionInput,
        calculateSourceOrder(molar, masses),
    );

    molar = testMolarTotals();
    masses = testMolarMasses();
    molar.ammonium_mol_n_per_step = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterElementalMassConversionInput,
        calculateSourceOrder(molar, masses),
    );

    molar = testMolarTotals();
    masses = testMolarMasses();
    molar.carbon_dioxide_mol_c_per_step = std.math.floatMax(f64);
    masses.carbon_g_c_per_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterElementalMassConversionResult,
        calculateSourceOrder(molar, masses),
    );
}
