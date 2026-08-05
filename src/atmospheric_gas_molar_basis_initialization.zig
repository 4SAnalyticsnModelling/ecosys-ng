const std = @import("std");

pub const MassConcentrations = struct {
    carbon_dioxide_g_c_per_m3: f64, // CCO2EI
    methane_g_c_per_m3: f64, // CCH4E
    oxygen_g_o2_per_m3: f64, // COXYE
    dinitrogen_g_n_per_m3: f64, // CZ2GE
    nitrous_oxide_g_n_per_m3: f64, // CZ2OE
};

pub const MolarMasses = struct {
    carbon_g_per_mol_c: f64, // source 12
    oxygen_g_per_mol_o2: f64, // source 32
    nitrogen_g_per_mol_n: f64, // source 14; N2 and N2O are N-basis
};

pub const MolarConcentrations = struct {
    carbon_dioxide_mol_c_per_m3: f64, // CCO2M
    methane_mol_c_per_m3: f64, // CCH4M
    oxygen_mol_o2_per_m3: f64, // COXYM
    dinitrogen_mol_n_per_m3: f64, // CZ2GM
    nitrous_oxide_mol_n_per_m3: f64, // CZ2OM
};

/// Direct translation of STARTE 112--116 for one runtime grid cell.
pub fn convert(mass: MassConcentrations, molar_mass: MolarMasses) !MolarConcentrations {
    inline for (@typeInfo(MassConcentrations).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidAtmosphericGasMassConcentration;
    }
    inline for (@typeInfo(MolarMasses).@"struct".fields) |field| {
        const value = @field(molar_mass, field.name);
        if (!std.math.isFinite(value) or value <= 0) return error.InvalidAtmosphericGasMolarMass;
    }
    const result: MolarConcentrations = .{
        .carbon_dioxide_mol_c_per_m3 = mass.carbon_dioxide_g_c_per_m3 / molar_mass.carbon_g_per_mol_c,
        .methane_mol_c_per_m3 = mass.methane_g_c_per_m3 / molar_mass.carbon_g_per_mol_c,
        .oxygen_mol_o2_per_m3 = mass.oxygen_g_o2_per_m3 / molar_mass.oxygen_g_per_mol_o2,
        .dinitrogen_mol_n_per_m3 = mass.dinitrogen_g_n_per_m3 / molar_mass.nitrogen_g_per_mol_n,
        .nitrous_oxide_mol_n_per_m3 = mass.nitrous_oxide_g_n_per_m3 / molar_mass.nitrogen_g_per_mol_n,
    };
    inline for (@typeInfo(MolarConcentrations).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteAtmosphericGasMolarConcentration;
    return result;
}

test "STARTE atmospheric gases preserve elemental molar basis and source order" {
    const result = try convert(.{
        .carbon_dioxide_g_c_per_m3 = 24,
        .methane_g_c_per_m3 = 6,
        .oxygen_g_o2_per_m3 = 64,
        .dinitrogen_g_n_per_m3 = 42,
        .nitrous_oxide_g_n_per_m3 = 7,
    }, .{
        .carbon_g_per_mol_c = 12,
        .oxygen_g_per_mol_o2 = 32,
        .nitrogen_g_per_mol_n = 14,
    });
    try std.testing.expectEqual(@as(f64, 2), result.carbon_dioxide_mol_c_per_m3);
    try std.testing.expectEqual(@as(f64, 0.5), result.methane_mol_c_per_m3);
    try std.testing.expectEqual(@as(f64, 2), result.oxygen_mol_o2_per_m3);
    try std.testing.expectEqual(@as(f64, 3), result.dinitrogen_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 0.5), result.nitrous_oxide_mol_n_per_m3);
}

test "STARTE atmospheric conversion uses runtime molar masses" {
    const result = try convert(.{
        .carbon_dioxide_g_c_per_m3 = 20,
        .methane_g_c_per_m3 = 10,
        .oxygen_g_o2_per_m3 = 30,
        .dinitrogen_g_n_per_m3 = 15,
        .nitrous_oxide_g_n_per_m3 = 5,
    }, .{
        .carbon_g_per_mol_c = 10,
        .oxygen_g_per_mol_o2 = 30,
        .nitrogen_g_per_mol_n = 5,
    });
    try std.testing.expectEqual(@as(f64, 2), result.carbon_dioxide_mol_c_per_m3);
    try std.testing.expectEqual(@as(f64, 3), result.dinitrogen_mol_n_per_m3);
}

test "STARTE atmospheric conversion rejects a late invalid input" {
    try std.testing.expectError(error.InvalidAtmosphericGasMassConcentration, convert(.{
        .carbon_dioxide_g_c_per_m3 = 1,
        .methane_g_c_per_m3 = 1,
        .oxygen_g_o2_per_m3 = 1,
        .dinitrogen_g_n_per_m3 = 1,
        .nitrous_oxide_g_n_per_m3 = std.math.nan(f64),
    }, .{
        .carbon_g_per_mol_c = 12,
        .oxygen_g_per_mol_o2 = 32,
        .nitrogen_g_per_mol_n = 14,
    }));
}
