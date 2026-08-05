const std = @import("std");

pub const Parameters = struct {
    nitrogen_molar_mass_g_per_mol: f64,
    phosphorus_molar_mass_g_per_mol: f64,
};

pub const WaterInputs = struct {
    /// FLQGQ + FLQRQ (m3 step-1).
    rain_m3: f64,
    /// FLQGI + FLQRI (m3 step-1).
    irrigation_m3: f64,
    /// PRECU (m3 h-1).
    subsurface_irrigation_m3_per_h: f64,
    /// XNFH (h step-1).
    timestep_h: f64,
};

pub const NitrogenConcentrations = struct {
    /// CN4R, CN3R, CNOR (mol N m-3).
    rain_nh4_mol_n_per_m3: f64,
    rain_nh3_mol_n_per_m3: f64,
    rain_no3_mol_n_per_m3: f64,
    /// CNNR, CN2R (g N m-3).
    rain_n2_g_n_per_m3: f64,
    rain_n2o_g_n_per_m3: f64,
    /// CN4Q, CN3Q, CNOQ (mol N m-3).
    irrigation_nh4_mol_n_per_m3: f64,
    irrigation_nh3_mol_n_per_m3: f64,
    irrigation_no3_mol_n_per_m3: f64,
    /// CNNQ, CN2Q (g N m-3).
    irrigation_n2_g_n_per_m3: f64,
    irrigation_n2o_g_n_per_m3: f64,
};

pub const PhosphorusConcentrations = struct {
    /// CPOR, CH1PR (mol P m-3).
    rain_h2po4_mol_p_per_m3: f64,
    rain_hpo4_mol_p_per_m3: f64,
    /// CPOQ, CH1PQ (mol P m-3).
    irrigation_h2po4_mol_p_per_m3: f64,
    irrigation_hpo4_mol_p_per_m3: f64,
};

/// Signed gas-N terms in g N step-1. Field order follows REDIST.F 4659--4664.
pub const GasNitrogenFluxes = struct {
    n2_soil_diffusion_g_n: f64,
    n2o_soil_diffusion_g_n: f64,
    nh3_soil_diffusion_g_n: f64,
    nh4_soil_diffusion_g_n: f64,
    n2_mineral_convection_g_n: f64,
    n2o_mineral_convection_g_n: f64,
    nh3_mineral_convection_g_n: f64,
    n2o_bulk_transfer_g_n: f64,
    nh3_bulk_transfer_g_n: f64,
    n2o_surface_dissolution_g_n: f64,
    n2_surface_dissolution_g_n: f64,
    nh3_surface_dissolution_g_n: f64,
    n2_litter_volatilization_g_n: f64,
    n2o_litter_volatilization_g_n: f64,
    nh3_litter_volatilization_g_n: f64,
};

/// Signed surface transformation terms in g element step-1.
pub const SurfaceMineralFluxes = struct {
    /// XNH4S, XNO3S, XNO2S, XN2GS at litter layer zero (g N step-1).
    nh4_g_n: f64,
    no3_g_n: f64,
    no2_g_n: f64,
    n2_g_n: f64,
    /// XH1PS, XH2PS at litter layer zero (g P step-1).
    hpo4_g_p: f64,
    h2po4_g_p: f64,
};

pub const Inputs = struct {
    parameters: Parameters,
    water: WaterInputs,
    nitrogen: NitrogenConcentrations,
    phosphorus: PhosphorusConcentrations,
    gas_nitrogen: GasNitrogenFluxes,
    surface_minerals: SurfaceMineralFluxes,
};

/// Source-weighted molar terms. The legacy source does not document a single
/// physical name for these coefficients, so their SIN/SGN/SIP/SNB/SPB/SNM0/
/// SPM0 traceability is retained explicitly here.
pub const Terms = struct {
    aqueous_nitrogen_weighted_mol: f64,
    gaseous_nitrogen_weighted_mol: f64,
    aqueous_phosphorus_weighted_mol: f64,
    subsurface_nitrogen_weighted_mol: f64,
    subsurface_phosphorus_weighted_mol: f64,
    surface_mineral_nitrogen_weighted_mol: f64,
    surface_mineral_phosphorus_weighted_mol: f64,
};

/// Direct translation of REDIST.F lines 4653--4674.
///
/// All additions retain source order. No algebraic regrouping is performed.
pub fn calculate(inputs: Inputs) !Terms {
    try validate(inputs);
    const water = inputs.water;
    const n = inputs.nitrogen;
    const p = inputs.phosphorus;
    const gas = inputs.gas_nitrogen;
    const mineral = inputs.surface_minerals;
    const n_mass = inputs.parameters.nitrogen_molar_mass_g_per_mol;
    const p_mass = inputs.parameters.phosphorus_molar_mass_g_per_mol;

    // SIN, lines 4653--4656.
    const aqueous_n = water.rain_m3 *
        (2.0 * n.rain_nh4_mol_n_per_m3 + n.rain_nh3_mol_n_per_m3 +
            n.rain_no3_mol_n_per_m3) +
        water.irrigation_m3 *
            (2.0 * n.irrigation_nh4_mol_n_per_m3 +
                n.irrigation_nh3_mol_n_per_m3 +
                n.irrigation_no3_mol_n_per_m3);

    // SGN, lines 4657--4664.
    const gaseous_n = (2.0 * water.rain_m3 *
        (n.rain_n2_g_n_per_m3 + n.rain_n2o_g_n_per_m3) +
        2.0 * water.irrigation_m3 *
            (n.irrigation_n2_g_n_per_m3 + n.irrigation_n2o_g_n_per_m3) +
        2.0 * (gas.n2_soil_diffusion_g_n + gas.n2o_soil_diffusion_g_n) +
        gas.nh3_soil_diffusion_g_n + gas.nh4_soil_diffusion_g_n +
        2.0 * (gas.n2_mineral_convection_g_n + gas.n2o_mineral_convection_g_n) +
        gas.nh3_mineral_convection_g_n + 2.0 * gas.n2o_bulk_transfer_g_n +
        gas.nh3_bulk_transfer_g_n +
        2.0 * (gas.n2o_surface_dissolution_g_n + gas.n2_surface_dissolution_g_n) +
        gas.nh3_surface_dissolution_g_n +
        2.0 * (gas.n2_litter_volatilization_g_n + gas.n2o_litter_volatilization_g_n) +
        gas.nh3_litter_volatilization_g_n) / n_mass;

    // SIP, lines 4665--4668.
    const aqueous_p = water.rain_m3 *
        (3.0 * p.rain_h2po4_mol_p_per_m3 +
            2.0 * p.rain_hpo4_mol_p_per_m3) +
        water.irrigation_m3 *
            (3.0 * p.irrigation_h2po4_mol_p_per_m3 +
                2.0 * p.irrigation_hpo4_mol_p_per_m3);

    // SNB and SPB, lines 4669--4671.
    const subsurface_n = (-water.subsurface_irrigation_m3_per_h *
        (n.irrigation_n2_g_n_per_m3 + n.irrigation_n2o_g_n_per_m3) -
        water.subsurface_irrigation_m3_per_h *
            (2.0 * n.irrigation_nh4_mol_n_per_m3 +
                n.irrigation_nh3_mol_n_per_m3 +
                n.irrigation_no3_mol_n_per_m3)) * water.timestep_h;
    const subsurface_p = (-water.subsurface_irrigation_m3_per_h *
        (3.0 * p.irrigation_h2po4_mol_p_per_m3 +
            2.0 * p.irrigation_hpo4_mol_p_per_m3)) * water.timestep_h;

    // SNM0 and SPM0, lines 4672--4674.
    const mineral_n = (2.0 * mineral.nh4_g_n + mineral.no3_g_n +
        mineral.no2_g_n - 2.0 * mineral.n2_g_n) / n_mass;
    const mineral_p = (2.0 * mineral.hpo4_g_p +
        3.0 * mineral.h2po4_g_p) / p_mass;

    const result = Terms{
        .aqueous_nitrogen_weighted_mol = aqueous_n,
        .gaseous_nitrogen_weighted_mol = gaseous_n,
        .aqueous_phosphorus_weighted_mol = aqueous_p,
        .subsurface_nitrogen_weighted_mol = subsurface_n,
        .subsurface_phosphorus_weighted_mol = subsurface_p,
        .surface_mineral_nitrogen_weighted_mol = mineral_n,
        .surface_mineral_phosphorus_weighted_mol = mineral_p,
    };
    inline for (@typeInfo(Terms).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceNutrientStoichiometricTerm;
    }
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(inputs.parameters, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidNutrientMolarMass;
    }
    inline for (@typeInfo(WaterInputs).@"struct".fields) |field| {
        const value = @field(inputs.water, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceNutrientInput;
        if (value < 0) return error.NegativeSurfaceNutrientWaterInput;
    }
    if (inputs.water.timestep_h == 0) return error.InvalidSurfaceNutrientTimestep;
    inline for (.{ NitrogenConcentrations, PhosphorusConcentrations }) |T| {
        const values = if (T == NitrogenConcentrations) inputs.nitrogen else inputs.phosphorus;
        inline for (@typeInfo(T).@"struct".fields) |field| {
            const value = @field(values, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteSurfaceNutrientInput;
            if (value < 0) return error.NegativeSurfaceNutrientConcentration;
        }
    }
    inline for (.{ GasNitrogenFluxes, SurfaceMineralFluxes }) |T| {
        const values = if (T == GasNitrogenFluxes) inputs.gas_nitrogen else inputs.surface_minerals;
        inline for (@typeInfo(T).@"struct".fields) |field| {
            if (!std.math.isFinite(@field(values, field.name)))
                return error.NonFiniteSurfaceNutrientInput;
        }
    }
}

fn testInputs() Inputs {
    return .{
        .parameters = .{
            .nitrogen_molar_mass_g_per_mol = 14,
            .phosphorus_molar_mass_g_per_mol = 31,
        },
        .water = .{
            .rain_m3 = 2,
            .irrigation_m3 = 3,
            .subsurface_irrigation_m3_per_h = 0.5,
            .timestep_h = 0.25,
        },
        .nitrogen = .{
            .rain_nh4_mol_n_per_m3 = 1,
            .rain_nh3_mol_n_per_m3 = 2,
            .rain_no3_mol_n_per_m3 = 3,
            .rain_n2_g_n_per_m3 = 4,
            .rain_n2o_g_n_per_m3 = 5,
            .irrigation_nh4_mol_n_per_m3 = 6,
            .irrigation_nh3_mol_n_per_m3 = 7,
            .irrigation_no3_mol_n_per_m3 = 8,
            .irrigation_n2_g_n_per_m3 = 9,
            .irrigation_n2o_g_n_per_m3 = 10,
        },
        .phosphorus = .{
            .rain_h2po4_mol_p_per_m3 = 1,
            .rain_hpo4_mol_p_per_m3 = 2,
            .irrigation_h2po4_mol_p_per_m3 = 3,
            .irrigation_hpo4_mol_p_per_m3 = 4,
        },
        .gas_nitrogen = std.mem.zeroes(GasNitrogenFluxes),
        .surface_minerals = .{
            .nh4_g_n = 14,
            .no3_g_n = 28,
            .no2_g_n = 42,
            .n2_g_n = 7,
            .hpo4_g_p = 31,
            .h2po4_g_p = 62,
        },
    };
}

test "REDIST nutrient terms preserve aqueous and subsurface source equations" {
    const result = try calculate(testInputs());
    try std.testing.expectEqual(
        @as(f64, 2 * (2 * 1 + 2 + 3) + 3 * (2 * 6 + 7 + 8)),
        result.aqueous_nitrogen_weighted_mol,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * (3 * 1 + 2 * 2) + 3 * (3 * 3 + 2 * 4)),
        result.aqueous_phosphorus_weighted_mol,
    );
    try std.testing.expectEqual(
        @as(f64, (-0.5 * (9 + 10) - 0.5 * (2 * 6 + 7 + 8)) * 0.25),
        result.subsurface_nitrogen_weighted_mol,
    );
    try std.testing.expectEqual(
        @as(f64, (-0.5 * (3 * 3 + 2 * 4)) * 0.25),
        result.subsurface_phosphorus_weighted_mol,
    );
}

test "REDIST nutrient terms preserve surface mineral signs and molar conversion" {
    const result = try calculate(testInputs());
    try std.testing.expectEqual(
        @as(f64, (2 * 14 + 28 + 42 - 2 * 7) / 14),
        result.surface_mineral_nitrogen_weighted_mol,
    );
    try std.testing.expectEqual(
        @as(f64, (2 * 31 + 3 * 62) / 31),
        result.surface_mineral_phosphorus_weighted_mol,
    );
}

test "REDIST gaseous nitrogen term preserves every source coefficient" {
    var inputs = testInputs();
    inline for (@typeInfo(GasNitrogenFluxes).@"struct".fields, 1..) |field, value| {
        @field(inputs.gas_nitrogen, field.name) = @floatFromInt(value);
    }
    const result = try calculate(inputs);
    const expected = (2.0 * 2.0 * (4.0 + 5.0) +
        2.0 * 3.0 * (9.0 + 10.0) +
        2.0 * (1.0 + 2.0) + 3.0 + 4.0 +
        2.0 * (5.0 + 6.0) + 7.0 + 2.0 * 8.0 + 9.0 +
        2.0 * (10.0 + 11.0) + 12.0 +
        2.0 * (13.0 + 14.0) + 15.0) / 14.0;
    try std.testing.expectEqual(expected, result.gaseous_nitrogen_weighted_mol);
}

test "REDIST nutrient terms reject negative aqueous concentrations" {
    var inputs = testInputs();
    inputs.nitrogen.rain_n2o_g_n_per_m3 = -1;
    try std.testing.expectError(
        error.NegativeSurfaceNutrientConcentration,
        calculate(inputs),
    );
}
