const std = @import("std");

pub const NutrientUptakeKinetics = struct {
    maximum_uptake_g_m2_h: f64,
    half_saturation_umol_l: f64,
    minimum_concentration_umol_l: f64,
};

pub const RootParameters = struct {
    porosity_fraction: f64,
    ammonium: NutrientUptakeKinetics,
    nitrate: NutrientUptakeKinetics,
    dihydrogen_phosphate: NutrientUptakeKinetics,
    radial_resistivity_m2_mpa_inv_h_inv: f64,
    axial_resistivity_m2_mpa_inv_h_inv: f64,
};

pub const RadiusParameters = struct {
    maximum_primary_radius_m: f64,
    maximum_secondary_radius_m: f64,
};

pub const MycorrhizalParameters = struct {
    maximum_primary_radius_m: f64,
    maximum_secondary_radius_m: f64,
    porosity_fraction: f64,
    ammonium: NutrientUptakeKinetics,
    nitrate: NutrientUptakeKinetics,
    dihydrogen_phosphate: NutrientUptakeKinetics,
    radial_resistivity_m2_mpa_inv_h_inv: f64,
    axial_resistivity_m2_mpa_inv_h_inv: f64,
};

pub const InitializationError = error{
    NonFiniteInput,
    InvalidPorosity,
    InvalidRadius,
    NegativeUptakeParameter,
    NegativeResistivity,
};

/// Translates STARTQ lines 390-403 for one plant species.
pub fn initialize(
    root: RootParameters,
    radii: RadiusParameters,
) InitializationError!MycorrhizalParameters {
    try validateRoot(root);
    if (!std.math.isFinite(radii.maximum_primary_radius_m) or
        !std.math.isFinite(radii.maximum_secondary_radius_m))
    {
        return error.NonFiniteInput;
    }
    if (radii.maximum_primary_radius_m <= 0.0 or
        radii.maximum_secondary_radius_m <= 0.0)
    {
        return error.InvalidRadius;
    }

    // Field order follows STARTQ: radii, porosity, NH4, NO3, H2PO4,
    // radial resistivity, then axial resistivity.
    return .{
        .maximum_primary_radius_m = radii.maximum_primary_radius_m,
        .maximum_secondary_radius_m = radii.maximum_secondary_radius_m,
        .porosity_fraction = root.porosity_fraction,
        .ammonium = root.ammonium,
        .nitrate = root.nitrate,
        .dihydrogen_phosphate = root.dihydrogen_phosphate,
        .radial_resistivity_m2_mpa_inv_h_inv = root.radial_resistivity_m2_mpa_inv_h_inv,
        .axial_resistivity_m2_mpa_inv_h_inv = root.axial_resistivity_m2_mpa_inv_h_inv,
    };
}

fn validateRoot(root: RootParameters) InitializationError!void {
    if (!std.math.isFinite(root.porosity_fraction)) return error.NonFiniteInput;
    if (root.porosity_fraction < 0.0 or root.porosity_fraction >= 1.0) {
        return error.InvalidPorosity;
    }
    inline for (.{ root.ammonium, root.nitrate, root.dihydrogen_phosphate }) |kinetics| {
        inline for (std.meta.fields(NutrientUptakeKinetics)) |field| {
            const value = @field(kinetics, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteInput;
            if (value < 0.0) return error.NegativeUptakeParameter;
        }
    }
    inline for (.{
        root.radial_resistivity_m2_mpa_inv_h_inv,
        root.axial_resistivity_m2_mpa_inv_h_inv,
    }) |resistivity| {
        if (!std.math.isFinite(resistivity)) return error.NonFiniteInput;
        if (resistivity < 0.0) return error.NegativeResistivity;
    }
}

test "mycorrhiza copies root uptake and resistance parameters in STARTQ order" {
    const root = RootParameters{
        .porosity_fraction = 0.35,
        .ammonium = .{
            .maximum_uptake_g_m2_h = 1.0,
            .half_saturation_umol_l = 2.0,
            .minimum_concentration_umol_l = 3.0,
        },
        .nitrate = .{
            .maximum_uptake_g_m2_h = 4.0,
            .half_saturation_umol_l = 5.0,
            .minimum_concentration_umol_l = 6.0,
        },
        .dihydrogen_phosphate = .{
            .maximum_uptake_g_m2_h = 7.0,
            .half_saturation_umol_l = 8.0,
            .minimum_concentration_umol_l = 9.0,
        },
        .radial_resistivity_m2_mpa_inv_h_inv = 10.0,
        .axial_resistivity_m2_mpa_inv_h_inv = 11.0,
    };
    const mycorrhiza = try initialize(root, .{
        .maximum_primary_radius_m = 2.5e-6,
        .maximum_secondary_radius_m = 2.5e-6,
    });

    try std.testing.expectEqual(@as(f64, 2.5e-6), mycorrhiza.maximum_primary_radius_m);
    try std.testing.expectEqual(root.porosity_fraction, mycorrhiza.porosity_fraction);
    try std.testing.expectEqual(root.ammonium, mycorrhiza.ammonium);
    try std.testing.expectEqual(root.nitrate, mycorrhiza.nitrate);
    try std.testing.expectEqual(
        root.dihydrogen_phosphate,
        mycorrhiza.dihydrogen_phosphate,
    );
    try std.testing.expectEqual(
        root.axial_resistivity_m2_mpa_inv_h_inv,
        mycorrhiza.axial_resistivity_m2_mpa_inv_h_inv,
    );
}

test "invalid root porosity fails explicitly" {
    var root = std.mem.zeroes(RootParameters);
    root.porosity_fraction = 1.0;
    try std.testing.expectError(error.InvalidPorosity, initialize(root, .{
        .maximum_primary_radius_m = 2.5e-6,
        .maximum_secondary_radius_m = 2.5e-6,
    }));
}
