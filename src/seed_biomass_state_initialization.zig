const std = @import("std");

pub const ElementalMass = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Inputs = struct {
    grain_nitrogen_to_carbon_g_n_g_c: f64,
    grain_phosphorus_to_carbon_g_p_g_c: f64,
    initial_leaf_carbon_g_c: f64,
    initial_sheath_carbon_g_c: f64,
    preceding_leaf_sheath_carbon_g_c: f64,
    stalk_carbon_g_c: f64,
    standing_stalk_area_m2: f64,
    total_canopy_water_potential_mpa: f64,
    canopy_nonstructural_carbon_g_c: f64,
    primary_root_layer_carbon_g_c: f64,
    total_primary_root_carbon_g_c: f64,
    maximum_root_protein_concentration_g_g_c: f64,
    root_nonstructural_carbon_g_c: f64,
};

pub const Parameters = struct {
    sapwood_thickness_m: f64,
    stalk_volume_per_carbon_m3_g_c: f64,
    fully_hydrated_dry_matter_fraction: f64,
    water_potential_dry_matter_numerator_mpa_inv: f64,
    water_potential_dry_matter_denominator_mpa_inv: f64,
    water_potential_dry_matter_denominator_offset: f64,
    cubic_metres_per_cubic_centimetre: f64,
};

pub const Result = struct {
    leaf_mass: ElementalMass,
    branch_leaf_sheath_carbon_g_c: f64,
    total_leaf_sheath_carbon_g_c: f64,
    water_equivalent_plant_carbon_g_c: f64,
    dry_matter_fraction: f64,
    internal_canopy_water_volume_m3: f64,
    surface_canopy_water_volume_m3: f64,
    canopy_nonstructural_pool: ElementalMass,
    primary_root_layer_mass: ElementalMass,
    total_primary_root_mass: ElementalMass,
    total_root_carbon_g_c: f64,
    structural_root_carbon_g_c: f64,
    root_protein_g: f64,
    root_nonstructural_pool: ElementalMass,
};

pub const InitializationError = error{
    NonFiniteInput,
    NegativeMassOrRatio,
    InvalidParameter,
    NonPositiveDryMatterFraction,
    NonFiniteResult,
};

/// Translates STARTQ lines 859-885 for fresh seed morphology and biomass.
pub fn initialize(
    inputs: Inputs,
    parameters: Parameters,
) InitializationError!Result {
    try validate(inputs, parameters);

    const leaf_mass = ElementalMass{
        .carbon_g_c = inputs.initial_leaf_carbon_g_c,
        .nitrogen_g_n = inputs.grain_nitrogen_to_carbon_g_n_g_c *
            inputs.initial_leaf_carbon_g_c,
        .phosphorus_g_p = inputs.grain_phosphorus_to_carbon_g_p_g_c *
            inputs.initial_leaf_carbon_g_c,
    };
    const branch_leaf_sheath_carbon_g_c =
        inputs.initial_leaf_carbon_g_c + inputs.initial_sheath_carbon_g_c;
    const total_leaf_sheath_carbon_g_c =
        inputs.preceding_leaf_sheath_carbon_g_c + branch_leaf_sheath_carbon_g_c;
    const sapwood_equivalent_carbon_g_c =
        parameters.sapwood_thickness_m * inputs.standing_stalk_area_m2 /
        parameters.stalk_volume_per_carbon_m3_g_c;
    const water_equivalent_plant_carbon_g_c = @max(
        0.0,
        total_leaf_sheath_carbon_g_c +
            @min(inputs.stalk_carbon_g_c, sapwood_equivalent_carbon_g_c),
    );
    const absolute_water_potential_mpa = @abs(inputs.total_canopy_water_potential_mpa);
    const dry_matter_fraction = parameters.fully_hydrated_dry_matter_fraction +
        parameters.water_potential_dry_matter_numerator_mpa_inv *
            absolute_water_potential_mpa /
            (parameters.water_potential_dry_matter_denominator_mpa_inv *
                absolute_water_potential_mpa +
                parameters.water_potential_dry_matter_denominator_offset);
    if (dry_matter_fraction <= 0.0) return error.NonPositiveDryMatterFraction;
    const internal_canopy_water_volume_m3 =
        parameters.cubic_metres_per_cubic_centimetre *
        water_equivalent_plant_carbon_g_c / dry_matter_fraction;

    const result = Result{
        .leaf_mass = leaf_mass,
        .branch_leaf_sheath_carbon_g_c = branch_leaf_sheath_carbon_g_c,
        .total_leaf_sheath_carbon_g_c = total_leaf_sheath_carbon_g_c,
        .water_equivalent_plant_carbon_g_c = water_equivalent_plant_carbon_g_c,
        .dry_matter_fraction = dry_matter_fraction,
        .internal_canopy_water_volume_m3 = internal_canopy_water_volume_m3,
        .surface_canopy_water_volume_m3 = 0.0,
        .canopy_nonstructural_pool = elementalFromCarbon(
            inputs.canopy_nonstructural_carbon_g_c,
            inputs,
        ),
        .primary_root_layer_mass = elementalFromCarbon(
            inputs.primary_root_layer_carbon_g_c,
            inputs,
        ),
        .total_primary_root_mass = elementalFromCarbon(
            inputs.total_primary_root_carbon_g_c,
            inputs,
        ),
        .total_root_carbon_g_c = inputs.primary_root_layer_carbon_g_c,
        .structural_root_carbon_g_c = inputs.primary_root_layer_carbon_g_c,
        .root_protein_g = inputs.primary_root_layer_carbon_g_c *
            inputs.maximum_root_protein_concentration_g_g_c,
        .root_nonstructural_pool = elementalFromCarbon(
            inputs.root_nonstructural_carbon_g_c,
            inputs,
        ),
    };
    inline for (std.meta.fields(Result)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(result, field.name))) {
            return error.NonFiniteResult;
        }
    }
    inline for (.{
        result.leaf_mass,
        result.canopy_nonstructural_pool,
        result.primary_root_layer_mass,
        result.total_primary_root_mass,
        result.root_nonstructural_pool,
    }) |mass| {
        inline for (std.meta.fields(ElementalMass)) |field| {
            if (!std.math.isFinite(@field(mass, field.name))) return error.NonFiniteResult;
        }
    }
    return result;
}

fn elementalFromCarbon(carbon_g_c: f64, inputs: Inputs) ElementalMass {
    return .{
        .carbon_g_c = carbon_g_c,
        .nitrogen_g_n = inputs.grain_nitrogen_to_carbon_g_n_g_c * carbon_g_c,
        .phosphorus_g_p = inputs.grain_phosphorus_to_carbon_g_p_g_c * carbon_g_c,
    };
}

fn validate(inputs: Inputs, parameters: Parameters) InitializationError!void {
    inline for (std.meta.fields(Inputs)) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (!std.mem.eql(u8, field.name, "total_canopy_water_potential_mpa") and
            value < 0.0)
        {
            return error.NegativeMassOrRatio;
        }
    }
    inline for (std.meta.fields(Parameters)) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value <= 0.0) return error.InvalidParameter;
    }
}

test "seed morphology preserves STARTQ C N P and water equations" {
    const inputs = Inputs{
        .grain_nitrogen_to_carbon_g_n_g_c = 0.04,
        .grain_phosphorus_to_carbon_g_p_g_c = 0.005,
        .initial_leaf_carbon_g_c = 2.0,
        .initial_sheath_carbon_g_c = 1.0,
        .preceding_leaf_sheath_carbon_g_c = 0.5,
        .stalk_carbon_g_c = 4.0,
        .standing_stalk_area_m2 = 10.0,
        .total_canopy_water_potential_mpa = -0.001,
        .canopy_nonstructural_carbon_g_c = 3.0,
        .primary_root_layer_carbon_g_c = 5.0,
        .total_primary_root_carbon_g_c = 5.0,
        .maximum_root_protein_concentration_g_g_c = 0.2,
        .root_nonstructural_carbon_g_c = 1.5,
    };
    const parameters = Parameters{
        .sapwood_thickness_m = 0.0025,
        .stalk_volume_per_carbon_m3_g_c = 4.0e-6,
        .fully_hydrated_dry_matter_fraction = 0.16,
        .water_potential_dry_matter_numerator_mpa_inv = 0.10,
        .water_potential_dry_matter_denominator_mpa_inv = 0.05,
        .water_potential_dry_matter_denominator_offset = 2.0,
        .cubic_metres_per_cubic_centimetre = 1.0e-6,
    };
    const result = try initialize(inputs, parameters);

    try std.testing.expectEqual(@as(f64, 0.08), result.leaf_mass.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 3.5), result.total_leaf_sheath_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 7.5), result.water_equivalent_plant_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 1.0), result.root_protein_g);
    try std.testing.expectEqual(@as(f64, 0.06), result.root_nonstructural_pool.nitrogen_g_n);
    try std.testing.expect(result.internal_canopy_water_volume_m3 > 0.0);
}

test "zero stalk volume coefficient fails explicitly" {
    var parameters = std.mem.zeroes(Parameters);
    parameters.sapwood_thickness_m = 0.0025;
    try std.testing.expectError(
        error.InvalidParameter,
        initialize(std.mem.zeroes(Inputs), parameters),
    );
}
