const std = @import("std");

pub const ResidueClass = struct {
    water_holding_capacity_m3_g: f64, // THETRX
    carbon_mass_g: f64, // RC0
    dry_bulk_density_megagrams_m3: f64, // BKRS
};

pub const Inputs = struct {
    residue_classes: []const ResidueClass,
    water_volume_m3: f64,
    ice_volume_m3: f64,
    underlying_soil_bulk_density_megagrams_m3: f64,
    organic_carbon_g: f64,
    charcoal_carbon_g: f64,
    horizontal_area_m2: f64,
    volume_threshold_m3: f64,
    positive_density_threshold_megagrams_m3: f64,
};

pub const ActiveProperties = struct {
    wet_volume_m3: f64,
    effective_volume_m3: f64,
    residue_mass_megagrams: f64,
    pore_capacity_m3: f64,
    air_filled_volume_m3: f64,
    porosity_m3_m3: f64,
    water_fraction_m3_m3: f64,
    ice_fraction_m3_m3: f64,
    air_fraction_m3_m3: f64,
    thickness_m: f64,
};

pub const Result = struct {
    water_holding_capacity_m3: f64,
    excess_water_and_ice_m3: f64,
    dry_residue_volume_m3: f64,
    active: ?ActiveProperties,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeInput,
    InvalidResidueBulkDensity,
    InvalidHorizontalArea,
    InvalidPorosity,
    NonFiniteResult,
};

/// Translates HOUR1 lines 4350-4382. Runtime residue classes replace the
/// fixed legacy class selection. Lines 4383 onward use the intentionally
/// replaced legacy retention formulation and are not included.
pub fn calculate(inputs: Inputs) CalculationError!Result {
    inline for (std.meta.fields(Inputs)) |field| {
        if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteInput;
            if (value < 0.0) return error.NegativeInput;
        }
    }
    for (inputs.residue_classes) |class| {
        inline for (std.meta.fields(ResidueClass)) |field| {
            const value = @field(class, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteInput;
            if (value < 0.0) return error.NegativeInput;
        }
        if (class.carbon_mass_g > 0.0 and class.dry_bulk_density_megagrams_m3 <= 0.0) {
            return error.InvalidResidueBulkDensity;
        }
    }

    var water_holding_capacity_m3: f64 = 0.0;
    for (inputs.residue_classes) |class| {
        water_holding_capacity_m3 = water_holding_capacity_m3 +
            class.water_holding_capacity_m3_g * class.carbon_mass_g;
    }
    water_holding_capacity_m3 = @max(0.0, water_holding_capacity_m3);
    const excess_water_and_ice_m3 = @max(
        0.0,
        inputs.water_volume_m3 + inputs.ice_volume_m3 - water_holding_capacity_m3,
    );
    var dry_density_sum: f64 = 0.0;
    for (inputs.residue_classes) |class| {
        if (class.carbon_mass_g > 0.0) {
            dry_density_sum =
                dry_density_sum + class.carbon_mass_g / class.dry_bulk_density_megagrams_m3;
        }
    }
    const dry_residue_volume_m3 = 1.0e-6 * @max(0.0, dry_density_sum);
    const wet_volume_m3 = excess_water_and_ice_m3 + dry_residue_volume_m3;
    if (wet_volume_m3 <= inputs.volume_threshold_m3) {
        return .{
            .water_holding_capacity_m3 = water_holding_capacity_m3,
            .excess_water_and_ice_m3 = excess_water_and_ice_m3,
            .dry_residue_volume_m3 = dry_residue_volume_m3,
            .active = null,
        };
    }
    if (inputs.horizontal_area_m2 <= 0.0) return error.InvalidHorizontalArea;

    const effective_volume_m3 =
        if (inputs.underlying_soil_bulk_density_megagrams_m3 <=
        inputs.positive_density_threshold_megagrams_m3)
            inputs.water_volume_m3 + inputs.ice_volume_m3
        else
            wet_volume_m3;
    const residue_mass_megagrams =
        1.82e-6 * (inputs.organic_carbon_g + inputs.charcoal_carbon_g);
    const pore_capacity_m3 =
        @max(0.0, dry_residue_volume_m3 - residue_mass_megagrams / 1.30);
    const air_filled_volume_m3 = @max(
        0.0,
        pore_capacity_m3 - inputs.water_volume_m3 - inputs.ice_volume_m3,
    );
    var porosity_m3_m3: f64 = 1.0;
    var water_fraction_m3_m3: f64 = 0.0;
    var ice_fraction_m3_m3: f64 = 0.0;
    var air_fraction_m3_m3: f64 = 0.0;
    if (dry_residue_volume_m3 > inputs.volume_threshold_m3) {
        porosity_m3_m3 = pore_capacity_m3 / dry_residue_volume_m3;
        water_fraction_m3_m3 =
            @max(0.0, @min(1.0, inputs.water_volume_m3 / dry_residue_volume_m3));
        ice_fraction_m3_m3 =
            @max(0.0, @min(1.0, inputs.ice_volume_m3 / dry_residue_volume_m3));
        air_fraction_m3_m3 =
            @max(0.0, @min(1.0, air_filled_volume_m3 / dry_residue_volume_m3));
    }
    if (porosity_m3_m3 < 0.0 or porosity_m3_m3 > 1.0) return error.InvalidPorosity;
    const active = ActiveProperties{
        .wet_volume_m3 = wet_volume_m3,
        .effective_volume_m3 = effective_volume_m3,
        .residue_mass_megagrams = residue_mass_megagrams,
        .pore_capacity_m3 = pore_capacity_m3,
        .air_filled_volume_m3 = air_filled_volume_m3,
        .porosity_m3_m3 = porosity_m3_m3,
        .water_fraction_m3_m3 = water_fraction_m3_m3,
        .ice_fraction_m3_m3 = ice_fraction_m3_m3,
        .air_fraction_m3_m3 = air_fraction_m3_m3,
        .thickness_m = wet_volume_m3 / inputs.horizontal_area_m2,
    };
    inline for (std.meta.fields(ActiveProperties)) |field| {
        if (!std.math.isFinite(@field(active, field.name))) return error.NonFiniteResult;
    }
    return .{
        .water_holding_capacity_m3 = water_holding_capacity_m3,
        .excess_water_and_ice_m3 = excess_water_and_ice_m3,
        .dry_residue_volume_m3 = dry_residue_volume_m3,
        .active = active,
    };
}

test "runtime residue classes determine volume porosity and thickness" {
    const classes = [_]ResidueClass{
        .{
            .water_holding_capacity_m3_g = 1.0e-6,
            .carbon_mass_g = 100_000.0,
            .dry_bulk_density_megagrams_m3 = 0.2,
        },
        .{
            .water_holding_capacity_m3_g = 2.0e-6,
            .carbon_mass_g = 50_000.0,
            .dry_bulk_density_megagrams_m3 = 0.1,
        },
    };
    const result = try calculate(.{
        .residue_classes = &classes,
        .water_volume_m3 = 0.3,
        .ice_volume_m3 = 0.0,
        .underlying_soil_bulk_density_megagrams_m3 = 1.2,
        .organic_carbon_g = 100_000.0,
        .charcoal_carbon_g = 0.0,
        .horizontal_area_m2 = 10.0,
        .volume_threshold_m3 = 1.0e-12,
        .positive_density_threshold_megagrams_m3 = 0.0,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.2),
        result.water_holding_capacity_m3,
        1.0e-15,
    );
    try std.testing.expect(result.active.?.thickness_m > 0.0);
    try std.testing.expect(result.active.?.porosity_m3_m3 <= 1.0);
}

test "empty residue returns no active properties" {
    const result = try calculate(.{
        .residue_classes = &.{},
        .water_volume_m3 = 0.0,
        .ice_volume_m3 = 0.0,
        .underlying_soil_bulk_density_megagrams_m3 = 0.0,
        .organic_carbon_g = 0.0,
        .charcoal_carbon_g = 0.0,
        .horizontal_area_m2 = 0.0,
        .volume_threshold_m3 = 0.0,
        .positive_density_threshold_megagrams_m3 = 0.0,
    });
    try std.testing.expect(result.active == null);
}
