const std = @import("std");

pub const Properties = struct {
    field_capacity_m3_m3: f64,
    wilting_point_m3_m3: f64,
    cation_exchange_capacity_mol_mg: f64,
    anion_exchange_capacity_mol_mg: f64,
};

pub const AdjustmentError = error{
    NonFiniteInput,
    NegativeInput,
    InvalidWaterRetention,
    NonFiniteResult,
};

/// Translates HOUR1 lines 3720-3725. `charcoal_carbon_g` is legacy DORGCC
/// and `effective_soil_volume_m3` is VOLY.
pub fn adjust(
    properties: Properties,
    charcoal_carbon_g: f64,
    effective_soil_volume_m3: f64,
    volume_threshold_m3: f64,
) AdjustmentError!Properties {
    inline for (std.meta.fields(Properties)) |field| {
        const value = @field(properties, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeInput;
    }
    const additional = [_]f64{
        charcoal_carbon_g,
        effective_soil_volume_m3,
        volume_threshold_m3,
    };
    for (additional) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeInput;
    }
    if (properties.field_capacity_m3_m3 > 1.0 or properties.wilting_point_m3_m3 > 1.0) {
        return error.InvalidWaterRetention;
    }
    if (effective_soil_volume_m3 <= volume_threshold_m3) return properties;

    const adjusted = Properties{
        .field_capacity_m3_m3 = properties.field_capacity_m3_m3 +
            1.0e-6 * charcoal_carbon_g / effective_soil_volume_m3,
        .wilting_point_m3_m3 = properties.wilting_point_m3_m3 +
            1.0e-6 * charcoal_carbon_g / effective_soil_volume_m3,
        .cation_exchange_capacity_mol_mg = properties.cation_exchange_capacity_mol_mg +
            1.0e-3 * charcoal_carbon_g / effective_soil_volume_m3,
        .anion_exchange_capacity_mol_mg = properties.anion_exchange_capacity_mol_mg +
            1.0e-3 * charcoal_carbon_g / effective_soil_volume_m3,
    };
    inline for (std.meta.fields(Properties)) |field| {
        if (!std.math.isFinite(@field(adjusted, field.name))) return error.NonFiniteResult;
    }
    if (adjusted.field_capacity_m3_m3 > 1.0 or
        adjusted.wilting_point_m3_m3 > adjusted.field_capacity_m3_m3)
    {
        return error.InvalidWaterRetention;
    }
    return adjusted;
}

test "charcoal increments retention and exchange capacities in source order" {
    const adjusted = try adjust(.{
        .field_capacity_m3_m3 = 0.30,
        .wilting_point_m3_m3 = 0.10,
        .cation_exchange_capacity_mol_mg = 20.0,
        .anion_exchange_capacity_mol_mg = 2.0,
    }, 100_000.0, 10.0, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0.31), adjusted.field_capacity_m3_m3);
    try std.testing.expectEqual(@as(f64, 0.11), adjusted.wilting_point_m3_m3);
    try std.testing.expectEqual(@as(f64, 30.0), adjusted.cation_exchange_capacity_mol_mg);
    try std.testing.expectEqual(@as(f64, 12.0), adjusted.anion_exchange_capacity_mol_mg);
}

test "volume at threshold preserves existing properties" {
    const properties = Properties{
        .field_capacity_m3_m3 = 0.30,
        .wilting_point_m3_m3 = 0.10,
        .cation_exchange_capacity_mol_mg = 20.0,
        .anion_exchange_capacity_mol_mg = 2.0,
    };
    try std.testing.expectEqual(properties, try adjust(properties, 100.0, 0.0, 0.0));
}
