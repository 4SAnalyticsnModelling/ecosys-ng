const std = @import("std");

pub const AllocationError = error{
    NonFiniteInput,
    InvalidClayMassFraction,
    InvalidOrganicCarbonConcentration,
    InvalidAllocationFraction,
};

/// Translates `hour1.f` lines 2910--2911 (EHUM).
///
/// `surface_clay_megagrams_megagrams` is Mg clay per Mg soil,
/// `soil_organic_carbon_g_per_megagram` is g C per Mg soil, and the result is the
/// dimensionless fraction of microbial decomposition product sent to humus.
pub fn humusAllocationFraction(
    surface_clay_megagrams_megagrams: f64,
    soil_organic_carbon_g_per_megagram: f64,
) AllocationError!f64 {
    if (!std.math.isFinite(surface_clay_megagrams_megagrams) or
        !std.math.isFinite(soil_organic_carbon_g_per_megagram))
    {
        return error.NonFiniteInput;
    }
    if (surface_clay_megagrams_megagrams < 0.0 or surface_clay_megagrams_megagrams > 1.0) {
        return error.InvalidClayMassFraction;
    }
    if (soil_organic_carbon_g_per_megagram < 0.0) {
        return error.InvalidOrganicCarbonConcentration;
    }

    const allocation_fraction = 0.150 +
        0.300 * @min(0.333, surface_clay_megagrams_megagrams) +
        0.182e-6 * soil_organic_carbon_g_per_megagram;
    if (!std.math.isFinite(allocation_fraction) or
        allocation_fraction < 0.0 or allocation_fraction > 1.0)
    {
        return error.InvalidAllocationFraction;
    }
    return allocation_fraction;
}

test "humus allocation preserves legacy operation order below clay cap" {
    const fraction = try humusAllocationFraction(0.2, 100_000.0);
    const expected = 0.150 + 0.300 * 0.2 + 0.182e-6 * 100_000.0;
    try std.testing.expectEqual(expected, fraction);
}

test "clay contribution is capped at legacy mass fraction" {
    const at_cap = try humusAllocationFraction(0.333, 0.0);
    const above_cap = try humusAllocationFraction(0.8, 0.0);
    try std.testing.expectEqual(at_cap, above_cap);
}

test "allocation exceeding a physical fraction fails explicitly" {
    try std.testing.expectError(
        error.InvalidAllocationFraction,
        humusAllocationFraction(0.2, 5_000_000.0),
    );
}
