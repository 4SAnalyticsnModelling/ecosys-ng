const std = @import("std");

pub const Inputs = struct {
    /// Runtime replacement for HOUR1's source default 0.55E+06 g Mg^-1.
    surface_litter_organic_carbon_g_per_megagram: f64,
};

/// HOUR1 line 2336. Publishes the required runtime surface-litter organic
/// carbon concentration without embedding a fixed parameter table.
pub fn compute(inputs: Inputs) !f64 {
    const concentration =
        inputs.surface_litter_organic_carbon_g_per_megagram;
    if (!std.math.isFinite(concentration))
        return error.NonFiniteSurfaceLitterOrganicCarbonConcentration;
    if (concentration < 0 or concentration > 1.0e6)
        return error.InvalidSurfaceLitterOrganicCarbonConcentration;
    return concentration;
}

test "runtime input can reproduce exact HOUR1 surface litter concentration" {
    try std.testing.expectEqual(
        @as(f64, 0.55e6),
        try compute(.{
            .surface_litter_organic_carbon_g_per_megagram = 0.55e6,
        }),
    );
}

test "concentration above one megagram per megagram fails" {
    try std.testing.expectError(
        error.InvalidSurfaceLitterOrganicCarbonConcentration,
        compute(.{
            .surface_litter_organic_carbon_g_per_megagram = 1.0e6 + 1,
        }),
    );
}
