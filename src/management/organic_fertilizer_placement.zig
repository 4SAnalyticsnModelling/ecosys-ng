const std = @import("std");

pub const Destination = union(enum) {
    surface_litter,
    soil_layer: usize,
};

/// `hour1.f` lines 623--635. The source gate considers organic carbon only;
/// positive N or P without residue/manure C does not admit an application.
pub fn determine(
    plant_residue_carbon_g_c_per_m2: f64,
    manure_carbon_g_c_per_m2: f64,
    application_depth_m: f64,
    soil_layer_bottom_depth_m: []const f64,
) !?Destination {
    inline for (.{
        plant_residue_carbon_g_c_per_m2,
        manure_carbon_g_c_per_m2,
        application_depth_m,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteOrganicFertilizerPlacementInput;
        if (value < 0) return error.InvalidOrganicFertilizerPlacementInput;
    }
    var previous_bottom_m: f64 = 0;
    for (soil_layer_bottom_depth_m) |bottom_m| {
        if (!std.math.isFinite(bottom_m) or bottom_m <= previous_bottom_m)
            return error.InvalidOrganicFertilizerLayerBoundary;
        previous_bottom_m = bottom_m;
    }

    if (plant_residue_carbon_g_c_per_m2 +
        manure_carbon_g_c_per_m2 <= 0)
        return null;
    if (application_depth_m <= 0) return .surface_litter;
    for (soil_layer_bottom_depth_m, 0..) |bottom_m, layer|
        if (bottom_m >= application_depth_m)
            return .{ .soil_layer = layer };
    return error.OrganicFertilizerApplicationBelowSoilProfile;
}

test "organic placement preserves carbon-only gate" {
    try std.testing.expect((try determine(0, 0, 0, &.{0.1})) == null);
    try std.testing.expect(
        (try determine(1, 0, 0, &.{0.1})).? == .surface_litter,
    );
}

test "organic placement selects first inclusive layer bottom" {
    const first = (try determine(0, 1, 0.1, &.{ 0.1, 0.3 })).?;
    try std.testing.expectEqual(@as(usize, 0), first.soil_layer);
    const second = (try determine(1, 0, 0.1001, &.{ 0.1, 0.3 })).?;
    try std.testing.expectEqual(@as(usize, 1), second.soil_layer);
}

test "invalid topology and below-profile depth fail explicitly" {
    try std.testing.expectError(
        error.InvalidOrganicFertilizerLayerBoundary,
        determine(1, 0, 0.2, &.{ 0.1, 0.1 }),
    );
    try std.testing.expectError(
        error.OrganicFertilizerApplicationBelowSoilProfile,
        determine(1, 0, 1, &.{ 0.1, 0.3 }),
    );
}
