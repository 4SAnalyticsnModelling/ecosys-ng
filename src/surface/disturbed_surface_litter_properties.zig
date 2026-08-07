const std = @import("std");

pub const DisturbanceStatus = enum {
    unchanged,
    disturbed,
};

pub const Inputs = struct {
    disturbance_status: DisturbanceStatus,
    surface_litter_volume_m3: f64,
    surface_litter_mass_megagrams: f64,
    negligible_volume_m3: f64,
    reference_fine_litter_bulk_density_megagrams_m3: f64,
};

/// `hour1.f` lines 1900--1905. Null means the source disturbance gate did not
/// execute and the caller-owned density must remain unchanged.
pub fn computeBulkDensity(inputs: Inputs) !?f64 {
    try validate(inputs);
    if (inputs.disturbance_status == .unchanged) return null;
    if (inputs.surface_litter_volume_m3 > inputs.negligible_volume_m3)
        return inputs.surface_litter_mass_megagrams /
            inputs.surface_litter_volume_m3;
    return inputs.reference_fine_litter_bulk_density_megagrams_m3;
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.surface_litter_volume_m3,
        inputs.surface_litter_mass_megagrams,
        inputs.negligible_volume_m3,
        inputs.reference_fine_litter_bulk_density_megagrams_m3,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDisturbedSurfaceLitterInput;
    if (inputs.reference_fine_litter_bulk_density_megagrams_m3 == 0)
        return error.InvalidDisturbedSurfaceLitterInput;
}

test "disturbed litter density uses mass divided by positive volume" {
    const density = (try computeBulkDensity(.{
        .disturbance_status = .disturbed,
        .surface_litter_volume_m3 = 2,
        .surface_litter_mass_megagrams = 6,
        .negligible_volume_m3 = 1e-12,
        .reference_fine_litter_bulk_density_megagrams_m3 = 0.1,
    })).?;
    try std.testing.expectEqual(@as(f64, 3), density);
}

test "negligible disturbed volume uses source reference density" {
    const density = (try computeBulkDensity(.{
        .disturbance_status = .disturbed,
        .surface_litter_volume_m3 = 0,
        .surface_litter_mass_megagrams = 0,
        .negligible_volume_m3 = 1e-12,
        .reference_fine_litter_bulk_density_megagrams_m3 = 0.2,
    })).?;
    try std.testing.expectEqual(@as(f64, 0.2), density);
}

test "undisturbed litter leaves caller density unchanged" {
    try std.testing.expect((try computeBulkDensity(.{
        .disturbance_status = .unchanged,
        .surface_litter_volume_m3 = 1,
        .surface_litter_mass_megagrams = 2,
        .negligible_volume_m3 = 1e-12,
        .reference_fine_litter_bulk_density_megagrams_m3 = 0.2,
    })) == null);
}
