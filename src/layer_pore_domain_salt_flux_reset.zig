const std = @import("std");

pub const SaltMode = enum { static_equilibrium, dynamic_transport };

/// Source-ordered micropore species from HOUR1 lines 2795-2844.
pub const MicroporeSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_3,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
    band_hydrogen_phosphate_0,
    band_hydrogen_phosphate_3,
    band_iron_phosphate_1,
    band_iron_phosphate_2,
    band_calcium_phosphate_0,
    band_calcium_phosphate_1,
    band_calcium_phosphate_2,
    band_magnesium_phosphate_1,
};

/// Source-ordered macropore species from HOUR1 lines 2845-2893.
pub const MacroporeSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_3,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
    band_hydrogen_phosphate_0,
    band_hydrogen_phosphate_3,
    band_iron_phosphate_1,
    band_iron_phosphate_2,
    band_calcium_phosphate_0,
    band_calcium_phosphate_1,
    band_calcium_phosphate_2,
    band_magnesium_phosphate_1,
};

/// Salt fluxes for one runtime-selected layer and lateral direction.
/// Values retain legacy species mass-per-timestep units.
pub const LayerFluxes = struct {
    micropore: []f64,
    macropore: []f64,
};

pub const ResetError = error{SpeciesCountMismatch};

/// Translates HOUR1 lines 2795-2893. Static-equilibrium mode is a no-op.
pub fn reset(mode: SaltMode, layers: []LayerFluxes) ResetError!void {
    if (mode == .static_equilibrium) return;
    const micropore_count = std.meta.fields(MicroporeSpecies).len;
    const macropore_count = std.meta.fields(MacroporeSpecies).len;
    for (layers) |layer| {
        if (layer.micropore.len != micropore_count or layer.macropore.len != macropore_count) {
            return error.SpeciesCountMismatch;
        }
    }
    for (layers) |layer| {
        for (layer.micropore) |*flux| flux.* = 0.0;
        for (layer.macropore) |*flux| flux.* = 0.0;
    }
}

test "dynamic mode clears all runtime layers and pore-domain species" {
    const allocator = std.testing.allocator;
    const layers = try allocator.alloc(LayerFluxes, 7);
    defer allocator.free(layers);
    for (layers) |*layer| {
        layer.micropore = try allocator.alloc(f64, std.meta.fields(MicroporeSpecies).len);
        layer.macropore = try allocator.alloc(f64, std.meta.fields(MacroporeSpecies).len);
        @memset(layer.micropore, 2.0);
        @memset(layer.macropore, 3.0);
    }
    defer for (layers) |layer| {
        allocator.free(layer.micropore);
        allocator.free(layer.macropore);
    };

    try reset(.dynamic_transport, layers);

    for (layers) |layer| {
        for (layer.micropore) |flux| try std.testing.expectEqual(@as(f64, 0.0), flux);
        for (layer.macropore) |flux| try std.testing.expectEqual(@as(f64, 0.0), flux);
    }
}

test "dimension mismatch is detected before any layer mutation" {
    var valid_micro = [_]f64{4.0} ** std.meta.fields(MicroporeSpecies).len;
    var valid_macro = [_]f64{5.0} ** std.meta.fields(MacroporeSpecies).len;
    var short_micro = [_]f64{6.0};
    var layers = [_]LayerFluxes{
        .{ .micropore = &valid_micro, .macropore = &valid_macro },
        .{ .micropore = &short_micro, .macropore = &valid_macro },
    };
    try std.testing.expectError(error.SpeciesCountMismatch, reset(.dynamic_transport, &layers));
    try std.testing.expectEqual(@as(f64, 4.0), valid_micro[0]);
    try std.testing.expectEqual(@as(f64, 5.0), valid_macro[0]);
}
