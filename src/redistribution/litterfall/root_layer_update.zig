const std = @import("std");

pub const Dimensions = struct {
    soil_layer_count: usize,
    root_litter_class_count: usize,
    biochemical_fraction_count: usize,

    pub fn elementCount(self: Dimensions) !usize {
        if (self.soil_layer_count == 0 or self.root_litter_class_count == 0 or
            self.biochemical_fraction_count == 0)
            return error.InvalidRootLitterfallDimensions;
        return std.math.mul(
            usize,
            try std.math.mul(usize, self.soil_layer_count, self.root_litter_class_count),
            self.biochemical_fraction_count,
        ) catch return error.RootLitterfallDimensionOverflow;
    }
};

pub const OrganicMatterState = struct {
    carbon_g: []f64,
    nitrogen_g: []f64,
    phosphorus_g: []f64,
};

pub const RootLitterfall = struct {
    carbon_g: []const f64,
    nitrogen_g: []const f64,
    phosphorus_g: []const f64,
};

fn index(dimensions: Dimensions, layer: usize, class: usize, fraction: usize) usize {
    return (layer * dimensions.root_litter_class_count + class) *
        dimensions.biochemical_fraction_count + fraction;
}

/// Runtime-dimension translation of redist.f lines 6057--6069.
///
/// Traversal remains layer-major, then root-litter class K, then biochemical
/// fraction M. The legacy run uses two K classes and five M fractions.
pub fn apply(
    dimensions: Dimensions,
    state: OrganicMatterState,
    litterfall: RootLitterfall,
) !void {
    const count = try dimensions.elementCount();
    inline for (.{
        state.carbon_g.len,
        state.nitrogen_g.len,
        state.phosphorus_g.len,
        litterfall.carbon_g.len,
        litterfall.nitrogen_g.len,
        litterfall.phosphorus_g.len,
    }) |length| if (length != count) return error.RootLitterfallDimensionMismatch;

    for (0..dimensions.soil_layer_count) |layer| {
        for (0..dimensions.root_litter_class_count) |class| {
            for (0..dimensions.biochemical_fraction_count) |fraction| {
                const i = index(dimensions, layer, class, fraction);
                inline for (.{
                    state.carbon_g[i],      state.nitrogen_g[i],      state.phosphorus_g[i],
                    litterfall.carbon_g[i], litterfall.nitrogen_g[i], litterfall.phosphorus_g[i],
                }) |value| if (!std.math.isFinite(value)) return error.InvalidRootLitterfallValue;
                state.carbon_g[i] = state.carbon_g[i] + litterfall.carbon_g[i];
                state.nitrogen_g[i] = state.nitrogen_g[i] + litterfall.nitrogen_g[i];
                state.phosphorus_g[i] = state.phosphorus_g[i] + litterfall.phosphorus_g[i];
                inline for (.{ state.carbon_g[i], state.nitrogen_g[i], state.phosphorus_g[i] }) |value|
                    if (!std.math.isFinite(value)) return error.NonFiniteRootLitterfallPool;
            }
        }
    }
}

test "REDIST root litterfall preserves legacy layer K M traversal dimensions" {
    const dimensions = Dimensions{ .soil_layer_count = 2, .root_litter_class_count = 2, .biochemical_fraction_count = 5 };
    var carbon = [_]f64{0} ** 20;
    var nitrogen = [_]f64{0} ** 20;
    var phosphorus = [_]f64{0} ** 20;
    var carbon_flux = [_]f64{0} ** 20;
    var nitrogen_flux = [_]f64{0} ** 20;
    var phosphorus_flux = [_]f64{0} ** 20;
    carbon_flux[index(dimensions, 0, 0, 0)] = 1;
    nitrogen_flux[index(dimensions, 1, 1, 4)] = 2;
    phosphorus_flux[index(dimensions, 1, 0, 3)] = 3;
    try apply(dimensions, .{ .carbon_g = &carbon, .nitrogen_g = &nitrogen, .phosphorus_g = &phosphorus }, .{ .carbon_g = &carbon_flux, .nitrogen_g = &nitrogen_flux, .phosphorus_g = &phosphorus_flux });
    try std.testing.expectEqual(@as(f64, 1), carbon[0]);
    try std.testing.expectEqual(@as(f64, 2), nitrogen[19]);
    try std.testing.expectEqual(@as(f64, 3), phosphorus[13]);
}

test "REDIST root litterfall supports runtime class and fraction dimensions" {
    const dimensions = Dimensions{ .soil_layer_count = 3, .root_litter_class_count = 4, .biochemical_fraction_count = 2 };
    var carbon = [_]f64{1} ** 24;
    var nitrogen = [_]f64{2} ** 24;
    var phosphorus = [_]f64{3} ** 24;
    const carbon_flux = [_]f64{0.5} ** 24;
    const nitrogen_flux = [_]f64{1} ** 24;
    const phosphorus_flux = [_]f64{-1} ** 24;
    try apply(dimensions, .{ .carbon_g = &carbon, .nitrogen_g = &nitrogen, .phosphorus_g = &phosphorus }, .{ .carbon_g = &carbon_flux, .nitrogen_g = &nitrogen_flux, .phosphorus_g = &phosphorus_flux });
    try std.testing.expectEqual(@as(f64, 1.5), carbon[23]);
    try std.testing.expectEqual(@as(f64, 3), nitrogen[23]);
    try std.testing.expectEqual(@as(f64, 2), phosphorus[23]);
}

test "REDIST root litterfall rejects dimensions invalid values and overflow" {
    try std.testing.expectError(error.InvalidRootLitterfallDimensions, (Dimensions{ .soil_layer_count = 0, .root_litter_class_count = 2, .biochemical_fraction_count = 5 }).elementCount());
    const dimensions = Dimensions{ .soil_layer_count = 1, .root_litter_class_count = 1, .biochemical_fraction_count = 1 };
    var carbon = [_]f64{std.math.floatMax(f64)};
    var nitrogen = [_]f64{0};
    var phosphorus = [_]f64{0};
    const carbon_flux = [_]f64{std.math.floatMax(f64)};
    const zero = [_]f64{0};
    try std.testing.expectError(error.NonFiniteRootLitterfallPool, apply(dimensions, .{ .carbon_g = &carbon, .nitrogen_g = &nitrogen, .phosphorus_g = &phosphorus }, .{ .carbon_g = &carbon_flux, .nitrogen_g = &zero, .phosphorus_g = &zero }));
}
