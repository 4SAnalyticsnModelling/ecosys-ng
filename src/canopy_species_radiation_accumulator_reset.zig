const std = @import("std");

pub const RadiationComponents = struct {
    leaf_shortwave: []f64,
    stalk_shortwave: []f64,
    standing_dead_shortwave: []f64,
    leaf_par: []f64,
    stalk_par: []f64,
    standing_dead_par: []f64,
};

pub const Accumulators = struct {
    direct_absorbed: RadiationComponents,
    diffuse_absorbed: RadiationComponents,
    direct_backscattered: RadiationComponents,
    diffuse_backscattered: RadiationComponents,
    direct_forward_scattered: RadiationComponents,
    diffuse_forward_scattered: RadiationComponents,
};

/// HOUR1 lines 1235--1271. Resets all 36 radiation accumulators for each
/// runtime species in the exact source group and component order.
pub fn apply(species_count: usize, accumulators: Accumulators) !void {
    try validate(species_count, accumulators);
    for (0..species_count) |species| {
        resetComponentsAt(accumulators.direct_absorbed, species);
        resetComponentsAt(accumulators.diffuse_absorbed, species);
        resetComponentsAt(accumulators.direct_backscattered, species);
        resetComponentsAt(accumulators.diffuse_backscattered, species);
        resetComponentsAt(accumulators.direct_forward_scattered, species);
        resetComponentsAt(accumulators.diffuse_forward_scattered, species);
    }
}

fn resetComponentsAt(components: RadiationComponents, species: usize) void {
    components.leaf_shortwave[species] = 0.0;
    components.stalk_shortwave[species] = 0.0;
    components.standing_dead_shortwave[species] = 0.0;
    components.leaf_par[species] = 0.0;
    components.stalk_par[species] = 0.0;
    components.standing_dead_par[species] = 0.0;
}

fn validate(species_count: usize, accumulators: Accumulators) !void {
    if (species_count == 0)
        return error.ZeroCanopyRadiationSpeciesExtent;
    inline for (@typeInfo(Accumulators).@"struct".fields) |group_field| {
        const components = @field(accumulators, group_field.name);
        inline for (@typeInfo(RadiationComponents).@"struct".fields) |field|
            if (@field(components, field.name).len != species_count)
                return error.CanopyRadiationAccumulatorDimensionMismatch;
    }
}

fn repeatedComponents(values: []f64) RadiationComponents {
    return .{
        .leaf_shortwave = values,
        .stalk_shortwave = values,
        .standing_dead_shortwave = values,
        .leaf_par = values,
        .stalk_par = values,
        .standing_dead_par = values,
    };
}

test "all thirty six accumulators reset for runtime species" {
    var direct_absorbed = [_]f64{1} ** 7;
    var diffuse_absorbed = [_]f64{2} ** 7;
    var direct_back = [_]f64{3} ** 7;
    var diffuse_back = [_]f64{4} ** 7;
    var direct_forward = [_]f64{5} ** 7;
    var diffuse_forward = [_]f64{6} ** 7;
    const accumulators: Accumulators = .{
        .direct_absorbed = repeatedComponents(&direct_absorbed),
        .diffuse_absorbed = repeatedComponents(&diffuse_absorbed),
        .direct_backscattered = repeatedComponents(&direct_back),
        .diffuse_backscattered = repeatedComponents(&diffuse_back),
        .direct_forward_scattered = repeatedComponents(&direct_forward),
        .diffuse_forward_scattered = repeatedComponents(&diffuse_forward),
    };
    try apply(7, accumulators);
    inline for (@typeInfo(Accumulators).@"struct".fields) |group_field| {
        const components = @field(accumulators, group_field.name);
        inline for (@typeInfo(RadiationComponents).@"struct".fields) |field|
            for (@field(components, field.name)) |value|
                try std.testing.expectEqual(@as(f64, 0), value);
    }
}

test "dimension mismatch leaves every accumulator unchanged" {
    var valid = [_]f64{ 41, 41 };
    var invalid = [_]f64{42};
    var accumulators: Accumulators = .{
        .direct_absorbed = repeatedComponents(&valid),
        .diffuse_absorbed = repeatedComponents(&valid),
        .direct_backscattered = repeatedComponents(&valid),
        .diffuse_backscattered = repeatedComponents(&valid),
        .direct_forward_scattered = repeatedComponents(&valid),
        .diffuse_forward_scattered = repeatedComponents(&valid),
    };
    accumulators.diffuse_forward_scattered.standing_dead_par = &invalid;
    try std.testing.expectError(
        error.CanopyRadiationAccumulatorDimensionMismatch,
        apply(2, accumulators),
    );
    try std.testing.expectEqualSlices(f64, &.{ 41, 41 }, &valid);
    try std.testing.expectEqual(@as(f64, 42), invalid[0]);
}
