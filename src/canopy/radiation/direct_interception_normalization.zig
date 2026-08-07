const std = @import("std");

pub const Components = struct {
    leaf_shortwave: []f64,
    stalk_shortwave: []f64,
    standing_dead_shortwave: []f64,
    leaf_par: []f64,
    stalk_par: []f64,
    standing_dead_par: []f64,
};

pub const SpeciesAccumulators = struct {
    absorbed: Components,
    backscattered: Components,
    forward_scattered: Components,
};

pub const State = struct {
    direct_transmittance_above: f64,
    current_layer_interception_fraction: f64,
};

/// `hour1.f` lines 1446--1474. Normalizes direct interception only when the
/// preceding plus current interception exceeds one, preserving every
/// species multiplication in source order.
pub fn apply(
    preceding_interception_fraction: f64,
    negligible_interception_fraction: f64,
    state: *State,
    accumulators: SpeciesAccumulators,
) !?f64 {
    const species_count = try validate(
        preceding_interception_fraction,
        negligible_interception_fraction,
        state.*,
        accumulators,
    );
    if (preceding_interception_fraction +
        state.current_layer_interception_fraction <= 1.0)
        return null;

    const normalization = if (state.current_layer_interception_fraction >
        negligible_interception_fraction)
        (1.0 - preceding_interception_fraction) /
            ((1.0 - preceding_interception_fraction) -
                (1.0 - preceding_interception_fraction -
                    state.current_layer_interception_fraction))
    else
        0.0;
    if (!std.math.isFinite(normalization) or normalization < 0)
        return error.InvalidDirectInterceptionNormalization;
    state.direct_transmittance_above =
        state.direct_transmittance_above * normalization;
    state.current_layer_interception_fraction =
        state.current_layer_interception_fraction * normalization;
    for (0..species_count) |species| {
        scaleAt(accumulators.absorbed, species, normalization);
        scaleAt(accumulators.backscattered, species, normalization);
        scaleAt(accumulators.forward_scattered, species, normalization);
    }
    return normalization;
}

fn scaleAt(components: Components, species: usize, scale: f64) void {
    components.leaf_shortwave[species] =
        components.leaf_shortwave[species] * scale;
    components.stalk_shortwave[species] =
        components.stalk_shortwave[species] * scale;
    components.standing_dead_shortwave[species] =
        components.standing_dead_shortwave[species] * scale;
    components.leaf_par[species] = components.leaf_par[species] * scale;
    components.stalk_par[species] = components.stalk_par[species] * scale;
    components.standing_dead_par[species] =
        components.standing_dead_par[species] * scale;
}

fn validate(
    preceding: f64,
    negligible: f64,
    state: State,
    accumulators: SpeciesAccumulators,
) !usize {
    inline for (.{
        preceding,
        negligible,
        state.direct_transmittance_above,
        state.current_layer_interception_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDirectInterceptionNormalizationInput;
    const species_count = accumulators.absorbed.leaf_shortwave.len;
    if (species_count == 0)
        return error.ZeroDirectInterceptionSpeciesExtent;
    inline for (@typeInfo(SpeciesAccumulators).@"struct".fields) |group_field| {
        const components = @field(accumulators, group_field.name);
        inline for (@typeInfo(Components).@"struct".fields) |field| {
            const values = @field(components, field.name);
            if (values.len != species_count)
                return error.DirectInterceptionNormalizationDimensionMismatch;
            for (values) |value| if (!std.math.isFinite(value))
                return error.InvalidDirectInterceptionNormalizationInput;
        }
    }
    return species_count;
}

fn partitioned(values: []f64, species_count: usize) Components {
    return .{
        .leaf_shortwave = values[0 * species_count .. 1 * species_count],
        .stalk_shortwave = values[1 * species_count .. 2 * species_count],
        .standing_dead_shortwave = values[2 * species_count .. 3 * species_count],
        .leaf_par = values[3 * species_count .. 4 * species_count],
        .stalk_par = values[4 * species_count .. 5 * species_count],
        .standing_dead_par = values[5 * species_count .. 6 * species_count],
    };
}

test "excess direct interception scales state and all species fields" {
    var absorbed = [_]f64{ 2, 4 } ** 6;
    var back = [_]f64{ 6, 8 } ** 6;
    var forward = [_]f64{ 10, 12 } ** 6;
    var state: State = .{
        .direct_transmittance_above = 0.9,
        .current_layer_interception_fraction = 0.4,
    };
    const normalization = (try apply(0.8, 1.0e-12, &state, .{
        .absorbed = partitioned(&absorbed, 2),
        .backscattered = partitioned(&back, 2),
        .forward_scattered = partitioned(&forward, 2),
    })).?;
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), normalization, 1e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.45),
        state.direct_transmittance_above,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.2),
        state.current_layer_interception_fraction,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 1), absorbed[0], 3e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), absorbed[1], 3e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), back[0], 3e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 4), back[1], 3e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 5), forward[0], 3e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6), forward[1], 3e-15);
}

test "interception not exceeding one leaves all state unchanged" {
    var values = [_]f64{42} ** 6;
    var state: State = .{
        .direct_transmittance_above = 0.9,
        .current_layer_interception_fraction = 0.2,
    };
    const result = try apply(0.7, 0, &state, .{
        .absorbed = partitioned(&values, 1),
        .backscattered = partitioned(&values, 1),
        .forward_scattered = partitioned(&values, 1),
    });
    try std.testing.expect(result == null);
    try std.testing.expectEqual(@as(f64, 0.9), state.direct_transmittance_above);
    try std.testing.expectEqual(@as(f64, 42), values[0]);
}
