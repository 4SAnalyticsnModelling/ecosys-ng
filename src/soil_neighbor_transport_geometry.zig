const std = @import("std");

pub const Inputs = struct {
    source_path_length_m: []const f64,
    destination_path_length_m: []const f64,
    face_area_m2: []const f64,
    boundary_kind: []const BoundaryKind,
    timestep_h: f64,
};

pub const BoundaryKind = enum {
    internal,
    top_surface,
};

pub const Parameters = struct {
    dispersivity_coefficient: f64,
    distance_exponent: f64,
};

pub const State = struct {
    mean_transport_distance_m: []f64,
    face_area_per_distance_m: []f64,
    dispersivity_m2_per_step: []f64,
};

/// Exact face-equation translation of legacy `STARTS` lines 741--748.
///
/// Internal faces use half the sum of source and destination path lengths.
/// The top boundary uses half the source layer thickness, matching lines
/// 745--748 rather than inventing a ghost-cell distance.
pub fn initialize(
    state: State,
    inputs: Inputs,
    parameters: Parameters,
) !void {
    const face_count = inputs.source_path_length_m.len;
    if (face_count == 0 or
        inputs.destination_path_length_m.len != face_count or
        inputs.face_area_m2.len != face_count or
        inputs.boundary_kind.len != face_count or
        state.mean_transport_distance_m.len != face_count or
        state.face_area_per_distance_m.len != face_count or
        state.dispersivity_m2_per_step.len != face_count)
        return error.SoilNeighborTransportDimensionMismatch;
    inline for (.{
        inputs.timestep_h,
        parameters.dispersivity_coefficient,
        parameters.distance_exponent,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSoilNeighborTransportParameter;
    }
    if (inputs.timestep_h <= 0 or
        parameters.dispersivity_coefficient < 0 or
        parameters.distance_exponent <= 0)
        return error.InvalidSoilNeighborTransportParameter;

    for (0..face_count) |face| {
        inline for (.{
            inputs.source_path_length_m[face],
            inputs.destination_path_length_m[face],
            inputs.face_area_m2[face],
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteSoilNeighborTransportInput;
        if (inputs.source_path_length_m[face] <= 0 or
            inputs.face_area_m2[face] <= 0 or
            (inputs.boundary_kind[face] == .internal and
                inputs.destination_path_length_m[face] <= 0))
            return error.InvalidSoilNeighborTransportInput;
        const distance_m =
            if (inputs.boundary_kind[face] == .top_surface)
                0.5 * inputs.source_path_length_m[face]
            else
                0.5 * (inputs.source_path_length_m[face] +
                    inputs.destination_path_length_m[face]);
        const area_per_distance_m =
            inputs.face_area_m2[face] / distance_m;
        const dispersivity_m2_per_step =
            parameters.dispersivity_coefficient *
            std.math.pow(f64, distance_m, parameters.distance_exponent) *
            inputs.timestep_h;
        inline for (.{
            distance_m,
            area_per_distance_m,
            dispersivity_m2_per_step,
        }) |candidate| if (!std.math.isFinite(candidate))
            return error.SoilNeighborTransportOverflow;
    }

    for (0..face_count) |face| {
        state.mean_transport_distance_m[face] =
            if (inputs.boundary_kind[face] == .top_surface)
                0.5 * inputs.source_path_length_m[face]
            else
                0.5 * (inputs.source_path_length_m[face] +
                    inputs.destination_path_length_m[face]);
        state.face_area_per_distance_m[face] =
            inputs.face_area_m2[face] /
            state.mean_transport_distance_m[face];
        state.dispersivity_m2_per_step[face] =
            parameters.dispersivity_coefficient *
            std.math.pow(
                f64,
                state.mean_transport_distance_m[face],
                parameters.distance_exponent,
            ) * inputs.timestep_h;
    }
}

test "STARTS internal and top faces preserve distinct distance equations" {
    var distance = [_]f64{0.0} ** 2;
    var area_per_distance = [_]f64{0.0} ** 2;
    var dispersivity = [_]f64{0.0} ** 2;
    try initialize(.{
        .mean_transport_distance_m = &distance,
        .face_area_per_distance_m = &area_per_distance,
        .dispersivity_m2_per_step = &dispersivity,
    }, .{
        .source_path_length_m = &.{ 0.2, 0.2 },
        .destination_path_length_m = &.{ 0.4, 0.0 },
        .face_area_m2 = &.{ 10, 10 },
        .boundary_kind = &.{ .internal, .top_surface },
        .timestep_h = 0.25,
    }, .{
        .dispersivity_coefficient = 0.20,
        .distance_exponent = 1.07,
    });

    try std.testing.expectApproxEqAbs(@as(f64, 0.3), distance[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), distance[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0 / 0.3), area_per_distance[0], 1e-14);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.20 * std.math.pow(f64, 0.3, 1.07) * 0.25),
        dispersivity[0],
        1e-15,
    );
}

test "late invalid face preserves all transport geometry" {
    var distance = [_]f64{ 7, 8 };
    var area_per_distance = [_]f64{ 9, 10 };
    var dispersivity = [_]f64{ 11, 12 };
    const before = distance;
    try std.testing.expectError(
        error.InvalidSoilNeighborTransportInput,
        initialize(.{
            .mean_transport_distance_m = &distance,
            .face_area_per_distance_m = &area_per_distance,
            .dispersivity_m2_per_step = &dispersivity,
        }, .{
            .source_path_length_m = &.{ 0.2, 0.2 },
            .destination_path_length_m = &.{ 0.4, -1 },
            .face_area_m2 = &.{ 10, 10 },
            .boundary_kind = &.{ .internal, .internal },
            .timestep_h = 1,
        }, .{
            .dispersivity_coefficient = 0.2,
            .distance_exponent = 1.07,
        }),
    );
    try std.testing.expectEqualSlices(f64, &before, &distance);
}
