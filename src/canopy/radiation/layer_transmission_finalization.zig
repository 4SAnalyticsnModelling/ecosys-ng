const std = @import("std");

pub const LayerExposure = enum {
    exposed,
    submerged,
};

pub const AboveBoundary = struct {
    forward_scattered_shortwave_megajoules_per_m2_h: f64,
    forward_scattered_par_umol_per_m2_s: f64,
    direct_transmittance: f64,
    diffuse_transmittance: f64,
};

pub const LayerState = struct {
    accumulated_direct_interception_fraction: f64,
    accumulated_diffuse_interception_fraction: f64,
    forward_scattered_shortwave_megajoules_per_m2_h: f64,
    forward_scattered_par_umol_per_m2_s: f64,
    direct_transmittance: f64,
    direct_interception_fraction: f64,
    diffuse_transmittance: f64,
};

/// `hour1.f` lines 1579--1591. Finalizes an exposed layer or propagates the
/// boundary above through a submerged layer in exact source assignment order.
pub fn apply(
    exposure: LayerExposure,
    current_direct_interception_fraction: f64,
    current_diffuse_interception_fraction: f64,
    above: AboveBoundary,
    state: *LayerState,
) !void {
    try validate(
        exposure,
        current_direct_interception_fraction,
        current_diffuse_interception_fraction,
        above,
        state.*,
    );
    switch (exposure) {
        .exposed => {
            state.accumulated_direct_interception_fraction +=
                current_direct_interception_fraction;
            state.accumulated_diffuse_interception_fraction +=
                current_diffuse_interception_fraction;
            state.direct_transmittance =
                1.0 - state.accumulated_direct_interception_fraction;
            state.direct_interception_fraction =
                1.0 - state.direct_transmittance;
            state.diffuse_transmittance =
                1.0 - state.accumulated_diffuse_interception_fraction;
        },
        .submerged => {
            state.forward_scattered_shortwave_megajoules_per_m2_h =
                above.forward_scattered_shortwave_megajoules_per_m2_h;
            state.forward_scattered_par_umol_per_m2_s =
                above.forward_scattered_par_umol_per_m2_s;
            state.direct_transmittance = above.direct_transmittance;
            state.direct_interception_fraction =
                1.0 - state.direct_transmittance;
            state.diffuse_transmittance = above.diffuse_transmittance;
        },
    }
}

fn validate(
    exposure: LayerExposure,
    direct_delta: f64,
    diffuse_delta: f64,
    above: AboveBoundary,
    state: LayerState,
) !void {
    inline for (.{
        direct_delta,
        diffuse_delta,
        above.forward_scattered_shortwave_megajoules_per_m2_h,
        above.forward_scattered_par_umol_per_m2_s,
        above.direct_transmittance,
        above.diffuse_transmittance,
        state.accumulated_direct_interception_fraction,
        state.accumulated_diffuse_interception_fraction,
        state.forward_scattered_shortwave_megajoules_per_m2_h,
        state.forward_scattered_par_umol_per_m2_s,
        state.direct_transmittance,
        state.direct_interception_fraction,
        state.diffuse_transmittance,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyLayerTransmissionInput;
    if (above.direct_transmittance > 1 or above.diffuse_transmittance > 1)
        return error.InvalidCanopyLayerTransmissionInput;
    if (exposure == .exposed and
        (state.accumulated_direct_interception_fraction + direct_delta > 1 or
            state.accumulated_diffuse_interception_fraction +
                diffuse_delta > 1))
        return error.CanopyLayerInterceptionExceedsOne;
}

test "exposed layer finalizes accumulated interception in source order" {
    var state: LayerState = .{
        .accumulated_direct_interception_fraction = 0.2,
        .accumulated_diffuse_interception_fraction = 0.1,
        .forward_scattered_shortwave_megajoules_per_m2_h = 7,
        .forward_scattered_par_umol_per_m2_s = 8,
        .direct_transmittance = 9,
        .direct_interception_fraction = 10,
        .diffuse_transmittance = 11,
    };
    try apply(.exposed, 0.3, 0.4, .{
        .forward_scattered_shortwave_megajoules_per_m2_h = 1,
        .forward_scattered_par_umol_per_m2_s = 2,
        .direct_transmittance = 0.8,
        .diffuse_transmittance = 0.7,
    }, &state);
    try std.testing.expectEqual(@as(f64, 0.5), state.accumulated_direct_interception_fraction);
    try std.testing.expectEqual(@as(f64, 0.5), state.accumulated_diffuse_interception_fraction);
    try std.testing.expectEqual(@as(f64, 0.5), state.direct_transmittance);
    try std.testing.expectEqual(@as(f64, 0.5), state.direct_interception_fraction);
    try std.testing.expectEqual(@as(f64, 0.5), state.diffuse_transmittance);
    try std.testing.expectEqual(@as(f64, 7), state.forward_scattered_shortwave_megajoules_per_m2_h);
}

test "submerged layer inherits the boundary above" {
    var state: LayerState = .{
        .accumulated_direct_interception_fraction = 0.2,
        .accumulated_diffuse_interception_fraction = 0.3,
        .forward_scattered_shortwave_megajoules_per_m2_h = 7,
        .forward_scattered_par_umol_per_m2_s = 8,
        .direct_transmittance = 0.1,
        .direct_interception_fraction = 0.9,
        .diffuse_transmittance = 0.2,
    };
    try apply(.submerged, 0, 0, .{
        .forward_scattered_shortwave_megajoules_per_m2_h = 1,
        .forward_scattered_par_umol_per_m2_s = 2,
        .direct_transmittance = 0.8,
        .diffuse_transmittance = 0.7,
    }, &state);
    try std.testing.expectEqual(@as(f64, 1), state.forward_scattered_shortwave_megajoules_per_m2_h);
    try std.testing.expectEqual(@as(f64, 2), state.forward_scattered_par_umol_per_m2_s);
    try std.testing.expectEqual(@as(f64, 0.8), state.direct_transmittance);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.direct_interception_fraction, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.7), state.diffuse_transmittance);
    try std.testing.expectEqual(@as(f64, 0.2), state.accumulated_direct_interception_fraction);
}
