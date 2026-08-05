const std = @import("std");

pub const Inputs = struct {
    cell_area_m2: f64,
    diffuse_shortwave_megajoules_m2_h: f64,
    diffuse_par_umol_m2_s: f64,
};

pub const InitialState = struct {
    inverse_cell_area_per_m2: f64,
    inverse_azimuth_area_per_m2: f64,
    working_diffuse_shortwave_megajoules_m2_h: f64,
    working_diffuse_par_umol_m2_s: f64,
    top_direct_transmittance: f64,
    top_diffuse_transmittance: f64,
    top_scattered_shortwave_megajoules_m2_h: f64,
    top_scattered_par_umol_m2_s: f64,
    interception_stop_fraction: f64,
};

/// HOUR1 lines 1128--1136. Initializes the canopy-top boundary and
/// interception accumulators in exact source assignment order.
pub fn initialize(inputs: Inputs) !InitialState {
    if (!std.math.isFinite(inputs.cell_area_m2) or inputs.cell_area_m2 <= 0)
        return error.InvalidCanopyInterceptionCellArea;
    if (!std.math.isFinite(inputs.diffuse_shortwave_megajoules_m2_h) or
        !std.math.isFinite(inputs.diffuse_par_umol_m2_s))
        return error.NonFiniteCanopyInterceptionInput;
    if (inputs.diffuse_shortwave_megajoules_m2_h < 0 or
        inputs.diffuse_par_umol_m2_s < 0)
        return error.InvalidCanopyInterceptionInput;

    const inverse_cell_area_per_m2 = 1.00 / inputs.cell_area_m2;
    const inverse_azimuth_area_per_m2 = 0.25 / inputs.cell_area_m2;
    const working_diffuse_shortwave_megajoules_m2_h =
        inputs.diffuse_shortwave_megajoules_m2_h;
    const working_diffuse_par_umol_m2_s = inputs.diffuse_par_umol_m2_s;
    const top_direct_transmittance: f64 = 1.0;
    const top_diffuse_transmittance: f64 = 1.0;
    const top_scattered_shortwave_megajoules_m2_h: f64 = 0.0;
    const top_scattered_par_umol_m2_s: f64 = 0.0;
    const interception_stop_fraction: f64 = 0.0;
    return .{
        .inverse_cell_area_per_m2 = inverse_cell_area_per_m2,
        .inverse_azimuth_area_per_m2 = inverse_azimuth_area_per_m2,
        .working_diffuse_shortwave_megajoules_m2_h = working_diffuse_shortwave_megajoules_m2_h,
        .working_diffuse_par_umol_m2_s = working_diffuse_par_umol_m2_s,
        .top_direct_transmittance = top_direct_transmittance,
        .top_diffuse_transmittance = top_diffuse_transmittance,
        .top_scattered_shortwave_megajoules_m2_h = top_scattered_shortwave_megajoules_m2_h,
        .top_scattered_par_umol_m2_s = top_scattered_par_umol_m2_s,
        .interception_stop_fraction = interception_stop_fraction,
    };
}

test "canopy interception accumulators preserve source initialization" {
    const state = try initialize(.{
        .cell_area_m2 = 20,
        .diffuse_shortwave_megajoules_m2_h = 1.25,
        .diffuse_par_umol_m2_s = 300,
    });
    try std.testing.expectEqual(@as(f64, 0.05), state.inverse_cell_area_per_m2);
    try std.testing.expectEqual(@as(f64, 0.0125), state.inverse_azimuth_area_per_m2);
    try std.testing.expectEqual(
        @as(f64, 1.25),
        state.working_diffuse_shortwave_megajoules_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 300),
        state.working_diffuse_par_umol_m2_s,
    );
    try std.testing.expectEqual(@as(f64, 1), state.top_direct_transmittance);
    try std.testing.expectEqual(@as(f64, 1), state.top_diffuse_transmittance);
    try std.testing.expectEqual(
        @as(f64, 0),
        state.top_scattered_shortwave_megajoules_m2_h,
    );
    try std.testing.expectEqual(@as(f64, 0), state.top_scattered_par_umol_m2_s);
    try std.testing.expectEqual(@as(f64, 0), state.interception_stop_fraction);
}

test "invalid cell area fails before reciprocal" {
    try std.testing.expectError(
        error.InvalidCanopyInterceptionCellArea,
        initialize(.{
            .cell_area_m2 = 0,
            .diffuse_shortwave_megajoules_m2_h = 1,
            .diffuse_par_umol_m2_s = 1,
        }),
    );
}
