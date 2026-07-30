const std = @import("std");

pub const Inputs = struct {
    cell_count: usize,
    layer_capacity: usize,
    active: []const bool,
    liquid_water_volume_m3: []const f64,
    solid_snow_volume_m3: []const f64,
    air_filled_volume_m3: []const f64,
    litter_cover_fraction: []const f64,
    micropore_fraction: []const f64,
    macropore_fraction: []const f64,
    topsoil_micropore_air_capacity_m3: []const f64,
    topsoil_macropore_air_capacity_m3: []const f64,
    other_micropore_water_input_m3: []const f64,
    other_macropore_water_input_m3: []const f64,
    /// Physical fraction of the hour represented by this converged call.
    step_fraction: f64,
};

pub const Outputs = struct {
    /// Destination-layer indexed, matching WATSUB `FLQWM(M,L2,...)`.
    downward_water_flux_m3: []f64,
    litter_water_flux_m3: []f64,
    soil_micropore_water_flux_m3: []f64,
    soil_macropore_water_flux_m3: []f64,
};

/// Builds only the converged WATSUB snowmelt carriers used by TRNSFR(S).
/// It intentionally does not repeat the ecosystem model at snow substeps.
pub fn calculate(inputs: Inputs, outputs: Outputs) !void {
    try validate(inputs, outputs);
    @memset(outputs.downward_water_flux_m3, 0);
    @memset(outputs.litter_water_flux_m3, 0);
    @memset(outputs.soil_micropore_water_flux_m3, 0);
    @memset(outputs.soil_macropore_water_flux_m3, 0);

    for (0..inputs.cell_count) |cell| {
        var lowest_active: ?usize = null;
        for (0..inputs.layer_capacity) |layer| {
            const index = cell * inputs.layer_capacity + layer;
            if (inputs.active[index]) lowest_active = layer;
        }
        const bottom = lowest_active orelse continue;
        for (0..bottom) |layer| {
            const source = cell * inputs.layer_capacity + layer;
            const destination = source + 1;
            if (!inputs.active[source] or !inputs.active[destination]) continue;
            const releasable = releasableWater(inputs.liquid_water_volume_m3[source], inputs.solid_snow_volume_m3[source], inputs.step_fraction);
            // WATSUB limits FLWQM by the receiving snow-layer air-filled
            // volume represented by THETP2 in its interim state.
            outputs.downward_water_flux_m3[destination] = @min(releasable, inputs.air_filled_volume_m3[destination]);
        }

        const bottom_index = cell * inputs.layer_capacity + bottom;
        const releasable = releasableWater(inputs.liquid_water_volume_m3[bottom_index], inputs.solid_snow_volume_m3[bottom_index], inputs.step_fraction);
        const bare_delivery = releasable * (1 - inputs.litter_cover_fraction[cell]);
        const micro_capacity = @max(0, inputs.topsoil_micropore_air_capacity_m3[cell] * inputs.step_fraction - inputs.other_micropore_water_input_m3[cell]);
        const macro_capacity = @max(0, inputs.topsoil_macropore_air_capacity_m3[cell] * inputs.step_fraction - inputs.other_macropore_water_input_m3[cell]);
        const micro = @min(micro_capacity, bare_delivery * inputs.micropore_fraction[cell]);
        const macro_candidate = @min(macro_capacity, bare_delivery * inputs.macropore_fraction[cell]);
        // Preserve an exact nonnegative remainder despite the final ulp in
        // user fractions whose validated sum is numerically close to one.
        const macro = @min(macro_candidate, @max(0, releasable - micro));
        outputs.soil_micropore_water_flux_m3[cell] = micro;
        outputs.soil_macropore_water_flux_m3[cell] = macro;
        outputs.litter_water_flux_m3[cell] = @max(0, releasable - micro - macro);
    }
}

fn releasableWater(liquid_m3: f64, solid_snow_m3: f64, step_fraction: f64) f64 {
    return @max(0, liquid_m3 - 0.05 * solid_snow_m3) * step_fraction;
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    if (inputs.cell_count == 0 or inputs.layer_capacity == 0) return error.InvalidSnowMeltRoutingDimensions;
    const layers = try std.math.mul(usize, inputs.cell_count, inputs.layer_capacity);
    inline for (.{ inputs.active.len, inputs.liquid_water_volume_m3.len, inputs.solid_snow_volume_m3.len, inputs.air_filled_volume_m3.len, outputs.downward_water_flux_m3.len }) |length| if (length != layers) return error.InvalidSnowMeltRoutingDimensions;
    inline for (.{ inputs.litter_cover_fraction.len, inputs.micropore_fraction.len, inputs.macropore_fraction.len, inputs.topsoil_micropore_air_capacity_m3.len, inputs.topsoil_macropore_air_capacity_m3.len, inputs.other_micropore_water_input_m3.len, inputs.other_macropore_water_input_m3.len, outputs.litter_water_flux_m3.len, outputs.soil_micropore_water_flux_m3.len, outputs.soil_macropore_water_flux_m3.len }) |length| if (length != inputs.cell_count) return error.InvalidSnowMeltRoutingDimensions;
    if (!std.math.isFinite(inputs.step_fraction) or inputs.step_fraction <= 0 or inputs.step_fraction > 1) return error.InvalidSnowMeltStepFraction;
    inline for (.{ inputs.liquid_water_volume_m3, inputs.solid_snow_volume_m3, inputs.air_filled_volume_m3, inputs.topsoil_micropore_air_capacity_m3, inputs.topsoil_macropore_air_capacity_m3, inputs.other_micropore_water_input_m3, inputs.other_macropore_water_input_m3 }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowMeltRoutingState;
    for (inputs.litter_cover_fraction, inputs.micropore_fraction, inputs.macropore_fraction) |litter, micro, macro| {
        if (!std.math.isFinite(litter) or litter < 0 or litter > 1 or !std.math.isFinite(micro) or micro < 0 or micro > 1 or !std.math.isFinite(macro) or macro < 0 or macro > 1 or @abs(micro + macro - 1) > 1e-10) return error.InvalidSnowMeltRoutingFraction;
    }
}

test "WATSUB snow liquid retention and lowest-layer routing are conservative" {
    var downward = [_]f64{0} ** 2;
    var litter = [_]f64{0};
    var micro = [_]f64{0};
    var macro = [_]f64{0};
    try calculate(.{
        .cell_count = 1,
        .layer_capacity = 2,
        .active = &.{ true, true },
        .liquid_water_volume_m3 = &.{ 0.3, 0.5 },
        .solid_snow_volume_m3 = &.{ 2, 4 },
        .air_filled_volume_m3 = &.{ 0.1, 0.15 },
        .litter_cover_fraction = &.{0.25},
        .micropore_fraction = &.{0.6},
        .macropore_fraction = &.{0.4},
        .topsoil_micropore_air_capacity_m3 = &.{0.1},
        .topsoil_macropore_air_capacity_m3 = &.{0.2},
        .other_micropore_water_input_m3 = &.{0},
        .other_macropore_water_input_m3 = &.{0},
        .step_fraction = 1,
    }, .{ .downward_water_flux_m3 = &downward, .litter_water_flux_m3 = &litter, .soil_micropore_water_flux_m3 = &micro, .soil_macropore_water_flux_m3 = &macro });
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), downward[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), micro[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), macro[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.11), litter[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), micro[0] + macro[0] + litter[0], 1e-12);
}

test "snow routing rejects non-finite state before publishing outputs" {
    var downward = [_]f64{7};
    var litter = [_]f64{7};
    var micro = [_]f64{7};
    var macro = [_]f64{7};
    try std.testing.expectError(error.InvalidSnowMeltRoutingState, calculate(.{ .cell_count = 1, .layer_capacity = 1, .active = &.{true}, .liquid_water_volume_m3 = &.{std.math.nan(f64)}, .solid_snow_volume_m3 = &.{0}, .air_filled_volume_m3 = &.{0}, .litter_cover_fraction = &.{0}, .micropore_fraction = &.{1}, .macropore_fraction = &.{0}, .topsoil_micropore_air_capacity_m3 = &.{1}, .topsoil_macropore_air_capacity_m3 = &.{0}, .other_micropore_water_input_m3 = &.{0}, .other_macropore_water_input_m3 = &.{0}, .step_fraction = 1 }, .{ .downward_water_flux_m3 = &downward, .litter_water_flux_m3 = &litter, .soil_micropore_water_flux_m3 = &micro, .soil_macropore_water_flux_m3 = &macro }));
    try std.testing.expectEqual(@as(f64, 7), downward[0]);
}
