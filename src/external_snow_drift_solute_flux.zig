const std = @import("std");

pub const species_count = 41;
pub const Side = enum { forward, reverse };
pub const Outcome = enum { zero_flux, loss_accumulated };

pub const Inputs = struct {
    side: Side,
    boundary_exchange_enabled: bool,
    boundary_charge_fraction: f64,
    total_snow_transfer_m3_per_step: f64,
    directional_snow_transfer_m3_per_step: f64,
    minimum_snow_transfer_m3_per_step: f64,
    /// 33 salt/complex species without H4SiO4, followed by eight P species.
    source_snow_solute_amount_per_step: []const f64,
};

pub const State = struct {
    boundary_flux_amount_per_step: []f64,
    accumulated_boundary_flux_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 7473--7707.
/// Boundary charge is a gate only. Forward loss is positive transfer and
/// reverse loss is negative; runon is explicitly assigned zero.
pub fn route(inputs: Inputs, state: State) !Outcome {
    try validateTopologyAndControls(inputs, state);
    const transport_enabled = inputs.boundary_exchange_enabled and inputs.boundary_charge_fraction != 0 and
        inputs.total_snow_transfer_m3_per_step > inputs.minimum_snow_transfer_m3_per_step;
    const outward_loss = switch (inputs.side) {
        .forward => inputs.directional_snow_transfer_m3_per_step > inputs.minimum_snow_transfer_m3_per_step,
        .reverse => inputs.directional_snow_transfer_m3_per_step < inputs.minimum_snow_transfer_m3_per_step,
    };
    if (!transport_enabled or !outward_loss) {
        @memset(state.boundary_flux_amount_per_step, 0);
        return .zero_flux;
    }

    try validateActiveFluxInputs(inputs, state);

    const directional_fraction = inputs.directional_snow_transfer_m3_per_step /
        inputs.total_snow_transfer_m3_per_step;
    for (0..species_count) |species| {
        const flux = inputs.source_snow_solute_amount_per_step[species] * directional_fraction;
        if (!std.math.isFinite(flux) or
            !std.math.isFinite(state.accumulated_boundary_flux_amount_per_step[species] + flux))
            return error.NonFiniteExternalSnowDriftSoluteResult;
    }
    for (0..species_count) |species| {
        const flux = inputs.source_snow_solute_amount_per_step[species] * directional_fraction;
        state.boundary_flux_amount_per_step[species] = flux;
        state.accumulated_boundary_flux_amount_per_step[species] =
            state.accumulated_boundary_flux_amount_per_step[species] + flux;
    }
    return .loss_accumulated;
}

fn validateTopologyAndControls(inputs: Inputs, state: State) !void {
    if (inputs.source_snow_solute_amount_per_step.len != species_count or
        state.boundary_flux_amount_per_step.len != species_count or
        state.accumulated_boundary_flux_amount_per_step.len != species_count)
        return error.ExternalSnowDriftSoluteDimensionMismatch;
    inline for (.{ inputs.boundary_charge_fraction, inputs.total_snow_transfer_m3_per_step, inputs.directional_snow_transfer_m3_per_step, inputs.minimum_snow_transfer_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSnowDriftSoluteInput;
    if (inputs.total_snow_transfer_m3_per_step < 0 or inputs.minimum_snow_transfer_m3_per_step < 0)
        return error.InvalidExternalSnowDriftSoluteInput;
}

fn validateActiveFluxInputs(inputs: Inputs, state: State) !void {
    for (inputs.source_snow_solute_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSnowDriftSoluteInput;
    for (state.accumulated_boundary_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSnowDriftSoluteInput;
}

fn fixture(source: []const f64, side: Side, directional_transfer: f64) Inputs {
    return .{ .side = side, .boundary_exchange_enabled = true, .boundary_charge_fraction = 0.25, .total_snow_transfer_m3_per_step = 4, .directional_snow_transfer_m3_per_step = directional_transfer, .minimum_snow_transfer_m3_per_step = 0.01, .source_snow_solute_amount_per_step = source };
}

test "TRNSFRS forward snow-drift loss partitions and accumulates 41 species" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    const outcome = try route(fixture(&source, .forward, 2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.loss_accumulated, outcome);
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 5), accumulated[40]);
}

test "TRNSFRS reverse snow-drift loss retains negative sign" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    const outcome = try route(fixture(&source, .reverse, -2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.loss_accumulated, outcome);
    try std.testing.expectEqual(@as(f64, -4), flux[0]);
    try std.testing.expectEqual(@as(f64, -3), accumulated[40]);
}

test "TRNSFRS snow runon zeros instantaneous flux without accumulation" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    const outcome = try route(fixture(&source, .forward, -2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.zero_flux, outcome);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{1} ** species_count), &accumulated);
}

test "TRNSFRS disabled charge gate zeros all 41 instantaneous species" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    var inputs = fixture(&source, .forward, 2);
    inputs.boundary_charge_fraction = 0;
    inputs.source_snow_solute_amount_per_step = &([_]f64{std.math.nan(f64)} ** species_count);
    accumulated[40] = std.math.inf(f64);
    try std.testing.expectEqual(Outcome.zero_flux, try route(inputs, .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated }));
    try std.testing.expectEqual(@as(f64, 0), flux[40]);
}

test "late invalid accumulator keeps snow flux atomic" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    accumulated[40] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteExternalSnowDriftSoluteInput, route(fixture(&source, .forward, 2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated }));
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** species_count), &flux);
}
