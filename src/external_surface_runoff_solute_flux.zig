const std = @import("std");

pub const species_count = 42;
pub const Direction = enum { east_west, north_south, vertical };
pub const Side = enum { forward, reverse };
pub const Outcome = enum { not_surface_horizontal, zero_flux, loss_accumulated };

pub const Inputs = struct {
    is_surface_layer: bool,
    direction: Direction,
    side: Side,
    boundary_exchange_enabled: bool,
    boundary_charge_fraction: f64,
    total_runoff_m3_per_step: f64,
    directional_runoff_m3_per_step: f64,
    minimum_runoff_m3_per_step: f64,
    source_runoff_solute_amount_per_step: []const f64,
};

pub const State = struct {
    boundary_flux_amount_per_step: []f64,
    accumulated_boundary_flux_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 7199--7466.
/// Boundary charge controls enablement only; it does not scale the flux.
/// Forward loss requires positive directional runoff, reverse loss negative.
pub fn route(inputs: Inputs, state: State) !Outcome {
    if (!inputs.is_surface_layer or inputs.direction == .vertical) return .not_surface_horizontal;
    try validateTopologyAndControls(inputs, state);

    const transport_enabled = inputs.boundary_exchange_enabled and inputs.boundary_charge_fraction != 0 and
        inputs.total_runoff_m3_per_step > inputs.minimum_runoff_m3_per_step;
    const outward_loss = switch (inputs.side) {
        .forward => inputs.directional_runoff_m3_per_step > inputs.minimum_runoff_m3_per_step,
        .reverse => inputs.directional_runoff_m3_per_step < inputs.minimum_runoff_m3_per_step,
    };
    if (!transport_enabled or !outward_loss) {
        @memset(state.boundary_flux_amount_per_step, 0);
        return .zero_flux;
    }

    try validateActiveFluxInputs(inputs, state);

    const directional_fraction = inputs.directional_runoff_m3_per_step / inputs.total_runoff_m3_per_step;
    for (0..species_count) |species| {
        const flux = inputs.source_runoff_solute_amount_per_step[species] * directional_fraction;
        if (!std.math.isFinite(flux) or
            !std.math.isFinite(state.accumulated_boundary_flux_amount_per_step[species] + flux))
            return error.NonFiniteExternalSurfaceRunoffSoluteResult;
    }
    for (0..species_count) |species| {
        const flux = inputs.source_runoff_solute_amount_per_step[species] * directional_fraction;
        state.boundary_flux_amount_per_step[species] = flux;
        state.accumulated_boundary_flux_amount_per_step[species] =
            state.accumulated_boundary_flux_amount_per_step[species] + flux;
    }
    return .loss_accumulated;
}

fn validateTopologyAndControls(inputs: Inputs, state: State) !void {
    if (inputs.source_runoff_solute_amount_per_step.len != species_count or
        state.boundary_flux_amount_per_step.len != species_count or
        state.accumulated_boundary_flux_amount_per_step.len != species_count)
        return error.ExternalSurfaceRunoffSoluteDimensionMismatch;
    inline for (.{ inputs.boundary_charge_fraction, inputs.total_runoff_m3_per_step, inputs.directional_runoff_m3_per_step, inputs.minimum_runoff_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSurfaceRunoffSoluteInput;
    if (inputs.total_runoff_m3_per_step < 0 or inputs.minimum_runoff_m3_per_step < 0)
        return error.InvalidExternalSurfaceRunoffSoluteInput;
}

fn validateActiveFluxInputs(inputs: Inputs, state: State) !void {
    for (inputs.source_runoff_solute_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSurfaceRunoffSoluteInput;
    for (state.accumulated_boundary_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSurfaceRunoffSoluteInput;
}

fn fixture(source: []const f64, side: Side, directional_runoff: f64) Inputs {
    return .{ .is_surface_layer = true, .direction = .east_west, .side = side, .boundary_exchange_enabled = true, .boundary_charge_fraction = 0.25, .total_runoff_m3_per_step = 4, .directional_runoff_m3_per_step = directional_runoff, .minimum_runoff_m3_per_step = 0.01, .source_runoff_solute_amount_per_step = source };
}

test "TRNSFRS forward external runoff loss partitions and accumulates 42 species" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    const outcome = try route(fixture(&source, .forward, 2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.loss_accumulated, outcome);
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 5), accumulated[41]);
}

test "TRNSFRS reverse external runoff loss retains negative sign" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    const outcome = try route(fixture(&source, .reverse, -2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.loss_accumulated, outcome);
    try std.testing.expectEqual(@as(f64, -4), flux[0]);
    try std.testing.expectEqual(@as(f64, -3), accumulated[41]);
}

test "TRNSFRS runon and disabled boundary zero instantaneous flux without accumulation" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    const outcome = try route(fixture(&source, .forward, -2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.zero_flux, outcome);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{1} ** species_count), &accumulated);
}

test "TRNSFRS nonsurface or vertical path leaves boundary state untouched" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    var inputs = fixture(&source, .forward, 2);
    inputs.direction = .vertical;
    const outcome = try route(inputs, .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.not_surface_horizontal, outcome);
    try std.testing.expectEqual(@as(f64, 9), flux[0]);
}

test "TRNSFRS nonsurface path leaves invalid dormant state unexamined" {
    const empty = [_]f64{};
    var empty_flux = [_]f64{};
    var empty_accumulated = [_]f64{};
    var inputs = fixture(&empty, .forward, std.math.nan(f64));
    inputs.is_surface_layer = false;
    inputs.total_runoff_m3_per_step = std.math.nan(f64);
    const outcome = try route(inputs, .{ .boundary_flux_amount_per_step = &empty_flux, .accumulated_boundary_flux_amount_per_step = &empty_accumulated });
    try std.testing.expectEqual(Outcome.not_surface_horizontal, outcome);
}

test "TRNSFRS disabled loss gate does not inspect dormant solute state" {
    var source = [_]f64{8} ** species_count;
    var flux = [_]f64{std.math.nan(f64)} ** species_count;
    var accumulated = [_]f64{std.math.inf(f64)} ** species_count;
    source[41] = std.math.nan(f64);
    var inputs = fixture(&source, .forward, 2);
    inputs.boundary_exchange_enabled = false;
    const outcome = try route(inputs, .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated });
    try std.testing.expectEqual(Outcome.zero_flux, outcome);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
    try std.testing.expect(std.math.isInf(accumulated[0]));
}

test "late invalid accumulator keeps instantaneous flux atomic" {
    const source = [_]f64{8} ** species_count;
    var flux = [_]f64{9} ** species_count;
    var accumulated = [_]f64{1} ** species_count;
    accumulated[41] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteExternalSurfaceRunoffSoluteInput, route(fixture(&source, .forward, 2), .{ .boundary_flux_amount_per_step = &flux, .accumulated_boundary_flux_amount_per_step = &accumulated }));
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** species_count), &flux);
}
