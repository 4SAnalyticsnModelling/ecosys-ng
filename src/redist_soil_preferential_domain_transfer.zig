const std = @import("std");

pub const Context = struct {
    source_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    redistribution_fraction: f64, // FX; REDIST assigns FHO=FX.
};

pub const Layer = struct {
    preferential_fraction: f64, // FHOL, dimensionless.
    water_volume_m3: f64, // VOLWH
    ice_volume_m3: f64, // VOLIH
    air_volume_m3: f64, // VOLAH
};

pub const State = struct {
    source: Layer,
    destination: Layer,
};

fn finiteStruct(value: anytype) bool {
    inline for (std.meta.fields(@TypeOf(value))) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(value, field.name))) return false;
    }
    return true;
}

fn transferVolume(source: f64, destination: f64, fraction: f64) !struct { source: f64, destination: f64 } {
    const moved = fraction * source;
    const next_destination = destination + moved;
    const next_source = source - moved;
    if (!std.math.isFinite(moved) or !std.math.isFinite(next_destination) or !std.math.isFinite(next_source))
        return error.NonFiniteSoilPreferentialDomainTransferResult;
    return .{ .source = next_source, .destination = next_destination };
}

/// Direct translation of REDIST 9611--9623 under its enclosing soil gates.
pub fn transfer(context: Context, state: *State) !void {
    if (!finiteStruct(context) or !finiteStruct(state.source) or !finiteStruct(state.destination) or
        context.redistribution_fraction < 0 or context.redistribution_fraction > 1)
        return error.InvalidSoilPreferentialDomainTransferInput;
    if (context.source_layer == 0 or context.source_bulk_density_megagrams_m3 <= 0 or
        context.destination_bulk_density_megagrams_m3 <= 0 or state.destination.preferential_fraction <= 0 or
        state.source.preferential_fraction <= 0) return;

    const preferential_fraction = context.redistribution_fraction; // FHO=FX.
    var next = state.*;
    next.destination.preferential_fraction = (1.0 - context.redistribution_fraction) *
        state.destination.preferential_fraction + context.redistribution_fraction * state.source.preferential_fraction;
    if (!std.math.isFinite(next.destination.preferential_fraction))
        return error.NonFiniteSoilPreferentialDomainTransferResult;
    inline for (.{ "water_volume_m3", "ice_volume_m3", "air_volume_m3" }) |field_name| {
        const result = try transferVolume(@field(state.source, field_name), @field(state.destination, field_name), preferential_fraction);
        @field(next.source, field_name) = result.source;
        @field(next.destination, field_name) = result.destination;
    }
    state.* = next;
}

test "REDIST preferential domain mixes FHOL then transfers water ice air" {
    var state = State{
        .source = .{ .preferential_fraction = 0.4, .water_volume_m3 = 2, .ice_volume_m3 = 4, .air_volume_m3 = 6 },
        .destination = .{ .preferential_fraction = 0.2, .water_volume_m3 = 10, .ice_volume_m3 = 20, .air_volume_m3 = 30 },
    };
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = 0.25 }, &state);
    try std.testing.expectEqual(@as(f64, 0.25), state.destination.preferential_fraction);
    try std.testing.expectEqual(@as(f64, 1.5), state.source.water_volume_m3);
    try std.testing.expectEqual(@as(f64, 10.5), state.destination.water_volume_m3);
    try std.testing.expectEqual(@as(f64, 3), state.source.ice_volume_m3);
    try std.testing.expectEqual(@as(f64, 21), state.destination.ice_volume_m3);
    try std.testing.expectEqual(@as(f64, 4.5), state.source.air_volume_m3);
    try std.testing.expectEqual(@as(f64, 31.5), state.destination.air_volume_m3);
}

test "REDIST preferential-domain transfer conserves each volume" {
    var state = State{
        .source = .{ .preferential_fraction = 0.4, .water_volume_m3 = 2, .ice_volume_m3 = 4, .air_volume_m3 = 6 },
        .destination = .{ .preferential_fraction = 0.2, .water_volume_m3 = 10, .ice_volume_m3 = 20, .air_volume_m3 = 30 },
    };
    const before = state;
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = 0.4 }, &state);
    inline for (.{ "water_volume_m3", "ice_volume_m3", "air_volume_m3" }) |field_name|
        try std.testing.expectEqual(@field(before.source, field_name) + @field(before.destination, field_name), @field(state.source, field_name) + @field(state.destination, field_name));
}

test "REDIST preferential-domain gates and invalid input are atomic" {
    var state = State{
        .source = .{ .preferential_fraction = 0, .water_volume_m3 = 2, .ice_volume_m3 = 4, .air_volume_m3 = 6 },
        .destination = .{ .preferential_fraction = 0.2, .water_volume_m3 = 10, .ice_volume_m3 = 20, .air_volume_m3 = 30 },
    };
    const before = state;
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = 1 }, &state);
    try std.testing.expectEqualDeep(before, state);
    try std.testing.expectError(error.InvalidSoilPreferentialDomainTransferInput, transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = std.math.nan(f64) }, &state));
    try std.testing.expectEqualDeep(before, state);
}
