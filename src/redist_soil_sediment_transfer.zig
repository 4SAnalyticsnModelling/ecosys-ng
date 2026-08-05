const std = @import("std");

pub const Context = struct {
    source_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    source_layer_thickness_m: f64,
    destination_layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
    redistribution_fraction: f64, // FX; REDIST assigns FBO=FX.
    organic_fraction: f64, // FO is retained for traceability but is not used by 9594--9607.
};

pub const Layer = struct {
    acidity_ph: f64,
    sand_megagrams: f64,
    silt_megagrams: f64,
    clay_megagrams: f64,
    rock_megagrams: f64,
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

fn transferPool(source: f64, destination: f64, bulk_fraction: f64) !struct { source: f64, destination: f64 } {
    const transferred = bulk_fraction * source;
    const next_destination = destination + transferred;
    const next_source = source - transferred;
    if (!std.math.isFinite(transferred) or !std.math.isFinite(next_destination) or !std.math.isFinite(next_source))
        return error.NonFiniteSoilSedimentTransferResult;
    return .{ .source = next_source, .destination = next_destination };
}

/// Direct translation of REDIST 9594--9607 under the enclosing soil-material gates.
pub fn transfer(context: Context, state: *State) !void {
    if (!finiteStruct(context) or !finiteStruct(state.source) or !finiteStruct(state.destination) or
        context.redistribution_fraction < 0 or context.redistribution_fraction > 1 or
        context.organic_fraction < 0 or context.organic_fraction > 1 or context.minimum_layer_thickness_m < 0)
        return error.InvalidSoilSedimentTransferInput;
    if (context.source_layer == 0 or context.source_bulk_density_megagrams_m3 <= 0 or
        context.destination_bulk_density_megagrams_m3 <= 0 or
        context.source_layer_thickness_m <= context.minimum_layer_thickness_m or
        context.destination_layer_thickness_m <= context.minimum_layer_thickness_m) return;

    const bulk_fraction = context.redistribution_fraction; // FBO=FX; FO is not used here.
    var next = state.*;
    next.destination.acidity_ph = (1.0 - context.redistribution_fraction) * state.destination.acidity_ph +
        context.redistribution_fraction * state.source.acidity_ph;
    if (!std.math.isFinite(next.destination.acidity_ph)) return error.NonFiniteSoilSedimentTransferResult;
    inline for (.{ "sand_megagrams", "silt_megagrams", "clay_megagrams", "rock_megagrams" }) |field_name| {
        const result = try transferPool(@field(state.source, field_name), @field(state.destination, field_name), bulk_fraction);
        @field(next.source, field_name) = result.source;
        @field(next.destination, field_name) = result.destination;
    }
    state.* = next;
}

test "REDIST soil sediment mixes pH then transfers exact pool order" {
    var state = State{
        .source = .{ .acidity_ph = 4, .sand_megagrams = 2, .silt_megagrams = 4, .clay_megagrams = 6, .rock_megagrams = 8 },
        .destination = .{ .acidity_ph = 8, .sand_megagrams = 10, .silt_megagrams = 20, .clay_megagrams = 30, .rock_megagrams = 40 },
    };
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .source_layer_thickness_m = 1, .destination_layer_thickness_m = 1, .minimum_layer_thickness_m = 0.01, .redistribution_fraction = 0.25, .organic_fraction = 0.9 }, &state);
    try std.testing.expectEqual(@as(f64, 7), state.destination.acidity_ph);
    try std.testing.expectEqual(@as(f64, 1.5), state.source.sand_megagrams);
    try std.testing.expectEqual(@as(f64, 10.5), state.destination.sand_megagrams);
    try std.testing.expectEqual(@as(f64, 6), state.source.rock_megagrams);
    try std.testing.expectEqual(@as(f64, 42), state.destination.rock_megagrams);
}

test "REDIST soil sediment conserves every extensive solid pool" {
    var state = State{
        .source = .{ .acidity_ph = 5, .sand_megagrams = 2, .silt_megagrams = 4, .clay_megagrams = 6, .rock_megagrams = 8 },
        .destination = .{ .acidity_ph = 7, .sand_megagrams = 10, .silt_megagrams = 20, .clay_megagrams = 30, .rock_megagrams = 40 },
    };
    const before = state;
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .source_layer_thickness_m = 1, .destination_layer_thickness_m = 1, .minimum_layer_thickness_m = 0, .redistribution_fraction = 0.4, .organic_fraction = 0 }, &state);
    inline for (.{ "sand_megagrams", "silt_megagrams", "clay_megagrams", "rock_megagrams" }) |field_name|
        try std.testing.expectEqual(@field(before.source, field_name) + @field(before.destination, field_name), @field(state.source, field_name) + @field(state.destination, field_name));
}

test "REDIST soil sediment gate and validation leave state atomic" {
    var state = State{
        .source = .{ .acidity_ph = 5, .sand_megagrams = 2, .silt_megagrams = 4, .clay_megagrams = 6, .rock_megagrams = 8 },
        .destination = .{ .acidity_ph = 7, .sand_megagrams = 10, .silt_megagrams = 20, .clay_megagrams = 30, .rock_megagrams = 40 },
    };
    const before = state;
    try transfer(.{ .source_layer = 0, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .source_layer_thickness_m = 1, .destination_layer_thickness_m = 1, .minimum_layer_thickness_m = 0, .redistribution_fraction = 1, .organic_fraction = 0 }, &state);
    try std.testing.expectEqualDeep(before, state);
    try std.testing.expectError(error.InvalidSoilSedimentTransferInput, transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .source_layer_thickness_m = 1, .destination_layer_thickness_m = 1, .minimum_layer_thickness_m = 0, .redistribution_fraction = std.math.nan(f64), .organic_fraction = 0 }, &state));
    try std.testing.expectEqualDeep(before, state);
}
