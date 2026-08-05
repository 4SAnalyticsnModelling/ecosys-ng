const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const AcceptedDirectLiquid = struct {
    rain_m3: f64,
    irrigation_m3: f64,
};

/// Allocates only the liquid volume accepted by the surface-water router to
/// its weather and irrigation sources. Unaccepted water must not carry solute
/// into either storage or the external-input ledger.
pub fn partitionAcceptedDirectLiquid(
    accepted_direct_m3: f64,
    available_rain_m3: f64,
    available_irrigation_m3: f64,
    tolerance_m3: f64,
) !AcceptedDirectLiquid {
    inline for (.{ accepted_direct_m3, available_rain_m3, available_irrigation_m3, tolerance_m3 }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidAtmosphericSoluteInputVolume;
    const available_m3 = available_rain_m3 + available_irrigation_m3;
    if (!std.math.isFinite(available_m3) or accepted_direct_m3 > available_m3 + tolerance_m3)
        return error.AtmosphericSoluteInputExceedsAvailableWater;
    if (accepted_direct_m3 == 0 or available_m3 == 0) return .{ .rain_m3 = 0, .irrigation_m3 = 0 };
    const accepted_m3 = @min(accepted_direct_m3, available_m3);
    const rain_m3 = accepted_m3 * available_rain_m3 / available_m3;
    return .{ .rain_m3 = rain_m3, .irrigation_m3 = accepted_m3 - rain_m3 };
}

/// Runtime-sized accepted precipitation and irrigation chemistry before it
/// branches into snow, litter, non-band soil, and band soil destinations.
/// Values retain each carrier's tracked-element grams.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    daily_input_g: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.EmptyAtmosphericSoluteInputGrid;
        const values = try allocator.alloc(
            f64,
            try std.math.mul(usize, cell_count, snow.species_count),
        );
        @memset(values, 0);
        return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .daily_input_g = values,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.daily_input_g);
        self.* = undefined;
    }

    pub fn resetDaily(self: *State) void {
        @memset(self.daily_input_g, 0);
    }

    /// Adds one fully assembled hour atomically. `snow_input_g` is the branch
    /// retained by snow; `direct_input` holds the mutually exclusive branch
    /// routed directly to litter and topsoil.
    pub fn accumulateAcceptedHour(
        self: *State,
        snow_input_g: []const f64,
        direct_input: []const snow.SurfaceDischarge,
    ) !void {
        if (snow_input_g.len != self.daily_input_g.len or
            direct_input.len != self.cell_count)
            return error.AtmosphericSoluteInputDimensionMismatch;

        for (0..self.cell_count) |cell| {
            const first = cell * snow.species_count;
            for (0..snow.species_count) |species| {
                const direct = direct_input[cell];
                const increment =
                    snow_input_g[first + species] +
                    direct.litter_g[species] +
                    direct.soil_nonband_g[species] +
                    direct.soil_band_g[species];
                const next = self.daily_input_g[first + species] + increment;
                if (!std.math.isFinite(increment) or increment < 0)
                    return error.InvalidAtmosphericSoluteInput;
                if (!std.math.isFinite(next))
                    return error.AtmosphericSoluteInputOverflow;
            }
        }
        for (0..self.cell_count) |cell| {
            const first = cell * snow.species_count;
            for (0..snow.species_count) |species| {
                const direct = direct_input[cell];
                self.daily_input_g[first + species] +=
                    snow_input_g[first + species] +
                    direct.litter_g[species] +
                    direct.soil_nonband_g[species] +
                    direct.soil_band_g[species];
            }
        }
    }

    pub fn speciesInputG(
        self: State,
        cell: usize,
        species: snow.Species,
    ) !f64 {
        if (cell >= self.cell_count)
            return error.AtmosphericSoluteInputCellOutOfBounds;
        return self.daily_input_g[
            cell * snow.species_count + @intFromEnum(species)
        ];
    }
};

test "accepted atmospheric input recombines snow litter and both soil zones" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var snow_input = [_]f64{0} ** (2 * snow.species_count);
    var direct = [_]snow.SurfaceDischarge{ .{}, .{} };
    const species = @intFromEnum(snow.Species.ammonium_nitrogen);
    snow_input[species] = 1;
    direct[0].litter_g[species] = 2;
    direct[0].soil_nonband_g[species] = 3;
    direct[0].soil_band_g[species] = 4;
    direct[1].soil_band_g[species] = 5;
    try state.accumulateAcceptedHour(&snow_input, &direct);
    try std.testing.expectEqual(
        @as(f64, 10),
        try state.speciesInputG(0, .ammonium_nitrogen),
    );
    try std.testing.expectEqual(
        @as(f64, 5),
        try state.speciesInputG(1, .ammonium_nitrogen),
    );
}

test "direct chemistry follows accepted water rather than nominal precipitation" {
    const accepted = try partitionAcceptedDirectLiquid(8, 9, 1, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7.2), accepted.rain_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), accepted.irrigation_m3, 1e-15);
    try std.testing.expectError(
        error.AtmosphericSoluteInputExceedsAvailableWater,
        partitionAcceptedDirectLiquid(11, 9, 1, 1e-12),
    );
}

test "invalid late atmospheric carrier leaves the daily ledger unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.daily_input_g[0] = 7;
    var snow_input = [_]f64{0} ** (2 * snow.species_count);
    snow_input[snow.species_count + 3] = std.math.nan(f64);
    const direct = [_]snow.SurfaceDischarge{ .{}, .{} };
    try std.testing.expectError(
        error.InvalidAtmosphericSoluteInput,
        state.accumulateAcceptedHour(&snow_input, &direct),
    );
    try std.testing.expectEqual(@as(f64, 7), state.daily_input_g[0]);
}
