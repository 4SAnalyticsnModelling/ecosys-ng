const std = @import("std");

/// Runtime-sized EXTRACT `TTRAN/CTRAN`, `TVOLWP/TVOLWC`, and
/// `TEVAPP/TEVAPC` publication. Plant arrays are cell-major then species-major.
/// Water fluxes retain the source sign: atmospheric losses are negative.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    hourly_water_source_m3_per_h_by_plant: []f64,
    cumulative_water_source_m3_by_plant: []f64,
    internal_water_m3_by_cell: []f64,
    surface_water_m3_by_cell: []f64,
    evapotranspiration_m3_per_h_by_cell: []f64,
    evaporation_m3_per_h_by_cell: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        species_count: usize,
    ) !State {
        if (cell_count == 0 or species_count == 0)
            return error.InvalidPlantWaterPublicationDimensions;
        const plant_count = try std.math.mul(usize, cell_count, species_count);
        const hourly = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(hourly);
        const cumulative = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(cumulative);
        const internal = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(internal);
        const surface = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(surface);
        const evapotranspiration = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(evapotranspiration);
        const evaporation = try allocator.alloc(f64, cell_count);
        @memset(hourly, 0);
        @memset(cumulative, 0);
        @memset(internal, 0);
        @memset(surface, 0);
        @memset(evapotranspiration, 0);
        @memset(evaporation, 0);
        return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .species_count = species_count,
            .hourly_water_source_m3_per_h_by_plant = hourly,
            .cumulative_water_source_m3_by_plant = cumulative,
            .internal_water_m3_by_cell = internal,
            .surface_water_m3_by_cell = surface,
            .evapotranspiration_m3_per_h_by_cell = evapotranspiration,
            .evaporation_m3_per_h_by_cell = evaporation,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.evaporation_m3_per_h_by_cell);
        self.allocator.free(self.evapotranspiration_m3_per_h_by_cell);
        self.allocator.free(self.surface_water_m3_by_cell);
        self.allocator.free(self.internal_water_m3_by_cell);
        self.allocator.free(self.cumulative_water_source_m3_by_plant);
        self.allocator.free(self.hourly_water_source_m3_per_h_by_plant);
        self.* = undefined;
    }

    /// DAY/HOUR1 reset of the legacy daily `CTRAN` carrier.
    pub fn resetDaily(self: *State) void {
        @memset(self.cumulative_water_source_m3_by_plant, 0);
    }
};

pub const Inputs = struct {
    cell_area_m2: []const f64,
    internal_water_depth_m_per_m2_by_plant: []const f64,
    living_surface_water_m3_by_plant: []const f64,
    standing_dead_surface_water_m3_by_plant: []const f64,
    transpiration_m3_per_h_by_plant: []const f64,
    living_evaporation_m3_per_h_by_plant: []const f64,
    standing_dead_evaporation_m3_per_h_by_plant: []const f64,
};

const Candidate = struct {
    water_source_m3_per_h: f64,
    cumulative_water_source_m3: f64,
    internal_water_m3: f64,
    surface_water_m3: f64,
    evaporation_m3_per_h: f64,
};

/// Exact EXTRACT lines 637–642 publication. The full domain is preflighted
/// before any hourly, cumulative, or cell aggregate is changed.
pub fn refresh(state: *State, inputs: Inputs) !void {
    const plant_count = try std.math.mul(
        usize,
        state.cell_count,
        state.species_count,
    );
    if (inputs.cell_area_m2.len != state.cell_count)
        return error.InvalidPlantWaterPublicationDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields[1..]) |field|
        if (@field(inputs, field.name).len != plant_count)
            return error.InvalidPlantWaterPublicationDimensions;

    for (0..plant_count) |plant| _ = try candidate(state, inputs, plant);

    @memset(state.internal_water_m3_by_cell, 0);
    @memset(state.surface_water_m3_by_cell, 0);
    @memset(state.evapotranspiration_m3_per_h_by_cell, 0);
    @memset(state.evaporation_m3_per_h_by_cell, 0);
    for (0..plant_count) |plant| {
        const next = candidate(state, inputs, plant) catch unreachable;
        const cell = plant / state.species_count;
        state.hourly_water_source_m3_per_h_by_plant[plant] =
            next.water_source_m3_per_h;
        state.cumulative_water_source_m3_by_plant[plant] =
            next.cumulative_water_source_m3;
        state.internal_water_m3_by_cell[cell] += next.internal_water_m3;
        state.surface_water_m3_by_cell[cell] += next.surface_water_m3;
        state.evapotranspiration_m3_per_h_by_cell[cell] +=
            next.water_source_m3_per_h;
        state.evaporation_m3_per_h_by_cell[cell] += next.evaporation_m3_per_h;
    }
}

fn candidate(state: *const State, inputs: Inputs, plant: usize) !Candidate {
    const cell = plant / state.species_count;
    const area = inputs.cell_area_m2[cell];
    inline for (@typeInfo(Inputs).@"struct".fields[1..]) |field| {
        const value = @field(inputs, field.name)[plant];
        if (!std.math.isFinite(value))
            return error.NonFinitePlantWaterPublicationInput;
    }
    if (!std.math.isFinite(area) or area <= 0 or
        inputs.internal_water_depth_m_per_m2_by_plant[plant] < 0 or
        inputs.living_surface_water_m3_by_plant[plant] < 0 or
        inputs.standing_dead_surface_water_m3_by_plant[plant] < 0)
        return error.InvalidPlantWaterPublicationInput;
    const evaporation =
        inputs.living_evaporation_m3_per_h_by_plant[plant] +
        inputs.standing_dead_evaporation_m3_per_h_by_plant[plant];
    const water_source =
        inputs.transpiration_m3_per_h_by_plant[plant] + evaporation;
    const result: Candidate = .{
        .water_source_m3_per_h = water_source,
        .cumulative_water_source_m3 = state.cumulative_water_source_m3_by_plant[plant] + water_source,
        .internal_water_m3 = inputs.internal_water_depth_m_per_m2_by_plant[plant] * area,
        .surface_water_m3 = inputs.living_surface_water_m3_by_plant[plant] +
            inputs.standing_dead_surface_water_m3_by_plant[plant],
        .evaporation_m3_per_h = evaporation,
    };
    inline for (@typeInfo(Candidate).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFinitePlantWaterPublication;
    return result;
}

test "EXTRACT water publication preserves source signs and runtime species order" {
    var state = try State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    try refresh(&state, .{
        .cell_area_m2 = &.{ 10, 20 },
        .internal_water_depth_m_per_m2_by_plant = &.{ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06 },
        .living_surface_water_m3_by_plant = &.{ 1, 2, 3, 4, 5, 6 },
        .standing_dead_surface_water_m3_by_plant = &.{ 6, 5, 4, 3, 2, 1 },
        .transpiration_m3_per_h_by_plant = &.{ -1, -2, -3, -4, -5, -6 },
        .living_evaporation_m3_per_h_by_plant = &.{ -0.1, -0.2, -0.3, -0.4, -0.5, -0.6 },
        .standing_dead_evaporation_m3_per_h_by_plant = &.{ -0.01, -0.02, -0.03, -0.04, -0.05, -0.06 },
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, -6.66),
        state.hourly_water_source_m3_per_h_by_plant[5],
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.6),
        state.internal_water_m3_by_cell[0],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 21.0),
        state.surface_water_m3_by_cell[1],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, -16.65),
        state.evapotranspiration_m3_per_h_by_cell[1],
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, -1.65),
        state.evaporation_m3_per_h_by_cell[1],
        1e-14,
    );
    try refresh(&state, .{
        .cell_area_m2 = &.{ 10, 20 },
        .internal_water_depth_m_per_m2_by_plant = &.{ 0, 0, 0, 0, 0, 0 },
        .living_surface_water_m3_by_plant = &.{ 0, 0, 0, 0, 0, 0 },
        .standing_dead_surface_water_m3_by_plant = &.{ 0, 0, 0, 0, 0, 0 },
        .transpiration_m3_per_h_by_plant = &.{ -1, -2, -3, -4, -5, -6 },
        .living_evaporation_m3_per_h_by_plant = &.{ 0, 0, 0, 0, 0, 0 },
        .standing_dead_evaporation_m3_per_h_by_plant = &.{ 0, 0, 0, 0, 0, 0 },
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, -12.66),
        state.cumulative_water_source_m3_by_plant[5],
        1e-14,
    );
    state.resetDaily();
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 0, 0, 0, 0, 0 },
        state.cumulative_water_source_m3_by_plant,
    );
}

test "late invalid plant leaves the complete water publication unchanged" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    @memset(state.hourly_water_source_m3_per_h_by_plant, 7);
    @memset(state.cumulative_water_source_m3_by_plant, 8);
    @memset(state.internal_water_m3_by_cell, 9);
    const invalid = [_]f64{ 0, std.math.nan(f64) };
    try std.testing.expectError(
        error.NonFinitePlantWaterPublicationInput,
        refresh(&state, .{
            .cell_area_m2 = &.{10},
            .internal_water_depth_m_per_m2_by_plant = &invalid,
            .living_surface_water_m3_by_plant = &.{ 0, 0 },
            .standing_dead_surface_water_m3_by_plant = &.{ 0, 0 },
            .transpiration_m3_per_h_by_plant = &.{ 0, 0 },
            .living_evaporation_m3_per_h_by_plant = &.{ 0, 0 },
            .standing_dead_evaporation_m3_per_h_by_plant = &.{ 0, 0 },
        }),
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 7, 7 },
        state.hourly_water_source_m3_per_h_by_plant,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 8, 8 },
        state.cumulative_water_source_m3_by_plant,
    );
    try std.testing.expectEqualSlices(f64, &.{9}, state.internal_water_m3_by_cell);
}
