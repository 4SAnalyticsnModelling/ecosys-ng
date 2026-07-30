const std = @import("std");
const LayerState = @import("canopy_layer_distribution.zig").State;
const CanopyState = @import("canopy_photosynthesis.zig").State;
const InterceptionState = @import("canopy_interception.zig").State;
const source_order = @import("canopy_precipitation_retention_source_order.zig");

pub const Parameters = struct {
    surface_water_capacity_m3_per_m2_by_root_profile: [4]f64,
    low_sun_extinction_per_area_index: f64,
    minimum_solar_angle_sine_for_radiation_shares: f64,

    pub fn validate(self: Parameters) !void {
        inline for (self.surface_water_capacity_m3_per_m2_by_root_profile) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidCanopyRetentionParameters;
        if (!std.math.isFinite(self.low_sun_extinction_per_area_index) or
            self.low_sun_extinction_per_area_index < 0 or
            !std.math.isFinite(self.minimum_solar_angle_sine_for_radiation_shares) or
            self.minimum_solar_angle_sine_for_radiation_shares < 0 or
            self.minimum_solar_angle_sine_for_radiation_shares > 1)
            return error.InvalidCanopyRetentionParameters;
    }
};

pub fn compatibilityParameters() Parameters {
    return .{
        .surface_water_capacity_m3_per_m2_by_root_profile = .{
            5.0e-4,
            2.5e-4,
            2.5e-4,
            2.5e-4,
        },
        .low_sun_extinction_per_area_index = 0.65,
        .minimum_solar_angle_sine_for_radiation_shares = 0.05,
    };
}

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    living_surface_area_m2: []f64,
    standing_dead_surface_area_m2: []f64,
    living_radiation_fraction: []f64,
    standing_dead_radiation_fraction: []f64,
    living_absorbed_shortwave_mj_per_m2: []f64,
    standing_dead_absorbed_shortwave_mj_per_m2: []f64,
    living_surface_water_m3: []f64,
    standing_dead_surface_water_m3: []f64,
    living_retention_m3_per_h: []f64,
    standing_dead_retention_m3_per_h: []f64,
    previous_water_energy_mj: []f64,
    cell_potential_interception_m3_per_h: []f64,
    cell_retention_m3_per_h: []f64,
    cell_throughfall_m3_per_h: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        species_count: usize,
    ) !State {
        if (cell_count == 0 or species_count == 0)
            return error.InvalidCanopyRetentionDimensions;
        const plant_count = try std.math.mul(usize, cell_count, species_count);
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        result.species_count = species_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64 and allocated > 0) {
                allocated -= 1;
                allocator.free(@field(result, field.name));
            };
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                const count = if (comptime std.mem.startsWith(u8, field.name, "cell_"))
                    cell_count
                else
                    plant_count;
                @field(result, field.name) = try allocator.alloc(f64, count);
                @memset(@field(result, field.name), 0);
                allocated += 1;
            }
        }
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

pub fn refreshFromModel(
    state: *State,
    layers: *const LayerState,
    canopy: *const CanopyState,
    interception: *const InterceptionState,
    rainfall_m_by_cell: []const f64,
    cell_area_m2: []const f64,
    root_profile_type_by_plant: []const u8,
    solar_angle_sine_by_cell: []const f64,
    incident_ground_shortwave_mj_per_m2: []const f64,
    parameters: Parameters,
) !void {
    try parameters.validate();
    const plant_count = try validateModelDimensions(
        state,
        layers,
        canopy,
        interception,
        rainfall_m_by_cell,
        cell_area_m2,
        root_profile_type_by_plant,
        solar_angle_sine_by_cell,
        incident_ground_shortwave_mj_per_m2,
    );
    var leaf_area_by_plant_m2 = try state.allocator.alloc(f64, plant_count);
    var stalk_area_by_plant_m2 = try state.allocator.alloc(f64, plant_count);
    defer {
        state.allocator.free(leaf_area_by_plant_m2);
        state.allocator.free(stalk_area_by_plant_m2);
    }
    for (0..plant_count) |plant| {
        const cell = plant / state.species_count;
        var living_leaf_area_m2: f64 = 0;
        var living_stalk_area_m2: f64 = 0;
        var dead_area_m2: f64 = 0;
        var dead_shortwave_mj_per_m2: f64 = 0;
        for (0..layers.layer_count) |layer| {
            const plant_layer = plant * layers.layer_count + layer;
            dead_area_m2 += layers.plant_standing_dead_area_m2[plant_layer];
            dead_shortwave_mj_per_m2 +=
                interception.standing_dead_absorbed_shortwave_by_layer_mj_per_m2[
                    plant_layer
                ];
        }
        const branch_first = canopy.plant_branch_offsets[plant];
        const branch_end = canopy.plant_branch_offsets[plant + 1];
        for (branch_first..branch_end) |branch| {
            const first = branch * layers.layer_count;
            for (0..layers.layer_count) |layer|
                living_stalk_area_m2 += layers.branch_stalk_area_m2[first + layer];
            const node_first = canopy.branch_node_offsets[branch];
            const node_end = canopy.branch_node_offsets[branch + 1];
            for (node_first..node_end) |node| {
                const node_layer_first = node * layers.layer_count;
                for (0..layers.layer_count) |layer|
                    living_leaf_area_m2 += layers.node_leaf_area_m2[node_layer_first + layer];
            }
        }
        const live_shortwave = interception.absorbed_shortwave_mj_per_m2[plant];
        const living_area_m2 = living_leaf_area_m2 + living_stalk_area_m2;
        const incident = incident_ground_shortwave_mj_per_m2[cell];
        const total_absorbed = live_shortwave + dead_shortwave_mj_per_m2;
        const low_sun_fraction = 1.0 -
            std.math.exp(-parameters.low_sun_extinction_per_area_index *
                (living_area_m2 + dead_area_m2) / cell_area_m2[cell]);
        const intercepted_fraction = if (solar_angle_sine_by_cell[cell] >
            parameters.minimum_solar_angle_sine_for_radiation_shares and incident > 0)
            std.math.clamp(total_absorbed / incident, 0, 1)
        else
            std.math.clamp(low_sun_fraction, 0, 1);
        const live_share = if (total_absorbed > 0)
            live_shortwave / total_absorbed
        else if (living_area_m2 + dead_area_m2 > 0)
            living_area_m2 / (living_area_m2 + dead_area_m2)
        else
            0;
        leaf_area_by_plant_m2[plant] = living_leaf_area_m2;
        stalk_area_by_plant_m2[plant] = living_stalk_area_m2;
        state.living_surface_area_m2[plant] = living_area_m2;
        state.standing_dead_surface_area_m2[plant] = dead_area_m2;
        state.living_absorbed_shortwave_mj_per_m2[plant] = live_shortwave;
        state.standing_dead_absorbed_shortwave_mj_per_m2[plant] =
            dead_shortwave_mj_per_m2;
        state.living_radiation_fraction[plant] =
            intercepted_fraction * live_share;
        state.standing_dead_radiation_fraction[plant] =
            intercepted_fraction * (1.0 - live_share);
    }
    for (0..state.cell_count) |cell| {
        const first = cell * state.species_count;
        const last = first + state.species_count;
        const totals = try source_order.compute(
            .{
                .precipitation_irrigation_m3_h =
                    rainfall_m_by_cell[cell] * cell_area_m2[cell],
                .retention_capacity_m3_per_m2_by_vegetation_type =
                    &parameters.surface_water_capacity_m3_per_m2_by_root_profile,
                .vegetation_type_by_species =
                    root_profile_type_by_plant[first..last],
                .leaf_area_m2 = leaf_area_by_plant_m2[first..last],
                .stalk_area_m2 = stalk_area_by_plant_m2[first..last],
                .standing_dead_area_m2 = state.standing_dead_surface_area_m2[first..last],
                .live_water_content_m3 = state.living_surface_water_m3[first..last],
                .standing_dead_water_content_m3 = state.standing_dead_surface_water_m3[first..last],
                .live_radiation_fraction = state.living_radiation_fraction[first..last],
                .standing_dead_radiation_fraction = state.standing_dead_radiation_fraction[first..last],
            },
            .{
                .live_retention_flux_m3_h = state.living_retention_m3_per_h[first..last],
                .standing_dead_retention_flux_m3_h = state.standing_dead_retention_m3_per_h[first..last],
            },
        );
        state.cell_potential_interception_m3_per_h[cell] =
            totals.potential_interception_m3_h;
        state.cell_retention_m3_per_h[cell] = totals.retention_flux_m3_h;
        state.cell_throughfall_m3_per_h[cell] =
            rainfall_m_by_cell[cell] * cell_area_m2[cell] - totals.retention_flux_m3_h;
    }
}

pub fn commitRetention(state: *State, timestep_h: f64) !void {
    if (!std.math.isFinite(timestep_h) or timestep_h < 0)
        return error.InvalidCanopyRetentionTimestep;
    for (state.living_surface_water_m3, state.standing_dead_surface_water_m3, state.living_retention_m3_per_h, state.standing_dead_retention_m3_per_h) |living, dead, live_rate, dead_rate| {
        const next_living = living + live_rate * timestep_h;
        const next_dead = dead + dead_rate * timestep_h;
        if (!std.math.isFinite(next_living) or !std.math.isFinite(next_dead) or
            next_living < -1.0e-14 or next_dead < -1.0e-14)
            return error.InvalidCanopyRetentionCommit;
    }
    for (state.living_surface_water_m3, state.standing_dead_surface_water_m3, state.living_retention_m3_per_h, state.standing_dead_retention_m3_per_h) |*living, *dead, live_rate, dead_rate| {
        living.* = @max(0, living.* + live_rate * timestep_h);
        dead.* = @max(0, dead.* + dead_rate * timestep_h);
    }
}

pub fn commitSurfaceWater(
    state: *State,
    living_water_change_m3_per_h: []const f64,
    standing_dead_water_change_m3_per_h: []const f64,
    timestep_h: f64,
) !void {
    const plant_count = state.living_surface_water_m3.len;
    if (living_water_change_m3_per_h.len != plant_count or
        standing_dead_water_change_m3_per_h.len != plant_count)
        return error.CanopyRetentionDimensionMismatch;
    for (0..plant_count) |plant| {
        const next_living = state.living_surface_water_m3[plant] +
            (state.living_retention_m3_per_h[plant] +
                living_water_change_m3_per_h[plant]) * timestep_h;
        const next_dead = state.standing_dead_surface_water_m3[plant] +
            (state.standing_dead_retention_m3_per_h[plant] +
                standing_dead_water_change_m3_per_h[plant]) * timestep_h;
        if (!std.math.isFinite(next_living) or !std.math.isFinite(next_dead) or
            next_living < -1.0e-14 or next_dead < -1.0e-14)
            return error.InvalidCanopySurfaceWaterCommit;
    }
    for (0..plant_count) |plant| {
        const next_living = state.living_surface_water_m3[plant] +
            (state.living_retention_m3_per_h[plant] +
                living_water_change_m3_per_h[plant]) * timestep_h;
        const next_dead = state.standing_dead_surface_water_m3[plant] +
            (state.standing_dead_retention_m3_per_h[plant] +
                standing_dead_water_change_m3_per_h[plant]) * timestep_h;
        state.living_surface_water_m3[plant] = @max(0, next_living);
        state.standing_dead_surface_water_m3[plant] = @max(0, next_dead);
    }
}

fn retentionFlux(
    precipitation_share_m3_per_h: f64,
    capacity_m3: f64,
    stored_m3: f64,
) f64 {
    return @max(0, @min(precipitation_share_m3_per_h, capacity_m3 - stored_m3)) -
        @max(0, stored_m3 - capacity_m3);
}

fn validateModelDimensions(
    state: *const State,
    layers: *const LayerState,
    canopy: *const CanopyState,
    interception: *const InterceptionState,
    rainfall: []const f64,
    area: []const f64,
    profiles: []const u8,
    solar_sine: []const f64,
    ground_shortwave: []const f64,
) !usize {
    const plants = try std.math.mul(usize, state.cell_count, state.species_count);
    if (layers.cell_count != state.cell_count or
        layers.species_count != state.species_count or
        canopy.cell_count != state.cell_count or
        canopy.species_count != state.species_count or
        interception.cell_count != state.cell_count or
        interception.species_count != state.species_count or
        rainfall.len != state.cell_count or area.len != state.cell_count or
        solar_sine.len != state.cell_count or
        ground_shortwave.len != state.cell_count or profiles.len != plants)
        return error.CanopyRetentionDimensionMismatch;
    for (0..state.cell_count) |cell| {
        inline for (.{ rainfall[cell], area[cell], solar_sine[cell], ground_shortwave[cell] }) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidCanopyRetentionInput;
        if (area[cell] == 0 or solar_sine[cell] > 1)
            return error.InvalidCanopyRetentionInput;
    }
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state, field.name)) |value|
            if (!std.math.isFinite(value) or value < -1.0e-14)
                return error.InvalidCanopyRetentionState;
    return plants;
}

test "runtime canopy retention fills storage and publishes throughfall" {
    const flux = retentionFlux(0.2, 0.3, 0.1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), flux, 1.0e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.1),
        retentionFlux(0, 0.1, 0.2),
        1.0e-15,
    );
}

test "runtime canopy retention state supports arbitrary species count" {
    var state = try State.init(std.testing.allocator, 2, 7);
    defer state.deinit();

    try std.testing.expectEqual(@as(usize, 14), state.living_surface_water_m3.len);
    try std.testing.expectEqual(@as(usize, 2), state.cell_throughfall_m3_per_h.len);
    for (state.previous_water_energy_mj) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
}

test "retention commit rejects a late invalid value atomically" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.living_surface_water_m3[0] = 0.2;
    state.living_surface_water_m3[1] = 0.3;
    state.living_retention_m3_per_h[0] = 0.1;
    state.living_retention_m3_per_h[1] = -1;

    try std.testing.expectError(
        error.InvalidCanopyRetentionCommit,
        commitRetention(&state, 1),
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.2, 0.3 },
        state.living_surface_water_m3,
    );
}

test "surface-water commit rejects non-finite late input atomically" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.living_surface_water_m3[0] = 0.2;
    state.living_surface_water_m3[1] = 0.3;

    try std.testing.expectError(
        error.InvalidCanopySurfaceWaterCommit,
        commitSurfaceWater(
            &state,
            &.{ 0.1, std.math.nan(f64) },
            &.{ 0, 0 },
            1,
        ),
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.2, 0.3 },
        state.living_surface_water_m3,
    );
}

test "production retention flux matches HOUR1 source-order kernel" {
    const precipitation_m3_per_h = 0.8;
    const capacity_m3_per_m2 = 0.1;
    const live_area_m2 = 3.0;
    const dead_area_m2 = 1.0;
    const live_water_m3 = 0.1;
    const dead_water_m3 = 0.2;
    const live_fraction = 0.2;
    const dead_fraction = 0.1;
    var source_live: [1]f64 = undefined;
    var source_dead: [1]f64 = undefined;
    const totals = try source_order.compute(.{
        .precipitation_irrigation_m3_h = precipitation_m3_per_h,
        .retention_capacity_m3_per_m2_by_vegetation_type = &.{capacity_m3_per_m2},
        .vegetation_type_by_species = &.{0},
        .leaf_area_m2 = &.{2},
        .stalk_area_m2 = &.{1},
        .standing_dead_area_m2 = &.{dead_area_m2},
        .live_water_content_m3 = &.{live_water_m3},
        .standing_dead_water_content_m3 = &.{dead_water_m3},
        .live_radiation_fraction = &.{live_fraction},
        .standing_dead_radiation_fraction = &.{dead_fraction},
    }, .{
        .live_retention_flux_m3_h = &source_live,
        .standing_dead_retention_flux_m3_h = &source_dead,
    });

    const production_live = retentionFlux(
        precipitation_m3_per_h * live_fraction,
        capacity_m3_per_m2 * live_area_m2,
        live_water_m3,
    );
    const production_dead = retentionFlux(
        precipitation_m3_per_h * dead_fraction,
        capacity_m3_per_m2 * dead_area_m2,
        dead_water_m3,
    );
    // Both source-order and production paths now use separate leaf and stalk
    // areas before aggregation.
    try std.testing.expectApproxEqAbs(source_live[0], production_live, 1.0e-15);
    try std.testing.expectApproxEqAbs(source_dead[0], production_dead, 1.0e-15);
    try std.testing.expectApproxEqAbs(
        totals.retention_flux_m3_h,
        production_live + production_dead,
        1.0e-15,
    );
}
