const std = @import("std");

pub const Inputs = struct {
    landscape_total_water_potential_mpa: []const f64,
    soil_bulk_density_megagrams_per_m3: []const f64,
    current_matrix_water_m3: []const f64,
    hygroscopic_water_content_m3_per_m3: []const f64,
    matrix_bulk_volume_m3: []const f64,
    total_pore_volume_m3: []const f64,
    total_water_m3: []const f64,
    total_ice_m3: []const f64,
    root_mass_g_c_by_layer_unit: []const f64,
    root_unit_count: usize,
    surface_elevation_m: f64,
    gravitational_water_potential_mpa_per_m: f64,
    negligible_bulk_density_megagrams_per_m3: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    root_referenced_total_water_potential_mpa: []f64,
    root_available_water_m3: []f64,
    air_volume_m3: []f64,
    total_biome_root_mass_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.InvalidRootZoneWaterDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.layer_count = layer_count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(state, field.name) = try allocator.alloc(f64, layer_count);
                @memset(@field(state, field.name), 0);
                allocated += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// UPTAKE.F 410--424. Initializes root-zone water, air, potential, and
/// biome-root-mass work arrays in source statement and summation order.
pub fn initialize(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.layer_count);
    errdefer staged.deinit();
    for (0..state.layer_count) |layer| {
        staged.root_referenced_total_water_potential_mpa[layer] =
            inputs.landscape_total_water_potential_mpa[layer] -
            inputs.gravitational_water_potential_mpa_per_m *
                inputs.surface_elevation_m;
        if (inputs.soil_bulk_density_megagrams_per_m3[layer] >
            inputs.negligible_bulk_density_megagrams_per_m3)
        {
            staged.root_available_water_m3[layer] =
                inputs.current_matrix_water_m3[layer] -
                inputs.hygroscopic_water_content_m3_per_m3[layer] *
                    inputs.matrix_bulk_volume_m3[layer];
            staged.air_volume_m3[layer] = @max(
                0,
                inputs.total_pore_volume_m3[layer] -
                    inputs.total_water_m3[layer] -
                    inputs.total_ice_m3[layer],
            );
        } else {
            staged.root_available_water_m3[layer] =
                inputs.current_matrix_water_m3[layer];
            staged.air_volume_m3[layer] = 0;
        }
        var root_mass_g_c: f64 = 0;
        const first = layer * inputs.root_unit_count;
        for (inputs.root_mass_g_c_by_layer_unit[first .. first + inputs.root_unit_count]) |mass|
            root_mass_g_c += @max(0, mass);
        staged.total_biome_root_mass_g_c[layer] = root_mass_g_c;
    }
    try validateResults(&staged);
    state.deinit();
    state.* = staged;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const count = state.layer_count;
    if (count == 0 or inputs.root_unit_count == 0)
        return error.InvalidRootZoneWaterDimensions;
    inline for (.{
        inputs.landscape_total_water_potential_mpa,
        inputs.soil_bulk_density_megagrams_per_m3,
        inputs.current_matrix_water_m3,
        inputs.hygroscopic_water_content_m3_per_m3,
        inputs.matrix_bulk_volume_m3,
        inputs.total_pore_volume_m3,
        inputs.total_water_m3,
        inputs.total_ice_m3,
    }) |values| if (values.len != count)
        return error.InvalidRootZoneWaterDimensions;
    const root_count = std.math.mul(usize, count, inputs.root_unit_count) catch
        return error.InvalidRootZoneWaterDimensions;
    if (inputs.root_mass_g_c_by_layer_unit.len != root_count)
        return error.InvalidRootZoneWaterDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == []const f64) for (@field(inputs, field.name)) |value|
            if (!std.math.isFinite(value)) return error.InvalidRootZoneWaterInput;
    }
    inline for (.{
        inputs.surface_elevation_m,
        inputs.gravitational_water_potential_mpa_per_m,
        inputs.negligible_bulk_density_megagrams_per_m3,
    }) |value| if (!std.math.isFinite(value))
        return error.InvalidRootZoneWaterInput;
    if (inputs.gravitational_water_potential_mpa_per_m < 0 or
        inputs.negligible_bulk_density_megagrams_per_m3 < 0)
        return error.InvalidRootZoneWaterInput;
}

fn validateResults(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state, field.name)) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteRootZoneWaterInitialization;
}

test "UPTAKE initializes mineral and non-mineral layers in source order" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try initialize(&state, .{
        .landscape_total_water_potential_mpa = &.{ 0.2, -0.4 },
        .soil_bulk_density_megagrams_per_m3 = &.{ 1.2, 0 },
        .current_matrix_water_m3 = &.{ 0.6, 0.3 },
        .hygroscopic_water_content_m3_per_m3 = &.{ 0.1, 0.2 },
        .matrix_bulk_volume_m3 = &.{ 2, 1 },
        .total_pore_volume_m3 = &.{ 1.5, 1 },
        .total_water_m3 = &.{ 0.7, 0.2 },
        .total_ice_m3 = &.{ 0.1, 0.1 },
        .root_mass_g_c_by_layer_unit = &.{ 1, -2, 3, 4, 5, -6 },
        .root_unit_count = 3,
        .surface_elevation_m = 10,
        .gravitational_water_potential_mpa_per_m = 0.0098,
        .negligible_bulk_density_megagrams_per_m3 = 1e-12,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.102), state.root_referenced_total_water_potential_mpa[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.498), state.root_referenced_total_water_potential_mpa[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.root_available_water_m3[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 0.3), state.root_available_water_m3[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), state.air_volume_m3[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), state.air_volume_m3[1]);
    try std.testing.expectEqual([2]f64{ 4, 9 }, state.total_biome_root_mass_g_c[0..2].*);
}

test "root-zone initialization preserves negative source available water" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try initialize(&state, .{
        .landscape_total_water_potential_mpa = &.{0},
        .soil_bulk_density_megagrams_per_m3 = &.{1},
        .current_matrix_water_m3 = &.{0.1},
        .hygroscopic_water_content_m3_per_m3 = &.{0.2},
        .matrix_bulk_volume_m3 = &.{1},
        .total_pore_volume_m3 = &.{1},
        .total_water_m3 = &.{2},
        .total_ice_m3 = &.{0},
        .root_mass_g_c_by_layer_unit = &.{0},
        .root_unit_count = 1,
        .surface_elevation_m = 0,
        .gravitational_water_potential_mpa_per_m = 0.0098,
        .negligible_bulk_density_megagrams_per_m3 = 0,
    });
    try std.testing.expectEqual(@as(f64, -0.1), state.root_available_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 0), state.air_volume_m3[0]);
}

test "late root-mass overflow leaves initialization state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.root_available_water_m3[0] = 7;
    try std.testing.expectError(
        error.NonFiniteRootZoneWaterInitialization,
        initialize(&state, .{
            .landscape_total_water_potential_mpa = &.{0},
            .soil_bulk_density_megagrams_per_m3 = &.{1},
            .current_matrix_water_m3 = &.{1},
            .hygroscopic_water_content_m3_per_m3 = &.{0},
            .matrix_bulk_volume_m3 = &.{1},
            .total_pore_volume_m3 = &.{1},
            .total_water_m3 = &.{0},
            .total_ice_m3 = &.{0},
            .root_mass_g_c_by_layer_unit = &.{ std.math.floatMax(f64), std.math.floatMax(f64) },
            .root_unit_count = 2,
            .surface_elevation_m = 0,
            .gravitational_water_potential_mpa_per_m = 0.0098,
            .negligible_bulk_density_megagrams_per_m3 = 0,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.root_available_water_m3[0]);
}
