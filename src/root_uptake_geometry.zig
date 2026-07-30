const std = @import("std");

pub const Inputs = struct {
    axis_count: usize,
    layer_count: usize,
    primary_root_depth_m: []const f64,
    layer_thickness_m: []const f64,
    depth_to_layer_top_m: []const f64,
    seeding_depth_m: f64,
    hypocotyl_height_m: f64,
    total_biome_root_mass_g_c_by_layer: []const f64,
    root_mass_g_c_by_axis_layer: []const f64,
    root_length_density_m_per_m3_by_axis_layer: []const f64,
    root_aqueous_volume_m3_by_axis_layer: []const f64,
    root_porosity_fraction_by_axis: []const f64,
    plant_population_count: f64,
    root_length_per_plant_m_by_axis_layer: []const f64,
    micropore_fraction_by_layer: []const f64,
    minimum_root_radius_m_by_axis: []const f64,
    fallback_root_radius_m_by_axis: []const f64,
    minimum_population_fraction_multiplier: f64,
    negligible_biome_root_mass_g_c: f64,
    negligible_root_length_density_m_per_m3: f64,
    minimum_soil_layer_thickness_m: f64,
    circle_area_coefficient: f64,
    circumference_coefficient: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    axis_count: usize,
    layer_count: usize,
    rooted_layer_fraction: []f64,
    plant_fraction_of_biome_root_mass: []f64,
    minimum_population_fraction: []f64,
    effective_root_radius_m: []f64,
    uptake_path_length_m: []f64,
    root_surface_area_per_radius_m: []f64,
    maximum_water_uptake_m3_per_step: []f64,

    pub fn init(allocator: std.mem.Allocator, axis_count: usize, layer_count: usize) !State {
        if (axis_count == 0 or layer_count == 0)
            return error.InvalidRootUptakeGeometryDimensions;
        const unit_count = std.math.mul(usize, axis_count, layer_count) catch
            return error.InvalidRootUptakeGeometryDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.axis_count = axis_count;
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
                const count = if (std.mem.eql(u8, field.name, "rooted_layer_fraction"))
                    layer_count
                else
                    unit_count;
                @field(state, field.name) = try allocator.alloc(f64, count);
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

/// UPTAKE.F 494--551. Axis-major, layer-minor traversal and all branch
/// operations retain source order.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.axis_count, state.layer_count);
    errdefer staged.deinit();

    for (0..inputs.axis_count) |axis| {
        for (0..inputs.layer_count) |layer| {
            if (axis == 0) {
                var maximum_primary_root_depth_m: f64 = 0;
                for (inputs.primary_root_depth_m) |depth|
                    maximum_primary_root_depth_m = @max(maximum_primary_root_depth_m, depth);
                if (layer == 0) {
                    staged.rooted_layer_fraction[layer] = 1;
                } else if (inputs.layer_thickness_m[layer] >
                    inputs.minimum_soil_layer_thickness_m)
                {
                    const layer_top_m = inputs.depth_to_layer_top_m[layer];
                    var rooted_depth_m = @max(
                        0,
                        maximum_primary_root_depth_m - layer_top_m,
                    );
                    rooted_depth_m = @max(
                        0,
                        @min(inputs.layer_thickness_m[layer], rooted_depth_m) -
                            @max(
                                0,
                                inputs.seeding_depth_m -
                                    layer_top_m -
                                    inputs.hypocotyl_height_m,
                            ),
                    );
                    staged.rooted_layer_fraction[layer] =
                        rooted_depth_m / inputs.layer_thickness_m[layer];
                } else {
                    staged.rooted_layer_fraction[layer] = 0;
                }
            }
            const unit = axis * inputs.layer_count + layer;
            staged.plant_fraction_of_biome_root_mass[unit] =
                if (inputs.total_biome_root_mass_g_c_by_layer[layer] >
                inputs.negligible_biome_root_mass_g_c)
                    @max(0, inputs.root_mass_g_c_by_axis_layer[unit]) /
                        inputs.total_biome_root_mass_g_c_by_layer[layer]
                else
                    1;
            staged.minimum_population_fraction[unit] =
                inputs.minimum_population_fraction_multiplier *
                staged.plant_fraction_of_biome_root_mass[unit];

            if (inputs.root_length_density_m_per_m3_by_axis_layer[unit] >
                inputs.negligible_root_length_density_m_per_m3 and
                staged.rooted_layer_fraction[layer] > 0)
            {
                const root_volume_per_length_m2 =
                    (inputs.root_aqueous_volume_m3_by_axis_layer[unit] /
                        (1 - inputs.root_porosity_fraction_by_axis[axis])) /
                    (inputs.circle_area_coefficient *
                        inputs.plant_population_count *
                        inputs.root_length_per_plant_m_by_axis_layer[unit]);
                staged.effective_root_radius_m[unit] = @max(
                    inputs.minimum_root_radius_m_by_axis[axis],
                    @sqrt(root_volume_per_length_m2),
                );
                staged.uptake_path_length_m[unit] = 1 /
                    @sqrt(
                        inputs.circle_area_coefficient *
                            (inputs.root_length_density_m_per_m3_by_axis_layer[unit] /
                                staged.rooted_layer_fraction[layer]) /
                            inputs.micropore_fraction_by_layer[layer],
                    );
                staged.root_surface_area_per_radius_m[unit] =
                    inputs.circumference_coefficient *
                    inputs.root_length_per_plant_m_by_axis_layer[unit] /
                    staged.rooted_layer_fraction[layer];
            } else {
                staged.effective_root_radius_m[unit] =
                    inputs.fallback_root_radius_m_by_axis[axis];
                staged.uptake_path_length_m[unit] =
                    inputs.layer_thickness_m[layer];
                staged.root_surface_area_per_radius_m[unit] =
                    inputs.circumference_coefficient *
                    inputs.root_length_per_plant_m_by_axis_layer[unit];
            }
            staged.maximum_water_uptake_m3_per_step[unit] = 0;
        }
    }
    try validateResults(&staged);
    state.deinit();
    state.* = staged;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.axis_count != state.axis_count or inputs.layer_count != state.layer_count)
        return error.InvalidRootUptakeGeometryDimensions;
    const units = std.math.mul(usize, inputs.axis_count, inputs.layer_count) catch
        return error.InvalidRootUptakeGeometryDimensions;
    if (inputs.primary_root_depth_m.len == 0 or
        inputs.layer_thickness_m.len != inputs.layer_count or
        inputs.depth_to_layer_top_m.len != inputs.layer_count or
        inputs.total_biome_root_mass_g_c_by_layer.len != inputs.layer_count or
        inputs.micropore_fraction_by_layer.len != inputs.layer_count or
        inputs.root_porosity_fraction_by_axis.len != inputs.axis_count or
        inputs.minimum_root_radius_m_by_axis.len != inputs.axis_count or
        inputs.fallback_root_radius_m_by_axis.len != inputs.axis_count)
        return error.InvalidRootUptakeGeometryDimensions;
    inline for (.{
        inputs.root_mass_g_c_by_axis_layer,
        inputs.root_length_density_m_per_m3_by_axis_layer,
        inputs.root_aqueous_volume_m3_by_axis_layer,
        inputs.root_length_per_plant_m_by_axis_layer,
    }) |values| if (values.len != units)
        return error.InvalidRootUptakeGeometryDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == []const f64) for (@field(inputs, field.name)) |value|
            if (!std.math.isFinite(value)) return error.InvalidRootUptakeGeometryInput;
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidRootUptakeGeometryInput;
    if (inputs.plant_population_count <= 0 or
        inputs.minimum_population_fraction_multiplier < 0 or
        inputs.negligible_biome_root_mass_g_c < 0 or
        inputs.negligible_root_length_density_m_per_m3 < 0 or
        inputs.minimum_soil_layer_thickness_m < 0 or
        inputs.circle_area_coefficient <= 0 or
        inputs.circumference_coefficient <= 0)
        return error.InvalidRootUptakeGeometryInput;
    for (inputs.root_porosity_fraction_by_axis) |porosity|
        if (porosity < 0 or porosity >= 1) return error.InvalidRootUptakeGeometryInput;
    for (inputs.micropore_fraction_by_layer) |fraction|
        if (fraction <= 0) return error.InvalidRootUptakeGeometryInput;
    for (inputs.root_length_per_plant_m_by_axis_layer) |length|
        if (length <= 0) return error.InvalidRootUptakeGeometryInput;
}

fn validateResults(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state, field.name)) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteRootUptakeGeometry;
}

test "UPTAKE root geometry preserves rooted fractions competition and branches" {
    var state = try State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    try calculate(&state, .{
        .axis_count = 2,
        .layer_count = 3,
        .primary_root_depth_m = &.{ 0.35, 0.5 },
        .layer_thickness_m = &.{ 0.1, 0.2, 0.3 },
        .depth_to_layer_top_m = &.{ 0, 0.1, 0.3 },
        .seeding_depth_m = 0.05,
        .hypocotyl_height_m = 0.01,
        .total_biome_root_mass_g_c_by_layer = &.{ 10, 0, 20 },
        .root_mass_g_c_by_axis_layer = &.{ 2, -1, 5, 3, 4, 10 },
        .root_length_density_m_per_m3_by_axis_layer = &.{ 4, 0, 3, 0, 2, 1 },
        .root_aqueous_volume_m3_by_axis_layer = &.{ 0.01, 0.01, 0.02, 0.01, 0.02, 0.01 },
        .root_porosity_fraction_by_axis = &.{ 0.2, 0.3 },
        .plant_population_count = 2,
        .root_length_per_plant_m_by_axis_layer = &.{ 2, 2, 3, 1, 2, 2 },
        .micropore_fraction_by_layer = &.{ 0.5, 0.5, 0.4 },
        .minimum_root_radius_m_by_axis = &.{ 0.001, 0.002 },
        .fallback_root_radius_m_by_axis = &.{ 0.003, 0.004 },
        .minimum_population_fraction_multiplier = 0.0001,
        .negligible_biome_root_mass_g_c = 1e-12,
        .negligible_root_length_density_m_per_m3 = 0,
        .minimum_soil_layer_thickness_m = 1e-6,
        .circle_area_coefficient = 3.1416,
        .circumference_coefficient = 6.283,
    });
    try std.testing.expectEqual(@as(f64, 1), state.rooted_layer_fraction[0]);
    try std.testing.expectEqual(@as(f64, 1), state.rooted_layer_fraction[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 3.0), state.rooted_layer_fraction[2], 1e-15);
    try std.testing.expectEqual(@as(f64, 0.2), state.plant_fraction_of_biome_root_mass[0]);
    try std.testing.expectEqual(@as(f64, 1), state.plant_fraction_of_biome_root_mass[1]);
    try std.testing.expectEqual(@as(f64, 0.003), state.effective_root_radius_m[1]);
    try std.testing.expectEqual(@as(f64, 0.2), state.uptake_path_length_m[1]);
    for (state.maximum_water_uptake_m3_per_step) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
}

test "thin unrooted layer uses exact fallback geometry" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, .{
        .axis_count = 1,
        .layer_count = 2,
        .primary_root_depth_m = &.{1},
        .layer_thickness_m = &.{ 0.1, 1e-7 },
        .depth_to_layer_top_m = &.{ 0, 0.1 },
        .seeding_depth_m = 0,
        .hypocotyl_height_m = 0,
        .total_biome_root_mass_g_c_by_layer = &.{ 1, 1 },
        .root_mass_g_c_by_axis_layer = &.{ 1, 1 },
        .root_length_density_m_per_m3_by_axis_layer = &.{ 1, 1 },
        .root_aqueous_volume_m3_by_axis_layer = &.{ 1, 1 },
        .root_porosity_fraction_by_axis = &.{0},
        .plant_population_count = 1,
        .root_length_per_plant_m_by_axis_layer = &.{ 1, 1 },
        .micropore_fraction_by_layer = &.{ 1, 1 },
        .minimum_root_radius_m_by_axis = &.{0.001},
        .fallback_root_radius_m_by_axis = &.{0.02},
        .minimum_population_fraction_multiplier = 0.0001,
        .negligible_biome_root_mass_g_c = 0,
        .negligible_root_length_density_m_per_m3 = 0,
        .minimum_soil_layer_thickness_m = 1e-6,
        .circle_area_coefficient = 3.1416,
        .circumference_coefficient = 6.283,
    });
    try std.testing.expectEqual(@as(f64, 0), state.rooted_layer_fraction[1]);
    try std.testing.expectEqual(@as(f64, 0.02), state.effective_root_radius_m[1]);
    try std.testing.expectEqual(@as(f64, 1e-7), state.uptake_path_length_m[1]);
    try std.testing.expectEqual(@as(f64, 6.283), state.root_surface_area_per_radius_m[1]);
}
