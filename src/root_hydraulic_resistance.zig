const std = @import("std");

pub const Inputs = struct {
    axis_count: usize,
    layer_count: usize,
    layer_thickness_m: []const f64,
    matrix_water_volume_m3: []const f64,
    root_length_density_m_per_m3: []const f64,
    soil_hydraulic_conductivity_m2_per_h_mpa: []const f64,
    primary_axis_count_per_plant: []const f64,
    secondary_axis_count_per_plant: []const f64,
    volumetric_water_content_m3_per_m3: []const f64,
    uptake_path_length_m: []const f64,
    effective_root_radius_m: []const f64,
    root_surface_area_per_radius_m: []const f64,
    secondary_root_radius_m: []const f64,
    root_length_per_plant_m: []const f64,
    radial_resistivity_mpa_h_per_m3_by_axis: []const f64,
    micropore_volume_m3_by_layer: []const f64,
    axial_resistivity_mpa_h_per_m2_by_axis: []const f64,
    primary_root_depth_m_by_layer: []const f64,
    primary_root_radius_m: []const f64,
    average_secondary_root_length_m: []const f64,
    plant_population_count: f64,
    hydraulic_canopy_height_m: f64,
    stalk_to_root_conducting_radius_factor: f64,
    minimum_layer_thickness_m: f64,
    negligible_matrix_water_m3: f64,
    negligible_root_density_m_per_m3: f64,
    negligible_axis_count_per_plant: f64,
    negligible_water_content_m3_per_m3: f64,
    circumference_coefficient: f64,
    reference_root_radius_m: f64,
    radius_conductance_exponent: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    rooted: []bool,
    soil_hydraulic_resistance_mpa_h_per_m: []f64,
    radial_root_resistance_mpa_h_per_m: []f64,
    primary_axial_resistance_mpa_h_per_m: []f64,
    secondary_axial_resistance_mpa_h_per_m: []f64,
    root_radial_plus_axial_resistance_mpa_h_per_m: []f64,
    total_soil_root_resistance_mpa_h_per_m: []f64,
    total_soil_root_conductance_m_per_h_mpa: f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidRootHydraulicResistanceDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.total_soil_root_conductance_m_per_h_mpa = 0;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if ((field.type == []bool or field.type == []f64) and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []bool) {
                @field(state, field.name) = try allocator.alloc(bool, unit_count);
                @memset(@field(state, field.name), false);
                allocated += 1;
            } else if (field.type == []f64) {
                @field(state, field.name) = try allocator.alloc(f64, unit_count);
                @memset(@field(state, field.name), 0);
                allocated += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []bool or field.type == []f64)
                self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// UPTAKE.F 665--775. Performs the source admission/resistance pass followed
/// by the separate total-resistance/conductance pass.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    errdefer staged.deinit();
    for (0..inputs.axis_count) |axis| {
        for (0..inputs.layer_count) |layer| {
            const unit = axis * inputs.layer_count + layer;
            if (inputs.layer_thickness_m[layer] > inputs.minimum_layer_thickness_m and
                inputs.matrix_water_volume_m3[layer] > inputs.negligible_matrix_water_m3 and
                inputs.root_length_density_m_per_m3[unit] > inputs.negligible_root_density_m_per_m3 and
                inputs.soil_hydraulic_conductivity_m2_per_h_mpa[layer] > 0 and
                inputs.primary_axis_count_per_plant[unit] > inputs.negligible_axis_count_per_plant and
                inputs.secondary_axis_count_per_plant[unit] > inputs.negligible_axis_count_per_plant and
                inputs.volumetric_water_content_m3_per_m3[layer] > inputs.negligible_water_content_m3_per_m3)
            {
                staged.rooted[unit] = true;
                const soil_density_resistance =
                    @log(
                        (inputs.uptake_path_length_m[unit] +
                            inputs.effective_root_radius_m[unit]) /
                            inputs.effective_root_radius_m[unit],
                    ) /
                    inputs.root_surface_area_per_radius_m[unit];
                staged.soil_hydraulic_resistance_mpa_h_per_m[unit] =
                    soil_density_resistance /
                    inputs.soil_hydraulic_conductivity_m2_per_h_mpa[layer];
                const secondary_area = inputs.circumference_coefficient *
                    inputs.secondary_root_radius_m[unit] *
                    inputs.root_length_per_plant_m[unit];
                staged.radial_root_resistance_mpa_h_per_m[unit] =
                    inputs.radial_resistivity_mpa_h_per_m3_by_axis[axis] /
                    secondary_area *
                    inputs.micropore_volume_m3_by_layer[layer] /
                    inputs.matrix_water_volume_m3[layer];
                const primary_radius_factor = std.math.pow(
                    f64,
                    inputs.primary_root_radius_m[unit] /
                        inputs.reference_root_radius_m,
                    inputs.radius_conductance_exponent,
                );
                staged.primary_axial_resistance_mpa_h_per_m[unit] =
                    inputs.axial_resistivity_mpa_h_per_m2_by_axis[axis] *
                    inputs.primary_root_depth_m_by_layer[layer] /
                    (primary_radius_factor *
                        inputs.primary_axis_count_per_plant[unit] /
                        inputs.plant_population_count) +
                    inputs.axial_resistivity_mpa_h_per_m2_by_axis[0] *
                        inputs.hydraulic_canopy_height_m /
                        (inputs.stalk_to_root_conducting_radius_factor *
                            inputs.primary_axis_count_per_plant[unit] /
                            inputs.plant_population_count);
                const secondary_radius_factor = std.math.pow(
                    f64,
                    inputs.secondary_root_radius_m[unit] /
                        inputs.reference_root_radius_m,
                    inputs.radius_conductance_exponent,
                );
                staged.secondary_axial_resistance_mpa_h_per_m[unit] =
                    inputs.axial_resistivity_mpa_h_per_m2_by_axis[axis] *
                    inputs.average_secondary_root_length_m[unit] /
                    (secondary_radius_factor *
                        inputs.secondary_axis_count_per_plant[unit] /
                        inputs.plant_population_count);
            }
        }
    }
    for (0..inputs.axis_count) |axis| {
        for (0..inputs.layer_count) |layer| {
            const unit = axis * inputs.layer_count + layer;
            if (!staged.rooted[unit]) continue;
            staged.root_radial_plus_axial_resistance_mpa_h_per_m[unit] =
                staged.radial_root_resistance_mpa_h_per_m[unit] +
                staged.primary_axial_resistance_mpa_h_per_m[unit] +
                staged.secondary_axial_resistance_mpa_h_per_m[unit];
            staged.total_soil_root_resistance_mpa_h_per_m[unit] =
                staged.soil_hydraulic_resistance_mpa_h_per_m[unit] +
                staged.root_radial_plus_axial_resistance_mpa_h_per_m[unit];
            staged.total_soil_root_conductance_m_per_h_mpa +=
                1 / staged.total_soil_root_resistance_mpa_h_per_m[unit];
        }
    }
    try validateResults(&staged);
    state.deinit();
    state.* = staged;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.axis_count == 0 or inputs.layer_count == 0) return error.InvalidRootHydraulicResistanceDimensions;
    const units = std.math.mul(usize, inputs.axis_count, inputs.layer_count) catch return error.InvalidRootHydraulicResistanceDimensions;
    if (state.unit_count != units or
        inputs.layer_thickness_m.len != inputs.layer_count or
        inputs.matrix_water_volume_m3.len != inputs.layer_count or
        inputs.soil_hydraulic_conductivity_m2_per_h_mpa.len != inputs.layer_count or
        inputs.volumetric_water_content_m3_per_m3.len != inputs.layer_count or
        inputs.micropore_volume_m3_by_layer.len != inputs.layer_count or
        inputs.primary_root_depth_m_by_layer.len != inputs.layer_count or
        inputs.radial_resistivity_mpa_h_per_m3_by_axis.len != inputs.axis_count or
        inputs.axial_resistivity_mpa_h_per_m2_by_axis.len != inputs.axis_count)
        return error.InvalidRootHydraulicResistanceDimensions;
    inline for (.{
        inputs.root_length_density_m_per_m3,
        inputs.primary_axis_count_per_plant,
        inputs.secondary_axis_count_per_plant,
        inputs.uptake_path_length_m,
        inputs.effective_root_radius_m,
        inputs.root_surface_area_per_radius_m,
        inputs.secondary_root_radius_m,
        inputs.root_length_per_plant_m,
        inputs.primary_root_radius_m,
        inputs.average_secondary_root_length_m,
    }) |values| if (values.len != units) return error.InvalidRootHydraulicResistanceDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == []const f64) for (@field(inputs, field.name)) |value|
            if (!std.math.isFinite(value)) return error.InvalidRootHydraulicResistanceInput;
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidRootHydraulicResistanceInput;
    if (inputs.plant_population_count <= 0 or
        inputs.hydraulic_canopy_height_m < 0 or
        inputs.stalk_to_root_conducting_radius_factor <= 0 or
        inputs.minimum_layer_thickness_m < 0 or
        inputs.negligible_matrix_water_m3 < 0 or
        inputs.negligible_root_density_m_per_m3 < 0 or
        inputs.negligible_axis_count_per_plant < 0 or
        inputs.negligible_water_content_m3_per_m3 < 0 or
        inputs.circumference_coefficient <= 0 or
        inputs.reference_root_radius_m <= 0 or
        inputs.radius_conductance_exponent <= 0)
        return error.InvalidRootHydraulicResistanceInput;
    inline for (.{
        inputs.uptake_path_length_m,
        inputs.effective_root_radius_m,
        inputs.root_surface_area_per_radius_m,
        inputs.secondary_root_radius_m,
        inputs.root_length_per_plant_m,
        inputs.primary_root_radius_m,
        inputs.average_secondary_root_length_m,
    }) |values| for (values) |value| if (value <= 0)
        return error.InvalidRootHydraulicResistanceInput;
}

fn validateResults(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64 and (!std.math.isFinite(@field(state, field.name)) or @field(state, field.name) < 0))
            return error.NonFiniteRootHydraulicResistance;
        if (field.type == []f64) for (@field(state, field.name)) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteRootHydraulicResistance;
    }
    for (state.rooted, state.total_soil_root_resistance_mpa_h_per_m) |rooted, resistance|
        if (rooted and resistance <= 0) return error.InvalidRootHydraulicResistance;
}

test "UPTAKE two-pass hydraulic resistance preserves active and inactive units" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, .{
        .axis_count = 1,
        .layer_count = 2,
        .layer_thickness_m = &.{ 0.2, 0.2 },
        .matrix_water_volume_m3 = &.{ 0.5, 0 },
        .root_length_density_m_per_m3 = &.{ 2, 2 },
        .soil_hydraulic_conductivity_m2_per_h_mpa = &.{ 0.1, 0.1 },
        .primary_axis_count_per_plant = &.{ 2, 2 },
        .secondary_axis_count_per_plant = &.{ 4, 4 },
        .volumetric_water_content_m3_per_m3 = &.{ 0.3, 0.3 },
        .uptake_path_length_m = &.{ 0.05, 0.05 },
        .effective_root_radius_m = &.{ 0.002, 0.002 },
        .root_surface_area_per_radius_m = &.{ 10, 10 },
        .secondary_root_radius_m = &.{ 0.001, 0.001 },
        .root_length_per_plant_m = &.{ 3, 3 },
        .radial_resistivity_mpa_h_per_m3_by_axis = &.{0.2},
        .micropore_volume_m3_by_layer = &.{ 0.8, 0.8 },
        .axial_resistivity_mpa_h_per_m2_by_axis = &.{0.4},
        .primary_root_depth_m_by_layer = &.{ 0.1, 0.3 },
        .primary_root_radius_m = &.{ 0.002, 0.002 },
        .average_secondary_root_length_m = &.{ 0.05, 0.05 },
        .plant_population_count = 5,
        .hydraulic_canopy_height_m = 1.6,
        .stalk_to_root_conducting_radius_factor = 1000,
        .minimum_layer_thickness_m = 1e-6,
        .negligible_matrix_water_m3 = 1e-12,
        .negligible_root_density_m_per_m3 = 0,
        .negligible_axis_count_per_plant = 0,
        .negligible_water_content_m3_per_m3 = 0,
        .circumference_coefficient = 6.283,
        .reference_root_radius_m = 0.1e-3,
        .radius_conductance_exponent = 4,
    });
    try std.testing.expect(state.rooted[0]);
    try std.testing.expect(!state.rooted[1]);
    try std.testing.expect(state.total_soil_root_resistance_mpa_h_per_m[0] > 0);
    try std.testing.expectEqual(
        1 / state.total_soil_root_resistance_mpa_h_per_m[0],
        state.total_soil_root_conductance_m_per_h_mpa,
    );
    try std.testing.expectEqual(@as(f64, 0), state.total_soil_root_resistance_mpa_h_per_m[1]);
}

test "hydraulic resistance supports runtime multiple axes and layers" {
    var state = try State.init(std.testing.allocator, 6);
    defer state.deinit();
    try std.testing.expectEqual(@as(usize, 6), state.rooted.len);
    try std.testing.expectEqual(@as(usize, 6), state.total_soil_root_resistance_mpa_h_per_m.len);
}
