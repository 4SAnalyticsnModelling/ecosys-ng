const std = @import("std");

pub const State = struct {
    root_length_per_plant_m: []f64,
    root_length_density_m_per_m3: []f64,
    root_air_volume_m3: []f64,
    root_water_volume_m3: []f64,
    root_surface_area_per_plant_m2: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    host_domain_index: usize,
    planting_layer_index: usize,
    seed_length_m: f64,
    seed_volume_m3_per_plant: f64,
    seed_surface_area_m2_per_plant: f64,
    population_count: f64,
    host_root_porosity_fraction: f64,
    planting_layer_thickness_m: f64,
    minimum_active_layer_thickness_m: f64,
};

fn validate(state: State, inputs: Inputs) !void {
    const count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.SeedRootGeometryDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.host_domain_index >= inputs.biological_domain_count or inputs.planting_layer_index >= inputs.soil_layer_count) return error.SeedRootGeometryDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != count) return error.SeedRootGeometryDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSeedRootGeometryState;
    }
    inline for (.{ inputs.seed_length_m, inputs.seed_volume_m3_per_plant, inputs.seed_surface_area_m2_per_plant, inputs.population_count, inputs.host_root_porosity_fraction, inputs.planting_layer_thickness_m, inputs.minimum_active_layer_thickness_m }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSeedRootGeometryInput;
    if (inputs.host_root_porosity_fraction > 1) return error.InvalidSeedRootGeometryInput;
}

/// Exact GROSUB 7510--7523 seed geometry augmentation. Only the host domain's
/// planting layer is changed. Length is m, volume m3, and surface area m2.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validate(state, inputs);
    try validate(workspace, inputs);
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(workspace, field.name), @field(state, field.name));
    const index = inputs.host_domain_index * inputs.soil_layer_count + inputs.planting_layer_index;
    workspace.root_length_per_plant_m[index] += inputs.seed_length_m;
    workspace.root_length_density_m_per_m3[index] = if (inputs.planting_layer_thickness_m > inputs.minimum_active_layer_thickness_m) workspace.root_length_per_plant_m[index] / inputs.planting_layer_thickness_m else 0;
    const total_volume_m3 = workspace.root_air_volume_m3[index] + workspace.root_water_volume_m3[index] + inputs.seed_volume_m3_per_plant * inputs.population_count;
    workspace.root_air_volume_m3[index] = inputs.host_root_porosity_fraction * total_volume_m3;
    workspace.root_water_volume_m3[index] = (1 - inputs.host_root_porosity_fraction) * total_volume_m3;
    workspace.root_surface_area_per_plant_m2[index] += inputs.seed_surface_area_m2_per_plant;
    try validate(workspace, inputs);
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(state, field.name), @field(workspace, field.name));
}

test "GROSUB seed geometry augments only host planting layer" {
    var values = [5][6]f64{ [_]f64{1} ** 6, [_]f64{2} ** 6, [_]f64{3} ** 6, [_]f64{4} ** 6, [_]f64{5} ** 6 };
    var work: [5][6]f64 = std.mem.zeroes([5][6]f64);
    const before = values;
    const state = State{ .root_length_per_plant_m = &values[0], .root_length_density_m_per_m3 = &values[1], .root_air_volume_m3 = &values[2], .root_water_volume_m3 = &values[3], .root_surface_area_per_plant_m2 = &values[4] };
    const workspace = State{ .root_length_per_plant_m = &work[0], .root_length_density_m_per_m3 = &work[1], .root_air_volume_m3 = &work[2], .root_water_volume_m3 = &work[3], .root_surface_area_per_plant_m2 = &work[4] };
    try apply(state, workspace, .{ .biological_domain_count = 2, .soil_layer_count = 3, .host_domain_index = 0, .planting_layer_index = 1, .seed_length_m = 0.5, .seed_volume_m3_per_plant = 0.1, .seed_surface_area_m2_per_plant = 0.25, .population_count = 2, .host_root_porosity_fraction = 0.25, .planting_layer_thickness_m = 0.5, .minimum_active_layer_thickness_m = 0.01 });
    try std.testing.expectApproxEqAbs(1.5, values[0][1], 1e-12);
    try std.testing.expectApproxEqAbs(3, values[1][1], 1e-12);
    try std.testing.expectApproxEqAbs(1.8, values[2][1], 1e-12);
    try std.testing.expectApproxEqAbs(5.4, values[3][1], 1e-12);
    try std.testing.expectApproxEqAbs(5.25, values[4][1], 1e-12);
    for (0..6) |index| if (index != 1) for (0..5) |field| try std.testing.expectEqual(before[field][index], values[field][index]);
}

test "GROSUB seed geometry rejects invalid late porosity atomically" {
    var values = [5][1]f64{ .{1}, .{2}, .{3}, .{4}, .{5} };
    const before = values;
    var work: [5][1]f64 = std.mem.zeroes([5][1]f64);
    try std.testing.expectError(error.InvalidSeedRootGeometryInput, apply(.{ .root_length_per_plant_m = &values[0], .root_length_density_m_per_m3 = &values[1], .root_air_volume_m3 = &values[2], .root_water_volume_m3 = &values[3], .root_surface_area_per_plant_m2 = &values[4] }, .{ .root_length_per_plant_m = &work[0], .root_length_density_m_per_m3 = &work[1], .root_air_volume_m3 = &work[2], .root_water_volume_m3 = &work[3], .root_surface_area_per_plant_m2 = &work[4] }, .{ .biological_domain_count = 1, .soil_layer_count = 1, .host_domain_index = 0, .planting_layer_index = 0, .seed_length_m = 1, .seed_volume_m3_per_plant = 1, .seed_surface_area_m2_per_plant = 1, .population_count = 1, .host_root_porosity_fraction = 1.1, .planting_layer_thickness_m = 1, .minimum_active_layer_thickness_m = 0 }));
    try std.testing.expectEqualDeep(before, values);
}
