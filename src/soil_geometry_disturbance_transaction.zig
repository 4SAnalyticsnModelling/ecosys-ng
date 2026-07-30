const std = @import("std");
const Geometry = @import("soil_layer_geometry.zig");
const Assembly = @import("soil_geometry_change_assembly.zig");

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    pond_boundary_change_m: []f64,
    freeze_thaw_boundary_change_m: []f64,
    erosion_boundary_change_m: []f64,
    organic_carbon_boundary_change_m: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, layer_capacity: usize) !Workspace {
        if (cell_count == 0 or layer_capacity == 0) return error.InvalidSoilGeometryDisturbanceDimensions;
        const count = try std.math.mul(usize, cell_count, try std.math.add(usize, layer_capacity, 1));
        var result: Workspace = undefined;
        result.allocator = allocator;
        result.pond_boundary_change_m = try zeroes(allocator, count);
        errdefer allocator.free(result.pond_boundary_change_m);
        result.freeze_thaw_boundary_change_m = try zeroes(allocator, count);
        errdefer allocator.free(result.freeze_thaw_boundary_change_m);
        result.erosion_boundary_change_m = try zeroes(allocator, count);
        errdefer allocator.free(result.erosion_boundary_change_m);
        result.organic_carbon_boundary_change_m = try zeroes(allocator, count);
        return result;
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.organic_carbon_boundary_change_m);
        self.allocator.free(self.erosion_boundary_change_m);
        self.allocator.free(self.freeze_thaw_boundary_change_m);
        self.allocator.free(self.pond_boundary_change_m);
        self.* = undefined;
    }

    pub fn changes(self: *const Workspace) Geometry.DisturbanceChanges {
        return .{
            .pond_m = self.pond_boundary_change_m,
            .freeze_thaw_m = self.freeze_thaw_boundary_change_m,
            .erosion_m = self.erosion_boundary_change_m,
            .organic_carbon_m = self.organic_carbon_boundary_change_m,
        };
    }
};

pub const Inputs = struct {
    total_ice_volume_change_m3: []const f64,
    soil_matrix_fraction: []const f64,
    net_sediment_Mg_by_cell: []const f64,
    snow_deposited_sediment_Mg_by_cell: []const f64,
    horizontal_area_m2_by_cell: []const f64,
    surface_soil_mass_Mg_by_cell: []const f64,
    surface_soil_volume_m3_by_cell: []const f64,
    receiving_soil_bulk_density_Mg_per_m3_by_cell: []const f64,
    organic_carbon_change_g_c: []const f64,
    macropore_fraction: []const f64,
    reference_bulk_density_Mg_per_m3: []const f64,
    initial_layer_thickness_m: []const f64,
    current_layer_thickness_m: []const f64,
    reset_organic_accumulation_by_layer: []const bool,
    ice_to_water_specific_volume_difference: f64,
    organic_carbon_specific_volume_m3_per_g: f64,
    freeze_thaw_enabled: bool,
    erosion_enabled: bool,
    organic_carbon_enabled: bool,
    negligible_ice_volume_change_m3: f64,
    negligible_sediment_Mg: f64,
    negligible_carbon_change_g_c: f64,
};

/// Stages every REDIST boundary driver without touching live geometry. A caller
/// commits `workspace.changes()` only after all extensive and intensive pool
/// remaps have also validated.
pub fn stage(workspace: *Workspace, geometry: *const Geometry.State, inputs: Inputs) !void {
    const expected = geometry.cell_count * (geometry.layer_capacity + 1);
    if (workspace.pond_boundary_change_m.len != expected or workspace.freeze_thaw_boundary_change_m.len != expected or workspace.erosion_boundary_change_m.len != expected or workspace.organic_carbon_boundary_change_m.len != expected) return error.SoilGeometryDisturbanceWorkspaceMismatch;
    @memset(workspace.pond_boundary_change_m, 0);
    if (inputs.freeze_thaw_enabled) {
        try Assembly.assembleFreezeThawBoundaryChangeM(workspace.freeze_thaw_boundary_change_m, geometry, inputs.total_ice_volume_change_m3, inputs.soil_matrix_fraction, inputs.horizontal_area_m2_by_cell, inputs.ice_to_water_specific_volume_difference, inputs.negligible_ice_volume_change_m3);
    } else {
        @memset(workspace.freeze_thaw_boundary_change_m, 0);
    }
    try Assembly.assembleErosionBoundaryChangeM(workspace.erosion_boundary_change_m, geometry, inputs.net_sediment_Mg_by_cell, inputs.snow_deposited_sediment_Mg_by_cell, inputs.horizontal_area_m2_by_cell, inputs.surface_soil_mass_Mg_by_cell, inputs.surface_soil_volume_m3_by_cell, inputs.receiving_soil_bulk_density_Mg_per_m3_by_cell, inputs.erosion_enabled, inputs.negligible_sediment_Mg);
    try Assembly.assembleOrganicCarbonBoundaryChangeM(workspace.organic_carbon_boundary_change_m, geometry, inputs.organic_carbon_change_g_c, inputs.macropore_fraction, inputs.reference_bulk_density_Mg_per_m3, inputs.initial_layer_thickness_m, inputs.current_layer_thickness_m, inputs.reset_organic_accumulation_by_layer, inputs.horizontal_area_m2_by_cell, inputs.organic_carbon_specific_volume_m3_per_g, inputs.organic_carbon_enabled, inputs.negligible_carbon_change_g_c);
}

fn zeroes(allocator: std.mem.Allocator, count: usize) ![]f64 {
    const values = try allocator.alloc(f64, count);
    @memset(values, 0);
    return values;
}

test "combined REDIST geometry drivers stage before one atomic geometry commit" {
    var geometry = try Geometry.State.init(std.testing.allocator, 1, 2);
    defer geometry.deinit();
    try Geometry.initializeCell(&geometry, 0, 0, &.{ 0.2, 0.3 }, 0, 1e-9);
    var workspace = try Workspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try stage(&workspace, &geometry, .{
        .total_ice_volume_change_m3 = &.{ 0.1, 0 },
        .soil_matrix_fraction = &.{ 1, 1 },
        .net_sediment_Mg_by_cell = &.{1},
        .snow_deposited_sediment_Mg_by_cell = &.{0},
        .horizontal_area_m2_by_cell = &.{10},
        .surface_soil_mass_Mg_by_cell = &.{20},
        .surface_soil_volume_m3_by_cell = &.{10},
        .receiving_soil_bulk_density_Mg_per_m3_by_cell = &.{2},
        .organic_carbon_change_g_c = &.{ -100, 0 },
        .macropore_fraction = &.{ 0, 0 },
        .reference_bulk_density_Mg_per_m3 = &.{ 1, 1 },
        .initial_layer_thickness_m = &.{ 0.2, 0.3 },
        .current_layer_thickness_m = &.{ 0.2, 0.3 },
        .reset_organic_accumulation_by_layer = &.{ false, true },
        .ice_to_water_specific_volume_difference = 0.083,
        .organic_carbon_specific_volume_m3_per_g = 1.82e-6,
        .freeze_thaw_enabled = true,
        .erosion_enabled = true,
        .organic_carbon_enabled = true,
        .negligible_ice_volume_change_m3 = 0,
        .negligible_sediment_Mg = 0,
        .negligible_carbon_change_g_c = 0,
    });
    const old_bottom = geometry.boundary_depth_m[2];
    try Geometry.applyDisturbances(&geometry, workspace.changes(), 1e-9);
    // Erosion moves every boundary equally; only freeze/SOC alter thickness.
    try std.testing.expect(geometry.boundary_depth_m[0] != 0);
    try std.testing.expectApproxEqAbs(old_bottom + 0.05, geometry.boundary_depth_m[2], 1e-14);
}
