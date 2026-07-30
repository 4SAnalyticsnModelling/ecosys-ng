const std = @import("std");
const inventory = @import("surface_pond_inventory_transfer.zig");
const chemistry = @import("surface_pond_chemistry_transfer.zig");
const surface_chemistry_module = @import("surface_litter_chemistry.zig");
const soil_chemistry_module = @import("solute_chemistry_state.zig");
const water_heat = @import("surface_pond_water_heat_transfer.zig");
const geometry = @import("soil_layer_geometry.zig");
const transition_module = @import("surface_pond_transition_step.zig");
const soil_properties_module = @import("soil_solver_properties.zig");
const face_geometry_module = @import("soil_face_geometry.zig");
const transport_module = @import("transport_hydrology.zig");

/// Heap-owned scratch storage for one atomic, whole-grid surface-pond
/// transaction. A complete preflight precedes every mutation.
pub const Workspace = struct {
    allocator: std.mem.Allocator,
    pond_boundary_change_m: []f64,
    zero_boundary_change_m: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, layer_capacity: usize) !Workspace {
        const count = try std.math.mul(usize, cell_count, try std.math.add(usize, layer_capacity, 1));
        const pond = try allocator.alloc(f64, count);
        errdefer allocator.free(pond);
        const zero = try allocator.alloc(f64, count);
        @memset(pond, 0);
        @memset(zero, 0);
        return .{ .allocator = allocator, .pond_boundary_change_m = pond, .zero_boundary_change_m = zero };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.zero_boundary_change_m);
        self.allocator.free(self.pond_boundary_change_m);
        self.* = undefined;
    }
};

pub const Owners = struct {
    inventories: inventory.Owners,
    surface_chemistry: *surface_chemistry_module.State,
    soil_chemistry: *soil_chemistry_module.State,
    water_heat: water_heat.Owners,
    soil_geometry: *geometry.State,
    soil_properties: *soil_properties_module.State,
    soil_faces: *const transport_module.SoilFaces,
    soil_face_geometry: *face_geometry_module.State,
};

pub const Inputs = struct {
    transitions: *const transition_module.State,
    dynamic_salts: bool,
    water_heat_parameters: water_heat.Parameters,
    minimum_heat_capacity_mj_per_k: []const f64,
    minimum_soil_layer_thickness_m: f64,
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
};

pub fn apply(workspace: *Workspace, owners: Owners, inputs: Inputs) !void {
    const grid = owners.water_heat.grid;
    const boundary_count = grid.cell_count * (grid.soil_layer_capacity + 1);
    if (inputs.transitions.cell_count != grid.cell_count or
        owners.soil_geometry.cell_count != grid.cell_count or
        owners.soil_geometry.layer_capacity != grid.soil_layer_capacity or
        workspace.pond_boundary_change_m.len != boundary_count or
        workspace.zero_boundary_change_m.len != boundary_count or
        owners.soil_properties.layer_count != grid.layer_count or
        owners.soil_properties.layer_thickness_m.len != grid.layer_count or
        owners.soil_properties.layer_midpoint_depth_m.len != grid.layer_count or
        owners.soil_properties.layer_bottom_depth_m.len != grid.layer_count or
        inputs.minimum_heat_capacity_mj_per_k.len != grid.cell_count)
        return error.SurfacePondDomainDimensionMismatch;

    @memset(workspace.pond_boundary_change_m, 0);
    @memset(workspace.zero_boundary_change_m, 0);

    // Domain-wide validation guarantees that a late invalid cell cannot leave
    // earlier cells partially transferred.
    for (0..grid.cell_count) |cell| {
        if (!inputs.transitions.active[cell]) continue;
        const layer = inputs.transitions.destination_soil_layer[cell];
        const fraction = inputs.transitions.transfer_fraction[cell];
        const destination = cell * grid.soil_layer_capacity + layer;
        const carriers = try carrierVolumes(owners, cell, destination, fraction);
        const organic_carbon_g_c = try owners.inventories.surface_organic.totalCarbon_g_c(cell);
        var heat_parameters = inputs.water_heat_parameters;
        heat_parameters.minimum_heat_capacity_mj_per_k = inputs.minimum_heat_capacity_mj_per_k[cell];
        try inventory.validateSurfaceFractionToSoil(owners.inventories, .{ .cell = cell, .destination_soil_layer = layer, .fraction = fraction });
        try chemistry.validateSurfaceFractionToSoil(owners.surface_chemistry, owners.soil_chemistry, cell, destination, carriers, inputs.dynamic_salts, fraction);
        try water_heat.validateSurfaceFractionToSoil(owners.water_heat, .{ .cell = cell, .destination_soil_layer = layer, .fraction = fraction, .surface_organic_carbon_g_c = organic_carbon_g_c, .parameters = heat_parameters });
        const boundary = cell * (grid.soil_layer_capacity + 1) + layer;
        const change = inputs.transitions.boundary_change_m[cell];
        if (!std.math.isFinite(change)) return error.InvalidSurfacePondBoundaryChange;
        workspace.pond_boundary_change_m[boundary] += change;
    }

    const changes = geometry.DisturbanceChanges{
        .pond_m = workspace.pond_boundary_change_m,
        .freeze_thaw_m = workspace.zero_boundary_change_m,
        .erosion_m = workspace.zero_boundary_change_m,
        .organic_carbon_m = workspace.zero_boundary_change_m,
    };
    try geometry.validateDisturbances(owners.soil_geometry, changes, inputs.minimum_soil_layer_thickness_m);
    try owners.soil_face_geometry.validateMapped(grid, owners.soil_faces, owners.soil_geometry.layer_thickness_m, inputs.horizontal_cell_width_m, inputs.vertical_cell_width_m);

    for (0..grid.cell_count) |cell| {
        if (!inputs.transitions.active[cell]) continue;
        const layer = inputs.transitions.destination_soil_layer[cell];
        const fraction = inputs.transitions.transfer_fraction[cell];
        const destination = cell * grid.soil_layer_capacity + layer;
        const carriers = carrierVolumes(owners, cell, destination, fraction) catch unreachable;
        const organic_carbon_g_c = owners.inventories.surface_organic.totalCarbon_g_c(cell) catch unreachable;
        var heat_parameters = inputs.water_heat_parameters;
        heat_parameters.minimum_heat_capacity_mj_per_k = inputs.minimum_heat_capacity_mj_per_k[cell];
        const moved_surface_volume_m3 = fraction * owners.water_heat.surface_geometry.expanded_total_volume_m3[cell];
        inventory.transferSurfaceFractionToSoil(owners.inventories, .{ .cell = cell, .destination_soil_layer = layer, .fraction = fraction }) catch unreachable;
        chemistry.transferSurfaceFractionToSoil(owners.surface_chemistry, owners.soil_chemistry, cell, destination, carriers, inputs.dynamic_salts, fraction) catch unreachable;
        water_heat.transferSurfaceFractionToSoil(owners.water_heat, .{ .cell = cell, .destination_soil_layer = layer, .fraction = fraction, .surface_organic_carbon_g_c = organic_carbon_g_c, .parameters = heat_parameters }) catch unreachable;
        owners.soil_properties.matrix_bulk_volume_m3[destination] += moved_surface_volume_m3;
        owners.soil_properties.layer_volume_m3[destination] = owners.water_heat.soil_thermal.layer_volume_m3[destination];
        owners.soil_properties.bulk_density_megagrams_per_m3[destination] =
            carriers.soil_dry_mass_after_Mg / owners.soil_properties.layer_volume_m3[destination];
        owners.soil_properties.porosity_fraction[destination] = owners.water_heat.soil_thermal.porosity_fraction[destination];
    }
    geometry.applyDisturbances(owners.soil_geometry, changes, inputs.minimum_soil_layer_thickness_m) catch unreachable;
    @memcpy(owners.soil_properties.layer_thickness_m, owners.soil_geometry.layer_thickness_m);
    @memcpy(owners.soil_properties.layer_midpoint_depth_m, owners.soil_geometry.layer_midpoint_depth_from_surface_m);
    @memcpy(owners.soil_properties.layer_bottom_depth_m, owners.soil_geometry.layer_bottom_depth_from_surface_m);
    owners.soil_face_geometry.refreshMapped(grid, owners.soil_faces, owners.soil_properties.layer_thickness_m, inputs.horizontal_cell_width_m, inputs.vertical_cell_width_m) catch unreachable;
}

fn carrierVolumes(owners: Owners, cell: usize, destination: usize, fraction: f64) !chemistry.CarrierVolumes {
    if (destination >= owners.soil_properties.layer_count) return error.SurfacePondDomainDimensionMismatch;
    const surface_water = owners.water_heat.surface_liquid_water_m3[cell];
    const soil_water = owners.water_heat.grid.matrix_liquid_water_m3[destination];
    const surface_dry_mass = owners.water_heat.surface_geometry.dry_mass_Mg[cell];
    const soil_dry_mass = owners.soil_properties.bulk_density_megagrams_per_m3[destination] *
        owners.water_heat.soil_thermal.layer_volume_m3[destination];
    return .{
        .surface_water_before_m3 = surface_water,
        .soil_shared_water_before_m3 = soil_water,
        .soil_phosphate_non_band_water_before_m3 = soil_water,
        .surface_water_after_m3 = (1 - fraction) * surface_water,
        .soil_shared_water_after_m3 = soil_water + fraction * surface_water,
        .soil_phosphate_non_band_water_after_m3 = soil_water + fraction * surface_water,
        .surface_dry_mass_before_Mg = surface_dry_mass,
        .soil_dry_mass_before_Mg = soil_dry_mass,
        .surface_dry_mass_after_Mg = (1 - fraction) * surface_dry_mass,
        .soil_dry_mass_after_Mg = soil_dry_mass + fraction * surface_dry_mass,
    };
}
