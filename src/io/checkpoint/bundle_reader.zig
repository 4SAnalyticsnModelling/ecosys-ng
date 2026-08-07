const std = @import("std");
const manifest_module = @import("manifest.zig");
const bundle_io = @import("bundle_io.zig");
const checkpoint = @import("checkpoint.zig");
const metadata_checkpoint = @import("plant_checkpoint_metadata.zig");
const development_checkpoint = @import("plant_development_checkpoint.zig");
const root_checkpoint = @import("plant_root_checkpoint.zig");
const canopy_checkpoint = @import("canopy_state_checkpoint.zig");
const biogeochemistry_checkpoint = @import("soil_biogeochemistry_checkpoint.zig");
const organic_checkpoint = @import("soil_organic_checkpoint.zig");
const transport_checkpoint = @import("transport_state_checkpoint.zig");
const geometry_checkpoint = @import("soil_geometry_checkpoint.zig");
const mass_balance_checkpoint = @import("landscape_mass_balance_checkpoint.zig");
const config_module = @import("../../core/config.zig");
const GridState = @import("../../state/grid.zig").GridState;
const PlantState = @import("../../state/grid.zig").PlantState;
const RootState = @import("../../plant/root/plant_root_system.zig").State;

pub const Limits = struct {
    manifest: manifest_module.Limits,
    plant_metadata: metadata_checkpoint.Limits,
    plant_development: development_checkpoint.Limits,
    plant_roots: root_checkpoint.Limits,
    plant_canopy: canopy_checkpoint.Limits,
    soil_biogeochemistry: biogeochemistry_checkpoint.Limits,
    soil_organic: organic_checkpoint.Limits,
    transport: transport_checkpoint.Limits,
    soil_geometry: geometry_checkpoint.Limits,
};

pub const Settings = struct {
    manifest_file_name: []const u8 = "restart.manifest",
    manifest_buffer_bytes: usize,
    section_read_buffer_bytes: usize,
    section_verify_buffer_bytes: usize,
};

pub const OwnedBundle = struct {
    manifest: manifest_module.Manifest,
    grid: GridState,
    plants: PlantState,
    plant_metadata: metadata_checkpoint.Metadata,
    plant_development: development_checkpoint.Owned,
    plant_roots: RootState,
    plant_canopy: canopy_checkpoint.Owned,
    soil_biogeochemistry: biogeochemistry_checkpoint.Owned,
    soil_organic_matter: organic_checkpoint.Owned,
    transport: transport_checkpoint.Owned,
    soil_geometry_and_hydrology: geometry_checkpoint.Owned,
    landscape_mass_balance: mass_balance_checkpoint.State,

    pub fn deinit(self: *OwnedBundle) void {
        self.soil_geometry_and_hydrology.deinit();
        self.transport.deinit();
        self.soil_organic_matter.deinit();
        self.soil_biogeochemistry.deinit();
        self.plant_canopy.deinit();
        self.plant_roots.deinit();
        self.plant_development.deinit();
        self.plant_metadata.deinit();
        self.plants.deinit();
        self.grid.deinit();
        self.manifest.deinit();
        self.* = undefined;
    }
};

pub const LiveTargets = struct {
    grid: *GridState,
    plants: *PlantState,
    plant_development: development_checkpoint.View,
    plant_roots: *RootState,
    plant_canopy: canopy_checkpoint.View,
    soil_biogeochemistry: biogeochemistry_checkpoint.View,
    soil_organic_matter: organic_checkpoint.View,
    transport: transport_checkpoint.View,
    soil_geometry_and_hydrology: geometry_checkpoint.View,
    landscape_mass_balance: *mass_balance_checkpoint.State,
};

/// The fallible validation happens before the first exchange. After it
/// succeeds every operation is a no-fail pointer/slice-owner swap. `bundle`
/// then owns the replaced states, allowing its normal `deinit` to release them.
pub fn swapIntoLive(bundle: *OwnedBundle, targets: LiveTargets) !void {
    try validateSwapTargets(bundle.*, targets);
    const runtime_targets = targets.soil_geometry_and_hydrology.runtime orelse
        return error.MissingSoilRuntimeCheckpointTarget;
    try bundle.soil_geometry_and_hydrology.runtime.restoreInto(
        @constCast(runtime_targets.soil_properties),
        @constCast(runtime_targets.soil_thermal),
    );
    const surface_boundary_targets =
        targets.soil_geometry_and_hydrology.surface_boundary orelse
        return error.MissingSurfaceBoundaryCheckpointTarget;
    try bundle.soil_geometry_and_hydrology.surface_boundary.restoreInto(
        @constCast(surface_boundary_targets.ground_air),
        @constCast(surface_boundary_targets.surface_aerodynamics),
    );
    std.mem.swap(GridState, &bundle.grid, targets.grid);
    std.mem.swap(PlantState, &bundle.plants, targets.plants);
    std.mem.swap(@TypeOf(bundle.plant_development.phenology), &bundle.plant_development.phenology, @constCast(targets.plant_development.phenology));
    std.mem.swap(@TypeOf(bundle.plant_development.growth), &bundle.plant_development.growth, @constCast(targets.plant_development.growth));
    std.mem.swap(@TypeOf(bundle.plant_development.dormancy), &bundle.plant_development.dormancy, @constCast(targets.plant_development.dormancy));
    std.mem.swap(@TypeOf(bundle.plant_development.branch_development), &bundle.plant_development.branch_development, @constCast(targets.plant_development.branch_development));
    std.mem.swap(RootState, &bundle.plant_roots, targets.plant_roots);
    std.mem.swap(@TypeOf(bundle.plant_canopy.canopy), &bundle.plant_canopy.canopy, @constCast(targets.plant_canopy.canopy));
    std.mem.swap(@TypeOf(bundle.plant_canopy.retention), &bundle.plant_canopy.retention, @constCast(targets.plant_canopy.retention));
    std.mem.swap(@TypeOf(bundle.plant_canopy.layer_distribution), &bundle.plant_canopy.layer_distribution, @constCast(targets.plant_canopy.layer_distribution));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.microbial), &bundle.soil_biogeochemistry.microbial, @constCast(targets.soil_biogeochemistry.microbial));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.chemistry), &bundle.soil_biogeochemistry.chemistry, @constCast(targets.soil_biogeochemistry.chemistry));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.available_nutrients), &bundle.soil_biogeochemistry.available_nutrients, @constCast(targets.soil_biogeochemistry.available_nutrients));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.fertilizer), &bundle.soil_biogeochemistry.fertilizer, @constCast(targets.soil_biogeochemistry.fertilizer));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.mineral_fertilizer), &bundle.soil_biogeochemistry.mineral_fertilizer, @constCast(targets.soil_biogeochemistry.mineral_fertilizer));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.fertilizer_band), &bundle.soil_biogeochemistry.fertilizer_band, @constCast(targets.soil_biogeochemistry.fertilizer_band));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.reactive_nitrogen), &bundle.soil_biogeochemistry.reactive_nitrogen, @constCast(targets.soil_biogeochemistry.reactive_nitrogen));
    std.mem.swap(@TypeOf(bundle.soil_biogeochemistry.microbial_phosphorus), &bundle.soil_biogeochemistry.microbial_phosphorus, @constCast(targets.soil_biogeochemistry.microbial_phosphorus));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.profile), &bundle.soil_organic_matter.profile, @constCast(targets.soil_organic_matter.profile));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.surface), &bundle.soil_organic_matter.surface, @constCast(targets.soil_organic_matter.surface));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.litter_chemistry), &bundle.soil_organic_matter.litter_chemistry, @constCast(targets.soil_organic_matter.litter_chemistry));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.litter_fertilizer), &bundle.soil_organic_matter.litter_fertilizer, @constCast(targets.soil_organic_matter.litter_fertilizer));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.surface_respiration), &bundle.soil_organic_matter.surface_respiration, @constCast(targets.soil_organic_matter.surface_respiration));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.surface_denitrification), &bundle.soil_organic_matter.surface_denitrification, @constCast(targets.soil_organic_matter.surface_denitrification));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.surface_fire_exchange), &bundle.soil_organic_matter.surface_fire_exchange, @constCast(targets.soil_organic_matter.surface_fire_exchange));
    std.mem.swap(@TypeOf(bundle.soil_organic_matter.litter_salt_ingress), &bundle.soil_organic_matter.litter_salt_ingress, @constCast(targets.soil_organic_matter.litter_salt_ingress));
    std.mem.swap(@TypeOf(bundle.transport.micropore), &bundle.transport.micropore, @constCast(targets.transport.micropore));
    std.mem.swap(@TypeOf(bundle.transport.macropore), &bundle.transport.macropore, @constCast(targets.transport.macropore));
    std.mem.swap(@TypeOf(bundle.transport.mineral_nitrogen), &bundle.transport.mineral_nitrogen, @constCast(targets.transport.mineral_nitrogen));
    std.mem.swap(@TypeOf(bundle.transport.organic), &bundle.transport.organic, @constCast(targets.transport.organic));
    std.mem.swap(@TypeOf(bundle.transport.gas), &bundle.transport.gas, @constCast(targets.transport.gas));
    std.mem.swap(@TypeOf(bundle.transport.litter_gas), &bundle.transport.litter_gas, @constCast(targets.transport.litter_gas));
    std.mem.swap(@TypeOf(bundle.transport.snow), &bundle.transport.snow, @constCast(targets.transport.snow));
    std.mem.swap(@TypeOf(bundle.transport.surface), &bundle.transport.surface, @constCast(targets.transport.surface));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.geometry), &bundle.soil_geometry_and_hydrology.geometry, @constCast(targets.soil_geometry_and_hydrology.geometry));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.hydrology), &bundle.soil_geometry_and_hydrology.hydrology, @constCast(targets.soil_geometry_and_hydrology.hydrology));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.surface), &bundle.soil_geometry_and_hydrology.surface, @constCast(targets.soil_geometry_and_hydrology.surface));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.surface_litter_geometry), &bundle.soil_geometry_and_hydrology.surface_litter_geometry, @constCast(targets.soil_geometry_and_hydrology.surface_litter_geometry));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.erosion), &bundle.soil_geometry_and_hydrology.erosion, @constCast(targets.soil_geometry_and_hydrology.erosion));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.climate), &bundle.soil_geometry_and_hydrology.climate, @constCast(targets.soil_geometry_and_hydrology.climate));
    std.mem.swap(@TypeOf(bundle.soil_geometry_and_hydrology.eroded_minerals), &bundle.soil_geometry_and_hydrology.eroded_minerals, @constCast(targets.soil_geometry_and_hydrology.eroded_minerals));
    @memcpy(@constCast(targets.soil_geometry_and_hydrology.delayed_live_canopy_combustion_heat_megajoules), bundle.soil_geometry_and_hydrology.delayed_live_canopy_combustion_heat_megajoules);
    @memcpy(@constCast(targets.soil_geometry_and_hydrology.delayed_standing_dead_combustion_heat_megajoules), bundle.soil_geometry_and_hydrology.delayed_standing_dead_combustion_heat_megajoules);
    @memcpy(@constCast(targets.soil_geometry_and_hydrology.delayed_subsurface_combustion_heat_megajoules), bundle.soil_geometry_and_hydrology.delayed_subsurface_combustion_heat_megajoules);
    @memcpy(@constCast(targets.soil_geometry_and_hydrology.delayed_surface_combustion_heat_megajoules), bundle.soil_geometry_and_hydrology.delayed_surface_combustion_heat_megajoules);
    @memcpy(@constCast(targets.soil_geometry_and_hydrology.surface_litter_ice_m3), bundle.soil_geometry_and_hydrology.surface_litter_ice_m3);
    targets.landscape_mass_balance.* = bundle.landscape_mass_balance;
}

/// Restores only into newly allocated owners. A caller may swap this bundle
/// into the executable after success; any error deinitializes every temporary
/// owner and leaves the live simulation unchanged.
pub fn read(
    allocator: std.mem.Allocator,
    io: std.Io,
    directory: std.Io.Dir,
    runtime_config: config_module.SimulationConfig,
    limits: Limits,
    settings: Settings,
) !OwnedBundle {
    if (settings.manifest_buffer_bytes == 0 or settings.section_read_buffer_bytes == 0 or settings.section_verify_buffer_bytes == 0) return error.InvalidCheckpointBundleBufferSize;
    var manifest = try readManifest(allocator, io, directory, settings.manifest_file_name, settings.manifest_buffer_bytes, limits.manifest);
    errdefer manifest.deinit();
    try validateConfigShape(runtime_config, manifest.shape);

    // Verify the complete generation before allocating or parsing model state.
    inline for (@typeInfo(manifest_module.Section).@"enum".fields) |field| try bundle_io.verifySectionFile(allocator, io, directory, manifest, @enumFromInt(field.value), settings.section_verify_buffer_bytes);

    var grid = try GridState.init(allocator, runtime_config);
    errdefer grid.deinit();
    var plants = try PlantState.init(allocator, runtime_config);
    errdefer plants.deinit();
    try readCoupled(io, directory, manifest, settings.section_read_buffer_bytes, &grid, &plants);

    var metadata = try readOwned(metadata_checkpoint.Metadata, allocator, io, directory, manifest, .plant_metadata, settings.section_read_buffer_bytes, limits.plant_metadata, metadata_checkpoint.read);
    errdefer metadata.deinit();
    if (metadata.year != manifest.instant.year or metadata.day_of_year != manifest.instant.day_of_year or metadata.cells.len != grid.cell_count) return error.CheckpointMetadataInstantMismatch;
    var development = try readOwned(development_checkpoint.Owned, allocator, io, directory, manifest, .plant_development, settings.section_read_buffer_bytes, limits.plant_development, development_checkpoint.read);
    errdefer development.deinit();
    var roots = try readOwned(RootState, allocator, io, directory, manifest, .plant_roots, settings.section_read_buffer_bytes, limits.plant_roots, root_checkpoint.read);
    errdefer roots.deinit();
    var canopy = try readOwned(canopy_checkpoint.Owned, allocator, io, directory, manifest, .plant_canopy, settings.section_read_buffer_bytes, limits.plant_canopy, canopy_checkpoint.read);
    errdefer canopy.deinit();
    var biogeochemistry = try readOwned(biogeochemistry_checkpoint.Owned, allocator, io, directory, manifest, .soil_biogeochemistry, settings.section_read_buffer_bytes, limits.soil_biogeochemistry, biogeochemistry_checkpoint.read);
    errdefer biogeochemistry.deinit();
    var organic = try readOwned(organic_checkpoint.Owned, allocator, io, directory, manifest, .soil_organic_matter, settings.section_read_buffer_bytes, limits.soil_organic, organic_checkpoint.read);
    errdefer organic.deinit();
    var transport = try readOwned(transport_checkpoint.Owned, allocator, io, directory, manifest, .solute_gas_and_snow_transport, settings.section_read_buffer_bytes, limits.transport, transport_checkpoint.read);
    errdefer transport.deinit();
    var geometry = try readOwned(geometry_checkpoint.Owned, allocator, io, directory, manifest, .soil_geometry_and_hydrology, settings.section_read_buffer_bytes, limits.soil_geometry, geometry_checkpoint.read);
    errdefer geometry.deinit();
    const landscape_mass_balance = try readMassBalance(
        allocator,
        io,
        directory,
        manifest,
        settings.section_read_buffer_bytes,
    );

    try validateOwnedShape(manifest.shape, grid, plants, development, roots, canopy, biogeochemistry, organic, transport, geometry);
    return .{ .manifest = manifest, .grid = grid, .plants = plants, .plant_metadata = metadata, .plant_development = development, .plant_roots = roots, .plant_canopy = canopy, .soil_biogeochemistry = biogeochemistry, .soil_organic_matter = organic, .transport = transport, .soil_geometry_and_hydrology = geometry, .landscape_mass_balance = landscape_mass_balance };
}

fn readManifest(allocator: std.mem.Allocator, io: std.Io, directory: std.Io.Dir, file_name: []const u8, buffer_bytes: usize, limits: manifest_module.Limits) !manifest_module.Manifest {
    if (!safeFileName(file_name)) return error.InvalidCheckpointFileName;
    const buffer = try allocator.alloc(u8, buffer_bytes);
    defer allocator.free(buffer);
    var file = try directory.openFile(io, file_name, .{});
    defer file.close(io);
    var file_reader = file.readerStreaming(io, buffer);
    return manifest_module.read(allocator, &file_reader.interface, limits);
}

fn safeFileName(name: []const u8) bool {
    if (name.len == 0 or std.mem.eql(u8, name, ".") or std.mem.eql(u8, name, ".."))
        return false;
    if (std.mem.indexOf(u8, name, "..") != null) return false;
    if (name[0] == ' ' or name[name.len - 1] == ' ' or name[name.len - 1] == '.')
        return false;
    for (name) |byte| {
        if (byte == 0 or byte < 0x20 or
            std.mem.indexOfScalar(u8, "/\\<>:\"|?*", byte) != null)
            return false;
    }
    return true;
}

fn readCoupled(io: std.Io, directory: std.Io.Dir, manifest: manifest_module.Manifest, buffer_bytes: usize, grid: *GridState, plants: *PlantState) !void {
    const entry = manifest.entry(.grid_and_plants) orelse return error.MissingCheckpointSection;
    const buffer = try grid.allocator.alloc(u8, buffer_bytes);
    defer grid.allocator.free(buffer);
    var file = try directory.openFile(io, entry.file_name, .{});
    defer file.close(io);
    var file_reader = file.readerStreaming(io, buffer);
    try checkpoint.readCoupledInto(&file_reader.interface, grid, plants);
}

fn readOwned(comptime T: type, allocator: std.mem.Allocator, io: std.Io, directory: std.Io.Dir, manifest: manifest_module.Manifest, section: manifest_module.Section, buffer_bytes: usize, limits: anytype, comptime readFunction: anytype) !T {
    const entry = manifest.entry(section) orelse return error.MissingCheckpointSection;
    const buffer = try allocator.alloc(u8, buffer_bytes);
    defer allocator.free(buffer);
    var file = try directory.openFile(io, entry.file_name, .{});
    defer file.close(io);
    var file_reader = file.readerStreaming(io, buffer);
    return try readFunction(allocator, &file_reader.interface, limits);
}

fn readMassBalance(
    allocator: std.mem.Allocator,
    io: std.Io,
    directory: std.Io.Dir,
    manifest: manifest_module.Manifest,
    buffer_bytes: usize,
) !mass_balance_checkpoint.State {
    const entry = manifest.entry(.landscape_mass_balance) orelse
        return error.MissingCheckpointSection;
    const buffer = try allocator.alloc(u8, buffer_bytes);
    defer allocator.free(buffer);
    var file = try directory.openFile(io, entry.file_name, .{});
    defer file.close(io);
    var file_reader = file.readerStreaming(io, buffer);
    return mass_balance_checkpoint.read(&file_reader.interface);
}

fn validateConfigShape(config: config_module.SimulationConfig, shape: manifest_module.RuntimeShape) !void {
    if (config.lon_count != shape.columns or config.lat_count != shape.rows or config.soil_layers != shape.soil_layers or config.plant_populations != shape.plant_species_per_cell) return error.CheckpointRuntimeConfigMismatch;
}

fn validateOwnedShape(shape: manifest_module.RuntimeShape, grid: GridState, plants: PlantState, development: development_checkpoint.Owned, roots: RootState, canopy: canopy_checkpoint.Owned, biogeochemistry: biogeochemistry_checkpoint.Owned, organic: organic_checkpoint.Owned, transport: transport_checkpoint.Owned, geometry: geometry_checkpoint.Owned) !void {
    const cells = try std.math.mul(usize, shape.columns, shape.rows);
    const plant_count = try std.math.mul(usize, cells, shape.plant_species_per_cell);
    const layer_cells = try std.math.mul(usize, cells, shape.soil_layers);
    if (grid.cell_count != cells or grid.soil_layer_capacity != shape.soil_layers or
        plants.cell_count != cells or plants.species_count != shape.plant_species_per_cell or
        development.phenology.cell_count != cells or development.phenology.species_count != shape.plant_species_per_cell or
        roots.plant_count != plant_count or roots.soil_layer_count != shape.soil_layers or roots.root_axis_count != shape.root_axes_per_plant or
        canopy.canopy.cell_count != cells or canopy.canopy.species_count != shape.plant_species_per_cell or canopy.retention.cell_count != cells or canopy.layer_distribution.cell_count != cells or canopy.layer_distribution.species_count != shape.plant_species_per_cell or
        biogeochemistry.chemistry.cell_count != layer_cells or biogeochemistry.available_nutrients.layer_count != layer_cells or biogeochemistry.fertilizer.cell_count != cells or biogeochemistry.fertilizer.layer_capacity != shape.soil_layers or biogeochemistry.fertilizer_band.cell_count != cells or biogeochemistry.fertilizer_band.layer_capacity != shape.soil_layers or
        organic.profile.layer_count != layer_cells or organic.surface.layer_count != cells or organic.litter_chemistry.cells.len != cells or organic.litter_fertilizer.cells.len != cells or organic.litter_fertilizer.formulation.len != cells or organic.surface_respiration.cell_count != cells or organic.surface_denitrification.cell_count != cells or organic.surface_fire_exchange.layer_count != cells or organic.litter_salt_ingress.cell_count != cells or organic.litter_salt_ingress.soil_layer_capacity != shape.soil_layers or
        transport.micropore.cell_count != layer_cells or transport.mineral_nitrogen.cell_count != layer_cells or transport.organic.layer_count != layer_cells or transport.litter_gas.cell_count != cells or transport.snow.cell_count != cells or transport.snow.layer_capacity != shape.snow_layers or
        transport.surface.columns != shape.columns or transport.surface.rows != shape.rows or
        geometry.geometry.cell_count != cells or geometry.geometry.layer_capacity != shape.soil_layers or
        geometry.hydrology.columns != shape.columns or geometry.hydrology.rows != shape.rows or geometry.hydrology.snow_layer_capacity != shape.snow_layers or geometry.surface.cell_count != cells or geometry.surface_litter_geometry.cell_count != cells or geometry.erosion.cell_count != cells or geometry.eroded_minerals.workspace.cell_count != cells or geometry.runtime.layer_count != layer_cells) return error.CheckpointOwnedShapeMismatch;
}

fn validateSwapTargets(bundle: OwnedBundle, targets: LiveTargets) !void {
    const shape = bundle.manifest.shape;
    const cells = try std.math.mul(usize, shape.columns, shape.rows);
    const plants = try std.math.mul(usize, cells, shape.plant_species_per_cell);
    const layer_cells = try std.math.mul(usize, cells, shape.soil_layers);
    if (targets.grid.cell_count != cells or targets.grid.soil_layer_capacity != shape.soil_layers or
        targets.plants.cell_count != cells or targets.plants.species_count != shape.plant_species_per_cell or
        targets.plant_development.phenology.cell_count != cells or targets.plant_development.phenology.species_count != shape.plant_species_per_cell or
        targets.plant_roots.plant_count != plants or targets.plant_roots.soil_layer_count != shape.soil_layers or targets.plant_roots.root_axis_count != shape.root_axes_per_plant or
        targets.plant_canopy.canopy.cell_count != cells or targets.plant_canopy.canopy.species_count != shape.plant_species_per_cell or targets.plant_canopy.retention.cell_count != cells or targets.plant_canopy.layer_distribution.cell_count != cells or targets.plant_canopy.layer_distribution.species_count != shape.plant_species_per_cell or
        targets.soil_biogeochemistry.chemistry.cell_count != layer_cells or targets.soil_biogeochemistry.available_nutrients.layer_count != layer_cells or targets.soil_biogeochemistry.fertilizer.cell_count != cells or targets.soil_biogeochemistry.fertilizer.layer_capacity != shape.soil_layers or targets.soil_biogeochemistry.mineral_fertilizer.cell_count != cells or targets.soil_biogeochemistry.mineral_fertilizer.layer_capacity != shape.soil_layers or targets.soil_biogeochemistry.fertilizer_band.cell_count != cells or targets.soil_biogeochemistry.fertilizer_band.layer_capacity != shape.soil_layers or targets.soil_biogeochemistry.reactive_nitrogen.layer_count != layer_cells or targets.soil_biogeochemistry.microbial_phosphorus.layer_count != layer_cells or
        targets.soil_organic_matter.profile.layer_count != layer_cells or targets.soil_organic_matter.surface.layer_count != cells or targets.soil_organic_matter.litter_chemistry.cells.len != cells or targets.soil_organic_matter.litter_fertilizer.cells.len != cells or targets.soil_organic_matter.litter_fertilizer.formulation.len != cells or targets.soil_organic_matter.surface_respiration.cell_count != cells or targets.soil_organic_matter.surface_denitrification.cell_count != cells or targets.soil_organic_matter.surface_fire_exchange.layer_count != cells or targets.soil_organic_matter.litter_salt_ingress.cell_count != cells or targets.soil_organic_matter.litter_salt_ingress.soil_layer_capacity != shape.soil_layers or
        targets.transport.micropore.cell_count != layer_cells or targets.transport.mineral_nitrogen.cell_count != layer_cells or targets.transport.organic.layer_count != layer_cells or targets.transport.litter_gas.cell_count != cells or targets.transport.snow.cell_count != cells or targets.transport.snow.layer_capacity != shape.snow_layers or
        targets.transport.surface.columns != shape.columns or targets.transport.surface.rows != shape.rows or
        targets.soil_geometry_and_hydrology.geometry.cell_count != cells or targets.soil_geometry_and_hydrology.geometry.layer_capacity != shape.soil_layers or
        targets.soil_geometry_and_hydrology.hydrology.columns != shape.columns or targets.soil_geometry_and_hydrology.hydrology.rows != shape.rows or targets.soil_geometry_and_hydrology.hydrology.snow_layer_capacity != shape.snow_layers or targets.soil_geometry_and_hydrology.surface.cell_count != cells or targets.soil_geometry_and_hydrology.surface_litter_geometry.cell_count != cells or targets.soil_geometry_and_hydrology.erosion.cell_count != cells or targets.soil_geometry_and_hydrology.eroded_minerals.workspace.cell_count != cells or targets.soil_geometry_and_hydrology.delayed_live_canopy_combustion_heat_megajoules.len != plants or targets.soil_geometry_and_hydrology.delayed_standing_dead_combustion_heat_megajoules.len != plants or targets.soil_geometry_and_hydrology.delayed_subsurface_combustion_heat_megajoules.len != layer_cells or targets.soil_geometry_and_hydrology.delayed_surface_combustion_heat_megajoules.len != cells or targets.soil_geometry_and_hydrology.surface_litter_ice_m3.len != cells) return error.CheckpointLiveSwapShapeMismatch;
    const runtime = targets.soil_geometry_and_hydrology.runtime orelse
        return error.MissingSoilRuntimeCheckpointTarget;
    try @import("soil_runtime_checkpoint.zig").validateView(runtime);
    try bundle.soil_geometry_and_hydrology.runtime.validate();
    const surface_boundary =
        targets.soil_geometry_and_hydrology.surface_boundary orelse
        return error.MissingSurfaceBoundaryCheckpointTarget;
    try @import("surface_boundary_checkpoint.zig").validateTargetDimensions(
        surface_boundary,
    );
    try bundle.soil_geometry_and_hydrology.surface_boundary.validate();
}

test "bundle reader rejects runtime configuration before model owner allocation" {
    const config = try config_module.SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 3 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    try std.testing.expectError(error.CheckpointRuntimeConfigMismatch, validateConfigShape(config, .{ .columns = 1, .rows = 1, .soil_layers = 2, .snow_layers = 2, .plant_species_per_cell = 4, .root_axes_per_plant = 10 }));
}

test "bundle reader requires an atomic manifest before any restore" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const config = try config_module.SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    const manifest_limits: manifest_module.Limits = .{ .maximum_columns = 1, .maximum_rows = 1, .maximum_soil_layers = 1, .maximum_snow_layers = 1, .maximum_plant_species_per_cell = 1, .maximum_root_axes_per_plant = 1 };
    const limits: Limits = .{
        .manifest = manifest_limits,
        .plant_metadata = .{ .maximum_cells = 1, .maximum_species_per_cell = 1, .maximum_species_name_bytes = 32 },
        .plant_development = .{ .maximum_cells = 1, .maximum_species = 1, .maximum_branches = 1 },
        .plant_roots = .{ .maximum_plants = 1, .maximum_soil_layers = 1, .maximum_root_axes = 1 },
        .plant_canopy = .{ .maximum_cells = 1, .maximum_species = 1, .maximum_branches = 1, .maximum_nodes = 1, .maximum_samples = 1, .maximum_layers = 1, .maximum_inclinations = 1, .maximum_azimuths = 1 },
        .soil_biogeochemistry = .{ .maximum_cells = 1, .maximum_layers = 1, .maximum_substrates = 1, .maximum_populations = 1 },
        .soil_organic = .{ .maximum_profile_layers = 1, .maximum_surface_cells = 1 },
        .transport = .{ .maximum_transport_cells = 1, .maximum_solute_species = 1, .maximum_snow_cells = 1, .maximum_snow_layers = 1 },
        .soil_geometry = .{ .maximum_columns = 1, .maximum_rows = 1, .maximum_soil_layers = 1, .maximum_snow_layers = 1, .maximum_plants = 1 },
    };
    try std.testing.expectError(error.FileNotFound, read(std.testing.allocator, std.testing.io, temporary.dir, config, limits, .{ .manifest_buffer_bytes = 64, .section_read_buffer_bytes = 64, .section_verify_buffer_bytes = 64 }));
}

test "bundle reader rejects unsafe manifest names before allocation or I/O" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const limits: manifest_module.Limits = .{
        .maximum_columns = 1,
        .maximum_rows = 1,
        .maximum_soil_layers = 1,
        .maximum_snow_layers = 1,
        .maximum_plant_species_per_cell = 1,
        .maximum_root_axes_per_plant = 1,
    };
    inline for (.{
        "",
        "../restart.manifest",
        "subdir/restart.manifest",
        "subdir\\restart.manifest",
        " restart.manifest",
        "restart.manifest ",
        "restart.",
        "restart:1.manifest",
        "restart|1.manifest",
        "restart?1.manifest",
    }) |name| try std.testing.expectError(
        error.InvalidCheckpointFileName,
        readManifest(
            std.testing.allocator,
            std.testing.io,
            temporary.dir,
            name,
            64,
            limits,
        ),
    );
}

test "bundle reader accepts portable manifest names at the I/O boundary" {
    try std.testing.expect(safeFileName("restart"));
    try std.testing.expect(safeFileName("restart-1.manifest"));
    try std.testing.expect(safeFileName("restart.1.manifest"));
}
