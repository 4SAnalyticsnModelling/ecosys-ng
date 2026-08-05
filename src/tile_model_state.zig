const std = @import("std");
const compute = @import("compute.zig");
const grid_module = @import("grid.zig");
const soil_thermal = @import("soil_thermal.zig");
const transport_hydrology = @import("transport_hydrology.zig");
const gas_transport = @import("gas_transport.zig");
const dissolved_gas_transport = @import("soil_dissolved_gas_transport.zig");
const solute_transport = @import("solute_transport.zig");
const spatial_grid = @import("spatial_grid.zig");
const tile_executor = @import("tile_executor.zig");
const tile_state_store = @import("tile_state_store.zig");

fn countSliceFields(comptime State: type, comptime Slice: type) usize {
    var count: usize = 0;
    for (@typeInfo(State).@"struct".fields) |declared|
        count += @intFromBool(declared.type == Slice);
    return count;
}

const base_field_count: usize = 30;
const hydrology_field_count: usize = countSliceFields(
    transport_hydrology.State,
    []f64,
);
const gas_field_count: usize = countSliceFields(gas_transport.State, []f64);
const solute_field_count: usize = countSliceFields(
    solute_transport.State,
    []f64,
);
pub const field_count: usize =
    base_field_count +
    hydrology_field_count +
    gas_field_count * 2 +
    1 +
    solute_field_count * 2;

/// Authoritative runtime-shaped fields owned by the core grid and plant
/// states. Field order is the on-disk schema and therefore intentionally
/// explicit. Layer and species counts remain runtime values.
pub fn fields(
    grid: *const grid_module.GridState,
    plants: *const grid_module.PlantState,
    thermal: *const soil_thermal.State,
    hydrology: *const transport_hydrology.State,
    soil_gas: *const gas_transport.State,
    litter_gas: *const gas_transport.State,
    dissolved_gas: *const dissolved_gas_transport.State,
    micropore_solutes: *const solute_transport.State,
    macropore_solutes: *const solute_transport.State,
) ![field_count]tile_state_store.Field {
    try validateShapes(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    const soil_layers = grid.soil_layer_capacity;
    const species = plants.species_count;
    const root_components = std.math.mul(
        usize,
        species,
        plants.soil_layer_count,
    ) catch return error.TileModelStateDimensionOverflow;
    const base = [_]tile_state_store.Field{
        unsignedField(
            "active_soil_layer_count",
            grid.active_soil_layer_count,
            1,
        ),
        field("soil_temperature_k", grid.soil_temperature_k, soil_layers),
        field("liquid_water_m3", grid.liquid_water_m3, soil_layers),
        field("ice_water_m3", grid.ice_water_m3, soil_layers),
        field("matrix_liquid_water_m3", grid.matrix_liquid_water_m3, soil_layers),
        field("macropore_liquid_water_m3", grid.macropore_liquid_water_m3, soil_layers),
        field("matrix_ice_water_m3", grid.matrix_ice_water_m3, soil_layers),
        field("macropore_ice_water_m3", grid.macropore_ice_water_m3, soil_layers),
        field("matrix_pore_capacity_m3", grid.matrix_pore_capacity_m3, soil_layers),
        field("macropore_pore_capacity_m3", grid.macropore_pore_capacity_m3, soil_layers),
        field("matrix_air_volume_m3", grid.matrix_air_volume_m3, soil_layers),
        field("macropore_air_volume_m3", grid.macropore_air_volume_m3, soil_layers),
        field("air_volume_m3", grid.air_volume_m3, soil_layers),
        field("water_vapor_volume_m3", grid.water_vapor_volume_m3, soil_layers),
        field("matric_potential_mpa", grid.matric_potential_mpa, soil_layers),
        field("surface_temperature_k", grid.surface_temperature_k, 1),
        field("canopy_temperature_k", plants.canopy_temperature_k, species),
        field("canopy_water_potential_mpa", plants.canopy_water_potential_mpa, species),
        field("canopy_water_storage_m_per_m2", plants.canopy_water_storage_m_per_m2, species),
        field("leaf_area_index_m2_m2", plants.leaf_area_index_m2_m2, species),
        field("shoot_carbon_g_m2", plants.shoot_carbon_g_m2, species),
        field("root_carbon_g_m2", plants.root_carbon_g_m2, root_components),
        field("soil_thermal_layer_volume_m3", thermal.layer_volume_m3, soil_layers),
        field("soil_thermal_layer_thickness_m", thermal.layer_thickness_m, soil_layers),
        field("soil_thermal_porosity_fraction", thermal.porosity_fraction, soil_layers),
        field("soil_thermal_dry_solid_heat_capacity_megajoules_per_m3_k", thermal.dry_solid_heat_capacity_megajoules_per_m3_k, soil_layers),
        field("soil_thermal_solid_conductivity_numerator_m_megajoules_per_h_k", thermal.solid_thermal_conductivity_numerator_m_megajoules_per_h_k, soil_layers),
        field("soil_thermal_solid_conductivity_denominator", thermal.solid_thermal_conductivity_denominator, soil_layers),
        field("soil_thermal_total_heat_capacity_megajoules_per_m3_k", thermal.total_heat_capacity_megajoules_per_m3_k, soil_layers),
        field("soil_thermal_conductivity_m_megajoules_per_h_k", thermal.thermal_conductivity_m_megajoules_per_h_k, soil_layers),
    };
    var result: [field_count]tile_state_store.Field = undefined;
    @memcpy(result[0..base_field_count], &base);
    var next = base_field_count;
    inline for (@typeInfo(transport_hydrology.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = field(
                "transport_hydrology_" ++ declared.name,
                @field(hydrology, declared.name),
                @field(hydrology, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    inline for (@typeInfo(gas_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = field(
                "soil_gas_" ++ declared.name,
                @field(soil_gas, declared.name),
                @field(soil_gas, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    inline for (@typeInfo(gas_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = field(
                "litter_gas_" ++ declared.name,
                @field(litter_gas, declared.name),
                @field(litter_gas, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    result[next] = field(
        "dissolved_gas_boundary_net_flux_g",
        dissolved_gas.boundary_net_flux_g,
        dissolved_gas.boundary_net_flux_g.len / grid.cell_count,
    );
    next += 1;
    inline for (@typeInfo(solute_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = field(
                "micropore_solute_" ++ declared.name,
                @field(micropore_solutes, declared.name),
                @field(micropore_solutes, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    inline for (@typeInfo(solute_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = field(
                "macropore_solute_" ++ declared.name,
                @field(macropore_solutes, declared.name),
                @field(macropore_solutes, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    std.debug.assert(next == field_count);
    return result;
}

pub fn mutableFields(
    grid: *grid_module.GridState,
    plants: *grid_module.PlantState,
    thermal: *soil_thermal.State,
    hydrology: *transport_hydrology.State,
    soil_gas: *gas_transport.State,
    litter_gas: *gas_transport.State,
    dissolved_gas: *dissolved_gas_transport.State,
    micropore_solutes: *solute_transport.State,
    macropore_solutes: *solute_transport.State,
) ![field_count]tile_state_store.MutableField {
    try validateShapes(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    const soil_layers = grid.soil_layer_capacity;
    const species = plants.species_count;
    const root_components = std.math.mul(
        usize,
        species,
        plants.soil_layer_count,
    ) catch return error.TileModelStateDimensionOverflow;
    const base = [_]tile_state_store.MutableField{
        mutableUnsignedField(
            "active_soil_layer_count",
            grid.active_soil_layer_count,
            1,
        ),
        mutableField("soil_temperature_k", grid.soil_temperature_k, soil_layers),
        mutableField("liquid_water_m3", grid.liquid_water_m3, soil_layers),
        mutableField("ice_water_m3", grid.ice_water_m3, soil_layers),
        mutableField("matrix_liquid_water_m3", grid.matrix_liquid_water_m3, soil_layers),
        mutableField("macropore_liquid_water_m3", grid.macropore_liquid_water_m3, soil_layers),
        mutableField("matrix_ice_water_m3", grid.matrix_ice_water_m3, soil_layers),
        mutableField("macropore_ice_water_m3", grid.macropore_ice_water_m3, soil_layers),
        mutableField("matrix_pore_capacity_m3", grid.matrix_pore_capacity_m3, soil_layers),
        mutableField("macropore_pore_capacity_m3", grid.macropore_pore_capacity_m3, soil_layers),
        mutableField("matrix_air_volume_m3", grid.matrix_air_volume_m3, soil_layers),
        mutableField("macropore_air_volume_m3", grid.macropore_air_volume_m3, soil_layers),
        mutableField("air_volume_m3", grid.air_volume_m3, soil_layers),
        mutableField("water_vapor_volume_m3", grid.water_vapor_volume_m3, soil_layers),
        mutableField("matric_potential_mpa", grid.matric_potential_mpa, soil_layers),
        mutableField("surface_temperature_k", grid.surface_temperature_k, 1),
        mutableField("canopy_temperature_k", plants.canopy_temperature_k, species),
        mutableField("canopy_water_potential_mpa", plants.canopy_water_potential_mpa, species),
        mutableField("canopy_water_storage_m_per_m2", plants.canopy_water_storage_m_per_m2, species),
        mutableField("leaf_area_index_m2_m2", plants.leaf_area_index_m2_m2, species),
        mutableField("shoot_carbon_g_m2", plants.shoot_carbon_g_m2, species),
        mutableField(
            "root_carbon_g_m2",
            plants.root_carbon_g_m2,
            root_components,
        ),
        mutableField("soil_thermal_layer_volume_m3", thermal.layer_volume_m3, soil_layers),
        mutableField("soil_thermal_layer_thickness_m", thermal.layer_thickness_m, soil_layers),
        mutableField("soil_thermal_porosity_fraction", thermal.porosity_fraction, soil_layers),
        mutableField("soil_thermal_dry_solid_heat_capacity_megajoules_per_m3_k", thermal.dry_solid_heat_capacity_megajoules_per_m3_k, soil_layers),
        mutableField("soil_thermal_solid_conductivity_numerator_m_megajoules_per_h_k", thermal.solid_thermal_conductivity_numerator_m_megajoules_per_h_k, soil_layers),
        mutableField("soil_thermal_solid_conductivity_denominator", thermal.solid_thermal_conductivity_denominator, soil_layers),
        mutableField("soil_thermal_total_heat_capacity_megajoules_per_m3_k", thermal.total_heat_capacity_megajoules_per_m3_k, soil_layers),
        mutableField("soil_thermal_conductivity_m_megajoules_per_h_k", thermal.thermal_conductivity_m_megajoules_per_h_k, soil_layers),
    };
    var result: [field_count]tile_state_store.MutableField = undefined;
    @memcpy(result[0..base_field_count], &base);
    var next = base_field_count;
    inline for (@typeInfo(transport_hydrology.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = mutableField(
                "transport_hydrology_" ++ declared.name,
                @field(hydrology, declared.name),
                @field(hydrology, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    inline for (@typeInfo(gas_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = mutableField(
                "soil_gas_" ++ declared.name,
                @field(soil_gas, declared.name),
                @field(soil_gas, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    inline for (@typeInfo(gas_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = mutableField(
                "litter_gas_" ++ declared.name,
                @field(litter_gas, declared.name),
                @field(litter_gas, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    result[next] = mutableField(
        "dissolved_gas_boundary_net_flux_g",
        dissolved_gas.boundary_net_flux_g,
        dissolved_gas.boundary_net_flux_g.len / grid.cell_count,
    );
    next += 1;
    inline for (@typeInfo(solute_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = mutableField(
                "micropore_solute_" ++ declared.name,
                @field(micropore_solutes, declared.name),
                @field(micropore_solutes, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    inline for (@typeInfo(solute_transport.State).@"struct".fields) |declared| {
        if (declared.type == []f64) {
            result[next] = mutableField(
                "macropore_solute_" ++ declared.name,
                @field(macropore_solutes, declared.name),
                @field(macropore_solutes, declared.name).len / grid.cell_count,
            );
            next += 1;
        }
    }
    std.debug.assert(next == field_count);
    return result;
}

pub fn saveTile(
    store: tile_state_store.FileStore,
    tile: spatial_grid.Tile,
    lat_count: usize,
    lon_count: usize,
    grid: *const grid_module.GridState,
    plants: *const grid_module.PlantState,
    thermal: *const soil_thermal.State,
    hydrology: *const transport_hydrology.State,
    soil_gas: *const gas_transport.State,
    litter_gas: *const gas_transport.State,
    dissolved_gas: *const dissolved_gas_transport.State,
    micropore_solutes: *const solute_transport.State,
    macropore_solutes: *const solute_transport.State,
) !void {
    const catalog = try fields(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    try store.saveOwnedFields(
        tile,
        lat_count,
        lon_count,
        &catalog,
    );
}

/// Creates the authoritative owned-interior file for every tile before the
/// in-memory domain is released or reused as a tile working set.
pub fn initializeStore(
    store: tile_state_store.FileStore,
    plan: spatial_grid.TilePlan,
    grid: *const grid_module.GridState,
    plants: *const grid_module.PlantState,
    thermal: *const soil_thermal.State,
    hydrology: *const transport_hydrology.State,
    soil_gas: *const gas_transport.State,
    litter_gas: *const gas_transport.State,
    dissolved_gas: *const dissolved_gas_transport.State,
    micropore_solutes: *const solute_transport.State,
    macropore_solutes: *const solute_transport.State,
) !void {
    try store.requireUnpublished();
    for (plan.tiles) |tile| try saveTile(
        store,
        tile,
        plan.lat_count,
        plan.lon_count,
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    const catalog = try fields(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    try store.publishGeneration(plan, &catalog);
}

pub fn loadTile(
    store: tile_state_store.FileStore,
    plan: spatial_grid.TilePlan,
    tile: spatial_grid.Tile,
    grid: *grid_module.GridState,
    plants: *grid_module.PlantState,
    thermal: *soil_thermal.State,
    hydrology: *transport_hydrology.State,
    soil_gas: *gas_transport.State,
    litter_gas: *gas_transport.State,
    dissolved_gas: *dissolved_gas_transport.State,
    micropore_solutes: *solute_transport.State,
    macropore_solutes: *solute_transport.State,
) !void {
    const catalog = try mutableFields(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    try store.loadTileFieldsInto(
        plan,
        tile,
        &catalog,
    );
    try grid.validateFinite();
    try plants.validateFinite();
    try thermal.validateFinite();
    try hydrology.validateFinite();
    try soil_gas.validateFinite();
    try litter_gas.validateFinite();
    try dissolved_gas.validateFinite();
    try micropore_solutes.validateFinite();
    try macropore_solutes.validateFinite();
}

/// Runs tiles strictly serially between two generations. Every tile loads its
/// authoritative interior and two-cell halo exclusively from the immutable
/// source generation. Its owned cells are split across CPU workers and, after
/// all workers join and validation succeeds, atomically committed to the
/// destination generation. The caller swaps generations only after this
/// function completes, preventing traversal-order-dependent same-step reads.
pub fn runSerialTiles(
    executor: compute.CpuExecutor,
    source_store: tile_state_store.FileStore,
    destination_store: tile_state_store.FileStore,
    plan: spatial_grid.TilePlan,
    grid: *grid_module.GridState,
    plants: *grid_module.PlantState,
    thermal: *soil_thermal.State,
    hydrology: *transport_hydrology.State,
    soil_gas: *gas_transport.State,
    litter_gas: *gas_transport.State,
    dissolved_gas: *dissolved_gas_transport.State,
    micropore_solutes: *solute_transport.State,
    macropore_solutes: *solute_transport.State,
    science_context: anytype,
    comptime cell_kernel: anytype,
) !void {
    if (source_store.generation_index == destination_store.generation_index)
        return error.TileGenerationMustBeDoubleBuffered;
    if (destination_store.generation_index !=
        std.math.add(
            u64,
            source_store.generation_index,
            1,
        ) catch return error.TileGenerationIndexOverflow)
        return error.TileGenerationIsNotConsecutive;
    const source_catalog = try mutableFields(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    try source_store.validateGeneration(plan, &source_catalog);
    try destination_store.requireUnpublished();
    const Context = struct {
        source_store: tile_state_store.FileStore,
        destination_store: tile_state_store.FileStore,
        plan: spatial_grid.TilePlan,
        grid: *grid_module.GridState,
        plants: *grid_module.PlantState,
        thermal: *soil_thermal.State,
        hydrology: *transport_hydrology.State,
        soil_gas: *gas_transport.State,
        litter_gas: *gas_transport.State,
        dissolved_gas: *dissolved_gas_transport.State,
        micropore_solutes: *solute_transport.State,
        macropore_solutes: *solute_transport.State,
        science: @TypeOf(science_context),
    };
    const Adapter = struct {
        fn load(context: *Context, tile: spatial_grid.Tile) !void {
            try loadTile(
                context.source_store,
                context.plan,
                tile,
                context.grid,
                context.plants,
                context.thermal,
                context.hydrology,
                context.soil_gas,
                context.litter_gas,
                context.dissolved_gas,
                context.micropore_solutes,
                context.macropore_solutes,
            );
        }

        fn computeCells(context: *Context, cells: []const usize) !void {
            try cell_kernel(context.science, cells);
        }

        fn commit(context: *Context, tile: spatial_grid.Tile) !void {
            // All workers have joined before commit is called, so validation
            // cannot race a worker updating another owned cell.
            try context.grid.validateFinite();
            try context.plants.validateFinite();
            try context.thermal.validateFinite();
            try context.hydrology.validateFinite();
            try context.soil_gas.validateFinite();
            try context.litter_gas.validateFinite();
            try context.dissolved_gas.validateFinite();
            try context.micropore_solutes.validateFinite();
            try context.macropore_solutes.validateFinite();
            try saveTile(
                context.destination_store,
                tile,
                context.plan.lat_count,
                context.plan.lon_count,
                context.grid,
                context.plants,
                context.thermal,
                context.hydrology,
                context.soil_gas,
                context.litter_gas,
                context.dissolved_gas,
                context.micropore_solutes,
                context.macropore_solutes,
            );
        }
    };
    var context = Context{
        .source_store = source_store,
        .destination_store = destination_store,
        .plan = plan,
        .grid = grid,
        .plants = plants,
        .thermal = thermal,
        .hydrology = hydrology,
        .soil_gas = soil_gas,
        .litter_gas = litter_gas,
        .dissolved_gas = dissolved_gas,
        .micropore_solutes = micropore_solutes,
        .macropore_solutes = macropore_solutes,
        .science = science_context,
    };
    try tile_executor.run(
        executor,
        plan,
        &context,
        Adapter.load,
        Adapter.computeCells,
        Adapter.commit,
    );
    const destination_catalog = try fields(
        grid,
        plants,
        thermal,
        hydrology,
        soil_gas,
        litter_gas,
        dissolved_gas,
        micropore_solutes,
        macropore_solutes,
    );
    try destination_store.publishGeneration(plan, &destination_catalog);
}

fn field(
    name: []const u8,
    values: []const f64,
    components_per_cell: usize,
) tile_state_store.Field {
    return .{
        .name = name,
        .values = values,
        .components_per_cell = components_per_cell,
    };
}

fn mutableField(
    name: []const u8,
    values: []f64,
    components_per_cell: usize,
) tile_state_store.MutableField {
    return .{
        .name = name,
        .values = values,
        .components_per_cell = components_per_cell,
    };
}

fn unsignedField(
    name: []const u8,
    values: []const usize,
    components_per_cell: usize,
) tile_state_store.Field {
    return .{
        .name = name,
        .unsigned_values = values,
        .components_per_cell = components_per_cell,
    };
}

fn mutableUnsignedField(
    name: []const u8,
    values: []usize,
    components_per_cell: usize,
) tile_state_store.MutableField {
    return .{
        .name = name,
        .unsigned_values = values,
        .components_per_cell = components_per_cell,
    };
}

fn validateShapes(
    grid: *const grid_module.GridState,
    plants: *const grid_module.PlantState,
    thermal: *const soil_thermal.State,
    hydrology: *const transport_hydrology.State,
    soil_gas: *const gas_transport.State,
    litter_gas: *const gas_transport.State,
    dissolved_gas: *const dissolved_gas_transport.State,
    micropore_solutes: *const solute_transport.State,
    macropore_solutes: *const solute_transport.State,
) !void {
    if (grid.cell_count == 0 or plants.cell_count != grid.cell_count or
        plants.soil_layer_count != grid.soil_layer_capacity or
        thermal.cell_count != grid.cell_count or
        thermal.soil_layer_capacity != grid.soil_layer_capacity or
        thermal.layer_volume_m3.len != grid.layer_count or
        hydrology.columns * hydrology.rows != grid.cell_count or
        hydrology.soil_layer_capacity != grid.soil_layer_capacity or
        soil_gas.cell_count != grid.layer_count or
        litter_gas.cell_count != grid.cell_count or
        dissolved_gas.layer_count != grid.layer_count or
        micropore_solutes.cell_count != grid.layer_count or
        macropore_solutes.cell_count != grid.layer_count or
        micropore_solutes.species_count != macropore_solutes.species_count)
        return error.TileModelStateDimensionMismatch;
    try micropore_solutes.validateFinite();
    try macropore_solutes.validateFinite();
}

fn initTestThermal(
    allocator: std.mem.Allocator,
    grid: grid_module.GridState,
) !soil_thermal.State {
    const layer_volume_m3 = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(layer_volume_m3);
    const layer_thickness_m = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(layer_thickness_m);
    const porosity_fraction = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(porosity_fraction);
    const dry_solid_heat_capacity = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(dry_solid_heat_capacity);
    const conductivity_numerator = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(conductivity_numerator);
    const conductivity_denominator = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(conductivity_denominator);
    const total_heat_capacity = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(total_heat_capacity);
    const thermal_conductivity = try allocator.alloc(f64, grid.layer_count);
    errdefer allocator.free(thermal_conductivity);
    @memset(layer_volume_m3, 1);
    @memset(layer_thickness_m, 0.1);
    @memset(porosity_fraction, 0.5);
    @memset(dry_solid_heat_capacity, 1);
    @memset(conductivity_numerator, 1);
    @memset(conductivity_denominator, 1);
    @memset(total_heat_capacity, 1);
    @memset(thermal_conductivity, 1);
    return .{
        .allocator = allocator,
        .cell_count = grid.cell_count,
        .soil_layer_capacity = grid.soil_layer_capacity,
        .layer_volume_m3 = layer_volume_m3,
        .layer_thickness_m = layer_thickness_m,
        .porosity_fraction = porosity_fraction,
        .dry_solid_heat_capacity_megajoules_per_m3_k = dry_solid_heat_capacity,
        .solid_thermal_conductivity_numerator_m_megajoules_per_h_k = conductivity_numerator,
        .solid_thermal_conductivity_denominator = conductivity_denominator,
        .total_heat_capacity_megajoules_per_m3_k = total_heat_capacity,
        .thermal_conductivity_m_megajoules_per_h_k = thermal_conductivity,
    };
}

test "production grid and runtime plant fields round trip through one Morton tile" {
    try std.testing.expectEqual(@as(usize, 81), field_count);
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 3,
            .lat_count = 2,
            .soil_layers = 4,
            .plant_populations = 7,
        },
        .{ .worker_threads = 2, .tile_cells = 6 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var plants = try grid_module.PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var thermal = try initTestThermal(std.testing.allocator, grid);
    defer thermal.deinit();
    var hydrology = try transport_hydrology.State.init(
        std.testing.allocator,
        config.lon_count,
        config.lat_count,
        config.soil_layers,
        3,
    );
    defer hydrology.deinit();
    var soil_gas = try gas_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer soil_gas.deinit();
    var litter_gas = try gas_transport.State.init(
        std.testing.allocator,
        grid.cell_count,
    );
    defer litter_gas.deinit();
    @memset(soil_gas.temperature_k, 273.15);
    @memset(litter_gas.temperature_k, 273.15);
    var dissolved_gas = try dissolved_gas_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer dissolved_gas.deinit();
    const runtime_aqueous_species_count: usize = 13;
    var micropore_solutes = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        runtime_aqueous_species_count,
    );
    defer micropore_solutes.deinit();
    var macropore_solutes = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        runtime_aqueous_species_count,
    );
    defer macropore_solutes.deinit();
    for (grid.matrix_liquid_water_m3, 0..) |*value, index|
        value.* = @floatFromInt(index + 100);
    for (plants.root_carbon_g_m2, 0..) |*value, index|
        value.* = @floatFromInt(index + 500);
    for (grid.active_soil_layer_count, 0..) |*value, cell|
        value.* = 1 + cell % config.soil_layers;
    for (thermal.thermal_conductivity_m_megajoules_per_h_k, 0..) |*value, index|
        value.* = @floatFromInt(index + 900);
    for (hydrology.micropore_face_flux_m3_per_step, 0..) |*value, index|
        value.* = -@as(f64, @floatFromInt(index + 1));
    for (soil_gas.dissolved_mass_g, 0..) |*value, index|
        value.* = @floatFromInt(index + 1000);
    for (litter_gas.gaseous_mass_g, 0..) |*value, index|
        value.* = @floatFromInt(index + 2000);
    for (dissolved_gas.boundary_net_flux_g, 0..) |*value, index|
        value.* = -@as(f64, @floatFromInt(index + 3000));
    for (micropore_solutes.water_volume_m3, 0..) |*value, index|
        value.* = @floatFromInt(index + 1);
    for (macropore_solutes.water_volume_m3, 0..) |*value, index|
        value.* = @floatFromInt(index + 2);
    for (micropore_solutes.amount_mol, 0..) |*value, index|
        value.* = @as(f64, @floatFromInt(index + 4000)) * 1e-6;
    for (macropore_solutes.amount_mol, 0..) |*value, index|
        value.* = @as(f64, @floatFromInt(index + 5000)) * 1e-6;
    const expected_thermal_conductivity = try std.testing.allocator.dupe(
        f64,
        thermal.thermal_conductivity_m_megajoules_per_h_k,
    );
    defer std.testing.allocator.free(expected_thermal_conductivity);
    const expected_micropore_face_flux = try std.testing.allocator.dupe(
        f64,
        hydrology.micropore_face_flux_m3_per_step,
    );
    defer std.testing.allocator.free(expected_micropore_face_flux);
    const expected_soil_dissolved_gas = try std.testing.allocator.dupe(
        f64,
        soil_gas.dissolved_mass_g,
    );
    defer std.testing.allocator.free(expected_soil_dissolved_gas);
    const expected_litter_gaseous_mass = try std.testing.allocator.dupe(
        f64,
        litter_gas.gaseous_mass_g,
    );
    defer std.testing.allocator.free(expected_litter_gaseous_mass);
    const expected_dissolved_boundary_flux = try std.testing.allocator.dupe(
        f64,
        dissolved_gas.boundary_net_flux_g,
    );
    defer std.testing.allocator.free(expected_dissolved_boundary_flux);
    const expected_micropore_solute_amount = try std.testing.allocator.dupe(
        f64,
        micropore_solutes.amount_mol,
    );
    defer std.testing.allocator.free(expected_micropore_solute_amount);
    const expected_macropore_solute_amount = try std.testing.allocator.dupe(
        f64,
        macropore_solutes.amount_mol,
    );
    defer std.testing.allocator.free(expected_macropore_solute_amount);
    const expected_active_layers = try std.testing.allocator.dupe(
        usize,
        grid.active_soil_layer_count,
    );
    defer std.testing.allocator.free(expected_active_layers);
    const expected_water =
        try std.testing.allocator.dupe(f64, grid.matrix_liquid_water_m3);
    defer std.testing.allocator.free(expected_water);
    const expected_roots =
        try std.testing.allocator.dupe(f64, plants.root_carbon_g_m2);
    defer std.testing.allocator.free(expected_roots);

    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const store = try tile_state_store.FileStore.init(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        4096,
    );
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        2,
        3,
        2,
        3,
        2,
    );
    defer plan.deinit();
    try initializeStore(
        store,
        plan,
        &grid,
        &plants,
        &thermal,
        &hydrology,
        &soil_gas,
        &litter_gas,
        &dissolved_gas,
        &micropore_solutes,
        &macropore_solutes,
    );
    @memset(grid.active_soil_layer_count, 0);
    @memset(grid.matrix_liquid_water_m3, 0);
    @memset(plants.root_carbon_g_m2, 0);
    @memset(thermal.thermal_conductivity_m_megajoules_per_h_k, 0);
    @memset(hydrology.micropore_face_flux_m3_per_step, 0);
    @memset(soil_gas.dissolved_mass_g, 0);
    @memset(litter_gas.gaseous_mass_g, 0);
    @memset(dissolved_gas.boundary_net_flux_g, 0);
    @memset(micropore_solutes.amount_mol, 0);
    @memset(macropore_solutes.amount_mol, 0);
    try loadTile(
        store,
        plan,
        plan.tiles[0],
        &grid,
        &plants,
        &thermal,
        &hydrology,
        &soil_gas,
        &litter_gas,
        &dissolved_gas,
        &micropore_solutes,
        &macropore_solutes,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_water,
        grid.matrix_liquid_water_m3,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_roots,
        plants.root_carbon_g_m2,
    );
    try std.testing.expectEqualSlices(
        usize,
        expected_active_layers,
        grid.active_soil_layer_count,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_thermal_conductivity,
        thermal.thermal_conductivity_m_megajoules_per_h_k,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_micropore_face_flux,
        hydrology.micropore_face_flux_m3_per_step,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_soil_dissolved_gas,
        soil_gas.dissolved_mass_g,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_litter_gaseous_mass,
        litter_gas.gaseous_mass_g,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_dissolved_boundary_flux,
        dissolved_gas.boundary_net_flux_g,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_micropore_solute_amount,
        micropore_solutes.amount_mol,
    );
    try std.testing.expectEqualSlices(
        f64,
        expected_macropore_solute_amount,
        macropore_solutes.amount_mol,
    );
}

const IncrementContext = struct {
    surface_temperature_k: []f64,
};

fn incrementSurfaceTemperature(
    context: *IncrementContext,
    cells: []const usize,
) !void {
    for (cells) |cell| context.surface_temperature_k[cell] += 1;
}

test "serial Morton transactions read an immutable source generation" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 6,
            .lat_count = 4,
            .soil_layers = 2,
            .plant_populations = 3,
        },
        .{ .worker_threads = 3, .tile_cells = 24 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var plants = try grid_module.PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var thermal = try initTestThermal(std.testing.allocator, grid);
    defer thermal.deinit();
    var hydrology = try transport_hydrology.State.init(
        std.testing.allocator,
        config.lon_count,
        config.lat_count,
        config.soil_layers,
        2,
    );
    defer hydrology.deinit();
    var soil_gas = try gas_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer soil_gas.deinit();
    var litter_gas = try gas_transport.State.init(
        std.testing.allocator,
        grid.cell_count,
    );
    defer litter_gas.deinit();
    @memset(soil_gas.temperature_k, 273.15);
    @memset(litter_gas.temperature_k, 273.15);
    var dissolved_gas = try dissolved_gas_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer dissolved_gas.deinit();
    var micropore_solutes = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        9,
    );
    defer micropore_solutes.deinit();
    var macropore_solutes = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        9,
    );
    defer macropore_solutes.deinit();
    for (grid.surface_temperature_k, 0..) |*value, cell|
        value.* = @floatFromInt(cell);
    var source_temporary = std.testing.tmpDir(.{});
    defer source_temporary.cleanup();
    var destination_temporary = std.testing.tmpDir(.{});
    defer destination_temporary.cleanup();
    const source_store = try tile_state_store.FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        source_temporary.dir,
        4096,
        17,
    );
    const destination_store = try tile_state_store.FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        destination_temporary.dir,
        4096,
        18,
    );
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        4,
        6,
        2,
        2,
        2,
    );
    defer plan.deinit();
    try initializeStore(
        source_store,
        plan,
        &grid,
        &plants,
        &thermal,
        &hydrology,
        &soil_gas,
        &litter_gas,
        &dissolved_gas,
        &micropore_solutes,
        &macropore_solutes,
    );
    var increment_context = IncrementContext{
        .surface_temperature_k = grid.surface_temperature_k,
    };
    const executor = try compute.CpuExecutor.init(
        std.testing.allocator,
        3,
        24,
    );
    try runSerialTiles(
        executor,
        source_store,
        destination_store,
        plan,
        &grid,
        &plants,
        &thermal,
        &hydrology,
        &soil_gas,
        &litter_gas,
        &dissolved_gas,
        &micropore_solutes,
        &macropore_solutes,
        &increment_context,
        incrementSurfaceTemperature,
    );
    @memset(grid.surface_temperature_k, 0);
    for (plan.tiles, 0..) |tile, tile_index| {
        try loadTile(
            destination_store,
            plan,
            tile,
            &grid,
            &plants,
            &thermal,
            &hydrology,
            &soil_gas,
            &litter_gas,
            &dissolved_gas,
            &micropore_solutes,
            &macropore_solutes,
        );
        const owned_cells = try plan.ownedCells(tile_index);
        for (owned_cells) |cell| {
            try std.testing.expectEqual(
                @as(f64, @floatFromInt(cell)) + 1,
                grid.surface_temperature_k[cell],
            );
        }
    }
}

test "serial Morton transactions reject one mutable generation" {
    const config = try @import("config.zig").SimulationConfig.init(
        .{
            .lon_count = 1,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 2,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var plants = try grid_module.PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var thermal = try initTestThermal(std.testing.allocator, grid);
    defer thermal.deinit();
    var hydrology = try transport_hydrology.State.init(
        std.testing.allocator,
        config.lon_count,
        config.lat_count,
        config.soil_layers,
        1,
    );
    defer hydrology.deinit();
    var soil_gas = try gas_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer soil_gas.deinit();
    var litter_gas = try gas_transport.State.init(
        std.testing.allocator,
        grid.cell_count,
    );
    defer litter_gas.deinit();
    var dissolved_gas = try dissolved_gas_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer dissolved_gas.deinit();
    var micropore_solutes = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        6,
    );
    defer micropore_solutes.deinit();
    var macropore_solutes = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        6,
    );
    defer macropore_solutes.deinit();
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const store = try tile_state_store.FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        4096,
        4,
    );
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        1,
        1,
        1,
        1,
        2,
    );
    defer plan.deinit();
    var increment_context = IncrementContext{
        .surface_temperature_k = grid.surface_temperature_k,
    };
    const executor = try compute.CpuExecutor.init(
        std.testing.allocator,
        1,
        1,
    );
    try std.testing.expectError(
        error.TileGenerationMustBeDoubleBuffered,
        runSerialTiles(
            executor,
            store,
            store,
            plan,
            &grid,
            &plants,
            &thermal,
            &hydrology,
            &soil_gas,
            &litter_gas,
            &dissolved_gas,
            &micropore_solutes,
            &macropore_solutes,
            &increment_context,
            incrementSurfaceTemperature,
        ),
    );
}
