const std = @import("std");
const CellRange = @import("compute.zig").CellRange;
const Topography = @import("topography.zig").Topography;
const SoilCatalog = @import("soil_catalog.zig").Catalog;
const RadiationState = @import("canopy_radiation.zig").State;
const InterceptionState = @import("canopy_interception.zig").State;
const TerrainState = @import("terrain_radiation.zig").State;

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    soil_albedo: []f64,
    initial_snow_depth_m: []f64,
    surface_albedo: []f64,
    incident_shortwave_megajoules_per_m2: []f64,
    absorbed_shortwave_megajoules_per_m2: []f64,
    reflected_shortwave_megajoules_per_m2: []f64,
    incident_par_micromol_per_m2_per_s: []f64,
    absorbed_par_micromol_per_m2_per_s: []f64,
    reflected_par_micromol_per_m2_per_s: []f64,

    pub fn initMapped(allocator: std.mem.Allocator, topography: Topography, topography_unit_by_cell: []const usize, soil_catalog: SoilCatalog, soil_catalog_index_by_cell: []const usize) !State {
        if (topography_unit_by_cell.len == 0 or topography_unit_by_cell.len != soil_catalog_index_by_cell.len) return error.InvalidGroundRadiationDimensions;
        const cell_count = topography_unit_by_cell.len;
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        var allocated: usize = 0;
        errdefer freeAllocated(&result, allocated);
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, cell_count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        for (topography_unit_by_cell, soil_catalog_index_by_cell, 0..) |topography_index, soil_index, cell| {
            if (topography_index >= topography.units.len or soil_index >= soil_catalog.entries.items.len) return error.GroundRadiationMappingOutOfBounds;
            const albedo = soil_catalog.entries.items[soil_index].profile.wet_soil_albedo;
            const snow_depth = topography.units[topography_index].initial_snowpack_depth_m;
            if (!std.math.isFinite(albedo) or albedo < 0 or albedo > 1 or !std.math.isFinite(snow_depth) or snow_depth < 0) return error.InvalidGroundSurfaceProperties;
            result.soil_albedo[cell] = albedo;
            result.initial_snow_depth_m[cell] = snow_depth;
            result.surface_albedo[cell] = initialSurfaceAlbedo(albedo, snow_depth);
        }
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| {
            if (!std.math.isFinite(value)) {
                std.log.err("non-finite ground radiation: field={s} index={d} value={e}", .{ field.name, index, value });
                return error.NonFiniteGroundRadiation;
            }
        };
    }
};

pub const ApplyContext = struct {
    result: *State,
    radiation: *const RadiationState,
    interception: ?*const InterceptionState,
    terrain: *const TerrainState,
};

pub fn applyTile(context: *ApplyContext, range: CellRange) !void {
    const result = context.result;
    if (range.end > result.cell_count or context.radiation.cellCount() != result.cell_count or context.terrain.cell_count != result.cell_count) return error.GroundRadiationDimensionMismatch;
    if (context.interception) |interception| if (interception.cell_count != result.cell_count) return error.GroundRadiationDimensionMismatch;
    for (range.first..range.end) |cell| {
        const direct_transmission = if (context.interception) |interception| interception.direct_transmission_fraction[cell] else 1.0;
        const diffuse_transmission = if (context.interception) |interception| interception.diffuse_transmission_fraction[cell] else 1.0;
        var diffuse_terrain_projection: f64 = 0;
        for (0..context.terrain.sky_sector_count) |sky| diffuse_terrain_projection += context.terrain.diffuse_sky_incidence_fraction[cell * context.terrain.sky_sector_count + sky];
        const direct_shortwave = context.radiation.direct_shortwave_megajoules_per_m2[cell] * direct_transmission * context.terrain.direct_solar_incidence_fraction[cell];
        var diffuse_shortwave = context.radiation.diffuse_shortwave_megajoules_per_m2[cell] * diffuse_transmission * diffuse_terrain_projection;
        const direct_par = context.radiation.direct_par_micromol_per_m2_per_s[cell] * direct_transmission * context.terrain.direct_solar_incidence_fraction[cell];
        var diffuse_par = context.radiation.diffuse_par_micromol_per_m2_per_s[cell] * diffuse_transmission * diffuse_terrain_projection;
        if (context.interception) |interception| {
            const bottom_boundary = cell * (interception.layer_count + 1);
            diffuse_shortwave += interception.downward_scattered_shortwave_by_boundary_megajoules_per_m2[bottom_boundary];
            diffuse_par += interception.downward_scattered_par_by_boundary_micromol_per_m2_per_s[bottom_boundary];
        }
        const shortwave = try partitionEnergy(direct_shortwave + diffuse_shortwave, result.surface_albedo[cell]);
        const par = try partitionEnergy(direct_par + diffuse_par, result.surface_albedo[cell]);
        result.incident_shortwave_megajoules_per_m2[cell] = shortwave.incident;
        result.absorbed_shortwave_megajoules_per_m2[cell] = shortwave.absorbed;
        result.reflected_shortwave_megajoules_per_m2[cell] = shortwave.reflected;
        result.incident_par_micromol_per_m2_per_s[cell] = par.incident;
        result.absorbed_par_micromol_per_m2_per_s[cell] = par.absorbed;
        result.reflected_par_micromol_per_m2_per_s[cell] = par.reflected;
    }
}

const EnergyPartition = struct { incident: f64, absorbed: f64, reflected: f64 };

fn partitionEnergy(incident: f64, albedo: f64) !EnergyPartition {
    if (!std.math.isFinite(incident) or !std.math.isFinite(albedo)) return error.NonFiniteGroundEnergy;
    if (incident < 0 or albedo < 0 or albedo > 1) return error.InvalidGroundEnergy;
    const reflected = incident * albedo;
    const absorbed = incident - reflected;
    if (@abs((absorbed + reflected) - incident) > 1.0e-12 * @max(1.0, incident)) return error.GroundEnergyImbalance;
    return .{ .incident = incident, .absorbed = absorbed, .reflected = reflected };
}

fn initialSurfaceAlbedo(soil_albedo: f64, snow_depth_m: f64) f64 {
    // Initially prescribed snow is solid snow in the available translated
    // state, matching ALBW=0.90. Later liquid/ice evolution updates this value.
    const snow_cover_fraction = @min(std.math.pow(f64, snow_depth_m / 0.07, 2), 1.0);
    return snow_cover_fraction * 0.90 + (1.0 - snow_cover_fraction) * soil_albedo;
}

fn freeAllocated(state: *State, count: usize) void {
    var visited: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (visited < count) state.allocator.free(@field(state, field.name));
        visited += 1;
    };
}

test "ground energy partition is conservative" {
    const result = try partitionEnergy(2.5, 0.2);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.absorbed, 1.0e-15);
    try std.testing.expectApproxEqAbs(result.incident, result.absorbed + result.reflected, 1.0e-15);
}

test "initial snow cover follows squared depth fraction" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), initialSurfaceAlbedo(0.2, 0), 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), initialSurfaceAlbedo(0.2, 0.07), 1.0e-15);
    try std.testing.expect(initialSurfaceAlbedo(0.2, 0.035) > 0.2);
}
