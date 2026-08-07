const std = @import("std");
const geometry = @import("litter_geometry.zig");
const organic = @import("../soil/organic/initialization.zig");
const compute = @import("../core/compute.zig");

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    water_retention_capacity_m3: []f64,
    dry_litter_volume_m3: []f64,
    expanded_total_volume_m3: []f64,
    dry_mass_megagrams: []f64,
    pore_volume_m3: []f64,
    air_volume_m3: []f64,
    porosity_m3_per_m3: []f64,
    field_capacity_m3_per_m3: []f64,
    wilting_point_m3_per_m3: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroSurfaceLitterGeometryCells;
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, cell_count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

pub const ApplyContext = struct {
    result: *State,
    surface_organic: *const organic.State,
    water_m3: []const f64,
    ice_m3: []const f64,
    charcoal_carbon_g_c: []const f64,
    parameters: geometry.Parameters,
};

pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    const cells = context.result.cell_count;
    if (range.first > range.end or range.end > cells or context.surface_organic.layer_count != cells or context.water_m3.len != cells or context.ice_m3.len != cells or context.charcoal_carbon_g_c.len != cells) return error.SurfaceLitterGeometryDimensionMismatch;
    for (range.first..range.end) |cell| {
        var carbon: [geometry.source_pool_count]f64 = undefined;
        for (&carbon, 0..) |*value, pool| value.* = try context.surface_organic.substrateCarbon_g_c(cell, pool);
        const value = try geometry.calculate(.{ .carbon_by_pool_g_c = carbon, .charcoal_carbon_g_c = context.charcoal_carbon_g_c[cell], .water_m3 = context.water_m3[cell], .ice_m3 = context.ice_m3[cell] }, context.parameters);
        inline for (@typeInfo(geometry.Result).@"struct".fields) |field| @field(context.result, field.name)[cell] = @field(value, field.name);
    }
}

test "surface litter geometry is runtime sized and tile independent" {
    var organic_state = try organic.State.init(std.testing.allocator, 3);
    defer organic_state.deinit();
    organic_state.dissolved[1 * organic.substrate_count].carbon_g_c = 10;
    organic_state.dissolved[2 * organic.substrate_count + 1].carbon_g_c = 20;
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    var context: ApplyContext = .{ .result = &state, .surface_organic = &organic_state, .water_m3 = &.{ 0, 0, 0 }, .ice_m3 = &.{ 0, 0, 0 }, .charcoal_carbon_g_c = &.{ 0, 0, 0 }, .parameters = .{ .water_retention_m3_per_g_c = .{ 2e-6, 5e-6, 5e-6, 5e-6, 5e-6 }, .dry_bulk_density_megagrams_per_m3 = .{ 0.1, 0.0125, 0.025, 0.025, 0.025 }, .dry_mass_megagrams_per_g_c = 1.82e-6, .particle_density_megagrams_per_m3 = 1.3, .field_capacity_fraction_of_porosity = 0.5, .wilting_point_fraction_of_porosity = 0.25 } };
    try applyTile(&context, .{ .first = 1, .end = 3 });
    try std.testing.expectEqual(@as(f64, 0), state.dry_litter_volume_m3[0]);
    try std.testing.expect(state.dry_litter_volume_m3[1] > 0);
    try std.testing.expect(state.dry_litter_volume_m3[2] > state.dry_litter_volume_m3[1]);
}
