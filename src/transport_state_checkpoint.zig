const std = @import("std");
const Solute = @import("solute_transport.zig").State;
const Gas = @import("gas_transport.zig").State;
const Snow = @import("snow_solute_transport.zig").State;
const Surface = @import("surface_solute_routing.zig").State;
const MineralNitrogen = @import("mineral_nitrogen_transport.zig").State;
const OrganicTransport = @import("soil_organic_transport.zig").State;
const magic = "ECOSTRNS";
// Version 8 persists the prognostic macropore organic-solute inventory.
const version: u32 = 8;
pub const View = struct { micropore: *const Solute, macropore: *const Solute, mineral_nitrogen: *const MineralNitrogen, organic: *const OrganicTransport, gas: *const Gas, litter_gas: *const Gas, snow: *const Snow, surface: *const Surface };
pub const Limits = struct { maximum_transport_cells: usize, maximum_solute_species: usize, maximum_snow_cells: usize, maximum_snow_layers: usize };
pub const Owned = struct {
    micropore: Solute,
    macropore: Solute,
    mineral_nitrogen: MineralNitrogen,
    organic: OrganicTransport,
    gas: Gas,
    litter_gas: Gas,
    snow: Snow,
    surface: Surface,
    pub fn deinit(self: *Owned) void {
        self.surface.deinit();
        self.snow.deinit();
        self.litter_gas.deinit();
        self.gas.deinit();
        self.organic.deinit();
        self.mineral_nitrogen.deinit();
        self.macropore.deinit();
        self.micropore.deinit();
        self.* = undefined;
    }
};

pub fn write(writer: anytype, view: View) !void {
    try validate(view);
    try writer.writeAll(magic);
    try writer.writeInt(u32, version, .little);
    try writer.writeInt(u64, @intCast(view.micropore.cell_count), .little);
    try writer.writeInt(u64, @intCast(view.micropore.species_count), .little);
    try writer.writeInt(u64, @intCast(view.snow.cell_count), .little);
    try writer.writeInt(u64, @intCast(view.snow.layer_capacity), .little);
    try writer.writeInt(u64, @intCast(view.surface.columns), .little);
    try writer.writeInt(u64, @intCast(view.surface.rows), .little);
    try writer.writeInt(u64, @intCast(view.surface.species_count), .little);
    inline for (.{ view.micropore.water_volume_m3, view.micropore.amount_mol, view.macropore.water_volume_m3, view.macropore.amount_mol, view.gas.air_volume_m3, view.gas.temperature_k, view.gas.water_vapor_mol, view.gas.gaseous_mass_g, view.gas.dissolved_mass_g, view.gas.macropore_dissolved_mass_g, view.gas.band_dissolved_mass_g }) |values| try writeF64Slice(writer, values);
    inline for (.{ view.mineral_nitrogen.matrix.water_volume_m3, view.mineral_nitrogen.matrix.amount_mol, view.mineral_nitrogen.macropore.water_volume_m3, view.mineral_nitrogen.macropore.amount_mol, view.mineral_nitrogen.boundary_export_g_n_per_step }) |values| try writeF64Slice(writer, values);
    inline for (.{ view.organic.micropore_amount_g, view.organic.macropore_amount_g, view.organic.boundary_net_flux_g }) |values| try writeF64Slice(writer, values);
    inline for (.{ view.litter_gas.air_volume_m3, view.litter_gas.temperature_k, view.litter_gas.water_vapor_mol, view.litter_gas.gaseous_mass_g, view.litter_gas.dissolved_mass_g, view.litter_gas.macropore_dissolved_mass_g, view.litter_gas.band_dissolved_mass_g }) |values| try writeF64Slice(writer, values);
    try writeBoolSlice(writer, view.snow.active);
    inline for (.{ view.snow.solid_snow_water_equivalent_m3, view.snow.liquid_water_volume_m3, view.snow.vapor_water_equivalent_m3, view.snow.ice_volume_m3, view.snow.air_filled_volume_m3, view.snow.total_layer_volume_m3, view.snow.target_layer_volume_m3, view.snow.layer_thickness_m, view.snow.cumulative_depth_m, view.snow.snow_density_megagrams_per_m3, view.snow.temperature_k, view.snow.heat_capacity_megajoules_per_k, view.snow.horizontal_area_m2 }) |values| try writeF64Slice(writer, values);
    try writeF64Slice(writer, view.snow.amount_g);
    try writeF64Slice(writer, view.surface.carrier_volume_m3);
    try writeF64Slice(writer, view.surface.amount_mol);
}

pub fn read(allocator: std.mem.Allocator, reader: *std.Io.Reader, limits: Limits) !Owned {
    if (limits.maximum_transport_cells == 0 or limits.maximum_solute_species == 0 or limits.maximum_snow_cells == 0 or limits.maximum_snow_layers == 0) return error.InvalidTransportCheckpointLimits;
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic)) return error.InvalidTransportCheckpointMagic;
    if (try reader.takeInt(u32, .little) != version) return error.UnsupportedTransportCheckpointVersion;
    const cells = try bounded(reader, limits.maximum_transport_cells, error.TransportCheckpointCellLimitExceeded);
    const species = try bounded(reader, limits.maximum_solute_species, error.TransportCheckpointSpeciesLimitExceeded);
    const snow_cells = try bounded(reader, limits.maximum_snow_cells, error.TransportCheckpointSnowCellLimitExceeded);
    const snow_layers = try bounded(reader, limits.maximum_snow_layers, error.TransportCheckpointSnowLayerLimitExceeded);
    const surface_columns = try bounded(reader, limits.maximum_transport_cells, error.TransportCheckpointCellLimitExceeded);
    const surface_rows = try bounded(reader, limits.maximum_transport_cells, error.TransportCheckpointCellLimitExceeded);
    const surface_species = try bounded(reader, limits.maximum_solute_species, error.TransportCheckpointSpeciesLimitExceeded);
    if (cells == 0 or species == 0 or snow_cells == 0 or snow_layers == 0) return error.InvalidTransportCheckpointDimensions;
    var micropore = try Solute.init(allocator, cells, species);
    errdefer micropore.deinit();
    var macropore = try Solute.init(allocator, cells, species);
    errdefer macropore.deinit();
    var mineral_nitrogen = try MineralNitrogen.init(allocator, cells);
    errdefer mineral_nitrogen.deinit();
    var organic_transport = try OrganicTransport.init(allocator, cells);
    errdefer organic_transport.deinit();
    var gas = try Gas.init(allocator, cells);
    errdefer gas.deinit();
    var litter_gas = try Gas.init(allocator, snow_cells);
    errdefer litter_gas.deinit();
    var snow = try Snow.init(allocator, snow_cells, snow_layers);
    errdefer snow.deinit();
    var surface = try Surface.init(allocator, surface_columns, surface_rows, surface_species);
    errdefer surface.deinit();
    inline for (.{ micropore.water_volume_m3, micropore.amount_mol, macropore.water_volume_m3, macropore.amount_mol, gas.air_volume_m3, gas.temperature_k, gas.water_vapor_mol, gas.gaseous_mass_g, gas.dissolved_mass_g, gas.macropore_dissolved_mass_g, gas.band_dissolved_mass_g }) |values| try readF64Slice(reader, values);
    inline for (.{ mineral_nitrogen.matrix.water_volume_m3, mineral_nitrogen.matrix.amount_mol, mineral_nitrogen.macropore.water_volume_m3, mineral_nitrogen.macropore.amount_mol, mineral_nitrogen.boundary_export_g_n_per_step }) |values| try readF64Slice(reader, values);
    inline for (.{ organic_transport.micropore_amount_g, organic_transport.macropore_amount_g, organic_transport.boundary_net_flux_g }) |values| try readF64Slice(reader, values);
    inline for (.{ litter_gas.air_volume_m3, litter_gas.temperature_k, litter_gas.water_vapor_mol, litter_gas.gaseous_mass_g, litter_gas.dissolved_mass_g, litter_gas.macropore_dissolved_mass_g, litter_gas.band_dissolved_mass_g }) |values| try readF64Slice(reader, values);
    try readBoolSlice(reader, snow.active);
    inline for (.{ snow.solid_snow_water_equivalent_m3, snow.liquid_water_volume_m3, snow.vapor_water_equivalent_m3, snow.ice_volume_m3, snow.air_filled_volume_m3, snow.total_layer_volume_m3, snow.target_layer_volume_m3, snow.layer_thickness_m, snow.cumulative_depth_m, snow.snow_density_megagrams_per_m3, snow.temperature_k, snow.heat_capacity_megajoules_per_k, snow.horizontal_area_m2 }) |values| try readF64Slice(reader, values);
    try readF64Slice(reader, snow.amount_g);
    try readF64Slice(reader, surface.carrier_volume_m3);
    try readF64Slice(reader, surface.amount_mol);
    if (reader.peekByte()) |_| return error.TrailingTransportCheckpointData else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
    const result = Owned{ .micropore = micropore, .macropore = macropore, .mineral_nitrogen = mineral_nitrogen, .organic = organic_transport, .gas = gas, .litter_gas = litter_gas, .snow = snow, .surface = surface };
    try validate(.{ .micropore = &result.micropore, .macropore = &result.macropore, .mineral_nitrogen = &result.mineral_nitrogen, .organic = &result.organic, .gas = &result.gas, .litter_gas = &result.litter_gas, .snow = &result.snow, .surface = &result.surface });
    return result;
}

fn validate(view: View) !void {
    if (view.micropore.cell_count == 0 or view.micropore.species_count == 0 or view.macropore.cell_count != view.micropore.cell_count or view.macropore.species_count != view.micropore.species_count or view.mineral_nitrogen.cell_count != view.micropore.cell_count or view.organic.layer_count != view.micropore.cell_count or view.gas.cell_count != view.micropore.cell_count or view.snow.cell_count == 0 or view.litter_gas.cell_count != view.snow.cell_count or view.snow.layer_capacity == 0 or view.surface.columns == 0 or view.surface.rows == 0 or view.surface.species_count != view.micropore.species_count or try std.math.mul(usize, view.surface.columns, view.surface.rows) != view.snow.cell_count) return error.InvalidTransportCheckpointDimensions;
    inline for (.{ view.micropore.water_volume_m3, view.micropore.amount_mol, view.macropore.water_volume_m3, view.macropore.amount_mol, view.gas.air_volume_m3, view.gas.water_vapor_mol, view.gas.gaseous_mass_g, view.gas.dissolved_mass_g, view.gas.macropore_dissolved_mass_g, view.gas.band_dissolved_mass_g, view.snow.solid_snow_water_equivalent_m3, view.snow.liquid_water_volume_m3, view.snow.vapor_water_equivalent_m3, view.snow.ice_volume_m3, view.snow.air_filled_volume_m3, view.snow.total_layer_volume_m3, view.snow.target_layer_volume_m3, view.snow.layer_thickness_m, view.snow.cumulative_depth_m, view.snow.snow_density_megagrams_per_m3, view.snow.temperature_k, view.snow.heat_capacity_megajoules_per_k, view.snow.horizontal_area_m2, view.snow.amount_g, view.surface.carrier_volume_m3, view.surface.amount_mol }) |values| try validateNonnegative(values);
    inline for (.{ view.litter_gas.air_volume_m3, view.litter_gas.water_vapor_mol, view.litter_gas.gaseous_mass_g, view.litter_gas.dissolved_mass_g, view.litter_gas.macropore_dissolved_mass_g, view.litter_gas.band_dissolved_mass_g }) |values| try validateNonnegative(values);
    inline for (.{ view.mineral_nitrogen.matrix.water_volume_m3, view.mineral_nitrogen.matrix.amount_mol, view.mineral_nitrogen.macropore.water_volume_m3, view.mineral_nitrogen.macropore.amount_mol, view.mineral_nitrogen.boundary_export_g_n_per_step }) |values| try validateNonnegative(values);
    inline for (.{ view.organic.micropore_amount_g, view.organic.macropore_amount_g }) |values| try validateNonnegative(values);
    for (view.organic.boundary_net_flux_g) |value| if (!std.math.isFinite(value)) return error.InvalidTransportCheckpointInventory;
    for (view.gas.temperature_k) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidTransportCheckpointTemperature;
    for (view.litter_gas.temperature_k) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidTransportCheckpointTemperature;
}
fn validateNonnegative(values: []const f64) !void {
    for (values) |value| if (!std.math.isFinite(value) or value < -1e-14) return error.InvalidTransportCheckpointInventory;
}
fn bounded(reader: *std.Io.Reader, limit: usize, comptime too_large: anyerror) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > limit or value > std.math.maxInt(usize)) return too_large;
    return @intCast(value);
}
fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteTransportCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    }
}
fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value| {
        value.* = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(value.*)) return error.NonFiniteTransportCheckpoint;
    }
}
fn writeBoolSlice(writer: anytype, values: []const bool) !void {
    for (values) |value| try writer.writeByte(@intFromBool(value));
}
fn readBoolSlice(reader: *std.Io.Reader, values: []bool) !void {
    for (values) |*value| value.* = switch (try reader.takeByte()) {
        0 => false,
        1 => true,
        else => return error.InvalidTransportCheckpointBoolean,
    };
}

test "transport checkpoint round trips pore gas and runtime snow inventories" {
    var micro = try Solute.init(std.testing.allocator, 6, 50);
    defer micro.deinit();
    var macro = try Solute.init(std.testing.allocator, 6, 50);
    defer macro.deinit();
    var mineral_nitrogen = try MineralNitrogen.init(std.testing.allocator, 6);
    defer mineral_nitrogen.deinit();
    var organic_transport = try OrganicTransport.init(std.testing.allocator, 6);
    defer organic_transport.deinit();
    var gas = try Gas.init(std.testing.allocator, 6);
    defer gas.deinit();
    var litter_gas = try Gas.init(std.testing.allocator, 2);
    defer litter_gas.deinit();
    var snow = try Snow.init(std.testing.allocator, 2, 4);
    defer snow.deinit();
    var surface = try Surface.init(std.testing.allocator, 2, 1, 50);
    defer surface.deinit();
    micro.water_volume_m3[5] = 2;
    micro.amount_mol[micro.amount_mol.len - 1] = 3;
    macro.amount_mol[macro.amount_mol.len - 1] = 4;
    mineral_nitrogen.matrix.water_volume_m3[5] = 0.25;
    mineral_nitrogen.macropore.amount_mol[mineral_nitrogen.macropore.amount_mol.len - 1] = 4.5;
    mineral_nitrogen.boundary_export_g_n_per_step[5] = 0.125;
    organic_transport.macropore_amount_g[organic_transport.macropore_amount_g.len - 1] = 4.75;
    organic_transport.boundary_net_flux_g[organic_transport.boundary_net_flux_g.len - 1] = -0.375;
    gas.temperature_k[5] = 280;
    gas.gaseous_mass_g[gas.gaseous_mass_g.len - 1] = 5;
    gas.macropore_dissolved_mass_g[gas.macropore_dissolved_mass_g.len - 1] = 5.25;
    litter_gas.dissolved_mass_g[litter_gas.dissolved_mass_g.len - 1] = 5.5;
    snow.active[7] = true;
    snow.liquid_water_volume_m3[7] = 0.1;
    snow.solid_snow_water_equivalent_m3[7] = 0.2;
    snow.temperature_k[7] = 268;
    snow.amount_g[snow.amount_g.len - 1] = 6;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    surface.amount_mol[surface.amount_mol.len - 1] = 7;
    try write(&bytes.writer, .{ .micropore = &micro, .macropore = &macro, .mineral_nitrogen = &mineral_nitrogen, .organic = &organic_transport, .gas = &gas, .litter_gas = &litter_gas, .snow = &snow, .surface = &surface });
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_transport_cells = 20, .maximum_solute_species = 100, .maximum_snow_cells = 10, .maximum_snow_layers = 20 });
    defer restored.deinit();
    try std.testing.expectEqualSlices(f64, micro.amount_mol, restored.micropore.amount_mol);
    try std.testing.expectEqualSlices(f64, macro.amount_mol, restored.macropore.amount_mol);
    try std.testing.expectEqualSlices(f64, mineral_nitrogen.matrix.water_volume_m3, restored.mineral_nitrogen.matrix.water_volume_m3);
    try std.testing.expectEqualSlices(f64, mineral_nitrogen.macropore.amount_mol, restored.mineral_nitrogen.macropore.amount_mol);
    try std.testing.expectEqualSlices(f64, mineral_nitrogen.boundary_export_g_n_per_step, restored.mineral_nitrogen.boundary_export_g_n_per_step);
    try std.testing.expectEqualSlices(f64, organic_transport.micropore_amount_g, restored.organic.micropore_amount_g);
    try std.testing.expectEqualSlices(f64, organic_transport.macropore_amount_g, restored.organic.macropore_amount_g);
    try std.testing.expectEqualSlices(f64, organic_transport.boundary_net_flux_g, restored.organic.boundary_net_flux_g);
    try std.testing.expectEqualSlices(f64, gas.gaseous_mass_g, restored.gas.gaseous_mass_g);
    try std.testing.expectEqualSlices(f64, gas.macropore_dissolved_mass_g, restored.gas.macropore_dissolved_mass_g);
    try std.testing.expectEqualSlices(f64, litter_gas.dissolved_mass_g, restored.litter_gas.dissolved_mass_g);
    try std.testing.expectEqualSlices(bool, snow.active, restored.snow.active);
    try std.testing.expectEqualSlices(f64, snow.solid_snow_water_equivalent_m3, restored.snow.solid_snow_water_equivalent_m3);
    try std.testing.expectEqualSlices(f64, snow.temperature_k, restored.snow.temperature_k);
    try std.testing.expectEqualSlices(f64, snow.amount_g, restored.snow.amount_g);
    try std.testing.expectEqualSlices(f64, surface.amount_mol, restored.surface.amount_mol);
}

test "transport checkpoint enforces snow layer limit before allocation" {
    var micro = try Solute.init(std.testing.allocator, 1, 1);
    defer micro.deinit();
    var macro = try Solute.init(std.testing.allocator, 1, 1);
    defer macro.deinit();
    var mineral_nitrogen = try MineralNitrogen.init(std.testing.allocator, 1);
    defer mineral_nitrogen.deinit();
    var organic_transport = try OrganicTransport.init(std.testing.allocator, 1);
    defer organic_transport.deinit();
    var gas = try Gas.init(std.testing.allocator, 1);
    defer gas.deinit();
    var litter_gas = try Gas.init(std.testing.allocator, 1);
    defer litter_gas.deinit();
    var snow = try Snow.init(std.testing.allocator, 1, 4);
    defer snow.deinit();
    var surface = try Surface.init(std.testing.allocator, 1, 1, 1);
    defer surface.deinit();
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, .{ .micropore = &micro, .macropore = &macro, .mineral_nitrogen = &mineral_nitrogen, .organic = &organic_transport, .gas = &gas, .litter_gas = &litter_gas, .snow = &snow, .surface = &surface });
    var reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(error.TransportCheckpointSnowLayerLimitExceeded, read(std.testing.allocator, &reader, .{ .maximum_transport_cells = 1, .maximum_solute_species = 1, .maximum_snow_cells = 1, .maximum_snow_layers = 3 }));
}

test "transport checkpoint rejects trailing corruption" {
    var micro = try Solute.init(std.testing.allocator, 1, 1);
    defer micro.deinit();
    var macro = try Solute.init(std.testing.allocator, 1, 1);
    defer macro.deinit();
    var mineral_nitrogen = try MineralNitrogen.init(std.testing.allocator, 1);
    defer mineral_nitrogen.deinit();
    var organic_transport = try OrganicTransport.init(std.testing.allocator, 1);
    defer organic_transport.deinit();
    var gas = try Gas.init(std.testing.allocator, 1);
    defer gas.deinit();
    var litter_gas = try Gas.init(std.testing.allocator, 1);
    defer litter_gas.deinit();
    var snow = try Snow.init(std.testing.allocator, 1, 1);
    defer snow.deinit();
    var surface = try Surface.init(std.testing.allocator, 1, 1, 1);
    defer surface.deinit();
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, .{ .micropore = &micro, .macropore = &macro, .mineral_nitrogen = &mineral_nitrogen, .organic = &organic_transport, .gas = &gas, .litter_gas = &litter_gas, .snow = &snow, .surface = &surface });
    try bytes.writer.writeByte(0xff);
    var reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(error.TrailingTransportCheckpointData, read(std.testing.allocator, &reader, .{ .maximum_transport_cells = 1, .maximum_solute_species = 1, .maximum_snow_cells = 1, .maximum_snow_layers = 1 }));
}
