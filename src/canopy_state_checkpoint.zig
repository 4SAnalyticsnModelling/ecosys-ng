const std = @import("std");
const State = @import("canopy_photosynthesis.zig").State;
const Retention = @import("canopy_precipitation_retention.zig").State;
const LayerDistribution = @import("canopy_layer_distribution.zig").State;
const magic = "ECOSCANP";
const version: u32 = 4;
pub const Limits = struct {
    maximum_cells: usize,
    maximum_species: usize,
    maximum_branches: usize,
    maximum_nodes: usize,
    maximum_samples: usize,
    maximum_layers: usize,
    maximum_inclinations: usize,
    maximum_azimuths: usize,
};
pub const View = struct {
    canopy: *const State,
    retention: *const Retention,
    layer_distribution: *const LayerDistribution,
};
pub const Owned = struct {
    canopy: State,
    retention: Retention,
    layer_distribution: LayerDistribution,
    pub fn deinit(self: *Owned) void {
        self.layer_distribution.deinit();
        self.retention.deinit();
        self.canopy.deinit();
        self.* = undefined;
    }
};

pub fn write(writer: anytype, view: View) !void {
    const state = view.canopy.*;
    try validate(view);
    try writer.writeAll(magic);
    try writer.writeInt(u32, version, .little);
    try writer.writeInt(u64, @intCast(state.cell_count), .little);
    try writer.writeInt(u64, @intCast(state.species_count), .little);
    try writer.writeInt(u64, @intCast(view.layer_distribution.layer_count), .little);
    try writer.writeInt(u64, @intCast(view.layer_distribution.inclination_count), .little);
    try writer.writeInt(u64, @intCast(view.layer_distribution.azimuth_count), .little);
    try writeOffsets(writer, state.plant_branch_offsets);
    try writeOffsets(writer, state.branch_node_offsets);
    try writeOffsets(writer, state.node_sample_offsets);
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) try writeF64Slice(writer, @field(state, field.name));
    inline for (@typeInfo(Retention).@"struct".fields) |field| if (field.type == []f64) try writeF64Slice(writer, @field(view.retention, field.name));
    inline for (@typeInfo(LayerDistribution).@"struct".fields) |field| if (field.type == []f64) try writeF64Slice(writer, @field(view.layer_distribution, field.name));
}

pub fn read(allocator: std.mem.Allocator, reader: *std.Io.Reader, limits: Limits) !Owned {
    if (limits.maximum_cells == 0 or limits.maximum_species == 0 or limits.maximum_branches == 0 or limits.maximum_nodes == 0 or limits.maximum_samples == 0 or limits.maximum_layers == 0 or limits.maximum_inclinations == 0 or limits.maximum_azimuths == 0) return error.InvalidCanopyCheckpointLimits;
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic)) return error.InvalidCanopyCheckpointMagic;
    if (try reader.takeInt(u32, .little) != version) return error.UnsupportedCanopyCheckpointVersion;
    const cells = try bounded(reader, limits.maximum_cells, error.CanopyCheckpointCellLimitExceeded);
    const species = try bounded(reader, limits.maximum_species, error.CanopyCheckpointSpeciesLimitExceeded);
    const layers = try bounded(reader, limits.maximum_layers, error.CanopyCheckpointLayerLimitExceeded);
    const inclinations = try bounded(reader, limits.maximum_inclinations, error.CanopyCheckpointInclinationLimitExceeded);
    const azimuths = try bounded(reader, limits.maximum_azimuths, error.CanopyCheckpointAzimuthLimitExceeded);
    if (layers == 0 or inclinations == 0 or azimuths == 0)
        return error.InvalidCanopyCheckpointDimensions;
    const plants = try std.math.mul(usize, cells, species);
    const branch_counts = try readCounts(allocator, reader, plants, limits.maximum_branches);
    defer allocator.free(branch_counts.counts);
    const node_counts = try readCounts(allocator, reader, branch_counts.total, limits.maximum_nodes);
    defer allocator.free(node_counts.counts);
    const sample_counts = try readCounts(allocator, reader, node_counts.total, limits.maximum_samples);
    defer allocator.free(sample_counts.counts);
    var state = try State.init(allocator, cells, species, branch_counts.counts, node_counts.counts, sample_counts.counts);
    errdefer state.deinit();
    var retention = try Retention.init(allocator, cells, species);
    errdefer retention.deinit();
    var layer_distribution = try LayerDistribution.init(
        allocator,
        cells,
        species,
        layers,
        inclinations,
        azimuths,
        &state,
    );
    errdefer layer_distribution.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) try readF64Slice(reader, @field(state, field.name));
    inline for (@typeInfo(Retention).@"struct".fields) |field| if (field.type == []f64) try readF64Slice(reader, @field(retention, field.name));
    inline for (@typeInfo(LayerDistribution).@"struct".fields) |field| if (field.type == []f64) try readF64Slice(reader, @field(layer_distribution, field.name));
    if (reader.peekByte()) |_| return error.TrailingCanopyCheckpointData else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
    const result = Owned{ .canopy = state, .retention = retention, .layer_distribution = layer_distribution };
    try validate(.{ .canopy = &result.canopy, .retention = &result.retention, .layer_distribution = &result.layer_distribution });
    return result;
}

const Counts = struct { counts: []usize, total: usize };
fn readCounts(allocator: std.mem.Allocator, reader: *std.Io.Reader, parent_count: usize, maximum_total: usize) !Counts {
    const stored_parent = try bounded(reader, parent_count, error.InvalidCanopyCheckpointOffsets);
    if (stored_parent != parent_count) return error.InvalidCanopyCheckpointOffsets;
    const total = try bounded(reader, maximum_total, error.CanopyCheckpointTopologyLimitExceeded);
    if (parent_count == 0 or total == 0) return error.InvalidCanopyCheckpointOffsets;
    const counts = try allocator.alloc(usize, parent_count);
    errdefer allocator.free(counts);
    var previous: usize = 0;
    for (0..parent_count + 1) |index| {
        const offset = try readUsize(reader);
        if (index == 0) {
            if (offset != 0) return error.InvalidCanopyCheckpointOffsets;
        } else {
            if (offset < previous or offset > total) return error.InvalidCanopyCheckpointOffsets;
            counts[index - 1] = offset - previous;
        }
        previous = offset;
    }
    if (previous != total) return error.InvalidCanopyCheckpointOffsets;
    return .{ .counts = counts, .total = total };
}
fn writeOffsets(writer: anytype, offsets: []const usize) !void {
    if (offsets.len < 2) return error.InvalidCanopyCheckpointOffsets;
    try writer.writeInt(u64, @intCast(offsets.len - 1), .little);
    try writer.writeInt(u64, @intCast(offsets[offsets.len - 1]), .little);
    for (offsets) |value| try writer.writeInt(u64, @intCast(value), .little);
}
fn validate(view: View) !void {
    const state = view.canopy.*;
    const plants = try std.math.mul(usize, state.cell_count, state.species_count);
    if (state.plant_branch_offsets.len != plants + 1) return error.InvalidCanopyCheckpointOffsets;
    inline for (.{ state.plant_branch_offsets, state.branch_node_offsets, state.node_sample_offsets }) |offsets| {
        if (offsets.len < 2 or offsets[0] != 0) return error.InvalidCanopyCheckpointOffsets;
        for (0..offsets.len - 1) |i| if (offsets[i] > offsets[i + 1]) return error.InvalidCanopyCheckpointOffsets;
    }
    try state.validateFinite();
    if (view.retention.cell_count != state.cell_count or view.retention.species_count != state.species_count) return error.CanopyCheckpointRetentionDimensionMismatch;
    inline for (@typeInfo(Retention).@"struct".fields) |field| if (field.type == []f64) for (@field(view.retention, field.name)) |value| if (!std.math.isFinite(value) or value < -1e-14) return error.InvalidCanopyCheckpointRetention;
    const layers = view.layer_distribution;
    if (layers.cell_count != state.cell_count or
        layers.species_count != state.species_count or
        layers.node_count != state.node_sample_offsets.len - 1 or
        layers.branch_count != state.branch_node_offsets.len - 1 or
        layers.layer_count == 0 or layers.inclination_count == 0 or
        layers.azimuth_count == 0)
        return error.CanopyCheckpointLayerDimensionMismatch;
    inline for (@typeInfo(LayerDistribution).@"struct".fields) |field|
        if (field.type == []f64) for (@field(layers, field.name)) |value|
            if (!std.math.isFinite(value) or value < -1e-14)
                return error.InvalidCanopyCheckpointLayerState;
}
fn bounded(reader: *std.Io.Reader, limit: usize, comptime too_large: anyerror) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > limit or value > std.math.maxInt(usize)) return too_large;
    return @intCast(value);
}
fn readUsize(reader: *std.Io.Reader) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > std.math.maxInt(usize)) return error.CanopyCheckpointIntegerOverflow;
    return @intCast(value);
}
fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteCanopyCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    }
}
fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value| {
        value.* = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(value.*)) return error.NonFiniteCanopyCheckpoint;
    }
}

test "canopy checkpoint reconstructs arbitrary species branch node sample topology" {
    const branch_counts = [_]usize{ 1, 2, 1, 1, 3, 1, 1 };
    const node_counts = [_]usize{ 1, 2, 1, 1, 1, 2, 1, 1, 1, 1 };
    const sample_counts = [_]usize{ 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1 };
    var source = try State.init(std.testing.allocator, 1, 7, &branch_counts, &node_counts, &sample_counts);
    defer source.deinit();
    var retention = try Retention.init(std.testing.allocator, 1, 7);
    defer retention.deinit();
    var layers = try LayerDistribution.init(
        std.testing.allocator,
        1,
        7,
        6,
        3,
        4,
        &source,
    );
    defer layers.deinit();
    source.node_leaf_carbon_g[source.node_leaf_carbon_g.len - 1] = 12.5;
    source.sample_leaf_nitrogen_g[source.sample_leaf_nitrogen_g.len - 1] = 0.4;
    source.plant_seed_storage_phosphorus_g[6] = 0.2;
    retention.living_surface_water_m3[6] = 0.03;
    retention.previous_water_energy_mj[6] = 7.25;
    layers.boundary_height_m[layers.boundary_height_m.len - 1] = 2.75;
    layers.node_leaf_area_m2[layers.node_leaf_area_m2.len - 1] = 0.125;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, .{ .canopy = &source, .retention = &retention, .layer_distribution = &layers });
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_cells = 2, .maximum_species = 20, .maximum_branches = 100, .maximum_nodes = 200, .maximum_samples = 300, .maximum_layers = 20, .maximum_inclinations = 10, .maximum_azimuths = 10 });
    defer restored.deinit();
    try std.testing.expectEqualSlices(usize, source.plant_branch_offsets, restored.canopy.plant_branch_offsets);
    try std.testing.expectEqualSlices(usize, source.branch_node_offsets, restored.canopy.branch_node_offsets);
    try std.testing.expectEqualSlices(usize, source.node_sample_offsets, restored.canopy.node_sample_offsets);
    try std.testing.expectEqualSlices(f64, source.node_leaf_carbon_g, restored.canopy.node_leaf_carbon_g);
    try std.testing.expectEqualSlices(f64, source.sample_leaf_nitrogen_g, restored.canopy.sample_leaf_nitrogen_g);
    try std.testing.expectEqualSlices(f64, source.plant_seed_storage_phosphorus_g, restored.canopy.plant_seed_storage_phosphorus_g);
    try std.testing.expectEqualSlices(f64, retention.living_surface_water_m3, restored.retention.living_surface_water_m3);
    try std.testing.expectEqualSlices(f64, retention.previous_water_energy_mj, restored.retention.previous_water_energy_mj);
    try std.testing.expectEqualSlices(f64, layers.boundary_height_m, restored.layer_distribution.boundary_height_m);
    try std.testing.expectEqualSlices(f64, layers.node_leaf_area_m2, restored.layer_distribution.node_leaf_area_m2);
}
