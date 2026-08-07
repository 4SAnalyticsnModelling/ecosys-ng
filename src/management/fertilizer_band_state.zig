const std = @import("std");
const geometry_module = @import("hourly_fertilizer_band_geometry.zig");
const phase_module = @import("fertilizer_band_phase_coordinator.zig");
const plant_available_nutrients = @import("../soil/nutrients/plant_available_nutrients.zig");
const charge_classification = @import("../soil/solute/charge_classification.zig");

pub const family_count = 3;
const checkpoint_magic = "ECOBAND1";
const checkpoint_version: u32 = 1;

pub const Family = phase_module.Family;
const LayerRange = struct { start: usize, end: usize };

pub const InitializationPolicy = enum {
    /// Reproduces the current ecosys-ng runscript fraction in every active
    /// layer without inventing a localized fertilizer application.
    uniform_profile_fraction,
};

pub const Initialization = struct {
    cell_count: usize,
    layer_capacity: usize,
    active_layer_count_by_cell: []const usize,
    layer_upper_depth_m: []const f64,
    layer_lower_depth_m: []const f64,
    layer_thickness_m: []const f64,
    initial_band_fraction_by_family: [family_count]f64,
    row_spacing_m_by_cell_family: []const f64,
    policy: InitializationPolicy = .uniform_profile_fraction,
};

/// Authoritative runtime fertilizer-band geometry and phase owner.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    layer_capacity: usize,
    active: []bool,
    row_spacing_m: []f64,
    upper_edge_depth_m: []f64,
    lower_edge_depth_m: []f64,
    band_depth_m: []f64,
    band_width_m: []f64,
    band_volume_fraction: []f64,
    non_band_volume_fraction: []f64,
    coordinators: []phase_module.Coordinator,

    pub fn init(
        allocator: std.mem.Allocator,
        initialization: Initialization,
    ) !State {
        try validateInitialization(initialization);
        const scalar_count = try std.math.mul(
            usize,
            initialization.cell_count,
            family_count,
        );
        const layer_count = try std.math.mul(
            usize,
            scalar_count,
            initialization.layer_capacity,
        );
        const active = try allocator.alloc(bool, scalar_count);
        errdefer allocator.free(active);
        const row = try allocator.alloc(f64, scalar_count);
        errdefer allocator.free(row);
        const upper = try allocator.alloc(f64, scalar_count);
        errdefer allocator.free(upper);
        const lower = try allocator.alloc(f64, scalar_count);
        errdefer allocator.free(lower);
        const depth = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(depth);
        const width = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(width);
        const band = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(band);
        const non_band = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(non_band);
        const coordinators = try allocator.alloc(
            phase_module.Coordinator,
            initialization.cell_count,
        );
        errdefer allocator.free(coordinators);
        var initialized_coordinators: usize = 0;
        errdefer for (coordinators[0..initialized_coordinators]) |*phase_state|
            phase_state.deinit();
        while (initialized_coordinators < coordinators.len) : (initialized_coordinators += 1) {
            coordinators[initialized_coordinators] =
                try phase_module.Coordinator.init(
                    allocator,
                    initialization.layer_capacity,
                );
        }

        var result: State = .{
            .allocator = allocator,
            .cell_count = initialization.cell_count,
            .layer_capacity = initialization.layer_capacity,
            .active = active,
            .row_spacing_m = row,
            .upper_edge_depth_m = upper,
            .lower_edge_depth_m = lower,
            .band_depth_m = depth,
            .band_width_m = width,
            .band_volume_fraction = band,
            .non_band_volume_fraction = non_band,
            .coordinators = coordinators,
        };
        result.initializeUniformProfile(initialization);
        return result;
    }

    /// Initializes the authoritative layer-aware geometry from the compulsory
    /// runscript plant-nutrient fractions and row spacings.
    pub fn initFromPlantNutrientParameters(
        allocator: std.mem.Allocator,
        cell_count: usize,
        layer_capacity: usize,
        active_layer_count_by_cell: []const usize,
        layer_upper_depth_m: []const f64,
        layer_lower_depth_m: []const f64,
        layer_thickness_m: []const f64,
        parameters: plant_available_nutrients.InitializationParameters,
    ) !State {
        try parameters.validate();
        const scalar_count = try std.math.mul(usize, cell_count, family_count);
        const row_spacing = try allocator.alloc(f64, scalar_count);
        defer allocator.free(row_spacing);
        const family_row_spacing_m = [family_count]f64{
            parameters.initial_ammonium_band_row_spacing_m,
            parameters.initial_nitrate_band_row_spacing_m,
            parameters.initial_phosphate_band_row_spacing_m,
        };
        for (0..cell_count) |cell| {
            inline for (family_row_spacing_m, 0..) |spacing_m, family_index| {
                row_spacing[cell * family_count + family_index] = spacing_m;
            }
        }
        return init(allocator, .{
            .cell_count = cell_count,
            .layer_capacity = layer_capacity,
            .active_layer_count_by_cell = active_layer_count_by_cell,
            .layer_upper_depth_m = layer_upper_depth_m,
            .layer_lower_depth_m = layer_lower_depth_m,
            .layer_thickness_m = layer_thickness_m,
            .initial_band_fraction_by_family = .{
                parameters.initial_ammonium_band_fraction,
                parameters.initial_nitrate_band_fraction,
                parameters.initial_phosphate_band_fraction,
            },
            .row_spacing_m_by_cell_family = row_spacing,
        });
    }

    pub fn deinit(self: *State) void {
        for (self.coordinators) |*phase_state| phase_state.deinit();
        self.allocator.free(self.coordinators);
        self.allocator.free(self.non_band_volume_fraction);
        self.allocator.free(self.band_volume_fraction);
        self.allocator.free(self.band_width_m);
        self.allocator.free(self.band_depth_m);
        self.allocator.free(self.lower_edge_depth_m);
        self.allocator.free(self.upper_edge_depth_m);
        self.allocator.free(self.row_spacing_m);
        self.allocator.free(self.active);
        self.* = undefined;
    }

    pub fn geometry(
        self: *State,
        cell: usize,
        family: Family,
    ) !geometry_module.State {
        const scalar = try self.scalarIndex(cell, family);
        const layers = try self.layerRange(cell, family);
        return .{
            .active = self.active[scalar],
            .row_spacing_m = self.row_spacing_m[scalar],
            .upper_edge_depth_m = self.upper_edge_depth_m[scalar],
            .lower_edge_depth_m = self.lower_edge_depth_m[scalar],
            .band_depth_m = self.band_depth_m[layers.start..layers.end],
            .band_width_m = self.band_width_m[layers.start..layers.end],
            .band_volume_fraction = self.band_volume_fraction[layers.start..layers.end],
            .non_band_volume_fraction = self.non_band_volume_fraction[layers.start..layers.end],
        };
    }

    /// Copies scalar fields changed through a geometry view back to its owner.
    pub fn commitGeometryScalars(
        self: *State,
        cell: usize,
        family: Family,
        view: geometry_module.State,
    ) !void {
        const scalar = try self.scalarIndex(cell, family);
        self.active[scalar] = view.active;
        self.row_spacing_m[scalar] = view.row_spacing_m;
        self.upper_edge_depth_m[scalar] = view.upper_edge_depth_m;
        self.lower_edge_depth_m[scalar] = view.lower_edge_depth_m;
    }

    pub fn coordinator(
        self: *State,
        cell: usize,
    ) !*phase_module.Coordinator {
        if (cell >= self.cell_count) return error.FertilizerBandCellOutOfBounds;
        return &self.coordinators[cell];
    }

    /// Returns the three band/non-band volume pairs for one runtime soil
    /// layer. `layer` is local to the cell, not a flattened storage index.
    pub fn zoneFractions(
        self: *const State,
        cell: usize,
        layer: usize,
    ) !charge_classification.ZoneFractions {
        if (cell >= self.cell_count)
            return error.FertilizerBandCellOutOfBounds;
        if (layer >= self.layer_capacity)
            return error.FertilizerBandLayerOutOfBounds;
        const ammonium = self.band_volume_fraction[
            self.layerRangeUnchecked(cell, .ammonium).start + layer
        ];
        const nitrate = self.band_volume_fraction[
            self.layerRangeUnchecked(cell, .nitrate).start + layer
        ];
        const phosphate = self.band_volume_fraction[
            self.layerRangeUnchecked(cell, .phosphate).start + layer
        ];
        return .{
            .ammonium_non_band = 1 - ammonium,
            .ammonium_band = ammonium,
            .nitrate_non_band = 1 - nitrate,
            .nitrate_band = nitrate,
            .phosphate_non_band = 1 - phosphate,
            .phosphate_band = phosphate,
        };
    }

    pub fn prepareFamily(
        self: *State,
        cell: usize,
        token: phase_module.HourToken,
        family: Family,
        layer_geometry: geometry_module.LayerGeometry,
        forcing: geometry_module.Forcing,
        workspace: *geometry_module.Workspace,
    ) !void {
        var view = try self.geometry(cell, family);
        try (try self.coordinator(cell)).prepareFamily(
            token,
            family,
            &view,
            layer_geometry,
            forcing,
            workspace,
        );
        try self.commitGeometryScalars(cell, family, view);
    }

    fn initializeUniformProfile(
        self: *State,
        initialization: Initialization,
    ) void {
        for (0..self.cell_count) |cell| {
            const active_layers =
                initialization.active_layer_count_by_cell[cell];
            const profile_base = cell * self.layer_capacity;
            inline for (std.enums.values(Family)) |family| {
                const family_index = @intFromEnum(family);
                const scalar = cell * family_count + family_index;
                const fraction =
                    initialization.initial_band_fraction_by_family[
                        family_index
                    ];
                const row_spacing =
                    initialization.row_spacing_m_by_cell_family[scalar];
                self.active[scalar] = fraction > 0;
                self.row_spacing_m[scalar] = row_spacing;
                self.upper_edge_depth_m[scalar] = if (fraction > 0)
                    initialization.layer_upper_depth_m[profile_base]
                else
                    0;
                self.lower_edge_depth_m[scalar] = if (fraction > 0)
                    initialization.layer_lower_depth_m[
                        profile_base + active_layers - 1
                    ]
                else
                    0;
                const layers = self.layerRangeUnchecked(cell, family);
                for (0..self.layer_capacity) |layer| {
                    const index = layers.start + layer;
                    const is_active = layer < active_layers;
                    self.band_depth_m[index] = if (fraction > 0 and is_active)
                        initialization.layer_thickness_m[
                            profile_base + layer
                        ]
                    else
                        0;
                    self.band_width_m[index] = if (fraction > 0 and is_active)
                        row_spacing * fraction
                    else
                        0;
                    self.band_volume_fraction[index] =
                        if (is_active) fraction else 0;
                    self.non_band_volume_fraction[index] =
                        if (is_active) 1 - fraction else 1;
                }
            }
        }
    }

    fn scalarIndex(self: *const State, cell: usize, family: Family) !usize {
        if (cell >= self.cell_count) return error.FertilizerBandCellOutOfBounds;
        return cell * family_count + @intFromEnum(family);
    }

    fn layerRange(
        self: *const State,
        cell: usize,
        family: Family,
    ) !LayerRange {
        if (cell >= self.cell_count) return error.FertilizerBandCellOutOfBounds;
        return self.layerRangeUnchecked(cell, family);
    }

    fn layerRangeUnchecked(
        self: *const State,
        cell: usize,
        family: Family,
    ) LayerRange {
        const start =
            (cell * family_count + @intFromEnum(family)) *
            self.layer_capacity;
        return .{ .start = start, .end = start + self.layer_capacity };
    }
};

pub const CheckpointLimits = struct {
    maximum_cells: usize,
    maximum_layers: usize,
};

pub fn writeCheckpoint(writer: anytype, state: *const State) !void {
    try validateState(state);
    try writer.writeAll(checkpoint_magic);
    try writer.writeInt(u32, checkpoint_version, .little);
    try writer.writeInt(u64, @intCast(state.cell_count), .little);
    try writer.writeInt(u64, @intCast(state.layer_capacity), .little);
    for (state.active) |value| try writer.writeByte(@intFromBool(value));
    inline for (.{
        state.row_spacing_m,
        state.upper_edge_depth_m,
        state.lower_edge_depth_m,
        state.band_depth_m,
        state.band_width_m,
        state.band_volume_fraction,
        state.non_band_volume_fraction,
    }) |values| try writeF64Slice(writer, values);
    for (state.coordinators) |*coordinator| {
        const metadata = coordinator.persistentMetadata();
        try writer.writeByte(@intFromEnum(metadata.phase));
        try writer.writeInt(u64, metadata.current_token, .little);
        try writer.writeInt(u64, metadata.last_completed_token, .little);
        try writer.writeInt(u64, @intCast(metadata.next_prepare_family), .little);
        try writer.writeInt(u64, @intCast(metadata.next_consume_family), .little);
        for (metadata.active_by_family) |value|
            try writer.writeByte(@intFromBool(value));
        for (metadata.first_active_layer_by_family) |value|
            try writer.writeInt(u64, @intCast(value), .little);
        for (metadata.last_active_layer_by_family) |value|
            try writer.writeInt(u64, @intCast(value), .little);
        try writeF64Slice(
            writer,
            coordinator.persistentRelativeChanges(),
        );
        for (coordinator.persistentDisappearanceFlags()) |value|
            try writer.writeByte(@intFromBool(value));
    }
}

pub fn readCheckpoint(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    limits: CheckpointLimits,
) !State {
    return readCheckpointInternal(allocator, reader, limits, true);
}

/// Reads one fertilizer-band payload embedded inside a larger versioned
/// checkpoint section and leaves the reader positioned at the next owner.
pub fn readEmbeddedCheckpoint(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    limits: CheckpointLimits,
) !State {
    return readCheckpointInternal(allocator, reader, limits, false);
}

fn readCheckpointInternal(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    limits: CheckpointLimits,
    require_end_of_stream: bool,
) !State {
    if (!std.mem.eql(
        u8,
        try reader.takeArray(checkpoint_magic.len),
        checkpoint_magic,
    )) return error.InvalidFertilizerBandCheckpointMagic;
    if (try reader.takeInt(u32, .little) != checkpoint_version)
        return error.UnsupportedFertilizerBandCheckpointVersion;
    const cells = try bounded(reader, limits.maximum_cells);
    const layers = try bounded(reader, limits.maximum_layers);
    if (cells == 0 or layers == 0)
        return error.InvalidFertilizerBandCheckpointDimensions;
    const active_counts = try allocator.alloc(usize, cells);
    defer allocator.free(active_counts);
    @memset(active_counts, layers);
    const profile_extent = try std.math.mul(usize, cells, layers);
    const upper = try allocator.alloc(f64, profile_extent);
    defer allocator.free(upper);
    const lower = try allocator.alloc(f64, profile_extent);
    defer allocator.free(lower);
    const thickness = try allocator.alloc(f64, profile_extent);
    defer allocator.free(thickness);
    for (0..cells) |cell| for (0..layers) |layer| {
        const index = cell * layers + layer;
        upper[index] = @floatFromInt(layer);
        lower[index] = @floatFromInt(layer + 1);
        thickness[index] = 1;
    };
    const rows = try allocator.alloc(f64, cells * family_count);
    defer allocator.free(rows);
    @memset(rows, 1);
    var state = try State.init(allocator, .{
        .cell_count = cells,
        .layer_capacity = layers,
        .active_layer_count_by_cell = active_counts,
        .layer_upper_depth_m = upper,
        .layer_lower_depth_m = lower,
        .layer_thickness_m = thickness,
        .initial_band_fraction_by_family = .{ 0, 0, 0 },
        .row_spacing_m_by_cell_family = rows,
    });
    errdefer state.deinit();
    for (state.active) |*value| value.* = try readBool(reader);
    inline for (.{
        state.row_spacing_m,
        state.upper_edge_depth_m,
        state.lower_edge_depth_m,
        state.band_depth_m,
        state.band_width_m,
        state.band_volume_fraction,
        state.non_band_volume_fraction,
    }) |values| try readF64Slice(reader, values);
    for (state.coordinators) |*coordinator| {
        const phase = std.enums.fromInt(
            phase_module.Phase,
            try reader.takeByte(),
        ) orelse return error.InvalidFertilizerBandCoordinatorCheckpoint;
        var metadata: phase_module.PersistentMetadata = .{
            .phase = phase,
            .current_token = try reader.takeInt(u64, .little),
            .last_completed_token = try reader.takeInt(u64, .little),
            .next_prepare_family = try bounded(reader, 3),
            .next_consume_family = try bounded(reader, 3),
            .active_by_family = undefined,
            .first_active_layer_by_family = undefined,
            .last_active_layer_by_family = undefined,
        };
        for (&metadata.active_by_family) |*value| value.* = try readBool(reader);
        for (&metadata.first_active_layer_by_family) |*value|
            value.* = try bounded(reader, layers - 1);
        for (&metadata.last_active_layer_by_family) |*value|
            value.* = try bounded(reader, layers - 1);
        const relative = try allocator.alloc(f64, family_count * layers);
        defer allocator.free(relative);
        const disappeared = try allocator.alloc(bool, family_count * layers);
        defer allocator.free(disappeared);
        try readF64Slice(reader, relative);
        for (disappeared) |*value| value.* = try readBool(reader);
        try coordinator.restorePersistent(metadata, relative, disappeared);
    }
    if (require_end_of_stream)
        if (reader.peekByte()) |_| return error.TrailingFertilizerBandCheckpointData else |err| switch (err) {
            error.EndOfStream => {},
            else => return err,
        };
    try validateState(&state);
    return state;
}

fn validateInitialization(initialization: Initialization) !void {
    if (initialization.cell_count == 0 or initialization.layer_capacity == 0 or
        initialization.active_layer_count_by_cell.len !=
            initialization.cell_count or
        initialization.row_spacing_m_by_cell_family.len !=
            initialization.cell_count * family_count)
        return error.InvalidFertilizerBandInitializationDimensions;
    const extent = try std.math.mul(
        usize,
        initialization.cell_count,
        initialization.layer_capacity,
    );
    if (initialization.layer_upper_depth_m.len != extent or
        initialization.layer_lower_depth_m.len != extent or
        initialization.layer_thickness_m.len != extent)
        return error.InvalidFertilizerBandInitializationDimensions;
    for (initialization.initial_band_fraction_by_family) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidInitialFertilizerBandFraction;
    for (0..initialization.cell_count) |cell| {
        const active = initialization.active_layer_count_by_cell[cell];
        if (active == 0 or active > initialization.layer_capacity)
            return error.InvalidFertilizerBandActiveLayerCount;
        for (0..initialization.layer_capacity) |layer| {
            const index = cell * initialization.layer_capacity + layer;
            const upper = initialization.layer_upper_depth_m[index];
            const lower = initialization.layer_lower_depth_m[index];
            const thickness = initialization.layer_thickness_m[index];
            if (!std.math.isFinite(upper) or !std.math.isFinite(lower) or
                !std.math.isFinite(thickness) or upper < 0 or lower <= upper or
                thickness <= 0 or @abs(lower - upper - thickness) > 1e-10)
                return error.InvalidFertilizerBandProfileGeometry;
        }
    }
    for (initialization.row_spacing_m_by_cell_family) |row|
        if (!std.math.isFinite(row) or row <= 0)
            return error.InvalidInitialFertilizerBandRowSpacing;
}

fn validateState(state: *const State) !void {
    const scalar_count = state.cell_count * family_count;
    const layer_count = scalar_count * state.layer_capacity;
    if (state.active.len != scalar_count or
        state.row_spacing_m.len != scalar_count or
        state.upper_edge_depth_m.len != scalar_count or
        state.lower_edge_depth_m.len != scalar_count or
        state.band_depth_m.len != layer_count or
        state.band_width_m.len != layer_count or
        state.band_volume_fraction.len != layer_count or
        state.non_band_volume_fraction.len != layer_count or
        state.coordinators.len != state.cell_count)
        return error.InvalidFertilizerBandStateDimensions;
    inline for (.{
        state.row_spacing_m,
        state.upper_edge_depth_m,
        state.lower_edge_depth_m,
        state.band_depth_m,
        state.band_width_m,
        state.band_volume_fraction,
        state.non_band_volume_fraction,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFertilizerBandState;
    for (state.band_volume_fraction, state.non_band_volume_fraction) |
        band,
        non_band,
    | if (band > 1 or non_band < 0 or non_band > 1 or
        @abs(band + non_band - 1) > 1e-12)
        return error.InvalidFertilizerBandState;
}

fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| try writer.writeInt(u64, @bitCast(value), .little);
}

fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value|
        value.* = @bitCast(try reader.takeInt(u64, .little));
}

fn readBool(reader: *std.Io.Reader) !bool {
    return switch (try reader.takeByte()) {
        0 => false,
        1 => true,
        else => error.InvalidFertilizerBandCheckpointBoolean,
    };
}

fn bounded(reader: *std.Io.Reader, maximum: usize) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > maximum or value > std.math.maxInt(usize))
        return error.FertilizerBandCheckpointLimitExceeded;
    return @intCast(value);
}

test "runtime owner initializes uniform fractions beyond legacy dimensions" {
    const cells = 7;
    const layers = 13;
    var active = [_]usize{layers} ** cells;
    var upper: [cells * layers]f64 = undefined;
    var lower: [cells * layers]f64 = undefined;
    var thickness: [cells * layers]f64 = undefined;
    for (0..cells) |cell| for (0..layers) |layer| {
        const index = cell * layers + layer;
        upper[index] = @as(f64, @floatFromInt(layer)) * 0.1;
        lower[index] = upper[index] + 0.1;
        thickness[index] = 0.1;
    };
    var rows = [_]f64{0.2} ** (cells * family_count);
    var state = try State.init(std.testing.allocator, .{
        .cell_count = cells,
        .layer_capacity = layers,
        .active_layer_count_by_cell = &active,
        .layer_upper_depth_m = &upper,
        .layer_lower_depth_m = &lower,
        .layer_thickness_m = &thickness,
        .initial_band_fraction_by_family = .{ 0.2, 0.3, 0.4 },
        .row_spacing_m_by_cell_family = &rows,
    });
    defer state.deinit();
    const phosphate = try state.geometry(6, .phosphate);
    try std.testing.expectEqual(@as(usize, layers), phosphate.band_depth_m.len);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), phosphate.band_volume_fraction[12], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.08), phosphate.band_width_m[12], 1e-15);
}

test "runscript nutrient parameters initialize family-specific row geometry" {
    const active = [_]usize{2};
    const upper = [_]f64{ 0, 0.1 };
    const lower = [_]f64{ 0.1, 0.3 };
    const thickness = [_]f64{ 0.1, 0.2 };
    var state = try State.initFromPlantNutrientParameters(
        std.testing.allocator,
        1,
        2,
        &active,
        &upper,
        &lower,
        &thickness,
        .{
            .initial_ammonium_band_fraction = 0.2,
            .initial_nitrate_band_fraction = 0.3,
            .initial_phosphate_band_fraction = 0.4,
            .initial_h2po4_fraction = 0.8,
            .initial_ammonium_band_row_spacing_m = 0.5,
            .initial_nitrate_band_row_spacing_m = 0.6,
            .initial_phosphate_band_row_spacing_m = 0.7,
        },
    );
    defer state.deinit();

    const ammonium = try state.geometry(0, .ammonium);
    const nitrate = try state.geometry(0, .nitrate);
    const phosphate = try state.geometry(0, .phosphate);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), ammonium.row_spacing_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.18), nitrate.band_width_m[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.28), phosphate.band_width_m[0], 1e-15);
    const fractions = try state.zoneFractions(0, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), fractions.nitrate_non_band, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), fractions.phosphate_band, 1e-15);
    try std.testing.expectError(
        error.FertilizerBandLayerOutOfBounds,
        state.zoneFractions(0, 2),
    );
}

test "checkpoint round trip preserves geometry and pending phase generation" {
    var active = [_]usize{2};
    const upper = [_]f64{ 0, 0.1 };
    const lower = [_]f64{ 0.1, 0.2 };
    const thickness = [_]f64{ 0.1, 0.1 };
    var rows = [_]f64{ 0.1, 0.1, 0.1 };
    var state = try State.init(std.testing.allocator, .{
        .cell_count = 1,
        .layer_capacity = 2,
        .active_layer_count_by_cell = &active,
        .layer_upper_depth_m = &upper,
        .layer_lower_depth_m = &lower,
        .layer_thickness_m = &thickness,
        .initial_band_fraction_by_family = .{ 0.2, 0.3, 0.4 },
        .row_spacing_m_by_cell_family = &rows,
    });
    defer state.deinit();
    const coordinator = try state.coordinator(0);
    try coordinator.beginHour(.{ .value = 99 });
    var geometry_workspace =
        try geometry_module.Workspace.init(std.testing.allocator, 2);
    defer geometry_workspace.deinit();
    try state.prepareFamily(
        0,
        .{ .value = 99 },
        .ammonium,
        .{
            .upper_depth_m = &upper,
            .lower_depth_m = &lower,
            .thickness_m = &thickness,
            .first_active_layer = 0,
            .last_active_layer = 1,
            .minimum_active_thickness_m = 0,
            .structural_presence_threshold = 0,
        },
        .{
            .diffusivity_m2_per_h = &.{ 0, 0 },
            .tortuosity = &.{ 1, 1 },
            .timestep_h = 1,
        },
        &geometry_workspace,
    );

    var bytes: [8192]u8 = undefined;
    var writer = std.Io.Writer.fixed(&bytes);
    try writeCheckpoint(&writer, &state);
    var reader = std.Io.Reader.fixed(writer.buffered());
    var restored = try readCheckpoint(std.testing.allocator, &reader, .{
        .maximum_cells = 10,
        .maximum_layers = 20,
    });
    defer restored.deinit();
    try std.testing.expectEqualSlices(
        f64,
        state.band_volume_fraction,
        restored.band_volume_fraction,
    );
    try std.testing.expectEqual(
        @as(u64, 99),
        (try restored.coordinator(0)).pendingToken().?.value,
    );
    try std.testing.expectEqual(
        @as(usize, 1),
        (try restored.coordinator(0)).persistentMetadata().next_prepare_family,
    );
}
