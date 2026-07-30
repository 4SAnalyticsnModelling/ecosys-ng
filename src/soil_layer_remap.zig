const std = @import("std");

pub const Unit = enum {
    cubic_meter,
    megagram,
    gram,
    mole,
    megajoule_per_kelvin,
    kelvin,
    mole_per_cubic_meter,
    gram_per_cubic_meter,
    unitless,
};

pub const Quantity = struct {
    name: []const u8,
    unit: Unit,
};

/// Reusable heap workspace for one cell. Counts are runtime values; chemistry,
/// organic-matter, gas, and plant-root sub-pools can all occupy the component
/// axes without compile-time Fortran dimensions.
pub const Workspace = struct {
    allocator: std.mem.Allocator,
    new_layer_capacity: usize,
    extensive_count: usize,
    intensive_count: usize,
    extensive_amounts: []f64,
    intensive_values: []f64,
    intensive_weights: []f64,
    weighted_intensive_sums: []f64,

    pub fn init(allocator: std.mem.Allocator, new_layer_capacity: usize, extensive_count: usize, intensive_count: usize) !Workspace {
        if (new_layer_capacity == 0 or extensive_count == 0) return error.InvalidLayerRemapDimensions;
        const extensive_len = try std.math.mul(usize, new_layer_capacity, extensive_count);
        const intensive_len = try std.math.mul(usize, new_layer_capacity, intensive_count);
        const extensive_amounts = try allocator.alloc(f64, extensive_len);
        errdefer allocator.free(extensive_amounts);
        const intensive_values = try allocator.alloc(f64, intensive_len);
        errdefer allocator.free(intensive_values);
        const intensive_weights = try allocator.alloc(f64, intensive_len);
        errdefer allocator.free(intensive_weights);
        const weighted_intensive_sums = try allocator.alloc(f64, intensive_len);
        errdefer allocator.free(weighted_intensive_sums);
        @memset(extensive_amounts, 0);
        @memset(intensive_values, 0);
        @memset(intensive_weights, 0);
        @memset(weighted_intensive_sums, 0);
        return .{ .allocator = allocator, .new_layer_capacity = new_layer_capacity, .extensive_count = extensive_count, .intensive_count = intensive_count, .extensive_amounts = extensive_amounts, .intensive_values = intensive_values, .intensive_weights = intensive_weights, .weighted_intensive_sums = weighted_intensive_sums };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.weighted_intensive_sums);
        self.allocator.free(self.intensive_weights);
        self.allocator.free(self.intensive_values);
        self.allocator.free(self.extensive_amounts);
        self.* = undefined;
    }
};

pub const Inputs = struct {
    old_boundary_depth_m: []const f64,
    new_boundary_depth_m: []const f64,
    extensive_quantities: []const Quantity,
    old_extensive_amounts: []const f64,
    intensive_quantities: []const Quantity,
    old_intensive_values: []const f64,
    old_intensive_weights: []const f64,
    minimum_layer_thickness_m: f64,
    conservation_relative_tolerance: f64,
    conservation_absolute_tolerance: f64,
};

/// Ports REDIST's repeated `destination += FX * source; source *= 1-FX`
/// operations using exact old/new layer overlap. Extensive state is assumed
/// homogeneous within each old layer, matching the Fortran fraction model.
/// Intensive state is remapped by its caller-supplied extensive weight (e.g.
/// temperature by heat capacity, concentration by water volume).
pub fn remap(workspace: *Workspace, inputs: Inputs) !void {
    if (inputs.old_boundary_depth_m.len < 2 or inputs.new_boundary_depth_m.len < 2) return error.InvalidLayerRemapDimensions;
    const old_layers = inputs.old_boundary_depth_m.len - 1;
    const new_layers = inputs.new_boundary_depth_m.len - 1;
    if (new_layers > workspace.new_layer_capacity or inputs.extensive_quantities.len != workspace.extensive_count or inputs.intensive_quantities.len != workspace.intensive_count or inputs.old_extensive_amounts.len != old_layers * workspace.extensive_count or inputs.old_intensive_values.len != old_layers * workspace.intensive_count or inputs.old_intensive_weights.len != old_layers * workspace.intensive_count) return error.LayerRemapDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_layer_thickness_m) or inputs.minimum_layer_thickness_m <= 0 or !std.math.isFinite(inputs.conservation_relative_tolerance) or inputs.conservation_relative_tolerance < 0 or !std.math.isFinite(inputs.conservation_absolute_tolerance) or inputs.conservation_absolute_tolerance < 0) return error.InvalidLayerRemapTolerance;
    try validateBoundaries(inputs.old_boundary_depth_m, inputs.minimum_layer_thickness_m);
    try validateBoundaries(inputs.new_boundary_depth_m, inputs.minimum_layer_thickness_m);
    const extent_tolerance = tolerance(inputs.old_boundary_depth_m[0], inputs.old_boundary_depth_m[old_layers], inputs.conservation_relative_tolerance, inputs.conservation_absolute_tolerance);
    if (@abs(inputs.old_boundary_depth_m[0] - inputs.new_boundary_depth_m[0]) > extent_tolerance or @abs(inputs.old_boundary_depth_m[old_layers] - inputs.new_boundary_depth_m[new_layers]) > extent_tolerance) return error.LayerRemapExtentMismatch;
    for (inputs.extensive_quantities) |quantity| if (quantity.name.len == 0) return error.EmptyLayerQuantityName;
    for (inputs.intensive_quantities) |quantity| if (quantity.name.len == 0) return error.EmptyLayerQuantityName;
    for (inputs.old_extensive_amounts) |amount| if (!std.math.isFinite(amount) or amount < 0) return error.InvalidLayerExtensiveAmount;
    for (inputs.old_intensive_values, inputs.old_intensive_weights) |value, weight| if (!std.math.isFinite(value) or !std.math.isFinite(weight) or weight < 0) return error.InvalidLayerIntensiveState;

    @memset(workspace.extensive_amounts, 0);
    @memset(workspace.intensive_values, 0);
    @memset(workspace.intensive_weights, 0);
    @memset(workspace.weighted_intensive_sums, 0);
    var old_layer: usize = 0;
    var new_layer: usize = 0;
    while (old_layer < old_layers and new_layer < new_layers) {
        const old_top = inputs.old_boundary_depth_m[old_layer];
        const old_bottom = inputs.old_boundary_depth_m[old_layer + 1];
        const new_top = inputs.new_boundary_depth_m[new_layer];
        const new_bottom = inputs.new_boundary_depth_m[new_layer + 1];
        const overlap_m = @max(0, @min(old_bottom, new_bottom) - @max(old_top, new_top));
        if (overlap_m > 0) {
            const old_fraction = overlap_m / (old_bottom - old_top);
            for (0..workspace.extensive_count) |component| workspace.extensive_amounts[new_layer * workspace.extensive_count + component] += old_fraction * inputs.old_extensive_amounts[old_layer * workspace.extensive_count + component];
            for (0..workspace.intensive_count) |component| {
                const old_index = old_layer * workspace.intensive_count + component;
                const new_index = new_layer * workspace.intensive_count + component;
                const transferred_weight = old_fraction * inputs.old_intensive_weights[old_index];
                workspace.intensive_weights[new_index] += transferred_weight;
                workspace.weighted_intensive_sums[new_index] += transferred_weight * inputs.old_intensive_values[old_index];
            }
        }
        if (old_bottom <= new_bottom) old_layer += 1;
        if (new_bottom <= old_bottom) new_layer += 1;
    }
    for (0..new_layers) |layer| for (0..workspace.intensive_count) |component| {
        const index = layer * workspace.intensive_count + component;
        if (workspace.intensive_weights[index] > 0) workspace.intensive_values[index] = workspace.weighted_intensive_sums[index] / workspace.intensive_weights[index];
    };
    try verifyConservation(workspace, inputs, old_layers, new_layers);
}

fn validateBoundaries(boundaries: []const f64, minimum_layer_thickness_m: f64) !void {
    for (boundaries) |depth| if (!std.math.isFinite(depth)) return error.NonFiniteLayerBoundary;
    for (0..boundaries.len - 1) |layer| if (boundaries[layer + 1] - boundaries[layer] < minimum_layer_thickness_m) return error.InvalidLayerBoundaryOrder;
}

fn verifyConservation(workspace: *const Workspace, inputs: Inputs, old_layers: usize, new_layers: usize) !void {
    for (0..workspace.extensive_count) |component| {
        var old_total: f64 = 0;
        var new_total: f64 = 0;
        for (0..old_layers) |layer| old_total += inputs.old_extensive_amounts[layer * workspace.extensive_count + component];
        for (0..new_layers) |layer| new_total += workspace.extensive_amounts[layer * workspace.extensive_count + component];
        if (@abs(new_total - old_total) > tolerance(0, old_total, inputs.conservation_relative_tolerance, inputs.conservation_absolute_tolerance)) return error.LayerRemapConservationFailure;
    }
    for (0..workspace.intensive_count) |component| {
        var old_weight: f64 = 0;
        var old_weighted_value: f64 = 0;
        var new_weight: f64 = 0;
        var new_weighted_value: f64 = 0;
        for (0..old_layers) |layer| {
            const index = layer * workspace.intensive_count + component;
            old_weight += inputs.old_intensive_weights[index];
            old_weighted_value += inputs.old_intensive_weights[index] * inputs.old_intensive_values[index];
        }
        for (0..new_layers) |layer| {
            const index = layer * workspace.intensive_count + component;
            new_weight += workspace.intensive_weights[index];
            new_weighted_value += workspace.intensive_weights[index] * workspace.intensive_values[index];
        }
        if (@abs(new_weight - old_weight) > tolerance(0, old_weight, inputs.conservation_relative_tolerance, inputs.conservation_absolute_tolerance) or @abs(new_weighted_value - old_weighted_value) > tolerance(0, old_weighted_value, inputs.conservation_relative_tolerance, inputs.conservation_absolute_tolerance)) return error.LayerRemapConservationFailure;
    }
}

fn tolerance(a: f64, b: f64, relative: f64, absolute: f64) f64 {
    return absolute + relative * @max(@abs(a), @abs(b));
}

test "splitting a layer conserves runtime extensive pools" {
    var workspace = try Workspace.init(std.testing.allocator, 2, 2, 0);
    defer workspace.deinit();
    const quantities = [_]Quantity{ .{ .name = "liquid_water", .unit = .cubic_meter }, .{ .name = "nitrate", .unit = .mole } };
    try remap(&workspace, .{ .old_boundary_depth_m = &.{ 0, 1 }, .new_boundary_depth_m = &.{ 0, 0.25, 1 }, .extensive_quantities = &quantities, .old_extensive_amounts = &.{ 8, 20 }, .intensive_quantities = &.{}, .old_intensive_values = &.{}, .old_intensive_weights = &.{}, .minimum_layer_thickness_m = 1e-9, .conservation_relative_tolerance = 1e-12, .conservation_absolute_tolerance = 1e-14 });
    try std.testing.expectApproxEqAbs(@as(f64, 2), workspace.extensive_amounts[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 5), workspace.extensive_amounts[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 6), workspace.extensive_amounts[2], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 15), workspace.extensive_amounts[3], 1e-14);
}

test "merging layers mixes temperature by heat capacity and conserves energy" {
    var workspace = try Workspace.init(std.testing.allocator, 1, 1, 1);
    defer workspace.deinit();
    const extensive = [_]Quantity{.{ .name = "soil_mass", .unit = .megagram }};
    const intensive = [_]Quantity{.{ .name = "temperature", .unit = .kelvin }};
    try remap(&workspace, .{ .old_boundary_depth_m = &.{ 0, 0.4, 1 }, .new_boundary_depth_m = &.{ 0, 1 }, .extensive_quantities = &extensive, .old_extensive_amounts = &.{ 4, 6 }, .intensive_quantities = &intensive, .old_intensive_values = &.{ 280, 300 }, .old_intensive_weights = &.{ 2, 3 }, .minimum_layer_thickness_m = 1e-9, .conservation_relative_tolerance = 1e-12, .conservation_absolute_tolerance = 1e-14 });
    try std.testing.expectApproxEqAbs(@as(f64, 10), workspace.extensive_amounts[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 292), workspace.intensive_values[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1460), workspace.intensive_weights[0] * workspace.intensive_values[0], 1e-12);
}

test "mismatched profile extents fail before remapping" {
    var workspace = try Workspace.init(std.testing.allocator, 1, 1, 0);
    defer workspace.deinit();
    const quantities = [_]Quantity{.{ .name = "clay", .unit = .megagram }};
    try std.testing.expectError(error.LayerRemapExtentMismatch, remap(&workspace, .{ .old_boundary_depth_m = &.{ 0, 1 }, .new_boundary_depth_m = &.{ 0, 0.9 }, .extensive_quantities = &quantities, .old_extensive_amounts = &.{1}, .intensive_quantities = &.{}, .old_intensive_values = &.{}, .old_intensive_weights = &.{}, .minimum_layer_thickness_m = 1e-9, .conservation_relative_tolerance = 1e-12, .conservation_absolute_tolerance = 1e-14 }));
}
