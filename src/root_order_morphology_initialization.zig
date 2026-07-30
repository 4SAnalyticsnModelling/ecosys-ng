const std = @import("std");

pub const ElementalMass = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const LayerOrderState = struct {
    secondary_root_number: f64,
    primary_root_length_m: f64,
    primary_root_mass: ElementalMass,
    secondary_root_length_m: f64,
    secondary_root_mass: ElementalMass,
};

pub const PopulationOrderState = struct {
    primary_root_depth_m: f64,
    total_primary_root_mass: ElementalMass,
};

pub const Dimensions = struct {
    population_count: usize,
    layer_count: usize,
    root_order_count: usize,
};

pub const InitializationError = error{
    ExtentOverflow,
    LayerOrderExtentMismatch,
    PopulationOrderExtentMismatch,
    NonFiniteSeedingDepth,
    NegativeSeedingDepth,
};

/// Translates STARTQ lines 812-826.
///
/// Layer-order storage is population-major, layer, then root order.
/// Population-order storage is population-major, then root order.
pub fn initialize(
    dimensions: Dimensions,
    seeding_depth_m: f64,
    layer_order_states: []LayerOrderState,
    population_order_states: []PopulationOrderState,
) InitializationError!void {
    if (!std.math.isFinite(seeding_depth_m)) return error.NonFiniteSeedingDepth;
    if (seeding_depth_m < 0.0) return error.NegativeSeedingDepth;
    const population_layers = std.math.mul(
        usize,
        dimensions.population_count,
        dimensions.layer_count,
    ) catch return error.ExtentOverflow;
    const expected_layer_orders = std.math.mul(
        usize,
        population_layers,
        dimensions.root_order_count,
    ) catch return error.ExtentOverflow;
    const expected_population_orders = std.math.mul(
        usize,
        dimensions.population_count,
        dimensions.root_order_count,
    ) catch return error.ExtentOverflow;
    if (layer_order_states.len != expected_layer_orders) {
        return error.LayerOrderExtentMismatch;
    }
    if (population_order_states.len != expected_population_orders) {
        return error.PopulationOrderExtentMismatch;
    }

    for (0..dimensions.population_count) |population| {
        for (0..dimensions.layer_count) |layer| {
            for (0..dimensions.root_order_count) |root_order| {
                const layer_order_index =
                    (population * dimensions.layer_count + layer) *
                    dimensions.root_order_count + root_order;
                layer_order_states[layer_order_index] = std.mem.zeroes(LayerOrderState);

                // STARTQ repeats these population-order assignments for every
                // soil layer; retaining that order keeps restart semantics exact.
                const population_order_index =
                    population * dimensions.root_order_count + root_order;
                population_order_states[population_order_index] = .{
                    .primary_root_depth_m = seeding_depth_m,
                    .total_primary_root_mass = std.mem.zeroes(ElementalMass),
                };
            }
        }
    }
}

test "runtime populations layers and root orders initialize source domains" {
    const dimensions = Dimensions{
        .population_count = 2,
        .layer_count = 3,
        .root_order_count = 4,
    };
    var layer_orders: [24]LayerOrderState = undefined;
    var population_orders: [8]PopulationOrderState = undefined;
    @memset(std.mem.asBytes(&layer_orders), 0xff);
    @memset(std.mem.asBytes(&population_orders), 0xff);

    try initialize(dimensions, 0.12, &layer_orders, &population_orders);

    for (layer_orders) |state| {
        try std.testing.expectEqual(std.mem.zeroes(LayerOrderState), state);
    }
    for (population_orders) |state| {
        try std.testing.expectEqual(@as(f64, 0.12), state.primary_root_depth_m);
        try std.testing.expectEqual(
            std.mem.zeroes(ElementalMass),
            state.total_primary_root_mass,
        );
    }
}

test "extent mismatch fails before mutation" {
    const dimensions = Dimensions{
        .population_count = 1,
        .layer_count = 2,
        .root_order_count = 2,
    };
    var layer_orders: [3]LayerOrderState = undefined;
    var population_orders: [2]PopulationOrderState = undefined;
    population_orders[0] = .{
        .primary_root_depth_m = 42.0,
        .total_primary_root_mass = std.mem.zeroes(ElementalMass),
    };
    try std.testing.expectError(error.LayerOrderExtentMismatch, initialize(
        dimensions,
        0.1,
        &layer_orders,
        &population_orders,
    ));
    try std.testing.expectEqual(@as(f64, 42.0), population_orders[0].primary_root_depth_m);
}
