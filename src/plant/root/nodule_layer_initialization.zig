const std = @import("std");

pub const ElementalMass = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const LayerState = struct {
    nonstructural_pool: ElementalMass,
    nodule_mass: ElementalMass,
    nitrogen_fixation_uptake_g_n_per_timestep: f64,
};

pub const KineticPools = struct {
    carbon_g_c: []f64,
    nitrogen_g_n: []f64,
    phosphorus_g_p: []f64,
};

pub const Dimensions = struct {
    layer_count: usize,
    compartment_count: usize,
    kinetic_pool_count: usize,
};

pub const InitializationError = error{
    ExtentOverflow,
    LayerExtentMismatch,
    KineticExtentMismatch,
};

/// Translates `startq.f` lines 827--841.
///
/// Kinetic storage is layer-major, then compartment, then kinetic pool.
/// This module is root-owned and has no mycorrhizal population dimension.
pub fn initialize(
    dimensions: Dimensions,
    layer_states: []LayerState,
    kinetic_pools: KineticPools,
) InitializationError!void {
    if (layer_states.len != dimensions.layer_count) {
        return error.LayerExtentMismatch;
    }
    const layer_compartments = std.math.mul(
        usize,
        dimensions.layer_count,
        dimensions.compartment_count,
    ) catch return error.ExtentOverflow;
    const expected_kinetic_count = std.math.mul(
        usize,
        layer_compartments,
        dimensions.kinetic_pool_count,
    ) catch return error.ExtentOverflow;
    inline for (std.meta.fields(KineticPools)) |field| {
        if (@field(kinetic_pools, field.name).len != expected_kinetic_count) {
            return error.KineticExtentMismatch;
        }
    }

    for (0..dimensions.layer_count) |layer| {
        for (0..dimensions.compartment_count) |compartment| {
            for (0..dimensions.kinetic_pool_count) |kinetic_pool| {
                const index =
                    (layer * dimensions.compartment_count + compartment) *
                    dimensions.kinetic_pool_count + kinetic_pool;
                kinetic_pools.carbon_g_c[index] = 0.0;
                kinetic_pools.nitrogen_g_n[index] = 0.0;
                kinetic_pools.phosphorus_g_p[index] = 0.0;
            }
        }
        layer_states[layer] = std.mem.zeroes(LayerState);
    }
}

test "runtime nodule layers compartments and kinetic pools reset" {
    const dimensions = Dimensions{
        .layer_count = 3,
        .compartment_count = 2,
        .kinetic_pool_count = 4,
    };
    var layers: [3]LayerState = undefined;
    var carbon: [24]f64 = undefined;
    var nitrogen: [24]f64 = undefined;
    var phosphorus: [24]f64 = undefined;
    @memset(std.mem.asBytes(&layers), 0xff);
    @memset(&carbon, 7.0);
    @memset(&nitrogen, 8.0);
    @memset(&phosphorus, 9.0);

    try initialize(dimensions, &layers, .{
        .carbon_g_c = &carbon,
        .nitrogen_g_n = &nitrogen,
        .phosphorus_g_p = &phosphorus,
    });

    for (layers) |layer| try std.testing.expectEqual(std.mem.zeroes(LayerState), layer);
    for (carbon) |value| try std.testing.expectEqual(@as(f64, 0.0), value);
    for (nitrogen) |value| try std.testing.expectEqual(@as(f64, 0.0), value);
    for (phosphorus) |value| try std.testing.expectEqual(@as(f64, 0.0), value);
}

test "kinetic extent mismatch fails before layer mutation" {
    const dimensions = Dimensions{
        .layer_count = 1,
        .compartment_count = 2,
        .kinetic_pool_count = 2,
    };
    var layers = [_]LayerState{std.mem.zeroes(LayerState)};
    layers[0].nitrogen_fixation_uptake_g_n_per_timestep = 42.0;
    var carbon: [3]f64 = undefined;
    var nitrogen: [4]f64 = undefined;
    var phosphorus: [4]f64 = undefined;
    try std.testing.expectError(error.KineticExtentMismatch, initialize(
        dimensions,
        &layers,
        .{
            .carbon_g_c = &carbon,
            .nitrogen_g_n = &nitrogen,
            .phosphorus_g_p = &phosphorus,
        },
    ));
    try std.testing.expectEqual(
        @as(f64, 42.0),
        layers[0].nitrogen_fixation_uptake_g_n_per_timestep,
    );
}
