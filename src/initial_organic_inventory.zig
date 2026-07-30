const std = @import("std");

pub const LayerKind = enum { surface_litter, mineral_soil };

pub const Dimensions = struct {
    organic_pool_count: usize,
    microbial_pool_count: usize,
    microbial_population_count: usize,
    microbial_component_count: usize,
    residue_component_count: usize,
    structural_component_count: usize,
    litter_pool_count: usize,
    detrital_pool_count: usize,
    noncharcoal_structural_component_count: usize,
};

pub const ElementPools = struct {
    carbon_g_c: []const f64,
    nitrogen_g_n: []const f64,
    phosphorus_g_p: []const f64,
};

pub const MobilePools = struct {
    micropore_dissolved: ElementPools,
    macropore_dissolved: ElementPools,
    adsorbed: ElementPools,
    micropore_acetate_g_c: []const f64,
    macropore_acetate_g_c: []const f64,
    adsorbed_acetate_g_c: []const f64,
};

pub const Inputs = struct {
    layer_kind: LayerKind,
    dimensions: Dimensions,
    microbial: ElementPools,
    microbial_residue: ElementPools,
    mobile: MobilePools,
    structural: ElementPools,
};

pub const InventoryState = struct {
    total_noncharcoal_carbon_g_c: *f64,
    total_noncharcoal_nitrogen_g_n: *f64,
    total_noncharcoal_phosphorus_g_p: *f64,
    charcoal_carbon_g_c: *f64,
    charcoal_nitrogen_g_n: *f64,
    charcoal_phosphorus_g_p: *f64,
    detrital_carbon_g_c: *f64,
    detrital_nitrogen_g_n: *f64,
    detrital_phosphorus_g_p: *f64,
    organic_carbon_snapshot_g_c: *f64,
    litter_carbon_g_c: *f64,
    charcoal_carbon_change_g_c: *f64,
    /// Persistent per-cell surface inventory; mineral layers add microbial C.
    surface_carbon_g_c_by_microbial_pool: []f64,
};

pub const TransformationDiagnostics = struct {
    oxygen_response: []f64,
    maximum_rate_4: []f64,
    maximum_rate_3: []f64,
    maximum_rate_2: []f64,
    maximum_rate_1: []f64,
    humus_nitrogen_inhibition: []f64,
    organic_nitrogen_inhibition: []f64,
    organic_phosphorus_inhibition: []f64,
    surface_humus_nitrogen_inhibition: []f64,
    surface_organic_nitrogen_inhibition: []f64,
    surface_organic_phosphorus_inhibition: []f64,
};

const Totals = struct {
    carbon_g_c: f64 = 0.0,
    nitrogen_g_n: f64 = 0.0,
    phosphorus_g_p: f64 = 0.0,
    charcoal_carbon_g_c: f64 = 0.0,
    charcoal_nitrogen_g_n: f64 = 0.0,
    charcoal_phosphorus_g_p: f64 = 0.0,
    detrital_carbon_g_c: f64 = 0.0,
    detrital_nitrogen_g_n: f64 = 0.0,
    detrital_phosphorus_g_p: f64 = 0.0,
};

fn mul(a: usize, b: usize) !usize {
    return std.math.mul(usize, a, b) catch
        return error.InitialOrganicInventoryDimensionOverflow;
}

fn dimensions(inputs: Inputs, state: InventoryState, diagnostics: TransformationDiagnostics) !void {
    const d = inputs.dimensions;
    if (d.organic_pool_count == 0 or d.microbial_pool_count == 0 or
        d.microbial_pool_count < d.organic_pool_count or
        d.microbial_population_count == 0 or d.microbial_component_count == 0 or
        d.residue_component_count == 0 or d.structural_component_count == 0 or
        d.litter_pool_count > d.organic_pool_count or
        d.detrital_pool_count > d.organic_pool_count or
        d.noncharcoal_structural_component_count >
            d.structural_component_count)
        return error.InvalidInitialOrganicInventoryDimensions;
    const microbial_count = try mul(
        try mul(d.microbial_pool_count, d.microbial_population_count),
        d.microbial_component_count,
    );
    const population_pool_count =
        try mul(d.microbial_pool_count, d.microbial_population_count);
    const residue_count =
        try mul(d.organic_pool_count, d.residue_component_count);
    const structural_count =
        try mul(d.organic_pool_count, d.structural_component_count);
    inline for (.{
        inputs.microbial.carbon_g_c,
        inputs.microbial.nitrogen_g_n,
        inputs.microbial.phosphorus_g_p,
    }) |values| if (values.len != microbial_count)
        return error.InitialOrganicInventoryDimensionMismatch;
    inline for (.{
        inputs.microbial_residue.carbon_g_c,
        inputs.microbial_residue.nitrogen_g_n,
        inputs.microbial_residue.phosphorus_g_p,
    }) |values| if (values.len != residue_count)
        return error.InitialOrganicInventoryDimensionMismatch;
    inline for (.{
        inputs.structural.carbon_g_c,
        inputs.structural.nitrogen_g_n,
        inputs.structural.phosphorus_g_p,
    }) |values| if (values.len != structural_count)
        return error.InitialOrganicInventoryDimensionMismatch;
    inline for (.{
        inputs.mobile.micropore_dissolved.carbon_g_c,
        inputs.mobile.micropore_dissolved.nitrogen_g_n,
        inputs.mobile.micropore_dissolved.phosphorus_g_p,
        inputs.mobile.macropore_dissolved.carbon_g_c,
        inputs.mobile.macropore_dissolved.nitrogen_g_n,
        inputs.mobile.macropore_dissolved.phosphorus_g_p,
        inputs.mobile.adsorbed.carbon_g_c,
        inputs.mobile.adsorbed.nitrogen_g_n,
        inputs.mobile.adsorbed.phosphorus_g_p,
        inputs.mobile.micropore_acetate_g_c,
        inputs.mobile.macropore_acetate_g_c,
        inputs.mobile.adsorbed_acetate_g_c,
    }) |values| if (values.len != d.organic_pool_count)
        return error.InitialOrganicInventoryDimensionMismatch;
    inline for (.{
        diagnostics.oxygen_response,
        diagnostics.maximum_rate_4,
        diagnostics.maximum_rate_3,
        diagnostics.maximum_rate_2,
        diagnostics.maximum_rate_1,
        diagnostics.humus_nitrogen_inhibition,
        diagnostics.organic_nitrogen_inhibition,
        diagnostics.organic_phosphorus_inhibition,
    }) |values| if (values.len != population_pool_count)
        return error.InitialOrganicInventoryDimensionMismatch;
    inline for (.{
        diagnostics.surface_humus_nitrogen_inhibition,
        diagnostics.surface_organic_nitrogen_inhibition,
        diagnostics.surface_organic_phosphorus_inhibition,
    }) |values| if (values.len != population_pool_count)
        return error.InitialOrganicInventoryDimensionMismatch;
    if (state.surface_carbon_g_c_by_microbial_pool.len !=
        d.microbial_pool_count)
        return error.InitialOrganicInventoryDimensionMismatch;
}

fn validatePools(inputs: Inputs, state: InventoryState) !void {
    inline for (.{
        inputs.microbial.carbon_g_c,
        inputs.microbial.nitrogen_g_n,
        inputs.microbial.phosphorus_g_p,
        inputs.microbial_residue.carbon_g_c,
        inputs.microbial_residue.nitrogen_g_n,
        inputs.microbial_residue.phosphorus_g_p,
        inputs.mobile.micropore_dissolved.carbon_g_c,
        inputs.mobile.micropore_dissolved.nitrogen_g_n,
        inputs.mobile.micropore_dissolved.phosphorus_g_p,
        inputs.mobile.macropore_dissolved.carbon_g_c,
        inputs.mobile.macropore_dissolved.nitrogen_g_n,
        inputs.mobile.macropore_dissolved.phosphorus_g_p,
        inputs.mobile.adsorbed.carbon_g_c,
        inputs.mobile.adsorbed.nitrogen_g_n,
        inputs.mobile.adsorbed.phosphorus_g_p,
        inputs.mobile.micropore_acetate_g_c,
        inputs.mobile.macropore_acetate_g_c,
        inputs.mobile.adsorbed_acetate_g_c,
        inputs.structural.carbon_g_c,
        inputs.structural.nitrogen_g_n,
        inputs.structural.phosphorus_g_p,
        state.surface_carbon_g_c_by_microbial_pool,
    }) |values| for (values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteInitialOrganicInventoryInput;
        if (value < 0.0) return error.InvalidInitialOrganicInventoryInput;
    };
}

fn microbialIndex(d: Dimensions, pool: usize, population: usize, component: usize) usize {
    return (pool * d.microbial_population_count + population) *
        d.microbial_component_count + component;
}

fn accumulate(
    inputs: Inputs,
    surface_carbon: []f64,
) !Totals {
    const d = inputs.dimensions;
    var totals: Totals = .{};
    if (inputs.layer_kind == .surface_litter) @memset(surface_carbon, 0.0);

    for (0..d.microbial_pool_count) |pool| {
        for (0..d.microbial_population_count) |population| {
            if (pool < d.litter_pool_count) {
                const final_component = d.microbial_component_count - 1;
                const index =
                    microbialIndex(d, pool, population, final_component);
                totals.detrital_carbon_g_c += inputs.microbial.carbon_g_c[index];
                totals.detrital_nitrogen_g_n +=
                    inputs.microbial.nitrogen_g_n[index];
                totals.detrital_phosphorus_g_p +=
                    inputs.microbial.phosphorus_g_p[index];
            }
            for (0..d.microbial_component_count) |component| {
                const index = microbialIndex(d, pool, population, component);
                totals.carbon_g_c += inputs.microbial.carbon_g_c[index];
                totals.nitrogen_g_n += inputs.microbial.nitrogen_g_n[index];
                totals.phosphorus_g_p += inputs.microbial.phosphorus_g_p[index];
                if (pool < d.detrital_pool_count)
                    totals.detrital_carbon_g_c +=
                        inputs.microbial.carbon_g_c[index];
                surface_carbon[pool] += inputs.microbial.carbon_g_c[index];
            }
        }
    }
    for (0..d.organic_pool_count) |pool| {
        for (0..d.residue_component_count) |component| {
            const index = pool * d.residue_component_count + component;
            totals.carbon_g_c += inputs.microbial_residue.carbon_g_c[index];
            totals.nitrogen_g_n += inputs.microbial_residue.nitrogen_g_n[index];
            totals.phosphorus_g_p +=
                inputs.microbial_residue.phosphorus_g_p[index];
            if (pool < d.detrital_pool_count) {
                totals.detrital_carbon_g_c +=
                    inputs.microbial_residue.carbon_g_c[index];
                totals.detrital_nitrogen_g_n +=
                    inputs.microbial_residue.nitrogen_g_n[index];
                totals.detrital_phosphorus_g_p +=
                    inputs.microbial_residue.phosphorus_g_p[index];
            }
            if (inputs.layer_kind == .surface_litter)
                surface_carbon[pool] +=
                    inputs.microbial_residue.carbon_g_c[index];
        }
        const dissolved = inputs.mobile.micropore_dissolved;
        const macropore = inputs.mobile.macropore_dissolved;
        const adsorbed = inputs.mobile.adsorbed;
        const acetate_c = inputs.mobile.micropore_acetate_g_c[pool] +
            inputs.mobile.macropore_acetate_g_c[pool] +
            inputs.mobile.adsorbed_acetate_g_c[pool];
        const mobile_c = dissolved.carbon_g_c[pool] +
            macropore.carbon_g_c[pool] + adsorbed.carbon_g_c[pool] + acetate_c;
        totals.carbon_g_c += mobile_c;
        totals.nitrogen_g_n += dissolved.nitrogen_g_n[pool] +
            macropore.nitrogen_g_n[pool] + adsorbed.nitrogen_g_n[pool];
        totals.phosphorus_g_p += dissolved.phosphorus_g_p[pool] +
            macropore.phosphorus_g_p[pool] + adsorbed.phosphorus_g_p[pool];
        // STARTS line 1542 adds micropore and macropore acetate a second time.
        totals.carbon_g_c += inputs.mobile.micropore_acetate_g_c[pool] +
            inputs.mobile.macropore_acetate_g_c[pool];
        if (pool < d.detrital_pool_count) {
            totals.detrital_carbon_g_c += mobile_c;
            totals.detrital_nitrogen_g_n += dissolved.nitrogen_g_n[pool] +
                macropore.nitrogen_g_n[pool] + adsorbed.nitrogen_g_n[pool];
            totals.detrital_phosphorus_g_p += dissolved.phosphorus_g_p[pool] +
                macropore.phosphorus_g_p[pool] +
                adsorbed.phosphorus_g_p[pool];
        }
        if (inputs.layer_kind == .surface_litter)
            surface_carbon[pool] += mobile_c;

        for (0..d.structural_component_count) |component| {
            const index = pool * d.structural_component_count + component;
            if (component < d.noncharcoal_structural_component_count) {
                totals.carbon_g_c += inputs.structural.carbon_g_c[index];
                totals.nitrogen_g_n += inputs.structural.nitrogen_g_n[index];
                totals.phosphorus_g_p += inputs.structural.phosphorus_g_p[index];
            } else {
                totals.charcoal_carbon_g_c +=
                    inputs.structural.carbon_g_c[index];
                totals.charcoal_nitrogen_g_n +=
                    inputs.structural.nitrogen_g_n[index];
                totals.charcoal_phosphorus_g_p +=
                    inputs.structural.phosphorus_g_p[index];
            }
            if (pool < d.detrital_pool_count) {
                totals.detrital_carbon_g_c +=
                    inputs.structural.carbon_g_c[index];
                totals.detrital_nitrogen_g_n +=
                    inputs.structural.nitrogen_g_n[index];
                totals.detrital_phosphorus_g_p +=
                    inputs.structural.phosphorus_g_p[index];
            }
            if (inputs.layer_kind == .surface_litter)
                surface_carbon[pool] += inputs.structural.carbon_g_c[index];
        }
    }
    inline for (@typeInfo(Totals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals, field.name)))
            return error.NonFiniteInitialOrganicInventoryResult;
    for (surface_carbon) |value| if (!std.math.isFinite(value))
        return error.NonFiniteInitialOrganicInventoryResult;
    return totals;
}

/// Exact source-order translation of legacy `STARTS` lines 1478--1584.
pub fn initialize(
    allocator: std.mem.Allocator,
    inventory: InventoryState,
    diagnostics: TransformationDiagnostics,
    inputs: Inputs,
) !void {
    try dimensions(inputs, inventory, diagnostics);
    try validatePools(inputs, inventory);
    const staged_surface = try allocator.dupe(
        f64,
        inventory.surface_carbon_g_c_by_microbial_pool,
    );
    defer allocator.free(staged_surface);
    const totals = try accumulate(inputs, staged_surface);

    inventory.total_noncharcoal_carbon_g_c.* = totals.carbon_g_c;
    inventory.total_noncharcoal_nitrogen_g_n.* = totals.nitrogen_g_n;
    inventory.total_noncharcoal_phosphorus_g_p.* = totals.phosphorus_g_p;
    inventory.charcoal_carbon_g_c.* = totals.charcoal_carbon_g_c;
    inventory.charcoal_nitrogen_g_n.* = totals.charcoal_nitrogen_g_n;
    inventory.charcoal_phosphorus_g_p.* = totals.charcoal_phosphorus_g_p;
    inventory.detrital_carbon_g_c.* = totals.detrital_carbon_g_c;
    inventory.detrital_nitrogen_g_n.* = totals.detrital_nitrogen_g_n;
    inventory.detrital_phosphorus_g_p.* = totals.detrital_phosphorus_g_p;
    inventory.organic_carbon_snapshot_g_c.* = totals.carbon_g_c;
    inventory.litter_carbon_g_c.* = totals.detrital_carbon_g_c;
    inventory.charcoal_carbon_change_g_c.* = 0.0;
    @memcpy(inventory.surface_carbon_g_c_by_microbial_pool, staged_surface);

    inline for (.{
        diagnostics.oxygen_response,
        diagnostics.maximum_rate_4,
        diagnostics.maximum_rate_3,
        diagnostics.maximum_rate_2,
        diagnostics.maximum_rate_1,
        diagnostics.humus_nitrogen_inhibition,
        diagnostics.organic_nitrogen_inhibition,
        diagnostics.organic_phosphorus_inhibition,
    }) |values| @memset(values, 0.0);
    if (inputs.layer_kind == .surface_litter) {
        inline for (.{
            diagnostics.surface_humus_nitrogen_inhibition,
            diagnostics.surface_organic_nitrogen_inhibition,
            diagnostics.surface_organic_phosphorus_inhibition,
        }) |values| @memset(values, 0.0);
    }
}

fn zeroPools(count: usize) ElementPools {
    return .{
        .carbon_g_c = (&[_]f64{0} ** 64)[0..count],
        .nitrogen_g_n = (&[_]f64{0} ** 64)[0..count],
        .phosphorus_g_p = (&[_]f64{0} ** 64)[0..count],
    };
}

test "STARTS aggregation preserves double acetate and detrital order" {
    const d: Dimensions = .{
        .organic_pool_count = 1,
        .microbial_pool_count = 1,
        .microbial_population_count = 1,
        .microbial_component_count = 2,
        .residue_component_count = 1,
        .structural_component_count = 2,
        .litter_pool_count = 1,
        .detrital_pool_count = 1,
        .noncharcoal_structural_component_count = 1,
    };
    const microbial: ElementPools = .{
        .carbon_g_c = &.{ 1, 2 },
        .nitrogen_g_n = &.{ 0.1, 0.2 },
        .phosphorus_g_p = &.{ 0.01, 0.02 },
    };
    const residue: ElementPools = .{
        .carbon_g_c = &.{3},
        .nitrogen_g_n = &.{0.3},
        .phosphorus_g_p = &.{0.03},
    };
    const structural: ElementPools = .{
        .carbon_g_c = &.{ 4, 5 },
        .nitrogen_g_n = &.{ 0.4, 0.5 },
        .phosphorus_g_p = &.{ 0.04, 0.05 },
    };
    const dissolved: ElementPools = .{
        .carbon_g_c = &.{6},
        .nitrogen_g_n = &.{0.6},
        .phosphorus_g_p = &.{0.06},
    };
    const zeros = zeroPools(1);
    var scalars = [_]f64{9.0} ** 12;
    var surface = [_]f64{9.0};
    var diagnostic = [_]f64{9.0} ** 11;
    try initialize(std.testing.allocator, .{
        .total_noncharcoal_carbon_g_c = &scalars[0],
        .total_noncharcoal_nitrogen_g_n = &scalars[1],
        .total_noncharcoal_phosphorus_g_p = &scalars[2],
        .charcoal_carbon_g_c = &scalars[3],
        .charcoal_nitrogen_g_n = &scalars[4],
        .charcoal_phosphorus_g_p = &scalars[5],
        .detrital_carbon_g_c = &scalars[6],
        .detrital_nitrogen_g_n = &scalars[7],
        .detrital_phosphorus_g_p = &scalars[8],
        .organic_carbon_snapshot_g_c = &scalars[9],
        .litter_carbon_g_c = &scalars[10],
        .charcoal_carbon_change_g_c = &scalars[11],
        .surface_carbon_g_c_by_microbial_pool = &surface,
    }, .{
        .oxygen_response = diagnostic[0..1],
        .maximum_rate_4 = diagnostic[1..2],
        .maximum_rate_3 = diagnostic[2..3],
        .maximum_rate_2 = diagnostic[3..4],
        .maximum_rate_1 = diagnostic[4..5],
        .humus_nitrogen_inhibition = diagnostic[5..6],
        .organic_nitrogen_inhibition = diagnostic[6..7],
        .organic_phosphorus_inhibition = diagnostic[7..8],
        .surface_humus_nitrogen_inhibition = diagnostic[8..9],
        .surface_organic_nitrogen_inhibition = diagnostic[9..10],
        .surface_organic_phosphorus_inhibition = diagnostic[10..11],
    }, .{
        .layer_kind = .surface_litter,
        .dimensions = d,
        .microbial = microbial,
        .microbial_residue = residue,
        .mobile = .{
            .micropore_dissolved = dissolved,
            .macropore_dissolved = zeros,
            .adsorbed = zeros,
            .micropore_acetate_g_c = &.{1},
            .macropore_acetate_g_c = &.{0},
            .adsorbed_acetate_g_c = &.{0},
        },
        .structural = structural,
    });
    try std.testing.expectEqual(@as(f64, 18), scalars[0]);
    try std.testing.expectEqual(@as(f64, 5), scalars[3]);
    try std.testing.expectEqual(@as(f64, 24), scalars[6]);
    try std.testing.expectEqual(@as(f64, 22), surface[0]);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** 11), &diagnostic);
}

test "dimension failure is atomic" {
    var scalar: f64 = 9;
    var surface = [_]f64{9};
    var diagnostic = [_]f64{9.0} ** 11;
    const inventory: InventoryState = .{
        .total_noncharcoal_carbon_g_c = &scalar,
        .total_noncharcoal_nitrogen_g_n = &scalar,
        .total_noncharcoal_phosphorus_g_p = &scalar,
        .charcoal_carbon_g_c = &scalar,
        .charcoal_nitrogen_g_n = &scalar,
        .charcoal_phosphorus_g_p = &scalar,
        .detrital_carbon_g_c = &scalar,
        .detrital_nitrogen_g_n = &scalar,
        .detrital_phosphorus_g_p = &scalar,
        .organic_carbon_snapshot_g_c = &scalar,
        .litter_carbon_g_c = &scalar,
        .charcoal_carbon_change_g_c = &scalar,
        .surface_carbon_g_c_by_microbial_pool = &surface,
    };
    const diag: TransformationDiagnostics = .{
        .oxygen_response = diagnostic[0..1],
        .maximum_rate_4 = diagnostic[1..2],
        .maximum_rate_3 = diagnostic[2..3],
        .maximum_rate_2 = diagnostic[3..4],
        .maximum_rate_1 = diagnostic[4..5],
        .humus_nitrogen_inhibition = diagnostic[5..6],
        .organic_nitrogen_inhibition = diagnostic[6..7],
        .organic_phosphorus_inhibition = diagnostic[7..8],
        .surface_humus_nitrogen_inhibition = diagnostic[8..9],
        .surface_organic_nitrogen_inhibition = diagnostic[9..10],
        .surface_organic_phosphorus_inhibition = diagnostic[10..11],
    };
    try std.testing.expectError(
        error.InitialOrganicInventoryDimensionMismatch,
        initialize(std.testing.allocator, inventory, diag, .{
            .layer_kind = .surface_litter,
            .dimensions = .{
                .organic_pool_count = 1,
                .microbial_pool_count = 1,
                .microbial_population_count = 1,
                .microbial_component_count = 1,
                .residue_component_count = 1,
                .structural_component_count = 1,
                .litter_pool_count = 1,
                .detrital_pool_count = 1,
                .noncharcoal_structural_component_count = 1,
            },
            .microbial = zeroPools(0),
            .microbial_residue = zeroPools(1),
            .mobile = .{
                .micropore_dissolved = zeroPools(1),
                .macropore_dissolved = zeroPools(1),
                .adsorbed = zeroPools(1),
                .micropore_acetate_g_c = &.{0},
                .macropore_acetate_g_c = &.{0},
                .adsorbed_acetate_g_c = &.{0},
            },
            .structural = zeroPools(1),
        }),
    );
    try std.testing.expectEqual(@as(f64, 9), scalar);
    try std.testing.expectEqual(@as(f64, 9), surface[0]);
}
