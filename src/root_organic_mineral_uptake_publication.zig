const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const MineralUptake = struct {
    ammonium_nonband_g_n: f64,
    ammonium_band_g_n: f64,
    nitrate_nonband_g_n: f64,
    nitrate_band_g_n: f64,
    dihydrogen_phosphate_nonband_g_p: f64,
    dihydrogen_phosphate_band_g_p: f64,
    hydrogen_phosphate_nonband_g_p: f64,
    hydrogen_phosphate_band_g_p: f64,
};

pub const RootLayerState = struct {
    mobile_carbon_g_c: []f64,
    mobile_nitrogen_g_n: []f64,
    mobile_phosphorus_g_p: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    organic_substrate_count: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    /// RDFOM*, flattened [domain][layer][substrate], signed positive uptake.
    organic_exchange: []const Pool,
    mineral_uptake_by_domain_layer: []const MineralUptake,
};

/// Exact GROSUB lines 5803--5822 root-side RDFOM C:N:P publication followed
/// by mineral NH4/NO3 and phosphate assimilation. Domains are outermost,
/// rooted layers next, organic substrate classes innermost. Runtime topology
/// replaces the source's compiled N/L/K extents.
pub fn apply(state: RootLayerState, inputs: Inputs) !void {
    const root_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.RootUptakePublicationDimensionOverflow;
    const exchange_count = std.math.mul(usize, root_layer_count, inputs.organic_substrate_count) catch
        return error.RootUptakePublicationDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or
        inputs.organic_substrate_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.organic_exchange.len != exchange_count or
        inputs.mineral_uptake_by_domain_layer.len != root_layer_count or
        state.mobile_carbon_g_c.len != root_layer_count or
        state.mobile_nitrogen_g_n.len != root_layer_count or
        state.mobile_phosphorus_g_p.len != root_layer_count)
        return error.RootUptakePublicationDimensionMismatch;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or
        inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidRootedLayerRange;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or
        inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidRootUptakePublicationInput;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidRootLayerThickness;

    for (0..inputs.biological_domain_count) |domain| {
        for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            _ = try nextPool(state, inputs, domain, layer);
        }
    }
    for (0..inputs.biological_domain_count) |domain| {
        for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const index = domain * inputs.soil_layer_count + layer;
            const next = try nextPool(state, inputs, domain, layer);
            state.mobile_carbon_g_c[index] = next.carbon_g_c;
            state.mobile_nitrogen_g_n[index] = next.nitrogen_g_n;
            state.mobile_phosphorus_g_p[index] = next.phosphorus_g_p;
        }
    }
}

fn nextPool(state: RootLayerState, inputs: Inputs, domain: usize, layer: usize) !Pool {
    const root_index = domain * inputs.soil_layer_count + layer;
    var next = Pool{
        .carbon_g_c = state.mobile_carbon_g_c[root_index],
        .nitrogen_g_n = state.mobile_nitrogen_g_n[root_index],
        .phosphorus_g_p = state.mobile_phosphorus_g_p[root_index],
    };
    try validateNonnegativePool(next);
    const exchange_base = root_index * inputs.organic_substrate_count;
    for (0..inputs.organic_substrate_count) |substrate| {
        const exchange = inputs.organic_exchange[exchange_base + substrate];
        try validateSignedPool(exchange);
        next.carbon_g_c += exchange.carbon_g_c;
        next.nitrogen_g_n += exchange.nitrogen_g_n;
        next.phosphorus_g_p += exchange.phosphorus_g_p;
        validateNonnegativePool(next) catch return error.RootOrganicExchangeWouldOverdraw;
    }

    const mineral = inputs.mineral_uptake_by_domain_layer[root_index];
    inline for (@typeInfo(MineralUptake).@"struct".fields) |field| {
        const value = @field(mineral, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootMineralUptake;
    }
    next.nitrogen_g_n = next.nitrogen_g_n + mineral.ammonium_nonband_g_n +
        mineral.ammonium_band_g_n + mineral.nitrate_nonband_g_n + mineral.nitrate_band_g_n;
    next.phosphorus_g_p = next.phosphorus_g_p + mineral.dihydrogen_phosphate_nonband_g_p +
        mineral.dihydrogen_phosphate_band_g_p + mineral.hydrogen_phosphate_nonband_g_p +
        mineral.hydrogen_phosphate_band_g_p;
    try validateNonnegativePool(next);
    return next;
}

fn validateNonnegativePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootMobilePool;
}

fn validateSignedPool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteRootOrganicExchange;
}

fn zeroMineral() MineralUptake {
    return std.mem.zeroes(MineralUptake);
}

test "GROSUB publishes signed organic classes before mineral N and P" {
    var carbon = [_]f64{10};
    var nitrogen = [_]f64{1};
    var phosphorus = [_]f64{0.1};
    const organic = [_]Pool{
        .{ .carbon_g_c = -1, .nitrogen_g_n = -0.1, .phosphorus_g_p = -0.01 },
        .{ .carbon_g_c = 0.5, .nitrogen_g_n = 0.05, .phosphorus_g_p = 0.005 },
        .{ .carbon_g_c = 0.4, .nitrogen_g_n = 0.04, .phosphorus_g_p = 0.004 },
        .{ .carbon_g_c = 0.3, .nitrogen_g_n = 0.03, .phosphorus_g_p = 0.003 },
        .{ .carbon_g_c = 0.2, .nitrogen_g_n = 0.02, .phosphorus_g_p = 0.002 },
    };
    var mineral = zeroMineral();
    mineral.ammonium_nonband_g_n = 0.2;
    mineral.nitrate_band_g_n = 0.3;
    mineral.dihydrogen_phosphate_band_g_p = 0.02;
    mineral.hydrogen_phosphate_nonband_g_p = 0.03;
    try apply(.{ .mobile_carbon_g_c = &carbon, .mobile_nitrogen_g_n = &nitrogen, .mobile_phosphorus_g_p = &phosphorus }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .organic_substrate_count = 5,
        .layer_thickness_m = &.{0.1},
        .minimum_active_layer_thickness_m = 0.001,
        .organic_exchange = &organic,
        .mineral_uptake_by_domain_layer = &.{mineral},
    });
    try std.testing.expectApproxEqAbs(@as(f64, 10.4), carbon[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1.54), nitrogen[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.154), phosphorus[0], 1e-14);
}

test "GROSUB root publication C N P delta equals organic plus mineral inputs" {
    var carbon = [_]f64{5};
    var nitrogen = [_]f64{2};
    var phosphorus = [_]f64{1};
    const organic = [_]Pool{.{ .carbon_g_c = 0.5, .nitrogen_g_n = -0.2, .phosphorus_g_p = 0.1 }};
    var mineral = zeroMineral();
    mineral.ammonium_band_g_n = 0.3;
    mineral.hydrogen_phosphate_band_g_p = 0.4;
    try apply(.{ .mobile_carbon_g_c = &carbon, .mobile_nitrogen_g_n = &nitrogen, .mobile_phosphorus_g_p = &phosphorus }, .{ .biological_domain_count = 1, .soil_layer_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .organic_substrate_count = 1, .layer_thickness_m = &.{1}, .minimum_active_layer_thickness_m = 0, .organic_exchange = &organic, .mineral_uptake_by_domain_layer = &.{mineral} });
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), carbon[0] - 5, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), nitrogen[0] - 2, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), phosphorus[0] - 1, 1e-15);
}

test "GROSUB runtime domain layer substrate sweep is atomic on late exudation overdraw" {
    const domains = 4;
    const layers = 9;
    const substrates = 7;
    var carbon = [_]f64{10} ** (domains * layers);
    var nitrogen = [_]f64{10} ** (domains * layers);
    var phosphorus = [_]f64{10} ** (domains * layers);
    var organic = [_]Pool{.{ .carbon_g_c = 0.1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.1 }} ** (domains * layers * substrates);
    const mineral = [_]MineralUptake{zeroMineral()} ** (domains * layers);
    const thickness = [_]f64{0.1} ** layers;
    const state = RootLayerState{ .mobile_carbon_g_c = &carbon, .mobile_nitrogen_g_n = &nitrogen, .mobile_phosphorus_g_p = &phosphorus };
    const inputs = Inputs{ .biological_domain_count = domains, .soil_layer_count = layers, .planting_layer_index = 1, .deepest_rooted_layer_index = 8, .organic_substrate_count = substrates, .layer_thickness_m = &thickness, .minimum_active_layer_thickness_m = 0.001, .organic_exchange = &organic, .mineral_uptake_by_domain_layer = &mineral };
    try apply(state, inputs);
    const before_carbon = carbon;
    const before_nitrogen = nitrogen;
    organic[((domains - 1) * layers + 8) * substrates + 6].carbon_g_c = -100;
    try std.testing.expectError(error.RootOrganicExchangeWouldOverdraw, apply(state, inputs));
    try std.testing.expectEqualDeep(before_carbon, carbon);
    try std.testing.expectEqualDeep(before_nitrogen, nitrogen);
}

test "GROSUB thin inactive layer does not read exchange or root pools" {
    var carbon = [_]f64{std.math.nan(f64)};
    var nitrogen = [_]f64{std.math.nan(f64)};
    var phosphorus = [_]f64{std.math.nan(f64)};
    const organic = [_]Pool{.{ .carbon_g_c = std.math.nan(f64), .nitrogen_g_n = 0, .phosphorus_g_p = 0 }};
    try apply(.{ .mobile_carbon_g_c = &carbon, .mobile_nitrogen_g_n = &nitrogen, .mobile_phosphorus_g_p = &phosphorus }, .{ .biological_domain_count = 1, .soil_layer_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .organic_substrate_count = 1, .layer_thickness_m = &.{0.001}, .minimum_active_layer_thickness_m = 0.001, .organic_exchange = &organic, .mineral_uptake_by_domain_layer = &.{zeroMineral()} });
    try std.testing.expect(std.math.isNan(carbon[0]));
}
