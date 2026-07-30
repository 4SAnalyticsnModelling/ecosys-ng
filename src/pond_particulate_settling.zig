const std = @import("std");

/// Unit carried by one extensive particulate inventory.
///
/// REDIST settles pools with different units using the same dimensionless
/// fraction. Values are never summed across unlike units.
pub const ExtensiveUnit = enum {
    megagrams,
    moles,
    grams_carbon,
    grams_nitrogen,
    grams_phosphorus,
};

pub const DiagnosticRole = enum {
    particulate_inventory,
    mineral_sediment_mass,
};

/// One runtime-layer extensive inventory. The caller supplies pools in the
/// legacy REDIST calculation order when intermediate diagnostics depend on it.
pub const Pool = struct {
    name: []const u8,
    unit: ExtensiveUnit,
    diagnostic_role: DiagnosticRole = .particulate_inventory,
    amount_by_layer: []f64,
};

pub const Geometry = struct {
    /// Bulk density by runtime layer (Mg m-3); zero denotes pond water.
    bulk_density_Mg_per_m3: []const f64,
    /// Current layer thickness by runtime layer (m).
    layer_thickness_m: []const f64,
    /// Fortran NL, converted to a zero-based inclusive layer index.
    last_layer: usize,
    /// Fortran NU, converted to a zero-based soil-surface layer index.
    surface_soil_layer: usize,
    /// Fortran DLYRM: minimum thickness for a receiving layer (m).
    minimum_receiver_thickness_m: f64,
};

pub const Options = struct {
    /// Fortran XNFH: current science-step duration (h).
    timestep_h: f64,
    /// Active compatibility formulation from REDIST line 352 (h-1).
    settling_rate_per_h: f64 = 0.001,
};

pub const Result = struct {
    /// Mineral sediment newly delivered into positive-density layers (Mg).
    deposited_mineral_sediment_Mg: f64,
    /// Deepest positive-density receiver that accepted sediment, if any.
    deepest_deposition_layer: ?usize,
};

/// Returns the REDIST compatibility settling fraction after applying the same
/// finite-value and physical-domain checks used by `settle`.
pub fn settlingFraction(options: Options) !f64 {
    if (!std.math.isFinite(options.timestep_h) or options.timestep_h < 0 or
        !std.math.isFinite(options.settling_rate_per_h) or
        options.settling_rate_per_h < 0)
        return error.InvalidPondSettlingControl;
    const fraction = options.settling_rate_per_h * options.timestep_h;
    if (!std.math.isFinite(fraction) or fraction > 1)
        return error.InvalidPondSettlingFraction;
    return fraction;
}

/// Settles pond particulates in the exact source traversal order.
///
/// Traceability: REDIST (`redist.f`) lines 333-614. Donors are visited from
/// `NL-1` through layer zero. The first deeper layer thicker than DLYRM is the
/// receiver. Descending traversal is scientifically significant: material
/// received by a deeper layer is not settled again during the same call.
pub fn settle(
    geometry: Geometry,
    pools: []const Pool,
    options: Options,
) !Result {
    try validate(geometry, pools, options);
    const settling_fraction = try settlingFraction(options);

    var result: Result = .{
        .deposited_mineral_sediment_Mg = 0,
        .deepest_deposition_layer = null,
    };
    if (geometry.last_layer == 0 or settling_fraction == 0) return result;

    var donor = geometry.last_layer;
    while (donor > 0) {
        donor -= 1;
        if (!isPondDonor(geometry, donor)) continue;
        const receiver = firstReceiver(geometry, donor) orelse
            return error.MissingPondSettlingReceiver;

        for (pools) |pool| {
            const transfer = settling_fraction * pool.amount_by_layer[donor];
            pool.amount_by_layer[donor] -= transfer;
            pool.amount_by_layer[receiver] += transfer;
            if (pool.diagnostic_role == .mineral_sediment_mass and
                geometry.bulk_density_Mg_per_m3[receiver] > 0)
            {
                result.deposited_mineral_sediment_Mg += transfer;
                result.deepest_deposition_layer = receiver;
            }
        }
    }
    return result;
}

fn isPondDonor(geometry: Geometry, layer: usize) bool {
    const water_layer = geometry.bulk_density_Mg_per_m3[layer] <= 0;
    const surface_over_pond =
        layer == 0 and
        geometry.bulk_density_Mg_per_m3[geometry.surface_soil_layer] <= 0;
    return (water_layer or surface_over_pond) and
        geometry.layer_thickness_m[layer] > 0;
}

fn firstReceiver(geometry: Geometry, donor: usize) ?usize {
    var receiver = donor + 1;
    while (receiver <= geometry.last_layer) : (receiver += 1) {
        if (geometry.layer_thickness_m[receiver] >
            geometry.minimum_receiver_thickness_m) return receiver;
    }
    return null;
}

fn validate(geometry: Geometry, pools: []const Pool, options: Options) !void {
    const layer_count = geometry.bulk_density_Mg_per_m3.len;
    if (layer_count == 0 or geometry.layer_thickness_m.len != layer_count or
        geometry.last_layer >= layer_count or
        geometry.surface_soil_layer >= layer_count)
        return error.PondSettlingGeometryDimensionMismatch;
    if (!std.math.isFinite(geometry.minimum_receiver_thickness_m) or
        geometry.minimum_receiver_thickness_m < 0)
        return error.InvalidPondSettlingControl;
    _ = try settlingFraction(options);

    for (geometry.bulk_density_Mg_per_m3, geometry.layer_thickness_m) |
        density,
        thickness,
    | {
        if (!std.math.isFinite(density) or density < 0 or
            !std.math.isFinite(thickness) or thickness < 0)
            return error.InvalidPondSettlingGeometry;
    }
    if (geometry.last_layer > 0) {
        var donor = geometry.last_layer;
        while (donor > 0) {
            donor -= 1;
            if (isPondDonor(geometry, donor) and
                firstReceiver(geometry, donor) == null)
                return error.MissingPondSettlingReceiver;
        }
    }
    for (pools) |pool| {
        if (pool.name.len == 0 or pool.amount_by_layer.len != layer_count)
            return error.PondSettlingPoolDimensionMismatch;
        var total: f64 = 0;
        for (pool.amount_by_layer) |amount| {
            if (!std.math.isFinite(amount) or amount < 0)
                return error.InvalidPondParticulateInventory;
            total += amount;
            if (!std.math.isFinite(total))
                return error.PondParticulateInventoryOverflow;
        }
    }
}

test "REDIST descending order prevents same-step settling cascade" {
    var sediment = [_]f64{ 100, 50, 0 };
    var carbon = [_]f64{ 20, 10, 0 };
    const pools = [_]Pool{
        .{
            .name = "sand",
            .unit = .megagrams,
            .diagnostic_role = .mineral_sediment_mass,
            .amount_by_layer = &sediment,
        },
        .{
            .name = "microbial carbon",
            .unit = .grams_carbon,
            .amount_by_layer = &carbon,
        },
    };

    const result = try settle(.{
        .bulk_density_Mg_per_m3 = &.{ 0, 0, 1.2 },
        .layer_thickness_m = &.{ 0.1, 0.1, 0.2 },
        .last_layer = 2,
        .surface_soil_layer = 2,
        .minimum_receiver_thickness_m = 1e-6,
    }, &pools, .{ .timestep_h = 1 });

    // Layer 1 settles first. Layer 0 then settles into layer 1, after its turn.
    try std.testing.expectApproxEqAbs(@as(f64, 99.9), sediment[0], 1e-13);
    try std.testing.expectApproxEqAbs(@as(f64, 50.05), sediment[1], 1e-13);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), sediment[2], 1e-13);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), result.deposited_mineral_sediment_Mg, 1e-13);
    try std.testing.expectEqual(@as(?usize, 2), result.deepest_deposition_layer);
}

test "mixed-unit particulate pools conserve each extensive inventory" {
    var mineral = [_]f64{ 2, 0 };
    var phosphorus = [_]f64{ 31, 0 };
    var exchange_sites = [_]f64{ 4, 0 };
    const pools = [_]Pool{
        .{ .name = "clay", .unit = .megagrams, .diagnostic_role = .mineral_sediment_mass, .amount_by_layer = &mineral },
        .{ .name = "organic phosphorus", .unit = .grams_phosphorus, .amount_by_layer = &phosphorus },
        .{ .name = "cation exchange capacity", .unit = .moles, .amount_by_layer = &exchange_sites },
    };

    _ = try settle(.{
        .bulk_density_Mg_per_m3 = &.{ 0, 1.1 },
        .layer_thickness_m = &.{ 0.2, 0.2 },
        .last_layer = 1,
        .surface_soil_layer = 1,
        .minimum_receiver_thickness_m = 0,
    }, &pools, .{ .timestep_h = 2 });

    try std.testing.expectApproxEqAbs(@as(f64, 2), mineral[0] + mineral[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 31), phosphorus[0] + phosphorus[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 4), exchange_sites[0] + exchange_sites[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.004), mineral[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.062), phosphorus[1], 1e-14);
}

test "zero-thickness layers are skipped when selecting receiver" {
    var nitrogen = [_]f64{ 10, 0, 0 };
    const pools = [_]Pool{
        .{ .name = "microbial nitrogen", .unit = .grams_nitrogen, .amount_by_layer = &nitrogen },
    };
    _ = try settle(.{
        .bulk_density_Mg_per_m3 = &.{ 0, 0, 1 },
        .layer_thickness_m = &.{ 0.1, 1e-7, 0.2 },
        .last_layer = 2,
        .surface_soil_layer = 2,
        .minimum_receiver_thickness_m = 1e-6,
    }, &pools, .{ .timestep_h = 1 });
    try std.testing.expectEqual(@as(f64, 0), nitrogen[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), nitrogen[2], 1e-15);
}

test "invalid late pool leaves all pools unchanged" {
    var valid = [_]f64{ 8, 0 };
    var invalid = [_]f64{ 2, -1 };
    const before = valid;
    const pools = [_]Pool{
        .{ .name = "valid", .unit = .moles, .amount_by_layer = &valid },
        .{ .name = "invalid", .unit = .moles, .amount_by_layer = &invalid },
    };
    try std.testing.expectError(error.InvalidPondParticulateInventory, settle(.{
        .bulk_density_Mg_per_m3 = &.{ 0, 1 },
        .layer_thickness_m = &.{ 0.1, 0.1 },
        .last_layer = 1,
        .surface_soil_layer = 1,
        .minimum_receiver_thickness_m = 0,
    }, &pools, .{ .timestep_h = 1 }));
    try std.testing.expectEqualSlices(f64, &before, &valid);
}

test "missing receiver is rejected before any pool mutation" {
    var inventory = [_]f64{ 5, 0 };
    const before = inventory;
    const pools = [_]Pool{
        .{ .name = "particulate", .unit = .moles, .amount_by_layer = &inventory },
    };
    try std.testing.expectError(error.MissingPondSettlingReceiver, settle(.{
        .bulk_density_Mg_per_m3 = &.{ 0, 0 },
        .layer_thickness_m = &.{ 0.1, 0 },
        .last_layer = 1,
        .surface_soil_layer = 1,
        .minimum_receiver_thickness_m = 0,
    }, &pools, .{ .timestep_h = 1 }));
    try std.testing.expectEqualSlices(f64, &before, &inventory);
}
