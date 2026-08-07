const std = @import("std");

pub const species_count = 49;
pub const Direction = enum { horizontal, vertical };

pub const Layer = struct {
    macropore_water_m3: f64,
    minimum_water_m3: f64,
    macropore_volume_m3: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    /// 33 salts/complexes (no H4SiO4), 8 non-band P, 8 band P.
    solute_inventory_amount: []const f64,
};

pub const Inputs = struct {
    current: Layer,
    adjacent: Layer,
    maximum_convective_fraction: f64,
    macropore_water_flux_m3_per_step: f64,
    direction: Direction,
    /// Current-cell macropore/micropore exchange, same 49-species order, amount step-1.
    current_exchange_amount_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 5431--5742.
/// Positive vertical flow conditionally subtracts only negative current-cell
/// exchange before convection, exactly retaining the legacy partition choices.
pub fn calculate(inputs: Inputs) ![species_count]f64 {
    try validateTopologyAndFlux(inputs);
    var flux = [_]f64{0} ** species_count;
    const water_flux = inputs.macropore_water_flux_m3_per_step;
    if (water_flux == 0) return flux;

    try validateLayer(inputs.current);
    try validateLayer(inputs.adjacent);

    const positive = water_flux > 0;
    const donor = if (positive) inputs.current else inputs.adjacent;
    const fraction = if (donor.macropore_water_m3 > donor.minimum_water_m3)
        if (positive)
            std.math.clamp(water_flux / donor.macropore_water_m3, 0, inputs.maximum_convective_fraction)
        else
            std.math.clamp(water_flux / donor.macropore_water_m3, -inputs.maximum_convective_fraction, 0)
    else if (positive)
        inputs.maximum_convective_fraction
    else
        -inputs.maximum_convective_fraction;

    const account_exchange = positive and inputs.direction == .vertical and
        inputs.adjacent.macropore_volume_m3 > inputs.adjacent.macropore_water_m3;
    if (account_exchange) for (inputs.current_exchange_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeConvectiveSoluteInput;
    for (0..species_count) |species| {
        const available = if (account_exchange)
            donor.solute_inventory_amount[species] - @min(0, inputs.current_exchange_amount_per_step[species])
        else
            donor.solute_inventory_amount[species];
        const partition_layer = if (account_exchange)
            inputs.current
        else if (positive)
            inputs.adjacent
        else
            inputs.adjacent;
        const partition = if (species >= 41)
            partition_layer.band_phosphate_fraction
        else if (species >= 33)
            partition_layer.nonband_phosphate_fraction
        else
            1;
        flux[species] = fraction * @max(0, available) * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteSoilBoundaryMacroporeConvectiveSoluteResult;
    }
    return flux;
}

fn validateTopologyAndFlux(inputs: Inputs) !void {
    if (inputs.current.solute_inventory_amount.len != species_count or
        inputs.adjacent.solute_inventory_amount.len != species_count or
        inputs.current_exchange_amount_per_step.len != species_count)
        return error.SoilBoundaryMacroporeConvectiveSoluteDimensionMismatch;
    inline for (.{ inputs.maximum_convective_fraction, inputs.macropore_water_flux_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeConvectiveSoluteInput;
    if (inputs.maximum_convective_fraction < 0) return error.InvalidSoilBoundaryMacroporeConvectiveSoluteInput;
}

fn validateLayer(layer: Layer) !void {
    inline for (.{ layer.macropore_water_m3, layer.minimum_water_m3, layer.macropore_volume_m3, layer.nonband_phosphate_fraction, layer.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeConvectiveSoluteInput;
    if (layer.macropore_water_m3 < 0 or layer.minimum_water_m3 < 0 or layer.macropore_volume_m3 < 0 or
        layer.nonband_phosphate_fraction < 0 or layer.nonband_phosphate_fraction > 1 or
        layer.band_phosphate_fraction < 0 or layer.band_phosphate_fraction > 1)
        return error.InvalidSoilBoundaryMacroporeConvectiveSoluteInput;
    for (layer.solute_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeConvectiveSoluteInput;
}

fn fixtureLayer(inventory: []const f64, nonband: f64, band: f64) Layer {
    return .{ .macropore_water_m3 = 4, .minimum_water_m3 = 0.01, .macropore_volume_m3 = 8, .nonband_phosphate_fraction = nonband, .band_phosphate_fraction = band, .solute_inventory_amount = inventory };
}

test "TRNSFRS positive nonvertical macropore flow uses current donor and adjacent P partition" {
    const current_inventory = [_]f64{8} ** species_count;
    const adjacent_inventory = [_]f64{20} ** species_count;
    const exchange = [_]f64{-4} ** species_count;
    const flux = try calculate(.{ .current = fixtureLayer(&current_inventory, 0.2, 0.3), .adjacent = fixtureLayer(&adjacent_inventory, 0.25, 0.5), .maximum_convective_fraction = 0.75, .macropore_water_flux_m3_per_step = 2, .direction = .horizontal, .current_exchange_amount_per_step = &exchange });
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 1), flux[33]);
    try std.testing.expectEqual(@as(f64, 2), flux[41]);
}

test "TRNSFRS positive unsaturated vertical flow includes negative exchange and current P partition" {
    const current_inventory = [_]f64{8} ** species_count;
    const adjacent_inventory = [_]f64{20} ** species_count;
    const exchange = [_]f64{-4} ** species_count;
    const flux = try calculate(.{ .current = fixtureLayer(&current_inventory, 0.2, 0.3), .adjacent = fixtureLayer(&adjacent_inventory, 0.25, 0.5), .maximum_convective_fraction = 0.75, .macropore_water_flux_m3_per_step = 2, .direction = .vertical, .current_exchange_amount_per_step = &exchange });
    try std.testing.expectEqual(@as(f64, 6), flux[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), flux[33], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.8), flux[41], 1e-15);
}

test "TRNSFRS negative flow uses adjacent inventory and adjacent partitions" {
    const current_inventory = [_]f64{99} ** species_count;
    const adjacent_inventory = [_]f64{8} ** species_count;
    const exchange = [_]f64{-99} ** species_count;
    const flux = try calculate(.{ .current = fixtureLayer(&current_inventory, 0.2, 0.3), .adjacent = fixtureLayer(&adjacent_inventory, 0.25, 0.5), .maximum_convective_fraction = 0.75, .macropore_water_flux_m3_per_step = -2, .direction = .vertical, .current_exchange_amount_per_step = &exchange });
    try std.testing.expectEqual(@as(f64, -4), flux[0]);
    try std.testing.expectEqual(@as(f64, -1), flux[33]);
    try std.testing.expectEqual(@as(f64, -2), flux[41]);
}

test "TRNSFRS zero macropore flow explicitly clears all species" {
    const inventory = [_]f64{8} ** species_count;
    const exchange = [_]f64{-4} ** species_count;
    const flux = try calculate(.{ .current = fixtureLayer(&inventory, 0.2, 0.3), .adjacent = fixtureLayer(&inventory, 0.25, 0.5), .maximum_convective_fraction = 0.75, .macropore_water_flux_m3_per_step = 0, .direction = .vertical, .current_exchange_amount_per_step = &exchange });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "TRNSFRS zero macropore flow does not evaluate dormant state" {
    var inventory = [_]f64{8} ** species_count;
    var exchange = [_]f64{-4} ** species_count;
    inventory[48] = std.math.nan(f64);
    exchange[48] = std.math.inf(f64);
    var current = fixtureLayer(&inventory, 0.2, 0.3);
    current.macropore_water_m3 = std.math.nan(f64);
    const flux = try calculate(.{ .current = current, .adjacent = fixtureLayer(&inventory, 0.25, 0.5), .maximum_convective_fraction = 0.75, .macropore_water_flux_m3_per_step = 0, .direction = .vertical, .current_exchange_amount_per_step = &exchange });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "late invalid exchange fails atomically" {
    const inventory = [_]f64{8} ** species_count;
    var exchange = [_]f64{0} ** species_count;
    exchange[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundaryMacroporeConvectiveSoluteInput, calculate(.{ .current = fixtureLayer(&inventory, 0.2, 0.3), .adjacent = fixtureLayer(&inventory, 0.25, 0.5), .maximum_convective_fraction = 0.75, .macropore_water_flux_m3_per_step = 1, .direction = .vertical, .current_exchange_amount_per_step = &exchange }));
}
