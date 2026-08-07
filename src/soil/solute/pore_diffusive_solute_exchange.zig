const std = @import("std");

pub const micropore_species_count = 50;
pub const exchange_species_count = 49;

pub const Inputs = struct {
    macropore_water_m3: f64,
    micropore_water_m3: f64,
    total_pore_volume_m3: f64,
    maximum_exchange_volume_fraction: f64,
    minimum_macropore_water_m3: f64,
    solute_exchange_timestep_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    micropore_inventory_amount: []const f64,
    /// Compact topology excludes micropore H4SiO4.
    macropore_inventory_amount: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 6570--6739.
/// The duplicate identical DFVNAC assignments at lines 6633--6640 are
/// represented once because the second assignment has no observable effect.
pub fn calculate(inputs: Inputs) ![exchange_species_count]f64 {
    try validateTopologyAndGateInputs(inputs);
    var flux = [_]f64{0} ** exchange_species_count;
    if (inputs.macropore_water_m3 <= inputs.minimum_macropore_water_m3) return flux;
    try validateActiveInputs(inputs);

    const exchanging_macropore_water_m3 = @min(
        inputs.maximum_exchange_volume_fraction * inputs.total_pore_volume_m3,
        inputs.macropore_water_m3,
    );
    const combined_water_m3 = inputs.micropore_water_m3 + exchanging_macropore_water_m3;
    if (!(combined_water_m3 > 0) or !std.math.isFinite(combined_water_m3))
        return error.InvalidSoilPoreDiffusiveSoluteExchangeWaterVolume;

    for (0..exchange_species_count) |species| {
        const micropore_species = if (species < 33) species else species + 1;
        const macro_inventory = @max(0, inputs.macropore_inventory_amount[species]);
        const micro_inventory = @max(0, inputs.micropore_inventory_amount[micropore_species]);
        const partition = if (species >= 41)
            inputs.band_phosphate_fraction
        else if (species >= 33)
            inputs.nonband_phosphate_fraction
        else
            1;
        flux[species] = inputs.solute_exchange_timestep_fraction *
            (macro_inventory * inputs.micropore_water_m3 - micro_inventory * exchanging_macropore_water_m3) /
            combined_water_m3 * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteSoilPoreDiffusiveSoluteExchangeResult;
    }
    return flux;
}

fn validateTopologyAndGateInputs(inputs: Inputs) !void {
    if (inputs.micropore_inventory_amount.len != micropore_species_count or
        inputs.macropore_inventory_amount.len != exchange_species_count)
        return error.SoilPoreDiffusiveSoluteExchangeDimensionMismatch;
    inline for (.{ inputs.macropore_water_m3, inputs.minimum_macropore_water_m3 }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreDiffusiveSoluteExchangeInput;
    if (inputs.macropore_water_m3 < 0 or inputs.minimum_macropore_water_m3 < 0)
        return error.InvalidSoilPoreDiffusiveSoluteExchangeInput;
}

fn validateActiveInputs(inputs: Inputs) !void {
    inline for (.{ inputs.micropore_water_m3, inputs.total_pore_volume_m3, inputs.maximum_exchange_volume_fraction, inputs.solute_exchange_timestep_fraction, inputs.nonband_phosphate_fraction, inputs.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreDiffusiveSoluteExchangeInput;
    if (inputs.micropore_water_m3 < 0 or inputs.total_pore_volume_m3 < 0 or
        inputs.maximum_exchange_volume_fraction < 0 or
        inputs.solute_exchange_timestep_fraction < 0 or inputs.nonband_phosphate_fraction < 0 or
        inputs.nonband_phosphate_fraction > 1 or inputs.band_phosphate_fraction < 0 or inputs.band_phosphate_fraction > 1)
        return error.InvalidSoilPoreDiffusiveSoluteExchangeInput;
    for (inputs.micropore_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreDiffusiveSoluteExchangeInput;
    for (inputs.macropore_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreDiffusiveSoluteExchangeInput;
}

fn fixture(micro: []const f64, macro: []const f64) Inputs {
    return .{ .macropore_water_m3 = 2, .micropore_water_m3 = 4, .total_pore_volume_m3 = 10, .maximum_exchange_volume_fraction = 0.1, .minimum_macropore_water_m3 = 0.01, .solute_exchange_timestep_fraction = 0.5, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .micropore_inventory_amount = micro, .macropore_inventory_amount = macro };
}

test "TRNSFRS pore diffusion uses capped macropore exchange water" {
    const micro = [_]f64{2} ** micropore_species_count;
    const macro = [_]f64{4} ** exchange_species_count;
    const flux = try calculate(fixture(&micro, &macro));
    try std.testing.expectEqual(@as(f64, 1.4), flux[0]);
    try std.testing.expectEqual(@as(f64, 0.35), flux[33]);
    try std.testing.expectEqual(@as(f64, 0.7), flux[41]);
}

test "TRNSFRS pore diffusion skips micropore H4SiO4 in compact mapping" {
    var micro = [_]f64{0} ** micropore_species_count;
    const macro = [_]f64{0} ** exchange_species_count;
    micro[33] = 1000;
    micro[34] = 10;
    const flux = try calculate(fixture(&micro, &macro));
    try std.testing.expectEqual(@as(f64, -0.25), flux[33]);
}

test "TRNSFRS equal weighted pore inventories have zero exchange" {
    const micro = [_]f64{8} ** micropore_species_count;
    const macro = [_]f64{2} ** exchange_species_count;
    const flux = try calculate(fixture(&micro, &macro));
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** exchange_species_count), &flux);
}

test "TRNSFRS dry macropore clears all diffusive exchange" {
    const micro = [_]f64{2} ** micropore_species_count;
    const macro = [_]f64{4} ** exchange_species_count;
    var inputs = fixture(&micro, &macro);
    inputs.macropore_water_m3 = inputs.minimum_macropore_water_m3;
    const flux = try calculate(inputs);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** exchange_species_count), &flux);
}

test "TRNSFRS dry macropore does not evaluate dormant exchange state" {
    var micro = [_]f64{2} ** micropore_species_count;
    var macro = [_]f64{4} ** exchange_species_count;
    micro[49] = std.math.nan(f64);
    macro[48] = std.math.inf(f64);
    var inputs = fixture(&micro, &macro);
    inputs.macropore_water_m3 = inputs.minimum_macropore_water_m3;
    inputs.micropore_water_m3 = std.math.nan(f64);
    inputs.solute_exchange_timestep_fraction = std.math.nan(f64);
    const flux = try calculate(inputs);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** exchange_species_count), &flux);
}

test "late invalid compact inventory fails before publishing a result" {
    const micro = [_]f64{2} ** micropore_species_count;
    var macro = [_]f64{4} ** exchange_species_count;
    macro[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilPoreDiffusiveSoluteExchangeInput, calculate(fixture(&micro, &macro)));
}
