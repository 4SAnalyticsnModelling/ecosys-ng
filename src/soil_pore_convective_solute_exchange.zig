const std = @import("std");

pub const micropore_species_count = 50;
pub const exchange_species_count = 49;
pub const Direction = enum { horizontal, vertical };

pub const Inputs = struct {
    direction: Direction,
    /// Positive is macropore to micropore, negative is micropore to macropore, m3 step-1.
    macropore_to_micropore_water_flux_m3_per_step: f64,
    micropore_water_m3: f64,
    macropore_water_m3: f64,
    minimum_water_m3: f64,
    maximum_convective_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    micropore_inventory_amount: []const f64,
    /// Compact topology excludes micropore H4SiO4.
    macropore_inventory_amount: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 6311--6542.
/// TRNSFRS lines 6278--6305 are an intentionally omitted diagnostic IF with
/// only commented WRITE statements. The legacy zero-flow branch fails to set
/// RFLM1B and instead sets unused RFLHYS; this safe translation clears all 49
/// published exchanges rather than propagating stale MgHPO4 band flux.
pub fn calculate(inputs: Inputs) !?[exchange_species_count]f64 {
    if (inputs.direction != .vertical) return null;
    try validateTopologyAndFlux(inputs);
    var flux = [_]f64{0} ** exchange_species_count;
    const water_flux = inputs.macropore_to_micropore_water_flux_m3_per_step;
    if (water_flux == 0) return flux;
    try validateActiveInputs(inputs);

    const macropore_donor = water_flux > 0;
    const donor_water_m3 = if (macropore_donor) inputs.macropore_water_m3 else inputs.micropore_water_m3;
    const fraction = if (donor_water_m3 > inputs.minimum_water_m3)
        if (macropore_donor)
            std.math.clamp(water_flux / donor_water_m3, 0, inputs.maximum_convective_fraction)
        else
            std.math.clamp(water_flux / donor_water_m3, -inputs.maximum_convective_fraction, 0)
    else if (macropore_donor)
        inputs.maximum_convective_fraction
    else
        -inputs.maximum_convective_fraction;

    for (0..exchange_species_count) |species| {
        const inventory = if (macropore_donor)
            inputs.macropore_inventory_amount[species]
        else
            inputs.micropore_inventory_amount[if (species < 33) species else species + 1];
        const partition = if (species >= 41)
            inputs.band_phosphate_fraction
        else if (species >= 33)
            inputs.nonband_phosphate_fraction
        else
            1;
        flux[species] = fraction * @max(0, inventory) * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteSoilPoreConvectiveSoluteExchangeResult;
    }
    return flux;
}

fn validateTopologyAndFlux(inputs: Inputs) !void {
    if (inputs.micropore_inventory_amount.len != micropore_species_count or
        inputs.macropore_inventory_amount.len != exchange_species_count)
        return error.SoilPoreConvectiveSoluteExchangeDimensionMismatch;
    if (!std.math.isFinite(inputs.macropore_to_micropore_water_flux_m3_per_step))
        return error.NonFiniteSoilPoreConvectiveSoluteExchangeInput;
}

fn validateActiveInputs(inputs: Inputs) !void {
    inline for (.{ inputs.micropore_water_m3, inputs.macropore_water_m3, inputs.minimum_water_m3, inputs.maximum_convective_fraction, inputs.nonband_phosphate_fraction, inputs.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreConvectiveSoluteExchangeInput;
    if (inputs.micropore_water_m3 < 0 or inputs.macropore_water_m3 < 0 or inputs.minimum_water_m3 < 0 or
        inputs.maximum_convective_fraction < 0 or inputs.nonband_phosphate_fraction < 0 or
        inputs.nonband_phosphate_fraction > 1 or inputs.band_phosphate_fraction < 0 or inputs.band_phosphate_fraction > 1)
        return error.InvalidSoilPoreConvectiveSoluteExchangeInput;
    for (inputs.micropore_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreConvectiveSoluteExchangeInput;
    for (inputs.macropore_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilPoreConvectiveSoluteExchangeInput;
}

fn fixture(micro: []const f64, macro: []const f64, water_flux: f64) Inputs {
    return .{ .direction = .vertical, .macropore_to_micropore_water_flux_m3_per_step = water_flux, .micropore_water_m3 = 4, .macropore_water_m3 = 2, .minimum_water_m3 = 0.01, .maximum_convective_fraction = 0.75, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .micropore_inventory_amount = micro, .macropore_inventory_amount = macro };
}

test "TRNSFRS positive pore exchange uses compact macropore donor and P fractions" {
    const micro = [_]f64{99} ** micropore_species_count;
    const macro = [_]f64{8} ** exchange_species_count;
    const flux = (try calculate(fixture(&micro, &macro, 1))).?;
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 1), flux[33]);
    try std.testing.expectEqual(@as(f64, 2), flux[41]);
}

test "TRNSFRS negative pore exchange skips micropore H4SiO4" {
    var micro = [_]f64{8} ** micropore_species_count;
    micro[33] = 1000;
    micro[34] = 12;
    const macro = [_]f64{99} ** exchange_species_count;
    const flux = (try calculate(fixture(&micro, &macro, -2))).?;
    try std.testing.expectEqual(@as(f64, -4), flux[0]);
    try std.testing.expectEqual(@as(f64, -1.5), flux[33]);
    try std.testing.expectEqual(@as(f64, -2), flux[41]);
}

test "TRNSFRS horizontal direction does not execute vertical exchange" {
    const micro = [_]f64{8} ** micropore_species_count;
    const macro = [_]f64{8} ** exchange_species_count;
    var inputs = fixture(&micro, &macro, 1);
    inputs.direction = .horizontal;
    try std.testing.expect((try calculate(inputs)) == null);
}

test "TRNSFRS horizontal direction leaves dormant invalid inputs unexamined" {
    const empty = [_]f64{};
    const inputs: Inputs = .{ .direction = .horizontal, .macropore_to_micropore_water_flux_m3_per_step = std.math.nan(f64), .micropore_water_m3 = std.math.nan(f64), .macropore_water_m3 = std.math.nan(f64), .minimum_water_m3 = -1, .maximum_convective_fraction = -1, .nonband_phosphate_fraction = -1, .band_phosphate_fraction = -1, .micropore_inventory_amount = &empty, .macropore_inventory_amount = &empty };
    try std.testing.expect((try calculate(inputs)) == null);
}

test "zero water exchange safely clears legacy stale final band species" {
    const micro = [_]f64{8} ** micropore_species_count;
    const macro = [_]f64{8} ** exchange_species_count;
    const flux = (try calculate(fixture(&micro, &macro, 0))).?;
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** exchange_species_count), &flux);
}

test "zero water exchange does not evaluate dormant inventories" {
    var micro = [_]f64{8} ** micropore_species_count;
    var macro = [_]f64{8} ** exchange_species_count;
    micro[49] = std.math.nan(f64);
    macro[48] = std.math.inf(f64);
    var inputs = fixture(&micro, &macro, 0);
    inputs.micropore_water_m3 = std.math.nan(f64);
    const flux = (try calculate(inputs)).?;
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** exchange_species_count), &flux);
}

test "late invalid inventory fails before publishing a result" {
    const micro = [_]f64{8} ** micropore_species_count;
    var macro = [_]f64{8} ** exchange_species_count;
    macro[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilPoreConvectiveSoluteExchangeInput, calculate(fixture(&micro, &macro, 1)));
}
