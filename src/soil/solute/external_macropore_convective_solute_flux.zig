const std = @import("std");

pub const species_count = 49;
pub const Side = enum { forward, reverse };

pub const Inputs = struct {
    side: Side,
    macropore_water_flux_m3_per_step: f64,
    micropore_water_flux_m3_per_step: f64,
    source_macropore_water_m3: f64,
    minimum_water_m3: f64,
    maximum_convective_fraction: f64,
    minimum_convective_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    /// Compact topology: 33 salts/complexes without H4SiO4, then 16 P fields.
    source_macropore_inventory_amount: []const f64,
    prescribed_recharge_concentration_amount_per_m3: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 8082--8415.
/// Preserves asymmetric legacy branches: discharge exists only for forward
/// positive FLWHM; recharge selection is gated by FLWM sign, while candidate
/// fluxes use FLWHM and the clamped FLWHM/VOLWHM fraction.
pub fn calculate(inputs: Inputs) ![species_count]f64 {
    try validateHydraulicControls(inputs);
    var flux = [_]f64{0} ** species_count;
    const fraction = if (inputs.source_macropore_water_m3 > inputs.minimum_water_m3)
        std.math.clamp(inputs.macropore_water_flux_m3_per_step / inputs.source_macropore_water_m3, -inputs.maximum_convective_fraction, inputs.maximum_convective_fraction)
    else
        0;
    if (@abs(fraction) <= inputs.minimum_convective_fraction) return flux;

    const discharge = inputs.side == .forward and inputs.macropore_water_flux_m3_per_step > 0;
    const forward_recharge = inputs.side == .forward and inputs.micropore_water_flux_m3_per_step < 0;
    const reverse_recharge = inputs.side == .reverse and inputs.micropore_water_flux_m3_per_step > 0;
    if (!discharge and !forward_recharge and !reverse_recharge) return flux;

    try validateActiveChemistry(inputs, discharge);

    for (0..species_count) |species| {
        const inventory_flux = fraction * @max(0, inputs.source_macropore_inventory_amount[species]);
        const recharge_flux = inputs.macropore_water_flux_m3_per_step *
            inputs.prescribed_recharge_concentration_amount_per_m3[species];
        const unpartitioned = if (discharge)
            inventory_flux
        else if (forward_recharge)
            @max(inventory_flux, recharge_flux)
        else
            @min(inventory_flux, recharge_flux);
        const partition = if (species >= 41)
            inputs.band_phosphate_fraction
        else if (species >= 33)
            inputs.nonband_phosphate_fraction
        else
            1;
        flux[species] = unpartitioned * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteExternalMacroporeConvectiveSoluteResult;
    }
    return flux;
}

fn validateHydraulicControls(inputs: Inputs) !void {
    inline for (.{ inputs.macropore_water_flux_m3_per_step, inputs.micropore_water_flux_m3_per_step, inputs.source_macropore_water_m3, inputs.minimum_water_m3, inputs.maximum_convective_fraction, inputs.minimum_convective_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalMacroporeConvectiveSoluteInput;
    if (inputs.source_macropore_water_m3 < 0 or inputs.minimum_water_m3 < 0 or
        inputs.maximum_convective_fraction < 0 or inputs.minimum_convective_fraction < 0)
        return error.InvalidExternalMacroporeConvectiveSoluteInput;
}

fn validateActiveChemistry(inputs: Inputs, discharge: bool) !void {
    inline for (.{ inputs.nonband_phosphate_fraction, inputs.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidExternalMacroporeConvectiveSoluteInput;
    const active_values = if (discharge)
        inputs.source_macropore_inventory_amount
    else
        inputs.prescribed_recharge_concentration_amount_per_m3;
    if (active_values.len != species_count)
        return error.ExternalMacroporeConvectiveSoluteDimensionMismatch;
    for (active_values) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalMacroporeConvectiveSoluteInput;
}

fn fixture(inventory: []const f64, concentration: []const f64, side: Side, macro_flux: f64, micro_flux: f64) Inputs {
    return .{ .side = side, .macropore_water_flux_m3_per_step = macro_flux, .micropore_water_flux_m3_per_step = micro_flux, .source_macropore_water_m3 = 4, .minimum_water_m3 = 0.01, .maximum_convective_fraction = 0.75, .minimum_convective_fraction = 1e-9, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .source_macropore_inventory_amount = inventory, .prescribed_recharge_concentration_amount_per_m3 = concentration };
}

test "TRNSFRS forward positive macropore discharge uses compact inventory" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{99} ** species_count;
    const flux = try calculate(fixture(&inventory, &concentration, .forward, 2, 0));
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 1), flux[33]);
    try std.testing.expectEqual(@as(f64, 2), flux[41]);
}

test "TRNSFRS forward recharge uses FLWM gate and maximum negative candidate" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    const flux = try calculate(fixture(&inventory, &concentration, .forward, -2, -1));
    try std.testing.expectEqual(@as(f64, -4), flux[0]);
    try std.testing.expectEqual(@as(f64, -1), flux[33]);
    try std.testing.expectEqual(@as(f64, -2), flux[41]);
}

test "TRNSFRS reverse recharge uses FLWM gate and minimum positive candidate" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    const flux = try calculate(fixture(&inventory, &concentration, .reverse, 2, 1));
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 1), flux[33]);
    try std.testing.expectEqual(@as(f64, 2), flux[41]);
}

test "TRNSFRS reverse negative macropore flow has no discharge branch" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    const flux = try calculate(fixture(&inventory, &concentration, .reverse, -2, -1));
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "late invalid prescribed concentration fails atomically" {
    const inventory = [_]f64{8} ** species_count;
    var concentration = [_]f64{3} ** species_count;
    concentration[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteExternalMacroporeConvectiveSoluteInput, calculate(fixture(&inventory, &concentration, .forward, -2, -1)));
}

test "TRNSFRS discharge does not inspect dormant recharge chemistry" {
    const inventory = [_]f64{8} ** species_count;
    var concentration = [_]f64{3} ** species_count;
    concentration[48] = std.math.inf(f64);
    const flux = try calculate(fixture(&inventory, &concentration, .forward, 2, 0));
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
}

test "TRNSFRS zero fraction returns before chemistry validation" {
    var inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    inventory[48] = std.math.inf(f64);
    var inputs = fixture(&inventory, &concentration, .forward, 2, 0);
    inputs.source_macropore_water_m3 = inputs.minimum_water_m3;
    inputs.nonband_phosphate_fraction = std.math.nan(f64);
    const flux = try calculate(inputs);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}
