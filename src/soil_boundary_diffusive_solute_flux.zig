const std = @import("std");

pub const species_count = 50;
pub const diffusivity_group_count = 14;

pub const Layer = struct {
    volumetric_water_content: f64,
    minimum_volumetric_water_content: f64,
    micropore_water_m3: f64,
    minimum_water_m3: f64,
    thickness_m: f64,
    tortuosity: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    band_nitrate_fraction: f64,
    solute_inventory_amount: []const f64,
};

pub const Inputs = struct {
    current: Layer,
    adjacent: Layer,
    minimum_length_m: f64,
    maximum_flux_velocity_m_per_step: f64,
    water_flux_m3_per_step: f64,
    interface_area_m2: f64,
    dispersivity: f64,
    /// PO4, Al, Fe, H, Ca, Mg, Na, K, OH, SO4, Cl, CO3, HCO3, H4SiO4; m2 step-1.
    aqueous_diffusivity_m2_per_step: [diffusivity_group_count]f64,
};

/// Compatibility translation of TRNSFRS.F lines 5074--5425.
/// Preserves the legacy adjacent-P denominators and band-P VLNOB weighting.
pub fn calculate(inputs: Inputs) ![species_count]f64 {
    try validateTopologyAndGateInputs(inputs);
    var flux = [_]f64{0} ** species_count;
    if (inputs.current.volumetric_water_content <= inputs.current.minimum_volumetric_water_content or
        inputs.adjacent.volumetric_water_content <= inputs.adjacent.minimum_volumetric_water_content or
        inputs.current.micropore_water_m3 <= inputs.current.minimum_water_m3 or
        inputs.adjacent.micropore_water_m3 <= inputs.adjacent.minimum_water_m3)
        return flux;
    try validateActiveInputs(inputs);

    const current_concentration = concentrations(inputs.current, inputs.current.minimum_water_m3, false);
    const adjacent_concentration = concentrations(inputs.adjacent, inputs.current.minimum_water_m3, true);
    const current_thickness = @max(inputs.minimum_length_m, inputs.current.thickness_m);
    const adjacent_thickness = @max(inputs.minimum_length_m, inputs.adjacent.thickness_m);
    const tortuosity_per_m = (inputs.current.tortuosity + inputs.adjacent.tortuosity) /
        (current_thickness + adjacent_thickness);
    const dispersion_m_per_step = inputs.dispersivity *
        @min(inputs.maximum_flux_velocity_m_per_step, @abs(inputs.water_flux_m3_per_step / inputs.interface_area_m2));

    for (0..species_count) |species| {
        const group = diffusivityGroup(species);
        const conductance_m3_per_step = (inputs.aqueous_diffusivity_m2_per_step[group] * tortuosity_per_m +
            dispersion_m_per_step) * inputs.interface_area_m2;
        const partition = if (species >= 42)
            @min(inputs.current.band_nitrate_fraction, inputs.adjacent.band_nitrate_fraction)
        else if (species >= 34)
            @min(inputs.current.nonband_phosphate_fraction, inputs.adjacent.nonband_phosphate_fraction)
        else
            1;
        flux[species] = conductance_m3_per_step *
            (current_concentration[species] - adjacent_concentration[species]) * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteSoilBoundaryDiffusiveSoluteResult;
    }
    return flux;
}

fn concentrations(layer: Layer, phosphate_zero_threshold_m3: f64, adjacent_legacy_denominator: bool) [species_count]f64 {
    var values: [species_count]f64 = undefined;
    for (0..34) |species|
        values[species] = @max(0, layer.solute_inventory_amount[species] / layer.micropore_water_m3);

    const nonband_water_m3 = layer.micropore_water_m3 * layer.nonband_phosphate_fraction;
    for (34..42) |species| {
        values[species] = if (nonband_water_m3 > phosphate_zero_threshold_m3)
            @max(0, layer.solute_inventory_amount[species] /
                (if (adjacent_legacy_denominator) layer.micropore_water_m3 else nonband_water_m3))
        else
            0;
    }
    const band_water_m3 = layer.micropore_water_m3 * layer.band_phosphate_fraction;
    for (42..50) |species| {
        values[species] = if (band_water_m3 > phosphate_zero_threshold_m3)
            @max(0, layer.solute_inventory_amount[species] /
                (if (adjacent_legacy_denominator) layer.micropore_water_m3 else band_water_m3))
        else
            values[species - 8];
    }
    return values;
}

fn diffusivityGroup(species: usize) usize {
    return switch (species) {
        0, 12...16 => 1,
        1, 17...21 => 2,
        2 => 3,
        3, 22...25 => 4,
        4, 26...29 => 5,
        5, 30, 31 => 6,
        6, 32 => 7,
        7 => 8,
        8 => 9,
        9 => 10,
        10 => 11,
        11 => 12,
        33 => 13,
        34...49 => 0,
        else => unreachable,
    };
}

fn validateTopologyAndGateInputs(inputs: Inputs) !void {
    if (inputs.current.solute_inventory_amount.len != species_count or inputs.adjacent.solute_inventory_amount.len != species_count)
        return error.SoilBoundaryDiffusiveSoluteDimensionMismatch;
    try validateLayerGate(inputs.current);
    try validateLayerGate(inputs.adjacent);
}

fn validateLayerGate(layer: Layer) !void {
    inline for (.{ layer.volumetric_water_content, layer.minimum_volumetric_water_content, layer.micropore_water_m3, layer.minimum_water_m3 }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryDiffusiveSoluteInput;
    if (layer.volumetric_water_content < 0 or layer.minimum_volumetric_water_content < 0 or
        layer.micropore_water_m3 < 0 or layer.minimum_water_m3 < 0)
        return error.InvalidSoilBoundaryDiffusiveSoluteInput;
}

fn validateActiveInputs(inputs: Inputs) !void {
    inline for (.{ inputs.minimum_length_m, inputs.maximum_flux_velocity_m_per_step, inputs.water_flux_m3_per_step, inputs.interface_area_m2, inputs.dispersivity }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryDiffusiveSoluteInput;
    if (inputs.minimum_length_m <= 0 or inputs.maximum_flux_velocity_m_per_step < 0 or inputs.interface_area_m2 <= 0 or inputs.dispersivity < 0)
        return error.InvalidSoilBoundaryDiffusiveSoluteInput;
    for (inputs.aqueous_diffusivity_m2_per_step) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilBoundaryDiffusiveSoluteInput;
    try validateActiveLayer(inputs.current);
    try validateActiveLayer(inputs.adjacent);
}

fn validateActiveLayer(layer: Layer) !void {
    inline for (.{ layer.thickness_m, layer.tortuosity, layer.nonband_phosphate_fraction, layer.band_phosphate_fraction, layer.band_nitrate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryDiffusiveSoluteInput;
    if (layer.thickness_m < 0 or layer.tortuosity < 0 or
        layer.nonband_phosphate_fraction < 0 or layer.nonband_phosphate_fraction > 1 or
        layer.band_phosphate_fraction < 0 or layer.band_phosphate_fraction > 1 or
        layer.band_nitrate_fraction < 0 or layer.band_nitrate_fraction > 1)
        return error.InvalidSoilBoundaryDiffusiveSoluteInput;
    for (layer.solute_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryDiffusiveSoluteInput;
}

fn fixtureLayer(inventory: []const f64) Layer {
    return .{ .volumetric_water_content = 0.3, .minimum_volumetric_water_content = 0.05, .micropore_water_m3 = 2, .minimum_water_m3 = 0.01, .thickness_m = 0.2, .tortuosity = 0.1, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .band_nitrate_fraction = 0.4, .solute_inventory_amount = inventory };
}

fn fixture(current: Layer, adjacent: Layer) Inputs {
    return .{ .current = current, .adjacent = adjacent, .minimum_length_m = 0.001, .maximum_flux_velocity_m_per_step = 3, .water_flux_m3_per_step = 0, .interface_area_m2 = 2, .dispersivity = 0, .aqueous_diffusivity_m2_per_step = [_]f64{1} ** diffusivity_group_count };
}

test "TRNSFRS micropore diffusion maps all diffusivity groups" {
    var current_inventory = [_]f64{4} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    current_inventory[34] = 1;
    current_inventory[42] = 2;
    const flux = try calculate(fixture(fixtureLayer(&current_inventory), fixtureLayer(&adjacent_inventory)));
    try std.testing.expectEqual(@as(f64, 2), flux[0]);
    try std.testing.expectEqual(@as(f64, 0.5), flux[34]);
    try std.testing.expectEqual(@as(f64, 0.8), flux[42]);
}

test "TRNSFRS preserves adjacent phosphate denominator asymmetry" {
    const current_inventory = [_]f64{0} ** species_count;
    var adjacent_inventory = [_]f64{0} ** species_count;
    adjacent_inventory[34] = 2;
    const flux = try calculate(fixture(fixtureLayer(&current_inventory), fixtureLayer(&adjacent_inventory)));
    try std.testing.expectEqual(@as(f64, -0.25), flux[34]);
}

test "TRNSFRS dry band falls back to corresponding nonband P concentration" {
    var current_inventory = [_]f64{0} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    current_inventory[34] = 1;
    var current = fixtureLayer(&current_inventory);
    current.band_phosphate_fraction = 0;
    const flux = try calculate(fixture(current, fixtureLayer(&adjacent_inventory)));
    try std.testing.expectEqual(@as(f64, 0.8), flux[42]);
}

test "TRNSFRS water threshold clears every diffusive flux" {
    const current_inventory = [_]f64{4} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    var current = fixtureLayer(&current_inventory);
    current.volumetric_water_content = current.minimum_volumetric_water_content;
    const flux = try calculate(fixture(current, fixtureLayer(&adjacent_inventory)));
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "TRNSFRS dry pair does not evaluate dormant diffusion state" {
    const current_inventory = [_]f64{4} ** species_count;
    var adjacent_inventory = [_]f64{0} ** species_count;
    adjacent_inventory[49] = std.math.nan(f64);
    var current = fixtureLayer(&current_inventory);
    current.volumetric_water_content = current.minimum_volumetric_water_content;
    var inputs = fixture(current, fixtureLayer(&adjacent_inventory));
    inputs.interface_area_m2 = 0;
    inputs.aqueous_diffusivity_m2_per_step[13] = std.math.nan(f64);
    const flux = try calculate(inputs);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "invalid late inventory fails before publishing a result" {
    const current_inventory = [_]f64{4} ** species_count;
    var adjacent_inventory = [_]f64{0} ** species_count;
    adjacent_inventory[49] = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundaryDiffusiveSoluteInput, calculate(fixture(fixtureLayer(&current_inventory), fixtureLayer(&adjacent_inventory))));
}
