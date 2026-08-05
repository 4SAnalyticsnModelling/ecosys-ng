const std = @import("std");

pub const species_count = 49;
pub const diffusivity_group_count = 14;

pub const Layer = struct {
    macropore_water_m3: f64,
    macropore_volume_m3: f64,
    minimum_volumetric_water_content: f64,
    thickness_m: f64,
    macropore_tortuosity: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    /// 33 salts/complexes (no H4SiO4), then 8 non-band and 8 band P species.
    solute_inventory_amount: []const f64,
};

pub const Inputs = struct {
    current: Layer,
    adjacent: Layer,
    minimum_length_m: f64,
    interface_area_m2: f64,
    /// PO4, Al, Fe, H, Ca, Mg, Na, K, OH, SO4, Cl, CO3, HCO3, H4SiO4; m2 step-1.
    aqueous_diffusivity_m2_per_step: [diffusivity_group_count]f64,
};

/// Compatibility translation of TRNSFRS.F lines 5751--6028.
/// Macropore phosphate concentrations divide by total macropore water and
/// their fluxes use the adjacent layer's partition fraction, not a minimum.
pub fn calculate(inputs: Inputs) ![species_count]f64 {
    try validateTopologyAndGateInputs(inputs);
    var flux = [_]f64{0} ** species_count;
    if (inputs.current.macropore_water_m3 <= inputs.current.minimum_volumetric_water_content * inputs.current.macropore_volume_m3 or
        inputs.adjacent.macropore_water_m3 <= inputs.adjacent.minimum_volumetric_water_content * inputs.adjacent.macropore_volume_m3)
        return flux;
    try validateActiveInputs(inputs);

    var current_concentration: [species_count]f64 = undefined;
    var adjacent_concentration: [species_count]f64 = undefined;
    for (0..species_count) |species| {
        current_concentration[species] = @max(0, inputs.current.solute_inventory_amount[species] / inputs.current.macropore_water_m3);
        adjacent_concentration[species] = @max(0, inputs.adjacent.solute_inventory_amount[species] / inputs.adjacent.macropore_water_m3);
    }

    const current_thickness = @max(inputs.minimum_length_m, inputs.current.thickness_m);
    const adjacent_thickness = @max(inputs.minimum_length_m, inputs.adjacent.thickness_m);
    const tortuosity_per_m = (inputs.current.macropore_tortuosity + inputs.adjacent.macropore_tortuosity) /
        (current_thickness + adjacent_thickness);
    var conductance_m3_per_step: [diffusivity_group_count]f64 = undefined;
    for (0..diffusivity_group_count) |group| {
        conductance_m3_per_step[group] = inputs.aqueous_diffusivity_m2_per_step[group] * tortuosity_per_m * inputs.interface_area_m2;
        if (!std.math.isFinite(conductance_m3_per_step[group])) return error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteResult;
    }

    for (0..species_count) |species| {
        const partition = if (species >= 41)
            inputs.adjacent.band_phosphate_fraction
        else if (species >= 33)
            inputs.adjacent.nonband_phosphate_fraction
        else
            1;
        flux[species] = conductance_m3_per_step[diffusivityGroup(species)] *
            (current_concentration[species] - adjacent_concentration[species]) * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteResult;
    }
    return flux;
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
        33...48 => 0,
        else => unreachable,
    };
}

fn validateTopologyAndGateInputs(inputs: Inputs) !void {
    if (inputs.current.solute_inventory_amount.len != species_count or inputs.adjacent.solute_inventory_amount.len != species_count)
        return error.SoilBoundaryMacroporeDiffusiveSoluteDimensionMismatch;
    try validateLayerGate(inputs.current);
    try validateLayerGate(inputs.adjacent);
}

fn validateActiveInputs(inputs: Inputs) !void {
    inline for (.{ inputs.minimum_length_m, inputs.interface_area_m2 }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteInput;
    if (inputs.minimum_length_m <= 0 or inputs.interface_area_m2 <= 0)
        return error.InvalidSoilBoundaryMacroporeDiffusiveSoluteInput;
    for (inputs.aqueous_diffusivity_m2_per_step) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilBoundaryMacroporeDiffusiveSoluteInput;
    try validateActiveLayer(inputs.current);
    try validateActiveLayer(inputs.adjacent);
}

fn validateLayerGate(layer: Layer) !void {
    inline for (.{ layer.macropore_water_m3, layer.macropore_volume_m3, layer.minimum_volumetric_water_content }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteInput;
    if (layer.macropore_water_m3 < 0 or layer.macropore_volume_m3 < 0 or layer.minimum_volumetric_water_content < 0)
        return error.InvalidSoilBoundaryMacroporeDiffusiveSoluteInput;
}

fn validateActiveLayer(layer: Layer) !void {
    inline for (.{ layer.thickness_m, layer.macropore_tortuosity, layer.nonband_phosphate_fraction, layer.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteInput;
    if (layer.thickness_m < 0 or layer.macropore_tortuosity < 0 or layer.nonband_phosphate_fraction < 0 or
        layer.nonband_phosphate_fraction > 1 or layer.band_phosphate_fraction < 0 or layer.band_phosphate_fraction > 1)
        return error.InvalidSoilBoundaryMacroporeDiffusiveSoluteInput;
    for (layer.solute_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteInput;
}

fn fixtureLayer(inventory: []const f64, nonband: f64, band: f64) Layer {
    return .{ .macropore_water_m3 = 2, .macropore_volume_m3 = 4, .minimum_volumetric_water_content = 0.05, .thickness_m = 0.2, .macropore_tortuosity = 0.1, .nonband_phosphate_fraction = nonband, .band_phosphate_fraction = band, .solute_inventory_amount = inventory };
}

fn fixture(current: Layer, adjacent: Layer) Inputs {
    return .{ .current = current, .adjacent = adjacent, .minimum_length_m = 0.001, .interface_area_m2 = 2, .aqueous_diffusivity_m2_per_step = [_]f64{1} ** diffusivity_group_count };
}

test "TRNSFRS macropore diffusion maps compact species and adjacent P fractions" {
    const current_inventory = [_]f64{4} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    const flux = try calculate(fixture(fixtureLayer(&current_inventory, 0.2, 0.3), fixtureLayer(&adjacent_inventory, 0.25, 0.5)));
    try std.testing.expectEqual(@as(f64, 2), flux[0]);
    try std.testing.expectEqual(@as(f64, 0.5), flux[33]);
    try std.testing.expectEqual(@as(f64, 1), flux[41]);
}

test "TRNSFRS macropore diffusion uses adjacent rather than minimum P partition" {
    const current_inventory = [_]f64{4} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    const flux = try calculate(fixture(fixtureLayer(&current_inventory, 0.01, 0.02), fixtureLayer(&adjacent_inventory, 0.5, 0.75)));
    try std.testing.expectEqual(@as(f64, 1), flux[33]);
    try std.testing.expectEqual(@as(f64, 1.5), flux[41]);
}

test "TRNSFRS saturation threshold clears all macropore diffusive fluxes" {
    const current_inventory = [_]f64{4} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    var current = fixtureLayer(&current_inventory, 0.2, 0.3);
    current.macropore_water_m3 = current.minimum_volumetric_water_content * current.macropore_volume_m3;
    const flux = try calculate(fixture(current, fixtureLayer(&adjacent_inventory, 0.25, 0.5)));
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "TRNSFRS dry macropore pair does not evaluate dormant diffusion state" {
    const current_inventory = [_]f64{4} ** species_count;
    var adjacent_inventory = [_]f64{0} ** species_count;
    adjacent_inventory[48] = std.math.nan(f64);
    var current = fixtureLayer(&current_inventory, 0.2, 0.3);
    current.macropore_water_m3 = current.minimum_volumetric_water_content * current.macropore_volume_m3;
    var inputs = fixture(current, fixtureLayer(&adjacent_inventory, 0.25, 0.5));
    inputs.interface_area_m2 = 0;
    inputs.aqueous_diffusivity_m2_per_step[13] = std.math.nan(f64);
    const flux = try calculate(inputs);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "TRNSFRS computes and validates the legacy unused H4SiO4 conductance" {
    const current_inventory = [_]f64{4} ** species_count;
    const adjacent_inventory = [_]f64{0} ** species_count;
    var inputs = fixture(fixtureLayer(&current_inventory, 0.2, 0.3), fixtureLayer(&adjacent_inventory, 0.25, 0.5));
    inputs.aqueous_diffusivity_m2_per_step[13] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilBoundaryMacroporeDiffusiveSoluteInput, calculate(inputs));
}

test "late invalid inventory fails before publishing a result" {
    const current_inventory = [_]f64{4} ** species_count;
    var adjacent_inventory = [_]f64{0} ** species_count;
    adjacent_inventory[48] = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundaryMacroporeDiffusiveSoluteInput, calculate(fixture(fixtureLayer(&current_inventory, 0.2, 0.3), fixtureLayer(&adjacent_inventory, 0.25, 0.5))));
}
