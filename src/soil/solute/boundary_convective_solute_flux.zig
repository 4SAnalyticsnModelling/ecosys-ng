const std = @import("std");

pub const species_count = 50;

pub const Layer = struct {
    thickness_m: f64,
    is_within_active_profile: bool,
    soil_volume_m3: f64,
    micropore_water_m3: f64,
    minimum_water_m3: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    /// Species 0...33 and 34...49 reproduce the exact legacy RFL ordering.
    /// Free non-band and band phosphate inventories retain their legacy amount units.
    solute_inventory_amount: []const f64,
};

pub const Inputs = struct {
    current: Layer,
    adjacent: Layer,
    minimum_transport_layer_thickness_m: f64,
    maximum_convective_fraction: f64,
    /// Positive means current to adjacent, non-positive means adjacent to current, m3 step-1.
    micropore_water_flux_m3_per_step: f64,
};

pub const Result = struct {
    current_volumetric_water_content: f64,
    adjacent_volumetric_water_content: f64,
    convective_fraction_per_step: f64,
    flux_amount_per_step: [species_count]f64,
};

/// Compatibility translation of TRNSFRS.F lines 4904--5068.
/// The phosphate multipliers intentionally come from the receiving layer.
pub fn calculate(inputs: Inputs) !?Result {
    try validateTopologyAndGateInputs(inputs);
    if (inputs.current.thickness_m <= inputs.minimum_transport_layer_thickness_m or
        inputs.adjacent.thickness_m <= inputs.minimum_transport_layer_thickness_m or
        !inputs.current.is_within_active_profile or !inputs.adjacent.is_within_active_profile)
        return null;
    try validateActiveLayer(inputs.current);
    try validateActiveLayer(inputs.adjacent);

    const current_theta = @max(0, inputs.current.micropore_water_m3 / inputs.current.soil_volume_m3);
    const adjacent_theta = @max(0, inputs.adjacent.micropore_water_m3 / inputs.adjacent.soil_volume_m3);
    const positive_flow = inputs.micropore_water_flux_m3_per_step > 0;
    const donor = if (positive_flow) inputs.current else inputs.adjacent;
    const receiver = if (positive_flow) inputs.adjacent else inputs.current;
    const fraction = if (positive_flow)
        if (donor.micropore_water_m3 > donor.minimum_water_m3)
            std.math.clamp(inputs.micropore_water_flux_m3_per_step / donor.micropore_water_m3, 0, inputs.maximum_convective_fraction)
        else
            inputs.maximum_convective_fraction
    else if (donor.micropore_water_m3 > donor.minimum_water_m3)
        std.math.clamp(inputs.micropore_water_flux_m3_per_step / donor.micropore_water_m3, -inputs.maximum_convective_fraction, 0)
    else
        -inputs.maximum_convective_fraction;

    var result = Result{
        .current_volumetric_water_content = current_theta,
        .adjacent_volumetric_water_content = adjacent_theta,
        .convective_fraction_per_step = fraction,
        .flux_amount_per_step = undefined,
    };
    for (0..species_count) |species| {
        const partition = if (species >= 42)
            receiver.band_phosphate_fraction
        else if (species >= 34)
            receiver.nonband_phosphate_fraction
        else
            1;
        result.flux_amount_per_step[species] = fraction * @max(0, donor.solute_inventory_amount[species]) * partition;
        if (!std.math.isFinite(result.flux_amount_per_step[species]))
            return error.NonFiniteSoilBoundaryConvectiveSoluteResult;
    }
    return result;
}

fn validateTopologyAndGateInputs(inputs: Inputs) !void {
    if (inputs.current.solute_inventory_amount.len != species_count or
        inputs.adjacent.solute_inventory_amount.len != species_count)
        return error.SoilBoundaryConvectiveSoluteDimensionMismatch;
    inline for (.{ inputs.minimum_transport_layer_thickness_m, inputs.maximum_convective_fraction, inputs.micropore_water_flux_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryConvectiveSoluteInput;
    if (inputs.minimum_transport_layer_thickness_m < 0 or inputs.maximum_convective_fraction < 0)
        return error.InvalidSoilBoundaryConvectiveSoluteInput;
    inline for (.{ inputs.current.thickness_m, inputs.adjacent.thickness_m }) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidSoilBoundaryConvectiveSoluteInput;
}

fn validateActiveLayer(layer: Layer) !void {
    inline for (.{ layer.soil_volume_m3, layer.micropore_water_m3, layer.minimum_water_m3, layer.nonband_phosphate_fraction, layer.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryConvectiveSoluteInput;
    if (layer.thickness_m < 0 or layer.soil_volume_m3 <= 0 or layer.micropore_water_m3 < 0 or
        layer.minimum_water_m3 < 0 or layer.nonband_phosphate_fraction < 0 or
        layer.nonband_phosphate_fraction > 1 or layer.band_phosphate_fraction < 0 or
        layer.band_phosphate_fraction > 1)
        return error.InvalidSoilBoundaryConvectiveSoluteInput;
    for (layer.solute_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBoundaryConvectiveSoluteInput;
}

fn fixtureLayer(inventory: []const f64) Layer {
    return .{ .thickness_m = 0.2, .is_within_active_profile = true, .soil_volume_m3 = 2, .micropore_water_m3 = 4, .minimum_water_m3 = 0.01, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .solute_inventory_amount = inventory };
}

test "TRNSFRS positive boundary flow uses current donor and adjacent P fractions" {
    const current_inventory = [_]f64{8} ** species_count;
    const adjacent_inventory = [_]f64{20} ** species_count;
    const result = (try calculate(.{ .current = fixtureLayer(&current_inventory), .adjacent = fixtureLayer(&adjacent_inventory), .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = 2 })).?;
    try std.testing.expectEqual(@as(f64, 2), result.current_volumetric_water_content);
    try std.testing.expectEqual(@as(f64, 0.5), result.convective_fraction_per_step);
    try std.testing.expectEqual(@as(f64, 4), result.flux_amount_per_step[0]);
    try std.testing.expectEqual(@as(f64, 1), result.flux_amount_per_step[34]);
    try std.testing.expectEqual(@as(f64, 2), result.flux_amount_per_step[42]);
}

test "TRNSFRS non-positive boundary flow uses adjacent donor and current P fractions" {
    const current_inventory = [_]f64{99} ** species_count;
    const adjacent_inventory = [_]f64{8} ** species_count;
    const result = (try calculate(.{ .current = fixtureLayer(&current_inventory), .adjacent = fixtureLayer(&adjacent_inventory), .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = -2 })).?;
    try std.testing.expectEqual(@as(f64, -0.5), result.convective_fraction_per_step);
    try std.testing.expectEqual(@as(f64, -4), result.flux_amount_per_step[0]);
    try std.testing.expectEqual(@as(f64, -1), result.flux_amount_per_step[34]);
    try std.testing.expectEqual(@as(f64, -2), result.flux_amount_per_step[42]);
}

test "TRNSFRS dry donor uses signed maximum fraction including zero flow" {
    const inventory = [_]f64{8} ** species_count;
    var current = fixtureLayer(&inventory);
    var adjacent = fixtureLayer(&inventory);
    current.micropore_water_m3 = 0;
    adjacent.micropore_water_m3 = 0;
    const positive = (try calculate(.{ .current = current, .adjacent = adjacent, .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = 1 })).?;
    const zero = (try calculate(.{ .current = current, .adjacent = adjacent, .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = 0 })).?;
    try std.testing.expectEqual(@as(f64, 0.75), positive.convective_fraction_per_step);
    try std.testing.expectEqual(@as(f64, -0.75), zero.convective_fraction_per_step);
}

test "TRNSFRS inactive or thin pair skips the complete boundary block" {
    const inventory = [_]f64{8} ** species_count;
    var adjacent = fixtureLayer(&inventory);
    adjacent.thickness_m = 0.01;
    try std.testing.expect((try calculate(.{ .current = fixtureLayer(&inventory), .adjacent = adjacent, .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = 1 })) == null);
}

test "TRNSFRS thin pair does not evaluate dormant layer state" {
    const inventory = [_]f64{8} ** species_count;
    var adjacent = fixtureLayer(&inventory);
    adjacent.thickness_m = 0;
    adjacent.soil_volume_m3 = 0;
    adjacent.micropore_water_m3 = std.math.nan(f64);
    try std.testing.expect((try calculate(.{ .current = fixtureLayer(&inventory), .adjacent = adjacent, .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = 1 })) == null);
}

test "invalid late species leaves no partially published result" {
    const current_inventory = [_]f64{8} ** species_count;
    var adjacent_inventory = [_]f64{8} ** species_count;
    adjacent_inventory[species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundaryConvectiveSoluteInput, calculate(.{ .current = fixtureLayer(&current_inventory), .adjacent = fixtureLayer(&adjacent_inventory), .minimum_transport_layer_thickness_m = 0.01, .maximum_convective_fraction = 0.75, .micropore_water_flux_m3_per_step = -1 }));
}
