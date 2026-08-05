const std = @import("std");

pub const species_count = 50;
pub const Direction = enum { east_west, north_south, vertical };
pub const Side = enum { forward, reverse };
pub const CellConnectivity = enum { interconnected, standalone };

pub const Inputs = struct {
    direction: Direction,
    side: Side,
    connectivity: CellConnectivity,
    source_layer_thickness_m: f64,
    minimum_transport_layer_thickness_m: f64,
    boundary_water_flux_m3_per_step: f64,
    source_micropore_water_m3: f64,
    minimum_water_m3: f64,
    maximum_convective_fraction: f64,
    minimum_convective_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    source_inventory_amount: []const f64,
    prescribed_recharge_concentration_amount_per_m3: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 7736--7947.
/// Recharge uses raw FLWM rather than clamped VFLW, forces H4SiO4 to zero,
/// and consumes caller-provided legacy concentration ordering (including the
/// RALSFS assignment from CALU at source line 7843).
pub fn calculate(inputs: Inputs) !?[species_count]f64 {
    try validateGuardInputs(inputs);
    if (inputs.source_layer_thickness_m <= inputs.minimum_transport_layer_thickness_m or
        (inputs.connectivity == .standalone and inputs.direction != .vertical))
        return null;

    try validateHydraulicControls(inputs);
    var flux = [_]f64{0} ** species_count;
    const convective_fraction = if (inputs.source_micropore_water_m3 > inputs.minimum_water_m3)
        std.math.clamp(inputs.boundary_water_flux_m3_per_step / inputs.source_micropore_water_m3, -inputs.maximum_convective_fraction, inputs.maximum_convective_fraction)
    else
        0;
    if (@abs(convective_fraction) <= inputs.minimum_convective_fraction) return flux;

    const discharge = switch (inputs.side) {
        .forward => inputs.boundary_water_flux_m3_per_step > 0,
        .reverse => inputs.boundary_water_flux_m3_per_step < 0,
    };
    const recharge = switch (inputs.side) {
        .forward => inputs.boundary_water_flux_m3_per_step < 0,
        .reverse => inputs.boundary_water_flux_m3_per_step > 0,
    };
    if (!discharge and !recharge) return flux;

    try validateActiveChemistry(inputs, discharge);

    for (0..species_count) |species| {
        const partition = if (species >= 42)
            inputs.band_phosphate_fraction
        else if (species >= 34)
            inputs.nonband_phosphate_fraction
        else
            1;
        flux[species] = if (discharge)
            convective_fraction * @max(0, inputs.source_inventory_amount[species]) * partition
        else if (species == 33)
            0
        else
            inputs.boundary_water_flux_m3_per_step * inputs.prescribed_recharge_concentration_amount_per_m3[species] * partition;
        if (!std.math.isFinite(flux[species])) return error.NonFiniteExternalMicroporeConvectiveSoluteResult;
    }
    return flux;
}

fn validateGuardInputs(inputs: Inputs) !void {
    inline for (.{ inputs.source_layer_thickness_m, inputs.minimum_transport_layer_thickness_m }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalMicroporeConvectiveSoluteInput;
    if (inputs.source_layer_thickness_m < 0 or inputs.minimum_transport_layer_thickness_m < 0)
        return error.InvalidExternalMicroporeConvectiveSoluteInput;
}

fn validateHydraulicControls(inputs: Inputs) !void {
    if (inputs.source_inventory_amount.len != species_count or
        inputs.prescribed_recharge_concentration_amount_per_m3.len != species_count)
        return error.ExternalMicroporeConvectiveSoluteDimensionMismatch;
    inline for (.{ inputs.boundary_water_flux_m3_per_step, inputs.source_micropore_water_m3, inputs.minimum_water_m3, inputs.maximum_convective_fraction, inputs.minimum_convective_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalMicroporeConvectiveSoluteInput;
    if (inputs.source_micropore_water_m3 < 0 or inputs.minimum_water_m3 < 0 or
        inputs.maximum_convective_fraction < 0 or inputs.minimum_convective_fraction < 0)
        return error.InvalidExternalMicroporeConvectiveSoluteInput;
}

fn validateActiveChemistry(inputs: Inputs, discharge: bool) !void {
    inline for (.{ inputs.nonband_phosphate_fraction, inputs.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidExternalMicroporeConvectiveSoluteInput;
    const active_values = if (discharge)
        inputs.source_inventory_amount
    else
        inputs.prescribed_recharge_concentration_amount_per_m3;
    for (active_values) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalMicroporeConvectiveSoluteInput;
}

fn fixture(inventory: []const f64, concentration: []const f64, side: Side, water_flux: f64) Inputs {
    return .{ .direction = .east_west, .side = side, .connectivity = .interconnected, .source_layer_thickness_m = 0.2, .minimum_transport_layer_thickness_m = 0.01, .boundary_water_flux_m3_per_step = water_flux, .source_micropore_water_m3 = 4, .minimum_water_m3 = 0.01, .maximum_convective_fraction = 0.75, .minimum_convective_fraction = 1e-9, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .source_inventory_amount = inventory, .prescribed_recharge_concentration_amount_per_m3 = concentration };
}

test "TRNSFRS forward discharge uses clamped inventory fraction and P partitions" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{99} ** species_count;
    const flux = (try calculate(fixture(&inventory, &concentration, .forward, 2))).?;
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 4), flux[33]);
    try std.testing.expectEqual(@as(f64, 1), flux[34]);
    try std.testing.expectEqual(@as(f64, 2), flux[42]);
}

test "TRNSFRS forward recharge uses raw water flux and forces H4SiO4 zero" {
    const inventory = [_]f64{99} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    const flux = (try calculate(fixture(&inventory, &concentration, .forward, -2))).?;
    try std.testing.expectEqual(@as(f64, -6), flux[0]);
    try std.testing.expectEqual(@as(f64, 0), flux[33]);
    try std.testing.expectEqual(@as(f64, -1.5), flux[34]);
    try std.testing.expectEqual(@as(f64, -3), flux[42]);
}

test "TRNSFRS reverse side swaps discharge and recharge signs" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    const discharge = (try calculate(fixture(&inventory, &concentration, .reverse, -2))).?;
    const recharge = (try calculate(fixture(&inventory, &concentration, .reverse, 2))).?;
    try std.testing.expectEqual(@as(f64, -4), discharge[0]);
    try std.testing.expectEqual(@as(f64, 6), recharge[0]);
}

test "TRNSFRS thin or standalone horizontal boundary skips publication" {
    const inventory = [_]f64{8} ** species_count;
    const concentration = [_]f64{3} ** species_count;
    var inputs = fixture(&inventory, &concentration, .forward, 2);
    inputs.connectivity = .standalone;
    inputs.source_inventory_amount = &[_]f64{};
    inputs.prescribed_recharge_concentration_amount_per_m3 = &[_]f64{};
    try std.testing.expect((try calculate(inputs)) == null);
    inputs.source_inventory_amount = &inventory;
    inputs.prescribed_recharge_concentration_amount_per_m3 = &concentration;
    inputs.direction = .vertical;
    try std.testing.expect((try calculate(inputs)) != null);
}

test "TRNSFRS discharge and recharge ignore the opposite dormant chemistry" {
    var inventory = [_]f64{8} ** species_count;
    var concentration = [_]f64{3} ** species_count;
    concentration[49] = std.math.nan(f64);
    _ = (try calculate(fixture(&inventory, &concentration, .forward, 2))).?;
    concentration[49] = 3;
    inventory[49] = std.math.nan(f64);
    _ = (try calculate(fixture(&inventory, &concentration, .forward, -2))).?;
}

test "late invalid prescribed concentration fails atomically" {
    const inventory = [_]f64{8} ** species_count;
    var concentration = [_]f64{3} ** species_count;
    concentration[49] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteExternalMicroporeConvectiveSoluteInput, calculate(fixture(&inventory, &concentration, .forward, -2)));
}
