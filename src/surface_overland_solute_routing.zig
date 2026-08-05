const std = @import("std");

/// Litter runoff includes H4SiO4: 34 salt/complex plus eight phosphate fields.
pub const runoff_species_count = 42;
/// Snow drift omits H4SiO4: 33 salt/complex plus eight phosphate fields.
pub const snow_drift_species_count = 41;

pub const Inputs = struct {
    runoff_water_m3_per_step: f64,
    snow_drift_water_m3_per_step: f64,
    minimum_overland_water_flux_m3_per_step: f64,
    litter_water_m3: f64,
    top_snow_water_m3: f64,
    minimum_donor_water_m3: f64,
    maximum_transport_fraction: f64,
    litter_inventory_mol: []const f64,
    top_snow_inventory_mol: []const f64,
};

/// Exact source-order translation of TRNSFRS.F lines 3829--4012.
/// Outputs are mol step-1. The caller owns runtime grid and substep axes; this
/// kernel handles one cell and validates both independent routes atomically.
pub fn calculate(
    inputs: Inputs,
    runoff_flux_mol_per_step: []f64,
    snow_drift_flux_mol_per_step: []f64,
) !void {
    if (inputs.litter_inventory_mol.len != runoff_species_count or
        inputs.top_snow_inventory_mol.len != snow_drift_species_count or
        runoff_flux_mol_per_step.len != runoff_species_count or
        snow_drift_flux_mol_per_step.len != snow_drift_species_count)
        return error.SurfaceOverlandSoluteRoutingDimensionMismatch;
    try validate(inputs);

    const runoff_fraction = routeFraction(
        inputs.runoff_water_m3_per_step,
        inputs.minimum_overland_water_flux_m3_per_step,
        inputs.litter_water_m3,
        inputs.minimum_donor_water_m3,
        inputs.maximum_transport_fraction,
    );
    const snow_fraction = routeFraction(
        inputs.snow_drift_water_m3_per_step,
        inputs.minimum_overland_water_flux_m3_per_step,
        inputs.top_snow_water_m3,
        inputs.minimum_donor_water_m3,
        inputs.maximum_transport_fraction,
    );
    for (inputs.litter_inventory_mol) |inventory_mol| {
        if (!std.math.isFinite(runoff_fraction * @max(0, inventory_mol)))
            return error.NonFiniteSurfaceOverlandSoluteRoutingResult;
    }
    for (inputs.top_snow_inventory_mol) |inventory_mol| {
        if (!std.math.isFinite(snow_fraction * @max(0, inventory_mol)))
            return error.NonFiniteSurfaceOverlandSoluteRoutingResult;
    }
    for (inputs.litter_inventory_mol, runoff_flux_mol_per_step) |inventory_mol, *flux_mol|
        flux_mol.* = runoff_fraction * @max(0, inventory_mol);
    for (inputs.top_snow_inventory_mol, snow_drift_flux_mol_per_step) |inventory_mol, *flux_mol|
        flux_mol.* = snow_fraction * @max(0, inventory_mol);
}

fn routeFraction(flow_m3: f64, minimum_flow_m3: f64, donor_water_m3: f64, minimum_water_m3: f64, maximum_fraction: f64) f64 {
    if (!(flow_m3 > minimum_flow_m3)) return 0;
    return if (donor_water_m3 > minimum_water_m3)
        @min(maximum_fraction, flow_m3 / donor_water_m3)
    else
        maximum_fraction;
}

fn validate(inputs: Inputs) !void {
    inline for (.{ inputs.runoff_water_m3_per_step, inputs.snow_drift_water_m3_per_step, inputs.minimum_overland_water_flux_m3_per_step, inputs.litter_water_m3, inputs.top_snow_water_m3, inputs.minimum_donor_water_m3, inputs.maximum_transport_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandSoluteRoutingInput;
    if (inputs.runoff_water_m3_per_step < 0 or inputs.snow_drift_water_m3_per_step < 0 or
        inputs.minimum_overland_water_flux_m3_per_step < 0 or inputs.litter_water_m3 < 0 or
        inputs.top_snow_water_m3 < 0 or inputs.minimum_donor_water_m3 < 0 or
        inputs.maximum_transport_fraction < 0)
        return error.InvalidSurfaceOverlandSoluteRoutingInput;
    for (inputs.litter_inventory_mol) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandSoluteRoutingInput;
    for (inputs.top_snow_inventory_mol) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandSoluteRoutingInput;
}

fn fixture(litter: []const f64, snow: []const f64) Inputs {
    return .{ .runoff_water_m3_per_step = 2, .snow_drift_water_m3_per_step = 1, .minimum_overland_water_flux_m3_per_step = 0.01, .litter_water_m3 = 4, .top_snow_water_m3 = 4, .minimum_donor_water_m3 = 0.1, .maximum_transport_fraction = 0.75, .litter_inventory_mol = litter, .top_snow_inventory_mol = snow };
}

test "TRNSFRS routes exact runoff and compact snow-drift topologies" {
    const litter = [_]f64{8} ** runoff_species_count;
    const snow = [_]f64{8} ** snow_drift_species_count;
    var runoff = [_]f64{0} ** runoff_species_count;
    var snow_drift = [_]f64{0} ** snow_drift_species_count;
    try calculate(fixture(&litter, &snow), &runoff, &snow_drift);
    try std.testing.expectEqual(@as(f64, 4), runoff[0]);
    try std.testing.expectEqual(@as(f64, 4), runoff[runoff_species_count - 1]);
    try std.testing.expectEqual(@as(f64, 2), snow_drift[0]);
    try std.testing.expectEqual(@as(f64, 2), snow_drift[snow_drift_species_count - 1]);
}

test "independent route gates zero only the inactive route" {
    const litter = [_]f64{8} ** runoff_species_count;
    const snow = [_]f64{8} ** snow_drift_species_count;
    var runoff = [_]f64{9} ** runoff_species_count;
    var snow_drift = [_]f64{9} ** snow_drift_species_count;
    var inputs = fixture(&litter, &snow);
    inputs.runoff_water_m3_per_step = 0;
    try calculate(inputs, &runoff, &snow_drift);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** runoff_species_count), &runoff);
    try std.testing.expectEqual(@as(f64, 2), snow_drift[0]);
}

test "dry donor uses maximum transport fraction and clamps negative inventory" {
    var litter = [_]f64{8} ** runoff_species_count;
    litter[0] = -2;
    const snow = [_]f64{8} ** snow_drift_species_count;
    var runoff = [_]f64{0} ** runoff_species_count;
    var snow_drift = [_]f64{0} ** snow_drift_species_count;
    var inputs = fixture(&litter, &snow);
    inputs.litter_water_m3 = 0;
    try calculate(inputs, &runoff, &snow_drift);
    try std.testing.expectEqual(@as(f64, 0), runoff[0]);
    try std.testing.expectEqual(@as(f64, 6), runoff[1]);
}

test "late invalid snow inventory leaves both routes atomic" {
    const litter = [_]f64{8} ** runoff_species_count;
    var snow = [_]f64{8} ** snow_drift_species_count;
    snow[snow_drift_species_count - 1] = std.math.inf(f64);
    var runoff = [_]f64{9} ** runoff_species_count;
    var snow_drift = [_]f64{7} ** snow_drift_species_count;
    try std.testing.expectError(error.NonFiniteSurfaceOverlandSoluteRoutingInput, calculate(fixture(&litter, &snow), &runoff, &snow_drift));
    try std.testing.expectEqual(@as(f64, 9), runoff[0]);
    try std.testing.expectEqual(@as(f64, 7), snow_drift[0]);
}
