const std = @import("std");

pub const snow_species_per_layer = 41;
pub const litter_species_count = 42;
pub const soil_species_count = 50;
const pre_phosphate_species_count = 33;
const phosphate_species_count = 8;

pub const SurfaceWaterRouting = struct {
    litter_water_flux_m3_per_step: f64,
    soil_micropore_water_flux_m3_per_step: f64,
    soil_macropore_water_flux_m3_per_step: f64,
    litter_cover_fraction: f64,
    bare_soil_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
};

pub const Discharge = struct {
    /// Named replacement for `ICHKL`; at most one snow layer may discharge.
    surface_discharge_claimed: *bool,
    heat_capacity_megajoules_per_k: []const f64,
    minimum_heat_capacity_megajoules_per_k: f64,
    liquid_water_m3: []const f64,
    /// Cell-specific numerical water threshold (`ZEROS2`), m3.
    minimum_liquid_water_m3: f64,
    inventory_mol_per_layer: []const f64,
    water: SurfaceWaterRouting,
    /// 34 salt/complex slots plus eight non-band phosphate slots, mol step-1.
    litter_flux_mol_per_step: []f64,
    /// Same first 42 slots plus eight band phosphate slots, mol step-1.
    soil_flux_mol_per_step: []f64,
};

/// Exact compatibility translation of TRNSFRS.F lines 2186--2293.
///
/// This preserves the first eligible bottom-layer claim, legacy dry-layer
/// cover fractions, explicit zero H4SiO4 outputs, and phosphate partitioning.
pub fn route(discharge: Discharge) !void {
    const layer_count = discharge.heat_capacity_megajoules_per_k.len;
    if (discharge.liquid_water_m3.len != layer_count)
        return error.SnowpackSurfacePhysicalDimensionMismatch;
    const inventory_len = std.math.mul(usize, layer_count, snow_species_per_layer) catch
        return error.SnowpackSurfaceSoluteDimensionOverflow;
    if (discharge.inventory_mol_per_layer.len != inventory_len or
        discharge.litter_flux_mol_per_step.len != litter_species_count or
        discharge.soil_flux_mol_per_step.len != soil_species_count)
        return error.SnowpackSurfaceSoluteDimensionMismatch;
    try validateScalarInputs(discharge);
    if (discharge.surface_discharge_claimed.*) return;

    const layer = eligibleLayer(discharge) orelse return;
    const water_m3 = discharge.liquid_water_m3[layer];
    const litter_fraction = if (water_m3 > discharge.minimum_liquid_water_m3)
        std.math.clamp(discharge.water.litter_water_flux_m3_per_step / water_m3, 0, 1)
    else
        discharge.water.litter_cover_fraction;
    const soil_fraction = if (water_m3 > discharge.minimum_liquid_water_m3)
        std.math.clamp((discharge.water.soil_micropore_water_flux_m3_per_step +
            discharge.water.soil_macropore_water_flux_m3_per_step) / water_m3, 0, 1)
    else
        discharge.water.bare_soil_fraction;
    const source_start = layer * snow_species_per_layer;
    for (0..snow_species_per_layer) |species| {
        const inventory_mol = discharge.inventory_mol_per_layer[source_start + species];
        if (!std.math.isFinite(inventory_mol)) return error.NonFiniteSnowpackSurfaceInput;
        const litter_mol = inventory_mol * litter_fraction;
        const soil_mol = inventory_mol * soil_fraction;
        if (!std.math.isFinite(litter_mol) or !std.math.isFinite(soil_mol))
            return error.NonFiniteSnowpackSurfaceResult;
        if (species >= pre_phosphate_species_count and
            (!std.math.isFinite(soil_mol * discharge.water.nonband_phosphate_fraction) or
                !std.math.isFinite(soil_mol * discharge.water.band_phosphate_fraction)))
            return error.NonFiniteSnowpackSurfaceResult;
    }

    for (0..pre_phosphate_species_count) |species| {
        const inventory_mol = discharge.inventory_mol_per_layer[source_start + species];
        discharge.litter_flux_mol_per_step[species] = inventory_mol * litter_fraction;
        discharge.soil_flux_mol_per_step[species] = inventory_mol * soil_fraction;
    }
    discharge.litter_flux_mol_per_step[pre_phosphate_species_count] = 0;
    discharge.soil_flux_mol_per_step[pre_phosphate_species_count] = 0;
    for (0..phosphate_species_count) |phosphate| {
        const inventory_mol = discharge.inventory_mol_per_layer[
            source_start + pre_phosphate_species_count + phosphate
        ];
        discharge.litter_flux_mol_per_step[34 + phosphate] = inventory_mol * litter_fraction;
        discharge.soil_flux_mol_per_step[34 + phosphate] =
            inventory_mol * soil_fraction * discharge.water.nonband_phosphate_fraction;
        discharge.soil_flux_mol_per_step[42 + phosphate] =
            inventory_mol * soil_fraction * discharge.water.band_phosphate_fraction;
    }
    discharge.surface_discharge_claimed.* = true;
}

fn eligibleLayer(discharge: Discharge) ?usize {
    for (0..discharge.heat_capacity_megajoules_per_k.len) |layer| {
        if (discharge.heat_capacity_megajoules_per_k[layer] <= discharge.minimum_heat_capacity_megajoules_per_k)
            continue;
        const lower_active = layer + 1 < discharge.heat_capacity_megajoules_per_k.len and
            discharge.heat_capacity_megajoules_per_k[layer + 1] > discharge.minimum_heat_capacity_megajoules_per_k;
        if (!lower_active) return layer;
    }
    return null;
}

fn validateScalarInputs(discharge: Discharge) !void {
    if (!std.math.isFinite(discharge.minimum_heat_capacity_megajoules_per_k) or
        discharge.minimum_heat_capacity_megajoules_per_k < 0 or
        !std.math.isFinite(discharge.minimum_liquid_water_m3) or
        discharge.minimum_liquid_water_m3 < 0)
        return error.NonFiniteSnowpackSurfaceInput;
    for (discharge.heat_capacity_megajoules_per_k) |value|
        if (!std.math.isFinite(value) or value < 0) return error.NonFiniteSnowpackSurfaceInput;
    for (discharge.liquid_water_m3) |value|
        if (!std.math.isFinite(value) or value < 0) return error.NonFiniteSnowpackSurfaceInput;
    inline for (@typeInfo(SurfaceWaterRouting).@"struct".fields) |field| {
        const value = @field(discharge.water, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSnowpackSurfaceInput;
    }
    inline for (.{ discharge.water.litter_cover_fraction, discharge.water.bare_soil_fraction, discharge.water.nonband_phosphate_fraction, discharge.water.band_phosphate_fraction }) |fraction|
        if (fraction < 0 or fraction > 1) return error.NonFiniteSnowpackSurfaceInput;
}

fn fixture(
    claimed: *bool,
    heat: []const f64,
    water_m3: []const f64,
    inventory: []const f64,
    litter: []f64,
    soil: []f64,
) Discharge {
    return .{
        .surface_discharge_claimed = claimed,
        .heat_capacity_megajoules_per_k = heat,
        .minimum_heat_capacity_megajoules_per_k = 1,
        .liquid_water_m3 = water_m3,
        .minimum_liquid_water_m3 = 0.1,
        .inventory_mol_per_layer = inventory,
        .water = .{ .litter_water_flux_m3_per_step = 2, .soil_micropore_water_flux_m3_per_step = 1, .soil_macropore_water_flux_m3_per_step = 1, .litter_cover_fraction = 0.6, .bare_soil_fraction = 0.4, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.75 },
        .litter_flux_mol_per_step = litter,
        .soil_flux_mol_per_step = soil,
    };
}

test "TRNSFRS first eligible bottom layer discharges with phosphate partitioning" {
    const heat = [_]f64{ 2, 2 };
    const water_m3 = [_]f64{ 9, 4 };
    var inventory = [_]f64{0} ** (2 * snow_species_per_layer);
    @memset(inventory[snow_species_per_layer..], 8);
    var litter = [_]f64{99} ** litter_species_count;
    var soil = [_]f64{99} ** soil_species_count;
    var claimed = false;
    try route(fixture(&claimed, &heat, &water_m3, &inventory, &litter, &soil));
    try std.testing.expect(claimed);
    try std.testing.expectEqual(@as(f64, 4), litter[0]);
    try std.testing.expectEqual(@as(f64, 0), litter[33]);
    try std.testing.expectEqual(@as(f64, 4), litter[34]);
    try std.testing.expectEqual(@as(f64, 1), soil[34]);
    try std.testing.expectEqual(@as(f64, 3), soil[42]);
}

test "TRNSFRS dry layer uses cover fractions" {
    const heat = [_]f64{2};
    const water_m3 = [_]f64{0};
    const inventory = [_]f64{10} ** snow_species_per_layer;
    var litter = [_]f64{0} ** litter_species_count;
    var soil = [_]f64{0} ** soil_species_count;
    var claimed = false;
    try route(fixture(&claimed, &heat, &water_m3, &inventory, &litter, &soil));
    try std.testing.expectEqual(@as(f64, 6), litter[0]);
    try std.testing.expectEqual(@as(f64, 4), soil[0]);
}

test "claimed discharge and zero layers preserve outputs" {
    const empty: [0]f64 = .{};
    var litter = [_]f64{9} ** litter_species_count;
    var soil = [_]f64{7} ** soil_species_count;
    var claimed = true;
    try route(fixture(&claimed, &empty, &empty, &empty, &litter, &soil));
    try std.testing.expectEqual(@as(f64, 9), litter[0]);
    try std.testing.expectEqual(@as(f64, 7), soil[0]);
}

test "late non-finite inventory failure is atomic" {
    const heat = [_]f64{2};
    const water_m3 = [_]f64{4};
    var inventory = [_]f64{8} ** snow_species_per_layer;
    inventory[snow_species_per_layer - 1] = std.math.inf(f64);
    var litter = [_]f64{9} ** litter_species_count;
    var soil = [_]f64{7} ** soil_species_count;
    var claimed = false;
    try std.testing.expectError(error.NonFiniteSnowpackSurfaceInput, route(fixture(&claimed, &heat, &water_m3, &inventory, &litter, &soil)));
    try std.testing.expect(!claimed);
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** litter_species_count), &litter);
    try std.testing.expectEqualSlices(f64, &([_]f64{7} ** soil_species_count), &soil);
}
