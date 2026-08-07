const std = @import("std");

pub const Inputs = struct {
    root_axis_count: usize,
    soil_layer_count: usize,
    first_soil_layer: usize,
    end_soil_layer: usize,
    adjusted_soil_total_water_potential_mpa_by_layer: []const f64,
    minimum_dry_matter_fraction_g_c_per_g: f64,
    root_nonstructural_carbon_g_per_g_c_by_root_layer: []const f64,
    root_nonstructural_nitrogen_g_per_g_c_by_root_layer: []const f64,
    root_nonstructural_phosphorus_g_per_g_c_by_root_layer: []const f64,
    osmotic_potential_at_zero_total_megapascal: f64,
    soil_temperature_k_by_layer: []const f64,
    root_salt_concentration_mol_per_g_c_by_root_layer: []const f64,
};

pub const Outputs = struct {
    total_water_potential_mpa_by_root_layer: []f64,
    osmotic_water_potential_mpa_by_root_layer: []f64,
    turgor_water_potential_mpa_by_root_layer: []f64,
    water_uptake_m3_per_step_by_root_layer: []f64,
};

/// UPTAKE.F 1517--1530. Initializes runtime root axes/layers when the canopy
/// is too small for coupled water/energy exchange.
pub fn initialize(inputs: Inputs, outputs: Outputs) !void {
    try validate(inputs, outputs);
    for (0..inputs.root_axis_count) |root_axis| {
        for (inputs.first_soil_layer..inputs.end_soil_layer) |layer| {
            const index = root_axis * inputs.soil_layer_count + layer;
            const total_potential =
                inputs.adjusted_soil_total_water_potential_mpa_by_layer[layer];
            const absolute_potential = @abs(total_potential);
            const dry_matter_fraction =
                inputs.minimum_dry_matter_fraction_g_c_per_g +
                0.10 * absolute_potential /
                    (0.05 * absolute_potential + 2.0);
            const nonstructural_solute_concentration =
                inputs.root_nonstructural_carbon_g_per_g_c_by_root_layer[index] +
                inputs.root_nonstructural_nitrogen_g_per_g_c_by_root_layer[index] +
                inputs.root_nonstructural_phosphorus_g_per_g_c_by_root_layer[index];
            const solute_molar_mass =
                36.0 + 840.0 * @max(0, nonstructural_solute_concentration);
            const osmotic_potential =
                dry_matter_fraction /
                inputs.minimum_dry_matter_fraction_g_c_per_g *
                inputs.osmotic_potential_at_zero_total_megapascal -
                8.3143 * inputs.soil_temperature_k_by_layer[layer] *
                    dry_matter_fraction *
                    (nonstructural_solute_concentration / solute_molar_mass +
                        inputs.root_salt_concentration_mol_per_g_c_by_root_layer[index]);
            const turgor_potential = @max(
                0,
                total_potential - osmotic_potential,
            );
            inline for (.{ total_potential, osmotic_potential, turgor_potential }) |value|
                if (!std.math.isFinite(value))
                    return error.NonFiniteInactiveRootWaterStateResult;
            outputs.total_water_potential_mpa_by_root_layer[index] = total_potential;
            outputs.osmotic_water_potential_mpa_by_root_layer[index] = osmotic_potential;
            outputs.turgor_water_potential_mpa_by_root_layer[index] = turgor_potential;
            outputs.water_uptake_m3_per_step_by_root_layer[index] = 0;
        }
    }
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    const count = std.math.mul(
        usize,
        inputs.root_axis_count,
        inputs.soil_layer_count,
    ) catch return error.InactiveRootWaterStateDimensionMismatch;
    if (inputs.first_soil_layer > inputs.end_soil_layer or
        inputs.end_soil_layer > inputs.soil_layer_count or
        inputs.adjusted_soil_total_water_potential_mpa_by_layer.len != inputs.soil_layer_count or
        inputs.root_nonstructural_carbon_g_per_g_c_by_root_layer.len != count or
        inputs.root_nonstructural_nitrogen_g_per_g_c_by_root_layer.len != count or
        inputs.root_nonstructural_phosphorus_g_per_g_c_by_root_layer.len != count or
        inputs.soil_temperature_k_by_layer.len != inputs.soil_layer_count or
        inputs.root_salt_concentration_mol_per_g_c_by_root_layer.len != count or
        outputs.total_water_potential_mpa_by_root_layer.len != count or
        outputs.osmotic_water_potential_mpa_by_root_layer.len != count or
        outputs.turgor_water_potential_mpa_by_root_layer.len != count or
        outputs.water_uptake_m3_per_step_by_root_layer.len != count)
        return error.InactiveRootWaterStateDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_dry_matter_fraction_g_c_per_g) or
        inputs.minimum_dry_matter_fraction_g_c_per_g <= 0 or
        !std.math.isFinite(inputs.osmotic_potential_at_zero_total_megapascal))
        return error.InvalidInactiveRootWaterStateInput;
    for (inputs.adjusted_soil_total_water_potential_mpa_by_layer) |value|
        if (!std.math.isFinite(value))
            return error.InvalidInactiveRootWaterStateInput;
    for (inputs.soil_temperature_k_by_layer) |value|
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidInactiveRootWaterStateInput;
    for (0..count) |index| {
        inline for (.{
            inputs.root_nonstructural_carbon_g_per_g_c_by_root_layer[index],
            inputs.root_nonstructural_nitrogen_g_per_g_c_by_root_layer[index],
            inputs.root_nonstructural_phosphorus_g_per_g_c_by_root_layer[index],
            inputs.root_salt_concentration_mol_per_g_c_by_root_layer[index],
        }) |value|
            if (!std.math.isFinite(value))
                return error.InvalidInactiveRootWaterStateInput;
        if (inputs.root_salt_concentration_mol_per_g_c_by_root_layer[index] < 0)
            return error.InvalidInactiveRootWaterStateInput;
    }
}

test "inactive roots copy soil potential and zero uptake in source order" {
    const inputs = Inputs{
        .root_axis_count = 2,
        .soil_layer_count = 2,
        .first_soil_layer = 0,
        .end_soil_layer = 2,
        .adjusted_soil_total_water_potential_mpa_by_layer = &.{ -0.5, -1 },
        .minimum_dry_matter_fraction_g_c_per_g = 0.16,
        .root_nonstructural_carbon_g_per_g_c_by_root_layer = &.{ 0.1, 0.2, 0.3, 0.4 },
        .root_nonstructural_nitrogen_g_per_g_c_by_root_layer = &.{ 0.01, 0.02, 0.03, 0.04 },
        .root_nonstructural_phosphorus_g_per_g_c_by_root_layer = &.{ 0.001, 0.002, 0.003, 0.004 },
        .osmotic_potential_at_zero_total_megapascal = -1.5,
        .soil_temperature_k_by_layer = &.{ 290, 295 },
        .root_salt_concentration_mol_per_g_c_by_root_layer = &.{ 0.001, 0.002, 0.003, 0.004 },
    };
    var total = [_]f64{0} ** 4;
    var osmotic = [_]f64{0} ** 4;
    var turgor = [_]f64{0} ** 4;
    var uptake = [_]f64{9} ** 4;
    try initialize(inputs, .{
        .total_water_potential_mpa_by_root_layer = &total,
        .osmotic_water_potential_mpa_by_root_layer = &osmotic,
        .turgor_water_potential_mpa_by_root_layer = &turgor,
        .water_uptake_m3_per_step_by_root_layer = &uptake,
    });
    try std.testing.expectEqualSlices(f64, &.{ -0.5, -1, -0.5, -1 }, &total);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0, 0 }, &uptake);
    for (0..4) |index| {
        try std.testing.expect(osmotic[index] < total[index]);
        try std.testing.expectEqual(
            @max(0, total[index] - osmotic[index]),
            turgor[index],
        );
    }
}

test "inactive root runtime range leaves outside entries untouched" {
    const inputs = Inputs{
        .root_axis_count = 1,
        .soil_layer_count = 2,
        .first_soil_layer = 1,
        .end_soil_layer = 2,
        .adjusted_soil_total_water_potential_mpa_by_layer = &.{ -0.5, -1 },
        .minimum_dry_matter_fraction_g_c_per_g = 0.16,
        .root_nonstructural_carbon_g_per_g_c_by_root_layer = &.{ 0.1, 0.2 },
        .root_nonstructural_nitrogen_g_per_g_c_by_root_layer = &.{ 0.01, 0.02 },
        .root_nonstructural_phosphorus_g_per_g_c_by_root_layer = &.{ 0.001, 0.002 },
        .osmotic_potential_at_zero_total_megapascal = -1.5,
        .soil_temperature_k_by_layer = &.{ 290, 295 },
        .root_salt_concentration_mol_per_g_c_by_root_layer = &.{ 0.001, 0.002 },
    };
    var total = [_]f64{ 9, 9 };
    var osmotic = [_]f64{ 9, 9 };
    var turgor = [_]f64{ 9, 9 };
    var uptake = [_]f64{ 9, 9 };
    try initialize(inputs, .{
        .total_water_potential_mpa_by_root_layer = &total,
        .osmotic_water_potential_mpa_by_root_layer = &osmotic,
        .turgor_water_potential_mpa_by_root_layer = &turgor,
        .water_uptake_m3_per_step_by_root_layer = &uptake,
    });
    try std.testing.expectEqual(@as(f64, 9), total[0]);
    try std.testing.expectEqual(@as(f64, 0), uptake[1]);
}
