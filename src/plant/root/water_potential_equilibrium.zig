const std = @import("std");

pub const Inputs = struct {
    root_axis_count: usize,
    soil_layer_count: usize,
    first_soil_layer: usize,
    end_soil_layer: usize,
    rooted_by_root_layer: []const bool,
    adjusted_soil_total_water_potential_mpa_by_layer: []const f64,
    canopy_total_water_potential_megapascal: f64,
    root_resistance_mpa_h_per_m_by_root_layer: []const f64,
    soil_resistance_mpa_h_per_m_by_root_layer: []const f64,
    soil_plus_root_resistance_mpa_h_per_m_by_root_layer: []const f64,
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
};

/// UPTAKE.F 1439--1474. Publishes converged root total, osmotic, and turgor
/// potentials over runtime root axes and the runtime active layer interval.
pub fn calculate(inputs: Inputs, outputs: Outputs) !void {
    try validate(inputs, outputs);
    for (0..inputs.root_axis_count) |root_axis| {
        for (inputs.first_soil_layer..inputs.end_soil_layer) |layer| {
            const index = root_axis * inputs.soil_layer_count + layer;
            const total_potential = if (inputs.rooted_by_root_layer[index])
                @min(
                    0,
                    (inputs.adjusted_soil_total_water_potential_mpa_by_layer[layer] *
                        inputs.root_resistance_mpa_h_per_m_by_root_layer[index] +
                        inputs.canopy_total_water_potential_megapascal *
                            inputs.soil_resistance_mpa_h_per_m_by_root_layer[index]) /
                        inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer[index],
                )
            else
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
                    return error.NonFiniteRootWaterPotentialResult;
            outputs.total_water_potential_mpa_by_root_layer[index] = total_potential;
            outputs.osmotic_water_potential_mpa_by_root_layer[index] = osmotic_potential;
            outputs.turgor_water_potential_mpa_by_root_layer[index] = turgor_potential;
        }
    }
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    const count = std.math.mul(
        usize,
        inputs.root_axis_count,
        inputs.soil_layer_count,
    ) catch return error.RootWaterPotentialDimensionMismatch;
    if (inputs.first_soil_layer > inputs.end_soil_layer or
        inputs.end_soil_layer > inputs.soil_layer_count or
        inputs.rooted_by_root_layer.len != count or
        inputs.adjusted_soil_total_water_potential_mpa_by_layer.len != inputs.soil_layer_count or
        inputs.root_resistance_mpa_h_per_m_by_root_layer.len != count or
        inputs.soil_resistance_mpa_h_per_m_by_root_layer.len != count or
        inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer.len != count or
        inputs.root_nonstructural_carbon_g_per_g_c_by_root_layer.len != count or
        inputs.root_nonstructural_nitrogen_g_per_g_c_by_root_layer.len != count or
        inputs.root_nonstructural_phosphorus_g_per_g_c_by_root_layer.len != count or
        inputs.soil_temperature_k_by_layer.len != inputs.soil_layer_count or
        inputs.root_salt_concentration_mol_per_g_c_by_root_layer.len != count or
        outputs.total_water_potential_mpa_by_root_layer.len != count or
        outputs.osmotic_water_potential_mpa_by_root_layer.len != count or
        outputs.turgor_water_potential_mpa_by_root_layer.len != count)
        return error.RootWaterPotentialDimensionMismatch;
    if (!std.math.isFinite(inputs.canopy_total_water_potential_megapascal) or
        !std.math.isFinite(inputs.minimum_dry_matter_fraction_g_c_per_g) or
        inputs.minimum_dry_matter_fraction_g_c_per_g <= 0 or
        !std.math.isFinite(inputs.osmotic_potential_at_zero_total_megapascal))
        return error.InvalidRootWaterPotentialInput;
    for (inputs.adjusted_soil_total_water_potential_mpa_by_layer) |value|
        if (!std.math.isFinite(value))
            return error.InvalidRootWaterPotentialInput;
    for (inputs.soil_temperature_k_by_layer) |value|
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidRootWaterPotentialInput;
    for (0..count) |index| {
        inline for (.{
            inputs.root_resistance_mpa_h_per_m_by_root_layer[index],
            inputs.soil_resistance_mpa_h_per_m_by_root_layer[index],
            inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer[index],
            inputs.root_nonstructural_carbon_g_per_g_c_by_root_layer[index],
            inputs.root_nonstructural_nitrogen_g_per_g_c_by_root_layer[index],
            inputs.root_nonstructural_phosphorus_g_per_g_c_by_root_layer[index],
            inputs.root_salt_concentration_mol_per_g_c_by_root_layer[index],
        }) |value|
            if (!std.math.isFinite(value))
                return error.InvalidRootWaterPotentialInput;
        if (inputs.root_resistance_mpa_h_per_m_by_root_layer[index] < 0 or
            inputs.soil_resistance_mpa_h_per_m_by_root_layer[index] < 0 or
            inputs.root_salt_concentration_mol_per_g_c_by_root_layer[index] < 0 or
            (inputs.rooted_by_root_layer[index] and
                inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer[index] <= 0))
            return error.InvalidRootWaterPotentialInput;
    }
}

fn sourceInputs() Inputs {
    return .{
        .root_axis_count = 1,
        .soil_layer_count = 2,
        .first_soil_layer = 0,
        .end_soil_layer = 2,
        .rooted_by_root_layer = &.{ true, false },
        .adjusted_soil_total_water_potential_mpa_by_layer = &.{ -0.5, -1.2 },
        .canopy_total_water_potential_megapascal = -2,
        .root_resistance_mpa_h_per_m_by_root_layer = &.{ 3, 3 },
        .soil_resistance_mpa_h_per_m_by_root_layer = &.{ 1, 1 },
        .soil_plus_root_resistance_mpa_h_per_m_by_root_layer = &.{ 4, 4 },
        .minimum_dry_matter_fraction_g_c_per_g = 0.16,
        .root_nonstructural_carbon_g_per_g_c_by_root_layer = &.{ 0.1, 0.2 },
        .root_nonstructural_nitrogen_g_per_g_c_by_root_layer = &.{ 0.01, 0.02 },
        .root_nonstructural_phosphorus_g_per_g_c_by_root_layer = &.{ 0.001, 0.002 },
        .osmotic_potential_at_zero_total_megapascal = -1.5,
        .soil_temperature_k_by_layer = &.{ 290, 295 },
        .root_salt_concentration_mol_per_g_c_by_root_layer = &.{ 0.001, 0.002 },
    };
}

test "UPTAKE converged rooted and unrooted total potentials preserve branches" {
    const inputs = sourceInputs();
    var total = [_]f64{ 9, 9 };
    var osmotic = [_]f64{ 9, 9 };
    var turgor = [_]f64{ 9, 9 };
    try calculate(inputs, .{
        .total_water_potential_mpa_by_root_layer = &total,
        .osmotic_water_potential_mpa_by_root_layer = &osmotic,
        .turgor_water_potential_mpa_by_root_layer = &turgor,
    });
    try std.testing.expectEqual(
        @as(f64, (-0.5 * 3.0 + -2.0 * 1.0) / 4.0),
        total[0],
    );
    try std.testing.expectEqual(@as(f64, -1.2), total[1]);
    try std.testing.expect(osmotic[0] < total[0]);
    try std.testing.expect(osmotic[1] < total[1]);
    try std.testing.expectEqual(@max(0, total[0] - osmotic[0]), turgor[0]);
    try std.testing.expectEqual(@max(0, total[1] - osmotic[1]), turgor[1]);
}

test "runtime active layer range leaves outside output untouched" {
    var inputs = sourceInputs();
    inputs.first_soil_layer = 1;
    var total = [_]f64{ 9, 9 };
    var osmotic = [_]f64{ 9, 9 };
    var turgor = [_]f64{ 9, 9 };
    try calculate(inputs, .{
        .total_water_potential_mpa_by_root_layer = &total,
        .osmotic_water_potential_mpa_by_root_layer = &osmotic,
        .turgor_water_potential_mpa_by_root_layer = &turgor,
    });
    try std.testing.expectEqual(@as(f64, 9), total[0]);
    try std.testing.expectEqual(@as(f64, -1.2), total[1]);
}

test "rooted zero combined resistance fails explicitly" {
    var combined = [_]f64{ 0, 4 };
    var inputs = sourceInputs();
    inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer = &combined;
    var total = [_]f64{0} ** 2;
    var osmotic = [_]f64{0} ** 2;
    var turgor = [_]f64{0} ** 2;
    try std.testing.expectError(
        error.InvalidRootWaterPotentialInput,
        calculate(inputs, .{
            .total_water_potential_mpa_by_root_layer = &total,
            .osmotic_water_potential_mpa_by_root_layer = &osmotic,
            .turgor_water_potential_mpa_by_root_layer = &turgor,
        }),
    );
}
