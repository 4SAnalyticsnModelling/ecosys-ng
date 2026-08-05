const std = @import("std");

pub const species_count = 41;
pub const phosphate_species_start = 33;
pub const phosphate_species_count = 8;
pub const Direction = enum { east_west, north_south, vertical };

pub const Inputs = struct {
    direction: Direction,
    /// Current volumetric snow heat capacity by layer (MJ K-1).
    snow_heat_capacity_megajoules_k: []const f64,
    minimum_snow_heat_capacity_megajoules_k: f64,
    /// Layer-major vertical snowpack fluxes `[layer][species]` (mol timestep-1).
    vertical_flux_mol_per_step: []const f64,
    /// Vertical fluxes at litter layer zero and the soil surface (mol timestep-1).
    litter_surface_flux_mol_per_step: []const f64,
    soil_surface_flux_mol_per_step: []const f64,
    /// Extra soil-surface BFB terms for species 33..40 (mol timestep-1).
    soil_surface_bubble_flux_mol_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 8819--9007.
/// Active snow layers receive their own downward flux minus either the next
/// active layer's flux or the litter/soil surface losses. The final eight
/// phosphate-family species also subtract the source BFB term at the surface.
pub fn accumulate(inputs: Inputs, totals_mol_per_step: []f64) !usize {
    if (inputs.direction != .vertical) return 0;
    const layer_count = inputs.snow_heat_capacity_megajoules_k.len;
    try validate(inputs, totals_mol_per_step, layer_count);

    var updated_layers: usize = 0;
    for (0..layer_count) |layer| {
        if (inputs.snow_heat_capacity_megajoules_k[layer] <= inputs.minimum_snow_heat_capacity_megajoules_k) continue;
        const lower_is_snow = layer + 1 < layer_count and
            inputs.snow_heat_capacity_megajoules_k[layer + 1] > inputs.minimum_snow_heat_capacity_megajoules_k;
        for (0..species_count) |species| {
            const index = layer * species_count + species;
            var result = totals_mol_per_step[index] + inputs.vertical_flux_mol_per_step[index];
            if (lower_is_snow) {
                result -= inputs.vertical_flux_mol_per_step[index + species_count];
            } else {
                result -= inputs.litter_surface_flux_mol_per_step[species];
                result -= inputs.soil_surface_flux_mol_per_step[species];
                if (species >= phosphate_species_start)
                    result -= inputs.soil_surface_bubble_flux_mol_per_step[species - phosphate_species_start];
            }
            if (!std.math.isFinite(result)) return error.NonFiniteSnowpackSoluteResult;
        }
    }

    for (0..layer_count) |layer| {
        if (inputs.snow_heat_capacity_megajoules_k[layer] <= inputs.minimum_snow_heat_capacity_megajoules_k) continue;
        const lower_is_snow = layer + 1 < layer_count and
            inputs.snow_heat_capacity_megajoules_k[layer + 1] > inputs.minimum_snow_heat_capacity_megajoules_k;
        for (0..species_count) |species| {
            const index = layer * species_count + species;
            totals_mol_per_step[index] += inputs.vertical_flux_mol_per_step[index];
            if (lower_is_snow) {
                totals_mol_per_step[index] -= inputs.vertical_flux_mol_per_step[index + species_count];
            } else {
                totals_mol_per_step[index] -= inputs.litter_surface_flux_mol_per_step[species];
                totals_mol_per_step[index] -= inputs.soil_surface_flux_mol_per_step[species];
                if (species >= phosphate_species_start)
                    totals_mol_per_step[index] -= inputs.soil_surface_bubble_flux_mol_per_step[species - phosphate_species_start];
            }
        }
        updated_layers += 1;
    }
    return updated_layers;
}

fn validate(inputs: Inputs, totals: []const f64, layer_count: usize) !void {
    if (!std.math.isFinite(inputs.minimum_snow_heat_capacity_megajoules_k))
        return error.NonFiniteSnowpackSoluteInput;
    const state_value_count = std.math.mul(usize, layer_count, species_count) catch
        return error.SnowpackSoluteDimensionMismatch;
    if (inputs.vertical_flux_mol_per_step.len != state_value_count or
        totals.len != state_value_count or
        inputs.litter_surface_flux_mol_per_step.len != species_count or
        inputs.soil_surface_flux_mol_per_step.len != species_count or
        inputs.soil_surface_bubble_flux_mol_per_step.len != phosphate_species_count)
        return error.SnowpackSoluteDimensionMismatch;
    const slices = [_][]const f64{
        inputs.snow_heat_capacity_megajoules_k,
        inputs.vertical_flux_mol_per_step,
        inputs.litter_surface_flux_mol_per_step,
        inputs.soil_surface_flux_mol_per_step,
        inputs.soil_surface_bubble_flux_mol_per_step,
        totals,
    };
    for (slices) |slice| for (slice) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSnowpackSoluteInput;
}

fn testInputs(heat: []const f64, vertical: []const f64, litter: []const f64, soil: []const f64, bubble: []const f64) Inputs {
    return .{
        .direction = .vertical,
        .snow_heat_capacity_megajoules_k = heat,
        .minimum_snow_heat_capacity_megajoules_k = 1,
        .vertical_flux_mol_per_step = vertical,
        .litter_surface_flux_mol_per_step = litter,
        .soil_surface_flux_mol_per_step = soil,
        .soil_surface_bubble_flux_mol_per_step = bubble,
    };
}

test "TRNSFRS active upper snow layer subtracts next active layer flux" {
    const heat = [_]f64{ 2, 2 };
    var vertical = [_]f64{0} ** (2 * species_count);
    @memset(vertical[0..species_count], 5);
    @memset(vertical[species_count..], 2);
    const litter = [_]f64{1} ** species_count;
    const soil = [_]f64{1} ** species_count;
    const bubble = [_]f64{1} ** phosphate_species_count;
    var totals = [_]f64{0} ** (2 * species_count);
    try std.testing.expectEqual(@as(usize, 2), try accumulate(testInputs(&heat, &vertical, &litter, &soil, &bubble), &totals));
    try std.testing.expectEqual(@as(f64, 3), totals[0]);
    try std.testing.expectEqual(@as(f64, 0), totals[species_count]);
    try std.testing.expectEqual(@as(f64, -1), totals[species_count + phosphate_species_start]);
}

test "inactive lower layer makes upper layer use surface boundary" {
    const heat = [_]f64{ 2, 1 };
    const vertical = [_]f64{5} ** (2 * species_count);
    const litter = [_]f64{1} ** species_count;
    const soil = [_]f64{2} ** species_count;
    const bubble = [_]f64{4} ** phosphate_species_count;
    var totals = [_]f64{0} ** (2 * species_count);
    try std.testing.expectEqual(@as(usize, 1), try accumulate(testInputs(&heat, &vertical, &litter, &soil, &bubble), &totals));
    try std.testing.expectEqual(@as(f64, 2), totals[0]);
    try std.testing.expectEqual(@as(f64, -2), totals[phosphate_species_start]);
    try std.testing.expectEqual(@as(f64, 0), totals[species_count]);
}

test "heat capacity equal to minimum is not snowpack" {
    const heat = [_]f64{1};
    const vertical = [_]f64{5} ** species_count;
    const litter = [_]f64{1} ** species_count;
    const soil = [_]f64{2} ** species_count;
    const bubble = [_]f64{4} ** phosphate_species_count;
    var totals = [_]f64{7} ** species_count;
    try std.testing.expectEqual(@as(usize, 0), try accumulate(testInputs(&heat, &vertical, &litter, &soil, &bubble), &totals));
    try std.testing.expectEqual(@as(f64, 7), totals[0]);
}

test "nonvertical direction leaves snowpack totals unchanged" {
    const heat = [_]f64{};
    const vertical = [_]f64{};
    const litter = [_]f64{};
    const soil = [_]f64{};
    const bubble = [_]f64{};
    var inputs = testInputs(&heat, &vertical, &litter, &soil, &bubble);
    inputs.direction = .east_west;
    var totals = [_]f64{};
    try std.testing.expectEqual(@as(usize, 0), try accumulate(inputs, &totals));
}

test "runtime snowpack topology mismatch fails atomically" {
    const heat = [_]f64{2};
    const short_vertical = [_]f64{5} ** (species_count - 1);
    const litter = [_]f64{1} ** species_count;
    const soil = [_]f64{2} ** species_count;
    const bubble = [_]f64{4} ** phosphate_species_count;
    var totals = [_]f64{7} ** species_count;
    try std.testing.expectError(error.SnowpackSoluteDimensionMismatch, accumulate(testInputs(&heat, &short_vertical, &litter, &soil, &bubble), &totals));
    try std.testing.expectEqual(@as(f64, 7), totals[0]);
}

test "late surface overflow leaves all snow totals atomic" {
    const heat = [_]f64{ 2, 1 };
    const vertical = [_]f64{0} ** (2 * species_count);
    const litter = [_]f64{0} ** species_count;
    var soil = [_]f64{0} ** species_count;
    soil[species_count - 1] = -std.math.floatMax(f64);
    const bubble = [_]f64{0} ** phosphate_species_count;
    var totals = [_]f64{0} ** (2 * species_count);
    totals[species_count - 1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSnowpackSoluteResult, accumulate(testInputs(&heat, &vertical, &litter, &soil, &bubble), &totals));
    try std.testing.expectEqual(@as(f64, 0), totals[0]);
}
