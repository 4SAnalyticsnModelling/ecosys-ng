const std = @import("std");

pub const Inputs = struct {
    leaf_protein_surface_density_g_per_m2: f64,
    rubisco_leaf_protein_fraction: f64,
    chlorophyll_leaf_protein_fraction: f64,
    rubisco_carboxylation_umol_per_g_s_25c: f64,
    carboxylation_temperature_factor: f64,
    rubisco_oxygenation_umol_per_g_s_25c: f64,
    oxygenation_temperature_factor: f64,
    dissolved_oxygen_umol_per_l: f64,
    rubisco_co2_half_saturation_umol_per_l: f64,
    rubisco_o2_half_saturation_umol_per_l: f64,
    intercellular_co2_umol_per_l: f64,
    rubisco_co2_half_saturation_with_o2_umol_per_l: f64,
    chlorophyll_electron_transport_umol_per_g_s_25c: f64,
    electron_transport_temperature_factor: f64,
    electron_requirement_umol_e_per_umol_co2: f64,
};

pub const Result = struct {
    rubisco_surface_density_g_per_m2: f64,
    chlorophyll_surface_density_g_per_m2: f64,
    co2_unlimited_carboxylation_umol_per_m2_s: f64,
    oxygenation_umol_per_m2_s: f64,
    co2_compensation_umol_per_l: f64,
    co2_limited_carboxylation_umol_per_m2_s: f64,
    light_saturated_electron_transport_umol_per_m2_s: f64,
    carboxylation_umol_co2_per_umol_electron: f64,
};

/// STOMATE.F lines 477--519 standalone mesophyll C3 capacity in exact
/// source assignment order.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var result: Result = undefined;
    result.rubisco_surface_density_g_per_m2 =
        inputs.rubisco_leaf_protein_fraction *
        inputs.leaf_protein_surface_density_g_per_m2;
    result.chlorophyll_surface_density_g_per_m2 =
        inputs.chlorophyll_leaf_protein_fraction *
        inputs.leaf_protein_surface_density_g_per_m2;
    result.co2_unlimited_carboxylation_umol_per_m2_s =
        inputs.rubisco_carboxylation_umol_per_g_s_25c *
        inputs.carboxylation_temperature_factor *
        result.rubisco_surface_density_g_per_m2;
    if (result.co2_unlimited_carboxylation_umol_per_m2_s <= 0)
        return error.ZeroMesophyllRubiscoCapacity;
    result.oxygenation_umol_per_m2_s =
        inputs.rubisco_oxygenation_umol_per_g_s_25c *
        inputs.oxygenation_temperature_factor *
        result.rubisco_surface_density_g_per_m2;
    result.co2_compensation_umol_per_l =
        0.5 *
        inputs.dissolved_oxygen_umol_per_l *
        result.oxygenation_umol_per_m2_s *
        inputs.rubisco_co2_half_saturation_umol_per_l /
        (result.co2_unlimited_carboxylation_umol_per_m2_s *
            inputs.rubisco_o2_half_saturation_umol_per_l);
    result.co2_limited_carboxylation_umol_per_m2_s = @max(
        0,
        result.co2_unlimited_carboxylation_umol_per_m2_s *
            (inputs.intercellular_co2_umol_per_l -
                result.co2_compensation_umol_per_l) /
            (inputs.intercellular_co2_umol_per_l +
                inputs.rubisco_co2_half_saturation_with_o2_umol_per_l),
    );
    result.light_saturated_electron_transport_umol_per_m2_s =
        inputs.chlorophyll_electron_transport_umol_per_g_s_25c *
        inputs.electron_transport_temperature_factor *
        result.chlorophyll_surface_density_g_per_m2;
    result.carboxylation_umol_co2_per_umol_electron = @max(
        0,
        (inputs.intercellular_co2_umol_per_l -
            result.co2_compensation_umol_per_l) /
            (inputs.electron_requirement_umol_e_per_umol_co2 *
                inputs.intercellular_co2_umol_per_l +
                10.5 * result.co2_compensation_umol_per_l),
    );
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteMesophyllC3CapacityResult;
    return result;
}

pub fn computeRuntimeNodes(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.MesophyllC3CapacityDimensionMismatch;
    for (inputs, scratch) |node_inputs, *candidate|
        candidate.* = try compute(node_inputs);
    @memcpy(destination, scratch);
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteMesophyllC3CapacityInput;
    if (inputs.leaf_protein_surface_density_g_per_m2 <= 0 or
        inputs.rubisco_leaf_protein_fraction < 0 or
        inputs.chlorophyll_leaf_protein_fraction < 0 or
        inputs.rubisco_carboxylation_umol_per_g_s_25c <= 0 or
        inputs.carboxylation_temperature_factor < 0 or
        inputs.rubisco_oxygenation_umol_per_g_s_25c < 0 or
        inputs.oxygenation_temperature_factor < 0 or
        inputs.dissolved_oxygen_umol_per_l < 0 or
        inputs.rubisco_co2_half_saturation_umol_per_l <= 0 or
        inputs.rubisco_o2_half_saturation_umol_per_l <= 0 or
        inputs.intercellular_co2_umol_per_l < 0 or
        inputs.rubisco_co2_half_saturation_with_o2_umol_per_l <= 0 or
        inputs.chlorophyll_electron_transport_umol_per_g_s_25c < 0 or
        inputs.electron_transport_temperature_factor < 0 or
        inputs.electron_requirement_umol_e_per_umol_co2 <= 0)
        return error.InvalidMesophyllC3CapacityInput;
}

fn fixture() Inputs {
    return .{
        .leaf_protein_surface_density_g_per_m2 = 2,
        .rubisco_leaf_protein_fraction = 0.2,
        .chlorophyll_leaf_protein_fraction = 0.1,
        .rubisco_carboxylation_umol_per_g_s_25c = 75,
        .carboxylation_temperature_factor = 0.8,
        .rubisco_oxygenation_umol_per_g_s_25c = 20,
        .oxygenation_temperature_factor = 0.7,
        .dissolved_oxygen_umol_per_l = 250,
        .rubisco_co2_half_saturation_umol_per_l = 30,
        .rubisco_o2_half_saturation_umol_per_l = 300,
        .intercellular_co2_umol_per_l = 25,
        .rubisco_co2_half_saturation_with_o2_umol_per_l = 55,
        .chlorophyll_electron_transport_umol_per_g_s_25c = 120,
        .electron_transport_temperature_factor = 0.9,
        .electron_requirement_umol_e_per_umol_co2 = 4.5,
    };
}

test "mesophyll C3 capacity publishes every source field" {
    const result = try compute(fixture());
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result.rubisco_surface_density_g_per_m2, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.chlorophyll_surface_density_g_per_m2, 1e-15);
    try std.testing.expect(result.co2_compensation_umol_per_l > 0);
    try std.testing.expect(result.light_saturated_electron_transport_umol_per_m2_s > 0);
}

test "mesophyll shared outputs reconcile with existing C3 capacity" {
    const stomatal = @import("../energy/stomatal_resistance.zig");
    const inputs = fixture();
    const gas: stomatal.GasEnvironment = .{
        .air_amount_mol_per_m3 = 1,
        .intercellular_co2_umol_per_mol = 1,
        .co2_solubility_umol_per_l_per_umol_per_mol = 1,
        .o2_solubility_umol_per_l_per_umol_per_mol = 1,
        .dissolved_co2_umol_per_l = inputs.intercellular_co2_umol_per_l,
        .dissolved_o2_umol_per_l = inputs.dissolved_oxygen_umol_per_l,
        .atmospheric_to_intercellular_co2_umol_per_m3 = 1,
        .rubisco_carboxylation_temperature_factor = inputs.carboxylation_temperature_factor,
        .rubisco_oxygenation_temperature_factor = inputs.oxygenation_temperature_factor,
        .electron_transport_temperature_factor = inputs.electron_transport_temperature_factor,
        .rubisco_co2_half_saturation_umol_per_l = inputs.rubisco_co2_half_saturation_umol_per_l,
        .rubisco_co2_half_saturation_with_o2_umol_per_l = inputs.rubisco_co2_half_saturation_with_o2_umol_per_l,
    };
    const existing = try stomatal.c3Capacity(
        inputs.leaf_protein_surface_density_g_per_m2,
        inputs.rubisco_leaf_protein_fraction,
        inputs.chlorophyll_leaf_protein_fraction,
        inputs.rubisco_carboxylation_umol_per_g_s_25c,
        inputs.rubisco_oxygenation_umol_per_g_s_25c,
        inputs.chlorophyll_electron_transport_umol_per_g_s_25c,
        gas,
        inputs.rubisco_o2_half_saturation_umol_per_l,
        inputs.intercellular_co2_umol_per_l,
    );
    const result = try compute(inputs);
    try std.testing.expectEqual(existing.co2_unlimited_carboxylation_umol_per_m2_s, result.co2_unlimited_carboxylation_umol_per_m2_s);
    try std.testing.expectEqual(existing.oxygenation_umol_per_m2_s, result.oxygenation_umol_per_m2_s);
    try std.testing.expectEqual(existing.co2_compensation_umol_per_l, result.co2_compensation_umol_per_l);
    try std.testing.expectEqual(existing.co2_limited_carboxylation_umol_per_m2_s, result.co2_limited_carboxylation_umol_per_m2_s);
    try std.testing.expectEqual(existing.light_saturated_electron_transport_umol_per_m2_s, result.light_saturated_electron_transport_umol_per_m2_s);
    try std.testing.expectEqual(existing.carboxylation_umol_co2_per_umol_electron, result.carboxylation_umol_co2_per_umol_electron);
}

test "later invalid mesophyll node leaves destination unchanged" {
    var invalid = fixture();
    invalid.electron_requirement_umol_e_per_umol_co2 = 0;
    const inputs = [_]Inputs{ fixture(), invalid };
    var scratch: [2]Result = undefined;
    var destination: [2]Result = undefined;
    destination[0] = try compute(fixture());
    destination[1] = destination[0];
    destination[0].co2_compensation_umol_per_l = 41;
    destination[1].co2_compensation_umol_per_l = 42;
    try std.testing.expectError(
        error.InvalidMesophyllC3CapacityInput,
        computeRuntimeNodes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].co2_compensation_umol_per_l);
    try std.testing.expectEqual(@as(f64, 42), destination[1].co2_compensation_umol_per_l);
}
