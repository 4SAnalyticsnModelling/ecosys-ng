const std = @import("std");

pub const Dimensions = struct {
    layer_count: usize,
    inclination_count: usize,
    azimuth_count: usize,
};

pub const Inputs = struct {
    dimensions: Dimensions,
    leaf_area_m2_by_layer: []const f64,
    unself_shaded_leaf_area_m2_by_inclination_layer: []const f64,
    direct_par_umol_per_m2_s: []const f64,
    diffuse_par_umol_per_m2_s: []const f64,
    direct_transmission_fraction_by_layer_boundary: []const f64,
    diffuse_transmission_fraction_by_layer_boundary: []const f64,
    negligible_leaf_area_m2: f64,
    quantum_efficiency_umol_e_per_umol_par: f64,
    light_saturated_electron_transport_umol_per_m2_s: f64,
    electron_transport_curvature_product: f64,
    electron_transport_curvature_divisor: f64,
    carboxylation_umol_co2_per_umol_electron: f64,
    co2_limited_carboxylation_umol_per_m2_s: f64,
    feedback_fraction: f64,
};

/// STOMATE lines 303--408. Runtime layers preserve the source's descending
/// traversal; inclination and azimuth samples ascend, with direct preceding
/// diffuse radiation at every sample.
pub fn integrate(inputs: Inputs) !f64 {
    try validate(inputs);

    var total_carboxylation_umol_per_s: f64 = 0;
    var layer = inputs.dimensions.layer_count;
    while (layer > 0) {
        layer -= 1;
        if (inputs.leaf_area_m2_by_layer[layer] >
            inputs.negligible_leaf_area_m2)
        {
            for (0..inputs.dimensions.inclination_count) |inclination| {
                const surface_index =
                    inclination * inputs.dimensions.layer_count + layer;
                const unself_shaded_leaf_area_m2 =
                    inputs.unself_shaded_leaf_area_m2_by_inclination_layer[surface_index];
                if (unself_shaded_leaf_area_m2 >
                    inputs.negligible_leaf_area_m2)
                {
                    for (0..inputs.dimensions.azimuth_count) |azimuth| {
                        const sample_index =
                            (inclination * inputs.dimensions.azimuth_count +
                                azimuth) *
                            inputs.dimensions.layer_count +
                            layer;
                        const direct_par =
                            inputs.direct_par_umol_per_m2_s[sample_index];
                        if (direct_par > 0) {
                            const rate = try limitedCarboxylation(inputs, direct_par);
                            total_carboxylation_umol_per_s +=
                                rate *
                                unself_shaded_leaf_area_m2 *
                                inputs.direct_transmission_fraction_by_layer_boundary[layer + 1];
                            if (!std.math.isFinite(total_carboxylation_umol_per_s))
                                return error.NonFiniteCanopyC4CarboxylationTotal;
                        }

                        const diffuse_par =
                            inputs.diffuse_par_umol_per_m2_s[sample_index];
                        if (diffuse_par > 0) {
                            const rate = try limitedCarboxylation(inputs, diffuse_par);
                            total_carboxylation_umol_per_s +=
                                rate *
                                unself_shaded_leaf_area_m2 *
                                inputs.diffuse_transmission_fraction_by_layer_boundary[layer + 1];
                            if (!std.math.isFinite(total_carboxylation_umol_per_s))
                                return error.NonFiniteCanopyC4CarboxylationTotal;
                        }
                    }
                }
            }
        }
    }
    return total_carboxylation_umol_per_s;
}

fn limitedCarboxylation(inputs: Inputs, par_umol_per_m2_s: f64) !f64 {
    const light_electrons =
        inputs.quantum_efficiency_umol_e_per_umol_par * par_umol_per_m2_s;
    const combined_electrons =
        light_electrons + inputs.light_saturated_electron_transport_umol_per_m2_s;
    const discriminant =
        combined_electrons * combined_electrons -
        inputs.electron_transport_curvature_product *
            light_electrons *
            inputs.light_saturated_electron_transport_umol_per_m2_s;
    if (!std.math.isFinite(discriminant) or discriminant < 0)
        return error.InvalidCanopyC4ElectronTransportDiscriminant;
    const light_limited_electron_transport =
        (combined_electrons - @sqrt(discriminant)) /
        inputs.electron_transport_curvature_divisor;
    const light_limited_carboxylation =
        light_limited_electron_transport *
        inputs.carboxylation_umol_co2_per_umol_electron;
    const rate = @min(
        inputs.co2_limited_carboxylation_umol_per_m2_s,
        light_limited_carboxylation,
    ) * inputs.feedback_fraction;
    if (!std.math.isFinite(rate) or rate < 0)
        return error.InvalidCanopyC4CarboxylationRate;
    return rate;
}

fn validate(inputs: Inputs) !void {
    const dimensions = inputs.dimensions;
    const surfaces = try std.math.mul(
        usize,
        dimensions.inclination_count,
        dimensions.layer_count,
    );
    const orientations = try std.math.mul(
        usize,
        dimensions.inclination_count,
        dimensions.azimuth_count,
    );
    const samples = try std.math.mul(usize, orientations, dimensions.layer_count);
    if (inputs.leaf_area_m2_by_layer.len != dimensions.layer_count or
        inputs.unself_shaded_leaf_area_m2_by_inclination_layer.len != surfaces or
        inputs.direct_par_umol_per_m2_s.len != samples or
        inputs.diffuse_par_umol_per_m2_s.len != samples or
        inputs.direct_transmission_fraction_by_layer_boundary.len !=
            dimensions.layer_count + 1 or
        inputs.diffuse_transmission_fraction_by_layer_boundary.len !=
            dimensions.layer_count + 1)
        return error.CanopyC4RadiationDimensionMismatch;

    const scalars = [_]f64{
        inputs.negligible_leaf_area_m2,
        inputs.quantum_efficiency_umol_e_per_umol_par,
        inputs.light_saturated_electron_transport_umol_per_m2_s,
        inputs.electron_transport_curvature_product,
        inputs.electron_transport_curvature_divisor,
        inputs.carboxylation_umol_co2_per_umol_electron,
        inputs.co2_limited_carboxylation_umol_per_m2_s,
        inputs.feedback_fraction,
    };
    for (scalars) |value| if (!std.math.isFinite(value))
        return error.NonFiniteCanopyC4RadiationInput;
    if (inputs.negligible_leaf_area_m2 < 0 or
        inputs.quantum_efficiency_umol_e_per_umol_par < 0 or
        inputs.light_saturated_electron_transport_umol_per_m2_s < 0 or
        inputs.electron_transport_curvature_product < 0 or
        inputs.electron_transport_curvature_divisor <= 0 or
        inputs.carboxylation_umol_co2_per_umol_electron < 0 or
        inputs.co2_limited_carboxylation_umol_per_m2_s < 0 or
        inputs.feedback_fraction < 0 or inputs.feedback_fraction > 1)
        return error.InvalidCanopyC4RadiationInput;

    for (inputs.leaf_area_m2_by_layer) |value| try validateNonnegative(value);
    for (inputs.unself_shaded_leaf_area_m2_by_inclination_layer) |value|
        try validateNonnegative(value);
    for (inputs.direct_par_umol_per_m2_s) |value| try validateNonnegative(value);
    for (inputs.diffuse_par_umol_per_m2_s) |value| try validateNonnegative(value);
    for (inputs.direct_transmission_fraction_by_layer_boundary) |value|
        try validateFraction(value);
    for (inputs.diffuse_transmission_fraction_by_layer_boundary) |value|
        try validateFraction(value);
}

fn validateNonnegative(value: f64) !void {
    if (!std.math.isFinite(value)) return error.NonFiniteCanopyC4RadiationInput;
    if (value < 0) return error.InvalidCanopyC4RadiationInput;
}

fn validateFraction(value: f64) !void {
    if (!std.math.isFinite(value)) return error.NonFiniteCanopyC4RadiationInput;
    if (value < 0 or value > 1) return error.InvalidCanopyC4RadiationInput;
}

fn fixture() Inputs {
    return .{
        .dimensions = .{ .layer_count = 2, .inclination_count = 1, .azimuth_count = 2 },
        .leaf_area_m2_by_layer = &.{ 1, 1 },
        .unself_shaded_leaf_area_m2_by_inclination_layer = &.{ 2, 3 },
        .direct_par_umol_per_m2_s = &.{ 10, 20, 30, 40 },
        .diffuse_par_umol_per_m2_s = &.{ 5, 10, 15, 20 },
        .direct_transmission_fraction_by_layer_boundary = &.{ 0, 0.5, 0.25 },
        .diffuse_transmission_fraction_by_layer_boundary = &.{ 0, 0.2, 0.1 },
        .negligible_leaf_area_m2 = 1.0e-12,
        .quantum_efficiency_umol_e_per_umol_par = 1,
        .light_saturated_electron_transport_umol_per_m2_s = 100,
        .electron_transport_curvature_product = 4,
        .electron_transport_curvature_divisor = 2,
        .carboxylation_umol_co2_per_umol_electron = 1,
        .co2_limited_carboxylation_umol_per_m2_s = 1.0e6,
        .feedback_fraction = 1,
    };
}

test "C4 radiation integration preserves descending layers and direct diffuse gates" {
    const inputs = fixture();
    // CURV4=4 and CURV2=2 make ETLF4=min(PARX, ETGR4). Source order is layer 2,
    // azimuths 1..2, then layer 1, with direct before diffuse.
    const expected =
        20 * 3 * 0.25 + 10 * 3 * 0.1 +
        40 * 3 * 0.25 + 20 * 3 * 0.1 +
        10 * 2 * 0.5 + 5 * 2 * 0.2 +
        30 * 2 * 0.5 + 15 * 2 * 0.2;
    try std.testing.expectEqual(@as(f64, expected), try integrate(inputs));
}

test "C4 radiation gates inactive layer surface and zero PAR" {
    var inputs = fixture();
    inputs.leaf_area_m2_by_layer = &.{ 0, 1 };
    inputs.unself_shaded_leaf_area_m2_by_inclination_layer = &.{ 2, 0 };
    try std.testing.expectEqual(@as(f64, 0), try integrate(inputs));
}

test "C4 radiation rejects malformed runtime topology and invalid samples" {
    var inputs = fixture();
    inputs.direct_par_umol_per_m2_s = &.{ 1, 2 };
    try std.testing.expectError(
        error.CanopyC4RadiationDimensionMismatch,
        integrate(inputs),
    );
    inputs = fixture();
    inputs.diffuse_par_umol_per_m2_s = &.{ 5, 10, -1, 20 };
    try std.testing.expectError(error.InvalidCanopyC4RadiationInput, integrate(inputs));
}
