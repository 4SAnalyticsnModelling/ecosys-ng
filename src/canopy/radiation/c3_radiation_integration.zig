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
    nutrient_feedback_fraction: f64,
};

/// STOMATE.F lines 526--623. Runtime layers descend; inclination and
/// azimuth classes ascend; direct precedes diffuse at every sample.
pub fn integrate(inputs: Inputs) !f64 {
    try validate(inputs);
    var total_umol_per_s: f64 = 0;
    var layer = inputs.dimensions.layer_count;
    while (layer > 0) {
        layer -= 1;
        if (inputs.leaf_area_m2_by_layer[layer] >
            inputs.negligible_leaf_area_m2)
        {
            for (0..inputs.dimensions.inclination_count) |inclination| {
                const surface_index =
                    inclination * inputs.dimensions.layer_count + layer;
                const surface_m2 =
                    inputs.unself_shaded_leaf_area_m2_by_inclination_layer[surface_index];
                if (surface_m2 > inputs.negligible_leaf_area_m2) {
                    for (0..inputs.dimensions.azimuth_count) |azimuth| {
                        const sample_index =
                            (inclination * inputs.dimensions.azimuth_count +
                                azimuth) *
                            inputs.dimensions.layer_count +
                            layer;
                        const direct_par =
                            inputs.direct_par_umol_per_m2_s[sample_index];
                        if (direct_par > 0) {
                            total_umol_per_s +=
                                try limitedRate(inputs, direct_par) *
                                surface_m2 *
                                inputs.direct_transmission_fraction_by_layer_boundary[layer + 1];
                            if (!std.math.isFinite(total_umol_per_s))
                                return error.NonFiniteCanopyC3CarboxylationTotal;
                        }
                        const diffuse_par =
                            inputs.diffuse_par_umol_per_m2_s[sample_index];
                        if (diffuse_par > 0) {
                            total_umol_per_s +=
                                try limitedRate(inputs, diffuse_par) *
                                surface_m2 *
                                inputs.diffuse_transmission_fraction_by_layer_boundary[layer + 1];
                            if (!std.math.isFinite(total_umol_per_s))
                                return error.NonFiniteCanopyC3CarboxylationTotal;
                        }
                    }
                }
            }
        }
    }
    return total_umol_per_s;
}

fn limitedRate(inputs: Inputs, par_umol_per_m2_s: f64) !f64 {
    const par_electrons =
        inputs.quantum_efficiency_umol_e_per_umol_par * par_umol_per_m2_s;
    const combined =
        par_electrons + inputs.light_saturated_electron_transport_umol_per_m2_s;
    const discriminant =
        combined * combined -
        inputs.electron_transport_curvature_product *
            par_electrons *
            inputs.light_saturated_electron_transport_umol_per_m2_s;
    if (!std.math.isFinite(discriminant) or discriminant < 0)
        return error.InvalidCanopyC3ElectronTransportDiscriminant;
    const electron_transport =
        (combined - @sqrt(discriminant)) /
        inputs.electron_transport_curvature_divisor;
    const light_limited =
        electron_transport *
        inputs.carboxylation_umol_co2_per_umol_electron;
    const rate = @min(
        inputs.co2_limited_carboxylation_umol_per_m2_s,
        light_limited,
    ) * inputs.nutrient_feedback_fraction;
    if (!std.math.isFinite(rate) or rate < 0)
        return error.InvalidCanopyC3CarboxylationRate;
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
        return error.CanopyC3RadiationDimensionMismatch;

    inline for (.{
        inputs.negligible_leaf_area_m2,
        inputs.quantum_efficiency_umol_e_per_umol_par,
        inputs.light_saturated_electron_transport_umol_per_m2_s,
        inputs.electron_transport_curvature_product,
        inputs.electron_transport_curvature_divisor,
        inputs.carboxylation_umol_co2_per_umol_electron,
        inputs.co2_limited_carboxylation_umol_per_m2_s,
        inputs.nutrient_feedback_fraction,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteCanopyC3RadiationInput;
    if (inputs.negligible_leaf_area_m2 < 0 or
        inputs.quantum_efficiency_umol_e_per_umol_par < 0 or
        inputs.light_saturated_electron_transport_umol_per_m2_s < 0 or
        inputs.electron_transport_curvature_product < 0 or
        inputs.electron_transport_curvature_divisor <= 0 or
        inputs.carboxylation_umol_co2_per_umol_electron < 0 or
        inputs.co2_limited_carboxylation_umol_per_m2_s < 0 or
        inputs.nutrient_feedback_fraction < 0 or
        inputs.nutrient_feedback_fraction > 1)
        return error.InvalidCanopyC3RadiationInput;
    for (inputs.leaf_area_m2_by_layer) |value| try nonnegative(value);
    for (inputs.unself_shaded_leaf_area_m2_by_inclination_layer) |value|
        try nonnegative(value);
    for (inputs.direct_par_umol_per_m2_s) |value| try nonnegative(value);
    for (inputs.diffuse_par_umol_per_m2_s) |value| try nonnegative(value);
    for (inputs.direct_transmission_fraction_by_layer_boundary) |value|
        try fraction(value);
    for (inputs.diffuse_transmission_fraction_by_layer_boundary) |value|
        try fraction(value);
}

fn nonnegative(value: f64) !void {
    if (!std.math.isFinite(value)) return error.NonFiniteCanopyC3RadiationInput;
    if (value < 0) return error.InvalidCanopyC3RadiationInput;
}

fn fraction(value: f64) !void {
    try nonnegative(value);
    if (value > 1) return error.InvalidCanopyC3RadiationInput;
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
        .nutrient_feedback_fraction = 0.5,
    };
}

test "C3 radiation preserves source topology order and feedback" {
    const inputs = fixture();
    const unfed =
        20 * 3 * 0.25 + 10 * 3 * 0.1 +
        40 * 3 * 0.25 + 20 * 3 * 0.1 +
        10 * 2 * 0.5 + 5 * 2 * 0.2 +
        30 * 2 * 0.5 + 15 * 2 * 0.2;
    try std.testing.expectEqual(@as(f64, unfed * 0.5), try integrate(inputs));
}

test "C3 radiation retains strict layer surface and PAR gates" {
    var inputs = fixture();
    inputs.leaf_area_m2_by_layer = &.{ 0, 1 };
    inputs.unself_shaded_leaf_area_m2_by_inclination_layer = &.{ 2, 0 };
    try std.testing.expectEqual(@as(f64, 0), try integrate(inputs));
}

test "C3 radiation rejects malformed topology and negative radiation" {
    var inputs = fixture();
    inputs.direct_par_umol_per_m2_s = &.{1};
    try std.testing.expectError(
        error.CanopyC3RadiationDimensionMismatch,
        integrate(inputs),
    );
    inputs = fixture();
    inputs.diffuse_par_umol_per_m2_s = &.{ 5, 10, -1, 20 };
    try std.testing.expectError(error.InvalidCanopyC3RadiationInput, integrate(inputs));
}
