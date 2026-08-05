const std = @import("std");

pub const Inputs = struct {
    layer_index: usize,
    top_mineral_layer_index: usize,
    layer_thickness_m: f64,
    minimum_active_layer_thickness_m: f64,
    matric_potential_mpa: f64,
    osmotic_potential_mpa: f64,
};

pub const Admission = struct {
    combined_water_potential_mpa: f64,
    /// Legacy KL=2 converted from one-based to zero-based indexing.
    surface_substrate_index: ?usize,
};

pub const Result = struct {
    /// Legacy RDOSL is reset before either layer-admission predicate.
    lignin_acidity_accumulator: f64,
    admission: ?Admission,
};

/// Exact source-order translation of NITRO.F lines 253--259.
///
/// Water potentials are MPa and layer thicknesses are m. The disabled VOLX
/// predicate at source line 255 is intentionally not introduced as a gate.
pub fn evaluate(inputs: Inputs) !Result {
    inline for (.{
        inputs.layer_thickness_m,
        inputs.minimum_active_layer_thickness_m,
        inputs.matric_potential_mpa,
        inputs.osmotic_potential_mpa,
    }) |value| if (!std.math.isFinite(value)) return error.NonFiniteMicrobialLayerAdmissionInput;
    if (inputs.layer_thickness_m < 0 or inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidMicrobialLayerThickness;
    if (inputs.top_mineral_layer_index == 0)
        return error.InvalidTopMineralLayerIndex;

    const combined_water_potential_mpa = inputs.matric_potential_mpa + inputs.osmotic_potential_mpa;
    if (!std.math.isFinite(combined_water_potential_mpa))
        return error.NonFiniteCombinedWaterPotential;

    var result: Result = .{
        .lignin_acidity_accumulator = 0,
        .admission = null,
    };
    if (inputs.layer_thickness_m > inputs.minimum_active_layer_thickness_m) {
        if (inputs.layer_index == 0 or inputs.layer_index >= inputs.top_mineral_layer_index) {
            result.admission = .{
                .combined_water_potential_mpa = combined_water_potential_mpa,
                .surface_substrate_index = if (inputs.layer_index == 0) 1 else null,
            };
        }
    }
    return result;
}

test "NITRO surface admission preserves reset potential sum and KL selection" {
    const result = try evaluate(.{
        .layer_index = 0,
        .top_mineral_layer_index = 2,
        .layer_thickness_m = 0.02,
        .minimum_active_layer_thickness_m = 0.001,
        .matric_potential_mpa = -0.4,
        .osmotic_potential_mpa = -0.1,
    });
    try std.testing.expectEqual(@as(f64, 0), result.lignin_acidity_accumulator);
    try std.testing.expectEqual(@as(f64, -0.5), result.admission.?.combined_water_potential_mpa);
    try std.testing.expectEqual(@as(?usize, 1), result.admission.?.surface_substrate_index);
}

test "NITRO mineral admission begins at the runtime top mineral layer" {
    const skipped = try evaluate(.{
        .layer_index = 1,
        .top_mineral_layer_index = 2,
        .layer_thickness_m = 0.1,
        .minimum_active_layer_thickness_m = 0.01,
        .matric_potential_mpa = -0.2,
        .osmotic_potential_mpa = 0,
    });
    try std.testing.expect(skipped.admission == null);

    const admitted = try evaluate(.{
        .layer_index = 2,
        .top_mineral_layer_index = 2,
        .layer_thickness_m = 0.1,
        .minimum_active_layer_thickness_m = 0.01,
        .matric_potential_mpa = -0.2,
        .osmotic_potential_mpa = 0,
    });
    try std.testing.expect(admitted.admission != null);
    try std.testing.expect(admitted.admission.?.surface_substrate_index == null);
}

test "NITRO layer thickness gate remains strict" {
    const result = try evaluate(.{
        .layer_index = 0,
        .top_mineral_layer_index = 1,
        .layer_thickness_m = 0.01,
        .minimum_active_layer_thickness_m = 0.01,
        .matric_potential_mpa = 0,
        .osmotic_potential_mpa = 0,
    });
    try std.testing.expect(result.admission == null);
}

test "microbial layer admission rejects invalid runtime domains" {
    try std.testing.expectError(error.InvalidTopMineralLayerIndex, evaluate(.{
        .layer_index = 0,
        .top_mineral_layer_index = 0,
        .layer_thickness_m = 0.01,
        .minimum_active_layer_thickness_m = 0,
        .matric_potential_mpa = 0,
        .osmotic_potential_mpa = 0,
    }));
}
