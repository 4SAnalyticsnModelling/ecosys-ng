const std = @import("std");

pub const micropore_species_count = 50;
pub const macropore_species_count = 49;

pub const BoundaryState = struct {
    current_volumetric_water_content: *f64,
    adjacent_volumetric_water_content: *f64,
    micropore_flux_amount_per_step: []f64,
    macropore_flux_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 6974--7076.
/// This is the ELSE of the strict dual `DLYR(3,...) > DLYRM` guard, so either
/// layer at or below the threshold clears both THETW1 values and all fluxes.
pub fn resetIfThin(
    current_layer_thickness_m: f64,
    adjacent_layer_thickness_m: f64,
    minimum_transport_layer_thickness_m: f64,
    state: BoundaryState,
) !bool {
    inline for (.{ current_layer_thickness_m, adjacent_layer_thickness_m, minimum_transport_layer_thickness_m }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteThinSoilBoundarySoluteFluxInput;
    if (current_layer_thickness_m < 0 or adjacent_layer_thickness_m < 0 or minimum_transport_layer_thickness_m < 0)
        return error.InvalidThinSoilBoundarySoluteFluxInput;
    if (current_layer_thickness_m > minimum_transport_layer_thickness_m and
        adjacent_layer_thickness_m > minimum_transport_layer_thickness_m)
        return false;
    if (state.micropore_flux_amount_per_step.len != micropore_species_count or
        state.macropore_flux_amount_per_step.len != macropore_species_count)
        return error.ThinSoilBoundarySoluteFluxDimensionMismatch;
    if (!std.math.isFinite(state.current_volumetric_water_content.*) or
        !std.math.isFinite(state.adjacent_volumetric_water_content.*))
        return error.NonFiniteThinSoilBoundarySoluteFluxInput;
    for (state.micropore_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteThinSoilBoundarySoluteFluxInput;
    for (state.macropore_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteThinSoilBoundarySoluteFluxInput;

    state.current_volumetric_water_content.* = 0;
    state.adjacent_volumetric_water_content.* = 0;
    @memset(state.micropore_flux_amount_per_step, 0);
    @memset(state.macropore_flux_amount_per_step, 0);
    return true;
}

test "TRNSFRS current layer at threshold triggers complete reset" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    const applied = try resetIfThin(0.01, 0.2, 0.01, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux });
    try std.testing.expect(applied);
    try std.testing.expectEqual(@as(f64, 0), current_theta);
    try std.testing.expectEqual(@as(f64, 0), adjacent_theta);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** micropore_species_count), &micropore_flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** macropore_species_count), &macropore_flux);
}

test "TRNSFRS adjacent layer below threshold triggers complete reset" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    try std.testing.expect(try resetIfThin(0.2, 0.005, 0.01, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux }));
    try std.testing.expectEqual(@as(f64, 0), micropore_flux[49]);
    try std.testing.expectEqual(@as(f64, 0), macropore_flux[48]);
}

test "TRNSFRS dual strict pass leaves state untouched" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    const applied = try resetIfThin(0.0101, 0.02, 0.01, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux });
    try std.testing.expect(!applied);
    try std.testing.expectEqual(@as(f64, 0.3), current_theta);
    try std.testing.expectEqual(@as(f64, 2), micropore_flux[49]);
}

test "thin-path dimension error leaves all targets unchanged" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var short_micro = [_]f64{2} ** (micropore_species_count - 1);
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    try std.testing.expectError(error.ThinSoilBoundarySoluteFluxDimensionMismatch, resetIfThin(0.01, 0.2, 0.01, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &short_micro, .macropore_flux_amount_per_step = &macropore_flux }));
    try std.testing.expectEqual(@as(f64, 0.3), current_theta);
    try std.testing.expectEqual(@as(f64, 0.4), adjacent_theta);
}

test "late invalid macropore value keeps thin reset atomic" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    macropore_flux[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteThinSoilBoundarySoluteFluxInput, resetIfThin(0.01, 0.2, 0.01, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux }));
    try std.testing.expectEqual(@as(f64, 0.3), current_theta);
    try std.testing.expectEqual(@as(f64, 2), micropore_flux[0]);
}
