const std = @import("std");

pub const micropore_species_count = 50;
pub const macropore_species_count = 49;
pub const Direction = enum { east, south, vertical };
pub const ProfileEligibility = enum { active_pair, inactive_pair };

pub const BoundaryState = struct {
    current_volumetric_water_content: *f64,
    adjacent_volumetric_water_content: *f64,
    micropore_flux_amount_per_step: []f64,
    macropore_flux_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 6871--6973.
/// This ELSEIF applies only when the preceding active-profile test failed and
/// the direction N is not vertical. It clears THETW1 and both flux families.
pub fn resetIfRequired(eligibility: ProfileEligibility, direction: Direction, state: BoundaryState) !bool {
    if (eligibility != .inactive_pair or direction == .vertical) return false;
    if (state.micropore_flux_amount_per_step.len != micropore_species_count or
        state.macropore_flux_amount_per_step.len != macropore_species_count)
        return error.InactiveSoilBoundarySoluteFluxDimensionMismatch;
    if (!std.math.isFinite(state.current_volumetric_water_content.*) or
        !std.math.isFinite(state.adjacent_volumetric_water_content.*))
        return error.NonFiniteInactiveSoilBoundarySoluteFluxInput;
    for (state.micropore_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteInactiveSoilBoundarySoluteFluxInput;
    for (state.macropore_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteInactiveSoilBoundarySoluteFluxInput;

    state.current_volumetric_water_content.* = 0;
    state.adjacent_volumetric_water_content.* = 0;
    @memset(state.micropore_flux_amount_per_step, 0);
    @memset(state.macropore_flux_amount_per_step, 0);
    return true;
}

test "TRNSFRS inactive horizontal pair clears THETW1 and both flux topologies" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    const applied = try resetIfRequired(.inactive_pair, .east, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux });
    try std.testing.expect(applied);
    try std.testing.expectEqual(@as(f64, 0), current_theta);
    try std.testing.expectEqual(@as(f64, 0), adjacent_theta);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** micropore_species_count), &micropore_flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** macropore_species_count), &macropore_flux);
}

test "TRNSFRS vertical direction bypasses inactive-pair reset" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    const applied = try resetIfRequired(.inactive_pair, .vertical, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux });
    try std.testing.expect(!applied);
    try std.testing.expectEqual(@as(f64, 0.3), current_theta);
    try std.testing.expectEqual(@as(f64, 2), micropore_flux[49]);
    try std.testing.expectEqual(@as(f64, 3), macropore_flux[48]);
}

test "TRNSFRS active horizontal pair bypasses inactive-pair reset" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    const applied = try resetIfRequired(.active_pair, .south, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux });
    try std.testing.expect(!applied);
    try std.testing.expectEqual(@as(f64, 0.4), adjacent_theta);
}

test "dimension error leaves all active reset targets unchanged" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var short_micro = [_]f64{2} ** (micropore_species_count - 1);
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    try std.testing.expectError(error.InactiveSoilBoundarySoluteFluxDimensionMismatch, resetIfRequired(.inactive_pair, .east, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &short_micro, .macropore_flux_amount_per_step = &macropore_flux }));
    try std.testing.expectEqual(@as(f64, 0.3), current_theta);
    try std.testing.expectEqual(@as(f64, 0.4), adjacent_theta);
}

test "late invalid macropore value leaves micropore state atomic" {
    var current_theta: f64 = 0.3;
    var adjacent_theta: f64 = 0.4;
    var micropore_flux = [_]f64{2} ** micropore_species_count;
    var macropore_flux = [_]f64{3} ** macropore_species_count;
    macropore_flux[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteInactiveSoilBoundarySoluteFluxInput, resetIfRequired(.inactive_pair, .south, .{ .current_volumetric_water_content = &current_theta, .adjacent_volumetric_water_content = &adjacent_theta, .micropore_flux_amount_per_step = &micropore_flux, .macropore_flux_amount_per_step = &macropore_flux }));
    try std.testing.expectEqual(@as(f64, 0.3), current_theta);
    try std.testing.expectEqual(@as(f64, 2), micropore_flux[0]);
}
