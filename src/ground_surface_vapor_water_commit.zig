const std = @import("std");

/// Commits the realized implicit ground-air vapor exchange to prognostic
/// liquid water. This is the storage counterpart of the surface vapor source
/// in the ground-air solve; positive exchange is evaporation to ground air.
pub const Context = struct {
    surface_vapor_conductance_m3_per_h: []const f64,
    surface_vapor_fraction: []const f64,
    accepted_ground_air_vapor_fraction: []const f64,
    litter_liquid_water_m3: []f64,
    soil_matrix_liquid_water_m3: []f64,
    active_soil_layer_count: []const usize,
    soil_layer_capacity: usize,
    evaporation_m3_per_h: []f64,
    condensation_m3_per_h: []f64,
    litter_liquid_water_change_m3: []f64,
    topsoil_liquid_water_change_m3: []f64,
};

/// Legacy traceability: WATSUB.F accumulates WFLVR/WFLVL and commits them to
/// VOLW1 at lines 6710 and 6811. ecosys-ng uses the accepted implicit
/// ground-air exchange and removes pure water from litter before topsoil.
pub fn commit(context: Context) !void {
    const cells = context.litter_liquid_water_m3.len;
    if (cells == 0 or context.soil_layer_capacity == 0 or
        context.surface_vapor_conductance_m3_per_h.len != cells or
        context.surface_vapor_fraction.len != cells or
        context.accepted_ground_air_vapor_fraction.len != cells or
        context.active_soil_layer_count.len != cells or
        context.evaporation_m3_per_h.len != cells or
        context.condensation_m3_per_h.len != cells or
        context.litter_liquid_water_change_m3.len != cells or
        context.topsoil_liquid_water_change_m3.len != cells or
        context.soil_matrix_liquid_water_m3.len != cells * context.soil_layer_capacity)
        return error.GroundSurfaceVaporWaterDimensionMismatch;

    // Validate the complete transaction before mutating any storage.
    for (0..cells) |cell| {
        if (context.active_soil_layer_count[cell] == 0 or
            context.active_soil_layer_count[cell] > context.soil_layer_capacity)
            return error.InvalidGroundSurfaceSoilLayerCount;
        const topsoil = cell * context.soil_layer_capacity;
        inline for (.{
            context.surface_vapor_conductance_m3_per_h[cell],
            context.surface_vapor_fraction[cell],
            context.accepted_ground_air_vapor_fraction[cell],
            context.litter_liquid_water_m3[cell],
            context.soil_matrix_liquid_water_m3[topsoil],
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidGroundSurfaceVaporWaterInput;
        const signed_evaporation = context.surface_vapor_conductance_m3_per_h[cell] *
            (context.surface_vapor_fraction[cell] - context.accepted_ground_air_vapor_fraction[cell]);
        if (!std.math.isFinite(signed_evaporation))
            return error.NonFiniteGroundSurfaceVaporExchange;
        if (signed_evaporation > context.litter_liquid_water_m3[cell] +
            context.soil_matrix_liquid_water_m3[topsoil])
            return error.InsufficientGroundSurfaceLiquidWater;
    }

    for (0..cells) |cell| {
        const topsoil = cell * context.soil_layer_capacity;
        const signed_evaporation = context.surface_vapor_conductance_m3_per_h[cell] *
            (context.surface_vapor_fraction[cell] - context.accepted_ground_air_vapor_fraction[cell]);
        const evaporation = @max(0, signed_evaporation);
        const condensation = @max(0, -signed_evaporation);
        const litter_evaporation = @min(evaporation, context.litter_liquid_water_m3[cell]);
        const topsoil_evaporation = evaporation - litter_evaporation;
        context.litter_liquid_water_change_m3[cell] = condensation - litter_evaporation;
        context.topsoil_liquid_water_change_m3[cell] = -topsoil_evaporation;
        context.litter_liquid_water_m3[cell] += context.litter_liquid_water_change_m3[cell];
        context.soil_matrix_liquid_water_m3[topsoil] += context.topsoil_liquid_water_change_m3[cell];
        context.evaporation_m3_per_h[cell] = evaporation;
        context.condensation_m3_per_h[cell] = condensation;
    }
}

test "accepted evaporation conservatively withdraws litter then topsoil" {
    var litter = [_]f64{2};
    var soil = [_]f64{ 5, 7 };
    var evaporation = [_]f64{0};
    var condensation = [_]f64{0};
    var litter_change = [_]f64{0};
    var topsoil_change = [_]f64{0};
    try commit(.{
        .surface_vapor_conductance_m3_per_h = &.{3},
        .surface_vapor_fraction = &.{2},
        .accepted_ground_air_vapor_fraction = &.{0},
        .litter_liquid_water_m3 = &litter,
        .soil_matrix_liquid_water_m3 = &soil,
        .active_soil_layer_count = &.{2},
        .soil_layer_capacity = 2,
        .evaporation_m3_per_h = &evaporation,
        .condensation_m3_per_h = &condensation,
        .litter_liquid_water_change_m3 = &litter_change,
        .topsoil_liquid_water_change_m3 = &topsoil_change,
    });
    try std.testing.expectEqual(@as(f64, 6), evaporation[0]);
    try std.testing.expectEqual(@as(f64, 0), litter[0]);
    try std.testing.expectEqual(@as(f64, 1), soil[0]);
    try std.testing.expectEqual(@as(f64, -2), litter_change[0]);
    try std.testing.expectEqual(@as(f64, -4), topsoil_change[0]);
}

test "accepted condensation enters litter storage" {
    var litter = [_]f64{2};
    var soil = [_]f64{5};
    var evaporation = [_]f64{0};
    var condensation = [_]f64{0};
    var litter_change = [_]f64{0};
    var topsoil_change = [_]f64{0};
    try commit(.{
        .surface_vapor_conductance_m3_per_h = &.{4},
        .surface_vapor_fraction = &.{0.25},
        .accepted_ground_air_vapor_fraction = &.{0.5},
        .litter_liquid_water_m3 = &litter,
        .soil_matrix_liquid_water_m3 = &soil,
        .active_soil_layer_count = &.{1},
        .soil_layer_capacity = 1,
        .evaporation_m3_per_h = &evaporation,
        .condensation_m3_per_h = &condensation,
        .litter_liquid_water_change_m3 = &litter_change,
        .topsoil_liquid_water_change_m3 = &topsoil_change,
    });
    try std.testing.expectEqual(@as(f64, 1), condensation[0]);
    try std.testing.expectEqual(@as(f64, 3), litter[0]);
    try std.testing.expectEqual(@as(f64, 5), soil[0]);
}

test "insufficient storage rejects the full transaction atomically" {
    var litter = [_]f64{ 1, 1 };
    var soil = [_]f64{ 2, 2 };
    var output = [_]f64{ 9, 9 };
    try std.testing.expectError(error.InsufficientGroundSurfaceLiquidWater, commit(.{
        .surface_vapor_conductance_m3_per_h = &.{ 1, 4 },
        .surface_vapor_fraction = &.{ 1, 1 },
        .accepted_ground_air_vapor_fraction = &.{ 0, 0 },
        .litter_liquid_water_m3 = &litter,
        .soil_matrix_liquid_water_m3 = &soil,
        .active_soil_layer_count = &.{ 1, 1 },
        .soil_layer_capacity = 1,
        .evaporation_m3_per_h = &output,
        .condensation_m3_per_h = &output,
        .litter_liquid_water_change_m3 = &output,
        .topsoil_liquid_water_change_m3 = &output,
    }));
    try std.testing.expectEqualSlices(f64, &.{ 1, 1 }, &litter);
    try std.testing.expectEqualSlices(f64, &.{ 2, 2 }, &soil);
}
