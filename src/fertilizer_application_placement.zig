const std = @import("std");

pub const SoilForcingAmounts = struct {
    banded_ammonium_g_n_per_m2: f64,
    banded_ammonia_g_n_per_m2: f64,
    banded_urea_g_n_per_m2: f64,
    banded_nitrate_g_n_per_m2: f64,
    banded_monocalcium_phosphate_g_p_per_m2: f64,
    calcium_carbonate_g_ca_per_m2: f64,
    calcium_sulfate_or_crushed_rock_g_per_m2: f64,
};

pub const Destination = union(enum) {
    surface_litter,
    soil_layer: usize,
};

pub const Placement = struct {
    destination: Destination,
    surface_litter_cover_fraction: f64,
    bare_soil_fraction: f64,
};

pub const Inputs = struct {
    application_depth_m: f64,
    soil_layer_bottom_depth_m: []const f64,
    surface_organic_carbon_g_c: f64,
    surface_horizontal_area_m2: f64,
    soil_forcing_amounts: SoilForcingAmounts,
    litter_cover_extinction_m2_per_g_c: f64 = 0.8e-2,
};

/// Exact HOUR1 LFDPTH/CVRDF/BAREF decision from hour1.f:273-288.
pub fn determine(inputs: Inputs) !Placement {
    inline for (.{
        inputs.application_depth_m,
        inputs.surface_organic_carbon_g_c,
        inputs.surface_horizontal_area_m2,
        inputs.litter_cover_extinction_m2_per_g_c,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteFertilizerPlacementInput;
    inline for (@typeInfo(SoilForcingAmounts).@"struct".fields) |field| {
        const value = @field(inputs.soil_forcing_amounts, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteFertilizerPlacementInput;
        if (value < 0) return error.InvalidFertilizerPlacementInput;
    }
    if (inputs.application_depth_m < 0 or
        inputs.surface_organic_carbon_g_c < 0 or
        inputs.surface_horizontal_area_m2 <= 0 or
        inputs.litter_cover_extinction_m2_per_g_c < 0)
        return error.InvalidFertilizerPlacementInput;
    var previous_bottom: f64 = 0;
    for (inputs.soil_layer_bottom_depth_m) |bottom| {
        if (!std.math.isFinite(bottom) or bottom <= previous_bottom)
            return error.InvalidFertilizerLayerBoundary;
        previous_bottom = bottom;
    }

    var soil_forcing_total: f64 = 0;
    inline for (@typeInfo(SoilForcingAmounts).@"struct".fields) |field|
        soil_forcing_total += @field(inputs.soil_forcing_amounts, field.name);
    if (!std.math.isFinite(soil_forcing_total))
        return error.FertilizerPlacementOverflow;

    if (inputs.application_depth_m <= 0 and soil_forcing_total == 0) {
        const cover = 1 - @exp(
            -inputs.litter_cover_extinction_m2_per_g_c *
                inputs.surface_organic_carbon_g_c /
                inputs.surface_horizontal_area_m2,
        );
        if (!std.math.isFinite(cover))
            return error.FertilizerPlacementOverflow;
        return .{
            .destination = .surface_litter,
            .surface_litter_cover_fraction = cover,
            .bare_soil_fraction = 1 - cover,
        };
    }

    for (inputs.soil_layer_bottom_depth_m, 0..) |bottom, layer|
        if (bottom >= inputs.application_depth_m)
            return .{
                .destination = .{ .soil_layer = layer },
                .surface_litter_cover_fraction = 1,
                .bare_soil_fraction = 0,
            };
    return error.FertilizerApplicationBelowSoilProfile;
}

fn noSoilForcing() SoilForcingAmounts {
    return std.mem.zeroes(SoilForcingAmounts);
}

test "zero-depth broadcast-only application partitions by litter cover" {
    const result = try determine(.{
        .application_depth_m = 0,
        .soil_layer_bottom_depth_m = &.{ 0.1, 0.3 },
        .surface_organic_carbon_g_c = 100,
        .surface_horizontal_area_m2 = 10,
        .soil_forcing_amounts = noSoilForcing(),
    });
    try std.testing.expect(result.destination == .surface_litter);
    const expected = 1 - @exp(-0.008 * 10.0);
    try std.testing.expectApproxEqAbs(
        expected,
        result.surface_litter_cover_fraction,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1) - expected,
        result.bare_soil_fraction,
        1e-15,
    );
}

test "banded nutrient forces zero-depth application into top soil" {
    var forcing = noSoilForcing();
    forcing.banded_ammonium_g_n_per_m2 = 1;
    const result = try determine(.{
        .application_depth_m = 0,
        .soil_layer_bottom_depth_m = &.{ 0.1, 0.3 },
        .surface_organic_carbon_g_c = 100,
        .surface_horizontal_area_m2 = 10,
        .soil_forcing_amounts = forcing,
    });
    try std.testing.expectEqual(@as(usize, 0), result.destination.soil_layer);
    try std.testing.expectEqual(
        @as(f64, 1),
        result.surface_litter_cover_fraction,
    );
}

test "application depth selects first inclusive layer bottom" {
    const at_boundary = try determine(.{
        .application_depth_m = 0.1,
        .soil_layer_bottom_depth_m = &.{ 0.1, 0.3, 0.8 },
        .surface_organic_carbon_g_c = 0,
        .surface_horizontal_area_m2 = 10,
        .soil_forcing_amounts = noSoilForcing(),
    });
    try std.testing.expectEqual(
        @as(usize, 0),
        at_boundary.destination.soil_layer,
    );
    const second = try determine(.{
        .application_depth_m = 0.1001,
        .soil_layer_bottom_depth_m = &.{ 0.1, 0.3, 0.8 },
        .surface_organic_carbon_g_c = 0,
        .surface_horizontal_area_m2 = 10,
        .soil_forcing_amounts = noSoilForcing(),
    });
    try std.testing.expectEqual(@as(usize, 1), second.destination.soil_layer);
}

test "runtime extinction coefficient controls only surface partition" {
    const result = try determine(.{
        .application_depth_m = 0,
        .soil_layer_bottom_depth_m = &.{0.1},
        .surface_organic_carbon_g_c = 100,
        .surface_horizontal_area_m2 = 10,
        .soil_forcing_amounts = noSoilForcing(),
        .litter_cover_extinction_m2_per_g_c = 0,
    });
    try std.testing.expectEqual(
        @as(f64, 0),
        result.surface_litter_cover_fraction,
    );
    try std.testing.expectEqual(@as(f64, 1), result.bare_soil_fraction);
}

test "invalid late boundary and below-profile depth fail immediately" {
    try std.testing.expectError(
        error.InvalidFertilizerLayerBoundary,
        determine(.{
            .application_depth_m = 0.2,
            .soil_layer_bottom_depth_m = &.{ 0.1, 0.1 },
            .surface_organic_carbon_g_c = 0,
            .surface_horizontal_area_m2 = 10,
            .soil_forcing_amounts = noSoilForcing(),
        }),
    );
    try std.testing.expectError(
        error.FertilizerApplicationBelowSoilProfile,
        determine(.{
            .application_depth_m = 1,
            .soil_layer_bottom_depth_m = &.{ 0.1, 0.3 },
            .surface_organic_carbon_g_c = 0,
            .surface_horizontal_area_m2 = 10,
            .soil_forcing_amounts = noSoilForcing(),
        }),
    );
}
