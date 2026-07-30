const std = @import("std");

pub const Inputs = struct {
    top_mineral_layer_sand_mass_Mg: f64,
    top_mineral_layer_silt_mass_Mg: f64,
    top_mineral_layer_clay_mass_Mg: f64,
    surface_noncharcoal_organic_carbon_g_c: f64,
    initial_litter_water_m3_per_g_c: f64,
};

pub const State = struct {
    top_mineral_layer_mass_Mg: *f64,
    surface_litter_water_m3: *f64,
};

/// Exact source-order translation of legacy `STARTS` lines 1646--1648.
///
/// The surface-water input is explicitly noncharcoal carbon because legacy
/// `ORGC` excludes the separately accumulated charcoal pool `ORGCC`.
pub fn initialize(state: State, inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceTopsoilInterfaceInput;
        if (value < 0.0) return error.InvalidSurfaceTopsoilInterfaceInput;
    }
    const top_mineral_layer_mass_Mg = @max(
        0.0,
        inputs.top_mineral_layer_sand_mass_Mg +
            inputs.top_mineral_layer_silt_mass_Mg +
            inputs.top_mineral_layer_clay_mass_Mg,
    );
    const surface_litter_water_m3 =
        inputs.initial_litter_water_m3_per_g_c *
        inputs.surface_noncharcoal_organic_carbon_g_c;
    if (!std.math.isFinite(top_mineral_layer_mass_Mg) or
        !std.math.isFinite(surface_litter_water_m3))
        return error.NonFiniteSurfaceTopsoilInterfaceResult;

    state.top_mineral_layer_mass_Mg.* = top_mineral_layer_mass_Mg;
    state.surface_litter_water_m3.* = surface_litter_water_m3;
}

test "STARTS final block refreshes mineral mass and noncharcoal litter water" {
    var mineral_mass_Mg: f64 = 0.0;
    var litter_water_m3: f64 = 0.0;
    try initialize(.{
        .top_mineral_layer_mass_Mg = &mineral_mass_Mg,
        .surface_litter_water_m3 = &litter_water_m3,
    }, .{
        .top_mineral_layer_sand_mass_Mg = 2.0,
        .top_mineral_layer_silt_mass_Mg = 3.0,
        .top_mineral_layer_clay_mass_Mg = 4.0,
        .surface_noncharcoal_organic_carbon_g_c = 100.0,
        .initial_litter_water_m3_per_g_c = 8.0e-6,
    });
    try std.testing.expectEqual(@as(f64, 9.0), mineral_mass_Mg);
    try std.testing.expectApproxEqAbs(
        @as(f64, 8.0e-4),
        litter_water_m3,
        1.0e-18,
    );
}

test "overflow fails before either state field mutates" {
    var mineral_mass_Mg: f64 = 7.0;
    var litter_water_m3: f64 = 8.0;
    try std.testing.expectError(
        error.NonFiniteSurfaceTopsoilInterfaceResult,
        initialize(.{
            .top_mineral_layer_mass_Mg = &mineral_mass_Mg,
            .surface_litter_water_m3 = &litter_water_m3,
        }, .{
            .top_mineral_layer_sand_mass_Mg = std.math.floatMax(f64),
            .top_mineral_layer_silt_mass_Mg = std.math.floatMax(f64),
            .top_mineral_layer_clay_mass_Mg = 0.0,
            .surface_noncharcoal_organic_carbon_g_c = 1.0,
            .initial_litter_water_m3_per_g_c = 1.0,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7.0), mineral_mass_Mg);
    try std.testing.expectEqual(@as(f64, 8.0), litter_water_m3);
}
