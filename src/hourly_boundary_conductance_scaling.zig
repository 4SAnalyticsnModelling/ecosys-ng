const std = @import("std");

pub const State = struct {
    latent_boundary_numerator_m2_h_per_iteration: f64,
    sensible_boundary_numerator_m2_megajoules_per_k_iteration: f64,
};

pub const Inputs = struct {
    first_subhourly_iteration: bool,
    horizontal_cell_area_m2: f64,
    water_heat_iteration_fraction_h: f64,
    volumetric_air_heat_capacity_megajoules_per_m3_k: f64 = 1.25e-3,
};

/// Exact HOUR1 PAREX/PARSX refresh from hour1.f:179-183.
///
/// HOUR1 executes these assignments only for NFZ==1. Later nonlinear
/// iterations retain the accepted hourly numerators.
pub fn refresh(state: *State, inputs: Inputs) !bool {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state.*, field.name)))
            return error.NonFiniteBoundaryConductanceScalingState;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == bool) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteBoundaryConductanceScalingInput;
    }
    if (inputs.horizontal_cell_area_m2 <= 0 or
        inputs.water_heat_iteration_fraction_h <= 0 or
        inputs.water_heat_iteration_fraction_h > 1 or
        inputs.volumetric_air_heat_capacity_megajoules_per_m3_k <= 0)
        return error.InvalidBoundaryConductanceScalingInput;
    if (!inputs.first_subhourly_iteration) return false;

    const next: State = .{
        .latent_boundary_numerator_m2_h_per_iteration = inputs.horizontal_cell_area_m2 *
            inputs.water_heat_iteration_fraction_h,
        .sensible_boundary_numerator_m2_megajoules_per_k_iteration = inputs.horizontal_cell_area_m2 *
            inputs.volumetric_air_heat_capacity_megajoules_per_m3_k *
            inputs.water_heat_iteration_fraction_h,
    };
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.BoundaryConductanceScalingOverflow;
    state.* = next;
    return true;
}

test "first subhour reproduces exact PAREX and PARSX equations" {
    var state: State = .{
        .latent_boundary_numerator_m2_h_per_iteration = 0,
        .sensible_boundary_numerator_m2_megajoules_per_k_iteration = 0,
    };
    try std.testing.expect(try refresh(&state, .{
        .first_subhourly_iteration = true,
        .horizontal_cell_area_m2 = 100,
        .water_heat_iteration_fraction_h = 0.05,
    }));
    try std.testing.expectEqual(
        @as(f64, 5),
        state.latent_boundary_numerator_m2_h_per_iteration,
    );
    try std.testing.expectEqual(
        @as(f64, 0.00625),
        state.sensible_boundary_numerator_m2_megajoules_per_k_iteration,
    );
}

test "later nonlinear iteration retains accepted hourly numerators" {
    var state: State = .{
        .latent_boundary_numerator_m2_h_per_iteration = 7,
        .sensible_boundary_numerator_m2_megajoules_per_k_iteration = 8,
    };
    const before = state;
    try std.testing.expect(!try refresh(&state, .{
        .first_subhourly_iteration = false,
        .horizontal_cell_area_m2 = 200,
        .water_heat_iteration_fraction_h = 0.1,
    }));
    try std.testing.expectEqualDeep(before, state);
}

test "runtime air heat capacity controls sensible numerator only" {
    var state: State = .{
        .latent_boundary_numerator_m2_h_per_iteration = 0,
        .sensible_boundary_numerator_m2_megajoules_per_k_iteration = 0,
    };
    _ = try refresh(&state, .{
        .first_subhourly_iteration = true,
        .horizontal_cell_area_m2 = 10,
        .water_heat_iteration_fraction_h = 0.2,
        .volumetric_air_heat_capacity_megajoules_per_m3_k = 0.002,
    });
    try std.testing.expectEqual(
        @as(f64, 2),
        state.latent_boundary_numerator_m2_h_per_iteration,
    );
    try std.testing.expectEqual(
        @as(f64, 0.004),
        state.sensible_boundary_numerator_m2_megajoules_per_k_iteration,
    );
}

test "invalid late input leaves state unchanged" {
    var state: State = .{
        .latent_boundary_numerator_m2_h_per_iteration = 7,
        .sensible_boundary_numerator_m2_megajoules_per_k_iteration = 8,
    };
    const before = state;
    try std.testing.expectError(
        error.NonFiniteBoundaryConductanceScalingInput,
        refresh(&state, .{
            .first_subhourly_iteration = true,
            .horizontal_cell_area_m2 = 100,
            .water_heat_iteration_fraction_h = std.math.nan(f64),
        }),
    );
    try std.testing.expectEqualDeep(before, state);
}

test "overflow leaves state unchanged" {
    var state: State = .{
        .latent_boundary_numerator_m2_h_per_iteration = 7,
        .sensible_boundary_numerator_m2_megajoules_per_k_iteration = 8,
    };
    const before = state;
    try std.testing.expectError(
        error.BoundaryConductanceScalingOverflow,
        refresh(&state, .{
            .first_subhourly_iteration = true,
            .horizontal_cell_area_m2 = std.math.floatMax(f64),
            .water_heat_iteration_fraction_h = 1,
            .volumetric_air_heat_capacity_megajoules_per_m3_k = 2,
        }),
    );
    try std.testing.expectEqualDeep(before, state);
}
