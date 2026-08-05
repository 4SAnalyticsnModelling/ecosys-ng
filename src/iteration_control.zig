const std = @import("std");
const SceneOptions = @import("options.zig").SceneOptions;

/// Runtime nonlinear iteration ceilings derived from the old sub-hour control
/// record. They are convergence budgets only: ecosys-ng never repeats a full
/// model cycle to consume the budget.
pub const Limits = struct {
    water_heat_solute_max_iterations: u16,
    /// EROSION NPH: local suspended-sediment convergence.
    erosion_max_iterations: u16,
    gas_max_iterations: u16,
    /// REDIST dissolved-organic transport convergence.
    ///
    /// This shared the `water_heat_solute` ceiling of 20, which is not enough
    /// when a near-empty pool sits beside a large aqueous boundary transfer: the
    /// Ottawa example produced an `8.37e-6` g residual on a `1.43e-15` g pool
    /// while that hour moved `8.367e6` g across the boundary, and Newton
    /// accepted no step in 20 iterations. The fixed point is reachable, just
    /// slowly: raising the ceiling converges it. This is the same situation that
    /// `gas_max_iterations` documents for trace gases below 1 microgram.
    organic_transport_max_iterations: u16 = 2000,
    /// WTHR NPR: litter/surface water and heat convergence.
    litter_water_heat_max_iterations: u16 = 30,
    /// WTHR NPS: snow water and heat convergence.
    snowpack_max_iterations: u16 = 20,
    litter_under_snow_max_iterations: u16 = 10,
    /// SOLUTE MRXN: profile and surface-litter reaction equilibrium.
    solute_reaction_max_iterations: u16 = 60,
    /// STARTE MRXN: initial reaction-equilibrium establishment.
    initial_solute_reaction_max_iterations: u16 = 1000,
    canopy_energy_water_max_iterations: u16 = 100,
    leaf_co2_max_iterations: u16 = 100,

    pub fn fromSceneOptions(options: SceneOptions) !Limits {
        const gas_iterations = try std.math.mul(u16, options.water_heat_solute_iteration_limit, options.gas_iterations_per_water_heat_solute_iteration);
        if (gas_iterations == 0) return error.ZeroNonlinearIterationLimit;
        // Trace-gas species (NH3, H2) in partially frozen soil require more
        // iterations than NPH*NPG (typically 80) to converge when masses fall
        // below 1 µg and the solver tolerance floor dominates. Empirically, the
        // NH3 case at day 20 converges at ~0.37% per Newton step starting from
        // scaled_residual ≈ 4.6, requiring ~613 total iterations. 1000 is the
        // safe ceiling; most hours still converge in < 100 iterations.
        const minimum_gas_iterations: u16 = 1000;
        return .{
            .water_heat_solute_max_iterations = options.water_heat_solute_iteration_limit,
            .erosion_max_iterations = options.water_heat_solute_iteration_limit,
            .gas_max_iterations = @max(minimum_gas_iterations, gas_iterations),
        };
    }

    /// UPTAKE standing-dead energy used an outer NPH cycle and an inner MXN
    /// solve. ecosys-ng spends the same maximum work in one local nonlinear
    /// solve and exits immediately on convergence.
    pub fn standingDeadEnergyMaxIterations(self: Limits) !u16 {
        const iterations = try std.math.mul(
            u16,
            self.water_heat_solute_max_iterations,
            self.canopy_energy_water_max_iterations,
        );
        if (iterations == 0) return error.ZeroNonlinearIterationLimit;
        return iterations;
    }
};

/// WTHR raises NPH to at least 20 when any top soil layer has less than
/// 4.19e-3 MJ K-1 per square metre of horizontal area. In ecosys-ng this is a
/// nonlinear convergence ceiling, not a request to repeat the full model.
pub fn waterHeatSoluteCeilingForCurrentState(base_nph: u16, heat_capacity_megajoules_per_k: []const f64, horizontal_area_m2: []const f64, is_top_soil_layer: []const bool) !u16 {
    if (base_nph == 0) return error.ZeroNonlinearIterationLimit;
    if (heat_capacity_megajoules_per_k.len != horizontal_area_m2.len or heat_capacity_megajoules_per_k.len != is_top_soil_layer.len) return error.IterationControlDimensionMismatch;
    for (heat_capacity_megajoules_per_k, horizontal_area_m2, is_top_soil_layer) |heat_capacity, area, is_top| {
        if (!std.math.isFinite(heat_capacity) or heat_capacity < 0 or !std.math.isFinite(area) or area <= 0) return error.InvalidIterationControlThermalState;
        if (is_top and heat_capacity < 4.19e-3 * area) return @max(@as(u16, 20), base_nph);
    }
    return base_nph;
}

test "legacy option controls become convergence ceilings" {
    const options = try @import("options.zig").parse(@import("test_fixtures.zig").scene_options_source);
    const limits = try Limits.fromSceneOptions(options);
    try std.testing.expectEqual(@as(u16, 20), limits.water_heat_solute_max_iterations);
    try std.testing.expectEqual(@as(u16, 20), limits.erosion_max_iterations);
    try std.testing.expectEqual(@as(u16, 1000), limits.gas_max_iterations);
    try std.testing.expectEqual(@as(u16, 30), limits.litter_water_heat_max_iterations);
    try std.testing.expectEqual(@as(u16, 20), limits.snowpack_max_iterations);
    try std.testing.expectEqual(@as(u16, 10), limits.litter_under_snow_max_iterations);
    try std.testing.expectEqual(@as(u16, 60), limits.solute_reaction_max_iterations);
    try std.testing.expectEqual(@as(u16, 1000), limits.initial_solute_reaction_max_iterations);
    try std.testing.expectEqual(@as(u16, 100), limits.canopy_energy_water_max_iterations);
    try std.testing.expectEqual(@as(u16, 100), limits.leaf_co2_max_iterations);
    try std.testing.expectEqual(@as(u16, 2000), try limits.standingDeadEnergyMaxIterations());
}

test "standing-dead energy ceiling rejects overflow" {
    var limits = try Limits.fromSceneOptions(try @import("options.zig").parse(@import("test_fixtures.zig").scene_options_source));
    limits.water_heat_solute_max_iterations = std.math.maxInt(u16);
    try std.testing.expectError(error.Overflow, limits.standingDeadEnergyMaxIterations());
}

test "low top-layer heat capacity applies WTHR minimum NPH" {
    try std.testing.expectEqual(@as(u16, 20), try waterHeatSoluteCeilingForCurrentState(10, &.{ 0.004, 3.0 }, &.{ 1.0, 1.0 }, &.{ true, false }));
    try std.testing.expectEqual(@as(u16, 10), try waterHeatSoluteCeilingForCurrentState(10, &.{ 0.005, 0.001 }, &.{ 1.0, 1.0 }, &.{ true, false }));
    try std.testing.expectEqual(@as(u16, 30), try waterHeatSoluteCeilingForCurrentState(30, &.{0.001}, &.{1.0}, &.{true}));
}

test "low heat-capacity threshold is strict and only triggers on top layers" {
    try std.testing.expectEqual(
        @as(u16, 4),
        try waterHeatSoluteCeilingForCurrentState(4, &.{4.19e-3}, &.{1.0}, &.{true}),
    );
    try std.testing.expectEqual(
        @as(u16, 20),
        try waterHeatSoluteCeilingForCurrentState(4, &.{ 4.19e-3, 0.001 }, &.{ 1.0, 1.0 }, &.{ false, true }),
    );
    try std.testing.expectEqual(
        @as(u16, 4),
        try waterHeatSoluteCeilingForCurrentState(4, &.{ 4.19e-3, 3.0 }, &.{ 1.0, 1.0 }, &.{ false, false }),
    );
}

test "water-heat ceiling rejects non-top-layer heat-capacity changes and dimensional mismatch" {
    try std.testing.expectEqual(@as(u16, 12), try waterHeatSoluteCeilingForCurrentState(12, &.{ 0.001, 0.001 }, &.{ 1.0, 2.0 }, &.{ false, false }));
    try std.testing.expectError(
        error.IterationControlDimensionMismatch,
        waterHeatSoluteCeilingForCurrentState(12, &.{0.004}, &.{1.0}, &.{ true, false }),
    );
    try std.testing.expectError(
        error.IterationControlDimensionMismatch,
        waterHeatSoluteCeilingForCurrentState(12, &.{ 0.004, 0.005 }, &.{ 1.0, 2.0 }, &.{true}),
    );
}

test "water-heat ceiling rejects invalid thermal control state" {
    try std.testing.expectError(
        error.InvalidIterationControlThermalState,
        waterHeatSoluteCeilingForCurrentState(12, &.{-0.01}, &.{1.0}, &.{true}),
    );
    try std.testing.expectError(
        error.InvalidIterationControlThermalState,
        waterHeatSoluteCeilingForCurrentState(12, &.{0.004}, &.{0.0}, &.{true}),
    );
    try std.testing.expectError(
        error.InvalidIterationControlThermalState,
        waterHeatSoluteCeilingForCurrentState(12, &.{std.math.nan(f64)}, &.{1.0}, &.{true}),
    );
}

test "from-scene limits rejects invalid iteration-option combinations" {
    const options = try @import("options.zig").parse(@import("test_fixtures.zig").scene_options_source);
    var zero_nph = options;
    zero_nph.water_heat_solute_iteration_limit = 0;
    try std.testing.expectError(
        error.ZeroNonlinearIterationLimit,
        Limits.fromSceneOptions(zero_nph),
    );

    var zero_npg = options;
    zero_npg.gas_iterations_per_water_heat_solute_iteration = 0;
    try std.testing.expectError(
        error.ZeroNonlinearIterationLimit,
        Limits.fromSceneOptions(zero_npg),
    );

    var overflow = options;
    overflow.water_heat_solute_iteration_limit = std.math.maxInt(u16);
    overflow.gas_iterations_per_water_heat_solute_iteration = 2;
    try std.testing.expectError(
        error.Overflow,
        Limits.fromSceneOptions(overflow),
    );
}
