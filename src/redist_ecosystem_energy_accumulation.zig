const std = @import("std");

/// Mutable per-cell energy diagnostics. Legacy owners: TRN/TLE/TSH/TGH and
/// ground-only TRNS/TLES/TSHS/TGHS. All values are MJ h-1.
pub const State = struct {
    ecosystem_net_radiation_megajoules_h: f64,
    ecosystem_latent_heat_megajoules_h: f64,
    ecosystem_sensible_heat_megajoules_h: f64,
    ecosystem_storage_heat_megajoules_h: f64,
    ground_net_radiation_megajoules_h: f64,
    ground_latent_heat_megajoules_h: f64,
    ground_sensible_heat_megajoules_h: f64,
    ground_storage_heat_megajoules_h: f64,
};

/// Converged surface flux increments owned by WATSUB/REDIST, MJ per hourly step.
pub const Inputs = struct {
    absorbed_surface_radiation_megajoules: f64, // HEATI
    surface_latent_heat_megajoules: f64, // HEATE
    surface_sensible_heat_megajoules: f64, // HEATS
    surface_storage_heat_megajoules: f64, // HEATH
    advected_water_heat_megajoules: f64, // HEATV
};

fn finiteStruct(value: anytype) bool {
    inline for (std.meta.fields(@TypeOf(value))) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 10624--10631 in exact assignment order.
pub fn accumulate(state: *State, inputs: Inputs) !void {
    if (!finiteStruct(state.*) or !finiteStruct(inputs)) return error.InvalidEcosystemEnergyAccumulationInput;
    var next = state.*;
    next.ecosystem_net_radiation_megajoules_h = state.ecosystem_net_radiation_megajoules_h + inputs.absorbed_surface_radiation_megajoules;
    next.ecosystem_latent_heat_megajoules_h = state.ecosystem_latent_heat_megajoules_h + inputs.surface_latent_heat_megajoules;
    next.ecosystem_sensible_heat_megajoules_h = state.ecosystem_sensible_heat_megajoules_h + inputs.surface_sensible_heat_megajoules;
    next.ecosystem_storage_heat_megajoules_h = state.ecosystem_storage_heat_megajoules_h -
        (inputs.surface_storage_heat_megajoules - inputs.advected_water_heat_megajoules);
    next.ground_net_radiation_megajoules_h = state.ground_net_radiation_megajoules_h + inputs.absorbed_surface_radiation_megajoules;
    next.ground_latent_heat_megajoules_h = state.ground_latent_heat_megajoules_h + inputs.surface_latent_heat_megajoules;
    next.ground_sensible_heat_megajoules_h = state.ground_sensible_heat_megajoules_h + inputs.surface_sensible_heat_megajoules;
    next.ground_storage_heat_megajoules_h = state.ground_storage_heat_megajoules_h -
        (inputs.surface_storage_heat_megajoules - inputs.advected_water_heat_megajoules);
    if (!finiteStruct(next)) return error.NonFiniteEcosystemEnergyAccumulationResult;
    state.* = next;
}

test "REDIST ecosystem energy accumulation preserves exact order and storage sign" {
    var state = State{
        .ecosystem_net_radiation_megajoules_h = 1,
        .ecosystem_latent_heat_megajoules_h = 2,
        .ecosystem_sensible_heat_megajoules_h = 3,
        .ecosystem_storage_heat_megajoules_h = 4,
        .ground_net_radiation_megajoules_h = 5,
        .ground_latent_heat_megajoules_h = 6,
        .ground_sensible_heat_megajoules_h = 7,
        .ground_storage_heat_megajoules_h = 8,
    };
    try accumulate(&state, .{ .absorbed_surface_radiation_megajoules = 10, .surface_latent_heat_megajoules = 20, .surface_sensible_heat_megajoules = 30, .surface_storage_heat_megajoules = 40, .advected_water_heat_megajoules = 4 });
    try std.testing.expectEqual(@as(f64, 11), state.ecosystem_net_radiation_megajoules_h);
    try std.testing.expectEqual(@as(f64, 22), state.ecosystem_latent_heat_megajoules_h);
    try std.testing.expectEqual(@as(f64, 33), state.ecosystem_sensible_heat_megajoules_h);
    try std.testing.expectEqual(@as(f64, -32), state.ecosystem_storage_heat_megajoules_h);
    try std.testing.expectEqual(@as(f64, 15), state.ground_net_radiation_megajoules_h);
    try std.testing.expectEqual(@as(f64, -28), state.ground_storage_heat_megajoules_h);
}

test "REDIST ecosystem energy accumulation validation is atomic" {
    var state = State{
        .ecosystem_net_radiation_megajoules_h = 1,
        .ecosystem_latent_heat_megajoules_h = 2,
        .ecosystem_sensible_heat_megajoules_h = 3,
        .ecosystem_storage_heat_megajoules_h = 4,
        .ground_net_radiation_megajoules_h = 5,
        .ground_latent_heat_megajoules_h = 6,
        .ground_sensible_heat_megajoules_h = 7,
        .ground_storage_heat_megajoules_h = 8,
    };
    const before = state;
    try std.testing.expectError(error.InvalidEcosystemEnergyAccumulationInput, accumulate(&state, .{
        .absorbed_surface_radiation_megajoules = 1,
        .surface_latent_heat_megajoules = 2,
        .surface_sensible_heat_megajoules = std.math.nan(f64),
        .surface_storage_heat_megajoules = 4,
        .advected_water_heat_megajoules = 5,
    }));
    try std.testing.expectEqualDeep(before, state);
}
