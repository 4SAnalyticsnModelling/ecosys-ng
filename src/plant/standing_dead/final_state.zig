const std = @import("std");

pub const Inputs = struct {
    solved_canopy_air_temperature_k: f64,
    solved_canopy_air_vapor_volume_fraction: f64,
    solved_surface_temperature_k: f64,
};

pub const State = struct {
    canopy_air_temperature_k: f64,
    canopy_air_vapor_volume_fraction: f64,
    surface_temperature_k: f64,
    surface_temperature_c: f64,
};

/// UPTAKE.F 4276--4279 final standing-dead temperature and vapor commit.
pub fn compute(inputs: Inputs) !State {
    try validate(inputs);
    var result: State = undefined;
    result.canopy_air_temperature_k = inputs.solved_canopy_air_temperature_k;
    result.canopy_air_vapor_volume_fraction =
        inputs.solved_canopy_air_vapor_volume_fraction;
    result.surface_temperature_k = inputs.solved_surface_temperature_k;
    result.surface_temperature_c = result.surface_temperature_k - 273.15;
    if (!std.math.isFinite(result.surface_temperature_c))
        return error.NonFiniteStandingDeadFinalState;
    return result;
}

pub fn computeRuntimeSpecies(
    inputs: []const Inputs,
    scratch: []State,
    destination: []State,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.StandingDeadFinalStateDimensionMismatch;
    for (inputs, scratch) |species_inputs, *candidate|
        candidate.* = try compute(species_inputs);
    @memcpy(destination, scratch);
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteStandingDeadFinalStateInput;
    if (inputs.solved_canopy_air_temperature_k <= 0 or
        inputs.solved_canopy_air_vapor_volume_fraction < 0 or
        inputs.solved_surface_temperature_k <= 0)
        return error.InvalidStandingDeadFinalStateInput;
}

test "final standing-dead state preserves exact source commit order" {
    const result = try compute(.{
        .solved_canopy_air_temperature_k = 281,
        .solved_canopy_air_vapor_volume_fraction = 0.01,
        .solved_surface_temperature_k = 275,
    });
    try std.testing.expectEqual(@as(f64, 281), result.canopy_air_temperature_k);
    try std.testing.expectEqual(@as(f64, 0.01), result.canopy_air_vapor_volume_fraction);
    try std.testing.expectEqual(@as(f64, 275), result.surface_temperature_k);
    try std.testing.expectApproxEqAbs(@as(f64, 1.85), result.surface_temperature_c, 1e-12);
}

test "runtime species dimensions are unrestricted" {
    const inputs = [_]Inputs{
        .{
            .solved_canopy_air_temperature_k = 280,
            .solved_canopy_air_vapor_volume_fraction = 0.01,
            .solved_surface_temperature_k = 274,
        },
        .{
            .solved_canopy_air_temperature_k = 290,
            .solved_canopy_air_vapor_volume_fraction = 0.02,
            .solved_surface_temperature_k = 284,
        },
    };
    var scratch: [2]State = undefined;
    var destination: [2]State = undefined;
    try computeRuntimeSpecies(&inputs, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 274), destination[0].surface_temperature_k);
    try std.testing.expectEqual(@as(f64, 284), destination[1].surface_temperature_k);
}

test "later nonfinite species leaves destination unchanged" {
    const inputs = [_]Inputs{
        .{
            .solved_canopy_air_temperature_k = 280,
            .solved_canopy_air_vapor_volume_fraction = 0.01,
            .solved_surface_temperature_k = 274,
        },
        .{
            .solved_canopy_air_temperature_k = 290,
            .solved_canopy_air_vapor_volume_fraction = std.math.nan(f64),
            .solved_surface_temperature_k = 284,
        },
    };
    var scratch: [2]State = undefined;
    var destination: [2]State = undefined;
    destination[0] = try compute(inputs[0]);
    destination[1] = destination[0];
    destination[0].surface_temperature_k = 241;
    destination[1].surface_temperature_k = 242;
    try std.testing.expectError(
        error.NonFiniteStandingDeadFinalStateInput,
        computeRuntimeSpecies(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 241), destination[0].surface_temperature_k);
    try std.testing.expectEqual(@as(f64, 242), destination[1].surface_temperature_k);
}
