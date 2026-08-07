const std = @import("std");

pub const RetentionMode = enum {
    supplied,
    estimated,
};

pub const Inputs = struct {
    is_first_simulation_day: bool,
    is_first_hour: bool,
    is_initialization_year: bool,
    retention_mode: RetentionMode,
    water_initialization_code: f64,
    ice_initialization_code: f64,
    layer_depth_m: f64,
    water_table_depth_m: f64,
    porosity_m3_m3: f64,
    field_capacity_m3_m3: f64,
    wilting_point_m3_m3: f64,
    micropore_volume_m3: f64,
    macropore_volume_m3: f64,
    mineral_heat_capacity_megajoules_k: f64,
};

pub const State = struct {
    water_content_m3_m3: f64,
    ice_content_m3_m3: f64,
    micropore_water_m3: f64,
    previous_micropore_water_m3: f64,
    macropore_water_m3: f64,
    micropore_ice_m3: f64,
    macropore_ice_m3: f64,
    heat_capacity_megajoules_k: f64,
    previous_water_content_m3_m3: f64,
    previous_ice_content_m3_m3: f64,
};

/// `hour1.f` lines 2131--2163 for one runtime soil layer. Sentinel codes retain
/// source equality branches; codes strictly between zero and one preserve
/// the caller's existing concentration.
pub fn apply(inputs: Inputs, state: *State) !void {
    try validate(inputs, state.*);
    if (!(inputs.is_first_simulation_day and inputs.is_first_hour and
        inputs.is_initialization_year)) return;
    if (inputs.water_initialization_code > 1.0 or
        inputs.layer_depth_m >= inputs.water_table_depth_m)
        state.water_content_m3_m3 = inputs.porosity_m3_m3
    else if (inputs.water_initialization_code == 1.0)
        state.water_content_m3_m3 = inputs.field_capacity_m3_m3
    else if (inputs.water_initialization_code == 0.0)
        state.water_content_m3_m3 = inputs.wilting_point_m3_m3
    else if (inputs.water_initialization_code < 0.0)
        state.water_content_m3_m3 = 0.0;

    const remaining_porosity =
        inputs.porosity_m3_m3 - state.water_content_m3_m3;
    if (inputs.ice_initialization_code > 1.0 or
        inputs.layer_depth_m >= inputs.water_table_depth_m)
        state.ice_content_m3_m3 =
            @max(0.0, @min(inputs.porosity_m3_m3, remaining_porosity))
    else if (inputs.ice_initialization_code == 1.0)
        state.ice_content_m3_m3 =
            @max(0.0, @min(inputs.field_capacity_m3_m3, remaining_porosity))
    else if (inputs.ice_initialization_code == 0.0)
        state.ice_content_m3_m3 =
            @max(0.0, @min(inputs.wilting_point_m3_m3, remaining_porosity))
    else if (inputs.ice_initialization_code < 0.0)
        state.ice_content_m3_m3 = 0.0;

    if (inputs.retention_mode == .estimated) {
        state.micropore_water_m3 =
            state.water_content_m3_m3 * inputs.micropore_volume_m3;
        state.previous_micropore_water_m3 = state.micropore_water_m3;
        state.macropore_water_m3 =
            state.water_content_m3_m3 * inputs.macropore_volume_m3;
        state.micropore_ice_m3 =
            state.ice_content_m3_m3 * inputs.micropore_volume_m3;
        state.macropore_ice_m3 =
            state.ice_content_m3_m3 * inputs.macropore_volume_m3;
        state.heat_capacity_megajoules_k = inputs.mineral_heat_capacity_megajoules_k +
            4.19 * (state.micropore_water_m3 + state.macropore_water_m3) +
            1.9274 * (state.micropore_ice_m3 + state.macropore_ice_m3);
        state.previous_water_content_m3_m3 = state.water_content_m3_m3;
        state.previous_ice_content_m3_m3 = state.ice_content_m3_m3;
    }
}

fn validate(inputs: Inputs, state: State) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSoilInitializationInput;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state, field.name)))
            return error.NonFiniteSoilInitializationState;
    if (inputs.porosity_m3_m3 < 0 or inputs.field_capacity_m3_m3 < 0 or
        inputs.wilting_point_m3_m3 < 0 or inputs.micropore_volume_m3 < 0 or
        inputs.macropore_volume_m3 < 0 or inputs.mineral_heat_capacity_megajoules_k < 0)
        return error.InvalidSoilInitializationInput;
}

test "first-run codes initialize water then ice and publish volumes" {
    var state: State = std.mem.zeroes(State);
    try apply(.{
        .is_first_simulation_day = true,
        .is_first_hour = true,
        .is_initialization_year = true,
        .retention_mode = .estimated,
        .water_initialization_code = 1,
        .ice_initialization_code = 1,
        .layer_depth_m = 0.1,
        .water_table_depth_m = 1,
        .porosity_m3_m3 = 0.5,
        .field_capacity_m3_m3 = 0.3,
        .wilting_point_m3_m3 = 0.1,
        .micropore_volume_m3 = 2,
        .macropore_volume_m3 = 1,
        .mineral_heat_capacity_megajoules_k = 5,
    }, &state);
    try std.testing.expectEqual(@as(f64, 0.3), state.water_content_m3_m3);
    try std.testing.expectEqual(@as(f64, 0.2), state.ice_content_m3_m3);
    try std.testing.expectEqual(@as(f64, 0.6), state.micropore_water_m3);
    try std.testing.expectApproxEqAbs(
        @as(f64, 5 + 4.19 * 0.9 + 1.9274 * 0.6),
        state.heat_capacity_megajoules_k,
        1e-14,
    );
}

test "noninitial execution leaves state unchanged" {
    var state: State = std.mem.zeroes(State);
    state.water_content_m3_m3 = 0.42;
    try apply(.{
        .is_first_simulation_day = false,
        .is_first_hour = true,
        .is_initialization_year = true,
        .retention_mode = .estimated,
        .water_initialization_code = 1,
        .ice_initialization_code = 1,
        .layer_depth_m = 0,
        .water_table_depth_m = 1,
        .porosity_m3_m3 = 0.5,
        .field_capacity_m3_m3 = 0.3,
        .wilting_point_m3_m3 = 0.1,
        .micropore_volume_m3 = 1,
        .macropore_volume_m3 = 1,
        .mineral_heat_capacity_megajoules_k = 1,
    }, &state);
    try std.testing.expectEqual(@as(f64, 0.42), state.water_content_m3_m3);
}
