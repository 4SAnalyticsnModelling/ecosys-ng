const std = @import("std");

pub const State = struct {
    carbon_dioxide_g_c: f64, // COX
    methane_g_c: f64, // CHX
    nitrogen_oxide_g_n: f64, // ZOX
    phosphorus_oxide_g_p: f64, // POX
    charcoal_g_c: f64, // COR
    ammonium_g_n: f64, // Z4M
    dihydrogen_phosphate_g_p: f64, // P4M
};

/// EXTRACT lines 450--454. Resets root-fire combustion scratch accumulators
/// in exact source assignment order before layer traversal.
pub fn apply(state: *State) void {
    state.carbon_dioxide_g_c = 0;
    state.methane_g_c = 0;
    state.nitrogen_oxide_g_n = 0;
    state.phosphorus_oxide_g_p = 0;
    state.charcoal_g_c = 0;
    state.ammonium_g_n = 0;
    state.dihydrogen_phosphate_g_p = 0;
}

fn filled(value: f64) State {
    var state: State = undefined;
    inline for (@typeInfo(State).@"struct".fields) |field|
        @field(state, field.name) = value;
    return state;
}

test "EXTRACT root fire scratch accumulators reset to exact zero" {
    var state = filled(7);
    apply(&state);
    try std.testing.expectEqualDeep(filled(0), state);
}

test "root fire reset overwrites invalid prior scratch values" {
    var state = filled(std.math.nan(f64));
    apply(&state);
    inline for (@typeInfo(State).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(state, field.name));
}
