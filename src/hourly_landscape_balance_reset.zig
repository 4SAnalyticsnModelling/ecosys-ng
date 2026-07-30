const std = @import("std");

pub const State = struct {
    soil_water_volume_m3: f64,
    soil_heat_mj: f64,
    soil_oxygen_mass_g_o: f64,
    total_hydrogen_gas_mass_g_h: f64,
    soil_sediment_mass_mg: f64,
    residue_carbon_mass_g_c: f64,
    soil_organic_carbon_mass_g_c: f64,
    soil_co2_carbon_mass_g_c: f64,
    residue_nitrogen_mass_g_n: f64,
    soil_organic_nitrogen_mass_g_n: f64,
    soil_dinitrogen_mass_g_n: f64,
    residue_phosphorus_mass_g_p: f64,
    soil_organic_phosphorus_mass_g_p: f64,
    soil_ammonium_mass_g_n: f64,
    soil_nitrate_mass_g_n: f64,
    soil_phosphate_mass_g_p: f64,
    soil_ion_amount_mol: f64,
    plant_carbon_balance_g_c: f64,
    plant_nitrogen_balance_g_n: f64,
    plant_phosphorus_balance_g_p: f64,
};

/// Exact HOUR1 lines 136--155 reset in source assignment order.
pub fn apply(state: *State) void {
    state.soil_water_volume_m3 = 0;
    state.soil_heat_mj = 0;
    state.soil_oxygen_mass_g_o = 0;
    state.total_hydrogen_gas_mass_g_h = 0;
    state.soil_sediment_mass_mg = 0;
    state.residue_carbon_mass_g_c = 0;
    state.soil_organic_carbon_mass_g_c = 0;
    state.soil_co2_carbon_mass_g_c = 0;
    state.residue_nitrogen_mass_g_n = 0;
    state.soil_organic_nitrogen_mass_g_n = 0;
    state.soil_dinitrogen_mass_g_n = 0;
    state.residue_phosphorus_mass_g_p = 0;
    state.soil_organic_phosphorus_mass_g_p = 0;
    state.soil_ammonium_mass_g_n = 0;
    state.soil_nitrate_mass_g_n = 0;
    state.soil_phosphate_mass_g_p = 0;
    state.soil_ion_amount_mol = 0;
    state.plant_carbon_balance_g_c = 0;
    state.plant_nitrogen_balance_g_n = 0;
    state.plant_phosphorus_balance_g_p = 0;
}

fn filled(value: f64) State {
    var state: State = undefined;
    inline for (@typeInfo(State).@"struct".fields) |field|
        @field(state, field.name) = value;
    return state;
}

test "HOUR1 resets all twenty landscape balances" {
    var state = filled(7);
    apply(&state);
    try std.testing.expectEqualDeep(filled(0), state);
}

test "reset deliberately replaces invalid prior accumulator values" {
    var state = filled(std.math.nan(f64));
    apply(&state);
    inline for (@typeInfo(State).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(state, field.name));
}
