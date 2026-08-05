const std = @import("std");

pub const Parameters = struct {
    minimum_combustion_temperature_k: f64, // TCMBX
    methane_half_saturation_umol_mol: f64, // CCH4GK
    oxygen_half_saturation_umol_mol: f64, // COXYGK
    oxygen_g_per_g_c: f64 = 2.667,
    methane_combustion_energy_megajoules_g_c: f64, // GCBC4
};
pub const Forcing = struct {
    ground_surface_air_temperature_k: []const f64, // TKQGX
    methane_concentration_umol_mol: []const f64, // CH4Q
    oxygen_concentration_umol_mol: []const f64, // OXYQ
    canopy_oxygen_content_g_o: []const f64, // OXYC
};
pub const State = struct {
    carbon_dioxide_net_input_g_c_step: []f64, // XCNET
    methane_net_input_g_c_step: []f64, // XHNET
    oxygen_net_input_g_o_step: []f64, // XONET
    combustion_heat_megajoules_step: []f64, // HCBFG
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 10946--10970 for one runtime-indexed cell.
pub fn combustCell(cell: usize, parameters: Parameters, forcing: Forcing, state: State) !void {
    const cell_count = forcing.ground_surface_air_temperature_k.len;
    if (cell_count == 0 or cell >= cell_count or
        forcing.methane_concentration_umol_mol.len != cell_count or
        forcing.oxygen_concentration_umol_mol.len != cell_count or
        forcing.canopy_oxygen_content_g_o.len != cell_count or
        state.carbon_dioxide_net_input_g_c_step.len != cell_count or
        state.methane_net_input_g_c_step.len != cell_count or
        state.oxygen_net_input_g_o_step.len != cell_count or
        state.combustion_heat_megajoules_step.len != cell_count)
        return error.CanopyMethaneCombustionDimensionMismatch;
    inline for (std.meta.fields(Parameters)) |field|
        if (!std.math.isFinite(@field(parameters, field.name))) return error.InvalidCanopyMethaneCombustionInput;
    inline for (.{ forcing.ground_surface_air_temperature_k, forcing.methane_concentration_umol_mol, forcing.oxygen_concentration_umol_mol, forcing.canopy_oxygen_content_g_o, state.carbon_dioxide_net_input_g_c_step, state.methane_net_input_g_c_step, state.oxygen_net_input_g_o_step, state.combustion_heat_megajoules_step }) |values|
        if (!finiteSlice(values)) return error.InvalidCanopyMethaneCombustionInput;
    if (parameters.minimum_combustion_temperature_k <= 0 or
        parameters.methane_half_saturation_umol_mol < 0 or
        parameters.oxygen_half_saturation_umol_mol < 0 or
        parameters.oxygen_g_per_g_c <= 0 or
        parameters.methane_combustion_energy_megajoules_g_c < 0 or
        forcing.ground_surface_air_temperature_k[cell] <= 0 or
        forcing.methane_concentration_umol_mol[cell] < 0 or
        forcing.oxygen_concentration_umol_mol[cell] < 0 or
        forcing.canopy_oxygen_content_g_o[cell] < 0)
        return error.InvalidCanopyMethaneCombustionInput;

    var next_co2 = state.carbon_dioxide_net_input_g_c_step[cell];
    var next_ch4 = state.methane_net_input_g_c_step[cell];
    var next_o2 = state.oxygen_net_input_g_o_step[cell];
    var next_heat: f64 = 0;
    if (forcing.ground_surface_air_temperature_k[cell] > parameters.minimum_combustion_temperature_k) {
        const methane_denominator = forcing.methane_concentration_umol_mol[cell] + parameters.methane_half_saturation_umol_mol;
        const oxygen_denominator = forcing.oxygen_concentration_umol_mol[cell] + parameters.oxygen_half_saturation_umol_mol;
        if (methane_denominator <= 0 or oxygen_denominator <= 0)
            return error.ZeroCanopyMethaneCombustionDenominator;
        const concentration_limited_g_c = -state.methane_net_input_g_c_step[cell] *
            forcing.methane_concentration_umol_mol[cell] / methane_denominator *
            forcing.oxygen_concentration_umol_mol[cell] / oxygen_denominator;
        const oxygen_available_g_o = @max(0.0, forcing.canopy_oxygen_content_g_o[cell] - state.oxygen_net_input_g_o_step[cell]);
        const oxygen_limited_g_c = oxygen_available_g_o / parameters.oxygen_g_per_g_c;
        const combusted_g_c = @min(concentration_limited_g_c, oxygen_limited_g_c);
        next_o2 = next_o2 + combusted_g_c * parameters.oxygen_g_per_g_c;
        next_co2 = next_co2 - combusted_g_c;
        next_ch4 = next_ch4 + combusted_g_c;
        next_heat = combusted_g_c * parameters.methane_combustion_energy_megajoules_g_c;
    }
    inline for (.{ next_co2, next_ch4, next_o2, next_heat }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteCanopyMethaneCombustionResult;
    state.oxygen_net_input_g_o_step[cell] = next_o2;
    state.carbon_dioxide_net_input_g_c_step[cell] = next_co2;
    state.methane_net_input_g_c_step[cell] = next_ch4;
    state.combustion_heat_megajoules_step[cell] = next_heat;
}

test "REDIST canopy methane combustion preserves gate calculation and signs" {
    const temperatures = [_]f64{ 300, 310 };
    const methane = [_]f64{ 2, 2 };
    const oxygen = [_]f64{ 8, 8 };
    const oxygen_content = [_]f64{ 20, 20 };
    var co2 = [_]f64{ 1, 1 };
    var methane_flux = [_]f64{ -12, -12 };
    var oxygen_flux = [_]f64{ 2, 2 };
    var heat = [_]f64{ 9, 9 };
    const parameters = Parameters{ .minimum_combustion_temperature_k = 305, .methane_half_saturation_umol_mol = 1, .oxygen_half_saturation_umol_mol = 2, .methane_combustion_energy_megajoules_g_c = 0.05 };
    const forcing = Forcing{ .ground_surface_air_temperature_k = &temperatures, .methane_concentration_umol_mol = &methane, .oxygen_concentration_umol_mol = &oxygen, .canopy_oxygen_content_g_o = &oxygen_content };
    const state = State{ .carbon_dioxide_net_input_g_c_step = &co2, .methane_net_input_g_c_step = &methane_flux, .oxygen_net_input_g_o_step = &oxygen_flux, .combustion_heat_megajoules_step = &heat };
    try combustCell(0, parameters, forcing, state);
    try std.testing.expectEqual(@as(f64, 1), co2[0]);
    try std.testing.expectEqual(@as(f64, -12), methane_flux[0]);
    try std.testing.expectEqual(@as(f64, 2), oxygen_flux[0]);
    try std.testing.expectEqual(@as(f64, 0), heat[0]);
    try combustCell(1, parameters, forcing, state);
    const combusted: f64 = 12 * 2.0 / 3.0 * 8.0 / 10.0;
    try std.testing.expectApproxEqAbs(1 - combusted, co2[1], 1e-12);
    try std.testing.expectApproxEqAbs(-12 + combusted, methane_flux[1], 1e-12);
    try std.testing.expectApproxEqAbs(2 + combusted * 2.667, oxygen_flux[1], 1e-12);
    try std.testing.expectApproxEqAbs(combusted * 0.05, heat[1], 1e-12);
}

test "REDIST canopy methane late overflow is atomic" {
    const temperatures = [_]f64{310};
    const methane = [_]f64{2};
    const oxygen = [_]f64{8};
    const oxygen_content = [_]f64{std.math.floatMax(f64)};
    var co2 = [_]f64{1};
    var methane_flux = [_]f64{-std.math.floatMax(f64)};
    var oxygen_flux = [_]f64{2};
    var heat = [_]f64{9};
    try std.testing.expectError(error.NonFiniteCanopyMethaneCombustionResult, combustCell(0, .{ .minimum_combustion_temperature_k = 305, .methane_half_saturation_umol_mol = 0, .oxygen_half_saturation_umol_mol = 0, .methane_combustion_energy_megajoules_g_c = std.math.floatMax(f64) }, .{ .ground_surface_air_temperature_k = &temperatures, .methane_concentration_umol_mol = &methane, .oxygen_concentration_umol_mol = &oxygen, .canopy_oxygen_content_g_o = &oxygen_content }, .{ .carbon_dioxide_net_input_g_c_step = &co2, .methane_net_input_g_c_step = &methane_flux, .oxygen_net_input_g_o_step = &oxygen_flux, .combustion_heat_megajoules_step = &heat }));
    try std.testing.expectEqual(@as(f64, 1), co2[0]);
    try std.testing.expectEqual(-std.math.floatMax(f64), methane_flux[0]);
    try std.testing.expectEqual(@as(f64, 2), oxygen_flux[0]);
    try std.testing.expectEqual(@as(f64, 9), heat[0]);
}
