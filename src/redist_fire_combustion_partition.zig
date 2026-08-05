const std = @import("std");

pub const Accumulators = struct {
    carbon_dioxide_emission_g_c: f64,
    methane_emission_g_c: f64,
    charcoal_carbon_g_c: f64,
    gaseous_nitrogen_g_n: f64,
    gaseous_phosphorus_g_p: f64,
    mineral_nitrogen_g_n: f64,
    mineral_phosphorus_g_p: f64,
};
pub const Parameters = struct {
    mineral_n_fraction_at_full_combustion: f64,
    mineral_n_fraction_at_zero_combustion: f64,
    mineral_p_fraction_at_full_combustion: f64,
    mineral_p_fraction_at_zero_combustion: f64,
    aerobic_energy_megajoules_g_c: f64,
    anaerobic_energy_megajoules_g_c: f64,
    methane_energy_megajoules_g_c: f64,
};
pub const Inputs = struct {
    total_soil_combustion_g_c: []const f64,
    oxygen_limited_combustion_g_c: []const f64,
    methane_unlimited_combustion_g_c: []const f64,
    methane_limited_combustion_g_c: []const f64,
    uncombusted_fraction: []const f64,
    minimum_combustion_g_c: f64,
};
pub const Workspace = struct {
    mineral_n_fraction: []f64,
    gaseous_n_fraction: []f64,
    mineral_p_fraction: []f64,
    gaseous_p_fraction: []f64,
    charcoal_fraction: []f64,
    aerobic_fraction: []f64,
    anaerobic_fraction: []f64,
};
pub const State = struct { heat_release_megajoules: []f64, accumulators: *Accumulators, workspace: Workspace };

fn finiteSlice(v: []const f64) bool {
    for (v) |x| if (!std.math.isFinite(x)) return false;
    return true;
}
fn finiteStruct(v: anytype) bool {
    inline for (std.meta.fields(@TypeOf(v))) |f| if (!std.math.isFinite(@field(v, f.name))) return false;
    return true;
}

/// Direct translation of REDIST 10684 and 10700--10719. Unsafe legacy divisions fail explicitly.
pub fn prepare(fire_flag: u8, parameters: Parameters, inputs: Inputs, state: State) !void {
    const n = inputs.total_soil_combustion_g_c.len;
    if (n == 0 or inputs.oxygen_limited_combustion_g_c.len != n or inputs.methane_unlimited_combustion_g_c.len != n or
        inputs.methane_limited_combustion_g_c.len != n or inputs.uncombusted_fraction.len != n or state.heat_release_megajoules.len != n)
        return error.FireCombustionPartitionDimensionMismatch;
    inline for (std.meta.fields(Workspace)) |f| if (@field(state.workspace, f.name).len != n) return error.FireCombustionPartitionDimensionMismatch;
    inline for (.{ inputs.total_soil_combustion_g_c, inputs.oxygen_limited_combustion_g_c, inputs.methane_unlimited_combustion_g_c, inputs.methane_limited_combustion_g_c, inputs.uncombusted_fraction, state.heat_release_megajoules }) |v| if (!finiteSlice(v)) return error.InvalidFireCombustionPartitionInput;
    if (!finiteStruct(parameters) or !finiteStruct(state.accumulators.*) or !std.math.isFinite(inputs.minimum_combustion_g_c)) return error.InvalidFireCombustionPartitionInput;
    if (inputs.minimum_combustion_g_c < 0 or
        parameters.mineral_n_fraction_at_full_combustion < 0 or parameters.mineral_n_fraction_at_full_combustion > 1 or
        parameters.mineral_n_fraction_at_zero_combustion < 0 or parameters.mineral_n_fraction_at_zero_combustion > 1 or
        parameters.mineral_p_fraction_at_full_combustion < 0 or parameters.mineral_p_fraction_at_full_combustion > 1 or
        parameters.mineral_p_fraction_at_zero_combustion < 0 or parameters.mineral_p_fraction_at_zero_combustion > 1)
        return error.InvalidFireCombustionPartitionInput;
    for (0..n) |layer| {
        if (inputs.total_soil_combustion_g_c[layer] < 0 or inputs.oxygen_limited_combustion_g_c[layer] < 0 or
            inputs.methane_unlimited_combustion_g_c[layer] < 0 or inputs.methane_limited_combustion_g_c[layer] < 0 or
            inputs.uncombusted_fraction[layer] < 0 or inputs.uncombusted_fraction[layer] > 1)
            return error.InvalidFireCombustionPartitionInput;
    }
    if (fire_flag != 1) return;

    for (0..n) |layer| if (inputs.total_soil_combustion_g_c[layer] > inputs.minimum_combustion_g_c) {
        const denominator = inputs.oxygen_limited_combustion_g_c[layer] + inputs.methane_unlimited_combustion_g_c[layer];
        if (denominator == 0) return error.ZeroFireCombustionGasPartition;
        const mineral_n = parameters.mineral_n_fraction_at_zero_combustion + (parameters.mineral_n_fraction_at_full_combustion - parameters.mineral_n_fraction_at_zero_combustion) * (1.0 - inputs.uncombusted_fraction[layer]);
        const mineral_p = parameters.mineral_p_fraction_at_zero_combustion + (parameters.mineral_p_fraction_at_full_combustion - parameters.mineral_p_fraction_at_zero_combustion) * (1.0 - inputs.uncombusted_fraction[layer]);
        const charcoal = (inputs.total_soil_combustion_g_c[layer] - inputs.oxygen_limited_combustion_g_c[layer] - inputs.methane_unlimited_combustion_g_c[layer]) / inputs.total_soil_combustion_g_c[layer];
        const aerobic = inputs.oxygen_limited_combustion_g_c[layer] / denominator;
        const anaerobic = inputs.methane_unlimited_combustion_g_c[layer] / denominator;
        const heat = state.heat_release_megajoules[layer] + inputs.oxygen_limited_combustion_g_c[layer] * parameters.aerobic_energy_megajoules_g_c + inputs.methane_unlimited_combustion_g_c[layer] * parameters.anaerobic_energy_megajoules_g_c + inputs.methane_limited_combustion_g_c[layer] * parameters.methane_energy_megajoules_g_c;
        inline for (.{ mineral_n, 1.0 - mineral_n, mineral_p, 1.0 - mineral_p, charcoal, aerobic, anaerobic, heat }) |x| if (!std.math.isFinite(x)) return error.NonFiniteFireCombustionPartitionResult;
        inline for (.{ mineral_n, mineral_p, charcoal, aerobic, anaerobic }) |fraction|
            if (fraction < 0 or fraction > 1) return error.InvalidFireCombustionPartitionResult;
    };
    state.accumulators.* = .{ .carbon_dioxide_emission_g_c = 0, .methane_emission_g_c = 0, .charcoal_carbon_g_c = 0, .gaseous_nitrogen_g_n = 0, .gaseous_phosphorus_g_p = 0, .mineral_nitrogen_g_n = 0, .mineral_phosphorus_g_p = 0 };
    for (0..n) |layer| if (inputs.total_soil_combustion_g_c[layer] > inputs.minimum_combustion_g_c) {
        const denominator = inputs.oxygen_limited_combustion_g_c[layer] + inputs.methane_unlimited_combustion_g_c[layer];
        const mineral_n = parameters.mineral_n_fraction_at_zero_combustion + (parameters.mineral_n_fraction_at_full_combustion - parameters.mineral_n_fraction_at_zero_combustion) * (1.0 - inputs.uncombusted_fraction[layer]);
        const mineral_p = parameters.mineral_p_fraction_at_zero_combustion + (parameters.mineral_p_fraction_at_full_combustion - parameters.mineral_p_fraction_at_zero_combustion) * (1.0 - inputs.uncombusted_fraction[layer]);
        state.workspace.mineral_n_fraction[layer] = mineral_n;
        state.workspace.gaseous_n_fraction[layer] = 1.0 - mineral_n;
        state.workspace.mineral_p_fraction[layer] = mineral_p;
        state.workspace.gaseous_p_fraction[layer] = 1.0 - mineral_p;
        state.workspace.charcoal_fraction[layer] = (inputs.total_soil_combustion_g_c[layer] - inputs.oxygen_limited_combustion_g_c[layer] - inputs.methane_unlimited_combustion_g_c[layer]) / inputs.total_soil_combustion_g_c[layer];
        state.workspace.aerobic_fraction[layer] = inputs.oxygen_limited_combustion_g_c[layer] / denominator;
        state.workspace.anaerobic_fraction[layer] = inputs.methane_unlimited_combustion_g_c[layer] / denominator;
        state.heat_release_megajoules[layer] = state.heat_release_megajoules[layer] + inputs.oxygen_limited_combustion_g_c[layer] * parameters.aerobic_energy_megajoules_g_c + inputs.methane_unlimited_combustion_g_c[layer] * parameters.anaerobic_energy_megajoules_g_c + inputs.methane_limited_combustion_g_c[layer] * parameters.methane_energy_megajoules_g_c;
    };
}

test "REDIST fire partition resets ledgers and conserves N P fractions" {
    const total = [_]f64{10};
    const aerobic = [_]f64{3};
    const anaerobic = [_]f64{2};
    const methane = [_]f64{1};
    const unburned = [_]f64{0.25};
    var heat = [_]f64{1};
    var storage: [7][1]f64 = .{.{0}} ** 7;
    var a = Accumulators{ .carbon_dioxide_emission_g_c = 1, .methane_emission_g_c = 1, .charcoal_carbon_g_c = 1, .gaseous_nitrogen_g_n = 1, .gaseous_phosphorus_g_p = 1, .mineral_nitrogen_g_n = 1, .mineral_phosphorus_g_p = 1 };
    try prepare(1, .{ .mineral_n_fraction_at_full_combustion = 0.8, .mineral_n_fraction_at_zero_combustion = 0.4, .mineral_p_fraction_at_full_combustion = 0.6, .mineral_p_fraction_at_zero_combustion = 0.2, .aerobic_energy_megajoules_g_c = 2, .anaerobic_energy_megajoules_g_c = 3, .methane_energy_megajoules_g_c = 4 }, .{ .total_soil_combustion_g_c = &total, .oxygen_limited_combustion_g_c = &aerobic, .methane_unlimited_combustion_g_c = &anaerobic, .methane_limited_combustion_g_c = &methane, .uncombusted_fraction = &unburned, .minimum_combustion_g_c = 0 }, .{ .heat_release_megajoules = &heat, .accumulators = &a, .workspace = .{ .mineral_n_fraction = &storage[0], .gaseous_n_fraction = &storage[1], .mineral_p_fraction = &storage[2], .gaseous_p_fraction = &storage[3], .charcoal_fraction = &storage[4], .aerobic_fraction = &storage[5], .anaerobic_fraction = &storage[6] } });
    try std.testing.expectEqual(@as(f64, 1), storage[0][0] + storage[1][0]);
    try std.testing.expectEqual(@as(f64, 1), storage[2][0] + storage[3][0]);
    try std.testing.expectEqual(@as(f64, 1), storage[5][0] + storage[6][0]);
    try std.testing.expectEqual(@as(f64, 17), heat[0]);
    try std.testing.expectEqual(@as(f64, 0), a.carbon_dioxide_emission_g_c);
}
test "REDIST fire partition failure is atomic" {
    const total = [_]f64{1};
    const zero = [_]f64{0};
    var heat = [_]f64{2};
    var storage: [7][1]f64 = .{.{9}} ** 7;
    var a = Accumulators{ .carbon_dioxide_emission_g_c = 1, .methane_emission_g_c = 1, .charcoal_carbon_g_c = 1, .gaseous_nitrogen_g_n = 1, .gaseous_phosphorus_g_p = 1, .mineral_nitrogen_g_n = 1, .mineral_phosphorus_g_p = 1 };
    try std.testing.expectError(error.ZeroFireCombustionGasPartition, prepare(1, .{ .mineral_n_fraction_at_full_combustion = 1, .mineral_n_fraction_at_zero_combustion = 0, .mineral_p_fraction_at_full_combustion = 1, .mineral_p_fraction_at_zero_combustion = 0, .aerobic_energy_megajoules_g_c = 1, .anaerobic_energy_megajoules_g_c = 1, .methane_energy_megajoules_g_c = 1 }, .{ .total_soil_combustion_g_c = &total, .oxygen_limited_combustion_g_c = &zero, .methane_unlimited_combustion_g_c = &zero, .methane_limited_combustion_g_c = &zero, .uncombusted_fraction = &zero, .minimum_combustion_g_c = 0 }, .{ .heat_release_megajoules = &heat, .accumulators = &a, .workspace = .{ .mineral_n_fraction = &storage[0], .gaseous_n_fraction = &storage[1], .mineral_p_fraction = &storage[2], .gaseous_p_fraction = &storage[3], .charcoal_fraction = &storage[4], .aerobic_fraction = &storage[5], .anaerobic_fraction = &storage[6] } }));
    try std.testing.expectEqual(@as(f64, 2), heat[0]);
    try std.testing.expectEqual(@as(f64, 1), a.carbon_dioxide_emission_g_c);
}
