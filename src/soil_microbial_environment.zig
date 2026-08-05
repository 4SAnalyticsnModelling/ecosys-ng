const std = @import("std");

pub const Inputs = struct {
    total_water_fraction: f64,
    porosity_fraction: f64,
    field_capacity_fraction: f64,
    inactive_water_fraction: f64,
    matrix_bulk_volume_m3: f64,
    soil_temperature_k: f64,
    thermal_adaptation_offset_k: f64,
    aqueous_oxygen_concentration_g_o_per_m3: f64,
    soil_organic_carbon_g_c: f64,
    soil_dry_mass_megagrams: f64,
};

pub const Result = struct {
    biologically_active_water_fraction: f64,
    biologically_active_water_m3: f64,
    growth_temperature_response: f64,
    maintenance_temperature_response: f64,
    fermentation_oxygen_inhibition_fraction: f64,
    density_limited_soil_organic_carbon_g_c: f64,
};

/// NITRO.F 287--326 for a mineral-soil layer. This retains the source
/// concentration transform and thermally adapted Arrhenius functions while
/// making every formerly global coordinate explicit.
pub fn calculate(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSoilMicrobialEnvironment;
    }
    if (inputs.total_water_fraction < 0 or
        inputs.porosity_fraction < 0 or inputs.porosity_fraction > 1 or
        inputs.field_capacity_fraction < 0 or inputs.field_capacity_fraction > 1 or
        inputs.inactive_water_fraction < 0 or
        inputs.matrix_bulk_volume_m3 < 0 or
        inputs.soil_temperature_k <= 0 or
        inputs.soil_temperature_k + inputs.thermal_adaptation_offset_k <= 0 or
        inputs.aqueous_oxygen_concentration_g_o_per_m3 < 0 or
        inputs.soil_organic_carbon_g_c < 0 or
        inputs.soil_dry_mass_megagrams < 0)
        return error.InvalidSoilMicrobialEnvironment;

    const retained_water_fraction = @min(
        @max(0.75 * inputs.porosity_fraction, inputs.field_capacity_fraction),
        inputs.total_water_fraction,
    );
    const active_water_fraction = @max(
        0,
        retained_water_fraction - inputs.inactive_water_fraction,
    );
    const active_water_m3 =
        active_water_fraction / (1 + active_water_fraction) *
        inputs.matrix_bulk_volume_m3;

    const adapted_temperature_k =
        inputs.soil_temperature_k + inputs.thermal_adaptation_offset_k;
    const rt_j_per_mol = 8.3143 * adapted_temperature_k;
    const entropy_temperature_j_per_mol = 710 * adapted_temperature_k;
    const growth_inactivation =
        1 +
        @exp((197500 - entropy_temperature_j_per_mol) / rt_j_per_mol) +
        @exp((entropy_temperature_j_per_mol - 222500) / rt_j_per_mol);
    const growth_response =
        @exp(25.229 - 62500 / rt_j_per_mol) / growth_inactivation;
    const maintenance_inactivation =
        1 + @exp((197500 - entropy_temperature_j_per_mol) / rt_j_per_mol);
    const maintenance_response = @min(
        1.0e3,
        @exp(25.216 - 62500 / rt_j_per_mol) / maintenance_inactivation,
    );
    const oxygen_inhibition =
        1 - 1 / (1 + @exp(-inputs.aqueous_oxygen_concentration_g_o_per_m3 + 2.5));
    const density_limited_carbon_g_c = @min(
        1.0e5 * inputs.soil_dry_mass_megagrams,
        inputs.soil_organic_carbon_g_c,
    );

    const result: Result = .{
        .biologically_active_water_fraction = active_water_fraction,
        .biologically_active_water_m3 = active_water_m3,
        .growth_temperature_response = growth_response,
        .maintenance_temperature_response = maintenance_response,
        .fermentation_oxygen_inhibition_fraction = oxygen_inhibition,
        .density_limited_soil_organic_carbon_g_c = density_limited_carbon_g_c,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSoilMicrobialEnvironment;
    }
    return result;
}

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    biologically_active_water_fraction: []f64,
    biologically_active_water_m3: []f64,
    growth_temperature_response: []f64,
    maintenance_temperature_response: []f64,
    fermentation_oxygen_inhibition_fraction: []f64,
    density_limited_soil_organic_carbon_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.InvalidSoilMicrobialLayerCount;
        var state: State = undefined;
        state.allocator = allocator;
        state.layer_count = layer_count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64) {
                    if (allocated > 0) {
                        allocated -= 1;
                        allocator.free(@field(state, field.name));
                    }
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(state, field.name) = try allocator.alloc(f64, layer_count);
                @memset(@field(state, field.name), 0);
                allocated += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn set(self: *State, layer: usize, result: Result) !void {
        if (layer >= self.layer_count) return error.SoilMicrobialLayerOutOfRange;
        inline for (@typeInfo(Result).@"struct".fields) |field|
            @field(self, field.name)[layer] = @field(result, field.name);
    }
};

test "mineral layer environment reproduces NITRO water temperature oxygen and carbon equations" {
    const inputs: Inputs = .{
        .total_water_fraction = 0.4,
        .porosity_fraction = 0.5,
        .field_capacity_fraction = 0.3,
        .inactive_water_fraction = 0.1,
        .matrix_bulk_volume_m3 = 2,
        .soil_temperature_k = 293.15,
        .thermal_adaptation_offset_k = 0,
        .aqueous_oxygen_concentration_g_o_per_m3 = 1.5,
        .soil_organic_carbon_g_c = 300_000,
        .soil_dry_mass_megagrams = 2,
    };
    const result = try calculate(inputs);
    const expected_active_fraction: f64 = 0.275;
    try std.testing.expectApproxEqAbs(expected_active_fraction, result.biologically_active_water_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(expected_active_fraction / (1 + expected_active_fraction) * 2, result.biologically_active_water_m3, 1e-15);
    const rt: f64 = 8.3143 * 293.15;
    const entropy_temperature: f64 = 710 * 293.15;
    try std.testing.expectApproxEqAbs(
        @exp(25.229 - 62500 / rt) /
            (1 + @exp((197500 - entropy_temperature) / rt) +
                @exp((entropy_temperature - 222500) / rt)),
        result.growth_temperature_response,
        5e-15,
    );
    try std.testing.expectApproxEqAbs(
        @min(1.0e3, @exp(25.216 - 62500 / rt) /
            (1 + @exp((197500 - entropy_temperature) / rt))),
        result.maintenance_temperature_response,
        5e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1) - @as(f64, 1) / (1 + @exp(1.0)),
        result.fermentation_oxygen_inhibition_fraction,
        1e-15,
    );
    try std.testing.expectEqual(@as(f64, 200_000), result.density_limited_soil_organic_carbon_g_c);
}

test "runtime state owns arbitrary layer count on heap" {
    var state = try State.init(std.testing.allocator, 37);
    defer state.deinit();
    const result = try calculate(.{
        .total_water_fraction = 0.2,
        .porosity_fraction = 0.4,
        .field_capacity_fraction = 0.25,
        .inactive_water_fraction = 0.05,
        .matrix_bulk_volume_m3 = 1,
        .soil_temperature_k = 280,
        .thermal_adaptation_offset_k = 2,
        .aqueous_oxygen_concentration_g_o_per_m3 = 0,
        .soil_organic_carbon_g_c = 10,
        .soil_dry_mass_megagrams = 1,
    });
    try state.set(36, result);
    try std.testing.expectEqual(result.biologically_active_water_m3, state.biologically_active_water_m3[36]);
    try std.testing.expectError(error.SoilMicrobialLayerOutOfRange, state.set(37, result));
}

test "invalid microbial environment fails explicitly" {
    try std.testing.expectError(error.InvalidSoilMicrobialEnvironment, calculate(.{
        .total_water_fraction = 0.2,
        .porosity_fraction = 0.4,
        .field_capacity_fraction = 0.25,
        .inactive_water_fraction = 0.05,
        .matrix_bulk_volume_m3 = 1,
        .soil_temperature_k = 280,
        .thermal_adaptation_offset_k = -280,
        .aqueous_oxygen_concentration_g_o_per_m3 = 0,
        .soil_organic_carbon_g_c = 10,
        .soil_dry_mass_megagrams = 1,
    }));
}
