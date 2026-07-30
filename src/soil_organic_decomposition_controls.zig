const std = @import("std");

pub const Inputs = struct {
    surface_layer: []const bool,
    total_active_biomass_carbon_g_c: []const f64,
    total_active_biomass_nitrogen_g_n: []const f64,
    total_active_biomass_phosphorus_g_p: []const f64,
    maximum_active_biomass_nitrogen_g_n: []const f64,
    maximum_active_biomass_phosphorus_g_p: []const f64,
    colonized_soil_organic_carbon_g_c: []const f64,
    microbial_activity_respiration_g_c: []const f64,
    biologically_active_water_m3: []const f64,
    soil_layer_mass_Mg: []const f64,
    soil_layer_volume_m3: []const f64,
    dissolved_activity_carbon_concentration_g_c_m3: []const f64,
    surface_zero_activity_half_saturation_g_c_m3: f64,
    surface_activity_inhibition_rate_g_c_m3_h: f64,
    soil_zero_activity_half_saturation_g_c_m3: f64,
    soil_activity_inhibition_rate_g_c_m3_h: f64,
    dissolved_carbon_product_inhibition_g_c_m3: f64,
    minimum_nutrient_limitation_fraction: f64 = 0.5,
    maximum_activity_concentration_g_c_m3_h: f64 = 100_000,
    negligible_biomass_g_c: f64,
    negligible_water_m3: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    complex_count: usize,
    relative_nitrogen_status: []f64,
    relative_phosphorus_status: []f64,
    nitrogen_limitation_fraction: []f64,
    phosphorus_limitation_fraction: []f64,
    microbial_activity_concentration_g_c_m3_h: []f64,
    decomposition_half_saturation_g_c_m3: []f64,
    effective_substrate_concentration_sqrt_g_c_per_Mg: []f64,
    substrate_density_response: []f64,
    dissolved_carbon_product_response: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize) !State {
        if (complex_count == 0) return error.InvalidDecompositionControlDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.complex_count = complex_count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(state, field.name) = try allocator.alloc(f64, complex_count);
            @memset(@field(state, field.name), 0);
            allocated += 1;
        };
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// Exact NITRO.F 3182--3240 nutrient and concentration controls on SOC decomposition.
/// Complex traversal and the nutrient-before-concentration order follow the source.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([9]f64, state.complex_count);
    defer state.allocator.free(temporary);
    for (0..state.complex_count) |complex| {
        var relative_nitrogen: f64 = 1;
        var relative_phosphorus: f64 = 1;
        var nitrogen_limitation: f64 = 1;
        var phosphorus_limitation: f64 = 1;
        if (inputs.total_active_biomass_carbon_g_c[complex] > inputs.negligible_biomass_g_c) {
            if (inputs.maximum_active_biomass_nitrogen_g_n[complex] == 0 or
                inputs.maximum_active_biomass_phosphorus_g_p[complex] == 0)
                return error.InvalidDecompositionNutrientCapacity;
            relative_nitrogen = inputs.total_active_biomass_nitrogen_g_n[complex] /
                inputs.maximum_active_biomass_nitrogen_g_n[complex];
            relative_phosphorus = inputs.total_active_biomass_phosphorus_g_p[complex] /
                inputs.maximum_active_biomass_phosphorus_g_p[complex];
            nitrogen_limitation = @min(1, @max(inputs.minimum_nutrient_limitation_fraction, relative_nitrogen));
            phosphorus_limitation = @min(1, @max(inputs.minimum_nutrient_limitation_fraction, relative_phosphorus));
        }
        var activity_concentration: f64 = 0;
        var half_saturation: f64 = 0;
        var effective_substrate: f64 = 0;
        var density_response: f64 = 0;
        var product_response: f64 = 0;
        if (inputs.colonized_soil_organic_carbon_g_c[complex] >
            inputs.negligible_biomass_g_c)
        {
            activity_concentration = if (inputs.biologically_active_water_m3[complex] >
                inputs.negligible_water_m3)
                @min(inputs.maximum_activity_concentration_g_c_m3_h, inputs.microbial_activity_respiration_g_c[complex] /
                    (inputs.biologically_active_water_m3[complex] *
                        inputs.biochemical_time_fraction_h))
            else
                inputs.maximum_activity_concentration_g_c_m3_h;
            half_saturation = if (inputs.surface_layer[complex])
                inputs.surface_zero_activity_half_saturation_g_c_m3 *
                    (1 + activity_concentration /
                        inputs.surface_activity_inhibition_rate_g_c_m3_h)
            else
                inputs.soil_zero_activity_half_saturation_g_c_m3 *
                    (1 + activity_concentration /
                        inputs.soil_activity_inhibition_rate_g_c_m3_h);
            const divisor = if (inputs.soil_layer_mass_Mg[complex] >
                inputs.negligible_biomass_g_c)
                inputs.soil_layer_mass_Mg[complex]
            else
                inputs.soil_layer_volume_m3[complex];
            if (divisor <= 0) return error.InvalidDecompositionLayerGeometry;
            effective_substrate = @sqrt(
                inputs.colonized_soil_organic_carbon_g_c[complex] / divisor,
            );
            density_response = effective_substrate / (effective_substrate + half_saturation);
            product_response = 1 / (1 +
                inputs.dissolved_activity_carbon_concentration_g_c_m3[complex] /
                    inputs.dissolved_carbon_product_inhibition_g_c_m3);
        }
        temporary[complex] = .{
            relative_nitrogen,     relative_phosphorus,    nitrogen_limitation,
            phosphorus_limitation, activity_concentration, half_saturation,
            effective_substrate,   density_response,       product_response,
        };
        for (temporary[complex]) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDecompositionControlResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, complex| @field(state, field.name)[complex] = values[index];
    };
}

fn stateFieldIndex(comptime name: []const u8) usize {
    var index: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (comptime std.mem.eql(u8, field.name, name)) return index;
        index += 1;
    };
    unreachable;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.complex_count;
    if (inputs.surface_layer.len != n) return error.InvalidDecompositionControlDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidDecompositionControlDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDecompositionControlInput;
    };
    inline for (.{
        inputs.surface_zero_activity_half_saturation_g_c_m3,
        inputs.surface_activity_inhibition_rate_g_c_m3_h,
        inputs.soil_zero_activity_half_saturation_g_c_m3,
        inputs.soil_activity_inhibition_rate_g_c_m3_h,
        inputs.dissolved_carbon_product_inhibition_g_c_m3,
        inputs.minimum_nutrient_limitation_fraction,
        inputs.maximum_activity_concentration_g_c_m3_h,
        inputs.negligible_biomass_g_c,
        inputs.negligible_water_m3,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDecompositionControlInput;
    if (inputs.surface_activity_inhibition_rate_g_c_m3_h == 0 or
        inputs.soil_activity_inhibition_rate_g_c_m3_h == 0 or
        inputs.dissolved_carbon_product_inhibition_g_c_m3 == 0 or
        inputs.biochemical_time_fraction_h == 0)
        return error.InvalidDecompositionControlInput;
}

fn fixture() Inputs {
    return .{
        .surface_layer = &.{ true, false },
        .total_active_biomass_carbon_g_c = &.{ 10, 10 },
        .total_active_biomass_nitrogen_g_n = &.{ 2, 1 },
        .total_active_biomass_phosphorus_g_p = &.{ 1, 0.25 },
        .maximum_active_biomass_nitrogen_g_n = &.{ 4, 4 },
        .maximum_active_biomass_phosphorus_g_p = &.{ 2, 2 },
        .colonized_soil_organic_carbon_g_c = &.{ 100, 100 },
        .microbial_activity_respiration_g_c = &.{ 4, 4 },
        .biologically_active_water_m3 = &.{ 2, 2 },
        .soil_layer_mass_Mg = &.{ 10, 10 },
        .soil_layer_volume_m3 = &.{ 20, 20 },
        .dissolved_activity_carbon_concentration_g_c_m3 = &.{ 2, 2 },
        .surface_zero_activity_half_saturation_g_c_m3 = 1,
        .surface_activity_inhibition_rate_g_c_m3_h = 2,
        .soil_zero_activity_half_saturation_g_c_m3 = 2,
        .soil_activity_inhibition_rate_g_c_m3_h = 4,
        .dissolved_carbon_product_inhibition_g_c_m3 = 2,
        .negligible_biomass_g_c = 1e-12,
        .negligible_water_m3 = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "nutrient controls preserve source half-floor and unit cap" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.5, state.nitrogen_limitation_fraction[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.5, state.phosphorus_limitation_fraction[1], 1e-12);
}

test "surface and soil use distinct activity inhibition parameters" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(2, state.decomposition_half_saturation_g_c_m3[0], 1e-12);
    try std.testing.expectApproxEqAbs(3, state.decomposition_half_saturation_g_c_m3[1], 1e-12);
    try std.testing.expect(state.substrate_density_response[0] > 0);
    try std.testing.expectApproxEqAbs(0.5, state.dissolved_carbon_product_response[0], 1e-12);
}

test "dry active substrate uses capped activity concentration" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.biologically_active_water_m3 = &.{ 0, 2 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(100_000, state.microbial_activity_concentration_g_c_m3_h[0]);
}

test "invalid late geometry leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.substrate_density_response[0] = 7;
    var inputs = fixture();
    inputs.soil_layer_mass_Mg = &.{ 0, 10 };
    inputs.soil_layer_volume_m3 = &.{ 0, 20 };
    try std.testing.expectError(error.InvalidDecompositionLayerGeometry, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.substrate_density_response[0]);
}

test "NITRO 3189 threshold equality initializes both limitations to one" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.total_active_biomass_carbon_g_c = &.{ 1e-12, 10 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(1, state.nitrogen_limitation_fraction[0]);
    try std.testing.expectEqual(1, state.phosphorus_limitation_fraction[0]);
}

test "NITRO 3182-3240 derived overflow preserves all decomposition controls" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.nitrogen_limitation_fraction[0] = 7;
    state.dissolved_carbon_product_response[1] = 11;
    var inputs = fixture();
    inputs.total_active_biomass_nitrogen_g_n =
        &.{ std.math.floatMax(f64), 1 };
    inputs.maximum_active_biomass_nitrogen_g_n =
        &.{ std.math.floatMin(f64), 4 };
    try std.testing.expectError(
        error.InvalidDecompositionControlResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.nitrogen_limitation_fraction[0]);
    try std.testing.expectEqual(11, state.dissolved_carbon_product_response[1]);
}
