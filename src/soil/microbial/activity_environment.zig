const std = @import("std");

pub const Inputs = struct {
    substrate_complex_count: usize,
    population_count: usize,
    enabled: []const bool,
    fungal_water_response: []const bool,
    use_total_colonized_carbon: []const bool,
    active_biomass_g_c: []const f64,
    active_biomass_by_complex_g_c: []const f64,
    colonized_structural_carbon_by_complex_g_c: []const f64,
    total_active_biomass_g_c: f64,
    density_reference_biomass_g_c: f64,
    total_colonized_structural_carbon_g_c: f64,
    growth_temperature_response: f64,
    maintenance_temperature_response: f64,
    matric_plus_osmotic_potential_megapascal: f64,
    fungal_water_sensitivity_per_megapascal: f64,
    other_microbe_water_sensitivity_per_megapascal: f64,
    decomposition_density_inhibition_ratio: f64,
    maintenance_density_inhibition_ratio: f64,
    negligible_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    water_response: []f64,
    temperature_water_response: []f64,
    maintenance_temperature_response: []f64,
    fraction_of_total_active_biomass: []f64,
    fraction_of_density_reference_biomass: []f64,
    fraction_of_complex_active_biomass: []f64,
    active_to_colonized_carbon_ratio: []f64,
    decomposition_density_response: []f64,
    maintenance_density_response: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidMicrobialActivityDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(state, field.name) = try allocator.alloc(f64, unit_count);
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
};

/// NITRO.F 579--664. This consumes TFNX/TFNY from the layer environment,
/// OMA/TOMA/TOMN/TOMK from the biomass environment, and OSAT/TOSA from the
/// substrate environment. Runtime masks retain the source fungal and
/// autotrophic-complex branches without fixed population or complex counts.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    calculateValidated(&staged, inputs);
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

fn calculateValidated(state: *State, inputs: Inputs) void {
    clear(state);
    for (0..inputs.substrate_complex_count) |complex| {
        for (0..inputs.population_count) |population| {
            const unit = complex * inputs.population_count + population;
            if (!inputs.enabled[unit]) continue;
            const water_coefficient_per_megapascal: f64 =
                if (inputs.fungal_water_response[unit])
                    inputs.fungal_water_sensitivity_per_megapascal
                else
                    inputs.other_microbe_water_sensitivity_per_megapascal;
            const water_response = @exp(
                water_coefficient_per_megapascal *
                    inputs.matric_plus_osmotic_potential_megapascal,
            );
            state.water_response[unit] = water_response;
            state.temperature_water_response[unit] =
                inputs.growth_temperature_response * water_response;
            state.maintenance_temperature_response[unit] =
                inputs.maintenance_temperature_response;

            const active = inputs.active_biomass_g_c[unit];
            if (active <= 0) continue;
            state.fraction_of_total_active_biomass[unit] =
                if (inputs.total_active_biomass_g_c >
                inputs.negligible_carbon_g_c)
                    active / inputs.total_active_biomass_g_c
                else
                    1;
            state.fraction_of_density_reference_biomass[unit] =
                if (inputs.density_reference_biomass_g_c >
                inputs.negligible_carbon_g_c)
                    active / inputs.density_reference_biomass_g_c
                else
                    1;
            state.fraction_of_complex_active_biomass[unit] =
                if (inputs.active_biomass_by_complex_g_c[complex] >
                inputs.negligible_carbon_g_c)
                    active / inputs.active_biomass_by_complex_g_c[complex]
                else
                    1;

            const colonized_carbon =
                if (inputs.use_total_colonized_carbon[unit])
                    inputs.total_colonized_structural_carbon_g_c
                else
                    inputs.colonized_structural_carbon_by_complex_g_c[complex];
            if (colonized_carbon > inputs.negligible_carbon_g_c) {
                const ratio = active / colonized_carbon;
                state.active_to_colonized_carbon_ratio[unit] = ratio;
                state.decomposition_density_response[unit] =
                    ratio /
                    (ratio + inputs.decomposition_density_inhibition_ratio);
                state.maintenance_density_response[unit] =
                    ratio /
                    (ratio + inputs.maintenance_density_inhibition_ratio);
            } else {
                state.decomposition_density_response[unit] = 1;
                state.maintenance_density_response[unit] = 1;
            }
        }
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteMicrobialActivityEnvironment;
        };
}

/// Exact source WFNG branch: only population N=3 in complexes K<=4 uses
/// the fungal 0.05 MPa^-1 coefficient.
pub fn sourceFungalWaterResponse(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return zero_based_complex <= 4 and zero_based_population == 2;
}

/// Exact source density denominator branch: K>4 uses layer-total TOSA.
pub fn sourceUsesTotalColonizedCarbon(zero_based_complex: usize) bool {
    return zero_based_complex > 4;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.substrate_complex_count == 0 or inputs.population_count == 0)
        return error.InvalidMicrobialActivityDimensions;
    const unit_count = try std.math.mul(
        usize,
        inputs.substrate_complex_count,
        inputs.population_count,
    );
    if (state.unit_count != unit_count or
        inputs.enabled.len != unit_count or
        inputs.fungal_water_response.len != unit_count or
        inputs.use_total_colonized_carbon.len != unit_count or
        inputs.active_biomass_g_c.len != unit_count or
        inputs.active_biomass_by_complex_g_c.len !=
            inputs.substrate_complex_count or
        inputs.colonized_structural_carbon_by_complex_g_c.len !=
            inputs.substrate_complex_count)
        return error.InvalidMicrobialActivityDimensions;
    inline for (.{
        inputs.active_biomass_g_c,
        inputs.active_biomass_by_complex_g_c,
        inputs.colonized_structural_carbon_by_complex_g_c,
    }) |values| for (values) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidMicrobialActivityInput;
    };
    inline for (.{
        inputs.total_active_biomass_g_c,
        inputs.density_reference_biomass_g_c,
        inputs.total_colonized_structural_carbon_g_c,
        inputs.growth_temperature_response,
        inputs.maintenance_temperature_response,
        inputs.fungal_water_sensitivity_per_megapascal,
        inputs.other_microbe_water_sensitivity_per_megapascal,
        inputs.decomposition_density_inhibition_ratio,
        inputs.maintenance_density_inhibition_ratio,
        inputs.negligible_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMicrobialActivityInput;
    if (!std.math.isFinite(inputs.matric_plus_osmotic_potential_megapascal) or
        inputs.fungal_water_sensitivity_per_megapascal <= 0 or
        inputs.other_microbe_water_sensitivity_per_megapascal <= 0 or
        inputs.decomposition_density_inhibition_ratio <= 0 or
        inputs.maintenance_density_inhibition_ratio <= 0)
        return error.InvalidMicrobialActivityInput;
    for (0..inputs.substrate_complex_count) |complex| {
        var summed_active_g_c: f64 = 0;
        for (0..inputs.population_count) |population| {
            const unit = complex * inputs.population_count + population;
            if (inputs.enabled[unit])
                summed_active_g_c += inputs.active_biomass_g_c[unit];
        }
        const tolerance = 64 * std.math.floatEps(f64) *
            @max(1, inputs.active_biomass_by_complex_g_c[complex]);
        if (@abs(
            summed_active_g_c -
                inputs.active_biomass_by_complex_g_c[complex],
        ) > tolerance) return error.InconsistentComplexActiveBiomass;
    }
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source fungal and density branches reproduce NITRO selectors" {
    try std.testing.expect(sourceFungalWaterResponse(0, 2));
    try std.testing.expect(sourceFungalWaterResponse(4, 2));
    try std.testing.expect(!sourceFungalWaterResponse(5, 2));
    try std.testing.expect(!sourceFungalWaterResponse(2, 1));
    try std.testing.expect(!sourceUsesTotalColonizedCarbon(4));
    try std.testing.expect(sourceUsesTotalColonizedCarbon(5));
}

test "temperature water allocation and density responses reproduce NITRO equations" {
    var state = try State.init(std.testing.allocator, 4);
    defer state.deinit();
    try calculate(&state, .{
        .substrate_complex_count = 2,
        .population_count = 2,
        .enabled = &.{ true, true, true, false },
        .fungal_water_response = &.{ false, true, false, false },
        .use_total_colonized_carbon = &.{ false, false, true, false },
        .active_biomass_g_c = &.{ 2, 3, 5, 0 },
        .active_biomass_by_complex_g_c = &.{ 5, 5 },
        .colonized_structural_carbon_by_complex_g_c = &.{ 20, 40 },
        .total_active_biomass_g_c = 10,
        .density_reference_biomass_g_c = 4,
        .total_colonized_structural_carbon_g_c = 100,
        .growth_temperature_response = 2,
        .maintenance_temperature_response = 3,
        .matric_plus_osmotic_potential_megapascal = -2,
        .fungal_water_sensitivity_per_megapascal = 0.05,
        .other_microbe_water_sensitivity_per_megapascal = 0.10,
        .decomposition_density_inhibition_ratio = 0.2,
        .maintenance_density_inhibition_ratio = 0.1,
        .negligible_carbon_g_c = 1e-12,
    });
    try std.testing.expectApproxEqAbs(@exp(-0.2), state.water_response[0], 1e-15);
    try std.testing.expectApproxEqAbs(@exp(-0.1), state.water_response[1], 1e-15);
    try std.testing.expectApproxEqAbs(2 * @exp(-0.2), state.temperature_water_response[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 3), state.maintenance_temperature_response[2]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.fraction_of_total_active_biomass[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.fraction_of_complex_active_biomass[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), state.active_to_colonized_carbon_ratio[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0 / 3.0), state.decomposition_density_response[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), state.maintenance_density_response[0], 1e-15);
    // Complex 1 is marked K>4-equivalent and therefore uses total TOSA.
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.active_to_colonized_carbon_ratio[2], 1e-15);
}

test "invalid final complex total leaves activity state unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.water_response[0] = 9;
    try std.testing.expectError(error.InconsistentComplexActiveBiomass, calculate(&state, .{
        .substrate_complex_count = 1,
        .population_count = 2,
        .enabled = &.{ true, true },
        .fungal_water_response = &.{ false, false },
        .use_total_colonized_carbon = &.{ false, false },
        .active_biomass_g_c = &.{ 1, 2 },
        .active_biomass_by_complex_g_c = &.{4},
        .colonized_structural_carbon_by_complex_g_c = &.{10},
        .total_active_biomass_g_c = 3,
        .density_reference_biomass_g_c = 1,
        .total_colonized_structural_carbon_g_c = 10,
        .growth_temperature_response = 1,
        .maintenance_temperature_response = 1,
        .matric_plus_osmotic_potential_megapascal = 0,
        .fungal_water_sensitivity_per_megapascal = 0.05,
        .other_microbe_water_sensitivity_per_megapascal = 0.10,
        .decomposition_density_inhibition_ratio = 0.1,
        .maintenance_density_inhibition_ratio = 0.1,
        .negligible_carbon_g_c = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 9), state.water_response[0]);
}

test "NITRO 579-664 derived overflow leaves activity state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.water_response[0] = 9;
    try std.testing.expectError(
        error.NonFiniteMicrobialActivityEnvironment,
        calculate(&state, .{
            .substrate_complex_count = 1,
            .population_count = 1,
            .enabled = &.{true},
            .fungal_water_response = &.{false},
            .use_total_colonized_carbon = &.{false},
            .active_biomass_g_c = &.{1},
            .active_biomass_by_complex_g_c = &.{1},
            .colonized_structural_carbon_by_complex_g_c = &.{1},
            .total_active_biomass_g_c = 1,
            .density_reference_biomass_g_c = 1,
            .total_colonized_structural_carbon_g_c = 1,
            .growth_temperature_response = 1,
            .maintenance_temperature_response = 1,
            .matric_plus_osmotic_potential_megapascal = 1.0e3,
            .fungal_water_sensitivity_per_megapascal = 1,
            .other_microbe_water_sensitivity_per_megapascal = 1,
            .decomposition_density_inhibition_ratio = 0.1,
            .maintenance_density_inhibition_ratio = 0.1,
            .negligible_carbon_g_c = 1e-12,
        }),
    );
    try std.testing.expectEqual(@as(f64, 9), state.water_response[0]);
}
