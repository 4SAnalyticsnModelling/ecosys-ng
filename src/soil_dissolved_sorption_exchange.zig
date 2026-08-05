const std = @import("std");

pub const Inputs = struct {
    surface_layer: []const bool,
    micropore_water_m3: []const f64,
    soil_mass_megagrams: []const f64,
    substrate_complex_fraction: []const f64,
    anion_exchange_capacity_mol_per_megagram: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    dissolved_acetate_g_c: []const f64,
    sorbed_organic_carbon_g_c: []const f64,
    sorbed_organic_nitrogen_g_n: []const f64,
    sorbed_organic_phosphorus_g_p: []const f64,
    sorbed_acetate_g_c: []const f64,
    dissolved_carbon_fraction: []const f64,
    acetate_fraction: []const f64,
    sorption_rate_h: f64,
    sorption_coefficient_megagrams_mol: f64,
    surface_anion_exchange_capacity_mol_per_megagram: f64 = 500,
    negligible_pool_mass: f64,
    negligible_water_m3: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    complex_count: usize,
    effective_anion_exchange_capacity_mol_per_megagram: []f64,
    exchange_volume: []f64,
    aqueous_volume_m3: []f64,
    dissolved_organic_carbon_sorption_g_c: []f64,
    dissolved_organic_nitrogen_sorption_g_n: []f64,
    dissolved_organic_phosphorus_sorption_g_p: []f64,
    acetate_sorption_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize) !State {
        if (complex_count == 0) return error.InvalidSorptionExchangeDimensions;
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

/// Exact NITRO.F 3460--3514 signed dissolved/sorbed organic exchange.
/// Positive flux adsorbs; negative flux desorbs.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([7]f64, state.complex_count);
    defer state.allocator.free(temporary);
    for (0..state.complex_count) |complex| {
        if (inputs.micropore_water_m3[complex] <= inputs.negligible_water_m3 or
            inputs.substrate_complex_fraction[complex] <= 0)
        {
            temporary[complex] = @splat(0);
            continue;
        }
        const exchange_capacity = if (inputs.surface_layer[complex])
            inputs.surface_anion_exchange_capacity_mol_per_megagram
        else
            inputs.anion_exchange_capacity_mol_per_megagram[complex];
        const exchange_volume = inputs.soil_mass_megagrams[complex] * exchange_capacity *
            inputs.sorption_coefficient_megagrams_mol *
            inputs.substrate_complex_fraction[complex];
        const aqueous_volume = inputs.micropore_water_m3[complex] *
            inputs.substrate_complex_fraction[complex];
        const dissolved_carbon = @max(inputs.negligible_pool_mass, inputs.dissolved_organic_carbon_g_c[complex]);
        const dissolved_nitrogen = @max(inputs.negligible_pool_mass, inputs.dissolved_organic_nitrogen_g_n[complex]);
        const dissolved_phosphorus = @max(inputs.negligible_pool_mass, inputs.dissolved_organic_phosphorus_g_p[complex]);
        const dissolved_acetate = @max(inputs.negligible_pool_mass, inputs.dissolved_acetate_g_c[complex]);
        const sorbed_carbon = @max(inputs.negligible_pool_mass, inputs.sorbed_organic_carbon_g_c[complex]);
        const sorbed_nitrogen = @max(inputs.negligible_pool_mass, inputs.sorbed_organic_nitrogen_g_n[complex]);
        const sorbed_phosphorus = @max(inputs.negligible_pool_mass, inputs.sorbed_organic_phosphorus_g_p[complex]);
        const sorbed_acetate = @max(inputs.negligible_pool_mass, inputs.sorbed_acetate_g_c[complex]);
        const carbon_flux = fractionatedFlux(
            inputs.sorption_rate_h,
            inputs.biochemical_time_fraction_h,
            dissolved_carbon,
            sorbed_carbon,
            exchange_volume,
            aqueous_volume,
            inputs.dissolved_carbon_fraction[complex],
        );
        const acetate_flux = fractionatedFlux(
            inputs.sorption_rate_h,
            inputs.biochemical_time_fraction_h,
            dissolved_acetate,
            sorbed_acetate,
            exchange_volume,
            aqueous_volume,
            inputs.acetate_fraction[complex],
        );
        const nitrogen_flux = exchangeFlux(
            inputs.sorption_rate_h,
            inputs.biochemical_time_fraction_h,
            dissolved_nitrogen,
            sorbed_nitrogen,
            exchange_volume,
            aqueous_volume,
        );
        const phosphorus_flux = exchangeFlux(
            inputs.sorption_rate_h,
            inputs.biochemical_time_fraction_h,
            dissolved_phosphorus,
            sorbed_phosphorus,
            exchange_volume,
            aqueous_volume,
        );
        temporary[complex] = .{
            exchange_capacity, exchange_volume, aqueous_volume,
            carbon_flux,       nitrogen_flux,   phosphorus_flux,
            acetate_flux,
        };
        for (temporary[complex]) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteSorptionExchangeResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, complex| @field(state, field.name)[complex] = values[index];
    };
}

fn fractionatedFlux(
    rate: f64,
    time: f64,
    dissolved: f64,
    sorbed: f64,
    exchange_volume: f64,
    aqueous_volume: f64,
    fraction: f64,
) f64 {
    if (fraction > 0) return exchangeFlux(rate, time, dissolved, sorbed, fraction * exchange_volume, fraction * aqueous_volume);
    return exchangeFlux(rate, time, dissolved, sorbed, exchange_volume, aqueous_volume);
}

fn exchangeFlux(
    rate: f64,
    time: f64,
    dissolved: f64,
    sorbed: f64,
    exchange_volume: f64,
    aqueous_volume: f64,
) f64 {
    const denominator = exchange_volume + aqueous_volume;
    if (denominator <= 0) return 0;
    return rate * (dissolved * exchange_volume - sorbed * aqueous_volume) /
        denominator * time;
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
    if (inputs.surface_layer.len != n) return error.InvalidSorptionExchangeDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidSorptionExchangeDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSorptionExchangeInput;
    };
    inline for (.{
        inputs.sorption_rate_h,                            inputs.sorption_coefficient_megagrams_mol,
        inputs.surface_anion_exchange_capacity_mol_per_megagram, inputs.negligible_pool_mass,
        inputs.negligible_water_m3,                        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSorptionExchangeInput;
}

fn fixture() Inputs {
    return .{
        .surface_layer = &.{ true, false },
        .micropore_water_m3 = &.{ 2, 2 },
        .soil_mass_megagrams = &.{ 1, 1 },
        .substrate_complex_fraction = &.{ 0.5, 0.5 },
        .anion_exchange_capacity_mol_per_megagram = &.{ 100, 100 },
        .dissolved_organic_carbon_g_c = &.{ 10, 10 },
        .dissolved_organic_nitrogen_g_n = &.{ 2, 2 },
        .dissolved_organic_phosphorus_g_p = &.{ 1, 1 },
        .dissolved_acetate_g_c = &.{ 4, 4 },
        .sorbed_organic_carbon_g_c = &.{ 1, 1 },
        .sorbed_organic_nitrogen_g_n = &.{ 0.2, 0.2 },
        .sorbed_organic_phosphorus_g_p = &.{ 0.1, 0.1 },
        .sorbed_acetate_g_c = &.{ 0.4, 0.4 },
        .dissolved_carbon_fraction = &.{ 0.8, 0.8 },
        .acetate_fraction = &.{ 0.2, 0.2 },
        .sorption_rate_h = 0.1,
        .sorption_coefficient_megagrams_mol = 0.01,
        .negligible_pool_mass = 1e-12,
        .negligible_water_m3 = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "surface layer uses fixed source anion exchange capacity" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(500, state.effective_anion_exchange_capacity_mol_per_megagram[0]);
    try std.testing.expectEqual(100, state.effective_anion_exchange_capacity_mol_per_megagram[1]);
}

test "signed exchange distinguishes adsorption and desorption" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.sorbed_organic_carbon_g_c = &.{ 1, 100 };
    try calculate(&state, inputs);
    try std.testing.expect(state.dissolved_organic_carbon_sorption_g_c[0] > 0);
    try std.testing.expect(state.dissolved_organic_carbon_sorption_g_c[1] < 0);
}

test "dry complex has zero sorption fluxes" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.micropore_water_m3 = &.{ 0, 2 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.dissolved_organic_carbon_sorption_g_c[0]);
}

test "NITRO wet and active threshold equality zero every exchange output" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.micropore_water_m3 = &.{ 1e-12, 2 };
    inputs.substrate_complex_fraction = &.{ 0.5, 0 };
    try calculate(&state, inputs);
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64)
        for (@field(state, field.name)) |value| try std.testing.expectEqual(0, value);
}

test "DOC and acetate retain distinct positive-fraction source controls" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.dissolved_carbon_fraction = &.{ 0.8, 0.8 };
    inputs.acetate_fraction = &.{ 0.2, 0.2 };
    try calculate(&state, inputs);
    try std.testing.expect(
        state.dissolved_organic_carbon_sorption_g_c[0] !=
            state.acetate_sorption_g_c[0],
    );
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_sorption_g_c[0] = 7;
    var inputs = fixture();
    inputs.soil_mass_megagrams = &.{ std.math.nan(f64), 1 };
    try std.testing.expectError(error.InvalidSorptionExchangeInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_sorption_g_c[0]);
}

test "NITRO 3460-3514 derived overflow preserves all sorption outputs" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_sorption_g_c[0] = 7;
    state.acetate_sorption_g_c[1] = 11;
    var inputs = fixture();
    inputs.soil_mass_megagrams = &.{ std.math.floatMax(f64), 1 };
    inputs.sorption_coefficient_megagrams_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSorptionExchangeResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_sorption_g_c[0]);
    try std.testing.expectEqual(11, state.acetate_sorption_g_c[1]);
}
