const std = @import("std");
const compute = @import("../core/compute.zig");
const organic = @import("../soil/organic/initialization.zig");
const priming = @import("../soil/biogeochemistry/organic_priming_exchange.zig");
const respiration = @import("microbial_respiration_step.zig");
const environment = @import("microbial_environment_step.zig");
const uptake = @import("microbial_substrate_uptake_step.zig");

pub const substrate_count: usize = respiration.litter_complex_count;

pub const State = struct {
    exchange: priming.State,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        return .{ .exchange = try priming.State.init(allocator, cell_count, substrate_count, respiration.source_population_count, organic.kinetic_fraction_count) };
    }

    pub fn deinit(self: *State) void {
        self.exchange.deinit();
        self.* = undefined;
    }

    pub fn cellCount(self: State) usize {
        return self.exchange.cell_count;
    }
};

pub const ApplyContext = struct {
    result: *State,
    surface_organic: *const organic.State,
    respiration: *const respiration.State,
    substrate_uptake: *const uptake.State,
    environment: *const environment.State,
    matric_plus_osmotic_potential_megapascal: []const f64,
    dissolved_priming_rate_per_h: f64,
    microbial_priming_rate_per_h: f64,
    timestep_h: f64,
    negligible_carbon_g_c: f64,
};

/// Surface specialization of NITRO's all-K priming redistribution.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    for (range.first..range.end) |cell| {
        var substrate_carbon_g_c: [substrate_count]f64 = undefined;
        var activity_g_c: [substrate_count]f64 = .{0} ** substrate_count;
        var microbial_response: [substrate_count * respiration.source_population_count]f64 = undefined;
        var dissolved_available_after_uptake: [substrate_count]organic.ElementPool = undefined;
        var acetate_available_after_uptake_g_c: [substrate_count]f64 = undefined;
        for (0..substrate_count) |complex| {
            var carbon_g_c: f64 = 0;
            const structural_base = (cell * organic.substrate_count + complex) * organic.structural_fraction_count;
            for (0..organic.structural_fraction_count) |fraction| carbon_g_c += context.surface_organic.structural[structural_base + fraction].carbon_g_c;
            const residue_base = (cell * organic.substrate_count + complex) * organic.residue_fraction_count;
            for (0..organic.residue_fraction_count) |fraction| carbon_g_c += context.surface_organic.residue[residue_base + fraction].carbon_g_c;
            const mobile = cell * organic.substrate_count + complex;
            dissolved_available_after_uptake[complex] = context.surface_organic.dissolved[mobile];
            acetate_available_after_uptake_g_c[complex] = context.surface_organic.dissolved_acetate_carbon_g_c[mobile];
            carbon_g_c += context.surface_organic.adsorbed[mobile].carbon_g_c + context.surface_organic.adsorbed_acetate_carbon_g_c[mobile];
            substrate_carbon_g_c[complex] = carbon_g_c;
            for (0..respiration.source_population_count) |population| {
                const unit = cell * respiration.unit_count_per_cell + complex * respiration.source_population_count + population;
                dissolved_available_after_uptake[complex].carbon_g_c -= context.substrate_uptake.doc_uptake_g_c[unit];
                dissolved_available_after_uptake[complex].nitrogen_g_n -= context.substrate_uptake.dissolved_organic_nitrogen_uptake_g_n[unit];
                dissolved_available_after_uptake[complex].phosphorus_g_p -= context.substrate_uptake.dissolved_organic_phosphorus_uptake_g_p[unit];
                acetate_available_after_uptake_g_c[complex] -= context.substrate_uptake.acetate_uptake_g_c[unit];
                activity_g_c[complex] += context.respiration.unlimited_respiration_g_c[unit];
                microbial_response[complex * respiration.source_population_count + population] = context.environment.growth_temperature_response[cell] * @exp((if (population == 2) @as(f64, 0.05) else 0.10) * context.matric_plus_osmotic_potential_megapascal[cell]);
            }
            inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| if (@field(dissolved_available_after_uptake[complex], field.name) < -context.negligible_carbon_g_c) return error.SurfaceOrganicPrimingUptakeOverdraw;
            dissolved_available_after_uptake[complex].carbon_g_c = @max(0, dissolved_available_after_uptake[complex].carbon_g_c);
            dissolved_available_after_uptake[complex].nitrogen_g_n = @max(0, dissolved_available_after_uptake[complex].nitrogen_g_n);
            dissolved_available_after_uptake[complex].phosphorus_g_p = @max(0, dissolved_available_after_uptake[complex].phosphorus_g_p);
            if (acetate_available_after_uptake_g_c[complex] < -context.negligible_carbon_g_c) return error.SurfaceOrganicPrimingUptakeOverdraw;
            acetate_available_after_uptake_g_c[complex] = @max(0, acetate_available_after_uptake_g_c[complex]);
        }
        const microbial_per_surface_cell = organic.microbial_substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count;
        const participating_microbial_count = substrate_count * respiration.source_population_count * organic.kinetic_fraction_count;
        try priming.deriveCell(&context.result.exchange, cell, .{
            .substrate_carbon_g_c = &substrate_carbon_g_c,
            .microbial_activity_g_c_per_step = &activity_g_c,
            .dissolved = context.surface_organic.dissolved[cell * organic.substrate_count ..][0..substrate_count],
            .dissolved_available_after_uptake = &dissolved_available_after_uptake,
            .dissolved_acetate_carbon_g_c = context.surface_organic.dissolved_acetate_carbon_g_c[cell * organic.substrate_count ..][0..substrate_count],
            .dissolved_acetate_available_after_uptake_g_c = &acetate_available_after_uptake_g_c,
            .microbial = context.surface_organic.microbial[cell * microbial_per_surface_cell ..][0..participating_microbial_count],
            .microbial_temperature_water_response = &microbial_response,
            .dissolved_priming_rate_per_h = context.dissolved_priming_rate_per_h,
            .microbial_priming_rate_per_h = context.microbial_priming_rate_per_h,
            .decomposition_temperature_response = context.environment.aqueous_diffusion_temperature_response[cell],
            .timestep_h = context.timestep_h,
            .negligible_carbon_g_c = context.negligible_carbon_g_c,
        });
    }
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const cells = context.result.cellCount();
    if (range.first > range.end or range.end > cells or context.surface_organic.layer_count != cells or context.respiration.cell_count != cells or context.substrate_uptake.cell_count != cells or context.environment.growth_temperature_response.len != cells or context.environment.aqueous_diffusion_temperature_response.len != cells or context.matric_plus_osmotic_potential_megapascal.len != cells) return error.SurfaceOrganicPrimingDimensionMismatch;
    inline for (.{ context.dissolved_priming_rate_per_h, context.microbial_priming_rate_per_h, context.timestep_h, context.negligible_carbon_g_c }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfaceOrganicPrimingParameter;
    if (context.timestep_h <= 0) return error.InvalidSurfaceOrganicPrimingParameter;
}

test "surface priming derives conservative independent cell deltas" {
    var organic_state = try organic.State.init(std.testing.allocator, 2);
    defer organic_state.deinit();
    organic_state.structural[0].carbon_g_c = 10;
    organic_state.structural[organic.structural_fraction_count].carbon_g_c = 20;
    organic_state.dissolved[0].carbon_g_c = 4;
    organic_state.dissolved[1].carbon_g_c = 1;
    var respiration_state = try respiration.State.init(std.testing.allocator, 2);
    defer respiration_state.deinit();
    respiration_state.unlimited_respiration_g_c[0] = 2;
    respiration_state.unlimited_respiration_g_c[respiration.source_population_count] = 0.5;
    var environment_state = try environment.State.init(std.testing.allocator, 2);
    defer environment_state.deinit();
    var uptake_state = try uptake.State.init(std.testing.allocator, 2);
    defer uptake_state.deinit();
    environment_state.growth_temperature_response[0] = 1;
    environment_state.aqueous_diffusion_temperature_response[0] = 1;
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var context: ApplyContext = .{ .result = &state, .surface_organic = &organic_state, .respiration = &respiration_state, .substrate_uptake = &uptake_state, .environment = &environment_state, .matric_plus_osmotic_potential_megapascal = &.{ 0, 0 }, .dissolved_priming_rate_per_h = 0.01, .microbial_priming_rate_per_h = 0.001, .timestep_h = 1, .negligible_carbon_g_c = 1e-12 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    var dissolved_carbon_change: f64 = 0;
    for (state.exchange.dissolved_change[0..substrate_count]) |pool| dissolved_carbon_change += pool.carbon_g_c;
    try std.testing.expectApproxEqAbs(@as(f64, 0), dissolved_carbon_change, 1e-15);
    try std.testing.expect(state.exchange.dissolved_change[0].carbon_g_c != 0);
    for (state.exchange.dissolved_change[substrate_count..]) |pool| try std.testing.expectEqual(@as(f64, 0), pool.carbon_g_c);
}
