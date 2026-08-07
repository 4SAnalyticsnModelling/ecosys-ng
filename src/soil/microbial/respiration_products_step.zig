const std = @import("std");
const compute = @import("../../core/compute.zig");
const microbial = @import("state.zig");
const fluxes = @import("../nutrients/nitrogen_flux_workspace.zig");
const respiration_activity = @import("respiration_activity.zig");

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    process_unit_count_per_layer: usize,
    carbon_dioxide_g_c: []f64,
    acetate_g_c: []f64,
    methane_g_c: []f64,
    hydrogen_g_h: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize, process_unit_count_per_layer: usize) !State {
        if (layer_count == 0 or process_unit_count_per_layer == 0) return error.InvalidSoilRespirationProductDimensions;
        const count = try std.math.mul(usize, layer_count, process_unit_count_per_layer);
        const co2 = try allocator.alloc(f64, count);
        errdefer allocator.free(co2);
        const acetate = try allocator.alloc(f64, count);
        errdefer allocator.free(acetate);
        const methane = try allocator.alloc(f64, count);
        errdefer allocator.free(methane);
        const hydrogen = try allocator.alloc(f64, count);
        @memset(co2, 0);
        @memset(acetate, 0);
        @memset(methane, 0);
        @memset(hydrogen, 0);
        return .{ .allocator = allocator, .layer_count = layer_count, .process_unit_count_per_layer = process_unit_count_per_layer, .carbon_dioxide_g_c = co2, .acetate_g_c = acetate, .methane_g_c = methane, .hydrogen_g_h = hydrogen };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.hydrogen_g_h);
        self.allocator.free(self.methane_g_c);
        self.allocator.free(self.acetate_g_c);
        self.allocator.free(self.carbon_dioxide_g_c);
        self.* = undefined;
    }
};

pub const ApplyContext = struct { result: *State, microbial_state: *const microbial.State, respiration_fluxes: *const fluxes.State };

/// NITRO RCO2X/RCH3X/RCH4X/RH2GX for K=0..4. Source populations N=4/7
/// ferment; N=5 is acetotrophic methanogenic; the others respire to CO2.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    const populations = context.microbial_state.population_count;
    for (range.first..range.end) |layer| for (0..context.microbial_state.substrate_count) |substrate| for (0..populations) |population| {
        const unit = layer * context.result.process_unit_count_per_layer + substrate * populations + population;
        const respiration_g_c = context.respiration_fluxes.actual_aerobic_respiration_g_c[unit];
        var co2_g_c: f64 = respiration_g_c;
        var acetate_g_c: f64 = 0;
        var methane_g_c: f64 = 0;
        var hydrogen_g_h: f64 = 0;
        switch (respiration_activity.sourceMetabolism(population)) {
            .fermenting_heterotroph => {
                co2_g_c = 0.333 * respiration_g_c;
                acetate_g_c = 0.667 * respiration_g_c;
                hydrogen_g_h = 0.111 * respiration_g_c;
            },
            .acetotrophic_methanogen => {
                co2_g_c = 0.5 * respiration_g_c;
                methane_g_c = 0.5 * respiration_g_c;
            },
            .aerobic_heterotroph => {},
        }
        context.result.carbon_dioxide_g_c[unit] = co2_g_c;
        context.result.acetate_g_c[unit] = acetate_g_c;
        context.result.methane_g_c[unit] = methane_g_c;
        context.result.hydrogen_g_h[unit] = hydrogen_g_h;
        inline for (.{ co2_g_c, acetate_g_c, methane_g_c, hydrogen_g_h }) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteSoilRespirationProduct;
        if (@abs(respiration_g_c - co2_g_c - acetate_g_c - methane_g_c) > 1e-12 * @max(1, respiration_g_c)) return error.SoilRespirationCarbonBalanceFailure;
    };
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const layers = context.result.layer_count;
    if (range.first > range.end or range.end > layers or context.microbial_state.cell_count * context.microbial_state.layer_count != layers or context.respiration_fluxes.layer_count != layers or context.result.process_unit_count_per_layer != context.microbial_state.substrate_count * context.microbial_state.population_count or context.respiration_fluxes.process_unit_count_per_layer != context.result.process_unit_count_per_layer) return error.InvalidSoilRespirationProductDimensions;
}

test "legacy soil populations route respiration to CO2 acetate methane and hydrogen" {
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 1, 7);
    defer microbial_state.deinit();
    var respiration = try fluxes.State.init(std.testing.allocator, 1, 7);
    defer respiration.deinit();
    @memset(respiration.actual_aerobic_respiration_g_c, 1);
    var state = try State.init(std.testing.allocator, 1, 7);
    defer state.deinit();
    var context: ApplyContext = .{ .result = &state, .microbial_state = &microbial_state, .respiration_fluxes = &respiration };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expectEqual(@as(f64, 1), state.carbon_dioxide_g_c[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.667), state.acetate_g_c[3], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.111), state.hydrogen_g_h[6], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), state.methane_g_c[4], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.carbon_dioxide_g_c[4] + state.methane_g_c[4], 1e-15);
}
