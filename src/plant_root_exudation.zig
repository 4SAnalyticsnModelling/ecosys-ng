const std = @import("std");
const RootState = @import("plant_root_system.zig").State;
const SoilState = @import("soil_organic_initialization.zig").State;

pub const substrate_count: usize = @import("plant_root_system.zig").organic_substrate_count;

pub const Parameters = struct {
    maximum_root_carbon_concentration_g_c_per_m3: f64,
    root_mobile_nitrogen_exchange_fraction: f64,
    root_mobile_phosphorus_exchange_fraction: f64,
    exchange_rate_per_h: f64,

    pub fn validate(self: Parameters) !void {
        inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(self, field.name))) return error.NonFiniteRootExudationParameter;
        if (self.maximum_root_carbon_concentration_g_c_per_m3 <= 0 or
            self.root_mobile_nitrogen_exchange_fraction < 0 or self.root_mobile_nitrogen_exchange_fraction > 1 or
            self.root_mobile_phosphorus_exchange_fraction < 0 or self.root_mobile_phosphorus_exchange_fraction > 1 or
            self.exchange_rate_per_h < 0)
            return error.InvalidRootExudationParameter;
    }
};

pub fn compatibilityParameters() Parameters {
    return .{
        .maximum_root_carbon_concentration_g_c_per_m3 = 1.0e3,
        .root_mobile_nitrogen_exchange_fraction = 0.1,
        .root_mobile_phosphorus_exchange_fraction = 0.1,
        .exchange_rate_per_h = 1.0e-3,
    };
}

pub const Input = struct {
    soil_water_volume_m3: f64,
    root_water_volume_m3: f64,
    soil_dissolved_carbon_g_c: f64,
    soil_dissolved_nitrogen_g_n: f64,
    soil_dissolved_phosphorus_g_p: f64,
    root_nonstructural_carbon_g_c: f64,
    root_nonstructural_nitrogen_g_n: f64,
    root_nonstructural_phosphorus_g_p: f64,
    maximum_root_carbon_concentration_g_c_per_m3: f64,
    root_mobile_nitrogen_exchange_fraction: f64,
    root_mobile_phosphorus_exchange_fraction: f64,
    exchange_rate_per_h: f64,
    timestep_h: f64,
    significance_threshold_g: f64,
};

/// Positive values are soil uptake by the root; negative values are exudation,
/// retaining the RDFOMC/RDFOMN/RDFOMP source sign convention.
pub const Result = struct { carbon_g_c: f64, nitrogen_g_n: f64, phosphorus_g_p: f64 };

pub const Competitor = struct { plant: usize, domain: usize, layer: usize };

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    competitor_capacity: usize,
    competitors: []Competitor,
    staged_results: []Result,

    pub fn init(allocator: std.mem.Allocator, competitor_capacity: usize) !Workspace {
        if (competitor_capacity == 0) return error.ZeroRootExudationCompetitorCapacity;
        const competitors = try allocator.alloc(Competitor, competitor_capacity);
        errdefer allocator.free(competitors);
        const staged = try allocator.alloc(Result, try std.math.mul(usize, competitor_capacity, substrate_count));
        @memset(staged, .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 });
        return .{ .allocator = allocator, .competitor_capacity = competitor_capacity, .competitors = competitors, .staged_results = staged };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.staged_results);
        self.allocator.free(self.competitors);
        self.* = undefined;
    }

    pub fn stageLayer(
        self: *Workspace,
        roots: *const RootState,
        soil: *const SoilState,
        soil_layer: usize,
        biologically_active_water_m3: f64,
        substrate_fraction: []const f64,
        parameters: Parameters,
        significance_threshold_g: f64,
        competitor_count: usize,
    ) !void {
        if (competitor_count > self.competitor_capacity) return error.RootExudationCompetitorOutOfBounds;
        if (substrate_fraction.len != substrate_count or !std.math.isFinite(biologically_active_water_m3) or biologically_active_water_m3 < 0) return error.InvalidRootExudationWorkspaceInput;
        try parameters.validate();
        for (self.competitors[0..competitor_count], 0..) |competitor, competitor_index| {
            const root = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
            for (0..substrate_count) |substrate| {
                const fraction = substrate_fraction[substrate];
                if (!std.math.isFinite(fraction) or fraction < 0) return error.InvalidRootExudationWorkspaceInput;
                const output = &self.staged_results[competitor_index * substrate_count + substrate];
                const soil_index = soil_layer * substrate_count + substrate;
                const dissolved = soil.dissolved[soil_index];
                const soil_water = biologically_active_water_m3 * fraction;
                output.* = if (soil_water > significance_threshold_g and roots.aqueous_volume_m3[root] > significance_threshold_g)
                    try calculate(.{
                        .soil_water_volume_m3 = soil_water,
                        .root_water_volume_m3 = roots.aqueous_volume_m3[root],
                        .soil_dissolved_carbon_g_c = dissolved.carbon_g_c,
                        .soil_dissolved_nitrogen_g_n = dissolved.nitrogen_g_n,
                        .soil_dissolved_phosphorus_g_p = dissolved.phosphorus_g_p,
                        .root_nonstructural_carbon_g_c = roots.mobile_carbon_g[root],
                        .root_nonstructural_nitrogen_g_n = roots.mobile_nitrogen_g[root],
                        .root_nonstructural_phosphorus_g_p = roots.mobile_phosphorus_g[root],
                        .maximum_root_carbon_concentration_g_c_per_m3 = parameters.maximum_root_carbon_concentration_g_c_per_m3,
                        .root_mobile_nitrogen_exchange_fraction = parameters.root_mobile_nitrogen_exchange_fraction,
                        .root_mobile_phosphorus_exchange_fraction = parameters.root_mobile_phosphorus_exchange_fraction,
                        .exchange_rate_per_h = parameters.exchange_rate_per_h,
                        .timestep_h = 1,
                        .significance_threshold_g = significance_threshold_g,
                    })
                else
                    .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
            }
        }
    }

    pub fn commitLayer(self: *Workspace, roots: *RootState, soil: *SoilState, soil_layer: usize, competitor_count: usize) !void {
        if (competitor_count > self.competitor_capacity) return error.RootExudationCompetitorOutOfBounds;
        // Every root/domain/layer competitor is unique. Validate all soil and
        // root destinations before publishing any of them.
        for (0..competitor_count) |first| for (first + 1..competitor_count) |second| {
            const a = self.competitors[first];
            const b = self.competitors[second];
            if (a.plant == b.plant and a.domain == b.domain and a.layer == b.layer) return error.DuplicateRootExudationCompetitor;
        };
        for (0..substrate_count) |substrate| {
            var total: Result = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
            for (0..competitor_count) |competitor_index| {
                const result = self.staged_results[competitor_index * substrate_count + substrate];
                total.carbon_g_c += result.carbon_g_c;
                total.nitrogen_g_n += result.nitrogen_g_n;
                total.phosphorus_g_p += result.phosphorus_g_p;
            }
            const soil_index = soil_layer * substrate_count + substrate;
            const dissolved = soil.dissolved[soil_index];
            inline for (.{
                dissolved.carbon_g_c - total.carbon_g_c,
                dissolved.nitrogen_g_n - total.nitrogen_g_n,
                dissolved.phosphorus_g_p - total.phosphorus_g_p,
            }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootExudationCommit;
            inline for (.{
                dissolved.carbon_g_c - total.carbon_g_c,
                dissolved.nitrogen_g_n - total.nitrogen_g_n,
                dissolved.phosphorus_g_p - total.phosphorus_g_p,
            }) |value| if (value < -1.0e-12) return error.InsufficientElementForRootExudation;
        }
        for (self.competitors[0..competitor_count], 0..) |competitor, competitor_index| {
            const root = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
            var root_delta: Result = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
            for (0..substrate_count) |substrate| {
                const result = self.staged_results[competitor_index * substrate_count + substrate];
                root_delta.carbon_g_c += result.carbon_g_c;
                root_delta.nitrogen_g_n += result.nitrogen_g_n;
                root_delta.phosphorus_g_p += result.phosphorus_g_p;
            }
            inline for (.{
                roots.mobile_carbon_g[root] + root_delta.carbon_g_c,
                roots.mobile_nitrogen_g[root] + root_delta.nitrogen_g_n,
                roots.mobile_phosphorus_g[root] + root_delta.phosphorus_g_p,
            }) |value| if (!std.math.isFinite(value) or value < -1.0e-12) return error.InsufficientElementForRootExudation;
        }
        for (self.competitors[0..competitor_count], 0..) |competitor, competitor_index| {
            const root = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
            for (0..substrate_count) |substrate| {
                const result = self.staged_results[competitor_index * substrate_count + substrate];
                const soil_index = soil_layer * substrate_count + substrate;
                const root_substrate = try roots.substrateIndex(competitor.plant, competitor.domain, competitor.layer, substrate);
                soil.dissolved[soil_index].carbon_g_c = @max(0, soil.dissolved[soil_index].carbon_g_c - result.carbon_g_c);
                soil.dissolved[soil_index].nitrogen_g_n = @max(0, soil.dissolved[soil_index].nitrogen_g_n - result.nitrogen_g_n);
                soil.dissolved[soil_index].phosphorus_g_p = @max(0, soil.dissolved[soil_index].phosphorus_g_p - result.phosphorus_g_p);
                roots.mobile_carbon_g[root] = @max(0, roots.mobile_carbon_g[root] + result.carbon_g_c);
                roots.mobile_nitrogen_g[root] = @max(0, roots.mobile_nitrogen_g[root] + result.nitrogen_g_n);
                roots.mobile_phosphorus_g[root] = @max(0, roots.mobile_phosphorus_g[root] + result.phosphorus_g_p);
                roots.exudate_carbon_exchange_g_c_per_h[root_substrate] += result.carbon_g_c;
                roots.exudate_nitrogen_exchange_g_n_per_h[root_substrate] += result.nitrogen_g_n;
                roots.exudate_phosphorus_exchange_g_p_per_h[root_substrate] += result.phosphorus_g_p;
            }
        }
    }
};

pub const GridWorkspace = struct {
    allocator: std.mem.Allocator,
    per_cell: []Workspace,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, competitor_capacity_per_cell: usize) !GridWorkspace {
        if (cell_count == 0) return error.ZeroRootExudationGridCells;
        const cells = try allocator.alloc(Workspace, cell_count);
        errdefer allocator.free(cells);
        var initialized: usize = 0;
        errdefer for (cells[0..initialized]) |*cell| cell.deinit();
        for (cells) |*cell| {
            cell.* = try Workspace.init(allocator, competitor_capacity_per_cell);
            initialized += 1;
        }
        return .{ .allocator = allocator, .per_cell = cells };
    }

    pub fn deinit(self: *GridWorkspace) void {
        for (self.per_cell) |*cell| cell.deinit();
        self.allocator.free(self.per_cell);
        self.* = undefined;
    }
};

pub fn calculate(input: Input) !Result {
    inline for (@typeInfo(Input).@"struct".fields) |field| if (!std.math.isFinite(@field(input, field.name))) return error.NonFiniteRootExudationInput;
    if (input.soil_water_volume_m3 <= 0 or input.root_water_volume_m3 <= 0 or input.soil_dissolved_carbon_g_c < 0 or input.soil_dissolved_nitrogen_g_n < 0 or input.soil_dissolved_phosphorus_g_p < 0 or input.root_nonstructural_carbon_g_c < 0 or input.root_nonstructural_nitrogen_g_n < 0 or input.root_nonstructural_phosphorus_g_p < 0 or input.maximum_root_carbon_concentration_g_c_per_m3 <= 0 or input.root_mobile_nitrogen_exchange_fraction < 0 or input.root_mobile_nitrogen_exchange_fraction > 1 or input.root_mobile_phosphorus_exchange_fraction < 0 or input.root_mobile_phosphorus_exchange_fraction > 1 or input.exchange_rate_per_h < 0 or input.timestep_h < 0 or input.significance_threshold_g < 0) return error.InvalidRootExudationInput;
    const total_water = input.soil_water_volume_m3 + input.root_water_volume_m3;
    const exchange_fraction = input.exchange_rate_per_h * input.timestep_h;
    const exchangeable_root_carbon = @min(input.maximum_root_carbon_concentration_g_c_per_m3 * input.root_water_volume_m3, input.root_nonstructural_carbon_g_c);
    const carbon = exchange_fraction * (input.soil_dissolved_carbon_g_c * input.root_water_volume_m3 - exchangeable_root_carbon * input.soil_water_volume_m3) / total_water;
    var nitrogen: f64 = 0;
    var phosphorus: f64 = 0;
    if (input.soil_dissolved_carbon_g_c > input.significance_threshold_g and input.root_nonstructural_carbon_g_c > input.significance_threshold_g) {
        const total_carbon = input.soil_dissolved_carbon_g_c + input.root_nonstructural_carbon_g_c;
        const exchangeable_nitrogen = input.root_mobile_nitrogen_exchange_fraction * input.root_nonstructural_nitrogen_g_n;
        const exchangeable_phosphorus = input.root_mobile_phosphorus_exchange_fraction * input.root_nonstructural_phosphorus_g_p;
        nitrogen = exchange_fraction * (input.soil_dissolved_nitrogen_g_n * input.root_nonstructural_carbon_g_c - exchangeable_nitrogen * input.soil_dissolved_carbon_g_c) / total_carbon;
        phosphorus = exchange_fraction * (input.soil_dissolved_phosphorus_g_p * input.root_nonstructural_carbon_g_c - exchangeable_phosphorus * input.soil_dissolved_carbon_g_c) / total_carbon;
    }
    const result = Result{ .carbon_g_c = carbon, .nitrogen_g_n = nitrogen, .phosphorus_g_p = phosphorus };
    inline for (@typeInfo(Result).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteRootExudationResult;
    return result;
}

/// Atomic RDFOM*/XOQ* publication. Positive exchange removes dissolved
/// material from soil and adds it to the root mobile pool; negative exchange
/// is exudation. The soil change ledger therefore receives the opposite sign.
pub fn commit(
    soil_dissolved_carbon_g_c: *f64,
    soil_dissolved_nitrogen_g_n: *f64,
    soil_dissolved_phosphorus_g_p: *f64,
    root_mobile_carbon_g_c: *f64,
    root_mobile_nitrogen_g_n: *f64,
    root_mobile_phosphorus_g_p: *f64,
    soil_carbon_change_g_c: *f64,
    soil_nitrogen_change_g_n: *f64,
    soil_phosphorus_change_g_p: *f64,
    exchange: Result,
) !void {
    const current = .{ soil_dissolved_carbon_g_c.*, soil_dissolved_nitrogen_g_n.*, soil_dissolved_phosphorus_g_p.*, root_mobile_carbon_g_c.*, root_mobile_nitrogen_g_n.*, root_mobile_phosphorus_g_p.*, soil_carbon_change_g_c.*, soil_nitrogen_change_g_n.*, soil_phosphorus_change_g_p.* };
    inline for (current) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootExudationCommitInput;
    inline for (@typeInfo(Result).@"struct".fields) |field| if (!std.math.isFinite(@field(exchange, field.name))) return error.NonFiniteRootExudationResult;
    inline for (.{ soil_dissolved_carbon_g_c.*, soil_dissolved_nitrogen_g_n.*, soil_dissolved_phosphorus_g_p.*, root_mobile_carbon_g_c.*, root_mobile_nitrogen_g_n.*, root_mobile_phosphorus_g_p.* }) |value| if (value < 0) return error.InvalidRootExudationCommitInput;

    const next_soil_c = soil_dissolved_carbon_g_c.* - exchange.carbon_g_c;
    const next_soil_n = soil_dissolved_nitrogen_g_n.* - exchange.nitrogen_g_n;
    const next_soil_p = soil_dissolved_phosphorus_g_p.* - exchange.phosphorus_g_p;
    const next_root_c = root_mobile_carbon_g_c.* + exchange.carbon_g_c;
    const next_root_n = root_mobile_nitrogen_g_n.* + exchange.nitrogen_g_n;
    const next_root_p = root_mobile_phosphorus_g_p.* + exchange.phosphorus_g_p;
    const next_change_c = soil_carbon_change_g_c.* - exchange.carbon_g_c;
    const next_change_n = soil_nitrogen_change_g_n.* - exchange.nitrogen_g_n;
    const next_change_p = soil_phosphorus_change_g_p.* - exchange.phosphorus_g_p;
    const next = .{ next_soil_c, next_soil_n, next_soil_p, next_root_c, next_root_n, next_root_p, next_change_c, next_change_n, next_change_p };
    inline for (next) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootExudationCommit;
    inline for (.{ next_soil_c, next_soil_n, next_soil_p, next_root_c, next_root_n, next_root_p }) |value| if (value < -1.0e-12) return error.InsufficientElementForRootExudation;

    soil_dissolved_carbon_g_c.* = @max(0.0, next_soil_c);
    soil_dissolved_nitrogen_g_n.* = @max(0.0, next_soil_n);
    soil_dissolved_phosphorus_g_p.* = @max(0.0, next_soil_p);
    root_mobile_carbon_g_c.* = @max(0.0, next_root_c);
    root_mobile_nitrogen_g_n.* = @max(0.0, next_root_n);
    root_mobile_phosphorus_g_p.* = @max(0.0, next_root_p);
    soil_carbon_change_g_c.* = next_change_c;
    soil_nitrogen_change_g_n.* = next_change_n;
    soil_phosphorus_change_g_p.* = next_change_p;
}

pub fn sumSubstrateExchange(carbon_g_c: []const f64, nitrogen_g_n: []const f64, phosphorus_g_p: []const f64) !Result {
    if (carbon_g_c.len != 5 or nitrogen_g_n.len != 5 or phosphorus_g_p.len != 5) return error.RootExudationSubstrateCountMismatch;
    var total: Result = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
    for (carbon_g_c, nitrogen_g_n, phosphorus_g_p) |carbon, nitrogen, phosphorus| {
        inline for (.{ carbon, nitrogen, phosphorus }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootExudationResult;
        total.carbon_g_c += carbon;
        total.nitrogen_g_n += nitrogen;
        total.phosphorus_g_p += phosphorus;
    }
    inline for (@typeInfo(Result).@"struct".fields) |field| if (!std.math.isFinite(@field(total, field.name))) return error.NonFiniteRootExudationResult;
    return total;
}

test "UPTAKE C N P exudation retains concentration and stoichiometric equations" {
    const input = Input{ .soil_water_volume_m3 = 2, .root_water_volume_m3 = 1, .soil_dissolved_carbon_g_c = 3, .soil_dissolved_nitrogen_g_n = 0.4, .soil_dissolved_phosphorus_g_p = 0.08, .root_nonstructural_carbon_g_c = 6, .root_nonstructural_nitrogen_g_n = 1, .root_nonstructural_phosphorus_g_p = 0.2, .maximum_root_carbon_concentration_g_c_per_m3 = 1000, .root_mobile_nitrogen_exchange_fraction = 0.1, .root_mobile_phosphorus_exchange_fraction = 0.1, .exchange_rate_per_h = 0.2, .timestep_h = 0.5, .significance_threshold_g = 1.0e-12 };
    const result = try calculate(input);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1) * (3.0 * 1.0 - 6.0 * 2.0) / 3.0, result.carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1) * (0.4 * 6.0 - 0.1 * 3.0) / 9.0, result.nitrogen_g_n, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1) * (0.08 * 6.0 - 0.02 * 3.0) / 9.0, result.phosphorus_g_p, 1.0e-12);
}

test "UPTAKE exudation commit conserves C N P and reverses soil ledger sign" {
    var soil_c: f64 = 2;
    var soil_n: f64 = 1;
    var soil_p: f64 = 0.5;
    var root_c: f64 = 3;
    var root_n: f64 = 0.4;
    var root_p: f64 = 0.2;
    var change_c: f64 = 0;
    var change_n: f64 = 0;
    var change_p: f64 = 0;
    const before = .{ soil_c + root_c, soil_n + root_n, soil_p + root_p };
    try commit(&soil_c, &soil_n, &soil_p, &root_c, &root_n, &root_p, &change_c, &change_n, &change_p, .{ .carbon_g_c = -0.5, .nitrogen_g_n = 0.2, .phosphorus_g_p = -0.1 });
    try std.testing.expectApproxEqAbs(before[0], soil_c + root_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(before[1], soil_n + root_n, 1.0e-12);
    try std.testing.expectApproxEqAbs(before[2], soil_p + root_p, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), change_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -0.2), change_n, 1.0e-12);
    const total = try sumSubstrateExchange(&.{ -0.5, 1, 2, 3, 4 }, &.{ 0.2, 0, 0, 0, 0 }, &.{ -0.1, 0, 0, 0, 0 });
    try std.testing.expectApproxEqAbs(@as(f64, 9.5), total.carbon_g_c, 1.0e-12);
}

test "five-substrate exudation workspace supports runtime competitors and rolls back atomically" {
    var roots = try RootState.init(std.testing.allocator, 7, 1, 1);
    defer roots.deinit();
    var soil = try SoilState.init(std.testing.allocator, 1);
    defer soil.deinit();
    var workspace = try Workspace.init(std.testing.allocator, 7);
    defer workspace.deinit();
    for (0..substrate_count) |substrate| {
        soil.dissolved[substrate] = .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 };
    }
    for (0..7) |plant| {
        const root = try roots.layerIndex(plant, 0, 0);
        roots.aqueous_volume_m3[root] = 0.1;
        roots.mobile_carbon_g[root] = 2;
        roots.mobile_nitrogen_g[root] = 0.2;
        roots.mobile_phosphorus_g[root] = 0.02;
        workspace.competitors[plant] = .{ .plant = plant, .domain = 0, .layer = 0 };
    }
    const fractions = [_]f64{0.2} ** substrate_count;
    var carbon_before: f64 = 0;
    for (soil.dissolved) |pool| carbon_before += pool.carbon_g_c;
    for (roots.mobile_carbon_g) |value| carbon_before += value;
    try workspace.stageLayer(&roots, &soil, 0, 1, &fractions, compatibilityParameters(), 1.0e-12, 7);
    try workspace.commitLayer(&roots, &soil, 0, 7);
    var carbon_after: f64 = 0;
    for (soil.dissolved) |pool| carbon_after += pool.carbon_g_c;
    for (roots.mobile_carbon_g) |value| carbon_after += value;
    try std.testing.expectApproxEqAbs(carbon_before, carbon_after, 1.0e-12);

    const soil_before_failure = soil.dissolved[0].carbon_g_c;
    const root_before_failure = roots.mobile_carbon_g[0];
    for (0..7) |competitor| workspace.staged_results[competitor * substrate_count] = .{ .carbon_g_c = 10, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
    try std.testing.expectError(error.InsufficientElementForRootExudation, workspace.commitLayer(&roots, &soil, 0, 7));
    try std.testing.expectEqual(soil_before_failure, soil.dissolved[0].carbon_g_c);
    try std.testing.expectEqual(root_before_failure, roots.mobile_carbon_g[0]);
}
