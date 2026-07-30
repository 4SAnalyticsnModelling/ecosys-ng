const std = @import("std");
const organic = @import("soil_organic_initialization.zig");
const parameters_module = @import("soil_organic_parameters.zig");

pub const Kind = enum { plant_residue, manure };

/// HOUR1 OFC/OFN/OFP application. The event is first partitioned into local
/// deltas, every destination is validated, and only then is live state changed.
pub fn apply(
    state: *organic.State,
    layer: usize,
    kind: Kind,
    material_type: u8,
    input: organic.ElementPool,
    parameters: *const parameters_module.OwnedParameters,
) !void {
    if (layer >= state.layer_count) return error.OrganicFertilizerLayerOutOfBounds;
    inline for (.{ input.carbon_g_c, input.nitrogen_g_n, input.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidOrganicFertilizerInput;
    if (input.carbon_g_c == 0 and input.nitrogen_g_n == 0 and input.phosphorus_g_p == 0) return;

    const target_substrate: usize = switch (kind) {
        .plant_residue => 1,
        .manure => 2,
    };
    if (kind == .plant_residue and material_type == 10) {
        const charcoal_index = (layer * organic.substrate_count + 3) * organic.structural_fraction_count + 4;
        try validateAddition(state.structural[charcoal_index], input);
        add(&state.structural[charcoal_index], input);
        return;
    }

    var microbial_delta = [_]organic.ElementPool{.{}} ** (organic.microbial_substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count);
    var allocated: organic.ElementPool = .{};
    for (0..organic.microbial_population_count) |population| for (0..organic.kinetic_fraction_count) |fraction| {
        const parameter_index = (target_substrate * organic.microbial_population_count + population) * organic.kinetic_fraction_count + fraction;
        const carbon = @max(0, @min(
            input.carbon_g_c * parameters.microbial_kinetic_fraction[target_substrate * organic.kinetic_fraction_count + fraction] * parameters.heterotroph_population_fraction[population],
            input.carbon_g_c - allocated.carbon_g_c,
        ));
        const pool: organic.ElementPool = .{
            .carbon_g_c = carbon,
            .nitrogen_g_n = @max(0, @min(carbon * parameters.microbial_nitrogen_to_carbon[parameter_index], input.nitrogen_g_n - allocated.nitrogen_g_n)),
            .phosphorus_g_p = @max(0, @min(carbon * parameters.microbial_phosphorus_to_carbon[parameter_index], input.phosphorus_g_p - allocated.phosphorus_g_p)),
        };
        add(&microbial_delta[microbialIndex(target_substrate, population, fraction)], pool);
        add(&allocated, pool);
        for (0..organic.microbial_population_count) |autotroph_population| {
            const autotroph = scale(pool, parameters.autotroph_population_fraction[autotroph_population]);
            add(&microbial_delta[microbialIndex(organic.autotrophic_substrate_index, autotroph_population, fraction)], autotroph);
            add(&allocated, autotroph);
        }
    };

    const dissolved: organic.ElementPool = .{
        .carbon_g_c = @max(0, @min(0.1 * allocated.carbon_g_c, input.carbon_g_c - allocated.carbon_g_c)),
        .nitrogen_g_n = @max(0, @min(0.1 * allocated.nitrogen_g_n, input.nitrogen_g_n - allocated.nitrogen_g_n)),
        .phosphorus_g_p = @max(0, @min(0.1 * allocated.phosphorus_g_p, input.phosphorus_g_p - allocated.phosphorus_g_p)),
    };
    add(&allocated, dissolved);

    const fractions = switch (kind) {
        .plant_residue => blk: {
            const type_index: usize = if (material_type >= 1 and material_type <= 10) material_type - 1 else 11;
            break :blk parameters.surface_plant_structural_fraction[type_index * organic.structural_fraction_count ..][0..organic.structural_fraction_count];
        },
        .manure => blk: {
            const type_index: usize = if (material_type == 1 or material_type == 3) 0 else 1;
            break :blk parameters.surface_manure_structural_fraction[type_index * organic.structural_fraction_count ..][0..organic.structural_fraction_count];
        },
    };
    const nutrient_weight_row: usize = if (kind == .plant_residue) 1 else 2;
    const nitrogen_weights = parameters.surface_residue_nitrogen_weight[nutrient_weight_row * organic.structural_fraction_count ..][0..organic.structural_fraction_count];
    const phosphorus_weights = parameters.surface_residue_phosphorus_weight[nutrient_weight_row * organic.structural_fraction_count ..][0..organic.structural_fraction_count];
    var weighted_nitrogen: f64 = 0;
    var weighted_phosphorus: f64 = 0;
    for (0..organic.structural_fraction_count) |fraction| {
        weighted_nitrogen += fractions[fraction] * nitrogen_weights[fraction];
        weighted_phosphorus += fractions[fraction] * phosphorus_weights[fraction];
    }
    const remaining: organic.ElementPool = .{
        .carbon_g_c = @max(0, input.carbon_g_c - allocated.carbon_g_c),
        .nitrogen_g_n = @max(0, input.nitrogen_g_n - allocated.nitrogen_g_n),
        .phosphorus_g_p = @max(0, input.phosphorus_g_p - allocated.phosphorus_g_p),
    };
    var structural_delta = [_]organic.ElementPool{.{}} ** organic.structural_fraction_count;
    for (0..organic.structural_fraction_count) |fraction| structural_delta[fraction] = .{
        .carbon_g_c = fractions[fraction] * remaining.carbon_g_c,
        .nitrogen_g_n = if (weighted_nitrogen > 0) fractions[fraction] * nitrogen_weights[fraction] / weighted_nitrogen * remaining.nitrogen_g_n else 0,
        .phosphorus_g_p = if (weighted_phosphorus > 0) fractions[fraction] * phosphorus_weights[fraction] / weighted_phosphorus * remaining.phosphorus_g_p else 0,
    };

    const microbial_first = layer * organic.microbial_substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count;
    for (microbial_delta, 0..) |delta, index| try validateAddition(state.microbial[microbial_first + index], delta);
    const dissolved_index = layer * organic.substrate_count + target_substrate;
    try validateAddition(state.dissolved[dissolved_index], dissolved);
    for (structural_delta, 0..) |delta, fraction| {
        const index = (layer * organic.substrate_count + target_substrate) * organic.structural_fraction_count + fraction;
        try validateAddition(state.structural[index], delta);
        const colonized = state.colonized_structural_carbon_g_c[index] + delta.carbon_g_c * parameters.microbial_kinetic_fraction[target_substrate * organic.kinetic_fraction_count];
        if (!std.math.isFinite(colonized)) return error.OrganicFertilizerApplicationOverflow;
    }
    for (microbial_delta, 0..) |delta, index| add(&state.microbial[microbial_first + index], delta);
    add(&state.dissolved[dissolved_index], dissolved);
    for (structural_delta, 0..) |delta, fraction| {
        const index = (layer * organic.substrate_count + target_substrate) * organic.structural_fraction_count + fraction;
        add(&state.structural[index], delta);
        state.colonized_structural_carbon_g_c[index] += delta.carbon_g_c * parameters.microbial_kinetic_fraction[target_substrate * organic.kinetic_fraction_count];
    }
}

fn microbialIndex(substrate: usize, population: usize, fraction: usize) usize {
    return (substrate * organic.microbial_population_count + population) * organic.kinetic_fraction_count + fraction;
}
fn add(target: *organic.ElementPool, value: organic.ElementPool) void {
    target.carbon_g_c += value.carbon_g_c;
    target.nitrogen_g_n += value.nitrogen_g_n;
    target.phosphorus_g_p += value.phosphorus_g_p;
}
fn scale(pool: organic.ElementPool, fraction: f64) organic.ElementPool {
    return .{ .carbon_g_c = pool.carbon_g_c * fraction, .nitrogen_g_n = pool.nitrogen_g_n * fraction, .phosphorus_g_p = pool.phosphorus_g_p * fraction };
}
fn validateAddition(current: organic.ElementPool, delta: organic.ElementPool) !void {
    inline for (.{ current.carbon_g_c + delta.carbon_g_c, current.nitrogen_g_n + delta.nitrogen_g_n, current.phosphorus_g_p + delta.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.OrganicFertilizerApplicationOverflow;
}

test "HOUR1 organic application conserves arbitrary residue C N P" {
    var state = try organic.State.init(std.testing.allocator, 2);
    defer state.deinit();
    var parameters = try parameters_module.sourceParameters(std.testing.allocator);
    defer parameters.deinit();
    const input: organic.ElementPool = .{ .carbon_g_c = 100, .nitrogen_g_n = 5, .phosphorus_g_p = 1 };
    try apply(&state, 1, .plant_residue, 2, input, &parameters);
    var carbon: f64 = 0;
    var nitrogen: f64 = 0;
    var phosphorus: f64 = 0;
    const microbial_first = organic.microbial_substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count;
    for (state.microbial[microbial_first .. microbial_first * 2]) |pool| {
        carbon += pool.carbon_g_c;
        nitrogen += pool.nitrogen_g_n;
        phosphorus += pool.phosphorus_g_p;
    }
    for (state.dissolved[organic.substrate_count .. organic.substrate_count * 2]) |pool| {
        carbon += pool.carbon_g_c;
        nitrogen += pool.nitrogen_g_n;
        phosphorus += pool.phosphorus_g_p;
    }
    for (state.structural[organic.substrate_count * organic.structural_fraction_count .. organic.substrate_count * organic.structural_fraction_count * 2]) |pool| {
        carbon += pool.carbon_g_c;
        nitrogen += pool.nitrogen_g_n;
        phosphorus += pool.phosphorus_g_p;
    }
    try std.testing.expectApproxEqAbs(input.carbon_g_c, carbon, 1e-10);
    try std.testing.expectApproxEqAbs(input.nitrogen_g_n, nitrogen, 1e-10);
    try std.testing.expectApproxEqAbs(input.phosphorus_g_p, phosphorus, 1e-10);
}

test "HOUR1 production parameters preserve plant material rows four five and eight" {
    var parameters = try parameters_module.sourceParameters(std.testing.allocator);
    defer parameters.deinit();

    const expected = @import("organic_fertilizer_material_fractions.zig").sourceParameters();
    inline for (.{ 4, 5, 8 }) |material_type| {
        const source = try @import("organic_fertilizer_material_fractions.zig").plant(material_type, expected);
        const offset = (material_type - 1) * organic.structural_fraction_count;
        const actual = parameters.surface_plant_structural_fraction[offset..][0..4];
        try std.testing.expectEqual(source.protein, actual[0]);
        try std.testing.expectEqual(source.soluble_carbohydrate, actual[1]);
        try std.testing.expectEqual(source.cellulose, actual[2]);
        try std.testing.expectEqual(source.lignin, actual[3]);
    }
}

test "HOUR1 manure type three selects runtime ruminant fractions" {
    var state = try organic.State.init(std.testing.allocator, 1);
    defer state.deinit();
    var parameters = try parameters_module.sourceParameters(std.testing.allocator);
    defer parameters.deinit();

    // Disable the preliminary microbial/dissolved allocation so structural
    // publication exposes only the material-fraction selection.
    @memset(parameters.microbial_kinetic_fraction, 0);
    const input: organic.ElementPool = .{ .carbon_g_c = 1 };
    try apply(&state, 0, .manure, 3, input, &parameters);

    const source = try @import("organic_fertilizer_material_fractions.zig").manure(
        3,
        @import("organic_fertilizer_material_fractions.zig").sourceParameters(),
    );
    const first = (2 * organic.structural_fraction_count);
    try std.testing.expectEqual(source.protein, state.structural[first].carbon_g_c);
    try std.testing.expectEqual(source.soluble_carbohydrate, state.structural[first + 1].carbon_g_c);
    try std.testing.expectEqual(source.cellulose, state.structural[first + 2].carbon_g_c);
    try std.testing.expectEqual(source.lignin, state.structural[first + 3].carbon_g_c);
}
