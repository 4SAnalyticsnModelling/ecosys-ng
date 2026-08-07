const std = @import("std");
const organic = @import("../organic/initialization.zig");
const constituents = @import("../../erosion/eroded_constituents.zig");

pub const ExportedElementsG = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub fn componentCount(state: *const organic.State) !usize {
    if (state.layer_count == 0 or state.microbial.len % state.layer_count != 0 or state.residue.len % state.layer_count != 0 or state.adsorbed.len % state.layer_count != 0 or state.adsorbed_acetate_carbon_g_c.len % state.layer_count != 0 or state.structural.len % state.layer_count != 0 or state.colonized_structural_carbon_g_c.len % state.layer_count != 0) return error.OrganicErosionDimensionMismatch;
    return 3 * (state.microbial.len / state.layer_count + state.residue.len / state.layer_count + state.adsorbed.len / state.layer_count + state.structural.len / state.layer_count) +
        state.adsorbed_acetate_carbon_g_c.len / state.layer_count +
        state.colonized_structural_carbon_g_c.len / state.layer_count;
}

pub fn route(
    columns: usize,
    rows: usize,
    soil_layer_capacity: usize,
    surface_soil_mass_megagrams: []const f64,
    state: *organic.State,
    sediment: constituents.DirectionalSediment,
    workspace: *constituents.PackedWorkspace,
) !void {
    const cells = try std.math.mul(usize, columns, rows);
    if (soil_layer_capacity == 0 or state.layer_count != try std.math.mul(usize, cells, soil_layer_capacity)) return error.OrganicErosionDimensionMismatch;
    const component_count = try componentCount(state);
    if (workspace.cell_count != cells or workspace.component_count != component_count) return error.OrganicErosionDimensionMismatch;
    try packSurface(state, cells, soil_layer_capacity, component_count, workspace.pools);
    try constituents.routePackedWorkspace(workspace, columns, rows, surface_soil_mass_megagrams, sediment);
    try unpackSurface(state, cells, soil_layer_capacity, component_count, workspace.pools);
}

/// Refreshes HOUR1's ORGC concentration from the authoritative extensive
/// organic pools. Colonized structural carbon is an annotation of structural
/// carbon and is therefore not counted a second time.
pub fn refreshSurfaceOrganicCarbonGPerMg(
    state: *const organic.State,
    soil_layer_capacity: usize,
    surface_soil_mass_megagrams: []const f64,
    total_organic_carbon_g_per_megagram: []f64,
) !void {
    const cells = surface_soil_mass_megagrams.len;
    if (soil_layer_capacity == 0 or state.layer_count != try std.math.mul(usize, cells, soil_layer_capacity) or total_organic_carbon_g_per_megagram.len != state.layer_count) return error.OrganicErosionDimensionMismatch;
    for (0..cells) |cell| {
        const layer = cell * soil_layer_capacity;
        const soil_mass_megagrams = surface_soil_mass_megagrams[cell];
        if (!std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams <= 0) return error.InvalidOrganicErosionState;
        var carbon_g: f64 = 0;
        inline for (.{ state.microbial, state.residue, state.adsorbed }) |pools| {
            const per_layer = pools.len / state.layer_count;
            for (pools[layer * per_layer ..][0..per_layer]) |pool| carbon_g += pool.carbon_g_c;
        }
        const acetate_per_layer = state.adsorbed_acetate_carbon_g_c.len / state.layer_count;
        for (state.adsorbed_acetate_carbon_g_c[layer * acetate_per_layer ..][0..acetate_per_layer]) |value| carbon_g += value;
        const structural_per_layer = state.structural.len / state.layer_count;
        for (state.structural[layer * structural_per_layer ..][0..structural_per_layer]) |pool| carbon_g += pool.carbon_g_c;
        if (!std.math.isFinite(carbon_g) or carbon_g < 0) return error.InvalidOrganicErosionState;
        total_organic_carbon_g_per_megagram[layer] = carbon_g / soil_mass_megagrams;
    }
}

/// REDIST `COE/ZOE/POE` external sediment loss. Colonized structural carbon
/// remains an annotation of its structural pool and is not counted twice.
pub fn exportedElements(
    state: *const organic.State,
    workspace: *const constituents.PackedWorkspace,
) !ExportedElementsG {
    const components = try componentCount(state);
    if (workspace.component_count != components or
        workspace.exported.len !=
            try std.math.mul(usize, workspace.cell_count, components))
        return error.OrganicErosionDimensionMismatch;
    var result: ExportedElementsG = .{};
    for (0..workspace.cell_count) |cell| {
        var cursor = cell * components;
        try sumElementExports(
            state.microbial.len / state.layer_count,
            workspace.exported,
            &cursor,
            &result,
        );
        try sumElementExports(
            state.residue.len / state.layer_count,
            workspace.exported,
            &cursor,
            &result,
        );
        try sumElementExports(
            state.adsorbed.len / state.layer_count,
            workspace.exported,
            &cursor,
            &result,
        );
        const acetate_count =
            state.adsorbed_acetate_carbon_g_c.len / state.layer_count;
        for (0..acetate_count) |_| {
            result.carbon_g_c = try addExport(
                result.carbon_g_c,
                workspace.exported[cursor],
            );
            cursor += 1;
        }
        try sumElementExports(
            state.structural.len / state.layer_count,
            workspace.exported,
            &cursor,
            &result,
        );
        cursor += state.colonized_structural_carbon_g_c.len / state.layer_count;
        if (cursor != (cell + 1) * components)
            return error.OrganicErosionDimensionMismatch;
    }
    return result;
}

fn sumElementExports(
    pool_count: usize,
    values: []const f64,
    cursor: *usize,
    result: *ExportedElementsG,
) !void {
    for (0..pool_count) |_| {
        result.carbon_g_c = try addExport(result.carbon_g_c, values[cursor.*]);
        cursor.* += 1;
        result.nitrogen_g_n = try addExport(result.nitrogen_g_n, values[cursor.*]);
        cursor.* += 1;
        result.phosphorus_g_p = try addExport(result.phosphorus_g_p, values[cursor.*]);
        cursor.* += 1;
    }
}

fn addExport(total: f64, value: f64) !f64 {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidOrganicErosionExport;
    const next = total + value;
    if (!std.math.isFinite(next)) return error.OrganicErosionExportOverflow;
    return next;
}

fn packSurface(state: *const organic.State, cells: usize, soil_layer_capacity: usize, component_count: usize, output: []f64) !void {
    for (0..cells) |cell| {
        const layer = cell * soil_layer_capacity;
        var cursor = cell * component_count;
        try packElementLayer(state.microbial, state.layer_count, layer, output, &cursor);
        try packElementLayer(state.residue, state.layer_count, layer, output, &cursor);
        try packElementLayer(state.adsorbed, state.layer_count, layer, output, &cursor);
        try packScalarLayer(state.adsorbed_acetate_carbon_g_c, state.layer_count, layer, output, &cursor);
        try packElementLayer(state.structural, state.layer_count, layer, output, &cursor);
        try packScalarLayer(state.colonized_structural_carbon_g_c, state.layer_count, layer, output, &cursor);
        if (cursor != (cell + 1) * component_count) return error.OrganicErosionDimensionMismatch;
    }
}

fn unpackSurface(state: *organic.State, cells: usize, soil_layer_capacity: usize, component_count: usize, input: []const f64) !void {
    for (input) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidOrganicErosionCandidate;
    for (0..cells) |cell| {
        const layer = cell * soil_layer_capacity;
        var cursor = cell * component_count;
        try unpackElementLayer(state.microbial, state.layer_count, layer, input, &cursor);
        try unpackElementLayer(state.residue, state.layer_count, layer, input, &cursor);
        try unpackElementLayer(state.adsorbed, state.layer_count, layer, input, &cursor);
        try unpackScalarLayer(state.adsorbed_acetate_carbon_g_c, state.layer_count, layer, input, &cursor);
        try unpackElementLayer(state.structural, state.layer_count, layer, input, &cursor);
        try unpackScalarLayer(state.colonized_structural_carbon_g_c, state.layer_count, layer, input, &cursor);
    }
}

fn packElementLayer(values: []const organic.ElementPool, layer_count: usize, layer: usize, output: []f64, cursor: *usize) !void {
    const per_layer = values.len / layer_count;
    for (values[layer * per_layer ..][0..per_layer]) |pool| inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| {
        const value = @field(pool, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidOrganicErosionState;
        output[cursor.*] = value;
        cursor.* += 1;
    };
}

fn unpackElementLayer(values: []organic.ElementPool, layer_count: usize, layer: usize, input: []const f64, cursor: *usize) !void {
    const per_layer = values.len / layer_count;
    for (values[layer * per_layer ..][0..per_layer]) |*pool| inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| {
        @field(pool, field.name) = input[cursor.*];
        cursor.* += 1;
    };
}

fn packScalarLayer(values: []const f64, layer_count: usize, layer: usize, output: []f64, cursor: *usize) !void {
    const per_layer = values.len / layer_count;
    for (values[layer * per_layer ..][0..per_layer]) |value| {
        if (!std.math.isFinite(value) or value < 0) return error.InvalidOrganicErosionState;
        output[cursor.*] = value;
        cursor.* += 1;
    }
}

fn unpackScalarLayer(values: []f64, layer_count: usize, layer: usize, input: []const f64, cursor: *usize) !void {
    const per_layer = values.len / layer_count;
    for (values[layer * per_layer ..][0..per_layer]) |*value| {
        value.* = input[cursor.*];
        cursor.* += 1;
    }
}

test "all solid organic families follow sediment conservatively" {
    var state = try organic.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.microbial[0].carbon_g_c = 10;
    state.residue[0].nitrogen_g_n = 20;
    state.adsorbed[0].phosphorus_g_p = 30;
    state.adsorbed_acetate_carbon_g_c[0] = 40;
    state.structural[0].carbon_g_c = 50;
    state.colonized_structural_carbon_g_c[0] = 60;
    const count = try componentCount(&state);
    var workspace = try constituents.PackedWorkspace.init(std.testing.allocator, 2, count);
    defer workspace.deinit();
    try route(2, 1, 1, &.{ 10, 10 }, &state, .{ .east_megagrams = &.{ 1, 0 }, .west_megagrams = &.{ 0, 0 }, .south_megagrams = &.{ 0, 0 }, .north_megagrams = &.{ 0, 0 } }, &workspace);
    try std.testing.expectApproxEqAbs(@as(f64, 9), state.microbial[0].carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.microbial[state.microbial.len / 2].carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 54), state.colonized_structural_carbon_g_c[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 6), state.colonized_structural_carbon_g_c[state.colonized_structural_carbon_g_c.len / 2], 1e-14);
    var total_organic_carbon_g_per_megagram = [_]f64{ 0, 0 };
    try refreshSurfaceOrganicCarbonGPerMg(&state, 1, &.{ 10, 10 }, &total_organic_carbon_g_per_megagram);
    try std.testing.expectApproxEqAbs(@as(f64, 9), total_organic_carbon_g_per_megagram[0], 1e-14);
}

test "external organic sediment export counts C N P without colonized double count" {
    var state = try organic.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const count = try componentCount(&state);
    var workspace = try constituents.PackedWorkspace.init(
        std.testing.allocator,
        1,
        count,
    );
    defer workspace.deinit();
    state.microbial[0] = .{
        .carbon_g_c = 10,
        .nitrogen_g_n = 20,
        .phosphorus_g_p = 30,
    };
    state.adsorbed_acetate_carbon_g_c[0] = 40;
    state.structural[0] = .{
        .carbon_g_c = 50,
        .nitrogen_g_n = 60,
        .phosphorus_g_p = 70,
    };
    state.colonized_structural_carbon_g_c[0] = 80;
    try route(1, 1, 1, &.{10}, &state, .{
        .east_megagrams = &.{1},
        .west_megagrams = &.{0},
        .south_megagrams = &.{0},
        .north_megagrams = &.{0},
    }, &workspace);
    const exported = try exportedElements(&state, &workspace);
    try std.testing.expectApproxEqAbs(@as(f64, 10), exported.carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 8), exported.nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10), exported.phosphorus_g_p, 1e-12);
}
