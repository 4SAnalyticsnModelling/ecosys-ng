const std = @import("std");
const inventory_module = @import("fertilizer_nitrogen_inventory.zig");
const fertilizer = @import("soil_fertilizer_dissolution.zig");
const constituents = @import("eroded_constituents.zig");

pub const component_count: usize = @typeInfo(fertilizer.FertilizerState).@"struct".fields.len;

pub const Exported = struct {
    nitrogen_g_n: f64,
    ion_mol: f64,
};

pub fn route(
    columns: usize,
    rows: usize,
    surface_soil_mass_Mg: []const f64,
    inventory: *inventory_module.State,
    sediment: constituents.DirectionalSediment,
    workspace: *constituents.PackedWorkspace,
) !void {
    const cells = try std.math.mul(usize, columns, rows);
    if (inventory.cell_count != cells or workspace.cell_count != cells or workspace.component_count != component_count) return error.FertilizerErosionDimensionMismatch;
    for (0..cells) |cell| {
        const top = try inventory.index(cell, 0);
        inline for (@typeInfo(fertilizer.FertilizerState).@"struct".fields, 0..) |field, component| {
            const value = @field(inventory.soil[top], field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidFertilizerErosionState;
            workspace.pools[cell * component_count + component] = value;
        }
    }
    try constituents.routePackedWorkspace(workspace, columns, rows, surface_soil_mass_Mg, sediment);
    for (workspace.pools) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidFertilizerErosionCandidate;
    for (0..cells) |cell| {
        const top = try inventory.index(cell, 0);
        inline for (@typeInfo(fertilizer.FertilizerState).@"struct".fields, 0..) |field, component| @field(inventory.soil[top], field.name) = workspace.pools[cell * component_count + component];
    }
}

/// REDIST `ZPE/SEF` fertilizer carried through external sediment faces.
pub fn exported(
    workspace: *const constituents.PackedWorkspace,
    nitrogen_g_per_mol: f64,
) !Exported {
    if (workspace.component_count != component_count or
        workspace.exported.len !=
            try std.math.mul(usize, workspace.cell_count, component_count))
        return error.FertilizerErosionDimensionMismatch;
    if (!std.math.isFinite(nitrogen_g_per_mol) or nitrogen_g_per_mol <= 0)
        return error.InvalidFertilizerErosionMolarMass;
    var nitrogen_mol: f64 = 0;
    var ion_mol: f64 = 0;
    for (0..workspace.cell_count) |cell| {
        const first = cell * component_count;
        inline for (@typeInfo(fertilizer.FertilizerState).@"struct".fields, 0..) |field, component| {
            const amount = workspace.exported[first + component];
            if (!std.math.isFinite(amount) or amount < 0)
                return error.InvalidFertilizerErosionExport;
            nitrogen_mol += amount;
            ion_mol += amount * if (comptime std.mem.indexOf(u8, field.name, "ammonium") != null) @as(f64, 2) else @as(f64, 1);
        }
    }
    const nitrogen_g_n = nitrogen_mol * nitrogen_g_per_mol;
    if (!std.math.isFinite(nitrogen_g_n) or !std.math.isFinite(ion_mol))
        return error.FertilizerErosionExportOverflow;
    return .{ .nitrogen_g_n = nitrogen_g_n, .ion_mol = ion_mol };
}

test "broadcast and banded fertilizer amounts follow sediment faces" {
    var inventory = try inventory_module.State.init(std.testing.allocator, 2, 1);
    defer inventory.deinit();
    inventory.soil[0].broadcast_ammonium_mol_n = 10;
    inventory.soil[0].banded_nitrate_mol_n = 20;
    var workspace = try constituents.PackedWorkspace.init(std.testing.allocator, 2, component_count);
    defer workspace.deinit();
    try route(2, 1, &.{ 10, 10 }, &inventory, .{ .east_Mg = &.{ 1, 0 }, .west_Mg = &.{ 0, 0 }, .south_Mg = &.{ 0, 0 }, .north_Mg = &.{ 0, 0 } }, &workspace);
    try std.testing.expectApproxEqAbs(@as(f64, 9), inventory.soil[0].broadcast_ammonium_mol_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), inventory.soil[1].broadcast_ammonium_mol_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 18), inventory.soil[0].banded_nitrate_mol_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2), inventory.soil[1].banded_nitrate_mol_n, 1e-14);
}

test "external fertilizer sediment export retains REDIST N and ion stoichiometry" {
    var inventory = try inventory_module.State.init(std.testing.allocator, 1, 1);
    defer inventory.deinit();
    inventory.soil[0].broadcast_ammonium_mol_n = 10;
    inventory.soil[0].banded_nitrate_mol_n = 20;
    var workspace = try constituents.PackedWorkspace.init(
        std.testing.allocator,
        1,
        component_count,
    );
    defer workspace.deinit();
    try route(1, 1, &.{10}, &inventory, .{
        .east_Mg = &.{1},
        .west_Mg = &.{0},
        .south_Mg = &.{0},
        .north_Mg = &.{0},
    }, &workspace);
    const loss = try exported(&workspace, 14);
    try std.testing.expectApproxEqAbs(@as(f64, 42), loss.nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4), loss.ion_mol, 1e-12);
}
