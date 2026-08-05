const std = @import("std");
const Chemistry = @import("surface_litter_chemistry.zig").State;
const runoff_carrier = @import("surface_runoff_carrier.zig");

pub const Directions = struct {
    east_m3: []const f64,
    west_m3: []const f64,
    south_m3: []const f64,
    north_m3: []const f64,
};

pub const Output = struct {
    inorganic_nitrogen_export_g_n_by_cell: []f64,
    inorganic_phosphorus_export_g_p_by_cell: []f64,
};

const Species = enum(u8) { ammonium, ammonia, nitrate, nitrite, hpo4, h2po4 };
const species_count = @typeInfo(Species).@"enum".fields.len;

/// Translates the TRNSFR surface `VFLW * pool` donor transaction without
/// repeating a full sub-hour model cycle. All directional transfers use the
/// converged hourly runoff, update simultaneously, and commit atomically.
pub fn advance(
    allocator: std.mem.Allocator,
    chemistry: *Chemistry,
    nitrite_g_n: []f64,
    columns: usize,
    rows: usize,
    post_runoff_water_m3: []const f64,
    runoff_water_change_m3: []const f64,
    directions: Directions,
    maximum_transport_fraction: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    phosphorus_molar_mass_g_per_mol: f64,
    output: Output,
) !void {
    const cells = try std.math.mul(usize, columns, rows);
    if (chemistry.cells.len != cells or nitrite_g_n.len != cells or post_runoff_water_m3.len != cells or runoff_water_change_m3.len != cells or output.inorganic_nitrogen_export_g_n_by_cell.len != cells or output.inorganic_phosphorus_export_g_p_by_cell.len != cells) return error.SurfaceMineralTransportDimensionMismatch;
    inline for (.{ directions.east_m3, directions.west_m3, directions.south_m3, directions.north_m3 }) |values| if (values.len != cells) return error.SurfaceMineralTransportDimensionMismatch;
    if (!std.math.isFinite(maximum_transport_fraction) or maximum_transport_fraction < 0 or maximum_transport_fraction > 1 or !std.math.isFinite(nitrogen_molar_mass_g_per_mol) or nitrogen_molar_mass_g_per_mol <= 0 or !std.math.isFinite(phosphorus_molar_mass_g_per_mol) or phosphorus_molar_mass_g_per_mol <= 0) return error.InvalidSurfaceMineralTransportParameter;

    const amount_count = try std.math.mul(usize, cells, species_count);
    const original = try allocator.alloc(f64, amount_count);
    defer allocator.free(original);
    const candidate = try allocator.alloc(f64, amount_count);
    defer allocator.free(candidate);
    const boundary_export_mol = try allocator.alloc(f64, amount_count);
    defer allocator.free(boundary_export_mol);
    const pre_runoff_water = try allocator.alloc(f64, cells);
    defer allocator.free(pre_runoff_water);
    @memset(output.inorganic_nitrogen_export_g_n_by_cell, 0);
    @memset(output.inorganic_phosphorus_export_g_p_by_cell, 0);

    for (0..cells) |cell| {
        const water_after = post_runoff_water_m3[cell];
        const water_before = water_after - runoff_water_change_m3[cell];
        if (!std.math.isFinite(water_after) or water_after < 0 or !std.math.isFinite(water_before) or water_before < -1e-12) return error.InvalidSurfaceMineralWaterState;
        pre_runoff_water[cell] = @max(0, water_before);
        const state = chemistry.cells[cell];
        const concentrations = [_]f64{ state.ammonium_mol_per_m3, state.ammonia_mol_per_m3, state.nitrate_mol_per_m3, 0, state.hpo4_mol_p_per_m3, state.h2po4_mol_p_per_m3 };
        for (concentrations, 0..) |concentration, species| {
            if (!std.math.isFinite(concentration) or concentration < 0) return error.InvalidSurfaceMineralChemistryState;
            original[cell * species_count + species] = if (species == @intFromEnum(Species.nitrite)) blk: {
                if (!std.math.isFinite(nitrite_g_n[cell]) or nitrite_g_n[cell] < 0) return error.InvalidSurfaceMineralChemistryState;
                break :blk nitrite_g_n[cell] / nitrogen_molar_mass_g_per_mol;
            } else concentration * pre_runoff_water[cell];
        }
    }
    try runoff_carrier.calculateChanges(
        columns,
        rows,
        species_count,
        original,
        pre_runoff_water,
        .{
            .east_m3 = directions.east_m3,
            .west_m3 = directions.west_m3,
            .south_m3 = directions.south_m3,
            .north_m3 = directions.north_m3,
        },
        maximum_transport_fraction,
        candidate,
        boundary_export_mol,
    );
    for (candidate, original) |*change, amount| change.* += amount;
    for (0..cells) |cell| for (0..species_count) |species| {
        const exported_mol =
            boundary_export_mol[cell * species_count + species];
        switch (@as(Species, @enumFromInt(species))) {
            .ammonium, .ammonia, .nitrate, .nitrite => output.inorganic_nitrogen_export_g_n_by_cell[cell] +=
                exported_mol * nitrogen_molar_mass_g_per_mol,
            .hpo4, .h2po4 => output.inorganic_phosphorus_export_g_p_by_cell[cell] +=
                exported_mol * phosphorus_molar_mass_g_per_mol,
        }
    };

    for (0..cells) |cell| {
        const water = post_runoff_water_m3[cell];
        for (0..species_count) |species| {
            const amount = candidate[cell * species_count + species];
            // A dry cell holding dissolved mineral mass is a legitimate state once
            // surface evaporation is active, and it is handled below by holding the
            // stored concentrations rather than dividing by zero. Rejecting it here
            // was the fourth blocker on enabling evaporation. See EXEC-004.
            _ = water;
            if (!std.math.isFinite(amount) or amount < -1e-12) return error.InvalidSurfaceMineralTransportCandidate;
        }
    }
    for (0..cells) |cell|
        try chemistry.renormalizeMinerals(cell, post_runoff_water_m3[cell]);
    for (0..cells) |cell| {
        const water = post_runoff_water_m3[cell];
        if (water <= 0) {
            // Dry cell: hold the stored concentrations. Writing zero here would
            // silently destroy the dissolved mineral mass, which is why this case
            // used to be rejected outright. The litter carrier rebase remembers the
            // carrier these concentrations refer to, so the extensive amount stays
            // recoverable on rewetting.
            continue;
        }
        const inverse_water = 1.0 / water;
        chemistry.cells[cell].ammonium_mol_per_m3 = @max(0, candidate[cell * species_count + @intFromEnum(Species.ammonium)]) * inverse_water;
        chemistry.cells[cell].ammonia_mol_per_m3 = @max(0, candidate[cell * species_count + @intFromEnum(Species.ammonia)]) * inverse_water;
        chemistry.cells[cell].nitrate_mol_per_m3 = @max(0, candidate[cell * species_count + @intFromEnum(Species.nitrate)]) * inverse_water;
        nitrite_g_n[cell] = @max(0, candidate[cell * species_count + @intFromEnum(Species.nitrite)]) * nitrogen_molar_mass_g_per_mol;
        chemistry.cells[cell].hpo4_mol_p_per_m3 = @max(0, candidate[cell * species_count + @intFromEnum(Species.hpo4)]) * inverse_water;
        chemistry.cells[cell].h2po4_mol_p_per_m3 = @max(0, candidate[cell * species_count + @intFromEnum(Species.h2po4)]) * inverse_water;
    }
}

test "surface mineral runoff conserves internal transfer and reports external N P" {
    var chemistry = try Chemistry.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    chemistry.cells[0].ammonium_mol_per_m3 = 2;
    chemistry.cells[0].ammonia_mol_per_m3 = 1;
    chemistry.cells[0].nitrate_mol_per_m3 = 3;
    chemistry.cells[0].hpo4_mol_p_per_m3 = 0.5;
    chemistry.cells[0].h2po4_mol_p_per_m3 = 1.5;
    chemistry.cells[1].ammonium_mol_per_m3 = 4;
    chemistry.cells[1].ammonia_mol_per_m3 = 2;
    chemistry.cells[1].nitrate_mol_per_m3 = 6;
    chemistry.cells[1].hpo4_mol_p_per_m3 = 1;
    chemistry.cells[1].h2po4_mol_p_per_m3 = 3;
    var nitrogen_export = [_]f64{ 0, 0 };
    var phosphorus_export = [_]f64{ 0, 0 };
    var nitrite = [_]f64{ 14, 0 };
    const zero = [_]f64{ 0, 0 };
    try advance(std.testing.allocator, &chemistry, &nitrite, 2, 1, &.{ 0.5, 1.25 }, &.{ -0.5, 0.25 }, .{
        .east_m3 = &.{ 0.5, 0.25 },
        .west_m3 = &zero,
        .south_m3 = &zero,
        .north_m3 = &zero,
    }, 1, 14, 31, .{ .inorganic_nitrogen_export_g_n_by_cell = &nitrogen_export, .inorganic_phosphorus_export_g_p_by_cell = &phosphorus_export });
    try std.testing.expectApproxEqAbs(@as(f64, 2), chemistry.cells[0].ammonium_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3.2), chemistry.cells[1].ammonium_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 42), nitrogen_export[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 7), nitrite[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 31), phosphorus_export[1], 1e-14);
}

test "a dry surface cell preserves its dissolved mineral concentrations" {
    // EXEC-004: this used to fail with InvalidSurfaceMineralTransportCandidate,
    // the fourth blocker on enabling surface evaporation. The guard existed for a
    // real reason: the writeback used `inverse_water = 0` for a dry cell, which
    // would have silently zeroed the concentrations and destroyed the mass. The fix
    // holds them instead, so both the error and the mass loss are gone.
    var chemistry = try Chemistry.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    chemistry.cells[0].ammonium_mol_per_m3 = 2;
    chemistry.cells[0].nitrate_mol_per_m3 = 3;
    chemistry.cells[0].h2po4_mol_p_per_m3 = 1.5;
    var nitrogen_export = [_]f64{0};
    var phosphorus_export = [_]f64{0};
    var nitrite = [_]f64{0};
    const zero = [_]f64{0};
    // Evaporate the carrier to exactly dry with no runoff.
    try advance(std.testing.allocator, &chemistry, &nitrite, 1, 1, &.{1}, &.{0}, .{
        .east_m3 = &zero,
        .west_m3 = &zero,
        .south_m3 = &zero,
        .north_m3 = &zero,
    }, 1, 14, 31, .{
        .inorganic_nitrogen_export_g_n_by_cell = &nitrogen_export,
        .inorganic_phosphorus_export_g_p_by_cell = &phosphorus_export,
    });
    // A wet baseline first: with a unit carrier the concentrations are unchanged.
    try std.testing.expectApproxEqAbs(@as(f64, 2), chemistry.cells[0].ammonium_mol_per_m3, 1e-15);

    // Now dry it out. This must not error and must not zero the pools.
    try advance(std.testing.allocator, &chemistry, &nitrite, 1, 1, &.{0}, &.{0}, .{
        .east_m3 = &zero,
        .west_m3 = &zero,
        .south_m3 = &zero,
        .north_m3 = &zero,
    }, 1, 14, 31, .{
        .inorganic_nitrogen_export_g_n_by_cell = &nitrogen_export,
        .inorganic_phosphorus_export_g_p_by_cell = &phosphorus_export,
    });
    try std.testing.expect(chemistry.cells[0].ammonium_mol_per_m3 > 0);
    try std.testing.expect(chemistry.cells[0].nitrate_mol_per_m3 > 0);
    try std.testing.expect(chemistry.cells[0].h2po4_mol_p_per_m3 > 0);
    // Specifically, they are held at their previous values rather than scaled.
    try std.testing.expectApproxEqAbs(@as(f64, 2), chemistry.cells[0].ammonium_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), chemistry.cells[0].nitrate_mol_per_m3, 1e-15);
}

test "failed surface mineral transport leaves chemistry unchanged" {
    var chemistry = try Chemistry.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    chemistry.cells[0].ammonium_mol_per_m3 = 2;
    var n = [_]f64{0};
    var p = [_]f64{0};
    var nitrite = [_]f64{0};
    try std.testing.expectError(error.InvalidSurfaceMineralWaterState, advance(std.testing.allocator, &chemistry, &nitrite, 1, 1, &.{0}, &.{1}, .{ .east_m3 = &.{0}, .west_m3 = &.{0}, .south_m3 = &.{0}, .north_m3 = &.{0} }, 1, 14, 31, .{ .inorganic_nitrogen_export_g_n_by_cell = &n, .inorganic_phosphorus_export_g_p_by_cell = &p }));
    try std.testing.expectEqual(@as(f64, 2), chemistry.cells[0].ammonium_mol_per_m3);
}

test "runoff water change preserves nonmobile solid mineral inventories" {
    var chemistry = try Chemistry.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    chemistry.cells[0].salt_minerals.calcite_mol_per_m3 = 2;
    chemistry.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3 = 3;
    try chemistry.bindMineralReferenceWater(&.{1});
    var n = [_]f64{0};
    var p = [_]f64{0};
    var nitrite = [_]f64{0};
    const zero = [_]f64{0};
    try advance(
        std.testing.allocator,
        &chemistry,
        &nitrite,
        1,
        1,
        &.{0.5},
        &.{-0.5},
        .{
            .east_m3 = &zero,
            .west_m3 = &zero,
            .south_m3 = &zero,
            .north_m3 = &zero,
        },
        1,
        14,
        31,
        .{
            .inorganic_nitrogen_export_g_n_by_cell = &n,
            .inorganic_phosphorus_export_g_p_by_cell = &p,
        },
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        chemistry.cells[0].salt_minerals.calcite_mol_per_m3 * 0.5,
    );
    try std.testing.expectEqual(
        @as(f64, 3),
        chemistry.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3 * 0.5,
    );
}
