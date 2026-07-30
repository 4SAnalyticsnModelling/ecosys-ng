const std = @import("std");
const chemistry_module = @import("solute_chemistry_state.zig");
const phosphate = @import("solute_phosphate_network.zig");
const cation = @import("solute_cation_exchange.zig");
const geochemistry = @import("solute_geochemistry_network.zig");
const constituents = @import("eroded_constituents.zig");

pub const component_count: usize = @typeInfo(cation.Cations).@"struct".fields.len +
    2 * phosphateErodibleFieldCount() +
    @typeInfo(geochemistry.SolidState).@"struct".fields.len + 1;

pub const Exported = struct {
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
    inorganic_carbon_g_c: f64 = 0,
    ion_mol: f64 = 0,
};

pub fn route(
    columns: usize,
    rows: usize,
    soil_layer_capacity: usize,
    surface_soil_mass_Mg: []const f64,
    topsoil_water_volume_m3: []const f64,
    chemistry: *chemistry_module.State,
    sediment: constituents.DirectionalSediment,
    workspace: *constituents.PackedWorkspace,
) !void {
    const cells = try std.math.mul(usize, columns, rows);
    if (soil_layer_capacity == 0 or chemistry.cell_count != try std.math.mul(usize, cells, soil_layer_capacity) or surface_soil_mass_Mg.len != cells or topsoil_water_volume_m3.len != chemistry.cell_count or workspace.cell_count != cells or workspace.component_count != component_count) return error.ChemistryErosionDimensionMismatch;
    try pack(cells, soil_layer_capacity, surface_soil_mass_Mg, topsoil_water_volume_m3, chemistry, workspace.pools);
    try constituents.routePackedWorkspace(workspace, columns, rows, surface_soil_mass_Mg, sediment);
    try unpack(cells, soil_layer_capacity, surface_soil_mass_Mg, topsoil_water_volume_m3, chemistry, workspace.pools);
}

/// REDIST `ZXE/PXE/PPE/CXE/SEX/SEP` external solid-chemistry loss.
pub fn exported(
    workspace: *const constituents.PackedWorkspace,
    carbon_g_per_mol: f64,
    nitrogen_g_per_mol: f64,
    phosphorus_g_per_mol: f64,
) !Exported {
    if (workspace.component_count != component_count or
        workspace.exported.len !=
            try std.math.mul(usize, workspace.cell_count, component_count))
        return error.ChemistryErosionDimensionMismatch;
    inline for (.{ carbon_g_per_mol, nitrogen_g_per_mol, phosphorus_g_per_mol }) |mass|
        if (!std.math.isFinite(mass) or mass <= 0)
            return error.InvalidChemistryErosionMolarMass;
    var result: Exported = .{};
    for (0..workspace.cell_count) |cell| {
        var cursor = cell * component_count;
        inline for (@typeInfo(cation.Cations).@"struct".fields) |field| {
            const amount = try exportedAmount(workspace.exported[cursor]);
            cursor += 1;
            if (comptime std.mem.startsWith(u8, field.name, "ammonium"))
                result.nitrogen_g_n += amount * nitrogen_g_per_mol;
            result.ion_mol += amount *
                if (comptime std.mem.startsWith(u8, field.name, "ammonium"))
                    2
                else
                    1;
        }
        // Carboxyl-bound hydrogen.
        result.ion_mol += try exportedAmount(workspace.exported[cursor]);
        cursor += 1;
        inline for (0..2) |_| inline for (@typeInfo(phosphate.State).@"struct".fields) |field| {
            if (comptime isErodiblePhosphateField(field.name)) {
                const amount = try exportedAmount(workspace.exported[cursor]);
                cursor += 1;
                result.phosphorus_g_p +=
                    amount * phosphateAtoms(field.name) * phosphorus_g_per_mol;
                result.ion_mol += amount * phosphateIonAtoms(field.name);
            }
        };
        inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field| {
            const amount = try exportedAmount(workspace.exported[cursor]);
            cursor += 1;
            result.ion_mol += amount * geochemistryIonAtoms(field.name);
            if (comptime std.mem.eql(u8, field.name, "calcite_solid_mol_per_m3"))
                result.inorganic_carbon_g_c += amount * carbon_g_per_mol;
        }
        if (cursor != (cell + 1) * component_count)
            return error.ChemistryErosionDimensionMismatch;
    }
    inline for (std.meta.fields(Exported)) |field|
        if (!std.math.isFinite(@field(result, field.name)) or
            @field(result, field.name) < 0)
            return error.ChemistryErosionExportOverflow;
    return result;
}

fn exportedAmount(value: f64) !f64 {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidChemistryErosionExport;
    return value;
}

fn phosphateAtoms(comptime name: []const u8) f64 {
    if (std.mem.indexOf(u8, name, "adsorbed_") != null) return 1;
    if (std.mem.indexOf(u8, name, "aluminum_phosphate") != null or
        std.mem.indexOf(u8, name, "iron_phosphate") != null or
        std.mem.indexOf(u8, name, "dicalcium_phosphate") != null)
        return 1;
    if (std.mem.indexOf(u8, name, "hydroxyapatite") != null) return 3;
    if (std.mem.indexOf(u8, name, "monocalcium_phosphate") != null) return 2;
    return 0;
}

fn phosphateIonAtoms(comptime name: []const u8) f64 {
    if (std.mem.indexOf(u8, name, "deprotonated_site") != null) return 1;
    if (std.mem.indexOf(u8, name, "hydroxyl_site") != null) return 2;
    if (std.mem.indexOf(u8, name, "protonated_site") != null) return 3;
    if (std.mem.indexOf(u8, name, "adsorbed_hpo4") != null) return 3;
    if (std.mem.indexOf(u8, name, "adsorbed_h2po4") != null) return 4;
    if (std.mem.indexOf(u8, name, "aluminum_phosphate") != null or
        std.mem.indexOf(u8, name, "iron_phosphate") != null)
        return 2;
    if (std.mem.indexOf(u8, name, "dicalcium_phosphate") != null) return 3;
    if (std.mem.indexOf(u8, name, "hydroxyapatite") != null) return 9;
    if (std.mem.indexOf(u8, name, "monocalcium_phosphate") != null) return 7;
    return 0;
}

fn geochemistryIonAtoms(comptime name: []const u8) f64 {
    if (std.mem.indexOf(u8, name, "gibbsite") != null or
        std.mem.indexOf(u8, name, "iron_hydroxide") != null)
        return 4;
    if (std.mem.indexOf(u8, name, "calcite") != null or
        std.mem.indexOf(u8, name, "gypsum") != null)
        return 2;
    return 1;
}

fn pack(cells: usize, layer_capacity: usize, soil_mass_Mg: []const f64, water_m3: []const f64, chemistry: *const chemistry_module.State, output: []f64) !void {
    for (0..cells) |cell| {
        const layer = cell * layer_capacity;
        const soil_mass = soil_mass_Mg[cell];
        const water = water_m3[layer];
        if (!std.math.isFinite(soil_mass) or soil_mass <= 0 or !std.math.isFinite(water) or water < 0) return error.InvalidChemistryErosionState;
        var cursor = cell * component_count;
        inline for (@typeInfo(cation.Cations).@"struct".fields) |field| {
            try put(@field(chemistry.cation_exchange_mol_per_Mg[layer], field.name) * soil_mass, output, &cursor);
        }
        try put(chemistry.carboxyl_bound_hydrogen_mol_per_Mg[layer] * soil_mass, output, &cursor);
        inline for (.{ chemistry.non_band_phosphate[layer], chemistry.band_phosphate[layer] }) |zone| inline for (@typeInfo(phosphate.State).@"struct".fields) |field| {
            if (comptime isErodiblePhosphateField(field.name)) {
                const scale = if (comptime std.mem.endsWith(u8, field.name, "_per_Mg")) soil_mass else water;
                try put(@field(zone, field.name) * scale, output, &cursor);
            }
        };
        inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field| {
            try put(@field(chemistry.geochemistry_solids[layer], field.name) * water, output, &cursor);
        }
        if (cursor != (cell + 1) * component_count) return error.ChemistryErosionDimensionMismatch;
    }
}

fn unpack(cells: usize, layer_capacity: usize, soil_mass_Mg: []const f64, water_m3: []const f64, chemistry: *chemistry_module.State, input: []const f64) !void {
    for (input) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryErosionCandidate;
    for (0..cells) |cell| {
        const layer = cell * layer_capacity;
        const soil_mass = soil_mass_Mg[cell];
        const water = water_m3[layer];
        var cursor = cell * component_count;
        inline for (@typeInfo(cation.Cations).@"struct".fields) |field| {
            @field(chemistry.cation_exchange_mol_per_Mg[layer], field.name) = input[cursor] / soil_mass;
            cursor += 1;
        }
        chemistry.carboxyl_bound_hydrogen_mol_per_Mg[layer] = input[cursor] / soil_mass;
        cursor += 1;
        inline for (.{ &chemistry.non_band_phosphate[layer], &chemistry.band_phosphate[layer] }) |zone| inline for (@typeInfo(phosphate.State).@"struct".fields) |field| {
            if (comptime isErodiblePhosphateField(field.name)) {
                const scale = if (comptime std.mem.endsWith(u8, field.name, "_per_Mg")) soil_mass else water;
                @field(zone, field.name) = try concentration(input[cursor], scale);
                cursor += 1;
            }
        };
        inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field| {
            @field(chemistry.geochemistry_solids[layer], field.name) = try concentration(input[cursor], water);
            cursor += 1;
        }
    }
}

fn put(value: f64, output: []f64, cursor: *usize) !void {
    if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryErosionState;
    output[cursor.*] = value;
    cursor.* += 1;
}

fn concentration(amount: f64, scale: f64) !f64 {
    if (scale > 0) return amount / scale;
    if (amount == 0) return 0;
    return error.ErodedSolidRequiresRecipientVolume;
}

fn isErodiblePhosphateField(comptime name: []const u8) bool {
    @setEvalBranchQuota(10_000);
    return std.mem.endsWith(u8, name, "_per_Mg") or std.mem.indexOf(u8, name, "_solid_mol_per_m3") != null;
}

fn phosphateErodibleFieldCount() usize {
    comptime var count: usize = 0;
    inline for (@typeInfo(phosphate.State).@"struct".fields) |field| if (isErodiblePhosphateField(field.name)) {
        count += 1;
    };
    return count;
}

test "solid chemistry amounts follow sediment while retaining native concentrations" {
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    chemistry.cation_exchange_mol_per_Mg[0].calcium = 2;
    chemistry.carboxyl_bound_hydrogen_mol_per_Mg[0] = 4;
    chemistry.non_band_phosphate[0].adsorbed_hpo4_mol_p_per_Mg = 3;
    chemistry.non_band_phosphate[0].aluminum_phosphate_solid_mol_per_m3 = 10;
    chemistry.geochemistry_solids[0].gibbsite_solid_mol_per_m3 = 20;
    var workspace = try constituents.PackedWorkspace.init(std.testing.allocator, 2, component_count);
    defer workspace.deinit();
    try route(2, 1, 1, &.{ 10, 10 }, &.{ 1, 1 }, &chemistry, .{ .east_Mg = &.{ 1, 0 }, .west_Mg = &.{ 0, 0 }, .south_Mg = &.{ 0, 0 }, .north_Mg = &.{ 0, 0 } }, &workspace);
    try std.testing.expectApproxEqAbs(@as(f64, 1.8), chemistry.cation_exchange_mol_per_Mg[0].calcium, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), chemistry.cation_exchange_mol_per_Mg[1].calcium, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 3.6), chemistry.carboxyl_bound_hydrogen_mol_per_Mg[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), chemistry.carboxyl_bound_hydrogen_mol_per_Mg[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.7), chemistry.non_band_phosphate[0].adsorbed_hpo4_mol_p_per_Mg, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), chemistry.non_band_phosphate[1].adsorbed_hpo4_mol_p_per_Mg, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 9), chemistry.non_band_phosphate[0].aluminum_phosphate_solid_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), chemistry.non_band_phosphate[1].aluminum_phosphate_solid_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 18), chemistry.geochemistry_solids[0].gibbsite_solid_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2), chemistry.geochemistry_solids[1].gibbsite_solid_mol_per_m3, 1e-14);
}

test "external solid chemistry export reproduces REDIST C N P ion counts" {
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    chemistry.cation_exchange_mol_per_Mg[0].ammonium_non_band = 10;
    chemistry.cation_exchange_mol_per_Mg[0].calcium = 20;
    chemistry.carboxyl_bound_hydrogen_mol_per_Mg[0] = 30;
    chemistry.non_band_phosphate[0].adsorbed_hpo4_mol_p_per_Mg = 40;
    chemistry.non_band_phosphate[0].hydroxyapatite_solid_mol_per_m3 = 50;
    chemistry.geochemistry_solids[0].calcite_solid_mol_per_m3 = 60;
    var workspace = try constituents.PackedWorkspace.init(
        std.testing.allocator,
        1,
        component_count,
    );
    defer workspace.deinit();
    try route(1, 1, 1, &.{10}, &.{1}, &chemistry, .{
        .east_Mg = &.{1},
        .west_Mg = &.{0},
        .south_Mg = &.{0},
        .north_Mg = &.{0},
    }, &workspace);
    const loss = try exported(&workspace, 12, 14, 31);
    try std.testing.expectApproxEqAbs(@as(f64, 140), loss.nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(
        @as(f64, 31 * (40 + 15)),
        loss.phosphorus_g_p,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 72), loss.inorganic_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(
        @as(f64, 20 + 20 + 30 + 120 + 45 + 12),
        loss.ion_mol,
        1e-12,
    );
}
