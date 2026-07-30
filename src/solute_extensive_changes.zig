const std = @import("std");
const aqueous_network = @import("solute_aqueous_network.zig");
const phosphate_network = @import("solute_phosphate_network.zig");
const cation_exchange = @import("solute_cation_exchange.zig");
const geochemistry = @import("solute_geochemistry_network.zig");
const chemistry = @import("solute_chemistry_state.zig");
const reaction_solver = @import("solute_reaction_solver.zig");

pub const Geometry = struct {
    shared_water_volume_m3: f64,
    ammonium_non_band_water_volume_m3: f64,
    ammonium_band_water_volume_m3: f64,
    phosphate_non_band_water_volume_m3: f64,
    phosphate_band_water_volume_m3: f64,
    shared_soil_mass_Mg: f64,
    ammonium_non_band_soil_mass_Mg: f64,
    ammonium_band_soil_mass_Mg: f64,
    phosphate_non_band_soil_mass_Mg: f64,
    phosphate_band_soil_mass_Mg: f64,
};

pub const Changes = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    component_count: usize,
    /// Component ordering is exactly `solute_chemistry_state.packCell`; every
    /// entry is an extensive molar change for REDIST.
    component_change_mol: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !Changes {
        if (cell_count == 0) return error.ZeroExtensiveChangeCellCount;
        const component_count = chemistry.State.packedComponentCount();
        const total = try std.math.mul(usize, cell_count, component_count);
        const values = try allocator.alloc(f64, total);
        @memset(values, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .component_count = component_count, .component_change_mol = values };
    }

    pub fn deinit(self: *Changes) void {
        self.allocator.free(self.component_change_mol);
        self.* = undefined;
    }

    pub fn cell(self: *Changes, cell_index: usize) ![]f64 {
        if (cell_index >= self.cell_count) return error.ExtensiveChangeCellIndexOutOfBounds;
        const start = cell_index * self.component_count;
        return self.component_change_mol[start .. start + self.component_count];
    }
};

/// Converts accepted before/after chemistry state to extensive molar changes.
///
/// This replaces SOLUTE.F lines 2557-2674. Numerical Newton/Picard iterates
/// are not physical sub-hour cycles: only the converged state difference is
/// published for redistribution.
pub fn calculate(before: []const f64, after: []const f64, geometry: Geometry, output_mol: []f64) !void {
    const count = chemistry.State.packedComponentCount();
    if (before.len != count or after.len != count or output_mol.len != count) return error.ChemistryVectorSizeMismatch;
    try validateGeometry(geometry);
    var cursor: usize = 0;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        const volume = aqueousVolume(field.name, geometry);
        output_mol[cursor] = try finiteChange(before[cursor], after[cursor], volume);
        cursor += 1;
    }
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
        const factor = if (isExchangeSite(field.name)) geometry.phosphate_non_band_soil_mass_Mg else geometry.phosphate_non_band_water_volume_m3;
        output_mol[cursor] = try finiteChange(before[cursor], after[cursor], factor);
        cursor += 1;
    }
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
        const factor = if (isExchangeSite(field.name)) geometry.phosphate_band_soil_mass_Mg else geometry.phosphate_band_water_volume_m3;
        output_mol[cursor] = try finiteChange(before[cursor], after[cursor], factor);
        cursor += 1;
    }
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field| {
        const mass = if (std.mem.eql(u8, field.name, "ammonium_non_band")) geometry.ammonium_non_band_soil_mass_Mg else if (std.mem.eql(u8, field.name, "ammonium_band")) geometry.ammonium_band_soil_mass_Mg else geometry.shared_soil_mass_Mg;
        output_mol[cursor] = try finiteChange(before[cursor], after[cursor], mass);
        cursor += 1;
    }
    output_mol[cursor] = try finiteChange(before[cursor], after[cursor], geometry.shared_soil_mass_Mg);
    cursor += 1;
    inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |_| {
        output_mol[cursor] = try finiteChange(before[cursor], after[cursor], geometry.shared_water_volume_m3);
        cursor += 1;
    }
    output_mol[cursor] = try finiteChange(before[cursor], after[cursor], geometry.shared_water_volume_m3);
}

/// Runs one coupled chemistry solve and records its extensive REDIST changes.
/// The before-state is retained until convergence; solver failure leaves both
/// chemistry and the caller's change record untouched.
pub fn solveCellAndCapture(allocator: std.mem.Allocator, state: *chemistry.State, cell_index: usize, parameters: chemistry.ReactionParameters, solver_options: reaction_solver.Options, geometry: Geometry, changes: *Changes) !reaction_solver.Result {
    if (cell_index >= state.cell_count or cell_index >= changes.cell_count) return error.ExtensiveChangeCellIndexOutOfBounds;
    try validateGeometry(geometry);
    const destination = try changes.cell(cell_index);
    const before = try allocator.alloc(f64, chemistry.State.packedComponentCount());
    defer allocator.free(before);
    const after = try allocator.alloc(f64, chemistry.State.packedComponentCount());
    defer allocator.free(after);
    const staged = try allocator.alloc(f64, destination.len);
    defer allocator.free(staged);
    try state.packCell(cell_index, before);
    const result = try reaction_solver.solveCell(allocator, state, cell_index, parameters, solver_options);
    try state.packCell(cell_index, after);
    try calculate(before, after, geometry, staged);
    @memcpy(destination, staged);
    return result;
}

fn aqueousVolume(comptime name: []const u8, geometry: Geometry) f64 {
    if (std.mem.endsWith(u8, name, "_non_band")) return geometry.ammonium_non_band_water_volume_m3;
    if (std.mem.endsWith(u8, name, "_band")) return geometry.ammonium_band_water_volume_m3;
    return geometry.shared_water_volume_m3;
}

fn isExchangeSite(comptime name: []const u8) bool {
    return std.mem.indexOf(u8, name, "site_mol_per_Mg") != null or std.mem.indexOf(u8, name, "adsorbed_") != null;
}

fn finiteChange(before: f64, after: f64, factor: f64) !f64 {
    if (!std.math.isFinite(before) or !std.math.isFinite(after)) return error.NonFiniteChemistryState;
    const result = (after - before) * factor;
    if (!std.math.isFinite(result)) return error.NonFiniteExtensiveChange;
    return result;
}

fn validateGeometry(geometry: Geometry) !void {
    inline for (@typeInfo(Geometry).@"struct".fields) |field| {
        const value = @field(geometry, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryGeometry;
    }
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "intensive chemistry differences use pool-specific REDIST geometry" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const count = chemistry.State.packedComponentCount();
    const before = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(before);
    const after = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(after);
    const changes = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(changes);
    try state.packCell(0, before);
    state.aqueous[0].calcium = 1;
    state.aqueous[0].ammonium_band = 1;
    state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 = 1;
    state.non_band_phosphate[0].adsorbed_h2po4_mol_p_per_Mg = 1;
    state.cation_exchange_mol_per_Mg[0] = filled(cation_exchange.Cations, 0);
    state.cation_exchange_mol_per_Mg[0].ammonium_band = 1;
    state.geochemistry_solids[0].calcite_solid_mol_per_m3 = 1;
    state.water_mol_per_m3[0] = 1;
    try state.packCell(0, after);
    try calculate(before, after, .{ .shared_water_volume_m3 = 2, .ammonium_non_band_water_volume_m3 = 3, .ammonium_band_water_volume_m3 = 4, .phosphate_non_band_water_volume_m3 = 5, .phosphate_band_water_volume_m3 = 6, .shared_soil_mass_Mg = 7, .ammonium_non_band_soil_mass_Mg = 8, .ammonium_band_soil_mass_Mg = 9, .phosphate_non_band_soil_mass_Mg = 10, .phosphate_band_soil_mass_Mg = 11 }, changes);
    var cursor: usize = 0;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        if (std.mem.eql(u8, field.name, "calcium")) try std.testing.expectApproxEqAbs(@as(f64, 2), changes[cursor], 1e-15);
        if (std.mem.eql(u8, field.name, "ammonium_band")) try std.testing.expectApproxEqAbs(@as(f64, 4), changes[cursor], 1e-15);
        cursor += 1;
    }
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
        if (std.mem.eql(u8, field.name, "dissolved_h2po4_mol_p_per_m3")) try std.testing.expectApproxEqAbs(@as(f64, 5), changes[cursor], 1e-15);
        if (std.mem.eql(u8, field.name, "adsorbed_h2po4_mol_p_per_Mg")) try std.testing.expectApproxEqAbs(@as(f64, 10), changes[cursor], 1e-15);
        cursor += 1;
    }
}

test "packed extensive changes retain carboxyl geochemistry and water order" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const count = chemistry.State.packedComponentCount();
    const before = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(before);
    const after = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(after);
    const changes = try std.testing.allocator.alloc(f64, count);
    defer std.testing.allocator.free(changes);

    try state.packCell(0, before);
    state.carboxyl_bound_hydrogen_mol_per_Mg[0] = 1;
    state.geochemistry_solids[0].calcite_solid_mol_per_m3 = 2;
    state.water_mol_per_m3[0] = 3;
    try state.packCell(0, after);
    try calculate(before, after, .{
        .shared_water_volume_m3 = 5,
        .ammonium_non_band_water_volume_m3 = 5,
        .ammonium_band_water_volume_m3 = 5,
        .phosphate_non_band_water_volume_m3 = 5,
        .phosphate_band_water_volume_m3 = 5,
        .shared_soil_mass_Mg = 7,
        .ammonium_non_band_soil_mass_Mg = 7,
        .ammonium_band_soil_mass_Mg = 7,
        .phosphate_non_band_soil_mass_Mg = 7,
        .phosphate_band_soil_mass_Mg = 7,
    }, changes);

    for (changes, 0..) |change, index| {
        const name = chemistry.State.packedComponentName(index).?;
        if (std.mem.eql(u8, name, "carboxyl_bound_hydrogen_mol_per_Mg")) {
            try std.testing.expectEqual(@as(f64, 7), change);
        } else if (std.mem.eql(u8, name, "geochemistry_solids.calcite_solid_mol_per_m3")) {
            try std.testing.expectEqual(@as(f64, 10), change);
        } else if (std.mem.eql(u8, name, "water_mol_per_m3")) {
            try std.testing.expectEqual(@as(f64, 15), change);
        } else {
            try std.testing.expectEqual(@as(f64, 0), change);
        }
    }
}
