const std = @import("std");
const compute = @import("compute.zig");
const chemistry = @import("surface_litter_chemistry.zig");
const extensive = @import("surface_litter_extensive_changes.zig");
const organic = @import("soil_organic_initialization.zig");
const chemistry_parameters = @import("soil_chemistry_parameters.zig");
const chemistry_initialization = @import("soil_chemistry_initialization.zig");
const cation_exchange = @import("solute_cation_exchange.zig");

pub const Diagnostics = struct {
    allocator: std.mem.Allocator,
    changes: []extensive.Changes,
    iterations: []u16,
    newton_raphson_steps: []u16,
    picard_steps: []u16,
    maximum_scaled_residual: []f64,
    solved: []bool,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !Diagnostics {
        if (cell_count == 0) return error.ZeroLitterChemistryCellCount;
        const changes = try allocator.alloc(extensive.Changes, cell_count);
        errdefer allocator.free(changes);
        const iterations = try allocator.alloc(u16, cell_count);
        errdefer allocator.free(iterations);
        const newton = try allocator.alloc(u16, cell_count);
        errdefer allocator.free(newton);
        const picard = try allocator.alloc(u16, cell_count);
        errdefer allocator.free(picard);
        const residual = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(residual);
        const solved = try allocator.alloc(bool, cell_count);
        var result = Diagnostics{ .allocator = allocator, .changes = changes, .iterations = iterations, .newton_raphson_steps = newton, .picard_steps = picard, .maximum_scaled_residual = residual, .solved = solved };
        result.reset();
        return result;
    }

    pub fn deinit(self: *Diagnostics) void {
        self.allocator.free(self.solved);
        self.allocator.free(self.maximum_scaled_residual);
        self.allocator.free(self.picard_steps);
        self.allocator.free(self.newton_raphson_steps);
        self.allocator.free(self.iterations);
        self.allocator.free(self.changes);
        self.* = undefined;
    }

    pub fn reset(self: *Diagnostics) void {
        @memset(self.changes, std.mem.zeroes(extensive.Changes));
        @memset(self.iterations, 0);
        @memset(self.newton_raphson_steps, 0);
        @memset(self.picard_steps, 0);
        @memset(self.maximum_scaled_residual, 0);
        @memset(self.solved, false);
    }
};

pub const ApplyContext = struct {
    state: *chemistry.State,
    surface_organic: *const organic.State,
    litter_water_m3: []const f64,
    chemistry_parameters: chemistry_parameters.Parameters,
    cation_selectivity_by_cell: []const cation_exchange.Selectivity,
    litter_dry_mass_Mg_per_g_c: f64,
    dynamic_salts: bool,
    solver_options: chemistry.Options,
    diagnostics: *Diagnostics,
};

pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    const count = context.state.cells.len;
    if (context.surface_organic.layer_count != count or context.litter_water_m3.len != count or context.cation_selectivity_by_cell.len != count or context.diagnostics.changes.len != count or range.first > range.end or range.end > count) return error.LitterChemistryStepDimensionMismatch;
    if (!std.math.isFinite(context.litter_dry_mass_Mg_per_g_c) or context.litter_dry_mass_Mg_per_g_c <= 0) return error.InvalidLitterDryMassCoefficient;
    for (range.first..range.end) |cell| try applyCell(context, cell);
}

fn applyCell(context: *ApplyContext, cell: usize) !void {
    const carbon_g_c = try context.surface_organic.totalCarbon_g_c(cell);
    const water_m3 = context.litter_water_m3[cell];
    if (!std.math.isFinite(water_m3) or water_m3 < 0) return error.InvalidLitterWaterVolume;
    try context.state.renormalizeMinerals(cell, water_m3);
    if (carbon_g_c == 0 or water_m3 == 0) return;
    const dry_mass_Mg = context.litter_dry_mass_Mg_per_g_c * carbon_g_c;
    if (!std.math.isFinite(dry_mass_Mg) or dry_mass_Mg <= 0) return error.InvalidLitterDryMass;
    const density = dry_mass_Mg / water_m3;
    const capacity = try chemistry_initialization.surfaceLitterCationExchangeCapacity_mol_charge_per_Mg_litter(carbon_g_c, dry_mass_Mg, context.chemistry_parameters.surface_litter.carboxyl_sites_mol_per_Mg_c);
    const before = context.state.cells[cell];
    const activity = try chemistry.activityCoefficients(before, water_m3);
    const rates_context = context.chemistry_parameters.forSurfaceLitter(activity, context.cation_selectivity_by_cell[cell], capacity, density, context.dynamic_salts);
    // STARTE initializes occupied carboxyl sites from H activity before its
    // first equilibrium cycle. Preserve a nonzero restored value thereafter.
    if (context.state.cells[cell].carboxyl_hydrogen_mol_per_Mg == 0 and capacity > 0) {
        const hydrogen_activity = before.hydrogen_mol_per_m3 * activity.monovalent_activity_coefficient;
        context.state.cells[cell].carboxyl_hydrogen_mol_per_Mg = capacity * @min(1.0, hydrogen_activity / context.chemistry_parameters.surface_litter.carboxyl_dissociation_constant);
    }
    const result = extensive.solveCellAndCapture(context.state, cell, &rates_context, context.solver_options, .{ .litter_water_volume_m3 = water_m3, .litter_dry_mass_Mg = dry_mass_Mg }, &context.diagnostics.changes[cell]) catch |err| {
        context.state.cells[cell] = before;
        return err;
    };
    context.diagnostics.iterations[cell] = result.iterations;
    context.diagnostics.newton_raphson_steps[cell] = result.newton_raphson_steps;
    context.diagnostics.picard_steps[cell] = result.picard_steps;
    context.diagnostics.maximum_scaled_residual[cell] = result.maximum_scaled_residual;
    context.diagnostics.solved[cell] = true;
}

test "dry litter cells skip without mutating chemistry" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    var surface = try organic.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var diagnostics = try Diagnostics.init(std.testing.allocator, 1);
    defer diagnostics.deinit();
    var context: ApplyContext = undefined;
    context.state = &state;
    context.surface_organic = &surface;
    context.litter_water_m3 = &.{0};
    context.cation_selectivity_by_cell = &.{.{ .calcium_ammonium = 1, .calcium_hydrogen = 1, .calcium_aluminum_and_iron = 1, .calcium_magnesium = 1, .calcium_sodium = 1, .calcium_potassium = 1 }};
    context.litter_dry_mass_Mg_per_g_c = 1.82e-6;
    context.dynamic_salts = false;
    context.solver_options = .{};
    context.diagnostics = &diagnostics;
    // Unused because zero carbon/water exits before parameter evaluation.
    context.chemistry_parameters = undefined;
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(!diagnostics.solved[0]);
    try std.testing.expectEqual(std.mem.zeroes(chemistry.Cell), state.cells[0]);
}

test "wet litter tile converges synthetic equilibrium conservatively" {
    const source =
        "aqueous_constants 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n" ++
        "aqueous_kinetics 0.2 0.2 0.1 0.1\n" ++
        "phosphate_constants 1 1 1 1 1 1 1 1\n" ++
        "phosphate_surface 1 1 1 1 1 1 0.1 0.2\n" ++
        "phosphate_minerals 1 1 1 1 1 1 0.1 0.1 0.1\n" ++
        "phosphate_kinetics 0.2 0.1\n" ++
        "cation_kinetics 0.2 0.1\n" ++
        "geochemistry_products 1 1 1 1 1 1 1 1 1 1\n" ++
        "geochemistry_kinetics 0.2 0.2 0.1 0.1 1 0 0\n" ++
        "water_equilibrium 1 1e-30 55555.555555555555\n" ++
        "surface_litter 1 0\n" ++
        "surface_fertilizer 1 1 1 0.05 50 0.03 0.05 0.01 0.005\n";
    const parameters = try chemistry_parameters.parse(source);
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].hydrogen_mol_per_m3 = 1;
    state.cells[0].hydroxide_mol_per_m3 = 1;
    state.cells[0].aluminum_mol_per_m3 = 1;
    state.cells[0].iron_mol_per_m3 = 1;
    state.cells[0].calcium_mol_per_m3 = 1;
    state.cells[0].magnesium_mol_per_m3 = 1;
    state.cells[0].sodium_mol_per_m3 = 1;
    state.cells[0].potassium_mol_per_m3 = 1;
    state.cells[0].hpo4_mol_p_per_m3 = 1;
    state.cells[0].h2po4_mol_p_per_m3 = 1;
    state.cells[0].phosphate_minerals = .{ .aluminum_phosphate_mol_per_m3 = 1, .iron_phosphate_mol_per_m3 = 1, .dicalcium_phosphate_mol_per_m3 = 1, .hydroxyapatite_mol_per_m3 = 1, .monocalcium_phosphate_mol_per_m3 = 1 };
    try state.bindMineralReferenceWater(&.{1});
    var surface = try organic.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    surface.structural[0].carbon_g_c = 1;
    const phosphorus_before = state.cells[0].hpo4_mol_p_per_m3 +
        state.cells[0].h2po4_mol_p_per_m3 +
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3 +
        state.cells[0].phosphate_minerals.iron_phosphate_mol_per_m3 +
        state.cells[0].phosphate_minerals.dicalcium_phosphate_mol_per_m3 +
        3 * state.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3 +
        2 * state.cells[0].phosphate_minerals.monocalcium_phosphate_mol_per_m3;
    var diagnostics = try Diagnostics.init(std.testing.allocator, 1);
    defer diagnostics.deinit();
    var context = ApplyContext{
        .state = &state,
        .surface_organic = &surface,
        .litter_water_m3 = &.{1},
        .chemistry_parameters = parameters,
        .cation_selectivity_by_cell = &.{.{ .calcium_ammonium = 1, .calcium_hydrogen = 1, .calcium_aluminum_and_iron = 1, .calcium_magnesium = 1, .calcium_sodium = 1, .calcium_potassium = 1 }},
        .litter_dry_mass_Mg_per_g_c = 1,
        .dynamic_salts = false,
        .solver_options = .{ .absolute_tolerance = 1e-10, .relative_tolerance = 1e-8, .picard_relaxation = 0.5, .max_iterations = 60 },
        .diagnostics = &diagnostics,
    };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    const phosphorus_after = state.cells[0].hpo4_mol_p_per_m3 +
        state.cells[0].h2po4_mol_p_per_m3 +
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3 +
        state.cells[0].phosphate_minerals.iron_phosphate_mol_per_m3 +
        state.cells[0].phosphate_minerals.dicalcium_phosphate_mol_per_m3 +
        3 * state.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3 +
        2 * state.cells[0].phosphate_minerals.monocalcium_phosphate_mol_per_m3;
    try std.testing.expectApproxEqAbs(phosphorus_before, phosphorus_after, 1e-12);
    try std.testing.expect(diagnostics.solved[0]);
}
