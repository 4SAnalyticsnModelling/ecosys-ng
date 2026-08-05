const std = @import("std");
const litter_chemistry = @import("surface_litter_chemistry.zig");

pub const Inventory = struct {
    ammonium_mol_n: f64,
    ammonia_mol_n: f64,
    urea_mol_n: f64,
    nitrate_mol_n: f64,
    initial_urease_inhibition_fraction: f64,
    current_urease_inhibition_fraction: f64,
};

pub const Environment = struct {
    litter_water_volume_m3: f64,
    litter_water_volume_fraction: f64,
    litter_dry_mass_megagrams: f64,
    biologically_active_water_volume_m3: f64,
    active_biomass_respiration_g_c_per_step: f64,
    microbial_temperature_factor: f64,
    step_duration_h: f64,
};

pub const Parameters = struct {
    ammonium_dissolution_fraction_per_step: f64,
    ammonia_dissolution_fraction_per_step: f64,
    nitrate_dissolution_fraction_per_step: f64,
    minimum_urea_half_saturation_mol_n_per_megagram: f64,
    microbial_activity_inhibition_g_c_per_m3_per_h: f64,
    specific_urea_hydrolysis_mol_n_per_g_c: f64,
    urease_inhibition_decline_fraction_per_step: f64,
};

pub const Fluxes = struct {
    ammonium_dissolution_mol_n: f64,
    ammonia_dissolution_mol_n: f64,
    urea_hydrolysis_mol_n: f64,
    nitrate_dissolution_mol_n: f64,
};

/// Runtime-sized surface-litter fertilizer inventories.  There is one entry
/// per horizontal model cell; no maximum grid extent is compiled in.
pub const State = struct {
    allocator: std.mem.Allocator,
    cells: []Inventory,
    formulation: []u8,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroLitterCellCount;
        const cells = try allocator.alloc(Inventory, cell_count);
        errdefer allocator.free(cells);
        const formulation = try allocator.alloc(u8, cell_count);
        @memset(cells, .{
            .ammonium_mol_n = 0,
            .ammonia_mol_n = 0,
            .urea_mol_n = 0,
            .nitrate_mol_n = 0,
            .initial_urease_inhibition_fraction = 0,
            .current_urease_inhibition_fraction = 0,
        });
        @memset(formulation, 0);
        return .{ .allocator = allocator, .cells = cells, .formulation = formulation };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.cells);
        self.allocator.free(self.formulation);
        self.* = undefined;
    }

    pub fn calculateAndCommit(self: *State, cell_index: usize, environment: Environment, parameters: Parameters) !Fluxes {
        if (cell_index >= self.cells.len) return error.LitterCellIndexOutOfBounds;
        try validate(self.cells[cell_index], environment, parameters);
        var staged = self.cells[cell_index];

        if (staged.initial_urease_inhibition_fraction > 0 and staged.current_urease_inhibition_fraction > 0) {
            const recovery = @max(parameters.urease_inhibition_decline_fraction_per_step, 1.0 - staged.current_urease_inhibition_fraction / staged.initial_urease_inhibition_fraction);
            staged.current_urease_inhibition_fraction = @max(0, staged.current_urease_inhibition_fraction - parameters.urease_inhibition_decline_fraction_per_step * staged.current_urease_inhibition_fraction * recovery);
        } else staged.current_urease_inhibition_fraction = 0;

        const activity_g_c_per_m3_per_h = if (environment.biologically_active_water_volume_m3 > 0)
            @min(100_000.0, environment.active_biomass_respiration_g_c_per_step / (environment.biologically_active_water_volume_m3 * environment.step_duration_h))
        else
            100_000.0;
        const effective_half_saturation = parameters.minimum_urea_half_saturation_mol_n_per_megagram *
            (1.0 + activity_g_c_per_m3_per_h / parameters.microbial_activity_inhibition_g_c_per_m3_per_h);
        const urea_concentration = if (environment.litter_dry_mass_megagrams > 0) staged.urea_mol_n / environment.litter_dry_mass_megagrams else 0;
        const substrate_factor = if (urea_concentration > 0) urea_concentration / (urea_concentration + effective_half_saturation) else 0;
        const hydrolysis = @min(staged.urea_mol_n, parameters.specific_urea_hydrolysis_mol_n_per_g_c * environment.active_biomass_respiration_g_c_per_step * substrate_factor * environment.microbial_temperature_factor * (1.0 - staged.current_urease_inhibition_fraction));
        const fluxes = Fluxes{
            .ammonium_dissolution_mol_n = parameters.ammonium_dissolution_fraction_per_step * staged.ammonium_mol_n * environment.litter_water_volume_fraction,
            .ammonia_dissolution_mol_n = parameters.ammonia_dissolution_fraction_per_step * staged.ammonia_mol_n,
            .urea_hydrolysis_mol_n = hydrolysis,
            .nitrate_dissolution_mol_n = parameters.nitrate_dissolution_fraction_per_step * staged.nitrate_mol_n * environment.litter_water_volume_fraction,
        };
        staged.ammonium_mol_n -= fluxes.ammonium_dissolution_mol_n;
        staged.ammonia_mol_n -= fluxes.ammonia_dissolution_mol_n;
        staged.urea_mol_n -= fluxes.urea_hydrolysis_mol_n;
        staged.nitrate_mol_n -= fluxes.nitrate_dissolution_mol_n;
        try validateInventory(staged);
        self.cells[cell_index] = staged;
        return fluxes;
    }
};

/// Commits fertilizer loss and its aqueous litter gain as one nitrogen-
/// conserving transaction. Urea hydrolysis enters NH3, matching SOLUTE.F
/// `TRN3S += RSN3AA + RSNUAA`.
pub fn calculateAndApplyToChemistry(state: *State, chemistry: *litter_chemistry.State, cell_index: usize, environment: Environment, parameters: Parameters) !Fluxes {
    if (cell_index >= state.cells.len or cell_index >= chemistry.cells.len) return error.LitterCellIndexOutOfBounds;
    const fertilizer_before = state.cells[cell_index];
    const chemistry_before = chemistry.cells[cell_index];
    const before_nitrogen_mol = inventoryNitrogen(fertilizer_before) + aqueousNitrogen(chemistry_before) * environment.litter_water_volume_m3;
    const fluxes = state.calculateAndCommit(cell_index, environment, parameters) catch |err| return err;
    if (environment.litter_water_volume_m3 <= 0 and (fluxes.ammonium_dissolution_mol_n > 0 or fluxes.ammonia_dissolution_mol_n > 0 or fluxes.urea_hydrolysis_mol_n > 0 or fluxes.nitrate_dissolution_mol_n > 0)) {
        state.cells[cell_index] = fertilizer_before;
        return error.FertilizerDissolutionRequiresLitterWater;
    }
    if (environment.litter_water_volume_m3 > 0) {
        const inverse_water = 1.0 / environment.litter_water_volume_m3;
        const cell = &chemistry.cells[cell_index];
        cell.ammonium_mol_per_m3 += fluxes.ammonium_dissolution_mol_n * inverse_water;
        cell.ammonia_mol_per_m3 += (fluxes.ammonia_dissolution_mol_n + fluxes.urea_hydrolysis_mol_n) * inverse_water;
        cell.nitrate_mol_per_m3 += fluxes.nitrate_dissolution_mol_n * inverse_water;
    }
    const after_nitrogen_mol = inventoryNitrogen(state.cells[cell_index]) + aqueousNitrogen(chemistry.cells[cell_index]) * environment.litter_water_volume_m3;
    const tolerance = 128 * std.math.floatEps(f64) * @max(1, @abs(before_nitrogen_mol));
    if (!std.math.isFinite(after_nitrogen_mol) or @abs(after_nitrogen_mol - before_nitrogen_mol) > tolerance) {
        state.cells[cell_index] = fertilizer_before;
        chemistry.cells[cell_index] = chemistry_before;
        return error.SurfaceFertilizerNitrogenImbalance;
    }
    return fluxes;
}

fn inventoryNitrogen(inventory: Inventory) f64 {
    return inventory.ammonium_mol_n + inventory.ammonia_mol_n + inventory.urea_mol_n + inventory.nitrate_mol_n;
}

fn aqueousNitrogen(cell: litter_chemistry.Cell) f64 {
    return cell.ammonium_mol_per_m3 + cell.ammonia_mol_per_m3 + cell.nitrate_mol_per_m3;
}

fn validate(inventory: Inventory, environment: Environment, parameters: Parameters) !void {
    try validateInventory(inventory);
    inline for (@typeInfo(Environment).@"struct".fields) |field|
        if (!std.math.isFinite(@field(environment, field.name)) or @field(environment, field.name) < 0) return error.InvalidLitterEnvironment;
    inline for (@typeInfo(Parameters).@"struct".fields) |field|
        if (!std.math.isFinite(@field(parameters, field.name)) or @field(parameters, field.name) < 0) return error.InvalidLitterFertilizerParameter;
    if (environment.step_duration_h <= 0 or environment.litter_water_volume_fraction > 1 or environment.microbial_temperature_factor > 1 or
        parameters.microbial_activity_inhibition_g_c_per_m3_per_h <= 0 or parameters.ammonium_dissolution_fraction_per_step > 1 or
        parameters.ammonia_dissolution_fraction_per_step > 1 or parameters.nitrate_dissolution_fraction_per_step > 1 or
        parameters.urease_inhibition_decline_fraction_per_step > 1) return error.InvalidLitterFertilizerParameter;
}

fn validateInventory(inventory: Inventory) !void {
    inline for (@typeInfo(Inventory).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inventory, field.name)) or @field(inventory, field.name) < -1e-12) return error.InvalidLitterFertilizerInventory;
    if (inventory.initial_urease_inhibition_fraction > 1 or inventory.current_urease_inhibition_fraction > 1) return error.InvalidLitterFertilizerInventory;
}

test "zero biologically active litter water saturates COQCK and preserves hydrolysis flow" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0] = .{ .ammonium_mol_n = 0, .ammonia_mol_n = 10, .urea_mol_n = 4, .nitrate_mol_n = 0, .initial_urease_inhibition_fraction = 0, .current_urease_inhibition_fraction = 0 };
    const fluxes = try state.calculateAndCommit(
        0,
        .{ .litter_water_volume_m3 = 2, .litter_water_volume_fraction = 0.5, .litter_dry_mass_megagrams = 2, .biologically_active_water_volume_m3 = 0, .active_biomass_respiration_g_c_per_step = 10, .microbial_temperature_factor = 1, .step_duration_h = 1 },
        .{
            .ammonium_dissolution_fraction_per_step = 0.2,
            .ammonia_dissolution_fraction_per_step = 0.25,
            .nitrate_dissolution_fraction_per_step = 0.5,
            .minimum_urea_half_saturation_mol_n_per_megagram = 1,
            .microbial_activity_inhibition_g_c_per_m3_per_h = 2,
            .specific_urea_hydrolysis_mol_n_per_g_c = 1,
            .urease_inhibition_decline_fraction_per_step = 0.1,
        },
    );
    const expected_urea_concentration = 4.0 / 2;
    const expected_activity = 100_000.0;
    const expected_half_saturation = 1.0 * (1 + expected_activity / 2);
    const expected_fraction = expected_urea_concentration / (expected_urea_concentration + expected_half_saturation);
    const expected_hydrolysis = @min(4.0, 10.0 * expected_fraction * 1 * (1 - 0));
    try std.testing.expect(fluxes.urea_hydrolysis_mol_n > 0);
    try std.testing.expectEqual(@as(f64, 4), state.cells[0].urea_mol_n + fluxes.urea_hydrolysis_mol_n);
    try std.testing.expectApproxEqAbs(expected_hydrolysis, fluxes.urea_hydrolysis_mol_n, 1e-12);
}

test "surface litter fertilizer fluxes commit atomically" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.cells[1] = .{ .ammonium_mol_n = 10, .ammonia_mol_n = 8, .urea_mol_n = 6, .nitrate_mol_n = 4, .initial_urease_inhibition_fraction = 0, .current_urease_inhibition_fraction = 0 };
    const fluxes = try state.calculateAndCommit(1, .{ .litter_water_volume_m3 = 2, .litter_water_volume_fraction = 0.5, .litter_dry_mass_megagrams = 2, .biologically_active_water_volume_m3 = 1, .active_biomass_respiration_g_c_per_step = 2, .microbial_temperature_factor = 1, .step_duration_h = 1 }, .{ .ammonium_dissolution_fraction_per_step = 0.2, .ammonia_dissolution_fraction_per_step = 0.25, .nitrate_dissolution_fraction_per_step = 0.5, .minimum_urea_half_saturation_mol_n_per_megagram = 1, .microbial_activity_inhibition_g_c_per_m3_per_h = 2, .specific_urea_hydrolysis_mol_n_per_g_c = 1, .urease_inhibition_decline_fraction_per_step = 0.1 });
    try std.testing.expectEqual(@as(f64, 1), fluxes.ammonium_dissolution_mol_n);
    try std.testing.expectEqual(@as(f64, 2), fluxes.ammonia_dissolution_mol_n);
    try std.testing.expectEqual(@as(f64, 1), fluxes.nitrate_dissolution_mol_n);
    try std.testing.expect(state.cells[1].urea_mol_n < 6);
    try std.testing.expectEqual(@as(f64, 0), state.cells[0].ammonium_mol_n);
}

test "invalid litter input leaves inventory unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].ammonium_mol_n = 2;
    const before = state.cells[0];
    try std.testing.expectError(error.InvalidLitterEnvironment, state.calculateAndCommit(0, .{ .litter_water_volume_m3 = std.math.nan(f64), .litter_water_volume_fraction = 0.5, .litter_dry_mass_megagrams = 1, .biologically_active_water_volume_m3 = 1, .active_biomass_respiration_g_c_per_step = 1, .microbial_temperature_factor = 1, .step_duration_h = 1 }, .{ .ammonium_dissolution_fraction_per_step = 0.1, .ammonia_dissolution_fraction_per_step = 0.1, .nitrate_dissolution_fraction_per_step = 0.1, .minimum_urea_half_saturation_mol_n_per_megagram = 1, .microbial_activity_inhibition_g_c_per_m3_per_h = 1, .specific_urea_hydrolysis_mol_n_per_g_c = 1, .urease_inhibition_decline_fraction_per_step = 0.1 }));
    try std.testing.expectEqualDeep(before, state.cells[0]);
}

test "surface fertilizer dissolution transfers nitrogen atomically into litter water" {
    var fertilizer = try State.init(std.testing.allocator, 1);
    defer fertilizer.deinit();
    fertilizer.cells[0] = .{ .ammonium_mol_n = 10, .ammonia_mol_n = 8, .urea_mol_n = 6, .nitrate_mol_n = 4, .initial_urease_inhibition_fraction = 0, .current_urease_inhibition_fraction = 0 };
    var chemistry = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    const fluxes = try calculateAndApplyToChemistry(&fertilizer, &chemistry, 0, .{ .litter_water_volume_m3 = 2, .litter_water_volume_fraction = 0.5, .litter_dry_mass_megagrams = 2, .biologically_active_water_volume_m3 = 1, .active_biomass_respiration_g_c_per_step = 2, .microbial_temperature_factor = 1, .step_duration_h = 1 }, .{ .ammonium_dissolution_fraction_per_step = 0.2, .ammonia_dissolution_fraction_per_step = 0.25, .nitrate_dissolution_fraction_per_step = 0.5, .minimum_urea_half_saturation_mol_n_per_megagram = 1, .microbial_activity_inhibition_g_c_per_m3_per_h = 2, .specific_urea_hydrolysis_mol_n_per_g_c = 1, .urease_inhibition_decline_fraction_per_step = 0.1 });
    try std.testing.expectApproxEqAbs(fluxes.ammonium_dissolution_mol_n / 2, chemistry.cells[0].ammonium_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs((fluxes.ammonia_dissolution_mol_n + fluxes.urea_hydrolysis_mol_n) / 2, chemistry.cells[0].ammonia_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(fluxes.nitrate_dissolution_mol_n / 2, chemistry.cells[0].nitrate_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 28), inventoryNitrogen(fertilizer.cells[0]) + aqueousNitrogen(chemistry.cells[0]) * 2, 1e-12);
}
