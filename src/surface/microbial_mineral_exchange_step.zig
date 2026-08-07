const std = @import("std");
const compute = @import("../core/compute.zig");
const organic = @import("../soil/organic/initialization.zig");
const chemistry = @import("litter_chemistry.zig");
const respiration = @import("microbial_respiration_step.zig");

pub const Parameters = struct {
    ammonium_maximum_uptake_g_n_per_m2_h: f64,
    ammonium_minimum_concentration_g_n_per_m3: f64,
    ammonium_half_saturation_g_n_per_m3: f64,
    nitrate_maximum_uptake_g_n_per_m2_h: f64,
    nitrate_minimum_concentration_g_n_per_m3: f64,
    nitrate_half_saturation_g_n_per_m3: f64,
    phosphate_maximum_uptake_g_p_per_m2_h: f64,
    phosphate_minimum_concentration_g_p_per_m3: f64,
    phosphate_half_saturation_g_p_per_m3: f64,
    phosphorus_molar_mass_g_per_mol: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    ammonium_exchange_g_n: []f64,
    nitrate_exchange_g_n: []f64,
    h2po4_exchange_g_p: []f64,
    hpo4_exchange_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroSurfaceMineralExchangeCells;
        const count = try std.math.mul(usize, cell_count, respiration.unit_count_per_cell);
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

pub const ApplyContext = struct {
    result: *State,
    surface_organic: *const organic.State,
    litter_chemistry: *const chemistry.State,
    litter_water_m3: []const f64,
    growth_temperature_response: []const f64,
    matric_plus_osmotic_potential_megapascal: []const f64,
    microbial_nitrogen_to_carbon_g_n_per_g_c: []const f64,
    microbial_phosphorus_to_carbon_g_p_per_g_c: []const f64,
    labile_biomass_fraction: f64,
    microbial_surface_area_m2_per_g_c: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    parameters: Parameters,
    timestep_h: f64,
};

/// NITRO L=0 RINH4/RINO3/RIPO4/RIP14 for litter mineral pools.
/// Positive values immobilize mineral nutrient into OMC3; negative values
/// mineralize microbial storage back to the corresponding aqueous pool.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    for (range.first..range.end) |cell| {
        const water = context.litter_water_m3[cell];
        const c = context.litter_chemistry.cells[cell];
        const ammonium_g_n = c.ammonium_mol_per_m3 * water * context.nitrogen_molar_mass_g_per_mol;
        const nitrate_g_n = c.nitrate_mol_per_m3 * water * context.nitrogen_molar_mass_g_per_mol;
        const h2po4_g_p = c.h2po4_mol_p_per_m3 * water * context.parameters.phosphorus_molar_mass_g_per_mol;
        const hpo4_g_p = c.hpo4_mol_p_per_m3 * water * context.parameters.phosphorus_molar_mass_g_per_mol;
        const ammonium_concentration = c.ammonium_mol_per_m3 * context.nitrogen_molar_mass_g_per_mol;
        const nitrate_concentration = c.nitrate_mol_per_m3 * context.nitrogen_molar_mass_g_per_mol;
        const h2po4_concentration = c.h2po4_mol_p_per_m3 * context.parameters.phosphorus_molar_mass_g_per_mol;
        const hpo4_concentration = c.hpo4_mol_p_per_m3 * context.parameters.phosphorus_molar_mass_g_per_mol;
        var total_active: f64 = 0;
        var active: [respiration.unit_count_per_cell]f64 = undefined;
        for (0..respiration.litter_complex_count) |complex| for (0..respiration.source_population_count) |population| {
            const local = complex * respiration.source_population_count + population;
            const microbial = ((cell * organic.microbial_substrate_count + complex) * organic.microbial_population_count + population) * organic.kinetic_fraction_count;
            active[local] = context.surface_organic.microbial[microbial].carbon_g_c / context.labile_biomass_fraction;
            total_active += active[local];
        };
        for (0..respiration.litter_complex_count) |complex| for (0..respiration.source_population_count) |population| {
            const local = complex * respiration.source_population_count + population;
            const unit = cell * respiration.unit_count_per_cell + local;
            const microbial = ((cell * organic.microbial_substrate_count + complex) * organic.microbial_population_count + population) * organic.kinetic_fraction_count;
            const nonstructural = context.surface_organic.microbial[microbial + 2];
            const ratio = (complex * organic.microbial_population_count + population) * organic.kinetic_fraction_count + 2;
            const nitrogen_demand = nonstructural.carbon_g_c * context.microbial_nitrogen_to_carbon_g_n_per_g_c[ratio] - nonstructural.nitrogen_g_n;
            const phosphorus_demand = nonstructural.carbon_g_c * context.microbial_phosphorus_to_carbon_g_p_per_g_c[ratio] - nonstructural.phosphorus_g_p;
            const share = if (total_active > 0) active[local] / total_active else 0;
            const activity = context.growth_temperature_response[cell] * @exp((if (population == 2) @as(f64, 0.05) else 0.10) * context.matric_plus_osmotic_potential_megapascal[cell]);
            const area_activity = context.microbial_surface_area_m2_per_g_c * active[local] * activity * context.timestep_h;
            const ammonium = exchange(nitrogen_demand, ammonium_g_n, ammonium_concentration, context.parameters.ammonium_minimum_concentration_g_n_per_m3, context.parameters.ammonium_half_saturation_g_n_per_m3, context.parameters.ammonium_maximum_uptake_g_n_per_m2_h * area_activity, share, water);
            const nitrate_demand = @max(0, nitrogen_demand - ammonium);
            const nitrate = if (nitrate_demand > 0) exchange(nitrate_demand, nitrate_g_n, nitrate_concentration, context.parameters.nitrate_minimum_concentration_g_n_per_m3, context.parameters.nitrate_half_saturation_g_n_per_m3, context.parameters.nitrate_maximum_uptake_g_n_per_m2_h * area_activity, share, water) else 0;
            const h2po4 = exchange(phosphorus_demand, h2po4_g_p, h2po4_concentration, context.parameters.phosphate_minimum_concentration_g_p_per_m3, context.parameters.phosphate_half_saturation_g_p_per_m3, context.parameters.phosphate_maximum_uptake_g_p_per_m2_h * area_activity, share, water);
            const hpo4_demand = @max(0, phosphorus_demand - h2po4);
            const hpo4 = if (hpo4_demand > 0) exchange(hpo4_demand, hpo4_g_p, hpo4_concentration, 0.25 * context.parameters.phosphate_minimum_concentration_g_p_per_m3, context.parameters.phosphate_half_saturation_g_p_per_m3, 0.25 * context.parameters.phosphate_maximum_uptake_g_p_per_m2_h * area_activity, share, water) else 0;
            context.result.ammonium_exchange_g_n[unit] = ammonium;
            context.result.nitrate_exchange_g_n[unit] = nitrate;
            context.result.h2po4_exchange_g_p[unit] = h2po4;
            context.result.hpo4_exchange_g_p[unit] = hpo4;
        };
    }
}

fn exchange(demand: f64, amount: f64, concentration: f64, minimum_concentration: f64, half_saturation: f64, uptake_capacity: f64, share: f64, water_m3: f64) f64 {
    if (demand <= 0) {
        // Mineralization (N/P surplus): cap release at the same enzymatic rate capacity
        // as immobilization so surplus does not flush into the aqueous pool all at once.
        return -@min(-demand, uptake_capacity);
    }
    const available_concentration = @max(0, concentration - minimum_concentration);
    const unlimited = @min(demand, uptake_capacity) * available_concentration / (available_concentration + half_saturation);
    const reserve = minimum_concentration * water_m3;
    return @max(0, @min(unlimited, share * @max(0, amount - reserve)));
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const ratios = organic.substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count;
    if (range.first > range.end or range.end > context.result.cell_count or context.surface_organic.layer_count != context.result.cell_count or context.litter_chemistry.cells.len != context.result.cell_count or context.litter_water_m3.len != context.result.cell_count or context.growth_temperature_response.len != context.result.cell_count or context.matric_plus_osmotic_potential_megapascal.len != context.result.cell_count or context.microbial_nitrogen_to_carbon_g_n_per_g_c.len != ratios or context.microbial_phosphorus_to_carbon_g_p_per_g_c.len != ratios) return error.SurfaceMineralExchangeDimensionMismatch;
    inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(context.parameters, field.name)) or @field(context.parameters, field.name) <= 0) return error.InvalidSurfaceMineralExchangeParameter;
    inline for (.{ context.labile_biomass_fraction, context.microbial_surface_area_m2_per_g_c, context.nitrogen_molar_mass_g_per_mol, context.timestep_h }) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSurfaceMineralExchangeParameter;
    if (context.labile_biomass_fraction > 1) return error.InvalidSurfaceMineralExchangeParameter;
}

test "surface mineral exchange immobilizes deficits and mineralizes surpluses" {
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    organic_state.microbial[0].carbon_g_c = 0.55;
    organic_state.microbial[2] = .{ .carbon_g_c = 1, .nitrogen_g_n = 0, .phosphorus_g_p = 0.2 };
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    chemistry_state.cells[0].ammonium_mol_per_m3 = 1;
    chemistry_state.cells[0].nitrate_mol_per_m3 = 1;
    chemistry_state.cells[0].h2po4_mol_p_per_m3 = 1;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const n = [_]f64{0.1} ** (organic.substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count);
    const p = [_]f64{0.01} ** (organic.substrate_count * organic.microbial_population_count * organic.kinetic_fraction_count);
    var context: ApplyContext = .{ .result = &state, .surface_organic = &organic_state, .litter_chemistry = &chemistry_state, .litter_water_m3 = &.{1}, .growth_temperature_response = &.{1}, .matric_plus_osmotic_potential_megapascal = &.{0}, .microbial_nitrogen_to_carbon_g_n_per_g_c = &n, .microbial_phosphorus_to_carbon_g_p_per_g_c = &p, .labile_biomass_fraction = 0.55, .microbial_surface_area_m2_per_g_c = 3, .nitrogen_molar_mass_g_per_mol = 14, .parameters = .{ .ammonium_maximum_uptake_g_n_per_m2_h = 0.014, .ammonium_minimum_concentration_g_n_per_m3 = 0.0125, .ammonium_half_saturation_g_n_per_m3 = 0.40, .nitrate_maximum_uptake_g_n_per_m2_h = 0.014, .nitrate_minimum_concentration_g_n_per_m3 = 0.03, .nitrate_half_saturation_g_n_per_m3 = 0.35, .phosphate_maximum_uptake_g_p_per_m2_h = 0.003, .phosphate_minimum_concentration_g_p_per_m3 = 0.009, .phosphate_half_saturation_g_p_per_m3 = 0.18, .phosphorus_molar_mass_g_per_mol = 31 }, .timestep_h = 1 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(state.ammonium_exchange_g_n[0] > 0);
    try std.testing.expect(state.h2po4_exchange_g_p[0] < 0);
}
