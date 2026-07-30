const std = @import("std");
const builtin = @import("builtin");
const surface_module = @import("surface_litter_chemistry.zig");
const soil_module = @import("solute_chemistry_state.zig");

pub const CarrierVolumes = struct {
    surface_water_before_m3: f64,
    soil_shared_water_before_m3: f64,
    soil_phosphate_non_band_water_before_m3: f64,
    surface_water_after_m3: f64,
    soil_shared_water_after_m3: f64,
    soil_phosphate_non_band_water_after_m3: f64,
    surface_dry_mass_before_Mg: f64,
    soil_dry_mass_before_Mg: f64,
    surface_dry_mass_after_Mg: f64,
    soil_dry_mass_after_Mg: f64,
};

const Candidate = struct {
    surface: surface_module.Cell,
    soil_aqueous: @import("solute_aqueous_network.zig").State,
    soil_phosphate: @import("solute_phosphate_network.zig").State,
    soil_cations: @import("solute_cation_exchange.zig").Cations,
    soil_carboxyl_hydrogen_mol_per_Mg: f64,
    soil_solids: @import("solute_geochemistry_network.zig").SolidState,
};

/// REDIST layer-zero chemistry transfer. Surface concentrations are converted
/// to extensive amounts on their native water or litter-mass carriers and
/// mixed into the destination soil non-band owners.
pub fn transferSurfaceFractionToSoil(
    surface: *surface_module.State,
    soil: *soil_module.State,
    cell: usize,
    destination: usize,
    carriers: CarrierVolumes,
    dynamic_salts: bool,
    fraction: f64,
) !void {
    const next = try calculate(surface, soil, cell, destination, carriers, dynamic_salts, fraction);
    surface.cells[cell] = next.surface;
    soil.aqueous[destination] = next.soil_aqueous;
    soil.non_band_phosphate[destination] = next.soil_phosphate;
    soil.cation_exchange_mol_per_Mg[destination] = next.soil_cations;
    soil.carboxyl_bound_hydrogen_mol_per_Mg[destination] = next.soil_carboxyl_hydrogen_mol_per_Mg;
    soil.geochemistry_solids[destination] = next.soil_solids;
}

pub fn validateSurfaceFractionToSoil(
    surface: *const surface_module.State,
    soil: *const soil_module.State,
    cell: usize,
    destination: usize,
    carriers: CarrierVolumes,
    dynamic_salts: bool,
    fraction: f64,
) !void {
    _ = try calculate(surface, soil, cell, destination, carriers, dynamic_salts, fraction);
}

fn calculate(
    surface: *const surface_module.State,
    soil: *const soil_module.State,
    cell: usize,
    destination: usize,
    carriers: CarrierVolumes,
    dynamic_salts: bool,
    fraction: f64,
) !Candidate {
    if (cell >= surface.cells.len or destination >= soil.cell_count) return error.SurfacePondChemistryIndexOutOfBounds;
    inline for (@typeInfo(CarrierVolumes).@"struct".fields) |field| {
        const value = @field(carriers, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfacePondChemistryCarrier;
    }
    if (carriers.surface_water_before_m3 <= 0 or carriers.soil_shared_water_before_m3 < 0 or carriers.soil_phosphate_non_band_water_before_m3 < 0 or carriers.soil_shared_water_after_m3 <= 0 or carriers.soil_phosphate_non_band_water_after_m3 <= 0 or carriers.soil_dry_mass_before_Mg <= 0 or carriers.soil_dry_mass_after_Mg <= 0 or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1) {
        if (!builtin.is_test) std.log.err(
            "invalid surface pond chemistry carriers: surface_water_before_m3={e} soil_shared_water_before_m3={e} soil_phosphate_water_before_m3={e} surface_water_after_m3={e} soil_shared_water_after_m3={e} soil_phosphate_water_after_m3={e} surface_dry_mass_before_Mg={e} soil_dry_mass_before_Mg={e} surface_dry_mass_after_Mg={e} soil_dry_mass_after_Mg={e} fraction={e}",
            .{
                carriers.surface_water_before_m3,
                carriers.soil_shared_water_before_m3,
                carriers.soil_phosphate_non_band_water_before_m3,
                carriers.surface_water_after_m3,
                carriers.soil_shared_water_after_m3,
                carriers.soil_phosphate_non_band_water_after_m3,
                carriers.surface_dry_mass_before_Mg,
                carriers.soil_dry_mass_before_Mg,
                carriers.surface_dry_mass_after_Mg,
                carriers.soil_dry_mass_after_Mg,
                fraction,
            },
        );
        return error.InvalidSurfacePondChemistryCarrier;
    }

    var result: Candidate = .{
        .surface = surface.cells[cell],
        .soil_aqueous = soil.aqueous[destination],
        .soil_phosphate = soil.non_band_phosphate[destination],
        .soil_cations = soil.cation_exchange_mol_per_Mg[destination],
        .soil_carboxyl_hydrogen_mol_per_Mg = soil.carboxyl_bound_hydrogen_mol_per_Mg[destination],
        .soil_solids = soil.geochemistry_solids[destination],
    };

    inline for (.{
        .{ "ammonium_mol_per_m3", "ammonium_non_band" },
        .{ "ammonia_mol_per_m3", "ammonia_non_band" },
        .{ "nitrate_mol_per_m3", "nitrate_non_band" },
        .{ "hydrogen_mol_per_m3", "hydrogen" },
        .{ "hydroxide_mol_per_m3", "hydroxide" },
        .{ "aluminum_mol_per_m3", "aluminum" },
        .{ "iron_mol_per_m3", "iron" },
        .{ "calcium_mol_per_m3", "calcium" },
        .{ "magnesium_mol_per_m3", "magnesium" },
        .{ "sodium_mol_per_m3", "sodium" },
        .{ "potassium_mol_per_m3", "potassium" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_aqueous, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3, fraction) catch return error.InvalidSurfacePondChemistryState;

    if (dynamic_salts) inline for (.{
        .{ "chloride_mol_per_m3", "chloride" },
        .{ "sulfate_mol_per_m3", "sulfate" },
        .{ "carbonate_mol_per_m3", "carbonate" },
        .{ "bicarbonate_mol_per_m3", "bicarbonate" },
        .{ "carbon_dioxide_mol_per_m3", "carbon_dioxide" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_aqueous, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3, fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "hpo4_mol_p_per_m3", "dissolved_hpo4_mol_p_per_m3" },
        .{ "h2po4_mol_p_per_m3", "dissolved_h2po4_mol_p_per_m3" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_phosphate, names[1]), carriers.surface_water_before_m3, carriers.soil_phosphate_non_band_water_before_m3, carriers.surface_water_after_m3, carriers.soil_phosphate_non_band_water_after_m3, fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "ammonium_mol_per_Mg", "ammonium_non_band" },
        .{ "hydrogen_mol_per_Mg", "hydrogen" },
        .{ "aluminum_mol_per_Mg", "aluminum" },
        .{ "iron_mol_per_Mg", "iron" },
        .{ "calcium_mol_per_Mg", "calcium" },
        .{ "magnesium_mol_per_Mg", "magnesium" },
        .{ "sodium_mol_per_Mg", "sodium" },
        .{ "potassium_mol_per_Mg", "potassium" },
    }) |names| transferMappedConcentration(&@field(result.surface.exchange, names[0]), &@field(result.soil_cations, names[1]), carriers.surface_dry_mass_before_Mg, carriers.soil_dry_mass_before_Mg, carriers.surface_dry_mass_after_Mg, carriers.soil_dry_mass_after_Mg, fraction) catch return error.InvalidSurfacePondChemistryState;

    transferMappedConcentration(&result.surface.carboxyl_hydrogen_mol_per_Mg, &result.soil_carboxyl_hydrogen_mol_per_Mg, carriers.surface_dry_mass_before_Mg, carriers.soil_dry_mass_before_Mg, carriers.surface_dry_mass_after_Mg, carriers.soil_dry_mass_after_Mg, fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "deprotonated_site_mol_per_Mg", "deprotonated_site_mol_per_Mg" },
        .{ "hydroxyl_site_mol_per_Mg", "hydroxyl_site_mol_per_Mg" },
        .{ "protonated_site_mol_per_Mg", "protonated_site_mol_per_Mg" },
        .{ "adsorbed_hpo4_mol_p_per_Mg", "adsorbed_hpo4_mol_p_per_Mg" },
        .{ "adsorbed_h2po4_mol_p_per_Mg", "adsorbed_h2po4_mol_p_per_Mg" },
    }) |names| transferMappedConcentration(&@field(result.surface.phosphate_surface, names[0]), &@field(result.soil_phosphate, names[1]), carriers.surface_dry_mass_before_Mg, carriers.soil_dry_mass_before_Mg, carriers.surface_dry_mass_after_Mg, carriers.soil_dry_mass_after_Mg, fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "aluminum_phosphate_mol_per_m3", "aluminum_phosphate_solid_mol_per_m3" },
        .{ "iron_phosphate_mol_per_m3", "iron_phosphate_solid_mol_per_m3" },
        .{ "dicalcium_phosphate_mol_per_m3", "dicalcium_phosphate_solid_mol_per_m3" },
        .{ "hydroxyapatite_mol_per_m3", "hydroxyapatite_solid_mol_per_m3" },
        .{ "monocalcium_phosphate_mol_per_m3", "monocalcium_phosphate_solid_mol_per_m3" },
    }) |names| transferMappedConcentration(&@field(result.surface.phosphate_minerals, names[0]), &@field(result.soil_phosphate, names[1]), carriers.surface_water_before_m3, carriers.soil_phosphate_non_band_water_before_m3, carriers.surface_water_after_m3, carriers.soil_phosphate_non_band_water_after_m3, fraction) catch return error.InvalidSurfacePondChemistryState;

    if (dynamic_salts) inline for (.{
        .{ "gibbsite_mol_per_m3", "gibbsite_solid_mol_per_m3" },
        .{ "iron_hydroxide_mol_per_m3", "iron_hydroxide_solid_mol_per_m3" },
        .{ "calcite_mol_per_m3", "calcite_solid_mol_per_m3" },
        .{ "gypsum_mol_per_m3", "gypsum_solid_mol_per_m3" },
    }) |names| transferMappedConcentration(&@field(result.surface.salt_minerals, names[0]), &@field(result.soil_solids, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3, fraction) catch return error.InvalidSurfacePondChemistryState;

    return result;
}

fn transferMappedConcentration(source: *f64, destination: *f64, source_base_before: f64, destination_base_before: f64, source_base_after: f64, destination_base_after: f64, fraction: f64) !void {
    inline for (.{ source.*, destination.*, source_base_before, destination_base_before, source_base_after, destination_base_after, fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfacePondChemistry;
    if (source.* < 0 or destination.* < 0 or source_base_before < 0 or destination_base_before < 0 or source_base_after < 0 or destination_base_after <= 0 or fraction < 0 or fraction > 1) return error.InvalidSurfacePondChemistry;
    const source_amount = source.* * source_base_before;
    const moved = fraction * source_amount;
    const destination_amount = destination.* * destination_base_before + moved;
    source.* = if (source_base_after > 0) (source_amount - moved) / source_base_after else 0;
    destination.* = destination_amount / destination_base_after;
    if (!std.math.isFinite(source.*) or source.* < 0 or !std.math.isFinite(destination.*) or destination.* < 0) return error.InvalidSurfacePondChemistry;
}

test "surface chemistry mixes into soil non-band owners on native carriers" {
    var surface = try surface_module.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try soil_module.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    surface.cells[0].ammonium_mol_per_m3 = 4;
    surface.cells[0].h2po4_mol_p_per_m3 = 2;
    surface.cells[0].exchange.calcium_mol_per_Mg = 3;
    soil.aqueous[0].ammonium_non_band = 1;
    const carriers: CarrierVolumes = .{ .surface_water_before_m3 = 2, .soil_shared_water_before_m3 = 2, .soil_phosphate_non_band_water_before_m3 = 2, .surface_water_after_m3 = 1, .soil_shared_water_after_m3 = 3, .soil_phosphate_non_band_water_after_m3 = 3, .surface_dry_mass_before_Mg = 2, .soil_dry_mass_before_Mg = 2, .surface_dry_mass_after_Mg = 1, .soil_dry_mass_after_Mg = 3 };
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    try std.testing.expectEqual(@as(f64, 4), surface.cells[0].ammonium_mol_per_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 2), soil.aqueous[0].ammonium_non_band, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 3.0), soil.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), soil.cation_exchange_mol_per_Mg[0].calcium, 1e-14);
}

test "fixed salt mode leaves optional surface salts in place" {
    var surface = try surface_module.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try soil_module.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    surface.cells[0].sulfate_mol_per_m3 = 5;
    const carriers: CarrierVolumes = .{ .surface_water_before_m3 = 1, .soil_shared_water_before_m3 = 1, .soil_phosphate_non_band_water_before_m3 = 1, .surface_water_after_m3 = 0.5, .soil_shared_water_after_m3 = 1.5, .soil_phosphate_non_band_water_after_m3 = 1.5, .surface_dry_mass_before_Mg = 1, .soil_dry_mass_before_Mg = 1, .surface_dry_mass_after_Mg = 0.5, .soil_dry_mass_after_Mg = 1.5 };
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    try std.testing.expectEqual(@as(f64, 5), surface.cells[0].sulfate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0), soil.aqueous[0].sulfate);
}

test "water-only pond chemistry transfers with zero dry carrier" {
    var surface = try surface_module.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try soil_module.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    surface.cells[0].nitrate_mol_per_m3 = 2;
    soil.aqueous[0].nitrate_non_band = 1;
    const carriers: CarrierVolumes = .{
        .surface_water_before_m3 = 3,
        .soil_shared_water_before_m3 = 1,
        .soil_phosphate_non_band_water_before_m3 = 1,
        .surface_water_after_m3 = 0,
        .soil_shared_water_after_m3 = 4,
        .soil_phosphate_non_band_water_after_m3 = 4,
        .surface_dry_mass_before_Mg = 0,
        .soil_dry_mass_before_Mg = 2,
        .surface_dry_mass_after_Mg = 0,
        .soil_dry_mass_after_Mg = 2,
    };
    try transferSurfaceFractionToSoil(
        &surface,
        &soil,
        0,
        0,
        carriers,
        false,
        1,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        surface.cells[0].nitrate_mol_per_m3,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.75),
        soil.aqueous[0].nitrate_non_band,
        1e-14,
    );
}
