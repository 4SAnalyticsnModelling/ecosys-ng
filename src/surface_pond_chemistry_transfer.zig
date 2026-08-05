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
    surface_dry_mass_before_megagrams: f64,
    soil_dry_mass_before_megagrams: f64,
    surface_dry_mass_after_megagrams: f64,
    soil_dry_mass_after_megagrams: f64,
    /// Non-band reaction concentrations occupy only these shares of matrix
    /// water. Surface pond nitrogen is incorporated into the non-band zone.
    ammonium_non_band_water_fraction: f64 = 1,
    nitrate_non_band_water_fraction: f64 = 1,
    // Fraction of dissolved (water-borne) chemistry to transfer. This may be
    // less than the dry-material fraction when water transfer is capped by
    // available pore capacity. Set equal to fraction when unconstrained.
    dissolved_chemistry_fraction: f64,
};

const Candidate = struct {
    surface: surface_module.Cell,
    soil_aqueous: @import("solute_aqueous_network.zig").State,
    soil_phosphate: @import("solute_phosphate_network.zig").State,
    soil_cations: @import("solute_cation_exchange.zig").Cations,
    soil_carboxyl_hydrogen_mol_per_megagram: f64,
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
    soil.cation_exchange_mol_per_megagram[destination] = next.soil_cations;
    soil.carboxyl_bound_hydrogen_mol_per_megagram[destination] = next.soil_carboxyl_hydrogen_mol_per_megagram;
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
    // A dry surface carrier (`surface_water_before_m3 == 0`) is a legitimate
    // state once surface evaporation is active: the pond can dry out while its
    // solid phase still transfers. `transferMappedConcentration` already handles
    // it correctly, because a zero source carrier makes the source amount zero so
    // nothing dissolved moves, and it guards its own division. Requiring a
    // positive surface carrier here therefore rejected a valid state and was the
    // second blocker on enabling evaporation. See EXEC-004. The soil-side
    // carriers are still required to be positive: they are the divisors.
    if (carriers.soil_shared_water_before_m3 < 0 or carriers.soil_phosphate_non_band_water_before_m3 < 0 or carriers.soil_shared_water_after_m3 <= 0 or carriers.soil_phosphate_non_band_water_after_m3 <= 0 or carriers.soil_dry_mass_before_megagrams <= 0 or carriers.soil_dry_mass_after_megagrams <= 0 or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1 or carriers.dissolved_chemistry_fraction > 1) {
        if (!builtin.is_test) std.log.err(
            "invalid surface pond chemistry carriers: surface_water_before_m3={e} soil_shared_water_before_m3={e} soil_phosphate_water_before_m3={e} surface_water_after_m3={e} soil_shared_water_after_m3={e} soil_phosphate_water_after_m3={e} surface_dry_mass_before_megagrams={e} soil_dry_mass_before_megagrams={e} surface_dry_mass_after_megagrams={e} soil_dry_mass_after_megagrams={e} fraction={e}",
            .{
                carriers.surface_water_before_m3,
                carriers.soil_shared_water_before_m3,
                carriers.soil_phosphate_non_band_water_before_m3,
                carriers.surface_water_after_m3,
                carriers.soil_shared_water_after_m3,
                carriers.soil_phosphate_non_band_water_after_m3,
                carriers.surface_dry_mass_before_megagrams,
                carriers.soil_dry_mass_before_megagrams,
                carriers.surface_dry_mass_after_megagrams,
                carriers.soil_dry_mass_after_megagrams,
                fraction,
            },
        );
        return error.InvalidSurfacePondChemistryCarrier;
    }

    var result: Candidate = .{
        .surface = surface.cells[cell],
        .soil_aqueous = soil.aqueous[destination],
        .soil_phosphate = soil.non_band_phosphate[destination],
        .soil_cations = soil.cation_exchange_mol_per_megagram[destination],
        .soil_carboxyl_hydrogen_mol_per_megagram = soil.carboxyl_bound_hydrogen_mol_per_megagram[destination],
        .soil_solids = soil.geochemistry_solids[destination],
    };

    inline for (.{
        .{ "hydrogen_mol_per_m3", "hydrogen" },
        .{ "hydroxide_mol_per_m3", "hydroxide" },
        .{ "aluminum_mol_per_m3", "aluminum" },
        .{ "iron_mol_per_m3", "iron" },
        .{ "calcium_mol_per_m3", "calcium" },
        .{ "magnesium_mol_per_m3", "magnesium" },
        .{ "sodium_mol_per_m3", "sodium" },
        .{ "potassium_mol_per_m3", "potassium" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_aqueous, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "ammonium_mol_per_m3", "ammonium_non_band" },
        .{ "ammonia_mol_per_m3", "ammonia_non_band" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_aqueous, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3 * carriers.ammonium_non_band_water_fraction, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3 * carriers.ammonium_non_band_water_fraction, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;
    transferMappedConcentration(&result.surface.nitrate_mol_per_m3, &result.soil_aqueous.nitrate_non_band, carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3 * carriers.nitrate_non_band_water_fraction, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3 * carriers.nitrate_non_band_water_fraction, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;

    // Pond material enters the non-band zone. Existing band concentrations
    // receive no solute, but their shared water carrier grows; dilute them so
    // the extensive band inventory remains unchanged.
    inline for (.{ "ammonium_band", "ammonia_band", "nitrate_band" }) |field_name|
        @field(result.soil_aqueous, field_name) *= carriers.soil_shared_water_before_m3 / carriers.soil_shared_water_after_m3;

    if (dynamic_salts) inline for (.{
        .{ "chloride_mol_per_m3", "chloride" },
        .{ "sulfate_mol_per_m3", "sulfate" },
        .{ "carbon_dioxide_mol_per_m3", "carbon_dioxide" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_aqueous, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;

    // Carbonate alkalinity follows liquid water even when optional salt
    // reactions are disabled. The feature flag controls salt chemistry, not
    // ownership of conserved carbon carriers.
    inline for (.{
        .{ "carbonate_mol_per_m3", "carbonate" },
        .{ "bicarbonate_mol_per_m3", "bicarbonate" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_aqueous, names[1]), carriers.surface_water_before_m3, carriers.soil_shared_water_before_m3, carriers.surface_water_after_m3, carriers.soil_shared_water_after_m3, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "hpo4_mol_p_per_m3", "dissolved_hpo4_mol_p_per_m3" },
        .{ "h2po4_mol_p_per_m3", "dissolved_h2po4_mol_p_per_m3" },
    }) |names| transferMappedConcentration(&@field(result.surface, names[0]), &@field(result.soil_phosphate, names[1]), carriers.surface_water_before_m3, carriers.soil_phosphate_non_band_water_before_m3, carriers.surface_water_after_m3, carriers.soil_phosphate_non_band_water_after_m3, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "ammonium_mol_per_megagram", "ammonium_non_band" },
        .{ "hydrogen_mol_per_megagram", "hydrogen" },
        .{ "aluminum_mol_per_megagram", "aluminum" },
        .{ "iron_mol_per_megagram", "iron" },
        .{ "calcium_mol_per_megagram", "calcium" },
        .{ "magnesium_mol_per_megagram", "magnesium" },
        .{ "sodium_mol_per_megagram", "sodium" },
        .{ "potassium_mol_per_megagram", "potassium" },
    }) |names| transferMappedConcentration(&@field(result.surface.exchange, names[0]), &@field(result.soil_cations, names[1]), carriers.surface_dry_mass_before_megagrams, carriers.soil_dry_mass_before_megagrams, carriers.surface_dry_mass_after_megagrams, carriers.soil_dry_mass_after_megagrams, fraction) catch return error.InvalidSurfacePondChemistryState;

    // As above, only non-band exchange receives surface material. Preserve
    // the pre-existing band ammonium amount on the enlarged dry-mass carrier.
    result.soil_cations.ammonium_band *= carriers.soil_dry_mass_before_megagrams / carriers.soil_dry_mass_after_megagrams;

    transferMappedConcentration(&result.surface.carboxyl_hydrogen_mol_per_megagram, &result.soil_carboxyl_hydrogen_mol_per_megagram, carriers.surface_dry_mass_before_megagrams, carriers.soil_dry_mass_before_megagrams, carriers.surface_dry_mass_after_megagrams, carriers.soil_dry_mass_after_megagrams, fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "deprotonated_site_mol_per_megagram", "deprotonated_site_mol_per_megagram" },
        .{ "hydroxyl_site_mol_per_megagram", "hydroxyl_site_mol_per_megagram" },
        .{ "protonated_site_mol_per_megagram", "protonated_site_mol_per_megagram" },
        .{ "adsorbed_hpo4_mol_p_per_megagram", "adsorbed_hpo4_mol_p_per_megagram" },
        .{ "adsorbed_h2po4_mol_p_per_megagram", "adsorbed_h2po4_mol_p_per_megagram" },
    }) |names| transferMappedConcentration(&@field(result.surface.phosphate_surface, names[0]), &@field(result.soil_phosphate, names[1]), carriers.surface_dry_mass_before_megagrams, carriers.soil_dry_mass_before_megagrams, carriers.surface_dry_mass_after_megagrams, carriers.soil_dry_mass_after_megagrams, fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "aluminum_phosphate_mol_per_m3", "aluminum_phosphate_solid_mol_per_m3" },
        .{ "iron_phosphate_mol_per_m3", "iron_phosphate_solid_mol_per_m3" },
        .{ "dicalcium_phosphate_mol_per_m3", "dicalcium_phosphate_solid_mol_per_m3" },
        .{ "hydroxyapatite_mol_per_m3", "hydroxyapatite_solid_mol_per_m3" },
        .{ "monocalcium_phosphate_mol_per_m3", "monocalcium_phosphate_solid_mol_per_m3" },
    }) |names| transferMappedConcentration(&@field(result.surface.phosphate_minerals, names[0]), &@field(result.soil_phosphate, names[1]), carriers.surface_water_before_m3, carriers.soil_phosphate_non_band_water_before_m3, carriers.surface_water_after_m3, carriers.soil_phosphate_non_band_water_after_m3, carriers.dissolved_chemistry_fraction) catch return error.InvalidSurfacePondChemistryState;

    inline for (.{
        .{ "gibbsite_mol_per_m3", "gibbsite_solid_mol_per_m3" },
        .{ "iron_hydroxide_mol_per_m3", "iron_hydroxide_solid_mol_per_m3" },
        .{ "calcite_mol_per_m3", "calcite_solid_mol_per_m3" },
        .{ "gypsum_mol_per_m3", "gypsum_solid_mol_per_m3" },
    }) |names| transferMappedConcentration(
        &@field(result.surface.salt_minerals, names[0]),
        &@field(result.soil_solids, names[1]),
        carriers.surface_water_before_m3,
        carriers.soil_shared_water_before_m3,
        carriers.surface_water_after_m3,
        carriers.soil_shared_water_after_m3,
        if (dynamic_salts) carriers.dissolved_chemistry_fraction else 0,
    ) catch return error.InvalidSurfacePondChemistryState;

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
    surface.cells[0].exchange.calcium_mol_per_megagram = 3;
    soil.aqueous[0].ammonium_non_band = 1;
    const carriers: CarrierVolumes = .{ .surface_water_before_m3 = 2, .soil_shared_water_before_m3 = 2, .soil_phosphate_non_band_water_before_m3 = 2, .surface_water_after_m3 = 1, .soil_shared_water_after_m3 = 3, .soil_phosphate_non_band_water_after_m3 = 3, .surface_dry_mass_before_megagrams = 2, .soil_dry_mass_before_megagrams = 2, .surface_dry_mass_after_megagrams = 1, .soil_dry_mass_after_megagrams = 3, .dissolved_chemistry_fraction = 0.5 };
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    try std.testing.expectEqual(@as(f64, 4), surface.cells[0].ammonium_mol_per_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 2), soil.aqueous[0].ammonium_non_band, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 3.0), soil.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), soil.cation_exchange_mol_per_megagram[0].calcium, 1e-14);
}

test "pond mineral nitrogen conserves amount on non-band zone water" {
    var surface = try surface_module.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try soil_module.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    surface.cells[0].ammonium_mol_per_m3 = 4;
    surface.cells[0].nitrate_mol_per_m3 = 2;
    soil.aqueous[0].ammonium_non_band = 1;
    soil.aqueous[0].nitrate_non_band = 3;
    soil.aqueous[0].ammonium_band = 7;
    soil.aqueous[0].nitrate_band = 11;
    soil.cation_exchange_mol_per_megagram[0].ammonium_band = 13;
    const carriers: CarrierVolumes = .{
        .surface_water_before_m3 = 2,
        .soil_shared_water_before_m3 = 5,
        .soil_phosphate_non_band_water_before_m3 = 5,
        .surface_water_after_m3 = 1,
        .soil_shared_water_after_m3 = 6,
        .soil_phosphate_non_band_water_after_m3 = 6,
        .surface_dry_mass_before_megagrams = 2,
        .soil_dry_mass_before_megagrams = 5,
        .surface_dry_mass_after_megagrams = 1,
        .soil_dry_mass_after_megagrams = 6,
        .ammonium_non_band_water_fraction = 0.8,
        .nitrate_non_band_water_fraction = 0.6,
        .dissolved_chemistry_fraction = 0.5,
    };
    const ammonium_before = 4.0 * 2.0 + 1.0 * 5.0 * 0.8;
    const nitrate_before = 2.0 * 2.0 + 3.0 * 5.0 * 0.6;
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    const ammonium_after = surface.cells[0].ammonium_mol_per_m3 * 1.0 + soil.aqueous[0].ammonium_non_band * 6.0 * 0.8;
    const nitrate_after = surface.cells[0].nitrate_mol_per_m3 * 1.0 + soil.aqueous[0].nitrate_non_band * 6.0 * 0.6;
    try std.testing.expectApproxEqAbs(ammonium_before, ammonium_after, 1e-14);
    try std.testing.expectApproxEqAbs(nitrate_before, nitrate_after, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 7 * 5), soil.aqueous[0].ammonium_band * 6, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 11 * 5), soil.aqueous[0].nitrate_band * 6, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 13 * 5), soil.cation_exchange_mol_per_megagram[0].ammonium_band * 6, 1e-14);
}

test "fixed salt mode preserves optional salt amounts while water carrier changes" {
    var surface = try surface_module.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try soil_module.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    surface.cells[0].sulfate_mol_per_m3 = 5;
    surface.cells[0].bicarbonate_mol_per_m3 = 4;
    const carriers: CarrierVolumes = .{ .surface_water_before_m3 = 1, .soil_shared_water_before_m3 = 1, .soil_phosphate_non_band_water_before_m3 = 1, .surface_water_after_m3 = 0.5, .soil_shared_water_after_m3 = 1.5, .soil_phosphate_non_band_water_after_m3 = 1.5, .surface_dry_mass_before_megagrams = 1, .soil_dry_mass_before_megagrams = 1, .surface_dry_mass_after_megagrams = 0.5, .soil_dry_mass_after_megagrams = 1.5, .dissolved_chemistry_fraction = 0.5 };
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    try std.testing.expectEqual(@as(f64, 5), surface.cells[0].sulfate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0), soil.aqueous[0].sulfate);
    try std.testing.expectEqual(@as(f64, 4), surface.cells[0].bicarbonate_mol_per_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0 / 3.0), soil.aqueous[0].bicarbonate, 1e-15);
    surface.cells[0].salt_minerals.calcite_mol_per_m3 = 4;
    soil.geochemistry_solids[0].calcite_solid_mol_per_m3 = 2;
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    try std.testing.expectEqual(@as(f64, 8), surface.cells[0].salt_minerals.calcite_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4.0 / 3.0), soil.geochemistry_solids[0].calcite_solid_mol_per_m3);
    const calcite_after_mol =
        surface.cells[0].salt_minerals.calcite_mol_per_m3 * carriers.surface_water_after_m3 +
        soil.geochemistry_solids[0].calcite_solid_mol_per_m3 * carriers.soil_shared_water_after_m3;
    try std.testing.expectEqual(@as(f64, 6), calcite_after_mol);
}

test "dissolved chemistry fraction zero leaves soil aqueous state unchanged while dry chemistry still transfers" {
    var surface = try surface_module.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try soil_module.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    surface.cells[0].ammonium_mol_per_m3 = 1_000_000;
    surface.cells[0].exchange.ammonium_mol_per_megagram = 100;
    soil.aqueous[0].ammonium_non_band = 5;
    soil.cation_exchange_mol_per_megagram[0].ammonium_non_band = 0;
    // Water transfer blocked (pore full): dissolved_chemistry_fraction = 0.
    // Dry mass transfers at the full fraction = 0.5.
    const carriers: CarrierVolumes = .{
        .surface_water_before_m3 = 10,
        .soil_shared_water_before_m3 = 2,
        .soil_phosphate_non_band_water_before_m3 = 2,
        .surface_water_after_m3 = 10,
        .soil_shared_water_after_m3 = 2,
        .soil_phosphate_non_band_water_after_m3 = 2,
        .surface_dry_mass_before_megagrams = 2,
        .soil_dry_mass_before_megagrams = 2,
        .surface_dry_mass_after_megagrams = 1,
        .soil_dry_mass_after_megagrams = 3,
        .dissolved_chemistry_fraction = 0,
    };
    try transferSurfaceFractionToSoil(&surface, &soil, 0, 0, carriers, false, 0.5);
    // Dissolved: no transfer, soil stays at 5
    try std.testing.expectEqual(@as(f64, 5), soil.aqueous[0].ammonium_non_band);
    // Dissolved: source concentration unchanged (no water moved)
    try std.testing.expectEqual(@as(f64, 1_000_000), surface.cells[0].ammonium_mol_per_m3);
    // Adsorbed: half of surface (0.5 * 100 * 2 / 3 Mg) moves to soil
    try std.testing.expectApproxEqAbs(@as(f64, 100.0 / 3.0), soil.cation_exchange_mol_per_megagram[0].ammonium_non_band, 1e-12);
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
        .surface_dry_mass_before_megagrams = 0,
        .soil_dry_mass_before_megagrams = 2,
        .surface_dry_mass_after_megagrams = 0,
        .soil_dry_mass_after_megagrams = 2,
        .dissolved_chemistry_fraction = 1,
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
