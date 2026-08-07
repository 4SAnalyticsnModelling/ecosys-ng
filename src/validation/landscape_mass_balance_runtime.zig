const std = @import("std");
const audit = @import("mass_balance_audit.zig");
const inventory = @import("landscape_mass_inventory.zig");
const boundary = @import("landscape_boundary_ledger.zig");
const grid_module = @import("../state/grid.zig");
const snow_module = @import("../soil/solute/snow_solute_transport.zig");
const thermal_module = @import("../soil/heat/thermal.zig");
const soil_properties_module = @import("../soil/water/solver_properties.zig");
const gas_module = @import("../soil/gas/transport.zig");
const organic_module = @import("../soil/organic/initialization.zig");
const organic_transport_module = @import("../soil/organic/transport.zig");
const mineral_nitrogen_module = @import("../soil/biogeochemistry/mineral_nitrogen_transport.zig");
const chemistry_module = @import("../soil/solute/chemistry_state.zig");
const nitrogen_fertilizer_module = @import("../management/fertilizer_nitrogen_inventory.zig");
const mineral_fertilizer_module = @import("../management/mineral_fertilizer_inventory.zig");
const solute_module = @import("../soil/solute/transport.zig");
const zone_module = @import("../soil/solute/charge_classification.zig");
const litter_chemistry_module = @import("../surface/litter_chemistry.zig");
const litter_fertilizer_module = @import("../surface/litter_fertilizer.zig");
const surface_module = @import("../surface/precipitation.zig");
const canopy_retention_module = @import("../canopy/energy/precipitation_retention.zig");
const plant_root_module = @import("../plant/root/plant_root_system.zig");

pub const Parameters = struct {
    snow_ice_density_megagrams_per_m3: f64,
    /// HEAT-001 resolution A. Latent heat of fusion used to convert the snow
    /// owner's frozen carriers from sensible heat to enthalpy. Authoritative
    /// source is `runscript.snow_latent_heat_of_fusion_megajoules_per_m3`.
    ///
    /// The default exists only because the composition root in
    /// `src/ecosys_ng.zig` is owned by the Integrator lane and cannot be
    /// edited here; see
    /// `docs/binding_requests/heat_001_landscape_enthalpy.md`. It equals the
    /// value every shipped runscript currently carries. Once the binding
    /// request lands the default must be deleted so a runscript that changes
    /// the constant cannot silently disagree with the census.
    snow_latent_heat_of_fusion_megajoules_per_m3: f64 = 333,
    /// As above, for the soil matrix and macropore ice carriers.
    /// Authoritative source is
    /// `runscript.soil_phase_heat_parameters.freeze_thaw.latent_heat_of_fusion_megajoules_per_m3`.
    soil_latent_heat_of_fusion_megajoules_per_m3: f64 = 333,
    carbon_g_per_mol: f64,
    nitrogen_g_per_mol: f64,
    phosphorus_g_per_mol: f64,
    phosphate_zone_fractions: zone_module.ZoneFractions,
    surface_physical: inventory.SurfacePhysicalParameters,
};

/// All authoritative owners needed by the seven EXEC storage equations.
/// `soil_mass_megagrams_scratch` is caller-owned runtime memory so reconstruction
/// performs no allocation and remains suitable at every daily audit boundary.
pub const Inputs = struct {
    grid: *const grid_module.GridState,
    plants: *const grid_module.PlantState,
    snow: *const snow_module.State,
    soil_thermal: *const thermal_module.State,
    soil_properties: *const soil_properties_module.State,
    soil_gas: *const gas_module.State,
    root_gas: ?*const plant_root_module.State,
    soil_organic: *const organic_module.State,
    soil_organic_transport: *const organic_transport_module.State,
    surface_organic: *const organic_module.State,
    mineral_nitrogen: *const mineral_nitrogen_module.State,
    soil_chemistry: *const chemistry_module.State,
    nitrogen_fertilizer: *const nitrogen_fertilizer_module.State,
    mineral_fertilizer: *const mineral_fertilizer_module.State,
    micropore_solutes: *const solute_module.State,
    macropore_solutes: *const solute_module.State,
    surface_chemistry: *const litter_chemistry_module.State,
    surface_fertilizer: *const litter_fertilizer_module.State,
    surface_denitrification_nitrite_g_n: []const f64,
    surface: *const surface_module.RuntimeState,
    surface_ice_water_equivalent_m3: []const f64,
    surface_gas: *const gas_module.State,
    surface_litter_dry_mass_megagrams: []const f64,
    canopy_retention: ?*const canopy_retention_module.State,
    cell_area_m2: []const f64,
    soil_mass_megagrams_scratch: []f64,
    parameters: Parameters,
};

/// Reconstructs one complete, finite EXEC snapshot without mutating any
/// scientific owner. Boundary history is published only after all storage
/// contributors have validated, preventing a partially assembled audit.
pub fn reconstruct(
    inputs: Inputs,
    boundary_ledger: *const boundary.State,
) !audit.Totals {
    try validateDimensions(inputs);
    var storage: inventory.Storage = .{};
    try storage.add(try inventory.aggregateSnowEnthalpy(
        inputs.snow,
        inputs.parameters.snow_ice_density_megagrams_per_m3,
        inputs.parameters.snow_latent_heat_of_fusion_megajoules_per_m3,
        // HEAT-001 second layer. The snow census re-bases its frozen carriers
        // onto the same ice-branch enthalpy definition the soil and surface
        // carriers use, so it needs both heat capacities. They come from the
        // same authoritative struct the soil aggregator reads below, which is
        // what keeps snow, surface, and soil on one definition.
        inputs.parameters.surface_physical.liquid_water_heat_capacity_megajoules_per_m3_k,
        inputs.parameters.surface_physical.ice_heat_capacity_megajoules_per_m3_k,
    ));
    try storage.add(try inventory.aggregateSoilPhysicalAndGas(
        inputs.grid,
        inputs.soil_thermal.dry_solid_heat_capacity_megajoules_per_m3_k,
        inputs.soil_thermal.layer_volume_m3,
        inputs.parameters.surface_physical.liquid_water_heat_capacity_megajoules_per_m3_k,
        inputs.parameters.surface_physical.ice_heat_capacity_megajoules_per_m3_k,
        inputs.parameters.soil_latent_heat_of_fusion_megajoules_per_m3,
        inputs.soil_gas,
    ));
    if (inputs.root_gas) |root_gas|
        try storage.add(try inventory.aggregateRootGas(root_gas));
    try storage.add(try inventory.aggregateSurfaceOrganic(
        inputs.surface_organic,
    ));
    try storage.add(try inventory.aggregateSoilOrganic(
        inputs.soil_organic,
        inputs.grid,
    ));
    try storage.add(try inventory.aggregateSoilOrganicTransportMacropore(
        inputs.soil_organic_transport,
        inputs.grid,
    ));
    try fillSoilMass(inputs);
    try storage.add(try inventory.aggregateProfileMineralNitrogen(
        inputs.grid,
        inputs.mineral_nitrogen,
        inputs.soil_chemistry,
        inputs.nitrogen_fertilizer,
        inputs.soil_mass_megagrams_scratch,
        inputs.parameters.nitrogen_g_per_mol,
    ));
    try storage.add(try inventory.aggregateProfilePhosphorusAndIons(
        inputs.grid,
        inputs.micropore_solutes,
        inputs.macropore_solutes,
        inputs.soil_chemistry,
        inputs.mineral_fertilizer,
        inputs.grid.matrix_liquid_water_m3,
        inputs.soil_mass_megagrams_scratch,
        inputs.parameters.phosphate_zone_fractions,
        inputs.parameters.carbon_g_per_mol,
        inputs.parameters.phosphorus_g_per_mol,
    ));
    try storage.add(try inventory.aggregatePendingSurfaceMinerals(
        inputs.mineral_fertilizer,
        inputs.parameters.carbon_g_per_mol,
        inputs.parameters.phosphorus_g_per_mol,
    ));
    try storage.add(try inventory.aggregateSurfaceChemistry(
        inputs.surface_chemistry,
        inputs.surface_fertilizer,
        inputs.surface_denitrification_nitrite_g_n,
        inputs.surface.litter_water_m3,
        inputs.surface_litter_dry_mass_megagrams,
        inputs.parameters.carbon_g_per_mol,
        inputs.parameters.nitrogen_g_per_mol,
        inputs.parameters.phosphorus_g_per_mol,
    ));
    try storage.add(try inventory.aggregateSurfacePhysicalAndGas(
        inputs.surface,
        inputs.surface_ice_water_equivalent_m3,
        inputs.grid,
        inputs.surface_gas,
        inputs.surface_organic,
        inputs.parameters.surface_physical,
    ));
    if (inputs.canopy_retention) |retention|
        try storage.add(try inventory.aggregateCanopyWaterAndHeat(
            inputs.plants,
            retention,
            inputs.cell_area_m2,
        ));

    var totals = std.mem.zeroes(audit.Totals);
    for (inputs.cell_area_m2) |area_m2| {
        if (!std.math.isFinite(area_m2) or area_m2 <= 0)
            return error.InvalidLandscapeCellArea;
        totals.landscape_area_m2 += area_m2;
        if (!std.math.isFinite(totals.landscape_area_m2))
            return error.LandscapeAreaOverflow;
    }
    try inventory.publishStorage(&totals, storage);
    try boundary_ledger.publish(&totals);
    std.log.debug("heat balance instrument: storage={e} cumulative_in={e} cumulative_out={e} balance={e}", .{
        totals.heat_storage_megajoules,
        totals.cumulative_heat_input_megajoules,
        totals.cumulative_heat_output_megajoules,
        totals.heat_storage_megajoules - totals.cumulative_heat_input_megajoules +
            totals.cumulative_heat_output_megajoules,
    });
    _ = try audit.balance(totals);
    return totals;
}

fn fillSoilMass(inputs: Inputs) !void {
    try deriveSoilMass(
        inputs.soil_properties.matrix_bulk_volume_m3,
        inputs.soil_properties.bulk_density_megagrams_per_m3,
        inputs.soil_mass_megagrams_scratch,
    );
}

pub fn deriveSoilMass(
    matrix_bulk_volume_m3: []const f64,
    bulk_density_megagrams_per_m3: []const f64,
    soil_mass_megagrams: []f64,
) !void {
    if (matrix_bulk_volume_m3.len == 0 or
        bulk_density_megagrams_per_m3.len != matrix_bulk_volume_m3.len or
        soil_mass_megagrams.len != matrix_bulk_volume_m3.len)
        return error.RuntimeSoilMassDimensionMismatch;
    for (matrix_bulk_volume_m3, bulk_density_megagrams_per_m3) |volume_m3, density_megagrams_per_m3| {
        if (!std.math.isFinite(volume_m3) or volume_m3 < 0 or
            !std.math.isFinite(density_megagrams_per_m3) or density_megagrams_per_m3 < 0)
            return error.InvalidRuntimeSoilMassInput;
        const mass_megagrams = volume_m3 * density_megagrams_per_m3;
        if (!std.math.isFinite(mass_megagrams))
            return error.RuntimeSoilMassOverflow;
    }
    for (
        matrix_bulk_volume_m3,
        bulk_density_megagrams_per_m3,
        soil_mass_megagrams,
    ) |volume_m3, density_megagrams_per_m3, *mass_megagrams| {
        mass_megagrams.* = volume_m3 * density_megagrams_per_m3;
    }
}

fn validateDimensions(inputs: Inputs) !void {
    const layers = inputs.grid.layer_count;
    const cells = inputs.grid.cell_count;
    if (layers == 0 or cells == 0 or
        inputs.soil_thermal.total_heat_capacity_megajoules_per_m3_k.len != layers or
        inputs.soil_thermal.layer_volume_m3.len != layers or
        inputs.soil_properties.matrix_bulk_volume_m3.len != layers or
        inputs.soil_properties.bulk_density_megagrams_per_m3.len != layers or
        inputs.soil_mass_megagrams_scratch.len != layers or
        inputs.surface_ice_water_equivalent_m3.len != cells or
        inputs.surface_litter_dry_mass_megagrams.len != cells or
        inputs.cell_area_m2.len != cells)
        return error.LandscapeMassBalanceRuntimeDimensionMismatch;
}

test "runtime soil mass derivation is explicit and rejects late invalid input" {
    var scratch = [_]f64{ 9, 9 };
    try std.testing.expectError(
        error.InvalidRuntimeSoilMassInput,
        deriveSoilMass(
            &.{ 2, 3 },
            &.{ 1.25, std.math.nan(f64) },
            &scratch,
        ),
    );
    try std.testing.expectEqualSlices(f64, &.{ 9, 9 }, &scratch);
    try deriveSoilMass(&.{ 2, 3 }, &.{ 1.25, 1.5 }, &scratch);
    try std.testing.expectEqualSlices(f64, &.{ 2.5, 4.5 }, &scratch);
}
