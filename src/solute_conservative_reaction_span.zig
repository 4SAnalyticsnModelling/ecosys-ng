const std = @import("std");
const chemistry = @import("solute_chemistry_state.zig");
const aqueous_network = @import("solute_aqueous_network.zig");
const aqueous_rates = @import("solute_aqueous_reaction_rates.zig");
const phosphate_network = @import("solute_phosphate_network.zig");
const phosphate_exchange = @import("solute_phosphate_exchange.zig");
const phosphate_rates = @import("solute_phosphate_reaction_rates.zig");
const cation_exchange = @import("solute_cation_exchange.zig");
const geochemistry = @import("solute_geochemistry_network.zig");
const geochemistry_rates = @import("solute_geochemistry_reaction_rates.zig");

pub const aqueous_reaction_count =
    @typeInfo(aqueous_network.Fluxes).@"struct".fields.len;
pub const phosphate_mineral_reaction_count =
    @typeInfo(phosphate_network.MineralFluxes).@"struct".fields.len;
pub const phosphate_surface_reaction_count =
    @typeInfo(phosphate_exchange.Flux).@"struct".fields.len;
pub const phosphate_aqueous_reaction_count =
    @typeInfo(phosphate_network.DissociationAndPairingFluxes).@"struct".fields.len;
pub const phosphate_zone_reaction_count =
    phosphate_mineral_reaction_count +
    phosphate_surface_reaction_count +
    phosphate_aqueous_reaction_count;
pub const equilibrium_mineral_reaction_count =
    @typeInfo(geochemistry.MineralExtents).@"struct".fields.len;

pub const aqueous_reaction_offset: usize = 0;
pub const non_band_phosphate_reaction_offset =
    aqueous_reaction_offset + aqueous_reaction_count;
pub const non_band_phosphate_mineral_offset =
    non_band_phosphate_reaction_offset;
pub const non_band_phosphate_surface_offset =
    non_band_phosphate_mineral_offset + phosphate_mineral_reaction_count;
pub const non_band_phosphate_aqueous_offset =
    non_band_phosphate_surface_offset + phosphate_surface_reaction_count;
pub const band_phosphate_reaction_offset =
    non_band_phosphate_reaction_offset + phosphate_zone_reaction_count;
pub const band_phosphate_mineral_offset = band_phosphate_reaction_offset;
pub const band_phosphate_surface_offset =
    band_phosphate_mineral_offset + phosphate_mineral_reaction_count;
pub const band_phosphate_aqueous_offset =
    band_phosphate_surface_offset + phosphate_surface_reaction_count;
pub const gapon_reaction_index =
    band_phosphate_reaction_offset + phosphate_zone_reaction_count;
pub const carboxyl_reaction_index = gapon_reaction_index + 1;
pub const equilibrium_mineral_reaction_offset = carboxyl_reaction_index + 1;

/// Complete conservative equilibrium-reaction axis count. Operator-split
/// silicate weathering is intentionally excluded from the equilibrium span.
pub const reaction_count =
    equilibrium_mineral_reaction_offset + equilibrium_mineral_reaction_count;

pub const ReactionDomain = enum {
    aqueous,
    non_band_phosphate_mineral,
    non_band_phosphate_surface,
    non_band_phosphate_aqueous,
    band_phosphate_mineral,
    band_phosphate_surface,
    band_phosphate_aqueous,
    cation_exchange,
    carboxyl,
    equilibrium_mineral,
};

pub const ReactionIdentity = struct {
    domain: ReactionDomain,
    name: []const u8,
};

/// Returns the scientific ledger identity for one reaction-span column.
pub fn reactionIdentity(column: usize) ?ReactionIdentity {
    inline for (
        @typeInfo(aqueous_network.Fluxes).@"struct".fields,
        0..,
    ) |field, index| {
        if (column == aqueous_reaction_offset + index)
            return .{ .domain = .aqueous, .name = field.name };
    }
    inline for (
        @typeInfo(phosphate_network.MineralFluxes).@"struct".fields,
        0..,
    ) |field, index| {
        if (column == non_band_phosphate_mineral_offset + index)
            return .{
                .domain = .non_band_phosphate_mineral,
                .name = field.name,
            };
        if (column == band_phosphate_mineral_offset + index)
            return .{
                .domain = .band_phosphate_mineral,
                .name = field.name,
            };
    }
    inline for (
        @typeInfo(phosphate_exchange.Flux).@"struct".fields,
        0..,
    ) |field, index| {
        if (column == non_band_phosphate_surface_offset + index)
            return .{
                .domain = .non_band_phosphate_surface,
                .name = field.name,
            };
        if (column == band_phosphate_surface_offset + index)
            return .{
                .domain = .band_phosphate_surface,
                .name = field.name,
            };
    }
    inline for (
        @typeInfo(phosphate_network.DissociationAndPairingFluxes).@"struct".fields,
        0..,
    ) |field, index| {
        if (column == non_band_phosphate_aqueous_offset + index)
            return .{
                .domain = .non_band_phosphate_aqueous,
                .name = field.name,
            };
        if (column == band_phosphate_aqueous_offset + index)
            return .{
                .domain = .band_phosphate_aqueous,
                .name = field.name,
            };
    }
    if (column == gapon_reaction_index)
        return .{
            .domain = .cation_exchange,
            .name = "gapon_charge_exchange",
        };
    if (column == carboxyl_reaction_index)
        return .{
            .domain = .carboxyl,
            .name = "hydrogen_protonation",
        };
    inline for (
        @typeInfo(geochemistry.MineralExtents).@"struct".fields,
        0..,
    ) |field, index| {
        if (column == equilibrium_mineral_reaction_offset + index)
            return .{
                .domain = .equilibrium_mineral,
                .name = field.name,
            };
    }
    return null;
}

/// Evaluates each native equilibrium rate without modifying `state`.
///
/// The caller must first project the cell's H+/OH- pair onto the configured
/// water equilibrium. This routine deliberately does not apply a second water
/// projection, so every rate is evaluated at exactly the solver's current
/// iterate. Phosphate surface rates are mol/Mg, other phosphate and aqueous
/// rates are mol/m3, the carboxyl rate is mol/Mg, and mineral rates are
/// mol/m3. The Gapon entry is one when its complete conservative vector is
/// active and zero otherwise; its native extent is therefore dimensionless.
pub fn evaluateRates(
    state: *const chemistry.State,
    cell_index: usize,
    parameters: chemistry.ReactionParameters,
    output: []f64,
) !void {
    if (output.len != reaction_count)
        return error.SoluteReactionRateVectorSizeMismatch;
    if (cell_index >= state.cell_count)
        return error.ChemistryCellIndexOutOfBounds;

    const coefficients =
        try state.activityCoefficients(cell_index, parameters.fractions);
    const aqueous = try aqueous_rates.calculate(
        state.aqueous[cell_index],
        coefficients,
        parameters.aqueous_constants,
        parameters.aqueous_kinetics,
    );
    writeStruct(aqueous_network.Fluxes, aqueous, output, aqueous_reaction_offset);

    try evaluatePhosphateZoneRates(
        state.aqueous[cell_index],
        state.non_band_phosphate[cell_index],
        coefficients,
        parameters.fractions.phosphate_non_band,
        parameters.non_band_phosphate_soil_mass_per_water_volume_Mg_per_m3,
        parameters,
        output[non_band_phosphate_reaction_offset..band_phosphate_reaction_offset],
    );
    try evaluatePhosphateZoneRates(
        state.aqueous[cell_index],
        state.band_phosphate[cell_index],
        coefficients,
        parameters.fractions.phosphate_band,
        parameters.band_phosphate_soil_mass_per_water_volume_Mg_per_m3,
        parameters,
        output[band_phosphate_reaction_offset..gapon_reaction_index],
    );

    const adsorption = try evaluateCationAdsorption(
        state,
        cell_index,
        coefficients,
        parameters,
    );
    output[gapon_reaction_index] =
        if (structHasNonzero(cation_exchange.Cations, adsorption)) 1 else 0;
    output[carboxyl_reaction_index] =
        try state.evaluateCarboxylHydrogenChange(
            cell_index,
            parameters.total_carboxyl_sites_mol_per_Mg,
            state.aqueous[cell_index].hydrogen *
                coefficients.monovalent_activity_coefficient,
            parameters.cation_exchange_water_ratios.shared_Mg_per_m3,
            parameters.carboxyl_exchange_parameters,
        );

    var equilibrium_kinetics = parameters.geochemistry_kinetics;
    equilibrium_kinetics.maximum_natural_weathering_mol_per_m3_step = 0;
    equilibrium_kinetics.maximum_ground_weathering_mol_per_m3_step = 0;
    const minerals = try geochemistry_rates.calculate(
        state.aqueous[cell_index],
        state.geochemistry_solids[cell_index],
        coefficients,
        parameters.geochemistry_products,
        equilibrium_kinetics,
    );
    output[equilibrium_mineral_reaction_offset + 0] =
        minerals.gibbsite_solid_mol_per_m3;
    output[equilibrium_mineral_reaction_offset + 1] =
        minerals.iron_hydroxide_solid_mol_per_m3;
    output[equilibrium_mineral_reaction_offset + 2] =
        minerals.calcite_solid_mol_per_m3;
    output[equilibrium_mineral_reaction_offset + 3] =
        minerals.gypsum_solid_mol_per_m3;
}

/// Creates an empty cell ledger with all runtime conversion metadata set.
pub fn zeroTransformations(
    parameters: chemistry.ReactionParameters,
) chemistry.CellTransformations {
    var result = std.mem.zeroes(chemistry.CellTransformations);
    result.non_band_phosphate_water_fraction =
        parameters.fractions.phosphate_non_band;
    result.band_phosphate_water_fraction =
        parameters.fractions.phosphate_band;
    result.cation_exchange_water_ratios =
        parameters.cation_exchange_water_ratios;
    result.carboxyl_soil_mass_per_water_volume_Mg_per_m3 =
        parameters.cation_exchange_water_ratios.shared_Mg_per_m3;
    return result;
}

/// Adds one conservative reaction axis to `target`.
///
/// `native_extent` uses the rate units documented by `evaluateRates`. For the
/// Gapon axis it is a dimensionless multiplier of the complete current
/// adsorption vector. Scaling that vector as one object is essential: its
/// charge-weighted sum is conservative, while individual cation fields are
/// not conservative reaction directions.
pub fn addReactionExtent(
    target: *chemistry.CellTransformations,
    column: usize,
    native_extent: f64,
    current_transformations: chemistry.CellTransformations,
    parameters: chemistry.ReactionParameters,
) !void {
    if (column >= reaction_count) return error.SoluteReactionColumnOutOfBounds;
    if (!std.math.isFinite(native_extent))
        return error.NonFiniteSoluteReactionExtent;

    if (column < non_band_phosphate_reaction_offset) {
        var fluxes = std.mem.zeroes(aqueous_network.Fluxes);
        setStructField(
            aqueous_network.Fluxes,
            &fluxes,
            column - aqueous_reaction_offset,
            native_extent,
        );
        const changes = try aqueous_network.assemble(fluxes, .{
            .non_band = parameters.fractions.ammonium_non_band,
            .band = parameters.fractions.ammonium_band,
        });
        try addStruct(aqueous_network.Transformations, &target.aqueous, changes);
        return;
    }
    if (column < band_phosphate_reaction_offset) {
        try addPhosphateExtent(
            &target.non_band_phosphate,
            column - non_band_phosphate_reaction_offset,
            native_extent,
            parameters.non_band_phosphate_soil_mass_per_water_volume_Mg_per_m3,
        );
        return;
    }
    if (column < gapon_reaction_index) {
        try addPhosphateExtent(
            &target.band_phosphate,
            column - band_phosphate_reaction_offset,
            native_extent,
            parameters.band_phosphate_soil_mass_per_water_volume_Mg_per_m3,
        );
        return;
    }
    if (column == gapon_reaction_index) {
        try addScaledCationDirection(
            &target.cation_adsorption_mol_per_Mg,
            current_transformations.cation_adsorption_mol_per_Mg,
            native_extent,
        );
        return;
    }
    if (column == carboxyl_reaction_index) {
        const next =
            target.carboxyl_hydrogen_change_mol_per_Mg + native_extent;
        if (!std.math.isFinite(next))
            return error.NonFiniteSoluteReactionTransformation;
        target.carboxyl_hydrogen_change_mol_per_Mg = next;
        return;
    }

    var extents = std.mem.zeroes(geochemistry.MineralExtents);
    setStructField(
        geochemistry.MineralExtents,
        &extents,
        column - equilibrium_mineral_reaction_offset,
        native_extent,
    );
    const changes = try geochemistry.assemble(
        extents,
        std.mem.zeroes(geochemistry.WeatheringExtents),
    );
    try addStruct(
        geochemistry.Transformations,
        &target.geochemistry,
        changes,
    );
}

fn evaluatePhosphateZoneRates(
    shared: aqueous_network.State,
    zone: phosphate_network.State,
    coefficients: @import("solute_activity_coefficients.zig").Result,
    water_fraction: f64,
    soil_mass_per_water_volume_Mg_per_m3: f64,
    parameters: chemistry.ReactionParameters,
    output: []f64,
) !void {
    if (output.len != phosphate_zone_reaction_count)
        return error.SolutePhosphateRateVectorSizeMismatch;
    if (water_fraction == 0) {
        @memset(output, 0);
        return;
    }
    const fluxes = try phosphate_rates.calculate(
        shared,
        zone,
        coefficients,
        soil_mass_per_water_volume_Mg_per_m3,
        parameters.phosphate_constants,
        parameters.phosphate_surface,
        parameters.phosphate_minerals,
        parameters.phosphate_kinetics,
    );
    writeStruct(phosphate_network.MineralFluxes, fluxes.minerals, output, 0);
    writeStruct(
        phosphate_exchange.Flux,
        fluxes.surface,
        output,
        phosphate_mineral_reaction_count,
    );
    writeStruct(
        phosphate_network.DissociationAndPairingFluxes,
        fluxes.aqueous,
        output,
        phosphate_mineral_reaction_count + phosphate_surface_reaction_count,
    );
}

fn evaluateCationAdsorption(
    state: *const chemistry.State,
    cell_index: usize,
    coefficients: @import("solute_activity_coefficients.zig").Result,
    parameters: chemistry.ReactionParameters,
) !cation_exchange.Cations {
    const shared = state.aqueous[cell_index];
    const concentrations = cation_exchange.Cations{
        .ammonium_non_band = shared.ammonium_non_band,
        .ammonium_band = shared.ammonium_band,
        .hydrogen = shared.hydrogen,
        .aluminum = shared.aluminum,
        .iron = shared.iron,
        .calcium = shared.calcium,
        .magnesium = shared.magnesium,
        .sodium = shared.sodium,
        .potassium = shared.potassium,
    };
    var activities = concentrations;
    activities.ammonium_non_band *= coefficients.monovalent_activity_coefficient;
    activities.ammonium_band *= coefficients.monovalent_activity_coefficient;
    activities.hydrogen *= coefficients.monovalent_activity_coefficient;
    activities.aluminum *= coefficients.trivalent_activity_coefficient;
    activities.iron *= coefficients.trivalent_activity_coefficient;
    activities.calcium *= coefficients.divalent_activity_coefficient;
    activities.magnesium *= coefficients.divalent_activity_coefficient;
    activities.sodium *= coefficients.monovalent_activity_coefficient;
    activities.potassium *= coefficients.monovalent_activity_coefficient;
    return cation_exchange.calculateSourceOrder(.{
        .cation_exchange_capacity_mol_charge_per_Mg = parameters.cation_exchange_capacity_mol_charge_per_Mg,
        .aqueous_concentration_mol_per_m3 = concentrations,
        .aqueous_activity_mol_per_m3 = activities,
        .exchange_concentration_mol_per_Mg = state.cation_exchange_mol_per_Mg[cell_index],
        .ammonium_non_band_fraction = parameters.fractions.ammonium_non_band,
        .ammonium_band_fraction = parameters.fractions.ammonium_band,
        .soil_mass_per_water_volume_Mg_per_m3 = parameters.cation_exchange_water_ratios.shared_Mg_per_m3,
    }, parameters.cation_exchange_parameters, .{
        .minimum_activity_mol_per_m3 = parameters.negligible_water_ion_concentration_mol_per_m3,
    });
}

fn addPhosphateExtent(
    target: *phosphate_network.Transformations,
    local_column: usize,
    native_extent: f64,
    soil_mass_per_water_volume_Mg_per_m3: f64,
) !void {
    var fluxes = std.mem.zeroes(phosphate_network.Fluxes);
    fluxes.soil_mass_per_water_volume_Mg_per_m3 =
        soil_mass_per_water_volume_Mg_per_m3;
    if (local_column < phosphate_mineral_reaction_count) {
        setStructField(
            phosphate_network.MineralFluxes,
            &fluxes.minerals,
            local_column,
            native_extent,
        );
    } else if (local_column <
        phosphate_mineral_reaction_count + phosphate_surface_reaction_count)
    {
        setStructField(
            phosphate_exchange.Flux,
            &fluxes.surface,
            local_column - phosphate_mineral_reaction_count,
            native_extent,
        );
    } else {
        setStructField(
            phosphate_network.DissociationAndPairingFluxes,
            &fluxes.aqueous,
            local_column -
                phosphate_mineral_reaction_count -
                phosphate_surface_reaction_count,
            native_extent,
        );
    }
    const changes = try phosphate_network.assemble(fluxes);
    try addStruct(phosphate_network.Transformations, target, changes);
}

fn addScaledCationDirection(
    target: *cation_exchange.Cations,
    current: cation_exchange.Cations,
    multiplier: f64,
) !void {
    const charge = current.ammonium_non_band +
        current.ammonium_band +
        current.hydrogen +
        3 * current.aluminum +
        3 * current.iron +
        2 * current.calcium +
        2 * current.magnesium +
        current.sodium +
        current.potassium;
    const magnitude = @abs(current.ammonium_non_band) +
        @abs(current.ammonium_band) +
        @abs(current.hydrogen) +
        3 * @abs(current.aluminum) +
        3 * @abs(current.iron) +
        2 * @abs(current.calcium) +
        2 * @abs(current.magnesium) +
        @abs(current.sodium) +
        @abs(current.potassium);
    if (!std.math.isFinite(charge) or
        @abs(charge) > 128 * std.math.floatEps(f64) * @max(1.0, magnitude))
    {
        return error.NonConservativeCationExchangeDirection;
    }
    var scaled = current;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        @field(scaled, field.name) *= multiplier;
    try addStruct(cation_exchange.Cations, target, scaled);
}

fn writeStruct(
    comptime T: type,
    value: T,
    output: []f64,
    offset: usize,
) void {
    inline for (@typeInfo(T).@"struct".fields, 0..) |field, index|
        output[offset + index] = @field(value, field.name);
}

fn setStructField(
    comptime T: type,
    target: *T,
    index: usize,
    value: f64,
) void {
    inline for (@typeInfo(T).@"struct".fields, 0..) |field, field_index| {
        if (index == field_index) @field(target.*, field.name) = value;
    }
}

fn addStruct(comptime T: type, target: *T, change: T) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        const next = @field(target.*, field.name) + @field(change, field.name);
        if (!std.math.isFinite(next))
            return error.NonFiniteSoluteReactionTransformation;
    }
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(target.*, field.name) += @field(change, field.name);
}

fn structHasNonzero(comptime T: type, value: T) bool {
    inline for (@typeInfo(T).@"struct".fields) |field|
        if (@field(value, field.name) != 0) return true;
    return false;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn testParameters() chemistry.ReactionParameters {
    return .{
        .fractions = .{
            .ammonium_non_band = 0.8,
            .ammonium_band = 0.2,
            .nitrate_non_band = 0.6,
            .nitrate_band = 0.4,
            .phosphate_non_band = 0.7,
            .phosphate_band = 0.3,
        },
        .non_band_phosphate_soil_mass_per_water_volume_Mg_per_m3 = 1.2,
        .band_phosphate_soil_mass_per_water_volume_Mg_per_m3 = 1.5,
        .cation_exchange_capacity_mol_charge_per_Mg = 10,
        .cation_exchange_water_ratios = .{
            .shared_Mg_per_m3 = 1.4,
            .ammonium_non_band_Mg_per_m3 = 1.1,
            .ammonium_band_Mg_per_m3 = 1.8,
        },
        .total_carboxyl_sites_mol_per_Mg = 2,
        .carboxyl_exchange_parameters = .{
            .dissociation_constant_mol_per_m3 = 0.01,
            .maximum_exchange_mol_per_m3_per_iteration = 0.01,
            .substrate_limit_fraction_per_iteration = 0.2,
        },
        .aqueous_constants = filled(aqueous_rates.EquilibriumConstants, 1),
        .aqueous_kinetics = .{
            .ammonium_substrate_limit_fraction = 0.2,
            .general_substrate_limit_fraction = 0.2,
            .maximum_fast_association_mol_per_m3_step = 0.01,
            .maximum_slow_association_mol_per_m3_step = 0.01,
        },
        .phosphate_constants = filled(phosphate_rates.EquilibriumConstants, 1),
        .phosphate_surface = .{
            .protonated_site_equilibrium_constant = 1,
            .hydroxyl_site_equilibrium_constant = 1,
            .h2po4_exchange_equilibrium_constant = 1,
            .hpo4_exchange_equilibrium_constant = 1,
            .water_activity_product_mol2_per_m6 = 1,
            .h2po4_dissociation_constant = 1,
            .maximum_exchange_mol_per_Mg_step = 0.01,
            .substrate_limit_fraction = 0.2,
        },
        .phosphate_minerals = filled(phosphate_rates.MineralParameters, 1),
        .phosphate_kinetics = .{
            .substrate_limit_fraction = 0.2,
            .maximum_pairing_mol_per_m3_step = 0.01,
        },
        .cation_exchange_parameters = .{
            .selectivity = .{
                .calcium_ammonium = 1,
                .calcium_hydrogen = 1,
                .calcium_aluminum_and_iron = 1,
                .calcium_magnesium = 1,
                .calcium_sodium = 1,
                .calcium_potassium = 1,
            },
            .substrate_limit_fraction = 0.2,
            .maximum_adsorption_mol_charge_per_m3_step = 0.01,
        },
        .geochemistry_products = filled(geochemistry_rates.SolubilityProducts, 1),
        .geochemistry_kinetics = .{
            .general_substrate_limit_fraction = 0.2,
            .hydrogen_coupled_substrate_limit_fraction = 0.2,
            .maximum_hydroxide_mineral_mol_per_m3_step = 0.01,
            .maximum_general_mineral_mol_per_m3_step = 0.01,
            .calcite_hydroxide_inhibition_constant_mol_per_m3 = 1,
            .maximum_natural_weathering_mol_per_m3_step = 0.01,
            .maximum_ground_weathering_mol_per_m3_step = 0.01,
        },
        .water_activity_product_mol2_per_m6 = 1,
        .negligible_water_ion_concentration_mol_per_m3 = 1e-32,
    };
}

fn resetPositiveTestState(state: *chemistry.State) void {
    state.aqueous[0] = filled(aqueous_network.State, 10);
    state.non_band_phosphate[0] = filled(phosphate_network.State, 10);
    state.band_phosphate[0] = filled(phosphate_network.State, 10);
    state.cation_exchange_mol_per_Mg[0] =
        filled(cation_exchange.Cations, 1);
    state.carboxyl_bound_hydrogen_mol_per_Mg[0] = 1;
    state.geochemistry_solids[0] = filled(geochemistry.SolidState, 10);
    state.water_mol_per_m3[0] = 10;
}

test "reaction span offsets cover every equilibrium ledger axis" {
    try std.testing.expectEqual(@as(usize, 25), aqueous_reaction_count);
    try std.testing.expectEqual(@as(usize, 5), phosphate_mineral_reaction_count);
    try std.testing.expectEqual(@as(usize, 5), phosphate_surface_reaction_count);
    try std.testing.expectEqual(@as(usize, 9), phosphate_aqueous_reaction_count);
    try std.testing.expectEqual(@as(usize, 25), non_band_phosphate_reaction_offset);
    try std.testing.expectEqual(@as(usize, 44), band_phosphate_reaction_offset);
    try std.testing.expectEqual(@as(usize, 63), gapon_reaction_index);
    try std.testing.expectEqual(@as(usize, 64), carboxyl_reaction_index);
    try std.testing.expectEqual(@as(usize, 65), equilibrium_mineral_reaction_offset);
    try std.testing.expectEqual(@as(usize, 69), reaction_count);
}

test "reaction identities cover every equilibrium ledger axis" {
    for (0..reaction_count) |column|
        try std.testing.expect(reactionIdentity(column) != null);
    try std.testing.expect(reactionIdentity(reaction_count) == null);

    const aqueous = reactionIdentity(aqueous_reaction_offset).?;
    try std.testing.expectEqual(ReactionDomain.aqueous, aqueous.domain);
    try std.testing.expectEqualStrings(
        "ammonium_non_band_association",
        aqueous.name,
    );
    const non_band_mineral =
        reactionIdentity(non_band_phosphate_mineral_offset).?;
    try std.testing.expectEqual(
        ReactionDomain.non_band_phosphate_mineral,
        non_band_mineral.domain,
    );
    try std.testing.expectEqualStrings(
        "aluminum_phosphate_mol_per_m3",
        non_band_mineral.name,
    );
    const band_surface =
        reactionIdentity(band_phosphate_surface_offset).?;
    try std.testing.expectEqual(
        ReactionDomain.band_phosphate_surface,
        band_surface.domain,
    );
    try std.testing.expectEqualStrings(
        "protonated_to_hydroxyl_site_mol_per_Mg",
        band_surface.name,
    );
    const band_aqueous =
        reactionIdentity(band_phosphate_aqueous_offset).?;
    try std.testing.expectEqual(
        ReactionDomain.band_phosphate_aqueous,
        band_aqueous.domain,
    );
    try std.testing.expectEqualStrings(
        "po4_hydrogen_association_mol_p_per_m3",
        band_aqueous.name,
    );
    try std.testing.expectEqualStrings(
        "gapon_charge_exchange",
        reactionIdentity(gapon_reaction_index).?.name,
    );
    try std.testing.expectEqualStrings(
        "hydrogen_protonation",
        reactionIdentity(carboxyl_reaction_index).?.name,
    );
    const final_mineral = reactionIdentity(reaction_count - 1).?;
    try std.testing.expectEqual(
        ReactionDomain.equilibrium_mineral,
        final_mineral.domain,
    );
    try std.testing.expectEqualStrings(
        "gypsum_precipitation_mol_per_m3",
        final_mineral.name,
    );
}

test "runtime rates retain disabled phosphate axes as zero" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    resetPositiveTestState(&state);
    const rates = try std.testing.allocator.alloc(f64, reaction_count);
    defer std.testing.allocator.free(rates);
    try evaluateRates(&state, 0, testParameters(), rates);

    const calcium_po4_local_index =
        phosphate_mineral_reaction_count +
        phosphate_surface_reaction_count + 5;
    try std.testing.expectEqual(
        @as(f64, 0),
        rates[non_band_phosphate_reaction_offset + calcium_po4_local_index],
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        rates[band_phosphate_reaction_offset + calcium_po4_local_index],
    );
    for (rates) |rate| try std.testing.expect(std.math.isFinite(rate));
}

test "every reaction axis assembles an admissible atomic cell transaction" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const parameters = testParameters();
    var current = zeroTransformations(parameters);
    current.cation_adsorption_mol_per_Mg.ammonium_non_band = 0.01;
    current.cation_adsorption_mol_per_Mg.sodium = -0.01;

    for (0..reaction_count) |column| {
        resetPositiveTestState(&state);
        var target = zeroTransformations(parameters);
        const extent: f64 = if (column == gapon_reaction_index)
            0.5
        else
            1e-4;
        try addReactionExtent(
            &target,
            column,
            extent,
            current,
            parameters,
        );
        try state.commitCell(0, target);
    }
}

test "Gapon axis scales the complete charge-conserving vector" {
    const parameters = testParameters();
    var current = zeroTransformations(parameters);
    current.cation_adsorption_mol_per_Mg = .{
        .ammonium_non_band = 0.1,
        .ammonium_band = 0,
        .hydrogen = 0,
        .aluminum = 0,
        .iron = 0,
        .calcium = -0.05,
        .magnesium = 0,
        .sodium = 0,
        .potassium = 0,
    };
    var target = zeroTransformations(parameters);
    try addReactionExtent(
        &target,
        gapon_reaction_index,
        0.25,
        current,
        parameters,
    );
    try std.testing.expectEqual(
        @as(f64, 0.025),
        target.cation_adsorption_mol_per_Mg.ammonium_non_band,
    );
    try std.testing.expectEqual(
        @as(f64, -0.0125),
        target.cation_adsorption_mol_per_Mg.calcium,
    );
    const charge =
        target.cation_adsorption_mol_per_Mg.ammonium_non_band +
        2 * target.cation_adsorption_mol_per_Mg.calcium;
    try std.testing.expectEqual(@as(f64, 0), charge);
}
