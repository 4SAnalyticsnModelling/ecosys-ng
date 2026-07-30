const std = @import("std");
const aqueous_network = @import("solute_aqueous_network.zig");
const phosphate_network = @import("solute_phosphate_network.zig");
const phosphate_exchange = @import("solute_phosphate_exchange.zig");
const ion_pairing = @import("solute_ion_pairing.zig");
const activity_coefficients = @import("solute_activity_coefficients.zig");
const precipitation = @import("solute_phosphate_precipitation.zig");

pub const EquilibriumConstants = struct {
    hpo4: f64,
    h2po4: f64,
    h3po4: f64,
    iron_hpo4: f64,
    iron_h2po4: f64,
    calcium_hpo4: f64,
    calcium_h2po4: f64,
    magnesium_hpo4: f64,
};

pub const Kinetics = struct {
    substrate_limit_fraction: f64,
    maximum_pairing_mol_per_m3_step: f64,
};

pub const MineralParameters = struct {
    aluminum_phosphate_solubility_product: f64,
    iron_phosphate_solubility_product: f64,
    dicalcium_phosphate_solubility_product: f64,
    hydroxyapatite_solubility_product: f64,
    monocalcium_phosphate_solubility_product: f64,
    water_activity_product_mol2_per_m6: f64,
    maximum_phosphate_precipitation_mol_per_m3_step: f64,
    maximum_apatite_precipitation_mol_per_m3_step: f64,
    maximum_mineral_dissolution_mol_per_m3_step: f64,
};

pub const ZoneIdentity = enum {
    non_band,
    band,
};

pub const SourceMineralInputs = struct {
    shared: aqueous_network.State,
    zone: phosphate_network.State,
    zone_identity: ZoneIdentity,
    coefficients: activity_coefficients.Result,
    constants: EquilibriumConstants,
    parameters: MineralParameters,
    substrate_limit_fraction: f64,
};

pub fn calculate(shared: aqueous_network.State, zone: phosphate_network.State, coefficients: activity_coefficients.Result, soil_mass_per_water_volume_Mg_per_m3: f64, constants: EquilibriumConstants, surface_parameters: phosphate_exchange.Parameters, mineral_parameters: ?MineralParameters, kinetics: Kinetics) !phosphate_network.Fluxes {
    try validate(shared, zone, coefficients, soil_mass_per_water_volume_Mg_per_m3, constants, kinetics);
    const g1 = coefficients.monovalent_activity_coefficient;
    const g2 = coefficients.divalent_activity_coefficient;
    const g3 = coefficients.trivalent_activity_coefficient;
    const limit = kinetics.substrate_limit_fraction;
    const maximum = kinetics.maximum_pairing_mol_per_m3_step;
    // Direct shared mapping for SOLUTE.F 1709--1893; the caller owns the
    // source `VOLWPO`/`VOLWPB` wet-zone gate.
    const aqueous = phosphate_network.DissociationAndPairingFluxes{
        .po4_hydrogen_association_mol_p_per_m3 = try reaction(zone.dissolved_po4_mol_p_per_m3, shared.hydrogen, zone.dissolved_hpo4_mol_p_per_m3, zone.dissolved_po4_mol_p_per_m3 * g3, shared.hydrogen * g1, zone.dissolved_hpo4_mol_p_per_m3 * g2, g3, constants.hpo4, limit, maximum),
        .hpo4_hydrogen_association_mol_p_per_m3 = try reaction(zone.dissolved_hpo4_mol_p_per_m3, shared.hydrogen, zone.dissolved_h2po4_mol_p_per_m3, zone.dissolved_hpo4_mol_p_per_m3 * g2, shared.hydrogen * g1, zone.dissolved_h2po4_mol_p_per_m3 * g1, g2, constants.h2po4, limit, maximum),
        .h2po4_hydrogen_association_mol_p_per_m3 = try reaction(zone.dissolved_h2po4_mol_p_per_m3, shared.hydrogen, zone.dissolved_h3po4_mol_p_per_m3, zone.dissolved_h2po4_mol_p_per_m3 * g1, shared.hydrogen * g1, zone.dissolved_h3po4_mol_p_per_m3, g1, constants.h3po4, limit, maximum),
        .iron_hpo4_pairing_mol_p_per_m3 = try reaction(shared.iron, zone.dissolved_hpo4_mol_p_per_m3, zone.iron_hpo4_pair_mol_per_m3, shared.iron * g3, zone.dissolved_hpo4_mol_p_per_m3 * g2, zone.iron_hpo4_pair_mol_per_m3 * g2, g3, constants.iron_hpo4, limit, maximum),
        .iron_h2po4_pairing_mol_p_per_m3 = try reaction(shared.iron, zone.dissolved_h2po4_mol_p_per_m3, zone.iron_h2po4_pair_mol_per_m3, shared.iron * g3, zone.dissolved_h2po4_mol_p_per_m3 * g1, zone.iron_h2po4_pair_mol_per_m3 * g1, g3, constants.iron_h2po4, limit, maximum),
        // CaPO4 pairing is explicitly disabled in both full SOLUTE branches.
        .calcium_po4_pairing_mol_p_per_m3 = 0,
        .calcium_hpo4_pairing_mol_p_per_m3 = try reaction(shared.calcium, zone.dissolved_hpo4_mol_p_per_m3, zone.calcium_hpo4_pair_mol_per_m3, shared.calcium * g2, zone.dissolved_hpo4_mol_p_per_m3 * g2, zone.calcium_hpo4_pair_mol_per_m3, g2, constants.calcium_hpo4, limit, maximum),
        .calcium_h2po4_pairing_mol_p_per_m3 = try reaction(shared.calcium, zone.dissolved_h2po4_mol_p_per_m3, zone.calcium_h2po4_pair_mol_per_m3, shared.calcium * g2, zone.dissolved_h2po4_mol_p_per_m3 * g1, zone.calcium_h2po4_pair_mol_per_m3 * g1, g2, constants.calcium_h2po4, limit, maximum),
        .magnesium_hpo4_pairing_mol_p_per_m3 = try reaction(shared.magnesium, zone.dissolved_hpo4_mol_p_per_m3, zone.magnesium_hpo4_pair_mol_per_m3, shared.magnesium * g2, zone.dissolved_hpo4_mol_p_per_m3 * g2, zone.magnesium_hpo4_pair_mol_per_m3, g2, constants.magnesium_hpo4, limit, maximum),
    };
    const surface = try phosphate_exchange.calculate(.{
        .hydrogen_concentration_mol_per_m3 = shared.hydrogen,
        .hydrogen_activity_mol_per_m3 = shared.hydrogen * g1,
        .hydroxide_activity_mol_per_m3 = shared.hydroxide * g1,
        .h2po4_concentration_mol_p_per_m3 = zone.dissolved_h2po4_mol_p_per_m3,
        .h2po4_activity_mol_p_per_m3 = zone.dissolved_h2po4_mol_p_per_m3 * g1,
        .hpo4_concentration_mol_p_per_m3 = zone.dissolved_hpo4_mol_p_per_m3,
        .hpo4_activity_mol_p_per_m3 = zone.dissolved_hpo4_mol_p_per_m3 * g2,
        .deprotonated_site_mol_per_Mg = zone.deprotonated_site_mol_per_Mg,
        .hydroxyl_site_mol_per_Mg = zone.hydroxyl_site_mol_per_Mg,
        .protonated_site_mol_per_Mg = zone.protonated_site_mol_per_Mg,
        .adsorbed_h2po4_mol_p_per_Mg = zone.adsorbed_h2po4_mol_p_per_Mg,
        .adsorbed_hpo4_mol_p_per_Mg = zone.adsorbed_hpo4_mol_p_per_Mg,
        .monovalent_activity_coefficient = g1,
        .divalent_activity_coefficient = g2,
    }, surface_parameters);
    const minerals = if (mineral_parameters) |parameters| try calculateMinerals(shared, zone, coefficients, constants, parameters, kinetics.substrate_limit_fraction) else phosphate_network.MineralFluxes{ .aluminum_phosphate_mol_per_m3 = 0, .iron_phosphate_mol_per_m3 = 0, .dicalcium_phosphate_mol_per_m3 = 0, .hydroxyapatite_mol_per_m3 = 0, .monocalcium_phosphate_mol_per_m3 = 0 };
    return .{
        .minerals = minerals,
        .surface = surface,
        .aqueous = aqueous,
        .soil_mass_per_water_volume_Mg_per_m3 = soil_mass_per_water_volume_Mg_per_m3,
    };
}

fn calculateMinerals(shared: aqueous_network.State, zone: phosphate_network.State, coefficients: activity_coefficients.Result, constants: EquilibriumConstants, parameters: MineralParameters, substrate_limit_fraction: f64) !phosphate_network.MineralFluxes {
    inline for (@typeInfo(MineralParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name)) or @field(parameters, field.name) <= 0) return error.InvalidPhosphateMineralParameter;
    const g1 = coefficients.monovalent_activity_coefficient;
    const g2 = coefficients.divalent_activity_coefficient;
    const g3 = coefficients.trivalent_activity_coefficient;
    const hydrogen_activity = shared.hydrogen * g1;
    const hydroxide_activity = shared.hydroxide * g1;
    const aluminum_activity = shared.aluminum * g3;
    const iron_activity = shared.iron * g3;
    const calcium_activity = shared.calcium * g2;
    const h2_activity = zone.dissolved_h2po4_mol_p_per_m3 * g1;
    const hpo4_activity = zone.dissolved_hpo4_mol_p_per_m3 * g2;
    const aluminum_target = try precipitation.aluminumOrIronPhosphateEquilibriumH2po4(parameters.aluminum_phosphate_solubility_product, hydrogen_activity, constants.h2po4, constants.hpo4, aluminum_activity);
    const iron_target = try precipitation.aluminumOrIronPhosphateEquilibriumH2po4(parameters.iron_phosphate_solubility_product, hydrogen_activity, constants.h2po4, constants.hpo4, iron_activity);
    const dicalcium_target = try precipitation.dicalciumPhosphateEquilibriumHpo4(parameters.dicalcium_phosphate_solubility_product, calcium_activity);
    const apatite_target = try precipitation.hydroxyapatiteEquilibriumH2po4(parameters.hydroxyapatite_solubility_product, hydrogen_activity, hydroxide_activity, calcium_activity, constants.h2po4, constants.hpo4);
    const monocalcium_target = try precipitation.monocalciumPhosphateEquilibriumH2po4(parameters.monocalcium_phosphate_solubility_product, calcium_activity);
    const standard = precipitation.Kinetics{ .substrate_limit_fraction = substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = parameters.maximum_phosphate_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = parameters.maximum_mineral_dissolution_mol_per_m3_step, .phosphate_activity_coefficient = g1 };
    return .{
        .aluminum_phosphate_mol_per_m3 = try precipitation.calculateExtent(.{ .dissolved_cation_mol_per_m3 = shared.aluminum, .dissolved_phosphate_mol_p_per_m3 = zone.dissolved_h2po4_mol_p_per_m3, .precipitate_mol_per_m3 = zone.aluminum_phosphate_solid_mol_per_m3 }, h2_activity, aluminum_target, .{ .cation_mol_per_mol_precipitate = 1, .phosphorus_mol_per_mol_precipitate = 1 }, standard),
        .iron_phosphate_mol_per_m3 = try precipitation.calculateExtent(.{ .dissolved_cation_mol_per_m3 = shared.iron, .dissolved_phosphate_mol_p_per_m3 = zone.dissolved_h2po4_mol_p_per_m3, .precipitate_mol_per_m3 = zone.iron_phosphate_solid_mol_per_m3 }, h2_activity, iron_target, .{ .cation_mol_per_mol_precipitate = 1, .phosphorus_mol_per_mol_precipitate = 1 }, standard),
        .dicalcium_phosphate_mol_per_m3 = try precipitation.calculateExtent(.{ .dissolved_cation_mol_per_m3 = shared.calcium, .dissolved_phosphate_mol_p_per_m3 = zone.dissolved_hpo4_mol_p_per_m3, .precipitate_mol_per_m3 = zone.dicalcium_phosphate_solid_mol_per_m3 }, hpo4_activity, dicalcium_target, .{ .cation_mol_per_mol_precipitate = 1, .phosphorus_mol_per_mol_precipitate = 1 }, .{ .substrate_limit_fraction = substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = parameters.maximum_phosphate_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = parameters.maximum_mineral_dissolution_mol_per_m3_step, .phosphate_activity_coefficient = g2 }),
        .hydroxyapatite_mol_per_m3 = try precipitation.calculateExtent(.{ .dissolved_cation_mol_per_m3 = shared.calcium, .dissolved_phosphate_mol_p_per_m3 = zone.dissolved_h2po4_mol_p_per_m3, .precipitate_mol_per_m3 = zone.hydroxyapatite_solid_mol_per_m3 }, h2_activity, apatite_target, .{ .cation_mol_per_mol_precipitate = 5, .phosphorus_mol_per_mol_precipitate = 3 }, .{ .substrate_limit_fraction = substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = parameters.maximum_apatite_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = parameters.maximum_mineral_dissolution_mol_per_m3_step, .phosphate_activity_coefficient = g1 }),
        .monocalcium_phosphate_mol_per_m3 = try precipitation.calculateExtent(.{ .dissolved_cation_mol_per_m3 = shared.calcium, .dissolved_phosphate_mol_p_per_m3 = zone.dissolved_h2po4_mol_p_per_m3, .precipitate_mol_per_m3 = zone.monocalcium_phosphate_solid_mol_per_m3 }, h2_activity, monocalcium_target, .{ .cation_mol_per_mol_precipitate = 1, .phosphorus_mol_per_mol_precipitate = 2 }, standard),
    };
}

/// Direct source-order mineral equations for SOLUTE.F lines 1973--2141.
///
/// This entry point retains the source's zone-specific dicalcium limiter,
/// monocalcium P limiter, `0.333` apatite exponent, and per-mineral ceiling
/// choices for controlled comparison. Production remains on the explicitly
/// stoichiometric, finite-safe path in `calculateMinerals`.
pub fn calculateMineralsSourceOrder(
    inputs: SourceMineralInputs,
) !phosphate_network.MineralFluxes {
    try validate(
        inputs.shared,
        inputs.zone,
        inputs.coefficients,
        1,
        inputs.constants,
        .{
            .substrate_limit_fraction = inputs.substrate_limit_fraction,
            .maximum_pairing_mol_per_m3_step = 0,
        },
    );
    inline for (@typeInfo(MineralParameters).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.parameters, field.name)) or
            @field(inputs.parameters, field.name) <= 0)
            return error.InvalidPhosphateMineralParameter;

    const shared = inputs.shared;
    const zone = inputs.zone;
    const p = inputs.parameters;
    const g1 = inputs.coefficients.monovalent_activity_coefficient;
    const g2 = inputs.coefficients.divalent_activity_coefficient;
    const g3 = inputs.coefficients.trivalent_activity_coefficient;
    const hydrogen_activity = shared.hydrogen * g1;
    const hydroxide_activity = shared.hydroxide * g1;
    const aluminum_activity = shared.aluminum * g3;
    const iron_activity = shared.iron * g3;
    const calcium_activity = shared.calcium * g2;
    const h2po4_activity = zone.dissolved_h2po4_mol_p_per_m3 * g1;
    const hpo4_activity = zone.dissolved_hpo4_mol_p_per_m3 * g2;
    const aluminum_target = try precipitation.aluminumOrIronPhosphateEquilibriumH2po4(p.aluminum_phosphate_solubility_product, hydrogen_activity, inputs.constants.h2po4, inputs.constants.hpo4, aluminum_activity);
    const iron_target = try precipitation.aluminumOrIronPhosphateEquilibriumH2po4(p.iron_phosphate_solubility_product, hydrogen_activity, inputs.constants.h2po4, inputs.constants.hpo4, iron_activity);
    const dicalcium_target = try precipitation.dicalciumPhosphateEquilibriumHpo4(p.dicalcium_phosphate_solubility_product, calcium_activity);
    const apatite_target = try sourceApatiteTarget(
        p.hydroxyapatite_solubility_product,
        hydrogen_activity,
        hydroxide_activity,
        calcium_activity,
        inputs.constants.h2po4,
        inputs.constants.hpo4,
    );
    const monocalcium_target = try precipitation.monocalciumPhosphateEquilibriumH2po4(p.monocalcium_phosphate_solubility_product, calcium_activity);
    const fraction = inputs.substrate_limit_fraction;
    const dicalcium_limiter_phosphate = switch (inputs.zone_identity) {
        .non_band => zone.dissolved_hpo4_mol_p_per_m3,
        .band => zone.dissolved_h2po4_mol_p_per_m3,
    };
    return .{
        .aluminum_phosphate_mol_per_m3 = sourceMineralExtent(zone.aluminum_phosphate_solid_mol_per_m3, h2po4_activity, aluminum_target, g1, fraction * @min(shared.aluminum, zone.dissolved_h2po4_mol_p_per_m3), p.maximum_phosphate_precipitation_mol_per_m3_step, p.maximum_phosphate_precipitation_mol_per_m3_step),
        .iron_phosphate_mol_per_m3 = sourceMineralExtent(zone.iron_phosphate_solid_mol_per_m3, h2po4_activity, iron_target, g1, fraction * @min(shared.iron, zone.dissolved_h2po4_mol_p_per_m3), p.maximum_phosphate_precipitation_mol_per_m3_step, p.maximum_phosphate_precipitation_mol_per_m3_step),
        .dicalcium_phosphate_mol_per_m3 = sourceMineralExtent(zone.dicalcium_phosphate_solid_mol_per_m3, hpo4_activity, dicalcium_target, g2, fraction * @min(shared.calcium, dicalcium_limiter_phosphate), p.maximum_phosphate_precipitation_mol_per_m3_step, p.maximum_phosphate_precipitation_mol_per_m3_step),
        .hydroxyapatite_mol_per_m3 = sourceMineralExtent(zone.hydroxyapatite_solid_mol_per_m3, h2po4_activity, apatite_target, g1, fraction * @min(shared.calcium / 5, zone.dissolved_h2po4_mol_p_per_m3 / 3), p.maximum_apatite_precipitation_mol_per_m3_step, p.maximum_apatite_precipitation_mol_per_m3_step),
        .monocalcium_phosphate_mol_per_m3 = sourceMineralExtent(zone.monocalcium_phosphate_solid_mol_per_m3, h2po4_activity, monocalcium_target, g1, fraction * @min(shared.calcium, zone.dissolved_h2po4_mol_p_per_m3), p.maximum_phosphate_precipitation_mol_per_m3_step, p.maximum_mineral_dissolution_mol_per_m3_step),
    };
}

fn sourceMineralExtent(
    solid_mol_per_m3: f64,
    phosphate_activity_mol_p_per_m3: f64,
    equilibrium_activity_mol_p_per_m3: f64,
    phosphate_activity_coefficient: f64,
    precipitation_inventory_limit_mol_per_m3: f64,
    precipitation_maximum_mol_per_m3_step: f64,
    dissolution_maximum_mol_per_m3_step: f64,
) f64 {
    return @max(
        -solid_mol_per_m3,
        -dissolution_maximum_mol_per_m3_step,
        @min(
            precipitation_maximum_mol_per_m3_step,
            precipitation_inventory_limit_mol_per_m3,
            (phosphate_activity_mol_p_per_m3 -
                equilibrium_activity_mol_p_per_m3) /
                phosphate_activity_coefficient,
        ),
    );
}

fn sourceApatiteTarget(
    solubility_product: f64,
    hydrogen_activity: f64,
    hydroxide_activity: f64,
    calcium_activity: f64,
    h2po4_dissociation_constant: f64,
    hpo4_dissociation_constant: f64,
) !f64 {
    if (calcium_activity == 0) return std.math.floatMax(f64);
    const ratio = solubility_product *
        std.math.pow(f64, hydrogen_activity, 6) /
        (std.math.pow(f64, calcium_activity, 5) * hydroxide_activity *
            std.math.pow(
                f64,
                h2po4_dissociation_constant * hpo4_dissociation_constant,
                3,
            ));
    if (!std.math.isFinite(ratio) or ratio < 0)
        return error.InvalidPhosphateEquilibrium;
    const target = std.math.pow(f64, ratio, 0.333);
    if (!std.math.isFinite(target)) return error.InvalidPhosphateEquilibrium;
    return target;
}

fn reaction(free_first: f64, free_second: f64, paired: f64, first_activity: f64, second_activity: f64, paired_activity: f64, first_coefficient: f64, dissociation_constant: f64, limit_fraction: f64, maximum: f64) !f64 {
    return ion_pairing.calculate(.{ .free_first_mol_per_m3 = free_first, .free_second_mol_per_m3 = free_second, .paired_mol_per_m3 = paired }, .{ .free_first_mol_per_m3 = first_activity, .free_second_mol_per_m3 = second_activity, .paired_mol_per_m3 = paired_activity, .free_first_activity_coefficient = first_coefficient }, .{ .dissociation_constant = dissociation_constant, .substrate_limit_fraction = limit_fraction, .maximum_association_mol_per_m3_step = maximum });
}

fn validate(shared: aqueous_network.State, zone: phosphate_network.State, coefficients: activity_coefficients.Result, density: f64, constants: EquilibriumConstants, kinetics: Kinetics) !void {
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| if (!std.math.isFinite(@field(shared, field.name)) or @field(shared, field.name) < 0) return error.InvalidSharedPhosphateState;
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| if (!std.math.isFinite(@field(zone, field.name)) or @field(zone, field.name) < 0) return error.InvalidZonePhosphateState;
    if (!std.math.isFinite(coefficients.monovalent_activity_coefficient) or coefficients.monovalent_activity_coefficient <= 0 or !std.math.isFinite(coefficients.divalent_activity_coefficient) or coefficients.divalent_activity_coefficient <= 0 or !std.math.isFinite(coefficients.trivalent_activity_coefficient) or coefficients.trivalent_activity_coefficient <= 0 or !std.math.isFinite(density) or density <= 0) return error.InvalidPhosphateReactionInput;
    inline for (@typeInfo(EquilibriumConstants).@"struct".fields) |field| if (!std.math.isFinite(@field(constants, field.name)) or @field(constants, field.name) < 0) return error.InvalidPhosphateEquilibriumConstant;
    if (!std.math.isFinite(kinetics.substrate_limit_fraction) or kinetics.substrate_limit_fraction < 0 or kinetics.substrate_limit_fraction > 1 or !std.math.isFinite(kinetics.maximum_pairing_mol_per_m3_step) or kinetics.maximum_pairing_mol_per_m3_step < 0) return error.InvalidPhosphateKinetics;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "phosphate aqueous equilibria feed a phosphorus-conserving ledger" {
    const shared = filled(aqueous_network.State, 1);
    const zone = filled(phosphate_network.State, 1);
    const coefficients = activity_coefficients.Result{ .ionic_strength_mol_per_l = 0, .monovalent_activity_coefficient = 1, .divalent_activity_coefficient = 1, .trivalent_activity_coefficient = 1, .total_ion_activity_mol_per_m3 = 1, .electrical_conductivity_dS_per_m = 0 };
    const fluxes = try calculate(shared, zone, coefficients, 1, filled(EquilibriumConstants, 1), .{ .protonated_site_equilibrium_constant = 1, .hydroxyl_site_equilibrium_constant = 1, .h2po4_exchange_equilibrium_constant = 1, .hpo4_exchange_equilibrium_constant = 1, .water_activity_product_mol2_per_m6 = 1, .h2po4_dissociation_constant = 1, .maximum_exchange_mol_per_Mg_step = 0.1, .substrate_limit_fraction = 0.2 }, null, .{ .substrate_limit_fraction = 0.2, .maximum_pairing_mol_per_m3_step = 0.1 });
    const changes = try phosphate_network.assemble(fluxes);
    const dissolved = changes.dissolved_po4_mol_p_per_m3 + changes.dissolved_hpo4_mol_p_per_m3 + changes.dissolved_h2po4_mol_p_per_m3 + changes.dissolved_h3po4_mol_p_per_m3;
    const adsorbed = changes.adsorbed_hpo4_mol_p_per_Mg + changes.adsorbed_h2po4_mol_p_per_Mg;
    const paired = changes.iron_hpo4_pair_mol_per_m3 + changes.iron_h2po4_pair_mol_per_m3 + changes.calcium_po4_pair_mol_per_m3 + changes.calcium_hpo4_pair_mol_per_m3 + changes.calcium_h2po4_pair_mol_per_m3 + changes.magnesium_hpo4_pair_mol_per_m3;
    try std.testing.expectApproxEqAbs(@as(f64, 0), dissolved + adsorbed + paired, 1e-14);
}

test "five phosphate minerals are evaluated with distinct stoichiometry and ceilings" {
    var shared = filled(aqueous_network.State, 1);
    var zone = filled(phosphate_network.State, 1);
    shared.aluminum = 2;
    shared.iron = 2;
    shared.calcium = 2;
    zone.dissolved_h2po4_mol_p_per_m3 = 2;
    zone.dissolved_hpo4_mol_p_per_m3 = 2;
    const coefficients = activity_coefficients.Result{ .ionic_strength_mol_per_l = 0, .monovalent_activity_coefficient = 1, .divalent_activity_coefficient = 1, .trivalent_activity_coefficient = 1, .total_ion_activity_mol_per_m3 = 1, .electrical_conductivity_dS_per_m = 0 };
    const mineral_parameters = filled(MineralParameters, 1);
    const fluxes = try calculate(shared, zone, coefficients, 1, filled(EquilibriumConstants, 1), .{ .protonated_site_equilibrium_constant = 1, .hydroxyl_site_equilibrium_constant = 1, .h2po4_exchange_equilibrium_constant = 1, .hpo4_exchange_equilibrium_constant = 1, .water_activity_product_mol2_per_m6 = 1, .h2po4_dissociation_constant = 1, .maximum_exchange_mol_per_Mg_step = 0.1, .substrate_limit_fraction = 0.2 }, mineral_parameters, .{ .substrate_limit_fraction = 0.2, .maximum_pairing_mol_per_m3_step = 0.1 });
    try std.testing.expect(fluxes.minerals.aluminum_phosphate_mol_per_m3 > 0);
    try std.testing.expect(fluxes.minerals.iron_phosphate_mol_per_m3 > 0);
    try std.testing.expect(fluxes.minerals.dicalcium_phosphate_mol_per_m3 > 0);
    try std.testing.expect(fluxes.minerals.hydroxyapatite_mol_per_m3 > 0);
    try std.testing.expect(fluxes.minerals.monocalcium_phosphate_mol_per_m3 > 0);
}

test "phosphate aqueous rates match every SOLUTE source equation" {
    var shared: aqueous_network.State = undefined;
    var shared_value: f64 = 0.2;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        @field(shared, field.name) = shared_value;
        shared_value += 0.031;
    }
    var zone: phosphate_network.State = undefined;
    var zone_value: f64 = 0.3;
    inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
        @field(zone, field.name) = zone_value;
        zone_value += 0.029;
    }
    var constants: EquilibriumConstants = undefined;
    var constant_value: f64 = 0.1;
    inline for (@typeInfo(EquilibriumConstants).@"struct".fields) |field| {
        @field(constants, field.name) = constant_value;
        constant_value += 0.017;
    }
    const coefficients = activity_coefficients.Result{
        .ionic_strength_mol_per_l = 0.1,
        .monovalent_activity_coefficient = 0.8,
        .divalent_activity_coefficient = 0.6,
        .trivalent_activity_coefficient = 0.4,
        .total_ion_activity_mol_per_m3 = 1,
        .electrical_conductivity_dS_per_m = 0.2,
    };
    const kinetics = Kinetics{
        .substrate_limit_fraction = 0.7,
        .maximum_pairing_mol_per_m3_step = 10,
    };
    const fluxes = try calculate(
        shared,
        zone,
        coefficients,
        1.4,
        constants,
        .{
            .protonated_site_equilibrium_constant = 1,
            .hydroxyl_site_equilibrium_constant = 1,
            .h2po4_exchange_equilibrium_constant = 1,
            .hpo4_exchange_equilibrium_constant = 1,
            .water_activity_product_mol2_per_m6 = 1,
            .h2po4_dissociation_constant = 1,
            .maximum_exchange_mol_per_Mg_step = 0.1,
            .substrate_limit_fraction = 0.2,
        },
        null,
        kinetics,
    );
    const a = fluxes.aqueous;
    const g1 = coefficients.monovalent_activity_coefficient;
    const g2 = coefficients.divalent_activity_coefficient;
    const g3 = coefficients.trivalent_activity_coefficient;
    const limit = kinetics.substrate_limit_fraction;
    const maximum = kinetics.maximum_pairing_mol_per_m3_step;

    try expectSourcePairingRate(a.po4_hydrogen_association_mol_p_per_m3, zone.dissolved_po4_mol_p_per_m3, shared.hydrogen, zone.dissolved_hpo4_mol_p_per_m3, zone.dissolved_po4_mol_p_per_m3 * g3, shared.hydrogen * g1, zone.dissolved_hpo4_mol_p_per_m3 * g2, g3, constants.hpo4, limit, maximum);
    try expectSourcePairingRate(a.hpo4_hydrogen_association_mol_p_per_m3, zone.dissolved_hpo4_mol_p_per_m3, shared.hydrogen, zone.dissolved_h2po4_mol_p_per_m3, zone.dissolved_hpo4_mol_p_per_m3 * g2, shared.hydrogen * g1, zone.dissolved_h2po4_mol_p_per_m3 * g1, g2, constants.h2po4, limit, maximum);
    try expectSourcePairingRate(a.h2po4_hydrogen_association_mol_p_per_m3, zone.dissolved_h2po4_mol_p_per_m3, shared.hydrogen, zone.dissolved_h3po4_mol_p_per_m3, zone.dissolved_h2po4_mol_p_per_m3 * g1, shared.hydrogen * g1, zone.dissolved_h3po4_mol_p_per_m3, g1, constants.h3po4, limit, maximum);
    try expectSourcePairingRate(a.iron_hpo4_pairing_mol_p_per_m3, shared.iron, zone.dissolved_hpo4_mol_p_per_m3, zone.iron_hpo4_pair_mol_per_m3, shared.iron * g3, zone.dissolved_hpo4_mol_p_per_m3 * g2, zone.iron_hpo4_pair_mol_per_m3 * g2, g3, constants.iron_hpo4, limit, maximum);
    try expectSourcePairingRate(a.iron_h2po4_pairing_mol_p_per_m3, shared.iron, zone.dissolved_h2po4_mol_p_per_m3, zone.iron_h2po4_pair_mol_per_m3, shared.iron * g3, zone.dissolved_h2po4_mol_p_per_m3 * g1, zone.iron_h2po4_pair_mol_per_m3 * g1, g3, constants.iron_h2po4, limit, maximum);
    try std.testing.expectEqual(@as(f64, 0), a.calcium_po4_pairing_mol_p_per_m3);
    try expectSourcePairingRate(a.calcium_hpo4_pairing_mol_p_per_m3, shared.calcium, zone.dissolved_hpo4_mol_p_per_m3, zone.calcium_hpo4_pair_mol_per_m3, shared.calcium * g2, zone.dissolved_hpo4_mol_p_per_m3 * g2, zone.calcium_hpo4_pair_mol_per_m3, g2, constants.calcium_hpo4, limit, maximum);
    try expectSourcePairingRate(a.calcium_h2po4_pairing_mol_p_per_m3, shared.calcium, zone.dissolved_h2po4_mol_p_per_m3, zone.calcium_h2po4_pair_mol_per_m3, shared.calcium * g2, zone.dissolved_h2po4_mol_p_per_m3 * g1, zone.calcium_h2po4_pair_mol_per_m3 * g1, g2, constants.calcium_h2po4, limit, maximum);
    try expectSourcePairingRate(a.magnesium_hpo4_pairing_mol_p_per_m3, shared.magnesium, zone.dissolved_hpo4_mol_p_per_m3, zone.magnesium_hpo4_pair_mol_per_m3, shared.magnesium * g2, zone.dissolved_hpo4_mol_p_per_m3 * g2, zone.magnesium_hpo4_pair_mol_per_m3, g2, constants.magnesium_hpo4, limit, maximum);
}

fn expectSourcePairingRate(
    actual: f64,
    free_first: f64,
    free_second: f64,
    paired: f64,
    first_activity: f64,
    second_activity: f64,
    paired_activity: f64,
    first_coefficient: f64,
    dissociation_constant: f64,
    substrate_limit_fraction: f64,
    maximum: f64,
) !void {
    const dissociation_limit = substrate_limit_fraction * paired;
    const association_limit =
        substrate_limit_fraction * @min(free_first, free_second);
    const equilibrium_first_activity =
        dissociation_constant * paired_activity / second_activity;
    const expected = @max(
        -maximum,
        -dissociation_limit,
        @min(
            maximum,
            association_limit,
            (first_activity - equilibrium_first_activity) / first_coefficient,
        ),
    );
    try std.testing.expectApproxEqAbs(expected, actual, 1e-15);
}

test "source-order phosphate mineral limiters remain explicit" {
    var shared = filled(aqueous_network.State, 1);
    var zone = filled(phosphate_network.State, 1);
    shared.aluminum = 1;
    shared.iron = 1;
    shared.calcium = 1;
    shared.hydrogen = 1;
    shared.hydroxide = 1;
    zone.dissolved_hpo4_mol_p_per_m3 = 0.8;
    zone.dissolved_h2po4_mol_p_per_m3 = 0.2;
    const coefficients = activity_coefficients.Result{
        .ionic_strength_mol_per_l = 0,
        .monovalent_activity_coefficient = 0.8,
        .divalent_activity_coefficient = 0.6,
        .trivalent_activity_coefficient = 0.4,
        .total_ion_activity_mol_per_m3 = 1,
        .electrical_conductivity_dS_per_m = 0,
    };
    var parameters = filled(MineralParameters, 1e-30);
    parameters.water_activity_product_mol2_per_m6 = 1;
    parameters.maximum_phosphate_precipitation_mol_per_m3_step = 10;
    parameters.maximum_apatite_precipitation_mol_per_m3_step = 10;
    parameters.maximum_mineral_dissolution_mol_per_m3_step = 10;
    const constants = filled(EquilibriumConstants, 1);
    const source = try calculateMineralsSourceOrder(.{
        .shared = shared,
        .zone = zone,
        .zone_identity = .band,
        .coefficients = coefficients,
        .constants = constants,
        .parameters = parameters,
        .substrate_limit_fraction = 0.5,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), source.aluminum_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), source.iron_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), source.dicalcium_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0 / 30.0), source.hydroxyapatite_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), source.monocalcium_phosphate_mol_per_m3, 1e-15);

    const production = try calculateMinerals(
        shared,
        zone,
        coefficients,
        constants,
        parameters,
        0.5,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), production.dicalcium_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), production.monocalcium_phosphate_mol_per_m3, 1e-15);
}

test "source-order phosphate minerals retain distinct dissolution ceilings" {
    const shared = filled(aqueous_network.State, 1);
    var zone = filled(phosphate_network.State, 1);
    zone.dissolved_hpo4_mol_p_per_m3 = 0;
    zone.dissolved_h2po4_mol_p_per_m3 = 0;
    const coefficients = activity_coefficients.Result{
        .ionic_strength_mol_per_l = 0,
        .monovalent_activity_coefficient = 1,
        .divalent_activity_coefficient = 1,
        .trivalent_activity_coefficient = 1,
        .total_ion_activity_mol_per_m3 = 1,
        .electrical_conductivity_dS_per_m = 0,
    };
    var parameters = filled(MineralParameters, 1);
    parameters.maximum_phosphate_precipitation_mol_per_m3_step = 0.2;
    parameters.maximum_apatite_precipitation_mol_per_m3_step = 0.3;
    parameters.maximum_mineral_dissolution_mol_per_m3_step = 0.04;
    const constants = filled(EquilibriumConstants, 1);
    const source = try calculateMineralsSourceOrder(.{
        .shared = shared,
        .zone = zone,
        .zone_identity = .non_band,
        .coefficients = coefficients,
        .constants = constants,
        .parameters = parameters,
        .substrate_limit_fraction = 0.5,
    });
    try std.testing.expectApproxEqAbs(@as(f64, -0.2), source.aluminum_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.2), source.iron_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.2), source.dicalcium_phosphate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.3), source.hydroxyapatite_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.04), source.monocalcium_phosphate_mol_per_m3, 1e-15);

    const production = try calculateMinerals(
        shared,
        zone,
        coefficients,
        constants,
        parameters,
        0.5,
    );
    inline for (@typeInfo(phosphate_network.MineralFluxes).@"struct".fields) |field|
        try std.testing.expectApproxEqAbs(
            @as(f64, -0.04),
            @field(production, field.name),
            1e-15,
        );
}
