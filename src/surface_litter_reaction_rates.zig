const std = @import("std");
const chemistry = @import("surface_litter_chemistry.zig");
const ledger = @import("surface_litter_reaction_ledger.zig");
const activity_coefficients = @import("solute_activity_coefficients.zig");
const ion_pairing = @import("solute_ion_pairing.zig");
const cation_exchange = @import("solute_cation_exchange.zig");
const phosphate_precipitation = @import("solute_phosphate_precipitation.zig");
const mineral_precipitation = @import("solute_mineral_precipitation.zig");
const water_equilibrium = @import("solute_water_equilibrium.zig");
const phosphate_exchange = @import("solute_phosphate_exchange.zig");

pub const DissociationConstants = struct {
    ammonium: f64,
    carbon_dioxide: f64,
    bicarbonate: f64,
    h2po4: f64,
    hpo4: f64,
    carboxyl: f64,
};

pub const MineralProducts = struct {
    aluminum_phosphate: f64,
    iron_phosphate: f64,
    dicalcium_phosphate: f64,
    hydroxyapatite: f64,
    monocalcium_phosphate: f64,
    gibbsite: f64,
    iron_hydroxide: f64,
    calcite: f64,
    gypsum: f64,
    fixed_ph_aluminum_h2po4: f64,
    fixed_ph_iron_h2po4: f64,
    fixed_ph_hydroxyapatite_h2po4: f64,
};

pub const Kinetics = struct {
    ammonium_substrate_limit_fraction: f64,
    general_substrate_limit_fraction: f64,
    maximum_ammonium_association_mol_per_m3_step: f64,
    maximum_association_mol_per_m3_step: f64,
    maximum_phosphate_precipitation_mol_per_m3_step: f64,
    maximum_apatite_precipitation_mol_per_m3_step: f64,
    maximum_monocalcium_dissolution_mol_per_m3_step: f64,
    maximum_cation_adsorption_mol_charge_per_m3_step: f64,
    calcite_hydroxide_inhibition_constant_mol_per_m3: f64,
};

pub const Parameters = struct {
    activity: activity_coefficients.Result,
    dissociation: DissociationConstants,
    minerals: MineralProducts,
    kinetics: Kinetics,
    cation_exchange_capacity_mol_charge_per_megagram: f64,
    cation_selectivity: cation_exchange.Selectivity,
    water_activity_product_mol2_per_m6: f64,
    negligible_water_ion_concentration_mol_per_m3: f64,
    external_hydrogen_mol_per_m3: f64 = 0,
    phosphate_surface: phosphate_exchange.Parameters = .{
        .protonated_site_equilibrium_constant = 0,
        .hydroxyl_site_equilibrium_constant = 0,
        .h2po4_exchange_equilibrium_constant = 0,
        .hpo4_exchange_equilibrium_constant = 0,
        .water_activity_product_mol2_per_m6 = 0,
        .h2po4_dissociation_constant = 0,
        .maximum_exchange_mol_per_megagram_step = 0,
        .substrate_limit_fraction = 0,
    },
};

pub const Context = struct {
    parameters: Parameters,
    litter_mass_per_water_volume_megagrams_per_m3: f64,
    dynamic_salts: bool,

    pub fn evaluator(self: *const Context) chemistry.Evaluator {
        return .{
            .context = self,
            .evaluate = evaluateOpaque,
            .equilibrate_cation_exchange = equilibrateCationExchangeOpaque,
            .phosphate_mineral_equilibrium_residuals = phosphateMineralEquilibriumResidualsOpaque,
        };
    }
};

pub fn evaluateOpaque(raw: *const anyopaque, cell: chemistry.Cell) !ledger.ReactionExtents {
    const context: *const Context = @ptrCast(@alignCast(raw));
    return calculate(cell, context.*);
}

fn phosphateMineralEquilibriumResidualsOpaque(
    raw: *const anyopaque,
    cell: chemistry.Cell,
) !ledger.PhosphateMineralExtents {
    const context: *const Context = @ptrCast(@alignCast(raw));
    try validate(cell, context.*);
    const p = context.parameters;
    const g1 = p.activity.monovalent_activity_coefficient;
    const g2 = p.activity.divalent_activity_coefficient;
    const water = try water_equilibrium.solve(.{
        .hydrogen_concentration_mol_per_m3 = cell.hydrogen_mol_per_m3,
        .hydroxide_concentration_mol_per_m3 = cell.hydroxide_mol_per_m3,
        .monovalent_activity_coefficient = g1,
        .water_activity_product_mol2_per_m6 = p.water_activity_product_mol2_per_m6,
        .negligible_concentration_mol_per_m3 = p.negligible_water_ion_concentration_mol_per_m3,
    });
    const hydrogen_activity =
        water.hydrogen_concentration_mol_per_m3 * g1;
    const hydroxide_activity =
        water.hydroxide_concentration_mol_per_m3 * g1;
    const aluminum_activity =
        cell.aluminum_mol_per_m3 *
        p.activity.trivalent_activity_coefficient;
    const iron_activity =
        cell.iron_mol_per_m3 *
        p.activity.trivalent_activity_coefficient;
    const calcium_activity = cell.calcium_mol_per_m3 * g2;
    if (aluminum_activity <= 0 or iron_activity <= 0 or
        calcium_activity <= 0 or hydrogen_activity <= 0 or
        hydroxide_activity <= 0)
        return error.InvalidLitterMineralActivity;
    const aluminum_target = if (context.dynamic_salts)
        try phosphate_precipitation.aluminumOrIronPhosphateEquilibriumH2po4(
            p.minerals.aluminum_phosphate,
            hydrogen_activity,
            p.dissociation.h2po4,
            p.dissociation.hpo4,
            aluminum_activity,
        )
    else
        p.minerals.fixed_ph_aluminum_h2po4 *
            hydrogen_activity * hydrogen_activity /
            (p.minerals.gibbsite /
                std.math.pow(f64, hydroxide_activity, 3));
    const iron_target = if (context.dynamic_salts)
        try phosphate_precipitation.aluminumOrIronPhosphateEquilibriumH2po4(
            p.minerals.iron_phosphate,
            hydrogen_activity,
            p.dissociation.h2po4,
            p.dissociation.hpo4,
            iron_activity,
        )
    else
        p.minerals.fixed_ph_iron_h2po4 *
            hydrogen_activity * hydrogen_activity /
            (p.minerals.iron_hydroxide /
                std.math.pow(f64, hydroxide_activity, 3));
    const dicalcium_target =
        try phosphate_precipitation.dicalciumPhosphateEquilibriumHpo4(
            p.minerals.dicalcium_phosphate,
            calcium_activity,
        );
    const apatite_target = if (context.dynamic_salts)
        try phosphate_precipitation.hydroxyapatiteEquilibriumH2po4(
            p.minerals.hydroxyapatite,
            hydrogen_activity,
            hydroxide_activity,
            calcium_activity,
            p.dissociation.h2po4,
            p.dissociation.hpo4,
        )
    else
        std.math.cbrt(
            p.minerals.fixed_ph_hydroxyapatite_h2po4 *
                std.math.pow(f64, hydrogen_activity, 7) /
                std.math.pow(f64, calcium_activity, 5),
        );
    const monocalcium_target =
        try phosphate_precipitation.monocalciumPhosphateEquilibriumH2po4(
            p.minerals.monocalcium_phosphate,
            calcium_activity,
        );
    return .{
        .aluminum_phosphate_mol_per_m3 = complementarityResidual(
            cell.phosphate_minerals.aluminum_phosphate_mol_per_m3,
            cell.h2po4_mol_p_per_m3 * g1 - aluminum_target,
        ),
        .iron_phosphate_mol_per_m3 = complementarityResidual(
            cell.phosphate_minerals.iron_phosphate_mol_per_m3,
            cell.h2po4_mol_p_per_m3 * g1 - iron_target,
        ),
        .dicalcium_phosphate_mol_per_m3 = complementarityResidual(
            cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3,
            cell.hpo4_mol_p_per_m3 * g2 - dicalcium_target,
        ),
        .hydroxyapatite_mol_per_m3 = complementarityResidual(
            cell.phosphate_minerals.hydroxyapatite_mol_per_m3,
            cell.h2po4_mol_p_per_m3 * g1 - apatite_target,
        ),
        .monocalcium_phosphate_mol_per_m3 = complementarityResidual(
            cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3,
            cell.h2po4_mol_p_per_m3 * g1 - monocalcium_target,
        ),
    };
}

fn complementarityResidual(solid_mol_per_m3: f64, saturation_residual: f64) f64 {
    return if (solid_mol_per_m3 <= 0)
        @max(0, saturation_residual)
    else
        saturation_residual;
}

fn equilibrateCationExchangeOpaque(raw: *const anyopaque, cell: chemistry.Cell) !chemistry.Cell {
    const context: *const Context = @ptrCast(@alignCast(raw));
    const density = context.litter_mass_per_water_volume_megagrams_per_m3;
    const p = context.parameters;
    if (p.cation_exchange_capacity_mol_charge_per_megagram == 0) return cell;
    const totals = ledger.ExchangeAdsorption{
        .ammonium_mol_per_megagram = (cell.ammonium_mol_per_m3 +
            cell.ammonia_mol_per_m3) / density +
            cell.exchange.ammonium_mol_per_megagram,
        .hydrogen_mol_per_megagram = 0,
        .aluminum_mol_per_megagram = cell.aluminum_mol_per_m3 / density +
            cell.exchange.aluminum_mol_per_megagram,
        .iron_mol_per_megagram = cell.iron_mol_per_m3 / density +
            cell.exchange.iron_mol_per_megagram,
        .calcium_mol_per_megagram = cell.calcium_mol_per_m3 / density +
            cell.exchange.calcium_mol_per_megagram,
        .magnesium_mol_per_megagram = cell.magnesium_mol_per_m3 / density +
            cell.exchange.magnesium_mol_per_megagram,
        .sodium_mol_per_megagram = cell.sodium_mol_per_m3 / density +
            cell.exchange.sodium_mol_per_megagram,
        .potassium_mol_per_megagram = cell.potassium_mol_per_m3 / density +
            cell.exchange.potassium_mol_per_megagram,
    };
    var result = cell;
    var iteration: u8 = 0;
    while (iteration < 96) : (iteration += 1) {
        const concentrations = cation_exchange.Cations{
            .ammonium_non_band = result.ammonium_mol_per_m3,
            .ammonium_band = 0,
            .hydrogen = result.hydrogen_mol_per_m3,
            .aluminum = result.aluminum_mol_per_m3,
            .iron = result.iron_mol_per_m3,
            .calcium = result.calcium_mol_per_m3,
            .magnesium = result.magnesium_mol_per_m3,
            .sodium = result.sodium_mol_per_m3,
            .potassium = result.potassium_mol_per_m3,
        };
        var activities = concentrations;
        activities.ammonium_non_band *= p.activity.monovalent_activity_coefficient;
        activities.hydrogen *= p.activity.monovalent_activity_coefficient;
        activities.aluminum *= p.activity.trivalent_activity_coefficient;
        activities.iron *= p.activity.trivalent_activity_coefficient;
        activities.calcium *= p.activity.divalent_activity_coefficient;
        activities.magnesium *= p.activity.divalent_activity_coefficient;
        activities.sodium *= p.activity.monovalent_activity_coefficient;
        activities.potassium *= p.activity.monovalent_activity_coefficient;
        const current_exchange = cation_exchange.Cations{
            .ammonium_non_band = result.exchange.ammonium_mol_per_megagram,
            .ammonium_band = 0,
            .hydrogen = result.exchange.hydrogen_mol_per_megagram,
            .aluminum = result.exchange.aluminum_mol_per_megagram,
            .iron = result.exchange.iron_mol_per_megagram,
            .calcium = result.exchange.calcium_mol_per_megagram,
            .magnesium = result.exchange.magnesium_mol_per_megagram,
            .sodium = result.exchange.sodium_mol_per_megagram,
            .potassium = result.exchange.potassium_mol_per_megagram,
        };
        const target = try cation_exchange.equilibriumIonConcentration(.{
            .cation_exchange_capacity_mol_charge_per_megagram = p.cation_exchange_capacity_mol_charge_per_megagram,
            .aqueous_concentration_mol_per_m3 = concentrations,
            .aqueous_activity_mol_per_m3 = activities,
            .exchange_concentration_mol_per_megagram = current_exchange,
            .ammonium_non_band_fraction = 1,
            .ammonium_band_fraction = 0,
            .soil_mass_per_water_volume_megagrams_per_m3 = density,
        }, p.cation_selectivity);
        var fraction: f64 = 0.5;
        fraction = admissibleExchangeFraction(fraction, result.exchange.ammonium_mol_per_megagram, target.ammonium_non_band, totals.ammonium_mol_per_megagram);
        fraction = admissibleExchangeFraction(fraction, result.exchange.aluminum_mol_per_megagram, target.aluminum, totals.aluminum_mol_per_megagram);
        fraction = admissibleExchangeFraction(fraction, result.exchange.iron_mol_per_megagram, target.iron, totals.iron_mol_per_megagram);
        fraction = admissibleExchangeFraction(fraction, result.exchange.calcium_mol_per_megagram, target.calcium, totals.calcium_mol_per_megagram);
        fraction = admissibleExchangeFraction(fraction, result.exchange.magnesium_mol_per_megagram, target.magnesium, totals.magnesium_mol_per_megagram);
        fraction = admissibleExchangeFraction(fraction, result.exchange.sodium_mol_per_megagram, target.sodium, totals.sodium_mol_per_megagram);
        fraction = admissibleExchangeFraction(fraction, result.exchange.potassium_mol_per_megagram, target.potassium, totals.potassium_mol_per_megagram);
        if (fraction <= std.math.floatEps(f64)) break;
        const before = result.exchange;
        result.exchange = .{
            .ammonium_mol_per_megagram = before.ammonium_mol_per_megagram + fraction * (target.ammonium_non_band - before.ammonium_mol_per_megagram),
            .hydrogen_mol_per_megagram = before.hydrogen_mol_per_megagram + fraction * (target.hydrogen - before.hydrogen_mol_per_megagram),
            .aluminum_mol_per_megagram = before.aluminum_mol_per_megagram + fraction * (target.aluminum - before.aluminum_mol_per_megagram),
            .iron_mol_per_megagram = before.iron_mol_per_megagram + fraction * (target.iron - before.iron_mol_per_megagram),
            .calcium_mol_per_megagram = before.calcium_mol_per_megagram + fraction * (target.calcium - before.calcium_mol_per_megagram),
            .magnesium_mol_per_megagram = before.magnesium_mol_per_megagram + fraction * (target.magnesium - before.magnesium_mol_per_megagram),
            .sodium_mol_per_megagram = before.sodium_mol_per_megagram + fraction * (target.sodium - before.sodium_mol_per_megagram),
            .potassium_mol_per_megagram = before.potassium_mol_per_megagram + fraction * (target.potassium - before.potassium_mol_per_megagram),
        };
        const aqueous_n_per_megagram =
            totals.ammonium_mol_per_megagram - result.exchange.ammonium_mol_per_megagram;
        const ammonia_to_ammonium =
            p.dissociation.ammonium / result.hydrogen_mol_per_m3;
        result.ammonium_mol_per_m3 =
            density * aqueous_n_per_megagram / (1 + ammonia_to_ammonium);
        result.ammonia_mol_per_m3 =
            density * aqueous_n_per_megagram - result.ammonium_mol_per_m3;
        result.aluminum_mol_per_m3 = density *
            (totals.aluminum_mol_per_megagram - result.exchange.aluminum_mol_per_megagram);
        result.iron_mol_per_m3 = density *
            (totals.iron_mol_per_megagram - result.exchange.iron_mol_per_megagram);
        result.calcium_mol_per_m3 = density *
            (totals.calcium_mol_per_megagram - result.exchange.calcium_mol_per_megagram);
        result.magnesium_mol_per_m3 = density *
            (totals.magnesium_mol_per_megagram - result.exchange.magnesium_mol_per_megagram);
        result.sodium_mol_per_m3 = density *
            (totals.sodium_mol_per_megagram - result.exchange.sodium_mol_per_megagram);
        result.potassium_mol_per_m3 = density *
            (totals.potassium_mol_per_megagram - result.exchange.potassium_mol_per_megagram);
        try normalizeExchangeReconstruction(&result, totals, density);
        if (maximumExchangeDifference(before, result.exchange) <=
            8 * std.math.floatEps(f64) *
                @max(1.0, maximumExchangeMagnitude(result.exchange)))
            break;
    }
    return result;
}

fn normalizeExchangeReconstruction(
    cell: *chemistry.Cell,
    totals: ledger.ExchangeAdsorption,
    density: f64,
) !void {
    inline for ([_][]const u8{
        "ammonium_mol_per_m3",
        "ammonia_mol_per_m3",
        "aluminum_mol_per_m3",
        "iron_mol_per_m3",
        "calcium_mol_per_m3",
        "magnesium_mol_per_m3",
        "sodium_mol_per_m3",
        "potassium_mol_per_m3",
    }) |field_name| {
        const value = &@field(cell.*, field_name);
        const total_per_megagram = if (comptime std.mem.eql(u8, field_name, "ammonium_mol_per_m3") or
            std.mem.eql(u8, field_name, "ammonia_mol_per_m3"))
            totals.ammonium_mol_per_megagram
        else if (comptime std.mem.eql(u8, field_name, "aluminum_mol_per_m3"))
            totals.aluminum_mol_per_megagram
        else if (comptime std.mem.eql(u8, field_name, "iron_mol_per_m3"))
            totals.iron_mol_per_megagram
        else if (comptime std.mem.eql(u8, field_name, "calcium_mol_per_m3"))
            totals.calcium_mol_per_megagram
        else if (comptime std.mem.eql(u8, field_name, "magnesium_mol_per_m3"))
            totals.magnesium_mol_per_megagram
        else if (comptime std.mem.eql(u8, field_name, "sodium_mol_per_m3"))
            totals.sodium_mol_per_megagram
        else
            totals.potassium_mol_per_megagram;
        const roundoff_limit = 1024 * std.math.floatEps(f64) *
            @max(1.0, @abs(density * total_per_megagram));
        if (!std.math.isFinite(value.*))
            return error.NonFiniteLitterExchangeReconstruction;
        if (value.* < 0) {
            if (value.* < -roundoff_limit)
                return error.NegativeLitterExchangeReconstruction;
            value.* = 0;
        }
    }
}

fn admissibleExchangeFraction(requested: f64, current: f64, target: f64, total: f64) f64 {
    if (target <= current) return requested;
    return @min(requested, @max(0, total - current) / (target - current));
}

fn maximumExchangeDifference(a: ledger.ExchangeAdsorption, b: ledger.ExchangeAdsorption) f64 {
    var maximum: f64 = 0;
    inline for (@typeInfo(ledger.ExchangeAdsorption).@"struct".fields) |field|
        maximum = @max(maximum, @abs(@field(a, field.name) - @field(b, field.name)));
    return maximum;
}

fn maximumExchangeMagnitude(value: ledger.ExchangeAdsorption) f64 {
    var maximum: f64 = 0;
    inline for (@typeInfo(ledger.ExchangeAdsorption).@"struct".fields) |field|
        maximum = @max(maximum, @abs(@field(value, field.name)));
    return maximum;
}

/// Runs the litter formulation carried by `context`; density and salt mode are
/// sourced once so the evaluator and conservative ledger cannot disagree.
pub fn solveCell(state: *chemistry.State, cell_index: usize, context: *const Context, options: chemistry.Options) !chemistry.Result {
    return chemistry.solveCell(state, cell_index, .{
        .litter_mass_per_water_volume_megagrams_per_m3 = context.litter_mass_per_water_volume_megagrams_per_m3,
        .dynamic_salts = context.dynamic_salts,
    }, context.evaluator(), options);
}

/// Litter-only equilibrium rates from SOLUTE.F. The litter has one
/// aqueous/exchange domain, phosphate surface sites, and five coprecipitates.
pub fn calculate(cell: chemistry.Cell, context: Context) !ledger.ReactionExtents {
    try validate(cell, context);
    const p = context.parameters;
    const g1 = p.activity.monovalent_activity_coefficient;
    const g2 = p.activity.divalent_activity_coefficient;
    const limit = p.kinetics.general_substrate_limit_fraction;
    const maximum = p.kinetics.maximum_association_mol_per_m3_step;

    const water = try water_equilibrium.solve(.{
        .hydrogen_concentration_mol_per_m3 = cell.hydrogen_mol_per_m3,
        .hydroxide_concentration_mol_per_m3 = cell.hydroxide_mol_per_m3,
        .monovalent_activity_coefficient = g1,
        .water_activity_product_mol2_per_m6 = p.water_activity_product_mol2_per_m6,
        .negligible_concentration_mol_per_m3 = p.negligible_water_ion_concentration_mol_per_m3,
    });
    const hydrogen = water.hydrogen_concentration_mol_per_m3;
    const hydroxide = water.hydroxide_concentration_mol_per_m3;
    const hydrogen_activity = hydrogen * g1;
    const hydroxide_activity = hydroxide * g1;
    const recombination = cell.hydrogen_mol_per_m3 - hydrogen;
    if (@abs(recombination - (cell.hydroxide_mol_per_m3 - hydroxide)) > 1e-10 * @max(1.0, @abs(recombination))) return error.LitterWaterEquilibriumImbalance;

    var result = zeroExtents();
    result.water_ion_recombination_mol_per_m3 = recombination;
    result.external_hydrogen_mol_per_m3 = p.external_hydrogen_mol_per_m3;
    result.ammonium_association_mol_per_m3 = try association(cell.ammonia_mol_per_m3, hydrogen, cell.ammonium_mol_per_m3, cell.ammonia_mol_per_m3, hydrogen_activity, cell.ammonium_mol_per_m3 * g1, 1, p.dissociation.ammonium, p.kinetics.ammonium_substrate_limit_fraction, p.kinetics.maximum_ammonium_association_mol_per_m3_step);
    result.h2po4_association_mol_p_per_m3 = try association(cell.hpo4_mol_p_per_m3, hydrogen, cell.h2po4_mol_p_per_m3, cell.hpo4_mol_p_per_m3 * g2, hydrogen_activity, cell.h2po4_mol_p_per_m3 * g1, g2, p.dissociation.h2po4, limit, maximum);

    if (context.dynamic_salts) {
        result.carbonate_hydrogen_association_mol_per_m3 = try association(cell.carbonate_mol_per_m3, hydrogen, cell.bicarbonate_mol_per_m3, cell.carbonate_mol_per_m3 * g2, hydrogen_activity, cell.bicarbonate_mol_per_m3 * g1, g2, p.dissociation.bicarbonate, limit, maximum);
        result.bicarbonate_hydrogen_association_mol_per_m3 = try association(cell.bicarbonate_mol_per_m3, hydrogen, cell.carbon_dioxide_mol_per_m3, cell.bicarbonate_mol_per_m3 * g1, hydrogen_activity, cell.carbon_dioxide_mol_per_m3, g1, p.dissociation.carbon_dioxide, limit, maximum);
    }

    result.exchange = try exchangeRates(cell, hydrogen, context);
    const surface_site_total = cell.phosphate_surface.deprotonated_site_mol_per_megagram + cell.phosphate_surface.hydroxyl_site_mol_per_megagram + cell.phosphate_surface.protonated_site_mol_per_megagram + cell.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram + cell.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram;
    if (surface_site_total > 0) result.phosphate_surface = try phosphate_exchange.calculate(.{
        .hydrogen_concentration_mol_per_m3 = hydrogen,
        .hydrogen_activity_mol_per_m3 = hydrogen_activity,
        .hydroxide_activity_mol_per_m3 = hydroxide_activity,
        .h2po4_concentration_mol_p_per_m3 = cell.h2po4_mol_p_per_m3,
        .h2po4_activity_mol_p_per_m3 = cell.h2po4_mol_p_per_m3 * g1,
        .hpo4_concentration_mol_p_per_m3 = cell.hpo4_mol_p_per_m3,
        .hpo4_activity_mol_p_per_m3 = cell.hpo4_mol_p_per_m3 * g2,
        .deprotonated_site_mol_per_megagram = cell.phosphate_surface.deprotonated_site_mol_per_megagram,
        .hydroxyl_site_mol_per_megagram = cell.phosphate_surface.hydroxyl_site_mol_per_megagram,
        .protonated_site_mol_per_megagram = cell.phosphate_surface.protonated_site_mol_per_megagram,
        .adsorbed_h2po4_mol_p_per_megagram = cell.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram,
        .adsorbed_hpo4_mol_p_per_megagram = cell.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram,
        .monovalent_activity_coefficient = g1,
        .divalent_activity_coefficient = g2,
    }, p.phosphate_surface);
    result.carboxyl_hydrogen_adsorption_mol_per_megagram = if (context.dynamic_salts) try carboxylRate(cell, hydrogen_activity, context) else 0;
    result.phosphate_minerals = try phosphateMinerals(cell, hydrogen_activity, hydroxide_activity, context);

    if (context.dynamic_salts) {
        const without_salt = try ledger.assemble(result, context.litter_mass_per_water_volume_megagrams_per_m3, true);
        result.salt_minerals = try saltMinerals(cell, without_salt, hydroxide_activity, context);
    }
    return result;
}

fn exchangeRates(cell: chemistry.Cell, hydrogen: f64, context: Context) !ledger.ExchangeAdsorption {
    const p = context.parameters;
    const concentrations = cation_exchange.Cations{
        .ammonium_non_band = cell.ammonium_mol_per_m3,
        .ammonium_band = 0,
        .hydrogen = hydrogen,
        .aluminum = cell.aluminum_mol_per_m3,
        .iron = cell.iron_mol_per_m3,
        .calcium = cell.calcium_mol_per_m3,
        .magnesium = cell.magnesium_mol_per_m3,
        .sodium = cell.sodium_mol_per_m3,
        .potassium = cell.potassium_mol_per_m3,
    };
    var activities = concentrations;
    activities.ammonium_non_band *= p.activity.monovalent_activity_coefficient;
    activities.hydrogen *= p.activity.monovalent_activity_coefficient;
    activities.aluminum *= p.activity.trivalent_activity_coefficient;
    activities.iron *= p.activity.trivalent_activity_coefficient;
    activities.calcium *= p.activity.divalent_activity_coefficient;
    activities.magnesium *= p.activity.divalent_activity_coefficient;
    activities.sodium *= p.activity.monovalent_activity_coefficient;
    activities.potassium *= p.activity.monovalent_activity_coefficient;
    const exchange_state = cation_exchange.Cations{
        .ammonium_non_band = cell.exchange.ammonium_mol_per_megagram,
        .ammonium_band = 0,
        .hydrogen = cell.exchange.hydrogen_mol_per_megagram,
        .aluminum = cell.exchange.aluminum_mol_per_megagram,
        .iron = cell.exchange.iron_mol_per_megagram,
        .calcium = cell.exchange.calcium_mol_per_megagram,
        .magnesium = cell.exchange.magnesium_mol_per_megagram,
        .sodium = cell.exchange.sodium_mol_per_megagram,
        .potassium = cell.exchange.potassium_mol_per_megagram,
    };
    const rates = try cation_exchange.calculateSourceOrder(.{
        .cation_exchange_capacity_mol_charge_per_megagram = p.cation_exchange_capacity_mol_charge_per_megagram,
        .aqueous_concentration_mol_per_m3 = concentrations,
        .aqueous_activity_mol_per_m3 = activities,
        .exchange_concentration_mol_per_megagram = exchange_state,
        .ammonium_non_band_fraction = 1,
        .ammonium_band_fraction = 0,
        .soil_mass_per_water_volume_megagrams_per_m3 = context.litter_mass_per_water_volume_megagrams_per_m3,
    }, .{ .selectivity = p.cation_selectivity, .substrate_limit_fraction = p.kinetics.general_substrate_limit_fraction, .maximum_adsorption_mol_charge_per_m3_step = p.kinetics.maximum_cation_adsorption_mol_charge_per_m3_step }, .{
        .minimum_activity_mol_per_m3 = p.negligible_water_ion_concentration_mol_per_m3,
    });
    return .{ .ammonium_mol_per_megagram = rates.ammonium_non_band, .hydrogen_mol_per_megagram = rates.hydrogen, .aluminum_mol_per_megagram = rates.aluminum, .iron_mol_per_megagram = rates.iron, .calcium_mol_per_megagram = rates.calcium, .magnesium_mol_per_megagram = rates.magnesium, .sodium_mol_per_megagram = rates.sodium, .potassium_mol_per_megagram = rates.potassium };
}

fn carboxylRate(cell: chemistry.Cell, hydrogen_activity: f64, context: Context) !f64 {
    if (context.parameters.cation_exchange_capacity_mol_charge_per_megagram == 0) return 0;
    if (hydrogen_activity <= 0) return error.InvalidLitterHydrogenActivity;
    const occupied = cell.carboxyl_hydrogen_mol_per_megagram;
    const open = @max(0, context.parameters.cation_exchange_capacity_mol_charge_per_megagram - occupied);
    const equilibrium_open = @min(context.parameters.cation_exchange_capacity_mol_charge_per_megagram, context.parameters.dissociation.carboxyl * occupied / hydrogen_activity);
    const density = context.litter_mass_per_water_volume_megagrams_per_m3;
    const maximum_per_megagram = context.parameters.kinetics.maximum_cation_adsorption_mol_charge_per_m3_step / density;
    const substrate_limit = context.parameters.kinetics.general_substrate_limit_fraction / density * occupied;
    return @max(-maximum_per_megagram, -substrate_limit, @min(maximum_per_megagram, substrate_limit, open - equilibrium_open));
}

fn phosphateMinerals(cell: chemistry.Cell, hydrogen_activity: f64, hydroxide_activity: f64, context: Context) !ledger.PhosphateMineralExtents {
    const p = context.parameters;
    const g1 = p.activity.monovalent_activity_coefficient;
    const g2 = p.activity.divalent_activity_coefficient;
    const aluminum_activity = cell.aluminum_mol_per_m3 * p.activity.trivalent_activity_coefficient;
    const iron_activity = cell.iron_mol_per_m3 * p.activity.trivalent_activity_coefficient;
    const calcium_activity = cell.calcium_mol_per_m3 * g2;
    if (aluminum_activity <= 0 or iron_activity <= 0 or calcium_activity <= 0 or hydrogen_activity <= 0 or hydroxide_activity <= 0) return error.InvalidLitterMineralActivity;
    const dynamic = context.dynamic_salts;
    const aluminum_target = if (dynamic)
        try phosphate_precipitation.aluminumOrIronPhosphateEquilibriumH2po4(p.minerals.aluminum_phosphate, hydrogen_activity, p.dissociation.h2po4, p.dissociation.hpo4, aluminum_activity)
    else
        p.minerals.fixed_ph_aluminum_h2po4 * hydrogen_activity * hydrogen_activity / (p.minerals.gibbsite / std.math.pow(f64, hydroxide_activity, 3));
    const iron_target = if (dynamic)
        try phosphate_precipitation.aluminumOrIronPhosphateEquilibriumH2po4(p.minerals.iron_phosphate, hydrogen_activity, p.dissociation.h2po4, p.dissociation.hpo4, iron_activity)
    else
        p.minerals.fixed_ph_iron_h2po4 * hydrogen_activity * hydrogen_activity / (p.minerals.iron_hydroxide / std.math.pow(f64, hydroxide_activity, 3));
    const dicalcium_target = try phosphate_precipitation.dicalciumPhosphateEquilibriumHpo4(p.minerals.dicalcium_phosphate, calcium_activity);
    const apatite_target = if (dynamic)
        try phosphate_precipitation.hydroxyapatiteEquilibriumH2po4(p.minerals.hydroxyapatite, hydrogen_activity, hydroxide_activity, calcium_activity, p.dissociation.h2po4, p.dissociation.hpo4)
    else
        std.math.cbrt(p.minerals.fixed_ph_hydroxyapatite_h2po4 * std.math.pow(f64, hydrogen_activity, 7) / std.math.pow(f64, calcium_activity, 5));
    const monocalcium_target = try phosphate_precipitation.monocalciumPhosphateEquilibriumH2po4(p.minerals.monocalcium_phosphate, calcium_activity);
    const standard = phosphate_precipitation.Kinetics{ .substrate_limit_fraction = p.kinetics.general_substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = p.kinetics.maximum_phosphate_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = p.kinetics.maximum_phosphate_precipitation_mol_per_m3_step, .phosphate_activity_coefficient = g1 };
    return .{
        .aluminum_phosphate_mol_per_m3 = try phosphateExtent(cell.aluminum_mol_per_m3, cell.h2po4_mol_p_per_m3, cell.phosphate_minerals.aluminum_phosphate_mol_per_m3, cell.h2po4_mol_p_per_m3 * g1, aluminum_target, 1, 1, standard, dynamic),
        .iron_phosphate_mol_per_m3 = try phosphateExtent(cell.iron_mol_per_m3, cell.h2po4_mol_p_per_m3, cell.phosphate_minerals.iron_phosphate_mol_per_m3, cell.h2po4_mol_p_per_m3 * g1, iron_target, 1, 1, standard, dynamic),
        .dicalcium_phosphate_mol_per_m3 = try phosphateExtent(cell.calcium_mol_per_m3, cell.hpo4_mol_p_per_m3, cell.phosphate_minerals.dicalcium_phosphate_mol_per_m3, cell.hpo4_mol_p_per_m3 * g2, dicalcium_target, 1, 1, .{ .substrate_limit_fraction = standard.substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = standard.maximum_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = standard.maximum_dissolution_mol_per_m3_step, .phosphate_activity_coefficient = g2 }, dynamic),
        .hydroxyapatite_mol_per_m3 = try phosphateExtent(cell.calcium_mol_per_m3, cell.h2po4_mol_p_per_m3, cell.phosphate_minerals.hydroxyapatite_mol_per_m3, cell.h2po4_mol_p_per_m3 * g1, apatite_target, 5, 3, .{ .substrate_limit_fraction = standard.substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = p.kinetics.maximum_apatite_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = p.kinetics.maximum_apatite_precipitation_mol_per_m3_step, .phosphate_activity_coefficient = g1 }, dynamic),
        .monocalcium_phosphate_mol_per_m3 = try phosphateExtent(cell.calcium_mol_per_m3, cell.h2po4_mol_p_per_m3, cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3, cell.h2po4_mol_p_per_m3 * g1, monocalcium_target, 1, 2, .{ .substrate_limit_fraction = standard.substrate_limit_fraction, .maximum_precipitation_mol_per_m3_step = standard.maximum_precipitation_mol_per_m3_step, .maximum_dissolution_mol_per_m3_step = p.kinetics.maximum_monocalcium_dissolution_mol_per_m3_step, .phosphate_activity_coefficient = g1 }, dynamic),
    };
}

fn phosphateExtent(cation: f64, phosphate: f64, solid: f64, activity: f64, target: f64, cation_count: f64, phosphorus_count: f64, kinetics: phosphate_precipitation.Kinetics, limit_by_cation: bool) !f64 {
    const state = phosphate_precipitation.State{ .dissolved_cation_mol_per_m3 = cation, .dissolved_phosphate_mol_p_per_m3 = phosphate, .precipitate_mol_per_m3 = solid };
    const stoichiometry = phosphate_precipitation.Stoichiometry{ .cation_mol_per_mol_precipitate = cation_count, .phosphorus_mol_per_mol_precipitate = phosphorus_count };
    return if (limit_by_cation)
        phosphate_precipitation.calculateExtent(state, activity, target, stoichiometry, kinetics)
    else
        phosphate_precipitation.calculateExtentWithPhosphateOnlyPrecipitationLimit(state, activity, target, stoichiometry, kinetics);
}

fn saltMinerals(cell: chemistry.Cell, base: chemistry.Cell, hydroxide_activity: f64, context: Context) !ledger.SaltMineralExtents {
    const p = context.parameters;
    const limit = p.kinetics.general_substrate_limit_fraction;
    const maximum = p.kinetics.maximum_phosphate_precipitation_mol_per_m3_step;
    const aluminum = try nonnegative(cell.aluminum_mol_per_m3 + base.aluminum_mol_per_m3);
    const iron = try nonnegative(cell.iron_mol_per_m3 + base.iron_mol_per_m3);
    const calcium = try nonnegative(cell.calcium_mol_per_m3 + base.calcium_mol_per_m3);
    const carbonate = try nonnegative(cell.carbonate_mol_per_m3 + base.carbonate_mol_per_m3);
    const sulfate = try nonnegative(cell.sulfate_mol_per_m3 + base.sulfate_mol_per_m3);
    const hydroxide = try nonnegative(cell.hydroxide_mol_per_m3 + base.hydroxide_mol_per_m3);
    const gibbsite = try binaryMineral(aluminum, hydroxide, cell.salt_minerals.gibbsite_mol_per_m3, aluminum * p.activity.trivalent_activity_coefficient, hydroxide_activity, 1, 3, p.minerals.gibbsite, p.activity.trivalent_activity_coefficient, limit, maximum, maximum);
    const iron_hydroxide = try binaryMineral(iron, hydroxide, cell.salt_minerals.iron_hydroxide_mol_per_m3, iron * p.activity.trivalent_activity_coefficient, hydroxide_activity, 1, 3, p.minerals.iron_hydroxide, p.activity.trivalent_activity_coefficient, limit, maximum, maximum);
    const calcite_dissolution = try mineral_precipitation.calciteDissolutionLimit(maximum, hydroxide_activity, p.kinetics.calcite_hydroxide_inhibition_constant_mol_per_m3);
    const calcite = try binaryMineral(calcium, carbonate, cell.salt_minerals.calcite_mol_per_m3, calcium * p.activity.divalent_activity_coefficient, carbonate * p.activity.divalent_activity_coefficient, 1, 1, p.minerals.calcite, p.activity.divalent_activity_coefficient, limit, maximum, calcite_dissolution);
    const calcium_after = try nonnegative(calcium - calcite);
    const gypsum = try binaryMineral(calcium_after, sulfate, cell.salt_minerals.gypsum_mol_per_m3, calcium_after * p.activity.divalent_activity_coefficient, sulfate * p.activity.divalent_activity_coefficient, 1, 1, p.minerals.gypsum, p.activity.divalent_activity_coefficient, limit, maximum, maximum);
    return .{ .gibbsite_mol_per_m3 = gibbsite, .iron_hydroxide_mol_per_m3 = iron_hydroxide, .calcite_mol_per_m3 = calcite, .gypsum_mol_per_m3 = gypsum };
}

fn binaryMineral(first: f64, second: f64, solid: f64, first_activity: f64, second_activity: f64, first_count: f64, second_count: f64, product: f64, first_coefficient: f64, limit: f64, maximum_precipitation: f64, maximum_dissolution: f64) !f64 {
    return mineral_precipitation.calculateExtent(.{ .dissolved_first_mol_per_m3 = first, .dissolved_second_mol_per_m3 = second, .solid_mol_per_m3 = solid }, first_activity, second_activity, .{ .first_mol_per_mol_solid = first_count, .second_mol_per_mol_solid = second_count }, .{ .solubility_product = product, .first_activity_coefficient = first_coefficient, .substrate_limit_fraction = limit, .maximum_precipitation_mol_per_m3_step = maximum_precipitation, .maximum_dissolution_mol_per_m3_step = maximum_dissolution });
}

fn association(free_first: f64, free_second: f64, paired: f64, first_activity: f64, second_activity: f64, paired_activity: f64, first_coefficient: f64, constant: f64, limit: f64, maximum: f64) !f64 {
    return ion_pairing.calculate(.{ .free_first_mol_per_m3 = free_first, .free_second_mol_per_m3 = free_second, .paired_mol_per_m3 = paired }, .{ .free_first_mol_per_m3 = first_activity, .free_second_mol_per_m3 = second_activity, .paired_mol_per_m3 = paired_activity, .free_first_activity_coefficient = first_coefficient }, .{ .dissociation_constant = constant, .substrate_limit_fraction = limit, .maximum_association_mol_per_m3_step = maximum });
}

fn nonnegative(value: f64) !f64 {
    if (!std.math.isFinite(value) or value < -1e-12) return error.NegativeIntermediateLitterChemistry;
    return @max(0, value);
}

fn validate(cell: chemistry.Cell, context: Context) !void {
    _ = cell;
    if (!std.math.isFinite(context.litter_mass_per_water_volume_megagrams_per_m3) or context.litter_mass_per_water_volume_megagrams_per_m3 <= 0) return error.InvalidLitterMassWaterRatio;
    inline for (@typeInfo(DissociationConstants).@"struct".fields) |field| if (!std.math.isFinite(@field(context.parameters.dissociation, field.name)) or @field(context.parameters.dissociation, field.name) <= 0) return error.InvalidLitterDissociationConstant;
    inline for (@typeInfo(MineralProducts).@"struct".fields) |field| if (!std.math.isFinite(@field(context.parameters.minerals, field.name)) or @field(context.parameters.minerals, field.name) <= 0) return error.InvalidLitterMineralProduct;
    inline for (@typeInfo(Kinetics).@"struct".fields) |field| if (!std.math.isFinite(@field(context.parameters.kinetics, field.name)) or @field(context.parameters.kinetics, field.name) < 0) return error.InvalidLitterKinetics;
    if (context.parameters.kinetics.ammonium_substrate_limit_fraction > 1 or context.parameters.kinetics.general_substrate_limit_fraction > 1 or context.parameters.kinetics.calcite_hydroxide_inhibition_constant_mol_per_m3 <= 0 or context.parameters.cation_exchange_capacity_mol_charge_per_megagram < 0) return error.InvalidLitterKinetics;
}

fn zeroExtents() ledger.ReactionExtents {
    var value: ledger.ReactionExtents = undefined;
    zeroStruct(ledger.ReactionExtents, &value);
    return value;
}

fn zeroStruct(comptime T: type, value: *T) void {
    inline for (@typeInfo(T).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .float => @field(value.*, field.name) = 0,
        .@"struct" => zeroStruct(field.type, &@field(value.*, field.name)),
        else => unreachable,
    };
}

fn unitParameters() Parameters {
    return .{
        .activity = .{ .ionic_strength_mol_per_l = 0, .monovalent_activity_coefficient = 1, .divalent_activity_coefficient = 1, .trivalent_activity_coefficient = 1, .total_ion_activity_mol_per_m3 = 0, .electrical_conductivity_dS_per_m = 0 },
        .dissociation = .{ .ammonium = 1, .carbon_dioxide = 1, .bicarbonate = 1, .h2po4 = 1, .hpo4 = 1, .carboxyl = 1 },
        .minerals = .{ .aluminum_phosphate = 1, .iron_phosphate = 1, .dicalcium_phosphate = 1, .hydroxyapatite = 1, .monocalcium_phosphate = 1, .gibbsite = 1, .iron_hydroxide = 1, .calcite = 1, .gypsum = 1, .fixed_ph_aluminum_h2po4 = 1, .fixed_ph_iron_h2po4 = 1, .fixed_ph_hydroxyapatite_h2po4 = 1 },
        .kinetics = .{ .ammonium_substrate_limit_fraction = 0.2, .general_substrate_limit_fraction = 0.2, .maximum_ammonium_association_mol_per_m3_step = 0.1, .maximum_association_mol_per_m3_step = 0.1, .maximum_phosphate_precipitation_mol_per_m3_step = 0.1, .maximum_apatite_precipitation_mol_per_m3_step = 0.1, .maximum_monocalcium_dissolution_mol_per_m3_step = 0.1, .maximum_cation_adsorption_mol_charge_per_m3_step = 0, .calcite_hydroxide_inhibition_constant_mol_per_m3 = 1 },
        .cation_exchange_capacity_mol_charge_per_megagram = 0,
        .cation_selectivity = .{ .calcium_ammonium = 1, .calcium_hydrogen = 1, .calcium_aluminum_and_iron = 1, .calcium_magnesium = 1, .calcium_sodium = 1, .calcium_potassium = 1 },
        .water_activity_product_mol2_per_m6 = 1,
        .negligible_water_ion_concentration_mol_per_m3 = 1e-32,
    };
}

test "litter rate evaluator supplies exact single-zone associations" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    var cell = state.cells[0];
    cell.hydrogen_mol_per_m3 = 1;
    cell.hydroxide_mol_per_m3 = 1;
    cell.ammonia_mol_per_m3 = 2;
    cell.ammonium_mol_per_m3 = 0.5;
    cell.hpo4_mol_p_per_m3 = 2;
    cell.h2po4_mol_p_per_m3 = 0.5;
    cell.aluminum_mol_per_m3 = 1;
    cell.iron_mol_per_m3 = 1;
    cell.calcium_mol_per_m3 = 1;
    cell.magnesium_mol_per_m3 = 1;
    cell.sodium_mol_per_m3 = 1;
    cell.potassium_mol_per_m3 = 1;
    cell.carbonate_mol_per_m3 = 1;
    cell.bicarbonate_mol_per_m3 = 1;
    cell.carbon_dioxide_mol_per_m3 = 1;
    cell.sulfate_mol_per_m3 = 1;
    const extents = try calculate(cell, .{ .parameters = unitParameters(), .litter_mass_per_water_volume_megagrams_per_m3 = 1, .dynamic_salts = false });
    try std.testing.expect(extents.ammonium_association_mol_per_m3 > 0);
    try std.testing.expect(extents.h2po4_association_mol_p_per_m3 > 0);
    try std.testing.expectEqual(@as(f64, 0), extents.bicarbonate_hydrogen_association_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0), extents.salt_minerals.calcite_mol_per_m3);
}

test "rate context binds directly to transactional litter solver" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].hydrogen_mol_per_m3 = 1;
    state.cells[0].hydroxide_mol_per_m3 = 1;
    state.cells[0].ammonia_mol_per_m3 = 1;
    state.cells[0].ammonium_mol_per_m3 = 1;
    state.cells[0].hpo4_mol_p_per_m3 = 1;
    state.cells[0].h2po4_mol_p_per_m3 = 1;
    state.cells[0].aluminum_mol_per_m3 = 1;
    state.cells[0].iron_mol_per_m3 = 1;
    state.cells[0].calcium_mol_per_m3 = 1;
    state.cells[0].magnesium_mol_per_m3 = 1;
    state.cells[0].sodium_mol_per_m3 = 1;
    state.cells[0].potassium_mol_per_m3 = 1;
    state.cells[0].carbonate_mol_per_m3 = 1;
    state.cells[0].bicarbonate_mol_per_m3 = 1;
    state.cells[0].carbon_dioxide_mol_per_m3 = 1;
    state.cells[0].sulfate_mol_per_m3 = 1;
    var parameters = unitParameters();
    parameters.kinetics.maximum_phosphate_precipitation_mol_per_m3_step = 0;
    parameters.kinetics.maximum_apatite_precipitation_mol_per_m3_step = 0;
    const context = Context{ .parameters = parameters, .litter_mass_per_water_volume_megagrams_per_m3 = 1, .dynamic_salts = false };
    const result = try solveCell(&state, 0, &context, .{});
    try std.testing.expect(result.iterations < 60);
}

test "dynamic-salt litter evaluates carbonate and sequential salt minerals" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    var cell = state.cells[0];
    cell.hydrogen_mol_per_m3 = 1;
    cell.hydroxide_mol_per_m3 = 1;
    cell.ammonia_mol_per_m3 = 1;
    cell.ammonium_mol_per_m3 = 1;
    cell.hpo4_mol_p_per_m3 = 1;
    cell.h2po4_mol_p_per_m3 = 1;
    cell.aluminum_mol_per_m3 = 2;
    cell.iron_mol_per_m3 = 2;
    cell.calcium_mol_per_m3 = 2;
    cell.magnesium_mol_per_m3 = 1;
    cell.sodium_mol_per_m3 = 1;
    cell.potassium_mol_per_m3 = 1;
    cell.carbonate_mol_per_m3 = 2;
    cell.bicarbonate_mol_per_m3 = 0.5;
    cell.carbon_dioxide_mol_per_m3 = 0.25;
    cell.sulfate_mol_per_m3 = 2;
    var parameters = unitParameters();
    parameters.minerals.gibbsite = 0.01;
    parameters.minerals.iron_hydroxide = 0.01;
    parameters.minerals.calcite = 0.01;
    parameters.minerals.gypsum = 0.01;
    parameters.kinetics.maximum_phosphate_precipitation_mol_per_m3_step = 0.01;
    parameters.kinetics.maximum_apatite_precipitation_mol_per_m3_step = 0.01;
    const extents = try calculate(cell, .{ .parameters = parameters, .litter_mass_per_water_volume_megagrams_per_m3 = 1, .dynamic_salts = true });
    try std.testing.expect(extents.carbonate_hydrogen_association_mol_per_m3 > 0);
    try std.testing.expect(extents.salt_minerals.gibbsite_mol_per_m3 > 0);
    try std.testing.expect(extents.salt_minerals.iron_hydroxide_mol_per_m3 > 0);
    try std.testing.expect(extents.salt_minerals.calcite_mol_per_m3 > 0);
    try std.testing.expect(extents.salt_minerals.gypsum_mol_per_m3 > 0);
}
