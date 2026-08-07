const std = @import("std");
const aqueous_network = @import("aqueous_network.zig");
const phosphate_network = @import("phosphate_network.zig");
const cation_exchange = @import("cation_exchange.zig");
const geochemistry = @import("geochemistry_network.zig");
const charge_classification = @import("charge_classification.zig");
const activity_coefficients = @import("activity_coefficients.zig");
const aqueous_reaction_rates = @import("aqueous_reaction_rates.zig");
const phosphate_reaction_rates = @import("phosphate_reaction_rates.zig");
const phosphate_exchange = @import("phosphate_exchange.zig");
const geochemistry_reaction_rates = @import("geochemistry_reaction_rates.zig");
const carboxyl_exchange = @import("carboxyl_exchange.zig");

pub const PhosphateTransformations = struct {
    non_band: phosphate_network.Transformations,
    band: phosphate_network.Transformations,
};

pub const ReactionParameters = struct {
    fractions: charge_classification.ZoneFractions,
    non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3: f64,
    band_phosphate_soil_mass_per_water_volume_megagrams_per_m3: f64,
    cation_exchange_capacity_mol_charge_per_megagram: f64,
    cation_exchange_water_ratios: CationExchangeWaterRatios,
    total_carboxyl_sites_mol_per_megagram: f64,
    carboxyl_exchange_parameters: carboxyl_exchange.Parameters,
    aqueous_constants: aqueous_reaction_rates.EquilibriumConstants,
    aqueous_kinetics: aqueous_reaction_rates.Kinetics,
    phosphate_constants: phosphate_reaction_rates.EquilibriumConstants,
    phosphate_surface: phosphate_exchange.Parameters,
    phosphate_minerals: ?phosphate_reaction_rates.MineralParameters,
    phosphate_kinetics: phosphate_reaction_rates.Kinetics,
    cation_exchange_parameters: cation_exchange.Parameters,
    geochemistry_products: geochemistry_reaction_rates.SolubilityProducts,
    geochemistry_kinetics: geochemistry_reaction_rates.Kinetics,
    water_activity_product_mol2_per_m6: f64,
    negligible_water_ion_concentration_mol_per_m3: f64,
};

pub const CationExchangeWaterRatios = struct {
    shared_megagrams_per_m3: f64,
    ammonium_non_band_megagrams_per_m3: f64,
    ammonium_band_megagrams_per_m3: f64,
};

pub const CellTransformations = struct {
    aqueous: aqueous_network.Transformations,
    non_band_phosphate: phosphate_network.Transformations,
    band_phosphate: phosphate_network.Transformations,
    non_band_phosphate_water_fraction: f64,
    band_phosphate_water_fraction: f64,
    cation_adsorption_mol_per_megagram: cation_exchange.Cations,
    cation_exchange_water_ratios: CationExchangeWaterRatios,
    geochemistry: geochemistry.Transformations,
    carboxyl_hydrogen_change_mol_per_megagram: f64 = 0,
    carboxyl_soil_mass_per_water_volume_megagrams_per_m3: f64 = 0,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    aqueous: []aqueous_network.State,
    non_band_phosphate: []phosphate_network.State,
    band_phosphate: []phosphate_network.State,
    water_mol_per_m3: []f64,
    cation_exchange_mol_per_megagram: []cation_exchange.Cations,
    carboxyl_bound_hydrogen_mol_per_megagram: []f64,
    geochemistry_solids: []geochemistry.SolidState,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroChemistryCellCount;
        const aqueous = try allocator.alloc(aqueous_network.State, cell_count);
        errdefer allocator.free(aqueous);
        const non_band = try allocator.alloc(phosphate_network.State, cell_count);
        errdefer allocator.free(non_band);
        const band = try allocator.alloc(phosphate_network.State, cell_count);
        errdefer allocator.free(band);
        const water = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(water);
        const exchange_state = try allocator.alloc(cation_exchange.Cations, cell_count);
        errdefer allocator.free(exchange_state);
        const carboxyl_bound_hydrogen = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(carboxyl_bound_hydrogen);
        const geochemistry_solids = try allocator.alloc(geochemistry.SolidState, cell_count);
        errdefer allocator.free(geochemistry_solids);
        for (aqueous) |*cell| zeroStruct(aqueous_network.State, cell);
        for (non_band) |*cell| zeroStruct(phosphate_network.State, cell);
        for (band) |*cell| zeroStruct(phosphate_network.State, cell);
        @memset(water, 0);
        for (exchange_state) |*cell| zeroStruct(cation_exchange.Cations, cell);
        @memset(carboxyl_bound_hydrogen, 0);
        for (geochemistry_solids) |*cell| zeroStruct(geochemistry.SolidState, cell);
        return .{ .allocator = allocator, .cell_count = cell_count, .aqueous = aqueous, .non_band_phosphate = non_band, .band_phosphate = band, .water_mol_per_m3 = water, .cation_exchange_mol_per_megagram = exchange_state, .carboxyl_bound_hydrogen_mol_per_megagram = carboxyl_bound_hydrogen, .geochemistry_solids = geochemistry_solids };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.geochemistry_solids);
        self.allocator.free(self.carboxyl_bound_hydrogen_mol_per_megagram);
        self.allocator.free(self.cation_exchange_mol_per_megagram);
        self.allocator.free(self.water_mol_per_m3);
        self.allocator.free(self.band_phosphate);
        self.allocator.free(self.non_band_phosphate);
        self.allocator.free(self.aqueous);
        self.* = undefined;
    }

    /// Dimensionally safe replacement for SOLUTE.F lines 2295--2452. A
    /// chemistry iterate is atomic across the shared aqueous system and both
    /// fertilizer zones; no partial state survives a failed invariant check.
    pub fn commitCell(self: *State, cell_index: usize, transformations: CellTransformations) !void {
        if (cell_index >= self.cell_count) return error.ChemistryCellIndexOutOfBounds;
        var staged_aqueous = self.aqueous[cell_index];
        var staged_non_band = self.non_band_phosphate[cell_index];
        var staged_band = self.band_phosphate[cell_index];
        var staged_exchange = self.cation_exchange_mol_per_megagram[cell_index];
        const staged_carboxyl_hydrogen = self.carboxyl_bound_hydrogen_mol_per_megagram[cell_index] + transformations.carboxyl_hydrogen_change_mol_per_megagram;
        var staged_geochemistry = self.geochemistry_solids[cell_index];
        if (!std.math.isFinite(staged_carboxyl_hydrogen) or staged_carboxyl_hydrogen < -1e-12) return error.InvalidCarboxylExchangeState;
        if (!std.math.isFinite(transformations.carboxyl_soil_mass_per_water_volume_megagrams_per_m3) or transformations.carboxyl_soil_mass_per_water_volume_megagrams_per_m3 < 0) return error.InvalidCarboxylExchangeSoilWaterRatio;
        if (!validFraction(transformations.non_band_phosphate_water_fraction) or !validFraction(transformations.band_phosphate_water_fraction)) return error.InvalidPhosphateWaterFraction;
        var aqueous_changes = transformations.aqueous;
        aqueous_changes.hydrogen -= transformations.carboxyl_hydrogen_change_mol_per_megagram * transformations.carboxyl_soil_mass_per_water_volume_megagrams_per_m3;
        try addCationExchange(&staged_exchange, &aqueous_changes, transformations.cation_adsorption_mol_per_megagram, transformations.cation_exchange_water_ratios);
        addGeochemistry(&aqueous_changes, transformations.geochemistry);
        addSharedPhosphate(&aqueous_changes, transformations.non_band_phosphate, transformations.non_band_phosphate_water_fraction);
        addSharedPhosphate(&aqueous_changes, transformations.band_phosphate, transformations.band_phosphate_water_fraction);
        const staged_water = self.water_mol_per_m3[cell_index] +
            transformations.non_band_phosphate.water_mol_per_m3 * transformations.non_band_phosphate_water_fraction +
            transformations.band_phosphate.water_mol_per_m3 * transformations.band_phosphate_water_fraction;
        if (!std.math.isFinite(staged_water) or staged_water < -1e-12) return error.InvalidChemistryWaterState;
        try aqueous_network.commit(&staged_aqueous, aqueous_changes);
        try phosphate_network.commit(&staged_non_band, transformations.non_band_phosphate);
        try phosphate_network.commit(&staged_band, transformations.band_phosphate);
        try geochemistry.commitSolids(&staged_geochemistry, transformations.geochemistry);
        self.aqueous[cell_index] = staged_aqueous;
        self.non_band_phosphate[cell_index] = staged_non_band;
        self.band_phosphate[cell_index] = staged_band;
        self.water_mol_per_m3[cell_index] = staged_water;
        self.cation_exchange_mol_per_megagram[cell_index] = staged_exchange;
        self.carboxyl_bound_hydrogen_mol_per_megagram[cell_index] = @max(0.0, staged_carboxyl_hydrogen);
        self.geochemistry_solids[cell_index] = staged_geochemistry;
    }

    /// Returns the complete shared-aqueous change assembled from every
    /// reaction family before the atomic cell commit.
    pub fn assembledAqueousChanges(
        transformations: CellTransformations,
        current_exchange: cation_exchange.Cations,
    ) !aqueous_network.Transformations {
        var aqueous_changes = transformations.aqueous;
        aqueous_changes.hydrogen -=
            transformations.carboxyl_hydrogen_change_mol_per_megagram *
            transformations.carboxyl_soil_mass_per_water_volume_megagrams_per_m3;
        var exchange = current_exchange;
        try addCationExchange(
            &exchange,
            &aqueous_changes,
            transformations.cation_adsorption_mol_per_megagram,
            transformations.cation_exchange_water_ratios,
        );
        addGeochemistry(&aqueous_changes, transformations.geochemistry);
        addSharedPhosphate(
            &aqueous_changes,
            transformations.non_band_phosphate,
            transformations.non_band_phosphate_water_fraction,
        );
        addSharedPhosphate(
            &aqueous_changes,
            transformations.band_phosphate,
            transformations.band_phosphate_water_fraction,
        );
        return aqueous_changes;
    }

    /// Direct source-order `RH2O` diagnostic for SOLUTE.F lines 2225--2226.
    /// Production state mutation remains unchanged until restart and protected
    /// coupled-solver comparisons include this additional water owner.
    pub fn sourceOrderWaterChangeMolPerM3(
        transformations: CellTransformations,
    ) !f64 {
        if (!validFraction(transformations.non_band_phosphate_water_fraction) or
            !validFraction(transformations.band_phosphate_water_fraction))
            return error.InvalidPhosphateWaterFraction;
        const change = transformations.aqueous.carbon_dioxide +
            transformations.non_band_phosphate.water_mol_per_m3 *
                transformations.non_band_phosphate_water_fraction +
            transformations.band_phosphate.water_mol_per_m3 *
                transformations.band_phosphate_water_fraction;
        if (!std.math.isFinite(change)) return error.InvalidChemistryWaterState;
        return change;
    }

    pub fn evaluateCarboxylHydrogenChange(
        self: *const State,
        cell_index: usize,
        total_carboxyl_sites_mol_per_megagram: f64,
        hydrogen_activity_mol_per_m3: f64,
        soil_mass_per_water_volume_megagrams_per_m3: f64,
        parameters: carboxyl_exchange.Parameters,
    ) !f64 {
        if (cell_index >= self.cell_count) return error.ChemistryCellIndexOutOfBounds;
        return carboxyl_exchange.calculateChangeMolPerMg(.{
            .total_carboxyl_sites_mol_per_megagram = total_carboxyl_sites_mol_per_megagram,
            .hydrogen_occupied_sites_mol_per_megagram = @min(
                self.carboxyl_bound_hydrogen_mol_per_megagram[cell_index],
                total_carboxyl_sites_mol_per_megagram,
            ),
            .hydrogen_activity_mol_per_m3 = hydrogen_activity_mol_per_m3,
            .soil_mass_per_water_volume_megagrams_per_m3 = soil_mass_per_water_volume_megagrams_per_m3,
        }, parameters);
    }

    pub fn packedComponentCount() usize {
        return @typeInfo(aqueous_network.State).@"struct".fields.len +
            2 * @typeInfo(phosphate_network.State).@"struct".fields.len +
            @typeInfo(cation_exchange.Cations).@"struct".fields.len +
            @typeInfo(geochemistry.SolidState).@"struct".fields.len + 2;
    }

    pub fn packedComponentName(index: usize) ?[]const u8 {
        var cursor: usize = 0;
        inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
            if (index == cursor) return "aqueous." ++ field.name;
            cursor += 1;
        }
        inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
            if (index == cursor) return "phosphate_non_band." ++ field.name;
            cursor += 1;
        }
        inline for (@typeInfo(phosphate_network.State).@"struct".fields) |field| {
            if (index == cursor) return "phosphate_band." ++ field.name;
            cursor += 1;
        }
        inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field| {
            if (index == cursor) return "cation_exchange." ++ field.name;
            cursor += 1;
        }
        if (index == cursor) return "carboxyl_bound_hydrogen_mol_per_megagram";
        cursor += 1;
        inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field| {
            if (index == cursor) return "geochemistry_solids." ++ field.name;
            cursor += 1;
        }
        if (index == cursor) return "water_mol_per_m3";
        return null;
    }

    /// Packs one runtime cell into the vector consumed by the coupled hybrid
    /// solver. The chemistry species list is scientific structure; grid size,
    /// layers, zones, and solver storage remain runtime allocated.
    pub fn packCell(self: *const State, cell_index: usize, output: []f64) !void {
        if (cell_index >= self.cell_count) return error.ChemistryCellIndexOutOfBounds;
        if (output.len != packedComponentCount()) return error.ChemistryVectorSizeMismatch;
        var cursor: usize = 0;
        packStruct(aqueous_network.State, self.aqueous[cell_index], output, &cursor);
        packStruct(phosphate_network.State, self.non_band_phosphate[cell_index], output, &cursor);
        packStruct(phosphate_network.State, self.band_phosphate[cell_index], output, &cursor);
        packStruct(cation_exchange.Cations, self.cation_exchange_mol_per_megagram[cell_index], output, &cursor);
        output[cursor] = self.carboxyl_bound_hydrogen_mol_per_megagram[cell_index];
        cursor += 1;
        packStruct(geochemistry.SolidState, self.geochemistry_solids[cell_index], output, &cursor);
        output[cursor] = self.water_mol_per_m3[cell_index];
    }

    pub fn unpackCell(self: *State, cell_index: usize, input: []const f64) !void {
        if (cell_index >= self.cell_count) return error.ChemistryCellIndexOutOfBounds;
        if (input.len != packedComponentCount()) return error.ChemistryVectorSizeMismatch;
        var cursor: usize = 0;
        const aqueous = try unpackStruct(aqueous_network.State, input, &cursor);
        const non_band = try unpackStruct(phosphate_network.State, input, &cursor);
        const band = try unpackStruct(phosphate_network.State, input, &cursor);
        const exchange_state = try unpackStruct(cation_exchange.Cations, input, &cursor);
        const carboxyl_bound_hydrogen = input[cursor];
        cursor += 1;
        if (!std.math.isFinite(carboxyl_bound_hydrogen) or carboxyl_bound_hydrogen < -1e-12) return error.InvalidCarboxylExchangeState;
        const geochemistry_solids = try unpackStruct(geochemistry.SolidState, input, &cursor);
        const water = input[cursor];
        if (!std.math.isFinite(water) or water < -1e-12) return error.InvalidChemistryWaterState;
        self.aqueous[cell_index] = aqueous;
        self.non_band_phosphate[cell_index] = non_band;
        self.band_phosphate[cell_index] = band;
        self.cation_exchange_mol_per_megagram[cell_index] = exchange_state;
        self.carboxyl_bound_hydrogen_mol_per_megagram[cell_index] = carboxyl_bound_hydrogen;
        self.geochemistry_solids[cell_index] = geochemistry_solids;
        self.water_mol_per_m3[cell_index] = water;
    }

    pub fn activityCoefficients(self: *const State, cell_index: usize, fractions: charge_classification.ZoneFractions) !activity_coefficients.Result {
        if (cell_index >= self.cell_count) return error.ChemistryCellIndexOutOfBounds;
        const totals_per_m3 = try charge_classification.classify(self.aqueous[cell_index], self.non_band_phosphate[cell_index], self.band_phosphate[cell_index], fractions);
        return activity_coefficients.calculate(totals_per_m3, 1);
    }

    pub fn evaluateAqueousTransformations(self: *const State, cell_index: usize, fractions: charge_classification.ZoneFractions, constants: aqueous_reaction_rates.EquilibriumConstants, kinetics: aqueous_reaction_rates.Kinetics) !aqueous_network.Transformations {
        const coefficients = try self.activityCoefficients(cell_index, fractions);
        const fluxes = try aqueous_reaction_rates.calculateSourceOrder(
            self.aqueous[cell_index],
            coefficients,
            constants,
            kinetics,
            .{
                .non_band = if (fractions.ammonium_non_band > 0) .wet else .dry,
                .band = if (fractions.ammonium_band > 0) .wet else .dry,
            },
        );
        return aqueous_network.assemble(fluxes, .{ .non_band = fractions.ammonium_non_band, .band = fractions.ammonium_band });
    }

    pub fn evaluatePhosphateTransformations(self: *const State, cell_index: usize, fractions: charge_classification.ZoneFractions, non_band_soil_mass_per_water_volume_megagrams_per_m3: f64, band_soil_mass_per_water_volume_megagrams_per_m3: f64, constants: phosphate_reaction_rates.EquilibriumConstants, surface_parameters: phosphate_exchange.Parameters, mineral_parameters: ?phosphate_reaction_rates.MineralParameters, kinetics: phosphate_reaction_rates.Kinetics) !PhosphateTransformations {
        const coefficients = try self.activityCoefficients(cell_index, fractions);
        const non_band = if (fractions.phosphate_non_band > 0) blk: {
            const fluxes = try phosphate_reaction_rates.calculate(self.aqueous[cell_index], self.non_band_phosphate[cell_index], coefficients, non_band_soil_mass_per_water_volume_megagrams_per_m3, constants, surface_parameters, mineral_parameters, kinetics);
            break :blk try phosphate_network.assemble(fluxes);
        } else std.mem.zeroes(phosphate_network.Transformations);
        const band = if (fractions.phosphate_band > 0) blk: {
            const fluxes = try phosphate_reaction_rates.calculate(self.aqueous[cell_index], self.band_phosphate[cell_index], coefficients, band_soil_mass_per_water_volume_megagrams_per_m3, constants, surface_parameters, mineral_parameters, kinetics);
            break :blk try phosphate_network.assemble(fluxes);
        } else std.mem.zeroes(phosphate_network.Transformations);
        return .{ .non_band = non_band, .band = band };
    }

    pub fn evaluateGeochemistryTransformations(self: *const State, cell_index: usize, fractions: charge_classification.ZoneFractions, products: geochemistry_reaction_rates.SolubilityProducts, kinetics: geochemistry_reaction_rates.Kinetics) !geochemistry.Transformations {
        const coefficients = try self.activityCoefficients(cell_index, fractions);
        return geochemistry_reaction_rates.calculate(self.aqueous[cell_index], self.geochemistry_solids[cell_index], coefficients, products, kinetics);
    }

    pub fn evaluateCell(self: *const State, cell_index: usize, parameters: ReactionParameters) !CellTransformations {
        const coefficients = try self.activityCoefficients(cell_index, parameters.fractions);
        const aqueous = try self.evaluateAqueousTransformations(cell_index, parameters.fractions, parameters.aqueous_constants, parameters.aqueous_kinetics);
        const phosphate = try self.evaluatePhosphateTransformations(cell_index, parameters.fractions, parameters.non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3, parameters.band_phosphate_soil_mass_per_water_volume_megagrams_per_m3, parameters.phosphate_constants, parameters.phosphate_surface, parameters.phosphate_minerals, parameters.phosphate_kinetics);
        const shared = self.aqueous[cell_index];
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
        const activities = cation_exchange.Cations{
            .ammonium_non_band = shared.ammonium_non_band * coefficients.monovalent_activity_coefficient,
            .ammonium_band = shared.ammonium_band * coefficients.monovalent_activity_coefficient,
            .hydrogen = shared.hydrogen * coefficients.monovalent_activity_coefficient,
            .aluminum = shared.aluminum * coefficients.trivalent_activity_coefficient,
            .iron = shared.iron * coefficients.trivalent_activity_coefficient,
            .calcium = shared.calcium * coefficients.divalent_activity_coefficient,
            .magnesium = shared.magnesium * coefficients.divalent_activity_coefficient,
            .sodium = shared.sodium * coefficients.monovalent_activity_coefficient,
            .potassium = shared.potassium * coefficients.monovalent_activity_coefficient,
        };
        const adsorption = try cation_exchange.calculateSourceOrder(.{
            .cation_exchange_capacity_mol_charge_per_megagram = parameters.cation_exchange_capacity_mol_charge_per_megagram,
            .aqueous_concentration_mol_per_m3 = concentrations,
            .aqueous_activity_mol_per_m3 = activities,
            .exchange_concentration_mol_per_megagram = self.cation_exchange_mol_per_megagram[cell_index],
            .ammonium_non_band_fraction = parameters.fractions.ammonium_non_band,
            .ammonium_band_fraction = parameters.fractions.ammonium_band,
            .soil_mass_per_water_volume_megagrams_per_m3 = parameters.cation_exchange_water_ratios.shared_megagrams_per_m3,
        }, parameters.cation_exchange_parameters, .{
            .minimum_activity_mol_per_m3 = parameters.negligible_water_ion_concentration_mol_per_m3,
        });
        const carboxyl_hydrogen_change = try self.evaluateCarboxylHydrogenChange(
            cell_index,
            parameters.total_carboxyl_sites_mol_per_megagram,
            activities.hydrogen,
            parameters.cation_exchange_water_ratios.shared_megagrams_per_m3,
            parameters.carboxyl_exchange_parameters,
        );
        return .{
            .aqueous = aqueous,
            .non_band_phosphate = phosphate.non_band,
            .band_phosphate = phosphate.band,
            .non_band_phosphate_water_fraction = parameters.fractions.phosphate_non_band,
            .band_phosphate_water_fraction = parameters.fractions.phosphate_band,
            .cation_adsorption_mol_per_megagram = adsorption,
            .cation_exchange_water_ratios = parameters.cation_exchange_water_ratios,
            .geochemistry = try self.evaluateGeochemistryTransformations(cell_index, parameters.fractions, parameters.geochemistry_products, parameters.geochemistry_kinetics),
            .carboxyl_hydrogen_change_mol_per_megagram = carboxyl_hydrogen_change,
            .carboxyl_soil_mass_per_water_volume_megagrams_per_m3 = parameters.cation_exchange_water_ratios.shared_megagrams_per_m3,
        };
    }
};

fn zeroStruct(comptime T: type, value: *T) void {
    inline for (@typeInfo(T).@"struct".fields) |field| @field(value.*, field.name) = 0;
}

fn filledStruct(comptime T: type, amount: f64) T {
    var value: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(value, field.name) = amount;
    return value;
}

test "inactive ammonium zone accepts only zero exchange flux" {
    var exchange_state = filledStruct(cation_exchange.Cations, 0);
    var aqueous = filledStruct(aqueous_network.Transformations, 0);
    var adsorption = filledStruct(cation_exchange.Cations, 0);
    const ratios: CationExchangeWaterRatios = .{
        .shared_megagrams_per_m3 = 1.5,
        .ammonium_non_band_megagrams_per_m3 = 1.5,
        .ammonium_band_megagrams_per_m3 = 0,
    };
    try addCationExchange(
        &exchange_state,
        &aqueous,
        adsorption,
        ratios,
    );
    adsorption.ammonium_band = 1.0e-6;
    try std.testing.expectError(
        error.InactiveAmmoniumZoneHasExchangeFlux,
        addCationExchange(&exchange_state, &aqueous, adsorption, ratios),
    );
}

fn packStruct(comptime T: type, value: T, output: []f64, cursor: *usize) void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        output[cursor.*] = @field(value, field.name);
        cursor.* += 1;
    }
}

fn unpackStruct(comptime T: type, input: []const f64, cursor: *usize) !T {
    var value: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| {
        const component = input[cursor.*];
        if (!std.math.isFinite(component)) return error.NonFiniteChemistryVector;
        if (component < -1e-12) return error.NegativeChemistryVector;
        @field(value, field.name) = component;
        cursor.* += 1;
    }
    return value;
}

fn validFraction(value: f64) bool {
    return std.math.isFinite(value) and value >= 0 and value <= 1;
}

fn addSharedPhosphate(aqueous: *aqueous_network.Transformations, phosphate: phosphate_network.Transformations, water_fraction: f64) void {
    aqueous.aluminum += phosphate.dissolved_aluminum_mol_per_m3 * water_fraction;
    aqueous.iron += phosphate.dissolved_iron_mol_per_m3 * water_fraction;
    aqueous.calcium += phosphate.dissolved_calcium_mol_per_m3 * water_fraction;
    aqueous.magnesium += phosphate.dissolved_magnesium_mol_per_m3 * water_fraction;
    aqueous.hydrogen += phosphate.dissolved_hydrogen_mol_per_m3 * water_fraction;
    aqueous.hydroxide += phosphate.dissolved_hydroxide_mol_per_m3 * water_fraction;
}

fn addCationExchange(exchange_state: *cation_exchange.Cations, aqueous: *aqueous_network.Transformations, adsorption: cation_exchange.Cations, ratios: CationExchangeWaterRatios) !void {
    if (!std.math.isFinite(ratios.shared_megagrams_per_m3) or
        ratios.shared_megagrams_per_m3 <= 0 or
        !std.math.isFinite(ratios.ammonium_non_band_megagrams_per_m3) or
        ratios.ammonium_non_band_megagrams_per_m3 < 0 or
        !std.math.isFinite(ratios.ammonium_band_megagrams_per_m3) or
        ratios.ammonium_band_megagrams_per_m3 < 0)
        return error.InvalidCationExchangeWaterRatio;
    // SOLUTE sets BKVLWH/BKVLWB to zero for an inactive zero-water zone.
    // Such a zone must also contribute exactly zero adsorption; accepting a
    // nonzero flux here would silently discard or create aqueous ammonium.
    if ((ratios.ammonium_non_band_megagrams_per_m3 == 0 and
        adsorption.ammonium_non_band != 0) or
        (ratios.ammonium_band_megagrams_per_m3 == 0 and
            adsorption.ammonium_band != 0))
        return error.InactiveAmmoniumZoneHasExchangeFlux;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field| {
        const flux = @field(adsorption, field.name);
        if (!std.math.isFinite(flux)) return error.NonFiniteCationExchangeTransformation;
        @field(exchange_state.*, field.name) += flux;
        if (@field(exchange_state.*, field.name) < -1e-12) return error.NegativeCationExchangeState;
    }
    aqueous.ammonium_non_band -= adsorption.ammonium_non_band * ratios.ammonium_non_band_megagrams_per_m3;
    aqueous.ammonium_band -= adsorption.ammonium_band * ratios.ammonium_band_megagrams_per_m3;
    aqueous.hydrogen -= adsorption.hydrogen * ratios.shared_megagrams_per_m3;
    aqueous.aluminum -= adsorption.aluminum * ratios.shared_megagrams_per_m3;
    aqueous.iron -= adsorption.iron * ratios.shared_megagrams_per_m3;
    aqueous.calcium -= adsorption.calcium * ratios.shared_megagrams_per_m3;
    aqueous.magnesium -= adsorption.magnesium * ratios.shared_megagrams_per_m3;
    aqueous.sodium -= adsorption.sodium * ratios.shared_megagrams_per_m3;
    aqueous.potassium -= adsorption.potassium * ratios.shared_megagrams_per_m3;
}

fn addGeochemistry(aqueous: *aqueous_network.Transformations, transformations: geochemistry.Transformations) void {
    aqueous.aluminum += transformations.dissolved_aluminum_mol_per_m3;
    aqueous.iron += transformations.dissolved_iron_mol_per_m3;
    aqueous.calcium += transformations.dissolved_calcium_mol_per_m3;
    aqueous.magnesium += transformations.dissolved_magnesium_mol_per_m3;
    aqueous.sodium += transformations.dissolved_sodium_mol_per_m3;
    aqueous.potassium += transformations.dissolved_potassium_mol_per_m3;
    aqueous.hydrogen += transformations.dissolved_hydrogen_mol_per_m3;
    aqueous.hydroxide += transformations.dissolved_hydroxide_mol_per_m3;
    aqueous.carbonate += transformations.dissolved_carbonate_mol_per_m3;
    aqueous.sulfate += transformations.dissolved_sulfate_mol_per_m3;
    aqueous.hydrogen_silicate += transformations.dissolved_hydrogen_silicate_mol_per_m3;
}

test "chemistry state uses runtime cell allocation" {
    var state = try State.init(std.testing.allocator, 7);
    defer state.deinit();
    try std.testing.expectEqual(@as(usize, 7), state.aqueous.len);
    try std.testing.expectEqual(@as(usize, 7), state.non_band_phosphate.len);
    try std.testing.expectEqual(@as(usize, 7), state.band_phosphate.len);
    try std.testing.expectEqual(@as(usize, 7), state.water_mol_per_m3.len);
    try std.testing.expectEqual(@as(usize, 7), state.cation_exchange_mol_per_megagram.len);
    try std.testing.expectEqual(@as(usize, 7), state.geochemistry_solids.len);
}

test "cell chemistry commit is atomic across aqueous and both phosphate zones" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.aqueous[1] = filledStruct(aqueous_network.State, 1);
    state.non_band_phosphate[1] = filledStruct(phosphate_network.State, 1);
    state.band_phosphate[1] = filledStruct(phosphate_network.State, 1);
    const before_aqueous = state.aqueous[1];
    const before_non_band = state.non_band_phosphate[1];
    const before_band = state.band_phosphate[1];
    var aqueous_change = filledStruct(aqueous_network.Transformations, 0);
    var non_band_change = filledStruct(phosphate_network.Transformations, 0);
    var band_change = filledStruct(phosphate_network.Transformations, 0);
    aqueous_change.calcium = 0.1;
    non_band_change.dissolved_h2po4_mol_p_per_m3 = 0.1;
    band_change.dissolved_h2po4_mol_p_per_m3 = -2;
    try std.testing.expectError(error.NegativePhosphateNetworkState, state.commitCell(1, .{ .aqueous = aqueous_change, .non_band_phosphate = non_band_change, .band_phosphate = band_change, .non_band_phosphate_water_fraction = 0.8, .band_phosphate_water_fraction = 0.2, .cation_adsorption_mol_per_megagram = filledStruct(cation_exchange.Cations, 0), .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 1, .ammonium_non_band_megagrams_per_m3 = 1, .ammonium_band_megagrams_per_m3 = 1 }, .geochemistry = filledStruct(geochemistry.Transformations, 0) }));
    try std.testing.expectEqualDeep(before_aqueous, state.aqueous[1]);
    try std.testing.expectEqualDeep(before_non_band, state.non_band_phosphate[1]);
    try std.testing.expectEqualDeep(before_band, state.band_phosphate[1]);
}

test "runtime chemistry cell round trips through hybrid solver vector" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    state.aqueous[2] = filledStruct(aqueous_network.State, 0.25);
    state.non_band_phosphate[2] = filledStruct(phosphate_network.State, 0.5);
    state.band_phosphate[2] = filledStruct(phosphate_network.State, 0.75);
    const vector = try std.testing.allocator.alloc(f64, State.packedComponentCount());
    defer std.testing.allocator.free(vector);
    try state.packCell(2, vector);
    state.aqueous[2] = filledStruct(aqueous_network.State, 0);
    state.non_band_phosphate[2] = filledStruct(phosphate_network.State, 0);
    state.band_phosphate[2] = filledStruct(phosphate_network.State, 0);
    try state.unpackCell(2, vector);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.aqueous[2].calcium, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), state.non_band_phosphate[2].dissolved_h2po4_mol_p_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), state.band_phosphate[2].hydroxyapatite_solid_mol_per_m3, 1e-15);
}

test "band and non-band phosphate reactions update one fraction-weighted shared state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 10);
    state.water_mol_per_m3[0] = 10;
    const aqueous_change = filledStruct(aqueous_network.Transformations, 0);
    var non_band_change = filledStruct(phosphate_network.Transformations, 0);
    var band_change = filledStruct(phosphate_network.Transformations, 0);
    non_band_change.dissolved_calcium_mol_per_m3 = -1;
    band_change.dissolved_calcium_mol_per_m3 = -1;
    non_band_change.water_mol_per_m3 = 2;
    band_change.water_mol_per_m3 = 4;
    try state.commitCell(0, .{
        .aqueous = aqueous_change,
        .non_band_phosphate = non_band_change,
        .band_phosphate = band_change,
        .non_band_phosphate_water_fraction = 0.75,
        .band_phosphate_water_fraction = 0.25,
        .cation_adsorption_mol_per_megagram = filledStruct(cation_exchange.Cations, 0),
        .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 1, .ammonium_non_band_megagrams_per_m3 = 1, .ammonium_band_megagrams_per_m3 = 1 },
        .geochemistry = filledStruct(geochemistry.Transformations, 0),
    });
    try std.testing.expectApproxEqAbs(@as(f64, 9), state.aqueous[0].calcium, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 12.5), state.water_mol_per_m3[0], 1e-15);
}

test "cation adsorption uses distinct shared and ammonium water ratios" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 10);
    var adsorption = filledStruct(cation_exchange.Cations, 0);
    adsorption.ammonium_non_band = 0.1;
    adsorption.ammonium_band = 0.2;
    adsorption.calcium = 0.3;
    try state.commitCell(0, .{
        .aqueous = filledStruct(aqueous_network.Transformations, 0),
        .non_band_phosphate = filledStruct(phosphate_network.Transformations, 0),
        .band_phosphate = filledStruct(phosphate_network.Transformations, 0),
        .non_band_phosphate_water_fraction = 0.8,
        .band_phosphate_water_fraction = 0.2,
        .cation_adsorption_mol_per_megagram = adsorption,
        .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 2, .ammonium_non_band_megagrams_per_m3 = 3, .ammonium_band_megagrams_per_m3 = 4 },
        .geochemistry = filledStruct(geochemistry.Transformations, 0),
    });
    try std.testing.expectApproxEqAbs(@as(f64, 9.7), state.aqueous[0].ammonium_non_band, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.2), state.aqueous[0].ammonium_band, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.4), state.aqueous[0].calcium, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), state.cation_exchange_mol_per_megagram[0].calcium, 1e-15);
}

test "carboxyl protonation conserves hydrogen across solid and aqueous pools" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].hydrogen = 1;
    const zero_aqueous = filledStruct(aqueous_network.Transformations, 0);
    const zero_phosphate = filledStruct(phosphate_network.Transformations, 0);
    try state.commitCell(0, .{
        .aqueous = zero_aqueous,
        .non_band_phosphate = zero_phosphate,
        .band_phosphate = zero_phosphate,
        .non_band_phosphate_water_fraction = 0.5,
        .band_phosphate_water_fraction = 0.5,
        .cation_adsorption_mol_per_megagram = filledStruct(cation_exchange.Cations, 0),
        .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 2, .ammonium_non_band_megagrams_per_m3 = 2, .ammonium_band_megagrams_per_m3 = 2 },
        .geochemistry = filledStruct(geochemistry.Transformations, 0),
        .carboxyl_hydrogen_change_mol_per_megagram = 0.2,
        .carboxyl_soil_mass_per_water_volume_megagrams_per_m3 = 2,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.carboxyl_bound_hydrogen_mol_per_megagram[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), state.aqueous[0].hydrogen, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.aqueous[0].hydrogen + 2 * state.carboxyl_bound_hydrogen_mol_per_megagram[0], 1e-15);
}

test "shared geochemistry updates aqueous and solid pools in one transaction" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 10);
    state.geochemistry_solids[0] = filledStruct(geochemistry.SolidState, 1);
    const changes = try geochemistry.assemble(
        .{ .gibbsite_precipitation_mol_per_m3 = 0.1, .iron_hydroxide_precipitation_mol_per_m3 = 0, .calcite_precipitation_mol_per_m3 = 0, .gypsum_precipitation_mol_per_m3 = 0 },
        .{ .aluminum_natural_mol_per_m3 = 0.02, .aluminum_ground_mol_per_m3 = 0, .iron_natural_mol_per_m3 = 0, .iron_ground_mol_per_m3 = 0, .calcium_natural_mol_per_m3 = 0, .calcium_ground_mol_per_m3 = 0, .magnesium_natural_mol_per_m3 = 0, .magnesium_ground_mol_per_m3 = 0, .sodium_natural_mol_per_m3 = 0, .sodium_ground_mol_per_m3 = 0, .potassium_natural_mol_per_m3 = 0, .potassium_ground_mol_per_m3 = 0 },
    );
    try state.commitCell(0, .{
        .aqueous = filledStruct(aqueous_network.Transformations, 0),
        .non_band_phosphate = filledStruct(phosphate_network.Transformations, 0),
        .band_phosphate = filledStruct(phosphate_network.Transformations, 0),
        .non_band_phosphate_water_fraction = 0.8,
        .band_phosphate_water_fraction = 0.2,
        .cation_adsorption_mol_per_megagram = filledStruct(cation_exchange.Cations, 0),
        .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 1, .ammonium_non_band_megagrams_per_m3 = 1, .ammonium_band_megagrams_per_m3 = 1 },
        .geochemistry = changes,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 9.92), state.aqueous[0].aluminum, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.1), state.geochemistry_solids[0].gibbsite_solid_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.98), state.geochemistry_solids[0].aluminum_natural_silicate_mol_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 10.015), state.aqueous[0].hydrogen_silicate, 1e-15);
}

test "runtime chemistry state feeds activity coefficients directly" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].calcium = 10;
    state.aqueous[0].chloride = 20;
    const result = try state.activityCoefficients(0, .{ .ammonium_non_band = 0.8, .ammonium_band = 0.2, .nitrate_non_band = 0.6, .nitrate_band = 0.4, .phosphate_non_band = 0.7, .phosphate_band = 0.3 });
    try std.testing.expect(result.ionic_strength_mol_per_l > 0);
    try std.testing.expect(result.divalent_activity_coefficient < result.monovalent_activity_coefficient);
}

test "runtime state evaluates named aqueous reactions into conservative changes" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 1);
    state.aqueous[0].calcium = 2;
    const transformations = try state.evaluateAqueousTransformations(0, .{ .ammonium_non_band = 0.8, .ammonium_band = 0.2, .nitrate_non_band = 0.6, .nitrate_band = 0.4, .phosphate_non_band = 0.7, .phosphate_band = 0.3 }, filledStruct(aqueous_reaction_rates.EquilibriumConstants, 1), .{ .ammonium_substrate_limit_fraction = 0.2, .general_substrate_limit_fraction = 0.2, .maximum_fast_association_mol_per_m3_step = 0.1, .maximum_slow_association_mol_per_m3_step = 0.1 });
    try std.testing.expectApproxEqAbs(@as(f64, 0), transformations.calcium + transformations.calcium_hydroxide + transformations.calcium_carbonate + transformations.calcium_bicarbonate + transformations.calcium_sulfate, 1e-14);
}

test "runtime aqueous evaluation honors independent ammonium water zones" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 1);
    var constants = filledStruct(aqueous_reaction_rates.EquilibriumConstants, 1);
    constants.ammonium = 0.5;
    const transformations = try state.evaluateAqueousTransformations(
        0,
        .{
            .ammonium_non_band = 0,
            .ammonium_band = 1,
            .nitrate_non_band = 1,
            .nitrate_band = 0,
            .phosphate_non_band = 1,
            .phosphate_band = 0,
        },
        constants,
        .{
            .ammonium_substrate_limit_fraction = 0.2,
            .general_substrate_limit_fraction = 0.2,
            .maximum_fast_association_mol_per_m3_step = 0.1,
            .maximum_slow_association_mol_per_m3_step = 0.1,
        },
    );
    try std.testing.expectEqual(@as(f64, 0), transformations.ammonium_non_band);
    try std.testing.expect(transformations.ammonium_band > 0);
}

test "runtime state evaluates both phosphate zones with shared activities" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 1);
    state.non_band_phosphate[0] = filledStruct(phosphate_network.State, 1);
    state.band_phosphate[0] = filledStruct(phosphate_network.State, 1);
    const changes = try state.evaluatePhosphateTransformations(0, .{ .ammonium_non_band = 0.8, .ammonium_band = 0.2, .nitrate_non_band = 0.6, .nitrate_band = 0.4, .phosphate_non_band = 0.7, .phosphate_band = 0.3 }, 1.2, 1.5, filledStruct(phosphate_reaction_rates.EquilibriumConstants, 1), .{ .protonated_site_equilibrium_constant = 1, .hydroxyl_site_equilibrium_constant = 1, .h2po4_exchange_equilibrium_constant = 1, .hpo4_exchange_equilibrium_constant = 1, .water_activity_product_mol2_per_m6 = 1, .h2po4_dissociation_constant = 1, .maximum_exchange_mol_per_megagram_step = 0.1, .substrate_limit_fraction = 0.2 }, null, .{ .substrate_limit_fraction = 0.2, .maximum_pairing_mol_per_m3_step = 0.1 });
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field| {
        try std.testing.expect(std.math.isFinite(@field(changes.non_band, field.name)));
        try std.testing.expect(std.math.isFinite(@field(changes.band, field.name)));
    }
}

test "runtime state evaluates shared mineral and weathering rates" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 1);
    state.geochemistry_solids[0] = filledStruct(geochemistry.SolidState, 1);
    const changes = try state.evaluateGeochemistryTransformations(0, .{ .ammonium_non_band = 0.8, .ammonium_band = 0.2, .nitrate_non_band = 0.6, .nitrate_band = 0.4, .phosphate_non_band = 0.7, .phosphate_band = 0.3 }, filledStruct(geochemistry_reaction_rates.SolubilityProducts, 1), .{ .general_substrate_limit_fraction = 0.2, .hydrogen_coupled_substrate_limit_fraction = 0.2, .maximum_hydroxide_mineral_mol_per_m3_step = 0.1, .maximum_general_mineral_mol_per_m3_step = 0.1, .calcite_hydroxide_inhibition_constant_mol_per_m3 = 1, .maximum_natural_weathering_mol_per_m3_step = 0.01, .maximum_ground_weathering_mol_per_m3_step = 0.02 });
    try std.testing.expectApproxEqAbs(@as(f64, 0), changes.dissolved_calcium_mol_per_m3 + changes.calcite_solid_mol_per_m3 + changes.gypsum_solid_mol_per_m3 + changes.calcium_natural_silicate_mol_per_m3 + changes.calcium_ground_silicate_mol_per_m3, 1e-14);
}

test "single cell evaluation assembles every active SOLUTE reaction family" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0] = filledStruct(aqueous_network.State, 1);
    state.non_band_phosphate[0] = filledStruct(phosphate_network.State, 1);
    state.band_phosphate[0] = filledStruct(phosphate_network.State, 1);
    state.cation_exchange_mol_per_megagram[0] = filledStruct(cation_exchange.Cations, 1);
    state.geochemistry_solids[0] = filledStruct(geochemistry.SolidState, 1);
    const transformations = try state.evaluateCell(0, .{
        .fractions = .{ .ammonium_non_band = 0.8, .ammonium_band = 0.2, .nitrate_non_band = 0.6, .nitrate_band = 0.4, .phosphate_non_band = 0.7, .phosphate_band = 0.3 },
        .non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = 1.2,
        .band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = 1.5,
        .cation_exchange_capacity_mol_charge_per_megagram = 10,
        .cation_exchange_water_ratios = .{ .shared_megagrams_per_m3 = 1.4, .ammonium_non_band_megagrams_per_m3 = 1.1, .ammonium_band_megagrams_per_m3 = 1.8 },
        .total_carboxyl_sites_mol_per_megagram = 0,
        .carboxyl_exchange_parameters = .{ .dissociation_constant_mol_per_m3 = 0.01, .maximum_exchange_mol_per_m3_per_iteration = 0.01, .substrate_limit_fraction_per_iteration = 0.2 },
        .aqueous_constants = filledStruct(aqueous_reaction_rates.EquilibriumConstants, 1),
        .aqueous_kinetics = .{ .ammonium_substrate_limit_fraction = 0.2, .general_substrate_limit_fraction = 0.2, .maximum_fast_association_mol_per_m3_step = 0.01, .maximum_slow_association_mol_per_m3_step = 0.01 },
        .phosphate_constants = filledStruct(phosphate_reaction_rates.EquilibriumConstants, 1),
        .phosphate_surface = .{ .protonated_site_equilibrium_constant = 1, .hydroxyl_site_equilibrium_constant = 1, .h2po4_exchange_equilibrium_constant = 1, .hpo4_exchange_equilibrium_constant = 1, .water_activity_product_mol2_per_m6 = 1, .h2po4_dissociation_constant = 1, .maximum_exchange_mol_per_megagram_step = 0.01, .substrate_limit_fraction = 0.2 },
        .phosphate_minerals = null,
        .phosphate_kinetics = .{ .substrate_limit_fraction = 0.2, .maximum_pairing_mol_per_m3_step = 0.01 },
        .cation_exchange_parameters = .{ .selectivity = .{ .calcium_ammonium = 1, .calcium_hydrogen = 1, .calcium_aluminum_and_iron = 1, .calcium_magnesium = 1, .calcium_sodium = 1, .calcium_potassium = 1 }, .substrate_limit_fraction = 0.2, .maximum_adsorption_mol_charge_per_m3_step = 0.01 },
        .geochemistry_products = filledStruct(geochemistry_reaction_rates.SolubilityProducts, 1),
        .geochemistry_kinetics = .{ .general_substrate_limit_fraction = 0.2, .hydrogen_coupled_substrate_limit_fraction = 0.2, .maximum_hydroxide_mineral_mol_per_m3_step = 0.01, .maximum_general_mineral_mol_per_m3_step = 0.01, .calcite_hydroxide_inhibition_constant_mol_per_m3 = 1, .maximum_natural_weathering_mol_per_m3_step = 0.001, .maximum_ground_weathering_mol_per_m3_step = 0.002 },
        .water_activity_product_mol2_per_m6 = 1,
        .negligible_water_ion_concentration_mol_per_m3 = 1e-32,
    });
    inline for (@typeInfo(aqueous_network.Transformations).@"struct".fields) |field| try std.testing.expect(std.math.isFinite(@field(transformations.aqueous, field.name)));
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field| try std.testing.expect(std.math.isFinite(@field(transformations.cation_adsorption_mol_per_megagram, field.name)));
}

test "source-order water change includes carbon dioxide hydration" {
    var transformations = std.mem.zeroes(CellTransformations);
    transformations.aqueous.carbon_dioxide = 0.2;
    transformations.non_band_phosphate.water_mol_per_m3 = 0.3;
    transformations.band_phosphate.water_mol_per_m3 = 0.4;
    transformations.non_band_phosphate_water_fraction = 0.75;
    transformations.band_phosphate_water_fraction = 0.25;
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.525),
        try State.sourceOrderWaterChangeMolPerM3(transformations),
        1e-15,
    );
}
