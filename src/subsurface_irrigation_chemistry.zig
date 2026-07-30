const std = @import("std");
const routing = @import("irrigation_layer_routing.zig");
const aqueous_transport = @import("solute_transport.zig");
const aqueous_species = @import("solute_transport_species.zig");
const chemistry = @import("solute_chemistry_state.zig");
const mineral_nitrogen = @import("mineral_nitrogen_transport.zig");
const nutrient_speciation = @import("precipitation_nutrient_speciation.zig");
const aqueous_rates = @import("solute_aqueous_reaction_rates.zig");
const phosphate_rates = @import("solute_phosphate_reaction_rates.zig");

pub const ElementMolarMassesGPerMol = struct {
    nitrogen: f64,
    phosphorus: f64,
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
    sulfur: f64,
    chloride: f64,
};

pub const EquilibriumConstants = struct {
    aqueous: aqueous_rates.EquilibriumConstants,
    phosphate: phosphate_rates.EquilibriumConstants,
};

pub const Parameters = struct {
    molar_mass_g_per_mol: ElementMolarMassesGPerMol,
    equilibrium: ?EquilibriumConstants,
    ammonium_band_fraction: f64,
    nitrate_band_fraction: f64,
    phosphate_band_fraction: f64,
};

pub const BoundaryInput = struct {
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
    ion_mol: f64,
};

/// Adds free hydrogen and the eight irrigation salt carriers to the selected
/// runtime soil layers before the same-hour 50-species aqueous solve.
pub fn addTransportedIons(
    loads: *const routing.Loads,
    transport: *aqueous_transport.State,
    parameters: Parameters,
) !void {
    try validate(loads, transport.cell_count, parameters);
    if (transport.species_count != aqueous_species.AqueousSpecies.count)
        return error.SubsurfaceIrrigationTransportDimensionMismatch;

    for (0..loads.subsurface_water_m3.len) |layer| {
        const first = layer * routing.dissolved_species_count;
        const mass = loads.subsurface_dissolved_mass_g[first..][0..routing.dissolved_species_count];
        const additions = ionAdditions(
            loads.subsurface_hydrogen_mol[layer],
            mass,
            parameters.molar_mass_g_per_mol,
        );
        inline for (transported_ion_species, additions) |species, addition| {
            const destination = layer * aqueous_species.AqueousSpecies.count +
                aqueous_species.index(species);
            const candidate = transport.amount_mol[destination] + addition;
            if (!std.math.isFinite(candidate) or candidate < 0)
                return error.InvalidSubsurfaceIrrigationChemistryTransaction;
        }
    }
    for (0..loads.subsurface_water_m3.len) |layer| {
        const first = layer * routing.dissolved_species_count;
        const additions = ionAdditions(
            loads.subsurface_hydrogen_mol[layer],
            loads.subsurface_dissolved_mass_g[first..][0..routing.dissolved_species_count],
            parameters.molar_mass_g_per_mol,
        );
        inline for (transported_ion_species, additions) |species, addition|
            transport.amount_mol[
                layer * aqueous_species.AqueousSpecies.count +
                    aqueous_species.index(species)
            ] += addition;
    }
}

/// HPO4 and H2PO4 are reaction-state species, not members of the conservative
/// 50-carrier vector. Publish them after aqueous transport has restored the
/// selected layer's chemistry concentrations.
pub fn addPhosphate(
    loads: *const routing.Loads,
    state: *chemistry.State,
    matrix_water_volume_m3: []const f64,
    parameters: Parameters,
) !void {
    try validate(loads, state.cell_count, parameters);
    if (matrix_water_volume_m3.len != state.cell_count)
        return error.SubsurfaceIrrigationChemistryDimensionMismatch;
    for (0..state.cell_count) |layer| {
        const species = try layerNutrients(loads, layer, parameters);
        const phosphorus_mol = species[3] + species[4];
        const water_m3 = matrix_water_volume_m3[layer];
        if (!std.math.isFinite(water_m3) or water_m3 < 0 or
            (phosphorus_mol > 0 and water_m3 == 0))
            return error.InvalidSubsurfaceIrrigationChemistryWater;
        if (phosphorus_mol == 0) continue;
        const non_band = 1 - parameters.phosphate_band_fraction;
        inline for (.{
            state.non_band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 +
                species[3] * non_band / water_m3,
            state.non_band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 +
                species[4] * non_band / water_m3,
            state.band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 +
                species[3] * parameters.phosphate_band_fraction / water_m3,
            state.band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 +
                species[4] * parameters.phosphate_band_fraction / water_m3,
        }) |candidate| if (!std.math.isFinite(candidate) or candidate < 0)
            return error.InvalidSubsurfaceIrrigationChemistryTransaction;
    }
    for (0..state.cell_count) |layer| {
        const species = try layerNutrients(loads, layer, parameters);
        if (species[3] + species[4] == 0) continue;
        const water_m3 = matrix_water_volume_m3[layer];
        const non_band = 1 - parameters.phosphate_band_fraction;
        state.non_band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 += species[3] * non_band / water_m3;
        state.non_band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 += species[4] * non_band / water_m3;
        state.band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 += species[3] * parameters.phosphate_band_fraction / water_m3;
        state.band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 += species[4] * parameters.phosphate_band_fraction / water_m3;
    }
}

/// Adds NH4, NH3, and NO3 amounts after the reaction concentrations have been
/// imported into the mineral-N transport owner and before its same-hour solve.
pub fn addMineralNitrogen(
    loads: *const routing.Loads,
    state: *mineral_nitrogen.State,
    parameters: Parameters,
) !void {
    try validate(loads, state.cell_count, parameters);
    if (state.matrix.species_count != mineral_nitrogen.species_count)
        return error.SubsurfaceIrrigationMineralNitrogenDimensionMismatch;
    for (0..state.cell_count) |layer| {
        const species = try layerNutrients(loads, layer, parameters);
        const additions = nitrogenAdditions(species, parameters);
        for (additions, 0..) |addition, component| {
            const candidate = state.matrix.amount_mol[
                layer * mineral_nitrogen.species_count + component
            ] + addition;
            if (!std.math.isFinite(candidate) or candidate < 0)
                return error.InvalidSubsurfaceIrrigationChemistryTransaction;
        }
    }
    for (0..state.cell_count) |layer| {
        const additions = nitrogenAdditions(
            try layerNutrients(loads, layer, parameters),
            parameters,
        );
        for (additions, 0..) |addition, component|
            state.matrix.amount_mol[
                layer * mineral_nitrogen.species_count + component
            ] += addition;
    }
}

/// Exact extensive boundary input corresponding to the committed subsurface
/// irrigation carriers. REDIST keeps aqueous N and P in their own balances;
/// its `SBU` ion boundary contains free H and the salt-system carriers only.
pub fn boundaryInput(
    loads: *const routing.Loads,
    parameters: Parameters,
) !BoundaryInput {
    try validate(loads, loads.subsurface_water_m3.len, parameters);
    var result: BoundaryInput = .{
        .nitrogen_g_n = 0,
        .phosphorus_g_p = 0,
        .ion_mol = 0,
    };
    for (0..loads.subsurface_water_m3.len) |layer| {
        const first = layer * routing.dissolved_species_count;
        const mass = loads.subsurface_dissolved_mass_g[first..][0..routing.dissolved_species_count];
        const ions = ionAdditions(
            loads.subsurface_hydrogen_mol[layer],
            mass,
            parameters.molar_mass_g_per_mol,
        );
        result.nitrogen_g_n += mass[0] + mass[1];
        result.phosphorus_g_p += mass[2];
        for (ions) |amount| result.ion_mol += amount;
        if (!std.math.isFinite(result.nitrogen_g_n) or
            !std.math.isFinite(result.phosphorus_g_p) or
            !std.math.isFinite(result.ion_mol))
            return error.InvalidSubsurfaceIrrigationChemistryTransaction;
    }
    return result;
}

const transported_ion_species = [_]aqueous_species.AqueousSpecies{
    .hydrogen, .aluminum,  .iron,    .calcium,  .magnesium,
    .sodium,   .potassium, .sulfate, .chloride,
};

fn ionAdditions(
    hydrogen_mol: f64,
    mass_g: []const f64,
    molar_mass: ElementMolarMassesGPerMol,
) [transported_ion_species.len]f64 {
    return .{
        hydrogen_mol,
        mass_g[3] / molar_mass.aluminum,
        mass_g[4] / molar_mass.iron,
        mass_g[5] / molar_mass.calcium,
        mass_g[6] / molar_mass.magnesium,
        mass_g[7] / molar_mass.sodium,
        mass_g[8] / molar_mass.potassium,
        mass_g[9] / molar_mass.sulfur,
        mass_g[10] / molar_mass.chloride,
    };
}

fn nitrogenAdditions(species: [5]f64, parameters: Parameters) [mineral_nitrogen.species_count]f64 {
    const ammonium_non_band = 1 - parameters.ammonium_band_fraction;
    const nitrate_non_band = 1 - parameters.nitrate_band_fraction;
    return .{
        species[0] * ammonium_non_band,
        species[0] * parameters.ammonium_band_fraction,
        species[1] * ammonium_non_band,
        species[1] * parameters.ammonium_band_fraction,
        species[2] * nitrate_non_band,
        species[2] * parameters.nitrate_band_fraction,
        0,
        0,
    };
}

fn layerNutrients(
    loads: *const routing.Loads,
    layer: usize,
    parameters: Parameters,
) ![5]f64 {
    const first = layer * routing.dissolved_species_count;
    const mass = loads.subsurface_dissolved_mass_g[first..][0..routing.dissolved_species_count];
    const water_m3 = loads.subsurface_water_m3[layer];
    if (water_m3 == 0) return .{ 0, 0, 0, 0, 0 };
    const hydrogen_mol_per_m3 = loads.subsurface_hydrogen_mol[layer] / water_m3;
    const ph = -std.math.log10(@max(hydrogen_mol_per_m3 / 1000.0, 1.0e-14));
    if (parameters.equilibrium) |equilibrium|
        return nutrient_speciation.calculate(.{
            .ph = ph,
            .ammonium_g_n_per_m3 = mass[0] / water_m3,
            .nitrate_g_n_per_m3 = mass[1] / water_m3,
            .phosphate_g_p_per_m3 = mass[2] / water_m3,
            .nitrogen_g_per_mol = parameters.molar_mass_g_per_mol.nitrogen,
            .phosphorus_g_per_mol = parameters.molar_mass_g_per_mol.phosphorus,
        }, equilibrium.aqueous, equilibrium.phosphate);
    return .{
        mass[0] / parameters.molar_mass_g_per_mol.nitrogen,
        0,
        mass[1] / parameters.molar_mass_g_per_mol.nitrogen,
        0,
        mass[2] / parameters.molar_mass_g_per_mol.phosphorus,
    };
}

fn validate(
    loads: *const routing.Loads,
    expected_layer_count: usize,
    parameters: Parameters,
) !void {
    const layer_count = try std.math.mul(
        usize,
        loads.cell_count,
        loads.soil_layer_capacity,
    );
    if (layer_count != expected_layer_count or
        loads.subsurface_water_m3.len != layer_count or
        loads.subsurface_hydrogen_mol.len != layer_count or
        loads.subsurface_dissolved_mass_g.len !=
            try std.math.mul(usize, layer_count, routing.dissolved_species_count))
        return error.SubsurfaceIrrigationChemistryDimensionMismatch;
    inline for (std.meta.fields(ElementMolarMassesGPerMol)) |field| {
        const value = @field(parameters.molar_mass_g_per_mol, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSubsurfaceIrrigationMolarMass;
    }
    inline for (.{
        parameters.ammonium_band_fraction,
        parameters.nitrate_band_fraction,
        parameters.phosphate_band_fraction,
    }) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidSubsurfaceIrrigationBandFraction;
    for (loads.subsurface_water_m3, loads.subsurface_hydrogen_mol) |water, hydrogen|
        if (!std.math.isFinite(water) or water < 0 or
            !std.math.isFinite(hydrogen) or hydrogen < 0 or
            (water == 0 and hydrogen != 0))
            return error.InvalidSubsurfaceIrrigationChemistryLoad;
    for (loads.subsurface_dissolved_mass_g) |mass|
        if (!std.math.isFinite(mass) or mass < 0)
            return error.InvalidSubsurfaceIrrigationChemistryLoad;
}

fn testParameters() Parameters {
    return .{
        .molar_mass_g_per_mol = .{
            .nitrogen = 14,
            .phosphorus = 31,
            .aluminum = 27,
            .iron = 56,
            .calcium = 40,
            .magnesium = 24,
            .sodium = 23,
            .potassium = 39,
            .sulfur = 32,
            .chloride = 35.5,
        },
        .equilibrium = null,
        .ammonium_band_fraction = 0.25,
        .nitrate_band_fraction = 0.4,
        .phosphate_band_fraction = 0.2,
    };
}

test "runtime-depth subsurface carriers enter only selected layers and conserve mass" {
    var loads = try routing.Loads.init(std.testing.allocator, 2, 7);
    defer loads.deinit();
    try loads.accumulate(1, 7, &.{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }, 10, 0.02, 0.45, .{
        .ph = 6,
        .ammonium_nitrogen = 14,
        .nitrate_nitrogen = 28,
        .phosphate_phosphorus = 31,
        .aluminum = 27,
        .iron = 56,
        .calcium = 40,
        .magnesium = 24,
        .sodium = 23,
        .potassium = 39,
        .sulfate_sulfur = 32,
        .chloride = 35.5,
    });
    var transport = try aqueous_transport.State.init(
        std.testing.allocator,
        14,
        aqueous_species.AqueousSpecies.count,
    );
    defer transport.deinit();
    try addTransportedIons(&loads, &transport, testParameters());
    const selected: usize = 7 + 4;
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.2),
        transport.amount_mol[
            selected * aqueous_species.AqueousSpecies.count +
                aqueous_species.index(.aluminum)
        ],
        1e-12,
    );
    for (0..14) |layer| if (layer != selected)
        try std.testing.expectEqual(
            @as(f64, 0),
            transport.amount_mol[
                layer * aqueous_species.AqueousSpecies.count +
                    aqueous_species.index(.aluminum)
            ],
        );
    const boundary = try boundaryInput(&loads, testParameters());
    try std.testing.expectApproxEqAbs(@as(f64, 8.4), boundary.nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 6.2), boundary.phosphorus_g_p, 1e-12);
    // 0.2 mol each of the eight supplied salt carriers plus pH-derived H.
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.6) + loads.subsurface_hydrogen_mol[selected],
        boundary.ion_mol,
        1e-12,
    );
}

test "late invalid load leaves transported ions unchanged" {
    var loads = try routing.Loads.init(std.testing.allocator, 1, 2);
    defer loads.deinit();
    loads.subsurface_water_m3[0] = 1;
    loads.subsurface_hydrogen_mol[0] = 0.001;
    loads.subsurface_dissolved_mass_g[3] = 27;
    loads.subsurface_dissolved_mass_g[2 * routing.dissolved_species_count - 1] = std.math.nan(f64);
    var transport = try aqueous_transport.State.init(
        std.testing.allocator,
        2,
        aqueous_species.AqueousSpecies.count,
    );
    defer transport.deinit();
    try std.testing.expectError(
        error.InvalidSubsurfaceIrrigationChemistryLoad,
        addTransportedIons(&loads, &transport, testParameters()),
    );
    for (transport.amount_mol) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "subsurface nutrient binding preserves runtime band splits and selected layer" {
    var loads = try routing.Loads.init(std.testing.allocator, 1, 7);
    defer loads.deinit();
    try loads.accumulate(0, 7, &.{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }, 1, 1, 0.61, .{
        .ph = 7,
        .ammonium_nitrogen = 14,
        .nitrate_nitrogen = 28,
        .phosphate_phosphorus = 31,
        .aluminum = 0,
        .iron = 0,
        .calcium = 0,
        .magnesium = 0,
        .sodium = 0,
        .potassium = 0,
        .sulfate_sulfur = 0,
        .chloride = 0,
    });
    var nitrogen = try mineral_nitrogen.State.init(std.testing.allocator, 7);
    defer nitrogen.deinit();
    try addMineralNitrogen(&loads, &nitrogen, testParameters());
    const selected: usize = 6;
    const first = selected * mineral_nitrogen.species_count;
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), nitrogen.matrix.amount_mol[first + @intFromEnum(mineral_nitrogen.Species.ammonium_non_band)], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), nitrogen.matrix.amount_mol[first + @intFromEnum(mineral_nitrogen.Species.ammonium_band)], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), nitrogen.matrix.amount_mol[first + @intFromEnum(mineral_nitrogen.Species.nitrate_non_band)], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), nitrogen.matrix.amount_mol[first + @intFromEnum(mineral_nitrogen.Species.nitrate_band)], 1e-12);

    var chemistry_state = try chemistry.State.init(std.testing.allocator, 7);
    defer chemistry_state.deinit();
    const water = [_]f64{1} ** 7;
    try addPhosphate(&loads, &chemistry_state, &water, testParameters());
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.8),
        chemistry_state.non_band_phosphate[selected].dissolved_h2po4_mol_p_per_m3,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.2),
        chemistry_state.band_phosphate[selected].dissolved_h2po4_mol_p_per_m3,
        1e-12,
    );
    for (0..selected) |layer| {
        try std.testing.expectEqual(
            @as(f64, 0),
            chemistry_state.non_band_phosphate[layer].dissolved_h2po4_mol_p_per_m3,
        );
        try std.testing.expectEqual(
            @as(f64, 0),
            nitrogen.matrix.amount_mol[
                layer * mineral_nitrogen.species_count +
                    @intFromEnum(mineral_nitrogen.Species.ammonium_non_band)
            ],
        );
    }
}

test "phosphate zero-water failure is atomic across runtime layers" {
    var loads = try routing.Loads.init(std.testing.allocator, 1, 2);
    defer loads.deinit();
    loads.subsurface_water_m3[1] = 1;
    loads.subsurface_hydrogen_mol[1] = 0.0001;
    loads.subsurface_dissolved_mass_g[
        routing.dissolved_species_count + 2
    ] = 31;
    var state = try chemistry.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 = 3;
    try std.testing.expectError(
        error.InvalidSubsurfaceIrrigationChemistryWater,
        addPhosphate(&loads, &state, &.{ 1, 0 }, testParameters()),
    );
    try std.testing.expectEqual(
        @as(f64, 3),
        state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        state.non_band_phosphate[1].dissolved_h2po4_mol_p_per_m3,
    );
}
