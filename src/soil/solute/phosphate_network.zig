const std = @import("std");
const exchange = @import("phosphate_exchange.zig");

pub const MineralFluxes = struct {
    aluminum_phosphate_mol_per_m3: f64,
    iron_phosphate_mol_per_m3: f64,
    dicalcium_phosphate_mol_per_m3: f64,
    hydroxyapatite_mol_per_m3: f64,
    monocalcium_phosphate_mol_per_m3: f64,
};

pub const DissociationAndPairingFluxes = struct {
    po4_hydrogen_association_mol_p_per_m3: f64,
    hpo4_hydrogen_association_mol_p_per_m3: f64,
    h2po4_hydrogen_association_mol_p_per_m3: f64,
    iron_hpo4_pairing_mol_p_per_m3: f64,
    iron_h2po4_pairing_mol_p_per_m3: f64,
    calcium_po4_pairing_mol_p_per_m3: f64,
    calcium_hpo4_pairing_mol_p_per_m3: f64,
    calcium_h2po4_pairing_mol_p_per_m3: f64,
    magnesium_hpo4_pairing_mol_p_per_m3: f64,
};

pub const Fluxes = struct {
    minerals: MineralFluxes,
    surface: exchange.Flux,
    aqueous: DissociationAndPairingFluxes,
    soil_mass_per_water_volume_megagrams_per_m3: f64,
};

pub const Transformations = struct {
    dissolved_po4_mol_p_per_m3: f64,
    dissolved_hpo4_mol_p_per_m3: f64,
    dissolved_h2po4_mol_p_per_m3: f64,
    dissolved_h3po4_mol_p_per_m3: f64,
    deprotonated_site_mol_per_megagram: f64,
    hydroxyl_site_mol_per_megagram: f64,
    protonated_site_mol_per_megagram: f64,
    adsorbed_hpo4_mol_p_per_megagram: f64,
    adsorbed_h2po4_mol_p_per_megagram: f64,
    dissolved_aluminum_mol_per_m3: f64,
    dissolved_iron_mol_per_m3: f64,
    dissolved_calcium_mol_per_m3: f64,
    dissolved_magnesium_mol_per_m3: f64,
    dissolved_hydrogen_mol_per_m3: f64,
    dissolved_hydroxide_mol_per_m3: f64,
    water_mol_per_m3: f64,
    aluminum_phosphate_solid_mol_per_m3: f64,
    iron_phosphate_solid_mol_per_m3: f64,
    dicalcium_phosphate_solid_mol_per_m3: f64,
    hydroxyapatite_solid_mol_per_m3: f64,
    monocalcium_phosphate_solid_mol_per_m3: f64,
    iron_hpo4_pair_mol_per_m3: f64,
    iron_h2po4_pair_mol_per_m3: f64,
    calcium_po4_pair_mol_per_m3: f64,
    calcium_hpo4_pair_mol_per_m3: f64,
    calcium_h2po4_pair_mol_per_m3: f64,
    magnesium_hpo4_pair_mol_per_m3: f64,
};

pub const State = struct {
    dissolved_po4_mol_p_per_m3: f64,
    dissolved_hpo4_mol_p_per_m3: f64,
    dissolved_h2po4_mol_p_per_m3: f64,
    dissolved_h3po4_mol_p_per_m3: f64,
    deprotonated_site_mol_per_megagram: f64,
    hydroxyl_site_mol_per_megagram: f64,
    protonated_site_mol_per_megagram: f64,
    adsorbed_hpo4_mol_p_per_megagram: f64,
    adsorbed_h2po4_mol_p_per_megagram: f64,
    aluminum_phosphate_solid_mol_per_m3: f64,
    iron_phosphate_solid_mol_per_m3: f64,
    dicalcium_phosphate_solid_mol_per_m3: f64,
    hydroxyapatite_solid_mol_per_m3: f64,
    monocalcium_phosphate_solid_mol_per_m3: f64,
    iron_hpo4_pair_mol_per_m3: f64,
    iron_h2po4_pair_mol_per_m3: f64,
    calcium_po4_pair_mol_per_m3: f64,
    calcium_hpo4_pair_mol_per_m3: f64,
    calcium_h2po4_pair_mol_per_m3: f64,
    magnesium_hpo4_pair_mol_per_m3: f64,
};

pub const SourceIterationStage = enum {
    before_iteration_ceiling,
    iteration_ceiling,
};

/// Applies every phosphate-network destination atomically. This is the state
/// transaction used by a chemistry iterate; failed validation changes nothing.
pub fn commit(state: *State, transformations: Transformations) !void {
    try validateState(state.*);
    var next = state.*;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const change = @field(transformations, field.name);
        if (!std.math.isFinite(change)) return error.NonFinitePhosphateNetworkTransformation;
        @field(next, field.name) += change;
    }
    try validateState(next);
    state.* = next;
}

/// Direct update for one zone in SOLUTE.F 2368--2387. The caller invokes it
/// for non-band then band state, preserving source order. Every listed pool is
/// floor-clamped and the entire block is skipped at `M == MRXN`.
pub fn applySourceOrderAqueousUpdate(
    current: State,
    transformations: Transformations,
    minimum_concentration_mol_per_m3: f64,
    stage: SourceIterationStage,
) !State {
    try validateState(current);
    if (!std.math.isFinite(minimum_concentration_mol_per_m3) or
        minimum_concentration_mol_per_m3 <= 0)
        return error.InvalidPhosphateConcentrationFloor;
    inline for (@typeInfo(Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(transformations, field.name)))
            return error.NonFinitePhosphateNetworkTransformation;
    if (stage == .iteration_ceiling) return current;

    var next = current;
    inline for (.{
        "dissolved_po4_mol_p_per_m3",
        "dissolved_hpo4_mol_p_per_m3",
        "dissolved_h2po4_mol_p_per_m3",
        "dissolved_h3po4_mol_p_per_m3",
        "iron_hpo4_pair_mol_per_m3",
        "iron_h2po4_pair_mol_per_m3",
        "calcium_po4_pair_mol_per_m3",
        "calcium_hpo4_pair_mol_per_m3",
        "calcium_h2po4_pair_mol_per_m3",
        "magnesium_hpo4_pair_mol_per_m3",
    }) |name| {
        @field(next, name) = @max(
            minimum_concentration_mol_per_m3,
            @field(current, name) + @field(transformations, name),
        );
    }
    return next;
}

/// Direct surface-site update for one zone in SOLUTE.F 2402--2411. Invoke for
/// non-band then band state. These five pools are not floor-clamped and the
/// block is skipped at the iteration ceiling.
pub fn applySourceOrderSurfaceUpdate(
    current: State,
    transformations: Transformations,
    stage: SourceIterationStage,
) !State {
    try validateState(current);
    inline for (@typeInfo(Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(transformations, field.name)))
            return error.NonFinitePhosphateNetworkTransformation;
    if (stage == .iteration_ceiling) return current;

    var next = current;
    inline for (.{
        "deprotonated_site_mol_per_megagram",
        "hydroxyl_site_mol_per_megagram",
        "protonated_site_mol_per_megagram",
        "adsorbed_hpo4_mol_p_per_megagram",
        "adsorbed_h2po4_mol_p_per_megagram",
    }) |name| {
        @field(next, name) = @field(current, name) +
            @field(transformations, name);
        if (!std.math.isFinite(@field(next, name)) or @field(next, name) < 0)
            return error.InvalidPhosphateSurfaceStateUpdate;
    }
    return next;
}

/// Direct precipitate update for one phosphate zone in SOLUTE.F 2442--2451.
/// Invoke for non-band then band state. No floor is applied and the whole
/// block is skipped at the iteration ceiling.
pub fn applySourceOrderMineralUpdate(
    current: State,
    transformations: Transformations,
    stage: SourceIterationStage,
) !State {
    try validateState(current);
    inline for (@typeInfo(Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(transformations, field.name)))
            return error.NonFinitePhosphateNetworkTransformation;
    if (stage == .iteration_ceiling) return current;

    var next = current;
    inline for (.{
        "aluminum_phosphate_solid_mol_per_m3",
        "iron_phosphate_solid_mol_per_m3",
        "dicalcium_phosphate_solid_mol_per_m3",
        "hydroxyapatite_solid_mol_per_m3",
        "monocalcium_phosphate_solid_mol_per_m3",
    }) |name| {
        @field(next, name) = @field(current, name) +
            @field(transformations, name);
        if (!std.math.isFinite(@field(next, name)) or @field(next, name) < 0)
            return error.InvalidPhosphateMineralStateUpdate;
    }
    return next;
}

/// Dimensionally corrected phosphate portion of SOLUTE.F lines 2143--2226.
/// Flux signs follow the reaction kernels: positive means association,
/// adsorption, or precipitation; negative means the reverse reaction. Surface
/// extents are mol/Mg and are explicitly converted to mol/m3; the source
/// transformation block omits that required conversion.
pub fn assemble(fluxes: Fluxes) !Transformations {
    try validate(fluxes);
    return assembleWithSiteFactor(
        fluxes,
        fluxes.soil_mass_per_water_volume_megagrams_per_m3,
    );
}

/// Exact phosphate-zone portion of SOLUTE.F 2205--2226. The source adds
/// mol Mg-1 surface extents directly to mol m-3 aqueous transformations;
/// factor 1 intentionally retains that dimensional defect for comparison.
pub fn assembleSourceOrder(fluxes: Fluxes) !Transformations {
    try validate(fluxes);
    return assembleWithSiteFactor(fluxes, 1);
}

fn assembleWithSiteFactor(
    fluxes: Fluxes,
    site_to_water: f64,
) Transformations {
    const m = fluxes.minerals;
    const s = fluxes.surface;
    const a = fluxes.aqueous;
    return .{
        .dissolved_po4_mol_p_per_m3 = -a.po4_hydrogen_association_mol_p_per_m3 - a.calcium_po4_pairing_mol_p_per_m3,
        .dissolved_hpo4_mol_p_per_m3 = -s.hpo4_with_hydroxyl_site_mol_p_per_megagram * site_to_water + a.po4_hydrogen_association_mol_p_per_m3 - a.hpo4_hydrogen_association_mol_p_per_m3 - a.iron_hpo4_pairing_mol_p_per_m3 - a.calcium_hpo4_pairing_mol_p_per_m3 - a.magnesium_hpo4_pairing_mol_p_per_m3,
        .dissolved_h2po4_mol_p_per_m3 = -m.aluminum_phosphate_mol_per_m3 - m.iron_phosphate_mol_per_m3 - m.dicalcium_phosphate_mol_per_m3 - 3 * m.hydroxyapatite_mol_per_m3 - 2 * m.monocalcium_phosphate_mol_per_m3 - (s.h2po4_with_protonated_site_mol_p_per_megagram + s.h2po4_with_hydroxyl_site_mol_p_per_megagram) * site_to_water + a.hpo4_hydrogen_association_mol_p_per_m3 - a.h2po4_hydrogen_association_mol_p_per_m3 - a.iron_h2po4_pairing_mol_p_per_m3 - a.calcium_h2po4_pairing_mol_p_per_m3,
        .dissolved_h3po4_mol_p_per_m3 = a.h2po4_hydrogen_association_mol_p_per_m3,
        .deprotonated_site_mol_per_megagram = -s.hydroxyl_to_deprotonated_site_mol_per_megagram,
        .hydroxyl_site_mol_per_megagram = s.hydroxyl_to_deprotonated_site_mol_per_megagram - s.protonated_to_hydroxyl_site_mol_per_megagram - s.h2po4_with_hydroxyl_site_mol_p_per_megagram - s.hpo4_with_hydroxyl_site_mol_p_per_megagram,
        .protonated_site_mol_per_megagram = s.protonated_to_hydroxyl_site_mol_per_megagram - s.h2po4_with_protonated_site_mol_p_per_megagram,
        .adsorbed_hpo4_mol_p_per_megagram = s.hpo4_with_hydroxyl_site_mol_p_per_megagram,
        .adsorbed_h2po4_mol_p_per_megagram = s.h2po4_with_protonated_site_mol_p_per_megagram + s.h2po4_with_hydroxyl_site_mol_p_per_megagram,
        .dissolved_aluminum_mol_per_m3 = -m.aluminum_phosphate_mol_per_m3,
        .dissolved_iron_mol_per_m3 = -m.iron_phosphate_mol_per_m3 - a.iron_hpo4_pairing_mol_p_per_m3 - a.iron_h2po4_pairing_mol_p_per_m3,
        .dissolved_calcium_mol_per_m3 = -m.dicalcium_phosphate_mol_per_m3 - m.monocalcium_phosphate_mol_per_m3 - 5 * m.hydroxyapatite_mol_per_m3 - a.calcium_po4_pairing_mol_p_per_m3 - a.calcium_hpo4_pairing_mol_p_per_m3 - a.calcium_h2po4_pairing_mol_p_per_m3,
        .dissolved_magnesium_mol_per_m3 = -a.magnesium_hpo4_pairing_mol_p_per_m3,
        .dissolved_hydrogen_mol_per_m3 = 2 * (m.aluminum_phosphate_mol_per_m3 + m.iron_phosphate_mol_per_m3) + 6 * m.hydroxyapatite_mol_per_m3 + m.dicalcium_phosphate_mol_per_m3 - (s.protonated_to_hydroxyl_site_mol_per_megagram + s.hydroxyl_to_deprotonated_site_mol_per_megagram) * site_to_water - a.po4_hydrogen_association_mol_p_per_m3 - a.hpo4_hydrogen_association_mol_p_per_m3 - a.h2po4_hydrogen_association_mol_p_per_m3,
        .dissolved_hydroxide_mol_per_m3 = -m.hydroxyapatite_mol_per_m3 + (s.h2po4_with_hydroxyl_site_mol_p_per_megagram + s.hpo4_with_hydroxyl_site_mol_p_per_megagram) * site_to_water,
        .water_mol_per_m3 = s.h2po4_with_protonated_site_mol_p_per_megagram * site_to_water,
        .aluminum_phosphate_solid_mol_per_m3 = m.aluminum_phosphate_mol_per_m3,
        .iron_phosphate_solid_mol_per_m3 = m.iron_phosphate_mol_per_m3,
        .dicalcium_phosphate_solid_mol_per_m3 = m.dicalcium_phosphate_mol_per_m3,
        .hydroxyapatite_solid_mol_per_m3 = m.hydroxyapatite_mol_per_m3,
        .monocalcium_phosphate_solid_mol_per_m3 = m.monocalcium_phosphate_mol_per_m3,
        .iron_hpo4_pair_mol_per_m3 = a.iron_hpo4_pairing_mol_p_per_m3,
        .iron_h2po4_pair_mol_per_m3 = a.iron_h2po4_pairing_mol_p_per_m3,
        .calcium_po4_pair_mol_per_m3 = a.calcium_po4_pairing_mol_p_per_m3,
        .calcium_hpo4_pair_mol_per_m3 = a.calcium_hpo4_pairing_mol_p_per_m3,
        .calcium_h2po4_pair_mol_per_m3 = a.calcium_h2po4_pairing_mol_p_per_m3,
        .magnesium_hpo4_pair_mol_per_m3 = a.magnesium_hpo4_pairing_mol_p_per_m3,
    };
}

fn validate(fluxes: Fluxes) !void {
    inline for (@typeInfo(MineralFluxes).@"struct".fields) |field| if (!std.math.isFinite(@field(fluxes.minerals, field.name))) return error.NonFinitePhosphateNetworkFlux;
    inline for (@typeInfo(exchange.Flux).@"struct".fields) |field| if (!std.math.isFinite(@field(fluxes.surface, field.name))) return error.NonFinitePhosphateNetworkFlux;
    inline for (@typeInfo(DissociationAndPairingFluxes).@"struct".fields) |field| if (!std.math.isFinite(@field(fluxes.aqueous, field.name))) return error.NonFinitePhosphateNetworkFlux;
    if (!std.math.isFinite(fluxes.soil_mass_per_water_volume_megagrams_per_m3) or fluxes.soil_mass_per_water_volume_megagrams_per_m3 <= 0) return error.InvalidPhosphateNetworkDensity;
}

fn validateState(state: State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const value = @field(state, field.name);
        if (!std.math.isFinite(value)) return error.NonFinitePhosphateNetworkState;
        if (value < -1e-12) return error.NegativePhosphateNetworkState;
    }
}

fn filledState(value: f64) State {
    var state: State = undefined;
    inline for (@typeInfo(State).@"struct".fields) |field| @field(state, field.name) = value;
    return state;
}

fn filledTransformations(value: f64) Transformations {
    var transformations: Transformations = undefined;
    inline for (@typeInfo(Transformations).@"struct".fields) |field| @field(transformations, field.name) = value;
    return transformations;
}

test "assembled phosphate network conserves phosphorus and exchange sites" {
    const density = 1.7;
    const result = try assemble(.{
        .soil_mass_per_water_volume_megagrams_per_m3 = density,
        .minerals = .{ .aluminum_phosphate_mol_per_m3 = 0.001, .iron_phosphate_mol_per_m3 = -0.002, .dicalcium_phosphate_mol_per_m3 = 0.003, .hydroxyapatite_mol_per_m3 = 0.004, .monocalcium_phosphate_mol_per_m3 = -0.001 },
        .surface = .{ .protonated_to_hydroxyl_site_mol_per_megagram = 0.002, .hydroxyl_to_deprotonated_site_mol_per_megagram = -0.001, .h2po4_with_protonated_site_mol_p_per_megagram = 0.003, .h2po4_with_hydroxyl_site_mol_p_per_megagram = -0.002, .hpo4_with_hydroxyl_site_mol_p_per_megagram = 0.001 },
        .aqueous = .{ .po4_hydrogen_association_mol_p_per_m3 = 0.003, .hpo4_hydrogen_association_mol_p_per_m3 = -0.002, .h2po4_hydrogen_association_mol_p_per_m3 = 0.001, .iron_hpo4_pairing_mol_p_per_m3 = 0.002, .iron_h2po4_pairing_mol_p_per_m3 = -0.001, .calcium_po4_pairing_mol_p_per_m3 = 0.002, .calcium_hpo4_pairing_mol_p_per_m3 = -0.003, .calcium_h2po4_pairing_mol_p_per_m3 = 0.001, .magnesium_hpo4_pairing_mol_p_per_m3 = 0.002 },
    });
    const dissolved_p = result.dissolved_po4_mol_p_per_m3 + result.dissolved_hpo4_mol_p_per_m3 + result.dissolved_h2po4_mol_p_per_m3 + result.dissolved_h3po4_mol_p_per_m3;
    const adsorbed_p = density * (result.adsorbed_hpo4_mol_p_per_megagram + result.adsorbed_h2po4_mol_p_per_megagram);
    const solid_p = result.aluminum_phosphate_solid_mol_per_m3 + result.iron_phosphate_solid_mol_per_m3 + result.dicalcium_phosphate_solid_mol_per_m3 + 3 * result.hydroxyapatite_solid_mol_per_m3 + 2 * result.monocalcium_phosphate_solid_mol_per_m3;
    const paired_p = result.iron_hpo4_pair_mol_per_m3 + result.iron_h2po4_pair_mol_per_m3 + result.calcium_po4_pair_mol_per_m3 + result.calcium_hpo4_pair_mol_per_m3 + result.calcium_h2po4_pair_mol_per_m3 + result.magnesium_hpo4_pair_mol_per_m3;
    try std.testing.expectApproxEqAbs(@as(f64, 0), dissolved_p + adsorbed_p + solid_p + paired_p, 1e-14);
    const sites = result.deprotonated_site_mol_per_megagram + result.hydroxyl_site_mol_per_megagram + result.protonated_site_mol_per_megagram + result.adsorbed_hpo4_mol_p_per_megagram + result.adsorbed_h2po4_mol_p_per_megagram;
    try std.testing.expectApproxEqAbs(@as(f64, 0), sites, 1e-14);
}

test "source-order phosphate assembly omits exchange density conversion" {
    const fluxes = Fluxes{
        .soil_mass_per_water_volume_megagrams_per_m3 = 3,
        .minerals = .{
            .aluminum_phosphate_mol_per_m3 = 0,
            .iron_phosphate_mol_per_m3 = 0,
            .dicalcium_phosphate_mol_per_m3 = 0,
            .hydroxyapatite_mol_per_m3 = 0,
            .monocalcium_phosphate_mol_per_m3 = 0,
        },
        .surface = .{
            .protonated_to_hydroxyl_site_mol_per_megagram = 0.2,
            .hydroxyl_to_deprotonated_site_mol_per_megagram = 0.1,
            .h2po4_with_protonated_site_mol_p_per_megagram = 0.4,
            .h2po4_with_hydroxyl_site_mol_p_per_megagram = 0.3,
            .hpo4_with_hydroxyl_site_mol_p_per_megagram = 0.5,
        },
        .aqueous = .{
            .po4_hydrogen_association_mol_p_per_m3 = 0,
            .hpo4_hydrogen_association_mol_p_per_m3 = 0,
            .h2po4_hydrogen_association_mol_p_per_m3 = 0,
            .iron_hpo4_pairing_mol_p_per_m3 = 0,
            .iron_h2po4_pairing_mol_p_per_m3 = 0,
            .calcium_po4_pairing_mol_p_per_m3 = 0,
            .calcium_hpo4_pairing_mol_p_per_m3 = 0,
            .calcium_h2po4_pairing_mol_p_per_m3 = 0,
            .magnesium_hpo4_pairing_mol_p_per_m3 = 0,
        },
    };
    const source = try assembleSourceOrder(fluxes);
    const corrected = try assemble(fluxes);

    try std.testing.expectEqual(@as(f64, -0.5), source.dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, -0.7), source.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectApproxEqAbs(@as(f64, -0.3), source.dissolved_hydrogen_mol_per_m3, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.8), source.dissolved_hydroxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0.4), source.water_mol_per_m3);
    try std.testing.expectEqual(3 * source.dissolved_hpo4_mol_p_per_m3, corrected.dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(3 * source.dissolved_h2po4_mol_p_per_m3, corrected.dissolved_h2po4_mol_p_per_m3);
}

test "phosphate network state commit is atomic" {
    var state = filledState(1);
    const before = state;
    var changes = filledTransformations(0);
    changes.dissolved_h2po4_mol_p_per_m3 = -2;
    changes.iron_phosphate_solid_mol_per_m3 = 0.2;
    try std.testing.expectError(error.NegativePhosphateNetworkState, commit(&state, changes));
    try std.testing.expectEqualDeep(before, state);

    changes.dissolved_h2po4_mol_p_per_m3 = -0.2;
    try commit(&state, changes);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.dissolved_h2po4_mol_p_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), state.iron_phosphate_solid_mol_per_m3, 1e-15);
}

test "source phosphate aqueous update floors before ceiling only" {
    const current = filledState(1);
    var changes = filledTransformations(0);
    changes.dissolved_po4_mol_p_per_m3 = -2;
    changes.dissolved_hpo4_mol_p_per_m3 = 0.25;
    changes.magnesium_hpo4_pair_mol_per_m3 = -0.5;

    const terminal = try applySourceOrderAqueousUpdate(
        current,
        changes,
        0.01,
        .iteration_ceiling,
    );
    try std.testing.expectEqualDeep(current, terminal);

    const continuing = try applySourceOrderAqueousUpdate(
        current,
        changes,
        0.01,
        .before_iteration_ceiling,
    );
    try std.testing.expectEqual(
        @as(f64, 0.01),
        continuing.dissolved_po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 1.25),
        continuing.dissolved_hpo4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        continuing.magnesium_hpo4_pair_mol_per_m3,
    );
}

test "source phosphate surface update is unfloored and ceiling gated" {
    const current = filledState(1);
    var changes = filledTransformations(0);
    changes.deprotonated_site_mol_per_megagram = -1;
    changes.hydroxyl_site_mol_per_megagram = 0.25;
    changes.adsorbed_h2po4_mol_p_per_megagram = -0.5;

    const terminal = try applySourceOrderSurfaceUpdate(
        current,
        changes,
        .iteration_ceiling,
    );
    try std.testing.expectEqualDeep(current, terminal);

    const continuing = try applySourceOrderSurfaceUpdate(
        current,
        changes,
        .before_iteration_ceiling,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        continuing.deprotonated_site_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @as(f64, 1.25),
        continuing.hydroxyl_site_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        continuing.adsorbed_h2po4_mol_p_per_megagram,
    );
}

test "source phosphate mineral update is unfloored and ceiling gated" {
    const current = filledState(1);
    var changes = filledTransformations(0);
    changes.aluminum_phosphate_solid_mol_per_m3 = -1;
    changes.hydroxyapatite_solid_mol_per_m3 = 0.25;
    changes.monocalcium_phosphate_solid_mol_per_m3 = -0.5;

    const terminal = try applySourceOrderMineralUpdate(
        current,
        changes,
        .iteration_ceiling,
    );
    try std.testing.expectEqualDeep(current, terminal);

    const continuing = try applySourceOrderMineralUpdate(
        current,
        changes,
        .before_iteration_ceiling,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        continuing.aluminum_phosphate_solid_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 1.25),
        continuing.hydroxyapatite_solid_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        continuing.monocalcium_phosphate_solid_mol_per_m3,
    );
}

test "surface extents use runtime soil mass per water volume" {
    const density = 2.0;
    const result = try assemble(.{
        .soil_mass_per_water_volume_megagrams_per_m3 = density,
        .minerals = std.mem.zeroes(MineralFluxes),
        .surface = .{
            .protonated_to_hydroxyl_site_mol_per_megagram = 0.2,
            .hydroxyl_to_deprotonated_site_mol_per_megagram = 0.1,
            .h2po4_with_protonated_site_mol_p_per_megagram = 0.3,
            .h2po4_with_hydroxyl_site_mol_p_per_megagram = 0,
            .hpo4_with_hydroxyl_site_mol_p_per_megagram = 0,
        },
        .aqueous = std.mem.zeroes(DissociationAndPairingFluxes),
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.6),
        result.dissolved_hydrogen_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.6),
        result.water_mol_per_m3,
        1e-15,
    );
}
