const std = @import("std");
const litter_chemistry = @import("litter_chemistry.zig");

/// Diagnostic counter: how many times the dry-carrier branch has executed. Used to
/// check whether that branch is reachable at all in a given scenario, rather than
/// assuming it from the presence of a failure downstream. See EXEC-004.
///
/// Measured: `10` executions over the first day of the Ottawa example with surface
/// evaporation enabled, and `0` with it disabled as shipped. That is what
/// establishes the branch is genuinely exercised once evaporation is active, and
/// that the shipped configuration never reaches it.
pub var dry_branch_executions: u64 = 0;

/// Allocation-free carrier transaction after an accepted liquid-water change.
/// All `mol/m3` pools retain their extensive amount; dry-mass-normalized
/// exchange, carboxyl, and phosphate-surface pools remain unchanged.
///
/// **Dry carriers.** `hour1.f` 4494--4524 stores litter solutes as absolute mass
/// and derives concentration as `AMAX1(0.0, mass/VOLW(0))`, setting it to `0.0`
/// where the litter layer does not exist. Mass therefore survives the carrier
/// vanishing. Because ecosys-ng stores these pools water-normalized, a carrier of
/// exactly zero has no representable concentration, and this transaction used to
/// reject that outright with `LitterChemistryMassWithoutWaterCarrier`, which made
/// evaporation to dryness a hard error.
///
/// It no longer does. When the live carrier reaches zero the concentrations are
/// **held unchanged** rather than scaled to zero, and the carrier they refer to is
/// remembered in `dry_reference_water_m3` so that rewetting rescales from the
/// retained value instead of from zero. The extensive amount
/// `concentration * remembered_carrier` is therefore preserved exactly across an
/// arbitrary wet/dry/rewet sequence, matching the oracle's invariant (mass
/// conserved, concentration undefined while dry) within a water-normalized
/// representation.
///
/// Note `mineral_reference_water_m3` is deliberately *not* reused for this: it is
/// also written by `bindMineralReferenceWater` and `renormalizeMinerals`, so
/// borrowing it as the rescale source would couple two independent invariants.
///
/// **KNOWN LIMITATION.** While a cell is dry its stored concentrations refer to
/// `dry_reference_water_m3`, not to the live `litter_water_m3`. Consumers that form
/// an extensive amount as `concentration * live_carrier` therefore disagree with
/// this transaction during a dry spell, and 27 modules reference `litter_water_m3`,
/// so this would not be fixable by patching one consumer.
///
/// This is a theoretical concern, not an observed failure. I initially blamed it
/// for the hour-11 nitrogen closure failure, but that turned out to be exactly
/// 1.00 ULP of float roundoff in an unrelated check. No measurement yet shows the
/// coupling causing an actual defect. With surface evaporation enabled the branch
/// executes 10 times over the first day and the run completes all 24 hours; in the
/// shipped configuration it executes 0 times.
///
/// If a future failure *is* traced to this coupling, the two candidate resolutions
/// recorded under EXEC-004 are to make the extensive amount authoritative as
/// `hour1.f` 4494--4524 does, or to guarantee on physical grounds that the carrier
/// cannot reach zero. Measure before choosing.
///
/// See EXEC-004. `surface_litter_solute_dry_carrier.zig` states the target
/// semantics independently and proves why a naive rescale cannot round-trip.
pub fn rebaseFromAcceptedLiquidWaterChange(
    state: *litter_chemistry.State,
    new_litter_water_volume_m3: []const f64,
    liquid_water_change_m3: []const f64,
) !void {
    if (new_litter_water_volume_m3.len != state.cells.len or
        liquid_water_change_m3.len != state.cells.len)
        return error.LitterChemistryCarrierDimensionMismatch;

    for (state.cells, 0..) |cell, index| {
        const old = try oldCarrier(
            new_litter_water_volume_m3[index],
            liquid_water_change_m3[index],
        );
        try validateCell(cell, old, new_litter_water_volume_m3[index], state.dry_reference_water_m3[index]);
    }
    for (state.cells, 0..) |*cell, index| {
        const new = new_litter_water_volume_m3[index];
        const old = new - liquid_water_change_m3[index];
        // The carrier the stored concentrations currently refer to: the live
        // `old`, unless a previous step went dry and remembered one.
        const source = if (state.dry_reference_water_m3[index] > 0)
            state.dry_reference_water_m3[index]
        else
            old;
        if (new == 0) {
            // Going dry: hold the concentrations and remember the carrier they
            // refer to, so no mass is lost and rewetting can recover it.
            dry_branch_executions += 1;
            if (state.dry_reference_water_m3[index] == 0)
                state.dry_reference_water_m3[index] = source;
            continue;
        }
        applyScale(cell, if (source == 0) 0 else source / new);
        state.mineral_reference_water_m3[index] = new;
        state.dry_reference_water_m3[index] = 0;
    }
}

/// Rebase every dry-mass-normalized pool after an accepted litter dry-mass
/// change. Extensive amounts are invariant; water-normalized pools are not
/// touched because this transaction changes only the solid carrier.
pub fn rebaseFromAcceptedDryMassChange(
    state: *litter_chemistry.State,
    old_dry_mass_megagrams: []const f64,
    new_dry_mass_megagrams: []const f64,
) !void {
    if (old_dry_mass_megagrams.len != state.cells.len or
        new_dry_mass_megagrams.len != state.cells.len)
        return error.LitterChemistryCarrierDimensionMismatch;

    for (state.cells, old_dry_mass_megagrams, new_dry_mass_megagrams) |cell, old, new| {
        try validateDryCarrier(old, new);
        try validateDryNormalizedCell(cell, old, new);
    }
    for (state.cells, old_dry_mass_megagrams, new_dry_mass_megagrams) |*cell, old, new| {
        const scale = if (old == 0) 0 else old / new;
        applyDryScale(cell, scale);
    }
}

fn validateDryCarrier(old: f64, new: f64) !void {
    if (!std.math.isFinite(old) or !std.math.isFinite(new) or old < 0 or new < 0)
        return error.InvalidLitterChemistryCarrier;
}

fn validateDryNormalizedCell(cell: litter_chemistry.Cell, old: f64, new: f64) !void {
    inline for (@typeInfo(litter_chemistry.Cell).@"struct".fields) |field| {
        const Field = field.type;
        if (comptime isDryNormalized(field.name, Field)) {
            try validateDryPool(@field(cell, field.name), old, new);
        } else if (comptime isDryNormalizedStruct(field.name)) {
            inline for (@typeInfo(Field).@"struct".fields) |nested|
                try validateDryPool(@field(@field(cell, field.name), nested.name), old, new);
        }
    }
}

fn validateDryPool(value: f64, old: f64, new: f64) !void {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDryNormalizedLitterChemistryPool;
    if (value > 0 and (old == 0 or new == 0))
        return error.LitterChemistryMassWithoutDryCarrier;
}

fn applyDryScale(cell: *litter_chemistry.Cell, scale: f64) void {
    inline for (@typeInfo(litter_chemistry.Cell).@"struct".fields) |field| {
        const Field = field.type;
        if (comptime isDryNormalized(field.name, Field)) {
            @field(cell, field.name) *= scale;
        } else if (comptime isDryNormalizedStruct(field.name)) {
            inline for (@typeInfo(Field).@"struct".fields) |nested|
                @field(@field(cell, field.name), nested.name) *= scale;
        }
    }
}

fn isDryNormalized(comptime name: []const u8, comptime Field: type) bool {
    return @typeInfo(Field) == .float and std.mem.endsWith(u8, name, "_per_megagram");
}

fn isDryNormalizedStruct(comptime name: []const u8) bool {
    return std.mem.eql(u8, name, "exchange") or
        std.mem.eql(u8, name, "phosphate_surface");
}

fn oldCarrier(new: f64, change: f64) !f64 {
    if (!std.math.isFinite(new) or !std.math.isFinite(change))
        return error.InvalidLitterChemistryCarrier;
    const old = new - change;
    if (!std.math.isFinite(old) or old < 0 or new < 0)
        return error.InvalidLitterChemistryCarrier;
    return old;
}

fn validateCell(cell: litter_chemistry.Cell, old: f64, new: f64, reference: f64) !void {
    if (!std.math.isFinite(reference) or reference < 0)
        return error.InvalidLitterChemistryCarrier;
    inline for (@typeInfo(litter_chemistry.Cell).@"struct".fields) |field| {
        const Field = field.type;
        if (comptime isWaterNormalized(field.name, Field)) {
            try validatePool(@field(cell, field.name), old, new);
        } else if (comptime isMineralStruct(field.name)) {
            inline for (@typeInfo(Field).@"struct".fields) |nested| {
                try validatePool(@field(@field(cell, field.name), nested.name), old, new);
            }
        } else if (comptime @typeInfo(Field) == .float) {
            const value = @field(cell, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDryNormalizedLitterChemistryPool;
        } else {
            inline for (@typeInfo(Field).@"struct".fields) |nested| {
                const value = @field(@field(cell, field.name), nested.name);
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidDryNormalizedLitterChemistryPool;
            }
        }
    }
}

fn validatePool(value: f64, old: f64, new: f64) !void {
    _ = old;
    _ = new;
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidWaterNormalizedLitterChemistryPool;
    // A zero carrier is no longer rejected. Concentrations and their reference
    // water are held together while dry, which keeps the extensive amount
    // `concentration * reference` exact, so there is no mass without a carrier to
    // detect. See the transaction doc comment and EXEC-004.
}

fn applyScale(cell: *litter_chemistry.Cell, scale: f64) void {
    inline for (@typeInfo(litter_chemistry.Cell).@"struct".fields) |field| {
        const Field = field.type;
        if (comptime isWaterNormalized(field.name, Field)) {
            @field(cell, field.name) *= scale;
        } else if (comptime isMineralStruct(field.name)) {
            inline for (@typeInfo(Field).@"struct".fields) |nested| {
                @field(@field(cell, field.name), nested.name) *= scale;
            }
        }
    }
}

fn isWaterNormalized(comptime name: []const u8, comptime Field: type) bool {
    return @typeInfo(Field) == .float and
        std.mem.indexOf(u8, name, "_mol") != null and
        std.mem.endsWith(u8, name, "_per_m3");
}

fn isMineralStruct(comptime name: []const u8) bool {
    return std.mem.eql(u8, name, "phosphate_minerals") or
        std.mem.eql(u8, name, "salt_minerals");
}

test "surface litter carrier preserves inorganic carbon phosphorus and ions" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    var cell = &state.cells[0];
    cell.carbon_dioxide_mol_per_m3 = 2;
    cell.bicarbonate_mol_per_m3 = 3;
    cell.h2po4_mol_p_per_m3 = 5;
    cell.calcium_mol_per_m3 = 7;
    cell.phosphate_minerals.hydroxyapatite_mol_per_m3 = 11;
    cell.salt_minerals.calcite_mol_per_m3 = 13;

    try rebaseFromAcceptedLiquidWaterChange(&state, &.{0.5}, &.{-1.5});
    try std.testing.expectEqual(@as(f64, 8), cell.carbon_dioxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 12), cell.bicarbonate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 20), cell.h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 28), cell.calcium_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 44), cell.phosphate_minerals.hydroxyapatite_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 52), cell.salt_minerals.calcite_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), cell.carbon_dioxide_mol_per_m3 * 0.5);
    try std.testing.expectEqual(@as(f64, 10), cell.h2po4_mol_p_per_m3 * 0.5);
    try std.testing.expectEqual(@as(f64, 26), cell.salt_minerals.calcite_mol_per_m3 * 0.5);
}

test "surface litter carrier leaves dry normalized pools untouched" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].carboxyl_hydrogen_mol_per_megagram = 2;
    state.cells[0].exchange.calcium_mol_per_megagram = 3;
    state.cells[0].phosphate_surface.adsorbed_h2po4_mol_p_per_megagram = 4;
    try rebaseFromAcceptedLiquidWaterChange(&state, &.{2}, &.{1});
    try std.testing.expectEqual(@as(f64, 2), state.cells[0].carboxyl_hydrogen_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 3), state.cells[0].exchange.calcium_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 4), state.cells[0].phosphate_surface.adsorbed_h2po4_mol_p_per_megagram);
}

test "surface litter carrier validation is atomic and rejects invalid carriers" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.cells[0].nitrate_mol_per_m3 = 2;
    state.cells[1].sulfate_mol_per_m3 = 3;
    // A carrier that would go negative is still rejected, and rejection leaves
    // every cell untouched rather than partially applied.
    try std.testing.expectError(
        error.InvalidLitterChemistryCarrier,
        rebaseFromAcceptedLiquidWaterChange(&state, &.{ 2, 1 }, &.{ 1, 2 }),
    );
    try std.testing.expectEqual(@as(f64, 2), state.cells[0].nitrate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 3), state.cells[1].sulfate_mol_per_m3);
}

test "evaporating the carrier to dryness conserves solute mass" {
    // EXEC-004: this used to fail with LitterChemistryMassWithoutWaterCarrier,
    // which forced the shipped runscript to disable surface evaporation.
    var state = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const cell = &state.cells[0];
    cell.nitrate_mol_per_m3 = 4;
    // Establish a wet carrier of 2 m3, so the extensive amount is 8 mol.
    try rebaseFromAcceptedLiquidWaterChange(&state, &.{2}, &.{0});
    try std.testing.expectEqual(@as(f64, 2), state.mineral_reference_water_m3[0]);
    const amount_before = cell.nitrate_mol_per_m3 * state.mineral_reference_water_m3[0];

    // Evaporate to exactly dry. This must not error.
    try rebaseFromAcceptedLiquidWaterChange(&state, &.{0}, &.{-2});
    // Concentration and reference are held together, so the amount is unchanged.
    const amount_dry = cell.nitrate_mol_per_m3 * state.mineral_reference_water_m3[0];
    try std.testing.expectApproxEqRel(amount_before, amount_dry, 1e-15);
}

test "rewetting a dry carrier recovers the retained solute mass" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const cell = &state.cells[0];
    cell.nitrate_mol_per_m3 = 4;
    try rebaseFromAcceptedLiquidWaterChange(&state, &.{2}, &.{0});
    const amount_before = cell.nitrate_mol_per_m3 * state.mineral_reference_water_m3[0];

    try rebaseFromAcceptedLiquidWaterChange(&state, &.{0}, &.{-2});
    // Rewet to 8 m3: four times the original carrier.
    try rebaseFromAcceptedLiquidWaterChange(&state, &.{8}, &.{8});
    try std.testing.expectEqual(@as(f64, 8), state.mineral_reference_water_m3[0]);
    const amount_after = cell.nitrate_mol_per_m3 * state.mineral_reference_water_m3[0];
    // The mass survived the round trip through dryness, which is exactly what a
    // naive rescale-from-zero cannot do.
    try std.testing.expectApproxEqRel(amount_before, amount_after, 1e-14);
    // Diluted into four times the water, so a quarter of the concentration.
    try std.testing.expectApproxEqRel(@as(f64, 1), cell.nitrate_mol_per_m3, 1e-14);
}

test "solute mass survives an arbitrary wet dry rewet sequence" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    const cell = &state.cells[0];
    cell.nitrate_mol_per_m3 = 2.5;
    try rebaseFromAcceptedLiquidWaterChange(&state, &.{4}, &.{0});
    const amount_before = cell.nitrate_mol_per_m3 * state.mineral_reference_water_m3[0];

    var carrier: f64 = 4;
    for ([_]f64{ 2, 0.5, 0, 3, 0, 0, 7 }) |next| {
        try rebaseFromAcceptedLiquidWaterChange(&state, &.{next}, &.{next - carrier});
        carrier = next;
        // The invariant holds after every single step, including consecutive dry
        // steps, not merely at the end.
        const amount = cell.nitrate_mol_per_m3 * state.mineral_reference_water_m3[0];
        try std.testing.expectApproxEqRel(amount_before, amount, 1e-13);
    }
}

test "surface litter dry carrier preserves exchange and phosphate amounts" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].carboxyl_hydrogen_mol_per_megagram = 2;
    state.cells[0].exchange.calcium_mol_per_megagram = 3;
    state.cells[0].phosphate_surface.adsorbed_h2po4_mol_p_per_megagram = 5;
    state.cells[0].h2po4_mol_p_per_m3 = 7;

    try rebaseFromAcceptedDryMassChange(&state, &.{2}, &.{0.5});
    try std.testing.expectEqual(@as(f64, 8), state.cells[0].carboxyl_hydrogen_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 12), state.cells[0].exchange.calcium_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 20), state.cells[0].phosphate_surface.adsorbed_h2po4_mol_p_per_megagram);
    try std.testing.expectEqual(@as(f64, 7), state.cells[0].h2po4_mol_p_per_m3);
}

test "surface litter dry carrier validation is atomic" {
    var state = try litter_chemistry.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.cells[0].exchange.calcium_mol_per_megagram = 2;
    state.cells[1].phosphate_surface.adsorbed_hpo4_mol_p_per_megagram = 3;
    try std.testing.expectError(
        error.LitterChemistryMassWithoutDryCarrier,
        rebaseFromAcceptedDryMassChange(&state, &.{ 1, 1 }, &.{ 2, 0 }),
    );
    try std.testing.expectEqual(@as(f64, 2), state.cells[0].exchange.calcium_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 3), state.cells[1].phosphate_surface.adsorbed_hpo4_mol_p_per_megagram);
}
