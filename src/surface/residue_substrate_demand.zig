//! **A5 DISPOSITION: never production bound.**
//!
//! Legacy source: `ecosys_f77/hour1.f` lines 4658--4678. This module is an exact,
//! tested, source-order translation of that range and is imported only by
//! `src/module_index.zig`, so it is reachable by `zig build test` and unreachable
//! from `executeHourlyScience`. That is intentional and must stay that way.
//!
//! Classification: diagnostic-only. This kernel resets legacy running totals that no
//! production module accumulates and no production module reads. Production
//! reconstructs the equivalent totals on demand in
//! `landscape_mass_balance_runtime.reconstruct`, which cannot drift from the
//! state it summarizes. Binding a reset for an accumulator that nothing
//! accumulates would add cost and no behaviour.
//!
//! Superseded by: the litter-side uptake owners, which carry their own previous-step demand.
//!
//! Field census: 11 fields, 2 shared, and both shared names
//! (`current_organic`, `previous_organic`) belong to
//! `substrate_uptake_state_rollover`, the other class D module. Nothing in
//! production reads a stored litter substrate demand.
//!
//! Do not bind this module, and do not delete it: the tests are a source-order
//! comparison oracle for the range above. Full argument and the census behind
//! it: `docs/traceability/hour1_2039_5200_binding_survey.md`.
const std = @import("std");

pub const MineralAndGasDemands = struct {
    oxygen_g_o_per_timestep: f64,
    ammonium_g_n_per_timestep: f64,
    nitrate_g_n_per_timestep: f64,
    nitrite_g_n_per_timestep: f64,
    nitrous_oxide_g_n_per_timestep: f64,
    hydrogen_phosphate_g_p_per_timestep: f64,
    dihydrogen_phosphate_g_p_per_timestep: f64,
};

pub const OrganicDemands = struct {
    dissolved_organic_carbon_g_c_per_timestep: []f64,
    acetate_g_c_per_timestep: []f64,
};

pub const DemandHistory = struct {
    previous: *MineralAndGasDemands,
    current: *MineralAndGasDemands,
    previous_organic: OrganicDemands,
    current_organic: OrganicDemands,
};

pub const RolloverError = error{
    OrganicClassCountMismatch,
    NonFiniteCurrentDemand,
};

/// Translates `hour1.f` lines 4658--4678 for one surface-residue layer.
///
/// Organic class count is runtime-defined; all four organic-demand slices must
/// have identical lengths. Previous values are overwritten and current
/// accumulators are reset only after every input has passed validation.
pub fn rollover(history: DemandHistory) RolloverError!void {
    const class_count = history.current_organic.dissolved_organic_carbon_g_c_per_timestep.len;
    if (history.current_organic.acetate_g_c_per_timestep.len != class_count or
        history.previous_organic.dissolved_organic_carbon_g_c_per_timestep.len != class_count or
        history.previous_organic.acetate_g_c_per_timestep.len != class_count)
    {
        return error.OrganicClassCountMismatch;
    }
    inline for (std.meta.fields(MineralAndGasDemands)) |field| {
        if (!std.math.isFinite(@field(history.current, field.name))) {
            return error.NonFiniteCurrentDemand;
        }
    }
    for (history.current_organic.dissolved_organic_carbon_g_c_per_timestep) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteCurrentDemand;
    }
    for (history.current_organic.acetate_g_c_per_timestep) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteCurrentDemand;
    }

    history.previous.oxygen_g_o_per_timestep = history.current.oxygen_g_o_per_timestep;
    history.previous.ammonium_g_n_per_timestep = history.current.ammonium_g_n_per_timestep;
    history.previous.nitrate_g_n_per_timestep = history.current.nitrate_g_n_per_timestep;
    history.previous.nitrite_g_n_per_timestep = history.current.nitrite_g_n_per_timestep;
    history.previous.nitrous_oxide_g_n_per_timestep = history.current.nitrous_oxide_g_n_per_timestep;
    history.previous.hydrogen_phosphate_g_p_per_timestep =
        history.current.hydrogen_phosphate_g_p_per_timestep;
    history.previous.dihydrogen_phosphate_g_p_per_timestep =
        history.current.dihydrogen_phosphate_g_p_per_timestep;

    history.current.oxygen_g_o_per_timestep = 0.0;
    history.current.ammonium_g_n_per_timestep = 0.0;
    history.current.nitrate_g_n_per_timestep = 0.0;
    history.current.nitrite_g_n_per_timestep = 0.0;
    history.current.nitrous_oxide_g_n_per_timestep = 0.0;
    history.current.hydrogen_phosphate_g_p_per_timestep = 0.0;
    history.current.dihydrogen_phosphate_g_p_per_timestep = 0.0;

    for (0..class_count) |class_index| {
        history.previous_organic.dissolved_organic_carbon_g_c_per_timestep[class_index] =
            history.current_organic.dissolved_organic_carbon_g_c_per_timestep[class_index];
        history.previous_organic.acetate_g_c_per_timestep[class_index] =
            history.current_organic.acetate_g_c_per_timestep[class_index];
        history.current_organic.dissolved_organic_carbon_g_c_per_timestep[class_index] = 0.0;
        history.current_organic.acetate_g_c_per_timestep[class_index] = 0.0;
    }
}

test "rollover copies demands then clears runtime-sized accumulators" {
    var current = MineralAndGasDemands{
        .oxygen_g_o_per_timestep = 1.0,
        .ammonium_g_n_per_timestep = 2.0,
        .nitrate_g_n_per_timestep = 3.0,
        .nitrite_g_n_per_timestep = 4.0,
        .nitrous_oxide_g_n_per_timestep = 5.0,
        .hydrogen_phosphate_g_p_per_timestep = 6.0,
        .dihydrogen_phosphate_g_p_per_timestep = 7.0,
    };
    var previous = std.mem.zeroes(MineralAndGasDemands);
    var current_doc = [_]f64{ 8.0, 9.0, 10.0 };
    var current_acetate = [_]f64{ 11.0, 12.0, 13.0 };
    var previous_doc = [_]f64{ 0.0, 0.0, 0.0 };
    var previous_acetate = [_]f64{ 0.0, 0.0, 0.0 };

    try rollover(.{
        .previous = &previous,
        .current = &current,
        .previous_organic = .{
            .dissolved_organic_carbon_g_c_per_timestep = &previous_doc,
            .acetate_g_c_per_timestep = &previous_acetate,
        },
        .current_organic = .{
            .dissolved_organic_carbon_g_c_per_timestep = &current_doc,
            .acetate_g_c_per_timestep = &current_acetate,
        },
    });

    try std.testing.expectEqual(@as(f64, 1.0), previous.oxygen_g_o_per_timestep);
    try std.testing.expectEqual(@as(f64, 7.0), previous.dihydrogen_phosphate_g_p_per_timestep);
    try std.testing.expectEqual(std.mem.zeroes(MineralAndGasDemands), current);
    try std.testing.expectEqualSlices(f64, &.{ 8.0, 9.0, 10.0 }, &previous_doc);
    try std.testing.expectEqualSlices(f64, &.{ 11.0, 12.0, 13.0 }, &previous_acetate);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0 }, &current_doc);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0 }, &current_acetate);
}

test "invalid demand fails before any rollover mutation" {
    var current = std.mem.zeroes(MineralAndGasDemands);
    current.nitrate_g_n_per_timestep = std.math.nan(f64);
    var previous = std.mem.zeroes(MineralAndGasDemands);
    previous.oxygen_g_o_per_timestep = 42.0;
    var current_doc = [_]f64{1.0};
    var current_acetate = [_]f64{2.0};
    var previous_doc = [_]f64{3.0};
    var previous_acetate = [_]f64{4.0};

    try std.testing.expectError(error.NonFiniteCurrentDemand, rollover(.{
        .previous = &previous,
        .current = &current,
        .previous_organic = .{
            .dissolved_organic_carbon_g_c_per_timestep = &previous_doc,
            .acetate_g_c_per_timestep = &previous_acetate,
        },
        .current_organic = .{
            .dissolved_organic_carbon_g_c_per_timestep = &current_doc,
            .acetate_g_c_per_timestep = &current_acetate,
        },
    }));
    try std.testing.expectEqual(@as(f64, 42.0), previous.oxygen_g_o_per_timestep);
    try std.testing.expectEqualSlices(f64, &.{3.0}, &previous_doc);
}
