const std = @import("std");

pub const Source = enum(u8) {
    precipitation = 1,
    irrigation = 2,
    soil = 3,
};

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const Inputs = struct {
    dissolved_calcium_mol_per_source_volume: f64, // CCA1
    dissolved_carbonate_mol_per_m3: f64, // CCO31
    calcium_activity_mol_per_source_volume: f64, // ACA1
    carbonate_activity_mol_per_m3: f64, // ACO31
    divalent_activity_coefficient: f64, // A2S
    calcite_solid_mol_per_megagram: f64, // PCACO1
    calcite_solubility_product: f64, // SPCAC
    substrate_limit_fraction: f64, // FION
    maximum_reaction_mol_per_iteration: f64, // TPD
};

pub const Result = struct {
    substrate_limit_mol_per_iteration: f64, // XMINP
    equilibrium_calcium_activity: f64, // ACA1Q
    bounded_unsuppressed_extent_mol_per_iteration: f64,
    calcite_extent_mol_per_iteration: f64, // RPCACX
};

/// Direct translation of STARTE.F lines 531--535 within the outer `K=3`
/// branch. STARTE evaluates the dynamic calcite limits, then suppresses the
/// initialization extent with a literal `0.0` multiplier. Production calcite
/// chemistry is owned separately by the enhanced mineral solver.
pub fn evaluate(source: Source, mode: SaltEquilibriumMode, inputs: Inputs) !?Result {
    if (source != .soil) return null;
    if (mode == .static) return .{
        .substrate_limit_mol_per_iteration = 0,
        .equilibrium_calcium_activity = 0,
        .bounded_unsuppressed_extent_mol_per_iteration = 0,
        .calcite_extent_mol_per_iteration = 0,
    };

    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteInitialCalciteInput;
    }
    if (inputs.dissolved_calcium_mol_per_source_volume < 0 or
        inputs.dissolved_carbonate_mol_per_m3 < 0 or
        inputs.calcium_activity_mol_per_source_volume < 0 or
        inputs.carbonate_activity_mol_per_m3 <= 0 or
        inputs.divalent_activity_coefficient <= 0 or
        inputs.calcite_solubility_product <= 0 or
        inputs.substrate_limit_fraction < 0 or
        inputs.maximum_reaction_mol_per_iteration < 0)
        return error.InvalidInitialCalciteInput;

    const substrate_limit = inputs.substrate_limit_fraction *
        @min(
            inputs.dissolved_calcium_mol_per_source_volume,
            inputs.dissolved_carbonate_mol_per_m3,
        );
    const equilibrium_calcium_activity = inputs.calcite_solubility_product /
        inputs.carbonate_activity_mol_per_m3;
    const bounded = @max(
        -@max(0.0, inputs.calcite_solid_mol_per_megagram),
        -inputs.maximum_reaction_mol_per_iteration,
        @min(
            inputs.maximum_reaction_mol_per_iteration,
            substrate_limit,
            (inputs.calcium_activity_mol_per_source_volume -
                equilibrium_calcium_activity) /
                inputs.divalent_activity_coefficient,
        ),
    );
    const result: Result = .{
        .substrate_limit_mol_per_iteration = substrate_limit,
        .equilibrium_calcium_activity = equilibrium_calcium_activity,
        .bounded_unsuppressed_extent_mol_per_iteration = bounded,
        .calcite_extent_mol_per_iteration = 0.0 * bounded,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteInitialCalciteResult;
    }
    return result;
}

fn dynamicInputs() Inputs {
    return .{
        .dissolved_calcium_mol_per_source_volume = 8,
        .dissolved_carbonate_mol_per_m3 = 6,
        .calcium_activity_mol_per_source_volume = 5,
        .carbonate_activity_mol_per_m3 = 2,
        .divalent_activity_coefficient = 0.5,
        .calcite_solid_mol_per_megagram = 4,
        .calcite_solubility_product = 2,
        .substrate_limit_fraction = 0.25,
        .maximum_reaction_mol_per_iteration = 3,
    };
}

test "STARTE dynamic calcite evaluates limits then applies literal suppression" {
    const result = (try evaluate(.soil, .dynamic, dynamicInputs())).?;
    try std.testing.expectEqual(@as(f64, 1.5), result.substrate_limit_mol_per_iteration);
    try std.testing.expectEqual(@as(f64, 1), result.equilibrium_calcium_activity);
    try std.testing.expectEqual(@as(f64, 1.5), result.bounded_unsuppressed_extent_mol_per_iteration);
    try std.testing.expectEqual(@as(f64, 0), result.calcite_extent_mol_per_iteration);
}

test "STARTE static calcite reset does not inspect dormant dynamic inputs" {
    var invalid = dynamicInputs();
    invalid.carbonate_activity_mol_per_m3 = std.math.nan(f64);
    const result = (try evaluate(.soil, .static, invalid)).?;
    try std.testing.expectEqual(@as(f64, 0), result.calcite_extent_mol_per_iteration);
}

test "STARTE non-soil calcite guard is dormant and dynamic validation is atomic" {
    var invalid = dynamicInputs();
    invalid.divalent_activity_coefficient = 0;
    try std.testing.expect((try evaluate(.irrigation, .dynamic, invalid)) == null);
    try std.testing.expectError(
        error.InvalidInitialCalciteInput,
        evaluate(.soil, .dynamic, invalid),
    );
}
