const std = @import("std");

pub const Inputs = struct {
    total_carboxyl_sites_mol_per_megagram: f64,
    protonated_carboxyl_sites_mol_per_megagram: f64,
    hydrogen_activity_mol_per_m3: f64,
};

pub const Parameters = struct {
    carboxyl_dissociation_mol_per_m3: f64,
    maximum_exchange_mol_per_megagram_iteration: f64,
    negligible_open_sites_mol_per_megagram: f64,
};

pub const Result = struct {
    open_carboxyl_sites_mol_per_megagram: f64,
    equilibrium_open_carboxyl_sites_mol_per_megagram: f64,
    protonation_change_mol_per_megagram_iteration: f64,
    protonated_carboxyl_sites_mol_per_megagram: f64,
};

/// Direct single-iteration translation of `starte.f` lines 1896--1900. This
/// initialization map intentionally differs from the later SOLUTE carboxyl
/// exchange rate and is retained for source-order compatibility validation.
pub fn iterate(inputs: Inputs, parameters: Parameters) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceCarboxylInitializationInput;
        if (value < 0) return error.InvalidSurfaceCarboxylInitializationInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceCarboxylInitializationParameter;
        if (value < 0)
            return error.InvalidSurfaceCarboxylInitializationParameter;
    }
    if (inputs.hydrogen_activity_mol_per_m3 <= 0 or
        parameters.carboxyl_dissociation_mol_per_m3 <= 0 or
        inputs.protonated_carboxyl_sites_mol_per_megagram >
            inputs.total_carboxyl_sites_mol_per_megagram)
        return error.InvalidSurfaceCarboxylInitializationInput;

    const open_sites = @max(
        parameters.negligible_open_sites_mol_per_megagram,
        inputs.total_carboxyl_sites_mol_per_megagram -
            inputs.protonated_carboxyl_sites_mol_per_megagram,
    );
    const equilibrium_open_sites = @min(
        inputs.total_carboxyl_sites_mol_per_megagram,
        parameters.carboxyl_dissociation_mol_per_m3 *
            inputs.protonated_carboxyl_sites_mol_per_megagram /
            inputs.hydrogen_activity_mol_per_m3,
    );
    const change = @max(
        -parameters.maximum_exchange_mol_per_megagram_iteration *
            inputs.protonated_carboxyl_sites_mol_per_megagram,
        -parameters.maximum_exchange_mol_per_megagram_iteration,
        @min(
            parameters.maximum_exchange_mol_per_megagram_iteration,
            open_sites - equilibrium_open_sites,
        ),
    );
    const protonated_sites =
        inputs.protonated_carboxyl_sites_mol_per_megagram + change;
    const result: Result = .{
        .open_carboxyl_sites_mol_per_megagram = open_sites,
        .equilibrium_open_carboxyl_sites_mol_per_megagram = equilibrium_open_sites,
        .protonation_change_mol_per_megagram_iteration = change,
        .protonated_carboxyl_sites_mol_per_megagram = protonated_sites,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceCarboxylInitializationResult;
    if (protonated_sites < 0 or
        protonated_sites > inputs.total_carboxyl_sites_mol_per_megagram)
        return error.InvalidSurfaceCarboxylInitializationResult;
    return result;
}

fn sourceParameters() Parameters {
    return .{
        .carboxyl_dissociation_mol_per_m3 = 1.0e-2,
        .maximum_exchange_mol_per_megagram_iteration = 1.0e-2,
        .negligible_open_sites_mol_per_megagram = 1.0e-48,
    };
}

test "STARTE litter carboxyl iteration preserves nested source bounds" {
    const result = try iterate(.{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .protonated_carboxyl_sites_mol_per_megagram = 0.4,
        .hydrogen_activity_mol_per_m3 = 0.1,
    }, sourceParameters());
    try std.testing.expectEqual(@as(f64, 0.6), result.open_carboxyl_sites_mol_per_megagram);
    try std.testing.expectEqual(
        @as(f64, 0.04),
        result.equilibrium_open_carboxyl_sites_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @as(f64, 0.01),
        result.protonation_change_mol_per_megagram_iteration,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.41),
        result.protonated_carboxyl_sites_mol_per_megagram,
        1.0e-15,
    );
}

test "STARTE litter carboxyl desorption is limited by occupied fraction" {
    const result = try iterate(.{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .protonated_carboxyl_sites_mol_per_megagram = 0.05,
        .hydrogen_activity_mol_per_m3 = 1.0e-6,
    }, .{
        .carboxyl_dissociation_mol_per_m3 = 1,
        .maximum_exchange_mol_per_megagram_iteration = 0.2,
        .negligible_open_sites_mol_per_megagram = 1.0e-48,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.01),
        result.protonation_change_mol_per_megagram_iteration,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.04),
        result.protonated_carboxyl_sites_mol_per_megagram,
        1.0e-15,
    );
}

test "STARTE litter carboxyl iteration rejects invalid occupancy" {
    try std.testing.expectError(
        error.InvalidSurfaceCarboxylInitializationInput,
        iterate(.{
            .total_carboxyl_sites_mol_per_megagram = 1,
            .protonated_carboxyl_sites_mol_per_megagram = 2,
            .hydrogen_activity_mol_per_m3 = 1,
        }, sourceParameters()),
    );
}
