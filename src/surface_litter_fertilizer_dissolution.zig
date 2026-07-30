const std = @import("std");

pub const Inputs = struct {
    ammonium_dissolution_fraction_per_step: f64,
    ammonia_dissolution_fraction_per_step: f64,
    nitrate_dissolution_fraction_per_step: f64,
    ammonium_fertilizer_mol_n: f64,
    ammonia_fertilizer_mol_n: f64,
    urea_fertilizer_mol_n: f64,
    nitrate_fertilizer_mol_n: f64,
    litter_water_concentration_m3_per_m3: f64,
    urea_hydrolysis_mol_n_per_step: f64,
};

pub const Fluxes = struct {
    ammonium_dissolution_mol_n_per_step: f64,
    ammonia_dissolution_mol_n_per_step: f64,
    urea_hydrolysis_mol_n_per_step: f64,
    nitrate_dissolution_mol_n_per_step: f64,
};

/// Direct source-order translation of SOLUTE.F lines 4106--4125.
///
/// The caller selects one runtime horizontal cell. This pure kernel computes
/// the four broadcast-fertilizer transfers without mutating their inventories.
pub fn calculate(inputs: Inputs) !Fluxes {
    try validateInputs(inputs);

    // SOLUTE.F 4122--4125. Preserve the source multiplication order.
    const fluxes: Fluxes = .{
        .ammonium_dissolution_mol_n_per_step = inputs.ammonium_dissolution_fraction_per_step *
            inputs.ammonium_fertilizer_mol_n *
            inputs.litter_water_concentration_m3_per_m3,
        .ammonia_dissolution_mol_n_per_step = inputs.ammonia_dissolution_fraction_per_step *
            inputs.ammonia_fertilizer_mol_n,
        .urea_hydrolysis_mol_n_per_step = inputs.urea_hydrolysis_mol_n_per_step,
        .nitrate_dissolution_mol_n_per_step = inputs.nitrate_dissolution_fraction_per_step *
            inputs.nitrate_fertilizer_mol_n *
            inputs.litter_water_concentration_m3_per_m3,
    };
    try validateFluxes(fluxes, inputs);
    return fluxes;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterFertilizerDissolutionInput;
    }
    if (inputs.litter_water_concentration_m3_per_m3 > 1 or
        inputs.urea_hydrolysis_mol_n_per_step >
            inputs.urea_fertilizer_mol_n)
    {
        return error.InvalidSurfaceLitterFertilizerDissolutionInput;
    }
}

fn validateFluxes(fluxes: Fluxes, inputs: Inputs) !void {
    inline for (@typeInfo(Fluxes).@"struct".fields) |field| {
        const value = @field(fluxes, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterFertilizerDissolution;
        if (value < 0)
            return error.InvalidSurfaceLitterFertilizerDissolutionFlux;
    }
    if (fluxes.ammonium_dissolution_mol_n_per_step >
        inputs.ammonium_fertilizer_mol_n or
        fluxes.ammonia_dissolution_mol_n_per_step >
            inputs.ammonia_fertilizer_mol_n or
        fluxes.urea_hydrolysis_mol_n_per_step >
            inputs.urea_fertilizer_mol_n or
        fluxes.nitrate_dissolution_mol_n_per_step >
            inputs.nitrate_fertilizer_mol_n)
    {
        return error.SurfaceLitterFertilizerDissolutionExceedsInventory;
    }
}

fn testInputs() Inputs {
    return .{
        .ammonium_dissolution_fraction_per_step = 0.2,
        .ammonia_dissolution_fraction_per_step = 0.25,
        .nitrate_dissolution_fraction_per_step = 0.4,
        .ammonium_fertilizer_mol_n = 10,
        .ammonia_fertilizer_mol_n = 8,
        .urea_fertilizer_mol_n = 6,
        .nitrate_fertilizer_mol_n = 5,
        .litter_water_concentration_m3_per_m3 = 0.5,
        .urea_hydrolysis_mol_n_per_step = 1.5,
    };
}

test "SOLUTE surface fertilizer dissolution preserves every source expression" {
    const inputs = testInputs();
    const fluxes = try calculate(inputs);

    const expected_ammonium =
        inputs.ammonium_dissolution_fraction_per_step *
        inputs.ammonium_fertilizer_mol_n *
        inputs.litter_water_concentration_m3_per_m3;
    const expected_ammonia =
        inputs.ammonia_dissolution_fraction_per_step *
        inputs.ammonia_fertilizer_mol_n;
    const expected_urea = inputs.urea_hydrolysis_mol_n_per_step;
    const expected_nitrate =
        inputs.nitrate_dissolution_fraction_per_step *
        inputs.nitrate_fertilizer_mol_n *
        inputs.litter_water_concentration_m3_per_m3;

    try std.testing.expectEqual(
        expected_ammonium,
        fluxes.ammonium_dissolution_mol_n_per_step,
    );
    try std.testing.expectEqual(
        expected_ammonia,
        fluxes.ammonia_dissolution_mol_n_per_step,
    );
    try std.testing.expectEqual(
        expected_urea,
        fluxes.urea_hydrolysis_mol_n_per_step,
    );
    try std.testing.expectEqual(
        expected_nitrate,
        fluxes.nitrate_dissolution_mol_n_per_step,
    );
}

test "surface fertilizer dissolution conserves each nitrogen donor" {
    const inputs = testInputs();
    const fluxes = try calculate(inputs);

    const initial_mol_n =
        inputs.ammonium_fertilizer_mol_n +
        inputs.ammonia_fertilizer_mol_n +
        inputs.urea_fertilizer_mol_n +
        inputs.nitrate_fertilizer_mol_n;
    const remaining_mol_n =
        inputs.ammonium_fertilizer_mol_n -
        fluxes.ammonium_dissolution_mol_n_per_step +
        inputs.ammonia_fertilizer_mol_n -
        fluxes.ammonia_dissolution_mol_n_per_step +
        inputs.urea_fertilizer_mol_n -
        fluxes.urea_hydrolysis_mol_n_per_step +
        inputs.nitrate_fertilizer_mol_n -
        fluxes.nitrate_dissolution_mol_n_per_step;
    const transferred_mol_n =
        fluxes.ammonium_dissolution_mol_n_per_step +
        fluxes.ammonia_dissolution_mol_n_per_step +
        fluxes.urea_hydrolysis_mol_n_per_step +
        fluxes.nitrate_dissolution_mol_n_per_step;

    try std.testing.expectApproxEqAbs(
        initial_mol_n,
        remaining_mol_n + transferred_mol_n,
        1.0e-14,
    );
}

test "dry litter suppresses water-mediated fertilizer dissolution only" {
    var inputs = testInputs();
    inputs.litter_water_concentration_m3_per_m3 = 0;
    const fluxes = try calculate(inputs);

    try std.testing.expectEqual(
        @as(f64, 0),
        fluxes.ammonium_dissolution_mol_n_per_step,
    );
    try std.testing.expectEqual(
        inputs.ammonia_dissolution_fraction_per_step *
            inputs.ammonia_fertilizer_mol_n,
        fluxes.ammonia_dissolution_mol_n_per_step,
    );
    try std.testing.expectEqual(
        inputs.urea_hydrolysis_mol_n_per_step,
        fluxes.urea_hydrolysis_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        fluxes.nitrate_dissolution_mol_n_per_step,
    );
}

test "surface fertilizer dissolution rejects invalid and unbounded transfers" {
    var inputs = testInputs();
    inputs.litter_water_concentration_m3_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterFertilizerDissolutionInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.urea_hydrolysis_mol_n_per_step =
        inputs.urea_fertilizer_mol_n + 1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterFertilizerDissolutionInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.ammonia_dissolution_fraction_per_step = 2;
    try std.testing.expectError(
        error.SurfaceLitterFertilizerDissolutionExceedsInventory,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.ammonium_dissolution_fraction_per_step =
        std.math.floatMax(f64);
    inputs.ammonium_fertilizer_mol_n = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterFertilizerDissolution,
        calculate(inputs),
    );
}
