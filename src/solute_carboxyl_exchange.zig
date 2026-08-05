const std = @import("std");

pub const Parameters = struct {
    dissociation_constant_mol_per_m3: f64,
    maximum_exchange_mol_per_m3_per_iteration: f64,
    substrate_limit_fraction_per_iteration: f64,
};

pub const Inputs = struct {
    total_carboxyl_sites_mol_per_megagram: f64,
    hydrogen_occupied_sites_mol_per_megagram: f64,
    hydrogen_activity_mol_per_m3: f64,
    soil_mass_per_water_volume_megagrams_per_m3: f64,
};

pub const SurfaceControls = struct {
    /// Runtime replacement for source `ZEROC` in `XCOO`.
    minimum_open_sites_mol_per_megagram: f64,
};

pub const SurfaceSourceOrderResult = struct {
    substrate_limit_mol_per_megagram_step: f64,
    maximum_exchange_mol_per_megagram_step: f64,
    total_carboxyl_sites_mol_per_megagram: f64,
    open_sites_mol_per_megagram: f64,
    equilibrium_open_sites_mol_per_megagram: f64,
    hydrogen_adsorption_mol_per_megagram_step: f64,
    open_site_floor_was_applied: bool,
};

pub const SourceOrderControls = SurfaceControls;
pub const SourceOrderResult = SurfaceSourceOrderResult;

pub const SourceIterationStage = enum {
    before_iteration_ceiling,
    iteration_ceiling,
};

/// Directly translates SOLUTE.F lines 1407--1423 (`RXHC`). A positive result
/// protonates a carboxyl site and consumes the same amount of aqueous hydrogen
/// after conversion from mol/Mg to mol/m3 with the soil-mass:water ratio.
pub fn calculateChangeMolPerMg(inputs: Inputs, parameters: Parameters) !f64 {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidCarboxylExchangeInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidCarboxylExchangeParameter;
    }
    if (inputs.soil_mass_per_water_volume_megagrams_per_m3 == 0) return error.ZeroCarboxylExchangeSoilWaterRatio;
    // Occupied sites are clipped here; float drift above capacity is handled
    // by the @min so no hard error is needed.
    const occupied = @min(inputs.hydrogen_occupied_sites_mol_per_megagram, inputs.total_carboxyl_sites_mol_per_megagram);
    const deprotonated = @max(0.0, inputs.total_carboxyl_sites_mol_per_megagram - occupied);
    const equilibrium_deprotonated = if (inputs.hydrogen_activity_mol_per_m3 > 0)
        @min(inputs.total_carboxyl_sites_mol_per_megagram, parameters.dissociation_constant_mol_per_m3 * occupied / inputs.hydrogen_activity_mol_per_m3)
    else
        inputs.total_carboxyl_sites_mol_per_megagram;
    const substrate_limit = parameters.substrate_limit_fraction_per_iteration / inputs.soil_mass_per_water_volume_megagrams_per_m3 * occupied;
    const kinetic_limit = parameters.maximum_exchange_mol_per_m3_per_iteration / inputs.soil_mass_per_water_volume_megagrams_per_m3;
    const change = @max(-kinetic_limit, @max(-substrate_limit, @min(kinetic_limit, @min(substrate_limit, deprotonated - equilibrium_deprotonated))));
    if (!std.math.isFinite(change)) return error.NonFiniteCarboxylExchangeChange;
    return change;
}

/// Direct state update for SOLUTE.F 2401. No source floor is applied, and the
/// occupied carboxyl pool is unchanged at `M == MRXN`.
pub fn applySourceOrderStateUpdate(
    current_hydrogen_occupied_sites_mol_per_megagram: f64,
    change_mol_per_megagram: f64,
    stage: SourceIterationStage,
) !f64 {
    if (!std.math.isFinite(current_hydrogen_occupied_sites_mol_per_megagram) or
        current_hydrogen_occupied_sites_mol_per_megagram < 0 or
        !std.math.isFinite(change_mol_per_megagram))
        return error.InvalidCarboxylExchangeStateUpdate;
    if (stage == .iteration_ceiling)
        return current_hydrogen_occupied_sites_mol_per_megagram;
    const next = current_hydrogen_occupied_sites_mol_per_megagram + change_mol_per_megagram;
    if (!std.math.isFinite(next) or next < 0)
        return error.InvalidCarboxylExchangeStateUpdate;
    return next;
}

/// Direct source-order translation of the soil carboxyl dissociation block,
/// SOLUTE.F 1407--1423. State is read only; the returned positive extent
/// protonates COO- and is in mol Mg-1 per reaction iteration.
pub fn calculateSourceOrder(
    inputs: Inputs,
    parameters: Parameters,
    controls: SourceOrderControls,
) !SourceOrderResult {
    return calculateSourceOrderSurface(inputs, parameters, controls);
}

/// Direct source-order translation of SOLUTE.F lines 4483--4498.
///
/// A positive result protonates a surface-litter carboxyl site. The pure
/// scalar kernel owns no grid dimensions and mutates no state.
pub fn calculateSourceOrderSurface(
    inputs: Inputs,
    parameters: Parameters,
    controls: SurfaceControls,
) !SurfaceSourceOrderResult {
    try validateSurfaceInputs(inputs, parameters, controls);
    const occupied = inputs.hydrogen_occupied_sites_mol_per_megagram;
    const density = inputs.soil_mass_per_water_volume_megagrams_per_m3;

    // SOLUTE.F 4492--4498. Preserve division, floor, and nested-bound order.
    const substrate_limit =
        parameters.substrate_limit_fraction_per_iteration / density * occupied;
    const maximum_exchange =
        parameters.maximum_exchange_mol_per_m3_per_iteration / density;
    const total = inputs.total_carboxyl_sites_mol_per_megagram;
    const unconstrained_open = total - occupied;
    const open = @max(
        controls.minimum_open_sites_mol_per_megagram,
        unconstrained_open,
    );
    const equilibrium_open = @min(
        total,
        parameters.dissociation_constant_mol_per_m3 * occupied /
            inputs.hydrogen_activity_mol_per_m3,
    );
    const change = @max(
        -maximum_exchange,
        -substrate_limit,
        @min(
            maximum_exchange,
            substrate_limit,
            open - equilibrium_open,
        ),
    );
    const result: SurfaceSourceOrderResult = .{
        .substrate_limit_mol_per_megagram_step = substrate_limit,
        .maximum_exchange_mol_per_megagram_step = maximum_exchange,
        .total_carboxyl_sites_mol_per_megagram = total,
        .open_sites_mol_per_megagram = open,
        .equilibrium_open_sites_mol_per_megagram = equilibrium_open,
        .hydrogen_adsorption_mol_per_megagram_step = change,
        .open_site_floor_was_applied = unconstrained_open < controls.minimum_open_sites_mol_per_megagram,
    };
    inline for (@typeInfo(SurfaceSourceOrderResult).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .bool) continue;
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceCarboxylExchangeResult;
    }
    return result;
}

fn validateSurfaceInputs(
    inputs: Inputs,
    parameters: Parameters,
    controls: SurfaceControls,
) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceCarboxylExchangeInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceCarboxylExchangeParameter;
    }
    if (!std.math.isFinite(controls.minimum_open_sites_mol_per_megagram) or
        controls.minimum_open_sites_mol_per_megagram <= 0 or
        inputs.hydrogen_activity_mol_per_m3 <= 0 or
        inputs.soil_mass_per_water_volume_megagrams_per_m3 <= 0 or
        inputs.hydrogen_occupied_sites_mol_per_megagram >
            inputs.total_carboxyl_sites_mol_per_megagram or
        parameters.substrate_limit_fraction_per_iteration > 1)
    {
        return error.InvalidSurfaceCarboxylExchangeInput;
    }
}

test "SOLUTE carboxyl exchange reproduces bounded proton adsorption" {
    const change = try calculateChangeMolPerMg(.{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .hydrogen_occupied_sites_mol_per_megagram = 0.4,
        .hydrogen_activity_mol_per_m3 = 0.1,
        .soil_mass_per_water_volume_megagrams_per_m3 = 2,
    }, .{
        .dissociation_constant_mol_per_m3 = 0.01,
        .maximum_exchange_mol_per_m3_per_iteration = 0.2,
        .substrate_limit_fraction_per_iteration = 0.2,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.04), change, 1e-15);
}

test "source carboxyl state update is unfloored and ceiling gated" {
    try std.testing.expectEqual(
        @as(f64, 0),
        try applySourceOrderStateUpdate(
            0.25,
            -0.25,
            .before_iteration_ceiling,
        ),
    );
    try std.testing.expectEqual(
        @as(f64, 0.25),
        try applySourceOrderStateUpdate(0.25, -0.25, .iteration_ceiling),
    );
    try std.testing.expectError(
        error.InvalidCarboxylExchangeStateUpdate,
        applySourceOrderStateUpdate(0.25, -0.3, .before_iteration_ceiling),
    );
}

test "SOLUTE carboxyl exchange desorption cannot exceed occupied sites" {
    const change = try calculateChangeMolPerMg(.{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .hydrogen_occupied_sites_mol_per_megagram = 0.1,
        .hydrogen_activity_mol_per_m3 = 0,
        .soil_mass_per_water_volume_megagrams_per_m3 = 1,
    }, .{
        .dissociation_constant_mol_per_m3 = 0.01,
        .maximum_exchange_mol_per_m3_per_iteration = 1,
        .substrate_limit_fraction_per_iteration = 1,
    });
    try std.testing.expectApproxEqAbs(@as(f64, -0.1), change, 1e-15);
}

test "SOLUTE carboxyl exchange matches source nested rate bounds" {
    const inputs = Inputs{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .hydrogen_occupied_sites_mol_per_megagram = 0.2,
        .hydrogen_activity_mol_per_m3 = 0.01,
        .soil_mass_per_water_volume_megagrams_per_m3 = 2,
    };
    const parameters = Parameters{
        .dissociation_constant_mol_per_m3 = 0.08,
        .maximum_exchange_mol_per_m3_per_iteration = 0.5,
        .substrate_limit_fraction_per_iteration = 0.5,
    };

    // SOLUTE.F: XMIN=0.05, TADCC=0.25, XCOO=0.8, XCOOQ=1.0,
    // so RXHC=max(-0.25,-0.05,min(0.25,0.05,-0.2))=-0.05 mol C/Mg.
    const change = try calculateChangeMolPerMg(inputs, parameters);
    try std.testing.expectApproxEqAbs(@as(f64, -0.05), change, 1e-15);
}

test "soil carboxyl source order retains runtime open-site floor" {
    const inputs: Inputs = .{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .hydrogen_occupied_sites_mol_per_megagram = 1,
        .hydrogen_activity_mol_per_m3 = 2,
        .soil_mass_per_water_volume_megagrams_per_m3 = 4,
    };
    const parameters: Parameters = .{
        .dissociation_constant_mol_per_m3 = 0.2,
        .maximum_exchange_mol_per_m3_per_iteration = 0.8,
        .substrate_limit_fraction_per_iteration = 0.5,
    };
    const result = try calculateSourceOrder(
        inputs,
        parameters,
        .{ .minimum_open_sites_mol_per_megagram = 0.03 },
    );

    // SOLUTE.F 1418--1423: XMIN=.125, TADCC=.2, XCOO=.03,
    // XCOOQ=.1, RXHC=max(-.2,-.125,min(.2,.125,-.07))=-.07.
    try std.testing.expectEqual(@as(f64, 0.125), result.substrate_limit_mol_per_megagram_step);
    try std.testing.expectEqual(@as(f64, 0.2), result.maximum_exchange_mol_per_megagram_step);
    try std.testing.expectEqual(@as(f64, 0.03), result.open_sites_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 0.1), result.equilibrium_open_sites_mol_per_megagram);
    try std.testing.expectApproxEqAbs(
        @as(f64, -0.07),
        result.hydrogen_adsorption_mol_per_megagram_step,
        1e-15,
    );
}

fn surfaceTestInputs() Inputs {
    return .{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .hydrogen_occupied_sites_mol_per_megagram = 0.2,
        .hydrogen_activity_mol_per_m3 = 0.01,
        .soil_mass_per_water_volume_megagrams_per_m3 = 2,
    };
}

fn surfaceTestParameters() Parameters {
    return .{
        .dissociation_constant_mol_per_m3 = 0.08,
        .maximum_exchange_mol_per_m3_per_iteration = 0.5,
        .substrate_limit_fraction_per_iteration = 0.5,
    };
}

test "SOLUTE surface carboxyl exchange preserves every source expression" {
    const inputs = surfaceTestInputs();
    const parameters = surfaceTestParameters();
    const controls: SurfaceControls = .{
        .minimum_open_sites_mol_per_megagram = 1.0e-32,
    };
    const result =
        try calculateSourceOrderSurface(inputs, parameters, controls);
    const expected_substrate =
        parameters.substrate_limit_fraction_per_iteration /
        inputs.soil_mass_per_water_volume_megagrams_per_m3 *
        inputs.hydrogen_occupied_sites_mol_per_megagram;
    const expected_maximum =
        parameters.maximum_exchange_mol_per_m3_per_iteration /
        inputs.soil_mass_per_water_volume_megagrams_per_m3;
    const expected_open = @max(
        controls.minimum_open_sites_mol_per_megagram,
        inputs.total_carboxyl_sites_mol_per_megagram -
            inputs.hydrogen_occupied_sites_mol_per_megagram,
    );
    const expected_equilibrium = @min(
        inputs.total_carboxyl_sites_mol_per_megagram,
        parameters.dissociation_constant_mol_per_m3 *
            inputs.hydrogen_occupied_sites_mol_per_megagram /
            inputs.hydrogen_activity_mol_per_m3,
    );
    const expected_change = @max(
        -expected_maximum,
        -expected_substrate,
        @min(
            expected_maximum,
            expected_substrate,
            expected_open - expected_equilibrium,
        ),
    );

    try std.testing.expectEqual(
        expected_substrate,
        result.substrate_limit_mol_per_megagram_step,
    );
    try std.testing.expectEqual(
        expected_maximum,
        result.maximum_exchange_mol_per_megagram_step,
    );
    try std.testing.expectEqual(
        inputs.total_carboxyl_sites_mol_per_megagram,
        result.total_carboxyl_sites_mol_per_megagram,
    );
    try std.testing.expectEqual(expected_open, result.open_sites_mol_per_megagram);
    try std.testing.expectEqual(
        expected_equilibrium,
        result.equilibrium_open_sites_mol_per_megagram,
    );
    try std.testing.expectEqual(
        expected_change,
        result.hydrogen_adsorption_mol_per_megagram_step,
    );
}

test "surface carboxyl source floor is explicit" {
    var inputs = surfaceTestInputs();
    inputs.hydrogen_occupied_sites_mol_per_megagram =
        inputs.total_carboxyl_sites_mol_per_megagram;
    const controls: SurfaceControls = .{
        .minimum_open_sites_mol_per_megagram = 1.0e-12,
    };
    const result =
        try calculateSourceOrderSurface(inputs, surfaceTestParameters(), controls);
    try std.testing.expect(result.open_site_floor_was_applied);
    try std.testing.expectEqual(
        controls.minimum_open_sites_mol_per_megagram,
        result.open_sites_mol_per_megagram,
    );
}

test "surface source open-site floor can exceed a fully occupied capacity" {
    const inputs: Inputs = .{
        .total_carboxyl_sites_mol_per_megagram = 1,
        .hydrogen_occupied_sites_mol_per_megagram = 1,
        .hydrogen_activity_mol_per_m3 = 1,
        .soil_mass_per_water_volume_megagrams_per_m3 = 1,
    };
    const parameters: Parameters = .{
        .dissociation_constant_mol_per_m3 = 0,
        .maximum_exchange_mol_per_m3_per_iteration = 1,
        .substrate_limit_fraction_per_iteration = 1,
    };
    const source = try calculateSourceOrderSurface(
        inputs,
        parameters,
        .{ .minimum_open_sites_mol_per_megagram = 1.0e-6 },
    );
    const production_safe = try calculateChangeMolPerMg(inputs, parameters);

    try std.testing.expectEqual(
        @as(f64, 1.0e-6),
        source.hydrogen_adsorption_mol_per_megagram_step,
    );
    try std.testing.expectEqual(@as(f64, 0), production_safe);
    try std.testing.expect(
        inputs.hydrogen_occupied_sites_mol_per_megagram +
            source.hydrogen_adsorption_mol_per_megagram_step >
            inputs.total_carboxyl_sites_mol_per_megagram,
    );
}

test "surface carboxyl rate cannot overdraw occupied or open sites" {
    var inputs = surfaceTestInputs();
    var parameters = surfaceTestParameters();
    parameters.maximum_exchange_mol_per_m3_per_iteration = 100;
    parameters.substrate_limit_fraction_per_iteration = 1;
    inputs.hydrogen_activity_mol_per_m3 = 1.0e-12;
    const desorption = try calculateSourceOrderSurface(
        inputs,
        parameters,
        .{ .minimum_open_sites_mol_per_megagram = 1.0e-32 },
    );
    try std.testing.expect(
        -desorption.hydrogen_adsorption_mol_per_megagram_step <=
            inputs.hydrogen_occupied_sites_mol_per_megagram /
                inputs.soil_mass_per_water_volume_megagrams_per_m3,
    );

    inputs.hydrogen_activity_mol_per_m3 = 1.0e12;
    const adsorption = try calculateSourceOrderSurface(
        inputs,
        parameters,
        .{ .minimum_open_sites_mol_per_megagram = 1.0e-32 },
    );
    try std.testing.expect(
        adsorption.hydrogen_adsorption_mol_per_megagram_step <=
            inputs.total_carboxyl_sites_mol_per_megagram -
                inputs.hydrogen_occupied_sites_mol_per_megagram,
    );
}

test "surface carboxyl exchange rejects invalid input and overflow" {
    var inputs = surfaceTestInputs();
    inputs.hydrogen_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceCarboxylExchangeInput,
        calculateSourceOrderSurface(
            inputs,
            surfaceTestParameters(),
            .{ .minimum_open_sites_mol_per_megagram = 1.0e-32 },
        ),
    );

    inputs = surfaceTestInputs();
    inputs.hydrogen_occupied_sites_mol_per_megagram = 2;
    try std.testing.expectError(
        error.InvalidSurfaceCarboxylExchangeInput,
        calculateSourceOrderSurface(
            inputs,
            surfaceTestParameters(),
            .{ .minimum_open_sites_mol_per_megagram = 1.0e-32 },
        ),
    );

    inputs = surfaceTestInputs();
    inputs.soil_mass_per_water_volume_megagrams_per_m3 =
        std.math.floatMin(f64);
    var parameters = surfaceTestParameters();
    parameters.maximum_exchange_mol_per_m3_per_iteration =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceCarboxylExchangeResult,
        calculateSourceOrderSurface(
            inputs,
            parameters,
            .{ .minimum_open_sites_mol_per_megagram = 1.0e-32 },
        ),
    );
}
