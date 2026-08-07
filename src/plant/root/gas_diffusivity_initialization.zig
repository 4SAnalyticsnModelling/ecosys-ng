const std = @import("std");

pub const SpeciesDiffusivity = struct {
    carbon_dioxide_m2_per_h: f64,
    oxygen_m2_per_h: f64,
    methane_m2_per_h: f64,
    nitrous_oxide_m2_per_h: f64,
    ammonia_m2_per_h: f64,
    hydrogen_m2_per_h: f64,
};

pub const SpeciesDiffusivityPerStep = struct {
    carbon_dioxide_m2_per_step: f64,
    oxygen_m2_per_step: f64,
    methane_m2_per_step: f64,
    nitrous_oxide_m2_per_step: f64,
    ammonia_m2_per_step: f64,
    hydrogen_m2_per_step: f64,
};

pub const GasFluxWorkspace = struct {
    carbon_dioxide: f64 = 0,
    oxygen: f64 = 0,
    methane: f64 = 0,
    nitrous_oxide: f64 = 0,
    ammonia: f64 = 0,
    hydrogen: f64 = 0,
};

pub const Inputs = struct {
    root_porosity_m3_per_m3: f64,
    porosity_tortuosity_exponent: f64,
    gas_transfer_timestep_h_per_step: f64,
    oxygen_competition_fraction: f64,
    gaseous_diffusivity: SpeciesDiffusivity,
    aqueous_diffusivity: SpeciesDiffusivity,
};

pub const Result = struct {
    porosity_tortuosity_factor: f64,
    gaseous_diffusivity: SpeciesDiffusivityPerStep,
    allocated_aqueous_diffusivity: SpeciesDiffusivityPerStep,
    unallocated_oxygen_aqueous_diffusivity_m2_per_step: f64,
    gas_flux_workspace: GasFluxWorkspace,
};

/// UPTAKE.F 1943--1962. Scales root/soil gas diffusivities and clears the
/// six local gas-flux workspaces in exact source order.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const porosity_factor = std.math.pow(
        f64,
        inputs.root_porosity_m3_per_m3,
        inputs.porosity_tortuosity_exponent,
    );
    var gaseous: SpeciesDiffusivityPerStep = undefined;
    gaseous.carbon_dioxide_m2_per_step =
        inputs.gaseous_diffusivity.carbon_dioxide_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step * porosity_factor;
    gaseous.oxygen_m2_per_step =
        inputs.gaseous_diffusivity.oxygen_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step * porosity_factor;
    gaseous.methane_m2_per_step =
        inputs.gaseous_diffusivity.methane_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step * porosity_factor;
    gaseous.nitrous_oxide_m2_per_step =
        inputs.gaseous_diffusivity.nitrous_oxide_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step * porosity_factor;
    gaseous.ammonia_m2_per_step =
        inputs.gaseous_diffusivity.ammonia_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step * porosity_factor;
    gaseous.hydrogen_m2_per_step =
        inputs.gaseous_diffusivity.hydrogen_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step * porosity_factor;

    var aqueous: SpeciesDiffusivityPerStep = undefined;
    aqueous.carbon_dioxide_m2_per_step =
        inputs.aqueous_diffusivity.carbon_dioxide_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step *
        inputs.oxygen_competition_fraction;
    aqueous.oxygen_m2_per_step =
        inputs.aqueous_diffusivity.oxygen_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step *
        inputs.oxygen_competition_fraction;
    aqueous.methane_m2_per_step =
        inputs.aqueous_diffusivity.methane_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step *
        inputs.oxygen_competition_fraction;
    aqueous.nitrous_oxide_m2_per_step =
        inputs.aqueous_diffusivity.nitrous_oxide_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step *
        inputs.oxygen_competition_fraction;
    aqueous.ammonia_m2_per_step =
        inputs.aqueous_diffusivity.ammonia_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step *
        inputs.oxygen_competition_fraction;
    aqueous.hydrogen_m2_per_step =
        inputs.aqueous_diffusivity.hydrogen_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step *
        inputs.oxygen_competition_fraction;
    const unallocated_oxygen =
        inputs.aqueous_diffusivity.oxygen_m2_per_h *
        inputs.gas_transfer_timestep_h_per_step;
    const result = Result{
        .porosity_tortuosity_factor = porosity_factor,
        .gaseous_diffusivity = gaseous,
        .allocated_aqueous_diffusivity = aqueous,
        .unallocated_oxygen_aqueous_diffusivity_m2_per_step = unallocated_oxygen,
        .gas_flux_workspace = .{},
    };
    try validateResult(result);
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.root_porosity_m3_per_m3,
        inputs.porosity_tortuosity_exponent,
        inputs.gas_transfer_timestep_h_per_step,
        inputs.oxygen_competition_fraction,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidRootGasDiffusivityInput;
    if (inputs.root_porosity_m3_per_m3 < 0 or
        inputs.root_porosity_m3_per_m3 > 1 or
        inputs.porosity_tortuosity_exponent <= 0 or
        inputs.gas_transfer_timestep_h_per_step < 0 or
        inputs.oxygen_competition_fraction < 0)
        return error.InvalidRootGasDiffusivityInput;
    inline for (@typeInfo(SpeciesDiffusivity).@"struct".fields) |field| {
        const gaseous = @field(inputs.gaseous_diffusivity, field.name);
        const aqueous = @field(inputs.aqueous_diffusivity, field.name);
        if (!std.math.isFinite(gaseous) or gaseous < 0 or
            !std.math.isFinite(aqueous) or aqueous < 0)
            return error.InvalidRootGasDiffusivityInput;
    }
}

fn validateResult(result: Result) !void {
    if (!std.math.isFinite(result.porosity_tortuosity_factor) or
        !std.math.isFinite(
            result.unallocated_oxygen_aqueous_diffusivity_m2_per_step,
        ))
        return error.NonFiniteRootGasDiffusivityResult;
    inline for (@typeInfo(SpeciesDiffusivityPerStep).@"struct".fields) |field|
        if (!std.math.isFinite(
            @field(result.gaseous_diffusivity, field.name),
        ) or !std.math.isFinite(
            @field(result.allocated_aqueous_diffusivity, field.name),
        ))
            return error.NonFiniteRootGasDiffusivityResult;
}

fn species(base: f64) SpeciesDiffusivity {
    return .{
        .carbon_dioxide_m2_per_h = base + 1,
        .oxygen_m2_per_h = base + 2,
        .methane_m2_per_h = base + 3,
        .nitrous_oxide_m2_per_h = base + 4,
        .ammonia_m2_per_h = base + 5,
        .hydrogen_m2_per_h = base + 6,
    };
}

test "UPTAKE gaseous and aqueous diffusivity scaling preserves source factors" {
    const inputs = Inputs{
        .root_porosity_m3_per_m3 = 0.5,
        .porosity_tortuosity_exponent = 1.33,
        .gas_transfer_timestep_h_per_step = 0.25,
        .oxygen_competition_fraction = 0.4,
        .gaseous_diffusivity = species(0),
        .aqueous_diffusivity = species(10),
    };
    const result = try calculate(inputs);
    const porosity_factor = std.math.pow(f64, 0.5, 1.33);
    try std.testing.expectEqual(
        1.0 * 0.25 * porosity_factor,
        result.gaseous_diffusivity.carbon_dioxide_m2_per_step,
    );
    try std.testing.expectApproxEqAbs(
        12.0 * 0.25 * 0.4,
        result.allocated_aqueous_diffusivity.oxygen_m2_per_step,
        3e-16,
    );
    try std.testing.expectEqual(
        @as(f64, 12.0 * 0.25),
        result.unallocated_oxygen_aqueous_diffusivity_m2_per_step,
    );
}

test "UPTAKE six gas flux workspaces are reset to zero" {
    const result = try calculate(.{
        .root_porosity_m3_per_m3 = 0,
        .porosity_tortuosity_exponent = 1.33,
        .gas_transfer_timestep_h_per_step = 0.25,
        .oxygen_competition_fraction = 0.4,
        .gaseous_diffusivity = species(0),
        .aqueous_diffusivity = species(10),
    });
    inline for (@typeInfo(GasFluxWorkspace).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(result.gas_flux_workspace, field.name),
        );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.gaseous_diffusivity.oxygen_m2_per_step,
    );
}

test "out of range root porosity fails explicitly" {
    try std.testing.expectError(
        error.InvalidRootGasDiffusivityInput,
        calculate(.{
            .root_porosity_m3_per_m3 = -0.1,
            .porosity_tortuosity_exponent = 1.33,
            .gas_transfer_timestep_h_per_step = 0.25,
            .oxygen_competition_fraction = 0.4,
            .gaseous_diffusivity = species(0),
            .aqueous_diffusivity = species(10),
        }),
    );
}
