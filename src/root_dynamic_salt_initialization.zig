const std = @import("std");

pub const IonValues = struct {
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
    sulfate: f64,
    chloride: f64,
};

pub const Inputs = struct {
    dynamic_salts_enabled: bool,
    soil_content_mol: IonValues,
    root_content_mol: IonValues,
    diffusivity_m2_per_h: IonValues,
    root_mass_fraction: f64,
    water_flux_timestep_h_per_step: f64,
};

pub const Active = struct {
    allocated_soil_content_mol: IonValues,
    root_content_mol: IonValues,
    diffusivity_m2_per_step: IonValues,
};

pub const Outcome = union(enum) {
    inactive,
    active: Active,
};

/// UPTAKE.F 1977--2002. Initializes the eight dynamically solved root/soil
/// salts in exact Al, Fe, Ca, Mg, Na, K, sulfate, chloride source order.
pub fn calculate(inputs: Inputs) !Outcome {
    if (!inputs.dynamic_salts_enabled) return .inactive;
    try validate(inputs);
    var allocated_soil: IonValues = undefined;
    var root_content: IonValues = undefined;
    var diffusivity: IonValues = undefined;
    inline for (@typeInfo(IonValues).@"struct".fields) |field| {
        @field(allocated_soil, field.name) =
            @field(inputs.soil_content_mol, field.name) *
            inputs.root_mass_fraction;
        @field(root_content, field.name) =
            @field(inputs.root_content_mol, field.name);
        @field(diffusivity, field.name) =
            @field(inputs.diffusivity_m2_per_h, field.name) *
            inputs.water_flux_timestep_h_per_step;
    }
    inline for (@typeInfo(IonValues).@"struct".fields) |field|
        if (!std.math.isFinite(@field(allocated_soil, field.name)) or
            !std.math.isFinite(@field(diffusivity, field.name)))
            return error.NonFiniteRootDynamicSaltResult;
    return .{ .active = .{
        .allocated_soil_content_mol = allocated_soil,
        .root_content_mol = root_content,
        .diffusivity_m2_per_step = diffusivity,
    } };
}

fn validate(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.root_mass_fraction) or
        !std.math.isFinite(inputs.water_flux_timestep_h_per_step) or
        inputs.root_mass_fraction < 0 or
        inputs.water_flux_timestep_h_per_step < 0)
        return error.InvalidRootDynamicSaltInput;
    inline for (@typeInfo(IonValues).@"struct".fields) |field| {
        const soil = @field(inputs.soil_content_mol, field.name);
        const root = @field(inputs.root_content_mol, field.name);
        const diffusivity = @field(inputs.diffusivity_m2_per_h, field.name);
        if (!std.math.isFinite(soil) or soil < 0 or
            !std.math.isFinite(root) or root < 0 or
            !std.math.isFinite(diffusivity) or diffusivity < 0)
            return error.InvalidRootDynamicSaltInput;
    }
}

fn ions(base: f64) IonValues {
    return .{
        .aluminum = base + 1,
        .iron = base + 2,
        .calcium = base + 3,
        .magnesium = base + 4,
        .sodium = base + 5,
        .potassium = base + 6,
        .sulfate = base + 7,
        .chloride = base + 8,
    };
}

test "UPTAKE dynamic salts preserve allocation copy and timestep scaling" {
    const outcome = try calculate(.{
        .dynamic_salts_enabled = true,
        .soil_content_mol = ions(0),
        .root_content_mol = ions(10),
        .diffusivity_m2_per_h = ions(20),
        .root_mass_fraction = 0.25,
        .water_flux_timestep_h_per_step = 0.5,
    });
    const active = outcome.active;
    try std.testing.expectEqual(@as(f64, 0.25), active.allocated_soil_content_mol.aluminum);
    try std.testing.expectEqual(@as(f64, 18), active.root_content_mol.chloride);
    try std.testing.expectEqual(@as(f64, 13.5), active.diffusivity_m2_per_step.sulfate);
}

test "disabled dynamic salts do not evaluate inactive values" {
    var soil = ions(0);
    soil.chloride = std.math.nan(f64);
    const outcome = try calculate(.{
        .dynamic_salts_enabled = false,
        .soil_content_mol = soil,
        .root_content_mol = ions(10),
        .diffusivity_m2_per_h = ions(20),
        .root_mass_fraction = -1,
        .water_flux_timestep_h_per_step = -1,
    });
    try std.testing.expect(outcome == .inactive);
}

test "negative active salt content fails explicitly" {
    var soil = ions(0);
    soil.aluminum = -1;
    try std.testing.expectError(
        error.InvalidRootDynamicSaltInput,
        calculate(.{
            .dynamic_salts_enabled = true,
            .soil_content_mol = soil,
            .root_content_mol = ions(10),
            .diffusivity_m2_per_h = ions(20),
            .root_mass_fraction = 0.25,
            .water_flux_timestep_h_per_step = 0.5,
        }),
    );
}
