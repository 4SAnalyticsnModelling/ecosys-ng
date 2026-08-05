const std = @import("std");

pub const RootProfile = enum {
    shallow,
    non_shallow,
};

pub const Inputs = struct {
    emerged: bool,
    biochemical_feedback_fraction: f64,
    solar_angle_sine: f64,
    absorbed_par_megajoules_per_timestep: f64,
    canopy_air_co2_umol_per_mol: f64,
    root_profile: RootProfile,
    stomatal_turgor_response: f64,
};

/// Exact GROSUB lines 949--953 branch-level admission to the canopy
/// carboxylation loops. Comparisons remain strict at zero.
pub fn admitted(inputs: Inputs) !bool {
    inline for (.{
        inputs.biochemical_feedback_fraction,
        inputs.solar_angle_sine,
        inputs.absorbed_par_megajoules_per_timestep,
        inputs.canopy_air_co2_umol_per_mol,
        inputs.stomatal_turgor_response,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteCanopyBranchCarboxylationAdmissionInput;
    if (inputs.biochemical_feedback_fraction < 0 or
        inputs.biochemical_feedback_fraction > 1 or
        inputs.solar_angle_sine < -1 or
        inputs.solar_angle_sine > 1 or
        inputs.absorbed_par_megajoules_per_timestep < 0 or
        inputs.canopy_air_co2_umol_per_mol < 0 or
        inputs.stomatal_turgor_response < 0)
        return error.InvalidCanopyBranchCarboxylationAdmissionInput;

    if (!inputs.emerged) return false;
    if (inputs.biochemical_feedback_fraction == 0) return false;
    if (inputs.solar_angle_sine <= 0) return false;
    if (inputs.absorbed_par_megajoules_per_timestep <= 0) return false;
    if (inputs.canopy_air_co2_umol_per_mol <= 0) return false;
    return inputs.root_profile != .shallow or
        inputs.stomatal_turgor_response > 0;
}

fn activeInputs() Inputs {
    return .{
        .emerged = true,
        .biochemical_feedback_fraction = 0.8,
        .solar_angle_sine = 0.5,
        .absorbed_par_megajoules_per_timestep = 1.2,
        .canopy_air_co2_umol_per_mol = 420,
        .root_profile = .non_shallow,
        .stomatal_turgor_response = 0.7,
    };
}

test "GROSUB admits a fully active emerged branch" {
    try std.testing.expect(try admitted(activeInputs()));
}

test "every source outer gate remains strict at zero" {
    inline for (.{
        "biochemical_feedback_fraction",
        "solar_angle_sine",
        "absorbed_par_megajoules_per_timestep",
        "canopy_air_co2_umol_per_mol",
    }) |field_name| {
        var inputs = activeInputs();
        @field(inputs, field_name) = 0;
        try std.testing.expect(!try admitted(inputs));
    }
    var pre_emergence = activeInputs();
    pre_emergence.emerged = false;
    try std.testing.expect(!try admitted(pre_emergence));
}

test "only shallow roots require positive stomatal turgor response" {
    var inputs = activeInputs();
    inputs.stomatal_turgor_response = 0;
    try std.testing.expect(try admitted(inputs));
    inputs.root_profile = .shallow;
    try std.testing.expect(!try admitted(inputs));
}

test "invalid admission inputs fail explicitly" {
    var inputs = activeInputs();
    inputs.canopy_air_co2_umol_per_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteCanopyBranchCarboxylationAdmissionInput,
        admitted(inputs),
    );
    inputs = activeInputs();
    inputs.biochemical_feedback_fraction = 1.1;
    try std.testing.expectError(
        error.InvalidCanopyBranchCarboxylationAdmissionInput,
        admitted(inputs),
    );
}
