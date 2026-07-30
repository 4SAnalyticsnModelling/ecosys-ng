const std = @import("std");

pub const PopulationInput = struct {
    initial_porosity_fraction: f64,
    maximum_primary_radius_m: f64,
    maximum_secondary_radius_m: f64,
};

pub const PopulationDimensions = struct {
    porosity_fraction: f64,
    radial_diffusion_path_log_factor: f64,
    volume_per_carbon_m3_g_c: f64,
    specific_primary_length_m_g_c: f64,
    specific_secondary_length_m_g_c: f64,
    primary_radius_m: f64,
    secondary_radius_m: f64,
    primary_cross_section_area_m2: f64,
    secondary_cross_section_area_m2: f64,
};

pub const Parameters = struct {
    minimum_porosity_for_diffusion: f64,
    cubic_metres_per_cubic_centimetre: f64,
    solid_carbon_density_g_c_cm3: f64,
    pi_approximation: f64,
};

pub const DerivationError = error{
    PopulationCountMismatch,
    NonFiniteInput,
    InvalidParameter,
    InvalidPorosity,
    InvalidRadius,
    NonFiniteResult,
};

/// Translates STARTQ lines 413-425 for runtime root/symbiont populations.
pub fn derive(
    inputs: []const PopulationInput,
    parameters: Parameters,
    outputs: []PopulationDimensions,
) DerivationError!void {
    if (inputs.len != outputs.len) return error.PopulationCountMismatch;
    try validate(inputs, parameters);

    for (inputs, outputs) |input, *output| {
        const porosity_fraction = input.initial_porosity_fraction;
        const radial_diffusion_path_log_factor =
            @log(1.0 / @sqrt(@max(
                parameters.minimum_porosity_for_diffusion,
                input.initial_porosity_fraction,
            )));
        const volume_per_carbon_m3_g_c = parameters.cubic_metres_per_cubic_centimetre /
            (parameters.solid_carbon_density_g_c_cm3 *
                (1.0 - input.initial_porosity_fraction));
        const specific_primary_length_m_g_c = volume_per_carbon_m3_g_c /
            (parameters.pi_approximation *
                input.maximum_primary_radius_m * input.maximum_primary_radius_m);
        const specific_secondary_length_m_g_c = volume_per_carbon_m3_g_c /
            (parameters.pi_approximation *
                input.maximum_secondary_radius_m * input.maximum_secondary_radius_m);
        const primary_radius_m = input.maximum_primary_radius_m;
        const secondary_radius_m = input.maximum_secondary_radius_m;
        const primary_cross_section_area_m2 = parameters.pi_approximation *
            primary_radius_m * primary_radius_m;
        const secondary_cross_section_area_m2 = parameters.pi_approximation *
            secondary_radius_m * secondary_radius_m;

        output.* = .{
            .porosity_fraction = porosity_fraction,
            .radial_diffusion_path_log_factor = radial_diffusion_path_log_factor,
            .volume_per_carbon_m3_g_c = volume_per_carbon_m3_g_c,
            .specific_primary_length_m_g_c = specific_primary_length_m_g_c,
            .specific_secondary_length_m_g_c = specific_secondary_length_m_g_c,
            .primary_radius_m = primary_radius_m,
            .secondary_radius_m = secondary_radius_m,
            .primary_cross_section_area_m2 = primary_cross_section_area_m2,
            .secondary_cross_section_area_m2 = secondary_cross_section_area_m2,
        };
        inline for (std.meta.fields(PopulationDimensions)) |field| {
            if (!std.math.isFinite(@field(output, field.name))) return error.NonFiniteResult;
        }
    }
}

fn validate(
    inputs: []const PopulationInput,
    parameters: Parameters,
) DerivationError!void {
    inline for (std.meta.fields(Parameters)) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value <= 0.0) return error.InvalidParameter;
    }
    if (parameters.minimum_porosity_for_diffusion > 1.0) {
        return error.InvalidParameter;
    }
    for (inputs) |input| {
        inline for (std.meta.fields(PopulationInput)) |field| {
            if (!std.math.isFinite(@field(input, field.name))) return error.NonFiniteInput;
        }
        if (input.initial_porosity_fraction < 0.0 or
            input.initial_porosity_fraction >= 1.0)
        {
            return error.InvalidPorosity;
        }
        if (input.maximum_primary_radius_m <= 0.0 or
            input.maximum_secondary_radius_m <= 0.0)
        {
            return error.InvalidRadius;
        }
    }
}

fn legacyParameters() Parameters {
    return .{
        .minimum_porosity_for_diffusion = 0.01,
        .cubic_metres_per_cubic_centimetre = 1.0e-6,
        .solid_carbon_density_g_c_cm3 = 0.05,
        .pi_approximation = 3.142,
    };
}

test "runtime root and mycorrhizal dimensions preserve STARTQ equations" {
    const inputs = [_]PopulationInput{
        .{
            .initial_porosity_fraction = 0.2,
            .maximum_primary_radius_m = 1.0e-4,
            .maximum_secondary_radius_m = 5.0e-5,
        },
        .{
            .initial_porosity_fraction = 0.4,
            .maximum_primary_radius_m = 2.5e-6,
            .maximum_secondary_radius_m = 2.5e-6,
        },
    };
    var outputs: [inputs.len]PopulationDimensions = undefined;
    try derive(&inputs, legacyParameters(), &outputs);

    const expected_volume = 1.0e-6 / (0.05 * (1.0 - 0.2));
    try std.testing.expectApproxEqRel(
        @log(1.0 / @sqrt(0.2)),
        outputs[0].radial_diffusion_path_log_factor,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        expected_volume / (3.142 * 1.0e-4 * 1.0e-4),
        outputs[0].specific_primary_length_m_g_c,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        3.142 * 2.5e-6 * 2.5e-6,
        outputs[1].primary_cross_section_area_m2,
        1.0e-14,
    );
}

test "porosity at one fails before output mutation" {
    const inputs = [_]PopulationInput{.{
        .initial_porosity_fraction = 1.0,
        .maximum_primary_radius_m = 1.0e-4,
        .maximum_secondary_radius_m = 1.0e-4,
    }};
    var outputs = [_]PopulationDimensions{std.mem.zeroes(PopulationDimensions)};
    outputs[0].primary_radius_m = 42.0;
    try std.testing.expectError(
        error.InvalidPorosity,
        derive(&inputs, legacyParameters(), &outputs),
    );
    try std.testing.expectEqual(@as(f64, 42.0), outputs[0].primary_radius_m);
}
