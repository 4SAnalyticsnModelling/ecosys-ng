const std = @import("std");

/// Original van Genuchten (1980) retention curve coupled to the Mualem
/// conductivity model. This intentionally does not implement the Ippisch
/// near-saturation modification.
pub const MualemVanGenuchtenParameters = struct {
    residual_water_content_m3_per_m3: f64,
    saturated_water_content_m3_per_m3: f64,
    alpha_per_m: f64,
    n: f64,
    pore_connectivity: f64 = 0.5,
    saturated_hydraulic_conductivity_m_per_h: f64,

    pub fn validate(self: MualemVanGenuchtenParameters) !void {
        inline for (@typeInfo(MualemVanGenuchtenParameters).@"struct".fields) |field| {
            if (!std.math.isFinite(@field(self, field.name)))
                return error.NonFiniteMualemVanGenuchtenParameter;
        }
        if (self.residual_water_content_m3_per_m3 < 0 or
            self.saturated_water_content_m3_per_m3 <= self.residual_water_content_m3_per_m3 or
            self.saturated_water_content_m3_per_m3 > 1 or
            self.alpha_per_m <= 0 or self.n <= 1 or
            self.saturated_hydraulic_conductivity_m_per_h < 0)
        {
            return error.InvalidMualemVanGenuchtenParameter;
        }
    }

    pub fn m(self: MualemVanGenuchtenParameters) f64 {
        return 1.0 - 1.0 / self.n;
    }

    pub fn effectiveSaturationAtPressureHead(
        self: MualemVanGenuchtenParameters,
        pressure_head_m: f64,
    ) !f64 {
        try self.validate();
        if (!std.math.isFinite(pressure_head_m))
            return error.NonFinitePressureHead;
        if (pressure_head_m >= 0) return 1;
        const scaled_head = self.alpha_per_m * -pressure_head_m;
        const effective_saturation =
            std.math.pow(f64, 1.0 + std.math.pow(f64, scaled_head, self.n), -self.m());
        if (!std.math.isFinite(effective_saturation))
            return error.NonFiniteEffectiveSaturation;
        return std.math.clamp(effective_saturation, 0, 1);
    }

    pub fn waterContentAtPressureHead(
        self: MualemVanGenuchtenParameters,
        pressure_head_m: f64,
    ) !f64 {
        const effective_saturation =
            try self.effectiveSaturationAtPressureHead(pressure_head_m);
        return self.residual_water_content_m3_per_m3 +
            effective_saturation *
                (self.saturated_water_content_m3_per_m3 -
                    self.residual_water_content_m3_per_m3);
    }

    pub fn pressureHeadAtWaterContent(
        self: MualemVanGenuchtenParameters,
        water_content_m3_per_m3: f64,
    ) !f64 {
        try self.validate();
        if (!std.math.isFinite(water_content_m3_per_m3))
            return error.NonFiniteWaterContent;
        if (water_content_m3_per_m3 < self.residual_water_content_m3_per_m3 or
            water_content_m3_per_m3 > self.saturated_water_content_m3_per_m3)
        {
            return error.WaterContentOutsideRetentionDomain;
        }
        if (water_content_m3_per_m3 == self.saturated_water_content_m3_per_m3)
            return 0;
        const effective_saturation = std.math.clamp(
            (water_content_m3_per_m3 - self.residual_water_content_m3_per_m3) /
                (self.saturated_water_content_m3_per_m3 -
                    self.residual_water_content_m3_per_m3),
            std.math.floatMin(f64),
            1,
        );
        // Evaluate the original van Genuchten inverse in the log domain.
        // At theta_r its mathematical limit is negative infinity; represent
        // that asymptote by the largest pressure magnitude for which both the
        // head and alpha*head remain finite. This avoids silently producing
        // infinity while retaining the unmodified constitutive curve.
        const inverse_saturation_log =
            -@log(effective_saturation) / self.m();
        const inner_log =
            if (inverse_saturation_log > 40)
                inverse_saturation_log
            else
                @log(@exp(inverse_saturation_log) - 1);
        const pressure_magnitude_log =
            inner_log / self.n - @log(self.alpha_per_m);
        const maximum_pressure_magnitude_log =
            @log(std.math.floatMax(f64)) -
            @max(0, @log(self.alpha_per_m));
        const pressure_head_m = -@exp(@min(
            pressure_magnitude_log,
            maximum_pressure_magnitude_log,
        ));
        if (!std.math.isFinite(pressure_head_m))
            return error.NonFinitePressureHead;
        return pressure_head_m;
    }

    pub fn waterCapacityPerM(
        self: MualemVanGenuchtenParameters,
        pressure_head_m: f64,
    ) !f64 {
        try self.validate();
        if (!std.math.isFinite(pressure_head_m))
            return error.NonFinitePressureHead;
        if (pressure_head_m >= 0) return 0;
        const absolute_head_m = -pressure_head_m;
        const scaled_head = self.alpha_per_m * absolute_head_m;
        const scaled_to_n = std.math.pow(f64, scaled_head, self.n);
        const capacity_per_m =
            (self.saturated_water_content_m3_per_m3 -
                self.residual_water_content_m3_per_m3) *
            self.m() * self.n * self.alpha_per_m *
            std.math.pow(f64, scaled_head, self.n - 1.0) *
            std.math.pow(f64, 1.0 + scaled_to_n, -self.m() - 1.0);
        if (!std.math.isFinite(capacity_per_m) or capacity_per_m < 0)
            return error.NonFiniteWaterCapacity;
        return capacity_per_m;
    }

    pub fn relativeHydraulicConductivityAtEffectiveSaturation(
        self: MualemVanGenuchtenParameters,
        effective_saturation: f64,
    ) !f64 {
        try self.validate();
        if (!std.math.isFinite(effective_saturation) or
            effective_saturation < 0 or effective_saturation > 1)
        {
            return error.InvalidEffectiveSaturation;
        }
        if (effective_saturation == 0) return 0;
        if (effective_saturation == 1) return 1;
        const saturation_to_inverse_m =
            std.math.pow(f64, effective_saturation, 1.0 / self.m());
        const mualem_integral =
            1.0 - std.math.pow(f64, 1.0 - saturation_to_inverse_m, self.m());
        const relative_conductivity =
            std.math.pow(f64, effective_saturation, self.pore_connectivity) *
            mualem_integral * mualem_integral;
        if (!std.math.isFinite(relative_conductivity) or
            relative_conductivity < 0 or relative_conductivity > 1)
        {
            return error.NonFiniteRelativeHydraulicConductivity;
        }
        return relative_conductivity;
    }

    pub fn hydraulicConductivityMPerH(
        self: MualemVanGenuchtenParameters,
        pressure_head_m: f64,
    ) !f64 {
        const effective_saturation =
            try self.effectiveSaturationAtPressureHead(pressure_head_m);
        return self.saturated_hydraulic_conductivity_m_per_h *
            try self.relativeHydraulicConductivityAtEffectiveSaturation(
                effective_saturation,
            );
    }
};

pub const MualemVanGenuchtenFitInputs = struct {
    saturated_water_content_m3_per_m3: f64,
    field_capacity_water_content_m3_per_m3: f64,
    field_capacity_pressure_head_m: f64,
    wilting_point_water_content_m3_per_m3: f64,
    wilting_point_pressure_head_m: f64,
    inflection_pressure_head_m: f64,
    saturated_hydraulic_conductivity_m_per_h: f64,
    pore_connectivity: f64 = 0.5,
};

pub const MualemVanGenuchtenFitOptions = struct {
    minimum_n: f64 = 1.01,
    maximum_n: f64 = 20,
    water_content_tolerance_m3_per_m3: f64 = 1.0e-10,
    maximum_iterations: u16,

    pub fn validate(self: MualemVanGenuchtenFitOptions) !void {
        inline for (@typeInfo(MualemVanGenuchtenFitOptions).@"struct".fields) |field| {
            if (field.type == f64 and !std.math.isFinite(@field(self, field.name)))
                return error.NonFiniteMualemVanGenuchtenFitOption;
        }
        if (self.minimum_n <= 1 or self.maximum_n <= self.minimum_n or
            self.water_content_tolerance_m3_per_m3 <= 0 or
            self.maximum_iterations == 0)
        {
            return error.InvalidMualemVanGenuchtenFitOption;
        }
    }
};

pub const MualemVanGenuchtenFitResult = struct {
    parameters: MualemVanGenuchtenParameters,
    iterations: u16,
    newton_raphson_steps: u16,
    bisection_steps: u16,
    residual_water_content_m3_per_m3: f64,
};

pub const SoilTextureClass = enum {
    sand,
    loamy_sand,
    sandy_loam,
    loam,
    silt_loam,
    silt,
    sandy_clay_loam,
    clay_loam,
    silty_clay_loam,
    sandy_clay,
    silty_clay,
    clay,
};

pub fn classifyUsdaSoilTexture(
    sand_fraction: f64,
    silt_fraction: f64,
    clay_fraction: f64,
) !SoilTextureClass {
    inline for (.{ sand_fraction, silt_fraction, clay_fraction }) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidSoilTextureFraction;
    }
    const mineral_fraction = sand_fraction + silt_fraction + clay_fraction;
    if (!std.math.isFinite(mineral_fraction) or mineral_fraction <= 0 or
        mineral_fraction > 1.0 + 1.0e-8)
        return error.SoilTextureFractionsDoNotSumToOne;
    // USDA texture is defined on the fine-earth mineral fraction. Soil input
    // mass fractions can sum below one when organic matter or coarse material
    // is present, so normalize only for classification.
    const normalized_sand_fraction = sand_fraction / mineral_fraction;
    const normalized_silt_fraction = silt_fraction / mineral_fraction;
    const normalized_clay_fraction = clay_fraction / mineral_fraction;
    if (normalized_clay_fraction >= 0.4) {
        if (normalized_silt_fraction >= 0.4) return .silty_clay;
        if (normalized_sand_fraction >= 0.45) return .sandy_clay;
        return .clay;
    }
    if (normalized_clay_fraction >= 0.27) {
        if (normalized_sand_fraction >= 0.45) return .sandy_clay;
        if (normalized_sand_fraction <= 0.2) return .silty_clay_loam;
        return .clay_loam;
    }
    if (normalized_clay_fraction >= 0.2 and normalized_sand_fraction >= 0.45)
        return .sandy_clay_loam;
    if (normalized_silt_fraction >= 0.8 and normalized_clay_fraction < 0.12) return .silt;
    if (normalized_silt_fraction >= 0.5 and normalized_sand_fraction < 0.5) return .silt_loam;
    if (normalized_sand_fraction >= 0.85 and normalized_clay_fraction < 0.1) return .sand;
    if (normalized_sand_fraction >= 0.7 and normalized_clay_fraction < 0.15) return .loamy_sand;
    if (normalized_sand_fraction >= 0.43 and normalized_clay_fraction < 0.2) return .sandy_loam;
    if (normalized_clay_fraction >= 0.27) return .clay_loam;
    return .loam;
}

const CarselParrishRow = struct {
    residual_water_content_m3_per_m3: f64,
    saturated_water_content_m3_per_m3: f64,
    alpha_per_cm: f64,
    n: f64,
    saturated_hydraulic_conductivity_cm_per_day: f64,
};

fn carselParrishRow(
    residual_water_content_m3_per_m3: f64,
    saturated_water_content_m3_per_m3: f64,
    alpha_per_cm: f64,
    n: f64,
    saturated_hydraulic_conductivity_cm_per_day: f64,
) CarselParrishRow {
    return .{
        .residual_water_content_m3_per_m3 = residual_water_content_m3_per_m3,
        .saturated_water_content_m3_per_m3 = saturated_water_content_m3_per_m3,
        .alpha_per_cm = alpha_per_cm,
        .n = n,
        .saturated_hydraulic_conductivity_cm_per_day = saturated_hydraulic_conductivity_cm_per_day,
    };
}

/// Carsel-Parrish-type mean parameters copied into ordinary runtime data.
/// Alpha is converted from cm^-1 to m^-1 and Ksat from cm day^-1 to m h^-1.
pub fn carselParrishDefault(
    texture: SoilTextureClass,
    saturated_water_content_m3_per_m3: ?f64,
) !MualemVanGenuchtenParameters {
    const row: CarselParrishRow = switch (texture) {
        .sand => carselParrishRow(0.045, 0.43, 0.145, 2.68, 712.8),
        .loamy_sand => carselParrishRow(0.057, 0.41, 0.124, 2.28, 350.2),
        .sandy_loam => carselParrishRow(0.065, 0.41, 0.075, 1.89, 106.1),
        .loam => carselParrishRow(0.078, 0.43, 0.036, 1.56, 24.96),
        .silt_loam => carselParrishRow(0.067, 0.45, 0.020, 1.41, 10.8),
        .silt => carselParrishRow(0.034, 0.46, 0.016, 1.37, 6.0),
        .sandy_clay_loam => carselParrishRow(0.100, 0.39, 0.059, 1.48, 31.44),
        .clay_loam => carselParrishRow(0.095, 0.41, 0.019, 1.31, 6.24),
        .silty_clay_loam => carselParrishRow(0.089, 0.43, 0.010, 1.23, 1.68),
        .sandy_clay => carselParrishRow(0.100, 0.38, 0.027, 1.23, 2.88),
        .silty_clay => carselParrishRow(0.070, 0.36, 0.005, 1.09, 0.48),
        .clay => carselParrishRow(0.068, 0.38, 0.008, 1.09, 4.80),
    };
    const theta_s = saturated_water_content_m3_per_m3 orelse
        row.saturated_water_content_m3_per_m3;
    const parameters: MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = row.residual_water_content_m3_per_m3,
        .saturated_water_content_m3_per_m3 = theta_s,
        .alpha_per_m = 100.0 * row.alpha_per_cm,
        .n = row.n,
        .saturated_hydraulic_conductivity_m_per_h = row.saturated_hydraulic_conductivity_cm_per_day / 2400.0,
    };
    try parameters.validate();
    return parameters;
}

const RetentionFitEvaluation = struct {
    residual: f64,
    residual_water_content_m3_per_m3: f64,
    alpha_per_m: f64,
};

/// Fits the original van Genuchten curve. For an inflection in water content
/// versus pressure head, `(alpha * |h_inflection|)^n = m`; this removes alpha
/// from the remaining scalar solve for n.
pub fn fitOriginalMualemVanGenuchten(
    inputs: MualemVanGenuchtenFitInputs,
    options: MualemVanGenuchtenFitOptions,
) !MualemVanGenuchtenFitResult {
    try options.validate();
    inline for (@typeInfo(MualemVanGenuchtenFitInputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteMualemVanGenuchtenFitInput;
    }
    if (inputs.saturated_water_content_m3_per_m3 <= 0 or
        inputs.saturated_water_content_m3_per_m3 > 1 or
        inputs.field_capacity_water_content_m3_per_m3 <= 0 or
        inputs.field_capacity_water_content_m3_per_m3 >= inputs.saturated_water_content_m3_per_m3 or
        inputs.wilting_point_water_content_m3_per_m3 <= 0 or
        inputs.wilting_point_water_content_m3_per_m3 >= inputs.field_capacity_water_content_m3_per_m3 or
        inputs.field_capacity_pressure_head_m >= 0 or
        inputs.wilting_point_pressure_head_m >= inputs.field_capacity_pressure_head_m or
        inputs.inflection_pressure_head_m >= 0 or
        inputs.saturated_hydraulic_conductivity_m_per_h < 0)
    {
        return error.InvalidMualemVanGenuchtenFitInput;
    }

    var lower_n = options.minimum_n;
    var upper_n = options.maximum_n;
    var lower = try evaluateRetentionFit(inputs, lower_n);
    var upper = try evaluateRetentionFit(inputs, upper_n);
    if (lower.residual == 0)
        return makeRetentionFitResult(inputs, lower_n, lower, 1, 0, 0);
    if (upper.residual == 0)
        return makeRetentionFitResult(inputs, upper_n, upper, 1, 0, 0);
    if (std.math.signbit(lower.residual) == std.math.signbit(upper.residual))
        return error.MualemVanGenuchtenFitNotBracketed;

    var n = 0.5 * (lower_n + upper_n);
    var evaluation = try evaluateRetentionFit(inputs, n);
    var newton_steps: u16 = 0;
    var bisection_steps: u16 = 0;
    var iteration: u16 = 0;
    while (iteration < options.maximum_iterations) : (iteration += 1) {
        if (@abs(evaluation.residual) <= options.water_content_tolerance_m3_per_m3)
            return makeRetentionFitResult(
                inputs,
                n,
                evaluation,
                iteration + 1,
                newton_steps,
                bisection_steps,
            );
        if (std.math.signbit(evaluation.residual) == std.math.signbit(lower.residual)) {
            lower_n = n;
            lower = evaluation;
        } else {
            upper_n = n;
            upper = evaluation;
        }

        const derivative_probe =
            std.math.sqrt(std.math.floatEps(f64)) * @max(1.0, @abs(n));
        const probe_n = @min(upper_n, n + derivative_probe);
        var next_n = 0.5 * (lower_n + upper_n);
        var used_newton = false;
        if (probe_n > n) {
            const probe = try evaluateRetentionFit(inputs, probe_n);
            const derivative = (probe.residual - evaluation.residual) / (probe_n - n);
            if (std.math.isFinite(derivative) and derivative != 0) {
                const newton_n = n - evaluation.residual / derivative;
                if (newton_n > lower_n and newton_n < upper_n and std.math.isFinite(newton_n)) {
                    next_n = newton_n;
                    used_newton = true;
                }
            }
        }
        if (used_newton) {
            newton_steps += 1;
        } else {
            bisection_steps += 1;
        }
        n = next_n;
        evaluation = try evaluateRetentionFit(inputs, n);
    }
    return error.MualemVanGenuchtenFitDidNotConverge;
}

fn evaluateRetentionFit(
    inputs: MualemVanGenuchtenFitInputs,
    n: f64,
) !RetentionFitEvaluation {
    const m = 1.0 - 1.0 / n;
    const alpha_per_m =
        std.math.pow(f64, m, 1.0 / n) / -inputs.inflection_pressure_head_m;
    const field_effective_saturation = std.math.pow(
        f64,
        1.0 + std.math.pow(
            f64,
            alpha_per_m * -inputs.field_capacity_pressure_head_m,
            n,
        ),
        -m,
    );
    const wilting_effective_saturation = std.math.pow(
        f64,
        1.0 + std.math.pow(
            f64,
            alpha_per_m * -inputs.wilting_point_pressure_head_m,
            n,
        ),
        -m,
    );
    const field_denominator = 1.0 - field_effective_saturation;
    if (!std.math.isFinite(alpha_per_m) or alpha_per_m <= 0 or
        !std.math.isFinite(field_denominator) or field_denominator <= 0)
    {
        return error.InvalidMualemVanGenuchtenFitCandidate;
    }
    const residual_water_content =
        (inputs.field_capacity_water_content_m3_per_m3 -
            inputs.saturated_water_content_m3_per_m3 *
                field_effective_saturation) /
        field_denominator;
    if (!std.math.isFinite(residual_water_content))
        return error.InvalidMualemVanGenuchtenFitCandidate;
    const predicted_wilting_water_content =
        residual_water_content +
        (inputs.saturated_water_content_m3_per_m3 - residual_water_content) *
            wilting_effective_saturation;
    return .{
        .residual = predicted_wilting_water_content -
            inputs.wilting_point_water_content_m3_per_m3,
        .residual_water_content_m3_per_m3 = residual_water_content,
        .alpha_per_m = alpha_per_m,
    };
}

fn makeRetentionFitResult(
    inputs: MualemVanGenuchtenFitInputs,
    n: f64,
    evaluation: RetentionFitEvaluation,
    iterations: u16,
    newton_steps: u16,
    bisection_steps: u16,
) !MualemVanGenuchtenFitResult {
    const parameters: MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = evaluation.residual_water_content_m3_per_m3,
        .saturated_water_content_m3_per_m3 = inputs.saturated_water_content_m3_per_m3,
        .alpha_per_m = evaluation.alpha_per_m,
        .n = n,
        .pore_connectivity = inputs.pore_connectivity,
        .saturated_hydraulic_conductivity_m_per_h = inputs.saturated_hydraulic_conductivity_m_per_h,
    };
    try parameters.validate();
    return .{
        .parameters = parameters,
        .iterations = iterations,
        .newton_raphson_steps = newton_steps,
        .bisection_steps = bisection_steps,
        .residual_water_content_m3_per_m3 = evaluation.residual_water_content_m3_per_m3,
    };
}

/// Runtime science parameters for the HOUR1/WATSUB water-retention functions.
/// Values are loaded by the caller; no parameters.h-sized global state exists.
pub const Parameters = struct {
    saturation_water_potential_mpa: f64,
    minimum_water_potential_mpa: f64,
    organic_soil_threshold_g_per_megagram: f64,
    mineral_field_capacity_intercept: f64,
    mineral_field_capacity_sand_coefficient: f64,
    mineral_field_capacity_clay_coefficient: f64,
    mineral_field_capacity_organic_coefficient_per_g_per_megagram: f64,
    mineral_wilting_point_intercept: f64,
    mineral_wilting_point_clay_coefficient: f64,
    mineral_wilting_point_organic_coefficient_per_g_per_megagram: f64,
    organic_bulk_density_threshold_1_megagrams_per_m3: f64,
    organic_bulk_density_threshold_2_megagrams_per_m3: f64,
    organic_field_capacity_1: f64,
    organic_field_capacity_2: f64,
    organic_field_capacity_3: f64,
    organic_wilting_point_1: f64,
    organic_wilting_point_2: f64,
    organic_wilting_point_3: f64,
    maximum_field_capacity_fraction_of_porosity: f64,
    maximum_wilting_point_fraction_of_field_capacity: f64,
    saturation_to_field_shape: f64,
    below_wilting_shape: f64,

    pub fn validate(self: Parameters) !void {
        inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(self, field.name))) return error.NonFiniteSoilRetentionParameter;
        if (self.saturation_water_potential_mpa >= 0 or self.minimum_water_potential_mpa >= self.saturation_water_potential_mpa or self.organic_soil_threshold_g_per_megagram < 0 or self.organic_bulk_density_threshold_1_megagrams_per_m3 <= 0 or self.organic_bulk_density_threshold_2_megagrams_per_m3 <= self.organic_bulk_density_threshold_1_megagrams_per_m3 or self.maximum_field_capacity_fraction_of_porosity <= 0 or self.maximum_field_capacity_fraction_of_porosity > 1 or self.maximum_wilting_point_fraction_of_field_capacity <= 0 or self.maximum_wilting_point_fraction_of_field_capacity > 1 or self.saturation_to_field_shape <= 0 or self.below_wilting_shape <= 0) return error.InvalidSoilRetentionParameter;
    }
};

/// HOUR1 coefficients used only when an older runscript omits the runtime
/// soil_solver record. The returned value is ordinary runtime data.
pub fn compatibilityParameters() Parameters {
    return .{
        .saturation_water_potential_mpa = -0.0005,
        .minimum_water_potential_mpa = -1.5e12,
        .organic_soil_threshold_g_per_megagram = 250_000,
        .mineral_field_capacity_intercept = 0.2576,
        .mineral_field_capacity_sand_coefficient = -0.20,
        .mineral_field_capacity_clay_coefficient = 0.36,
        .mineral_field_capacity_organic_coefficient_per_g_per_megagram = 0.60e-6,
        .mineral_wilting_point_intercept = 0.0260,
        .mineral_wilting_point_clay_coefficient = 0.50,
        .mineral_wilting_point_organic_coefficient_per_g_per_megagram = 0.32e-6,
        .organic_bulk_density_threshold_1_megagrams_per_m3 = 0.075,
        .organic_bulk_density_threshold_2_megagrams_per_m3 = 0.195,
        .organic_field_capacity_1 = 0.27,
        .organic_field_capacity_2 = 0.62,
        .organic_field_capacity_3 = 0.71,
        .organic_wilting_point_1 = 0.04,
        .organic_wilting_point_2 = 0.15,
        .organic_wilting_point_3 = 0.22,
        .maximum_field_capacity_fraction_of_porosity = 0.75,
        .maximum_wilting_point_fraction_of_field_capacity = 0.75,
        .saturation_to_field_shape = 0.5,
        .below_wilting_shape = 0.5,
    };
}

pub const LayerInputs = struct {
    porosity_fraction: f64,
    macropore_fraction: f64,
    sand_fraction: f64,
    clay_fraction: f64,
    organic_carbon_g_per_megagram: f64,
    bulk_density_megagrams_per_m3: f64,
    supplied_field_capacity_fraction: ?f64,
    supplied_wilting_point_fraction: ?f64,
};

pub const Curve = struct {
    field_capacity_fraction: f64,
    wilting_point_fraction: f64,
    saturation_water_potential_mpa: f64,
    field_capacity_water_potential_mpa: f64,
    wilting_point_water_potential_mpa: f64,
    minimum_water_potential_mpa: f64,
    saturation_to_field_shape: f64,
    below_wilting_shape: f64,
};

pub const ResolvedCurve = struct {
    porosity_fraction: f64,
    curve: Curve,

    pub fn waterPotentialMpa(self: ResolvedCurve, water_fraction: f64) !f64 {
        if (!std.math.isFinite(water_fraction) or water_fraction <= 0) return error.InvalidSoilWaterFraction;
        const water = @min(water_fraction, self.porosity_fraction);
        if (water >= self.porosity_fraction) return self.curve.saturation_water_potential_mpa;
        const log_water = @log(water);
        const log_porosity = @log(self.porosity_fraction);
        const log_field_capacity = @log(self.curve.field_capacity_fraction);
        const log_wilting_point = @log(self.curve.wilting_point_fraction);
        const log_saturation_potential = @log(-self.curve.saturation_water_potential_mpa);
        const log_field_potential = @log(-self.curve.field_capacity_water_potential_mpa);
        const log_wilting_potential = @log(-self.curve.wilting_point_water_potential_mpa);
        const potential = if (water < self.curve.wilting_point_fraction)
            -@exp(log_wilting_potential + self.curve.below_wilting_shape * ((log_wilting_point - log_water) / (log_field_capacity - log_wilting_point) * (log_wilting_potential - log_field_potential)))
        else if (water < self.curve.field_capacity_fraction)
            -@exp(log_field_potential + ((log_field_capacity - log_water) / (log_field_capacity - log_wilting_point) * (log_wilting_potential - log_field_potential)))
        else
            -@exp(log_saturation_potential + std.math.pow(f64, @max(0.0, (log_porosity - log_water) / (log_porosity - log_field_capacity)), self.curve.saturation_to_field_shape) * (log_field_potential - log_saturation_potential));
        return @max(self.curve.minimum_water_potential_mpa, potential);
    }

    /// Exact inverse of the HOUR1 log-water branches. NITRO evaluates this at
    /// PSIHY=-1.5e4 MPa to exclude hygroscopic water from active water.
    pub fn waterFractionAtPotentialMpa(self: ResolvedCurve, target_potential_mpa: f64) !f64 {
        if (!std.math.isFinite(target_potential_mpa) or target_potential_mpa >= 0) return error.InvalidTargetWaterPotential;
        const c = self.curve;
        const log_target = @log(-target_potential_mpa);
        const log_saturation_potential = @log(-c.saturation_water_potential_mpa);
        const log_field_potential = @log(-c.field_capacity_water_potential_mpa);
        const log_wilting_potential = @log(-c.wilting_point_water_potential_mpa);
        const log_porosity = @log(self.porosity_fraction);
        const log_field = @log(c.field_capacity_fraction);
        const log_wilting = @log(c.wilting_point_fraction);
        const log_water = if (target_potential_mpa < c.wilting_point_water_potential_mpa)
            log_wilting - (log_target - log_wilting_potential) * (log_field - log_wilting) / (c.below_wilting_shape * (log_wilting_potential - log_field_potential))
        else if (target_potential_mpa < c.field_capacity_water_potential_mpa)
            log_field - (log_target - log_field_potential) * (log_field - log_wilting) / (log_wilting_potential - log_field_potential)
        else
            log_porosity - std.math.pow(f64, std.math.clamp((log_target - log_saturation_potential) / (log_field_potential - log_saturation_potential), 0, 1), 1 / c.saturation_to_field_shape) * (log_porosity - log_field);
        const water = std.math.clamp(@exp(log_water), 0, self.porosity_fraction);
        if (!std.math.isFinite(water)) return error.NonFiniteWaterFraction;
        return water;
    }
};

pub fn resolve(parameters: Parameters, inputs: LayerInputs, field_capacity_water_potential_mpa: f64, wilting_point_water_potential_mpa: f64) !ResolvedCurve {
    try parameters.validate();
    inline for (.{ inputs.porosity_fraction, inputs.macropore_fraction, inputs.sand_fraction, inputs.clay_fraction, inputs.organic_carbon_g_per_megagram, inputs.bulk_density_megagrams_per_m3, field_capacity_water_potential_mpa, wilting_point_water_potential_mpa }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSoilRetentionInput;
    if (inputs.porosity_fraction <= 0 or inputs.porosity_fraction > 1 or inputs.macropore_fraction < 0 or inputs.macropore_fraction >= 1 or inputs.sand_fraction < 0 or inputs.clay_fraction < 0 or inputs.sand_fraction + inputs.clay_fraction > 1 or inputs.organic_carbon_g_per_megagram < 0 or inputs.bulk_density_megagrams_per_m3 < 0 or field_capacity_water_potential_mpa >= 0 or wilting_point_water_potential_mpa >= field_capacity_water_potential_mpa) return error.InvalidSoilRetentionInput;

    const both_supplied = inputs.supplied_field_capacity_fraction != null and inputs.supplied_wilting_point_fraction != null;
    var field_capacity = if (both_supplied) inputs.supplied_field_capacity_fraction.? else estimateFieldCapacity(parameters, inputs);
    var wilting_point = if (both_supplied) inputs.supplied_wilting_point_fraction.? else estimateWiltingPoint(parameters, inputs);
    // HOUR1 estimates both properties when either ISOIL flag says unknown;
    // only that estimated branch applies the non-macropore correction and
    // the 0.75*POROS / 0.75*FC ceilings.
    if (!both_supplied) {
        field_capacity = @min(parameters.maximum_field_capacity_fraction_of_porosity * inputs.porosity_fraction, field_capacity / (1.0 - inputs.macropore_fraction));
        wilting_point = @min(parameters.maximum_wilting_point_fraction_of_field_capacity * field_capacity, wilting_point / (1.0 - inputs.macropore_fraction));
    }
    if (!std.math.isFinite(field_capacity) or !std.math.isFinite(wilting_point) or wilting_point <= 0 or field_capacity <= wilting_point or field_capacity >= inputs.porosity_fraction) return error.InvalidResolvedSoilRetentionCurve;
    return .{ .porosity_fraction = inputs.porosity_fraction, .curve = .{ .field_capacity_fraction = field_capacity, .wilting_point_fraction = wilting_point, .saturation_water_potential_mpa = parameters.saturation_water_potential_mpa, .field_capacity_water_potential_mpa = field_capacity_water_potential_mpa, .wilting_point_water_potential_mpa = wilting_point_water_potential_mpa, .minimum_water_potential_mpa = parameters.minimum_water_potential_mpa, .saturation_to_field_shape = parameters.saturation_to_field_shape, .below_wilting_shape = parameters.below_wilting_shape } };
}

fn estimateFieldCapacity(parameters: Parameters, inputs: LayerInputs) f64 {
    if (inputs.organic_carbon_g_per_megagram < parameters.organic_soil_threshold_g_per_megagram) return parameters.mineral_field_capacity_intercept + parameters.mineral_field_capacity_sand_coefficient * inputs.sand_fraction + parameters.mineral_field_capacity_clay_coefficient * inputs.clay_fraction + parameters.mineral_field_capacity_organic_coefficient_per_g_per_megagram * inputs.organic_carbon_g_per_megagram;
    if (inputs.bulk_density_megagrams_per_m3 < parameters.organic_bulk_density_threshold_1_megagrams_per_m3) return parameters.organic_field_capacity_1;
    if (inputs.bulk_density_megagrams_per_m3 < parameters.organic_bulk_density_threshold_2_megagrams_per_m3) return parameters.organic_field_capacity_2;
    return parameters.organic_field_capacity_3;
}

fn estimateWiltingPoint(parameters: Parameters, inputs: LayerInputs) f64 {
    if (inputs.organic_carbon_g_per_megagram < parameters.organic_soil_threshold_g_per_megagram) return parameters.mineral_wilting_point_intercept + parameters.mineral_wilting_point_clay_coefficient * inputs.clay_fraction + parameters.mineral_wilting_point_organic_coefficient_per_g_per_megagram * inputs.organic_carbon_g_per_megagram;
    if (inputs.bulk_density_megagrams_per_m3 < parameters.organic_bulk_density_threshold_1_megagrams_per_m3) return parameters.organic_wilting_point_1;
    if (inputs.bulk_density_megagrams_per_m3 < parameters.organic_bulk_density_threshold_2_megagrams_per_m3) return parameters.organic_wilting_point_2;
    return parameters.organic_wilting_point_3;
}

test "original Mualem van Genuchten retention is invertible and monotone" {
    const loam: MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.078,
        .saturated_water_content_m3_per_m3 = 0.43,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.0104,
    };
    try loam.validate();
    try std.testing.expectApproxEqAbs(@as(f64, 0.43), try loam.waterContentAtPressureHead(0), 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0104), try loam.hydraulicConductivityMPerH(0), 1.0e-15);
    const pressure_heads_m = [_]f64{ -0.1, -0.33, -1.5, -15.0 };
    var previous_water_content = loam.saturated_water_content_m3_per_m3;
    var previous_conductivity = loam.saturated_hydraulic_conductivity_m_per_h;
    for (pressure_heads_m) |pressure_head_m| {
        const water_content = try loam.waterContentAtPressureHead(pressure_head_m);
        const reconstructed_head =
            try loam.pressureHeadAtWaterContent(water_content);
        const conductivity =
            try loam.hydraulicConductivityMPerH(pressure_head_m);
        try std.testing.expectApproxEqRel(pressure_head_m, reconstructed_head, 1.0e-12);
        try std.testing.expect(water_content < previous_water_content);
        try std.testing.expect(conductivity < previous_conductivity);
        try std.testing.expect(try loam.waterCapacityPerM(pressure_head_m) > 0);
        previous_water_content = water_content;
        previous_conductivity = conductivity;
    }
}

test "original van Genuchten residual endpoint has finite asymptotic head" {
    const parameters: MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.45,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.01,
    };
    const pressure_head_m = try parameters.pressureHeadAtWaterContent(
        parameters.residual_water_content_m3_per_m3,
    );
    try std.testing.expect(std.math.isFinite(pressure_head_m));
    try std.testing.expect(pressure_head_m < 0);
    const recovered = try parameters.waterContentAtPressureHead(pressure_head_m);
    try std.testing.expectApproxEqAbs(
        parameters.residual_water_content_m3_per_m3,
        recovered,
        std.math.floatEps(f64),
    );
}

test "original Mualem van Genuchten excludes nonphysical parameters" {
    var parameters: MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.1,
        .saturated_water_content_m3_per_m3 = 0.5,
        .alpha_per_m = 2,
        .n = 1.5,
        .saturated_hydraulic_conductivity_m_per_h = 0.01,
    };
    parameters.n = 1;
    try std.testing.expectError(error.InvalidMualemVanGenuchtenParameter, parameters.validate());
    parameters.n = 1.5;
    parameters.saturated_water_content_m3_per_m3 = 0.1;
    try std.testing.expectError(error.InvalidMualemVanGenuchtenParameter, parameters.validate());
}

test "runtime inflection fit recovers original Mualem van Genuchten parameters" {
    const expected: MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.078,
        .saturated_water_content_m3_per_m3 = 0.43,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.0104,
    };
    const field_capacity_pressure_head_m = -0.33;
    const wilting_point_pressure_head_m = -15;
    const inflection_pressure_head_m =
        -std.math.pow(f64, expected.m(), 1.0 / expected.n) /
        expected.alpha_per_m;
    const fitted = try fitOriginalMualemVanGenuchten(.{
        .saturated_water_content_m3_per_m3 = expected.saturated_water_content_m3_per_m3,
        .field_capacity_water_content_m3_per_m3 = try expected.waterContentAtPressureHead(field_capacity_pressure_head_m),
        .field_capacity_pressure_head_m = field_capacity_pressure_head_m,
        .wilting_point_water_content_m3_per_m3 = try expected.waterContentAtPressureHead(wilting_point_pressure_head_m),
        .wilting_point_pressure_head_m = wilting_point_pressure_head_m,
        .inflection_pressure_head_m = inflection_pressure_head_m,
        .saturated_hydraulic_conductivity_m_per_h = expected.saturated_hydraulic_conductivity_m_per_h,
    }, .{ .maximum_iterations = 60 });
    try std.testing.expect(fitted.iterations < 60);
    try std.testing.expectApproxEqRel(expected.n, fitted.parameters.n, 1.0e-8);
    try std.testing.expectApproxEqRel(expected.alpha_per_m, fitted.parameters.alpha_per_m, 1.0e-8);
    try std.testing.expectApproxEqAbs(
        expected.residual_water_content_m3_per_m3,
        fitted.parameters.residual_water_content_m3_per_m3,
        1.0e-9,
    );
}

test "runtime Carsel Parrish fallback follows USDA texture and supplied porosity" {
    const texture = try classifyUsdaSoilTexture(0.40, 0.40, 0.20);
    try std.testing.expectEqual(SoilTextureClass.loam, texture);
    const loam = try carselParrishDefault(texture, 0.47);
    try std.testing.expectApproxEqAbs(@as(f64, 0.078), loam.residual_water_content_m3_per_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.47), loam.saturated_water_content_m3_per_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3.6), loam.alpha_per_m, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.56), loam.n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0104), loam.saturated_hydraulic_conductivity_m_per_h, 1.0e-15);
}

test "runtime retention parameters reproduce mineral threshold equations" {
    const parameters: Parameters = .{ .saturation_water_potential_mpa = -0.0005, .minimum_water_potential_mpa = -1.5e12, .organic_soil_threshold_g_per_megagram = 250_000, .mineral_field_capacity_intercept = 0.2576, .mineral_field_capacity_sand_coefficient = -0.20, .mineral_field_capacity_clay_coefficient = 0.36, .mineral_field_capacity_organic_coefficient_per_g_per_megagram = 0.60e-6, .mineral_wilting_point_intercept = 0.0260, .mineral_wilting_point_clay_coefficient = 0.50, .mineral_wilting_point_organic_coefficient_per_g_per_megagram = 0.32e-6, .organic_bulk_density_threshold_1_megagrams_per_m3 = 0.075, .organic_bulk_density_threshold_2_megagrams_per_m3 = 0.195, .organic_field_capacity_1 = 0.27, .organic_field_capacity_2 = 0.62, .organic_field_capacity_3 = 0.71, .organic_wilting_point_1 = 0.04, .organic_wilting_point_2 = 0.15, .organic_wilting_point_3 = 0.22, .maximum_field_capacity_fraction_of_porosity = 0.75, .maximum_wilting_point_fraction_of_field_capacity = 0.75, .saturation_to_field_shape = 0.5, .below_wilting_shape = 0.5 };
    const resolved = try resolve(parameters, .{ .porosity_fraction = 0.5, .macropore_fraction = 0, .sand_fraction = 0.4, .clay_fraction = 0.2, .organic_carbon_g_per_megagram = 10_000, .bulk_density_megagrams_per_m3 = 1.3, .supplied_field_capacity_fraction = null, .supplied_wilting_point_fraction = null }, -0.01, -1.5);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2556), resolved.curve.field_capacity_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1292), resolved.curve.wilting_point_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -0.01), try resolved.waterPotentialMpa(resolved.curve.field_capacity_fraction), 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -1.5), try resolved.waterPotentialMpa(resolved.curve.wilting_point_fraction), 1.0e-12);
    const hygroscopic = try resolved.waterFractionAtPotentialMpa(-1.5e4);
    try std.testing.expectApproxEqRel(@as(f64, -1.5e4), try resolved.waterPotentialMpa(hygroscopic), 1.0e-12);
}

test "one missing hydraulic property makes HOUR1 estimate both properties" {
    const parameters = compatibilityParameters();
    const inputs: LayerInputs = .{ .porosity_fraction = 0.5, .macropore_fraction = 0.05, .sand_fraction = 0.4, .clay_fraction = 0.2, .organic_carbon_g_per_megagram = 10_000, .bulk_density_megagrams_per_m3 = 1.3, .supplied_field_capacity_fraction = 0.4, .supplied_wilting_point_fraction = null };
    const resolved = try resolve(parameters, inputs, -0.033, -1.5);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2556 / 0.95), resolved.curve.field_capacity_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1292 / 0.95), resolved.curve.wilting_point_fraction, 1.0e-12);
}
