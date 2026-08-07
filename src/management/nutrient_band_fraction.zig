const std = @import("std");

pub const PartitionInput = struct {
    non_band_mass_g: f64,
    band_mass_g: f64,
    fallback_non_band_fraction: f64,
    fallback_band_fraction: f64,
};

pub const Inputs = struct {
    ammonium: PartitionInput,
    nitrate: PartitionInput,
    dihydrogen_phosphate: PartitionInput,
    hydrogen_phosphate: PartitionInput,
    mass_threshold_g: f64,
};

pub const Fractions = struct {
    non_band: f64,
    band: f64,
};

pub const Result = struct {
    ammonium: Fractions,
    nitrate: Fractions,
    dihydrogen_phosphate: Fractions,
    hydrogen_phosphate: Fractions,
};

pub const PartitionError = error{
    NonFiniteInput,
    InvalidMassThreshold,
    InvalidFallbackFraction,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 3296--3335. Nutrient masses are g element and
/// returned partitions are dimensionless.
pub fn calculate(inputs: Inputs) PartitionError!Result {
    if (!std.math.isFinite(inputs.mass_threshold_g)) return error.NonFiniteInput;
    if (inputs.mass_threshold_g < 0.0) return error.InvalidMassThreshold;
    inline for ([_][]const u8{
        "ammonium",
        "nitrate",
        "dihydrogen_phosphate",
        "hydrogen_phosphate",
    }) |name| {
        const partition_input = @field(inputs, name);
        inline for (std.meta.fields(PartitionInput)) |field| {
            if (!std.math.isFinite(@field(partition_input, field.name))) {
                return error.NonFiniteInput;
            }
        }
        if (partition_input.fallback_non_band_fraction < 0.0 or
            partition_input.fallback_non_band_fraction > 1.0 or
            partition_input.fallback_band_fraction < 0.0 or
            partition_input.fallback_band_fraction > 1.0)
        {
            return error.InvalidFallbackFraction;
        }
    }

    const result = Result{
        .ammonium = partition(inputs.ammonium, inputs.mass_threshold_g),
        .nitrate = partition(inputs.nitrate, inputs.mass_threshold_g),
        .dihydrogen_phosphate = partition(inputs.dihydrogen_phosphate, inputs.mass_threshold_g),
        .hydrogen_phosphate = partition(inputs.hydrogen_phosphate, inputs.mass_threshold_g),
    };
    inline for (std.meta.fields(Result)) |field| {
        const fractions = @field(result, field.name);
        if (!std.math.isFinite(fractions.non_band) or !std.math.isFinite(fractions.band)) {
            return error.NonFiniteResult;
        }
    }
    return result;
}

fn partition(input: PartitionInput, mass_threshold_g: f64) Fractions {
    const non_band_mass_g = @max(0.0, input.non_band_mass_g);
    const band_mass_g = @max(0.0, input.band_mass_g);
    const total_mass_g = non_band_mass_g + band_mass_g;
    if (total_mass_g > mass_threshold_g) {
        return .{
            .non_band = non_band_mass_g / total_mass_g,
            .band = band_mass_g / total_mass_g,
        };
    }
    return .{
        .non_band = input.fallback_non_band_fraction,
        .band = input.fallback_band_fraction,
    };
}

test "positive nutrient masses produce normalized band partitions" {
    const common = PartitionInput{
        .non_band_mass_g = 3.0,
        .band_mass_g = 1.0,
        .fallback_non_band_fraction = 0.2,
        .fallback_band_fraction = 0.8,
    };
    const result = try calculate(.{
        .ammonium = common,
        .nitrate = common,
        .dihydrogen_phosphate = common,
        .hydrogen_phosphate = common,
        .mass_threshold_g = 1.0e-12,
    });
    inline for (std.meta.fields(Result)) |field| {
        try std.testing.expectEqual(@as(f64, 0.75), @field(result, field.name).non_band);
        try std.testing.expectEqual(@as(f64, 0.25), @field(result, field.name).band);
    }
}

test "negative masses clamp to zero and low totals use fallback" {
    const low = PartitionInput{
        .non_band_mass_g = -1.0,
        .band_mass_g = 0.0,
        .fallback_non_band_fraction = 0.6,
        .fallback_band_fraction = 0.4,
    };
    const result = try calculate(.{
        .ammonium = low,
        .nitrate = low,
        .dihydrogen_phosphate = low,
        .hydrogen_phosphate = low,
        .mass_threshold_g = 0.0,
    });
    try std.testing.expectEqual(@as(f64, 0.6), result.ammonium.non_band);
    try std.testing.expectEqual(@as(f64, 0.4), result.ammonium.band);
}
