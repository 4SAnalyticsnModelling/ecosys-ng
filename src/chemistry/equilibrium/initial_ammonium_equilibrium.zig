const std = @import("std");

pub const Source = enum(u8) {
    precipitation = 1,
    irrigation = 2,
    soil = 3,
};

pub const Inputs = struct {
    source: Source,
    total_ammoniacal_n_mol_per_source_volume: f64, // CN4Z
    hydrogen_mol_per_m3: f64, // CHY1
    ammonium_dissociation_constant_mol_per_m3: f64, // DPN4
    soil_ammonium_seed_fraction: f64, // source 0.10
};

pub const Result = struct {
    ammonium_mol_n_per_source_volume: f64, // CN41
    ammonia_mol_n_per_source_volume: f64, // CN31
};

/// Direct translation of `starte.f` lines 258--264. Precipitation and
/// irrigation concentrations are in `mol m-3`; soil inputs retain the
/// current soil source-volume basis established by the enclosing branch.
pub fn calculate(inputs: Inputs) !Result {
    if (!std.math.isFinite(inputs.total_ammoniacal_n_mol_per_source_volume))
        return error.NonFiniteInitialAmmoniumInput;

    const result: Result = switch (inputs.source) {
        .precipitation, .irrigation => blk: {
            if (!std.math.isFinite(inputs.hydrogen_mol_per_m3) or
                !std.math.isFinite(inputs.ammonium_dissociation_constant_mol_per_m3))
                return error.NonFiniteInitialAmmoniumInput;
            if (inputs.hydrogen_mol_per_m3 <= 0 or
                inputs.ammonium_dissociation_constant_mol_per_m3 <= 0)
                return error.InvalidInitialAmmoniumInput;
            const ammonium = inputs.total_ammoniacal_n_mol_per_source_volume /
                (1.0 + inputs.ammonium_dissociation_constant_mol_per_m3 /
                    inputs.hydrogen_mol_per_m3);
            break :blk .{
                .ammonium_mol_n_per_source_volume = ammonium,
                .ammonia_mol_n_per_source_volume = ammonium *
                    inputs.ammonium_dissociation_constant_mol_per_m3 /
                    inputs.hydrogen_mol_per_m3,
            };
        },
        .soil => blk: {
            if (!std.math.isFinite(inputs.soil_ammonium_seed_fraction))
                return error.NonFiniteInitialAmmoniumInput;
            if (inputs.soil_ammonium_seed_fraction < 0)
                return error.InvalidInitialAmmoniumInput;
            break :blk .{
                .ammonium_mol_n_per_source_volume = inputs.total_ammoniacal_n_mol_per_source_volume *
                    inputs.soil_ammonium_seed_fraction,
                .ammonia_mol_n_per_source_volume = 0.0,
            };
        },
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteInitialAmmoniumResult;
    }
    return result;
}

fn fixtureInputs(source: Source) Inputs {
    return .{
        .source = source,
        .total_ammoniacal_n_mol_per_source_volume = 12,
        .hydrogen_mol_per_m3 = 2,
        .ammonium_dissociation_constant_mol_per_m3 = 4,
        .soil_ammonium_seed_fraction = 0.10,
    };
}

test "STARTE precipitation and irrigation share exact ammonium equilibrium order" {
    const precipitation = try calculate(fixtureInputs(.precipitation));
    const irrigation = try calculate(fixtureInputs(.irrigation));
    try std.testing.expectEqualDeep(precipitation, irrigation);
    try std.testing.expectEqual(@as(f64, 4), precipitation.ammonium_mol_n_per_source_volume);
    try std.testing.expectEqual(@as(f64, 8), precipitation.ammonia_mol_n_per_source_volume);
    try std.testing.expectEqual(
        @as(f64, 12),
        precipitation.ammonium_mol_n_per_source_volume +
            precipitation.ammonia_mol_n_per_source_volume,
    );
}

test "STARTE soil ammonium seed preserves source fraction and zero ammonia" {
    const result = try calculate(fixtureInputs(.soil));
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.2),
        result.ammonium_mol_n_per_source_volume,
        1.0e-15,
    );
    try std.testing.expectEqual(@as(f64, 0), result.ammonia_mol_n_per_source_volume);
}

test "STARTE ammonium calculation validates before publishing a result" {
    var invalid = fixtureInputs(.irrigation);
    invalid.ammonium_dissociation_constant_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteInitialAmmoniumInput, calculate(invalid));

    invalid = fixtureInputs(.irrigation);
    invalid.hydrogen_mol_per_m3 = 0;
    try std.testing.expectError(error.InvalidInitialAmmoniumInput, calculate(invalid));
}

test "STARTE ammonium branches do not inspect dormant parameters" {
    var soil = fixtureInputs(.soil);
    soil.hydrogen_mol_per_m3 = std.math.nan(f64);
    soil.ammonium_dissociation_constant_mol_per_m3 = std.math.nan(f64);
    _ = try calculate(soil);

    var rain = fixtureInputs(.precipitation);
    rain.soil_ammonium_seed_fraction = std.math.nan(f64);
    _ = try calculate(rain);
}
