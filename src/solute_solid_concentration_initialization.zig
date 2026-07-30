const std = @import("std");
const geochemistry = @import("solute_geochemistry_network.zig");

pub const PrecipitateAmounts = struct {
    gibbsite_mol: f64,
    iron_hydroxide_mol: f64,
    calcite_mol: f64,
    gypsum_mol: f64,
};

pub const SilicateAmounts = struct {
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
};

pub const Inputs = struct {
    precipitates: PrecipitateAmounts,
    natural_silicates: SilicateAmounts,
    ground_silicates: SilicateAmounts,
    water_volume_m3: f64,
    minimum_water_volume_m3: f64,
    ground_silicate_specific_surface_area_m2_per_mol: f64,
};

pub const Result = struct {
    ground_silicate_surface_area_m2_per_m3: f64,
};

/// Direct translation of SOLUTE lines 779--818. The source floors solid
/// inventories at zero, unlike aqueous `ZEROC` concentration floors.
pub fn applyConcentrations(
    state: *geochemistry.SolidState,
    inputs: Inputs,
) !Result {
    try validateInputs(inputs);
    var staged = state.*;
    if (inputs.water_volume_m3 > inputs.minimum_water_volume_m3) {
        const volume = inputs.water_volume_m3;
        const precipitates = inputs.precipitates;
        staged.gibbsite_solid_mol_per_m3 =
            nonnegativeConcentration(precipitates.gibbsite_mol, volume);
        staged.iron_hydroxide_solid_mol_per_m3 =
            nonnegativeConcentration(precipitates.iron_hydroxide_mol, volume);
        staged.calcite_solid_mol_per_m3 =
            nonnegativeConcentration(precipitates.calcite_mol, volume);
        staged.gypsum_solid_mol_per_m3 =
            nonnegativeConcentration(precipitates.gypsum_mol, volume);
        applySilicateConcentrations(
            &staged,
            inputs.natural_silicates,
            volume,
            .natural,
        );
        applySilicateConcentrations(
            &staged,
            inputs.ground_silicates,
            volume,
            .ground,
        );
    } else {
        inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field|
            @field(staged, field.name) = 0;
    }
    try validateState(staged);
    const ground_total_mol_per_m3 =
        staged.aluminum_ground_silicate_mol_per_m3 +
        staged.iron_ground_silicate_mol_per_m3 +
        staged.calcium_ground_silicate_mol_per_m3 +
        staged.magnesium_ground_silicate_mol_per_m3 +
        staged.sodium_ground_silicate_mol_per_m3 +
        staged.potassium_ground_silicate_mol_per_m3;
    const surface_area_m2_per_m3 =
        inputs.ground_silicate_specific_surface_area_m2_per_mol *
        ground_total_mol_per_m3;
    if (!std.math.isFinite(surface_area_m2_per_m3))
        return error.NonFiniteGroundSilicateSurfaceArea;
    state.* = staged;
    return .{
        .ground_silicate_surface_area_m2_per_m3 = surface_area_m2_per_m3,
    };
}

const SilicateDomain = enum { natural, ground };

fn applySilicateConcentrations(
    state: *geochemistry.SolidState,
    amounts: SilicateAmounts,
    water_volume_m3: f64,
    comptime domain: SilicateDomain,
) void {
    inline for (.{
        .{ "aluminum", "aluminum_mol" },
        .{ "iron", "iron_mol" },
        .{ "calcium", "calcium_mol" },
        .{ "magnesium", "magnesium_mol" },
        .{ "sodium", "sodium_mol" },
        .{ "potassium", "potassium_mol" },
    }) |mapping| {
        const field_name = mapping[0] ++ "_" ++
            @tagName(domain) ++ "_silicate_mol_per_m3";
        @field(state, field_name) = nonnegativeConcentration(
            @field(amounts, mapping[1]),
            water_volume_m3,
        );
    }
}

fn nonnegativeConcentration(amount_mol: f64, water_volume_m3: f64) f64 {
    return @max(0, amount_mol) / water_volume_m3;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(PrecipitateAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.precipitates, field.name)))
            return error.InvalidSolidConcentrationInput;
    }
    inline for (@typeInfo(SilicateAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.natural_silicates, field.name)) or
            !std.math.isFinite(@field(inputs.ground_silicates, field.name)))
        {
            return error.InvalidSolidConcentrationInput;
        }
    }
    if (!std.math.isFinite(inputs.water_volume_m3) or
        inputs.water_volume_m3 < 0 or
        !std.math.isFinite(inputs.minimum_water_volume_m3) or
        inputs.minimum_water_volume_m3 < 0 or
        !std.math.isFinite(inputs.ground_silicate_specific_surface_area_m2_per_mol) or
        inputs.ground_silicate_specific_surface_area_m2_per_mol < 0)
    {
        return error.InvalidSolidConcentrationInput;
    }
}

fn validateState(state: geochemistry.SolidState) !void {
    inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field| {
        const value = @field(state, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSolidConcentrationState;
    }
}

fn silicates(scale: f64) SilicateAmounts {
    return .{
        .aluminum_mol = 1 * scale,
        .iron_mol = 2 * scale,
        .calcium_mol = 3 * scale,
        .magnesium_mol = 4 * scale,
        .sodium_mol = 5 * scale,
        .potassium_mol = 6 * scale,
    };
}

fn sourceInputs() Inputs {
    return .{
        .precipitates = .{
            .gibbsite_mol = 2,
            .iron_hydroxide_mol = 4,
            .calcite_mol = 6,
            .gypsum_mol = 8,
        },
        .natural_silicates = silicates(2),
        .ground_silicates = silicates(4),
        .water_volume_m3 = 2,
        .minimum_water_volume_m3 = 0,
        .ground_silicate_specific_surface_area_m2_per_mol = 1000,
    };
}

test "SOLUTE solid concentrations preserve source order and ground surface area" {
    var state = std.mem.zeroes(geochemistry.SolidState);
    var source = sourceInputs();
    source.precipitates.iron_hydroxide_mol = -1;
    const result = try applyConcentrations(&state, source);

    try std.testing.expectEqual(@as(f64, 1), state.gibbsite_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0), state.iron_hydroxide_solid_mol_per_m3);
    try std.testing.expectEqual(
        @as(f64, 2),
        state.aluminum_ground_silicate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 12),
        state.potassium_ground_silicate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 42_000),
        result.ground_silicate_surface_area_m2_per_m3,
    );
}

test "SOLUTE dry solid concentrations clear all geochemistry solids" {
    var state: geochemistry.SolidState = undefined;
    inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field|
        @field(state, field.name) = 7;
    var source = sourceInputs();
    source.water_volume_m3 = 1.0e-9;
    source.minimum_water_volume_m3 = 1.0e-9;
    const result = try applyConcentrations(&state, source);
    inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(state, field.name));
    try std.testing.expectEqual(
        @as(f64, 0),
        result.ground_silicate_surface_area_m2_per_m3,
    );
}

test "SOLUTE invalid late solid input leaves state unchanged" {
    var state: geochemistry.SolidState = undefined;
    inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field|
        @field(state, field.name) = 3;
    const before = state;
    var source = sourceInputs();
    source.ground_silicates.potassium_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSolidConcentrationInput,
        applyConcentrations(&state, source),
    );
    try std.testing.expectEqualDeep(before, state);
}
