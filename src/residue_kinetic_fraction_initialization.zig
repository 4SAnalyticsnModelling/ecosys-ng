const std = @import("std");

pub const LayerKind = enum {
    surface_litter,
    mineral_soil,
};

pub const PlantResidueKind = enum {
    maize,
    wheat,
    soybean,
    new_straw,
    old_straw,
    compost,
    green_manure,
    new_deciduous_forest,
    new_coniferous_forest,
    old_deciduous_forest,
    old_coniferous_forest,
    other,
};

pub const ManureKind = enum {
    ruminant,
    non_ruminant,
    other,
};

/// Fractions are ordered as protein, nonstructural carbohydrate,
/// cellulose, and lignin. Each slice must contain exactly four entries.
pub const State = struct {
    previous_coarse_residue_fraction: []f64,
    plant_or_fine_root_residue_fraction: []f64,
    manure_fraction: []f64,
    particulate_organic_matter_fraction: []f64,
};

pub const Inputs = struct {
    layer_kind: LayerKind,
    plant_residue_kind: PlantResidueKind,
    manure_kind: ManureKind,
};

const FractionSet = [4]f64;

fn surfacePlantFractions(kind: PlantResidueKind) FractionSet {
    return switch (kind) {
        .maize => .{ 0.080, 0.245, 0.613, 0.062 },
        .wheat => .{ 0.125, 0.171, 0.560, 0.144 },
        .soybean => .{ 0.138, 0.426, 0.316, 0.120 },
        .new_straw => .{ 0.036, 0.044, 0.767, 0.153 },
        .old_straw => .{ 0.075, 0.125, 0.550, 0.250 },
        .compost => .{ 0.143, 0.015, 0.640, 0.202 },
        .green_manure => .{ 0.202, 0.013, 0.560, 0.225 },
        .new_deciduous_forest => .{ 0.07, 0.41, 0.36, 0.16 },
        .new_coniferous_forest => .{ 0.07, 0.25, 0.38, 0.30 },
        .old_deciduous_forest, .old_coniferous_forest => .{
            0.02,
            0.06,
            0.34,
            0.58,
        },
        .other => .{ 0.075, 0.125, 0.550, 0.250 },
    };
}

fn manureFractions(kind: ManureKind) FractionSet {
    return switch (kind) {
        .ruminant => .{ 0.036, 0.044, 0.630, 0.290 },
        .non_ruminant, .other => .{ 0.138, 0.401, 0.316, 0.145 },
    };
}

fn validateFractions(fractions: FractionSet) !void {
    var total: f64 = 0.0;
    for (fractions) |fraction| {
        if (!std.math.isFinite(fraction))
            return error.NonFiniteResidueKineticFraction;
        if (fraction < 0.0 or fraction > 1.0)
            return error.InvalidResidueKineticFraction;
        total += fraction;
    }
    if (@abs(total - 1.0) > 1.0e-12)
        return error.InvalidResidueKineticFractionSum;
}

fn copyFractions(destination: []f64, source: FractionSet) void {
    @memcpy(destination, source[0..]);
}

/// Exact source-order translation of legacy `STARTS` lines 840--990.
///
/// Fractions are dimensionless. Surface particulate organic matter is
/// deliberately left unchanged, matching the legacy `L.NE.0` branch.
pub fn initialize(state: State, inputs: Inputs) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (@field(state, field.name).len != 4)
            return error.ResidueKineticFractionDimensionMismatch;
    }

    const previous_coarse: FractionSet = switch (inputs.layer_kind) {
        .surface_litter => .{ 0.000, 0.045, 0.660, 0.295 },
        .mineral_soil => .{ 0.00, 0.00, 0.20, 0.80 },
    };
    const plant_or_fine_root: FractionSet = switch (inputs.layer_kind) {
        .surface_litter => surfacePlantFractions(inputs.plant_residue_kind),
        .mineral_soil => .{ 0.02, 0.06, 0.34, 0.58 },
    };
    const manure = manureFractions(inputs.manure_kind);
    const particulate: FractionSet = .{ 1.00, 0.00, 0.00, 0.00 };

    try validateFractions(previous_coarse);
    try validateFractions(plant_or_fine_root);
    try validateFractions(manure);
    if (inputs.layer_kind == .mineral_soil)
        try validateFractions(particulate);

    copyFractions(state.previous_coarse_residue_fraction, previous_coarse);
    copyFractions(
        state.plant_or_fine_root_residue_fraction,
        plant_or_fine_root,
    );
    copyFractions(state.manure_fraction, manure);
    if (inputs.layer_kind == .mineral_soil)
        copyFractions(state.particulate_organic_matter_fraction, particulate);
}

test "STARTS surface fractions preserve selected litter and POM state" {
    var previous_coarse = [_]f64{9.0} ** 4;
    var plant = [_]f64{9.0} ** 4;
    var manure = [_]f64{9.0} ** 4;
    var particulate = [_]f64{9.0} ** 4;
    try initialize(.{
        .previous_coarse_residue_fraction = &previous_coarse,
        .plant_or_fine_root_residue_fraction = &plant,
        .manure_fraction = &manure,
        .particulate_organic_matter_fraction = &particulate,
    }, .{
        .layer_kind = .surface_litter,
        .plant_residue_kind = .soybean,
        .manure_kind = .ruminant,
    });

    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.000, 0.045, 0.660, 0.295 },
        &previous_coarse,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.138, 0.426, 0.316, 0.120 },
        &plant,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.036, 0.044, 0.630, 0.290 },
        &manure,
    );
    try std.testing.expectEqualSlices(f64, &.{ 9, 9, 9, 9 }, &particulate);
}

test "STARTS mineral fractions initialize roots manure and POM" {
    var previous_coarse = [_]f64{9.0} ** 4;
    var plant = [_]f64{9.0} ** 4;
    var manure = [_]f64{9.0} ** 4;
    var particulate = [_]f64{9.0} ** 4;
    try initialize(.{
        .previous_coarse_residue_fraction = &previous_coarse,
        .plant_or_fine_root_residue_fraction = &plant,
        .manure_fraction = &manure,
        .particulate_organic_matter_fraction = &particulate,
    }, .{
        .layer_kind = .mineral_soil,
        .plant_residue_kind = .maize,
        .manure_kind = .other,
    });

    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.00, 0.00, 0.20, 0.80 },
        &previous_coarse,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.02, 0.06, 0.34, 0.58 },
        &plant,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0.138, 0.401, 0.316, 0.145 },
        &manure,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 1.00, 0.00, 0.00, 0.00 },
        &particulate,
    );
}

test "dimension failure is atomic" {
    var short = [_]f64{ 7, 7, 7 };
    var plant = [_]f64{8.0} ** 4;
    var manure = [_]f64{8.0} ** 4;
    var particulate = [_]f64{8.0} ** 4;
    try std.testing.expectError(
        error.ResidueKineticFractionDimensionMismatch,
        initialize(.{
            .previous_coarse_residue_fraction = &short,
            .plant_or_fine_root_residue_fraction = &plant,
            .manure_fraction = &manure,
            .particulate_organic_matter_fraction = &particulate,
        }, .{
            .layer_kind = .mineral_soil,
            .plant_residue_kind = .wheat,
            .manure_kind = .non_ruminant,
        }),
    );
    try std.testing.expectEqualSlices(f64, &.{ 8, 8, 8, 8 }, &plant);
}
