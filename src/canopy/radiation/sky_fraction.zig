const std = @import("std");

pub const SpeciesGeometry = struct {
    standing_dead_area_m2_by_layer: []const f64,
    leaf_area_m2_by_layer_branch_orientation: []const f64,
    stem_area_m2_by_layer_branch: []const f64,
    branch_count: usize,
    leaf_orientation_count: usize,
};

pub const SpeciesRadiation = struct {
    live_canopy_radiation_megajoules_h: f64,
    standing_dead_radiation_megajoules_h: f64,
};

pub const SpeciesThermodynamics = struct {
    live_surface_temperature_k: f64,
    standing_dead_surface_temperature_k: f64,
    live_air_temperature_k: f64,
    standing_dead_air_temperature_k: f64,
    live_vapor_pressure_kpa: f64,
    standing_dead_vapor_pressure_kpa: f64,
};

pub const Inputs = struct {
    canopy_layer_bottom_height_m: []const f64,
    snow_depth_m: f64,
    surface_water_depth_m: f64,
    height_tolerance_m: f64,
    sine_solar_elevation: f64,
    daylight_threshold: f64,
    total_canopy_radiation_megajoules_h: f64,
    ground_radiation_megajoules_h: f64,
    cell_surface_area_m2: f64,
    division_threshold: f64,
    species_geometry: []const SpeciesGeometry,
    species_radiation: []const SpeciesRadiation,
    species_thermodynamics: []const SpeciesThermodynamics,
    initial_bulk_surface_temperature_k: f64,
    initial_bulk_air_temperature_k: f64,
    initial_bulk_vapor_pressure_kpa: f64,
};

pub const SpeciesOutputs = struct {
    exposed_live_area_m2: []f64,
    exposed_standing_dead_area_m2: []f64,
    live_radiation_fraction: []f64,
    standing_dead_radiation_fraction: []f64,
    live_canopy_fraction: []f64,
    standing_dead_canopy_fraction: []f64,
};

pub const Result = struct {
    total_exposed_canopy_area_m2: f64,
    canopy_radiation_fraction: f64,
    ground_radiation_fraction: f64,
    bulk_surface_temperature_k: f64,
    bulk_air_temperature_k: f64,
    bulk_vapor_pressure_kpa: f64,
};

pub const CalculationError = error{
    SpeciesCountMismatch,
    OutputCountMismatch,
    LayerCountMismatch,
    GeometryExtentOverflow,
    GeometryExtentMismatch,
    NonFiniteInput,
    NegativeArea,
    NonPositiveCellArea,
    NegativeThreshold,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 4713--4819 for one grid cell.
///
/// Species, canopy layers, branches, and leaf-orientation classes all have
/// runtime extents. Flattened leaf indexing is layer-major, then branch, then
/// orientation; flattened stem indexing is layer-major, then branch.
pub fn calculate(inputs: Inputs, outputs: SpeciesOutputs) CalculationError!Result {
    const species_count = inputs.species_geometry.len;
    if (inputs.species_radiation.len != species_count or
        inputs.species_thermodynamics.len != species_count)
    {
        return error.SpeciesCountMismatch;
    }
    inline for (std.meta.fields(SpeciesOutputs)) |field| {
        if (@field(outputs, field.name).len != species_count) {
            return error.OutputCountMismatch;
        }
    }
    try validateInputs(inputs);

    var total_exposed_area_m2: f64 = 0.0;
    for (inputs.species_geometry, 0..) |geometry, species_index| {
        outputs.exposed_live_area_m2[species_index] = 0.0;
        outputs.exposed_standing_dead_area_m2[species_index] = 0.0;
        for (inputs.canopy_layer_bottom_height_m, 0..) |layer_bottom_m, layer_index| {
            if (layer_bottom_m >= inputs.snow_depth_m - inputs.height_tolerance_m and
                layer_bottom_m >= inputs.surface_water_depth_m - inputs.height_tolerance_m)
            {
                const standing_dead_area_m2 =
                    geometry.standing_dead_area_m2_by_layer[layer_index];
                outputs.exposed_standing_dead_area_m2[species_index] += standing_dead_area_m2;
                total_exposed_area_m2 += standing_dead_area_m2;

                for (0..geometry.branch_count) |branch_index| {
                    for (0..geometry.leaf_orientation_count) |orientation_index| {
                        const leaf_index = (layer_index * geometry.branch_count + branch_index) *
                            geometry.leaf_orientation_count + orientation_index;
                        const leaf_area_m2 =
                            geometry.leaf_area_m2_by_layer_branch_orientation[leaf_index];
                        outputs.exposed_live_area_m2[species_index] += leaf_area_m2;
                        total_exposed_area_m2 += leaf_area_m2;
                    }
                    const stem_index = layer_index * geometry.branch_count + branch_index;
                    const stem_area_m2 = geometry.stem_area_m2_by_layer_branch[stem_index];
                    outputs.exposed_live_area_m2[species_index] += stem_area_m2;
                    total_exposed_area_m2 += stem_area_m2;
                }
            }
        }
    }

    var total_radiation_megajoules_h: f64 = undefined;
    var canopy_fraction: f64 = 0.0;
    var ground_fraction: f64 = 1.0;
    if (inputs.sine_solar_elevation > inputs.daylight_threshold) {
        total_radiation_megajoules_h =
            inputs.total_canopy_radiation_megajoules_h + inputs.ground_radiation_megajoules_h;
        if (total_radiation_megajoules_h > inputs.division_threshold) {
            for (inputs.species_radiation, 0..) |radiation, species_index| {
                outputs.live_radiation_fraction[species_index] =
                    radiation.live_canopy_radiation_megajoules_h / total_radiation_megajoules_h;
                outputs.standing_dead_radiation_fraction[species_index] =
                    radiation.standing_dead_radiation_megajoules_h / total_radiation_megajoules_h;
                canopy_fraction += outputs.live_radiation_fraction[species_index] +
                    outputs.standing_dead_radiation_fraction[species_index];
                ground_fraction -= outputs.live_radiation_fraction[species_index] +
                    outputs.standing_dead_radiation_fraction[species_index];
            }
        } else {
            total_radiation_megajoules_h = 0.0;
            zeroRadiationFractions(outputs);
        }
    } else if (total_exposed_area_m2 > inputs.division_threshold) {
        total_radiation_megajoules_h =
            1.0 - @exp(-0.65 * total_exposed_area_m2 / inputs.cell_surface_area_m2);
        for (0..species_count) |species_index| {
            outputs.live_radiation_fraction[species_index] = total_radiation_megajoules_h *
                outputs.exposed_live_area_m2[species_index] / total_exposed_area_m2;
            outputs.standing_dead_radiation_fraction[species_index] = total_radiation_megajoules_h *
                outputs.exposed_standing_dead_area_m2[species_index] / total_exposed_area_m2;
            canopy_fraction += outputs.live_radiation_fraction[species_index] +
                outputs.standing_dead_radiation_fraction[species_index];
            ground_fraction -= outputs.live_radiation_fraction[species_index] +
                outputs.standing_dead_radiation_fraction[species_index];
        }
    } else {
        total_radiation_megajoules_h = 0.0;
        zeroRadiationFractions(outputs);
    }

    var bulk_surface_temperature_k = inputs.initial_bulk_surface_temperature_k;
    var bulk_air_temperature_k = inputs.initial_bulk_air_temperature_k;
    var bulk_vapor_pressure_kpa = inputs.initial_bulk_vapor_pressure_kpa;
    for (inputs.species_thermodynamics, 0..) |thermodynamics, species_index| {
        if (canopy_fraction > inputs.height_tolerance_m) {
            outputs.live_canopy_fraction[species_index] =
                outputs.live_radiation_fraction[species_index] / canopy_fraction;
            outputs.standing_dead_canopy_fraction[species_index] =
                outputs.standing_dead_radiation_fraction[species_index] / canopy_fraction;
        } else {
            outputs.live_canopy_fraction[species_index] = 0.0;
            outputs.standing_dead_canopy_fraction[species_index] = 0.0;
        }
        bulk_surface_temperature_k += thermodynamics.live_surface_temperature_k *
            outputs.live_canopy_fraction[species_index] +
            thermodynamics.standing_dead_surface_temperature_k *
                outputs.standing_dead_canopy_fraction[species_index];
        bulk_air_temperature_k += thermodynamics.live_air_temperature_k *
            outputs.live_canopy_fraction[species_index] +
            thermodynamics.standing_dead_air_temperature_k *
                outputs.standing_dead_canopy_fraction[species_index];
        bulk_vapor_pressure_kpa += thermodynamics.live_vapor_pressure_kpa *
            outputs.live_canopy_fraction[species_index] +
            thermodynamics.standing_dead_vapor_pressure_kpa *
                outputs.standing_dead_canopy_fraction[species_index];
    }

    const result = Result{
        .total_exposed_canopy_area_m2 = total_exposed_area_m2,
        .canopy_radiation_fraction = canopy_fraction,
        .ground_radiation_fraction = ground_fraction,
        .bulk_surface_temperature_k = bulk_surface_temperature_k,
        .bulk_air_temperature_k = bulk_air_temperature_k,
        .bulk_vapor_pressure_kpa = bulk_vapor_pressure_kpa,
    };
    inline for (std.meta.fields(Result)) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteResult;
    }
    return result;
}

fn validateInputs(inputs: Inputs) CalculationError!void {
    const scalar_values = [_]f64{
        inputs.snow_depth_m,
        inputs.surface_water_depth_m,
        inputs.height_tolerance_m,
        inputs.sine_solar_elevation,
        inputs.daylight_threshold,
        inputs.total_canopy_radiation_megajoules_h,
        inputs.ground_radiation_megajoules_h,
        inputs.cell_surface_area_m2,
        inputs.division_threshold,
        inputs.initial_bulk_surface_temperature_k,
        inputs.initial_bulk_air_temperature_k,
        inputs.initial_bulk_vapor_pressure_kpa,
    };
    for (scalar_values) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
    }
    if (inputs.cell_surface_area_m2 <= 0.0) return error.NonPositiveCellArea;
    if (inputs.height_tolerance_m < 0.0 or inputs.division_threshold < 0.0) {
        return error.NegativeThreshold;
    }
    for (inputs.canopy_layer_bottom_height_m) |height_m| {
        if (!std.math.isFinite(height_m)) return error.NonFiniteInput;
    }
    for (inputs.species_geometry) |geometry| {
        try validateGeometry(geometry, inputs.canopy_layer_bottom_height_m.len);
    }
    for (inputs.species_radiation) |radiation| {
        inline for (std.meta.fields(SpeciesRadiation)) |field| {
            if (!std.math.isFinite(@field(radiation, field.name))) return error.NonFiniteInput;
        }
    }
    for (inputs.species_thermodynamics) |thermodynamics| {
        inline for (std.meta.fields(SpeciesThermodynamics)) |field| {
            if (!std.math.isFinite(@field(thermodynamics, field.name))) {
                return error.NonFiniteInput;
            }
        }
    }
}

fn validateGeometry(geometry: SpeciesGeometry, layer_count: usize) CalculationError!void {
    if (geometry.standing_dead_area_m2_by_layer.len != layer_count) {
        return error.LayerCountMismatch;
    }
    const stem_count = std.math.mul(usize, layer_count, geometry.branch_count) catch
        return error.GeometryExtentOverflow;
    const leaf_count = std.math.mul(
        usize,
        stem_count,
        geometry.leaf_orientation_count,
    ) catch return error.GeometryExtentOverflow;
    if (geometry.stem_area_m2_by_layer_branch.len != stem_count or
        geometry.leaf_area_m2_by_layer_branch_orientation.len != leaf_count)
    {
        return error.GeometryExtentMismatch;
    }
    const area_groups = [_][]const f64{
        geometry.standing_dead_area_m2_by_layer,
        geometry.stem_area_m2_by_layer_branch,
        geometry.leaf_area_m2_by_layer_branch_orientation,
    };
    for (area_groups) |areas| {
        for (areas) |area_m2| {
            if (!std.math.isFinite(area_m2)) return error.NonFiniteInput;
            if (area_m2 < 0.0) return error.NegativeArea;
        }
    }
}

fn zeroRadiationFractions(outputs: SpeciesOutputs) void {
    @memset(outputs.live_radiation_fraction, 0.0);
    @memset(outputs.standing_dead_radiation_fraction, 0.0);
}

test "daylight fractions use radiation and aggregate canopy thermodynamics" {
    const layer_bottoms_m = [_]f64{ 0.0, 1.0 };
    const geometry = [_]SpeciesGeometry{.{
        .standing_dead_area_m2_by_layer = &.{ 1.0, 2.0 },
        .leaf_area_m2_by_layer_branch_orientation = &.{ 3.0, 4.0 },
        .stem_area_m2_by_layer_branch = &.{ 5.0, 6.0 },
        .branch_count = 1,
        .leaf_orientation_count = 1,
    }};
    var live_area: [1]f64 = undefined;
    var dead_area: [1]f64 = undefined;
    var live_radiation_fraction: [1]f64 = undefined;
    var dead_radiation_fraction: [1]f64 = undefined;
    var live_canopy_fraction: [1]f64 = undefined;
    var dead_canopy_fraction: [1]f64 = undefined;
    const result = try calculate(.{
        .canopy_layer_bottom_height_m = &layer_bottoms_m,
        .snow_depth_m = 0.5,
        .surface_water_depth_m = 0.0,
        .height_tolerance_m = 1.0e-12,
        .sine_solar_elevation = 0.5,
        .daylight_threshold = 0.05,
        .total_canopy_radiation_megajoules_h = 7.0,
        .ground_radiation_megajoules_h = 3.0,
        .cell_surface_area_m2 = 10.0,
        .division_threshold = 1.0e-12,
        .species_geometry = &geometry,
        .species_radiation = &.{.{
            .live_canopy_radiation_megajoules_h = 4.0,
            .standing_dead_radiation_megajoules_h = 2.0,
        }},
        .species_thermodynamics = &.{.{
            .live_surface_temperature_k = 300.0,
            .standing_dead_surface_temperature_k = 280.0,
            .live_air_temperature_k = 295.0,
            .standing_dead_air_temperature_k = 285.0,
            .live_vapor_pressure_kpa = 2.0,
            .standing_dead_vapor_pressure_kpa = 1.0,
        }},
        .initial_bulk_surface_temperature_k = 0.0,
        .initial_bulk_air_temperature_k = 0.0,
        .initial_bulk_vapor_pressure_kpa = 0.0,
    }, .{
        .exposed_live_area_m2 = &live_area,
        .exposed_standing_dead_area_m2 = &dead_area,
        .live_radiation_fraction = &live_radiation_fraction,
        .standing_dead_radiation_fraction = &dead_radiation_fraction,
        .live_canopy_fraction = &live_canopy_fraction,
        .standing_dead_canopy_fraction = &dead_canopy_fraction,
    });

    try std.testing.expectEqual(@as(f64, 12.0), result.total_exposed_canopy_area_m2);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result.canopy_radiation_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result.ground_radiation_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, 293.3333333333333),
        result.bulk_surface_temperature_k,
        1.0e-12,
    );
}

test "night fractions use exposed area with runtime species count" {
    const geometry = [_]SpeciesGeometry{
        .{
            .standing_dead_area_m2_by_layer = &.{1.0},
            .leaf_area_m2_by_layer_branch_orientation = &.{2.0},
            .stem_area_m2_by_layer_branch = &.{1.0},
            .branch_count = 1,
            .leaf_orientation_count = 1,
        },
        .{
            .standing_dead_area_m2_by_layer = &.{2.0},
            .leaf_area_m2_by_layer_branch_orientation = &.{3.0},
            .stem_area_m2_by_layer_branch = &.{1.0},
            .branch_count = 1,
            .leaf_orientation_count = 1,
        },
    };
    var output_storage: [12]f64 = undefined;
    const zero_thermodynamics = SpeciesThermodynamics{
        .live_surface_temperature_k = 0.0,
        .standing_dead_surface_temperature_k = 0.0,
        .live_air_temperature_k = 0.0,
        .standing_dead_air_temperature_k = 0.0,
        .live_vapor_pressure_kpa = 0.0,
        .standing_dead_vapor_pressure_kpa = 0.0,
    };
    const result = try calculate(.{
        .canopy_layer_bottom_height_m = &.{0.0},
        .snow_depth_m = 0.0,
        .surface_water_depth_m = 0.0,
        .height_tolerance_m = 0.0,
        .sine_solar_elevation = 0.0,
        .daylight_threshold = 0.05,
        .total_canopy_radiation_megajoules_h = 0.0,
        .ground_radiation_megajoules_h = 0.0,
        .cell_surface_area_m2 = 20.0,
        .division_threshold = 0.0,
        .species_geometry = &geometry,
        .species_radiation = &.{ std.mem.zeroes(SpeciesRadiation), std.mem.zeroes(SpeciesRadiation) },
        .species_thermodynamics = &.{ zero_thermodynamics, zero_thermodynamics },
        .initial_bulk_surface_temperature_k = 0.0,
        .initial_bulk_air_temperature_k = 0.0,
        .initial_bulk_vapor_pressure_kpa = 0.0,
    }, .{
        .exposed_live_area_m2 = output_storage[0..2],
        .exposed_standing_dead_area_m2 = output_storage[2..4],
        .live_radiation_fraction = output_storage[4..6],
        .standing_dead_radiation_fraction = output_storage[6..8],
        .live_canopy_fraction = output_storage[8..10],
        .standing_dead_canopy_fraction = output_storage[10..12],
    });

    const expected_total_fraction = 1.0 - @exp(-0.65 * 10.0 / 20.0);
    try std.testing.expectApproxEqRel(
        expected_total_fraction,
        result.canopy_radiation_fraction,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        expected_total_fraction * 3.0 / 10.0,
        output_storage[4],
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        expected_total_fraction * 2.0 / 10.0,
        output_storage[7],
        1.0e-14,
    );
}
