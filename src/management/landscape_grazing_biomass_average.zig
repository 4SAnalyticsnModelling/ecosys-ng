const std = @import("std");

pub const HarvestKind = enum {
    other,
    animal_grazing,
    insect_grazing,
};

pub const Inputs = struct {
    lon_count: usize,
    lat_count: usize,
    species_count: usize,
    target_cell: usize,
    target_species: usize,
    harvest_kind: HarvestKind,
    landscape_section_by_plant: []const i32,
    active_by_plant: []const bool,
    shoot_carbon_g_c_by_plant: []const f64,
};

/// Exact grosub.f lines 247--264 landscape-section grazing average.
///
/// Plant arrays use cell-major, species-minor runtime layout. The result is
/// `null` outside the source animal/insect grazing gate. Within that gate,
/// only active plants of the target species and landscape section contribute;
/// an empty section retains the target plant's shoot carbon (g C).
pub fn calculate(inputs: Inputs) !?f64 {
    if (inputs.harvest_kind != .animal_grazing and
        inputs.harvest_kind != .insect_grazing)
        return null;
    if (inputs.lon_count == 0 or inputs.lat_count == 0 or
        inputs.species_count == 0)
        return error.EmptyGrazingLandscape;
    const cell_count = try std.math.mul(
        usize,
        inputs.lon_count,
        inputs.lat_count,
    );
    const plant_count = try std.math.mul(
        usize,
        cell_count,
        inputs.species_count,
    );
    if (inputs.landscape_section_by_plant.len != plant_count or
        inputs.active_by_plant.len != plant_count or
        inputs.shoot_carbon_g_c_by_plant.len != plant_count)
        return error.GrazingLandscapeDimensionMismatch;
    if (inputs.target_cell >= cell_count or
        inputs.target_species >= inputs.species_count)
        return error.GrazingTargetOutOfBounds;
    const target = inputs.target_cell * inputs.species_count +
        inputs.target_species;
    if (!std.math.isFinite(inputs.shoot_carbon_g_c_by_plant[target]))
        return error.NonFiniteGrazingShootCarbon;
    const target_section = inputs.landscape_section_by_plant[target];
    var total_shoot_carbon_g_c: f64 = 0;
    var active_cell_count: usize = 0;
    // Preserve GROSUB's NX-outer, NY-inner accumulation order while using
    // the production row-major cell layout.
    for (0..inputs.lon_count) |column| {
        for (0..inputs.lat_count) |row| {
            const cell = row * inputs.lon_count + column;
            const plant = cell * inputs.species_count + inputs.target_species;
            if (inputs.landscape_section_by_plant[plant] == target_section and
                inputs.active_by_plant[plant])
            {
                if (!std.math.isFinite(inputs.shoot_carbon_g_c_by_plant[plant]))
                    return error.NonFiniteGrazingShootCarbon;
                total_shoot_carbon_g_c +=
                    inputs.shoot_carbon_g_c_by_plant[plant];
                active_cell_count += 1;
            }
        }
    }
    if (!std.math.isFinite(total_shoot_carbon_g_c))
        return error.GrazingShootCarbonOverflow;
    return if (active_cell_count > 0)
        total_shoot_carbon_g_c / @as(f64, @floatFromInt(active_cell_count))
    else
        inputs.shoot_carbon_g_c_by_plant[target];
}

test "GROSUB averages target species across active cells in its section" {
    const result = (try calculate(.{
        .lon_count = 2,
        .lat_count = 2,
        .species_count = 7,
        .target_cell = 1,
        .target_species = 6,
        .harvest_kind = .animal_grazing,
        .landscape_section_by_plant = &sections,
        .active_by_plant = &active,
        .shoot_carbon_g_c_by_plant = &shoot,
    })).?;
    try std.testing.expectEqual(@as(f64, 30), result);
}

test "empty active section falls back to target shoot carbon" {
    var sections_local = [_]i32{2} ** 6;
    var active_local = [_]bool{false} ** 6;
    var shoot_local = [_]f64{ 3, 4, 5, 6, 7, 8 };
    const result = (try calculate(.{
        .lon_count = 3,
        .lat_count = 1,
        .species_count = 2,
        .target_cell = 2,
        .target_species = 1,
        .harvest_kind = .insect_grazing,
        .landscape_section_by_plant = &sections_local,
        .active_by_plant = &active_local,
        .shoot_carbon_g_c_by_plant = &shoot_local,
    })).?;
    try std.testing.expectEqual(@as(f64, 8), result);
}

test "non-grazing event leaves source average unpublished" {
    try std.testing.expectEqual(@as(?f64, null), try calculate(.{
        .lon_count = 1,
        .lat_count = 1,
        .species_count = 1,
        .target_cell = 0,
        .target_species = 0,
        .harvest_kind = .other,
        .landscape_section_by_plant = &.{1},
        .active_by_plant = &.{true},
        .shoot_carbon_g_c_by_plant = &.{5},
    }));
}

test "invalid late shoot value fails before publishing an average" {
    try std.testing.expectError(
        error.NonFiniteGrazingShootCarbon,
        calculate(.{
            .lon_count = 2,
            .lat_count = 1,
            .species_count = 1,
            .target_cell = 0,
            .target_species = 0,
            .harvest_kind = .animal_grazing,
            .landscape_section_by_plant = &.{ 1, 1 },
            .active_by_plant = &.{ true, true },
            .shoot_carbon_g_c_by_plant = &.{ 5, std.math.nan(f64) },
        }),
    );
}

test "irrelevant species carbon cannot contaminate the source traversal" {
    const result = (try calculate(.{
        .lon_count = 1,
        .lat_count = 1,
        .species_count = 2,
        .target_cell = 0,
        .target_species = 0,
        .harvest_kind = .animal_grazing,
        .landscape_section_by_plant = &.{ 1, 1 },
        .active_by_plant = &.{ true, true },
        .shoot_carbon_g_c_by_plant = &.{ 5, std.math.nan(f64) },
    })).?;
    try std.testing.expectEqual(@as(f64, 5), result);
}

test "GROSUB column outer row inner accumulation order is retained" {
    const result = (try calculate(.{
        .lon_count = 2,
        .lat_count = 2,
        .species_count = 1,
        .target_cell = 0,
        .target_species = 0,
        .harvest_kind = .animal_grazing,
        .landscape_section_by_plant = &.{ 1, 1, 1, 1 },
        .active_by_plant = &.{ true, true, true, true },
        .shoot_carbon_g_c_by_plant = &.{ 1.0e16, 1, -1.0e16, 1 },
    })).?;
    try std.testing.expectEqual(@as(f64, 0.5), result);
}

const sections = [_]i32{
    9, 9, 9, 9, 9, 9, 4,
    9, 9, 9, 9, 9, 9, 4,
    9, 9, 9, 9, 9, 9, 4,
    9, 9, 9, 9, 9, 9, 8,
};

const active = [_]bool{
    false, false, false, false, false, false, true,
    false, false, false, false, false, false, false,
    false, false, false, false, false, false, true,
    false, false, false, false, false, false, true,
};

const shoot = [_]f64{
    1, 1, 1, 1, 1, 1, 20,
    1, 1, 1, 1, 1, 1, 99,
    1, 1, 1, 1, 1, 1, 40,
    1, 1, 1, 1, 1, 1, 80,
};
