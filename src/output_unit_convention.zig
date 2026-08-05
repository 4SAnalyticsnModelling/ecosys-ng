const std = @import("std");
const soil_output_catalog = @import("soil_output_catalog.zig");
const plant_output_catalog = @import("plant_output_catalog.zig");

/// Dimensionless qualifiers that are deliberately words rather than SI symbols.
const dimensionless_labels = [_][]const u8{
    "count",
    "enum_code",
    "fraction",
    "model_concentration",
    "stage",
};

fn isDimensionless(unit: []const u8) bool {
    for (dimensionless_labels) |label| if (std.mem.eql(u8, unit, label)) return true;
    return false;
}

/// Rejects `per`-style identifier spellings inside a heading. Case-carrying SI
/// symbols such as `MJ` and `Mg` are correct in a heading, because a heading is
/// a unit string rather than a snake_case identifier and its case is preserved
/// verbatim. The ambiguity rule that forces `megajoule`/`megagram` applies to
/// identifiers, which is checked separately by
/// `tools/inventory_ambiguous_unit_identifiers.py`. What a heading must never
/// contain is a lowercase `mg`/`mj`, since those read as milli- rather than
/// mega- prefixes.
fn checkUnit(unit: []const u8) !void {
    if (unit.len == 0) return error.EmptyOutputUnit;
    if (isDimensionless(unit)) return;
    if (std.mem.indexOf(u8, unit, "_per_") != null) return error.NonSiOutputUnit;
    if (std.mem.indexOfScalar(u8, unit, '_') != null) return error.NonSiOutputUnit;
    var symbols = std.mem.splitScalar(u8, unit, ' ');
    while (symbols.next()) |symbol| {
        if (symbol.len == 0) return error.NonSiOutputUnit;
        // A lowercase `mg`/`mj` is a milligram/millijoule reading of a quantity
        // that is actually mega-scaled, so require the correct capitalization.
        for ([_][]const u8{ "mg", "mj", "Mj", "MG" }) |ambiguous|
            if (std.mem.eql(u8, symbol, ambiguous)) return error.AmbiguousOutputUnit;
    }
}

fn checkCatalog(catalog: anytype) !void {
    for (catalog.variables) |variable| {
        checkUnit(variable.unit) catch |err| {
            std.log.err("invalid output unit: name={s} unit={s}", .{ variable.name, variable.unit });
            return err;
        };
    }
}

test "every soil output catalog heading uses SI or dimensionless unit labels" {
    const allocator = std.testing.allocator;
    var soil_water = try soil_output_catalog.water(allocator, 11);
    defer soil_water.deinit();
    try checkCatalog(soil_water);
    var soil_heat = try soil_output_catalog.heat(allocator, 11);
    defer soil_heat.deinit();
    try checkCatalog(soil_heat);
    var soil_carbon = try soil_output_catalog.carbon(allocator, 11, 11, 11);
    defer soil_carbon.deinit();
    try checkCatalog(soil_carbon);
    var soil_nitrogen = try soil_output_catalog.nitrogen(allocator, 11, 11);
    defer soil_nitrogen.deinit();
    try checkCatalog(soil_nitrogen);
    var soil_phosphorus = try soil_output_catalog.phosphorus(allocator);
    defer soil_phosphorus.deinit();
    try checkCatalog(soil_phosphorus);
    var daily_carbon = try soil_output_catalog.dailyCarbon(allocator, 11);
    defer daily_carbon.deinit();
    try checkCatalog(daily_carbon);
    var daily_water = try soil_output_catalog.dailyWater(allocator, 11, 11, 11);
    defer daily_water.deinit();
    try checkCatalog(daily_water);
    var daily_nitrogen = try soil_output_catalog.dailyNitrogen(allocator, 11);
    defer daily_nitrogen.deinit();
    try checkCatalog(daily_nitrogen);
    var daily_phosphorus = try soil_output_catalog.dailyPhosphorus(allocator, 11);
    defer daily_phosphorus.deinit();
    try checkCatalog(daily_phosphorus);
    var daily_heat = try soil_output_catalog.dailyHeat(allocator, 11, 11);
    defer daily_heat.deinit();
    try checkCatalog(daily_heat);
}

test "every plant output catalog heading uses SI or dimensionless unit labels" {
    const allocator = std.testing.allocator;
    var plant_carbon = try plant_output_catalog.carbon(allocator);
    defer plant_carbon.deinit();
    try checkCatalog(plant_carbon);
    var plant_water = try plant_output_catalog.water(allocator, 11);
    defer plant_water.deinit();
    try checkCatalog(plant_water);
    var plant_nitrogen = try plant_output_catalog.nitrogen(allocator, 11);
    defer plant_nitrogen.deinit();
    try checkCatalog(plant_nitrogen);
    var plant_phosphorus = try plant_output_catalog.phosphorus(allocator, 11);
    defer plant_phosphorus.deinit();
    try checkCatalog(plant_phosphorus);
    var plant_heat = try plant_output_catalog.heat(allocator);
    defer plant_heat.deinit();
    try checkCatalog(plant_heat);
    var plant_daily_carbon = try plant_output_catalog.dailyCarbon(allocator, 11);
    defer plant_daily_carbon.deinit();
    try checkCatalog(plant_daily_carbon);
    var plant_daily_water = try plant_output_catalog.dailyWater(allocator);
    defer plant_daily_water.deinit();
    try checkCatalog(plant_daily_water);
    var plant_daily_nitrogen = try plant_output_catalog.dailyNitrogen(allocator);
    defer plant_daily_nitrogen.deinit();
    try checkCatalog(plant_daily_nitrogen);
    var plant_daily_phosphorus = try plant_output_catalog.dailyPhosphorus(allocator);
    defer plant_daily_phosphorus.deinit();
    try checkCatalog(plant_daily_phosphorus);
    var plant_daily_development = try plant_output_catalog.dailyDevelopment(allocator);
    defer plant_daily_development.deinit();
    try checkCatalog(plant_daily_development);
}

test "unit checker rejects underscore-per spellings and mis-cased prefixes" {
    try checkUnit("g m-3");
    try checkUnit("g N m-2 h-1");
    try checkUnit("umol mol-1");
    try checkUnit("fraction");
    // Correctly capitalized SI prefixes belong in a heading.
    try checkUnit("MJ m-2");
    try checkUnit("Mg m-3");
    try std.testing.expectError(error.NonSiOutputUnit, checkUnit("g_per_m3"));
    try std.testing.expectError(error.NonSiOutputUnit, checkUnit("W_per_m2"));
    try std.testing.expectError(error.AmbiguousOutputUnit, checkUnit("mg m-3"));
    try std.testing.expectError(error.AmbiguousOutputUnit, checkUnit("mj m-2 h-1"));
    try std.testing.expectError(error.EmptyOutputUnit, checkUnit(""));
}
