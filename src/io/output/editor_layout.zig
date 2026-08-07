const std = @import("std");

pub const source_domain_choice_count: usize = 50;
pub const source_total_choice_count: usize = 2 * source_domain_choice_count;

pub const Resolved = struct {
    allocator: std.mem.Allocator,
    soil_enabled: []bool,
    plant_enabled: []bool,

    pub fn deinit(self: *Resolved) void {
        self.allocator.free(self.plant_enabled);
        self.allocator.free(self.soil_enabled);
        self.* = undefined;
    }
};

/// Resolves the source FOUTS/FOUTP 50+50 editor layout or a modern compact
/// soil-count + plant-count layout. Runtime variables beyond the historical
/// fifty slots default off in a source editor and remain explicitly selectable
/// in a compact editor.
pub fn resolve(allocator: std.mem.Allocator, choices: []const bool, soil_variable_count: usize, plant_variable_count: usize) !Resolved {
    if (soil_variable_count == 0 or plant_variable_count == 0) return error.EmptyOutputCatalog;
    const compact_count = try std.math.add(usize, soil_variable_count, plant_variable_count);
    if (choices.len != source_total_choice_count and choices.len != compact_count) return error.OutputEditorChoiceCountMismatch;
    const soil = try allocator.alloc(bool, soil_variable_count);
    errdefer allocator.free(soil);
    const plant = try allocator.alloc(bool, plant_variable_count);
    errdefer allocator.free(plant);
    @memset(soil, false);
    @memset(plant, false);
    if (choices.len == source_total_choice_count) {
        @memcpy(soil[0..@min(soil.len, source_domain_choice_count)], choices[0..@min(soil.len, source_domain_choice_count)]);
        @memcpy(plant[0..@min(plant.len, source_domain_choice_count)], choices[source_domain_choice_count..][0..@min(plant.len, source_domain_choice_count)]);
    } else {
        @memcpy(soil, choices[0..soil.len]);
        @memcpy(plant, choices[soil.len..]);
    }
    return .{ .allocator = allocator, .soil_enabled = soil, .plant_enabled = plant };
}

test "source editor maps independent fifty-choice soil and plant domains" {
    var choices = [_]bool{false} ** source_total_choice_count;
    choices[2] = true;
    choices[source_domain_choice_count + 3] = true;
    var resolved = try resolve(std.testing.allocator, &choices, 60, 7);
    defer resolved.deinit();
    try std.testing.expect(resolved.soil_enabled[2]);
    try std.testing.expect(!resolved.soil_enabled[52]);
    try std.testing.expect(resolved.plant_enabled[3]);
}

test "compact editor selects every runtime-expanded variable" {
    const choices = [_]bool{ true, false, true, true, false };
    var resolved = try resolve(std.testing.allocator, &choices, 3, 2);
    defer resolved.deinit();
    try std.testing.expectEqualSlices(bool, &.{ true, false, true }, resolved.soil_enabled);
    try std.testing.expectEqualSlices(bool, &.{ true, false }, resolved.plant_enabled);
}
