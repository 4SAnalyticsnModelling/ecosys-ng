const std = @import("std");
const topography_module = @import("topography.zig");

pub const Entry = struct {
    name: []const u8,
    topography: topography_module.Topography,
};

/// Owns each distinct topography file once. A per-cell map may point many
/// cells at the same entry or select a different entry for every cell.
pub const Catalog = struct {
    allocator: std.mem.Allocator,
    entries: std.ArrayList(Entry) = .empty,

    pub fn init(allocator: std.mem.Allocator) Catalog {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *Catalog) void {
        for (self.entries.items) |*entry| {
            entry.topography.deinit();
            self.allocator.free(entry.name);
        }
        self.entries.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn find(self: Catalog, name: []const u8) ?usize {
        for (self.entries.items, 0..) |entry, index| {
            if (std.mem.eql(u8, entry.name, name)) return index;
        }
        return null;
    }

    pub fn appendFromSource(self: *Catalog, name: []const u8, source: []const u8) !usize {
        try validateInputName(name);
        if (self.find(name)) |index| return index;
        const owned_name = try self.allocator.dupe(u8, name);
        errdefer self.allocator.free(owned_name);
        var parsed_topography = try topography_module.parse(self.allocator, source);
        errdefer parsed_topography.deinit();
        const index = self.entries.items.len;
        try self.entries.append(self.allocator, .{
            .name = owned_name,
            .topography = parsed_topography,
        });
        return index;
    }
};

fn validateInputName(name: []const u8) !void {
    if (name.len == 0 or
        std.mem.trim(u8, name, " \t\r\n").len != name.len or
        std.mem.indexOfAny(u8, name, " \t\r\n,|#") != null)
        return error.InvalidTopographyInputName;
}

const test_topography_source = "1 1 1 1 180 2 0 0\nsoil_profile\n";

test "topography catalog caches repeated filenames while retaining distinct files" {
    var catalog = Catalog.init(std.testing.allocator);
    defer catalog.deinit();

    const first = try catalog.appendFromSource("shared_topography", test_topography_source);
    const repeated = try catalog.appendFromSource("shared_topography", "not parsed");
    const distinct = try catalog.appendFromSource("other_topography", test_topography_source);

    try std.testing.expectEqual(first, repeated);
    try std.testing.expect(first != distinct);
    try std.testing.expectEqual(@as(usize, 2), catalog.entries.items.len);
}

test "topography catalog rejects malformed compulsory input names" {
    inline for (.{
        "",
        " ",
        " topography_a",
        "topography_a ",
        "topography,a",
        "topography|a",
        "topography\ta",
        "topography#a",
    }) |name| {
        var catalog = Catalog.init(std.testing.allocator);
        defer catalog.deinit();
        try std.testing.expectError(
            error.InvalidTopographyInputName,
            catalog.appendFromSource(name, test_topography_source),
        );
        try std.testing.expectEqual(@as(usize, 0), catalog.entries.items.len);
    }
}

test "topography catalog accepts extensionless and conventional file names" {
    inline for (.{ "topography_a", "topography-a.txt", "topography.a.csv" }) |name| {
        var catalog = Catalog.init(std.testing.allocator);
        defer catalog.deinit();
        _ = try catalog.appendFromSource(name, test_topography_source);
        try std.testing.expectEqualStrings(name, catalog.entries.items[0].name);
    }
}
