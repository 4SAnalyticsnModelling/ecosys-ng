const std = @import("std");

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    process_unit_count_per_layer: usize,
    previous_total_non_band_h2po4_demand_g_p: []f64,
    previous_total_band_h2po4_demand_g_p: []f64,
    previous_total_non_band_hpo4_demand_g_p: []f64,
    previous_total_band_hpo4_demand_g_p: []f64,
    previous_non_band_h2po4_capacity_g_p: []f64,
    previous_band_h2po4_capacity_g_p: []f64,
    previous_non_band_hpo4_capacity_g_p: []f64,
    previous_band_hpo4_capacity_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize, process_unit_count_per_layer: usize) !State {
        if (layer_count == 0 or process_unit_count_per_layer == 0) return error.InvalidSoilMicrobialPhosphorusDimensions;
        var result: State = undefined;
        result.allocator = allocator;
        result.layer_count = layer_count;
        result.process_unit_count_per_layer = process_unit_count_per_layer;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const count = if (std.mem.startsWith(u8, field.name, "previous_total_")) layer_count else try std.math.mul(usize, layer_count, process_unit_count_per_layer);
            @field(result, field.name) = try allocator.alloc(f64, count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

test "phosphorus competition history is runtime sized" {
    var state = try State.init(std.testing.allocator, 3, 17);
    defer state.deinit();
    try std.testing.expectEqual(@as(usize, 3), state.previous_total_band_hpo4_demand_g_p.len);
    try std.testing.expectEqual(@as(usize, 51), state.previous_non_band_h2po4_capacity_g_p.len);
}
