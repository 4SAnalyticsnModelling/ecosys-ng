const std = @import("std");

pub const Destination = enum { none, snowpack, soil_surface };

pub const Routing = struct {
    snowfall_m3_per_step: f64,
    rainfall_m3_per_step: f64,
    irrigation_m3_per_step: f64,
    snowpack_heat_capacity_mj_per_k: f64,
    minimum_snowpack_heat_capacity_mj_per_k: f64,
};

pub const Flux = struct {
    destination: Destination,
    amount_mol: []f64,

    pub fn deinit(self: *Flux, allocator: std.mem.Allocator) void {
        allocator.free(self.amount_mol);
        self.* = undefined;
    }
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    snowpack_amount_mol: []f64,
    soil_surface_amount_mol: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize) !State {
        if (cell_count == 0 or species_count == 0) return error.ZeroBoundaryStateExtent;
        const count = try std.math.mul(usize, cell_count, species_count);
        const snow = try allocator.alloc(f64, count);
        errdefer allocator.free(snow);
        const soil = try allocator.alloc(f64, count);
        errdefer allocator.free(soil);
        @memset(snow, 0);
        @memset(soil, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .species_count = species_count, .snowpack_amount_mol = snow, .soil_surface_amount_mol = soil };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.soil_surface_amount_mol);
        self.allocator.free(self.snowpack_amount_mol);
        self.* = undefined;
    }

    pub fn commit(self: *State, cell_index: usize, flux: Flux) !void {
        if (cell_index >= self.cell_count) return error.BoundaryCellIndexOutOfBounds;
        if (flux.amount_mol.len != self.species_count) return error.TransportSpeciesCountMismatch;
        const start = cell_index * self.species_count;
        const destination: ?[]f64 = switch (flux.destination) {
            .none => null,
            .snowpack => self.snowpack_amount_mol[start .. start + self.species_count],
            .soil_surface => self.soil_surface_amount_mol[start .. start + self.species_count],
        };
        for (flux.amount_mol) |amount| if (!std.math.isFinite(amount) or amount < 0) return error.InvalidSoluteBoundaryFlux;
        if (destination) |values| {
            for (values, flux.amount_mol) |value, amount| if (!std.math.isFinite(value + amount)) return error.NonFiniteSoluteBoundaryState;
            for (values, flux.amount_mol) |*value, amount| value.* += amount;
        }
    }
};

/// Exact TRNSFRS precipitation/irrigation routing. Rain enters an existing
/// snowpack; otherwise rain and irrigation enter the soil surface. Snowfall
/// always enters the snowpack.
pub fn calculate(allocator: std.mem.Allocator, routing: Routing, precipitation_concentration_mol_per_m3: []const f64, irrigation_concentration_mol_per_m3: []const f64) !Flux {
    if (precipitation_concentration_mol_per_m3.len == 0 or precipitation_concentration_mol_per_m3.len != irrigation_concentration_mol_per_m3.len) return error.TransportSpeciesCountMismatch;
    try validateRouting(routing);
    const snowpack_present = routing.snowpack_heat_capacity_mj_per_k > routing.minimum_snowpack_heat_capacity_mj_per_k;
    const destination: Destination = if (routing.snowfall_m3_per_step > 0 or (routing.rainfall_m3_per_step > 0 and snowpack_present))
        .snowpack
    else if ((routing.rainfall_m3_per_step > 0 or routing.irrigation_m3_per_step > 0) and !snowpack_present)
        .soil_surface
    else
        .none;
    const output = try allocator.alloc(f64, precipitation_concentration_mol_per_m3.len);
    errdefer allocator.free(output);
    const precipitation_water = routing.snowfall_m3_per_step + routing.rainfall_m3_per_step;
    for (output, precipitation_concentration_mol_per_m3, irrigation_concentration_mol_per_m3) |*amount, precipitation, irrigation| {
        if (!std.math.isFinite(precipitation) or precipitation < 0 or !std.math.isFinite(irrigation) or irrigation < 0) return error.InvalidBoundaryConcentration;
        amount.* = if (destination == .none) 0 else precipitation_water * precipitation + routing.irrigation_m3_per_step * irrigation;
        if (!std.math.isFinite(amount.*)) return error.NonFiniteSoluteBoundaryFlux;
    }
    return .{ .destination = destination, .amount_mol = output };
}

fn validateRouting(routing: Routing) !void {
    inline for (@typeInfo(Routing).@"struct".fields) |field| if (!std.math.isFinite(@field(routing, field.name)) or @field(routing, field.name) < 0) return error.InvalidSoluteBoundaryRouting;
}

test "rain routes through existing snowpack" {
    var flux = try calculate(std.testing.allocator, .{ .snowfall_m3_per_step = 0, .rainfall_m3_per_step = 2, .irrigation_m3_per_step = 1, .snowpack_heat_capacity_mj_per_k = 2, .minimum_snowpack_heat_capacity_mj_per_k = 1 }, &[_]f64{ 3, 4 }, &[_]f64{ 5, 6 });
    defer flux.deinit(std.testing.allocator);
    try std.testing.expectEqual(Destination.snowpack, flux.destination);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 11, 14 }, flux.amount_mol);
}

test "rain and irrigation route to bare soil surface" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    var flux = try calculate(std.testing.allocator, .{ .snowfall_m3_per_step = 0, .rainfall_m3_per_step = 2, .irrigation_m3_per_step = 1, .snowpack_heat_capacity_mj_per_k = 0, .minimum_snowpack_heat_capacity_mj_per_k = 1 }, &[_]f64{ 3, 4 }, &[_]f64{ 5, 6 });
    defer flux.deinit(std.testing.allocator);
    try std.testing.expectEqual(Destination.soil_surface, flux.destination);
    try state.commit(1, flux);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 0, 0, 11, 14 }, state.soil_surface_amount_mol);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 0, 0, 0, 0 }, state.snowpack_amount_mol);
}
