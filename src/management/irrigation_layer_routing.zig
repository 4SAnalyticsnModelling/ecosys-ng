const std = @import("std");
const irrigation = @import("irrigation_schedule.zig");

pub const dissolved_species_count = 11;

/// Heap-owned hourly irrigation loads. Surface carriers remain indexed by
/// cell; subsurface carriers use the runtime cell×soil-layer layout. All
/// chemistry is extensive, so simultaneous events at different depths retain
/// their own destination and can be committed without concentration averaging.
pub const Loads = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    soil_layer_capacity: usize,
    surface_water_m3: []f64,
    subsurface_water_m3: []f64,
    surface_dissolved_mass_g: []f64,
    subsurface_dissolved_mass_g: []f64,
    surface_hydrogen_mol: []f64,
    subsurface_hydrogen_mol: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        soil_layer_capacity: usize,
    ) !Loads {
        if (cell_count == 0 or soil_layer_capacity == 0)
            return error.EmptyIrrigationRoutingDomain;
        const layer_count = try std.math.mul(
            usize,
            cell_count,
            soil_layer_capacity,
        );
        const surface_species_count = try std.math.mul(
            usize,
            cell_count,
            dissolved_species_count,
        );
        const subsurface_species_count = try std.math.mul(
            usize,
            layer_count,
            dissolved_species_count,
        );
        var result: Loads = .{
            .allocator = allocator,
            .cell_count = cell_count,
            .soil_layer_capacity = soil_layer_capacity,
            .surface_water_m3 = try allocator.alloc(f64, cell_count),
            .subsurface_water_m3 = undefined,
            .surface_dissolved_mass_g = undefined,
            .subsurface_dissolved_mass_g = undefined,
            .surface_hydrogen_mol = undefined,
            .subsurface_hydrogen_mol = undefined,
        };
        errdefer allocator.free(result.surface_water_m3);
        result.subsurface_water_m3 = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(result.subsurface_water_m3);
        result.surface_dissolved_mass_g =
            try allocator.alloc(f64, surface_species_count);
        errdefer allocator.free(result.surface_dissolved_mass_g);
        result.subsurface_dissolved_mass_g =
            try allocator.alloc(f64, subsurface_species_count);
        errdefer allocator.free(result.subsurface_dissolved_mass_g);
        result.surface_hydrogen_mol = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.surface_hydrogen_mol);
        result.subsurface_hydrogen_mol =
            try allocator.alloc(f64, layer_count);
        result.reset();
        return result;
    }

    pub fn deinit(self: *Loads) void {
        self.allocator.free(self.subsurface_hydrogen_mol);
        self.allocator.free(self.surface_hydrogen_mol);
        self.allocator.free(self.subsurface_dissolved_mass_g);
        self.allocator.free(self.surface_dissolved_mass_g);
        self.allocator.free(self.subsurface_water_m3);
        self.allocator.free(self.surface_water_m3);
        self.* = undefined;
    }

    pub fn reset(self: *Loads) void {
        @memset(self.surface_water_m3, 0);
        @memset(self.subsurface_water_m3, 0);
        @memset(self.surface_dissolved_mass_g, 0);
        @memset(self.subsurface_dissolved_mass_g, 0);
        @memset(self.surface_hydrogen_mol, 0);
        @memset(self.subsurface_hydrogen_mol, 0);
    }

    pub fn accumulate(
        self: *Loads,
        cell: usize,
        active_layer_count: usize,
        layer_thickness_m: []const f64,
        cell_area_m2: f64,
        water_depth_m: f64,
        application_depth_m: f64,
        chemistry: irrigation.WaterChemistry_g_per_m3,
    ) !void {
        if (cell >= self.cell_count or
            layer_thickness_m.len != self.soil_layer_capacity)
            return error.IrrigationRoutingDimensionMismatch;
        if (active_layer_count == 0 or
            active_layer_count > self.soil_layer_capacity)
            return error.InvalidActiveIrrigationLayerCount;
        inline for (.{ cell_area_m2, water_depth_m, application_depth_m }) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteIrrigationRoutingInput;
        if (cell_area_m2 <= 0 or water_depth_m < 0 or application_depth_m < 0)
            return error.InvalidIrrigationRoutingInput;
        const concentrations = chemistryValues(chemistry);
        for (concentrations) |concentration|
            if (!std.math.isFinite(concentration) or concentration < 0)
                return error.InvalidIrrigationChemistry;
        if (!std.math.isFinite(chemistry.ph) or
            chemistry.ph < 0 or chemistry.ph > 14)
            return error.InvalidIrrigationChemistry;
        var cumulative_depth_m: f64 = 0;
        for (layer_thickness_m[0..active_layer_count]) |thickness_m| {
            if (!std.math.isFinite(thickness_m) or thickness_m <= 0)
                return error.InvalidIrrigationLayerThickness;
            cumulative_depth_m += thickness_m;
            if (!std.math.isFinite(cumulative_depth_m))
                return error.IrrigationRoutingOverflow;
        }
        if (water_depth_m == 0) return;
        const water_volume_m3 = water_depth_m * cell_area_m2;
        if (!std.math.isFinite(water_volume_m3))
            return error.IrrigationRoutingOverflow;
        const hydrogen_mol =
            water_volume_m3 * 1000.0 * std.math.pow(f64, 10.0, -chemistry.ph);
        if (!std.math.isFinite(hydrogen_mol))
            return error.IrrigationRoutingOverflow;

        if (application_depth_m == 0) {
            const next_water = self.surface_water_m3[cell] + water_volume_m3;
            const next_hydrogen =
                self.surface_hydrogen_mol[cell] + hydrogen_mol;
            if (!std.math.isFinite(next_water) or
                !std.math.isFinite(next_hydrogen))
                return error.IrrigationRoutingOverflow;
            const first = cell * dissolved_species_count;
            var candidate = [_]f64{0} ** dissolved_species_count;
            for (&candidate, concentrations, 0..) |*next, concentration, species| {
                next.* = self.surface_dissolved_mass_g[first + species] +
                    water_volume_m3 * concentration;
                if (!std.math.isFinite(next.*))
                    return error.IrrigationRoutingOverflow;
            }
            self.surface_water_m3[cell] = next_water;
            self.surface_hydrogen_mol[cell] = next_hydrogen;
            @memcpy(
                self.surface_dissolved_mass_g[first..][0..dissolved_species_count],
                &candidate,
            );
            return;
        }

        const local_layer = try layerAtDepth(
            layer_thickness_m[0..active_layer_count],
            application_depth_m,
        );
        const layer = cell * self.soil_layer_capacity + local_layer;
        const next_water = self.subsurface_water_m3[layer] + water_volume_m3;
        const next_hydrogen =
            self.subsurface_hydrogen_mol[layer] + hydrogen_mol;
        if (!std.math.isFinite(next_water) or !std.math.isFinite(next_hydrogen))
            return error.IrrigationRoutingOverflow;
        const first = layer * dissolved_species_count;
        var candidate = [_]f64{0} ** dissolved_species_count;
        for (&candidate, concentrations, 0..) |*next, concentration, species| {
            next.* = self.subsurface_dissolved_mass_g[first + species] +
                water_volume_m3 * concentration;
            if (!std.math.isFinite(next.*))
                return error.IrrigationRoutingOverflow;
        }
        self.subsurface_water_m3[layer] = next_water;
        self.subsurface_hydrogen_mol[layer] = next_hydrogen;
        @memcpy(
            self.subsurface_dissolved_mass_g[first..][0..dissolved_species_count],
            &candidate,
        );
    }
};

pub fn layerAtDepth(
    active_layer_thickness_m: []const f64,
    application_depth_m: f64,
) !usize {
    if (active_layer_thickness_m.len == 0)
        return error.EmptyActiveIrrigationProfile;
    if (!std.math.isFinite(application_depth_m) or application_depth_m <= 0)
        return error.InvalidSubsurfaceIrrigationDepth;

    var cumulative_depth_m: f64 = 0;
    for (active_layer_thickness_m, 0..) |thickness_m, layer| {
        if (!std.math.isFinite(thickness_m) or thickness_m <= 0)
            return error.InvalidIrrigationLayerThickness;
        cumulative_depth_m += thickness_m;
        if (!std.math.isFinite(cumulative_depth_m))
            return error.IrrigationRoutingOverflow;
        if (cumulative_depth_m >= application_depth_m) return layer;
    }
    return error.SubsurfaceIrrigationDepthBelowProfile;
}

fn chemistryValues(
    water: irrigation.WaterChemistry_g_per_m3,
) [dissolved_species_count]f64 {
    return .{
        water.ammonium_nitrogen,
        water.nitrate_nitrogen,
        water.phosphate_phosphorus,
        water.aluminum,
        water.iron,
        water.calcium,
        water.magnesium,
        water.sodium,
        water.potassium,
        water.sulfate_sulfur,
        water.chloride,
    };
}

test "runtime irrigation depths retain separate surface and soil-layer loads" {
    var loads = try Loads.init(std.testing.allocator, 1, 7);
    defer loads.deinit();
    const chemistry: irrigation.WaterChemistry_g_per_m3 = .{
        .ph = 7,
        .ammonium_nitrogen = 1,
        .nitrate_nitrogen = 2,
        .phosphate_phosphorus = 3,
        .aluminum = 4,
        .iron = 5,
        .calcium = 6,
        .magnesium = 7,
        .sodium = 8,
        .potassium = 9,
        .sulfate_sulfur = 10,
        .chloride = 11,
    };
    const thickness = [_]f64{ 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35 };
    try loads.accumulate(0, 7, &thickness, 10, 0.001, 0, chemistry);
    try loads.accumulate(0, 7, &thickness, 10, 0.002, 0.12, chemistry);
    try loads.accumulate(0, 7, &thickness, 10, 0.003, 0.90, chemistry);
    try std.testing.expectEqual(@as(f64, 0.01), loads.surface_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 0.02), loads.subsurface_water_m3[1]);
    try std.testing.expectEqual(@as(f64, 0.03), loads.subsurface_water_m3[5]);
    try std.testing.expectEqual(
        @as(f64, 0.04),
        loads.surface_dissolved_mass_g[3],
    );
    try std.testing.expectEqual(
        @as(f64, 0.10),
        loads.subsurface_dissolved_mass_g[
            1 * dissolved_species_count + 4
        ],
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.33),
        loads.subsurface_dissolved_mass_g[
            5 * dissolved_species_count + 10
        ],
        1.0e-14,
    );
}

test "subsurface irrigation depth below profile fails without mutation" {
    var loads = try Loads.init(std.testing.allocator, 2, 3);
    defer loads.deinit();
    const chemistry = std.mem.zeroInit(
        irrigation.WaterChemistry_g_per_m3,
        .{ .ph = 6 },
    );
    try std.testing.expectError(
        error.SubsurfaceIrrigationDepthBelowProfile,
        loads.accumulate(
            1,
            2,
            &.{ 0.1, 0.2, 9.0 },
            5,
            0.004,
            5,
            chemistry,
        ),
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        loads.subsurface_water_m3[1 * 3 + 1],
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        loads.subsurface_hydrogen_mol[1 * 3 + 1],
    );
}

test "subsurface irrigation layer scan uses active cumulative bottoms" {
    try std.testing.expectEqual(
        @as(usize, 0),
        try layerAtDepth(&.{ 0.05, 0.10, 0.20 }, 0.05),
    );
    try std.testing.expectEqual(
        @as(usize, 1),
        try layerAtDepth(&.{ 0.05, 0.10, 0.20 }, 0.0500001),
    );
    try std.testing.expectEqual(
        @as(usize, 2),
        try layerAtDepth(&.{ 0.05, 0.10, 0.20 }, 0.35),
    );
    try std.testing.expectError(
        error.SubsurfaceIrrigationDepthBelowProfile,
        layerAtDepth(&.{ 0.05, 0.10, 0.20 }, 0.350001),
    );
}

test "failed irrigation event cannot partially alter its destination" {
    var loads = try Loads.init(std.testing.allocator, 1, 1);
    defer loads.deinit();
    var chemistry = std.mem.zeroInit(
        irrigation.WaterChemistry_g_per_m3,
        .{ .ph = 7 },
    );
    try loads.accumulate(0, 1, &.{0.1}, 10, 0.001, 0, chemistry);
    const before_water = loads.surface_water_m3[0];
    const before_hydrogen = loads.surface_hydrogen_mol[0];
    chemistry.nitrate_nitrogen = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidIrrigationChemistry,
        loads.accumulate(0, 1, &.{0.1}, 10, 0.002, 0, chemistry),
    );
    try std.testing.expectEqual(before_water, loads.surface_water_m3[0]);
    try std.testing.expectEqual(
        before_hydrogen,
        loads.surface_hydrogen_mol[0],
    );
}
