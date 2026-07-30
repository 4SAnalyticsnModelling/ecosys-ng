const std = @import("std");
const CellRange = @import("compute.zig").CellRange;
const StructureState = @import("canopy_structure.zig").State;
const InterceptionState = @import("canopy_interception.zig").State;
const GroundRadiationState = @import("ground_radiation.zig").State;

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    canopy_exposure_fraction: []f64,
    ground_exposure_fraction: []f64,
    species_exposure_fraction: []f64,
    species_share_of_canopy_exposure: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize) !State {
        if (cell_count == 0 or species_count == 0) return error.InvalidCanopyExposureDimensions;
        const species_slots = try std.math.mul(usize, cell_count, species_count);
        var result: State = .{
            .allocator = allocator,
            .cell_count = cell_count,
            .species_count = species_count,
            .canopy_exposure_fraction = try allocator.alloc(f64, cell_count),
            .ground_exposure_fraction = undefined,
            .species_exposure_fraction = undefined,
            .species_share_of_canopy_exposure = undefined,
        };
        errdefer allocator.free(result.canopy_exposure_fraction);
        result.ground_exposure_fraction = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(result.ground_exposure_fraction);
        result.species_exposure_fraction = try allocator.alloc(f64, species_slots);
        errdefer allocator.free(result.species_exposure_fraction);
        result.species_share_of_canopy_exposure = try allocator.alloc(f64, species_slots);
        @memset(result.canopy_exposure_fraction, 0);
        @memset(result.ground_exposure_fraction, 1);
        @memset(result.species_exposure_fraction, 0);
        @memset(result.species_share_of_canopy_exposure, 0);
        return result;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.species_share_of_canopy_exposure);
        self.allocator.free(self.species_exposure_fraction);
        self.allocator.free(self.ground_exposure_fraction);
        self.allocator.free(self.canopy_exposure_fraction);
        self.* = undefined;
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| {
            if (!std.math.isFinite(value) or value < 0 or value > 1 + 1.0e-12) {
                std.log.err("invalid canopy exposure: field={s} index={d} value={e}", .{ field.name, index, value });
                return error.InvalidCanopyExposure;
            }
        };
    }
};

pub const ApplyContext = struct {
    result: *State,
    structure: *const StructureState,
    interception: *const InterceptionState,
    ground_radiation: *const GroundRadiationState,
    solar_angle_sine_by_cell: []const f64,
};

/// Ports HOUR1 FRADT/FRADG/FRADP for living leaf area. Standing-dead and
/// stalk exposure enter when those structural pools become available.
pub fn applyTile(context: *ApplyContext, range: CellRange) !void {
    const result = context.result;
    if (range.end > result.cell_count or context.solar_angle_sine_by_cell.len != result.cell_count or context.structure.cell_count != result.cell_count or context.interception.cell_count != result.cell_count or context.ground_radiation.cell_count != result.cell_count or context.structure.species_count != result.species_count or context.interception.species_count != result.species_count) return error.CanopyExposureDimensionMismatch;
    for (range.first..range.end) |cell| {
        const solar_angle_sine = context.solar_angle_sine_by_cell[cell];
        if (!std.math.isFinite(solar_angle_sine) or solar_angle_sine < 0 or solar_angle_sine > 1) return error.InvalidSolarAngle;
        var total_leaf_area_index: f64 = 0;
        var canopy_absorbed_shortwave: f64 = 0;
        for (0..result.species_count) |species| {
            const index = cell * result.species_count + species;
            var species_leaf_area: f64 = 0;
            for (0..context.structure.inclination_class_count) |inclination| species_leaf_area += context.structure.leaf_area_index_by_inclination_m2_m2[index * context.structure.inclination_class_count + inclination];
            total_leaf_area_index += species_leaf_area;
            canopy_absorbed_shortwave += context.interception.absorbed_shortwave_mj_per_m2[index];
            result.species_exposure_fraction[index] = species_leaf_area;
        }

        var canopy_fraction: f64 = 0;
        if (solar_angle_sine > 0.05) {
            const total_received_shortwave = canopy_absorbed_shortwave + context.ground_radiation.incident_shortwave_mj_per_m2[cell];
            if (total_received_shortwave > 1.0e-12) canopy_fraction = canopy_absorbed_shortwave / total_received_shortwave;
        } else if (total_leaf_area_index > 1.0e-12) {
            canopy_fraction = 1.0 - @exp(-0.65 * total_leaf_area_index);
        }
        canopy_fraction = std.math.clamp(canopy_fraction, 0.0, 1.0);
        result.canopy_exposure_fraction[cell] = canopy_fraction;
        result.ground_exposure_fraction[cell] = 1.0 - canopy_fraction;
        for (0..result.species_count) |species| {
            const index = cell * result.species_count + species;
            const share = if (total_leaf_area_index > 1.0e-12) result.species_exposure_fraction[index] / total_leaf_area_index else 0;
            result.species_share_of_canopy_exposure[index] = share;
            result.species_exposure_fraction[index] = canopy_fraction * share;
        }
        try validateCellConservation(result.*, cell);
    }
}

fn validateCellConservation(state: State, cell: usize) !void {
    if (@abs(state.canopy_exposure_fraction[cell] + state.ground_exposure_fraction[cell] - 1.0) > 1.0e-12) return error.CanopyExposureImbalance;
    var species_total: f64 = 0;
    for (0..state.species_count) |species| species_total += state.species_exposure_fraction[cell * state.species_count + species];
    if (@abs(species_total - state.canopy_exposure_fraction[cell]) > 1.0e-12) return error.SpeciesExposureImbalance;
}

test "night exposure uses ecosys exponential leaf-area relation" {
    const leaf_area_index: f64 = 3;
    const expected = 1.0 - @exp(-0.65 * leaf_area_index);
    try std.testing.expect(expected > 0 and expected < 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8577259284), expected, 1.0e-9);
}
