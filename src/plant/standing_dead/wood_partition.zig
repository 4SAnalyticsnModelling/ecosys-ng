const std = @import("std");

pub const litter_position_count: usize = 2;

pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    carbon_fraction_by_plant_and_position: []f64,
    nitrogen_fraction_by_plant_and_position: []f64,
    phosphorus_fraction_by_plant_and_position: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize) !State {
        if (plant_count == 0)
            return error.InvalidStandingDeadWoodPartitionDimensions;
        const element_count = try std.math.mul(
            usize,
            plant_count,
            litter_position_count,
        );
        const values = try allocator.alloc(
            f64,
            try std.math.mul(usize, element_count, 3),
        );
        @memset(values, 0);
        return .{
            .allocator = allocator,
            .plant_count = plant_count,
            .carbon_fraction_by_plant_and_position = values[0..element_count],
            .nitrogen_fraction_by_plant_and_position = values[element_count .. 2 * element_count],
            .phosphorus_fraction_by_plant_and_position = values[2 * element_count ..],
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(
            self.carbon_fraction_by_plant_and_position.ptr[0 .. 3 * self.plant_count * litter_position_count],
        );
        self.* = undefined;
    }
};

pub const Inputs = struct {
    biomass_turnover_type_by_plant: []const u8,
    root_profile_type_by_plant: []const u8,
    structural_presence_threshold_g_c_by_plant: []const f64,
    plant_branch_offsets: []const usize,
    stalk_carbon_g_c_by_branch: []const f64,
    sapwood_carbon_g_c_by_branch: []const f64,
};

/// grosub.f lines 422–434, 453–454, and 461–462. Position zero is woody and
/// position one is nonwoody. N and P copy the C fractions exactly.
pub fn refresh(state: *State, inputs: Inputs) !void {
    if (inputs.biomass_turnover_type_by_plant.len != state.plant_count or
        inputs.root_profile_type_by_plant.len != state.plant_count or
        inputs.structural_presence_threshold_g_c_by_plant.len != state.plant_count or
        inputs.plant_branch_offsets.len != state.plant_count + 1 or
        inputs.plant_branch_offsets[0] != 0)
        return error.InvalidStandingDeadWoodPartitionDimensions;
    const branch_count = inputs.plant_branch_offsets[state.plant_count];
    if (inputs.stalk_carbon_g_c_by_branch.len != branch_count or
        inputs.sapwood_carbon_g_c_by_branch.len != branch_count)
        return error.InvalidStandingDeadWoodPartitionDimensions;
    try validateOffsets(inputs.plant_branch_offsets, branch_count);

    for (0..state.plant_count) |plant| _ = try fractionsForPlant(inputs, plant);

    for (0..state.plant_count) |plant| {
        const fractions = fractionsForPlant(inputs, plant) catch unreachable;
        const first = plant * litter_position_count;
        state.carbon_fraction_by_plant_and_position[first..][0..2].* = fractions;
        state.nitrogen_fraction_by_plant_and_position[first..][0..2].* = fractions;
        state.phosphorus_fraction_by_plant_and_position[first..][0..2].* = fractions;
    }
}

fn fractionsForPlant(inputs: Inputs, plant: usize) ![2]f64 {
    const threshold = inputs.structural_presence_threshold_g_c_by_plant[plant];
    if (!std.math.isFinite(threshold) or threshold < 0)
        return error.InvalidStandingDeadWoodPartitionInput;
    var stalk_carbon_g_c: f64 = 0;
    var sapwood_carbon_g_c: f64 = 0;
    for (inputs.plant_branch_offsets[plant]..inputs.plant_branch_offsets[plant + 1]) |branch| {
        const stalk = inputs.stalk_carbon_g_c_by_branch[branch];
        const sapwood = inputs.sapwood_carbon_g_c_by_branch[branch];
        if (!std.math.isFinite(stalk) or stalk < 0 or
            !std.math.isFinite(sapwood) or sapwood < 0)
            return error.InvalidStandingDeadWoodPartitionInput;
        stalk_carbon_g_c += stalk;
        sapwood_carbon_g_c += sapwood;
        if (!std.math.isFinite(stalk_carbon_g_c) or
            !std.math.isFinite(sapwood_carbon_g_c))
            return error.NonFiniteStandingDeadWoodPartition;
    }
    if (sapwood_carbon_g_c > stalk_carbon_g_c)
        return error.SapwoodExceedsStalkCarbon;
    const nonwoody_fraction =
        if (inputs.biomass_turnover_type_by_plant[plant] == 0 or
        inputs.root_profile_type_by_plant[plant] <= 1 or
        stalk_carbon_g_c <= threshold)
            1.0
        else
            sapwood_carbon_g_c / stalk_carbon_g_c;
    if (!std.math.isFinite(nonwoody_fraction) or
        nonwoody_fraction < 0 or nonwoody_fraction > 1)
        return error.NonFiniteStandingDeadWoodPartition;
    return .{ 1.0 - nonwoody_fraction, nonwoody_fraction };
}

fn validateOffsets(offsets: []const usize, branch_count: usize) !void {
    var previous: usize = 0;
    for (offsets) |offset| {
        if (offset < previous or offset > branch_count)
            return error.InvalidStandingDeadWoodPartitionTopology;
        previous = offset;
    }
}

test "dynamic FWOOD preserves branch order and copies C fractions to N and P" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try refresh(&state, .{
        .biomass_turnover_type_by_plant = &.{ 2, 2 },
        .root_profile_type_by_plant = &.{ 2, 3 },
        .structural_presence_threshold_g_c_by_plant = &.{ 0, 0 },
        .plant_branch_offsets = &.{ 0, 2, 3 },
        .stalk_carbon_g_c_by_branch = &.{ 3, 7, 8 },
        .sapwood_carbon_g_c_by_branch = &.{ 1, 3, 2 },
    });
    const expected = [_]f64{ 0.6, 0.4, 0.75, 0.25 };
    try std.testing.expectEqualSlices(
        f64,
        &expected,
        state.carbon_fraction_by_plant_and_position,
    );
    try std.testing.expectEqualSlices(
        f64,
        &expected,
        state.nitrogen_fraction_by_plant_and_position,
    );
    try std.testing.expectEqualSlices(
        f64,
        &expected,
        state.phosphorus_fraction_by_plant_and_position,
    );
}

test "turnover root profile and threshold branches force nonwoody" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    try refresh(&state, .{
        .biomass_turnover_type_by_plant = &.{ 0, 2, 2 },
        .root_profile_type_by_plant = &.{ 3, 1, 3 },
        .structural_presence_threshold_g_c_by_plant = &.{ 0, 0, 2 },
        .plant_branch_offsets = &.{ 0, 1, 2, 3 },
        .stalk_carbon_g_c_by_branch = &.{ 10, 10, 2 },
        .sapwood_carbon_g_c_by_branch = &.{ 2, 2, 1 },
    });
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 1, 0, 1, 0, 1 },
        state.carbon_fraction_by_plant_and_position,
    );
}

test "late invalid plant preserves all element publications" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    @memset(state.carbon_fraction_by_plant_and_position, 7);
    @memset(state.nitrogen_fraction_by_plant_and_position, 7);
    @memset(state.phosphorus_fraction_by_plant_and_position, 7);
    try std.testing.expectError(
        error.SapwoodExceedsStalkCarbon,
        refresh(&state, .{
            .biomass_turnover_type_by_plant = &.{ 2, 2 },
            .root_profile_type_by_plant = &.{ 2, 2 },
            .structural_presence_threshold_g_c_by_plant = &.{ 0, 0 },
            .plant_branch_offsets = &.{ 0, 1, 2 },
            .stalk_carbon_g_c_by_branch = &.{ 10, 1 },
            .sapwood_carbon_g_c_by_branch = &.{ 5, 2 },
        }),
    );
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            for (@field(state, field.name)) |value|
                try std.testing.expectEqual(@as(f64, 7), value);
}
