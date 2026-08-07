const std = @import("std");

pub const Dimensions = struct {
    cell_count: usize,
    plant_species_per_cell: usize,
    layer_count: usize,
    litter_pool_count: usize = 5,
    litter_position_count: usize = 2,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    dimensions: Dimensions,
    layered_litter_carbon_g_c_per_h: []f64,
    layered_litter_nitrogen_g_n_per_h: []f64,
    layered_litter_phosphorus_g_p_per_h: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        dimensions: Dimensions,
    ) !State {
        if (dimensions.cell_count == 0 or
            dimensions.plant_species_per_cell == 0 or
            dimensions.layer_count == 0 or
            dimensions.litter_pool_count == 0 or
            dimensions.litter_position_count == 0)
            return error.InvalidLayeredLitterPublicationDimensions;
        const layered_count = try std.math.mul(
            usize,
            dimensions.cell_count,
            try perCellValueCount(dimensions),
        );
        const value_count = try std.math.mul(usize, layered_count, 3);
        const values = try allocator.alloc(f64, value_count);
        @memset(values, 0);
        return .{
            .allocator = allocator,
            .dimensions = dimensions,
            .layered_litter_carbon_g_c_per_h = values[0..layered_count],
            .layered_litter_nitrogen_g_n_per_h = values[layered_count .. 2 * layered_count],
            .layered_litter_phosphorus_g_p_per_h = values[2 * layered_count .. 3 * layered_count],
        };
    }

    pub fn deinit(self: *State) void {
        const layered_count = self.layered_litter_carbon_g_c_per_h.len;
        self.allocator.free(
            self.layered_litter_carbon_g_c_per_h.ptr[0 .. 3 * layered_count],
        );
        self.* = undefined;
    }
};

pub const Inputs = struct {
    layered_litter_carbon_g_c_per_h_by_plant: []const f64,
    layered_litter_nitrogen_g_n_per_h_by_plant: []const f64,
    layered_litter_phosphorus_g_p_per_h_by_plant: []const f64,
};

/// Exact EXTRACT lines 64--69. Runtime cell and plant dimensions replace
/// source `NY/NX/NZ`, while preserving source litter-loop order L, K, then M.
pub fn refresh(state: *State, inputs: Inputs) !void {
    const plant_count = try std.math.mul(
        usize,
        state.dimensions.cell_count,
        state.dimensions.plant_species_per_cell,
    );
    const per_plant = try perPlantValueCount(state.dimensions);
    const input_count = try std.math.mul(usize, plant_count, per_plant);
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (@field(inputs, field.name).len != input_count)
            return error.InvalidLayeredLitterPublicationDimensions;
    const output_count = state.layered_litter_carbon_g_c_per_h.len;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64 and
            @field(state.*, field.name).len != output_count)
            return error.InvalidLayeredLitterPublicationDimensions;

    for (0..state.dimensions.cell_count) |cell|
        for (0..per_plant) |index| {
            _ = try summedValues(state.dimensions, inputs, cell, index);
        };

    for (0..state.dimensions.cell_count) |cell| {
        const first = cell * per_plant;
        for (0..per_plant) |index| {
            const summed = summedValues(state.dimensions, inputs, cell, index) catch unreachable;
            state.layered_litter_carbon_g_c_per_h[first + index] = summed.carbon;
            state.layered_litter_nitrogen_g_n_per_h[first + index] = summed.nitrogen;
            state.layered_litter_phosphorus_g_p_per_h[first + index] = summed.phosphorus;
        }
    }
}

const SummedValues = struct {
    carbon: f64,
    nitrogen: f64,
    phosphorus: f64,
};

fn summedValues(
    dimensions: Dimensions,
    inputs: Inputs,
    cell: usize,
    index: usize,
) !SummedValues {
    const per_plant = try perPlantValueCount(dimensions);
    var totals: SummedValues = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 };
    const first_plant = cell * dimensions.plant_species_per_cell;
    for (0..dimensions.plant_species_per_cell) |local_plant| {
        const plant = first_plant + local_plant;
        const plant_first = plant * per_plant;
        const source_carbon =
            inputs.layered_litter_carbon_g_c_per_h_by_plant[plant_first + index];
        const source_nitrogen =
            inputs.layered_litter_nitrogen_g_n_per_h_by_plant[plant_first + index];
        const source_phosphorus =
            inputs.layered_litter_phosphorus_g_p_per_h_by_plant[plant_first + index];
        inline for (.{
            source_carbon,
            source_nitrogen,
            source_phosphorus,
        }) |value| if (!std.math.isFinite(value))
            return error.InvalidLayeredLitterPublicationInput;

        totals.carbon += source_carbon;
        totals.nitrogen += source_nitrogen;
        totals.phosphorus += source_phosphorus;
        inline for (.{
            totals.carbon,
            totals.nitrogen,
            totals.phosphorus,
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteLayeredLitterPublication;
    }
    return totals;
}

fn perCellValueCount(dimensions: Dimensions) !usize {
    return try std.math.mul(
        usize,
        dimensions.layer_count,
        try std.math.mul(
            usize,
            dimensions.litter_position_count,
            dimensions.litter_pool_count,
        ),
    );
}

fn perPlantValueCount(dimensions: Dimensions) !usize {
    return perCellValueCount(dimensions);
}

test "EXTRACT layered litter publication preserves L then K then M order" {
    const dimensions: Dimensions = .{
        .cell_count = 1,
        .plant_species_per_cell = 2,
        .layer_count = 2,
        .litter_pool_count = 3,
    };
    var state = try State.init(std.testing.allocator, dimensions);
    defer state.deinit();

    try refresh(&state, .{
        .layered_litter_carbon_g_c_per_h_by_plant = &.{
            1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,
            10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
        },
        .layered_litter_nitrogen_g_n_per_h_by_plant = &.{
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
            1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,
        },
        .layered_litter_phosphorus_g_p_per_h_by_plant = &.{
            0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
            0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,
        },
    });

    try std.testing.expectEqualSlices(
        f64,
        &.{ 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 121, 132 },
        state.layered_litter_carbon_g_c_per_h,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.1),
        state.layered_litter_nitrogen_g_n_per_h[0],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 13.2),
        state.layered_litter_nitrogen_g_n_per_h[11],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.32),
        state.layered_litter_phosphorus_g_p_per_h[11],
        1e-15,
    );
}

test "late invalid plant value preserves all publication arrays" {
    const dimensions: Dimensions = .{
        .cell_count = 1,
        .plant_species_per_cell = 2,
        .layer_count = 1,
    };
    var state = try State.init(std.testing.allocator, dimensions);
    defer state.deinit();
    @memset(state.layered_litter_carbon_g_c_per_h, 9);
    @memset(state.layered_litter_nitrogen_g_n_per_h, 9);
    @memset(state.layered_litter_phosphorus_g_p_per_h, 9);
    try std.testing.expectError(
        error.InvalidLayeredLitterPublicationInput,
        refresh(&state, .{
            .layered_litter_carbon_g_c_per_h_by_plant = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, std.math.nan(f64), 12, 13, 14, 15, 16, 17, 18, 19, 20 },
            .layered_litter_nitrogen_g_n_per_h_by_plant = &([_]f64{0} ** 20),
            .layered_litter_phosphorus_g_p_per_h_by_plant = &([_]f64{0} ** 20),
        }),
    );
    for (state.layered_litter_carbon_g_c_per_h) |value|
        try std.testing.expectEqual(@as(f64, 9), value);
    for (state.layered_litter_nitrogen_g_n_per_h) |value|
        try std.testing.expectEqual(@as(f64, 9), value);
    for (state.layered_litter_phosphorus_g_p_per_h) |value|
        try std.testing.expectEqual(@as(f64, 9), value);
}
