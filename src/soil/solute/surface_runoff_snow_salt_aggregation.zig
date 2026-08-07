const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const NeighborFluxMode = enum {
    blocked,
    subtract,
};

/// Runoff species follow the REDIST source update order. Each value indexes a
/// runtime-owned slice whose entries are mol per model step.
pub const RunoffSaltSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_monohydroxide,
    aluminum_dihydroxide,
    aluminum_trihydroxide,
    aluminum_tetrahydroxide,
    aluminum_sulfate,
    iron_monohydroxide,
    iron_dihydroxide,
    iron_trihydroxide,
    iron_tetrahydroxide,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_sulfate,
    phosphate,
    phosphoric_acid,
    iron_hydrogen_phosphate,
    iron_dihydrogen_phosphate,
    calcium_phosphate,
    calcium_hydrogen_phosphate,
    calcium_dihydrogen_phosphate,
    magnesium_hydrogen_phosphate,
};

/// Snow-drift species follow the REDIST source update order. Snow omits the
/// runoff hydrogen-sulfate entry, so its phosphorus indexes differ by one.
pub const SnowSaltSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_monohydroxide,
    aluminum_dihydroxide,
    aluminum_trihydroxide,
    aluminum_tetrahydroxide,
    aluminum_sulfate,
    iron_monohydroxide,
    iron_dihydroxide,
    iron_trihydroxide,
    iron_tetrahydroxide,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    phosphate,
    phosphoric_acid,
    iron_hydrogen_phosphate,
    iron_dihydrogen_phosphate,
    calcium_phosphate,
    calcium_hydrogen_phosphate,
    calcium_dihydrogen_phosphate,
    magnesium_hydrogen_phosphate,
};

pub const runoff_species_count =
    @typeInfo(RunoffSaltSpecies).@"enum".fields.len;
pub const snow_species_count =
    @typeInfo(SnowSaltSpecies).@"enum".fields.len;

pub const SaltFlux = struct {
    runoff_mol_per_step_by_species: []const f64,
    snow_mol_per_step_by_species: []const f64,
};

pub const PositiveNeighbor = struct {
    flux: ?SaltFlux,
    runoff_mode: NeighborFluxMode,
    snow_mode: NeighborFluxMode,
};

pub const Inputs = struct {
    local_flux_by_boundary_side: []const SaltFlux,
    positive_neighbor_by_boundary_side: []const PositiveNeighbor,
    /// REDIST subtracts this neighbor only during its first `NN` side.
    opposite_neighbor_first_side: ?SaltFlux,
};

pub const State = struct {
    runoff_mol_per_step_by_species: []f64,
    snow_mol_per_step_by_species: []f64,
};

/// Caller-owned runtime scratch keeps a failed candidate update from changing
/// the cell state and avoids allocation inside the transport kernel.
pub const Workspace = struct {
    runoff_mol_per_step_by_species: []f64,
    snow_mol_per_step_by_species: []f64,
};

const Operation = enum {
    add,
    subtract,
};

/// Aggregates dynamic salt transport for one surface cell and horizontal axis.
///
/// Traceability: REDIST.F lines 2139--2397 (`TQRAL`--`TQSM1P`). The runtime
/// boundary-side extent replaces the source's two fixed `NN` positions.
/// Within each side, updates retain source order: local runoff, local snow,
/// gated positive-neighbor runoff, gated positive-neighbor snow, then the
/// first-side opposite neighbor. Static equilibrium bypasses the entire source
/// block. All state is committed only after finite candidate evaluation.
pub fn aggregate(
    equilibrium_mode: SaltEquilibriumMode,
    inputs: Inputs,
    state: *State,
    workspace: Workspace,
) !void {
    if (equilibrium_mode == .static) return;

    try validateDimensions(inputs, state.*, workspace);
    try validateInputs(inputs, state.*, workspace);
    @memcpy(
        workspace.runoff_mol_per_step_by_species,
        state.runoff_mol_per_step_by_species,
    );
    @memcpy(
        workspace.snow_mol_per_step_by_species,
        state.snow_mol_per_step_by_species,
    );

    for (inputs.local_flux_by_boundary_side, 0..) |local, side| {
        try update(
            workspace.runoff_mol_per_step_by_species,
            local.runoff_mol_per_step_by_species,
            .add,
        );
        try update(
            workspace.snow_mol_per_step_by_species,
            local.snow_mol_per_step_by_species,
            .add,
        );

        const positive = inputs.positive_neighbor_by_boundary_side[side];
        if (positive.runoff_mode == .subtract) {
            try update(
                workspace.runoff_mol_per_step_by_species,
                positive.flux.?.runoff_mol_per_step_by_species,
                .subtract,
            );
        }
        if (positive.snow_mode == .subtract) {
            try update(
                workspace.snow_mol_per_step_by_species,
                positive.flux.?.snow_mol_per_step_by_species,
                .subtract,
            );
        }

        if (side == 0) {
            if (inputs.opposite_neighbor_first_side) |opposite| {
                try update(
                    workspace.runoff_mol_per_step_by_species,
                    opposite.runoff_mol_per_step_by_species,
                    .subtract,
                );
                try update(
                    workspace.snow_mol_per_step_by_species,
                    opposite.snow_mol_per_step_by_species,
                    .subtract,
                );
            }
        }
    }

    @memcpy(
        state.runoff_mol_per_step_by_species,
        workspace.runoff_mol_per_step_by_species,
    );
    @memcpy(
        state.snow_mol_per_step_by_species,
        workspace.snow_mol_per_step_by_species,
    );
}

fn validateDimensions(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (inputs.local_flux_by_boundary_side.len == 0)
        return error.InvalidSurfaceSaltAggregationDimensions;
    if (inputs.positive_neighbor_by_boundary_side.len !=
        inputs.local_flux_by_boundary_side.len)
    {
        return error.SurfaceSaltAggregationDimensionMismatch;
    }
    try validateSpeciesDimensions(.{
        .runoff_mol_per_step_by_species = state.runoff_mol_per_step_by_species,
        .snow_mol_per_step_by_species = state.snow_mol_per_step_by_species,
    });
    try validateSpeciesDimensions(.{
        .runoff_mol_per_step_by_species = workspace.runoff_mol_per_step_by_species,
        .snow_mol_per_step_by_species = workspace.snow_mol_per_step_by_species,
    });
    for (inputs.local_flux_by_boundary_side) |flux|
        try validateSpeciesDimensions(flux);
    for (inputs.positive_neighbor_by_boundary_side) |neighbor| {
        if (neighbor.flux) |flux| {
            try validateSpeciesDimensions(flux);
        } else if (neighbor.runoff_mode == .subtract or
            neighbor.snow_mode == .subtract)
        {
            return error.MissingPositiveNeighborSaltFlux;
        }
    }
    if (inputs.opposite_neighbor_first_side) |flux|
        try validateSpeciesDimensions(flux);
}

fn validateSpeciesDimensions(flux: SaltFlux) !void {
    if (flux.runoff_mol_per_step_by_species.len != runoff_species_count or
        flux.snow_mol_per_step_by_species.len != snow_species_count)
    {
        return error.SurfaceSaltAggregationDimensionMismatch;
    }
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    try validateFinite(state.runoff_mol_per_step_by_species);
    try validateFinite(state.snow_mol_per_step_by_species);

    const workspace_runoff = workspace.runoff_mol_per_step_by_species;
    const workspace_snow = workspace.snow_mol_per_step_by_species;
    if (overlap(workspace_runoff, state.runoff_mol_per_step_by_species) or
        overlap(workspace_runoff, state.snow_mol_per_step_by_species) or
        overlap(workspace_snow, state.runoff_mol_per_step_by_species) or
        overlap(workspace_snow, state.snow_mol_per_step_by_species) or
        overlap(workspace_runoff, workspace_snow))
    {
        return error.SurfaceSaltAggregationWorkspaceOverlap;
    }

    for (inputs.local_flux_by_boundary_side) |flux|
        try validateFlux(flux, workspace);
    for (inputs.positive_neighbor_by_boundary_side) |neighbor|
        if (neighbor.flux) |flux| try validateFlux(flux, workspace);
    if (inputs.opposite_neighbor_first_side) |flux|
        try validateFlux(flux, workspace);
}

fn validateFlux(flux: SaltFlux, workspace: Workspace) !void {
    try validateFinite(flux.runoff_mol_per_step_by_species);
    try validateFinite(flux.snow_mol_per_step_by_species);
    if (overlap(
        workspace.runoff_mol_per_step_by_species,
        flux.runoff_mol_per_step_by_species,
    ) or overlap(
        workspace.runoff_mol_per_step_by_species,
        flux.snow_mol_per_step_by_species,
    ) or overlap(
        workspace.snow_mol_per_step_by_species,
        flux.runoff_mol_per_step_by_species,
    ) or overlap(
        workspace.snow_mol_per_step_by_species,
        flux.snow_mol_per_step_by_species,
    )) {
        return error.SurfaceSaltAggregationWorkspaceOverlap;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceSaltAggregationInput;
}

fn update(candidate: []f64, contribution: []const f64, operation: Operation) !void {
    for (candidate, contribution) |*current, change| {
        const result = switch (operation) {
            .add => current.* + change,
            .subtract => current.* - change,
        };
        if (!std.math.isFinite(result))
            return error.NonFiniteSurfaceSaltAggregationResult;
        current.* = result;
    }
}

fn overlap(left: []const f64, right: []const f64) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(f64);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn filledFlux(
    runoff_mol_per_step_by_species: []const f64,
    snow_mol_per_step_by_species: []const f64,
) SaltFlux {
    return .{
        .runoff_mol_per_step_by_species = runoff_mol_per_step_by_species,
        .snow_mol_per_step_by_species = snow_mol_per_step_by_species,
    };
}

fn expectAllEqual(values: []const f64, expected: f64) !void {
    for (values) |value| try std.testing.expectEqual(expected, value);
}

test "source side order and independent runoff snow gates are exact" {
    const runoff_one = [_]f64{1} ** runoff_species_count;
    const runoff_two = [_]f64{2} ** runoff_species_count;
    const runoff_three = [_]f64{3} ** runoff_species_count;
    const runoff_four = [_]f64{4} ** runoff_species_count;
    const runoff_five = [_]f64{5} ** runoff_species_count;
    const snow_one = [_]f64{1} ** snow_species_count;
    const snow_two = [_]f64{2} ** snow_species_count;
    const snow_three = [_]f64{3} ** snow_species_count;
    const snow_four = [_]f64{4} ** snow_species_count;
    const snow_five = [_]f64{5} ** snow_species_count;
    const local = [_]SaltFlux{
        filledFlux(&runoff_one, &snow_one),
        filledFlux(&runoff_four, &snow_four),
    };
    const positive = [_]PositiveNeighbor{
        .{
            .flux = filledFlux(&runoff_two, &snow_two),
            .runoff_mode = .subtract,
            .snow_mode = .blocked,
        },
        .{
            .flux = filledFlux(&runoff_five, &snow_five),
            .runoff_mode = .blocked,
            .snow_mode = .subtract,
        },
    };
    var runoff_state = [_]f64{100} ** runoff_species_count;
    var snow_state = [_]f64{100} ** snow_species_count;
    var runoff_scratch: [runoff_species_count]f64 = undefined;
    var snow_scratch: [snow_species_count]f64 = undefined;
    var state = State{
        .runoff_mol_per_step_by_species = &runoff_state,
        .snow_mol_per_step_by_species = &snow_state,
    };

    try aggregate(.dynamic, .{
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_by_boundary_side = &positive,
        .opposite_neighbor_first_side = filledFlux(&runoff_three, &snow_three),
    }, &state, .{
        .runoff_mol_per_step_by_species = &runoff_scratch,
        .snow_mol_per_step_by_species = &snow_scratch,
    });

    try expectAllEqual(state.runoff_mol_per_step_by_species, 100);
    try expectAllEqual(state.snow_mol_per_step_by_species, 97);
}

test "shared horizontal face conserves every salt species exactly" {
    const runoff_seven = [_]f64{7} ** runoff_species_count;
    const snow_seven = [_]f64{7} ** snow_species_count;
    const runoff_zero = [_]f64{0} ** runoff_species_count;
    const snow_zero = [_]f64{0} ** snow_species_count;
    const first_local = [_]SaltFlux{filledFlux(&runoff_seven, &snow_seven)};
    const second_local = [_]SaltFlux{filledFlux(&runoff_zero, &snow_zero)};
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    var first_runoff = runoff_zero;
    var first_snow = snow_zero;
    var second_runoff = runoff_zero;
    var second_snow = snow_zero;
    var first_runoff_scratch: [runoff_species_count]f64 = undefined;
    var first_snow_scratch: [snow_species_count]f64 = undefined;
    var second_runoff_scratch: [runoff_species_count]f64 = undefined;
    var second_snow_scratch: [snow_species_count]f64 = undefined;
    var first = State{
        .runoff_mol_per_step_by_species = &first_runoff,
        .snow_mol_per_step_by_species = &first_snow,
    };
    var second = State{
        .runoff_mol_per_step_by_species = &second_runoff,
        .snow_mol_per_step_by_species = &second_snow,
    };

    try aggregate(.dynamic, .{
        .local_flux_by_boundary_side = &first_local,
        .positive_neighbor_by_boundary_side = &blocked,
        .opposite_neighbor_first_side = null,
    }, &first, .{
        .runoff_mol_per_step_by_species = &first_runoff_scratch,
        .snow_mol_per_step_by_species = &first_snow_scratch,
    });
    try aggregate(.dynamic, .{
        .local_flux_by_boundary_side = &second_local,
        .positive_neighbor_by_boundary_side = &blocked,
        .opposite_neighbor_first_side = filledFlux(&runoff_seven, &snow_seven),
    }, &second, .{
        .runoff_mol_per_step_by_species = &second_runoff_scratch,
        .snow_mol_per_step_by_species = &second_snow_scratch,
    });

    for (first.runoff_mol_per_step_by_species, second.runoff_mol_per_step_by_species) |a, b|
        try std.testing.expectEqual(@as(f64, 0), a + b);
    for (first.snow_mol_per_step_by_species, second.snow_mol_per_step_by_species) |a, b|
        try std.testing.expectEqual(@as(f64, 0), a + b);
}

test "species layouts and runtime side extent are explicit" {
    try std.testing.expectEqual(@as(usize, 42), runoff_species_count);
    try std.testing.expectEqual(@as(usize, 41), snow_species_count);
    try std.testing.expectEqual(
        @as(usize, 33),
        @intFromEnum(RunoffSaltSpecies.hydrogen_sulfate),
    );
    try std.testing.expectEqual(
        @as(usize, 34),
        @intFromEnum(RunoffSaltSpecies.phosphate),
    );
    try std.testing.expectEqual(
        @as(usize, 33),
        @intFromEnum(SnowSaltSpecies.phosphate),
    );

    const runoff_one = [_]f64{1} ** runoff_species_count;
    const snow_one = [_]f64{1} ** snow_species_count;
    const local = [_]SaltFlux{filledFlux(&runoff_one, &snow_one)} ** 4;
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }} ** 4;
    var runoff_state = [_]f64{0} ** runoff_species_count;
    var snow_state = [_]f64{0} ** snow_species_count;
    var runoff_scratch: [runoff_species_count]f64 = undefined;
    var snow_scratch: [snow_species_count]f64 = undefined;
    var state = State{
        .runoff_mol_per_step_by_species = &runoff_state,
        .snow_mol_per_step_by_species = &snow_state,
    };

    try aggregate(.dynamic, .{
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_by_boundary_side = &blocked,
        .opposite_neighbor_first_side = null,
    }, &state, .{
        .runoff_mol_per_step_by_species = &runoff_scratch,
        .snow_mol_per_step_by_species = &snow_scratch,
    });
    try expectAllEqual(state.runoff_mol_per_step_by_species, 4);
    try expectAllEqual(state.snow_mol_per_step_by_species, 4);
}

test "static equilibrium bypasses the dynamic aggregation block" {
    var runoff_state = [_]f64{11};
    var snow_state = [_]f64{13};
    var runoff_scratch = [_]f64{17};
    var snow_scratch = [_]f64{19};
    var state = State{
        .runoff_mol_per_step_by_species = &runoff_state,
        .snow_mol_per_step_by_species = &snow_state,
    };

    try aggregate(.static, .{
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_by_boundary_side = &.{},
        .opposite_neighbor_first_side = null,
    }, &state, .{
        .runoff_mol_per_step_by_species = &runoff_scratch,
        .snow_mol_per_step_by_species = &snow_scratch,
    });

    try std.testing.expectEqual(@as(f64, 11), runoff_state[0]);
    try std.testing.expectEqual(@as(f64, 13), snow_state[0]);
    try std.testing.expectEqual(@as(f64, 17), runoff_scratch[0]);
    try std.testing.expectEqual(@as(f64, 19), snow_scratch[0]);
}

test "late invalid input and arithmetic overflow preserve state atomically" {
    const runoff_one = [_]f64{1} ** runoff_species_count;
    var snow_invalid = [_]f64{1} ** snow_species_count;
    snow_invalid[snow_species_count - 1] = std.math.nan(f64);
    const invalid_local = [_]SaltFlux{filledFlux(&runoff_one, &snow_invalid)};
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    var runoff_state = [_]f64{3} ** runoff_species_count;
    var snow_state = [_]f64{3} ** snow_species_count;
    var runoff_scratch: [runoff_species_count]f64 = undefined;
    var snow_scratch: [snow_species_count]f64 = undefined;
    var state = State{
        .runoff_mol_per_step_by_species = &runoff_state,
        .snow_mol_per_step_by_species = &snow_state,
    };
    const workspace = Workspace{
        .runoff_mol_per_step_by_species = &runoff_scratch,
        .snow_mol_per_step_by_species = &snow_scratch,
    };

    try std.testing.expectError(
        error.NonFiniteSurfaceSaltAggregationInput,
        aggregate(.dynamic, .{
            .local_flux_by_boundary_side = &invalid_local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, workspace),
    );
    try expectAllEqual(state.runoff_mol_per_step_by_species, 3);
    try expectAllEqual(state.snow_mol_per_step_by_species, 3);

    const runoff_max = [_]f64{std.math.floatMax(f64)} ** runoff_species_count;
    const snow_zero = [_]f64{0} ** snow_species_count;
    const overflowing_local = [_]SaltFlux{filledFlux(&runoff_max, &snow_zero)};
    @memset(runoff_state[0..], std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteSurfaceSaltAggregationResult,
        aggregate(.dynamic, .{
            .local_flux_by_boundary_side = &overflowing_local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, workspace),
    );
    try expectAllEqual(state.runoff_mol_per_step_by_species, std.math.floatMax(f64));
    try expectAllEqual(state.snow_mol_per_step_by_species, 3);
}

test "shape missing neighbor and workspace alias errors precede mutation" {
    const runoff_one = [_]f64{1} ** runoff_species_count;
    const snow_one = [_]f64{1} ** snow_species_count;
    const local = [_]SaltFlux{filledFlux(&runoff_one, &snow_one)};
    const missing = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .subtract,
        .snow_mode = .blocked,
    }};
    var runoff_state = [_]f64{5} ** runoff_species_count;
    var snow_state = [_]f64{5} ** snow_species_count;
    var snow_scratch: [snow_species_count]f64 = undefined;
    var state = State{
        .runoff_mol_per_step_by_species = &runoff_state,
        .snow_mol_per_step_by_species = &snow_state,
    };

    try std.testing.expectError(
        error.MissingPositiveNeighborSaltFlux,
        aggregate(.dynamic, .{
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_by_boundary_side = &missing,
            .opposite_neighbor_first_side = null,
        }, &state, .{
            .runoff_mol_per_step_by_species = &runoff_state,
            .snow_mol_per_step_by_species = &snow_scratch,
        }),
    );
    try expectAllEqual(state.runoff_mol_per_step_by_species, 5);

    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    try std.testing.expectError(
        error.SurfaceSaltAggregationWorkspaceOverlap,
        aggregate(.dynamic, .{
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, .{
            .runoff_mol_per_step_by_species = &runoff_state,
            .snow_mol_per_step_by_species = &snow_scratch,
        }),
    );
    try expectAllEqual(state.runoff_mol_per_step_by_species, 5);
}
