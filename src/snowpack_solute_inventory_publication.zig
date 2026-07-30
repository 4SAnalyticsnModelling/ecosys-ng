const std = @import("std");

pub const SaltTransportMode = enum {
    static_equilibrium,
    dynamic_equilibrium,
};

/// Snow solutes stored on the elemental mass basis used by ecosys.
///
/// Field order preserves REDIST lines 4143--4152.
pub const PrimarySolutes = struct {
    carbon_dioxide_g_c: f64 = 0,
    methane_g_c: f64 = 0,
    oxygen_g_o: f64 = 0,
    dinitrogen_g_n: f64 = 0,
    nitrous_oxide_g_n: f64 = 0,
    ammonium_g_n: f64 = 0,
    ammonia_g_n: f64 = 0,
    nitrate_g_n: f64 = 0,
    hydrogen_phosphate_g_p: f64 = 0,
    dihydrogen_phosphate_g_p: f64 = 0,
};

pub const AcceptedIncrements = struct {
    /// Signed accepted transport increment [snow_layer].
    primary_by_layer: []const PrimarySolutes,
    /// Signed accepted increment [snow_layer][salt_species], mol.
    salt_mol_by_layer_species: []const f64,
    salt_species_count: usize,
};

pub const State = struct {
    primary_inventory_by_layer: []PrimarySolutes,
    salt_inventory_mol_by_layer_species: []f64,
};

pub const Workspace = struct {
    primary_inventory_by_layer: []PrimarySolutes,
    salt_inventory_mol_by_layer_species: []f64,
};

/// Publishes accepted snow-transport increments into chemical inventory.
///
/// Traceability: REDIST.F lines 4143--4215. The ten primary carriers are added
/// in source order for every runtime layer. Salt inventories are added in the
/// caller-declared species order only when dynamic equilibrium is enabled.
/// Static equilibrium retains salt inventory; a nonzero accepted salt increment
/// then fails because silently discarding transported mass is unsafe.
///
/// The caller owns one cell's layer-major slices and scratch. Independent cells
/// can run in parallel. Candidate construction finishes before state commit, so
/// a non-finite or negative inventory leaves the complete cell unchanged.
pub fn publishAcceptedIncrements(
    mode: SaltTransportMode,
    increments: AcceptedIncrements,
    state: *State,
    workspace: Workspace,
) !void {
    const layer_count = try validateDimensions(increments, state.*, workspace);
    try validateInputs(mode, increments, state.*, workspace);

    @memcpy(
        workspace.primary_inventory_by_layer,
        state.primary_inventory_by_layer,
    );
    @memcpy(
        workspace.salt_inventory_mol_by_layer_species,
        state.salt_inventory_mol_by_layer_species,
    );

    for (0..layer_count) |layer| {
        try addPrimaryInSourceOrder(
            &workspace.primary_inventory_by_layer[layer],
            increments.primary_by_layer[layer],
        );
    }
    if (mode == .dynamic_equilibrium) {
        for (
            workspace.salt_inventory_mol_by_layer_species,
            increments.salt_mol_by_layer_species,
        ) |*candidate, increment| {
            candidate.* = try checkedInventoryAddition(candidate.*, increment);
        }
    }

    @memcpy(
        state.primary_inventory_by_layer,
        workspace.primary_inventory_by_layer,
    );
    if (mode == .dynamic_equilibrium) {
        @memcpy(
            state.salt_inventory_mol_by_layer_species,
            workspace.salt_inventory_mol_by_layer_species,
        );
    }
}

fn addPrimaryInSourceOrder(
    candidate: *PrimarySolutes,
    increment: PrimarySolutes,
) !void {
    inline for (@typeInfo(PrimarySolutes).@"struct".fields) |field| {
        @field(candidate.*, field.name) = try checkedInventoryAddition(
            @field(candidate.*, field.name),
            @field(increment, field.name),
        );
    }
}

fn checkedInventoryAddition(inventory: f64, increment: f64) !f64 {
    const result = inventory + increment;
    if (!std.math.isFinite(result))
        return error.NonFiniteSnowpackSolutePublicationResult;
    if (result < 0) return error.NegativeSnowpackSoluteInventory;
    return result;
}

fn validateDimensions(
    increments: AcceptedIncrements,
    state: State,
    workspace: Workspace,
) !usize {
    const layer_count = state.primary_inventory_by_layer.len;
    if (layer_count == 0 or increments.salt_species_count == 0)
        return error.InvalidSnowpackSolutePublicationDimensions;
    const salt_value_count = std.math.mul(
        usize,
        layer_count,
        increments.salt_species_count,
    ) catch return error.SnowpackSolutePublicationDimensionOverflow;
    if (increments.primary_by_layer.len != layer_count or
        workspace.primary_inventory_by_layer.len != layer_count or
        increments.salt_mol_by_layer_species.len != salt_value_count or
        state.salt_inventory_mol_by_layer_species.len != salt_value_count or
        workspace.salt_inventory_mol_by_layer_species.len != salt_value_count)
    {
        return error.SnowpackSolutePublicationDimensionMismatch;
    }
    return layer_count;
}

fn validateInputs(
    mode: SaltTransportMode,
    increments: AcceptedIncrements,
    state: State,
    workspace: Workspace,
) !void {
    try validatePrimarySlice(state.primary_inventory_by_layer, .inventory);
    try validatePrimarySlice(increments.primary_by_layer, .signed_increment);
    try validateScalarSlice(
        state.salt_inventory_mol_by_layer_species,
        .inventory,
    );
    try validateScalarSlice(
        increments.salt_mol_by_layer_species,
        .signed_increment,
    );
    if (mode == .static_equilibrium) {
        for (increments.salt_mol_by_layer_species) |increment| {
            if (increment != 0)
                return error.StaticSnowSaltTransportIncrement;
        }
    }
    if (overlapPrimary(
        workspace.primary_inventory_by_layer,
        state.primary_inventory_by_layer,
    ) or overlapPrimary(
        workspace.primary_inventory_by_layer,
        increments.primary_by_layer,
    ) or overlapScalar(
        workspace.salt_inventory_mol_by_layer_species,
        state.salt_inventory_mol_by_layer_species,
    ) or overlapScalar(
        workspace.salt_inventory_mol_by_layer_species,
        increments.salt_mol_by_layer_species,
    )) {
        return error.SnowpackSolutePublicationWorkspaceOverlap;
    }
}

const ValueDomain = enum {
    inventory,
    signed_increment,
};

fn validatePrimarySlice(
    values: []const PrimarySolutes,
    domain: ValueDomain,
) !void {
    for (values) |value| {
        inline for (@typeInfo(PrimarySolutes).@"struct".fields) |field| {
            try validateScalar(@field(value, field.name), domain);
        }
    }
}

fn validateScalarSlice(values: []const f64, domain: ValueDomain) !void {
    for (values) |value| try validateScalar(value, domain);
}

fn validateScalar(value: f64, domain: ValueDomain) !void {
    if (!std.math.isFinite(value))
        return error.NonFiniteSnowpackSolutePublicationInput;
    if (domain == .inventory and value < 0)
        return error.NegativeSnowpackSoluteInventory;
}

fn overlapPrimary(
    left: []const PrimarySolutes,
    right: []const PrimarySolutes,
) bool {
    return rangesOverlap(
        @intFromPtr(left.ptr),
        extentBytes(left.len, @sizeOf(PrimarySolutes)) catch return true,
        @intFromPtr(right.ptr),
        extentBytes(right.len, @sizeOf(PrimarySolutes)) catch return true,
    );
}

fn overlapScalar(left: []const f64, right: []const f64) bool {
    return rangesOverlap(
        @intFromPtr(left.ptr),
        extentBytes(left.len, @sizeOf(f64)) catch return true,
        @intFromPtr(right.ptr),
        extentBytes(right.len, @sizeOf(f64)) catch return true,
    );
}

fn extentBytes(count: usize, item_size: usize) !usize {
    return std.math.mul(usize, count, item_size) catch
        error.SnowpackSolutePublicationDimensionOverflow;
}

fn rangesOverlap(
    left_start: usize,
    left_bytes: usize,
    right_start: usize,
    right_bytes: usize,
) bool {
    if (left_bytes == 0 or right_bytes == 0) return false;
    const left_end = std.math.add(usize, left_start, left_bytes) catch
        return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch
        return true;
    return left_start < right_end and right_start < left_end;
}

fn filledPrimary(value: f64) PrimarySolutes {
    var result: PrimarySolutes = undefined;
    inline for (@typeInfo(PrimarySolutes).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn sumPrimary(value: PrimarySolutes) f64 {
    var result: f64 = 0;
    inline for (@typeInfo(PrimarySolutes).@"struct".fields) |field|
        result += @field(value, field.name);
    return result;
}

test "REDIST primary solutes accumulate in source field order" {
    const allocator = std.testing.allocator;
    const primary_state = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_state);
    const primary_increment = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_increment);
    const primary_workspace = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_workspace);
    const salt_state = try allocator.alloc(f64, 3);
    defer allocator.free(salt_state);
    const salt_increment = try allocator.alloc(f64, 3);
    defer allocator.free(salt_increment);
    const salt_workspace = try allocator.alloc(f64, 3);
    defer allocator.free(salt_workspace);

    primary_state[0] = filledPrimary(10);
    primary_increment[0] = .{
        .carbon_dioxide_g_c = 1,
        .methane_g_c = 2,
        .oxygen_g_o = 3,
        .dinitrogen_g_n = 4,
        .nitrous_oxide_g_n = 5,
        .ammonium_g_n = 6,
        .ammonia_g_n = 7,
        .nitrate_g_n = 8,
        .hydrogen_phosphate_g_p = 9,
        .dihydrogen_phosphate_g_p = 10,
    };
    @memset(salt_state, 2);
    @memset(salt_increment, 0);
    var state = State{
        .primary_inventory_by_layer = primary_state,
        .salt_inventory_mol_by_layer_species = salt_state,
    };
    try publishAcceptedIncrements(
        .static_equilibrium,
        .{
            .primary_by_layer = primary_increment,
            .salt_mol_by_layer_species = salt_increment,
            .salt_species_count = 3,
        },
        &state,
        .{
            .primary_inventory_by_layer = primary_workspace,
            .salt_inventory_mol_by_layer_species = salt_workspace,
        },
    );

    try std.testing.expectEqual(
        @as(f64, 11),
        state.primary_inventory_by_layer[0].carbon_dioxide_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 15),
        state.primary_inventory_by_layer[0].nitrous_oxide_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, 20),
        state.primary_inventory_by_layer[0].dihydrogen_phosphate_g_p,
    );
    try std.testing.expectEqual(@as(f64, 155), sumPrimary(primary_state[0]));
}

test "dynamic salt publication supports runtime layers and species" {
    const allocator = std.testing.allocator;
    const layer_count: usize = 3;
    const salt_species_count: usize = 17;
    const salt_value_count = layer_count * salt_species_count;
    const primary_state = try allocator.alloc(PrimarySolutes, layer_count);
    defer allocator.free(primary_state);
    const primary_increment = try allocator.alloc(PrimarySolutes, layer_count);
    defer allocator.free(primary_increment);
    const primary_workspace = try allocator.alloc(PrimarySolutes, layer_count);
    defer allocator.free(primary_workspace);
    const salt_state = try allocator.alloc(f64, salt_value_count);
    defer allocator.free(salt_state);
    const salt_increment = try allocator.alloc(f64, salt_value_count);
    defer allocator.free(salt_increment);
    const salt_workspace = try allocator.alloc(f64, salt_value_count);
    defer allocator.free(salt_workspace);
    @memset(primary_state, filledPrimary(5));
    @memset(primary_increment, filledPrimary(-1));
    @memset(salt_state, 4);
    @memset(salt_increment, 0.25);
    var state = State{
        .primary_inventory_by_layer = primary_state,
        .salt_inventory_mol_by_layer_species = salt_state,
    };

    try publishAcceptedIncrements(
        .dynamic_equilibrium,
        .{
            .primary_by_layer = primary_increment,
            .salt_mol_by_layer_species = salt_increment,
            .salt_species_count = salt_species_count,
        },
        &state,
        .{
            .primary_inventory_by_layer = primary_workspace,
            .salt_inventory_mol_by_layer_species = salt_workspace,
        },
    );

    for (primary_state) |primary|
        try std.testing.expectEqual(@as(f64, 40), sumPrimary(primary));
    var salt_total_mol: f64 = 0;
    for (salt_state) |inventory_mol| {
        try std.testing.expectEqual(@as(f64, 4.25), inventory_mol);
        salt_total_mol += inventory_mol;
    }
    try std.testing.expectEqual(
        @as(f64, @floatFromInt(salt_value_count)) * 4.25,
        salt_total_mol,
    );
}

test "static equilibrium rejects rather than discards salt transport" {
    const allocator = std.testing.allocator;
    const primary_state = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_state);
    const primary_increment = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_increment);
    const primary_workspace = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_workspace);
    const salt_state = try allocator.alloc(f64, 2);
    defer allocator.free(salt_state);
    const salt_increment = try allocator.alloc(f64, 2);
    defer allocator.free(salt_increment);
    const salt_workspace = try allocator.alloc(f64, 2);
    defer allocator.free(salt_workspace);
    @memset(primary_state, filledPrimary(1));
    @memset(primary_increment, filledPrimary(1));
    @memset(salt_state, 7);
    @memset(salt_increment, 0);
    salt_increment[1] = 0.5;
    var state = State{
        .primary_inventory_by_layer = primary_state,
        .salt_inventory_mol_by_layer_species = salt_state,
    };

    try std.testing.expectError(
        error.StaticSnowSaltTransportIncrement,
        publishAcceptedIncrements(
            .static_equilibrium,
            .{
                .primary_by_layer = primary_increment,
                .salt_mol_by_layer_species = salt_increment,
                .salt_species_count = 2,
            },
            &state,
            .{
                .primary_inventory_by_layer = primary_workspace,
                .salt_inventory_mol_by_layer_species = salt_workspace,
            },
        ),
    );
    try std.testing.expectEqual(@as(f64, 10), sumPrimary(primary_state[0]));
    try std.testing.expectEqual(@as(f64, 7), salt_state[1]);
}

test "negative dynamic salt result leaves every inventory unchanged" {
    const allocator = std.testing.allocator;
    const primary_state = try allocator.alloc(PrimarySolutes, 2);
    defer allocator.free(primary_state);
    const primary_increment = try allocator.alloc(PrimarySolutes, 2);
    defer allocator.free(primary_increment);
    const primary_workspace = try allocator.alloc(PrimarySolutes, 2);
    defer allocator.free(primary_workspace);
    const salt_state = try allocator.alloc(f64, 6);
    defer allocator.free(salt_state);
    const salt_increment = try allocator.alloc(f64, 6);
    defer allocator.free(salt_increment);
    const salt_workspace = try allocator.alloc(f64, 6);
    defer allocator.free(salt_workspace);
    @memset(primary_state, filledPrimary(3));
    @memset(primary_increment, filledPrimary(2));
    @memset(salt_state, 1);
    @memset(salt_increment, 0);
    salt_increment[5] = -2;
    var state = State{
        .primary_inventory_by_layer = primary_state,
        .salt_inventory_mol_by_layer_species = salt_state,
    };

    try std.testing.expectError(
        error.NegativeSnowpackSoluteInventory,
        publishAcceptedIncrements(
            .dynamic_equilibrium,
            .{
                .primary_by_layer = primary_increment,
                .salt_mol_by_layer_species = salt_increment,
                .salt_species_count = 3,
            },
            &state,
            .{
                .primary_inventory_by_layer = primary_workspace,
                .salt_inventory_mol_by_layer_species = salt_workspace,
            },
        ),
    );
    for (primary_state) |primary|
        try std.testing.expectEqual(@as(f64, 30), sumPrimary(primary));
    for (salt_state) |inventory_mol|
        try std.testing.expectEqual(@as(f64, 1), inventory_mol);
}

test "dimension non-finite and workspace overlap failures are explicit" {
    const allocator = std.testing.allocator;
    const primary_state = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_state);
    const primary_increment = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_increment);
    const primary_workspace = try allocator.alloc(PrimarySolutes, 1);
    defer allocator.free(primary_workspace);
    const salt_state = try allocator.alloc(f64, 2);
    defer allocator.free(salt_state);
    const salt_increment = try allocator.alloc(f64, 2);
    defer allocator.free(salt_increment);
    const salt_workspace = try allocator.alloc(f64, 2);
    defer allocator.free(salt_workspace);
    @memset(primary_state, filledPrimary(1));
    @memset(primary_increment, filledPrimary(0));
    @memset(salt_state, 1);
    @memset(salt_increment, 0);
    var state = State{
        .primary_inventory_by_layer = primary_state,
        .salt_inventory_mol_by_layer_species = salt_state,
    };
    const valid_increments = AcceptedIncrements{
        .primary_by_layer = primary_increment,
        .salt_mol_by_layer_species = salt_increment,
        .salt_species_count = 2,
    };

    try std.testing.expectError(
        error.SnowpackSolutePublicationDimensionMismatch,
        publishAcceptedIncrements(
            .dynamic_equilibrium,
            .{
                .primary_by_layer = primary_increment,
                .salt_mol_by_layer_species = salt_increment[0..1],
                .salt_species_count = 2,
            },
            &state,
            .{
                .primary_inventory_by_layer = primary_workspace,
                .salt_inventory_mol_by_layer_species = salt_workspace,
            },
        ),
    );
    salt_increment[1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSnowpackSolutePublicationInput,
        publishAcceptedIncrements(
            .dynamic_equilibrium,
            valid_increments,
            &state,
            .{
                .primary_inventory_by_layer = primary_workspace,
                .salt_inventory_mol_by_layer_species = salt_workspace,
            },
        ),
    );
    salt_increment[1] = 0;
    try std.testing.expectError(
        error.SnowpackSolutePublicationWorkspaceOverlap,
        publishAcceptedIncrements(
            .dynamic_equilibrium,
            valid_increments,
            &state,
            .{
                .primary_inventory_by_layer = primary_state,
                .salt_inventory_mol_by_layer_species = salt_workspace,
            },
        ),
    );
}
