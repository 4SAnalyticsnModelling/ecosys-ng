const std = @import("std");

pub const TransportAxis = enum {
    column,
    row,
    vertical,
};

/// Site-file NCNG semantics for horizontal exchange.
pub const HorizontalExchange = enum {
    enabled,
    standalone,
};

/// Dissolved organic pools for one runtime organic-matter class.
pub const OrganicClassFlux = struct {
    micropore_dissolved_organic_carbon_g_per_step: f64 = 0,
    micropore_dissolved_organic_nitrogen_g_per_step: f64 = 0,
    micropore_dissolved_organic_phosphorus_g_per_step: f64 = 0,
    micropore_acetate_g_per_step: f64 = 0,
    macropore_dissolved_organic_carbon_g_per_step: f64 = 0,
    macropore_dissolved_organic_nitrogen_g_per_step: f64 = 0,
    macropore_dissolved_organic_phosphorus_g_per_step: f64 = 0,
    macropore_acetate_g_per_step: f64 = 0,
};

pub const DissolvedGasFlux = struct {
    carbon_dioxide_g_per_step: f64 = 0,
    methane_g_per_step: f64 = 0,
    oxygen_g_per_step: f64 = 0,
    dinitrogen_g_per_step: f64 = 0,
    nitrous_oxide_g_per_step: f64 = 0,
    hydrogen_g_per_step: f64 = 0,
};

pub const NutrientFlux = struct {
    ammonium_g_per_step: f64 = 0,
    ammonia_g_per_step: f64 = 0,
    nitrate_g_per_step: f64 = 0,
    nitrite_g_per_step: f64 = 0,
    hydrogen_phosphate_g_per_step: f64 = 0,
    dihydrogen_phosphate_g_per_step: f64 = 0,
};

pub const PoreDomainSoluteFlux = struct {
    dissolved_gases: DissolvedGasFlux = .{},
    nonband_nutrients: NutrientFlux = .{},
    band_nutrients: NutrientFlux = .{},
};

pub const SoluteFlux = struct {
    /// [organic_class], with every class extent supplied at runtime.
    organic_by_class: []const OrganicClassFlux,
    micropore: PoreDomainSoluteFlux,
    macropore: PoreDomainSoluteFlux,
};

pub const Inputs = struct {
    transport_axis: TransportAxis,
    horizontal_exchange: HorizontalExchange,
    current_layer_thickness_m: f64,
    layer_activity_threshold_m: f64,
    organic_class_count: usize,
    current_cell_flux: SoluteFlux,
    positive_neighbor_flux: SoluteFlux,
};

pub const State = struct {
    net_organic_by_class: []OrganicClassFlux,
    net_micropore: PoreDomainSoluteFlux,
    net_macropore: PoreDomainSoluteFlux,
};

pub const Workspace = struct {
    net_organic_by_class: []OrganicClassFlux,
};

/// Aggregates net dissolved-solute transfer between adjacent cells.
///
/// Traceability: REDIST.F lines 3414--3506 under enclosing gates 3350 and
/// 3363. Runtime organic classes replace fixed `K=0,4`. For each class the
/// source order is micropore DOC, DON, DOP, acetate, then the corresponding
/// macropore pools. Scalar order then remains micropore gases, non-band
/// nutrients, band nutrients, followed by those three macropore groups.
/// Every inventory is g per model step. Net transfer is current-cell face flux
/// minus the already-selected positive-neighbor `N6` face flux. Candidate
/// state commits atomically after source-ordered finite evaluation.
pub fn aggregate(inputs: Inputs, state: *State, workspace: Workspace) !void {
    if (inputs.transport_axis != .vertical and
        inputs.horizontal_exchange == .standalone)
    {
        return;
    }
    try validateInputs(inputs, state.*, workspace);
    if (inputs.current_layer_thickness_m <= inputs.layer_activity_threshold_m)
        return;

    @memcpy(workspace.net_organic_by_class, state.net_organic_by_class);
    for (
        workspace.net_organic_by_class,
        inputs.current_cell_flux.organic_by_class,
        inputs.positive_neighbor_flux.organic_by_class,
    ) |*candidate, current, positive| {
        try addDifference(candidate, current, positive);
    }

    var micropore_candidate = state.net_micropore;
    try addPoreDomainDifference(
        &micropore_candidate,
        inputs.current_cell_flux.micropore,
        inputs.positive_neighbor_flux.micropore,
    );
    var macropore_candidate = state.net_macropore;
    try addPoreDomainDifference(
        &macropore_candidate,
        inputs.current_cell_flux.macropore,
        inputs.positive_neighbor_flux.macropore,
    );

    @memcpy(state.net_organic_by_class, workspace.net_organic_by_class);
    state.net_micropore = micropore_candidate;
    state.net_macropore = macropore_candidate;
}

fn addPoreDomainDifference(
    candidate: *PoreDomainSoluteFlux,
    current: PoreDomainSoluteFlux,
    positive: PoreDomainSoluteFlux,
) !void {
    try addDifference(
        &candidate.dissolved_gases,
        current.dissolved_gases,
        positive.dissolved_gases,
    );
    try addDifference(
        &candidate.nonband_nutrients,
        current.nonband_nutrients,
        positive.nonband_nutrients,
    );
    try addDifference(
        &candidate.band_nutrients,
        current.band_nutrients,
        positive.band_nutrients,
    );
}

fn addDifference(candidate: anytype, current: anytype, positive: anytype) !void {
    inline for (@typeInfo(@TypeOf(candidate.*)).@"struct".fields) |field| {
        const with_current = @field(candidate.*, field.name) +
            @field(current, field.name);
        if (!std.math.isFinite(with_current))
            return error.NonFiniteAdjacentSoluteResult;
        const result = with_current - @field(positive, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteAdjacentSoluteResult;
        @field(candidate.*, field.name) = result;
    }
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (inputs.organic_class_count == 0)
        return error.InvalidAdjacentSoluteDimensions;
    if (inputs.current_cell_flux.organic_by_class.len !=
        inputs.organic_class_count or
        inputs.positive_neighbor_flux.organic_by_class.len !=
            inputs.organic_class_count or
        state.net_organic_by_class.len != inputs.organic_class_count or
        workspace.net_organic_by_class.len != inputs.organic_class_count)
    {
        return error.AdjacentSoluteDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.current_layer_thickness_m) or
        !std.math.isFinite(inputs.layer_activity_threshold_m))
    {
        return error.NonFiniteAdjacentSoluteInput;
    }
    if (inputs.current_layer_thickness_m < 0 or
        inputs.layer_activity_threshold_m < 0)
    {
        return error.InvalidAdjacentSoluteThickness;
    }
    try validateOrganicFinite(inputs.current_cell_flux.organic_by_class);
    try validateOrganicFinite(inputs.positive_neighbor_flux.organic_by_class);
    try validateOrganicFinite(state.net_organic_by_class);
    try validatePoreDomain(inputs.current_cell_flux.micropore);
    try validatePoreDomain(inputs.current_cell_flux.macropore);
    try validatePoreDomain(inputs.positive_neighbor_flux.micropore);
    try validatePoreDomain(inputs.positive_neighbor_flux.macropore);
    try validatePoreDomain(state.net_micropore);
    try validatePoreDomain(state.net_macropore);

    const scratch = workspace.net_organic_by_class;
    if (overlap(scratch, state.net_organic_by_class) or
        overlap(scratch, inputs.current_cell_flux.organic_by_class) or
        overlap(scratch, inputs.positive_neighbor_flux.organic_by_class))
    {
        return error.AdjacentSoluteWorkspaceOverlap;
    }
}

fn validateOrganicFinite(values: []const OrganicClassFlux) !void {
    for (values) |value|
        try validateStruct(value);
}

fn validatePoreDomain(value: PoreDomainSoluteFlux) !void {
    try validateStruct(value.dissolved_gases);
    try validateStruct(value.nonband_nutrients);
    try validateStruct(value.band_nutrients);
}

fn validateStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteAdjacentSoluteInput;
}

fn overlap(
    left: []const OrganicClassFlux,
    right: []const OrganicClassFlux,
) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(OrganicClassFlux);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .@"struct")
            @field(result, field.name) = filled(field.type, value)
        else
            @field(result, field.name) = value;
    }
    return result;
}

fn expectStruct(actual: anytype, expected: f64) !void {
    inline for (@typeInfo(@TypeOf(actual)).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .@"struct")
            try expectStruct(@field(actual, field.name), expected)
        else
            try std.testing.expectEqual(expected, @field(actual, field.name));
    }
}

fn expectOrganic(values: []const OrganicClassFlux, expected: f64) !void {
    for (values) |value| try expectStruct(value, expected);
}

fn baseInputs(
    current_organic: []const OrganicClassFlux,
    positive_organic: []const OrganicClassFlux,
    current_scalar: f64,
    positive_scalar: f64,
) Inputs {
    return .{
        .transport_axis = .column,
        .horizontal_exchange = .enabled,
        .current_layer_thickness_m = 0.2,
        .layer_activity_threshold_m = 0.1,
        .organic_class_count = current_organic.len,
        .current_cell_flux = .{
            .organic_by_class = current_organic,
            .micropore = filled(PoreDomainSoluteFlux, current_scalar),
            .macropore = filled(PoreDomainSoluteFlux, current_scalar),
        },
        .positive_neighbor_flux = .{
            .organic_by_class = positive_organic,
            .micropore = filled(PoreDomainSoluteFlux, positive_scalar),
            .macropore = filled(PoreDomainSoluteFlux, positive_scalar),
        },
    };
}

test "runtime organic classes and all scalar solutes follow face difference" {
    const current = [_]OrganicClassFlux{
        filled(OrganicClassFlux, 5),
    } ** 3;
    const positive = [_]OrganicClassFlux{
        filled(OrganicClassFlux, 2),
    } ** 3;
    var organic_state = [_]OrganicClassFlux{
        filled(OrganicClassFlux, 100),
    } ** 3;
    var organic_scratch: [3]OrganicClassFlux = undefined;
    var state = State{
        .net_organic_by_class = &organic_state,
        .net_micropore = filled(PoreDomainSoluteFlux, 100),
        .net_macropore = filled(PoreDomainSoluteFlux, 100),
    };
    try aggregate(
        baseInputs(&current, &positive, 5, 2),
        &state,
        .{ .net_organic_by_class = &organic_scratch },
    );
    try expectOrganic(&organic_state, 103);
    try expectStruct(state.net_micropore, 103);
    try expectStruct(state.net_macropore, 103);
}

test "reversed adjacent faces are exactly antisymmetric for every pool" {
    const first_organic = [_]OrganicClassFlux{
        filled(OrganicClassFlux, 7),
    } ** 2;
    const second_organic = [_]OrganicClassFlux{
        filled(OrganicClassFlux, 3),
    } ** 2;
    var first_state_organic = [_]OrganicClassFlux{.{}} ** 2;
    var second_state_organic = [_]OrganicClassFlux{.{}} ** 2;
    var first_scratch: [2]OrganicClassFlux = undefined;
    var second_scratch: [2]OrganicClassFlux = undefined;
    var first_state = State{
        .net_organic_by_class = &first_state_organic,
        .net_micropore = .{},
        .net_macropore = .{},
    };
    var second_state = State{
        .net_organic_by_class = &second_state_organic,
        .net_micropore = .{},
        .net_macropore = .{},
    };
    try aggregate(
        baseInputs(&first_organic, &second_organic, 7, 3),
        &first_state,
        .{ .net_organic_by_class = &first_scratch },
    );
    try aggregate(
        baseInputs(&second_organic, &first_organic, 3, 7),
        &second_state,
        .{ .net_organic_by_class = &second_scratch },
    );

    for (first_state_organic, second_state_organic) |first, second| {
        inline for (@typeInfo(OrganicClassFlux).@"struct".fields) |field|
            try std.testing.expectEqual(
                @as(f64, 0),
                @field(first, field.name) + @field(second, field.name),
            );
    }
    inline for (@typeInfo(PoreDomainSoluteFlux).@"struct".fields) |group| {
        inline for (@typeInfo(group.type).@"struct".fields) |field| {
            try std.testing.expectEqual(
                @as(f64, 0),
                @field(@field(first_state.net_micropore, group.name), field.name) +
                    @field(@field(second_state.net_micropore, group.name), field.name),
            );
            try std.testing.expectEqual(
                @as(f64, 0),
                @field(@field(first_state.net_macropore, group.name), field.name) +
                    @field(@field(second_state.net_macropore, group.name), field.name),
            );
        }
    }
}

test "equal face values cause exact zero increment in all solute pools" {
    const shared = [_]OrganicClassFlux{
        filled(OrganicClassFlux, 9),
    } ** 4;
    var organic_state = [_]OrganicClassFlux{.{}} ** 4;
    var scratch: [4]OrganicClassFlux = undefined;
    var state = State{
        .net_organic_by_class = &organic_state,
        .net_micropore = .{},
        .net_macropore = .{},
    };
    try aggregate(
        baseInputs(&shared, &shared, 9, 9),
        &state,
        .{ .net_organic_by_class = &scratch },
    );
    try expectOrganic(&organic_state, 0);
    try expectStruct(state.net_micropore, 0);
    try expectStruct(state.net_macropore, 0);
}

test "standalone horizontal and inactive layer gates do not mutate" {
    var state = State{
        .net_organic_by_class = &.{},
        .net_micropore = filled(PoreDomainSoluteFlux, 9),
        .net_macropore = filled(PoreDomainSoluteFlux, 9),
    };
    var inputs = Inputs{
        .transport_axis = .row,
        .horizontal_exchange = .standalone,
        .current_layer_thickness_m = std.math.nan(f64),
        .layer_activity_threshold_m = std.math.nan(f64),
        .organic_class_count = 0,
        .current_cell_flux = .{
            .organic_by_class = &.{},
            .micropore = .{},
            .macropore = .{},
        },
        .positive_neighbor_flux = .{
            .organic_by_class = &.{},
            .micropore = .{},
            .macropore = .{},
        },
    };
    try aggregate(inputs, &state, .{ .net_organic_by_class = &.{} });
    try expectStruct(state.net_micropore, 9);

    const zero = [_]OrganicClassFlux{.{}};
    var state_organic = [_]OrganicClassFlux{filled(OrganicClassFlux, 9)};
    var scratch: [1]OrganicClassFlux = undefined;
    state.net_organic_by_class = &state_organic;
    inputs = baseInputs(&zero, &zero, 100, 0);
    inputs.current_layer_thickness_m = 0.1;
    try aggregate(
        inputs,
        &state,
        .{ .net_organic_by_class = &scratch },
    );
    try expectOrganic(&state_organic, 9);
    try expectStruct(state.net_micropore, 9);
}

test "dimension finite alias and overflow failures preserve state atomically" {
    const current = [_]OrganicClassFlux{
        filled(OrganicClassFlux, std.math.floatMax(f64)),
    } ** 2;
    const positive = [_]OrganicClassFlux{
        filled(OrganicClassFlux, -std.math.floatMax(f64)),
    } ** 2;
    var organic_state = [_]OrganicClassFlux{
        filled(OrganicClassFlux, std.math.floatMax(f64)),
    } ** 2;
    var scratch: [2]OrganicClassFlux = undefined;
    var state = State{
        .net_organic_by_class = &organic_state,
        .net_micropore = filled(PoreDomainSoluteFlux, 5),
        .net_macropore = filled(PoreDomainSoluteFlux, 5),
    };
    var inputs = baseInputs(
        &current,
        &positive,
        std.math.floatMax(f64),
        -std.math.floatMax(f64),
    );

    inputs.organic_class_count = 3;
    try std.testing.expectError(
        error.AdjacentSoluteDimensionMismatch,
        aggregate(
            inputs,
            &state,
            .{ .net_organic_by_class = &scratch },
        ),
    );
    try expectOrganic(&organic_state, std.math.floatMax(f64));

    inputs.organic_class_count = 2;
    try std.testing.expectError(
        error.AdjacentSoluteWorkspaceOverlap,
        aggregate(
            inputs,
            &state,
            .{ .net_organic_by_class = &organic_state },
        ),
    );
    try expectOrganic(&organic_state, std.math.floatMax(f64));

    try std.testing.expectError(
        error.NonFiniteAdjacentSoluteResult,
        aggregate(
            inputs,
            &state,
            .{ .net_organic_by_class = &scratch },
        ),
    );
    try expectOrganic(&organic_state, std.math.floatMax(f64));
    try expectStruct(state.net_micropore, 5);
    try expectStruct(state.net_macropore, 5);
}
