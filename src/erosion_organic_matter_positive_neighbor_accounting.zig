const std = @import("std");

pub const DisturbanceMode = enum {
    no_profile_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter,
    freeze_thaw_erosion_and_organic_matter,
};

pub const TransportAxis = enum {
    east_west,
    north_south,
    vertical,
};

pub const NeighborConnection = enum {
    blocked,
    connected,
};

pub const ElementalOrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
};

pub const AdsorbedOrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
    acetate_carbon_g_c_per_step: f64 = 0,
};

pub const SoilOrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    colonized_carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
};

pub const OrganicDimensions = struct {
    microbial_substrate_class_count: usize,
    microbial_group_count: usize,
    microbial_component_count: usize,
    stable_organic_class_count: usize,
    residue_component_count: usize,
    soil_organic_fraction_count: usize,
};

pub const OrganicMatterFlux = struct {
    /// [microbial class][microbial group][component], component fastest.
    microbial_by_class_group_component: []const ElementalOrganicFlux,
    /// [stable class][residue component], component fastest.
    residue_by_class_component: []const ElementalOrganicFlux,
    /// [stable class].
    adsorbed_by_class: []const AdsorbedOrganicFlux,
    /// [stable class][soil organic fraction], fraction fastest.
    soil_by_class_fraction: []const SoilOrganicFlux,
};

pub const NeighborSide = struct {
    local_total_sediment_Mg_per_step: f64,
    positive_neighbor_total_sediment_Mg_per_step: f64,
    positive_neighbor_organic_matter: OrganicMatterFlux,
    connection: NeighborConnection,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_Mg_per_step: f64,
    dimensions: OrganicDimensions,
    neighbor_by_boundary_side: []const NeighborSide,
};

pub const State = struct {
    microbial_by_class_group_component: []ElementalOrganicFlux,
    residue_by_class_component: []ElementalOrganicFlux,
    adsorbed_by_class: []AdsorbedOrganicFlux,
    soil_by_class_fraction: []SoilOrganicFlux,
};

pub const Workspace = State;

const Extents = struct {
    microbial: usize,
    residue: usize,
    adsorbed: usize,
    soil: usize,
};

const ByteRegion = struct {
    start: usize,
    byte_len: usize,
};

/// Subtracts connected positive-neighbor organic-matter erosion.
///
/// Traceability: REDIST.F lines 3160--3186, within gates 2888, 2894--2895,
/// and 3051. All pool dimensions are runtime values. Flattening retains
/// Fortran traversal with component/fraction fastest; stable pools are updated
/// class-by-class in residue, adsorbed, then soil order. Elemental inventories
/// remain g C, g N, or g P per model step. Candidate state commits atomically.
pub fn account(inputs: Inputs, state: *State, workspace: Workspace) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical)
    {
        return;
    }
    const extents = try validateDimensions(inputs, state.*, workspace);
    try validateInputs(inputs, state.*, workspace);
    copyState(workspace, state.*);

    for (inputs.neighbor_by_boundary_side) |side| {
        if (@abs(side.local_total_sediment_Mg_per_step) <=
            inputs.sediment_activity_threshold_Mg_per_step and
            @abs(side.positive_neighbor_total_sediment_Mg_per_step) <=
                inputs.sediment_activity_threshold_Mg_per_step)
        {
            continue;
        }
        if (side.connection == .blocked) continue;
        try subtractOrganicMatter(
            inputs.dimensions,
            extents,
            workspace,
            side.positive_neighbor_organic_matter,
        );
    }
    copyState(state.*, workspace);
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateDimensions(
    inputs: Inputs,
    state: State,
    workspace: Workspace,
) !Extents {
    inline for (@typeInfo(OrganicDimensions).@"struct".fields) |field|
        if (@field(inputs.dimensions, field.name) == 0)
            return error.InvalidOrganicNeighborErosionDimensions;
    if (inputs.neighbor_by_boundary_side.len == 0)
        return error.InvalidOrganicNeighborErosionDimensions;

    const dimensions = inputs.dimensions;
    const extents = Extents{
        .microbial = try checkedProduct3(
            dimensions.microbial_substrate_class_count,
            dimensions.microbial_group_count,
            dimensions.microbial_component_count,
        ),
        .residue = try checkedProduct(
            dimensions.stable_organic_class_count,
            dimensions.residue_component_count,
        ),
        .adsorbed = dimensions.stable_organic_class_count,
        .soil = try checkedProduct(
            dimensions.stable_organic_class_count,
            dimensions.soil_organic_fraction_count,
        ),
    };
    try validateStateDimensions(state, extents);
    try validateStateDimensions(workspace, extents);
    for (inputs.neighbor_by_boundary_side) |side|
        try validateFluxDimensions(
            side.positive_neighbor_organic_matter,
            extents,
        );
    return extents;
}

fn checkedProduct(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        error.InvalidOrganicNeighborErosionDimensions;
}

fn checkedProduct3(first: usize, second: usize, third: usize) !usize {
    return checkedProduct(try checkedProduct(first, second), third);
}

fn validateStateDimensions(state: State, extents: Extents) !void {
    if (state.microbial_by_class_group_component.len != extents.microbial or
        state.residue_by_class_component.len != extents.residue or
        state.adsorbed_by_class.len != extents.adsorbed or
        state.soil_by_class_fraction.len != extents.soil)
    {
        return error.OrganicNeighborErosionDimensionMismatch;
    }
}

fn validateFluxDimensions(flux: OrganicMatterFlux, extents: Extents) !void {
    if (flux.microbial_by_class_group_component.len != extents.microbial or
        flux.residue_by_class_component.len != extents.residue or
        flux.adsorbed_by_class.len != extents.adsorbed or
        flux.soil_by_class_fraction.len != extents.soil)
    {
        return error.OrganicNeighborErosionDimensionMismatch;
    }
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (!std.math.isFinite(inputs.sediment_activity_threshold_Mg_per_step))
        return error.NonFiniteOrganicNeighborErosionInput;
    if (inputs.sediment_activity_threshold_Mg_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateStateFinite(state);
    try validateWorkspaceOwnership(workspace, state, inputs);
    for (inputs.neighbor_by_boundary_side) |side| {
        if (!std.math.isFinite(side.local_total_sediment_Mg_per_step) or
            !std.math.isFinite(
                side.positive_neighbor_total_sediment_Mg_per_step,
            ))
        {
            return error.NonFiniteOrganicNeighborErosionInput;
        }
        try validateFluxFinite(side.positive_neighbor_organic_matter);
    }
}

fn validateStateFinite(state: State) !void {
    try validateFinite(state.microbial_by_class_group_component);
    try validateFinite(state.residue_by_class_component);
    try validateFinite(state.adsorbed_by_class);
    try validateFinite(state.soil_by_class_fraction);
}

fn validateFluxFinite(flux: OrganicMatterFlux) !void {
    try validateFinite(flux.microbial_by_class_group_component);
    try validateFinite(flux.residue_by_class_component);
    try validateFinite(flux.adsorbed_by_class);
    try validateFinite(flux.soil_by_class_fraction);
}

fn validateFinite(values: anytype) !void {
    for (values) |value| {
        inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
            if (!std.math.isFinite(@field(value, field.name)))
                return error.NonFiniteOrganicNeighborErosionInput;
    }
}

fn validateWorkspaceOwnership(
    workspace: Workspace,
    state: State,
    inputs: Inputs,
) !void {
    const workspace_regions = stateRegions(workspace);
    const state_regions = stateRegions(state);
    try requireDisjoint(workspace_regions[0..], workspace_regions[0..], true);
    try requireDisjoint(workspace_regions[0..], state_regions[0..], false);
    for (inputs.neighbor_by_boundary_side) |side| {
        const flux_regions = fluxRegions(side.positive_neighbor_organic_matter);
        try requireDisjoint(workspace_regions[0..], flux_regions[0..], false);
    }
}

fn stateRegions(state: State) [4]ByteRegion {
    return .{
        byteRegion(state.microbial_by_class_group_component),
        byteRegion(state.residue_by_class_component),
        byteRegion(state.adsorbed_by_class),
        byteRegion(state.soil_by_class_fraction),
    };
}

fn fluxRegions(flux: OrganicMatterFlux) [4]ByteRegion {
    return .{
        byteRegion(flux.microbial_by_class_group_component),
        byteRegion(flux.residue_by_class_component),
        byteRegion(flux.adsorbed_by_class),
        byteRegion(flux.soil_by_class_fraction),
    };
}

fn byteRegion(values: anytype) ByteRegion {
    const Pointer = @typeInfo(@TypeOf(values)).pointer;
    return .{
        .start = @intFromPtr(values.ptr),
        .byte_len = std.math.mul(usize, values.len, @sizeOf(Pointer.child)) catch std.math.maxInt(usize),
    };
}

fn requireDisjoint(
    left: []const ByteRegion,
    right: []const ByteRegion,
    same_collection: bool,
) !void {
    for (left, 0..) |left_region, left_index| {
        for (right, 0..) |right_region, right_index| {
            if (same_collection and right_index <= left_index) continue;
            if (regionsOverlap(left_region, right_region))
                return error.OrganicNeighborErosionWorkspaceOverlap;
        }
    }
}

fn regionsOverlap(left: ByteRegion, right: ByteRegion) bool {
    if (left.byte_len == 0 or right.byte_len == 0) return false;
    const left_end = std.math.add(usize, left.start, left.byte_len) catch return true;
    const right_end = std.math.add(usize, right.start, right.byte_len) catch return true;
    return left.start < right_end and right.start < left_end;
}

fn copyState(destination: State, source: State) void {
    @memcpy(
        destination.microbial_by_class_group_component,
        source.microbial_by_class_group_component,
    );
    @memcpy(
        destination.residue_by_class_component,
        source.residue_by_class_component,
    );
    @memcpy(destination.adsorbed_by_class, source.adsorbed_by_class);
    @memcpy(destination.soil_by_class_fraction, source.soil_by_class_fraction);
}

fn subtractOrganicMatter(
    dimensions: OrganicDimensions,
    extents: Extents,
    candidate: State,
    contribution: OrganicMatterFlux,
) !void {
    try subtractSlice(
        candidate.microbial_by_class_group_component,
        contribution.microbial_by_class_group_component,
    );
    for (0..dimensions.stable_organic_class_count) |class_index| {
        const residue_start = class_index *
            dimensions.residue_component_count;
        const residue_end = residue_start +
            dimensions.residue_component_count;
        try subtractSlice(
            candidate.residue_by_class_component[residue_start..residue_end],
            contribution.residue_by_class_component[residue_start..residue_end],
        );
        try subtractValue(
            &candidate.adsorbed_by_class[class_index],
            contribution.adsorbed_by_class[class_index],
        );
        const soil_start = class_index *
            dimensions.soil_organic_fraction_count;
        const soil_end = soil_start +
            dimensions.soil_organic_fraction_count;
        try subtractSlice(
            candidate.soil_by_class_fraction[soil_start..soil_end],
            contribution.soil_by_class_fraction[soil_start..soil_end],
        );
    }
    std.debug.assert(extents.adsorbed == dimensions.stable_organic_class_count);
}

fn subtractSlice(candidate: anytype, contribution: anytype) !void {
    for (candidate, contribution) |*current, change|
        try subtractValue(current, change);
}

fn subtractValue(candidate: anytype, contribution: @TypeOf(candidate.*)) !void {
    inline for (@typeInfo(@TypeOf(candidate.*)).@"struct".fields) |field| {
        const result = @field(candidate.*, field.name) -
            @field(contribution, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteOrganicNeighborErosionResult;
        @field(candidate.*, field.name) = result;
    }
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn expectSlice(values: anytype, expected: f64) !void {
    for (values) |value| {
        inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
            try std.testing.expectEqual(expected, @field(value, field.name));
    }
}

const test_dimensions = OrganicDimensions{
    .microbial_substrate_class_count = 2,
    .microbial_group_count = 2,
    .microbial_component_count = 2,
    .stable_organic_class_count = 2,
    .residue_component_count = 2,
    .soil_organic_fraction_count = 2,
};

test "runtime dimensions subtract every connected organic inventory" {
    const microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 2),
    } ** 8;
    const residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 2),
    } ** 4;
    const adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 2),
    } ** 2;
    const soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 2)} ** 4;
    const sides = [_]NeighborSide{.{
        .local_total_sediment_Mg_per_step = 2,
        .positive_neighbor_total_sediment_Mg_per_step = 0,
        .positive_neighbor_organic_matter = .{
            .microbial_by_class_group_component = &microbial,
            .residue_by_class_component = &residue,
            .adsorbed_by_class = &adsorbed,
            .soil_by_class_fraction = &soil,
        },
        .connection = .connected,
    }};
    var state_microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 10),
    } ** 8;
    var state_residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 10),
    } ** 4;
    var state_adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 10),
    } ** 2;
    var state_soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 10)} ** 4;
    var scratch_microbial: [8]ElementalOrganicFlux = undefined;
    var scratch_residue: [4]ElementalOrganicFlux = undefined;
    var scratch_adsorbed: [2]AdsorbedOrganicFlux = undefined;
    var scratch_soil: [4]SoilOrganicFlux = undefined;
    var state = State{
        .microbial_by_class_group_component = &state_microbial,
        .residue_by_class_component = &state_residue,
        .adsorbed_by_class = &state_adsorbed,
        .soil_by_class_fraction = &state_soil,
    };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .dimensions = test_dimensions,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{
        .microbial_by_class_group_component = &scratch_microbial,
        .residue_by_class_component = &scratch_residue,
        .adsorbed_by_class = &scratch_adsorbed,
        .soil_by_class_fraction = &scratch_soil,
    });
    try expectSlice(&state_microbial, 8);
    try expectSlice(&state_residue, 8);
    try expectSlice(&state_adsorbed, 8);
    try expectSlice(&state_soil, 8);
}

test "strict activity and connection gates are independent" {
    const microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 4),
    } ** 8;
    const residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 4),
    } ** 4;
    const adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 4),
    } ** 2;
    const soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 4)} ** 4;
    const flux = OrganicMatterFlux{
        .microbial_by_class_group_component = &microbial,
        .residue_by_class_component = &residue,
        .adsorbed_by_class = &adsorbed,
        .soil_by_class_fraction = &soil,
    };
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_total_sediment_Mg_per_step = 2,
            .positive_neighbor_organic_matter = flux,
            .connection = .blocked,
        },
        .{
            .local_total_sediment_Mg_per_step = 1,
            .positive_neighbor_total_sediment_Mg_per_step = -1,
            .positive_neighbor_organic_matter = flux,
            .connection = .connected,
        },
        .{
            .local_total_sediment_Mg_per_step = 0,
            .positive_neighbor_total_sediment_Mg_per_step = -2,
            .positive_neighbor_organic_matter = flux,
            .connection = .connected,
        },
    };
    var state_microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 10),
    } ** 8;
    var state_residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 10),
    } ** 4;
    var state_adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 10),
    } ** 2;
    var state_soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 10)} ** 4;
    var scratch_microbial: [8]ElementalOrganicFlux = undefined;
    var scratch_residue: [4]ElementalOrganicFlux = undefined;
    var scratch_adsorbed: [2]AdsorbedOrganicFlux = undefined;
    var scratch_soil: [4]SoilOrganicFlux = undefined;
    var state = State{
        .microbial_by_class_group_component = &state_microbial,
        .residue_by_class_component = &state_residue,
        .adsorbed_by_class = &state_adsorbed,
        .soil_by_class_fraction = &state_soil,
    };

    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_Mg_per_step = 1,
        .dimensions = test_dimensions,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{
        .microbial_by_class_group_component = &scratch_microbial,
        .residue_by_class_component = &scratch_residue,
        .adsorbed_by_class = &scratch_adsorbed,
        .soil_by_class_fraction = &scratch_soil,
    });
    try expectSlice(&state_microbial, 6);
    try expectSlice(&state_residue, 6);
    try expectSlice(&state_adsorbed, 6);
    try expectSlice(&state_soil, 6);
}

test "shared face transfer conserves carbon nitrogen and phosphorus exactly" {
    const microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 3),
    } ** 8;
    const residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 3),
    } ** 4;
    const adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 3),
    } ** 2;
    const soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 3)} ** 4;
    var source_microbial = microbial;
    var source_residue = residue;
    var source_adsorbed = adsorbed;
    var source_soil = soil;
    const sides = [_]NeighborSide{.{
        .local_total_sediment_Mg_per_step = 0,
        .positive_neighbor_total_sediment_Mg_per_step = 2,
        .positive_neighbor_organic_matter = .{
            .microbial_by_class_group_component = &microbial,
            .residue_by_class_component = &residue,
            .adsorbed_by_class = &adsorbed,
            .soil_by_class_fraction = &soil,
        },
        .connection = .connected,
    }};
    var scratch_microbial: [8]ElementalOrganicFlux = undefined;
    var scratch_residue: [4]ElementalOrganicFlux = undefined;
    var scratch_adsorbed: [2]AdsorbedOrganicFlux = undefined;
    var scratch_soil: [4]SoilOrganicFlux = undefined;
    var state = State{
        .microbial_by_class_group_component = &source_microbial,
        .residue_by_class_component = &source_residue,
        .adsorbed_by_class = &source_adsorbed,
        .soil_by_class_fraction = &source_soil,
    };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .dimensions = test_dimensions,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{
        .microbial_by_class_group_component = &scratch_microbial,
        .residue_by_class_component = &scratch_residue,
        .adsorbed_by_class = &scratch_adsorbed,
        .soil_by_class_fraction = &scratch_soil,
    });
    try expectSlice(&source_microbial, 0);
    try expectSlice(&source_residue, 0);
    try expectSlice(&source_adsorbed, 0);
    try expectSlice(&source_soil, 0);
}

test "disabled and vertical modes bypass unused organic storage" {
    var state = State{
        .microbial_by_class_group_component = &.{},
        .residue_by_class_component = &.{},
        .adsorbed_by_class = &.{},
        .soil_by_class_fraction = &.{},
    };
    const workspace = Workspace{
        .microbial_by_class_group_component = &.{},
        .residue_by_class_component = &.{},
        .adsorbed_by_class = &.{},
        .soil_by_class_fraction = &.{},
    };
    try account(.{
        .disturbance_mode = .freeze_thaw,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .dimensions = std.mem.zeroes(OrganicDimensions),
        .neighbor_by_boundary_side = &.{},
    }, &state, workspace);
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .dimensions = std.mem.zeroes(OrganicDimensions),
        .neighbor_by_boundary_side = &.{},
    }, &state, workspace);
}

test "invalid dimensions alias and overflow preserve organic state atomically" {
    const microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, -std.math.floatMax(f64)),
    } ** 8;
    const residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 1),
    } ** 4;
    const adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 1),
    } ** 2;
    const soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 1)} ** 4;
    const sides = [_]NeighborSide{.{
        .local_total_sediment_Mg_per_step = 2,
        .positive_neighbor_total_sediment_Mg_per_step = 2,
        .positive_neighbor_organic_matter = .{
            .microbial_by_class_group_component = &microbial,
            .residue_by_class_component = &residue,
            .adsorbed_by_class = &adsorbed,
            .soil_by_class_fraction = &soil,
        },
        .connection = .connected,
    }};
    var state_microbial = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, std.math.floatMax(f64)),
    } ** 8;
    var state_residue = [_]ElementalOrganicFlux{
        filled(ElementalOrganicFlux, 5),
    } ** 4;
    var state_adsorbed = [_]AdsorbedOrganicFlux{
        filled(AdsorbedOrganicFlux, 5),
    } ** 2;
    var state_soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 5)} ** 4;
    var scratch_microbial: [8]ElementalOrganicFlux = undefined;
    var scratch_residue: [4]ElementalOrganicFlux = undefined;
    var scratch_adsorbed: [2]AdsorbedOrganicFlux = undefined;
    var scratch_soil: [4]SoilOrganicFlux = undefined;
    var state = State{
        .microbial_by_class_group_component = &state_microbial,
        .residue_by_class_component = &state_residue,
        .adsorbed_by_class = &state_adsorbed,
        .soil_by_class_fraction = &state_soil,
    };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .dimensions = test_dimensions,
        .neighbor_by_boundary_side = &sides,
    };
    const workspace = Workspace{
        .microbial_by_class_group_component = &scratch_microbial,
        .residue_by_class_component = &scratch_residue,
        .adsorbed_by_class = &scratch_adsorbed,
        .soil_by_class_fraction = &scratch_soil,
    };

    try std.testing.expectError(
        error.NonFiniteOrganicNeighborErosionResult,
        account(inputs, &state, workspace),
    );
    try expectSlice(&state_microbial, std.math.floatMax(f64));
    try expectSlice(&state_residue, 5);

    var invalid_inputs = inputs;
    invalid_inputs.dimensions.microbial_component_count = 3;
    try std.testing.expectError(
        error.OrganicNeighborErosionDimensionMismatch,
        account(invalid_inputs, &state, workspace),
    );
    try expectSlice(&state_microbial, std.math.floatMax(f64));

    const alias_workspace = Workspace{
        .microbial_by_class_group_component = &state_microbial,
        .residue_by_class_component = &scratch_residue,
        .adsorbed_by_class = &scratch_adsorbed,
        .soil_by_class_fraction = &scratch_soil,
    };
    try std.testing.expectError(
        error.OrganicNeighborErosionWorkspaceOverlap,
        account(inputs, &state, alias_workspace),
    );
    try expectSlice(&state_microbial, std.math.floatMax(f64));
}
