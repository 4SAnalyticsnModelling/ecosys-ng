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

pub const BoundarySide = enum {
    first,
    second,
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

pub const OppositeNeighborFlux = struct {
    total_sediment_megagrams_per_step: f64,
    organic_matter: OrganicMatterFlux,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    boundary_side: BoundarySide,
    sediment_activity_threshold_megagrams_per_step: f64,
    dimensions: OrganicDimensions,
    /// Null when the geometry-derived opposite-neighbor coordinate is absent.
    opposite_neighbor_first_side_flux: ?OppositeNeighborFlux,
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

/// Subtracts organic erosion from the first side's opposite neighbor.
///
/// Traceability: REDIST.F lines 3308--3334 under gates 3192--3193. Optional
/// geometry replaces `N4B,N5B`; `.first` replaces `NN == 1`. Every class,
/// group, component, and fraction extent is runtime configured. Flattening
/// preserves the Fortran traversal, including stable-class residue, adsorbed,
/// then soil order. Elemental inventories remain g C, g N, or g P per step.
pub fn account(inputs: Inputs, state: *State, workspace: Workspace) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical or
        inputs.boundary_side != .first)
    {
        return;
    }
    const flux = inputs.opposite_neighbor_first_side_flux orelse return;
    try validateDimensions(inputs.dimensions, flux.organic_matter, state.*, workspace);
    try validateInputs(
        inputs.sediment_activity_threshold_megagrams_per_step,
        flux,
        state.*,
        workspace,
    );
    if (@abs(flux.total_sediment_megagrams_per_step) <=
        inputs.sediment_activity_threshold_megagrams_per_step)
    {
        return;
    }

    copyState(workspace, state.*);
    try subtractOrganicMatter(
        inputs.dimensions,
        workspace,
        flux.organic_matter,
    );
    copyState(state.*, workspace);
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateDimensions(
    dimensions: OrganicDimensions,
    flux: OrganicMatterFlux,
    state: State,
    workspace: Workspace,
) !void {
    inline for (@typeInfo(OrganicDimensions).@"struct".fields) |field|
        if (@field(dimensions, field.name) == 0)
            return error.InvalidOppositeNeighborOrganicErosionDimensions;
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
    if (flux.microbial_by_class_group_component.len != extents.microbial or
        flux.residue_by_class_component.len != extents.residue or
        flux.adsorbed_by_class.len != extents.adsorbed or
        flux.soil_by_class_fraction.len != extents.soil)
    {
        return error.OppositeNeighborOrganicErosionDimensionMismatch;
    }
}

fn checkedProduct(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        error.InvalidOppositeNeighborOrganicErosionDimensions;
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
        return error.OppositeNeighborOrganicErosionDimensionMismatch;
    }
}

fn validateInputs(
    threshold: f64,
    flux: OppositeNeighborFlux,
    state: State,
    workspace: Workspace,
) !void {
    if (!std.math.isFinite(threshold) or
        !std.math.isFinite(flux.total_sediment_megagrams_per_step))
    {
        return error.NonFiniteOppositeNeighborOrganicErosionInput;
    }
    if (threshold < 0) return error.InvalidErosionActivityThreshold;
    try validateStateFinite(state);
    try validateFluxFinite(flux.organic_matter);
    try validateWorkspaceOwnership(workspace, state, flux.organic_matter);
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
                return error.NonFiniteOppositeNeighborOrganicErosionInput;
    }
}

fn validateWorkspaceOwnership(
    workspace: Workspace,
    state: State,
    flux: OrganicMatterFlux,
) !void {
    const workspace_regions = stateRegions(workspace);
    const state_regions = stateRegions(state);
    const flux_regions = fluxRegions(flux);
    try requireDisjoint(workspace_regions[0..], workspace_regions[0..], true);
    try requireDisjoint(workspace_regions[0..], state_regions[0..], false);
    try requireDisjoint(workspace_regions[0..], flux_regions[0..], false);
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
                return error.OppositeNeighborOrganicErosionWorkspaceOverlap;
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
            return error.NonFiniteOppositeNeighborOrganicErosionResult;
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

const TestStorage = struct {
    microbial: [8]ElementalOrganicFlux,
    residue: [4]ElementalOrganicFlux,
    adsorbed: [2]AdsorbedOrganicFlux,
    soil: [4]SoilOrganicFlux,

    fn init(value: f64) TestStorage {
        return .{
            .microbial = [_]ElementalOrganicFlux{
                filled(ElementalOrganicFlux, value),
            } ** 8,
            .residue = [_]ElementalOrganicFlux{
                filled(ElementalOrganicFlux, value),
            } ** 4,
            .adsorbed = [_]AdsorbedOrganicFlux{
                filled(AdsorbedOrganicFlux, value),
            } ** 2,
            .soil = [_]SoilOrganicFlux{
                filled(SoilOrganicFlux, value),
            } ** 4,
        };
    }

    fn state(self: *TestStorage) State {
        return .{
            .microbial_by_class_group_component = &self.microbial,
            .residue_by_class_component = &self.residue,
            .adsorbed_by_class = &self.adsorbed,
            .soil_by_class_fraction = &self.soil,
        };
    }

    fn flux(self: *const TestStorage) OrganicMatterFlux {
        return .{
            .microbial_by_class_group_component = &self.microbial,
            .residue_by_class_component = &self.residue,
            .adsorbed_by_class = &self.adsorbed,
            .soil_by_class_fraction = &self.soil,
        };
    }

    fn expectAll(self: *const TestStorage, expected: f64) !void {
        try expectSlice(&self.microbial, expected);
        try expectSlice(&self.residue, expected);
        try expectSlice(&self.adsorbed, expected);
        try expectSlice(&self.soil, expected);
    }
};

test "runtime dimensions subtract every opposite-neighbor organic inventory" {
    var contribution = TestStorage.init(2);
    var storage = TestStorage.init(10);
    var scratch = TestStorage.init(0);
    var state = storage.state();
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .dimensions = test_dimensions,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = -2,
            .organic_matter = contribution.flux(),
        },
    }, &state, scratch.state());
    try storage.expectAll(8);
}

test "shared face organic transfer conserves every elemental inventory" {
    var contribution = TestStorage.init(3);
    var source = TestStorage.init(3);
    var scratch = TestStorage.init(0);
    var state = source.state();
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .north_south,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .dimensions = test_dimensions,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = 2,
            .organic_matter = contribution.flux(),
        },
    }, &state, scratch.state());
    try source.expectAll(0);
}

test "strict sediment threshold bypasses large organic pool values" {
    var contribution = TestStorage.init(100);
    var storage = TestStorage.init(9);
    var scratch = TestStorage.init(0);
    var state = storage.state();
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .dimensions = test_dimensions,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = -1,
            .organic_matter = contribution.flux(),
        },
    }, &state, scratch.state());
    try storage.expectAll(9);
}

test "outer side axis and absent-geometry gates bypass unused storage" {
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
    const bypass = [_]Inputs{
        .{
            .disturbance_mode = .freeze_thaw,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .dimensions = std.mem.zeroes(OrganicDimensions),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .vertical,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .dimensions = std.mem.zeroes(OrganicDimensions),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .second,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .dimensions = std.mem.zeroes(OrganicDimensions),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .dimensions = std.mem.zeroes(OrganicDimensions),
            .opposite_neighbor_first_side_flux = null,
        },
    };
    for (bypass) |inputs| try account(inputs, &state, workspace);
}

test "dimension alias and overflow failures preserve state atomically" {
    var contribution = TestStorage.init(-std.math.floatMax(f64));
    var storage = TestStorage.init(std.math.floatMax(f64));
    var scratch = TestStorage.init(0);
    var state = storage.state();
    var inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .dimensions = test_dimensions,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = 2,
            .organic_matter = contribution.flux(),
        },
    };
    try std.testing.expectError(
        error.NonFiniteOppositeNeighborOrganicErosionResult,
        account(inputs, &state, scratch.state()),
    );
    try storage.expectAll(std.math.floatMax(f64));

    inputs.dimensions.microbial_component_count = 3;
    try std.testing.expectError(
        error.OppositeNeighborOrganicErosionDimensionMismatch,
        account(inputs, &state, scratch.state()),
    );
    try storage.expectAll(std.math.floatMax(f64));

    inputs.dimensions = test_dimensions;
    const alias_workspace = storage.state();
    try std.testing.expectError(
        error.OppositeNeighborOrganicErosionWorkspaceOverlap,
        account(inputs, &state, alias_workspace),
    );
    try storage.expectAll(std.math.floatMax(f64));
}
