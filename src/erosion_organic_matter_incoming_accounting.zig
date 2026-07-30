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
    /// [microbial_substrate_class][microbial_group][microbial_component].
    microbial_by_class_group_component: []const ElementalOrganicFlux,
    /// [stable_organic_class][residue_component].
    residue_by_class_component: []const ElementalOrganicFlux,
    /// [stable_organic_class].
    adsorbed_by_class: []const AdsorbedOrganicFlux,
    /// [stable_organic_class][soil_organic_fraction].
    soil_by_class_fraction: []const SoilOrganicFlux,
};

pub const BoundaryFlux = struct {
    total_sediment_Mg_per_step: f64,
    organic_matter: OrganicMatterFlux,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_Mg_per_step: f64,
    dimensions: OrganicDimensions,
    local_flux_by_boundary_side: []const BoundaryFlux,
    positive_neighbor_total_sediment_Mg_per_step_by_boundary_side: []const f64,
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

/// Accounts organic-matter erosion inflow on one horizontal axis.
///
/// Traceability: REDIST.F lines 3012--3038 under gates 2888--2895. All model
/// extents are runtime inputs. Flattened storage retains Fortran traversal:
/// microbial component `M` fastest inside group `NO` and class `K`; residue
/// component and soil fraction `M` fastest inside class `K`. Carbon, nitrogen,
/// phosphorus, colonized carbon, and acetate remain on their elemental gram
/// basis. State commits atomically after source-ordered finite evaluation.
pub fn account(
    inputs: Inputs,
    state: *State,
    workspace: Workspace,
) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical)
    {
        return;
    }
    const extents = try validateDimensions(inputs, state.*, workspace);
    try validateInputs(inputs, state.*, workspace, extents);
    copyState(workspace, state.*);

    for (
        inputs.local_flux_by_boundary_side,
        inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side,
    ) |local, positive_neighbor_total| {
        if (@abs(local.total_sediment_Mg_per_step) <=
            inputs.sediment_activity_threshold_Mg_per_step and
            @abs(positive_neighbor_total) <=
                inputs.sediment_activity_threshold_Mg_per_step)
        {
            continue;
        }
        try addOrganicMatter(workspace, local.organic_matter);
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
    const dimensions = inputs.dimensions;
    inline for (@typeInfo(OrganicDimensions).@"struct".fields) |field|
        if (@field(dimensions, field.name) == 0)
            return error.InvalidOrganicErosionDimensions;
    if (inputs.local_flux_by_boundary_side.len == 0)
        return error.InvalidOrganicErosionDimensions;
    if (inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side.len !=
        inputs.local_flux_by_boundary_side.len)
    {
        return error.OrganicErosionDimensionMismatch;
    }
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
    for (inputs.local_flux_by_boundary_side) |flux|
        try validateFluxDimensions(flux.organic_matter, extents);
    return extents;
}

fn checkedProduct(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        error.InvalidOrganicErosionDimensions;
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
        return error.OrganicErosionDimensionMismatch;
    }
}

fn validateFluxDimensions(flux: OrganicMatterFlux, extents: Extents) !void {
    if (flux.microbial_by_class_group_component.len != extents.microbial or
        flux.residue_by_class_component.len != extents.residue or
        flux.adsorbed_by_class.len != extents.adsorbed or
        flux.soil_by_class_fraction.len != extents.soil)
    {
        return error.OrganicErosionDimensionMismatch;
    }
}

fn validateInputs(
    inputs: Inputs,
    state: State,
    workspace: Workspace,
    extents: Extents,
) !void {
    _ = extents;
    if (!std.math.isFinite(inputs.sediment_activity_threshold_Mg_per_step))
        return error.NonFiniteOrganicErosionInput;
    if (inputs.sediment_activity_threshold_Mg_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateStateFinite(state);
    try validateWorkspaceOwnership(workspace, state, inputs);
    for (inputs.local_flux_by_boundary_side) |flux| {
        if (!std.math.isFinite(flux.total_sediment_Mg_per_step))
            return error.NonFiniteOrganicErosionInput;
        try validateFluxFinite(flux.organic_matter);
    }
    for (inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteOrganicErosionInput;
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
                return error.NonFiniteOrganicErosionInput;
    }
}

fn validateWorkspaceOwnership(
    workspace: Workspace,
    state: State,
    inputs: Inputs,
) !void {
    if (overlapElemental(
        workspace.microbial_by_class_group_component,
        state.microbial_by_class_group_component,
    ) or overlapElemental(
        workspace.microbial_by_class_group_component,
        state.residue_by_class_component,
    ) or overlapElemental(
        workspace.residue_by_class_component,
        state.microbial_by_class_group_component,
    ) or overlapElemental(
        workspace.residue_by_class_component,
        state.residue_by_class_component,
    ) or overlapElemental(
        workspace.microbial_by_class_group_component,
        workspace.residue_by_class_component,
    ) or overlapAdsorbed(
        workspace.adsorbed_by_class,
        state.adsorbed_by_class,
    ) or overlapSoil(
        workspace.soil_by_class_fraction,
        state.soil_by_class_fraction,
    )) return error.OrganicErosionWorkspaceOverlap;

    for (inputs.local_flux_by_boundary_side) |flux| {
        if (overlapElemental(
            workspace.microbial_by_class_group_component,
            flux.organic_matter.microbial_by_class_group_component,
        ) or overlapElemental(
            workspace.microbial_by_class_group_component,
            flux.organic_matter.residue_by_class_component,
        ) or overlapElemental(
            workspace.residue_by_class_component,
            flux.organic_matter.microbial_by_class_group_component,
        ) or overlapElemental(
            workspace.residue_by_class_component,
            flux.organic_matter.residue_by_class_component,
        ) or overlapAdsorbed(
            workspace.adsorbed_by_class,
            flux.organic_matter.adsorbed_by_class,
        ) or overlapSoil(
            workspace.soil_by_class_fraction,
            flux.organic_matter.soil_by_class_fraction,
        )) return error.OrganicErosionWorkspaceOverlap;
    }
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

fn addOrganicMatter(candidate: State, contribution: OrganicMatterFlux) !void {
    try addSlice(
        candidate.microbial_by_class_group_component,
        contribution.microbial_by_class_group_component,
    );
    try addSlice(
        candidate.residue_by_class_component,
        contribution.residue_by_class_component,
    );
    try addSlice(candidate.adsorbed_by_class, contribution.adsorbed_by_class);
    try addSlice(
        candidate.soil_by_class_fraction,
        contribution.soil_by_class_fraction,
    );
}

fn addSlice(candidate: anytype, contribution: anytype) !void {
    for (candidate, contribution) |*current, change| {
        inline for (@typeInfo(@TypeOf(current.*)).@"struct".fields) |field| {
            const result = @field(current.*, field.name) +
                @field(change, field.name);
            if (!std.math.isFinite(result))
                return error.NonFiniteOrganicErosionResult;
            @field(current.*, field.name) = result;
        }
    }
}

fn overlapElemental(
    left: []const ElementalOrganicFlux,
    right: []const ElementalOrganicFlux,
) bool {
    return overlapBytes(
        @intFromPtr(left.ptr),
        left.len,
        @sizeOf(ElementalOrganicFlux),
        @intFromPtr(right.ptr),
        right.len,
        @sizeOf(ElementalOrganicFlux),
    );
}

fn overlapAdsorbed(
    left: []const AdsorbedOrganicFlux,
    right: []const AdsorbedOrganicFlux,
) bool {
    return overlapBytes(
        @intFromPtr(left.ptr),
        left.len,
        @sizeOf(AdsorbedOrganicFlux),
        @intFromPtr(right.ptr),
        right.len,
        @sizeOf(AdsorbedOrganicFlux),
    );
}

fn overlapSoil(
    left: []const SoilOrganicFlux,
    right: []const SoilOrganicFlux,
) bool {
    return overlapBytes(
        @intFromPtr(left.ptr),
        left.len,
        @sizeOf(SoilOrganicFlux),
        @intFromPtr(right.ptr),
        right.len,
        @sizeOf(SoilOrganicFlux),
    );
}

fn overlapBytes(
    left_start: usize,
    left_len: usize,
    left_item_size: usize,
    right_start: usize,
    right_len: usize,
    right_item_size: usize,
) bool {
    if (left_len == 0 or right_len == 0) return false;
    const left_bytes = std.math.mul(usize, left_len, left_item_size) catch return true;
    const right_bytes = std.math.mul(usize, right_len, right_item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
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

test "source default extents add every organic inventory" {
    const dimensions = OrganicDimensions{
        .microbial_substrate_class_count = 6,
        .microbial_group_count = 7,
        .microbial_component_count = 3,
        .stable_organic_class_count = 5,
        .residue_component_count = 2,
        .soil_organic_fraction_count = 5,
    };
    const microbial = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 2)} ** 126;
    const residue = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 2)} ** 10;
    const adsorbed = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, 2)} ** 5;
    const soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 2)} ** 25;
    const local = [_]BoundaryFlux{.{
        .total_sediment_Mg_per_step = 2,
        .organic_matter = .{
            .microbial_by_class_group_component = &microbial,
            .residue_by_class_component = &residue,
            .adsorbed_by_class = &adsorbed,
            .soil_by_class_fraction = &soil,
        },
    }};
    const positive = [_]f64{0};
    var state_microbial = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 100)} ** 126;
    var state_residue = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 100)} ** 10;
    var state_adsorbed = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, 100)} ** 5;
    var state_soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 100)} ** 25;
    var scratch_microbial: [126]ElementalOrganicFlux = undefined;
    var scratch_residue: [10]ElementalOrganicFlux = undefined;
    var scratch_adsorbed: [5]AdsorbedOrganicFlux = undefined;
    var scratch_soil: [25]SoilOrganicFlux = undefined;
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
        .dimensions = dimensions,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    }, &state, .{
        .microbial_by_class_group_component = &scratch_microbial,
        .residue_by_class_component = &scratch_residue,
        .adsorbed_by_class = &scratch_adsorbed,
        .soil_by_class_fraction = &scratch_soil,
    });
    try expectSlice(state.microbial_by_class_group_component, 102);
    try expectSlice(state.residue_by_class_component, 102);
    try expectSlice(state.adsorbed_by_class, 102);
    try expectSlice(state.soil_by_class_fraction, 102);
}

test "nonlegacy runtime extents conserve every accepted side exactly" {
    const dimensions = OrganicDimensions{
        .microbial_substrate_class_count = 2,
        .microbial_group_count = 4,
        .microbial_component_count = 2,
        .stable_organic_class_count = 3,
        .residue_component_count = 4,
        .soil_organic_fraction_count = 6,
    };
    const microbial_count = 16;
    const residue_count = 12;
    const adsorbed_count = 3;
    const soil_count = 18;
    const two_m = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 2)} ** microbial_count;
    const two_r = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 2)} ** residue_count;
    const two_a = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, 2)} ** adsorbed_count;
    const two_s = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 2)} ** soil_count;
    const three_m = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, -3)} ** microbial_count;
    const three_r = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, -3)} ** residue_count;
    const three_a = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, -3)} ** adsorbed_count;
    const three_s = [_]SoilOrganicFlux{filled(SoilOrganicFlux, -3)} ** soil_count;
    const four_m = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 4)} ** microbial_count;
    const four_r = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 4)} ** residue_count;
    const four_a = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, 4)} ** adsorbed_count;
    const four_s = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 4)} ** soil_count;
    const local = [_]BoundaryFlux{
        makeBoundary(1, &two_m, &two_r, &two_a, &two_s),
        makeBoundary(-3, &three_m, &three_r, &three_a, &three_s),
        makeBoundary(4, &four_m, &four_r, &four_a, &four_s),
    };
    const positive = [_]f64{ 1, 0, 0 };
    var state_m = [_]ElementalOrganicFlux{.{}} ** microbial_count;
    var state_r = [_]ElementalOrganicFlux{.{}} ** residue_count;
    var state_a = [_]AdsorbedOrganicFlux{.{}} ** adsorbed_count;
    var state_s = [_]SoilOrganicFlux{.{}} ** soil_count;
    var scratch_m: [microbial_count]ElementalOrganicFlux = undefined;
    var scratch_r: [residue_count]ElementalOrganicFlux = undefined;
    var scratch_a: [adsorbed_count]AdsorbedOrganicFlux = undefined;
    var scratch_s: [soil_count]SoilOrganicFlux = undefined;
    var state = State{
        .microbial_by_class_group_component = &state_m,
        .residue_by_class_component = &state_r,
        .adsorbed_by_class = &state_a,
        .soil_by_class_fraction = &state_s,
    };
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_Mg_per_step = 1,
        .dimensions = dimensions,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    }, &state, .{
        .microbial_by_class_group_component = &scratch_m,
        .residue_by_class_component = &scratch_r,
        .adsorbed_by_class = &scratch_a,
        .soil_by_class_fraction = &scratch_s,
    });
    try expectSlice(state.microbial_by_class_group_component, 1);
    try expectSlice(state.residue_by_class_component, 1);
    try expectSlice(state.adsorbed_by_class, 1);
    try expectSlice(state.soil_by_class_fraction, 1);
}

fn makeBoundary(
    sediment: f64,
    microbial: []const ElementalOrganicFlux,
    residue: []const ElementalOrganicFlux,
    adsorbed: []const AdsorbedOrganicFlux,
    soil: []const SoilOrganicFlux,
) BoundaryFlux {
    return .{
        .total_sediment_Mg_per_step = sediment,
        .organic_matter = .{
            .microbial_by_class_group_component = microbial,
            .residue_by_class_component = residue,
            .adsorbed_by_class = adsorbed,
            .soil_by_class_fraction = soil,
        },
    };
}

test "disabled and vertical modes bypass unused runtime layouts" {
    var state = State{
        .microbial_by_class_group_component = &.{},
        .residue_by_class_component = &.{},
        .adsorbed_by_class = &.{},
        .soil_by_class_fraction = &.{},
    };
    const empty_workspace = Workspace{
        .microbial_by_class_group_component = &.{},
        .residue_by_class_component = &.{},
        .adsorbed_by_class = &.{},
        .soil_by_class_fraction = &.{},
    };
    const empty_dimensions = OrganicDimensions{
        .microbial_substrate_class_count = 0,
        .microbial_group_count = 0,
        .microbial_component_count = 0,
        .stable_organic_class_count = 0,
        .residue_component_count = 0,
        .soil_organic_fraction_count = 0,
    };
    try account(.{
        .disturbance_mode = .freeze_thaw,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .dimensions = empty_dimensions,
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &.{},
    }, &state, empty_workspace);
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .dimensions = empty_dimensions,
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &.{},
    }, &state, empty_workspace);
}

test "late invalid input overflow dimensions and alias preserve state" {
    const dimensions = OrganicDimensions{
        .microbial_substrate_class_count = 1,
        .microbial_group_count = 1,
        .microbial_component_count = 1,
        .stable_organic_class_count = 1,
        .residue_component_count = 1,
        .soil_organic_fraction_count = 1,
    };
    var microbial = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 2)};
    var residue = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 2)};
    var adsorbed = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, 2)};
    var soil = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 2)};
    soil[0].phosphorus_g_p_per_step = std.math.nan(f64);
    const local = [_]BoundaryFlux{
        makeBoundary(2, &microbial, &residue, &adsorbed, &soil),
    };
    const positive = [_]f64{0};
    var state_m = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 5)};
    var state_r = [_]ElementalOrganicFlux{filled(ElementalOrganicFlux, 5)};
    var state_a = [_]AdsorbedOrganicFlux{filled(AdsorbedOrganicFlux, 5)};
    var state_s = [_]SoilOrganicFlux{filled(SoilOrganicFlux, 5)};
    var scratch_m: [1]ElementalOrganicFlux = undefined;
    var scratch_r: [1]ElementalOrganicFlux = undefined;
    var scratch_a: [1]AdsorbedOrganicFlux = undefined;
    var scratch_s: [1]SoilOrganicFlux = undefined;
    var state = State{
        .microbial_by_class_group_component = &state_m,
        .residue_by_class_component = &state_r,
        .adsorbed_by_class = &state_a,
        .soil_by_class_fraction = &state_s,
    };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .dimensions = dimensions,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    };
    const workspace = Workspace{
        .microbial_by_class_group_component = &scratch_m,
        .residue_by_class_component = &scratch_r,
        .adsorbed_by_class = &scratch_a,
        .soil_by_class_fraction = &scratch_s,
    };
    try std.testing.expectError(
        error.NonFiniteOrganicErosionInput,
        account(inputs, &state, workspace),
    );
    try expectSlice(state.soil_by_class_fraction, 5);

    soil[0] = filled(SoilOrganicFlux, std.math.floatMax(f64));
    state_s[0] = filled(SoilOrganicFlux, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteOrganicErosionResult,
        account(inputs, &state, workspace),
    );
    try expectSlice(state.soil_by_class_fraction, std.math.floatMax(f64));

    try std.testing.expectError(
        error.OrganicErosionWorkspaceOverlap,
        account(inputs, &state, .{
            .microbial_by_class_group_component = &state_m,
            .residue_by_class_component = &scratch_r,
            .adsorbed_by_class = &scratch_a,
            .soil_by_class_fraction = &scratch_s,
        }),
    );
    var short_soil = [_]SoilOrganicFlux{};
    try std.testing.expectError(
        error.OrganicErosionDimensionMismatch,
        account(.{
            .disturbance_mode = inputs.disturbance_mode,
            .transport_axis = inputs.transport_axis,
            .sediment_activity_threshold_Mg_per_step = inputs.sediment_activity_threshold_Mg_per_step,
            .dimensions = .{
                .microbial_substrate_class_count = 1,
                .microbial_group_count = 1,
                .microbial_component_count = 1,
                .stable_organic_class_count = 2,
                .residue_component_count = 1,
                .soil_organic_fraction_count = 1,
            },
            .local_flux_by_boundary_side = inputs.local_flux_by_boundary_side,
            .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side,
        }, &state, .{
            .microbial_by_class_group_component = &scratch_m,
            .residue_by_class_component = &scratch_r,
            .adsorbed_by_class = &scratch_a,
            .soil_by_class_fraction = &short_soil,
        }),
    );
    try expectSlice(state.microbial_by_class_group_component, 5);
}
