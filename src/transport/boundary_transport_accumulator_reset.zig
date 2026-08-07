const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const DisturbanceMode = enum {
    no_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter_change,
    freeze_thaw_erosion_and_organic_matter_change,

    fn includesErosion(self: DisturbanceMode) bool {
        return self == .freeze_thaw_and_erosion or
            self == .freeze_thaw_erosion_and_organic_matter_change;
    }
};

/// Runtime extents replace REDIST's `JS`, `K`, `NO`, and `M` bounds.
pub const Dimensions = struct {
    grid_cell_count: usize,
    snow_layer_count: usize,
    runoff_organic_class_count: usize,
    microbial_class_count: usize,
    microbial_group_count: usize,
    microbial_component_count: usize,
    residue_and_soil_organic_class_count: usize,
    residue_component_count: usize,
    soil_organic_component_count: usize,
};

pub const SurfaceWaterHeatFlux = struct {
    runoff_water_m3_per_step: f64 = 0,
    runoff_heat_megajoules_per_step: f64 = 0,
    drifting_snow_m3_per_step: f64 = 0,
    drifting_snow_liquid_water_m3_per_step: f64 = 0,
    drifting_snow_ice_m3_per_step: f64 = 0,
    drifting_snow_heat_megajoules_per_step: f64 = 0,
};

pub const OrganicElementFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
};

pub const RunoffOrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
    acetate_carbon_g_c_per_step: f64 = 0,
};

pub const RunoffElementGasFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    hydrogen_g_h_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    other_nitrogen_g_n_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    snow_carbon_dioxide_g_c_per_step: f64 = 0,
    snow_methane_g_c_per_step: f64 = 0,
    snow_oxygen_g_o_per_step: f64 = 0,
    snow_dinitrogen_g_n_per_step: f64 = 0,
    snow_nitrous_oxide_g_n_per_step: f64 = 0,
    snow_ammonium_g_n_per_step: f64 = 0,
    snow_ammonia_g_n_per_step: f64 = 0,
    snow_nitrate_g_n_per_step: f64 = 0,
    snow_h2po4_g_p_per_step: f64 = 0,
    snow_hpo4_g_p_per_step: f64 = 0,
};

pub const SnowLayerElementGasFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
};

/// Source order of redist.f lines 1491--1532. Each flattened state entry is
/// mol/step. The runoff-only hydrogen sulfate entry is retained explicitly.
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
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
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
    iron_hpo4,
    iron_h2po4,
    calcium_po4,
    calcium_hpo4,
    calcium_h2po4,
    magnesium_hpo4,
};

/// Source order shared by drifting-snow and snow-layer salt resets at REDIST
/// lines 1533--1618. Each flattened state entry is mol/step.
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
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
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
    iron_hpo4,
    iron_h2po4,
    calcium_po4,
    calcium_hpo4,
    calcium_h2po4,
    magnesium_hpo4,
};

pub const ErodedMineralFlux = struct {
    sediment_megagrams_per_step: f64 = 0,
    sand_megagrams_per_step: f64 = 0,
    silt_megagrams_per_step: f64 = 0,
    clay_megagrams_per_step: f64 = 0,
    cation_exchange_capacity_mol_per_step: f64 = 0,
    anion_exchange_capacity_mol_per_step: f64 = 0,
};

pub const ErodedNitrogenFlux = struct {
    fertilizer_ammonium_non_band_mol_n_per_step: f64 = 0,
    fertilizer_ammonia_non_band_mol_n_per_step: f64 = 0,
    fertilizer_urea_non_band_mol_n_per_step: f64 = 0,
    fertilizer_nitrate_non_band_mol_n_per_step: f64 = 0,
    fertilizer_ammonium_band_mol_n_per_step: f64 = 0,
    fertilizer_ammonia_band_mol_n_per_step: f64 = 0,
    fertilizer_urea_band_mol_n_per_step: f64 = 0,
    fertilizer_nitrate_band_mol_n_per_step: f64 = 0,
    adsorbed_ammonium_non_band_mol_n_per_step: f64 = 0,
    adsorbed_ammonium_band_mol_n_per_step: f64 = 0,
};

pub const ErodedAdsorbedSaltFlux = struct {
    hydrogen_mol_per_step: f64 = 0,
    aluminum_mol_per_step: f64 = 0,
    iron_mol_per_step: f64 = 0,
    calcium_mol_per_step: f64 = 0,
    magnesium_mol_per_step: f64 = 0,
    sodium_mol_per_step: f64 = 0,
    potassium_mol_per_step: f64 = 0,
    bicarbonate_mol_per_step: f64 = 0,
    aluminum_hydroxide_2_mol_per_step: f64 = 0,
    iron_hydroxide_2_mol_per_step: f64 = 0,
    deprotonated_site_non_band_mol_per_step: f64 = 0,
    hydroxyl_site_non_band_mol_per_step: f64 = 0,
    protonated_site_non_band_mol_per_step: f64 = 0,
    deprotonated_site_band_mol_per_step: f64 = 0,
    hydroxyl_site_band_mol_per_step: f64 = 0,
    protonated_site_band_mol_per_step: f64 = 0,
};

pub const ErodedMatrixPrecipitateFlux = struct {
    aluminum_hydroxide_mol_per_step: f64 = 0,
    iron_hydroxide_mol_per_step: f64 = 0,
    calcium_carbonate_mol_per_step: f64 = 0,
    calcium_sulfate_mol_per_step: f64 = 0,
};

pub const ErodedSulfateComplexFlux = struct {
    aluminum_sulfate_soil_mol_per_step: f64 = 0,
    iron_sulfate_soil_mol_per_step: f64 = 0,
    calcium_sulfate_soil_mol_per_step: f64 = 0,
    magnesium_sulfate_soil_mol_per_step: f64 = 0,
    sodium_sulfate_soil_mol_per_step: f64 = 0,
    potassium_sulfate_soil_mol_per_step: f64 = 0,
    aluminum_sulfate_fertilizer_mol_per_step: f64 = 0,
    iron_sulfate_fertilizer_mol_per_step: f64 = 0,
    calcium_sulfate_fertilizer_mol_per_step: f64 = 0,
    magnesium_sulfate_fertilizer_mol_per_step: f64 = 0,
    sodium_sulfate_fertilizer_mol_per_step: f64 = 0,
    potassium_sulfate_fertilizer_mol_per_step: f64 = 0,
};

pub const ErodedPhosphorusFlux = struct {
    adsorbed_hpo4_non_band_mol_p_per_step: f64 = 0,
    adsorbed_h2po4_non_band_mol_p_per_step: f64 = 0,
    adsorbed_hpo4_band_mol_p_per_step: f64 = 0,
    adsorbed_h2po4_band_mol_p_per_step: f64 = 0,
    aluminum_phosphate_non_band_mol_p_per_step: f64 = 0,
    iron_phosphate_non_band_mol_p_per_step: f64 = 0,
    calcium_hpo4_non_band_mol_p_per_step: f64 = 0,
    calcium_h2po4_non_band_mol_p_per_step: f64 = 0,
    apatite_non_band_mol_p_per_step: f64 = 0,
    aluminum_phosphate_band_mol_p_per_step: f64 = 0,
    iron_phosphate_band_mol_p_per_step: f64 = 0,
    calcium_hpo4_band_mol_p_per_step: f64 = 0,
    calcium_h2po4_band_mol_p_per_step: f64 = 0,
    apatite_band_mol_p_per_step: f64 = 0,
};

pub const ErosionScalarFlux = struct {
    mineral: ErodedMineralFlux = .{},
    nitrogen: ErodedNitrogenFlux = .{},
    adsorbed_salt: ErodedAdsorbedSaltFlux = .{},
    matrix_precipitate: ErodedMatrixPrecipitateFlux = .{},
    sulfate_complex: ErodedSulfateComplexFlux = .{},
    phosphorus: ErodedPhosphorusFlux = .{},
};

pub const AdsorbedOrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
    acetate_carbon_g_c_per_step: f64 = 0,
};

pub const SoilOrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    ash_g_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
};

pub const SnowpackWaterHeatFlux = struct {
    snow_m3_per_step: f64 = 0,
    liquid_water_m3_per_step: f64 = 0,
    water_vapor_m3_per_step: f64 = 0,
    ice_m3_per_step: f64 = 0,
    convective_heat_megajoules_per_step: f64 = 0,
};

/// All storage is caller-owned and runtime allocated. Flattened indexes are
/// cell-major; scientific subindexes retain the REDIST loop order.
pub const State = struct {
    surface_water_heat_by_cell: []SurfaceWaterHeatFlux,
    runoff_organic_by_cell_class: []RunoffOrganicFlux,
    runoff_element_gas_by_cell: []RunoffElementGasFlux,
    snow_layer_element_gas_by_cell_layer: []SnowLayerElementGasFlux,
    runoff_salt_mol_by_cell_species: []f64,
    snow_drift_salt_mol_by_cell_species: []f64,
    snow_layer_salt_mol_by_cell_layer_species: []f64,
    erosion_scalar_by_cell: []ErosionScalarFlux,
    erosion_microbial_by_cell_class_group_component: []OrganicElementFlux,
    erosion_residue_by_cell_class_component: []OrganicElementFlux,
    erosion_adsorbed_by_cell_class: []AdsorbedOrganicFlux,
    erosion_soil_organic_by_cell_class_component: []SoilOrganicFlux,
    snowpack_water_heat_by_cell_layer: []SnowpackWaterHeatFlux,
};

const Lengths = struct {
    cell: usize,
    runoff_organic: usize,
    snow_layer: usize,
    runoff_salt: usize,
    snow_drift_salt: usize,
    snow_layer_salt: usize,
    microbial: usize,
    residue: usize,
    adsorbed: usize,
    soil_organic: usize,
};

/// Resets REDIST's per-cell net transport accumulators.
///
/// Traceability: REDIST.F lines 1441--1721. Static salt equilibrium preserves
/// salt storage, and disturbance modes without erosion preserve all erosion
/// storage, matching the two source gates. Dimension and finite-value checks
/// complete before the first mutation, so failure is atomic.
pub fn reset(
    dimensions: Dimensions,
    salt_equilibrium_mode: SaltEquilibriumMode,
    disturbance_mode: DisturbanceMode,
    state: *State,
) !void {
    const lengths = try calculateLengths(dimensions);
    try validateLengths(state.*, lengths);
    try validateTargetValues(state.*, salt_equilibrium_mode, disturbance_mode);

    zeroSlice(state.surface_water_heat_by_cell);
    zeroSlice(state.runoff_organic_by_cell_class);
    zeroSlice(state.runoff_element_gas_by_cell);
    zeroSlice(state.snow_layer_element_gas_by_cell_layer);

    if (salt_equilibrium_mode == .dynamic) {
        @memset(state.runoff_salt_mol_by_cell_species, 0);
        @memset(state.snow_drift_salt_mol_by_cell_species, 0);
        @memset(state.snow_layer_salt_mol_by_cell_layer_species, 0);
    }

    if (disturbance_mode.includesErosion()) {
        zeroSlice(state.erosion_scalar_by_cell);
        zeroSlice(state.erosion_microbial_by_cell_class_group_component);
        zeroSlice(state.erosion_residue_by_cell_class_component);
        zeroSlice(state.erosion_adsorbed_by_cell_class);
        zeroSlice(state.erosion_soil_organic_by_cell_class_component);
    }

    zeroSlice(state.snowpack_water_heat_by_cell_layer);
}

fn calculateLengths(dimensions: Dimensions) !Lengths {
    inline for (@typeInfo(Dimensions).@"struct".fields) |field| {
        if (@field(dimensions, field.name) == 0) return error.InvalidAccumulatorDimensions;
    }
    const cell_snow = try multiply(dimensions.grid_cell_count, dimensions.snow_layer_count);
    const cell_organic = try multiply(
        dimensions.grid_cell_count,
        dimensions.runoff_organic_class_count,
    );
    const microbial_class_group = try multiply(
        dimensions.microbial_class_count,
        dimensions.microbial_group_count,
    );
    const microbial_per_cell = try multiply(
        microbial_class_group,
        dimensions.microbial_component_count,
    );
    const organic_class_count = dimensions.residue_and_soil_organic_class_count;
    return .{
        .cell = dimensions.grid_cell_count,
        .runoff_organic = cell_organic,
        .snow_layer = cell_snow,
        .runoff_salt = try multiply(
            dimensions.grid_cell_count,
            @typeInfo(RunoffSaltSpecies).@"enum".fields.len,
        ),
        .snow_drift_salt = try multiply(
            dimensions.grid_cell_count,
            @typeInfo(SnowSaltSpecies).@"enum".fields.len,
        ),
        .snow_layer_salt = try multiply(
            cell_snow,
            @typeInfo(SnowSaltSpecies).@"enum".fields.len,
        ),
        .microbial = try multiply(dimensions.grid_cell_count, microbial_per_cell),
        .residue = try multiply(
            dimensions.grid_cell_count,
            try multiply(organic_class_count, dimensions.residue_component_count),
        ),
        .adsorbed = try multiply(dimensions.grid_cell_count, organic_class_count),
        .soil_organic = try multiply(
            dimensions.grid_cell_count,
            try multiply(organic_class_count, dimensions.soil_organic_component_count),
        ),
    };
}

fn multiply(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        return error.AccumulatorDimensionOverflow;
}

fn validateLengths(state: State, lengths: Lengths) !void {
    if (state.surface_water_heat_by_cell.len != lengths.cell or
        state.runoff_organic_by_cell_class.len != lengths.runoff_organic or
        state.runoff_element_gas_by_cell.len != lengths.cell or
        state.snow_layer_element_gas_by_cell_layer.len != lengths.snow_layer or
        state.runoff_salt_mol_by_cell_species.len != lengths.runoff_salt or
        state.snow_drift_salt_mol_by_cell_species.len != lengths.snow_drift_salt or
        state.snow_layer_salt_mol_by_cell_layer_species.len != lengths.snow_layer_salt or
        state.erosion_scalar_by_cell.len != lengths.cell or
        state.erosion_microbial_by_cell_class_group_component.len != lengths.microbial or
        state.erosion_residue_by_cell_class_component.len != lengths.residue or
        state.erosion_adsorbed_by_cell_class.len != lengths.adsorbed or
        state.erosion_soil_organic_by_cell_class_component.len != lengths.soil_organic or
        state.snowpack_water_heat_by_cell_layer.len != lengths.snow_layer)
    {
        return error.AccumulatorDimensionMismatch;
    }
}

fn validateTargetValues(
    state: State,
    salt_equilibrium_mode: SaltEquilibriumMode,
    disturbance_mode: DisturbanceMode,
) !void {
    try validateFiniteSlice(state.surface_water_heat_by_cell);
    try validateFiniteSlice(state.runoff_organic_by_cell_class);
    try validateFiniteSlice(state.runoff_element_gas_by_cell);
    try validateFiniteSlice(state.snow_layer_element_gas_by_cell_layer);
    if (salt_equilibrium_mode == .dynamic) {
        try validateFiniteSlice(state.runoff_salt_mol_by_cell_species);
        try validateFiniteSlice(state.snow_drift_salt_mol_by_cell_species);
        try validateFiniteSlice(state.snow_layer_salt_mol_by_cell_layer_species);
    }
    if (disturbance_mode.includesErosion()) {
        try validateFiniteSlice(state.erosion_scalar_by_cell);
        try validateFiniteSlice(state.erosion_microbial_by_cell_class_group_component);
        try validateFiniteSlice(state.erosion_residue_by_cell_class_component);
        try validateFiniteSlice(state.erosion_adsorbed_by_cell_class);
        try validateFiniteSlice(state.erosion_soil_organic_by_cell_class_component);
    }
    try validateFiniteSlice(state.snowpack_water_heat_by_cell_layer);
}

fn validateFiniteSlice(values: anytype) !void {
    for (values) |value| try validateFinite(value);
}

fn validateFinite(value: anytype) !void {
    switch (@typeInfo(@TypeOf(value))) {
        .float => if (!std.math.isFinite(value))
            return error.NonFiniteTransportAccumulator,
        .@"struct" => inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
            try validateFinite(@field(value, field.name)),
        else => @compileError("transport accumulator must contain only floats or structs"),
    }
}

fn zeroSlice(values: anytype) void {
    for (values) |*value| value.* = std.mem.zeroes(@TypeOf(value.*));
}

fn setValue(target: anytype, value: f64) void {
    const Child = @typeInfo(@TypeOf(target)).pointer.child;
    switch (@typeInfo(Child)) {
        .float => target.* = value,
        .@"struct" => inline for (@typeInfo(Child).@"struct".fields) |field|
            setValue(&@field(target.*, field.name), value),
        else => @compileError("test accumulator must contain only floats or structs"),
    }
}

fn allocateFilled(
    comptime T: type,
    allocator: std.mem.Allocator,
    length: usize,
    value: f64,
) ![]T {
    const result = try allocator.alloc(T, length);
    for (result) |*entry| setValue(entry, value);
    return result;
}

fn allocateState(
    allocator: std.mem.Allocator,
    dimensions: Dimensions,
    value: f64,
) !State {
    const lengths = try calculateLengths(dimensions);
    return .{
        .surface_water_heat_by_cell = try allocateFilled(
            SurfaceWaterHeatFlux,
            allocator,
            lengths.cell,
            value,
        ),
        .runoff_organic_by_cell_class = try allocateFilled(
            RunoffOrganicFlux,
            allocator,
            lengths.runoff_organic,
            value,
        ),
        .runoff_element_gas_by_cell = try allocateFilled(
            RunoffElementGasFlux,
            allocator,
            lengths.cell,
            value,
        ),
        .snow_layer_element_gas_by_cell_layer = try allocateFilled(
            SnowLayerElementGasFlux,
            allocator,
            lengths.snow_layer,
            value,
        ),
        .runoff_salt_mol_by_cell_species = try allocateFilled(
            f64,
            allocator,
            lengths.runoff_salt,
            value,
        ),
        .snow_drift_salt_mol_by_cell_species = try allocateFilled(
            f64,
            allocator,
            lengths.snow_drift_salt,
            value,
        ),
        .snow_layer_salt_mol_by_cell_layer_species = try allocateFilled(
            f64,
            allocator,
            lengths.snow_layer_salt,
            value,
        ),
        .erosion_scalar_by_cell = try allocateFilled(
            ErosionScalarFlux,
            allocator,
            lengths.cell,
            value,
        ),
        .erosion_microbial_by_cell_class_group_component = try allocateFilled(
            OrganicElementFlux,
            allocator,
            lengths.microbial,
            value,
        ),
        .erosion_residue_by_cell_class_component = try allocateFilled(
            OrganicElementFlux,
            allocator,
            lengths.residue,
            value,
        ),
        .erosion_adsorbed_by_cell_class = try allocateFilled(
            AdsorbedOrganicFlux,
            allocator,
            lengths.adsorbed,
            value,
        ),
        .erosion_soil_organic_by_cell_class_component = try allocateFilled(
            SoilOrganicFlux,
            allocator,
            lengths.soil_organic,
            value,
        ),
        .snowpack_water_heat_by_cell_layer = try allocateFilled(
            SnowpackWaterHeatFlux,
            allocator,
            lengths.snow_layer,
            value,
        ),
    };
}

fn expectSliceValue(values: anytype, expected: f64) !void {
    for (values) |value| try expectValue(value, expected);
}

fn expectValue(value: anytype, expected: f64) !void {
    switch (@typeInfo(@TypeOf(value))) {
        .float => try std.testing.expectEqual(expected, value),
        .@"struct" => inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
            try expectValue(@field(value, field.name), expected),
        else => @compileError("test accumulator must contain only floats or structs"),
    }
}

fn floatFieldCount(comptime T: type) usize {
    return switch (@typeInfo(T)) {
        .float => 1,
        .@"struct" => |structure| count: {
            var total: usize = 0;
            inline for (structure.fields) |field| total += floatFieldCount(field.type);
            break :count total;
        },
        else => @compileError("accumulator must contain only floats or structs"),
    };
}

const nonlegacy_dimensions = Dimensions{
    .grid_cell_count = 2,
    .snow_layer_count = 3,
    .runoff_organic_class_count = 4,
    .microbial_class_count = 2,
    .microbial_group_count = 3,
    .microbial_component_count = 4,
    .residue_and_soil_organic_class_count = 3,
    .residue_component_count = 4,
    .soil_organic_component_count = 6,
};

test "dynamic salt and erosion reset every runtime-sized inventory" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 7);

    try reset(
        nonlegacy_dimensions,
        .dynamic,
        .freeze_thaw_erosion_and_organic_matter_change,
        &state,
    );

    inline for (@typeInfo(State).@"struct".fields) |field| {
        try expectSliceValue(@field(state, field.name), 0);
    }
}

test "static salt and no-erosion modes preserve their gated inventories" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 5);

    try reset(nonlegacy_dimensions, .static, .freeze_thaw, &state);

    try expectSliceValue(state.surface_water_heat_by_cell, 0);
    try expectSliceValue(state.runoff_organic_by_cell_class, 0);
    try expectSliceValue(state.runoff_element_gas_by_cell, 0);
    try expectSliceValue(state.snow_layer_element_gas_by_cell_layer, 0);
    try expectSliceValue(state.snowpack_water_heat_by_cell_layer, 0);
    try expectSliceValue(state.runoff_salt_mol_by_cell_species, 5);
    try expectSliceValue(state.snow_drift_salt_mol_by_cell_species, 5);
    try expectSliceValue(state.snow_layer_salt_mol_by_cell_layer_species, 5);
    try expectSliceValue(state.erosion_scalar_by_cell, 5);
    try expectSliceValue(state.erosion_microbial_by_cell_class_group_component, 5);
    try expectSliceValue(state.erosion_residue_by_cell_class_component, 5);
    try expectSliceValue(state.erosion_adsorbed_by_cell_class, 5);
    try expectSliceValue(state.erosion_soil_organic_by_cell_class_component, 5);
}

test "source species order and default loop extents are represented exactly" {
    try std.testing.expectEqual(@as(usize, 6), floatFieldCount(SurfaceWaterHeatFlux));
    try std.testing.expectEqual(@as(usize, 22), floatFieldCount(RunoffElementGasFlux));
    try std.testing.expectEqual(@as(usize, 10), floatFieldCount(SnowLayerElementGasFlux));
    try std.testing.expectEqual(@as(usize, 62), floatFieldCount(ErosionScalarFlux));
    try std.testing.expectEqual(@as(usize, 5), floatFieldCount(SnowpackWaterHeatFlux));
    try std.testing.expectEqual(
        @as(usize, 42),
        @typeInfo(RunoffSaltSpecies).@"enum".fields.len,
    );
    try std.testing.expectEqual(
        @as(usize, 41),
        @typeInfo(SnowSaltSpecies).@"enum".fields.len,
    );
    try std.testing.expectEqual(@as(usize, 33), @intFromEnum(RunoffSaltSpecies.hydrogen_sulfate));

    const source_dimensions = Dimensions{
        .grid_cell_count = 1,
        .snow_layer_count = 10,
        .runoff_organic_class_count = 3,
        .microbial_class_count = 6,
        .microbial_group_count = 7,
        .microbial_component_count = 3,
        .residue_and_soil_organic_class_count = 5,
        .residue_component_count = 2,
        .soil_organic_component_count = 5,
    };
    const lengths = try calculateLengths(source_dimensions);
    try std.testing.expectEqual(@as(usize, 3), lengths.runoff_organic);
    try std.testing.expectEqual(@as(usize, 126), lengths.microbial);
    try std.testing.expectEqual(@as(usize, 10), lengths.residue);
    try std.testing.expectEqual(@as(usize, 5), lengths.adsorbed);
    try std.testing.expectEqual(@as(usize, 25), lengths.soil_organic);
    try std.testing.expectEqual(@as(usize, 410), lengths.snow_layer_salt);
}

test "late non-finite accumulator fails before any reset mutation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 3);
    state.snowpack_water_heat_by_cell_layer[
        state.snowpack_water_heat_by_cell_layer.len - 1
    ].convective_heat_megajoules_per_step = std.math.nan(f64);

    try std.testing.expectError(
        error.NonFiniteTransportAccumulator,
        reset(nonlegacy_dimensions, .dynamic, .freeze_thaw_and_erosion, &state),
    );
    try expectSliceValue(state.surface_water_heat_by_cell, 3);
    try std.testing.expect(std.math.isNan(
        state.snowpack_water_heat_by_cell_layer[
            state.snowpack_water_heat_by_cell_layer.len - 1
        ].convective_heat_megajoules_per_step,
    ));
}

test "short storage and overflowing dimensions fail atomically" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 2);
    state.snowpack_water_heat_by_cell_layer =
        state.snowpack_water_heat_by_cell_layer[0 .. state.snowpack_water_heat_by_cell_layer.len - 1];
    try std.testing.expectError(
        error.AccumulatorDimensionMismatch,
        reset(nonlegacy_dimensions, .dynamic, .freeze_thaw_and_erosion, &state),
    );
    try expectSliceValue(state.surface_water_heat_by_cell, 2);

    var overflowing = nonlegacy_dimensions;
    overflowing.grid_cell_count = std.math.maxInt(usize);
    try std.testing.expectError(
        error.AccumulatorDimensionOverflow,
        reset(overflowing, .dynamic, .freeze_thaw_and_erosion, &state),
    );
    try expectSliceValue(state.surface_water_heat_by_cell, 2);
}
