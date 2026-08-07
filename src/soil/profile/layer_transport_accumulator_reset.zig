const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const Dimensions = struct {
    grid_cell_count: usize,
    soil_layer_capacity: usize,
    organic_class_count: usize,
};

pub const ActiveLayerRange = struct {
    first_layer: usize,
    layer_count: usize,
};

pub const WaterHeatFlux = struct {
    micropore_water_m3_per_step: f64 = 0,
    water_vapor_m3_per_step: f64 = 0,
    micropore_wetting_front_water_m3_per_step: f64 = 0,
    macropore_water_m3_per_step: f64 = 0,
    conductive_convective_heat_megajoules_per_step: f64 = 0,
    micropore_freeze_thaw_water_m3_per_step: f64 = 0,
    macropore_freeze_thaw_water_m3_per_step: f64 = 0,
    freeze_thaw_latent_heat_megajoules_per_step: f64 = 0,
    evaporation_condensation_water_m3_per_step: f64 = 0,
    evaporation_condensation_latent_heat_megajoules_per_step: f64 = 0,
};

pub const OrganicPoreFlux = struct {
    micropore_carbon_g_c_per_step: f64 = 0,
    micropore_nitrogen_g_n_per_step: f64 = 0,
    micropore_phosphorus_g_p_per_step: f64 = 0,
    micropore_acetate_carbon_g_c_per_step: f64 = 0,
    macropore_carbon_g_c_per_step: f64 = 0,
    macropore_nitrogen_g_n_per_step: f64 = 0,
    macropore_phosphorus_g_p_per_step: f64 = 0,
    macropore_acetate_carbon_g_c_per_step: f64 = 0,
};

pub const PoreElementFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    hydrogen_g_h_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    nitrite_g_n_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
};

pub const BandElementFlux = struct {
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    nitrite_g_n_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
};

pub const GaseousFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    hydrogen_g_h_per_step: f64 = 0,
};

pub const ElementGasFlux = struct {
    micropore: PoreElementFlux = .{},
    micropore_band: BandElementFlux = .{},
    macropore: PoreElementFlux = .{},
    macropore_band: BandElementFlux = .{},
    gaseous: GaseousFlux = .{},
};

/// Source order of redist.f lines 1827--1868. Values are mol/step.
pub const MicroporeSaltSpecies = enum {
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

/// Source order of redist.f lines 1877--1917. Values are mol/step.
pub const MacroporeSaltSpecies = enum {
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

/// Source order of both fertilizer-band phosphate blocks. Values are mol/step.
pub const BandPhosphateSpecies = enum {
    phosphate,
    phosphoric_acid,
    iron_hpo4,
    iron_h2po4,
    calcium_po4,
    calcium_hpo4,
    calcium_h2po4,
    magnesium_hpo4,
};

/// All slices are caller-owned, runtime allocated, and cell-layer-major.
pub const State = struct {
    water_heat_by_cell_layer: []WaterHeatFlux,
    organic_by_cell_layer_class: []OrganicPoreFlux,
    element_gas_by_cell_layer: []ElementGasFlux,
    micropore_salt_mol_per_step_by_cell_layer_species: []f64,
    micropore_band_phosphate_mol_per_step_by_cell_layer_species: []f64,
    macropore_salt_mol_per_step_by_cell_layer_species: []f64,
    macropore_band_phosphate_mol_per_step_by_cell_layer_species: []f64,
};

const Lengths = struct {
    cell_layer: usize,
    organic: usize,
    micropore_salt: usize,
    band_phosphate: usize,
    macropore_salt: usize,
};

/// Resets active-layer within-soil transport accumulators.
///
/// Traceability: REDIST.F lines 1760--1926. Assignment and group order are
/// water/heat, runtime organic classes, element/gas, then dynamic salt.
/// Static salt equilibrium leaves all salt arrays untouched. Validation of
/// every active target completes before the first mutation.
pub fn reset(
    dimensions: Dimensions,
    active_layers_by_cell: []const ActiveLayerRange,
    salt_equilibrium_mode: SaltEquilibriumMode,
    state: *State,
) !void {
    const lengths = try calculateLengths(dimensions);
    try validateLengths(
        dimensions,
        active_layers_by_cell,
        state.*,
        lengths,
    );
    try preflightActiveValues(
        dimensions,
        active_layers_by_cell,
        salt_equilibrium_mode,
        state.*,
    );

    const micropore_salt_count = speciesCount(MicroporeSaltSpecies);
    const macropore_salt_count = speciesCount(MacroporeSaltSpecies);
    const band_phosphate_count = speciesCount(BandPhosphateSpecies);
    for (active_layers_by_cell, 0..) |active, cell| {
        const cell_start = cell * dimensions.soil_layer_capacity;
        for (active.first_layer..active.first_layer + active.layer_count) |local_layer| {
            const layer = cell_start + local_layer;
            state.water_heat_by_cell_layer[layer] = .{};
            const organic_start = layer * dimensions.organic_class_count;
            for (0..dimensions.organic_class_count) |class| {
                state.organic_by_cell_layer_class[organic_start + class] = .{};
            }
            state.element_gas_by_cell_layer[layer] = .{};
            if (salt_equilibrium_mode == .dynamic) {
                zeroSpecies(
                    state.micropore_salt_mol_per_step_by_cell_layer_species,
                    layer,
                    micropore_salt_count,
                );
                zeroSpecies(
                    state.micropore_band_phosphate_mol_per_step_by_cell_layer_species,
                    layer,
                    band_phosphate_count,
                );
                zeroSpecies(
                    state.macropore_salt_mol_per_step_by_cell_layer_species,
                    layer,
                    macropore_salt_count,
                );
                zeroSpecies(
                    state.macropore_band_phosphate_mol_per_step_by_cell_layer_species,
                    layer,
                    band_phosphate_count,
                );
            }
        }
    }
}

fn speciesCount(comptime Species: type) usize {
    return @typeInfo(Species).@"enum".fields.len;
}

fn calculateLengths(dimensions: Dimensions) !Lengths {
    if (dimensions.grid_cell_count == 0 or
        dimensions.soil_layer_capacity == 0 or
        dimensions.organic_class_count == 0)
    {
        return error.InvalidSoilTransportAccumulatorDimensions;
    }
    const cell_layer = try multiply(
        dimensions.grid_cell_count,
        dimensions.soil_layer_capacity,
    );
    return .{
        .cell_layer = cell_layer,
        .organic = try multiply(cell_layer, dimensions.organic_class_count),
        .micropore_salt = try multiply(
            cell_layer,
            speciesCount(MicroporeSaltSpecies),
        ),
        .band_phosphate = try multiply(
            cell_layer,
            speciesCount(BandPhosphateSpecies),
        ),
        .macropore_salt = try multiply(
            cell_layer,
            speciesCount(MacroporeSaltSpecies),
        ),
    };
}

fn multiply(left: usize, right: usize) !usize {
    return std.math.mul(usize, left, right) catch
        return error.SoilTransportAccumulatorDimensionOverflow;
}

fn validateLengths(
    dimensions: Dimensions,
    active_layers_by_cell: []const ActiveLayerRange,
    state: State,
    lengths: Lengths,
) !void {
    if (active_layers_by_cell.len != dimensions.grid_cell_count or
        state.water_heat_by_cell_layer.len != lengths.cell_layer or
        state.organic_by_cell_layer_class.len != lengths.organic or
        state.element_gas_by_cell_layer.len != lengths.cell_layer or
        state.micropore_salt_mol_per_step_by_cell_layer_species.len !=
            lengths.micropore_salt or
        state.micropore_band_phosphate_mol_per_step_by_cell_layer_species.len !=
            lengths.band_phosphate or
        state.macropore_salt_mol_per_step_by_cell_layer_species.len !=
            lengths.macropore_salt or
        state.macropore_band_phosphate_mol_per_step_by_cell_layer_species.len !=
            lengths.band_phosphate)
    {
        return error.SoilTransportAccumulatorDimensionMismatch;
    }
    for (active_layers_by_cell) |active| {
        if (active.layer_count == 0 or
            active.first_layer > dimensions.soil_layer_capacity or
            active.layer_count > dimensions.soil_layer_capacity - active.first_layer)
        {
            return error.InvalidActiveSoilLayerRange;
        }
    }
}

fn preflightActiveValues(
    dimensions: Dimensions,
    active_layers_by_cell: []const ActiveLayerRange,
    salt_equilibrium_mode: SaltEquilibriumMode,
    state: State,
) !void {
    const micropore_salt_count = speciesCount(MicroporeSaltSpecies);
    const macropore_salt_count = speciesCount(MacroporeSaltSpecies);
    const band_phosphate_count = speciesCount(BandPhosphateSpecies);
    for (active_layers_by_cell, 0..) |active, cell| {
        const cell_start = cell * dimensions.soil_layer_capacity;
        for (active.first_layer..active.first_layer + active.layer_count) |local_layer| {
            const layer = cell_start + local_layer;
            try validateFinite(state.water_heat_by_cell_layer[layer]);
            const organic_start = layer * dimensions.organic_class_count;
            for (0..dimensions.organic_class_count) |class| {
                try validateFinite(state.organic_by_cell_layer_class[organic_start + class]);
            }
            try validateFinite(state.element_gas_by_cell_layer[layer]);
            if (salt_equilibrium_mode == .dynamic) {
                try validateSpecies(
                    state.micropore_salt_mol_per_step_by_cell_layer_species,
                    layer,
                    micropore_salt_count,
                );
                try validateSpecies(
                    state.micropore_band_phosphate_mol_per_step_by_cell_layer_species,
                    layer,
                    band_phosphate_count,
                );
                try validateSpecies(
                    state.macropore_salt_mol_per_step_by_cell_layer_species,
                    layer,
                    macropore_salt_count,
                );
                try validateSpecies(
                    state.macropore_band_phosphate_mol_per_step_by_cell_layer_species,
                    layer,
                    band_phosphate_count,
                );
            }
        }
    }
}

fn validateSpecies(values: []const f64, layer: usize, count: usize) !void {
    const start = layer * count;
    for (values[start .. start + count]) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSoilTransportAccumulator;
    }
}

fn zeroSpecies(values: []f64, layer: usize, count: usize) void {
    const start = layer * count;
    @memset(values[start .. start + count], 0);
}

fn validateFinite(value: anytype) !void {
    switch (@typeInfo(@TypeOf(value))) {
        .float => if (!std.math.isFinite(value))
            return error.NonFiniteSoilTransportAccumulator,
        .@"struct" => inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
            try validateFinite(@field(value, field.name)),
        else => @compileError("soil transport accumulator must contain only floats or structs"),
    }
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
    const values = try allocator.alloc(T, length);
    for (values) |*entry| setValue(entry, value);
    return values;
}

fn allocateState(
    allocator: std.mem.Allocator,
    dimensions: Dimensions,
    value: f64,
) !State {
    const lengths = try calculateLengths(dimensions);
    return .{
        .water_heat_by_cell_layer = try allocateFilled(
            WaterHeatFlux,
            allocator,
            lengths.cell_layer,
            value,
        ),
        .organic_by_cell_layer_class = try allocateFilled(
            OrganicPoreFlux,
            allocator,
            lengths.organic,
            value,
        ),
        .element_gas_by_cell_layer = try allocateFilled(
            ElementGasFlux,
            allocator,
            lengths.cell_layer,
            value,
        ),
        .micropore_salt_mol_per_step_by_cell_layer_species = try allocateFilled(
            f64,
            allocator,
            lengths.micropore_salt,
            value,
        ),
        .micropore_band_phosphate_mol_per_step_by_cell_layer_species = try allocateFilled(
            f64,
            allocator,
            lengths.band_phosphate,
            value,
        ),
        .macropore_salt_mol_per_step_by_cell_layer_species = try allocateFilled(
            f64,
            allocator,
            lengths.macropore_salt,
            value,
        ),
        .macropore_band_phosphate_mol_per_step_by_cell_layer_species = try allocateFilled(
            f64,
            allocator,
            lengths.band_phosphate,
            value,
        ),
    };
}

fn expectValue(actual: anytype, expected: f64) !void {
    switch (@typeInfo(@TypeOf(actual))) {
        .float => try std.testing.expectEqual(expected, actual),
        .@"struct" => inline for (@typeInfo(@TypeOf(actual)).@"struct".fields) |field|
            try expectValue(@field(actual, field.name), expected),
        else => @compileError("test accumulator must contain only floats or structs"),
    }
}

fn expectLayerValue(
    state: State,
    dimensions: Dimensions,
    layer: usize,
    expected: f64,
) !void {
    try expectValue(state.water_heat_by_cell_layer[layer], expected);
    try expectValue(state.element_gas_by_cell_layer[layer], expected);
    const organic_start = layer * dimensions.organic_class_count;
    for (state.organic_by_cell_layer_class[organic_start .. organic_start + dimensions.organic_class_count]) |value| try expectValue(value, expected);
    try expectSaltLayerValue(state, layer, expected);
}

fn expectSaltLayerValue(state: State, layer: usize, expected: f64) !void {
    inline for (.{
        .{ state.micropore_salt_mol_per_step_by_cell_layer_species, speciesCount(MicroporeSaltSpecies) },
        .{ state.micropore_band_phosphate_mol_per_step_by_cell_layer_species, speciesCount(BandPhosphateSpecies) },
        .{ state.macropore_salt_mol_per_step_by_cell_layer_species, speciesCount(MacroporeSaltSpecies) },
        .{ state.macropore_band_phosphate_mol_per_step_by_cell_layer_species, speciesCount(BandPhosphateSpecies) },
    }) |entry| {
        const start = layer * entry[1];
        for (entry[0][start .. start + entry[1]]) |value|
            try std.testing.expectEqual(expected, value);
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
    .soil_layer_capacity = 4,
    .organic_class_count = 3,
};

const nonlegacy_active = [_]ActiveLayerRange{
    .{ .first_layer = 1, .layer_count = 2 },
    .{ .first_layer = 0, .layer_count = 4 },
};

test "dynamic reset covers runtime active layers and preserves inactive inventories" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 7);

    try reset(nonlegacy_dimensions, &nonlegacy_active, .dynamic, &state);

    for (0..nonlegacy_dimensions.grid_cell_count) |cell| {
        const active = nonlegacy_active[cell];
        for (0..nonlegacy_dimensions.soil_layer_capacity) |local_layer| {
            const layer = cell * nonlegacy_dimensions.soil_layer_capacity + local_layer;
            const is_active = local_layer >= active.first_layer and
                local_layer < active.first_layer + active.layer_count;
            try expectLayerValue(
                state,
                nonlegacy_dimensions,
                layer,
                if (is_active) 0 else 7,
            );
        }
    }
}

test "static salt mode preserves active salt while resetting other groups" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 5);

    try reset(nonlegacy_dimensions, &nonlegacy_active, .static, &state);

    const active_layer = 1;
    try expectValue(state.water_heat_by_cell_layer[active_layer], 0);
    try expectValue(state.element_gas_by_cell_layer[active_layer], 0);
    for (state.organic_by_cell_layer_class[3..6]) |value|
        try expectValue(value, 0);
    try expectSaltLayerValue(state, active_layer, 5);
}

test "source scalar counts and default organic extent remain exact" {
    try std.testing.expectEqual(@as(usize, 10), floatFieldCount(WaterHeatFlux));
    try std.testing.expectEqual(@as(usize, 8), floatFieldCount(OrganicPoreFlux));
    try std.testing.expectEqual(@as(usize, 43), floatFieldCount(ElementGasFlux));
    try std.testing.expectEqual(@as(usize, 42), speciesCount(MicroporeSaltSpecies));
    try std.testing.expectEqual(@as(usize, 41), speciesCount(MacroporeSaltSpecies));
    try std.testing.expectEqual(@as(usize, 8), speciesCount(BandPhosphateSpecies));
    const source = try calculateLengths(.{
        .grid_cell_count = 1,
        .soil_layer_capacity = 10,
        .organic_class_count = 5,
    });
    try std.testing.expectEqual(@as(usize, 50), source.organic);
    try std.testing.expectEqual(@as(usize, 420), source.micropore_salt);
    try std.testing.expectEqual(@as(usize, 410), source.macropore_salt);
    try std.testing.expectEqual(@as(usize, 80), source.band_phosphate);
}

test "late non-finite active salt fails before any mutation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 3);
    const final_active_layer = 7;
    const final_salt =
        (final_active_layer + 1) * speciesCount(MacroporeSaltSpecies) - 1;
    state.macropore_salt_mol_per_step_by_cell_layer_species[final_salt] =
        std.math.nan(f64);

    try std.testing.expectError(
        error.NonFiniteSoilTransportAccumulator,
        reset(nonlegacy_dimensions, &nonlegacy_active, .dynamic, &state),
    );
    try expectLayerValue(state, nonlegacy_dimensions, 1, 3);
    try std.testing.expect(std.math.isNan(
        state.macropore_salt_mol_per_step_by_cell_layer_species[final_salt],
    ));
}

test "short storage and overflowing runtime dimensions fail atomically" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var state = try allocateState(arena.allocator(), nonlegacy_dimensions, 2);
    state.element_gas_by_cell_layer =
        state.element_gas_by_cell_layer[0 .. state.element_gas_by_cell_layer.len - 1];
    try std.testing.expectError(
        error.SoilTransportAccumulatorDimensionMismatch,
        reset(nonlegacy_dimensions, &nonlegacy_active, .dynamic, &state),
    );
    try expectValue(state.water_heat_by_cell_layer[1], 2);

    var overflowing = nonlegacy_dimensions;
    overflowing.grid_cell_count = std.math.maxInt(usize);
    try std.testing.expectError(
        error.SoilTransportAccumulatorDimensionOverflow,
        reset(overflowing, &nonlegacy_active, .dynamic, &state),
    );
    try expectValue(state.water_heat_by_cell_layer[1], 2);
}
