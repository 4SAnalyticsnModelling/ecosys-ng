const std = @import("std");

pub const Dimensions = struct {
    grid_cell_count: usize,
    soil_layer_capacity: usize,
};

/// Cell-local, zero-based active soil-layer window.
pub const ActiveLayerRange = struct {
    first_layer: usize,
    layer_count: usize,
};

/// Tracked-element masses corresponding to REDIST's seven gaseous species.
pub const GaseousMasses = struct {
    carbon_dioxide_g_c: []const f64,
    methane_g_c: []const f64,
    oxygen_g_o: []const f64,
    dinitrogen_g_n: []const f64,
    nitrous_oxide_g_n: []const f64,
    ammonia_g_n: []const f64,
    hydrogen_g_h: []const f64,
};

/// Runtime values reproduce REDIST's divisors and ideal-gas coefficient when
/// configured as 12, 12, 32, 28, 28, 14, 2, and 1.2194e4 respectively.
pub const RuntimeParameters = struct {
    carbon_dioxide_carbon_g_c_per_mol: f64,
    methane_carbon_g_c_per_mol: f64,
    oxygen_g_o_per_mol: f64,
    dinitrogen_g_n_per_mol: f64,
    nitrous_oxide_g_n_per_mol: f64,
    ammonia_g_n_per_mol: f64,
    hydrogen_g_h_per_mol: f64,
    atmospheric_pressure_over_gas_constant_mol_k_per_m3: f64,
    minimum_air_filled_porosity_m3_per_m3: f64,
};

pub const Inputs = struct {
    dimensions: Dimensions,
    active_layers_by_cell: []const ActiveLayerRange,
    air_filled_porosity_m3_per_m3: []const f64,
    soil_air_volume_m3: []const f64,
    soil_temperature_k: []const f64,
    gaseous_mass: GaseousMasses,
    matrix_liquid_water_m3: []const f64,
    matrix_ice_water_m3: []const f64,
    macropore_liquid_water_m3: []const f64,
    macropore_ice_water_m3: []const f64,
    parameters: RuntimeParameters,
};

pub const Outputs = struct {
    /// Null replaces REDIST's ambiguous `LG=0` sentinel. A non-null value is
    /// the zero-based layer index local to its grid cell.
    lowest_bubble_release_layer_by_cell: []?usize,
    matrix_liquid_water_m3_before_transport: []f64,
    matrix_ice_water_m3_before_transport: []f64,
    macropore_liquid_water_m3_before_transport: []f64,
    macropore_ice_water_m3_before_transport: []f64,
};

const LayerAssessment = struct {
    total_gas_mol: f64,
    atmospheric_capacity_mol: f64,
};

/// Selects the bubble-release layer and snapshots soil water inventories.
///
/// Traceability: REDIST.F lines 1722--1756 (`LG`, `LX`, `V*G2`, `VTATM`,
/// `VTGAS`, and `VOL*1`). The source's strict `<` porosity and `>` gas
/// capacity barriers are retained. Once encountered, a barrier remains
/// latched for all deeper layers in that cell. All inputs are validated before
/// output mutation, and inactive layer slots remain untouched.
pub fn selectAndSnapshot(inputs: Inputs, outputs: *Outputs) !void {
    try validateDimensions(inputs, outputs.*);
    try validateParameters(inputs.parameters);
    try preflightActiveLayers(inputs);

    for (inputs.active_layers_by_cell, 0..) |active, cell| {
        const cell_start = cell * inputs.dimensions.soil_layer_capacity;
        var barrier_encountered = false;
        var lowest_layer: ?usize = null;
        for (active.first_layer..active.first_layer + active.layer_count) |local_layer| {
            const layer = cell_start + local_layer;
            const assessment = assessLayer(inputs, layer) catch unreachable;
            const porosity = inputs.air_filled_porosity_m3_per_m3[layer];
            if (porosity < inputs.parameters.minimum_air_filled_porosity_m3_per_m3 or
                assessment.total_gas_mol > assessment.atmospheric_capacity_mol)
            {
                barrier_encountered = true;
            }
            if (porosity >= inputs.parameters.minimum_air_filled_porosity_m3_per_m3 and
                !barrier_encountered)
            {
                lowest_layer = local_layer;
            }

            outputs.matrix_liquid_water_m3_before_transport[layer] =
                inputs.matrix_liquid_water_m3[layer];
            outputs.matrix_ice_water_m3_before_transport[layer] =
                inputs.matrix_ice_water_m3[layer];
            outputs.macropore_liquid_water_m3_before_transport[layer] =
                inputs.macropore_liquid_water_m3[layer];
            outputs.macropore_ice_water_m3_before_transport[layer] =
                inputs.macropore_ice_water_m3[layer];
        }
        outputs.lowest_bubble_release_layer_by_cell[cell] = lowest_layer;
    }
}

fn validateDimensions(inputs: Inputs, outputs: Outputs) !void {
    const dimensions = inputs.dimensions;
    if (dimensions.grid_cell_count == 0 or dimensions.soil_layer_capacity == 0)
        return error.InvalidBubbleLayerDimensions;
    const total = std.math.mul(
        usize,
        dimensions.grid_cell_count,
        dimensions.soil_layer_capacity,
    ) catch return error.BubbleLayerDimensionOverflow;
    if (inputs.active_layers_by_cell.len != dimensions.grid_cell_count or
        outputs.lowest_bubble_release_layer_by_cell.len != dimensions.grid_cell_count)
    {
        return error.BubbleLayerDimensionMismatch;
    }
    inline for (.{
        inputs.air_filled_porosity_m3_per_m3,
        inputs.soil_air_volume_m3,
        inputs.soil_temperature_k,
        inputs.gaseous_mass.carbon_dioxide_g_c,
        inputs.gaseous_mass.methane_g_c,
        inputs.gaseous_mass.oxygen_g_o,
        inputs.gaseous_mass.dinitrogen_g_n,
        inputs.gaseous_mass.nitrous_oxide_g_n,
        inputs.gaseous_mass.ammonia_g_n,
        inputs.gaseous_mass.hydrogen_g_h,
        inputs.matrix_liquid_water_m3,
        inputs.matrix_ice_water_m3,
        inputs.macropore_liquid_water_m3,
        inputs.macropore_ice_water_m3,
    }) |values| if (values.len != total)
        return error.BubbleLayerDimensionMismatch;
    inline for (.{
        outputs.matrix_liquid_water_m3_before_transport,
        outputs.matrix_ice_water_m3_before_transport,
        outputs.macropore_liquid_water_m3_before_transport,
        outputs.macropore_ice_water_m3_before_transport,
    }) |values| if (values.len != total)
        return error.BubbleLayerDimensionMismatch;
}

fn validateParameters(parameters: RuntimeParameters) !void {
    inline for (@typeInfo(RuntimeParameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteBubbleLayerParameter;
    }
    inline for (.{
        parameters.carbon_dioxide_carbon_g_c_per_mol,
        parameters.methane_carbon_g_c_per_mol,
        parameters.oxygen_g_o_per_mol,
        parameters.dinitrogen_g_n_per_mol,
        parameters.nitrous_oxide_g_n_per_mol,
        parameters.ammonia_g_n_per_mol,
        parameters.hydrogen_g_h_per_mol,
        parameters.atmospheric_pressure_over_gas_constant_mol_k_per_m3,
    }) |value| if (value <= 0)
        return error.InvalidBubbleLayerParameter;
    if (parameters.minimum_air_filled_porosity_m3_per_m3 < 0 or
        parameters.minimum_air_filled_porosity_m3_per_m3 > 1)
    {
        return error.InvalidBubbleLayerParameter;
    }
}

fn preflightActiveLayers(inputs: Inputs) !void {
    for (inputs.active_layers_by_cell, 0..) |active, cell| {
        if (active.layer_count == 0 or
            active.first_layer > inputs.dimensions.soil_layer_capacity or
            active.layer_count >
                inputs.dimensions.soil_layer_capacity - active.first_layer)
        {
            return error.InvalidActiveSoilLayerRange;
        }
        const cell_start = cell * inputs.dimensions.soil_layer_capacity;
        for (active.first_layer..active.first_layer + active.layer_count) |local_layer| {
            const layer = cell_start + local_layer;
            const porosity = inputs.air_filled_porosity_m3_per_m3[layer];
            if (!std.math.isFinite(porosity) or porosity < 0 or porosity > 1)
                return error.InvalidBubbleLayerState;
            const air_volume = inputs.soil_air_volume_m3[layer];
            const temperature = inputs.soil_temperature_k[layer];
            if (!std.math.isFinite(air_volume) or air_volume < 0 or
                !std.math.isFinite(temperature) or temperature <= 0)
            {
                return error.InvalidBubbleLayerState;
            }
            inline for (.{
                inputs.gaseous_mass.carbon_dioxide_g_c[layer],
                inputs.gaseous_mass.methane_g_c[layer],
                inputs.gaseous_mass.oxygen_g_o[layer],
                inputs.gaseous_mass.dinitrogen_g_n[layer],
                inputs.gaseous_mass.nitrous_oxide_g_n[layer],
                inputs.gaseous_mass.ammonia_g_n[layer],
                inputs.gaseous_mass.hydrogen_g_h[layer],
                inputs.matrix_liquid_water_m3[layer],
                inputs.matrix_ice_water_m3[layer],
                inputs.macropore_liquid_water_m3[layer],
                inputs.macropore_ice_water_m3[layer],
            }) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidBubbleLayerState;
            _ = try assessLayer(inputs, layer);
        }
    }
}

fn assessLayer(inputs: Inputs, layer: usize) !LayerAssessment {
    const parameters = inputs.parameters;
    const carbon_dioxide_mol =
        inputs.gaseous_mass.carbon_dioxide_g_c[layer] /
        parameters.carbon_dioxide_carbon_g_c_per_mol;
    const methane_mol =
        inputs.gaseous_mass.methane_g_c[layer] /
        parameters.methane_carbon_g_c_per_mol;
    const oxygen_mol =
        inputs.gaseous_mass.oxygen_g_o[layer] / parameters.oxygen_g_o_per_mol;
    const dinitrogen_mol =
        inputs.gaseous_mass.dinitrogen_g_n[layer] / parameters.dinitrogen_g_n_per_mol;
    const nitrous_oxide_mol =
        inputs.gaseous_mass.nitrous_oxide_g_n[layer] /
        parameters.nitrous_oxide_g_n_per_mol;
    const ammonia_mol =
        inputs.gaseous_mass.ammonia_g_n[layer] / parameters.ammonia_g_n_per_mol;
    const hydrogen_mol =
        inputs.gaseous_mass.hydrogen_g_h[layer] / parameters.hydrogen_g_h_per_mol;

    var total_gas_mol = carbon_dioxide_mol;
    total_gas_mol += methane_mol;
    total_gas_mol += oxygen_mol;
    total_gas_mol += dinitrogen_mol;
    total_gas_mol += nitrous_oxide_mol;
    total_gas_mol += ammonia_mol;
    total_gas_mol += hydrogen_mol;
    const scaled_air_volume =
        parameters.atmospheric_pressure_over_gas_constant_mol_k_per_m3 *
        inputs.soil_air_volume_m3[layer];
    const atmospheric_capacity_mol = @max(
        0,
        scaled_air_volume / inputs.soil_temperature_k[layer],
    );
    if (!std.math.isFinite(total_gas_mol) or
        !std.math.isFinite(atmospheric_capacity_mol))
    {
        return error.NonFiniteBubbleLayerCalculation;
    }
    return .{
        .total_gas_mol = total_gas_mol,
        .atmospheric_capacity_mol = atmospheric_capacity_mol,
    };
}

const compatibility_parameters = RuntimeParameters{
    .carbon_dioxide_carbon_g_c_per_mol = 12,
    .methane_carbon_g_c_per_mol = 12,
    .oxygen_g_o_per_mol = 32,
    .dinitrogen_g_n_per_mol = 28,
    .nitrous_oxide_g_n_per_mol = 28,
    .ammonia_g_n_per_mol = 14,
    .hydrogen_g_h_per_mol = 2,
    .atmospheric_pressure_over_gas_constant_mol_k_per_m3 = 1.2194e4,
    .minimum_air_filled_porosity_m3_per_m3 = 1.0e-3,
};

fn zeroMasses(values: []const f64) GaseousMasses {
    return .{
        .carbon_dioxide_g_c = values,
        .methane_g_c = values,
        .oxygen_g_o = values,
        .dinitrogen_g_n = values,
        .nitrous_oxide_g_n = values,
        .ammonia_g_n = values,
        .hydrogen_g_h = values,
    };
}

test "deepest continuous eligible layer is selected and active water is snapshotted" {
    const dimensions = Dimensions{ .grid_cell_count = 3, .soil_layer_capacity = 4 };
    const active = [_]ActiveLayerRange{
        .{ .first_layer = 1, .layer_count = 3 },
        .{ .first_layer = 0, .layer_count = 4 },
        .{ .first_layer = 0, .layer_count = 3 },
    };
    var porosity = [_]f64{0.2} ** 12;
    porosity[5] = 0.5e-3;
    const air_volume = [_]f64{1} ** 12;
    const temperature = [_]f64{300} ** 12;
    const zero = [_]f64{0} ** 12;
    var carbon_dioxide = zero;
    carbon_dioxide[8] = 600;
    const water = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    const ice = [_]f64{0.5} ** 12;
    const macropore_water = [_]f64{0.25} ** 12;
    const macropore_ice = [_]f64{0.125} ** 12;
    var selected = [_]?usize{99} ** 3;
    var water_before = [_]f64{-1} ** 12;
    var ice_before = [_]f64{-1} ** 12;
    var macropore_water_before = [_]f64{-1} ** 12;
    var macropore_ice_before = [_]f64{-1} ** 12;
    var gases = zeroMasses(&zero);
    gases.carbon_dioxide_g_c = &carbon_dioxide;
    const inputs = Inputs{
        .dimensions = dimensions,
        .active_layers_by_cell = &active,
        .air_filled_porosity_m3_per_m3 = &porosity,
        .soil_air_volume_m3 = &air_volume,
        .soil_temperature_k = &temperature,
        .gaseous_mass = gases,
        .matrix_liquid_water_m3 = &water,
        .matrix_ice_water_m3 = &ice,
        .macropore_liquid_water_m3 = &macropore_water,
        .macropore_ice_water_m3 = &macropore_ice,
        .parameters = compatibility_parameters,
    };
    var outputs = Outputs{
        .lowest_bubble_release_layer_by_cell = &selected,
        .matrix_liquid_water_m3_before_transport = &water_before,
        .matrix_ice_water_m3_before_transport = &ice_before,
        .macropore_liquid_water_m3_before_transport = &macropore_water_before,
        .macropore_ice_water_m3_before_transport = &macropore_ice_before,
    };

    try selectAndSnapshot(inputs, &outputs);

    try std.testing.expectEqual(@as(?usize, 3), selected[0]);
    try std.testing.expectEqual(@as(?usize, 0), selected[1]);
    try std.testing.expectEqual(@as(?usize, null), selected[2]);
    var source_water_total: f64 = 0;
    var snapshot_water_total: f64 = 0;
    for (active, 0..) |range, cell| {
        const base = cell * dimensions.soil_layer_capacity;
        for (range.first_layer..range.first_layer + range.layer_count) |local| {
            const layer = base + local;
            try std.testing.expectEqual(water[layer], water_before[layer]);
            try std.testing.expectEqual(ice[layer], ice_before[layer]);
            try std.testing.expectEqual(macropore_water[layer], macropore_water_before[layer]);
            try std.testing.expectEqual(macropore_ice[layer], macropore_ice_before[layer]);
            source_water_total += water[layer] + ice[layer] +
                macropore_water[layer] + macropore_ice[layer];
            snapshot_water_total += water_before[layer] + ice_before[layer] +
                macropore_water_before[layer] + macropore_ice_before[layer];
        }
    }
    try std.testing.expectEqual(source_water_total, snapshot_water_total);
    try std.testing.expectEqual(@as(f64, -1), water_before[0]);
    try std.testing.expectEqual(@as(f64, -1), water_before[11]);
}

test "source strict thresholds accept equality" {
    const dimensions = Dimensions{ .grid_cell_count = 1, .soil_layer_capacity = 2 };
    const active = [_]ActiveLayerRange{.{ .first_layer = 0, .layer_count = 2 }};
    const porosity = [_]f64{ 0.1, 0.1 };
    const air_volume = [_]f64{ 1, 1 };
    const temperature = [_]f64{ 300, 300 };
    const carbon_dioxide = [_]f64{ 1, 1 };
    const zero = [_]f64{ 0, 0 };
    var gases = zeroMasses(&zero);
    gases.carbon_dioxide_g_c = &carbon_dioxide;
    const water = [_]f64{ 1, 2 };
    var selected = [_]?usize{null};
    var snapshot = [_]f64{ 9, 9 };
    var snapshot_two = snapshot;
    var snapshot_three = snapshot;
    var snapshot_four = snapshot;
    var parameters = compatibility_parameters;
    parameters.carbon_dioxide_carbon_g_c_per_mol = 1;
    parameters.atmospheric_pressure_over_gas_constant_mol_k_per_m3 = 300;
    parameters.minimum_air_filled_porosity_m3_per_m3 = 0.1;
    const inputs = Inputs{
        .dimensions = dimensions,
        .active_layers_by_cell = &active,
        .air_filled_porosity_m3_per_m3 = &porosity,
        .soil_air_volume_m3 = &air_volume,
        .soil_temperature_k = &temperature,
        .gaseous_mass = gases,
        .matrix_liquid_water_m3 = &water,
        .matrix_ice_water_m3 = &zero,
        .macropore_liquid_water_m3 = &zero,
        .macropore_ice_water_m3 = &zero,
        .parameters = parameters,
    };
    var outputs = Outputs{
        .lowest_bubble_release_layer_by_cell = &selected,
        .matrix_liquid_water_m3_before_transport = &snapshot,
        .matrix_ice_water_m3_before_transport = &snapshot_two,
        .macropore_liquid_water_m3_before_transport = &snapshot_three,
        .macropore_ice_water_m3_before_transport = &snapshot_four,
    };

    try selectAndSnapshot(inputs, &outputs);
    try std.testing.expectEqual(@as(?usize, 1), selected[0]);
}

test "tracked gas mole inventory retains source conversion and addition order" {
    const one = [_]f64{1};
    const inputs = Inputs{
        .dimensions = .{ .grid_cell_count = 1, .soil_layer_capacity = 1 },
        .active_layers_by_cell = &.{.{ .first_layer = 0, .layer_count = 1 }},
        .air_filled_porosity_m3_per_m3 = &one,
        .soil_air_volume_m3 = &one,
        .soil_temperature_k = &.{300},
        .gaseous_mass = .{
            .carbon_dioxide_g_c = &.{12},
            .methane_g_c = &.{24},
            .oxygen_g_o = &.{32},
            .dinitrogen_g_n = &.{56},
            .nitrous_oxide_g_n = &.{84},
            .ammonia_g_n = &.{28},
            .hydrogen_g_h = &.{4},
        },
        .matrix_liquid_water_m3 = &one,
        .matrix_ice_water_m3 = &one,
        .macropore_liquid_water_m3 = &one,
        .macropore_ice_water_m3 = &one,
        .parameters = compatibility_parameters,
    };
    const assessment = try assessLayer(inputs, 0);
    try std.testing.expectEqual(@as(f64, 13), assessment.total_gas_mol);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.2194e4) / 300.0,
        assessment.atmospheric_capacity_mol,
        1e-15,
    );
}

test "late invalid active state fails before outputs mutate" {
    const dimensions = Dimensions{ .grid_cell_count = 1, .soil_layer_capacity = 3 };
    const active = [_]ActiveLayerRange{.{ .first_layer = 0, .layer_count = 3 }};
    const porosity = [_]f64{0.2} ** 3;
    const air_volume = [_]f64{1} ** 3;
    const temperature = [_]f64{300} ** 3;
    const zero = [_]f64{0} ** 3;
    const water = [_]f64{ 1, 2, std.math.nan(f64) };
    var selected = [_]?usize{2};
    var snapshot = [_]f64{9} ** 3;
    var snapshot_two = snapshot;
    var snapshot_three = snapshot;
    var snapshot_four = snapshot;
    const inputs = Inputs{
        .dimensions = dimensions,
        .active_layers_by_cell = &active,
        .air_filled_porosity_m3_per_m3 = &porosity,
        .soil_air_volume_m3 = &air_volume,
        .soil_temperature_k = &temperature,
        .gaseous_mass = zeroMasses(&zero),
        .matrix_liquid_water_m3 = &water,
        .matrix_ice_water_m3 = &zero,
        .macropore_liquid_water_m3 = &zero,
        .macropore_ice_water_m3 = &zero,
        .parameters = compatibility_parameters,
    };
    var outputs = Outputs{
        .lowest_bubble_release_layer_by_cell = &selected,
        .matrix_liquid_water_m3_before_transport = &snapshot,
        .matrix_ice_water_m3_before_transport = &snapshot_two,
        .macropore_liquid_water_m3_before_transport = &snapshot_three,
        .macropore_ice_water_m3_before_transport = &snapshot_four,
    };

    try std.testing.expectError(
        error.InvalidBubbleLayerState,
        selectAndSnapshot(inputs, &outputs),
    );
    try std.testing.expectEqual(@as(?usize, 2), selected[0]);
    for (snapshot) |value| try std.testing.expectEqual(@as(f64, 9), value);
}

test "invalid range and dimension overflow fail explicitly" {
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const active = [_]ActiveLayerRange{.{ .first_layer = 1, .layer_count = 1 }};
    var selected = [_]?usize{null};
    var snapshot = [_]f64{0};
    const inputs = Inputs{
        .dimensions = .{ .grid_cell_count = 1, .soil_layer_capacity = 1 },
        .active_layers_by_cell = &active,
        .air_filled_porosity_m3_per_m3 = &one,
        .soil_air_volume_m3 = &one,
        .soil_temperature_k = &one,
        .gaseous_mass = zeroMasses(&zero),
        .matrix_liquid_water_m3 = &one,
        .matrix_ice_water_m3 = &one,
        .macropore_liquid_water_m3 = &one,
        .macropore_ice_water_m3 = &one,
        .parameters = compatibility_parameters,
    };
    var outputs = Outputs{
        .lowest_bubble_release_layer_by_cell = &selected,
        .matrix_liquid_water_m3_before_transport = &snapshot,
        .matrix_ice_water_m3_before_transport = &snapshot,
        .macropore_liquid_water_m3_before_transport = &snapshot,
        .macropore_ice_water_m3_before_transport = &snapshot,
    };
    try std.testing.expectError(
        error.InvalidActiveSoilLayerRange,
        selectAndSnapshot(inputs, &outputs),
    );

    var overflowing = inputs;
    overflowing.dimensions = .{
        .grid_cell_count = std.math.maxInt(usize),
        .soil_layer_capacity = 2,
    };
    try std.testing.expectError(
        error.BubbleLayerDimensionOverflow,
        selectAndSnapshot(overflowing, &outputs),
    );
}
