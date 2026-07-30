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

pub const WaterHeatFlux = struct {
    micropore_water_m3_per_step: f64 = 0,
    water_vapor_m3_per_step: f64 = 0,
    /// Micropore water including the infiltration wetting-front adjustment.
    wetting_front_micropore_water_m3_per_step: f64 = 0,
    macropore_water_m3_per_step: f64 = 0,
    convective_heat_MJ_per_step: f64 = 0,
};

pub const PositiveNeighborColumn = struct {
    /// Initial mapped layer (`N6`) before the vertical active-layer scan.
    initial_layer_index: usize,
    deepest_layer_index_inclusive: usize,
    /// Vertical soil thickness by neighbor layer, m.
    layer_thickness_m: []const f64,
    /// Axis-specific neighbor flux by layer.
    flux_by_layer: []const WaterHeatFlux,
    /// Stored flux used if a neighboring pond surface disappeared.
    lake_surface_disappearance_flux: WaterHeatFlux,
};

pub const Inputs = struct {
    transport_axis: TransportAxis,
    horizontal_exchange: HorizontalExchange,
    layer_activity_threshold_m: f64,
    current_layer_index: usize,
    current_top_active_layer_index: usize,
    current_layer_thickness_m: f64,
    current_cell_flux: WaterHeatFlux,
    positive_neighbor: PositiveNeighborColumn,
};

pub const State = struct {
    net_flux: WaterHeatFlux,
};

/// Aggregates adjacent-cell micropore, vapor, macropore, and heat fluxes.
///
/// Traceability: REDIST.F lines 3350--3388. Horizontal standalone cells are
/// skipped, while vertical transfer remains enabled (`NCN == 3`, `N == 3`).
/// On the vertical axis, the positive column is scanned from mapped `N6`
/// through `NL` for the first layer thicker than `DLYRM`; if none is found,
/// source behavior retains the initial layer. A current vertical surface layer
/// subtracts the stored pond-disappearance flux; every other active layer
/// subtracts the selected positive-neighbor layer flux. All five source-order
/// updates commit atomically after finite evaluation.
pub fn aggregate(inputs: Inputs, state: *State) !void {
    if (inputs.transport_axis != .vertical and
        inputs.horizontal_exchange == .standalone)
    {
        return;
    }
    try validateInputs(inputs, state.*);
    const positive_layer = selectPositiveLayer(inputs);
    if (inputs.current_layer_thickness_m <= inputs.layer_activity_threshold_m)
        return;

    const outgoing = if (inputs.transport_axis == .vertical and
        inputs.current_layer_index == inputs.current_top_active_layer_index)
        inputs.positive_neighbor.lake_surface_disappearance_flux
    else
        inputs.positive_neighbor.flux_by_layer[positive_layer];

    var candidate = state.net_flux;
    inline for (@typeInfo(WaterHeatFlux).@"struct".fields) |field| {
        const with_current = @field(candidate, field.name) +
            @field(inputs.current_cell_flux, field.name);
        if (!std.math.isFinite(with_current))
            return error.NonFiniteAdjacentWaterHeatResult;
        const result = with_current - @field(outgoing, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteAdjacentWaterHeatResult;
        @field(candidate, field.name) = result;
    }
    state.net_flux = candidate;
}

fn selectPositiveLayer(inputs: Inputs) usize {
    const neighbor = inputs.positive_neighbor;
    if (inputs.transport_axis != .vertical) return neighbor.initial_layer_index;
    for (
        neighbor.initial_layer_index..neighbor.deepest_layer_index_inclusive + 1,
    ) |layer| {
        if (neighbor.layer_thickness_m[layer] >
            inputs.layer_activity_threshold_m)
        {
            return layer;
        }
    }
    return neighbor.initial_layer_index;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    const neighbor = inputs.positive_neighbor;
    if (neighbor.layer_thickness_m.len == 0 or
        neighbor.flux_by_layer.len != neighbor.layer_thickness_m.len or
        neighbor.initial_layer_index > neighbor.deepest_layer_index_inclusive or
        neighbor.deepest_layer_index_inclusive >= neighbor.layer_thickness_m.len or
        inputs.current_layer_index < inputs.current_top_active_layer_index)
    {
        return error.InvalidAdjacentWaterHeatDimensions;
    }
    if (!std.math.isFinite(inputs.layer_activity_threshold_m) or
        !std.math.isFinite(inputs.current_layer_thickness_m))
    {
        return error.NonFiniteAdjacentWaterHeatInput;
    }
    if (inputs.layer_activity_threshold_m < 0 or
        inputs.current_layer_thickness_m < 0)
    {
        return error.InvalidAdjacentWaterHeatThickness;
    }
    for (neighbor.layer_thickness_m) |thickness_m| {
        if (!std.math.isFinite(thickness_m))
            return error.NonFiniteAdjacentWaterHeatInput;
        if (thickness_m < 0) return error.InvalidAdjacentWaterHeatThickness;
    }
    try validateFlux(inputs.current_cell_flux);
    try validateFlux(neighbor.lake_surface_disappearance_flux);
    for (neighbor.flux_by_layer) |flux| try validateFlux(flux);
    try validateFlux(state.net_flux);
}

fn validateFlux(flux: WaterHeatFlux) !void {
    inline for (@typeInfo(WaterHeatFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteAdjacentWaterHeatInput;
}

fn filledFlux(value: f64) WaterHeatFlux {
    var result: WaterHeatFlux = undefined;
    inline for (@typeInfo(WaterHeatFlux).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn expectFlux(actual: WaterHeatFlux, expected: f64) !void {
    inline for (@typeInfo(WaterHeatFlux).@"struct".fields) |field|
        try std.testing.expectEqual(expected, @field(actual, field.name));
}

fn baseInputs(
    axis: TransportAxis,
    thicknesses: []const f64,
    fluxes: []const WaterHeatFlux,
) Inputs {
    return .{
        .transport_axis = axis,
        .horizontal_exchange = .enabled,
        .layer_activity_threshold_m = 0.1,
        .current_layer_index = 1,
        .current_top_active_layer_index = 0,
        .current_layer_thickness_m = 0.2,
        .current_cell_flux = filledFlux(10),
        .positive_neighbor = .{
            .initial_layer_index = 1,
            .deepest_layer_index_inclusive = thicknesses.len - 1,
            .layer_thickness_m = thicknesses,
            .flux_by_layer = fluxes,
            .lake_surface_disappearance_flux = filledFlux(4),
        },
    };
}

test "horizontal aggregation adds current and subtracts mapped neighbor" {
    const thicknesses = [_]f64{ 0, 0.2, 0.3 };
    const fluxes = [_]WaterHeatFlux{
        filledFlux(99),
        filledFlux(3),
        filledFlux(7),
    };
    var state = State{ .net_flux = filledFlux(100) };
    try aggregate(baseInputs(.column, &thicknesses, &fluxes), &state);
    try expectFlux(state.net_flux, 107);
}

test "vertical scan selects first neighbor layer thicker than threshold" {
    const thicknesses = [_]f64{ 0, 0.05, 0.2, 0.3 };
    const fluxes = [_]WaterHeatFlux{
        filledFlux(99),
        filledFlux(50),
        filledFlux(3),
        filledFlux(7),
    };
    var inputs = baseInputs(.vertical, &thicknesses, &fluxes);
    inputs.positive_neighbor.deepest_layer_index_inclusive = 3;
    var state = State{ .net_flux = filledFlux(100) };
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 107);
}

test "vertical scan without active layer retains mapped neighbor layer" {
    const thicknesses = [_]f64{ 0, 0.05, 0.1 };
    const fluxes = [_]WaterHeatFlux{
        filledFlux(99),
        filledFlux(3),
        filledFlux(7),
    };
    var state = State{ .net_flux = filledFlux(100) };
    try aggregate(baseInputs(.vertical, &thicknesses, &fluxes), &state);
    try expectFlux(state.net_flux, 107);
}

test "vertical current surface uses lake disappearance flux" {
    const thicknesses = [_]f64{ 0, 0.05, 0.2 };
    const fluxes = [_]WaterHeatFlux{
        filledFlux(99),
        filledFlux(50),
        filledFlux(80),
    };
    var inputs = baseInputs(.vertical, &thicknesses, &fluxes);
    inputs.current_layer_index = 0;
    inputs.current_top_active_layer_index = 0;
    var state = State{ .net_flux = filledFlux(100) };
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 106);
}

test "standalone horizontal and inactive current layers do not mutate" {
    const thicknesses = [_]f64{0.2};
    const fluxes = [_]WaterHeatFlux{filledFlux(3)};
    var inputs = baseInputs(.row, &thicknesses, &fluxes);
    inputs.positive_neighbor.initial_layer_index = 0;
    inputs.horizontal_exchange = .standalone;
    inputs.layer_activity_threshold_m = std.math.nan(f64);
    var state = State{ .net_flux = filledFlux(9) };
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 9);

    inputs.horizontal_exchange = .enabled;
    inputs.layer_activity_threshold_m = 0.1;
    inputs.current_layer_thickness_m = 0.1;
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 9);
}

test "equal face fluxes produce exact zero net increment" {
    const shared = filledFlux(7);
    const thicknesses = [_]f64{0.2};
    const fluxes = [_]WaterHeatFlux{shared};
    var inputs = baseInputs(.column, &thicknesses, &fluxes);
    inputs.positive_neighbor.initial_layer_index = 0;
    inputs.current_cell_flux = shared;
    var state = State{ .net_flux = .{} };
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 0);
}

test "invalid input and arithmetic overflow preserve state atomically" {
    const thicknesses = [_]f64{0.2};
    const fluxes = [_]WaterHeatFlux{filledFlux(3)};
    var inputs = baseInputs(.column, &thicknesses, &fluxes);
    inputs.positive_neighbor.initial_layer_index = 0;
    var state = State{ .net_flux = filledFlux(5) };

    inputs.current_layer_thickness_m = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteAdjacentWaterHeatInput,
        aggregate(inputs, &state),
    );
    try expectFlux(state.net_flux, 5);

    inputs.current_layer_thickness_m = 0.2;
    inputs.current_cell_flux = filledFlux(std.math.floatMax(f64));
    state.net_flux = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteAdjacentWaterHeatResult,
        aggregate(inputs, &state),
    );
    try expectFlux(state.net_flux, std.math.floatMax(f64));
}
