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

/// Convective plus diffusive gas transfer, in REDIST source order.
pub const GasFlux = struct {
    carbon_dioxide_g_per_step: f64 = 0,
    methane_g_per_step: f64 = 0,
    oxygen_g_per_step: f64 = 0,
    dinitrogen_g_per_step: f64 = 0,
    nitrous_oxide_g_per_step: f64 = 0,
    ammonia_g_per_step: f64 = 0,
    hydrogen_g_per_step: f64 = 0,
};

pub const Inputs = struct {
    transport_axis: TransportAxis,
    horizontal_exchange: HorizontalExchange,
    current_layer_thickness_m: f64,
    layer_activity_threshold_m: f64,
    current_cell_flux: GasFlux,
    positive_neighbor_flux: GasFlux,
};

pub const State = struct {
    net_flux: GasFlux,
};

/// Aggregates net gas transfer between adjacent grid cells.
///
/// Traceability: REDIST.F lines 3516--3529 under enclosing gates 3350 and
/// 3363. Each gas remains g per model step. In source order, current-cell CO2,
/// CH4, O2, N2, N2O, NH3, and H2 face fluxes are added and the already
/// selected positive-neighbor `N6` face fluxes are subtracted. Horizontal
/// standalone cells are skipped while vertical transfer remains enabled.
/// Candidate state commits only after every result is finite.
pub fn aggregate(inputs: Inputs, state: *State) !void {
    if (inputs.transport_axis != .vertical and
        inputs.horizontal_exchange == .standalone)
    {
        return;
    }
    try validateInputs(inputs, state.*);
    if (inputs.current_layer_thickness_m <= inputs.layer_activity_threshold_m)
        return;

    var candidate = state.net_flux;
    inline for (@typeInfo(GasFlux).@"struct".fields) |field| {
        const with_current = @field(candidate, field.name) +
            @field(inputs.current_cell_flux, field.name);
        if (!std.math.isFinite(with_current))
            return error.NonFiniteAdjacentGasResult;
        const result = with_current -
            @field(inputs.positive_neighbor_flux, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteAdjacentGasResult;
        @field(candidate, field.name) = result;
    }
    state.net_flux = candidate;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    if (!std.math.isFinite(inputs.current_layer_thickness_m) or
        !std.math.isFinite(inputs.layer_activity_threshold_m))
    {
        return error.NonFiniteAdjacentGasInput;
    }
    if (inputs.current_layer_thickness_m < 0 or
        inputs.layer_activity_threshold_m < 0)
    {
        return error.InvalidAdjacentGasThickness;
    }
    try validateFlux(inputs.current_cell_flux);
    try validateFlux(inputs.positive_neighbor_flux);
    try validateFlux(state.net_flux);
}

fn validateFlux(flux: GasFlux) !void {
    inline for (@typeInfo(GasFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteAdjacentGasInput;
}

fn filledFlux(value: f64) GasFlux {
    var result: GasFlux = undefined;
    inline for (@typeInfo(GasFlux).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn expectFlux(actual: GasFlux, expected: f64) !void {
    inline for (@typeInfo(GasFlux).@"struct".fields) |field|
        try std.testing.expectEqual(expected, @field(actual, field.name));
}

fn baseInputs(current: f64, positive: f64) Inputs {
    return .{
        .transport_axis = .column,
        .horizontal_exchange = .enabled,
        .current_layer_thickness_m = 0.2,
        .layer_activity_threshold_m = 0.1,
        .current_cell_flux = filledFlux(current),
        .positive_neighbor_flux = filledFlux(positive),
    };
}

test "all seven gases follow current minus positive-neighbor flux" {
    var state = State{ .net_flux = filledFlux(100) };
    try aggregate(baseInputs(5, 2), &state);
    try expectFlux(state.net_flux, 103);
}

test "reversed adjacent faces are exactly antisymmetric" {
    var first = State{ .net_flux = .{} };
    var second = State{ .net_flux = .{} };
    try aggregate(baseInputs(7, 3), &first);
    try aggregate(baseInputs(3, 7), &second);
    inline for (@typeInfo(GasFlux).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(first.net_flux, field.name) +
                @field(second.net_flux, field.name),
        );
}

test "equal face fluxes produce exact zero increment" {
    var state = State{ .net_flux = .{} };
    try aggregate(baseInputs(9, 9), &state);
    try expectFlux(state.net_flux, 0);
}

test "standalone horizontal and inactive layers do not mutate" {
    var state = State{ .net_flux = filledFlux(9) };
    var inputs = baseInputs(100, 0);
    inputs.horizontal_exchange = .standalone;
    inputs.current_layer_thickness_m = std.math.nan(f64);
    inputs.layer_activity_threshold_m = std.math.nan(f64);
    inputs.current_cell_flux = filledFlux(std.math.nan(f64));
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 9);

    inputs.horizontal_exchange = .enabled;
    inputs.current_layer_thickness_m = 0.1;
    inputs.layer_activity_threshold_m = 0.1;
    inputs.current_cell_flux = filledFlux(100);
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 9);

    inputs.transport_axis = .vertical;
    inputs.horizontal_exchange = .standalone;
    inputs.current_layer_thickness_m = 0.2;
    try aggregate(inputs, &state);
    try expectFlux(state.net_flux, 109);
}

test "invalid input and arithmetic overflow preserve state atomically" {
    var state = State{ .net_flux = filledFlux(5) };
    var inputs = baseInputs(3, 1);
    inputs.current_cell_flux.hydrogen_g_per_step = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteAdjacentGasInput,
        aggregate(inputs, &state),
    );
    try expectFlux(state.net_flux, 5);

    inputs.current_cell_flux = filledFlux(std.math.floatMax(f64));
    inputs.positive_neighbor_flux = filledFlux(-std.math.floatMax(f64));
    state.net_flux = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteAdjacentGasResult,
        aggregate(inputs, &state),
    );
    try expectFlux(state.net_flux, std.math.floatMax(f64));
}
