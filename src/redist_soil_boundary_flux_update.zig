const std = @import("std");

/// Gas pools indexed by the original runtime soil-layer number, including 0.
pub const GasPools = struct {
    co2_g_c: f64,
    ch4_g_c: f64,
    o2_g_o: f64,
    n2_g_n: f64,
    n2o_g_n: f64,
    nh3_g_n: f64,
    h2_g_h: f64,
};

pub const BubbleFluxes = struct {
    co2_g_c: f64, // XCOBBL
    ch4_g_c: f64, // XCHBBL
    o2_g_o: f64, // XOXBBL
    n2_g_n: f64, // XNGBBL
    n2o_g_n: f64, // XN2BBL
    nh3_nonband_g_n: f64, // XN3BBL
    nh3_band_g_n: f64, // XNBBBL
    h2_g_h: f64, // XHGBBL
};

pub const LayerInputs = struct {
    heat_conduction_megajoules: f64, // THFLFL
    heat_convection_megajoules: f64, // THFLVL
    root_water_heat_megajoules: f64, // TUPHT
    root_co2_atmosphere_g_c: f64, // TCOFLA
    root_ch4_atmosphere_g_c: f64, // TCHFLA
    root_o2_atmosphere_g_o: f64, // TOXFLA
    root_n2o_atmosphere_g_n: f64, // TN2FLA
    root_nh3_atmosphere_g_n: f64, // TNHFLA
    bubbles: BubbleFluxes,
    root_co2_exchange_g_c: f64, // TCO2P
    soil_co2_exchange_g_c: f64, // TCO2S
    co2_equilibrium_change_g_c: f64, // TRCO2
    root_o2_uptake_g_o: f64, // RUPOXO
    root_primary_o2_uptake_g_o: f64, // TUPOXP
    root_secondary_o2_uptake_g_o: f64, // TUPOXS
    microbial_o2_uptake_g_o: f64, // ROGOX
    methane_oxidation_g_c: f64, // RC4OX
    hydrogen_output_g_h: f64, // RH2GO
    root_hydrogen_uptake_g_h: f64, // TUPHGS
};

/// Grid-cell and run accumulators updated in literal REDIST source order.
pub const Accumulators = struct {
    heat_input_megajoules: f64, // HEATIN
    total_lost_co2_g_c: f64, // TLCO2G
    subsurface_co2_g_c: f64, // UCO2S
    soil_o2_g_o: f64, // OXYGSO
    total_lost_n_g_n: f64, // TLN2G
    total_lost_h2_g_h: f64, // TLH2G
    co2_input_g_c: f64, // CO2GIN
    co2_output_g_c: f64, // TCOU
    co2_net_g_c: f64, // XCNET
    ch4_net_g_c: f64, // XHNET
    hourly_co2_g_c: f64, // HCO2G
    cumulative_co2_g_c: f64, // UCO2G
    hourly_ch4_g_c: f64, // HCH4G
    cumulative_ch4_g_c: f64, // UCH4G
    root_soil_co2_g_c: f64, // UCOP
    o2_net_g_o: f64, // XONET
    o2_input_g_o: f64, // OXYGIN
    o2_output_g_o: f64, // OXYGOU
    cumulative_o2_g_o: f64, // UOXYG
    hourly_o2_g_o: f64, // HOXYG
    h2_input_g_h: f64, // H2GIN
    h2_output_g_h: f64, // H2GOU
    nitrogen_gas_input_g_n: f64, // ZN2GIN
    cumulative_n2o_g_n: f64, // UN2OG
    hourly_n2o_g_n: f64, // HN2OG
    cumulative_nh3_g_n: f64, // UNH3G
    hourly_nh3_g_n: f64, // HNH3G
    cumulative_h2_g_h: f64, // UH2GG
    cumulative_co2_equilibrium_g_c: f64, // TXCO2
};

fn allFinite(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (@TypeOf(field_value) == f64 and !std.math.isFinite(field_value)) return false;
    }
    return true;
}

/// Direct translation of REDIST 6529--6594 for one source layer.
pub fn updateLayer(
    accumulators: *Accumulators,
    gas_pools_by_fortran_layer: []GasPools,
    source_layer_number: usize,
    gas_boundary_layer_number: usize,
    input: LayerInputs,
) !void {
    if (source_layer_number >= gas_pools_by_fortran_layer.len or
        gas_boundary_layer_number >= gas_pools_by_fortran_layer.len)
        return error.SoilBoundaryLayerOutOfRange;
    if (!allFinite(accumulators.*) or !allFinite(input) or !allFinite(input.bubbles))
        return error.InvalidSoilBoundaryFlux;
    for (gas_pools_by_fortran_layer) |pools|
        if (!allFinite(pools)) return error.InvalidSoilBoundaryGasPool;

    var next = accumulators.*;
    var destination = gas_pools_by_fortran_layer[@min(source_layer_number, gas_boundary_layer_number)];
    next.heat_input_megajoules = next.heat_input_megajoules + input.heat_conduction_megajoules + input.heat_convection_megajoules + input.root_water_heat_megajoules;
    var co2_input = input.root_co2_atmosphere_g_c;
    var ch4_input = input.root_ch4_atmosphere_g_c;
    var o2_input = input.root_o2_atmosphere_g_o;
    var h2_input: f64 = 0.0;
    var n2_input: f64 = 0.0;
    var n2o_input = input.root_n2o_atmosphere_g_n;
    var nh3_input = input.root_nh3_atmosphere_g_n;

    if (gas_boundary_layer_number == 0) {
        co2_input = co2_input + input.bubbles.co2_g_c;
        ch4_input = ch4_input + input.bubbles.ch4_g_c;
        o2_input = o2_input + input.bubbles.o2_g_o;
        n2_input = n2_input + input.bubbles.n2_g_n;
        n2o_input = n2o_input + input.bubbles.n2o_g_n;
        nh3_input = nh3_input + input.bubbles.nh3_nonband_g_n + input.bubbles.nh3_band_g_n;
        h2_input = h2_input + input.bubbles.h2_g_h;
    } else {
        destination.co2_g_c = destination.co2_g_c - input.bubbles.co2_g_c;
        destination.ch4_g_c = destination.ch4_g_c - input.bubbles.ch4_g_c;
        destination.o2_g_o = destination.o2_g_o - input.bubbles.o2_g_o;
        destination.n2_g_n = destination.n2_g_n - input.bubbles.n2_g_n;
        destination.n2o_g_n = destination.n2o_g_n - input.bubbles.n2o_g_n;
        destination.nh3_g_n = destination.nh3_g_n - input.bubbles.nh3_nonband_g_n - input.bubbles.nh3_band_g_n;
        destination.h2_g_h = destination.h2_g_h - input.bubbles.h2_g_h;
        if (gas_boundary_layer_number < source_layer_number) {
            next.total_lost_co2_g_c = next.total_lost_co2_g_c - input.bubbles.co2_g_c - input.bubbles.ch4_g_c;
            next.subsurface_co2_g_c = next.subsurface_co2_g_c - input.bubbles.co2_g_c - input.bubbles.ch4_g_c;
            next.soil_o2_g_o = next.soil_o2_g_o - input.bubbles.o2_g_o;
            next.total_lost_n_g_n = next.total_lost_n_g_n - input.bubbles.n2_g_n - input.bubbles.n2o_g_n - input.bubbles.nh3_nonband_g_n - input.bubbles.nh3_band_g_n;
            next.total_lost_h2_g_h = next.total_lost_h2_g_h - input.bubbles.h2_g_h;
        }
    }

    next.co2_input_g_c = next.co2_input_g_c + co2_input + ch4_input;
    const co2_output = input.root_co2_exchange_g_c + input.soil_co2_exchange_g_c - input.co2_equilibrium_change_g_c;
    next.co2_output_g_c = next.co2_output_g_c + co2_output;
    next.co2_net_g_c = next.co2_net_g_c + co2_input;
    next.ch4_net_g_c = next.ch4_net_g_c + ch4_input;
    next.hourly_co2_g_c = next.hourly_co2_g_c + co2_input;
    next.cumulative_co2_g_c = next.cumulative_co2_g_c + co2_input;
    next.hourly_ch4_g_c = next.hourly_ch4_g_c + ch4_input;
    next.cumulative_ch4_g_c = next.cumulative_ch4_g_c + ch4_input;
    next.root_soil_co2_g_c = next.root_soil_co2_g_c + input.root_co2_exchange_g_c + input.soil_co2_exchange_g_c;
    next.o2_net_g_o = next.o2_net_g_o + o2_input;
    next.o2_input_g_o = next.o2_input_g_o + o2_input;
    const o2_output = input.root_o2_uptake_g_o + input.root_primary_o2_uptake_g_o + input.root_secondary_o2_uptake_g_o + input.microbial_o2_uptake_g_o + input.methane_oxidation_g_c * 2.667;
    next.o2_output_g_o = next.o2_output_g_o + o2_output;
    next.cumulative_o2_g_o = next.cumulative_o2_g_o + o2_input;
    next.hourly_o2_g_o = next.hourly_o2_g_o + o2_input;
    next.h2_input_g_h = next.h2_input_g_h + h2_input;
    const h2_output = input.hydrogen_output_g_h + input.root_hydrogen_uptake_g_h;
    next.h2_output_g_h = next.h2_output_g_h + h2_output;
    next.nitrogen_gas_input_g_n = next.nitrogen_gas_input_g_n + n2_input + n2o_input + nh3_input;
    next.cumulative_n2o_g_n = next.cumulative_n2o_g_n + n2o_input;
    next.hourly_n2o_g_n = next.hourly_n2o_g_n + n2o_input;
    next.cumulative_nh3_g_n = next.cumulative_nh3_g_n + nh3_input;
    next.hourly_nh3_g_n = next.hourly_nh3_g_n + nh3_input;
    next.cumulative_h2_g_h = next.cumulative_h2_g_h + h2_input;
    next.cumulative_co2_equilibrium_g_c = next.cumulative_co2_equilibrium_g_c + input.co2_equilibrium_change_g_c;

    if (!allFinite(next) or !allFinite(destination)) return error.NonFiniteSoilBoundaryFlux;
    accumulators.* = next;
    if (gas_boundary_layer_number != 0)
        gas_pools_by_fortran_layer[@min(source_layer_number, gas_boundary_layer_number)] = destination;
}

/// Executes consecutive runtime layers serially because each layer mutates
/// shared cell accumulators and may target the same `MIN(L,LG)` gas pool.
pub fn updateLayers(
    accumulators: *Accumulators,
    gas_pools_by_fortran_layer: []GasPools,
    first_source_layer_number: usize,
    gas_boundary_layer_number: usize,
    inputs_by_layer: []const LayerInputs,
) !void {
    if (inputs_by_layer.len == 0 or
        first_source_layer_number >= gas_pools_by_fortran_layer.len or
        inputs_by_layer.len > gas_pools_by_fortran_layer.len - first_source_layer_number)
        return error.SoilBoundaryDimensionMismatch;

    for (inputs_by_layer, first_source_layer_number..) |input, source_layer_number|
        try updateLayer(
            accumulators,
            gas_pools_by_fortran_layer,
            source_layer_number,
            gas_boundary_layer_number,
            input,
        );
}

test "REDIST boundary bubbling with LG zero becomes atmospheric input" {
    var accumulators = std.mem.zeroes(Accumulators);
    var pools = [_]GasPools{ std.mem.zeroes(GasPools), std.mem.zeroes(GasPools) };
    var input = std.mem.zeroes(LayerInputs);
    input.root_co2_atmosphere_g_c = 1.0;
    input.bubbles.co2_g_c = 2.0;
    input.bubbles.ch4_g_c = 3.0;
    input.bubbles.n2_g_n = 4.0;
    input.bubbles.nh3_nonband_g_n = 5.0;
    input.bubbles.nh3_band_g_n = 6.0;
    try updateLayer(&accumulators, &pools, 1, 0, input);
    try std.testing.expectEqual(@as(f64, 6.0), accumulators.co2_input_g_c);
    try std.testing.expectEqual(@as(f64, 15.0), accumulators.nitrogen_gas_input_g_n);
    try std.testing.expectEqual(std.mem.zeroes(GasPools), pools[0]);
}

test "REDIST boundary bubbling with buried LG subtracts destination and ledgers" {
    var accumulators = std.mem.zeroes(Accumulators);
    var pools = [_]GasPools{ std.mem.zeroes(GasPools), std.mem.zeroes(GasPools), std.mem.zeroes(GasPools) };
    pools[1] = .{ .co2_g_c = 10, .ch4_g_c = 10, .o2_g_o = 10, .n2_g_n = 10, .n2o_g_n = 10, .nh3_g_n = 10, .h2_g_h = 10 };
    var input = std.mem.zeroes(LayerInputs);
    input.bubbles = .{ .co2_g_c = 1, .ch4_g_c = 2, .o2_g_o = 3, .n2_g_n = 4, .n2o_g_n = 5, .nh3_nonband_g_n = 1, .nh3_band_g_n = 2, .h2_g_h = 6 };
    try updateLayer(&accumulators, &pools, 2, 1, input);
    try std.testing.expectEqual(@as(f64, 9.0), pools[1].co2_g_c);
    try std.testing.expectEqual(@as(f64, 7.0), pools[1].nh3_g_n);
    try std.testing.expectEqual(@as(f64, -3.0), accumulators.total_lost_co2_g_c);
    try std.testing.expectEqual(@as(f64, -12.0), accumulators.total_lost_n_g_n);
}

test "REDIST boundary flux preserves diagnostic terms and source signs" {
    var accumulators = std.mem.zeroes(Accumulators);
    var pools = [_]GasPools{ std.mem.zeroes(GasPools), std.mem.zeroes(GasPools) };
    var input = std.mem.zeroes(LayerInputs);
    input.heat_conduction_megajoules = 1;
    input.heat_convection_megajoules = 2;
    input.root_water_heat_megajoules = 3;
    input.root_co2_exchange_g_c = 4;
    input.soil_co2_exchange_g_c = 5;
    input.co2_equilibrium_change_g_c = 6;
    input.root_o2_uptake_g_o = 1;
    input.root_primary_o2_uptake_g_o = 2;
    input.root_secondary_o2_uptake_g_o = 3;
    input.microbial_o2_uptake_g_o = 4;
    input.methane_oxidation_g_c = 5;
    try updateLayer(&accumulators, &pools, 1, 0, input);
    try std.testing.expectEqual(@as(f64, 6.0), accumulators.heat_input_megajoules);
    try std.testing.expectEqual(@as(f64, 3.0), accumulators.co2_output_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 23.335), accumulators.o2_output_g_o, 1e-12);
    try std.testing.expectEqual(@as(f64, 6.0), accumulators.cumulative_co2_equilibrium_g_c);
}

test "REDIST boundary flux rejects indexes invalid values and overflow" {
    var accumulators = std.mem.zeroes(Accumulators);
    var pools = [_]GasPools{std.mem.zeroes(GasPools)};
    const input = std.mem.zeroes(LayerInputs);
    try std.testing.expectError(error.SoilBoundaryLayerOutOfRange, updateLayer(&accumulators, &pools, 1, 0, input));
    var bad = input;
    bad.bubbles.co2_g_c = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilBoundaryFlux, updateLayer(&accumulators, &pools, 0, 0, bad));
    accumulators.heat_input_megajoules = std.math.floatMax(f64);
    var overflow = input;
    overflow.heat_conduction_megajoules = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundaryFlux, updateLayer(&accumulators, &pools, 0, 0, overflow));
}

test "REDIST boundary runtime layers execute in ascending source order" {
    var accumulators = std.mem.zeroes(Accumulators);
    var pools = [_]GasPools{
        std.mem.zeroes(GasPools),
        std.mem.zeroes(GasPools),
        std.mem.zeroes(GasPools),
    };
    var inputs = [_]LayerInputs{ std.mem.zeroes(LayerInputs), std.mem.zeroes(LayerInputs) };
    inputs[0].heat_conduction_megajoules = 1.0e16;
    inputs[1].heat_conduction_megajoules = -1.0e16;
    accumulators.heat_input_megajoules = 1.0;
    try updateLayers(&accumulators, &pools, 1, 0, &inputs);
    try std.testing.expectEqual(@as(f64, 0.0), accumulators.heat_input_megajoules);

    const no_inputs: [0]LayerInputs = .{};
    try std.testing.expectError(
        error.SoilBoundaryDimensionMismatch,
        updateLayers(&accumulators, &pools, 1, 0, &no_inputs),
    );
}
