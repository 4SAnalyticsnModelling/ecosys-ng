const std = @import("std");

/// Soil gas pools for one runtime soil layer (g of the named element).
pub const SoilGasPhasePools = struct {
    co2_g_c: f64, // CO2G
    ch4_g_c: f64, // CH4G
    o2_g_o: f64, // OXYG
    n2_g_n: f64, // Z2GG
    n2o_g_n: f64, // Z2OG
    nh3_g_n: f64, // ZNH3G
    h2_g_h: f64, // H2GG
};

/// Diagnostic fluxes overwritten for each layer by REDIST 6484--6490.
pub const GasPhaseAuxFluxes = struct {
    gaseous_o2_flux_g_o: f64, // ROXYF = TOXFLG
    gaseous_co2_flux_g_c: f64, // RCO2F = TCOFLG
    gaseous_ch4_flux_g_c: f64, // RCH4F = TCHFLG
    aqueous_o2_flux_g_o: f64, // ROXYL
    aqueous_ch4_flux_g_c: f64, // RCH4L
};

/// Gas-pool and diagnostic source terms for one timestep.
pub const SoilGasPhaseFluxes = struct {
    co2_gas_transport_g_c: f64, // TCOFLG
    co2_water_air_transfer_g_c: f64, // XCODFG, subtracted
    co2_aerobic_oxidation_g_c: f64, // RCGOX
    co2_methane_oxidation_g_c: f64, // RC4OX
    ch4_gas_transport_g_c: f64, // TCHFLG
    ch4_water_air_transfer_g_c: f64, // XCHDFG, subtracted
    ch4_production_g_c: f64, // RCHOX
    o2_gas_transport_g_o: f64, // TOXFLG
    o2_water_air_transfer_g_o: f64, // XOXDFG, subtracted
    o2_aerobic_oxidation_g_o: f64, // ROGOX, subtracted
    n2_gas_transport_g_n: f64, // TNGFLG
    n2_water_air_transfer_g_n: f64, // XNGDFG, subtracted
    n2o_gas_transport_g_n: f64, // TN2FLG
    n2o_water_air_transfer_g_n: f64, // XN2DFG, subtracted
    nh3_gas_transport_g_n: f64, // TNHFLG
    nh3_nonband_dissolution_g_n: f64, // XN3DFG, subtracted
    nh3_band_dissolution_g_n: f64, // XNBDFG, subtracted
    nh3_transformation_g_n: f64, // TRN3G
    h2_gas_transport_g_h: f64, // THGFLG
    h2_water_air_transfer_g_h: f64, // XHGDFG, subtracted
    o2_aqueous_transport_g_o: f64, // TOXFLS
    o2_subsurface_flux_g_o: f64, // ROXFLU
    o2_pore_exchange_g_o: f64, // XOXFXS
    o2_bubble_flux_g_o: f64, // XOXBBL
    ch4_aqueous_transport_g_c: f64, // TCHFLS
    ch4_subsurface_flux_g_c: f64, // RCHFLU
    ch4_pore_exchange_g_c: f64, // XCHFXS
    ch4_bubble_flux_g_c: f64, // XCHBBL
};

pub const Result = struct {
    pools: SoilGasPhasePools,
    aux: GasPhaseAuxFluxes,
};

/// Direct translation of REDIST 6473--6490, preserving source operation order.
pub fn update(pools: SoilGasPhasePools, fluxes: SoilGasPhaseFluxes) !Result {
    inline for (@typeInfo(SoilGasPhasePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidSoilGasPhasePool;
    inline for (@typeInfo(SoilGasPhaseFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fluxes, field.name)))
            return error.InvalidSoilGasPhaseFlux;

    const new_pools = SoilGasPhasePools{
        .co2_g_c = pools.co2_g_c + fluxes.co2_gas_transport_g_c - fluxes.co2_water_air_transfer_g_c + fluxes.co2_aerobic_oxidation_g_c + fluxes.co2_methane_oxidation_g_c,
        .ch4_g_c = pools.ch4_g_c + fluxes.ch4_gas_transport_g_c - fluxes.ch4_water_air_transfer_g_c + fluxes.ch4_production_g_c - fluxes.co2_methane_oxidation_g_c,
        .o2_g_o = pools.o2_g_o + fluxes.o2_gas_transport_g_o - fluxes.o2_water_air_transfer_g_o - fluxes.o2_aerobic_oxidation_g_o - fluxes.co2_methane_oxidation_g_c * 2.667,
        .n2_g_n = pools.n2_g_n + fluxes.n2_gas_transport_g_n - fluxes.n2_water_air_transfer_g_n,
        .n2o_g_n = pools.n2o_g_n + fluxes.n2o_gas_transport_g_n - fluxes.n2o_water_air_transfer_g_n,
        .nh3_g_n = pools.nh3_g_n + fluxes.nh3_gas_transport_g_n - fluxes.nh3_nonband_dissolution_g_n - fluxes.nh3_band_dissolution_g_n + fluxes.nh3_transformation_g_n,
        .h2_g_h = pools.h2_g_h + fluxes.h2_gas_transport_g_h - fluxes.h2_water_air_transfer_g_h,
    };
    inline for (@typeInfo(SoilGasPhasePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(new_pools, field.name)))
            return error.NonFiniteSoilGasPhasePool;

    const aux = GasPhaseAuxFluxes{
        .gaseous_o2_flux_g_o = fluxes.o2_gas_transport_g_o,
        .gaseous_co2_flux_g_c = fluxes.co2_gas_transport_g_c,
        .gaseous_ch4_flux_g_c = fluxes.ch4_gas_transport_g_c,
        .aqueous_o2_flux_g_o = fluxes.o2_aqueous_transport_g_o + fluxes.o2_subsurface_flux_g_o + fluxes.o2_pore_exchange_g_o + fluxes.o2_bubble_flux_g_o,
        .aqueous_ch4_flux_g_c = fluxes.ch4_aqueous_transport_g_c + fluxes.ch4_subsurface_flux_g_c + fluxes.ch4_pore_exchange_g_c + fluxes.ch4_bubble_flux_g_c,
    };
    inline for (@typeInfo(GasPhaseAuxFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(aux, field.name)))
            return error.NonFiniteSoilGasAuxFlux;
    return .{ .pools = new_pools, .aux = aux };
}

/// Applies the REDIST DO 125 body to the runtime-configured soil layer span.
pub fn updateLayers(
    pools_by_layer: []SoilGasPhasePools,
    fluxes_by_layer: []const SoilGasPhaseFluxes,
    aux_by_layer: []GasPhaseAuxFluxes,
) !void {
    if (pools_by_layer.len == 0 or
        fluxes_by_layer.len != pools_by_layer.len or
        aux_by_layer.len != pools_by_layer.len)
        return error.SoilGasPhaseDimensionMismatch;

    for (pools_by_layer, fluxes_by_layer, aux_by_layer) |*pools, fluxes, *aux| {
        const result = try update(pools.*, fluxes);
        pools.* = result.pools;
        aux.* = result.aux;
    }
}

test "REDIST soil gas pools preserve every source term and sign" {
    var pools = std.mem.zeroes(SoilGasPhasePools);
    pools.co2_g_c = 10.0;
    var fluxes = std.mem.zeroes(SoilGasPhaseFluxes);
    fluxes.co2_gas_transport_g_c = 1.0;
    fluxes.co2_water_air_transfer_g_c = 2.0;
    fluxes.co2_aerobic_oxidation_g_c = 3.0;
    fluxes.co2_methane_oxidation_g_c = 4.0;
    fluxes.o2_gas_transport_g_o = 20.0;
    fluxes.o2_water_air_transfer_g_o = 1.0;
    fluxes.o2_aerobic_oxidation_g_o = 2.0;
    fluxes.nh3_gas_transport_g_n = 8.0;
    fluxes.nh3_nonband_dissolution_g_n = 1.0;
    fluxes.nh3_band_dissolution_g_n = 2.0;
    fluxes.nh3_transformation_g_n = 3.0;

    const result = try update(pools, fluxes);
    try std.testing.expectEqual(@as(f64, 16.0), result.pools.co2_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 6.332), result.pools.o2_g_o, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 8.0), result.pools.nh3_g_n);
}

test "REDIST soil gas auxiliary fluxes use literal transport sums" {
    var fluxes = std.mem.zeroes(SoilGasPhaseFluxes);
    fluxes.o2_gas_transport_g_o = 1.0;
    fluxes.co2_gas_transport_g_c = 2.0;
    fluxes.ch4_gas_transport_g_c = 3.0;
    fluxes.o2_aqueous_transport_g_o = 4.0;
    fluxes.o2_subsurface_flux_g_o = 5.0;
    fluxes.o2_pore_exchange_g_o = 6.0;
    fluxes.o2_bubble_flux_g_o = 7.0;
    fluxes.ch4_aqueous_transport_g_c = 8.0;
    fluxes.ch4_subsurface_flux_g_c = 9.0;
    fluxes.ch4_pore_exchange_g_c = 10.0;
    fluxes.ch4_bubble_flux_g_c = 11.0;
    const result = try update(std.mem.zeroes(SoilGasPhasePools), fluxes);
    try std.testing.expectEqual(@as(f64, 1.0), result.aux.gaseous_o2_flux_g_o);
    try std.testing.expectEqual(@as(f64, 2.0), result.aux.gaseous_co2_flux_g_c);
    try std.testing.expectEqual(@as(f64, 3.0), result.aux.gaseous_ch4_flux_g_c);
    try std.testing.expectEqual(@as(f64, 22.0), result.aux.aqueous_o2_flux_g_o);
    try std.testing.expectEqual(@as(f64, 38.0), result.aux.aqueous_ch4_flux_g_c);
}

test "REDIST soil gas runtime layer traversal remains independent" {
    var pools = [_]SoilGasPhasePools{ std.mem.zeroes(SoilGasPhasePools), std.mem.zeroes(SoilGasPhasePools) };
    var fluxes = [_]SoilGasPhaseFluxes{ std.mem.zeroes(SoilGasPhaseFluxes), std.mem.zeroes(SoilGasPhaseFluxes) };
    var aux: [2]GasPhaseAuxFluxes = undefined;
    fluxes[0].n2_gas_transport_g_n = 2.0;
    fluxes[1].h2_gas_transport_g_h = 3.0;
    try updateLayers(&pools, &fluxes, &aux);
    try std.testing.expectEqual(@as(f64, 2.0), pools[0].n2_g_n);
    try std.testing.expectEqual(@as(f64, 0.0), pools[1].n2_g_n);
    try std.testing.expectEqual(@as(f64, 3.0), pools[1].h2_g_h);
}

test "REDIST soil gas rejects dimensions non-finite input and overflow" {
    var pools = [_]SoilGasPhasePools{std.mem.zeroes(SoilGasPhasePools)};
    const no_fluxes: [0]SoilGasPhaseFluxes = .{};
    var aux: [1]GasPhaseAuxFluxes = undefined;
    try std.testing.expectError(error.SoilGasPhaseDimensionMismatch, updateLayers(&pools, &no_fluxes, &aux));

    var invalid = std.mem.zeroes(SoilGasPhaseFluxes);
    invalid.n2_gas_transport_g_n = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilGasPhaseFlux, update(pools[0], invalid));

    pools[0].n2_g_n = std.math.floatMax(f64);
    var overflowing = std.mem.zeroes(SoilGasPhaseFluxes);
    overflowing.n2_gas_transport_g_n = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilGasPhasePool, update(pools[0], overflowing));
}
