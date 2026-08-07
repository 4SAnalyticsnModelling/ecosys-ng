const std = @import("std");

pub const SpeciesTerms = struct {
    transport_g: f64 = 0,
    dissolution_g: f64 = 0,
    transformation_loss_g: f64 = 0,
    root_exchange_loss_g: f64 = 0,
    subsurface_input_g: f64 = 0,
    micropore_macropore_exchange_g: f64 = 0,
    additional_loss_g: f64 = 0,
    chemical_transformation_g: f64 = 0,
    bubble_flux_g: f64 = 0,
};
pub const Pools = struct { co2_c_g: f64, ch4_c_g: f64, oxygen_g: f64, n2_n_g: f64, n2o_n_g: f64, h2_h_g: f64 };
pub const Fluxes = struct { co2: SpeciesTerms, ch4: SpeciesTerms, oxygen: SpeciesTerms, n2: SpeciesTerms, n2o: SpeciesTerms, h2: SpeciesTerms };

fn update(pool: f64, t: SpeciesTerms) !f64 {
    if (!std.math.isFinite(pool)) return error.InvalidSoilGasSolutePool;
    inline for (@typeInfo(SpeciesTerms).@"struct".fields) |f| if (!std.math.isFinite(@field(t, f.name))) return error.InvalidSoilGasSoluteFlux;
    const result = pool + t.transport_g + t.dissolution_g - t.transformation_loss_g -
        t.root_exchange_loss_g + t.subsurface_input_g + t.micropore_macropore_exchange_g -
        t.additional_loss_g + t.chemical_transformation_g + t.bubble_flux_g;
    if (!std.math.isFinite(result)) return error.NonFiniteSoilGasSolutePool;
    return result;
}

/// Exact source-order translation of redist.f lines 6140--6157 for one layer.
/// Callers leave terms absent from a legacy species equation at zero.
pub fn apply(pools: Pools, fluxes: Fluxes) !Pools {
    return .{
        .co2_c_g = try update(pools.co2_c_g, fluxes.co2),
        .ch4_c_g = try update(pools.ch4_c_g, fluxes.ch4),
        .oxygen_g = try update(pools.oxygen_g, fluxes.oxygen),
        .n2_n_g = try update(pools.n2_n_g, fluxes.n2),
        .n2o_n_g = try update(pools.n2o_n_g, fluxes.n2o),
        .h2_h_g = try update(pools.h2_h_g, fluxes.h2),
    };
}

/// Runtime NU..NL traversal in source layer order.
pub fn applyLayers(pools: []Pools, fluxes: []const Fluxes) !void {
    if (pools.len == 0 or pools.len != fluxes.len) return error.SoilGasSoluteDimensionMismatch;
    for (pools, fluxes) |*pool, layer_fluxes| pool.* = try apply(pool.*, layer_fluxes);
}

test "REDIST soil gas solute preserves every source sign" {
    const terms = SpeciesTerms{ .transport_g = 1, .dissolution_g = 2, .transformation_loss_g = 3, .root_exchange_loss_g = 4, .subsurface_input_g = 5, .micropore_macropore_exchange_g = 6, .additional_loss_g = 7, .chemical_transformation_g = 8, .bubble_flux_g = 9 };
    const result = try apply(std.mem.zeroes(Pools), .{ .co2 = terms, .ch4 = terms, .oxygen = terms, .n2 = terms, .n2o = terms, .h2 = terms });
    inline for (@typeInfo(Pools).@"struct".fields) |f| try std.testing.expectEqual(@as(f64, 17), @field(result, f.name));
}

test "REDIST soil gas solute runtime layers retain species independence" {
    var pools = [_]Pools{ std.mem.zeroes(Pools), std.mem.zeroes(Pools) };
    var fluxes = [_]Fluxes{ std.mem.zeroes(Fluxes), std.mem.zeroes(Fluxes) };
    fluxes[0].co2.chemical_transformation_g = 2;
    fluxes[1].n2.additional_loss_g = 3;
    try applyLayers(&pools, &fluxes);
    try std.testing.expectEqual(@as(f64, 2), pools[0].co2_c_g);
    try std.testing.expectEqual(@as(f64, -3), pools[1].n2_n_g);
    try std.testing.expectEqual(@as(f64, 0), pools[1].co2_c_g);
}

test "REDIST soil gas solute rejects invalid and overflow" {
    var fluxes = std.mem.zeroes(Fluxes);
    fluxes.h2.bubble_flux_g = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilGasSoluteFlux, apply(std.mem.zeroes(Pools), fluxes));
    var pools = std.mem.zeroes(Pools);
    pools.oxygen_g = std.math.floatMax(f64);
    fluxes = std.mem.zeroes(Fluxes);
    fluxes.oxygen.transport_g = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilGasSolutePool, apply(pools, fluxes));
}
