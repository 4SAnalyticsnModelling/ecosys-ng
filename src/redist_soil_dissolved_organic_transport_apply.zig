const std = @import("std");

/// Dissolved organic micropore and macropore pools for one layer, K=0-4.
pub const DissolvedOrganicLayer = struct {
    /// OQC[K=0-4] micropore DOC (g C).
    oqc: [5]f64,
    /// OQN[K=0-4] micropore DON (g N).
    oqn: [5]f64,
    /// OQP[K=0-4] micropore DOP (g P).
    oqp: [5]f64,
    /// OQA[K=0-4] micropore DOC acetate (g C).
    oqa: [5]f64,
    /// OQCH[K=0-4] macropore DOC (g C).
    oqch: [5]f64,
    /// OQNH[K=0-4] macropore DON (g N).
    oqnh: [5]f64,
    /// OQPH[K=0-4] macropore DOP (g P).
    oqph: [5]f64,
    /// OQAH[K=0-4] macropore DOC acetate (g C).
    oqah: [5]f64,
};

/// Transport fluxes for dissolved organics at one layer (g step-1).
/// Micropore net flux: TOCFLS + XOCFXS (micro receives +XOCFXS, macro receives -XOCFXS).
pub const DissolvedOrganicFluxes = struct {
    /// TOCFLS[K=0-4]. Net micropore DOC flux.
    tocfls: [5]f64,
    /// TONFLS[K=0-4]. Net micropore DON flux.
    tonfls: [5]f64,
    /// TOPFLS[K=0-4]. Net micropore DOP flux.
    topfls: [5]f64,
    /// TOAFLS[K=0-4]. Net micropore DOC acetate flux.
    toafls: [5]f64,
    /// TOCFHS[K=0-4]. Net macropore DOC flux.
    tocfhs: [5]f64,
    /// TONFHS[K=0-4]. Net macropore DON flux.
    tonfhs: [5]f64,
    /// TOPFHS[K=0-4]. Net macropore DOP flux.
    topfhs: [5]f64,
    /// TOAFHS[K=0-4]. Net macropore DOC acetate flux.
    toafhs: [5]f64,
    /// XOCFXS[K=0-4]. Micro-macropore exchange DOC (+micro, -macro).
    xocfxs: [5]f64,
    /// XONFXS[K=0-4]. Micro-macropore exchange DON.
    xonfxs: [5]f64,
    /// XOPFXS[K=0-4]. Micro-macropore exchange DOP.
    xopfxs: [5]f64,
    /// XOAFXS[K=0-4]. Micro-macropore exchange acetate.
    xoafxs: [5]f64,
};

pub const RuntimeState = struct {
    micropore_doc_g: []f64,
    micropore_don_g: []f64,
    micropore_dop_g: []f64,
    micropore_acetate_g: []f64,
    macropore_doc_g: []f64,
    macropore_don_g: []f64,
    macropore_dop_g: []f64,
    macropore_acetate_g: []f64,
};
pub const RuntimeFluxes = struct {
    micropore_doc_g: []const f64,
    micropore_don_g: []const f64,
    micropore_dop_g: []const f64,
    micropore_acetate_g: []const f64,
    macropore_doc_g: []const f64,
    macropore_don_g: []const f64,
    macropore_dop_g: []const f64,
    macropore_acetate_g: []const f64,
    exchange_doc_g: []const f64,
    exchange_don_g: []const f64,
    exchange_dop_g: []const f64,
    exchange_acetate_g: []const f64,
};

/// Runtime-dimension form of REDIST lines 6080--6096. Flat storage is
/// layer-major then biochemical class K, preserving the enclosing L/K order.
pub fn applyLayers(layer_count: usize, class_count: usize, state: RuntimeState, fluxes: RuntimeFluxes) !void {
    if (layer_count == 0 or class_count == 0) return error.InvalidDissolvedOrganicDimensions;
    const count = std.math.mul(usize, layer_count, class_count) catch return error.DissolvedOrganicDimensionOverflow;
    inline for (@typeInfo(RuntimeState).@"struct".fields) |field|
        if (@field(state, field.name).len != count) return error.DissolvedOrganicDimensionMismatch;
    inline for (@typeInfo(RuntimeFluxes).@"struct".fields) |field|
        if (@field(fluxes, field.name).len != count) return error.DissolvedOrganicDimensionMismatch;
    for (0..layer_count) |layer_index| for (0..class_count) |class_index| {
        const i = layer_index * class_count + class_index;
        inline for (@typeInfo(RuntimeState).@"struct".fields) |field|
            if (!std.math.isFinite(@field(state, field.name)[i])) return error.InvalidDissolvedOrganicPool;
        inline for (@typeInfo(RuntimeFluxes).@"struct".fields) |field|
            if (!std.math.isFinite(@field(fluxes, field.name)[i])) return error.InvalidDissolvedOrganicFlux;
        state.micropore_doc_g[i] = state.micropore_doc_g[i] + fluxes.micropore_doc_g[i] + fluxes.exchange_doc_g[i];
        state.micropore_don_g[i] = state.micropore_don_g[i] + fluxes.micropore_don_g[i] + fluxes.exchange_don_g[i];
        state.micropore_dop_g[i] = state.micropore_dop_g[i] + fluxes.micropore_dop_g[i] + fluxes.exchange_dop_g[i];
        state.micropore_acetate_g[i] = state.micropore_acetate_g[i] + fluxes.micropore_acetate_g[i] + fluxes.exchange_acetate_g[i];
        state.macropore_doc_g[i] = state.macropore_doc_g[i] + fluxes.macropore_doc_g[i] - fluxes.exchange_doc_g[i];
        state.macropore_don_g[i] = state.macropore_don_g[i] + fluxes.macropore_don_g[i] - fluxes.exchange_don_g[i];
        state.macropore_dop_g[i] = state.macropore_dop_g[i] + fluxes.macropore_dop_g[i] - fluxes.exchange_dop_g[i];
        state.macropore_acetate_g[i] = state.macropore_acetate_g[i] + fluxes.macropore_acetate_g[i] - fluxes.exchange_acetate_g[i];
        inline for (@typeInfo(RuntimeState).@"struct".fields) |field|
            if (!std.math.isFinite(@field(state, field.name)[i])) return error.NonFiniteDissolvedOrganicPool;
    };
}

/// Direct translation of REDIST lines 6080--6096 (inner body of DO 125 L loop).
pub fn apply(layer: DissolvedOrganicLayer, fluxes: DissolvedOrganicFluxes) !DissolvedOrganicLayer {
    inline for (@typeInfo(DissolvedOrganicLayer).@"struct".fields) |field|
        for (@field(layer, field.name)) |v| if (!std.math.isFinite(v)) return error.InvalidDissolvedOrganicPool;
    inline for (@typeInfo(DissolvedOrganicFluxes).@"struct".fields) |field| {
        const arr: [5]f64 = @field(fluxes, field.name);
        for (arr) |v| if (!std.math.isFinite(v)) return error.InvalidDissolvedOrganicFlux;
    }

    var new = layer;
    for (0..5) |k| {
        new.oqc[k] = layer.oqc[k] + fluxes.tocfls[k] + fluxes.xocfxs[k];
        new.oqn[k] = layer.oqn[k] + fluxes.tonfls[k] + fluxes.xonfxs[k];
        new.oqp[k] = layer.oqp[k] + fluxes.topfls[k] + fluxes.xopfxs[k];
        new.oqa[k] = layer.oqa[k] + fluxes.toafls[k] + fluxes.xoafxs[k];
        new.oqch[k] = layer.oqch[k] + fluxes.tocfhs[k] - fluxes.xocfxs[k];
        new.oqnh[k] = layer.oqnh[k] + fluxes.tonfhs[k] - fluxes.xonfxs[k];
        new.oqph[k] = layer.oqph[k] + fluxes.topfhs[k] - fluxes.xopfxs[k];
        new.oqah[k] = layer.oqah[k] + fluxes.toafhs[k] - fluxes.xoafxs[k];
    }
    inline for (@typeInfo(DissolvedOrganicLayer).@"struct".fields) |field|
        for (@field(new, field.name)) |v| if (!std.math.isFinite(v)) return error.NonFiniteDissolvedOrganicPool;
    return new;
}

test "REDIST dissolved organic transport adds micropore and exchange fluxes to OQC" {
    var fluxes = std.mem.zeroes(DissolvedOrganicFluxes);
    fluxes.tocfls[2] = 1.0;
    fluxes.xocfxs[2] = 0.5;
    const result = try apply(std.mem.zeroes(DissolvedOrganicLayer), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 1.5), result.oqc[2], 1.0e-15);
}

test "REDIST dissolved organic transport subtracts exchange from macropore OQCH" {
    var fluxes = std.mem.zeroes(DissolvedOrganicFluxes);
    fluxes.tocfhs[1] = 2.0;
    fluxes.xocfxs[1] = 0.8;
    const result = try apply(std.mem.zeroes(DissolvedOrganicLayer), fluxes);
    // OQCH = TOCFHS - XOCFXS = 2.0 - 0.8 = 1.2
    try std.testing.expectApproxEqRel(@as(f64, 1.2), result.oqch[1], 1.0e-15);
}

test "REDIST dissolved organic micropore and macropore exchange are opposite signs" {
    var fluxes = std.mem.zeroes(DissolvedOrganicFluxes);
    fluxes.xocfxs[0] = 3.0;
    const result = try apply(std.mem.zeroes(DissolvedOrganicLayer), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.oqc[0], 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, -3.0), result.oqch[0], 1.0e-15);
}

test "REDIST dissolved organic rejects non-finite transport flux" {
    var bad = std.mem.zeroes(DissolvedOrganicFluxes);
    bad.tonfls[3] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidDissolvedOrganicFlux,
        apply(std.mem.zeroes(DissolvedOrganicLayer), bad),
    );
}

test "REDIST dissolved organic runtime layers preserve opposite exchange signs" {
    var micro_doc = [_]f64{0} ** 6;
    var micro_don = [_]f64{0} ** 6;
    var micro_dop = [_]f64{0} ** 6;
    var micro_a = [_]f64{0} ** 6;
    var macro_doc = [_]f64{0} ** 6;
    var macro_don = [_]f64{0} ** 6;
    var macro_dop = [_]f64{0} ** 6;
    var macro_a = [_]f64{0} ** 6;
    const zero = [_]f64{0} ** 6;
    var exchange_doc = zero;
    var exchange_don = zero;
    var exchange_dop = zero;
    var exchange_a = zero;
    exchange_doc[5] = 2;
    exchange_don[5] = 3;
    exchange_dop[5] = 4;
    exchange_a[5] = 5;
    try applyLayers(2, 3, .{ .micropore_doc_g = &micro_doc, .micropore_don_g = &micro_don, .micropore_dop_g = &micro_dop, .micropore_acetate_g = &micro_a, .macropore_doc_g = &macro_doc, .macropore_don_g = &macro_don, .macropore_dop_g = &macro_dop, .macropore_acetate_g = &macro_a }, .{ .micropore_doc_g = &zero, .micropore_don_g = &zero, .micropore_dop_g = &zero, .micropore_acetate_g = &zero, .macropore_doc_g = &zero, .macropore_don_g = &zero, .macropore_dop_g = &zero, .macropore_acetate_g = &zero, .exchange_doc_g = &exchange_doc, .exchange_don_g = &exchange_don, .exchange_dop_g = &exchange_dop, .exchange_acetate_g = &exchange_a });
    try std.testing.expectEqual(@as(f64, 2), micro_doc[5]);
    try std.testing.expectEqual(@as(f64, -2), macro_doc[5]);
    try std.testing.expectEqual(@as(f64, 3), micro_don[5]);
    try std.testing.expectEqual(@as(f64, -4), macro_dop[5]);
    try std.testing.expectEqual(@as(f64, 5), micro_a[5]);
    try std.testing.expectEqual(@as(f64, -5), macro_a[5]);
}

test "REDIST dissolved organic rejects non-DOC invalid pool and overflow" {
    var layer = std.mem.zeroes(DissolvedOrganicLayer);
    layer.oqnh[2] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidDissolvedOrganicPool, apply(layer, std.mem.zeroes(DissolvedOrganicFluxes)));
    layer = std.mem.zeroes(DissolvedOrganicLayer);
    layer.oqp[1] = std.math.floatMax(f64);
    var flux = std.mem.zeroes(DissolvedOrganicFluxes);
    flux.topfls[1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteDissolvedOrganicPool, apply(layer, flux));
}
