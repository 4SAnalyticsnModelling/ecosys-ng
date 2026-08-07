const std = @import("std");

pub const Pools = struct {
    layer_count: usize,
    microbial_c_n_p: [3][]const f64, // each 6*7*3*layers
    residue_c_n_p: [3][]const f64, // each 5*2*layers
    /// OQC,OQCH,OHC,OQA,OQAH,OHA,OQN,OQNH,OHN,OQP,OQPH,OHP.
    soluble_and_adsorbed: [12][]const f64, // each 5*layers
    /// OSC,OSA,OSN,OSP; each 5 M * 5 K * layers.
    som_c_a_n_p: [4][]const f64,
};
pub const Totals = struct {
    excluding_fraction_4_carbon_g_c: f64 = 0,
    excluding_fraction_4_nitrogen_g_n: f64 = 0,
    excluding_fraction_4_phosphorus_g_p: f64 = 0, // DC/DN/DP
    all_fractions_carbon_g_c: f64 = 0,
    all_fractions_nitrogen_g_n: f64 = 0,
    all_fractions_phosphorus_g_p: f64 = 0, // OC/ON/OP
    excluding_fraction_4_charcoal_carbon_g_c: f64 = 0,
    excluding_fraction_4_charcoal_nitrogen_g_n: f64 = 0,
    excluding_fraction_4_charcoal_phosphorus_g_p: f64 = 0, // DCC/DNC/DPC
    fraction_4_charcoal_carbon_g_c: f64 = 0,
    fraction_4_charcoal_nitrogen_g_n: f64 = 0,
    fraction_4_charcoal_phosphorus_g_p: f64 = 0, // OCC/ONC/OPC
};
pub const Ledgers = struct {
    organic_carbon_g_c: f64,
    organic_nitrogen_g_n: f64,
    charcoal_organic_carbon_g_c: f64,
    charcoal_organic_nitrogen_g_n: f64,
    residue_organic_carbon_g_c: f64,
};
pub const Result = struct { totals: Totals, ledgers: Ledgers };

fn index(layers: usize, pool: usize, layer: usize) usize {
    return pool * layers + layer;
}
fn finiteNonnegative(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value) or value < 0) return false;
    return true;
}
fn addCnp(t: *Totals, c: f64, n: f64, p: f64, plant: bool) void {
    if (plant) {
        t.all_fractions_carbon_g_c += c;
        t.all_fractions_nitrogen_g_n += n;
        t.all_fractions_phosphorus_g_p += p;
    } else {
        t.excluding_fraction_4_carbon_g_c += c;
        t.excluding_fraction_4_nitrogen_g_n += n;
        t.excluding_fraction_4_phosphorus_g_p += p;
    }
}

/// Direct translation of REDIST 12651--12724 for one runtime-indexed layer.
pub fn recalculate(layer: usize, pools: Pools) !Result {
    const layers = pools.layer_count;
    if (layers == 0 or layer >= layers) return error.TillageOrganicLedgerDimensionMismatch;
    inline for (0..3) |i| {
        if (pools.microbial_c_n_p[i].len != 126 * layers or pools.residue_c_n_p[i].len != 10 * layers)
            return error.TillageOrganicLedgerDimensionMismatch;
        if (!finiteNonnegative(pools.microbial_c_n_p[i]) or !finiteNonnegative(pools.residue_c_n_p[i]))
            return error.InvalidTillageOrganicLedgerInput;
    }
    inline for (0..12) |i| {
        if (pools.soluble_and_adsorbed[i].len != 5 * layers) return error.TillageOrganicLedgerDimensionMismatch;
        if (!finiteNonnegative(pools.soluble_and_adsorbed[i])) return error.InvalidTillageOrganicLedgerInput;
    }
    inline for (0..4) |i| {
        if (pools.som_c_a_n_p[i].len != 25 * layers) return error.TillageOrganicLedgerDimensionMismatch;
        if (!finiteNonnegative(pools.som_c_a_n_p[i])) return error.InvalidTillageOrganicLedgerInput;
    }
    var t: Totals = .{};
    for (0..6) |k| for (0..7) |n| for (0..3) |m| {
        const pool = (k * 7 + n) * 3 + m;
        const c = pools.microbial_c_n_p[0][index(layers, pool, layer)];
        const nitrogen = pools.microbial_c_n_p[1][index(layers, pool, layer)];
        const phosphorus = pools.microbial_c_n_p[2][index(layers, pool, layer)];
        addCnp(&t, c, nitrogen, phosphorus, true);
        if (k != 4) addCnp(&t, c, nitrogen, phosphorus, false);
    };
    for (0..5) |k| {
        for (0..2) |m| {
            const pool = k * 2 + m;
            const c = pools.residue_c_n_p[0][index(layers, pool, layer)];
            const n = pools.residue_c_n_p[1][index(layers, pool, layer)];
            const p = pools.residue_c_n_p[2][index(layers, pool, layer)];
            addCnp(&t, c, n, p, true);
            if (k != 4) addCnp(&t, c, n, p, false);
        }
        const c = pools.soluble_and_adsorbed[0][index(layers, k, layer)] + pools.soluble_and_adsorbed[1][index(layers, k, layer)] + pools.soluble_and_adsorbed[2][index(layers, k, layer)] + pools.soluble_and_adsorbed[3][index(layers, k, layer)] + pools.soluble_and_adsorbed[4][index(layers, k, layer)] + pools.soluble_and_adsorbed[5][index(layers, k, layer)];
        const n = pools.soluble_and_adsorbed[6][index(layers, k, layer)] + pools.soluble_and_adsorbed[7][index(layers, k, layer)] + pools.soluble_and_adsorbed[8][index(layers, k, layer)];
        const p = pools.soluble_and_adsorbed[9][index(layers, k, layer)] + pools.soluble_and_adsorbed[10][index(layers, k, layer)] + pools.soluble_and_adsorbed[11][index(layers, k, layer)];
        addCnp(&t, c, n, p, true);
        if (k != 4) addCnp(&t, c, n, p, false);
        for (0..5) |m| {
            const pool = k * 5 + m;
            const som_c = pools.som_c_a_n_p[0][index(layers, pool, layer)];
            const som_n = pools.som_c_a_n_p[2][index(layers, pool, layer)];
            const som_p = pools.som_c_a_n_p[3][index(layers, pool, layer)];
            if (k != 4) {
                if (m < 4) addCnp(&t, som_c, som_n, som_p, false) else {
                    t.excluding_fraction_4_charcoal_carbon_g_c += som_c;
                    t.excluding_fraction_4_charcoal_nitrogen_g_n += som_n;
                    t.excluding_fraction_4_charcoal_phosphorus_g_p = t.excluding_fraction_4_phosphorus_g_p + som_p;
                }
            } else if (m < 4) addCnp(&t, som_c, som_n, som_p, true) else {
                t.fraction_4_charcoal_carbon_g_c += som_c;
                t.fraction_4_charcoal_nitrogen_g_n += som_n;
                t.fraction_4_charcoal_phosphorus_g_p += som_p;
            }
        }
    }
    inline for (std.meta.fields(Totals)) |field| if (!std.math.isFinite(@field(t, field.name))) return error.NonFiniteTillageOrganicLedgerResult;
    return .{ .totals = t, .ledgers = .{ .organic_carbon_g_c = t.excluding_fraction_4_carbon_g_c + t.all_fractions_carbon_g_c, .organic_nitrogen_g_n = t.excluding_fraction_4_nitrogen_g_n + t.all_fractions_nitrogen_g_n, .charcoal_organic_carbon_g_c = t.excluding_fraction_4_charcoal_carbon_g_c + t.fraction_4_charcoal_carbon_g_c, .charcoal_organic_nitrogen_g_n = t.excluding_fraction_4_charcoal_nitrogen_g_n + t.fraction_4_charcoal_nitrogen_g_n, .residue_organic_carbon_g_c = t.excluding_fraction_4_carbon_g_c } };
}

test "REDIST organic ledger preserves K M nesting and literal DPC replacement" {
    var microbial: [3][126]f64 = @splat(@splat(1));
    var residue: [3][10]f64 = @splat(@splat(1));
    var soluble: [12][5]f64 = @splat(@splat(1));
    var som: [4][25]f64 = @splat(@splat(1));
    var microbial_s: [3][]const f64 = undefined;
    var residue_s: [3][]const f64 = undefined;
    var soluble_s: [12][]const f64 = undefined;
    var som_s: [4][]const f64 = undefined;
    inline for (0..3) |i| {
        microbial_s[i] = &microbial[i];
        residue_s[i] = &residue[i];
    }
    inline for (0..12) |i| soluble_s[i] = &soluble[i];
    inline for (0..4) |i| som_s[i] = &som[i];
    const result = try recalculate(0, .{ .layer_count = 1, .microbial_c_n_p = microbial_s, .residue_c_n_p = residue_s, .soluble_and_adsorbed = soluble_s, .som_c_a_n_p = som_s });
    try std.testing.expectEqual(@as(f64, 153), result.totals.excluding_fraction_4_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 170), result.totals.all_fractions_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 142), result.totals.excluding_fraction_4_charcoal_phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 323), result.ledgers.organic_carbon_g_c);
}
