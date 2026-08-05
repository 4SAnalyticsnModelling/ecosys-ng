const std = @import("std");
const organic = @import("redist_erosion_organic_matter_apply.zig");

pub const DissolvedPools = struct {
    doc_matrix_g: [3]f64,
    doc_macropore_g: [3]f64,
    adsorbed_carbon_g: [3]f64,
    acetate_matrix_g: [3]f64,
    acetate_macropore_g: [3]f64,
    adsorbed_acetate_g: [3]f64,
    don_matrix_g: [3]f64,
    don_macropore_g: [3]f64,
    adsorbed_nitrogen_g: [3]f64,
    dop_matrix_g: [3]f64,
    dop_macropore_g: [3]f64,
    adsorbed_phosphorus_g: [3]f64,
};

pub const MicrobialDiagnostics = struct {
    carbon_g: f64 = 0.0,
    nitrogen_g: f64 = 0.0,
    phosphorus_g: f64 = 0.0,
};

pub const Result = struct {
    organic_carbon_g: f64,
    organic_nitrogen_g: f64,
    organic_phosphorus_g: f64,
    colonized_carbon_g: f64,
    colonized_nitrogen_g: f64,
    colonized_phosphorus_g: f64,
    carbon_by_class_g: [6]f64,
    microbial_carbon_layer_g: f64,
    microbial_nitrogen_layer_g: f64,
    microbial_diagnostics: MicrobialDiagnostics,
};

fn addChecked(total: *f64, value: f64) !void {
    if (!std.math.isFinite(value)) return error.InvalidSurfaceLitterOrganicPool;
    total.* += value;
    if (!std.math.isFinite(total.*)) return error.NonFiniteSurfaceLitterOrganicInventory;
}

/// Exact standalone translation of REDIST lines 5522--5630.
pub fn calculate(
    microbial: organic.MicrobialPools,
    residue: organic.ResiduePools,
    dissolved: DissolvedPools,
    soc: organic.SocPools,
    initial_diagnostics: MicrobialDiagnostics,
) !Result {
    var dc: f64 = 0.0;
    var dn: f64 = 0.0;
    var dp: f64 = 0.0;
    var dcc: f64 = 0.0;
    var dnc: f64 = 0.0;
    var dpc: f64 = 0.0;
    var rc0 = [_]f64{0.0} ** 6;
    var omcl: f64 = 0.0;
    var omnl: f64 = 0.0;
    var diagnostics = initial_diagnostics;
    inline for (@typeInfo(MicrobialDiagnostics).@"struct".fields) |field|
        if (!std.math.isFinite(@field(diagnostics, field.name))) return error.InvalidSurfaceLitterOrganicPool;

    for (0..6) |k| {
        if (k != 4) {
            for (0..7) |n| for (0..3) |m| {
                const c = microbial.omc[k][n][m];
                const nitrogen = microbial.omn[k][n][m];
                const phosphorus = microbial.omp[k][n][m];
                try addChecked(&dc, c);
                try addChecked(&dn, nitrogen);
                try addChecked(&dp, phosphorus);
                try addChecked(&rc0[k], c);
                try addChecked(&diagnostics.carbon_g, c);
                try addChecked(&diagnostics.nitrogen_g, nitrogen);
                try addChecked(&diagnostics.phosphorus_g, phosphorus);
                try addChecked(&omcl, c);
                try addChecked(&omnl, nitrogen);
            };
        }
    }
    for (0..3) |k| {
        for (0..2) |m| {
            try addChecked(&dc, residue.orc[k][m]);
            try addChecked(&dn, residue.orn[k][m]);
            try addChecked(&dp, residue.orp[k][m]);
            try addChecked(&rc0[k], residue.orc[k][m]);
        }
        const carbon_terms = .{
            dissolved.doc_matrix_g[k],     dissolved.doc_macropore_g[k],     dissolved.adsorbed_carbon_g[k],
            dissolved.acetate_matrix_g[k], dissolved.acetate_macropore_g[k], dissolved.adsorbed_acetate_g[k],
        };
        inline for (carbon_terms) |value| {
            try addChecked(&dc, value);
            try addChecked(&rc0[k], value);
        }
        inline for (.{ dissolved.don_matrix_g[k], dissolved.don_macropore_g[k], dissolved.adsorbed_nitrogen_g[k] }) |value|
            try addChecked(&dn, value);
        inline for (.{ dissolved.dop_matrix_g[k], dissolved.dop_macropore_g[k], dissolved.adsorbed_phosphorus_g[k] }) |value|
            try addChecked(&dp, value);
    }
    for (0..5) |k| for (0..5) |m| {
        const c = soc.osc[k][m];
        const nitrogen = soc.osn[k][m];
        const phosphorus = soc.osp[k][m];
        if (m <= 3) {
            try addChecked(&dc, c);
            try addChecked(&dn, nitrogen);
            try addChecked(&dp, phosphorus);
        } else {
            try addChecked(&dcc, c);
            try addChecked(&dnc, nitrogen);
            try addChecked(&dpc, phosphorus);
        }
        try addChecked(&rc0[k], c);
    };
    return .{
        .organic_carbon_g = dc,
        .organic_nitrogen_g = dn,
        .organic_phosphorus_g = dp,
        .colonized_carbon_g = dcc,
        .colonized_nitrogen_g = dnc,
        .colonized_phosphorus_g = dpc,
        .carbon_by_class_g = rc0,
        .microbial_carbon_layer_g = omcl,
        .microbial_nitrogen_layer_g = omnl,
        .microbial_diagnostics = diagnostics,
    };
}

test "REDIST litter organic inventory preserves class and loop boundaries" {
    var microbial = std.mem.zeroes(organic.MicrobialPools);
    microbial.omc[0][0][0] = 1.0;
    microbial.omn[5][6][2] = 2.0;
    microbial.omp[4][0][0] = 99.0; // K=4 excluded from microbial totals.
    var residue = std.mem.zeroes(organic.ResiduePools);
    residue.orc[2][1] = 3.0;
    var soc = std.mem.zeroes(organic.SocPools);
    soc.osc[4][3] = 4.0;
    soc.osc[4][4] = 5.0;
    soc.osn[4][4] = 6.0;
    soc.osp[4][4] = 7.0;
    const result = try calculate(microbial, residue, std.mem.zeroes(DissolvedPools), soc, .{});
    try std.testing.expectEqual(@as(f64, 8.0), result.organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 5.0), result.colonized_carbon_g);
    try std.testing.expectEqual(@as(f64, 9.0), result.carbon_by_class_g[4]);
    try std.testing.expectEqual(@as(f64, 2.0), result.organic_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 6.0), result.colonized_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 7.0), result.colonized_phosphorus_g);
}

test "REDIST litter organic inventory includes all dissolved carriers" {
    var dissolved = std.mem.zeroes(DissolvedPools);
    inline for (@typeInfo(DissolvedPools).@"struct".fields) |field| @field(dissolved, field.name)[1] = 1.0;
    const result = try calculate(
        std.mem.zeroes(organic.MicrobialPools),
        std.mem.zeroes(organic.ResiduePools),
        dissolved,
        std.mem.zeroes(organic.SocPools),
        .{},
    );
    try std.testing.expectEqual(@as(f64, 6.0), result.organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 3.0), result.organic_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 3.0), result.organic_phosphorus_g);
    try std.testing.expectEqual(@as(f64, 6.0), result.carbon_by_class_g[1]);
}

test "REDIST litter organic inventory rejects invalid and overflow" {
    var soc = std.mem.zeroes(organic.SocPools);
    soc.osn[0][0] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSurfaceLitterOrganicPool, calculate(
        std.mem.zeroes(organic.MicrobialPools),
        std.mem.zeroes(organic.ResiduePools),
        std.mem.zeroes(DissolvedPools),
        soc,
        .{},
    ));
    var microbial = std.mem.zeroes(organic.MicrobialPools);
    microbial.omc[0][0][0] = std.math.floatMax(f64);
    microbial.omc[0][0][1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSurfaceLitterOrganicInventory, calculate(
        microbial,
        std.mem.zeroes(organic.ResiduePools),
        std.mem.zeroes(DissolvedPools),
        std.mem.zeroes(organic.SocPools),
        .{},
    ));
}
