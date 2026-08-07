const std = @import("std");

/// Per-layer, per-direction flux accumulators. Element names identify the
/// transported species; values retain legacy mass-per-timestep units.
pub const FluxAccumulators = struct {
    carbon_dioxide_g_timestep: f64, // XCOFLG
    methane_g_timestep: f64, // XCHFLG
    oxygen_g_timestep: f64, // XOXFLG
    nitrogen_g_timestep: f64, // XNGFLG
    nitrous_oxide_g_timestep: f64, // XN2FLG
    ammonia_g_timestep: f64, // XN3FLG
    hydrogen_g_timestep: f64, // XHGFLG
    dissolved_carbon_dioxide_g_timestep: f64, // XCOFLS
    dissolved_methane_g_timestep: f64, // XCHFLS
    dissolved_oxygen_g_timestep: f64, // XOXFLS
    dissolved_nitrogen_g_timestep: f64, // XNGFLS
    dissolved_nitrous_oxide_g_timestep: f64, // XN2FLS
    dissolved_hydrogen_g_timestep: f64, // XHGFLS
    ammonium_g_n_timestep: f64, // XN4FLW
    ammonia_g_n_timestep: f64, // XN3FLW
    nitrate_g_n_timestep: f64, // XNOFLW
    other_nitrogen_g_n_timestep: f64, // XNXFLS
    phosphate_h2po4_g_p_timestep: f64, // XH1PFS
    phosphate_hpo4_g_p_timestep: f64, // XH2PFS
    organic_carbon_g_timestep: []f64, // XOCFLS
    organic_nitrogen_g_timestep: []f64, // XONFLS
    organic_phosphorus_g_timestep: []f64, // XOPFLS
    organic_acetate_g_timestep: []f64, // XOAFLS
};

pub const ResetError = error{OrganicPoolLengthMismatch};

/// Translates `hour1.f` lines 2620--2644 in source assignment order.
pub fn reset(
    accumulators: *FluxAccumulators,
    organic_pool_count: usize,
) ResetError!void {
    if (accumulators.organic_carbon_g_timestep.len != organic_pool_count or
        accumulators.organic_nitrogen_g_timestep.len != organic_pool_count or
        accumulators.organic_phosphorus_g_timestep.len != organic_pool_count or
        accumulators.organic_acetate_g_timestep.len != organic_pool_count)
    {
        return error.OrganicPoolLengthMismatch;
    }

    accumulators.carbon_dioxide_g_timestep = 0.0;
    accumulators.methane_g_timestep = 0.0;
    accumulators.oxygen_g_timestep = 0.0;
    accumulators.nitrogen_g_timestep = 0.0;
    accumulators.nitrous_oxide_g_timestep = 0.0;
    accumulators.ammonia_g_timestep = 0.0;
    accumulators.hydrogen_g_timestep = 0.0;
    accumulators.dissolved_carbon_dioxide_g_timestep = 0.0;
    accumulators.dissolved_methane_g_timestep = 0.0;
    accumulators.dissolved_oxygen_g_timestep = 0.0;
    accumulators.dissolved_nitrogen_g_timestep = 0.0;
    accumulators.dissolved_nitrous_oxide_g_timestep = 0.0;
    accumulators.dissolved_hydrogen_g_timestep = 0.0;
    accumulators.ammonium_g_n_timestep = 0.0;
    accumulators.ammonia_g_n_timestep = 0.0;
    accumulators.nitrate_g_n_timestep = 0.0;
    accumulators.other_nitrogen_g_n_timestep = 0.0;
    accumulators.phosphate_h2po4_g_p_timestep = 0.0;
    accumulators.phosphate_hpo4_g_p_timestep = 0.0;
    for (0..organic_pool_count) |pool| {
        accumulators.organic_carbon_g_timestep[pool] = 0.0;
        accumulators.organic_nitrogen_g_timestep[pool] = 0.0;
        accumulators.organic_phosphorus_g_timestep[pool] = 0.0;
        accumulators.organic_acetate_g_timestep[pool] = 0.0;
    }
}

test "layer gas and solute reset covers runtime organic pools" {
    var organic_carbon = [_]f64{ 3.0, 4.0, 5.0 };
    var organic_nitrogen = [_]f64{ 3.0, 4.0, 5.0 };
    var organic_phosphorus = [_]f64{ 3.0, 4.0, 5.0 };
    var organic_acetate = [_]f64{ 3.0, 4.0, 5.0 };
    var accumulators: FluxAccumulators = undefined;
    inline for (std.meta.fields(FluxAccumulators)) |field| {
        if (field.type == f64) @field(accumulators, field.name) = 8.0;
    }
    accumulators.organic_carbon_g_timestep = &organic_carbon;
    accumulators.organic_nitrogen_g_timestep = &organic_nitrogen;
    accumulators.organic_phosphorus_g_timestep = &organic_phosphorus;
    accumulators.organic_acetate_g_timestep = &organic_acetate;

    try reset(&accumulators, 3);

    inline for (std.meta.fields(FluxAccumulators)) |field| {
        if (field.type == f64) {
            try std.testing.expectEqual(@as(f64, 0.0), @field(accumulators, field.name));
        } else {
            for (@field(accumulators, field.name)) |value| {
                try std.testing.expectEqual(@as(f64, 0.0), value);
            }
        }
    }
}

test "organic pool length mismatch fails before mutation" {
    var short = [_]f64{1.0};
    var accumulators: FluxAccumulators = undefined;
    accumulators.carbon_dioxide_g_timestep = 9.0;
    accumulators.organic_carbon_g_timestep = &short;
    accumulators.organic_nitrogen_g_timestep = &short;
    accumulators.organic_phosphorus_g_timestep = &short;
    accumulators.organic_acetate_g_timestep = &short;
    try std.testing.expectError(error.OrganicPoolLengthMismatch, reset(&accumulators, 2));
    try std.testing.expectEqual(@as(f64, 9.0), accumulators.carbon_dioxide_g_timestep);
}
