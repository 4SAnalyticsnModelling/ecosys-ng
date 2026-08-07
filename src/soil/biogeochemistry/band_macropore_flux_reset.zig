const std = @import("std");

pub const LayerPosition = enum { surface, subsurface };

pub const HydraulicFluxes = struct {
    water_m3_timestep: f64, // FLW
    vapor_m3_timestep: f64, // FLV
    water_x_m3_timestep: f64, // FLWX
    water_horizontal_m3_timestep: f64, // FLWH
    water_y_m3_timestep: f64, // FLWY
    water_horizontal_y_m3_timestep: f64, // FLWHY
    heat_j_timestep: f64, // HFLW
};

pub const BandSoluteFluxes = struct {
    ammonium_g_n_timestep: f64, // XN4FLB
    ammonia_g_n_timestep: f64, // XN3FLB
    nitrate_g_n_timestep: f64, // XNOFLB
    other_nitrogen_g_n_timestep: f64, // XNXFLB
    phosphate_h2po4_g_p_timestep: f64, // XH1BFB
    phosphate_hpo4_g_p_timestep: f64, // XH2BFB
};

pub const MacroporeFluxes = struct {
    carbon_dioxide_g_timestep: f64, // XCOFHS
    methane_g_timestep: f64, // XCHFHS
    oxygen_g_timestep: f64, // XOXFHS
    nitrogen_g_timestep: f64, // XNGFHS
    nitrous_oxide_g_timestep: f64, // XN2FHS
    hydrogen_g_timestep: f64, // XHGFHS
    ammonium_g_n_timestep: f64, // XN4FHW
    ammonia_g_n_timestep: f64, // XN3FHW
    nitrate_g_n_timestep: f64, // XNOFHW
    other_nitrogen_g_n_timestep: f64, // XNXFHS
    phosphate_h2po4_g_p_timestep: f64, // XH1PHS
    phosphate_hpo4_g_p_timestep: f64, // XH2PHS
    band_ammonium_g_n_timestep: f64, // XN4FHB
    band_ammonia_g_n_timestep: f64, // XN3FHB
    band_nitrate_g_n_timestep: f64, // XNOFHB
    band_other_nitrogen_g_n_timestep: f64, // XNXFHB
    band_phosphate_h2po4_g_p_timestep: f64, // XH1BHB
    band_phosphate_hpo4_g_p_timestep: f64, // XH2BHB
};

pub const OrganicFluxes = struct {
    carbon_g_timestep: []f64, // XOCFHS
    nitrogen_g_timestep: []f64, // XONFHS
    phosphorus_g_timestep: []f64, // XOPFHS
    acetate_g_timestep: []f64, // XOAFHS
};

pub const Accumulators = struct {
    hydraulic: HydraulicFluxes,
    band_solutes: BandSoluteFluxes,
    macropore: MacroporeFluxes,
    organic: OrganicFluxes,
};

pub const ResetError = error{OrganicPoolLengthMismatch};

/// Translates `hour1.f` lines 2649--2689. Surface-layer calls are a no-op,
/// matching the legacy `L.NE.0` guard.
pub fn reset(
    layer_position: LayerPosition,
    organic_pool_count: usize,
    accumulators: *Accumulators,
) ResetError!void {
    if (layer_position == .surface) return;
    if (accumulators.organic.carbon_g_timestep.len != organic_pool_count or
        accumulators.organic.nitrogen_g_timestep.len != organic_pool_count or
        accumulators.organic.phosphorus_g_timestep.len != organic_pool_count or
        accumulators.organic.acetate_g_timestep.len != organic_pool_count)
    {
        return error.OrganicPoolLengthMismatch;
    }

    inline for (std.meta.fields(HydraulicFluxes)) |field| {
        @field(accumulators.hydraulic, field.name) = 0.0;
    }
    inline for (std.meta.fields(BandSoluteFluxes)) |field| {
        @field(accumulators.band_solutes, field.name) = 0.0;
    }
    inline for (std.meta.fields(MacroporeFluxes)) |field| {
        @field(accumulators.macropore, field.name) = 0.0;
    }
    for (0..organic_pool_count) |pool| {
        accumulators.organic.carbon_g_timestep[pool] = 0.0;
        accumulators.organic.nitrogen_g_timestep[pool] = 0.0;
        accumulators.organic.phosphorus_g_timestep[pool] = 0.0;
        accumulators.organic.acetate_g_timestep[pool] = 0.0;
    }
}

test "subsurface reset clears every flux and runtime organic pool" {
    var carbon = [_]f64{ 4.0, 5.0, 6.0 };
    var nitrogen = carbon;
    var phosphorus = carbon;
    var acetate = carbon;
    var accumulators = Accumulators{
        .hydraulic = undefined,
        .band_solutes = undefined,
        .macropore = undefined,
        .organic = .{
            .carbon_g_timestep = &carbon,
            .nitrogen_g_timestep = &nitrogen,
            .phosphorus_g_timestep = &phosphorus,
            .acetate_g_timestep = &acetate,
        },
    };
    inline for (std.meta.fields(HydraulicFluxes)) |field| {
        @field(accumulators.hydraulic, field.name) = 8.0;
    }
    inline for (std.meta.fields(BandSoluteFluxes)) |field| {
        @field(accumulators.band_solutes, field.name) = 8.0;
    }
    inline for (std.meta.fields(MacroporeFluxes)) |field| {
        @field(accumulators.macropore, field.name) = 8.0;
    }

    try reset(.subsurface, 3, &accumulators);

    inline for (std.meta.fields(HydraulicFluxes)) |field| {
        try std.testing.expectEqual(@as(f64, 0.0), @field(accumulators.hydraulic, field.name));
    }
    inline for (std.meta.fields(BandSoluteFluxes)) |field| {
        try std.testing.expectEqual(@as(f64, 0.0), @field(accumulators.band_solutes, field.name));
    }
    inline for (std.meta.fields(MacroporeFluxes)) |field| {
        try std.testing.expectEqual(@as(f64, 0.0), @field(accumulators.macropore, field.name));
    }
    inline for (std.meta.fields(OrganicFluxes)) |field| {
        for (@field(accumulators.organic, field.name)) |value| {
            try std.testing.expectEqual(@as(f64, 0.0), value);
        }
    }
}

test "surface layer preserves flux state" {
    var accumulators: Accumulators = undefined;
    accumulators.hydraulic.water_m3_timestep = 12.0;
    try reset(.surface, 99, &accumulators);
    try std.testing.expectEqual(@as(f64, 12.0), accumulators.hydraulic.water_m3_timestep);
}
