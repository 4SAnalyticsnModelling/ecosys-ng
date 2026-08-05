const std = @import("std");

pub const LitterState = struct {
    /// VHCP(0). Volumetric heat capacity (MJ K-1).
    heat_capacity_megajoules_per_k: f64,
    /// TKS(0). Temperature (K).
    temperature_k: f64,
    /// VOLW(0). Liquid water (m3).
    liquid_m3: f64,
    /// VOLV(0). Vapor (m3).
    vapor_m3: f64,
    /// VOLI(0). Ice (m3).
    ice_m3: f64,
    /// ORGC(0). Organic carbon (g C).
    organic_carbon_g: f64,
    /// ORGCC(0). Charcoal (g C).
    charcoal_g: f64,
};

/// Water fluxes from watsub.f to surface litter (m3 step-1).
pub const LitterWaterFluxes = struct {
    /// FLWR. Liquid water input.
    liquid_in_m3: f64,
    /// FLVR. Vapor input.
    vapor_in_m3: f64,
    /// XWFLVR. Evaporation-condensation (positive = condensation).
    evap_condensation_m3: f64,
    /// XWFLFR. Freeze-thaw (positive = freeze).
    freeze_thaw_m3: f64,
    /// TQR. Runoff.
    runoff_m3: f64,
};

/// Heat fluxes from watsub.f to surface litter (MJ step-1).
pub const LitterHeatFluxes = struct {
    /// HFLWR. Liquid water heat input.
    liquid_megajoules: f64,
    /// XHFLFR. Freeze-thaw heat flux.
    freeze_thaw_megajoules: f64,
    /// XHFLVR. Evaporation-condensation heat flux.
    evap_condensation_megajoules: f64,
    /// THQR. Runoff heat flux.
    runoff_megajoules: f64,
};

pub const Parameters = struct {
    /// DENSI. Ice density (Mg m-3); converts freeze-thaw water to ice volume.
    ice_density_megagrams_per_m3: f64,
    /// VHCPRX. Minimum litter heat capacity for temperature update (MJ K-1).
    min_heat_capacity_megajoules_per_k: f64,
    /// TKS(NUM). Mineral layer temperature fallback (K).
    mineral_layer_temperature_k: f64,
};

pub const HeatDiagnostics = struct {
    /// HEATIN. Cumulative surface heat exchange before this block (MJ).
    cumulative_surface_heat_megajoules: f64,
    /// HEATSO. Total heat content before this block (MJ).
    total_heat_content_megajoules: f64,
    /// TSMX(0). Running maximum litter temperature (°C).
    max_temperature_c: f64,
    /// TSMN(0). Running minimum litter temperature (°C).
    min_temperature_c: f64,
};

pub const Inputs = struct {
    litter: LitterState,
    water_fluxes: LitterWaterFluxes,
    heat_fluxes: LitterHeatFluxes,
    /// HCBFX(0). Heat released by surface litter combustion (MJ step-1).
    combustion_heat_megajoules: f64,
    parameters: Parameters,
    diagnostics: HeatDiagnostics,
};

pub const UpdatedLitter = struct {
    liquid_m3: f64,
    vapor_m3: f64,
    ice_m3: f64,
    heat_capacity_megajoules_per_k: f64,
    temperature_k: f64,
    /// TCS(0). Celsius temperature.
    temperature_c: f64,
};

pub const UpdatedDiagnostics = struct {
    cumulative_surface_heat_megajoules: f64,
    total_heat_content_megajoules: f64,
    max_temperature_c: f64,
    min_temperature_c: f64,
};

pub const Result = struct {
    litter: UpdatedLitter,
    diagnostics: UpdatedDiagnostics,
};

// Heat capacity coefficients (MJ m-3 K-1, or MJ g-1 K-1 for organic).
const liquid_water_heat_capacity: f64 = 4.19;
const ice_heat_capacity: f64 = 1.9274;
const organic_carbon_heat_capacity: f64 = 2.496e-6;

/// Direct translation of REDIST lines 4318--4347.
///
/// Computes the surface litter temperature update from watsub.f water and
/// heat fluxes. Source-order preserved:
///   1. Save VHCPZ=VHCP(0), compute capacity from current volumes (VHCPY),
///      derive heat-capacity-change flux HFLXO = (VHCPY-VHCPZ)*TKS(0).
///   2. Apply water fluxes to VOLW(0), VOLV(0), VOLI(0).
///   3. Save old energy ENGYZ = VHCPZ*TKS(0).
///   4. Recompute VHCP(0) from updated volumes.
///   5. Update TKS(0) via energy conservation or fallback to TKS(NUM).
///   6. Accumulate HEATIN and HEATSO diagnostics.
///   7. Derive TCS(0), update TSMX/TSMN.
pub fn update(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const p = inputs.parameters;
    const litter = inputs.litter;
    const wf = inputs.water_fluxes;
    const hf = inputs.heat_fluxes;

    // Lines 4318--4319: save current capacity.
    const vhcpz = litter.heat_capacity_megajoules_per_k;

    // Lines 4320--4324: capacity from current volumes, then capacity-change flux.
    const vhcpy = organic_carbon_heat_capacity * (litter.organic_carbon_g + litter.charcoal_g) +
        liquid_water_heat_capacity * (litter.liquid_m3 + litter.vapor_m3) +
        ice_heat_capacity * litter.ice_m3;
    const hflxo = (vhcpy - vhcpz) * litter.temperature_k;

    // Lines 4325--4328: apply water fluxes.
    const volw = litter.liquid_m3 + wf.liquid_in_m3 + wf.evap_condensation_m3 +
        wf.freeze_thaw_m3 + wf.runoff_m3;
    const volv = litter.vapor_m3 + wf.vapor_in_m3 - wf.evap_condensation_m3;
    const voli = litter.ice_m3 - wf.freeze_thaw_m3 / p.ice_density_megagrams_per_m3;
    inline for (.{ volw, volv, voli }) |volume_m3| {
        if (!std.math.isFinite(volume_m3))
            return error.NonFiniteLitterHeatStorageResult;
        if (volume_m3 < 0)
            return error.NegativeLitterWaterPhaseInventory;
    }

    // Line 4329: old energy using saved vhcpz and current (old) temperature.
    const engyz = vhcpz * litter.temperature_k;

    // Lines 4330--4332: new heat capacity from updated volumes.
    const vhcp_new = organic_carbon_heat_capacity * (litter.organic_carbon_g + litter.charcoal_g) +
        liquid_water_heat_capacity * (volw + volv) +
        ice_heat_capacity * voli;

    // Lines 4333--4341: temperature update with energy conservation or fallback.
    var tks_new: f64 = undefined;
    var heatin_after_temperature_update = inputs.diagnostics.cumulative_surface_heat_megajoules +
        hflxo + inputs.combustion_heat_megajoules;
    if (vhcp_new > p.min_heat_capacity_megajoules_per_k) {
        // Lines 4334--4336. Keep the source expression tree; each addition is
        // intentionally written in the same order as the fixed-form statement.
        tks_new = (engyz + hf.liquid_megajoules + hf.freeze_thaw_megajoules +
            hf.evap_condensation_megajoules + hf.runoff_megajoules + hflxo +
            inputs.combustion_heat_megajoules) / vhcp_new;
    } else {
        // Lines 4338-4340.
        heatin_after_temperature_update = heatin_after_temperature_update +
            (p.mineral_layer_temperature_k - litter.temperature_k) * vhcp_new;
        tks_new = p.mineral_layer_temperature_k;
    }

    // Lines 4342--4344: energy content and heat diagnostics.
    const engyr = vhcp_new * tks_new;
    const heatso_new = inputs.diagnostics.total_heat_content_megajoules + engyr;
    // Line 4344: add evaporation and freeze-thaw heat terms.
    const heatin_new = heatin_after_temperature_update +
        hf.freeze_thaw_megajoules + hf.evap_condensation_megajoules;

    // Lines 4345--4347: Celsius temperature and running extrema.
    const tcs_new = tks_new - 273.15;
    const tsmx_new = @max(inputs.diagnostics.max_temperature_c, tcs_new);
    const tsmn_new = @min(inputs.diagnostics.min_temperature_c, tcs_new);

    const result = Result{
        .litter = .{
            .liquid_m3 = volw,
            .vapor_m3 = volv,
            .ice_m3 = voli,
            .heat_capacity_megajoules_per_k = vhcp_new,
            .temperature_k = tks_new,
            .temperature_c = tcs_new,
        },
        .diagnostics = .{
            .cumulative_surface_heat_megajoules = heatin_new,
            .total_heat_content_megajoules = heatso_new,
            .max_temperature_c = tsmx_new,
            .min_temperature_c = tsmn_new,
        },
    };
    inline for (@typeInfo(UpdatedLitter).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.litter, field.name)))
            return error.NonFiniteLitterHeatStorageResult;
    inline for (@typeInfo(UpdatedDiagnostics).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.diagnostics, field.name)))
            return error.NonFiniteLitterHeatStorageResult;
    if (result.litter.heat_capacity_megajoules_per_k < 0 or
        result.litter.temperature_k <= 0)
    {
        return error.InvalidLitterHeatStorageResult;
    }
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(LitterState).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.litter, field.name)))
            return error.NonFiniteLitterHeatStorageInput;
    }
    inline for (@typeInfo(LitterWaterFluxes).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.water_fluxes, field.name)))
            return error.NonFiniteLitterHeatStorageInput;
    }
    inline for (@typeInfo(LitterHeatFluxes).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.heat_fluxes, field.name)))
            return error.NonFiniteLitterHeatStorageInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.parameters, field.name)))
            return error.NonFiniteLitterHeatStorageInput;
    }
    inline for (@typeInfo(HeatDiagnostics).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.diagnostics, field.name)))
            return error.NonFiniteLitterHeatStorageInput;
    }
    if (!std.math.isFinite(inputs.combustion_heat_megajoules))
        return error.NonFiniteLitterHeatStorageInput;

    const litter = inputs.litter;
    if (litter.heat_capacity_megajoules_per_k < 0 or litter.temperature_k <= 0 or
        litter.liquid_m3 < 0 or litter.vapor_m3 < 0 or litter.ice_m3 < 0 or
        litter.organic_carbon_g < 0 or litter.charcoal_g < 0)
    {
        return error.InvalidLitterHeatStorageState;
    }
    const p = inputs.parameters;
    if (p.ice_density_megagrams_per_m3 <= 0 or
        p.min_heat_capacity_megajoules_per_k < 0 or
        p.mineral_layer_temperature_k <= 0)
    {
        return error.InvalidLitterHeatStorageParameter;
    }
}

fn defaultInputs() Inputs {
    const litter = LitterState{
        .heat_capacity_megajoules_per_k = 0.05,
        .temperature_k = 278.0,
        .liquid_m3 = 0.01,
        .vapor_m3 = 0.002,
        .ice_m3 = 0.0,
        .organic_carbon_g = 800.0,
        .charcoal_g = 20.0,
    };
    return .{
        .litter = litter,
        .water_fluxes = .{
            .liquid_in_m3 = 0.005,
            .vapor_in_m3 = 0.001,
            .evap_condensation_m3 = 0.0,
            .freeze_thaw_m3 = 0.0,
            .runoff_m3 = 0.0,
        },
        .heat_fluxes = .{
            .liquid_megajoules = 0.01,
            .freeze_thaw_megajoules = 0.0,
            .evap_condensation_megajoules = 0.0,
            .runoff_megajoules = 0.0,
        },
        .combustion_heat_megajoules = 0.0,
        .parameters = .{
            .ice_density_megagrams_per_m3 = 0.917,
            .min_heat_capacity_megajoules_per_k = 1.0e-6,
            .mineral_layer_temperature_k = 280.0,
        },
        .diagnostics = .{
            .cumulative_surface_heat_megajoules = 0.0,
            .total_heat_content_megajoules = 0.0,
            .max_temperature_c = -50.0,
            .min_temperature_c = 50.0,
        },
    };
}

test "REDIST litter heat storage preserves source-order water flux application" {
    const inp = defaultInputs();
    const result = try update(inp);

    const wf = inp.water_fluxes;
    try std.testing.expectApproxEqRel(
        inp.litter.liquid_m3 + wf.liquid_in_m3 + wf.evap_condensation_m3 +
            wf.freeze_thaw_m3 + wf.runoff_m3,
        result.litter.liquid_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        inp.litter.vapor_m3 + wf.vapor_in_m3 - wf.evap_condensation_m3,
        result.litter.vapor_m3,
        1.0e-15,
    );
}

test "REDIST litter heat storage applies energy conservation when capacity sufficient" {
    const inp = defaultInputs();
    const result = try update(inp);

    const vhcpy = organic_carbon_heat_capacity * (inp.litter.organic_carbon_g + inp.litter.charcoal_g) +
        liquid_water_heat_capacity * (inp.litter.liquid_m3 + inp.litter.vapor_m3) +
        ice_heat_capacity * inp.litter.ice_m3;
    const hflxo = (vhcpy - inp.litter.heat_capacity_megajoules_per_k) * inp.litter.temperature_k;
    const engyz = inp.litter.heat_capacity_megajoules_per_k * inp.litter.temperature_k;
    const volw_new = inp.litter.liquid_m3 + inp.water_fluxes.liquid_in_m3;
    const volv_new = inp.litter.vapor_m3 + inp.water_fluxes.vapor_in_m3;
    const vhcp_new = organic_carbon_heat_capacity * (inp.litter.organic_carbon_g + inp.litter.charcoal_g) +
        liquid_water_heat_capacity * (volw_new + volv_new) +
        ice_heat_capacity * inp.litter.ice_m3;
    const heat_sum = inp.heat_fluxes.liquid_megajoules;
    const expected_tks = (engyz + heat_sum + hflxo + inp.combustion_heat_megajoules) / vhcp_new;
    try std.testing.expectApproxEqRel(
        expected_tks,
        result.litter.temperature_k,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        expected_tks - 273.15,
        result.litter.temperature_c,
        1.0e-14,
    );
}

test "REDIST litter heat storage falls back to mineral temperature when capacity below threshold" {
    var inp = defaultInputs();
    inp.parameters.min_heat_capacity_megajoules_per_k = 1.0;
    const result = try update(inp);
    try std.testing.expectEqual(
        inp.parameters.mineral_layer_temperature_k,
        result.litter.temperature_k,
    );
}

test "REDIST litter heat storage updates extrema diagnostics" {
    const inp = defaultInputs();
    const result = try update(inp);
    // Extrema follow AMAX1/AMIN1 of prior value and new temperature.
    try std.testing.expectEqual(
        @max(inp.diagnostics.max_temperature_c, result.litter.temperature_c),
        result.diagnostics.max_temperature_c,
    );
    try std.testing.expectEqual(
        @min(inp.diagnostics.min_temperature_c, result.litter.temperature_c),
        result.diagnostics.min_temperature_c,
    );
}

test "REDIST litter heat storage rejects a negative post-flux phase inventory" {
    var inp = defaultInputs();
    inp.water_fluxes.evap_condensation_m3 = inp.litter.vapor_m3 +
        inp.water_fluxes.vapor_in_m3 + 1.0e-6;
    try std.testing.expectError(
        error.NegativeLitterWaterPhaseInventory,
        update(inp),
    );
}

test "REDIST litter heat storage rejects non-finite diagnostics" {
    var inp = defaultInputs();
    inp.diagnostics.total_heat_content_megajoules = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteLitterHeatStorageInput,
        update(inp),
    );
}
