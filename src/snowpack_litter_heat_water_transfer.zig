const std = @import("std");

pub const SnowLayer = struct {
    /// VOLSSL(1). Solid snow volume (m3).
    solid_m3: f64,
    /// VOLWSL(1). Liquid water volume (m3).
    liquid_m3: f64,
    /// VOLVSL(1). Vapor volume (m3).
    vapor_m3: f64,
    /// VOLISL(1). Ice volume (m3).
    ice_m3: f64,
    /// TKW(1). Temperature (K).
    temperature_k: f64,
    /// VHCPW(1). Volumetric heat capacity (MJ K-1).
    heat_capacity_megajoules_per_k: f64,
};

pub const SurfaceLitter = struct {
    /// VOLW(0). Liquid water volume (m3).
    liquid_m3: f64,
    /// VOLV(0). Vapor volume (m3).
    vapor_m3: f64,
    /// VOLI(0). Ice volume (m3).
    ice_m3: f64,
    /// TKS(0). Temperature (K).
    temperature_k: f64,
    /// VHCP(0). Volumetric heat capacity (MJ K-1).
    heat_capacity_megajoules_per_k: f64,
};

pub const Fluxes = struct {
    /// XFLWSX. Solid snow flux from watsub.f to litter (m3 step-1).
    solid_m3: f64,
    /// XFLWWX. Liquid water flux (m3 step-1).
    liquid_m3: f64,
    /// XFLWVX. Vapor flux (m3 step-1).
    vapor_m3: f64,
    /// XFLWIX. Ice flux (m3 step-1).
    ice_m3: f64,
    /// XHFLWX. Heat flux from watsub.f to litter (MJ step-1). Must be > 0.
    heat_megajoules: f64,
};

pub const Parameters = struct {
    /// DENSI. Ice density (Mg m-3); converts solid snow volume to ice volume.
    ice_density_megagrams_per_m3: f64,
    /// DENS0. Initial snow density applied to all layers after full melt
    /// (Mg m-3). Also used to compute VOLS from VOLSSL.
    reset_snow_density_megagrams_per_m3: f64,
    /// VHCPWX. Minimum snow layer heat capacity for temperature update (MJ K-1).
    min_snow_heat_capacity_megajoules_per_k: f64,
    /// TKQ. Ambient temperature fallback for snow layer when heat capacity
    /// falls below threshold (K).
    ambient_temperature_k: f64,
    /// VHCPRX. Minimum litter heat capacity for temperature update (MJ K-1).
    min_litter_heat_capacity_megajoules_per_k: f64,
    /// TKS(NUM). Mineral layer temperature fallback for litter when litter
    /// heat capacity falls below threshold (K).
    mineral_layer_temperature_k: f64,
};

pub const LitterOrganic = struct {
    /// ORGC(0). Surface litter organic carbon (g C).
    organic_carbon_g: f64,
    /// ORGCC(0). Surface litter charcoal carbon (g C).
    charcoal_g: f64,
};

pub const Inputs = struct {
    snow_layer: SnowLayer,
    litter: SurfaceLitter,
    litter_organic: LitterOrganic,
    fluxes: Fluxes,
    parameters: Parameters,
};

pub const UpdatedSnowLayer = struct {
    solid_m3: f64,
    liquid_m3: f64,
    vapor_m3: f64,
    ice_m3: f64,
    heat_capacity_megajoules_per_k: f64,
    temperature_k: f64,
};

pub const UpdatedLitter = struct {
    liquid_m3: f64,
    vapor_m3: f64,
    ice_m3: f64,
    heat_capacity_megajoules_per_k: f64,
    temperature_k: f64,
};

pub const SnowpackTotals = struct {
    /// VOLSS. Total solid snow volume (m3).
    solid_m3: f64,
    /// VOLWS. Total liquid water (m3).
    liquid_m3: f64,
    /// VOLIS. Total ice (m3).
    ice_m3: f64,
    /// VOLS. Total snowpack volume (m3).
    total_m3: f64,
    /// DPTHS. Set to 0 after complete melt-to-litter transfer.
    depth_m: f64,
    /// DENSS value assigned to all snow layers (= DENS0).
    layer_density_megagrams_per_m3: f64,
};

pub const Result = struct {
    snow_layer: UpdatedSnowLayer,
    litter: UpdatedLitter,
    snowpack: SnowpackTotals,
};

// Heat capacity coefficients (MJ m-3 K-1 or MJ g-1 K-1).
const solid_snow_heat_capacity: f64 = 2.095;
const liquid_water_heat_capacity: f64 = 4.19;
const ice_heat_capacity: f64 = 1.9274;
const organic_carbon_heat_capacity: f64 = 2.496e-6; // MJ g-1 K-1

/// Direct translation of REDIST lines 4259--4300 (body of the
/// `IF(XHFLWX > ZEROS)` gate).
///
/// Transfers snow phase volumes from the bottom snow layer (layer 1) to
/// the surface litter (layer 0), updates both heat capacities from the new
/// volumes using source-default phase coefficients, and updates temperatures
/// via energy conservation. Falls back to ambient/mineral temperatures when
/// heat capacity drops below the respective minimum thresholds.
/// Resets all snow layer densities to DENS0 and zeroes DPTHS.
pub fn transfer(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const fluxes = inputs.fluxes;
    const p = inputs.parameters;

    const snow = inputs.snow_layer;
    const litter = inputs.litter;
    const org = inputs.litter_organic;

    // Lines 4260--4263: remove fluxes from snow layer 1.
    const snow_solid = snow.solid_m3 - fluxes.solid_m3;
    const snow_liquid = snow.liquid_m3 - fluxes.liquid_m3;
    const snow_vapor = snow.vapor_m3 - fluxes.vapor_m3;
    const snow_ice = snow.ice_m3 - fluxes.ice_m3;
    inline for (.{ snow_solid, snow_liquid, snow_vapor, snow_ice }) |volume_m3| {
        if (!std.math.isFinite(volume_m3))
            return error.NonFiniteSnowpackLitterTransferResult;
        if (volume_m3 < 0)
            return error.SnowpackLitterFluxExceedsInventory;
    }

    // Lines 4264--4266: add fluxes to surface litter.
    const litter_liquid = litter.liquid_m3 + fluxes.liquid_m3;
    const litter_vapor = litter.vapor_m3 + fluxes.vapor_m3;
    const litter_ice = litter.ice_m3 + fluxes.ice_m3 + fluxes.solid_m3 / p.ice_density_megagrams_per_m3;

    // Lines 4267--4275: snow layer heat capacity and temperature.
    const snow_energy = snow.temperature_k * snow.heat_capacity_megajoules_per_k;
    const snow_vhcp = solid_snow_heat_capacity * snow_solid +
        liquid_water_heat_capacity * (snow_liquid + snow_vapor) +
        ice_heat_capacity * snow_ice;
    const snow_temperature = if (snow_vhcp > p.min_snow_heat_capacity_megajoules_per_k)
        (snow_energy - fluxes.heat_megajoules) / snow_vhcp
    else
        p.ambient_temperature_k;

    // Lines 4276--4284: litter heat capacity and temperature.
    const litter_energy = litter.temperature_k * litter.heat_capacity_megajoules_per_k;
    const litter_vhcp = organic_carbon_heat_capacity * (org.organic_carbon_g + org.charcoal_g) +
        liquid_water_heat_capacity * (litter_liquid + litter_vapor) +
        ice_heat_capacity * litter_ice;
    const litter_temperature = if (litter_vhcp > p.min_litter_heat_capacity_megajoules_per_k)
        (litter_energy + fluxes.heat_megajoules) / litter_vhcp
    else
        p.mineral_layer_temperature_k;

    // Lines 4285--4293: snowpack totals. DENSS(L) = DENS0 for all L;
    // VOLS uses DENSS(1) = DENS0. DPTHS = 0.
    const total_snow_m3 = snow_solid / p.reset_snow_density_megagrams_per_m3 +
        snow_liquid + snow_ice;

    const result = Result{
        .snow_layer = .{
            .solid_m3 = snow_solid,
            .liquid_m3 = snow_liquid,
            .vapor_m3 = snow_vapor,
            .ice_m3 = snow_ice,
            .heat_capacity_megajoules_per_k = snow_vhcp,
            .temperature_k = snow_temperature,
        },
        .litter = .{
            .liquid_m3 = litter_liquid,
            .vapor_m3 = litter_vapor,
            .ice_m3 = litter_ice,
            .heat_capacity_megajoules_per_k = litter_vhcp,
            .temperature_k = litter_temperature,
        },
        .snowpack = .{
            .solid_m3 = snow_solid,
            .liquid_m3 = snow_liquid,
            .ice_m3 = snow_ice,
            .total_m3 = total_snow_m3,
            .depth_m = 0.0,
            .layer_density_megagrams_per_m3 = p.reset_snow_density_megagrams_per_m3,
        },
    };
    inline for (@typeInfo(UpdatedSnowLayer).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.snow_layer, field.name)))
            return error.NonFiniteSnowpackLitterTransferResult;
    inline for (@typeInfo(UpdatedLitter).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.litter, field.name)))
            return error.NonFiniteSnowpackLitterTransferResult;
    inline for (@typeInfo(SnowpackTotals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.snowpack, field.name)))
            return error.NonFiniteSnowpackLitterTransferResult;
    if (result.snow_layer.heat_capacity_megajoules_per_k < 0 or
        result.snow_layer.temperature_k <= 0 or
        result.litter.liquid_m3 < 0 or
        result.litter.vapor_m3 < 0 or
        result.litter.ice_m3 < 0 or
        result.litter.heat_capacity_megajoules_per_k < 0 or
        result.litter.temperature_k <= 0 or
        result.snowpack.total_m3 < 0)
    {
        return error.InvalidSnowpackLitterTransferResult;
    }
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(SnowLayer).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.snow_layer, field.name)))
            return error.NonFiniteSnowpackLitterTransferInput;
    }
    inline for (@typeInfo(SurfaceLitter).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.litter, field.name)))
            return error.NonFiniteSnowpackLitterTransferInput;
    }
    inline for (@typeInfo(LitterOrganic).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.litter_organic, field.name)))
            return error.NonFiniteSnowpackLitterTransferInput;
    }
    inline for (@typeInfo(Fluxes).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.fluxes, field.name)))
            return error.NonFiniteSnowpackLitterTransferInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.parameters, field.name)))
            return error.NonFiniteSnowpackLitterTransferInput;
    }

    const snow = inputs.snow_layer;
    const litter = inputs.litter;
    const organic = inputs.litter_organic;
    const fluxes = inputs.fluxes;
    const p = inputs.parameters;
    if (snow.solid_m3 < 0 or snow.liquid_m3 < 0 or snow.vapor_m3 < 0 or
        snow.ice_m3 < 0 or snow.heat_capacity_megajoules_per_k < 0 or
        snow.temperature_k <= 0 or litter.liquid_m3 < 0 or
        litter.vapor_m3 < 0 or litter.ice_m3 < 0 or
        litter.heat_capacity_megajoules_per_k < 0 or litter.temperature_k <= 0 or
        organic.organic_carbon_g < 0 or organic.charcoal_g < 0)
    {
        return error.InvalidSnowpackLitterTransferState;
    }
    if (fluxes.heat_megajoules <= 0)
        return error.InvalidSnowpackLitterHeatFlux;
    if (fluxes.solid_m3 < 0 or fluxes.liquid_m3 < 0 or
        fluxes.vapor_m3 < 0 or fluxes.ice_m3 < 0)
    {
        return error.InvalidSnowpackLitterFlux;
    }
    if (p.ice_density_megagrams_per_m3 <= 0 or
        p.reset_snow_density_megagrams_per_m3 <= 0 or
        p.min_snow_heat_capacity_megajoules_per_k < 0 or
        p.ambient_temperature_k <= 0 or
        p.min_litter_heat_capacity_megajoules_per_k < 0 or
        p.mineral_layer_temperature_k <= 0)
    {
        return error.InvalidSnowpackLitterTransferParameter;
    }
}

fn defaultInputs() Inputs {
    return .{
        .snow_layer = .{
            .solid_m3 = 0.10,
            .liquid_m3 = 0.05,
            .vapor_m3 = 0.01,
            .ice_m3 = 0.02,
            .temperature_k = 273.0,
            .heat_capacity_megajoules_per_k = 2.095 * 0.10 + 4.19 * (0.05 + 0.01) + 1.9274 * 0.02,
        },
        .litter = .{
            .liquid_m3 = 0.005,
            .vapor_m3 = 0.001,
            .ice_m3 = 0.0,
            .temperature_k = 275.0,
            .heat_capacity_megajoules_per_k = 0.01,
        },
        .litter_organic = .{
            .organic_carbon_g = 500.0,
            .charcoal_g = 10.0,
        },
        .fluxes = .{
            .solid_m3 = 0.05,
            .liquid_m3 = 0.02,
            .vapor_m3 = 0.005,
            .ice_m3 = 0.01,
            .heat_megajoules = 0.1,
        },
        .parameters = .{
            .ice_density_megagrams_per_m3 = 0.917,
            .reset_snow_density_megagrams_per_m3 = 0.10,
            .min_snow_heat_capacity_megajoules_per_k = 1.0e-5,
            .ambient_temperature_k = 273.15,
            .min_litter_heat_capacity_megajoules_per_k = 1.0e-5,
            .mineral_layer_temperature_k = 278.0,
        },
    };
}

test "REDIST snowpack-litter transfer preserves source-order volume balance" {
    const inp = defaultInputs();
    const result = try transfer(inp);

    // Snow layer volumes decrease by fluxes.
    try std.testing.expectApproxEqRel(
        inp.snow_layer.solid_m3 - inp.fluxes.solid_m3,
        result.snow_layer.solid_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqRel(
        inp.snow_layer.liquid_m3 - inp.fluxes.liquid_m3,
        result.snow_layer.liquid_m3,
        1.0e-15,
    );

    // Litter liquid and vapor increase by fluxes.
    try std.testing.expectApproxEqRel(
        inp.litter.liquid_m3 + inp.fluxes.liquid_m3,
        result.litter.liquid_m3,
        1.0e-15,
    );

    // Litter ice includes solid snow converted by ice density.
    const expected_litter_ice =
        inp.litter.ice_m3 + inp.fluxes.ice_m3 +
        inp.fluxes.solid_m3 / inp.parameters.ice_density_megagrams_per_m3;
    try std.testing.expectApproxEqRel(
        expected_litter_ice,
        result.litter.ice_m3,
        1.0e-15,
    );
}

test "REDIST snowpack-litter transfer applies energy conservation when heat capacity sufficient" {
    const inp = defaultInputs();
    const result = try transfer(inp);

    // Snow temperature: (old energy - heat flux) / new heat capacity.
    const snow_energy = inp.snow_layer.temperature_k * inp.snow_layer.heat_capacity_megajoules_per_k;
    const snow_solid_new = inp.snow_layer.solid_m3 - inp.fluxes.solid_m3;
    const snow_liquid_new = inp.snow_layer.liquid_m3 - inp.fluxes.liquid_m3;
    const snow_vapor_new = inp.snow_layer.vapor_m3 - inp.fluxes.vapor_m3;
    const snow_ice_new = inp.snow_layer.ice_m3 - inp.fluxes.ice_m3;
    const snow_vhcp_new = solid_snow_heat_capacity * snow_solid_new +
        liquid_water_heat_capacity * (snow_liquid_new + snow_vapor_new) +
        ice_heat_capacity * snow_ice_new;
    const expected_snow_temp = (snow_energy - inp.fluxes.heat_megajoules) / snow_vhcp_new;
    try std.testing.expectApproxEqRel(
        expected_snow_temp,
        result.snow_layer.temperature_k,
        1.0e-14,
    );
}

test "REDIST snowpack-litter transfer falls back to ambient when snow heat capacity below threshold" {
    var inp = defaultInputs();
    // Force heat capacity below threshold by making snow layer near-empty.
    inp.snow_layer.solid_m3 = 1.0e-10;
    inp.snow_layer.liquid_m3 = 1.0e-10;
    inp.snow_layer.vapor_m3 = 1.0e-10;
    inp.snow_layer.ice_m3 = 1.0e-10;
    inp.snow_layer.heat_capacity_megajoules_per_k = 1.0e-10;
    inp.fluxes.solid_m3 = 1.0e-11;
    inp.fluxes.liquid_m3 = 1.0e-11;
    inp.fluxes.vapor_m3 = 1.0e-11;
    inp.fluxes.ice_m3 = 1.0e-11;
    // VHCPWX set high so new snow VHCP < threshold.
    inp.parameters.min_snow_heat_capacity_megajoules_per_k = 1.0;
    const result = try transfer(inp);
    try std.testing.expectEqual(
        inp.parameters.ambient_temperature_k,
        result.snow_layer.temperature_k,
    );
}

test "REDIST snowpack-litter transfer rejects non-positive heat flux" {
    var inp = defaultInputs();
    inp.fluxes.heat_megajoules = 0.0;
    try std.testing.expectError(
        error.InvalidSnowpackLitterHeatFlux,
        transfer(inp),
    );
    inp.fluxes.heat_megajoules = -1.0;
    try std.testing.expectError(
        error.InvalidSnowpackLitterHeatFlux,
        transfer(inp),
    );
}

test "REDIST snowpack-litter transfer rejects a phase flux larger than its inventory" {
    var inp = defaultInputs();
    inp.fluxes.vapor_m3 = inp.snow_layer.vapor_m3 + 1.0e-6;
    try std.testing.expectError(
        error.SnowpackLitterFluxExceedsInventory,
        transfer(inp),
    );
}

test "REDIST snowpack-litter transfer rejects non-finite state before mutation" {
    var inp = defaultInputs();
    inp.litter.temperature_k = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSnowpackLitterTransferInput,
        transfer(inp),
    );
}
