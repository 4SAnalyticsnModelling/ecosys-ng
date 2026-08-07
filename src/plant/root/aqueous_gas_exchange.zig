const std = @import("std");
const transfer = @import("gas_transfer_transaction.zig");

pub const Inputs = struct {
    prepared: transfer.PreparedTransfer,
    concentrations: transfer.AqueousConcentrationPair,
    convective_mass_flow: transfer.ConvectiveMassFlow,
    aqueous_masses: transfer.AqueousMassPools,
    root_aqueous_volume_m3: f64,
    root_aqueous_non_band_volume_m3: f64,
    root_aqueous_band_volume_m3: f64,
    root_oxygen_soil_flux_g_o_per_plant_step: f64,
    root_oxygen_internal_flux_g_o_per_plant_step: f64,
    population: f64,
    shorter_time_fraction: f64,
    negligible_volume_m3: f64,
    ammonium_non_band_fraction: f64,
    ammonium_band_fraction: f64,
    root_axis_index: usize,
};

pub const Result = struct {
    oxygen_from_soil_g_o_per_step: f64,
    oxygen_from_root_g_o_per_step: f64,
    carbon_dioxide_diffusive_flux_g_c_per_plant_step: f64,
    carbon_dioxide_exchange_g_c_per_step: f64,
    methane_diffusive_flux_g_c_per_plant_step: f64,
    methane_exchange_g_c_per_step: f64,
    nitrous_oxide_diffusive_flux_g_n_per_plant_step: f64,
    nitrous_oxide_exchange_g_n_per_step: f64,
    ammonia_non_band_diffusive_flux_g_n_per_plant_step: f64,
    ammonia_non_band_exchange_g_n_per_step: f64,
    ammonia_band_diffusive_flux_g_n_per_plant_step: f64,
    ammonia_band_exchange_g_n_per_step: f64,
    hydrogen_diffusive_flux_g_h_per_plant_step: f64,
    hydrogen_exchange_g_h_per_step: f64,
};

/// UPTAKE.F 2270--2337. Runtime cells are computed into scratch before the
/// destination is committed, preserving failure atomicity.
pub fn computeCells(
    cells: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (cells.len != scratch.len or cells.len != destination.len)
        return error.RootAqueousExchangeDimensionMismatch;
    for (cells, scratch) |inputs, *candidate|
        candidate.* = try compute(inputs);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const concentration = inputs.concentrations;
    const prepared = inputs.prepared;
    const masses = inputs.aqueous_masses;
    const carbon_dioxide_flux =
        inputs.convective_mass_flow.carbon_dioxide_g_c_per_step +
        prepared.soil_to_root_methane_diffusivity_m3_per_step *
            (concentration.soil.carbon_dioxide_g_c_m3 -
                concentration.root.carbon_dioxide_g_c_m3);
    const carbon_dioxide_gradient = gradientExchange(
        inputs.root_aqueous_volume_m3,
        masses.soil.carbon_dioxide_g_c,
        prepared.root_allocated_water_m3,
        masses.root.carbon_dioxide_g_c,
        prepared.combined_root_soil_water_m3,
        inputs.shorter_time_fraction,
    );
    var result = std.mem.zeroes(Result);
    result.oxygen_from_soil_g_o_per_step =
        inputs.root_oxygen_soil_flux_g_o_per_plant_step * inputs.population;
    result.oxygen_from_root_g_o_per_step =
        inputs.root_oxygen_internal_flux_g_o_per_plant_step * inputs.population;
    result.carbon_dioxide_diffusive_flux_g_c_per_plant_step = carbon_dioxide_flux;
    result.carbon_dioxide_exchange_g_c_per_step =
        signedLimitedExchange(carbon_dioxide_flux, carbon_dioxide_gradient, inputs.population);
    if (inputs.root_axis_index != 1) {
        try validateResult(result);
        return result;
    }

    const methane_flux =
        inputs.convective_mass_flow.methane_g_c_per_step +
        prepared.soil_to_root_methane_diffusivity_m3_per_step *
            (concentration.soil.methane_g_c_m3 -
                concentration.root.methane_g_c_m3);
    const methane_gradient = gradientExchange(
        inputs.root_aqueous_volume_m3,
        masses.soil.methane_g_c,
        prepared.root_allocated_water_m3,
        masses.root.methane_g_c,
        prepared.combined_root_soil_water_m3,
        inputs.shorter_time_fraction,
    );
    result.methane_diffusive_flux_g_c_per_plant_step = methane_flux;
    result.methane_exchange_g_c_per_step =
        signedLimitedExchange(methane_flux, methane_gradient, inputs.population);

    const nitrous_oxide_flux =
        inputs.convective_mass_flow.nitrous_oxide_g_n_per_step +
        prepared.soil_to_root_nitrous_oxide_diffusivity_m3_per_step *
            (concentration.soil.nitrous_oxide_g_n_m3 -
                concentration.root.nitrous_oxide_g_n_m3);
    const nitrous_oxide_gradient = gradientExchange(
        inputs.root_aqueous_volume_m3,
        masses.soil.nitrous_oxide_g_n,
        prepared.root_allocated_water_m3,
        masses.root.nitrous_oxide_g_n,
        prepared.combined_root_soil_water_m3,
        inputs.shorter_time_fraction,
    );
    result.nitrous_oxide_diffusive_flux_g_n_per_plant_step = nitrous_oxide_flux;
    result.nitrous_oxide_exchange_g_n_per_step =
        signedLimitedExchange(nitrous_oxide_flux, nitrous_oxide_gradient, inputs.population);

    const ammonia_root_mass_non_band =
        masses.root.ammonia_non_band_g_n * inputs.ammonium_non_band_fraction;
    const ammonia_non_band_flux =
        inputs.convective_mass_flow.ammonia_non_band_g_n_per_step +
        prepared.soil_to_root_ammonia_non_band_diffusivity_m3_per_step *
            (concentration.soil.ammonia_non_band_g_n_m3 -
                concentration.root.ammonia_non_band_g_n_m3);
    const ammonia_non_band_gradient = if (prepared.non_band_water_m3 > inputs.negligible_volume_m3)
        gradientExchange(
            inputs.root_aqueous_non_band_volume_m3,
            masses.soil.ammonia_non_band_g_n,
            prepared.root_allocated_water_m3 * inputs.ammonium_non_band_fraction,
            ammonia_root_mass_non_band,
            prepared.non_band_water_m3,
            inputs.shorter_time_fraction,
        )
    else
        0;
    result.ammonia_non_band_diffusive_flux_g_n_per_plant_step = ammonia_non_band_flux;
    result.ammonia_non_band_exchange_g_n_per_step =
        signedLimitedExchange(ammonia_non_band_flux, ammonia_non_band_gradient, inputs.population);

    const ammonia_root_mass_band =
        masses.root.ammonia_non_band_g_n * inputs.ammonium_band_fraction;
    const ammonia_band_flux =
        inputs.convective_mass_flow.ammonia_band_g_n_per_step +
        prepared.soil_to_root_ammonia_band_diffusivity_m3_per_step *
            (concentration.soil.ammonia_band_g_n_m3 -
                concentration.root.ammonia_non_band_g_n_m3);
    const ammonia_band_gradient = if (prepared.band_water_m3 > inputs.negligible_volume_m3)
        gradientExchange(
            inputs.root_aqueous_band_volume_m3,
            masses.soil.ammonia_band_g_n,
            prepared.root_allocated_water_m3 * inputs.ammonium_band_fraction,
            ammonia_root_mass_band,
            prepared.band_water_m3,
            inputs.shorter_time_fraction,
        )
    else
        0;
    result.ammonia_band_diffusive_flux_g_n_per_plant_step = ammonia_band_flux;
    result.ammonia_band_exchange_g_n_per_step =
        signedLimitedExchange(ammonia_band_flux, ammonia_band_gradient, inputs.population);

    const hydrogen_flux =
        inputs.convective_mass_flow.hydrogen_g_h_per_step +
        prepared.soil_to_root_hydrogen_diffusivity_m3_per_step *
            (concentration.soil.hydrogen_g_h_m3 -
                concentration.root.hydrogen_g_h_m3);
    const hydrogen_gradient = gradientExchange(
        inputs.root_aqueous_volume_m3,
        masses.soil.hydrogen_g_h,
        prepared.root_allocated_water_m3,
        masses.root.hydrogen_g_h,
        prepared.combined_root_soil_water_m3,
        inputs.shorter_time_fraction,
    );
    result.hydrogen_diffusive_flux_g_h_per_plant_step = hydrogen_flux;
    result.hydrogen_exchange_g_h_per_step =
        signedLimitedExchange(hydrogen_flux, hydrogen_gradient, inputs.population);
    try validateResult(result);
    return result;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteRootAqueousExchangeResult;
}

fn gradientExchange(
    root_volume_m3: f64,
    soil_mass_g: f64,
    allocated_soil_water_m3: f64,
    root_mass_g: f64,
    combined_water_m3: f64,
    shorter_time_fraction: f64,
) f64 {
    return (root_volume_m3 * @max(0, soil_mass_g) -
        allocated_soil_water_m3 * @max(0, root_mass_g)) /
        combined_water_m3 * shorter_time_fraction;
}

fn signedLimitedExchange(flux_per_plant: f64, gradient: f64, population: f64) f64 {
    const population_flux = flux_per_plant * population;
    return if (flux_per_plant > 0)
        @min(@max(0, gradient), population_flux)
    else
        @max(@min(0, gradient), population_flux);
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.root_aqueous_volume_m3,
        inputs.root_aqueous_non_band_volume_m3,
        inputs.root_aqueous_band_volume_m3,
        inputs.population,
        inputs.shorter_time_fraction,
        inputs.negligible_volume_m3,
        inputs.ammonium_non_band_fraction,
        inputs.ammonium_band_fraction,
    }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootAqueousExchangeInput;
    if (!std.math.isFinite(inputs.root_oxygen_soil_flux_g_o_per_plant_step) or
        !std.math.isFinite(inputs.root_oxygen_internal_flux_g_o_per_plant_step))
        return error.InvalidRootAqueousExchangeInput;
    if (inputs.root_aqueous_volume_m3 <= 0 or
        inputs.population <= 0 or
        inputs.shorter_time_fraction <= 0 or
        inputs.ammonium_non_band_fraction > 1 or
        inputs.ammonium_band_fraction > 1 or
        inputs.prepared.combined_root_soil_water_m3 <= 0)
        return error.InvalidRootAqueousExchangeInput;
}

fn testInputs() Inputs {
    var prepared = std.mem.zeroes(transfer.PreparedTransfer);
    prepared.root_allocated_water_m3 = 1;
    prepared.combined_root_soil_water_m3 = 2;
    prepared.non_band_water_m3 = 1;
    prepared.band_water_m3 = 1;
    prepared.soil_to_root_methane_diffusivity_m3_per_step = 0.2;
    prepared.soil_to_root_nitrous_oxide_diffusivity_m3_per_step = 0.2;
    prepared.soil_to_root_ammonia_non_band_diffusivity_m3_per_step = 0.2;
    prepared.soil_to_root_ammonia_band_diffusivity_m3_per_step = 0.2;
    prepared.soil_to_root_hydrogen_diffusivity_m3_per_step = 0.2;
    return .{
        .prepared = prepared,
        .concentrations = .{
            .soil = .{
                .carbon_dioxide_g_c_m3 = 2,
                .oxygen_g_o_m3 = 2,
                .methane_g_c_m3 = 2,
                .nitrous_oxide_g_n_m3 = 2,
                .ammonia_non_band_g_n_m3 = 2,
                .ammonia_band_g_n_m3 = 2,
                .hydrogen_g_h_m3 = 2,
            },
            .root = .{
                .carbon_dioxide_g_c_m3 = 1,
                .oxygen_g_o_m3 = 1,
                .methane_g_c_m3 = 1,
                .nitrous_oxide_g_n_m3 = 1,
                .ammonia_non_band_g_n_m3 = 1,
                .ammonia_band_g_n_m3 = 0,
                .hydrogen_g_h_m3 = 1,
            },
        },
        .convective_mass_flow = std.mem.zeroes(transfer.ConvectiveMassFlow),
        .aqueous_masses = .{
            .soil = .{
                .carbon_dioxide_g_c = 2,
                .oxygen_g_o = 2,
                .methane_g_c = 2,
                .nitrous_oxide_g_n = 2,
                .ammonia_non_band_g_n = 2,
                .ammonia_band_g_n = 2,
                .hydrogen_g_h = 2,
            },
            .root = .{
                .carbon_dioxide_g_c = 1,
                .oxygen_g_o = 1,
                .methane_g_c = 1,
                .nitrous_oxide_g_n = 1,
                .ammonia_non_band_g_n = 1,
                .ammonia_band_g_n = 0,
                .hydrogen_g_h = 1,
            },
        },
        .root_aqueous_volume_m3 = 1,
        .root_aqueous_non_band_volume_m3 = 0.5,
        .root_aqueous_band_volume_m3 = 0.5,
        .root_oxygen_soil_flux_g_o_per_plant_step = 0.3,
        .root_oxygen_internal_flux_g_o_per_plant_step = 0.1,
        .population = 2,
        .shorter_time_fraction = 0.5,
        .negligible_volume_m3 = 1e-12,
        .ammonium_non_band_fraction = 0.6,
        .ammonium_band_fraction = 0.4,
        .root_axis_index = 1,
    };
}

test "positive exchange is bounded by population-scaled diffusive flux" {
    const result = try compute(testInputs());
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result.oxygen_from_soil_g_o_per_step, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.oxygen_from_root_g_o_per_step, 1e-12);
    try std.testing.expect(result.carbon_dioxide_exchange_g_c_per_step >= 0);
    try std.testing.expect(result.carbon_dioxide_exchange_g_c_per_step <=
        result.carbon_dioxide_diffusive_flux_g_c_per_plant_step * 2);
}

test "negative exchange keeps the legacy signed limiter" {
    var inputs = testInputs();
    inputs.concentrations.soil.carbon_dioxide_g_c_m3 = 0;
    inputs.concentrations.root.carbon_dioxide_g_c_m3 = 2;
    inputs.aqueous_masses.soil.carbon_dioxide_g_c = 0;
    inputs.aqueous_masses.root.carbon_dioxide_g_c = 2;
    const result = try compute(inputs);
    try std.testing.expect(result.carbon_dioxide_diffusive_flux_g_c_per_plant_step < 0);
    try std.testing.expect(result.carbon_dioxide_exchange_g_c_per_step <= 0);
    try std.testing.expect(result.carbon_dioxide_exchange_g_c_per_step >=
        result.carbon_dioxide_diffusive_flux_g_c_per_plant_step * inputs.population);
}

test "nonprimary root axis suppresses five legacy gas exchanges" {
    var inputs = testInputs();
    inputs.root_axis_index = 2;
    const result = try compute(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.methane_exchange_g_c_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.nitrous_oxide_exchange_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.ammonia_non_band_exchange_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.ammonia_band_exchange_g_n_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.hydrogen_exchange_g_h_per_step);
}

test "runtime batch failure is atomic" {
    var cells = [_]Inputs{ testInputs(), testInputs() };
    cells[1].population = std.math.nan(f64);
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].oxygen_from_soil_g_o_per_step = 41;
    destination[1].oxygen_from_soil_g_o_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootAqueousExchangeInput,
        computeCells(&cells, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].oxygen_from_soil_g_o_per_step);
    try std.testing.expectEqual(@as(f64, 42), destination[1].oxygen_from_soil_g_o_per_step);
}
