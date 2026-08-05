const std = @import("std");
const aqueous_network = @import("solute_aqueous_network.zig");

pub const LayerReactionAdmissionInputs = struct {
    /// Legacy DLYR(3,L,NY,NX), soil-layer thickness (m).
    layer_thickness_m: f64,
    /// Legacy DLYRM, minimum active soil-layer thickness (m).
    minimum_layer_thickness_m: f64,
    /// Legacy VOLW(L,NY,NX), water volume in the layer (m3).
    water_volume_m3: f64,
    /// Legacy ZEROS2(NY,NX), cell-specific minimum water volume (m3).
    minimum_water_volume_m3: f64,
};

/// Direct translation of SOLUTE.F lines 158--161, the admission gate around
/// the layer preparation and reactions beginning at line 163. Equality with
/// either threshold is inactive because the source uses strict `.GT.` tests.
pub fn admitsLayerReactions(inputs: LayerReactionAdmissionInputs) !bool {
    inline for (@typeInfo(LayerReactionAdmissionInputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLayerReactionAdmissionInput;
    }
    return inputs.layer_thickness_m > inputs.minimum_layer_thickness_m and
        inputs.water_volume_m3 > inputs.minimum_water_volume_m3;
}

test "SOLUTE layer admission preserves strict thickness and water gates" {
    const active: LayerReactionAdmissionInputs = .{
        .layer_thickness_m = 0.1,
        .minimum_layer_thickness_m = 1.0e-6,
        .water_volume_m3 = 0.02,
        .minimum_water_volume_m3 = 1.0e-12,
    };
    try std.testing.expect(try admitsLayerReactions(active));

    var boundary = active;
    boundary.layer_thickness_m = boundary.minimum_layer_thickness_m;
    try std.testing.expect(!try admitsLayerReactions(boundary));
    boundary = active;
    boundary.water_volume_m3 = boundary.minimum_water_volume_m3;
    try std.testing.expect(!try admitsLayerReactions(boundary));
}

test "SOLUTE layer admission rejects invalid physical inputs" {
    const invalid: LayerReactionAdmissionInputs = .{
        .layer_thickness_m = 0.1,
        .minimum_layer_thickness_m = 0,
        .water_volume_m3 = std.math.nan(f64),
        .minimum_water_volume_m3 = 0,
    };
    try std.testing.expectError(
        error.InvalidLayerReactionAdmissionInput,
        admitsLayerReactions(invalid),
    );
}

pub const LayerZoneFractions = struct {
    ammonium_non_band: f64,
    ammonium_band: f64,
    nitrate_non_band: f64,
    nitrate_band: f64,
    phosphate_non_band: f64,
    phosphate_band: f64,
};

pub const LayerZonePreparationInputs = struct {
    /// Legacy VOLW (m3 layer-1).
    water_volume_m3: f64,
    /// Legacy BKVL (Mg layer-1).
    soil_mass_megagrams: f64,
    /// Legacy VOLA (m3 layer-1), used by SOLUTE when soil mass is absent.
    soil_volume_m3: f64,
    fractions: LayerZoneFractions,
    /// Legacy ZEROS threshold applied to BKVL (Mg layer-1).
    positive_soil_mass_threshold_megagrams: f64,
};

/// SOLUTE's zone water volumes and normalization bases. A normalization base
/// has units of Mg when soil mass is present, but intentionally has units of
/// m3 in the source's zero-soil-mass fallback branch.
pub const LayerZonePreparation = struct {
    ammonium_non_band_water_m3: f64,
    ammonium_band_water_m3: f64,
    nitrate_non_band_water_m3: f64,
    nitrate_band_water_m3: f64,
    phosphate_non_band_water_m3: f64,
    phosphate_band_water_m3: f64,
    whole_layer_normalization_basis: f64,
    ammonium_non_band_normalization_basis: f64,
    ammonium_band_normalization_basis: f64,
    nitrate_non_band_normalization_basis: f64,
    nitrate_band_normalization_basis: f64,
    phosphate_non_band_normalization_basis: f64,
    phosphate_band_normalization_basis: f64,
};

/// Direct source-order translation of SOLUTE.F lines 163--193 (`VOLWNH`
/// through `BKVLPB`). The mixed-unit fallback is retained deliberately for
/// translation fidelity and surfaced in the result type's documentation.
pub fn prepareLayerZones(inputs: LayerZonePreparationInputs) !LayerZonePreparation {
    try validateLayerZonePreparationInputs(inputs);

    const ammonium_non_band_water_m3 = inputs.water_volume_m3 * inputs.fractions.ammonium_non_band;
    const ammonium_band_water_m3 = inputs.water_volume_m3 * inputs.fractions.ammonium_band;
    const nitrate_non_band_water_m3 = inputs.water_volume_m3 * inputs.fractions.nitrate_non_band;
    const nitrate_band_water_m3 = inputs.water_volume_m3 * inputs.fractions.nitrate_band;
    const phosphate_non_band_water_m3 = inputs.water_volume_m3 * inputs.fractions.phosphate_non_band;
    const phosphate_band_water_m3 = inputs.water_volume_m3 * inputs.fractions.phosphate_band;

    if (inputs.soil_mass_megagrams > inputs.positive_soil_mass_threshold_megagrams) {
        return .{
            .ammonium_non_band_water_m3 = ammonium_non_band_water_m3,
            .ammonium_band_water_m3 = ammonium_band_water_m3,
            .nitrate_non_band_water_m3 = nitrate_non_band_water_m3,
            .nitrate_band_water_m3 = nitrate_band_water_m3,
            .phosphate_non_band_water_m3 = phosphate_non_band_water_m3,
            .phosphate_band_water_m3 = phosphate_band_water_m3,
            .whole_layer_normalization_basis = inputs.soil_mass_megagrams,
            .ammonium_non_band_normalization_basis = inputs.soil_mass_megagrams * inputs.fractions.ammonium_non_band,
            .ammonium_band_normalization_basis = inputs.soil_mass_megagrams * inputs.fractions.ammonium_band,
            .nitrate_non_band_normalization_basis = inputs.soil_mass_megagrams * inputs.fractions.nitrate_non_band,
            .nitrate_band_normalization_basis = inputs.soil_mass_megagrams * inputs.fractions.nitrate_band,
            .phosphate_non_band_normalization_basis = inputs.soil_mass_megagrams * inputs.fractions.phosphate_non_band,
            .phosphate_band_normalization_basis = inputs.soil_mass_megagrams * inputs.fractions.phosphate_band,
        };
    }

    return .{
        .ammonium_non_band_water_m3 = ammonium_non_band_water_m3,
        .ammonium_band_water_m3 = ammonium_band_water_m3,
        .nitrate_non_band_water_m3 = nitrate_non_band_water_m3,
        .nitrate_band_water_m3 = nitrate_band_water_m3,
        .phosphate_non_band_water_m3 = phosphate_non_band_water_m3,
        .phosphate_band_water_m3 = phosphate_band_water_m3,
        .whole_layer_normalization_basis = inputs.soil_volume_m3,
        .ammonium_non_band_normalization_basis = ammonium_non_band_water_m3,
        .ammonium_band_normalization_basis = ammonium_band_water_m3,
        .nitrate_non_band_normalization_basis = nitrate_non_band_water_m3,
        .nitrate_band_normalization_basis = nitrate_band_water_m3,
        .phosphate_non_band_normalization_basis = phosphate_non_band_water_m3,
        .phosphate_band_normalization_basis = phosphate_band_water_m3,
    };
}

/// Converts one source normalization basis (`BKVL*`) to its corresponding
/// water-volume ratio. Dry or zero-width zones are inactive in SOLUTE and
/// therefore publish zero rather than inventing a denominator.
pub fn normalizationBasisPerWaterVolume(
    normalization_basis: f64,
    water_volume_m3: f64,
) !f64 {
    if (!std.math.isFinite(normalization_basis) or normalization_basis < 0 or
        !std.math.isFinite(water_volume_m3) or water_volume_m3 < 0)
        return error.InvalidLayerZoneNormalization;
    if (water_volume_m3 == 0) return 0;
    const ratio = normalization_basis / water_volume_m3;
    if (!std.math.isFinite(ratio))
        return error.NonFiniteLayerZoneNormalization;
    return ratio;
}

fn validateLayerZonePreparationInputs(inputs: LayerZonePreparationInputs) !void {
    if (!std.math.isFinite(inputs.water_volume_m3) or inputs.water_volume_m3 < 0 or
        !std.math.isFinite(inputs.soil_mass_megagrams) or inputs.soil_mass_megagrams < 0 or
        !std.math.isFinite(inputs.soil_volume_m3) or inputs.soil_volume_m3 < 0 or
        !std.math.isFinite(inputs.positive_soil_mass_threshold_megagrams) or inputs.positive_soil_mass_threshold_megagrams < 0)
        return error.InvalidLayerZonePreparationInput;
    inline for (@typeInfo(LayerZoneFractions).@"struct".fields) |field| {
        const fraction = @field(inputs.fractions, field.name);
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidLayerZoneFraction;
    }
}

test "SOLUTE layer zones use soil mass normalization when BKVL is positive" {
    const result = try prepareLayerZones(.{
        .water_volume_m3 = 8,
        .soil_mass_megagrams = 12,
        .soil_volume_m3 = 10,
        .fractions = .{
            .ammonium_non_band = 0.75,
            .ammonium_band = 0.25,
            .nitrate_non_band = 0.6,
            .nitrate_band = 0.4,
            .phosphate_non_band = 0.8,
            .phosphate_band = 0.2,
        },
        .positive_soil_mass_threshold_megagrams = 1.0e-12,
    });

    try std.testing.expectEqual(@as(f64, 6), result.ammonium_non_band_water_m3);
    try std.testing.expectEqual(@as(f64, 2), result.ammonium_band_water_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 4.8), result.nitrate_non_band_water_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3.2), result.nitrate_band_water_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6.4), result.phosphate_non_band_water_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.6), result.phosphate_band_water_m3, 1.0e-15);
    try std.testing.expectEqual(@as(f64, 12), result.whole_layer_normalization_basis);
    try std.testing.expectEqual(@as(f64, 9), result.ammonium_non_band_normalization_basis);
    try std.testing.expectEqual(@as(f64, 3), result.ammonium_band_normalization_basis);
    try std.testing.expectApproxEqAbs(@as(f64, 7.2), result.nitrate_non_band_normalization_basis, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 4.8), result.nitrate_band_normalization_basis, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.6), result.phosphate_non_band_normalization_basis, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.4), result.phosphate_band_normalization_basis, 1.0e-15);
}

test "SOLUTE layer zones preserve zero-BKVL water-volume fallback" {
    const result = try prepareLayerZones(.{
        .water_volume_m3 = 8,
        .soil_mass_megagrams = 0,
        .soil_volume_m3 = 10,
        .fractions = .{
            .ammonium_non_band = 0.75,
            .ammonium_band = 0.25,
            .nitrate_non_band = 0.6,
            .nitrate_band = 0.4,
            .phosphate_non_band = 0.8,
            .phosphate_band = 0.2,
        },
        .positive_soil_mass_threshold_megagrams = 1.0e-12,
    });

    try std.testing.expectEqual(@as(f64, 10), result.whole_layer_normalization_basis);
    try std.testing.expectEqual(result.ammonium_non_band_water_m3, result.ammonium_non_band_normalization_basis);
    try std.testing.expectEqual(result.ammonium_band_water_m3, result.ammonium_band_normalization_basis);
    try std.testing.expectEqual(result.nitrate_non_band_water_m3, result.nitrate_non_band_normalization_basis);
    try std.testing.expectEqual(result.nitrate_band_water_m3, result.nitrate_band_normalization_basis);
    try std.testing.expectEqual(result.phosphate_non_band_water_m3, result.phosphate_non_band_normalization_basis);
    try std.testing.expectEqual(result.phosphate_band_water_m3, result.phosphate_band_normalization_basis);
}

test "SOLUTE zone mass-water ratios cancel matched source fractions" {
    const prepared = try prepareLayerZones(.{
        .water_volume_m3 = 8,
        .soil_mass_megagrams = 12,
        .soil_volume_m3 = 10,
        .fractions = .{
            .ammonium_non_band = 0.75,
            .ammonium_band = 0.25,
            .nitrate_non_band = 0.6,
            .nitrate_band = 0.4,
            .phosphate_non_band = 0.8,
            .phosphate_band = 0.2,
        },
        .positive_soil_mass_threshold_megagrams = 1.0e-12,
    });
    const expected_ammonium_non_band_ratio =
        prepared.ammonium_non_band_normalization_basis /
        prepared.ammonium_non_band_water_m3;
    const expected_ammonium_band_ratio =
        prepared.ammonium_band_normalization_basis /
        prepared.ammonium_band_water_m3;
    const expected_phosphate_non_band_ratio =
        prepared.phosphate_non_band_normalization_basis /
        prepared.phosphate_non_band_water_m3;
    try std.testing.expectEqual(
        expected_ammonium_non_band_ratio,
        try normalizationBasisPerWaterVolume(
            prepared.ammonium_non_band_normalization_basis,
            prepared.ammonium_non_band_water_m3,
        ),
    );
    try std.testing.expectEqual(
        expected_ammonium_band_ratio,
        try normalizationBasisPerWaterVolume(
            prepared.ammonium_band_normalization_basis,
            prepared.ammonium_band_water_m3,
        ),
    );
    try std.testing.expectEqual(
        expected_phosphate_non_band_ratio,
        try normalizationBasisPerWaterVolume(
            prepared.phosphate_non_band_normalization_basis,
            prepared.phosphate_non_band_water_m3,
        ),
    );
}

test "SOLUTE zero-width zone publishes zero normalization ratio" {
    try std.testing.expectEqual(
        @as(f64, 0),
        try normalizationBasisPerWaterVolume(0, 0),
    );
    try std.testing.expectError(
        error.InvalidLayerZoneNormalization,
        normalizationBasisPerWaterVolume(1, std.math.nan(f64)),
    );
}

pub const UreaInputs = struct {
    broadcast_urea_mol_n: f64,
    banded_urea_mol_n: f64,
    soil_mass_megagrams: f64,
    water_volume_m3: f64,
    biologically_active_water_volume_m3: f64,
    total_microbial_respiration_activity_g_c_per_step: f64,
    temperature_response: f64,
    initial_inhibitor_activity: f64,
    current_inhibitor_activity: f64,
    timestep_h: f64,
};

pub const UreaParameters = struct {
    minimum_half_saturation_mol_n_per_megagram: f64,
    microbial_activity_inhibition_g_c_per_m3_h: f64,
    specific_hydrolysis_mol_n_per_g_c_h: f64,
    inhibitor_decline_rate_per_h: f64,
    negligible_amount: f64,
};

pub const UreaResult = struct {
    broadcast_hydrolysis_mol_n: f64,
    banded_hydrolysis_mol_n: f64,
    effective_half_saturation_mol_n_per_megagram: f64,
    next_inhibitor_activity: f64,
};

pub fn ureaHydrolysis(inputs: UreaInputs, parameters: UreaParameters) !UreaResult {
    try validateUrea(inputs, parameters);
    // Legacy mapping: TOQCK / VOLQ / XNFH in FORTRAN SOLUTE.F.
    const coqck = if (inputs.biologically_active_water_volume_m3 > parameters.negligible_amount)
        @min(0.1e6, inputs.total_microbial_respiration_activity_g_c_per_step / (inputs.biologically_active_water_volume_m3 * inputs.timestep_h))
    else
        0.1e6;
    const effective_half_saturation = parameters.minimum_half_saturation_mol_n_per_megagram * (1 + coqck / parameters.microbial_activity_inhibition_g_c_per_m3_h);
    var next_inhibitor: f64 = 0;
    if (inputs.initial_inhibitor_activity > parameters.negligible_amount and inputs.current_inhibitor_activity > parameters.negligible_amount) {
        const decline_per_step = parameters.inhibitor_decline_rate_per_h * inputs.timestep_h;
        next_inhibitor = inputs.current_inhibitor_activity - decline_per_step * inputs.current_inhibitor_activity * @max(decline_per_step, 1 - inputs.current_inhibitor_activity / inputs.initial_inhibitor_activity);
        if (next_inhibitor < 0 or next_inhibitor > 1) return error.InvalidUreaInhibitorEvolution;
    }
    const specific_hydrolysis = parameters.specific_hydrolysis_mol_n_per_g_c_h * inputs.timestep_h;
    const broadcast_concentration = fertilizerConcentration(inputs.broadcast_urea_mol_n, inputs.soil_mass_megagrams, inputs.water_volume_m3, parameters.negligible_amount);
    const banded_concentration = fertilizerConcentration(inputs.banded_urea_mol_n, inputs.soil_mass_megagrams, inputs.water_volume_m3, parameters.negligible_amount);
    const broadcast_limitation = broadcast_concentration / (broadcast_concentration + effective_half_saturation);
    const banded_limitation = banded_concentration / (banded_concentration + effective_half_saturation);
    const common_capacity = specific_hydrolysis * inputs.total_microbial_respiration_activity_g_c_per_step * inputs.temperature_response * (1 - next_inhibitor);
    return .{
        .broadcast_hydrolysis_mol_n = @min(inputs.broadcast_urea_mol_n, common_capacity * broadcast_limitation),
        .banded_hydrolysis_mol_n = @min(inputs.banded_urea_mol_n, common_capacity * banded_limitation),
        .effective_half_saturation_mol_n_per_megagram = effective_half_saturation,
        .next_inhibitor_activity = next_inhibitor,
    };
}

pub const FertilizerState = struct {
    broadcast_ammonium_mol_n: f64,
    broadcast_ammonia_mol_n: f64,
    broadcast_urea_mol_n: f64,
    broadcast_nitrate_mol_n: f64,
    banded_ammonium_mol_n: f64,
    banded_ammonia_mol_n: f64,
    banded_urea_mol_n: f64,
    banded_nitrate_mol_n: f64,
};

pub const ZoneFractions = struct { ammonium_non_band: f64, ammonium_band: f64, nitrate_non_band: f64, nitrate_band: f64 };
pub const DissolutionRates = struct { ammonium_per_h: f64, ammonia_per_h: f64, nitrate_per_h: f64 };

pub const DissolutionFlux = struct {
    broadcast_ammonium_non_band_mol_n: f64,
    broadcast_ammonia_non_band_mol_n: f64,
    broadcast_urea_non_band_mol_n: f64,
    broadcast_nitrate_non_band_mol_n: f64,
    broadcast_ammonium_band_mol_n: f64,
    broadcast_ammonia_band_mol_n: f64,
    broadcast_urea_band_mol_n: f64,
    broadcast_nitrate_band_mol_n: f64,
    banded_ammonium_mol_n: f64,
    banded_ammonia_mol_n: f64,
    banded_urea_mol_n: f64,
    banded_nitrate_mol_n: f64,
};

/// Direct source-order translation of SOLUTE.F lines 326--338.
pub fn dissolution(state: FertilizerState, hydrolysis: UreaResult, fractions: ZoneFractions, rates: DissolutionRates, water_content_m3_per_m3: f64, timestep_h: f64) !DissolutionFlux {
    try validateFertilizerState(state);
    inline for (@typeInfo(ZoneFractions).@"struct".fields) |field| if (!std.math.isFinite(@field(fractions, field.name)) or @field(fractions, field.name) < 0 or @field(fractions, field.name) > 1) return error.InvalidFertilizerZoneFraction;
    inline for (@typeInfo(DissolutionRates).@"struct".fields) |field| if (!std.math.isFinite(@field(rates, field.name)) or @field(rates, field.name) < 0) return error.InvalidFertilizerDissolutionRate;
    if (!std.math.isFinite(water_content_m3_per_m3) or water_content_m3_per_m3 < 0 or !std.math.isFinite(timestep_h) or timestep_h <= 0) return error.InvalidFertilizerDissolutionInput;
    const ammonium_rate = rates.ammonium_per_h * timestep_h;
    const ammonia_rate = rates.ammonia_per_h * timestep_h;
    const nitrate_rate = rates.nitrate_per_h * timestep_h;
    const result: DissolutionFlux = .{
        .broadcast_ammonium_non_band_mol_n = ammonium_rate * state.broadcast_ammonium_mol_n * fractions.ammonium_non_band * water_content_m3_per_m3,
        .broadcast_ammonia_non_band_mol_n = ammonia_rate * state.broadcast_ammonia_mol_n * fractions.ammonium_non_band,
        .broadcast_urea_non_band_mol_n = hydrolysis.broadcast_hydrolysis_mol_n * fractions.ammonium_non_band,
        .broadcast_nitrate_non_band_mol_n = nitrate_rate * state.broadcast_nitrate_mol_n * fractions.nitrate_non_band * water_content_m3_per_m3,
        .broadcast_ammonium_band_mol_n = ammonium_rate * state.broadcast_ammonium_mol_n * fractions.ammonium_band * water_content_m3_per_m3,
        .broadcast_ammonia_band_mol_n = ammonia_rate * state.broadcast_ammonia_mol_n * fractions.ammonium_band,
        .broadcast_urea_band_mol_n = hydrolysis.broadcast_hydrolysis_mol_n * fractions.ammonium_band,
        .broadcast_nitrate_band_mol_n = nitrate_rate * state.broadcast_nitrate_mol_n * fractions.nitrate_band * water_content_m3_per_m3,
        .banded_ammonium_mol_n = ammonium_rate * state.banded_ammonium_mol_n * water_content_m3_per_m3,
        .banded_ammonia_mol_n = ammonia_rate * state.banded_ammonia_mol_n,
        .banded_urea_mol_n = hydrolysis.banded_hydrolysis_mol_n * fractions.ammonium_band,
        .banded_nitrate_mol_n = nitrate_rate * state.banded_nitrate_mol_n * water_content_m3_per_m3,
    };
    try validateFluxAgainstState(result, state);
    return result;
}

pub fn commit(state: *FertilizerState, flux: DissolutionFlux) !void {
    try validateFertilizerState(state.*);
    try validateFluxAgainstState(flux, state.*);
    var next = state.*;
    next.broadcast_ammonium_mol_n -= flux.broadcast_ammonium_non_band_mol_n + flux.broadcast_ammonium_band_mol_n;
    next.broadcast_ammonia_mol_n -= flux.broadcast_ammonia_non_band_mol_n + flux.broadcast_ammonia_band_mol_n;
    next.broadcast_urea_mol_n -= flux.broadcast_urea_non_band_mol_n + flux.broadcast_urea_band_mol_n;
    next.broadcast_nitrate_mol_n -= flux.broadcast_nitrate_non_band_mol_n + flux.broadcast_nitrate_band_mol_n;
    next.banded_ammonium_mol_n -= flux.banded_ammonium_mol_n;
    next.banded_ammonia_mol_n -= flux.banded_ammonia_mol_n;
    next.banded_urea_mol_n -= flux.banded_urea_mol_n;
    next.banded_nitrate_mol_n -= flux.banded_nitrate_mol_n;
    try validateFertilizerState(next);
    state.* = next;
}

/// Publishes the SOLUTE fertilizer recipients atomically. NH3 fertilizer
/// enters soil gas, while hydrolyzed urea enters aqueous NH3.
pub fn commitToRecipients(
    state: *FertilizerState,
    aqueous: *aqueous_network.State,
    gaseous_ammonia_g_n: *f64,
    flux: DissolutionFlux,
    fractions: ZoneFractions,
    water_volume_m3: f64,
    nitrogen_molar_mass_g_per_mol: f64,
) !void {
    try validateFertilizerState(state.*);
    try validateFluxAgainstState(flux, state.*);
    if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or
        !std.math.isFinite(nitrogen_molar_mass_g_per_mol) or
        nitrogen_molar_mass_g_per_mol <= 0 or
        !std.math.isFinite(gaseous_ammonia_g_n.*) or gaseous_ammonia_g_n.* < 0)
        return error.InvalidFertilizerRecipientWaterVolume;
    inline for (@typeInfo(ZoneFractions).@"struct".fields) |field| {
        const fraction = @field(fractions, field.name);
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidFertilizerZoneFraction;
    }

    const ammonium_non_band_water =
        water_volume_m3 * fractions.ammonium_non_band;
    const ammonium_band_water =
        water_volume_m3 * fractions.ammonium_band;
    const nitrate_non_band_water =
        water_volume_m3 * fractions.nitrate_non_band;
    const nitrate_band_water =
        water_volume_m3 * fractions.nitrate_band;
    try requireRecipientVolume(
        flux.broadcast_ammonium_non_band_mol_n,
        ammonium_non_band_water,
    );
    try requireRecipientVolume(
        flux.broadcast_ammonium_band_mol_n +
            flux.banded_ammonium_mol_n,
        ammonium_band_water,
    );
    try requireRecipientVolume(
        flux.broadcast_urea_non_band_mol_n,
        ammonium_non_band_water,
    );
    try requireRecipientVolume(
        flux.broadcast_urea_band_mol_n +
            flux.banded_urea_mol_n,
        ammonium_band_water,
    );
    try requireRecipientVolume(
        flux.broadcast_nitrate_non_band_mol_n,
        nitrate_non_band_water,
    );
    try requireRecipientVolume(
        flux.broadcast_nitrate_band_mol_n +
            flux.banded_nitrate_mol_n,
        nitrate_band_water,
    );

    var next_state = state.*;
    var next_aqueous = aqueous.*;
    var next_gaseous_ammonia_g_n = gaseous_ammonia_g_n.*;
    try commit(&next_state, flux);
    next_aqueous.ammonium_non_band +=
        flux.broadcast_ammonium_non_band_mol_n /
        nonzeroOrOne(ammonium_non_band_water);
    next_aqueous.ammonium_band +=
        (flux.broadcast_ammonium_band_mol_n +
            flux.banded_ammonium_mol_n) /
        nonzeroOrOne(ammonium_band_water);
    next_aqueous.ammonia_non_band +=
        flux.broadcast_urea_non_band_mol_n /
        nonzeroOrOne(ammonium_non_band_water);
    next_aqueous.ammonia_band +=
        (flux.broadcast_urea_band_mol_n +
            flux.banded_urea_mol_n) /
        nonzeroOrOne(ammonium_band_water);
    next_gaseous_ammonia_g_n +=
        (flux.broadcast_ammonia_non_band_mol_n +
            flux.broadcast_ammonia_band_mol_n +
            flux.banded_ammonia_mol_n) *
        nitrogen_molar_mass_g_per_mol;
    next_aqueous.nitrate_non_band +=
        flux.broadcast_nitrate_non_band_mol_n /
        nonzeroOrOne(nitrate_non_band_water);
    next_aqueous.nitrate_band +=
        (flux.broadcast_nitrate_band_mol_n +
            flux.banded_nitrate_mol_n) /
        nonzeroOrOne(nitrate_band_water);
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        const value = @field(next_aqueous, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFertilizerRecipientAqueousState;
    }
    if (!std.math.isFinite(next_gaseous_ammonia_g_n) or
        next_gaseous_ammonia_g_n < 0)
        return error.InvalidFertilizerRecipientGasState;
    state.* = next_state;
    aqueous.* = next_aqueous;
    gaseous_ammonia_g_n.* = next_gaseous_ammonia_g_n;
}

fn requireRecipientVolume(amount_mol_n: f64, volume_m3: f64) !void {
    if (amount_mol_n > 0 and volume_m3 <= 0)
        return error.MissingFertilizerRecipientWaterVolume;
}

fn nonzeroOrOne(value: f64) f64 {
    return if (value > 0) value else 1;
}

fn fertilizerConcentration(amount_mol_n: f64, soil_mass_megagrams: f64, water_volume_m3: f64, negligible: f64) f64 {
    if (amount_mol_n > negligible and soil_mass_megagrams > negligible) return amount_mol_n / soil_mass_megagrams;
    if (water_volume_m3 > negligible) return amount_mol_n / water_volume_m3;
    return 0;
}

fn validateUrea(inputs: UreaInputs, parameters: UreaParameters) !void {
    inline for (@typeInfo(UreaInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name)) or @field(inputs, field.name) < 0) return error.InvalidUreaHydrolysisInput;
    inline for (@typeInfo(UreaParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name)) or @field(parameters, field.name) < 0) return error.InvalidUreaHydrolysisParameter;
    if (inputs.current_inhibitor_activity > 1 or inputs.timestep_h <= 0 or parameters.minimum_half_saturation_mol_n_per_megagram <= 0 or parameters.microbial_activity_inhibition_g_c_per_m3_h <= 0) return error.InvalidUreaHydrolysisInput;
}

fn validateFertilizerState(state: FertilizerState) !void {
    inline for (@typeInfo(FertilizerState).@"struct".fields) |field| if (!std.math.isFinite(@field(state, field.name)) or @field(state, field.name) < -1e-14) return error.InvalidFertilizerState;
}

fn validateFluxAgainstState(flux: DissolutionFlux, state: FertilizerState) !void {
    inline for (@typeInfo(DissolutionFlux).@"struct".fields) |field| if (!std.math.isFinite(@field(flux, field.name)) or @field(flux, field.name) < 0) return error.InvalidFertilizerDissolutionFlux;
    if (flux.broadcast_ammonium_non_band_mol_n + flux.broadcast_ammonium_band_mol_n > state.broadcast_ammonium_mol_n or flux.broadcast_ammonia_non_band_mol_n + flux.broadcast_ammonia_band_mol_n > state.broadcast_ammonia_mol_n or flux.broadcast_urea_non_band_mol_n + flux.broadcast_urea_band_mol_n > state.broadcast_urea_mol_n or flux.broadcast_nitrate_non_band_mol_n + flux.broadcast_nitrate_band_mol_n > state.broadcast_nitrate_mol_n or flux.banded_ammonium_mol_n > state.banded_ammonium_mol_n or flux.banded_ammonia_mol_n > state.banded_ammonia_mol_n or flux.banded_urea_mol_n > state.banded_urea_mol_n or flux.banded_nitrate_mol_n > state.banded_nitrate_mol_n) return error.FertilizerDissolutionExceedsPool;
}

test "urea hydrolysis and fertilizer dissolution conserve fertilizer nitrogen" {
    var state: FertilizerState = .{ .broadcast_ammonium_mol_n = 1, .broadcast_ammonia_mol_n = 1, .broadcast_urea_mol_n = 1, .broadcast_nitrate_mol_n = 1, .banded_ammonium_mol_n = 1, .banded_ammonia_mol_n = 1, .banded_urea_mol_n = 1, .banded_nitrate_mol_n = 1 };
    const hydrolysis = try ureaHydrolysis(.{ .broadcast_urea_mol_n = 1, .banded_urea_mol_n = 1, .soil_mass_megagrams = 1, .water_volume_m3 = 1, .biologically_active_water_volume_m3 = 1, .total_microbial_respiration_activity_g_c_per_step = 0.1, .temperature_response = 1, .initial_inhibitor_activity = 1, .current_inhibitor_activity = 0.5, .timestep_h = 1 }, .{ .minimum_half_saturation_mol_n_per_megagram = 0.05, .microbial_activity_inhibition_g_c_per_m3_h = 50, .specific_hydrolysis_mol_n_per_g_c_h = 0.03, .inhibitor_decline_rate_per_h = 0.01, .negligible_amount = 1e-12 });
    const flux = try dissolution(state, hydrolysis, .{ .ammonium_non_band = 0.7, .ammonium_band = 0.3, .nitrate_non_band = 0.6, .nitrate_band = 0.4 }, .{ .ammonium_per_h = 0.1, .ammonia_per_h = 0.1, .nitrate_per_h = 0.1 }, 0.5, 1);
    const before = sumState(state);
    const dissolved = sumFlux(flux);
    try commit(&state, flux);
    try std.testing.expectApproxEqAbs(before, sumState(state) + dissolved, 1e-13);
}

test "SOLUTE lines 327-338 retain fertilizer flux assignment order" {
    const expected_names = [_][]const u8{
        "broadcast_ammonium_non_band_mol_n",
        "broadcast_ammonia_non_band_mol_n",
        "broadcast_urea_non_band_mol_n",
        "broadcast_nitrate_non_band_mol_n",
        "broadcast_ammonium_band_mol_n",
        "broadcast_ammonia_band_mol_n",
        "broadcast_urea_band_mol_n",
        "broadcast_nitrate_band_mol_n",
        "banded_ammonium_mol_n",
        "banded_ammonia_mol_n",
        "banded_urea_mol_n",
        "banded_nitrate_mol_n",
    };
    inline for (@typeInfo(DissolutionFlux).@"struct".fields, 0..) |field, index|
        try std.testing.expectEqualStrings(expected_names[index], field.name);
}

test "urea hydrolysis uses legacy TOQCK / VOLQ / TFNQ operands" {
    const toqck_g_c_per_step = 45.0;
    const volq_m3 = 3.0;
    const temperature_factor = 1.5;
    const soil_mass_megagrams = 2.0;
    const water_volume_m3 = 6.0;
    const broadcast_urea_mol_n = 10.0;
    const banded_urea_mol_n = 4.0;
    const specific_hydrolysis = 0.2;
    const microbial_activity_inhibition = 15.0;
    const minimum_half_saturation = 0.5;
    const initial_inhibitor = 0.8;
    const current_inhibitor = 0.4;
    const inhibitor_decline = 0.05;
    const timestep_h = 1.0;

    const coqck = @min(0.1e6, toqck_g_c_per_step / (volq_m3 * timestep_h));
    const expected_effective_half_saturation = minimum_half_saturation * (1 + coqck / microbial_activity_inhibition);
    const broadcast_fraction = (broadcast_urea_mol_n / soil_mass_megagrams) / (broadcast_urea_mol_n / soil_mass_megagrams + expected_effective_half_saturation);
    const banded_fraction = (banded_urea_mol_n / soil_mass_megagrams) / (banded_urea_mol_n / soil_mass_megagrams + expected_effective_half_saturation);
    const expected_decline_per_step = inhibitor_decline * timestep_h;
    const expected_inhibitor_decline_rate = expected_decline_per_step * current_inhibitor *
        @max(expected_decline_per_step, 1 - current_inhibitor / initial_inhibitor);
    const expected_next_inhibitor = current_inhibitor - expected_inhibitor_decline_rate;
    const expected_common_hydrolysis_capacity = specific_hydrolysis * timestep_h *
        toqck_g_c_per_step * temperature_factor * (1 - expected_next_inhibitor);
    const result = try ureaHydrolysis(
        .{
            .broadcast_urea_mol_n = broadcast_urea_mol_n,
            .banded_urea_mol_n = banded_urea_mol_n,
            .soil_mass_megagrams = soil_mass_megagrams,
            .water_volume_m3 = water_volume_m3,
            .biologically_active_water_volume_m3 = volq_m3,
            .total_microbial_respiration_activity_g_c_per_step = toqck_g_c_per_step,
            .temperature_response = temperature_factor,
            .initial_inhibitor_activity = initial_inhibitor,
            .current_inhibitor_activity = current_inhibitor,
            .timestep_h = timestep_h,
        },
        .{
            .minimum_half_saturation_mol_n_per_megagram = minimum_half_saturation,
            .microbial_activity_inhibition_g_c_per_m3_h = microbial_activity_inhibition,
            .specific_hydrolysis_mol_n_per_g_c_h = specific_hydrolysis,
            .inhibitor_decline_rate_per_h = inhibitor_decline,
            .negligible_amount = 1e-12,
        },
    );
    try std.testing.expectApproxEqAbs(@min(broadcast_urea_mol_n, expected_common_hydrolysis_capacity * broadcast_fraction), result.broadcast_hydrolysis_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(@min(banded_urea_mol_n, expected_common_hydrolysis_capacity * banded_fraction), result.banded_hydrolysis_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_effective_half_saturation, result.effective_half_saturation_mol_n_per_megagram, 1e-12);
    try std.testing.expectApproxEqAbs(expected_next_inhibitor, result.next_inhibitor_activity, 1e-12);
}

test "first fertilizer boundary commits legacy-formula hydrolysis and dissolution to aqueous state" {
    const toqck_g_c_per_step = 90.0;
    const volq_m3 = 4.0;
    const tfnq = 1.25;
    const timestep_h = 1.0;
    const water_volume_m3 = 12.0;
    const soil_mass_megagrams = 8.0;
    const water_content_m3_per_m3 = 0.6;
    const broadcast_urea_mol_n = 10.0;
    const banded_urea_mol_n = 4.0;
    const broadcast_ammonium_mol_n = 6.0;
    const banded_ammonium_mol_n = 2.0;
    const broadcast_ammonia_mol_n = 3.0;
    const banded_ammonia_mol_n = 1.0;
    const broadcast_nitrate_mol_n = 5.0;
    const banded_nitrate_mol_n = 7.0;

    const fractions: ZoneFractions = .{
        .ammonium_non_band = 0.3,
        .ammonium_band = 0.7,
        .nitrate_non_band = 0.6,
        .nitrate_band = 0.4,
    };
    const rates: DissolutionRates = .{
        .ammonium_per_h = 0.2,
        .ammonia_per_h = 0.1,
        .nitrate_per_h = 0.4,
    };
    const initial_state: FertilizerState = .{
        .broadcast_ammonium_mol_n = broadcast_ammonium_mol_n,
        .broadcast_ammonia_mol_n = broadcast_ammonia_mol_n,
        .broadcast_urea_mol_n = broadcast_urea_mol_n,
        .broadcast_nitrate_mol_n = broadcast_nitrate_mol_n,
        .banded_ammonium_mol_n = banded_ammonium_mol_n,
        .banded_ammonia_mol_n = banded_ammonia_mol_n,
        .banded_urea_mol_n = banded_urea_mol_n,
        .banded_nitrate_mol_n = banded_nitrate_mol_n,
    };
    const fertilizer_parameters: UreaParameters = .{
        .minimum_half_saturation_mol_n_per_megagram = 1.5,
        .microbial_activity_inhibition_g_c_per_m3_h = 30,
        .specific_hydrolysis_mol_n_per_g_c_h = 0.018,
        .inhibitor_decline_rate_per_h = 0.0,
        .negligible_amount = 1e-12,
    };

    const hydrolysis = try ureaHydrolysis(
        .{
            .broadcast_urea_mol_n = broadcast_urea_mol_n,
            .banded_urea_mol_n = banded_urea_mol_n,
            .soil_mass_megagrams = soil_mass_megagrams,
            .water_volume_m3 = water_volume_m3,
            .biologically_active_water_volume_m3 = volq_m3,
            .total_microbial_respiration_activity_g_c_per_step = toqck_g_c_per_step,
            .temperature_response = tfnq,
            .initial_inhibitor_activity = 1.0,
            .current_inhibitor_activity = 0.2,
            .timestep_h = timestep_h,
        },
        fertilizer_parameters,
    );
    const coqck = @min(1.0e6, toqck_g_c_per_step / (volq_m3 * timestep_h));
    const expected_effective_half_saturation = fertilizer_parameters.minimum_half_saturation_mol_n_per_megagram * (1 + coqck / fertilizer_parameters.microbial_activity_inhibition_g_c_per_m3_h);
    const broadcast_urea_concentration = broadcast_urea_mol_n / soil_mass_megagrams;
    const banded_urea_concentration = banded_urea_mol_n / soil_mass_megagrams;
    const broadcast_limitation = broadcast_urea_concentration / (broadcast_urea_concentration + expected_effective_half_saturation);
    const banded_limitation = banded_urea_concentration / (banded_urea_concentration + expected_effective_half_saturation);
    const expected_common_hydrolysis = fertilizer_parameters.specific_hydrolysis_mol_n_per_g_c_h * toqck_g_c_per_step * tfnq * (1.0 - 0.2);
    const expected_broadcast_hydrolysis_mol_n = @min(broadcast_urea_mol_n, expected_common_hydrolysis * broadcast_limitation);
    const expected_banded_hydrolysis_mol_n = @min(banded_urea_mol_n, expected_common_hydrolysis * banded_limitation);
    try std.testing.expectApproxEqAbs(expected_broadcast_hydrolysis_mol_n, hydrolysis.broadcast_hydrolysis_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_banded_hydrolysis_mol_n, hydrolysis.banded_hydrolysis_mol_n, 1e-12);

    const dissolved = try dissolution(
        initial_state,
        hydrolysis,
        fractions,
        rates,
        water_content_m3_per_m3,
        timestep_h,
    );
    const expected_dissolution_broadcast_non_band_ammonium = rates.ammonium_per_h * broadcast_ammonium_mol_n * fractions.ammonium_non_band * water_content_m3_per_m3;
    const expected_dissolution_broadcast_band_ammonium = rates.ammonium_per_h * broadcast_ammonium_mol_n * fractions.ammonium_band * water_content_m3_per_m3;
    const expected_dissolution_banded_ammonium = rates.ammonium_per_h * banded_ammonium_mol_n * water_content_m3_per_m3;
    const expected_dissolution_broadcast_non_band_ammonia = rates.ammonia_per_h * broadcast_ammonia_mol_n * fractions.ammonium_non_band;
    const expected_dissolution_broadcast_band_ammonia = rates.ammonia_per_h * broadcast_ammonia_mol_n * fractions.ammonium_band;
    const expected_dissolution_banded_ammonia = rates.ammonia_per_h * banded_ammonia_mol_n;
    const expected_dissolution_non_band_nitrate = rates.nitrate_per_h * broadcast_nitrate_mol_n * fractions.nitrate_non_band * water_content_m3_per_m3;
    const expected_dissolution_broadcast_band_nitrate = rates.nitrate_per_h * broadcast_nitrate_mol_n * fractions.nitrate_band * water_content_m3_per_m3;
    const expected_dissolution_banded_nitrate = rates.nitrate_per_h * banded_nitrate_mol_n * water_content_m3_per_m3;
    const expected_broadcast_urea_to_non_band = expected_broadcast_hydrolysis_mol_n * fractions.ammonium_non_band;
    const expected_broadcast_urea_to_band = expected_broadcast_hydrolysis_mol_n * fractions.ammonium_band;
    const expected_banded_urea_to_band = expected_banded_hydrolysis_mol_n * fractions.ammonium_band;

    try std.testing.expectApproxEqAbs(expected_dissolution_broadcast_non_band_ammonium, dissolved.broadcast_ammonium_non_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_broadcast_band_ammonium, dissolved.broadcast_ammonium_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_banded_ammonium, dissolved.banded_ammonium_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_broadcast_non_band_ammonia, dissolved.broadcast_ammonia_non_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_broadcast_band_ammonia, dissolved.broadcast_ammonia_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_banded_ammonia, dissolved.banded_ammonia_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_broadcast_urea_to_non_band, dissolved.broadcast_urea_non_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_broadcast_urea_to_band, dissolved.broadcast_urea_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_banded_urea_to_band, dissolved.banded_urea_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_non_band_nitrate, dissolved.broadcast_nitrate_non_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_broadcast_band_nitrate, dissolved.broadcast_nitrate_band_mol_n, 1e-12);
    try std.testing.expectApproxEqAbs(expected_dissolution_banded_nitrate, dissolved.banded_nitrate_mol_n, 1e-12);

    var state = initial_state;
    var aqueous = std.mem.zeroes(aqueous_network.State);
    var gaseous_ammonia_g_n: f64 = 0;
    try commitToRecipients(&state, &aqueous, &gaseous_ammonia_g_n, dissolved, fractions, water_volume_m3, 14);

    try std.testing.expectApproxEqAbs(
        broadcast_urea_mol_n - expected_broadcast_hydrolysis_mol_n,
        state.broadcast_urea_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        banded_urea_mol_n - expected_banded_urea_to_band,
        state.banded_urea_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        broadcast_ammonium_mol_n -
            expected_dissolution_broadcast_non_band_ammonium -
            expected_dissolution_broadcast_band_ammonium,
        state.broadcast_ammonium_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        banded_ammonium_mol_n - expected_dissolution_banded_ammonium,
        state.banded_ammonium_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        broadcast_ammonia_mol_n -
            expected_dissolution_broadcast_non_band_ammonia -
            expected_dissolution_broadcast_band_ammonia,
        state.broadcast_ammonia_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        banded_ammonia_mol_n - expected_dissolution_banded_ammonia,
        state.banded_ammonia_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        broadcast_nitrate_mol_n -
            expected_dissolution_non_band_nitrate -
            expected_dissolution_broadcast_band_nitrate,
        state.broadcast_nitrate_mol_n,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        banded_nitrate_mol_n - expected_dissolution_banded_nitrate,
        state.banded_nitrate_mol_n,
        1e-12,
    );

    const non_band_water_m3 = water_volume_m3 * fractions.ammonium_non_band;
    const band_water_m3 = water_volume_m3 * fractions.ammonium_band;
    const nitrate_non_band_water_m3 = water_volume_m3 * fractions.nitrate_non_band;
    const nitrate_band_water_m3 = water_volume_m3 * fractions.nitrate_band;
    const expected_aqueous_ammonium_non_band = expected_dissolution_broadcast_non_band_ammonium / non_band_water_m3;
    const expected_aqueous_ammonium_band = (expected_dissolution_broadcast_band_ammonium +
        expected_dissolution_banded_ammonium) / band_water_m3;
    const expected_aqueous_ammonia_non_band = expected_broadcast_urea_to_non_band / non_band_water_m3;
    const expected_aqueous_ammonia_band = (expected_broadcast_urea_to_band +
        expected_banded_urea_to_band) / band_water_m3;
    const expected_aqueous_nitrate_non_band = expected_dissolution_non_band_nitrate / nitrate_non_band_water_m3;
    const expected_aqueous_nitrate_band = (expected_dissolution_broadcast_band_nitrate +
        expected_dissolution_banded_nitrate) / nitrate_band_water_m3;

    try std.testing.expectApproxEqAbs(expected_aqueous_ammonium_non_band, aqueous.ammonium_non_band, 1e-12);
    try std.testing.expectApproxEqAbs(expected_aqueous_ammonium_band, aqueous.ammonium_band, 1e-12);
    try std.testing.expectApproxEqAbs(expected_aqueous_ammonia_non_band, aqueous.ammonia_non_band, 1e-12);
    try std.testing.expectApproxEqAbs(expected_aqueous_ammonia_band, aqueous.ammonia_band, 1e-12);
    try std.testing.expectApproxEqAbs(expected_aqueous_nitrate_non_band, aqueous.nitrate_non_band, 1e-12);
    try std.testing.expectApproxEqAbs(expected_aqueous_nitrate_band, aqueous.nitrate_band, 1e-12);
    try std.testing.expectApproxEqAbs(
        14 * (expected_dissolution_broadcast_non_band_ammonia +
            expected_dissolution_broadcast_band_ammonia +
            expected_dissolution_banded_ammonia),
        gaseous_ammonia_g_n,
        1e-12,
    );
}

test "zero biologically active water uses COQCK saturation cap" {
    const toqck_g_c_per_step = 10.0;
    const volq_m3 = 0.0;
    const temperature_factor = 0.8;
    const soil_mass_megagrams = 1.0;
    const water_volume_m3 = 3.0;
    const broadcast_urea_mol_n = 5.0;
    const expected_coqck = 0.1e6;
    const expected_effective_half_saturation = 0.2 * (1 + expected_coqck / 10.0);
    const expected_fraction = (broadcast_urea_mol_n / soil_mass_megagrams) / (broadcast_urea_mol_n / soil_mass_megagrams + expected_effective_half_saturation);
    const expected_next_inhibitor = 0.2 - 0.05 * 0.2 * 0.5;
    const expected_common_capacity = 0.4 * toqck_g_c_per_step *
        temperature_factor * (1 - expected_next_inhibitor);
    const expected_hydrolysis = @min(broadcast_urea_mol_n, expected_common_capacity * expected_fraction);
    const result = try ureaHydrolysis(
        .{
            .broadcast_urea_mol_n = broadcast_urea_mol_n,
            .banded_urea_mol_n = 0,
            .soil_mass_megagrams = soil_mass_megagrams,
            .water_volume_m3 = water_volume_m3,
            .biologically_active_water_volume_m3 = volq_m3,
            .total_microbial_respiration_activity_g_c_per_step = toqck_g_c_per_step,
            .temperature_response = temperature_factor,
            .initial_inhibitor_activity = 0.4,
            .current_inhibitor_activity = 0.2,
            .timestep_h = 1,
        },
        .{
            .minimum_half_saturation_mol_n_per_megagram = 0.2,
            .microbial_activity_inhibition_g_c_per_m3_h = 10.0,
            .specific_hydrolysis_mol_n_per_g_c_h = 0.4,
            .inhibitor_decline_rate_per_h = 0.05,
            .negligible_amount = 1e-12,
        },
    );
    try std.testing.expectApproxEqAbs(expected_hydrolysis, result.broadcast_hydrolysis_mol_n, 1e-12);
}

test "SOLUTE fertilizer publication conserves nitrogen and is transactional" {
    var state: FertilizerState = .{
        .broadcast_ammonium_mol_n = 1,
        .broadcast_ammonia_mol_n = 2,
        .broadcast_urea_mol_n = 3,
        .broadcast_nitrate_mol_n = 4,
        .banded_ammonium_mol_n = 5,
        .banded_ammonia_mol_n = 6,
        .banded_urea_mol_n = 7,
        .banded_nitrate_mol_n = 8,
    };
    var aqueous = std.mem.zeroes(aqueous_network.State);
    const fractions: ZoneFractions = .{
        .ammonium_non_band = 0.75,
        .ammonium_band = 0.25,
        .nitrate_non_band = 0.6,
        .nitrate_band = 0.4,
    };
    const flux: DissolutionFlux = .{
        .broadcast_ammonium_non_band_mol_n = 0.1,
        .broadcast_ammonium_band_mol_n = 0.2,
        .broadcast_ammonia_non_band_mol_n = 0.3,
        .broadcast_ammonia_band_mol_n = 0.4,
        .broadcast_urea_non_band_mol_n = 0.5,
        .broadcast_urea_band_mol_n = 0.6,
        .broadcast_nitrate_non_band_mol_n = 0.7,
        .broadcast_nitrate_band_mol_n = 0.8,
        .banded_ammonium_mol_n = 0.9,
        .banded_ammonia_mol_n = 1,
        .banded_urea_mol_n = 1.1,
        .banded_nitrate_mol_n = 1.2,
    };
    const before = sumState(state);
    var gaseous_ammonia_g_n: f64 = 0;
    try commitToRecipients(&state, &aqueous, &gaseous_ammonia_g_n, flux, fractions, 2, 14);
    const aqueous_n =
        aqueous.ammonium_non_band * 1.5 +
        aqueous.ammonium_band * 0.5 +
        aqueous.ammonia_non_band * 1.5 +
        aqueous.ammonia_band * 0.5 +
        aqueous.nitrate_non_band * 1.2 +
        aqueous.nitrate_band * 0.8;
    try std.testing.expectApproxEqAbs(
        before,
        sumState(state) + aqueous_n + gaseous_ammonia_g_n / 14,
        1e-13,
    );

    const donor_before_failure = state;
    const aqueous_before_failure = aqueous;
    const gaseous_before_failure = gaseous_ammonia_g_n;
    var invalid_flux = std.mem.zeroes(DissolutionFlux);
    invalid_flux.banded_ammonium_mol_n = 0.1;
    try std.testing.expectError(
        error.MissingFertilizerRecipientWaterVolume,
        commitToRecipients(
            &state,
            &aqueous,
            &gaseous_ammonia_g_n,
            invalid_flux,
            .{
                .ammonium_non_band = 1,
                .ammonium_band = 0,
                .nitrate_non_band = 1,
                .nitrate_band = 0,
            },
            2,
            14,
        ),
    );
    try std.testing.expectEqual(donor_before_failure, state);
    try std.testing.expectEqual(aqueous_before_failure, aqueous);
    try std.testing.expectEqual(gaseous_before_failure, gaseous_ammonia_g_n);
}

fn sumState(state: FertilizerState) f64 {
    var total: f64 = 0;
    inline for (@typeInfo(FertilizerState).@"struct".fields) |field| total += @field(state, field.name);
    return total;
}

fn sumFlux(flux: DissolutionFlux) f64 {
    var total: f64 = 0;
    inline for (@typeInfo(DissolutionFlux).@"struct".fields) |field| total += @field(flux, field.name);
    return total;
}
