const std = @import("std");
const delimited_input = @import("delimited_input.zig");
const nitrification = @import("soil_nitrification.zig");
const denitrification = @import("soil_denitrification.zig");
const chemodenitrification = @import("soil_chemodenitrification.zig");

pub const AutotrophicDenitrificationParameters = struct {
    anaerobic_growth_respiration_fraction: f64,
    additional_ammonium_oxidation_per_nitrite_reduction: f64,
};

pub const NitrifierIndices = struct {
    autotrophic_substrate_index: usize,
    ammonia_oxidizer_population_index: usize,
    nitrite_oxidizer_population_index: usize,
    heterotrophic_denitrifier_population_index: usize,
};

pub const NitrifierEnvironmentParameters = struct {
    labile_biomass_fraction: f64,
    ammonia_oxidizer_target_nitrogen_per_carbon_g_n_per_g_c: f64,
    nitrite_oxidizer_target_nitrogen_per_carbon_g_n_per_g_c: f64,
    ammonia_oxidizer_target_phosphorus_per_carbon_g_p_per_g_c: f64,
    nitrite_oxidizer_target_phosphorus_per_carbon_g_p_per_g_c: f64,
    aqueous_co2_half_saturation_g_c_per_m3: f64,
    water_potential_sensitivity_per_mpa: f64,
};

pub const OxygenUptakeParameters = struct {
    microbial_radius_m: f64,
    microbial_count_per_g_c: f64,
    oxygen_half_saturation_g_o_per_m3: f64,
    hygroscopic_water_potential_mpa: f64,
    air_water_exchange_reference_time_h: f64,
    wet_exchange_exponent: f64,
    dry_exchange_exponent: f64,
    minimum_transition_water_fraction: f64,
    aqueous_tortuosity_coefficient: f64,
    minimum_allocation_fraction: f64,
    negligible_oxygen_demand_g_o: f64,
};

pub const HeterotrophicRespirationParameters = struct {
    substrate_unlimited_respiration_per_h: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
    doc_half_saturation_g_c_per_m3: f64,
    acetate_half_saturation_g_c_per_m3: f64,
    doc_respiration_requirement_g_c_per_g_c: f64,
    acetate_respiration_requirement_g_c_per_g_c: f64,
    water_potential_sensitivity_per_mpa: f64,
    oxygen_per_respired_carbon_g_o_per_g_c: f64,
    specific_maintenance_respiration_g_c_per_g_n_per_h: f64,
    decomposition_density_half_saturation_g_c_per_g_c: f64,
    maintenance_density_half_saturation_g_c_per_g_c: f64,
    acidity_half_response_mol_per_m3: f64,
    denitrification_growth_respiration_fraction_g_c_per_g_c: f64,
};

pub const MicrobialMineralExchangeParameters = struct {
    ammonium_maximum_uptake_g_n_per_m2_h: f64,
    ammonium_minimum_concentration_g_n_per_m3: f64,
    ammonium_half_saturation_g_n_per_m3: f64,
    nitrate_maximum_uptake_g_n_per_m2_h: f64,
    nitrate_minimum_concentration_g_n_per_m3: f64,
    nitrate_half_saturation_g_n_per_m3: f64,
    phosphate_maximum_uptake_g_p_per_m2_h: f64,
    phosphate_minimum_concentration_g_p_per_m3: f64,
    phosphate_half_saturation_g_p_per_m3: f64,
    phosphorus_molar_mass_g_per_mol: f64,
};

pub const NonsymbioticNitrogenFixationParameters = struct {
    aerobic_diazotroph_population_index: usize,
    anaerobic_diazotroph_population_index: usize,
    aerobic_yield_g_n_per_g_c: f64,
    anaerobic_yield_g_n_per_g_c: f64,
    dinitrogen_half_saturation_g_n_per_m3: f64,
    nonstructural_to_structural_rate_per_h: f64,
};

pub const FermenterRespirationParameters = struct {
    fermenter_population_index: usize,
    anaerobic_diazotroph_population_index: usize,
    specific_oxidation_rate_g_c_per_g_c_h: f64,
    dissolved_organic_carbon_half_saturation_g_c_per_m3: f64,
    acetate_product_inhibition_g_c_per_m3: f64,
    reference_energy_yield_kj_per_g_c: f64,
    growth_energy_requirement_kj_per_g_c: f64,
    diazotroph_growth_energy_requirement_kj_per_g_c: f64,
    minimum_respiration_requirement_g_c_per_g_c: f64,
    diazotroph_minimum_respiration_requirement_g_c_per_g_c: f64,
    gas_constant_kj_per_mol_k: f64,
    acetate_feedback_stoichiometric_exponent: f64,
    feedback_carbon_conversion_g_c_per_mol: f64,
};

pub const AcetotrophicMethanogenesisParameters = struct {
    population_index: usize,
    acetate_product_inhibition_g_c_per_m3: f64,
    acetate_half_saturation_g_c_per_m3: f64,
    specific_respiration_rate_g_c_per_g_c_h: f64,
    reference_energy_yield_kj_per_g_c: f64,
    growth_energy_requirement_kj_per_g_c: f64,
    minimum_growth_respiration_fraction: f64,
    gas_constant_kj_per_mol_k: f64,
    feedback_carbon_conversion_g_c_per_mol: f64,
    methane_carbon_yield_g_c_per_g_c_oxidized: f64,
};

pub const AnaerobicEnergyParameters = struct {
    fermenter: FermenterRespirationParameters,
    acetotrophic_methanogenesis: AcetotrophicMethanogenesisParameters,
};

pub const MicrobialTurnoverParameters = struct {
    labile_basal_decomposition_rate_per_h: f64,
    resistant_basal_decomposition_rate_per_h: f64,
    minimum_carbon_recycling_fraction: f64,
    carbon_recycling_range_fraction: f64,
    maximum_nitrogen_recycling_fraction: f64,
    maximum_phosphorus_recycling_fraction: f64,
    humification_intercept: f64,
    humification_clay_coefficient: f64,
    humification_maximum_clay_fraction: f64,
    woody_colonization_per_g_respired_carbon: f64,
    fine_litter_colonization_per_g_respired_carbon: f64,
    manure_colonization_per_g_respired_carbon: f64,
    particulate_colonization_per_g_respired_carbon: f64,
    humus_colonization_per_g_respired_carbon: f64,
};

pub const MethaneParameters = struct {
    hydrogenotroph_population_index: usize,
    methanotroph_population_index: usize,
    hydrogen_product_inhibition_g_h_per_m3: f64,
    hydrogen_half_saturation_g_h_per_m3: f64,
    hydrogenotrophic_specific_co2_reduction_g_c_per_g_c_h: f64,
    hydrogenotrophic_reference_energy_yield_kj_per_g_c: f64,
    methanogen_growth_energy_requirement_kj_per_g_c: f64,
    minimum_growth_respiration_fraction: f64,
    hydrogen_supply_conversion_g_c_per_g_h: f64,
    fermentation_hydrogen_to_pool_fraction: f64,
    methane_half_saturation_g_c_per_m3: f64,
    methane_solubility_water_to_air: f64,
    gas_exchange_rate_per_step: f64,
    methanotroph_biomass_conversion_efficiency_g_c_per_g_c: f64,
    methanotroph_growth_respiration_g_c_per_g_c: f64,
    methanotroph_specific_oxidation_per_h: f64,
};

pub const Parameters = struct {
    nitrification: nitrification.Parameters,
    denitrification: denitrification.Parameters,
    autotrophic_denitrification: AutotrophicDenitrificationParameters,
    chemodenitrification: chemodenitrification.Parameters,
    nitrous_acid_dissociation_mol_per_m3: f64,
    microbial_thermal_adaptation_offset_k: f64,
    nitrifier_indices: NitrifierIndices,
    nitrifier_environment: NitrifierEnvironmentParameters,
    oxygen_uptake: OxygenUptakeParameters,
    heterotrophic_respiration: HeterotrophicRespirationParameters,
    microbial_mineral_exchange: MicrobialMineralExchangeParameters,
    nonsymbiotic_nitrogen_fixation: NonsymbioticNitrogenFixationParameters,
    microbial_turnover: MicrobialTurnoverParameters,
    microbial_layer_mixing_rate_per_h: f64 = 0.001,
    methane: ?MethaneParameters = null,
};

/// Runtime compatibility defaults transcribed from the source NITRO parameter
/// set. They are returned as ordinary runtime state and remain fully
/// replaceable by a user-supplied parameter file.
pub fn sourceParameters() !Parameters {
    return parse(source_parameter_text);
}

/// Runtime compatibility defaults for NITRO.F 900--1048. This isolated
/// parser is intentionally not production-bound until the two kernels have
/// authoritative state-array owners.
pub fn sourceAnaerobicEnergyParameters() !AnaerobicEnergyParameters {
    return parseAnaerobicEnergyParameters(source_anaerobic_energy_parameter_text);
}

const source_anaerobic_energy_parameter_text =
    "soil_fermenter_respiration,3,6,0.125,12,12,3,37.5,37.5,0.4,0.5,0.0083143,2,72\n" ++
    "soil_acetotrophic_methanogenesis,4,12,12,0.125,1.5,37.5,0.4,0.0083143,24,0.5";

const source_parameter_text =
    "soil_nitrification,0.001,0.0002,7000,14,1.4,1.4,0.125,0.125,0.3,0.1,0.5,2.667,3.429,1.143\n" ++
    "soil_denitrification,0.001,1.4,1.4,0.014,1,0.429,0.429,0.214,0.875\n" ++
    "soil_autotrophic_denitrification,0.5,0.333\n" ++
    "soil_chemodenitrification,0.0005,0.001,1e-12,0.5,0,0.5\n" ++
    "nitrous_acid_dissociation_mol_per_m3,0.45\n" ++
    "soil_microbial_thermal_adaptation_offset_k,0\n" ++
    "soil_nitrifier_indices,5,0,1,1\n" ++
    "soil_nitrifier_environment,0.55,0.1,0.1,0.01,0.01,12,0.1\n" ++
    "soil_oxygen_uptake,1e-6,2.3866348449e11,0.064,-1.5e4,0.5,12,12,0.5,0.7,0.001,1e-12\n" ++
    "soil_heterotrophic_respiration,0.125,0.1,0.01,12,12,0.5,0.42016806722689076,0.1,2.667,0.01,0.01,1e-6,1,0.7142857142857143\n" ++
    "soil_microbial_mineral_exchange,0.014,0.0125,0.40,0.014,0.03,0.35,0.003,0.009,0.18,31\n" ++
    "soil_nonsymbiotic_nitrogen_fixation,5,6,0.25,0.02,0.14,0.25\n" ++
    "soil_microbial_turnover,0.01,0.001,0.167,0.333,0.333,0.333,0.150,0.300,0.333,0.25,2.0,5.0,1.0,0.5\n" ++
    "soil_microbial_layer_mixing,0.001\n" ++
    "soil_methane,4,2,0.001,0.01,0.1,1,10,0.05,1.5,0.111,0.2,0.03,0.5,0.4,0.5,0.125";

/// Runtime soil NITRO parameters. Record names are ASCII case-insensitive and
/// values accept comma, tab, space, or pipe delimiters with any file extension.
pub fn parse(source: []const u8) !Parameters {
    try validateRecordArities(source);
    var tokens = delimited_input.tokens(source);
    var result: Parameters = undefined;
    try expectRecord(&tokens, "soil_nitrification");
    try fillStruct(nitrification.Parameters, &result.nitrification, &tokens);
    try expectRecord(&tokens, "soil_denitrification");
    try fillStruct(denitrification.Parameters, &result.denitrification, &tokens);
    try expectRecord(&tokens, "soil_autotrophic_denitrification");
    try fillStruct(AutotrophicDenitrificationParameters, &result.autotrophic_denitrification, &tokens);
    try expectRecord(&tokens, "soil_chemodenitrification");
    try fillStruct(chemodenitrification.Parameters, &result.chemodenitrification, &tokens);
    try expectRecord(&tokens, "nitrous_acid_dissociation_mol_per_m3");
    result.nitrous_acid_dissociation_mol_per_m3 = try nextFloat(&tokens);
    try expectRecord(&tokens, "soil_microbial_thermal_adaptation_offset_k");
    result.microbial_thermal_adaptation_offset_k = try nextFloat(&tokens);
    try expectRecord(&tokens, "soil_nitrifier_indices");
    result.nitrifier_indices = .{ .autotrophic_substrate_index = try nextUnsigned(&tokens), .ammonia_oxidizer_population_index = try nextUnsigned(&tokens), .nitrite_oxidizer_population_index = try nextUnsigned(&tokens), .heterotrophic_denitrifier_population_index = try nextUnsigned(&tokens) };
    try expectRecord(&tokens, "soil_nitrifier_environment");
    try fillStruct(NitrifierEnvironmentParameters, &result.nitrifier_environment, &tokens);
    try expectRecord(&tokens, "soil_oxygen_uptake");
    try fillStruct(OxygenUptakeParameters, &result.oxygen_uptake, &tokens);
    try expectRecord(&tokens, "soil_heterotrophic_respiration");
    try fillStruct(HeterotrophicRespirationParameters, &result.heterotrophic_respiration, &tokens);
    try expectRecord(&tokens, "soil_microbial_mineral_exchange");
    try fillStruct(MicrobialMineralExchangeParameters, &result.microbial_mineral_exchange, &tokens);
    try expectRecord(&tokens, "soil_nonsymbiotic_nitrogen_fixation");
    result.nonsymbiotic_nitrogen_fixation = .{
        .aerobic_diazotroph_population_index = try nextUnsigned(&tokens),
        .anaerobic_diazotroph_population_index = try nextUnsigned(&tokens),
        .aerobic_yield_g_n_per_g_c = try nextFloat(&tokens),
        .anaerobic_yield_g_n_per_g_c = try nextFloat(&tokens),
        .dinitrogen_half_saturation_g_n_per_m3 = try nextFloat(&tokens),
        .nonstructural_to_structural_rate_per_h = try nextFloat(&tokens),
    };
    try expectRecord(&tokens, "soil_microbial_turnover");
    try fillStruct(MicrobialTurnoverParameters, &result.microbial_turnover, &tokens);
    result.microbial_layer_mixing_rate_per_h = 0.001;
    result.methane = null;
    var has_mixing = false;
    while (tokens.next()) |record| {
        if (std.ascii.eqlIgnoreCase(record, "soil_microbial_layer_mixing")) {
            if (has_mixing) return error.DuplicateSoilNitrogenParameter;
            result.microbial_layer_mixing_rate_per_h = try nextFloat(&tokens);
            has_mixing = true;
        } else if (std.ascii.eqlIgnoreCase(record, "soil_methane")) {
            if (result.methane != null) return error.DuplicateSoilNitrogenParameter;
            var methane: MethaneParameters = undefined;
            methane.hydrogenotroph_population_index = try nextUnsigned(&tokens);
            methane.methanotroph_population_index = try nextUnsigned(&tokens);
            inline for (@typeInfo(MethaneParameters).@"struct".fields[2..]) |field| @field(methane, field.name) = try nextFloat(&tokens);
            result.methane = methane;
        } else return error.UnexpectedSoilNitrogenParameter;
    }
    try validate(result);
    return result;
}

/// Parses the two compulsory anaerobic-energy records with exact physical
/// line arity. Record names are ASCII case-insensitive.
pub fn parseAnaerobicEnergyParameters(source: []const u8) !AnaerobicEnergyParameters {
    try validateAnaerobicEnergyRecordArities(source);
    var tokens = delimited_input.tokens(source);
    try expectRecord(&tokens, "soil_fermenter_respiration");
    const fermenter = FermenterRespirationParameters{
        .fermenter_population_index = try nextUnsigned(&tokens),
        .anaerobic_diazotroph_population_index = try nextUnsigned(&tokens),
        .specific_oxidation_rate_g_c_per_g_c_h = try nextFloat(&tokens),
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = try nextFloat(&tokens),
        .acetate_product_inhibition_g_c_per_m3 = try nextFloat(&tokens),
        .reference_energy_yield_kj_per_g_c = try nextFloat(&tokens),
        .growth_energy_requirement_kj_per_g_c = try nextFloat(&tokens),
        .diazotroph_growth_energy_requirement_kj_per_g_c = try nextFloat(&tokens),
        .minimum_respiration_requirement_g_c_per_g_c = try nextFloat(&tokens),
        .diazotroph_minimum_respiration_requirement_g_c_per_g_c = try nextFloat(&tokens),
        .gas_constant_kj_per_mol_k = try nextFloat(&tokens),
        .acetate_feedback_stoichiometric_exponent = try nextFloat(&tokens),
        .feedback_carbon_conversion_g_c_per_mol = try nextFloat(&tokens),
    };
    try expectRecord(&tokens, "soil_acetotrophic_methanogenesis");
    const acetotrophic = AcetotrophicMethanogenesisParameters{
        .population_index = try nextUnsigned(&tokens),
        .acetate_product_inhibition_g_c_per_m3 = try nextFloat(&tokens),
        .acetate_half_saturation_g_c_per_m3 = try nextFloat(&tokens),
        .specific_respiration_rate_g_c_per_g_c_h = try nextFloat(&tokens),
        .reference_energy_yield_kj_per_g_c = try nextFloat(&tokens),
        .growth_energy_requirement_kj_per_g_c = try nextFloat(&tokens),
        .minimum_growth_respiration_fraction = try nextFloat(&tokens),
        .gas_constant_kj_per_mol_k = try nextFloat(&tokens),
        .feedback_carbon_conversion_g_c_per_mol = try nextFloat(&tokens),
        .methane_carbon_yield_g_c_per_g_c_oxidized = try nextFloat(&tokens),
    };
    if (tokens.next() != null) return error.UnexpectedSoilNitrogenParameter;
    const result = AnaerobicEnergyParameters{
        .fermenter = fermenter,
        .acetotrophic_methanogenesis = acetotrophic,
    };
    try validateAnaerobicEnergyParameters(result);
    return result;
}

fn validateAnaerobicEnergyRecordArities(source: []const u8) !void {
    var records = delimited_input.records(source);
    var record_index: usize = 0;
    const names = [_][]const u8{
        "soil_fermenter_respiration",
        "soil_acetotrophic_methanogenesis",
    };
    const counts = [_]usize{
        structFieldCount(FermenterRespirationParameters),
        structFieldCount(AcetotrophicMethanogenesisParameters),
    };
    while (records.next()) |record| : (record_index += 1) {
        if (record_index >= names.len) return error.UnexpectedSoilNitrogenParameterRecord;
        if (hasEmptyExplicitField(record)) return error.EmptySoilNitrogenParameterValue;
        var tokens = delimited_input.recordTokens(record);
        const label = tokens.next() orelse continue;
        if (!std.ascii.eqlIgnoreCase(label, names[record_index]))
            return error.UnexpectedSoilNitrogenParameterRecord;
        var actual_count: usize = 0;
        while (tokens.next() != null) actual_count += 1;
        if (actual_count != counts[record_index])
            return error.InvalidSoilNitrogenParameterRecordArity;
    }
    if (record_index != names.len) return error.MissingSoilNitrogenParameterRecord;
}

fn validateAnaerobicEnergyParameters(parameters: AnaerobicEnergyParameters) !void {
    const fermenter = parameters.fermenter;
    if (fermenter.fermenter_population_index == fermenter.anaerobic_diazotroph_population_index)
        return error.InvalidSoilNitrogenParameter;
    inline for (@typeInfo(FermenterRespirationParameters).@"struct".fields[2..]) |field| {
        const value = @field(fermenter, field.name);
        if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilNitrogenParameter;
    }
    if (fermenter.minimum_respiration_requirement_g_c_per_g_c > 1 or
        fermenter.diazotroph_minimum_respiration_requirement_g_c_per_g_c > 1)
        return error.InvalidSoilNitrogenParameter;
    const acetotrophic = parameters.acetotrophic_methanogenesis;
    inline for (@typeInfo(AcetotrophicMethanogenesisParameters).@"struct".fields[1..]) |field| {
        const value = @field(acetotrophic, field.name);
        if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilNitrogenParameter;
    }
    if (acetotrophic.minimum_growth_respiration_fraction > 1 or
        acetotrophic.methane_carbon_yield_g_c_per_g_c_oxidized > 1)
        return error.InvalidSoilNitrogenParameter;
}

fn validateRecordArities(source: []const u8) !void {
    var records = delimited_input.records(source);
    while (records.next()) |record| {
        if (hasEmptyExplicitField(record))
            return error.EmptySoilNitrogenParameterValue;
        var tokens = delimited_input.recordTokens(record);
        const label = tokens.next() orelse continue;
        if (label[0] == '#') continue;
        const expected_count = expectedRecordValueCount(label) orelse
            return error.UnexpectedSoilNitrogenParameterRecord;
        var actual_count: usize = 0;
        while (tokens.next() != null) actual_count += 1;
        if (actual_count != expected_count)
            return error.InvalidSoilNitrogenParameterRecordArity;
    }
}

fn hasEmptyExplicitField(line: []const u8) bool {
    const content = if (std.mem.indexOfScalar(u8, line, '#')) |comment| line[0..comment] else line;
    const trimmed = std.mem.trim(u8, content, " \r");
    if (trimmed.len == 0) return false;
    var field_start: usize = 0;
    var saw_explicit_delimiter = false;
    for (trimmed, 0..) |byte, index| {
        if (byte != ',' and byte != '|' and byte != '\t') continue;
        const field = std.mem.trim(u8, trimmed[field_start..index], " \r");
        if (field.len == 0) return true;
        field_start = index + 1;
        saw_explicit_delimiter = true;
    }
    return saw_explicit_delimiter and std.mem.trim(u8, trimmed[field_start..], " \r").len == 0;
}

fn expectedRecordValueCount(label: []const u8) ?usize {
    const Entry = struct {
        name: []const u8,
        count: usize,
    };
    const entries = [_]Entry{
        .{ .name = "soil_nitrification", .count = structFieldCount(nitrification.Parameters) },
        .{ .name = "soil_denitrification", .count = structFieldCount(denitrification.Parameters) },
        .{ .name = "soil_autotrophic_denitrification", .count = structFieldCount(AutotrophicDenitrificationParameters) },
        .{ .name = "soil_chemodenitrification", .count = structFieldCount(chemodenitrification.Parameters) },
        .{ .name = "nitrous_acid_dissociation_mol_per_m3", .count = 1 },
        .{ .name = "soil_microbial_thermal_adaptation_offset_k", .count = 1 },
        .{ .name = "soil_nitrifier_indices", .count = structFieldCount(NitrifierIndices) },
        .{ .name = "soil_nitrifier_environment", .count = structFieldCount(NitrifierEnvironmentParameters) },
        .{ .name = "soil_oxygen_uptake", .count = structFieldCount(OxygenUptakeParameters) },
        .{ .name = "soil_heterotrophic_respiration", .count = structFieldCount(HeterotrophicRespirationParameters) },
        .{ .name = "soil_microbial_mineral_exchange", .count = structFieldCount(MicrobialMineralExchangeParameters) },
        .{ .name = "soil_nonsymbiotic_nitrogen_fixation", .count = structFieldCount(NonsymbioticNitrogenFixationParameters) },
        .{ .name = "soil_microbial_turnover", .count = structFieldCount(MicrobialTurnoverParameters) },
        .{ .name = "soil_microbial_layer_mixing", .count = 1 },
        .{ .name = "soil_methane", .count = structFieldCount(MethaneParameters) },
    };
    for (entries) |entry| {
        if (std.ascii.eqlIgnoreCase(label, entry.name)) return entry.count;
    }
    return null;
}

fn structFieldCount(comptime T: type) usize {
    return @typeInfo(T).@"struct".fields.len;
}

pub fn validate(parameters: Parameters) !void {
    try parameters.nitrification.validate();
    try parameters.denitrification.validate();
    const autotrophic = parameters.autotrophic_denitrification;
    if (!std.math.isFinite(autotrophic.anaerobic_growth_respiration_fraction) or autotrophic.anaerobic_growth_respiration_fraction < 0 or autotrophic.anaerobic_growth_respiration_fraction > 1 or !std.math.isFinite(autotrophic.additional_ammonium_oxidation_per_nitrite_reduction) or autotrophic.additional_ammonium_oxidation_per_nitrite_reduction <= 0 or autotrophic.additional_ammonium_oxidation_per_nitrite_reduction >= 1) return error.InvalidSoilNitrogenParameter;
    const chemo = parameters.chemodenitrification;
    inline for (@typeInfo(chemodenitrification.Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(chemo, field.name)) or @field(chemo, field.name) < 0) return error.InvalidSoilNitrogenParameter;
    if (chemo.negligible_demand_g_n <= 0 or @abs(chemo.nitrous_oxide_product_fraction + chemo.dinitrogen_product_fraction + chemo.dissolved_organic_nitrogen_product_fraction - 1) > 1e-12) return error.InvalidSoilNitrogenParameter;
    if (!std.math.isFinite(parameters.nitrous_acid_dissociation_mol_per_m3) or parameters.nitrous_acid_dissociation_mol_per_m3 <= 0) return error.InvalidSoilNitrogenParameter;
    if (!std.math.isFinite(parameters.microbial_thermal_adaptation_offset_k)) return error.InvalidSoilNitrogenParameter;
    if (parameters.nitrifier_indices.ammonia_oxidizer_population_index == parameters.nitrifier_indices.nitrite_oxidizer_population_index) return error.InvalidSoilNitrogenParameter;
    inline for (@typeInfo(NitrifierEnvironmentParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters.nitrifier_environment, field.name)) or @field(parameters.nitrifier_environment, field.name) <= 0) return error.InvalidSoilNitrogenParameter;
    if (parameters.nitrifier_environment.labile_biomass_fraction > 1) return error.InvalidSoilNitrogenParameter;
    inline for (@typeInfo(OxygenUptakeParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters.oxygen_uptake, field.name))) return error.InvalidSoilNitrogenParameter;
    const oxygen = parameters.oxygen_uptake;
    if (oxygen.microbial_radius_m <= 0 or oxygen.microbial_count_per_g_c < 0 or oxygen.oxygen_half_saturation_g_o_per_m3 <= 0 or oxygen.hygroscopic_water_potential_mpa >= 0 or oxygen.air_water_exchange_reference_time_h <= 0 or oxygen.minimum_transition_water_fraction < 0 or oxygen.minimum_transition_water_fraction > 1 or oxygen.aqueous_tortuosity_coefficient < 0 or oxygen.minimum_allocation_fraction < 0 or oxygen.negligible_oxygen_demand_g_o < 0) return error.InvalidSoilNitrogenParameter;
    inline for (@typeInfo(HeterotrophicRespirationParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters.heterotrophic_respiration, field.name)) or @field(parameters.heterotrophic_respiration, field.name) <= 0) return error.InvalidSoilNitrogenParameter;
    inline for (@typeInfo(MicrobialMineralExchangeParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters.microbial_mineral_exchange, field.name)) or @field(parameters.microbial_mineral_exchange, field.name) <= 0) return error.InvalidSoilNitrogenParameter;
    const fixation = parameters.nonsymbiotic_nitrogen_fixation;
    if (fixation.aerobic_diazotroph_population_index == fixation.anaerobic_diazotroph_population_index) return error.InvalidSoilNitrogenParameter;
    inline for (.{ fixation.aerobic_yield_g_n_per_g_c, fixation.anaerobic_yield_g_n_per_g_c, fixation.dinitrogen_half_saturation_g_n_per_m3, fixation.nonstructural_to_structural_rate_per_h }) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilNitrogenParameter;
    const turnover = parameters.microbial_turnover;
    inline for (@typeInfo(MicrobialTurnoverParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(turnover, field.name)) or @field(turnover, field.name) < 0) return error.InvalidSoilNitrogenParameter;
    if (!std.math.isFinite(parameters.microbial_layer_mixing_rate_per_h) or parameters.microbial_layer_mixing_rate_per_h < 0) return error.InvalidSoilNitrogenParameter;
    if (parameters.methane) |methane| {
        if (methane.hydrogenotroph_population_index == methane.methanotroph_population_index) return error.InvalidSoilNitrogenParameter;
        inline for (@typeInfo(MethaneParameters).@"struct".fields[2..]) |field| if (!std.math.isFinite(@field(methane, field.name)) or @field(methane, field.name) < 0) return error.InvalidSoilNitrogenParameter;
        if (methane.hydrogen_half_saturation_g_h_per_m3 <= 0 or methane.methanogen_growth_energy_requirement_kj_per_g_c <= 0 or methane.hydrogen_supply_conversion_g_c_per_g_h <= 0 or methane.methane_half_saturation_g_c_per_m3 <= 0 or methane.minimum_growth_respiration_fraction > 1) return error.InvalidSoilNitrogenParameter;
    }
    if (turnover.minimum_carbon_recycling_fraction + turnover.carbon_recycling_range_fraction > 1 or turnover.maximum_nitrogen_recycling_fraction > 1 or turnover.maximum_phosphorus_recycling_fraction > 1 or turnover.humification_maximum_clay_fraction > 1) return error.InvalidSoilNitrogenParameter;
}

fn fillStruct(comptime T: type, value: *T, tokens: anytype) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| @field(value, field.name) = try nextFloat(tokens);
}

fn expectRecord(tokens: anytype, expected: []const u8) !void {
    const actual = tokens.next() orelse return error.MissingSoilNitrogenParameterRecord;
    if (!std.ascii.eqlIgnoreCase(actual, expected)) return error.UnexpectedSoilNitrogenParameterRecord;
}

fn nextFloat(tokens: anytype) !f64 {
    const text = tokens.next() orelse return error.MissingSoilNitrogenParameter;
    const value = std.fmt.parseFloat(f64, text) catch return error.InvalidSoilNitrogenParameter;
    if (!std.math.isFinite(value)) return error.NonFiniteSoilNitrogenParameter;
    return value;
}

fn nextUnsigned(tokens: anytype) !usize {
    const text = tokens.next() orelse return error.MissingSoilNitrogenParameter;
    return std.fmt.parseUnsigned(usize, text, 10) catch return error.InvalidSoilNitrogenParameter;
}

test "soil nitrogen parameters are self-contained flexible and case-insensitive" {
    const source =
        "SOIL_NITRIFICATION,0.001,0.0002,7000,14,1.4,1.4,0.125,0.125,0.3,0.1,0.5,2.667,3.429,1.143\n" ++
        "Soil_Denitrification|0.001|1.4|1.4|0.014|1|0.429|0.429|0.214|0.875\n" ++
        "soil_autotrophic_denitrification 0.5 0.333\n" ++
        "soil_chemodenitrification\t0.0005\t0.001\t1e-12\t0.5\t0\t0.5\n" ++
        "nitrous_acid_dissociation_mol_per_m3 0.45\n" ++
        "soil_microbial_thermal_adaptation_offset_k 0\n" ++
        "soil_nitrifier_indices 5 0 1 1\n" ++
        "soil_nitrifier_environment 0.55 0.1 0.1 0.01 0.01 12 0.1\n" ++
        "soil_oxygen_uptake 1e-6 2.3866348449e11 0.064 -1.5e4 0.5 12 12 0.5 0.7 0.001 1e-12\n" ++
        "soil_heterotrophic_respiration 0.125 0.1 0.01 12 12 0.5 0.42016806722689076 0.1 2.667 0.01 0.01 1e-6 1 0.7142857142857143\n" ++
        "soil_microbial_mineral_exchange 0.014 0.0125 0.40 0.014 0.03 0.35 0.003 0.009 0.18 31\n" ++
        "soil_nonsymbiotic_nitrogen_fixation 5 6 0.25 0.02 0.14 0.25\n" ++
        "soil_microbial_turnover 0.01 0.001 0.167 0.333 0.333 0.333 0.150 0.300 0.333 0.25 2.0 5.0 1.0 0.5\n" ++
        "Soil_Microbial_Layer_Mixing|0.0025\n" ++
        "Soil_Methane|4|2|0.001|0.01|0.1|1|10|0.05|1.5|0.111|0.2|0.03|0.5|0.4|0.5|0.125";
    const parameters = try parse(source);
    try std.testing.expectEqual(@as(f64, 0.125), parameters.nitrification.ammonia_oxidation_rate_g_n_per_g_c_h);
    try std.testing.expectEqual(@as(f64, 0.45), parameters.nitrous_acid_dissociation_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0.064), parameters.oxygen_uptake.oxygen_half_saturation_g_o_per_m3);
    try std.testing.expectEqual(@as(f64, 0.0025), parameters.microbial_layer_mixing_rate_per_h);
    try std.testing.expectEqual(@as(usize, 4), parameters.methane.?.hydrogenotroph_population_index);
    try std.testing.expectEqual(@as(f64, 0.125), parameters.methane.?.methanotroph_specific_oxidation_per_h);
}

test "source NITRO defaults are runtime validated and keep methane active" {
    const parameters = try sourceParameters();
    try validate(parameters);
    try std.testing.expect(parameters.methane != null);
    try std.testing.expectEqual(@as(f64, 0.125), parameters.nitrification.ammonia_oxidation_rate_g_n_per_g_c_h);
    try std.testing.expectEqual(@as(f64, 31), parameters.microbial_mineral_exchange.phosphorus_molar_mass_g_per_mol);
    try std.testing.expectEqual(@as(f64, 0.001), parameters.microbial_layer_mixing_rate_per_h);
}

test "source NITRO derived geometry phosphate and EN2F defaults reach runtime owners" {
    const parameters = try sourceParameters();
    const radius_m = parameters.oxygen_uptake.microbial_radius_m;
    const source_count_per_g_c = 1.0e-6 / (4.19 * radius_m * radius_m * radius_m);
    try std.testing.expectApproxEqRel(
        source_count_per_g_c,
        parameters.oxygen_uptake.microbial_count_per_g_c,
        2e-11,
    );
    const source_surface_area_m2_per_g_c =
        source_count_per_g_c * 12.57 * radius_m * radius_m;
    const kernel_surface_area_m2_per_g_c =
        parameters.oxygen_uptake.microbial_count_per_g_c *
        4 * std.math.pi * radius_m * radius_m;
    // The kernel uses exact 4*pi while Fortran used rounded 12.57.
    try std.testing.expectApproxEqRel(
        source_surface_area_m2_per_g_c,
        kernel_surface_area_m2_per_g_c,
        3e-4,
    );

    const phosphate = parameters.microbial_mineral_exchange;
    try std.testing.expectEqual(
        @as(f64, 0.25 * 0.003),
        0.25 * phosphate.phosphate_maximum_uptake_g_p_per_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.25 * 0.009),
        0.25 * phosphate.phosphate_minimum_concentration_g_p_per_m3,
    );

    var fixation_yield_by_population: [7]f64 = @splat(0);
    const fixation = parameters.nonsymbiotic_nitrogen_fixation;
    fixation_yield_by_population[fixation.aerobic_diazotroph_population_index] =
        fixation.aerobic_yield_g_n_per_g_c;
    fixation_yield_by_population[fixation.anaerobic_diazotroph_population_index] =
        fixation.anaerobic_yield_g_n_per_g_c;
    try std.testing.expectEqual(
        [7]f64{ 0, 0, 0, 0, 0, 0.25, 0.02 },
        fixation_yield_by_population,
    );
    try std.testing.expectEqual(@as(f64, 37.5) / 150.0, fixation.aerobic_yield_g_n_per_g_c);
    try std.testing.expectEqual(@as(f64, 3.0) / 150.0, fixation.anaerobic_yield_g_n_per_g_c);
}

test "anaerobic energy owners remain outside production parameter aggregate" {
    // NITRO EOMF/EOMY/GCHX and EOMH/GC4X reach the isolated kernels only as
    // caller-supplied arrays. Isolated runtime owners now exist, but keep
    // this negative proof until production bindings are added deliberately.
    try std.testing.expect(!@hasField(Parameters, "fermenter_respiration"));
    try std.testing.expect(!@hasField(Parameters, "acetotrophic_methanogenesis"));

    const methane = (try sourceParameters()).methane.?;
    // The existing methane record belongs to the hydrogenotrophic path and
    // must not be mistaken for source acetotrophic GC4X=1.5/EOMH=37.5.
    try std.testing.expect(methane.hydrogenotrophic_reference_energy_yield_kj_per_g_c != 1.5);
    try std.testing.expect(methane.methanogen_growth_energy_requirement_kj_per_g_c != 37.5);
}

test "NITRO anaerobic energy defaults have isolated runtime owners" {
    const parameters = try sourceAnaerobicEnergyParameters();
    const fermenter = parameters.fermenter;
    try std.testing.expectEqual(@as(usize, 3), fermenter.fermenter_population_index);
    try std.testing.expectEqual(@as(usize, 6), fermenter.anaerobic_diazotroph_population_index);
    try std.testing.expectEqual(@as(f64, 3), fermenter.reference_energy_yield_kj_per_g_c);
    try std.testing.expectEqual(@as(f64, 37.5), fermenter.growth_energy_requirement_kj_per_g_c);
    try std.testing.expectEqual(
        1.0 / (1.0 + 37.5 / 25.0),
        fermenter.minimum_respiration_requirement_g_c_per_g_c,
    );
    try std.testing.expectEqual(
        1.0 / (1.0 + 37.5 / 37.5),
        fermenter.diazotroph_minimum_respiration_requirement_g_c_per_g_c,
    );

    const acetotrophic = parameters.acetotrophic_methanogenesis;
    try std.testing.expectEqual(@as(usize, 4), acetotrophic.population_index);
    try std.testing.expectEqual(@as(f64, 1.5), acetotrophic.reference_energy_yield_kj_per_g_c);
    try std.testing.expectEqual(@as(f64, 37.5), acetotrophic.growth_energy_requirement_kj_per_g_c);
    try std.testing.expectEqual(
        1.0 / (1.0 + 37.5 / 25.0),
        acetotrophic.minimum_growth_respiration_fraction,
    );
    try std.testing.expectEqual(@as(f64, 0.5), acetotrophic.methane_carbon_yield_g_c_per_g_c_oxidized);
}

test "anaerobic energy records are strict flexible and case-insensitive" {
    const parameters = try parseAnaerobicEnergyParameters(
        "Soil_Fermenter_Respiration|8|9|0.2|13|14|4|38|39|0.41|0.51|0.008|2.1|73 # runtime roles\n" ++
            "SOIL_ACETOTROPHIC_METHANOGENESIS,10,15,16,0.3,1.6,40,0.42,0.009,25,0.6",
    );
    try std.testing.expectEqual(@as(usize, 8), parameters.fermenter.fermenter_population_index);
    try std.testing.expectEqual(@as(f64, 4), parameters.fermenter.reference_energy_yield_kj_per_g_c);
    try std.testing.expectEqual(@as(usize, 10), parameters.acetotrophic_methanogenesis.population_index);
    try std.testing.expectEqual(@as(f64, 40), parameters.acetotrophic_methanogenesis.growth_energy_requirement_kj_per_g_c);
    try std.testing.expectError(
        error.InvalidSoilNitrogenParameterRecordArity,
        parseAnaerobicEnergyParameters(
            "soil_fermenter_respiration 3 6 0.125\n" ++
                "soil_acetotrophic_methanogenesis 4 12 12 0.125 1.5 37.5 0.4 0.0083143 24 0.5",
        ),
    );
}

test "soil nitrogen records accept comments and retain case-insensitive labels" {
    const parameters = try parse(
        "# Every parameter record has exact physical-line arity.\n" ++
            source_parameter_text ++
            "\n# End of soil nitrogen parameters.\n",
    );
    try std.testing.expect(parameters.methane != null);
    try std.testing.expectEqual(
        @as(f64, 0.001),
        parameters.microbial_layer_mixing_rate_per_h,
    );
}

test "soil nitrogen preflight rejects short and long physical records" {
    try std.testing.expectError(
        error.InvalidSoilNitrogenParameterRecordArity,
        parse("soil_nitrification 0.001 # remaining values cannot come from another line\n"),
    );
    try std.testing.expectError(
        error.InvalidSoilNitrogenParameterRecordArity,
        parse(source_parameter_text ++ " 99 # extra methane value\n"),
    );
}

test "soil nitrogen records reject explicit empty values" {
    try std.testing.expectError(
        error.EmptySoilNitrogenParameterValue,
        parse("soil_nitrification 0.001,,0.0002 7000 14 1.4 1.4 0.125 0.125 0.3 0.1 0.5 2.667 3.429\n"),
    );
}
