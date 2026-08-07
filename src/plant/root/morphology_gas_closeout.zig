const std = @import("std");

pub const State = struct {
    root_length_per_plant_m: []f64,
    root_length_density_m_per_m3: []f64,
    root_air_volume_m3: []f64,
    root_water_volume_m3: []f64,
    primary_radius_m: []f64,
    secondary_radius_m: []f64,
    root_surface_area_per_plant_m2: []f64,
    average_secondary_length_m: []f64,
    gaseous_mass_g_by_domain_layer_species: []f64,
    aqueous_mass_g_by_domain_layer_species: []f64,
    withdrawal_ledger_g_by_species: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    gas_species_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    population_count: f64,
    primary_length_per_plant_m_by_domain_layer: []const f64,
    primary_carbon_g_c_by_domain_layer: []const f64,
    secondary_length_m_by_domain_layer: []const f64,
    secondary_carbon_g_c_by_domain_layer: []const f64,
    nonwoody_carbon_fraction_by_domain: []const f64,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    secondary_cross_section_m2_by_domain: []const f64,
    volume_per_carbon_m3_per_g_c_by_domain: []const f64,
    root_turgor_water_potential_mpa_by_domain_layer: []const f64,
    root_growth_water_potential_mpa_by_domain_layer: []const f64,
    porosity_fraction_by_domain: []const f64,
    minimum_primary_radius_m_by_domain: []const f64,
    minimum_secondary_radius_m_by_domain: []const f64,
    base_primary_radius_m_by_domain: []const f64,
    base_secondary_radius_m_by_domain: []const f64,
    root_elastic_modulus_megapascal: f64,
    secondary_axis_count_by_domain_layer: []const f64,
    minimum_average_secondary_length_m: f64,
    root_profile_type: i32,
    presence_threshold: f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, domain_layers: usize, gas_values: usize, gases: usize) !void {
    inline for (.{ state.root_length_per_plant_m, state.root_length_density_m_per_m3, state.root_air_volume_m3, state.root_water_volume_m3, state.primary_radius_m, state.secondary_radius_m, state.root_surface_area_per_plant_m2, state.average_secondary_length_m }) |values| if (values.len != domain_layers) return error.RootMorphologyDimensionMismatch;
    inline for (.{ state.gaseous_mass_g_by_domain_layer_species, state.aqueous_mass_g_by_domain_layer_species }) |values| if (values.len != gas_values) return error.RootMorphologyDimensionMismatch;
    if (state.withdrawal_ledger_g_by_species.len != gases) return error.RootMorphologyDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| for (@field(state, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootMorphologyState;
}

fn validateInputs(inputs: Inputs, domain_layers: usize) !void {
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.gas_species_count == 0 or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count) return error.RootMorphologyDimensionMismatch;
    inline for (.{ inputs.primary_length_per_plant_m_by_domain_layer, inputs.primary_carbon_g_c_by_domain_layer, inputs.secondary_length_m_by_domain_layer, inputs.secondary_carbon_g_c_by_domain_layer, inputs.root_turgor_water_potential_mpa_by_domain_layer, inputs.root_growth_water_potential_mpa_by_domain_layer, inputs.secondary_axis_count_by_domain_layer }) |values| if (values.len != domain_layers) return error.RootMorphologyDimensionMismatch;
    inline for (.{ inputs.nonwoody_carbon_fraction_by_domain, inputs.secondary_cross_section_m2_by_domain, inputs.volume_per_carbon_m3_per_g_c_by_domain, inputs.porosity_fraction_by_domain, inputs.minimum_primary_radius_m_by_domain, inputs.minimum_secondary_radius_m_by_domain, inputs.base_primary_radius_m_by_domain, inputs.base_secondary_radius_m_by_domain }) |values| if (values.len != inputs.biological_domain_count) return error.RootMorphologyDimensionMismatch;
    if (inputs.layer_thickness_m.len != inputs.soil_layer_count) return error.RootMorphologyDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| switch (field.type) {
        f64 => if (!std.math.isFinite(@field(inputs, field.name))) return error.InvalidRootMorphologyInput,
        []const f64 => for (@field(inputs, field.name)) |value| if (!std.math.isFinite(value)) return error.InvalidRootMorphologyInput,
        else => {},
    };
    if (inputs.population_count < 0 or inputs.minimum_active_layer_thickness_m < 0 or inputs.root_elastic_modulus_megapascal <= 0 or inputs.minimum_average_secondary_length_m < 0 or inputs.presence_threshold < 0) return error.InvalidRootMorphologyInput;
    inline for (.{ inputs.primary_length_per_plant_m_by_domain_layer, inputs.primary_carbon_g_c_by_domain_layer, inputs.secondary_length_m_by_domain_layer, inputs.secondary_carbon_g_c_by_domain_layer, inputs.secondary_axis_count_by_domain_layer, inputs.layer_thickness_m, inputs.nonwoody_carbon_fraction_by_domain, inputs.secondary_cross_section_m2_by_domain, inputs.volume_per_carbon_m3_per_g_c_by_domain, inputs.minimum_primary_radius_m_by_domain, inputs.minimum_secondary_radius_m_by_domain, inputs.base_primary_radius_m_by_domain, inputs.base_secondary_radius_m_by_domain }) |values| for (values) |value| if (value < 0) return error.InvalidRootMorphologyInput;
    for (inputs.nonwoody_carbon_fraction_by_domain) |value| if (value > 1) return error.InvalidRootMorphologyInput;
    for (inputs.porosity_fraction_by_domain) |value| if (value < 0 or value > 1) return error.InvalidRootMorphologyInput;
}

/// Exact GROSUB 7425--7488 live-root morphology publication and empty-root
/// gas closeout. Traversal is N then L; gas species replace the six fixed
/// Fortran arrays while retaining their elemental-gram and signed-ledger units.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.RootMorphologyDimensionOverflow;
    const gas_values = std.math.mul(usize, domain_layers, inputs.gas_species_count) catch return error.RootMorphologyDimensionOverflow;
    try validateState(state, domain_layers, gas_values, inputs.gas_species_count);
    try validateState(workspace, domain_layers, gas_values, inputs.gas_species_count);
    try validateInputs(inputs, domain_layers);
    copyState(workspace, state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        const index = domain * inputs.soil_layer_count + layer;
        const primary_population_length_m = inputs.primary_length_per_plant_m_by_domain_layer[index] * inputs.population_count;
        const total_length_m = inputs.secondary_length_m_by_domain_layer[index] + primary_population_length_m;
        const total_carbon_g_c = inputs.secondary_carbon_g_c_by_domain_layer[index] + inputs.primary_carbon_g_c_by_domain_layer[index];
        if (total_length_m > inputs.presence_threshold and total_carbon_g_c > inputs.presence_threshold and inputs.population_count > inputs.presence_threshold) {
            workspace.root_length_per_plant_m[index] = total_length_m / inputs.population_count * inputs.nonwoody_carbon_fraction_by_domain[domain];
            workspace.root_length_density_m_per_m3[index] = if (inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m) workspace.root_length_per_plant_m[index] / inputs.layer_thickness_m[layer] else 0;
            const volume_m3 = @max(inputs.secondary_cross_section_m2_by_domain[domain] * inputs.secondary_length_m_by_domain_layer[index], inputs.secondary_carbon_g_c_by_domain_layer[index] * inputs.volume_per_carbon_m3_per_g_c_by_domain[domain] * inputs.root_growth_water_potential_mpa_by_domain_layer[index]);
            workspace.root_air_volume_m3[index] = inputs.porosity_fraction_by_domain[domain] * volume_m3;
            workspace.root_water_volume_m3[index] = (1 - inputs.porosity_fraction_by_domain[domain]) * volume_m3;
            workspace.primary_radius_m[index] = @max(inputs.minimum_primary_radius_m_by_domain[domain], (1 + inputs.root_turgor_water_potential_mpa_by_domain_layer[index] / inputs.root_elastic_modulus_megapascal) * inputs.base_primary_radius_m_by_domain[domain]);
            workspace.secondary_radius_m[index] = @max(inputs.minimum_secondary_radius_m_by_domain[domain], (1 + inputs.root_turgor_water_potential_mpa_by_domain_layer[index] / inputs.root_elastic_modulus_megapascal) * inputs.base_secondary_radius_m_by_domain[domain]);
            workspace.average_secondary_length_m[index] = if (inputs.secondary_axis_count_by_domain_layer[index] > inputs.presence_threshold) @max(inputs.minimum_average_secondary_length_m, inputs.secondary_length_m_by_domain_layer[index] / inputs.secondary_axis_count_by_domain_layer[index]) else inputs.minimum_average_secondary_length_m;
            var area_m2 = 6.283 * workspace.primary_radius_m[index] * primary_population_length_m + 6.283 * workspace.secondary_radius_m[index] * inputs.secondary_length_m_by_domain_layer[index];
            if (inputs.root_profile_type != 0) area_m2 *= inputs.minimum_average_secondary_length_m / workspace.average_secondary_length_m[index];
            workspace.root_surface_area_per_plant_m2[index] = area_m2 / inputs.population_count * inputs.nonwoody_carbon_fraction_by_domain[domain];
        } else {
            workspace.root_length_per_plant_m[index] = 0;
            workspace.root_length_density_m_per_m3[index] = 0;
            workspace.root_air_volume_m3[index] = 0;
            workspace.root_water_volume_m3[index] = 0;
            workspace.primary_radius_m[index] = inputs.base_primary_radius_m_by_domain[domain];
            workspace.secondary_radius_m[index] = inputs.base_secondary_radius_m_by_domain[domain];
            workspace.root_surface_area_per_plant_m2[index] = 0;
            workspace.average_secondary_length_m[index] = inputs.minimum_average_secondary_length_m;
            for (0..inputs.gas_species_count) |species| {
                const gas_index = index * inputs.gas_species_count + species;
                workspace.withdrawal_ledger_g_by_species[species] -= workspace.gaseous_mass_g_by_domain_layer_species[gas_index] + workspace.aqueous_mass_g_by_domain_layer_species[gas_index];
                workspace.gaseous_mass_g_by_domain_layer_species[gas_index] = 0;
                workspace.aqueous_mass_g_by_domain_layer_species[gas_index] = 0;
            }
        }
    };
    inline for (@typeInfo(State).@"struct".fields) |field| for (@field(workspace, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootMorphologyResult;
    copyState(state, workspace);
}

test "GROSUB live root morphology preserves source equations" {
    var out: [8][1]f64 = std.mem.zeroes([8][1]f64);
    var gas = [_]f64{0};
    var aq = [_]f64{0};
    var ledger = [_]f64{0};
    var wo: [8][1]f64 = std.mem.zeroes([8][1]f64);
    var wg = [_]f64{0};
    var wa = [_]f64{0};
    var wl = [_]f64{0};
    const state = State{ .root_length_per_plant_m = &out[0], .root_length_density_m_per_m3 = &out[1], .root_air_volume_m3 = &out[2], .root_water_volume_m3 = &out[3], .primary_radius_m = &out[4], .secondary_radius_m = &out[5], .root_surface_area_per_plant_m2 = &out[6], .average_secondary_length_m = &out[7], .gaseous_mass_g_by_domain_layer_species = &gas, .aqueous_mass_g_by_domain_layer_species = &aq, .withdrawal_ledger_g_by_species = &ledger };
    const work = State{ .root_length_per_plant_m = &wo[0], .root_length_density_m_per_m3 = &wo[1], .root_air_volume_m3 = &wo[2], .root_water_volume_m3 = &wo[3], .primary_radius_m = &wo[4], .secondary_radius_m = &wo[5], .root_surface_area_per_plant_m2 = &wo[6], .average_secondary_length_m = &wo[7], .gaseous_mass_g_by_domain_layer_species = &wg, .aqueous_mass_g_by_domain_layer_species = &wa, .withdrawal_ledger_g_by_species = &wl };
    try apply(state, work, .{ .biological_domain_count = 1, .soil_layer_count = 1, .gas_species_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .population_count = 2, .primary_length_per_plant_m_by_domain_layer = &.{3}, .primary_carbon_g_c_by_domain_layer = &.{2}, .secondary_length_m_by_domain_layer = &.{4}, .secondary_carbon_g_c_by_domain_layer = &.{3}, .nonwoody_carbon_fraction_by_domain = &.{0.5}, .layer_thickness_m = &.{0.25}, .minimum_active_layer_thickness_m = 0.01, .secondary_cross_section_m2_by_domain = &.{0.001}, .volume_per_carbon_m3_per_g_c_by_domain = &.{0.002}, .root_turgor_water_potential_mpa_by_domain_layer = &.{0}, .root_growth_water_potential_mpa_by_domain_layer = &.{1}, .porosity_fraction_by_domain = &.{0.25}, .minimum_primary_radius_m_by_domain = &.{0.001}, .minimum_secondary_radius_m_by_domain = &.{0.0005}, .base_primary_radius_m_by_domain = &.{0.002}, .base_secondary_radius_m_by_domain = &.{0.001}, .root_elastic_modulus_megapascal = 10, .secondary_axis_count_by_domain_layer = &.{2}, .minimum_average_secondary_length_m = 0.1, .root_profile_type = 0, .presence_threshold = 1e-12 });
    try std.testing.expectApproxEqAbs(2.5, out[0][0], 1e-12);
    try std.testing.expectApproxEqAbs(10, out[1][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.0015, out[2][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.0045, out[3][0], 1e-12);
    try std.testing.expectApproxEqAbs(2, out[7][0], 1e-12);
}

test "GROSUB empty roots release every runtime gas species conservatively" {
    var out: [8][1]f64 = std.mem.zeroes([8][1]f64);
    var gas = [_]f64{ 1, 2, 3 };
    var aq = [_]f64{ 0.5, 1, 1.5 };
    var ledger = [_]f64{ 0, 0, 0 };
    var wo: [8][1]f64 = std.mem.zeroes([8][1]f64);
    var wg = [_]f64{0} ** 3;
    var wa = wg;
    var wl = wg;
    const state = State{ .root_length_per_plant_m = &out[0], .root_length_density_m_per_m3 = &out[1], .root_air_volume_m3 = &out[2], .root_water_volume_m3 = &out[3], .primary_radius_m = &out[4], .secondary_radius_m = &out[5], .root_surface_area_per_plant_m2 = &out[6], .average_secondary_length_m = &out[7], .gaseous_mass_g_by_domain_layer_species = &gas, .aqueous_mass_g_by_domain_layer_species = &aq, .withdrawal_ledger_g_by_species = &ledger };
    const work = State{ .root_length_per_plant_m = &wo[0], .root_length_density_m_per_m3 = &wo[1], .root_air_volume_m3 = &wo[2], .root_water_volume_m3 = &wo[3], .primary_radius_m = &wo[4], .secondary_radius_m = &wo[5], .root_surface_area_per_plant_m2 = &wo[6], .average_secondary_length_m = &wo[7], .gaseous_mass_g_by_domain_layer_species = &wg, .aqueous_mass_g_by_domain_layer_species = &wa, .withdrawal_ledger_g_by_species = &wl };
    try apply(state, work, .{ .biological_domain_count = 1, .soil_layer_count = 1, .gas_species_count = 3, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .population_count = 1, .primary_length_per_plant_m_by_domain_layer = &.{0}, .primary_carbon_g_c_by_domain_layer = &.{0}, .secondary_length_m_by_domain_layer = &.{0}, .secondary_carbon_g_c_by_domain_layer = &.{0}, .nonwoody_carbon_fraction_by_domain = &.{1}, .layer_thickness_m = &.{0.2}, .minimum_active_layer_thickness_m = 0, .secondary_cross_section_m2_by_domain = &.{0}, .volume_per_carbon_m3_per_g_c_by_domain = &.{0}, .root_turgor_water_potential_mpa_by_domain_layer = &.{0}, .root_growth_water_potential_mpa_by_domain_layer = &.{0}, .porosity_fraction_by_domain = &.{0.5}, .minimum_primary_radius_m_by_domain = &.{0.001}, .minimum_secondary_radius_m_by_domain = &.{0.001}, .base_primary_radius_m_by_domain = &.{0.002}, .base_secondary_radius_m_by_domain = &.{0.002}, .root_elastic_modulus_megapascal = 1, .secondary_axis_count_by_domain_layer = &.{0}, .minimum_average_secondary_length_m = 0.1, .root_profile_type = 0, .presence_threshold = 0 });
    for (0..3) |species| try std.testing.expectApproxEqAbs(-ledger[species], @as(f64, @floatFromInt(species + 1)) * 1.5, 1e-12);
}
