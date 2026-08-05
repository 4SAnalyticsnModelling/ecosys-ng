const std = @import("std");

pub const State = struct {
    total_phytomass_removal_g_c_h_by_species: []f64,
    leaf_structural_removal_g_c_h_by_species: []f64,
    petiole_structural_removal_g_c_h_by_species: []f64,
    husk_removal_g_c_h_by_species: []f64,
    ear_removal_g_c_h_by_species: []f64,
    grain_removal_g_c_h_by_species: []f64,
    nonstructural_carbon_removal_g_c_h_by_species: []f64,
    bacterial_carbon_removal_g_c_h_by_species: []f64,
    stalk_removal_g_c_h_by_species: []f64,
    reserve_removal_g_c_h_by_species: []f64,
    unmet_removal_demand_g_c_h_by_species: []f64,
};

pub const Inputs = struct {
    plant_species_count: usize,
    harvest_material_pool_count: usize,
    harvest_type_by_species: []const i32,
    harvest_height_or_grazer_biomass_g_lm_m2_by_species: []const f64,
    thinning_or_specific_herbivory_rate_d_inv_by_species: []const f64,
    harvest_efficiency_by_species_destination_material_pool: []const f64,
    hour_index: i32,
    local_solar_noon_h: f64,
    planting_layer_area_m2_by_species: []const f64,
    average_shoot_carbon_g_c_by_species: []const f64,
    shoot_carbon_g_c_by_species: []const f64,
    projected_leaf_area_m2_by_species: []const f64,
    projected_stalk_area_m2_by_species: []const f64,
    canopy_temperature_response_by_species: []const f64,
    nonstructural_carbon_ratio_by_species: []const f64,
    bacterial_carbon_ratio_by_species: []const f64,
    leaf_carbon_g_c_by_species: []const f64,
    petiole_carbon_g_c_by_species: []const f64,
    husk_carbon_g_c_by_species: []const f64,
    ear_carbon_g_c_by_species: []const f64,
    grain_carbon_g_c_by_species: []const f64,
    total_nonfoliar_carbon_g_c_by_species: []const f64,
    stalk_carbon_g_c_by_species: []const f64,
    reserve_carbon_g_c_by_species: []const f64,
    total_stalk_reserve_carbon_g_c_by_species: []const f64,
    minimum_pool_g_c: f64,
};

fn efficiencyIndex(inputs: Inputs, species: usize, destination: usize, pool: usize) usize {
    return (species * 2 + destination) * inputs.harvest_material_pool_count + pool;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, species_count: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != species_count) return error.GrazingDemandDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidGrazingDemandState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.plant_species_count == 0 or inputs.harvest_material_pool_count < 3) return error.GrazingDemandDimensionMismatch;
    const efficiency_count = std.math.mul(usize, inputs.plant_species_count, inputs.harvest_material_pool_count) catch return error.GrazingDemandDimensionOverflow;
    const two_destination_count = std.math.mul(usize, efficiency_count, 2) catch return error.GrazingDemandDimensionOverflow;
    if (inputs.harvest_efficiency_by_species_destination_material_pool.len != two_destination_count) return error.GrazingDemandDimensionMismatch;
    for (inputs.harvest_efficiency_by_species_destination_material_pool) |value| if (!std.math.isFinite(value) or value < 0.0 or value > 1.0) return error.InvalidGrazingDemandInput;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64 and !std.mem.startsWith(u8, field.name, "harvest_efficiency")) {
        const values = @field(inputs, field.name);
        if (values.len != inputs.plant_species_count) return error.GrazingDemandDimensionMismatch;
        const signed = std.mem.startsWith(u8, field.name, "harvest_height_or_grazer");
        for (values) |value| if (!std.math.isFinite(value) or (!signed and value < 0.0)) return error.InvalidGrazingDemandInput;
    };
    if (inputs.harvest_type_by_species.len != inputs.plant_species_count or !std.math.isFinite(inputs.local_solar_noon_h) or !std.math.isFinite(inputs.minimum_pool_g_c) or inputs.minimum_pool_g_c < 0.0) return error.InvalidGrazingDemandInput;
    for (inputs.harvest_type_by_species, 0..) |harvest_type, species| if ((harvest_type == 4 or harvest_type == 6) and inputs.harvest_height_or_grazer_biomass_g_lm_m2_by_species[species] < 0.0) return error.InvalidGrazingDemandInput;
}

fn resetRemoval(workspace: State, species: usize) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @field(workspace, field.name)[species] = 0.0;
}

fn eventEnabled(inputs: Inputs, harvest_type: i32) bool {
    return (harvest_type >= 0 and inputs.hour_index == @as(i32, @intFromFloat(@trunc(inputs.local_solar_noon_h))) and harvest_type != 4 and harvest_type != 6) or harvest_type == 4 or harvest_type == 6;
}

/// Exact GROSUB 8626--8785 removal accumulator initialization and animal or
/// insect grazing-demand allocation. Runtime topology is plant species x
/// harvest material pools. Removal is g C h-1; areas are m2; grazer biomass is
/// g live mass m-2; herbivory input is d-1; the source 0.5/24 conversion remains.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs.plant_species_count);
    try validateState(workspace, inputs.plant_species_count);
    copyState(workspace, state);
    for (0..inputs.plant_species_count) |species| {
        const harvest_type = inputs.harvest_type_by_species[species];
        if (!eventEnabled(inputs, harvest_type)) continue;
        resetRemoval(workspace, species);
        if (harvest_type != 4 and harvest_type != 6) continue;

        const total_removal_g_c_h = if (inputs.average_shoot_carbon_g_c_by_species[species] > inputs.minimum_pool_g_c)
            if (harvest_type == 4)
                inputs.harvest_height_or_grazer_biomass_g_lm_m2_by_species[species] * inputs.thinning_or_specific_herbivory_rate_d_inv_by_species[species] * inputs.planting_layer_area_m2_by_species[species] * 0.5 / 24.0 * inputs.shoot_carbon_g_c_by_species[species] / inputs.average_shoot_carbon_g_c_by_species[species]
            else
                inputs.harvest_height_or_grazer_biomass_g_lm_m2_by_species[species] * inputs.thinning_or_specific_herbivory_rate_d_inv_by_species[species] * (inputs.projected_leaf_area_m2_by_species[species] + inputs.projected_stalk_area_m2_by_species[species]) * 0.5 / 24.0 * inputs.canopy_temperature_response_by_species[species] * inputs.shoot_carbon_g_c_by_species[species] / inputs.average_shoot_carbon_g_c_by_species[species]
        else
            0.0;
        workspace.total_phytomass_removal_g_c_h_by_species[species] = total_removal_g_c_h;
        const nonstructural_fraction = inputs.nonstructural_carbon_ratio_by_species[species] / (1.0 + inputs.nonstructural_carbon_ratio_by_species[species]);
        const bacterial_fraction = inputs.bacterial_carbon_ratio_by_species[species] / (1.0 + inputs.bacterial_carbon_ratio_by_species[species]);

        const requested_leaf_g_c_h = total_removal_g_c_h * inputs.harvest_efficiency_by_species_destination_material_pool[efficiencyIndex(inputs, species, 0, 0)];
        const actual_leaf_g_c_h = @min(inputs.leaf_carbon_g_c_by_species[species], requested_leaf_g_c_h);
        var leaf_structural_g_c_h = actual_leaf_g_c_h * (1.0 - nonstructural_fraction);
        var leaf_nonstructural_g_c_h = actual_leaf_g_c_h * nonstructural_fraction;
        var leaf_bacterial_g_c_h = actual_leaf_g_c_h * bacterial_fraction;
        var unmet_g_c_h = @max(0.0, requested_leaf_g_c_h - actual_leaf_g_c_h);
        const requested_nonfoliar_g_c_h = total_removal_g_c_h * inputs.harvest_efficiency_by_species_destination_material_pool[efficiencyIndex(inputs, species, 0, 1)];

        var petiole_structural_g_c_h: f64 = 0.0;
        var petiole_nonstructural_g_c_h: f64 = 0.0;
        var petiole_bacterial_g_c_h: f64 = 0.0;
        var husk_g_c_h: f64 = 0.0;
        var ear_g_c_h: f64 = 0.0;
        var grain_g_c_h: f64 = 0.0;
        const total_nonfoliar_g_c = inputs.total_nonfoliar_carbon_g_c_by_species[species];
        if (total_nonfoliar_g_c > inputs.minimum_pool_g_c) {
            const requested_petiole = requested_nonfoliar_g_c_h * inputs.petiole_carbon_g_c_by_species[species] / total_nonfoliar_g_c + unmet_g_c_h;
            const actual_petiole = @min(inputs.petiole_carbon_g_c_by_species[species], requested_petiole);
            petiole_structural_g_c_h = actual_petiole * (1.0 - nonstructural_fraction);
            petiole_nonstructural_g_c_h = actual_petiole * nonstructural_fraction;
            petiole_bacterial_g_c_h = actual_petiole * bacterial_fraction;
            unmet_g_c_h = @max(0.0, requested_petiole - actual_petiole);
            const requested_husk = requested_nonfoliar_g_c_h * inputs.husk_carbon_g_c_by_species[species] / total_nonfoliar_g_c + unmet_g_c_h;
            husk_g_c_h = @min(inputs.husk_carbon_g_c_by_species[species], requested_husk);
            unmet_g_c_h = @max(0.0, requested_husk - husk_g_c_h);
            const requested_ear = requested_nonfoliar_g_c_h * inputs.ear_carbon_g_c_by_species[species] / total_nonfoliar_g_c + unmet_g_c_h;
            ear_g_c_h = @min(inputs.ear_carbon_g_c_by_species[species], requested_ear);
            unmet_g_c_h = @max(0.0, requested_ear - ear_g_c_h);
            const requested_grain = requested_nonfoliar_g_c_h * inputs.grain_carbon_g_c_by_species[species] / total_nonfoliar_g_c + unmet_g_c_h;
            grain_g_c_h = @min(inputs.grain_carbon_g_c_by_species[species], requested_grain);
            unmet_g_c_h = @max(0.0, requested_grain - grain_g_c_h);
        } else unmet_g_c_h += requested_nonfoliar_g_c_h;

        var stalk_g_c_h: f64 = 0.0;
        var reserve_g_c_h: f64 = 0.0;
        const requested_woody_g_c_h = total_removal_g_c_h * inputs.harvest_efficiency_by_species_destination_material_pool[efficiencyIndex(inputs, species, 0, 2)];
        const total_stalk_reserve_g_c = inputs.total_stalk_reserve_carbon_g_c_by_species[species];
        const stalk_reserve_sufficient = total_stalk_reserve_g_c > requested_woody_g_c_h + unmet_g_c_h;
        if (stalk_reserve_sufficient) {
            const requested_stalk = requested_woody_g_c_h * inputs.stalk_carbon_g_c_by_species[species] / total_stalk_reserve_g_c + unmet_g_c_h;
            stalk_g_c_h = @min(inputs.stalk_carbon_g_c_by_species[species], requested_stalk);
            unmet_g_c_h = @max(0.0, requested_stalk - stalk_g_c_h);
            const requested_reserve = requested_woody_g_c_h * inputs.reserve_carbon_g_c_by_species[species] / total_stalk_reserve_g_c + unmet_g_c_h;
            reserve_g_c_h = @min(inputs.reserve_carbon_g_c_by_species[species], requested_reserve);
            unmet_g_c_h = @max(0.0, requested_reserve - reserve_g_c_h);
        } else {
            stalk_g_c_h = 0.0;
            reserve_g_c_h = 0.0;
            unmet_g_c_h = 0.0;
        }

        // GROSUB 8758 remains inside the 8742 ELSE branch. That branch sets
        // WHVXXX=0 immediately beforehand, so this legacy redistribution is
        // intentionally unreachable without changing the source order.
        if (!stalk_reserve_sufficient and unmet_g_c_h > 0.0) {
            const additional_leaf = @min(inputs.leaf_carbon_g_c_by_species[species] - leaf_structural_g_c_h - leaf_nonstructural_g_c_h, unmet_g_c_h);
            leaf_structural_g_c_h += additional_leaf * (1.0 - nonstructural_fraction);
            leaf_nonstructural_g_c_h += additional_leaf * nonstructural_fraction;
            leaf_bacterial_g_c_h += additional_leaf * bacterial_fraction;
            unmet_g_c_h = @max(0.0, unmet_g_c_h - additional_leaf);
            if (total_nonfoliar_g_c > inputs.minimum_pool_g_c) {
                const requested_petiole = unmet_g_c_h * inputs.petiole_carbon_g_c_by_species[species] / total_nonfoliar_g_c;
                const actual_petiole = @min(inputs.petiole_carbon_g_c_by_species[species], requested_petiole);
                petiole_structural_g_c_h += actual_petiole * (1.0 - nonstructural_fraction);
                petiole_nonstructural_g_c_h += actual_petiole * nonstructural_fraction;
                petiole_bacterial_g_c_h += actual_petiole * bacterial_fraction;
                unmet_g_c_h = @max(0.0, unmet_g_c_h - actual_petiole);
                const requested_husk = unmet_g_c_h * inputs.husk_carbon_g_c_by_species[species] / total_nonfoliar_g_c;
                const actual_husk = @min(inputs.husk_carbon_g_c_by_species[species], requested_husk);
                husk_g_c_h += actual_husk;
                unmet_g_c_h = @max(0.0, unmet_g_c_h - actual_husk);
                const requested_ear = unmet_g_c_h * inputs.ear_carbon_g_c_by_species[species] / total_nonfoliar_g_c;
                const actual_ear = @min(inputs.ear_carbon_g_c_by_species[species], requested_ear);
                ear_g_c_h += actual_ear;
                unmet_g_c_h = @max(0.0, requested_ear - actual_ear);
                const requested_grain = unmet_g_c_h * inputs.grain_carbon_g_c_by_species[species] / total_nonfoliar_g_c;
                const actual_grain = @min(inputs.grain_carbon_g_c_by_species[species], requested_grain);
                grain_g_c_h += actual_grain;
                unmet_g_c_h = @max(0.0, requested_grain - actual_grain);
            }
        }
        workspace.leaf_structural_removal_g_c_h_by_species[species] = leaf_structural_g_c_h;
        workspace.petiole_structural_removal_g_c_h_by_species[species] = petiole_structural_g_c_h;
        workspace.husk_removal_g_c_h_by_species[species] = husk_g_c_h;
        workspace.ear_removal_g_c_h_by_species[species] = ear_g_c_h;
        workspace.grain_removal_g_c_h_by_species[species] = grain_g_c_h;
        workspace.nonstructural_carbon_removal_g_c_h_by_species[species] = leaf_nonstructural_g_c_h + petiole_nonstructural_g_c_h;
        workspace.bacterial_carbon_removal_g_c_h_by_species[species] = leaf_bacterial_g_c_h + petiole_bacterial_g_c_h;
        workspace.stalk_removal_g_c_h_by_species[species] = stalk_g_c_h;
        workspace.reserve_removal_g_c_h_by_species[species] = reserve_g_c_h;
        workspace.unmet_removal_demand_g_c_h_by_species[species] = unmet_g_c_h;
    }
    try validateState(workspace, inputs.plant_species_count);
    copyState(state, workspace);
}

fn makeState(storage: *[11]f64) State {
    return .{ .total_phytomass_removal_g_c_h_by_species = storage[0..1], .leaf_structural_removal_g_c_h_by_species = storage[1..2], .petiole_structural_removal_g_c_h_by_species = storage[2..3], .husk_removal_g_c_h_by_species = storage[3..4], .ear_removal_g_c_h_by_species = storage[4..5], .grain_removal_g_c_h_by_species = storage[5..6], .nonstructural_carbon_removal_g_c_h_by_species = storage[6..7], .bacterial_carbon_removal_g_c_h_by_species = storage[7..8], .stalk_removal_g_c_h_by_species = storage[8..9], .reserve_removal_g_c_h_by_species = storage[9..10], .unmet_removal_demand_g_c_h_by_species = storage[10..11] };
}

test "GROSUB animal grazing demand partitions leaf carbon exactly" {
    var values = [_]f64{0} ** 11;
    var work = [_]f64{0} ** 11;
    const inputs: Inputs = .{ .plant_species_count = 1, .harvest_material_pool_count = 4, .harvest_type_by_species = &.{4}, .harvest_height_or_grazer_biomass_g_lm_m2_by_species = &.{24}, .thinning_or_specific_herbivory_rate_d_inv_by_species = &.{1}, .harvest_efficiency_by_species_destination_material_pool = &.{ 1, 0, 0, 0, 0, 0, 0, 0 }, .hour_index = 1, .local_solar_noon_h = 12, .planting_layer_area_m2_by_species = &.{1}, .average_shoot_carbon_g_c_by_species = &.{100}, .shoot_carbon_g_c_by_species = &.{100}, .projected_leaf_area_m2_by_species = &.{1}, .projected_stalk_area_m2_by_species = &.{1}, .canopy_temperature_response_by_species = &.{1}, .nonstructural_carbon_ratio_by_species = &.{0.25}, .bacterial_carbon_ratio_by_species = &.{0.1111111111111111}, .leaf_carbon_g_c_by_species = &.{10}, .petiole_carbon_g_c_by_species = &.{0}, .husk_carbon_g_c_by_species = &.{0}, .ear_carbon_g_c_by_species = &.{0}, .grain_carbon_g_c_by_species = &.{0}, .total_nonfoliar_carbon_g_c_by_species = &.{0}, .stalk_carbon_g_c_by_species = &.{1}, .reserve_carbon_g_c_by_species = &.{1}, .total_stalk_reserve_carbon_g_c_by_species = &.{2}, .minimum_pool_g_c = 1e-12 };
    try apply(makeState(&values), makeState(&work), inputs);
    try std.testing.expectApproxEqAbs(0.5, values[0], 1e-14);
    try std.testing.expectApproxEqAbs(0.4, values[1], 1e-14);
    try std.testing.expectApproxEqAbs(0.1, values[6], 1e-14);
    try std.testing.expectApproxEqAbs(0.05, values[7], 1e-14);
    try std.testing.expectApproxEqAbs(0.0, values[10], 1e-14);
}
