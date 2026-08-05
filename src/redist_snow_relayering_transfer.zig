const std = @import("std");
const planning = @import("redist_snow_relayering_plan.zig");

pub const SaltMode = enum { static, dynamic };

pub const State = struct {
    layer_count: usize,
    solid_water_equivalent_m3: []f64,
    liquid_water_m3: []f64,
    vapour_water_equivalent_m3: []f64,
    ice_m3: []f64,
    current_volume_m3: []f64,
    snow_density_megagrams_per_m3: []const f64,
    heat_capacity_megajoules_per_k: []f64,
    temperature_k: []f64,
    temperature_c: []f64,
    solute_species_count: usize,
    solute_amount_g: []f64,
    salt_species_count: usize,
    salt_amount_mol: []f64,
    minimum_heat_capacity_megajoules_per_k: f64,
};

fn validateState(state: State) !void {
    if (state.layer_count < 2 or state.solid_water_equivalent_m3.len != state.layer_count or
        state.liquid_water_m3.len != state.layer_count or state.vapour_water_equivalent_m3.len != state.layer_count or
        state.ice_m3.len != state.layer_count or state.current_volume_m3.len != state.layer_count or
        state.snow_density_megagrams_per_m3.len != state.layer_count or state.heat_capacity_megajoules_per_k.len != state.layer_count or
        state.temperature_k.len != state.layer_count or state.temperature_c.len != state.layer_count or
        state.solute_amount_g.len != state.layer_count * state.solute_species_count or
        state.salt_amount_mol.len != state.layer_count * state.salt_species_count)
        return error.SnowRelayeringTransferDimensionMismatch;
    if (!std.math.isFinite(state.minimum_heat_capacity_megajoules_per_k) or state.minimum_heat_capacity_megajoules_per_k < 0)
        return error.InvalidSnowRelayeringTransferState;
    inline for (.{ state.solid_water_equivalent_m3, state.liquid_water_m3, state.vapour_water_equivalent_m3, state.ice_m3, state.current_volume_m3, state.snow_density_megagrams_per_m3, state.heat_capacity_megajoules_per_k, state.temperature_k, state.temperature_c, state.solute_amount_g, state.salt_amount_mol }) |values|
        for (values) |value| if (!std.math.isFinite(value)) return error.InvalidSnowRelayeringTransferState;
    for (state.snow_density_megagrams_per_m3) |density|
        if (density <= 0) return error.InvalidSnowDensity;
}

/// Direct translation of REDIST 7487--7712 for one plan. Invoke immediately
/// after planning each layer to preserve the source's interleaved loop state.
pub fn applyPlan(state: State, plan: planning.TransferPlan, salt_mode: SaltMode) !void {
    try validateState(state);
    if (!std.math.isFinite(plan.fraction) or plan.fraction < 0 or plan.fraction > 1)
        return error.InvalidSnowRelayeringFraction;
    if (plan.fraction == 0) return;
    const source = plan.source_layer;
    const destination = plan.destination_layer;
    if (source >= state.layer_count or destination >= state.layer_count or source == destination)
        return error.InvalidSnowRelayeringTransferLayers;

    const fraction = plan.fraction;
    const retained = 1.0 - fraction;
    const source_energy_megajoules = state.heat_capacity_megajoules_per_k[source] * state.temperature_k[source];
    const destination_energy_megajoules = state.heat_capacity_megajoules_per_k[destination] * state.temperature_k[destination];
    const destination_solid = state.solid_water_equivalent_m3[destination] + fraction * state.solid_water_equivalent_m3[source];
    const destination_liquid = state.liquid_water_m3[destination] + fraction * state.liquid_water_m3[source];
    const destination_vapour = state.vapour_water_equivalent_m3[destination] + fraction * state.vapour_water_equivalent_m3[source];
    const destination_ice = state.ice_m3[destination] + fraction * state.ice_m3[source];
    const destination_volume = destination_solid / state.snow_density_megagrams_per_m3[destination] + destination_liquid + destination_ice;
    const destination_energy = destination_energy_megajoules + fraction * source_energy_megajoules;
    const destination_capacity = 2.095 * destination_solid + 4.19 * (destination_liquid + destination_vapour) + 1.9274 * destination_ice;
    const destination_temperature = if (destination_capacity > state.minimum_heat_capacity_megajoules_per_k)
        destination_energy / destination_capacity
    else
        state.temperature_k[source];

    const source_solid = retained * state.solid_water_equivalent_m3[source];
    const source_liquid = retained * state.liquid_water_m3[source];
    const source_vapour = retained * state.vapour_water_equivalent_m3[source];
    const source_ice = retained * state.ice_m3[source];
    const source_volume = source_solid / state.snow_density_megagrams_per_m3[source] + source_liquid + source_ice;
    const remaining_source_energy = retained * source_energy_megajoules;
    const source_capacity = 2.095 * source_solid + 4.19 * (source_liquid + source_vapour) + 1.9274 * source_ice;
    const source_temperature = if (source_capacity > state.minimum_heat_capacity_megajoules_per_k)
        remaining_source_energy / source_capacity
    else
        destination_temperature;

    inline for (.{ destination_solid, destination_liquid, destination_vapour, destination_ice, destination_volume, destination_capacity, destination_temperature, source_solid, source_liquid, source_vapour, source_ice, source_volume, source_capacity, source_temperature }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSnowRelayeringTransfer;
    for (0..state.solute_species_count) |species| {
        const source_index = source * state.solute_species_count + species;
        const destination_index = destination * state.solute_species_count + species;
        const destination_amount = state.solute_amount_g[destination_index] + fraction * state.solute_amount_g[source_index];
        const source_amount = retained * state.solute_amount_g[source_index];
        if (!std.math.isFinite(destination_amount) or !std.math.isFinite(source_amount))
            return error.NonFiniteSnowRelayeringTransfer;
    }
    if (salt_mode == .dynamic) for (0..state.salt_species_count) |species| {
        const source_index = source * state.salt_species_count + species;
        const destination_index = destination * state.salt_species_count + species;
        const destination_amount = state.salt_amount_mol[destination_index] + fraction * state.salt_amount_mol[source_index];
        const source_amount = retained * state.salt_amount_mol[source_index];
        if (!std.math.isFinite(destination_amount) or !std.math.isFinite(source_amount))
            return error.NonFiniteSnowRelayeringTransfer;
    };

    state.solid_water_equivalent_m3[destination] = destination_solid;
    state.liquid_water_m3[destination] = destination_liquid;
    state.vapour_water_equivalent_m3[destination] = destination_vapour;
    state.ice_m3[destination] = destination_ice;
    state.current_volume_m3[destination] = destination_volume;
    state.heat_capacity_megajoules_per_k[destination] = destination_capacity;
    state.temperature_k[destination] = destination_temperature;
    state.temperature_c[destination] = destination_temperature - 273.15;
    for (0..state.solute_species_count) |species| {
        const source_index = source * state.solute_species_count + species;
        const destination_index = destination * state.solute_species_count + species;
        state.solute_amount_g[destination_index] = state.solute_amount_g[destination_index] + fraction * state.solute_amount_g[source_index];
    }
    if (salt_mode == .dynamic) for (0..state.salt_species_count) |species| {
        const source_index = source * state.salt_species_count + species;
        const destination_index = destination * state.salt_species_count + species;
        state.salt_amount_mol[destination_index] = state.salt_amount_mol[destination_index] + fraction * state.salt_amount_mol[source_index];
    };

    state.solid_water_equivalent_m3[source] = source_solid;
    state.liquid_water_m3[source] = source_liquid;
    state.vapour_water_equivalent_m3[source] = source_vapour;
    state.ice_m3[source] = source_ice;
    state.current_volume_m3[source] = source_volume;
    state.heat_capacity_megajoules_per_k[source] = source_capacity;
    state.temperature_k[source] = source_temperature;
    state.temperature_c[source] = source_temperature - 273.15;
    for (0..state.solute_species_count) |species| {
        const source_index = source * state.solute_species_count + species;
        state.solute_amount_g[source_index] = retained * state.solute_amount_g[source_index];
    }
    if (salt_mode == .dynamic) for (0..state.salt_species_count) |species| {
        const source_index = source * state.salt_species_count + species;
        state.salt_amount_mol[source_index] = retained * state.salt_amount_mol[source_index];
    };
}

const Fixture = struct {
    solid: [2]f64 = .{ 1, 3 },
    liquid: [2]f64 = .{ 2, 4 },
    vapour: [2]f64 = .{ 0.5, 1.5 },
    ice: [2]f64 = .{ 1, 2 },
    volume: [2]f64 = .{ 4, 9 },
    density: [2]f64 = .{ 1, 1 },
    capacity: [2]f64 = .{ 12, 24 },
    temperature_k: [2]f64 = .{ 270, 260 },
    temperature_c: [2]f64 = .{ -3.15, -13.15 },
    solute: [4]f64 = .{ 1, 2, 3, 5 },
    salt: [4]f64 = .{ 7, 11, 13, 17 },
    fn state(self: *Fixture) State {
        return .{ .layer_count = 2, .solid_water_equivalent_m3 = &self.solid, .liquid_water_m3 = &self.liquid, .vapour_water_equivalent_m3 = &self.vapour, .ice_m3 = &self.ice, .current_volume_m3 = &self.volume, .snow_density_megagrams_per_m3 = &self.density, .heat_capacity_megajoules_per_k = &self.capacity, .temperature_k = &self.temperature_k, .temperature_c = &self.temperature_c, .solute_species_count = 2, .solute_amount_g = &self.solute, .salt_species_count = 2, .salt_amount_mol = &self.salt, .minimum_heat_capacity_megajoules_per_k = 0 };
    }
};

fn transferPlan(fraction: f64) planning.TransferPlan {
    return .{ .depth_change_m = 0, .unconstrained_depth_change_m = 0, .source_layer = 1, .destination_layer = 0, .fraction = fraction, .adjustment_mode = .direct };
}

test "REDIST snow transfer conserves physical contents solutes salts and energy" {
    var fixture = Fixture{};
    const water_before = fixture.solid[0] + fixture.solid[1] + fixture.liquid[0] + fixture.liquid[1] + fixture.vapour[0] + fixture.vapour[1] + fixture.ice[0] + fixture.ice[1];
    const energy_before = fixture.capacity[0] * fixture.temperature_k[0] + fixture.capacity[1] * fixture.temperature_k[1];
    try applyPlan(fixture.state(), transferPlan(0.25), .dynamic);
    const water_after = fixture.solid[0] + fixture.solid[1] + fixture.liquid[0] + fixture.liquid[1] + fixture.vapour[0] + fixture.vapour[1] + fixture.ice[0] + fixture.ice[1];
    const energy_after = fixture.capacity[0] * fixture.temperature_k[0] + fixture.capacity[1] * fixture.temperature_k[1];
    try std.testing.expectApproxEqAbs(water_before, water_after, 1e-14);
    try std.testing.expectApproxEqAbs(energy_before, energy_after, 1e-12);
    try std.testing.expectEqual(@as(f64, 4), fixture.solute[0] + fixture.solute[2]);
    try std.testing.expectEqual(@as(f64, 20), fixture.salt[0] + fixture.salt[2]);
}

test "REDIST snow transfer static salt mode leaves salt unmoved" {
    var fixture = Fixture{};
    const before = fixture.salt;
    try applyPlan(fixture.state(), transferPlan(0.5), .static);
    try std.testing.expectEqual(before, fixture.salt);
    try std.testing.expectEqual(@as(f64, 4), fixture.solute[0] + fixture.solute[2]);
}

test "REDIST snow transfer applies fallback temperatures at low capacity" {
    var fixture = Fixture{};
    var state = fixture.state();
    state.minimum_heat_capacity_megajoules_per_k = 1000;
    try applyPlan(state, transferPlan(1.0), .dynamic);
    try std.testing.expectEqual(@as(f64, 260), fixture.temperature_k[0]);
    try std.testing.expectEqual(@as(f64, 260), fixture.temperature_k[1]);
}

test "REDIST snow transfer rejects dimensions fractions and overflow" {
    var fixture = Fixture{};
    var state = fixture.state();
    state.solute_species_count = 3;
    try std.testing.expectError(error.SnowRelayeringTransferDimensionMismatch, applyPlan(state, transferPlan(0.5), .dynamic));
    try std.testing.expectError(error.InvalidSnowRelayeringFraction, applyPlan(fixture.state(), transferPlan(1.1), .dynamic));
    fixture.solute[0] = std.math.floatMax(f64);
    fixture.solute[2] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSnowRelayeringTransfer, applyPlan(fixture.state(), transferPlan(1.0), .dynamic));
}
