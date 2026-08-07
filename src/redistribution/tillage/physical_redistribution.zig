const std = @import("std");

pub const intensive_family_count = 10;
pub const inventory_family_count = 5;
pub const Family = struct { layer_values: []f64, mixed_total: f64 };
pub const Inputs = struct {
    first_soil_layer: usize,
    last_mixed_layer: usize,
    mixing_depth_m: f64, // DCORPZ
    mixing_fraction: f64, // CORP
    cumulative_layer_bottom_m: []const f64,
    layer_thickness_m: []const f64,
    minimum_layer_thickness_m: f64,
    intensive_families: []const Family, // BKDSI,FC,WP,SCNV,SCNH,GKC4,A,M,N,K
    inventory_families: []const Family, // SAND,SILT,CLAY,XCEC,XAEC
    mixed_water_m3: f64,
    mixed_vapor_m3: f64,
    mixed_ice_m3: f64,
    incorporated_surface_water_m3: f64,
    mixed_thermal_energy_megajoules: f64,
    incorporated_surface_energy_megajoules: f64,
};
pub const State = struct {
    water_m3: []f64,
    vapor_m3: []f64,
    ice_m3: []f64,
    water_snapshot_m3: []f64,
    bound_water_m3: []const f64,
    bound_ice_m3: []const f64,
    mineral_heat_capacity_megajoules_k: []const f64,
    heat_capacity_megajoules_k: []f64,
    temperature_k: []f64,
    temperature_c: []f64,
};

fn finite(v: []const f64) bool {
    for (v) |x| if (!std.math.isFinite(x)) return false;
    return true;
}

/// Direct translation of REDIST 12137--12203 for one runtime soil column.
pub fn redistribute(allocator: std.mem.Allocator, inputs: Inputs, state: State) !void {
    const layers = inputs.cumulative_layer_bottom_m.len;
    if (layers == 0 or inputs.first_soil_layer > inputs.last_mixed_layer or inputs.last_mixed_layer >= layers or inputs.layer_thickness_m.len != layers or inputs.intensive_families.len != intensive_family_count or inputs.inventory_families.len != inventory_family_count) return error.TillagePhysicalRedistributionDimensionMismatch;
    inline for (std.meta.fields(State)) |field| {
        if (@field(state, field.name).len != layers) return error.TillagePhysicalRedistributionDimensionMismatch;
        if (!finite(@field(state, field.name))) return error.InvalidTillagePhysicalRedistributionInput;
    }
    for (inputs.intensive_families) |family| if (family.layer_values.len != layers or !finite(family.layer_values) or !std.math.isFinite(family.mixed_total)) return error.InvalidTillagePhysicalRedistributionInput;
    for (inputs.inventory_families) |family| if (family.layer_values.len != layers or !finite(family.layer_values) or !std.math.isFinite(family.mixed_total)) return error.InvalidTillagePhysicalRedistributionInput;
    inline for (.{ inputs.mixing_depth_m, inputs.mixing_fraction, inputs.minimum_layer_thickness_m, inputs.mixed_water_m3, inputs.mixed_vapor_m3, inputs.mixed_ice_m3, inputs.incorporated_surface_water_m3, inputs.mixed_thermal_energy_megajoules, inputs.incorporated_surface_energy_megajoules }) |x| if (!std.math.isFinite(x)) return error.InvalidTillagePhysicalRedistributionInput;
    if (!finite(inputs.cumulative_layer_bottom_m) or !finite(inputs.layer_thickness_m) or inputs.mixing_depth_m <= 0 or inputs.mixing_fraction < 0 or inputs.mixing_fraction > 1 or inputs.minimum_layer_thickness_m < 0) return error.InvalidTillagePhysicalRedistributionInput;

    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();
    const intensive = try arena.allocator().alloc(f64, intensive_family_count * layers);
    const inventory = try arena.allocator().alloc(f64, inventory_family_count * layers);
    const water = try arena.allocator().dupe(f64, state.water_m3);
    const vapor = try arena.allocator().dupe(f64, state.vapor_m3);
    const ice = try arena.allocator().dupe(f64, state.ice_m3);
    const snapshot = try arena.allocator().dupe(f64, state.water_snapshot_m3);
    const heat_capacity = try arena.allocator().dupe(f64, state.heat_capacity_megajoules_k);
    const temperature_k = try arena.allocator().dupe(f64, state.temperature_k);
    const temperature_c = try arena.allocator().dupe(f64, state.temperature_c);
    for (inputs.intensive_families, 0..) |family, i| @memcpy(intensive[i * layers ..][0..layers], family.layer_values);
    for (inputs.inventory_families, 0..) |family, i| @memcpy(inventory[i * layers ..][0..layers], family.layer_values);
    for (inputs.first_soil_layer..inputs.last_mixed_layer + 1) |layer| {
        const thickness = inputs.layer_thickness_m[layer];
        if (thickness <= inputs.minimum_layer_thickness_m) continue;
        const overlap = @min(thickness, inputs.mixing_depth_m - (inputs.cumulative_layer_bottom_m[layer] - thickness));
        const fi = overlap / inputs.mixing_depth_m;
        const ti = overlap / thickness;
        const tx = 1.0 - ti;
        for (inputs.intensive_families, 0..) |family, i| {
            const old = family.layer_values[layer];
            intensive[i * layers + layer] = ti * (old + inputs.mixing_fraction * (family.mixed_total - old)) + tx * old;
        }
        for (inputs.inventory_families, 0..) |family, i| {
            const old = family.layer_values[layer];
            inventory[i * layers + layer] = ti * old + inputs.mixing_fraction * (fi * family.mixed_total - ti * old) + tx * old;
        }
        const mineral_energy = state.mineral_heat_capacity_megajoules_k[layer] * state.temperature_k[layer];
        const volume_energy = (4.19 * (state.water_m3[layer] + state.bound_water_m3[layer]) + 1.9274 * (state.ice_m3[layer] + state.bound_ice_m3[layer])) * state.temperature_k[layer];
        water[layer] = ti * state.water_m3[layer] + inputs.mixing_fraction * (fi * inputs.mixed_water_m3 - ti * state.water_m3[layer]) + tx * state.water_m3[layer] + fi * inputs.incorporated_surface_water_m3;
        vapor[layer] = ti * state.vapor_m3[layer] + inputs.mixing_fraction * (fi * inputs.mixed_vapor_m3 - ti * state.vapor_m3[layer]) + tx * state.vapor_m3[layer];
        ice[layer] = ti * state.ice_m3[layer] + inputs.mixing_fraction * (fi * inputs.mixed_ice_m3 - ti * state.ice_m3[layer]) + tx * state.ice_m3[layer];
        snapshot[layer] = water[layer];
        const layer_energy = ti * volume_energy + inputs.mixing_fraction * (fi * inputs.mixed_thermal_energy_megajoules - ti * volume_energy) + tx * volume_energy + fi * inputs.incorporated_surface_energy_megajoules;
        heat_capacity[layer] = state.mineral_heat_capacity_megajoules_k[layer] + 4.19 * (water[layer] + vapor[layer] + state.bound_water_m3[layer]) + 1.9274 * (ice[layer] + state.bound_ice_m3[layer]);
        if (heat_capacity[layer] <= 0) return error.InvalidTillageHeatCapacityResult;
        temperature_k[layer] = (mineral_energy + layer_energy) / heat_capacity[layer];
        temperature_c[layer] = temperature_k[layer] - 273.15;
        inline for (.{ water[layer], vapor[layer], ice[layer], heat_capacity[layer], temperature_k[layer], temperature_c[layer] }) |x| if (!std.math.isFinite(x)) return error.NonFiniteTillagePhysicalRedistributionResult;
        for (0..intensive_family_count) |i| if (!std.math.isFinite(intensive[i * layers + layer])) return error.NonFiniteTillagePhysicalRedistributionResult;
        for (0..inventory_family_count) |i| if (!std.math.isFinite(inventory[i * layers + layer])) return error.NonFiniteTillagePhysicalRedistributionResult;
    }
    for (inputs.intensive_families, 0..) |family, i| @memcpy(family.layer_values, intensive[i * layers ..][0..layers]);
    for (inputs.inventory_families, 0..) |family, i| @memcpy(family.layer_values, inventory[i * layers ..][0..layers]);
    @memcpy(state.water_m3, water);
    @memcpy(state.vapor_m3, vapor);
    @memcpy(state.ice_m3, ice);
    @memcpy(state.water_snapshot_m3, snapshot);
    @memcpy(state.heat_capacity_megajoules_k, heat_capacity);
    @memcpy(state.temperature_k, temperature_k);
    @memcpy(state.temperature_c, temperature_c);
}

test "REDIST tillage physical redistribution preserves overlap equations and atomic layers" {
    const bottoms = [_]f64{ 0.1, 0.3 };
    const thickness = [_]f64{ 0.1, 0.2 };
    var intensive_values: [intensive_family_count][2]f64 = @splat(@splat(2));
    var inventory_values: [inventory_family_count][2]f64 = @splat(@splat(2));
    var intensive: [intensive_family_count]Family = undefined;
    var inventory: [inventory_family_count]Family = undefined;
    for (0..intensive_family_count) |i| intensive[i] = .{ .layer_values = &intensive_values[i], .mixed_total = 4 };
    for (0..inventory_family_count) |i| inventory[i] = .{ .layer_values = &inventory_values[i], .mixed_total = 4 };
    var water = [_]f64{ 1, 1 };
    var vapor = [_]f64{ 1, 1 };
    var ice = [_]f64{ 1, 1 };
    var snapshot = [_]f64{ 0, 0 };
    const bound = [_]f64{ 0, 0 };
    const mineral = [_]f64{ 10, 10 };
    var capacity = [_]f64{ 0, 0 };
    var tk = [_]f64{ 300, 300 };
    var tc = [_]f64{ 0, 0 };
    try redistribute(std.testing.allocator, .{ .first_soil_layer = 0, .last_mixed_layer = 1, .mixing_depth_m = 0.2, .mixing_fraction = 0.5, .cumulative_layer_bottom_m = &bottoms, .layer_thickness_m = &thickness, .minimum_layer_thickness_m = 0.001, .intensive_families = &intensive, .inventory_families = &inventory, .mixed_water_m3 = 2, .mixed_vapor_m3 = 2, .mixed_ice_m3 = 2, .incorporated_surface_water_m3 = 0.5, .mixed_thermal_energy_megajoules = 1000, .incorporated_surface_energy_megajoules = 10 }, .{ .water_m3 = &water, .vapor_m3 = &vapor, .ice_m3 = &ice, .water_snapshot_m3 = &snapshot, .bound_water_m3 = &bound, .bound_ice_m3 = &bound, .mineral_heat_capacity_megajoules_k = &mineral, .heat_capacity_megajoules_k = &capacity, .temperature_k = &tk, .temperature_c = &tc });
    try std.testing.expectEqual(@as(f64, 3), intensive_values[0][0]);
    try std.testing.expectApproxEqAbs(@as(f64, 2), inventory_values[0][0], 1e-12);
    try std.testing.expectEqual(water[0], snapshot[0]);
    try std.testing.expect(std.math.isFinite(tk[1]));
}
