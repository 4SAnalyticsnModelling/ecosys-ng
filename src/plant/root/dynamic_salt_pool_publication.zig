const std = @import("std");

pub const Ion = enum(u8) {
    aluminum,
    iron,
    calcium,
    magnesium,
    sodium,
    potassium,
    sulfate,
    chloride,
};

pub const ion_count: usize = @typeInfo(Ion).@"enum".fields.len;

pub const Inputs = struct {
    dynamic_salts_enabled: bool,
    biological_domain_count: usize,
    soil_layer_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    /// RUPZ*, flattened [domain][layer][ion], in Ion declaration order (mol).
    uptake_mol_by_domain_layer_ion: []const f64,
};

/// Exact GROSUB dynamic-salt root publication immediately following mineral
/// uptake: N domain outer, L rooted layer next, then Al, Fe, Ca, Mg, Na, K,
/// SO4, and Cl assignments. Signed uptake is added to root inventory. Charge
/// equivalents use the same source ion valences and are diagnostic only.
pub fn apply(
    root_content_mol_by_domain_layer_ion: []f64,
    charge_equivalent_mol_by_domain_layer: []f64,
    inputs: Inputs,
) !void {
    const root_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.RootSaltPublicationDimensionOverflow;
    const value_count = std.math.mul(usize, root_layer_count, ion_count) catch
        return error.RootSaltPublicationDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.uptake_mol_by_domain_layer_ion.len != value_count or
        root_content_mol_by_domain_layer_ion.len != value_count or
        charge_equivalent_mol_by_domain_layer.len != root_layer_count)
        return error.RootSaltPublicationDimensionMismatch;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or
        inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidRootedLayerRange;
    if (!inputs.dynamic_salts_enabled) return;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or
        inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidRootSaltPublicationInput;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidRootLayerThickness;

    for (0..inputs.biological_domain_count) |domain| {
        for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            _ = try calculateNext(root_content_mol_by_domain_layer_ion, inputs, domain, layer);
        }
    }
    for (0..inputs.biological_domain_count) |domain| {
        for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const result = try calculateNext(root_content_mol_by_domain_layer_ion, inputs, domain, layer);
            const root_layer = domain * inputs.soil_layer_count + layer;
            const base = root_layer * ion_count;
            inline for (0..ion_count) |ion| root_content_mol_by_domain_layer_ion[base + ion] = result.content[ion];
            charge_equivalent_mol_by_domain_layer[root_layer] = result.charge_equivalent_mol;
        }
    }
}

const Next = struct {
    content: [ion_count]f64,
    charge_equivalent_mol: f64,
};

fn calculateNext(current: []const f64, inputs: Inputs, domain: usize, layer: usize) !Next {
    const root_layer = domain * inputs.soil_layer_count + layer;
    const base = root_layer * ion_count;
    var next: [ion_count]f64 = undefined;
    inline for (0..ion_count) |ion| {
        const content = current[base + ion];
        const uptake = inputs.uptake_mol_by_domain_layer_ion[base + ion];
        if (!std.math.isFinite(content) or content < 0 or !std.math.isFinite(uptake))
            return error.InvalidRootSaltPublicationState;
        next[ion] = content + uptake;
        if (!std.math.isFinite(next[ion]) or next[ion] < 0)
            return error.RootSaltPublicationWouldOverdraw;
    }
    const uptake = inputs.uptake_mol_by_domain_layer_ion[base .. base + ion_count];
    const charge = 3.0 * (uptake[@intFromEnum(Ion.aluminum)] + uptake[@intFromEnum(Ion.iron)]) +
        2.0 * (uptake[@intFromEnum(Ion.calcium)] + uptake[@intFromEnum(Ion.magnesium)]) +
        uptake[@intFromEnum(Ion.sodium)] + uptake[@intFromEnum(Ion.potassium)] -
        2.0 * uptake[@intFromEnum(Ion.sulfate)] - uptake[@intFromEnum(Ion.chloride)];
    if (!std.math.isFinite(charge)) return error.NonFiniteRootSaltChargeEquivalent;
    return .{ .content = next, .charge_equivalent_mol = charge };
}

test "GROSUB dynamic salt publication preserves exact ion assignment order" {
    var root = [_]f64{10} ** ion_count;
    var charge = [_]f64{0};
    const uptake = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8 };
    try apply(&root, &charge, .{
        .dynamic_salts_enabled = true,
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.1},
        .minimum_active_layer_thickness_m = 0.001,
        .uptake_mol_by_domain_layer_ion = &uptake,
    });
    try std.testing.expectEqualSlices(f64, &.{ 11, 12, 13, 14, 15, 16, 17, 18 }, &root);
    const expected_charge = 3.0 * (1 + 2) + 2.0 * (3 + 4) + 5 + 6 - 2.0 * 7 - 8;
    try std.testing.expectEqual(expected_charge, charge[0]);
}

test "GROSUB signed salt uptake closes root inventory and charge accounting" {
    var root = [_]f64{10} ** ion_count;
    const before = root;
    var charge = [_]f64{0};
    const uptake = [_]f64{ 0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8 };
    try apply(&root, &charge, .{ .dynamic_salts_enabled = true, .biological_domain_count = 1, .soil_layer_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .layer_thickness_m = &.{1}, .minimum_active_layer_thickness_m = 0, .uptake_mol_by_domain_layer_ion = &uptake });
    inline for (0..ion_count) |ion|
        try std.testing.expectApproxEqAbs(uptake[ion], root[ion] - before[ion], 1e-15);
    try std.testing.expectApproxEqAbs(3 * (uptake[0] + uptake[1]) + 2 * (uptake[2] + uptake[3]) + uptake[4] + uptake[5] - 2 * uptake[6] - uptake[7], charge[0], 1e-15);
}

test "GROSUB disabled dynamic salts do not read or modify root ion state" {
    var root = [_]f64{std.math.nan(f64)} ** ion_count;
    var charge = [_]f64{std.math.nan(f64)};
    const uptake = [_]f64{std.math.nan(f64)} ** ion_count;
    try apply(&root, &charge, .{ .dynamic_salts_enabled = false, .biological_domain_count = 1, .soil_layer_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .layer_thickness_m = &.{std.math.nan(f64)}, .minimum_active_layer_thickness_m = std.math.nan(f64), .uptake_mol_by_domain_layer_ion = &uptake });
    try std.testing.expect(std.math.isNan(root[0]));
    try std.testing.expect(std.math.isNan(charge[0]));
}

test "GROSUB runtime root salt sweep is atomic on late ion overdraw" {
    const domains = 5;
    const layers = 10;
    var root = [_]f64{10} ** (domains * layers * ion_count);
    var charge = [_]f64{0} ** (domains * layers);
    var uptake = [_]f64{0.1} ** (domains * layers * ion_count);
    const thickness = [_]f64{0.1} ** layers;
    const inputs = Inputs{ .dynamic_salts_enabled = true, .biological_domain_count = domains, .soil_layer_count = layers, .planting_layer_index = 2, .deepest_rooted_layer_index = 9, .layer_thickness_m = &thickness, .minimum_active_layer_thickness_m = 0.001, .uptake_mol_by_domain_layer_ion = &uptake };
    try apply(&root, &charge, inputs);
    const before_root = root;
    const before_charge = charge;
    uptake[((domains - 1) * layers + 9) * ion_count + @intFromEnum(Ion.chloride)] = -100;
    try std.testing.expectError(error.RootSaltPublicationWouldOverdraw, apply(&root, &charge, inputs));
    try std.testing.expectEqualDeep(before_root, root);
    try std.testing.expectEqualDeep(before_charge, charge);
}
