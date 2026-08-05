const std = @import("std");
const growth_stages = @import("plant_growth_stages.zig");

pub const ElementMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const State = struct {
    seasonal_storage: *ElementMass,
    branch_mobile_by_branch: []ElementMass,
    root_structural_carbon_g_c_by_layer: []const f64,
    root_mobile_by_layer: []ElementMass,
};

pub const Inputs = struct {
    current_branch: usize,
    growth_habit: growth_stages.GrowthHabit,
    oxidized_storage_carbon_g_c: f64,
    nutrient_equilibration_fraction_per_h: f64,
    biological_timestep_h: f64,
    shoot_partition_fraction: f64,
    root_partition_fraction: f64,
    carbon_presence_threshold_g_c: f64,
    planting_layer: usize,
};

pub const Result = struct {
    shoot_nitrogen_g_n: f64,
    root_nitrogen_g_n: f64,
    shoot_phosphorus_g_p: f64,
    root_phosphorus_g_p: f64,
    used_root_mobile_carbon_distribution: bool,
};

/// Exact GROSUB lines 4767--4869 seasonal-storage N/P transfer. CH2OH is the
/// carbon oxidized by the retained-DATRP step. Perennial gradients and annual
/// carbon-coupled transfers retain source operation order; root N/P uses
/// runtime CPOOLR fractions or the planting-layer fallback. Units are g C,
/// g N, g P, and h.
pub fn apply(state: State, inputs: Inputs) !Result {
    try validate(state, inputs);
    const storage = state.seasonal_storage.*;
    const branch = state.branch_mobile_by_branch[inputs.current_branch];
    const storage_has_carbon =
        storage.carbon_g_c > inputs.carbon_presence_threshold_g_c;

    var shoot_nitrogen_g_n: f64 = 0;
    var shoot_phosphorus_g_p: f64 = 0;
    if (storage_has_carbon) {
        if (inputs.growth_habit == .perennial) {
            const combined_carbon_g_c = @max(0.0, storage.carbon_g_c + branch.carbon_g_c);
            const nitrogen_gradient_g_n =
                (storage.nitrogen_g_n * branch.carbon_g_c -
                    branch.nitrogen_g_n * storage.carbon_g_c) /
                combined_carbon_g_c;
            const phosphorus_gradient_g_p =
                (storage.phosphorus_g_p * branch.carbon_g_c -
                    branch.phosphorus_g_p * storage.carbon_g_c) /
                combined_carbon_g_c;
            shoot_nitrogen_g_n = @max(0.0, inputs.nutrient_equilibration_fraction_per_h * nitrogen_gradient_g_n) *
                inputs.biological_timestep_h;
            shoot_phosphorus_g_p = @max(0.0, inputs.nutrient_equilibration_fraction_per_h * phosphorus_gradient_g_p) *
                inputs.biological_timestep_h;
        } else {
            shoot_nitrogen_g_n = @max(0.0, inputs.shoot_partition_fraction * inputs.oxidized_storage_carbon_g_c *
                storage.nitrogen_g_n / storage.carbon_g_c);
            shoot_phosphorus_g_p = @max(0.0, inputs.shoot_partition_fraction * inputs.oxidized_storage_carbon_g_c *
                storage.phosphorus_g_p / storage.carbon_g_c);
        }
    } else {
        shoot_nitrogen_g_n = @max(0.0, inputs.shoot_partition_fraction * storage.nitrogen_g_n);
        shoot_phosphorus_g_p = @max(0.0, inputs.shoot_partition_fraction * storage.phosphorus_g_p);
    }

    var root_structural_carbon_g_c: f64 = 0;
    var root_mobile_total: ElementMass = .{};
    for (state.root_structural_carbon_g_c_by_layer, state.root_mobile_by_layer) |structural, root| {
        root_structural_carbon_g_c += structural;
        root_mobile_total.carbon_g_c += root.carbon_g_c;
        root_mobile_total.nitrogen_g_n += root.nitrogen_g_n;
        root_mobile_total.phosphorus_g_p += root.phosphorus_g_p;
        inline for (.{
            root_structural_carbon_g_c,
            root_mobile_total.carbon_g_c,
            root_mobile_total.nitrogen_g_n,
            root_mobile_total.phosphorus_g_p,
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteSeasonalStorageNutrientResult;
    }

    var root_nitrogen_g_n: f64 = 0;
    var root_phosphorus_g_p: f64 = 0;
    if (storage_has_carbon) {
        if (inputs.growth_habit == .perennial) {
            const combined_carbon_g_c = @max(inputs.carbon_presence_threshold_g_c, storage.carbon_g_c + root_mobile_total.carbon_g_c);
            const nitrogen_gradient_g_n =
                (storage.nitrogen_g_n * root_mobile_total.carbon_g_c -
                    root_mobile_total.nitrogen_g_n * storage.carbon_g_c) /
                combined_carbon_g_c;
            const phosphorus_gradient_g_p =
                (storage.phosphorus_g_p * root_mobile_total.carbon_g_c -
                    root_mobile_total.phosphorus_g_p * storage.carbon_g_c) /
                combined_carbon_g_c;
            root_nitrogen_g_n = @max(0.0, inputs.nutrient_equilibration_fraction_per_h * nitrogen_gradient_g_n) *
                inputs.biological_timestep_h;
            root_phosphorus_g_p = @max(0.0, inputs.nutrient_equilibration_fraction_per_h * phosphorus_gradient_g_p) *
                inputs.biological_timestep_h;
        } else {
            root_nitrogen_g_n = @max(0.0, inputs.root_partition_fraction * inputs.oxidized_storage_carbon_g_c *
                storage.nitrogen_g_n / storage.carbon_g_c);
            root_phosphorus_g_p = @max(0.0, inputs.root_partition_fraction * inputs.oxidized_storage_carbon_g_c *
                storage.phosphorus_g_p / storage.carbon_g_c);
        }
    } else {
        root_nitrogen_g_n = @max(0.0, inputs.root_partition_fraction * storage.nitrogen_g_n);
        root_phosphorus_g_p = @max(0.0, inputs.root_partition_fraction * storage.phosphorus_g_p);
    }
    inline for (.{
        shoot_nitrogen_g_n,
        shoot_phosphorus_g_p,
        root_nitrogen_g_n,
        root_phosphorus_g_p,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.NonFiniteSeasonalStorageNutrientResult;

    const next_storage_nitrogen_g_n =
        storage.nitrogen_g_n - shoot_nitrogen_g_n - root_nitrogen_g_n;
    const next_storage_phosphorus_g_p =
        storage.phosphorus_g_p - shoot_phosphorus_g_p - root_phosphorus_g_p;
    const next_branch_nitrogen_g_n = branch.nitrogen_g_n + shoot_nitrogen_g_n;
    const next_branch_phosphorus_g_p = branch.phosphorus_g_p + shoot_phosphorus_g_p;
    inline for (.{
        next_storage_nitrogen_g_n,
        next_storage_phosphorus_g_p,
        next_branch_nitrogen_g_n,
        next_branch_phosphorus_g_p,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.SeasonalStorageNutrientOverdraw;

    const distributed = root_structural_carbon_g_c > inputs.carbon_presence_threshold_g_c and
        root_mobile_total.carbon_g_c > inputs.carbon_presence_threshold_g_c;
    for (state.root_mobile_by_layer, 0..) |root, layer| {
        const fraction = if (distributed)
            @max(0.0, root.carbon_g_c) / root_mobile_total.carbon_g_c
        else
            @as(f64, @floatFromInt(@intFromBool(layer == inputs.planting_layer)));
        inline for (.{
            root.nitrogen_g_n + fraction * root_nitrogen_g_n,
            root.phosphorus_g_p + fraction * root_phosphorus_g_p,
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteSeasonalStorageNutrientResult;
    }

    state.seasonal_storage.nitrogen_g_n = next_storage_nitrogen_g_n;
    state.seasonal_storage.phosphorus_g_p = next_storage_phosphorus_g_p;
    state.branch_mobile_by_branch[inputs.current_branch].nitrogen_g_n =
        next_branch_nitrogen_g_n;
    state.branch_mobile_by_branch[inputs.current_branch].phosphorus_g_p =
        next_branch_phosphorus_g_p;
    for (state.root_mobile_by_layer, 0..) |*root, layer| {
        const fraction = if (distributed)
            @max(0.0, root.carbon_g_c) / root_mobile_total.carbon_g_c
        else
            @as(f64, @floatFromInt(@intFromBool(layer == inputs.planting_layer)));
        root.nitrogen_g_n += fraction * root_nitrogen_g_n;
        root.phosphorus_g_p += fraction * root_phosphorus_g_p;
    }
    return .{
        .shoot_nitrogen_g_n = shoot_nitrogen_g_n,
        .root_nitrogen_g_n = root_nitrogen_g_n,
        .shoot_phosphorus_g_p = shoot_phosphorus_g_p,
        .root_phosphorus_g_p = root_phosphorus_g_p,
        .used_root_mobile_carbon_distribution = distributed,
    };
}

fn validate(state: State, inputs: Inputs) !void {
    if (state.branch_mobile_by_branch.len == 0 or
        inputs.current_branch >= state.branch_mobile_by_branch.len or
        state.root_structural_carbon_g_c_by_layer.len == 0 or
        state.root_structural_carbon_g_c_by_layer.len != state.root_mobile_by_layer.len or
        inputs.planting_layer >= state.root_mobile_by_layer.len)
        return error.SeasonalStorageNutrientDimensionMismatch;
    inline for (.{
        inputs.oxidized_storage_carbon_g_c,
        inputs.nutrient_equilibration_fraction_per_h,
        inputs.biological_timestep_h,
        inputs.shoot_partition_fraction,
        inputs.root_partition_fraction,
        inputs.carbon_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSeasonalStorageNutrientInput;
    if (inputs.biological_timestep_h == 0 or
        @abs(inputs.shoot_partition_fraction + inputs.root_partition_fraction - 1) > 1e-12)
        return error.InvalidSeasonalStorageNutrientInput;
    inline for (.{ state.seasonal_storage.*, state.branch_mobile_by_branch[inputs.current_branch] }) |mass|
        try validateMass(mass);
    for (state.root_structural_carbon_g_c_by_layer, state.root_mobile_by_layer) |structural, root| {
        if (!std.math.isFinite(structural) or structural < 0)
            return error.InvalidSeasonalStorageNutrientState;
        try validateMass(root);
    }
}

fn validateMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSeasonalStorageNutrientState;
}

const Fixture = struct {
    storage: ElementMass = .{ .carbon_g_c = 98, .nitrogen_g_n = 10, .phosphorus_g_p = 1 },
    branches: [3]ElementMass = .{
        .{},
        .{ .carbon_g_c = 20, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .{},
    },
    structural_g_c: [3]f64 = .{ 1, 2, 3 },
    roots: [3]ElementMass = .{
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 },
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 },
        .{ .carbon_g_c = 1, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03 },
    },

    fn state(self: *Fixture) State {
        return .{
            .seasonal_storage = &self.storage,
            .branch_mobile_by_branch = &self.branches,
            .root_structural_carbon_g_c_by_layer = &self.structural_g_c,
            .root_mobile_by_layer = &self.roots,
        };
    }
};

fn testInputs() Inputs {
    return .{
        .current_branch = 1,
        .growth_habit = .annual,
        .oxidized_storage_carbon_g_c = 2,
        .nutrient_equilibration_fraction_per_h = 0.1,
        .biological_timestep_h = 1,
        .shoot_partition_fraction = 0.25,
        .root_partition_fraction = 0.75,
        .carbon_presence_threshold_g_c = 1e-12,
        .planting_layer = 1,
    };
}

fn totals(fixture: *const Fixture) ElementMass {
    var result = fixture.storage;
    result.carbon_g_c += fixture.branches[1].carbon_g_c;
    result.nitrogen_g_n += fixture.branches[1].nitrogen_g_n;
    result.phosphorus_g_p += fixture.branches[1].phosphorus_g_p;
    for (fixture.roots) |root| {
        result.carbon_g_c += root.carbon_g_c;
        result.nitrogen_g_n += root.nitrogen_g_n;
        result.phosphorus_g_p += root.phosphorus_g_p;
    }
    return result;
}

test "GROSUB annual retained-DATRP nutrient transfer conserves C N P" {
    var fixture: Fixture = .{};
    const before = totals(&fixture);
    const result = try apply(fixture.state(), testInputs());
    const after = totals(&fixture);
    try std.testing.expectEqual(before.carbon_g_c, after.carbon_g_c);
    try std.testing.expectApproxEqAbs(before.nitrogen_g_n, after.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before.phosphorus_g_p, after.phosphorus_g_p, 1e-14);
    try std.testing.expect(result.shoot_nitrogen_g_n > 0);
    try std.testing.expect(result.root_nitrogen_g_n > result.shoot_nitrogen_g_n);
    try std.testing.expect(result.used_root_mobile_carbon_distribution);
}

test "carbon transaction publishes post-CH2OH state consumed by annual nutrients" {
    const carbon = @import("seasonal_storage_carbon_remobilization.zig");
    var times_h = [_]f64{1};
    var workspace: carbon.TimeIncrementWorkspace = .{};
    var storage_carbon_g_c: f64 = 100;
    var branch_carbon_g_c = [_]f64{10};
    const structural_carbon_g_c = [_]f64{ 1, 2 };
    var root_carbon_g_c = [_]f64{ 3, 1 };
    const carbon_result = try carbon.apply(.{
        .time_since_germination_h_by_branch = &times_h,
        .time_increment_workspace = &workspace,
        .seasonal_storage_carbon_g_c = &storage_carbon_g_c,
        .branch_mobile_carbon_g_c = &branch_carbon_g_c,
        .root_structural_carbon_g_c_by_layer = &structural_carbon_g_c,
        .root_mobile_carbon_g_c_by_layer = &root_carbon_g_c,
    }, .{
        .current_branch = 0,
        .main_branch = 0,
        .photoperiod_type = .insensitive,
        .phenology = .winter_deciduous,
        .growth_habit = .annual,
        .critical_photoperiod_h = 12,
        .photoperiod_induction_difference_h = 0,
        .daylength_h = 12,
        .canopy_growth_temperature_response = 0.5,
        .canopy_growth_water_fraction = 0.4,
        .biological_timestep_h = 1,
        .remobilization_duration_h = 100,
        .storage_carbon_oxidation_fraction_per_h = 0.1,
        .shoot_partition_fraction = 0.25,
        .root_partition_fraction = 0.75,
        .carbon_presence_threshold_g_c = 1e-12,
        .planting_layer = 0,
    });
    var storage: ElementMass = .{
        .carbon_g_c = storage_carbon_g_c,
        .nitrogen_g_n = 10,
        .phosphorus_g_p = 1,
    };
    var branches = [_]ElementMass{.{
        .carbon_g_c = branch_carbon_g_c[0],
        .nitrogen_g_n = 1,
        .phosphorus_g_p = 0.1,
    }};
    var roots = [_]ElementMass{
        .{ .carbon_g_c = root_carbon_g_c[0], .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 },
        .{ .carbon_g_c = root_carbon_g_c[1], .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 },
    };
    const nutrient_result = try apply(.{
        .seasonal_storage = &storage,
        .branch_mobile_by_branch = &branches,
        .root_structural_carbon_g_c_by_layer = &structural_carbon_g_c,
        .root_mobile_by_layer = &roots,
    }, .{
        .current_branch = 0,
        .growth_habit = .annual,
        .oxidized_storage_carbon_g_c = carbon_result.oxidized_storage_carbon_g_c,
        .nutrient_equilibration_fraction_per_h = 0.1,
        .biological_timestep_h = 1,
        .shoot_partition_fraction = 0.25,
        .root_partition_fraction = 0.75,
        .carbon_presence_threshold_g_c = 1e-12,
        .planting_layer = 0,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.25) * 2.0 * 10.0 / 98.0,
        nutrient_result.shoot_nitrogen_g_n,
        1e-15,
    );
    try std.testing.expectEqual(storage_carbon_g_c, storage.carbon_g_c);
    try std.testing.expectEqual(branch_carbon_g_c[0], branches[0].carbon_g_c);
    try std.testing.expectEqual(root_carbon_g_c[0], roots[0].carbon_g_c);
    try std.testing.expectEqual(root_carbon_g_c[1], roots[1].carbon_g_c);
}

test "GROSUB perennial shoot and root gradients preserve source signs" {
    var fixture: Fixture = .{};
    var inputs = testInputs();
    inputs.growth_habit = .perennial;
    const result = try apply(fixture.state(), inputs);
    try std.testing.expect(result.shoot_nitrogen_g_n > 0);
    try std.testing.expect(result.root_nitrogen_g_n > 0);
    try std.testing.expect(result.shoot_phosphorus_g_p > 0);
    try std.testing.expect(result.root_phosphorus_g_p > 0);
}

test "depleted storage carbon partitions all remaining N and P" {
    var fixture: Fixture = .{};
    fixture.storage.carbon_g_c = 0;
    const before_n = fixture.storage.nitrogen_g_n;
    const before_p = fixture.storage.phosphorus_g_p;
    const result = try apply(fixture.state(), testInputs());
    try std.testing.expectEqual(@as(f64, 0), fixture.storage.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 0), fixture.storage.phosphorus_g_p);
    try std.testing.expectApproxEqAbs(before_n, result.shoot_nitrogen_g_n + result.root_nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(before_p, result.shoot_phosphorus_g_p + result.root_phosphorus_g_p, 1e-15);
}

test "root fallback publishes nutrients only to planting layer" {
    var fixture: Fixture = .{};
    fixture.structural_g_c = .{ 0, 0, 0 };
    const before = fixture.roots;
    const result = try apply(fixture.state(), testInputs());
    try std.testing.expect(!result.used_root_mobile_carbon_distribution);
    try std.testing.expectEqual(before[0], fixture.roots[0]);
    try std.testing.expect(fixture.roots[1].nitrogen_g_n > before[1].nitrogen_g_n);
    try std.testing.expectEqual(before[2], fixture.roots[2]);
}

test "runtime root layers preserve CPOOLR allocation without source ceiling" {
    const allocator = std.testing.allocator;
    const count = 41;
    const structural = try allocator.alloc(f64, count);
    defer allocator.free(structural);
    const roots = try allocator.alloc(ElementMass, count);
    defer allocator.free(roots);
    for (structural, roots, 0..) |*structural_g_c, *root, layer| {
        structural_g_c.* = @floatFromInt(layer + 1);
        root.* = .{ .carbon_g_c = @floatFromInt(layer + 1), .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
    }
    var fixture: Fixture = .{};
    var state = fixture.state();
    state.root_structural_carbon_g_c_by_layer = structural;
    state.root_mobile_by_layer = roots;
    _ = try apply(state, testInputs());
    try std.testing.expect(roots[count - 1].nitrogen_g_n > roots[0].nitrogen_g_n);
}

test "late invalid root and nutrient overdraw leave transaction atomic" {
    var fixture: Fixture = .{};
    fixture.roots[2].phosphorus_g_p = std.math.nan(f64);
    const before = fixture;
    try std.testing.expectError(
        error.InvalidSeasonalStorageNutrientState,
        apply(fixture.state(), testInputs()),
    );
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&before), std.mem.asBytes(&fixture));

    fixture = .{};
    var inputs = testInputs();
    inputs.nutrient_equilibration_fraction_per_h = 1e6;
    inputs.growth_habit = .perennial;
    try std.testing.expectError(
        error.SeasonalStorageNutrientOverdraw,
        apply(fixture.state(), inputs),
    );
}
