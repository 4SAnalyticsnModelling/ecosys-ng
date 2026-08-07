const std = @import("std");
const node_window = @import("../../plant/growth/shoot_growing_node_window.zig");

pub const State = struct {
    leaf_carbon_g_c: []const f64,
    sheath_carbon_g_c: []f64,
    sheath_nitrogen_g_n: []f64,
    sheath_phosphorus_g_p: []f64,
    sheath_protein_g: []f64,
    sheath_height_m: []f64,
};

pub const Inputs = struct {
    newest_growing_node: usize,
    maximum_concurrently_growing_nodes: usize,
    is_main_branch: bool,
    canopy_height_m: f64,
    sowing_depth_m: f64,
    total_carbon_growth_g_c_per_timestep: f64,
    total_nitrogen_growth_g_n_per_timestep: f64,
    total_phosphorus_growth_g_p_per_timestep: f64,
    protein_per_nitrogen_g_per_g_n: f64,
    protein_per_phosphorus_g_per_g_p: f64,
    etiolation_factor: f64,
    base_specific_sheath_length_m_per_g_c: f64,
    minimum_sheath_carbon_g_c: f64,
    plant_population: f64,
    sheath_mass_exponent: f64,
    turgor_extension_fraction: f64,
    sheath_vertical_projection_fraction: f64,
};

/// grosub.f lines 2342--2395: distribute and publish sheath or petiole growth.
/// Runtime nodes replace the fixed 25-node ring; all source assignments and the
/// leaf-presence gate remain in their original order.
pub fn publish(state: State, inputs: Inputs) !void {
    inline for (.{
        inputs.total_carbon_growth_g_c_per_timestep,
        inputs.total_nitrogen_growth_g_n_per_timestep,
        inputs.total_phosphorus_growth_g_p_per_timestep,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSheathNodeGrowth;
    if (inputs.total_carbon_growth_g_c_per_timestep == 0) {
        if (inputs.total_nitrogen_growth_g_n_per_timestep != 0 or
            inputs.total_phosphorus_growth_g_p_per_timestep != 0)
            return error.InconsistentZeroSheathCarbonGrowth;
        return;
    }

    const runtime_node_count = state.sheath_carbon_g_c.len;
    inline for (.{
        state.leaf_carbon_g_c,
        state.sheath_nitrogen_g_n,
        state.sheath_phosphorus_g_p,
        state.sheath_protein_g,
        state.sheath_height_m,
    }) |values| if (values.len != runtime_node_count)
        return error.SheathNodeGrowthDimensionMismatch;
    inline for (.{
        inputs.protein_per_nitrogen_g_per_g_n,
        inputs.protein_per_phosphorus_g_per_g_p,
        inputs.etiolation_factor,
        inputs.base_specific_sheath_length_m_per_g_c,
        inputs.minimum_sheath_carbon_g_c,
        inputs.plant_population,
        inputs.sheath_mass_exponent,
        inputs.turgor_extension_fraction,
        inputs.sheath_vertical_projection_fraction,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSheathNodeGrowthInput;
    if (inputs.protein_per_nitrogen_g_per_g_n < 0 or
        inputs.protein_per_phosphorus_g_per_g_p < 0 or
        inputs.etiolation_factor < 0 or
        inputs.base_specific_sheath_length_m_per_g_c < 0 or
        inputs.minimum_sheath_carbon_g_c < 0 or
        inputs.plant_population <= 0 or
        inputs.turgor_extension_fraction < 0 or
        inputs.sheath_vertical_projection_fraction < 0 or
        inputs.sheath_vertical_projection_fraction > 1)
        return error.InvalidSheathNodeGrowth;

    const window = try node_window.calculate(.{
        .newest_growing_node = inputs.newest_growing_node,
        .maximum_concurrently_growing_nodes = inputs.maximum_concurrently_growing_nodes,
        .runtime_node_count = runtime_node_count,
        .is_main_branch = inputs.is_main_branch,
        .canopy_height_m = inputs.canopy_height_m,
        .sowing_depth_m = inputs.sowing_depth_m,
    });
    const carbon_growth_g_c = window.allocation_fraction *
        inputs.total_carbon_growth_g_c_per_timestep;
    const nitrogen_growth_g_n = window.allocation_fraction *
        inputs.total_nitrogen_growth_g_n_per_timestep;
    const phosphorus_growth_g_p = window.allocation_fraction *
        inputs.total_phosphorus_growth_g_p_per_timestep;

    for (window.first_node..window.last_node + 1) |node|
        try previewNode(state, inputs, node, carbon_growth_g_c, nitrogen_growth_g_n, phosphorus_growth_g_p);

    for (window.first_node..window.last_node + 1) |node| {
        state.sheath_carbon_g_c[node] += carbon_growth_g_c;
        state.sheath_nitrogen_g_n[node] += nitrogen_growth_g_n;
        state.sheath_phosphorus_g_p[node] += phosphorus_growth_g_p;
        state.sheath_protein_g[node] += @min(
            nitrogen_growth_g_n * inputs.protein_per_nitrogen_g_per_g_n,
            phosphorus_growth_g_p * inputs.protein_per_phosphorus_g_per_g_p,
        );
        if (state.leaf_carbon_g_c[node] > 0) {
            const specific_length_m_per_g_c = inputs.etiolation_factor *
                inputs.base_specific_sheath_length_m_per_g_c *
                std.math.pow(
                    f64,
                    @max(inputs.minimum_sheath_carbon_g_c, state.sheath_carbon_g_c[node]) /
                        inputs.plant_population,
                    inputs.sheath_mass_exponent,
                ) * inputs.turgor_extension_fraction;
            const height_growth_m = carbon_growth_g_c / inputs.plant_population *
                specific_length_m_per_g_c;
            state.sheath_height_m[node] += height_growth_m *
                inputs.sheath_vertical_projection_fraction;
        }
    }
}

fn previewNode(
    state: State,
    inputs: Inputs,
    node: usize,
    carbon_growth_g_c: f64,
    nitrogen_growth_g_n: f64,
    phosphorus_growth_g_p: f64,
) !void {
    const carbon_g_c = state.sheath_carbon_g_c[node] + carbon_growth_g_c;
    const nitrogen_g_n = state.sheath_nitrogen_g_n[node] + nitrogen_growth_g_n;
    const phosphorus_g_p = state.sheath_phosphorus_g_p[node] + phosphorus_growth_g_p;
    const protein_g = state.sheath_protein_g[node] + @min(
        nitrogen_growth_g_n * inputs.protein_per_nitrogen_g_per_g_n,
        phosphorus_growth_g_p * inputs.protein_per_phosphorus_g_per_g_p,
    );
    var height_m = state.sheath_height_m[node];
    if (state.leaf_carbon_g_c[node] > 0) {
        const specific_length_m_per_g_c = inputs.etiolation_factor *
            inputs.base_specific_sheath_length_m_per_g_c *
            std.math.pow(
                f64,
                @max(inputs.minimum_sheath_carbon_g_c, carbon_g_c) /
                    inputs.plant_population,
                inputs.sheath_mass_exponent,
            ) * inputs.turgor_extension_fraction;
        const height_growth_m = carbon_growth_g_c / inputs.plant_population *
            specific_length_m_per_g_c;
        height_m += height_growth_m * inputs.sheath_vertical_projection_fraction;
        if (!std.math.isFinite(specific_length_m_per_g_c) or specific_length_m_per_g_c < 0 or
            !std.math.isFinite(height_growth_m) or height_growth_m < 0)
            return error.NonFiniteSheathNodeGrowthResult;
    }
    inline for (.{ carbon_g_c, nitrogen_g_n, phosphorus_g_p, protein_g, height_m }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteSheathNodeGrowthResult;
}

fn sampleInputs() Inputs {
    return .{
        .newest_growing_node = 3,
        .maximum_concurrently_growing_nodes = 3,
        .is_main_branch = false,
        .canopy_height_m = 0.2,
        .sowing_depth_m = 0.03,
        .total_carbon_growth_g_c_per_timestep = 1.5,
        .total_nitrogen_growth_g_n_per_timestep = 0.15,
        .total_phosphorus_growth_g_p_per_timestep = 0.03,
        .protein_per_nitrogen_g_per_g_n = 2,
        .protein_per_phosphorus_g_per_g_p = 10,
        .etiolation_factor = 1.5,
        .base_specific_sheath_length_m_per_g_c = 0.4,
        .minimum_sheath_carbon_g_c = 0.1,
        .plant_population = 2,
        .sheath_mass_exponent = 0.5,
        .turgor_extension_fraction = 0.8,
        .sheath_vertical_projection_fraction = 0.5,
    };
}

test "GROSUB sheath publication preserves window protein and height equations" {
    var leaf = [_]f64{ 1, 1, 0, 2 };
    var carbon = [_]f64{ 1, 1, 1, 1 };
    var nitrogen = [_]f64{ 0, 0, 0, 0 };
    var phosphorus = [_]f64{ 0, 0, 0, 0 };
    var protein = [_]f64{ 0, 0, 0, 0 };
    var height = [_]f64{ 0, 0, 0, 0 };
    const state: State = .{
        .leaf_carbon_g_c = &leaf,
        .sheath_carbon_g_c = &carbon,
        .sheath_nitrogen_g_n = &nitrogen,
        .sheath_phosphorus_g_p = &phosphorus,
        .sheath_protein_g = &protein,
        .sheath_height_m = &height,
    };
    const inputs = sampleInputs();
    try publish(state, inputs);
    try std.testing.expectEqual(@as(f64, 1), carbon[0]);
    try std.testing.expectEqual(@as(f64, 1.5), carbon[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), protein[1], 1.0e-15);
    try std.testing.expectEqual(@as(f64, 0), height[2]);
    const specific_length = 1.5 * 0.4 * std.math.pow(f64, 1.5 / 2.0, 0.5) * 0.8;
    try std.testing.expectApproxEqAbs(0.5 / 2.0 * specific_length * 0.5, height[1], 1.0e-15);
}

test "zero sheath carbon growth returns before node-window evaluation" {
    var none: [0]f64 = .{};
    const state: State = .{
        .leaf_carbon_g_c = &none,
        .sheath_carbon_g_c = &none,
        .sheath_nitrogen_g_n = &none,
        .sheath_phosphorus_g_p = &none,
        .sheath_protein_g = &none,
        .sheath_height_m = &none,
    };
    var inputs = sampleInputs();
    inputs.total_carbon_growth_g_c_per_timestep = 0;
    inputs.total_nitrogen_growth_g_n_per_timestep = 0;
    inputs.total_phosphorus_growth_g_p_per_timestep = 0;
    try publish(state, inputs);
}

test "late invalid node leaves all prior nodes unchanged" {
    var leaf = [_]f64{ 1, 1, 1 };
    var carbon = [_]f64{ 1, 1, std.math.floatMax(f64) };
    var nitrogen = [_]f64{ 0, 0, 0 };
    var phosphorus = [_]f64{ 0, 0, 0 };
    var protein = [_]f64{ 0, 0, 0 };
    var height = [_]f64{ 0, 0, 0 };
    const state: State = .{
        .leaf_carbon_g_c = &leaf,
        .sheath_carbon_g_c = &carbon,
        .sheath_nitrogen_g_n = &nitrogen,
        .sheath_phosphorus_g_p = &phosphorus,
        .sheath_protein_g = &protein,
        .sheath_height_m = &height,
    };
    var inputs = sampleInputs();
    inputs.newest_growing_node = 2;
    inputs.maximum_concurrently_growing_nodes = 2;
    inputs.total_carbon_growth_g_c_per_timestep = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSheathNodeGrowthResult, publish(state, inputs));
    try std.testing.expectEqual(@as(f64, 1), carbon[1]);
    try std.testing.expectEqual(@as(f64, 0), height[1]);
}

test "dimension and inconsistent zero-growth errors are explicit" {
    var leaf = [_]f64{ 1, 1 };
    var one = [_]f64{0};
    var two = [_]f64{ 0, 0 };
    const state: State = .{
        .leaf_carbon_g_c = &leaf,
        .sheath_carbon_g_c = &two,
        .sheath_nitrogen_g_n = &one,
        .sheath_phosphorus_g_p = &two,
        .sheath_protein_g = &two,
        .sheath_height_m = &two,
    };
    try std.testing.expectError(error.SheathNodeGrowthDimensionMismatch, publish(state, sampleInputs()));
    var inputs = sampleInputs();
    inputs.total_carbon_growth_g_c_per_timestep = 0;
    try std.testing.expectError(error.InconsistentZeroSheathCarbonGrowth, publish(state, inputs));
}
