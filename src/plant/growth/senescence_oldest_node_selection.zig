const std = @import("std");

pub const Inputs = struct {
    lowest_node_remobilization_enabled: bool,
    newest_growing_node: usize,
    maximum_previous_nodes_retained: usize,
    runtime_node_count: usize,
    canopy_growth_temperature_response: f64,
    leaf_appearance_rate_per_h: f64,
    timestep_h: f64,
};

pub const Selection = struct {
    oldest_node: usize,
    newest_node: usize,
    requested_remobilization_fraction: f64,
};

/// grosub.f lines 2542--2549. Logical runtime node ordinals replace K/KX's
/// modulo-25 storage slots. `maximum_previous_nodes_retained=24` reproduces
/// KVSTGX=MAX(0,KVSTG-24) while allowing runtime-configured history depth.
pub fn select(inputs: Inputs) !?Selection {
    if (!inputs.lowest_node_remobilization_enabled) return null;
    if (inputs.runtime_node_count == 0 or
        inputs.newest_growing_node >= inputs.runtime_node_count)
        return error.SenescenceNodeIndexOutOfBounds;

    const oldest_node = inputs.newest_growing_node -|
        inputs.maximum_previous_nodes_retained;
    // The source uses KVSTGX.GT.0, so logical node zero never enters this path.
    if (oldest_node == 0) return null;

    inline for (.{
        inputs.canopy_growth_temperature_response,
        inputs.leaf_appearance_rate_per_h,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSenescenceNodeSelectionInput;
    if (inputs.canopy_growth_temperature_response < 0 or
        inputs.leaf_appearance_rate_per_h < 0 or
        inputs.timestep_h <= 0)
        return error.InvalidSenescenceNodeSelectionInput;

    const requested_fraction = inputs.canopy_growth_temperature_response *
        inputs.leaf_appearance_rate_per_h * inputs.timestep_h;
    if (!std.math.isFinite(requested_fraction) or requested_fraction < 0)
        return error.InvalidSenescenceNodeSelectionResult;
    return .{
        .oldest_node = oldest_node,
        .newest_node = inputs.newest_growing_node,
        .requested_remobilization_fraction = requested_fraction,
    };
}

test "GROSUB source history selects logical ordinal without modulo" {
    const selection = (try select(.{
        .lowest_node_remobilization_enabled = true,
        .newest_growing_node = 50,
        .maximum_previous_nodes_retained = 24,
        .runtime_node_count = 64,
        .canopy_growth_temperature_response = 0.8,
        .leaf_appearance_rate_per_h = 0.04,
        .timestep_h = 0.5,
    })).?;
    try std.testing.expectEqual(@as(usize, 26), selection.oldest_node);
    try std.testing.expectEqual(@as(usize, 50), selection.newest_node);
    try std.testing.expectApproxEqAbs(@as(f64, 0.016), selection.requested_remobilization_fraction, 1.0e-15);
}

test "turnover begins only when source oldest ordinal is greater than zero" {
    try std.testing.expectEqual(@as(?Selection, null), try select(.{
        .lowest_node_remobilization_enabled = true,
        .newest_growing_node = 24,
        .maximum_previous_nodes_retained = 24,
        .runtime_node_count = 32,
        .canopy_growth_temperature_response = 1,
        .leaf_appearance_rate_per_h = 0.1,
        .timestep_h = 1,
    }));
    const first = (try select(.{
        .lowest_node_remobilization_enabled = true,
        .newest_growing_node = 25,
        .maximum_previous_nodes_retained = 24,
        .runtime_node_count = 32,
        .canopy_growth_temperature_response = 1,
        .leaf_appearance_rate_per_h = 0.1,
        .timestep_h = 1,
    })).?;
    try std.testing.expectEqual(@as(usize, 1), first.oldest_node);
}

test "runtime history depth is user controlled" {
    const selection = (try select(.{
        .lowest_node_remobilization_enabled = true,
        .newest_growing_node = 30,
        .maximum_previous_nodes_retained = 10,
        .runtime_node_count = 40,
        .canopy_growth_temperature_response = 1,
        .leaf_appearance_rate_per_h = 0.05,
        .timestep_h = 1,
    })).?;
    try std.testing.expectEqual(@as(usize, 20), selection.oldest_node);
}

test "disabled turnover exits before inactive inputs are evaluated" {
    try std.testing.expectEqual(@as(?Selection, null), try select(.{
        .lowest_node_remobilization_enabled = false,
        .newest_growing_node = 99,
        .maximum_previous_nodes_retained = 24,
        .runtime_node_count = 0,
        .canopy_growth_temperature_response = std.math.nan(f64),
        .leaf_appearance_rate_per_h = -1,
        .timestep_h = 0,
    }));
}

test "active invalid topology and flux controls fail explicitly" {
    try std.testing.expectError(error.SenescenceNodeIndexOutOfBounds, select(.{
        .lowest_node_remobilization_enabled = true,
        .newest_growing_node = 32,
        .maximum_previous_nodes_retained = 24,
        .runtime_node_count = 32,
        .canopy_growth_temperature_response = 1,
        .leaf_appearance_rate_per_h = 0.1,
        .timestep_h = 1,
    }));
    try std.testing.expectError(error.NonFiniteSenescenceNodeSelectionInput, select(.{
        .lowest_node_remobilization_enabled = true,
        .newest_growing_node = 25,
        .maximum_previous_nodes_retained = 24,
        .runtime_node_count = 32,
        .canopy_growth_temperature_response = std.math.inf(f64),
        .leaf_appearance_rate_per_h = 0.1,
        .timestep_h = 1,
    }));
}
