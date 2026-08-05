const std = @import("std");

pub const SaltMode = enum { static_equilibrium, dynamic_transport };

pub const LayerState = struct {
    carbon_transfer_g: []f64, // CSNT
    nitrogen_transfer_g: []f64, // ZSNT
    phosphorus_transfer_g: []f64, // PSNT
    dissolved_organic_carbon_g: []f64, // XOQCS
    dissolved_organic_nitrogen_g: []f64, // XOQNS
    dissolved_organic_phosphorus_g: []f64, // XOQPS
    dissolved_acetate_g: []f64, // XOQAS
    previous_combustion_heat_megajoules_timestep: f64, // HCBFX
    current_combustion_heat_megajoules_timestep: f64, // HCBFL
    root_combustion_g_c_timestep: f64, // RCGSK
    oxygen_consumption_g_o_timestep: f64, // ROGOX
    carbon_dioxide_g_c_timestep: f64, // RCGOX
    carbon_monoxide_g_c_timestep: f64, // RCHOX
    methane_g_c_timestep: f64, // RC4OX
    plant_water_uptake_m3_timestep: f64, // TUPWTR
    plant_heat_uptake_j_timestep: f64, // TUPHT
    carbon_dioxide_diffusion_g_timestep: f64, // XCODFG
    methane_diffusion_g_timestep: f64, // XCHDFG
    oxygen_diffusion_g_timestep: f64, // XOXDFG
    nitrogen_diffusion_g_timestep: f64, // XNGDFG
    nitrous_oxide_diffusion_g_timestep: f64, // XN2DFG
    ammonia_diffusion_g_timestep: f64, // XN3DFG
    nitrogen_balance_diffusion_g_timestep: f64, // XNBDFG
    hydrogen_diffusion_g_timestep: f64, // XHGDFG
    salt_transfer_g: []f64, // ALSNT through CLSNT
    mineral_organic_carbon_flux_g: []f64, // TDFOMC
    mineral_organic_nitrogen_flux_g: []f64, // TDFOMN
    mineral_organic_phosphorus_flux_g: []f64, // TDFOMP
    respiration_iteration_workspace: []f64, // ROXSK
};

pub const Dimensions = struct {
    transfer_pool_count: usize,
    dissolved_organic_pool_count: usize,
    mineral_organic_pool_count: usize,
    salt_species_count: usize,
    maximum_solver_iterations: usize,
    top_mineral_layer_index: usize,
};

pub const GridCombustionHeat = struct {
    previous_megajoules_timestep: f64, // HCBFH
    current_megajoules_timestep: f64, // HCBFG
};

pub const ResetError = error{
    DimensionMismatch,
    NonFiniteCombustionHeat,
};

/// Translates HOUR1 lines 3157-3212. All formerly fixed loop extents are
/// supplied at runtime, including the nonlinear workspace (`NPH`).
pub fn reset(
    salt_mode: SaltMode,
    dimensions: Dimensions,
    grid_combustion_heat: *GridCombustionHeat,
    layers: []LayerState,
) ResetError!void {
    if (!std.math.isFinite(grid_combustion_heat.current_megajoules_timestep)) {
        return error.NonFiniteCombustionHeat;
    }
    for (layers, 0..) |layer, layer_index| {
        if (layer.carbon_transfer_g.len != dimensions.transfer_pool_count or
            layer.nitrogen_transfer_g.len != dimensions.transfer_pool_count or
            layer.phosphorus_transfer_g.len != dimensions.transfer_pool_count or
            layer.dissolved_organic_carbon_g.len != dimensions.dissolved_organic_pool_count or
            layer.dissolved_organic_nitrogen_g.len != dimensions.dissolved_organic_pool_count or
            layer.dissolved_organic_phosphorus_g.len != dimensions.dissolved_organic_pool_count or
            layer.dissolved_acetate_g.len != dimensions.dissolved_organic_pool_count or
            layer.respiration_iteration_workspace.len != dimensions.maximum_solver_iterations or
            (salt_mode == .dynamic_transport and
                layer.salt_transfer_g.len != dimensions.salt_species_count) or
            (layer_index >= dimensions.top_mineral_layer_index and
                (layer.mineral_organic_carbon_flux_g.len != dimensions.mineral_organic_pool_count or
                    layer.mineral_organic_nitrogen_flux_g.len != dimensions.mineral_organic_pool_count or
                    layer.mineral_organic_phosphorus_flux_g.len != dimensions.mineral_organic_pool_count)))
        {
            return error.DimensionMismatch;
        }
        if (!std.math.isFinite(layer.current_combustion_heat_megajoules_timestep)) {
            return error.NonFiniteCombustionHeat;
        }
    }

    grid_combustion_heat.previous_megajoules_timestep = grid_combustion_heat.current_megajoules_timestep;
    grid_combustion_heat.current_megajoules_timestep = 0.0;
    for (layers, 0..) |*layer, layer_index| {
        @memset(layer.carbon_transfer_g, 0.0);
        @memset(layer.nitrogen_transfer_g, 0.0);
        @memset(layer.phosphorus_transfer_g, 0.0);
        @memset(layer.dissolved_organic_carbon_g, 0.0);
        @memset(layer.dissolved_organic_nitrogen_g, 0.0);
        @memset(layer.dissolved_organic_phosphorus_g, 0.0);
        @memset(layer.dissolved_acetate_g, 0.0);
        layer.previous_combustion_heat_megajoules_timestep = layer.current_combustion_heat_megajoules_timestep;
        layer.current_combustion_heat_megajoules_timestep = 0.0;
        inline for (std.meta.fields(LayerState)[9..24]) |field| {
            @field(layer, field.name) = 0.0;
        }
        if (salt_mode == .dynamic_transport) @memset(layer.salt_transfer_g, 0.0);
        if (layer_index >= dimensions.top_mineral_layer_index) {
            @memset(layer.mineral_organic_carbon_flux_g, 0.0);
            @memset(layer.mineral_organic_nitrogen_flux_g, 0.0);
            @memset(layer.mineral_organic_phosphorus_flux_g, 0.0);
        }
        @memset(layer.respiration_iteration_workspace, 0.0);
    }
}

test "runtime transfer pools layers and solver iterations are reset" {
    const allocator = std.testing.allocator;
    const dimensions = Dimensions{
        .transfer_pool_count = 6,
        .dissolved_organic_pool_count = 4,
        .mineral_organic_pool_count = 3,
        .salt_species_count = 8,
        .maximum_solver_iterations = 11,
        .top_mineral_layer_index = 1,
    };
    const layers = try allocator.alloc(LayerState, 3);
    defer allocator.free(layers);
    for (layers) |*layer| {
        inline for (std.meta.fields(LayerState)) |field| {
            if (field.type == f64) @field(layer, field.name) = 5.0;
        }
        layer.carbon_transfer_g = try allocator.alloc(f64, 6);
        layer.nitrogen_transfer_g = try allocator.alloc(f64, 6);
        layer.phosphorus_transfer_g = try allocator.alloc(f64, 6);
        layer.dissolved_organic_carbon_g = try allocator.alloc(f64, 4);
        layer.dissolved_organic_nitrogen_g = try allocator.alloc(f64, 4);
        layer.dissolved_organic_phosphorus_g = try allocator.alloc(f64, 4);
        layer.dissolved_acetate_g = try allocator.alloc(f64, 4);
        layer.salt_transfer_g = try allocator.alloc(f64, 8);
        layer.mineral_organic_carbon_flux_g = try allocator.alloc(f64, 3);
        layer.mineral_organic_nitrogen_flux_g = try allocator.alloc(f64, 3);
        layer.mineral_organic_phosphorus_flux_g = try allocator.alloc(f64, 3);
        layer.respiration_iteration_workspace = try allocator.alloc(f64, 11);
        inline for (std.meta.fields(LayerState)) |field| {
            if (field.type == []f64) @memset(@field(layer, field.name), 5.0);
        }
    }
    defer for (layers) |layer| inline for (std.meta.fields(LayerState)) |field| {
        if (field.type == []f64) allocator.free(@field(layer, field.name));
    };
    var heat = GridCombustionHeat{ .previous_megajoules_timestep = 1.0, .current_megajoules_timestep = 7.0 };

    try reset(.dynamic_transport, dimensions, &heat, layers);

    try std.testing.expectEqual(@as(f64, 7.0), heat.previous_megajoules_timestep);
    try std.testing.expectEqual(@as(f64, 0.0), heat.current_megajoules_timestep);
    for (layers) |layer| {
        try std.testing.expectEqual(@as(f64, 5.0), layer.previous_combustion_heat_megajoules_timestep);
        try std.testing.expectEqual(@as(f64, 0.0), layer.current_combustion_heat_megajoules_timestep);
        for (layer.respiration_iteration_workspace) |value| {
            try std.testing.expectEqual(@as(f64, 0.0), value);
        }
    }
}
