const std = @import("std");
const grid_module = @import("grid.zig");
const catalog_module = @import("soil_catalog.zig");
const retention = @import("soil_water_retention.zig");
const conductivity = @import("soil_hydraulic_conductivity.zig");
const profile_derivation = @import("soil_profile_derivation.zig");

pub const RuntimeParameters = struct {
    retention: retention.Parameters,
    mualem_van_genuchten_fit_max_iterations: u16,
    hydraulic_conductivity_class_count: usize,
    pore_interaction_exponent: f64,
    air_entry_fraction_of_vertical_saturated_conductivity: f64,
    mineral_saturated_conductivity_scale_m2_per_h_mpa: f64,
    mineral_reference_water_potential_mpa_magnitude: f64,
    organic_saturated_conductivity_intercept_m2_per_h_mpa: f64,
    organic_saturated_conductivity_scale_m2_per_h_mpa: f64,
    organic_saturated_conductivity_bulk_density_base: f64,
    profile_derivation: profile_derivation.Parameters,
};

/// Values encoded by HOUR1 in the historical source. This constructs a runtime
/// value only for old runscripts that predate the explicit soil_solver record.
pub fn compatibilityParameters() RuntimeParameters {
    return .{
        .retention = retention.compatibilityParameters(),
        .mualem_van_genuchten_fit_max_iterations = 80,
        .hydraulic_conductivity_class_count = 100,
        .pore_interaction_exponent = 1.33,
        .air_entry_fraction_of_vertical_saturated_conductivity = 0.1,
        .mineral_saturated_conductivity_scale_m2_per_h_mpa = 1.54,
        .mineral_reference_water_potential_mpa_magnitude = 0.033,
        .organic_saturated_conductivity_intercept_m2_per_h_mpa = 0.10,
        .organic_saturated_conductivity_scale_m2_per_h_mpa = 75,
        .organic_saturated_conductivity_bulk_density_base = 1.0e-15,
        .profile_derivation = profile_derivation.compatibilityParameters(),
    };
}

/// Heap-owned, cell-layer-resolved HOUR1 properties consumed by the coupled
/// soil solvers. Every scientific coefficient enters through RuntimeParameters
/// or a user soil profile; none of the extents or coefficients are comptime.
pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    hydraulic_conductivity_class_count: usize,
    retention_curve: []retention.ResolvedCurve,
    mualem_van_genuchten_parameters: []retention.MualemVanGenuchtenParameters = @constCast(&.{}),
    matrix_bulk_volume_m3: []f64,
    layer_volume_m3: []f64,
    layer_thickness_m: []f64,
    layer_midpoint_depth_m: []f64,
    layer_bottom_depth_m: []f64,
    bulk_density_megagrams_per_m3: []f64,
    sand_mass_fraction: []f64,
    clay_mass_fraction: []f64,
    /// STARTS/REDIST SAND, SILT, CLAY authoritative extensive inventories.
    sand_mass_Mg: []f64,
    silt_mass_Mg: []f64,
    clay_mass_Mg: []f64,
    total_organic_carbon_g_per_megagram: []f64,
    cation_exchange_capacity_mol_per_Mg: []f64,
    anion_exchange_capacity_mol_per_Mg: []f64,
    cation_exchange_capacity_mol: []f64,
    anion_exchange_capacity_mol: []f64,
    porosity_fraction: []f64,
    matrix_air_entry_water_fraction: []f64,
    saturation_water_potential_mpa: []f64,
    /// Rainfall-impact reduction applied during lookup; the base table remains immutable.
    rainfall_conductivity_multiplier: []f64,
    /// Cell-major, then x/y/z, then wettest-to-driest water class.
    matrix_hydraulic_conductivity_m2_per_h_mpa: []f64,

    pub fn initMapped(
        allocator: std.mem.Allocator,
        grid: *const grid_module.GridState,
        catalog_entries: []const catalog_module.Entry,
        catalog_index_by_cell: []const usize,
        horizontal_cell_width_m: []const f64,
        vertical_cell_width_m: []const f64,
        parameters: RuntimeParameters,
    ) !State {
        try parameters.retention.validate();
        try parameters.profile_derivation.validate();
        if (parameters.mualem_van_genuchten_fit_max_iterations == 0 or parameters.hydraulic_conductivity_class_count == 0 or !std.math.isFinite(parameters.pore_interaction_exponent) or parameters.pore_interaction_exponent <= 0 or !std.math.isFinite(parameters.air_entry_fraction_of_vertical_saturated_conductivity) or parameters.air_entry_fraction_of_vertical_saturated_conductivity < 0 or parameters.air_entry_fraction_of_vertical_saturated_conductivity > 1 or !std.math.isFinite(parameters.mineral_saturated_conductivity_scale_m2_per_h_mpa) or parameters.mineral_saturated_conductivity_scale_m2_per_h_mpa < 0 or !std.math.isFinite(parameters.mineral_reference_water_potential_mpa_magnitude) or parameters.mineral_reference_water_potential_mpa_magnitude <= 0 or !std.math.isFinite(parameters.organic_saturated_conductivity_intercept_m2_per_h_mpa) or parameters.organic_saturated_conductivity_intercept_m2_per_h_mpa < 0 or !std.math.isFinite(parameters.organic_saturated_conductivity_scale_m2_per_h_mpa) or parameters.organic_saturated_conductivity_scale_m2_per_h_mpa < 0 or !std.math.isFinite(parameters.organic_saturated_conductivity_bulk_density_base) or parameters.organic_saturated_conductivity_bulk_density_base <= 0) return error.InvalidSoilSolverRuntimeParameters;
        if (catalog_index_by_cell.len != grid.cell_count or
            horizontal_cell_width_m.len != grid.cell_count or
            vertical_cell_width_m.len != grid.cell_count)
            return error.SoilSolverPropertyDimensionMismatch;
        var result: State = undefined;
        result.allocator = allocator;
        result.layer_count = grid.layer_count;
        result.hydraulic_conductivity_class_count = parameters.hydraulic_conductivity_class_count;
        var allocated: usize = 0;
        errdefer result.freeAllocated(allocated);
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                const length = if (comptime std.mem.eql(u8, field.name, "matrix_hydraulic_conductivity_m2_per_h_mpa")) try std.math.mul(usize, grid.layer_count, try std.math.mul(usize, 3, parameters.hydraulic_conductivity_class_count)) else grid.layer_count;
                @field(result, field.name) = try allocator.alloc(f64, length);
                @memset(@field(result, field.name), 0);
                allocated += 1;
            } else if (field.type == []retention.ResolvedCurve) {
                @field(result, field.name) = try allocator.alloc(retention.ResolvedCurve, grid.layer_count);
                allocated += 1;
            } else if (field.type == []retention.MualemVanGenuchtenParameters) {
                @field(result, field.name) = try allocator.alloc(retention.MualemVanGenuchtenParameters, grid.layer_count);
                allocated += 1;
            }
        }
        const inactive_curve = try retention.resolve(parameters.retention, .{ .porosity_fraction = 0.5, .macropore_fraction = 0, .sand_fraction = 0.5, .clay_fraction = 0.25, .organic_carbon_g_per_megagram = 0, .bulk_density_megagrams_per_m3 = 1, .supplied_field_capacity_fraction = null, .supplied_wilting_point_fraction = null }, -0.01, -1.5);
        @memset(result.retention_curve, inactive_curve);
        const inactive_mualem_van_genuchten = try retention.carselParrishDefault(.loam, 0.5);
        @memset(result.mualem_van_genuchten_parameters, inactive_mualem_van_genuchten);
        @memset(result.matrix_bulk_volume_m3, 1);
        @memset(result.layer_volume_m3, 1);
        @memset(result.layer_thickness_m, 1);
        @memset(result.layer_midpoint_depth_m, 0.5);
        @memset(result.layer_bottom_depth_m, 1);
        @memset(result.bulk_density_megagrams_per_m3, 1);
        @memset(result.sand_mass_fraction, 0.5);
        @memset(result.clay_mass_fraction, 0.25);
        @memset(result.sand_mass_Mg, 0.5);
        @memset(result.silt_mass_Mg, 0.25);
        @memset(result.clay_mass_Mg, 0.25);
        @memset(result.total_organic_carbon_g_per_megagram, 0);
        @memset(result.cation_exchange_capacity_mol_per_Mg, 0);
        @memset(result.anion_exchange_capacity_mol_per_Mg, 0);
        @memset(result.cation_exchange_capacity_mol, 0);
        @memset(result.anion_exchange_capacity_mol, 0);
        @memset(result.porosity_fraction, 0.5);
        @memset(result.matrix_air_entry_water_fraction, inactive_curve.porosity_fraction);
        @memset(result.saturation_water_potential_mpa, inactive_curve.curve.saturation_water_potential_mpa);
        @memset(result.rainfall_conductivity_multiplier, 1);

        for (catalog_index_by_cell, 0..) |catalog_index, cell| {
            if (catalog_index >= catalog_entries.len) return error.SoilCatalogMapOutOfBounds;
            const entry = &catalog_entries[catalog_index];
            if (entry.profile.total_layer_count != grid.active_soil_layer_count[cell]) return error.SoilLayerCountMismatch;
            const area_m2 = horizontal_cell_width_m[cell] * vertical_cell_width_m[cell];
            var layer_top_depth_m: f64 = 0;
            for (0..entry.profile.total_layer_count) |layer| {
                const index = cell * grid.soil_layer_capacity + layer;
                const curve = try retention.resolve(parameters.retention, .{
                    .porosity_fraction = entry.hydrology_per_m2.total_porosity_fraction[layer],
                    .macropore_fraction = entry.profile.macropore_fraction[layer],
                    .sand_fraction = entry.material.sand_mass_fraction[layer],
                    .clay_fraction = entry.material.clay_mass_fraction[layer],
                    .organic_carbon_g_per_megagram = entry.material.total_organic_carbon_g_per_megagram[layer],
                    .bulk_density_megagrams_per_m3 = entry.material.bulk_density_megagrams_per_m3[layer],
                    .supplied_field_capacity_fraction = supplied(entry.profile.field_capacity_m3_m3[layer]),
                    .supplied_wilting_point_fraction = supplied(entry.profile.wilting_point_m3_m3[layer]),
                }, entry.profile.field_capacity_potential_mpa, entry.profile.wilting_point_potential_mpa);
                const texture = try retention.classifyUsdaSoilTexture(
                    entry.material.sand_mass_fraction[layer],
                    entry.material.silt_mass_fraction[layer],
                    entry.material.clay_mass_fraction[layer],
                );
                var mualem_van_genuchten = try retention.carselParrishDefault(
                    texture,
                    curve.porosity_fraction,
                );
                if (entry.profile.vertical_saturated_conductivity_mm_h[layer] >= 0)
                    mualem_van_genuchten.saturated_hydraulic_conductivity_m_per_h =
                        entry.profile.vertical_saturated_conductivity_mm_h[layer] / 1000.0;
                const inflection_pressure_head_m =
                    entry.profile.van_genuchten_inflection_pressure_head_m[layer];
                if (inflection_pressure_head_m < 0) {
                    const fitted = try retention.fitOriginalMualemVanGenuchten(.{
                        .saturated_water_content_m3_per_m3 = curve.porosity_fraction,
                        .field_capacity_water_content_m3_per_m3 = curve.curve.field_capacity_fraction,
                        .field_capacity_pressure_head_m = try mpaToPressureHeadM(
                            entry.profile.field_capacity_potential_mpa,
                        ),
                        .wilting_point_water_content_m3_per_m3 = curve.curve.wilting_point_fraction,
                        .wilting_point_pressure_head_m = try mpaToPressureHeadM(
                            entry.profile.wilting_point_potential_mpa,
                        ),
                        .inflection_pressure_head_m = inflection_pressure_head_m,
                        .saturated_hydraulic_conductivity_m_per_h = mualem_van_genuchten.saturated_hydraulic_conductivity_m_per_h,
                    }, .{
                        .maximum_iterations = parameters.mualem_van_genuchten_fit_max_iterations,
                    });
                    mualem_van_genuchten = fitted.parameters;
                }
                // Conductivity remains an intrinsic m2 MPa-1 h-1 property.
                // Face area is applied exactly once by the flux kernel.
                const estimated_saturated_conductivity = try estimateSaturatedConductivity(parameters, curve, entry.material.total_organic_carbon_g_per_megagram[layer], entry.material.bulk_density_megagrams_per_m3[layer], entry.material.micropore_fraction[layer]);
                const vertical_saturated_conductivity = if (entry.profile.vertical_saturated_conductivity_mm_h[layer] < 0)
                    estimated_saturated_conductivity
                else
                    entry.material.vertical_hydraulic_conductivity_m2_per_mpa_h[layer];
                const lateral_saturated_conductivity = if (entry.profile.lateral_saturated_conductivity_mm_h[layer] < 0)
                    estimated_saturated_conductivity
                else
                    entry.material.lateral_hydraulic_conductivity_m2_per_mpa_h[layer];
                var table = try conductivity.build(allocator, curve, lateral_saturated_conductivity, vertical_saturated_conductivity, .{ .class_count = parameters.hydraulic_conductivity_class_count, .pore_interaction_exponent = parameters.pore_interaction_exponent, .air_entry_fraction_of_vertical_saturated_conductivity = parameters.air_entry_fraction_of_vertical_saturated_conductivity });
                defer table.deinit();
                result.retention_curve[index] = curve;
                result.mualem_van_genuchten_parameters[index] =
                    mualem_van_genuchten;
                result.clay_mass_fraction[index] = entry.material.clay_mass_fraction[layer];
                result.sand_mass_fraction[index] = entry.material.sand_mass_fraction[layer];
                const bulk_soil_mass_Mg = entry.material.bulk_density_megagrams_per_m3[layer] * entry.hydrology_per_m2.total_layer_volume_m3[layer] * area_m2;
                result.sand_mass_Mg[index] = entry.material.sand_mass_fraction[layer] * bulk_soil_mass_Mg;
                result.silt_mass_Mg[index] = entry.material.silt_mass_fraction[layer] * bulk_soil_mass_Mg;
                result.clay_mass_Mg[index] = entry.material.clay_mass_fraction[layer] * bulk_soil_mass_Mg;
                result.total_organic_carbon_g_per_megagram[index] = entry.material.total_organic_carbon_g_per_megagram[layer];
                result.cation_exchange_capacity_mol_per_Mg[index] = entry.material.cation_exchange_capacity_mol_per_megagram[layer];
                result.anion_exchange_capacity_mol_per_Mg[index] = 10.0 * entry.profile.anion_exchange_capacity_cmol_kg[layer];
                result.cation_exchange_capacity_mol[index] = result.cation_exchange_capacity_mol_per_Mg[index] * bulk_soil_mass_Mg;
                result.anion_exchange_capacity_mol[index] = result.anion_exchange_capacity_mol_per_Mg[index] * bulk_soil_mass_Mg;
                result.matrix_bulk_volume_m3[index] = entry.hydrology_per_m2.matrix_volume_m3[layer] * area_m2;
                result.layer_volume_m3[index] = entry.hydrology_per_m2.total_layer_volume_m3[layer] * area_m2;
                result.layer_thickness_m[index] = entry.hydrology_per_m2.layer_thickness_m[layer];
                result.layer_midpoint_depth_m[index] = layer_top_depth_m + 0.5 * entry.hydrology_per_m2.layer_thickness_m[layer];
                layer_top_depth_m += entry.hydrology_per_m2.layer_thickness_m[layer];
                result.layer_bottom_depth_m[index] = layer_top_depth_m;
                result.bulk_density_megagrams_per_m3[index] = entry.material.bulk_density_megagrams_per_m3[layer];
                result.porosity_fraction[index] = curve.porosity_fraction;
                result.matrix_air_entry_water_fraction[index] = table.air_entry_water_fraction;
                result.saturation_water_potential_mpa[index] = curve.curve.saturation_water_potential_mpa;
                const classes_per_cell = 3 * parameters.hydraulic_conductivity_class_count;
                @memcpy(result.matrix_hydraulic_conductivity_m2_per_h_mpa[index * classes_per_cell ..][0..classes_per_cell], table.conductivity_m2_per_h_mpa);
            }
        }
        try result.validateFinite();
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 or field.type == []retention.ResolvedCurve or field.type == []retention.MualemVanGenuchtenParameters) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn validateFinite(self: *const State) !void {
        for (self.retention_curve) |curve| {
            if (!std.math.isFinite(curve.porosity_fraction) or curve.porosity_fraction <= 0) return error.NonFiniteSoilSolverProperty;
        }
        for (self.mualem_van_genuchten_parameters) |parameters|
            try parameters.validate();
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFiniteSoilSolverProperty;
    }

    fn freeAllocated(self: *State, count: usize) void {
        var visited: usize = 0;
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 or field.type == []retention.ResolvedCurve or field.type == []retention.MualemVanGenuchtenParameters) {
            if (visited < count) self.allocator.free(@field(self, field.name));
            visited += 1;
        };
    }
};

fn estimateSaturatedConductivity(parameters: RuntimeParameters, curve: retention.ResolvedCurve, organic_carbon_g_per_megagram: f64, bulk_density_megagrams_per_m3: f64, micropore_fraction: f64) !f64 {
    inline for (.{ parameters.mineral_saturated_conductivity_scale_m2_per_h_mpa, parameters.mineral_reference_water_potential_mpa_magnitude, parameters.organic_saturated_conductivity_intercept_m2_per_h_mpa, parameters.organic_saturated_conductivity_scale_m2_per_h_mpa, parameters.organic_saturated_conductivity_bulk_density_base, organic_carbon_g_per_megagram, bulk_density_megagrams_per_m3, micropore_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSaturatedConductivityInput;
    if (parameters.mineral_saturated_conductivity_scale_m2_per_h_mpa < 0 or parameters.mineral_reference_water_potential_mpa_magnitude <= 0 or parameters.organic_saturated_conductivity_intercept_m2_per_h_mpa < 0 or parameters.organic_saturated_conductivity_scale_m2_per_h_mpa < 0 or parameters.organic_saturated_conductivity_bulk_density_base <= 0 or organic_carbon_g_per_megagram < 0 or bulk_density_megagrams_per_m3 < 0 or micropore_fraction < 0 or micropore_fraction > 1) return error.InvalidSaturatedConductivityInput;
    const value = if (organic_carbon_g_per_megagram < parameters.retention.organic_soil_threshold_g_per_megagram) mineral: {
        const porosity = curve.porosity_fraction;
        const field_capacity = curve.curve.field_capacity_fraction;
        const log_saturation_potential = @log(-curve.curve.saturation_water_potential_mpa);
        const log_field_potential = @log(-curve.curve.field_capacity_water_potential_mpa);
        const water_at_reference = @min(porosity, @exp(
            (log_saturation_potential - @log(parameters.mineral_reference_water_potential_mpa_magnitude)) *
                (@log(porosity) - @log(field_capacity)) /
                (log_field_potential - log_saturation_potential) +
                @log(porosity),
        ));
        break :mineral parameters.mineral_saturated_conductivity_scale_m2_per_h_mpa * std.math.pow(f64, (porosity - water_at_reference) / water_at_reference, 2);
    } else (parameters.organic_saturated_conductivity_intercept_m2_per_h_mpa +
        parameters.organic_saturated_conductivity_scale_m2_per_h_mpa *
            std.math.pow(f64, parameters.organic_saturated_conductivity_bulk_density_base, bulk_density_megagrams_per_m3)) *
        micropore_fraction;
    if (!std.math.isFinite(value) or value < 0) return error.InvalidResolvedSaturatedConductivity;
    return value;
}

fn supplied(value: f64) ?f64 {
    return if (value >= 0) value else null;
}

fn mpaToPressureHeadM(water_potential_mpa: f64) !f64 {
    if (!std.math.isFinite(water_potential_mpa) or water_potential_mpa >= 0)
        return error.InvalidSoilWaterPotential;
    const water_density_kg_per_m3 = 1000.0;
    const gravitational_acceleration_m_per_s2 = 9.80665;
    const pressure_head_m =
        water_potential_mpa * 1_000_000.0 /
        (water_density_kg_per_m3 * gravitational_acceleration_m_per_s2);
    if (!std.math.isFinite(pressure_head_m) or pressure_head_m >= 0)
        return error.InvalidSoilWaterPotential;
    return pressure_head_m;
}

test "mapped solver properties use runtime dimensions classes and profile science" {
    const allocator = std.testing.allocator;
    const source = try @import("test_fixtures.zig").soilProfileSource(allocator, @typeInfo(@import("soil_profile.zig").LayerProperty).@"enum".fields.len);
    defer allocator.free(source);
    var catalog = catalog_module.Catalog.init(allocator);
    defer catalog.deinit();
    _ = try catalog.appendFromSource("soil", source, compatibilityParameters().retention, compatibilityParameters().profile_derivation);
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = catalog.entries.items[0].profile.total_layer_count, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(allocator, cfg);
    defer grid.deinit();
    try @import("model_initialization.zig").initializeCellHydrology(&grid, 0, catalog.entries.items[0].hydrology_per_m2);
    var parameters = compatibilityParameters();
    parameters.hydraulic_conductivity_class_count = 37;
    catalog.entries.items[0].profile.vertical_saturated_conductivity_mm_h[0] = -1;
    catalog.entries.items[0].material.vertical_hydraulic_conductivity_m2_per_mpa_h[0] = -1;
    var state = try State.initMapped(allocator, &grid, catalog.entries.items, &.{0}, &.{1}, &.{1}, parameters);
    defer state.deinit();
    try std.testing.expectEqual(@as(usize, 37), state.hydraulic_conductivity_class_count);
    try std.testing.expectEqual(grid.layer_count * 3 * 37, state.matrix_hydraulic_conductivity_m2_per_h_mpa.len);
    try std.testing.expect(state.matrix_bulk_volume_m3[0] > 0);
    const entry = catalog.entries.items[0];
    const bulk_soil_mass_Mg = entry.material.bulk_density_megagrams_per_m3[0] * entry.hydrology_per_m2.total_layer_volume_m3[0];
    try std.testing.expectApproxEqAbs(entry.material.sand_mass_fraction[0] * bulk_soil_mass_Mg, state.sand_mass_Mg[0], 1e-14);
    try std.testing.expectApproxEqAbs(entry.material.silt_mass_fraction[0] * bulk_soil_mass_Mg, state.silt_mass_Mg[0], 1e-14);
    try std.testing.expectApproxEqAbs(entry.material.clay_mass_fraction[0] * bulk_soil_mass_Mg, state.clay_mass_Mg[0], 1e-14);
    try std.testing.expectApproxEqAbs(state.cation_exchange_capacity_mol_per_Mg[0] * bulk_soil_mass_Mg, state.cation_exchange_capacity_mol[0], 1e-14);
    try std.testing.expectApproxEqAbs(state.anion_exchange_capacity_mol_per_Mg[0] * bulk_soil_mass_Mg, state.anion_exchange_capacity_mol[0], 1e-14);
    try std.testing.expect(state.matrix_air_entry_water_fraction[0] > 0);
    try std.testing.expect(state.matrix_hydraulic_conductivity_m2_per_h_mpa[2 * parameters.hydraulic_conductivity_class_count] > 0);
    try std.testing.expectApproxEqAbs(0.5 * state.layer_thickness_m[0], state.layer_midpoint_depth_m[0], 1e-15);
    for (1..grid.active_soil_layer_count[0]) |layer| {
        try std.testing.expectApproxEqAbs(state.layer_bottom_depth_m[layer - 1] + 0.5 * state.layer_thickness_m[layer], state.layer_midpoint_depth_m[layer], 1e-15);
        try std.testing.expectApproxEqAbs(state.layer_bottom_depth_m[layer - 1] + state.layer_thickness_m[layer], state.layer_bottom_depth_m[layer], 1e-15);
    }
}

test "HOUR1 missing saturated conductivity uses mineral and organic equations" {
    const parameters = compatibilityParameters();
    const mineral_curve = try retention.resolve(parameters.retention, .{
        .porosity_fraction = 0.5,
        .macropore_fraction = 0,
        .sand_fraction = 0.4,
        .clay_fraction = 0.2,
        .organic_carbon_g_per_megagram = 10_000,
        .bulk_density_megagrams_per_m3 = 1.3,
        .supplied_field_capacity_fraction = 0.25,
        .supplied_wilting_point_fraction = 0.1,
    }, -0.033, -1.5);
    const mineral = try estimateSaturatedConductivity(parameters, mineral_curve, 10_000, 1.3, 0.9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.54), mineral, 1.0e-12);

    const organic = try estimateSaturatedConductivity(parameters, mineral_curve, 300_000, 0.1, 0.8);
    const expected = (0.10 + 75.0 * std.math.pow(f64, 1.0e-15, 0.1)) * 0.8;
    try std.testing.expectApproxEqAbs(expected, organic, 1.0e-12);
}
