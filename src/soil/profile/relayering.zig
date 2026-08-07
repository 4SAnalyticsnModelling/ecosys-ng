const std = @import("std");
const SimulationConfig = @import("../../core/config.zig").SimulationConfig;
const GridState = @import("../../state/grid.zig").GridState;
const thermal_module = @import("../heat/thermal.zig");
const gas = @import("../gas/transport.zig");
const organic = @import("../organic/initialization.zig");
const chemistry_module = @import("../solute/chemistry_state.zig");
const reactive = @import("../nutrients/reactive_nitrogen_state.zig");
const nitrogen_fertilizer = @import("../../management/fertilizer_nitrogen_inventory.zig");
const mineral_fertilizer = @import("../../management/mineral_fertilizer_inventory.zig");
const properties = @import("../water/solver_properties.zig");
const RootSystem = @import("../../plant/root/plant_root_system.zig");
const Geometry = @import("layer_geometry.zig");
const face_geometry_module = @import("../water/face_geometry.zig");
const transport_module = @import("../../transport/hydrology.zig");
const geometry_assembly = @import("geometry_change_assembly.zig");

const water_heat_remap = @import("../water/heat_layer_remap.zig");
const gas_remap = @import("../gas/layer_remap.zig");
const organic_remap = @import("../organic/layer_remap.zig");
const chemistry_remap = @import("../chemistry/layer_remap.zig");
const nitrite_remap = @import("../microbial/nitrite_layer_remap.zig");
const fertilizer_remap = @import("../nutrients/fertilizer_layer_remap.zig");
const mineral_remap = @import("../nutrients/mineral_layer_remap.zig");
const root_remap = @import("../../plant/root/plant_root_layer_remap.zig");

/// Heap-owned scratch for one hourly geometry change transaction.
pub const Workspace = struct {
    allocator: std.mem.Allocator,
    /// Boundary change array for the freeze-thaw driver: cell_count * (layer_capacity + 1).
    freeze_thaw_boundary_change_m: []f64,
    /// Reusable zero-filled boundary change array for inactive drivers.
    zero_boundary_change_m: []f64,
    /// Scratch: per-layer ice volume delta extracted from the water-heat solver result.
    /// Size: cell_count * layer_capacity.
    ice_volume_delta_m3: []f64,
    /// Scratch: matrix zone fraction (matrix_bulk_volume / layer_volume) per layer.
    /// Size: cell_count * layer_capacity.
    matrix_zone_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, layer_capacity: usize) !Workspace {
        if (cell_count == 0 or layer_capacity == 0) return error.InvalidSoilRelayeringDimensions;
        const bnd_count = try std.math.mul(usize, cell_count, try std.math.add(usize, layer_capacity, 1));
        const layer_count = try std.math.mul(usize, cell_count, layer_capacity);
        const ft = try allocator.alloc(f64, bnd_count);
        errdefer allocator.free(ft);
        @memset(ft, 0);
        const zero = try allocator.alloc(f64, bnd_count);
        errdefer allocator.free(zero);
        @memset(zero, 0);
        const ice_delta = try allocator.alloc(f64, layer_count);
        errdefer allocator.free(ice_delta);
        @memset(ice_delta, 0);
        const mzf = try allocator.alloc(f64, layer_count);
        @memset(mzf, 1);
        return .{
            .allocator = allocator,
            .freeze_thaw_boundary_change_m = ft,
            .zero_boundary_change_m = zero,
            .ice_volume_delta_m3 = ice_delta,
            .matrix_zone_fraction = mzf,
        };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.matrix_zone_fraction);
        self.allocator.free(self.ice_volume_delta_m3);
        self.allocator.free(self.zero_boundary_change_m);
        self.allocator.free(self.freeze_thaw_boundary_change_m);
        self.* = undefined;
    }

    pub fn changes(self: *const Workspace) Geometry.DisturbanceChanges {
        return .{
            .pond_m = self.zero_boundary_change_m,
            .freeze_thaw_m = self.freeze_thaw_boundary_change_m,
            .erosion_m = self.zero_boundary_change_m,
            .organic_carbon_m = self.zero_boundary_change_m,
        };
    }
};

pub const Context = struct {
    grid: *GridState,
    soil_thermal: *thermal_module.State,
    gas_transport: *gas.State,
    soil_organic: *organic.State,
    soil_chemistry: *chemistry_module.State,
    reactive_nitrogen: *reactive.State,
    soil_fertilizer_inventory: *nitrogen_fertilizer.State,
    mineral_fertilizer_inventory: *mineral_fertilizer.State,
    /// Optional: when non-null, solid chemistry and mineral transfer are applied
    /// and geometry thickness is synced into properties after applyDisturbances.
    soil_properties: ?*properties.State,
    plant_roots: ?*RootSystem.State,
    soil_geometry: *Geometry.State,
    /// Optional: when non-null (along with soil_transport_faces), face geometry
    /// is refreshed after applyDisturbances.
    soil_face_geometry: ?*face_geometry_module.State,
    soil_transport_faces: ?*const transport_module.SoilFaces,
    water_heat_parameters: water_heat_remap.Parameters,
    dynamic_salts: bool,
    nutrient_zone_fractions: NutrientZoneFractions,
    plant_populations: usize,
    minimum_layer_thickness_m: f64,
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
    /// EXEC-002 remediation switch. When true, the soil-thermal layer volume is
    /// restated onto the committed geometry after the transfers, which improves
    /// hourly heat conservation but currently stagnates the phase-enthalpy
    /// solver. Defaults to false so production behaviour is unchanged.
    rebase_thermal_volume_to_geometry: bool = false,
};

pub const NutrientZoneFractions = struct {
    ammonium_non_band: f64,
    ammonium_band: f64,
    nitrate_non_band: f64,
    nitrate_band: f64,
    phosphate_non_band: f64,
    phosphate_band: f64,
};

/// REDIST DO 245 (NN=1): applies inter-soil-layer pool redistribution for all
/// active cells, then commits the geometry change. Pool remaps run in source
/// order matching Fortran REDIST: water/heat → gas → organic → chemistry
/// (aqueous then solid) → nitrite → fertilizer → mineral → plant roots.
/// The geometry transaction is applied only after all pool transfers complete.
///
/// geometry_changes must cover freeze-thaw, SOC, erosion, and pond boundary
/// changes for the current hour as assembled by soil_geometry_disturbance_transaction
/// or soil_geometry_change_assembly helpers. Zero arrays are acceptable when a
/// driver is inactive.
pub fn applyLayerRedistribution(ctx: Context, geometry_changes: Geometry.DisturbanceChanges) !void {
    if (ctx.grid.cell_count == 0) return;

    // Validate geometry changes before touching any pool state.
    try Geometry.validateDisturbances(ctx.soil_geometry, geometry_changes, ctx.minimum_layer_thickness_m);

    const cap = ctx.grid.soil_layer_capacity;

    for (0..ctx.grid.cell_count) |cell| {
        const active = ctx.grid.active_soil_layer_count[cell];
        if (active < 2) continue;
        const first = ctx.soil_geometry.first_active_layer[cell];
        const bnd_base = cell * (cap + 1);

        for (0..active - 1) |offset| {
            const layer = first + offset;
            const g = cell * cap + layer;

            // DDLYRX: total change at the bottom boundary of layer `layer`.
            const bnd_idx = bnd_base + layer + 1;
            const ddlyrx = geometry_changes.pond_m[bnd_idx] + geometry_changes.freeze_thaw_m[bnd_idx] + geometry_changes.erosion_m[bnd_idx] + geometry_changes.organic_carbon_m[bnd_idx];

            if (@abs(ddlyrx) < ctx.minimum_layer_thickness_m) continue;

            // DDLYRX > 0: layer expands upward, pull material from the layer below.
            // DDLYRX < 0: layer contracts, push material to the layer below.
            const src_global: usize = if (ddlyrx > 0) g + 1 else g;
            const dst_global: usize = if (ddlyrx > 0) g else g + 1;
            const src_layer: usize = if (ddlyrx > 0) layer + 1 else layer;
            const dst_layer: usize = if (ddlyrx > 0) layer else layer + 1;

            const src_thickness = ctx.soil_geometry.layer_thickness_m[src_global];
            if (src_thickness <= ctx.minimum_layer_thickness_m) continue;

            // DISC-WATSUB-002. `redist.f:8468`/`:8498` gate the transfer on
            // `IF(DLYR(3,L0,NY,NX).GT.DLYRM)`, and `DLYR` is mutated *inside* the
            // Fortran relayering loop, so each boundary sees the LIVE thickness
            // left by the previous boundary. In this translation the geometry
            // commit is deferred to `applyDisturbances` after the loop, so
            // `layer_thickness_m` above is the thickness as it stood at the start
            // of the hour for every iteration. The live carrier is
            // `soil_thermal.layer_volume_m3`, which the transfer itself mutates.
            //
            // Measured consequence on `Arctic Fen CH` cell 0, hour 1 (A6 probe):
            // boundary `layer=0` transfers all of layer 1 into layer 0 at
            // `fx=1.0`, leaving layer 1 with `layer_volume_m3=0`,
            // `total_heat_capacity=0` and a meaningless `920.1 K`. Boundary
            // `layer=1` then selects that emptied layer 1 as its DESTINATION,
            // with `geom_dst_thick=2e-2` still reading healthy, and the kernel
            // correctly refuses the pair with `InvalidWaterHeatLayerRemapState`.
            //
            // So this is the stale-carrier half of the defect A7 recorded, and it
            // resolves A7's two open readings in favour of the second: correcting
            // the energy split let an earlier transfer proceed that previously
            // aborted, which is not a regression, and this pair is reached for the
            // first time. Gating on the live carrier is what the source does.
            //
            // Skipping is the conserving choice, not the deferring one: an emptied
            // layer holds no matrix water, ice or vapor to move, so there is
            // nothing to lose. Its retained macropore stores stay attached to it
            // exactly as `redist.f` 9609--9622 requires.
            if (ctx.soil_thermal.layer_volume_m3[src_global] <= 0 or
                ctx.soil_thermal.layer_volume_m3[dst_global] <= 0) continue;

            const fx = @min(1.0, @abs(ddlyrx) / src_thickness);
            if (fx == 0) continue;

            // Save water volumes before water/heat transfer; chemistry needs both.
            const src_water_before = ctx.grid.matrix_liquid_water_m3[src_global];
            const dst_water_before = ctx.grid.matrix_liquid_water_m3[dst_global];

            // 1. Water / heat.
            try water_heat_remap.transferLayerFraction(
                ctx.grid,
                ctx.soil_thermal,
                src_global,
                dst_global,
                fx,
                ctx.water_heat_parameters,
            );

            // 2. Gas (gaseous, dissolved, band-dissolved species).
            //    Band carrier: use updated destination water as conservative proxy;
            //    band fertilizer is not redistributed by freeze-thaw in practice.
            try gas_remap.transferLayerFraction(
                ctx.gas_transport,
                src_global,
                dst_global,
                ctx.grid.matrix_liquid_water_m3[dst_global],
                fx,
            );

            // 3. Organic matter.
            try organic_remap.transferLayerFraction(
                ctx.soil_organic,
                src_global,
                dst_global,
                fx,
            );

            // 4a. Aqueous chemistry (concentration-based, before/after volumes).
            const src_water_after = ctx.grid.matrix_liquid_water_m3[src_global];
            const dst_water_after = ctx.grid.matrix_liquid_water_m3[dst_global];
            const src_zones_before = waterZones(src_water_before, ctx.nutrient_zone_fractions);
            const dst_zones_before = waterZones(dst_water_before, ctx.nutrient_zone_fractions);
            const src_zones_after = waterZones(src_water_after, ctx.nutrient_zone_fractions);
            const dst_zones_after = waterZones(dst_water_after, ctx.nutrient_zone_fractions);
            try chemistry_remap.transferAqueousLayerFraction(
                ctx.soil_chemistry,
                src_global,
                dst_global,
                src_zones_before,
                dst_zones_before,
                src_zones_after,
                dst_zones_after,
                ctx.dynamic_salts,
                fx,
            );

            // 4b. Solid chemistry and 7. Mineral: both need soil mass, which
            // requires soil_properties. Skip when not provided (e.g. in tests).
            if (ctx.soil_properties) |props| {
                // Solid chemistry is defined per Mg of mineral matrix. Total
                // layer volume includes the separate macropore domain and is
                // not its carrier; using it creates adsorbed mass whenever
                // freeze-thaw relayering changes the matrix share.
                const src_mass_before = props.bulk_density_megagrams_per_m3[src_global] * props.matrix_bulk_volume_m3[src_global];
                const dst_mass_before = props.bulk_density_megagrams_per_m3[dst_global] * props.matrix_bulk_volume_m3[dst_global];

                // 4b. Solid chemistry (adsorbed cations, phosphate surfaces, precipitates).
                if (src_mass_before > 0 and dst_mass_before > 0) {
                    const moved_mass = fx * src_mass_before;
                    try chemistry_remap.transferSolidLayerFraction(
                        ctx.soil_chemistry,
                        src_global,
                        dst_global,
                        src_mass_before,
                        dst_mass_before,
                        src_water_before,
                        dst_water_before,
                        src_mass_before - moved_mass,
                        dst_mass_before + moved_mass,
                        src_water_after,
                        dst_water_after,
                        fx,
                    );
                }

                // 7. Mineral / sediment and its physical carrier. The legacy
                // soil branch transfers FX of SAND/SILT/CLAY for every FX > 0.
                const src_mass_after = src_mass_before - fx * src_mass_before;
                const dst_mass_after = dst_mass_before + fx * src_mass_before;
                try mineral_remap.transferLayerFraction(
                    props,
                    src_global,
                    dst_global,
                    fx,
                    src_mass_after,
                    dst_mass_after,
                );
                const moved_matrix_volume_m3 = fx * props.matrix_bulk_volume_m3[src_global];
                props.matrix_bulk_volume_m3[src_global] -= moved_matrix_volume_m3;
                props.matrix_bulk_volume_m3[dst_global] += moved_matrix_volume_m3;
                if (props.matrix_bulk_volume_m3[src_global] > 0)
                    props.bulk_density_megagrams_per_m3[src_global] = src_mass_after / props.matrix_bulk_volume_m3[src_global];
                props.bulk_density_megagrams_per_m3[dst_global] = dst_mass_after / props.matrix_bulk_volume_m3[dst_global];
            }

            // 5. Nitrite.
            try nitrite_remap.transferLayerFraction(
                ctx.reactive_nitrogen,
                src_global,
                dst_global,
                ctx.grid.matrix_liquid_water_m3[dst_global],
                fx,
            );

            // 6. Fertilizer (cell-relative layer indices).
            try fertilizer_remap.transferCellLayerFraction(
                ctx.soil_fertilizer_inventory,
                ctx.mineral_fertilizer_inventory,
                cell,
                src_layer,
                dst_layer,
                fx,
            );

            // 8. Plant roots (loop over all plants in this cell).
            if (ctx.plant_roots) |roots| {
                for (0..ctx.plant_populations) |species| {
                    const plant = cell * ctx.plant_populations + species;
                    try root_remap.transferLayerFraction(roots, plant, src_layer, dst_layer, fx);
                }
            }
        }
    }

    // Commit geometry after all pool transfers succeed.
    try Geometry.applyDisturbances(ctx.soil_geometry, geometry_changes, ctx.minimum_layer_thickness_m);
    // EXEC-002: `soil_thermal.layer_volume_m3` drifts from the committed
    // `layer_thickness_m * cell_area_m2` here. The drift is confined to the top
    // active layer, whose upper face is a free surface with no paired transfer
    // (2.8e-2 relative there, 5.3e-16 for every deeper layer). Because the EXEC
    // heat census forms extensive capacity as `dry_capacity_per_m3 * volume`,
    // that drift reads as created heat.
    //
    // `rebaseThermalVolumeToGeometry` corrects it and is tested, but enabling it
    // measured *worse* end to end: the day-1 audit deviation rises from
    // 7.075e-2 to 1.636e-1 per m2 even with the companion
    // `soil_properties.layer_volume_m3` update below. Correcting the volume
    // without also re-deriving the dry solid capacity from soil material moves
    // the census and the water/heat commit onto different footings. It therefore
    // stays off. See the EXEC-002 entry in `docs/discrepancy_register.md`.
    if (ctx.rebase_thermal_volume_to_geometry) try rebaseThermalVolumeToGeometry(ctx, cap);
    if (ctx.soil_properties) |props| {
        @memcpy(props.layer_thickness_m, ctx.soil_geometry.layer_thickness_m);
        @memcpy(props.layer_midpoint_depth_m, ctx.soil_geometry.layer_midpoint_depth_from_surface_m);
        @memcpy(props.layer_bottom_depth_m, ctx.soil_geometry.layer_bottom_depth_from_surface_m);
        // `matrix_bulk_volume_m3` was moved between layers above, but
        // `layer_volume_m3` was left at its initialization value. Consumers
        // divide one by the other: `matrix_zone_fraction` and
        // `soil_matrix_fraction` in `ecosys_ng.zig` both do, and the latter
        // feeds the next hour's freeze-thaw boundary change. Leaving them
        // inconsistent lets a transfer silently change a fraction that is
        // supposed to describe the same layer. Restate the total volume from the
        // committed geometry so the ratio stays meaningful.
        if (ctx.rebase_thermal_volume_to_geometry) {
            for (0..ctx.grid.cell_count) |cell| {
                const area_m2 = ctx.horizontal_cell_width_m[cell] * ctx.vertical_cell_width_m[cell];
                if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidRelayeringCellArea;
                const active = ctx.grid.active_soil_layer_count[cell];
                const first = ctx.soil_geometry.first_active_layer[cell];
                for (0..active) |offset| {
                    const g = cell * cap + first + offset;
                    const volume_m3 = ctx.soil_geometry.layer_thickness_m[g] * area_m2;
                    if (!std.math.isFinite(volume_m3) or volume_m3 <= 0) return error.InvalidRelayeringLayerVolume;
                    props.layer_volume_m3[g] = volume_m3;
                }
            }
        }
    }
    if (ctx.soil_face_geometry) |fg| if (ctx.soil_transport_faces) |tf| {
        fg.refreshMapped(
            ctx.grid,
            tf,
            ctx.soil_geometry.layer_thickness_m,
            ctx.horizontal_cell_width_m,
            ctx.vertical_cell_width_m,
        ) catch unreachable;
    };
}

/// Restates the soil-thermal layer volume onto the committed geometry without
/// moving energy. The extensive dry and total heat capacities are held fixed, so
/// the EXEC heat census sees the same energy before and after; only the
/// per-volume densities and the porosity fraction are restated.
///
/// Measured scope: after `applyLayerRedistribution` only the top active layer
/// drifts materially from the committed geometry. Every deeper layer agrees to
/// within f64 rounding (worst observed 5.3e-16 relative, versus 2.8e-2 for the
/// top layer). That is structural, not incidental: the transfer loop pairs each
/// layer with the one below it through the boundary at `layer + 1`, and the
/// freeze-thaw assembly anchors the bottom datum and accumulates upward, so the
/// top layer's own upper face is a free surface with no paired transfer. This
/// matches REDIST, which at `LX.EQ.NU` propagates the cumulative change to the
/// boundary above the top layer (`redist.f` lines 7973--7977).
///
/// Only that free surface is rebased. Restating the deeper layers as well was
/// tried and measured indistinguishable, confirming the top layer carries the
/// whole effect.
///
/// This is not enabled in production: see the call site and the EXEC-002 entry
/// in `docs/discrepancy_register.md` for the end-to-end measurement that keeps
/// it off.
fn rebaseThermalVolumeToGeometry(ctx: Context, cap: usize) !void {
    for (0..ctx.grid.cell_count) |cell| {
        const area_m2 = ctx.horizontal_cell_width_m[cell] * ctx.vertical_cell_width_m[cell];
        if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidRelayeringCellArea;
        const active = ctx.grid.active_soil_layer_count[cell];
        if (active == 0) continue;
        const first = ctx.soil_geometry.first_active_layer[cell];
        // The top active layer is the only free surface in this transaction.
        const g = cell * cap + first;
        const geometry_volume_m3 = ctx.soil_geometry.layer_thickness_m[g] * area_m2;
        if (!std.math.isFinite(geometry_volume_m3) or geometry_volume_m3 <= 0) return error.InvalidRelayeringLayerVolume;
        const previous_volume_m3 = ctx.soil_thermal.layer_volume_m3[g];
        if (!std.math.isFinite(previous_volume_m3) or previous_volume_m3 <= 0) return error.InvalidRelayeringLayerVolume;
        if (previous_volume_m3 == geometry_volume_m3) continue;
        const extensive_dry = ctx.soil_thermal.dry_solid_heat_capacity_megajoules_per_m3_k[g] * previous_volume_m3;
        const extensive_total = ctx.soil_thermal.total_heat_capacity_megajoules_per_m3_k[g] * previous_volume_m3;
        ctx.soil_thermal.layer_volume_m3[g] = geometry_volume_m3;
        // Both capacities preserve their extensive value, so the census energy is
        // unchanged by the correction. Holding the dry density constant instead,
        // on the theory that it is a material property of unchanged solids, was
        // tried and measured clearly worse (day-1 deviation 4.120e-1 versus
        // 1.636e-1 per m2), so the extensive reading is the one kept.
        ctx.soil_thermal.dry_solid_heat_capacity_megajoules_per_m3_k[g] = extensive_dry / geometry_volume_m3;
        ctx.soil_thermal.total_heat_capacity_megajoules_per_m3_k[g] = extensive_total / geometry_volume_m3;
        // Porosity is not an independent store: recompute it from the live pore
        // capacities against the new volume, exactly as
        // `soil_water_heat_layer_remap.commitDerived` and `soil_thermal` do.
        ctx.soil_thermal.porosity_fraction[g] =
            (ctx.grid.matrix_pore_capacity_m3[g] + ctx.grid.macropore_pore_capacity_m3[g]) / geometry_volume_m3;
        inline for (.{
            ctx.soil_thermal.dry_solid_heat_capacity_megajoules_per_m3_k[g],
            ctx.soil_thermal.total_heat_capacity_megajoules_per_m3_k[g],
            ctx.soil_thermal.porosity_fraction[g],
        }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRelayeringThermalRebase;
    }
}

/// Build the distinct runtime nutrient carriers used by concentration state.
fn waterZones(shared_m3: f64, fractions: NutrientZoneFractions) chemistry_remap.ZoneWaterVolumes {
    return .{
        .shared_m3 = shared_m3,
        .ammonium_non_band_m3 = shared_m3 * fractions.ammonium_non_band,
        .ammonium_band_m3 = shared_m3 * fractions.ammonium_band,
        .nitrate_non_band_m3 = shared_m3 * fractions.nitrate_non_band,
        .nitrate_band_m3 = shared_m3 * fractions.nitrate_band,
        .phosphate_non_band_m3 = shared_m3 * fractions.phosphate_non_band,
        .phosphate_band_m3 = shared_m3 * fractions.phosphate_band,
    };
}

/// REDIST DO 225 freeze-thaw geometry: assembles boundary changes from the
/// accepted per-layer ice volume delta (DVOLI). Wraps
/// soil_geometry_change_assembly.assembleFreezeThawBoundaryChangeM with a
/// workspace already sized to the geometry.
pub fn assembleFreezeThawChanges(
    workspace: *Workspace,
    geometry: *const Geometry.State,
    total_ice_volume_change_m3: []const f64,
    soil_matrix_fraction: []const f64,
    horizontal_area_m2_by_cell: []const f64,
    ice_to_water_specific_volume_difference: f64,
    negligible_ice_volume_change_m3: f64,
) !void {
    @memset(workspace.freeze_thaw_boundary_change_m, 0);
    try geometry_assembly.assembleFreezeThawBoundaryChangeM(
        workspace.freeze_thaw_boundary_change_m,
        geometry,
        total_ice_volume_change_m3,
        soil_matrix_fraction,
        horizontal_area_m2_by_cell,
        ice_to_water_specific_volume_difference,
        negligible_ice_volume_change_m3,
    );
}

test "layer redistribution conserves water between adjacent soil layers" {
    const allocator = std.testing.allocator;

    var geometry = try Geometry.State.init(allocator, 1, 3);
    defer geometry.deinit();
    try Geometry.initializeCell(&geometry, 0, 0, &.{ 0.1, 0.2, 0.3 }, 0, 1e-9);

    const cfg = try SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 3, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-6, .absolute_tolerance = 1e-6, .max_nonlinear_iterations = 10 },
    );
    var grid = try GridState.init(allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.05;
    grid.matrix_liquid_water_m3[1] = 0.10;
    grid.matrix_liquid_water_m3[2] = 0.15;
    grid.soil_temperature_k[0] = 280;
    grid.soil_temperature_k[1] = 278;
    grid.soil_temperature_k[2] = 275;
    // active_soil_layer_count already set to 3 by GridState.init

    // Thermal state — use undefined + stack slices (no allocator needed).
    var thermal_lv = [_]f64{ 0.05, 0.10, 0.15 };
    var thermal_dry = [_]f64{ 1.0, 1.0, 1.0 };
    var thermal_total = [_]f64{ 2.0, 2.0, 2.0 };
    var thermal_porosity = [_]f64{ 0.4, 0.4, 0.4 };
    var thermal_state: thermal_module.State = undefined;
    thermal_state.layer_volume_m3 = &thermal_lv;
    thermal_state.dry_solid_heat_capacity_megajoules_per_m3_k = &thermal_dry;
    thermal_state.total_heat_capacity_megajoules_per_m3_k = &thermal_total;
    thermal_state.porosity_fraction = &thermal_porosity;

    const total_water_before = grid.matrix_liquid_water_m3[0] + grid.matrix_liquid_water_m3[1] + grid.matrix_liquid_water_m3[2];

    // Simulate layer 0 expanding (bottom boundary pushes down 0.02 m):
    // DDLYRX > 0 at boundary between layer 0 and layer 1 → pull from layer 1.
    var ft_change = [_]f64{0} ** 4; // 1 cell × (3+1) boundaries
    ft_change[1] = 0.02; // boundary between layer 0 and layer 1 moves down

    var zero_change = [_]f64{0} ** 4;
    const geometry_changes = Geometry.DisturbanceChanges{
        .pond_m = &zero_change,
        .freeze_thaw_m = &ft_change,
        .erosion_m = &zero_change,
        .organic_carbon_m = &zero_change,
    };

    var gas_state = try gas.State.init(allocator, 3);
    defer gas_state.deinit();
    var organic_state = try organic.State.init(allocator, 3);
    defer organic_state.deinit();
    var chem_state = try chemistry_module.State.init(allocator, 3);
    defer chem_state.deinit();
    var reactive_state = try reactive.State.init(allocator, 3, 1);
    defer reactive_state.deinit();
    var nfert_state = try nitrogen_fertilizer.State.init(allocator, 1, 3);
    defer nfert_state.deinit();
    var mfert_state = try mineral_fertilizer.State.init(allocator, 1, 3);
    defer mfert_state.deinit();

    try applyLayerRedistribution(.{
        .grid = &grid,
        .soil_thermal = &thermal_state,
        .gas_transport = &gas_state,
        .soil_organic = &organic_state,
        .soil_chemistry = &chem_state,
        .reactive_nitrogen = &reactive_state,
        .soil_fertilizer_inventory = &nfert_state,
        .mineral_fertilizer_inventory = &mfert_state,
        .soil_properties = null,
        .plant_roots = null,
        .soil_geometry = &geometry,
        .soil_face_geometry = null,
        .soil_transport_faces = null,
        .water_heat_parameters = .{
            .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
            .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
            .minimum_heat_capacity_megajoules_per_k = 1e-6,
        },
        .dynamic_salts = false,
        .nutrient_zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 },
        .plant_populations = 0,
        .minimum_layer_thickness_m = 1e-9,
        .horizontal_cell_width_m = &.{1},
        .vertical_cell_width_m = &.{1},
    }, geometry_changes);

    const total_water_after = grid.matrix_liquid_water_m3[0] + grid.matrix_liquid_water_m3[1] + grid.matrix_liquid_water_m3[2];
    try std.testing.expectApproxEqAbs(total_water_before, total_water_after, 1e-14);
    // Layer 0 expanded (pulled from layer 1), so layer 1 loses water and layer 0 gains.
    try std.testing.expect(grid.matrix_liquid_water_m3[0] > 0.05);
    try std.testing.expect(grid.matrix_liquid_water_m3[1] < 0.10);

    // Without the EXEC-002 rebase the thermal volume no longer matches the
    // committed geometry, which is exactly the drift the census sees.
    var worst_relative_mismatch: f64 = 0;
    for (0..3) |layer| {
        const geometry_volume_m3 = geometry.layer_thickness_m[layer] * 1.0;
        const scale = @max(@abs(geometry_volume_m3), @abs(thermal_state.layer_volume_m3[layer]));
        if (scale > 0)
            worst_relative_mismatch = @max(worst_relative_mismatch, @abs(thermal_state.layer_volume_m3[layer] - geometry_volume_m3) / scale);
    }
    try std.testing.expect(worst_relative_mismatch > 1e-6);
}

test "EXEC-002 rebase aligns thermal volume with geometry and preserves census energy" {
    const allocator = std.testing.allocator;
    const config = try @import("../../core/config.zig").SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 3, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 4 },
    );
    var grid = try GridState.init(allocator, config);
    defer grid.deinit();

    var geometry = try Geometry.State.init(allocator, 1, 3);
    defer geometry.deinit();
    try Geometry.initializeCell(&geometry, 0, 0, &.{ 0.1, 0.2, 0.3 }, 0, 1e-9);

    var thermal_state: thermal_module.State = undefined;
    var thermal_volume = [_]f64{ 0.1, 0.2, 0.3 };
    var thermal_dry = [_]f64{ 1.0, 1.1, 1.2 };
    var thermal_total = [_]f64{ 2.0, 2.1, 2.2 };
    var thermal_porosity = [_]f64{ 0.4, 0.4, 0.4 };
    thermal_state.layer_volume_m3 = &thermal_volume;
    thermal_state.dry_solid_heat_capacity_megajoules_per_m3_k = &thermal_dry;
    thermal_state.total_heat_capacity_megajoules_per_m3_k = &thermal_total;
    thermal_state.porosity_fraction = &thermal_porosity;
    for (0..3) |layer| {
        grid.soil_temperature_k[layer] = 275 + @as(f64, @floatFromInt(layer));
        grid.matrix_pore_capacity_m3[layer] = 0.4 * thermal_volume[layer];
        grid.macropore_pore_capacity_m3[layer] = 0;
    }

    // Detach the top layer's volume from geometry the way the free surface does.
    // Only that layer is rebased, so the census energy of layer 0 must be
    // preserved while layers 1 and 2 are left untouched.
    thermal_volume[0] += 0.02;
    const untouched_volume_1 = thermal_volume[1];
    const untouched_volume_2 = thermal_volume[2];
    var census_energy_before: f64 = 0;
    for (0..3) |layer| census_energy_before += thermal_total[layer] * thermal_volume[layer] * grid.soil_temperature_k[layer];

    try rebaseThermalVolumeToGeometry(.{
        .grid = &grid,
        .soil_thermal = &thermal_state,
        .gas_transport = undefined,
        .soil_organic = undefined,
        .soil_chemistry = undefined,
        .reactive_nitrogen = undefined,
        .soil_fertilizer_inventory = undefined,
        .mineral_fertilizer_inventory = undefined,
        .soil_properties = null,
        .plant_roots = null,
        .soil_geometry = &geometry,
        .soil_face_geometry = null,
        .soil_transport_faces = null,
        .water_heat_parameters = .{
            .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
            .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
            .minimum_heat_capacity_megajoules_per_k = 1e-6,
        },
        .dynamic_salts = false,
        .nutrient_zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 },
        .plant_populations = 0,
        .minimum_layer_thickness_m = 1e-9,
        .horizontal_cell_width_m = &.{1},
        .vertical_cell_width_m = &.{1},
    }, 3);

    // The top layer is realigned with geometry.
    try std.testing.expectApproxEqRel(geometry.layer_thickness_m[0], thermal_volume[0], 1e-15);
    // Deeper layers are deliberately untouched: their volume already matched.
    try std.testing.expectEqual(untouched_volume_1, thermal_volume[1]);
    try std.testing.expectEqual(untouched_volume_2, thermal_volume[2]);
    var census_energy_after: f64 = 0;
    for (0..3) |layer| census_energy_after += thermal_total[layer] * thermal_volume[layer] * grid.soil_temperature_k[layer];
    try std.testing.expectApproxEqRel(census_energy_before, census_energy_after, 1e-12);
    // Porosity of the rebased layer is recomputed from live pore capacity.
    try std.testing.expectApproxEqRel(grid.matrix_pore_capacity_m3[0] / thermal_volume[0], thermal_porosity[0], 1e-12);
}

test "DISC-WATSUB-002 relayering skips a boundary whose layer was emptied earlier in the same hour" {
    // Measured on `Arctic Fen CH` cell 0, hour 1: boundary `layer=0` transferred
    // ALL of layer 1 into layer 0 at `fx=1.0`, and boundary `layer=1` then chose
    // that emptied layer 1 as its DESTINATION. The geometry thickness still read
    // a healthy `2e-2` because the geometry commit is deferred to after the loop,
    // while the live `soil_thermal.layer_volume_m3` was already `0` and the layer
    // temperature had become a meaningless `9.2e2 K`.
    //
    // `redist.f:8468`/`:8498` gate on `DLYR(3,L0)`, which the Fortran mutates
    // inside the loop, so the source never presents such a pair. This test pins
    // the live-carrier gate: with an emptied layer in the chain the transaction
    // must complete rather than fail, and must not resurrect the emptied layer.
    const allocator = std.testing.allocator;

    var geometry = try Geometry.State.init(allocator, 1, 3);
    defer geometry.deinit();
    try Geometry.initializeCell(&geometry, 0, 0, &.{ 0.02, 0.025, 0.3 }, 0, 1e-9);

    const cfg = try SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 3, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-6, .absolute_tolerance = 1e-6, .max_nonlinear_iterations = 10 },
    );
    var grid = try GridState.init(allocator, cfg);
    defer grid.deinit();

    var thermal_lv = [_]f64{ 1.0, 0.0, 1.5 }; // layer 1 already emptied
    var thermal_dry = [_]f64{ 1.0, 0.0, 1.0 };
    var thermal_total = [_]f64{ 2.0, 0.0, 2.0 };
    var thermal_porosity = [_]f64{ 0.4, 0.0, 0.4 };
    var thermal_state: thermal_module.State = undefined;
    thermal_state.layer_volume_m3 = &thermal_lv;
    thermal_state.dry_solid_heat_capacity_megajoules_per_m3_k = &thermal_dry;
    thermal_state.total_heat_capacity_megajoules_per_m3_k = &thermal_total;
    thermal_state.porosity_fraction = &thermal_porosity;

    grid.matrix_liquid_water_m3[0] = 0.05;
    grid.matrix_liquid_water_m3[1] = 0; // emptied
    grid.matrix_liquid_water_m3[2] = 0.15;
    grid.soil_temperature_k[0] = 275;
    // The nonsense temperature the emptied layer actually carried. It must not be
    // consumed, and it must not be overwritten by a transfer either.
    grid.soil_temperature_k[1] = 920.1111879438225;
    grid.soil_temperature_k[2] = 270;

    // A boundary change at the layer1/layer2 face, which would select the emptied
    // layer 1 as one endpoint.
    var ft_change = [_]f64{0} ** 4;
    ft_change[2] = -0.01;
    var zero_change = [_]f64{0} ** 4;

    var gas_state = try gas.State.init(allocator, 3);
    defer gas_state.deinit();
    var organic_state = try organic.State.init(allocator, 3);
    defer organic_state.deinit();
    var chem_state = try chemistry_module.State.init(allocator, 3);
    defer chem_state.deinit();
    var reactive_state = try reactive.State.init(allocator, 3, 1);
    defer reactive_state.deinit();
    var nfert_state = try nitrogen_fertilizer.State.init(allocator, 1, 3);
    defer nfert_state.deinit();
    var mfert_state = try mineral_fertilizer.State.init(allocator, 1, 3);
    defer mfert_state.deinit();

    const water_before = grid.matrix_liquid_water_m3[0] + grid.matrix_liquid_water_m3[1] + grid.matrix_liquid_water_m3[2];

    // Before the live-carrier gate this returned InvalidWaterHeatLayerRemapState.
    try applyLayerRedistribution(.{
        .grid = &grid,
        .soil_thermal = &thermal_state,
        .gas_transport = &gas_state,
        .soil_organic = &organic_state,
        .soil_chemistry = &chem_state,
        .reactive_nitrogen = &reactive_state,
        .soil_fertilizer_inventory = &nfert_state,
        .mineral_fertilizer_inventory = &mfert_state,
        .soil_properties = null,
        .plant_roots = null,
        .soil_geometry = &geometry,
        .soil_face_geometry = null,
        .soil_transport_faces = null,
        .water_heat_parameters = .{
            .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
            .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
            .minimum_heat_capacity_megajoules_per_k = 1e-6,
        },
        .dynamic_salts = false,
        .nutrient_zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 },
        .plant_populations = 0,
        .minimum_layer_thickness_m = 1e-9,
        .horizontal_cell_width_m = &.{1},
        .vertical_cell_width_m = &.{1},
    }, .{
        .pond_m = &zero_change,
        .freeze_thaw_m = &ft_change,
        .erosion_m = &zero_change,
        .organic_carbon_m = &zero_change,
    });

    // Skipping conserves: an emptied layer has no matrix water to move, so
    // nothing is lost by declining the pair.
    const water_after = grid.matrix_liquid_water_m3[0] + grid.matrix_liquid_water_m3[1] + grid.matrix_liquid_water_m3[2];
    try std.testing.expectApproxEqAbs(water_before, water_after, 1e-14);
    // The emptied layer stays empty rather than being resurrected with a share of
    // a neighbour's water against a zero heat capacity.
    try std.testing.expectEqual(@as(f64, 0), grid.matrix_liquid_water_m3[1]);
    try std.testing.expectEqual(@as(f64, 0), thermal_state.layer_volume_m3[1]);
}
