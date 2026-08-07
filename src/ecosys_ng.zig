//! ecosys-ng entry point.
//!
//! This file holds `main` and nothing else. The science helpers it
//! drives live in `src/stages/`, grouped by role; they are aliased back
//! into scope below so call sites inside `main` read unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");

const tile_kernels = @import("stages/tile_kernels.zig");
const plant_daily = @import("stages/plant_daily.zig");
const root_processes = @import("stages/root_processes.zig");
const soil_chemistry_convergence = @import("stages/soil_chemistry_convergence.zig");
const biogeochemistry_batches = @import("stages/biogeochemistry_batches.zig");
const surface_litter_convergence = @import("stages/surface_litter_convergence.zig");
const hourly_science = @import("stages/hourly_science.zig");
const diagnostics = @import("stages/diagnostics.zig");
const run_support = @import("stages/run_support.zig");

const DailyPlantElement = plant_daily.DailyPlantElement;
const applyRootMetabolism = root_processes.applyRootMetabolism;
const assemblePlantBalanceInputs = plant_daily.assemblePlantBalanceInputs;
const calculateDailyPlantElementPools = plant_daily.calculateDailyPlantElementPools;
const convergeSurfaceLitterChemistry = surface_litter_convergence.convergeSurfaceLitterChemistry;
const dailySoilWaterPotentialLayerCount = run_support.dailySoilWaterPotentialLayerCount;
const dayOfYearFromTimestamp = run_support.dayOfYearFromTimestamp;
const executeHourlyScience = hourly_science.executeHourlyScience;
const openOsTempDir = run_support.openOsTempDir;
const outputSpeciesLabel = run_support.outputSpeciesLabel;
const plantOutputCatalog = run_support.plantOutputCatalog;
const reconstructLandscapeMassBalance = diagnostics.reconstructLandscapeMassBalance;
const resolveInputPath = run_support.resolveInputPath;
const runKernelAcrossSerialTilePlan = tile_kernels.runKernelAcrossSerialTilePlan;
const sameCalendarDay = run_support.sameCalendarDay;
const sameWeatherTimestamp = run_support.sameWeatherTimestamp;
const soilOutputCatalog = run_support.soilOutputCatalog;

/// Seeds the initial soil profile from the site's mean annual soil temperature.
///
/// `starts.f` lines 397--400 and 1222 seed the profile isothermally at
/// `ATKS = ATCS + 273.15`, where `ATCS` comes from the site's mean annual air
/// temperature. `GridState` allocates `soil_temperature_k` at 273.15 K as a
/// placeholder, so without this every layer would begin exactly on the
/// pure-water freezing point, the worst possible start for a phase-enthalpy
/// solver.
///
/// Enabling this reduced the day-1 heat imbalance from `-3.913e8` MJ to
/// `-2.695e7` MJ, a 14.5x improvement, and cut the share carried by the first
/// three hours from 67.3% to 8.3%. `profile_relayering` also drops to exactly
/// zero for all 24 hours, confirming its earlier large deltas were an artifact
/// of the freezing-point start rather than a redistribution defect. See
/// EXEC-002.
const seed_soil_profile_from_mean_annual = true;

pub fn main(init: std.process.Init) !void {
    var debug_allocator: std.heap.DebugAllocator(.{}) = .init;
    defer _ = debug_allocator.deinit();
    const allocator = debug_allocator.allocator();

    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.log.err("usage: ecosys_ng <runscript>", .{});
        return error.MissingRunscriptPath;
    }
    const source = try std.Io.Dir.cwd().readFileAlloc(init.io, args[1], allocator, .limited(16 * 1024 * 1024));
    defer allocator.free(source);
    var runscript = try ecosys.runscript.parse(allocator, source);
    defer runscript.deinit();
    if (runscript.geospatial_bounds == null) {
        std.log.err("missing compulsory geospatial_grid record: expected minimum/maximum latitude, minimum/maximum longitude, and latitude/longitude intervals", .{});
        return error.MissingGeospatialGridRecord;
    }
    std.log.warn("runtime max_nonlinear_iterations is recorded, but process solvers use option-derived NPH/NPG or legacy NPR/NPS/NPRS iteration ceilings", .{});
    const runscript_directory = std.fs.path.dirname(args[1]) orelse ".";
    const organic_parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.organic_initialization_file);
    defer allocator.free(organic_parameter_path);
    const organic_parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, organic_parameter_path, allocator, .limited(16 * 1024 * 1024));
    defer allocator.free(organic_parameter_source);
    var organic_parameters = try ecosys.soil_organic_parameters.parse(allocator, organic_parameter_source);
    defer organic_parameters.deinit();
    const surface_gas_parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.surface_gas_parameter_file);
    defer allocator.free(surface_gas_parameter_path);
    const surface_gas_parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, surface_gas_parameter_path, allocator, .limited(1024 * 1024));
    defer allocator.free(surface_gas_parameter_source);
    var surface_gas_parameters = try ecosys.surface_gas_parameters.parse(surface_gas_parameter_source);
    std.log.info("loaded runtime surface gas and microbial oxygen parameters", .{});
    const soil_nitrogen_parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.soil_nitrogen_parameter_file);
    defer allocator.free(soil_nitrogen_parameter_path);
    const soil_nitrogen_parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, soil_nitrogen_parameter_path, allocator, .limited(1024 * 1024));
    defer allocator.free(soil_nitrogen_parameter_source);
    const soil_nitrogen_parameters = try ecosys.soil_nitrogen_parameters.parse(soil_nitrogen_parameter_source);
    std.log.info("loaded runtime layer-resolved soil nitrogen parameters", .{});
    const chemistry_reaction_parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.chemistry_reaction_file);
    defer allocator.free(chemistry_reaction_parameter_path);
    const chemistry_reaction_parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, chemistry_reaction_parameter_path, allocator, .limited(4 * 1024 * 1024));
    defer allocator.free(chemistry_reaction_parameter_source);
    const chemistry_reaction_parameters = try ecosys.soil_chemistry_parameters.parse(chemistry_reaction_parameter_source);
    const grid_input_path = try resolveInputPath(
        allocator,
        init.io,
        runscript_directory,
        runscript.grid_input_file,
    );
    defer allocator.free(grid_input_path);
    const grid_input_source = try std.Io.Dir.cwd().readFileAlloc(
        init.io,
        grid_input_path,
        allocator,
        .limited(16 * 1024 * 1024),
    );
    defer allocator.free(grid_input_source);
    var grid_input_files = try ecosys.grid_input_files.CellFiles.parse(
        allocator,
        grid_input_source,
        try runscript.domain.columns(),
        try runscript.domain.rows(),
    );
    defer grid_input_files.deinit();

    var site_catalog = ecosys.site_catalog.Catalog.init(allocator);
    defer site_catalog.deinit();
    var topography_catalog = ecosys.topography_catalog.Catalog.init(allocator);
    defer topography_catalog.deinit();
    for (grid_input_files.site_file_by_cell) |site_name| {
        if (site_catalog.find(site_name) != null) continue;
        const site_path = try resolveInputPath(allocator, init.io, runscript_directory, site_name);
        defer allocator.free(site_path);
        const site_source = try std.Io.Dir.cwd().readFileAlloc(
            init.io,
            site_path,
            allocator,
            .limited(16 * 1024 * 1024),
        );
        defer allocator.free(site_source);
        _ = site_catalog.appendFromSource(
            site_name,
            site_source,
            1,
            1,
        ) catch |err| {
            std.log.err("site input validation failed: file='{s}' error={s}", .{ site_name, @errorName(err) });
            return err;
        };
    }
    for (grid_input_files.topography_file_by_cell) |topography_name| {
        if (topography_catalog.find(topography_name) != null) continue;
        const topography_path = try resolveInputPath(allocator, init.io, runscript_directory, topography_name);
        defer allocator.free(topography_path);
        const topography_source = try std.Io.Dir.cwd().readFileAlloc(
            init.io,
            topography_path,
            allocator,
            .limited(64 * 1024 * 1024),
        );
        defer allocator.free(topography_source);
        _ = topography_catalog.appendFromSource(topography_name, topography_source) catch |err| {
            std.log.err("topography input validation failed: file='{s}' error={s}", .{ topography_name, @errorName(err) });
            return err;
        };
    }
    var grid_environment = try ecosys.grid_environment.Assignments.init(
        allocator,
        runscript.domain,
        grid_input_files,
        site_catalog,
        topography_catalog,
    );
    defer grid_environment.deinit();
    const lateral_connection_mode_by_cell = try allocator.alloc(
        u8,
        grid_environment.site_catalog_index_by_cell.len,
    );
    defer allocator.free(lateral_connection_mode_by_cell);
    for (lateral_connection_mode_by_cell, 0..) |*mode, cell| {
        mode.* = @intFromEnum(site_catalog.entries.items[
            grid_environment.site_catalog_index_by_cell[cell]
        ].site.lateral_connection_mode);
    }
    var geospatial_grid: ?ecosys.spatial_grid.RegularGrid = if (runscript.geospatial_bounds) |bounds|
        try ecosys.spatial_grid.RegularGrid.init(allocator, bounds)
    else
        null;
    defer if (geospatial_grid) |*grid| grid.deinit();
    var tile_plan: ?ecosys.spatial_grid.TilePlan = null;
    defer if (tile_plan) |*plan| plan.deinit();
    if (geospatial_grid) |grid| {
        const cell_count = try grid.cellCount();
        if (cell_count != grid_environment.horizontal_cell_width_m.len)
            return error.GeospatialGridEnvironmentDimensionMismatch;
        for (0..cell_count) |cell| {
            const selected_site = site_catalog.entries.items[
                grid_environment.site_catalog_index_by_cell[cell]
            ].site;
            try grid.validateSiteCoordinate(
                cell,
                selected_site.latitude_degrees_north,
                selected_site.longitude_degrees_east,
            );
            // Geospatial dimensions are authoritative. Site widths remain
            // available only to exterior-boundary fallbacks.
            grid_environment.horizontal_cell_width_m[cell] =
                grid.east_west_cell_width_m[cell];
            grid_environment.vertical_cell_width_m[cell] =
                grid.north_south_cell_width_m[cell];
        }
        try ecosys.spatial_grid.validateNeighborExchangeControls(
            grid.row_count,
            grid.column_count,
            lateral_connection_mode_by_cell,
        );
        tile_plan = try ecosys.spatial_grid.TilePlan.init(
            allocator,
            grid.row_count,
            grid.column_count,
            runscript.tile_row_count,
            runscript.tile_column_count,
            runscript.lateral_flow_halo_cell_count,
        );
    }
    if (tile_plan == null) {
        tile_plan = try ecosys.spatial_grid.TilePlan.init(
            allocator,
            try runscript.domain.rows(),
            try runscript.domain.columns(),
            runscript.tile_row_count,
            runscript.tile_column_count,
            runscript.lateral_flow_halo_cell_count,
        );
    }
    const site = site_catalog.entries.items[grid_environment.site_catalog_index_by_cell[0]].site;
    // Atmospheric composition is environmental state owned by the site
    // input, never a process-parameter default. Supplying the compulsory
    // surface-gas coefficient file must not disable this conversion.
    surface_gas_parameters.atmospheric_concentration_g_per_m3 = try ecosys.surface_gas_parameters.atmosphericConcentrationsFromUmolPerMol(.{
        site.atmospheric_co2_umol_mol,
        site.atmospheric_methane_umol_mol,
        site.atmospheric_oxygen_umol_mol,
        site.atmospheric_nitrogen_umol_mol,
        site.atmospheric_nitrous_oxide_umol_mol,
        site.atmospheric_ammonia_umol_mol,
        0,
    }, site.mean_annual_air_temperature_c + 273.15);

    var topography = try grid_environment.buildDomainTopography(
        runscript.domain,
        grid_input_files,
        topography_catalog,
    );
    defer topography.deinit();
    const unit_by_cell = try topography.buildCellUnitMap(runscript.domain, allocator);
    defer allocator.free(unit_by_cell);
    var terrain_hydrology_state = try ecosys.terrain_hydrology.State.initMapped(
        allocator,
        topography,
        unit_by_cell,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
        try runscript.domain.columns(),
        try runscript.domain.rows(),
    );
    defer terrain_hydrology_state.deinit();
    var surface_runoff_state = try ecosys.surface_runoff.State.init(allocator, unit_by_cell.len);
    defer surface_runoff_state.deinit();
    var surface_pond_transition_state = try ecosys.surface_pond_transition_step.State.init(allocator, unit_by_cell.len);
    defer surface_pond_transition_state.deinit();
    const surface_inorganic_nitrogen_export_g_n_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_inorganic_nitrogen_export_g_n_per_h);
    @memset(surface_inorganic_nitrogen_export_g_n_per_h, 0);
    const surface_inorganic_phosphorus_export_g_p_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_inorganic_phosphorus_export_g_p_per_h);
    @memset(surface_inorganic_phosphorus_export_g_p_per_h, 0);
    const surface_organic_carbon_export_g_c_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_organic_carbon_export_g_c_per_h);
    @memset(surface_organic_carbon_export_g_c_per_h, 0);
    const surface_inorganic_carbon_export_g_c_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_inorganic_carbon_export_g_c_per_h);
    @memset(surface_inorganic_carbon_export_g_c_per_h, 0);
    const surface_organic_nitrogen_export_g_n_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_organic_nitrogen_export_g_n_per_h);
    @memset(surface_organic_nitrogen_export_g_n_per_h, 0);
    const surface_organic_phosphorus_export_g_p_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_organic_phosphorus_export_g_p_per_h);
    @memset(surface_organic_phosphorus_export_g_p_per_h, 0);
    var surface_erosion_state = try ecosys.soil_erosion.RuntimeState.init(allocator, try runscript.domain.columns(), try runscript.domain.rows());
    defer surface_erosion_state.deinit();
    const surface_soil_mass_at_erosion_start_megagrams = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_soil_mass_at_erosion_start_megagrams);
    const net_sediment_megagrams_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(net_sediment_megagrams_per_h);
    @memset(surface_soil_mass_at_erosion_start_megagrams, 0);
    @memset(net_sediment_megagrams_per_h, 0);
    var eroded_mineral_state = try ecosys.soil_erosion_mineral_bridge.State.init(allocator, unit_by_cell.len);
    defer eroded_mineral_state.deinit();
    var soil_catalog = ecosys.soil_catalog.Catalog.init(allocator);
    defer soil_catalog.deinit();
    const catalog_index_by_cell = try allocator.alloc(usize, unit_by_cell.len);
    defer allocator.free(catalog_index_by_cell);
    for (unit_by_cell, 0..) |unit_index, cell_index| {
        const profile_name = topography.units[unit_index].soil_profile_file;
        const catalog_index = soil_catalog.find(profile_name) orelse load: {
            const profile_path = try resolveInputPath(allocator, init.io, runscript_directory, profile_name);
            defer allocator.free(profile_path);
            const profile_source = try std.Io.Dir.cwd().readFileAlloc(init.io, profile_path, allocator, .limited(64 * 1024 * 1024));
            defer allocator.free(profile_source);
            break :load try soil_catalog.appendFromSource(profile_name, profile_source, runscript.soil_solver_parameters.retention, runscript.soil_solver_parameters.profile_derivation);
        };
        catalog_index_by_cell[cell_index] = catalog_index;
    }

    const scene_options = try allocator.alloc(ecosys.options.SceneOptions, runscript.scenes.len);
    defer allocator.free(scene_options);
    const weather_forcing_geometry_by_cell = try allocator.alloc(
        ecosys.weather_grid.ForcingGeometry,
        grid_environment.site_catalog_index_by_cell.len,
    );
    defer allocator.free(weather_forcing_geometry_by_cell);
    for (weather_forcing_geometry_by_cell, 0..) |*geometry, cell| {
        const selected_site = site_catalog.entries.items[
            grid_environment.site_catalog_index_by_cell[cell]
        ].site;
        geometry.* = .{
            .altitude_m = selected_site.elevation_m,
            .latitude_degrees_north = selected_site.latitude_degrees_north,
            .longitude_degrees_east = selected_site.longitude_degrees_east,
            .phytotron = selected_site.ecosystem_type == -2,
        };
    }
    const scene_weather_assignments = try allocator.alloc(
        ecosys.weather_grid.Assignments,
        runscript.scenes.len,
    );
    defer allocator.free(scene_weather_assignments);
    var initialized_scene_weather_assignments: usize = 0;
    defer for (scene_weather_assignments[0..initialized_scene_weather_assignments]) |*assignments|
        assignments.deinit();

    var total_weather_records: usize = 0;
    for (runscript.scenes, 0..) |scene, scene_index| {
        const options_path = try resolveInputPath(allocator, init.io, runscript_directory, scene.options);
        defer allocator.free(options_path);
        const options_source = try std.Io.Dir.cwd().readFileAlloc(init.io, options_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(options_source);
        scene_options[scene_index] = ecosys.options.parse(options_source) catch |err| {
            std.log.err("scene {d} options validation failed for '{s}'", .{ scene_index + 1, scene.options });
            return err;
        };
        _ = try ecosys.iteration_control.Limits.fromSceneOptions(scene_options[scene_index]);
        const weather_grid_path = try resolveInputPath(
            allocator,
            init.io,
            runscript_directory,
            scene.weather_grid_file,
        );
        defer allocator.free(weather_grid_path);
        const weather_grid_source = try std.Io.Dir.cwd().readFileAlloc(
            init.io,
            weather_grid_path,
            allocator,
            .limited(16 * 1024 * 1024),
        );
        defer allocator.free(weather_grid_source);
        var weather_files = ecosys.grid_input_files.WeatherFiles.parse(
            allocator,
            weather_grid_source,
            try runscript.domain.columns(),
            try runscript.domain.rows(),
        ) catch |err| {
            std.log.err("scene {d} weather-grid validation failed for '{s}': {s}", .{
                scene_index + 1,
                scene.weather_grid_file,
                @errorName(err),
            });
            return err;
        };
        defer weather_files.deinit();
        scene_weather_assignments[scene_index] = try ecosys.weather_grid.Assignments.init(
            allocator,
            weather_files,
            weather_forcing_geometry_by_cell,
        );
        initialized_scene_weather_assignments += 1;
        for (scene_weather_assignments[scene_index].unique_file_names) |weather_name| {
            const weather_path = try resolveInputPath(allocator, init.io, runscript_directory, weather_name);
            defer allocator.free(weather_path);
            var weather_summary = ecosys.weather.scanFile(allocator, init.io, weather_path, 64 * 1024) catch |err| {
                std.log.err("scene {d} weather validation failed for '{s}'", .{ scene_index + 1, weather_name });
                return err;
            };
            defer weather_summary.deinit();
            total_weather_records = try std.math.add(
                usize,
                total_weather_records,
                weather_summary.observation_count,
            );
        }
    }
    const iteration_limits = try ecosys.iteration_control.Limits.fromSceneOptions(scene_options[0]);
    const timeline = try ecosys.simulation_timeline.summarize(scene_options, runscript.scenarios, runscript.execution_repeat_count);
    var plant_assignments: ?ecosys.plant_assignment.Assignments = null;
    defer if (plant_assignments) |*assignments| assignments.deinit();
    var plant_unit_by_cell: ?[]usize = null;
    defer if (plant_unit_by_cell) |map| allocator.free(map);
    var plant_catalog = ecosys.plant_catalog.Catalog.init(allocator);
    defer plant_catalog.deinit();
    var plant_management_catalog = ecosys.plant_management.Catalog.init(allocator);
    defer plant_management_catalog.deinit();
    var land_management_assignments: ?ecosys.land_management.Assignments = null;
    defer if (land_management_assignments) |*assignments| assignments.deinit();
    var land_management_unit_by_cell: ?[]usize = null;
    defer if (land_management_unit_by_cell) |map| allocator.free(map);
    var disturbance_catalog = ecosys.disturbance_schedule.Catalog.init(allocator);
    defer disturbance_catalog.deinit();
    const disturbance_schedule_maps = try allocator.alloc(?ecosys.disturbance_management_dispatch.ScheduleMap, runscript.scenes.len);
    @memset(disturbance_schedule_maps, null);
    defer {
        for (disturbance_schedule_maps) |*maybe_map| if (maybe_map.*) |*map| map.deinit();
        allocator.free(disturbance_schedule_maps);
    }
    var fertilizer_catalog = ecosys.fertilizer_schedule.Catalog.init(allocator);
    defer fertilizer_catalog.deinit();
    const fertilizer_schedule_maps = try allocator.alloc(?ecosys.fertilizer_management_dispatch.ScheduleMap, runscript.scenes.len);
    @memset(fertilizer_schedule_maps, null);
    defer {
        for (fertilizer_schedule_maps) |*maybe_map| if (maybe_map.*) |*map| map.deinit();
        allocator.free(fertilizer_schedule_maps);
    }
    var irrigation_catalog = ecosys.irrigation_schedule.Catalog.init(allocator);
    defer irrigation_catalog.deinit();
    const irrigation_schedule_maps = try allocator.alloc(?ecosys.irrigation_management_dispatch.ScheduleMap, runscript.scenes.len);
    @memset(irrigation_schedule_maps, null);
    defer {
        for (irrigation_schedule_maps) |*maybe_map| if (maybe_map.*) |*map| map.deinit();
        allocator.free(irrigation_schedule_maps);
    }
    var output_selection_catalog = ecosys.output_selection.Catalog.init(allocator);
    defer output_selection_catalog.deinit();
    for (runscript.scenes, 0..) |scene, scene_index| for (scene.output_editors) |editor_name| {
        if (ecosys.delimited_input.isNo(editor_name) or output_selection_catalog.find(editor_name) != null) continue;
        const path = try resolveInputPath(allocator, init.io, runscript_directory, editor_name);
        defer allocator.free(path);
        const editor_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(editor_source);
        _ = output_selection_catalog.appendFromSource(editor_name, editor_source) catch |err| {
            std.log.err("scene {d} output-selection validation failed for '{s}'", .{ scene_index + 1, editor_name });
            return err;
        };
    };
    if (!ecosys.delimited_input.isNo(runscript.scenes[0].land_management)) {
        const management_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.scenes[0].land_management);
        defer allocator.free(management_path);
        const management_source = try std.Io.Dir.cwd().readFileAlloc(init.io, management_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(management_source);
        land_management_assignments = try ecosys.land_management.parse(allocator, management_source);
        land_management_unit_by_cell = try land_management_assignments.?.buildCellUnitMap(allocator, runscript.domain);
        for (land_management_assignments.?.units) |unit| {
            if (!ecosys.delimited_input.isNo(unit.tillage_file) and disturbance_catalog.find(unit.tillage_file) == null) {
                const path = try resolveInputPath(allocator, init.io, runscript_directory, unit.tillage_file);
                defer allocator.free(path);
                const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(64 * 1024 * 1024));
                defer allocator.free(schedule_source);
                _ = try disturbance_catalog.appendFromSource(unit.tillage_file, schedule_source);
            }
            if (!ecosys.delimited_input.isNo(unit.fertilizer_file) and fertilizer_catalog.find(unit.fertilizer_file) == null) {
                const path = try resolveInputPath(allocator, init.io, runscript_directory, unit.fertilizer_file);
                defer allocator.free(path);
                const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(64 * 1024 * 1024));
                defer allocator.free(schedule_source);
                _ = try fertilizer_catalog.appendFromSource(unit.fertilizer_file, schedule_source);
            }
            if (!ecosys.delimited_input.isNo(unit.irrigation_file) and irrigation_catalog.find(unit.irrigation_file) == null) {
                const path = try resolveInputPath(allocator, init.io, runscript_directory, unit.irrigation_file);
                defer allocator.free(path);
                const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(64 * 1024 * 1024));
                defer allocator.free(schedule_source);
                _ = try irrigation_catalog.appendFromSource(unit.irrigation_file, schedule_source);
            }
        }
        fertilizer_schedule_maps[0] = try ecosys.fertilizer_management_dispatch.ScheduleMap.init(allocator, land_management_assignments.?, land_management_unit_by_cell.?, fertilizer_catalog);
        disturbance_schedule_maps[0] = try ecosys.disturbance_management_dispatch.ScheduleMap.init(allocator, land_management_assignments.?, land_management_unit_by_cell.?, disturbance_catalog);
        irrigation_schedule_maps[0] = try ecosys.irrigation_management_dispatch.ScheduleMap.init(allocator, land_management_assignments.?, land_management_unit_by_cell.?, irrigation_catalog);
    }
    if (!ecosys.delimited_input.isNo(runscript.scenes[0].plant_management)) {
        const plant_assignment_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.scenes[0].plant_management);
        defer allocator.free(plant_assignment_path);
        const plant_assignment_source = try std.Io.Dir.cwd().readFileAlloc(init.io, plant_assignment_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(plant_assignment_source);
        plant_assignments = try ecosys.plant_assignment.parse(allocator, plant_assignment_source, site.ecosystem_type);
        plant_unit_by_cell = try plant_assignments.?.buildCellUnitMap(allocator, runscript.domain, runscript.plant_species_count);
        for (plant_assignments.?.units) |unit| for (unit.species) |species| {
            if (plant_catalog.find(species.species_file) == null) {
                const traits_path = try resolveInputPath(allocator, init.io, runscript_directory, species.species_file);
                defer allocator.free(traits_path);
                const traits_source = try std.Io.Dir.cwd().readFileAlloc(init.io, traits_path, allocator, .limited(16 * 1024 * 1024));
                defer allocator.free(traits_source);
                _ = try plant_catalog.appendFromSource(species.species_file, traits_source);
            }
            if (plant_management_catalog.find(species.management_file) == null and !ecosys.delimited_input.isNo(species.management_file)) {
                const management_path = try resolveInputPath(allocator, init.io, runscript_directory, species.management_file);
                defer allocator.free(management_path);
                const management_source = try std.Io.Dir.cwd().readFileAlloc(init.io, management_path, allocator, .limited(16 * 1024 * 1024));
                defer allocator.free(management_source);
                _ = try plant_management_catalog.appendFromSource(species.management_file, management_source);
            }
        };
    }

    // Validate later scene assignments and populate the shared input catalogs.
    // The first scene maps remain alive for state initialization; subsequent
    // maps are checked and released immediately.
    for (runscript.scenes[1..], 1..) |scene, scene_index| {
        if (!ecosys.delimited_input.isNo(scene.land_management)) {
            const path = try resolveInputPath(allocator, init.io, runscript_directory, scene.land_management);
            defer allocator.free(path);
            const assignment_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(16 * 1024 * 1024));
            defer allocator.free(assignment_source);
            var assignments = ecosys.land_management.parse(allocator, assignment_source) catch |err| {
                std.log.err("scene {d} land-management validation failed for '{s}'", .{ scene_index + 1, scene.land_management });
                return err;
            };
            defer assignments.deinit();
            const map = try assignments.buildCellUnitMap(allocator, runscript.domain);
            defer allocator.free(map);
            for (assignments.units) |unit| {
                if (!ecosys.delimited_input.isNo(unit.tillage_file) and disturbance_catalog.find(unit.tillage_file) == null) {
                    const schedule_path = try resolveInputPath(allocator, init.io, runscript_directory, unit.tillage_file);
                    defer allocator.free(schedule_path);
                    const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, schedule_path, allocator, .limited(64 * 1024 * 1024));
                    defer allocator.free(schedule_source);
                    _ = try disturbance_catalog.appendFromSource(unit.tillage_file, schedule_source);
                }
                if (!ecosys.delimited_input.isNo(unit.fertilizer_file) and fertilizer_catalog.find(unit.fertilizer_file) == null) {
                    const schedule_path = try resolveInputPath(allocator, init.io, runscript_directory, unit.fertilizer_file);
                    defer allocator.free(schedule_path);
                    const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, schedule_path, allocator, .limited(64 * 1024 * 1024));
                    defer allocator.free(schedule_source);
                    _ = try fertilizer_catalog.appendFromSource(unit.fertilizer_file, schedule_source);
                }
                if (!ecosys.delimited_input.isNo(unit.irrigation_file) and irrigation_catalog.find(unit.irrigation_file) == null) {
                    const schedule_path = try resolveInputPath(allocator, init.io, runscript_directory, unit.irrigation_file);
                    defer allocator.free(schedule_path);
                    const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, schedule_path, allocator, .limited(64 * 1024 * 1024));
                    defer allocator.free(schedule_source);
                    _ = try irrigation_catalog.appendFromSource(unit.irrigation_file, schedule_source);
                }
            }
            fertilizer_schedule_maps[scene_index] = try ecosys.fertilizer_management_dispatch.ScheduleMap.init(allocator, assignments, map, fertilizer_catalog);
            disturbance_schedule_maps[scene_index] = try ecosys.disturbance_management_dispatch.ScheduleMap.init(allocator, assignments, map, disturbance_catalog);
            irrigation_schedule_maps[scene_index] = try ecosys.irrigation_management_dispatch.ScheduleMap.init(allocator, assignments, map, irrigation_catalog);
        }
        if (!ecosys.delimited_input.isNo(scene.plant_management)) {
            const path = try resolveInputPath(allocator, init.io, runscript_directory, scene.plant_management);
            defer allocator.free(path);
            const assignment_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(16 * 1024 * 1024));
            defer allocator.free(assignment_source);
            var assignments = ecosys.plant_assignment.parse(allocator, assignment_source, site.ecosystem_type) catch |err| {
                std.log.err("scene {d} plant-management validation failed for '{s}'", .{ scene_index + 1, scene.plant_management });
                return err;
            };
            defer assignments.deinit();
            const map = try assignments.buildCellUnitMap(allocator, runscript.domain, runscript.plant_species_count);
            defer allocator.free(map);
            for (assignments.units) |unit| for (unit.species) |species| {
                if (plant_catalog.find(species.species_file) == null) {
                    const traits_path = try resolveInputPath(allocator, init.io, runscript_directory, species.species_file);
                    defer allocator.free(traits_path);
                    const traits_source = try std.Io.Dir.cwd().readFileAlloc(init.io, traits_path, allocator, .limited(16 * 1024 * 1024));
                    defer allocator.free(traits_source);
                    _ = try plant_catalog.appendFromSource(species.species_file, traits_source);
                }
                if (!ecosys.delimited_input.isNo(species.management_file) and plant_management_catalog.find(species.management_file) == null) {
                    const management_path = try resolveInputPath(allocator, init.io, runscript_directory, species.management_file);
                    defer allocator.free(management_path);
                    const management_source = try std.Io.Dir.cwd().readFileAlloc(init.io, management_path, allocator, .limited(16 * 1024 * 1024));
                    defer allocator.free(management_source);
                    _ = try plant_management_catalog.appendFromSource(species.management_file, management_source);
                }
            };
        }
    }

    const config = try ecosys.config.SimulationConfig.init(
        .{
            .lon_count = try runscript.domain.columns(),
            .lat_count = try runscript.domain.rows(),
            .soil_layers = try soil_catalog.maximumLayerCount(),
            .plant_populations = runscript.plant_species_count,
        },
        .{ .worker_threads = runscript.worker_count, .tile_cells = runscript.tile_cell_count },
        .{
            .relative_tolerance = runscript.relative_tolerance,
            .absolute_tolerance = runscript.absolute_tolerance,
            .mass_balance_tolerance = runscript.mass_balance_tolerance,
            .negligible_quantity_threshold = runscript.negligible_quantity_threshold,
            .max_nonlinear_iterations = iteration_limits.water_heat_solute_max_iterations,
            .picard_relaxation = runscript.picard_relaxation,
        },
    );
    const resolved_output_editors = try allocator.alloc(?ecosys.output_editor_layout.Resolved, try std.math.mul(usize, runscript.scenes.len, 10));
    @memset(resolved_output_editors, null);
    defer {
        for (resolved_output_editors) |*maybe_editor| if (maybe_editor.*) |*editor| editor.deinit();
        allocator.free(resolved_output_editors);
    }
    for (runscript.scenes, 0..) |scene, scene_index| for (scene.output_editors, 0..) |editor_name, editor_index| {
        if (ecosys.delimited_input.isNo(editor_name)) continue;
        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
        var soil_catalog_for_editor = try soilOutputCatalog(allocator, editor_index, config.soil_layers);
        defer soil_catalog_for_editor.deinit();
        var plant_catalog_for_editor = try plantOutputCatalog(allocator, editor_index, config.soil_layers);
        defer plant_catalog_for_editor.deinit();
        resolved_output_editors[scene_index * 10 + editor_index] = try ecosys.output_editor_layout.resolve(
            allocator,
            output_selection_catalog.entries.items[selection_index].selection.enabled_variables,
            soil_catalog_for_editor.variables.len,
            plant_catalog_for_editor.variables.len,
        );
    };
    var state = try ecosys.grid.GridState.init(allocator, config);
    defer state.deinit();
    var any_visualization_enabled = false;
    for (scene_options) |options|
        any_visualization_enabled =
            any_visualization_enabled or options.visualization_enabled;
    var visualization_streams = try ecosys.visualization_output.Streams.init(
        allocator,
        init.io,
        std.Io.Dir.cwd(),
        any_visualization_enabled,
        state.cell_count,
        config.plant_populations,
        config.soil_layers,
        scene_options[0].visualization_start_year,
        scene_options[0].visualization_end_year,
        64 * 1024,
    );
    defer visualization_streams.deinit();
    const visualization_species_respired_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_species_respired_carbon_g);
    const visualization_leaf_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_leaf_carbon_g);
    const visualization_stalk_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_stalk_carbon_g);
    const visualization_root_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_root_carbon_g);
    const visualization_root_profile_carbon_g =
        try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(visualization_root_profile_carbon_g);
    const visualization_organic_carbon_g =
        try allocator.alloc(f64, config.soil_layers + 1);
    defer allocator.free(visualization_organic_carbon_g);
    var daily_soil_gas_flux = try ecosys.soil_daily_gas_flux.State.init(allocator, state.cell_count);
    defer daily_soil_gas_flux.deinit();
    // HOUR1 lines 2408-2432: per-cell hourly accumulator array (XCNET, XHNET, XONET, etc.).
    const hourly_accumulators = try allocator.alloc(
        ecosys.hourly_grid_cell_accumulator_reset.GridCellAccumulators,
        state.cell_count,
    );
    defer allocator.free(hourly_accumulators);
    for (hourly_accumulators) |*acc| ecosys.hourly_grid_cell_accumulator_reset.reset(acc);
    var daily_canopy_gas_exchange = try ecosys.daily_canopy_gas_exchange.State.init(allocator, state.cell_count);
    defer daily_canopy_gas_exchange.deinit();
    var daily_heterotrophic_respiration = try ecosys.soil_daily_heterotrophic_respiration.State.init(allocator, state.cell_count, state.layer_count);
    defer daily_heterotrophic_respiration.deinit();
    var daily_carbon_export = try ecosys.soil_daily_carbon_export.State.init(allocator, state.cell_count);
    defer daily_carbon_export.deinit();
    var daily_phosphorus_export = try ecosys.soil_daily_phosphorus_export.State.init(allocator, state.cell_count);
    defer daily_phosphorus_export.deinit();
    var daily_nitrogen_export = try ecosys.soil_daily_nitrogen_export.State.init(allocator, state.cell_count);
    defer daily_nitrogen_export.deinit();
    var daily_water_ledger = try ecosys.soil_daily_water_ledger.State.init(allocator, state.cell_count);
    defer daily_water_ledger.deinit();
    var daily_heat_ledger = try ecosys.soil_daily_heat_ledger.State.init(allocator, state.cell_count, config.soil_layers);
    defer daily_heat_ledger.deinit();
    var landscape_mass_balance_state: ecosys.landscape_mass_balance_checkpoint.State = .{
        .boundary_ledger = .{},
        .monitor = null,
    };
    const landscape_soil_mass_megagrams_scratch = try allocator.alloc(
        f64,
        state.layer_count,
    );
    defer allocator.free(landscape_soil_mass_megagrams_scratch);
    @memset(landscape_soil_mass_megagrams_scratch, 0);
    var daily_ecosystem_carbon = try ecosys.soil_daily_ecosystem_carbon.State.init(allocator, state.cell_count);
    defer daily_ecosystem_carbon.deinit();
    const daily_soil_fire_charcoal_production_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_charcoal_production_g_c);
    @memset(daily_soil_fire_charcoal_production_g_c, 0);
    const daily_soil_fire_carbon_dioxide_emission_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_carbon_dioxide_emission_g_c);
    @memset(daily_soil_fire_carbon_dioxide_emission_g_c, 0);
    const daily_soil_fire_methane_emission_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_methane_emission_g_c);
    @memset(daily_soil_fire_methane_emission_g_c, 0);
    const daily_soil_fire_phosphorus_flux_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_phosphorus_flux_g_p);
    @memset(daily_soil_fire_phosphorus_flux_g_p, 0);
    const daily_soil_fire_nitrogen_flux_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_nitrogen_flux_g_n);
    @memset(daily_soil_fire_nitrogen_flux_g_n, 0);
    const daily_root_fire_carbon_dioxide_emission_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_root_fire_carbon_dioxide_emission_g_c);
    @memset(daily_root_fire_carbon_dioxide_emission_g_c, 0);
    const daily_root_fire_methane_emission_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_root_fire_methane_emission_g_c);
    @memset(daily_root_fire_methane_emission_g_c, 0);
    const daily_root_fire_nitrogen_flux_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_root_fire_nitrogen_flux_g_n);
    @memset(daily_root_fire_nitrogen_flux_g_n, 0);
    const daily_root_fire_phosphorus_flux_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_root_fire_phosphorus_flux_g_p);
    @memset(daily_root_fire_phosphorus_flux_g_p, 0);
    const daily_root_fire_nitrogen_loss_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_root_fire_nitrogen_loss_g_n);
    @memset(daily_root_fire_nitrogen_loss_g_n, 0);
    const daily_root_fire_phosphorus_loss_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_root_fire_phosphorus_loss_g_p);
    @memset(daily_root_fire_phosphorus_loss_g_p, 0);
    const daily_soil_combusted_nitrogen_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_combusted_nitrogen_g_n);
    @memset(daily_soil_combusted_nitrogen_g_n, 0);
    const daily_soil_combusted_phosphorus_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_combusted_phosphorus_g_p);
    @memset(daily_soil_combusted_phosphorus_g_p, 0);
    const daily_microbial_phosphate_mineralization_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_microbial_phosphate_mineralization_g_p);
    @memset(daily_microbial_phosphate_mineralization_g_p, 0);
    const daily_microbial_nitrogen_mineralization_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_microbial_nitrogen_mineralization_g_n);
    @memset(daily_microbial_nitrogen_mineralization_g_n, 0);
    const daily_organic_fertilizer_carbon_input_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_organic_fertilizer_carbon_input_g_c);
    @memset(daily_organic_fertilizer_carbon_input_g_c, 0);
    const daily_biome_organic_carbon_input_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_biome_organic_carbon_input_g_c);
    @memset(daily_biome_organic_carbon_input_g_c, 0);
    const daily_organic_fertilizer_phosphorus_input_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_organic_fertilizer_phosphorus_input_g_p);
    @memset(daily_organic_fertilizer_phosphorus_input_g_p, 0);
    const daily_organic_fertilizer_nitrogen_input_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_organic_fertilizer_nitrogen_input_g_n);
    @memset(daily_organic_fertilizer_nitrogen_input_g_n, 0);
    var plant_state = try ecosys.grid.PlantState.init(allocator, config);
    defer plant_state.deinit();
    var hourly_soil_carbon_catalog = try ecosys.soil_output_catalog.carbon(allocator, config.soil_layers, config.soil_layers, config.soil_layers);
    defer hourly_soil_carbon_catalog.deinit();
    var hourly_soil_carbon_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_carbon_enabled = hourly_soil_carbon_enabled or !ecosys.delimited_input.isNo(scene.output_editors[0]);
    var hourly_soil_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_carbon_enabled, state.cell_count, hourly_soil_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_carbon_bank.deinit();
    var hourly_plant_carbon_catalog = try ecosys.plant_output_catalog.carbon(allocator);
    defer hourly_plant_carbon_catalog.deinit();
    var hourly_plant_carbon_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[0])) continue;
        const resolved = resolved_output_editors[scene_index * 10] orelse continue;
        hourly_plant_carbon_enabled = hourly_plant_carbon_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    const output_plant_count = try std.math.mul(usize, state.cell_count, config.plant_populations);
    var hourly_plant_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_carbon_enabled, output_plant_count, hourly_plant_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_carbon_bank.deinit();
    var hourly_soil_water_catalog = try ecosys.soil_output_catalog.water(allocator, config.soil_layers);
    defer hourly_soil_water_catalog.deinit();
    var hourly_soil_water_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_water_enabled = hourly_soil_water_enabled or !ecosys.delimited_input.isNo(scene.output_editors[1]);
    var hourly_soil_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_water_enabled, state.cell_count, hourly_soil_water_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_water_bank.deinit();
    var hourly_plant_water_catalog = try ecosys.plant_output_catalog.water(allocator, config.soil_layers);
    defer hourly_plant_water_catalog.deinit();
    var hourly_plant_water_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[1])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 1] orelse continue;
        hourly_plant_water_enabled = hourly_plant_water_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_water_enabled, output_plant_count, hourly_plant_water_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_water_bank.deinit();
    var daily_plant_water_catalog = try ecosys.plant_output_catalog.dailyWater(allocator);
    defer daily_plant_water_catalog.deinit();
    var daily_plant_water_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[6])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 6] orelse continue;
        daily_plant_water_enabled = daily_plant_water_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_water_enabled, output_plant_count, daily_plant_water_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_water_bank.deinit();
    var daily_plant_carbon_catalog = try ecosys.plant_output_catalog.dailyCarbon(allocator, config.soil_layers);
    defer daily_plant_carbon_catalog.deinit();
    var daily_plant_carbon_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[5])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 5] orelse continue;
        daily_plant_carbon_enabled = daily_plant_carbon_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_carbon_enabled, output_plant_count, daily_plant_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_carbon_bank.deinit();
    var daily_soil_carbon_catalog = try ecosys.soil_output_catalog.dailyCarbon(allocator, config.soil_layers);
    defer daily_soil_carbon_catalog.deinit();
    var daily_soil_carbon_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[5])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 5] orelse continue;
        daily_soil_carbon_enabled = daily_soil_carbon_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_carbon_enabled, state.cell_count, daily_soil_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_carbon_bank.deinit();
    const daily_soil_water_potential_layers = dailySoilWaterPotentialLayerCount(config.soil_layers);
    var daily_soil_water_catalog = try ecosys.soil_output_catalog.dailyWater(allocator, config.soil_layers, config.soil_layers, daily_soil_water_potential_layers);
    defer daily_soil_water_catalog.deinit();
    var daily_soil_water_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[6])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 6] orelse continue;
        daily_soil_water_enabled = daily_soil_water_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_water_enabled, state.cell_count, daily_soil_water_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_water_bank.deinit();
    const daily_water_profile_workspace = try allocator.alloc(f64, 3 * config.soil_layers);
    defer allocator.free(daily_water_profile_workspace);
    var daily_soil_phosphorus_catalog = try ecosys.soil_output_catalog.dailyPhosphorus(allocator, config.soil_layers);
    defer daily_soil_phosphorus_catalog.deinit();
    var daily_soil_phosphorus_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[8])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 8] orelse continue;
        daily_soil_phosphorus_enabled = daily_soil_phosphorus_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_phosphorus_enabled, state.cell_count, daily_soil_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_phosphorus_bank.deinit();
    var daily_soil_nitrogen_catalog = try ecosys.soil_output_catalog.dailyNitrogen(allocator, config.soil_layers);
    defer daily_soil_nitrogen_catalog.deinit();
    var daily_soil_nitrogen_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[7])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 7] orelse continue;
        daily_soil_nitrogen_enabled = daily_soil_nitrogen_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_nitrogen_enabled, state.cell_count, daily_soil_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_nitrogen_bank.deinit();
    var daily_soil_heat_catalog = try ecosys.soil_output_catalog.dailyHeat(allocator, config.soil_layers, config.soil_layers);
    defer daily_soil_heat_catalog.deinit();
    var daily_soil_heat_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[9])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 9] orelse continue;
        daily_soil_heat_enabled = daily_soil_heat_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_heat_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_heat_enabled, state.cell_count, daily_soil_heat_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_heat_bank.deinit();
    const daily_soil_heat_workspace = try allocator.alloc(f64, 3 * config.soil_layers);
    defer allocator.free(daily_soil_heat_workspace);
    const daily_soil_organic_carbon_by_layer_g_c = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_organic_carbon_by_layer_g_c);
    @memset(daily_soil_organic_carbon_by_layer_g_c, 0);
    const daily_soil_carbonate_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_carbonate_mol_per_m3);
    const daily_soil_bicarbonate_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_bicarbonate_mol_per_m3);
    const daily_soil_carbonate_complexes_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_carbonate_complexes_mol_per_m3);
    const daily_soil_calcite_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_calcite_mol_per_m3);
    var daily_plant_nitrogen_catalog = try ecosys.plant_output_catalog.dailyNitrogen(allocator);
    defer daily_plant_nitrogen_catalog.deinit();
    var daily_plant_nitrogen_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[7])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 7] orelse continue;
        daily_plant_nitrogen_enabled = daily_plant_nitrogen_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_nitrogen_enabled, output_plant_count, daily_plant_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_nitrogen_bank.deinit();
    var daily_plant_phosphorus_catalog = try ecosys.plant_output_catalog.dailyPhosphorus(allocator);
    defer daily_plant_phosphorus_catalog.deinit();
    var daily_plant_phosphorus_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[8])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 8] orelse continue;
        daily_plant_phosphorus_enabled = daily_plant_phosphorus_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_phosphorus_enabled, output_plant_count, daily_plant_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_phosphorus_bank.deinit();
    const daily_element_root_workspace = try allocator.alloc(f64, try std.math.mul(usize, output_plant_count, config.soil_layers));
    defer allocator.free(daily_element_root_workspace);
    @memset(daily_element_root_workspace, 0);
    const zero_root_layer_flux_g_per_h = try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(zero_root_layer_flux_g_per_h);
    @memset(zero_root_layer_flux_g_per_h, 0);
    var daily_plant_development_catalog = try ecosys.plant_output_catalog.dailyDevelopment(allocator);
    defer daily_plant_development_catalog.deinit();
    var daily_plant_development_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[9])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 9] orelse continue;
        daily_plant_development_enabled = daily_plant_development_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_development_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_development_enabled, output_plant_count, daily_plant_development_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_development_bank.deinit();
    var hourly_soil_nitrogen_catalog = try ecosys.soil_output_catalog.nitrogen(allocator, config.soil_layers, config.soil_layers);
    defer hourly_soil_nitrogen_catalog.deinit();
    var hourly_soil_nitrogen_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_nitrogen_enabled = hourly_soil_nitrogen_enabled or !ecosys.delimited_input.isNo(scene.output_editors[2]);
    var hourly_soil_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_nitrogen_enabled, state.cell_count, hourly_soil_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_nitrogen_bank.deinit();
    var hourly_plant_nitrogen_catalog = try ecosys.plant_output_catalog.nitrogen(allocator, config.soil_layers);
    defer hourly_plant_nitrogen_catalog.deinit();
    var hourly_plant_nitrogen_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[2])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 2] orelse continue;
        hourly_plant_nitrogen_enabled = hourly_plant_nitrogen_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_nitrogen_enabled, output_plant_count, hourly_plant_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_nitrogen_bank.deinit();
    var hourly_soil_phosphorus_catalog = try ecosys.soil_output_catalog.phosphorus(allocator);
    defer hourly_soil_phosphorus_catalog.deinit();
    var hourly_soil_phosphorus_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_phosphorus_enabled = hourly_soil_phosphorus_enabled or !ecosys.delimited_input.isNo(scene.output_editors[3]);
    var hourly_soil_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_phosphorus_enabled, state.cell_count, hourly_soil_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_phosphorus_bank.deinit();
    var hourly_plant_phosphorus_catalog = try ecosys.plant_output_catalog.phosphorus(allocator, config.soil_layers);
    defer hourly_plant_phosphorus_catalog.deinit();
    var hourly_plant_phosphorus_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[3])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 3] orelse continue;
        hourly_plant_phosphorus_enabled = hourly_plant_phosphorus_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_phosphorus_enabled, output_plant_count, hourly_plant_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_phosphorus_bank.deinit();
    var hourly_soil_heat_catalog = try ecosys.soil_output_catalog.heat(allocator, config.soil_layers);
    defer hourly_soil_heat_catalog.deinit();
    var hourly_soil_heat_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_heat_enabled = hourly_soil_heat_enabled or !ecosys.delimited_input.isNo(scene.output_editors[4]);
    var hourly_soil_heat_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_heat_enabled, state.cell_count, hourly_soil_heat_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_heat_bank.deinit();
    var hourly_plant_heat_catalog = try ecosys.plant_output_catalog.heat(allocator);
    defer hourly_plant_heat_catalog.deinit();
    var hourly_plant_heat_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[4])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 4] orelse continue;
        hourly_plant_heat_enabled = hourly_plant_heat_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_heat_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_heat_enabled, output_plant_count, hourly_plant_heat_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_heat_bank.deinit();
    try state.validateFinite();
    try plant_state.validateFinite();
    const executor = try ecosys.compute.CpuExecutor.init(allocator, config.worker_threads, config.tile_cells);
    var soil_thermal_state = try ecosys.soil_thermal.State.initMapped(
        allocator,
        state,
        soil_catalog.entries.items,
        catalog_index_by_cell,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
    );
    defer soil_thermal_state.deinit();
    var initialization_context = ecosys.model_initialization.MappedHydrologyContext{
        .grid = &state,
        .catalog_entries = soil_catalog.entries.items,
        .catalog_index_by_cell = catalog_index_by_cell,
        .horizontal_cell_width_m = grid_environment.horizontal_cell_width_m,
        .vertical_cell_width_m = grid_environment.vertical_cell_width_m,
    };
    try runKernelAcrossSerialTilePlan(
        executor,
        &tile_plan.?,
        &initialization_context,
        ecosys.model_initialization.initializeMappedHydrologyTile,
    );
    const site_by_cell = try allocator.alloc(ecosys.site.Site, state.cell_count);
    defer allocator.free(site_by_cell);
    for (site_by_cell, grid_environment.site_catalog_index_by_cell) |*selected_site, site_index|
        selected_site.* = site_catalog.entries.items[site_index].site;
    const initial_surface_elevation_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_surface_elevation_m);
    for (site_by_cell, initial_surface_elevation_m) |selected_site, *elevation_m|
        elevation_m.* = selected_site.elevation_m;
    try terrain_hydrology_state.bindInitialSurfaceElevations(initial_surface_elevation_m);
    const surface_runoff_boundary_storage = try allocator.alloc(
        f64,
        try std.math.mul(usize, state.cell_count, 4),
    );
    defer allocator.free(surface_runoff_boundary_storage);
    const mean_annual_temperature_c_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(mean_annual_temperature_c_by_cell);
    const mean_annual_temperature_k_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(mean_annual_temperature_k_by_cell);
    const geothermal_enabled_by_cell = try allocator.alloc(bool, state.cell_count);
    defer allocator.free(geothermal_enabled_by_cell);
    const surface_runoff_boundary_fraction_by_direction = [4][]f64{
        surface_runoff_boundary_storage[0 * state.cell_count ..][0..state.cell_count],
        surface_runoff_boundary_storage[1 * state.cell_count ..][0..state.cell_count],
        surface_runoff_boundary_storage[2 * state.cell_count ..][0..state.cell_count],
        surface_runoff_boundary_storage[3 * state.cell_count ..][0..state.cell_count],
    };
    for (site_by_cell, 0..) |selected_site, cell| {
        mean_annual_temperature_c_by_cell[cell] =
            selected_site.mean_annual_air_temperature_c;
        mean_annual_temperature_k_by_cell[cell] =
            selected_site.mean_annual_air_temperature_c + 273.15;
        geothermal_enabled_by_cell[cell] =
            runscript.geothermal_controls.enabled and selected_site.ecosystem_type != -2;
        for (0..4) |direction| {
            surface_runoff_boundary_fraction_by_direction[direction][cell] =
                selected_site.surface_runoff_boundary_fraction[direction];
        }
        // STARTS lines 397--400 and 1222: `TKS(L,NY,NX)=ATKS(NY,NX)` seeds the
        // profile isothermally at the mean annual soil temperature. See
        // `seed_soil_profile_from_mean_annual` at the top of this file for the
        // measured effect on the EXEC heat audit.
        if (seed_soil_profile_from_mean_annual) {
            const initial_soil_temperature_k = selected_site.mean_annual_air_temperature_c + 273.15;
            if (!std.math.isFinite(initial_soil_temperature_k) or initial_soil_temperature_k <= 0)
                return error.InvalidInitialSoilTemperature;
            const layer_base = cell * state.soil_layer_capacity;
            for (0..state.soil_layer_capacity) |layer|
                state.soil_temperature_k[layer_base + layer] = initial_soil_temperature_k;
            // BIND-STARTS-654, source `starts.f` 659: `TKS(0,NY,NX)=ATKS(NY,NX)`.
            // This is the LAYER-0 companion of the `starts.f` 1222 seed
            // immediately above, and it was left on the `273.15 K` placeholder
            // that `grid.zig:79` installs at allocation. `273.15 K` is the
            // pure-water melting point, so every example began with its surface
            // litter sitting exactly on the ice/liquid phase boundary
            // regardless of site climate, which is the worst possible start for
            // a phase-enthalpy solver for the same reason recorded at the
            // `seed_soil_profile_from_mean_annual` declaration.
            //
            // A8 measured the placeholder error across the whole migrated
            // corpus and it FLIPS SIGN, so it cannot be corrected by an offset
            // and Ottawa is not the worst case: Ottawa `-5.40 K`, Boreal Crop
            // AB `-2.00 K`, Arctic Fen CH `+6.90 K`, Arctic Tundra IQ
            // `+9.80 K`. On the warm sites the surface started too cold; on the
            // frozen sites it held a `-9.8 C` site at the melting point, the
            // wrong side of the phase boundary for a permafrost column.
            //
            // Bound here rather than in `grid.zig` deliberately. The value is
            // `lower_soil_heat_boundary_initialization.initialize`'s
            // `surface_litter_temperature_k`, which that module sets to exactly
            // `inputs.mean_annual_soil_temperature_k` (`:108--109`), the same
            // per-cell quantity already in hand here. Seeding it at the one
            // place that already owns the mean-annual seed keeps a single
            // writer, needs no claim on the shared `src/grid.zig`, and leaves
            // `grid.zig:79`'s `@memset` as the correct allocation-time default
            // for callers that never reach site selection. Per A8's
            // recommendation this binds ONLY the surface temperature and does
            // not wire `lower_heat_source_depth_m` or `deep_source_temperature_k`,
            // which `solver.zig:498--504` and `:698--699` already
            // recompute per step from the same runscript controls, so `INIT-003`
            // is deferred entirely rather than half-answered.
            //
            // No double mutation: `surface_temperature_k` is written here before
            // any solver runs, and the hourly surface-temperature solver is its
            // only other writer, which advances this value rather than
            // establishing it. This runs before the EXEC heat baseline is
            // captured, which is the whole point: the baseline must see the
            // physical initial condition, not the placeholder.
            state.surface_temperature_k[cell] = initial_soil_temperature_k;
        }
    }
    var soil_boundary_topology_state = try ecosys.soil_boundary_topology.State.initMapped(
        allocator,
        &state,
        &terrain_hydrology_state,
        config.lon_count,
        config.lat_count,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
        site_by_cell,
    );
    defer soil_boundary_topology_state.deinit();
    var soil_solver_property_state = try ecosys.soil_solver_properties.State.initMapped(
        allocator,
        &state,
        soil_catalog.entries.items,
        catalog_index_by_cell,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
        runscript.soil_solver_parameters,
    );
    defer soil_solver_property_state.deinit();
    var soil_geometry_state = try ecosys.soil_layer_geometry.State.init(allocator, state.cell_count, config.soil_layers);
    defer soil_geometry_state.deinit();
    var surface_pond_domain_workspace = try ecosys.surface_pond_domain_transaction.Workspace.init(allocator, state.cell_count, config.soil_layers);
    defer surface_pond_domain_workspace.deinit();
    var soil_profile_relayering_workspace = try ecosys.soil_profile_relayering.Workspace.init(allocator, state.cell_count, config.soil_layers);
    defer soil_profile_relayering_workspace.deinit();
    for (0..state.cell_count) |cell| {
        const first_layer: usize = 0;
        const active_count = state.active_soil_layer_count[cell];
        const layer_base = cell * config.soil_layers;
        try ecosys.soil_layer_geometry.initializeCell(&soil_geometry_state, cell, first_layer, soil_solver_property_state.layer_thickness_m[layer_base .. layer_base + active_count], 0, runscript.soil_geometry_parameters.minimum_layer_thickness_m);
    }
    const fertilizer_band_layer_upper_depth_m =
        try allocator.alloc(f64, state.layer_count);
    defer allocator.free(fertilizer_band_layer_upper_depth_m);
    const fertilizer_band_layer_lower_depth_m =
        try allocator.alloc(f64, state.layer_count);
    defer allocator.free(fertilizer_band_layer_lower_depth_m);
    for (0..state.cell_count) |cell| {
        const boundary_base = cell * (config.soil_layers + 1);
        const layer_base = cell * config.soil_layers;
        const surface_depth_m =
            soil_geometry_state.boundary_depth_m[boundary_base];
        for (0..config.soil_layers) |layer| {
            fertilizer_band_layer_upper_depth_m[layer_base + layer] =
                soil_geometry_state.boundary_depth_m[
                    boundary_base + layer
                ] - surface_depth_m;
            fertilizer_band_layer_lower_depth_m[layer_base + layer] =
                soil_geometry_state.boundary_depth_m[
                    boundary_base + layer + 1
                ] - surface_depth_m;
        }
    }
    var fertilizer_band_state =
        try ecosys.fertilizer_band_state.State.initFromPlantNutrientParameters(
            allocator,
            state.cell_count,
            config.soil_layers,
            soil_geometry_state.active_layer_count,
            fertilizer_band_layer_upper_depth_m,
            fertilizer_band_layer_lower_depth_m,
            soil_geometry_state.layer_thickness_m,
            runscript.plant_nutrient_initialization,
        );
    defer fertilizer_band_state.deinit();
    var soil_microbial_state = try ecosys.soil_microbial_state.State.init(allocator, state.cell_count, config.soil_layers, runscript.microbial_substrate_count, runscript.microbial_population_count);
    defer soil_microbial_state.deinit();
    const soil_nitrogen_process_units_per_layer = try std.math.mul(usize, runscript.microbial_substrate_count, runscript.microbial_population_count);
    var soil_reactive_nitrogen_state = try ecosys.soil_reactive_nitrogen_state.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_reactive_nitrogen_state.deinit();
    var soil_microbial_phosphorus_state = try ecosys.soil_microbial_phosphorus_state.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_microbial_phosphorus_state.deinit();
    var soil_microbial_turnover_state = try ecosys.soil_microbial_turnover_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_microbial_turnover_state.deinit();
    var soil_litter_colonization_state = try ecosys.soil_litter_colonization_step.State.init(allocator, state.layer_count);
    defer soil_litter_colonization_state.deinit();
    var soil_organic_sorption_state = try ecosys.soil_organic_sorption_step.State.init(allocator, state.layer_count);
    defer soil_organic_sorption_state.deinit();
    var soil_organic_decomposition_state = try ecosys.soil_organic_decomposition_step.State.init(allocator, state.layer_count);
    defer soil_organic_decomposition_state.deinit();
    var soil_organic_priming_state = try ecosys.soil_organic_priming_step.State.init(allocator, state.layer_count, runscript.microbial_population_count);
    defer soil_organic_priming_state.deinit();
    var soil_respiration_products_state = try ecosys.soil_respiration_products_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_respiration_products_state.deinit();
    var soil_autotrophic_carbon_state = try ecosys.soil_autotrophic_carbon_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_autotrophic_carbon_state.deinit();
    var soil_microbial_layer_mixing_state = try ecosys.soil_microbial_layer_mixing.State.init(allocator, state.layer_count);
    defer soil_microbial_layer_mixing_state.deinit();
    var soil_biogeochemical_gas_flux_state = try ecosys.soil_biogeochemical_gas_aggregation.State.init(allocator, state.layer_count);
    defer soil_biogeochemical_gas_flux_state.deinit();
    var soil_methane_state = try ecosys.soil_methane_step.State.init(allocator, state.layer_count);
    defer soil_methane_state.deinit();
    var organic_matter_fire_exchange_state = try ecosys.organic_matter_fire_exchange.State.init(allocator, state.layer_count, ecosys.soil_organic_initialization.microbial_substrate_count);
    defer organic_matter_fire_exchange_state.deinit();
    var surface_fire_exchange_state = try ecosys.organic_matter_fire_exchange.State.init(allocator, state.cell_count, ecosys.soil_organic_initialization.microbial_substrate_count);
    defer surface_fire_exchange_state.deinit();
    var soil_microbial_oxygen_state = try ecosys.soil_oxygen_allocation.State.init(allocator, state.cell_count, config.soil_layers, soil_nitrogen_process_units_per_layer);
    defer soil_microbial_oxygen_state.deinit();
    var soil_oxygen_staging_state = try ecosys.soil_oxygen_step.State.init(allocator, state.cell_count, config.soil_layers, soil_nitrogen_process_units_per_layer);
    defer soil_oxygen_staging_state.deinit();
    var soil_nitrogen_flux_workspace = try ecosys.soil_nitrogen_flux_workspace.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_nitrogen_flux_workspace.deinit();
    var soil_nitrifier_environment_state = try ecosys.soil_nitrifier_environment_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_nitrifier_environment_state.deinit();
    const soil_field_capacity_fraction = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_field_capacity_fraction);
    for (soil_solver_property_state.retention_curve, soil_field_capacity_fraction) |curve, *field_capacity| field_capacity.* = curve.curve.field_capacity_fraction;
    const soil_redox_satisfaction_fraction = try allocator.alloc(f64, try std.math.mul(usize, state.layer_count, soil_nitrogen_process_units_per_layer));
    defer allocator.free(soil_redox_satisfaction_fraction);
    @memset(soil_redox_satisfaction_fraction, 1);
    // Runtime transport domains. Capacities derive only from parsed grid and
    // profile inputs; there is no PARAMETERS.H-equivalent extent.
    const aqueous_species_count = ecosys.solute_transport_species.AqueousSpecies.count;
    var micropore_solute_state = try ecosys.solute_transport.State.init(allocator, state.layer_count, aqueous_species_count);
    defer micropore_solute_state.deinit();
    var macropore_solute_state = try ecosys.solute_transport.State.init(allocator, state.layer_count, aqueous_species_count);
    defer macropore_solute_state.deinit();
    const soil_recharge_concentration_mol_per_m3 = try allocator.alloc(f64, try std.math.mul(usize, state.layer_count, aqueous_species_count));
    defer allocator.free(soil_recharge_concentration_mol_per_m3);
    const soil_solute_boundary_net_flux_mol = try allocator.alloc(f64, try std.math.mul(usize, state.layer_count, aqueous_species_count));
    defer allocator.free(soil_solute_boundary_net_flux_mol);
    @memset(soil_solute_boundary_net_flux_mol, 0);
    var gas_transport_state = try ecosys.gas_transport.State.init(allocator, state.layer_count);
    defer gas_transport_state.deinit();
    var soil_gas_transport_state = try ecosys.soil_gas_transport_step.State.init(allocator, state.layer_count);
    defer soil_gas_transport_state.deinit();
    var surface_litter_gas_transport_state = try ecosys.surface_litter_gas_transport_step.State.init(allocator, state.cell_count);
    defer surface_litter_gas_transport_state.deinit();
    var litter_gas_transport_state = try ecosys.gas_transport.State.init(allocator, state.cell_count);
    defer litter_gas_transport_state.deinit();
    var plant_available_nutrient_state = try ecosys.soil_plant_available_nutrients.State.init(allocator, state.layer_count);
    defer plant_available_nutrient_state.deinit();
    try plant_available_nutrient_state.initializeMapped(&state, &soil_solver_property_state, soil_catalog.entries.items, catalog_index_by_cell, runscript.plant_nutrient_initialization);
    var soil_organic_state = try ecosys.soil_organic_initialization.State.init(allocator, state.layer_count);
    defer soil_organic_state.deinit();
    const soil_organic_carbon_at_hour_start_g_c = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_organic_carbon_at_hour_start_g_c);
    const soil_organic_carbon_change_g_c_per_h = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_organic_carbon_change_g_c_per_h);
    @memset(soil_organic_carbon_at_hour_start_g_c, 0);
    @memset(soil_organic_carbon_change_g_c_per_h, 0);
    var surface_organic_state = try ecosys.soil_organic_initialization.State.init(allocator, state.cell_count);
    defer surface_organic_state.deinit();
    const shoot_senescence_products_by_plant = try allocator.alloc(ecosys.canopy_photosynthesis.SenescenceProducts, output_plant_count);
    defer allocator.free(shoot_senescence_products_by_plant);
    @memset(shoot_senescence_products_by_plant, .{});
    const seasonal_turnover_event_by_plant = try allocator.alloc(bool, output_plant_count);
    defer allocator.free(seasonal_turnover_event_by_plant);
    @memset(seasonal_turnover_event_by_plant, false);
    const automatic_harvest_date_by_plant = try allocator.alloc(ecosys.plant_management.PackedDate, output_plant_count);
    defer allocator.free(automatic_harvest_date_by_plant);
    @memset(automatic_harvest_date_by_plant, .{ .day = 1, .month = 1, .year = 0 });
    const root_litter_products_by_plant = try allocator.alloc(ecosys.plant_root_metabolism.RootLitter, output_plant_count);
    defer allocator.free(root_litter_products_by_plant);
    @memset(root_litter_products_by_plant, std.mem.zeroes(ecosys.plant_root_metabolism.RootLitter));
    var root_litter_carbon_ledger = try ecosys.plant_root_litter_ledger.State.init(
        allocator,
        output_plant_count,
        ecosys.plant_root_system.biological_domain_count,
        config.soil_layers,
    );
    defer root_litter_carbon_ledger.deinit();
    var plant_salt_harvest_workspace = try ecosys.plant_salt_harvest_adapter.Workspace.init(
        allocator,
        config.soil_layers,
    );
    defer plant_salt_harvest_workspace.deinit();
    const cumulative_harvest_salt_mol_by_plant = try allocator.alloc(
        f64,
        output_plant_count * ecosys.plant_salt_harvest.salt_count,
    );
    defer allocator.free(cumulative_harvest_salt_mol_by_plant);
    @memset(cumulative_harvest_salt_mol_by_plant, 0);
    const litter_salt_mol_by_plant_exchange_layer = try allocator.alloc(
        f64,
        output_plant_count * (config.soil_layers + 1) * ecosys.plant_salt_harvest.salt_count,
    );
    defer allocator.free(litter_salt_mol_by_plant_exchange_layer);
    @memset(litter_salt_mol_by_plant_exchange_layer, 0);
    var plant_litter_salt_ingress_state = try ecosys.plant_litter_salt_ingress.State.init(
        allocator,
        state.cell_count,
        config.soil_layers,
    );
    defer plant_litter_salt_ingress_state.deinit();
    var litter_salt_publication_state = try ecosys.litter_salt_publication.State.init(
        allocator,
        state.cell_count,
        config.plant_populations,
        config.soil_layers,
    );
    defer litter_salt_publication_state.deinit();
    const root_litter_carbon_g_c_by_layer = try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(root_litter_carbon_g_c_by_layer);
    @memset(root_litter_carbon_g_c_by_layer, 0);
    const shoot_harvest_litter_carbon_g_c_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(shoot_harvest_litter_carbon_g_c_by_plant);
    @memset(shoot_harvest_litter_carbon_g_c_by_plant, 0);
    const hourly_manure_products_by_plant = try allocator.alloc(
        ecosys.grazing_manure.Products,
        output_plant_count,
    );
    defer allocator.free(hourly_manure_products_by_plant);
    @memset(hourly_manure_products_by_plant, .{});
    const shoot_harvest_litter_nitrogen_g_n_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(shoot_harvest_litter_nitrogen_g_n_by_plant);
    @memset(shoot_harvest_litter_nitrogen_g_n_by_plant, 0);
    const shoot_harvest_litter_phosphorus_g_p_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(shoot_harvest_litter_phosphorus_g_p_by_plant);
    @memset(shoot_harvest_litter_phosphorus_g_p_by_plant, 0);
    const harvest_carbon_at_hour_start_g_c_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(harvest_carbon_at_hour_start_g_c_by_plant);
    @memset(harvest_carbon_at_hour_start_g_c_by_plant, 0);
    const eroded_organic_component_count = try ecosys.soil_erosion_organic_bridge.componentCount(&soil_organic_state);
    var eroded_organic_workspace = try ecosys.eroded_constituents.PackedWorkspace.init(allocator, state.cell_count, eroded_organic_component_count);
    defer eroded_organic_workspace.deinit();
    var surface_litter_geometry_state = try ecosys.surface_litter_geometry_step.State.init(allocator, state.cell_count);
    defer surface_litter_geometry_state.deinit();
    var surface_litter_water_environment_state = try ecosys.surface_litter_water_environment.State.init(allocator, state.cell_count);
    defer surface_litter_water_environment_state.deinit();
    var surface_microbial_environment_state = try ecosys.surface_microbial_environment_step.State.init(allocator, state.cell_count);
    defer surface_microbial_environment_state.deinit();
    const surface_field_capacity_potential_megapascal = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_field_capacity_potential_megapascal);
    const surface_wilting_point_potential_megapascal = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_wilting_point_potential_megapascal);
    for (0..state.cell_count) |cell| {
        const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
        surface_field_capacity_potential_megapascal[cell] = profile.field_capacity_potential_megapascal;
        surface_wilting_point_potential_megapascal[cell] = profile.wilting_point_potential_megapascal;
    }
    const surface_litter_ice_m3 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_litter_ice_m3);
    @memset(surface_litter_ice_m3, 0);
    const surface_charcoal_carbon_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_charcoal_carbon_g_c);
    @memset(surface_charcoal_carbon_g_c, 0);
    var surface_litter_chemistry_state = try ecosys.surface_litter_chemistry.State.init(allocator, state.cell_count);
    defer surface_litter_chemistry_state.deinit();
    var surface_litter_fertilizer_state = try ecosys.surface_litter_fertilizer.State.init(allocator, state.cell_count);
    defer surface_litter_fertilizer_state.deinit();
    var surface_microbial_respiration_state = try ecosys.surface_microbial_respiration_step.State.init(allocator, state.cell_count);
    defer surface_microbial_respiration_state.deinit();
    var surface_microbial_oxygen_state = try ecosys.surface_microbial_oxygen_driver.State.init(allocator, state.cell_count);
    defer surface_microbial_oxygen_state.deinit();
    var surface_microbial_maintenance_state = try ecosys.surface_microbial_maintenance_step.State.init(allocator, state.cell_count);
    defer surface_microbial_maintenance_state.deinit();
    var surface_nonsymbiotic_nitrogen_fixation_state = try ecosys.surface_nonsymbiotic_nitrogen_fixation_step.State.init(allocator, state.cell_count);
    defer surface_nonsymbiotic_nitrogen_fixation_state.deinit();
    var surface_microbial_substrate_uptake_state = try ecosys.surface_microbial_substrate_uptake_step.State.init(allocator, state.cell_count);
    defer surface_microbial_substrate_uptake_state.deinit();
    var surface_denitrification_state = try ecosys.surface_denitrification_step.State.init(allocator, state.cell_count);
    defer surface_denitrification_state.deinit();
    var surface_microbial_assimilation_state = try ecosys.surface_microbial_assimilation_step.State.init(allocator, state.cell_count);
    defer surface_microbial_assimilation_state.deinit();
    var surface_microbial_mineral_exchange_state = try ecosys.surface_microbial_mineral_exchange_step.State.init(allocator, state.cell_count);
    defer surface_microbial_mineral_exchange_state.deinit();
    var surface_topsoil_mineral_exchange_state = try ecosys.surface_topsoil_mineral_exchange_step.State.init(allocator, state.cell_count);
    defer surface_topsoil_mineral_exchange_state.deinit();
    var surface_microbial_turnover_state = try ecosys.surface_microbial_turnover_step.State.init(allocator, state.cell_count);
    defer surface_microbial_turnover_state.deinit();
    var surface_organic_priming_state = try ecosys.surface_organic_priming_step.State.init(allocator, state.cell_count);
    defer surface_organic_priming_state.deinit();
    var surface_organic_decomposition_state = try ecosys.surface_organic_decomposition_step.State.init(allocator, state.cell_count);
    defer surface_organic_decomposition_state.deinit();
    var surface_organic_sorption_state = try ecosys.surface_organic_sorption_step.State.init(allocator, state.cell_count);
    defer surface_organic_sorption_state.deinit();
    var surface_litter_colonization_state = try ecosys.surface_litter_colonization_step.State.init(allocator, state.cell_count);
    defer surface_litter_colonization_state.deinit();
    const surface_humification_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_humification_fraction);
    const topsoil_humus_partition = try allocator.alloc([2]f64, state.cell_count);
    defer allocator.free(topsoil_humus_partition);
    for (0..state.cell_count) |cell| {
        const topsoil_material = soil_catalog.entries.items[catalog_index_by_cell[cell]].material;
        surface_humification_fraction[cell] = 0.150 + 0.300 * @min(0.333, topsoil_material.clay_mass_fraction[0]);
        // STARTS FCY=0.20 when no initialized humus inventory is available.
        // NITRO converts that humus fraction to its microbial-residue routing
        // fraction through 5*FC/(4*FC+1).
        const source_less_resistant_humus_fraction: f64 = 0.20;
        const source_microbial_litter_partition =
            5 * source_less_resistant_humus_fraction /
            (4 * source_less_resistant_humus_fraction + 1);
        topsoil_humus_partition[cell] = .{
            source_microbial_litter_partition,
            1 - source_microbial_litter_partition,
        };
    }
    var surface_litter_fertilizer_diagnostics = try ecosys.surface_litter_fertilizer_step.Diagnostics.init(allocator, state.cell_count);
    defer surface_litter_fertilizer_diagnostics.deinit();
    var soil_fertilizer_inventory = try ecosys.fertilizer_nitrogen_inventory.State.init(allocator, state.cell_count, config.soil_layers);
    defer soil_fertilizer_inventory.deinit();
    var mineral_fertilizer_inventory = try ecosys.mineral_fertilizer_inventory.State.init(allocator, state.cell_count, config.soil_layers);
    defer mineral_fertilizer_inventory.deinit();
    var eroded_fertilizer_workspace = try ecosys.eroded_constituents.PackedWorkspace.init(allocator, state.cell_count, ecosys.soil_erosion_fertilizer_bridge.component_count);
    defer eroded_fertilizer_workspace.deinit();
    var surface_litter_chemistry_diagnostics = try ecosys.surface_litter_chemistry_step.Diagnostics.init(allocator, state.cell_count);
    defer surface_litter_chemistry_diagnostics.deinit();
    const surface_litter_cation_selectivity = try allocator.alloc(ecosys.solute_cation_exchange.Selectivity, state.cell_count);
    defer allocator.free(surface_litter_cation_selectivity);
    for (surface_litter_cation_selectivity, 0..) |*selectivity, cell| {
        const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
        selectivity.* = .{
            .calcium_ammonium = profile.property(.gapon_calcium_ammonium)[0],
            .calcium_hydrogen = profile.property(.gapon_calcium_hydrogen)[0],
            .calcium_aluminum_and_iron = profile.property(.gapon_calcium_aluminum)[0],
            .calcium_magnesium = profile.property(.gapon_calcium_magnesium)[0],
            .calcium_sodium = profile.property(.gapon_calcium_sodium)[0],
            .calcium_potassium = profile.property(.gapon_calcium_potassium)[0],
        };
    }
    {
        const parameters = &organic_parameters;
        const depth_partition_factor = try allocator.alloc(f64, state.layer_count);
        defer allocator.free(depth_partition_factor);
        try ecosys.soil_organic_initialization.deriveMappedDepthPartition(
            allocator,
            &state,
            &soil_solver_property_state,
            soil_catalog.entries.items,
            catalog_index_by_cell,
            grid_environment.horizontal_cell_width_m,
            grid_environment.vertical_cell_width_m,
            terrain_hydrology_state.relative_surface_elevation_m,
            terrain_hydrology_state.minimum_surface_elevation_m,
            grid_environment.initial_water_table_depth_m,
            grid_environment.natural_water_table_surface_slope,
            parameters.humusDepthParameters(),
            depth_partition_factor,
        );
        try soil_organic_state.initializeMapped(&state, &soil_solver_property_state, soil_catalog.entries.items, catalog_index_by_cell, parameters.mapped(depth_partition_factor));
        try plant_available_nutrient_state.bindInitializedDissolvedOrganic(&soil_organic_state);
        try surface_organic_state.initializeMappedSurface(
            soil_catalog.entries.items,
            catalog_index_by_cell,
            grid_environment.horizontal_cell_width_m,
            grid_environment.vertical_cell_width_m,
            parameters.surface(),
        );
        for (0..state.cell_count) |cell| {
            const top = try state.layerIndex(cell, 0);
            const humus_base = (top * ecosys.soil_organic_initialization.substrate_count + 4) * ecosys.soil_organic_initialization.structural_fraction_count;
            const less_resistant_carbon_g_c = soil_organic_state.structural[humus_base].carbon_g_c;
            const more_resistant_carbon_g_c = soil_organic_state.structural[humus_base + 1].carbon_g_c;
            const total_humus_carbon_g_c = less_resistant_carbon_g_c + more_resistant_carbon_g_c;
            const less_resistant_fraction = if (total_humus_carbon_g_c > 0) less_resistant_carbon_g_c / total_humus_carbon_g_c else parameters.less_resistant_humus_fraction_at_surface;
            const microbial_litter_to_less_resistant = 5 * less_resistant_fraction / (4 * less_resistant_fraction + 1);
            topsoil_humus_partition[cell] = .{ microbial_litter_to_less_resistant, 1 - microbial_litter_to_less_resistant };
        }
    }
    for (0..state.layer_count) |layer| for (0..@min(runscript.microbial_substrate_count, ecosys.soil_organic_initialization.microbial_substrate_count)) |substrate| for (0..@min(runscript.microbial_population_count, ecosys.soil_organic_initialization.microbial_population_count)) |population| {
        const runtime_index = try soil_microbial_state.populationIndex(layer / config.soil_layers, layer % config.soil_layers, substrate, population);
        const organic_index = ((layer * ecosys.soil_organic_initialization.microbial_substrate_count + substrate) * ecosys.soil_organic_initialization.microbial_population_count + population) * ecosys.soil_organic_initialization.kinetic_fraction_count;
        const initialized_labile = soil_organic_state.microbial[organic_index];
        const initialized_resistant = soil_organic_state.microbial[organic_index + 1];
        const initialized_nonstructural = soil_organic_state.microbial[organic_index + 2];
        soil_microbial_state.structural[runtime_index * 2] = .{ .carbon_g_c = initialized_labile.carbon_g_c, .nitrogen_g_n = initialized_labile.nitrogen_g_n, .phosphorus_g_p = initialized_labile.phosphorus_g_p };
        soil_microbial_state.structural[runtime_index * 2 + 1] = .{ .carbon_g_c = initialized_resistant.carbon_g_c, .nitrogen_g_n = initialized_resistant.nitrogen_g_n, .phosphorus_g_p = initialized_resistant.phosphorus_g_p };
        soil_microbial_state.nonstructural[runtime_index] = .{ .carbon_g_c = initialized_nonstructural.carbon_g_c, .nitrogen_g_n = initialized_nonstructural.nitrogen_g_n, .phosphorus_g_p = initialized_nonstructural.phosphorus_g_p };
    };
    var initial_chemistry_state = try ecosys.solute_chemistry_state.State.init(allocator, state.layer_count);
    defer initial_chemistry_state.deinit();
    const soil_chemistry_layer_parameters = try allocator.alloc(ecosys.solute_chemistry_state.ReactionParameters, state.layer_count);
    defer allocator.free(soil_chemistry_layer_parameters);
    var soil_chemistry_solver_workspace = try ecosys.solute_reaction_solver.Workspace.init(allocator);
    defer soil_chemistry_solver_workspace.deinit();
    var eroded_chemistry_workspace = try ecosys.eroded_constituents.PackedWorkspace.init(allocator, state.cell_count, ecosys.soil_erosion_chemistry_bridge.component_count);
    defer eroded_chemistry_workspace.deinit();
    const soil_numerical_scales_state = try ecosys.soil_initial_numerical_scales.initialize(.{
        .calculation_floor = 1.0e-15,
        .geometry_floor_m = 1.0e-6,
        .initial_water_content_m3_per_m3 = 1.0e-3,
        .initial_ice_porosity_m3_per_m3 = 0.0,
        .pure_ice_density_megagrams_per_m3 = 0.92,
    });
    var initial_balance_ledger_state: ecosys.initial_balance_ledger.State = undefined;
    ecosys.initial_balance_ledger.initialize(&initial_balance_ledger_state);
    const chemistry_parameters = runscript.chemistry_initialization;
    {
        for (0..state.cell_count) |cell| {
            const soil_entry = soil_catalog.entries.items[catalog_index_by_cell[cell]];
            const profile = soil_entry.profile;
            for (0..profile.total_layer_count) |layer| {
                const layer_cell = try state.layerIndex(cell, layer);
                {
                    const primary_parameters = runscript.chemistry_primary_initialization;
                    const profile_chemistry_inputs = try ecosys.soil_chemistry_initialization.resolveProfileEquilibriumSentinels(.{
                        .soil_ph = profile.ph[layer],
                        .ammonium_g_n_per_megagram = soil_entry.material.initial_ammonium_g_per_megagram[layer],
                        .nitrate_g_n_per_megagram = profile.property(.nitrate_g_per_megagram)[layer],
                        .phosphate_g_p_per_megagram = profile.property(.phosphate_g_per_megagram)[layer],
                        .aluminum_g_per_megagram = profile.property(.aluminum_g_per_megagram)[layer],
                        .iron_g_per_megagram = profile.property(.iron_g_per_megagram)[layer],
                        .calcium_g_per_megagram = soil_entry.material.initial_calcium_g_per_megagram[layer],
                        .magnesium_g_per_megagram = profile.property(.magnesium_g_per_megagram)[layer],
                        .sodium_g_per_megagram = profile.property(.sodium_g_per_megagram)[layer],
                        .potassium_g_per_megagram = profile.property(.potassium_g_per_megagram)[layer],
                        .sulfate_sulfur_g_s_per_megagram = profile.property(.sulfate_sulfur_g_per_megagram)[layer],
                        .chloride_g_per_megagram = profile.property(.chloride_g_per_megagram)[layer],
                        .aluminum_phosphate_g_p_per_megagram = profile.property(.aluminum_phosphate_p_g_per_megagram)[layer],
                        .iron_phosphate_g_p_per_megagram = profile.property(.iron_phosphate_p_g_per_megagram)[layer],
                        .dicalcium_phosphate_g_p_per_megagram = profile.property(.calcium_hydrogen_phosphate_p_g_per_megagram)[layer],
                        .apatite_g_p_per_megagram = profile.property(.apatite_phosphorus_g_per_megagram)[layer],
                        .aluminum_hydroxide_g_al_per_megagram = profile.property(.aluminum_hydroxide_al_g_per_megagram)[layer],
                        .iron_hydroxide_g_fe_per_megagram = profile.property(.iron_hydroxide_fe_g_per_megagram)[layer],
                        .calcium_carbonate_g_ca_per_megagram = profile.property(.calcium_carbonate_ca_g_per_megagram)[layer],
                        .calcium_sulfate_g_ca_per_megagram = profile.property(.calcium_sulfate_ca_g_per_megagram)[layer],
                    }, primary_parameters);
                    try ecosys.soil_chemistry_initialization.seedProfilePrimaryState(&initial_chemistry_state, layer_cell, profile_chemistry_inputs, soil_solver_property_state.matrix_bulk_volume_m3[layer_cell] * soil_solver_property_state.bulk_density_megagrams_per_m3[layer_cell], state.matrix_liquid_water_m3[layer_cell], primary_parameters);
                }
            }
        }
        {
            const reaction_parameters = chemistry_reaction_parameters;
            const zone_volume_storage = try allocator.alloc(f64, 6 * state.layer_count);
            defer allocator.free(zone_volume_storage);
            @memset(zone_volume_storage, 0);
            const ammonium_non_band_volume_m3 = zone_volume_storage[0 * state.layer_count ..][0..state.layer_count];
            const ammonium_band_volume_m3 = zone_volume_storage[1 * state.layer_count ..][0..state.layer_count];
            const nitrate_non_band_volume_m3 = zone_volume_storage[2 * state.layer_count ..][0..state.layer_count];
            const nitrate_band_volume_m3 = zone_volume_storage[3 * state.layer_count ..][0..state.layer_count];
            const phosphate_non_band_volume_m3 = zone_volume_storage[4 * state.layer_count ..][0..state.layer_count];
            const phosphate_band_volume_m3 = zone_volume_storage[5 * state.layer_count ..][0..state.layer_count];
            const nutrient_zones = runscript.plant_nutrient_initialization;
            const fractions = ecosys.solute_charge_classification.ZoneFractions{
                .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
                .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
                .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - nutrient_zones.initial_phosphate_band_fraction,
                .phosphate_band = nutrient_zones.initial_phosphate_band_fraction,
            };
            for (0..state.cell_count) |cell| {
                const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
                for (0..profile.total_layer_count) |layer| {
                    const layer_cell = try state.layerIndex(cell, layer);
                    const water_volume_m3 = state.matrix_liquid_water_m3[layer_cell];
                    // REDIST-POND-ORGANIC-MASS, source `starts.f` 1297--1304.
                    // The carrier that intensive concentrations are made
                    // extensive against is chosen PER LAYER: mineral mass
                    // `BKVL` where a matrix exists, total layer volume `VOLT`
                    // where it does not. This call site was the third open-coded
                    // copy of the mineral-mass branch with no copy of the
                    // source's `ELSE`, so a ponded layer of open water (no
                    // mineral matrix, large water volume) was rejected as
                    // malformed input rather than carried on its volume. A6
                    // landed `selectOrganicCarrier` as the single owner and
                    // delegated the two in-module copies; this is the third.
                    // The acceptance predicate is NOT loosened: the owner still
                    // rejects non-finite and negative inputs, and a layer with
                    // neither mineral mass nor volume is still a hard error, so
                    // no check is lost. Passing `0` as the negligible-mass
                    // threshold preserves today's strict `> 0` reading and
                    // deliberately does not pre-empt the runtime `ZEROS`
                    // question open as `INIT-006`.
                    // Publish set: NONE. `soil_mass_megagrams` is a
                    // function-local `const` consumed by `prepareLayerZones` on
                    // the next lines, so no owner is replaced and there is no
                    // double-mutation risk. `layer_volume_m3` is a new READ
                    // only, written once at `solver_properties.zig:213`
                    // during initialization and not mutated before this point.
                    const soil_mass_megagrams = blk: {
                        const carrier = ecosys.soil_organic_initialization.selectOrganicCarrier(
                            soil_solver_property_state.matrix_bulk_volume_m3[layer_cell] * soil_solver_property_state.bulk_density_megagrams_per_m3[layer_cell],
                            soil_solver_property_state.layer_volume_m3[layer_cell],
                            0,
                        ) catch |err| {
                            std.log.err("STARTE requires a positive organic carrier (mineral mass or total layer volume): cell={d} layer={d} mineral_mass_megagrams={e} layer_volume_m3={e} water_volume_m3={e}", .{ cell, layer, soil_solver_property_state.matrix_bulk_volume_m3[layer_cell] * soil_solver_property_state.bulk_density_megagrams_per_m3[layer_cell], soil_solver_property_state.layer_volume_m3[layer_cell], water_volume_m3 });
                            return err;
                        };
                        break :blk carrier.amount;
                    };
                    if (!std.math.isFinite(water_volume_m3) or water_volume_m3 <= 0 or !std.math.isFinite(soil_mass_megagrams)) {
                        std.log.err("STARTE requires positive runtime soil and matrix-water volumes: cell={d} layer={d} soil_mass_megagrams={e} water_volume_m3={e}", .{ cell, layer, soil_mass_megagrams, water_volume_m3 });
                        return error.InvalidStarteLayerGeometry;
                    }
                    const prepared_zones = try ecosys.soil_fertilizer_dissolution.prepareLayerZones(.{
                        .water_volume_m3 = water_volume_m3,
                        .soil_mass_megagrams = soil_mass_megagrams,
                        .soil_volume_m3 = soil_solver_property_state.matrix_bulk_volume_m3[layer_cell],
                        .fractions = .{
                            .ammonium_non_band = fractions.ammonium_non_band,
                            .ammonium_band = fractions.ammonium_band,
                            .nitrate_non_band = fractions.nitrate_non_band,
                            .nitrate_band = fractions.nitrate_band,
                            .phosphate_non_band = fractions.phosphate_non_band,
                            .phosphate_band = fractions.phosphate_band,
                        },
                        .positive_soil_mass_threshold_megagrams = runscript.absolute_tolerance,
                    });
                    ammonium_non_band_volume_m3[layer_cell] = prepared_zones.ammonium_non_band_water_m3;
                    ammonium_band_volume_m3[layer_cell] = prepared_zones.ammonium_band_water_m3;
                    nitrate_non_band_volume_m3[layer_cell] = prepared_zones.nitrate_non_band_water_m3;
                    nitrate_band_volume_m3[layer_cell] = prepared_zones.nitrate_band_water_m3;
                    phosphate_non_band_volume_m3[layer_cell] = prepared_zones.phosphate_non_band_water_m3;
                    phosphate_band_volume_m3[layer_cell] = prepared_zones.phosphate_band_water_m3;
                    const shared_soil_mass_per_water_volume_megagrams_per_m3 = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(prepared_zones.whole_layer_normalization_basis, water_volume_m3);
                    const phosphate_non_band_ratio = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(prepared_zones.phosphate_non_band_normalization_basis, prepared_zones.phosphate_non_band_water_m3);
                    const phosphate_band_ratio = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(prepared_zones.phosphate_band_normalization_basis, prepared_zones.phosphate_band_water_m3);
                    const ammonium_non_band_ratio = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(prepared_zones.ammonium_non_band_normalization_basis, prepared_zones.ammonium_non_band_water_m3);
                    const ammonium_band_ratio = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(prepared_zones.ammonium_band_normalization_basis, prepared_zones.ammonium_band_water_m3);
                    const selectivity = ecosys.solute_cation_exchange.Selectivity{
                        .calcium_ammonium = profile.property(.gapon_calcium_ammonium)[layer],
                        .calcium_hydrogen = profile.property(.gapon_calcium_hydrogen)[layer],
                        .calcium_aluminum_and_iron = profile.property(.gapon_calcium_aluminum)[layer],
                        .calcium_magnesium = profile.property(.gapon_calcium_magnesium)[layer],
                        .calcium_sodium = profile.property(.gapon_calcium_sodium)[layer],
                        .calcium_potassium = profile.property(.gapon_calcium_potassium)[layer],
                    };
                    const cation_exchange_capacity_mol_charge_per_megagram = soil_solver_property_state.cation_exchange_capacity_mol_per_megagram[layer_cell];
                    const anion_exchange_capacity_mol_charge_per_megagram = 10 * profile.anion_exchange_capacity_cmol_kg[layer];
                    const total_carboxyl_sites_mol_per_megagram = reaction_parameters.surface_litter.carboxyl_sites_mol_per_megagram_c *
                        1.0e-6 * soil_solver_property_state.total_organic_carbon_g_per_megagram[layer_cell];
                    try ecosys.soil_chemistry_initialization.seedProfilePhosphateSurfaceSites(&initial_chemistry_state, layer_cell, profile.property(.phosphate_g_per_megagram)[layer] / 31.0, anion_exchange_capacity_mol_charge_per_megagram, reaction_parameters.phosphate_surface, reaction_parameters.phosphate_constants.h2po4);
                    try ecosys.soil_chemistry_initialization.seedProfileCationExchange(&initial_chemistry_state, layer_cell, cation_exchange_capacity_mol_charge_per_megagram, selectivity, fractions);
                    // SOLUTE-CARBOXYL-INIT-SEED, source `starte.f` 403,
                    // `XHC1=XCOOH*AMIN1(1.0,CHY1/DPCOH)`. Establishes the
                    // hydrogen-occupied carboxyl pool, which `State.init`
                    // correctly zeroes at allocation and which nothing then
                    // seeded. Zero occupancy is an ABSORBING state rather than
                    // merely a small one: the carboxyl reaction's substrate
                    // bound is the source's `XMIN=FIONX/BKVLW*XHC1`, which is
                    // proportional to the occupied pool, so a zero pool clamps
                    // the extent to exactly zero every iteration for any pH and
                    // any budget, and the organic proton buffer contributes
                    // identically nothing. That is the `carboxyl=0e0` term A9
                    // measured in four chemically dissimilar wave-1 examples.
                    // Same defect class as `SOLUTE-PHOSPHATE-EQUIL`: a term
                    // with no descent direction reports as a plateau and looks
                    // like an iteration-budget problem while not being one. No
                    // ceiling is raised and no tolerance is relaxed.
                    // Publish set is one element of one field, with no other
                    // writer on the initialization path; the layer-remap,
                    // erosion-bridge, pond-transfer and solver `commit` writers
                    // all run later and TRANSFORM an existing value, so this
                    // fills the hole they currently propagate rather than
                    // double-mutating anything. Ordering: must follow
                    // `seedProfilePrimaryState` (`:7700`), which sets the
                    // `aqueous[].hydrogen` that is the source's `CHY1`; this
                    // position satisfies that and matches source order, where
                    // STARTE 400--403 precede the `DO 1000 M=1,MRXN` loop.
                    try ecosys.soil_chemistry_initialization.seedProfileCarboxylOccupancy(
                        &initial_chemistry_state,
                        layer_cell,
                        total_carboxyl_sites_mol_per_megagram,
                        reaction_parameters.surface_litter.carboxyl_dissociation_constant,
                    );
                    initial_chemistry_state.water_mol_per_m3[layer_cell] = reaction_parameters.water_concentration_mol_per_m3;
                    const layer_parameters = reaction_parameters.forLayer(fractions, phosphate_non_band_ratio, phosphate_band_ratio, cation_exchange_capacity_mol_charge_per_megagram, total_carboxyl_sites_mol_per_megagram, .{ .shared_megagrams_per_m3 = shared_soil_mass_per_water_volume_megagrams_per_m3, .ammonium_non_band_megagrams_per_m3 = ammonium_non_band_ratio, .ammonium_band_megagrams_per_m3 = ammonium_band_ratio }, selectivity);
                    soil_chemistry_layer_parameters[layer_cell] = layer_parameters;
                    // STARTE lines 234-248: seed dissolved CO2 from atmospheric
                    // concentration and temperature-dependent Henry solubility
                    // before equilibration so bicarbonate/carbonate buffers are
                    // present for the initial charge-balance solve.
                    {
                        const gas_params = surface_gas_parameters;
                        const co2 = @intFromEnum(ecosys.gas_transport.Species.carbon_dioxide);
                        // Use mean annual air temperature (ATCA) like Fortran STARTE line 238.
                        // Soil temperature is 273.15 K here (thermal solver hasn't run yet).
                        const sol = try ecosys.gas_transport.surfaceSolubilityWaterToAir(mean_annual_temperature_k_by_cell[cell], gas_params.solubility);
                        initial_chemistry_state.aqueous[layer_cell].carbon_dioxide = @max(0, gas_params.atmospheric_concentration_g_per_m3[co2] * sol[co2] / ecosys.gas_transport.g_per_mol_tracked[co2]);
                    }
                    // SOLUTE-036 REVERTED by A0 2026-08-05. Kept as a comment so the
                    // next lane does not re-derive the same wrong binding.
                    //
                    // The request was to make the STARTE initial equilibrium use
                    // `starte.f:813`'s carboxyl rate limit `AMIN1(XHC1,AHY1)` rather
                    // than `solute.f:1418`'s. The rate limit reading is correct. The
                    // BINDING was not, for a structural reason:
                    //
                    // `starte.f:813` sits INSIDE `DO 1000 M=1,MRXN` (opened at
                    // `starte.f:408`, `MRXN=1000`), i.e. it is one rate-limited step
                    // of the iteration that converges toward equilibrium, applied
                    // once per iteration alongside every other reaction. The binding
                    // hoisted it to run ONCE, BEFORE the solver, then handed the
                    // mutated pool to the solver as an initial condition.
                    //
                    // Two consequences, both fatal:
                    //  1. Double mutation. `carboxyl_bound_hydrogen_mol_per_megagram`
                    //     is owned by the reaction solver, which equilibrates it. A
                    //     pre-pass writing the same field makes two writers of one
                    //     carrier, the failure mode `agent_workflow.md` section 8
                    //     names as worst.
                    //  2. It applies a per-iteration RATE LIMIT (`TADC`, and the
                    //     0.75 mol/Mg per-iteration cap above) as if it were a
                    //     one-shot seed, so the pool starts displaced from, not
                    //     closer to, the equilibrium the solver then has to find.
                    //
                    // Measured: Ottawa regressed from reaching
                    // `year1998-day1-hour24` to failing during INITIALIZATION with
                    // `SoluteReactionSolverDidNotConverge`, cell=11,
                    // `max_scaled_residual=2.94e7`, `newton_steps=1000
                    // picard_steps=0`. Bisect: `692b3ce` reaches hour 24
                    // (`deviation_per_m2=3.8914823262074827e-1`); `1c97ba0` does not
                    // reach hour 1. That all-Newton/no-Picard split is the same
                    // signature A4 traced to dead escape code in
                    // `SOLUTE-PHOSPHATE-EQUIL`, which is why it looked like a solver
                    // defect rather than a binding defect.
                    //
                    // If the `starte.f:813` limit is to be honoured, it belongs
                    // INSIDE the solver's iteration as one more rate-limited
                    // reaction, with the solver remaining the sole writer of the
                    // pool. That is a change to `solute_reaction_solver`, not a
                    // call-site addition here. SOLUTE-036 reopened.
                    const initialization_solver_options: ecosys.solute_reaction_solver.Options = .{
                        .absolute_tolerance = runscript.absolute_tolerance,
                        .picard_relaxation = runscript.picard_relaxation,
                        .max_iterations = iteration_limits.initial_solute_reaction_max_iterations,
                    };
                    _ = ecosys.solute_reaction_solver.solveCell(
                        allocator,
                        &initial_chemistry_state,
                        layer_cell,
                        layer_parameters,
                        initialization_solver_options,
                    ) catch |err| return ecosys.solute_failure_reporter
                        .reportPreservingSolverError(
                        allocator,
                        .{
                            .io = init.io,
                            .directory = std.Io.Dir.cwd(),
                            .file_path = "ecosys-ng-solute-initialization-failure.bin",
                            .context = .{
                                .global_cell_id = @intCast(cell),
                                .soil_layer_id = @intCast(layer),
                                .packed_cell_index = @intCast(layer_cell),
                            },
                        },
                        &initial_chemistry_state,
                        layer_cell,
                        layer_parameters,
                        initialization_solver_options,
                        err,
                    );
                }
            }
            try plant_available_nutrient_state.bindEquilibratedMineralPools(&initial_chemistry_state, .{
                .ammonium_non_band_m3 = ammonium_non_band_volume_m3,
                .ammonium_band_m3 = ammonium_band_volume_m3,
                .nitrate_non_band_m3 = nitrate_non_band_volume_m3,
                .nitrate_band_m3 = nitrate_band_volume_m3,
                .phosphate_non_band_m3 = phosphate_non_band_volume_m3,
                .phosphate_band_m3 = phosphate_band_volume_m3,
            });
        }
        for (0..state.cell_count) |cell| {
            const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
            try ecosys.soil_chemistry_initialization.seedSurfaceLitterFromTopsoil(
                &surface_litter_chemistry_state,
                cell,
                &initial_chemistry_state,
                try state.layerIndex(cell, 0),
                profile.surface_litter_ph,
                chemistry_parameters.water_activity_product_mol2_per_m6,
                runscript.dynamic_plant_salts,
            );
        }
    }
    for (0..state.layer_count) |layer_cell| {
        // STARTE lines 1418-1422 gate dissolved oxygen on depth: only layers
        // whose upper face lies above the water table (`CDPTH(L-1) < DTBLZ`)
        // receive `OXYS`; saturated soil below it starts anoxic. The oracle's
        // `CDPTH(L-1)` is this layer's top face, which the solver geometry
        // already stores as bottom depth minus thickness.
        const depth_cell = layer_cell / state.soil_layer_capacity;
        const layer_top_depth_m = soil_solver_property_state.layer_bottom_depth_m[layer_cell] -
            soil_solver_property_state.layer_thickness_m[layer_cell];
        micropore_solute_state.water_volume_m3[layer_cell] = state.matrix_liquid_water_m3[layer_cell];
        macropore_solute_state.water_volume_m3[layer_cell] = state.macropore_liquid_water_m3[layer_cell];
        gas_transport_state.temperature_k[layer_cell] = state.soil_temperature_k[layer_cell];
        gas_transport_state.air_volume_m3[layer_cell] = state.matrix_air_volume_m3[layer_cell];
        {
            const parameters = surface_gas_parameters;
            // STARTE lines 1418-1433: seed dissolved gas from atmospheric
            // concentrations at mean annual air temperature (ATCA), not current
            // soil temperature. Soil temperature is 273.15 K before the thermal
            // solver first runs; ATCA is a site-level constant already available.
            const cell = layer_cell / state.soil_layer_capacity;
            // STARTE 1419--1433 divide each dissolved concentration by
            // `exp(A*CSTR1)`. `CSTR1` is the layer's salt-derived ionic
            // strength, so this is exactly 1 for non-saline soil.
            const seed_zone_fractions: ecosys.solute_charge_classification.ZoneFractions = .{
                .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                .phosphate_band = runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
            };
            const layer_ionic_strength =
                (try initial_chemistry_state.activityCoefficients(layer_cell, seed_zone_fractions)).ionic_strength_mol_per_l;
            try ecosys.gas_transport.initializeSoilLayerCell(
                &gas_transport_state,
                layer_cell,
                state.matrix_air_volume_m3[layer_cell],
                state.matrix_liquid_water_m3[layer_cell],
                layer_top_depth_m,
                site_by_cell[depth_cell].initial_water_table_depth_m,
                mean_annual_temperature_k_by_cell[cell],
                parameters.atmospheric_concentration_g_per_m3,
                parameters.solubility,
                layer_ionic_strength,
                ecosys.gas_transport.starte_activity_coefficient,
            );
        }
    }
    try ecosys.soil_aqueous_transport_bridge.exportConcentrations(&initial_chemistry_state, soil_recharge_concentration_mol_per_m3);
    var snow_transport_state = try ecosys.snow_solute_transport.State.init(allocator, state.cell_count, runscript.snow_layer_bottom_depth_m.len);
    defer snow_transport_state.deinit();
    const snow_atmospheric_input_g = try allocator.alloc(f64, state.cell_count * ecosys.snow_solute_transport.species_count);
    defer allocator.free(snow_atmospheric_input_g);
    @memset(snow_atmospheric_input_g, 0);
    var atmospheric_solute_input_ledger_state =
        try ecosys.atmospheric_solute_input_ledger.State.init(
            allocator,
            state.cell_count,
        );
    defer atmospheric_solute_input_ledger_state.deinit();
    const snow_surface_partitions = try allocator.alloc(ecosys.snow_solute_transport.SurfacePartition, state.cell_count);
    defer allocator.free(snow_surface_partitions);
    const snow_surface_discharge = try allocator.alloc(ecosys.snow_solute_transport.SurfaceDischarge, state.cell_count);
    defer allocator.free(snow_surface_discharge);
    const direct_surface_solute_input = try allocator.alloc(ecosys.snow_solute_transport.SurfaceDischarge, state.cell_count);
    defer allocator.free(direct_surface_solute_input);
    @memset(direct_surface_solute_input, .{});
    const snowpack_internal_solute_flux_by_layer = try allocator.alloc(
        ecosys.snowpack_internal_solute_aggregation.SoluteFlux,
        try std.math.mul(
            usize,
            state.cell_count,
            snow_transport_state.layer_capacity,
        ),
    );
    defer allocator.free(snowpack_internal_solute_flux_by_layer);
    @memset(snowpack_internal_solute_flux_by_layer, .{});
    const snowpack_internal_solute_flux_workspace = try allocator.alloc(
        ecosys.snowpack_internal_solute_aggregation.SoluteFlux,
        try std.math.mul(
            usize,
            state.cell_count,
            snow_transport_state.layer_capacity,
        ),
    );
    defer allocator.free(snowpack_internal_solute_flux_workspace);
    @memset(snowpack_internal_solute_flux_workspace, .{});
    const snowpack_internal_salt_flux_mol_by_layer_species = try allocator.alloc(
        f64,
        try std.math.mul(
            usize,
            try std.math.mul(
                usize,
                state.cell_count,
                snow_transport_state.layer_capacity,
            ),
            ecosys.snowpack_internal_salt_aggregation.salt_species_count,
        ),
    );
    defer allocator.free(snowpack_internal_salt_flux_mol_by_layer_species);
    @memset(snowpack_internal_salt_flux_mol_by_layer_species, 0);
    const snowpack_internal_salt_flux_workspace_by_layer_species = try allocator.alloc(
        f64,
        try std.math.mul(
            usize,
            try std.math.mul(
                usize,
                state.cell_count,
                snow_transport_state.layer_capacity,
            ),
            ecosys.snowpack_internal_salt_aggregation.salt_species_count,
        ),
    );
    defer allocator.free(snowpack_internal_salt_flux_workspace_by_layer_species);
    @memset(snowpack_internal_salt_flux_workspace_by_layer_species, 0);
    var transport_hydrology_state = try ecosys.transport_hydrology.State.init(allocator, config.lon_count, config.lat_count, config.soil_layers, snow_transport_state.layer_capacity);
    defer transport_hydrology_state.deinit();
    try transport_hydrology_state.syncStorage(&state, &snow_transport_state);
    var soil_transport_faces = try ecosys.transport_hydrology.buildSoilFacesMapped(allocator, &transport_hydrology_state, &state, lateral_connection_mode_by_cell);
    defer soil_transport_faces.deinit();
    var soil_face_geometry_state = try ecosys.soil_face_geometry.State.initMapped(
        allocator,
        &state,
        &soil_transport_faces,
        soil_solver_property_state.layer_thickness_m,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
    );
    defer soil_face_geometry_state.deinit();
    var soil_solute_face_parameter_state = try ecosys.soil_solute_face_parameters.State.init(allocator, soil_transport_faces.micropore_faces.len);
    defer soil_solute_face_parameter_state.deinit();
    const soil_solute_face_component_count = try std.math.mul(
        usize,
        soil_transport_faces.micropore_faces.len,
        aqueous_species_count,
    );
    const micropore_solute_face_flux_mol = try allocator.alloc(
        f64,
        soil_solute_face_component_count,
    );
    defer allocator.free(micropore_solute_face_flux_mol);
    const macropore_solute_face_flux_mol = try allocator.alloc(
        f64,
        soil_solute_face_component_count,
    );
    defer allocator.free(macropore_solute_face_flux_mol);
    @memset(micropore_solute_face_flux_mol, 0);
    @memset(macropore_solute_face_flux_mol, 0);
    var soil_organic_face_parameter_state = try ecosys.soil_organic_face_parameters.State.init(allocator, soil_transport_faces.micropore_faces.len);
    defer soil_organic_face_parameter_state.deinit();
    var soil_dissolved_gas_face_parameter_state = try ecosys.soil_dissolved_gas_face_parameters.State.init(allocator, soil_transport_faces.micropore_faces.len);
    defer soil_dissolved_gas_face_parameter_state.deinit();
    var soil_dissolved_gas_transport_state = try ecosys.soil_dissolved_gas_transport.State.init(allocator, state.layer_count);
    defer soil_dissolved_gas_transport_state.deinit();
    const soil_dissolved_gas_recharge_concentration_g_per_m3 = try allocator.alloc(f64, state.layer_count * ecosys.gas_transport.species_count);
    defer allocator.free(soil_dissolved_gas_recharge_concentration_g_per_m3);
    @memset(soil_dissolved_gas_recharge_concentration_g_per_m3, 0);
    var soil_organic_transport_state = try ecosys.soil_organic_transport.State.init(allocator, state.layer_count);
    defer soil_organic_transport_state.deinit();
    try soil_organic_transport_state.initializeFromProfile(&soil_organic_state);
    const soil_organic_recharge_concentration_g_per_m3 = try allocator.alloc(f64, state.layer_count * ecosys.soil_organic_transport.component_count);
    defer allocator.free(soil_organic_recharge_concentration_g_per_m3);
    @memset(soil_organic_recharge_concentration_g_per_m3, 0);
    var mineral_nitrogen_transport_state = try ecosys.mineral_nitrogen_transport.State.init(allocator, state.layer_count);
    defer mineral_nitrogen_transport_state.deinit();
    var mineral_nitrogen_face_parameters = try ecosys.mineral_nitrogen_transport.FaceParameters.init(allocator, soil_transport_faces.micropore_faces.len);
    defer mineral_nitrogen_face_parameters.deinit();
    const mineral_nitrogen_zone_fractions: ecosys.mineral_nitrogen_transport.ZoneFractions = .{
        .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
        .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
    };
    try mineral_nitrogen_transport_state.initializeMatrix(
        &initial_chemistry_state,
        &soil_reactive_nitrogen_state,
        state.matrix_liquid_water_m3,
        mineral_nitrogen_zone_fractions,
        runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    var surface_transport_state = try ecosys.surface_solute_routing.State.init(allocator, config.lon_count, config.lat_count, aqueous_species_count);
    defer surface_transport_state.deinit();
    var redistribution_ledger = try ecosys.transport_redistribution.Ledger.init(
        allocator,
        micropore_solute_state.amount_mol.len,
        gas_transport_state.gaseous_mass_g.len,
        snow_transport_state.amount_g.len,
        surface_transport_state.amount_mol.len,
        micropore_solute_state.species_count,
        ecosys.gas_transport.species_count,
        ecosys.snow_solute_transport.species_count,
        aqueous_species_count,
    );
    defer redistribution_ledger.deinit();
    var soil_thermal_context: ecosys.soil_thermal.UpdateContext = .{ .thermal = &soil_thermal_state, .grid = &state, .liquid_water_heat_capacity_megajoules_per_m3_k = runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k, .ice_heat_capacity_megajoules_per_m3_k = runscript.soil_phase_heat_parameters.ice_heat_capacity_megajoules_per_m3_k };
    try runKernelAcrossSerialTilePlan(
        executor,
        &tile_plan.?,
        &soil_thermal_context,
        ecosys.soil_thermal.updateTile,
    );
    try soil_thermal_state.validateFinite();
    var soil_hourly_workspace = try ecosys.soil_hourly_workspace.State.init(allocator, state.layer_count);
    defer soil_hourly_workspace.deinit();
    try soil_hourly_workspace.refresh(&state, &soil_solver_property_state, &soil_thermal_state, &terrain_hydrology_state, runscript.soil_process_parameters);
    try soil_hourly_workspace.fillMacroporeFaceConductance(&soil_transport_faces, &soil_face_geometry_state, soil_face_geometry_state.macropore_hydraulic_conductance_m_per_h_megapascal);
    var soil_heat_solver_workspace = try ecosys.soil_heat_solver.Workspace.init(
        allocator,
        state.layer_count,
        soil_transport_faces.micropore_faces.len,
        config.tile_cells,
    );
    defer soil_heat_solver_workspace.deinit();
    var atmospheric_state = try ecosys.atmospheric_forcing.State.init(allocator, state.cell_count);
    defer atmospheric_state.deinit();
    const hourly_forcing_by_cell = try allocator.alloc(ecosys.weather.HourlyForcing, state.cell_count);
    defer allocator.free(hourly_forcing_by_cell);
    const hourly_radiation_by_cell = try allocator.alloc(ecosys.atmospheric_radiation.Result, state.cell_count);
    defer allocator.free(hourly_radiation_by_cell);
    const hourly_weather_header_by_cell = try allocator.alloc(ecosys.weather.Header, state.cell_count);
    defer allocator.free(hourly_weather_header_by_cell);
    const plant_calendar_by_cell = try allocator.alloc(
        ecosys.plant_development.Calendar,
        state.cell_count,
    );
    defer allocator.free(plant_calendar_by_cell);
    const hourly_weather_reference_height_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_weather_reference_height_m);
    const hourly_adjusted_shortwave_megajoules_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_adjusted_shortwave_megajoules_per_m2);
    const hourly_extraterrestrial_shortwave_megajoules_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_extraterrestrial_shortwave_megajoules_per_m2);
    const hourly_solar_angle_sine = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_solar_angle_sine);
    const hourly_solar_azimuth_radians = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_solar_azimuth_radians);
    const irrigation_water_depth_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(irrigation_water_depth_m);
    @memset(irrigation_water_depth_m, 0);
    const irrigation_dissolved_mass_g_per_m2 = try allocator.alloc(f64, state.cell_count * ecosys.irrigation_management_dispatch.dissolved_species_count);
    defer allocator.free(irrigation_dissolved_mass_g_per_m2);
    @memset(irrigation_dissolved_mass_g_per_m2, 0);
    const irrigation_hydrogen_mol_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(irrigation_hydrogen_mol_per_m2);
    @memset(irrigation_hydrogen_mol_per_m2, 0);
    var irrigation_loads = try ecosys.irrigation_layer_routing.Loads.init(
        allocator,
        state.cell_count,
        state.soil_layer_capacity,
    );
    defer irrigation_loads.deinit();
    const automatic_irrigation_depth_m_by_cell_hour = try allocator.alloc(f64, state.cell_count * 24);
    defer allocator.free(automatic_irrigation_depth_m_by_cell_hour);
    @memset(automatic_irrigation_depth_m_by_cell_hour, 0);
    const minimum_canopy_water_potential_mpa_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(minimum_canopy_water_potential_mpa_by_cell);
    @memset(minimum_canopy_water_potential_mpa_by_cell, 0);
    const soil_wilting_point_fraction = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_wilting_point_fraction);
    for (soil_solver_property_state.retention_curve, soil_wilting_point_fraction) |curve, *wilting| wilting.* = curve.curve.wilting_point_fraction;
    var canopy_radiation_state = try ecosys.canopy_radiation.State.init(allocator, state.cell_count);
    defer canopy_radiation_state.deinit();
    var canopy_optics_state: ?ecosys.canopy_optics.State = null;
    defer if (canopy_optics_state) |*optics_state| optics_state.deinit();
    if (plant_assignments) |assignments| canopy_optics_state = try ecosys.canopy_optics.State.initMapped(
        allocator,
        state.cell_count,
        config.plant_populations,
        assignments,
        plant_unit_by_cell.?,
        plant_catalog,
    );
    const canopy_cell_area_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(canopy_cell_area_m2);
    const surface_pond_minimum_heat_capacity_megajoules_per_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_pond_minimum_heat_capacity_megajoules_per_k);
    var total_grid_area_m2: f64 = 0;
    for (canopy_cell_area_m2, 0..) |*area_m2, cell| {
        area_m2.* = grid_environment.horizontal_cell_width_m[cell] *
            grid_environment.vertical_cell_width_m[cell];
        if (!std.math.isFinite(area_m2.*) or area_m2.* <= 0) return error.InvalidCanopyCellArea;
        surface_pond_minimum_heat_capacity_megajoules_per_k[cell] = runscript.surface_pond_activation_heat_capacity_megajoules_per_m2_k * area_m2.*;
        total_grid_area_m2 += area_m2.*;
    }
    if (!std.math.isFinite(total_grid_area_m2) or total_grid_area_m2 <= 0) return error.InvalidOutputGridArea;
    const initial_ground_air_temperature_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_ground_air_temperature_k);
    @memset(initial_ground_air_temperature_k, site.mean_annual_air_temperature_c + 273.15);
    const initial_ground_air_vapor_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_ground_air_vapor_fraction);
    const initial_ground_air_temperature = site.mean_annual_air_temperature_c + 273.15;
    const initial_ground_air_saturation = runscript.ground_air_parameters.saturation_vapor_prefactor_k / initial_ground_air_temperature * runscript.ground_air_parameters.saturation_relative_humidity * @exp(runscript.ground_air_parameters.saturation_temperature_k * (runscript.ground_air_parameters.saturation_reference_inverse_temperature_per_k - 1 / initial_ground_air_temperature));
    @memset(initial_ground_air_vapor_fraction, initial_ground_air_saturation);
    var surface_aerodynamic_state = try ecosys.surface_aerodynamics.State.init(allocator, state.cell_count, runscript.surface_aerodynamic_parameters.soil_roughness_height_m);
    defer surface_aerodynamic_state.deinit();
    var ground_air_state = try ecosys.ground_air_exchange.State.init(allocator, initial_ground_air_temperature_k, initial_ground_air_vapor_fraction, canopy_cell_area_m2, runscript.ground_air_parameters.minimum_air_column_height_m, runscript.ground_air_parameters);
    defer ground_air_state.deinit();
    const surface_total_canopy_area_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_total_canopy_area_m2);
    @memset(surface_total_canopy_area_m2, 0);
    const surface_canopy_height_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_canopy_height_m);
    @memset(surface_canopy_height_m, 0);
    const ground_air_vapor_pressure_kpa = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_vapor_pressure_kpa);
    @memset(ground_air_vapor_pressure_kpa, 0);
    const atmospheric_vapor_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(atmospheric_vapor_fraction);
    @memset(atmospheric_vapor_fraction, 0);
    const ground_air_canopy_resistance_h_per_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_canopy_resistance_h_per_m);
    @memset(ground_air_canopy_resistance_h_per_m, 0);
    const ground_air_sensible_source_megajoules_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_sensible_source_megajoules_per_h);
    @memset(ground_air_sensible_source_megajoules_per_h, 0);
    const runtime_plant_count = try std.math.mul(usize, state.cell_count, config.plant_populations);
    const grazing_average_shoot_carbon_g_c = try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(grazing_average_shoot_carbon_g_c);
    @memset(grazing_average_shoot_carbon_g_c, 0);
    const grazing_active_cell_count = try allocator.alloc(usize, config.plant_populations);
    defer allocator.free(grazing_active_cell_count);
    @memset(grazing_active_cell_count, 0);
    const grazing_total_shoot_carbon_g_c_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(grazing_total_shoot_carbon_g_c_by_cell);
    const grazing_total_grazer_biomass_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(grazing_total_grazer_biomass_by_cell);
    var plant_daily_flux_ledger = try ecosys.plant_daily_flux_ledger.State.init(allocator, runtime_plant_count);
    defer plant_daily_flux_ledger.deinit();
    var plant_root_soil_exchange_state = try ecosys.plant_root_soil_exchange_accumulation.State.init(allocator, runtime_plant_count);
    defer plant_root_soil_exchange_state.deinit();
    const delayed_live_canopy_combustion_heat_megajoules = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(delayed_live_canopy_combustion_heat_megajoules);
    @memset(delayed_live_canopy_combustion_heat_megajoules, 0);
    const delayed_standing_dead_combustion_heat_megajoules = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(delayed_standing_dead_combustion_heat_megajoules);
    @memset(delayed_standing_dead_combustion_heat_megajoules, 0);
    const delayed_subsurface_combustion_heat_megajoules = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(delayed_subsurface_combustion_heat_megajoules);
    @memset(delayed_subsurface_combustion_heat_megajoules, 0);
    const delayed_surface_combustion_heat_megajoules = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(delayed_surface_combustion_heat_megajoules);
    @memset(delayed_surface_combustion_heat_megajoules, 0);
    const surface_combustion_heat_megajoules_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_combustion_heat_megajoules_per_m2);
    @memset(surface_combustion_heat_megajoules_per_m2, 0);
    const ground_air_vapor_source_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_vapor_source_m3_per_h);
    @memset(ground_air_vapor_source_m3_per_h, 0);
    const ground_air_surface_sensible_conductance_megajoules_per_h_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_surface_sensible_conductance_megajoules_per_h_k);
    @memset(ground_air_surface_sensible_conductance_megajoules_per_h_k, 0);
    const ground_air_surface_vapor_conductance_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_surface_vapor_conductance_m3_per_h);
    @memset(ground_air_surface_vapor_conductance_m3_per_h, 0);
    const ground_air_surface_vapor_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_surface_vapor_fraction);
    @memset(ground_air_surface_vapor_fraction, 0);
    const ground_surface_evaporation_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_surface_evaporation_m3_per_h);
    @memset(ground_surface_evaporation_m3_per_h, 0);
    const ground_surface_condensation_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_surface_condensation_m3_per_h);
    @memset(ground_surface_condensation_m3_per_h, 0);
    const ground_surface_litter_water_change_m3 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_surface_litter_water_change_m3);
    @memset(ground_surface_litter_water_change_m3, 0);
    const ground_surface_topsoil_water_change_m3 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_surface_topsoil_water_change_m3);
    @memset(ground_surface_topsoil_water_change_m3, 0);
    var canopy_geometry = try ecosys.canopy_geometry.Geometry.init(allocator, runscript.canopy_discretization);
    defer canopy_geometry.deinit();
    var canopy_irradiance_interception_geometry_state = try ecosys.canopy_irradiance_interception_geometry.initialize(allocator, .{
        .leaf_inclination_sine = canopy_geometry.leaf_inclination_sine,
        .leaf_inclination_cosine = canopy_geometry.leaf_inclination_cosine,
        .sky_azimuth_class_count = runscript.canopy_discretization.diffuse_sky_sector_count,
        .leaf_azimuth_class_count = runscript.canopy_discretization.leaf_azimuth_class_count,
        .sky_elevation_rad = runscript.canopy_discretization.diffuse_sky_elevation_radians,
        .pi = std.math.pi,
        .reflected_angle_threshold_rad = -std.math.pi / 2.0,
    });
    defer canopy_irradiance_interception_geometry_state.deinit();
    var canopy_structure_state: ?ecosys.canopy_structure.State = null;
    defer if (canopy_structure_state) |*structure_state| structure_state.deinit();
    if (plant_assignments) |assignments| {
        canopy_structure_state = try ecosys.canopy_structure.State.initMapped(
            allocator,
            state.cell_count,
            config.plant_populations,
            canopy_geometry.leaf_inclination_sine.len,
            assignments,
            plant_unit_by_cell.?,
            plant_catalog,
        );
        var structure_context: ecosys.canopy_structure.ApplyContext = .{ .structure = &canopy_structure_state.?, .plants = &plant_state };
        try runKernelAcrossSerialTilePlan(
            executor,
            &tile_plan.?,
            &structure_context,
            ecosys.canopy_structure.applyLeafAreaTile,
        );
    }
    var detailed_canopy_state: ?ecosys.canopy_photosynthesis.State = null;
    defer if (detailed_canopy_state) |*canopy_state| canopy_state.deinit();
    var canopy_carbon_exchange_state: ?ecosys.canopy_carbon_exchange.State = null;
    defer if (canopy_carbon_exchange_state) |*ledger| ledger.deinit();
    var canopy_layer_distribution_state: ?ecosys.canopy_layer_distribution.State = null;
    defer if (canopy_layer_distribution_state) |*layer_state| layer_state.deinit();
    var canopy_precipitation_retention_state: ?ecosys.canopy_precipitation_retention.State = null;
    defer if (canopy_precipitation_retention_state) |*retention_state| retention_state.deinit();
    var surface_precipitation_state = try ecosys.surface_precipitation.RuntimeState.init(allocator, state.cell_count);
    defer surface_precipitation_state.deinit();
    {
        const parameters = surface_gas_parameters;
        for (0..state.cell_count) |cell| surface_precipitation_state.litter_water_m3[cell] = parameters.initial_litter_water_m3_per_g_c * try surface_organic_state.totalCarbon_g_c(cell);
        try surface_litter_chemistry_state.bindMineralReferenceWater(
            surface_precipitation_state.litter_water_m3,
        );
        var geometry_context: ecosys.surface_litter_geometry_step.ApplyContext = .{ .result = &surface_litter_geometry_state, .surface_organic = &surface_organic_state, .water_m3 = surface_precipitation_state.litter_water_m3, .ice_m3 = surface_litter_ice_m3, .charcoal_carbon_g_c = surface_charcoal_carbon_g_c, .parameters = parameters.litter_geometry };
        try runKernelAcrossSerialTilePlan(
            executor,
            &tile_plan.?,
            &geometry_context,
            ecosys.surface_litter_geometry_step.applyTile,
        );
        @memcpy(surface_precipitation_state.litter_water_capacity_m3, surface_litter_geometry_state.water_retention_capacity_m3);
        {
            // EXEC-002: STARTS line 1648 seeds `VOLW(0)=8.0E-06*ORGC(0)`, which
            // this reproduces faithfully, but that seed exceeds the litter
            // retention capacity derived from the litter geometry. That is not
            // itself a defect: the excess freezes on the first cold step and the
            // water is conserved. It matters because the surface freeze exposes
            // the latent-heat accounting defect described under EXEC-002 in the
            // discrepancy register, so a nonzero excess here predicts an early
            // freezing event and hence an early audit failure.
            var seeded_m3: f64 = 0;
            var capacity_m3: f64 = 0;
            for (0..state.cell_count) |cell| {
                seeded_m3 += surface_precipitation_state.litter_water_m3[cell];
                capacity_m3 += surface_precipitation_state.litter_water_capacity_m3[cell];
            }
            std.log.debug("litter water seed vs capacity: seeded_m3={e} capacity_m3={e} excess_m3={e}", .{ seeded_m3, capacity_m3, seeded_m3 - capacity_m3 });
        }
        for (0..state.cell_count) |cell| try ecosys.gas_transport.initializeSurfaceCell(&litter_gas_transport_state, cell, surface_litter_geometry_state.air_volume_m3[cell], surface_precipitation_state.litter_water_m3[cell], state.surface_temperature_k[cell], parameters.atmospheric_concentration_g_per_m3, parameters.solubility);
    }
    var plant_root_state: ?ecosys.plant_root_system.State = null;
    defer if (plant_root_state) |*root_state| root_state.deinit();
    var plant_root_nutrient_workspace: ?ecosys.plant_root_nutrient_uptake.GridWorkspace = null;
    defer if (plant_root_nutrient_workspace) |*workspace| workspace.deinit();
    var plant_root_salt_workspace: ?ecosys.plant_root_salt_exchange.GridWorkspace = null;
    defer if (plant_root_salt_workspace) |*workspace| workspace.deinit();
    var plant_root_exudation_workspace: ?ecosys.plant_root_exudation.GridWorkspace = null;
    defer if (plant_root_exudation_workspace) |*workspace| workspace.deinit();
    var plant_root_metabolism_workspace: ?ecosys.plant_root_metabolism.GridWorkspace = null;
    defer if (plant_root_metabolism_workspace) |*workspace| workspace.deinit();
    var plant_storage_remobilization_workspace: ?ecosys.plant_storage_remobilization.Workspace = null;
    defer if (plant_storage_remobilization_workspace) |*workspace| workspace.deinit();
    var plant_litter_partition_state: ?ecosys.plant_litter_partition.State = null;
    defer if (plant_litter_partition_state) |*partition_state| partition_state.deinit();
    var plant_phenology_state: ?ecosys.plant_phenology.State = null;
    defer if (plant_phenology_state) |*phenology_state| phenology_state.deinit();
    const phenology_root_oxygen_fraction = try allocator.alloc(f64, try std.math.mul(usize, state.cell_count, config.plant_populations));
    defer allocator.free(phenology_root_oxygen_fraction);
    @memset(phenology_root_oxygen_fraction, 0);
    var branch_development_state: ?ecosys.plant_phenology.BranchDevelopmentState = null;
    defer if (branch_development_state) |*development_state| development_state.deinit();
    var plant_growth_stage_state: ?ecosys.plant_growth_stages.State = null;
    defer if (plant_growth_stage_state) |*growth_stage_state| growth_stage_state.deinit();
    var plant_dormancy_state: ?ecosys.plant_dormancy.RuntimeState = null;
    defer if (plant_dormancy_state) |*dormancy_state| dormancy_state.deinit();
    const root_nutrient_traits = try allocator.alloc([3]ecosys.plant_traits.NutrientUptake, runtime_plant_count);
    defer allocator.free(root_nutrient_traits);
    const canopy_biochemistry_parameters = try allocator.alloc(ecosys.canopy_biochemistry.Parameters, runtime_plant_count);
    defer allocator.free(canopy_biochemistry_parameters);
    const shoot_growth_plant_parameters = try allocator.alloc(ecosys.shoot_growth_runtime.PlantParameters, runtime_plant_count);
    defer allocator.free(shoot_growth_plant_parameters);
    @memset(shoot_growth_plant_parameters, ecosys.shoot_growth_runtime.inactivePlantParameters());
    @memset(canopy_biochemistry_parameters, .{
        .pathway = .c3,
        .growth_habit = 0,
        .phenology_type = 0,
        .aboveground_turnover_type = 0,
        .rubisco_carboxylation_umol_per_g_protein_s = 1,
        .rubisco_oxygenation_umol_per_g_protein_s = 0,
        .pep_carboxylation_umol_per_g_protein_s = 0,
        .rubisco_co2_half_saturation_umol_per_l = 1,
        .rubisco_o2_half_saturation_umol_per_l = 1,
        .pep_co2_half_saturation_umol_per_l = 1,
        .rubisco_leaf_protein_fraction = 0,
        .pep_leaf_protein_fraction = 0,
        .chlorophyll_electron_transport_umol_per_g_protein_s = 0,
        .c3_chlorophyll_leaf_protein_fraction = 0,
        .c4_chlorophyll_leaf_protein_fraction = 0,
        .intercellular_to_atmospheric_co2_ratio = 0.7,
    });
    const inactive_nutrient_trait: ecosys.plant_traits.NutrientUptake = .{ .maximum_rate_g_per_m2_h = 0, .half_saturation_umol_per_l = 0, .minimum_concentration_umol_per_l = 0 };
    @memset(root_nutrient_traits, .{inactive_nutrient_trait} ** 3);
    const root_nutrient_feedback_enabled = try allocator.alloc(bool, runtime_plant_count);
    defer allocator.free(root_nutrient_feedback_enabled);
    @memset(root_nutrient_feedback_enabled, false);
    const root_metabolism_plant_parameters = try allocator.alloc(ecosys.plant_root_metabolism.RuntimePlantParameters, runtime_plant_count);
    defer allocator.free(root_metabolism_plant_parameters);
    const root_biological_domain_count_by_plant = try allocator.alloc(u8, runtime_plant_count);
    defer allocator.free(root_biological_domain_count_by_plant);
    @memset(root_biological_domain_count_by_plant, 1);
    @memset(root_metabolism_plant_parameters, .{
        .root_profile_type = 0,
        .mycorrhizal_type = 0,
        .growth_habit = 0,
        .leaf_phenology_type = 0,
        .root_growth_yield_g_c_per_g_c = 1,
        .root_nitrogen_to_carbon_g_n_per_g_c = 1,
        .root_phosphorus_to_carbon_g_p_per_g_c = 1,
        .stalk_nitrogen_to_carbon_g_n_per_g_c = 1,
        .stalk_phosphorus_to_carbon_g_p_per_g_c = 1,
        .primary_root_radius_m = 1,
        .secondary_root_radius_m = 1,
        .primary_specific_length_m_per_g_c = 1,
        .secondary_specific_length_m_per_g_c = 1,
        .secondary_root_branching_per_m = 1,
        .shoot_root_equilibration_fraction_per_h = 0,
    });
    const checkpoint_species_names = try allocator.alloc([]const u8, runtime_plant_count);
    defer allocator.free(checkpoint_species_names);
    @memset(checkpoint_species_names, "inactive");
    const checkpoint_metadata_cells = try allocator.alloc(ecosys.plant_checkpoint_metadata.CellView, state.cell_count);
    defer allocator.free(checkpoint_metadata_cells);
    @memset(checkpoint_metadata_cells, .{ .species_names = &.{}, .species_alive = &.{} });
    var plant_topology_controls = try ecosys.plant_topology.Controls.init(allocator, runtime_plant_count);
    defer plant_topology_controls.deinit();
    var plant_reproduction_controls = try ecosys.plant_reproduction.Controls.init(allocator, runtime_plant_count);
    defer plant_reproduction_controls.deinit();
    var canopy_layer_controls = try ecosys.canopy_layer_distribution.Controls.init(allocator, runtime_plant_count);
    defer canopy_layer_controls.deinit();
    const development_species_parameters = try allocator.alloc(ecosys.plant_growth_stages.SpeciesParameters, runtime_plant_count);
    defer allocator.free(development_species_parameters);
    @memset(development_species_parameters, .{ .growth_habit = .annual, .phenology_type = .evergreen, .photoperiod_type = .insensitive, .maturity_group_node_count = 1, .branch_floral_node_requirement = 1, .critical_photoperiod_h = 12, .photoperiod_sensitivity_h = 0, .determinate = true, .vegetative_stage_duration = runscript.phenology_parameters.vegetative_stage_duration, .reproductive_stage_duration = runscript.phenology_parameters.reproductive_stage_duration });
    const development_dormancy_parameters = try allocator.alloc(ecosys.plant_dormancy.Parameters, runtime_plant_count);
    defer allocator.free(development_dormancy_parameters);
    @memset(development_dormancy_parameters, .{ .required_leafout_h = 0, .required_leafoff_h = 0, .leafout_temperature_threshold_c = 0, .leafoff_temperature_threshold_c = 0, .chilling_temperature_c = 0, .drought_leafout_total_water_potential_megapascal = runscript.phenology_parameters.drought_leafout_total_water_potential_megapascal, .combined_leafout_turgor_potential_megapascal = runscript.phenology_parameters.minimum_turgor_potential_megapascal, .leafoff_total_water_potential_megapascal = runscript.phenology_parameters.vascular_leafoff_total_water_potential_megapascal, .maximum_photoperiod_counter_h = runscript.phenology_parameters.maximum_photoperiod_counter_h, .evergreen_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.evergreen_leafoff_remobilization_start_fraction, .deciduous_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.deciduous_leafoff_remobilization_start_fraction, .full_senescence_duration_h = runscript.root_metabolism_parameters.full_senescence_duration_h });
    const development_planting_dates = try allocator.alloc(ecosys.plant_management.PackedDate, runtime_plant_count);
    defer allocator.free(development_planting_dates);
    @memset(development_planting_dates, .{ .day = 1, .month = 1, .year = 9999 });
    const development_planting_day_of_year = try allocator.alloc(u16, runtime_plant_count);
    defer allocator.free(development_planting_day_of_year);
    @memset(development_planting_day_of_year, 1);
    const development_planting_year = try allocator.alloc(i32, runtime_plant_count);
    defer allocator.free(development_planting_year);
    @memset(development_planting_year, 0);
    const development_canopy_height_m = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(development_canopy_height_m);
    @memset(development_canopy_height_m, 0);
    const development_surface_water_potential_megapascal = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(development_surface_water_potential_megapascal);
    @memset(development_surface_water_potential_megapascal, 0);
    const development_seed_layer_water_potential_megapascal = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(development_seed_layer_water_potential_megapascal);
    @memset(development_seed_layer_water_potential_megapascal, 0);
    const newly_activated_plants = try allocator.alloc(bool, runtime_plant_count);
    defer allocator.free(newly_activated_plants);
    @memset(newly_activated_plants, false);
    const reseed_population_per_m2_by_plant = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(reseed_population_per_m2_by_plant);
    @memset(reseed_population_per_m2_by_plant, 0);
    const replant_layer_bottom_depth_m = try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(replant_layer_bottom_depth_m);
    var management_schedule_map: ?ecosys.plant_management_dispatch.ScheduleMap = null;
    defer if (management_schedule_map) |*schedule_map| schedule_map.deinit();
    var harvest_science_by_plant: ?[]ecosys.plant_harvest_runtime.ScienceParameters = null;
    defer if (harvest_science_by_plant) |science| allocator.free(science);
    var harvest_products_by_plant: ?[]ecosys.plant_harvest_runtime.ProductLedger = null;
    defer if (harvest_products_by_plant) |products| allocator.free(products);
    var harvest_context: ?ecosys.plant_harvest_runtime.Context = null;
    if (plant_assignments) |assignments| {
        const plant_count = try std.math.mul(usize, state.cell_count, config.plant_populations);
        plant_root_state = try ecosys.plant_root_system.State.init(allocator, plant_count, config.soil_layers, runscript.root_axes_per_plant);
        plant_root_nutrient_workspace = try ecosys.plant_root_nutrient_uptake.GridWorkspace.initWithAdmissionCapacity(
            allocator,
            state.cell_count,
            try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count),
            try std.math.mul(usize, try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count), config.soil_layers),
        );
        plant_root_salt_workspace = try ecosys.plant_root_salt_exchange.GridWorkspace.initWithAdmissionCapacity(
            allocator,
            state.cell_count,
            try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count),
            try std.math.mul(usize, try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count), config.soil_layers),
        );
        plant_root_exudation_workspace = try ecosys.plant_root_exudation.GridWorkspace.initWithAdmissionCapacity(
            allocator,
            state.cell_count,
            try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count),
            try std.math.mul(usize, try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count), config.soil_layers),
        );
        plant_root_metabolism_workspace = try ecosys.plant_root_metabolism.GridWorkspace.init(allocator, state.cell_count, runscript.root_axes_per_plant, config.soil_layers);
        plant_storage_remobilization_workspace = try ecosys.plant_storage_remobilization.Workspace.init(allocator, plant_count, config.soil_layers);
        plant_litter_partition_state = try ecosys.plant_litter_partition.State.init(allocator, plant_count);
        plant_phenology_state = try ecosys.plant_phenology.State.init(allocator, state.cell_count, config.plant_populations);
        const branch_count_by_plant = try allocator.alloc(usize, plant_count);
        defer allocator.free(branch_count_by_plant);
        @memset(branch_count_by_plant, 0);
        for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
            const active_species_count = assignments.units[unit_index].species.len;
            for (0..active_species_count) |species| branch_count_by_plant[cell * config.plant_populations + species] = 1;
        }
        var branch_count: usize = 0;
        for (branch_count_by_plant) |count| branch_count = try std.math.add(usize, branch_count, count);
        if (branch_count > 0) {
            plant_growth_stage_state = try ecosys.plant_growth_stages.State.init(allocator, branch_count_by_plant);
            plant_dormancy_state = try ecosys.plant_dormancy.RuntimeState.init(allocator, branch_count);
            const node_count_by_branch = try allocator.alloc(usize, branch_count);
            defer allocator.free(node_count_by_branch);
            @memset(node_count_by_branch, 1);
            const sample_count_by_node = try allocator.alloc(usize, branch_count);
            defer allocator.free(sample_count_by_node);
            const angular_sample_count = try std.math.mul(usize, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len);
            const radiation_sample_count_per_node = try std.math.mul(usize, runscript.canopy_layer_count, angular_sample_count);
            @memset(sample_count_by_node, radiation_sample_count_per_node);
            detailed_canopy_state = try ecosys.canopy_photosynthesis.State.init(allocator, state.cell_count, config.plant_populations, branch_count_by_plant, node_count_by_branch, sample_count_by_node);
            canopy_carbon_exchange_state = try ecosys.canopy_carbon_exchange.State.init(allocator, branch_count);
            canopy_layer_distribution_state = try ecosys.canopy_layer_distribution.State.init(allocator, state.cell_count, config.plant_populations, runscript.canopy_layer_count, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len, &detailed_canopy_state.?);
            canopy_precipitation_retention_state = try ecosys.canopy_precipitation_retention.State.init(allocator, state.cell_count, config.plant_populations);
            for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
                const cell_area_m2 = grid_environment.horizontal_cell_width_m[cell] *
                    grid_environment.vertical_cell_width_m[cell];
                const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
                const layer_bottom_depth_m = try allocator.alloc(f64, profile.total_layer_count);
                defer allocator.free(layer_bottom_depth_m);
                var cumulative_depth_m: f64 = 0;
                for (layer_bottom_depth_m, 0..) |*bottom, layer| {
                    cumulative_depth_m += soil_solver_property_state.layer_thickness_m[cell * config.soil_layers + layer];
                    bottom.* = cumulative_depth_m;
                }
                for (assignments.units[unit_index].species, 0..) |assignment, species| {
                    if (ecosys.delimited_input.isNo(assignment.management_file)) continue;
                    const trait_index = plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile;
                    const management_index = plant_management_catalog.find(assignment.management_file) orelse return error.MissingPlantManagementSchedule;
                    const traits = plant_catalog.entries.items[trait_index].traits;
                    const planting = plant_management_catalog.entries.items[management_index].schedule.planting;
                    const planting_layer = try ecosys.plant_initialization.plantingLayer(planting.seed_depth_m, layer_bottom_depth_m);
                    const storage = try ecosys.plant_initialization.seedStorage(
                        traits.morphology.seed_mass_at_planting_g,
                        planting.population_per_m2,
                        cell_area_m2,
                        traits.organ_nitrogen_to_carbon_ratio.grain,
                        traits.organ_phosphorus_to_carbon_ratio.grain,
                    );
                    const plant = try detailed_canopy_state.?.plantIndex(cell, species);
                    reseed_population_per_m2_by_plant[plant] = planting.population_per_m2;
                    root_nutrient_traits[plant] = .{ traits.ammonium_uptake, traits.nitrate_uptake, traits.phosphate_uptake };
                    canopy_biochemistry_parameters[plant] = .{
                        .pathway = if (traits.functional_type.photosynthesis_pathway == 3) .c3 else .c4,
                        .growth_habit = traits.functional_type.growth_habit,
                        .phenology_type = traits.functional_type.leaf_phenology_type,
                        .aboveground_turnover_type = traits.functional_type.aboveground_turnover_type,
                        .rubisco_carboxylation_umol_per_g_protein_s = traits.photosynthesis.rubisco_carboxylase_umol_c_per_g_enzyme_s,
                        .rubisco_oxygenation_umol_per_g_protein_s = traits.photosynthesis.rubisco_oxygenase_umol_o_per_g_enzyme_s,
                        .pep_carboxylation_umol_per_g_protein_s = traits.photosynthesis.pep_carboxylase_umol_per_g_enzyme_s,
                        .rubisco_co2_half_saturation_umol_per_l = traits.photosynthesis.co2_half_saturation_umol_per_l,
                        .rubisco_o2_half_saturation_umol_per_l = traits.photosynthesis.oxygen_half_saturation_umol_per_l,
                        .pep_co2_half_saturation_umol_per_l = traits.photosynthesis.pep_co2_half_saturation_umol_per_l,
                        .rubisco_leaf_protein_fraction = traits.photosynthesis.rubisco_leaf_protein_fraction,
                        .pep_leaf_protein_fraction = traits.photosynthesis.pep_carboxylase_leaf_protein_fraction,
                        .chlorophyll_electron_transport_umol_per_g_protein_s = traits.photosynthesis.chlorophyll_electron_activity_umol_per_g_enzyme_s,
                        .c3_chlorophyll_leaf_protein_fraction = traits.photosynthesis.mesophyll_chlorophyll_leaf_protein_fraction,
                        .c4_chlorophyll_leaf_protein_fraction = traits.photosynthesis.c4_mesophyll_chlorophyll_leaf_protein_fraction,
                        .intercellular_to_atmospheric_co2_ratio = traits.photosynthesis.intercellular_to_atmospheric_co2_ratio,
                    };
                    try canopy_biochemistry_parameters[plant].validate();
                    shoot_growth_plant_parameters[plant] = ecosys.shoot_growth_runtime.parametersFromTraits(traits);
                    root_nutrient_feedback_enabled[plant] = traits.functional_type.growth_habit != 0;
                    const root_geometry = try ecosys.plant_initialization.rootGeometry(traits.roots.root_porosity_fraction, traits.roots.primary_root_radius_m, traits.roots.secondary_root_radius_m, runscript.plant_geometry_parameters);
                    root_metabolism_plant_parameters[plant] = .{
                        .root_profile_type = traits.functional_type.root_profile_type,
                        .mycorrhizal_type = traits.functional_type.mycorrhizal_type,
                        .growth_habit = if (traits.functional_type.growth_habit == 0) 0 else 1,
                        .leaf_phenology_type = traits.functional_type.leaf_phenology_type,
                        .root_growth_yield_g_c_per_g_c = traits.organ_growth_yield_g_c_per_g_c.root,
                        .root_nitrogen_to_carbon_g_n_per_g_c = traits.organ_nitrogen_to_carbon_ratio.root,
                        .root_phosphorus_to_carbon_g_p_per_g_c = traits.organ_phosphorus_to_carbon_ratio.root,
                        .stalk_nitrogen_to_carbon_g_n_per_g_c = traits.organ_nitrogen_to_carbon_ratio.stalk,
                        .stalk_phosphorus_to_carbon_g_p_per_g_c = traits.organ_phosphorus_to_carbon_ratio.stalk,
                        .primary_root_radius_m = traits.roots.primary_root_radius_m,
                        .secondary_root_radius_m = traits.roots.secondary_root_radius_m,
                        .primary_specific_length_m_per_g_c = root_geometry.primary_specific_length_m_per_g_c,
                        .secondary_specific_length_m_per_g_c = root_geometry.secondary_specific_length_m_per_g_c,
                        .secondary_root_branching_per_m = traits.roots.secondary_root_branching_per_m,
                        .shoot_root_equilibration_fraction_per_h = traits.roots.shoot_root_carbon_equilibration_fraction_per_h,
                    };
                    try root_metabolism_plant_parameters[plant].validate();
                    root_biological_domain_count_by_plant[plant] = @intCast(root_metabolism_plant_parameters[plant].biologicalDomainCount());
                    const adjusted_phenology = try ecosys.plant_initialization.adjustedPhenology(traits, runscript.phenology_initialization_parameters);
                    try plant_topology_controls.setPlant(plant, traits.phenology.branching_nonstructural_carbon_fraction, traits.roots.branching_nonstructural_carbon_fraction, traits.functional_type.growth_habit, traits.functional_type.leaf_phenology_type, adjusted_phenology.floral_initiation_node_count_after_seed);
                    try plant_reproduction_controls.setPlant(plant, traits.morphology.maximum_seed_sites, traits.morphology.maximum_seeds_per_site, traits.morphology.maximum_seed_mass_g, traits.phenology.chilling_temperature_c, traits.water_relations.stomatal_turgor_shape, traits.functional_type.root_profile_type == 0);
                    try canopy_layer_controls.setPlant(plant, traits.phenology.leaf_length_to_width_ratio, runscript.stalk_volume_m3_per_g_c, traits.functional_type.aboveground_turnover_type, traits.functional_type.root_profile_type, traits.functional_type.growth_habit == 0, traits.morphology.specific_internode_length_m_per_g_c, @sin(traits.morphology.stem_angle_degrees * std.math.pi / 180.0));
                    const thermal_acclimation = try ecosys.plant_initialization.thermalAcclimation(
                        assignment.species_file,
                        traits.functional_type.photosynthesis_pathway,
                        traits.functional_type.thermal_adaptation_zone,
                        runscript.thermal_acclimation_parameters,
                    );
                    development_species_parameters[plant] = .{
                        .growth_habit = ecosys.plant_growth_stages.growthHabitFromReadq(traits.functional_type.growth_habit),
                        .phenology_type = try ecosys.plant_growth_stages.phenologyTypeFromReadq(traits.functional_type.leaf_phenology_type),
                        .photoperiod_type = try ecosys.plant_growth_stages.photoperiodTypeFromReadq(traits.functional_type.photoperiod_type),
                        .maturity_group_node_count = adjusted_phenology.floral_initiation_node_count_after_seed,
                        .branch_floral_node_requirement = adjusted_phenology.floral_initiation_node_count_after_seed,
                        .critical_photoperiod_h = traits.phenology.critical_photoperiod_h,
                        .photoperiod_sensitivity_h = traits.phenology.floral_induction_photoperiod_difference_h,
                        .determinate = traits.functional_type.determinacy_type == 0,
                        .vegetative_stage_duration = runscript.phenology_parameters.vegetative_stage_duration,
                        .reproductive_stage_duration = runscript.phenology_parameters.reproductive_stage_duration,
                    };
                    development_dormancy_parameters[plant] = .{
                        .required_leafout_h = traits.phenology.spring_leafout_requirement_h,
                        .required_leafoff_h = traits.phenology.autumn_leafoff_requirement_h,
                        .leafout_temperature_threshold_c = thermal_acclimation.leafout_threshold_c,
                        .leafoff_temperature_threshold_c = thermal_acclimation.leafoff_threshold_c,
                        .chilling_temperature_c = traits.phenology.chilling_temperature_c,
                        .drought_leafout_total_water_potential_megapascal = runscript.phenology_parameters.drought_leafout_total_water_potential_megapascal,
                        .combined_leafout_turgor_potential_megapascal = runscript.phenology_parameters.minimum_turgor_potential_megapascal,
                        .leafoff_total_water_potential_megapascal = if (traits.functional_type.root_profile_type == 0) runscript.phenology_parameters.nonvascular_leafoff_total_water_potential_megapascal else runscript.phenology_parameters.vascular_leafoff_total_water_potential_megapascal,
                        .maximum_photoperiod_counter_h = runscript.phenology_parameters.maximum_photoperiod_counter_h,
                        .evergreen_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.evergreen_leafoff_remobilization_start_fraction,
                        .deciduous_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.deciduous_leafoff_remobilization_start_fraction,
                        .full_senescence_duration_h = runscript.root_metabolism_parameters.full_senescence_duration_h,
                    };
                    development_planting_dates[plant] = planting.date;
                    plant_phenology_state.?.active[plant] = false;
                    plant_phenology_state.?.annual_growth_habit[plant] = traits.functional_type.growth_habit == 0;
                    plant_phenology_state.?.node_initiation_rate_at_25c_per_h[plant] = adjusted_phenology.node_initiation_per_h;
                    plant_phenology_state.?.leaf_appearance_rate_at_25c_per_h[plant] = adjusted_phenology.leaf_appearance_per_h;
                    plant_phenology_state.?.initiated_node_count[plant] = adjusted_phenology.seed_initial_node_count;
                    plant_phenology_state.?.appeared_leaf_count[plant] = adjusted_phenology.seed_initial_node_count;
                    const growth_stage_range = try plant_growth_stage_state.?.branchRange(plant);
                    if (growth_stage_range.first < growth_stage_range.end) {
                        plant_growth_stage_state.?.branches[growth_stage_range.first].initiated_node_count = adjusted_phenology.seed_initial_node_count;
                        plant_growth_stage_state.?.branches[growth_stage_range.first].appeared_leaf_count = adjusted_phenology.seed_initial_node_count;
                    }
                    try plant_root_state.?.initializePlant(plant, traits, planting_layer, planting.seed_depth_m, runscript.root_initialization_parameters);
                    const seed_geometry = try ecosys.plant_initialization.seedGeometry(traits.morphology.seed_mass_at_planting_g, runscript.plant_geometry_parameters);
                    try plant_root_state.?.setSeedGeometry(plant, seed_geometry.volume_m3, seed_geometry.length_m, seed_geometry.surface_area_m2);
                    try plant_litter_partition_state.?.initializePlant(plant, traits, runscript.standing_dead_partition_parameters);
                    try ecosys.plant_initialization.initializeShootControls(
                        &detailed_canopy_state.?,
                        plant,
                        planting.population_per_m2,
                        cell_area_m2,
                        traits.water_relations.cuticular_resistance_s_per_m,
                        traits.functional_type.photosynthesis_pathway,
                        runscript.shoot_control_parameters,
                    );
                    try ecosys.plant_initialization.initializeThermalAcclimation(
                        &detailed_canopy_state.?,
                        plant,
                        thermal_acclimation,
                    );
                    try ecosys.plant_initialization.initializeSeedStorage(&detailed_canopy_state.?, plant, storage);
                    const standing_dead = try ecosys.plant_initialization.standingDeadStorage(
                        traits.morphology.standing_dead_carbon_g_per_m2,
                        cell_area_m2,
                        traits.organ_nitrogen_to_carbon_ratio.stalk,
                        traits.organ_phosphorus_to_carbon_ratio.stalk,
                        runscript.standing_dead_partition_parameters,
                    );
                    try ecosys.plant_initialization.initializeStandingDeadStorage(&detailed_canopy_state.?, plant, standing_dead);
                    try ecosys.plant_initialization.initializePlantHeatAndWater(
                        &plant_state,
                        &detailed_canopy_state.?,
                        plant,
                        site.mean_annual_air_temperature_c,
                        traits.water_relations.osmotic_potential_megapascal,
                        runscript.plant_heat_water_parameters,
                    );
                }
            }
            branch_development_state = try ecosys.plant_phenology.BranchDevelopmentState.init(allocator, branch_count);
            for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
                for (assignments.units[unit_index].species, 0..) |assignment, species| {
                    if (ecosys.delimited_input.isNo(assignment.management_file)) continue;
                    const trait_index = plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile;
                    const traits = plant_catalog.entries.items[trait_index].traits;
                    const plant = try detailed_canopy_state.?.plantIndex(cell, species);
                    const branch_range = try detailed_canopy_state.?.branchRange(plant);
                    const adjusted = try ecosys.plant_initialization.adjustedPhenology(traits, runscript.phenology_initialization_parameters);
                    const concurrent = try ecosys.plant_initialization.concurrentNodeSettings(traits.functional_type.aboveground_turnover_type, adjusted.leaf_appearance_per_h, adjusted.floral_initiation_node_count_after_seed, runscript.phenology_initialization_parameters);
                    try branch_development_state.?.initializeRange(
                        branch_range.first,
                        branch_range.end,
                        adjusted.floral_initiation_node_count_after_seed,
                        adjusted.seed_initial_node_count,
                        traits.functional_type.growth_habit == 0 and traits.functional_type.leaf_phenology_type == 1,
                        concurrent.perennial_node_scaling,
                        concurrent.maximum_concurrently_growing_nodes,
                    );
                }
            }
            harvest_science_by_plant = try allocator.alloc(ecosys.plant_harvest_runtime.ScienceParameters, plant_count);
            @memset(harvest_science_by_plant.?, .{
                .carbon_woody_fraction = .{ 0, 1 },
                .leaf_nitrogen_woody_fraction = .{ 0, 1 },
                .sheath_nitrogen_woody_fraction = .{ 0, 1 },
                .leaf_phosphorus_woody_fraction = .{ 0, 1 },
                .sheath_phosphorus_woody_fraction = .{ 0, 1 },
            });
            for (harvest_science_by_plant.?, shoot_growth_plant_parameters) |*harvest_science, plant_parameters|
                harvest_science.nitrogen_fixation_type = plant_parameters.nitrogen_fixation_type;
            harvest_products_by_plant = try allocator.alloc(ecosys.plant_harvest_runtime.ProductLedger, plant_count);
            @memset(harvest_products_by_plant.?, .{});
            management_schedule_map = try ecosys.plant_management_dispatch.ScheduleMap.initMapped(allocator, state.cell_count, config.plant_populations, assignments, plant_unit_by_cell.?, plant_management_catalog);
            harvest_context = .{
                .canopy_state = &detailed_canopy_state.?,
                .canopy_structure_state = &canopy_structure_state.?,
                .canopy_layer_state = &canopy_layer_distribution_state.?,
                .branch_development = &branch_development_state.?,
                .science_by_plant = harvest_science_by_plant.?,
                .products_by_plant = harvest_products_by_plant.?,
                .leaf_area_presence_tolerance_m2 = config.negligible_quantity_threshold,
                .plant_structural_presence_threshold_g_per_plant = runscript.plant_pool_parameters.branch_structural_presence_g_per_plant,
                .plant_tissue_presence_threshold_g_per_plant = runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
                .canopy_biochemistry_parameters_by_plant = canopy_biochemistry_parameters,
                .plant_phenology = &plant_phenology_state.?,
                .growth_stages = &plant_growth_stage_state.?,
                .emerged_by_plant = plant_phenology_state.?.emerged,
                .root_state = &plant_root_state.?,
                .root_litter_partition = &plant_litter_partition_state.?,
                .root_litter_carbon_ledger = &root_litter_carbon_ledger,
                .shoot_litter_carbon_g_c_by_plant = shoot_harvest_litter_carbon_g_c_by_plant,
                .shoot_litter_nitrogen_g_n_by_plant = shoot_harvest_litter_nitrogen_g_n_by_plant,
                .shoot_litter_phosphorus_g_p_by_plant = shoot_harvest_litter_phosphorus_g_p_by_plant,
                .soil_organic_state = &soil_organic_state,
                .surface_organic_state = &surface_organic_state,
                .surface_nutrient_state = &surface_fire_exchange_state,
                .daily_manure_carbon_input_g_c = daily_organic_fertilizer_carbon_input_g_c,
                .daily_manure_nitrogen_input_g_n = daily_organic_fertilizer_nitrogen_input_g_n,
                .daily_manure_phosphorus_input_g_p = daily_organic_fertilizer_phosphorus_input_g_p,
                .hourly_manure_products_by_plant = hourly_manure_products_by_plant,
                .grid = &state,
                .reseed_population_per_m2_by_plant = reseed_population_per_m2_by_plant,
                .cell_area_m2_by_cell = canopy_cell_area_m2,
                .automatic_harvest_date_by_plant = automatic_harvest_date_by_plant,
            };
            for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
                const active_species = assignments.units[unit_index].species;
                const first = cell * config.plant_populations;
                for (active_species, 0..) |assignment, species| checkpoint_species_names[first + species] = assignment.species_file;
                checkpoint_metadata_cells[cell] = .{
                    .species_names = checkpoint_species_names[first .. first + active_species.len],
                    .species_alive = plant_phenology_state.?.active[first .. first + active_species.len],
                };
            }
        }
    }
    var terrain_radiation_state = try ecosys.terrain_radiation.State.init(allocator, topography, unit_by_cell, canopy_geometry);
    defer terrain_radiation_state.deinit();
    const direct_incidence_count = try std.math.mul(usize, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len);
    const direct_incidence_buffer_count =
        try std.math.mul(usize, state.cell_count, direct_incidence_count);
    const direct_incidence_fraction =
        try allocator.alloc(f64, direct_incidence_buffer_count);
    defer allocator.free(direct_incidence_fraction);
    @memset(direct_incidence_fraction, 0);
    const direct_incidence_per_horizontal_area =
        try allocator.alloc(f64, direct_incidence_buffer_count);
    defer allocator.free(direct_incidence_per_horizontal_area);
    @memset(direct_incidence_per_horizontal_area, 0);
    const direct_scattering_direction =
        try allocator.alloc(
            ecosys.canopy_geometry.ScatteringDirection,
            direct_incidence_buffer_count,
        );
    defer allocator.free(direct_scattering_direction);
    @memset(direct_scattering_direction, .forward);
    var canopy_interception_state: ?ecosys.canopy_interception.State = null;
    defer if (canopy_interception_state) |*interception_state| interception_state.deinit();
    if (canopy_structure_state != null and canopy_optics_state != null) canopy_interception_state = try ecosys.canopy_interception.State.init(allocator, state.cell_count, config.plant_populations, runscript.canopy_layer_count);
    var ground_radiation_state = try ecosys.ground_radiation.State.initMapped(allocator, topography, unit_by_cell, soil_catalog, catalog_index_by_cell);
    defer ground_radiation_state.deinit();
    const initial_snow_temperature_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_snow_temperature_k);
    @memset(initial_snow_temperature_k, site.mean_annual_air_temperature_c + 273.15);
    try snow_transport_state.initializePhysicalState(ground_radiation_state.initial_snow_depth_m, canopy_cell_area_m2, initial_snow_temperature_k, runscript.snow_layer_bottom_depth_m, runscript.initial_snow_density_megagrams_per_m3);
    const snow_depth_m = try allocator.dupe(f64, ground_radiation_state.initial_snow_depth_m);
    defer allocator.free(snow_depth_m);
    @memset(surface_precipitation_state.solid_snow_water_equivalent_m3, 0);
    for (0..state.cell_count) |cell| {
        for (0..snow_transport_state.layer_capacity) |layer| surface_precipitation_state.solid_snow_water_equivalent_m3[cell] += snow_transport_state.solid_snow_water_equivalent_m3[cell * snow_transport_state.layer_capacity + layer];
    }
    try transport_hydrology_state.syncStorage(&state, &snow_transport_state);
    var canopy_exposure_state: ?ecosys.canopy_exposure.State = null;
    defer if (canopy_exposure_state) |*exposure_state| exposure_state.deinit();
    if (canopy_interception_state != null) canopy_exposure_state = try ecosys.canopy_exposure.State.init(allocator, state.cell_count, config.plant_populations);
    var canopy_airflow_state: ?ecosys.canopy_airflow.State = null;
    defer if (canopy_airflow_state) |*airflow_state| airflow_state.deinit();
    var canopy_air_exchange_state: ?ecosys.canopy_air_exchange.State = null;
    defer if (canopy_air_exchange_state) |*air_state| air_state.deinit();
    var standing_dead_air_exchange_state: ?ecosys.canopy_air_exchange.State = null;
    defer if (standing_dead_air_exchange_state) |*air_state| air_state.deinit();
    var canopy_surface_exchange_state: ?ecosys.canopy_surface_exchange.State = null;
    defer if (canopy_surface_exchange_state) |*exchange_state| exchange_state.deinit();
    var standing_dead_surface_exchange_state: ?ecosys.standing_dead_surface_exchange.State = null;
    defer if (standing_dead_surface_exchange_state) |*exchange_state| exchange_state.deinit();
    var canopy_surface_input_workspace: ?ecosys.canopy_surface_exchange.SurfaceInputWorkspace = null;
    defer if (canopy_surface_input_workspace) |*workspace| workspace.deinit();
    if (detailed_canopy_state != null and canopy_exposure_state != null and canopy_precipitation_retention_state != null) {
        canopy_airflow_state = try ecosys.canopy_airflow.State.init(allocator, state.cell_count, config.plant_populations);
        canopy_air_exchange_state = try ecosys.canopy_air_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        standing_dead_air_exchange_state = try ecosys.canopy_air_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        canopy_surface_exchange_state = try ecosys.canopy_surface_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        standing_dead_surface_exchange_state = try ecosys.standing_dead_surface_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        canopy_surface_input_workspace = try ecosys.canopy_surface_exchange.SurfaceInputWorkspace.init(allocator, runtime_plant_count);
    }
    const canopy_atmospheric_vapor_diffusivity_m2_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(canopy_atmospheric_vapor_diffusivity_m2_per_h);
    @memset(canopy_atmospheric_vapor_diffusivity_m2_per_h, 0);
    const soil_atmospheric_gas_conductance_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(soil_atmospheric_gas_conductance_m3_per_h);
    @memset(soil_atmospheric_gas_conductance_m3_per_h, 0);
    const litter_atmospheric_gas_conductance_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(litter_atmospheric_gas_conductance_m3_per_h);
    @memset(litter_atmospheric_gas_conductance_m3_per_h, 0);
    const snow_layer_gas_diffusivity_m2_per_h = try allocator.alloc(f64, snow_transport_state.active.len);
    defer allocator.free(snow_layer_gas_diffusivity_m2_per_h);
    @memset(snow_layer_gas_diffusivity_m2_per_h, 0);
    const canopy_surface_roughness_height_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(canopy_surface_roughness_height_m);
    @memset(canopy_surface_roughness_height_m, runscript.surface_aerodynamic_parameters.soil_roughness_height_m);
    const canopy_available_intercepted_water_m3 = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_available_intercepted_water_m3);
    @memset(canopy_available_intercepted_water_m3, 0);
    const standing_dead_evaporation_m3_per_h = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(standing_dead_evaporation_m3_per_h);
    @memset(standing_dead_evaporation_m3_per_h, 0);
    var surface_energy_state = try ecosys.surface_energy.State.init(allocator, state.cell_count);
    defer surface_energy_state.deinit();
    var ecosystem_energy_ledger_state = try ecosys.ecosystem_energy_ledger.State.init(allocator, state.cell_count);
    defer ecosystem_energy_ledger_state.deinit();
    var plant_energy_publication_state = try ecosys.plant_energy_publication.State.init(
        allocator,
        runtime_plant_count,
    );
    defer plant_energy_publication_state.deinit();
    var plant_water_publication_state = try ecosys.plant_water_publication.State.init(
        allocator,
        state.cell_count,
        config.plant_populations,
    );
    defer plant_water_publication_state.deinit();
    var canopy_water_energy_publication_state =
        try ecosys.canopy_water_energy_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer canopy_water_energy_publication_state.deinit();
    var living_canopy_layer_publication_state =
        try ecosys.living_canopy_layer_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            runscript.canopy_layer_count,
        );
    defer living_canopy_layer_publication_state.deinit();
    var standing_dead_area_publication_state =
        try ecosys.standing_dead_area_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            runscript.canopy_layer_count,
        );
    defer standing_dead_area_publication_state.deinit();
    var root_uptake_ledger_state = try ecosys.root_uptake_ledger.State.init(allocator, state.layer_count);
    defer root_uptake_ledger_state.deinit();
    var root_water_uptake_publication_state =
        try ecosys.root_water_uptake_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_water_uptake_publication_state.deinit();
    var root_gas_content_publication_state =
        try ecosys.root_gas_content_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_gas_content_publication_state.deinit();
    var root_atmosphere_gas_publication_state =
        try ecosys.root_atmosphere_gas_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_atmosphere_gas_publication_state.deinit();
    var root_internal_gas_publication_state =
        try ecosys.root_internal_gas_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_internal_gas_publication_state.deinit();
    var root_soil_gas_publication_state =
        try ecosys.root_soil_gas_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_soil_gas_publication_state.deinit();
    var root_nutrient_uptake_publication_state =
        try ecosys.root_nutrient_uptake_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_nutrient_uptake_publication_state.deinit();
    var root_salt_uptake_publication_state =
        try ecosys.root_salt_uptake_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_salt_uptake_publication_state.deinit();
    var root_soil_ammonia_exchange_publication_state =
        try ecosys.root_soil_ammonia_exchange_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_soil_ammonia_exchange_publication_state.deinit();
    var root_exudate_publication_state =
        try ecosys.root_exudate_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_exudate_publication_state.deinit();
    var root_competition_demand_publication_state =
        try ecosys.root_competition_demand_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_competition_demand_publication_state.deinit();
    var root_nitrogen_fixation_publication_state =
        try ecosys.root_nitrogen_fixation_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_nitrogen_fixation_publication_state.deinit();
    var canopy_carbon_publication_state =
        try ecosys.canopy_carbon_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer canopy_carbon_publication_state.deinit();
    const canopy_net_fixation_g_c_per_h_by_plant =
        try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_net_fixation_g_c_per_h_by_plant);
    @memset(canopy_net_fixation_g_c_per_h_by_plant, 0);
    var root_soil_element_exchange_publication_state =
        try ecosys.root_soil_element_exchange_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer root_soil_element_exchange_publication_state.deinit();
    var plant_balance_publication_state =
        try ecosys.plant_balance_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer plant_balance_publication_state.deinit();
    const plant_balance_carbon_inputs_workspace =
        try allocator.alloc(ecosys.plant_daily_output.CarbonBalanceInputs, runtime_plant_count);
    defer allocator.free(plant_balance_carbon_inputs_workspace);
    @memset(std.mem.sliceAsBytes(plant_balance_carbon_inputs_workspace), 0);
    const plant_balance_nitrogen_inputs_workspace =
        try allocator.alloc(ecosys.plant_daily_output.NutrientBalanceInputs, runtime_plant_count);
    defer allocator.free(plant_balance_nitrogen_inputs_workspace);
    @memset(std.mem.sliceAsBytes(plant_balance_nitrogen_inputs_workspace), 0);
    const plant_balance_phosphorus_inputs_workspace =
        try allocator.alloc(ecosys.plant_daily_output.NutrientBalanceInputs, runtime_plant_count);
    defer allocator.free(plant_balance_phosphorus_inputs_workspace);
    @memset(std.mem.sliceAsBytes(plant_balance_phosphorus_inputs_workspace), 0);
    const plant_balance_carbon_root_by_layer =
        try allocator.alloc(f64, state.soil_layer_capacity);
    defer allocator.free(plant_balance_carbon_root_by_layer);
    @memset(plant_balance_carbon_root_by_layer, 0);
    const plant_balance_carbon_root_length_density_by_layer =
        try allocator.alloc(f64, state.soil_layer_capacity);
    defer allocator.free(plant_balance_carbon_root_length_density_by_layer);
    @memset(plant_balance_carbon_root_length_density_by_layer, 0);
    const plant_balance_nutrient_root_by_layer =
        try allocator.alloc(f64, state.soil_layer_capacity);
    defer allocator.free(plant_balance_nutrient_root_by_layer);
    @memset(plant_balance_nutrient_root_by_layer, 0);
    const root_soil_element_exchange_workspace =
        try allocator.alloc(f64, 9 * runtime_plant_count);
    defer allocator.free(root_soil_element_exchange_workspace);
    @memset(root_soil_element_exchange_workspace, 0);
    const root_soil_carbon_exchange_g_c_per_h_by_plant =
        root_soil_element_exchange_workspace[0..runtime_plant_count];
    const root_soil_nitrogen_exchange_g_n_per_h_by_plant =
        root_soil_element_exchange_workspace[runtime_plant_count .. 2 * runtime_plant_count];
    const root_soil_phosphorus_exchange_g_p_per_h_by_plant =
        root_soil_element_exchange_workspace[2 * runtime_plant_count .. 3 * runtime_plant_count];
    const rse_organic_carbon_by_plant =
        root_soil_element_exchange_workspace[3 * runtime_plant_count .. 4 * runtime_plant_count];
    const rse_organic_nitrogen_by_plant =
        root_soil_element_exchange_workspace[4 * runtime_plant_count .. 5 * runtime_plant_count];
    const rse_organic_phosphorus_by_plant =
        root_soil_element_exchange_workspace[5 * runtime_plant_count .. 6 * runtime_plant_count];
    const rse_phosphate_h2_by_plant =
        root_soil_element_exchange_workspace[6 * runtime_plant_count .. 7 * runtime_plant_count];
    const rse_phosphate_h_by_plant =
        root_soil_element_exchange_workspace[7 * runtime_plant_count .. 8 * runtime_plant_count];
    const rse_canopy_fixation_g_n_per_h_by_plant =
        root_soil_element_exchange_workspace[8 * runtime_plant_count .. 9 * runtime_plant_count];
    var root_gas_withdrawal_publication_state =
        try ecosys.root_gas_withdrawal_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer root_gas_withdrawal_publication_state.deinit();
    // The five UPTAKE/EXTRACT publication owners above advance as one atomic
    // hour, so a partially advanced hour can never be observed by REDIST,
    // NITRO, or the output writers. The workspace is the transaction's
    // rollback scratch and is allocated once here, next to the owners it
    // protects, so `apply` stays allocation-free inside the hourly loop.
    var uptake_coupled_transaction_workspace =
        try ecosys.uptake_coupled_transaction.Workspace.init(allocator, .{
            .cell_count = state.cell_count,
            .species_count = config.plant_populations,
            .soil_layer_capacity = state.soil_layer_capacity,
            .root_domain_capacity = ecosys.plant_root_system.biological_domain_count,
        });
    defer uptake_coupled_transaction_workspace.deinit();
    var canopy_ammonia_publication_state =
        try ecosys.canopy_ammonia_publication.State.init(
            allocator,
            runtime_plant_count,
        );
    defer canopy_ammonia_publication_state.deinit();
    var cell_litter_standing_dead_publication_state =
        try ecosys.grid_cell_litter_standing_dead_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer cell_litter_standing_dead_publication_state.deinit();
    var manure_deposition_publication_state =
        try ecosys.manure_deposition_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer manure_deposition_publication_state.deinit();
    var canopy_fire_publication_state =
        try ecosys.canopy_fire_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer canopy_fire_publication_state.deinit();
    var plant_combustion_publication_state =
        try ecosys.plant_combustion_publication.State.init(
            allocator,
            runtime_plant_count,
        );
    defer plant_combustion_publication_state.deinit();
    var root_combustion_boundary_ledger = std.mem.zeroes(
        ecosys.root_combustion_boundary_publication.Ledger,
    );
    var root_combustion_salt_publication_state =
        try ecosys.root_combustion_salt_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            config.soil_layers,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_combustion_salt_publication_state.deinit();
    const aboveground_litter_publication_workspace =
        try allocator.alloc(f64, 3 * runtime_plant_count);
    defer allocator.free(aboveground_litter_publication_workspace);
    @memset(aboveground_litter_publication_workspace, 0);
    const aboveground_litter_carbon_g_c_per_h_by_plant =
        aboveground_litter_publication_workspace[0..runtime_plant_count];
    const aboveground_litter_nitrogen_g_n_per_h_by_plant =
        aboveground_litter_publication_workspace[runtime_plant_count .. 2 * runtime_plant_count];
    const aboveground_litter_phosphorus_g_p_per_h_by_plant =
        aboveground_litter_publication_workspace[2 * runtime_plant_count .. 3 * runtime_plant_count];
    const canopy_net_radiation_megajoules = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_net_radiation_megajoules);
    const canopy_storage_heat_megajoules = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_storage_heat_megajoules);
    const zero_plant_energy_megajoules = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(zero_plant_energy_megajoules);
    @memset(canopy_net_radiation_megajoules, 0);
    @memset(canopy_storage_heat_megajoules, 0);
    @memset(zero_plant_energy_megajoules, 0);
    const fire_active_this_hour = try allocator.alloc(bool, state.cell_count);
    defer allocator.free(fire_active_this_hour);
    @memset(fire_active_this_hour, false);
    var surface_temperature_solver_state = try ecosys.surface_temperature_solver.State.init(allocator, state.cell_count);
    defer surface_temperature_solver_state.deinit();
    const surface_heat_capacity_megajoules_per_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_heat_capacity_megajoules_per_k);
    @memset(surface_heat_capacity_megajoules_per_k, 0);
    var canopy_energy_state: ?ecosys.canopy_energy.State = null;
    defer if (canopy_energy_state) |*energy_state| energy_state.deinit();
    if (canopy_exposure_state != null) canopy_energy_state = try ecosys.canopy_energy.State.init(allocator, state.cell_count, config.plant_populations);
    var plant_water_balance_state: ?ecosys.plant_water_balance.State = null;
    defer if (plant_water_balance_state) |*balance| balance.deinit();
    var plant_water_workspace: ?ecosys.plant_water_balance.Workspace = null;
    defer if (plant_water_workspace) |*workspace| workspace.deinit();
    if (plant_root_state != null and detailed_canopy_state != null) {
        plant_water_balance_state = try ecosys.plant_water_balance.State.init(allocator, state.cell_count, config.plant_populations, config.soil_layers);
        plant_water_workspace = try ecosys.plant_water_balance.Workspace.init(allocator, state.cell_count, config.plant_populations, config.soil_layers);
        for (0..state.cell_count) |cell|
            plant_water_workspace.?.cell_area_m2[cell] =
                grid_environment.horizontal_cell_width_m[cell] *
                grid_environment.vertical_cell_width_m[cell];
        if (plant_assignments) |assignments| for (plant_unit_by_cell.?, 0..) |unit_index, cell| for (assignments.units[unit_index].species, 0..) |assignment, species| {
            const trait_index = plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile;
            const management_index = plant_management_catalog.find(assignment.management_file) orelse continue;
            const plant = cell * config.plant_populations + species;
            const traits = plant_catalog.entries.items[trait_index].traits;
            const planting = plant_management_catalog.entries.items[management_index].schedule.planting;
            plant_water_workspace.?.leaf_osmotic_potential_at_zero_total_megapascal[plant] = traits.water_relations.osmotic_potential_megapascal;
            plant_water_workspace.?.plant_population_count[plant] = planting.population_per_m2 * plant_water_workspace.?.cell_area_m2[cell];
            plant_water_workspace.?.root_porosity_fraction[plant] = traits.roots.root_porosity_fraction;
            plant_water_workspace.?.root_radial_resistivity_mpa_h_per_m3[plant] = traits.roots.radial_resistivity_mpa_h_per_m3;
            plant_water_workspace.?.root_axial_resistivity_mpa_h_per_m2[plant] = traits.roots.axial_resistivity_mpa_h_per_m2;
            plant_water_workspace.?.primary_root_radius_m[plant] = traits.roots.primary_root_radius_m;
            plant_water_workspace.?.secondary_root_radius_m[plant] = traits.roots.secondary_root_radius_m;
            plant_water_workspace.?.seeding_depth_m[plant] = planting.seed_depth_m;
            plant_water_workspace.?.woody_root_fraction[plant] = 1;
            plant_water_workspace.?.dynamically_woody[plant] = traits.functional_type.root_profile_type != 0 and traits.functional_type.growth_habit > 1;
            plant_water_workspace.?.vascular_growth_habit[plant] = traits.functional_type.growth_habit != 0;
        };
    }
    if (harvest_context) |*context| {
        context.root_woody_fraction_by_plant = if (plant_water_workspace) |*workspace| workspace.woody_root_fraction else null;
        context.carbon_exchange_state = if (canopy_carbon_exchange_state) |*exchange| exchange else null;
    }
    const surface_energy_settings: ecosys.surface_energy.Settings = .{
        .soil_longwave_emissivity = runscript.soil_longwave_emissivity,
        .snow_longwave_emissivity = runscript.snow_longwave_emissivity,
        .snow_full_cover_depth_m = runscript.snow_full_cover_depth_m,
    };
    const nonlinear_solver_options: ecosys.numerics.SolverOptions = .{
        .absolute_tolerance = config.absolute_tolerance,
        .relative_tolerance = config.relative_tolerance,
        .picard_relaxation = config.picard_relaxation,
        .max_iterations = config.max_nonlinear_iterations,
    };
    var current_atmospheric_co2_umol_per_mol = site.atmospheric_co2_umol_mol;
    const current_atmospheric_gas_concentration_g_per_m3 = try allocator.create([ecosys.gas_transport.species_count]f64);
    defer allocator.destroy(current_atmospheric_gas_concentration_g_per_m3);
    current_atmospheric_gas_concentration_g_per_m3.* = surface_gas_parameters.atmospheric_concentration_g_per_m3;
    var active_tile_cells: ?[]const usize = null;
    const runtime_tile_plan: *const ecosys.spatial_grid.TilePlan =
        &tile_plan.?;
    var executed_weather_hours: usize = 0;
    const lateral_contribution_path = try std.fmt.allocPrint(
        allocator,
        "{s}.ecosys-ng-tile-io.{d}",
        // Workspace accepts one safe directory name relative to its parent;
        // input paths may be absolute or contain separators.
        //
        // The suffix makes the directory unique per run. This tree is
        // per-process out-of-core scratch, not a durable artifact: the
        // generation-indexed manifests inside it are double buffered and are
        // meaningful only to the run that wrote them. Naming it from the
        // runscript alone meant two concurrent runs of the same runscript, which
        // is routine on a shared machine and during parallel validation, shared
        // one directory and read each other's manifests. That surfaced as a
        // spurious `LateralContributionFileGenerationMismatch` partway through
        // an otherwise healthy run, so a run could silently corrupt another
        // run's state. The name is made unique rather than the generation check
        // relaxed, because that check is correct. The clock comes from the
        // injected `Io` so no new platform dependency is introduced.
        .{
            std.fs.path.basename(args[1]),
            std.Io.Clock.now(.boot, init.io).nanoseconds,
        },
    );
    defer allocator.free(lateral_contribution_path);
    // Use the OS temp directory so that cloud sync tools (e.g. OneDrive) do not
    // hold file locks on tile manifests during atomic rename operations. The
    // workspace only needs the parent directory handle during init; it keeps its
    // own handles to the subdirectories it creates.
    const tile_io_parent_temp_dir: ?std.Io.Dir = try openOsTempDir(allocator, init.io);
    defer if (tile_io_parent_temp_dir) |d| d.close(init.io);
    // A unique name would otherwise accumulate one tree per invocation.
    // Registered before the workspace so that, defers being LIFO, removal runs
    // after the workspace has released its directory handles. A failure to
    // remove is not fatal: the run's scientific result is already complete and
    // the OS reclaims its own temp space.
    defer if (tile_io_parent_temp_dir) |parent|
        parent.deleteTree(init.io, lateral_contribution_path) catch {};
    var lateral_contribution_workspace =
        try ecosys.hourly_lateral_contribution_io.Workspace.init(
            init.io,
            tile_io_parent_temp_dir orelse std.Io.Dir.cwd(),
            lateral_contribution_path,
        );
    defer lateral_contribution_workspace.deinit();
    var hourly_science_context = .{
        .allocator = allocator,
        .executor = executor,
        .active_tile_cells = &active_tile_cells,
        .tile_plan = runtime_tile_plan,
        .executed_weather_hours = &executed_weather_hours,
        .restoring_checkpoint = scene_options[0].resume_from_checkpoint,
        .lateral_contribution_workspace = &lateral_contribution_workspace,
        .landscape_boundary_ledger = &landscape_mass_balance_state.boundary_ledger,
        .landscape_soil_mass_megagrams_scratch = landscape_soil_mass_megagrams_scratch,
        .grid = &state,
        .lateral_connection_mode_by_cell = lateral_connection_mode_by_cell,
        .site_by_cell = site_by_cell,
        .mean_annual_temperature_c_by_cell = mean_annual_temperature_c_by_cell,
        .mean_annual_temperature_k_by_cell = mean_annual_temperature_k_by_cell,
        .geothermal_enabled_by_cell = geothermal_enabled_by_cell,
        .surface_runoff_boundary_fraction_by_direction = surface_runoff_boundary_fraction_by_direction,
        .atmosphere = &atmospheric_state,
        .hourly_adjusted_shortwave_megajoules_per_m2 = hourly_adjusted_shortwave_megajoules_per_m2,
        .hourly_weather_reference_height_m = hourly_weather_reference_height_m,
        .hourly_extraterrestrial_shortwave_megajoules_per_m2 = hourly_extraterrestrial_shortwave_megajoules_per_m2,
        .hourly_solar_angle_sine = hourly_solar_angle_sine,
        .hourly_solar_azimuth_radians = hourly_solar_azimuth_radians,
        .irrigation_water_depth_m = irrigation_water_depth_m,
        .irrigation_dissolved_mass_g_per_m2 = irrigation_dissolved_mass_g_per_m2,
        .irrigation_hydrogen_mol_per_m2 = irrigation_hydrogen_mol_per_m2,
        .subsurface_irrigation_water_m3 = irrigation_loads.subsurface_water_m3,
        .irrigation_loads = &irrigation_loads,
        .canopy_radiation = &canopy_radiation_state,
        .canopy_optics = &canopy_optics_state,
        .canopy_geometry = &canopy_geometry,
        .canopy_irradiance_interception_geometry = &canopy_irradiance_interception_geometry_state,
        .canopy_cell_area_m2 = canopy_cell_area_m2,
        .surface_aerodynamics = &surface_aerodynamic_state,
        .ground_air = &ground_air_state,
        .surface_total_canopy_area_m2 = surface_total_canopy_area_m2,
        .surface_canopy_height_m = surface_canopy_height_m,
        .ground_air_vapor_pressure_kpa = ground_air_vapor_pressure_kpa,
        .atmospheric_vapor_fraction = atmospheric_vapor_fraction,
        .ground_air_canopy_resistance_h_per_m = ground_air_canopy_resistance_h_per_m,
        .ground_air_sensible_source_megajoules_per_h = ground_air_sensible_source_megajoules_per_h,
        .delayed_live_canopy_combustion_heat_megajoules = delayed_live_canopy_combustion_heat_megajoules,
        .delayed_standing_dead_combustion_heat_megajoules = delayed_standing_dead_combustion_heat_megajoules,
        .delayed_subsurface_combustion_heat_megajoules = delayed_subsurface_combustion_heat_megajoules,
        .delayed_surface_combustion_heat_megajoules = delayed_surface_combustion_heat_megajoules,
        .surface_combustion_heat_megajoules_per_m2 = surface_combustion_heat_megajoules_per_m2,
        .ground_air_vapor_source_m3_per_h = ground_air_vapor_source_m3_per_h,
        .ground_air_surface_sensible_conductance_megajoules_per_h_k = ground_air_surface_sensible_conductance_megajoules_per_h_k,
        .ground_air_surface_vapor_conductance_m3_per_h = ground_air_surface_vapor_conductance_m3_per_h,
        .ground_air_surface_vapor_fraction = ground_air_surface_vapor_fraction,
        .ground_surface_evaporation_m3_per_h = ground_surface_evaporation_m3_per_h,
        .ground_surface_condensation_m3_per_h = ground_surface_condensation_m3_per_h,
        .ground_surface_litter_water_change_m3 = ground_surface_litter_water_change_m3,
        .ground_surface_topsoil_water_change_m3 = ground_surface_topsoil_water_change_m3,
        .direct_incidence_fraction = direct_incidence_fraction,
        .direct_incidence_per_horizontal_area = direct_incidence_per_horizontal_area,
        .direct_scattering_direction = direct_scattering_direction,
        .canopy_interception = &canopy_interception_state,
        .canopy_structure = &canopy_structure_state,
        .terrain_radiation = &terrain_radiation_state,
        .ground_radiation = &ground_radiation_state,
        .snow_depth_m = snow_depth_m,
        .canopy_exposure = &canopy_exposure_state,
        .surface_energy = &surface_energy_state,
        .surface_energy_settings = surface_energy_settings,
        .surface_temperature = &surface_temperature_solver_state,
        .surface_heat_capacity_megajoules_per_k = surface_heat_capacity_megajoules_per_k,
        .nonlinear_solver_options = nonlinear_solver_options,
        .soil_thermal = &soil_thermal_state,
        .soil_thermal_context = &soil_thermal_context,
        .soil_hourly_workspace = &soil_hourly_workspace,
        .soil_heat_solver_workspace = &soil_heat_solver_workspace,
        .soil_solver_properties = &soil_solver_property_state,
        .soil_geometry = &soil_geometry_state,
        .fertilizer_band = &fertilizer_band_state,
        .surface_pond_domain_workspace = &surface_pond_domain_workspace,
        .soil_profile_relayering_workspace = &soil_profile_relayering_workspace,
        .terrain_hydrology = &terrain_hydrology_state,
        .surface_runoff = &surface_runoff_state,
        .soil_numerical_scales = soil_numerical_scales_state,
        .initial_balance_ledger = &initial_balance_ledger_state,
        .surface_inorganic_nitrogen_export_g_n_per_h = surface_inorganic_nitrogen_export_g_n_per_h,
        .surface_inorganic_phosphorus_export_g_p_per_h = surface_inorganic_phosphorus_export_g_p_per_h,
        .surface_organic_carbon_export_g_c_per_h = surface_organic_carbon_export_g_c_per_h,
        .surface_inorganic_carbon_export_g_c_per_h = surface_inorganic_carbon_export_g_c_per_h,
        .surface_organic_nitrogen_export_g_n_per_h = surface_organic_nitrogen_export_g_n_per_h,
        .surface_organic_phosphorus_export_g_p_per_h = surface_organic_phosphorus_export_g_p_per_h,
        .surface_erosion = &surface_erosion_state,
        .surface_soil_mass_at_erosion_start_megagrams = surface_soil_mass_at_erosion_start_megagrams,
        .net_sediment_megagrams_per_h = net_sediment_megagrams_per_h,
        .eroded_mineral_state = &eroded_mineral_state,
        .soil_transport_faces = &soil_transport_faces,
        .soil_face_geometry = &soil_face_geometry_state,
        .soil_solute_face_parameters = &soil_solute_face_parameter_state,
        .micropore_solute_face_flux_mol = micropore_solute_face_flux_mol,
        .macropore_solute_face_flux_mol = macropore_solute_face_flux_mol,
        .soil_organic_face_parameters = &soil_organic_face_parameter_state,
        .soil_organic_transport = &soil_organic_transport_state,
        .soil_organic_recharge_concentration_g_per_m3 = soil_organic_recharge_concentration_g_per_m3,
        .soil_dissolved_gas_face_parameters = &soil_dissolved_gas_face_parameter_state,
        .soil_dissolved_gas_transport = &soil_dissolved_gas_transport_state,
        .soil_dissolved_gas_recharge_concentration_g_per_m3 = soil_dissolved_gas_recharge_concentration_g_per_m3,
        .mineral_nitrogen_transport = &mineral_nitrogen_transport_state,
        .mineral_nitrogen_face_parameters = &mineral_nitrogen_face_parameters,
        .mineral_nitrogen_zone_fractions = mineral_nitrogen_zone_fractions,
        .soil_boundary_topology = &soil_boundary_topology_state,
        .transport_hydrology = &transport_hydrology_state,
        .snow_transport = &snow_transport_state,
        .snow_atmospheric_input_g = snow_atmospheric_input_g,
        .atmospheric_solute_input_ledger = &atmospheric_solute_input_ledger_state,
        .snow_surface_partitions = snow_surface_partitions,
        .snow_surface_discharge = snow_surface_discharge,
        .direct_surface_solute_input = direct_surface_solute_input,
        .snowpack_internal_solute_flux_by_layer = snowpack_internal_solute_flux_by_layer,
        .snowpack_internal_solute_flux_workspace = snowpack_internal_solute_flux_workspace,
        .snowpack_internal_salt_flux_mol_by_layer_species = snowpack_internal_salt_flux_mol_by_layer_species,
        .snowpack_internal_salt_flux_workspace_by_layer_species = snowpack_internal_salt_flux_workspace_by_layer_species,
        .micropore_solute_state = &micropore_solute_state,
        .macropore_solute_state = &macropore_solute_state,
        .soil_recharge_concentration_mol_per_m3 = soil_recharge_concentration_mol_per_m3,
        .soil_solute_boundary_net_flux_mol = soil_solute_boundary_net_flux_mol,
        .gas_transport = &gas_transport_state,
        .soil_gas_transport = &soil_gas_transport_state,
        .surface_litter_gas_transport = &surface_litter_gas_transport_state,
        .soil_reactive_nitrogen = &soil_reactive_nitrogen_state,
        .soil_microbial_phosphorus = &soil_microbial_phosphorus_state,
        .soil_microbial_turnover = &soil_microbial_turnover_state,
        .soil_litter_colonization = &soil_litter_colonization_state,
        .soil_organic_sorption = &soil_organic_sorption_state,
        .soil_organic_decomposition = &soil_organic_decomposition_state,
        .soil_organic_priming = &soil_organic_priming_state,
        .soil_respiration_products = &soil_respiration_products_state,
        .daily_heterotrophic_respiration = &daily_heterotrophic_respiration,
        .soil_autotrophic_carbon = &soil_autotrophic_carbon_state,
        .soil_microbial_layer_mixing = &soil_microbial_layer_mixing_state,
        .soil_biogeochemical_gas_fluxes = &soil_biogeochemical_gas_flux_state,
        .soil_methane = &soil_methane_state,
        .soil_microbial_oxygen = &soil_microbial_oxygen_state,
        .soil_oxygen_staging = &soil_oxygen_staging_state,
        .soil_field_capacity_fraction = soil_field_capacity_fraction,
        .soil_nitrogen_flux_workspace = &soil_nitrogen_flux_workspace,
        .soil_nitrifier_environment = &soil_nitrifier_environment_state,
        .soil_microbial = &soil_microbial_state,
        .soil_redox_satisfaction_fraction = soil_redox_satisfaction_fraction,
        .soil_nitrogen_parameters = &soil_nitrogen_parameters,
        .plant_available_nutrients = &plant_available_nutrient_state,
        .iteration_limits = iteration_limits,
        .config = config,
        .runscript = &runscript,
        .site = &site,
        .horizontal_cell_width_m = grid_environment.horizontal_cell_width_m,
        .vertical_cell_width_m = grid_environment.vertical_cell_width_m,
        .current_atmospheric_co2_umol_per_mol = &current_atmospheric_co2_umol_per_mol,
        .current_atmospheric_gas_concentration_g_per_m3 = current_atmospheric_gas_concentration_g_per_m3,
        .canopy_energy = &canopy_energy_state,
        .canopy_airflow = &canopy_airflow_state,
        .canopy_air_exchange = &canopy_air_exchange_state,
        .standing_dead_air_exchange = &standing_dead_air_exchange_state,
        .canopy_surface_exchange = &canopy_surface_exchange_state,
        .standing_dead_surface_exchange = &standing_dead_surface_exchange_state,
        .canopy_surface_input_workspace = &canopy_surface_input_workspace,
        .canopy_atmospheric_vapor_diffusivity_m2_per_h = canopy_atmospheric_vapor_diffusivity_m2_per_h,
        .soil_atmospheric_gas_conductance_m3_per_h = soil_atmospheric_gas_conductance_m3_per_h,
        .litter_atmospheric_gas_conductance_m3_per_h = litter_atmospheric_gas_conductance_m3_per_h,
        .snow_layer_gas_diffusivity_m2_per_h = snow_layer_gas_diffusivity_m2_per_h,
        .canopy_surface_roughness_height_m = canopy_surface_roughness_height_m,
        .canopy_available_intercepted_water_m3 = canopy_available_intercepted_water_m3,
        .standing_dead_evaporation_m3_per_h = standing_dead_evaporation_m3_per_h,
        .plants = &plant_state,
        .plant_water_balance = &plant_water_balance_state,
        .plant_water_workspace = &plant_water_workspace,
        .plant_roots = &plant_root_state,
        .plant_root_nutrient_workspace = &plant_root_nutrient_workspace,
        .plant_root_salt_workspace = &plant_root_salt_workspace,
        .plant_root_exudation_workspace = &plant_root_exudation_workspace,
        .plant_root_metabolism_workspace = &plant_root_metabolism_workspace,
        .plant_litter_partition = &plant_litter_partition_state,
        .plant_storage_remobilization_workspace = &plant_storage_remobilization_workspace,
        .plant_harvest = if (harvest_context) |*value| value else null,
        .root_nutrient_traits = root_nutrient_traits,
        .root_nutrient_feedback_enabled = root_nutrient_feedback_enabled,
        .root_metabolism_plant_parameters = root_metabolism_plant_parameters,
        .root_biological_domain_count_by_plant = root_biological_domain_count_by_plant,
        .plant_phenology = &plant_phenology_state,
        .plant_topology_controls = &plant_topology_controls,
        .plant_reproduction_controls = &plant_reproduction_controls,
        .canopy_layer_controls = &canopy_layer_controls,
        .branch_development = &branch_development_state,
        .plant_growth_stages = &plant_growth_stage_state,
        .plant_dormancy = &plant_dormancy_state,
        .phenology_root_oxygen_fraction = phenology_root_oxygen_fraction,
        .development_species_parameters = development_species_parameters,
        .development_dormancy_parameters = development_dormancy_parameters,
        .development_planting_day_of_year = development_planting_day_of_year,
        .development_planting_year = development_planting_year,
        .development_canopy_height_m = development_canopy_height_m,
        .development_surface_water_potential_megapascal = development_surface_water_potential_megapascal,
        .development_seed_layer_water_potential_megapascal = development_seed_layer_water_potential_megapascal,
        .development_emerged = plant_phenology_state.?.emerged,
        .root_gas_parameters = runscript.root_gas_parameters,
        .detailed_canopy = &detailed_canopy_state,
        .canopy_carbon_exchange = &canopy_carbon_exchange_state,
        .canopy_biochemistry_parameters = canopy_biochemistry_parameters,
        .shoot_growth_plant_parameters = shoot_growth_plant_parameters,
        .canopy_layer_distribution = &canopy_layer_distribution_state,
        .canopy_precipitation_retention = &canopy_precipitation_retention_state,
        .surface_precipitation = &surface_precipitation_state,
        .surface_pond_transition = &surface_pond_transition_state,
        .surface_pond_minimum_heat_capacity_megajoules_per_k = surface_pond_minimum_heat_capacity_megajoules_per_k,
        .surface_litter_geometry = &surface_litter_geometry_state,
        .surface_litter_water_environment = &surface_litter_water_environment_state,
        .surface_microbial_environment = &surface_microbial_environment_state,
        .surface_microbial_respiration = &surface_microbial_respiration_state,
        .surface_microbial_oxygen = &surface_microbial_oxygen_state,
        .surface_microbial_maintenance = &surface_microbial_maintenance_state,
        .surface_nonsymbiotic_nitrogen_fixation = &surface_nonsymbiotic_nitrogen_fixation_state,
        .surface_microbial_substrate_uptake = &surface_microbial_substrate_uptake_state,
        .surface_denitrification = &surface_denitrification_state,
        .surface_microbial_assimilation = &surface_microbial_assimilation_state,
        .surface_microbial_mineral_exchange = &surface_microbial_mineral_exchange_state,
        .surface_topsoil_mineral_exchange = &surface_topsoil_mineral_exchange_state,
        .surface_microbial_turnover = &surface_microbial_turnover_state,
        .surface_organic_priming = &surface_organic_priming_state,
        .surface_organic_decomposition = &surface_organic_decomposition_state,
        .surface_organic_sorption = &surface_organic_sorption_state,
        .surface_litter_colonization = &surface_litter_colonization_state,
        .surface_humification_fraction = surface_humification_fraction,
        .soil_chemistry = &initial_chemistry_state,
        .soil_chemistry_layer_parameters = soil_chemistry_layer_parameters,
        .soil_chemistry_solver_workspace = &soil_chemistry_solver_workspace,
        .eroded_chemistry_workspace = &eroded_chemistry_workspace,
        .soil_organic = &soil_organic_state,
        .soil_organic_carbon_at_hour_start_g_c = soil_organic_carbon_at_hour_start_g_c,
        .soil_organic_carbon_change_g_c_per_h = soil_organic_carbon_change_g_c_per_h,
        .topsoil_humus_partition = topsoil_humus_partition,
        .surface_litter_fertilizer = &surface_litter_fertilizer_state,
        .soil_fertilizer_inventory = &soil_fertilizer_inventory,
        .mineral_fertilizer_inventory = &mineral_fertilizer_inventory,
        .eroded_fertilizer_workspace = &eroded_fertilizer_workspace,
        .surface_litter_fertilizer_diagnostics = &surface_litter_fertilizer_diagnostics,
        .surface_field_capacity_potential_megapascal = surface_field_capacity_potential_megapascal,
        .surface_wilting_point_potential_megapascal = surface_wilting_point_potential_megapascal,
        .surface_litter_ice_m3 = surface_litter_ice_m3,
        .surface_charcoal_carbon_g_c = surface_charcoal_carbon_g_c,
        .litter_gas_transport = &litter_gas_transport_state,
        .surface_gas_parameters = &surface_gas_parameters,
        .surface_litter_chemistry = &surface_litter_chemistry_state,
        .surface_litter_chemistry_diagnostics = &surface_litter_chemistry_diagnostics,
        .surface_litter_cation_selectivity = surface_litter_cation_selectivity,
        .surface_organic = &surface_organic_state,
        .shoot_senescence_products_by_plant = shoot_senescence_products_by_plant,
        .seasonal_turnover_event_by_plant = seasonal_turnover_event_by_plant,
        .root_litter_products_by_plant = root_litter_products_by_plant,
        .root_litter_carbon_ledger = &root_litter_carbon_ledger,
        .eroded_organic_workspace = &eroded_organic_workspace,
        .organic_parameters = &organic_parameters,
        .chemistry_reaction_parameters = &chemistry_reaction_parameters,
        .root_soil_ammonia_exchange_publication_state = &root_soil_ammonia_exchange_publication_state,
        .root_nutrient_uptake_publication_state = &root_nutrient_uptake_publication_state,
    };
    var climate_state: ecosys.climate_change.State = .{};
    var executed_scene_passes: usize = 0;
    var resume_position: ?ecosys.checkpoint_manifest.SimulationInstant = null;
    var first_pass_preview = try ecosys.simulation_timeline.PassIterator.init(runscript.scenarios, runscript.execution_repeat_count);
    const first_requested_pass = first_pass_preview.next() orelse return error.EmptySimulationTimeline;
    hourly_science_context.restoring_checkpoint =
        scene_options[first_requested_pass.scene_index].resume_from_checkpoint;
    // Legacy multi-scene runs mark later scenes as continuation scenes. They
    // already inherit live state in this process; only a resume flag on the
    // first executed scene requests an external bundle restore.
    if (scene_options[first_requested_pass.scene_index].resume_from_checkpoint) {
        if (plant_root_state == null or detailed_canopy_state == null or canopy_precipitation_retention_state == null or plant_phenology_state == null or plant_growth_stage_state == null or plant_dormancy_state == null or branch_development_state == null) return error.CheckpointResumeRequiresCompletePlantState;
        var maximum_species_name_bytes: usize = 1;
        for (checkpoint_species_names) |name| maximum_species_name_bytes = @max(maximum_species_name_bytes, name.len);
        const canopy = &detailed_canopy_state.?;
        const maximum_branches = canopy.plant_branch_offsets[canopy.plant_branch_offsets.len - 1];
        const maximum_nodes = canopy.branch_node_offsets[canopy.branch_node_offsets.len - 1];
        const maximum_samples = canopy.node_sample_offsets[canopy.node_sample_offsets.len - 1];
        var restored = try ecosys.checkpoint_bundle_reader.read(
            allocator,
            init.io,
            std.Io.Dir.cwd(),
            config,
            .{
                .manifest = .{ .maximum_columns = config.lon_count, .maximum_rows = config.lat_count, .maximum_soil_layers = config.soil_layers, .maximum_snow_layers = snow_transport_state.layer_capacity, .maximum_plant_species_per_cell = config.plant_populations, .maximum_root_axes_per_plant = runscript.root_axes_per_plant },
                .plant_metadata = .{ .maximum_cells = state.cell_count, .maximum_species_per_cell = config.plant_populations, .maximum_species_name_bytes = maximum_species_name_bytes },
                .plant_development = .{ .maximum_cells = state.cell_count, .maximum_species = config.plant_populations, .maximum_branches = plant_growth_stage_state.?.branches.len },
                .plant_roots = .{ .maximum_plants = plant_root_state.?.plant_count, .maximum_soil_layers = config.soil_layers, .maximum_root_axes = runscript.root_axes_per_plant },
                .plant_canopy = .{ .maximum_cells = state.cell_count, .maximum_species = config.plant_populations, .maximum_branches = maximum_branches, .maximum_nodes = maximum_nodes, .maximum_samples = maximum_samples, .maximum_layers = runscript.canopy_layer_count, .maximum_inclinations = canopy_geometry.leaf_inclination_sine.len, .maximum_azimuths = canopy_geometry.leaf_azimuth_radians.len },
                .soil_biogeochemistry = .{ .maximum_cells = state.cell_count, .maximum_layers = config.soil_layers, .maximum_substrates = runscript.microbial_substrate_count, .maximum_populations = runscript.microbial_population_count },
                .soil_organic = .{ .maximum_profile_layers = state.layer_count, .maximum_surface_cells = state.cell_count },
                .transport = .{ .maximum_transport_cells = state.layer_count, .maximum_solute_species = micropore_solute_state.species_count, .maximum_snow_cells = state.cell_count, .maximum_snow_layers = snow_transport_state.layer_capacity },
                .soil_geometry = .{ .maximum_columns = config.lon_count, .maximum_rows = config.lat_count, .maximum_soil_layers = config.soil_layers, .maximum_snow_layers = snow_transport_state.layer_capacity, .maximum_plants = runtime_plant_count },
            },
            .{ .manifest_buffer_bytes = 16 * 1024, .section_read_buffer_bytes = 64 * 1024, .section_verify_buffer_bytes = 64 * 1024 },
        );
        const restored_generation = restored.manifest.generation;
        resume_position = restored.manifest.instant;
        try ecosys.checkpoint_bundle_reader.swapIntoLive(&restored, .{
            .grid = &state,
            .plants = &plant_state,
            .plant_development = .{ .phenology = &plant_phenology_state.?, .growth = &plant_growth_stage_state.?, .dormancy = &plant_dormancy_state.?, .branch_development = &branch_development_state.? },
            .plant_roots = &plant_root_state.?,
            .plant_canopy = .{ .canopy = &detailed_canopy_state.?, .retention = &canopy_precipitation_retention_state.?, .layer_distribution = &canopy_layer_distribution_state.? },
            .soil_biogeochemistry = .{ .microbial = &soil_microbial_state, .chemistry = &initial_chemistry_state, .available_nutrients = &plant_available_nutrient_state, .fertilizer = &soil_fertilizer_inventory, .mineral_fertilizer = &mineral_fertilizer_inventory, .fertilizer_band = &fertilizer_band_state, .reactive_nitrogen = &soil_reactive_nitrogen_state, .microbial_phosphorus = &soil_microbial_phosphorus_state },
            .soil_organic_matter = .{ .profile = &soil_organic_state, .surface = &surface_organic_state, .litter_chemistry = &surface_litter_chemistry_state, .litter_fertilizer = &surface_litter_fertilizer_state, .surface_respiration = &surface_microbial_respiration_state, .surface_denitrification = &surface_denitrification_state, .surface_fire_exchange = &surface_fire_exchange_state, .litter_salt_ingress = &plant_litter_salt_ingress_state },
            .transport = .{ .micropore = &micropore_solute_state, .macropore = &macropore_solute_state, .mineral_nitrogen = &mineral_nitrogen_transport_state, .organic = &soil_organic_transport_state, .gas = &gas_transport_state, .litter_gas = &litter_gas_transport_state, .snow = &snow_transport_state, .surface = &surface_transport_state },
            .soil_geometry_and_hydrology = .{ .geometry = &soil_geometry_state, .hydrology = &transport_hydrology_state, .surface = &surface_precipitation_state, .erosion = &surface_erosion_state, .climate = &climate_state, .eroded_minerals = &eroded_mineral_state, .runtime = .{ .soil_properties = &soil_solver_property_state, .soil_thermal = &soil_thermal_state }, .surface_boundary = .{ .ground_air = &ground_air_state, .surface_aerodynamics = &surface_aerodynamic_state }, .surface_litter_geometry = &surface_litter_geometry_state, .surface_litter_ice_m3 = surface_litter_ice_m3, .delayed_live_canopy_combustion_heat_megajoules = delayed_live_canopy_combustion_heat_megajoules, .delayed_standing_dead_combustion_heat_megajoules = delayed_standing_dead_combustion_heat_megajoules, .delayed_subsurface_combustion_heat_megajoules = delayed_subsurface_combustion_heat_megajoules, .delayed_surface_combustion_heat_megajoules = delayed_surface_combustion_heat_megajoules },
            .landscape_mass_balance = &landscape_mass_balance_state,
        });
        // The owner swap replaces these slices. Rebind consumers before the
        // restored bundle releases the former live allocation.
        hourly_science_context.development_emerged = plant_phenology_state.?.emerged;
        hourly_science_context.restoring_checkpoint = true;
        harvest_context.?.emerged_by_plant = plant_phenology_state.?.emerged;
        if (plant_assignments) |assignments| for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
            const active_species = assignments.units[unit_index].species.len;
            const first = cell * config.plant_populations;
            checkpoint_metadata_cells[cell].species_alive = plant_phenology_state.?.active[first .. first + active_species];
        };
        restored.deinit();
        // The checkpointed thermal owner is authoritative. Its accepted
        // layer thickness and volume can intentionally differ from the soil
        // geometry carrier after freeze/thaw and remapping; rebuilding either
        // value here changes the next surface heat-storage equation.
        // The same rule applies to the checkpointed runtime soil-property
        // owner. Erosion texture, organic-carbon concentration, exchange
        // capacities, and rainfall conductivity damage are accepted mutable
        // state, not restart-time derivations.
        try ecosys.transport_hydrology.refreshSoilFacesFromHydrology(
            &soil_transport_faces,
            &transport_hydrology_state,
        );
        try soil_face_geometry_state.refreshMapped(
            &state,
            &soil_transport_faces,
            soil_solver_property_state.layer_thickness_m,
            grid_environment.horizontal_cell_width_m,
            grid_environment.vertical_cell_width_m,
        );
        executed_weather_hours = std.math.cast(usize, restored_generation) orelse return error.CheckpointGenerationExceedsRuntimeCounter;
    }
    var pass_iterator = try ecosys.simulation_timeline.PassIterator.init(runscript.scenarios, runscript.execution_repeat_count);
    while (pass_iterator.next()) |pass| {
        {
            const scene = runscript.scenes[pass.scene_index];
            const active_options = scene_options[pass.scene_index];
            try visualization_streams.configureScene(
                active_options.visualization_enabled,
                active_options.visualization_start_year,
                active_options.visualization_end_year,
            );
            if (resume_position) |position| {
                const before_target = pass.execution_iteration < position.execution_iteration or
                    (pass.execution_iteration == position.execution_iteration and pass.scenario_index < position.scenario_index) or
                    (pass.execution_iteration == position.execution_iteration and pass.scenario_index == position.scenario_index and pass.scenario_iteration < position.scenario_iteration) or
                    (pass.execution_iteration == position.execution_iteration and pass.scenario_index == position.scenario_index and pass.scenario_iteration == position.scenario_iteration and pass.scene_index < position.scene_index);
                if (before_target) continue;
                if (pass.execution_iteration != position.execution_iteration or pass.scenario_index != position.scenario_index or pass.scenario_iteration != position.scenario_iteration or pass.scene_index != position.scene_index) return error.CheckpointTimelinePositionNotInRunscript;
            }
            const active_iteration_limits = try ecosys.iteration_control.Limits.fromSceneOptions(active_options);
            hourly_science_context.iteration_limits = active_iteration_limits;
            hourly_science_context.nonlinear_solver_options.max_iterations = active_iteration_limits.water_heat_solute_max_iterations;

            const weather_assignment = scene_weather_assignments[pass.scene_index];
            const weather_streams = try allocator.alloc(
                *ecosys.weather.WeatherStream,
                weather_assignment.unique_file_names.len,
            );
            defer allocator.free(weather_streams);
            var opened_weather_stream_count: usize = 0;
            defer for (weather_streams[0..opened_weather_stream_count]) |stream| stream.deinit();
            const hourly_weather_streams = try allocator.alloc(
                ecosys.hourly_weather_stream.Stream,
                weather_assignment.unique_file_names.len,
            );
            defer allocator.free(hourly_weather_streams);
            for (
                weather_assignment.unique_file_names,
                weather_assignment.forcing_geometry_by_stream,
                0..,
            ) |weather_name, geometry, stream_index| {
                const weather_path = try resolveInputPath(
                    allocator,
                    init.io,
                    runscript_directory,
                    weather_name,
                );
                defer allocator.free(weather_path);
                weather_streams[stream_index] = ecosys.weather.WeatherStream.init(
                    allocator,
                    init.io,
                    weather_path,
                    64 * 1024,
                ) catch |err| {
                    std.log.err("weather stream open failed: execution={d} scenario={d} scenario_repeat={d} scene={d} file='{s}' error={s}", .{ pass.execution_iteration + 1, pass.scenario_index + 1, pass.scenario_iteration + 1, pass.scene_index + 1, weather_name, @errorName(err) });
                    return err;
                };
                hourly_science_context.restoring_checkpoint = false;
                opened_weather_stream_count += 1;
                hourly_weather_streams[stream_index] = try ecosys.hourly_weather_stream.Stream.init(
                    weather_streams[stream_index],
                    geometry.altitude_m,
                    geometry.latitude_degrees_north,
                    geometry.phytotron,
                    .{
                        .snowfall_temperature_threshold_c = active_options.snowfall_temperature_threshold_c,
                        .minimum_snowfall_water_equivalent_m = active_options.minimum_snowfall_water_equivalent_m,
                    },
                );
            }
            const weather_stream = weather_streams[0];
            const fertilizer_solar_noon_hour_by_cell = try allocator.alloc(u8, state.cell_count);
            defer allocator.free(fertilizer_solar_noon_hour_by_cell);
            for (fertilizer_solar_noon_hour_by_cell, 0..) |*solar_noon, cell| {
                const stream_index = weather_assignment.stream_index_by_cell[cell];
                if (stream_index >= weather_streams.len)
                    return error.WeatherGridStreamIndexOutOfRange;
                const value = weather_streams[stream_index].header.solar_noon_hour;
                if (!std.math.isFinite(value) or value < 0 or value >= 25)
                    return error.InvalidHourlyFertilizerSchedule;
                solar_noon.* = @intFromFloat(@trunc(value));
            }
            const weather_hour_by_stream = try allocator.alloc(
                ecosys.hourly_weather_stream.Observation,
                hourly_weather_streams.len,
            );
            defer allocator.free(weather_hour_by_stream);
            var previous_weather_timestamp: ?ecosys.weather.Timestamp = null;
            var resolved_development_year: ?u16 = null;
            var dormant_seed_branches: std.ArrayList(ecosys.plant_harvest_runtime.SourceOrderDormantSeedBranch) = .empty;
            defer dormant_seed_branches.deinit(allocator);
            var scene_weather_hours: usize = 0;
            if (resume_position) |position| {
                const completed_scene_hours = std.math.cast(usize, position.completed_scene_hours) orelse return error.CheckpointSceneHourExceedsRuntimeCounter;
                for (0..completed_scene_hours) |_| {
                    for (hourly_weather_streams, 0..) |*stream, stream_index| {
                        const skipped = try stream.next() orelse return error.CheckpointSceneHourExceedsWeatherStream;
                        if (stream_index == 0) {
                            previous_weather_timestamp = skipped.timestamp;
                        } else if (!sameWeatherTimestamp(previous_weather_timestamp.?, skipped.timestamp)) {
                            return error.WeatherGridTimestampMismatch;
                        }
                    }
                }
                const checkpoint_timestamp = previous_weather_timestamp orelse return error.CheckpointHasNoCompletedWeatherHour;
                const checkpoint_year = checkpoint_timestamp.year orelse active_options.start_date.year;
                const checkpoint_day = checkpoint_timestamp.day_of_year orelse return error.CheckpointWeatherRequiresDayOfYear;
                const checkpoint_hour = if (checkpoint_timestamp.hour == 24) 23 else checkpoint_timestamp.hour;
                if (checkpoint_year != position.year or checkpoint_day != position.day_of_year or checkpoint_hour != position.hour) return error.CheckpointWeatherTimestampMismatch;
                scene_weather_hours = completed_scene_hours;
                resume_position = null;
            }
            // EXEC IBEGIN/ISTART reset: establish a new seven-domain
            // reference only at a fresh scene boundary. A resumed scene keeps
            // the checkpointed monitor and cumulative boundary history.
            if (scene_weather_hours == 0) {
                // Sync gas_transport dissolved CO2 and micropore solutes from
                // STARTE-equilibrated chemistry before capturing the EXEC baseline,
                // so the baseline counts the same carbonate inventory that every
                // subsequent hour counts after exportChemistry runs at line 2708.
                for (0..state.cell_count) |_sync_cell| {
                    const _sync_first = _sync_cell * state.soil_layer_capacity;
                    for (0..state.active_soil_layer_count[_sync_cell]) |_sync_local| {
                        const _sync_layer = _sync_first + _sync_local;
                        const _sync_water = state.matrix_liquid_water_m3[_sync_layer];
                        if (_sync_water <= config.negligible_quantity_threshold) continue;
                        const _sync_co2_idx = try ecosys.gas_transport.massIndex(_sync_layer, .carbon_dioxide, gas_transport_state.cell_count);
                        gas_transport_state.dissolved_mass_g[_sync_co2_idx] = initial_chemistry_state.aqueous[_sync_layer].carbon_dioxide * 12.0 * _sync_water;
                    }
                }
                try ecosys.soil_aqueous_transport_bridge.exportChemistry(&initial_chemistry_state, &micropore_solute_state);
                const totals = try reconstructLandscapeMassBalance(
                    hourly_science_context,
                );
                {
                    const area = totals.landscape_area_m2;
                    const b = try ecosys.mass_balance_audit.balance(totals);
                    std.log.debug("carbon baseline: residue={e} organic={e} co2={e} co2in={e} co2out={e} fert={e} sink={e} total_g_m2={e}", .{ totals.residue_carbon_g / area, totals.organic_carbon_g / area, totals.carbon_dioxide_carbon_g / area, totals.cumulative_carbon_dioxide_input_g / area, totals.cumulative_carbon_output_g / area, totals.cumulative_organic_fertilizer_carbon_g / area, totals.cumulative_carbon_sink_g / area, b.carbon_g / area });
                    std.log.debug("nitrogen baseline: residue={e} organic={e} n2={e} nh4={e} no3={e} n2in={e} nin={e} nout={e} fert={e} sink={e} root_uptake={e} root_exudate={e} total_g_m2={e}", .{ totals.residue_nitrogen_g / area, totals.organic_nitrogen_g / area, totals.dinitrogen_nitrogen_g / area, totals.ammonium_nitrogen_g / area, totals.nitrate_nitrogen_g / area, totals.cumulative_dinitrogen_input_g / area, totals.cumulative_nitrogen_input_g / area, totals.cumulative_nitrogen_output_g / area, totals.cumulative_organic_fertilizer_nitrogen_g / area, totals.cumulative_nitrogen_sink_g / area, totals.cumulative_plant_root_organic_nitrogen_uptake_g / area, totals.cumulative_plant_root_organic_nitrogen_exudate_g / area, b.nitrogen_g / area });
                    // EXEC-003: the oxygen baseline is captured after STARTE
                    // gas seeding, so a nonzero value here confirms seeding ran
                    // and any later climb is genuine accumulation, not fill-in.
                    std.log.debug("oxygen baseline: storage={e} storage_g_m2={e} area_m2={e}", .{ totals.oxygen_storage_g, totals.oxygen_storage_g / area, area });
                }
                if (landscape_mass_balance_state.monitor) |*monitor|
                    try monitor.reset(totals)
                else
                    landscape_mass_balance_state.monitor =
                        try ecosys.mass_balance_audit.Monitor.init(
                            totals,
                            config.mass_balance_tolerance,
                        );
            } else if (landscape_mass_balance_state.monitor == null) {
                return error.CheckpointResumeMissingMassBalanceMonitor;
            }
            const scene_expected_hours = try std.math.mul(usize, try ecosys.simulation_timeline.inclusiveDays(active_options.start_date, active_options.end_date), 24);
            while (scene_weather_hours < scene_expected_hours) {
                @memset(fire_active_this_hour, false);
                for (hourly_weather_streams, 0..) |*stream, stream_index| {
                    weather_hour_by_stream[stream_index] = try stream.next() orelse {
                        std.log.err("weather stream ended before scene: scene={d} expected_hours={d} received_hours={d} file='{s}'", .{
                            pass.scene_index + 1,
                            scene_expected_hours,
                            scene_weather_hours,
                            weather_assignment.unique_file_names[stream_index],
                        });
                        return error.WeatherStreamEndedBeforeScene;
                    };
                    if (stream_index > 0 and !sameWeatherTimestamp(
                        weather_hour_by_stream[0].timestamp,
                        weather_hour_by_stream[stream_index].timestamp,
                    )) {
                        std.log.err("weather-grid timestamp mismatch: scene={d} hour_index={d} first_file='{s}' mismatched_file='{s}'", .{
                            pass.scene_index + 1,
                            scene_weather_hours,
                            weather_assignment.unique_file_names[0],
                            weather_assignment.unique_file_names[stream_index],
                        });
                        return error.WeatherGridTimestampMismatch;
                    }
                }
                const weather_hour = weather_hour_by_stream[0];
                var timestamp = weather_hour.timestamp;
                if (timestamp.year == null) timestamp.year = active_options.start_date.year;
                // HOUR1 clears hourly plant ledgers before GROSUB management.
                // Management dispatch below may publish disturbance, harvest,
                // or combustion fluxes that must survive into this hour's
                // HCNET and output aggregation.
                if (plant_root_state) |*roots| roots.resetHourlyFluxes();
                if (canopy_carbon_exchange_state) |*ledger| ledger.reset();
                for (hourly_accumulators) |*acc| ecosys.hourly_grid_cell_accumulator_reset.reset(acc);
                const current_year = timestamp.year.?;
                const current_day_of_year = try dayOfYearFromTimestamp(timestamp);
                if (resolved_development_year == null or resolved_development_year.? != current_year) {
                    for (development_planting_dates, 0..) |date, plant| {
                        const resolved = try date.resolve(current_year);
                        development_planting_day_of_year[plant] = try date.dayOfYear(current_year);
                        development_planting_year[plant] = resolved.year;
                    }
                    resolved_development_year = current_year;
                }
                const begins_day = if (previous_weather_timestamp) |previous| !sameCalendarDay(previous, timestamp) else true;
                // HFUNC/GROSUB evaluates IFLGC/IFLGI at NFZ=1, the first
                // source substep of every hour. The converged hourly model
                // therefore evaluates activation once per hour, while the
                // active/lifecycle flags keep reconstruction single-shot.
                if (plant_phenology_state != null and plant_growth_stage_state != null and plant_dormancy_state != null) {
                    for (0..runtime_plant_count) |plant| {
                        if (!plant_phenology_state.?.reseed_pending[plant]) continue;
                        dormant_seed_branches.items.len = 0;
                        const branch_range = try plant_growth_stage_state.?.branchRange(plant);
                        const cell = plant / config.plant_populations;
                        const boundary_base = cell * (soil_geometry_state.layer_capacity + 1);
                        const soil_surface_boundary_depth_m = soil_geometry_state.boundary_depth_m[boundary_base + soil_geometry_state.first_active_layer[cell]];
                        for (branch_range.first..branch_range.end) |branch| {
                            const branch_dormancy = plant_dormancy_state.?.branches[branch];
                            try dormant_seed_branches.append(allocator, .{
                                .leafout_disabled = branch_dormancy.leafout_disabled,
                                .accumulated_leafout_h = branch_dormancy.accumulated_leafout_h,
                                .required_leafout_h = development_dormancy_parameters[plant].required_leafout_h,
                            });
                        }
                        const activation = try ecosys.plant_harvest_runtime.sourceOrderDormantSeedActivation(
                            true,
                            plant_phenology_state.?.reseed_pending[plant],
                            dormant_seed_branches.items,
                            current_day_of_year,
                            current_year,
                            soil_surface_boundary_depth_m,
                        ) orelse continue;
                        plant_phenology_state.?.reseed_pending[plant] = activation.initialization_pending;
                        development_planting_day_of_year[plant] = activation.planting_day_of_year;
                        development_planting_year[plant] = activation.planting_year;
                        if (plant_water_workspace != null) plant_water_workspace.?.seeding_depth_m[plant] = activation.seeding_depth_m;
                    }
                }
                if (management_schedule_map != null) {
                    const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    if (plant_phenology_state != null and plant_growth_stage_state != null) {
                        const activation_count = try ecosys.plant_management_dispatch.refreshPlantActivity(management_schedule_map.?, plant_management_catalog, management_date, &plant_phenology_state.?, &plant_growth_stage_state.?, plant_phenology_state.?.emerged, plant_phenology_state.?.lifecycle_initialized, newly_activated_plants);
                        if (activation_count > 0) {
                            const assignments = plant_assignments orelse return error.PlantActivationRequiresAssignments;
                            var topology_changed = false;
                            for (newly_activated_plants, 0..) |activated, plant| {
                                if (!activated) continue;
                                const cell = plant / config.plant_populations;
                                const species = plant % config.plant_populations;
                                const unit_index = plant_unit_by_cell.?[cell];
                                if (species >= assignments.units[unit_index].species.len) return error.PlantActivationSpeciesOutOfBounds;
                                const assignment = assignments.units[unit_index].species[species];
                                const traits = plant_catalog.entries.items[plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile].traits;
                                const planting = plant_management_catalog.entries.items[plant_management_catalog.find(assignment.management_file) orelse return error.MissingPlantManagementSchedule].schedule.planting;
                                var cumulative_depth_m: f64 = 0;
                                for (replant_layer_bottom_depth_m, 0..) |*bottom, layer| {
                                    cumulative_depth_m += soil_solver_property_state.layer_thickness_m[cell * config.soil_layers + layer];
                                    bottom.* = cumulative_depth_m;
                                }
                                const planting_layer = try ecosys.plant_initialization.plantingLayer(planting.seed_depth_m, replant_layer_bottom_depth_m);
                                const adjusted = try ecosys.plant_initialization.adjustedPhenology(traits, runscript.phenology_initialization_parameters);
                                const concurrent = try ecosys.plant_initialization.concurrentNodeSettings(traits.functional_type.aboveground_turnover_type, adjusted.leaf_appearance_per_h, adjusted.floral_initiation_node_count_after_seed, runscript.phenology_initialization_parameters);
                                const storage = try ecosys.plant_initialization.seedStorage(traits.morphology.seed_mass_at_planting_g, planting.population_per_m2, canopy_cell_area_m2[cell], traits.organ_nitrogen_to_carbon_ratio.grain, traits.organ_phosphorus_to_carbon_ratio.grain);
                                const thermal = try ecosys.plant_initialization.thermalAcclimation(assignment.species_file, traits.functional_type.photosynthesis_pathway, traits.functional_type.thermal_adaptation_zone, runscript.thermal_acclimation_parameters);
                                const preserve_reseed_inventories = plant_phenology_state.?.reseed_pending[plant];
                                const preserve_death_standing_dead = plant_phenology_state.?.death_replant_pending[plant];
                                const preserve_persistent_inventories = preserve_reseed_inventories or preserve_death_standing_dead;
                                const persistent_reseed_inventories = if (preserve_persistent_inventories)
                                    try ecosys.canopy_photosynthesis.capturePersistentReseedInventories(&detailed_canopy_state.?, plant)
                                else
                                    undefined;
                                const carried_standing_dead_surface_water_m3 = canopy_precipitation_retention_state.?.standing_dead_surface_water_m3[plant];

                                const old_branch_range = try detailed_canopy_state.?.branchRange(plant);
                                const old_node_range = try detailed_canopy_state.?.nodeRange(old_branch_range.first);
                                topology_changed = topology_changed or old_branch_range.end - old_branch_range.first != 1 or old_node_range.end - old_node_range.first != 1;
                                try ecosys.plant_topology.compactPlantToInitialTopology(.{
                                    .canopy = &detailed_canopy_state.?,
                                    .growth_stages = &plant_growth_stage_state.?,
                                    .dormancy = &plant_dormancy_state.?,
                                    .branch_development = &branch_development_state.?,
                                }, plant);
                                try detailed_canopy_state.?.clearPlantForReconstruction(plant);
                                try plant_root_state.?.reconstructPlant(plant, traits, planting_layer, planting.seed_depth_m, runscript.root_initialization_parameters);
                                const seed_geometry = try ecosys.plant_initialization.seedGeometry(traits.morphology.seed_mass_at_planting_g, runscript.plant_geometry_parameters);
                                try plant_root_state.?.setSeedGeometry(plant, seed_geometry.volume_m3, seed_geometry.length_m, seed_geometry.surface_area_m2);
                                try plant_litter_partition_state.?.initializePlant(plant, traits, runscript.standing_dead_partition_parameters);
                                const branch_range = try plant_growth_stage_state.?.clearPlantForReconstruction(plant);
                                try plant_dormancy_state.?.clearRangeForReconstruction(branch_range.first, branch_range.end);
                                try branch_development_state.?.clearRangeForReconstruction(branch_range.first, branch_range.end);
                                const initial_branch_end = @min(branch_range.end, branch_range.first + 1);
                                try branch_development_state.?.initializeRange(branch_range.first, initial_branch_end, adjusted.floral_initiation_node_count_after_seed, adjusted.seed_initial_node_count, traits.functional_type.growth_habit == 0 and traits.functional_type.leaf_phenology_type == 1, concurrent.perennial_node_scaling, concurrent.maximum_concurrently_growing_nodes);
                                for (plant_growth_stage_state.?.branches[initial_branch_end..branch_range.end]) |*branch| branch.dead = true;
                                plant_phenology_state.?.floral_initiated[plant] = false;
                                plant_phenology_state.?.node_initiation_rate_at_25c_per_h[plant] = adjusted.node_initiation_per_h;
                                plant_phenology_state.?.leaf_appearance_rate_at_25c_per_h[plant] = adjusted.leaf_appearance_per_h;
                                plant_phenology_state.?.initiated_node_count[plant] = adjusted.seed_initial_node_count;
                                plant_phenology_state.?.appeared_leaf_count[plant] = adjusted.seed_initial_node_count;
                                plant_phenology_state.?.node_initiation_per_timestep[plant] = 0;
                                plant_phenology_state.?.leaf_appearance_per_timestep[plant] = 0;
                                if (branch_range.first < branch_range.end) {
                                    plant_growth_stage_state.?.branches[branch_range.first].initiated_node_count = adjusted.seed_initial_node_count;
                                    plant_growth_stage_state.?.branches[branch_range.first].appeared_leaf_count = adjusted.seed_initial_node_count;
                                }
                                try ecosys.plant_initialization.initializeShootControls(&detailed_canopy_state.?, plant, planting.population_per_m2, canopy_cell_area_m2[cell], traits.water_relations.cuticular_resistance_s_per_m, traits.functional_type.photosynthesis_pathway, runscript.shoot_control_parameters);
                                try ecosys.plant_initialization.initializeThermalAcclimation(&detailed_canopy_state.?, plant, thermal);
                                if (!preserve_reseed_inventories) {
                                    try ecosys.plant_initialization.initializeSeedStorage(&detailed_canopy_state.?, plant, storage);
                                }
                                const initialized_seed_storage = ecosys.canopy_photosynthesis.ElementalMass{
                                    .carbon_g = detailed_canopy_state.?.plant_seed_storage_carbon_g[plant],
                                    .nitrogen_g = detailed_canopy_state.?.plant_seed_storage_nitrogen_g[plant],
                                    .phosphorus_g = detailed_canopy_state.?.plant_seed_storage_phosphorus_g[plant],
                                };
                                const standing_dead = try ecosys.plant_initialization.standingDeadStorage(traits.morphology.standing_dead_carbon_g_per_m2, canopy_cell_area_m2[cell], traits.organ_nitrogen_to_carbon_ratio.stalk, traits.organ_phosphorus_to_carbon_ratio.stalk, runscript.standing_dead_partition_parameters);
                                if (!preserve_reseed_inventories) {
                                    try ecosys.plant_initialization.initializeStandingDeadStorage(&detailed_canopy_state.?, plant, standing_dead);
                                }
                                try ecosys.plant_initialization.initializePlantHeatAndWater(&plant_state, &detailed_canopy_state.?, plant, site.mean_annual_air_temperature_c, traits.water_relations.osmotic_potential_megapascal, runscript.plant_heat_water_parameters);
                                if (preserve_persistent_inventories) {
                                    try ecosys.canopy_photosynthesis.restorePersistentReseedInventories(&detailed_canopy_state.?, plant, persistent_reseed_inventories);
                                    if (preserve_death_standing_dead and !preserve_reseed_inventories) {
                                        detailed_canopy_state.?.plant_seed_storage_carbon_g[plant] = initialized_seed_storage.carbon_g;
                                        detailed_canopy_state.?.plant_seed_storage_nitrogen_g[plant] = initialized_seed_storage.nitrogen_g;
                                        detailed_canopy_state.?.plant_seed_storage_phosphorus_g[plant] = initialized_seed_storage.phosphorus_g;
                                    }
                                }
                                plant_phenology_state.?.reseed_pending[plant] = false;
                                plant_phenology_state.?.death_replant_pending[plant] = false;
                                plant_state.canopy_water_storage_m_per_m2[plant] = 0;
                                plant_state.leaf_area_index_m2_m2[plant] = 0;
                                plant_state.shoot_carbon_g_m2[plant] = 0;
                                @memset(plant_state.root_carbon_g_m2[plant * config.soil_layers .. (plant + 1) * config.soil_layers], 0);
                                const retention_state = &canopy_precipitation_retention_state.?;
                                inline for (@typeInfo(ecosys.canopy_precipitation_retention.State).@"struct".fields) |field| {
                                    if (field.type == []f64 and !std.mem.startsWith(u8, field.name, "cell_")) @field(retention_state, field.name)[plant] = 0;
                                }
                                if (preserve_reseed_inventories)
                                    retention_state.standing_dead_surface_water_m3[plant] = carried_standing_dead_surface_water_m3;
                            }
                            if (topology_changed) {
                                var replacement_layers = try ecosys.canopy_layer_distribution.State.init(allocator, state.cell_count, config.plant_populations, runscript.canopy_layer_count, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len, &detailed_canopy_state.?);
                                const replacement_carbon_exchange = ecosys.canopy_carbon_exchange.State.init(allocator, detailed_canopy_state.?.branch_node_offsets.len - 1) catch |err| {
                                    replacement_layers.deinit();
                                    return err;
                                };
                                var previous_carbon_exchange = canopy_carbon_exchange_state.?;
                                canopy_carbon_exchange_state.? = replacement_carbon_exchange;
                                previous_carbon_exchange.deinit();
                                var previous_layers = canopy_layer_distribution_state.?;
                                canopy_layer_distribution_state.? = replacement_layers;
                                previous_layers.deinit();
                            }
                            if (harvest_context != null) for (0..runtime_plant_count) |plant| {
                                const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                                try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                            };
                        }
                    }
                }
                if (fertilizer_schedule_maps[pass.scene_index] != null) {
                    const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    var fertilizer_context: ecosys.fertilizer_management_dispatch.NitrogenApplyContext = .{
                        .soil = &soil_fertilizer_inventory,
                        .surface = &surface_litter_fertilizer_state,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active_soil_layer_count = state.active_soil_layer_count,
                        .soil_layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                        .nitrogen_molar_mass_g_per_mol = runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                        .surface_organic = &surface_organic_state,
                        .source_hour_one_through_twenty_four = timestamp.hour,
                        .solar_noon_hour_by_cell = fertilizer_solar_noon_hour_by_cell,
                    };
                    _ = try ecosys.fertilizer_management_dispatch.dispatchDate(fertilizer_schedule_maps[pass.scene_index].?, fertilizer_catalog, management_date, &fertilizer_context, ecosys.fertilizer_management_dispatch.applyNitrogen);
                    var mineral_fertilizer_context: ecosys.fertilizer_management_dispatch.MineralApplyContext = .{
                        .inventory = &mineral_fertilizer_inventory,
                        .surface_organic = &surface_organic_state,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active_soil_layer_count = state.active_soil_layer_count,
                        .soil_layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                        .source_hour_one_through_twenty_four = timestamp.hour,
                        .solar_noon_hour_by_cell = fertilizer_solar_noon_hour_by_cell,
                    };
                    _ = try ecosys.fertilizer_management_dispatch.dispatchDate(fertilizer_schedule_maps[pass.scene_index].?, fertilizer_catalog, management_date, &mineral_fertilizer_context, ecosys.fertilizer_management_dispatch.applyMinerals);
                    var organic_fertilizer_context: ecosys.fertilizer_management_dispatch.OrganicApplyContext = .{
                        .soil = &soil_organic_state,
                        .surface = &surface_organic_state,
                        .parameters = &organic_parameters,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active_soil_layer_count = state.active_soil_layer_count,
                        .soil_layer_capacity = config.soil_layers,
                        .soil_layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                        .daily_organic_carbon_input_g_c = daily_organic_fertilizer_carbon_input_g_c,
                        .daily_biome_carbon_input_g_c = daily_biome_organic_carbon_input_g_c,
                        .daily_organic_phosphorus_input_g_p = daily_organic_fertilizer_phosphorus_input_g_p,
                        .daily_organic_nitrogen_input_g_n = daily_organic_fertilizer_nitrogen_input_g_n,
                        .source_hour_one_through_twenty_four = timestamp.hour,
                        .solar_noon_hour_by_cell = fertilizer_solar_noon_hour_by_cell,
                    };
                    _ = try ecosys.fertilizer_management_dispatch.dispatchDate(fertilizer_schedule_maps[pass.scene_index].?, fertilizer_catalog, management_date, &organic_fertilizer_context, ecosys.fertilizer_management_dispatch.applyOrganic);
                }
                // GROSUB applies tillage at integer solar noon, not at the
                // daily management initialization boundary.
                if (@as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour) and disturbance_schedule_maps[pass.scene_index] != null) {
                    const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    var disturbance_context: ecosys.disturbance_management_dispatch.ApplyContext = .{
                        .roots = &plant_root_state.?,
                        .litter_partition = &plant_litter_partition_state.?,
                        .soil_organic = &soil_organic_state,
                        .grid = &state,
                        .species_count = config.plant_populations,
                        .biological_domain_count_by_plant = root_biological_domain_count_by_plant,
                        .root_nonwoody_fraction_by_plant = plant_water_workspace.?.woody_root_fraction,
                        .biomass_turnover_type_by_plant = canopy_layer_controls.biomass_turnover_type,
                        .root_profile_type_by_plant = canopy_layer_controls.root_profile_type,
                        .growth_habit_by_plant = plant_topology_controls.growth_habit_code,
                        .leaf_phenology_type_by_plant = plant_topology_controls.leaf_phenology_code,
                        .planting_day_of_year_by_plant = development_planting_day_of_year,
                        .planting_year_by_plant = development_planting_year,
                        .current_day_of_year = current_day_of_year,
                        .current_year = @intCast(current_year),
                        .plant_harvest = if (harvest_context) |*value| value else null,
                        .root_litter_carbon_ledger = &root_litter_carbon_ledger,
                        .surface_energy = &surface_energy_state,
                        .fire_active_this_hour = fire_active_this_hour,
                    };
                    _ = try ecosys.disturbance_management_dispatch.dispatchDatePhase(
                        disturbance_schedule_maps[pass.scene_index].?,
                        disturbance_catalog,
                        management_date,
                        &disturbance_context,
                        .pre_science,
                    );
                }
                if (begins_day and active_options.climate_change_mode == 2) {
                    try climate_state.advanceDay(active_options, ecosys.climate_change.daysInYear(timestamp.year.?));
                }
                for (0..state.cell_count) |cell| {
                    const stream_index = weather_assignment.stream_index_by_cell[cell];
                    if (stream_index >= weather_hour_by_stream.len)
                        return error.WeatherGridStreamIndexOutOfRange;
                    const selected_site = site_by_cell[cell];
                    const cell_weather_stream = weather_streams[stream_index];
                    hourly_weather_header_by_cell[cell] = cell_weather_stream.header;
                    var cell_radiation = try ecosys.atmospheric_radiation.prepare(
                        weather_hour_by_stream[stream_index].forcing,
                        timestamp,
                        selected_site.latitude_degrees_north,
                        cell_weather_stream.header.solar_noon_hour,
                        selected_site.ecosystem_type == -2,
                    );
                    hourly_forcing_by_cell[cell] = try climate_state.apply(
                        cell_radiation.forcing,
                        active_options,
                        timestamp,
                        cell_weather_stream.header.solar_noon_hour,
                        selected_site.elevation_m,
                    );
                    cell_radiation.forcing = hourly_forcing_by_cell[cell];
                    hourly_radiation_by_cell[cell] = cell_radiation;
                }
                current_atmospheric_co2_umol_per_mol = site.atmospheric_co2_umol_mol * try climate_state.atmosphericCo2Multiplier(active_options, current_day_of_year);
                if (!std.math.isFinite(current_atmospheric_co2_umol_per_mol) or current_atmospheric_co2_umol_per_mol <= 0) return error.InvalidAtmosphericCarbonDioxide;
                {
                    const parameters = surface_gas_parameters;
                    current_atmospheric_gas_concentration_g_per_m3.* = parameters.atmospheric_concentration_g_per_m3;
                    current_atmospheric_gas_concentration_g_per_m3[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] *= current_atmospheric_co2_umol_per_mol / site.atmospheric_co2_umol_mol;
                }
                const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                if (begins_day and irrigation_schedule_maps[pass.scene_index] != null) {
                    @memset(minimum_canopy_water_potential_mpa_by_cell, 0);
                    try ecosys.plant_water_balance.publishFirstPlantDailyMinimumByCell(&detailed_canopy_state.?, minimum_canopy_water_potential_mpa_by_cell);
                    _ = try ecosys.irrigation_management_dispatch.planAutomatedDay(
                        irrigation_schedule_maps[pass.scene_index].?,
                        irrigation_catalog,
                        .{
                            .date = management_date,
                            .cell_area_m2 = canopy_cell_area_m2,
                            .active_layer_count = state.active_soil_layer_count,
                            .soil_layer_capacity = state.soil_layer_capacity,
                            .layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                            .porosity_fraction = soil_solver_property_state.porosity_fraction,
                            .field_capacity_fraction = soil_field_capacity_fraction,
                            .wilting_point_fraction = soil_wilting_point_fraction,
                            .liquid_water_m3 = state.liquid_water_m3,
                            .ice_water_m3 = state.ice_water_m3,
                            .minimum_canopy_water_potential_megapascal = minimum_canopy_water_potential_mpa_by_cell,
                        },
                        automatic_irrigation_depth_m_by_cell_hour,
                    );
                }
                if (begins_day and detailed_canopy_state != null) ecosys.plant_water_balance.resetDailyMinimumCanopyWaterPotential(&detailed_canopy_state.?);
                @memset(irrigation_water_depth_m, 0);
                @memset(irrigation_dissolved_mass_g_per_m2, 0);
                @memset(irrigation_hydrogen_mol_per_m2, 0);
                irrigation_loads.reset();
                const irrigation_source_hour: u8 = if (timestamp.hour == 0) 1 else timestamp.hour;
                if (irrigation_schedule_maps[pass.scene_index]) |irrigation_map| {
                    _ = try ecosys.irrigation_management_dispatch.accumulateScheduledHourRouted(
                        irrigation_map,
                        irrigation_catalog,
                        management_date,
                        irrigation_source_hour,
                        .{
                            .cell_area_m2 = canopy_cell_area_m2,
                            .active_layer_count = state.active_soil_layer_count,
                            .layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                            .loads = &irrigation_loads,
                        },
                    );
                    _ = try ecosys.irrigation_management_dispatch.accumulatePlannedAutomatedHourRouted(
                        irrigation_map,
                        irrigation_catalog,
                        irrigation_source_hour,
                        automatic_irrigation_depth_m_by_cell_hour,
                        .{
                            .cell_area_m2 = canopy_cell_area_m2,
                            .active_layer_count = state.active_soil_layer_count,
                            .layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                            .loads = &irrigation_loads,
                        },
                    );
                }
                for (0..state.cell_count) |cell| {
                    const area_m2 = canopy_cell_area_m2[cell];
                    if (!std.math.isFinite(area_m2) or area_m2 <= 0)
                        return error.InvalidIrrigationCellArea;
                    irrigation_water_depth_m[cell] =
                        irrigation_loads.surface_water_m3[cell] / area_m2;
                    irrigation_hydrogen_mol_per_m2[cell] =
                        irrigation_loads.surface_hydrogen_mol[cell] / area_m2;
                    const first = cell *
                        ecosys.irrigation_management_dispatch.dissolved_species_count;
                    for (0..ecosys.irrigation_management_dispatch.dissolved_species_count) |species|
                        irrigation_dissolved_mass_g_per_m2[first + species] =
                            irrigation_loads.surface_dissolved_mass_g[first + species] /
                            area_m2;
                }
                const previous_day_of_year: u16 = if (current_day_of_year > 1) current_day_of_year - 1 else ecosys.climate_change.daysInYear(current_year -| 1);
                for (site_by_cell, plant_calendar_by_cell) |cell_site, *calendar| {
                    const phytotron = cell_site.ecosystem_type == -2;
                    const current_daylength_h = if (phytotron) 24 else ecosys.daily_weather_disaggregation.calculateDaylength(current_day_of_year, cell_site.latitude_degrees_north);
                    const previous_daylength_h = if (phytotron) 24 else ecosys.daily_weather_disaggregation.calculateDaylength(previous_day_of_year, cell_site.latitude_degrees_north);
                    const solstice_day: u16 = if (cell_site.latitude_degrees_north >= 0) 173 else 355;
                    const maximum_daylength_h = if (phytotron) 24 else ecosys.daily_weather_disaggregation.calculateDaylength(solstice_day, cell_site.latitude_degrees_north);
                    calendar.* = .{
                        .day_of_year = current_day_of_year,
                        .current_year = current_year,
                        .current_daylength_h = current_daylength_h,
                        .previous_daylength_h = previous_daylength_h,
                        .maximum_seasonal_daylength_h = maximum_daylength_h,
                        .latitude_degrees_north = cell_site.latitude_degrees_north,
                    };
                }
                organic_matter_fire_exchange_state.resetHourly();
                surface_fire_exchange_state.resetHourly();
                @memset(shoot_harvest_litter_carbon_g_c_by_plant, 0);
                @memset(shoot_harvest_litter_nitrogen_g_n_by_plant, 0);
                @memset(shoot_harvest_litter_phosphorus_g_p_by_plant, 0);
                @memset(hourly_manure_products_by_plant, .{});
                @memcpy(
                    harvest_carbon_at_hour_start_g_c_by_plant,
                    plant_daily_flux_ledger.harvested_carbon_g,
                );
                const gas_failure_file_name = try std.fmt.allocPrint(
                    allocator,
                    "ecosys-ng-gas-failure-ex{d}-scenario{d}-repeat{d}-scene{d}-year{d}-day{d}-hour{d}-cells0-{d}.bin",
                    .{
                        pass.execution_iteration + 1,
                        pass.scenario_index + 1,
                        pass.scenario_iteration + 1,
                        pass.scene_index + 1,
                        current_year,
                        current_day_of_year,
                        timestamp.hour,
                        state.cell_count - 1,
                    },
                );
                defer allocator.free(gas_failure_file_name);
                const gas_failure_path = try std.fs.path.join(
                    allocator,
                    &.{ runscript_directory, gas_failure_file_name },
                );
                defer allocator.free(gas_failure_path);
                const solute_failure_file_name = try std.fmt.allocPrint(
                    allocator,
                    "ecosys-ng-solute-failure-ex{d}-scenario{d}-repeat{d}-scene{d}-year{d}-day{d}-hour{d}.bin",
                    .{
                        pass.execution_iteration + 1,
                        pass.scenario_index + 1,
                        pass.scenario_iteration + 1,
                        pass.scene_index + 1,
                        current_year,
                        current_day_of_year,
                        timestamp.hour,
                    },
                );
                defer allocator.free(solute_failure_file_name);
                const solute_failure_path = try std.fs.path.join(
                    allocator,
                    &.{ runscript_directory, solute_failure_file_name },
                );
                defer allocator.free(solute_failure_path);
                executeHourlyScience(
                    hourly_science_context,
                    timestamp.hour,
                    hourly_radiation_by_cell,
                    hourly_forcing_by_cell,
                    hourly_weather_header_by_cell,
                    plant_calendar_by_cell,
                    .{ .value = executed_weather_hours + 1 },
                    .{
                        .io = init.io,
                        .directory = std.Io.Dir.cwd(),
                        .file_path = gas_failure_path,
                    },
                    .{
                        .io = init.io,
                        .directory = std.Io.Dir.cwd(),
                        .file_path = solute_failure_path,
                        .context = .{
                            .execution_id = @intCast(
                                pass.execution_iteration + 1,
                            ),
                            .scenario_id = @intCast(
                                pass.scenario_index + 1,
                            ),
                            .repeat_id = @intCast(
                                pass.scenario_iteration + 1,
                            ),
                            .scene_id = @intCast(pass.scene_index + 1),
                            .scene_hour = @intCast(scene_weather_hours + 1),
                            .year = @intCast(current_year),
                            .day_of_year = @intCast(current_day_of_year),
                            .hour = @intCast(timestamp.hour),
                        },
                    },
                ) catch |err| {
                    std.log.err("hourly science failed: execution={d} scenario={d} scenario_repeat={d} scene={d} scene_hour={d} total_hour={d} year={?} day_of_year={?} month={?} day={?} hour={d} error={s}", .{ pass.execution_iteration + 1, pass.scenario_index + 1, pass.scenario_iteration + 1, pass.scene_index + 1, scene_weather_hours + 1, executed_weather_hours + 1, timestamp.year, timestamp.day_of_year, timestamp.month, timestamp.day_of_month, timestamp.hour, @errorName(err) });
                    return err;
                };
                if (harvest_context != null and plant_phenology_state != null) {
                    var transitioned = false;
                    for (plant_phenology_state.?.leafout_transition_this_step, 0..) |leafout, plant| {
                        if (!leafout or !plant_phenology_state.?.active[plant] or plant_topology_controls.growth_habit_code[plant] == 0) continue;
                        transitioned = true;
                    }
                    if (transitioned) for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    };
                }
                // GROSUB creates an immediate grain/JHVST=2 event when the
                // main branch of a self-seeding deciduous annual enters
                // end-of-season reproductive turnover.
                if (harvest_context != null) {
                    const automatic_harvest_date = try ecosys.execution_calendar_date.fromDayOfYear(
                        current_day_of_year,
                        @intCast(current_year),
                    );
                    const automatic_self_seed_count = try ecosys.plant_harvest_runtime.applyAutomaticSelfSeedingHarvests(
                        &harvest_context.?,
                        seasonal_turnover_event_by_plant,
                        plant_topology_controls.growth_habit_code,
                        plant_topology_controls.leaf_phenology_code,
                        .{ .day = automatic_harvest_date.day, .month = automatic_harvest_date.month, .year = automatic_harvest_date.year },
                    );
                    if (automatic_self_seed_count > 0) for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    };
                }
                // Tillage is a late-GROSUB transaction. Applying it here keeps
                // current-hour growth and FWODR composition, and prevents the
                // removed roots from regrowing before the next hour.
                if (@as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour) and disturbance_schedule_maps[pass.scene_index] != null) {
                    var disturbance_context: ecosys.disturbance_management_dispatch.ApplyContext = .{
                        .roots = &plant_root_state.?,
                        .litter_partition = &plant_litter_partition_state.?,
                        .soil_organic = &soil_organic_state,
                        .grid = &state,
                        .species_count = config.plant_populations,
                        .biological_domain_count_by_plant = root_biological_domain_count_by_plant,
                        .root_nonwoody_fraction_by_plant = plant_water_workspace.?.woody_root_fraction,
                        .biomass_turnover_type_by_plant = canopy_layer_controls.biomass_turnover_type,
                        .root_profile_type_by_plant = canopy_layer_controls.root_profile_type,
                        .growth_habit_by_plant = plant_topology_controls.growth_habit_code,
                        .leaf_phenology_type_by_plant = plant_topology_controls.leaf_phenology_code,
                        .planting_day_of_year_by_plant = development_planting_day_of_year,
                        .planting_year_by_plant = development_planting_year,
                        .current_day_of_year = current_day_of_year,
                        .current_year = @intCast(current_year),
                        .plant_harvest = if (harvest_context) |*value| value else null,
                        .surface_energy = &surface_energy_state,
                        .fire_active_this_hour = fire_active_this_hour,
                    };
                    _ = try ecosys.disturbance_management_dispatch.dispatchDatePhase(
                        disturbance_schedule_maps[pass.scene_index].?,
                        disturbance_catalog,
                        try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp),
                        &disturbance_context,
                        .post_science,
                    );
                    if (harvest_context != null) for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    };
                }
                if (harvest_context != null and management_schedule_map != null and plant_phenology_state != null) {
                    @memset(grazing_average_shoot_carbon_g_c, 0);
                    @memset(grazing_active_cell_count, 0);
                    const canopy = &detailed_canopy_state.?;
                    for (0..runtime_plant_count) |plant| {
                        if (!plant_phenology_state.?.active[plant]) continue;
                        const species = plant % config.plant_populations;
                        grazing_average_shoot_carbon_g_c[species] += canopy.plant_total_shoot_carbon_g[plant];
                        grazing_active_cell_count[species] += 1;
                    }
                    for (grazing_average_shoot_carbon_g_c, grazing_active_cell_count) |*average, count| {
                        if (count > 0) average.* /= @floatFromInt(count);
                    }
                    try ecosys.plant_management_dispatch.accumulateGrazingRebalanceTotals(
                        management_schedule_map.?,
                        plant_management_catalog,
                        management_date,
                        canopy.plant_total_shoot_carbon_g,
                        plant_phenology_state.?.active,
                        grazing_total_shoot_carbon_g_c_by_cell,
                        grazing_total_grazer_biomass_by_cell,
                    );
                    const GrazingDispatch = struct {
                        harvest: *ecosys.plant_harvest_runtime.Context,
                        average_shoot_carbon_g_c: []const f64,
                        cell_area_m2: []const f64,
                        active: []const bool,
                        turnover_type: []const u8,
                        root_profile_type: []const u8,
                        total_shoot_carbon_g_c_by_cell: []const f64,
                        total_grazer_biomass_by_cell: []const f64,
                        species_count: usize,
                        monthly_self_thinning_hour: bool,

                        fn apply(self: *@This(), plant: usize, event: ecosys.plant_management.HarvestEvent) !void {
                            if (!self.active[plant]) return;
                            if (self.monthly_self_thinning_hour and self.turnover_type[plant] != 0 and self.root_profile_type[plant] > 1) return;
                            const species = plant % self.species_count;
                            const cell = plant / self.species_count;
                            var effective_event = event;
                            if (event.kind == .animal_grazing and self.total_shoot_carbon_g_c_by_cell[cell] > 0) {
                                effective_event.cutting_height_m_or_lai_fraction =
                                    self.total_grazer_biomass_by_cell[cell] *
                                    self.harvest.canopy_state.plant_total_shoot_carbon_g[plant] /
                                    self.total_shoot_carbon_g_c_by_cell[cell];
                            }
                            _ = try ecosys.plant_harvest_runtime.applyGrazingEvent(
                                self.harvest,
                                plant,
                                effective_event,
                                self.average_shoot_carbon_g_c[species],
                                self.cell_area_m2[cell],
                            );
                        }
                    };
                    var grazing_dispatch: GrazingDispatch = .{
                        .harvest = &harvest_context.?,
                        .average_shoot_carbon_g_c = grazing_average_shoot_carbon_g_c,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active = plant_phenology_state.?.active,
                        .turnover_type = canopy_layer_controls.biomass_turnover_type,
                        .root_profile_type = canopy_layer_controls.root_profile_type,
                        .total_shoot_carbon_g_c_by_cell = grazing_total_shoot_carbon_g_c_by_cell,
                        .total_grazer_biomass_by_cell = grazing_total_grazer_biomass_by_cell,
                        .species_count = config.plant_populations,
                        .monthly_self_thinning_hour = current_day_of_year % 30 == 0 and
                            @as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour),
                    };
                    _ = try ecosys.plant_management_dispatch.dispatchDateGrazing(
                        management_schedule_map.?,
                        plant_management_catalog,
                        management_date,
                        &grazing_dispatch,
                        GrazingDispatch.apply,
                    );
                    for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    }
                }
                // GROSUB monthly forest self-thinning is evaluated at local
                // solar noon after the current-hour stem diameter refresh.
                if (harvest_context != null and current_day_of_year % 30 == 0 and
                    @as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour))
                {
                    for (0..runtime_plant_count) |plant| {
                        if (canopy_layer_controls.biomass_turnover_type[plant] == 0 or
                            canopy_layer_controls.root_profile_type[plant] <= 1 or
                            plant_phenology_state == null or !plant_phenology_state.?.active[plant]) continue;
                        if (management_schedule_map != null and
                            !try ecosys.plant_management_dispatch.allowsForestSelfThinning(
                                management_schedule_map.?,
                                plant_management_catalog,
                                try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp),
                                plant,
                            )) continue;
                        const thinned_fraction = try ecosys.plant_harvest_runtime.applyForestSelfThinning(&harvest_context.?, plant);
                        if (thinned_fraction > 0) {
                            const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                            try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                        }
                    }
                }
                if (harvest_context != null and management_schedule_map != null and
                    @as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour))
                {
                    const harvest_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    _ = try ecosys.plant_management_dispatch.dispatchDateNonGrazing(
                        management_schedule_map.?,
                        plant_management_catalog,
                        harvest_date,
                        &harvest_context.?,
                        ecosys.plant_harvest_runtime.applyEvent,
                    );
                    for (harvest_products_by_plant.?, 0..) |_, plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    }
                }
                try daily_soil_gas_flux.accumulateHour(
                    if (plant_root_state) |*roots| roots else null,
                    config.plant_populations,
                    config.soil_layers,
                    soil_gas_transport_state.atmospheric_flux_g_per_h,
                    surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                );
                try daily_soil_gas_flux.accumulateSubsurfaceNitrogenBoundaryHour(
                    config.soil_layers,
                    state.active_soil_layer_count,
                    soil_gas_transport_state.subsurface_flux_g_per_h,
                );
                // REDIST lines 4465--4507: accumulate soil/litter gas boundary
                // fluxes and precipitation/irrigation dissolved gases into the
                // source-signed hourly element ledgers.
                for (0..state.cell_count) |cell| {
                    const first_layer = cell * config.soil_layers;
                    const layer_count = state.active_soil_layer_count[cell];
                    const co2_species = @intFromEnum(ecosys.gas_transport.Species.carbon_dioxide);
                    const ch4_species = @intFromEnum(ecosys.gas_transport.Species.methane);
                    const o2_species = @intFromEnum(ecosys.gas_transport.Species.oxygen);
                    var co2_soil_surface_flux: f64 = 0;
                    var ch4_soil_surface_flux: f64 = 0;
                    var o2_soil_surface_flux: f64 = 0;
                    const co2_litter_flux = surface_litter_gas_transport_state.atmospheric_flux_g_per_h[cell * ecosys.gas_transport.species_count + co2_species];
                    const ch4_litter_flux = surface_litter_gas_transport_state.atmospheric_flux_g_per_h[cell * ecosys.gas_transport.species_count + ch4_species];
                    const o2_litter_flux = surface_litter_gas_transport_state.atmospheric_flux_g_per_h[cell * ecosys.gas_transport.species_count + o2_species];
                    var subsurface_irrigation_m3_per_h: f64 = 0;
                    for (first_layer..first_layer + layer_count) |layer| {
                        subsurface_irrigation_m3_per_h += irrigation_loads.subsurface_water_m3[layer];
                    }
                    for (first_layer..first_layer + layer_count) |layer| {
                        co2_soil_surface_flux += soil_gas_transport_state.atmospheric_flux_g_per_h[layer * ecosys.gas_transport.species_count + co2_species];
                        ch4_soil_surface_flux += soil_gas_transport_state.atmospheric_flux_g_per_h[layer * ecosys.gas_transport.species_count + ch4_species];
                        o2_soil_surface_flux += soil_gas_transport_state.atmospheric_flux_g_per_h[layer * ecosys.gas_transport.species_count + o2_species];
                    }
                    const carbon_fluxes = ecosys.redist_surface_gas_flux_accounting.SurfaceGasFluxes{
                        .soil_surface_exchange_g = co2_soil_surface_flux,
                        .mineral_layer_convective_g = 0,
                        .bulk_mass_transfer_g = 0,
                        .litter_layer_convective_g = co2_litter_flux,
                        .litter_volatilization_g = 0,
                    };
                    const methane_fluxes = ecosys.redist_surface_gas_flux_accounting.SurfaceGasFluxes{
                        .soil_surface_exchange_g = ch4_soil_surface_flux,
                        .mineral_layer_convective_g = 0,
                        .bulk_mass_transfer_g = 0,
                        .litter_layer_convective_g = ch4_litter_flux,
                        .litter_volatilization_g = 0,
                    };
                    const oxygen_fluxes = ecosys.redist_surface_gas_flux_accounting.SurfaceGasFluxes{
                        .soil_surface_exchange_g = o2_soil_surface_flux,
                        .mineral_layer_convective_g = 0,
                        .bulk_mass_transfer_g = 0,
                        .litter_layer_convective_g = o2_litter_flux,
                        .litter_volatilization_g = 0,
                    };
                    const aqueous = ecosys.redist_surface_gas_flux_accounting.AqueousFluxes{
                        .rain_to_surface_m3 = surface_precipitation_state.water_to_matrix_m3_per_h[cell] + surface_precipitation_state.water_to_macropore_m3_per_h[cell],
                        .rain_to_litter_m3 = surface_precipitation_state.water_to_litter_m3_per_h[cell],
                        .irrigation_to_surface_m3 = irrigation_loads.surface_water_m3[cell],
                        .irrigation_to_litter_m3 = 0,
                    };
                    var precipitation_gas = std.mem.zeroes(
                        ecosys.precipitation_irrigation_dissolved_gases.DissolvedGasConcentrations,
                    );
                    var irrigation_gas = std.mem.zeroes(
                        ecosys.precipitation_irrigation_dissolved_gases.DissolvedGasConcentrations,
                    );
                    {
                        const parameters = surface_gas_parameters;
                        const dissolved = try ecosys.precipitation_irrigation_dissolved_gases.calculate(
                            .{
                                .carbon_dioxide_g_m3 = current_atmospheric_gas_concentration_g_per_m3.*[co2_species],
                                .methane_g_m3 = current_atmospheric_gas_concentration_g_per_m3.*[ch4_species],
                                .oxygen_g_m3 = current_atmospheric_gas_concentration_g_per_m3.*[o2_species],
                                .nitrogen_g_m3 = current_atmospheric_gas_concentration_g_per_m3.*[@intFromEnum(ecosys.gas_transport.Species.nitrogen)],
                                .nitrous_oxide_g_m3 = current_atmospheric_gas_concentration_g_per_m3.*[@intFromEnum(ecosys.gas_transport.Species.nitrous_oxide)],
                            },
                            .{
                                .carbon_dioxide_ratio = parameters.solubility.reference_water_to_air[co2_species],
                                .methane_ratio = parameters.solubility.reference_water_to_air[ch4_species],
                                .oxygen_ratio = parameters.solubility.reference_water_to_air[o2_species],
                                .nitrogen_ratio = parameters.solubility.reference_water_to_air[@intFromEnum(ecosys.gas_transport.Species.nitrogen)],
                                .nitrous_oxide_ratio = parameters.solubility.reference_water_to_air[@intFromEnum(ecosys.gas_transport.Species.nitrous_oxide)],
                            },
                            .{
                                .carbon_dioxide = parameters.precipitation_activity_log[co2_species],
                                .methane = parameters.precipitation_activity_log[ch4_species],
                                .oxygen = parameters.precipitation_activity_log[o2_species],
                                .nitrogen = parameters.precipitation_activity_log[@intFromEnum(ecosys.gas_transport.Species.nitrogen)],
                                .nitrous_oxide = parameters.precipitation_activity_log[@intFromEnum(ecosys.gas_transport.Species.nitrous_oxide)],
                            },
                            atmospheric_state.air_temperature_k[cell] - 273.15,
                        );
                        precipitation_gas = dissolved.precipitation;
                        irrigation_gas = dissolved.irrigation;
                    }
                    const carbon_flux = try ecosys.redist_surface_gas_flux_accounting.computeCarbonGas(.{
                        .co2_fluxes = carbon_fluxes,
                        .ch4_fluxes = methane_fluxes,
                        .aqueous = aqueous,
                        .co2_rain_concentration_g_per_m3 = precipitation_gas.carbon_dioxide_g_m3,
                        .co2_irrigation_concentration_g_per_m3 = irrigation_gas.carbon_dioxide_g_m3,
                        .ch4_rain_concentration_g_per_m3 = precipitation_gas.methane_g_m3,
                        .ch4_irrigation_concentration_g_per_m3 = irrigation_gas.methane_g_m3,
                        .subsurface_irrigation_m3_per_h = subsurface_irrigation_m3_per_h,
                        .litter_subsurface_co2_transfer_g = 0,
                        .timestep_h = 1,
                    });
                    const oxygen_flux = try ecosys.redist_surface_gas_flux_accounting.computeOxygenHydrogen(.{
                        .o2_fluxes = oxygen_fluxes,
                        .h2_fluxes = std.mem.zeroes(ecosys.redist_surface_gas_flux_accounting.SurfaceGasFluxes),
                        .aqueous = aqueous,
                        .o2_rain_concentration_g_per_m3 = precipitation_gas.oxygen_g_m3,
                        .o2_irrigation_concentration_g_per_m3 = irrigation_gas.oxygen_g_m3,
                        .subsurface_irrigation_m3_per_h = subsurface_irrigation_m3_per_h,
                        // RUPOXO is driven by surface microbial O2 uptake in this port.
                        // ROGOX/RC4OX remain intentionally zero pending a dedicated surface
                        // methane-oxidation/combustion and oxygen-limited uptake export.
                        .litter_microbial_o2_uptake_g = blk: {
                            var microbial_o2_uptake_g: f64 = 0;
                            const first = cell * ecosys.surface_microbial_respiration_step.unit_count_per_cell;
                            for (0..ecosys.surface_microbial_respiration_step.unit_count_per_cell) |unit| {
                                microbial_o2_uptake_g +=
                                    surface_microbial_oxygen_state.allocation.oxygen_uptake_g_o[first + unit];
                            }
                            break :blk microbial_o2_uptake_g;
                        },
                        .litter_o2_limited_uptake_g = 0,
                        .litter_ch4_combustion_g_c = 0,
                        .litter_h2_output_g = -surface_microbial_oxygen_state.respiration_hydrogen_g_h_per_step[cell],
                        .timestep_h = 1,
                    });
                    hourly_accumulators[cell].canopy_co2_exchange_g_c_timestep += carbon_flux.co2_surface_input_g;
                    hourly_accumulators[cell].canopy_ch4_exchange_g_c_timestep += carbon_flux.ch4_surface_input_g;
                    hourly_accumulators[cell].canopy_o2_exchange_g_o_timestep += oxygen_flux.o2_surface_input_g;
                    hourly_accumulators[cell].redist_carbon_surface_input_g_c_timestep +=
                        carbon_flux.co2_surface_input_g +
                        carbon_flux.ch4_surface_input_g;
                    hourly_accumulators[cell].redist_carbon_subsurface_output_g_c_timestep +=
                        carbon_flux.total_subsurface_carbon_flux_g;
                    hourly_accumulators[cell].redist_oxygen_surface_input_g_o_timestep +=
                        oxygen_flux.o2_surface_input_g;
                    hourly_accumulators[cell].redist_oxygen_subsurface_output_g_o_timestep +=
                        oxygen_flux.o2_subsurface_out_g;
                    hourly_accumulators[cell].redist_hydrogen_surface_input_g_h_timestep +=
                        oxygen_flux.h2_surface_input_g;
                    hourly_accumulators[cell].redist_hydrogen_subsurface_output_g_h_timestep +=
                        oxygen_flux.h2_subsurface_out_g;
                    try daily_soil_gas_flux.accumulateRedistSurfaceGasHour(
                        cell,
                        carbon_flux.co2_surface_input_g + carbon_flux.ch4_surface_input_g,
                        carbon_flux.total_subsurface_carbon_flux_g,
                        oxygen_flux.o2_surface_input_g,
                        oxygen_flux.o2_subsurface_out_g,
                        oxygen_flux.h2_surface_input_g,
                        oxygen_flux.h2_subsurface_out_g,
                    );
                }
                for (0..state.cell_count) |cell| {
                    var evaporation_m3: f64 = ground_surface_evaporation_m3_per_h[cell];
                    for (0..config.plant_populations) |species| {
                        const plant = cell * config.plant_populations + species;
                        if (canopy_surface_exchange_state) |exchange| {
                            evaporation_m3 += @max(0, -exchange.transpiration_m3_per_h[plant]);
                            evaporation_m3 += @max(0, -exchange.intercepted_water_change_m3_per_h[plant]);
                        }
                        evaporation_m3 += @max(0, -standing_dead_evaporation_m3_per_h[plant]);
                    }
                    var external_water_outflow_m3: f64 = 0;
                    var external_water_inflow_m3: f64 = 0;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        external_water_outflow_m3 += @max(0, transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer]);
                        external_water_outflow_m3 += @max(0, transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]);
                        external_water_inflow_m3 += @max(0, -transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer]);
                        external_water_inflow_m3 += @max(0, -transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]);
                    }
                    const top_layer = try state.layerIndex(cell, 0);
                    const bulk_density_megagrams_per_m3 = soil_solver_property_state.bulk_density_megagrams_per_m3[top_layer];
                    if (!std.math.isFinite(bulk_density_megagrams_per_m3) or bulk_density_megagrams_per_m3 <= 0) return error.InvalidOutputSoilBulkDensity;
                    try daily_water_ledger.accumulateCell(cell, .{
                        .rainfall_m3 = atmospheric_state.precipitation_m[cell] * canopy_cell_area_m2[cell] + ground_surface_condensation_m3_per_h[cell],
                        .boundary_water_inflow_m3 = external_water_inflow_m3,
                        .evaporation_m3 = evaporation_m3,
                        .runoff_m3 = @max(0, -surface_runoff_state.exported_water_m3[cell]),
                        .water_outflow_m3 = external_water_outflow_m3,
                        .lateral_water_outflow_m3 = transport_hydrology_state.artificial_drainage_outflow_m3_per_step[cell],
                        .sediment_outflow_m3 = surface_erosion_state.routing.sediment_export_megagrams[cell] / bulk_density_megagrams_per_m3,
                    });
                    var ionic_outflow_mol: f64 = 0;
                    const active_layers = state.active_soil_layer_count[cell];
                    const first_layer = try state.layerIndex(cell, 0);
                    for (0..active_layers) |local_layer| {
                        const base = (first_layer + local_layer) * aqueous_species_count;
                        for (soil_solute_boundary_net_flux_mol[base .. base + aqueous_species_count]) |flux_mol| ionic_outflow_mol += @max(0, -flux_mol);
                    }
                    try daily_heat_ledger.accumulateHour(
                        cell,
                        active_layers,
                        atmospheric_state.shortwave_radiation_megajoules_per_m2[cell],
                        atmospheric_state.air_temperature_k[cell],
                        atmospheric_state.vapor_pressure_kpa[cell],
                        atmospheric_state.wind_speed_m_per_h[cell],
                        atmospheric_state.precipitation_m[cell],
                        state.soil_temperature_k[first_layer .. first_layer + active_layers],
                        state.surface_temperature_k[cell],
                        ionic_outflow_mol,
                    );
                }
                try daily_heterotrophic_respiration.accumulateHour(&state);
                try daily_carbon_export.accumulateOrganicDrainageHour(&state, &soil_organic_transport_state);
                try daily_carbon_export.accumulateOrganicRunoffHour(surface_organic_carbon_export_g_c_per_h);
                try daily_carbon_export.accumulateInorganicRunoffHour(surface_inorganic_carbon_export_g_c_per_h);
                try daily_carbon_export.accumulateDissolvedInorganicDrainageHour(&state, &soil_dissolved_gas_transport_state);
                try daily_carbon_export.accumulateGaseousInorganicDrainageHour(&state, &soil_gas_transport_state);
                try daily_carbon_export.accumulateSoluteCarbonateDrainageHour(
                    &state,
                    soil_solute_boundary_net_flux_mol,
                    aqueous_species_count,
                    12,
                );
                try daily_phosphorus_export.accumulateOrganicDrainageHour(&state, &soil_organic_transport_state);
                try daily_phosphorus_export.accumulateInorganicDrainageHour(
                    &state,
                    soil_solute_boundary_net_flux_mol,
                    runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                );
                try daily_phosphorus_export.accumulateRunoffHour(
                    surface_organic_phosphorus_export_g_p_per_h,
                    surface_inorganic_phosphorus_export_g_p_per_h,
                );
                try daily_nitrogen_export.accumulateOrganicDrainageHour(&state, &soil_organic_transport_state);
                try daily_nitrogen_export.accumulateInorganicDrainageHour(&state, &mineral_nitrogen_transport_state);
                try daily_nitrogen_export.accumulateDissolvedGasDrainageHour(
                    &state,
                    &soil_dissolved_gas_transport_state,
                );
                try daily_nitrogen_export.accumulateRunoffHour(
                    surface_organic_nitrogen_export_g_n_per_h,
                    surface_inorganic_nitrogen_export_g_n_per_h,
                );
                if (scene_weather_hours < 24) {
                    var diagnostic_gas_input_g_n: f64 = 0;
                    var diagnostic_aqueous_input_g_n: f64 = 0;
                    var diagnostic_export_g_n: f64 = 0;
                    var diagnostic_aqueous_input_g_p: f64 = 0;
                    var diagnostic_export_g_p: f64 = 0;
                    var diagnostic_input_g_p: f64 = 0;
                    var diagnostic_rain_m3: f64 = 0;
                    var diagnostic_water_output_m3: f64 = 0;
                    var diagnostic_ground_evaporation_m3: f64 = 0;
                    var diagnostic_ground_condensation_m3: f64 = 0;
                    var diagnostic_transpiration_m3: f64 = 0;
                    var diagnostic_living_interception_evaporation_m3: f64 = 0;
                    var diagnostic_dead_evaporation_m3: f64 = 0;
                    var diagnostic_runoff_m3: f64 = 0;
                    var diagnostic_external_outflow_m3: f64 = 0;
                    var diagnostic_lateral_outflow_m3: f64 = 0;
                    for (0..state.cell_count) |cell| {
                        inline for ([_]ecosys.gas_transport.Species{ .nitrogen, .nitrous_oxide, .ammonia }) |species|
                            diagnostic_gas_input_g_n += try daily_soil_gas_flux.get(cell, species);
                        inline for ([_]ecosys.snow_solute_transport.Species{ .ammonium_nitrogen, .ammonia_nitrogen, .nitrate_nitrogen }) |species|
                            diagnostic_aqueous_input_g_n += try atmospheric_solute_input_ledger_state.speciesInputG(cell, species);
                        inline for ([_]ecosys.snow_solute_transport.Species{ .dinitrogen_nitrogen, .nitrous_oxide_nitrogen }) |species|
                            diagnostic_gas_input_g_n += try atmospheric_solute_input_ledger_state.speciesInputG(cell, species);
                        diagnostic_export_g_n += daily_nitrogen_export.dissolved_organic_nitrogen_runoff_g_n[cell] +
                            daily_nitrogen_export.dissolved_organic_nitrogen_drainage_g_n[cell] +
                            daily_nitrogen_export.dissolved_inorganic_nitrogen_runoff_g_n[cell] +
                            daily_nitrogen_export.dissolved_inorganic_nitrogen_drainage_g_n[cell];
                        inline for ([_]ecosys.snow_solute_transport.Species{ .hydrogen_phosphate_phosphorus, .dihydrogen_phosphate_phosphorus }) |species|
                            diagnostic_aqueous_input_g_p += try atmospheric_solute_input_ledger_state.speciesInputG(cell, species);
                        diagnostic_input_g_p += daily_phosphorus_export.dissolved_organic_phosphorus_input_g_p[cell] +
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_input_g_p[cell];
                        diagnostic_export_g_p += daily_phosphorus_export.dissolved_organic_phosphorus_runoff_g_p[cell] +
                            daily_phosphorus_export.dissolved_organic_phosphorus_drainage_g_p[cell] +
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_runoff_g_p[cell] +
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_drainage_g_p[cell];
                        diagnostic_rain_m3 += daily_water_ledger.rainfall_m3[cell] + daily_water_ledger.boundary_water_inflow_m3[cell];
                        diagnostic_ground_evaporation_m3 += ground_surface_evaporation_m3_per_h[cell];
                        diagnostic_ground_condensation_m3 += ground_surface_condensation_m3_per_h[cell];
                        for (0..config.plant_populations) |species| {
                            const plant = cell * config.plant_populations + species;
                            if (canopy_surface_exchange_state) |exchange| {
                                diagnostic_transpiration_m3 += @max(0, -exchange.transpiration_m3_per_h[plant]);
                                diagnostic_living_interception_evaporation_m3 += @max(0, -exchange.intercepted_water_change_m3_per_h[plant]);
                            }
                            diagnostic_dead_evaporation_m3 += @max(0, -standing_dead_evaporation_m3_per_h[plant]);
                        }
                        diagnostic_water_output_m3 += daily_water_ledger.runoff_m3[cell] +
                            daily_water_ledger.evaporation_m3[cell] +
                            daily_water_ledger.water_outflow_m3[cell];
                        diagnostic_runoff_m3 += daily_water_ledger.runoff_m3[cell];
                        diagnostic_external_outflow_m3 += daily_water_ledger.water_outflow_m3[cell];
                        diagnostic_lateral_outflow_m3 += daily_water_ledger.lateral_water_outflow_m3[cell];
                    }
                    const diagnostic_totals = try reconstructLandscapeMassBalance(hourly_science_context);
                    std.log.debug("hourly water closure: hour={d} residual_m3={e} storage_change_m3={e} rain_m3={e} output_m3={e} current_ground_evaporation_m3={e} current_ground_condensation_m3={e}", .{ scene_weather_hours + 1, diagnostic_totals.water_storage_m3 - landscape_mass_balance_state.monitor.?.baseline.water_m3 - diagnostic_rain_m3 + diagnostic_water_output_m3, diagnostic_totals.water_storage_m3 - landscape_mass_balance_state.monitor.?.baseline.water_m3, diagnostic_rain_m3, diagnostic_water_output_m3, diagnostic_ground_evaporation_m3, diagnostic_ground_condensation_m3 });
                    const diagnostic_boundary = landscape_mass_balance_state.boundary_ledger.cumulative;
                    // NOTE: this read precedes the surface/canopy heat booking
                    // (`accumulateAcceptedSurfaceAndCanopyHeat`, ~700 lines below),
                    // so `heat_output_megajoules` here lags the storage term by one
                    // hour and this residual is NOT the audited quantity. Correcting
                    // the phase by hand makes the cumulative residual converge on the
                    // day-1 audit total (ratio 0.969 -> 0.992 -> 1.0, hours 20--23).
                    std.log.debug("hourly heat closure: hour={d} residual_megajoules={e} storage_change_megajoules={e} heat_input_megajoules={e} heat_output_megajoules={e}", .{ scene_weather_hours + 1, diagnostic_totals.heat_storage_megajoules - landscape_mass_balance_state.monitor.?.baseline.heat_megajoules - diagnostic_boundary.heat_input_megajoules + diagnostic_boundary.heat_output_megajoules, diagnostic_totals.heat_storage_megajoules - landscape_mass_balance_state.monitor.?.baseline.heat_megajoules, diagnostic_boundary.heat_input_megajoules, diagnostic_boundary.heat_output_megajoules });
                    // EXEC-002 follow-on: the oxygen audit fails on day 1 and
                    // blocks any multi-day heat measurement, so report its hourly
                    // closure on the same footing as heat.
                    std.log.debug("hourly oxygen closure: hour={d} residual_g={e} storage_change_g={e} oxygen_input_g={e} oxygen_output_g={e}", .{ scene_weather_hours + 1, diagnostic_totals.oxygen_storage_g - landscape_mass_balance_state.monitor.?.baseline.oxygen_g - diagnostic_boundary.oxygen_input_g + diagnostic_boundary.oxygen_output_g, diagnostic_totals.oxygen_storage_g - landscape_mass_balance_state.monitor.?.baseline.oxygen_g, diagnostic_boundary.oxygen_input_g, diagnostic_boundary.oxygen_output_g });
                    // EXEC-004: carbon closes exactly (`0e0`) with surface
                    // evaporation disabled but fails the day-1 audit once it is
                    // enabled, so report its hourly closure to localise the hour.
                    // Mirrors the `mass_balance_audit.balance` carbon expression.
                    {
                        // Use the audit's own balance rather than re-summing the
                        // terms here. `mass_balance_audit.balance` applies a
                        // compensated (Kahan) sum, and these terms are of order
                        // 4e11 g C, so a naive sum differs from the audited value by
                        // ~1.6e8 g C of accumulated rounding: enough to invent an
                        // apparent defect that the audit does not see.
                        //
                        // IMPORTANT: this hourly residual is NOT comparable to the
                        // daily audited carbon deviation. The plant carbon sink and
                        // the CO2 boundary terms accumulate once per day, after the
                        // hourly diagnostics run, so `cumulative_carbon_sink_g` and
                        // `cumulative_carbon_dioxide_input_g` read zero here for the
                        // whole day. In the shipped configuration this diagnostic
                        // therefore reports `-1.589e8` g C at hour 24 while the daily
                        // audit reports carbon closing at exactly `0e0`. Use it only
                        // to compare hours WITHIN a day, or to compare the same hour
                        // across two configurations. The same caveat applies to the
                        // hourly oxygen closure above, as recorded under EXEC-003.
                        const carbon_balance_g = (try ecosys.mass_balance_audit.balance(diagnostic_totals)).carbon_g;
                        std.log.debug("hourly carbon closure: hour={d} residual_g={e} residue_g={e} organic_g={e} co2_g={e} co2_input_g={e} output_g={e}", .{ scene_weather_hours + 1, carbon_balance_g - landscape_mass_balance_state.monitor.?.baseline.carbon_g, diagnostic_totals.residue_carbon_g, diagnostic_totals.organic_carbon_g, diagnostic_totals.carbon_dioxide_carbon_g, diagnostic_totals.cumulative_carbon_dioxide_input_g, diagnostic_totals.cumulative_carbon_output_g });
                        // EXEC-004: break the surface organic carbon down by pool
                        // family so the family losing the unaccounted carbon can be
                        // named rather than guessed at.
                        var surface_structural_g_c: f64 = 0;
                        var surface_residue_g_c: f64 = 0;
                        var surface_dissolved_g_c: f64 = 0;
                        var surface_adsorbed_g_c: f64 = 0;
                        var surface_microbial_g_c: f64 = 0;
                        var surface_acetate_g_c: f64 = 0;
                        for (surface_organic_state.structural) |pool| surface_structural_g_c += pool.carbon_g_c;
                        for (surface_organic_state.residue) |pool| surface_residue_g_c += pool.carbon_g_c;
                        for (surface_organic_state.dissolved) |pool| surface_dissolved_g_c += pool.carbon_g_c;
                        for (surface_organic_state.adsorbed) |pool| surface_adsorbed_g_c += pool.carbon_g_c;
                        for (surface_organic_state.microbial) |pool| surface_microbial_g_c += pool.carbon_g_c;
                        for (surface_organic_state.dissolved_acetate_carbon_g_c) |value| surface_acetate_g_c += value;
                        for (surface_organic_state.adsorbed_acetate_carbon_g_c) |value| surface_acetate_g_c += value;
                        std.log.debug("surface organic carbon families: hour={d} structural_g_c={e} residue_g_c={e} dissolved_g_c={e} adsorbed_g_c={e} microbial_g_c={e} acetate_g_c={e}", .{ scene_weather_hours + 1, surface_structural_g_c, surface_residue_g_c, surface_dissolved_g_c, surface_adsorbed_g_c, surface_microbial_g_c, surface_acetate_g_c });
                        // EXEC-004: `aggregateSurfaceOrganic` skips microbial
                        // substrate 4 (humus) with `if (substrate == 4) continue`,
                        // so carbon growing there is counted by neither the residue
                        // nor the organic term. Split the microbial total by whether
                        // the census counts it, to test that directly.
                        var counted_microbial_g_c: f64 = 0;
                        var skipped_microbial_g_c: f64 = 0;
                        for (0..state.cell_count) |cell| {
                            for (0..ecosys.soil_organic_initialization.microbial_substrate_count) |substrate| {
                                const first = (cell * ecosys.soil_organic_initialization.microbial_substrate_count + substrate) *
                                    ecosys.soil_organic_initialization.microbial_population_count *
                                    ecosys.soil_organic_initialization.kinetic_fraction_count;
                                const count = ecosys.soil_organic_initialization.microbial_population_count *
                                    ecosys.soil_organic_initialization.kinetic_fraction_count;
                                for (surface_organic_state.microbial[first .. first + count]) |pool| {
                                    if (substrate == 4) skipped_microbial_g_c += pool.carbon_g_c else counted_microbial_g_c += pool.carbon_g_c;
                                }
                            }
                        }
                        std.log.debug("surface microbial census split: hour={d} counted_g_c={e} skipped_humus_g_c={e}", .{ scene_weather_hours + 1, counted_microbial_g_c, skipped_microbial_g_c });
                        // EXEC-005: the litter inorganic carbon term is counted as
                        // `carbon_g_per_mol * live_water * (carbonate + bicarbonate +
                        // calcite)`. If the live-carrier coupling is the mechanism,
                        // this term must change when evaporation changes the water
                        // while the rebase rescales the concentrations. Report it and
                        // its two factors separately so the product can be checked.
                        var litter_inorganic_carbon_g_c: f64 = 0;
                        var litter_carbon_concentration_sum: f64 = 0;
                        var litter_carbon_water_m3: f64 = 0;
                        for (surface_litter_chemistry_state.cells, 0..) |litter_cell, litter_index| {
                            const water = surface_precipitation_state.litter_water_m3[litter_index];
                            const concentration = litter_cell.carbonate_mol_per_m3 +
                                litter_cell.bicarbonate_mol_per_m3 +
                                litter_cell.salt_minerals.calcite_mol_per_m3;
                            litter_carbon_water_m3 += water;
                            litter_carbon_concentration_sum += concentration;
                            litter_inorganic_carbon_g_c += 12.0 * water * concentration;
                        }
                        std.log.debug("litter inorganic carbon: hour={d} carbon_g_c={e} water_m3={e} concentration_sum_mol_per_m3={e}", .{ scene_weather_hours + 1, litter_inorganic_carbon_g_c, litter_carbon_water_m3, litter_carbon_concentration_sum });
                    }
                    std.log.debug("current evap components: transpiration_m3={e} living_interception_m3={e} standing_dead_m3={e}", .{ diagnostic_transpiration_m3, diagnostic_living_interception_evaporation_m3, diagnostic_dead_evaporation_m3 });
                    std.log.debug("cumulative water output components: runoff_m3={e} evaporation_m3={e} external_outflow_m3={e} lateral_outflow_m3={e}", .{ diagnostic_runoff_m3, diagnostic_water_output_m3 - diagnostic_runoff_m3 - diagnostic_external_outflow_m3 - diagnostic_lateral_outflow_m3, diagnostic_external_outflow_m3, diagnostic_lateral_outflow_m3 });
                    const diagnostic_storage_g_n = diagnostic_totals.residue_nitrogen_g + diagnostic_totals.organic_nitrogen_g + diagnostic_totals.dinitrogen_nitrogen_g + diagnostic_totals.ammonium_nitrogen_g + diagnostic_totals.nitrate_nitrogen_g;
                    const diagnostic_baseline_g_n = landscape_mass_balance_state.monitor.?.baseline.nitrogen_g;
                    std.log.debug("hourly nitrogen closure: hour={d} residual_g={e} storage_change_g={e} gas_input_g={e} aqueous_input_g={e} export_g={e}", .{ scene_weather_hours + 1, diagnostic_storage_g_n - diagnostic_baseline_g_n - diagnostic_gas_input_g_n - diagnostic_aqueous_input_g_n + diagnostic_export_g_n, diagnostic_storage_g_n - diagnostic_baseline_g_n, diagnostic_gas_input_g_n, diagnostic_aqueous_input_g_n, diagnostic_export_g_n });
                    const diagnostic_storage_g_p = diagnostic_totals.residue_phosphorus_g + diagnostic_totals.organic_phosphorus_g + diagnostic_totals.phosphate_phosphorus_g;
                    const diagnostic_baseline_g_p = landscape_mass_balance_state.monitor.?.baseline.phosphorus_g;
                    std.log.debug(
                        "hourly phosphorus closure: hour={d} residual_g={e} storage_change_g={e} atmospheric_input_g={e} runoff_drainage_input_g={e} export_g={e} residue_g={e} organic_g={e} phosphate_g={e}",
                        .{
                            scene_weather_hours + 1,
                            diagnostic_storage_g_p - diagnostic_baseline_g_p - diagnostic_aqueous_input_g_p - diagnostic_input_g_p + diagnostic_export_g_p,
                            diagnostic_storage_g_p - diagnostic_baseline_g_p,
                            diagnostic_aqueous_input_g_p,
                            diagnostic_input_g_p,
                            diagnostic_export_g_p,
                            diagnostic_totals.residue_phosphorus_g,
                            diagnostic_totals.organic_phosphorus_g,
                            diagnostic_totals.phosphate_phosphorus_g,
                        },
                    );
                }
                const soil_process_units = soil_nitrogen_flux_workspace.process_unit_count_per_layer;
                const surface_process_units = surface_microbial_mineral_exchange_state.h2po4_exchange_g_p.len / state.cell_count;
                for (0..state.cell_count) |cell| {
                    var immobilization_g_p: f64 = 0;
                    var immobilization_g_n: f64 = 0;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        const first_unit = layer * soil_process_units;
                        for (first_unit..first_unit + soil_process_units) |unit| immobilization_g_p +=
                            soil_nitrogen_flux_workspace.non_band_microbial_h2po4_exchange_g_p[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_h2po4_exchange_g_p[unit] +
                            soil_nitrogen_flux_workspace.non_band_microbial_hpo4_exchange_g_p[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_hpo4_exchange_g_p[unit];
                        for (first_unit..first_unit + soil_process_units) |unit| immobilization_g_n +=
                            soil_nitrogen_flux_workspace.non_band_microbial_ammonium_exchange_g_n[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_ammonium_exchange_g_n[unit] +
                            soil_nitrogen_flux_workspace.non_band_microbial_nitrate_exchange_g_n[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_nitrate_exchange_g_n[unit];
                    }
                    const surface_first = cell * surface_process_units;
                    for (surface_first..surface_first + surface_process_units) |unit| immobilization_g_p +=
                        surface_microbial_mineral_exchange_state.h2po4_exchange_g_p[unit] +
                        surface_microbial_mineral_exchange_state.hpo4_exchange_g_p[unit];
                    for (surface_first..surface_first + surface_process_units) |unit| immobilization_g_n +=
                        surface_microbial_mineral_exchange_state.ammonium_exchange_g_n[unit] +
                        surface_microbial_mineral_exchange_state.nitrate_exchange_g_n[unit];
                    const next = daily_microbial_phosphate_mineralization_g_p[cell] - immobilization_g_p;
                    if (!std.math.isFinite(next)) return error.NonFiniteDailyMicrobialPhosphateMineralization;
                    daily_microbial_phosphate_mineralization_g_p[cell] = next;
                    const nitrogen_next = daily_microbial_nitrogen_mineralization_g_n[cell] - immobilization_g_n;
                    if (!std.math.isFinite(nitrogen_next)) return error.NonFiniteDailyMicrobialNitrogenMineralization;
                    daily_microbial_nitrogen_mineralization_g_n[cell] = nitrogen_next;
                }
                // Late GROSUB salt removal follows all natural turnover,
                // management harvest, grazing, and tillage publications.
                // This is the single serial pre-EXTRACT ownership boundary.
                if (runscript.dynamic_plant_salts) {
                    const canopy = if (detailed_canopy_state) |*value| value else return error.MissingCanopyForPlantSaltHarvest;
                    const roots = if (plant_root_state) |*value| value else return error.MissingPlantRootsForPlantSaltHarvest;
                    for (0..runtime_plant_count) |plant| {
                        const active_domain_count = root_metabolism_plant_parameters[plant].biologicalDomainCount();
                        try root_litter_carbon_ledger.carbonByPlantLayer(
                            plant,
                            active_domain_count,
                            root_litter_carbon_g_c_by_layer,
                        );
                        var natural_shoot_litter_carbon_g_c: f64 = 0;
                        const natural_litter = shoot_senescence_products_by_plant[plant];
                        for (natural_litter.woody_carbon_g, natural_litter.nonwoody_carbon_g) |woody, nonwoody| {
                            natural_shoot_litter_carbon_g_c += woody + nonwoody;
                            if (!std.math.isFinite(natural_shoot_litter_carbon_g_c))
                                return error.NonFinitePlantSaltLitterCarbon;
                        }
                        try ecosys.plant_salt_harvest_adapter.applyPlant(
                            &plant_salt_harvest_workspace,
                            canopy,
                            roots,
                            .{
                                .cumulative_harvest_salt_mol_by_plant = cumulative_harvest_salt_mol_by_plant,
                                .litter_salt_mol_by_plant_exchange_layer = litter_salt_mol_by_plant_exchange_layer,
                            },
                            plant,
                            .{
                                .shoot_carbon_g_c = canopy.plant_total_shoot_carbon_g[plant],
                                .current_harvest_carbon_g_c = plant_daily_flux_ledger.harvested_carbon_g[plant],
                                .previous_harvest_carbon_g_c = harvest_carbon_at_hour_start_g_c_by_plant[plant],
                                .shoot_litterfall_carbon_g_c = natural_shoot_litter_carbon_g_c +
                                    shoot_harvest_litter_carbon_g_c_by_plant[plant],
                                .root_litterfall_carbon_g_c_by_layer = root_litter_carbon_g_c_by_layer,
                                .active_root_domain_count = active_domain_count,
                                .minimum_plant_mass_g_c = config.negligible_quantity_threshold,
                            },
                        );
                    }
                    try ecosys.litter_salt_publication.refresh(
                        &litter_salt_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .mol_by_plant_exchange_layer_salt = litter_salt_mol_by_plant_exchange_layer,
                        },
                    );
                    try ecosys.plant_litter_salt_ingress.apply(
                        &plant_litter_salt_ingress_state,
                        .{
                            .species_count = config.plant_populations,
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .litter_salt_mol_by_plant_exchange_layer = litter_salt_mol_by_plant_exchange_layer,
                            .surface_water_m3 = surface_precipitation_state.litter_water_m3,
                            .soil_water_m3 = state.matrix_liquid_water_m3,
                            .minimum_water_m3 = config.negligible_quantity_threshold,
                        },
                        &surface_litter_chemistry_state,
                        &initial_chemistry_state,
                    );
                }
                {
                    const checkpoint = try reconstructLandscapeMassBalance(hourly_science_context);
                    std.log.debug("phosphorus after plant salt ingress: hour={d} stored_g={e}", .{ scene_weather_hours + 1, checkpoint.residue_phosphorus_g + checkpoint.organic_phosphorus_g + checkpoint.phosphate_phosphorus_g });
                }
                try ecosys.manure_deposition_publication.refresh(
                    &manure_deposition_publication_state,
                    hourly_manure_products_by_plant,
                );
                {
                    const checkpoint = try reconstructLandscapeMassBalance(hourly_science_context);
                    std.log.debug("phosphorus after manure publication: hour={d} stored_g={e}", .{ scene_weather_hours + 1, checkpoint.residue_phosphorus_g + checkpoint.organic_phosphorus_g + checkpoint.phosphate_phosphorus_g });
                }
                // EXTRACT accumulates these source-signed hourly exchanges
                // for OUTPD before DAY carries balances and clears them.
                @memset(root_soil_element_exchange_workspace, 0);
                if (detailed_canopy_state) |*canopy| {
                    const carbon_exchange = if (canopy_carbon_exchange_state) |*value| value else return error.PlantDailyLedgerRequiresCarbonExchangeState;
                    try ecosys.canopy_ammonia_publication.refresh(
                        &canopy_ammonia_publication_state,
                        .{
                            .plant_branch_offsets = canopy.plant_branch_offsets,
                            .branch_exchange_g_n_per_h = canopy.branch_canopy_ammonia_exchange_g_n_per_h,
                            .preceding_cumulative_exchange_g_n_by_plant = plant_daily_flux_ledger.ammonia_exchange_g_n,
                        },
                    );
                    for (0..runtime_plant_count) |plant| {
                        const branches = try canopy.branchRange(plant);
                        const canopy_fixed_nitrogen_g_n_per_h =
                            try ecosys.canopy_symbiotic_respiration_fixation.sumFixedNitrogenPerHour(
                                canopy.branch_symbiotic_fixed_nitrogen_g_n_per_h[branches.first..branches.end],
                            );
                        rse_canopy_fixation_g_n_per_h_by_plant[plant] = canopy_fixed_nitrogen_g_n_per_h;
                        const carbon_flux = try carbon_exchange.fluxForPlant(canopy, plant);
                        const living_water_source_m3 = if (canopy_surface_exchange_state) |exchange|
                            exchange.transpiration_m3_per_h[plant] + exchange.intercepted_water_change_m3_per_h[plant]
                        else
                            0;
                        var signed_total_respiration_carbon_g = carbon_flux.signed_aboveground_respiration_g_c_per_h;
                        if (plant_root_state) |roots| {
                            for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..config.soil_layers) |layer| {
                                const root = try roots.layerIndex(plant, domain, layer);
                                signed_total_respiration_carbon_g -= roots.actual_respiration_g_c_per_h[root] +
                                    roots.symbiotic_respiration_actual_g_c_per_h[root];
                            };
                        }
                        try plant_daily_flux_ledger.accumulateHourlyExchange(plant, .{
                            .net_canopy_carbon_g = carbon_flux.net_co2_g_c_per_h,
                            .gross_primary_productivity_g = canopy.plant_gross_primary_productivity_g_c_per_h[plant],
                            .signed_total_respiration_carbon_g = signed_total_respiration_carbon_g,
                            .signed_aboveground_respiration_carbon_g = carbon_flux.signed_aboveground_respiration_g_c_per_h,
                            .canopy_and_standing_dead_water_source_m3 = living_water_source_m3 + standing_dead_evaporation_m3_per_h[plant],
                            .canopy_ammonia_exchange_g_n = 0,
                        });
                        if (plant_root_state) |*roots| {
                            var root_soil_carbon_exchange_g: f64 = 0;
                            var root_soil_nitrogen_exchange_g: f64 = 0;
                            var root_soil_phosphorus_exchange_g: f64 = 0;
                            for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..config.soil_layers) |layer| {
                                const root = try roots.layerIndex(plant, domain, layer);
                                for (0..ecosys.plant_root_system.organic_substrate_count) |fraction| {
                                    const substrate =
                                        root * ecosys.plant_root_system.organic_substrate_count +
                                        fraction;
                                    root_soil_carbon_exchange_g +=
                                        roots.exudate_carbon_exchange_g_c_per_h[substrate];
                                    root_soil_nitrogen_exchange_g +=
                                        roots.exudate_nitrogen_exchange_g_n_per_h[substrate];
                                    root_soil_phosphorus_exchange_g +=
                                        roots.exudate_phosphorus_exchange_g_p_per_h[substrate];
                                }
                                rse_phosphate_h2_by_plant[plant] +=
                                    roots.phosphate_h2_uptake_nonband_g_p_per_h[root] +
                                    roots.phosphate_h2_uptake_band_g_p_per_h[root];
                                rse_phosphate_h_by_plant[plant] +=
                                    roots.phosphate_h_uptake_nonband_g_p_per_h[root] +
                                    roots.phosphate_h_uptake_band_g_p_per_h[root];
                            };
                            rse_organic_carbon_by_plant[plant] = root_soil_carbon_exchange_g;
                            rse_organic_nitrogen_by_plant[plant] = root_soil_nitrogen_exchange_g;
                            rse_organic_phosphorus_by_plant[plant] = root_soil_phosphorus_exchange_g;
                            root_soil_nitrogen_exchange_g +=
                                roots.ammonium_uptake_g_n_per_h[plant] +
                                roots.nitrate_uptake_g_n_per_h[plant];
                            root_soil_phosphorus_exchange_g += roots.phosphate_uptake_g_p_per_h[plant];
                            root_soil_carbon_exchange_g_c_per_h_by_plant[plant] =
                                root_soil_carbon_exchange_g;
                            root_soil_nitrogen_exchange_g_n_per_h_by_plant[plant] =
                                root_soil_nitrogen_exchange_g;
                            root_soil_phosphorus_exchange_g_p_per_h_by_plant[plant] =
                                root_soil_phosphorus_exchange_g;
                            try plant_daily_flux_ledger.accumulateHourlyRootSoilExchange(
                                plant,
                                root_soil_carbon_exchange_g,
                                root_soil_nitrogen_exchange_g,
                                root_soil_phosphorus_exchange_g,
                                roots.fixation_uptake_g_n_per_h[plant] + canopy_fixed_nitrogen_g_n_per_h,
                            );
                        } else {
                            try plant_daily_flux_ledger.accumulateHourlyRootSoilExchange(
                                plant,
                                0,
                                0,
                                0,
                                canopy_fixed_nitrogen_g_n_per_h,
                            );
                        }
                        const litterfall = try ecosys.plant_litterfall_publication.totals(
                            shoot_senescence_products_by_plant[plant],
                            root_litter_products_by_plant[plant],
                            .{
                                .carbon_g_c = shoot_harvest_litter_carbon_g_c_by_plant[plant],
                                .nitrogen_g_n = shoot_harvest_litter_nitrogen_g_n_by_plant[plant],
                                .phosphorus_g_p = shoot_harvest_litter_phosphorus_g_p_by_plant[plant],
                            },
                        );
                        aboveground_litter_carbon_g_c_per_h_by_plant[plant] =
                            litterfall.aboveground.carbon_g_c;
                        aboveground_litter_nitrogen_g_n_per_h_by_plant[plant] =
                            litterfall.aboveground.nitrogen_g_n;
                        aboveground_litter_phosphorus_g_p_per_h_by_plant[plant] =
                            litterfall.aboveground.phosphorus_g_p;
                        try plant_daily_flux_ledger.accumulateLitterPublication(
                            plant,
                            litterfall.aboveground.carbon_g_c,
                            litterfall.aboveground.nitrogen_g_n,
                            litterfall.aboveground.phosphorus_g_p,
                            litterfall.belowground.carbon_g_c,
                            litterfall.belowground.nitrogen_g_n,
                            litterfall.belowground.phosphorus_g_p,
                        );
                    }
                    if (plant_root_state) |*roots| {
                        try ecosys.plant_root_soil_exchange_accumulation.accumulate(
                            &plant_root_soil_exchange_state,
                            .{
                                .organic_carbon_exchange_g_c = rse_organic_carbon_by_plant,
                                .organic_nitrogen_exchange_g_n = rse_organic_nitrogen_by_plant,
                                .organic_phosphorus_exchange_g_p = rse_organic_phosphorus_by_plant,
                                .ammonium_uptake_g_n = roots.ammonium_uptake_g_n_per_h,
                                .nitrate_uptake_g_n = roots.nitrate_uptake_g_n_per_h,
                                .phosphate_h2_uptake_g_p = rse_phosphate_h2_by_plant,
                                .phosphate_h_uptake_g_p = rse_phosphate_h_by_plant,
                                .symbiotic_fixation_g_n = roots.fixation_uptake_g_n_per_h,
                                .additional_fixation_g_n = rse_canopy_fixation_g_n_per_h_by_plant,
                                .carbon_balance_g_c = plant_daily_flux_ledger.net_carbon_change_g,
                                .cumulative_respiration_g_c = plant_daily_flux_ledger.signed_total_respiration_carbon_g,
                            },
                        );
                    }
                    try ecosys.grid_cell_litter_standing_dead_publication.refresh(
                        &cell_litter_standing_dead_publication_state,
                        .{
                            .aboveground_litter_carbon_g_c_per_h_by_plant = aboveground_litter_carbon_g_c_per_h_by_plant,
                            .aboveground_litter_nitrogen_g_n_per_h_by_plant = aboveground_litter_nitrogen_g_n_per_h_by_plant,
                            .aboveground_litter_phosphorus_g_p_per_h_by_plant = aboveground_litter_phosphorus_g_p_per_h_by_plant,
                            .standing_dead_carbon_g_c_by_plant = canopy.plant_standing_dead_carbon_g,
                        },
                    );
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.RootSoilElementPublicationRequiresPhenologyState;
                    try ecosys.root_soil_element_exchange_publication.refresh(
                        &root_soil_element_exchange_publication_state,
                        .{
                            .active_by_plant = phenology.active,
                            .exchange_g_element_per_h_by_element_and_plant = .{
                                root_soil_carbon_exchange_g_c_per_h_by_plant,
                                root_soil_nitrogen_exchange_g_n_per_h_by_plant,
                                root_soil_phosphorus_exchange_g_p_per_h_by_plant,
                            },
                        },
                    );
                    @memcpy(
                        plant_daily_flux_ledger.ammonia_exchange_g_n,
                        canopy_ammonia_publication_state
                            .cumulative_exchange_g_n_by_plant,
                    );
                } else {
                    inline for (
                        root_soil_element_exchange_publication_state
                            .change_g_element_per_h_by_element_and_cell,
                    ) |values| @memset(values, 0);
                    @memset(
                        canopy_ammonia_publication_state
                            .current_exchange_g_n_per_h_by_plant,
                        0,
                    );
                    @memcpy(
                        canopy_ammonia_publication_state
                            .cumulative_exchange_g_n_by_plant,
                        plant_daily_flux_ledger.ammonia_exchange_g_n,
                    );
                    inline for (@typeInfo(
                        ecosys.grid_cell_litter_standing_dead_publication.State,
                    ).@"struct".fields) |field|
                        if (field.type == []f64)
                            @memset(
                                @field(
                                    cell_litter_standing_dead_publication_state,
                                    field.name,
                                ),
                                0,
                            );
                }
                @memset(canopy_net_radiation_megajoules, 0);
                @memset(canopy_storage_heat_megajoules, 0);
                if (canopy_energy_state) |energy| if (canopy_surface_exchange_state) |exchange| {
                    for (0..runtime_plant_count) |plant| {
                        const cell = plant / config.plant_populations;
                        canopy_net_radiation_megajoules[plant] = energy.net_radiation_megajoules_per_m2[plant] * canopy_cell_area_m2[cell];
                        // HFLXC is the converged canopy energy residual sum;
                        // EXTRACT later subtracts its convective VFLXC term.
                        canopy_storage_heat_megajoules[plant] = canopy_net_radiation_megajoules[plant] +
                            exchange.latent_heat_flux_megajoules_per_h[plant] +
                            exchange.sensible_heat_flux_megajoules_per_h[plant] +
                            exchange.vapor_sensible_heat_flux_megajoules_per_h[plant];
                    }
                };
                const living_exchange = if (canopy_surface_exchange_state) |*exchange| exchange else null;
                const dead_exchange = if (standing_dead_surface_exchange_state) |*exchange| exchange else null;
                try ecosys.plant_energy_publication.refresh(
                    &plant_energy_publication_state,
                    .{
                        .living_net_radiation_megajoules = canopy_net_radiation_megajoules,
                        .living_latent_heat_megajoules = if (living_exchange) |exchange| exchange.latent_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                        .living_sensible_heat_megajoules = if (living_exchange) |exchange| exchange.sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                        .living_storage_heat_megajoules = canopy_storage_heat_megajoules,
                        .living_convective_water_heat_megajoules = if (living_exchange) |exchange| exchange.vapor_sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                        .standing_dead_net_radiation_megajoules = if (dead_exchange) |exchange| exchange.net_radiation_megajoules_per_h else zero_plant_energy_megajoules,
                        .standing_dead_latent_heat_megajoules = if (dead_exchange) |exchange| exchange.latent_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                        .standing_dead_sensible_heat_megajoules = if (dead_exchange) |exchange| exchange.sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                        .standing_dead_storage_heat_megajoules = if (dead_exchange) |exchange| exchange.storage_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                        .standing_dead_convective_water_heat_megajoules = if (dead_exchange) |exchange| exchange.vapor_sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    },
                );
                if (canopy_precipitation_retention_state) |*retention| {
                    try ecosys.plant_water_publication.refresh(
                        &plant_water_publication_state,
                        .{
                            .cell_area_m2 = canopy_cell_area_m2,
                            .internal_water_depth_m_per_m2_by_plant = plant_state.canopy_water_storage_m_per_m2,
                            .living_surface_water_m3_by_plant = retention.living_surface_water_m3,
                            .standing_dead_surface_water_m3_by_plant = retention.standing_dead_surface_water_m3,
                            .transpiration_m3_per_h_by_plant = if (living_exchange) |exchange|
                                exchange.transpiration_m3_per_h
                            else
                                zero_plant_energy_megajoules,
                            .living_evaporation_m3_per_h_by_plant = if (living_exchange) |exchange|
                                exchange.intercepted_water_change_m3_per_h
                            else
                                zero_plant_energy_megajoules,
                            .standing_dead_evaporation_m3_per_h_by_plant = if (dead_exchange) |exchange|
                                exchange.intercepted_water_change_m3_per_h
                            else
                                zero_plant_energy_megajoules,
                        },
                    );
                }
                try ecosys.ecosystem_energy_ledger.refresh(&ecosystem_energy_ledger_state, .{
                    .cell_area_m2 = canopy_cell_area_m2,
                    .ground_net_radiation_megajoules_per_m2 = surface_energy_state.net_radiation_megajoules_per_m2,
                    .ground_latent_heat_megajoules_per_m2 = surface_temperature_solver_state.latent_heat_flux_megajoules_per_m2,
                    .ground_sensible_heat_megajoules_per_m2 = surface_temperature_solver_state.sensible_heat_flux_megajoules_per_m2,
                    .ground_storage_heat_megajoules_per_m2 = surface_temperature_solver_state.storage_heat_flux_megajoules_per_m2,
                    .species_count = config.plant_populations,
                    .canopy_net_radiation_megajoules = canopy_net_radiation_megajoules,
                    .canopy_latent_heat_megajoules = if (living_exchange) |exchange| exchange.latent_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    .canopy_sensible_heat_megajoules = if (living_exchange) |exchange| exchange.sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    .canopy_storage_heat_megajoules = canopy_storage_heat_megajoules,
                    .canopy_convective_water_heat_megajoules = if (living_exchange) |exchange| exchange.vapor_sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    .standing_dead_net_radiation_megajoules = if (dead_exchange) |exchange| exchange.net_radiation_megajoules_per_h else zero_plant_energy_megajoules,
                    .standing_dead_latent_heat_megajoules = if (dead_exchange) |exchange| exchange.latent_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    .standing_dead_sensible_heat_megajoules = if (dead_exchange) |exchange| exchange.sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    .standing_dead_storage_heat_megajoules = if (dead_exchange) |exchange| exchange.storage_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                    .standing_dead_convective_water_heat_megajoules = if (dead_exchange) |exchange| exchange.vapor_sensible_heat_flux_megajoules_per_h else zero_plant_energy_megajoules,
                });
                if (canopy_precipitation_retention_state) |*retention| {
                    // `uptake.f` closes canopy energy/water, root water uptake,
                    // root gas content, root gas withdrawal, and root nutrient
                    // uptake within one hourly pass. Publishing them from five
                    // separate call sites let a half-advanced hour be observed
                    // by REDIST, NITRO, and the output writers, and the EXTRACT
                    // `ENGYX` carry made a repeated publish non-idempotent. One
                    // transaction snapshots every destination, runs the same
                    // five already validated owners with the same inputs,
                    // checks the cross-owner invariants no single owner can
                    // see, and rolls every array back on any failure.
                    //
                    // Placement is the earliest of the five former sites. The
                    // root gas withdrawal ledger now published ~540 lines early
                    // is unchanged over that span: its writers are
                    // `plant_root_system.resetHourlyFluxes` (line 9618),
                    // `plant_root_disturbance.releaseRootGasFraction` (reached
                    // through harvest and tillage dispatch, all earlier), and
                    // `plant_root_disturbance.withdrawRootAxisLayer` (reached
                    // through `applyRootMetabolism`, line 5565). The root phase
                    // inventories the gas content owner reads are next mutated
                    // by `applyRootFireCombustion`, which is later still.
                    // `plant_root_system.refreshLayerMorphology` is the one
                    // remaining potential writer and is currently unbound; if a
                    // lane binds it between here and the former late site, this
                    // placement must be re-derived.
                    //
                    // Guard equivalence: a non-null retention state implies a
                    // non-null root, phenology, and water workspace state,
                    // because those are constructed together, and the former
                    // root block already returned an error in the mixed case.
                    const uptake_roots = if (plant_root_state) |*value|
                        value
                    else
                        return error.UptakeCoupledTransactionRequiresRootState;
                    const uptake_phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.UptakeCoupledTransactionRequiresPhenologyState;
                    const uptake_water_workspace = if (plant_water_workspace) |*value|
                        value
                    else
                        return error.UptakeCoupledTransactionRequiresPlantWaterWorkspace;
                    try ecosys.uptake_coupled_transaction.apply(
                        &uptake_coupled_transaction_workspace,
                        .{
                            .canopy_energy = &canopy_water_energy_publication_state,
                            .root_water = &root_water_uptake_publication_state,
                            .root_gas_content = &root_gas_content_publication_state,
                            .root_gas_withdrawal = &root_gas_withdrawal_publication_state,
                            .root_nutrient = &root_nutrient_uptake_publication_state,
                        },
                        .{
                            .shape = .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = uptake_phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            },
                            .canopy_energy = .{
                                .air_temperature_k_by_cell = atmospheric_state.air_temperature_k,
                                .canopy_temperature_k_by_plant = plant_state.canopy_temperature_k,
                                .living_surface_water_m3_by_plant = retention.living_surface_water_m3,
                                .standing_dead_surface_water_m3_by_plant = retention.standing_dead_surface_water_m3,
                                .living_retention_m3_per_h_by_plant = retention.living_retention_m3_per_h,
                                .standing_dead_retention_m3_per_h_by_plant = retention.standing_dead_retention_m3_per_h,
                            },
                            .previous_canopy_water_energy_megajoules_by_plant = retention.previous_water_energy_megajoules,
                            .root_water = .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = uptake_phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                                .plant_population_count = uptake_water_workspace.plant_population_count,
                                .cell_area_m2 = canopy_cell_area_m2,
                                .soil_temperature_k_by_layer = state.soil_temperature_k,
                                .root_length_density_m_per_m3 = uptake_roots.root_length_density_m_per_m3,
                                .water_uptake_m3_per_h = uptake_roots.water_uptake_m3_per_h,
                            },
                            .root_gas_content = .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = uptake_phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                                .gaseous_g_by_gas = .{
                                    uptake_roots.gaseous_carbon_dioxide_g_c,
                                    uptake_roots.gaseous_oxygen_g_o,
                                    uptake_roots.gaseous_methane_g_c,
                                    uptake_roots.gaseous_nitrous_oxide_g_n,
                                    uptake_roots.gaseous_ammonia_g_n,
                                    uptake_roots.gaseous_hydrogen_g_h,
                                },
                                .aqueous_g_by_gas = .{
                                    uptake_roots.aqueous_carbon_dioxide_g_c,
                                    uptake_roots.aqueous_oxygen_g_o,
                                    uptake_roots.aqueous_methane_g_c,
                                    uptake_roots.aqueous_nitrous_oxide_g_n,
                                    uptake_roots.aqueous_ammonia_g_n,
                                    uptake_roots.aqueous_hydrogen_g_h,
                                },
                            },
                            .root_gas_withdrawal = .{
                                .loss_g_element_per_h_by_gas_and_plant = .{
                                    uptake_roots.withdrawal_carbon_dioxide_loss_g_c_per_h,
                                    uptake_roots.withdrawal_oxygen_loss_g_o_per_h,
                                    uptake_roots.withdrawal_methane_loss_g_c_per_h,
                                    uptake_roots.withdrawal_nitrous_oxide_loss_g_n_per_h,
                                    uptake_roots.withdrawal_ammonia_loss_g_n_per_h,
                                    uptake_roots.withdrawal_hydrogen_loss_g_h_per_h,
                                },
                            },
                            .root_nutrient = .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = uptake_phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                                .uptake_g_element_per_h_by_nutrient_and_root = .{
                                    uptake_roots.ammonium_uptake_nonband_g_n_per_h,
                                    uptake_roots.nitrate_uptake_nonband_g_n_per_h,
                                    uptake_roots.phosphate_h2_uptake_nonband_g_p_per_h,
                                    uptake_roots.phosphate_h_uptake_nonband_g_p_per_h,
                                    uptake_roots.ammonium_uptake_band_g_n_per_h,
                                    uptake_roots.nitrate_uptake_band_g_n_per_h,
                                    uptake_roots.phosphate_h2_uptake_band_g_p_per_h,
                                    uptake_roots.phosphate_h_uptake_band_g_p_per_h,
                                },
                            },
                        },
                    );
                    @memcpy(
                        ecosystem_energy_ledger_state.canopy_water_energy_megajoules,
                        canopy_water_energy_publication_state
                            .water_energy_megajoules_by_cell,
                    );
                    @memcpy(
                        ecosystem_energy_ledger_state.canopy_water_energy_change_megajoules_per_h,
                        canopy_water_energy_publication_state
                            .water_energy_change_megajoules_per_h_by_cell,
                    );
                } else {
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_megajoules_by_plant,
                        0,
                    );
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_change_megajoules_per_h_by_plant,
                        0,
                    );
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_megajoules_by_cell,
                        0,
                    );
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_change_megajoules_per_h_by_cell,
                        0,
                    );
                    @memset(ecosystem_energy_ledger_state.canopy_water_energy_megajoules, 0);
                    @memset(ecosystem_energy_ledger_state.canopy_water_energy_change_megajoules_per_h, 0);
                }
                if (canopy_layer_distribution_state) |*layers| {
                    const canopy = if (detailed_canopy_state) |*value|
                        value
                    else
                        return error.LivingCanopyPublicationRequiresCanopyState;
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.LivingCanopyPublicationRequiresPhenologyState;
                    try ecosys.living_canopy_layer_publication.refresh(
                        &living_canopy_layer_publication_state,
                        .{
                            .active_by_plant = phenology.active,
                            .plant_branch_offsets = canopy.plant_branch_offsets,
                            .branch_node_offsets = canopy.branch_node_offsets,
                            .node_leaf_area_m2_by_layer = layers.node_leaf_area_m2,
                            .node_leaf_carbon_g_c_by_layer = layers.node_leaf_carbon_g,
                            .branch_stalk_area_m2_by_layer = layers.branch_stalk_area_m2,
                        },
                    );
                    try ecosys.standing_dead_area_publication.refresh(
                        &standing_dead_area_publication_state,
                        layers.plant_standing_dead_area_m2,
                    );
                } else {
                    @memset(
                        living_canopy_layer_publication_state
                            .leaf_area_m2_by_cell_layer,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .leaf_carbon_g_c_by_cell_layer,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .stalk_area_m2_by_cell_layer,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .leaf_area_m2_by_cell,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .stalk_area_m2_by_cell,
                        0,
                    );
                    @memset(
                        standing_dead_area_publication_state
                            .area_m2_by_cell_layer,
                        0,
                    );
                    @memset(
                        standing_dead_area_publication_state.area_m2_by_cell,
                        0,
                    );
                }
                const heat_before = landscape_mass_balance_state.boundary_ledger.cumulative;
                try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedSurfaceAndCanopyHeat(
                    surface_energy_state.net_radiation_megajoules_per_m2,
                    surface_temperature_solver_state.sensible_heat_flux_megajoules_per_m2,
                    surface_temperature_solver_state.latent_heat_flux_megajoules_per_m2,
                    surface_temperature_solver_state.vapor_sensible_heat_flux_megajoules_per_m2,
                    canopy_cell_area_m2,
                    ecosystem_energy_ledger_state.canopy_water_energy_change_megajoules_per_h,
                );
                try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedSurfacePhaseSensibleAdjustment(
                    surface_temperature_solver_state.phase_heat_flux_megajoules_per_m2,
                    surface_temperature_solver_state.ice_water_equivalent_change_m3,
                    state.surface_temperature_k,
                    canopy_cell_area_m2,
                    runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
                    runscript.soil_phase_heat_parameters.ice_heat_capacity_megajoules_per_m3_k,
                );
                // Convective heat carried by water crossing the soil's
                // external (lateral/drainage) boundary. This term was
                // previously computed only inside the `scene_weather_hours < 24`
                // debug block below, so it was never booked into the ledger.
                // Leaving it unbooked under-counts boundary heat and shows up
                // as the `richards` audit deviation.
                {
                    var boundary_water_heat_outward_megajoules: f64 = 0;
                    for (0..state.cell_count) |cell| {
                        const first_layer = cell * state.soil_layer_capacity;
                        for (0..state.active_soil_layer_count[cell]) |local_layer| {
                            const layer = first_layer + local_layer;
                            boundary_water_heat_outward_megajoules += (transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer] + transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]) * runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k * state.soil_temperature_k[layer];
                        }
                    }
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedSignedHeat(
                        -boundary_water_heat_outward_megajoules,
                    );
                }
                if (scene_weather_hours < 24) {
                    var radiation_megajoules: f64 = 0;
                    var sensible_megajoules: f64 = 0;
                    var latent_megajoules: f64 = 0;
                    var vapor_sensible_megajoules: f64 = 0;
                    var phase_megajoules: f64 = 0;
                    var ground_storage_megajoules: f64 = 0;
                    var ground_conduction_megajoules: f64 = 0;
                    var canopy_water_change_megajoules: f64 = 0;
                    var direct_litter_precipitation_heat_megajoules: f64 = 0;
                    var external_soil_water_heat_outward_megajoules: f64 = 0;
                    for (0..state.cell_count) |cell| {
                        const area_m2 = canopy_cell_area_m2[cell];
                        radiation_megajoules += surface_energy_state.net_radiation_megajoules_per_m2[cell] * area_m2;
                        sensible_megajoules += surface_temperature_solver_state.sensible_heat_flux_megajoules_per_m2[cell] * area_m2;
                        latent_megajoules += surface_temperature_solver_state.latent_heat_flux_megajoules_per_m2[cell] * area_m2;
                        vapor_sensible_megajoules += surface_temperature_solver_state.vapor_sensible_heat_flux_megajoules_per_m2[cell] * area_m2;
                        phase_megajoules += surface_temperature_solver_state.phase_heat_flux_megajoules_per_m2[cell] * area_m2;
                        ground_storage_megajoules += surface_temperature_solver_state.storage_heat_flux_megajoules_per_m2[cell] * area_m2;
                        ground_conduction_megajoules += surface_temperature_solver_state.conductive_heat_flux_megajoules_per_m2[cell] * area_m2;
                        canopy_water_change_megajoules += ecosystem_energy_ledger_state.canopy_water_energy_change_megajoules_per_h[cell];
                        const direct_litter_water_m3 = @max(0, surface_precipitation_state.water_to_litter_m3_per_h[cell] - transport_hydrology_state.snow_to_litter_water_flux_m3_per_step[cell]);
                        direct_litter_precipitation_heat_megajoules += direct_litter_water_m3 * runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k * atmospheric_state.air_temperature_k[cell];
                        const first_layer = cell * state.soil_layer_capacity;
                        for (0..state.active_soil_layer_count[cell]) |local_layer| {
                            const layer = first_layer + local_layer;
                            external_soil_water_heat_outward_megajoules += (transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer] + transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]) * runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k * state.soil_temperature_k[layer];
                        }
                    }
                    std.log.debug("current surface/canopy heat boundary: input_megajoules={e} output_megajoules={e} radiation_megajoules={e} sensible_megajoules={e} latent_megajoules={e} vapor_sensible_megajoules={e} phase_megajoules={e} storage_term_megajoules={e} conduction_megajoules={e} canopy_water_change_megajoules={e} direct_litter_precipitation_heat_megajoules={e} external_soil_water_heat_outward_megajoules={e}", .{ landscape_mass_balance_state.boundary_ledger.cumulative.heat_input_megajoules - heat_before.heat_input_megajoules, landscape_mass_balance_state.boundary_ledger.cumulative.heat_output_megajoules - heat_before.heat_output_megajoules, radiation_megajoules, sensible_megajoules, latent_megajoules, vapor_sensible_megajoules, phase_megajoules, ground_storage_megajoules, ground_conduction_megajoules, canopy_water_change_megajoules, direct_litter_precipitation_heat_megajoules, external_soil_water_heat_outward_megajoules });
                }
                if (plant_root_state) |*roots| {
                    const water_workspace = if (plant_water_workspace) |*workspace| workspace else return error.RootUptakeLedgerRequiresPlantWaterWorkspace;
                    try ecosys.root_uptake_ledger.refresh(
                        &root_uptake_ledger_state,
                        &state,
                        roots,
                        config.plant_populations,
                        root_biological_domain_count_by_plant,
                        water_workspace.plant_population_count,
                        canopy_cell_area_m2,
                        runscript.dynamic_plant_salts,
                    );
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.RootWaterPublicationRequiresPhenologyState;
                    // When a canopy precipitation retention state exists, root
                    // water uptake, root gas content, and root nutrient uptake
                    // were already published above by
                    // `uptake_coupled_transaction.apply`, atomically with
                    // canopy energy and root gas withdrawal; publishing them
                    // again here would be a double mutation. These three remain
                    // the owners for the configuration that has roots but no
                    // retention state, which is reachable when no cell carries
                    // a branch, and where the transaction does not run.
                    if (canopy_precipitation_retention_state == null) {
                        try ecosys.root_water_uptake_publication.refresh(
                            &root_water_uptake_publication_state,
                            .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                                .plant_population_count = water_workspace.plant_population_count,
                                .cell_area_m2 = canopy_cell_area_m2,
                                .soil_temperature_k_by_layer = state.soil_temperature_k,
                                .root_length_density_m_per_m3 = roots.root_length_density_m_per_m3,
                                .water_uptake_m3_per_h = roots.water_uptake_m3_per_h,
                            },
                        );
                        try ecosys.root_gas_content_publication.refresh(
                            &root_gas_content_publication_state,
                            .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                                .gaseous_g_by_gas = .{
                                    roots.gaseous_carbon_dioxide_g_c,
                                    roots.gaseous_oxygen_g_o,
                                    roots.gaseous_methane_g_c,
                                    roots.gaseous_nitrous_oxide_g_n,
                                    roots.gaseous_ammonia_g_n,
                                    roots.gaseous_hydrogen_g_h,
                                },
                                .aqueous_g_by_gas = .{
                                    roots.aqueous_carbon_dioxide_g_c,
                                    roots.aqueous_oxygen_g_o,
                                    roots.aqueous_methane_g_c,
                                    roots.aqueous_nitrous_oxide_g_n,
                                    roots.aqueous_ammonia_g_n,
                                    roots.aqueous_hydrogen_g_h,
                                },
                            },
                        );
                        try ecosys.root_nutrient_uptake_publication.refresh(
                            &root_nutrient_uptake_publication_state,
                            .{
                                .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                                .active_by_plant = phenology.active,
                                .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                                .uptake_g_element_per_h_by_nutrient_and_root = .{
                                    roots.ammonium_uptake_nonband_g_n_per_h,
                                    roots.nitrate_uptake_nonband_g_n_per_h,
                                    roots.phosphate_h2_uptake_nonband_g_p_per_h,
                                    roots.phosphate_h_uptake_nonband_g_p_per_h,
                                    roots.ammonium_uptake_band_g_n_per_h,
                                    roots.nitrate_uptake_band_g_n_per_h,
                                    roots.phosphate_h2_uptake_band_g_p_per_h,
                                    roots.phosphate_h_uptake_band_g_p_per_h,
                                },
                            },
                        );
                    }
                    @memcpy(
                        root_uptake_ledger_state
                            .total_root_length_density_m_per_m3,
                        root_water_uptake_publication_state
                            .root_length_density_m_per_m3,
                    );
                    @memcpy(
                        root_uptake_ledger_state.total_water_uptake_m3_per_h,
                        root_water_uptake_publication_state
                            .water_uptake_m3_per_h,
                    );
                    @memcpy(
                        root_uptake_ledger_state.convective_water_heat_megajoules_per_h,
                        root_water_uptake_publication_state
                            .convective_water_heat_megajoules_per_h,
                    );
                    inline for (
                        root_uptake_ledger_state.total_root_gas_content_g,
                        root_gas_content_publication_state
                            .total_g_by_gas_and_layer,
                    ) |destination, source_values| @memcpy(destination, source_values);
                    try ecosys.root_atmosphere_gas_publication.refresh(
                        &root_atmosphere_gas_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .exchange_g_per_h_by_root_and_transport_gas = roots.atmosphere_to_root_gas_exchange_g_per_h,
                        },
                    );
                    inline for (
                        root_uptake_ledger_state
                            .total_atmosphere_to_root_gas_exchange_g_per_h,
                        root_atmosphere_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer,
                    ) |destination, source_values| @memcpy(destination, source_values);
                    try ecosys.root_internal_gas_publication.refresh(
                        &root_internal_gas_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .aqueous_carbon_dioxide_reaction_g_c_per_h_by_root = roots.aqueous_carbon_dioxide_reaction_g_c_per_h,
                            .oxygen_uptake_from_root_pool_g_o_per_h_by_root = roots.oxygen_uptake_from_root_pool_g_o_per_h,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .root_pool_oxygen_uptake_g_o_per_h,
                        root_internal_gas_publication_state
                            .oxygen_uptake_g_o_per_h_by_layer,
                    );
                    try ecosys.root_soil_gas_publication.refresh(
                        &root_soil_gas_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .soil_to_root_exchange_g_per_h_by_root_and_transport_gas = roots.soil_to_root_gas_exchange_g_per_h,
                            .oxygen_uptake_from_soil_g_o_per_h_by_root = roots.oxygen_uptake_from_soil_g_o_per_h,
                        },
                    );
                    for (
                        root_uptake_ledger_state
                            .total_soil_to_root_gas_exchange_g_per_h[0..4],
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer[0..4],
                    ) |destination, source_values| @memcpy(destination, source_values);
                    @memcpy(
                        root_uptake_ledger_state
                            .total_soil_to_root_gas_exchange_g_per_h[5],
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer[4],
                    );
                    @memcpy(
                        root_uptake_ledger_state.soil_oxygen_uptake_g_o_per_h,
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer[1],
                    );
                    for (
                        root_uptake_ledger_state.total_nutrient_uptake_by_layer,
                        root_nutrient_uptake_publication_state
                            .uptake_g_element_per_h_by_nutrient_and_layer,
                    ) |destination, source_values| @memcpy(
                        destination,
                        source_values,
                    );
                    try ecosys.root_salt_uptake_publication.refresh(
                        &root_salt_uptake_publication_state,
                        .{
                            .dynamic_salts_enabled = runscript.dynamic_plant_salts,
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .uptake_mol_per_h_by_root_and_salt = roots.salt_uptake_mol_per_h,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state.total_salt_uptake_mol_per_h,
                        root_salt_uptake_publication_state
                            .uptake_mol_per_h_by_layer_and_salt,
                    );
                    try ecosys.root_soil_ammonia_exchange_publication.refresh(
                        &root_soil_ammonia_exchange_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .non_band_exchange_g_n_per_h_by_root = roots.ammonia_nonband_soil_exchange_g_n_per_h,
                            .band_exchange_g_n_per_h_by_root = roots.ammonia_band_soil_exchange_g_n_per_h,
                        },
                    );
                    for (
                        root_uptake_ledger_state
                            .total_soil_to_root_gas_exchange_g_per_h[4],
                        root_soil_ammonia_exchange_publication_state
                            .non_band_exchange_g_n_per_h_by_layer,
                        root_soil_ammonia_exchange_publication_state
                            .band_exchange_g_n_per_h_by_layer,
                    ) |*dest, nonband, band| dest.* = nonband + band;
                    try ecosys.root_exudate_publication.refresh(
                        &root_exudate_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .exchange_g_element_per_h_by_element_root_and_fraction = .{
                                roots.exudate_carbon_exchange_g_c_per_h,
                                roots.exudate_nitrogen_exchange_g_n_per_h,
                                roots.exudate_phosphorus_exchange_g_p_per_h,
                            },
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_exudate_carbon_change_g_c_per_h,
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction[0],
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_exudate_nitrogen_change_g_n_per_h,
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction[1],
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_exudate_phosphorus_change_g_p_per_h,
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction[2],
                    );
                    try ecosys.root_competition_demand_publication.refresh(
                        &root_competition_demand_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .demand_g_element_per_h_by_kind_and_root = .{
                                roots.oxygen_demand_g_o_per_h,
                                roots.ammonium_demand_nonband_g_n_per_h,
                                roots.nitrate_demand_nonband_g_n_per_h,
                                roots.phosphate_h2_demand_nonband_g_p_per_h,
                                roots.phosphate_h_demand_nonband_g_p_per_h,
                                roots.ammonium_demand_band_g_n_per_h,
                                roots.nitrate_demand_band_g_n_per_h,
                                roots.phosphate_h2_demand_band_g_p_per_h,
                                roots.phosphate_h_demand_band_g_p_per_h,
                            },
                        },
                    );
                    for (
                        root_uptake_ledger_state
                            .total_competition_demand_by_layer,
                        root_competition_demand_publication_state
                            .demand_g_element_per_h_by_kind_and_layer,
                    ) |destination, source_values| @memcpy(
                        destination,
                        source_values,
                    );
                    try ecosys.root_nitrogen_fixation_publication.refresh(
                        &root_nitrogen_fixation_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .fixation_g_n_per_h_by_root = roots.fixation_uptake_g_n_per_h_by_layer,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_nitrogen_fixation_g_n_per_h,
                        root_nitrogen_fixation_publication_state
                            .fixation_g_n_per_h_by_layer,
                    );
                } else {
                    @memset(
                        root_water_uptake_publication_state
                            .root_length_density_m_per_m3,
                        0,
                    );
                    @memset(
                        root_water_uptake_publication_state
                            .water_uptake_m3_per_h,
                        0,
                    );
                    @memset(
                        root_water_uptake_publication_state
                            .convective_water_heat_megajoules_per_h,
                        0,
                    );
                    inline for (
                        root_gas_content_publication_state
                            .total_g_by_gas_and_layer,
                    ) |values| @memset(values, 0);
                    inline for (
                        root_atmosphere_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer,
                    ) |values| @memset(values, 0);
                    @memset(
                        root_internal_gas_publication_state
                            .carbon_dioxide_source_g_c_per_h_by_layer,
                        0,
                    );
                    inline for (
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer,
                    ) |values| @memset(values, 0);
                    inline for (
                        root_nutrient_uptake_publication_state
                            .uptake_g_element_per_h_by_nutrient_and_layer,
                    ) |values| @memset(values, 0);
                    @memset(
                        root_salt_uptake_publication_state
                            .uptake_mol_per_h_by_layer_and_salt,
                        0,
                    );
                    @memset(
                        root_soil_ammonia_exchange_publication_state
                            .non_band_exchange_g_n_per_h_by_layer,
                        0,
                    );
                    @memset(
                        root_soil_ammonia_exchange_publication_state
                            .band_exchange_g_n_per_h_by_layer,
                        0,
                    );
                    inline for (
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction,
                    ) |values| @memset(values, 0);
                    inline for (
                        root_competition_demand_publication_state
                            .demand_g_element_per_h_by_kind_and_layer,
                    ) |values| @memset(values, 0);
                    @memset(
                        root_nitrogen_fixation_publication_state
                            .fixation_g_n_per_h_by_layer,
                        0,
                    );
                    @memset(
                        root_internal_gas_publication_state
                            .oxygen_uptake_g_o_per_h_by_layer,
                        0,
                    );
                    @memset(root_uptake_ledger_state.total_root_length_density_m_per_m3, 0);
                    @memset(root_uptake_ledger_state.total_water_uptake_m3_per_h, 0);
                    @memset(root_uptake_ledger_state.convective_water_heat_megajoules_per_h, 0);
                    inline for (root_uptake_ledger_state.total_root_gas_content_g) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_soil_to_root_gas_exchange_g_per_h) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_aqueous_to_gaseous_root_exchange_g_per_h) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_atmosphere_to_root_gas_exchange_g_per_h) |values| @memset(values, 0);
                    @memset(root_uptake_ledger_state.total_salt_uptake_mol_per_h, 0);
                    @memset(root_uptake_ledger_state.total_exudate_carbon_change_g_c_per_h, 0);
                    @memset(root_uptake_ledger_state.total_exudate_nitrogen_change_g_n_per_h, 0);
                    @memset(root_uptake_ledger_state.total_exudate_phosphorus_change_g_p_per_h, 0);
                    inline for (root_uptake_ledger_state.total_competition_demand_by_layer) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_nutrient_uptake_by_layer) |values| @memset(values, 0);
                    @memset(root_uptake_ledger_state.soil_oxygen_uptake_g_o_per_h, 0);
                    @memset(root_uptake_ledger_state.root_pool_oxygen_uptake_g_o_per_h, 0);
                    @memset(root_uptake_ledger_state.total_nitrogen_fixation_g_n_per_h, 0);
                }
                var soil_fire = runscript.soil_fire_combustion_parameters;
                {
                    const checkpoint = try reconstructLandscapeMassBalance(hourly_science_context);
                    std.log.debug("phosphorus before soil fire: hour={d} stored_g={e}", .{ scene_weather_hours + 1, checkpoint.residue_phosphorus_g + checkpoint.organic_phosphorus_g + checkpoint.phosphate_phosphorus_g });
                }
                soil_fire.negligible_carbon_g_c = config.negligible_quantity_threshold;
                const fire_product_parameters: ecosys.plant_soil_exchange.SubsurfaceFireParameters = .{
                    .oxygen_g_per_g_carbon = soil_fire.oxygen_g_per_g_combusted_carbon,
                    .maximum_aerobic_charcoal_fraction = soil_fire.maximum_aerobic_charcoal_fraction,
                    .maximum_anaerobic_charcoal_fraction = soil_fire.maximum_anaerobic_charcoal_fraction,
                    .oxygen_half_saturation_g_o_per_m3 = soil_fire.oxygen_half_saturation_g_o_per_m3,
                    .methane_half_saturation_g_c_per_m3 = soil_fire.methane_half_saturation_g_c_per_m3,
                    .aerobic_combustion_energy_megajoules_per_g_carbon = soil_fire.aerobic_combustion_energy_megajoules_per_g_c,
                    .anaerobic_combustion_energy_megajoules_per_g_carbon = soil_fire.anaerobic_combustion_energy_megajoules_per_g_c,
                    .methane_combustion_energy_megajoules_per_g_carbon = soil_fire.methane_combustion_energy_megajoules_per_g_c,
                };
                for (0..state.cell_count) |cell| {
                    _ = ecosys.soil_combustion.burnSurfaceOrganicStateCell(
                        &surface_organic_state,
                        &surface_fire_exchange_state,
                        cell,
                        fire_active_this_hour[cell],
                        state.surface_temperature_k[cell],
                        canopy_cell_area_m2[cell],
                        1,
                        soil_fire,
                    ) catch |err| {
                        std.log.err("surface organic matter combustion failed: scene={d} scene_hour={d} cell={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, @errorName(err) });
                        return err;
                    };
                    surface_fire_exchange_state.finalizeSurfaceCell(
                        cell,
                        &litter_gas_transport_state,
                        &surface_litter_chemistry_state,
                        &surface_organic_state,
                        surface_precipitation_state.litter_water_m3,
                        runscript.dynamic_plant_salts,
                        delayed_surface_combustion_heat_megajoules,
                        config.negligible_quantity_threshold,
                        runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                        fire_product_parameters,
                    ) catch |err| {
                        std.log.err("surface combustion products failed: scene={d} scene_hour={d} cell={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, @errorName(err) });
                        return err;
                    };
                }
                // Fire mineral products enter the surface solution above.
                // Re-converge only when fire changed that solution; repeating
                // the nonlinear projection without a source is not idempotent
                // and would manufacture phosphate from stale carrier totals.
                // Never repeat the full hourly or legacy sub-hourly cycle.
                const any_surface_fire_active =
                    std.mem.indexOfScalar(bool, fire_active_this_hour, true) != null;
                if (any_surface_fire_active)
                    try convergeSurfaceLitterChemistry(hourly_science_context);
                {
                    const checkpoint = try reconstructLandscapeMassBalance(hourly_science_context);
                    std.log.debug("phosphorus after soil fire chemistry: hour={d} stored_g={e}", .{ scene_weather_hours + 1, checkpoint.residue_phosphorus_g + checkpoint.organic_phosphorus_g + checkpoint.phosphate_phosphorus_g });
                }
                for (0..state.cell_count) |cell| {
                    if (!fire_active_this_hour[cell]) continue;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        _ = ecosys.soil_combustion.burnOrganicStateLayer(
                            &soil_organic_state,
                            &organic_matter_fire_exchange_state,
                            layer,
                            true,
                            state.soil_temperature_k[layer],
                            canopy_cell_area_m2[cell],
                            1,
                            soil_fire,
                        ) catch |err| {
                            std.log.err("soil organic matter combustion failed: scene={d} scene_hour={d} cell={d} layer={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, local_layer, @errorName(err) });
                            return err;
                        };
                    }
                }
                if (plant_root_state) |*roots| {
                    ecosys.disturbance_management_dispatch.applyRootFireCombustion(
                        roots,
                        if (detailed_canopy_state) |*canopy| canopy else null,
                        &state,
                        config.plant_populations,
                        root_biological_domain_count_by_plant,
                        canopy_cell_area_m2,
                        fire_active_this_hour,
                        1,
                        runscript.dynamic_plant_salts,
                        runscript.plant_fire_combustion_parameters,
                        &organic_matter_fire_exchange_state,
                    ) catch |err| {
                        std.log.err("root nodule combustion failed: scene={d} scene_hour={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, @errorName(err) });
                        return err;
                    };
                    try ecosys.root_combustion_salt_publication.refresh(
                        &root_combustion_salt_publication_state,
                        .{
                            .dynamic_salts_enabled = runscript.dynamic_plant_salts,
                            .fire_active_by_cell = fire_active_this_hour,
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .released_mol_per_h_by_root_and_salt = roots.combustion_salt_loss_mol_per_h,
                        },
                    );
                } else {
                    @memset(
                        root_combustion_salt_publication_state
                            .released_mol_per_h_by_layer_and_salt,
                        0,
                    );
                    @memset(
                        root_combustion_salt_publication_state
                            .total_released_mol_per_h_by_cell,
                        0,
                    );
                }
                // Soil organic fire products exist independently of whether
                // the run has plant roots; roots merely add to this ledger.
                for (0..state.cell_count) |cell| {
                    if (!fire_active_this_hour[cell]) continue;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        organic_matter_fire_exchange_state.finalizeLayer(
                            layer,
                            &gas_transport_state,
                            &plant_available_nutrient_state,
                            &soil_organic_state,
                            &micropore_solute_state,
                            runscript.dynamic_plant_salts,
                            delayed_subsurface_combustion_heat_megajoules,
                            config.negligible_quantity_threshold,
                            fire_product_parameters,
                        ) catch |err| {
                            std.log.err("subsurface combustion products failed: scene={d} scene_hour={d} cell={d} layer={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, local_layer, @errorName(err) });
                            return err;
                        };
                    }
                    var emitted_carbon_dioxide_g_c =
                        surface_fire_exchange_state.carbon_dioxide_emission_g_c[cell];
                    var emitted_methane_g_c =
                        surface_fire_exchange_state.methane_emission_g_c[cell];
                    var produced_charcoal_g_c =
                        surface_fire_exchange_state.charcoal_production_g_c[cell];
                    var emitted_phosphorus_g_p = surface_fire_exchange_state.gaseous_phosphorus_emission_g_p[cell];
                    var combusted_phosphorus_g_p = surface_fire_exchange_state.combusted_phosphorus_g_p[cell];
                    var emitted_nitrogen_g_n = surface_fire_exchange_state.gaseous_nitrogen_emission_g_n[cell];
                    var combusted_nitrogen_g_n = surface_fire_exchange_state.combusted_nitrogen_g_n[cell];
                    var ammonium_production_g_n = surface_fire_exchange_state.ammonium_production_g_n[cell];
                    var phosphate_production_g_p = surface_fire_exchange_state.phosphate_production_g_p[cell];
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        emitted_carbon_dioxide_g_c +=
                            organic_matter_fire_exchange_state.carbon_dioxide_emission_g_c[layer];
                        emitted_methane_g_c +=
                            organic_matter_fire_exchange_state.methane_emission_g_c[layer];
                        produced_charcoal_g_c +=
                            organic_matter_fire_exchange_state.charcoal_production_g_c[layer];
                        emitted_phosphorus_g_p += organic_matter_fire_exchange_state.gaseous_phosphorus_emission_g_p[layer];
                        combusted_phosphorus_g_p += organic_matter_fire_exchange_state.combusted_phosphorus_g_p[layer];
                        emitted_nitrogen_g_n += organic_matter_fire_exchange_state.gaseous_nitrogen_emission_g_n[layer];
                        combusted_nitrogen_g_n += organic_matter_fire_exchange_state.combusted_nitrogen_g_n[layer];
                        ammonium_production_g_n += organic_matter_fire_exchange_state.ammonium_production_g_n[layer];
                        phosphate_production_g_p += organic_matter_fire_exchange_state.phosphate_production_g_p[layer];
                    }
                    var root_ratio: f64 = 0;
                    var root_combustion_carbon_abs: f64 = 0;
                    if (plant_root_state) |*roots| {
                        const first_plant = cell * config.plant_populations;
                        const last_plant = try std.math.add(
                            usize,
                            first_plant,
                            config.plant_populations,
                        );
                        for (first_plant..last_plant) |plant| {
                            const carbon_loss = roots.combustion_carbon_loss_g_c_per_h[plant];
                            if (carbon_loss > 0) return error.InvalidSoilFireCarbonEmission;
                            root_combustion_carbon_abs += -carbon_loss;
                        }
                        const total_carbon_g_c = emitted_carbon_dioxide_g_c +
                            emitted_methane_g_c +
                            produced_charcoal_g_c;
                        if (total_carbon_g_c > 0 and root_combustion_carbon_abs > 0 and std.math.isFinite(total_carbon_g_c)) {
                            root_ratio = @min(
                                1.0,
                                root_combustion_carbon_abs / total_carbon_g_c,
                            );
                        }
                    }
                    if (!std.math.isFinite(emitted_carbon_dioxide_g_c) or emitted_carbon_dioxide_g_c < 0 or
                        !std.math.isFinite(emitted_methane_g_c) or emitted_methane_g_c < 0)
                        return error.InvalidSoilFireCarbonEmission;
                    if (!std.math.isFinite(produced_charcoal_g_c) or produced_charcoal_g_c < 0) return error.InvalidSoilFireCharcoalProduction;
                    if (!std.math.isFinite(emitted_phosphorus_g_p) or emitted_phosphorus_g_p < 0) return error.InvalidSoilFirePhosphorusEmission;
                    if (!std.math.isFinite(combusted_phosphorus_g_p) or combusted_phosphorus_g_p < 0) return error.InvalidSoilFirePhosphorusCombustion;
                    if (!std.math.isFinite(emitted_nitrogen_g_n) or emitted_nitrogen_g_n < 0) return error.InvalidSoilFireNitrogenEmission;
                    if (!std.math.isFinite(combusted_nitrogen_g_n) or combusted_nitrogen_g_n < 0) return error.InvalidSoilFireNitrogenCombustion;
                    // REDIST fire ledgers retain their source signs:
                    // VCO2G/VCH4G are negative emissions, while VCOXFS is
                    // positive charcoal transferred back to soil organic C.
                    const root_carbon_dioxide_emission_g_c =
                        emitted_carbon_dioxide_g_c * root_ratio;
                    const root_methane_emission_g_c =
                        emitted_methane_g_c * root_ratio;
                    const root_charcoal_production_g_c =
                        produced_charcoal_g_c * root_ratio;
                    const root_nitrogen_emission_g_n =
                        emitted_nitrogen_g_n * root_ratio;
                    const root_phosphorus_emission_g_p =
                        emitted_phosphorus_g_p * root_ratio;
                    const root_ammonium_g_n =
                        ammonium_production_g_n * root_ratio;
                    const root_dihydrogen_phosphate_g_p =
                        phosphate_production_g_p * root_ratio;
                    const root_combusted_nitrogen_g_n =
                        combusted_nitrogen_g_n * root_ratio;
                    const root_combusted_phosphorus_g_p =
                        combusted_phosphorus_g_p * root_ratio;
                    daily_soil_fire_carbon_dioxide_emission_g_c[cell] -= emitted_carbon_dioxide_g_c - root_carbon_dioxide_emission_g_c;
                    daily_soil_fire_methane_emission_g_c[cell] -= emitted_methane_g_c - root_methane_emission_g_c;
                    daily_soil_fire_charcoal_production_g_c[cell] += produced_charcoal_g_c - root_charcoal_production_g_c;
                    daily_soil_fire_phosphorus_flux_g_p[cell] -= emitted_phosphorus_g_p - root_phosphorus_emission_g_p;
                    daily_soil_combusted_phosphorus_g_p[cell] -= combusted_phosphorus_g_p - root_combusted_phosphorus_g_p;
                    daily_soil_fire_nitrogen_flux_g_n[cell] -= emitted_nitrogen_g_n - root_nitrogen_emission_g_n;
                    daily_soil_combusted_nitrogen_g_n[cell] -= combusted_nitrogen_g_n - root_combusted_nitrogen_g_n;
                    daily_root_fire_carbon_dioxide_emission_g_c[cell] -= root_carbon_dioxide_emission_g_c;
                    daily_root_fire_methane_emission_g_c[cell] -= root_methane_emission_g_c;
                    daily_root_fire_nitrogen_flux_g_n[cell] -= root_nitrogen_emission_g_n;
                    daily_root_fire_phosphorus_flux_g_p[cell] -= root_phosphorus_emission_g_p;
                    daily_root_fire_nitrogen_loss_g_n[cell] -= root_combusted_nitrogen_g_n;
                    daily_root_fire_phosphorus_loss_g_p[cell] -= root_combusted_phosphorus_g_p;
                    if (root_combustion_carbon_abs > 0) {
                        try ecosys.root_combustion_boundary_publication.publish(
                            &root_combustion_boundary_ledger,
                            .{
                                .carbon_dioxide_g_c = root_carbon_dioxide_emission_g_c,
                                .methane_g_c = root_methane_emission_g_c,
                                .charcoal_g_c = root_charcoal_production_g_c,
                                .nitrogen_oxide_g_n = root_nitrogen_emission_g_n,
                                .ammonium_g_n = root_ammonium_g_n,
                                .phosphorus_oxide_g_p = root_phosphorus_emission_g_p,
                                .dihydrogen_phosphate_g_p = root_dihydrogen_phosphate_g_p,
                            },
                        );
                    }
                }
                if (plant_root_state) |*roots| {
                    const any_fire_active = std.mem.indexOfScalar(bool, fire_active_this_hour, true) != null;
                    if (any_fire_active) if (detailed_canopy_state) |*canopy| {
                        const canopy_fire_gas_parameters = surface_gas_parameters;
                        ecosys.plant_shoot_fire.apply(
                            canopy,
                            &canopy_layer_controls,
                            roots,
                            canopy_cell_area_m2,
                            fire_active_this_hour,
                            1,
                            runscript.dynamic_plant_salts,
                            runscript.plant_fire_combustion_parameters,
                            ground_air_state.temperature_k,
                            ground_air_state.air_volume_m3,
                            site.atmospheric_oxygen_umol_mol,
                            site.atmospheric_methane_umol_mol,
                            canopy_fire_gas_parameters.atmospheric_concentration_g_per_m3[@intFromEnum(ecosys.gas_transport.Species.oxygen)],
                            delayed_live_canopy_combustion_heat_megajoules,
                            delayed_standing_dead_combustion_heat_megajoules,
                            &surface_organic_state,
                        ) catch |err| {
                            std.log.err("shoot combustion failed: scene={d} scene_hour={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, @errorName(err) });
                            return err;
                        };
                    };
                    if (detailed_canopy_state) |*canopy| {
                        try ecosys.plant_combustion_publication.refresh(
                            &plant_combustion_publication_state,
                            .{
                                .plant_species_per_cell = config.plant_populations,
                                .fire_active_by_cell = fire_active_this_hour,
                                .signed_carbon_loss_g_c_per_h_by_plant = canopy.plant_combustion_carbon_loss_g_c_per_h,
                                .charcoal_return_g_c_per_h_by_plant = canopy.plant_fire_charcoal_production_g_c_per_h,
                                .signed_nitrogen_loss_g_n_per_h_by_plant = canopy.plant_combustion_nitrogen_loss_g_n_per_h,
                                .signed_phosphorus_loss_g_p_per_h_by_plant = canopy.plant_combustion_phosphorus_loss_g_p_per_h,
                            },
                        );
                        for (0..runtime_plant_count) |plant| {
                            const cell = plant / config.plant_populations;
                            if (!fire_active_this_hour[cell]) continue;
                            try plant_daily_flux_ledger.accumulateOxidation(
                                plant,
                                plant_combustion_publication_state.signed_carbon_balance_g_c_per_h_by_plant[plant],
                                plant_combustion_publication_state.signed_nitrogen_loss_g_n_per_h_by_plant[plant],
                                plant_combustion_publication_state.signed_phosphorus_loss_g_p_per_h_by_plant[plant],
                            );
                        }
                    }
                }
                if (detailed_canopy_state) |*canopy| {
                    try ecosys.canopy_fire_publication.refresh(
                        &canopy_fire_publication_state,
                        .{
                            .fire_active_by_cell = fire_active_this_hour,
                            .carbon_dioxide_emission_g_c_per_h_by_plant = canopy.plant_fire_carbon_dioxide_emission_g_c_per_h,
                            .methane_emission_g_c_per_h_by_plant = canopy.plant_fire_methane_emission_g_c_per_h,
                            .oxygen_consumption_g_o_per_h_by_plant = canopy.plant_fire_oxygen_consumption_g_o_per_h,
                            .charcoal_production_g_c_per_h_by_plant = canopy.plant_fire_charcoal_production_g_c_per_h,
                            .heat_release_megajoules_per_h_by_plant = canopy.plant_fire_heat_release_megajoules_per_h,
                        },
                    );
                } else {
                    inline for (@typeInfo(
                        ecosys.canopy_fire_publication.State,
                    ).@"struct".fields) |field|
                        if (field.type == []f64)
                            @memset(
                                @field(canopy_fire_publication_state, field.name),
                                0,
                            );
                }
                // The coupled UPTAKE transaction above already published this
                // owner, atomically with the other four, whenever a canopy
                // precipitation retention state exists. This remains the owner
                // for the configuration that has roots but no retention state,
                // where the transaction does not run at all; publishing here in
                // the retention case as well would advance the ledger twice.
                if (plant_root_state) |*roots| {
                    if (canopy_precipitation_retention_state == null)
                        try ecosys.root_gas_withdrawal_publication.refresh(
                            &root_gas_withdrawal_publication_state,
                            .{
                                .loss_g_element_per_h_by_gas_and_plant = .{
                                    roots.withdrawal_carbon_dioxide_loss_g_c_per_h,
                                    roots.withdrawal_oxygen_loss_g_o_per_h,
                                    roots.withdrawal_methane_loss_g_c_per_h,
                                    roots.withdrawal_nitrous_oxide_loss_g_n_per_h,
                                    roots.withdrawal_ammonia_loss_g_n_per_h,
                                    roots.withdrawal_hydrogen_loss_g_h_per_h,
                                },
                            },
                        );
                } else {
                    inline for (
                        root_gas_withdrawal_publication_state
                            .loss_g_element_per_h_by_gas_and_cell,
                    ) |values| @memset(values, 0);
                }
                if (detailed_canopy_state) |*canopy| {
                    const exchange = if (canopy_carbon_exchange_state) |*value|
                        value
                    else
                        return error.CanopyCarbonPublicationRequiresExchangeState;
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.CanopyCarbonPublicationRequiresPhenologyState;
                    for (0..runtime_plant_count) |plant|
                        canopy_net_fixation_g_c_per_h_by_plant[plant] =
                            (try exchange.fluxForPlant(canopy, plant))
                                .net_co2_g_c_per_h;
                    try ecosys.canopy_carbon_publication.refresh(
                        &canopy_carbon_publication_state,
                        .{
                            .active_by_plant = phenology.active,
                            .net_fixation_g_c_per_h_by_plant = canopy_net_fixation_g_c_per_h_by_plant,
                            .oxygen_g_o_per_g_c = runscript.root_gas_parameters
                                .oxygen_to_carbon_respiration_ratio_g_o_per_g_c,
                        },
                    );
                } else {
                    @memset(
                        canopy_carbon_publication_state
                            .hourly_net_fixation_g_c_per_h_by_cell,
                        0,
                    );
                    @memset(
                        canopy_carbon_publication_state
                            .canopy_carbon_exchange_g_c_per_h_by_cell,
                        0,
                    );
                    @memset(
                        canopy_carbon_publication_state
                            .canopy_oxygen_exchange_g_o_per_h_by_cell,
                        0,
                    );
                }
                for (0..state.cell_count) |cell| {
                    const canopy_net_fixation_g_c = if (canopy_carbon_exchange_state) |ledger|
                        try ledger.netFixationForCell(&detailed_canopy_state.?, cell)
                    else
                        0;
                    try daily_canopy_gas_exchange.accumulateHour(
                        cell,
                        canopy_net_fixation_g_c,
                        canopy_fire_publication_state
                            .carbon_dioxide_emission_g_c_per_h_by_cell[cell],
                        canopy_fire_publication_state
                            .methane_emission_g_c_per_h_by_cell[cell],
                        canopy_fire_publication_state
                            .oxygen_consumption_g_o_per_h_by_cell[cell],
                        runscript.root_gas_parameters.oxygen_to_carbon_respiration_ratio_g_o_per_g_c,
                    );
                    // EXTRACT lines 939-943, REDIST line 10991: add plant photosynthesis
                    // contribution to hourly XCNET/XONET (CNET per cell from canopy carbon).
                    hourly_accumulators[cell].canopy_co2_exchange_g_c_timestep +=
                        canopy_carbon_publication_state.canopy_carbon_exchange_g_c_per_h_by_cell[cell];
                    hourly_accumulators[cell].canopy_o2_exchange_g_o_timestep +=
                        canopy_carbon_publication_state.canopy_oxygen_exchange_g_o_per_h_by_cell[cell];
                }
                const completed_hour_count = try std.math.add(usize, scene_weather_hours, 1);
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 0;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilCarbonEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const root_gas_withdrawal = if (plant_root_state) |*roots|
                                try ecosys.plant_root_disturbance.rootGasWithdrawalForCell(roots, cell, config.plant_populations)
                            else
                                ecosys.plant_root_disturbance.CellRootGasWithdrawal{};
                            const values = try hourly_soil_carbon_bank.row(cell);
                            const carbon_dioxide_concentration = values[4 .. 4 + config.soil_layers];
                            const methane_first = 5 + config.soil_layers;
                            const methane_concentration = values[methane_first .. methane_first + config.soil_layers];
                            const oxygen_first = methane_first + config.soil_layers;
                            const oxygen_concentration = values[oxygen_first .. oxygen_first + config.soil_layers];
                            for (0..state.active_soil_layer_count[cell]) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const water_m3 = state.liquid_water_m3[layer];
                                inline for (.{ ecosys.gas_transport.Species.carbon_dioxide, ecosys.gas_transport.Species.methane, ecosys.gas_transport.Species.oxygen }) |species| {
                                    const component = layer * ecosys.gas_transport.species_count + @intFromEnum(species);
                                    const concentration = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(gas_transport_state.dissolved_mass_g[component], water_m3);
                                    switch (species) {
                                        .carbon_dioxide => carbon_dioxide_concentration[local_layer] = concentration,
                                        .methane => methane_concentration[local_layer] = concentration,
                                        .oxygen => oxygen_concentration[local_layer] = concentration,
                                        else => unreachable,
                                    }
                                }
                            }
                            var ground_flux_g: [ecosys.gas_transport.species_count]f64 = @splat(0);
                            inline for (.{ ecosys.gas_transport.Species.carbon_dioxide, ecosys.gas_transport.Species.methane, ecosys.gas_transport.Species.oxygen }) |species|
                                ground_flux_g[@intFromEnum(species)] = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                    soil_gas_transport_state.atmospheric_flux_g_per_h,
                                    surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                    try state.layerIndex(cell, 0),
                                    config.soil_layers,
                                    cell,
                                    species,
                                );
                            const litter_oxygen = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(
                                litter_gas_transport_state.dissolved_mass_g[cell * ecosys.gas_transport.species_count + @intFromEnum(ecosys.gas_transport.Species.oxygen)],
                                surface_litter_gas_transport_state.water_volume_m3[cell],
                            );
                            const canopy_net_fixation_g_c_per_h = if (canopy_carbon_exchange_state) |ledger|
                                try ledger.netFixationForCell(&detailed_canopy_state.?, cell)
                            else
                                0;
                            const root_atmospheric_co2_g_c_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 0);
                            const root_atmospheric_oxygen_g_o_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 1);
                            const root_atmospheric_methane_g_c_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 2);
                            try ecosys.soil_biogeochemistry_output.calculateCarbonInto(.{
                                .carbon_dioxide_emission_g_c_per_h = ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] + root_gas_withdrawal.carbon_dioxide_g_c_per_h - root_atmospheric_co2_g_c_per_h,
                                .net_carbon_exchange_g_c_per_h = canopy_net_fixation_g_c_per_h + ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] + root_gas_withdrawal.carbon_dioxide_g_c_per_h - root_atmospheric_co2_g_c_per_h,
                                .methane_emission_g_c_per_h = ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.methane)] + root_gas_withdrawal.methane_g_c_per_h - root_atmospheric_methane_g_c_per_h,
                                .oxygen_exchange_g_o2_per_h = ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.oxygen)] + root_gas_withdrawal.oxygen_g_o_per_h - root_atmospheric_oxygen_g_o_per_h,
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .carbon_dioxide_concentration_by_layer = carbon_dioxide_concentration,
                                .canopy_air_carbon_dioxide_umol_per_mol = current_atmospheric_co2_umol_per_mol,
                                .methane_concentration_by_layer = methane_concentration,
                                .oxygen_concentration_by_layer = oxygen_concentration,
                                .litter_oxygen_concentration = litter_oxygen,
                            }, values);
                            const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_carbon_bank.streams[cell].write(file_name, hourly_soil_carbon_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                // The file name keeps the FOUTS grid-index stem,
                                // but the row itself carries the physical site
                                // coordinates the user entered, so a row is
                                // locatable without knowing the grid layout.
                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 0;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantCarbonEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresCanopyState;
                            const exchange = if (canopy_carbon_exchange_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresCarbonExchangeState;
                            const surface = if (canopy_surface_exchange_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresSurfaceExchangeState;
                            const surface_workspace = if (canopy_surface_input_workspace) |*value| value else return error.HourlyPlantCarbonOutputRequiresSurfaceWorkspace;
                            const layers = if (canopy_layer_distribution_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresLayerDistribution;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const flux = try exchange.fluxForPlant(canopy, plant);
                                const output = try ecosys.plant_hourly_output.carbon(
                                    flux.net_co2_g_c_per_h,
                                    canopy.plant_gross_primary_productivity_g_c_per_h[plant],
                                    flux.signed_aboveground_respiration_g_c_per_h,
                                    canopy.plant_mobile_carbon_concentration_g_per_g[plant],
                                    surface_workspace.stomatal_resistance_h_per_m[plant],
                                    surface.boundary_layer_resistance_h_per_m[plant],
                                    try layers.plantLeafAreaM2(canopy, plant),
                                    canopy_cell_area_m2[cell],
                                );
                                const values = try hourly_plant_carbon_bank.row(plant);
                                values[0..7].* = output.values();
                                // Only assigned species produce output; an unassigned population
                                // would emit all-zero rows, so no file is created for it.
                                const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_carbon_bank.streams[plant].write(file_name, hourly_plant_carbon_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                    .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 1;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilWaterEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const values = try hourly_soil_water_bank.row(cell);
                            const liquid_fraction = values[6 .. 6 + config.soil_layers];
                            const surface_liquid_index = 6 + config.soil_layers;
                            const ice_fraction = values[surface_liquid_index + 1 ..][0..config.soil_layers];
                            var evapotranspiration_m3: f64 = 0;
                            var root_water_uptake_m3: f64 = 0;
                            for (0..config.plant_populations) |species| {
                                const plant = cell * config.plant_populations + species;
                                if (canopy_surface_exchange_state) |exchange| {
                                    evapotranspiration_m3 += @max(0, -exchange.transpiration_m3_per_h[plant]);
                                    evapotranspiration_m3 += @max(0, -exchange.intercepted_water_change_m3_per_h[plant]);
                                }
                                evapotranspiration_m3 += @max(0, -standing_dead_evaporation_m3_per_h[plant]);
                                if (plant_root_state) |roots| for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..config.soil_layers) |local_layer| {
                                    root_water_uptake_m3 += @max(0, -roots.water_uptake_m3_per_h[try roots.layerIndex(plant, domain, local_layer)]);
                                };
                            }
                            var external_water_outflow_m3: f64 = 0;
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const layer_volume_m3 = soil_solver_property_state.layer_volume_m3[layer];
                                if (!std.math.isFinite(layer_volume_m3) or layer_volume_m3 <= 0) return error.InvalidOutputSoilLayerVolume;
                                liquid_fraction[local_layer] = state.liquid_water_m3[layer] / layer_volume_m3;
                                ice_fraction[local_layer] = state.ice_water_m3[layer] / layer_volume_m3;
                                external_water_outflow_m3 += @max(0, transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer]);
                                external_water_outflow_m3 += @max(0, transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]);
                            }
                            const litter_capacity_m3 = surface_precipitation_state.litter_water_capacity_m3[cell];
                            const surface_liquid_fraction = if (litter_capacity_m3 > 0) surface_precipitation_state.litter_water_m3[cell] / litter_capacity_m3 else 0;
                            const surface_ice_fraction = if (litter_capacity_m3 > 0) surface_litter_ice_m3[cell] / litter_capacity_m3 else 0;
                            const top_layer = try state.layerIndex(cell, 0);
                            const bulk_density_megagrams_per_m3 = soil_solver_property_state.bulk_density_megagrams_per_m3[top_layer];
                            if (!std.math.isFinite(bulk_density_megagrams_per_m3) or bulk_density_megagrams_per_m3 <= 0) return error.InvalidOutputSoilBulkDensity;
                            try ecosys.soil_water_output.calculateInto(.{
                                .evapotranspiration_m3 = evapotranspiration_m3,
                                .runoff_m3 = -surface_runoff_state.exported_water_m3[cell],
                                .sediment_discharge_water_m3 = surface_erosion_state.routing.sediment_export_megagrams[cell] / bulk_density_megagrams_per_m3,
                                .root_water_uptake_m3 = root_water_uptake_m3,
                                .external_water_outflow_m3 = external_water_outflow_m3,
                                .surface_snow_volume_m3 = surface_precipitation_state.solid_snow_water_equivalent_m3[cell],
                                .surface_ice_volume_m3 = surface_litter_ice_m3[cell],
                                .surface_liquid_water_m3 = surface_precipitation_state.litter_water_m3[cell],
                                // Dall'Amico soil/surface ice is already
                                // stored as water-equivalent volume.
                                .ice_density_megagrams_per_m3 = 1,
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .total_grid_area_m2 = total_grid_area_m2,
                                .volumetric_liquid_water_fraction_by_layer = liquid_fraction,
                                .surface_volumetric_liquid_water_fraction = surface_liquid_fraction,
                                .volumetric_ice_fraction_by_layer = ice_fraction,
                                .surface_volumetric_ice_fraction = surface_ice_fraction,
                                .active_layer_boundary_depth_m = -soil_boundary_topology_state.active_layer_depth_m[cell],
                                .water_table_boundary_depth_m = -soil_boundary_topology_state.internal_water_table_depth_m[cell],
                                .soil_surface_reference_depth_m = 0,
                            }, values);
                            const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_water_bank.streams[cell].write(file_name, hourly_soil_water_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 1;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantWaterEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantWaterOutputRequiresCanopyState;
                            const roots = if (plant_root_state) |*value| value else return error.HourlyPlantWaterOutputRequiresRootState;
                            const surface = if (canopy_surface_exchange_state) |*value| value else return error.HourlyPlantWaterOutputRequiresSurfaceExchangeState;
                            const surface_workspace = if (canopy_surface_input_workspace) |*value| value else return error.HourlyPlantWaterOutputRequiresSurfaceWorkspace;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const values = try hourly_plant_water_bank.row(plant);
                                const primary_root_first = try roots.layerIndex(plant, 0, 0);
                                const primary_root_potential = roots.total_water_potential_megapascal[primary_root_first .. primary_root_first + config.soil_layers];
                                const outward_water_flux_m3_per_h = -(surface.transpiration_m3_per_h[plant] +
                                    surface.intercepted_water_change_m3_per_h[plant] +
                                    standing_dead_evaporation_m3_per_h[plant]);
                                try ecosys.plant_hourly_output.calculateWaterInto(
                                    canopy.plant_canopy_osmotic_potential_megapascal[plant] + canopy.plant_canopy_turgor_potential_megapascal[plant],
                                    canopy.plant_canopy_turgor_potential_megapascal[plant],
                                    surface_workspace.stomatal_resistance_h_per_m[plant],
                                    surface.boundary_layer_resistance_h_per_m[plant],
                                    outward_water_flux_m3_per_h,
                                    phenology_root_oxygen_fraction[plant],
                                    primary_root_potential,
                                    canopy_cell_area_m2[cell],
                                    values,
                                );
                                // Only assigned species produce output; an unassigned population
                                // would emit all-zero rows, so no file is created for it.
                                const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_water_bank.streams[plant].write(file_name, hourly_plant_water_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                    .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 2;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilNitrogenEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const root_gas_withdrawal = if (plant_root_state) |*roots|
                                try ecosys.plant_root_disturbance.rootGasWithdrawalForCell(roots, cell, config.plant_populations)
                            else
                                ecosys.plant_root_disturbance.CellRootGasWithdrawal{};
                            const values = try hourly_soil_nitrogen_bank.row(cell);
                            const nitrous_oxide_concentration = values[5 .. 5 + config.soil_layers];
                            const litter_nitrous_oxide_index = 5 + config.soil_layers;
                            const ammonia_concentration = values[litter_nitrous_oxide_index + 1 ..][0..config.soil_layers];
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const water_m3 = state.liquid_water_m3[layer];
                                const n2o_component = layer * ecosys.gas_transport.species_count + @intFromEnum(ecosys.gas_transport.Species.nitrous_oxide);
                                nitrous_oxide_concentration[local_layer] = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(gas_transport_state.dissolved_mass_g[n2o_component], water_m3);
                            }
                            const first_soil_layer = try state.layerIndex(cell, 0);
                            try ecosys.soil_hourly_output_binding.writeMineralAmmoniaNitrogenProfile(
                                initial_chemistry_state.aqueous,
                                first_soil_layer,
                                config.soil_layers,
                                runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                                ammonia_concentration,
                            );
                            const nitrous_oxide_exchange_g_n_per_h = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                soil_gas_transport_state.atmospheric_flux_g_per_h,
                                surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                first_soil_layer,
                                config.soil_layers,
                                cell,
                                .nitrous_oxide,
                            );
                            const dinitrogen_exchange_g_n_per_h = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                soil_gas_transport_state.atmospheric_flux_g_per_h,
                                surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                first_soil_layer,
                                config.soil_layers,
                                cell,
                                .nitrogen,
                            );
                            const ammonia_exchange_g_n_per_h = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                soil_gas_transport_state.atmospheric_flux_g_per_h,
                                surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                first_soil_layer,
                                config.soil_layers,
                                cell,
                                .ammonia,
                            );
                            const litter_base = cell * ecosys.gas_transport.species_count;
                            const litter_nitrous_oxide = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(
                                litter_gas_transport_state.dissolved_mass_g[litter_base + @intFromEnum(ecosys.gas_transport.Species.nitrous_oxide)],
                                surface_litter_gas_transport_state.water_volume_m3[cell],
                            );
                            const litter_ammonia_mol_n_per_m3 = surface_litter_chemistry_state.cells[cell].ammonia_mol_per_m3;
                            if (!std.math.isFinite(litter_ammonia_mol_n_per_m3) or litter_ammonia_mol_n_per_m3 < 0) return error.InvalidOutputLitterAmmoniaConcentration;
                            var dissolved_inorganic_nitrogen_drainage_g_n_per_h: f64 = 0;
                            for (0..config.soil_layers) |local_layer| {
                                dissolved_inorganic_nitrogen_drainage_g_n_per_h += mineral_nitrogen_transport_state.boundary_export_g_n_per_step[try state.layerIndex(cell, local_layer)];
                            }
                            const root_atmospheric_nitrous_oxide_g_n_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 3);
                            const root_atmospheric_ammonia_g_n_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 4);
                            try ecosys.soil_biogeochemistry_output.calculateNitrogenInto(.{
                                .nitrous_oxide_emission_g_n_per_h = nitrous_oxide_exchange_g_n_per_h + root_gas_withdrawal.nitrous_oxide_g_n_per_h - root_atmospheric_nitrous_oxide_g_n_per_h,
                                .dinitrogen_emission_g_n_per_h = dinitrogen_exchange_g_n_per_h,
                                .ammonia_emission_g_n_per_h = ammonia_exchange_g_n_per_h + root_gas_withdrawal.ammonia_g_n_per_h - root_atmospheric_ammonia_g_n_per_h,
                                .dissolved_inorganic_nitrogen_runoff_g_n_per_h = surface_inorganic_nitrogen_export_g_n_per_h[cell],
                                .dissolved_inorganic_nitrogen_drainage_g_n_per_h = dissolved_inorganic_nitrogen_drainage_g_n_per_h,
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .total_grid_area_m2 = total_grid_area_m2,
                                .nitrous_oxide_concentration_by_layer = nitrous_oxide_concentration,
                                .litter_nitrous_oxide_concentration = litter_nitrous_oxide,
                                .ammonia_concentration_by_layer = ammonia_concentration,
                                .litter_ammonia_concentration = litter_ammonia_mol_n_per_m3 * 14,
                            }, values);
                            const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_nitrogen_bank.streams[cell].write(file_name, hourly_soil_nitrogen_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 2;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantNitrogenEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantNitrogenOutputRequiresCanopyState;
                            const roots = if (plant_root_state) |*value| value else return error.HourlyPlantNitrogenOutputRequiresRootState;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const root_first = try roots.layerIndex(plant, 0, 0);
                                const mycorrhiza_first = try roots.layerIndex(plant, 1, 0);
                                const has_mycorrhiza = root_biological_domain_count_by_plant[plant] == 2;
                                const first_soil_layer = try state.layerIndex(cell, 0);
                                const values = try hourly_plant_nitrogen_bank.row(plant);
                                try ecosys.plant_hourly_output.calculateNitrogenInto(
                                    roots.ammonium_uptake_g_n_per_h[plant],
                                    roots.nitrate_uptake_g_n_per_h[plant],
                                    roots.fixation_uptake_g_n_per_h[plant],
                                    canopy.plant_mobile_nitrogen_concentration_g_per_g[plant],
                                    canopy_ammonia_publication_state
                                        .current_exchange_g_n_per_h_by_plant[plant],
                                    canopy_cell_area_m2[cell],
                                    .{
                                        .ammonium_root_non_band_g_n_per_h = roots.ammonium_uptake_nonband_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .ammonium_mycorrhiza_non_band_g_n_per_h = if (has_mycorrhiza) roots.ammonium_uptake_nonband_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .ammonium_root_band_g_n_per_h = roots.ammonium_uptake_band_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .ammonium_mycorrhiza_band_g_n_per_h = if (has_mycorrhiza) roots.ammonium_uptake_band_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .nitrate_root_non_band_g_n_per_h = roots.nitrate_uptake_nonband_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .nitrate_mycorrhiza_non_band_g_n_per_h = if (has_mycorrhiza) roots.nitrate_uptake_nonband_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .nitrate_root_band_g_n_per_h = roots.nitrate_uptake_band_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .nitrate_mycorrhiza_band_g_n_per_h = if (has_mycorrhiza) roots.nitrate_uptake_band_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .layer_area_m2 = soil_hourly_workspace.horizontal_face_area_m2[first_soil_layer .. first_soil_layer + config.soil_layers],
                                    },
                                    values,
                                );
                                // Only assigned species produce output; an unassigned population
                                // would emit all-zero rows, so no file is created for it.
                                const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_nitrogen_bank.streams[plant].write(file_name, hourly_plant_nitrogen_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                    .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 4;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilHeatEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const values = try hourly_soil_heat_bank.row(cell);
                            const soil_temperature_c = values[13 .. 13 + config.soil_layers];
                            for (soil_temperature_c, 0..) |*temperature_c, local_layer| temperature_c.* = state.soil_temperature_k[try state.layerIndex(cell, local_layer)] - 273.15;
                            const litter_air_volume_m3 = litter_gas_transport_state.air_volume_m3[cell];
                            const litter_vapor_pressure_kpa = if (litter_air_volume_m3 > 0)
                                litter_gas_transport_state.water_vapor_mol[cell] * 8.31446261815324e-3 * litter_gas_transport_state.temperature_k[cell] / litter_air_volume_m3
                            else
                                0;
                            try ecosys.soil_heat_output.calculateInto(.{
                                .incoming_shortwave_radiation_megajoules_per_m2_h = atmospheric_state.shortwave_radiation_megajoules_per_m2[cell],
                                .air_temperature_c = atmospheric_state.air_temperature_k[cell] - 273.15,
                                .atmospheric_vapor_pressure_kpa = atmospheric_state.vapor_pressure_kpa[cell],
                                .wind_travel_m_per_h = atmospheric_state.wind_speed_m_per_h[cell],
                                .rainfall_m3 = atmospheric_state.rainfall_m[cell] * canopy_cell_area_m2[cell],
                                .irrigation_m3 = irrigation_water_depth_m[cell] * canopy_cell_area_m2[cell],
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .ground_surface_net_radiation_megajoules = ecosystem_energy_ledger_state.ground_surface_net_radiation_megajoules[cell],
                                .ground_surface_latent_heat_megajoules = ecosystem_energy_ledger_state.ground_surface_latent_heat_megajoules[cell],
                                .ground_surface_sensible_heat_megajoules = ecosystem_energy_ledger_state.ground_surface_sensible_heat_megajoules[cell],
                                .ground_surface_storage_heat_megajoules = ecosystem_energy_ledger_state.ground_surface_storage_heat_megajoules[cell],
                                .ecosystem_net_radiation_megajoules = ecosystem_energy_ledger_state.ecosystem_net_radiation_megajoules[cell],
                                .ecosystem_latent_heat_megajoules = ecosystem_energy_ledger_state.ecosystem_latent_heat_megajoules[cell],
                                .ecosystem_sensible_heat_megajoules = ecosystem_energy_ledger_state.ecosystem_sensible_heat_megajoules[cell],
                                .ecosystem_storage_heat_megajoules = ecosystem_energy_ledger_state.ecosystem_storage_heat_megajoules[cell],
                                .soil_temperature_c_by_layer = soil_temperature_c,
                                .surface_soil_temperature_c = state.soil_temperature_k[try state.layerIndex(cell, 0)] - 273.15,
                                .surface_water_temperature_c = state.surface_temperature_k[cell] - 273.15,
                                .litter_temperature_c = litter_gas_transport_state.temperature_k[cell] - 273.15,
                                .litter_water_vapor_partial_pressure_kpa = litter_vapor_pressure_kpa,
                                .litter_absolute_temperature_k = litter_gas_transport_state.temperature_k[cell],
                            }, values);
                            const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_heat_bank.streams[cell].write(file_name, hourly_soil_heat_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 4;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantHeatEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantHeatOutputRequiresCanopyState;
                            const living = if (canopy_surface_exchange_state) |*value| value else return error.HourlyPlantHeatOutputRequiresSurfaceExchangeState;
                            const dead = if (standing_dead_surface_exchange_state) |*value| value else return error.HourlyPlantHeatOutputRequiresStandingDeadExchangeState;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const temperature_function = canopy.plant_uptake_growth_temperature_response[plant];
                                const output = try ecosys.plant_hourly_output.heat(
                                    canopy_net_radiation_megajoules[plant] + dead.net_radiation_megajoules_per_h[plant],
                                    living.latent_heat_flux_megajoules_per_h[plant] + dead.latent_heat_flux_megajoules_per_h[plant],
                                    living.sensible_heat_flux_megajoules_per_h[plant] + dead.sensible_heat_flux_megajoules_per_h[plant],
                                    -canopy_storage_heat_megajoules[plant] + living.vapor_sensible_heat_flux_megajoules_per_h[plant] - dead.storage_heat_flux_megajoules_per_h[plant] + dead.vapor_sensible_heat_flux_megajoules_per_h[plant],
                                    plant_state.canopy_temperature_k[plant] - 273.15,
                                    temperature_function,
                                    canopy.plant_standing_dead_surface_temperature_k[plant] - 273.15,
                                    canopy_cell_area_m2[cell],
                                );
                                const values = try hourly_plant_heat_bank.row(plant);
                                values[0..7].* = output.values();
                                // Only assigned species produce output; an unassigned population
                                // would emit all-zero rows, so no file is created for it.
                                const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_heat_bank.streams[plant].write(file_name, hourly_plant_heat_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                    .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 3;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilPhosphorusEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const values = try hourly_soil_phosphorus_bank.row(cell);
                            var dissolved_inorganic_phosphorus_drainage_g_p_per_h: f64 = 0;
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const base = layer * aqueous_species_count;
                                inline for (@typeInfo(ecosys.solute_transport_species.AqueousSpecies).@"enum".fields) |field| {
                                    const species: ecosys.solute_transport_species.AqueousSpecies = @enumFromInt(field.value);
                                    if (ecosys.solute_transport_species.diffusivityClass(species) == .phosphate) {
                                        dissolved_inorganic_phosphorus_drainage_g_p_per_h += @max(0, -soil_solute_boundary_net_flux_mol[base + field.value]) * runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol;
                                    }
                                }
                            }
                            const phosphorus = try ecosys.soil_biogeochemistry_output.calculatePhosphorus(
                                surface_inorganic_phosphorus_export_g_p_per_h[cell],
                                dissolved_inorganic_phosphorus_drainage_g_p_per_h,
                                total_grid_area_m2,
                            );
                            values[0..2].* = ecosys.soil_biogeochemistry_output.phosphorusValues(phosphorus);
                            const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_phosphorus_bank.streams[cell].write(file_name, hourly_soil_phosphorus_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 3;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantPhosphorusEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantPhosphorusOutputRequiresCanopyState;
                            const roots = if (plant_root_state) |*value| value else return error.HourlyPlantPhosphorusOutputRequiresRootState;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const root_first = try roots.layerIndex(plant, 0, 0);
                                const mycorrhiza_first = try roots.layerIndex(plant, 1, 0);
                                const has_mycorrhiza = root_biological_domain_count_by_plant[plant] == 2;
                                const first_soil_layer = try state.layerIndex(cell, 0);
                                const root_non_band = roots.phosphate_h2_uptake_nonband_g_p_per_h[root_first .. root_first + config.soil_layers];
                                const mycorrhiza_non_band = if (has_mycorrhiza) roots.phosphate_h2_uptake_nonband_g_p_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h;
                                const root_band = roots.phosphate_h2_uptake_band_g_p_per_h[root_first .. root_first + config.soil_layers];
                                const mycorrhiza_band = if (has_mycorrhiza) roots.phosphate_h2_uptake_band_g_p_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h;
                                var h2_phosphate_uptake_g_p_per_h: f64 = 0;
                                for (0..config.soil_layers) |layer| h2_phosphate_uptake_g_p_per_h += root_non_band[layer] + mycorrhiza_non_band[layer] + root_band[layer] + mycorrhiza_band[layer];
                                const values = try hourly_plant_phosphorus_bank.row(plant);
                                try ecosys.plant_hourly_output.calculatePhosphorusInto(
                                    h2_phosphate_uptake_g_p_per_h,
                                    canopy.plant_mobile_phosphorus_concentration_g_per_g[plant],
                                    root_non_band,
                                    mycorrhiza_non_band,
                                    root_band,
                                    mycorrhiza_band,
                                    soil_hourly_workspace.horizontal_face_area_m2[first_soil_layer .. first_soil_layer + config.soil_layers],
                                    canopy_cell_area_m2[cell],
                                    values,
                                );
                                // Only assigned species produce output; an unassigned population
                                // would emit all-zero rows, so no file is created for it.
                                const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_phosphorus_bank.streams[plant].write(file_name, hourly_plant_phosphorus_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                    .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (active_options.visualizationIncludesYear(current_year)) {
                    for (0..state.cell_count) |cell| {
                        const first_plant = cell * config.plant_populations;
                        const plant_end =
                            first_plant + config.plant_populations;
                        var ecosystem_gross_primary_productivity_g: f64 = 0;
                        var ecosystem_signed_autotrophic_respiration_g: f64 = 0;
                        var harvested_carbon_g: f64 = 0;
                        for (0..config.plant_populations) |species| {
                            const plant = first_plant + species;
                            ecosystem_gross_primary_productivity_g +=
                                plant_daily_flux_ledger.gross_primary_productivity_g[plant];
                            ecosystem_signed_autotrophic_respiration_g +=
                                plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
                            harvested_carbon_g +=
                                plant_daily_flux_ledger.harvested_carbon_g[plant];
                            visualization_species_respired_carbon_g[species] =
                                -plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
                        }
                        const ecosystem_net_primary_productivity_g =
                            ecosystem_gross_primary_productivity_g +
                            ecosystem_signed_autotrophic_respiration_g;
                        const ecosystem_carbon = try ecosys.soil_daily_ecosystem_carbon.calculate(.{
                            .net_primary_productivity_g_c = ecosystem_net_primary_productivity_g,
                            .signed_heterotrophic_respiration_g_c = daily_heterotrophic_respiration.daily_signed_g_c[cell],
                            .dissolved_organic_carbon_runoff_g_c = daily_carbon_export.dissolved_organic_carbon_runoff_g[cell],
                            .dissolved_inorganic_carbon_runoff_g_c = daily_carbon_export.dissolved_inorganic_carbon_runoff_g[cell],
                            .dissolved_organic_carbon_drainage_g_c = daily_carbon_export.dissolved_organic_carbon_drainage_g[cell],
                            .dissolved_inorganic_carbon_drainage_g_c = daily_carbon_export.dissolved_inorganic_carbon_drainage_g[cell],
                            .harvested_carbon_g_c = harvested_carbon_g,
                            .organic_fertilizer_carbon_input_g_c = daily_biome_organic_carbon_input_g_c[cell],
                        });
                        const first_layer = try state.layerIndex(cell, 0);
                        const layer_end = first_layer + config.soil_layers;
                        _ = try visualization_streams.calculateAndWriteHourly(
                            cell,
                            site_by_cell[cell].longitude_degrees_east,
                            site_by_cell[cell].latitude_degrees_north,
                            .{
                                .year = current_year,
                                .day_of_year = current_day_of_year,
                                .month = management_date.month,
                                .day = management_date.day,
                                .hour = if (timestamp.hour == 24)
                                    23
                                else
                                    timestamp.hour,
                            },
                            .{
                                .cell_area_m2 = canopy_cell_area_m2[cell],
                                .species_cumulative_gross_carbon_g = plant_daily_flux_ledger.gross_primary_productivity_g[first_plant..plant_end],
                                .species_cumulative_respired_carbon_g = visualization_species_respired_carbon_g,
                                .ecosystem_gross_primary_productivity_g = ecosystem_gross_primary_productivity_g,
                                .ecosystem_net_primary_productivity_g = ecosystem_net_primary_productivity_g,
                                .ecosystem_autotrophic_respiration_g = -ecosystem_signed_autotrophic_respiration_g,
                                .ecosystem_heterotrophic_respiration_g = -daily_heterotrophic_respiration.daily_signed_g_c[cell],
                                .net_biome_productivity_g = ecosystem_carbon.net_biome_productivity_g_c,
                                .carbon_dioxide_flux_g = (try daily_soil_gas_flux.get(cell, .carbon_dioxide)) +
                                    daily_canopy_gas_exchange.net_carbon_dioxide_uptake_g_c[cell],
                                .methane_flux_g = (try daily_soil_gas_flux.get(cell, .methane)) +
                                    daily_canopy_gas_exchange.net_methane_uptake_g_c[cell],
                                .net_radiation_megajoules_h = ecosystem_energy_ledger_state.ecosystem_net_radiation_megajoules[cell],
                                .latent_heat_flux_megajoules_h = ecosystem_energy_ledger_state.ecosystem_latent_heat_megajoules[cell],
                                .sensible_heat_flux_megajoules_h = ecosystem_energy_ledger_state.ecosystem_sensible_heat_megajoules[cell],
                                .liquid_water_m3_by_layer = state.matrix_liquid_water_m3[first_layer..layer_end],
                                .air_volume_m3_by_layer = state.matrix_air_volume_m3[first_layer..layer_end],
                                .macropore_water_m3_by_layer = state.macropore_liquid_water_m3[first_layer..layer_end],
                                .total_volume_m3_by_layer = soil_solver_property_state.layer_volume_m3[first_layer..layer_end],
                                .soil_temperature_k_by_layer = state.soil_temperature_k[first_layer..layer_end],
                            },
                        );
                    }
                }
                if ((scene_weather_hours + 1) % 24 == 0) {
                    // The midnight step belongs to the completed day. Use the
                    // previous step's timestamp to label this output correctly
                    // regardless of whether the weather file uses hour 0000 or
                    // hour 2400 for the day transition.
                    const completed_day_of_year = try dayOfYearFromTimestamp(previous_weather_timestamp.?);
                    const completed_management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(previous_weather_timestamp.?);
                    try daily_ecosystem_carbon.finalizeDay(
                        &state,
                        config.plant_populations,
                        plant_daily_flux_ledger.gross_primary_productivity_g,
                        plant_daily_flux_ledger.signed_total_respiration_carbon_g,
                        plant_daily_flux_ledger.harvested_carbon_g,
                        daily_heterotrophic_respiration.daily_signed_g_c,
                        daily_biome_organic_carbon_input_g_c,
                        &daily_carbon_export,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedWater(
                        daily_water_ledger.rainfall_m3,
                        daily_water_ledger.boundary_water_inflow_m3,
                        daily_water_ledger.runoff_m3,
                        daily_water_ledger.evaporation_m3,
                        daily_water_ledger.water_outflow_m3,
                        daily_water_ledger.lateral_water_outflow_m3,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedFertilizer(
                        soil_fertilizer_inventory.daily_nitrogen_input_g_n,
                        daily_organic_fertilizer_nitrogen_input_g_n,
                        mineral_fertilizer_inventory.daily_phosphorus_input_g_p,
                        daily_organic_fertilizer_phosphorus_input_g_p,
                        daily_organic_fertilizer_carbon_input_g_c,
                        daily_biome_organic_carbon_input_g_c,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedExports(
                        .{
                            daily_carbon_export.dissolved_organic_carbon_runoff_g,
                            daily_carbon_export.dissolved_inorganic_carbon_runoff_g,
                            daily_carbon_export.dissolved_organic_carbon_drainage_g,
                            daily_carbon_export.dissolved_inorganic_carbon_drainage_g,
                        },
                        .{
                            daily_nitrogen_export.dissolved_organic_nitrogen_runoff_g_n,
                            daily_nitrogen_export.dissolved_inorganic_nitrogen_runoff_g_n,
                            daily_nitrogen_export.dissolved_organic_nitrogen_drainage_g_n,
                            daily_nitrogen_export.dissolved_inorganic_nitrogen_drainage_g_n,
                        },
                        .{
                            daily_phosphorus_export.dissolved_organic_phosphorus_runoff_g_p,
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_runoff_g_p,
                            daily_phosphorus_export.dissolved_organic_phosphorus_drainage_g_p,
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_drainage_g_p,
                        },
                        daily_heat_ledger.ionic_outflow_mol,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedCarbonTransportInputs(
                        daily_carbon_export.dissolved_organic_carbon_input_g,
                        daily_carbon_export.dissolved_inorganic_carbon_input_g,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedNitrogenTransportInputs(
                        daily_nitrogen_export.dissolved_organic_nitrogen_input_g_n,
                        daily_nitrogen_export.dissolved_inorganic_nitrogen_input_g_n,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedPhosphorusTransportInputs(
                        daily_phosphorus_export.dissolved_organic_phosphorus_input_g_p,
                        daily_phosphorus_export.dissolved_inorganic_phosphorus_input_g_p,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedPlantLitter(
                        plant_daily_flux_ledger.carbon_sink_g,
                        plant_daily_flux_ledger.nitrogen_sink_g,
                        plant_daily_flux_ledger.phosphorus_sink_g,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedRootSoilOrganicExchange(
                        plant_daily_flux_ledger.root_soil_carbon_exchange_g,
                        plant_daily_flux_ledger.root_soil_nitrogen_exchange_g,
                        plant_daily_flux_ledger.root_soil_phosphorus_exchange_g,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedAtmosphericGas(
                        &daily_soil_gas_flux,
                        &daily_canopy_gas_exchange,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedRedistSurfaceGas(
                        try daily_soil_gas_flux.getRedistSurfaceGasTotals(),
                    );
                    const atmospheric_parameters = runscript.chemistry_primary_initialization.molar_mass_g_per_mol;
                    const atmospheric_ion_molar_masses: ecosys.snow_surface_discharge.IonMolarMassesGPerMol = .{
                        .aluminum = atmospheric_parameters.aluminum,
                        .iron = atmospheric_parameters.iron,
                        .calcium = atmospheric_parameters.calcium,
                        .magnesium = atmospheric_parameters.magnesium,
                        .sodium = atmospheric_parameters.sodium,
                        .potassium = atmospheric_parameters.potassium,
                        .sulfur = atmospheric_parameters.sulfur,
                        .chloride = atmospheric_parameters.chloride,
                    };
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedAtmosphericSolutes(
                        &atmospheric_solute_input_ledger_state,
                        atmospheric_ion_molar_masses,
                    );
                    const completed_scene_day = (scene_weather_hours + 1) / 24;
                    if (try ecosys.mass_balance_audit.shouldAudit(
                        completed_scene_day,
                        active_options.daily_output_interval_days,
                    )) {
                        const totals = try reconstructLandscapeMassBalance(
                            hourly_science_context,
                        );
                        const monitor = landscape_mass_balance_state.monitor orelse
                            return error.MassBalanceMonitorNotInitialized;
                        {
                            const area = totals.landscape_area_m2;
                            const b = try ecosys.mass_balance_audit.balance(totals);
                            std.log.debug("phosphorus day-end: residue={e} organic={e} phosphate={e} input={e} output={e} fertilizer={e} sink={e} total_g_m2={e}", .{ totals.residue_phosphorus_g / area, totals.organic_phosphorus_g / area, totals.phosphate_phosphorus_g / area, totals.cumulative_phosphorus_input_g / area, totals.cumulative_phosphorus_output_g / area, totals.cumulative_organic_fertilizer_phosphorus_g / area, totals.cumulative_phosphorus_sink_g / area, b.phosphorus_g / area });
                            std.log.debug("nitrogen day-end: residue={e} organic={e} n2={e} nh4={e} no3={e} n2in={e} nin={e} nout={e} fert={e} sink={e} root_uptake={e} root_exudate={e} total_g_m2={e}", .{ totals.residue_nitrogen_g / area, totals.organic_nitrogen_g / area, totals.dinitrogen_nitrogen_g / area, totals.ammonium_nitrogen_g / area, totals.nitrate_nitrogen_g / area, totals.cumulative_dinitrogen_input_g / area, totals.cumulative_nitrogen_input_g / area, totals.cumulative_nitrogen_output_g / area, totals.cumulative_organic_fertilizer_nitrogen_g / area, totals.cumulative_nitrogen_sink_g / area, totals.cumulative_plant_root_organic_nitrogen_uptake_g / area, totals.cumulative_plant_root_organic_nitrogen_exudate_g / area, b.nitrogen_g / area });
                            var diagnostic_gas_n_g: f64 = 0;
                            var diagnostic_boundary_gas_n_g: f64 = 0;
                            var diagnostic_don_runoff_g_n: f64 = 0;
                            var diagnostic_din_runoff_g_n: f64 = 0;
                            var diagnostic_don_drainage_g_n: f64 = 0;
                            var diagnostic_din_drainage_g_n: f64 = 0;
                            for (0..state.cell_count) |cell| {
                                inline for (.{ ecosys.gas_transport.Species.nitrogen, .nitrous_oxide, .ammonia }) |species| {
                                    diagnostic_gas_n_g += try daily_soil_gas_flux.get(cell, species);
                                    diagnostic_boundary_gas_n_g += try daily_soil_gas_flux.getSoilLitterBoundary(cell, species);
                                }
                                diagnostic_don_runoff_g_n += daily_nitrogen_export.dissolved_organic_nitrogen_runoff_g_n[cell];
                                diagnostic_din_runoff_g_n += daily_nitrogen_export.dissolved_inorganic_nitrogen_runoff_g_n[cell];
                                diagnostic_don_drainage_g_n += daily_nitrogen_export.dissolved_organic_nitrogen_drainage_g_n[cell];
                                diagnostic_din_drainage_g_n += daily_nitrogen_export.dissolved_inorganic_nitrogen_drainage_g_n[cell];
                            }
                            std.log.debug("nitrogen boundary probe: gas_combined={e} gas_soil_litter={e} don_runoff={e} din_runoff={e} don_drainage={e} din_drainage={e}", .{ diagnostic_gas_n_g / area, diagnostic_boundary_gas_n_g / area, diagnostic_don_runoff_g_n / area, diagnostic_din_runoff_g_n / area, diagnostic_don_drainage_g_n / area, diagnostic_din_drainage_g_n / area });
                        }
                        _ = try monitor.check(
                            completed_scene_day,
                            current_year,
                            totals,
                        );
                    }
                    if (@mod(completed_scene_day, active_options.daily_output_interval_days) == 0) {
                        const carbon_editor_index: usize = 5;
                        const carbon_editor_name = scene.output_editors[carbon_editor_index];
                        if (!ecosys.delimited_input.isNo(carbon_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + carbon_editor_index] orelse return error.UnresolvedDailyPlantCarbonEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(carbon_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                for (0..state.cell_count) |cell| {
                                    const values = try daily_soil_carbon_bank.row(cell);
                                    @memset(values, std.math.nan(f64));
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const layer_carbon = daily_soil_organic_carbon_by_layer_g_c[first_layer .. first_layer + config.soil_layers];
                                    @memset(layer_carbon, 0);
                                    const pools = try ecosys.soil_daily_carbon_pools.calculate(
                                        &soil_organic_state,
                                        &surface_organic_state,
                                        first_layer,
                                        active_layers,
                                        cell,
                                        layer_carbon[0..active_layers],
                                    );
                                    const plant_flux = try ecosys.soil_daily_plant_carbon_flux.calculate(&plant_daily_flux_ledger, cell, config.plant_populations);
                                    var gas_inventory_g: [ecosys.gas_transport.species_count]f64 = @splat(0);
                                    for (0..active_layers) |local_layer| {
                                        const layer = first_layer + local_layer;
                                        for (0..ecosys.gas_transport.species_count) |species| {
                                            const component = layer * ecosys.gas_transport.species_count + species;
                                            gas_inventory_g[species] +=
                                                gas_transport_state.gaseous_mass_g[component] +
                                                gas_transport_state.dissolved_mass_g[component] +
                                                gas_transport_state.macropore_dissolved_mass_g[component] +
                                                gas_transport_state.band_dissolved_mass_g[component];
                                        }
                                    }
                                    for (0..ecosys.gas_transport.species_count) |species| {
                                        const component = cell * ecosys.gas_transport.species_count + species;
                                        gas_inventory_g[species] +=
                                            litter_gas_transport_state.gaseous_mass_g[component] +
                                            litter_gas_transport_state.dissolved_mass_g[component] +
                                            litter_gas_transport_state.macropore_dissolved_mass_g[component] +
                                            litter_gas_transport_state.band_dissolved_mass_g[component];
                                    }
                                    const area_inverse = 1.0 / canopy_cell_area_m2[cell];
                                    values[0] = pools.residue_carbon_g_c * area_inverse;
                                    values[1] = pools.humus_carbon_g_c * area_inverse;
                                    values[2] = daily_organic_fertilizer_carbon_input_g_c[cell] * area_inverse;
                                    values[3] = plant_flux.plant_litterfall_g_c * area_inverse;
                                    values[4] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .carbon_dioxide)) * area_inverse;
                                    values[5] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .oxygen)) * area_inverse;
                                    values[6] = daily_heterotrophic_respiration.daily_carbon_dioxide_production_g_c[cell] * area_inverse;
                                    values[7] = pools.microbial_carbon_g_c * area_inverse;
                                    values[8] = pools.surface_noncharcoal_organic_carbon_g_c * area_inverse;
                                    values[9] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .methane)) * area_inverse;
                                    values[10] = daily_carbon_export.dissolved_organic_carbon_runoff_g[cell] / total_grid_area_m2;
                                    values[11] = daily_carbon_export.dissolved_organic_carbon_drainage_g[cell] / total_grid_area_m2;
                                    values[12] = daily_carbon_export.dissolved_inorganic_carbon_runoff_g[cell] / total_grid_area_m2;
                                    values[13] = daily_carbon_export.dissolved_inorganic_carbon_drainage_g[cell] / total_grid_area_m2;
                                    values[14] = current_atmospheric_co2_umol_per_mol;
                                    values[15] = daily_ecosystem_carbon.net_biome_productivity_g_c[cell] * area_inverse;
                                    values[16] = (daily_soil_fire_carbon_dioxide_emission_g_c[cell] + daily_root_fire_carbon_dioxide_emission_g_c[cell]) * area_inverse;
                                    for (0..config.soil_layers) |layer| values[17 + layer] = layer_carbon[layer] * area_inverse;
                                    const after_layers = 17 + config.soil_layers;
                                    values[after_layers] = daily_soil_fire_charcoal_production_g_c[cell] * area_inverse;
                                    values[after_layers + 1] =
                                        ((try daily_soil_gas_flux.get(cell, .carbon_dioxide)) +
                                            daily_canopy_gas_exchange.net_carbon_dioxide_uptake_g_c[cell]) * area_inverse;
                                    values[after_layers + 2] =
                                        ((try daily_soil_gas_flux.get(cell, .methane)) +
                                            daily_canopy_gas_exchange.net_methane_uptake_g_c[cell]) * area_inverse;
                                    values[after_layers + 3] =
                                        ((try daily_soil_gas_flux.get(cell, .oxygen)) +
                                            daily_canopy_gas_exchange.net_oxygen_uptake_g_o[cell]) * area_inverse;
                                    // OUTSD choices 36..40 are source-reserved zeros.
                                    @memset(values[after_layers + 4 .. after_layers + 9], 0);
                                    const tail = after_layers + 9;
                                    // OUTSD choice 41 is UH2GG: accepted daily
                                    // ecosystem-atmosphere H2 exchange, not
                                    // the end-of-day H2 storage inventory.
                                    values[tail] = try ecosys.soil_daily_output.dailyHydrogenFluxGPerM2(
                                        &daily_soil_gas_flux,
                                        cell,
                                        canopy_cell_area_m2[cell],
                                    );
                                    values[tail + 1] = plant_flux.harvested_carbon_g_c * area_inverse;
                                    // With no active plant branches, OUTSD ARLFC and WTSTGT
                                    // are both authoritative zero-valued ecosystem pools.
                                    values[tail + 2] = 0;
                                    values[tail + 9] = 0;
                                    if (detailed_canopy_state) |canopy| {
                                        var total_leaf_area_m2: f64 = 0;
                                        for (0..config.plant_populations) |species| {
                                            const plant = cell * config.plant_populations + species;
                                            const branches = try canopy.branchRange(plant);
                                            for (canopy.branch_leaf_area_m2[branches.first..branches.end]) |leaf_area_m2| total_leaf_area_m2 += leaf_area_m2;
                                        }
                                        values[tail + 2] = total_leaf_area_m2 * area_inverse;
                                        values[tail + 9] =
                                            cell_litter_standing_dead_publication_state
                                                .standing_dead_carbon_g_c_by_cell[cell] *
                                            area_inverse;
                                    }
                                    values[tail + 3] = plant_flux.gross_primary_productivity_g_c * area_inverse;
                                    values[tail + 4] = plant_flux.autotrophic_respiration_g_c * area_inverse;
                                    values[tail + 5] = plant_flux.net_primary_productivity_g_c * area_inverse;
                                    values[tail + 6] = daily_heterotrophic_respiration.daily_signed_g_c[cell] * area_inverse;
                                    values[tail + 7] = (daily_soil_fire_methane_emission_g_c[cell] + daily_root_fire_methane_emission_g_c[cell]) * area_inverse;
                                    const active_carbonate = daily_soil_carbonate_mol_per_m3[first_layer .. first_layer + active_layers];
                                    const active_bicarbonate = daily_soil_bicarbonate_mol_per_m3[first_layer .. first_layer + active_layers];
                                    const active_carbonate_complexes = daily_soil_carbonate_complexes_mol_per_m3[first_layer .. first_layer + active_layers];
                                    const active_calcite = daily_soil_calcite_mol_per_m3[first_layer .. first_layer + active_layers];
                                    for (0..active_layers) |local_layer| {
                                        const layer = first_layer + local_layer;
                                        const aqueous = initial_chemistry_state.aqueous[layer];
                                        active_carbonate[local_layer] = aqueous.carbonate;
                                        active_bicarbonate[local_layer] = aqueous.bicarbonate;
                                        active_carbonate_complexes[local_layer] =
                                            aqueous.calcium_carbonate + aqueous.calcium_bicarbonate +
                                            aqueous.magnesium_carbonate + aqueous.magnesium_bicarbonate +
                                            aqueous.sodium_carbonate;
                                        active_calcite[local_layer] = initial_chemistry_state.geochemistry_solids[layer].calcite_solid_mol_per_m3;
                                    }
                                    const litter_aqueous = surface_litter_chemistry_state.cells[cell];
                                    values[tail + 8] = try ecosys.soil_inorganic_carbon_storage.calculate(.{
                                        .carbon_dioxide_and_methane_g_c = gas_inventory_g[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] +
                                            gas_inventory_g[@intFromEnum(ecosys.gas_transport.Species.methane)],
                                        .soil_water_volume_m3 = state.matrix_liquid_water_m3[first_layer .. first_layer + active_layers],
                                        .soil_bulk_volume_m3 = soil_solver_property_state.matrix_bulk_volume_m3[first_layer .. first_layer + active_layers],
                                        .carbonate_mol_per_m3 = active_carbonate,
                                        .bicarbonate_mol_per_m3 = active_bicarbonate,
                                        .dissolved_carbonate_complexes_mol_per_m3 = active_carbonate_complexes,
                                        .calcite_mol_per_m3 = active_calcite,
                                        .litter_water_volume_m3 = surface_precipitation_state.litter_water_m3[cell],
                                        .litter_carbonate_mol_per_m3 = litter_aqueous.carbonate_mol_per_m3,
                                        .litter_bicarbonate_mol_per_m3 = litter_aqueous.bicarbonate_mol_per_m3,
                                        .litter_dissolved_carbonate_complexes_mol_per_m3 = 0,
                                        .litter_calcite_mol_per_m3 = litter_aqueous.salt_minerals.calcite_mol_per_m3,
                                        .litter_bulk_volume_m3 = surface_litter_geometry_state.expanded_total_volume_m3[cell],
                                    }) * area_inverse;
                                    for (resolved.soil_enabled, values, 0..) |enabled, value, choice| if (enabled and !std.math.isFinite(value)) {
                                        const name = daily_soil_carbon_catalog.variables[choice].name;
                                        if (std.mem.startsWith(u8, name, "reserved_zero_")) {
                                            values[choice] = 0;
                                        } else {
                                            std.log.err("daily soil carbon choice is not yet backed by an authoritative translated owner: choice={d} name={s}", .{ choice + 1, name });
                                            return error.UntranslatedDailySoilCarbonChoice;
                                        }
                                    };
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, carbon_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_carbon_bank.streams[cell].write(file_name, daily_soil_carbon_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                }
                            }
                            if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(carbon_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const canopy = if (detailed_canopy_state) |*value| value else return error.DailyPlantCarbonOutputRequiresCanopyState;
                                const roots = if (plant_root_state) |*value| value else return error.DailyPlantCarbonOutputRequiresRootState;
                                for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                    const plant = try canopy.plantIndex(cell, species);
                                    const branches = try canopy.branchRange(plant);
                                    const values = try daily_plant_carbon_bank.row(plant);
                                    const root_profile = values[20 .. 20 + config.soil_layers];
                                    @memset(root_profile, 0);
                                    var root_symbiont_carbon_g: f64 = 0;
                                    for (0..config.soil_layers) |layer| {
                                        for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                                            const root = try roots.layerIndex(plant, domain, layer);
                                            root_profile[layer] += roots.mobile_carbon_g[root];
                                            root_symbiont_carbon_g +=
                                                roots.symbiont_structural_carbon_g_c[root] +
                                                roots.symbiont_mobile_carbon_g_c[root];
                                            for (0..roots.active_root_axis_count[plant]) |axis| {
                                                const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                                                root_profile[layer] += roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
                                            }
                                        }
                                    }
                                    var leaf_intermediate_carbon_g: f64 = 0;
                                    var canopy_symbiont_carbon_g: f64 = 0;
                                    var projected_leaf_area_m2: f64 = 0;
                                    for (branches.first..branches.end) |branch| {
                                        canopy_symbiont_carbon_g +=
                                            canopy.branch_symbiont_structural_carbon_g[branch] +
                                            canopy.branch_symbiont_mobile_carbon_g[branch];
                                        projected_leaf_area_m2 += canopy.branch_leaf_area_m2[branch];
                                        const nodes = try canopy.nodeRange(branch);
                                        for (nodes.first..nodes.end) |node| {
                                            leaf_intermediate_carbon_g +=
                                                canopy.node_c3_nonstructural_carbon_g[node] +
                                                canopy.node_c4_mesophyll_nonstructural_carbon_g[node] +
                                                canopy.node_bundle_sheath_co2_carbon_g[node] +
                                                canopy.node_bundle_sheath_bicarbonate_carbon_g[node];
                                        }
                                    }
                                    const primary_root_first = try roots.layerIndex(plant, 0, 0);
                                    const pools = try ecosys.plant_daily_pool_aggregation.calculateCarbonInto(.{
                                        .branch_leaf_carbon_g = canopy.branch_leaf_carbon_g[branches.first..branches.end],
                                        .branch_sheath_carbon_g = canopy.branch_sheath_carbon_g[branches.first..branches.end],
                                        .branch_stalk_carbon_g = canopy.branch_stalk_carbon_g[branches.first..branches.end],
                                        .branch_reserve_carbon_g = canopy.branch_reserve_carbon_g[branches.first..branches.end],
                                        .branch_husk_carbon_g = canopy.branch_husk_carbon_g[branches.first..branches.end],
                                        .branch_ear_carbon_g = canopy.branch_ear_carbon_g[branches.first..branches.end],
                                        .branch_grain_carbon_g = canopy.branch_grain_carbon_g[branches.first..branches.end],
                                        .branch_mobile_carbon_g = canopy.branch_mobile_carbon_g[branches.first..branches.end],
                                        .branch_seed_count = canopy.branch_seed_count[branches.first..branches.end],
                                        .c4_intermediate_carbon_g = leaf_intermediate_carbon_g,
                                        .canopy_symbiont_carbon_g = canopy_symbiont_carbon_g,
                                        .root_carbon_g_by_layer = root_profile,
                                        .primary_root_length_density_m_per_m3_by_layer = roots.root_length_density_m_per_m3[primary_root_first .. primary_root_first + config.soil_layers],
                                        .root_symbiont_carbon_g = root_symbiont_carbon_g,
                                        .standing_dead_carbon_g = canopy.plant_standing_dead_carbon_g[plant],
                                        .seed_storage_carbon_g = canopy.plant_seed_storage_carbon_g[plant],
                                        .projected_leaf_area_m2 = projected_leaf_area_m2,
                                        .plant_population_count = canopy.plant_population_count[plant],
                                    }, root_profile);
                                    const net_primary_productivity_g =
                                        plant_daily_flux_ledger.net_carbon_change_g[plant] +
                                        plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
                                    const carbon_balance_inputs: ecosys.plant_daily_output.CarbonBalanceInputs = .{
                                        .shoot_carbon_g = pools.shoot_carbon_g,
                                        .root_carbon_g = pools.root_carbon_g,
                                        .nodule_carbon_g = pools.nodule_carbon_g,
                                        .storage_carbon_g = pools.storage_carbon_g,
                                        .standing_dead_carbon_g = pools.vegetative_residue_carbon_g,
                                        .cumulative_carbon_sink_g = plant_daily_flux_ledger.carbon_sink_g[plant],
                                        .cumulative_root_soil_carbon_exchange_g = plant_daily_flux_ledger.root_soil_carbon_exchange_g[plant],
                                        .cumulative_carbon_balance_g = plant_daily_flux_ledger.cumulative_carbon_balance_g[plant],
                                        .cumulative_harvested_carbon_g = plant_daily_flux_ledger.cumulative_harvested_carbon_g[plant],
                                        .harvested_carbon_g = plant_daily_flux_ledger.harvested_carbon_g[plant],
                                        .carbon_oxidation_g = plant_daily_flux_ledger.carbon_oxidation_g[plant],
                                        .cumulative_net_primary_productivity_g = net_primary_productivity_g,
                                    };
                                    plant_balance_carbon_inputs_workspace[plant] = carbon_balance_inputs;
                                    const balance_carbon_g = try ecosys.plant_daily_output.carbonBalance(carbon_balance_inputs);
                                    try ecosys.plant_daily_output.calculateCarbonInto(.{
                                        .cell_area_m2 = canopy_cell_area_m2[cell],
                                        .shoot_carbon_g = pools.shoot_carbon_g,
                                        .leaf_carbon_g = pools.leaf_carbon_g,
                                        .sheath_carbon_g = pools.sheath_carbon_g,
                                        .stalk_carbon_g = pools.stalk_carbon_g,
                                        .reserve_carbon_g = pools.reserve_carbon_g,
                                        .husk_carbon_g = pools.husk_carbon_g,
                                        .ear_carbon_g = pools.ear_carbon_g,
                                        .grain_carbon_g = pools.grain_carbon_g,
                                        .root_carbon_g = pools.root_carbon_g,
                                        .nodule_carbon_g = pools.nodule_carbon_g,
                                        .vegetative_residue_carbon_g = pools.vegetative_residue_carbon_g,
                                        .grain_number = pools.grain_number,
                                        .projected_leaf_area_m2 = pools.projected_leaf_area_m2,
                                        .daily_net_carbon_change_g = plant_daily_flux_ledger.net_carbon_change_g[plant],
                                        .cumulative_carbon_uptake_g = plant_daily_flux_ledger.root_soil_carbon_exchange_g[plant],
                                        .cumulative_carbon_sink_g = plant_daily_flux_ledger.carbon_sink_g[plant],
                                        .initial_cumulative_carbon_sink_g = plant_daily_flux_ledger.initial_carbon_sink_g[plant],
                                        .signed_total_respiration_carbon_g = plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant],
                                        .signed_aboveground_respiration_carbon_g = plant_daily_flux_ledger.signed_aboveground_respiration_carbon_g[plant],
                                        .carbon_pollination_factor = canopy.plant_mobile_carbon_concentration_g_per_g[plant],
                                        .harvested_carbon_g = plant_daily_flux_ledger.harvested_carbon_g[plant],
                                        .root_length_density_m_per_m3_by_layer = root_profile,
                                        .plant_population = canopy.plant_population_count[plant],
                                        .balance_carbon_g = balance_carbon_g,
                                        .storage_carbon_g = pools.storage_carbon_g,
                                        .carbon_oxidation_flux_g = plant_daily_flux_ledger.carbon_oxidation_g[plant],
                                        .net_primary_productivity_g = net_primary_productivity_g,
                                        .canopy_height_m = development_canopy_height_m[plant],
                                    }, values);
                                    // Only assigned species produce output; an unassigned population
                                    // would emit all-zero rows, so no file is created for it.
                                    const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, carbon_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_plant_carbon_bank.streams[plant].write(file_name, daily_plant_carbon_catalog.variables, selection, resolved.plant_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                };
                            }
                        }
                        const editor_index: usize = 6;
                        const editor_name = scene.output_editors[editor_index];
                        if (!ecosys.delimited_input.isNo(editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedDailySoilWaterEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const liquid_fraction = daily_water_profile_workspace[0..config.soil_layers];
                                const ice_fraction = daily_water_profile_workspace[config.soil_layers .. 2 * config.soil_layers];
                                const potential_count = dailySoilWaterPotentialLayerCount(config.soil_layers);
                                const matric_plus_osmotic_potential = daily_water_profile_workspace[2 * config.soil_layers ..][0..potential_count];
                                for (0..state.cell_count) |cell| {
                                    var soil_water_storage_m3: f64 = 0;
                                    for (0..config.soil_layers) |local_layer| {
                                        const layer = try state.layerIndex(cell, local_layer);
                                        const layer_volume_m3 = soil_solver_property_state.layer_volume_m3[layer];
                                        if (!std.math.isFinite(layer_volume_m3) or layer_volume_m3 <= 0) return error.InvalidOutputSoilLayerVolume;
                                        liquid_fraction[local_layer] = state.liquid_water_m3[layer] / layer_volume_m3;
                                        ice_fraction[local_layer] = state.ice_water_m3[layer] / layer_volume_m3;
                                        soil_water_storage_m3 += state.liquid_water_m3[layer] +
                                            state.water_vapor_volume_m3[layer] +
                                            state.ice_water_m3[layer];
                                        if (local_layer < potential_count) matric_plus_osmotic_potential[local_layer] = soil_hourly_workspace.matric_plus_osmotic_potential_megapascal[layer];
                                    }
                                    const snow_base = cell * snow_transport_state.layer_capacity;
                                    for (0..snow_transport_state.layer_capacity) |snow_layer| {
                                        const snow = snow_base + snow_layer;
                                        soil_water_storage_m3 += snow_transport_state.solid_snow_water_equivalent_m3[snow] +
                                            snow_transport_state.liquid_water_volume_m3[snow] +
                                            snow_transport_state.vapor_water_equivalent_m3[snow] +
                                            snow_transport_state.ice_volume_m3[snow] * runscript.snow_ice_density_megagrams_per_m3;
                                    }
                                    soil_water_storage_m3 +=
                                        surface_precipitation_state.litter_water_m3[cell] +
                                        surface_litter_ice_m3[cell];
                                    for (0..config.plant_populations) |species| {
                                        const plant = cell * config.plant_populations + species;
                                        soil_water_storage_m3 += plant_state.canopy_water_storage_m_per_m2[plant] * canopy_cell_area_m2[cell];
                                        if (canopy_precipitation_retention_state) |retention| {
                                            soil_water_storage_m3 += retention.living_surface_water_m3[plant] + retention.standing_dead_surface_water_m3[plant];
                                        }
                                    }
                                    const litter_capacity_m3 = surface_precipitation_state.litter_water_capacity_m3[cell];
                                    const surface_litter_thickness_m =
                                        surface_litter_geometry_state.expanded_total_volume_m3[cell] /
                                        canopy_cell_area_m2[cell];
                                    if (!std.math.isFinite(surface_litter_thickness_m) or
                                        surface_litter_thickness_m < 0)
                                        return error.InvalidOutputSurfaceLitterThickness;
                                    const active_boundary_base = cell * (soil_geometry_state.layer_capacity + 1);
                                    const active_surface_depth_m = soil_geometry_state.boundary_depth_m[active_boundary_base + soil_geometry_state.first_active_layer[cell]];
                                    const values = try daily_soil_water_bank.row(cell);
                                    try ecosys.soil_daily_output.calculateWaterInto(.{
                                        .cell_area_m2 = canopy_cell_area_m2[cell],
                                        .grid_area_m2 = total_grid_area_m2,
                                        .rainfall_m3 = daily_water_ledger.rainfall_m3[cell],
                                        .evaporation_m3 = daily_water_ledger.evaporation_m3[cell],
                                        .runoff_m3 = daily_water_ledger.runoff_m3[cell],
                                        .soil_water_storage_m3 = soil_water_storage_m3,
                                        .outflow_m3 = daily_water_ledger.water_outflow_m3[cell],
                                        .snow_depth_m = snow_depth_m[cell],
                                        .liquid_water_fraction_by_layer = liquid_fraction,
                                        .surface_excess_liquid_water_depth_m = @max(
                                            0,
                                            (surface_precipitation_state.litter_water_m3[cell] -
                                                litter_capacity_m3) /
                                                canopy_cell_area_m2[cell],
                                        ),
                                        .ice_fraction_by_layer = ice_fraction,
                                        .surface_excess_ice_water_depth_m = @max(
                                            0,
                                            (surface_litter_ice_m3[cell] -
                                                litter_capacity_m3) /
                                                canopy_cell_area_m2[cell],
                                        ),
                                        .matric_plus_osmotic_potential_mpa_by_layer = matric_plus_osmotic_potential,
                                        .surface_matric_potential_megapascal = surface_litter_water_environment_state.matric_water_potential_megapascal[cell],
                                        .lateral_water_outflow_m3 = daily_water_ledger.lateral_water_outflow_m3[cell],
                                        .sediment_outflow_m3 = daily_water_ledger.sediment_outflow_m3[cell],
                                        .mineral_soil_surface_depth_m = active_surface_depth_m,
                                        .surface_litter_thickness_m = surface_litter_thickness_m,
                                        .active_layer_bottom_depth_m = soil_boundary_topology_state.active_layer_depth_m[cell],
                                        .water_table_depth_m = soil_boundary_topology_state.internal_water_table_depth_m[cell],
                                    }, values);
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_water_bank.streams[cell].write(file_name, daily_soil_water_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        if (!ecosys.delimited_input.isNo(editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedDailyPlantWaterEditor;
                            if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const growth = if (plant_growth_stage_state) |*value| value else return error.DailyPlantWaterOutputRequiresGrowthStageState;
                                const dormancy = if (plant_dormancy_state) |*value| value else return error.DailyPlantWaterOutputRequiresDormancyState;
                                for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                    const plant = cell * config.plant_populations + species;
                                    const output = try ecosys.plant_daily_output.water(
                                        plant_daily_flux_ledger.transpiration_source_m3[plant],
                                        try ecosys.plant_development.coldOrWaterStressHours(growth, dormancy, plant),
                                        phenology_root_oxygen_fraction[plant],
                                        plant_state.canopy_water_storage_m_per_m2[plant] * canopy_cell_area_m2[cell],
                                        canopy_cell_area_m2[cell],
                                    );
                                    const values = try daily_plant_water_bank.row(plant);
                                    values[0] = output.transpiration_mm;
                                    values[1] = output.cold_or_water_stress_h;
                                    values[2] = output.oxygen_stress_factor;
                                    values[3] = output.plant_water_storage_mm;
                                    // Only assigned species produce output; an unassigned population
                                    // would emit all-zero rows, so no file is created for it.
                                    const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_plant_water_bank.streams[plant].write(file_name, daily_plant_water_catalog.variables, selection, resolved.plant_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                };
                            }
                        }
                        const soil_nitrogen_editor_index: usize = 7;
                        const soil_nitrogen_editor_name = scene.output_editors[soil_nitrogen_editor_index];
                        if (!ecosys.delimited_input.isNo(soil_nitrogen_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + soil_nitrogen_editor_index] orelse return error.UnresolvedDailySoilNitrogenEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(soil_nitrogen_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                for (0..state.cell_count) |cell| {
                                    const values = try daily_soil_nitrogen_bank.row(cell);
                                    @memset(values, std.math.nan(f64));
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const pools = try ecosys.soil_daily_nitrogen_pools.calculate(
                                        &soil_organic_state,
                                        &surface_organic_state,
                                        first_layer,
                                        active_layers,
                                        cell,
                                    );
                                    const area_inverse = 1.0 / canopy_cell_area_m2[cell];
                                    values[0] = pools.residue_nitrogen_g_n * area_inverse;
                                    values[1] = pools.humus_nitrogen_g_n * area_inverse;
                                    values[2] = (soil_fertilizer_inventory.daily_nitrogen_input_g_n[cell] + daily_organic_fertilizer_nitrogen_input_g_n[cell]) * area_inverse;
                                    var plant_nitrogen_sink_g_n: f64 = 0;
                                    var harvested_nitrogen_g_n: f64 = 0;
                                    for (0..config.plant_populations) |species| {
                                        const plant = cell * config.plant_populations + species;
                                        plant_nitrogen_sink_g_n += plant_daily_flux_ledger.nitrogen_sink_g[plant];
                                        harvested_nitrogen_g_n += plant_daily_flux_ledger.harvested_nitrogen_g[plant];
                                    }
                                    values[3] = plant_nitrogen_sink_g_n * area_inverse;
                                    values[6] = daily_nitrogen_export.dissolved_organic_nitrogen_runoff_g_n[cell] / total_grid_area_m2;
                                    values[7] = daily_nitrogen_export.dissolved_organic_nitrogen_drainage_g_n[cell] / total_grid_area_m2;
                                    values[8] = daily_nitrogen_export.dissolved_inorganic_nitrogen_runoff_g_n[cell] / total_grid_area_m2;
                                    values[9] = daily_nitrogen_export.dissolved_inorganic_nitrogen_drainage_g_n[cell] / total_grid_area_m2;
                                    var dissolved_dinitrogen_g_n: f64 = 0;
                                    for (0..active_layers) |local_layer| {
                                        const layer = first_layer + local_layer;
                                        for (0..ecosys.gas_transport.species_count) |species| {
                                            const component = layer * ecosys.gas_transport.species_count + species;
                                            if (species == @intFromEnum(ecosys.gas_transport.Species.nitrogen)) {
                                                dissolved_dinitrogen_g_n +=
                                                    gas_transport_state.dissolved_mass_g[component] +
                                                    gas_transport_state.macropore_dissolved_mass_g[component] +
                                                    gas_transport_state.band_dissolved_mass_g[component];
                                            }
                                        }
                                    }
                                    for (0..ecosys.gas_transport.species_count) |species| {
                                        const component = cell * ecosys.gas_transport.species_count + species;
                                        if (species == @intFromEnum(ecosys.gas_transport.Species.nitrogen)) {
                                            dissolved_dinitrogen_g_n +=
                                                litter_gas_transport_state.dissolved_mass_g[component] +
                                                litter_gas_transport_state.macropore_dissolved_mass_g[component] +
                                                litter_gas_transport_state.band_dissolved_mass_g[component];
                                        }
                                    }
                                    values[10] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .nitrous_oxide)) * area_inverse;
                                    values[11] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .ammonia)) * area_inverse;
                                    values[12] = dissolved_dinitrogen_g_n * area_inverse;
                                    values[13] = (pools.residue_nitrogen_g_n + pools.humus_nitrogen_g_n) * area_inverse;
                                    const ammonium_profile = values[14 .. 14 + config.soil_layers];
                                    const oxidized_profile = values[14 + config.soil_layers .. 14 + 2 * config.soil_layers];
                                    const mineral = try ecosys.soil_daily_mineral_nitrogen.calculateInto(
                                        &initial_chemistry_state,
                                        &soil_reactive_nitrogen_state,
                                        first_layer,
                                        active_layers,
                                        state.matrix_liquid_water_m3,
                                        soil_solver_property_state.matrix_bulk_volume_m3,
                                        soil_solver_property_state.bulk_density_megagrams_per_m3,
                                        .{
                                            .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                            .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                            .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                            .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                        },
                                        runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                                        config.negligible_quantity_threshold,
                                        ammonium_profile,
                                        oxidized_profile,
                                    );
                                    values[4] = mineral.ammonium_nitrogen_g_n * area_inverse;
                                    values[5] = mineral.nitrate_plus_nitrite_nitrogen_g_n * area_inverse;
                                    const tail = 14 + 2 * config.soil_layers;
                                    const litter = surface_litter_chemistry_state.cells[cell];
                                    const litter_bulk_volume_m3 = surface_litter_geometry_state.expanded_total_volume_m3[cell];
                                    const litter_ammonium_g_n = runscript.fertilizer_nitrogen_molar_mass_g_per_mol *
                                        (litter.ammonium_mol_per_m3 * surface_precipitation_state.litter_water_m3[cell] +
                                            litter.exchange.ammonium_mol_per_megagram * surface_litter_geometry_state.dry_mass_megagrams[cell]);
                                    values[tail] = if (litter_bulk_volume_m3 > config.negligible_quantity_threshold) litter_ammonium_g_n / litter_bulk_volume_m3 else 0;
                                    values[tail + 1] = daily_soil_combusted_nitrogen_g_n[cell] * area_inverse;
                                    values[tail + 2] = harvested_nitrogen_g_n * area_inverse;
                                    values[tail + 3] = daily_microbial_nitrogen_mineralization_g_n[cell] * area_inverse;
                                    values[tail + 4] = (daily_soil_fire_nitrogen_flux_g_n[cell] + daily_root_fire_nitrogen_flux_g_n[cell]) * area_inverse;
                                    values[tail + 5] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .nitrogen)) * area_inverse;
                                    for (resolved.soil_enabled, values, 0..) |enabled, value, choice| if (enabled and !std.math.isFinite(value)) {
                                        std.log.err("daily soil nitrogen choice has non-finite authoritative state: choice={d} name={s}", .{ choice + 1, daily_soil_nitrogen_catalog.variables[choice].name });
                                        return error.NonFiniteDailySoilNitrogenChoice;
                                    };
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, soil_nitrogen_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_nitrogen_bank.streams[cell].write(file_name, daily_soil_nitrogen_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        const soil_phosphorus_editor_index: usize = 8;
                        const soil_phosphorus_editor_name = scene.output_editors[soil_phosphorus_editor_index];
                        if (!ecosys.delimited_input.isNo(soil_phosphorus_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + soil_phosphorus_editor_index] orelse return error.UnresolvedDailySoilPhosphorusEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(soil_phosphorus_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                for (0..state.cell_count) |cell| {
                                    const values = try daily_soil_phosphorus_bank.row(cell);
                                    @memset(values, std.math.nan(f64));
                                    const area_inverse = 1.0 / canopy_cell_area_m2[cell];
                                    values[2] = (mineral_fertilizer_inventory.daily_phosphorus_input_g_p[cell] + daily_organic_fertilizer_phosphorus_input_g_p[cell]) * area_inverse;
                                    var plant_litterfall_phosphorus_g_p: f64 = 0;
                                    for (0..config.plant_populations) |species| plant_litterfall_phosphorus_g_p += plant_daily_flux_ledger.phosphorus_sink_g[cell * config.plant_populations + species];
                                    values[3] = plant_litterfall_phosphorus_g_p * area_inverse;
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const organic_pools = try ecosys.soil_daily_phosphorus_pools.calculate(
                                        &soil_organic_state,
                                        &surface_organic_state,
                                        first_layer,
                                        active_layers,
                                        cell,
                                    );
                                    values[0] = organic_pools.residue_phosphorus_g_p * area_inverse;
                                    values[1] = organic_pools.humus_phosphorus_g_p * area_inverse;
                                    values[5] = daily_phosphorus_export.dissolved_organic_phosphorus_runoff_g_p[cell] / total_grid_area_m2;
                                    values[6] = daily_phosphorus_export.dissolved_organic_phosphorus_drainage_g_p[cell] / total_grid_area_m2;
                                    values[7] = daily_phosphorus_export.dissolved_inorganic_phosphorus_runoff_g_p[cell] / total_grid_area_m2;
                                    values[8] = daily_phosphorus_export.dissolved_inorganic_phosphorus_drainage_g_p[cell] / total_grid_area_m2;
                                    const phosphate_band_fraction = runscript.plant_nutrient_initialization.initial_phosphate_band_fraction;
                                    const phosphate_non_band_fraction = 1.0 - phosphate_band_fraction;
                                    values[4] = try ecosys.soil_phosphate_inventory.exchangeablePhosphorus_g_p(
                                        &state,
                                        &initial_chemistry_state,
                                        &surface_litter_chemistry_state,
                                        soil_solver_property_state.matrix_bulk_volume_m3,
                                        soil_solver_property_state.bulk_density_megagrams_per_m3,
                                        surface_litter_geometry_state.dry_mass_megagrams,
                                        cell,
                                        phosphate_non_band_fraction,
                                        phosphate_band_fraction,
                                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                    ) * area_inverse;
                                    values[9] = try ecosys.soil_phosphate_inventory.precipitatedPhosphorus_g_p(
                                        &state,
                                        &initial_chemistry_state,
                                        &surface_litter_chemistry_state,
                                        state.matrix_liquid_water_m3,
                                        surface_precipitation_state.litter_water_m3,
                                        cell,
                                        phosphate_non_band_fraction,
                                        phosphate_band_fraction,
                                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                    ) * area_inverse;
                                    values[10] = organic_pools.microbial_phosphorus_g_p * area_inverse;
                                    values[11] = (try ecosys.soil_daily_fire_phosphorus.sourceSignedFlux_g_p(
                                        cell,
                                        daily_soil_fire_phosphorus_flux_g_p,
                                        plant_daily_flux_ledger.phosphorus_oxidation_g,
                                        config.plant_populations,
                                    ) + daily_root_fire_phosphorus_flux_g_p[cell]) * area_inverse;
                                    for (0..config.soil_layers) |local_layer| {
                                        const aqueous_index = 12 + local_layer;
                                        const sorbed_index = 12 + config.soil_layers + local_layer;
                                        if (local_layer >= active_layers) {
                                            values[aqueous_index] = 0;
                                            values[sorbed_index] = 0;
                                            continue;
                                        }
                                        const layer = first_layer + local_layer;
                                        const non_band = initial_chemistry_state.non_band_phosphate[layer];
                                        const band = initial_chemistry_state.band_phosphate[layer];
                                        values[aqueous_index] = 31.0 * (phosphate_non_band_fraction * (non_band.dissolved_hpo4_mol_p_per_m3 + non_band.dissolved_h2po4_mol_p_per_m3) + phosphate_band_fraction * (band.dissolved_hpo4_mol_p_per_m3 + band.dissolved_h2po4_mol_p_per_m3));
                                        values[sorbed_index] = try ecosys.soil_daily_output.sorbedPhosphorusConcentrationGPerM3(
                                            phosphate_non_band_fraction *
                                                (non_band.adsorbed_hpo4_mol_p_per_megagram +
                                                    non_band.adsorbed_h2po4_mol_p_per_megagram) +
                                                phosphate_band_fraction *
                                                    (band.adsorbed_hpo4_mol_p_per_megagram +
                                                        band.adsorbed_h2po4_mol_p_per_megagram),
                                            soil_solver_property_state.bulk_density_megagrams_per_m3[layer],
                                            runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                        );
                                    }
                                    const after_layers = 12 + 2 * config.soil_layers;
                                    const litter = surface_litter_chemistry_state.cells[cell];
                                    values[after_layers] = 31.0 * (litter.hpo4_mol_p_per_m3 + litter.h2po4_mol_p_per_m3);
                                    values[after_layers + 1] = if (surface_litter_geometry_state.expanded_total_volume_m3[cell] >
                                        config.negligible_quantity_threshold)
                                        try ecosys.soil_daily_output.sorbedPhosphorusConcentrationGPerM3(
                                            litter.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram +
                                                litter.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram,
                                            surface_litter_geometry_state.dry_mass_megagrams[cell] /
                                                surface_litter_geometry_state.expanded_total_volume_m3[cell],
                                            runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                        )
                                    else
                                        0;
                                    const final_fixed = after_layers + 2;
                                    values[final_fixed] = daily_soil_combusted_phosphorus_g_p[cell] * area_inverse;
                                    values[final_fixed + 1] = try ecosys.soil_phosphate_inventory.solublePhosphorus_g_p(
                                        &state,
                                        &initial_chemistry_state,
                                        &surface_litter_chemistry_state,
                                        state.matrix_liquid_water_m3,
                                        surface_precipitation_state.litter_water_m3,
                                        cell,
                                        phosphate_non_band_fraction,
                                        phosphate_band_fraction,
                                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                    ) * area_inverse;
                                    var harvested_phosphorus_g_p: f64 = 0;
                                    for (0..config.plant_populations) |species| harvested_phosphorus_g_p += plant_daily_flux_ledger.harvested_phosphorus_g[cell * config.plant_populations + species];
                                    values[final_fixed + 2] = harvested_phosphorus_g_p * area_inverse;
                                    values[final_fixed + 3] = daily_microbial_phosphate_mineralization_g_p[cell] * area_inverse;
                                    // OUTSD phosphorus choices 49 and 50 are source-reserved.
                                    values[values.len - 2] = 0;
                                    values[values.len - 1] = 0;
                                    for (resolved.soil_enabled, values, 0..) |enabled, value, choice| if (enabled and !std.math.isFinite(value)) {
                                        const name = daily_soil_phosphorus_catalog.variables[choice].name;
                                        if (std.mem.startsWith(u8, name, "reserved_zero_")) {
                                            values[choice] = 0;
                                        } else {
                                            std.log.err("daily soil phosphorus choice is not yet backed by an authoritative translated owner: choice={d} name={s}", .{ choice + 1, name });
                                            return error.UntranslatedDailySoilPhosphorusChoice;
                                        }
                                    };
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, soil_phosphorus_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_phosphorus_bank.streams[cell].write(file_name, daily_soil_phosphorus_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        const soil_heat_editor_index: usize = 9;
                        const soil_heat_editor_name = scene.output_editors[soil_heat_editor_index];
                        if (!ecosys.delimited_input.isNo(soil_heat_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + soil_heat_editor_index] orelse return error.UnresolvedDailySoilHeatEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(soil_heat_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const maximum_temperature_c = daily_soil_heat_workspace[0..config.soil_layers];
                                const minimum_temperature_c = daily_soil_heat_workspace[config.soil_layers .. 2 * config.soil_layers];
                                const electrical_conductivity_dS_per_m = daily_soil_heat_workspace[2 * config.soil_layers .. 3 * config.soil_layers];
                                const charge_fractions = ecosys.solute_charge_classification.ZoneFractions{
                                    .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                    .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                    .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                    .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                    .phosphate_non_band = 1 - runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                                    .phosphate_band = runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                                };
                                for (0..state.cell_count) |cell| {
                                    try daily_heat_ledger.validateCell(cell);
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const ledger_first = cell * config.soil_layers;
                                    @memset(maximum_temperature_c, 0);
                                    @memset(minimum_temperature_c, 0);
                                    @memset(electrical_conductivity_dS_per_m, 0);
                                    @memcpy(maximum_temperature_c[0..active_layers], daily_heat_ledger.maximum_soil_temperature_c[ledger_first .. ledger_first + active_layers]);
                                    @memcpy(minimum_temperature_c[0..active_layers], daily_heat_ledger.minimum_soil_temperature_c[ledger_first .. ledger_first + active_layers]);
                                    for (0..active_layers) |local_layer| {
                                        const coefficients = try initial_chemistry_state.activityCoefficients(first_layer + local_layer, charge_fractions);
                                        electrical_conductivity_dS_per_m[local_layer] = coefficients.electrical_conductivity_dS_per_m;
                                    }
                                    const values = try daily_soil_heat_bank.row(cell);
                                    try ecosys.soil_daily_output.calculateHeatInto(.{
                                        .total_radiation_megajoules_m2 = daily_heat_ledger.total_radiation_megajoules_per_m2[cell],
                                        .maximum_air_temperature_c = daily_heat_ledger.maximum_air_temperature_c[cell],
                                        .minimum_air_temperature_c = daily_heat_ledger.minimum_air_temperature_c[cell],
                                        .maximum_atmospheric_vapor_pressure_kpa = daily_heat_ledger.maximum_vapor_pressure_kpa[cell],
                                        .minimum_atmospheric_vapor_pressure_kpa = daily_heat_ledger.minimum_vapor_pressure_kpa[cell],
                                        .cumulative_wind_m = daily_heat_ledger.cumulative_wind_distance_m[cell],
                                        .total_precipitation_mm = daily_heat_ledger.total_precipitation_mm[cell],
                                        .maximum_soil_temperature_c_by_layer = maximum_temperature_c,
                                        .minimum_soil_temperature_c_by_layer = minimum_temperature_c,
                                        .surface_maximum_soil_temperature_c = daily_heat_ledger.surface_maximum_temperature_c[cell],
                                        .surface_minimum_soil_temperature_c = daily_heat_ledger.surface_minimum_temperature_c[cell],
                                        .electrical_conductivity_by_layer = electrical_conductivity_dS_per_m,
                                        .ionic_outflow_mol = daily_heat_ledger.ionic_outflow_mol[cell],
                                        .grid_area_m2 = total_grid_area_m2,
                                    }, values);
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .soil_or_eco, current_year, soil_heat_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_heat_bank.streams[cell].write(file_name, daily_soil_heat_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        inline for (.{ DailyPlantElement.nitrogen, DailyPlantElement.phosphorus }) |element| {
                            const nutrient_editor_index: usize = switch (element) {
                                .nitrogen => 7,
                                .phosphorus => 8,
                            };
                            const nutrient_editor_name = scene.output_editors[nutrient_editor_index];
                            if (!ecosys.delimited_input.isNo(nutrient_editor_name)) {
                                const resolved = resolved_output_editors[pass.scene_index * 10 + nutrient_editor_index] orelse return error.UnresolvedDailyPlantNutrientEditor;
                                if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                    const selection_index = output_selection_catalog.find(nutrient_editor_name) orelse return error.MissingOutputSelection;
                                    const selection = output_selection_catalog.entries.items[selection_index].selection;
                                    const canopy = if (detailed_canopy_state) |*value| value else return error.DailyPlantNutrientOutputRequiresCanopyState;
                                    const roots = if (plant_root_state) |*value| value else return error.DailyPlantNutrientOutputRequiresRootState;
                                    for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                        const plant = try canopy.plantIndex(cell, species);
                                        const branches = try canopy.branchRange(plant);
                                        const root_workspace = daily_element_root_workspace[plant * config.soil_layers .. (plant + 1) * config.soil_layers];
                                        const pools = try calculateDailyPlantElementPools(
                                            canopy,
                                            roots,
                                            plant,
                                            branches,
                                            element,
                                            root_metabolism_plant_parameters[plant].biologicalDomainCount(),
                                            root_workspace,
                                        );
                                        var leaf_carbon_g: f64 = 0;
                                        var leaf_nitrogen_g: f64 = 0;
                                        var leaf_phosphorus_g: f64 = 0;
                                        for (branches.first..branches.end) |branch| {
                                            leaf_carbon_g += canopy.branch_leaf_carbon_g[branch];
                                            leaf_nitrogen_g += canopy.branch_leaf_nitrogen_g[branch];
                                            leaf_phosphorus_g += canopy.branch_leaf_phosphorus_g[branch];
                                        }
                                        const balance_inputs: ecosys.plant_daily_output.NutrientBalanceInputs = switch (element) {
                                            .nitrogen => .{
                                                .shoot_g = pools.shoot_g,
                                                .root_g = pools.root_g,
                                                .nodule_g = pools.nodule_g,
                                                .storage_g = pools.storage_g,
                                                .standing_dead_g = pools.vegetative_residue_g,
                                                .cumulative_sink_g = plant_daily_flux_ledger.nitrogen_sink_g[plant],
                                                .cumulative_root_soil_exchange_g = plant_daily_flux_ledger.root_soil_nitrogen_exchange_g[plant],
                                                .cumulative_balance_g = plant_daily_flux_ledger.cumulative_nitrogen_balance_g[plant],
                                                .cumulative_harvested_g = plant_daily_flux_ledger.cumulative_harvested_nitrogen_g[plant],
                                                .harvested_g = plant_daily_flux_ledger.harvested_nitrogen_g[plant],
                                                .oxidation_g = plant_daily_flux_ledger.nitrogen_oxidation_g[plant],
                                                .atmospheric_exchange_g = plant_daily_flux_ledger.ammonia_exchange_g_n[plant],
                                                .biological_fixation_g = plant_daily_flux_ledger.symbiotic_nitrogen_fixation_g[plant],
                                            },
                                            .phosphorus => .{
                                                .shoot_g = pools.shoot_g,
                                                .root_g = pools.root_g,
                                                .nodule_g = pools.nodule_g,
                                                .storage_g = pools.storage_g,
                                                .standing_dead_g = pools.vegetative_residue_g,
                                                .cumulative_sink_g = plant_daily_flux_ledger.phosphorus_sink_g[plant],
                                                .cumulative_root_soil_exchange_g = plant_daily_flux_ledger.root_soil_phosphorus_exchange_g[plant],
                                                .cumulative_balance_g = plant_daily_flux_ledger.cumulative_phosphorus_balance_g[plant],
                                                .cumulative_harvested_g = plant_daily_flux_ledger.cumulative_harvested_phosphorus_g[plant],
                                                .harvested_g = plant_daily_flux_ledger.harvested_phosphorus_g[plant],
                                                .oxidation_g = plant_daily_flux_ledger.phosphorus_oxidation_g[plant],
                                            },
                                        };
                                        switch (element) {
                                            .nitrogen => plant_balance_nitrogen_inputs_workspace[plant] = balance_inputs,
                                            .phosphorus => plant_balance_phosphorus_inputs_workspace[plant] = balance_inputs,
                                        }
                                        const nutrient_inputs: ecosys.plant_daily_output.NutrientInputs = .{
                                            .cell_area_m2 = canopy_cell_area_m2[cell],
                                            .shoot_g = pools.shoot_g,
                                            .leaf_g = pools.leaf_g,
                                            .sheath_g = pools.sheath_g,
                                            .stalk_g = pools.stalk_g,
                                            .reserve_g = pools.reserve_g,
                                            .husk_g = pools.husk_g,
                                            .ear_g = pools.ear_g,
                                            .grain_g = pools.grain_g,
                                            .root_g = pools.root_g,
                                            .nodule_g = pools.nodule_g,
                                            .vegetative_residue_g = pools.vegetative_residue_g,
                                            .cumulative_uptake_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.root_soil_nitrogen_exchange_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.root_soil_phosphorus_exchange_g[plant],
                                            },
                                            .cumulative_sink_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.nitrogen_sink_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.phosphorus_sink_g[plant],
                                            },
                                            .nitrogen_fixation_g = plant_daily_flux_ledger.symbiotic_nitrogen_fixation_g[plant],
                                            .aboveground_litter_sink_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.initial_nitrogen_sink_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.initial_phosphorus_sink_g[plant],
                                            },
                                            .pollination_factor = switch (element) {
                                                .nitrogen => canopy.plant_mobile_nitrogen_concentration_g_per_g[plant],
                                                .phosphorus => canopy.plant_mobile_phosphorus_concentration_g_per_g[plant],
                                            },
                                            .leaf_carbon_g = leaf_carbon_g,
                                            .leaf_nitrogen_g = leaf_nitrogen_g,
                                            .leaf_phosphorus_g = leaf_phosphorus_g,
                                            .nonstructural_carbon_g = canopy.plant_mobile_carbon_g[plant],
                                            .nonstructural_nitrogen_g = canopy.plant_mobile_nitrogen_g[plant],
                                            .nonstructural_phosphorus_g = canopy.plant_mobile_phosphorus_g[plant],
                                            .minimum_leaf_carbon_g = runscript.plant_pool_parameters.branch_structural_presence_g_per_plant * canopy.plant_population_count[plant],
                                            .other_pollination_factor = canopy.plant_mobile_phosphorus_concentration_g_per_g[plant],
                                            .gaseous_loss_g = plant_daily_flux_ledger.ammonia_exchange_g_n[plant],
                                            .harvested_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.harvested_nitrogen_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.harvested_phosphorus_g[plant],
                                            },
                                            .balance_g = switch (element) {
                                                .nitrogen => try ecosys.plant_daily_output.nitrogenBalance(balance_inputs),
                                                .phosphorus => try ecosys.plant_daily_output.phosphorusBalance(balance_inputs),
                                            },
                                            .storage_g = pools.storage_g,
                                            .oxidation_flux_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.nitrogen_oxidation_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.phosphorus_oxidation_g[plant],
                                            },
                                        };
                                        const values = switch (element) {
                                            .nitrogen => try daily_plant_nitrogen_bank.row(plant),
                                            .phosphorus => try daily_plant_phosphorus_bank.row(plant),
                                        };
                                        switch (element) {
                                            .nitrogen => @memcpy(values, &(try ecosys.plant_daily_output.nitrogen(nutrient_inputs)).values),
                                            .phosphorus => @memcpy(values, &(try ecosys.plant_daily_output.phosphorus(nutrient_inputs)).values),
                                        }
                                        // Only assigned species produce output; an unassigned population
                                        // would emit all-zero rows, so no file is created for it.
                                        const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                        const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, nutrient_editor_name);
                                        defer allocator.free(file_name);
                                        _ = switch (element) {
                                            .nitrogen => try daily_plant_nitrogen_bank.streams[plant].write(file_name, daily_plant_nitrogen_catalog.variables, selection, resolved.plant_enabled, .{
                                                .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                                .values = values,
                                            }),
                                            .phosphorus => try daily_plant_phosphorus_bank.streams[plant].write(file_name, daily_plant_phosphorus_catalog.variables, selection, resolved.plant_enabled, .{
                                                .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                                .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                                .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                                .values = values,
                                            }),
                                        };
                                    };
                                }
                            }
                        }
                        // GROSUB 13071–13165 / EXTRACT 949–951: per-plant C/N/P balance tuple
                        // and cell totals.
                        if (detailed_canopy_state != null and
                            plant_root_state != null and
                            plant_phenology_state != null)
                        {
                            const canopy = &detailed_canopy_state.?;
                            const roots = &plant_root_state.?;
                            try assemblePlantBalanceInputs(
                                canopy,
                                roots,
                                root_metabolism_plant_parameters,
                                &plant_daily_flux_ledger,
                                plant_phenology_state.?.active,
                                plant_balance_carbon_inputs_workspace,
                                plant_balance_nitrogen_inputs_workspace,
                                plant_balance_phosphorus_inputs_workspace,
                                plant_balance_carbon_root_by_layer,
                                plant_balance_carbon_root_length_density_by_layer,
                                plant_balance_nutrient_root_by_layer,
                            );
                            try ecosys.plant_balance_publication.refresh(
                                &plant_balance_publication_state,
                                .{
                                    .active_by_plant = plant_phenology_state.?.active,
                                    .carbon_by_plant = plant_balance_carbon_inputs_workspace,
                                    .nitrogen_by_plant = plant_balance_nitrogen_inputs_workspace,
                                    .phosphorus_by_plant = plant_balance_phosphorus_inputs_workspace,
                                },
                            );
                        }
                        const development_editor_index: usize = 9;
                        const development_editor_name = scene.output_editors[development_editor_index];
                        if (!ecosys.delimited_input.isNo(development_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + development_editor_index] orelse return error.UnresolvedDailyPlantDevelopmentEditor;
                            if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(development_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const canopy = if (detailed_canopy_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresCanopyState;
                                const growth = if (plant_growth_stage_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresGrowthStageState;
                                const development_state = if (branch_development_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresBranchDevelopmentState;
                                const phenology = if (plant_phenology_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresPhenologyState;
                                for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                    const plant = cell * config.plant_populations + species;
                                    const branches = try growth.branchRange(plant);
                                    const has_branch = branches.first < branches.end;
                                    const main_branch = (try growth.mainLivingBranch(plant)) orelse branches.first;
                                    if (has_branch and (main_branch >= canopy.branch_leaf_carbon_g.len or main_branch >= development_state.branch_count)) return error.DailyPlantDevelopmentBranchDimensionMismatch;
                                    const stage: ecosys.plant_growth_stages.BranchState = if (has_branch) growth.branches[main_branch] else .{};
                                    const milestones = [_]u32{
                                        stage.emergence_day,
                                        stage.floral_initiation_day,
                                        stage.stem_elongation_start_day,
                                        stage.stem_elongation_midpoint_day,
                                        stage.stem_elongation_end_day,
                                        stage.anthesis_day,
                                        stage.grain_fill_start_day,
                                        stage.seed_number_set_end_day,
                                        stage.seed_size_set_end_day,
                                        if (has_branch) development_state.stage_day[main_branch * 10 + 9] else 0,
                                    };
                                    var leaf_carbon_g: f64 = 0;
                                    var leaf_nitrogen_g: f64 = 0;
                                    var leaf_phosphorus_g: f64 = 0;
                                    for (branches.first..branches.end) |branch| {
                                        if (growth.branches[branch].dead) continue;
                                        leaf_carbon_g += canopy.branch_leaf_carbon_g[branch];
                                        leaf_nitrogen_g += canopy.branch_leaf_nitrogen_g[branch];
                                        leaf_phosphorus_g += canopy.branch_leaf_phosphorus_g[branch];
                                    }
                                    const output = try ecosys.plant_daily_output.development(.{
                                        .alive = has_branch and phenology.active[plant],
                                        .milestone_day_by_stage = &milestones,
                                        .branch_count = branches.end - branches.first,
                                        .main_branch_stage = stage.accumulated_vegetative_stage,
                                        .development_feedback = canopy.plant_nitrogen_phosphorus_fixation_constraint_fraction[plant],
                                        .leaf_structural_carbon_g = leaf_carbon_g,
                                        .leaf_structural_nitrogen_g = leaf_nitrogen_g,
                                        .leaf_structural_phosphorus_g = leaf_phosphorus_g,
                                        .plant_nonstructural_carbon_g = canopy.plant_mobile_carbon_g[plant],
                                        .plant_nonstructural_nitrogen_g = canopy.plant_mobile_nitrogen_g[plant],
                                        .plant_nonstructural_phosphorus_g = canopy.plant_mobile_phosphorus_g[plant],
                                        .minimum_leaf_carbon_g = runscript.plant_pool_parameters.branch_structural_presence_g_per_plant * canopy.plant_population_count[plant],
                                        .minimum_daily_canopy_water_potential_megapascal = canopy.plant_minimum_daily_canopy_water_potential_megapascal[plant],
                                        .oxygen_stress_factor = phenology_root_oxygen_fraction[plant],
                                        .temperature_function = canopy.plant_uptake_growth_temperature_response[plant],
                                    });
                                    const values = try daily_plant_development_bank.row(plant);
                                    values[0] = @floatFromInt(@intFromEnum(output.phase));
                                    values[1] = @floatFromInt(output.branch_count);
                                    values[2] = output.main_branch_stage;
                                    values[3] = output.development_feedback;
                                    values[4] = output.leaf_nitrogen_to_carbon_ratio;
                                    values[5] = output.leaf_phosphorus_to_carbon_ratio;
                                    values[6] = output.minimum_daily_canopy_water_potential_megapascal;
                                    values[7] = output.oxygen_stress_factor;
                                    values[8] = output.temperature_function;
                                    // Only assigned species produce output; an unassigned population
                                    // would emit all-zero rows, so no file is created for it.
                                    const species_label = outputSpeciesLabel(plant_assignments, plant_unit_by_cell, cell, species) orelse continue;
                                    const file_name = try ecosys.output_record.buildOutputFileName(allocator, site_by_cell[cell].latitude_degrees_north, site_by_cell[cell].longitude_degrees_east, .{ .species = species_label }, current_year, development_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_plant_development_bank.streams[plant].write(file_name, daily_plant_development_catalog.variables, selection, resolved.plant_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = completed_day_of_year, .month = completed_management_date.month, .day = completed_management_date.day, .hour = 23 },
                                        .longitude_degrees_east = site_by_cell[cell].longitude_degrees_east,
                                        .latitude_degrees_north = site_by_cell[cell].latitude_degrees_north,
                                        .values = values,
                                    });
                                };
                            }
                        }
                    }
                    if (active_options.visualizationIncludesYear(current_year)) {
                        const canopy = if (detailed_canopy_state) |*value|
                            value
                        else
                            return error.VisualizationDailyOutputRequiresCanopyState;
                        const roots = if (plant_root_state) |*value|
                            value
                        else
                            return error.VisualizationDailyOutputRequiresRootState;
                        for (0..state.cell_count) |cell| {
                            var total_leaf_area_m2: f64 = 0;
                            for (0..config.plant_populations) |species| {
                                const plant =
                                    cell * config.plant_populations + species;
                                const branches = try canopy.branchRange(plant);
                                var leaf_carbon_g: f64 = 0;
                                var stalk_carbon_g: f64 = 0;
                                for (branches.first..branches.end) |branch| {
                                    leaf_carbon_g +=
                                        canopy.branch_leaf_carbon_g[branch];
                                    stalk_carbon_g +=
                                        canopy.branch_stalk_carbon_g[branch];
                                    total_leaf_area_m2 +=
                                        canopy.branch_leaf_area_m2[branch];
                                }
                                @memset(
                                    visualization_root_profile_carbon_g,
                                    0,
                                );
                                var root_carbon_g: f64 = 0;
                                for (0..config.soil_layers) |local_layer| {
                                    for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                                        const root = try roots.layerIndex(
                                            plant,
                                            domain,
                                            local_layer,
                                        );
                                        visualization_root_profile_carbon_g[local_layer] +=
                                            roots.mobile_carbon_g[root];
                                        for (0..roots.active_root_axis_count[plant]) |axis| {
                                            const axis_layer =
                                                try roots.layerAxisIndex(
                                                    plant,
                                                    domain,
                                                    local_layer,
                                                    axis,
                                                );
                                            visualization_root_profile_carbon_g[local_layer] +=
                                                roots.axis_primary_carbon_g[axis_layer] +
                                                roots.axis_secondary_carbon_g[axis_layer];
                                        }
                                    }
                                    root_carbon_g +=
                                        visualization_root_profile_carbon_g[local_layer];
                                }
                                visualization_leaf_carbon_g[species] =
                                    leaf_carbon_g;
                                visualization_stalk_carbon_g[species] =
                                    stalk_carbon_g;
                                visualization_root_carbon_g[species] =
                                    root_carbon_g;
                            }

                            const first_layer = try state.layerIndex(cell, 0);
                            const active_layers =
                                state.active_soil_layer_count[cell];
                            const carbon_pools =
                                try ecosys.soil_daily_carbon_pools.calculate(
                                    &soil_organic_state,
                                    &surface_organic_state,
                                    first_layer,
                                    active_layers,
                                    cell,
                                    visualization_organic_carbon_g[1 .. active_layers + 1],
                                );
                            visualization_organic_carbon_g[0] =
                                carbon_pools.surface_noncharcoal_organic_carbon_g_c;
                            @memset(
                                visualization_organic_carbon_g[active_layers + 1 ..],
                                0,
                            );
                            const active_boundary_base =
                                cell * (soil_geometry_state.layer_capacity + 1);
                            const active_surface_depth_m =
                                soil_geometry_state.boundary_depth_m[
                                    active_boundary_base +
                                        soil_geometry_state.first_active_layer[cell]
                                ];
                            const first_plant =
                                cell * config.plant_populations;
                            const plant_end =
                                first_plant + config.plant_populations;
                            _ = try visualization_streams.calculateAndWriteDaily(
                                cell,
                                site_by_cell[cell].longitude_degrees_east,
                                site_by_cell[cell].latitude_degrees_north,
                                .{
                                    .year = current_year,
                                    .day_of_year = completed_day_of_year,
                                    .month = completed_management_date.month,
                                    .day = completed_management_date.day,
                                    .hour = 23,
                                },
                                .{
                                    .cell_area_m2 = canopy_cell_area_m2[cell],
                                    .cumulative_transpiration_source_m3_by_species = plant_daily_flux_ledger.transpiration_source_m3[first_plant..plant_end],
                                    .cumulative_evapotranspiration_m3 = daily_water_ledger.evaporation_m3[cell],
                                    .cumulative_runoff_m3 = daily_water_ledger.runoff_m3[cell],
                                    .cumulative_outflow_m3 = daily_water_ledger.water_outflow_m3[cell] +
                                        daily_water_ledger.lateral_water_outflow_m3[cell],
                                    .water_table_bottom_depth_m = soil_boundary_topology_state.internal_water_table_depth_m[cell],
                                    .active_surface_depth_m = active_surface_depth_m,
                                    .snow_depth_m = snow_depth_m[cell],
                                    .total_leaf_area_m2 = total_leaf_area_m2,
                                    .leaf_carbon_g_by_species = visualization_leaf_carbon_g,
                                    .stalk_carbon_g_by_species = visualization_stalk_carbon_g,
                                    .root_carbon_g_by_species = visualization_root_carbon_g,
                                    .organic_carbon_g_by_layer_including_surface = visualization_organic_carbon_g,
                                },
                            );
                        }
                    }
                    // Legacy SOIL writes OUTPD after hour 24; DAY carries the
                    // balances and resets these fields before the next hour 1.
                    try plant_daily_flux_ledger.closeDayAfterOutput();
                    daily_soil_gas_flux.reset();
                    daily_canopy_gas_exchange.reset();
                    plant_water_publication_state.resetDaily();
                    daily_heterotrophic_respiration.resetDaily();
                    daily_carbon_export.resetDaily();
                    daily_phosphorus_export.resetDaily();
                    daily_nitrogen_export.resetDaily();
                    daily_water_ledger.resetDaily();
                    daily_heat_ledger.reset();
                    atmospheric_solute_input_ledger_state.resetDaily();
                    @memset(daily_soil_fire_charcoal_production_g_c, 0);
                    @memset(daily_soil_fire_carbon_dioxide_emission_g_c, 0);
                    @memset(daily_soil_fire_methane_emission_g_c, 0);
                    @memset(daily_soil_fire_phosphorus_flux_g_p, 0);
                    @memset(daily_soil_fire_nitrogen_flux_g_n, 0);
                    @memset(daily_root_fire_carbon_dioxide_emission_g_c, 0);
                    @memset(daily_root_fire_methane_emission_g_c, 0);
                    @memset(daily_root_fire_nitrogen_flux_g_n, 0);
                    @memset(daily_root_fire_phosphorus_flux_g_p, 0);
                    @memset(daily_root_fire_nitrogen_loss_g_n, 0);
                    @memset(daily_root_fire_phosphorus_loss_g_p, 0);
                    @memset(daily_soil_combusted_nitrogen_g_n, 0);
                    @memset(daily_soil_combusted_phosphorus_g_p, 0);
                    @memset(daily_microbial_phosphate_mineralization_g_p, 0);
                    @memset(daily_microbial_nitrogen_mineralization_g_n, 0);
                    @memset(daily_organic_fertilizer_carbon_input_g_c, 0);
                    @memset(daily_organic_fertilizer_phosphorus_input_g_p, 0);
                    @memset(daily_organic_fertilizer_nitrogen_input_g_n, 0);
                    soil_fertilizer_inventory.resetDaily();
                    mineral_fertilizer_inventory.resetDaily();
                    @memset(daily_biome_organic_carbon_input_g_c, 0);
                }
                if (scene_weather_hours < 24) {
                    const diagnostic_totals = try reconstructLandscapeMassBalance(hourly_science_context);
                    std.log.debug("phosphorus end hourly loop: hour={d} stored_g={e}", .{ scene_weather_hours + 1, diagnostic_totals.residue_phosphorus_g + diagnostic_totals.organic_phosphorus_g + diagnostic_totals.phosphate_phosphorus_g });
                }
                previous_weather_timestamp = timestamp;
                scene_weather_hours = try std.math.add(usize, scene_weather_hours, 1);
                executed_weather_hours = try std.math.add(usize, executed_weather_hours, 1);
                const scene_final_hour = scene_weather_hours == scene_expected_hours;
                if (try ecosys.checkpoint_schedule.shouldPublish(active_options.checkpoint_output_enabled, active_options.checkpoint_interval_days, current_day_of_year, timestamp.hour, scene_final_hour)) {
                    if (plant_root_state == null or detailed_canopy_state == null or canopy_precipitation_retention_state == null or plant_phenology_state == null or plant_growth_stage_state == null or plant_dormancy_state == null or branch_development_state == null) return error.CheckpointRequestedWithoutCompletePlantState;
                    try ecosys.checkpoint_bundle_writer.publish(
                        allocator,
                        init.io,
                        std.Io.Dir.cwd(),
                        @intCast(executed_weather_hours),
                        .{
                            .year = current_year,
                            .day_of_year = current_day_of_year,
                            .hour = if (timestamp.hour == 24) 23 else timestamp.hour,
                            .execution_iteration = @intCast(pass.execution_iteration),
                            .scenario_index = @intCast(pass.scenario_index),
                            .scenario_iteration = @intCast(pass.scenario_iteration),
                            .scene_index = @intCast(pass.scene_index),
                            .completed_scene_hours = @intCast(scene_weather_hours),
                        },
                        .{ .columns = config.lon_count, .rows = config.lat_count, .soil_layers = config.soil_layers, .snow_layers = snow_transport_state.layer_capacity, .plant_species_per_cell = config.plant_populations, .root_axes_per_plant = runscript.root_axes_per_plant },
                        .{
                            .grid = &state,
                            .plants = &plant_state,
                            .plant_metadata_cells = checkpoint_metadata_cells,
                            .plant_development = .{ .phenology = &plant_phenology_state.?, .growth = &plant_growth_stage_state.?, .dormancy = &plant_dormancy_state.?, .branch_development = &branch_development_state.? },
                            .plant_roots = &plant_root_state.?,
                            .plant_canopy = .{ .canopy = &detailed_canopy_state.?, .retention = &canopy_precipitation_retention_state.?, .layer_distribution = &canopy_layer_distribution_state.? },
                            .soil_biogeochemistry = .{ .microbial = &soil_microbial_state, .chemistry = &initial_chemistry_state, .available_nutrients = &plant_available_nutrient_state, .fertilizer = &soil_fertilizer_inventory, .mineral_fertilizer = &mineral_fertilizer_inventory, .fertilizer_band = &fertilizer_band_state, .reactive_nitrogen = &soil_reactive_nitrogen_state, .microbial_phosphorus = &soil_microbial_phosphorus_state },
                            .soil_organic_matter = .{ .profile = &soil_organic_state, .surface = &surface_organic_state, .litter_chemistry = &surface_litter_chemistry_state, .litter_fertilizer = &surface_litter_fertilizer_state, .surface_respiration = &surface_microbial_respiration_state, .surface_denitrification = &surface_denitrification_state, .surface_fire_exchange = &surface_fire_exchange_state, .litter_salt_ingress = &plant_litter_salt_ingress_state },
                            .transport = .{ .micropore = &micropore_solute_state, .macropore = &macropore_solute_state, .mineral_nitrogen = &mineral_nitrogen_transport_state, .organic = &soil_organic_transport_state, .gas = &gas_transport_state, .litter_gas = &litter_gas_transport_state, .snow = &snow_transport_state, .surface = &surface_transport_state },
                            .soil_geometry_and_hydrology = .{ .geometry = &soil_geometry_state, .hydrology = &transport_hydrology_state, .surface = &surface_precipitation_state, .erosion = &surface_erosion_state, .climate = &climate_state, .eroded_minerals = &eroded_mineral_state, .runtime = .{ .soil_properties = &soil_solver_property_state, .soil_thermal = &soil_thermal_state }, .surface_boundary = .{ .ground_air = &ground_air_state, .surface_aerodynamics = &surface_aerodynamic_state }, .surface_litter_geometry = &surface_litter_geometry_state, .surface_litter_ice_m3 = surface_litter_ice_m3, .delayed_live_canopy_combustion_heat_megajoules = delayed_live_canopy_combustion_heat_megajoules, .delayed_standing_dead_combustion_heat_megajoules = delayed_standing_dead_combustion_heat_megajoules, .delayed_subsurface_combustion_heat_megajoules = delayed_subsurface_combustion_heat_megajoules, .delayed_surface_combustion_heat_megajoules = delayed_surface_combustion_heat_megajoules },
                            .landscape_mass_balance = landscape_mass_balance_state,
                        },
                        .{ .section_write_buffer_bytes = 64 * 1024, .section_verify_buffer_bytes = 64 * 1024, .manifest_write_buffer_bytes = 16 * 1024 },
                    );
                }
            }
            if (scene_weather_hours == 0) return error.EmptyWeatherStream;
            executed_scene_passes = try std.math.add(usize, executed_scene_passes, 1);
        }
    }
    try hourly_soil_carbon_bank.finish();
    try hourly_plant_carbon_bank.finish();
    try hourly_soil_water_bank.finish();
    try hourly_plant_water_bank.finish();
    try hourly_soil_nitrogen_bank.finish();
    try hourly_plant_nitrogen_bank.finish();
    try hourly_soil_phosphorus_bank.finish();
    try hourly_plant_phosphorus_bank.finish();
    try hourly_soil_heat_bank.finish();
    try daily_soil_carbon_bank.finish();
    try daily_soil_water_bank.finish();
    try daily_soil_phosphorus_bank.finish();
    try daily_soil_nitrogen_bank.finish();
    try daily_soil_heat_bank.finish();
    try hourly_plant_heat_bank.finish();
    try daily_plant_carbon_bank.finish();
    try daily_plant_water_bank.finish();
    try daily_plant_nitrogen_bank.finish();
    try daily_plant_phosphorus_bank.finish();
    try daily_plant_development_bank.finish();
    try visualization_streams.finish();
    if (executed_weather_hours == 0) return error.EmptyWeatherStream;
    try atmospheric_state.validateFinite();
    if (canopy_optics_state) |optics_state| try optics_state.validateFinite();
    if (canopy_structure_state) |structure_state| try structure_state.validateFinite();
    if (canopy_interception_state) |interception_state| try interception_state.validateFinite();
    try ground_radiation_state.validateFinite();
    if (canopy_exposure_state) |exposure_state| try exposure_state.validateFinite();
    try surface_energy_state.validateFinite();
    try surface_temperature_solver_state.validateFinite();
    const surface_solver_diagnostics = try surface_temperature_solver_state.validateConvergence();
    std.log.info("surface Newton-Raphson/Picard diagnostics: cells={d}, iterations={d}, newton_raphson_steps={d}, picard_steps={d}, cells_using_picard={d}, maximum_absolute_residual_megajoules_per_m2={e}", .{ state.cell_count, surface_solver_diagnostics.total_iterations, surface_solver_diagnostics.total_newton_raphson_steps, surface_solver_diagnostics.total_picard_steps, surface_solver_diagnostics.cells_using_picard, surface_solver_diagnostics.maximum_absolute_residual_megajoules_per_m2 });
    if (canopy_energy_state) |energy_state| try energy_state.validateFinite();
    if (detailed_canopy_state) |canopy_state| try canopy_state.validateFinite();
    if (plant_root_state) |root_state| try root_state.validateFinite();
    if (branch_development_state) |development_state| try development_state.validateFinite();
    try terrain_radiation_state.validateFinite();
    std.log.info("ecosys-ng inputs validated: cells={d}, soil_layer_capacity={d}, soil_profiles={d}, plant_species_capacity={d}, plant_trait_profiles={d}, plant_management_schedules={d}, plant_assignment_units={d}, land_management_units={d}, disturbance_schedules={d}, fertilizer_schedules={d}, irrigation_schedules={d}, output_selections={d}, landscape_units={d}, scenarios={d}, scenes={d}, planned_days={d}, planned_hours={d}, total_weather_records={d}, executed_scene_passes={d}, executed_weather_hours={d}, latitude_deg_n={d:.3}, first_scene_water_heat_solute_max_iterations={d}, first_scene_gas_max_iterations={d}", .{
        state.cell_count,
        config.soil_layers,
        soil_catalog.entries.items.len,
        config.plant_populations,
        plant_catalog.entries.items.len,
        plant_management_catalog.entries.items.len,
        if (plant_assignments) |assignments| assignments.units.len else 0,
        if (land_management_assignments) |assignments| assignments.units.len else 0,
        disturbance_catalog.entries.items.len,
        fertilizer_catalog.entries.items.len,
        irrigation_catalog.entries.items.len,
        output_selection_catalog.entries.items.len,
        topography.units.len,
        runscript.scenarios.len,
        runscript.scenes.len,
        timeline.execution_weighted_days,
        timeline.execution_weighted_hours,
        total_weather_records,
        executed_scene_passes,
        executed_weather_hours,
        site.latitude_degrees_north,
        iteration_limits.water_heat_solute_max_iterations,
        iteration_limits.gas_max_iterations,
    });
}
