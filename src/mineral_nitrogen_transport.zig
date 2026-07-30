const std = @import("std");
const transport = @import("solute_transport.zig");
const solver = @import("solute_transport_solver.zig");
const chemistry_module = @import("solute_chemistry_state.zig");
const reactive_module = @import("soil_reactive_nitrogen_state.zig");
const grid_module = @import("grid.zig");
const hydrology_module = @import("transport_hydrology.zig");
const geometry_module = @import("soil_face_geometry.zig");
const nutrient_parameters_module = @import("plant_root_nutrient_uptake.zig");

pub const Species = enum(u8) {
    ammonium_non_band,
    ammonium_band,
    ammonia_non_band,
    ammonia_band,
    nitrate_non_band,
    nitrate_band,
    nitrite_non_band,
    nitrite_band,
};

pub const species_count = @typeInfo(Species).@"enum".fields.len;

pub const ZoneFractions = struct {
    ammonium_non_band: f64,
    ammonium_band: f64,
    nitrate_non_band: f64,
    nitrate_band: f64,
};

/// Runtime-owned aqueous mineral-N inventories. Matrix and macropore amounts
/// remain distinct, matching the separate Z...S and Z...SH stores in TRNSFR.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    matrix: transport.State,
    macropore: transport.State,
    /// Positive values are nitrogen lost through an external boundary.
    boundary_export_g_n_per_step: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        var matrix = try transport.State.init(allocator, cell_count, species_count);
        errdefer matrix.deinit();
        var macropore = try transport.State.init(allocator, cell_count, species_count);
        errdefer macropore.deinit();
        const boundary_export = try allocator.alloc(f64, cell_count);
        @memset(boundary_export, 0);
        return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .matrix = matrix,
            .macropore = macropore,
            .boundary_export_g_n_per_step = boundary_export,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.boundary_export_g_n_per_step);
        self.macropore.deinit();
        self.matrix.deinit();
        self.* = undefined;
    }

    /// Imports authoritative reaction-state concentrations into the matrix
    /// domain. Macropore inventories are deliberately not overwritten: STARTE
    /// initializes them to zero and subsequent transport owns their history.
    pub fn initializeMatrix(
        self: *State,
        chemistry: *const chemistry_module.State,
        reactive: *const reactive_module.State,
        matrix_water_volume_m3: []const f64,
        fractions: ZoneFractions,
        nitrogen_molar_mass_g_per_mol: f64,
    ) !void {
        try validateBinding(self, chemistry, reactive, matrix_water_volume_m3, fractions, nitrogen_molar_mass_g_per_mol);
        @memcpy(self.matrix.water_volume_m3, matrix_water_volume_m3);
        for (0..self.cell_count) |cell| {
            const aqueous = chemistry.aqueous[cell];
            const water = matrix_water_volume_m3[cell];
            const amounts = try self.matrix.cellAmounts(cell);
            amounts[index(.ammonium_non_band)] = aqueous.ammonium_non_band * water * fractions.ammonium_non_band;
            amounts[index(.ammonium_band)] = aqueous.ammonium_band * water * fractions.ammonium_band;
            amounts[index(.ammonia_non_band)] = aqueous.ammonia_non_band * water * fractions.ammonium_non_band;
            amounts[index(.ammonia_band)] = aqueous.ammonia_band * water * fractions.ammonium_band;
            amounts[index(.nitrate_non_band)] = aqueous.nitrate_non_band * water * fractions.nitrate_non_band;
            amounts[index(.nitrate_band)] = aqueous.nitrate_band * water * fractions.nitrate_band;
            amounts[index(.nitrite_non_band)] = reactive.non_band_nitrite_g_n[cell] / nitrogen_molar_mass_g_per_mol;
            amounts[index(.nitrite_band)] = reactive.band_nitrite_g_n[cell] / nitrogen_molar_mass_g_per_mol;
        }
        try self.validate();
    }

    /// Captures the authoritative matrix inventory before hourly water
    /// movement. This is the TRNSFR `ZNH4S -> ZNH4S2` ownership boundary.
    pub fn captureHourStartMatrix(
        self: *State,
        chemistry: *const chemistry_module.State,
        reactive: *const reactive_module.State,
        matrix_water_volume_m3: []const f64,
        fractions: ZoneFractions,
        nitrogen_molar_mass_g_per_mol: f64,
    ) !void {
        try self.initializeMatrix(
            chemistry,
            reactive,
            matrix_water_volume_m3,
            fractions,
            nitrogen_molar_mass_g_per_mol,
        );
    }

    pub fn validate(self: *const State) !void {
        if (self.matrix.cell_count != self.cell_count or self.macropore.cell_count != self.cell_count or self.matrix.species_count != species_count or self.macropore.species_count != species_count) return error.InvalidMineralNitrogenTransportDimensions;
        for (self.matrix.amount_mol, self.macropore.amount_mol) |matrix_amount, macropore_amount| {
            if (!std.math.isFinite(matrix_amount) or matrix_amount < 0 or !std.math.isFinite(macropore_amount) or macropore_amount < 0) return error.InvalidMineralNitrogenTransportState;
        }
    }

    /// Publishes only the matrix domain back to reaction state. Macropore
    /// inventory stays transport-owned and cannot be silently collapsed.
    pub fn publishMatrix(
        self: *const State,
        chemistry: *chemistry_module.State,
        reactive: *reactive_module.State,
        fractions: ZoneFractions,
        nitrogen_molar_mass_g_per_mol: f64,
    ) !void {
        try validateBinding(self, chemistry, reactive, self.matrix.water_volume_m3, fractions, nitrogen_molar_mass_g_per_mol);
        for (0..self.cell_count) |cell| {
            const water = self.matrix.water_volume_m3[cell];
            const amounts = try self.matrix.cellAmountsConst(cell);
            chemistry.aqueous[cell].ammonium_non_band = try concentration(amounts[index(.ammonium_non_band)], water, fractions.ammonium_non_band);
            chemistry.aqueous[cell].ammonium_band = try concentration(amounts[index(.ammonium_band)], water, fractions.ammonium_band);
            chemistry.aqueous[cell].ammonia_non_band = try concentration(amounts[index(.ammonia_non_band)], water, fractions.ammonium_non_band);
            chemistry.aqueous[cell].ammonia_band = try concentration(amounts[index(.ammonia_band)], water, fractions.ammonium_band);
            chemistry.aqueous[cell].nitrate_non_band = try concentration(amounts[index(.nitrate_non_band)], water, fractions.nitrate_non_band);
            chemistry.aqueous[cell].nitrate_band = try concentration(amounts[index(.nitrate_band)], water, fractions.nitrate_band);
            reactive.non_band_nitrite_g_n[cell] = amounts[index(.nitrite_non_band)] * nitrogen_molar_mass_g_per_mol;
            reactive.band_nitrite_g_n[cell] = amounts[index(.nitrite_band)] * nitrogen_molar_mass_g_per_mol;
        }
    }
};

pub const FaceParameters = struct {
    allocator: std.mem.Allocator,
    matrix_conductance_m3_per_step: []f64,
    macropore_conductance_m3_per_step: []f64,
    mobility_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, face_count: usize) !FaceParameters {
        const count = try std.math.mul(usize, face_count, species_count);
        const matrix = try allocator.alloc(f64, count);
        errdefer allocator.free(matrix);
        const macropore = try allocator.alloc(f64, count);
        errdefer allocator.free(macropore);
        const mobility_values = try allocator.alloc(f64, count);
        @memset(matrix, 0);
        @memset(macropore, 0);
        @memset(mobility_values, 1);
        return .{ .allocator = allocator, .matrix_conductance_m3_per_step = matrix, .macropore_conductance_m3_per_step = macropore, .mobility_fraction = mobility_values };
    }

    pub fn deinit(self: *FaceParameters) void {
        self.allocator.free(self.mobility_fraction);
        self.allocator.free(self.macropore_conductance_m3_per_step);
        self.allocator.free(self.matrix_conductance_m3_per_step);
        self.* = undefined;
    }

    pub fn refresh(
        self: *FaceParameters,
        grid: *const grid_module.GridState,
        faces: *const hydrology_module.SoilFaces,
        geometry: *const geometry_module.State,
        matrix_bulk_volume_m3: []const f64,
        fractions: ZoneFractions,
        nutrient_parameters: nutrient_parameters_module.RuntimeParameters,
        step_h: f64,
    ) !void {
        if (matrix_bulk_volume_m3.len != grid.layer_count or faces.micropore_faces.len != geometry.face_area_m2.len or self.matrix_conductance_m3_per_step.len != faces.micropore_faces.len * species_count) return error.MineralNitrogenFaceDimensionMismatch;
        try validateFractions(fractions);
        if (!std.math.isFinite(step_h) or step_h <= 0) return error.InvalidMineralNitrogenFaceInput;
        for (faces.micropore_faces, 0..) |face, face_index| {
            const first = face.first_cell;
            const second = face.second_cell;
            const path_m = geometry.source_path_length_m[face_index] + geometry.destination_path_length_m[face_index];
            const area_m2 = geometry.face_area_m2[face_index];
            if (!std.math.isFinite(path_m) or path_m <= 0 or !std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidMineralNitrogenFaceGeometry;
            const first_macro_fraction = macroporeFraction(grid, first);
            const second_macro_fraction = macroporeFraction(grid, second);
            const first_theta = std.math.clamp(grid.matrix_liquid_water_m3[first] / matrix_bulk_volume_m3[first], 0, 1);
            const second_theta = std.math.clamp(grid.matrix_liquid_water_m3[second] / matrix_bulk_volume_m3[second], 0, 1);
            const matrix_tortuosity_per_m = (0.7 * first_theta * first_theta * (1 - first_macro_fraction) + 0.7 * second_theta * second_theta * (1 - second_macro_fraction)) / path_m;
            const first_macro_theta = if (grid.macropore_pore_capacity_m3[first] > 0) std.math.clamp(grid.macropore_liquid_water_m3[first] / grid.macropore_pore_capacity_m3[first], 0, 1) else 0;
            const second_macro_theta = if (grid.macropore_pore_capacity_m3[second] > 0) std.math.clamp(grid.macropore_liquid_water_m3[second] / grid.macropore_pore_capacity_m3[second], 0, 1) else 0;
            const macropore_tortuosity_per_m = (@min(1.0, 2.8 * first_macro_theta * first_macro_theta * first_macro_theta) * first_macro_fraction + @min(1.0, 2.8 * second_macro_theta * second_macro_theta * second_macro_theta) * second_macro_fraction) / path_m;
            const water_velocity_m_per_step = @abs(faces.micropore_water_flux_m3_per_step[face_index]) / area_m2;
            const mean_distance_m = 0.5 * path_m;
            const dispersion_m2_per_step = 0.20 * std.math.pow(f64, mean_distance_m, 1.07) * step_h * @min(step_h, water_velocity_m_per_step);
            for (0..species_count) |species_index| {
                const nutrient_index: usize = switch (@as(Species, @enumFromInt(species_index))) {
                    .ammonium_non_band, .ammonium_band, .ammonia_non_band, .ammonia_band => 0,
                    .nitrate_non_band, .nitrate_band, .nitrite_non_band, .nitrite_band => 1,
                };
                const diffusivity = try nutrient_parameters.diffusivityM2PerH(nutrient_index, grid.soil_temperature_k[second]) * step_h;
                if (!std.math.isFinite(diffusivity) or diffusivity < 0) return error.InvalidMineralNitrogenDiffusivity;
                const component = face_index * species_count + species_index;
                self.matrix_conductance_m3_per_step[component] = (diffusivity * matrix_tortuosity_per_m + dispersion_m2_per_step) * area_m2;
                self.macropore_conductance_m3_per_step[component] = diffusivity * macropore_tortuosity_per_m * area_m2;
                // Inventories are already extensive within their band/non-band
                // water share, so applying the zone fraction again would
                // incorrectly square its effect on convection and diffusion.
                self.mobility_fraction[component] = 1;
            }
        }
    }
};

pub const AdvanceInputs = struct {
    matrix_water_volume_m3: []const f64,
    macropore_water_volume_m3: []const f64,
    layer_volume_m3: []const f64,
    matrix_faces: []const transport.Face,
    macropore_faces: []const transport.Face,
    matrix_conductance_m3_per_step: []const f64,
    macropore_conductance_m3_per_step: []const f64,
    mobility_fraction: []const f64,
    matrix_external_water_flux_m3_per_step: []const f64,
    macropore_external_water_flux_m3_per_step: []const f64,
    maximum_convective_fraction: f64,
    pore_exchange_step_fraction: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    solver_options: solver.Options,
};

pub const AdvanceResult = struct {
    matrix: solver.Result,
    macropore: solver.Result,
};

/// Runs only converged mineral-N transport kernels. It never repeats the full
/// model cycle; NPH is solely the hybrid Newton/Picard iteration ceiling.
pub fn advance(allocator: std.mem.Allocator, state: *State, inputs: AdvanceInputs) !AdvanceResult {
    try validateAdvanceInputs(state, inputs);
    const matrix_before = try allocator.dupe(f64, state.matrix.amount_mol);
    defer allocator.free(matrix_before);
    const macropore_before = try allocator.dupe(f64, state.macropore.amount_mol);
    defer allocator.free(macropore_before);
    var committed = false;
    defer if (!committed) {
        @memcpy(state.matrix.amount_mol, matrix_before);
        @memcpy(state.macropore.amount_mol, macropore_before);
    };

    @memcpy(state.matrix.water_volume_m3, inputs.matrix_water_volume_m3);
    @memcpy(state.macropore.water_volume_m3, inputs.macropore_water_volume_m3);
    @memset(state.boundary_export_g_n_per_step, 0);
    const matrix_result = try solver.solve(allocator, &state.matrix, inputs.matrix_faces, inputs.matrix_conductance_m3_per_step, inputs.mobility_fraction, .{ .maximum_convective_fraction = inputs.maximum_convective_fraction }, inputs.solver_options);
    const macropore_result = try solver.solve(allocator, &state.macropore, inputs.macropore_faces, inputs.macropore_conductance_m3_per_step, inputs.mobility_fraction, .{ .maximum_convective_fraction = inputs.maximum_convective_fraction }, inputs.solver_options);

    for (0..state.cell_count) |cell| {
        try applyOutwardBoundary(state, &state.matrix, cell, inputs.matrix_external_water_flux_m3_per_step[cell], inputs.maximum_convective_fraction, inputs.nitrogen_molar_mass_g_per_mol);
        try applyOutwardBoundary(state, &state.macropore, cell, inputs.macropore_external_water_flux_m3_per_step[cell], inputs.maximum_convective_fraction, inputs.nitrogen_molar_mass_g_per_mol);
        for (0..species_count) |species_index| {
            const component = cell * species_count + species_index;
            const exchange = try transport.calculatePoreExchangeFlux(
                state.matrix.amount_mol[component],
                state.macropore.amount_mol[component],
                state.matrix.water_volume_m3[cell],
                state.macropore.water_volume_m3[cell],
                inputs.layer_volume_m3[cell],
                inputs.pore_exchange_step_fraction,
            );
            try transport.commitPoreExchange(&state.matrix.amount_mol[component], &state.macropore.amount_mol[component], exchange);
        }
    }
    try state.validate();
    committed = true;
    return .{ .matrix = matrix_result, .macropore = macropore_result };
}

fn applyOutwardBoundary(state: *State, domain: *transport.State, cell: usize, outward_water_m3: f64, maximum_fraction: f64, nitrogen_molar_mass: f64) !void {
    if (outward_water_m3 <= 0) return;
    const water = domain.water_volume_m3[cell];
    const fraction = if (water > 0) @min(maximum_fraction, outward_water_m3 / water) else 0;
    const amounts = try domain.cellAmounts(cell);
    for (amounts) |*amount| {
        const exported = amount.* * fraction;
        amount.* -= exported;
        state.boundary_export_g_n_per_step[cell] += exported * nitrogen_molar_mass;
    }
}

fn validateBinding(state: *const State, chemistry: *const chemistry_module.State, reactive: *const reactive_module.State, water: []const f64, fractions: ZoneFractions, molar_mass: f64) !void {
    if (chemistry.cell_count != state.cell_count or reactive.layer_count != state.cell_count or water.len != state.cell_count) return error.MineralNitrogenBindingDimensionMismatch;
    try validateFractions(fractions);
    if (!std.math.isFinite(molar_mass) or molar_mass <= 0) return error.InvalidMineralNitrogenBinding;
    for (water) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidMineralNitrogenBinding;
}

fn validateAdvanceInputs(state: *const State, inputs: AdvanceInputs) !void {
    try state.validate();
    const cells = state.cell_count;
    if (inputs.matrix_water_volume_m3.len != cells or inputs.macropore_water_volume_m3.len != cells or inputs.layer_volume_m3.len != cells or inputs.matrix_external_water_flux_m3_per_step.len != cells or inputs.macropore_external_water_flux_m3_per_step.len != cells) return error.MineralNitrogenTransportDimensionMismatch;
    if (inputs.matrix_faces.len != inputs.macropore_faces.len) return error.MineralNitrogenTransportDimensionMismatch;
    const components = try std.math.mul(usize, inputs.matrix_faces.len, species_count);
    if (inputs.matrix_conductance_m3_per_step.len != components or inputs.macropore_conductance_m3_per_step.len != components or inputs.mobility_fraction.len != components) return error.MineralNitrogenTransportDimensionMismatch;
    if (!std.math.isFinite(inputs.maximum_convective_fraction) or inputs.maximum_convective_fraction < 0 or inputs.maximum_convective_fraction > 1 or !std.math.isFinite(inputs.pore_exchange_step_fraction) or inputs.pore_exchange_step_fraction < 0 or inputs.pore_exchange_step_fraction > 1 or !std.math.isFinite(inputs.nitrogen_molar_mass_g_per_mol) or inputs.nitrogen_molar_mass_g_per_mol <= 0) return error.InvalidMineralNitrogenTransportInput;
}

fn index(species: Species) usize {
    return @intFromEnum(species);
}

fn concentration(amount_mol: f64, water_m3: f64, fraction: f64) !f64 {
    if (fraction == 0 or water_m3 == 0) {
        if (amount_mol > 1e-12) return error.MineralNitrogenInZeroWaterDomain;
        return 0;
    }
    const value = amount_mol / (water_m3 * fraction);
    if (!std.math.isFinite(value) or value < 0) return error.InvalidMineralNitrogenConcentration;
    return value;
}

fn validateFractions(fractions: ZoneFractions) !void {
    inline for (@typeInfo(ZoneFractions).@"struct".fields) |field| {
        const value = @field(fractions, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidMineralNitrogenZoneFraction;
    }
    if (@abs(fractions.ammonium_non_band + fractions.ammonium_band - 1) > 1e-12 or @abs(fractions.nitrate_non_band + fractions.nitrate_band - 1) > 1e-12) return error.InvalidMineralNitrogenZoneFraction;
}

fn macroporeFraction(grid: *const grid_module.GridState, cell: usize) f64 {
    const total = grid.matrix_pore_capacity_m3[cell] + grid.macropore_pore_capacity_m3[cell];
    return if (total > 0) std.math.clamp(grid.macropore_pore_capacity_m3[cell] / total, 0, 1) else 0;
}

test "matrix and macropore mineral nitrogen remain distinct and conservative" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.matrix.amount_mol[index(.nitrate_non_band)] = 2;
    const face = [_]transport.Face{.{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0 }};
    const conductance = [_]f64{0.1} ** species_count;
    const mobility = [_]f64{1} ** species_count;
    const result = try advance(std.testing.allocator, &state, .{
        .matrix_water_volume_m3 = &.{ 1, 1 },
        .macropore_water_volume_m3 = &.{ 0, 0 },
        .layer_volume_m3 = &.{ 1, 1 },
        .matrix_faces = &face,
        .macropore_faces = &face,
        .matrix_conductance_m3_per_step = &conductance,
        .macropore_conductance_m3_per_step = &conductance,
        .mobility_fraction = &mobility,
        .matrix_external_water_flux_m3_per_step = &.{ 0, 0 },
        .macropore_external_water_flux_m3_per_step = &.{ 0, 0 },
        .maximum_convective_fraction = 1,
        .pore_exchange_step_fraction = 1,
        .nitrogen_molar_mass_g_per_mol = 14,
        .solver_options = .{ .max_iterations = 40 },
    });
    try std.testing.expect(result.matrix.iterations < 40);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.matrix.amount_mol[index(.nitrate_non_band)] + state.matrix.amount_mol[species_count + index(.nitrate_non_band)], 1e-10);
    for (state.macropore.amount_mol) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "boundary export ledger is positive nitrogen loss" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.matrix.amount_mol[index(.ammonium_non_band)] = 2;
    const result = try advance(std.testing.allocator, &state, .{
        .matrix_water_volume_m3 = &.{1},
        .macropore_water_volume_m3 = &.{0},
        .layer_volume_m3 = &.{1},
        .matrix_faces = &.{},
        .macropore_faces = &.{},
        .matrix_conductance_m3_per_step = &.{},
        .macropore_conductance_m3_per_step = &.{},
        .mobility_fraction = &.{},
        .matrix_external_water_flux_m3_per_step = &.{0.25},
        .macropore_external_water_flux_m3_per_step = &.{0},
        .maximum_convective_fraction = 1,
        .pore_exchange_step_fraction = 1,
        .nitrogen_molar_mass_g_per_mol = 14,
        .solver_options = .{ .max_iterations = 4 },
    });
    _ = result;
    try std.testing.expectApproxEqAbs(@as(f64, 7), state.boundary_export_g_n_per_step[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), state.matrix.amount_mol[index(.ammonium_non_band)], 1e-12);
}

test "hour-start capture conserves zoned ammonium through wetting and drying" {
    const fractions: ZoneFractions = .{
        .ammonium_non_band = 0.75,
        .ammonium_band = 0.25,
        .nitrate_non_band = 0.6,
        .nitrate_band = 0.4,
    };
    for ([_]f64{ 2, 0.5 }) |water_after_m3| {
        var state = try State.init(std.testing.allocator, 1);
        defer state.deinit();
        var chemistry = try chemistry_module.State.init(std.testing.allocator, 1);
        defer chemistry.deinit();
        var reactive = try reactive_module.State.init(std.testing.allocator, 1, 1);
        defer reactive.deinit();
        chemistry.aqueous[0].ammonium_non_band = 2;
        chemistry.aqueous[0].ammonium_band = 4;
        state.macropore.amount_mol[index(.ammonium_non_band)] = 0.125;
        try state.captureHourStartMatrix(&chemistry, &reactive, &.{1}, fractions, 14);
        _ = try advance(std.testing.allocator, &state, .{
            .matrix_water_volume_m3 = &.{water_after_m3},
            .macropore_water_volume_m3 = &.{0.1},
            .layer_volume_m3 = &.{2},
            .matrix_faces = &.{},
            .macropore_faces = &.{},
            .matrix_conductance_m3_per_step = &.{},
            .macropore_conductance_m3_per_step = &.{},
            .mobility_fraction = &.{},
            .matrix_external_water_flux_m3_per_step = &.{0},
            .macropore_external_water_flux_m3_per_step = &.{0},
            .maximum_convective_fraction = 1,
            .pore_exchange_step_fraction = 0,
            .nitrogen_molar_mass_g_per_mol = 14,
            .solver_options = .{ .max_iterations = 4 },
        });
        try state.publishMatrix(&chemistry, &reactive, fractions, 14);
        try std.testing.expectApproxEqAbs(@as(f64, 1.5), state.matrix.amount_mol[index(.ammonium_non_band)], 1e-15);
        try std.testing.expectApproxEqAbs(@as(f64, 1), state.matrix.amount_mol[index(.ammonium_band)], 1e-15);
        try std.testing.expectApproxEqAbs(@as(f64, 0.125), state.macropore.amount_mol[index(.ammonium_non_band)], 1e-15);
        try std.testing.expectApproxEqAbs(2 / water_after_m3, chemistry.aqueous[0].ammonium_non_band, 1e-15);
        try std.testing.expectApproxEqAbs(4 / water_after_m3, chemistry.aqueous[0].ammonium_band, 1e-15);
    }
}
