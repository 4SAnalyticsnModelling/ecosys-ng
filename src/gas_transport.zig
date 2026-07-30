const std = @import("std");

pub const Species = enum(u8) {
    carbon_dioxide,
    methane,
    oxygen,
    nitrogen,
    nitrous_oxide,
    ammonia,
    hydrogen,
};

pub const species_count = @typeInfo(Species).@"enum".fields.len;

pub fn massIndex(cell: usize, species: Species, cell_count: usize) !usize {
    if (cell >= cell_count) return error.GasTransportCellIndexOutOfBounds;
    return try std.math.add(usize, try std.math.mul(usize, cell, species_count), @intFromEnum(species));
}

/// ecosys stores each gas as the mass of its tracked element: C for CO2/CH4,
/// O2 as molecular oxygen, N for N2/N2O/NH3, and H2 as molecular hydrogen.
pub const g_per_mol_tracked = [species_count]f64{ 12, 12, 32, 28, 28, 14, 2 };
pub const atmospheric_boundary_multiplier = [species_count]f64{ 0.74, 1.04, 0.83, 0.86, 0.74, 1.02, 2.08 };

pub const SurfaceSolubilityParameters = struct {
    reference_water_to_air: [species_count]f64,
    log_intercept: [species_count]f64,
    temperature_coefficient_per_c: [species_count]f64,
};

/// HOUR1 L=0 temperature-dependent gas solubilities. Array order follows
/// `Species`, retaining the distinct NH3 coefficient even though STARTE does
/// not seed aqueous litter NH3.
pub fn surfaceSolubilityWaterToAir(temperature_k: f64, parameters: SurfaceSolubilityParameters) ![species_count]f64 {
    if (!std.math.isFinite(temperature_k) or temperature_k <= 0) return error.InvalidSurfaceGasTemperature;
    const temperature_c = temperature_k - 273.15;
    var result: [species_count]f64 = undefined;
    for (&result, parameters.reference_water_to_air, parameters.log_intercept, parameters.temperature_coefficient_per_c) |*value, reference, intercept, coefficient| {
        if (!std.math.isFinite(reference) or reference < 0 or !std.math.isFinite(intercept) or !std.math.isFinite(coefficient)) return error.InvalidSurfaceGasSolubilityParameter;
        value.* = reference * @exp(intercept - coefficient * temperature_c);
        if (!std.math.isFinite(value.*) or value.* < 0) return error.NonFiniteSurfaceGasSolubility;
    }
    return result;
}

/// STARTE L=0 gas initialization from runtime atmospheric concentrations.
/// Initial aqueous volume is passed explicitly because the source uses FC(0)
/// rather than assuming current liquid storage. NH3 aqueous mass starts at 0.
pub fn initializeSurfaceCell(state: *State, cell: usize, air_volume_m3: f64, initial_water_volume_m3: f64, temperature_k: f64, atmospheric_concentration_g_per_m3: [species_count]f64, solubility_parameters: SurfaceSolubilityParameters) !void {
    if (cell >= state.cell_count or !std.math.isFinite(air_volume_m3) or air_volume_m3 < 0 or !std.math.isFinite(initial_water_volume_m3) or initial_water_volume_m3 < 0) return error.InvalidSurfaceGasInitialization;
    const solubility = try surfaceSolubilityWaterToAir(temperature_k, solubility_parameters);
    var gaseous: [species_count]f64 = undefined;
    var dissolved: [species_count]f64 = undefined;
    for (&gaseous, &dissolved, atmospheric_concentration_g_per_m3, solubility, 0..) |*gas_mass, *water_mass, concentration, ratio, species_index| {
        if (!std.math.isFinite(concentration) or concentration < 0) return error.InvalidSurfaceGasInitialization;
        gas_mass.* = concentration * air_volume_m3;
        water_mass.* = if (species_index == @intFromEnum(Species.ammonia)) 0 else concentration * ratio * initial_water_volume_m3;
        if (!std.math.isFinite(gas_mass.*) or !std.math.isFinite(water_mass.*)) return error.NonFiniteSurfaceGasInitialization;
    }
    const first = cell * species_count;
    state.air_volume_m3[cell] = air_volume_m3;
    state.temperature_k[cell] = temperature_k;
    @memcpy(state.gaseous_mass_g[first .. first + species_count], &gaseous);
    @memcpy(state.dissolved_mass_g[first .. first + species_count], &dissolved);
    @memset(state.macropore_dissolved_mass_g[first .. first + species_count], 0);
    @memset(state.band_dissolved_mass_g[first .. first + species_count], 0);
}

pub const Face = struct {
    first_cell: usize,
    second_cell: usize,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    air_volume_m3: []f64,
    temperature_k: []f64,
    water_vapor_mol: []f64,
    /// Cell-major gas inventories in grams of the tracked element.
    gaseous_mass_g: []f64,
    dissolved_mass_g: []f64,
    /// Dissolved inventory carried by runtime soil macropore water.
    macropore_dissolved_mass_g: []f64,
    /// Band-water inventory; currently used by ammonia, retained as a full
    /// species-major buffer so kernels remain uniform and GPU-portable.
    band_dissolved_mass_g: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroGasTransportCellCount;
        const air = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(air);
        const temperature = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(temperature);
        const vapor = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(vapor);
        const n = try std.math.mul(usize, cell_count, species_count);
        const gaseous = try allocator.alloc(f64, n);
        errdefer allocator.free(gaseous);
        const dissolved = try allocator.alloc(f64, n);
        errdefer allocator.free(dissolved);
        const macropore_dissolved = try allocator.alloc(f64, n);
        errdefer allocator.free(macropore_dissolved);
        const band_dissolved = try allocator.alloc(f64, n);
        errdefer allocator.free(band_dissolved);
        @memset(air, 0);
        @memset(temperature, 0);
        @memset(vapor, 0);
        @memset(gaseous, 0);
        @memset(dissolved, 0);
        @memset(macropore_dissolved, 0);
        @memset(band_dissolved, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .air_volume_m3 = air, .temperature_k = temperature, .water_vapor_mol = vapor, .gaseous_mass_g = gaseous, .dissolved_mass_g = dissolved, .macropore_dissolved_mass_g = macropore_dissolved, .band_dissolved_mass_g = band_dissolved };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.band_dissolved_mass_g);
        self.allocator.free(self.macropore_dissolved_mass_g);
        self.allocator.free(self.dissolved_mass_g);
        self.allocator.free(self.gaseous_mass_g);
        self.allocator.free(self.water_vapor_mol);
        self.allocator.free(self.temperature_k);
        self.allocator.free(self.air_volume_m3);
        self.* = undefined;
    }

    pub fn clone(
        self: *const State,
        allocator: std.mem.Allocator,
    ) !State {
        var result = try State.init(allocator, self.cell_count);
        errdefer result.deinit();
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64)
                @memcpy(@field(&result, field.name), @field(self, field.name));
        }
        return result;
    }

    pub fn gaseousMasses(self: *State, cell: usize) ![]f64 {
        if (cell >= self.cell_count) return error.GasTransportCellIndexOutOfBounds;
        return self.gaseous_mass_g[cell * species_count .. (cell + 1) * species_count];
    }

    pub fn gaseousMassesConst(self: *const State, cell: usize) ![]const f64 {
        if (cell >= self.cell_count) return error.GasTransportCellIndexOutOfBounds;
        return self.gaseous_mass_g[cell * species_count .. (cell + 1) * species_count];
    }

    pub fn dissolvedMass(self: *State, cell: usize, species: Species) !*f64 {
        return &self.dissolved_mass_g[try massIndex(cell, species, self.cell_count)];
    }

    pub fn dissolvedMassConst(self: *const State, cell: usize, species: Species) !f64 {
        return self.dissolved_mass_g[try massIndex(cell, species, self.cell_count)];
    }

    /// Validates allocation dimensions before any field is indexed.
    pub fn validateShape(self: *const State) !void {
        const mass_count = std.math.mul(usize, self.cell_count, species_count) catch
            return error.CoupledGasStateSizeMismatch;
        if (self.cell_count == 0 or
            self.air_volume_m3.len != self.cell_count or
            self.temperature_k.len != self.cell_count or
            self.water_vapor_mol.len != self.cell_count or
            self.gaseous_mass_g.len != mass_count or
            self.dissolved_mass_g.len != mass_count or
            self.macropore_dissolved_mass_g.len != mass_count or
            self.band_dissolved_mass_g.len != mass_count)
        {
            return error.CoupledGasStateSizeMismatch;
        }
    }

    pub fn validateFinite(self: *const State) !void {
        try self.validateShape();
        inline for (@typeInfo(State).@"struct".fields) |declared| {
            if (declared.type == []f64) {
                for (@field(self, declared.name), 0..) |value, index| {
                    if (!std.math.isFinite(value)) {
                        std.log.warn(
                            "non-finite gas state: field={s} index={d} value={e}",
                            .{ declared.name, index, value },
                        );
                        return error.NonFiniteGasTransportState;
                    }
                    if (value < 0) {
                        std.log.warn(
                            "negative gas state: field={s} index={d} value={e}",
                            .{ declared.name, index, value },
                        );
                        return error.NegativeGasTransportState;
                    }
                }
            }
        }
        for (self.temperature_k) |temperature_k| {
            if (temperature_k <= 0) return error.InvalidGasTransportTemperature;
        }
    }
};

/// Penman linear-reduction expression used for litter and soil gas diffusion.
pub fn airFilledDiffusionGeometry(air_filled_porosity_m3_per_m3: f64, tortuosity: f64, total_porosity_m3_per_m3: f64, face_area_m2: f64, path_length_m: f64) !f64 {
    const v = [_]f64{ air_filled_porosity_m3_per_m3, tortuosity, total_porosity_m3_per_m3, face_area_m2, path_length_m };
    for (v) |x| if (!std.math.isFinite(x) or x < 0) return error.InvalidGasDiffusionGeometry;
    if (total_porosity_m3_per_m3 == 0 or path_length_m == 0) return error.InvalidGasDiffusionGeometry;
    return air_filled_porosity_m3_per_m3 * tortuosity * air_filled_porosity_m3_per_m3 / total_porosity_m3_per_m3 * face_area_m2 / path_length_m;
}

pub fn seriesConductance(interior_m3_per_step: f64, boundary_m3_per_step: f64) !f64 {
    if (!std.math.isFinite(interior_m3_per_step) or interior_m3_per_step < 0 or !std.math.isFinite(boundary_m3_per_step) or boundary_m3_per_step < 0) return error.InvalidGasConductance;
    const lower = @min(interior_m3_per_step, boundary_m3_per_step);
    const upper = @max(interior_m3_per_step, boundary_m3_per_step);
    return if (lower > 0) lower / (1 + lower / upper) else 0;
}

test "series conductance remains finite for an effectively open outer boundary" {
    const result = try seriesConductance(0.25, std.math.floatMax(f64));
    try std.testing.expectEqual(@as(f64, 0.25), result);
}

/// Conservative gaseous diffusion across one grid face. Positive flux moves
/// first -> second. Face coloring makes this kernel parallel without atomics.
pub fn calculateFaceDiffusiveFluxesG(state: *const State, face: Face, conductance_m3_per_step: []const f64, output_flux_g: []f64) !void {
    if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell) return error.InvalidGasTransportFace;
    if (conductance_m3_per_step.len != species_count or output_flux_g.len != species_count) return error.GasSpeciesCountMismatch;
    const first = try state.gaseousMassesConst(face.first_cell);
    const second = try state.gaseousMassesConst(face.second_cell);
    const first_air = state.air_volume_m3[face.first_cell];
    const second_air = state.air_volume_m3[face.second_cell];
    if (!std.math.isFinite(first_air) or first_air < 0 or !std.math.isFinite(second_air) or second_air < 0) return error.InvalidGasTransportState;
    for (first, second, conductance_m3_per_step, output_flux_g) |first_mass, second_mass, conductance, *flux| {
        if (!std.math.isFinite(first_mass) or first_mass < 0 or !std.math.isFinite(second_mass) or second_mass < 0 or !std.math.isFinite(conductance) or conductance < 0) return error.InvalidGasTransportState;
        const first_concentration = if (first_air > 0) first_mass / first_air else 0;
        const second_concentration = if (second_air > 0) second_mass / second_air else 0;
        flux.* = std.math.clamp(conductance * (first_concentration - second_concentration), -second_mass, first_mass);
        if (!std.math.isFinite(flux.*)) return error.NonFiniteGasTransportFlux;
    }
}

pub fn commitFaceFluxesG(state: *State, face: Face, flux_g: []const f64) !void {
    if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell) return error.InvalidGasTransportFace;
    if (flux_g.len != species_count) return error.GasSpeciesCountMismatch;
    const first = try state.gaseousMasses(face.first_cell);
    const second = try state.gaseousMasses(face.second_cell);
    for (first, second, flux_g) |a, b, flux| if (!std.math.isFinite(flux) or a - flux < -1e-12 or b + flux < -1e-12) return error.InsufficientGasForTransport;
    for (first, second, flux_g) |*a, *b, flux| {
        a.* = @max(0, a.* - flux);
        b.* = @max(0, b.* + flux);
    }
}

/// Exact bounded TRNSFRS atmospheric diffusion convention. Positive is into
/// the modeled cell. `iteration_fraction` is XNPG (= 1/NPG).
pub fn atmosphericDiffusiveFluxG(cell_mass_g: f64, cell_air_volume_m3: f64, atmospheric_concentration_g_per_m3: f64, conductance_m3_per_step: f64, iteration_fraction: f64) !f64 {
    const v = [_]f64{ cell_mass_g, cell_air_volume_m3, atmospheric_concentration_g_per_m3, conductance_m3_per_step, iteration_fraction };
    for (v) |x| if (!std.math.isFinite(x) or x < 0) return error.InvalidAtmosphericGasInput;
    if (iteration_fraction > 1) return error.InvalidAtmosphericGasInput;
    const requested = conductance_m3_per_step * (atmospheric_concentration_g_per_m3 - (if (cell_air_volume_m3 > 0) cell_mass_g / cell_air_volume_m3 else 0));
    return @min((atmospheric_concentration_g_per_m3 * cell_air_volume_m3 - cell_mass_g) * iteration_fraction, @max(-cell_mass_g * iteration_fraction, requested));
}

test "surface gas initialization reproduces STARTE and HOUR1 solubility" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    const reference = [species_count]f64{ 0.7391, 0.03156, 0.02925, 0.01510, 0.5241, 285.2, 0.03156 };
    const intercept = [species_count]f64{ 0.843, 0.597, 0.516, 0.456, 0.897, 0.513, 0.597 };
    const coefficient = [species_count]f64{ 0.0281, 0.0199, 0.0172, 0.0152, 0.0299, 0.0171, 0.0199 };
    const atmosphere = [species_count]f64{ 0.2, 0.001, 0.3, 0.8, 0.0001, 0.00001, 0.000001 };
    try initializeSurfaceCell(&state, 1, 2, 0.5, 293.15, atmosphere, .{ .reference_water_to_air = reference, .log_intercept = intercept, .temperature_coefficient_per_c = coefficient });
    const first = species_count;
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), state.gaseous_mass_g[first + @intFromEnum(Species.oxygen)], 1e-15);
    const oxygen_solubility = reference[@intFromEnum(Species.oxygen)] * @exp(0.516 - 0.0172 * 20);
    try std.testing.expectApproxEqAbs(atmosphere[@intFromEnum(Species.oxygen)] * oxygen_solubility * 0.5, state.dissolved_mass_g[first + @intFromEnum(Species.oxygen)], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), state.dissolved_mass_g[first + @intFromEnum(Species.ammonia)]);
}

/// Pressure/temperature correction from TRNSFRS. Fluxes are positive into the
/// cell and retain the existing gas mixture composition.
pub fn pressureDrivenFluxesG(air_volume_m3: f64, temperature_k: f64, water_vapor_mol: f64, gaseous_mass_g: []const f64, iteration_fraction: f64, output_flux_g: []f64) !void {
    if (gaseous_mass_g.len != species_count or output_flux_g.len != species_count) return error.GasSpeciesCountMismatch;
    if (!std.math.isFinite(air_volume_m3) or air_volume_m3 < 0 or !std.math.isFinite(temperature_k) or temperature_k <= 0 or !std.math.isFinite(water_vapor_mol) or water_vapor_mol < 0 or !std.math.isFinite(iteration_fraction) or iteration_fraction < 0 or iteration_fraction > 1) return error.InvalidGasPressureInput;
    var actual_mol = water_vapor_mol;
    for (gaseous_mass_g, g_per_mol_tracked) |mass, molar_mass| {
        if (!std.math.isFinite(mass) or mass < 0) return error.InvalidGasTransportState;
        actual_mol += mass / molar_mass;
    }
    if (actual_mol == 0) {
        @memset(output_flux_g, 0);
        return;
    }
    const capacity_mol = @max(0, 1.2194e4 * air_volume_m3 / temperature_k);
    const bulk_flux_mol = (capacity_mol - actual_mol) * iteration_fraction;
    for (output_flux_g, gaseous_mass_g, g_per_mol_tracked) |*flux, mass, molar_mass| {
        flux.* = bulk_flux_mol * (mass / molar_mass) / actual_mol * molar_mass;
        flux.* = @max(-mass * iteration_fraction, flux.*);
        if (!std.math.isFinite(flux.*)) return error.NonFiniteGasTransportFlux;
    }
}

/// TRNSFRS intercell pressure displacement. The receiving (`second`) cell's
/// temperature, air volume, vapor, and mixture define the bulk displacement;
/// positive flux is first -> second. Inventory clipping replaces the source's
/// dimensionally inconsistent mol-vs-g `AMIN1` bounds while retaining the
/// calculated scientific flux direction and mixture allocation.
pub fn adjacentPressureDrivenFluxesG(first_mass_g: []const f64, second_mass_g: []const f64, second_air_volume_m3: f64, second_temperature_k: f64, second_water_vapor_mol: f64, iteration_fraction: f64, output_flux_g: []f64) !void {
    if (first_mass_g.len != species_count or second_mass_g.len != species_count or output_flux_g.len != species_count) return error.GasSpeciesCountMismatch;
    if (!std.math.isFinite(second_air_volume_m3) or second_air_volume_m3 < 0 or !std.math.isFinite(second_temperature_k) or second_temperature_k <= 0 or !std.math.isFinite(second_water_vapor_mol) or second_water_vapor_mol < 0 or !std.math.isFinite(iteration_fraction) or iteration_fraction < 0 or iteration_fraction > 1) return error.InvalidGasPressureInput;
    var second_total_mol = second_water_vapor_mol;
    for (first_mass_g, second_mass_g, g_per_mol_tracked) |first, second, molar_mass| {
        if (!std.math.isFinite(first) or first < 0 or !std.math.isFinite(second) or second < 0) return error.InvalidGasTransportState;
        second_total_mol += second / molar_mass;
    }
    if (second_total_mol == 0) {
        @memset(output_flux_g, 0);
        return;
    }
    const capacity_mol = @max(0, 1.2194e4 * second_air_volume_m3 / second_temperature_k);
    const bulk_flux_mol = (capacity_mol - second_total_mol) * iteration_fraction;
    for (output_flux_g, first_mass_g, second_mass_g, g_per_mol_tracked) |*flux, first, second, molar_mass| {
        const requested_g = bulk_flux_mol * (second / molar_mass) / second_total_mol * molar_mass;
        // TRNSFR `AMIN1(V*G1,RFL*G)` uses the upstream molar inventory as
        // the numerical upper bound after RFL*G has been converted to grams.
        // Preserve that source behavior exactly; using the full gram mass
        // overstates convection by the tracked molar mass (12--32 times).
        flux.* = std.math.clamp(requested_g, -second, first / molar_mass);
        if (!std.math.isFinite(flux.*)) return error.NonFiniteGasTransportFlux;
    }
}

/// Air-water equilibration equation shared by all seven gases. Positive flux
/// dissolves gaseous mass into water; negative flux volatilizes aqueous mass.
pub fn phaseExchangeFluxG(gaseous_mass_g: f64, dissolved_mass_g: f64, air_volume_m3: f64, water_volume_m3: f64, mass_solubility_ratio: f64, exchange_rate_per_step: f64) !f64 {
    const v = [_]f64{ gaseous_mass_g, dissolved_mass_g, air_volume_m3, water_volume_m3, mass_solubility_ratio, exchange_rate_per_step };
    for (v) |x| if (!std.math.isFinite(x) or x < 0) return error.InvalidGasPhaseExchangeInput;
    const equivalent_water_volume = water_volume_m3 * mass_solubility_ratio;
    const total_equivalent_volume = equivalent_water_volume + air_volume_m3;
    if (total_equivalent_volume == 0) return 0;
    const flux = exchange_rate_per_step * (gaseous_mass_g * equivalent_water_volume - dissolved_mass_g * air_volume_m3) / total_equivalent_volume;
    if (!std.math.isFinite(flux)) return error.NonFiniteGasPhaseExchangeFlux;
    return std.math.clamp(flux, -dissolved_mass_g, gaseous_mass_g);
}

/// Supersaturation bubbling calculation from TRNSFRS. Dissolved masses are
/// expressed as tracked-element mass. Negative output means loss from water.
pub fn bubblingFluxesG(water_volume_m3: f64, temperature_k: f64, dissolved_mass_g: []const f64, mass_solubility_ratio: []const f64, iteration_fraction: f64, output_flux_g: []f64) !void {
    if (dissolved_mass_g.len != species_count or mass_solubility_ratio.len != species_count or output_flux_g.len != species_count) return error.GasSpeciesCountMismatch;
    if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or !std.math.isFinite(temperature_k) or temperature_k <= 0 or !std.math.isFinite(iteration_fraction) or iteration_fraction < 0 or iteration_fraction > 1) return error.InvalidGasBubblingInput;
    var equivalent_gas_mol: f64 = 0;
    for (dissolved_mass_g, mass_solubility_ratio, g_per_mol_tracked) |mass, solubility, molar_mass| {
        if (!std.math.isFinite(mass) or mass < 0 or !std.math.isFinite(solubility) or solubility <= 0) return error.InvalidGasBubblingInput;
        equivalent_gas_mol += mass / (molar_mass * solubility);
    }
    const capacity_mol = @max(0, 1.2194e4 * water_volume_m3 / temperature_k);
    if (equivalent_gas_mol <= capacity_mol or equivalent_gas_mol == 0) {
        @memset(output_flux_g, 0);
        return;
    }
    const bubble_volume_mol = (capacity_mol - equivalent_gas_mol) * iteration_fraction;
    for (output_flux_g, dissolved_mass_g, mass_solubility_ratio, g_per_mol_tracked) |*flux, mass, solubility, molar_mass| {
        const equivalent_species_mol = mass / (molar_mass * solubility);
        flux.* = @max(-mass * iteration_fraction, @min(0, bubble_volume_mol * equivalent_species_mol / equivalent_gas_mol * molar_mass * solubility));
        if (!std.math.isFinite(flux.*)) return error.NonFiniteGasBubblingFlux;
    }
}

test "gas diffusion geometry and series boundary preserve source equations" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.4 * 0.7 * 0.4 / 0.5 * 2 / 0.1), try airFilledDiffusionGeometry(0.4, 0.7, 0.5, 2, 0.1), 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), try seriesConductance(2, 3), 1e-14);
}

test "atmospheric and pressure gas fluxes are bounded" {
    const diffusion = try atmosphericDiffusiveFluxG(2, 1, 1, 10, 0.25);
    // MIN(equilibrium extent, MAX(inventory extent, gradient extent)).
    try std.testing.expectEqual(@as(f64, -0.5), diffusion);
    var flux: [species_count]f64 = undefined;
    try pressureDrivenFluxesG(0, 300, 0, &[_]f64{ 12, 0, 0, 0, 0, 0, 0 }, 0.25, &flux);
    try std.testing.expectEqual(@as(f64, -3), flux[0]);
}

test "phase exchange conserves gas plus dissolved inventory" {
    const flux = try phaseExchangeFluxG(4, 6, 2, 3, 0.5, 0.4);
    try std.testing.expectApproxEqAbs(@as(f64, -0.6857142857142857), flux, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 10), (4 + flux) + (6 - flux), 1e-14);
}

test "internal face diffusion conserves all seven gases" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.air_volume_m3[0] = 1;
    state.air_volume_m3[1] = 2;
    @memcpy(try state.gaseousMasses(0), &[_]f64{ 2, 3, 4, 5, 6, 7, 8 });
    var flux: [species_count]f64 = undefined;
    try calculateFaceDiffusiveFluxesG(&state, .{ .first_cell = 0, .second_cell = 1 }, &([_]f64{0.1} ** species_count), &flux);
    try commitFaceFluxesG(&state, .{ .first_cell = 0, .second_cell = 1 }, &flux);
    const first = try state.gaseousMassesConst(0);
    const second = try state.gaseousMassesConst(1);
    for (first, second, [_]f64{ 2, 3, 4, 5, 6, 7, 8 }) |a, b, total| try std.testing.expectApproxEqAbs(total, a + b, 1e-14);
}

test "bubbling only removes supersaturated dissolved gas" {
    var flux: [species_count]f64 = undefined;
    try bubblingFluxesG(0, 300, &[_]f64{ 12, 12, 32, 28, 28, 14, 2 }, &([_]f64{1} ** species_count), 0.25, &flux);
    for (flux, [_]f64{ 12, 12, 32, 28, 28, 14, 2 }) |value, mass| try std.testing.expectEqual(-0.25 * mass, value);
}

test "adjacent pressure displacement is conservative and donor bounded" {
    var flux: [species_count]f64 = undefined;
    try adjacentPressureDrivenFluxesG(&[_]f64{ 1, 0, 0, 0, 0, 0, 0 }, &[_]f64{ 12, 0, 0, 0, 0, 0, 0 }, 1, 300, 0, 1, &flux);
    try std.testing.expect(flux[0] <= 1);
    try std.testing.expect(flux[0] >= -12);
    for (flux) |value| try std.testing.expect(std.math.isFinite(value));
}
