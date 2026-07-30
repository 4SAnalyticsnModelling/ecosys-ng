const std = @import("std");

pub const Species = enum(u8) {
    carbon_dioxide_carbon,
    methane_carbon,
    oxygen,
    dinitrogen_nitrogen,
    nitrous_oxide_nitrogen,
    ammonium_nitrogen,
    ammonia_nitrogen,
    nitrate_nitrogen,
    hydrogen_phosphate_phosphorus,
    dihydrogen_phosphate_phosphorus,
    aluminum,
    iron,
    calcium,
    magnesium,
    sodium,
    potassium,
    sulfate_sulfur,
    chloride,
};

pub const species_count = @typeInfo(Species).@"enum".fields.len;
pub const nitrogen_g_per_mol: f64 = 14;
pub const phosphorus_g_per_mol: f64 = 31;

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    layer_capacity: usize,
    active: []bool,
    solid_snow_water_equivalent_m3: []f64,
    liquid_water_volume_m3: []f64,
    vapor_water_equivalent_m3: []f64,
    ice_volume_m3: []f64,
    air_filled_volume_m3: []f64,
    total_layer_volume_m3: []f64,
    target_layer_volume_m3: []f64,
    layer_thickness_m: []f64,
    cumulative_depth_m: []f64,
    snow_density_Mg_per_m3: []f64,
    temperature_k: []f64,
    heat_capacity_mj_per_k: []f64,
    horizontal_area_m2: []f64,
    amount_g: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, layer_capacity: usize) !State {
        if (cell_count == 0 or layer_capacity == 0) return error.ZeroSnowTransportDimension;
        const layer_count = try std.math.mul(usize, cell_count, layer_capacity);
        const active = try allocator.alloc(bool, layer_count);
        errdefer allocator.free(active);
        var physical: [13][]f64 = undefined;
        var physical_count: usize = 0;
        errdefer for (physical[0..physical_count]) |values| allocator.free(values);
        for (&physical) |*values| {
            values.* = try allocator.alloc(f64, layer_count);
            @memset(values.*, 0);
            physical_count += 1;
        }
        const amount = try allocator.alloc(f64, try std.math.mul(usize, layer_count, species_count));
        errdefer allocator.free(amount);
        @memset(active, false);
        @memset(amount, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .layer_capacity = layer_capacity, .active = active, .solid_snow_water_equivalent_m3 = physical[0], .liquid_water_volume_m3 = physical[1], .vapor_water_equivalent_m3 = physical[2], .ice_volume_m3 = physical[3], .air_filled_volume_m3 = physical[4], .total_layer_volume_m3 = physical[5], .target_layer_volume_m3 = physical[6], .layer_thickness_m = physical[7], .cumulative_depth_m = physical[8], .snow_density_Mg_per_m3 = physical[9], .temperature_k = physical[10], .heat_capacity_mj_per_k = physical[11], .horizontal_area_m2 = physical[12], .amount_g = amount };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.amount_g);
        inline for (.{ self.solid_snow_water_equivalent_m3, self.liquid_water_volume_m3, self.vapor_water_equivalent_m3, self.ice_volume_m3, self.air_filled_volume_m3, self.total_layer_volume_m3, self.target_layer_volume_m3, self.layer_thickness_m, self.cumulative_depth_m, self.snow_density_Mg_per_m3, self.temperature_k, self.heat_capacity_mj_per_k, self.horizontal_area_m2 }) |values| self.allocator.free(values);
        self.allocator.free(self.active);
        self.* = undefined;
    }

    /// STARTS snow-layer initialization with runtime layer boundaries.
    pub fn initializePhysicalState(self: *State, snow_depth_m: []const f64, cell_area_m2: []const f64, atmospheric_temperature_k: []const f64, layer_bottom_depth_m: []const f64, initial_snow_density_Mg_per_m3: f64) !void {
        if (snow_depth_m.len != self.cell_count or cell_area_m2.len != self.cell_count or atmospheric_temperature_k.len != self.cell_count or layer_bottom_depth_m.len != self.layer_capacity or !std.math.isFinite(initial_snow_density_Mg_per_m3) or initial_snow_density_Mg_per_m3 <= 0) return error.InvalidSnowPhysicalInitialization;
        var previous_bottom: f64 = 0;
        for (layer_bottom_depth_m) |bottom| {
            if (!std.math.isFinite(bottom) or bottom <= previous_bottom) return error.InvalidSnowLayerBoundary;
            previous_bottom = bottom;
        }
        for (0..self.cell_count) |cell| {
            const depth = snow_depth_m[cell];
            const area = cell_area_m2[cell];
            const air_temperature = atmospheric_temperature_k[cell];
            if (!std.math.isFinite(depth) or depth < 0 or !std.math.isFinite(area) or area <= 0 or !std.math.isFinite(air_temperature) or air_temperature <= 0) return error.InvalidSnowPhysicalInitialization;
            var cumulative: f64 = 0;
            for (0..self.layer_capacity) |layer| {
                const index = cell * self.layer_capacity + layer;
                const top = if (layer == 0) 0 else layer_bottom_depth_m[layer - 1];
                const nominal = layer_bottom_depth_m[layer] - top;
                const thickness = @min(nominal, @max(0, depth - top));
                const solid = thickness * initial_snow_density_Mg_per_m3 * area;
                const total = if (solid > 0) solid / initial_snow_density_Mg_per_m3 else 0;
                self.active[index] = thickness > 0;
                self.solid_snow_water_equivalent_m3[index] = solid;
                self.liquid_water_volume_m3[index] = 0;
                self.vapor_water_equivalent_m3[index] = 0;
                self.ice_volume_m3[index] = 0;
                self.total_layer_volume_m3[index] = total;
                self.target_layer_volume_m3[index] = nominal * area;
                self.air_filled_volume_m3[index] = @max(0, total - solid);
                self.layer_thickness_m[index] = thickness;
                cumulative += thickness;
                self.cumulative_depth_m[index] = cumulative;
                self.snow_density_Mg_per_m3[index] = initial_snow_density_Mg_per_m3;
                self.temperature_k[index] = @min(273.15, air_temperature);
                self.heat_capacity_mj_per_k[index] = 2.095 * solid;
                self.horizontal_area_m2[index] = area;
            }
        }
    }

    /// REDIST top-layer `TQS/TQW/THQS` update, committed only after every
    /// runtime cell has valid mass and energy inputs.
    pub fn commitAtmosphericWater(self: *State, solid_snow_input_m3: []const f64, liquid_water_input_m3: []const f64, heat_input_mj: []const f64, fallback_temperature_k: []const f64, initial_snow_density_Mg_per_m3: f64) !void {
        inline for (.{ solid_snow_input_m3.len, liquid_water_input_m3.len, heat_input_mj.len, fallback_temperature_k.len }) |length| if (length != self.cell_count) return error.SnowAtmosphericInputDimensionMismatch;
        for (0..self.cell_count) |cell| {
            inline for (.{ solid_snow_input_m3[cell], liquid_water_input_m3[cell], heat_input_mj[cell], fallback_temperature_k[cell] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowAtmosphericInput;
            if (solid_snow_input_m3[cell] < 0 or liquid_water_input_m3[cell] < 0 or fallback_temperature_k[cell] <= 0) return error.InvalidSnowAtmosphericInput;
            const top = cell * self.layer_capacity;
            const solid = self.solid_snow_water_equivalent_m3[top] + solid_snow_input_m3[cell];
            const liquid = self.liquid_water_volume_m3[top] + liquid_water_input_m3[cell];
            const heat_capacity = 2.095 * solid + 4.19 * (liquid + self.vapor_water_equivalent_m3[top]) + 1.9274 * self.ice_volume_m3[top];
            const old_energy = self.heat_capacity_mj_per_k[top] * self.temperature_k[top];
            const temperature = if (heat_capacity > 0) (old_energy + heat_input_mj[cell]) / heat_capacity else fallback_temperature_k[cell];
            if (!std.math.isFinite(temperature) or temperature <= 0) return error.InvalidSnowAtmosphericEnergy;
        }
        for (0..self.cell_count) |cell| {
            const top = cell * self.layer_capacity;
            const old_energy = self.heat_capacity_mj_per_k[top] * self.temperature_k[top];
            self.solid_snow_water_equivalent_m3[top] += solid_snow_input_m3[cell];
            self.liquid_water_volume_m3[top] += liquid_water_input_m3[cell];
            self.heat_capacity_mj_per_k[top] = 2.095 * self.solid_snow_water_equivalent_m3[top] + 4.19 * (self.liquid_water_volume_m3[top] + self.vapor_water_equivalent_m3[top]) + 1.9274 * self.ice_volume_m3[top];
            if (self.heat_capacity_mj_per_k[top] > 0) self.temperature_k[top] = (old_energy + heat_input_mj[cell]) / self.heat_capacity_mj_per_k[top];
            if (self.solid_snow_water_equivalent_m3[top] + self.liquid_water_volume_m3[top] + self.ice_volume_m3[top] > 0) {
                self.active[top] = true;
                if (self.snow_density_Mg_per_m3[top] <= 0) self.snow_density_Mg_per_m3[top] = initial_snow_density_Mg_per_m3;
            }
            self.refreshCellGeometry(cell);
        }
    }

    pub fn commitMeltWater(self: *State, downward_water_flux_m3: []const f64, litter_water_flux_m3: []const f64, soil_micropore_water_flux_m3: []const f64, soil_macropore_water_flux_m3: []const f64) !void {
        const layers = self.cell_count * self.layer_capacity;
        if (downward_water_flux_m3.len != layers or litter_water_flux_m3.len != self.cell_count or soil_micropore_water_flux_m3.len != self.cell_count or soil_macropore_water_flux_m3.len != self.cell_count) return error.SnowMeltCommitDimensionMismatch;
        const candidate = try self.allocator.dupe(f64, self.liquid_water_volume_m3);
        defer self.allocator.free(candidate);
        for (0..self.cell_count) |cell| {
            var bottom: ?usize = null;
            for (0..self.layer_capacity) |layer| {
                if (self.active[cell * self.layer_capacity + layer]) bottom = layer;
            }
            for (1..self.layer_capacity) |layer| {
                const destination = cell * self.layer_capacity + layer;
                const flux = downward_water_flux_m3[destination];
                if (!std.math.isFinite(flux) or flux < 0) return error.InvalidSnowMeltCommit;
                candidate[destination - 1] -= flux;
                candidate[destination] += flux;
            }
            if (bottom) |layer| {
                const discharge = litter_water_flux_m3[cell] + soil_micropore_water_flux_m3[cell] + soil_macropore_water_flux_m3[cell];
                if (!std.math.isFinite(discharge) or discharge < 0) return error.InvalidSnowMeltCommit;
                candidate[cell * self.layer_capacity + layer] -= discharge;
            }
        }
        for (candidate) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.SnowMeltExceedsLiquidInventory;
        for (candidate, self.liquid_water_volume_m3) |value, *destination| destination.* = @max(0, value);
        for (0..self.cell_count) |cell| self.refreshCellGeometry(cell);
    }

    fn refreshCellGeometry(self: *State, cell: usize) void {
        var cumulative: f64 = 0;
        for (0..self.layer_capacity) |layer| {
            const index = cell * self.layer_capacity + layer;
            const density = self.snow_density_Mg_per_m3[index];
            const solid_volume = if (density > 0) self.solid_snow_water_equivalent_m3[index] / density else 0;
            const total = solid_volume + self.liquid_water_volume_m3[index] + self.ice_volume_m3[index];
            self.total_layer_volume_m3[index] = total;
            self.air_filled_volume_m3[index] = @max(0, total - self.solid_snow_water_equivalent_m3[index] - self.liquid_water_volume_m3[index] - self.ice_volume_m3[index]);
            self.layer_thickness_m[index] = if (self.horizontal_area_m2[index] > 0) total / self.horizontal_area_m2[index] else 0;
            cumulative += self.layer_thickness_m[index];
            self.cumulative_depth_m[index] = cumulative;
        }
    }

    pub fn refreshAllGeometry(self: *State) void {
        for (0..self.cell_count) |cell| self.refreshCellGeometry(cell);
    }

    pub fn layerIndex(self: *const State, cell: usize, layer: usize) !usize {
        if (cell >= self.cell_count or layer >= self.layer_capacity) return error.SnowTransportIndexOutOfBounds;
        return cell * self.layer_capacity + layer;
    }

    pub fn amounts(self: *State, cell: usize, layer: usize) ![]f64 {
        const index = try self.layerIndex(cell, layer);
        return self.amount_g[index * species_count .. (index + 1) * species_count];
    }

    pub fn amountsConst(self: *const State, cell: usize, layer: usize) ![]const f64 {
        const index = try self.layerIndex(cell, layer);
        return self.amount_g[index * species_count .. (index + 1) * species_count];
    }
};

pub const SurfacePartition = struct {
    litter_cover_fraction: f64,
    bare_soil_fraction: f64,
    nonband_ammonium_fraction: f64,
    band_ammonium_fraction: f64,
    nonband_nitrate_fraction: f64,
    band_nitrate_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
};

pub const SurfaceDischarge = struct {
    litter_g: [species_count]f64 = [_]f64{0} ** species_count,
    soil_nonband_g: [species_count]f64 = [_]f64{0} ** species_count,
    soil_band_g: [species_count]f64 = [_]f64{0} ** species_count,
};

pub const Fluxes = struct {
    allocator: std.mem.Allocator,
    /// Layer-major flux from a layer to the next layer; zero at discharge layer.
    downward_g: []f64,
    surface_discharge: []SurfaceDischarge,

    pub fn deinit(self: *Fluxes) void {
        self.allocator.free(self.surface_discharge);
        self.allocator.free(self.downward_g);
        self.* = undefined;
    }
};

/// Converts precipitation/irrigation concentrations into tracked snow
/// inventories. The first five concentrations are already tracked-mass g/m3;
/// nutrient inputs are mol/m3 and receive the exact 14 or 31 multipliers.
pub fn atmosphericInputG(rain_water_m3: f64, irrigation_water_m3: f64, rain_first_five_g_per_m3: [5]f64, irrigation_first_five_g_per_m3: [5]f64, rain_nutrients_mol_per_m3: [5]f64, irrigation_nutrients_mol_per_m3: [5]f64, rain_ions_g_per_m3: [8]f64, irrigation_ions_g_per_m3: [8]f64) ![species_count]f64 {
    if (!std.math.isFinite(rain_water_m3) or rain_water_m3 < 0 or !std.math.isFinite(irrigation_water_m3) or irrigation_water_m3 < 0) return error.InvalidSnowAtmosphericInput;
    var output = [_]f64{0} ** species_count;
    for (rain_first_five_g_per_m3, irrigation_first_five_g_per_m3, 0..) |rain, irrigation, species| {
        if (!std.math.isFinite(rain) or rain < 0 or !std.math.isFinite(irrigation) or irrigation < 0) return error.InvalidSnowAtmosphericInput;
        output[species] = rain_water_m3 * rain + irrigation_water_m3 * irrigation;
    }
    for (rain_nutrients_mol_per_m3, irrigation_nutrients_mol_per_m3, 0..) |rain, irrigation, nutrient| {
        if (!std.math.isFinite(rain) or rain < 0 or !std.math.isFinite(irrigation) or irrigation < 0) return error.InvalidSnowAtmosphericInput;
        const factor = if (nutrient < 3) nitrogen_g_per_mol else phosphorus_g_per_mol;
        output[5 + nutrient] = (rain_water_m3 * rain + irrigation_water_m3 * irrigation) * factor;
    }
    for (rain_ions_g_per_m3, irrigation_ions_g_per_m3, 0..) |rain, irrigation, ion| {
        if (!std.math.isFinite(rain) or rain < 0 or !std.math.isFinite(irrigation) or irrigation < 0) return error.InvalidSnowAtmosphericInput;
        output[10 + ion] = rain_water_m3 * rain + irrigation_water_m3 * irrigation;
    }
    return output;
}

/// Evaluates one snow-water transport residual from a supplied trial state.
/// `water_flux_to_lower_m3` is layer-major. Surface fluxes apply only to the
/// first active layer whose lower neighbor is absent, matching `ICHKL`.
pub fn calculateFluxes(allocator: std.mem.Allocator, state: *const State, water_flux_to_lower_m3: []const f64, litter_water_flux_m3: []const f64, soil_micropore_water_flux_m3: []const f64, soil_macropore_water_flux_m3: []const f64, partitions: []const SurfacePartition) !Fluxes {
    const layer_count = try std.math.mul(usize, state.cell_count, state.layer_capacity);
    if (water_flux_to_lower_m3.len != layer_count or litter_water_flux_m3.len != state.cell_count or soil_micropore_water_flux_m3.len != state.cell_count or soil_macropore_water_flux_m3.len != state.cell_count or partitions.len != state.cell_count) return error.SnowTransportInputSizeMismatch;
    const downward = try allocator.alloc(f64, state.amount_g.len);
    errdefer allocator.free(downward);
    @memset(downward, 0);
    const discharge = try allocator.alloc(SurfaceDischarge, state.cell_count);
    errdefer allocator.free(discharge);
    @memset(discharge, .{});
    try validateState(state);
    for (0..state.cell_count) |cell| {
        try validatePartition(partitions[cell]);
        var discharged = false;
        for (0..state.layer_capacity) |layer| {
            const layer_index = try state.layerIndex(cell, layer);
            if (!state.active[layer_index]) continue;
            const amounts_g = try state.amountsConst(cell, layer);
            const has_active_lower = layer + 1 < state.layer_capacity and state.active[try state.layerIndex(cell, layer + 1)];
            if (has_active_lower) {
                const fraction = if (state.liquid_water_volume_m3[layer_index] > 0) std.math.clamp(water_flux_to_lower_m3[try state.layerIndex(cell, layer + 1)] / state.liquid_water_volume_m3[layer_index], 0, 1) else 1;
                for (amounts_g, 0..) |amount, species| downward[layer_index * species_count + species] = amount * fraction;
            } else if (!discharged) {
                const water = state.liquid_water_volume_m3[layer_index];
                const litter_fraction = if (water > 0) std.math.clamp(litter_water_flux_m3[cell] / water, 0, 1) else partitions[cell].litter_cover_fraction;
                const soil_fraction = if (water > 0) std.math.clamp((soil_micropore_water_flux_m3[cell] + soil_macropore_water_flux_m3[cell]) / water, 0, 1) else partitions[cell].bare_soil_fraction;
                if (litter_fraction + soil_fraction > 1 + 1e-12) return error.SnowSurfaceFluxExceedsInventory;
                routeSurface(amounts_g, litter_fraction, soil_fraction, partitions[cell], &discharge[cell]);
                discharged = true;
            }
        }
    }
    return .{ .allocator = allocator, .downward_g = downward, .surface_discharge = discharge };
}

pub fn commit(state: *State, atmospheric_top_input_g: []const f64, fluxes: *const Fluxes) !void {
    if (atmospheric_top_input_g.len != try std.math.mul(usize, state.cell_count, species_count) or fluxes.downward_g.len != state.amount_g.len or fluxes.surface_discharge.len != state.cell_count) return error.SnowTransportInputSizeMismatch;
    const candidate = try state.allocator.dupe(f64, state.amount_g);
    defer state.allocator.free(candidate);
    for (0..state.cell_count) |cell| {
        const top = try state.layerIndex(cell, 0);
        for (0..species_count) |species| candidate[top * species_count + species] += atmospheric_top_input_g[cell * species_count + species];
        for (0..state.layer_capacity) |layer| {
            const index = try state.layerIndex(cell, layer);
            if (!state.active[index]) continue;
            for (0..species_count) |species| {
                const component = index * species_count + species;
                const outgoing = fluxes.downward_g[component];
                candidate[component] -= outgoing;
                if (layer + 1 < state.layer_capacity and state.active[try state.layerIndex(cell, layer + 1)]) candidate[(index + 1) * species_count + species] += outgoing;
                if (layer + 1 >= state.layer_capacity or !state.active[try state.layerIndex(cell, layer + 1)]) candidate[component] -= fluxes.surface_discharge[cell].litter_g[species] + fluxes.surface_discharge[cell].soil_nonband_g[species] + fluxes.surface_discharge[cell].soil_band_g[species];
            }
        }
    }
    for (candidate) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.InvalidSnowTransportCandidate;
    for (candidate, state.amount_g) |value, *amount| amount.* = @max(0, value);
}

fn routeSurface(amounts_g: []const f64, litter_fraction: f64, soil_fraction: f64, partition: SurfacePartition, output: *SurfaceDischarge) void {
    for (amounts_g, 0..) |amount, species| {
        output.litter_g[species] = amount * litter_fraction;
        switch (@as(Species, @enumFromInt(species))) {
            .ammonium_nitrogen, .ammonia_nitrogen => {
                output.soil_nonband_g[species] = amount * soil_fraction * partition.nonband_ammonium_fraction;
                output.soil_band_g[species] = amount * soil_fraction * partition.band_ammonium_fraction;
            },
            .nitrate_nitrogen => {
                output.soil_nonband_g[species] = amount * soil_fraction * partition.nonband_nitrate_fraction;
                output.soil_band_g[species] = amount * soil_fraction * partition.band_nitrate_fraction;
            },
            .hydrogen_phosphate_phosphorus, .dihydrogen_phosphate_phosphorus => {
                output.soil_nonband_g[species] = amount * soil_fraction * partition.nonband_phosphate_fraction;
                output.soil_band_g[species] = amount * soil_fraction * partition.band_phosphate_fraction;
            },
            else => output.soil_nonband_g[species] = amount * soil_fraction,
        }
    }
}

fn validateState(state: *const State) !void {
    inline for (.{ state.solid_snow_water_equivalent_m3, state.liquid_water_volume_m3, state.vapor_water_equivalent_m3, state.ice_volume_m3, state.air_filled_volume_m3, state.total_layer_volume_m3, state.target_layer_volume_m3, state.layer_thickness_m, state.cumulative_depth_m, state.snow_density_Mg_per_m3, state.temperature_k, state.heat_capacity_mj_per_k, state.horizontal_area_m2 }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowTransportState;
    for (state.amount_g) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowTransportState;
}

fn validatePartition(partition: SurfacePartition) !void {
    inline for (@typeInfo(SurfacePartition).@"struct".fields) |field| {
        const value = @field(partition, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSnowSurfacePartition;
    }
    if (@abs(partition.litter_cover_fraction + partition.bare_soil_fraction - 1) > 1e-10 or @abs(partition.nonband_ammonium_fraction + partition.band_ammonium_fraction - 1) > 1e-10 or @abs(partition.nonband_nitrate_fraction + partition.band_nitrate_fraction - 1) > 1e-10 or @abs(partition.nonband_phosphate_fraction + partition.band_phosphate_fraction - 1) > 1e-10) return error.InvalidSnowSurfacePartition;
}

test "snow atmospheric nutrients and irrigation ions retain units" {
    const input = try atmosphericInputG(2, 1, [_]f64{ 1, 2, 3, 4, 5 }, [_]f64{ 2, 3, 4, 5, 6 }, [_]f64{1} ** 5, [_]f64{2} ** 5, [_]f64{0} ** 8, [_]f64{1} ** 8);
    try std.testing.expectEqual(@as(f64, 4), input[0]);
    try std.testing.expectEqual(@as(f64, 56), input[5]);
    try std.testing.expectEqual(@as(f64, 124), input[8]);
    try std.testing.expectEqual(@as(f64, 1), input[@intFromEnum(Species.aluminum)]);
}

test "STARTS snow physical state uses runtime layer boundaries" {
    var state = try State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    try state.initializePhysicalState(&.{ 0.2, 0 }, &.{ 10, 20 }, &.{ 278.0, 268.0 }, &.{ 0.05, 0.125, 0.25 }, 0.05);
    try std.testing.expectEqualSlices(bool, &.{ true, true, true, false, false, false }, state.active);
    try std.testing.expectApproxEqAbs(@as(f64, 0.025), state.solid_snow_water_equivalent_m3[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0375), state.solid_snow_water_equivalent_m3[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0375), state.solid_snow_water_equivalent_m3[2], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.cumulative_depth_m[2], 1e-12);
    try std.testing.expectEqual(@as(f64, 273.15), state.temperature_k[0]);
    try std.testing.expectEqual(@as(f64, 268), state.temperature_k[3]);
}

test "REDIST atmospheric input and WATSUB melt commit conserve snow water" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try state.initializePhysicalState(&.{0}, &.{10}, &.{270}, &.{ 0.05, 0.125 }, 0.05);
    try state.commitAtmosphericWater(&.{0.2}, &.{0.1}, &.{2.095 * 0.2 * 268 + 4.19 * 0.1 * 268}, &.{268}, 0.05);
    try std.testing.expect(state.active[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 268), state.temperature_k[0], 1e-10);
    const before = state.liquid_water_volume_m3[0] + state.liquid_water_volume_m3[1];
    try state.commitMeltWater(&.{ 0, 0.04 }, &.{0.02}, &.{0.01}, &.{0.01});
    const after = state.liquid_water_volume_m3[0] + state.liquid_water_volume_m3[1];
    try std.testing.expectApproxEqAbs(before - 0.04, after, 1e-12);
}

test "failed snow melt commit leaves liquid state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.01}, &.{1}, &.{270}, &.{0.05}, 0.05);
    state.liquid_water_volume_m3[0] = 0.01;
    try std.testing.expectError(error.SnowMeltExceedsLiquidInventory, state.commitMeltWater(&.{0}, &.{0.02}, &.{0}, &.{0}));
    try std.testing.expectEqual(@as(f64, 0.01), state.liquid_water_volume_m3[0]);
}

test "runtime snow layers route melt conservatively to lower layer and surface" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    state.active[0] = true;
    state.active[1] = true;
    state.liquid_water_volume_m3[0] = 2;
    state.liquid_water_volume_m3[1] = 2;
    @memset(try state.amounts(0, 0), 10);
    @memset(try state.amounts(0, 1), 4);
    const downward_water = [_]f64{ 0, 1, 0 };
    const partition = SurfacePartition{ .litter_cover_fraction = 0.25, .bare_soil_fraction = 0.75, .nonband_ammonium_fraction = 0.6, .band_ammonium_fraction = 0.4, .nonband_nitrate_fraction = 0.7, .band_nitrate_fraction = 0.3, .nonband_phosphate_fraction = 0.8, .band_phosphate_fraction = 0.2 };
    var fluxes = try calculateFluxes(std.testing.allocator, &state, &downward_water, &[_]f64{0.5}, &[_]f64{0.5}, &[_]f64{0}, &[_]SurfacePartition{partition});
    defer fluxes.deinit();
    const zero_input = [_]f64{0} ** species_count;
    try commit(&state, &zero_input, &fluxes);
    var remaining: f64 = 0;
    for (state.amount_g) |amount| remaining += amount;
    var discharged: f64 = 0;
    for (fluxes.surface_discharge[0].litter_g, fluxes.surface_discharge[0].soil_nonband_g, fluxes.surface_discharge[0].soil_band_g) |a, b, c| discharged += a + b + c;
    try std.testing.expectApproxEqAbs(@as(f64, 14 * species_count), remaining + discharged, 1e-12);
}

test "failed snow commit does not modify state" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    state.active[0] = true;
    state.amount_g[0] = 1;
    const downward = try std.testing.allocator.alloc(f64, species_count);
    defer std.testing.allocator.free(downward);
    @memset(downward, 0);
    const discharge = try std.testing.allocator.alloc(SurfaceDischarge, 1);
    defer std.testing.allocator.free(discharge);
    discharge[0] = .{};
    discharge[0].litter_g[0] = 2;
    const fluxes = Fluxes{ .allocator = std.testing.allocator, .downward_g = downward, .surface_discharge = discharge };
    const zero_input = [_]f64{0} ** species_count;
    try std.testing.expectError(error.InvalidSnowTransportCandidate, commit(&state, &zero_input, &fluxes));
    try std.testing.expectEqual(@as(f64, 1), state.amount_g[0]);
}
