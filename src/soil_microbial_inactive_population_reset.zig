const std = @import("std");

pub const SurfaceLayer = enum { subsurface, surface };
pub const MethanogenicComplex = enum { other, fifth };
pub const HydrogenotrophPopulation = enum { other, fifth };

pub const MetabolicFluxes = struct {
    oxygen_uptake_g_o: []f64,
    aerobic_respiration_g_c: []f64,
    carbon_dioxide_production_g_c: []f64,
    methane_reduction_g_c: []f64,
    methane_oxidation_g_c: []f64,
    oxygen_unlimited_respiration_g_c: []f64,
    denitrification_respiration_g_c: []f64,
    total_carbon_uptake_g_c: []f64,
    dissolved_organic_nitrogen_uptake_g_n: []f64,
    dissolved_organic_phosphorus_uptake_g_p: []f64,
    primary_carbon_substrate_uptake_g_c: []f64,
    acetate_uptake_g_c: []f64,
};

pub const NitrogenFluxes = struct {
    nitrate_reduction_g_n: []f64,
    band_nitrate_reduction_g_n: []f64,
    nitrite_reduction_g_n: []f64,
    band_nitrite_reduction_g_n: []f64,
    nitrous_oxide_reduction_g_n: []f64,
    fixed_dinitrogen_g_n: []f64,
};

pub const MineralFluxes = struct {
    ammonium_exchange_g_n: []f64,
    nitrate_exchange_g_n: []f64,
    dihydrogen_phosphate_exchange_g_p: []f64,
    hydrogen_phosphate_exchange_g_p: []f64,
    band_ammonium_exchange_g_n: []f64,
    band_nitrate_exchange_g_n: []f64,
    band_dihydrogen_phosphate_exchange_g_p: []f64,
    band_hydrogen_phosphate_exchange_g_p: []f64,
};

pub const SurfaceFluxes = struct {
    ammonium_exchange_g_n: []f64,
    nitrate_exchange_g_n: []f64,
    dihydrogen_phosphate_exchange_g_p: []f64,
    hydrogen_phosphate_exchange_g_p: []f64,
    ammonium_competition_fraction: []f64,
    nitrate_competition_fraction: []f64,
    dihydrogen_phosphate_competition_fraction: []f64,
    hydrogen_phosphate_competition_fraction: []f64,
};

pub const ComponentAssimilationFluxes = struct {
    carbon_assimilation_g_c: []f64,
    nitrogen_assimilation_g_n: []f64,
    phosphorus_assimilation_g_p: []f64,
    maintenance_respiration_g_c: []f64,
};

pub const ComponentSenescenceFluxes = struct {
    senesced_carbon_g_c: []f64,
    senesced_nitrogen_g_n: []f64,
    senesced_phosphorus_g_p: []f64,
    senescence_litterfall_carbon_g_c: []f64,
    senescence_litterfall_nitrogen_g_n: []f64,
    senescence_litterfall_phosphorus_g_p: []f64,
    senescence_recycled_carbon_g_c: []f64,
    senescence_recycled_nitrogen_g_n: []f64,
    senescence_recycled_phosphorus_g_p: []f64,
    senescence_humus_carbon_g_c: []f64,
    senescence_humus_nitrogen_g_n: []f64,
    senescence_humus_phosphorus_g_p: []f64,
    senescence_residue_carbon_g_c: []f64,
    senescence_residue_nitrogen_g_n: []f64,
    senescence_residue_phosphorus_g_p: []f64,
};

pub const ComponentDecompositionFluxes = struct {
    decomposed_carbon_g_c: []f64,
    decomposed_nitrogen_g_n: []f64,
    decomposed_phosphorus_g_p: []f64,
    decomposition_litterfall_carbon_g_c: []f64,
    decomposition_litterfall_nitrogen_g_n: []f64,
    decomposition_litterfall_phosphorus_g_p: []f64,
    decomposition_recycled_carbon_g_c: []f64,
    decomposition_recycled_nitrogen_g_n: []f64,
    decomposition_recycled_phosphorus_g_p: []f64,
    decomposition_humus_carbon_g_c: []f64,
    decomposition_humus_nitrogen_g_n: []f64,
    decomposition_humus_phosphorus_g_p: []f64,
    decomposition_residue_carbon_g_c: []f64,
    decomposition_residue_nitrogen_g_n: []f64,
    decomposition_residue_phosphorus_g_p: []f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    component_count: usize,
    population_count: usize,
    metabolic: MetabolicFluxes,
    nitrogen: NitrogenFluxes,
    mineral: MineralFluxes,
    surface: SurfaceFluxes,
    component_assimilation: ComponentAssimilationFluxes,
    component_senescence: ComponentSenescenceFluxes,
    component_decomposition: ComponentDecompositionFluxes,
    hydrogen_production_g_h: []f64,
    fifth_complex_non_band_oxidation_g_n: []f64,
    fifth_complex_band_oxidation_g_n: []f64,
    fifth_population_hydrogen_g_h: f64,

    pub fn init(
        allocator: std.mem.Allocator,
        unit_count: usize,
        component_count: usize,
        population_count: usize,
    ) !State {
        if (unit_count == 0 or component_count == 0 or population_count == 0)
            return error.InvalidInactivePopulationResetDimensions;
        const item_count = std.math.mul(usize, unit_count, component_count) catch
            return error.InvalidInactivePopulationResetDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.component_count = component_count;
        state.population_count = population_count;
        state.metabolic = try allocateGroup(MetabolicFluxes, allocator, unit_count);
        errdefer freeGroup(MetabolicFluxes, allocator, &state.metabolic);
        state.nitrogen = try allocateGroup(NitrogenFluxes, allocator, unit_count);
        errdefer freeGroup(NitrogenFluxes, allocator, &state.nitrogen);
        state.mineral = try allocateGroup(MineralFluxes, allocator, unit_count);
        errdefer freeGroup(MineralFluxes, allocator, &state.mineral);
        state.surface = try allocateGroup(SurfaceFluxes, allocator, unit_count);
        errdefer freeGroup(SurfaceFluxes, allocator, &state.surface);
        state.component_assimilation =
            try allocateGroup(ComponentAssimilationFluxes, allocator, item_count);
        errdefer freeGroup(
            ComponentAssimilationFluxes,
            allocator,
            &state.component_assimilation,
        );
        state.component_senescence =
            try allocateGroup(ComponentSenescenceFluxes, allocator, item_count);
        errdefer freeGroup(
            ComponentSenescenceFluxes,
            allocator,
            &state.component_senescence,
        );
        state.component_decomposition =
            try allocateGroup(ComponentDecompositionFluxes, allocator, item_count);
        errdefer freeGroup(
            ComponentDecompositionFluxes,
            allocator,
            &state.component_decomposition,
        );
        state.hydrogen_production_g_h = try allocator.alloc(f64, unit_count);
        errdefer allocator.free(state.hydrogen_production_g_h);
        state.fifth_complex_non_band_oxidation_g_n =
            try allocator.alloc(f64, population_count);
        errdefer allocator.free(state.fifth_complex_non_band_oxidation_g_n);
        state.fifth_complex_band_oxidation_g_n =
            try allocator.alloc(f64, population_count);
        @memset(state.hydrogen_production_g_h, 0);
        @memset(state.fifth_complex_non_band_oxidation_g_n, 0);
        @memset(state.fifth_complex_band_oxidation_g_n, 0);
        state.fifth_population_hydrogen_g_h = 0;
        return state;
    }

    pub fn deinit(self: *State) void {
        freeGroup(MetabolicFluxes, self.allocator, &self.metabolic);
        freeGroup(NitrogenFluxes, self.allocator, &self.nitrogen);
        freeGroup(MineralFluxes, self.allocator, &self.mineral);
        freeGroup(SurfaceFluxes, self.allocator, &self.surface);
        freeGroup(
            ComponentAssimilationFluxes,
            self.allocator,
            &self.component_assimilation,
        );
        freeGroup(
            ComponentSenescenceFluxes,
            self.allocator,
            &self.component_senescence,
        );
        freeGroup(
            ComponentDecompositionFluxes,
            self.allocator,
            &self.component_decomposition,
        );
        self.allocator.free(self.hydrogen_production_g_h);
        self.allocator.free(self.fifth_complex_non_band_oxidation_g_n);
        self.allocator.free(self.fifth_complex_band_oxidation_g_n);
        self.* = undefined;
    }
};

pub const Target = struct {
    unit: usize,
    population: usize,
    layer: SurfaceLayer,
    complex: MethanogenicComplex,
    hydrogenotroph: HydrogenotrophPopulation,
};

/// Exact NITRO.F 2833--2918 inactive microbial-population zero publication.
pub fn reset(state: *State, target: Target) !void {
    if (target.unit >= state.unit_count or target.population >= state.population_count)
        return error.InactivePopulationResetIndexOutOfBounds;
    zeroGroupItem(MetabolicFluxes, &state.metabolic, target.unit);
    zeroGroupItem(NitrogenFluxes, &state.nitrogen, target.unit);
    zeroGroupItem(MineralFluxes, &state.mineral, target.unit);
    if (target.layer == .surface)
        zeroGroupItem(SurfaceFluxes, &state.surface, target.unit);
    for (0..state.component_count) |component| {
        const item = target.unit * state.component_count + component;
        zeroGroupItem(
            ComponentAssimilationFluxes,
            &state.component_assimilation,
            item,
        );
        zeroGroupItem(ComponentSenescenceFluxes, &state.component_senescence, item);
        zeroGroupItem(
            ComponentDecompositionFluxes,
            &state.component_decomposition,
            item,
        );
    }
    state.hydrogen_production_g_h[target.unit] = 0;
    if (target.complex == .fifth) {
        state.fifth_complex_non_band_oxidation_g_n[target.population] = 0;
        state.fifth_complex_band_oxidation_g_n[target.population] = 0;
        if (target.hydrogenotroph == .fifth)
            state.fifth_population_hydrogen_g_h = 0;
    }
}

fn allocateGroup(
    comptime T: type,
    allocator: std.mem.Allocator,
    count: usize,
) !T {
    var result: T = undefined;
    var allocated: usize = 0;
    errdefer {
        inline for (@typeInfo(T).@"struct".fields) |field| {
            if (allocated > 0) {
                allocated -= 1;
                allocator.free(@field(result, field.name));
            }
        }
    }
    inline for (@typeInfo(T).@"struct".fields) |field| {
        @field(result, field.name) = try allocator.alloc(f64, count);
        @memset(@field(result, field.name), 0);
        allocated += 1;
    }
    return result;
}

fn freeGroup(comptime T: type, allocator: std.mem.Allocator, group: *T) void {
    inline for (@typeInfo(T).@"struct".fields) |field|
        allocator.free(@field(group, field.name));
}

fn zeroGroupItem(comptime T: type, group: *T, index: usize) void {
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(group, field.name)[index] = 0;
}

fn fillGroup(comptime T: type, group: *T, value: f64) void {
    inline for (@typeInfo(T).@"struct".fields) |field|
        @memset(@field(group, field.name), value);
}

test "NITRO 2833-2918 resets complete inactive surface population" {
    var state = try State.init(std.testing.allocator, 2, 2, 7);
    defer state.deinit();
    fillGroup(MetabolicFluxes, &state.metabolic, 1);
    fillGroup(NitrogenFluxes, &state.nitrogen, 1);
    fillGroup(MineralFluxes, &state.mineral, 1);
    fillGroup(SurfaceFluxes, &state.surface, 1);
    fillGroup(ComponentAssimilationFluxes, &state.component_assimilation, 1);
    fillGroup(ComponentSenescenceFluxes, &state.component_senescence, 1);
    fillGroup(ComponentDecompositionFluxes, &state.component_decomposition, 1);
    @memset(state.hydrogen_production_g_h, 1);
    @memset(state.fifth_complex_non_band_oxidation_g_n, 1);
    @memset(state.fifth_complex_band_oxidation_g_n, 1);
    state.fifth_population_hydrogen_g_h = 1;
    try reset(&state, .{
        .unit = 0,
        .population = 4,
        .layer = .surface,
        .complex = .fifth,
        .hydrogenotroph = .fifth,
    });
    inline for (@typeInfo(MetabolicFluxes).@"struct".fields) |field|
        try std.testing.expectEqual(0, @field(state.metabolic, field.name)[0]);
    inline for (@typeInfo(SurfaceFluxes).@"struct".fields) |field|
        try std.testing.expectEqual(0, @field(state.surface, field.name)[0]);
    inline for (@typeInfo(ComponentAssimilationFluxes).@"struct".fields) |field|
        for (0..2) |item| try std.testing.expectEqual(
            0,
            @field(state.component_assimilation, field.name)[item],
        );
    inline for (@typeInfo(ComponentSenescenceFluxes).@"struct".fields) |field|
        for (0..2) |item| try std.testing.expectEqual(
            0,
            @field(state.component_senescence, field.name)[item],
        );
    inline for (@typeInfo(ComponentDecompositionFluxes).@"struct".fields) |field|
        for (0..2) |item| try std.testing.expectEqual(
            0,
            @field(state.component_decomposition, field.name)[item],
        );
    try std.testing.expectEqual(0, state.fifth_population_hydrogen_g_h);
    try std.testing.expectEqual(1, state.hydrogen_production_g_h[1]);
}

test "NITRO inactive subsurface reset preserves surface-only fields" {
    var state = try State.init(std.testing.allocator, 1, 2, 7);
    defer state.deinit();
    fillGroup(SurfaceFluxes, &state.surface, 3);
    try reset(&state, .{
        .unit = 0,
        .population = 0,
        .layer = .subsurface,
        .complex = .other,
        .hydrogenotroph = .other,
    });
    inline for (@typeInfo(SurfaceFluxes).@"struct".fields) |field|
        try std.testing.expectEqual(3, @field(state.surface, field.name)[0]);
}

test "NITRO inactive reset invalid target is atomic" {
    var state = try State.init(std.testing.allocator, 1, 2, 7);
    defer state.deinit();
    state.metabolic.oxygen_uptake_g_o[0] = 7;
    try std.testing.expectError(
        error.InactivePopulationResetIndexOutOfBounds,
        reset(&state, .{
            .unit = 1,
            .population = 0,
            .layer = .surface,
            .complex = .fifth,
            .hydrogenotroph = .fifth,
        }),
    );
    try std.testing.expectEqual(7, state.metabolic.oxygen_uptake_g_o[0]);
}
