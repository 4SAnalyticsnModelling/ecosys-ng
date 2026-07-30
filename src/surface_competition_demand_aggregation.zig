const std = @import("std");

pub const organic_substrate_count: usize = 5;
pub const microbial_substrate_count: usize = 6;
pub const surface_litter_substrate_count: usize = 3;
pub const autotrophic_substrate_index: usize = 5;

pub const PopulationDemand = struct {
    oxygen_g_o_per_h: f64 = 0,
    ammonium_oxidation_g_n_per_h: f64 = 0,
    nitrate_reduction_g_n_per_h: f64 = 0,
    nitrite_oxidation_g_n_per_h: f64 = 0,
    nitrous_oxide_reduction_g_n_per_h: f64 = 0,
    surface_ammonium_immobilization_g_n_per_h: f64 = 0,
    surface_nitrate_immobilization_g_n_per_h: f64 = 0,
    surface_h2po4_immobilization_g_p_per_h: f64 = 0,
    surface_hpo4_immobilization_g_p_per_h: f64 = 0,
    topsoil_ammonium_immobilization_g_n_per_h: f64 = 0,
    topsoil_nitrate_immobilization_g_n_per_h: f64 = 0,
    topsoil_h2po4_immobilization_g_p_per_h: f64 = 0,
    topsoil_hpo4_immobilization_g_p_per_h: f64 = 0,
    dissolved_organic_carbon_demand_g_c_per_h: f64 = 0,
    dissolved_acetate_carbon_demand_g_c_per_h: f64 = 0,
};

pub const SurfaceDemand = struct {
    oxygen_g_o_per_h: f64 = 0,
    ammonium_g_n_per_h: f64 = 0,
    nitrate_g_n_per_h: f64 = 0,
    nitrite_g_n_per_h: f64 = 0,
    nitrous_oxide_g_n_per_h: f64 = 0,
    h2po4_g_p_per_h: f64 = 0,
    hpo4_g_p_per_h: f64 = 0,
};

pub const TopsoilDemand = struct {
    ammonium_g_n_per_h: f64 = 0,
    nitrate_g_n_per_h: f64 = 0,
    h2po4_g_p_per_h: f64 = 0,
    hpo4_g_p_per_h: f64 = 0,
};

pub const DissolvedCarbonDemand = struct {
    organic_carbon_g_c_per_h: f64 = 0,
    acetate_carbon_g_c_per_h: f64 = 0,
};

pub const Inputs = struct {
    column_count: usize,
    row_count: usize,
    soil_layer_capacity: usize,
    population_count: usize,
    top_soil_layer_by_cell: []const usize,
    population_demand_by_cell_substrate: []const PopulationDemand,
    chemodenitrification_nitrite_demand_g_n_per_h_by_cell: []const f64,
};

pub const State = struct {
    surface_demand_by_cell: []SurfaceDemand,
    topsoil_demand_by_soil_layer: []TopsoilDemand,
    dissolved_carbon_demand_by_cell_substrate: []DissolvedCarbonDemand,
};

/// Aggregates surface microbial and autotrophic competition demands.
///
/// Traceability: REDIST.F lines 308--331 (`ROXYX`, `RNH4X`, `RNO3X`,
/// `RNO2X`, `RN2OX`, `RPO4X`, `RP14X`, `ROQCX`, and `ROQAX`). The source
/// admits the three litter substrates and the autotrophic substrate while
/// skipping particulate and humus substrate slots. Runtime traversal retains
/// column, row, substrate, population, then source statement order. All
/// accepted-hour demand fields are extensive g element h-1.
pub fn apply(inputs: Inputs, state: *State) !void {
    const cell_count = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, cell_count);
    try preflightUpdates(inputs, state.*);

    for (0..inputs.column_count) |column| {
        for (0..inputs.row_count) |row| {
            const cell = row * inputs.column_count + column;
            const topsoil_index =
                cell * inputs.soil_layer_capacity + inputs.top_soil_layer_by_cell[cell];
            for (0..microbial_substrate_count) |substrate| {
                if (!admittedSubstrate(substrate)) continue;
                for (0..inputs.population_count) |population| {
                    const demand = inputs.population_demand_by_cell_substrate[
                        demandIndex(cell, substrate, population, inputs.population_count)
                    ];
                    commitPopulation(
                        &state.surface_demand_by_cell[cell],
                        &state.topsoil_demand_by_soil_layer[topsoil_index],
                        if (substrate < organic_substrate_count)
                            &state.dissolved_carbon_demand_by_cell_substrate[
                                cell * organic_substrate_count + substrate
                            ]
                        else
                            null,
                        demand,
                    );
                }
            }
            state.surface_demand_by_cell[cell].nitrite_g_n_per_h +=
                inputs.chemodenitrification_nitrite_demand_g_n_per_h_by_cell[cell];
        }
    }
}

fn validateDimensions(inputs: Inputs, state: State) !usize {
    if (inputs.column_count == 0 or
        inputs.row_count == 0 or
        inputs.soil_layer_capacity == 0 or
        inputs.population_count == 0)
        return error.InvalidSurfaceCompetitionDimensions;
    const cell_count = checkedProduct(
        &.{ inputs.column_count, inputs.row_count },
    ) catch return error.InvalidSurfaceCompetitionDimensions;
    const demand_count = checkedProduct(
        &.{ cell_count, microbial_substrate_count, inputs.population_count },
    ) catch return error.InvalidSurfaceCompetitionDimensions;
    const soil_value_count = checkedProduct(
        &.{ cell_count, inputs.soil_layer_capacity },
    ) catch return error.InvalidSurfaceCompetitionDimensions;
    const dissolved_value_count = checkedProduct(
        &.{ cell_count, organic_substrate_count },
    ) catch return error.InvalidSurfaceCompetitionDimensions;
    if (inputs.top_soil_layer_by_cell.len != cell_count or
        inputs.population_demand_by_cell_substrate.len != demand_count or
        inputs.chemodenitrification_nitrite_demand_g_n_per_h_by_cell.len != cell_count or
        state.surface_demand_by_cell.len != cell_count or
        state.topsoil_demand_by_soil_layer.len != soil_value_count or
        state.dissolved_carbon_demand_by_cell_substrate.len != dissolved_value_count)
        return error.InvalidSurfaceCompetitionDimensions;
    return cell_count;
}

fn validateInputsAndState(inputs: Inputs, state: State, cell_count: usize) !void {
    for (inputs.population_demand_by_cell_substrate) |demand|
        try validateStruct(demand, error.InvalidSurfaceCompetitionInput);
    for (inputs.chemodenitrification_nitrite_demand_g_n_per_h_by_cell) |demand|
        if (!nonnegativeFinite(demand))
            return error.InvalidSurfaceCompetitionInput;
    for (state.surface_demand_by_cell) |demand|
        try validateStruct(demand, error.InvalidSurfaceCompetitionState);
    for (state.topsoil_demand_by_soil_layer) |demand|
        try validateStruct(demand, error.InvalidSurfaceCompetitionState);
    for (state.dissolved_carbon_demand_by_cell_substrate) |demand|
        try validateStruct(demand, error.InvalidSurfaceCompetitionState);
    for (0..cell_count) |cell|
        if (inputs.top_soil_layer_by_cell[cell] >= inputs.soil_layer_capacity)
            return error.InvalidSurfaceCompetitionInput;
}

fn preflightUpdates(inputs: Inputs, state: State) !void {
    for (0..inputs.column_count) |column| {
        for (0..inputs.row_count) |row| {
            const cell = row * inputs.column_count + column;
            const topsoil_index =
                cell * inputs.soil_layer_capacity + inputs.top_soil_layer_by_cell[cell];
            var surface = state.surface_demand_by_cell[cell];
            var topsoil = state.topsoil_demand_by_soil_layer[topsoil_index];
            for (0..microbial_substrate_count) |substrate| {
                if (!admittedSubstrate(substrate)) continue;
                var dissolved = if (substrate < organic_substrate_count)
                    state.dissolved_carbon_demand_by_cell_substrate[
                        cell * organic_substrate_count + substrate
                    ]
                else
                    DissolvedCarbonDemand{};
                for (0..inputs.population_count) |population| {
                    const demand = inputs.population_demand_by_cell_substrate[
                        demandIndex(cell, substrate, population, inputs.population_count)
                    ];
                    try accumulatePopulation(
                        &surface,
                        &topsoil,
                        if (substrate < organic_substrate_count) &dissolved else null,
                        demand,
                    );
                }
            }
            surface.nitrite_g_n_per_h = try checkedSum(
                surface.nitrite_g_n_per_h,
                inputs.chemodenitrification_nitrite_demand_g_n_per_h_by_cell[cell],
            );
        }
    }
}

fn accumulatePopulation(
    surface: *SurfaceDemand,
    topsoil: *TopsoilDemand,
    dissolved: ?*DissolvedCarbonDemand,
    demand: PopulationDemand,
) !void {
    surface.oxygen_g_o_per_h =
        try checkedSum(surface.oxygen_g_o_per_h, demand.oxygen_g_o_per_h);
    surface.ammonium_g_n_per_h = try checkedSum(
        surface.ammonium_g_n_per_h,
        demand.ammonium_oxidation_g_n_per_h,
    );
    surface.nitrate_g_n_per_h = try checkedSum(
        surface.nitrate_g_n_per_h,
        demand.nitrate_reduction_g_n_per_h,
    );
    surface.nitrite_g_n_per_h = try checkedSum(
        surface.nitrite_g_n_per_h,
        demand.nitrite_oxidation_g_n_per_h,
    );
    surface.nitrous_oxide_g_n_per_h = try checkedSum(
        surface.nitrous_oxide_g_n_per_h,
        demand.nitrous_oxide_reduction_g_n_per_h,
    );
    surface.ammonium_g_n_per_h = try checkedSum(
        surface.ammonium_g_n_per_h,
        demand.surface_ammonium_immobilization_g_n_per_h,
    );
    surface.nitrate_g_n_per_h = try checkedSum(
        surface.nitrate_g_n_per_h,
        demand.surface_nitrate_immobilization_g_n_per_h,
    );
    surface.h2po4_g_p_per_h = try checkedSum(
        surface.h2po4_g_p_per_h,
        demand.surface_h2po4_immobilization_g_p_per_h,
    );
    surface.hpo4_g_p_per_h = try checkedSum(
        surface.hpo4_g_p_per_h,
        demand.surface_hpo4_immobilization_g_p_per_h,
    );
    topsoil.ammonium_g_n_per_h = try checkedSum(
        topsoil.ammonium_g_n_per_h,
        demand.topsoil_ammonium_immobilization_g_n_per_h,
    );
    topsoil.nitrate_g_n_per_h = try checkedSum(
        topsoil.nitrate_g_n_per_h,
        demand.topsoil_nitrate_immobilization_g_n_per_h,
    );
    topsoil.h2po4_g_p_per_h = try checkedSum(
        topsoil.h2po4_g_p_per_h,
        demand.topsoil_h2po4_immobilization_g_p_per_h,
    );
    topsoil.hpo4_g_p_per_h = try checkedSum(
        topsoil.hpo4_g_p_per_h,
        demand.topsoil_hpo4_immobilization_g_p_per_h,
    );
    if (dissolved) |value| {
        value.organic_carbon_g_c_per_h = try checkedSum(
            value.organic_carbon_g_c_per_h,
            demand.dissolved_organic_carbon_demand_g_c_per_h,
        );
        value.acetate_carbon_g_c_per_h = try checkedSum(
            value.acetate_carbon_g_c_per_h,
            demand.dissolved_acetate_carbon_demand_g_c_per_h,
        );
    }
}

fn commitPopulation(
    surface: *SurfaceDemand,
    topsoil: *TopsoilDemand,
    dissolved: ?*DissolvedCarbonDemand,
    demand: PopulationDemand,
) void {
    surface.oxygen_g_o_per_h += demand.oxygen_g_o_per_h;
    surface.ammonium_g_n_per_h += demand.ammonium_oxidation_g_n_per_h;
    surface.nitrate_g_n_per_h += demand.nitrate_reduction_g_n_per_h;
    surface.nitrite_g_n_per_h += demand.nitrite_oxidation_g_n_per_h;
    surface.nitrous_oxide_g_n_per_h += demand.nitrous_oxide_reduction_g_n_per_h;
    surface.ammonium_g_n_per_h +=
        demand.surface_ammonium_immobilization_g_n_per_h;
    surface.nitrate_g_n_per_h += demand.surface_nitrate_immobilization_g_n_per_h;
    surface.h2po4_g_p_per_h += demand.surface_h2po4_immobilization_g_p_per_h;
    surface.hpo4_g_p_per_h += demand.surface_hpo4_immobilization_g_p_per_h;
    topsoil.ammonium_g_n_per_h += demand.topsoil_ammonium_immobilization_g_n_per_h;
    topsoil.nitrate_g_n_per_h += demand.topsoil_nitrate_immobilization_g_n_per_h;
    topsoil.h2po4_g_p_per_h += demand.topsoil_h2po4_immobilization_g_p_per_h;
    topsoil.hpo4_g_p_per_h += demand.topsoil_hpo4_immobilization_g_p_per_h;
    if (dissolved) |value| {
        value.organic_carbon_g_c_per_h +=
            demand.dissolved_organic_carbon_demand_g_c_per_h;
        value.acetate_carbon_g_c_per_h +=
            demand.dissolved_acetate_carbon_demand_g_c_per_h;
    }
}

fn admittedSubstrate(substrate: usize) bool {
    return substrate < surface_litter_substrate_count or
        substrate == autotrophic_substrate_index;
}

fn demandIndex(
    cell: usize,
    substrate: usize,
    population: usize,
    population_count: usize,
) usize {
    return (cell * microbial_substrate_count + substrate) * population_count +
        population;
}

fn checkedProduct(values: []const usize) !usize {
    var result: usize = 1;
    for (values) |value|
        result = std.math.mul(usize, result, value) catch return error.Overflow;
    return result;
}

fn validateStruct(value: anytype, comptime failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!nonnegativeFinite(@field(value, field.name))) return failure;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!nonnegativeFinite(result))
        return error.NonFiniteSurfaceCompetitionResult;
    return result;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

fn demandWithScale(scale: f64) PopulationDemand {
    return .{
        .oxygen_g_o_per_h = 1 * scale,
        .ammonium_oxidation_g_n_per_h = 2 * scale,
        .nitrate_reduction_g_n_per_h = 3 * scale,
        .nitrite_oxidation_g_n_per_h = 4 * scale,
        .nitrous_oxide_reduction_g_n_per_h = 5 * scale,
        .surface_ammonium_immobilization_g_n_per_h = 6 * scale,
        .surface_nitrate_immobilization_g_n_per_h = 7 * scale,
        .surface_h2po4_immobilization_g_p_per_h = 8 * scale,
        .surface_hpo4_immobilization_g_p_per_h = 9 * scale,
        .topsoil_ammonium_immobilization_g_n_per_h = 10 * scale,
        .topsoil_nitrate_immobilization_g_n_per_h = 11 * scale,
        .topsoil_h2po4_immobilization_g_p_per_h = 12 * scale,
        .topsoil_hpo4_immobilization_g_p_per_h = 13 * scale,
        .dissolved_organic_carbon_demand_g_c_per_h = 14 * scale,
        .dissolved_acetate_carbon_demand_g_c_per_h = 15 * scale,
    };
}

test "REDIST surface competition preserves every source addition and skip" {
    const population_count = 2;
    const demand_count = microbial_substrate_count * population_count;
    var demands = [_]PopulationDemand{demandWithScale(1)} ** demand_count;
    for (0..population_count) |population| {
        demands[demandIndex(0, 3, population, population_count)] =
            demandWithScale(100);
        demands[demandIndex(0, 4, population, population_count)] =
            demandWithScale(100);
    }
    var surface = [_]SurfaceDemand{.{}};
    var topsoil = [_]TopsoilDemand{ .{}, .{} };
    var dissolved = [_]DissolvedCarbonDemand{.{}} ** organic_substrate_count;
    var state: State = .{
        .surface_demand_by_cell = &surface,
        .topsoil_demand_by_soil_layer = &topsoil,
        .dissolved_carbon_demand_by_cell_substrate = &dissolved,
    };
    try apply(.{
        .column_count = 1,
        .row_count = 1,
        .soil_layer_capacity = 2,
        .population_count = population_count,
        .top_soil_layer_by_cell = &.{1},
        .population_demand_by_cell_substrate = &demands,
        .chemodenitrification_nitrite_demand_g_n_per_h_by_cell = &.{5},
    }, &state);

    try std.testing.expectEqual(@as(f64, 8), surface[0].oxygen_g_o_per_h);
    try std.testing.expectEqual(@as(f64, 64), surface[0].ammonium_g_n_per_h);
    try std.testing.expectEqual(@as(f64, 80), surface[0].nitrate_g_n_per_h);
    try std.testing.expectEqual(@as(f64, 37), surface[0].nitrite_g_n_per_h);
    try std.testing.expectEqual(@as(f64, 40), surface[0].nitrous_oxide_g_n_per_h);
    try std.testing.expectEqual(@as(f64, 64), surface[0].h2po4_g_p_per_h);
    try std.testing.expectEqual(@as(f64, 72), surface[0].hpo4_g_p_per_h);
    try std.testing.expectEqual(TopsoilDemand{}, topsoil[0]);
    try std.testing.expectEqual(@as(f64, 80), topsoil[1].ammonium_g_n_per_h);
    try std.testing.expectEqual(@as(f64, 88), topsoil[1].nitrate_g_n_per_h);
    try std.testing.expectEqual(@as(f64, 96), topsoil[1].h2po4_g_p_per_h);
    try std.testing.expectEqual(@as(f64, 104), topsoil[1].hpo4_g_p_per_h);
    for (0..surface_litter_substrate_count) |substrate| {
        try std.testing.expectEqual(
            @as(f64, 28),
            dissolved[substrate].organic_carbon_g_c_per_h,
        );
        try std.testing.expectEqual(
            @as(f64, 30),
            dissolved[substrate].acetate_carbon_g_c_per_h,
        );
    }
    try std.testing.expectEqual(DissolvedCarbonDemand{}, dissolved[3]);
    try std.testing.expectEqual(DissolvedCarbonDemand{}, dissolved[4]);
}

test "runtime population and grid extents route each topsoil cell" {
    const cell_count = 4;
    const population_count = 3;
    const demand_count =
        cell_count * microbial_substrate_count * population_count;
    var demands = [_]PopulationDemand{.{}} ** demand_count;
    demands[demandIndex(3, autotrophic_substrate_index, 2, population_count)] =
        demandWithScale(2);
    var surface = [_]SurfaceDemand{.{}} ** cell_count;
    var topsoil = [_]TopsoilDemand{.{}} ** (cell_count * 3);
    var dissolved =
        [_]DissolvedCarbonDemand{.{}} ** (cell_count * organic_substrate_count);
    var state: State = .{
        .surface_demand_by_cell = &surface,
        .topsoil_demand_by_soil_layer = &topsoil,
        .dissolved_carbon_demand_by_cell_substrate = &dissolved,
    };
    try apply(.{
        .column_count = 2,
        .row_count = 2,
        .soil_layer_capacity = 3,
        .population_count = population_count,
        .top_soil_layer_by_cell = &.{ 2, 1, 0, 2 },
        .population_demand_by_cell_substrate = &demands,
        .chemodenitrification_nitrite_demand_g_n_per_h_by_cell = &.{ 0, 0, 0, 0 },
    }, &state);

    try std.testing.expectEqual(@as(f64, 2), surface[3].oxygen_g_o_per_h);
    try std.testing.expectEqual(
        @as(f64, 20),
        topsoil[3 * 3 + 2].ammonium_g_n_per_h,
    );
    for (dissolved) |demand|
        try std.testing.expectEqual(DissolvedCarbonDemand{}, demand);
}

test "late invalid competition demand leaves all state unchanged" {
    const cell_count = 2;
    const population_count = 2;
    const demand_count =
        cell_count * microbial_substrate_count * population_count;
    var demands = [_]PopulationDemand{.{}} ** demand_count;
    demands[0] = demandWithScale(1);
    demands[demand_count - 1].topsoil_hpo4_immobilization_g_p_per_h =
        std.math.nan(f64);
    var surface = [_]SurfaceDemand{.{ .oxygen_g_o_per_h = 1 }} ** cell_count;
    var topsoil = [_]TopsoilDemand{.{ .ammonium_g_n_per_h = 2 }} ** cell_count;
    var dissolved =
        [_]DissolvedCarbonDemand{.{ .organic_carbon_g_c_per_h = 3 }} **
        (cell_count * organic_substrate_count);
    var state: State = .{
        .surface_demand_by_cell = &surface,
        .topsoil_demand_by_soil_layer = &topsoil,
        .dissolved_carbon_demand_by_cell_substrate = &dissolved,
    };
    try std.testing.expectError(error.InvalidSurfaceCompetitionInput, apply(.{
        .column_count = 2,
        .row_count = 1,
        .soil_layer_capacity = 1,
        .population_count = population_count,
        .top_soil_layer_by_cell = &.{ 0, 0 },
        .population_demand_by_cell_substrate = &demands,
        .chemodenitrification_nitrite_demand_g_n_per_h_by_cell = &.{ 0, 0 },
    }, &state));
    for (surface) |demand|
        try std.testing.expectEqual(@as(f64, 1), demand.oxygen_g_o_per_h);
    for (topsoil) |demand|
        try std.testing.expectEqual(@as(f64, 2), demand.ammonium_g_n_per_h);
    for (dissolved) |demand|
        try std.testing.expectEqual(
            @as(f64, 3),
            demand.organic_carbon_g_c_per_h,
        );
}

test "invalid runtime dimensions and top-layer index fail explicitly" {
    var surface = [_]SurfaceDemand{.{}};
    var topsoil = [_]TopsoilDemand{.{}};
    var dissolved = [_]DissolvedCarbonDemand{.{}} ** organic_substrate_count;
    var state: State = .{
        .surface_demand_by_cell = &surface,
        .topsoil_demand_by_soil_layer = &topsoil,
        .dissolved_carbon_demand_by_cell_substrate = &dissolved,
    };
    const demands = [_]PopulationDemand{.{}} ** microbial_substrate_count;
    try std.testing.expectError(error.InvalidSurfaceCompetitionDimensions, apply(.{
        .column_count = 1,
        .row_count = 1,
        .soil_layer_capacity = 1,
        .population_count = 0,
        .top_soil_layer_by_cell = &.{0},
        .population_demand_by_cell_substrate = &demands,
        .chemodenitrification_nitrite_demand_g_n_per_h_by_cell = &.{0},
    }, &state));
    try std.testing.expectError(error.InvalidSurfaceCompetitionInput, apply(.{
        .column_count = 1,
        .row_count = 1,
        .soil_layer_capacity = 1,
        .population_count = 1,
        .top_soil_layer_by_cell = &.{1},
        .population_demand_by_cell_substrate = &demands,
        .chemodenitrification_nitrite_demand_g_n_per_h_by_cell = &.{0},
    }, &state));
}
