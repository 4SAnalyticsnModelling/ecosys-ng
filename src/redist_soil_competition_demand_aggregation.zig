const std = @import("std");

/// One `(N,K,L)` contribution to competition demand (g timestep-1).
pub const PopulationDemand = struct {
    o2_g_o: f64, // ROXYS
    nh4_maximum_g_n: f64, // RVMX4
    nh4_immobilization_g_n: f64, // RINHO
    no3_maximum_g_n: f64, // RVMX3
    no3_immobilization_g_n: f64, // RINOO
    no2_g_n: f64, // RVMX2
    n2o_g_n: f64, // RVMX1
    h2po4_g_p: f64, // RIPOO
    hpo4_g_p: f64, // RIPO1
    band_nh4_maximum_g_n: f64, // RVMB4
    band_nh4_immobilization_g_n: f64, // RINHB
    band_no3_maximum_g_n: f64, // RVMB3
    band_no3_immobilization_g_n: f64, // RINOB
    band_no2_g_n: f64, // RVMB2
    band_h2po4_g_p: f64, // RIPBO
    band_hpo4_g_p: f64, // RIPB1
    doc_oxidation_g_c: f64, // ROQCS, only K <= 4
    acetate_oxidation_g_c: f64, // ROQAS, only K <= 4
};

/// Corrections added after both legacy nested loops finish.
pub const LayerCorrections = struct {
    no2_g_n: f64, // RVMXC
    band_no2_g_n: f64, // RVMBC
};

/// Total substrate demands for one soil layer.
pub const LayerDemand = struct {
    o2_g_o: f64, // ROXYX
    nh4_g_n: f64, // RNH4X
    no3_g_n: f64, // RNO3X
    no2_g_n: f64, // RNO2X
    n2o_g_n: f64, // RN2OX
    h2po4_g_p: f64, // RPO4X
    hpo4_g_p: f64, // RP14X
    band_nh4_g_n: f64, // RNHBX
    band_no3_g_n: f64, // RN3BX
    band_no2_g_n: f64, // RN2BX
    band_h2po4_g_p: f64, // RPOBX
    band_hpo4_g_p: f64, // RP1BX
};

pub const OrganicDemand = struct {
    doc_oxidation_g_c: f64, // ROQCX(K)
    acetate_oxidation_g_c: f64, // ROQAX(K)
};

fn finiteStruct(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 6626--6651.
///
/// Runtime slices use `[layer][K][N]` order. The legacy scientific dimensions
/// remain six K classes and seven N populations, but no storage is static.
pub fn aggregateLayers(
    demands_by_layer: []LayerDemand,
    organic_by_layer_and_class: []OrganicDemand,
    population_by_layer_class_population: []const PopulationDemand,
    corrections_by_layer: []const LayerCorrections,
) !void {
    const class_count: usize = 6;
    const population_count: usize = 7;
    const organic_class_count: usize = 5;
    if (demands_by_layer.len == 0 or
        corrections_by_layer.len != demands_by_layer.len or
        organic_by_layer_and_class.len != demands_by_layer.len * organic_class_count or
        population_by_layer_class_population.len != demands_by_layer.len * class_count * population_count)
        return error.SoilCompetitionDemandDimensionMismatch;

    for (demands_by_layer, corrections_by_layer, 0..) |*demand, correction, layer_index| {
        if (!finiteStruct(demand.*) or !finiteStruct(correction))
            return error.InvalidSoilCompetitionDemand;
        const organic_start = layer_index * organic_class_count;
        for (organic_by_layer_and_class[organic_start .. organic_start + organic_class_count]) |organic|
            if (!finiteStruct(organic)) return error.InvalidSoilCompetitionDemand;

        const population_start = layer_index * class_count * population_count;
        var next_demand = demand.*;
        var next_organic: [organic_class_count]OrganicDemand = undefined;
        for (0..organic_class_count) |class_index|
            next_organic[class_index] = organic_by_layer_and_class[organic_start + class_index];
        for (0..class_count) |class_index| {
            for (0..population_count) |population_index| {
                const index = population_start + class_index * population_count + population_index;
                const contribution = population_by_layer_class_population[index];
                if (!finiteStruct(contribution)) return error.InvalidSoilCompetitionDemand;
                next_demand.o2_g_o = next_demand.o2_g_o + contribution.o2_g_o;
                next_demand.nh4_g_n = next_demand.nh4_g_n + contribution.nh4_maximum_g_n + contribution.nh4_immobilization_g_n;
                next_demand.no3_g_n = next_demand.no3_g_n + contribution.no3_maximum_g_n + contribution.no3_immobilization_g_n;
                next_demand.no2_g_n = next_demand.no2_g_n + contribution.no2_g_n;
                next_demand.n2o_g_n = next_demand.n2o_g_n + contribution.n2o_g_n;
                next_demand.h2po4_g_p = next_demand.h2po4_g_p + contribution.h2po4_g_p;
                next_demand.hpo4_g_p = next_demand.hpo4_g_p + contribution.hpo4_g_p;
                next_demand.band_nh4_g_n = next_demand.band_nh4_g_n + contribution.band_nh4_maximum_g_n + contribution.band_nh4_immobilization_g_n;
                next_demand.band_no3_g_n = next_demand.band_no3_g_n + contribution.band_no3_maximum_g_n + contribution.band_no3_immobilization_g_n;
                next_demand.band_no2_g_n = next_demand.band_no2_g_n + contribution.band_no2_g_n;
                next_demand.band_h2po4_g_p = next_demand.band_h2po4_g_p + contribution.band_h2po4_g_p;
                next_demand.band_hpo4_g_p = next_demand.band_hpo4_g_p + contribution.band_hpo4_g_p;
                if (class_index <= 4) {
                    const organic = &next_organic[class_index];
                    organic.doc_oxidation_g_c = organic.doc_oxidation_g_c + contribution.doc_oxidation_g_c;
                    organic.acetate_oxidation_g_c = organic.acetate_oxidation_g_c + contribution.acetate_oxidation_g_c;
                }
            }
        }
        next_demand.no2_g_n = next_demand.no2_g_n + correction.no2_g_n;
        next_demand.band_no2_g_n = next_demand.band_no2_g_n + correction.band_no2_g_n;
        if (!finiteStruct(next_demand)) return error.NonFiniteSoilCompetitionDemand;
        for (next_organic) |organic|
            if (!finiteStruct(organic)) return error.NonFiniteSoilCompetitionDemand;
        demand.* = next_demand;
        for (0..organic_class_count) |class_index|
            organic_by_layer_and_class[organic_start + class_index] = next_organic[class_index];
    }
}

test "REDIST soil competition aggregates all 42 populations and post-loop corrections" {
    var demand = [_]LayerDemand{std.mem.zeroes(LayerDemand)};
    var organic = [_]OrganicDemand{std.mem.zeroes(OrganicDemand)} ** 5;
    var populations = [_]PopulationDemand{std.mem.zeroes(PopulationDemand)} ** 42;
    for (&populations) |*population| {
        population.o2_g_o = 1.0;
        population.nh4_maximum_g_n = 2.0;
        population.band_hpo4_g_p = 3.0;
        population.doc_oxidation_g_c = 4.0;
    }
    const corrections = [_]LayerCorrections{.{ .no2_g_n = 5.0, .band_no2_g_n = 6.0 }};
    try aggregateLayers(&demand, &organic, &populations, &corrections);
    try std.testing.expectEqual(@as(f64, 42.0), demand[0].o2_g_o);
    try std.testing.expectEqual(@as(f64, 84.0), demand[0].nh4_g_n);
    try std.testing.expectEqual(@as(f64, 126.0), demand[0].band_hpo4_g_p);
    try std.testing.expectEqual(@as(f64, 5.0), demand[0].no2_g_n);
    try std.testing.expectEqual(@as(f64, 6.0), demand[0].band_no2_g_n);
    for (organic) |value| try std.testing.expectEqual(@as(f64, 28.0), value.doc_oxidation_g_c);
}

test "REDIST soil competition excludes K five from organic arrays" {
    var demand = [_]LayerDemand{std.mem.zeroes(LayerDemand)};
    var organic = [_]OrganicDemand{std.mem.zeroes(OrganicDemand)} ** 5;
    var populations = [_]PopulationDemand{std.mem.zeroes(PopulationDemand)} ** 42;
    populations[5 * 7].o2_g_o = 2.0;
    populations[5 * 7].doc_oxidation_g_c = 99.0;
    const corrections = [_]LayerCorrections{std.mem.zeroes(LayerCorrections)};
    try aggregateLayers(&demand, &organic, &populations, &corrections);
    try std.testing.expectEqual(@as(f64, 2.0), demand[0].o2_g_o);
    for (organic) |value| try std.testing.expectEqual(@as(f64, 0.0), value.doc_oxidation_g_c);
}

test "REDIST soil competition runtime layers remain independent" {
    var demands = [_]LayerDemand{ std.mem.zeroes(LayerDemand), std.mem.zeroes(LayerDemand) };
    var organic = [_]OrganicDemand{std.mem.zeroes(OrganicDemand)} ** 10;
    var populations = [_]PopulationDemand{std.mem.zeroes(PopulationDemand)} ** 84;
    populations[0].no3_maximum_g_n = 2.0;
    populations[42].no3_maximum_g_n = 3.0;
    const corrections = [_]LayerCorrections{ std.mem.zeroes(LayerCorrections), std.mem.zeroes(LayerCorrections) };
    try aggregateLayers(&demands, &organic, &populations, &corrections);
    try std.testing.expectEqual(@as(f64, 2.0), demands[0].no3_g_n);
    try std.testing.expectEqual(@as(f64, 3.0), demands[1].no3_g_n);
}

test "REDIST soil competition rejects dimensions invalid values and overflow" {
    var demand = [_]LayerDemand{std.mem.zeroes(LayerDemand)};
    var organic = [_]OrganicDemand{std.mem.zeroes(OrganicDemand)} ** 5;
    var populations = [_]PopulationDemand{std.mem.zeroes(PopulationDemand)} ** 42;
    const no_corrections: [0]LayerCorrections = .{};
    try std.testing.expectError(error.SoilCompetitionDemandDimensionMismatch, aggregateLayers(&demand, &organic, &populations, &no_corrections));
    const corrections = [_]LayerCorrections{std.mem.zeroes(LayerCorrections)};
    populations[0].n2o_g_n = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilCompetitionDemand, aggregateLayers(&demand, &organic, &populations, &corrections));
    populations[0].n2o_g_n = std.math.floatMax(f64);
    populations[1].n2o_g_n = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilCompetitionDemand, aggregateLayers(&demand, &organic, &populations, &corrections));
    try std.testing.expectEqual(std.mem.zeroes(LayerDemand), demand[0]);
    for (organic) |value| try std.testing.expectEqual(std.mem.zeroes(OrganicDemand), value);
}
