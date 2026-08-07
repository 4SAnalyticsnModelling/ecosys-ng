const std = @import("std");
const snow = @import("../solute/snow_solute_transport.zig");
const gas = @import("../gas/transport.zig");
const litter = @import("../../surface/litter_chemistry.zig");
const chemistry = @import("../solute/chemistry_state.zig");

pub const Inputs = struct {
    discharge: []const snow.SurfaceDischarge,
    litter_water_volume_m3: []const f64,
    topsoil_water_volume_m3: []const f64,
    soil_layer_capacity: usize,
    ion_molar_mass_g_per_mol: IonMolarMassesGPerMol,
};

pub const IonMolarMassesGPerMol = struct {
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
    sulfur: f64,
    chloride: f64,
};

/// Conservatively receives TRNSFR snow discharge into the authoritative
/// litter and topsoil gas/mineral pools. All destinations publish together.
pub fn commit(allocator: std.mem.Allocator, inputs: Inputs, litter_gas: *gas.State, soil_gas: *gas.State, litter_chemistry: *litter.State, soil_chemistry: *chemistry.State) !void {
    const cells = inputs.discharge.len;
    if (cells == 0 or inputs.litter_water_volume_m3.len != cells or inputs.soil_layer_capacity == 0 or litter_gas.cell_count != cells or litter_chemistry.cells.len != cells or soil_gas.cell_count != soil_chemistry.cell_count or soil_chemistry.cell_count != cells * inputs.soil_layer_capacity or inputs.topsoil_water_volume_m3.len != soil_chemistry.cell_count) return error.SnowDischargeDimensionMismatch;
    const staged_litter_gas = try allocator.dupe(f64, litter_gas.dissolved_mass_g);
    defer allocator.free(staged_litter_gas);
    const staged_soil_gas = try allocator.dupe(f64, soil_gas.dissolved_mass_g);
    defer allocator.free(staged_soil_gas);
    const staged_litter = try allocator.dupe(litter.Cell, litter_chemistry.cells);
    defer allocator.free(staged_litter);
    const staged_aqueous = try allocator.dupe(@TypeOf(soil_chemistry.aqueous[0]), soil_chemistry.aqueous);
    defer allocator.free(staged_aqueous);
    const staged_non_band = try allocator.dupe(@TypeOf(soil_chemistry.non_band_phosphate[0]), soil_chemistry.non_band_phosphate);
    defer allocator.free(staged_non_band);
    const staged_band = try allocator.dupe(@TypeOf(soil_chemistry.band_phosphate[0]), soil_chemistry.band_phosphate);
    defer allocator.free(staged_band);

    for (0..cells) |cell| {
        const litter_water = inputs.litter_water_volume_m3[cell];
        const topsoil = cell * inputs.soil_layer_capacity;
        const soil_water = inputs.topsoil_water_volume_m3[topsoil];
        if (!std.math.isFinite(litter_water) or litter_water < 0 or !std.math.isFinite(soil_water) or soil_water < 0) return error.InvalidSnowDischargeWaterVolume;
        for (0..snow.species_count) |species| {
            const litter_g = inputs.discharge[cell].litter_g[species];
            const non_band_g = inputs.discharge[cell].soil_nonband_g[species];
            const band_g = inputs.discharge[cell].soil_band_g[species];
            inline for (.{ litter_g, non_band_g, band_g }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowSurfaceDischarge;
            if (species < 5) {
                staged_litter_gas[cell * gas.species_count + species] += litter_g;
                staged_soil_gas[topsoil * gas.species_count + species] += non_band_g + band_g;
                continue;
            }
            if (litter_g > 0 and litter_water <= 0 or non_band_g + band_g > 0 and soil_water <= 0) return error.SnowDischargeRequiresWater;
            const molar_mass = molarMass(inputs.ion_molar_mass_g_per_mol, species);
            if (!std.math.isFinite(molar_mass) or molar_mass <= 0) return error.InvalidSnowDischargeMolarMass;
            const litter_mol_per_m3 = if (litter_water > 0) litter_g / molar_mass / litter_water else 0;
            const non_band_mol_per_m3 = if (soil_water > 0) non_band_g / molar_mass / soil_water else 0;
            const band_mol_per_m3 = if (soil_water > 0) band_g / molar_mass / soil_water else 0;
            switch (@as(snow.Species, @enumFromInt(species))) {
                .ammonium_nitrogen => {
                    staged_litter[cell].ammonium_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].ammonium_non_band += non_band_mol_per_m3;
                    staged_aqueous[topsoil].ammonium_band += band_mol_per_m3;
                },
                .ammonia_nitrogen => {
                    staged_litter[cell].ammonia_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].ammonia_non_band += non_band_mol_per_m3;
                    staged_aqueous[topsoil].ammonia_band += band_mol_per_m3;
                },
                .nitrate_nitrogen => {
                    staged_litter[cell].nitrate_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].nitrate_non_band += non_band_mol_per_m3;
                    staged_aqueous[topsoil].nitrate_band += band_mol_per_m3;
                },
                .hydrogen_phosphate_phosphorus => {
                    staged_litter[cell].hpo4_mol_p_per_m3 += litter_mol_per_m3;
                    staged_non_band[topsoil].dissolved_hpo4_mol_p_per_m3 += non_band_mol_per_m3;
                    staged_band[topsoil].dissolved_hpo4_mol_p_per_m3 += band_mol_per_m3;
                },
                .dihydrogen_phosphate_phosphorus => {
                    staged_litter[cell].h2po4_mol_p_per_m3 += litter_mol_per_m3;
                    staged_non_band[topsoil].dissolved_h2po4_mol_p_per_m3 += non_band_mol_per_m3;
                    staged_band[topsoil].dissolved_h2po4_mol_p_per_m3 += band_mol_per_m3;
                },
                .aluminum => {
                    staged_litter[cell].aluminum_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].aluminum += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .iron => {
                    staged_litter[cell].iron_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].iron += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .calcium => {
                    staged_litter[cell].calcium_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].calcium += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .magnesium => {
                    staged_litter[cell].magnesium_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].magnesium += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .sodium => {
                    staged_litter[cell].sodium_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].sodium += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .potassium => {
                    staged_litter[cell].potassium_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].potassium += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .sulfate_sulfur => {
                    staged_litter[cell].sulfate_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].sulfate += non_band_mol_per_m3 + band_mol_per_m3;
                },
                .chloride => {
                    staged_litter[cell].chloride_mol_per_m3 += litter_mol_per_m3;
                    staged_aqueous[topsoil].chloride += non_band_mol_per_m3 + band_mol_per_m3;
                },
                else => unreachable,
            }
        }
    }
    inline for (.{ staged_litter_gas, staged_soil_gas }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowDischargeCandidate;
    @memcpy(litter_gas.dissolved_mass_g, staged_litter_gas);
    @memcpy(soil_gas.dissolved_mass_g, staged_soil_gas);
    @memcpy(litter_chemistry.cells, staged_litter);
    @memcpy(soil_chemistry.aqueous, staged_aqueous);
    @memcpy(soil_chemistry.non_band_phosphate, staged_non_band);
    @memcpy(soil_chemistry.band_phosphate, staged_band);
}

fn molarMass(ions: IonMolarMassesGPerMol, species: usize) f64 {
    return switch (@as(snow.Species, @enumFromInt(species))) {
        .ammonium_nitrogen, .ammonia_nitrogen, .nitrate_nitrogen => snow.nitrogen_g_per_mol,
        .hydrogen_phosphate_phosphorus, .dihydrogen_phosphate_phosphorus => snow.phosphorus_g_per_mol,
        .aluminum => ions.aluminum,
        .iron => ions.iron,
        .calcium => ions.calcium,
        .magnesium => ions.magnesium,
        .sodium => ions.sodium,
        .potassium => ions.potassium,
        .sulfate_sulfur => ions.sulfur,
        .chloride => ions.chloride,
        else => unreachable,
    };
}

test "snow discharge conserves tracked gas nitrogen and phosphorus into runtime recipients" {
    var litter_gas = try gas.State.init(std.testing.allocator, 1);
    defer litter_gas.deinit();
    var soil_gas = try gas.State.init(std.testing.allocator, 2);
    defer soil_gas.deinit();
    var litter_chemistry = try litter.State.init(std.testing.allocator, 1);
    defer litter_chemistry.deinit();
    var soil_chemistry = try chemistry.State.init(std.testing.allocator, 2);
    defer soil_chemistry.deinit();
    var discharge = [_]snow.SurfaceDischarge{.{}};
    discharge[0].litter_g[@intFromEnum(snow.Species.carbon_dioxide_carbon)] = 2;
    discharge[0].soil_nonband_g[@intFromEnum(snow.Species.ammonium_nitrogen)] = 14;
    discharge[0].soil_band_g[@intFromEnum(snow.Species.hydrogen_phosphate_phosphorus)] = 31;
    discharge[0].soil_nonband_g[@intFromEnum(snow.Species.calcium)] = 40;
    try commit(std.testing.allocator, .{ .discharge = &discharge, .litter_water_volume_m3 = &.{1}, .topsoil_water_volume_m3 = &.{ 2, 0 }, .soil_layer_capacity = 2, .ion_molar_mass_g_per_mol = .{ .aluminum = 27, .iron = 55.8, .calcium = 40, .magnesium = 24.3, .sodium = 23, .potassium = 39.1, .sulfur = 32, .chloride = 35.5 } }, &litter_gas, &soil_gas, &litter_chemistry, &soil_chemistry);
    try std.testing.expectEqual(@as(f64, 2), litter_gas.dissolved_mass_g[0]);
    try std.testing.expectEqual(@as(f64, 0.5), soil_chemistry.aqueous[0].ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 0.5), soil_chemistry.band_phosphate[0].dissolved_hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 0.5), soil_chemistry.aqueous[0].calcium);
}
