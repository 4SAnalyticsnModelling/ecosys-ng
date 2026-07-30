const std = @import("std");
const chemistry = @import("solute_chemistry_state.zig");
const transport = @import("solute_transport.zig");
const species_module = @import("solute_transport_species.zig");

const Species = species_module.AqueousSpecies;

/// Copies dissolved chemistry concentrations into the conservative transport
/// inventory. Chemistry uses mol m-3; transport uses mol per runtime layer.
pub fn exportChemistry(chemistry_state: *const chemistry.State, transport_state: *transport.State) !void {
    try validateDimensions(chemistry_state, transport_state);
    // Validate the complete transaction before publishing any amount.
    for (0..chemistry_state.cell_count) |cell| {
        const water_volume_m3 = transport_state.water_volume_m3[cell];
        if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0) return error.InvalidAqueousTransportWaterVolume;
        inline for (@typeInfo(Species).@"enum".fields) |field| {
            const dissolved_mol_per_m3 = concentration(chemistry_state, cell, @enumFromInt(field.value));
            if (!std.math.isFinite(dissolved_mol_per_m3) or dissolved_mol_per_m3 < 0) return error.InvalidAqueousChemistryConcentration;
            const amount_mol = dissolved_mol_per_m3 * water_volume_m3;
            if (!std.math.isFinite(amount_mol)) return error.InvalidAqueousTransportAmount;
        }
    }
    for (0..chemistry_state.cell_count) |cell| {
        const water_volume_m3 = transport_state.water_volume_m3[cell];
        inline for (@typeInfo(Species).@"enum".fields) |field| {
            transport_state.amount_mol[cell * Species.count + field.value] =
                concentration(chemistry_state, cell, @enumFromInt(field.value)) * water_volume_m3;
        }
    }
}

/// Captures STARTE-equilibrated concentrations used by TRNSFRS for later
/// external water-table recharge (`C*U`). Output is cell-major mol m-3.
pub fn exportConcentrations(chemistry_state: *const chemistry.State, output_mol_per_m3: []f64) !void {
    if (output_mol_per_m3.len != try std.math.mul(usize, chemistry_state.cell_count, Species.count)) return error.AqueousTransportDimensionMismatch;
    for (0..chemistry_state.cell_count) |cell| inline for (@typeInfo(Species).@"enum".fields) |field| {
        const value = concentration(chemistry_state, cell, @enumFromInt(field.value));
        if (!std.math.isFinite(value) or value < 0) return error.InvalidAqueousChemistryConcentration;
        output_mol_per_m3[cell * Species.count + field.value] = value;
    };
}

/// Publishes a transported micropore inventory back to the chemistry state.
/// The update is staged so zero water or a non-finite value cannot partially
/// overwrite the reaction state.
pub fn importChemistry(transport_state: *const transport.State, chemistry_state: *chemistry.State) !void {
    try validateDimensions(chemistry_state, transport_state);
    // Validate the complete transaction before changing any concentration.
    for (0..chemistry_state.cell_count) |cell| {
        const water_volume_m3 = transport_state.water_volume_m3[cell];
        if (!std.math.isFinite(water_volume_m3) or water_volume_m3 <= 0) return error.AqueousTransportRequiresPositiveWaterVolume;
        inline for (@typeInfo(Species).@"enum".fields) |field| {
            const amount_mol = transport_state.amount_mol[cell * Species.count + field.value];
            if (!std.math.isFinite(amount_mol) or amount_mol < 0) return error.InvalidAqueousTransportAmount;
        }
    }
    for (0..chemistry_state.cell_count) |cell| {
        const water_volume_m3 = transport_state.water_volume_m3[cell];
        inline for (@typeInfo(Species).@"enum".fields) |field| {
            const amount_mol = transport_state.amount_mol[cell * Species.count + field.value];
            setConcentration(&chemistry_state.aqueous[cell], &chemistry_state.non_band_phosphate[cell], &chemistry_state.band_phosphate[cell], @enumFromInt(field.value), amount_mol / water_volume_m3);
        }
    }
}

fn validateDimensions(chemistry_state: *const chemistry.State, transport_state: *const transport.State) !void {
    if (chemistry_state.cell_count != transport_state.cell_count or transport_state.species_count != Species.count) return error.AqueousTransportDimensionMismatch;
}

fn concentration(state: *const chemistry.State, cell: usize, species: Species) f64 {
    const aqueous = state.aqueous[cell];
    const non_band = state.non_band_phosphate[cell];
    const band = state.band_phosphate[cell];
    return switch (species) {
        .non_band_phosphate => non_band.dissolved_po4_mol_p_per_m3,
        .non_band_phosphoric_acid => non_band.dissolved_h3po4_mol_p_per_m3,
        .non_band_iron_hpo4 => non_band.iron_hpo4_pair_mol_per_m3,
        .non_band_iron_h2po4 => non_band.iron_h2po4_pair_mol_per_m3,
        .non_band_calcium_phosphate => non_band.calcium_po4_pair_mol_per_m3,
        .non_band_calcium_hpo4 => non_band.calcium_hpo4_pair_mol_per_m3,
        .non_band_calcium_h2po4 => non_band.calcium_h2po4_pair_mol_per_m3,
        .non_band_magnesium_hpo4 => non_band.magnesium_hpo4_pair_mol_per_m3,
        .band_phosphate => band.dissolved_po4_mol_p_per_m3,
        .band_phosphoric_acid => band.dissolved_h3po4_mol_p_per_m3,
        .band_iron_hpo4 => band.iron_hpo4_pair_mol_per_m3,
        .band_iron_h2po4 => band.iron_h2po4_pair_mol_per_m3,
        .band_calcium_phosphate => band.calcium_po4_pair_mol_per_m3,
        .band_calcium_hpo4 => band.calcium_hpo4_pair_mol_per_m3,
        .band_calcium_h2po4 => band.calcium_h2po4_pair_mol_per_m3,
        .band_magnesium_hpo4 => band.magnesium_hpo4_pair_mol_per_m3,
        .aluminum => aqueous.aluminum,
        .iron => aqueous.iron,
        .hydrogen => aqueous.hydrogen,
        .calcium => aqueous.calcium,
        .magnesium => aqueous.magnesium,
        .sodium => aqueous.sodium,
        .potassium => aqueous.potassium,
        .hydroxide => aqueous.hydroxide,
        .sulfate => aqueous.sulfate,
        .chloride => aqueous.chloride,
        .carbonate => aqueous.carbonate,
        .bicarbonate => aqueous.bicarbonate,
        .aluminum_hydroxide_1 => aqueous.aluminum_hydroxide_1,
        .aluminum_hydroxide_2 => aqueous.aluminum_hydroxide_2,
        .aluminum_hydroxide_3 => aqueous.aluminum_hydroxide_3,
        .aluminum_hydroxide_4 => aqueous.aluminum_hydroxide_4,
        .aluminum_sulfate => aqueous.aluminum_sulfate,
        .iron_hydroxide_1 => aqueous.iron_hydroxide_1,
        .iron_hydroxide_2 => aqueous.iron_hydroxide_2,
        .iron_hydroxide_3 => aqueous.iron_hydroxide_3,
        .iron_hydroxide_4 => aqueous.iron_hydroxide_4,
        .iron_sulfate => aqueous.iron_sulfate,
        .calcium_hydroxide => aqueous.calcium_hydroxide,
        .calcium_carbonate => aqueous.calcium_carbonate,
        .calcium_bicarbonate => aqueous.calcium_bicarbonate,
        .calcium_sulfate => aqueous.calcium_sulfate,
        .magnesium_hydroxide => aqueous.magnesium_hydroxide,
        .magnesium_carbonate => aqueous.magnesium_carbonate,
        .magnesium_bicarbonate => aqueous.magnesium_bicarbonate,
        .magnesium_sulfate => aqueous.magnesium_sulfate,
        .sodium_carbonate => aqueous.sodium_carbonate,
        .sodium_sulfate => aqueous.sodium_sulfate,
        .potassium_sulfate => aqueous.potassium_sulfate,
        .hydrogen_silicate => aqueous.hydrogen_silicate,
    };
}

fn setConcentration(aqueous: anytype, non_band: anytype, band: anytype, species: Species, value: f64) void {
    switch (species) {
        .non_band_phosphate => non_band.dissolved_po4_mol_p_per_m3 = value,
        .non_band_phosphoric_acid => non_band.dissolved_h3po4_mol_p_per_m3 = value,
        .non_band_iron_hpo4 => non_band.iron_hpo4_pair_mol_per_m3 = value,
        .non_band_iron_h2po4 => non_band.iron_h2po4_pair_mol_per_m3 = value,
        .non_band_calcium_phosphate => non_band.calcium_po4_pair_mol_per_m3 = value,
        .non_band_calcium_hpo4 => non_band.calcium_hpo4_pair_mol_per_m3 = value,
        .non_band_calcium_h2po4 => non_band.calcium_h2po4_pair_mol_per_m3 = value,
        .non_band_magnesium_hpo4 => non_band.magnesium_hpo4_pair_mol_per_m3 = value,
        .band_phosphate => band.dissolved_po4_mol_p_per_m3 = value,
        .band_phosphoric_acid => band.dissolved_h3po4_mol_p_per_m3 = value,
        .band_iron_hpo4 => band.iron_hpo4_pair_mol_per_m3 = value,
        .band_iron_h2po4 => band.iron_h2po4_pair_mol_per_m3 = value,
        .band_calcium_phosphate => band.calcium_po4_pair_mol_per_m3 = value,
        .band_calcium_hpo4 => band.calcium_hpo4_pair_mol_per_m3 = value,
        .band_calcium_h2po4 => band.calcium_h2po4_pair_mol_per_m3 = value,
        .band_magnesium_hpo4 => band.magnesium_hpo4_pair_mol_per_m3 = value,
        .aluminum => aqueous.aluminum = value,
        .iron => aqueous.iron = value,
        .hydrogen => aqueous.hydrogen = value,
        .calcium => aqueous.calcium = value,
        .magnesium => aqueous.magnesium = value,
        .sodium => aqueous.sodium = value,
        .potassium => aqueous.potassium = value,
        .hydroxide => aqueous.hydroxide = value,
        .sulfate => aqueous.sulfate = value,
        .chloride => aqueous.chloride = value,
        .carbonate => aqueous.carbonate = value,
        .bicarbonate => aqueous.bicarbonate = value,
        .aluminum_hydroxide_1 => aqueous.aluminum_hydroxide_1 = value,
        .aluminum_hydroxide_2 => aqueous.aluminum_hydroxide_2 = value,
        .aluminum_hydroxide_3 => aqueous.aluminum_hydroxide_3 = value,
        .aluminum_hydroxide_4 => aqueous.aluminum_hydroxide_4 = value,
        .aluminum_sulfate => aqueous.aluminum_sulfate = value,
        .iron_hydroxide_1 => aqueous.iron_hydroxide_1 = value,
        .iron_hydroxide_2 => aqueous.iron_hydroxide_2 = value,
        .iron_hydroxide_3 => aqueous.iron_hydroxide_3 = value,
        .iron_hydroxide_4 => aqueous.iron_hydroxide_4 = value,
        .iron_sulfate => aqueous.iron_sulfate = value,
        .calcium_hydroxide => aqueous.calcium_hydroxide = value,
        .calcium_carbonate => aqueous.calcium_carbonate = value,
        .calcium_bicarbonate => aqueous.calcium_bicarbonate = value,
        .calcium_sulfate => aqueous.calcium_sulfate = value,
        .magnesium_hydroxide => aqueous.magnesium_hydroxide = value,
        .magnesium_carbonate => aqueous.magnesium_carbonate = value,
        .magnesium_bicarbonate => aqueous.magnesium_bicarbonate = value,
        .magnesium_sulfate => aqueous.magnesium_sulfate = value,
        .sodium_carbonate => aqueous.sodium_carbonate = value,
        .sodium_sulfate => aqueous.sodium_sulfate = value,
        .potassium_sulfate => aqueous.potassium_sulfate = value,
        .hydrogen_silicate => aqueous.hydrogen_silicate = value,
    }
}

test "all TRNSFRS species round trip between concentration and runtime amount" {
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 2);
    defer chemistry_state.deinit();
    var transport_state = try transport.State.init(std.testing.allocator, 2, Species.count);
    defer transport_state.deinit();
    transport_state.water_volume_m3[0] = 2;
    transport_state.water_volume_m3[1] = 0.5;
    inline for (@typeInfo(Species).@"enum".fields) |field| setConcentration(&chemistry_state.aqueous[0], &chemistry_state.non_band_phosphate[0], &chemistry_state.band_phosphate[0], @enumFromInt(field.value), @as(f64, @floatFromInt(field.value + 1)));
    try exportChemistry(&chemistry_state, &transport_state);
    inline for (@typeInfo(Species).@"enum".fields) |field| try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (field.value + 1))), transport_state.amount_mol[field.value]);
    @memcpy(transport_state.amount_mol[Species.count .. 2 * Species.count], transport_state.amount_mol[0..Species.count]);
    try importChemistry(&transport_state, &chemistry_state);
    inline for (@typeInfo(Species).@"enum".fields) |field| try std.testing.expectEqual(@as(f64, @floatFromInt(4 * (field.value + 1))), concentration(&chemistry_state, 1, @enumFromInt(field.value)));
}

test "failed import cannot partially overwrite chemistry" {
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    chemistry_state.aqueous[0].calcium = 3;
    var transport_state = try transport.State.init(std.testing.allocator, 1, Species.count);
    defer transport_state.deinit();
    transport_state.water_volume_m3[0] = 0;
    try std.testing.expectError(error.AqueousTransportRequiresPositiveWaterVolume, importChemistry(&transport_state, &chemistry_state));
    try std.testing.expectEqual(@as(f64, 3), chemistry_state.aqueous[0].calcium);
}
