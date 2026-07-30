const std = @import("std");
const aqueous_network = @import("solute_aqueous_network.zig");

const minimum_concentration_mol_per_m3 = 1.0e-32;

pub const NitrateZoneInput = struct {
    nitrate_g_n: f64,
    water_volume_m3: f64,
};

pub const NitrateConcentrations = struct {
    non_band_mol_n_per_m3: f64,
    band_mol_n_per_m3: f64,
};

/// SOLUTE lines 526--542. The source stores nitrate as g N and divides by
/// 14 g N mol-1 before converting each runtime water-volume zone.
pub fn initializeNitrateConcentrations(
    non_band: NitrateZoneInput,
    band: NitrateZoneInput,
    minimum_water_volume_m3: f64,
) !NitrateConcentrations {
    try validateFiniteNonnegative(minimum_water_volume_m3);
    return .{
        .non_band_mol_n_per_m3 = try nitrateConcentration(
            non_band,
            minimum_water_volume_m3,
        ),
        .band_mol_n_per_m3 = try nitrateConcentration(
            band,
            minimum_water_volume_m3,
        ),
    };
}

fn nitrateConcentration(
    input: NitrateZoneInput,
    minimum_water_volume_m3: f64,
) !f64 {
    if (!std.math.isFinite(input.nitrate_g_n) or
        !std.math.isFinite(input.water_volume_m3) or
        input.water_volume_m3 < 0)
        return error.InvalidNitrateConcentrationInput;
    if (input.water_volume_m3 <= minimum_water_volume_m3) return 0;
    return @max(
        minimum_concentration_mol_per_m3,
        input.nitrate_g_n / (14.0 * input.water_volume_m3),
    );
}

pub const CationExchangeAmounts = struct {
    capacity_mol_charge: f64,
    hydrogen_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    carboxyl_bound_hydrogen_mol: f64,
};

pub const CationExchangeConcentrations = struct {
    capacity_mol_charge_per_Mg: f64,
    hydrogen_mol_per_Mg: f64,
    aluminum_mol_per_Mg: f64,
    iron_mol_per_Mg: f64,
    calcium_mol_per_Mg: f64,
    magnesium_mol_per_Mg: f64,
    sodium_mol_per_Mg: f64,
    potassium_mol_per_Mg: f64,
    carboxyl_bound_hydrogen_mol_per_Mg: f64,
    total_carboxyl_sites_mol_per_Mg: f64,
    soil_mass_per_water_volume_Mg_per_m3: f64,
};

/// SOLUTE lines 544--578. `organic_carbon_g_c` is converted to Mg C by
/// 1e-6 before applying the carboxyl-site density in mol Mg-C-1.
pub fn initializeCationExchangeConcentrations(
    amounts: CationExchangeAmounts,
    soil_mass_Mg: f64,
    water_volume_m3: f64,
    carboxyl_sites_mol_per_Mg_c: f64,
    organic_carbon_g_c: f64,
    minimum_active_extent: f64,
) !CationExchangeConcentrations {
    inline for (@typeInfo(CationExchangeAmounts).@"struct".fields) |field|
        if (!std.math.isFinite(@field(amounts, field.name)))
            return error.InvalidCationExchangeInitialization;
    inline for (.{
        soil_mass_Mg,
        water_volume_m3,
        carboxyl_sites_mol_per_Mg_c,
        organic_carbon_g_c,
        minimum_active_extent,
    }) |value| try validateFiniteNonnegative(value);

    // Legacy ZEROS(NY,NX) is one cell-specific numerical threshold applied
    // contextually to both the Mg soil mass and m3 water volume.
    if (soil_mass_Mg <= minimum_active_extent or
        water_volume_m3 <= minimum_active_extent)
        return zeroCationExchangeConcentrations();

    return .{
        .capacity_mol_charge_per_Mg = concentration(amounts.capacity_mol_charge, soil_mass_Mg),
        .hydrogen_mol_per_Mg = concentration(amounts.hydrogen_mol, soil_mass_Mg),
        .aluminum_mol_per_Mg = concentration(amounts.aluminum_mol, soil_mass_Mg),
        .iron_mol_per_Mg = concentration(amounts.iron_mol, soil_mass_Mg),
        .calcium_mol_per_Mg = concentration(amounts.calcium_mol, soil_mass_Mg),
        .magnesium_mol_per_Mg = concentration(amounts.magnesium_mol, soil_mass_Mg),
        .sodium_mol_per_Mg = concentration(amounts.sodium_mol, soil_mass_Mg),
        .potassium_mol_per_Mg = concentration(amounts.potassium_mol, soil_mass_Mg),
        .carboxyl_bound_hydrogen_mol_per_Mg = concentration(
            amounts.carboxyl_bound_hydrogen_mol,
            soil_mass_Mg,
        ),
        .total_carboxyl_sites_mol_per_Mg = concentration(
            carboxyl_sites_mol_per_Mg_c * 1.0e-6 * organic_carbon_g_c,
            soil_mass_Mg,
        ),
        .soil_mass_per_water_volume_Mg_per_m3 = soil_mass_Mg / water_volume_m3,
    };
}

fn concentration(amount_mol: f64, soil_mass_Mg: f64) f64 {
    return @max(minimum_concentration_mol_per_m3, amount_mol / soil_mass_Mg);
}

fn zeroCationExchangeConcentrations() CationExchangeConcentrations {
    return .{
        .capacity_mol_charge_per_Mg = 0,
        .hydrogen_mol_per_Mg = 0,
        .aluminum_mol_per_Mg = 0,
        .iron_mol_per_Mg = 0,
        .calcium_mol_per_Mg = 0,
        .magnesium_mol_per_Mg = 0,
        .sodium_mol_per_Mg = 0,
        .potassium_mol_per_Mg = 0,
        .carboxyl_bound_hydrogen_mol_per_Mg = 0,
        .total_carboxyl_sites_mol_per_Mg = 0,
        .soil_mass_per_water_volume_Mg_per_m3 = 1,
    };
}

fn validateFiniteNonnegative(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSoluteConcentrationInitialization;
}

pub const FreeIonAmounts = struct {
    hydrogen_mol: f64,
    hydroxide_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
    carbonate_mol: f64,
    bicarbonate_mol: f64,
    carbon_dioxide_g_c: f64,
};

pub const DissolvedComplexAmounts = struct {
    aluminum_hydroxide_1_mol: f64,
    aluminum_hydroxide_2_mol: f64,
    aluminum_hydroxide_3_mol: f64,
    aluminum_hydroxide_4_mol: f64,
    aluminum_sulfate_mol: f64,
    iron_hydroxide_1_mol: f64,
    iron_hydroxide_2_mol: f64,
    iron_hydroxide_3_mol: f64,
    iron_hydroxide_4_mol: f64,
    iron_sulfate_mol: f64,
    calcium_hydroxide_mol: f64,
    calcium_carbonate_mol: f64,
    calcium_bicarbonate_mol: f64,
    calcium_sulfate_mol: f64,
    magnesium_hydroxide_mol: f64,
    magnesium_carbonate_mol: f64,
    magnesium_bicarbonate_mol: f64,
    magnesium_sulfate_mol: f64,
    sodium_carbonate_mol: f64,
    sodium_sulfate_mol: f64,
    potassium_sulfate_mol: f64,
    hydrogen_silicate_mol: f64,
};

pub const DynamicSaltAqueousInputs = struct {
    free_ions: FreeIonAmounts,
    dissolved_complexes: DissolvedComplexAmounts,
    water_volume_m3: f64,
    minimum_water_volume_m3: f64,
    dry_soil_ph: f64,
    water_activity_product_mol2_per_m6: f64,
};

/// Direct translation of SOLUTE lines 646--719. Nutrient-zone ammonium,
/// ammonia, and nitrate fields are separate owners and remain unchanged.
pub fn applyDynamicSaltAqueousConcentrations(
    state: *aqueous_network.State,
    inputs: DynamicSaltAqueousInputs,
) !void {
    try validateDynamicSaltAqueousInputs(inputs);
    var staged = state.*;
    if (inputs.water_volume_m3 > inputs.minimum_water_volume_m3) {
        const volume = inputs.water_volume_m3;
        const free = inputs.free_ions;
        staged.hydrogen = amountPerWaterVolume(free.hydrogen_mol, volume);
        staged.hydroxide = amountPerWaterVolume(free.hydroxide_mol, volume);
        staged.aluminum = amountPerWaterVolume(free.aluminum_mol, volume);
        staged.iron = amountPerWaterVolume(free.iron_mol, volume);
        staged.calcium = amountPerWaterVolume(free.calcium_mol, volume);
        staged.magnesium = amountPerWaterVolume(free.magnesium_mol, volume);
        staged.sodium = amountPerWaterVolume(free.sodium_mol, volume);
        staged.potassium = amountPerWaterVolume(free.potassium_mol, volume);
        staged.sulfate = amountPerWaterVolume(free.sulfate_mol, volume);
        staged.chloride = amountPerWaterVolume(free.chloride_mol, volume);
        staged.carbonate = amountPerWaterVolume(free.carbonate_mol, volume);
        staged.bicarbonate = amountPerWaterVolume(free.bicarbonate_mol, volume);
        staged.carbon_dioxide = amountPerWaterVolume(
            free.carbon_dioxide_g_c / 12.0,
            volume,
        );

        const complexes = inputs.dissolved_complexes;
        staged.aluminum_hydroxide_1 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_1_mol,
            volume,
        );
        staged.aluminum_hydroxide_2 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_2_mol,
            volume,
        );
        staged.aluminum_hydroxide_3 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_3_mol,
            volume,
        );
        staged.aluminum_hydroxide_4 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_4_mol,
            volume,
        );
        staged.aluminum_sulfate = amountPerWaterVolume(
            complexes.aluminum_sulfate_mol,
            volume,
        );
        staged.iron_hydroxide_1 = amountPerWaterVolume(
            complexes.iron_hydroxide_1_mol,
            volume,
        );
        staged.iron_hydroxide_2 = amountPerWaterVolume(
            complexes.iron_hydroxide_2_mol,
            volume,
        );
        staged.iron_hydroxide_3 = amountPerWaterVolume(
            complexes.iron_hydroxide_3_mol,
            volume,
        );
        staged.iron_hydroxide_4 = amountPerWaterVolume(
            complexes.iron_hydroxide_4_mol,
            volume,
        );
        staged.iron_sulfate = amountPerWaterVolume(
            complexes.iron_sulfate_mol,
            volume,
        );
        staged.calcium_hydroxide = amountPerWaterVolume(
            complexes.calcium_hydroxide_mol,
            volume,
        );
        staged.calcium_carbonate = amountPerWaterVolume(
            complexes.calcium_carbonate_mol,
            volume,
        );
        staged.calcium_bicarbonate = amountPerWaterVolume(
            complexes.calcium_bicarbonate_mol,
            volume,
        );
        staged.calcium_sulfate = amountPerWaterVolume(
            complexes.calcium_sulfate_mol,
            volume,
        );
        staged.magnesium_hydroxide = amountPerWaterVolume(
            complexes.magnesium_hydroxide_mol,
            volume,
        );
        staged.magnesium_carbonate = amountPerWaterVolume(
            complexes.magnesium_carbonate_mol,
            volume,
        );
        staged.magnesium_bicarbonate = amountPerWaterVolume(
            complexes.magnesium_bicarbonate_mol,
            volume,
        );
        staged.magnesium_sulfate = amountPerWaterVolume(
            complexes.magnesium_sulfate_mol,
            volume,
        );
        staged.sodium_carbonate = amountPerWaterVolume(
            complexes.sodium_carbonate_mol,
            volume,
        );
        staged.sodium_sulfate = amountPerWaterVolume(
            complexes.sodium_sulfate_mol,
            volume,
        );
        staged.potassium_sulfate = amountPerWaterVolume(
            complexes.potassium_sulfate_mol,
            volume,
        );
        staged.hydrogen_silicate = amountPerWaterVolume(
            complexes.hydrogen_silicate_mol,
            volume,
        );
    } else {
        staged.hydrogen =
            std.math.pow(f64, 10, -inputs.dry_soil_ph) * 1.0e3;
        staged.hydroxide =
            inputs.water_activity_product_mol2_per_m6 / staged.hydrogen;
        inline for (.{
            "aluminum",
            "iron",
            "calcium",
            "magnesium",
            "sodium",
            "potassium",
            "sulfate",
            "chloride",
            "carbonate",
            "bicarbonate",
            "carbon_dioxide",
            "aluminum_hydroxide_1",
            "aluminum_hydroxide_2",
            "aluminum_hydroxide_3",
            "aluminum_hydroxide_4",
            "aluminum_sulfate",
            "iron_hydroxide_1",
            "iron_hydroxide_2",
            "iron_hydroxide_3",
            "iron_hydroxide_4",
            "iron_sulfate",
            "calcium_hydroxide",
            "calcium_carbonate",
            "calcium_bicarbonate",
            "calcium_sulfate",
            "magnesium_hydroxide",
            "magnesium_carbonate",
            "magnesium_bicarbonate",
            "magnesium_sulfate",
            "sodium_carbonate",
            "sodium_sulfate",
            "potassium_sulfate",
            "hydrogen_silicate",
        }) |field_name| @field(staged, field_name) = 0;
    }
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        const value = @field(staged, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDynamicSaltAqueousConcentration;
    }
    state.* = staged;
}

fn validateDynamicSaltAqueousInputs(inputs: DynamicSaltAqueousInputs) !void {
    inline for (@typeInfo(FreeIonAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.free_ions, field.name)))
            return error.InvalidDynamicSaltAqueousInput;
    }
    inline for (@typeInfo(DissolvedComplexAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.dissolved_complexes, field.name)))
            return error.InvalidDynamicSaltAqueousInput;
    }
    if (!std.math.isFinite(inputs.water_volume_m3) or
        inputs.water_volume_m3 < 0 or
        !std.math.isFinite(inputs.minimum_water_volume_m3) or
        inputs.minimum_water_volume_m3 < 0 or
        !std.math.isFinite(inputs.dry_soil_ph) or
        !std.math.isFinite(inputs.water_activity_product_mol2_per_m6) or
        inputs.water_activity_product_mol2_per_m6 <= 0)
    {
        return error.InvalidDynamicSaltAqueousInput;
    }
}

fn amountPerWaterVolume(amount_mol: f64, water_volume_m3: f64) f64 {
    return @max(
        minimum_concentration_mol_per_m3,
        amount_mol / water_volume_m3,
    );
}

fn dynamicSaltInputs() DynamicSaltAqueousInputs {
    return .{
        .free_ions = .{
            .hydrogen_mol = 2,
            .hydroxide_mol = 4,
            .aluminum_mol = 6,
            .iron_mol = 8,
            .calcium_mol = 10,
            .magnesium_mol = 12,
            .sodium_mol = 14,
            .potassium_mol = 16,
            .sulfate_mol = 18,
            .chloride_mol = 20,
            .carbonate_mol = 22,
            .bicarbonate_mol = 24,
            .carbon_dioxide_g_c = 48,
        },
        .dissolved_complexes = .{
            .aluminum_hydroxide_1_mol = 2,
            .aluminum_hydroxide_2_mol = 4,
            .aluminum_hydroxide_3_mol = 6,
            .aluminum_hydroxide_4_mol = 8,
            .aluminum_sulfate_mol = 10,
            .iron_hydroxide_1_mol = 12,
            .iron_hydroxide_2_mol = 14,
            .iron_hydroxide_3_mol = 16,
            .iron_hydroxide_4_mol = 18,
            .iron_sulfate_mol = 20,
            .calcium_hydroxide_mol = 22,
            .calcium_carbonate_mol = 24,
            .calcium_bicarbonate_mol = 26,
            .calcium_sulfate_mol = 28,
            .magnesium_hydroxide_mol = 30,
            .magnesium_carbonate_mol = 32,
            .magnesium_bicarbonate_mol = 34,
            .magnesium_sulfate_mol = 36,
            .sodium_carbonate_mol = 38,
            .sodium_sulfate_mol = 40,
            .potassium_sulfate_mol = 42,
            .hydrogen_silicate_mol = 44,
        },
        .water_volume_m3 = 2,
        .minimum_water_volume_m3 = 1.0e-12,
        .dry_soil_ph = 7,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
    };
}

test "SOLUTE nitrate concentrations preserve source conversion and dry branch" {
    const result = try initializeNitrateConcentrations(
        .{ .nitrate_g_n = 28, .water_volume_m3 = 2 },
        .{ .nitrate_g_n = -7, .water_volume_m3 = 1 },
        0,
    );
    try std.testing.expectEqual(@as(f64, 1), result.non_band_mol_n_per_m3);
    try std.testing.expectEqual(
        minimum_concentration_mol_per_m3,
        result.band_mol_n_per_m3,
    );
    const dry = try initializeNitrateConcentrations(
        .{ .nitrate_g_n = 28, .water_volume_m3 = 0 },
        .{ .nitrate_g_n = 28, .water_volume_m3 = 1.0e-9 },
        1.0e-9,
    );
    try std.testing.expectEqual(@as(f64, 0), dry.non_band_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 0), dry.band_mol_n_per_m3);
}

test "SOLUTE cation exchange setup preserves source order and units" {
    const result = try initializeCationExchangeConcentrations(
        .{
            .capacity_mol_charge = 20,
            .hydrogen_mol = 2,
            .aluminum_mol = 4,
            .iron_mol = 6,
            .calcium_mol = 8,
            .magnesium_mol = 10,
            .sodium_mol = 12,
            .potassium_mol = 14,
            .carboxyl_bound_hydrogen_mol = 16,
        },
        2,
        0.5,
        5,
        4.0e8,
        0,
    );
    try std.testing.expectEqual(@as(f64, 10), result.capacity_mol_charge_per_Mg);
    try std.testing.expectEqual(@as(f64, 1), result.hydrogen_mol_per_Mg);
    try std.testing.expectEqual(@as(f64, 8), result.carboxyl_bound_hydrogen_mol_per_Mg);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1000),
        result.total_carboxyl_sites_mol_per_Mg,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 4), result.soil_mass_per_water_volume_Mg_per_m3);
}

test "SOLUTE cation exchange dry branch returns source defaults" {
    const zero = std.mem.zeroes(CationExchangeAmounts);
    const result = try initializeCationExchangeConcentrations(
        zero,
        0,
        0,
        0,
        0,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), result.capacity_mol_charge_per_Mg);
    try std.testing.expectEqual(@as(f64, 1), result.soil_mass_per_water_volume_Mg_per_m3);
}

test "SOLUTE wet dynamic salt initialization preserves units and nutrient owners" {
    var state = std.mem.zeroes(aqueous_network.State);
    state.ammonium_non_band = 101;
    state.ammonia_non_band = 102;
    state.ammonium_band = 103;
    state.ammonia_band = 104;
    state.nitrate_non_band = 105;
    state.nitrate_band = 106;
    var inputs = dynamicSaltInputs();
    inputs.free_ions.iron_mol = -1;
    try applyDynamicSaltAqueousConcentrations(&state, inputs);

    try std.testing.expectEqual(@as(f64, 1), state.hydrogen);
    try std.testing.expectEqual(
        minimum_concentration_mol_per_m3,
        state.iron,
    );
    try std.testing.expectEqual(@as(f64, 2), state.carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 1), state.aluminum_hydroxide_1);
    try std.testing.expectEqual(@as(f64, 22), state.hydrogen_silicate);
    try std.testing.expectEqual(@as(f64, 101), state.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 106), state.nitrate_band);
}

test "SOLUTE dry dynamic salt initialization uses pH and clears only salts" {
    var state: aqueous_network.State = undefined;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field|
        @field(state, field.name) = 9;
    var inputs = dynamicSaltInputs();
    inputs.water_volume_m3 = inputs.minimum_water_volume_m3;
    try applyDynamicSaltAqueousConcentrations(&state, inputs);

    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydrogen, 1.0e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydroxide, 1.0e-18);
    try std.testing.expectEqual(@as(f64, 0), state.aluminum);
    try std.testing.expectEqual(@as(f64, 0), state.hydrogen_silicate);
    try std.testing.expectEqual(@as(f64, 9), state.ammonium_band);
    try std.testing.expectEqual(@as(f64, 9), state.nitrate_non_band);
}

test "SOLUTE dynamic salt initialization rejects late invalid input atomically" {
    var state: aqueous_network.State = undefined;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field|
        @field(state, field.name) = 3;
    const before = state;
    var inputs = dynamicSaltInputs();
    inputs.dissolved_complexes.potassium_sulfate_mol =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidDynamicSaltAqueousInput,
        applyDynamicSaltAqueousConcentrations(&state, inputs),
    );
    try std.testing.expectEqualDeep(before, state);
}
