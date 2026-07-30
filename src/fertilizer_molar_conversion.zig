const std = @import("std");

pub const Amounts = struct {
    broadcast_ammonium_g_n_per_m2: f64,
    broadcast_ammonia_g_n_per_m2: f64,
    broadcast_urea_g_n_per_m2: f64,
    broadcast_nitrate_g_n_per_m2: f64,
    banded_ammonium_g_n_per_m2: f64,
    banded_ammonia_g_n_per_m2: f64,
    banded_urea_g_n_per_m2: f64,
    banded_nitrate_g_n_per_m2: f64,
    broadcast_monocalcium_phosphate_g_p_per_m2: f64,
    banded_monocalcium_phosphate_g_p_per_m2: f64,
    broadcast_hydroxyapatite_g_p_per_m2: f64,
    calcium_carbonate_g_ca_per_m2: f64,
    calcium_sulfate_or_ground_rock_g_per_m2: f64,
};

pub const Parameters = struct {
    nitrogen_g_per_mol: f64,
    monocalcium_phosphate_g_p_per_mol: f64,
    hydroxyapatite_g_p_per_mol: f64,
    calcium_g_per_mol: f64,
    ground_rock_g_per_mol_si: f64,
    ground_rock_formulation_minimum_code: u8,
};

pub const Result = struct {
    broadcast_ammonium_mol_n: f64,
    broadcast_ammonia_mol_n: f64,
    broadcast_urea_mol_n: f64,
    broadcast_nitrate_mol_n: f64,
    banded_ammonium_mol_n: f64,
    banded_ammonia_mol_n: f64,
    banded_urea_mol_n: f64,
    banded_nitrate_mol_n: f64,
    broadcast_monocalcium_phosphate_mol_p: f64,
    banded_monocalcium_phosphate_mol_p: f64,
    broadcast_hydroxyapatite_mol_p: f64,
    calcium_carbonate_mol_ca: f64,
    calcium_sulfate_or_ground_rock_mol: f64,
    is_ground_rock: bool,
    /// Source `TZIN`, summed broadcast groups before banded groups.
    total_nitrogen_input_g_n: f64,
    /// Source `UFERTN`, summed by chemical species across placement modes.
    daily_nitrogen_fertilizer_g_n: f64,
    total_phosphorus_input_g_p: f64,
    dynamic_salt_ionic_equivalents_mol: f64,
};

/// Exact HOUR1 fertilizer conversion/input ledgers from hour1.f:522-603.
pub fn convert(
    amounts: Amounts,
    cell_area_m2: f64,
    formulation_code: u8,
    dynamic_salt_chemistry: bool,
    parameters: Parameters,
) !Result {
    if (!std.math.isFinite(cell_area_m2) or cell_area_m2 <= 0)
        return error.InvalidFertilizerConversionArea;
    inline for (@typeInfo(Amounts).@"struct".fields) |field| {
        const value = @field(amounts, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteFertilizerConversionAmount;
        if (value < 0) return error.InvalidFertilizerConversionAmount;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (field.type == u8) continue;
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFertilizerConversionParameter;
    }
    const is_ground_rock =
        formulation_code >= parameters.ground_rock_formulation_minimum_code;

    var result: Result = .{
        .broadcast_ammonium_mol_n = amounts.broadcast_ammonium_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .broadcast_ammonia_mol_n = amounts.broadcast_ammonia_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .broadcast_urea_mol_n = amounts.broadcast_urea_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .broadcast_nitrate_mol_n = amounts.broadcast_nitrate_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .banded_ammonium_mol_n = amounts.banded_ammonium_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .banded_ammonia_mol_n = amounts.banded_ammonia_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .banded_urea_mol_n = amounts.banded_urea_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .banded_nitrate_mol_n = amounts.banded_nitrate_g_n_per_m2 * cell_area_m2 / parameters.nitrogen_g_per_mol,
        .broadcast_monocalcium_phosphate_mol_p = amounts.broadcast_monocalcium_phosphate_g_p_per_m2 *
            cell_area_m2 / parameters.monocalcium_phosphate_g_p_per_mol,
        .banded_monocalcium_phosphate_mol_p = amounts.banded_monocalcium_phosphate_g_p_per_m2 *
            cell_area_m2 / parameters.monocalcium_phosphate_g_p_per_mol,
        .broadcast_hydroxyapatite_mol_p = amounts.broadcast_hydroxyapatite_g_p_per_m2 *
            cell_area_m2 / parameters.hydroxyapatite_g_p_per_mol,
        .calcium_carbonate_mol_ca = amounts.calcium_carbonate_g_ca_per_m2 *
            cell_area_m2 / parameters.calcium_g_per_mol,
        .calcium_sulfate_or_ground_rock_mol = amounts.calcium_sulfate_or_ground_rock_g_per_m2 *
            cell_area_m2 /
            (if (is_ground_rock)
                parameters.ground_rock_g_per_mol_si
            else
                parameters.calcium_g_per_mol),
        .is_ground_rock = is_ground_rock,
        .total_nitrogen_input_g_n = 0,
        .daily_nitrogen_fertilizer_g_n = 0,
        .total_phosphorus_input_g_p = 0,
        .dynamic_salt_ionic_equivalents_mol = 0,
    };
    result.total_nitrogen_input_g_n = parameters.nitrogen_g_per_mol *
        (result.broadcast_ammonium_mol_n +
            result.broadcast_ammonia_mol_n +
            result.broadcast_urea_mol_n +
            result.broadcast_nitrate_mol_n +
            result.banded_ammonium_mol_n +
            result.banded_ammonia_mol_n +
            result.banded_urea_mol_n +
            result.banded_nitrate_mol_n);
    result.daily_nitrogen_fertilizer_g_n = parameters.nitrogen_g_per_mol *
        (result.broadcast_ammonium_mol_n +
            result.banded_ammonium_mol_n +
            result.broadcast_ammonia_mol_n +
            result.banded_ammonia_mol_n +
            result.broadcast_urea_mol_n +
            result.banded_urea_mol_n +
            result.broadcast_nitrate_mol_n +
            result.banded_nitrate_mol_n);
    result.total_phosphorus_input_g_p =
        parameters.monocalcium_phosphate_g_p_per_mol *
        (result.broadcast_monocalcium_phosphate_mol_p +
            result.banded_monocalcium_phosphate_mol_p) +
        parameters.hydroxyapatite_g_p_per_mol *
            result.broadcast_hydroxyapatite_mol_p;
    if (dynamic_salt_chemistry) {
        result.dynamic_salt_ionic_equivalents_mol = if (!is_ground_rock)
            2 * (result.calcium_carbonate_mol_ca +
                result.calcium_sulfate_or_ground_rock_mol +
                result.broadcast_ammonium_mol_n +
                result.banded_ammonium_mol_n) +
                result.broadcast_urea_mol_n +
                result.broadcast_nitrate_mol_n +
                result.banded_urea_mol_n +
                result.banded_nitrate_mol_n +
                7 * (result.broadcast_monocalcium_phosphate_mol_p +
                    result.banded_monocalcium_phosphate_mol_p) +
                9 * result.broadcast_hydroxyapatite_mol_p
        else
            2 * (result.calcium_carbonate_mol_ca +
                result.broadcast_ammonium_mol_n +
                result.banded_ammonium_mol_n) +
                result.calcium_sulfate_or_ground_rock_mol +
                result.broadcast_urea_mol_n +
                result.broadcast_nitrate_mol_n +
                result.banded_urea_mol_n +
                result.banded_nitrate_mol_n +
                7 * (result.broadcast_monocalcium_phosphate_mol_p +
                    result.banded_monocalcium_phosphate_mol_p) +
                9 * result.broadcast_hydroxyapatite_mol_p;
    }
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(result, field.name)))
            return error.FertilizerConversionOverflow;
    return result;
}

pub const sourceParameters: Parameters = .{
    .nitrogen_g_per_mol = 14,
    .monocalcium_phosphate_g_p_per_mol = 62,
    .hydroxyapatite_g_p_per_mol = 93,
    .calcium_g_per_mol = 40,
    .ground_rock_g_per_mol_si = 92,
    .ground_rock_formulation_minimum_code = 10,
};

fn exampleAmounts() Amounts {
    var result: Amounts = undefined;
    inline for (@typeInfo(Amounts).@"struct".fields) |field|
        @field(result, field.name) = 1;
    return result;
}

test "source conversion conserves runtime area-scaled N and P mass" {
    const result = try convert(exampleAmounts(), 14, 0, true, sourceParameters);
    try std.testing.expectEqual(@as(f64, 112), result.total_nitrogen_input_g_n);
    try std.testing.expectEqual(@as(f64, 112), result.daily_nitrogen_fertilizer_g_n);
    try std.testing.expectEqual(@as(f64, 42), result.total_phosphorus_input_g_p);
    try std.testing.expect(!result.is_ground_rock);
}

test "distinct nitrogen ledgers preserve both source summation orders" {
    var amounts = exampleAmounts();
    amounts.broadcast_ammonium_g_n_per_m2 = 1.0e16;
    amounts.broadcast_ammonia_g_n_per_m2 = 1;
    amounts.broadcast_urea_g_n_per_m2 = 1;
    amounts.broadcast_nitrate_g_n_per_m2 = 1;
    amounts.banded_ammonium_g_n_per_m2 = 1;
    amounts.banded_ammonia_g_n_per_m2 = 1;
    amounts.banded_urea_g_n_per_m2 = 1;
    amounts.banded_nitrate_g_n_per_m2 = 1;
    const result = try convert(amounts, 14, 0, false, sourceParameters);
    var tzin_sum =
        result.broadcast_ammonium_mol_n +
        result.broadcast_ammonia_mol_n +
        result.broadcast_urea_mol_n +
        result.broadcast_nitrate_mol_n +
        result.banded_ammonium_mol_n +
        result.banded_ammonia_mol_n +
        result.banded_urea_mol_n +
        result.banded_nitrate_mol_n;
    tzin_sum *= 14;
    var ufertn_sum =
        result.broadcast_ammonium_mol_n +
        result.banded_ammonium_mol_n +
        result.broadcast_ammonia_mol_n +
        result.banded_ammonia_mol_n +
        result.broadcast_urea_mol_n +
        result.banded_urea_mol_n +
        result.broadcast_nitrate_mol_n +
        result.banded_nitrate_mol_n;
    ufertn_sum *= 14;
    try std.testing.expectEqual(tzin_sum, result.total_nitrogen_input_g_n);
    try std.testing.expectEqual(ufertn_sum, result.daily_nitrogen_fertilizer_g_n);
}

test "gypsum and ground rock use distinct exact denominators and charge" {
    var input = exampleAmounts();
    inline for (@typeInfo(Amounts).@"struct".fields) |field|
        @field(input, field.name) = 0;
    input.calcium_sulfate_or_ground_rock_g_per_m2 = 40;
    const gypsum = try convert(input, 1, 9, true, sourceParameters);
    const rock = try convert(input, 1, 10, true, sourceParameters);
    try std.testing.expectEqual(@as(f64, 1), gypsum.calcium_sulfate_or_ground_rock_mol);
    try std.testing.expectEqual(@as(f64, 2), gypsum.dynamic_salt_ionic_equivalents_mol);
    try std.testing.expectApproxEqAbs(
        @as(f64, 40.0 / 92.0),
        rock.calcium_sulfate_or_ground_rock_mol,
        1e-15,
    );
    try std.testing.expectEqual(
        rock.calcium_sulfate_or_ground_rock_mol,
        rock.dynamic_salt_ionic_equivalents_mol,
    );
}

test "static salt mode reports no dynamic ionic-equivalent input" {
    const result = try convert(exampleAmounts(), 1, 0, false, sourceParameters);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.dynamic_salt_ionic_equivalents_mol,
    );
}

test "invalid late amount fails before returning conversion" {
    var input = exampleAmounts();
    input.calcium_sulfate_or_ground_rock_g_per_m2 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteFertilizerConversionAmount,
        convert(input, 1, 0, true, sourceParameters),
    );
}
