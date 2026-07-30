const std = @import("std");

pub const PerAreaAmounts = struct {
    carbon_g_c_per_m2: f64,
    nitrogen_g_n_per_m2: f64,
    phosphorus_g_p_per_m2: f64,
};

pub const ExtensiveAmounts = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Result = struct {
    input: ExtensiveAmounts,
    allocated_microbial_and_dissolved: ExtensiveAmounts,
};

/// HOUR1 lines 762--773. The caller invokes this in source K=plant,manure
/// order; multiplication and zero assignments retain source order.
pub fn prepare(amounts: PerAreaAmounts, cell_area_m2: f64) !Result {
    inline for (.{
        amounts.carbon_g_c_per_m2,
        amounts.nitrogen_g_n_per_m2,
        amounts.phosphorus_g_p_per_m2,
        cell_area_m2,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteOrganicAmendmentPreparationInput;
    if (amounts.carbon_g_c_per_m2 < 0 or
        amounts.nitrogen_g_n_per_m2 < 0 or
        amounts.phosphorus_g_p_per_m2 < 0 or
        cell_area_m2 <= 0)
        return error.InvalidOrganicAmendmentPreparationInput;

    var result: Result = undefined;
    result.input.carbon_g_c = amounts.carbon_g_c_per_m2 * cell_area_m2;
    result.input.nitrogen_g_n = amounts.nitrogen_g_n_per_m2 * cell_area_m2;
    result.input.phosphorus_g_p = amounts.phosphorus_g_p_per_m2 * cell_area_m2;
    result.allocated_microbial_and_dissolved.carbon_g_c = 0;
    result.allocated_microbial_and_dissolved.nitrogen_g_n = 0;
    result.allocated_microbial_and_dissolved.phosphorus_g_p = 0;
    inline for (@typeInfo(Result).@"struct".fields) |group_field|
        inline for (@typeInfo(ExtensiveAmounts).@"struct".fields) |field|
            if (!std.math.isFinite(@field(@field(result, group_field.name), field.name)))
                return error.OrganicAmendmentPreparationOverflow;
    return result;
}

test "organic amendment preparation preserves source assignments" {
    const result = try prepare(.{
        .carbon_g_c_per_m2 = 2,
        .nitrogen_g_n_per_m2 = 3,
        .phosphorus_g_p_per_m2 = 4,
    }, 5);
    try std.testing.expectEqual(@as(f64, 10), result.input.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 15), result.input.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 20), result.input.phosphorus_g_p);
    try std.testing.expectEqualDeep(
        ExtensiveAmounts{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        result.allocated_microbial_and_dissolved,
    );
}

test "late nonfinite amount fails before publishing preparation" {
    try std.testing.expectError(
        error.NonFiniteOrganicAmendmentPreparationInput,
        prepare(.{
            .carbon_g_c_per_m2 = 1,
            .nitrogen_g_n_per_m2 = 2,
            .phosphorus_g_p_per_m2 = std.math.nan(f64),
        }, 5),
    );
}
