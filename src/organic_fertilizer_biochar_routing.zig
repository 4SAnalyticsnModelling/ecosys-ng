const std = @import("std");

pub const AmendmentKind = enum {
    plant_residue,
    animal_manure,
};

pub const Amounts = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

/// HOUR1 lines 788--796. Plant-residue type ten routes the complete
/// amendment to biochar; manure and other plant types continue downstream.
pub fn routedAmounts(
    kind: AmendmentKind,
    material_type: u8,
    input: Amounts,
) !?Amounts {
    inline for (@typeInfo(Amounts).@"struct".fields) |field| {
        const value = @field(input, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteOrganicBiocharInput;
        if (value < 0) return error.InvalidOrganicBiocharInput;
    }
    return if (kind == .plant_residue and material_type == 10)
        input
    else
        null;
}

test "plant type ten routes complete C N P to biochar" {
    const input: Amounts = .{
        .carbon_g_c = 10,
        .nitrogen_g_n = 2,
        .phosphorus_g_p = 1,
    };
    try std.testing.expectEqual(
        input,
        (try routedAmounts(.plant_residue, 10, input)).?,
    );
}

test "manure type ten and ordinary plant residue continue downstream" {
    const input: Amounts = .{
        .carbon_g_c = 10,
        .nitrogen_g_n = 2,
        .phosphorus_g_p = 1,
    };
    try std.testing.expect((try routedAmounts(.animal_manure, 10, input)) == null);
    try std.testing.expect((try routedAmounts(.plant_residue, 9, input)) == null);
}

test "invalid biochar input fails before routing" {
    try std.testing.expectError(
        error.NonFiniteOrganicBiocharInput,
        routedAmounts(.plant_residue, 10, .{
            .carbon_g_c = 1,
            .nitrogen_g_n = std.math.nan(f64),
            .phosphorus_g_p = 0,
        }),
    );
}
