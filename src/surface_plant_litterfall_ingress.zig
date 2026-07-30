const std = @import("std");

pub const LitterComplex = enum(u8) {
    woody,
    nonwoody,

    pub const count: usize = @typeInfo(LitterComplex).@"enum".fields.len;
};

pub const structural_fraction_count: usize = 5;
pub const noncharcoal_fraction_count: usize = 4;

pub const ElementalMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const Inputs = struct {
    column_count: usize,
    row_count: usize,
    litterfall_g_by_cell_complex_fraction: []const ElementalMass,
};

pub const State = struct {
    surface_residue_g_by_cell_complex_fraction: []ElementalMass,
    noncharcoal_organic_carbon_g_c_by_cell: []f64,
    charcoal_carbon_g_c_by_cell: []f64,
    charcoal_carbon_before_ingress_g_c_by_cell: []f64,
};

/// Commits accepted plant litterfall to the two surface residue complexes.
///
/// Traceability: REDIST.F lines 252--270 (`ORGCCX`, `CSNT`, `ZSNT`, `PSNT`,
/// `OSC`, `OSN`, `OSP`, `ORGC`, and `ORGCC`). The runtime storage order is
/// cell, litter complex, then structural fraction. The source traversal is
/// retained as column, row, complex, then fraction. Fractions zero through
/// three contribute to non-charcoal organic C; fraction four contributes to
/// charcoal C. Validation completes before any caller-owned state changes.
pub fn apply(inputs: Inputs, state: *State) !void {
    const cell_count = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, cell_count);
    try preflightUpdates(inputs, state.*);

    for (0..inputs.column_count) |column| {
        for (0..inputs.row_count) |row| {
            const cell = row * inputs.column_count + column;
            state.charcoal_carbon_before_ingress_g_c_by_cell[cell] =
                state.charcoal_carbon_g_c_by_cell[cell];
            for (0..LitterComplex.count) |complex| {
                for (0..structural_fraction_count) |fraction| {
                    const index = residueIndex(cell, complex, fraction);
                    const litterfall = inputs.litterfall_g_by_cell_complex_fraction[index];
                    state.surface_residue_g_by_cell_complex_fraction[index].carbon_g_c +=
                        litterfall.carbon_g_c;
                    state.surface_residue_g_by_cell_complex_fraction[index].nitrogen_g_n +=
                        litterfall.nitrogen_g_n;
                    state.surface_residue_g_by_cell_complex_fraction[index].phosphorus_g_p +=
                        litterfall.phosphorus_g_p;
                    if (fraction < noncharcoal_fraction_count) {
                        state.noncharcoal_organic_carbon_g_c_by_cell[cell] +=
                            litterfall.carbon_g_c;
                    } else {
                        state.charcoal_carbon_g_c_by_cell[cell] +=
                            litterfall.carbon_g_c;
                    }
                }
            }
        }
    }
}

fn validateDimensions(inputs: Inputs, state: State) !usize {
    if (inputs.column_count == 0 or inputs.row_count == 0)
        return error.InvalidPlantLitterfallIngressDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.column_count,
        inputs.row_count,
    ) catch return error.InvalidPlantLitterfallIngressDimensions;
    const complex_value_count = std.math.mul(
        usize,
        cell_count,
        LitterComplex.count,
    ) catch return error.InvalidPlantLitterfallIngressDimensions;
    const residue_value_count = std.math.mul(
        usize,
        complex_value_count,
        structural_fraction_count,
    ) catch return error.InvalidPlantLitterfallIngressDimensions;
    if (inputs.litterfall_g_by_cell_complex_fraction.len != residue_value_count or
        state.surface_residue_g_by_cell_complex_fraction.len != residue_value_count or
        state.noncharcoal_organic_carbon_g_c_by_cell.len != cell_count or
        state.charcoal_carbon_g_c_by_cell.len != cell_count or
        state.charcoal_carbon_before_ingress_g_c_by_cell.len != cell_count)
        return error.InvalidPlantLitterfallIngressDimensions;
    return cell_count;
}

fn validateInputsAndState(inputs: Inputs, state: State, cell_count: usize) !void {
    for (inputs.litterfall_g_by_cell_complex_fraction) |mass|
        try validateMass(mass, error.InvalidPlantLitterfallIngressInput);
    for (state.surface_residue_g_by_cell_complex_fraction) |mass|
        try validateMass(mass, error.InvalidPlantLitterfallIngressState);
    for (0..cell_count) |cell| {
        if (!nonnegativeFinite(state.noncharcoal_organic_carbon_g_c_by_cell[cell]) or
            !nonnegativeFinite(state.charcoal_carbon_g_c_by_cell[cell]) or
            !nonnegativeFinite(state.charcoal_carbon_before_ingress_g_c_by_cell[cell]))
            return error.InvalidPlantLitterfallIngressState;
    }
}

fn preflightUpdates(inputs: Inputs, state: State) !void {
    for (0..inputs.column_count) |column| {
        for (0..inputs.row_count) |row| {
            const cell = row * inputs.column_count + column;
            var next_noncharcoal_g_c =
                state.noncharcoal_organic_carbon_g_c_by_cell[cell];
            var next_charcoal_g_c = state.charcoal_carbon_g_c_by_cell[cell];
            for (0..LitterComplex.count) |complex| {
                for (0..structural_fraction_count) |fraction| {
                    const index = residueIndex(cell, complex, fraction);
                    const current = state.surface_residue_g_by_cell_complex_fraction[index];
                    const litterfall = inputs.litterfall_g_by_cell_complex_fraction[index];
                    _ = try checkedSum(current.carbon_g_c, litterfall.carbon_g_c);
                    _ = try checkedSum(current.nitrogen_g_n, litterfall.nitrogen_g_n);
                    _ = try checkedSum(
                        current.phosphorus_g_p,
                        litterfall.phosphorus_g_p,
                    );
                    if (fraction < noncharcoal_fraction_count) {
                        next_noncharcoal_g_c = try checkedSum(
                            next_noncharcoal_g_c,
                            litterfall.carbon_g_c,
                        );
                    } else {
                        next_charcoal_g_c = try checkedSum(
                            next_charcoal_g_c,
                            litterfall.carbon_g_c,
                        );
                    }
                }
            }
        }
    }
}

fn residueIndex(cell: usize, complex: usize, fraction: usize) usize {
    return (cell * LitterComplex.count + complex) * structural_fraction_count +
        fraction;
}

fn validateMass(mass: ElementalMass, comptime failure: anyerror) !void {
    inline for (@typeInfo(ElementalMass).@"struct".fields) |field|
        if (!nonnegativeFinite(@field(mass, field.name))) return failure;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!nonnegativeFinite(result))
        return error.NonFinitePlantLitterfallIngressResult;
    return result;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

test "REDIST plant litterfall conserves C N P and separates charcoal" {
    const value_count = LitterComplex.count * structural_fraction_count;
    var litterfall: [value_count]ElementalMass = undefined;
    for (&litterfall, 0..) |*mass, index| {
        const carbon_g_c: f64 = @floatFromInt(index + 1);
        mass.* = .{
            .carbon_g_c = carbon_g_c,
            .nitrogen_g_n = 0.1 * carbon_g_c,
            .phosphorus_g_p = 0.01 * carbon_g_c,
        };
    }
    var residue = [_]ElementalMass{.{ .carbon_g_c = 1, .nitrogen_g_n = 2, .phosphorus_g_p = 3 }} ** value_count;
    var noncharcoal = [_]f64{100};
    var charcoal = [_]f64{10};
    var charcoal_before = [_]f64{999};
    var state: State = .{
        .surface_residue_g_by_cell_complex_fraction = &residue,
        .noncharcoal_organic_carbon_g_c_by_cell = &noncharcoal,
        .charcoal_carbon_g_c_by_cell = &charcoal,
        .charcoal_carbon_before_ingress_g_c_by_cell = &charcoal_before,
    };
    try apply(.{
        .column_count = 1,
        .row_count = 1,
        .litterfall_g_by_cell_complex_fraction = &litterfall,
    }, &state);

    var residue_carbon_increment_g_c: f64 = 0;
    var residue_nitrogen_increment_g_n: f64 = 0;
    var residue_phosphorus_increment_g_p: f64 = 0;
    for (residue) |mass| {
        residue_carbon_increment_g_c += mass.carbon_g_c - 1;
        residue_nitrogen_increment_g_n += mass.nitrogen_g_n - 2;
        residue_phosphorus_increment_g_p += mass.phosphorus_g_p - 3;
    }
    try std.testing.expectApproxEqAbs(
        @as(f64, 55),
        residue_carbon_increment_g_c,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 5.5),
        residue_nitrogen_increment_g_n,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.55),
        residue_phosphorus_increment_g_p,
        1e-14,
    );
    try std.testing.expectEqual(@as(f64, 140), noncharcoal[0]);
    try std.testing.expectEqual(@as(f64, 25), charcoal[0]);
    try std.testing.expectEqual(@as(f64, 10), charcoal_before[0]);
    try std.testing.expectEqual(
        residue_carbon_increment_g_c,
        (noncharcoal[0] - 100) + (charcoal[0] - 10),
    );
}

test "runtime grid snapshots charcoal and maps every cell independently" {
    const cell_count = 6;
    const value_count =
        cell_count * LitterComplex.count * structural_fraction_count;
    var litterfall = [_]ElementalMass{.{}} ** value_count;
    litterfall[residueIndex(1, @intFromEnum(LitterComplex.nonwoody), 2)] =
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03 };
    litterfall[residueIndex(5, @intFromEnum(LitterComplex.woody), 4)] =
        .{ .carbon_g_c = 7, .nitrogen_g_n = 0.7, .phosphorus_g_p = 0.07 };
    var residue = [_]ElementalMass{.{}} ** value_count;
    var noncharcoal = [_]f64{0} ** cell_count;
    var charcoal = [_]f64{ 10, 11, 12, 13, 14, 15 };
    var charcoal_before = [_]f64{999} ** cell_count;
    var state: State = .{
        .surface_residue_g_by_cell_complex_fraction = &residue,
        .noncharcoal_organic_carbon_g_c_by_cell = &noncharcoal,
        .charcoal_carbon_g_c_by_cell = &charcoal,
        .charcoal_carbon_before_ingress_g_c_by_cell = &charcoal_before,
    };
    try apply(.{
        .column_count = 3,
        .row_count = 2,
        .litterfall_g_by_cell_complex_fraction = &litterfall,
    }, &state);

    try std.testing.expectEqualSlices(
        f64,
        &.{ 10, 11, 12, 13, 14, 15 },
        &charcoal_before,
    );
    try std.testing.expectEqual(@as(f64, 3), noncharcoal[1]);
    try std.testing.expectEqual(@as(f64, 22), charcoal[5]);
    try std.testing.expectEqual(
        @as(f64, 3),
        residue[residueIndex(1, @intFromEnum(LitterComplex.nonwoody), 2)].carbon_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 7),
        residue[residueIndex(5, @intFromEnum(LitterComplex.woody), 4)].carbon_g_c,
    );
}

test "late invalid litterfall leaves all state unchanged" {
    const cell_count = 2;
    const value_count =
        cell_count * LitterComplex.count * structural_fraction_count;
    var litterfall = [_]ElementalMass{.{}} ** value_count;
    litterfall[0].carbon_g_c = 1;
    litterfall[value_count - 1].phosphorus_g_p = std.math.nan(f64);
    var residue = [_]ElementalMass{.{ .carbon_g_c = 2, .nitrogen_g_n = 3, .phosphorus_g_p = 4 }} ** value_count;
    var noncharcoal = [_]f64{ 5, 6 };
    var charcoal = [_]f64{ 7, 8 };
    var charcoal_before = [_]f64{ 9, 10 };
    var state: State = .{
        .surface_residue_g_by_cell_complex_fraction = &residue,
        .noncharcoal_organic_carbon_g_c_by_cell = &noncharcoal,
        .charcoal_carbon_g_c_by_cell = &charcoal,
        .charcoal_carbon_before_ingress_g_c_by_cell = &charcoal_before,
    };
    try std.testing.expectError(
        error.InvalidPlantLitterfallIngressInput,
        apply(.{
            .column_count = 2,
            .row_count = 1,
            .litterfall_g_by_cell_complex_fraction = &litterfall,
        }, &state),
    );

    for (residue) |mass| try std.testing.expectEqual(
        ElementalMass{ .carbon_g_c = 2, .nitrogen_g_n = 3, .phosphorus_g_p = 4 },
        mass,
    );
    try std.testing.expectEqualSlices(f64, &.{ 5, 6 }, &noncharcoal);
    try std.testing.expectEqualSlices(f64, &.{ 7, 8 }, &charcoal);
    try std.testing.expectEqualSlices(f64, &.{ 9, 10 }, &charcoal_before);
}

test "invalid runtime dimensions and nonfinite state fail explicitly" {
    const value_count = LitterComplex.count * structural_fraction_count;
    const litterfall = [_]ElementalMass{.{}} ** value_count;
    var residue = [_]ElementalMass{.{}} ** value_count;
    var noncharcoal = [_]f64{0};
    var charcoal = [_]f64{0};
    var charcoal_before = [_]f64{0};
    var state: State = .{
        .surface_residue_g_by_cell_complex_fraction = &residue,
        .noncharcoal_organic_carbon_g_c_by_cell = &noncharcoal,
        .charcoal_carbon_g_c_by_cell = &charcoal,
        .charcoal_carbon_before_ingress_g_c_by_cell = &charcoal_before,
    };
    try std.testing.expectError(
        error.InvalidPlantLitterfallIngressDimensions,
        apply(.{
            .column_count = 0,
            .row_count = 1,
            .litterfall_g_by_cell_complex_fraction = &litterfall,
        }, &state),
    );
    charcoal[0] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidPlantLitterfallIngressState,
        apply(.{
            .column_count = 1,
            .row_count = 1,
            .litterfall_g_by_cell_complex_fraction = &litterfall,
        }, &state),
    );
}
