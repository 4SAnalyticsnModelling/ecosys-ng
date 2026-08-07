const std = @import("std");

pub const biochemical_fraction_count: usize = 4;

pub const ElementalMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const UreaFormulation = enum(u8) {
    fast_release,
    normal_release,
    slow_release,
};

pub const ConversionParameters = struct {
    nitrogen_g_n_per_mol_n: f64,
    phosphorus_g_p_per_mol_monocalcium_phosphate: f64,
    ion_mol_per_mol_monocalcium_phosphate: f64,
};

pub const Inputs = struct {
    column_count: usize,
    row_count: usize,
    soil_layer_capacity: usize,
    top_soil_layer_by_cell: []const usize,
    negligible_carbon_g_c_by_cell: []const f64,
    deposited_organic_g_by_cell_fraction: []const ElementalMass,
    deposited_inorganic_nitrogen_g_n_by_cell: []const f64,
    deposited_inorganic_phosphorus_g_p_by_cell: []const f64,
    salt_equilibrium_mode: SaltEquilibriumMode,
    conversion: ConversionParameters,
};

pub const State = struct {
    surface_residue_g_by_cell_fraction: []ElementalMass,
    surface_active_carbon_g_c_by_cell_fraction: []f64,
    broadcast_urea_mol_n_by_soil_layer: []f64,
    broadcast_monocalcium_phosphate_mol_by_soil_layer: []f64,
    urea_formulation_by_cell: []UreaFormulation,
    total_dynamic_ion_input_mol: f64,
};

/// Transfers accepted manure deposition into surface residue and the current
/// top mineral layer.
///
/// Traceability: REDIST.F lines 224--243 (`CSNMT`, `ZSNMT`, `PSNMT`,
/// `ZSNIT`, `PSNIT`, `ZNHUFA`, `PCAPM`, `IUTYP`, and `TIONIN`). The removed
/// sub-hour full-model cycle makes this one accepted-hour transaction
/// equivalent to the source `NFZ == 1` gate. Cells are traversed in the
/// source's column-outer, row-inner order. All extensive inventories are
/// validated before any caller-owned state is mutated.
pub fn apply(inputs: Inputs, state: *State) !void {
    const cell_count = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, cell_count);
    try preflightUpdates(inputs, state.*, cell_count);

    for (0..inputs.column_count) |column| {
        for (0..inputs.row_count) |row| {
            const cell = row * inputs.column_count + column;
            const first_fraction = cell * biochemical_fraction_count;
            if (inputs.deposited_organic_g_by_cell_fraction[first_fraction].carbon_g_c <=
                inputs.negligible_carbon_g_c_by_cell[cell]) continue;

            for (0..biochemical_fraction_count) |fraction| {
                const index = first_fraction + fraction;
                const deposited = inputs.deposited_organic_g_by_cell_fraction[index];
                state.surface_residue_g_by_cell_fraction[index].carbon_g_c +=
                    deposited.carbon_g_c;
                state.surface_active_carbon_g_c_by_cell_fraction[index] +=
                    deposited.carbon_g_c;
                state.surface_residue_g_by_cell_fraction[index].nitrogen_g_n +=
                    deposited.nitrogen_g_n;
                state.surface_residue_g_by_cell_fraction[index].phosphorus_g_p +=
                    deposited.phosphorus_g_p;
            }

            const soil_index =
                cell * inputs.soil_layer_capacity + inputs.top_soil_layer_by_cell[cell];
            const inorganic_nitrogen_g_n =
                inputs.deposited_inorganic_nitrogen_g_n_by_cell[cell];
            const inorganic_phosphorus_g_p =
                inputs.deposited_inorganic_phosphorus_g_p_by_cell[cell];
            state.broadcast_urea_mol_n_by_soil_layer[soil_index] +=
                inorganic_nitrogen_g_n /
                inputs.conversion.nitrogen_g_n_per_mol_n;
            state.broadcast_monocalcium_phosphate_mol_by_soil_layer[soil_index] +=
                inorganic_phosphorus_g_p /
                inputs.conversion.phosphorus_g_p_per_mol_monocalcium_phosphate;
            state.urea_formulation_by_cell[cell] = .fast_release;
            if (inputs.salt_equilibrium_mode == .dynamic) {
                state.total_dynamic_ion_input_mol +=
                    inorganic_nitrogen_g_n /
                    inputs.conversion.nitrogen_g_n_per_mol_n;
                state.total_dynamic_ion_input_mol +=
                    inputs.conversion.ion_mol_per_mol_monocalcium_phosphate *
                    inorganic_phosphorus_g_p /
                    inputs.conversion.phosphorus_g_p_per_mol_monocalcium_phosphate;
            }
        }
    }
}

fn validateDimensions(inputs: Inputs, state: State) !usize {
    if (inputs.column_count == 0 or
        inputs.row_count == 0 or
        inputs.soil_layer_capacity == 0)
        return error.InvalidManureIngressDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.column_count,
        inputs.row_count,
    ) catch return error.InvalidManureIngressDimensions;
    const fraction_value_count = std.math.mul(
        usize,
        cell_count,
        biochemical_fraction_count,
    ) catch return error.InvalidManureIngressDimensions;
    const soil_value_count = std.math.mul(
        usize,
        cell_count,
        inputs.soil_layer_capacity,
    ) catch return error.InvalidManureIngressDimensions;
    if (inputs.top_soil_layer_by_cell.len != cell_count or
        inputs.negligible_carbon_g_c_by_cell.len != cell_count or
        inputs.deposited_organic_g_by_cell_fraction.len != fraction_value_count or
        inputs.deposited_inorganic_nitrogen_g_n_by_cell.len != cell_count or
        inputs.deposited_inorganic_phosphorus_g_p_by_cell.len != cell_count or
        state.surface_residue_g_by_cell_fraction.len != fraction_value_count or
        state.surface_active_carbon_g_c_by_cell_fraction.len != fraction_value_count or
        state.broadcast_urea_mol_n_by_soil_layer.len != soil_value_count or
        state.broadcast_monocalcium_phosphate_mol_by_soil_layer.len != soil_value_count or
        state.urea_formulation_by_cell.len != cell_count)
        return error.InvalidManureIngressDimensions;
    return cell_count;
}

fn validateInputsAndState(inputs: Inputs, state: State, cell_count: usize) !void {
    const conversion = inputs.conversion;
    if (!positiveFinite(conversion.nitrogen_g_n_per_mol_n) or
        !positiveFinite(conversion.phosphorus_g_p_per_mol_monocalcium_phosphate) or
        !nonnegativeFinite(conversion.ion_mol_per_mol_monocalcium_phosphate) or
        !nonnegativeFinite(state.total_dynamic_ion_input_mol))
        return error.InvalidManureIngressState;

    for (0..cell_count) |cell| {
        if (inputs.top_soil_layer_by_cell[cell] >= inputs.soil_layer_capacity or
            !nonnegativeFinite(inputs.negligible_carbon_g_c_by_cell[cell]) or
            !nonnegativeFinite(inputs.deposited_inorganic_nitrogen_g_n_by_cell[cell]) or
            !nonnegativeFinite(inputs.deposited_inorganic_phosphorus_g_p_by_cell[cell]))
            return error.InvalidManureIngressInput;
    }
    for (inputs.deposited_organic_g_by_cell_fraction) |mass|
        try validateMass(mass, error.InvalidManureIngressInput);
    for (state.surface_residue_g_by_cell_fraction) |mass|
        try validateMass(mass, error.InvalidManureIngressState);
    for (state.surface_active_carbon_g_c_by_cell_fraction) |carbon_g_c|
        if (!nonnegativeFinite(carbon_g_c))
            return error.InvalidManureIngressState;
    for (state.broadcast_urea_mol_n_by_soil_layer) |amount_mol_n|
        if (!nonnegativeFinite(amount_mol_n))
            return error.InvalidManureIngressState;
    for (state.broadcast_monocalcium_phosphate_mol_by_soil_layer) |amount_mol|
        if (!nonnegativeFinite(amount_mol))
            return error.InvalidManureIngressState;
}

fn preflightUpdates(inputs: Inputs, state: State, cell_count: usize) !void {
    _ = cell_count;
    var next_total_dynamic_ion_input_mol = state.total_dynamic_ion_input_mol;
    for (0..inputs.column_count) |column| {
        for (0..inputs.row_count) |row| {
            const cell = row * inputs.column_count + column;
            const first_fraction = cell * biochemical_fraction_count;
            if (inputs.deposited_organic_g_by_cell_fraction[first_fraction].carbon_g_c <=
                inputs.negligible_carbon_g_c_by_cell[cell]) continue;

            for (0..biochemical_fraction_count) |fraction| {
                const index = first_fraction + fraction;
                const deposited = inputs.deposited_organic_g_by_cell_fraction[index];
                const current = state.surface_residue_g_by_cell_fraction[index];
                _ = try checkedSum(current.carbon_g_c, deposited.carbon_g_c);
                _ = try checkedSum(
                    state.surface_active_carbon_g_c_by_cell_fraction[index],
                    deposited.carbon_g_c,
                );
                _ = try checkedSum(current.nitrogen_g_n, deposited.nitrogen_g_n);
                _ = try checkedSum(current.phosphorus_g_p, deposited.phosphorus_g_p);
            }

            const soil_index =
                cell * inputs.soil_layer_capacity + inputs.top_soil_layer_by_cell[cell];
            const nitrogen_mol_n =
                inputs.deposited_inorganic_nitrogen_g_n_by_cell[cell] /
                inputs.conversion.nitrogen_g_n_per_mol_n;
            const phosphate_mol =
                inputs.deposited_inorganic_phosphorus_g_p_by_cell[cell] /
                inputs.conversion.phosphorus_g_p_per_mol_monocalcium_phosphate;
            _ = try checkedSum(
                state.broadcast_urea_mol_n_by_soil_layer[soil_index],
                nitrogen_mol_n,
            );
            _ = try checkedSum(
                state.broadcast_monocalcium_phosphate_mol_by_soil_layer[soil_index],
                phosphate_mol,
            );
            if (inputs.salt_equilibrium_mode == .dynamic) {
                next_total_dynamic_ion_input_mol = try checkedSum(
                    next_total_dynamic_ion_input_mol,
                    nitrogen_mol_n,
                );
                next_total_dynamic_ion_input_mol = try checkedSum(
                    next_total_dynamic_ion_input_mol,
                    inputs.conversion.ion_mol_per_mol_monocalcium_phosphate *
                        inputs.deposited_inorganic_phosphorus_g_p_by_cell[cell] /
                        inputs.conversion.phosphorus_g_p_per_mol_monocalcium_phosphate,
                );
            }
        }
    }
}

fn validateMass(mass: ElementalMass, comptime failure: anyerror) !void {
    inline for (@typeInfo(ElementalMass).@"struct".fields) |field|
        if (!nonnegativeFinite(@field(mass, field.name))) return failure;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!nonnegativeFinite(result)) return error.NonFiniteManureIngressResult;
    return result;
}

fn positiveFinite(value: f64) bool {
    return std.math.isFinite(value) and value > 0;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

test "REDIST manure ingress preserves source arithmetic and elemental mass" {
    const deposited = [_]ElementalMass{
        .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 },
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 },
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03 },
        .{ .carbon_g_c = 4, .nitrogen_g_n = 0.4, .phosphorus_g_p = 0.04 },
    };
    var residue = [_]ElementalMass{.{}} ** biochemical_fraction_count;
    var active_carbon = [_]f64{0} ** biochemical_fraction_count;
    var urea = [_]f64{ 0, 5 };
    var monocalcium_phosphate = [_]f64{ 0, 7 };
    var formulation = [_]UreaFormulation{.slow_release};
    var state: State = .{
        .surface_residue_g_by_cell_fraction = &residue,
        .surface_active_carbon_g_c_by_cell_fraction = &active_carbon,
        .broadcast_urea_mol_n_by_soil_layer = &urea,
        .broadcast_monocalcium_phosphate_mol_by_soil_layer = &monocalcium_phosphate,
        .urea_formulation_by_cell = &formulation,
        .total_dynamic_ion_input_mol = 11,
    };
    try apply(.{
        .column_count = 1,
        .row_count = 1,
        .soil_layer_capacity = 2,
        .top_soil_layer_by_cell = &.{1},
        .negligible_carbon_g_c_by_cell = &.{1e-12},
        .deposited_organic_g_by_cell_fraction = &deposited,
        .deposited_inorganic_nitrogen_g_n_by_cell = &.{28},
        .deposited_inorganic_phosphorus_g_p_by_cell = &.{62},
        .salt_equilibrium_mode = .dynamic,
        .conversion = .{
            .nitrogen_g_n_per_mol_n = 14,
            .phosphorus_g_p_per_mol_monocalcium_phosphate = 62,
            .ion_mol_per_mol_monocalcium_phosphate = 7,
        },
    }, &state);

    var carbon_g_c: f64 = 0;
    var nitrogen_g_n: f64 = 0;
    var phosphorus_g_p: f64 = 0;
    for (state.surface_residue_g_by_cell_fraction) |mass| {
        carbon_g_c += mass.carbon_g_c;
        nitrogen_g_n += mass.nitrogen_g_n;
        phosphorus_g_p += mass.phosphorus_g_p;
    }
    try std.testing.expectApproxEqAbs(@as(f64, 10), carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), phosphorus_g_p, 1e-15);
    try std.testing.expectEqual(@as(f64, 7), urea[1]);
    try std.testing.expectEqual(@as(f64, 8), monocalcium_phosphate[1]);
    try std.testing.expectEqual(UreaFormulation.fast_release, formulation[0]);
    try std.testing.expectEqual(@as(f64, 20), state.total_dynamic_ion_input_mol);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, &active_carbon);
    try std.testing.expectApproxEqAbs(
        @as(f64, 28),
        (urea[1] - 5) * 14,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 62),
        (monocalcium_phosphate[1] - 7) * 62,
        1e-15,
    );
}

test "runtime grid and top-layer mapping apply only source-gated cells" {
    const cell_count = 6;
    var deposited = [_]ElementalMass{.{}} ** (cell_count * biochemical_fraction_count);
    deposited[0] = .{ .carbon_g_c = 0.5, .nitrogen_g_n = 0.05, .phosphorus_g_p = 0.005 };
    deposited[5 * biochemical_fraction_count] =
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 };
    var residue = [_]ElementalMass{.{}} ** (cell_count * biochemical_fraction_count);
    var active_carbon = [_]f64{0} ** (cell_count * biochemical_fraction_count);
    var urea = [_]f64{0} ** (cell_count * 3);
    var phosphate = [_]f64{0} ** (cell_count * 3);
    var formulation = [_]UreaFormulation{.normal_release} ** cell_count;
    var state: State = .{
        .surface_residue_g_by_cell_fraction = &residue,
        .surface_active_carbon_g_c_by_cell_fraction = &active_carbon,
        .broadcast_urea_mol_n_by_soil_layer = &urea,
        .broadcast_monocalcium_phosphate_mol_by_soil_layer = &phosphate,
        .urea_formulation_by_cell = &formulation,
        .total_dynamic_ion_input_mol = 0,
    };
    try apply(.{
        .column_count = 3,
        .row_count = 2,
        .soil_layer_capacity = 3,
        .top_soil_layer_by_cell = &.{ 2, 1, 0, 1, 2, 2 },
        .negligible_carbon_g_c_by_cell = &.{ 1, 0, 0, 0, 0, 1 },
        .deposited_organic_g_by_cell_fraction = &deposited,
        .deposited_inorganic_nitrogen_g_n_by_cell = &.{ 14, 0, 0, 0, 0, 28 },
        .deposited_inorganic_phosphorus_g_p_by_cell = &.{ 62, 0, 0, 0, 0, 124 },
        .salt_equilibrium_mode = .static,
        .conversion = .{
            .nitrogen_g_n_per_mol_n = 14,
            .phosphorus_g_p_per_mol_monocalcium_phosphate = 62,
            .ion_mol_per_mol_monocalcium_phosphate = 7,
        },
    }, &state);

    try std.testing.expectEqual(@as(f64, 0), residue[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), residue[5 * biochemical_fraction_count].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), urea[2]);
    try std.testing.expectEqual(@as(f64, 2), urea[5 * 3 + 2]);
    try std.testing.expectEqual(@as(f64, 2), phosphate[5 * 3 + 2]);
    try std.testing.expectEqual(UreaFormulation.normal_release, formulation[0]);
    try std.testing.expectEqual(UreaFormulation.fast_release, formulation[5]);
    try std.testing.expectEqual(@as(f64, 0), state.total_dynamic_ion_input_mol);
}

test "late invalid manure input leaves every state field unchanged" {
    var deposited = [_]ElementalMass{.{}} ** (2 * biochemical_fraction_count);
    deposited[0].carbon_g_c = 1;
    deposited[biochemical_fraction_count].carbon_g_c = std.math.nan(f64);
    var residue = [_]ElementalMass{.{ .carbon_g_c = 3 }} ** (2 * biochemical_fraction_count);
    var active_carbon = [_]f64{4} ** (2 * biochemical_fraction_count);
    var urea = [_]f64{ 5, 6 };
    var phosphate = [_]f64{ 7, 8 };
    var formulation = [_]UreaFormulation{ .slow_release, .normal_release };
    var state: State = .{
        .surface_residue_g_by_cell_fraction = &residue,
        .surface_active_carbon_g_c_by_cell_fraction = &active_carbon,
        .broadcast_urea_mol_n_by_soil_layer = &urea,
        .broadcast_monocalcium_phosphate_mol_by_soil_layer = &phosphate,
        .urea_formulation_by_cell = &formulation,
        .total_dynamic_ion_input_mol = 9,
    };
    try std.testing.expectError(error.InvalidManureIngressInput, apply(.{
        .column_count = 2,
        .row_count = 1,
        .soil_layer_capacity = 1,
        .top_soil_layer_by_cell = &.{ 0, 0 },
        .negligible_carbon_g_c_by_cell = &.{ 0, 0 },
        .deposited_organic_g_by_cell_fraction = &deposited,
        .deposited_inorganic_nitrogen_g_n_by_cell = &.{ 14, 14 },
        .deposited_inorganic_phosphorus_g_p_by_cell = &.{ 62, 62 },
        .salt_equilibrium_mode = .dynamic,
        .conversion = .{
            .nitrogen_g_n_per_mol_n = 14,
            .phosphorus_g_p_per_mol_monocalcium_phosphate = 62,
            .ion_mol_per_mol_monocalcium_phosphate = 7,
        },
    }, &state));
    for (residue) |mass| try std.testing.expectEqual(@as(f64, 3), mass.carbon_g_c);
    try std.testing.expectEqualSlices(f64, &.{ 4, 4, 4, 4, 4, 4, 4, 4 }, &active_carbon);
    try std.testing.expectEqualSlices(f64, &.{ 5, 6 }, &urea);
    try std.testing.expectEqualSlices(f64, &.{ 7, 8 }, &phosphate);
    try std.testing.expectEqualSlices(
        UreaFormulation,
        &.{ .slow_release, .normal_release },
        &formulation,
    );
    try std.testing.expectEqual(@as(f64, 9), state.total_dynamic_ion_input_mol);
}

test "invalid runtime extents and conversion factors fail explicitly" {
    var residue = [_]ElementalMass{.{}} ** biochemical_fraction_count;
    var active_carbon = [_]f64{0} ** biochemical_fraction_count;
    var urea = [_]f64{0};
    var phosphate = [_]f64{0};
    var formulation = [_]UreaFormulation{.normal_release};
    var state: State = .{
        .surface_residue_g_by_cell_fraction = &residue,
        .surface_active_carbon_g_c_by_cell_fraction = &active_carbon,
        .broadcast_urea_mol_n_by_soil_layer = &urea,
        .broadcast_monocalcium_phosphate_mol_by_soil_layer = &phosphate,
        .urea_formulation_by_cell = &formulation,
        .total_dynamic_ion_input_mol = 0,
    };
    const deposited = [_]ElementalMass{.{}} ** biochemical_fraction_count;
    try std.testing.expectError(error.InvalidManureIngressInput, apply(.{
        .column_count = 1,
        .row_count = 1,
        .soil_layer_capacity = 1,
        .top_soil_layer_by_cell = &.{1},
        .negligible_carbon_g_c_by_cell = &.{0},
        .deposited_organic_g_by_cell_fraction = &deposited,
        .deposited_inorganic_nitrogen_g_n_by_cell = &.{0},
        .deposited_inorganic_phosphorus_g_p_by_cell = &.{0},
        .salt_equilibrium_mode = .static,
        .conversion = .{
            .nitrogen_g_n_per_mol_n = 14,
            .phosphorus_g_p_per_mol_monocalcium_phosphate = 62,
            .ion_mol_per_mol_monocalcium_phosphate = 7,
        },
    }, &state));
    try std.testing.expectError(error.InvalidManureIngressState, apply(.{
        .column_count = 1,
        .row_count = 1,
        .soil_layer_capacity = 1,
        .top_soil_layer_by_cell = &.{0},
        .negligible_carbon_g_c_by_cell = &.{0},
        .deposited_organic_g_by_cell_fraction = &deposited,
        .deposited_inorganic_nitrogen_g_n_by_cell = &.{0},
        .deposited_inorganic_phosphorus_g_p_by_cell = &.{0},
        .salt_equilibrium_mode = .static,
        .conversion = .{
            .nitrogen_g_n_per_mol_n = 0,
            .phosphorus_g_p_per_mol_monocalcium_phosphate = 62,
            .ion_mol_per_mol_monocalcium_phosphate = 7,
        },
    }, &state));
}
