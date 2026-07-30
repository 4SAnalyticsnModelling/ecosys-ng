const std = @import("std");

pub const MultivalentActivities = struct {
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
};

pub const ActivityRoots = struct {
    aluminum_mol_per_m3_cuberoot: f64,
    iron_mol_per_m3_cuberoot: f64,
    calcium_mol_per_m3_squareroot: f64,
    magnesium_mol_per_m3_squareroot: f64,
};

pub const FloorFlags = struct {
    aluminum: bool,
    iron: bool,
    calcium: bool,
    magnesium: bool,
};

pub const Inputs = struct {
    activities: MultivalentActivities,
    minimum_activity_mol_per_m3: f64,
};

pub const Result = struct {
    roots: ActivityRoots,
    floors_applied: FloorFlags,
};

/// Direct source-order translation of SOLUTE.F lines 4746--4749.
///
/// The caller resolves surface-litter cell index `(0, NY, NX)`. These roots
/// are pure one-cell intermediates for the subsequent Gapon exchange kernel.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const activities = inputs.activities;
    const floor = inputs.minimum_activity_mol_per_m3;

    // Retain source exponents exactly: 0.333 is not replaced by 1 / 3.
    const roots: ActivityRoots = .{
        .aluminum_mol_per_m3_cuberoot = std.math.pow(
            f64,
            @max(floor, activities.aluminum_mol_per_m3),
            0.333,
        ),
        .iron_mol_per_m3_cuberoot = std.math.pow(
            f64,
            @max(floor, activities.iron_mol_per_m3),
            0.333,
        ),
        .calcium_mol_per_m3_squareroot = std.math.pow(
            f64,
            @max(floor, activities.calcium_mol_per_m3),
            0.500,
        ),
        .magnesium_mol_per_m3_squareroot = std.math.pow(
            f64,
            @max(floor, activities.magnesium_mol_per_m3),
            0.500,
        ),
    };
    const result: Result = .{
        .roots = roots,
        .floors_applied = .{
            .aluminum = activities.aluminum_mol_per_m3 < floor,
            .iron = activities.iron_mol_per_m3 < floor,
            .calcium = activities.calcium_mol_per_m3 < floor,
            .magnesium = activities.magnesium_mol_per_m3 < floor,
        },
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.minimum_activity_mol_per_m3) or
        inputs.minimum_activity_mol_per_m3 <= 0)
    {
        return error.InvalidSurfaceLitterStaticActivityRootInput;
    }
    inline for (@typeInfo(MultivalentActivities).@"struct".fields) |field| {
        const value = @field(inputs.activities, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterStaticActivityRootInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(ActivityRoots).@"struct".fields) |field| {
        const value = @field(result.roots, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.NonFiniteSurfaceLitterStaticActivityRootResult;
    }
}

fn testInputs() Inputs {
    return .{
        .activities = .{
            .aluminum_mol_per_m3 = 8,
            .iron_mol_per_m3 = 27,
            .calcium_mol_per_m3 = 16,
            .magnesium_mol_per_m3 = 25,
        },
        .minimum_activity_mol_per_m3 = 1.0e-32,
    };
}

test "SOLUTE static multivalent roots preserve every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const activity = inputs.activities;
    const floor = inputs.minimum_activity_mol_per_m3;

    try std.testing.expectEqual(
        std.math.pow(
            f64,
            @max(floor, activity.aluminum_mol_per_m3),
            0.333,
        ),
        result.roots.aluminum_mol_per_m3_cuberoot,
    );
    try std.testing.expectEqual(
        std.math.pow(
            f64,
            @max(floor, activity.iron_mol_per_m3),
            0.333,
        ),
        result.roots.iron_mol_per_m3_cuberoot,
    );
    try std.testing.expectEqual(
        std.math.pow(
            f64,
            @max(floor, activity.calcium_mol_per_m3),
            0.500,
        ),
        result.roots.calcium_mol_per_m3_squareroot,
    );
    try std.testing.expectEqual(
        std.math.pow(
            f64,
            @max(floor, activity.magnesium_mol_per_m3),
            0.500,
        ),
        result.roots.magnesium_mol_per_m3_squareroot,
    );
}

test "static multivalent activity floors are explicit" {
    var inputs = testInputs();
    inputs.activities = .{
        .aluminum_mol_per_m3 = 0,
        .iron_mol_per_m3 = 0,
        .calcium_mol_per_m3 = 0,
        .magnesium_mol_per_m3 = 0,
    };
    const result = try calculateSourceOrder(inputs);

    try std.testing.expect(result.floors_applied.aluminum);
    try std.testing.expect(result.floors_applied.iron);
    try std.testing.expect(result.floors_applied.calcium);
    try std.testing.expect(result.floors_applied.magnesium);
    try std.testing.expectEqual(
        std.math.pow(f64, inputs.minimum_activity_mol_per_m3, 0.333),
        result.roots.aluminum_mol_per_m3_cuberoot,
    );
    try std.testing.expectEqual(
        std.math.pow(f64, inputs.minimum_activity_mol_per_m3, 0.500),
        result.roots.calcium_mol_per_m3_squareroot,
    );
}

test "static aluminum and iron retain literal 0.333 exponent" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);

    try std.testing.expect(
        result.roots.aluminum_mol_per_m3_cuberoot !=
            std.math.pow(
                f64,
                inputs.activities.aluminum_mol_per_m3,
                1.0 / 3.0,
            ),
    );
    try std.testing.expect(
        result.roots.iron_mol_per_m3_cuberoot !=
            std.math.pow(
                f64,
                inputs.activities.iron_mol_per_m3,
                1.0 / 3.0,
            ),
    );
}

test "static multivalent roots reject invalid physical domains" {
    var inputs = testInputs();
    inputs.minimum_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticActivityRootInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.activities.aluminum_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticActivityRootInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.activities.calcium_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticActivityRootInput,
        calculateSourceOrder(inputs),
    );
}
