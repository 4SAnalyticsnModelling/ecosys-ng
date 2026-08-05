const std = @import("std");

pub const ExchangeInventories = struct {
    hydrogen_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    carboxyl_hydrogen_mol: f64,
};

pub const ExchangeConcentrations = struct {
    hydrogen_mol_per_megagram: f64,
    aluminum_mol_per_megagram: f64,
    iron_mol_per_megagram: f64,
    calcium_mol_per_megagram: f64,
    magnesium_mol_per_megagram: f64,
    sodium_mol_per_megagram: f64,
    potassium_mol_per_megagram: f64,
    carboxyl_hydrogen_mol_per_megagram: f64,
};

pub const MultivalentActivities = struct {
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
};

pub const MultivalentActivityRoots = struct {
    aluminum_mol_per_m3_cuberoot: f64,
    iron_mol_per_m3_cuberoot: f64,
    calcium_mol_per_m3_squareroot: f64,
    magnesium_mol_per_m3_squareroot: f64,
};

pub const FloorFlags = struct {
    hydrogen: bool,
    aluminum: bool,
    iron: bool,
    calcium: bool,
    magnesium: bool,
    sodium: bool,
    potassium: bool,
    carboxyl_hydrogen: bool,
    aluminum_activity: bool,
    iron_activity: bool,
    calcium_activity: bool,
    magnesium_activity: bool,
};

pub const Inputs = struct {
    litter_dry_mass_megagrams: f64,
    exchange_inventories: ExchangeInventories,
    multivalent_activities: MultivalentActivities,
    minimum_exchange_concentration_mol_per_megagram: f64,
    minimum_aqueous_activity_mol_per_m3: f64,
};

pub const Result = struct {
    exchange_concentrations: ExchangeConcentrations,
    activity_roots: MultivalentActivityRoots,
    floors_applied: FloorFlags,
};

/// Direct source-order translation of SOLUTE.F lines 4323--4343.
///
/// The source reaches these divisions through its wet-litter branch but does
/// not guard zero litter mass. ecosys-ng rejects that unsafe state instead of
/// allowing an infinite exchange concentration.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const mass = inputs.litter_dry_mass_megagrams;
    const floor = inputs.minimum_exchange_concentration_mol_per_megagram;
    const inventory = inputs.exchange_inventories;

    // SOLUTE.F 4332--4339. Preserve field and division/floor order.
    const unconstrained: ExchangeConcentrations = .{
        .hydrogen_mol_per_megagram = inventory.hydrogen_mol / mass,
        .aluminum_mol_per_megagram = inventory.aluminum_mol / mass,
        .iron_mol_per_megagram = inventory.iron_mol / mass,
        .calcium_mol_per_megagram = inventory.calcium_mol / mass,
        .magnesium_mol_per_megagram = inventory.magnesium_mol / mass,
        .sodium_mol_per_megagram = inventory.sodium_mol / mass,
        .potassium_mol_per_megagram = inventory.potassium_mol / mass,
        .carboxyl_hydrogen_mol_per_megagram = inventory.carboxyl_hydrogen_mol / mass,
    };
    const concentrations: ExchangeConcentrations = .{
        .hydrogen_mol_per_megagram = @max(floor, unconstrained.hydrogen_mol_per_megagram),
        .aluminum_mol_per_megagram = @max(floor, unconstrained.aluminum_mol_per_megagram),
        .iron_mol_per_megagram = @max(floor, unconstrained.iron_mol_per_megagram),
        .calcium_mol_per_megagram = @max(floor, unconstrained.calcium_mol_per_megagram),
        .magnesium_mol_per_megagram = @max(floor, unconstrained.magnesium_mol_per_megagram),
        .sodium_mol_per_megagram = @max(floor, unconstrained.sodium_mol_per_megagram),
        .potassium_mol_per_megagram = @max(floor, unconstrained.potassium_mol_per_megagram),
        .carboxyl_hydrogen_mol_per_megagram = @max(floor, unconstrained.carboxyl_hydrogen_mol_per_megagram),
    };

    const activity = inputs.multivalent_activities;
    const activity_floor = inputs.minimum_aqueous_activity_mol_per_m3;
    // SOLUTE.F 4340--4343 retains the source exponents, including 0.333.
    const roots: MultivalentActivityRoots = .{
        .aluminum_mol_per_m3_cuberoot = std.math.pow(f64, @max(activity_floor, activity.aluminum_mol_per_m3), 0.333),
        .iron_mol_per_m3_cuberoot = std.math.pow(f64, @max(activity_floor, activity.iron_mol_per_m3), 0.333),
        .calcium_mol_per_m3_squareroot = std.math.pow(f64, @max(activity_floor, activity.calcium_mol_per_m3), 0.500),
        .magnesium_mol_per_m3_squareroot = std.math.pow(f64, @max(activity_floor, activity.magnesium_mol_per_m3), 0.500),
    };
    const result: Result = .{
        .exchange_concentrations = concentrations,
        .activity_roots = roots,
        .floors_applied = .{
            .hydrogen = unconstrained.hydrogen_mol_per_megagram < floor,
            .aluminum = unconstrained.aluminum_mol_per_megagram < floor,
            .iron = unconstrained.iron_mol_per_megagram < floor,
            .calcium = unconstrained.calcium_mol_per_megagram < floor,
            .magnesium = unconstrained.magnesium_mol_per_megagram < floor,
            .sodium = unconstrained.sodium_mol_per_megagram < floor,
            .potassium = unconstrained.potassium_mol_per_megagram < floor,
            .carboxyl_hydrogen = unconstrained.carboxyl_hydrogen_mol_per_megagram < floor,
            .aluminum_activity = activity.aluminum_mol_per_m3 < activity_floor,
            .iron_activity = activity.iron_mol_per_m3 < activity_floor,
            .calcium_activity = activity.calcium_mol_per_m3 < activity_floor,
            .magnesium_activity = activity.magnesium_mol_per_m3 < activity_floor,
        },
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    if (!std.math.isFinite(inputs.litter_dry_mass_megagrams) or
        inputs.litter_dry_mass_megagrams <= 0 or
        !std.math.isFinite(inputs.minimum_exchange_concentration_mol_per_megagram) or
        inputs.minimum_exchange_concentration_mol_per_megagram <= 0 or
        !std.math.isFinite(inputs.minimum_aqueous_activity_mol_per_m3) or
        inputs.minimum_aqueous_activity_mol_per_m3 <= 0)
    {
        return error.InvalidSurfaceLitterCationExchangeInitializationInput;
    }
    inline for (@typeInfo(ExchangeInventories).@"struct".fields) |field| {
        const value = @field(inputs.exchange_inventories, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterCationExchangeInitializationInput;
    }
    inline for (@typeInfo(MultivalentActivities).@"struct".fields) |field| {
        const value = @field(inputs.multivalent_activities, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterCationExchangeInitializationInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(ExchangeConcentrations).@"struct".fields) |field| {
        const value = @field(result.exchange_concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteSurfaceLitterCationExchangeInitialization;
    }
    inline for (@typeInfo(MultivalentActivityRoots).@"struct".fields) |field| {
        const value = @field(result.activity_roots, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.NonFiniteSurfaceLitterCationExchangeInitialization;
    }
}

fn testInputs() Inputs {
    return .{
        .litter_dry_mass_megagrams = 2,
        .exchange_inventories = .{
            .hydrogen_mol = 2,
            .aluminum_mol = 4,
            .iron_mol = 6,
            .calcium_mol = 8,
            .magnesium_mol = 10,
            .sodium_mol = 12,
            .potassium_mol = 14,
            .carboxyl_hydrogen_mol = 16,
        },
        .multivalent_activities = .{
            .aluminum_mol_per_m3 = 8,
            .iron_mol_per_m3 = 27,
            .calcium_mol_per_m3 = 16,
            .magnesium_mol_per_m3 = 25,
        },
        .minimum_exchange_concentration_mol_per_megagram = 1.0e-32,
        .minimum_aqueous_activity_mol_per_m3 = 1.0e-32,
    };
}

test "SOLUTE surface exchange initialization preserves every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);
    const inventory = inputs.exchange_inventories;
    const mass = inputs.litter_dry_mass_megagrams;
    const floor = inputs.minimum_exchange_concentration_mol_per_megagram;
    const activity = inputs.multivalent_activities;
    const activity_floor = inputs.minimum_aqueous_activity_mol_per_m3;

    try std.testing.expectEqual(
        @max(floor, inventory.hydrogen_mol / mass),
        result.exchange_concentrations.hydrogen_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.aluminum_mol / mass),
        result.exchange_concentrations.aluminum_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.iron_mol / mass),
        result.exchange_concentrations.iron_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.calcium_mol / mass),
        result.exchange_concentrations.calcium_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.magnesium_mol / mass),
        result.exchange_concentrations.magnesium_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.sodium_mol / mass),
        result.exchange_concentrations.sodium_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.potassium_mol / mass),
        result.exchange_concentrations.potassium_mol_per_megagram,
    );
    try std.testing.expectEqual(
        @max(floor, inventory.carboxyl_hydrogen_mol / mass),
        result.exchange_concentrations.carboxyl_hydrogen_mol_per_megagram,
    );
    try std.testing.expectEqual(
        std.math.pow(f64, @max(activity_floor, activity.aluminum_mol_per_m3), 0.333),
        result.activity_roots.aluminum_mol_per_m3_cuberoot,
    );
    try std.testing.expectEqual(
        std.math.pow(f64, @max(activity_floor, activity.iron_mol_per_m3), 0.333),
        result.activity_roots.iron_mol_per_m3_cuberoot,
    );
    try std.testing.expectEqual(
        std.math.pow(f64, @max(activity_floor, activity.calcium_mol_per_m3), 0.500),
        result.activity_roots.calcium_mol_per_m3_squareroot,
    );
    try std.testing.expectEqual(
        std.math.pow(f64, @max(activity_floor, activity.magnesium_mol_per_m3), 0.500),
        result.activity_roots.magnesium_mol_per_m3_squareroot,
    );
}

test "surface exchange concentrations reconstruct extensive inventories" {
    const inputs = testInputs();
    const concentrations = (try calculate(inputs)).exchange_concentrations;
    inline for (@typeInfo(ExchangeInventories).@"struct".fields, 0..) |field, index| {
        const concentration_field =
            @typeInfo(ExchangeConcentrations).@"struct".fields[index];
        try std.testing.expectApproxEqAbs(
            @field(inputs.exchange_inventories, field.name),
            @field(concentrations, concentration_field.name) *
                inputs.litter_dry_mass_megagrams,
            1.0e-14,
        );
    }
}

test "surface exchange reconstruction conserves inventory as litter mass changes" {
    const inputs = testInputs();
    const initial = (try calculate(inputs)).exchange_concentrations;
    var reduced_mass = inputs;
    reduced_mass.litter_dry_mass_megagrams = 0.5;
    const concentrated = (try calculate(reduced_mass)).exchange_concentrations;

    try std.testing.expectEqual(
        initial.calcium_mol_per_megagram * inputs.litter_dry_mass_megagrams,
        concentrated.calcium_mol_per_megagram * reduced_mass.litter_dry_mass_megagrams,
    );
    try std.testing.expectEqual(
        initial.carboxyl_hydrogen_mol_per_megagram * inputs.litter_dry_mass_megagrams,
        concentrated.carboxyl_hydrogen_mol_per_megagram *
            reduced_mass.litter_dry_mass_megagrams,
    );
}

test "surface exchange and activity floors are explicit" {
    var inputs = testInputs();
    inputs.exchange_inventories = std.mem.zeroes(ExchangeInventories);
    inputs.multivalent_activities = std.mem.zeroes(MultivalentActivities);
    const result = try calculate(inputs);

    inline for (@typeInfo(FloorFlags).@"struct".fields) |field|
        try std.testing.expect(@field(result.floors_applied, field.name));
    try std.testing.expectEqual(
        inputs.minimum_exchange_concentration_mol_per_megagram,
        result.exchange_concentrations.calcium_mol_per_megagram,
    );
    try std.testing.expectEqual(
        std.math.pow(
            f64,
            inputs.minimum_aqueous_activity_mol_per_m3,
            0.333,
        ),
        result.activity_roots.aluminum_mol_per_m3_cuberoot,
    );
}

test "surface exchange initialization rejects unsafe input and overflow" {
    var inputs = testInputs();
    inputs.litter_dry_mass_megagrams = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterCationExchangeInitializationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.exchange_inventories.calcium_mol = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterCationExchangeInitializationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.multivalent_activities.aluminum_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterCationExchangeInitializationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.litter_dry_mass_megagrams = std.math.floatMin(f64);
    inputs.exchange_inventories.hydrogen_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterCationExchangeInitialization,
        calculate(inputs),
    );
}
