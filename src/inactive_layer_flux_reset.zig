const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static_concentrations,
    dynamic_transport,
};

/// REDIST lines 3767--3776 reset these values in this exact order.
pub const WaterHeatPhasePool = enum {
    micropore_water,
    water_vapor,
    wetting_front_micropore_water,
    macropore_water,
    convective_heat,
    micropore_freeze_thaw_water,
    macropore_freeze_thaw_water,
    freeze_thaw_latent_heat,
    evaporation_condensation_water,
    evaporation_condensation_latent_heat,
};

pub const water_heat_phase_pool_count: usize =
    @typeInfo(WaterHeatPhasePool).@"enum".fields.len;

/// Pool order within each runtime organic-matter class.
pub const OrganicClassPool = enum {
    micropore_dissolved_organic_carbon,
    micropore_dissolved_organic_nitrogen,
    micropore_dissolved_organic_phosphorus,
    micropore_acetate,
    macropore_dissolved_organic_carbon,
    macropore_dissolved_organic_nitrogen,
    macropore_dissolved_organic_phosphorus,
    macropore_acetate,
};

pub const organic_class_pool_count: usize =
    @typeInfo(OrganicClassPool).@"enum".fields.len;

/// Scalar dissolved-solute pools in REDIST lines 3787--3822.
pub const SolutePool = enum {
    micropore_carbon_dioxide,
    micropore_methane,
    micropore_oxygen,
    micropore_dinitrogen,
    micropore_nitrous_oxide,
    micropore_hydrogen,
    micropore_nonband_ammonium,
    micropore_nonband_ammonia,
    micropore_nonband_nitrate,
    micropore_nonband_nitrite,
    micropore_nonband_hydrogen_phosphate,
    micropore_nonband_dihydrogen_phosphate,
    micropore_band_ammonium,
    micropore_band_ammonia,
    micropore_band_nitrate,
    micropore_band_nitrite,
    micropore_band_hydrogen_phosphate,
    micropore_band_dihydrogen_phosphate,
    macropore_carbon_dioxide,
    macropore_methane,
    macropore_oxygen,
    macropore_dinitrogen,
    macropore_nitrous_oxide,
    macropore_hydrogen,
    macropore_nonband_ammonium,
    macropore_nonband_ammonia,
    macropore_nonband_nitrate,
    macropore_nonband_nitrite,
    macropore_nonband_hydrogen_phosphate,
    macropore_nonband_dihydrogen_phosphate,
    macropore_band_ammonium,
    macropore_band_ammonia,
    macropore_band_nitrate,
    macropore_band_nitrite,
    macropore_band_hydrogen_phosphate,
    macropore_band_dihydrogen_phosphate,
};

pub const solute_pool_count: usize =
    @typeInfo(SolutePool).@"enum".fields.len;

pub const GasPool = enum {
    carbon_dioxide,
    methane,
    oxygen,
    dinitrogen,
    nitrous_oxide,
    ammonia,
    hydrogen,
};

pub const gas_pool_count: usize = @typeInfo(GasPool).@"enum".fields.len;

pub const PrimaryIonPool = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
};

pub const MetalComplexPool = enum {
    aluminum_monohydroxide,
    aluminum_dihydroxide,
    aluminum_trihydroxide,
    aluminum_tetrahydroxide,
    aluminum_sulfate,
    iron_monohydroxide,
    iron_dihydroxide,
    iron_trihydroxide,
    iron_tetrahydroxide,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
};

pub const PhosphorusComplexPool = enum {
    phosphate,
    phosphoric_acid,
    iron_hydrogen_phosphate,
    iron_dihydrogen_phosphate,
    calcium_phosphate,
    calcium_hydrogen_phosphate,
    calcium_dihydrogen_phosphate,
    magnesium_hydrogen_phosphate,
};

const primary_ion_pool_count =
    @typeInfo(PrimaryIonPool).@"enum".fields.len;
const metal_complex_pool_count =
    @typeInfo(MetalComplexPool).@"enum".fields.len;
const phosphorus_complex_pool_count =
    @typeInfo(PhosphorusComplexPool).@"enum".fields.len;

/// Flattened order is primary ions, metal complexes, H-silicate, non-band P,
/// then band P.
pub const micropore_salt_pool_count: usize =
    primary_ion_pool_count +
    metal_complex_pool_count +
    1 +
    2 * phosphorus_complex_pool_count;

/// Flattened order matches micropores but omits H-silicate.
pub const macropore_salt_pool_count: usize =
    primary_ion_pool_count +
    metal_complex_pool_count +
    2 * phosphorus_complex_pool_count;

pub const Inputs = struct {
    salt_equilibrium_mode: SaltEquilibriumMode,
    organic_class_count: usize,
};

pub const State = struct {
    /// [water_heat_phase_pool]; water m3/step and heat MJ/step.
    water_heat_phase_by_pool: []f64,
    /// [organic_class][organic_class_pool], class-major, g/step.
    organic_g_per_step_by_class_pool: []f64,
    /// [solute_pool], g/step.
    solute_g_per_step_by_pool: []f64,
    /// [gas_pool], g/step.
    gas_g_per_step_by_pool: []f64,
    /// [micropore salt layout], mol/step.
    micropore_salt_mol_per_step_by_pool: []f64,
    /// [macropore salt layout], mol/step.
    macropore_salt_mol_per_step_by_pool: []f64,
};

/// Clears flux publication for a layer whose thickness is not active.
///
/// Traceability: REDIST.F lines 3767--3930, the `DLYR <= DLYRM` branch.
/// Water/heat/phase, all runtime organic classes, 36 dissolved-solute pools,
/// and seven gas pools are always reset in source order. The 50 micropore and
/// 49 macropore salt pools are reset only in dynamic salt mode, preserving the
/// nested `ISALTG != 0` branch. All dimensions, ownership, and values that will
/// be touched are validated before the first mutation.
pub fn reset(inputs: Inputs, state: State) !void {
    try validateDimensions(inputs, state);
    try validateOwnership(state);
    try validateFinite(state.water_heat_phase_by_pool);
    try validateFinite(state.organic_g_per_step_by_class_pool);
    try validateFinite(state.solute_g_per_step_by_pool);
    try validateFinite(state.gas_g_per_step_by_pool);
    if (inputs.salt_equilibrium_mode == .dynamic_transport) {
        try validateFinite(state.micropore_salt_mol_per_step_by_pool);
        try validateFinite(state.macropore_salt_mol_per_step_by_pool);
    }

    @memset(state.water_heat_phase_by_pool, 0);
    @memset(state.organic_g_per_step_by_class_pool, 0);
    @memset(state.solute_g_per_step_by_pool, 0);
    @memset(state.gas_g_per_step_by_pool, 0);
    if (inputs.salt_equilibrium_mode == .dynamic_transport) {
        @memset(state.micropore_salt_mol_per_step_by_pool, 0);
        @memset(state.macropore_salt_mol_per_step_by_pool, 0);
    }
}

fn validateDimensions(inputs: Inputs, state: State) !void {
    if (inputs.organic_class_count == 0)
        return error.InvalidInactiveLayerFluxDimensions;
    const organic_extent = std.math.mul(
        usize,
        inputs.organic_class_count,
        organic_class_pool_count,
    ) catch return error.InvalidInactiveLayerFluxDimensions;
    if (state.water_heat_phase_by_pool.len != water_heat_phase_pool_count or
        state.organic_g_per_step_by_class_pool.len != organic_extent or
        state.solute_g_per_step_by_pool.len != solute_pool_count or
        state.gas_g_per_step_by_pool.len != gas_pool_count or
        state.micropore_salt_mol_per_step_by_pool.len !=
            micropore_salt_pool_count or
        state.macropore_salt_mol_per_step_by_pool.len !=
            macropore_salt_pool_count)
    {
        return error.InactiveLayerFluxDimensionMismatch;
    }
}

fn validateOwnership(state: State) !void {
    const slices = [_][]const f64{
        state.water_heat_phase_by_pool,
        state.organic_g_per_step_by_class_pool,
        state.solute_g_per_step_by_pool,
        state.gas_g_per_step_by_pool,
        state.micropore_salt_mol_per_step_by_pool,
        state.macropore_salt_mol_per_step_by_pool,
    };
    for (slices, 0..) |left, left_index|
        for (slices[left_index + 1 ..]) |right|
            if (overlap(left, right))
                return error.InactiveLayerFluxStorageOverlap;
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteInactiveLayerFluxState;
}

fn overlap(left: []const f64, right: []const f64) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(f64);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn expectAll(values: []const f64, expected: f64) !void {
    for (values) |value| try std.testing.expectEqual(expected, value);
}

const TestStorage = struct {
    water: [water_heat_phase_pool_count]f64,
    organic: [3 * organic_class_pool_count]f64,
    solute: [solute_pool_count]f64,
    gas: [gas_pool_count]f64,
    micropore_salt: [micropore_salt_pool_count]f64,
    macropore_salt: [macropore_salt_pool_count]f64,

    fn init(value: f64) TestStorage {
        return .{
            .water = [_]f64{value} ** water_heat_phase_pool_count,
            .organic = [_]f64{value} ** (3 * organic_class_pool_count),
            .solute = [_]f64{value} ** solute_pool_count,
            .gas = [_]f64{value} ** gas_pool_count,
            .micropore_salt = [_]f64{value} ** micropore_salt_pool_count,
            .macropore_salt = [_]f64{value} ** macropore_salt_pool_count,
        };
    }

    fn state(self: *TestStorage) State {
        return .{
            .water_heat_phase_by_pool = &self.water,
            .organic_g_per_step_by_class_pool = &self.organic,
            .solute_g_per_step_by_pool = &self.solute,
            .gas_g_per_step_by_pool = &self.gas,
            .micropore_salt_mol_per_step_by_pool = &self.micropore_salt,
            .macropore_salt_mol_per_step_by_pool = &self.macropore_salt,
        };
    }
};

test "source layouts retain every inactive-layer reset position" {
    try std.testing.expectEqual(@as(usize, 10), water_heat_phase_pool_count);
    try std.testing.expectEqual(@as(usize, 8), organic_class_pool_count);
    try std.testing.expectEqual(@as(usize, 36), solute_pool_count);
    try std.testing.expectEqual(@as(usize, 7), gas_pool_count);
    try std.testing.expectEqual(@as(usize, 50), micropore_salt_pool_count);
    try std.testing.expectEqual(@as(usize, 49), macropore_salt_pool_count);
    try std.testing.expectEqual(
        @as(usize, 0),
        @intFromEnum(WaterHeatPhasePool.micropore_water),
    );
    try std.testing.expectEqual(
        @as(usize, 9),
        @intFromEnum(
            WaterHeatPhasePool.evaporation_condensation_latent_heat,
        ),
    );
    try std.testing.expectEqual(
        @as(usize, 35),
        @intFromEnum(SolutePool.macropore_band_dihydrogen_phosphate),
    );
}

test "dynamic mode resets every runtime flux inventory" {
    var storage = TestStorage.init(-7);
    try reset(.{
        .salt_equilibrium_mode = .dynamic_transport,
        .organic_class_count = 3,
    }, storage.state());
    try expectAll(&storage.water, 0);
    try expectAll(&storage.organic, 0);
    try expectAll(&storage.solute, 0);
    try expectAll(&storage.gas, 0);
    try expectAll(&storage.micropore_salt, 0);
    try expectAll(&storage.macropore_salt, 0);
}

test "static salt mode resets other fluxes and retains salt state" {
    var storage = TestStorage.init(11);
    try reset(.{
        .salt_equilibrium_mode = .static_concentrations,
        .organic_class_count = 3,
    }, storage.state());
    try expectAll(&storage.water, 0);
    try expectAll(&storage.organic, 0);
    try expectAll(&storage.solute, 0);
    try expectAll(&storage.gas, 0);
    try expectAll(&storage.micropore_salt, 11);
    try expectAll(&storage.macropore_salt, 11);
}

test "invalid dimension and nonfinite state fail before any reset" {
    var storage = TestStorage.init(5);
    var state = storage.state();
    state.organic_g_per_step_by_class_pool =
        state.organic_g_per_step_by_class_pool[0..23];
    try std.testing.expectError(
        error.InactiveLayerFluxDimensionMismatch,
        reset(.{
            .salt_equilibrium_mode = .dynamic_transport,
            .organic_class_count = 3,
        }, state),
    );
    try expectAll(&storage.water, 5);
    try expectAll(&storage.organic, 5);

    state = storage.state();
    storage.gas[gas_pool_count - 1] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteInactiveLayerFluxState,
        reset(.{
            .salt_equilibrium_mode = .dynamic_transport,
            .organic_class_count = 3,
        }, state),
    );
    try expectAll(&storage.water, 5);
    try expectAll(&storage.organic, 5);
}

test "overlapping storage fails before mutation including static salts" {
    var storage = TestStorage.init(13);
    var state = storage.state();
    state.water_heat_phase_by_pool = storage.micropore_salt[0..10];
    try std.testing.expectError(
        error.InactiveLayerFluxStorageOverlap,
        reset(.{
            .salt_equilibrium_mode = .static_concentrations,
            .organic_class_count = 3,
        }, state),
    );
    try expectAll(&storage.organic, 13);
    try expectAll(&storage.micropore_salt, 13);
}
