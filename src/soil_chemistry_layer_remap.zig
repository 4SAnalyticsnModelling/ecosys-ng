const std = @import("std");
const chemistry_module = @import("solute_chemistry_state.zig");
const phosphate = @import("solute_phosphate_network.zig");
const cation = @import("solute_cation_exchange.zig");
const geochemistry = @import("solute_geochemistry_network.zig");
const aqueous = @import("solute_aqueous_network.zig");

pub const ZoneWaterVolumes = struct {
    shared_m3: f64,
    ammonium_non_band_m3: f64,
    ammonium_band_m3: f64,
    nitrate_non_band_m3: f64,
    nitrate_band_m3: f64,
    phosphate_non_band_m3: f64,
    phosphate_band_m3: f64,
};

/// REDIST ponding for soluble N/P and the complete optional dynamic-salt
/// network. Band pools retain the source gate: no destination band carrier
/// means neither destination addition nor source removal.
pub fn transferAqueousLayerFraction(
    chemistry: *chemistry_module.State,
    source: usize,
    destination: usize,
    source_water: ZoneWaterVolumes,
    destination_water: ZoneWaterVolumes,
    source_water_after: ZoneWaterVolumes,
    destination_water_after: ZoneWaterVolumes,
    dynamic_salts: bool,
    fraction: f64,
) !void {
    if (source >= chemistry.cell_count or destination >= chemistry.cell_count or source == destination) return error.ChemistryLayerRemapIndexOutOfBounds;
    try validateZoneWater(source_water);
    try validateZoneWater(destination_water);
    try validateZoneWater(source_water_after);
    try validateZoneWater(destination_water_after);
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidChemistryLayerRemapInput;
    var next_source_aqueous = chemistry.aqueous[source];
    var next_destination_aqueous = chemistry.aqueous[destination];
    inline for (@typeInfo(aqueous.State).@"struct".fields) |field| {
        if ((comptime isBasePondedAqueousField(field.name)) or dynamic_salts) {
            const source_scale = aqueousFieldVolume(field.name, source_water);
            const destination_scale = aqueousFieldVolume(field.name, destination_water);
            const source_scale_after = aqueousFieldVolume(field.name, source_water_after);
            const destination_scale_after = aqueousFieldVolume(field.name, destination_water_after);
            const is_band = comptime std.mem.endsWith(u8, field.name, "_band");
            if (!is_band or destination_scale > 0) {
                const next = try transferConcentration(@field(next_source_aqueous, field.name), @field(next_destination_aqueous, field.name), source_scale, destination_scale, source_scale_after, destination_scale_after, fraction);
                @field(next_source_aqueous, field.name) = next.source;
                @field(next_destination_aqueous, field.name) = next.destination;
            }
        }
    }
    const next_reaction_water = try transferConcentration(
        chemistry.water_mol_per_m3[source],
        chemistry.water_mol_per_m3[destination],
        source_water.shared_m3,
        destination_water.shared_m3,
        source_water_after.shared_m3,
        destination_water_after.shared_m3,
        fraction,
    );
    var next_source_non_band = chemistry.non_band_phosphate[source];
    var next_destination_non_band = chemistry.non_band_phosphate[destination];
    var next_source_band = chemistry.band_phosphate[source];
    var next_destination_band = chemistry.band_phosphate[destination];
    inline for (@typeInfo(phosphate.State).@"struct".fields) |field| if ((comptime isBaseAqueousPhosphateField(field.name)) or (dynamic_salts and (comptime isAqueousPhosphateField(field.name)))) {
        const non_band = try transferConcentration(@field(next_source_non_band, field.name), @field(next_destination_non_band, field.name), source_water.phosphate_non_band_m3, destination_water.phosphate_non_band_m3, source_water_after.phosphate_non_band_m3, destination_water_after.phosphate_non_band_m3, fraction);
        @field(next_source_non_band, field.name) = non_band.source;
        @field(next_destination_non_band, field.name) = non_band.destination;
        if (destination_water.phosphate_band_m3 > 0) {
            const band = try transferConcentration(@field(next_source_band, field.name), @field(next_destination_band, field.name), source_water.phosphate_band_m3, destination_water.phosphate_band_m3, source_water_after.phosphate_band_m3, destination_water_after.phosphate_band_m3, fraction);
            @field(next_source_band, field.name) = band.source;
            @field(next_destination_band, field.name) = band.destination;
        }
    };
    chemistry.aqueous[source] = next_source_aqueous;
    chemistry.aqueous[destination] = next_destination_aqueous;
    chemistry.water_mol_per_m3[source] = next_reaction_water.source;
    chemistry.water_mol_per_m3[destination] = next_reaction_water.destination;
    chemistry.non_band_phosphate[source] = next_source_non_band;
    chemistry.non_band_phosphate[destination] = next_destination_non_band;
    chemistry.band_phosphate[source] = next_source_band;
    chemistry.band_phosphate[destination] = next_destination_band;
}

pub fn validateAqueousLayerFraction(
    chemistry: *const chemistry_module.State,
    source: usize,
    destination: usize,
    source_water: ZoneWaterVolumes,
    destination_water: ZoneWaterVolumes,
    source_water_after: ZoneWaterVolumes,
    destination_water_after: ZoneWaterVolumes,
    dynamic_salts: bool,
    fraction: f64,
) !void {
    if (source >= chemistry.cell_count or destination >= chemistry.cell_count) return error.ChemistryLayerRemapIndexOutOfBounds;
    var view = twoCellView(chemistry, source, destination);
    var state = view.bindState();
    try transferAqueousLayerFraction(&state, 0, 1, source_water, destination_water, source_water_after, destination_water_after, dynamic_salts, fraction);
}

/// REDIST ponding of adsorbed cations, carboxyl H, phosphate surfaces and
/// precipitates, and geochemical solids. Native concentration units are
/// converted to extensive amounts on their exact soil-mass or water-volume
/// basis before the conservative transfer.
pub fn transferSolidLayerFraction(
    chemistry: *chemistry_module.State,
    source: usize,
    destination: usize,
    source_soil_mass_Mg: f64,
    destination_soil_mass_Mg: f64,
    source_water_m3: f64,
    destination_water_m3: f64,
    source_soil_mass_after_Mg: f64,
    destination_soil_mass_after_Mg: f64,
    source_water_after_m3: f64,
    destination_water_after_m3: f64,
    fraction: f64,
) !void {
    if (source >= chemistry.cell_count or destination >= chemistry.cell_count or source == destination) return error.ChemistryLayerRemapIndexOutOfBounds;
    inline for (.{ source_soil_mass_Mg, destination_soil_mass_Mg, source_water_m3, destination_water_m3, source_soil_mass_after_Mg, destination_soil_mass_after_Mg, source_water_after_m3, destination_water_after_m3, fraction }) |value|
        if (!std.math.isFinite(value)) return error.InvalidChemistryLayerRemapInput;
    if (source_soil_mass_Mg <= 0 or destination_soil_mass_Mg <= 0 or source_water_m3 < 0 or destination_water_m3 < 0 or source_soil_mass_after_Mg < 0 or destination_soil_mass_after_Mg <= 0 or source_water_after_m3 < 0 or destination_water_after_m3 < 0 or fraction < 0 or fraction > 1) return error.InvalidChemistryLayerRemapInput;

    var next_source_cations = chemistry.cation_exchange_mol_per_Mg[source];
    var next_destination_cations = chemistry.cation_exchange_mol_per_Mg[destination];
    inline for (@typeInfo(cation.Cations).@"struct".fields) |field| {
        const next = try transferConcentration(@field(next_source_cations, field.name), @field(next_destination_cations, field.name), source_soil_mass_Mg, destination_soil_mass_Mg, source_soil_mass_after_Mg, destination_soil_mass_after_Mg, fraction);
        @field(next_source_cations, field.name) = next.source;
        @field(next_destination_cations, field.name) = next.destination;
    }
    const next_carboxyl = try transferConcentration(chemistry.carboxyl_bound_hydrogen_mol_per_Mg[source], chemistry.carboxyl_bound_hydrogen_mol_per_Mg[destination], source_soil_mass_Mg, destination_soil_mass_Mg, source_soil_mass_after_Mg, destination_soil_mass_after_Mg, fraction);

    var next_source_non_band = chemistry.non_band_phosphate[source];
    var next_destination_non_band = chemistry.non_band_phosphate[destination];
    var next_source_band = chemistry.band_phosphate[source];
    var next_destination_band = chemistry.band_phosphate[destination];
    inline for (.{
        .{ &next_source_non_band, &next_destination_non_band },
        .{ &next_source_band, &next_destination_band },
    }) |zones| inline for (@typeInfo(phosphate.State).@"struct".fields) |field| {
        if (comptime isPondedPhosphateField(field.name)) {
            const source_scale = if (comptime std.mem.endsWith(u8, field.name, "_per_Mg")) source_soil_mass_Mg else source_water_m3;
            const destination_scale = if (comptime std.mem.endsWith(u8, field.name, "_per_Mg")) destination_soil_mass_Mg else destination_water_m3;
            const source_scale_after = if (comptime std.mem.endsWith(u8, field.name, "_per_Mg")) source_soil_mass_after_Mg else source_water_after_m3;
            const destination_scale_after = if (comptime std.mem.endsWith(u8, field.name, "_per_Mg")) destination_soil_mass_after_Mg else destination_water_after_m3;
            const next = try transferConcentration(@field(zones[0].*, field.name), @field(zones[1].*, field.name), source_scale, destination_scale, source_scale_after, destination_scale_after, fraction);
            @field(zones[0].*, field.name) = next.source;
            @field(zones[1].*, field.name) = next.destination;
        }
    };

    var next_source_solids = chemistry.geochemistry_solids[source];
    var next_destination_solids = chemistry.geochemistry_solids[destination];
    inline for (@typeInfo(geochemistry.SolidState).@"struct".fields) |field| {
        const next = try transferConcentration(@field(next_source_solids, field.name), @field(next_destination_solids, field.name), source_water_m3, destination_water_m3, source_water_after_m3, destination_water_after_m3, fraction);
        @field(next_source_solids, field.name) = next.source;
        @field(next_destination_solids, field.name) = next.destination;
    }

    chemistry.cation_exchange_mol_per_Mg[source] = next_source_cations;
    chemistry.cation_exchange_mol_per_Mg[destination] = next_destination_cations;
    chemistry.carboxyl_bound_hydrogen_mol_per_Mg[source] = next_carboxyl.source;
    chemistry.carboxyl_bound_hydrogen_mol_per_Mg[destination] = next_carboxyl.destination;
    chemistry.non_band_phosphate[source] = next_source_non_band;
    chemistry.non_band_phosphate[destination] = next_destination_non_band;
    chemistry.band_phosphate[source] = next_source_band;
    chemistry.band_phosphate[destination] = next_destination_band;
    chemistry.geochemistry_solids[source] = next_source_solids;
    chemistry.geochemistry_solids[destination] = next_destination_solids;
}

pub fn validateSolidLayerFraction(
    chemistry: *const chemistry_module.State,
    source: usize,
    destination: usize,
    source_soil_mass_Mg: f64,
    destination_soil_mass_Mg: f64,
    source_water_m3: f64,
    destination_water_m3: f64,
    source_soil_mass_after_Mg: f64,
    destination_soil_mass_after_Mg: f64,
    source_water_after_m3: f64,
    destination_water_after_m3: f64,
    fraction: f64,
) !void {
    if (source >= chemistry.cell_count or destination >= chemistry.cell_count) return error.ChemistryLayerRemapIndexOutOfBounds;
    var view = twoCellView(chemistry, source, destination);
    var state = view.bindState();
    try transferSolidLayerFraction(&state, 0, 1, source_soil_mass_Mg, destination_soil_mass_Mg, source_water_m3, destination_water_m3, source_soil_mass_after_Mg, destination_soil_mass_after_Mg, source_water_after_m3, destination_water_after_m3, fraction);
}

const TwoCellView = struct {
    aqueous_values: [2]aqueous.State,
    non_band_values: [2]phosphate.State,
    band_values: [2]phosphate.State,
    water_values: [2]f64,
    cation_values: [2]cation.Cations,
    carboxyl_values: [2]f64,
    solid_values: [2]geochemistry.SolidState,

    fn bindState(self: *TwoCellView) chemistry_module.State {
        return .{
            .allocator = undefined,
            .cell_count = 2,
            .aqueous = &self.aqueous_values,
            .non_band_phosphate = &self.non_band_values,
            .band_phosphate = &self.band_values,
            .water_mol_per_m3 = &self.water_values,
            .cation_exchange_mol_per_Mg = &self.cation_values,
            .carboxyl_bound_hydrogen_mol_per_Mg = &self.carboxyl_values,
            .geochemistry_solids = &self.solid_values,
        };
    }
};

fn twoCellView(chemistry: *const chemistry_module.State, source: usize, destination: usize) TwoCellView {
    var result: TwoCellView = undefined;
    result.aqueous_values = .{ chemistry.aqueous[source], chemistry.aqueous[destination] };
    result.non_band_values = .{ chemistry.non_band_phosphate[source], chemistry.non_band_phosphate[destination] };
    result.band_values = .{ chemistry.band_phosphate[source], chemistry.band_phosphate[destination] };
    result.water_values = .{ chemistry.water_mol_per_m3[source], chemistry.water_mol_per_m3[destination] };
    result.cation_values = .{ chemistry.cation_exchange_mol_per_Mg[source], chemistry.cation_exchange_mol_per_Mg[destination] };
    result.carboxyl_values = .{ chemistry.carboxyl_bound_hydrogen_mol_per_Mg[source], chemistry.carboxyl_bound_hydrogen_mol_per_Mg[destination] };
    result.solid_values = .{ chemistry.geochemistry_solids[source], chemistry.geochemistry_solids[destination] };
    return result;
}

const ConcentrationPair = struct { source: f64, destination: f64 };

fn transferConcentration(source: f64, destination: f64, source_scale: f64, destination_scale: f64, source_scale_after: f64, destination_scale_after: f64, fraction: f64) !ConcentrationPair {
    inline for (.{ source, destination, source_scale, destination_scale, source_scale_after, destination_scale_after }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryLayerRemapState;
    const source_amount = source * source_scale;
    const destination_amount = destination * destination_scale;
    const moved = fraction * source_amount;
    const next_source_amount = source_amount - moved;
    const next_destination_amount = destination_amount + moved;
    inline for (.{ source_amount, destination_amount, moved, next_source_amount, next_destination_amount }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryLayerRemapState;
    return .{
        .source = try concentration(next_source_amount, source_scale_after),
        .destination = try concentration(next_destination_amount, destination_scale_after),
    };
}

fn concentration(amount: f64, scale: f64) !f64 {
    if (scale > 0) return amount / scale;
    if (amount == 0) return 0;
    return error.ChemistryLayerRemapRequiresRecipientVolume;
}

fn isPondedPhosphateField(comptime name: []const u8) bool {
    @setEvalBranchQuota(10_000);
    return std.mem.endsWith(u8, name, "_per_Mg") or std.mem.indexOf(u8, name, "_solid_mol_per_m3") != null;
}

fn isAqueousPhosphateField(comptime name: []const u8) bool {
    @setEvalBranchQuota(10_000);
    return std.mem.startsWith(u8, name, "dissolved_") or std.mem.indexOf(u8, name, "_pair_") != null;
}

fn isBaseAqueousPhosphateField(comptime name: []const u8) bool {
    return std.mem.eql(u8, name, "dissolved_hpo4_mol_p_per_m3") or std.mem.eql(u8, name, "dissolved_h2po4_mol_p_per_m3");
}

fn isBasePondedAqueousField(comptime name: []const u8) bool {
    @setEvalBranchQuota(10_000);
    if (std.mem.startsWith(u8, name, "ammonium_") or std.mem.startsWith(u8, name, "ammonia_") or std.mem.startsWith(u8, name, "nitrate_")) return true;
    if (std.mem.eql(u8, name, "hydrogen") or std.mem.eql(u8, name, "hydroxide") or std.mem.eql(u8, name, "aluminum") or std.mem.eql(u8, name, "iron") or std.mem.eql(u8, name, "calcium") or std.mem.eql(u8, name, "magnesium") or std.mem.eql(u8, name, "sodium") or std.mem.eql(u8, name, "potassium")) return true;
    return false;
}

fn aqueousFieldVolume(comptime name: []const u8, volumes: ZoneWaterVolumes) f64 {
    if (std.mem.startsWith(u8, name, "ammonium_non_band") or std.mem.startsWith(u8, name, "ammonia_non_band")) return volumes.ammonium_non_band_m3;
    if (std.mem.startsWith(u8, name, "ammonium_band") or std.mem.startsWith(u8, name, "ammonia_band")) return volumes.ammonium_band_m3;
    if (std.mem.startsWith(u8, name, "nitrate_non_band")) return volumes.nitrate_non_band_m3;
    if (std.mem.startsWith(u8, name, "nitrate_band")) return volumes.nitrate_band_m3;
    return volumes.shared_m3;
}

fn validateZoneWater(volumes: ZoneWaterVolumes) !void {
    inline for (@typeInfo(ZoneWaterVolumes).@"struct".fields) |field| {
        const value = @field(volumes, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryLayerRemapInput;
    }
}

test "REDIST solid chemistry remap conserves native amounts across unequal layer bases" {
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    chemistry.cation_exchange_mol_per_Mg[0].calcium = 2;
    chemistry.carboxyl_bound_hydrogen_mol_per_Mg[0] = 3;
    chemistry.non_band_phosphate[0].adsorbed_hpo4_mol_p_per_Mg = 4;
    chemistry.band_phosphate[0].aluminum_phosphate_solid_mol_per_m3 = 5;
    chemistry.geochemistry_solids[0].potassium_ground_silicate_mol_per_m3 = 6;
    try transferSolidLayerFraction(&chemistry, 0, 1, 10, 20, 2, 4, 10, 20, 2, 4, 0.25);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), chemistry.cation_exchange_mol_per_Mg[0].calcium, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), chemistry.cation_exchange_mol_per_Mg[1].calcium, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), chemistry.non_band_phosphate[1].adsorbed_hpo4_mol_p_per_Mg, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.625), chemistry.band_phosphate[1].aluminum_phosphate_solid_mol_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), chemistry.geochemistry_solids[1].potassium_ground_silicate_mol_per_m3, 1e-14);
}

test "REDIST solid chemistry remap rejects dry recipient before mutation" {
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    chemistry.cation_exchange_mol_per_Mg[0].calcium = 2;
    chemistry.geochemistry_solids[0].calcite_solid_mol_per_m3 = 3;
    try std.testing.expectError(error.ChemistryLayerRemapRequiresRecipientVolume, transferSolidLayerFraction(&chemistry, 0, 1, 10, 10, 1, 0, 10, 10, 1, 0, 0.5));
    try std.testing.expectEqual(@as(f64, 2), chemistry.cation_exchange_mol_per_Mg[0].calcium);
    try std.testing.expectEqual(@as(f64, 0), chemistry.cation_exchange_mol_per_Mg[1].calcium);
}

test "REDIST aqueous chemistry respects independent runtime band carriers" {
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    chemistry.aqueous[0].calcium = 4;
    chemistry.aqueous[0].ammonium_band = 8;
    chemistry.aqueous[0].nitrate_band = 10;
    chemistry.water_mol_per_m3[0] = 16;
    chemistry.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 = 12;
    chemistry.band_phosphate[0].dissolved_hpo4_mol_p_per_m3 = 14;
    const source: ZoneWaterVolumes = .{ .shared_m3 = 2, .ammonium_non_band_m3 = 2, .ammonium_band_m3 = 1, .nitrate_non_band_m3 = 2, .nitrate_band_m3 = 1, .phosphate_non_band_m3 = 2, .phosphate_band_m3 = 1 };
    const destination: ZoneWaterVolumes = .{ .shared_m3 = 4, .ammonium_non_band_m3 = 4, .ammonium_band_m3 = 2, .nitrate_non_band_m3 = 4, .nitrate_band_m3 = 0, .phosphate_non_band_m3 = 4, .phosphate_band_m3 = 2 };
    try transferAqueousLayerFraction(&chemistry, 0, 1, source, destination, source, destination, true, 0.25);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), chemistry.aqueous[1].calcium, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1), chemistry.aqueous[1].ammonium_band, 1e-14);
    try std.testing.expectEqual(@as(f64, 10), chemistry.aqueous[0].nitrate_band);
    try std.testing.expectEqual(@as(f64, 0), chemistry.aqueous[1].nitrate_band);
    try std.testing.expectApproxEqAbs(
        @as(f64, 12),
        chemistry.water_mol_per_m3[0],
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 2),
        chemistry.water_mol_per_m3[1],
        1e-14,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), chemistry.non_band_phosphate[1].dissolved_h2po4_mol_p_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), chemistry.band_phosphate[1].dissolved_hpo4_mol_p_per_m3, 1e-14);
}

test "REDIST aqueous concentration uses pre-transfer amount and post-transfer carrier" {
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    chemistry.aqueous[0].calcium = 4;
    const source: ZoneWaterVolumes = .{ .shared_m3 = 2, .ammonium_non_band_m3 = 2, .ammonium_band_m3 = 0, .nitrate_non_band_m3 = 2, .nitrate_band_m3 = 0, .phosphate_non_band_m3 = 2, .phosphate_band_m3 = 0 };
    const destination: ZoneWaterVolumes = .{ .shared_m3 = 4, .ammonium_non_band_m3 = 4, .ammonium_band_m3 = 0, .nitrate_non_band_m3 = 4, .nitrate_band_m3 = 0, .phosphate_non_band_m3 = 4, .phosphate_band_m3 = 0 };
    const source_after: ZoneWaterVolumes = .{ .shared_m3 = 1.5, .ammonium_non_band_m3 = 1.5, .ammonium_band_m3 = 0, .nitrate_non_band_m3 = 1.5, .nitrate_band_m3 = 0, .phosphate_non_band_m3 = 1.5, .phosphate_band_m3 = 0 };
    const destination_after: ZoneWaterVolumes = .{ .shared_m3 = 4.5, .ammonium_non_band_m3 = 4.5, .ammonium_band_m3 = 0, .nitrate_non_band_m3 = 4.5, .nitrate_band_m3 = 0, .phosphate_non_band_m3 = 4.5, .phosphate_band_m3 = 0 };
    try transferAqueousLayerFraction(&chemistry, 0, 1, source, destination, source_after, destination_after, false, 0.25);
    try std.testing.expectApproxEqAbs(@as(f64, 4), chemistry.aqueous[0].calcium, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 4.5), chemistry.aqueous[1].calcium, 1e-14);
}
