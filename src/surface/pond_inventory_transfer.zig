const std = @import("std");
const organic_module = @import("../soil/organic/initialization.zig");
const gas_module = @import("../soil/gas/transport.zig");
const surface_fertilizer_module = @import("litter_fertilizer.zig");
const soil_fertilizer_module = @import("../management/fertilizer_nitrogen_inventory.zig");
const mineral_fertilizer_module = @import("../management/mineral_fertilizer_inventory.zig");

pub const Owners = struct {
    surface_organic: *organic_module.State,
    soil_organic: *organic_module.State,
    surface_gas: *gas_module.State,
    soil_gas: *gas_module.State,
    surface_nitrogen_fertilizer: *surface_fertilizer_module.State,
    soil_nitrogen_fertilizer: *soil_fertilizer_module.State,
    mineral_fertilizer: *mineral_fertilizer_module.State,
};

pub const Inputs = struct {
    cell: usize,
    destination_soil_layer: usize,
    fraction: f64,
};

/// Cross-domain portion of REDIST L0=0 pond transfer. All surface owners are
/// validated before any soil or surface inventory changes.
pub fn transferSurfaceFractionToSoil(owners: Owners, inputs: Inputs) !void {
    const destination = try validate(owners, inputs);
    transferOrganic(owners.surface_organic, inputs.cell, owners.soil_organic, destination, inputs.fraction);
    transferGas(owners.surface_gas, inputs.cell, owners.soil_gas, destination, inputs.fraction);
    transferNitrogenFertilizer(owners.surface_nitrogen_fertilizer, inputs.cell, owners.soil_nitrogen_fertilizer, destination, inputs.fraction);
    transferMineralFertilizer(owners.mineral_fertilizer, inputs.cell, destination, inputs.fraction);
}

pub fn validateSurfaceFractionToSoil(owners: Owners, inputs: Inputs) !void {
    _ = try validate(owners, inputs);
}

/// REDIST line-333 settling subset represented by separated runtime owners.
/// Dissolved organic matter and gases are intentionally excluded.
pub fn validateParticulateFractionToSoil(owners: Owners, inputs: Inputs) !void {
    const destination = try validateIndicesAndFraction(owners, inputs);
    try validateParticulateOrganic(owners.surface_organic, inputs.cell, owners.soil_organic, destination, inputs.fraction);
    try validateNitrogenFertilizer(owners.surface_nitrogen_fertilizer.cells[inputs.cell], owners.soil_nitrogen_fertilizer.soil[destination], inputs.fraction);
    try validateStructPair(mineral_fertilizer_module.Inventory, owners.mineral_fertilizer.surface[inputs.cell], owners.mineral_fertilizer.soil[destination], inputs.fraction);
}

pub fn transferParticulateFractionToSoil(owners: Owners, inputs: Inputs) !void {
    try validateParticulateFractionToSoil(owners, inputs);
    const destination = inputs.cell * owners.soil_nitrogen_fertilizer.layer_capacity + inputs.destination_soil_layer;
    transferParticulateOrganic(owners.surface_organic, inputs.cell, owners.soil_organic, destination, inputs.fraction);
    transferNitrogenFertilizer(owners.surface_nitrogen_fertilizer, inputs.cell, owners.soil_nitrogen_fertilizer, destination, inputs.fraction);
    transferMineralFertilizer(owners.mineral_fertilizer, inputs.cell, destination, inputs.fraction);
}

fn validate(owners: Owners, inputs: Inputs) !usize {
    const destination = try validateIndicesAndFraction(owners, inputs);
    try validateOrganic(owners.surface_organic, inputs.cell, owners.soil_organic, destination, inputs.fraction);
    try validateGas(owners.surface_gas, inputs.cell, owners.soil_gas, destination, inputs.fraction);
    try validateNitrogenFertilizer(owners.surface_nitrogen_fertilizer.cells[inputs.cell], owners.soil_nitrogen_fertilizer.soil[destination], inputs.fraction);
    try validateStructPair(mineral_fertilizer_module.Inventory, owners.mineral_fertilizer.surface[inputs.cell], owners.mineral_fertilizer.soil[destination], inputs.fraction);
    return destination;
}

fn validateIndicesAndFraction(owners: Owners, inputs: Inputs) !usize {
    if (!std.math.isFinite(inputs.fraction) or inputs.fraction < 0 or inputs.fraction > 1) return error.InvalidSurfacePondTransferFraction;
    if (inputs.cell >= owners.surface_organic.layer_count or inputs.cell >= owners.surface_gas.cell_count or inputs.cell >= owners.surface_nitrogen_fertilizer.cells.len or inputs.cell >= owners.mineral_fertilizer.cell_count or inputs.destination_soil_layer >= owners.soil_nitrogen_fertilizer.layer_capacity) return error.SurfacePondTransferIndexOutOfBounds;
    const destination = inputs.cell * owners.soil_nitrogen_fertilizer.layer_capacity + inputs.destination_soil_layer;
    if (destination >= owners.soil_organic.layer_count or destination >= owners.soil_gas.cell_count or destination >= owners.soil_nitrogen_fertilizer.soil.len or destination >= owners.mineral_fertilizer.soil.len) return error.SurfacePondTransferDimensionMismatch;
    return destination;
}

fn validateParticulateOrganic(source: *const organic_module.State, source_layer: usize, destination: *const organic_module.State, destination_layer: usize, fraction: f64) !void {
    inline for (.{
        .{ "microbial", organic_module.microbial_substrate_count * organic_module.microbial_population_count * organic_module.kinetic_fraction_count },
        .{ "residue", organic_module.substrate_count * organic_module.residue_fraction_count },
        .{ "adsorbed", organic_module.substrate_count },
        .{ "structural", organic_module.substrate_count * organic_module.structural_fraction_count },
    }) |descriptor| {
        const source_values = @field(source, descriptor[0]);
        const destination_values = @field(destination, descriptor[0]);
        for (0..descriptor[1]) |offset| try validateStructPair(organic_module.ElementPool, source_values[source_layer * descriptor[1] + offset], destination_values[destination_layer * descriptor[1] + offset], fraction);
    }
    for (0..organic_module.substrate_count) |offset| try validateNumberPair(source.adsorbed_acetate_carbon_g_c[source_layer * organic_module.substrate_count + offset], destination.adsorbed_acetate_carbon_g_c[destination_layer * organic_module.substrate_count + offset], fraction);
    for (0..organic_module.substrate_count * organic_module.structural_fraction_count) |offset| try validateNumberPair(source.colonized_structural_carbon_g_c[source_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset], destination.colonized_structural_carbon_g_c[destination_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset], fraction);
}

fn transferParticulateOrganic(source: *organic_module.State, source_layer: usize, destination: *organic_module.State, destination_layer: usize, fraction: f64) void {
    inline for (.{
        .{ "microbial", organic_module.microbial_substrate_count * organic_module.microbial_population_count * organic_module.kinetic_fraction_count },
        .{ "residue", organic_module.substrate_count * organic_module.residue_fraction_count },
        .{ "adsorbed", organic_module.substrate_count },
        .{ "structural", organic_module.substrate_count * organic_module.structural_fraction_count },
    }) |descriptor| {
        const source_values = @field(source, descriptor[0]);
        const destination_values = @field(destination, descriptor[0]);
        for (0..descriptor[1]) |offset| transferStruct(organic_module.ElementPool, &source_values[source_layer * descriptor[1] + offset], &destination_values[destination_layer * descriptor[1] + offset], fraction);
    }
    for (0..organic_module.substrate_count) |offset| transferNumber(&source.adsorbed_acetate_carbon_g_c[source_layer * organic_module.substrate_count + offset], &destination.adsorbed_acetate_carbon_g_c[destination_layer * organic_module.substrate_count + offset], fraction);
    for (0..organic_module.substrate_count * organic_module.structural_fraction_count) |offset| transferNumber(&source.colonized_structural_carbon_g_c[source_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset], &destination.colonized_structural_carbon_g_c[destination_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset], fraction);
}

fn validateOrganic(source: *const organic_module.State, source_layer: usize, destination: *const organic_module.State, destination_layer: usize, fraction: f64) !void {
    inline for (.{
        .{ "microbial", organic_module.microbial_substrate_count * organic_module.microbial_population_count * organic_module.kinetic_fraction_count },
        .{ "residue", organic_module.substrate_count * organic_module.residue_fraction_count },
        .{ "dissolved", organic_module.substrate_count },
        .{ "adsorbed", organic_module.substrate_count },
        .{ "structural", organic_module.substrate_count * organic_module.structural_fraction_count },
    }) |descriptor| {
        const name = descriptor[0];
        const stride: usize = descriptor[1];
        const source_values = @field(source, name);
        const destination_values = @field(destination, name);
        for (0..stride) |offset| try validateStructPair(organic_module.ElementPool, source_values[source_layer * stride + offset], destination_values[destination_layer * stride + offset], fraction);
    }
    inline for (.{ "dissolved_acetate_carbon_g_c", "adsorbed_acetate_carbon_g_c" }) |name| {
        const source_values = @field(source, name);
        const destination_values = @field(destination, name);
        for (0..organic_module.substrate_count) |offset| try validateNumberPair(source_values[source_layer * organic_module.substrate_count + offset], destination_values[destination_layer * organic_module.substrate_count + offset], fraction);
    }
    for (0..organic_module.substrate_count * organic_module.structural_fraction_count) |offset| try validateNumberPair(
        source.colonized_structural_carbon_g_c[source_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset],
        destination.colonized_structural_carbon_g_c[destination_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset],
        fraction,
    );
}

fn transferOrganic(source: *organic_module.State, source_layer: usize, destination: *organic_module.State, destination_layer: usize, fraction: f64) void {
    inline for (.{
        .{ "microbial", organic_module.microbial_substrate_count * organic_module.microbial_population_count * organic_module.kinetic_fraction_count },
        .{ "residue", organic_module.substrate_count * organic_module.residue_fraction_count },
        .{ "dissolved", organic_module.substrate_count },
        .{ "adsorbed", organic_module.substrate_count },
        .{ "structural", organic_module.substrate_count * organic_module.structural_fraction_count },
    }) |descriptor| {
        const name = descriptor[0];
        const stride: usize = descriptor[1];
        const source_values = @field(source, name);
        const destination_values = @field(destination, name);
        for (0..stride) |offset| transferStruct(organic_module.ElementPool, &source_values[source_layer * stride + offset], &destination_values[destination_layer * stride + offset], fraction);
    }
    inline for (.{ "dissolved_acetate_carbon_g_c", "adsorbed_acetate_carbon_g_c" }) |name| {
        const source_values = @field(source, name);
        const destination_values = @field(destination, name);
        for (0..organic_module.substrate_count) |offset| transferNumber(&source_values[source_layer * organic_module.substrate_count + offset], &destination_values[destination_layer * organic_module.substrate_count + offset], fraction);
    }
    for (0..organic_module.substrate_count * organic_module.structural_fraction_count) |offset| transferNumber(
        &source.colonized_structural_carbon_g_c[source_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset],
        &destination.colonized_structural_carbon_g_c[destination_layer * organic_module.substrate_count * organic_module.structural_fraction_count + offset],
        fraction,
    );
}

fn validateGas(source: *const gas_module.State, source_cell: usize, destination: *const gas_module.State, destination_cell: usize, fraction: f64) !void {
    try validateNumberPair(source.water_vapor_mol[source_cell], destination.water_vapor_mol[destination_cell], fraction);
    for (0..gas_module.species_count) |species| {
        const source_index = source_cell * gas_module.species_count + species;
        const destination_index = destination_cell * gas_module.species_count + species;
        try validateNumberPair(source.gaseous_mass_g[source_index], destination.gaseous_mass_g[destination_index], fraction);
        try validateNumberPair(source.dissolved_mass_g[source_index], destination.dissolved_mass_g[destination_index], fraction);
    }
}

fn transferGas(source: *gas_module.State, source_cell: usize, destination: *gas_module.State, destination_cell: usize, fraction: f64) void {
    transferNumber(&source.water_vapor_mol[source_cell], &destination.water_vapor_mol[destination_cell], fraction);
    for (0..gas_module.species_count) |species| {
        const source_index = source_cell * gas_module.species_count + species;
        const destination_index = destination_cell * gas_module.species_count + species;
        transferNumber(&source.gaseous_mass_g[source_index], &destination.gaseous_mass_g[destination_index], fraction);
        transferNumber(&source.dissolved_mass_g[source_index], &destination.dissolved_mass_g[destination_index], fraction);
    }
}

fn validateNitrogenFertilizer(source: surface_fertilizer_module.Inventory, destination: @import("../soil/nutrients/fertilizer_dissolution.zig").FertilizerState, fraction: f64) !void {
    inline for (.{
        .{ "ammonium_mol_n", "broadcast_ammonium_mol_n" },
        .{ "ammonia_mol_n", "broadcast_ammonia_mol_n" },
        .{ "urea_mol_n", "broadcast_urea_mol_n" },
        .{ "nitrate_mol_n", "broadcast_nitrate_mol_n" },
    }) |names| try validateNumberPair(@field(source, names[0]), @field(destination, names[1]), fraction);
}

fn transferNitrogenFertilizer(source_state: *surface_fertilizer_module.State, cell: usize, destination_state: *soil_fertilizer_module.State, destination: usize, fraction: f64) void {
    inline for (.{
        .{ "ammonium_mol_n", "broadcast_ammonium_mol_n" },
        .{ "ammonia_mol_n", "broadcast_ammonia_mol_n" },
        .{ "urea_mol_n", "broadcast_urea_mol_n" },
        .{ "nitrate_mol_n", "broadcast_nitrate_mol_n" },
    }) |names| transferNumber(&@field(source_state.cells[cell], names[0]), &@field(destination_state.soil[destination], names[1]), fraction);
}

fn transferMineralFertilizer(state: *mineral_fertilizer_module.State, cell: usize, destination: usize, fraction: f64) void {
    transferStruct(mineral_fertilizer_module.Inventory, &state.surface[cell], &state.soil[destination], fraction);
}

fn validateStructPair(comptime T: type, source: T, destination: T, fraction: f64) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| try validateNumberPair(@field(source, field.name), @field(destination, field.name), fraction);
}

fn validateNumberPair(source: f64, destination: f64, fraction: f64) !void {
    if (!std.math.isFinite(source) or source < 0 or !std.math.isFinite(destination) or destination < 0 or !std.math.isFinite(destination + fraction * source)) return error.InvalidSurfacePondInventory;
}

fn transferStruct(comptime T: type, source: *T, destination: *T, fraction: f64) void {
    inline for (@typeInfo(T).@"struct".fields) |field| transferNumber(&@field(source.*, field.name), &@field(destination.*, field.name), fraction);
}

fn transferNumber(source: *f64, destination: *f64, fraction: f64) void {
    const moved = fraction * source.*;
    destination.* += moved;
    source.* -= moved;
}

test "surface pond inventories transfer atomically into selected soil layer" {
    var surface_organic = try organic_module.State.init(std.testing.allocator, 1);
    defer surface_organic.deinit();
    var soil_organic = try organic_module.State.init(std.testing.allocator, 2);
    defer soil_organic.deinit();
    var surface_gas = try gas_module.State.init(std.testing.allocator, 1);
    defer surface_gas.deinit();
    var soil_gas = try gas_module.State.init(std.testing.allocator, 2);
    defer soil_gas.deinit();
    var surface_n = try surface_fertilizer_module.State.init(std.testing.allocator, 1);
    defer surface_n.deinit();
    var soil_n = try soil_fertilizer_module.State.init(std.testing.allocator, 1, 2);
    defer soil_n.deinit();
    var mineral = try mineral_fertilizer_module.State.init(std.testing.allocator, 1, 2);
    defer mineral.deinit();
    surface_organic.microbial[0].carbon_g_c = 8;
    surface_gas.gaseous_mass_g[0] = 6;
    surface_n.cells[0].ammonium_mol_n = 4;
    mineral.surface[0].gypsum_mol = 2;
    const owners: Owners = .{ .surface_organic = &surface_organic, .soil_organic = &soil_organic, .surface_gas = &surface_gas, .soil_gas = &soil_gas, .surface_nitrogen_fertilizer = &surface_n, .soil_nitrogen_fertilizer = &soil_n, .mineral_fertilizer = &mineral };
    try transferSurfaceFractionToSoil(owners, .{ .cell = 0, .destination_soil_layer = 1, .fraction = 0.25 });
    try std.testing.expectEqual(@as(f64, 6), surface_organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), soil_organic.microbial[organic_module.microbial_substrate_count * organic_module.microbial_population_count * organic_module.kinetic_fraction_count].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 1.5), soil_gas.gaseous_mass_g[gas_module.species_count]);
    try std.testing.expectEqual(@as(f64, 1), soil_n.soil[1].broadcast_ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 0.5), mineral.soil[1].gypsum_mol);
}

test "late invalid particulate inventory leaves every surface and soil owner unchanged" {
    var surface_organic = try organic_module.State.init(std.testing.allocator, 1);
    defer surface_organic.deinit();
    var soil_organic = try organic_module.State.init(std.testing.allocator, 1);
    defer soil_organic.deinit();
    var surface_gas = try gas_module.State.init(std.testing.allocator, 1);
    defer surface_gas.deinit();
    var soil_gas = try gas_module.State.init(std.testing.allocator, 1);
    defer soil_gas.deinit();
    var surface_n = try surface_fertilizer_module.State.init(std.testing.allocator, 1);
    defer surface_n.deinit();
    var soil_n = try soil_fertilizer_module.State.init(std.testing.allocator, 1, 1);
    defer soil_n.deinit();
    var mineral = try mineral_fertilizer_module.State.init(std.testing.allocator, 1, 1);
    defer mineral.deinit();
    surface_organic.microbial[0].carbon_g_c = 8;
    surface_gas.gaseous_mass_g[0] = 6;
    mineral.surface[0].potassium_ground_silicate_mol = std.math.nan(f64);
    const owners: Owners = .{ .surface_organic = &surface_organic, .soil_organic = &soil_organic, .surface_gas = &surface_gas, .soil_gas = &soil_gas, .surface_nitrogen_fertilizer = &surface_n, .soil_nitrogen_fertilizer = &soil_n, .mineral_fertilizer = &mineral };
    try std.testing.expectError(error.InvalidSurfacePondInventory, transferParticulateFractionToSoil(owners, .{ .cell = 0, .destination_soil_layer = 0, .fraction = 0.5 }));
    try std.testing.expectEqual(@as(f64, 8), surface_organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), soil_organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 6), surface_gas.gaseous_mass_g[0]);
}
