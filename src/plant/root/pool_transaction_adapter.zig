const std = @import("std");
const admission = @import("rooted_layer_admission.zig");
const nutrient = @import("plant_root_nutrient_uptake.zig");
const organic = @import("plant_root_exudation.zig");
const replay = @import("pool_transaction_replay.zig");
const salt = @import("plant_root_salt_exchange.zig");

pub const Selection = struct { nutrient: []const bool, organic: []const bool, salt: []const bool };

/// Caller-owned adapter storage. Organic slices in `tuples` point into
/// `organic_storage` and remain valid until the next build.
pub const Workspace = struct {
    tuples: []replay.StagedTuple,
    organic_storage: []replay.CarbonNitrogenPhosphorus,
    requested_charge_mol: []f64,

    pub fn init(allocator: std.mem.Allocator, admission_capacity: usize) !Workspace {
        if (admission_capacity == 0) return error.ZeroRootTransactionAdmissionCapacity;
        const tuples = try allocator.alloc(replay.StagedTuple, admission_capacity);
        errdefer allocator.free(tuples);
        const organic_storage = try allocator.alloc(replay.CarbonNitrogenPhosphorus, try std.math.mul(usize, admission_capacity, organic.substrate_count));
        errdefer allocator.free(organic_storage);
        const charge = try allocator.alloc(f64, admission_capacity);
        return .{ .tuples = tuples, .organic_storage = organic_storage, .requested_charge_mol = charge };
    }

    pub fn deinit(self: *Workspace, allocator: std.mem.Allocator) void {
        allocator.free(self.requested_charge_mol);
        allocator.free(self.organic_storage);
        allocator.free(self.tuples);
        self.* = undefined;
    }

    pub fn build(
        self: *Workspace,
        coordinates: []const admission.Coordinate,
        selection: Selection,
        nutrient_by_tuple: []const nutrient.TransactionMapping,
        organic_by_tuple: []const organic.TransactionMapping,
        salt_by_tuple: []const salt.TransactionMapping,
        nitrogen_molar_mass_g_per_mol: f64,
        phosphorus_molar_mass_g_per_mol: f64,
    ) ![]const replay.StagedTuple {
        const count = coordinates.len;
        if (count > self.tuples.len or selection.nutrient.len != count or selection.organic.len != count or selection.salt.len != count or nutrient_by_tuple.len != count or organic_by_tuple.len != count or salt_by_tuple.len != count) return error.RootTransactionAdapterDimensionMismatch;
        if (!std.math.isFinite(nitrogen_molar_mass_g_per_mol) or nitrogen_molar_mass_g_per_mol <= 0 or !std.math.isFinite(phosphorus_molar_mass_g_per_mol) or phosphorus_molar_mass_g_per_mol <= 0) return error.InvalidRootTransactionMolarMass;
        for (coordinates, 0..) |coordinate, tuple| {
            const organic_start = tuple * organic.substrate_count;
            for (organic_by_tuple[tuple], 0..) |value, substrate| self.organic_storage[organic_start + substrate] = .{ .carbon_g_c = value.carbon_g_c, .nitrogen_g_n = value.nitrogen_g_n, .phosphorus_g_p = value.phosphorus_g_p };
            const charge = try requestedCharge(nutrient_by_tuple[tuple], salt_by_tuple[tuple], selection.nutrient[tuple], selection.salt[tuple], nitrogen_molar_mass_g_per_mol, phosphorus_molar_mass_g_per_mol);
            self.requested_charge_mol[tuple] = charge;
            self.tuples[tuple] = .{
                .coordinate = coordinate,
                .selected = .{ .nutrient = selection.nutrient[tuple], .organic = selection.organic[tuple], .dynamic_salt = selection.salt[tuple] },
                .nutrient_uptake_g_element = nutrient_by_tuple[tuple].uptake_g_element,
                .nutrient_respiration_cost = nutrient_by_tuple[tuple].respiration_cost,
                .organic_exchange_by_substrate = self.organic_storage[organic_start..][0..organic.substrate_count],
                .salt_exchange_mol = salt_by_tuple[tuple],
                .hydrogen_charge_mol = charge,
            };
        }
        return self.tuples[0..count];
    }

    /// Applies the intentional GROSUB-024 aggregate available-H+ floor. The
    /// correction is assigned to the final selected tuple in each layer so
    /// tuple charges sum exactly to the production-applied layer charge.
    pub fn boundAggregateCharge(self: *Workspace, tuples: []replay.StagedTuple, soil_layer_count: usize, available_hydrogen_mol_by_layer: []const f64) !void {
        _ = self;
        if (available_hydrogen_mol_by_layer.len != soil_layer_count) return error.RootTransactionAdapterDimensionMismatch;
        for (0..soil_layer_count) |layer| {
            if (!std.math.isFinite(available_hydrogen_mol_by_layer[layer]) or available_hydrogen_mol_by_layer[layer] < 0) return error.InvalidRootTransactionHydrogenInventory;
            var requested: f64 = 0;
            var final_selected: ?usize = null;
            for (tuples, 0..) |tuple, index| if (tuple.coordinate.soil_layer == layer and (tuple.selected.nutrient or tuple.selected.dynamic_salt)) {
                requested += tuple.hydrogen_charge_mol;
                final_selected = index;
            };
            if (!std.math.isFinite(requested)) return error.NonFiniteRootHydrogenCharge;
            if (final_selected) |index| tuples[index].hydrogen_charge_mol += @max(-available_hydrogen_mol_by_layer[layer], requested) - requested;
        }
    }
};

pub const GridWorkspace = struct {
    allocator: std.mem.Allocator,
    per_cell: []Workspace,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, admission_capacity_per_cell: usize) !GridWorkspace {
        if (cell_count == 0) return error.ZeroRootTransactionGridCells;
        const cells = try allocator.alloc(Workspace, cell_count);
        errdefer allocator.free(cells);
        var initialized: usize = 0;
        errdefer for (cells[0..initialized]) |*cell| cell.deinit(allocator);
        for (cells) |*cell| {
            cell.* = try Workspace.init(allocator, admission_capacity_per_cell);
            initialized += 1;
        }
        return .{ .allocator = allocator, .per_cell = cells };
    }

    pub fn deinit(self: *GridWorkspace) void {
        for (self.per_cell) |*cell| cell.deinit(self.allocator);
        self.allocator.free(self.per_cell);
        self.* = undefined;
    }
};

fn requestedCharge(nutrients: nutrient.TransactionMapping, salts: salt.TransactionMapping, nutrient_selected: bool, salt_selected: bool, nitrogen_mass: f64, phosphorus_mass: f64) !f64 {
    var charge: f64 = 0;
    if (nutrient_selected) {
        charge += (nutrients.uptake_g_element[0] + nutrients.uptake_g_element[1] - nutrients.uptake_g_element[2] - nutrients.uptake_g_element[3]) / nitrogen_mass;
        charge -= (2 * (nutrients.uptake_g_element[4] + nutrients.uptake_g_element[5]) + nutrients.uptake_g_element[6] + nutrients.uptake_g_element[7]) / phosphorus_mass;
    }
    if (salt_selected) charge += 3 * (salts[0] + salts[1]) + 2 * (salts[2] + salts[3]) + salts[4] + salts[5] - 2 * salts[6] - salts[7];
    if (!std.math.isFinite(charge)) return error.NonFiniteRootHydrogenCharge;
    return charge;
}

/// Replay currently lacks these compatibility destinations, so production
/// cutover must remain disabled until publication is extended atomically.
pub const MissingPublication = enum {
    root_exudate_exchange_ledgers,
    root_salt_uptake_ledger,
    plant_available_organic_mirrors,
    plant_available_mineral_mirrors,
    aqueous_concentration_publication,
    positive_respiration_ledger_adapter,
};

pub const production_cutover_blockers = [_]MissingPublication{
    .root_exudate_exchange_ledgers,
    .root_salt_uptake_ledger,
    .plant_available_organic_mirrors,
    .plant_available_mineral_mirrors,
    .aqueous_concentration_publication,
    .positive_respiration_ledger_adapter,
};

test "mixed transaction adapter preserves unequal NI order and aggregate charge floor" {
    var workspace = try Workspace.init(std.testing.allocator, 3);
    defer workspace.deinit(std.testing.allocator);
    const coordinates = [_]admission.Coordinate{ .{ .plant = 0, .biological_domain = 0, .soil_layer = 0 }, .{ .plant = 0, .biological_domain = 0, .soil_layer = 1 }, .{ .plant = 1, .biological_domain = 0, .soil_layer = 0 } };
    var nutrients = [_]nutrient.TransactionMapping{.{ .uptake_g_element = .{0} ** 8, .respiration_cost = .{} }} ** 3;
    nutrients[0].uptake_g_element[0] = 14;
    var organics = [_]organic.TransactionMapping{.{organic.Result{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }} ** organic.substrate_count} ** 3;
    organics[2][0].carbon_g_c = -0.25;
    var salts = [_]salt.TransactionMapping{.{0} ** 8} ** 3;
    salts[2][7] = 2;
    const tuples = try workspace.build(&coordinates, .{ .nutrient = &.{ true, false, false }, .organic = &.{ false, false, true }, .salt = &.{ false, false, true } }, &nutrients, &organics, &salts, 14, 31);
    try workspace.boundAggregateCharge(@constCast(tuples), 2, &.{ 0.5, 0.5 });
    try std.testing.expectEqual(coordinates[2], tuples[2].coordinate);
    try std.testing.expectEqual(@as(f64, -0.5), tuples[0].hydrogen_charge_mol + tuples[2].hydrogen_charge_mol);
    try std.testing.expectEqual(@as(f64, -0.25), tuples[2].organic_exchange_by_substrate[0].carbon_g_c);
}

test "adapter reports every missing atomic production destination" {
    try std.testing.expectEqual(@as(usize, 6), production_cutover_blockers.len);
}
