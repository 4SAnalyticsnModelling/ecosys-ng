const std = @import("std");

/// Hourly extensive energy ledgers corresponding to OUTSH TRNS/TLES/TSHS/TGHS
/// (ground surface) and TRN/TLE/TSH/TGH (whole ecosystem). All storage is
/// runtime-sized and reset/refreshed once after the converged hourly kernels.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    ground_surface_net_radiation_mj: []f64,
    ground_surface_latent_heat_mj: []f64,
    ground_surface_sensible_heat_mj: []f64,
    ground_surface_storage_heat_mj: []f64,
    ecosystem_net_radiation_mj: []f64,
    ecosystem_latent_heat_mj: []f64,
    ecosystem_sensible_heat_mj: []f64,
    ecosystem_storage_heat_mj: []f64,
    canopy_water_energy_mj: []f64,
    canopy_water_energy_change_mj_per_h: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.EmptyEcosystemEnergyLedger;
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, cell_count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| {
            if (!std.math.isFinite(value)) {
                std.log.err("non-finite ecosystem energy ledger: field={s} cell={d} value={e}", .{ field.name, index, value });
                return error.NonFiniteEcosystemEnergyLedger;
            }
        };
    }
};

pub const Inputs = struct {
    cell_area_m2: []const f64,
    ground_net_radiation_mj_per_m2: []const f64,
    ground_latent_heat_mj_per_m2: []const f64,
    ground_sensible_heat_mj_per_m2: []const f64,
    ground_storage_heat_mj_per_m2: []const f64,
    species_count: usize,
    canopy_net_radiation_mj: []const f64,
    canopy_latent_heat_mj: []const f64,
    canopy_sensible_heat_mj: []const f64,
    canopy_storage_heat_mj: []const f64,
    canopy_convective_water_heat_mj: []const f64,
    standing_dead_net_radiation_mj: []const f64,
    standing_dead_latent_heat_mj: []const f64,
    standing_dead_sensible_heat_mj: []const f64,
    standing_dead_storage_heat_mj: []const f64,
    standing_dead_convective_water_heat_mj: []const f64,
};

pub fn refresh(state: *State, inputs: Inputs) !void {
    if (inputs.species_count == 0) return error.InvalidEcosystemEnergySpeciesCount;
    inline for (.{ inputs.cell_area_m2, inputs.ground_net_radiation_mj_per_m2, inputs.ground_latent_heat_mj_per_m2, inputs.ground_sensible_heat_mj_per_m2, inputs.ground_storage_heat_mj_per_m2 }) |values| if (values.len != state.cell_count) return error.EcosystemEnergyCellDimensionMismatch;
    const plants = try std.math.mul(usize, state.cell_count, inputs.species_count);
    inline for (.{ inputs.canopy_net_radiation_mj, inputs.canopy_latent_heat_mj, inputs.canopy_sensible_heat_mj, inputs.canopy_storage_heat_mj, inputs.canopy_convective_water_heat_mj, inputs.standing_dead_net_radiation_mj, inputs.standing_dead_latent_heat_mj, inputs.standing_dead_sensible_heat_mj, inputs.standing_dead_storage_heat_mj, inputs.standing_dead_convective_water_heat_mj }) |values| if (values.len != plants) return error.EcosystemEnergyPlantDimensionMismatch;

    for (0..state.cell_count) |cell| {
        const area = inputs.cell_area_m2[cell];
        if (!std.math.isFinite(area) or area <= 0) return error.InvalidEcosystemEnergyCellArea;
        const ground_net = inputs.ground_net_radiation_mj_per_m2[cell] * area;
        const ground_latent = inputs.ground_latent_heat_mj_per_m2[cell] * area;
        const ground_sensible = inputs.ground_sensible_heat_mj_per_m2[cell] * area;
        const ground_storage = inputs.ground_storage_heat_mj_per_m2[cell] * area;
        state.ground_surface_net_radiation_mj[cell] = ground_net;
        state.ground_surface_latent_heat_mj[cell] = ground_latent;
        state.ground_surface_sensible_heat_mj[cell] = ground_sensible;
        state.ground_surface_storage_heat_mj[cell] = ground_storage;
        var ecosystem_net = ground_net;
        var ecosystem_latent = ground_latent;
        var ecosystem_sensible = ground_sensible;
        var ecosystem_storage = ground_storage;
        for (cell * inputs.species_count..(cell + 1) * inputs.species_count) |plant| {
            inline for (.{ inputs.canopy_net_radiation_mj[plant], inputs.canopy_latent_heat_mj[plant], inputs.canopy_sensible_heat_mj[plant], inputs.canopy_storage_heat_mj[plant], inputs.canopy_convective_water_heat_mj[plant], inputs.standing_dead_net_radiation_mj[plant], inputs.standing_dead_latent_heat_mj[plant], inputs.standing_dead_sensible_heat_mj[plant], inputs.standing_dead_storage_heat_mj[plant], inputs.standing_dead_convective_water_heat_mj[plant] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteEcosystemEnergyInput;
            ecosystem_net += inputs.canopy_net_radiation_mj[plant] + inputs.standing_dead_net_radiation_mj[plant];
            ecosystem_latent += inputs.canopy_latent_heat_mj[plant] + inputs.standing_dead_latent_heat_mj[plant];
            ecosystem_sensible += inputs.canopy_sensible_heat_mj[plant] + inputs.standing_dead_sensible_heat_mj[plant];
            // Exact EXTRACT convention: TGH -= HFLXC - VFLXC.
            ecosystem_storage -= inputs.canopy_storage_heat_mj[plant] - inputs.canopy_convective_water_heat_mj[plant];
            ecosystem_storage -= inputs.standing_dead_storage_heat_mj[plant] - inputs.standing_dead_convective_water_heat_mj[plant];
        }
        state.ecosystem_net_radiation_mj[cell] = ecosystem_net;
        state.ecosystem_latent_heat_mj[cell] = ecosystem_latent;
        state.ecosystem_sensible_heat_mj[cell] = ecosystem_sensible;
        state.ecosystem_storage_heat_mj[cell] = ecosystem_storage;
    }
    try state.validateFinite();
}

test "energy ledger reproduces REDIST ground and EXTRACT ecosystem accumulation" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try refresh(&state, .{
        .cell_area_m2 = &.{10},
        .ground_net_radiation_mj_per_m2 = &.{1},
        .ground_latent_heat_mj_per_m2 = &.{2},
        .ground_sensible_heat_mj_per_m2 = &.{3},
        .ground_storage_heat_mj_per_m2 = &.{4},
        .species_count = 2,
        .canopy_net_radiation_mj = &.{ 1, 2 },
        .canopy_latent_heat_mj = &.{ 3, 4 },
        .canopy_sensible_heat_mj = &.{ 5, 6 },
        .canopy_storage_heat_mj = &.{ 7, 8 },
        .canopy_convective_water_heat_mj = &.{ 0.5, 1 },
        .standing_dead_net_radiation_mj = &.{ 0.1, 0.2 },
        .standing_dead_latent_heat_mj = &.{ 0.3, 0.4 },
        .standing_dead_sensible_heat_mj = &.{ 0.5, 0.6 },
        .standing_dead_storage_heat_mj = &.{ 0.7, 0.8 },
        .standing_dead_convective_water_heat_mj = &.{ 0.05, 0.1 },
    });
    try std.testing.expectEqual(@as(f64, 10), state.ground_surface_net_radiation_mj[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 13.3), state.ecosystem_net_radiation_mj[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 27.7), state.ecosystem_latent_heat_mj[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 42.1), state.ecosystem_sensible_heat_mj[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 25.15), state.ecosystem_storage_heat_mj[0], 1e-14);
}
