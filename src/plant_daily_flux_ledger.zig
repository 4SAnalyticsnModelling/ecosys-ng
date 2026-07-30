const std = @import("std");
const builtin = @import("builtin");

/// Runtime-sized equivalent of the per-plant daily fields carried in BLK14
/// and reset by DAY after the preceding day's OUTPD record has been written.
pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,

    net_carbon_change_g: []f64,
    gross_primary_productivity_g: []f64,
    root_soil_carbon_exchange_g: []f64,
    carbon_sink_g: []f64,
    initial_carbon_sink_g: []f64,
    signed_total_respiration_carbon_g: []f64,
    signed_aboveground_respiration_carbon_g: []f64,
    transpiration_source_m3: []f64,
    ammonia_exchange_g_n: []f64,
    root_soil_nitrogen_exchange_g: []f64,
    symbiotic_nitrogen_fixation_g: []f64,
    nitrogen_sink_g: []f64,
    initial_nitrogen_sink_g: []f64,
    root_soil_phosphorus_exchange_g: []f64,
    phosphorus_sink_g: []f64,
    initial_phosphorus_sink_g: []f64,
    carbon_oxidation_g: []f64,
    nitrogen_oxidation_g: []f64,
    phosphorus_oxidation_g: []f64,
    harvested_carbon_g: []f64,
    harvested_nitrogen_g: []f64,
    harvested_phosphorus_g: []f64,

    cumulative_carbon_balance_g: []f64,
    cumulative_nitrogen_balance_g: []f64,
    cumulative_phosphorus_balance_g: []f64,
    cumulative_harvested_carbon_g: []f64,
    cumulative_harvested_nitrogen_g: []f64,
    cumulative_harvested_phosphorus_g: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize) !State {
        if (plant_count == 0) return error.InvalidPlantDailyFluxLedgerDimensions;
        var result: State = undefined;
        result.allocator = allocator;
        result.plant_count = plant_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64 and allocated > 0) {
                allocated -= 1;
                allocator.free(@field(result, field.name));
            }
        };
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(result, field.name) = try allocator.alloc(f64, plant_count);
                @memset(@field(result, field.name), 0);
                allocated += 1;
            }
        }
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub const HourlyExchange = struct {
        net_canopy_carbon_g: f64,
        gross_primary_productivity_g: f64,
        signed_total_respiration_carbon_g: f64,
        signed_aboveground_respiration_carbon_g: f64,
        canopy_and_standing_dead_water_source_m3: f64,
        canopy_ammonia_exchange_g_n: f64,
    };

    /// EXTRACT accumulates these three hourly quantities directly into
    /// CARBN, CTRAN, and TNH3C. Their source signs are retained here.
    pub fn accumulateHourlyExchange(self: *State, plant: usize, exchange: HourlyExchange) !void {
        if (plant >= self.plant_count) return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (std.meta.fields(HourlyExchange)) |field| {
            if (!std.math.isFinite(@field(exchange, field.name))) return error.NonFinitePlantDailyFluxLedgerInput;
        }
        self.net_carbon_change_g[plant] += exchange.net_canopy_carbon_g;
        self.gross_primary_productivity_g[plant] += exchange.gross_primary_productivity_g;
        self.signed_total_respiration_carbon_g[plant] += exchange.signed_total_respiration_carbon_g;
        self.signed_aboveground_respiration_carbon_g[plant] += exchange.signed_aboveground_respiration_carbon_g;
        self.transpiration_source_m3[plant] += exchange.canopy_and_standing_dead_water_source_m3;
        self.ammonia_exchange_g_n[plant] += exchange.canopy_ammonia_exchange_g_n;
        inline for (.{
            self.net_carbon_change_g[plant],
            self.gross_primary_productivity_g[plant],
            self.signed_total_respiration_carbon_g[plant],
            self.signed_aboveground_respiration_carbon_g[plant],
            self.transpiration_source_m3[plant],
            self.ammonia_exchange_g_n[plant],
        }) |value| if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedger;
    }

    pub fn accumulateHourlyRootSoilExchange(
        self: *State,
        plant: usize,
        carbon_g: f64,
        nitrogen_g: f64,
        phosphorus_g: f64,
        symbiotic_nitrogen_fixation_g: f64,
    ) !void {
        if (plant >= self.plant_count) return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (.{ carbon_g, nitrogen_g, phosphorus_g, symbiotic_nitrogen_fixation_g }) |value|
            if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedgerInput;
        self.root_soil_carbon_exchange_g[plant] += carbon_g;
        self.root_soil_nitrogen_exchange_g[plant] += nitrogen_g;
        self.root_soil_phosphorus_exchange_g[plant] += phosphorus_g;
        self.symbiotic_nitrogen_fixation_g[plant] += symbiotic_nitrogen_fixation_g;
        inline for (.{
            self.root_soil_carbon_exchange_g[plant],
            self.root_soil_nitrogen_exchange_g[plant],
            self.root_soil_phosphorus_exchange_g[plant],
            self.symbiotic_nitrogen_fixation_g[plant],
        }) |value| if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedger;
    }

    pub fn accumulateHarvest(self: *State, plant: usize, carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64) !void {
        if (plant >= self.plant_count) return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (.{ carbon_g, nitrogen_g, phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantDailyHarvest;
        self.harvested_carbon_g[plant] += carbon_g;
        self.harvested_nitrogen_g[plant] += nitrogen_g;
        self.harvested_phosphorus_g[plant] += phosphorus_g;
        inline for (.{ self.harvested_carbon_g[plant], self.harvested_nitrogen_g[plant], self.harvested_phosphorus_g[plant] }) |value|
            if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedger;
    }

    pub fn accumulateLitterSink(self: *State, plant: usize, carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64) !void {
        if (plant >= self.plant_count) return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (.{ carbon_g, nitrogen_g, phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantDailyLitterSink;
        self.carbon_sink_g[plant] += carbon_g;
        self.nitrogen_sink_g[plant] += nitrogen_g;
        self.phosphorus_sink_g[plant] += phosphorus_g;
        inline for (.{ self.carbon_sink_g[plant], self.nitrogen_sink_g[plant], self.phosphorus_sink_g[plant] }) |value|
            if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedger;
    }

    /// Above-ground GROSUB litter contributes to both total TCSNC/TZSNC/TPSNC
    /// and the surface-only TCSN0/TZSN0/TPSN0 ledgers.
    pub fn accumulateAbovegroundLitterSink(self: *State, plant: usize, carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64) !void {
        if (plant >= self.plant_count) return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (.{ carbon_g, nitrogen_g, phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantDailyLitterSink;
        const next_total_carbon = self.carbon_sink_g[plant] + carbon_g;
        const next_total_nitrogen = self.nitrogen_sink_g[plant] + nitrogen_g;
        const next_total_phosphorus = self.phosphorus_sink_g[plant] + phosphorus_g;
        const next_initial_carbon = self.initial_carbon_sink_g[plant] + carbon_g;
        const next_initial_nitrogen = self.initial_nitrogen_sink_g[plant] + nitrogen_g;
        const next_initial_phosphorus = self.initial_phosphorus_sink_g[plant] + phosphorus_g;
        inline for (.{ next_total_carbon, next_total_nitrogen, next_total_phosphorus, next_initial_carbon, next_initial_nitrogen, next_initial_phosphorus }) |value|
            if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedger;
        self.carbon_sink_g[plant] = next_total_carbon;
        self.nitrogen_sink_g[plant] = next_total_nitrogen;
        self.phosphorus_sink_g[plant] = next_total_phosphorus;
        self.initial_carbon_sink_g[plant] = next_initial_carbon;
        self.initial_nitrogen_sink_g[plant] = next_initial_nitrogen;
        self.initial_phosphorus_sink_g[plant] = next_initial_phosphorus;
    }

    /// Commits the complete EXTRACT plant litter publication atomically:
    /// aboveground material contributes to total and surface-only ledgers,
    /// while belowground material contributes only to total litter.
    pub fn accumulateLitterPublication(
        self: *State,
        plant: usize,
        aboveground_carbon_g_c: f64,
        aboveground_nitrogen_g_n: f64,
        aboveground_phosphorus_g_p: f64,
        belowground_carbon_g_c: f64,
        belowground_nitrogen_g_n: f64,
        belowground_phosphorus_g_p: f64,
    ) !void {
        if (plant >= self.plant_count)
            return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (.{
            aboveground_carbon_g_c,
            aboveground_nitrogen_g_n,
            aboveground_phosphorus_g_p,
            belowground_carbon_g_c,
            belowground_nitrogen_g_n,
            belowground_phosphorus_g_p,
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPlantDailyLitterSink;
        const next_total = .{
            self.carbon_sink_g[plant] + aboveground_carbon_g_c + belowground_carbon_g_c,
            self.nitrogen_sink_g[plant] + aboveground_nitrogen_g_n + belowground_nitrogen_g_n,
            self.phosphorus_sink_g[plant] + aboveground_phosphorus_g_p + belowground_phosphorus_g_p,
        };
        const next_surface = .{
            self.initial_carbon_sink_g[plant] + aboveground_carbon_g_c,
            self.initial_nitrogen_sink_g[plant] + aboveground_nitrogen_g_n,
            self.initial_phosphorus_sink_g[plant] + aboveground_phosphorus_g_p,
        };
        inline for (.{
            next_total[0],
            next_total[1],
            next_total[2],
            next_surface[0],
            next_surface[1],
            next_surface[2],
        }) |value|
            if (!std.math.isFinite(value))
                return error.NonFinitePlantDailyFluxLedger;
        self.carbon_sink_g[plant] = next_total[0];
        self.nitrogen_sink_g[plant] = next_total[1];
        self.phosphorus_sink_g[plant] = next_total[2];
        self.initial_carbon_sink_g[plant] = next_surface[0];
        self.initial_nitrogen_sink_g[plant] = next_surface[1];
        self.initial_phosphorus_sink_g[plant] = next_surface[2];
    }

    /// GROSUB/EXTRACT VCOXF/VNOXF/VPOXF retain signed plant-pool removals.
    /// Carbon includes the positive COR charcoal return supplied by the caller.
    pub fn accumulateOxidation(self: *State, plant: usize, carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64) !void {
        if (plant >= self.plant_count) return error.PlantDailyFluxLedgerIndexOutOfBounds;
        inline for (.{ carbon_g, nitrogen_g, phosphorus_g }) |value|
            if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedgerInput;
        const next_carbon = self.carbon_oxidation_g[plant] + carbon_g;
        const next_nitrogen = self.nitrogen_oxidation_g[plant] + nitrogen_g;
        const next_phosphorus = self.phosphorus_oxidation_g[plant] + phosphorus_g;
        inline for (.{ next_carbon, next_nitrogen, next_phosphorus }) |value|
            if (!std.math.isFinite(value)) return error.NonFinitePlantDailyFluxLedger;
        self.carbon_oxidation_g[plant] = next_carbon;
        self.nitrogen_oxidation_g[plant] = next_nitrogen;
        self.phosphorus_oxidation_g[plant] = next_phosphorus;
    }

    /// Applies the exact DAY carry equations, then clears only daily fields.
    /// Call this after OUTPD has consumed the completed day's values.
    pub fn closeDayAfterOutput(self: *State) !void {
        try self.validateFinite();
        for (0..self.plant_count) |plant| {
            self.cumulative_carbon_balance_g[plant] +=
                self.net_carbon_change_g[plant] +
                self.root_soil_carbon_exchange_g[plant] -
                self.carbon_sink_g[plant] +
                self.signed_total_respiration_carbon_g[plant] +
                self.carbon_oxidation_g[plant];
            self.cumulative_nitrogen_balance_g[plant] +=
                self.root_soil_nitrogen_exchange_g[plant] +
                self.ammonia_exchange_g_n[plant] -
                self.nitrogen_sink_g[plant] +
                self.nitrogen_oxidation_g[plant] +
                self.symbiotic_nitrogen_fixation_g[plant];
            self.cumulative_phosphorus_balance_g[plant] +=
                self.root_soil_phosphorus_exchange_g[plant] -
                self.phosphorus_sink_g[plant] +
                self.phosphorus_oxidation_g[plant];
            self.cumulative_harvested_carbon_g[plant] += self.harvested_carbon_g[plant];
            self.cumulative_harvested_nitrogen_g[plant] += self.harvested_nitrogen_g[plant];
            self.cumulative_harvested_phosphorus_g[plant] += self.harvested_phosphorus_g[plant];
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64 and !std.mem.startsWith(u8, field.name, "cumulative_"))
                @memset(@field(self, field.name), 0);
        }
        try self.validateFinite();
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) for (@field(self, field.name), 0..) |value, plant| {
                if (!std.math.isFinite(value)) {
                    if (!builtin.is_test) std.log.err("non-finite daily plant ledger: field={s} plant={d} value={e}", .{ field.name, plant, value });
                    return error.NonFinitePlantDailyFluxLedger;
                }
            };
        }
    }
};

test "DAY carry equations execute after OUTPD and reset daily fields" {
    var state = try State.init(std.testing.allocator, 7);
    defer state.deinit();
    const plant = 6;
    state.net_carbon_change_g[plant] = 10;
    state.root_soil_carbon_exchange_g[plant] = 2;
    state.carbon_sink_g[plant] = 3;
    state.signed_total_respiration_carbon_g[plant] = -4;
    state.carbon_oxidation_g[plant] = 1;
    state.root_soil_nitrogen_exchange_g[plant] = 5;
    state.ammonia_exchange_g_n[plant] = -1;
    state.nitrogen_sink_g[plant] = 2;
    state.nitrogen_oxidation_g[plant] = 0.5;
    state.symbiotic_nitrogen_fixation_g[plant] = 1.5;
    state.root_soil_phosphorus_exchange_g[plant] = 3;
    state.phosphorus_sink_g[plant] = 1;
    state.phosphorus_oxidation_g[plant] = 0.25;
    state.harvested_carbon_g[plant] = 8;
    state.harvested_nitrogen_g[plant] = 0.8;
    state.harvested_phosphorus_g[plant] = 0.08;

    // Values remain available to OUTPD until the close operation.
    try std.testing.expectEqual(@as(f64, 10), state.net_carbon_change_g[plant]);
    try state.closeDayAfterOutput();
    try std.testing.expectEqual(@as(f64, 6), state.cumulative_carbon_balance_g[plant]);
    try std.testing.expectEqual(@as(f64, 4), state.cumulative_nitrogen_balance_g[plant]);
    try std.testing.expectEqual(@as(f64, 2.25), state.cumulative_phosphorus_balance_g[plant]);
    try std.testing.expectEqual(@as(f64, 8), state.cumulative_harvested_carbon_g[plant]);
    try std.testing.expectEqual(@as(f64, 0), state.net_carbon_change_g[plant]);
    try std.testing.expectEqual(@as(f64, 0), state.transpiration_source_m3[plant]);
}

test "daily plant ledger fails immediately on non-finite accumulation" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.transpiration_source_m3[0] = std.math.nan(f64);
    try std.testing.expectError(error.NonFinitePlantDailyFluxLedger, state.closeDayAfterOutput());
}

test "EXTRACT hourly exchange retains source signs and arbitrary plant index" {
    var state = try State.init(std.testing.allocator, 12);
    defer state.deinit();
    try state.accumulateHourlyExchange(11, .{
        .net_canopy_carbon_g = 2.5,
        .gross_primary_productivity_g = 3,
        .signed_total_respiration_carbon_g = -0.4,
        .signed_aboveground_respiration_carbon_g = -0.3,
        .canopy_and_standing_dead_water_source_m3 = -0.003,
        .canopy_ammonia_exchange_g_n = -0.04,
    });
    try state.accumulateHourlyExchange(11, .{
        .net_canopy_carbon_g = -0.5,
        .gross_primary_productivity_g = 1,
        .signed_total_respiration_carbon_g = -0.2,
        .signed_aboveground_respiration_carbon_g = -0.1,
        .canopy_and_standing_dead_water_source_m3 = -0.002,
        .canopy_ammonia_exchange_g_n = 0.01,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.net_carbon_change_g[11], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.6), state.signed_total_respiration_carbon_g[11], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.4), state.signed_aboveground_respiration_carbon_g[11], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.005), state.transpiration_source_m3[11], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.03), state.ammonia_exchange_g_n[11], 1e-15);
}

test "GROSUB root-soil C N P exchange and fixation accumulate with source signs" {
    var state = try State.init(std.testing.allocator, 9);
    defer state.deinit();
    try state.accumulateHourlyRootSoilExchange(8, -0.2, 0.03, -0.004, 0.01);
    try state.accumulateHourlyRootSoilExchange(8, 0.1, 0.02, 0.001, 0.02);
    try std.testing.expectApproxEqAbs(@as(f64, -0.1), state.root_soil_carbon_exchange_g[8], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.root_soil_nitrogen_exchange_g[8], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.003), state.root_soil_phosphorus_exchange_g[8], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.03), state.symbiotic_nitrogen_fixation_g[8], 1e-15);
}

test "management exports accumulate in daily HVST elemental ledgers" {
    var state = try State.init(std.testing.allocator, 6);
    defer state.deinit();
    try state.accumulateHarvest(5, 3, 0.3, 0.03);
    try state.accumulateHarvest(5, 2, 0.2, 0.02);
    try std.testing.expectEqual(@as(f64, 5), state.harvested_carbon_g[5]);
    try std.testing.expectEqual(@as(f64, 0.5), state.harvested_nitrogen_g[5]);
    try std.testing.expectEqual(@as(f64, 0.05), state.harvested_phosphorus_g[5]);
}

test "senesced litter accumulates TCSNC TZSNC TPSNC without recycled pools" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try state.accumulateLitterSink(1, 4, 0.4, 0.04);
    try state.accumulateLitterSink(1, 1, 0.1, 0.01);
    try std.testing.expectEqual(@as(f64, 5), state.carbon_sink_g[1]);
    try std.testing.expectEqual(@as(f64, 0.5), state.nitrogen_sink_g[1]);
    try std.testing.expectEqual(@as(f64, 0.05), state.phosphorus_sink_g[1]);
}

test "aboveground litter also accumulates source TCSN0 TZSN0 TPSN0" {
    var state = try State.init(std.testing.allocator, 8);
    defer state.deinit();
    try state.accumulateAbovegroundLitterSink(7, 2, 0.2, 0.02);
    try state.accumulateLitterSink(7, 3, 0.3, 0.03);
    try std.testing.expectEqual(@as(f64, 5), state.carbon_sink_g[7]);
    try std.testing.expectEqual(@as(f64, 2), state.initial_carbon_sink_g[7]);
    try std.testing.expectEqual(@as(f64, 0.2), state.initial_nitrogen_sink_g[7]);
    try std.testing.expectEqual(@as(f64, 0.02), state.initial_phosphorus_sink_g[7]);
}

test "complete EXTRACT litter publication is atomic across surface and root totals" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.carbon_sink_g[0] = 4;
    state.initial_carbon_sink_g[0] = 2;
    state.phosphorus_sink_g[0] = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFinitePlantDailyFluxLedger,
        state.accumulateLitterPublication(0, 3, 0.3, 1, 5, 0.5, std.math.floatMax(f64)),
    );
    try std.testing.expectEqual(@as(f64, 4), state.carbon_sink_g[0]);
    try std.testing.expectEqual(@as(f64, 2), state.initial_carbon_sink_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.initial_phosphorus_sink_g[0]);
}

test "disturbance oxidation retains source signs and carbon charcoal return" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    try state.accumulateOxidation(2, -5.5 + 1.25, -0.4, -0.06);
    try std.testing.expectEqual(@as(f64, -4.25), state.carbon_oxidation_g[2]);
    try std.testing.expectEqual(@as(f64, -0.4), state.nitrogen_oxidation_g[2]);
    try std.testing.expectEqual(@as(f64, -0.06), state.phosphorus_oxidation_g[2]);
}
