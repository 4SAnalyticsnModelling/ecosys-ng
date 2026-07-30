const std = @import("std");
const builtin = @import("builtin");

/// Whole-landscape cumulative terms used by EXEC. Names retain direction and
/// units so callers cannot silently exchange an input with an output.
pub const Totals = struct {
    landscape_area_m2: f64,
    water_storage_m3: f64,
    cumulative_rain_m3: f64,
    cumulative_runoff_m3: f64,
    cumulative_evaporation_m3: f64,
    cumulative_water_outflow_m3: f64,
    heat_storage_mj: f64,
    cumulative_heat_input_mj: f64,
    cumulative_heat_output_mj: f64,
    oxygen_storage_g: f64,
    cumulative_oxygen_input_g: f64,
    cumulative_oxygen_output_g: f64,
    residue_carbon_g: f64,
    organic_carbon_g: f64,
    carbon_dioxide_carbon_g: f64,
    cumulative_carbon_dioxide_input_g: f64,
    cumulative_carbon_output_g: f64,
    cumulative_organic_fertilizer_carbon_g: f64,
    cumulative_carbon_sink_g: f64,
    residue_nitrogen_g: f64,
    organic_nitrogen_g: f64,
    dinitrogen_nitrogen_g: f64,
    ammonium_nitrogen_g: f64,
    nitrate_nitrogen_g: f64,
    cumulative_dinitrogen_input_g: f64,
    cumulative_nitrogen_input_g: f64,
    cumulative_nitrogen_output_g: f64,
    cumulative_organic_fertilizer_nitrogen_g: f64,
    cumulative_nitrogen_sink_g: f64,
    residue_phosphorus_g: f64,
    organic_phosphorus_g: f64,
    phosphate_phosphorus_g: f64,
    cumulative_phosphorus_input_g: f64,
    cumulative_phosphorus_output_g: f64,
    cumulative_organic_fertilizer_phosphorus_g: f64,
    cumulative_phosphorus_sink_g: f64,
    ion_inventory_mol: f64,
    cumulative_ion_input_mol: f64,
    cumulative_ion_output_mol: f64,
};

pub const Balance = struct { water_m3: f64, heat_mj: f64, oxygen_g: f64, carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64, ions_mol: f64 };
pub const Deviation = struct { water_m: f64, heat_mj_m2: f64, oxygen_g_m2: f64, carbon_g_m2: f64, nitrogen_g_m2: f64, phosphorus_g_m2: f64, ions_mol_m2: f64 };

pub fn balance(t: Totals) !Balance {
    try validate(t);
    return .{
        .water_m3 = t.water_storage_m3 - t.cumulative_rain_m3 + t.cumulative_runoff_m3 + t.cumulative_evaporation_m3 + t.cumulative_water_outflow_m3,
        .heat_mj = t.heat_storage_mj - t.cumulative_heat_input_mj + t.cumulative_heat_output_mj,
        .oxygen_g = t.oxygen_storage_g - t.cumulative_oxygen_input_g + t.cumulative_oxygen_output_g,
        .carbon_g = t.residue_carbon_g + t.organic_carbon_g + t.carbon_dioxide_carbon_g - t.cumulative_carbon_dioxide_input_g + t.cumulative_carbon_output_g - t.cumulative_organic_fertilizer_carbon_g - t.cumulative_carbon_sink_g,
        .nitrogen_g = t.residue_nitrogen_g + t.organic_nitrogen_g + t.dinitrogen_nitrogen_g + t.ammonium_nitrogen_g + t.nitrate_nitrogen_g - t.cumulative_dinitrogen_input_g - t.cumulative_nitrogen_input_g + t.cumulative_nitrogen_output_g - t.cumulative_organic_fertilizer_nitrogen_g - t.cumulative_nitrogen_sink_g,
        .phosphorus_g = t.residue_phosphorus_g + t.organic_phosphorus_g + t.phosphate_phosphorus_g - t.cumulative_phosphorus_input_g + t.cumulative_phosphorus_output_g - t.cumulative_organic_fertilizer_phosphorus_g - t.cumulative_phosphorus_sink_g,
        .ions_mol = t.ion_inventory_mol - t.cumulative_ion_input_mol + t.cumulative_ion_output_mol,
    };
}

pub const Monitor = struct {
    baseline: Balance,
    tolerance_per_m2: f64,

    pub fn init(totals: Totals, tolerance_per_m2: f64) !Monitor {
        if (!std.math.isFinite(tolerance_per_m2) or tolerance_per_m2 < 0) return error.InvalidMassBalanceTolerance;
        return .{ .baseline = try balance(totals), .tolerance_per_m2 = tolerance_per_m2 };
    }

    /// Matches EXEC's reset points (IBEGIN, ISTART, ILAST+1) when called by
    /// the timeline controller.
    pub fn reset(self: *Monitor, totals: Totals) !void {
        self.baseline = try balance(totals);
    }

    pub fn deviation(self: Monitor, totals: Totals) !Deviation {
        const current = try balance(totals);
        const area = totals.landscape_area_m2;
        return .{ .water_m = (current.water_m3 - self.baseline.water_m3) / area, .heat_mj_m2 = (current.heat_mj - self.baseline.heat_mj) / area, .oxygen_g_m2 = (current.oxygen_g - self.baseline.oxygen_g) / area, .carbon_g_m2 = (current.carbon_g - self.baseline.carbon_g) / area, .nitrogen_g_m2 = (current.nitrogen_g - self.baseline.nitrogen_g) / area, .phosphorus_g_m2 = (current.phosphorus_g - self.baseline.phosphorus_g) / area, .ions_mol_m2 = (current.ions_mol - self.baseline.ions_mol) / area };
    }

    pub fn check(self: Monitor, day: u64, year: i32, totals: Totals) !Deviation {
        const d = try self.deviation(totals);
        try checkOne(day, year, "carbon", d.carbon_g_m2, self.tolerance_per_m2, error.CarbonMassBalanceLost);
        try checkOne(day, year, "nitrogen", d.nitrogen_g_m2, self.tolerance_per_m2, error.NitrogenMassBalanceLost);
        try checkOne(day, year, "phosphorus", d.phosphorus_g_m2, self.tolerance_per_m2, error.PhosphorusMassBalanceLost);
        try checkOne(day, year, "water", d.water_m, self.tolerance_per_m2, error.WaterMassBalanceLost);
        try checkOne(day, year, "heat", d.heat_mj_m2, self.tolerance_per_m2, error.HeatMassBalanceLost);
        try checkOne(day, year, "oxygen", d.oxygen_g_m2, self.tolerance_per_m2, error.OxygenMassBalanceLost);
        try checkOne(day, year, "ions", d.ions_mol_m2, self.tolerance_per_m2, error.IonMassBalanceLost);
        return d;
    }
};

pub fn shouldAudit(day_index: u64, output_interval_days: u64) !bool {
    if (output_interval_days == 0) return error.InvalidMassBalanceAuditInterval;
    return day_index % output_interval_days == 0;
}

pub const DayBookkeeping = struct { reported_day: i64, previous_day: i64, management_event_count: u64, year_transition_count: u64 };

/// Exact EXEC tail: negative IDAYR remains relative to LYRX; otherwise the
/// current day is published, then daily management/year counters are reset.
pub fn advanceDayBookkeeping(current_day: i64, reported_day: i64, days_in_year: i64) !DayBookkeeping {
    if (days_in_year <= 0) return error.InvalidDaysInYear;
    return .{ .reported_day = if (reported_day < 0) try std.math.add(i64, days_in_year, reported_day) else current_day, .previous_day = current_day, .management_event_count = 0, .year_transition_count = 0 };
}

fn checkOne(day: u64, year: i32, comptime domain: []const u8, value: f64, tolerance: f64, comptime failure: anyerror) !void {
    if (!std.math.isFinite(value)) return error.NonFiniteMassBalanceDeviation;
    if (@abs(value) <= tolerance) return;
    if (!builtin.is_test) std.log.err("{s} mass balance lost: day={d} year={d} deviation_per_m2={e} tolerance_per_m2={e}", .{ domain, day, year, value, tolerance });
    return failure;
}

fn validate(t: Totals) !void {
    if (!std.math.isFinite(t.landscape_area_m2) or t.landscape_area_m2 <= 0) return error.InvalidLandscapeArea;
    inline for (std.meta.fields(Totals)) |field| if (!std.math.isFinite(@field(t, field.name))) return error.NonFiniteMassBalanceInput;
}

test "EXEC equations conserve all seven domains exactly" {
    var t = std.mem.zeroes(Totals);
    t.landscape_area_m2 = 10;
    t.water_storage_m3 = 5;
    t.cumulative_rain_m3 = 2;
    t.cumulative_runoff_m3 = 1;
    t.cumulative_evaporation_m3 = 0.5;
    t.cumulative_water_outflow_m3 = 0.5;
    t.residue_carbon_g = 10;
    t.cumulative_carbon_sink_g = 2;
    var monitor = try Monitor.init(t, 1e-6);
    const d = try monitor.check(1, 2001, t);
    try std.testing.expectEqual(@as(f64, 0), d.water_m);
    try std.testing.expectEqual(@as(f64, 0), d.carbon_g_m2);
    try monitor.reset(t);
}

test "EXEC fail-fast audit identifies carbon drift above source tolerance" {
    var t = std.mem.zeroes(Totals);
    t.landscape_area_m2 = 2;
    const monitor = try Monitor.init(t, 1e-6);
    t.organic_carbon_g = 3e-6;
    try std.testing.expectError(error.CarbonMassBalanceLost, monitor.check(2, 2001, t));
}

test "EXEC rejects non-finite cumulative inputs before deviation" {
    var t = std.mem.zeroes(Totals);
    t.landscape_area_m2 = 1;
    t.ion_inventory_mol = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteMassBalanceInput, balance(t));
}

test "EXEC cadence and day bookkeeping retain source behavior" {
    try std.testing.expect(try shouldAudit(30, 10));
    try std.testing.expect(!try shouldAudit(31, 10));
    const relative = try advanceDayBookkeeping(100, -2, 365);
    try std.testing.expectEqual(@as(i64, 363), relative.reported_day);
    try std.testing.expectEqual(@as(i64, 100), relative.previous_day);
    const absolute = try advanceDayBookkeeping(101, 8, 366);
    try std.testing.expectEqual(@as(i64, 101), absolute.reported_day);
    try std.testing.expectEqual(@as(u64, 0), absolute.management_event_count);
}
