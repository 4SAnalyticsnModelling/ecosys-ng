const std = @import("std");
const retention = @import("soil_water_retention.zig");

pub const Options = struct {
    class_count: usize,
    pore_interaction_exponent: f64 = 1.33,
    air_entry_fraction_of_vertical_saturated_conductivity: f64,
};

pub const Table = struct {
    allocator: std.mem.Allocator,
    class_count: usize,
    /// Direction-major x/y/z, then water class from wettest to driest.
    conductivity_m2_per_h_mpa: []f64,
    air_entry_water_potential_mpa: f64,
    air_entry_water_fraction: f64,

    pub fn deinit(self: *Table) void {
        self.allocator.free(self.conductivity_m2_per_h_mpa);
        self.* = undefined;
    }

    pub fn at(self: *const Table, axis: usize, water_fraction: f64, porosity_fraction: f64) !f64 {
        if (axis >= 3) return error.InvalidHydraulicConductivityAxis;
        return self.conductivity_m2_per_h_mpa[axis * self.class_count + try classIndex(self.class_count, water_fraction, porosity_fraction)];
    }
};

/// HOUR1's Mualem-style HCND table, including its K-0.01*JK class offset.
/// The number of classes is runtime input rather than a PARAMETER constant.
pub fn build(allocator: std.mem.Allocator, curve: retention.ResolvedCurve, lateral_saturated_conductivity_m2_per_h_mpa: f64, vertical_saturated_conductivity_m2_per_h_mpa: f64, options: Options) !Table {
    if (options.class_count == 0 or !std.math.isFinite(lateral_saturated_conductivity_m2_per_h_mpa) or lateral_saturated_conductivity_m2_per_h_mpa < 0 or !std.math.isFinite(vertical_saturated_conductivity_m2_per_h_mpa) or vertical_saturated_conductivity_m2_per_h_mpa < 0 or !std.math.isFinite(options.pore_interaction_exponent) or options.pore_interaction_exponent <= 0 or !std.math.isFinite(options.air_entry_fraction_of_vertical_saturated_conductivity) or options.air_entry_fraction_of_vertical_saturated_conductivity < 0 or options.air_entry_fraction_of_vertical_saturated_conductivity > 1) return error.InvalidHydraulicConductivityOptions;
    const class_count_f64: f64 = @floatFromInt(options.class_count);
    const water_fraction = try allocator.alloc(f64, options.class_count);
    defer allocator.free(water_fraction);
    const potential_mpa = try allocator.alloc(f64, options.class_count);
    defer allocator.free(potential_mpa);
    var denominator: f64 = 0;
    for (0..options.class_count) |class| {
        const fortran_k: f64 = @floatFromInt(class + 1);
        const offset = @max(0.0, fortran_k - 0.01 * class_count_f64);
        water_fraction[class] = curve.porosity_fraction - offset / class_count_f64 * curve.porosity_fraction;
        potential_mpa[class] = try curve.waterPotentialMpa(@max(water_fraction[class], std.math.floatMin(f64)));
        denominator += (2.0 * fortran_k - 1.0) / (potential_mpa[class] * potential_mpa[class]);
    }
    if (!std.math.isFinite(denominator) or denominator <= 0) return error.InvalidHydraulicConductivityIntegral;
    const values = try allocator.alloc(f64, try std.math.mul(usize, 3, options.class_count));
    errdefer allocator.free(values);
    var air_entry_potential = curve.curve.saturation_water_potential_mpa;
    var air_entry_fraction = curve.porosity_fraction;
    var previous_vertical = vertical_saturated_conductivity_m2_per_h_mpa;
    for (0..options.class_count) |class| {
        const fortran_k: f64 = @floatFromInt(class + 1);
        const class_offset = @max(0.0, fortran_k - 0.01 * class_count_f64);
        const pore_factor = std.math.pow(f64, (class_count_f64 - class_offset) / class_count_f64, options.pore_interaction_exponent);
        var numerator: f64 = 0;
        for (class..options.class_count) |integral_class| {
            const fortran_m: f64 = @floatFromInt(integral_class + 1);
            const integral_offset = @max(0.0, fortran_m - 0.01 * class_count_f64);
            numerator += (2.0 * integral_offset + 1.0 - 2.0 * class_offset) / (potential_mpa[integral_class] * potential_mpa[integral_class]);
        }
        const relative = pore_factor * numerator / denominator;
        values[class] = lateral_saturated_conductivity_m2_per_h_mpa * relative;
        values[options.class_count + class] = lateral_saturated_conductivity_m2_per_h_mpa * relative;
        const vertical = vertical_saturated_conductivity_m2_per_h_mpa * relative;
        values[2 * options.class_count + class] = vertical;
        const threshold = options.air_entry_fraction_of_vertical_saturated_conductivity * vertical_saturated_conductivity_m2_per_h_mpa;
        if (class > 0 and vertical < threshold and previous_vertical >= threshold) {
            air_entry_potential = potential_mpa[class - 1];
            air_entry_fraction = water_fraction[class - 1];
        }
        previous_vertical = vertical;
    }
    for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteHydraulicConductivity;
    return .{ .allocator = allocator, .class_count = options.class_count, .conductivity_m2_per_h_mpa = values, .air_entry_water_potential_mpa = air_entry_potential, .air_entry_water_fraction = air_entry_fraction };
}

/// WATSUB K=MAX(1,MIN(JK,INT(JK*(POROS-THETA)/POROS)+1)),
/// returned as a zero-based index.
pub fn classIndex(class_count: usize, water_fraction: f64, porosity_fraction: f64) !usize {
    if (class_count == 0 or !std.math.isFinite(water_fraction) or water_fraction < 0 or !std.math.isFinite(porosity_fraction) or porosity_fraction <= 0) return error.InvalidHydraulicConductivityClassInput;
    const dry_fraction = std.math.clamp((porosity_fraction - water_fraction) / porosity_fraction, 0.0, 1.0);
    const raw: usize = @intFromFloat(@floor(@as(f64, @floatFromInt(class_count)) * dry_fraction));
    return @min(class_count - 1, raw);
}

fn testCurve() retention.ResolvedCurve {
    return .{ .porosity_fraction = 0.5, .curve = .{ .field_capacity_fraction = 0.3, .wilting_point_fraction = 0.1, .saturation_water_potential_mpa = -0.0005, .field_capacity_water_potential_mpa = -0.01, .wilting_point_water_potential_mpa = -1.5, .minimum_water_potential_mpa = -1.5e12, .saturation_to_field_shape = 0.5, .below_wilting_shape = 0.5 } };
}

test "runtime HCND table is saturated at wet class and declines when dry" {
    var table = try build(std.testing.allocator, testCurve(), 2, 1, .{ .class_count = 100, .air_entry_fraction_of_vertical_saturated_conductivity = 0.1 });
    defer table.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 2), try table.at(0, 0.5, 0.5), 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), try table.at(2, 0.5, 0.5), 1e-12);
    try std.testing.expect(try table.at(2, 0.1, 0.5) < try table.at(2, 0.4, 0.5));
    try std.testing.expect(table.air_entry_water_fraction > 0 and table.air_entry_water_fraction <= 0.5);
}

test "water class selection reproduces one-based Fortran clamp" {
    try std.testing.expectEqual(@as(usize, 0), try classIndex(100, 0.5, 0.5));
    try std.testing.expectEqual(@as(usize, 50), try classIndex(100, 0.25, 0.5));
    try std.testing.expectEqual(@as(usize, 99), try classIndex(100, 0, 0.5));
}
