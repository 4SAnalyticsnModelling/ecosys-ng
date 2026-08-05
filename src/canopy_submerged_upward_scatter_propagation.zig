const std = @import("std");

pub const ScatterFlux = struct {
    forward_shortwave_megajoules_m2_h: f64,
    forward_par_umol_m2_s: f64,
    backscattered_shortwave_megajoules_m2_h: f64,
    backscattered_par_umol_m2_s: f64,
};

/// HOUR1 lines 1735--1738. Copies the boundary below through a submerged
/// canopy layer in exact RAFSL, RAFPL, RABSL, RABPL order.
pub fn propagate(below: ScatterFlux) !ScatterFlux {
    inline for (@typeInfo(ScatterFlux).@"struct".fields) |field| {
        const value = @field(below, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubmergedUpwardScatterFlux;
    }
    const forward_shortwave_megajoules_m2_h = below.forward_shortwave_megajoules_m2_h;
    const forward_par_umol_m2_s = below.forward_par_umol_m2_s;
    const backscattered_shortwave_megajoules_m2_h =
        below.backscattered_shortwave_megajoules_m2_h;
    const backscattered_par_umol_m2_s = below.backscattered_par_umol_m2_s;
    return .{
        .forward_shortwave_megajoules_m2_h = forward_shortwave_megajoules_m2_h,
        .forward_par_umol_m2_s = forward_par_umol_m2_s,
        .backscattered_shortwave_megajoules_m2_h = backscattered_shortwave_megajoules_m2_h,
        .backscattered_par_umol_m2_s = backscattered_par_umol_m2_s,
    };
}

test "submerged layer copies all upward scatter fluxes" {
    const result = try propagate(.{
        .forward_shortwave_megajoules_m2_h = 1,
        .forward_par_umol_m2_s = 2,
        .backscattered_shortwave_megajoules_m2_h = 3,
        .backscattered_par_umol_m2_s = 4,
    });
    try std.testing.expectEqual(@as(f64, 1), result.forward_shortwave_megajoules_m2_h);
    try std.testing.expectEqual(@as(f64, 2), result.forward_par_umol_m2_s);
    try std.testing.expectEqual(
        @as(f64, 3),
        result.backscattered_shortwave_megajoules_m2_h,
    );
    try std.testing.expectEqual(@as(f64, 4), result.backscattered_par_umol_m2_s);
}

test "nonfinite scatter flux fails explicitly" {
    try std.testing.expectError(
        error.InvalidSubmergedUpwardScatterFlux,
        propagate(.{
            .forward_shortwave_megajoules_m2_h = std.math.nan(f64),
            .forward_par_umol_m2_s = 0,
            .backscattered_shortwave_megajoules_m2_h = 0,
            .backscattered_par_umol_m2_s = 0,
        }),
    );
}
