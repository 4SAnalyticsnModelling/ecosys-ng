const std = @import("std");
const retention = @import("soil_water_retention.zig");

pub const FaceDirection = enum { horizontal, vertical };

/// Runtime inputs for one WATSUB micropore face. Positive flux is source to
/// destination. Volumes are m3, potentials MPa, conductivity m2 h-1 MPa-1,
/// face area m2, and path lengths m.
pub const MatrixFaceInputs = struct {
    direction: FaceDirection,
    source_water_m3: f64,
    destination_water_m3: f64,
    source_air_m3: f64,
    destination_air_m3: f64,
    source_micropore_volume_m3: f64,
    destination_micropore_volume_m3: f64,
    source_water_fraction: f64,
    destination_water_fraction: f64,
    source_air_entry_water_fraction: f64,
    destination_air_entry_water_fraction: f64,
    source_total_water_potential_mpa: f64,
    destination_total_water_potential_mpa: f64,
    source_hydraulic_conductivity_m2_per_h_mpa: f64,
    destination_hydraulic_conductivity_m2_per_h_mpa: f64,
    source_path_length_m: f64,
    destination_path_length_m: f64,
    face_area_m2: f64,
    wetting_depth_factor: f64,
    time_fraction: f64,
    source_incoming_vertical_water_m3: f64 = 0,
    source_excess_pore_volume_m3: f64 = 0,
    destination_excess_pore_volume_m3: f64 = 0,
};

pub const MatrixFaceFlux = struct {
    /// FLQX: conductive flux before water/air availability limits.
    unlimited_water_m3: f64,
    /// FLQZ: conductive flux plus air-entry wetting-front enhancement.
    wetting_enhanced_water_m3: f64,
    /// FLQL: storage-changing flux.
    limited_water_m3: f64,
    /// FLQ2: unenhanced, availability-limited flux consumed by transport.
    transport_water_m3: f64,
};

pub fn calculateMatrixFaceFlux(inputs: MatrixFaceInputs) !MatrixFaceFlux {
    try validateMatrixInputs(inputs);
    const denominator = inputs.source_hydraulic_conductivity_m2_per_h_mpa * inputs.destination_path_length_m +
        inputs.destination_hydraulic_conductivity_m2_per_h_mpa * inputs.source_path_length_m;
    const conductance_m_per_h_mpa = if (denominator > 0)
        2.0 * inputs.source_hydraulic_conductivity_m2_per_h_mpa * inputs.destination_hydraulic_conductivity_m2_per_h_mpa / denominator
    else
        0.0;
    const unlimited = conductance_m_per_h_mpa *
        (inputs.source_total_water_potential_mpa - inputs.destination_total_water_potential_mpa) *
        inputs.face_area_m2 * inputs.time_fraction;
    const depth_factor = @min(1.0, inputs.wetting_depth_factor);
    var enhanced = unlimited;
    if (unlimited >= 0 and inputs.source_water_fraction > inputs.source_air_entry_water_fraction) {
        const source_excess = (inputs.source_water_fraction - inputs.source_air_entry_water_fraction) * inputs.source_micropore_volume_m3;
        const destination_deficit = @max(0.0, (inputs.destination_air_entry_water_fraction - inputs.destination_water_fraction) * inputs.destination_micropore_volume_m3);
        enhanced += @min(source_excess, destination_deficit) * depth_factor;
    } else if (unlimited < 0 and inputs.destination_water_fraction > inputs.destination_air_entry_water_fraction) {
        const destination_excess = (inputs.destination_air_entry_water_fraction - inputs.destination_water_fraction) * inputs.destination_micropore_volume_m3;
        const source_deficit = @min(0.0, (inputs.source_water_fraction - inputs.source_air_entry_water_fraction) * inputs.source_micropore_volume_m3);
        enhanced += @max(destination_excess, source_deficit) * depth_factor;
    }

    const source_available = @max(0.0, inputs.source_water_m3 + if (inputs.direction == .vertical) inputs.source_incoming_vertical_water_m3 else 0.0);
    var limited: f64 = undefined;
    var transport: f64 = undefined;
    if (unlimited >= 0) {
        limited = @max(0.0, @min(enhanced, @min(source_available * inputs.time_fraction, inputs.destination_air_m3 * inputs.time_fraction)));
        transport = @max(0.0, @min(unlimited, @min(source_available * inputs.time_fraction, inputs.destination_air_m3 * inputs.time_fraction)));
    } else {
        const source_air = @max(0.0, inputs.source_air_m3 - if (inputs.direction == .vertical) inputs.source_incoming_vertical_water_m3 else 0.0);
        limited = @min(0.0, @max(enhanced, @max(-inputs.destination_water_m3 * inputs.time_fraction, -source_air * inputs.time_fraction)));
        transport = @min(0.0, @max(unlimited, @max(-inputs.destination_water_m3 * inputs.time_fraction, -source_air * inputs.time_fraction)));
    }
    if (inputs.direction == .vertical and inputs.destination_excess_pore_volume_m3 < 0) {
        const freezing_displacement = @min(0.0, @max(-inputs.destination_water_m3 * inputs.time_fraction, inputs.destination_excess_pore_volume_m3));
        limited += freezing_displacement;
        transport += freezing_displacement;
    }
    if (!std.math.isFinite(limited) or !std.math.isFinite(transport)) return error.NonFiniteSoilWaterFlux;
    return .{ .unlimited_water_m3 = unlimited, .wetting_enhanced_water_m3 = enhanced, .limited_water_m3 = limited, .transport_water_m3 = transport };
}

pub const MacroporeFaceInputs = struct {
    direction: FaceDirection,
    source_water_m3: f64,
    destination_water_m3: f64,
    source_air_m3: f64,
    destination_air_m3: f64,
    source_total_volume_m3: f64,
    destination_total_volume_m3: f64,
    source_gravitational_potential_mpa: f64,
    destination_gravitational_potential_mpa: f64,
    source_vertical_thickness_m: f64,
    destination_vertical_thickness_m: f64,
    hydraulic_conductance_m_per_h_mpa: f64,
    face_area_m2: f64,
    time_fraction: f64,
    incoming_vertical_water_m3: f64 = 0,
    destination_excess_pore_volume_m3: f64 = 0,
};

pub const MacroporeFaceFlux = struct { unlimited_water_m3: f64, limited_water_m3: f64 };

pub const DualDomainExchangeInputs = struct {
    matrix_parameters: retention.MualemVanGenuchtenParameters,
    matrix_pressure_head_m: f64,
    macropore_pressure_head_m: f64,
    matrix_bulk_volume_m3: f64,
    characteristic_matrix_length_m: f64,
    geometry_factor: f64,
    scaling_coefficient: f64,
    time_fraction: f64,
    matrix_water_m3: f64,
    matrix_air_m3: f64,
    macropore_water_m3: f64,
    macropore_air_m3: f64,
};

/// Gerke–van Genuchten first-order dual-permeability transfer. Positive
/// transfer is macropore to matrix. The interface conductivity is the
/// published arithmetic mean evaluated at both domain pressure heads.
pub fn calculateDualDomainExchange(inputs: DualDomainExchangeInputs) !f64 {
    inline for (@typeInfo(DualDomainExchangeInputs).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteDualDomainExchangeInput;
    }
    try inputs.matrix_parameters.validate();
    if (inputs.matrix_bulk_volume_m3 <= 0 or
        inputs.characteristic_matrix_length_m <= 0 or
        inputs.geometry_factor <= 0 or inputs.scaling_coefficient <= 0 or
        inputs.time_fraction <= 0 or inputs.time_fraction > 1 or
        inputs.matrix_water_m3 < 0 or inputs.matrix_air_m3 < 0 or
        inputs.macropore_water_m3 < 0 or inputs.macropore_air_m3 < 0)
    {
        return error.InvalidDualDomainExchangeInput;
    }
    const matrix_interface_conductivity_m_per_h =
        try inputs.matrix_parameters.hydraulicConductivityMPerH(
            inputs.matrix_pressure_head_m,
        );
    const macropore_interface_conductivity_m_per_h =
        try inputs.matrix_parameters.hydraulicConductivityMPerH(
            inputs.macropore_pressure_head_m,
        );
    const interface_conductivity_m_per_h =
        0.5 * (matrix_interface_conductivity_m_per_h +
            macropore_interface_conductivity_m_per_h);
    const transfer_coefficient_per_m_h =
        inputs.geometry_factor * inputs.scaling_coefficient *
        interface_conductivity_m_per_h /
        (inputs.characteristic_matrix_length_m *
            inputs.characteristic_matrix_length_m);
    const unlimited_m3 =
        transfer_coefficient_per_m_h *
        (inputs.macropore_pressure_head_m - inputs.matrix_pressure_head_m) *
        inputs.matrix_bulk_volume_m3 * inputs.time_fraction;
    const limited_m3 = if (unlimited_m3 >= 0)
        @min(unlimited_m3, @min(inputs.macropore_water_m3, inputs.matrix_air_m3))
    else
        @max(unlimited_m3, -@min(inputs.matrix_water_m3, inputs.macropore_air_m3));
    if (!std.math.isFinite(limited_m3))
        return error.NonFiniteDualDomainExchange;
    return limited_m3;
}

pub fn calculateMacroporeFaceFlux(inputs: MacroporeFaceInputs) !MacroporeFaceFlux {
    inline for (@typeInfo(MacroporeFaceInputs).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSoilWaterFluxInput;
    }
    if (inputs.source_water_m3 < 0 or inputs.destination_water_m3 < 0 or inputs.source_air_m3 < 0 or inputs.destination_air_m3 < 0 or inputs.source_total_volume_m3 <= 0 or inputs.destination_total_volume_m3 <= 0 or inputs.source_vertical_thickness_m <= 0 or inputs.destination_vertical_thickness_m <= 0 or inputs.hydraulic_conductance_m_per_h_mpa < 0 or inputs.face_area_m2 < 0 or inputs.time_fraction <= 0 or inputs.time_fraction > 1) return error.InvalidSoilWaterFluxInput;
    const source_potential = inputs.source_gravitational_potential_mpa + 0.0098 * inputs.source_vertical_thickness_m * (@min(1.0, @max(0.0, inputs.source_water_m3 / inputs.source_total_volume_m3)) - 0.5);
    const destination_potential = inputs.destination_gravitational_potential_mpa + 0.0098 * inputs.destination_vertical_thickness_m * (@min(1.0, @max(0.0, inputs.destination_water_m3 / inputs.destination_total_volume_m3)) - 0.5);
    const unlimited = inputs.hydraulic_conductance_m_per_h_mpa * (source_potential - destination_potential) * inputs.face_area_m2 * inputs.time_fraction;
    var limited: f64 = 0;
    if (inputs.direction == .horizontal) {
        if (source_potential > destination_potential) limited = @max(0.0, @min(unlimited, @min(inputs.source_water_m3, inputs.destination_air_m3) * inputs.time_fraction));
        if (source_potential < destination_potential) limited = @min(0.0, @max(unlimited, @max(-inputs.destination_water_m3, -inputs.source_air_m3) * inputs.time_fraction));
    } else {
        limited = @max(0.0, @min(unlimited, @min(inputs.source_water_m3 * inputs.time_fraction + inputs.incoming_vertical_water_m3, inputs.destination_air_m3 * inputs.time_fraction)));
        if (inputs.destination_excess_pore_volume_m3 < 0) limited += @min(0.0, @max(-inputs.destination_water_m3 * inputs.time_fraction, inputs.destination_excess_pore_volume_m3));
    }
    if (!std.math.isFinite(unlimited) or !std.math.isFinite(limited)) return error.NonFiniteSoilWaterFlux;
    return .{ .unlimited_water_m3 = unlimited, .limited_water_m3 = limited };
}

pub const VaporFaceInputs = struct {
    source_air_volume_m3: f64,
    destination_air_volume_m3: f64,
    source_vapor_volume_m3: f64,
    destination_vapor_volume_m3: f64,
    source_vapor_diffusivity_m2_per_h: f64,
    destination_vapor_diffusivity_m2_per_h: f64,
    source_air_fraction: f64,
    destination_air_fraction: f64,
    source_porosity_fraction: f64,
    destination_porosity_fraction: f64,
    tortuosity: f64,
    source_path_length_m: f64,
    destination_path_length_m: f64,
    face_area_m2: f64,
    time_fraction: f64,
};

pub const VaporFaceFlux = struct {
    unlimited_vapor_m3: f64,
    limited_vapor_m3: f64,
    conductance_m_per_h: f64,
};

/// Exact WATSUB air-filled-pore vapor diffusion term. Unlike liquid flux,
/// XNPHX is already contained in the unlimited flux and is not reapplied to
/// the donor-vapor availability bound in the source model.
pub fn calculateVaporFaceFlux(inputs: VaporFaceInputs) !VaporFaceFlux {
    inline for (@typeInfo(VaporFaceInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSoilWaterFluxInput;
    if (inputs.source_air_volume_m3 < 0 or inputs.destination_air_volume_m3 < 0 or inputs.source_vapor_volume_m3 < 0 or inputs.destination_vapor_volume_m3 < 0 or inputs.source_vapor_diffusivity_m2_per_h < 0 or inputs.destination_vapor_diffusivity_m2_per_h < 0 or inputs.source_air_fraction < 0 or inputs.destination_air_fraction < 0 or inputs.source_porosity_fraction <= 0 or inputs.destination_porosity_fraction <= 0 or inputs.tortuosity < 0 or inputs.source_path_length_m <= 0 or inputs.destination_path_length_m <= 0 or inputs.face_area_m2 < 0 or inputs.time_fraction <= 0 or inputs.time_fraction > 1) return error.InvalidSoilWaterFluxInput;
    if (inputs.source_air_volume_m3 == 0 or inputs.destination_air_volume_m3 == 0) return .{ .unlimited_vapor_m3 = 0, .limited_vapor_m3 = 0, .conductance_m_per_h = 0 };
    const source_concentration = @max(0.0, inputs.source_vapor_volume_m3 / inputs.source_air_volume_m3);
    const destination_concentration = @max(0.0, inputs.destination_vapor_volume_m3 / inputs.destination_air_volume_m3);
    const source_conductivity = inputs.source_vapor_diffusivity_m2_per_h * inputs.tortuosity * inputs.source_air_fraction * inputs.source_air_fraction / inputs.source_porosity_fraction;
    const destination_conductivity = inputs.destination_vapor_diffusivity_m2_per_h * inputs.tortuosity * inputs.destination_air_fraction * inputs.destination_air_fraction / inputs.destination_porosity_fraction;
    const denominator = source_conductivity * inputs.destination_path_length_m + destination_conductivity * inputs.source_path_length_m;
    const conductance = if (denominator > 0) 2.0 * source_conductivity * destination_conductivity / denominator else 0.0;
    const unlimited = conductance * (source_concentration - destination_concentration) * inputs.face_area_m2 * inputs.time_fraction;
    const limited = if (unlimited >= 0)
        @max(0.0, @min(unlimited, inputs.source_vapor_volume_m3))
    else
        @min(0.0, @max(unlimited, -inputs.destination_vapor_volume_m3));
    if (!std.math.isFinite(unlimited) or !std.math.isFinite(limited)) return error.NonFiniteSoilWaterFlux;
    return .{ .unlimited_vapor_m3 = unlimited, .limited_vapor_m3 = limited, .conductance_m_per_h = conductance };
}

fn validateMatrixInputs(inputs: MatrixFaceInputs) !void {
    inline for (@typeInfo(MatrixFaceInputs).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSoilWaterFluxInput;
    }
    if (inputs.source_water_m3 < 0 or inputs.destination_water_m3 < 0 or inputs.source_air_m3 < 0 or inputs.destination_air_m3 < 0 or inputs.source_micropore_volume_m3 <= 0 or inputs.destination_micropore_volume_m3 <= 0 or inputs.source_water_fraction < 0 or inputs.destination_water_fraction < 0 or inputs.source_air_entry_water_fraction < 0 or inputs.destination_air_entry_water_fraction < 0 or inputs.source_hydraulic_conductivity_m2_per_h_mpa < 0 or inputs.destination_hydraulic_conductivity_m2_per_h_mpa < 0 or inputs.source_path_length_m <= 0 or inputs.destination_path_length_m <= 0 or inputs.face_area_m2 < 0 or inputs.wetting_depth_factor < 0 or inputs.time_fraction <= 0 or inputs.time_fraction > 1) return error.InvalidSoilWaterFluxInput;
}

test "matrix face reproduces WATSUB conductance enhancement and bounds" {
    const flux = try calculateMatrixFaceFlux(.{ .direction = .horizontal, .source_water_m3 = 0.30, .destination_water_m3 = 0.10, .source_air_m3 = 0.20, .destination_air_m3 = 0.40, .source_micropore_volume_m3 = 1, .destination_micropore_volume_m3 = 1, .source_water_fraction = 0.30, .destination_water_fraction = 0.10, .source_air_entry_water_fraction = 0.20, .destination_air_entry_water_fraction = 0.20, .source_total_water_potential_mpa = -0.1, .destination_total_water_potential_mpa = -0.2, .source_hydraulic_conductivity_m2_per_h_mpa = 1, .destination_hydraulic_conductivity_m2_per_h_mpa = 1, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .wetting_depth_factor = 0.5, .time_fraction = 1 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), flux.unlimited_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), flux.wetting_enhanced_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), flux.limited_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), flux.transport_water_m3, 1e-12);
}

test "vertical macropore flux is downward only" {
    const common: MacroporeFaceInputs = .{ .direction = .vertical, .source_water_m3 = 0.1, .destination_water_m3 = 0.1, .source_air_m3 = 0.1, .destination_air_m3 = 0.1, .source_total_volume_m3 = 0.2, .destination_total_volume_m3 = 0.2, .source_gravitational_potential_mpa = 0, .destination_gravitational_potential_mpa = 1, .source_vertical_thickness_m = 0.1, .destination_vertical_thickness_m = 0.1, .hydraulic_conductance_m_per_h_mpa = 1, .face_area_m2 = 1, .time_fraction = 1 };
    const flux = try calculateMacroporeFaceFlux(common);
    try std.testing.expect(flux.unlimited_water_m3 < 0);
    try std.testing.expectEqual(@as(f64, 0), flux.limited_water_m3);
}

test "vapor diffusion is donor bounded and conserves volume" {
    const flux = try calculateVaporFaceFlux(.{ .source_air_volume_m3 = 1, .destination_air_volume_m3 = 1, .source_vapor_volume_m3 = 0.01, .destination_vapor_volume_m3 = 0, .source_vapor_diffusivity_m2_per_h = 10, .destination_vapor_diffusivity_m2_per_h = 10, .source_air_fraction = 0.5, .destination_air_fraction = 0.5, .source_porosity_fraction = 0.5, .destination_porosity_fraction = 0.5, .tortuosity = 1, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .time_fraction = 1 });
    try std.testing.expect(flux.unlimited_vapor_m3 > 0.01);
    try std.testing.expectEqual(@as(f64, 0.01), flux.limited_vapor_m3);
}

test "Gerke van Genuchten exchange follows pressure head and donor bounds" {
    const parameters: retention.MualemVanGenuchtenParameters = .{
        .residual_water_content_m3_per_m3 = 0.05,
        .saturated_water_content_m3_per_m3 = 0.45,
        .alpha_per_m = 3.6,
        .n = 1.56,
        .saturated_hydraulic_conductivity_m_per_h = 0.01,
    };
    const matrix_head_m = -2;
    const macropore_head_m = -0.1;
    const interface_conductivity_m_per_h = 0.5 *
        (try parameters.hydraulicConductivityMPerH(matrix_head_m) +
            try parameters.hydraulicConductivityMPerH(macropore_head_m));
    const expected_unbounded_m3 = 3.0 * 0.4 *
        interface_conductivity_m_per_h / (0.2 * 0.2) *
        (macropore_head_m - matrix_head_m) * 1.5;
    const actual = try calculateDualDomainExchange(.{
        .matrix_parameters = parameters,
        .matrix_pressure_head_m = matrix_head_m,
        .macropore_pressure_head_m = macropore_head_m,
        .matrix_bulk_volume_m3 = 1.5,
        .characteristic_matrix_length_m = 0.2,
        .geometry_factor = 3,
        .scaling_coefficient = 0.4,
        .time_fraction = 1,
        .matrix_water_m3 = 0.2,
        .matrix_air_m3 = 0.03,
        .macropore_water_m3 = 0.04,
        .macropore_air_m3 = 0.01,
    });
    try std.testing.expect(actual > 0);
    try std.testing.expectApproxEqAbs(
        @min(expected_unbounded_m3, 0.03),
        actual,
        1e-14,
    );
}
