const std = @import("std");

pub const Face = struct {
    first_cell: usize,
    second_cell: usize,
    /// Positive water moves first -> second; negative moves second -> first.
    water_flux_m3_per_step: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    water_volume_m3: []f64,
    /// Cell-major extensive aqueous inventories in mol.
    amount_mol: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize) !State {
        if (cell_count == 0) return error.ZeroTransportCellCount;
        if (species_count == 0) return error.ZeroTransportSpeciesCount;
        const water = try allocator.alloc(f64, cell_count);
        errdefer allocator.free(water);
        const amount = try allocator.alloc(f64, try std.math.mul(usize, cell_count, species_count));
        errdefer allocator.free(amount);
        @memset(water, 0);
        @memset(amount, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .species_count = species_count, .water_volume_m3 = water, .amount_mol = amount };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.amount_mol);
        self.allocator.free(self.water_volume_m3);
        self.* = undefined;
    }

    pub fn cellAmounts(self: *State, cell_index: usize) ![]f64 {
        if (cell_index >= self.cell_count) return error.TransportCellIndexOutOfBounds;
        const start = cell_index * self.species_count;
        return self.amount_mol[start .. start + self.species_count];
    }

    pub fn cellAmountsConst(self: *const State, cell_index: usize) ![]const f64 {
        if (cell_index >= self.cell_count) return error.TransportCellIndexOutOfBounds;
        const start = cell_index * self.species_count;
        return self.amount_mol[start .. start + self.species_count];
    }

    /// Rejects corrupted runtime-shaped aqueous state before it can silently
    /// enter a transport solve or an authoritative tile generation.
    pub fn validateFinite(self: *const State) !void {
        if (self.cell_count == 0) return error.ZeroTransportCellCount;
        if (self.species_count == 0) return error.ZeroTransportSpeciesCount;
        if (self.water_volume_m3.len != self.cell_count)
            return error.TransportWaterVolumeDimensionMismatch;
        const expected_amount_count = std.math.mul(
            usize,
            self.cell_count,
            self.species_count,
        ) catch return error.TransportStateDimensionOverflow;
        if (self.amount_mol.len != expected_amount_count)
            return error.TransportAmountDimensionMismatch;

        for (self.water_volume_m3, 0..) |water_volume_m3, cell_index| {
            if (!std.math.isFinite(water_volume_m3))
                return error.NonFiniteTransportWaterVolume;
            if (water_volume_m3 < 0)
                return error.NegativeTransportWaterVolume;
            const first_species = cell_index * self.species_count;
            for (
                self.amount_mol[first_species..][0..self.species_count],
            ) |amount_mol| {
                if (!std.math.isFinite(amount_mol))
                    return error.NonFiniteTransportAmount;
                if (amount_mol < 0)
                    return error.NegativeTransportAmount;
            }
        }
    }
};

pub const FaceParameters = struct {
    /// Legacy `VFLWX`: maximum donor-water fraction transported per solve step.
    maximum_convective_fraction: f64,
};

/// Computes positive first->second molar fluxes for every runtime species.
/// `diffusive_conductance_m3_per_step` includes diffusivity, tortuosity,
/// dispersivity, interface area, and distance. `mobility_fraction` represents
/// band/non-band water fractions and is normally one for shared solutes.
pub fn calculateFaceFluxes(state: *const State, face: Face, diffusive_conductance_m3_per_step: []const f64, mobility_fraction: []const f64, parameters: FaceParameters, output_flux_mol: []f64) !void {
    try validateFaceInputs(state, face, diffusive_conductance_m3_per_step, mobility_fraction, parameters, output_flux_mol);
    const first = try state.cellAmountsConst(face.first_cell);
    const second = try state.cellAmountsConst(face.second_cell);
    const first_water = state.water_volume_m3[face.first_cell];
    const second_water = state.water_volume_m3[face.second_cell];

    for (output_flux_mol, first, second, diffusive_conductance_m3_per_step, mobility_fraction) |*flux, first_amount, second_amount, conductance, mobility| {
        const first_mobile = first_amount * mobility;
        const second_mobile = second_amount * mobility;
        const first_concentration = if (first_water > 0) first_amount / first_water else 0;
        const second_concentration = if (second_water > 0) second_amount / second_water else 0;
        const donor_fraction = if (face.water_flux_m3_per_step >= 0)
            if (first_water > 0) @min(parameters.maximum_convective_fraction, face.water_flux_m3_per_step / first_water) else parameters.maximum_convective_fraction
        else if (second_water > 0)
            @min(parameters.maximum_convective_fraction, -face.water_flux_m3_per_step / second_water)
        else
            parameters.maximum_convective_fraction;
        const convection = if (face.water_flux_m3_per_step >= 0) donor_fraction * first_mobile else -donor_fraction * second_mobile;
        const diffusion = conductance * (first_concentration - second_concentration) * mobility;
        // The old NPH loop limited each small transfer indirectly. The direct
        // solve applies the identical direction but clips the combined extent
        // to the available mobile inventory, preventing negative progression.
        flux.* = std.math.clamp(convection + diffusion, -second_mobile, first_mobile);
        if (!std.math.isFinite(flux.*)) return error.NonFiniteSoluteTransportFlux;
    }
}

/// Applies a previously calculated face vector atomically and conserves every
/// species exactly. Suitable face coloring lets callers execute independent
/// grid faces in parallel without atomics.
pub fn commitFaceFluxes(state: *State, face: Face, flux_mol: []const f64) !void {
    if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell) return error.InvalidTransportFace;
    if (flux_mol.len != state.species_count) return error.TransportSpeciesCountMismatch;
    const first = try state.cellAmounts(face.first_cell);
    const second = try state.cellAmounts(face.second_cell);
    for (flux_mol, first, second) |flux, first_amount, second_amount| {
        if (!std.math.isFinite(flux)) return error.NonFiniteSoluteTransportFlux;
        if (first_amount - flux < -1e-12 or second_amount + flux < -1e-12) return error.InsufficientSoluteForTransport;
    }
    for (flux_mol, first, second) |flux, *first_amount, *second_amount| {
        first_amount.* = @max(0, first_amount.* - flux);
        second_amount.* = @max(0, second_amount.* + flux);
    }
}

/// Exact TRNSFRS micropore/macropore exchange equation. Positive flux moves
/// macropore -> micropore. The exchanging macropore water is capped at 5% of
/// geometric layer volume per transport step (`XFRS=0.05`, legacy `VOLT`).
pub fn calculatePoreExchangeFlux(micropore_amount_mol: f64, macropore_amount_mol: f64, micropore_water_m3: f64, macropore_water_m3: f64, layer_volume_m3: f64, step_fraction: f64) !f64 {
    const values = [_]f64{ micropore_amount_mol, macropore_amount_mol, micropore_water_m3, macropore_water_m3, layer_volume_m3, step_fraction };
    for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPoreExchangeInput;
    if (step_fraction > 1) return error.InvalidPoreExchangeInput;
    if (macropore_water_m3 == 0) return 0;
    const exchanging_macropore_water_m3 = @min(0.05 * layer_volume_m3, macropore_water_m3);
    const combined_water_m3 = micropore_water_m3 + exchanging_macropore_water_m3;
    if (combined_water_m3 == 0) return 0;
    const flux = step_fraction * (macropore_amount_mol * micropore_water_m3 - micropore_amount_mol * exchanging_macropore_water_m3) / combined_water_m3;
    if (!std.math.isFinite(flux)) return error.NonFinitePoreExchangeFlux;
    return std.math.clamp(flux, -micropore_amount_mol, macropore_amount_mol);
}

pub fn commitPoreExchange(micropore_amount_mol: *f64, macropore_amount_mol: *f64, flux_macropore_to_micropore_mol: f64) !void {
    if (!std.math.isFinite(micropore_amount_mol.*) or micropore_amount_mol.* < 0 or !std.math.isFinite(macropore_amount_mol.*) or macropore_amount_mol.* < 0 or !std.math.isFinite(flux_macropore_to_micropore_mol)) return error.InvalidPoreExchangeInput;
    if (micropore_amount_mol.* + flux_macropore_to_micropore_mol < -1e-12 or macropore_amount_mol.* - flux_macropore_to_micropore_mol < -1e-12) return error.InsufficientSoluteForPoreExchange;
    micropore_amount_mol.* = @max(0, micropore_amount_mol.* + flux_macropore_to_micropore_mol);
    macropore_amount_mol.* = @max(0, macropore_amount_mol.* - flux_macropore_to_micropore_mol);
}

fn validateFaceInputs(state: *const State, face: Face, conductance: []const f64, mobility: []const f64, parameters: FaceParameters, output: []f64) !void {
    try state.validateFinite();
    if (face.first_cell >= state.cell_count or face.second_cell >= state.cell_count or face.first_cell == face.second_cell) return error.InvalidTransportFace;
    if (conductance.len != state.species_count or mobility.len != state.species_count or output.len != state.species_count) return error.TransportSpeciesCountMismatch;
    if (!std.math.isFinite(face.water_flux_m3_per_step) or !std.math.isFinite(parameters.maximum_convective_fraction) or parameters.maximum_convective_fraction < 0 or parameters.maximum_convective_fraction > 1) return error.InvalidTransportParameter;
    if (!std.math.isFinite(state.water_volume_m3[face.first_cell]) or state.water_volume_m3[face.first_cell] < 0 or !std.math.isFinite(state.water_volume_m3[face.second_cell]) or state.water_volume_m3[face.second_cell] < 0) return error.InvalidTransportWaterVolume;
    const first = try state.cellAmountsConst(face.first_cell);
    const second = try state.cellAmountsConst(face.second_cell);
    for (first, second, conductance, mobility) |a, b, d, mobile| {
        if (!std.math.isFinite(a) or a < 0 or !std.math.isFinite(b) or b < 0) return error.InvalidSoluteTransportState;
        if (!std.math.isFinite(d) or d < 0 or !std.math.isFinite(mobile) or mobile < 0 or mobile > 1) return error.InvalidTransportParameter;
    }
}

test "upwind convection and diffusion conserve every runtime species" {
    var state = try State.init(std.testing.allocator, 2, 3);
    defer state.deinit();
    state.water_volume_m3[0] = 2;
    state.water_volume_m3[1] = 1;
    @memcpy(try state.cellAmounts(0), &[_]f64{ 4, 2, 8 });
    @memcpy(try state.cellAmounts(1), &[_]f64{ 1, 3, 0 });
    const before = [_]f64{ 5, 5, 8 };
    var flux: [3]f64 = undefined;
    try calculateFaceFluxes(&state, .{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0.2 }, &[_]f64{ 0.1, 0.1, 0.1 }, &[_]f64{ 1, 0.5, 1 }, .{ .maximum_convective_fraction = 0.2 }, &flux);
    try commitFaceFluxes(&state, .{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 0.2 }, &flux);
    const first = try state.cellAmountsConst(0);
    const second = try state.cellAmountsConst(1);
    for (before, first, second) |total, a, b| try std.testing.expectApproxEqAbs(total, a + b, 1e-14);
}

test "negative water flux uses the second cell as donor" {
    var state = try State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    state.water_volume_m3[0] = 1;
    state.water_volume_m3[1] = 2;
    (try state.cellAmounts(0))[0] = 0;
    (try state.cellAmounts(1))[0] = 4;
    var flux: [1]f64 = undefined;
    try calculateFaceFluxes(&state, .{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = -0.2 }, &[_]f64{0}, &[_]f64{1}, .{ .maximum_convective_fraction = 1 }, &flux);
    try std.testing.expectApproxEqAbs(@as(f64, -0.4), flux[0], 1e-15);
}

test "micropore macropore exchange matches XFRS equation" {
    const flux = try calculatePoreExchangeFlux(2, 6, 4, 3, 20, 0.25);
    const expected: f64 = 0.25 * (6.0 * 4.0 - 2.0 * 1.0) / 5.0;
    try std.testing.expectApproxEqAbs(expected, flux, 1e-15);
    var micro: f64 = 2;
    var macro: f64 = 6;
    try commitPoreExchange(&micro, &macro, flux);
    try std.testing.expectApproxEqAbs(@as(f64, 8), micro + macro, 1e-15);
}

test "combined transport cannot overdraw a donor" {
    var state = try State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    state.water_volume_m3[0] = 1;
    state.water_volume_m3[1] = 1;
    (try state.cellAmounts(0))[0] = 0.1;
    var flux: [1]f64 = undefined;
    try calculateFaceFluxes(&state, .{ .first_cell = 0, .second_cell = 1, .water_flux_m3_per_step = 100 }, &[_]f64{100}, &[_]f64{1}, .{ .maximum_convective_fraction = 1 }, &flux);
    try std.testing.expectEqual(@as(f64, 0.1), flux[0]);
}

test "aqueous state supports more than five runtime species and validates all inventories" {
    const runtime_species_count: usize = 13;
    var state = try State.init(
        std.testing.allocator,
        4,
        runtime_species_count,
    );
    defer state.deinit();

    for (state.water_volume_m3, 0..) |*water_volume_m3, cell_index|
        water_volume_m3.* = @as(f64, @floatFromInt(cell_index + 1));
    for (state.amount_mol, 0..) |*amount_mol, component_index|
        amount_mol.* = @as(f64, @floatFromInt(component_index + 1)) * 1e-6;
    try state.validateFinite();

    state.amount_mol[runtime_species_count + 7] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteTransportAmount,
        state.validateFinite(),
    );
}
