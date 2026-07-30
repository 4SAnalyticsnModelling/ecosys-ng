const std = @import("std");
const audit = @import("mass_balance_audit.zig");
const gas = @import("gas_transport.zig");
const soil_daily_gas = @import("soil_daily_gas_flux.zig");
const canopy_daily_gas = @import("daily_canopy_gas_exchange.zig");
const atmospheric_solutes = @import("atmospheric_solute_input_ledger.zig");
const snow = @import("snow_solute_transport.zig");
const snow_discharge = @import("snow_surface_discharge.zig");

/// Direction-separated whole-landscape boundary fluxes for one accepted
/// model interval. Every value is nonnegative; the EXEC sign convention is
/// applied only by `mass_balance_audit.balance`.
pub const Fluxes = struct {
    rain_m3: f64 = 0,
    runoff_m3: f64 = 0,
    evaporation_m3: f64 = 0,
    water_outflow_m3: f64 = 0,
    heat_input_mj: f64 = 0,
    heat_output_mj: f64 = 0,
    oxygen_input_g: f64 = 0,
    oxygen_output_g: f64 = 0,
    carbon_dioxide_input_g_c: f64 = 0,
    carbon_output_g_c: f64 = 0,
    organic_fertilizer_carbon_g_c: f64 = 0,
    carbon_sink_g_c: f64 = 0,
    dinitrogen_input_g_n: f64 = 0,
    nitrogen_input_g_n: f64 = 0,
    nitrogen_output_g_n: f64 = 0,
    organic_fertilizer_nitrogen_g_n: f64 = 0,
    nitrogen_sink_g_n: f64 = 0,
    phosphorus_input_g_p: f64 = 0,
    phosphorus_output_g_p: f64 = 0,
    organic_fertilizer_phosphorus_g_p: f64 = 0,
    phosphorus_sink_g_p: f64 = 0,
    ion_input_mol: f64 = 0,
    ion_output_mol: f64 = 0,
};

// Deliberately no biological-fixation boundary field exists here.
//
// EXEC excludes living plant biomass from its reconstructed storage. GROSUB
// symbiotic fixation is therefore confined to the separate plant balance;
// fixed N entering REDIST residue is paired with the XZSN plant-litter sink.
// NITRO nonsymbiotic fixation is an internal transfer: the accepted soil and
// surface commits remove exactly the fixed mass from dissolved N2 while
// adding it to microbial organic N. Treating either process as an additional
// atmospheric boundary input would double-count nitrogen.

/// Accepted cumulative EXEC boundary history. The transaction is staged so a
/// non-finite or negative process flux cannot partially advance the ledger.
pub const State = struct {
    cumulative: Fluxes = .{},

    pub fn accumulateAccepted(self: *State, fluxes: Fluxes) !void {
        try validate(fluxes);
        var next = self.cumulative;
        inline for (std.meta.fields(Fluxes)) |field| {
            const value = @field(next, field.name) + @field(fluxes, field.name);
            if (!std.math.isFinite(value))
                return error.LandscapeBoundaryLedgerOverflow;
            @field(next, field.name) = value;
        }
        self.cumulative = next;
    }

    pub fn reset(self: *State) void {
        self.cumulative = .{};
    }

    pub fn accumulateAcceptedWater(
        self: *State,
        rainfall_m3: []const f64,
        runoff_m3: []const f64,
        evaporation_m3: []const f64,
        water_outflow_m3: []const f64,
        artificial_drainage_outflow_m3: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            rainfall_m3,
            runoff_m3,
            evaporation_m3,
            water_outflow_m3,
            artificial_drainage_outflow_m3,
        });
        try self.accumulateAccepted(.{
            .rain_m3 = try sumNonnegative(rainfall_m3),
            .runoff_m3 = try sumNonnegative(runoff_m3),
            .evaporation_m3 = try sumNonnegative(evaporation_m3),
            .water_outflow_m3 = try addFinite(
                try sumNonnegative(water_outflow_m3),
                try sumNonnegative(artificial_drainage_outflow_m3),
            ),
        });
    }

    /// Atmospheric sensible heat carried into the landscape by liquid
    /// precipitation (including irrigation already merged with rain) and
    /// snowfall water equivalent. This is evaluated before canopy/snow/
    /// litter/soil routing so intercepted water is included exactly once.
    pub fn accumulateAcceptedPrecipitationHeat(
        self: *State,
        liquid_water_depth_m: []const f64,
        snowfall_water_equivalent_depth_m: []const f64,
        cell_area_m2: []const f64,
        atmospheric_temperature_k: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            liquid_water_depth_m,
            snowfall_water_equivalent_depth_m,
            cell_area_m2,
            atmospheric_temperature_k,
        });
        var heat_input_mj: f64 = 0;
        for (
            liquid_water_depth_m,
            snowfall_water_equivalent_depth_m,
            cell_area_m2,
            atmospheric_temperature_k,
        ) |liquid_depth_m, snow_depth_m, area_m2, temperature_k| {
            inline for (.{ liquid_depth_m, snow_depth_m, area_m2, temperature_k }) |value|
                if (!std.math.isFinite(value))
                    return error.NonFiniteLandscapeBoundaryFlux;
            if (liquid_depth_m < 0 or snow_depth_m < 0 or area_m2 <= 0 or temperature_k <= 0)
                return error.InvalidPrecipitationHeatBoundaryInput;
            const cell_heat_mj = temperature_k * area_m2 *
                (4.19 * liquid_depth_m + 2.095 * snow_depth_m);
            if (!std.math.isFinite(cell_heat_mj))
                return error.LandscapeBoundaryLedgerOverflow;
            heat_input_mj = try addFinite(heat_input_mj, cell_heat_mj);
        }
        try self.accumulateAccepted(.{ .heat_input_mj = heat_input_mj });
    }

    /// Source HEATH plus THFLXC boundary accounting. Ground flux components
    /// are positive into the surface and expressed per horizontal area.
    /// Canopy water-energy changes are already extensive signed MJ and retain
    /// EXTRACT's exact current-minus-previous-minus-retained-rain equation.
    pub fn accumulateAcceptedSurfaceAndCanopyHeat(
        self: *State,
        ground_net_radiation_mj_per_m2: []const f64,
        ground_sensible_heat_mj_per_m2: []const f64,
        ground_latent_heat_mj_per_m2: []const f64,
        ground_vapor_sensible_heat_mj_per_m2: []const f64,
        cell_area_m2: []const f64,
        canopy_water_energy_change_mj: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            ground_net_radiation_mj_per_m2,
            ground_sensible_heat_mj_per_m2,
            ground_latent_heat_mj_per_m2,
            ground_vapor_sensible_heat_mj_per_m2,
            cell_area_m2,
        });
        var signed_heat_into_landscape_mj: f64 = 0;
        for (
            ground_net_radiation_mj_per_m2,
            ground_sensible_heat_mj_per_m2,
            ground_latent_heat_mj_per_m2,
            ground_vapor_sensible_heat_mj_per_m2,
            cell_area_m2,
        ) |net_radiation, sensible, latent, vapor_sensible, area_m2| {
            inline for (.{ net_radiation, sensible, latent, vapor_sensible, area_m2 }) |value|
                if (!std.math.isFinite(value))
                    return error.NonFiniteLandscapeBoundaryFlux;
            if (area_m2 <= 0) return error.InvalidSurfaceHeatBoundaryArea;
            signed_heat_into_landscape_mj = try addFinite(
                signed_heat_into_landscape_mj,
                (net_radiation + sensible + latent + vapor_sensible) * area_m2,
            );
        }
        for (canopy_water_energy_change_mj) |change_mj| {
            if (!std.math.isFinite(change_mj))
                return error.NonFiniteLandscapeBoundaryFlux;
            signed_heat_into_landscape_mj =
                try addFinite(signed_heat_into_landscape_mj, change_mj);
        }
        try self.accumulateAccepted(.{
            .heat_input_mj = @max(0, signed_heat_into_landscape_mj),
            .heat_output_mj = @max(0, -signed_heat_into_landscape_mj),
        });
    }

    pub fn accumulateAcceptedSignedHeat(
        self: *State,
        signed_heat_into_landscape_mj: f64,
    ) !void {
        if (!std.math.isFinite(signed_heat_into_landscape_mj))
            return error.NonFiniteLandscapeBoundaryFlux;
        try self.accumulateAccepted(.{
            .heat_input_mj = @max(0, signed_heat_into_landscape_mj),
            .heat_output_mj = @max(0, -signed_heat_into_landscape_mj),
        });
    }

    pub fn accumulateAcceptedFertilizer(
        self: *State,
        mineral_nitrogen_g_n: []const f64,
        organic_nitrogen_g_n: []const f64,
        mineral_phosphorus_g_p: []const f64,
        organic_phosphorus_g_p: []const f64,
        organic_carbon_g_c: []const f64,
        biome_organic_carbon_g_c: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            mineral_nitrogen_g_n,
            organic_nitrogen_g_n,
            mineral_phosphorus_g_p,
            organic_phosphorus_g_p,
            organic_carbon_g_c,
            biome_organic_carbon_g_c,
        });
        try self.accumulateAccepted(.{
            .nitrogen_input_g_n = try sumNonnegative(mineral_nitrogen_g_n),
            .organic_fertilizer_nitrogen_g_n = try sumNonnegative(organic_nitrogen_g_n),
            .phosphorus_input_g_p = try sumNonnegative(mineral_phosphorus_g_p),
            .organic_fertilizer_phosphorus_g_p = try sumNonnegative(organic_phosphorus_g_p),
            .organic_fertilizer_carbon_g_c = try addFinite(
                try sumNonnegative(organic_carbon_g_c),
                try sumNonnegative(biome_organic_carbon_g_c),
            ),
        });
    }

    pub fn accumulateAcceptedExports(
        self: *State,
        carbon_exports_g_c: [4][]const f64,
        nitrogen_exports_g_n: [4][]const f64,
        phosphorus_exports_g_p: [4][]const f64,
        ion_outflow_mol: []const f64,
    ) !void {
        const cell_count = ion_outflow_mol.len;
        if (cell_count == 0) return error.EmptyLandscapeBoundaryGrid;
        inline for (carbon_exports_g_c ++ nitrogen_exports_g_n ++
            phosphorus_exports_g_p) |values|
            if (values.len != cell_count)
                return error.LandscapeBoundaryGridDimensionMismatch;

        var carbon_output_g_c: f64 = 0;
        var nitrogen_output_g_n: f64 = 0;
        var phosphorus_output_g_p: f64 = 0;
        inline for (carbon_exports_g_c) |values|
            carbon_output_g_c = try addFinite(
                carbon_output_g_c,
                try sumNonnegative(values),
            );
        inline for (nitrogen_exports_g_n) |values|
            nitrogen_output_g_n = try addFinite(
                nitrogen_output_g_n,
                try sumNonnegative(values),
            );
        inline for (phosphorus_exports_g_p) |values|
            phosphorus_output_g_p = try addFinite(
                phosphorus_output_g_p,
                try sumNonnegative(values),
            );
        try self.accumulateAccepted(.{
            .carbon_output_g_c = carbon_output_g_c,
            .nitrogen_output_g_n = nitrogen_output_g_n,
            .phosphorus_output_g_p = phosphorus_output_g_p,
            .ion_output_mol = try sumNonnegative(ion_outflow_mol),
        });
    }

    /// EXEC XCSN/XZSN/XPSN: plant litter transferred into the reconstructed
    /// residue stores. These cumulative sinks remove that internal transfer
    /// from the whole-landscape equation because living plant biomass is not
    /// included in REDIST storage.
    pub fn accumulateAcceptedPlantLitter(
        self: *State,
        carbon_litter_g_c: []const f64,
        nitrogen_litter_g_n: []const f64,
        phosphorus_litter_g_p: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            carbon_litter_g_c,
            nitrogen_litter_g_n,
            phosphorus_litter_g_p,
        });
        try self.accumulateAccepted(.{
            .carbon_sink_g_c = try sumNonnegative(carbon_litter_g_c),
            .nitrogen_sink_g_n = try sumNonnegative(nitrogen_litter_g_n),
            .phosphorus_sink_g_p = try sumNonnegative(phosphorus_litter_g_p),
        });
    }

    /// REDIST atmospheric gas boundary terms. Soil/litter DAY flux is
    /// positive ecosystem -> atmosphere; canopy EXTRACT flux is positive
    /// atmosphere -> ecosystem. Each species is split independently so
    /// simultaneous uptake and emission are retained instead of netted away.
    pub fn accumulateAcceptedAtmosphericGas(
        self: *State,
        soil: *const soil_daily_gas.State,
        canopy: *const canopy_daily_gas.State,
    ) !void {
        if (soil.cell_count == 0 or
            soil.tracked_element_mass_g_by_cell_and_species.len !=
                soil.cell_count * gas.species_count or
            canopy.net_carbon_dioxide_uptake_g_c.len != soil.cell_count or
            canopy.net_methane_uptake_g_c.len != soil.cell_count or
            canopy.net_oxygen_uptake_g_o.len != soil.cell_count)
            return error.LandscapeGasBoundaryDimensionMismatch;

        var carbon_input_g_c: f64 = 0;
        var carbon_output_g_c: f64 = 0;
        var oxygen_input_g: f64 = 0;
        var oxygen_output_g: f64 = 0;
        var gaseous_nitrogen_input_g_n: f64 = 0;
        var gaseous_nitrogen_output_g_n: f64 = 0;
        for (0..soil.cell_count) |cell| {
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                canopy.net_carbon_dioxide_uptake_g_c[cell],
            );
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                canopy.net_methane_uptake_g_c[cell],
            );
            try splitSignedAtmosphereFlux(
                &oxygen_input_g,
                &oxygen_output_g,
                canopy.net_oxygen_uptake_g_o[cell],
            );
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                -(try soil.get(cell, .carbon_dioxide)),
            );
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                -(try soil.get(cell, .methane)),
            );
            try splitSignedAtmosphereFlux(
                &oxygen_input_g,
                &oxygen_output_g,
                -(try soil.get(cell, .oxygen)),
            );
            inline for (.{
                gas.Species.nitrogen,
                gas.Species.nitrous_oxide,
                gas.Species.ammonia,
            }) |species|
                try splitSignedAtmosphereFlux(
                    &gaseous_nitrogen_input_g_n,
                    &gaseous_nitrogen_output_g_n,
                    -(try soil.get(cell, species)),
                );
        }
        try self.accumulateAccepted(.{
            .carbon_dioxide_input_g_c = carbon_input_g_c,
            .carbon_output_g_c = carbon_output_g_c,
            .oxygen_input_g = oxygen_input_g,
            .oxygen_output_g = oxygen_output_g,
            .dinitrogen_input_g_n = gaseous_nitrogen_input_g_n,
            .nitrogen_output_g_n = gaseous_nitrogen_output_g_n,
        });
    }

    /// REDIST precipitation/irrigation chemistry after snow/direct routing
    /// has been recombined. Nutrients retain tracked-element grams; the eight
    /// salt carriers are converted to mol using runtime chemistry constants.
    pub fn accumulateAcceptedAtmosphericSolutes(
        self: *State,
        inputs: *const atmospheric_solutes.State,
        ion_molar_mass_g_per_mol: snow_discharge.IonMolarMassesGPerMol,
    ) !void {
        if (inputs.cell_count == 0 or inputs.daily_input_g.len !=
            inputs.cell_count * snow.species_count)
            return error.AtmosphericSoluteInputDimensionMismatch;
        inline for (std.meta.fields(snow_discharge.IonMolarMassesGPerMol)) |field| {
            const value = @field(ion_molar_mass_g_per_mol, field.name);
            if (!std.math.isFinite(value) or value <= 0)
                return error.InvalidIonMolarMass;
        }

        var carbon_input_g_c: f64 = 0;
        var oxygen_input_g: f64 = 0;
        var gaseous_nitrogen_input_g_n: f64 = 0;
        var aqueous_nitrogen_input_g_n: f64 = 0;
        var phosphorus_input_g_p: f64 = 0;
        var ion_input_mol: f64 = 0;
        for (0..inputs.cell_count) |cell| {
            const first = cell * snow.species_count;
            const amounts = inputs.daily_input_g[first .. first + snow.species_count];
            for (amounts) |amount|
                if (!std.math.isFinite(amount) or amount < 0)
                    return error.InvalidAtmosphericSoluteInput;
            carbon_input_g_c = try addFinite(
                carbon_input_g_c,
                amounts[@intFromEnum(snow.Species.carbon_dioxide_carbon)] +
                    amounts[@intFromEnum(snow.Species.methane_carbon)],
            );
            oxygen_input_g = try addFinite(
                oxygen_input_g,
                amounts[@intFromEnum(snow.Species.oxygen)],
            );
            gaseous_nitrogen_input_g_n = try addFinite(
                gaseous_nitrogen_input_g_n,
                amounts[@intFromEnum(snow.Species.dinitrogen_nitrogen)] +
                    amounts[
                        @intFromEnum(
                            snow.Species.nitrous_oxide_nitrogen,
                        )
                    ],
            );
            aqueous_nitrogen_input_g_n = try addFinite(
                aqueous_nitrogen_input_g_n,
                amounts[@intFromEnum(snow.Species.ammonium_nitrogen)] +
                    amounts[@intFromEnum(snow.Species.ammonia_nitrogen)] +
                    amounts[@intFromEnum(snow.Species.nitrate_nitrogen)],
            );
            const hpo4_g_p = amounts[
                @intFromEnum(snow.Species.hydrogen_phosphate_phosphorus)
            ];
            const h2po4_g_p = amounts[
                @intFromEnum(snow.Species.dihydrogen_phosphate_phosphorus)
            ];
            phosphorus_input_g_p = try addFinite(
                phosphorus_input_g_p,
                hpo4_g_p + h2po4_g_p,
            );
            ion_input_mol = try addFinite(
                ion_input_mol,
                3 * hpo4_g_p / snow.phosphorus_g_per_mol +
                    4 * h2po4_g_p / snow.phosphorus_g_per_mol +
                    amounts[@intFromEnum(snow.Species.aluminum)] /
                        ion_molar_mass_g_per_mol.aluminum +
                    amounts[@intFromEnum(snow.Species.iron)] /
                        ion_molar_mass_g_per_mol.iron +
                    amounts[@intFromEnum(snow.Species.calcium)] /
                        ion_molar_mass_g_per_mol.calcium +
                    amounts[@intFromEnum(snow.Species.magnesium)] /
                        ion_molar_mass_g_per_mol.magnesium +
                    amounts[@intFromEnum(snow.Species.sodium)] /
                        ion_molar_mass_g_per_mol.sodium +
                    amounts[@intFromEnum(snow.Species.potassium)] /
                        ion_molar_mass_g_per_mol.potassium +
                    amounts[@intFromEnum(snow.Species.sulfate_sulfur)] /
                        ion_molar_mass_g_per_mol.sulfur +
                    amounts[@intFromEnum(snow.Species.chloride)] /
                        ion_molar_mass_g_per_mol.chloride,
            );
        }
        try self.accumulateAccepted(.{
            .carbon_dioxide_input_g_c = carbon_input_g_c,
            .oxygen_input_g = oxygen_input_g,
            .dinitrogen_input_g_n = gaseous_nitrogen_input_g_n,
            .nitrogen_input_g_n = aqueous_nitrogen_input_g_n,
            .phosphorus_input_g_p = phosphorus_input_g_p,
            .ion_input_mol = ion_input_mol,
        });
    }

    /// Publishes boundary fields only. Area and reconstructed storage remain
    /// owned by the runtime landscape inventory.
    pub fn publish(self: State, totals: *audit.Totals) !void {
        try validate(self.cumulative);
        const f = self.cumulative;
        totals.cumulative_rain_m3 = f.rain_m3;
        totals.cumulative_runoff_m3 = f.runoff_m3;
        totals.cumulative_evaporation_m3 = f.evaporation_m3;
        totals.cumulative_water_outflow_m3 = f.water_outflow_m3;
        totals.cumulative_heat_input_mj = f.heat_input_mj;
        totals.cumulative_heat_output_mj = f.heat_output_mj;
        totals.cumulative_oxygen_input_g = f.oxygen_input_g;
        totals.cumulative_oxygen_output_g = f.oxygen_output_g;
        totals.cumulative_carbon_dioxide_input_g =
            f.carbon_dioxide_input_g_c;
        totals.cumulative_carbon_output_g = f.carbon_output_g_c;
        totals.cumulative_organic_fertilizer_carbon_g =
            f.organic_fertilizer_carbon_g_c;
        totals.cumulative_carbon_sink_g = f.carbon_sink_g_c;
        totals.cumulative_dinitrogen_input_g = f.dinitrogen_input_g_n;
        totals.cumulative_nitrogen_input_g = f.nitrogen_input_g_n;
        totals.cumulative_nitrogen_output_g = f.nitrogen_output_g_n;
        totals.cumulative_organic_fertilizer_nitrogen_g =
            f.organic_fertilizer_nitrogen_g_n;
        totals.cumulative_nitrogen_sink_g = f.nitrogen_sink_g_n;
        totals.cumulative_phosphorus_input_g = f.phosphorus_input_g_p;
        totals.cumulative_phosphorus_output_g = f.phosphorus_output_g_p;
        totals.cumulative_organic_fertilizer_phosphorus_g =
            f.organic_fertilizer_phosphorus_g_p;
        totals.cumulative_phosphorus_sink_g = f.phosphorus_sink_g_p;
        totals.cumulative_ion_input_mol = f.ion_input_mol;
        totals.cumulative_ion_output_mol = f.ion_output_mol;
    }
};

fn validate(fluxes: Fluxes) !void {
    inline for (std.meta.fields(Fluxes)) |field| {
        const value = @field(fluxes, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteLandscapeBoundaryFlux;
        if (value < 0) return error.NegativeLandscapeBoundaryFlux;
    }
}

fn requireSameNonzeroLength(slices: anytype) !void {
    const length = slices[0].len;
    if (length == 0) return error.EmptyLandscapeBoundaryGrid;
    inline for (slices) |values|
        if (values.len != length)
            return error.LandscapeBoundaryGridDimensionMismatch;
}

fn sumNonnegative(values: []const f64) !f64 {
    var result: f64 = 0;
    for (values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteLandscapeBoundaryFlux;
        if (value < 0) return error.NegativeLandscapeBoundaryFlux;
        result = try addFinite(result, value);
    }
    return result;
}

fn addFinite(first: f64, second: f64) !f64 {
    const result = first + second;
    if (!std.math.isFinite(result))
        return error.LandscapeBoundaryLedgerOverflow;
    return result;
}

fn splitSignedAtmosphereFlux(
    input: *f64,
    output: *f64,
    atmosphere_to_ecosystem_g: f64,
) !void {
    if (!std.math.isFinite(atmosphere_to_ecosystem_g))
        return error.NonFiniteLandscapeBoundaryFlux;
    if (atmosphere_to_ecosystem_g >= 0)
        input.* = try addFinite(input.*, atmosphere_to_ecosystem_g)
    else
        output.* = try addFinite(output.*, -atmosphere_to_ecosystem_g);
}

test "accepted landscape boundary transactions publish every EXEC ledger" {
    var state: State = .{};
    var first = std.mem.zeroes(Fluxes);
    inline for (std.meta.fields(Fluxes), 0..) |field, index|
        @field(first, field.name) = @floatFromInt(index + 1);
    try state.accumulateAccepted(first);
    try state.accumulateAccepted(first);

    var totals = std.mem.zeroes(audit.Totals);
    totals.landscape_area_m2 = 99;
    totals.water_storage_m3 = 77;
    try state.publish(&totals);
    try std.testing.expectEqual(@as(f64, 2), totals.cumulative_rain_m3);
    try std.testing.expectEqual(@as(f64, 46), totals.cumulative_ion_output_mol);
    try std.testing.expectEqual(@as(f64, 99), totals.landscape_area_m2);
    try std.testing.expectEqual(@as(f64, 77), totals.water_storage_m3);
}

test "invalid boundary transaction cannot partially advance cumulative state" {
    var state: State = .{};
    try state.accumulateAccepted(.{ .rain_m3 = 3, .nitrogen_input_g_n = 4 });
    const before = state.cumulative;
    try std.testing.expectError(
        error.NonFiniteLandscapeBoundaryFlux,
        state.accumulateAccepted(.{
            .rain_m3 = 5,
            .phosphorus_output_g_p = std.math.nan(f64),
        }),
    );
    try std.testing.expectEqualDeep(before, state.cumulative);
}

test "accepted runtime grid ledgers aggregate water fertilizer and exports" {
    var state: State = .{};
    const first = [_]f64{ 1, 2 };
    const second = [_]f64{ 3, 4 };
    const third = [_]f64{ 5, 6 };
    const fourth = [_]f64{ 7, 8 };
    try state.accumulateAcceptedWater(
        &first,
        &second,
        &third,
        &fourth,
        &first,
    );
    try state.accumulateAcceptedFertilizer(
        &first,
        &second,
        &third,
        &fourth,
        &first,
        &second,
    );
    try state.accumulateAcceptedExports(
        .{ &first, &second, &third, &fourth },
        .{ &first, &second, &third, &fourth },
        .{ &first, &second, &third, &fourth },
        &second,
    );
    try state.accumulateAcceptedPlantLitter(&first, &second, &third);
    try std.testing.expectEqual(@as(f64, 3), state.cumulative.rain_m3);
    try std.testing.expectEqual(
        @as(f64, 36),
        state.cumulative.carbon_output_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 7),
        state.cumulative.ion_output_mol,
    );
    try std.testing.expectEqual(
        @as(f64, 3),
        state.cumulative.carbon_sink_g_c,
    );
}

test "precipitation heat counts raw liquid and snow once across runtime cells" {
    var state: State = .{};
    try state.accumulateAcceptedPrecipitationHeat(
        &.{ 0.001, 0.002 },
        &.{ 0.003, 0.004 },
        &.{ 10, 20 },
        &.{ 280, 285 },
    );
    const expected =
        280 * 10 * (4.19 * 0.001 + 2.095 * 0.003) +
        285 * 20 * (4.19 * 0.002 + 2.095 * 0.004);
    try std.testing.expectApproxEqAbs(
        expected,
        state.cumulative.heat_input_mj,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), state.cumulative.heat_output_mj);
}

test "invalid precipitation heat transaction cannot advance ledger" {
    var state: State = .{};
    try state.accumulateAccepted(.{ .heat_input_mj = 7 });
    const before = state.cumulative;
    try std.testing.expectError(
        error.InvalidPrecipitationHeatBoundaryInput,
        state.accumulateAcceptedPrecipitationHeat(
            &.{0.001},
            &.{0.002},
            &.{0},
            &.{280},
        ),
    );
    try std.testing.expectEqualDeep(before, state.cumulative);
}

test "surface HEATH and canopy THFLXC retain signed boundary direction" {
    var state: State = .{};
    try state.accumulateAcceptedSurfaceAndCanopyHeat(
        &.{ 2, -1 },
        &.{ 0.5, -0.5 },
        &.{ -0.25, 0.25 },
        &.{ -0.05, 0.05 },
        &.{ 10, 20 },
        &.{ -2, 1 },
    );
    // Ground: 22 MJ + (-24 MJ); canopy: -1 MJ.
    try std.testing.expectEqual(@as(f64, 0), state.cumulative.heat_input_mj);
    try std.testing.expectApproxEqAbs(
        @as(f64, 3),
        state.cumulative.heat_output_mj,
        1.0e-12,
    );
}

test "invalid late canopy heat leaves boundary ledger unchanged" {
    var state: State = .{};
    try state.accumulateAccepted(.{ .heat_input_mj = 4 });
    const before = state.cumulative;
    try std.testing.expectError(
        error.NonFiniteLandscapeBoundaryFlux,
        state.accumulateAcceptedSurfaceAndCanopyHeat(
            &.{1},
            &.{2},
            &.{3},
            &.{4},
            &.{5},
            &.{std.math.nan(f64)},
        ),
    );
    try std.testing.expectEqualDeep(before, state.cumulative);
}

test "atmospheric gases retain simultaneous species uptake and emission" {
    var soil = try soil_daily_gas.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    var canopy = try canopy_daily_gas.State.init(std.testing.allocator, 1);
    defer canopy.deinit();
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.carbon_dioxide)
    ] = 3;
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.methane)
    ] = -2;
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.oxygen)
    ] = 4;
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.nitrogen)
    ] = -7;
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.nitrous_oxide)
    ] = 8;
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.ammonia)
    ] = -9;
    canopy.net_carbon_dioxide_uptake_g_c[0] = 5;
    canopy.net_methane_uptake_g_c[0] = -1;
    canopy.net_oxygen_uptake_g_o[0] = 6;

    var state: State = .{};
    try state.accumulateAcceptedAtmosphericGas(&soil, &canopy);
    try std.testing.expectEqual(
        @as(f64, 7),
        state.cumulative.carbon_dioxide_input_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 4),
        state.cumulative.carbon_output_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 6),
        state.cumulative.oxygen_input_g,
    );
    try std.testing.expectEqual(
        @as(f64, 4),
        state.cumulative.oxygen_output_g,
    );
    try std.testing.expectEqual(
        @as(f64, 16),
        state.cumulative.dinitrogen_input_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, 8),
        state.cumulative.nitrogen_output_g_n,
    );
}

test "atmospheric solutes map tracked elements and phosphate atom counts" {
    var inputs = try atmospheric_solutes.State.init(std.testing.allocator, 1);
    defer inputs.deinit();
    const amounts = inputs.daily_input_g[0..snow.species_count];
    amounts[@intFromEnum(snow.Species.carbon_dioxide_carbon)] = 1;
    amounts[@intFromEnum(snow.Species.methane_carbon)] = 2;
    amounts[@intFromEnum(snow.Species.oxygen)] = 3;
    amounts[@intFromEnum(snow.Species.dinitrogen_nitrogen)] = 4;
    amounts[@intFromEnum(snow.Species.nitrous_oxide_nitrogen)] = 5;
    amounts[@intFromEnum(snow.Species.ammonium_nitrogen)] = 6;
    amounts[@intFromEnum(snow.Species.ammonia_nitrogen)] = 7;
    amounts[@intFromEnum(snow.Species.nitrate_nitrogen)] = 8;
    amounts[@intFromEnum(snow.Species.hydrogen_phosphate_phosphorus)] = 31;
    amounts[@intFromEnum(snow.Species.dihydrogen_phosphate_phosphorus)] = 62;
    inline for (0..8) |ion|
        amounts[10 + ion] += @as(f64, @floatFromInt(ion + 1));

    var state: State = .{};
    try state.accumulateAcceptedAtmosphericSolutes(&inputs, .{
        .aluminum = 1,
        .iron = 2,
        .calcium = 3,
        .magnesium = 4,
        .sodium = 5,
        .potassium = 6,
        .sulfur = 7,
        .chloride = 8,
    });
    try std.testing.expectEqual(
        @as(f64, 3),
        state.cumulative.carbon_dioxide_input_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 9),
        state.cumulative.dinitrogen_input_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, 21),
        state.cumulative.nitrogen_input_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, 93),
        state.cumulative.phosphorus_input_g_p,
    );
    // Eight free-ion carriers contribute 1 mol each; HPO4 contributes 3
    // atoms and two mol H2PO4 contribute 8.
    try std.testing.expectEqual(
        @as(f64, 19),
        state.cumulative.ion_input_mol,
    );
}
