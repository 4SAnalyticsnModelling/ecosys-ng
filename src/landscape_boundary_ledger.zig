const std = @import("std");
const audit = @import("mass_balance_audit.zig");
const gas = @import("gas_transport.zig");
const soil_daily_gas = @import("soil_daily_gas_flux.zig");
const canopy_daily_gas = @import("daily_canopy_gas_exchange.zig");
const atmospheric_solutes = @import("atmospheric_solute_input_ledger.zig");
const snow = @import("snow_solute_transport.zig");
const snow_discharge = @import("snow_surface_discharge.zig");
/// HEAT-001 third layer. Imported for `frozenWaterEnthalpyPerM3` only, so the
/// boundary credit for snowfall and the census's stored snow enthalpy share
/// one definition rather than two agreeing expressions.
const inventory = @import("landscape_mass_inventory.zig");

/// Direction-separated whole-landscape boundary fluxes for one accepted
/// model interval. Every value is nonnegative; the EXEC sign convention is
/// applied only by `mass_balance_audit.balance`.
pub const Fluxes = struct {
    rain_m3: f64 = 0,
    boundary_water_inflow_m3: f64 = 0,
    runoff_m3: f64 = 0,
    evaporation_m3: f64 = 0,
    water_outflow_m3: f64 = 0,
    heat_input_megajoules: f64 = 0,
    heat_output_megajoules: f64 = 0,
    oxygen_input_g: f64 = 0,
    oxygen_output_g: f64 = 0,
    redist_carbon_surface_input_g_c: f64 = 0,
    redist_carbon_subsurface_output_g_c: f64 = 0,
    redist_oxygen_surface_input_g_o: f64 = 0,
    redist_oxygen_subsurface_output_g_o: f64 = 0,
    redist_hydrogen_surface_input_g_h: f64 = 0,
    redist_hydrogen_subsurface_output_g_h: f64 = 0,
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
    plant_root_organic_uptake_g_c: f64 = 0,
    plant_root_organic_exudate_g_c: f64 = 0,
    plant_root_organic_uptake_g_n: f64 = 0,
    plant_root_organic_exudate_g_n: f64 = 0,
    plant_root_organic_uptake_g_p: f64 = 0,
    plant_root_organic_exudate_g_p: f64 = 0,
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
/// HEAT-001 resolution A. Latent heat of fusion used to place boundary
/// snowfall on the same enthalpy reference state as the landscape census.
///
/// This is a temporary bridge. Zig has no default function parameters and the
/// only caller of `accumulateAcceptedPrecipitationHeat` lives in
/// `src/ecosys_ng.zig`, which is owned by the Integrator lane, so the value
/// cannot yet be threaded from
/// `runscript.snow_latent_heat_of_fusion_megajoules_per_m3`. It equals the
/// value every shipped runscript carries. See
/// `docs/binding_requests/heat_001_landscape_enthalpy.md`; once that binding
/// lands this constant must become a parameter so a runscript that changes
/// the constant cannot silently disagree with the census.
pub const snowfall_latent_heat_of_fusion_megajoules_per_m3: f64 = 333;

/// HEAT-001 third layer. Heat capacities used to place boundary snowfall on
/// the *ice branch* of the enthalpy curve, identically to
/// `landscape_mass_inventory.aggregateSnowEnthalpy`.
///
/// `2.095` is the solid-snow volumetric heat capacity, `redist.f:4072`
/// (`VHCPW = 2.095*VOLSSL + ...`), and it is what the snow owner publishes in
/// `snow_solute_transport` for solid water equivalent. `4.19` and `1.9274` are
/// the liquid-water and ice capacities from the same REDIST expression, and
/// they are the pair that defines the ice-branch reference state. They are
/// constants here for the same reason the latent heat above is: the sole
/// caller lives in `src/ecosys_ng.zig`, which lane A6 may not edit. See
/// `docs/binding_requests/A6_snowfall_boundary_ice_branch.md`.
pub const solid_snowfall_heat_capacity_megajoules_per_m3_k: f64 = 2.095;
pub const snowfall_liquid_water_heat_capacity_megajoules_per_m3_k: f64 = 4.19;
pub const snowfall_ice_heat_capacity_megajoules_per_m3_k: f64 = 1.9274;

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
        boundary_water_inflow_m3: []const f64,
        runoff_m3: []const f64,
        evaporation_m3: []const f64,
        water_outflow_m3: []const f64,
        artificial_drainage_outflow_m3: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            rainfall_m3,
            boundary_water_inflow_m3,
            runoff_m3,
            evaporation_m3,
            water_outflow_m3,
            artificial_drainage_outflow_m3,
        });
        try self.accumulateAccepted(.{
            .rain_m3 = try sumNonnegative(rainfall_m3),
            .boundary_water_inflow_m3 = try sumNonnegative(boundary_water_inflow_m3),
            .runoff_m3 = try sumNonnegative(runoff_m3),
            .evaporation_m3 = try sumNonnegative(evaporation_m3),
            // The accepted Richards external boundary residual already
            // contains artificial drainage. Keep its process diagnostic in
            // the daily ledger, but do not count it a second time here.
            .water_outflow_m3 = try sumNonnegative(water_outflow_m3),
        });
        _ = try sumNonnegative(artificial_drainage_outflow_m3);
    }

    /// Atmospheric sensible heat carried into the landscape by liquid
    /// precipitation (including irrigation already merged with rain) and
    /// snowfall water equivalent. This is evaluated before canopy/snow/
    /// litter/soil routing so intercepted water is included exactly once.
    ///
    /// HEAT-001 resolution A. The landscape census is an *enthalpy* on the
    /// reference state "liquid water at 0 K", so frozen water is carried one
    /// latent heat of fusion below liquid water at the same temperature.
    /// Snowfall crosses the boundary already frozen, so the enthalpy it
    /// actually delivers is `C_ice*T - L*swe`, not `C_ice*T`. Omitting the
    /// latent term made the ledger credit the landscape with energy the
    /// inventory never stored, and it is what made the day-one Ottawa
    /// deviation grow rather than close after the census became an enthalpy.
    ///
    /// This is *not* the forbidden Resolution B. No internal freeze/thaw
    /// conversion is booked here: this term is the reference-state offset of
    /// a boundary *mass* inflow that arrives in the solid phase, exactly as
    /// the 4.19 and 2.095 heat capacities below are the sensible part of the
    /// same inflow.
    ///
    /// HEAT-001 third layer, and the correction that closes the day-one ramp.
    /// The census reference state is liquid water at 0 K, on which a cubic
    /// metre of frozen water carries the *ice-branch* enthalpy
    /// `C_l*Tm - L + C_i*(T - Tm)`, not `C*T - L`. The two differ by
    /// `(C_l - C_i)*Tm = 617.5 MJ m-3`. This routine booked the superseded
    /// form while `landscape_mass_inventory.aggregateSnowEnthalpy` stored the
    /// corrected one, so every cubic metre of snowfall was credited
    /// `617.5 MJ` less than the snowpack then stored, and the audit read that
    /// shortfall as a positive deviation.
    ///
    /// Measured on Ottawa day one: `5.753063130900103e5 m3` of snowfall water
    /// equivalent arrives in hours 16--23, `617.5 * 5.753e5 / 8.717e7 =
    /// +4.0754 MJ m-2` against a recorded deviation of `+4.4684 MJ m-2`. The
    /// hourly residual ramp is confined to exactly those hours.
    ///
    /// The credited value is composed the same way the snow census composes
    /// its stored value: the owner's published solid-snow sensible product
    /// `C_solid*T`, plus the shared re-basing correction
    /// `frozenWaterEnthalpyPerM3(T) - C_i*T`. It calls the census's own
    /// exported definition so the two sides cannot drift.
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
        var heat_input_megajoules: f64 = 0;
        var diagnostic_snowfall_water_equivalent_m3: f64 = 0;
        var diagnostic_rainfall_m3: f64 = 0;
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
            // Liquid precipitation is pure sensible `C_l*T` on this reference
            // state, unchanged, and it is what the liquid carriers store.
            const cell_liquid_heat_megajoules = temperature_k * area_m2 *
                snowfall_liquid_water_heat_capacity_megajoules_per_m3_k *
                liquid_depth_m;
            // Snowfall arrives frozen, so it is credited the same ice-branch
            // enthalpy the snow census will store for it.
            const snowfall_enthalpy_megajoules_per_m3 =
                solid_snowfall_heat_capacity_megajoules_per_m3_k * temperature_k +
                (inventory.frozenWaterEnthalpyPerM3(
                    temperature_k,
                    snowfall_liquid_water_heat_capacity_megajoules_per_m3_k,
                    snowfall_ice_heat_capacity_megajoules_per_m3_k,
                    snowfall_latent_heat_of_fusion_megajoules_per_m3,
                    inventory.default_pure_water_melting_temperature_k,
                ) - snowfall_ice_heat_capacity_megajoules_per_m3_k * temperature_k);
            const cell_snowfall_heat_megajoules =
                snowfall_enthalpy_megajoules_per_m3 * snow_depth_m * area_m2;
            const cell_heat_megajoules =
                cell_liquid_heat_megajoules + cell_snowfall_heat_megajoules;
            if (!std.math.isFinite(cell_heat_megajoules))
                return error.LandscapeBoundaryLedgerOverflow;
            heat_input_megajoules = try addFinite(
                heat_input_megajoules,
                cell_heat_megajoules,
            );
            diagnostic_snowfall_water_equivalent_m3 += snow_depth_m * area_m2;
            diagnostic_rainfall_m3 += liquid_depth_m * area_m2;
        }
        std.log.debug("heat latent instrument: precipitation snowfall_swe_m3={e} rain_m3={e} cumulative_snow_swe_m3={e}", .{ diagnostic_snowfall_water_equivalent_m3, diagnostic_rainfall_m3, diagnostic_snowfall_water_equivalent_m3 });
        try self.accumulateAcceptedSignedHeat(heat_input_megajoules);
    }

    /// Source HEATH plus THFLXC boundary accounting. Ground flux components
    /// are positive into the surface and expressed per horizontal area.
    /// Canopy water-energy changes are already extensive signed MJ and retain
    /// EXTRACT's exact current-minus-previous-minus-retained-rain equation.
    pub fn accumulateAcceptedSurfaceAndCanopyHeat(
        self: *State,
        ground_net_radiation_megajoules_per_m2: []const f64,
        ground_sensible_heat_megajoules_per_m2: []const f64,
        ground_latent_heat_megajoules_per_m2: []const f64,
        ground_vapor_sensible_heat_megajoules_per_m2: []const f64,
        cell_area_m2: []const f64,
        canopy_water_energy_change_megajoules: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{
            ground_net_radiation_megajoules_per_m2,
            ground_sensible_heat_megajoules_per_m2,
            ground_latent_heat_megajoules_per_m2,
            ground_vapor_sensible_heat_megajoules_per_m2,
            cell_area_m2,
        });
        var signed_heat_into_landscape_megajoules: f64 = 0;
        for (
            ground_net_radiation_megajoules_per_m2,
            ground_sensible_heat_megajoules_per_m2,
            ground_latent_heat_megajoules_per_m2,
            ground_vapor_sensible_heat_megajoules_per_m2,
            cell_area_m2,
        ) |net_radiation, sensible, latent, vapor_sensible, area_m2| {
            inline for (.{ net_radiation, sensible, latent, vapor_sensible, area_m2 }) |value|
                if (!std.math.isFinite(value))
                    return error.NonFiniteLandscapeBoundaryFlux;
            if (area_m2 <= 0) return error.InvalidSurfaceHeatBoundaryArea;
            signed_heat_into_landscape_megajoules = try addFinite(
                signed_heat_into_landscape_megajoules,
                (net_radiation + sensible + latent + vapor_sensible) * area_m2,
            );
        }
        for (canopy_water_energy_change_megajoules) |change_megajoules| {
            if (!std.math.isFinite(change_megajoules))
                return error.NonFiniteLandscapeBoundaryFlux;
            signed_heat_into_landscape_megajoules =
                try addFinite(signed_heat_into_landscape_megajoules, change_megajoules);
        }
        try self.accumulateAccepted(.{
            .heat_input_megajoules = @max(0, signed_heat_into_landscape_megajoules),
            .heat_output_megajoules = @max(0, -signed_heat_into_landscape_megajoules),
        });
    }

    pub fn accumulateAcceptedSignedHeat(
        self: *State,
        signed_heat_into_landscape_megajoules: f64,
    ) !void {
        if (!std.math.isFinite(signed_heat_into_landscape_megajoules))
            return error.NonFiniteLandscapeBoundaryFlux;
        try self.accumulateAccepted(.{
            .heat_input_megajoules = @max(0, signed_heat_into_landscape_megajoules),
            .heat_output_megajoules = @max(0, -signed_heat_into_landscape_megajoules),
        });
    }

    /// REDIST XHFLFR plus the phase part of HFLXO.
    ///
    /// HEAT-001 resolution A. EXEC now inventories *enthalpy*, so the latent
    /// heat of fusion released or absorbed by a surface liquid/ice
    /// repartition is already carried by the surface storage census in
    /// `landscape_mass_inventory.aggregateSurfacePhysicalAndGas` and must not
    /// also be booked here. Booking it would count the same joules twice and,
    /// worse, would put an internal phase conversion into the landscape
    /// boundary stream, corrupting the tier work in `docs/validation.md`.
    ///
    /// What remains a genuine audit adjustment is the *sensible* consequence
    /// of the repartition: liquid water and ice have different volumetric
    /// heat capacities, so moving water between them at fixed temperature
    /// changes `C*T` with no corresponding energy flow. `phase_heat_flux`
    /// is still accepted so the caller contract and dimension checks are
    /// unchanged and the argument stays available for validation.
    /// Positive values increase the enthalpy inventory.
    pub fn accumulateAcceptedSurfacePhaseSensibleAdjustment(
        self: *State,
        phase_heat_megajoules_per_m2: []const f64,
        ice_water_equivalent_change_m3: []const f64,
        accepted_temperature_k: []const f64,
        cell_area_m2: []const f64,
        liquid_heat_capacity_megajoules_per_m3_k: f64,
        ice_heat_capacity_megajoules_per_m3_k: f64,
    ) !void {
        try requireSameNonzeroLength(.{
            phase_heat_megajoules_per_m2,
            ice_water_equivalent_change_m3,
            accepted_temperature_k,
            cell_area_m2,
        });
        inline for (.{ liquid_heat_capacity_megajoules_per_m3_k, ice_heat_capacity_megajoules_per_m3_k }) |value|
            if (!std.math.isFinite(value) or value <= 0)
                return error.InvalidSurfacePhaseHeatCapacity;
        var signed_adjustment_megajoules: f64 = 0;
        for (
            phase_heat_megajoules_per_m2,
            ice_water_equivalent_change_m3,
            accepted_temperature_k,
            cell_area_m2,
        ) |phase_megajoules_per_m2, ice_change_m3, temperature_k, area_m2| {
            inline for (.{ phase_megajoules_per_m2, ice_change_m3, temperature_k, area_m2 }) |value|
                if (!std.math.isFinite(value))
                    return error.NonFiniteLandscapeBoundaryFlux;
            if (temperature_k <= 0 or area_m2 <= 0)
                return error.InvalidSurfacePhaseAuditState;
            signed_adjustment_megajoules = try addFinite(
                signed_adjustment_megajoules,
                (ice_heat_capacity_megajoules_per_m3_k -
                    liquid_heat_capacity_megajoules_per_m3_k) *
                    ice_change_m3 * temperature_k,
            );
        }
        try self.accumulateAcceptedSignedHeat(signed_adjustment_megajoules);
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

    /// Positive source-signed organic and gas recharge through external
    /// subsurface boundaries. Precipitation and irrigation chemistry retain
    /// separate authoritative input ledgers so recharge is not counted twice.
    pub fn accumulateAcceptedCarbonTransportInputs(
        self: *State,
        organic_carbon_input_g_c: []const f64,
        inorganic_carbon_input_g_c: []const f64,
    ) !void {
        if (organic_carbon_input_g_c.len == 0 or organic_carbon_input_g_c.len != inorganic_carbon_input_g_c.len)
            return error.LandscapeBoundaryGridDimensionMismatch;
        try self.accumulateAccepted(.{
            .carbon_dioxide_input_g_c = try addFinite(
                try sumNonnegative(organic_carbon_input_g_c),
                try sumNonnegative(inorganic_carbon_input_g_c),
            ),
        });
    }

    pub fn accumulateAcceptedNitrogenTransportInputs(
        self: *State,
        dissolved_organic_nitrogen_g_n: []const f64,
        dissolved_inorganic_nitrogen_g_n: []const f64,
    ) !void {
        try self.accumulateAccepted(.{
            .nitrogen_input_g_n = try addFinite(
                try sumNonnegative(dissolved_organic_nitrogen_g_n),
                try sumNonnegative(dissolved_inorganic_nitrogen_g_n),
            ),
        });
    }

    pub fn accumulateAcceptedPhosphorusTransportInputs(
        self: *State,
        dissolved_organic_phosphorus_g_p: []const f64,
        dissolved_inorganic_phosphorus_g_p: []const f64,
    ) !void {
        try self.accumulateAccepted(.{
            .phosphorus_input_g_p = try addFinite(
                try sumNonnegative(dissolved_organic_phosphorus_g_p),
                try sumNonnegative(dissolved_inorganic_phosphorus_g_p),
            ),
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

    /// Root-soil dissolved organic exchange: carbon that crossed the plant/soil
    /// boundary via RDFOMC/RDFOMN/RDFOMP. Positive per-plant values indicate
    /// root uptake from soil (tracked organic decreases); negative values
    /// indicate root exudation to soil (tracked organic increases). Both
    /// directions are split into nonneg fields so the EXEC sign convention is
    /// applied only in `mass_balance_audit.balance`.
    pub fn accumulateAcceptedRootSoilOrganicExchange(
        self: *State,
        carbon_g_c: []const f64,
        nitrogen_g_n: []const f64,
        phosphorus_g_p: []const f64,
    ) !void {
        try requireSameNonzeroLength(.{ carbon_g_c, nitrogen_g_n, phosphorus_g_p });
        var uptake_c: f64 = 0;
        var exudate_c: f64 = 0;
        var uptake_n: f64 = 0;
        var exudate_n: f64 = 0;
        var uptake_p: f64 = 0;
        var exudate_p: f64 = 0;
        for (carbon_g_c, nitrogen_g_n, phosphorus_g_p) |c, n, p| {
            inline for (.{ c, n, p }) |v|
                if (!std.math.isFinite(v)) return error.NonFiniteLandscapeBoundaryFlux;
            if (c > 0) uptake_c = try addFinite(uptake_c, c) else if (c < 0) exudate_c = try addFinite(exudate_c, -c);
            if (n > 0) uptake_n = try addFinite(uptake_n, n) else if (n < 0) exudate_n = try addFinite(exudate_n, -n);
            if (p > 0) uptake_p = try addFinite(uptake_p, p) else if (p < 0) exudate_p = try addFinite(exudate_p, -p);
        }
        try self.accumulateAccepted(.{
            .plant_root_organic_uptake_g_c = uptake_c,
            .plant_root_organic_exudate_g_c = exudate_c,
            .plant_root_organic_uptake_g_n = uptake_n,
            .plant_root_organic_exudate_g_n = exudate_n,
            .plant_root_organic_uptake_g_p = uptake_p,
            .plant_root_organic_exudate_g_p = exudate_p,
        });
    }

    /// REDIST atmospheric gas boundary terms entering the EXEC balances.
    /// Soil/litter DAY flux is positive atmosphere -> ecosystem. EXTRACT
    /// canopy XCNET/XHNET/XONET is deliberately excluded: living plant
    /// biomass is absent from EXEC storage, and legacy CO2GIN/OXYGIN never
    /// receive canopy photosynthesis or respiration. Canopy exchange remains
    /// an ecosystem-flux diagnostic owned by `daily_canopy_gas_exchange`.
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
            canopy.net_oxygen_uptake_g_o.len != soil.cell_count or
            canopy.fire_carbon_dioxide_emission_g_c.len != soil.cell_count or
            canopy.fire_methane_emission_g_c.len != soil.cell_count or
            canopy.fire_oxygen_consumption_g_o.len != soil.cell_count)
            return error.LandscapeGasBoundaryDimensionMismatch;

        var carbon_input_g_c: f64 = 0;
        var carbon_output_g_c: f64 = 0;
        var oxygen_input_g: f64 = 0;
        var oxygen_output_g: f64 = 0;
        var gaseous_nitrogen_input_g_n: f64 = 0;
        var gaseous_nitrogen_output_g_n: f64 = 0;
        for (0..soil.cell_count) |cell| {
            // EXTRACT adds shoot/root fire products directly to CO2GIN and
            // OXYGIN. They share the canopy diagnostic accumulator but are
            // true external boundary fluxes, unlike photosynthesis.
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                -canopy.fire_carbon_dioxide_emission_g_c[cell],
            );
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                -canopy.fire_methane_emission_g_c[cell],
            );
            try splitSignedAtmosphereFlux(
                &oxygen_input_g,
                &oxygen_output_g,
                canopy.fire_oxygen_consumption_g_o[cell],
            );
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                try soil.get(cell, .carbon_dioxide),
            );
            try splitSignedAtmosphereFlux(
                &carbon_input_g_c,
                &carbon_output_g_c,
                try soil.get(cell, .methane),
            );
            try splitSignedAtmosphereFlux(
                &oxygen_input_g,
                &oxygen_output_g,
                try soil.get(cell, .oxygen),
            );
            inline for (.{
                gas.Species.nitrogen,
                gas.Species.nitrous_oxide,
                gas.Species.ammonia,
            }) |species|
                try splitSignedAtmosphereFlux(
                    &gaseous_nitrogen_input_g_n,
                    &gaseous_nitrogen_output_g_n,
                    try soil.get(cell, species),
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

    pub fn accumulateAcceptedRedistSurfaceGas(
        self: *State,
        redist: soil_daily_gas.RedistSurfaceGasIncrements,
    ) !void {
        inline for (std.meta.fields(soil_daily_gas.RedistSurfaceGasIncrements)) |field| {
            if (!std.math.isFinite(@field(redist, field.name)))
                return error.NonFiniteLandscapeBoundaryFlux;
        }
        var next = self.cumulative;
        next.redist_carbon_surface_input_g_c = try addFinite(
            next.redist_carbon_surface_input_g_c,
            redist.carbon_surface_input_g_c,
        );
        next.redist_carbon_subsurface_output_g_c = try addFinite(
            next.redist_carbon_subsurface_output_g_c,
            redist.carbon_subsurface_output_g_c,
        );
        next.redist_oxygen_surface_input_g_o = try addFinite(
            next.redist_oxygen_surface_input_g_o,
            redist.oxygen_surface_input_g_o,
        );
        next.redist_oxygen_subsurface_output_g_o = try addFinite(
            next.redist_oxygen_subsurface_output_g_o,
            redist.oxygen_subsurface_output_g_o,
        );
        next.redist_hydrogen_surface_input_g_h = try addFinite(
            next.redist_hydrogen_surface_input_g_h,
            redist.hydrogen_surface_input_g_h,
        );
        next.redist_hydrogen_subsurface_output_g_h = try addFinite(
            next.redist_hydrogen_subsurface_output_g_h,
            redist.hydrogen_subsurface_output_g_h,
        );
        self.cumulative = next;
    }

    /// Publishes boundary fields only. Area and reconstructed storage remain
    /// owned by the runtime landscape inventory.
    pub fn publish(self: State, totals: *audit.Totals) !void {
        try validate(self.cumulative);
        const f = self.cumulative;
        totals.cumulative_rain_m3 = try addFinite(f.rain_m3, f.boundary_water_inflow_m3);
        totals.cumulative_runoff_m3 = f.runoff_m3;
        totals.cumulative_evaporation_m3 = f.evaporation_m3;
        totals.cumulative_water_outflow_m3 = f.water_outflow_m3;
        totals.cumulative_heat_input_megajoules = f.heat_input_megajoules;
        totals.cumulative_heat_output_megajoules = f.heat_output_megajoules;
        totals.cumulative_oxygen_input_g = f.oxygen_input_g;
        totals.cumulative_oxygen_output_g = f.oxygen_output_g;
        totals.cumulative_redist_carbon_surface_input_g_c = f.redist_carbon_surface_input_g_c;
        totals.cumulative_redist_carbon_subsurface_output_g_c = f.redist_carbon_subsurface_output_g_c;
        totals.cumulative_redist_oxygen_surface_input_g_o = f.redist_oxygen_surface_input_g_o;
        totals.cumulative_redist_oxygen_subsurface_output_g_o = f.redist_oxygen_subsurface_output_g_o;
        totals.cumulative_redist_hydrogen_surface_input_g_h = f.redist_hydrogen_surface_input_g_h;
        totals.cumulative_redist_hydrogen_subsurface_output_g_h = f.redist_hydrogen_subsurface_output_g_h;
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
        totals.cumulative_plant_root_organic_carbon_uptake_g = f.plant_root_organic_uptake_g_c;
        totals.cumulative_plant_root_organic_carbon_exudate_g = f.plant_root_organic_exudate_g_c;
        totals.cumulative_plant_root_organic_nitrogen_uptake_g = f.plant_root_organic_uptake_g_n;
        totals.cumulative_plant_root_organic_nitrogen_exudate_g = f.plant_root_organic_exudate_g_n;
        totals.cumulative_plant_root_organic_phosphorus_uptake_g = f.plant_root_organic_uptake_g_p;
        totals.cumulative_plant_root_organic_phosphorus_exudate_g = f.plant_root_organic_exudate_g_p;
    }
};

fn validate(fluxes: Fluxes) !void {
    inline for (std.meta.fields(Fluxes)) |field| {
        const value = @field(fluxes, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteLandscapeBoundaryFlux;
        const is_redist_term = std.mem.eql(
            u8,
            field.name,
            "redist_carbon_surface_input_g_c",
        ) or std.mem.eql(u8, field.name, "redist_carbon_subsurface_output_g_c") or
            std.mem.eql(u8, field.name, "redist_oxygen_surface_input_g_o") or
            std.mem.eql(u8, field.name, "redist_oxygen_subsurface_output_g_o") or
            std.mem.eql(u8, field.name, "redist_hydrogen_surface_input_g_h") or
            std.mem.eql(u8, field.name, "redist_hydrogen_subsurface_output_g_h");
        if (!is_redist_term and value < 0) return error.NegativeLandscapeBoundaryFlux;
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
    try std.testing.expectEqual(@as(f64, 6), totals.cumulative_rain_m3);
    try std.testing.expectEqual(@as(f64, 2) * first.ion_output_mol, totals.cumulative_ion_output_mol);
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
        &second,
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
    try std.testing.expectEqual(@as(f64, 7), state.cumulative.boundary_water_inflow_m3);
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
    // HEAT-001 third layer. Snowfall arrives frozen, so the enthalpy it
    // delivers on the census reference state (liquid water at 0 K) is the
    // ice-branch value `C_l*Tm - L + C_i*(T - Tm)` re-based onto the snow
    // owner's published solid capacity, NOT `C*T - L`. Rain carries no latent
    // offset. Written here as the closed form rather than by calling the
    // production helper, so the test is an independent statement of the
    // reference state and not a restatement of the code under test.
    const latent = snowfall_latent_heat_of_fusion_megajoules_per_m3;
    const melting = inventory.default_pure_water_melting_temperature_k;
    const liquid_capacity = snowfall_liquid_water_heat_capacity_megajoules_per_m3_k;
    const ice_capacity = snowfall_ice_heat_capacity_megajoules_per_m3_k;
    const solid_capacity = solid_snowfall_heat_capacity_megajoules_per_m3_k;
    // Per cubic metre of snowfall water equivalent at temperature T.
    const snowfallEnthalpy = struct {
        fn at(temperature: f64) f64 {
            return solid_capacity * temperature +
                (liquid_capacity * melting - latent +
                    ice_capacity * (temperature - melting)) -
                ice_capacity * temperature;
        }
    }.at;
    const expected =
        280 * 10 * liquid_capacity * 0.001 + snowfallEnthalpy(280) * 0.003 * 10 +
        285 * 20 * liquid_capacity * 0.002 + snowfallEnthalpy(285) * 0.004 * 20;
    // The ice-branch offset `(C_l - C_i)*Tm = 617.5 MJ m-3` is nearly twice the
    // latent heat, so the corrected snowfall enthalpy is *positive* where the
    // superseded `C*T - L` form was negative at these depths. That sign change
    // is the whole day-one Ottawa correction, and pinning it here is what stops
    // the boundary from silently drifting away from the census again.
    try std.testing.expect(snowfallEnthalpy(280) > 0);
    try std.testing.expectApproxEqAbs(
        (liquid_capacity - ice_capacity) * melting,
        snowfallEnthalpy(280) -
            (solid_capacity * 280 - latent + ice_capacity * 280 - ice_capacity * 280),
        1.0e-9,
    );
    try std.testing.expect(expected > 0);
    try std.testing.expectApproxEqAbs(
        expected,
        state.cumulative.heat_input_megajoules,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), state.cumulative.heat_output_megajoules);
}

test "invalid precipitation heat transaction cannot advance ledger" {
    var state: State = .{};
    try state.accumulateAccepted(.{ .heat_input_megajoules = 7 });
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
    try std.testing.expectEqual(@as(f64, 0), state.cumulative.heat_input_megajoules);
    try std.testing.expectApproxEqAbs(
        @as(f64, 3),
        state.cumulative.heat_output_megajoules,
        1.0e-12,
    );
}

test "invalid late canopy heat leaves boundary ledger unchanged" {
    var state: State = .{};
    try state.accumulateAccepted(.{ .heat_input_megajoules = 4 });
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

test "surface phase audit books only the sensible capacity rebase, not latent" {
    var state: State = .{};
    try state.accumulateAcceptedSurfacePhaseSensibleAdjustment(
        &.{3},
        &.{0.25},
        &.{280},
        &.{2},
        4,
        2,
    );
    // HEAT-001 resolution A: the `3*2` latent release is now carried by the
    // surface enthalpy inventory and must not be booked here as well, or the
    // same joules are counted twice and an internal phase conversion enters
    // the boundary stream. Only the capacity rebase remains:
    // (2-4)*0.25*280 = -140 MJ.
    try std.testing.expectEqual(@as(f64, 0), state.cumulative.heat_input_megajoules);
    try std.testing.expectEqual(@as(f64, 140), state.cumulative.heat_output_megajoules);
}

test "EXEC atmospheric gases exclude canopy exchange and retain soil species signs" {
    var soil = try soil_daily_gas.State.init(std.testing.allocator, 1);
    defer soil.deinit();
    var canopy = try canopy_daily_gas.State.init(std.testing.allocator, 1);
    defer canopy.deinit();
    // Soil gas convention matches the gas transport solver: positive = atmosphere
    // entering the ecosystem, negative = ecosystem emitting to atmosphere.
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.carbon_dioxide)
    ] = -3; // CO2 emitted from soil
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.methane)
    ] = 2; // CH4 absorbed from atmosphere into soil
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.oxygen)
    ] = -4; // O2 emitted from soil
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.nitrogen)
    ] = 7; // N2 entering soil from atmosphere
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.nitrous_oxide)
    ] = -8; // N2O emitted from soil
    soil.tracked_element_mass_g_by_cell_and_species[
        @intFromEnum(gas.Species.ammonia)
    ] = 9; // NH3 absorbed from atmosphere into soil
    canopy.net_carbon_dioxide_uptake_g_c[0] = 5;
    canopy.net_methane_uptake_g_c[0] = -1;
    canopy.net_oxygen_uptake_g_o[0] = 6;
    canopy.fire_carbon_dioxide_emission_g_c[0] = 11;
    canopy.fire_methane_emission_g_c[0] = 13;
    canopy.fire_oxygen_consumption_g_o[0] = 17;

    var state: State = .{};
    try state.accumulateAcceptedAtmosphericGas(&soil, &canopy);
    try std.testing.expectEqual(
        @as(f64, 2),
        state.cumulative.carbon_dioxide_input_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 27),
        state.cumulative.carbon_output_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 17),
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
