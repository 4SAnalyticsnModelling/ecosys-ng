const std = @import("std");

pub const HorizontalSubsurfaceExchange = enum {
    enabled,
    standalone,
};

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const ExternalBoundaryWindow = struct {
    first_column: usize,
    last_column_inclusive: usize,
    first_row: usize,
    last_row_inclusive: usize,
};

pub const Face = enum(u8) {
    east,
    west,
    south,
    north,
    lower,

    pub const count: usize = @typeInfo(Face).@"enum".fields.len;
};

pub const OrganicFlux = struct {
    micropore_carbon_g_c_per_step: []const f64,
    micropore_acetate_carbon_g_c_per_step: []const f64,
    macropore_carbon_g_c_per_step: []const f64,
    macropore_acetate_carbon_g_c_per_step: []const f64,
    micropore_nitrogen_g_n_per_step: []const f64,
    macropore_nitrogen_g_n_per_step: []const f64,
    micropore_phosphorus_g_p_per_step: []const f64,
    macropore_phosphorus_g_p_per_step: []const f64,
};

pub const CarbonFlux = struct {
    micropore_carbon_dioxide_g_c_per_step: f64 = 0,
    macropore_carbon_dioxide_g_c_per_step: f64 = 0,
    gas_carbon_dioxide_g_c_per_step: f64 = 0,
    micropore_methane_g_c_per_step: f64 = 0,
    macropore_methane_g_c_per_step: f64 = 0,
    gas_methane_g_c_per_step: f64 = 0,
};

pub const NitrogenFlux = struct {
    micropore_ammonium_non_band_g_n_per_step: f64 = 0,
    micropore_ammonia_non_band_g_n_per_step: f64 = 0,
    micropore_nitrate_non_band_g_n_per_step: f64 = 0,
    micropore_ammonium_band_g_n_per_step: f64 = 0,
    micropore_ammonia_band_g_n_per_step: f64 = 0,
    micropore_nitrate_band_g_n_per_step: f64 = 0,
    micropore_nitrite_non_band_g_n_per_step: f64 = 0,
    micropore_nitrite_band_g_n_per_step: f64 = 0,
    macropore_ammonium_non_band_g_n_per_step: f64 = 0,
    macropore_ammonia_non_band_g_n_per_step: f64 = 0,
    macropore_nitrate_non_band_g_n_per_step: f64 = 0,
    macropore_ammonium_band_g_n_per_step: f64 = 0,
    macropore_ammonia_band_g_n_per_step: f64 = 0,
    macropore_nitrate_band_g_n_per_step: f64 = 0,
    macropore_nitrite_non_band_g_n_per_step: f64 = 0,
    macropore_nitrite_band_g_n_per_step: f64 = 0,
    micropore_dinitrogen_g_n_per_step: f64 = 0,
    gas_dinitrogen_g_n_per_step: f64 = 0,
    macropore_dinitrogen_g_n_per_step: f64 = 0,
    micropore_nitrous_oxide_g_n_per_step: f64 = 0,
    gas_nitrous_oxide_g_n_per_step: f64 = 0,
    macropore_nitrous_oxide_g_n_per_step: f64 = 0,
    gas_ammonia_g_n_per_step: f64 = 0,
};

pub const PhosphorusElementFlux = struct {
    micropore_h2po4_non_band_g_p_per_step: f64 = 0,
    micropore_h2po4_band_g_p_per_step: f64 = 0,
    macropore_h2po4_non_band_g_p_per_step: f64 = 0,
    macropore_h2po4_band_g_p_per_step: f64 = 0,
    micropore_hpo4_non_band_g_p_per_step: f64 = 0,
    micropore_hpo4_band_g_p_per_step: f64 = 0,
    macropore_hpo4_non_band_g_p_per_step: f64 = 0,
    macropore_hpo4_band_g_p_per_step: f64 = 0,
};

pub const GasFlux = struct {
    micropore_oxygen_g_o_per_step: f64 = 0,
    macropore_oxygen_g_o_per_step: f64 = 0,
    gas_oxygen_g_o_per_step: f64 = 0,
    micropore_hydrogen_g_h_per_step: f64 = 0,
    macropore_hydrogen_g_h_per_step: f64 = 0,
    gas_hydrogen_g_h_per_step: f64 = 0,
};

pub const PrimaryIonFlux = struct {
    aluminum_mol_per_step: f64 = 0,
    iron_mol_per_step: f64 = 0,
    hydrogen_mol_per_step: f64 = 0,
    calcium_mol_per_step: f64 = 0,
    magnesium_mol_per_step: f64 = 0,
    sodium_mol_per_step: f64 = 0,
    potassium_mol_per_step: f64 = 0,
    hydroxide_mol_per_step: f64 = 0,
    sulfate_mol_per_step: f64 = 0,
    chloride_mol_per_step: f64 = 0,
    carbonate_mol_per_step: f64 = 0,
};

pub const SecondaryIonFlux = struct {
    bicarbonate_mol_per_step: f64 = 0,
    aluminum_hydroxide_1_mol_per_step: f64 = 0,
    aluminum_sulfate_mol_per_step: f64 = 0,
    iron_hydroxide_1_mol_per_step: f64 = 0,
    iron_sulfate_mol_per_step: f64 = 0,
    calcium_hydroxide_mol_per_step: f64 = 0,
    calcium_carbonate_mol_per_step: f64 = 0,
    calcium_sulfate_mol_per_step: f64 = 0,
    magnesium_hydroxide_mol_per_step: f64 = 0,
    magnesium_carbonate_mol_per_step: f64 = 0,
    magnesium_sulfate_mol_per_step: f64 = 0,
    sodium_carbonate_mol_per_step: f64 = 0,
    sodium_sulfate_mol_per_step: f64 = 0,
    potassium_sulfate_mol_per_step: f64 = 0,
};

pub const TertiaryIonFlux = struct {
    aluminum_hydroxide_2_mol_per_step: f64 = 0,
    iron_hydroxide_2_mol_per_step: f64 = 0,
    calcium_bicarbonate_mol_per_step: f64 = 0,
    magnesium_bicarbonate_mol_per_step: f64 = 0,
};

pub const QuaternaryIonFlux = struct {
    aluminum_hydroxide_3_mol_per_step: f64 = 0,
    iron_hydroxide_3_mol_per_step: f64 = 0,
};

pub const QuinaryIonFlux = struct {
    aluminum_hydroxide_4_mol_per_step: f64 = 0,
    iron_hydroxide_4_mol_per_step: f64 = 0,
};

pub const PhosphateComplexFlux = struct {
    po4_mol_p_per_step: f64 = 0,
    calcium_po4_mol_p_per_step: f64 = 0,
    iron_hpo4_mol_p_per_step: f64 = 0,
    calcium_hpo4_mol_p_per_step: f64 = 0,
    magnesium_hpo4_mol_p_per_step: f64 = 0,
    phosphoric_acid_mol_p_per_step: f64 = 0,
    iron_h2po4_mol_p_per_step: f64 = 0,
    calcium_h2po4_mol_p_per_step: f64 = 0,
};

pub const PoreSaltFlux = struct {
    primary: PrimaryIonFlux = .{},
    secondary: SecondaryIonFlux = .{},
    tertiary: TertiaryIonFlux = .{},
    quaternary: QuaternaryIonFlux = .{},
    quinary: QuinaryIonFlux = .{},
    phosphate_non_band: PhosphateComplexFlux = .{},
    phosphate_band: PhosphateComplexFlux = .{},
};

pub const BoundaryFlux = struct {
    micropore_water_m3_per_step: f64 = 0,
    macropore_water_m3_per_step: f64 = 0,
    organic: OrganicFlux,
    carbon: CarbonFlux = .{},
    nitrogen: NitrogenFlux = .{},
    phosphorus: PhosphorusElementFlux = .{},
    gases: GasFlux = .{},
    micropore_salts: PoreSaltFlux = .{},
    macropore_salts: PoreSaltFlux = .{},
    micropore_hydrogen_exchange_mol_per_step: f64 = 0,
};

pub const ConductivityCoefficients = struct {
    hydrogen_dS_m2_per_mol: f64,
    hydroxide_dS_m2_per_mol: f64,
    aluminum_dS_m2_per_equivalent: f64,
    iron_dS_m2_per_equivalent: f64,
    calcium_dS_m2_per_equivalent: f64,
    magnesium_dS_m2_per_equivalent: f64,
    sodium_dS_m2_per_mol: f64,
    potassium_dS_m2_per_mol: f64,
    carbonate_dS_m2_per_equivalent: f64,
    bicarbonate_dS_m2_per_mol: f64,
    sulfate_dS_m2_per_equivalent: f64,
    chloride_dS_m2_per_mol: f64,
    nitrate_dS_m2_per_mol: f64,
};

pub const CumulativeOutwardLedger = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
    oxygen_g_o: f64 = 0,
    hydrogen_g_h: f64 = 0,
    ion_components_mol: f64 = 0,
};

pub const CellOutwardLedger = struct {
    organic_carbon_g_c: f64 = 0,
    inorganic_carbon_g_c: f64 = 0,
    organic_nitrogen_g_n: f64 = 0,
    inorganic_nitrogen_g_n: f64 = 0,
    organic_phosphorus_g_p: f64 = 0,
    inorganic_phosphorus_g_p: f64 = 0,
    ion_components_mol: f64 = 0,
};

pub const Inputs = struct {
    lon_count: usize,
    lat_count: usize,
    soil_layer_count: usize,
    organic_fraction_count: usize,
    external_boundary_window: ExternalBoundaryWindow,
    horizontal_exchange_by_cell: []const HorizontalSubsurfaceExchange,
    salt_equilibrium_mode: SaltEquilibriumMode,
    flux_by_face: [Face.count][]const BoundaryFlux,
    negligible_water_m3_by_cell: []const f64,
    nitrogen_g_n_per_mol_n: f64,
    phosphorus_g_p_per_mol_p: f64,
    conductivity: ConductivityCoefficients,
};

pub const State = struct {
    cumulative: CumulativeOutwardLedger,
    outward_by_cell: []CellOutwardLedger,
    last_electrical_conductivity_dS_per_m_by_cell: []f64,
    /// Source `SG`: intentionally unsigned and excludes macropore H2.
    source_hydrogen_diagnostic_g_h: f64,
};

/// Accounts for subsurface C, N, P, O2, H2, salts, and conductivity.
///
/// Traceability: REDIST.F lines 1191--1432 (`COD`, `ZOD`, `POD`, `CXD`,
/// `ZXD`, `ZGD`, `PXD`, `OOD`, `HOD`, `PQD`, `PHD`, `SSD`, `SHD`, `SO`,
/// `ECNDX`, and `SG`). The enclosing exact nonzero ordinary-water gate and
/// source column/row/layer/face order are retained.
///
/// Two symmetric source defects are corrected: line 1345 repeats macropore
/// AlOH4 and omits FeOH4; lines 1399--1400 use macropore Ca in the aluminum
/// conductivity term. The corresponding micropore equations and neighboring
/// ion equations establish the intended macropore FeOH4 and Al terms.
pub fn apply(inputs: Inputs, state: *State) !void {
    const dimensions = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, dimensions);
    try preflightUpdates(inputs, state.*);

    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            for (0..inputs.soil_layer_count) |layer| {
                const record = cell * inputs.soil_layer_count + layer;
                for (std.meta.tags(Face)) |face| {
                    if (!isActiveBoundary(inputs, face, column, row, cell, layer))
                        continue;
                    const flux = inputs.flux_by_face[@intFromEnum(face)][record];
                    const water_m3 = tryWater(flux) catch unreachable;
                    if (water_m3 == 0) continue;
                    const transfer = calculateTransfer(
                        flux,
                        inwardSign(face),
                        inputs,
                    ) catch unreachable;
                    state.cumulative =
                        advanceCumulative(state.cumulative, transfer) catch unreachable;
                    state.outward_by_cell[cell] =
                        advanceCell(state.outward_by_cell[cell], transfer) catch unreachable;
                    if (inputs.salt_equilibrium_mode == .dynamic)
                        state.last_electrical_conductivity_dS_per_m_by_cell[cell] =
                            electricalConductivityForCell(
                                flux,
                                water_m3,
                                inputs,
                                cell,
                            ) catch unreachable;
                    state.source_hydrogen_diagnostic_g_h = checkedSum(
                        state.source_hydrogen_diagnostic_g_h,
                        transfer.source_hydrogen_diagnostic_increment_g_h,
                    ) catch unreachable;
                }
            }
        }
    }
}

const SignedTransfer = struct {
    organic_carbon_g_c: f64,
    inorganic_carbon_g_c: f64,
    organic_nitrogen_g_n: f64,
    inorganic_nitrogen_g_n: f64,
    gaseous_nitrogen_g_n: f64,
    organic_phosphorus_g_p: f64,
    inorganic_phosphorus_g_p: f64,
    salt_phosphorus_g_p: f64,
    oxygen_g_o: f64,
    hydrogen_g_h: f64,
    ion_components_mol: f64,
    source_hydrogen_diagnostic_increment_g_h: f64,
};

fn calculateTransfer(
    flux: BoundaryFlux,
    direction: f64,
    inputs: Inputs,
) !SignedTransfer {
    const organic = try organicElements(
        flux.organic,
        inputs.organic_fraction_count,
        direction,
    );
    const inorganic_carbon_g_c =
        try signedStructSum(direction, flux.carbon);
    const inorganic_nitrogen_g_n =
        try mineralNitrogen(flux.nitrogen, direction);
    const gaseous_nitrogen_g_n =
        try gaseousNitrogen(flux.nitrogen, direction);
    const inorganic_phosphorus_g_p =
        try signedStructSum(direction, flux.phosphorus);
    const oxygen_g_o = try signedSourceSum(direction, &.{
        flux.gases.micropore_oxygen_g_o_per_step,
        flux.gases.macropore_oxygen_g_o_per_step,
        flux.gases.gas_oxygen_g_o_per_step,
    });
    const hydrogen_g_h = try signedSourceSum(direction, &.{
        flux.gases.micropore_hydrogen_g_h_per_step,
        flux.gases.macropore_hydrogen_g_h_per_step,
        flux.gases.gas_hydrogen_g_h_per_step,
    });
    var salt_phosphorus_g_p: f64 = 0;
    var ion_components_mol: f64 = 0;
    if (inputs.salt_equilibrium_mode == .dynamic) {
        salt_phosphorus_g_p = try saltPhosphorus(
            flux,
            direction,
            inputs.phosphorus_g_p_per_mol_p,
        );
        ion_components_mol = try saltIonComponents(flux, direction);
    }
    return .{
        .organic_carbon_g_c = organic.carbon_g_c,
        .inorganic_carbon_g_c = inorganic_carbon_g_c,
        .organic_nitrogen_g_n = organic.nitrogen_g_n,
        .inorganic_nitrogen_g_n = inorganic_nitrogen_g_n,
        .gaseous_nitrogen_g_n = gaseous_nitrogen_g_n,
        .organic_phosphorus_g_p = organic.phosphorus_g_p,
        .inorganic_phosphorus_g_p = inorganic_phosphorus_g_p,
        .salt_phosphorus_g_p = salt_phosphorus_g_p,
        .oxygen_g_o = oxygen_g_o,
        .hydrogen_g_h = hydrogen_g_h,
        .ion_components_mol = ion_components_mol,
        .source_hydrogen_diagnostic_increment_g_h = try checkedSum(
            flux.gases.micropore_hydrogen_g_h_per_step,
            flux.gases.gas_hydrogen_g_h_per_step,
        ),
    };
}

const OrganicElements = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

fn organicElements(
    flux: OrganicFlux,
    fraction_count: usize,
    direction: f64,
) !OrganicElements {
    var result: OrganicElements = .{
        .carbon_g_c = 0,
        .nitrogen_g_n = 0,
        .phosphorus_g_p = 0,
    };
    for (0..fraction_count) |fraction| {
        result.carbon_g_c = try checkedSum(
            result.carbon_g_c,
            try signedSourceSum(direction, &.{
                flux.micropore_carbon_g_c_per_step[fraction],
                flux.micropore_acetate_carbon_g_c_per_step[fraction],
                flux.macropore_carbon_g_c_per_step[fraction],
                flux.macropore_acetate_carbon_g_c_per_step[fraction],
            }),
        );
        result.nitrogen_g_n = try checkedSum(
            result.nitrogen_g_n,
            try signedSourceSum(direction, &.{
                flux.micropore_nitrogen_g_n_per_step[fraction],
                flux.macropore_nitrogen_g_n_per_step[fraction],
            }),
        );
        result.phosphorus_g_p = try checkedSum(
            result.phosphorus_g_p,
            try signedSourceSum(direction, &.{
                flux.micropore_phosphorus_g_p_per_step[fraction],
                flux.macropore_phosphorus_g_p_per_step[fraction],
            }),
        );
    }
    return result;
}

fn mineralNitrogen(flux: NitrogenFlux, direction: f64) !f64 {
    return signedSourceSum(direction, &.{
        flux.micropore_ammonium_non_band_g_n_per_step,
        flux.micropore_ammonia_non_band_g_n_per_step,
        flux.micropore_nitrate_non_band_g_n_per_step,
        flux.micropore_ammonium_band_g_n_per_step,
        flux.micropore_ammonia_band_g_n_per_step,
        flux.micropore_nitrate_band_g_n_per_step,
        flux.micropore_nitrite_non_band_g_n_per_step,
        flux.micropore_nitrite_band_g_n_per_step,
        flux.macropore_ammonium_non_band_g_n_per_step,
        flux.macropore_ammonia_non_band_g_n_per_step,
        flux.macropore_nitrate_non_band_g_n_per_step,
        flux.macropore_ammonium_band_g_n_per_step,
        flux.macropore_ammonia_band_g_n_per_step,
        flux.macropore_nitrate_band_g_n_per_step,
        flux.macropore_nitrite_non_band_g_n_per_step,
        flux.macropore_nitrite_band_g_n_per_step,
    });
}

fn gaseousNitrogen(flux: NitrogenFlux, direction: f64) !f64 {
    return signedSourceSum(direction, &.{
        flux.micropore_dinitrogen_g_n_per_step,
        flux.gas_dinitrogen_g_n_per_step,
        flux.macropore_dinitrogen_g_n_per_step,
        flux.micropore_nitrous_oxide_g_n_per_step,
        flux.gas_nitrous_oxide_g_n_per_step,
        flux.macropore_nitrous_oxide_g_n_per_step,
        flux.gas_ammonia_g_n_per_step,
    });
}

fn saltPhosphorus(
    flux: BoundaryFlux,
    direction: f64,
    phosphorus_g_p_per_mol_p: f64,
) !f64 {
    const micropore_g_p = try checkedProduct(
        try checkedProduct(direction, phosphorus_g_p_per_mol_p),
        try porePhosphorusMoles(flux.micropore_salts),
    );
    const macropore_g_p = try checkedProduct(
        try checkedProduct(direction, phosphorus_g_p_per_mol_p),
        try porePhosphorusMoles(flux.macropore_salts),
    );
    return checkedSum(micropore_g_p, macropore_g_p);
}

fn porePhosphorusMoles(pore: PoreSaltFlux) !f64 {
    const non_band = pore.phosphate_non_band;
    const band = pore.phosphate_band;
    return sourceSum(&.{
        non_band.po4_mol_p_per_step,
        band.po4_mol_p_per_step,
        non_band.calcium_po4_mol_p_per_step,
        band.calcium_po4_mol_p_per_step,
        non_band.iron_hpo4_mol_p_per_step,
        non_band.calcium_hpo4_mol_p_per_step,
        non_band.magnesium_hpo4_mol_p_per_step,
        band.iron_hpo4_mol_p_per_step,
        band.calcium_hpo4_mol_p_per_step,
        band.magnesium_hpo4_mol_p_per_step,
        non_band.phosphoric_acid_mol_p_per_step,
        non_band.iron_h2po4_mol_p_per_step,
        non_band.calcium_h2po4_mol_p_per_step,
        band.phosphoric_acid_mol_p_per_step,
        band.iron_h2po4_mol_p_per_step,
        band.calcium_h2po4_mol_p_per_step,
    });
}

fn saltIonComponents(flux: BoundaryFlux, direction: f64) !f64 {
    const micropore = try poreIonComponents(
        flux.micropore_salts,
        flux.micropore_hydrogen_exchange_mol_per_step,
    );
    const macropore = try poreIonComponents(flux.macropore_salts, 0);
    return checkedProduct(
        direction,
        try checkedSum(micropore, macropore),
    );
}

fn poreIonComponents(
    pore: PoreSaltFlux,
    quaternary_extra_mol_per_step: f64,
) !f64 {
    const non_band = pore.phosphate_non_band;
    const band = pore.phosphate_band;
    var total = try checkedSum(
        try sourceStructSum(pore.primary),
        try sourceSum(&.{
            non_band.po4_mol_p_per_step,
            band.po4_mol_p_per_step,
        }),
    );
    total = try checkedSum(
        total,
        try checkedProduct(
            2,
            try checkedSum(
                try sourceStructSum(pore.secondary),
                try sourceSum(&.{
                    non_band.calcium_po4_mol_p_per_step,
                    band.calcium_po4_mol_p_per_step,
                }),
            ),
        ),
    );
    total = try checkedSum(
        total,
        try checkedProduct(
            3,
            try checkedSum(
                try sourceStructSum(pore.tertiary),
                try sourceSum(&.{
                    non_band.iron_hpo4_mol_p_per_step,
                    non_band.calcium_hpo4_mol_p_per_step,
                    non_band.magnesium_hpo4_mol_p_per_step,
                    band.iron_hpo4_mol_p_per_step,
                    band.calcium_hpo4_mol_p_per_step,
                    band.magnesium_hpo4_mol_p_per_step,
                }),
            ),
        ),
    );
    var quaternary = try checkedSum(
        try sourceStructSum(pore.quaternary),
        try sourceSum(&.{
            non_band.phosphoric_acid_mol_p_per_step,
            non_band.iron_h2po4_mol_p_per_step,
            non_band.calcium_h2po4_mol_p_per_step,
            band.phosphoric_acid_mol_p_per_step,
            band.iron_h2po4_mol_p_per_step,
            band.calcium_h2po4_mol_p_per_step,
        }),
    );
    quaternary = try checkedSum(quaternary, quaternary_extra_mol_per_step);
    total = try checkedSum(
        total,
        try checkedProduct(4, quaternary),
    );
    total = try checkedSum(
        total,
        try checkedProduct(5, try sourceStructSum(pore.quinary)),
    );
    return total;
}

fn electricalConductivityAboveThreshold(
    flux: BoundaryFlux,
    water_m3: f64,
    inputs: Inputs,
) !f64 {
    const micro = flux.micropore_salts;
    const macro = flux.macropore_salts;
    const coefficient = inputs.conductivity;
    var total: f64 = 0;
    total = try addConductivity(total, coefficient.hydrogen_dS_m2_per_mol, micro.primary.hydrogen_mol_per_step, macro.primary.hydrogen_mol_per_step, 1, water_m3);
    total = try addConductivity(total, coefficient.hydroxide_dS_m2_per_mol, micro.primary.hydroxide_mol_per_step, macro.primary.hydroxide_mol_per_step, 1, water_m3);
    total = try addConductivity(total, coefficient.aluminum_dS_m2_per_equivalent, micro.primary.aluminum_mol_per_step, macro.primary.aluminum_mol_per_step, 3, water_m3);
    total = try addConductivity(total, coefficient.iron_dS_m2_per_equivalent, micro.primary.iron_mol_per_step, macro.primary.iron_mol_per_step, 3, water_m3);
    total = try addConductivity(total, coefficient.calcium_dS_m2_per_equivalent, micro.primary.calcium_mol_per_step, macro.primary.calcium_mol_per_step, 2, water_m3);
    total = try addConductivity(total, coefficient.magnesium_dS_m2_per_equivalent, micro.primary.magnesium_mol_per_step, macro.primary.magnesium_mol_per_step, 2, water_m3);
    total = try addConductivity(total, coefficient.sodium_dS_m2_per_mol, micro.primary.sodium_mol_per_step, macro.primary.sodium_mol_per_step, 1, water_m3);
    total = try addConductivity(total, coefficient.potassium_dS_m2_per_mol, micro.primary.potassium_mol_per_step, macro.primary.potassium_mol_per_step, 1, water_m3);
    total = try addConductivity(total, coefficient.carbonate_dS_m2_per_equivalent, micro.primary.carbonate_mol_per_step, macro.primary.carbonate_mol_per_step, 2, water_m3);
    total = try addConductivity(total, coefficient.bicarbonate_dS_m2_per_mol, micro.secondary.bicarbonate_mol_per_step, macro.secondary.bicarbonate_mol_per_step, 1, water_m3);
    total = try addConductivity(total, coefficient.sulfate_dS_m2_per_equivalent, micro.primary.sulfate_mol_per_step, macro.primary.sulfate_mol_per_step, 2, water_m3);
    total = try addConductivity(total, coefficient.chloride_dS_m2_per_mol, micro.primary.chloride_mol_per_step, macro.primary.chloride_mol_per_step, 1, water_m3);
    const nitrate_g_n = try checkedSum(
        flux.nitrogen.micropore_nitrate_non_band_g_n_per_step,
        flux.nitrogen.macropore_nitrate_non_band_g_n_per_step,
    );
    const nitrate_denominator = try checkedProduct(
        water_m3,
        inputs.nitrogen_g_n_per_mol_n,
    );
    const nitrate_mol_per_m3 = nitrate_g_n / nitrate_denominator;
    if (!std.math.isFinite(nitrate_mol_per_m3))
        return error.NonFiniteSubsurfaceSoluteResult;
    total = try checkedSum(
        total,
        try checkedProduct(
            coefficient.nitrate_dS_m2_per_mol,
            @max(0, nitrate_mol_per_m3),
        ),
    );
    return total;
}

fn addConductivity(
    current: f64,
    coefficient: f64,
    micropore_mol: f64,
    macropore_mol: f64,
    valence: f64,
    water_m3: f64,
) !f64 {
    const pore_sum = try checkedSum(micropore_mol, macropore_mol);
    const concentration =
        try checkedProduct(pore_sum, valence) / water_m3;
    if (!std.math.isFinite(concentration))
        return error.NonFiniteSubsurfaceSoluteResult;
    return checkedSum(
        current,
        try checkedProduct(coefficient, @max(0, concentration)),
    );
}

fn advanceCumulative(
    current: CumulativeOutwardLedger,
    transfer: SignedTransfer,
) !CumulativeOutwardLedger {
    var next = current;
    next.carbon_g_c =
        try checkedSum(next.carbon_g_c, -transfer.organic_carbon_g_c);
    next.carbon_g_c =
        try checkedSum(next.carbon_g_c, -transfer.inorganic_carbon_g_c);
    next.nitrogen_g_n =
        try checkedSum(next.nitrogen_g_n, -transfer.organic_nitrogen_g_n);
    next.nitrogen_g_n =
        try checkedSum(next.nitrogen_g_n, -transfer.inorganic_nitrogen_g_n);
    next.nitrogen_g_n =
        try checkedSum(next.nitrogen_g_n, -transfer.gaseous_nitrogen_g_n);
    next.phosphorus_g_p =
        try checkedSum(next.phosphorus_g_p, -transfer.organic_phosphorus_g_p);
    next.phosphorus_g_p =
        try checkedSum(next.phosphorus_g_p, -transfer.inorganic_phosphorus_g_p);
    next.phosphorus_g_p =
        try checkedSum(next.phosphorus_g_p, -transfer.salt_phosphorus_g_p);
    next.oxygen_g_o =
        try checkedSum(next.oxygen_g_o, -transfer.oxygen_g_o);
    next.hydrogen_g_h =
        try checkedSum(next.hydrogen_g_h, -transfer.hydrogen_g_h);
    next.ion_components_mol =
        try checkedSum(next.ion_components_mol, -transfer.ion_components_mol);
    return next;
}

fn advanceCell(
    current: CellOutwardLedger,
    transfer: SignedTransfer,
) !CellOutwardLedger {
    var next = current;
    next.organic_carbon_g_c =
        try checkedSum(next.organic_carbon_g_c, -transfer.organic_carbon_g_c);
    next.inorganic_carbon_g_c =
        try checkedSum(next.inorganic_carbon_g_c, -transfer.inorganic_carbon_g_c);
    next.organic_nitrogen_g_n =
        try checkedSum(next.organic_nitrogen_g_n, -transfer.organic_nitrogen_g_n);
    next.inorganic_nitrogen_g_n = try checkedSum(
        next.inorganic_nitrogen_g_n,
        -transfer.inorganic_nitrogen_g_n,
    );
    next.organic_phosphorus_g_p = try checkedSum(
        next.organic_phosphorus_g_p,
        -transfer.organic_phosphorus_g_p,
    );
    next.inorganic_phosphorus_g_p = try checkedSum(
        next.inorganic_phosphorus_g_p,
        -transfer.inorganic_phosphorus_g_p,
    );
    next.ion_components_mol =
        try checkedSum(next.ion_components_mol, -transfer.ion_components_mol);
    return next;
}

fn preflightUpdates(inputs: Inputs, state: State) !void {
    var cumulative = state.cumulative;
    var source_hydrogen_diagnostic_g_h =
        state.source_hydrogen_diagnostic_g_h;
    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            var cell_ledger = state.outward_by_cell[cell];
            for (0..inputs.soil_layer_count) |layer| {
                const record = cell * inputs.soil_layer_count + layer;
                for (std.meta.tags(Face)) |face| {
                    if (!isActiveBoundary(inputs, face, column, row, cell, layer))
                        continue;
                    const flux = inputs.flux_by_face[@intFromEnum(face)][record];
                    const water_m3 = try tryWater(flux);
                    if (water_m3 == 0) continue;
                    const transfer = try calculateTransfer(
                        flux,
                        inwardSign(face),
                        inputs,
                    );
                    cumulative = try advanceCumulative(cumulative, transfer);
                    cell_ledger = try advanceCell(cell_ledger, transfer);
                    if (inputs.salt_equilibrium_mode == .dynamic)
                        _ = try electricalConductivityForCell(
                            flux,
                            water_m3,
                            inputs,
                            cell,
                        );
                    source_hydrogen_diagnostic_g_h = try checkedSum(
                        source_hydrogen_diagnostic_g_h,
                        transfer.source_hydrogen_diagnostic_increment_g_h,
                    );
                }
            }
        }
    }
}

fn electricalConductivityForCell(
    flux: BoundaryFlux,
    water_m3: f64,
    inputs: Inputs,
    cell: usize,
) !f64 {
    if (@abs(water_m3) <= inputs.negligible_water_m3_by_cell[cell])
        return 0;
    return electricalConductivityAboveThreshold(flux, water_m3, inputs);
}

fn tryWater(flux: BoundaryFlux) !f64 {
    return checkedSum(
        flux.micropore_water_m3_per_step,
        flux.macropore_water_m3_per_step,
    );
}

const Dimensions = struct {
    cell_count: usize,
    record_count: usize,
};

fn validateDimensions(inputs: Inputs, state: State) !Dimensions {
    if (inputs.lon_count == 0 or
        inputs.lat_count == 0 or
        inputs.soil_layer_count == 0)
        return error.InvalidSubsurfaceSoluteDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.lon_count,
        inputs.lat_count,
    ) catch return error.InvalidSubsurfaceSoluteDimensions;
    const record_count = std.math.mul(
        usize,
        cell_count,
        inputs.soil_layer_count,
    ) catch return error.InvalidSubsurfaceSoluteDimensions;
    inline for (inputs.flux_by_face) |values|
        if (values.len != record_count)
            return error.InvalidSubsurfaceSoluteDimensions;
    if (inputs.horizontal_exchange_by_cell.len != cell_count or
        inputs.negligible_water_m3_by_cell.len != cell_count or
        state.outward_by_cell.len != cell_count or
        state.last_electrical_conductivity_dS_per_m_by_cell.len != cell_count)
        return error.InvalidSubsurfaceSoluteDimensions;
    const window = inputs.external_boundary_window;
    if (window.first_column > window.last_column_inclusive or
        window.first_row > window.last_row_inclusive or
        window.last_column_inclusive >= inputs.lon_count or
        window.last_row_inclusive >= inputs.lat_count)
        return error.InvalidSubsurfaceSoluteWindow;
    return .{ .cell_count = cell_count, .record_count = record_count };
}

fn validateInputsAndState(
    inputs: Inputs,
    state: State,
    dimensions: Dimensions,
) !void {
    if (!positiveFinite(inputs.nitrogen_g_n_per_mol_n) or
        !positiveFinite(inputs.phosphorus_g_p_per_mol_p))
        return error.InvalidSubsurfaceSoluteInput;
    try validatePositiveStruct(inputs.conductivity);
    inline for (inputs.flux_by_face) |records| {
        for (records) |flux| {
            inline for (.{
                flux.micropore_water_m3_per_step,
                flux.macropore_water_m3_per_step,
                flux.micropore_hydrogen_exchange_mol_per_step,
            }) |value| if (!std.math.isFinite(value))
                return error.InvalidSubsurfaceSoluteInput;
            try validateOrganic(flux.organic, inputs.organic_fraction_count);
            try validateFiniteStruct(flux.carbon);
            try validateFiniteStruct(flux.nitrogen);
            try validateFiniteStruct(flux.phosphorus);
            try validateFiniteStruct(flux.gases);
            try validatePoreSalt(flux.micropore_salts);
            try validatePoreSalt(flux.macropore_salts);
        }
    }
    for (0..dimensions.cell_count) |cell| {
        if (!nonnegativeFinite(inputs.negligible_water_m3_by_cell[cell]))
            return error.InvalidSubsurfaceSoluteInput;
        validateFiniteStruct(state.outward_by_cell[cell]) catch
            return error.InvalidSubsurfaceSoluteState;
        if (!nonnegativeFinite(
            state.last_electrical_conductivity_dS_per_m_by_cell[cell],
        ))
            return error.InvalidSubsurfaceSoluteState;
    }
    validateFiniteStruct(state.cumulative) catch
        return error.InvalidSubsurfaceSoluteState;
    if (!std.math.isFinite(state.source_hydrogen_diagnostic_g_h))
        return error.InvalidSubsurfaceSoluteState;
}

fn validateOrganic(flux: OrganicFlux, fraction_count: usize) !void {
    inline for (@typeInfo(OrganicFlux).@"struct".fields) |field| {
        const values = @field(flux, field.name);
        if (values.len != fraction_count)
            return error.InvalidSubsurfaceSoluteDimensions;
        for (values) |value|
            if (!std.math.isFinite(value))
                return error.InvalidSubsurfaceSoluteInput;
    }
}

fn validatePoreSalt(pore: PoreSaltFlux) !void {
    try validateFiniteStruct(pore.primary);
    try validateFiniteStruct(pore.secondary);
    try validateFiniteStruct(pore.tertiary);
    try validateFiniteStruct(pore.quaternary);
    try validateFiniteStruct(pore.quinary);
    try validateFiniteStruct(pore.phosphate_non_band);
    try validateFiniteStruct(pore.phosphate_band);
}

fn isActiveBoundary(
    inputs: Inputs,
    face: Face,
    column: usize,
    row: usize,
    cell: usize,
    layer: usize,
) bool {
    if (face == .lower)
        return layer + 1 == inputs.soil_layer_count;
    if (inputs.horizontal_exchange_by_cell[cell] == .standalone)
        return false;
    const window = inputs.external_boundary_window;
    return switch (face) {
        .east => column == window.last_column_inclusive,
        .west => column == window.first_column,
        .south => row == window.last_row_inclusive,
        .north => row == window.first_row,
        .lower => unreachable,
    };
}

fn inwardSign(face: Face) f64 {
    return switch (face) {
        .east, .south, .lower => -1,
        .west, .north => 1,
    };
}

fn signedStructSum(direction: f64, value: anytype) !f64 {
    return checkedProduct(direction, try sourceStructSum(value));
}

fn signedSourceSum(direction: f64, values: []const f64) !f64 {
    return checkedProduct(direction, try sourceSum(values));
}

fn sourceStructSum(value: anytype) !f64 {
    var total: f64 = 0;
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        total = try checkedSum(total, @field(value, field.name));
    return total;
}

fn sourceSum(values: []const f64) !f64 {
    var total: f64 = 0;
    for (values) |value| total = try checkedSum(total, value);
    return total;
}

fn checkedProduct(left: f64, right: f64) !f64 {
    const result = left * right;
    if (!std.math.isFinite(result))
        return error.NonFiniteSubsurfaceSoluteResult;
    return result;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!std.math.isFinite(result))
        return error.NonFiniteSubsurfaceSoluteResult;
    return result;
}

fn validateFiniteStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name)))
            return error.InvalidSubsurfaceSoluteInput;
}

fn validatePositiveStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!positiveFinite(@field(value, field.name)))
            return error.InvalidSubsurfaceSoluteInput;
}

fn positiveFinite(value: f64) bool {
    return std.math.isFinite(value) and value > 0;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn filledPoreSalt(value: f64) PoreSaltFlux {
    return .{
        .primary = filled(PrimaryIonFlux, value),
        .secondary = filled(SecondaryIonFlux, value),
        .tertiary = filled(TertiaryIonFlux, value),
        .quaternary = filled(QuaternaryIonFlux, value),
        .quinary = filled(QuinaryIonFlux, value),
        .phosphate_non_band = filled(PhosphateComplexFlux, value),
        .phosphate_band = filled(PhosphateComplexFlux, value),
    };
}

fn testFlux(value: f64, organic: []const f64) BoundaryFlux {
    return .{
        .micropore_water_m3_per_step = value,
        .macropore_water_m3_per_step = value,
        .organic = .{
            .micropore_carbon_g_c_per_step = organic,
            .micropore_acetate_carbon_g_c_per_step = organic,
            .macropore_carbon_g_c_per_step = organic,
            .macropore_acetate_carbon_g_c_per_step = organic,
            .micropore_nitrogen_g_n_per_step = organic,
            .macropore_nitrogen_g_n_per_step = organic,
            .micropore_phosphorus_g_p_per_step = organic,
            .macropore_phosphorus_g_p_per_step = organic,
        },
        .carbon = filled(CarbonFlux, value),
        .nitrogen = filled(NitrogenFlux, value),
        .phosphorus = filled(PhosphorusElementFlux, value),
        .gases = filled(GasFlux, value),
        .micropore_salts = filledPoreSalt(value),
        .macropore_salts = filledPoreSalt(value),
        .micropore_hydrogen_exchange_mol_per_step = value,
    };
}

fn unitConductivity() ConductivityCoefficients {
    return filled(ConductivityCoefficients, 1);
}

fn oneCellInputs(
    east: []const BoundaryFlux,
    empty: []const BoundaryFlux,
) Inputs {
    return .{
        .lon_count = 1,
        .lat_count = 1,
        .soil_layer_count = 1,
        .organic_fraction_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{.enabled},
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ east, empty, empty, empty, empty },
        .negligible_water_m3_by_cell = &.{0},
        .nitrogen_g_n_per_mol_n = 14,
        .phosphorus_g_p_per_mol_p = 31,
        .conductivity = unitConductivity(),
    };
}

fn zeroState(cells: []CellOutwardLedger, conductivity: []f64) State {
    return .{
        .cumulative = .{},
        .outward_by_cell = cells,
        .last_electrical_conductivity_dS_per_m_by_cell = conductivity,
        .source_hydrogen_diagnostic_g_h = 0,
    };
}

test "REDIST subsurface element and salt fixture preserves source weights" {
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const east = [_]BoundaryFlux{testFlux(1, &one)};
    const empty = [_]BoundaryFlux{testFlux(0, &zero)};
    var cells = [_]CellOutwardLedger{.{}};
    var conductivity = [_]f64{0};
    var state = zeroState(&cells, &conductivity);

    try apply(oneCellInputs(&east, &empty), &state);

    try std.testing.expectEqual(@as(f64, 10), state.cumulative.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 25), state.cumulative.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 1002), state.cumulative.phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 3), state.cumulative.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 3), state.cumulative.hydrogen_g_h);
    try std.testing.expectEqual(@as(f64, 238), state.cumulative.ion_components_mol);
    try std.testing.expectEqual(@as(f64, 4), cells[0].organic_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 6), cells[0].inorganic_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), cells[0].organic_nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 16), cells[0].inorganic_nitrogen_g_n);
    try std.testing.expectEqual(
        state.cumulative.nitrogen_g_n,
        cells[0].organic_nitrogen_g_n +
            cells[0].inorganic_nitrogen_g_n + 7,
    );
    try std.testing.expectEqual(
        state.cumulative.phosphorus_g_p,
        cells[0].organic_phosphorus_g_p +
            cells[0].inorganic_phosphorus_g_p + 992,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 20 + 1.0 / 14.0),
        conductivity[0],
        1.0e-14,
    );
    try std.testing.expectEqual(@as(f64, 2), state.source_hydrogen_diagnostic_g_h);
}

test "static salt mode leaves salt phosphorus ions and conductivity untouched" {
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const east = [_]BoundaryFlux{testFlux(1, &one)};
    const empty = [_]BoundaryFlux{testFlux(0, &zero)};
    var cells = [_]CellOutwardLedger{.{}};
    var conductivity = [_]f64{7};
    var state = zeroState(&cells, &conductivity);
    var inputs = oneCellInputs(&east, &empty);
    inputs.salt_equilibrium_mode = .static;
    try apply(inputs, &state);
    try std.testing.expectEqual(@as(f64, 10), state.cumulative.phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 0), state.cumulative.ion_components_mol);
    try std.testing.expectEqual(@as(f64, 7), conductivity[0]);
}

test "exact zero water gate suppresses every element and diagnostic" {
    const one = [_]f64{1};
    const zero = [_]f64{0};
    var cancelling_flux = testFlux(1, &one);
    cancelling_flux.macropore_water_m3_per_step = -1;
    const east = [_]BoundaryFlux{cancelling_flux};
    const empty = [_]BoundaryFlux{testFlux(0, &zero)};
    var cells = [_]CellOutwardLedger{.{}};
    var conductivity = [_]f64{0};
    var state = zeroState(&cells, &conductivity);
    try apply(oneCellInputs(&east, &empty), &state);
    try std.testing.expectEqualDeep(CumulativeOutwardLedger{}, state.cumulative);
    try std.testing.expectEqualDeep(CellOutwardLedger{}, cells[0]);
    try std.testing.expectEqual(@as(f64, 0), state.source_hydrogen_diagnostic_g_h);
}

test "corrected macropore FeOH4 and aluminum conductivity remain distinct" {
    const zero = [_]f64{0};
    var corrected = testFlux(0, &zero);
    corrected.micropore_water_m3_per_step = 1;
    corrected.macropore_salts.quinary.iron_hydroxide_4_mol_per_step = 2;
    corrected.macropore_salts.primary.aluminum_mol_per_step = 4;
    const east = [_]BoundaryFlux{corrected};
    const empty = [_]BoundaryFlux{testFlux(0, &zero)};
    var cells = [_]CellOutwardLedger{.{}};
    var conductivity = [_]f64{0};
    var state = zeroState(&cells, &conductivity);
    try apply(oneCellInputs(&east, &empty), &state);
    try std.testing.expectEqual(@as(f64, 14), state.cumulative.ion_components_mol);
    try std.testing.expectEqual(@as(f64, 12), conductivity[0]);
}

test "conductivity threshold zeros EC without suppressing element ledgers" {
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const east = [_]BoundaryFlux{testFlux(1, &one)};
    const empty = [_]BoundaryFlux{testFlux(0, &zero)};
    var cells = [_]CellOutwardLedger{.{}};
    var conductivity = [_]f64{9};
    var state = zeroState(&cells, &conductivity);
    var inputs = oneCellInputs(&east, &empty);
    inputs.negligible_water_m3_by_cell = &.{2};
    try apply(inputs, &state);
    try std.testing.expectEqual(@as(f64, 10), state.cumulative.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), conductivity[0]);
}

test "runtime grid closes carbon and ion ledgers against per-cell exports" {
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const active = testFlux(1, &one);
    const inactive = testFlux(0, &zero);
    const east = [_]BoundaryFlux{
        inactive, inactive, inactive, active,
        inactive, inactive, active,   inactive,
    };
    const empty = [_]BoundaryFlux{inactive} ** 8;
    var cells = [_]CellOutwardLedger{ .{}, .{}, .{}, .{} };
    var conductivity = [_]f64{ 0, 0, 0, 0 };
    var state = zeroState(&cells, &conductivity);
    try apply(.{
        .lon_count = 2,
        .lat_count = 2,
        .soil_layer_count = 2,
        .organic_fraction_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 1,
        },
        .horizontal_exchange_by_cell = &.{
            .enabled, .enabled, .enabled, .enabled,
        },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &east, &empty, &empty, &empty, &empty },
        .negligible_water_m3_by_cell = &.{ 0, 0, 0, 0 },
        .nitrogen_g_n_per_mol_n = 14,
        .phosphorus_g_p_per_mol_p = 31,
        .conductivity = unitConductivity(),
    }, &state);
    var carbon_g_c: f64 = 0;
    var ion_components_mol: f64 = 0;
    for (cells) |cell| {
        carbon_g_c += cell.organic_carbon_g_c + cell.inorganic_carbon_g_c;
        ion_components_mol += cell.ion_components_mol;
    }
    try std.testing.expectEqual(state.cumulative.carbon_g_c, carbon_g_c);
    try std.testing.expectEqual(state.cumulative.ion_components_mol, ion_components_mol);
}

test "late non-finite salt flux leaves all ledgers unchanged" {
    const one = [_]f64{1};
    var east = [_]BoundaryFlux{ testFlux(1, &one), testFlux(1, &one) };
    east[1].macropore_salts.quinary.iron_hydroxide_4_mol_per_step =
        std.math.nan(f64);
    const empty = [_]BoundaryFlux{ testFlux(0, &one), testFlux(0, &one) };
    var cells = [_]CellOutwardLedger{ .{ .organic_carbon_g_c = 2 }, .{} };
    var conductivity = [_]f64{ 3, 4 };
    var state: State = .{
        .cumulative = .{ .carbon_g_c = 5 },
        .outward_by_cell = &cells,
        .last_electrical_conductivity_dS_per_m_by_cell = &conductivity,
        .source_hydrogen_diagnostic_g_h = 6,
    };
    try std.testing.expectError(error.InvalidSubsurfaceSoluteInput, apply(.{
        .lon_count = 2,
        .lat_count = 1,
        .soil_layer_count = 1,
        .organic_fraction_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{ .enabled, .enabled },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &east, &empty, &empty, &empty, &empty },
        .negligible_water_m3_by_cell = &.{ 0, 0 },
        .nitrogen_g_n_per_mol_n = 14,
        .phosphorus_g_p_per_mol_p = 31,
        .conductivity = unitConductivity(),
    }, &state));
    try std.testing.expectEqual(@as(f64, 5), state.cumulative.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), cells[0].organic_carbon_g_c);
    try std.testing.expectEqualSlices(f64, &.{ 3, 4 }, &conductivity);
    try std.testing.expectEqual(@as(f64, 6), state.source_hydrogen_diagnostic_g_h);
}

test "invalid runtime organic extent fails before mutation" {
    const one = [_]f64{1};
    const east = [_]BoundaryFlux{testFlux(1, &one)};
    const empty = [_]BoundaryFlux{testFlux(0, &one)};
    var cells = [_]CellOutwardLedger{.{}};
    var conductivity = [_]f64{0};
    var state = zeroState(&cells, &conductivity);
    var inputs = oneCellInputs(&east, &empty);
    inputs.organic_fraction_count = 2;
    try std.testing.expectError(
        error.InvalidSubsurfaceSoluteDimensions,
        apply(inputs, &state),
    );
    try std.testing.expectEqualDeep(CumulativeOutwardLedger{}, state.cumulative);
}
