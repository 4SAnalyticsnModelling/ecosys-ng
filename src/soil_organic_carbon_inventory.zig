const std = @import("std");

pub const DisturbanceMode = enum {
    no_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter_change,
    freeze_thaw_erosion_and_organic_matter_change,

    fn changesOrganicMatter(self: DisturbanceMode) bool {
        return self == .freeze_thaw_and_organic_matter_change or
            self == .freeze_thaw_erosion_and_organic_matter_change;
    }
};

pub const Dimensions = struct {
    microbial_substrate_count: usize,
    microbial_group_count: usize,
    microbial_component_count: usize,
    residue_class_count: usize,
    residue_component_count: usize,
    surface_organic_component_count: usize,
};

/// Carbon masses for one runtime soil layer, all in g C.
pub const CarbonPools = struct {
    microbial_biomass_g_c: []const f64, // OMC, order K then N then M
    microbial_residue_g_c: []const f64, // ORC, order K then M
    micropore_dissolved_organic_g_c: []const f64, // OQC
    macropore_dissolved_organic_g_c: []const f64, // OQCH
    adsorbed_organic_g_c: []const f64, // OHC
    micropore_acetate_g_c: []const f64, // OQA
    macropore_acetate_g_c: []const f64, // OQAH
    adsorbed_acetate_g_c: []const f64, // OHA
    surface_organic_g_c: []const f64, // OSC, order K then M
};

pub const InventoryError = error{
    DimensionMismatch,
    NonFiniteCarbonMass,
    NegativeCarbonMass,
    NonFiniteInventory,
};

/// Translates HOUR1 lines 3216-3247. Returns null when the selected
/// disturbance mode does not require an SOC-change baseline.
pub fn calculate(
    mode: DisturbanceMode,
    dimensions: Dimensions,
    pools: CarbonPools,
) InventoryError!?f64 {
    if (!mode.changesOrganicMatter()) return null;

    const microbial_count = dimensions.microbial_substrate_count *
        dimensions.microbial_group_count * dimensions.microbial_component_count;
    const residue_count = dimensions.residue_class_count *
        dimensions.residue_component_count;
    const surface_count = dimensions.residue_class_count *
        dimensions.surface_organic_component_count;
    if (pools.microbial_biomass_g_c.len != microbial_count or
        pools.microbial_residue_g_c.len != residue_count or
        pools.micropore_dissolved_organic_g_c.len != dimensions.residue_class_count or
        pools.macropore_dissolved_organic_g_c.len != dimensions.residue_class_count or
        pools.adsorbed_organic_g_c.len != dimensions.residue_class_count or
        pools.micropore_acetate_g_c.len != dimensions.residue_class_count or
        pools.macropore_acetate_g_c.len != dimensions.residue_class_count or
        pools.adsorbed_acetate_g_c.len != dimensions.residue_class_count or
        pools.surface_organic_g_c.len != surface_count)
    {
        return error.DimensionMismatch;
    }
    inline for (std.meta.fields(CarbonPools)) |field| {
        for (@field(pools, field.name)) |mass_g_c| {
            if (!std.math.isFinite(mass_g_c)) return error.NonFiniteCarbonMass;
            if (mass_g_c < 0.0) return error.NegativeCarbonMass;
        }
    }

    var total_organic_carbon_g: f64 = 0.0;
    for (0..dimensions.microbial_substrate_count) |substrate| {
        for (0..dimensions.microbial_group_count) |group| {
            for (0..dimensions.microbial_component_count) |component| {
                const index = (substrate * dimensions.microbial_group_count + group) *
                    dimensions.microbial_component_count + component;
                total_organic_carbon_g =
                    total_organic_carbon_g + pools.microbial_biomass_g_c[index];
            }
        }
    }
    for (0..dimensions.residue_class_count) |class| {
        for (0..dimensions.residue_component_count) |component| {
            const index = class * dimensions.residue_component_count + component;
            total_organic_carbon_g =
                total_organic_carbon_g + pools.microbial_residue_g_c[index];
        }
        total_organic_carbon_g = total_organic_carbon_g +
            pools.micropore_dissolved_organic_g_c[class] +
            pools.macropore_dissolved_organic_g_c[class] +
            pools.adsorbed_organic_g_c[class] +
            pools.micropore_acetate_g_c[class] +
            pools.macropore_acetate_g_c[class] +
            pools.adsorbed_acetate_g_c[class];
        for (0..dimensions.surface_organic_component_count) |component| {
            const index = class * dimensions.surface_organic_component_count + component;
            total_organic_carbon_g =
                total_organic_carbon_g + pools.surface_organic_g_c[index];
        }
    }
    if (!std.math.isFinite(total_organic_carbon_g)) return error.NonFiniteInventory;
    return total_organic_carbon_g;
}

test "SOC inventory preserves legacy nested accumulation order" {
    const dimensions = Dimensions{
        .microbial_substrate_count = 2,
        .microbial_group_count = 2,
        .microbial_component_count = 2,
        .residue_class_count = 2,
        .residue_component_count = 2,
        .surface_organic_component_count = 3,
    };
    const microbial = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8 };
    const residue = [_]f64{ 9, 10, 11, 12 };
    const per_class = [_]f64{ 1, 2 };
    const surface = [_]f64{ 1, 2, 3, 4, 5, 6 };
    const total = (try calculate(
        .freeze_thaw_and_organic_matter_change,
        dimensions,
        .{
            .microbial_biomass_g_c = &microbial,
            .microbial_residue_g_c = &residue,
            .micropore_dissolved_organic_g_c = &per_class,
            .macropore_dissolved_organic_g_c = &per_class,
            .adsorbed_organic_g_c = &per_class,
            .micropore_acetate_g_c = &per_class,
            .macropore_acetate_g_c = &per_class,
            .adsorbed_acetate_g_c = &per_class,
            .surface_organic_g_c = &surface,
        },
    )).?;
    try std.testing.expectEqual(@as(f64, 117.0), total);
}

test "mode without organic matter change skips inventory" {
    const result = try calculate(.freeze_thaw_and_erosion, undefined, undefined);
    try std.testing.expect(result == null);
}
