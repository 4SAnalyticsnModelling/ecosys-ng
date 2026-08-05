const std = @import("std");

pub const Source = enum(u8) {
    precipitation = 1,
    irrigation = 2,
    soil = 3,
};

pub const ReactionRates = struct {
    aluminum_hydroxide_precipitation_mol_per_iteration: f64 = 0, // RPALOX
    iron_hydroxide_precipitation_mol_per_iteration: f64 = 0, // RPFEOX
    calcite_precipitation_mol_per_iteration: f64 = 0, // RPCACX
    gypsum_precipitation_mol_per_iteration: f64 = 0, // RPCASO
    aluminum_phosphate_precipitation_mol_per_iteration: f64 = 0, // RPALPX
    iron_phosphate_precipitation_mol_per_iteration: f64 = 0, // RPFEPX
    dicalcium_phosphate_precipitation_mol_per_iteration: f64 = 0, // RPCADX
    hydroxyapatite_precipitation_mol_per_iteration: f64 = 0, // RPCAHX
    protonated_to_hydroxyl_site_mol_per_iteration: f64 = 0, // RXOH2
    hydroxyl_to_deprotonated_site_mol_per_iteration: f64 = 0, // RXOH1
    h2po4_protonated_site_exchange_mol_per_iteration: f64 = 0, // RXH2P
    h2po4_hydroxyl_site_exchange_mol_per_iteration: f64 = 0, // RYH2P
    hpo4_hydroxyl_site_exchange_mol_per_iteration: f64 = 0, // RXH1P
};

pub const ExchangeableCations = struct {
    ammonium_mol_n_per_megagram: f64 = 0, // XN41
    hydrogen_mol_per_megagram: f64 = 0, // XHY1
    aluminum_mol_per_megagram: f64 = 0, // XAL1
    iron_mol_per_megagram: f64 = 0, // XFE1
    calcium_mol_per_megagram: f64 = 0, // XCA1
    magnesium_mol_per_megagram: f64 = 0, // XMG1
    sodium_mol_per_megagram: f64 = 0, // XNA1
    potassium_mol_per_megagram: f64 = 0, // XKA1
};

pub const Workspace = struct {
    rates: ReactionRates,
    exchangeable_cations: ExchangeableCations,
};

/// Direct translation of STARTE.F lines 780--802. Precipitation and
/// irrigation iterations must not retain soil-only mineral or exchange state.
/// The soil branch does not inspect or mutate this workspace here.
pub fn reset(source: Source, workspace: *Workspace) bool {
    if (source == .soil) return false;
    workspace.* = .{
        .rates = .{},
        .exchangeable_cations = .{},
    };
    return true;
}

fn filledWorkspace(value: f64) Workspace {
    var rates: ReactionRates = undefined;
    inline for (@typeInfo(ReactionRates).@"struct".fields) |field|
        @field(rates, field.name) = value;
    var cations: ExchangeableCations = undefined;
    inline for (@typeInfo(ExchangeableCations).@"struct".fields) |field|
        @field(cations, field.name) = value;
    return .{ .rates = rates, .exchangeable_cations = cations };
}

test "STARTE precipitation clears all soil-only equilibrium workspace in source order" {
    var workspace = filledWorkspace(7);
    try std.testing.expect(reset(.precipitation, &workspace));
    try std.testing.expectEqualDeep(Workspace{
        .rates = .{},
        .exchangeable_cations = .{},
    }, workspace);
}

test "STARTE irrigation clears the same soil-only workspace" {
    var workspace = filledWorkspace(-3);
    try std.testing.expect(reset(.irrigation, &workspace));
    inline for (@typeInfo(ReactionRates).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(workspace.rates, field.name));
    inline for (@typeInfo(ExchangeableCations).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(workspace.exchangeable_cations, field.name));
}

test "STARTE soil branch leaves equilibrium workspace untouched" {
    var workspace = filledWorkspace(std.math.nan(f64));
    const before = workspace;
    try std.testing.expect(!reset(.soil, &workspace));
    inline for (@typeInfo(ReactionRates).@"struct".fields) |field|
        try std.testing.expect(std.math.isNan(@field(workspace.rates, field.name)) and
            std.math.isNan(@field(before.rates, field.name)));
    inline for (@typeInfo(ExchangeableCations).@"struct".fields) |field|
        try std.testing.expect(std.math.isNan(@field(workspace.exchangeable_cations, field.name)) and
            std.math.isNan(@field(before.exchangeable_cations, field.name)));
}
