const std = @import("std");
const nitrogen_parameters = @import("../nutrients/nitrogen_parameters.zig");
const fermenter = @import("fermenter_respiration.zig");
const acetotrophic = @import("acetotrophic_methanogenesis.zig");

/// Runtime-owned parameter expansion for NITRO.F 900--1048.
///
/// One unit represents one model-owned microbial population entry. This
/// workspace contains no scientific state and can be reused across timesteps
/// while population roles and parsed parameters remain unchanged.
pub const Workspace = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    population_index: []usize,
    fermenter_enabled: []bool,
    fermenter_reference_energy_yield_kj_per_g_c: []f64,
    fermenter_growth_energy_requirement_kj_per_g_c: []f64,
    fermenter_minimum_respiration_requirement_g_c_per_g_c: []f64,
    fermenter_specific_oxidation_rate_g_c_per_g_c_h: []f64,
    fermenter_doc_half_saturation_g_c_per_m3: []f64,
    fermenter_acetate_product_inhibition_g_c_per_m3: []f64,
    acetotrophic_enabled: []bool,
    acetotrophic_acetate_product_inhibition_g_c_per_m3: []f64,
    acetotrophic_acetate_half_saturation_g_c_per_m3: []f64,
    acetotrophic_specific_respiration_rate_g_c_per_g_c_h: []f64,
    acetotrophic_reference_energy_yield_kj_per_g_c: []f64,
    acetotrophic_growth_energy_requirement_kj_per_g_c: []f64,
    acetotrophic_minimum_growth_respiration_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, population_index: []const usize) !Workspace {
        if (population_index.len == 0) return error.InvalidAnaerobicEnergyWorkspaceDimensions;
        var workspace: Workspace = undefined;
        workspace.allocator = allocator;
        workspace.unit_count = population_index.len;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(Workspace).@"struct".fields) |field| {
                if ((field.type == []usize or field.type == []bool or field.type == []f64) and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(workspace, field.name));
                }
            }
        }
        inline for (@typeInfo(Workspace).@"struct".fields) |field| {
            if (field.type == []usize) {
                @field(workspace, field.name) = try allocator.dupe(usize, population_index);
                allocated += 1;
            } else if (field.type == []bool) {
                @field(workspace, field.name) = try allocator.alloc(bool, population_index.len);
                @memset(@field(workspace, field.name), false);
                allocated += 1;
            } else if (field.type == []f64) {
                @field(workspace, field.name) = try allocator.alloc(f64, population_index.len);
                @memset(@field(workspace, field.name), 0);
                allocated += 1;
            }
        }
        return workspace;
    }

    pub fn deinit(self: *Workspace) void {
        inline for (@typeInfo(Workspace).@"struct".fields) |field|
            if (field.type == []usize or field.type == []bool or field.type == []f64)
                self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn populate(
        self: *Workspace,
        parameters: nitrogen_parameters.AnaerobicEnergyParameters,
    ) !void {
        try fermenter.fillParameterArrays(
            parameters.fermenter,
            self.population_index,
            self.fermenterArrays(),
        );
        try acetotrophic.fillParameterArrays(
            parameters.acetotrophic_methanogenesis,
            self.population_index,
            self.acetotrophicArrays(),
        );
    }

    pub fn fermenterArrays(self: *Workspace) fermenter.ParameterArrays {
        return .{
            .enabled = self.fermenter_enabled,
            .reference_energy_yield_kj_per_g_c = self.fermenter_reference_energy_yield_kj_per_g_c,
            .growth_energy_requirement_kj_per_g_c = self.fermenter_growth_energy_requirement_kj_per_g_c,
            .minimum_respiration_requirement_g_c_per_g_c = self.fermenter_minimum_respiration_requirement_g_c_per_g_c,
            .specific_oxidation_rate_g_c_per_g_c_h = self.fermenter_specific_oxidation_rate_g_c_per_g_c_h,
            .dissolved_organic_carbon_half_saturation_g_c_per_m3 = self.fermenter_doc_half_saturation_g_c_per_m3,
            .acetate_product_inhibition_g_c_per_m3 = self.fermenter_acetate_product_inhibition_g_c_per_m3,
        };
    }

    pub fn acetotrophicArrays(self: *Workspace) acetotrophic.ParameterArrays {
        return .{
            .enabled = self.acetotrophic_enabled,
            .acetate_product_inhibition_g_c_per_m3 = self.acetotrophic_acetate_product_inhibition_g_c_per_m3,
            .acetate_half_saturation_g_c_per_m3 = self.acetotrophic_acetate_half_saturation_g_c_per_m3,
            .specific_respiration_rate_g_c_per_g_c_h = self.acetotrophic_specific_respiration_rate_g_c_per_g_c_h,
            .reference_energy_yield_kj_per_g_c = self.acetotrophic_reference_energy_yield_kj_per_g_c,
            .growth_energy_requirement_kj_per_g_c = self.acetotrophic_growth_energy_requirement_kj_per_g_c,
            .minimum_growth_respiration_fraction = self.acetotrophic_minimum_growth_respiration_fraction,
        };
    }
};

test "runtime workspace expands arbitrary population catalog without fixed ceiling" {
    const populations = [_]usize{ 0, 3, 4, 6, 8, 15, 23, 42, 64, 99 };
    var workspace = try Workspace.init(std.testing.allocator, &populations);
    defer workspace.deinit();
    try workspace.populate(try nitrogen_parameters.sourceAnaerobicEnergyParameters());

    try std.testing.expectEqual(populations.len, workspace.unit_count);
    try std.testing.expectEqualSlices(usize, &populations, workspace.population_index);
    try std.testing.expect(workspace.fermenter_enabled[1]);
    try std.testing.expect(workspace.fermenter_enabled[3]);
    try std.testing.expect(!workspace.fermenter_enabled[2]);
    try std.testing.expect(workspace.acetotrophic_enabled[2]);
    try std.testing.expect(!workspace.acetotrophic_enabled[1]);
    try std.testing.expectEqual(@as(f64, 0.5), workspace.fermenter_minimum_respiration_requirement_g_c_per_g_c[3]);
    try std.testing.expectEqual(@as(f64, 37.5), workspace.acetotrophic_growth_energy_requirement_kj_per_g_c[2]);
}

test "runtime workspace can be repopulated from changed user roles and values" {
    var workspace = try Workspace.init(std.testing.allocator, &.{ 2, 7, 12 });
    defer workspace.deinit();
    const source =
        "soil_fermenter_respiration 7 12 0.2 13 14 4 38 39 0.41 0.51 0.008 2.1 73\n" ++
        "soil_acetotrophic_methanogenesis 2 15 16 0.3 1.6 40 0.42 0.009 25 0.6";
    try workspace.populate(try nitrogen_parameters.parseAnaerobicEnergyParameters(source));

    try std.testing.expectEqualSlices(bool, &.{ false, true, true }, workspace.fermenter_enabled);
    try std.testing.expectEqualSlices(bool, &.{ true, false, false }, workspace.acetotrophic_enabled);
    try std.testing.expectEqual(@as(f64, 38), workspace.fermenter_growth_energy_requirement_kj_per_g_c[1]);
    try std.testing.expectEqual(@as(f64, 39), workspace.fermenter_growth_energy_requirement_kj_per_g_c[2]);
    try std.testing.expectEqual(@as(f64, 1.6), workspace.acetotrophic_reference_energy_yield_kj_per_g_c[0]);
}

test "runtime workspace rejects empty population catalog" {
    try std.testing.expectError(
        error.InvalidAnaerobicEnergyWorkspaceDimensions,
        Workspace.init(std.testing.allocator, &.{}),
    );
}
