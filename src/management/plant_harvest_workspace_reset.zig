const std = @import("std");

pub const ElementMass = struct {
    carbon_g_c_per_h: f64,
    nitrogen_g_n_per_h: f64,
    phosphorus_g_p_per_h: f64,
};

pub const Workspace = struct {
    remaining_plant_fraction: f64,
    branch_leaf_harvest_threshold_g_c: f64,
    shoot_harvest_by_class: []ElementMass,
    root_harvest_by_class: []ElementMass,
    auxiliary_harvest_by_class: []ElementMass,
    grain_harvest: ElementMass,
};

/// Exact grosub.f lines 311--357 per-plant harvest scratch initialization.
///
/// Callers allocate every class array at runtime. The three element fields
/// preserve each source C, N, P assignment group and their extensive hourly
/// units. This workspace is private to one plant kernel and contains no global
/// mutable state.
pub fn reset(workspace: *Workspace) !void {
    if (workspace.shoot_harvest_by_class.len == 0 or
        workspace.root_harvest_by_class.len == 0 or
        workspace.auxiliary_harvest_by_class.len == 0)
        return error.EmptyPlantHarvestWorkspace;

    workspace.remaining_plant_fraction = 1;
    workspace.branch_leaf_harvest_threshold_g_c = 0;
    for (workspace.shoot_harvest_by_class) |*mass| resetMass(mass);
    for (workspace.root_harvest_by_class) |*mass| resetMass(mass);
    for (workspace.auxiliary_harvest_by_class) |*mass| resetMass(mass);
    resetMass(&workspace.grain_harvest);
}

fn resetMass(mass: *ElementMass) void {
    mass.carbon_g_c_per_h = 0;
    mass.nitrogen_g_n_per_h = 0;
    mass.phosphorus_g_p_per_h = 0;
}

test "GROSUB resets every harvest class in source element order" {
    var shoot = [_]ElementMass{filled(2)} ** 5;
    var root = [_]ElementMass{filled(3)} ** 4;
    var auxiliary = [_]ElementMass{filled(4)} ** 5;
    var workspace: Workspace = .{
        .remaining_plant_fraction = 0.25,
        .branch_leaf_harvest_threshold_g_c = 9,
        .shoot_harvest_by_class = &shoot,
        .root_harvest_by_class = &root,
        .auxiliary_harvest_by_class = &auxiliary,
        .grain_harvest = filled(5),
    };

    try reset(&workspace);

    try std.testing.expectEqual(@as(f64, 1), workspace.remaining_plant_fraction);
    try std.testing.expectEqual(@as(f64, 0), workspace.branch_leaf_harvest_threshold_g_c);
    for (workspace.shoot_harvest_by_class) |mass| try expectZero(mass);
    for (workspace.root_harvest_by_class) |mass| try expectZero(mass);
    for (workspace.auxiliary_harvest_by_class) |mass| try expectZero(mass);
    try expectZero(workspace.grain_harvest);
}

test "runtime harvest class counts are not compile-time ceilings" {
    const allocator = std.testing.allocator;
    const shoot = try allocator.alloc(ElementMass, 17);
    defer allocator.free(shoot);
    const root = try allocator.alloc(ElementMass, 23);
    defer allocator.free(root);
    const auxiliary = try allocator.alloc(ElementMass, 19);
    defer allocator.free(auxiliary);
    @memset(shoot, filled(1));
    @memset(root, filled(1));
    @memset(auxiliary, filled(1));
    var workspace: Workspace = .{
        .remaining_plant_fraction = 0,
        .branch_leaf_harvest_threshold_g_c = 1,
        .shoot_harvest_by_class = shoot,
        .root_harvest_by_class = root,
        .auxiliary_harvest_by_class = auxiliary,
        .grain_harvest = filled(1),
    };
    try reset(&workspace);
    try expectZero(shoot[16]);
    try expectZero(root[22]);
    try expectZero(auxiliary[18]);
}

test "empty late class domain leaves workspace unchanged" {
    var shoot = [_]ElementMass{filled(2)};
    var root = [_]ElementMass{filled(3)};
    var workspace: Workspace = .{
        .remaining_plant_fraction = 0.5,
        .branch_leaf_harvest_threshold_g_c = 7,
        .shoot_harvest_by_class = &shoot,
        .root_harvest_by_class = &root,
        .auxiliary_harvest_by_class = &.{},
        .grain_harvest = filled(4),
    };
    const shoot_before = shoot;
    try std.testing.expectError(error.EmptyPlantHarvestWorkspace, reset(&workspace));
    try std.testing.expectEqual(@as(f64, 0.5), workspace.remaining_plant_fraction);
    try std.testing.expectEqualSlices(ElementMass, &shoot_before, &shoot);
}

fn filled(value: f64) ElementMass {
    return .{
        .carbon_g_c_per_h = value,
        .nitrogen_g_n_per_h = value,
        .phosphorus_g_p_per_h = value,
    };
}

fn expectZero(mass: ElementMass) !void {
    inline for (@typeInfo(ElementMass).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(mass, field.name));
}
