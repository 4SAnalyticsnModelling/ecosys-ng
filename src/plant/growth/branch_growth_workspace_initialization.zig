const std = @import("std");

pub const Remobilization = enum(u1) {
    disabled,
    enabled,
};

/// Branch-local scalar workspace reused during the source NB loop.
pub const Workspace = struct {
    shoot_remobilization_enabled: Remobilization,
    shoot_remobilization_active: Remobilization,
    root_remobilization_enabled: Remobilization,
    root_remobilization_active: Remobilization,
};

/// Exact grosub.f lines 551--557 initialization for one NB iteration.
///
/// The four source flags are scalar workspace, not branch-indexed state. The
/// caller invokes this operation in runtime NB order immediately before the
/// same branch's partition calculations. Carbon inputs and result are g C.
pub fn initializeBranch(
    workspace: *Workspace,
    leaf_carbon_g_c: f64,
    petiole_carbon_g_c: f64,
) !f64 {
    if (!std.math.isFinite(leaf_carbon_g_c) or
        !std.math.isFinite(petiole_carbon_g_c))
        return error.NonFiniteBranchCarbon;
    const total_carbon_g_c = leaf_carbon_g_c + petiole_carbon_g_c;
    if (!std.math.isFinite(total_carbon_g_c))
        return error.NonFiniteBranchCarbonTotal;

    workspace.* = .{
        .shoot_remobilization_enabled = .disabled,
        .shoot_remobilization_active = .disabled,
        .root_remobilization_enabled = .disabled,
        .root_remobilization_active = .disabled,
    };
    return @max(0.0, total_carbon_g_c);
}

test "GROSUB resets scalar flags for every runtime NB iteration" {
    var workspace: Workspace = undefined;
    const leaf = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    const petiole = [_]f64{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
    for (leaf, petiole, 0..) |leaf_g_c, petiole_g_c, branch| {
        workspace = .{
            .shoot_remobilization_enabled = .enabled,
            .shoot_remobilization_active = .enabled,
            .root_remobilization_enabled = .enabled,
            .root_remobilization_active = .enabled,
        };
        try std.testing.expectEqual(
            leaf_g_c + petiole_g_c,
            try initializeBranch(&workspace, leaf_g_c, petiole_g_c),
        );
        _ = branch;
        inline for (.{
            workspace.shoot_remobilization_enabled,
            workspace.shoot_remobilization_active,
            workspace.root_remobilization_enabled,
            workspace.root_remobilization_active,
        }) |flag| try std.testing.expectEqual(Remobilization.disabled, flag);
    }
}

test "GROSUB clamps only combined branch carbon" {
    var workspace: Workspace = undefined;
    try std.testing.expectEqual(
        @as(f64, 0),
        try initializeBranch(&workspace, -3.0, 1.0),
    );
    try std.testing.expectEqual(
        @as(f64, 1.5),
        try initializeBranch(&workspace, 2.0, -0.5),
    );
}

test "invalid branch carbon leaves scalar workspace unchanged" {
    var workspace: Workspace = .{
        .shoot_remobilization_enabled = .enabled,
        .shoot_remobilization_active = .enabled,
        .root_remobilization_enabled = .enabled,
        .root_remobilization_active = .enabled,
    };
    const before = workspace;
    try std.testing.expectError(
        error.NonFiniteBranchCarbon,
        initializeBranch(&workspace, std.math.nan(f64), 1.0),
    );
    try std.testing.expectEqualDeep(before, workspace);
}
