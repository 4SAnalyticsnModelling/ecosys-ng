const std = @import("std");

pub const StateSource = enum {
    new_simulation,
    checkpoint,
};

pub const ExecutionPass = enum {
    first,
    subsequent,
};

pub const CarbonLedger = struct {
    net_balance_g_c: f64,
    initial_sink_g_c: f64,
    cumulative_carbon_dioxide_exchange_g_c: f64,
    cumulative_aboveground_carbon_dioxide_exchange_g_c: f64,
    cumulative_uptake_g_c: f64,
    cumulative_sink_g_c: f64,
    oxidation_g_c: f64,
};

pub const NitrogenLedger = struct {
    initial_sink_g_n: f64,
    cumulative_uptake_g_n: f64,
    cumulative_sink_g_n: f64,
    cumulative_fixation_g_n: f64,
    current_ammonia_exchange_g_n: f64,
    cumulative_ammonia_exchange_g_n: f64,
    oxidation_g_n: f64,
};

pub const PhosphorusLedger = struct {
    initial_sink_g_p: f64,
    cumulative_uptake_g_p: f64,
    cumulative_sink_g_p: f64,
    oxidation_g_p: f64,
};

pub const HarvestLedger = struct {
    cumulative_carbon_g_c: f64,
    cumulative_nitrogen_g_n: f64,
    cumulative_phosphorus_g_p: f64,
    current_carbon_g_c: f64,
    current_nitrogen_g_n: f64,
    current_phosphorus_g_p: f64,
};

pub const ReseedingLedger = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const PlantLedger = struct {
    carbon: CarbonLedger,
    nitrogen: NitrogenLedger,
    phosphorus: PhosphorusLedger,
    harvest: HarvestLedger,
    reseeding: ReseedingLedger,
    cumulative_transpiration_m3: f64,
};

/// Direct translation of STARTQ lines 637--665. Balance history is cleared
/// only for a new simulation's first execution pass. Checkpoint restoration
/// and later passes retain every caller-provided ledger value.
pub fn initialize(
    ledgers: []PlantLedger,
    state_source: StateSource,
    execution_pass: ExecutionPass,
) !bool {
    if (ledgers.len == 0) return error.InvalidPlantBalanceLedgerDimensions;
    if (state_source != .new_simulation or execution_pass != .first)
        return false;
    @memset(ledgers, std.mem.zeroes(PlantLedger));
    return true;
}

fn filledLedger(value: f64) PlantLedger {
    var result: PlantLedger = undefined;
    inline for (std.meta.fields(PlantLedger)) |outer| {
        if (outer.type == f64) {
            @field(result, outer.name) = value;
        } else {
            inline for (std.meta.fields(outer.type)) |inner|
                @field(@field(result, outer.name), inner.name) = value;
        }
    }
    return result;
}

fn expectAll(value: PlantLedger, expected: f64) !void {
    inline for (std.meta.fields(PlantLedger)) |outer| {
        if (outer.type == f64) {
            try std.testing.expectEqual(expected, @field(value, outer.name));
        } else {
            inline for (std.meta.fields(outer.type)) |inner|
                try std.testing.expectEqual(
                    expected,
                    @field(@field(value, outer.name), inner.name),
                );
        }
    }
}

test "STARTQ new first pass zeros every runtime plant balance ledger" {
    var ledgers = [_]PlantLedger{
        filledLedger(1),
        filledLedger(2),
        filledLedger(3),
        filledLedger(4),
        filledLedger(5),
        filledLedger(6),
        filledLedger(7),
    };
    try std.testing.expect(try initialize(&ledgers, .new_simulation, .first));
    for (ledgers) |ledger| try expectAll(ledger, 0);
}

test "STARTQ checkpoint initialization preserves restored balances" {
    var ledgers = [_]PlantLedger{filledLedger(9)};
    try std.testing.expect(!try initialize(&ledgers, .checkpoint, .first));
    try expectAll(ledgers[0], 9);
}

test "STARTQ later execution pass preserves accumulated balances" {
    var ledgers = [_]PlantLedger{filledLedger(11)};
    try std.testing.expect(!try initialize(
        &ledgers,
        .new_simulation,
        .subsequent,
    ));
    try expectAll(ledgers[0], 11);
}
