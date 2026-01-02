const std = @import("std");
const print = std.debug.print;
///Returns model run completion time with success message if the run encounters no error(s), otherwise it returns model run failure message
pub const CompletionTime = struct {
    start_time_us: i64,
    runfile: []const u8,
    err_log: *std.Io.Writer,
    ///This method initalize the CompletionTime struct
    pub fn init(start_time_us: i64, runfile: []const u8, err_log: *std.Io.Writer) CompletionTime {
        return .{
            .start_time_us = start_time_us,
            .runfile = runfile,
            .err_log = err_log,
        };
    }
    ///This method confirms successful completion of a model run with time required to complete the run
    pub fn success(self: *CompletionTime) !void {
        const end = std.time.microTimestamp();
        const elapsed = self.elapsedTime(end);
        try self.err_log.print("success: Ecosys model run in {s} completed in {d} {s}!\n", .{ self.runfile, elapsed.value, elapsed.unit });
        try self.err_log.flush();
        print("\x1b[1;32msuccess:\x1b[0m Ecosys model run in {s} completed in {d} {s}!\n", .{ self.runfile, elapsed.value, elapsed.unit });
    }
    ///This method releases model failure errors
    pub fn fail(self: *CompletionTime) void {
        print("\x1b[1;31merror:\x1b[0m Ecosys model run in {s} failed!\n", .{self.runfile});
        self.err_log.print("error: Ecosys model run in {s} failed!\n", .{self.runfile}) catch {};
        self.err_log.flush() catch {};
    }
    ///This method calculates elapsed time from start to finish of a model run
    fn elapsedTime(self: *CompletionTime, end_time_us: i64) struct { value: f64, unit: []const u8 } {
        const stdt = std.time;
        const elapsed_us: i64 = end_time_us - self.start_time_us;
        const e = @as(f64, @floatFromInt(elapsed_us));
        if (elapsed_us < stdt.us_per_ms) return .{ .value = e, .unit = "us" };
        if (elapsed_us < stdt.us_per_s) return .{ .value = e / stdt.us_per_ms, .unit = "ms" };
        if (elapsed_us < stdt.us_per_min) return .{ .value = e / stdt.us_per_s, .unit = "s" };
        if (elapsed_us < stdt.us_per_hour) return .{ .value = e / stdt.us_per_min, .unit = "min" };
        if (elapsed_us < stdt.us_per_day) return .{ .value = e / stdt.us_per_hour, .unit = "hr" };
        if (elapsed_us < stdt.us_per_week) return .{ .value = e / stdt.us_per_day, .unit = "d" };
        return .{ .value = e / stdt.us_per_week, .unit = "wk" };
    }
};
