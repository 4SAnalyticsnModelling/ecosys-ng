///This module contains input parser helper methods to be used throughout other moduels
const std = @import("std");
const utils = @import("utils");
const print = std.debug.print;
const max_path_len = 1024;
const max_io_buf = 10 * 1024;
const max_tok_num = 512;

///This struct helps check ecosys run submission arguments and parse the runfile name
pub const RunArg = struct {
    runfile_name: [max_path_len]u8 = undefined,
    file_reader: utils.FileReader = utils.FileReader{},

    ///This method gets the runfile name from run submission args
    pub fn getRunfile(self: *RunArg) ![]const u8 {
        const fba = std.heap.FixedBufferAllocator;
        var args_buf: [max_io_buf]u8 = undefined;
        var fba_args = fba.init(&args_buf);
        const allocator_args = fba_args.allocator();
        const args = try std.process.argsAlloc(allocator_args); //use of allocator is required for windows os. So allocator less option is not feasible in this case
        defer std.process.argsFree(allocator_args, args);
        if (args.len < 2) {
            const err = error.MissingArguments;
            print(
                "\x1b[1;31merror:\x1b[0m {s} occured during ecosys job submission. Correct format is: <path/to/ecosys/binary> <runfile>\n",
                .{@errorName(err)},
            );
            return err;
        }
        if (args[1].len >= self.runfile_name.len) {
            return error.RunfilePathTooLong;
        }
        @memcpy(self.runfile_name[0..args[1].len], args[1]);
        return self.runfile_name[0..args[1].len];
    }
};
///Logs an token count mismatch error to both err_log and stdout
fn logMismatch(comptime err: TokenBoundsErrors, err_log: *std.Io.Writer, context: []const u8, file_name: []const u8) TokenBoundsErrors!void {
    err_log.print("error: {s} while reading {s} in {s}\n", .{ @errorName(err), context, file_name }) catch {
        return error.PrintFailed;
    };
    err_log.flush() catch {
        return error.PrintFailed;
    };
    print("\x1b[1;31merror:\x1b[0m {s} while reading {s} in {s}\n", .{ @errorName(err), context, file_name });
}
///This method checks min max bounds for grids, plants, scenes etc.
pub fn boundsCheck(comptime err: TokenBoundsErrors, conds: anytype, context: []const u8, file_name: []const u8, err_log: *std.Io.Writer) TokenBoundsErrors!void {
    inline for (conds) |ok| { // comptile unrolling, so don't use this method if there's a lot of conditions
        if (ok) {
            try logMismatch(err, err_log, context, file_name);
            return err;
        }
    }
}
///This method logs parsing of strings => integer errors.
pub fn parseTokToInt(comptime T: type, tok: []const u8, context: []const u8, file_name: []const u8, err_log: *std.Io.Writer) TokenBoundsErrors!T {
    return std.fmt.parseInt(T, tok, 10) catch {
        const err = error.InvalidCharacter;
        try logMismatch(err, err_log, context, file_name);
        return err;
    };
}
test "parseTokToInt parses tokens to integers" {
    var buf: [255]u8 = undefined;
    var err_log = std.Io.Writer.fixed(&buf);
    var tok: []const u8 = "10";
    const tok_to_int = try parseTokToInt(usize, tok, "test string like number", "test file", &err_log);
    try std.testing.expect(tok_to_int == 10);
    tok = "10.0";
    try std.testing.expectError(error.InvalidCharacter, parseTokToInt(usize, tok, "test string like number", "test file", &err_log));
    try std.testing.expect(std.mem.containsAtLeast(u8, err_log.buffered(), 1, "error: InvalidCharacter while reading test string like number in test file"));
    tok = "10e2";
    try std.testing.expectError(error.InvalidCharacter, parseTokToInt(usize, tok, "test string like number", "test file", &err_log));
    tok = "10a";
    try std.testing.expectError(error.InvalidCharacter, parseTokToInt(usize, tok, "test string like number", "test file", &err_log));
}
///This method logs parsing of strings => float errors.
pub fn parseTokToFloat(comptime T: type, tok: []const u8, context: []const u8, file_name: []const u8, err_log: *std.Io.Writer) TokenBoundsErrors!T {
    return std.fmt.parseFloat(T, tok) catch {
        const err = error.InvalidCharacter;
        try logMismatch(err, err_log, context, file_name);
        return err;
    };
}
test "parseTokToFloat parses tokens to floating point numbers" {
    var buf: [255]u8 = undefined;
    var err_log = std.Io.Writer.fixed(&buf);
    var tok: []const u8 = "10";
    var tok_to_float = try parseTokToFloat(f32, tok, "test string like number", "test file", &err_log);
    try std.testing.expectApproxEqAbs(tok_to_float, 10.0, 0.0001);
    tok = "10.0";
    tok_to_float = try parseTokToFloat(f32, tok, "test string like number", "test file", &err_log);
    try std.testing.expectApproxEqAbs(tok_to_float, 10.0, 0.0001);
    tok = "10e2";
    tok_to_float = try parseTokToFloat(f32, tok, "test string like number", "test file", &err_log);
    try std.testing.expectApproxEqAbs(tok_to_float, 1000.0, 0.0001);
    tok = "10e-2";
    tok_to_float = try parseTokToFloat(f32, tok, "test string like number", "test file", &err_log);
    try std.testing.expectApproxEqAbs(tok_to_float, 0.1, 0.0001);
    tok = "0.001";
    const tok_to_f64 = try parseTokToFloat(f64, tok, "test string like number", "test file", &err_log);
    try std.testing.expectApproxEqAbs(tok_to_f64, 0.001, 0.0001);
    tok = "10a";
    try std.testing.expectError(error.InvalidCharacter, parseTokToInt(usize, tok, "test string like number", "test file", &err_log));
    try std.testing.expect(std.mem.containsAtLeast(u8, err_log.buffered(), 1, "error: InvalidCharacter while reading test string like number in test file"));
}
const TokenBoundsErrors = error{
    TooManyTokens,
    TokenCountMismatch,
    PrintFailed,
    OutOfBounds,
    InvalidCharacter,
    FilePathTooLong,
};
///This is a helper struct to tokenize items in a read line
pub const Tokens = struct {
    items: [max_tok_num][]const u8 = undefined, //using []const u8 since it's a pointer to actual line, not a storage, short-lived but memory efficient
    len: usize = 0,
    ///This method resets tokens
    pub fn reset(self: *Tokens) void {
        self.len = 0;
        self.items = undefined;
    }
    ///This method parses a line into a list of numbers/strings (called tokens hereafter) by eliminating spaces, tabs, or commas between tokens
    pub fn tokenizeLine(self: *Tokens, line: []const u8, expected: usize, context: []const u8, file_name: []const u8, err_log: *std.Io.Writer) TokenBoundsErrors!void {
        self.reset();
        const hash_idx = std.mem.indexOfScalar(u8, line, '#');
        const data_before_hash = if (hash_idx) |i| line[0..i] else line;
        var it = std.mem.tokenizeAny(u8, data_before_hash, " ,\t\r\n");
        while (it.next()) |tok| {
            try boundsCheck(error.TooManyTokens, .{self.len >= self.items.len}, context, file_name, err_log);
            self.items[self.len] = tok;
            self.len += 1;
        }
        try boundsCheck(error.TokenCountMismatch, .{self.len != expected}, context, file_name, err_log);
    }
};
test "tokenizeLine splits lines into tokens" {
    var buf: [255]u8 = undefined;
    var err_log = std.Io.Writer.fixed(&buf);
    var tokens = Tokens{};
    const fake_line1 = "0 0.1 1 10 200.0 5e-2 10e8";
    try tokens.tokenizeLine(fake_line1, 7, "fake line #1", "test file", &err_log);
    try std.testing.expect(tokens.len == 7);
    try std.testing.expectEqualStrings(tokens.items[3], "10");
    try std.testing.expectEqualStrings(tokens.items[5], "5e-2");
    const fake_line2 = "0 0.1 1 10 \r\n200.0 5e-2 10e8";
    try tokens.tokenizeLine(fake_line2, 7, "fake line #2", "test file", &err_log);
    try std.testing.expect(tokens.len == 7);
    try std.testing.expectEqualStrings(tokens.items[3], "10");
    try std.testing.expectEqualStrings(tokens.items[5], "5e-2");
    try std.testing.expectError(error.TokenCountMismatch, tokens.tokenizeLine(fake_line2, 8, "fake line #2", "test file", &err_log));
    try std.testing.expect(std.mem.containsAtLeast(u8, err_log.buffered(), 1, "error: TokenCountMismatch"));
    try std.testing.expect(std.mem.containsAtLeast(u8, err_log.buffered(), 1, "in test file"));
}
///Skip line starts with # or blank line while reading only the data lines. Also break out if end of the file reaches
pub fn readNextDataLine(reader: *std.Io.Reader) ![]const u8 {
    while (true) {
        const line = reader.takeDelimiterInclusive('\n') catch |err| switch (err) {
            error.EndOfStream => return "EndOfStream",
            else => return err,
        };
        const trimmed = std.mem.trimLeft(u8, line, " ,\t\r\n");
        if (trimmed.len == 0 or trimmed[0] == '#') continue;
        return line;
    }
}
test "readNextDataLine skips comments and blanks" {
    const input =
        "# full line comment\n" ++
        "   # indented comment\r\n" ++
        "\n" ++
        "\t \r\n" ++
        "1 2 3\n";
    var reader = std.Io.Reader.fixed(input);
    var line = try readNextDataLine(&reader);
    try std.testing.expectEqualStrings("1 2 3\n", line);
    line = try readNextDataLine(&reader);
    try std.testing.expectEqualStrings("EndOfStream", line);
}
