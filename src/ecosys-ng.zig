const std = @import("std");
const utils = @import("utils");
const parser = @import("input_parser");
const err_handler = @import("err_handler");
const load_run = @import("load_run");
///ecosys-ng main function
pub fn main() !void {
    const start_time_us: i64 = std.time.microTimestamp();
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    //Create directory tree for saving the outputs
    var out: utils.OutDir = utils.OutDir{};
    try out.mkOutDirs();
    //Create the error log file in the outputs folder
    var err_log: utils.ErrorLog = utils.ErrorLog{};
    try err_log.init(&out);
    defer err_log.file_writer.close();
    err_log.file_writer.writer();
    //Initialize run args to store the runfile name from the run submission
    var run: parser.RunArg = parser.RunArg{};
    //Read and check ecosys run submission arguments and read ecosys runfile/runscript name from run submission arguments
    const runfile = try run.getRunfile();
    //Open runfile to read.
    try run.file_reader.open(err_log.file_writer.buf_writer, runfile);
    //Close the runfile later before the main function returns
    defer run.file_reader.close();
    //Open a buffered reader to read runfile
    run.file_reader.reader();
    //Initialize model runtime and completion tracking
    var runtime = err_handler.CompletionTime.init(start_time_us, runfile, err_log.file_writer.buf_writer);
    //Generate run failure message, if the run fails before completion
    errdefer runtime.fail();
    var load_runfile = load_run.LoadRun{
        .allocator = allocator,
        .runfile = runfile,
        .err_log = err_log.file_writer.buf_writer,
        .buf_reader = run.file_reader.buf_reader,
    };
    try load_runfile.load();
    try runtime.success();
}
