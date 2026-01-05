const std = @import("std");
const utils = @import("utils");
const parser = @import("input_parser");
const max_path_len = 1024;
///Lat-lon ranges and tile specifications
pub const LatLonRangeAndTileSpecs = packed struct {
    reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,
    lat_min_ud: i32 = undefined, //minimum latitude of the domain rounded to the nearest micro degree
    lat_max_ud: i32 = undefined, //maximum latitude of the domain rounded to the nearest micro degree
    lon_min_ud: i32 = undefined, //minimum longitude of the domain rounded to the nearest micro degree
    lon_max_ud: i32 = undefined, //maximum longitude of the domain rounded to the nearest micro degree
    dlat_ud: i32 = undefined, //latitude spacing to be simulated rounded to the nearest micro degree
    dlon_ud: i32 = undefined, //longitude spacing to be simulated rounded to the nearest micro degree
    ntx: usize = undefined, //number of west-east grids in a tile
    nty: usize = undefined, //number of north-south grids in a tile

    pub fn load(self: *LatLonRangeAndTileSpecs, file_name: []const u8) !void {
        var line = try parser.readNextDataLine(self.reader);
        var tokens = parser.Tokens{};
        try tokens.tokenizeLine(line, 6, "domain lat-lon range, and lat-lon spacing", file_name, self.err_log);
        const fields = [_]*i32{ &self.lat_min_ud, &self.lat_max_ud, &self.lon_min_ud, &self.lon_max_ud, &self.dlat_ud, &self.dlon_ud };
        for (tokens.items[0..fields.len], 0..) |tok, i| {
            const geo_loc_d = try parser.parseTokToFloat(f32, tok, "domain lat-lon range", file_name, self.err_log);
            const geo_loc_ud: i32 = @intFromFloat(@round(geo_loc_d * 1e6));
            fields[i].* = geo_loc_ud;
        }
        line = try parser.readNextDataLine(self.reader);
        try tokens.tokenizeLine(line, 2, "tile specifications", file_name, self.err_log);
        const fields_tile = [_]*usize{ &self.ntx, &self.nty };
        for (tokens.items[0..fields_tile.len], 0..) |tok, i| {
            fields_tile[i].* = try parser.parseTokToInt(usize, tok, "tile specifications", file_name, self.err_log);
            //check if tile specifications are valid: power-of-two
            try utils.requirePowerOfTwo(fields_tile[i].*);
        }
    }
};
test "LatLonRangeAndTileSpecs.load parses ranges and tile specs" {
    // Two lines: 6 floats, then 2 ints (power-of-two).
    const input = "10.5,12.25,-120.5,-118.25,0.5,0.25\n" ++
        "8,4\n";
    var input_reader = std.Io.Reader.fixed(input);

    var err_buf: [256]u8 = undefined;
    var err_log = std.Io.Writer.fixed(&err_buf);

    var specs = LatLonRangeAndTileSpecs{
        .reader = &input_reader,
        .err_log = &err_log,
    };

    try specs.load("test-input");

    try std.testing.expectEqual(@as(i32, 10_500_000), specs.lat_min_ud);
    try std.testing.expectEqual(@as(i32, 12_250_000), specs.lat_max_ud);
    try std.testing.expectEqual(@as(i32, -120_500_000), specs.lon_min_ud);
    try std.testing.expectEqual(@as(i32, -118_250_000), specs.lon_max_ud);
    try std.testing.expectEqual(@as(i32, 500_000), specs.dlat_ud);
    try std.testing.expectEqual(@as(i32, 250_000), specs.dlon_ud);
    try std.testing.expectEqual(@as(usize, 8), specs.ntx);
    try std.testing.expectEqual(@as(usize, 4), specs.nty);
}
///Simulation start year and initialization
pub const SimInit = packed struct {
    reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,
    start_yr: usize = undefined, //simulation start year
    sim_from_prev_run: usize = undefined, //does simulation start from a previous run? 0=no, 1=yes

    pub fn load(self: *SimInit, file_name: []const u8) !void {
        const line = try parser.readNextDataLine(self.reader);
        var tokens = parser.Tokens{};
        try tokens.tokenizeLine(line, 2, "simulation start year", file_name, self.err_log);
        const fields = [_]*usize{ &self.start_yr, &self.sim_from_prev_run };
        for (tokens.items[0..fields.len], 0..) |tok, i| {
            fields[i].* = try parser.parseTokToInt(usize, tok, "simulation start year", file_name, self.err_log);
        }
        try parser.boundsCheck(error.OutOfBounds, .{self.sim_from_prev_run > 1}, "whether the simulation starts from a previous run", file_name, self.err_log);
    }
};
test "SimInit.load parses simulation initialization params" {
    var input = "2001,0\n";
    var input_reader = std.Io.Reader.fixed(input);

    var err_buf: [256]u8 = undefined;
    var err_log = std.Io.Writer.fixed(&err_buf);

    var sim_init_params = SimInit{
        .reader = &input_reader,
        .err_log = &err_log,
    };

    try sim_init_params.load("test-input");

    try std.testing.expectEqual(@as(usize, 2001), sim_init_params.start_yr);
    try std.testing.expectEqual(@as(usize, 0), sim_init_params.sim_from_prev_run);

    input = "2001,2\n";
    input_reader = std.Io.Reader.fixed(input);

    try std.testing.expectError(error.OutOfBounds, sim_init_params.load("test-input"));
}
///Geographical attributes
pub const GeoAttr = packed struct {
    reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,
    bin_writer: *std.Io.Writer = undefined,
    lat_ud: i32 = undefined, //latitude rounded to the nearest micro degree
    lon_ud: i32 = undefined, //longitude rounded to the nearest micro degree
    elevation: f32 = undefined, //elevation/altitude (m)
    matc: f32 = undefined, //mean annual temperature (⁰C)
    ix: usize = 0, //longitude snapped to X coordinate id
    iy: usize = 0, //latitude snapped to Y coordinate id
    ix_old: usize = 0, //ix of the previous read
    iy_old: usize = 0, //iy of the previous read
    nx: usize = 0, //max domain ix
    ny: usize = 0, //max domain iy
    tile_ix: usize = 0, //tile id along X coordinate
    tile_iy: usize = 0, //tile id along Y coordinate
    local_ix: usize = 0, //local ix within a tile
    local_iy: usize = 0, //local iy within a tile

    pub fn loadAndWriteBin(self: *GeoAttr, lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs, file_name: []const u8) !void {
        var line = try parser.readNextDataLine(self.reader);
        var tokens = parser.Tokens{};
        try tokens.tokenizeLine(line, 1, "geographical attributes file path", file_name, self.err_log);
        var geo_buf: [max_path_len]u8 = undefined;
        @memcpy(geo_buf[0..tokens.items[0].len], tokens.items[0]);
        const geo_attr_filename: []const u8 = geo_buf[0..tokens.items[0].len];
        var geo_reader: utils.FileReader = utils.FileReader{};
        try geo_reader.open(self.err_log, geo_attr_filename);
        defer geo_reader.close();
        geo_reader.reader();

        try self.writeHeader();

        while (true) {
            line = try parser.readNextDataLine(geo_reader.buf_reader);
            if (std.mem.eql(u8, line, "EndOfStream")) break;
            try tokens.tokenizeLine(line, 4, "geographical attributes", geo_attr_filename, self.err_log);
            const fields = [_]*i32{ &self.lat_ud, &self.lon_ud };
            for (tokens.items[0..fields.len], 0..) |tok, i| {
                const geo_loc_d = try parser.parseTokToFloat(f32, tok, "latitude and longitude", geo_attr_filename, self.err_log);
                const geo_loc_ud: i32 = @intFromFloat(@round(geo_loc_d * 1e6));
                fields[i].* = geo_loc_ud;
            }
            const fields_elev_mat = [_]*f32{ &self.elevation, &self.matc };
            for (tokens.items[fields.len .. fields.len + fields_elev_mat.len], 0..) |tok, i| {
                fields_elev_mat[i].* = try parser.parseTokToFloat(f32, tok, "elevation and MAT", geo_attr_filename, self.err_log);
            }
            try self.assignId(lat_lon_rng_n_tile_specs);

            const morton: u64 = @intCast(utils.morton2D(self.ix, self.iy));
            try self.writeRecord(morton);
        }
    }
    fn assignId(self: *GeoAttr, lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs) !void {
        self.ix_old = self.ix;
        self.iy_old = self.iy;
        self.ix = utils.toXY(self.lon_ud, lat_lon_rng_n_tile_specs.lon_min_ud, lat_lon_rng_n_tile_specs.dlon_ud);
        self.iy = utils.toXY(self.lat_ud, lat_lon_rng_n_tile_specs.lat_min_ud, lat_lon_rng_n_tile_specs.dlat_ud);
        self.nx = @max(self.ix_old, self.ix);
        self.ny = @max(self.iy_old, self.iy);
        self.tile_ix = @divFloor(self.nx, lat_lon_rng_n_tile_specs.ntx);
        self.tile_iy = @divFloor(self.ny, lat_lon_rng_n_tile_specs.nty);
    }
    fn writeHeader(self: *GeoAttr) !void {
        try self.bin_writer.writeAll("geo_attr");
        try self.writeIntLittle(u16, 1); // version
        try self.writeIntLittle(u16, 1); // flags: little endian
    }
    fn writeIntLittle(self: *GeoAttr, comptime T: type, value: T) !void {
        var buf: [@sizeOf(T)]u8 = undefined;
        std.mem.writeInt(T, &buf, value, .little);
        try self.bin_writer.writeAll(&buf);
    }
    fn writeRecord(self: *GeoAttr, morton: u64) !void {
        try self.writeIntLittle(u64, morton);
        try self.writeIntLittle(i32, self.lat_ud);
        try self.writeIntLittle(i32, self.lon_ud);
        try self.writeIntLittle(u32, @as(u32, @bitCast(self.elevation)));
        try self.writeIntLittle(u32, @as(u32, @bitCast(self.matc)));
    }
};

pub const GeoAttrBinData = struct {
    arena: std.heap.ArenaAllocator = undefined,
    err_log: *std.Io.Writer = undefined,
    lat_ud: []i32 = undefined,
    lon_ud: []i32 = undefined,
    elev: []f32 = undefined,
    matc: []f32 = undefined,

    pub fn deinit(self: *GeoAttrBinData) void {
        self.arena.deinit();
    }

    pub fn readBinToArraysWithArena(self: *GeoAttrBinData, geo_attr: *const GeoAttr, bin_path: []const u8) !void {
        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        errdefer arena.deinit();
        const allocator = arena.allocator();

        // const nx_i64 = @divFloor(@as(i64, specs.lon_max_ud) - @as(i64, specs.lon_min_ud), @as(i64, specs.dlon_ud)) + 1;
        // const ny_i64 = @divFloor(@as(i64, specs.lat_max_ud) - @as(i64, specs.lat_min_ud), @as(i64, specs.dlat_ud)) + 1;
        // if (nx_i64 <= 0 or ny_i64 <= 0) return error.OutOfBounds;
        //
        // self.nx = @as(usize, @intCast(nx_i64));
        // self.ny = @as(usize, @intCast(ny_i64));
        const n = geo_attr.nx * geo_attr.ny;

        self.lat_ud = try allocator.alloc(i32, n);
        self.lon_ud = try allocator.alloc(i32, n);
        self.elev = try allocator.alloc(f32, n);
        self.matc = try allocator.alloc(f32, n);

        const lat_sentinel = std.math.minInt(i32);
        const lon_sentinel = std.math.minInt(i32);
        const elev_nan = std.math.nan(f32);
        const mat_nan = std.math.nan(f32);
        for (self.lat_ud, self.lon_ud, self.elev, self.matc) |*lat_p, *lon_p, *elev_p, *mat_p| {
            lat_p.* = lat_sentinel;
            lon_p.* = lon_sentinel;
            elev_p.* = elev_nan;
            mat_p.* = mat_nan;
        }

        var bin_reader: utils.FileReader = utils.FileReader{};
        try bin_reader.open(self.err_log, bin_path);
        defer bin_reader.close();
        bin_reader.reader();
        var buf_header: [12]u8 = undefined;
        try bin_reader.buf_reader.readSliceAll(&buf_header);
        std.debug.print("test bin header: {s}\n", .{buf_header[0..8]});

        if (!std.mem.eql(u8, buf_header[0..8], "geo_attr")) return error.BadHeader;
        const version = std.mem.readInt(u16, buf_header[8..10], .little);
        const flags = std.mem.readInt(u16, buf_header[10..12], .little);
        if (version != 1 or flags != 1) return error.BadHeader;
        var rec_buf: [24]u8 = undefined;
        while (true) {
            try bin_reader.buf_reader.readSliceAll(&rec_buf);
            std.debug.print("test bin lat_ud: {s}\n", .{rec_buf[0..4]});
        }
        //     if (nread == 0) break;
        //     if (nread != buf.len) return error.TruncatedRecord;
        //
        //     const lat_ud_ = std.mem.readInt(i32, buf[8..12], .little);
        //     const lon_ud_ = std.mem.readInt(i32, buf[12..16], .little);
        //     const elev_bits = std.mem.readInt(u32, buf[16..20], .little);
        //     const mat_bits = std.mem.readInt(u32, buf[20..24], .little);
        //
        //     const ix_i64 = @divFloor(@as(i64, lon_ud_) - @as(i64, specs.lon_min_ud), @as(i64, specs.dlon_ud));
        //     const iy_i64 = @divFloor(@as(i64, lat_ud_) - @as(i64, specs.lat_min_ud), @as(i64, specs.dlat_ud));
        //     if (ix_i64 < 0 or iy_i64 < 0) return error.OutOfBounds;
        //
        //     const ix = @as(usize, @intCast(ix_i64));
        //     const iy = @as(usize, @intCast(iy_i64));
        //     if (ix >= nx or iy >= ny) return error.OutOfBounds;
        //
        //     const flat_id = ix + iy * nx;
        //     lat_ud[flat_id] = lat_ud;
        //     lon_ud[flat_id] = lon_ud;
        //     elev[flat_id] = @as(f32, @bitCast(elev_bits));
        //     matc[flat_id] = @as(f32, @bitCast(mat_bits));
        // }
    }
};

// test "readBinToArraysWithArena loads lat/lon micro-degrees into row-major arrays" {
//     var tmp = std.testing.tmpDir(.{});
//     defer tmp.cleanup();
//
//     var orig_dir = try std.fs.cwd().openDir(".", .{});
//     defer orig_dir.close();
//     try tmp.dir.setAsCwd();
//     defer orig_dir.setAsCwd() catch {};
//
//     const geo_contents =
//         "1.0,2.0,100.5,5.0\n" ++
//         "3.0,4.0,200.25,-1.5\n";
//
//     var geo_file = try tmp.dir.createFile("geo_attr.txt", .{});
//     defer geo_file.close();
//     try geo_file.writeAll(geo_contents);
//
//     const geo_path = "geo_attr.txt";
//     const input = try std.fmt.allocPrint(std.testing.allocator, "{s}\n", .{geo_path});
//     defer std.testing.allocator.free(input);
//
//     var input_reader = std.Io.Reader.fixed(input);
//     var err_buf: [256]u8 = undefined;
//     var err_log = std.Io.Writer.fixed(&err_buf);
//
//     var bin_writer = utils.FileWriter{ .err_log = &err_log, .is_err_log = false };
//     try bin_writer.create("geo_attr.bin");
//     defer bin_writer.close();
//     bin_writer.writer();
//
//     var specs = LatLonRangeAndTileSpecs{
//         .reader = &input_reader,
//         .err_log = &err_log,
//         .lat_min_ud = 0,
//         .lat_max_ud = 3_000_000,
//         .lon_min_ud = 0,
//         .lon_max_ud = 4_000_000,
//         .dlat_ud = 1_000_000,
//         .dlon_ud = 1_000_000,
//         .ntx = 2,
//         .nty = 2,
//     };
//
//     var geo = GeoAttr{
//         .reader = &input_reader,
//         .err_log = &err_log,
//         .bin_writer = bin_writer.buf_writer,
//     };
//     try geo.loadAndWriteBin(&specs, "test-input");
//
//     var data = try readBinToArraysWithArena(&specs, "geo_attr.bin");
//     defer data.deinit();
//
//     const flat_id0 = 2 + 1 * data.nx;
//     try std.testing.expectEqual(@as(i32, 1_000_000), data.lat_ud[flat_id0]);
//     try std.testing.expectEqual(@as(i32, 2_000_000), data.lon_ud[flat_id0]);
//     try std.testing.expectApproxEqAbs(@as(f32, 100.5), data.elev[flat_id0], 1e-5);
//     try std.testing.expectApproxEqAbs(@as(f32, 5.0), data.matc[flat_id0], 1e-5);
// }
// test "GeoAttr.loadAndWriteBin writes header and records" {
//     var tmp = std.testing.tmpDir(.{});
//     defer tmp.cleanup();
//
//     var orig_dir = try std.fs.cwd().openDir(".", .{});
//     defer orig_dir.close();
//     try tmp.dir.setAsCwd();
//     defer orig_dir.setAsCwd() catch {};
//
//     const geo_contents =
//         "1.0,2.0,100.5,5.0\n" ++
//         "3.0,4.0,200.25,-1.5\n";
//
//     var geo_file = try tmp.dir.createFile("geo_attr.txt", .{});
//     defer geo_file.close();
//     try geo_file.writeAll(geo_contents);
//
//     const geo_path = "geo_attr.txt";
//
//     const input = try std.fmt.allocPrint(std.testing.allocator, "{s}\n", .{geo_path});
//     defer std.testing.allocator.free(input);
//
//     var input_reader = std.Io.Reader.fixed(input);
//
//     var err_buf: [256]u8 = undefined;
//     var err_log = std.Io.Writer.fixed(&err_buf);
//
//     var bin_buf: [128]u8 = undefined;
//     var bin_writer = std.Io.Writer.fixed(&bin_buf);
//
//     var specs = LatLonRangeAndTileSpecs{
//         .lat_min_ud = 0,
//         .lon_min_ud = 0,
//         .dlat_ud = 1_000_000,
//         .dlon_ud = 1_000_000,
//         .ntx = 2,
//         .nty = 2,
//     };
//
//     var geo = GeoAttr{
//         .reader = &input_reader,
//         .err_log = &err_log,
//         .bin_writer = &bin_writer,
//     };
//
//     try geo.loadAndWriteBin(&specs, "test-input");
//
//     const out = bin_writer.buffered();
//     try std.testing.expectEqual(@as(usize, 59), out.len);
//     try std.testing.expectEqualStrings("GEOATTR", out[0..7]);
//     try std.testing.expectEqual(@as(u16, 1), std.mem.readInt(u16, out[7..9], .little));
//     try std.testing.expectEqual(@as(u16, 1), std.mem.readInt(u16, out[9..11], .little));
//
//     const morton0 = std.mem.readInt(u64, out[11..19], .little);
//     const lat0 = std.mem.readInt(i32, out[19..23], .little);
//     const lon0 = std.mem.readInt(i32, out[23..27], .little);
//     try std.testing.expectEqual(@as(u64, @intCast(utils.morton2D(2, 1))), morton0);
//     try std.testing.expectEqual(@as(i32, 1_000_000), lat0);
//     try std.testing.expectEqual(@as(i32, 2_000_000), lon0);
// }
